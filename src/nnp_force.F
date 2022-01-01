!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with Neural Network potentials
!> \author Christoph Schran (christoph.schran@rub.de)
!> \date   2020-10-10
! **************************************************************************************************
MODULE nnp_force

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_sum
   USE nnp_acsf,                        ONLY: nnp_calc_acsf
   USE nnp_environment_types,           ONLY: nnp_env_get,&
                                              nnp_type
   USE nnp_model,                       ONLY: nnp_gradients,&
                                              nnp_predict
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE physcon,                         ONLY: angstrom
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nnp_force'

   PUBLIC :: nnp_calc_energy_force

CONTAINS

! **************************************************************************************************
!> \brief Calculate the energy and force for a given configuration with the NNP
!> \param nnp ...
!> \param calc_forces ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_calc_energy_force(nnp, calc_forces)
      TYPE(nnp_type), INTENT(INOUT), POINTER             :: nnp
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_calc_energy_force', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, i_com, ig, ind, istart, j, k, &
                                                            m, mecalc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: allcalc
      LOGICAL                                            :: calc_stress
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: denergydsym
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: dsymdxyz, stress
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (particle_set, logger, local_particles, subsys, &
               atomic_kind_set)
      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(nnp))
      CALL nnp_env_get(nnp_env=nnp, particle_set=particle_set, &
                       subsys=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set)

      CALL cp_subsys_get(subsys, &
                         virial=virial)

      calc_stress = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (calc_stress .AND. .NOT. calc_forces) CPABORT('Stress cannot be calculated without forces')

      ! Initialize energy and gradient:
      nnp%atomic_energy(:, :) = 0.0_dp
      IF (calc_forces) nnp%myforce(:, :, :) = 0.0_dp
      IF (calc_forces) nnp%committee_forces(:, :, :) = 0.0_dp
      IF (calc_stress) nnp%committee_stress(:, :, :) = 0.0_dp

      !fill coord array
      ig = 1
      DO i = 1, nnp%n_ele
         DO j = 1, nnp%num_atoms
            IF (nnp%ele(i) == particle_set(j)%atomic_kind%element_symbol) THEN
               DO m = 1, 3
                  nnp%coord(m, ig) = particle_set(j)%r(m)
               END DO
               nnp%atoms(ig) = nnp%ele(i)
               CALL get_ptable_info(nnp%atoms(ig), number=nnp%nuc_atoms(ig))
               nnp%ele_ind(ig) = i
               nnp%sort(ig) = j
               nnp%sort_inv(j) = ig
               ig = ig + 1
            END IF
         END DO
      END DO

      ! parallization:
      mecalc = nnp%num_atoms/logger%para_env%num_pe + &
               MIN(MOD(nnp%num_atoms, logger%para_env%num_pe)/ &
                   (logger%para_env%mepos + 1), 1)
      ALLOCATE (allcalc(logger%para_env%num_pe))
      allcalc(:) = 0
      CALL mp_allgather(mecalc, allcalc, logger%para_env%group)
      istart = 1
      DO i = 2, logger%para_env%mepos + 1
         istart = istart + allcalc(i - 1)
      END DO

      ! calc atomic contribution to energy and force
      DO i = istart, istart + mecalc - 1

         ! determine index of atom type and offset
         ind = nnp%ele_ind(i)

         ! reset input nodes of ele(ind):
         nnp%arc(ind)%layer(1)%node(:) = 0.0_dp

         ! compute sym fnct values
         IF (calc_forces) THEN
            !reset input grads of ele(ind):
            nnp%arc(ind)%layer(1)%node_grad(:) = 0.0_dp
            ALLOCATE (dsymdxyz(3, nnp%arc(ind)%n_nodes(1), nnp%num_atoms))
            dsymdxyz(:, :, :) = 0.0_dp
            IF (calc_stress) THEN
               ALLOCATE (stress(3, 3, nnp%arc(ind)%n_nodes(1)))
               stress(:, :, :) = 0.0_dp
               CALL nnp_calc_acsf(nnp, i, dsymdxyz, stress)
            ELSE
               CALL nnp_calc_acsf(nnp, i, dsymdxyz)
            END IF
         ELSE
            CALL nnp_calc_acsf(nnp, i)
         END IF

         DO i_com = 1, nnp%n_committee
            ! predict energy
            CALL nnp_predict(nnp%arc(ind), nnp, i_com)
            nnp%atomic_energy(i, i_com) = nnp%arc(ind)%layer(nnp%n_layer)%node(1) + &
                                          nnp%atom_energies(ind)

            ! predict forces
            IF (calc_forces) THEN
               ALLOCATE (denergydsym(nnp%arc(ind)%n_nodes(1)))
               denergydsym(:) = 0.0_dp
               CALL nnp_gradients(nnp%arc(ind), nnp, i_com, denergydsym)
               DO j = 1, nnp%arc(ind)%n_nodes(1)
                  DO k = 1, nnp%num_atoms
                     DO m = 1, 3
                        nnp%myforce(m, k, i_com) = nnp%myforce(m, k, i_com) - denergydsym(j)*dsymdxyz(m, j, k)
                     END DO
                  END DO
                  IF (calc_stress) THEN
                     nnp%committee_stress(:, :, i_com) = nnp%committee_stress(:, :, i_com) - &
                                                         denergydsym(j)*stress(:, :, j)
                  END IF
               END DO
               DEALLOCATE (denergydsym)
            END IF
         END DO

         !deallocate memory
         IF (calc_forces) THEN
            DEALLOCATE (dsymdxyz)
            IF (calc_stress) THEN
               DEALLOCATE (stress)
            END IF
         END IF

      END DO ! loop over num_atoms

      ! calculate energy:
      CALL mp_sum(nnp%atomic_energy(:, :), logger%para_env%group)
      nnp%committee_energy(:) = SUM(nnp%atomic_energy, 1)
      nnp%nnp_potential_energy = SUM(nnp%committee_energy)/REAL(nnp%n_committee, dp)

      IF (calc_forces) THEN
         ! bring myforce to force array
         DO j = 1, nnp%num_atoms
            DO k = 1, 3
               nnp%committee_forces(k, (nnp%sort(j)), :) = nnp%myforce(k, j, :)
            END DO
         END DO
         CALL mp_sum(nnp%committee_forces, logger%para_env%group)
         nnp%nnp_forces(:, :) = SUM(nnp%committee_forces, DIM=3)/REAL(nnp%n_committee, dp)
         DO j = 1, nnp%num_atoms
            particle_set(j)%f(:) = nnp%nnp_forces(:, j)
         END DO
      END IF

      IF (calc_stress) THEN
         CALL mp_sum(nnp%committee_stress, logger%para_env%group)
         virial%pv_virial = SUM(nnp%committee_stress, DIM=3)/REAL(nnp%n_committee, dp)
      END IF

      ! Bias the standard deviation of committee disagreement
      IF (nnp%bias) THEN
         CALL nnp_bias_sigma(nnp, calc_forces)
         nnp%nnp_potential_energy = nnp%nnp_potential_energy + nnp%bias_energy
         IF (calc_forces) THEN
            DO j = 1, nnp%num_atoms
               particle_set(j)%f(:) = particle_set(j)%f(:) + nnp%bias_forces(:, j)
            END DO
         END IF
         ! print properties if requested
         print_section => section_vals_get_subs_vals(nnp%nnp_input, "BIAS%PRINT")
         CALL nnp_print_bias(nnp, print_section)
      END IF

      ! print properties if requested
      print_section => section_vals_get_subs_vals(nnp%nnp_input, "PRINT")
      CALL nnp_print(nnp, print_section)

      DEALLOCATE (allcalc)

      CALL timestop(handle)

   END SUBROUTINE nnp_calc_energy_force

! **************************************************************************************************
!> \brief Calculate bias potential and force based on standard deviation of committee disagreement
!> \param nnp ...
!> \param calc_forces ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_bias_sigma(nnp, calc_forces)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'nnp_bias_sigma', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: avrg, pref, sigma

      CALL timeset(routineN, handle)

      ! init
      sigma = 0.0_dp
      nnp%bias_energy = 0.0_dp
      IF (calc_forces) nnp%bias_forces = 0.0_dp

      ! Subtract reference energy of each committee member, if requested
      IF (nnp%bias_align) THEN
         ! committee energy not used afterward, therefore overwritten
         nnp%committee_energy(:) = nnp%committee_energy(:) - nnp%bias_e_avrg(:)
      END IF

      ! <E> = 1/n sum(E_i)
      ! sigma = sqrt(1/n sum((E_i - <E>)**2))
      !       = sqrt(1/n sum(dE_i**2))
      avrg = SUM(nnp%committee_energy)/REAL(nnp%n_committee, dp)
      DO i = 1, nnp%n_committee
         sigma = sigma + (nnp%committee_energy(i) - avrg)**2
      END DO
      sigma = SQRT(sigma/REAL(nnp%n_committee, dp))
      nnp%bias_sigma = sigma

      IF (sigma > nnp%bias_sigma0) THEN
         ! E_b = 0.5 * kb * (sigma - sigma_0)**2
         nnp%bias_energy = 0.5_dp*nnp%bias_kb*(sigma - nnp%bias_sigma0)**2

         IF (calc_forces) THEN
            ! nabla(E_b) = kb*(sigma - sigma_0)*nabla(sigma)
            ! nabla(sigma) = 1/sigma * 1/n sum(dE_i* nabla(dE_i))
            ! nabla(dE_i) = nabla(E_i) - nabla(<E>)
            pref = nnp%bias_kb*(1.0_dp - nnp%bias_sigma0/sigma)
            DO i = 1, nnp%n_committee
               nnp%bias_forces(:, :) = nnp%bias_forces(:, :) + &
                                       (nnp%committee_energy(i) - avrg)* &
                                       (nnp%committee_forces(:, :, i) - nnp%nnp_forces(:, :))
            END DO
            pref = pref/REAL(nnp%n_committee, dp)
            nnp%bias_forces(:, :) = nnp%bias_forces(:, :)*pref
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE nnp_bias_sigma

! **************************************************************************************************
!> \brief Print properties according to the requests in input file
!> \param nnp ...
!> \param print_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print(nnp, print_section)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section

      INTEGER                                            :: unit_nr
      LOGICAL                                            :: explicit, file_is_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(print_section, "ENERGIES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".data", &
                                        middle_name="nnp-energies", is_new_file=file_is_new)
         IF (unit_nr > 0) CALL nnp_print_energies(nnp, unit_nr, file_is_new)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "FORCES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (unit_nr > 0) CALL nnp_print_forces(nnp, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "FORCES_SIGMA")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                        middle_name="nnp-forces-std")
         IF (unit_nr > 0) CALL nnp_print_force_sigma(nnp, unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      ! Output structures with extrapolation warning on any processor
      CALL mp_sum(nnp%output_expol, logger%para_env%group)
      IF (nnp%output_expol) THEN
         print_key => section_vals_get_subs_vals(print_section, "EXTRAPOLATION")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                           middle_name="nnp-extrapolation")
            IF (unit_nr > 0) CALL nnp_print_expol(nnp, unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)
         END IF
      END IF

      print_key => section_vals_get_subs_vals(print_section, "SUM_FORCE")

      CALL section_vals_val_get(print_section, "SUM_FORCE%ATOM_LIST", &
                                explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".dat", &
                                           middle_name="nnp-sumforce", is_new_file=file_is_new)
            IF (unit_nr > 0) CALL nnp_print_sumforces(nnp, print_section, unit_nr, file_is_new)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)
         END IF
      END IF

   END SUBROUTINE nnp_print

! **************************************************************************************************
!> \brief Print NNP energies and standard deviation sigma
!> \param nnp ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_energies(nnp, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: std

      IF (file_is_new) THEN
         WRITE (unit_nr, "(A1,1X,A20)", ADVANCE='no') "#", "NNP Average [a.u.],"
         WRITE (unit_nr, "(A20)", ADVANCE='no') "NNP sigma [a.u.]"
         DO i = 1, nnp%n_committee
            WRITE (unit_nr, "(A17,I3)", ADVANCE='no') "NNP", i
         END DO
         WRITE (unit_nr, "(A)") ""
      END IF

      fmt_string = "(2X,2(F20.9))"
      WRITE (fmt_string, "(A,I3,A)") "(2X", nnp%n_committee + 2, "(F20.9))"
      std = SUM((SUM(nnp%atomic_energy, 1) - nnp%nnp_potential_energy)**2)
      std = std/REAL(nnp%n_committee, dp)
      std = SQRT(std)
      WRITE (unit_nr, fmt_string) nnp%nnp_potential_energy, std, SUM(nnp%atomic_energy, 1)

   END SUBROUTINE nnp_print_energies

! **************************************************************************************************
!> \brief Print nnp forces
!> \param nnp ...
!> \param print_key ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_forces(nnp, print_key)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_key

      CHARACTER(len=default_path_length)                 :: fmt_string, middle_name
      INTEGER                                            :: i, j, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      ! TODO use write_particle_coordinates from particle_methods.F
      DO i = 1, nnp%n_committee
         WRITE (fmt_string, *) i
         WRITE (middle_name, "(A,A)") "nnp-forces-", ADJUSTL(TRIM(fmt_string))
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                        middle_name=TRIM(middle_name))
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) nnp%num_atoms
            WRITE (unit_nr, "(A,1X,A,A,F20.9)") "NNP forces [a.u.] of committee member", &
               ADJUSTL(TRIM(fmt_string)), "energy [a.u.]=", nnp%committee_energy(i)

            fmt_string = "(A4,1X,3F20.10)"
            DO j = 1, nnp%num_atoms
               WRITE (unit_nr, fmt_string) nnp%atoms(nnp%sort_inv(j)), nnp%committee_forces(:, j, i)
            END DO
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END DO

   END SUBROUTINE nnp_print_forces

! **************************************************************************************************
!> \brief Print standard deviation sigma of NNP forces
!> \param nnp ...
!> \param unit_nr ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_force_sigma(nnp, unit_nr)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr

      INTEGER                                            :: i, j
      REAL(KIND=dp), DIMENSION(3)                        :: var

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) nnp%num_atoms
         WRITE (unit_nr, "(A,1X,A)") "NNP sigma of forces [a.u.]"

         DO i = 1, nnp%num_atoms
            var = 0.0_dp
            DO j = 1, nnp%n_committee
               var = var + (nnp%committee_forces(:, i, j) - nnp%nnp_forces(:, i))**2
            END DO
            var = var/REAL(nnp%n_committee)
            var = SQRT(var)
            WRITE (unit_nr, "(A4,1X,3F20.10)") nnp%atoms(nnp%sort_inv(i)), var
         END DO
      END IF

   END SUBROUTINE nnp_print_force_sigma

! **************************************************************************************************
!> \brief Print structures with extrapolation warning
!> \param nnp ...
!> \param unit_nr ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_expol(nnp, unit_nr)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: mass, unit_conv
      REAL(KIND=dp), DIMENSION(3)                        :: com

      nnp%expol = nnp%expol + 1
      WRITE (unit_nr, *) nnp%num_atoms
      WRITE (unit_nr, "(A,1X,I6)") "NNP extrapolation point N =", nnp%expol

      ! move to COM of solute and wrap the box
      ! coord not needed afterwards, therefore manipulation ok
      com = 0.0_dp
      mass = 0.0_dp
      DO i = 1, nnp%num_atoms
         CALL get_ptable_info(nnp%atoms(i), amass=unit_conv)
         com(:) = com(:) + nnp%coord(:, i)*unit_conv
         mass = mass + unit_conv
      END DO
      com(:) = com(:)/mass

      DO i = 1, nnp%num_atoms
         nnp%coord(:, i) = nnp%coord(:, i) - com(:)
      END DO

      ! write out coordinates
      unit_conv = cp_unit_from_cp2k(1.0_dp, TRIM("angstrom"))
      fmt_string = "(A4,1X,3F20.10)"
      DO i = 1, nnp%num_atoms
         WRITE (unit_nr, fmt_string) &
            nnp%atoms(nnp%sort_inv(i)), &
            nnp%coord(1, nnp%sort_inv(i))*unit_conv, &
            nnp%coord(2, nnp%sort_inv(i))*unit_conv, &
            nnp%coord(3, nnp%sort_inv(i))*unit_conv
      END DO

   END SUBROUTINE nnp_print_expol

! **************************************************************************************************
!> \brief Print properties number according the requests in input file
!> \param nnp ...
!> \param print_section ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_bias(nnp, print_section)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section

      INTEGER                                            :: unit_nr
      LOGICAL                                            :: file_is_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(print_section, "BIAS_ENERGY")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".data", &
                                        middle_name="nnp-bias-energy", is_new_file=file_is_new)
         IF (unit_nr > 0) CALL nnp_print_bias_energy(nnp, unit_nr, file_is_new)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "BIAS_FORCES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                        middle_name="nnp-bias-forces")
         IF (unit_nr > 0) CALL nnp_print_bias_forces(nnp, unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

   END SUBROUTINE nnp_print_bias

! **************************************************************************************************
!> \brief Print NNP bias energies
!> \param nnp ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_bias_energy(nnp, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i

      IF (file_is_new) THEN
         WRITE (unit_nr, "(A1)", ADVANCE='no') "#"
         WRITE (unit_nr, "(2(2X,A19))", ADVANCE='no') "Sigma [a.u.]", "Bias energy [a.u.]"
         DO i = 1, nnp%n_committee
            IF (nnp%bias_align) THEN
               WRITE (unit_nr, "(2X,A16,I3)", ADVANCE='no') "shifted E_NNP", i
            ELSE
               WRITE (unit_nr, "(2X,A16,I3)", ADVANCE='no') "E_NNP", i
            END IF
         END DO
         WRITE (unit_nr, "(A)") ""

      END IF

      WRITE (fmt_string, "(A,I3,A)") "(2X,", nnp%n_committee + 2, "(F20.9,1X))"
      WRITE (unit_nr, fmt_string) nnp%bias_sigma, nnp%bias_energy, nnp%committee_energy

   END SUBROUTINE nnp_print_bias_energy

! **************************************************************************************************
!> \brief Print NNP bias forces
!> \param nnp ...
!> \param unit_nr ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_bias_forces(nnp, unit_nr)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i

      ! TODO use write_particle_coordinates from particle_methods.F
      WRITE (unit_nr, *) nnp%num_atoms
      WRITE (unit_nr, "(A,F20.9)") "NNP bias forces [a.u.] for bias energy [a.u]=", nnp%bias_energy

      fmt_string = "(A4,1X,3F20.10)"
      DO i = 1, nnp%num_atoms
         WRITE (unit_nr, fmt_string) nnp%atoms(nnp%sort_inv(i)), nnp%bias_forces(:, i)
      END DO

   END SUBROUTINE nnp_print_bias_forces

! **************************************************************************************************
!> \brief Print NNP summed forces
!> \param nnp ...
!> \param print_section ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2020-10-10
!> \author Christoph Schran (christoph.schran@rub.de)
! **************************************************************************************************
   SUBROUTINE nnp_print_sumforces(nnp, print_section, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      CHARACTER(len=default_path_length)                 :: fmt_string
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: atomlist
      INTEGER                                            :: i, ig, j, n
      REAL(KIND=dp), DIMENSION(3)                        :: rvec

      NULLIFY (atomlist)
      IF (file_is_new) THEN
         WRITE (unit_nr, "(A)") "# Summed forces [a.u.]"
      END IF

      rvec = 0.0_dp

      ! get atoms to sum over:
      CALL section_vals_val_get(print_section, "SUM_FORCE%ATOM_LIST", &
                                c_vals=atomlist)
      IF (ASSOCIATED(atomlist)) THEN
         n = SIZE(atomlist)
         DO i = 1, nnp%num_atoms
            DO j = 1, n
               ig = nnp%sort_inv(i)
               IF (TRIM(ADJUSTL(atomlist(j))) == TRIM(ADJUSTL(nnp%atoms(ig)))) THEN
                  rvec(:) = rvec(:) + nnp%nnp_forces(:, i)
               END IF
            END DO
         END DO
      END IF

      fmt_string = "(3(F20.10,1X))"
      WRITE (unit_nr, fmt_string) rvec

   END SUBROUTINE nnp_print_sumforces

END MODULE nnp_force
