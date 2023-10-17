!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods dealing with Neural Network interaction potential
!> \author Laura Duran
!> \date   2023-02-17
! **************************************************************************************************
MODULE helium_nnp

   USE bibliography,                    ONLY: Behler2007,&
                                              Behler2011,&
                                              Schran2020a,&
                                              Schran2020b,&
                                              cite_reference
   USE cell_methods,                    ONLY: cell_create
   USE cell_types,                      ONLY: cell_release,&
                                              cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE helium_types,                    ONLY: helium_solvent_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE nnp_environment,                 ONLY: nnp_init_model
   USE nnp_environment_types,           ONLY: nnp_env_get,&
                                              nnp_env_set,&
                                              nnp_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE physcon,                         ONLY: angstrom
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_nnp'

   PUBLIC :: helium_init_nnp, &
             helium_nnp_print

CONTAINS

! ***************************************************************************
!> \brief  Read and initialize all the information for neural network potentials
!> \param helium ...
!> \param nnp ...
!> \param input ...
!> \date   2023-02-21
!> \author lduran
! **************************************************************************************************
   SUBROUTINE helium_init_nnp(helium, nnp, input)
      TYPE(helium_solvent_type), INTENT(INOUT)           :: helium
      TYPE(nnp_type), POINTER                            :: nnp
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(len=default_path_length)                 :: msg_str
      CHARACTER(len=default_string_length)               :: elem
      INTEGER                                            :: i, ig, is, j
      INTEGER, DIMENSION(3)                              :: periodicity
      LOGICAL                                            :: found
      TYPE(cell_type), POINTER                           :: he_cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: sr_cutoff_section

      CALL cite_reference(Behler2007)
      CALL cite_reference(Behler2011)
      CALL cite_reference(Schran2020a)
      CALL cite_reference(Schran2020b)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      CALL nnp_env_set(nnp_env=nnp, nnp_input=input)

      nnp%num_atoms = helium%solute_atoms + 1

      CALL nnp_init_model(nnp, "HELIUM NNP")

      periodicity = 0
      IF (helium%periodic) periodicity = 1
      NULLIFY (he_cell)
      CALL cell_create(he_cell, hmat=helium%cell_m, &
                       periodic=periodicity, tag="HELIUM NNP")
      CALL nnp_env_set(nnp, cell=he_cell)
      CALL cell_release(he_cell)

      ! Set up arrays for calculation:
      ALLOCATE (nnp%ele_ind(nnp%num_atoms))
      ALLOCATE (nnp%nuc_atoms(nnp%num_atoms))
      ALLOCATE (nnp%coord(3, nnp%num_atoms))
      ALLOCATE (nnp%nnp_forces(3, nnp%num_atoms))
      ALLOCATE (nnp%atoms(nnp%num_atoms))
      ALLOCATE (nnp%sort(nnp%num_atoms - 1))

      !fill arrays, assume that order will not change during simulation
      ig = 1
      is = 1
      DO i = 1, nnp%n_ele
         IF (nnp%ele(i) == 'He') THEN
            nnp%atoms(ig) = 'He'
            CALL get_ptable_info(nnp%atoms(ig), number=nnp%nuc_atoms(ig))
            nnp%ele_ind(ig) = i
            ig = ig + 1
         END IF
         DO j = 1, helium%solute_atoms
            IF (nnp%ele(i) == helium%solute_element(j)) THEN
               nnp%atoms(ig) = nnp%ele(i)
               CALL get_ptable_info(nnp%atoms(ig), number=nnp%nuc_atoms(ig))
               nnp%ele_ind(ig) = i
               nnp%sort(is) = j
               ig = ig + 1
               is = is + 1
            END IF
         END DO
      END DO

      ALLOCATE (helium%nnp_sr_cut(nnp%n_ele))
      helium%nnp_sr_cut = 0.0_dp

      sr_cutoff_section => section_vals_get_subs_vals(nnp%nnp_input, "SR_CUTOFF")
      CALL section_vals_get(sr_cutoff_section, n_repetition=is)
      DO i = 1, is
         CALL section_vals_val_get(sr_cutoff_section, "ELEMENT", c_val=elem, i_rep_section=i)
         found = .FALSE.
         DO ig = 1, nnp%n_ele
            IF (TRIM(nnp%ele(ig)) == TRIM(elem)) THEN
               found = .TRUE.
               CALL section_vals_val_get(sr_cutoff_section, "RADIUS", r_val=helium%nnp_sr_cut(ig), &
                                         i_rep_section=i)
            END IF
         END DO
         IF (.NOT. found) THEN
            msg_str = "SR_CUTOFF for element "//TRIM(elem)//" defined but not found in NNP"
            CPWARN(msg_str)
         END IF
      END DO
      helium%nnp_sr_cut(:) = helium%nnp_sr_cut(:)**2

      RETURN

   END SUBROUTINE helium_init_nnp

! **************************************************************************************************
!> \brief Print properties according to the requests in input file
!> \param nnp ...
!> \param print_section ...
!> \param ind_he ...
!> \date   2023-07-31
!> \author Laura Duran
! **************************************************************************************************
   SUBROUTINE helium_nnp_print(nnp, print_section, ind_he)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      TYPE(section_vals_type), INTENT(IN), POINTER       :: print_section
      INTEGER, INTENT(IN)                                :: ind_he

      INTEGER                                            :: unit_nr
      LOGICAL                                            :: file_is_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, print_key)
      logger => cp_get_default_logger()

      print_key => section_vals_get_subs_vals(print_section, "ENERGIES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".data", &
                                        middle_name="helium-nnp-energies", is_new_file=file_is_new)
         IF (unit_nr > 0) CALL helium_nnp_print_energies(nnp, unit_nr, file_is_new)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      print_key => section_vals_get_subs_vals(print_section, "FORCES_SIGMA")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".data", &
                                        middle_name="helium-nnp-forces-std", is_new_file=file_is_new)
         IF (unit_nr > 0) CALL helium_nnp_print_force_sigma(nnp, unit_nr, file_is_new)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      CALL logger%para_env%sum(nnp%output_expol)
      IF (nnp%output_expol) THEN
         print_key => section_vals_get_subs_vals(print_section, "EXTRAPOLATION")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            unit_nr = cp_print_key_unit_nr(logger, print_key, extension=".xyz", &
                                           middle_name="-NNP-He-extrapolation")
            IF (unit_nr > 0) CALL helium_nnp_print_expol(nnp, unit_nr, ind_he)
            CALL cp_print_key_finished_output(unit_nr, logger, print_key)
         END IF
      END IF

   END SUBROUTINE helium_nnp_print

! **************************************************************************************************
!> \brief Print NNP energies and standard deviation sigma
!> \param nnp ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2023-07-31
!> \author Laura Duran
! **************************************************************************************************
   SUBROUTINE helium_nnp_print_energies(nnp, unit_nr, file_is_new)
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

   END SUBROUTINE helium_nnp_print_energies

! **************************************************************************************************
!> \brief Print standard deviation sigma of NNP forces
!> \param nnp ...
!> \param unit_nr ...
!> \param file_is_new ...
!> \date   2023-07-31
!> \author Laura Duran
! **************************************************************************************************
   SUBROUTINE helium_nnp_print_force_sigma(nnp, unit_nr, file_is_new)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: file_is_new

      INTEGER                                            :: i, ig, j
      REAL(KIND=dp), DIMENSION(3)                        :: var

      IF (unit_nr > 0) THEN
         IF (file_is_new) THEN
            WRITE (unit_nr, "(A,1X,A)") "#   NNP sigma of forces [a.u.]    x, y, z coordinates"
         END IF

         ig = 1
         DO i = 1, nnp%num_atoms
            IF (nnp%ele(i) == 'He') THEN
               var = 0.0_dp
               DO j = 1, nnp%n_committee
                  var = var + (nnp%committee_forces(:, i, j) - nnp%nnp_forces(:, i))**2
               END DO
               WRITE (unit_nr, "(A4,1X,3E20.10)") nnp%ele(i), var
            END IF
            ig = ig + 1
         END DO
      END IF

   END SUBROUTINE helium_nnp_print_force_sigma

! **************************************************************************************************
!> \brief Print structures with extrapolation warning
!> \param nnp ...
!> \param unit_nr ...
!> \param ind_he ...
!> \date   2023-10-11
!> \author Harald Forbert (harald.forbert@rub.de)
! **************************************************************************************************
   SUBROUTINE helium_nnp_print_expol(nnp, unit_nr, ind_he)
      TYPE(nnp_type), INTENT(INOUT)                      :: nnp
      INTEGER, INTENT(IN)                                :: unit_nr, ind_he

      CHARACTER(len=default_path_length)                 :: fmt_string
      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: mass, unit_conv
      REAL(KIND=dp), DIMENSION(3)                        :: com
      TYPE(cell_type), POINTER                           :: cell

      NULLIFY (cell)
      CALL nnp_env_get(nnp_env=nnp, cell=cell)

      nnp%expol = nnp%expol + 1
      WRITE (unit_nr, *) nnp%num_atoms
      WRITE (unit_nr, "(A,1X,I6)") "HELIUM-NNP extrapolation point N =", nnp%expol

      ! move to COM of solute and wrap the box
      ! coord not needed afterwards, therefore manipulation ok
      com = 0.0_dp
      mass = 0.0_dp
      DO i = 1, nnp%num_atoms
         IF (i == ind_he) CYCLE
         CALL get_ptable_info(nnp%atoms(i), amass=unit_conv)
         com(:) = com(:) + nnp%coord(:, i)*unit_conv
         mass = mass + unit_conv
      END DO
      com(:) = com(:)/mass

      DO i = 1, nnp%num_atoms
         nnp%coord(:, i) = nnp%coord(:, i) - com(:)
         nnp%coord(:, i) = pbc(nnp%coord(:, i), cell)
      END DO

      unit_conv = cp_unit_from_cp2k(1.0_dp, TRIM("angstrom"))
      fmt_string = "(A4,1X,3F20.10)"
      DO i = 1, nnp%num_atoms
         WRITE (unit_nr, fmt_string) &
            nnp%atoms(i), &
            nnp%coord(1, i)*unit_conv, &
            nnp%coord(2, i)*unit_conv, &
            nnp%coord(3, i)*unit_conv
      END DO

   END SUBROUTINE helium_nnp_print_expol

END MODULE helium_nnp
