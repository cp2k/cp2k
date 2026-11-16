!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Fractional occupation number weighted density (Grimme and Hansen).
!> \author Falk William Seidel (initial implementation)
! **************************************************************************************************
MODULE qs_fod
   USE bibliography,                    ONLY: Grimme2015FOD,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_constants,                 ONLY: smear_fermi_dirac
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              mo_set_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fod'
   PUBLIC :: qs_scf_post_fod, qs_fod_validate, fod_weights

CONTAINS

! **************************************************************************************************
!> \brief Hole weights below the chemical potential, particle weights above it.
!> \param occupation Fermi-Dirac occupations (0 to maxocc).
!> \param eigenvalues Orbital energies in the same units as mu.
!> \param mu Chemical potential for this spin channel.
!> \param maxocc Maximum occupation, two for RKS and one for UKS.
!> \return Nonnegative FOD weights; both branches agree at half occupation.
! **************************************************************************************************
   PURE FUNCTION fod_weights(occupation, eigenvalues, mu, maxocc) RESULT(weights)
      REAL(dp), DIMENSION(:), INTENT(IN)                 :: occupation, eigenvalues
      REAL(dp), INTENT(IN)                               :: mu, maxocc
      REAL(dp), DIMENSION(SIZE(occupation))              :: weights

      WHERE (eigenvalues <= mu)
         weights = maxocc - occupation
      ELSE WHERE
         weights = occupation
      END WHERE
   END FUNCTION fod_weights

! **************************************************************************************************
!> \brief Reject unsupported FOD settings before starting an expensive SCF calculation.
!> \param input FORCE_EVAL input section.
!> \param logger Output logger.
!> \param qs_env Quickstep environment.
! **************************************************************************************************
   SUBROUTINE qs_fod_validate(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: fod_section

      fod_section => section_vals_get_subs_vals(input, "DFT%PRINT%FOD")
      IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, fod_section, ""), cp_p_file)) RETURN
      CALL get_qs_env(qs_env, dft_control=dft_control, scf_control=scf_control, do_kpoints=do_kpoints)
      IF (do_kpoints) CPABORT("FOD currently supports Gamma-point calculations only.")
      IF (.NOT. dft_control%qs_control%gpw) CPABORT("FOD currently requires the GPW method.")
      IF (qs_env%run_rtp) CPABORT("FOD is not supported for real-time propagation.")
      IF (.NOT. scf_control%smear%do_smear .OR. scf_control%smear%method /= smear_fermi_dirac) THEN
         CPABORT("FOD requires SCF%SMEAR with METHOD FERMI_DIRAC.")
      END IF
      IF (scf_control%smear%electronic_temperature <= 0.0_dp) THEN
         CPABORT("FOD requires a positive electronic temperature.")
      END IF
   END SUBROUTINE qs_fod_validate

! **************************************************************************************************
!> \brief Print the orbital FOD sum and, optionally, its real-space density.
!> \param input FORCE_EVAL input section.
!> \param logger Output logger.
!> \param qs_env Quickstep environment; the SCF orbitals and density are not modified.
!> \param output_unit Main output unit, nonpositive on non-writing ranks.
! **************************************************************************************************
   SUBROUTINE qs_scf_post_fod(input, logger, qs_env, output_unit)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_scf_post_fod'

      INTEGER                                            :: handle, iatom, ikind, ispin, unit_nr
      LOGICAL                                            :: write_cube
      REAL(dp)                                           :: grid_integral, nfod, spin_integral
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: weights, zcharge
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_type), POINTER                          :: matrix_fod
      TYPE(mo_set_type)                                  :: mo_fod
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_c1d_gs_type)                               :: tmp_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: fod_r, tmp_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: fod_section

      fod_section => section_vals_get_subs_vals(input, "DFT%PRINT%FOD")
      IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, fod_section, ""), cp_p_file)) RETURN
      CALL timeset(routineN, handle)
      CALL qs_fod_validate(input, logger, qs_env)
      CALL get_qs_env(qs_env, mos=mos, scf_control=scf_control)
      CPASSERT(ASSOCIATED(mos))
      CALL cite_reference(Grimme2015FOD)
      CALL section_vals_val_get(fod_section, "CUBE", l_val=write_cube)

      NULLIFY (rho, rho_ao, ks_env, pw_env, pool, subsys, particles, matrix_fod)
      IF (write_cube) THEN
         CALL get_qs_env(qs_env, rho=rho, ks_env=ks_env, pw_env=pw_env, subsys=subsys, qs_kind_set=qs_kind_set)
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         CALL qs_subsys_get(subsys, particles=particles)
         ALLOCATE (zcharge(particles%n_els))
         DO iatom = 1, particles%n_els
            ikind = particles%els(iatom)%atomic_kind%kind_number
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zcharge(iatom))
         END DO
         CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
         CALL pool%create_pw(fod_r)
         CALL pool%create_pw(tmp_r)
         CALL pool%create_pw(tmp_g)
         fod_r%array = 0.0_dp
      END IF

      nfod = 0.0_dp
      grid_integral = 0.0_dp
      DO ispin = 1, SIZE(mos)
         CPASSERT(mos(ispin)%nmo > 0)
         IF (mos(ispin)%occupation_numbers(mos(ispin)%nmo) > &
             mos(ispin)%maxocc*scf_control%smear%eps_fermi_dirac) THEN
            CPWARN("FOD: the highest available orbital is occupied; increase ADDED_MOS and check convergence.")
         END IF
         ALLOCATE (weights(mos(ispin)%nmo))
         weights(:) = fod_weights(mos(ispin)%occupation_numbers, mos(ispin)%eigenvalues, &
                                  mos(ispin)%mu, mos(ispin)%maxocc)
         nfod = nfod + SUM(weights)
         IF (write_cube) THEN
            CALL duplicate_mo_set(mo_fod, mos(ispin))
            mo_fod%occupation_numbers = weights
            mo_fod%uniform_occupation = .FALSE.
            ! Include every supplied orbital, including the fractional virtual tail.
            mo_fod%homo = mo_fod%nmo
            ALLOCATE (matrix_fod)
            CALL dbcsr_copy(matrix_fod, rho_ao(ispin)%matrix)
            CALL calculate_density_matrix(mo_fod, matrix_fod)
            CALL calculate_rho_elec(matrix_p=matrix_fod, rho=tmp_r, rho_gspace=tmp_g, &
                                    total_rho=spin_integral, ks_env=ks_env)
            fod_r%array = fod_r%array + tmp_r%array
            ! calculate_rho_elec returns electronic charge, not electron count.
            grid_integral = grid_integral - spin_integral
            CALL dbcsr_release(matrix_fod)
            DEALLOCATE (matrix_fod)
            CALL deallocate_mo_set(mo_fod)
         END IF
         DEALLOCATE (weights)
      END DO

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T2,A,T61,F20.10)') "FOD| N_FOD (orbital sum)", nfod
         IF (write_cube) WRITE (output_unit, '(T2,A,T61,F20.10)') "FOD| Grid integral", grid_integral
      END IF
      IF (write_cube) THEN
         unit_nr = cp_print_key_unit_nr(logger, fod_section, "", extension=".cube", &
                                        middle_name="FOD", file_position="REWIND", log_filename=.FALSE.)
         ! Collective call: non-writing ranks must participate even with unit_nr=-1.
         CALL cp_pw_to_cube(fod_r, unit_nr, "FRACTIONAL OCCUPATION DENSITY [e/bohr^3]", &
                            particles=particles, zeff=zcharge, stride=section_get_ivals(fod_section, "STRIDE"))
         CALL cp_print_key_finished_output(unit_nr, logger, fod_section, "")
         CALL pool%give_back_pw(fod_r)
         CALL pool%give_back_pw(tmp_r)
         CALL pool%give_back_pw(tmp_g)
         DEALLOCATE (zcharge)
      END IF
      CALL timestop(handle)
   END SUBROUTINE qs_scf_post_fod
END MODULE qs_fod
