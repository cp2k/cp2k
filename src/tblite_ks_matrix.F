!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief tblite matrix build
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_ks_matrix

   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_dot,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE tblite_interface,                ONLY: tb_get_energy,&
                                              tb_ham_add_coulomb,&
                                              tb_update_charges
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_ks_matrix'

   PUBLIC :: build_tblite_ks_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param ext_ks_matrix ...
! **************************************************************************************************
   SUBROUTINE build_tblite_ks_matrix(qs_env, calculate_forces, just_energy, ext_ks_matrix)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: ext_ks_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'build_tblite_ks_matrix'

      INTEGER                                            :: handle, img, iounit, ispin, nimg, ns, &
                                                            nspins
      LOGICAL                                            :: do_efield
      REAL(KIND=dp)                                      :: pc_ener, qmmm_el
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: charges
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p1, mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_h
      TYPE(dbcsr_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: scf_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, logger, scf_section, ks_env, &
               ks_matrix, rho, energy)
      CPASSERT(ASSOCIATED(qs_env))

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      ks_env=ks_env, &
                      matrix_ks_kp=ks_matrix, &
                      rho=rho, &
                      energy=energy)

      IF (PRESENT(ext_ks_matrix)) THEN
         ! remap pointer to allow for non-kpoint external ks matrix
         ! ext_ks_matrix is used in linear response code
         ns = SIZE(ext_ks_matrix)
         ks_matrix(1:ns, 1:1) => ext_ks_matrix(1:ns)
      END IF

      energy%qmmm_el = 0.0_dp

      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")
      nspins = dft_control%nspins
      nimg = dft_control%nimages
      CPASSERT(ASSOCIATED(matrix_h))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(SIZE(ks_matrix) > 0)

      !WRITE (*, *) "::::::::::::::::: building GFN1 matrix, ks_matrix, in:"
      DO ispin = 1, nspins
         DO img = 1, nimg
            ! copy the core matrix into the fock matrix
            CALL dbcsr_copy(ks_matrix(ispin, img)%matrix, matrix_h(1, img)%matrix)
            !CALL dbcsr_print(ks_matrix(ispin, img)%matrix)
         END DO
      END DO

      IF (dft_control%apply_period_efield .OR. dft_control%apply_efield .OR. &
          dft_control%apply_efield_field) THEN
         do_efield = .TRUE.
         CPABORT("Not implemented yet. Use CP2K routines for GFN1")
      ELSE
         do_efield = .FALSE.
      END IF

      CALL tb_update_charges(qs_env, para_env, dft_control, qs_env%tb_tblite, charges, calculate_forces)

      CALL tb_ham_add_coulomb(qs_env, qs_env%tb_tblite, ks_matrix, dft_control, charges, calculate_forces)

      IF (qs_env%qmmm) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         DO ispin = 1, nspins
            ! If QM/MM sumup the 1el Hamiltonian
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           1.0_dp, 1.0_dp)
            CALL qs_rho_get(rho, rho_ao=matrix_p1)
            ! Compute QM/MM Energy
            CALL dbcsr_dot(qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           matrix_p1(ispin)%matrix, qmmm_el)
            energy%qmmm_el = energy%qmmm_el + qmmm_el
         END DO
         pc_ener = qs_env%ks_qmmm_env%pc_ener
         energy%qmmm_el = energy%qmmm_el + pc_ener
      END IF

      CALL tb_get_energy(qs_env%tb_tblite, energy)

      energy%total = energy%core + energy%hartree + energy%efield + energy%qmmm_el + &
                     energy%repulsive + energy%dispersion + energy%dftb3 + energy%kTS

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Charge fluctuation energy:                     ", energy%hartree, &
            "London dispersion energy:                      ", energy%dispersion
         IF (dft_control%qs_control%xtb_control%xb_interaction) &
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for halogen bonding:                ", energy%xtb_xb_inter
         IF (dft_control%qs_control%xtb_control%do_nonbonded) &
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for nonbonded interactions:         ", energy%xtb_nonbonded
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "DFTB3 3rd Order Energy Correction              ", energy%dftb3
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")
      ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         BLOCK
            TYPE(mo_set_type), DIMENSION(:), POINTER         :: mo_array
            CALL get_qs_env(qs_env, mo_derivs=mo_derivs, mos=mo_array)
            DO ispin = 1, SIZE(mo_derivs)
               CALL get_mo_set(mo_set=mo_array(ispin), mo_coeff_b=mo_coeff)
               IF (.NOT. mo_array(ispin)%use_mo_coeff_b) THEN
                  CPABORT("")
               END IF
               CALL dbcsr_multiply('n', 'n', 1.0_dp, ks_matrix(ispin, 1)%matrix, mo_coeff, &
                                   0.0_dp, mo_derivs(ispin)%matrix)
            END DO
         END BLOCK
      END IF

      !WRITE (*, *) "::::::::::::::::: built GFN1 matrix, ks_matrix, out:"
      !DO ispin = 1, nspins
      !   CALL dbcsr_print(ks_matrix(ispin, 1)%matrix)
      !END DO

      CALL timestop(handle)

   END SUBROUTINE build_tblite_ks_matrix

END MODULE tblite_ks_matrix
