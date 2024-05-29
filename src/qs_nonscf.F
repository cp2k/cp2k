!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for Quickstep NON-SCF run.
!> \par History
!>      - initial setup [JGH, 2024]
!> \author JGH (13.05.2024)
! **************************************************************************************************
MODULE qs_nonscf
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_copy,&
                                              dbcsr_dot,&
                                              dbcsr_p_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf,                          ONLY: init_scf_loop
   USE qs_scf_initialization,           ONLY: qs_scf_env_initialize
   USE qs_scf_loop_utils,               ONLY: qs_scf_new_mos,&
                                              qs_scf_new_mos_kp
   USE qs_scf_output,                   ONLY: qs_scf_loop_print,&
                                              qs_scf_print_summary,&
                                              qs_scf_write_mos
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_wf_history_methods,           ONLY: wfi_update
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_nonscf'

   PUBLIC :: nonscf

CONTAINS

! **************************************************************************************************
!> \brief Find solution to HC=SCE
!> \param qs_env the qs_environment where to perform the scf procedure
!> \par History
!>      none
!> \author JGH
!> \note
! **************************************************************************************************
   SUBROUTINE nonscf(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL get_qs_env(qs_env, scf_env=scf_env, scf_control=scf_control, &
                      dft_control=dft_control)

      IF (dft_control%qs_control%do_ls_scf) THEN
         ! Density matrix based solver

         CPABORT("NOT AVAILABLE")
         !CALL ls_scf(qs_env=qs_env)

      ELSE
         ! Wavefunction based solver

         IF (.NOT. ASSOCIATED(scf_env)) THEN
            CALL qs_scf_env_initialize(qs_env, scf_env)
            CALL set_qs_env(qs_env, scf_env=scf_env)
         ELSE
            CALL qs_scf_env_initialize(qs_env, scf_env)
         END IF

         CALL do_nonscf(qs_env, scf_env, scf_control)

         ! add the converged wavefunction to the wavefunction history
         IF (ASSOCIATED(qs_env%wf_history)) THEN
            CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
         END IF

         ! compute properties that depend on the wavefunction
         CALL qs_scf_compute_properties(qs_env)

      END IF

   END SUBROUTINE nonscf

! **************************************************************************************************
!> \brief Solve KS equation for fixed potential
!> \param qs_env ...
!> \param scf_env the scf_env where to perform the scf procedure
!> \param scf_control ...
!> \par History
!>      none
!> \author JGH
!> \note
! **************************************************************************************************
   SUBROUTINE do_nonscf(qs_env, scf_env, scf_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'do_nonscf'

      INTEGER                                            :: handle, img, ispin, output_unit
      LOGICAL                                            :: diis_step, do_kpoints
      REAL(KIND=dp)                                      :: pc_ener, qmmm_el, t1, t2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_ks, rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env=qs_env, &
                      energy=energy, &
                      ks_env=ks_env, &
                      rho=rho, &
                      mos=mos, &
                      input=input, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      para_env=para_env)

      DO ispin = 1, dft_control%nspins
         CPASSERT(.NOT. mos(ispin)%use_mo_coeff_b)
      END DO

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      CALL init_scf_loop(scf_env=scf_env, qs_env=qs_env, scf_section=scf_section)

      ! Calculate KS matrix
      CALL qs_ks_update_qs_env(qs_env, just_energy=.FALSE., calculate_forces=.FALSE.)

      ! print 'heavy weight' or relatively expensive quantities
      CALL qs_scf_loop_print(qs_env, scf_env, para_env)

      ! Diagonalization
      IF (do_kpoints) THEN
         ! kpoints
         CALL qs_scf_new_mos_kp(qs_env, scf_env, scf_control, diis_step)
      ELSE
         ! Gamma points only
         CALL qs_scf_new_mos(qs_env, scf_env, scf_control, scf_section, diis_step, .FALSE.)
      END IF

      ! Print requested MO information (can be computationally expensive with OT)
      CALL qs_scf_write_mos(qs_env, scf_env, final_mos=.TRUE.)

      ! copy density matrix
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      DO ispin = 1, dft_control%nspins
         DO img = 1, SIZE(rho_ao_kp, 2)
            CALL dbcsr_copy(rho_ao_kp(ispin, img)%matrix, scf_env%p_mix_new(ispin, img)%matrix)
         END DO
      END DO
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE., potential_changed=.TRUE.)

      ! core energy : Tr(PH)
      CALL get_qs_env(qs_env, matrix_ks_kp=matrixkp_ks)
      IF (qs_env%qmmm) THEN
         ! Compute QM/MM Energy
         CPASSERT(SIZE(matrixkp_ks, 2) == 1)
         DO ispin = 1, dft_control%nspins
            CALL dbcsr_dot(qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           matrixkp_ks(ispin, 1)%matrix, qmmm_el)
            energy%qmmm_el = energy%qmmm_el + qmmm_el
         END DO
         pc_ener = qs_env%ks_qmmm_env%pc_ener
         energy%qmmm_el = energy%qmmm_el + pc_ener
      ELSE
         energy%qmmm_el = 0.0_dp
      END IF
      energy%total = energy%total - energy%core
      CALL calculate_ptrace(matrixkp_ks, rho_ao_kp, energy%core, dft_control%nspins)
      energy%total = energy%total + energy%core + energy%qmmm_el

      t2 = m_walltime()

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, &
                FMT="(T2,A,T40,A,F10.2,T61,F20.10)") &
            "Diagonalization", "Time:", t2 - t1, energy%total
      END IF

      CALL qs_scf_print_summary(output_unit, qs_env)

      CALL timestop(handle)

   END SUBROUTINE do_nonscf

END MODULE qs_nonscf
