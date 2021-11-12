!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate EXX in RPA
!> \par History
!>      07.2020 separated from mp2.F [F. Stein, code by Jan Wilhelm]
!> \author Jan Wilhelm, Frederick Stein
! **************************************************************************************************
MODULE rpa_hfx
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_set
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_ri,                          ONLY: hfx_ri_update_ks
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_hfx'

   PUBLIC :: calculate_exx

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param do_gw ...
!> \param do_admm ...
!> \param E_ex_from_GW ...
!> \param t3 ...
! **************************************************************************************************
   SUBROUTINE calculate_exx(qs_env, unit_nr, do_gw, do_admm, E_ex_from_GW, t3)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: do_gw, do_admm
      REAL(KIND=dp), INTENT(IN)                          :: E_ex_from_GW, t3

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_exx', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, irep, n_rep_hf, ns
      REAL(KIND=dp)                                      :: eh1, ehfx, t1, t2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_2d, rho_ao_2d
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      NULLIFY (hfx_sections, input, para_env, matrix_ks, matrix_ks_aux_fit, rho, rho_ao)

      ! if ADMM for exchange energy, get rho_aux_fit and matrix_ks_aux_fit
      IF (do_admm) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         input=input, &
                         para_env=para_env, &
                         energy=energy, &
                         rho_aux_fit=rho, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_aux_fit=matrix_ks_aux_fit)
      ELSE
         CALL get_qs_env(qs_env=qs_env, &
                         input=input, &
                         para_env=para_env, &
                         energy=energy, &
                         rho=rho, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_aux_fit=matrix_ks_aux_fit)
      END IF

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

      ! put matrix_ks to zero
      DO i = 1, SIZE(matrix_ks)
         CALL dbcsr_set(matrix_ks(i)%matrix, 0.0_dp)
         IF (do_admm) THEN
            CALL dbcsr_set(matrix_ks_aux_fit(i)%matrix, 0.0_dp)
         END IF
      END DO

      ! take the exact exchange energy from GW or calculate it
      IF (do_gw) THEN

         IF (qs_env%mp2_env%ri_g0w0%update_xc_energy) THEN
            CALL remove_exc_energy(energy)
            energy%total = energy%total + E_ex_from_GW
            energy%ex = E_ex_from_GW
            t2 = m_walltime()

            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total EXX Time=', t2 - t1 + t3
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'EXX energy  =   ', energy%ex
         END IF

      ELSE

         CALL remove_exc_energy(energy)

         ehfx = 0.0_dp
         DO irep = 1, n_rep_hf
            ns = SIZE(rho_ao)
            rho_ao_2d(1:ns, 1:1) => rho_ao(1:ns)
            ns = SIZE(matrix_ks)
            IF (do_admm) THEN
               matrix_ks_2d(1:ns, 1:1) => matrix_ks_aux_fit(1:ns)
            ELSE
               matrix_ks_2d(1:ns, 1:1) => matrix_ks(1:ns)
            END IF

            IF (qs_env%mp2_env%ri_rpa%x_data(irep, 1)%do_hfx_ri) THEN
               CALL hfx_ri_update_ks(qs_env, qs_env%mp2_env%ri_rpa%x_data(irep, 1)%ri_data, matrix_ks_2d, ehfx, &
                                     rho_ao=rho_ao_2d, geometry_did_change=.TRUE., nspins=ns, &
                                     hf_fraction=qs_env%mp2_env%ri_rpa%x_data(irep, 1)%general_parameter%fraction)
            ELSE

               CALL integrate_four_center(qs_env, qs_env%mp2_env%ri_rpa%x_data, matrix_ks_2d, eh1, &
                                          rho_ao_2d, hfx_sections, &
                                          para_env, .TRUE., irep, .TRUE., &
                                          ispin=1)
               ehfx = ehfx + eh1
            END IF
         END DO

         ! include the EXX contribution to the total energy
         energy%ex = ehfx
         energy%total = energy%total + energy%ex

         t2 = m_walltime()

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total EXX Time=', t2 - t1 + t3
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'EXX energy  =   ', energy%ex

      END IF

      ! reset to zero the Hartree-Fock energy
      energy%ex = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE calculate_exx

! **************************************************************************************************
!> \brief ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE remove_exc_energy(energy)
      TYPE(qs_energy_type), POINTER                      :: energy

      CHARACTER(len=*), PARAMETER :: routineN = 'remove_exc_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Remove the Exchange-correlation energy contributions from the total energy
      energy%total = energy%total - (energy%exc + energy%exc1 + energy%ex + &
                                     energy%exc_aux_fit + energy%exc1_aux_fit)

      energy%exc = 0.0_dp
      energy%exc1 = 0.0_dp
      energy%exc_aux_fit = 0.0_dp
      energy%exc1_aux_fit = 0.0_dp
      energy%ex = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE

END MODULE rpa_hfx

