!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RPA correlation energy and exact exchange energy with the real-space resolution of the
!>        identity (RI-RS) and k-point sampling, input section WF_CORRELATION%RI_RS_RPA
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE rpa_ri_rs
   USE gw_ri_rs_small_cell_full_kp_types,ONLY: gwkp_type
   USE gw_ri_rs_small_cell_full_kp_utils,ONLY: compute_W0,&
                                              compute_chi,&
                                              compute_exchange_energy_ri_rs,&
                                              prepare_ri_rs,&
                                              release_gw_work
   USE gw_utils,                        ONLY: create_and_init_bs_env_for_gw,&
                                              de_init_bs_env
   USE hfx_exx,                         ONLY: remove_exc_energy
   USE input_constants,                 ONLY: ri_rs_small_cell_full_kp
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE post_scf_bandstructure_types,    ONLY: bs_env_release,&
                                              post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: create_and_init_bs_env
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_ri_rs'

   PUBLIC :: run_ri_rs_rpa

CONTAINS

! **************************************************************************************************
!> \brief RPA total energy E = E^DFT - E^xc + E^EXX + E_c^RPA from the converged SCF orbitals:
!>        replaces E^xc by E^EXX in the energies of qs_env and returns E_c^RPA
!> \param qs_env Quickstep environment containing the converged SCF state.
!> \param e_c_rpa RPA correlation energy per unit cell
! **************************************************************************************************
   SUBROUTINE run_ri_rs_rpa(qs_env, e_c_rpa)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: e_c_rpa

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'run_ri_rs_rpa'

      INTEGER                                            :: handle, u
      REAL(KIND=dp)                                      :: e_exx
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: bs_sec, input, rpa_sec

      CALL timeset(routineN, handle)

      NULLIFY (energy, input)
      CALL get_qs_env(qs_env, energy=energy, input=input)
      bs_sec => section_vals_get_subs_vals(input, "PROPERTIES%BANDSTRUCTURE")
      rpa_sec => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RS_RPA")

      CALL create_and_init_bs_env(qs_env, qs_env%bs_env, bs_sec, rpa_sec)
      CALL create_and_init_bs_env_for_gw(qs_env, qs_env%bs_env, bs_sec, rpa_sec)

      SELECT CASE (qs_env%bs_env%gw_implementation)
      CASE (ri_rs_small_cell_full_kp)
         CALL rpa_calc_ri_rs_small_cell_full_kp(qs_env, qs_env%bs_env, e_exx, e_c_rpa)
      CASE DEFAULT
         CPABORT("RI_RS_RPA is implemented only with DFT%KPOINTS")
      END SELECT

      CALL remove_exc_energy(energy)
      energy%ex = e_exx
      energy%total = energy%total + e_exx

      u = qs_env%bs_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,T55,F26.10)') 'RI-RS-RPA| Exact exchange energy E^EXX:', e_exx
         WRITE (u, '(T2,A,T55,F26.10)') 'RI-RS-RPA| RPA correlation energy E_c^RPA:', e_c_rpa
         WRITE (u, '(T2,A,T55,F26.10)') 'RI-RS-RPA| Total energy E^DFT-E^xc+E^EXX+E_c^RPA:', &
            energy%total + e_c_rpa
         WRITE (u, '(T2,A)') ''
      END IF

      CALL bs_env_release(qs_env%bs_env)

      CALL timestop(handle)

   END SUBROUTINE run_ri_rs_rpa

! **************************************************************************************************
!> \brief RPA with periodic RI-RS and full k-point sampling, the GW stages without Σᶜ and QP equation.
!>        E_c^RPA = 1/(2π N_k) Σ_q ∫_0^∞ dω Tr[ln ε(q,iω) - ε(q,iω) + I] is summed in compute_W0,
!>        E^EXX from compute_exchange_energy_ri_rs.
!>        No q = 0 head correction of ε is applied to E_c^RPA.
!> \param qs_env Quickstep environment containing the converged SCF state.
!> \param bs_env Band-structure settings, orbitals, and RI-RS work data.
!> \param e_exx exact-exchange energy per unit cell
!> \param e_c_rpa RPA correlation energy per unit cell
! **************************************************************************************************
   SUBROUTINE rpa_calc_ri_rs_small_cell_full_kp(qs_env, bs_env, e_exx, e_c_rpa)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(OUT)                         :: e_exx, e_c_rpa

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_calc_ri_rs_small_cell_full_kp'

      INTEGER                                            :: handle
      TYPE(gwkp_type), POINTER                           :: gwkp_env

      CALL timeset(routineN, handle)

      gwkp_env => bs_env%gwkp_env

      ! Z_lP(q) and ψ_nk(r_l)
      CALL prepare_ri_rs(bs_env, gwkp_env)

      ! χ_PQ(q,iω)
      CALL compute_chi(bs_env, gwkp_env)

      ! ε(q,iω), Wᶜ and E_c^RPA = 1/(2π N_k) Σ_q ∫_0^∞ dω Tr[ln ε(q,iω) - ε(q,iω) + I],
      ! stored in gwkp_env%e_c_rpa
      CALL compute_W0(qs_env, bs_env, gwkp_env)

      ! E^EXX = -s/2 Σ_σ Σ_R Σ_ll' V_ll'(R) D^σ_ll'(R)² + correction of the 1/q² singularity
      CALL compute_exchange_energy_ri_rs(bs_env, gwkp_env, e_exx)

      e_c_rpa = gwkp_env%e_c_rpa

      CALL release_gw_work(bs_env)
      NULLIFY (gwkp_env)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE rpa_calc_ri_rs_small_cell_full_kp

END MODULE rpa_ri_rs
