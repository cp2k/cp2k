!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using the RI-RS approximation for large periodic cells at Gamma
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_ri_rs_large_cell_Gamma
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE gw_compute_V_and_W,              ONLY: get_W_MIC
   USE gw_ri_rs_chi_sigma_kernels,      ONLY: compute_Sigma_c,&
                                              compute_Sigma_x,&
                                              get_mat_chi_Gamma_tau
   USE gw_tensor_large_cell_Gamma,      ONLY: compute_QP_energies
   USE gw_utils,                        ONLY: de_init_bs_env
   USE kinds,                           ONLY: dp
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_setup,                     ONLY: ensure_ri_rs_representation
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_large_cell_Gamma'

   PUBLIC :: gw_calc_ri_rs_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for large periodic cells at Gamma
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************

   SUBROUTINE gw_calc_ri_rs_large_cell_Gamma(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_large_cell_Gamma'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      ! Build the common grid, phi and Z representation with the system's boundary conditions.
      CALL ensure_ri_rs_representation(qs_env, bs_env)

       !!========================================================================
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|)  = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^vir_µλ(i|τ|)  = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^occ_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')
       !!    G^vir_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^vir_µν ϕ_ν(r_l')
       !!    χ_ll'(iτ)       = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
       !!    χ_PQ(iτ)        = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, cutoff=0.0_dp)

       !!========================================================================
       !! 5. Compute Screened Interaction (W^MIC)
       !!    χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W^MIC_PQ(iτ)
       !!========================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn C_νn
       !!    D_ll'         = sum_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma, cutoff=0.0_dp)

       !!========================================================================
       !! 7. Compute Correlation Self-Energy (Σ^c)
       !!    W^MIC_ll'(iτ) =  sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!    Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W^MIC_ll'(iτ), for τ < 0
       !!    Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W^MIC_ll'(iτ), for τ > 0
       !!    Σ^c_λσ(iτ)    =  sum_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_c(bs_env, fm_W_time, fm_Sigma_c_Gamma_time, cutoff=0.0_dp)

       !!========================================================================
       !! 8. Compute Quasiparticle Energies
       !!    Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
       !!    ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ) + Σ^x_nn - v^xc_nn
       !!========================================================================
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_large_cell_Gamma

END MODULE gw_ri_rs_large_cell_Gamma
