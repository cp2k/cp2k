!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Work data of periodic RI-RS GW with full k-point sampling.
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE gw_ri_rs_small_cell_full_kp_types
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE kpoint_transform,                ONLY: kp_transform_plan_type
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_small_cell_full_kp_types'

   PUBLIC :: gamma_atom_projection_type, gamma_correction_type, gwkp_type, &
             ri_q_matrix_block_type

   ! Gamma correction of W per real-space cell R, kept factored from the grid vectors:
   ! head_cell_tau(R,τ) multiplies Z_l0 Z_l'0*, and the wing of R at τ is
   ! Σ_ω cos(τ,ω) wing_cartesian(:,:,ω) wing_cell_frequency(:,ω,R).
   TYPE :: gamma_correction_type
      COMPLEX(KIND=dp), ALLOCATABLE :: head_cell_tau(:, :), wing_cartesian(:, :, :), &
                                       wing_cell_frequency(:, :, :)
   END TYPE gamma_correction_type

   TYPE gamma_atom_projection_type
      INTEGER :: first = 0
      INTEGER, ALLOCATABLE :: positions(:)
      REAL(KIND=dp), ALLOCATABLE :: values(:, :)
   END TYPE gamma_atom_projection_type

   !> values(P,Q,j) stores one q point, with j an imaginary frequency or time.
   !> Each q allocation can be released after its last use independently.
   TYPE ri_q_matrix_block_type
      COMPLEX(KIND=dp), ALLOCATABLE :: values(:, :, :)
   END TYPE ri_q_matrix_block_type

   ! Work data shared by the periodic GW stages; input settings remain in bs_env.
   TYPE :: gwkp_type
      REAL(KIND=dp), ALLOCATABLE                 :: eigenvalues_G(:, :, :)
      COMPLEX(dp), ALLOCATABLE                  :: metric_coulomb_factor(:, :, :), psi_lm_k(:, :, :, :), &
                                                   Z_lP_q(:, :, :), Z_lP_zero(:, :, :)
      REAL(dp), ALLOCATABLE                     :: sigma_c_bands(:, :, :, :, :), sigma_x_bands(:, :, :)
      ! Band Bloch orbitals on the RI-RS grid, only of the band k points owned in the Σ ranks, MOD(k-1,P).
      COMPLEX(dp), ALLOCATABLE                  :: band_orbitals(:, :, :, :)
      ! Charge neutralisation of the fitted pair densities in χ(q), see compute_W0: e(q) = M(q) c/|c| and
      ! X(q,iω) = M(q) (s/N_k) Σ_k,ia d_ia m_ia* w_ia(iω) per owned q, L0(iω) = (s/N_k) Σ_k,ia |m_ia|² w_ia(iω).
      COMPLEX(dp), ALLOCATABLE                  :: leak_e(:, :), leak_x(:, :, :)
      REAL(dp), ALLOCATABLE                     :: leak_l0(:)
      TYPE(mp_para_env_type), POINTER           :: group_env => NULL()
      INTEGER                                   :: ngroups = 1, my_group = 0
      TYPE(gamma_correction_type)               :: gamma_correction
      TYPE(ri_q_matrix_block_type), ALLOCATABLE :: chi_q_frequency(:), W_PQ_tau(:)
      ! AO image operators retained only until the Gamma-correction contraction.
      TYPE(dbcsr_p_type), POINTER               :: gamma_corr_d(:, :) => NULL(), &
                                                   gamma_corr_h(:, :) => NULL(), gamma_corr_s(:, :) => NULL()
      TYPE(kp_transform_plan_type)              :: gamma_corr_plan, gamma_corr_scf_plan
      INTEGER, POINTER                          :: gamma_corr_index_to_cell(:, :) => NULL(), &
                                                   gamma_corr_cell_to_index(:, :, :) => NULL()
   END TYPE gwkp_type

END MODULE gw_ri_rs_small_cell_full_kp_types
