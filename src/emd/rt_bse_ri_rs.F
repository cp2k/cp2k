!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RT-BSE RI-RS kernels: SEX and Hartree evaluated by collocation on grid points r_l.
!>        Once-built grid objects:
!>          φ_µ(r_l)  [grid×AO],   Z_lP  [grid×RI],
!>          V^aux_PQ = [M^-1 V^tr M^-1]_PQ  [RI×RI],  W^0_ll' = sum_PQ Z_lP (V + W^c(ω=0))_PQ Z_l'Q.
!>        Per-call kernels, collocation X = φ_µ(r_l) (AO domain, written φ_lµ below):
!>          SEX:      ρ^grid_ll' = sum_µν φ_lµ Δρ_µν φ_l'ν ;
!>                    Σ_µν = pref * sum_ll' φ_lµ [ρ^grid ∘ W^0]_ll' φ_l'ν
!>          Hartree:  n_l = sum_µν φ_lµ Δρ_µν φ_lν ;
!>                    v_l = sum_PQl' Z_lP V^aux_PQ Z_l'Q n_l'   (applied factorized, stage 2) ;
!>                    V^H_µν = sum_l φ_lµ v_l φ_lν      (diagonal-only, no grid×grid)
!>        Independent of which GW variant produced bs_env%fm_W_MIC_freq_zero.
!> \author Maximilian Graml (05.26)
! **************************************************************************************************
MODULE rt_bse_ri_rs
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_distribution_type, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_release, dbcsr_scale, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_get_diag
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE gw_large_cell_Gamma_ri_rs,       ONLY: contract_A_B_A,&
                                              hadamard_product_inplace,&
                                              release_dbcsr_topology_and_matrices,&
                                              setup_square_topology
   USE gw_large_cell_gamma,             ONLY: multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_non_periodic_ri_rs,           ONLY: atomic_basis_at_grid_point,&
                                              compute_coeff_Z_lP,&
                                              precompute_ri_rs_radii,&
                                              ri_rs_grid_assembler
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_bse_ri_rs'

   PUBLIC :: rt_bse_ri_rs_ensure_grid, &
             rt_bse_ri_rs_ensure_V_grid, &
             rt_bse_ri_rs_ensure_W0_grid, &
             compute_sigma_ri_rs, &
             compute_sigma_ri_rs_complex, &
             compute_hartree_ri_rs, &
             compute_hartree_ri_rs_complex, &
             compute_hartree_ri_rs_from_diag, &
             hartree_potential_from_diag_ri_rs

CONTAINS

! **************************************************************************************************
!> \brief Make sure the AO collocation φ_µ(r_l) (mat_phi_mu_l) and the RI fit coefficients
!>        Z_lP (mat_Z_lP) are populated in memory.
!>        If GW was run with RI-RS the grid is already built; otherwise build it here so the
!>        AO-RI GW + RI-RS RT-BSE combination is possible.
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_grid'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      t1 = m_walltime()

      CALL ri_rs_grid_assembler(qs_env, bs_env, bs_env%ri_rs%grid_points)
      ! Per-atom AO/RI screening radii required by the screened grid-fill and d_lP routines; the
      ! GW RI-RS driver populates these, but the standalone RT-BSE grid build must do so itself.
      IF (.NOT. ALLOCATED(bs_env%ri_rs%radius_ao_per_atom)) THEN
         CALL precompute_ri_rs_radii(qs_env, bs_env)
      END IF
      CALL atomic_basis_at_grid_point(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                                      bs_env%ri_rs%mat_phi_mu_l)
      CALL compute_coeff_Z_lP(qs_env, bs_env, bs_env%ri_rs%grid_points, &
                              bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

      bs_env%ri_rs%grid_built = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Built RI-RS grid for RT-BSE (no GW_RI_RS used),', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_grid

! **************************************************************************************************
!> \brief Build V^aux_PQ = [M^-1 V^tr M^-1]_PQ (truncated Coulomb in the RI basis, M^-1-sandwiched
!>        to match the W^MIC convention). The grid Coulomb V_ll' = sum_PQ Z_lP V^aux_PQ Z_l'Q is
!>        never materialized -- the Hartree stage applies it factorized:
!>          v_l = sum_PQl' Z_lP V^aux_PQ Z_l'Q n_l'
!>        Needed by RT-BSE RI-RS Hartree.
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_V_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_V_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%V_grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL rt_bse_ri_rs_ensure_grid(bs_env, qs_env)
      t1 = m_walltime()

      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)
      ! Apply M^-1 sandwich to match the W^MIC convention; same scale as W^c when both used.
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      ! Store the M^-1-sandwiched RI-basis Coulomb; the grid kernel Z V Z^T is applied factorized.
      CALL dbcsr_create(bs_env%ri_rs%mat_V_aux_rtbse, "V_aux_rtbse", dist_aux_aux, &
                        dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), bs_env%ri_rs%mat_V_aux_rtbse, &
                            keep_sparsity=.FALSE.)

      bs_env%ri_rs%V_grid_built = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
            'Precomputed RT-BSE RI-RS V_aux kernel,', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_V_grid

! **************************************************************************************************
!> \brief Build W^0_ll' = sum_PQ Z_lP (V + W^c(ω=0))_PQ Z_l'Q (statically screened W on the grid).
!>        Needed by RT-BSE RI-RS SEX/COH. Reuses bs_env%fm_W_MIC_freq_zero which must already
!>        contain M^-1 W^c(ω=0) M^-1 (built by either GW path under the BSE rtp_method gate).
!>        W^0 enters only through the Hadamard ρ^grid ∘ W^0 -- its elements are needed, so it is
!>        the one persistent grid×grid object of the kernel layer.
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_W0_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_W0_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_V_aux, matrix_W_aux

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%W0_grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      ! W(w=0) is built by the GW step only under its RTBSE rtp_method gate; reaching this consumer
      ! without it means gate and consumer disagree. Abort rather than read a never-created cp_fm.
      IF (.NOT. ASSOCIATED(bs_env%fm_W_MIC_freq_zero%matrix_struct)) THEN
         CALL cp_abort(__LOCATION__, &
                       "RT-BSE RI-RS kernel needs the screened interaction W(w=0), which the GW "// &
                       "step did not build. Select the RT-BSE propagator with '&RTBSE' or "// &
                       "'&RTBSE RTBSE', not '&RTBSE TDDFT'.")
      END IF

      CALL rt_bse_ri_rs_ensure_grid(bs_env, qs_env)
      t1 = m_walltime()

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)
      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux_rtbse_W0", dist_aux_aux, dbcsr_type_no_symmetry, &
                        blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)

      CALL dbcsr_create(matrix_W_aux, "W_aux_rtbse", dist_aux_aux, dbcsr_type_no_symmetry, &
                        blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(bs_env%fm_W_MIC_freq_zero, matrix_W_aux, keep_sparsity=.FALSE.)
      CALL dbcsr_add(matrix_W_aux, matrix_V_aux, 1.0_dp, 1.0_dp)

      CALL dbcsr_create(bs_env%ri_rs%mat_W0_grid_rtbse, "W0_grid_rtbse", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_Z_lP, matrix_W_aux, &
                          bs_env%ri_rs%mat_W0_grid_rtbse, bs_env%eps_filter)

      bs_env%ri_rs%W0_grid_built = .TRUE.
      bs_env%ri_rs%rtbse_kernels_ready = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
            'Precomputed RT-BSE RI-RS W0_grid kernel,', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux, &
                                               m1=matrix_V_aux, m2=matrix_W_aux)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_W0_grid

! **************************************************************************************************
!> \brief AO-domain SEX:  Σ_µν = pref * sum_ll' φ_lµ [ρ^grid ∘ W^0]_ll' φ_l'ν,
!>        ρ^grid_ll' = sum_µν φ_lµ Δρ_µν φ_l'ν.
!>        The grid×grid ρ^grid is intrinsic to SEX -- the Hadamard needs W^0's elements, so no
!>        factorized application exists (unlike the Hartree V_ll'). Real input, real output;
!>        used for COH (input S^-1) and the init reference (ρ^0); dynamic Δρ goes through the
!>        complex variant. Mirrors the AO-RI get_sigma(rtbse_env, sigma_fm, prefactor, rho_fm) API.
!> \param bs_env ...
!> \param sigma_AO_fm result, AO x AO
!> \param prefactor scaling applied to the final result
!> \param rho_AO_fm input density-like matrix, AO x AO
!> \param grid_diag_accum ...
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs(bs_env, sigma_AO_fm, prefactor, rho_AO_fm, &
                                  grid_diag_accum)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: sigma_AO_fm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_AO_fm
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: grid_diag_accum(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs'

      INTEGER                                            :: handle, n_grid
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, blk_grid, dist_col_grid, &
                                                            dist_row_ao
      REAL(KIND=dp), ALLOCATABLE                         :: diag_local(:)
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_rho_AO, matrix_rho_grid, &
                                                            matrix_Sigma_AO

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%W0_grid_built)

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)

      CALL dbcsr_create(matrix_rho_AO, "rho_AO_ri_rs", dist_ao_ao, dbcsr_type_no_symmetry, &
                        blk_ao, blk_ao)
      CALL copy_fm_to_dbcsr(rho_AO_fm, matrix_rho_AO, keep_sparsity=.FALSE.)

      ! ρ^grid_ll' = sum_µν φ_lµ Δρ_µν φ_l'ν   (the SEX-intrinsic grid×grid transient)
      CALL dbcsr_create(matrix_rho_grid, "rho_grid_ri_rs", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      ! No eps_filter on the SEX projections: a Δρ-dependent block drop breaks kernel
      ! self-adjointness (L non-Hermitian). The W0 Hadamard below supplies the symmetric,
      ! Δρ-independent grid sparsity, so unfiltered projections stay sparse via W0.
      ! PERF(rirs-selfadjoint): a retain_sparsity-to-W0 build (dbcsr_multiply) gives identical
      ! numbers with W0-sparse memory; NOT YET APPLIED, needs a perf test.
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_phi_mu_l, matrix_rho_AO, &
                          matrix_rho_grid, 0.0_dp)

      ! Harvest n_l = diag(ρ^grid) BEFORE the in-place Hadamard destroys it; the Hartree reuses it
      ! (compute_hartree_ri_rs_from_diag) instead of rebuilding φρφ^T. Bare accumulate (caller
      ! pre-zeroes): the cross-spin Hartree density is the spin SUM of these, and spin_degeneracy is
      ! applied to the V_H OUTPUT (post-filter, bit-identical) -- never to the diagonal, which would
      ! shift the stage-3 eps_filter cut (coarse-filter sensitivity).
      IF (PRESENT(grid_diag_accum)) THEN
         n_grid = SIZE(grid_diag_accum)
         ALLOCATE (diag_local(n_grid))
         diag_local = 0.0_dp
         CALL dbcsr_get_diag(matrix_rho_grid, diag_local)
         CALL bs_env%para_env%sum(diag_local)
         grid_diag_accum(:) = grid_diag_accum(:) + diag_local(:)
         DEALLOCATE (diag_local)
      END IF

      ! ρ^grid_ll' <- ρ^grid_ll' * W^0_ll'   (in place; blocks without a W^0 partner zeroed)
      CALL hadamard_product_inplace(matrix_rho_grid, bs_env%ri_rs%mat_W0_grid_rtbse, &
                                    1.0_dp)

      ! Σ_µν = sum_ll' φ_lµ [ρ^grid ∘ W^0]_ll' φ_l'ν
      CALL dbcsr_create(matrix_Sigma_AO, template=matrix_rho_AO)
      ! Unfiltered for self-adjointness (see the forward projection above; PERF note there).
      CALL contract_A_B_A("T", "N", bs_env%ri_rs%mat_phi_mu_l, matrix_rho_grid, &
                          matrix_Sigma_AO, 0.0_dp)

      CALL dbcsr_scale(matrix_Sigma_AO, prefactor)
      CALL copy_dbcsr_to_fm(matrix_Sigma_AO, sigma_AO_fm)

      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, &
                                               m1=matrix_rho_AO, m2=matrix_Sigma_AO)
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=matrix_rho_grid)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs

! **************************************************************************************************
!> \brief Complex-input AO SEX via Re/Im split: the kernel is real, so complex linearity holds as
!>        Σ[Δρ] = Σ[Re Δρ] + i Σ[Im Δρ].  Required for non-Hermitian Δρ inputs
!>        (TDA OV-only / ABBA OV+VO).
!> \param bs_env ...
!> \param sigma_AO_cfm result, AO x AO (complex)
!> \param prefactor scaling applied to the final result
!> \param rho_AO_cfm input AO x AO complex matrix
!> \param grid_diag_re_accum optional: accumulate diag(φ.Re(ρ).φ^T) (bare; for the Hartree reuse)
!> \param grid_diag_im_accum optional: accumulate diag(φ.Im(ρ).φ^T) (bare; for the Hartree reuse)
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs_complex(bs_env, sigma_AO_cfm, prefactor, rho_AO_cfm, &
                                          grid_diag_re_accum, grid_diag_im_accum)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: sigma_AO_cfm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_AO_cfm
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: grid_diag_re_accum(:), &
                                                            grid_diag_im_accum(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs_complex'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_rho, fm_sigma

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_rho, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_sigma, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_real_part, bs_env%fm_s_Gamma%matrix_struct)

      ! Re/Im each harvest into their own accumulator; absent optionals propagate as absent.
      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargetr=fm_rho)
      CALL compute_sigma_ri_rs(bs_env, fm_sigma, prefactor, fm_rho, &
                               grid_diag_accum=grid_diag_re_accum)
      CALL cp_fm_to_cfm(msourcer=fm_sigma, mtarget=cfm_real_part)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargeti=fm_rho)
      CALL compute_sigma_ri_rs(bs_env, fm_sigma, prefactor, fm_rho, &
                               grid_diag_accum=grid_diag_im_accum)
      CALL cp_fm_to_cfm(msourcei=fm_sigma, mtarget=sigma_AO_cfm)

      CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), sigma_AO_cfm, &
                                CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)

      CALL cp_fm_release(fm_rho)
      CALL cp_fm_release(fm_sigma)
      CALL cp_cfm_release(cfm_real_part)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs_complex

! **************************************************************************************************
!> \brief AO-domain Hartree via RI-RS:
!>          n_l = sum_µν φ_lµ Δρ_µν φ_lν = (φ ρ φ^T)_ll  (diagonal of materialized grid×grid) ;
!>          v_l = sum_PQl' Z_lP V^aux_PQ Z_l'Q n_l'   (factorized, stage 2) ;
!>          V^H_µν = sum_l φ_lµ v_l φ_lν   (diagonal-only row-scale, stage 3)
!>        Real input, real output; complex inputs go through compute_hartree_ri_rs_complex.
!> \param bs_env ...
!> \param rho_AO_fm input AO x AO density matrix
!> \param V_H_AO_fm output AO x AO Hartree potential
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs(bs_env, rho_AO_fm, V_H_AO_fm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_AO_fm
      TYPE(cp_fm_type), INTENT(INOUT)                    :: V_H_AO_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs'

      INTEGER                                            :: handle, n_grid
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, blk_grid, dist_col_grid, &
                                                            dist_row_ao
      REAL(KIND=dp), ALLOCATABLE                         :: n_vec(:)
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_rho_AO, matrix_rho_grid

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%V_grid_built)

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)

      CALL dbcsr_create(matrix_rho_AO, "rho_AO_hartree", dist_ao_ao, dbcsr_type_no_symmetry, &
                        blk_ao, blk_ao)
      CALL copy_fm_to_dbcsr(rho_AO_fm, matrix_rho_AO, keep_sparsity=.FALSE.)

      n_grid = SUM(blk_grid)
      ALLOCATE (n_vec(n_grid))

      ! stage 1: n_l = (φ ρ φ^T)_ll. A grid×AO row-dot (iterate φρ, look up φ) silently drops
      ! off-rank pairs at ≥2 ranks — φ is sparse and the product is not co-located with φ — so
      ! materialize φρφ^T and take its diagonal, co-location-safe like the SEX grid kernel.
      CALL dbcsr_create(matrix_rho_grid, "rho_grid_hartree", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      ! No eps_filter: a Δρ-dependent diagonal-block drop zeroes n_l inconsistently between OV
      ! pairs and breaks Hartree self-adjointness (L non-Hermitian).
      ! PERF(rirs-selfadjoint): this full grid×grid is built only for its diagonal n_l; a
      ! diagonal-only build would avoid it. NOT YET APPLIED, needs a perf test.
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_phi_mu_l, matrix_rho_AO, &
                          matrix_rho_grid, 0.0_dp)
      n_vec = 0.0_dp
      CALL dbcsr_get_diag(matrix_rho_grid, n_vec)
      CALL bs_env%para_env%sum(n_vec)
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=matrix_rho_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, &
                                               m1=matrix_rho_AO)

      ! stages 2-3: factorized Coulomb v = Z V^aux Z^T n, then V^H = φ^T diag(v) φ.
      CALL hartree_potential_from_diag_ri_rs(bs_env, n_vec, V_H_AO_fm)

      DEALLOCATE (n_vec)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs

! **************************************************************************************************
!> \brief Hartree stages 2-3 from a precomputed grid density n_l (skips the stage-1 φρφ^T build):
!>          v_l = sum_PQl' Z_lP V^aux_PQ Z_l'Q n_l'   (factorized Coulomb) ;
!>          V^H_µν = sum_l φ_lµ v_l φ_lν   (Φ_lν = v_l φ_lν rowscale, then V^H = φ^T Φ).
!>        n_l is harvested as diag(φρφ^T) inside compute_sigma_ri_rs (the SEX grid kernel), so the
!>        Hartree never rebuilds the grid×grid product. Real in/out.
!> \param bs_env ...
!> \param n_vec grid density n_l (length n_grid, replicated)
!> \param V_H_AO_fm output AO x AO Hartree potential
! **************************************************************************************************
   SUBROUTINE hartree_potential_from_diag_ri_rs(bs_env, n_vec, V_H_AO_fm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: n_vec(:)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: V_H_AO_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hartree_potential_from_diag_ri_rs'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      REAL(KIND=dp), ALLOCATABLE                         :: u_ri(:), v_vec(:), w_ri(:)
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao
      TYPE(dbcsr_type)                                   :: matrix_Phi, matrix_V_H_AO

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%V_grid_built)

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      ALLOCATE (v_vec(SIZE(n_vec)))

      ! stage 2: v_l = sum_PQl' Z_lP V^aux_PQ Z_l'Q n_l'  (factorized; helper zeroes each output)
      ALLOCATE (w_ri(bs_env%n_RI), u_ri(bs_env%n_RI))
      ! w_Q = sum_l Z_lQ n_l
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_Z_lP, n_vec, w_ri, bs_env%para_env, &
                                   transposed=.TRUE.)
      ! u_P = sum_Q V^aux_PQ w_Q
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_V_aux_rtbse, w_ri, u_ri, bs_env%para_env)
      ! v_l = sum_P Z_lP u_P
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_Z_lP, u_ri, v_vec, bs_env%para_env)
      DEALLOCATE (w_ri, u_ri)

      ! stage 3: V^H_µν = sum_l φ_lµ v_l φ_lν  (Φ_lν = v_l φ_lν rowscale, then V^H = φ^T Φ)
      CALL dbcsr_create(matrix_Phi, template=bs_env%ri_rs%mat_phi_mu_l)
      CALL dbcsr_copy(matrix_Phi, bs_env%ri_rs%mat_phi_mu_l)
      CALL dbcsr_scale_rows_replicated(matrix_Phi, v_vec)

      CALL dbcsr_create(matrix_V_H_AO, "V_H_AO_hartree", dist_ao_ao, dbcsr_type_no_symmetry, &
                        blk_ao, blk_ao)
      ! No eps_filter: a Δρ-dependent block drop breaks Hartree kernel self-adjointness.
      CALL dbcsr_multiply("T", "N", 1.0_dp, bs_env%ri_rs%mat_phi_mu_l, matrix_Phi, &
                          0.0_dp, matrix_V_H_AO, filter_eps=0.0_dp)
      CALL dbcsr_release(matrix_Phi)

      CALL copy_dbcsr_to_fm(matrix_V_H_AO, V_H_AO_fm)

      DEALLOCATE (v_vec)

      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, &
                                               m1=matrix_V_H_AO)

      CALL timestop(handle)

   END SUBROUTINE hartree_potential_from_diag_ri_rs

! **************************************************************************************************
!> \brief Complex Hartree from precomputed grid diagonals: V^H = V^H[n_re] + i V^H[n_im], each via
!>        hartree_potential_from_diag_ri_rs (stages 2-3 only). n_re/n_im are the spin-summed grid
!>        densities harvested in the SEX kernel; this is the cross-spin / TDA complex consumer that
!>        replaces compute_hartree_ri_rs_complex when SEX already built the grid. n_im optional: when
!>        absent the result is purely real (matches the Re-only real-input Hartree).
!> \param bs_env ...
!> \param n_re grid density Re part (length n_grid, replicated)
!> \param V_H_AO_cfm output AO x AO complex Hartree potential
!> \param n_im optional grid density Im part (length n_grid, replicated)
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs_from_diag(bs_env, n_re, V_H_AO_cfm, n_im)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: n_re(:)
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: V_H_AO_cfm
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: n_im(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs_from_diag'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_v

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_v, bs_env%fm_s_Gamma%matrix_struct)

      CALL hartree_potential_from_diag_ri_rs(bs_env, n_re, fm_v)
      IF (PRESENT(n_im)) THEN
         CALL cp_cfm_create(cfm_real_part, bs_env%fm_s_Gamma%matrix_struct)
         CALL cp_fm_to_cfm(msourcer=fm_v, mtarget=cfm_real_part)
         CALL hartree_potential_from_diag_ri_rs(bs_env, n_im, fm_v)
         CALL cp_fm_to_cfm(msourcei=fm_v, mtarget=V_H_AO_cfm)
         CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), V_H_AO_cfm, &
                                   CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)
         CALL cp_cfm_release(cfm_real_part)
      ELSE
         CALL cp_fm_to_cfm(msourcer=fm_v, mtarget=V_H_AO_cfm)
      END IF

      CALL cp_fm_release(fm_v)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs_from_diag

! **************************************************************************************************
!> \brief Complex-input Hartree potential via RI-RS. Re/Im split: feed each part to the real
!>        compute_hartree_ri_rs and reassemble. Real-input Hartree on a non-Hermitian input
!>        would silently drop Im and break Hermitian conjugacy of OV+VO contributions in TDA.
!> \param bs_env ...
!> \param rho_AO_cfm input AO x AO complex density-like matrix
!> \param V_H_AO_cfm output AO x AO complex Hartree potential
! **************************************************************************************************
   SUBROUTINE compute_hartree_ri_rs_complex(bs_env, rho_AO_cfm, V_H_AO_cfm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_AO_cfm
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: V_H_AO_cfm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_hartree_ri_rs_complex'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_rho, fm_v

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_rho, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_v, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_real_part, bs_env%fm_s_Gamma%matrix_struct)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargetr=fm_rho)
      CALL compute_hartree_ri_rs(bs_env, fm_rho, fm_v)
      CALL cp_fm_to_cfm(msourcer=fm_v, mtarget=cfm_real_part)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargeti=fm_rho)
      CALL compute_hartree_ri_rs(bs_env, fm_rho, fm_v)
      CALL cp_fm_to_cfm(msourcei=fm_v, mtarget=V_H_AO_cfm)

      CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), V_H_AO_cfm, &
                                CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)

      CALL cp_fm_release(fm_rho)
      CALL cp_fm_release(fm_v)
      CALL cp_cfm_release(cfm_real_part)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs_complex

! **************************************************************************************************
!> \brief Scale each row of a dbcsr matrix by a replicated full-length vector:
!>        block(ir,ic) <- vec(global_row(ir)) * block(ir,ic). Value mutation only.
!> \param matrix ...
!> \param vec full-length replicated row-scaling vector
! **************************************************************************************************
   SUBROUTINE dbcsr_scale_rows_replicated(matrix, vec)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: vec(:)

      INTEGER                                            :: col_blk, ic, ir, nblkrows_total, &
                                                            row_blk, row_off
      INTEGER, ALLOCATABLE                               :: row_offset(:)
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows_total, row_blk_size=row_blk_size)
      ALLOCATE (row_offset(nblkrows_total + 1))
      row_offset(1) = 0
      DO ir = 1, nblkrows_total
         row_offset(ir + 1) = row_offset(ir) + row_blk_size(ir)
      END DO

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row_blk, col_blk, blk)
         row_off = row_offset(row_blk)
         DO ic = 1, SIZE(blk, 2)
            DO ir = 1, SIZE(blk, 1)
               blk(ir, ic) = vec(row_off + ir)*blk(ir, ic)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (row_offset)

   END SUBROUTINE dbcsr_scale_rows_replicated

! **************************************************************************************************
!> \brief Replicated matvec: vec_out = matrix * vec_in (or matrix^T * vec_in if transposed) for a
!>        possibly rectangular distributed dbcsr matrix and replicated full-length vectors.
!>        Iterates over local blocks and reduces.
!> \param matrix distributed dbcsr matrix (may be rectangular)
!> \param vec_in full input vector, replicated on all ranks (column length, or row length if transposed)
!> \param vec_out full output vector, replicated on all ranks (zeroed on entry; sum-reduced on exit)
!> \param para_env ...
!> \param transposed if .TRUE. compute vec_out = matrix^T * vec_in
! **************************************************************************************************
   SUBROUTINE dbcsr_matvec_replicated(matrix, vec_in, vec_out, para_env, transposed)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: vec_in(:)
      REAL(KIND=dp), INTENT(INOUT)                       :: vec_out(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      INTEGER                                            :: col_blk, col_off, ic, ir, &
                                                            nblkcols_total, nblkrows_total, &
                                                            row_blk, row_off
      INTEGER, ALLOCATABLE                               :: col_offset(:), row_offset(:)
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      LOGICAL                                            :: my_trans
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

      my_trans = .FALSE.
      IF (PRESENT(transposed)) my_trans = transposed

      CALL dbcsr_get_info(matrix, nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size)

      ALLOCATE (row_offset(nblkrows_total + 1), col_offset(nblkcols_total + 1))
      row_offset(1) = 0
      DO ir = 1, nblkrows_total
         row_offset(ir + 1) = row_offset(ir) + row_blk_size(ir)
      END DO
      col_offset(1) = 0
      DO ic = 1, nblkcols_total
         col_offset(ic + 1) = col_offset(ic) + col_blk_size(ic)
      END DO

      vec_out(:) = 0.0_dp

      CALL dbcsr_iterator_start(iter, matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row_blk, col_blk, block)
         row_off = row_offset(row_blk)
         col_off = col_offset(col_blk)
         IF (my_trans) THEN
            DO ic = 1, SIZE(block, 2)
               DO ir = 1, SIZE(block, 1)
                  vec_out(col_off + ic) = vec_out(col_off + ic) + &
                                          block(ir, ic)*vec_in(row_off + ir)
               END DO
            END DO
         ELSE
            DO ic = 1, SIZE(block, 2)
               DO ir = 1, SIZE(block, 1)
                  vec_out(row_off + ir) = vec_out(row_off + ir) + &
                                          block(ir, ic)*vec_in(col_off + ic)
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (row_offset, col_offset)

      CALL para_env%sum(vec_out)

   END SUBROUTINE dbcsr_matvec_replicated

END MODULE rt_bse_ri_rs
