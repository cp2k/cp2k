!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RT-BSE specific RI-RS kernels (V_grid, W0_grid) and contraction routines.
!>        Independent of which GW variant produced bs_env%fm_W_MIC_freq_zero: this module
!>        builds the RI-RS grid on demand (mat_phi_mu_l, mat_Z_lP) if missing, and assembles
!>        the V_grid and W0_grid kernels that the linearized RT-BSE propagation contracts
!>        against rho/S^-1 to obtain Sigma and the Hartree potential.
!> \author Maximilian Graml
!> \date 05.2026
! **************************************************************************************************
MODULE rt_bse_ri_rs
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create, &
                                              cp_cfm_release, &
                                              cp_cfm_to_fm, &
                                              cp_cfm_type, &
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: dbcsr_add, &
                                              dbcsr_create, &
                                              dbcsr_distribution_type, &
                                              dbcsr_get_info, &
                                              dbcsr_iterator_blocks_left, &
                                              dbcsr_iterator_next_block, &
                                              dbcsr_iterator_start, &
                                              dbcsr_iterator_stop, &
                                              dbcsr_iterator_type, &
                                              dbcsr_release, &
                                              dbcsr_scale, &
                                              dbcsr_type, &
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_add_on_diag, &
                                              dbcsr_get_diag, &
                                              dbcsr_set_diag
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, &
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_create, &
                                              cp_fm_release, &
                                              cp_fm_type
   USE gw_large_cell_gamma,             ONLY: multiply_fm_W_MIC_time_with_Minv_Gamma
   USE gw_large_cell_Gamma_ri_rs,       ONLY: atomic_basis_at_grid_point, &
                                              compute_coeff_Z_lP, &
                                              contract_A_B_A, &
                                              hadamard_product, &
                                              release_dbcsr_topology_and_matrices, &
                                              ri_rs_grid_assembler, &
                                              setup_square_topology
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
             compute_hartree_ri_rs_complex

CONTAINS

! **************************************************************************************************
!> \brief Make sure mat_phi_mu_l and mat_Z_lP are populated in memory.
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
!> \brief Build V_grid(l,l') = sum_PQ Z_lP V_PQ Z_l'Q (RI Coulomb projected on the grid).
!>        Needed by RT-BSE RI-RS Hartree.
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rt_bse_ri_rs_ensure_V_grid(bs_env, qs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_bse_ri_rs_ensure_V_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_V_aux

      CALL timeset(routineN, handle)

      IF (bs_env%ri_rs%V_grid_built) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL rt_bse_ri_rs_ensure_grid(bs_env, qs_env)
      t1 = m_walltime()

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)
      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)
      ! Apply M^-1 sandwich to match the W^MIC convention; same scale as W^c when both used.
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux_rtbse", dist_aux_aux, dbcsr_type_no_symmetry, &
                        blk_aux, blk_aux)
      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)

      CALL dbcsr_create(bs_env%ri_rs%mat_V_grid_rtbse, "V_grid_rtbse", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_Z_lP, matrix_V_aux, &
                          bs_env%ri_rs%mat_V_grid_rtbse, bs_env%eps_filter)

      bs_env%ri_rs%V_grid_built = .TRUE.

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Precomputed RT-BSE RI-RS V_grid kernel,', ' Execution time', &
            m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux, &
                                               m1=matrix_V_aux)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE rt_bse_ri_rs_ensure_V_grid

! **************************************************************************************************
!> \brief Build W0_grid(l,l') = sum_PQ Z_lP (V+W^c)_PQ Z_l'Q (full W at omega=0 on the grid).
!>        Needed by RT-BSE RI-RS SEX/COH. Reuses bs_env%fm_W_MIC_freq_zero which must already
!>        contain M^-1 W^c(w=0) M^-1 (built by either GW path under the BSE rtp_method gate).
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
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
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
!> \brief Real-input self-energy via RI-RS grid kernel:
!>        Sigma_AO(mu,nu) = prefactor * sum_l phi(l,mu) [ rho_grid(l,l') . W0_grid(l,l') ] phi(l',nu)
!>        where rho_grid = phi * rho_AO * phi^T (l x l).
!>        Mirrors the AO-RI get_sigma(rtbse_env, sigma_fm, prefactor, rho_fm) API.
!> \param bs_env ...
!> \param sigma_AO_fm result, AO x AO
!> \param prefactor scaling applied to the final result
!> \param rho_AO_fm input density-like matrix, AO x AO
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs(bs_env, sigma_AO_fm, prefactor, rho_AO_fm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: sigma_AO_fm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_fm_type), INTENT(IN)                       :: rho_AO_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, blk_grid, dist_col_grid, &
                                                            dist_row_ao
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_rho_AO, matrix_rho_grid, &
                                                            matrix_Sigma_AO, matrix_Sigma_grid

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%W0_grid_built)

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)

      CALL dbcsr_create(matrix_rho_AO, "rho_AO_ri_rs", dist_ao_ao, dbcsr_type_no_symmetry, &
                        blk_ao, blk_ao)
      CALL copy_fm_to_dbcsr(rho_AO_fm, matrix_rho_AO, keep_sparsity=.FALSE.)

      CALL dbcsr_create(matrix_rho_grid, "rho_grid_ri_rs", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_phi_mu_l, matrix_rho_AO, &
                          matrix_rho_grid, bs_env%eps_filter)

      CALL dbcsr_create(matrix_Sigma_grid, template=matrix_rho_grid)
      CALL hadamard_product(matrix_rho_grid, bs_env%ri_rs%mat_W0_grid_rtbse, &
                            matrix_Sigma_grid, 1.0_dp)

      CALL dbcsr_create(matrix_Sigma_AO, template=matrix_rho_AO)
      CALL contract_A_B_A("T", "N", bs_env%ri_rs%mat_phi_mu_l, matrix_Sigma_grid, &
                          matrix_Sigma_AO, bs_env%eps_filter)

      CALL dbcsr_scale(matrix_Sigma_AO, prefactor)
      CALL copy_dbcsr_to_fm(matrix_Sigma_AO, sigma_AO_fm)

      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, &
                                               m1=matrix_rho_AO, m2=matrix_Sigma_AO)
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=matrix_rho_grid, m2=matrix_Sigma_grid)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs

! **************************************************************************************************
!> \brief Complex-input self-energy via RI-RS: applies compute_sigma_ri_rs to Re and Im parts
!>        independently and reassembles. Required for non-Hermitian Delta-rho inputs in TDA
!>        (kernel must be complex-linear, see rt_bse_linearized notes).
!> \param bs_env ...
!> \param sigma_AO_cfm result, AO x AO (complex)
!> \param prefactor scaling applied to the final result
!> \param rho_AO_cfm input AO x AO complex matrix
! **************************************************************************************************
   SUBROUTINE compute_sigma_ri_rs_complex(bs_env, sigma_AO_cfm, prefactor, rho_AO_cfm)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: sigma_AO_cfm
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_AO_cfm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_sigma_ri_rs_complex'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_real_part
      TYPE(cp_fm_type)                                   :: fm_rho, fm_sigma

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_rho, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_sigma, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_real_part, bs_env%fm_s_Gamma%matrix_struct)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargetr=fm_rho)
      CALL compute_sigma_ri_rs(bs_env, fm_sigma, prefactor, fm_rho)
      CALL cp_fm_to_cfm(msourcer=fm_sigma, mtarget=cfm_real_part)

      CALL cp_cfm_to_fm(msource=rho_AO_cfm, mtargeti=fm_rho)
      CALL compute_sigma_ri_rs(bs_env, fm_sigma, prefactor, fm_rho)
      CALL cp_fm_to_cfm(msourcei=fm_sigma, mtarget=sigma_AO_cfm)

      CALL cp_cfm_scale_and_add(CMPLX(1.0_dp, 0.0_dp, kind=dp), sigma_AO_cfm, &
                                CMPLX(1.0_dp, 0.0_dp, kind=dp), cfm_real_part)

      CALL cp_fm_release(fm_rho)
      CALL cp_fm_release(fm_sigma)
      CALL cp_cfm_release(cfm_real_part)

      CALL timestop(handle)

   END SUBROUTINE compute_sigma_ri_rs_complex

! **************************************************************************************************
!> \brief Hartree potential via RI-RS grid kernel:
!>        V_H_AO(la,si) = sum_l phi(l,la) v_l phi(l,si),  v_l = sum_l' V_grid(l,l') n_l',
!>        n_l = (phi * rho_AO * phi^T)_{ll} = density on grid point l.
!>        Real input, real output. Used as building block for the complex variant.
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
      REAL(KIND=dp), ALLOCATABLE                         :: n_vec(:), v_vec(:)
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao, dist_grid_grid
      TYPE(dbcsr_type)                                   :: matrix_D_grid, matrix_rho_AO, &
                                                            matrix_rho_grid, matrix_V_H_AO

      CALL timeset(routineN, handle)

      CPASSERT(bs_env%ri_rs%V_grid_built)

      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL setup_square_topology(bs_env%ri_rs%mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, &
                                 dist_col_grid)

      CALL dbcsr_create(matrix_rho_AO, "rho_AO_hartree", dist_ao_ao, dbcsr_type_no_symmetry, &
                        blk_ao, blk_ao)
      CALL copy_fm_to_dbcsr(rho_AO_fm, matrix_rho_AO, keep_sparsity=.FALSE.)

      CALL dbcsr_create(matrix_rho_grid, "rho_grid_hartree", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL contract_A_B_A("N", "T", bs_env%ri_rs%mat_phi_mu_l, matrix_rho_AO, &
                          matrix_rho_grid, bs_env%eps_filter)

      n_grid = SUM(blk_grid)
      ALLOCATE (n_vec(n_grid), v_vec(n_grid))
      n_vec = 0.0_dp
      CALL dbcsr_get_diag(matrix_rho_grid, n_vec)
      CALL bs_env%para_env%sum(n_vec)

      v_vec = 0.0_dp
      CALL dbcsr_matvec_replicated(bs_env%ri_rs%mat_V_grid_rtbse, n_vec, v_vec, bs_env%para_env)

      CALL dbcsr_create(matrix_D_grid, "D_grid_hartree", dist_grid_grid, &
                        dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_add_on_diag(matrix_D_grid, 1.0_dp)
      CALL dbcsr_set_diag(matrix_D_grid, v_vec)

      CALL dbcsr_create(matrix_V_H_AO, template=matrix_rho_AO)
      CALL contract_A_B_A("T", "N", bs_env%ri_rs%mat_phi_mu_l, matrix_D_grid, &
                          matrix_V_H_AO, bs_env%eps_filter)

      CALL copy_dbcsr_to_fm(matrix_V_H_AO, V_H_AO_fm)

      DEALLOCATE (n_vec, v_vec)

      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, &
                                               m1=matrix_rho_AO, m2=matrix_V_H_AO)
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=matrix_rho_grid, m2=matrix_D_grid)

      CALL timestop(handle)

   END SUBROUTINE compute_hartree_ri_rs

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
!> \brief Replicated matvec: vec_out(i) = sum_j matrix(i,j) vec_in(j) for a distributed dbcsr
!>        matrix and replicated full-length vectors. Iterates over local blocks and reduces.
!> \param matrix square distributed dbcsr matrix
!> \param vec_in full input vector, replicated on all ranks
!> \param vec_out full output vector, replicated on all ranks (zeroed on entry; sum-reduced on exit)
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE dbcsr_matvec_replicated(matrix, vec_in, vec_out, para_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      REAL(KIND=dp), INTENT(IN)                          :: vec_in(:)
      REAL(KIND=dp), INTENT(INOUT)                       :: vec_out(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: col_blk, col_off, ic, ir, &
                                                            nblkcols_total, nblkrows_total, &
                                                            row_blk, row_off
      INTEGER, ALLOCATABLE                               :: col_offset(:), row_offset(:)
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbcsr_iterator_type)                          :: iter

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
         DO ic = 1, SIZE(block, 2)
            DO ir = 1, SIZE(block, 1)
               vec_out(row_off + ir) = vec_out(row_off + ir) + &
                                       block(ir, ic)*vec_in(col_off + ic)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (row_offset, col_offset)

      CALL para_env%sum(vec_out)

   END SUBROUTINE dbcsr_matvec_replicated

END MODULE rt_bse_ri_rs
