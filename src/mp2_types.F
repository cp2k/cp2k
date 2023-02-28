!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! **************************************************************************************************
MODULE mp2_types
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_finalize,&
                                              cp_eri_mme_param
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE hfx_types,                       ONLY: hfx_release,&
                                              hfx_type,&
                                              pair_list_element_type
   USE input_constants,                 ONLY: do_eri_mme,&
                                              mp2_method_direct,&
                                              mp2_method_gpw,&
                                              mp2_method_none,&
                                              mp2_ri_optimize_basis,&
                                              ri_mp2_laplace,&
                                              ri_mp2_method_gpw,&
                                              ri_rpa_method_gpw
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE iso_c_binding,                   ONLY: c_ptr
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU,&
                                              local_gemm_create,&
                                              local_gemm_destroy,&
                                              local_gemm_set_op_threshold_gpu
   USE message_passing,                 ONLY: mp_request_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_types'

   PUBLIC :: mp2_type, &
             integ_mat_buffer_type, &
             integ_mat_buffer_type_2D, &
             mp2_method_none, &
             mp2_method_direct, &
             mp2_method_gpw, &
             mp2_ri_optimize_basis, &
             ri_mp2_method_gpw, &
             ri_rpa_method_gpw, &
             ri_mp2_laplace, &
             init_TShPSC_lmax

   PUBLIC :: mp2_env_create, &
             mp2_env_release, &
             mp2_biel_type, &
             pair_list_type_mp2, &
             one_dim_int_array, &
             two_dim_int_array, &
             one_dim_real_array, &
             two_dim_real_array, &
             three_dim_real_array

   INTEGER, SAVE                                         :: init_TShPSC_lmax = -1

! TYPE definitions

   TYPE one_dim_int_array
      INTEGER, DIMENSION(:), ALLOCATABLE    :: array
   END TYPE

   TYPE two_dim_int_array
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: array
   END TYPE

   TYPE one_dim_real_array
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: array
   END TYPE

   TYPE two_dim_real_array
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: array
   END TYPE

   TYPE three_dim_real_array
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE :: array
   END TYPE

   TYPE mp2_biel_type
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: index_table
   END TYPE mp2_biel_type

   TYPE mp2_laplace_type
      INTEGER       :: n_quadrature
      INTEGER       :: num_integ_groups
   END TYPE

   TYPE mp2_direct_type
      LOGICAL  :: big_send
   END TYPE

   TYPE mp2_gpw_type
      REAL(KIND=dp)            :: eps_grid, eps_filter, eps_pgf_orb_S
      INTEGER                  :: print_level
      REAL(KIND=dp)            :: cutoff
      REAL(KIND=dp)            :: relative_cutoff
      INTEGER                  :: size_lattice_sum
   END TYPE mp2_gpw_type

   TYPE ri_mp2_type
      INTEGER                  :: block_size, number_integration_groups
      LOGICAL                  :: print_dgemm_info
   END TYPE

   TYPE ri_rpa_type
      INTEGER                  :: rpa_num_quad_points
      INTEGER                  :: rpa_num_integ_groups
      INTEGER                  :: mm_style
      TYPE(hfx_type), DIMENSION(:, :), POINTER &
         :: x_data => NULL()
      TYPE(section_vals_type), POINTER         :: xc_section_primary => Null(), &
                                                  xc_section_aux => Null()
      LOGICAL                  :: reuse_hfx
      LOGICAL                  :: minimax_quad
      LOGICAL                  :: do_ri_g0w0
      LOGICAL                  :: do_admm
      LOGICAL                  :: do_ri_axk
      LOGICAL                  :: do_rse
      LOGICAL                  :: print_dgemm_info
      TYPE(dbcsr_type), POINTER             :: mo_coeff_o, &
                                               mo_coeff_v
      REAL(KIND=dp)            :: ener_axk
      REAL(KIND=dp)            :: rse_corr_diag
      REAL(KIND=dp)            :: rse_corr
      REAL(KIND=dp)            :: scale_rpa
   END TYPE

   TYPE ri_rpa_im_time_type
      INTEGER                  :: cut_memory
      LOGICAL                  :: memory_info, make_chi_pos_definite, make_overlap_mat_ao_pos_definite, &
                                  trunc_coulomb_ri_x, keep_quad, &
                                  do_kpoints_from_Gamma, do_extrapolate_kpoints
      REAL(KIND=dp)            :: eps_filter, &
                                  eps_filter_factor, eps_compress, exp_tailored_weights, regularization_RI, &
                                  eps_eigval_S, eps_eigval_S_Gamma, rel_cutoff_trunc_coulomb_ri_x
      REAL(KIND=dp), DIMENSION(:), POINTER :: tau_tj, tau_wj, tj, wj
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: weights_cos_tf_t_to_w, weights_cos_tf_w_to_t
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: Eigenval_Gamma, wkp_V
      INTEGER                  :: group_size_P, group_size_3c, kpoint_weights_W_method, k_mesh_g_factor
      INTEGER, DIMENSION(:), POINTER     :: kp_grid
      INTEGER, DIMENSION(3) ::  kp_grid_extra
      LOGICAL                  :: do_im_time_kpoints
      INTEGER                  :: min_bsize, min_bsize_mo, nkp_orig, nkp_extra
      TYPE(kpoint_type), POINTER :: kpoints_G, kpoints_Sigma, kpoints_Sigma_no_xc
      INTEGER, ALLOCATABLE, DIMENSION(:)      :: starts_array_mc_RI, ends_array_mc_RI, &
                                                 starts_array_mc_block_RI, &
                                                 ends_array_mc_block_RI, &
                                                 starts_array_mc, ends_array_mc, &
                                                 starts_array_mc_block, &
                                                 ends_array_mc_block

   END TYPE

   TYPE ri_g0w0_type
      INTEGER                  :: corr_mos_occ
      INTEGER                  :: corr_mos_virt
      INTEGER                  :: corr_mos_occ_beta
      INTEGER                  :: corr_mos_virt_beta
      INTEGER                  :: num_poles
      INTEGER                  :: nparam_pade
      INTEGER                  :: analytic_continuation
      REAL(KIND=dp)            :: omega_max_fit
      INTEGER                  :: crossing_search
      REAL(KIND=dp)            :: fermi_level_offset
      INTEGER                  :: iter_evGW, iter_sc_GW0
      REAL(KIND=dp)            :: eps_iter
      LOGICAL                  :: do_ri_Sigma_x, &
                                  do_periodic, &
                                  print_self_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: vec_Sigma_x_minus_vxc_gw
      INTEGER, DIMENSION(:), POINTER    :: kp_grid, kp_grid_Sigma
      INTEGER                  :: num_kp_grids
      REAL(KIND=dp)            :: eps_kpoint
      LOGICAL                  :: do_mo_coeff_gamma, do_average_deg_levels
      REAL(KIND=dp)            :: eps_eigenval
      LOGICAL                  :: do_extra_kpoints, do_aux_bas_gw
      REAL(KIND=dp)            :: frac_aux_mos
      INTEGER                  :: num_omega_points
      LOGICAL                  :: do_bse
      INTEGER                  :: num_z_vectors, max_iter_bse
      REAL(KIND=dp)            :: eps_min_trans
      LOGICAL                  :: do_ic_model, print_ic_values
      REAL(KIND=dp)            :: eps_dist
      TYPE(one_dim_real_array), DIMENSION(2) :: ic_corr_list
      INTEGER :: print_exx
      LOGICAL :: do_gamma_only_sigma
      LOGICAL :: update_xc_energy, do_kpoints_Sigma, print_local_bandgap
      INTEGER :: n_kp_in_kp_line, n_special_kp, nkp_self_energy, &
                 nkp_self_energy_special_kp, nkp_self_energy_monkh_pack
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: xkp_special_kp
      TYPE(dbcsr_p_type), DIMENSION(:), ALLOCATABLE :: &
         matrix_sigma_x_minus_vxc, matrix_ks
      REAL(KIND=dp)            :: broadening_print_loc_bandgap, energy_window_print_loc_bandgap, &
                                  ldos_thresh_print_loc_bandgap, energy_spacing_print_loc_bandgap
      INTEGER, DIMENSION(:), POINTER :: stride_loc_bandgap
   END TYPE

   TYPE ri_basis_opt
      REAL(KIND=dp)            :: DI_rel
      REAL(KIND=dp)            :: DRI
      REAL(KIND=dp)            :: eps_step
      INTEGER                  :: max_num_iter
      INTEGER                  :: basis_quality
      INTEGER, DIMENSION(:), ALLOCATABLE :: RI_nset_per_l
   END TYPE

   TYPE grad_util
      TYPE(two_dim_real_array), DIMENSION(2) :: P_ij, P_ab
      TYPE(three_dim_real_array), DIMENSION(2) :: Gamma_P_ia
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: operator_half, PQ_half, Gamma_PQ, Gamma_PQ_2
      TYPE(dbcsr_p_type), DIMENSION(:, :), ALLOCATABLE :: G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:), ALLOCATABLE :: mo_coeff_o, mo_coeff_v
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: P_mo, W_mo, L_jb
      REAL(KIND=dp) :: cphf_eps_conv, scale_step_size
      INTEGER :: cphf_max_num_iter, z_solver_method, cphf_restart
      LOGICAL :: enforce_decrease, recalc_residual, polak_ribiere
      TYPE(qs_p_env_type), POINTER :: p_env => NULL()
      TYPE(qs_force_type), DIMENSION(:), POINTER :: mp2_force => NULL()
      REAL(KIND=dp), DIMENSION(3, 3) :: mp2_virial
      REAL(dp)                 :: eps_canonical
      LOGICAL                  :: free_hfx_buffer
      LOGICAL                  :: use_old_grad
      INTEGER :: dot_blksize
      INTEGER :: max_parallel_comm
   END TYPE

   TYPE mp2_type
      INTEGER                  :: method
      TYPE(mp2_laplace_type)   :: ri_laplace
      TYPE(mp2_direct_type)    :: direct_canonical
      TYPE(libint_potential_type) :: potential_parameter
      TYPE(mp2_gpw_type)       :: mp2_gpw
      TYPE(ri_mp2_type)        :: ri_mp2
      TYPE(ri_rpa_type)        :: ri_rpa
      TYPE(ri_rpa_im_time_type) &
         :: ri_rpa_im_time
      TYPE(ri_g0w0_type)       :: ri_g0w0
      TYPE(ri_basis_opt)       :: ri_opt_param
      TYPE(grad_util)          :: ri_grad
      REAL(dp) :: mp2_memory
      REAL(dp) :: scale_S
      REAL(dp) :: scale_T
      INTEGER  :: mp2_num_proc
      INTEGER  :: block_size_row
      INTEGER  :: block_size_col
      LOGICAL  :: calc_PQ_cond_num
      LOGICAL  :: hf_fail
      LOGICAL  :: p_screen
      LOGICAL  :: not_last_hfx
      LOGICAL  :: do_im_time
      INTEGER  :: eri_method
      TYPE(cp_eri_mme_param), POINTER  :: eri_mme_param
      INTEGER, DIMENSION(:), POINTER  :: eri_blksize => NULL()
      LOGICAL  :: do_svd
      REAL(KIND=dp) :: eps_range
      TYPE(libint_potential_type) :: ri_metric
      TYPE(C_ptr)                 :: local_gemm_ctx
      REAL(dp) :: e_gap, e_range
   END TYPE

   TYPE integ_mat_buffer_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: msg
      INTEGER, DIMENSION(:), ALLOCATABLE  :: sizes
      INTEGER, DIMENSION(:, :), ALLOCATABLE  :: indx
      INTEGER :: proc
      TYPE(mp_request_type) :: msg_req
   END TYPE

   TYPE integ_mat_buffer_type_2D
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE  :: msg
      INTEGER :: proc
      TYPE(mp_request_type) :: msg_req
   END TYPE

   TYPE pair_list_type_mp2
      TYPE(pair_list_element_type), DIMENSION(:), ALLOCATABLE :: elements
      INTEGER :: n_element
   END TYPE pair_list_type_mp2

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE mp2_env_release(mp2_env)
      TYPE(mp2_type)                                     :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mp2_env_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! release the HFX section for the EXX calculation
      IF (.NOT. mp2_env%ri_rpa%reuse_hfx) THEN
         IF (ASSOCIATED(mp2_env%ri_rpa%x_data)) CALL hfx_release(mp2_env%ri_rpa%x_data)
      END IF
      IF (ASSOCIATED(mp2_env%ri_rpa%xc_section_aux)) CALL section_vals_release(mp2_env%ri_rpa%xc_section_aux)
      IF (ASSOCIATED(mp2_env%ri_rpa%xc_section_primary)) CALL section_vals_release(mp2_env%ri_rpa%xc_section_primary)

      IF (mp2_env%eri_method .EQ. do_eri_mme) CALL cp_eri_mme_finalize(mp2_env%eri_mme_param)
      IF (ASSOCIATED(mp2_env%eri_mme_param)) DEALLOCATE (mp2_env%eri_mme_param)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%tau_tj)) DEALLOCATE (mp2_env%ri_rpa_im_time%tau_tj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%tau_wj)) DEALLOCATE (mp2_env%ri_rpa_im_time%tau_wj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%tj)) DEALLOCATE (mp2_env%ri_rpa_im_time%tj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%wj)) DEALLOCATE (mp2_env%ri_rpa_im_time%wj)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w)) DEALLOCATE (mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w)
      IF (ASSOCIATED(mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t)) DEALLOCATE (mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t)

      CALL local_gemm_destroy(mp2_env%local_gemm_ctx)

      CALL timestop(handle)

   END SUBROUTINE mp2_env_release

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE mp2_env_create(mp2_env)
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mp2_env_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(.NOT. ASSOCIATED(mp2_env))

      ALLOCATE (mp2_env)

      ! these two functions are empty if spla is build without gpu support and
      ! OFFLOAD_GEMM is not given at compilation time

      CALL local_gemm_create(mp2_env%local_gemm_ctx, LOCAL_GEMM_PU_GPU)
      CALL local_gemm_set_op_threshold_gpu(mp2_env%local_gemm_ctx, 128*128*128*2)

      NULLIFY (mp2_env%ri_rpa%x_data)

      CALL timestop(handle)

   END SUBROUTINE mp2_env_create

END MODULE mp2_types
