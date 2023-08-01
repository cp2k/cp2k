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
   USE input_constants,                 ONLY: &
        do_eri_mme, eri_default, gw_pade_approx, kp_weights_W_auto, mp2_method_direct, &
        mp2_method_gpw, mp2_method_none, mp2_ri_optimize_basis, ri_mp2_laplace, ri_mp2_method_gpw, &
        ri_rpa_g0w0_crossing_z_shot, ri_rpa_method_gpw, soc_none, wfc_mm_style_gemm
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_type
   USE iso_c_binding,                   ONLY: c_null_ptr,&
                                              c_ptr
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
      INTEGER       :: n_quadrature = -1
      INTEGER       :: num_integ_groups = -1
   END TYPE

   TYPE mp2_direct_type
      LOGICAL  :: big_send = .FALSE.
   END TYPE

   TYPE mp2_gpw_type
      REAL(KIND=dp)            :: eps_grid = 0.0_dp, eps_filter = 0.0_dp, eps_pgf_orb_S = 0.0_dp
      INTEGER                  :: print_level = 0
      REAL(KIND=dp)            :: cutoff = 0.0_dp
      REAL(KIND=dp)            :: relative_cutoff = 0.0_dp
      INTEGER                  :: size_lattice_sum = 0
   END TYPE mp2_gpw_type

   TYPE ri_mp2_type
      INTEGER                  :: block_size = 0, number_integration_groups = 0
      LOGICAL                  :: print_dgemm_info = .FALSE.
   END TYPE

   TYPE ri_rpa_type
      INTEGER                  :: rpa_num_quad_points = -1
      INTEGER                  :: rpa_num_integ_groups = -1
      INTEGER                  :: mm_style = wfc_mm_style_gemm
      TYPE(hfx_type), DIMENSION(:, :), POINTER &
         :: x_data => NULL()
      TYPE(section_vals_type), POINTER         :: xc_section_primary => Null(), &
                                                  xc_section_aux => Null()
      LOGICAL                  :: reuse_hfx = .FALSE.
      LOGICAL                  :: minimax_quad = .FALSE.
      LOGICAL                  :: do_ri_g0w0 = .FALSE.
      LOGICAL                  :: do_admm = .FALSE.
      LOGICAL                  :: do_ri_axk = .FALSE.
      LOGICAL                  :: do_rse = .FALSE.
      LOGICAL                  :: print_dgemm_info = .FALSE.
      TYPE(dbcsr_type), POINTER             :: mo_coeff_o => NULL(), &
                                               mo_coeff_v => NULL()
      REAL(KIND=dp)            :: ener_axk = 0.0_dp
      REAL(KIND=dp)            :: rse_corr_diag = 0.0_dp
      REAL(KIND=dp)            :: rse_corr = 0.0_dp
      REAL(KIND=dp)            :: scale_rpa = 0.0_dp
   END TYPE

   TYPE ri_rpa_im_time_type
      INTEGER                  :: cut_memory = 0
      LOGICAL                  :: memory_info = .FALSE., make_chi_pos_definite = .FALSE., &
                                  make_overlap_mat_ao_pos_definite = .FALSE., &
                                  trunc_coulomb_ri_x = .FALSE., keep_quad = .FALSE., &
                                  do_kpoints_from_Gamma = .FALSE., do_extrapolate_kpoints = .FALSE.
      REAL(KIND=dp)            :: eps_filter = 0.0_dp, &
                     eps_filter_factor = 0.0_dp, eps_compress = 0.0_dp, exp_tailored_weights = 0.0_dp, regularization_RI = 0.0_dp, &
                                  eps_eigval_S = 0.0_dp, eps_eigval_S_Gamma = 0.0_dp, rel_cutoff_trunc_coulomb_ri_x = 0.0_dp
      REAL(KIND=dp), DIMENSION(:), POINTER :: tau_tj => NULL(), tau_wj => NULL(), tj => NULL(), wj => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: weights_cos_tf_t_to_w => NULL(), weights_cos_tf_w_to_t => NULL()
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: Eigenval_Gamma, wkp_V
   INTEGER                  :: group_size_P = 0, group_size_3c = 0, kpoint_weights_W_method = kp_weights_W_auto, k_mesh_g_factor = 0
      INTEGER, DIMENSION(:), POINTER     :: kp_grid => NULL()
      INTEGER, DIMENSION(3) ::  kp_grid_extra = -1
      LOGICAL                  :: do_im_time_kpoints = .FALSE.
      INTEGER                  :: min_bsize = 0, min_bsize_mo = 0, nkp_orig = 0, nkp_extra = 0
      TYPE(kpoint_type), POINTER :: kpoints_G => NULL(), kpoints_Sigma => NULL(), kpoints_Sigma_no_xc => NULL()
      INTEGER, ALLOCATABLE, DIMENSION(:)      :: starts_array_mc_RI, ends_array_mc_RI, &
                                                 starts_array_mc_block_RI, &
                                                 ends_array_mc_block_RI, &
                                                 starts_array_mc, ends_array_mc, &
                                                 starts_array_mc_block, &
                                                 ends_array_mc_block

   END TYPE

   TYPE ri_g0w0_type
      INTEGER                  :: corr_mos_occ = 0
      INTEGER                  :: corr_mos_virt = 0
      INTEGER                  :: corr_mos_occ_beta = 0
      INTEGER                  :: corr_mos_virt_beta = 0
      INTEGER                  :: num_poles = 0
      INTEGER                  :: nparam_pade = 0
      INTEGER                  :: analytic_continuation = gw_pade_approx
      REAL(KIND=dp)            :: omega_max_fit = 0.0_dp
      INTEGER                  :: crossing_search = ri_rpa_g0w0_crossing_z_shot
      REAL(KIND=dp)            :: fermi_level_offset = 0.0_dp
      INTEGER                  :: iter_evGW = 0, iter_sc_GW0 = 0
      REAL(KIND=dp)            :: eps_iter = 0.0_dp
      LOGICAL                  :: do_ri_Sigma_x = .FALSE., &
                                  do_periodic = .FALSE., &
                                  print_self_energy = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: vec_Sigma_x_minus_vxc_gw
      INTEGER, DIMENSION(:), POINTER    :: kp_grid => NULL(), kp_grid_Sigma => NULL()
      INTEGER                  :: num_kp_grids = 0
      REAL(KIND=dp)            :: eps_kpoint = 0.0_dp
      LOGICAL                  :: do_mo_coeff_gamma = .FALSE., do_average_deg_levels = .FALSE.
      REAL(KIND=dp)            :: eps_eigenval = 0.0_dp
      LOGICAL                  :: do_extra_kpoints = .FALSE., do_aux_bas_gw = .FALSE.
      REAL(KIND=dp)            :: frac_aux_mos = 0.0_dp
      INTEGER                  :: num_omega_points = 0
      LOGICAL                  :: do_bse = .FALSE.
      INTEGER                  :: num_z_vectors = 0, max_iter_bse = 0
      REAL(KIND=dp)            :: eps_min_trans = 0.0_dp
      LOGICAL                  :: do_ic_model = .FALSE., print_ic_values = .FALSE.
      REAL(KIND=dp)            :: eps_dist = 0.0_dp
      TYPE(one_dim_real_array), DIMENSION(2) :: ic_corr_list
      INTEGER :: print_exx = 0
      LOGICAL :: do_gamma_only_sigma = .FALSE.
      LOGICAL :: update_xc_energy = .FALSE., do_kpoints_Sigma = .FALSE., print_local_bandgap = .FALSE.
      INTEGER :: n_kp_in_kp_line = 0, n_special_kp = 0, nkp_self_energy = 0, &
                 nkp_self_energy_special_kp = 0, nkp_self_energy_monkh_pack = 0, &
                 soc_type = soc_none
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: xkp_special_kp
      TYPE(dbcsr_p_type), DIMENSION(:), ALLOCATABLE :: &
         matrix_sigma_x_minus_vxc, matrix_ks
      REAL(KIND=dp)            :: broadening_print_loc_bandgap = 0.0_dp, energy_window_print_loc_bandgap = 0.0_dp, &
                                  ldos_thresh_print_loc_bandgap = 0.0_dp, energy_spacing_print_loc_bandgap = 0.0_dp, &
                                  regularization_minimax = 0.0_dp, soc_energy_window = 0.0_dp
      INTEGER, DIMENSION(:), POINTER :: stride_loc_bandgap => NULL()
   END TYPE

   TYPE ri_basis_opt
      REAL(KIND=dp)            :: DI_rel = 0.0_dp
      REAL(KIND=dp)            :: DRI = 0.0_dp
      REAL(KIND=dp)            :: eps_step = 0.0_dp
      INTEGER                  :: max_num_iter = 0
      INTEGER                  :: basis_quality = 0
      INTEGER, DIMENSION(:), ALLOCATABLE :: RI_nset_per_l
   END TYPE

   TYPE grad_util
      TYPE(two_dim_real_array), DIMENSION(2) :: P_ij, P_ab
      TYPE(three_dim_real_array), DIMENSION(2) :: Gamma_P_ia
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: operator_half, PQ_half, Gamma_PQ, Gamma_PQ_2
      TYPE(dbcsr_p_type), DIMENSION(:, :), ALLOCATABLE :: G_P_ia
      TYPE(dbcsr_p_type), DIMENSION(:), ALLOCATABLE :: mo_coeff_o, mo_coeff_v
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:) :: P_mo, W_mo, L_jb
      REAL(KIND=dp) :: cphf_eps_conv = 0.0_dp, scale_step_size = 0.0_dp
      INTEGER :: cphf_max_num_iter = 0, z_solver_method = 0, cphf_restart = 0
      LOGICAL :: enforce_decrease = .FALSE., recalc_residual = .FALSE., polak_ribiere = .FALSE.
      TYPE(qs_p_env_type), POINTER :: p_env => NULL()
      TYPE(qs_force_type), DIMENSION(:), POINTER :: mp2_force => NULL()
      REAL(KIND=dp), DIMENSION(3, 3) :: mp2_virial = 0.0_dp
      REAL(dp)                 :: eps_canonical = 0.0_dp
      LOGICAL                  :: free_hfx_buffer = .FALSE.
      LOGICAL                  :: use_old_grad = .FALSE.
      INTEGER :: dot_blksize = 0
      INTEGER :: max_parallel_comm = 0
   END TYPE

   TYPE mp2_type
      INTEGER                  :: method = mp2_method_none
      TYPE(mp2_laplace_type)   :: ri_laplace = mp2_laplace_type()
      TYPE(mp2_direct_type)    :: direct_canonical = mp2_direct_type()
      TYPE(libint_potential_type) :: potential_parameter = libint_potential_type()
      TYPE(mp2_gpw_type)       :: mp2_gpw = mp2_gpw_type()
      TYPE(ri_mp2_type)        :: ri_mp2 = ri_mp2_type()
      TYPE(ri_rpa_type)        :: ri_rpa = ri_rpa_type()
      TYPE(ri_rpa_im_time_type) &
         :: ri_rpa_im_time
      TYPE(ri_g0w0_type)       :: ri_g0w0
      TYPE(ri_basis_opt)       :: ri_opt_param
      TYPE(grad_util)          :: ri_grad
      REAL(dp) :: mp2_memory = 0.0_dp
      REAL(dp) :: scale_S = 0.0_dp
      REAL(dp) :: scale_T = 0.0_dp
      INTEGER  :: mp2_num_proc = 0
      INTEGER  :: block_size_row = 0
      INTEGER  :: block_size_col = 0
      LOGICAL  :: calc_PQ_cond_num = .FALSE.
      LOGICAL  :: hf_fail = .FALSE.
      LOGICAL  :: p_screen = .FALSE.
      LOGICAL  :: not_last_hfx = .FALSE.
      LOGICAL  :: do_im_time = .FALSE.
      INTEGER  :: eri_method = eri_default
      TYPE(cp_eri_mme_param), POINTER  :: eri_mme_param => NULL()
      INTEGER, DIMENSION(:), POINTER  :: eri_blksize => NULL()
      LOGICAL  :: do_svd = .FALSE.
      REAL(KIND=dp) :: eps_range = 0.0_dp
      TYPE(libint_potential_type) :: ri_metric = libint_potential_type()
      TYPE(C_ptr)                 :: local_gemm_ctx = C_NULL_PTR
      REAL(dp) :: e_gap = 0.0_dp, e_range = 0.0_dp
   END TYPE

   TYPE integ_mat_buffer_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE  :: msg
      INTEGER, DIMENSION(:), ALLOCATABLE  :: sizes
      INTEGER, DIMENSION(:, :), ALLOCATABLE  :: indx
      INTEGER :: proc = -1
      TYPE(mp_request_type) :: msg_req = mp_request_type()
   END TYPE

   TYPE integ_mat_buffer_type_2D
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE  :: msg
      INTEGER :: proc = -1
      TYPE(mp_request_type) :: msg_req = mp_request_type()
   END TYPE

   TYPE pair_list_type_mp2
      TYPE(pair_list_element_type), DIMENSION(:), ALLOCATABLE :: elements
      INTEGER :: n_element = 0
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
