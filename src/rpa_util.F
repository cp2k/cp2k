!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility functions for RPA calculations
!> \par History
!>      06.2019 Moved from rpa_ri_gpw.F [Frederick Stein]
! **************************************************************************************************
MODULE rpa_util
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_syrk,&
                                              cp_fm_transpose
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_filter, dbcsr_get_info, &
        dbcsr_get_occupation, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_release_p, dbcsr_set, dbcsr_type
   USE dbcsr_tensor_api,                ONLY: dbcsr_t_destroy,&
                                              dbcsr_t_type
   USE input_constants,                 ONLY: wfc_mm_style_gemm,&
                                              wfc_mm_style_syrk
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: z_zero
   USE message_passing,                 ONLY: mp_comm_split_direct,&
                                              mp_sum
   USE mp2_laplace,                     ONLY: calc_fm_mat_S_laplace
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type,&
                                              two_dim_int_array
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_communication,               ONLY: fm_redistribute
   USE rpa_gw_kpoints,                  ONLY: compute_wkp_W
   USE rpa_im_time,                     ONLY: setup_mat_for_mem_cut_3c
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_util'

   PUBLIC :: RPA_postprocessing_nokp, alloc_im_time, calc_mat_Q, RPA_postprocessing_start, get_mat_3c_overl_int_cut, &
             dealloc_im_time, contract_P_omega_with_mat_L

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param num_integ_points ...
!> \param fm_mat_Q ...
!> \param do_mao ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_matrix_L_RI_metric ...
!> \param mat_munu ...
!> \param mat_dm_occ_local ...
!> \param mat_dm_virt_local ...
!> \param mat_P_global ...
!> \param mat_M ...
!> \param mat_3c_overl_int ...
!> \param mat_3c_overl_int_mao_for_occ ...
!> \param mat_3c_overl_int_mao_for_virt ...
!> \param do_dbcsr_t ...
!> \param t_3c_O ...
!> \param matrix_s ...
!> \param kpoints ...
!> \param eps_filter ...
!> \param eps_filter_im_time ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param cut_RI ...
!> \param cut_memory ...
!> \param nkp ...
!> \param num_cells_dm ...
!> \param num_3c_repl ...
!> \param size_P ...
!> \param n_group_col ...
!> \param n_group_row ...
!> \param ikp_local ...
!> \param mepos_P_from_RI_row ...
!> \param my_group_L_sizes_im_time ...
!> \param my_group_L_starts_im_time ...
!> \param row_from_LLL ...
!> \param ends_array_prim_col ...
!> \param ends_array_prim_fullcol ...
!> \param ends_array_prim_fullrow ...
!> \param ends_array_prim_row ...
!> \param index_to_cell_3c ...
!> \param sizes_array_prim_col ...
!> \param sizes_array_prim_row ...
!> \param starts_array_prim_col ...
!> \param cell_to_index_3c ...
!> \param non_zero_blocks_3c ...
!> \param starts_array_prim_fullrow ...
!> \param starts_array_prim_row ...
!> \param starts_array_prim_fullcol ...
!> \param col_blk_size ...
!> \param ends_array_cm ...
!> \param ends_array_cm_mao_occ ...
!> \param ends_array_cm_mao_virt ...
!> \param row_blk_offset ...
!> \param row_blk_size ...
!> \param starts_array_cm ...
!> \param starts_array_cm_mao_occ ...
!> \param starts_array_cm_mao_virt ...
!> \param do_ic_model ...
!> \param do_kpoints_cubic_RPA ...
!> \param do_kpoints_from_Gamma ...
!> \param do_ri_Sigma_x ...
!> \param my_open_shell ...
!> \param cycle_due_to_sparse_dm ...
!> \param multiply_needed_occ ...
!> \param multiply_needed_virt ...
!> \param has_mat_P_blocks ...
!> \param buffer_mat_M ...
!> \param wkp_W ...
!> \param cfm_mat_Q ...
!> \param fm_mat_L ...
!> \param fm_mat_RI_global_work ...
!> \param fm_mat_work ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param mat_dm ...
!> \param mat_L ...
!> \param mat_M_P_munu_occ ...
!> \param mat_M_P_munu_virt ...
!> \param mat_P_global_copy ...
!> \param mat_SinvVSinv ...
!> \param mat_M_mu_Pnu_occ ...
!> \param mat_M_mu_Pnu_virt ...
!> \param mat_P_omega ...
!> \param mat_P_omega_kp ...
!> \param mat_dm_loc_occ_cut ...
!> \param mat_dm_loc_virt_cut ...
!> \param mat_dm_loc_occ ...
!> \param mat_dm_loc_virt ...
!> \param mat_work ...
!> \param offset_combi_block ...
!> \param mat_P_omega_beta ...
! **************************************************************************************************
   SUBROUTINE alloc_im_time(qs_env, mp2_env, para_env, para_env_sub, dimen_RI, dimen_RI_red, num_integ_points, &
                            fm_mat_Q, do_mao, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                            fm_matrix_L_RI_metric, mat_munu, mat_dm_occ_local, mat_dm_virt_local, mat_P_global, mat_M, &
                            mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, do_dbcsr_t, &
                            t_3c_O, matrix_s, kpoints, eps_filter, eps_filter_im_time, &
                            do_ri_sos_laplace_mp2, cut_RI, cut_memory, nkp, num_cells_dm, num_3c_repl, &
                            size_P, n_group_col, n_group_row, ikp_local, mepos_P_from_RI_row, &
                            my_group_L_sizes_im_time, my_group_L_starts_im_time, row_from_LLL, ends_array_prim_col, &
                            ends_array_prim_fullcol, ends_array_prim_fullrow, ends_array_prim_row, index_to_cell_3c, &
                            sizes_array_prim_col, sizes_array_prim_row, starts_array_prim_col, cell_to_index_3c, &
                            non_zero_blocks_3c, starts_array_prim_fullrow, starts_array_prim_row, &
                            starts_array_prim_fullcol, col_blk_size, ends_array_cm, ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                            row_blk_offset, row_blk_size, starts_array_cm, starts_array_cm_mao_occ, &
                            starts_array_cm_mao_virt, do_ic_model, do_kpoints_cubic_RPA, &
                            do_kpoints_from_Gamma, do_ri_Sigma_x, my_open_shell, cycle_due_to_sparse_dm, multiply_needed_occ, &
                            multiply_needed_virt, has_mat_P_blocks, buffer_mat_M, wkp_W, &
                            cfm_mat_Q, fm_mat_L, fm_mat_RI_global_work, fm_mat_work, fm_mo_coeff_occ_scaled, &
                            fm_mo_coeff_virt_scaled, mat_dm, mat_L, mat_M_P_munu_occ, mat_M_P_munu_virt, mat_P_global_copy, &
                            mat_SinvVSinv, mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, mat_P_omega, mat_P_omega_kp, mat_dm_loc_occ_cut, &
                            mat_dm_loc_virt_cut, mat_dm_loc_occ, mat_dm_loc_virt, mat_work, offset_combi_block, &
                            mat_P_omega_beta)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red, num_integ_points
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      LOGICAL, INTENT(IN)                                :: do_mao
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu, mat_dm_occ_local, &
                                                            mat_dm_virt_local, mat_P_global, mat_M
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      LOGICAL, INTENT(IN)                                :: do_dbcsr_t
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: t_3c_O
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, eps_filter_im_time
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2
      INTEGER, INTENT(IN)                                :: cut_RI, cut_memory
      INTEGER, INTENT(OUT)                               :: nkp, num_cells_dm, num_3c_repl, size_P, &
                                                            n_group_col, n_group_row
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ikp_local, mepos_P_from_RI_row, &
                                                            my_group_L_sizes_im_time, &
                                                            my_group_L_starts_im_time, row_from_LLL
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: ends_array_prim_col, &
         ends_array_prim_fullcol, ends_array_prim_fullrow, ends_array_prim_row, index_to_cell_3c, &
         sizes_array_prim_col, sizes_array_prim_row, starts_array_prim_col
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: cell_to_index_3c, non_zero_blocks_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: starts_array_prim_fullrow, &
                                                            starts_array_prim_row, &
                                                            starts_array_prim_fullcol
      INTEGER, DIMENSION(:), POINTER :: col_blk_size, ends_array_cm, ends_array_cm_mao_occ, &
         ends_array_cm_mao_virt, row_blk_offset, row_blk_size, starts_array_cm, &
         starts_array_cm_mao_occ, starts_array_cm_mao_virt
      LOGICAL, INTENT(IN)                                :: do_ic_model, do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma, do_ri_Sigma_x, &
                                                            my_open_shell
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: cycle_due_to_sparse_dm, &
                                                            multiply_needed_occ, &
                                                            multiply_needed_virt
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :, :), &
         INTENT(OUT)                                     :: has_mat_P_blocks
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: buffer_mat_M
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp_W
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_type), POINTER                          :: fm_mat_RI_global_work, fm_mat_work, &
                                                            fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_dm, mat_L, mat_M_P_munu_occ, &
                                                            mat_M_P_munu_virt, mat_P_global_copy, &
                                                            mat_SinvVSinv
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega, mat_P_omega_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_dm_loc_occ_cut, mat_dm_loc_virt_cut
      TYPE(dbcsr_type), POINTER                          :: mat_dm_loc_occ, mat_dm_loc_virt, mat_work
      TYPE(two_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: offset_combi_block
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_im_time', routineP = moduleN//':'//routineN

      INTEGER :: cell_grid_dm(3), col_start_local, color_sub_col, color_sub_row, first_ikp_local, &
         handle, i_cell, i_cut_RI, i_dim, i_mem, j_cell, j_mem, LLL, n_local_col, n_local_row, &
         nblkrows_total, periodic(3), row, row_start_local
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_sub_kp

      CALL timeset(routineN, handle)

      ALLOCATE (my_group_L_starts_im_time(cut_RI))
      my_group_L_starts_im_time(:) = mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time
      ALLOCATE (my_group_L_sizes_im_time(cut_RI))
      my_group_L_sizes_im_time(:) = mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time

      IF (.NOT. do_dbcsr_t) THEN
         num_3c_repl = SIZE(mat_3c_overl_int, 3)
      ELSE
         num_3c_repl = SIZE(t_3c_O, 2)
      ENDIF

      IF (.NOT. do_dbcsr_t) THEN

         DO i_cut_RI = 1, cut_RI

            DO i_cell = 1, num_3c_repl
               DO j_cell = 1, num_3c_repl

                  CALL dbcsr_filter(mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix, &
                                    eps_filter)

               END DO
            END DO

         END DO

         CALL get_non_zero_blocks_3c(mat_3c_overl_int, para_env_sub, cut_RI, non_zero_blocks_3c)

      ENDIF

      IF (do_kpoints_cubic_RPA) THEN
         ! we always use an odd number of image cells
         ! CAUTION: also at another point, cell_grid_dm is defined, these definitions have to be identical
         DO i_dim = 1, 3
            cell_grid_dm(i_dim) = (kpoints%nkp_grid(i_dim)/2)*2 - 1
         END DO
         num_cells_dm = cell_grid_dm(1)*cell_grid_dm(2)*cell_grid_dm(3)
         ALLOCATE (index_to_cell_3c(3, SIZE(kpoints%index_to_cell, 2)))
         CPASSERT(SIZE(kpoints%index_to_cell, 1) == 3)
         index_to_cell_3c(:, :) = kpoints%index_to_cell(:, :)
         ALLOCATE (cell_to_index_3c(LBOUND(kpoints%cell_to_index, 1):UBOUND(kpoints%cell_to_index, 1), &
                                    LBOUND(kpoints%cell_to_index, 2):UBOUND(kpoints%cell_to_index, 2), &
                                    LBOUND(kpoints%cell_to_index, 3):UBOUND(kpoints%cell_to_index, 3)))
         cell_to_index_3c(:, :, :) = kpoints%cell_to_index(:, :, :)

      ELSE
         ALLOCATE (index_to_cell_3c(3, 1))
         index_to_cell_3c(:, 1) = 0
         ALLOCATE (cell_to_index_3c(0:0, 0:0, 0:0))
         cell_to_index_3c(0, 0, 0) = 1
         num_cells_dm = 1
      END IF

      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN

         CALL get_sub_para_kp(fm_struct_sub_kp, para_env, kpoints%nkp, &
                              dimen_RI, ikp_local, first_ikp_local, do_kpoints_cubic_RPA)

         NULLIFY (cfm_mat_Q)
         CALL cp_cfm_create(cfm_mat_Q, fm_struct_sub_kp)
         CALL cp_cfm_set_all(cfm_mat_Q, z_zero)
      ELSE
         first_ikp_local = 1
      END IF

      IF (.NOT. do_dbcsr_t) THEN
         CALL allocate_dm_loc(mat_dm_loc_occ, mat_dm_occ_local, mat_dm_loc_occ_cut, cut_RI, cut_memory, num_cells_dm)

         CALL allocate_dm_loc(mat_dm_loc_virt, mat_dm_virt_local, mat_dm_loc_virt_cut, cut_RI, cut_memory, num_cells_dm)

         CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
         CALL dbcsr_filter(mat_munu%matrix, 1.0_dp)

         IF (.NOT. do_mao) THEN
            mat_3c_overl_int_mao_for_occ => mat_3c_overl_int
            mat_3c_overl_int_mao_for_virt => mat_3c_overl_int
         END IF

         CALL allocate_mat_M_P_munu(mat_M_P_munu_occ, mat_M, mat_M_mu_Pnu_occ, cut_RI, mat_3c_overl_int_mao_for_occ)

         CALL allocate_mat_M_P_munu(mat_M_P_munu_virt, mat_M, mat_M_mu_Pnu_virt, cut_RI, mat_3c_overl_int_mao_for_virt)

      ENDIF

      ! if we do kpoints, mat_P has a kpoint and mat_P_omega has the inted
      ! mat_P(tau, kpoint)
      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN

         NULLIFY (cell)
         CALL get_qs_env(qs_env, cell=cell)
         CALL get_cell(cell=cell, periodic=periodic)

         CALL get_kpoint_info(kpoints, nkp=nkp)
         ! compute k-point weights such that functions 1/k^2, 1/k and const function are
         ! integrated correctly
         CALL compute_wkp_W(wkp_W, kpoints, cell%hmat, cell%h_inv, &
                            qs_env%mp2_env%ri_rpa_im_time%exp_kpoints, periodic)
      ELSE
         nkp = 1
      END IF

      IF (do_kpoints_cubic_RPA) THEN
         size_P = MAX(num_cells_dm/2 + 1, nkp)
      ELSE IF (do_kpoints_from_Gamma) THEN
         size_P = MAX(3**(periodic(1) + periodic(2) + periodic(3)), nkp)
      ELSE
         size_P = 1
      END IF

      CALL alloc_mat_P_omega(mat_P_omega, num_integ_points, size_P, mat_P_global%matrix)

      IF (my_open_shell .AND. do_ri_sos_laplace_mp2) THEN
         CALL alloc_mat_P_omega(mat_P_omega_beta, num_integ_points, size_P, mat_P_global%matrix)
      END IF

      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
         CALL alloc_mat_P_omega(mat_P_omega_kp, 2, size_P, mat_P_global%matrix)
      END IF

      IF (.NOT. do_dbcsr_t) THEN
         NULLIFY (mat_P_global_copy%matrix)
         ALLOCATE (mat_P_global_copy%matrix)
         CALL dbcsr_create(mat_P_global_copy%matrix, template=mat_P_global%matrix)
         CALL dbcsr_copy(mat_P_global_copy%matrix, mat_P_global%matrix)

         n_group_row = mp2_env%ri_rpa_im_time_util(1)%n_group_row
         ALLOCATE (sizes_array_prim_row(0:n_group_row - 1, cut_memory))
         DO i_mem = 1, cut_memory
            sizes_array_prim_row(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(:)
         END DO
         ALLOCATE (starts_array_prim_row(0:n_group_row - 1, cut_memory))
         DO i_mem = 1, cut_memory
            starts_array_prim_row(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(:)
         END DO
         ALLOCATE (ends_array_prim_row(0:n_group_row - 1, cut_memory))
         DO i_mem = 1, cut_memory
            ends_array_prim_row(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(:)
         END DO
         ALLOCATE (starts_array_prim_fullrow(0:n_group_row - 1, cut_memory))
         DO i_mem = 1, cut_memory
            starts_array_prim_fullrow(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(:)
         END DO
         ALLOCATE (ends_array_prim_fullrow(0:n_group_row - 1, cut_memory))
         DO i_mem = 1, cut_memory
            ends_array_prim_fullrow(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(:)
         END DO

         n_group_col = mp2_env%ri_rpa_im_time_util(1)%n_group_col
         ALLOCATE (sizes_array_prim_col(0:n_group_col - 1, cut_memory))
         DO j_mem = 1, cut_memory
            sizes_array_prim_col(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(:)
         END DO
         ALLOCATE (starts_array_prim_col(0:n_group_col - 1, cut_memory))
         DO j_mem = 1, cut_memory
            starts_array_prim_col(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(:)
         END DO
         ALLOCATE (ends_array_prim_col(0:n_group_col - 1, cut_memory))
         DO j_mem = 1, cut_memory
            ends_array_prim_col(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(:)
         END DO

         ALLOCATE (starts_array_prim_fullcol(0:n_group_col - 1, cut_memory))
         DO j_mem = 1, cut_memory
            starts_array_prim_fullcol(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(:)
         END DO
         ALLOCATE (ends_array_prim_fullcol(0:n_group_col - 1, cut_memory))
         DO j_mem = 1, cut_memory
            ends_array_prim_fullcol(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(:)
         END DO

         ALLOCATE (offset_combi_block(cut_memory, cut_memory))

         color_sub_row = mp2_env%ri_rpa_im_time_util(1)%color_sub_row
         color_sub_col = mp2_env%ri_rpa_im_time_util(1)%color_sub_col

         DO i_mem = 1, cut_memory
            DO j_mem = 1, cut_memory

               n_local_row = sizes_array_prim_row(color_sub_row, i_mem)
               row_start_local = starts_array_prim_row(color_sub_row, i_mem)

               n_local_col = sizes_array_prim_col(color_sub_col, j_mem)
               col_start_local = starts_array_prim_col(color_sub_col, j_mem)

               ALLOCATE (offset_combi_block(i_mem, j_mem)%array(row_start_local:row_start_local + n_local_row - 1, &
                                                                col_start_local:col_start_local + n_local_col - 1))
               offset_combi_block(i_mem, j_mem)%array(:, :) = &
                  mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block(:, :)

            END DO
         END DO

         NULLIFY (starts_array_cm, ends_array_cm)
         starts_array_cm => mp2_env%ri_rpa_im_time%starts_array_cm
         ends_array_cm => mp2_env%ri_rpa_im_time%ends_array_cm

         NULLIFY (starts_array_cm_mao_occ, starts_array_cm_mao_virt, ends_array_cm_mao_occ, ends_array_cm_mao_virt)
         IF (do_mao) THEN
            starts_array_cm_mao_occ => mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ
            starts_array_cm_mao_virt => mp2_env%ri_rpa_im_time%starts_array_cm_mao_virt
            ends_array_cm_mao_occ => mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ
            ends_array_cm_mao_virt => mp2_env%ri_rpa_im_time%ends_array_cm_mao_virt
         ELSE
            starts_array_cm_mao_occ => starts_array_cm
            starts_array_cm_mao_virt => starts_array_cm
            ends_array_cm_mao_occ => ends_array_cm
            ends_array_cm_mao_virt => ends_array_cm
         END IF

         ! allocatable arrays for fast filling of dbcsr matrices
         CALL dbcsr_get_info(mat_M_P_munu_occ%matrix, row_blk_size=row_blk_size, &
                             col_blk_size=col_blk_size, nblkrows_total=nblkrows_total)
         ALLOCATE (buffer_mat_M(MAXVAL(row_blk_size), MAXVAL(col_blk_size)))

         ALLOCATE (mepos_P_from_RI_row(nblkrows_total))
         mepos_P_from_RI_row(:) = mp2_env%ri_rpa_im_time_util(1)%mepos_P_from_RI_row(:)

         ! allocatable arrays for fast filling of dbcsr matrices
         CALL dbcsr_get_info(mat_P_global%matrix, row_blk_size=row_blk_size, col_blk_size=col_blk_size)

         ALLOCATE (row_from_LLL(dimen_RI))
         row_from_LLL = 0

         CALL dbcsr_get_info(mat_M_P_munu_occ%matrix, &
                             nblkrows_total=nblkrows_total, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_size)

         DO LLL = 1, dimen_RI
            DO row = 1, nblkrows_total
               IF (row_blk_offset(row) <= LLL .AND. LLL < row_blk_offset(row) + row_blk_size(row)) THEN
                  row_from_LLL(LLL) = row
               END IF
            END DO
         END DO
      ENDIF

      NULLIFY (fm_mo_coeff_occ_scaled)
      CALL cp_fm_create(fm_mo_coeff_occ_scaled, fm_mo_coeff_occ%matrix_struct)
      CALL cp_fm_to_fm(fm_mo_coeff_occ, fm_mo_coeff_occ_scaled)
      CALL cp_fm_set_all(matrix=fm_mo_coeff_occ_scaled, alpha=0.0_dp)

      NULLIFY (fm_mo_coeff_virt_scaled)
      CALL cp_fm_create(fm_mo_coeff_virt_scaled, fm_mo_coeff_virt%matrix_struct)
      CALL cp_fm_to_fm(fm_mo_coeff_virt, fm_mo_coeff_virt_scaled)
      CALL cp_fm_set_all(matrix=fm_mo_coeff_virt_scaled, alpha=0.0_dp)

      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
         NULLIFY (fm_mat_RI_global_work)
         CALL cp_fm_create(fm_mat_RI_global_work, fm_matrix_L_RI_metric(1, 1)%matrix%matrix_struct)
         CALL cp_fm_to_fm(fm_matrix_L_RI_metric(1, 1)%matrix, fm_mat_RI_global_work)
         CALL cp_fm_set_all(fm_mat_RI_global_work, 0.0_dp)
      END IF

      IF (.NOT. do_dbcsr_t) THEN
         ALLOCATE (cycle_due_to_sparse_dm(cut_memory, cut_memory, num_integ_points))
         cycle_due_to_sparse_dm = .FALSE.

         ALLOCATE (multiply_needed_occ(cut_memory, cut_RI, num_cells_dm))
         multiply_needed_occ = .TRUE.

         ALLOCATE (multiply_needed_virt(cut_memory, cut_RI, num_cells_dm))
         multiply_needed_virt = .TRUE.
      ENDIF

      ALLOCATE (has_mat_P_blocks(num_cells_dm/2 + 1, cut_memory, cut_memory, num_3c_repl, num_3c_repl))
      has_mat_P_blocks = .TRUE.

      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
         CALL reorder_mat_L(fm_mat_L, fm_matrix_L_RI_metric, fm_mat_Q%matrix_struct, para_env, mat_L, &
                            mat_P_global%matrix, dimen_RI, dimen_RI_red, first_ikp_local, ikp_local, fm_struct_sub_kp)

         CALL cp_fm_struct_release(fm_struct_sub_kp)
      ELSE
         CALL reorder_mat_L(fm_mat_L, fm_matrix_L_RI_metric, fm_mat_Q%matrix_struct, para_env, mat_L, &
                            mat_P_global%matrix, dimen_RI, dimen_RI_red, first_ikp_local)
      END IF

      ! Create Scalapack working matrix for the contraction with the metric
      IF (dimen_RI == dimen_RI_red) THEN
         NULLIFY (fm_mat_work)
         CALL cp_fm_create(fm_mat_work, fm_mat_Q%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_Q, fm_mat_work)
         CALL cp_fm_set_all(matrix=fm_mat_work, alpha=0.0_dp)

      ELSE
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, template_fmstruct=fm_mat_Q%matrix_struct, &
                                  ncol_global=dimen_RI_red, nrow_global=dimen_RI)

         NULLIFY (fm_mat_work)
         CALL cp_fm_create(fm_mat_work, fm_struct)
         CALL cp_fm_set_all(matrix=fm_mat_work, alpha=0.0_dp)

         CALL cp_fm_struct_release(fm_struct)

      END IF

      ! Then its DBCSR counter part
      IF (.NOT. (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma)) THEN
         CALL dbcsr_get_info(mat_L%matrix, col_blk_size=col_blk_size, row_blk_size=row_blk_size)

         ! Create mat_work having the shape of the transposed of mat_L (compare with contract_P_omega_with_mat_L)
         NULLIFY (mat_work)
         ALLOCATE (mat_work)
         CALL dbcsr_create(mat_work, template=mat_L%matrix, row_blk_size=col_blk_size, col_blk_size=row_blk_size)
      END IF

      IF (do_ri_Sigma_x .OR. do_ic_model) THEN

         NULLIFY (mat_SinvVSinv%matrix)
         ALLOCATE (mat_SinvVSinv%matrix)
         CALL dbcsr_create(mat_SinvVSinv%matrix, template=mat_P_global%matrix)
         CALL dbcsr_set(mat_SinvVSinv%matrix, 0.0_dp)

         ! for kpoints we compute SinvVSinv later with kpoints
         IF (.NOT. do_kpoints_from_Gamma) THEN

            !  get the Coulomb matrix for Sigma_x = G*V
            CALL dbcsr_multiply("T", "N", 1.0_dp, mat_L%matrix, mat_L%matrix, &
                                0.0_dp, mat_SinvVSinv%matrix, filter_eps=eps_filter_im_time)

         END IF

      END IF

      IF (do_ri_Sigma_x) THEN

         NULLIFY (mat_dm%matrix)
         ALLOCATE (mat_dm%matrix)
         CALL dbcsr_create(mat_dm%matrix, template=matrix_s(1)%matrix)

      END IF

      CALL timestop(handle)

   END SUBROUTINE alloc_im_time

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_loc ...
!> \param mat_dm_local ...
!> \param mat_dm_loc_cut ...
!> \param cut_RI ...
!> \param cut_memory ...
!> \param num_cells_dm ...
! **************************************************************************************************
   SUBROUTINE allocate_dm_loc(mat_dm_loc, mat_dm_local, mat_dm_loc_cut, cut_RI, cut_memory, num_cells_dm)
      TYPE(dbcsr_type), POINTER                          :: mat_dm_loc
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_dm_local
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_dm_loc_cut
      INTEGER, INTENT(IN)                                :: cut_RI, cut_memory, num_cells_dm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_dm_loc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cell, i_cut_RI, i_mem

      CALL timeset(routineN, handle)

      NULLIFY (mat_dm_loc)
      CALL dbcsr_init_p(mat_dm_loc)
      CALL dbcsr_desymmetrize(mat_dm_local%matrix, mat_dm_loc)

      NULLIFY (mat_dm_loc_cut)
      CALL dbcsr_allocate_matrix_set(mat_dm_loc_cut, cut_RI, cut_memory, num_cells_dm)

      DO i_mem = 1, cut_memory
         DO i_cut_RI = 1, cut_RI
            DO i_cell = 1, num_cells_dm

               ALLOCATE (mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix)
               CALL dbcsr_create(matrix=mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix, &
                                 template=mat_dm_loc)

            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_dm_loc

! **************************************************************************************************
!> \brief ...
!> \param mat_M_P_munu ...
!> \param mat_M ...
!> \param mat_M_mu_Pnu ...
!> \param cut_RI ...
!> \param mat_3c_overl_int_mao ...
! **************************************************************************************************
   SUBROUTINE allocate_mat_M_P_munu(mat_M_P_munu, mat_M, mat_M_mu_Pnu, cut_RI, mat_3c_overl_int_mao)
      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_M_P_munu
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_M
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_M_mu_Pnu
      INTEGER, INTENT(IN)                                :: cut_RI
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int_mao

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_mat_M_P_munu', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cut_RI

      CALL timeset(routineN, handle)

      NULLIFY (mat_M_P_munu%matrix)
      ALLOCATE (mat_M_P_munu%matrix)
      CALL dbcsr_create(mat_M_P_munu%matrix, template=mat_M%matrix)

      NULLIFY (mat_M_mu_Pnu)
      CALL dbcsr_allocate_matrix_set(mat_M_mu_Pnu, cut_RI)
      DO i_cut_RI = 1, cut_RI
         ALLOCATE (mat_M_mu_Pnu(i_cut_RI)%matrix)
         CALL dbcsr_create(matrix=mat_M_mu_Pnu(i_cut_RI)%matrix, &
                           template=mat_3c_overl_int_mao(i_cut_RI, 1, 1)%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE allocate_mat_M_P_munu

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param num_integ_points ...
!> \param size_P ...
!> \param template ...
! **************************************************************************************************
   SUBROUTINE alloc_mat_P_omega(mat_P_omega, num_integ_points, size_P, template)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega
      INTEGER, INTENT(IN)                                :: num_integ_points, size_P
      TYPE(dbcsr_type), POINTER                          :: template

      CHARACTER(LEN=*), PARAMETER :: routineN = 'alloc_mat_P_omega', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_kp, jquad

      CALL timeset(routineN, handle)

      NULLIFY (mat_P_omega)
      CALL dbcsr_allocate_matrix_set(mat_P_omega, num_integ_points, size_P)
      DO jquad = 1, num_integ_points
         DO i_kp = 1, size_P
            ALLOCATE (mat_P_omega(jquad, i_kp)%matrix)
            CALL dbcsr_create(matrix=mat_P_omega(jquad, i_kp)%matrix, &
                              template=template)
            CALL dbcsr_set(mat_P_omega(jquad, i_kp)%matrix, 0.0_dp)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE alloc_mat_P_omega

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_L ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_struct_template ...
!> \param para_env ...
!> \param mat_L ...
!> \param mat_template ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param first_ikp_local ...
!> \param ikp_local ...
!> \param fm_struct_sub_kp ...
! **************************************************************************************************
   SUBROUTINE reorder_mat_L(fm_mat_L, fm_matrix_L_RI_metric, fm_struct_template, para_env, mat_L, mat_template, &
                            dimen_RI, dimen_RI_red, first_ikp_local, ikp_local, fm_struct_sub_kp)
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L, fm_matrix_L_RI_metric
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_template
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_L
      TYPE(dbcsr_type), INTENT(IN)                       :: mat_template
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red, first_ikp_local
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: ikp_local
      TYPE(cp_fm_struct_type), OPTIONAL, POINTER         :: fm_struct_sub_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reorder_mat_L', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_size, j_size, nblk
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_size
      LOGICAL                                            :: do_kpoints
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L_transposed, fmdummy

      CALL timeset(routineN, handle)

      do_kpoints = .FALSE.
      IF (PRESENT(ikp_local) .AND. PRESENT(fm_struct_sub_kp)) THEN
         do_kpoints = .TRUE.
      END IF

      ! Get the fm_struct for fm_mat_L
      NULLIFY (fm_struct)
      IF (dimen_RI == dimen_RI_red) THEN
         fm_struct => fm_struct_template
      ELSE
         ! The template is assumed to be square such that we need a new fm_struct if dimensions are not equal
         CALL cp_fm_struct_create(fm_struct, nrow_global=dimen_RI_red, ncol_global=dimen_RI, template_fmstruct=fm_struct_template)
      END IF

      ! Start to allocate the new full matrix
      ALLOCATE (fm_mat_L(SIZE(fm_matrix_L_RI_metric, 1), SIZE(fm_matrix_L_RI_metric, 2)))
      DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
         DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 2)
            IF (do_kpoints) THEN
               IF (ANY(ikp_local(:) == i_size)) THEN
                  CALL cp_fm_create(fm_mat_L(i_size, j_size)%matrix, fm_struct_sub_kp)
                  CALL cp_fm_set_all(fm_mat_L(i_size, j_size)%matrix, 0.0_dp)
               END IF
            ELSE
               CALL cp_fm_create(fm_mat_L(i_size, j_size)%matrix, fm_struct)
               CALL cp_fm_set_all(fm_mat_L(i_size, j_size)%matrix, 0.0_dp)
            END IF
         END DO
      END DO

      ! For the transposed matric we need a different fm_struct
      IF (dimen_RI == dimen_RI_red) THEN
         fm_struct => fm_mat_L(first_ikp_local, 1)%matrix%matrix_struct
      ELSE
         CALL cp_fm_struct_release(fm_struct)

         ! Create a fm_struct with transposed sizes
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, nrow_global=dimen_RI, ncol_global=dimen_RI_red, &
                                  template_fmstruct=fm_mat_L(first_ikp_local, 1)%matrix%matrix_struct) !, force_block=.TRUE.)
      END IF

      ! Allocate buffer matrix
      NULLIFY (fm_mat_L_transposed)
      CALL cp_fm_create(fm_mat_L_transposed, fm_struct)
      CALL cp_fm_set_all(matrix=fm_mat_L_transposed, alpha=0.0_dp)

      IF (dimen_RI /= dimen_RI_red) CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(fm_mat_L_transposed, context=blacs_env)

      ! For k-points copy matrices of your group
      ! Without kpoints, transpose matrix
      ! without kpoints, the size of fm_mat_L is 1x1. with kpoints, the size is N_kpoints x 2 (2 for real/complex)
      DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
      DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 2)
         IF (do_kpoints) THEN
            IF (ANY(ikp_local(:) == i_size)) THEN
               CALL cp_fm_copy_general(fm_matrix_L_RI_metric(i_size, j_size)%matrix, fm_mat_L_transposed, para_env)
               CALL cp_fm_to_fm(fm_mat_L_transposed, fm_mat_L(i_size, j_size)%matrix)
            ELSE
               NULLIFY (fmdummy)
               CALL cp_fm_copy_general(fm_matrix_L_RI_metric(i_size, j_size)%matrix, fmdummy, para_env)
            END IF
         ELSE
            CALL cp_fm_copy_general(fm_matrix_L_RI_metric(i_size, j_size)%matrix, fm_mat_L_transposed, blacs_env%para_env)
            CALL cp_fm_transpose(fm_mat_L_transposed, fm_mat_L(i_size, j_size)%matrix)
         END IF
      END DO
      END DO

      ! Release old matrix
      DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
      DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 2)
         CALL cp_fm_release(fm_matrix_L_RI_metric(i_size, j_size)%matrix)
      END DO
      END DO
      DEALLOCATE (fm_matrix_L_RI_metric)
      ! Release buffer
      CALL cp_fm_release(fm_mat_L_transposed)

      ! Create sparse variant of L
      NULLIFY (mat_L%matrix)
      ALLOCATE (mat_L%matrix)
      IF (dimen_RI == dimen_RI_red) THEN
         CALL dbcsr_create(mat_L%matrix, template=mat_template)
      ELSE
         CALL dbcsr_get_info(mat_template, nblkrows_total=nblk)

         CALL calculate_equal_blk_size(row_blk_size, dimen_RI_red, nblk)

         CALL dbcsr_create(mat_L%matrix, template=mat_template, row_blk_size=row_blk_size)

         DEALLOCATE (row_blk_size)
      END IF
      IF (.NOT. (do_kpoints)) THEN
         CALL copy_fm_to_dbcsr(fm_mat_L(1, 1)%matrix, mat_L%matrix)
      END IF

      CALL timestop(handle)

   END SUBROUTINE reorder_mat_L

! **************************************************************************************************
!> \brief ...
!> \param blk_size_new ...
!> \param dimen_RI_red ...
!> \param nblk ...
! **************************************************************************************************
   SUBROUTINE calculate_equal_blk_size(blk_size_new, dimen_RI_red, nblk)
      INTEGER, DIMENSION(:), POINTER                     :: blk_size_new
      INTEGER, INTENT(IN)                                :: dimen_RI_red, nblk

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_equal_blk_size', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col_per_blk, remainder

      NULLIFY (blk_size_new)
      ALLOCATE (blk_size_new(nblk))

      remainder = MOD(dimen_RI_red, nblk)
      col_per_blk = dimen_RI_red/nblk

      ! Determine a new distribution for the columns (corresponding to the number of columns)
      IF (remainder > 0) blk_size_new(1:remainder) = col_per_blk + 1
      blk_size_new(remainder + 1:nblk) = col_per_blk

   END SUBROUTINE calculate_equal_blk_size

! **************************************************************************************************
!> \brief ...
!> \param para_env_sub ...
!> \param do_mao ...
!> \param do_gw_im_time ...
!> \param mat_3c_overl_int ...
!> \param mat_3c_overl_int_mao_for_occ ...
!> \param mat_3c_overl_int_mao_for_virt ...
!> \param do_dbcsr_t ...
!> \param eps_filter ...
!> \param cut_RI ...
!> \param cut_memory ...
!> \param num_3c_repl ...
!> \param my_group_L_sizes_im_time ...
!> \param non_zero_blocks_3c_cut_col ...
!> \param ends_array_cm ...
!> \param ends_array_cm_mao_occ ...
!> \param ends_array_cm_mao_virt ...
!> \param starts_array_cm ...
!> \param starts_array_cm_mao_occ ...
!> \param starts_array_cm_mao_virt ...
!> \param do_kpoints_cubic_RPA ...
!> \param needed_cutRI_mem_R1vec_R2vec_for_kp ...
!> \param mat_3c_overl_int_cut ...
!> \param mat_3c_overl_int_mao_for_occ_cut ...
!> \param mat_3c_overl_int_mao_for_virt_cut ...
! **************************************************************************************************
   SUBROUTINE get_mat_3c_overl_int_cut(para_env_sub, do_mao, do_gw_im_time, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                                       mat_3c_overl_int_mao_for_virt, do_dbcsr_t, eps_filter, cut_RI, &
                                       cut_memory, num_3c_repl, my_group_L_sizes_im_time, &
                                       non_zero_blocks_3c_cut_col, ends_array_cm, ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                       starts_array_cm, starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                       do_kpoints_cubic_RPA, needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                       mat_3c_overl_int_cut, mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut)
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      LOGICAL, INTENT(IN)                                :: do_mao, do_gw_im_time
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      LOGICAL, INTENT(IN)                                :: do_dbcsr_t
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, INTENT(IN)                                :: cut_RI, cut_memory
      INTEGER, INTENT(OUT)                               :: num_3c_repl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_group_L_sizes_im_time
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :), &
         INTENT(OUT)                                     :: non_zero_blocks_3c_cut_col
      INTEGER, DIMENSION(:), POINTER :: ends_array_cm, ends_array_cm_mao_occ, &
         ends_array_cm_mao_virt, starts_array_cm, starts_array_cm_mao_occ, starts_array_cm_mao_virt
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(OUT) :: &
         needed_cutRI_mem_R1vec_R2vec_for_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut, &
         mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_3c_overl_int_cut', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! This splitting allows to keep allocation for gw im.time outside from remaining allocation routines retaining slim interfaces
      IF (.NOT. do_dbcsr_t) THEN
         CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_cut, mat_3c_overl_int, cut_memory, cut_RI, &
                                       starts_array_cm, ends_array_cm, my_group_L_sizes_im_time, eps_filter, &
                                       do_kpoints_cubic_RPA, do_gw_im_time, num_3c_repl)
      ENDIF

      IF (do_mao) THEN

         ! in setup_mat_for_mem_cut_3c, one deallocates mat_3c_overl_int, therefore for the beginning, deallocate
         ! the mao 3c overlap integrals here
         CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_occ, &
                                       cut_memory, cut_RI, starts_array_cm_mao_virt, ends_array_cm_mao_virt, &
                                       my_group_L_sizes_im_time, eps_filter, do_kpoints_cubic_RPA, do_gw_im_time, &
                                       num_3c_repl)
         CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_mao_for_virt_cut, mat_3c_overl_int_mao_for_virt, &
                                       cut_memory, cut_RI, starts_array_cm_mao_occ, ends_array_cm_mao_occ, &
                                       my_group_L_sizes_im_time, eps_filter, do_kpoints_cubic_RPA, do_gw_im_time, &
                                       num_3c_repl)

      ELSE
         IF (.NOT. do_dbcsr_t) THEN

            mat_3c_overl_int_mao_for_occ_cut => mat_3c_overl_int_cut
            mat_3c_overl_int_mao_for_virt_cut => mat_3c_overl_int_cut
         ENDIF
      END IF

      IF (.NOT. do_dbcsr_t) THEN
         CALL get_non_zero_blocks_3c_cut_col(mat_3c_overl_int_cut, para_env_sub, cut_RI, &
                                             cut_memory, non_zero_blocks_3c_cut_col)

         CALL check_sparsity_arrays_for_kp(needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                           mat_3c_overl_int_cut, cut_RI, cut_memory, para_env_sub, &
                                           do_kpoints_cubic_RPA)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE get_mat_3c_overl_int_cut

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param first_cycle ...
!> \param count_ev_sc_GW ...
!> \param virtual ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param homo ...
!> \param omega ...
!> \param omega_old ...
!> \param jquad ...
!> \param mm_style ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param alpha ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param para_env_RPA ...
!> \param do_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
!> \param RPA_proc_map ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param number_of_send ...
!> \param map_send_size ...
!> \param map_rec_size ...
!> \param local_size_source ...
!> \param my_num_dgemm_call ...
!> \param my_flop_rate ...
! **************************************************************************************************
   SUBROUTINE calc_mat_Q(fm_mat_S, do_ri_sos_laplace_mp2, first_cycle, count_ev_sc_GW, virtual, Eigenval, Eigenval_last, &
                         homo, omega, omega_old, jquad, mm_style, dimen_RI, dimen_ia, alpha, fm_mat_Q, fm_mat_Q_gemm, &
                         para_env_RPA, do_bse, fm_mat_Q_static_bse_gemm, RPA_proc_map, buffer_rec, buffer_send, number_of_send, &
                         map_send_size, map_rec_size, local_size_source, my_num_dgemm_call, my_flop_rate)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2, first_cycle
      INTEGER, INTENT(IN)                                :: count_ev_sc_GW, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval, Eigenval_last
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(IN)                          :: omega, omega_old
      INTEGER, INTENT(IN)                                :: jquad, mm_style, dimen_RI, dimen_ia
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, fm_mat_Q_gemm
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      LOGICAL, INTENT(IN)                                :: do_bse
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q_static_bse_gemm
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: RPA_proc_map
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: buffer_rec, buffer_send
      INTEGER, INTENT(IN)                                :: number_of_send
      INTEGER, DIMENSION(0:para_env_RPA%num_pe-1), &
         INTENT(IN)                                      :: map_send_size, map_rec_size
      INTEGER, DIMENSION(2, 0:para_env_RPA%num_pe-1), &
         INTENT(IN)                                      :: local_size_source
      INTEGER, INTENT(INOUT)                             :: my_num_dgemm_call
      REAL(KIND=dp), INTENT(INOUT)                       :: my_flop_rate

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mat_Q', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (do_ri_sos_laplace_mp2) THEN
         ! the first index of tau_tj starts with 0 (see mp2_weights)
         CALL calc_fm_mat_S_laplace(fm_mat_S, first_cycle, homo, virtual, Eigenval, omega, omega_old)

      ELSE
         CALL calc_fm_mat_S_rpa(fm_mat_S, first_cycle, count_ev_sc_GW, virtual, Eigenval, Eigenval_last, &
                                homo, omega, omega_old)
      END IF

      CALL contract_S_to_Q(mm_style, dimen_RI, dimen_ia, alpha, fm_mat_S, fm_mat_Q_gemm, para_env_RPA, &
                           my_num_dgemm_call, fm_mat_Q, RPA_proc_map, buffer_rec, buffer_send, number_of_send, &
                           map_send_size, map_rec_size, local_size_source, my_flop_rate)

      IF (do_bse .AND. jquad == 1) THEN
         CALL cp_fm_to_fm(fm_mat_Q_gemm, fm_mat_Q_static_bse_gemm)
      END IF
      CALL timestop(handle)

   END SUBROUTINE calc_mat_Q

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param first_cycle ...
!> \param count_ev_sc_GW ...
!> \param virtual ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param homo ...
!> \param omega ...
!> \param omega_old ...
! **************************************************************************************************
   SUBROUTINE calc_fm_mat_S_rpa(fm_mat_S, first_cycle, count_ev_sc_GW, virtual, Eigenval, Eigenval_last, homo, &
                                omega, omega_old)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: first_cycle
      INTEGER, INTENT(IN)                                :: count_ev_sc_GW, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval, Eigenval_last
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(IN)                          :: omega, omega_old

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_fm_mat_S_rpa', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: avirt, handle, i_global, iiB, iocc, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigen_diff

      CALL timeset(routineN, handle)

      ! get info of fm_mat_S
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! remove eigenvalue part of S matrix from the last eigenvalue self-c. cycle
      IF (first_cycle .AND. count_ev_sc_GW > 1) THEN
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval_last,fm_mat_S,virtual,homo,omega_old)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)

               iocc = MAX(1, i_global - 1)/virtual + 1
               avirt = i_global - (iocc - 1)*virtual
               eigen_diff = Eigenval_last(avirt + homo) - Eigenval_last(iocc)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)/ &
                                               SQRT(eigen_diff/(eigen_diff**2 + omega_old**2))

            END DO
         END DO

      END IF

      ! update G matrix with the new value of omega
      IF (first_cycle) THEN
         ! In this case just update the matrix (symmetric form) with
         ! SQRT((epsi_a-epsi_i)/((epsi_a-epsi_i)**2+omega**2))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,fm_mat_S,virtual,homo,omega)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)

               iocc = MAX(1, i_global - 1)/virtual + 1
               avirt = i_global - (iocc - 1)*virtual
               eigen_diff = Eigenval(avirt + homo) - Eigenval(iocc)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)* &
                                               SQRT(eigen_diff/(eigen_diff**2 + omega**2))

            END DO
         END DO
      ELSE
         ! In this case the update has to remove the old omega component thus
         ! SQRT(((epsi_a-epsi_i)**2+omega_old**2)/((epsi_a-epsi_i)**2+omega**2))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,&
!$OMP                    fm_mat_S,virtual,homo,omega,omega_old)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)

               iocc = MAX(1, i_global - 1)/virtual + 1
               avirt = i_global - (iocc - 1)*virtual
               eigen_diff = Eigenval(avirt + homo) - Eigenval(iocc)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)* &
                                               SQRT((eigen_diff**2 + omega_old**2)/(eigen_diff**2 + omega**2))

            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_fm_mat_S_rpa

! **************************************************************************************************
!> \brief ...
!> \param mm_style ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param alpha ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_gemm ...
!> \param para_env_RPA ...
!> \param my_num_dgemm_call ...
!> \param fm_mat_Q ...
!> \param RPA_proc_map ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param number_of_send ...
!> \param map_send_size ...
!> \param map_rec_size ...
!> \param local_size_source ...
!> \param my_flop_rate ...
! **************************************************************************************************
   SUBROUTINE contract_S_to_Q(mm_style, dimen_RI, dimen_ia, alpha, fm_mat_S, fm_mat_Q_gemm, para_env_RPA, my_num_dgemm_call, &
                              fm_mat_Q, RPA_proc_map, buffer_rec, buffer_send, number_of_send, &
                              map_send_size, map_rec_size, local_size_source, my_flop_rate)

      INTEGER, INTENT(IN)                                :: mm_style, dimen_RI, dimen_ia
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S, fm_mat_Q_gemm
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      INTEGER, INTENT(INOUT)                             :: my_num_dgemm_call
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: RPA_proc_map
      TYPE(integ_mat_buffer_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: buffer_rec, buffer_send
      INTEGER, INTENT(IN)                                :: number_of_send
      INTEGER, DIMENSION(0:para_env_RPA%num_pe-1), &
         INTENT(IN)                                      :: map_send_size, map_rec_size
      INTEGER, DIMENSION(2, 0:para_env_RPA%num_pe-1), &
         INTENT(IN)                                      :: local_size_source
      REAL(KIND=dp), INTENT(INOUT)                       :: my_flop_rate

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_S_to_Q', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: actual_flop_rate, t_end, t_start

      CALL timeset(routineN, handle)

      t_start = m_walltime()
      SELECT CASE (mm_style)
      CASE (wfc_mm_style_gemm)
         ! waste-fully computes the full symmetrix matrix, but maybe faster than cp_fm_syrk for optimized cp_fm_gemm
         CALL cp_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia, alpha=alpha, &
                      matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                      matrix_c=fm_mat_Q_gemm)
      CASE (wfc_mm_style_syrk)
         ! will only compute the upper half of the matrix, which is fine, since we only use it for cholesky later
         CALL cp_fm_syrk(uplo='U', trans='T', k=dimen_ia, alpha=alpha, matrix_a=fm_mat_S, &
                         ia=1, ja=1, beta=0.0_dp, matrix_c=fm_mat_Q_gemm)
      CASE DEFAULT
         CPABORT("")
      END SELECT
      t_end = m_walltime()
      actual_flop_rate = 2.0_dp*REAL(dimen_ia, KIND=dp)*dimen_RI*REAL(dimen_RI, KIND=dp)/(MAX(0.01_dp, t_end - t_start))
      IF (para_env_RPA%mepos == 0) my_flop_rate = my_flop_rate + actual_flop_rate
      my_num_dgemm_call = my_num_dgemm_call + 1

      ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
      CALL cp_fm_set_all(matrix=fm_mat_Q, alpha=0.0_dp)
      CALL fm_redistribute(fm_mat_Q_gemm, fm_mat_Q, RPA_proc_map, buffer_rec, buffer_send, &
                           number_of_send, map_send_size, map_rec_size, local_size_source, para_env_RPA)

      CALL timestop(handle)

   END SUBROUTINE contract_S_to_Q

! **************************************************************************************************
!> \brief ...
!> \param dimen_RI ...
!> \param trace_Qomega ...
!> \param fm_mat_Q ...
!> \param para_env_RPA ...
! **************************************************************************************************
   SUBROUTINE RPA_postprocessing_start(dimen_RI, trace_Qomega, fm_mat_Q, para_env_RPA)

      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), DIMENSION(dimen_RI), INTENT(OUT)    :: trace_Qomega
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'RPA_postprocessing_start', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! calculate the trace of Q and add 1 on the diagonal
      trace_Qomega = 0.0_dp
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                       SHARED(ncol_local,nrow_local,col_indices,row_indices,trace_Qomega,fm_mat_Q,dimen_RI)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               trace_Qomega(i_global) = fm_mat_Q%local_data(iiB, jjB)
               fm_mat_Q%local_data(iiB, jjB) = fm_mat_Q%local_data(iiB, jjB) + 1.0_dp
            END IF
         END DO
      END DO
      CALL mp_sum(trace_Qomega, para_env_RPA%group)

      CALL timestop(handle)

   END SUBROUTINE RPA_postprocessing_start

! **************************************************************************************************
!> \brief ...
!> \param dimen_RI ...
!> \param trace_Qomega ...
!> \param fm_mat_Q ...
!> \param para_env_RPA ...
!> \param Erpa ...
!> \param wjquad ...
! **************************************************************************************************
   SUBROUTINE RPA_postprocessing_nokp(dimen_RI, trace_Qomega, fm_mat_Q, para_env_RPA, Erpa, wjquad)

      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), DIMENSION(dimen_RI), INTENT(IN)     :: trace_Qomega
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      REAL(KIND=dp), INTENT(IN)                          :: wjquad

      CHARACTER(LEN=*), PARAMETER :: routineN = 'RPA_postprocessing_nokp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: FComega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Q_log

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! calculate Trace(Log(Matrix)) as Log(DET(Matrix)) via cholesky decomposition
      CALL cp_fm_cholesky_decompose(matrix=fm_mat_Q, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

      ALLOCATE (Q_log(dimen_RI))
      Q_log = 0.0_dp
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                         SHARED(ncol_local,nrow_local,col_indices,row_indices,Q_log,fm_mat_Q,dimen_RI)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               Q_log(i_global) = 2.0_dp*LOG(fm_mat_Q%local_data(iiB, jjB))
            END IF
         END DO
      END DO
      CALL mp_sum(Q_log, para_env_RPA%group)

      FComega = 0.0_dp
      DO iiB = 1, dimen_RI
         IF (MODULO(iiB, para_env_RPA%num_pe) /= para_env_RPA%mepos) CYCLE
         FComega = FComega + (Q_log(iiB) - trace_Qomega(iiB))/2.0_dp
      END DO
      Erpa = Erpa + FComega*wjquad

      DEALLOCATE (Q_log)

      CALL timestop(handle)

   END SUBROUTINE RPA_postprocessing_nokp

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int ...
!> \param para_env_sub ...
!> \param cut_RI ...
!> \param non_zero_blocks_3c ...
! **************************************************************************************************
   SUBROUTINE get_non_zero_blocks_3c(mat_3c_overl_int, para_env_sub, cut_RI, non_zero_blocks_3c)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, INTENT(IN)                                :: cut_RI
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: non_zero_blocks_3c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_non_zero_blocks_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, block_counter, col, handle, i_cell, i_cut_RI, iblk, imepos, j_cell, &
         maxlength, maxlength_tmp, nblkrows_total, num_cells_3c, row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c_tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      num_cells_3c = SIZE(mat_3c_overl_int, 2)
      CPASSERT(num_cells_3c == SIZE(mat_3c_overl_int, 3))

      CALL dbcsr_get_info(mat_3c_overl_int(1, 1, 1)%matrix, nblkrows_total=nblkrows_total)

      ALLOCATE (non_zero_blocks_3c_tmp(1:cut_RI, 1:nblkrows_total, 0:(para_env_sub%num_pe - 1)))
      non_zero_blocks_3c_tmp = 0

      DO i_cut_RI = 1, cut_RI

         DO i_cell = 1, num_cells_3c

            DO j_cell = 1, num_cells_3c

               CALL dbcsr_iterator_start(iter, mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix)
               DO WHILE (dbcsr_iterator_blocks_left(iter))
                  CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk)

                  non_zero_blocks_3c_tmp(i_cut_RI, row, para_env_sub%mepos) = 1

               ENDDO

               CALL dbcsr_iterator_stop(iter)

            END DO

         END DO

      END DO

      CALL mp_sum(non_zero_blocks_3c_tmp, para_env_sub%group)

      maxlength = 0

      DO imepos = 0, para_env_sub%num_pe - 1
         DO i_cut_RI = 1, cut_RI
            maxlength_tmp = 0
            DO iblk = 1, nblkrows_total
               IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                  maxlength_tmp = maxlength_tmp + 1
               END IF
            END DO
            IF (maxlength_tmp > maxlength) THEN
               maxlength = maxlength_tmp
            END IF
         END DO
      END DO

      ! save memory with smaller non_zero_blocks_3c
      ALLOCATE (non_zero_blocks_3c(1:cut_RI, 1:maxlength, 0:(para_env_sub%num_pe - 1)))
      non_zero_blocks_3c = 0

      DO imepos = 0, para_env_sub%num_pe - 1
         DO i_cut_RI = 1, cut_RI
            block_counter = 0
            DO iblk = 1, nblkrows_total
               IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                  block_counter = block_counter + 1
                  non_zero_blocks_3c(i_cut_RI, block_counter, imepos) = iblk
               END IF
            END DO
         END DO
      END DO

      DEALLOCATE (non_zero_blocks_3c_tmp)

      CALL timestop(handle)

   END SUBROUTINE get_non_zero_blocks_3c

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int_cut_col ...
!> \param para_env_sub ...
!> \param cut_RI ...
!> \param cut_memory ...
!> \param non_zero_blocks_3c_cut ...
! **************************************************************************************************
   SUBROUTINE get_non_zero_blocks_3c_cut_col(mat_3c_overl_int_cut_col, para_env_sub, cut_RI, cut_memory, &
                                             non_zero_blocks_3c_cut)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut_col
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, INTENT(IN)                                :: cut_RI, cut_memory
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :), &
         INTENT(OUT)                                     :: non_zero_blocks_3c_cut

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_non_zero_blocks_3c_cut_col', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, block_counter, col, handle, i_cell, i_cut_RI, i_mem, iblk, imepos, j_cell, &
         maxlength, maxlength_tmp, nblkrows_total, nblkrows_total_max, num_3c_repl, row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c_tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nblkrows_total_max = 0

      num_3c_repl = SIZE(mat_3c_overl_int_cut_col, 3)
      ! quickly check that replica of fourth component matches the one of third component
      CPASSERT(num_3c_repl == SIZE(mat_3c_overl_int_cut_col, 4))

      DO i_cut_RI = 1, cut_RI
         DO i_mem = 1, cut_memory
            CALL dbcsr_get_info(mat_3c_overl_int_cut_col(i_cut_RI, i_mem, 1, 1)%matrix, &
                                nblkrows_total=nblkrows_total)
            IF (nblkrows_total > nblkrows_total_max) THEN
               nblkrows_total_max = nblkrows_total
            END IF
         END DO
      END DO

      ALLOCATE (non_zero_blocks_3c_tmp(1:cut_RI, 1:nblkrows_total_max, 0:(para_env_sub%num_pe - 1)))
      non_zero_blocks_3c_tmp = 0

      maxlength = 0

      ! first, determine maxlength
      DO i_mem = 1, cut_memory

         DO i_cut_RI = 1, cut_RI

            DO i_cell = 1, num_3c_repl

               DO j_cell = 1, num_3c_repl

                  CALL dbcsr_iterator_start(iter, mat_3c_overl_int_cut_col(i_cut_RI, i_mem, i_cell, j_cell)%matrix)
                  DO WHILE (dbcsr_iterator_blocks_left(iter))
                     CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk)

                     non_zero_blocks_3c_tmp(i_cut_RI, col, para_env_sub%mepos) = 1

                  ENDDO

                  CALL dbcsr_iterator_stop(iter)

               END DO

            END DO

         END DO ! cut_RI

         CALL mp_sum(non_zero_blocks_3c_tmp, para_env_sub%group)

         maxlength_tmp = 0

         DO imepos = 0, para_env_sub%num_pe - 1
            DO i_cut_RI = 1, cut_RI
               maxlength_tmp = 0
               DO iblk = 1, nblkrows_total
                  IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                     maxlength_tmp = maxlength_tmp + 1
                  END IF
               END DO
               IF (maxlength_tmp > maxlength) THEN
                  maxlength = maxlength_tmp
               END IF
            END DO
         END DO

         non_zero_blocks_3c_tmp = 0

      END DO ! i_mem
      ! end determine maxlength

      ! save memory with a smaller non_zero_blocks_3c_cut
      ALLOCATE (non_zero_blocks_3c_cut(1:cut_RI, 1:maxlength, 0:(para_env_sub%num_pe - 1), 1:cut_memory))
      non_zero_blocks_3c_cut = 0

      ! now, fill non_zero_blocks_3c_cut
      DO i_mem = 1, cut_memory

         DO i_cut_RI = 1, cut_RI

            DO i_cell = 1, num_3c_repl

               DO j_cell = 1, num_3c_repl

                  CALL dbcsr_iterator_start(iter, mat_3c_overl_int_cut_col(i_cut_RI, i_mem, i_cell, j_cell)%matrix)
                  DO WHILE (dbcsr_iterator_blocks_left(iter))
                     CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk)

                     non_zero_blocks_3c_tmp(i_cut_RI, col, para_env_sub%mepos) = 1

                  ENDDO

                  CALL dbcsr_iterator_stop(iter)

               END DO
            END DO

         END DO ! cut_RI

         CALL mp_sum(non_zero_blocks_3c_tmp, para_env_sub%group)

         DO imepos = 0, para_env_sub%num_pe - 1
            DO i_cut_RI = 1, cut_RI
               block_counter = 0
               DO iblk = 1, nblkrows_total
                  IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                     block_counter = block_counter + 1
                     non_zero_blocks_3c_cut(i_cut_RI, block_counter, imepos, i_mem) = iblk
                  END IF
               END DO
            END DO
         END DO

         non_zero_blocks_3c_tmp = 0

      END DO ! i_mem
      ! end fill non_zero_blocks_3c_cut

      DEALLOCATE (non_zero_blocks_3c_tmp)

      CALL timestop(handle)

   END SUBROUTINE get_non_zero_blocks_3c_cut_col

! **************************************************************************************************
!> \brief ...
!> \param needed_cutRI_mem_R1vec_R2vec_for_kp ...
!> \param mat_3c_overl_int_cut ...
!> \param cut_RI ...
!> \param cut_memory ...
!> \param para_env_sub ...
!> \param do_kpoints_cubic_RPA ...
! **************************************************************************************************
   SUBROUTINE check_sparsity_arrays_for_kp(needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                           mat_3c_overl_int_cut, cut_RI, cut_memory, para_env_sub, &
                                           do_kpoints_cubic_RPA)

      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(OUT) :: &
         needed_cutRI_mem_R1vec_R2vec_for_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut
      INTEGER, INTENT(IN)                                :: cut_RI, cut_memory
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_sparsity_arrays_for_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cell, i_cut_RI, i_mem, j_cell, &
                                                            num_cells_3c
      REAL(KIND=dp)                                      :: occ_local, occ_local_sum_j_cell

      CALL timeset(routineN, handle)

      num_cells_3c = SIZE(mat_3c_overl_int_cut, 3)

      ALLOCATE (needed_cutRI_mem_R1vec_R2vec_for_kp(cut_RI, cut_memory, 0:num_cells_3c, 0:num_cells_3c))
      needed_cutRI_mem_R1vec_R2vec_for_kp(:, :, :, :) = .TRUE.

      IF (do_kpoints_cubic_RPA) THEN

         DO i_mem = 1, cut_memory

            DO i_cell = 1, num_cells_3c

               occ_local_sum_j_cell = 0.0_dp

               DO i_cut_RI = 1, cut_RI

                  DO j_cell = 1, num_cells_3c

                     occ_local = dbcsr_get_occupation(mat_3c_overl_int_cut(i_cut_RI, i_mem, &
                                                                           j_cell, i_cell)%matrix)

                     CALL mp_sum(occ_local, para_env_sub%group)

                     IF (occ_local < 1E-100_dp) THEN

                        needed_cutRI_mem_R1vec_R2vec_for_kp(i_cut_RI, i_mem, j_cell, i_cell) = .FALSE.

                     END IF

                     occ_local_sum_j_cell = occ_local_sum_j_cell + occ_local

                  END DO ! j_cell

               END DO ! i_cut_RI

            END DO ! i_cell
         END DO ! i_mem

      END IF

      CALL timestop(handle)

   END SUBROUTINE check_sparsity_arrays_for_kp

! **************************************************************************************************
!> \brief ...
!> \param fm_struct_sub_kp ...
!> \param para_env ...
!> \param nkp ...
!> \param dimen_RI ...
!> \param ikp_local ...
!> \param first_ikp_local ...
!> \param do_kpoints_cubic_RPA ...
! **************************************************************************************************
   SUBROUTINE get_sub_para_kp(fm_struct_sub_kp, para_env, nkp, dimen_RI, &
                              ikp_local, first_ikp_local, do_kpoints_cubic_RPA)
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_sub_kp
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nkp, dimen_RI
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: ikp_local
      INTEGER, INTENT(OUT)                               :: first_ikp_local
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA

      CHARACTER(len=*), PARAMETER :: routineN = 'get_sub_para_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: color_sub_kp, comm_sub_kp, handle, ikp, &
                                                            num_proc_per_kp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub_kp
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_kp

      CALL timeset(routineN, handle)

      IF (nkp > para_env%num_pe .OR. do_kpoints_cubic_RPA) THEN
         ! we have all kpoints on every processpr
         num_proc_per_kp = para_env%num_pe
      ELSE
         ! we have only one kpoint per group
         num_proc_per_kp = para_env%num_pe/nkp
      END IF

      color_sub_kp = MOD(para_env%mepos/num_proc_per_kp, nkp)
      CALL mp_comm_split_direct(para_env%group, comm_sub_kp, color_sub_kp)
      NULLIFY (para_env_sub_kp)
      CALL cp_para_env_create(para_env_sub_kp, comm_sub_kp)
      NULLIFY (blacs_env_sub_kp)
      CALL cp_blacs_env_create(blacs_env=blacs_env_sub_kp, para_env=para_env_sub_kp)

      NULLIFY (fm_struct_sub_kp)
      CALL cp_fm_struct_create(fm_struct_sub_kp, context=blacs_env_sub_kp, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env_sub_kp)

      CALL cp_para_env_release(para_env_sub_kp)

      CALL cp_blacs_env_release(blacs_env_sub_kp)

      ALLOCATE (ikp_local(nkp))
      ikp_local = 0
      first_ikp_local = 1
      DO ikp = 1, nkp
         IF (nkp > para_env%num_pe .OR. do_kpoints_cubic_RPA .OR. ikp == color_sub_kp + 1) THEN
            ikp_local(ikp) = ikp
            first_ikp_local = ikp
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_sub_para_kp

! **************************************************************************************************
!> \brief ...
!> \param do_mao ...
!> \param do_dbcsr_t ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param ikp_local ...
!> \param row_from_LLL ...
!> \param index_to_cell_3c ...
!> \param cell_to_index_3c ...
!> \param non_zero_blocks_3c ...
!> \param non_zero_blocks_3c_cut_col ...
!> \param do_ic_model ...
!> \param do_kpoints_cubic_RPA ...
!> \param do_kpoints_from_Gamma ...
!> \param do_ri_Sigma_x ...
!> \param cycle_due_to_sparse_dm ...
!> \param multiply_needed_occ ...
!> \param multiply_needed_virt ...
!> \param needed_cutRI_mem_R1vec_R2vec_for_kp ...
!> \param has_mat_P_blocks ...
!> \param buffer_mat_M ...
!> \param wkp_W ...
!> \param cfm_mat_Q ...
!> \param fm_mat_L ...
!> \param fm_mat_RI_global_work ...
!> \param fm_mat_work ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param mat_dm ...
!> \param mat_L ...
!> \param mat_M_P_munu_occ ...
!> \param mat_M_P_munu_virt ...
!> \param mat_P_global_copy ...
!> \param mat_SinvVSinv ...
!> \param mat_M_mu_Pnu_occ ...
!> \param mat_M_mu_Pnu_virt ...
!> \param mat_P_omega ...
!> \param mat_P_omega_kp ...
!> \param mat_dm_loc_occ_cut ...
!> \param mat_dm_loc_virt_cut ...
!> \param mat_3c_overl_int_cut ...
!> \param mat_3c_overl_int_mao_for_occ_cut ...
!> \param mat_3c_overl_int_mao_for_virt_cut ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param mat_dm_loc_occ ...
!> \param mat_dm_loc_virt ...
!> \param mat_work ...
!> \param fm_mo_coeff_occ_beta ...
!> \param fm_mo_coeff_virt_beta ...
!> \param mat_P_omega_beta ...
! **************************************************************************************************
   SUBROUTINE dealloc_im_time(do_mao, do_dbcsr_t, do_ri_sos_laplace_mp2, fm_mo_coeff_occ, fm_mo_coeff_virt, fm_scaled_dm_occ_tau, &
                              fm_scaled_dm_virt_tau, ikp_local, row_from_LLL, index_to_cell_3c, &
                              cell_to_index_3c, non_zero_blocks_3c, non_zero_blocks_3c_cut_col, do_ic_model, &
                              do_kpoints_cubic_RPA, do_kpoints_from_Gamma, do_ri_Sigma_x, cycle_due_to_sparse_dm, &
                              multiply_needed_occ, multiply_needed_virt, needed_cutRI_mem_R1vec_R2vec_for_kp, has_mat_P_blocks, &
                              buffer_mat_M, wkp_W, cfm_mat_Q, fm_mat_L, fm_mat_RI_global_work, fm_mat_work, &
                              fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, mat_dm, mat_L, mat_M_P_munu_occ, mat_M_P_munu_virt, &
                              mat_P_global_copy, mat_SinvVSinv, mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, mat_P_omega, mat_P_omega_kp, &
                              mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, mat_3c_overl_int_cut, mat_3c_overl_int_mao_for_occ_cut, &
                              mat_3c_overl_int_mao_for_virt_cut, t_3c_M, t_3c_O, mat_dm_loc_occ, mat_dm_loc_virt, mat_work, &
                              fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, mat_P_omega_beta)

      LOGICAL, INTENT(IN)                                :: do_mao, do_dbcsr_t, do_ri_sos_laplace_mp2
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ikp_local
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: row_from_LLL
      INTEGER, ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: cell_to_index_3c, non_zero_blocks_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :), &
         INTENT(INOUT)                                   :: non_zero_blocks_3c_cut_col
      LOGICAL, INTENT(IN)                                :: do_ic_model, do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma, do_ri_Sigma_x
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: cycle_due_to_sparse_dm, &
                                                            multiply_needed_occ, &
                                                            multiply_needed_virt
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :), INTENT(INOUT) :: &
         needed_cutRI_mem_R1vec_R2vec_for_kp
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :, :), &
         INTENT(INOUT)                                   :: has_mat_P_blocks
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: buffer_mat_M
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp_W
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_type), POINTER                          :: fm_mat_RI_global_work, fm_mat_work, &
                                                            fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_dm, mat_L, mat_M_P_munu_occ, &
                                                            mat_M_P_munu_virt, mat_P_global_copy, &
                                                            mat_SinvVSinv
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega, mat_P_omega_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_dm_loc_occ_cut, mat_dm_loc_virt_cut
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut, &
         mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut
      TYPE(dbcsr_t_type)                                 :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_O
      TYPE(dbcsr_type), POINTER                          :: mat_dm_loc_occ, mat_dm_loc_virt, mat_work
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mo_coeff_occ_beta, &
                                                            fm_mo_coeff_virt_beta
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: mat_P_omega_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dealloc_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_size, j_size
      LOGICAL                                            :: my_open_shell

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(fm_mo_coeff_occ_beta) .AND. PRESENT(fm_mo_coeff_virt_beta) .AND. PRESENT(mat_P_omega_beta)) THEN
         my_open_shell = .TRUE.
      END IF

      CALL cp_fm_release(fm_scaled_dm_occ_tau)
      CALL cp_fm_release(fm_scaled_dm_virt_tau)
      CALL cp_fm_release(fm_mo_coeff_occ)
      CALL cp_fm_release(fm_mo_coeff_virt)
      CALL cp_fm_release(fm_mo_coeff_occ_scaled)
      CALL cp_fm_release(fm_mo_coeff_virt_scaled)

      IF (my_open_shell) THEN
         CALL cp_fm_release(fm_mo_coeff_occ_beta)
         CALL cp_fm_release(fm_mo_coeff_virt_beta)
      END IF

      DO i_size = 1, SIZE(fm_mat_L, 1)
         DO j_size = 1, SIZE(fm_mat_L, 2)
            IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
               IF (ANY(ikp_local(:) == i_size)) THEN
                  CALL cp_fm_release(fm_mat_L(i_size, j_size)%matrix)
               END IF
            ELSE
               CALL cp_fm_release(fm_mat_L(i_size, j_size)%matrix)
            END IF
         END DO
      END DO
      DEALLOCATE (fm_mat_L)
      CALL cp_fm_release(fm_mat_work)

      IF (.NOT. do_dbcsr_t) THEN
         CALL dbcsr_release_p(mat_dm_loc_occ)
         CALL dbcsr_release_p(mat_dm_loc_virt)

         CALL dbcsr_deallocate_matrix_set(mat_dm_loc_occ_cut)
         CALL dbcsr_deallocate_matrix_set(mat_dm_loc_virt_cut)

         CALL dbcsr_release(mat_M_P_munu_occ%matrix)
         DEALLOCATE (mat_M_P_munu_occ%matrix)

         CALL dbcsr_release(mat_M_P_munu_virt%matrix)
         DEALLOCATE (mat_M_P_munu_virt%matrix)

         CALL dbcsr_release(mat_P_global_copy%matrix)
         DEALLOCATE (mat_P_global_copy%matrix)
      ENDIF

      IF (.NOT. (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma)) THEN
         CALL dbcsr_release(mat_work)
         DEALLOCATE (mat_work)
      END IF

      CALL dbcsr_release(mat_L%matrix)
      DEALLOCATE (mat_L%matrix)
      IF (do_ri_Sigma_x .OR. do_ic_model) THEN
         CALL dbcsr_release(mat_SinvVSinv%matrix)
         DEALLOCATE (mat_SinvVSinv%matrix)
      END IF
      IF (do_ri_Sigma_x) THEN
         CALL dbcsr_release(mat_dm%matrix)
         DEALLOCATE (mat_dm%matrix)
      END IF

      DEALLOCATE (index_to_cell_3c, cell_to_index_3c)

      IF (.NOT. do_dbcsr_t) THEN
         CALL dbcsr_deallocate_matrix_set(mat_M_mu_Pnu_occ)
         CALL dbcsr_deallocate_matrix_set(mat_M_mu_Pnu_virt)
      ENDIF

      CALL dbcsr_deallocate_matrix_set(mat_P_omega)
      IF (my_open_shell .AND. do_ri_sos_laplace_mp2) CALL dbcsr_deallocate_matrix_set(mat_P_omega_beta)

      IF (.NOT. do_dbcsr_t) THEN
         CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_cut)

         IF (do_mao) THEN
            CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_mao_for_occ_cut)
            CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_mao_for_virt_cut)
         END IF
      ELSE
         DO i_size = 1, SIZE(t_3c_O, 1)
            DO j_size = 1, SIZE(t_3c_O, 2)
               CALL dbcsr_t_destroy(t_3c_O(i_size, j_size))
            ENDDO
         ENDDO

         DEALLOCATE (t_3c_O)
         CALL dbcsr_t_destroy(t_3c_M)
      ENDIF

      IF (.NOT. do_dbcsr_t) THEN
         DEALLOCATE (buffer_mat_M)
         DEALLOCATE (non_zero_blocks_3c)
         DEALLOCATE (non_zero_blocks_3c_cut_col)
         DEALLOCATE (cycle_due_to_sparse_dm, multiply_needed_occ, multiply_needed_virt)
         DEALLOCATE (row_from_LLL)
         DEALLOCATE (needed_cutRI_mem_R1vec_R2vec_for_kp)
      ENDIF
      DEALLOCATE (has_mat_P_blocks)

      IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
         CALL cp_cfm_release(cfm_mat_Q)
         CALL cp_fm_release(fm_mat_RI_global_work)
         CALL dbcsr_deallocate_matrix_set(mat_P_omega_kp)
         DEALLOCATE (wkp_W)
      END IF

      CALL timestop(handle)

   END SUBROUTINE dealloc_im_time

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param mat_L ...
!> \param mat_work ...
!> \param eps_filter_im_time ...
!> \param fm_mat_work ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param fm_mat_L ...
!> \param fm_mat_Q ...
! **************************************************************************************************
   SUBROUTINE contract_P_omega_with_mat_L(mat_P_omega, mat_L, mat_work, eps_filter_im_time, fm_mat_work, dimen_RI, &
                                          dimen_RI_red, fm_mat_L, fm_mat_Q)

      TYPE(dbcsr_type), POINTER                          :: mat_P_omega, mat_L, mat_work
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time
      TYPE(cp_fm_type), POINTER                          :: fm_mat_work
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L, fm_mat_Q

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_P_omega_with_mat_L', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! multiplication with RI metric/Coulomb operator
      CALL dbcsr_multiply("N", "T", 1.0_dp, mat_P_omega, mat_L, &
                          0.0_dp, mat_work, filter_eps=eps_filter_im_time)

      CALL copy_dbcsr_to_fm(mat_work, fm_mat_work)

      CALL cp_gemm('N', 'N', dimen_RI_red, dimen_RI_red, dimen_RI, 1.0_dp, fm_mat_L, fm_mat_work, &
                   0.0_dp, fm_mat_Q)

      ! Reset mat_work to save memory
      CALL dbcsr_set(mat_work, 0.0_dp)
      CALL dbcsr_filter(mat_work, 1.0_dp)

      CALL timestop(handle)

   END SUBROUTINE contract_P_omega_with_mat_L

END MODULE rpa_util
