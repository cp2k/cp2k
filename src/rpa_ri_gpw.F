!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate RI-RPA energy
!> \par History
!>      06.2012 created [Mauro Del Ben]
!>      04.2015 GW routines added [Jan Wilhelm]
!>      10.2015 Cubic-scaling RPA routines added [Jan Wilhelm]
!>      10.2018 Cubic-scaling SOS-MP2 added [Frederick Stein]
!>      03.2019 Refactoring [Frederick Stein]
! **************************************************************************************************
MODULE rpa_ri_gpw
   USE bibliography,                    ONLY: Bates2013,&
                                              DelBen2013,&
                                              DelBen2015,&
                                              cite_reference
   USE bse,                             ONLY: do_subspace_iterations,&
                                              mult_B_with_W_and_fill_local_3c_arrays
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_row_template,&
                                              cp_dbcsr_m_by_n_from_template,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_syrk,&
                                              cp_fm_transpose
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_indxg2l, cp_fm_indxg2p, &
        cp_fm_p_type, cp_fm_release, cp_fm_set_all, cp_fm_set_element, cp_fm_to_fm, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_filter, dbcsr_get_info, &
        dbcsr_get_occupation, dbcsr_get_stored_coordinates, dbcsr_init_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_release_p, dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_trace, dbcsr_transposed, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_default
   USE input_constants,                 ONLY: wfc_mm_style_gemm,&
                                              wfc_mm_style_syrk
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_memory,&
                                              m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: &
        mp_alltoall, mp_bcast, mp_comm_split_direct, mp_irecv, mp_isend, mp_min, mp_sendrecv, &
        mp_sum, mp_wait, mp_waitall
   USE minimax_exp,                     ONLY: check_exp_minimax_range,&
                                              get_exp_minimax_coeff
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff
   USE mp2_laplace,                     ONLY: laplace_minimax_approx
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type,&
                                              two_dim_int_array
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_axk,                         ONLY: compute_axk_ener
   USE rpa_communication,               ONLY: fm_redistribute,&
                                              initialize_buffer,&
                                              release_buffer
   USE rpa_gw,                          ONLY: GW_matrix_operations,&
                                              allocate_matrices_gw,&
                                              allocate_matrices_gw_im_time,&
                                              get_weights_gw,&
                                              gw_cleanup,&
                                              gw_postprocessing
   USE rpa_im_time,                     ONLY: communicate_buffer,&
                                              compute_mat_P_omega,&
                                              gap_and_max_eig_diff_kpoints,&
                                              setup_mat_for_mem_cut_3c,&
                                              zero_mat_P_omega
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_ri_gpw'

   PUBLIC :: rpa_ri_compute_en, test_least_square_ft

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param Erpa ...
!> \param mp2_env ...
!> \param BIb_C ...
!> \param BIb_C_gw ...
!> \param BIb_C_bse_ij ...
!> \param BIb_C_bse_ab ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param ends_array ...
!> \param ends_B_virtual ...
!> \param ends_B_all ...
!> \param sizes_array ...
!> \param sizes_B_virtual ...
!> \param sizes_B_all ...
!> \param starts_array ...
!> \param starts_B_virtual ...
!> \param starts_B_all ...
!> \param starts_B_occ_bse ...
!> \param sizes_B_occ_bse ...
!> \param ends_B_occ_bse ...
!> \param starts_B_virt_bse ...
!> \param sizes_B_virt_bse ...
!> \param ends_B_virt_bse ...
!> \param mo_coeff ...
!> \param fm_matrix_L_RI_metric ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param unit_nr ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param my_do_gw ...
!> \param do_im_time ...
!> \param do_mao ...
!> \param do_bse ...
!> \param matrix_s ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param mao_coeff_occ_A ...
!> \param mao_coeff_virt_A ...
!> \param mat_munu ...
!> \param mat_dm_occ_local ...
!> \param mat_dm_virt_local ...
!> \param mat_P_local ...
!> \param mat_P_global ...
!> \param mat_M ...
!> \param mat_3c_overl_int ...
!> \param mat_3c_overl_int_mao_for_occ ...
!> \param mat_3c_overl_int_mao_for_virt ...
!> \param eps_filter ...
!> \param BIb_C_beta ...
!> \param homo_beta ...
!> \param Eigenval_beta ...
!> \param ends_B_virtual_beta ...
!> \param sizes_B_virtual_beta ...
!> \param starts_B_virtual_beta ...
!> \param mo_coeff_beta ...
!> \param BIb_C_gw_beta ...
!> \param gw_corr_lev_occ_beta ...
!> \param gw_corr_lev_virt_beta ...
! **************************************************************************************************
   SUBROUTINE rpa_ri_compute_en(qs_env, Erpa, mp2_env, BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, &
                                para_env, para_env_sub, color_sub, &
                                ends_array, ends_B_virtual, ends_B_all, sizes_array, sizes_B_virtual, sizes_B_all, &
                                starts_array, starts_B_virtual, starts_B_all, starts_B_occ_bse, sizes_B_occ_bse, &
                                ends_B_occ_bse, starts_B_virt_bse, sizes_B_virt_bse, ends_B_virt_bse, &
                                mo_coeff, fm_matrix_L_RI_metric, &
                                Eigenval, nmo, homo, dimen_RI, gw_corr_lev_occ, gw_corr_lev_virt, &
                                unit_nr, do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_mao, do_bse, matrix_s, &
                                mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                mat_munu, mat_dm_occ_local, mat_dm_virt_local, &
                                mat_P_local, mat_P_global, &
                                mat_M, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                                mat_3c_overl_int_mao_for_virt, &
                                eps_filter, BIb_C_beta, homo_beta, Eigenval_beta, &
                                ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, &
                                mo_coeff_beta, BIb_C_gw_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: Erpa
      TYPE(mp2_type), POINTER                            :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: BIb_C, BIb_C_gw, BIb_C_bse_ij, &
                                                            BIb_C_bse_ab
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(INOUT)                             :: color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: ends_array, ends_B_virtual, ends_B_all, &
         sizes_array, sizes_B_virtual, sizes_B_all, starts_array, starts_B_virtual, starts_B_all, &
         starts_B_occ_bse, sizes_B_occ_bse, ends_B_occ_bse, starts_B_virt_bse, sizes_B_virt_bse, &
         ends_B_virt_bse
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      INTEGER, INTENT(IN)                                :: nmo, homo, dimen_RI, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, unit_nr
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2, my_do_gw, &
                                                            do_im_time, do_mao, do_bse
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, mao_coeff_occ, mao_coeff_virt, &
                                                            mao_coeff_occ_A, mao_coeff_virt_A
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu, mat_dm_occ_local, &
                                                            mat_dm_virt_local, mat_P_local, &
                                                            mat_P_global, mat_M
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT), OPTIONAL                         :: BIb_C_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: Eigenval_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: ends_B_virtual_beta, &
                                                            sizes_B_virtual_beta, &
                                                            starts_B_virtual_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mo_coeff_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT), OPTIONAL                         :: BIb_C_gw_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: gw_corr_lev_occ_beta, &
                                                            gw_corr_lev_virt_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_ri_compute_en', &
         routineP = moduleN//':'//routineN

      INTEGER :: best_integ_group_size, best_num_integ_point, color_rpa_group, dimen_homo_square, &
         dimen_ia, dimen_ia_beta, dimen_nm_gw, dimen_virt_square, handle, handle2, handle3, i, &
         ierr, iiB, input_integ_group_size, integ_group_size, jjB, min_integ_group_size, &
         my_ab_comb_bse_end, my_ab_comb_bse_size, my_ab_comb_bse_start, my_group_L_end, &
         my_group_L_size, my_group_L_start, my_homo_beta, my_ia_end, my_ia_end_beta, my_ia_size, &
         my_ia_size_beta, my_ia_start, my_ia_start_beta, my_ij_comb_bse_end, my_ij_comb_bse_size, &
         my_ij_comb_bse_start, my_nm_gw_end, my_nm_gw_size, my_nm_gw_start, ncol_block_mat
      INTEGER :: ngroup, nrow_block_mat, num_integ_group, num_integ_points, pos_integ_group, &
         virtual, virtual_beta
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sub_proc_map
      LOGICAL                                            :: do_minimax_quad, my_open_shell, &
                                                            skip_integ_group_opt
      REAL(KIND=dp) :: allowed_memory, avail_mem, E_Range, Emax, Emax_beta, Emin, Emin_beta, &
         mem_for_iaK, mem_for_QK, mem_min, mem_per_group, mem_real, needed_mem
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D, BIb_C_2D_beta, &
                                                            BIb_C_2D_bse_ab, BIb_C_2D_bse_ij, &
                                                            BIb_C_2D_gw, BIb_C_2D_gw_beta
      TYPE(cp_fm_type), POINTER :: fm_mat_dummy_1, fm_mat_dummy_2, fm_mat_dummy_3, fm_mat_dummy_4, &
         fm_mat_Q, fm_mat_Q_beta, fm_mat_Q_gemm, fm_mat_Q_gemm_beta, fm_mat_R_gw, &
         fm_mat_R_gw_beta, fm_mat_R_gw_gemm, fm_mat_R_gw_gemm_beta, fm_mat_S, fm_mat_S_ab_bse, &
         fm_mat_S_beta, fm_mat_S_gw, fm_mat_S_gw_beta, fm_mat_S_ij_bse, fm_mo_coeff_occ, &
         fm_mo_coeff_occ_beta, fm_mo_coeff_virt, fm_mo_coeff_virt_beta, fm_scaled_dm_occ_tau, &
         fm_scaled_dm_virt_tau
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA

      CALL timeset(routineN, handle)

      CALL cite_reference(DelBen2013)
      CALL cite_reference(DelBen2015)

      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         CALL cite_reference(Bates2013)

      ENDIF

      my_open_shell = .FALSE.
      IF (PRESENT(BIb_C_beta) .AND. &
          PRESENT(ends_B_virtual_beta) .AND. &
          PRESENT(sizes_B_virtual_beta) .AND. &
          PRESENT(starts_B_virtual_beta) .AND. &
          PRESENT(homo_beta) .AND. &
          PRESENT(Eigenval_beta)) my_open_shell = .TRUE.

      virtual = nmo-homo
      IF (my_open_shell) THEN
         virtual_beta = nmo-homo_beta
      END IF

      IF (do_ri_sos_laplace_mp2) THEN
         num_integ_points = mp2_env%ri_laplace%n_quadrature
         input_integ_group_size = mp2_env%ri_laplace%integ_group_size

         ! check the range for the minimax approximation
         Emin = 2.0_dp*(Eigenval(homo+1)-Eigenval(homo))
         Emax = 2.0_dp*(MAXVAL(Eigenval)-MINVAL(Eigenval))
         IF (my_open_shell) THEN
            IF (homo_beta > 0) THEN
               Emin_beta = 2.0_dp*(Eigenval_beta(homo_beta+1)-Eigenval_beta(homo_beta))
               Emax_beta = 2.0_dp*(MAXVAL(Eigenval_beta)-MINVAL(Eigenval_beta))
               Emin = MIN(Emin, Emin_beta)
               Emax = MAX(Emax, Emax_beta)
            END IF
         END IF
         E_Range = Emax/Emin
         IF (E_Range < 2.0_dp) E_Range = 2.0_dp
         ierr = 0
         CALL check_exp_minimax_range(num_integ_points, E_Range, ierr)
         IF (ierr /= 0) THEN
            jjB = num_integ_points-1
            DO iiB = 1, jjB
               num_integ_points = num_integ_points-1
               ierr = 0
               CALL check_exp_minimax_range(num_integ_points, E_Range, ierr)
               IF (ierr == 0) EXIT
            END DO
         END IF
         CPASSERT(num_integ_points >= 1)
      ELSE
         num_integ_points = mp2_env%ri_rpa%rpa_num_quad_points
         input_integ_group_size = mp2_env%ri_rpa%rpa_integ_group_size
         do_minimax_quad = mp2_env%ri_rpa%minimax_quad
         IF (do_minimax_quad .AND. num_integ_points > 20) THEN
            CALL cp_warn(__LOCATION__, &
                         "The required number of quadrature point exceeds the maximum possible in the "// &
                         "Minimax quadrature scheme. The number of quadrature point has been reset to 20.")
            num_integ_points = 20
         END IF
      END IF
      allowed_memory = mp2_env%mp2_memory

      my_group_L_size = sizes_array(color_sub)
      my_group_L_start = starts_array(color_sub)
      my_group_L_end = ends_array(color_sub)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! for imaginary time or periodic GW or BSE, we use all processors for a single frequency/time point
      IF (do_im_time .OR. mp2_env%ri_g0w0%do_periodic .OR. do_bse) THEN

         IF (do_im_time) color_sub = para_env%mepos/mp2_env%ri_rpa_im_time%group_size_3c

         integ_group_size = ngroup
         best_num_integ_point = num_integ_points

      ELSE

         ! Calculate available memory and create integral group according to that
         ! mem_for_iaK is the memory needed for storing the 3 centre integrals
         mem_for_iaK = REAL(homo, KIND=dp)*virtual*dimen_RI*8.0_dp/(1024_dp**2)
         mem_for_QK = REAL(dimen_RI, KIND=dp)*dimen_RI*8.0_dp/(1024_dp**2)

         IF (my_open_shell) THEN
            mem_for_iaK = mem_for_iaK+REAL(homo_beta, KIND=dp)*virtual_beta*dimen_RI*8.0_dp/(1024_dp**2)
            mem_for_QK = mem_for_QK*2.0_dp
         END IF

         CALL m_memory(mem)
         mem_real = (mem+1024*1024-1)/(1024*1024)
         ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
         ! has not been given back to the OS yet.
         CALL mp_min(mem_real, para_env%group)

         mem_min = 2.0_dp*REAL(homo, KIND=dp)*MAXVAL(sizes_B_virtual)*MAXVAL(sizes_array)*8.0_dp/(1024**2)

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Minimum required memory per MPI process:', &
            mem_min, ' MiB'

         mem_real = allowed_memory-mem_real
         mem_real = MAX(mem_real, mem_min)

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Available memory per MPI process:', &
            mem_real, ' MiB'

         mem_per_group = mem_real*para_env_sub%num_pe

         needed_mem = mem_for_iaK*2.0_dp+mem_for_QK*3.0_dp

         ! here we try to find the best rpa/laplace group size
         skip_integ_group_opt = .FALSE.

         IF (input_integ_group_size > 0) THEN
            IF (MOD(input_integ_group_size, para_env_sub%num_pe) == 0) THEN
               best_integ_group_size = input_integ_group_size/para_env_sub%num_pe
               IF (MOD(ngroup, best_integ_group_size) == 0) THEN
                  num_integ_group = ngroup/best_integ_group_size
                  IF ((num_integ_points > num_integ_group) .AND. MOD(num_integ_points, num_integ_group) == 0) THEN
                     best_num_integ_point = num_integ_points/num_integ_group
                     skip_integ_group_opt = .TRUE.
                  ELSE
                     IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'NUM_QUAD_POINTS not multiple of the number of INTEG_GROUP'
                  END IF
               ELSE
                  IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Total number of groups not multiple of SIZE_INTEG_GROUP'
               END IF
            ELSE
               IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'SIZE_INTEG_GROUP not multiple of GROUP_SIZE'
            END IF
         END IF

         IF (.NOT. skip_integ_group_opt) THEN
            best_integ_group_size = ngroup
            best_num_integ_point = num_integ_points

            min_integ_group_size = MAX(1, ngroup/num_integ_points)

            integ_group_size = min_integ_group_size-1
            DO iiB = min_integ_group_size+1, ngroup
               integ_group_size = integ_group_size+1

               ! check that the ngroup is a multiple of  integ_group_size
               IF (MOD(ngroup, integ_group_size) /= 0) CYCLE

               ! check for memory
               avail_mem = integ_group_size*mem_per_group
               IF (avail_mem < needed_mem) CYCLE

               ! check the number of integration points is a multiple of the  number of integ_group
               num_integ_group = ngroup/integ_group_size
               IF (num_integ_points < num_integ_group) CYCLE
               IF (MOD(num_integ_points, num_integ_group) /= 0) CYCLE

               ! if all the test passed then decide
               IF ((num_integ_points/num_integ_group) < best_num_integ_point) THEN
                  best_num_integ_point = num_integ_points/num_integ_group
                  best_integ_group_size = integ_group_size
               END IF

            END DO
         END IF

         integ_group_size = best_integ_group_size

      END IF

      IF (unit_nr > 0) THEN
         IF (do_ri_sos_laplace_mp2) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Group size for laplace numerical integration:", integ_group_size*para_env_sub%num_pe
            WRITE (UNIT=unit_nr, FMT="(T3,A)") &
               "INTEG_INFO| MINIMAX approximation"
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "INTEG_INFO| Number of integration points:", num_integ_points
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "INTEG_INFO| Number of integration points per Laplace group:", best_num_integ_point
         ELSE
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Group size for frequency integration:", integ_group_size*para_env_sub%num_pe
            IF (do_minimax_quad) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A)") &
                  "INTEG_INFO| MINIMAX quadrature"
            ELSE
               WRITE (UNIT=unit_nr, FMT="(T3,A)") &
                  "INTEG_INFO| Clenshaw-Curtius quadrature"
            END IF
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "INTEG_INFO| Number of integration points:", num_integ_points
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "INTEG_INFO| Number of integration points per RPA group:", best_num_integ_point
         END IF
         CALL m_flush(unit_nr)
      END IF

      num_integ_group = ngroup/integ_group_size

      pos_integ_group = MOD(color_sub, integ_group_size)
      color_rpa_group = color_sub/integ_group_size

      IF (.NOT. do_im_time) THEN
         IF (my_open_shell) THEN
            my_homo_beta = homo_beta
         ELSE
            my_homo_beta = homo
         END IF

      END IF ! not imaginary time

      CALL timeset(routineN//"_reorder", handle2)

      ! create the sub_proc_map
      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
      DO i = 0, para_env_sub%num_pe-1
         sub_proc_map(i) = i
         sub_proc_map(-i-1) = para_env_sub%num_pe-i-1
         sub_proc_map(para_env_sub%num_pe+i) = i
      END DO

      ! not necessary for imaginary time

      IF (do_im_time) THEN

         dimen_ia = homo*virtual

      ELSE

         ! reorder the local data in such a way to help the next stage of matrix creation
         ! now the data inside the group are divided into a ia x K matrix
         CALL calculate_BIb_C_2D(BIb_C_2D, BIb_C, para_env_sub, dimen_ia, homo, virtual, &
                                 sizes_B_virtual, starts_B_virtual, ends_B_virtual, &
                                 sub_proc_map, my_ia_size, my_ia_start, my_ia_end, my_group_L_size)
         DEALLOCATE (BIb_C)
         DEALLOCATE (starts_B_virtual)
         DEALLOCATE (ends_B_virtual)
         DEALLOCATE (sizes_B_virtual)

         ! The same for open shell
         IF (my_open_shell) THEN
            CALL calculate_BIb_C_2D(BIb_C_2D_beta, BIb_C_beta, para_env_sub, &
                                    dimen_ia_beta, homo_beta, virtual_beta, &
                                    sizes_B_virtual_beta, starts_B_virtual_beta, ends_B_virtual_beta, &
                                    sub_proc_map, my_ia_size_beta, my_ia_start_beta, my_ia_end_beta, my_group_L_size)

            DEALLOCATE (BIb_C_beta)
            DEALLOCATE (starts_B_virtual_beta)
            DEALLOCATE (ends_B_virtual_beta)
            DEALLOCATE (sizes_B_virtual_beta)

         END IF

      END IF

      ! in the GW case, BIb_C_2D_gw is an nm x K matrix, with n: number of corr GW levels, m=nmo
      IF (my_do_gw) THEN

         CALL timeset(routineN//"_reorder_gw", handle3)

         CALL calculate_BIb_C_2D(BIb_C_2D_gw, BIb_C_gw, para_env_sub, dimen_nm_gw, &
                                 gw_corr_lev_occ+gw_corr_lev_virt, nmo, &
                                 sizes_B_all, starts_B_all, ends_B_all, &
                                 sub_proc_map, my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, my_group_L_size)

         ! The same for open shell
         IF (my_open_shell) THEN
            CALL calculate_BIb_C_2D(BIb_C_2D_gw_beta, BIb_C_gw_beta, para_env_sub, dimen_nm_gw, &
                                    gw_corr_lev_occ+gw_corr_lev_virt, nmo, &
                                    sizes_B_all, starts_B_all, ends_B_all, &
                                    sub_proc_map, my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, my_group_L_size)
            DEALLOCATE (BIb_C_gw_beta)
         END IF

         DEALLOCATE (BIb_C_gw)
         DEALLOCATE (starts_B_all)
         DEALLOCATE (ends_B_all)
         DEALLOCATE (sizes_B_all)

         CALL timestop(handle3)

      END IF

      IF (do_bse) THEN

         CALL timeset(routineN//"_reorder_bse1", handle3)

         CALL calculate_BIb_C_2D(BIb_C_2D_bse_ij, BIb_C_bse_ij, para_env_sub, dimen_homo_square, &
                                 homo, homo, &
                                 sizes_B_occ_bse, starts_B_occ_bse, ends_B_occ_bse, &
                                 sub_proc_map, my_ij_comb_bse_size, my_ij_comb_bse_start, my_ij_comb_bse_end, my_group_L_size)

         DEALLOCATE (BIb_C_bse_ij)
         DEALLOCATE (starts_B_occ_bse)
         DEALLOCATE (ends_B_occ_bse)
         DEALLOCATE (sizes_B_occ_bse)

         CALL timestop(handle3)

         CALL timeset(routineN//"_reorder_bse2", handle3)

         CALL calculate_BIb_C_2D(BIb_C_2D_bse_ab, BIb_C_bse_ab, para_env_sub, dimen_virt_square, &
                                 virtual, virtual, &
                                 sizes_B_virt_bse, starts_B_virt_bse, ends_B_virt_bse, &
                                 sub_proc_map, my_ab_comb_bse_size, my_ab_comb_bse_start, my_ab_comb_bse_end, my_group_L_size)

         DEALLOCATE (BIb_C_bse_ab)
         DEALLOCATE (starts_B_virt_bse)
         DEALLOCATE (ends_B_virt_bse)
         DEALLOCATE (sizes_B_virt_bse)

         CALL timestop(handle3)

      END IF

      CALL timestop(handle2)

      ! now create the matrices needed for the calculation, Q, S and G
      ! Q and G will have omega dependence

      IF (my_open_shell .AND. do_im_time .AND. do_ri_sos_laplace_mp2) THEN
         CALL create_integ_mat(BIb_C_2D, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_ia, dimen_ia, color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ia_size, my_ia_start, my_ia_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                               nrow_block_mat, ncol_block_mat, &
                               do_im_time=do_im_time, fm_scaled_dm_occ_tau=fm_scaled_dm_occ_tau, &
                               fm_scaled_dm_virt_tau=fm_scaled_dm_virt_tau, mo_coeff=mo_coeff, &
                               fm_mo_coeff_occ=fm_mo_coeff_occ, fm_mo_coeff_virt=fm_mo_coeff_virt, &
                               nmo=nmo, homo=homo, do_mao=do_mao, mao_coeff_occ_A=mao_coeff_occ_A, &
                               mao_coeff_virt_A=mao_coeff_virt_A, matrix_s=matrix_s, fm_mat_Q_beta=fm_mat_Q_beta)
      ELSE
         CALL create_integ_mat(BIb_C_2D, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_ia, dimen_ia, color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ia_size, my_ia_start, my_ia_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                               nrow_block_mat, ncol_block_mat, &
                               do_im_time=do_im_time, fm_scaled_dm_occ_tau=fm_scaled_dm_occ_tau, &
                               fm_scaled_dm_virt_tau=fm_scaled_dm_virt_tau, mo_coeff=mo_coeff, &
                               fm_mo_coeff_occ=fm_mo_coeff_occ, fm_mo_coeff_virt=fm_mo_coeff_virt, &
                               nmo=nmo, homo=homo, do_mao=do_mao, mao_coeff_occ_A=mao_coeff_occ_A, &
                               mao_coeff_virt_A=mao_coeff_virt_A, matrix_s=matrix_s)
      END IF

      IF (my_open_shell) THEN

         ! for imaginary time, we only have to build beta mo coefficients
         IF (do_im_time) THEN

            CALL create_occ_virt_mo_coeffs(fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, mo_coeff_beta, &
                                           nmo, homo_beta, do_mao, mao_coeff_occ_A, &
                                           mao_coeff_virt_A, matrix_s, 2)

            ! for RPA with imaginary frequency, we have to build the same matrices for beta as for alpha
         ELSE

            CALL create_integ_mat(BIb_C_2D_beta, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                                  dimen_RI, dimen_ia_beta, dimen_ia_beta, color_rpa_group, &
                                  mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                                  my_ia_size_beta, my_ia_start_beta, my_ia_end_beta, &
                                  my_group_L_size, my_group_L_start, my_group_L_end, &
                                  para_env_RPA, fm_mat_S_beta, fm_mat_Q_gemm_beta, fm_mat_Q_beta, &
                                  nrow_block_mat, ncol_block_mat, &
                                  .TRUE., fm_mat_Q%matrix_struct%context)

         END IF

      END IF

      ! for GW, we need other matrix fm_mat_S
      IF (my_do_gw) THEN

         CALL create_integ_mat(BIb_C_2D_gw, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_nm_gw, dimen_ia, color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_gw, fm_mat_R_gw_gemm, fm_mat_R_gw, &
                               nrow_block_mat, ncol_block_mat, &
                               .TRUE., fm_mat_Q%matrix_struct%context, fm_mat_Q%matrix_struct%context)

         ! for GW, we don't need fm_mat_R_gw_gemm (in contrast to RPA)
         CALL cp_fm_release(fm_mat_R_gw_gemm)

         IF (my_open_shell) THEN
            CALL create_integ_mat(BIb_C_2D_gw_beta, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                                  dimen_RI, dimen_nm_gw, dimen_ia, color_rpa_group, &
                                  mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                                  my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, &
                                  my_group_L_size, my_group_L_start, my_group_L_end, &
                                  para_env_RPA, fm_mat_S_gw_beta, fm_mat_R_gw_gemm_beta, fm_mat_R_gw_beta, &
                                  nrow_block_mat, ncol_block_mat, &
                                  .TRUE., fm_mat_Q%matrix_struct%context, fm_mat_Q%matrix_struct%context, &
                                  do_gw_open_shell=.TRUE.)

         END IF

      END IF

      ! for Bethe-Salpeter, we need other matrix fm_mat_S
      IF (do_bse) THEN

         CALL create_integ_mat(BIb_C_2D_bse_ij, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_homo_square, dimen_ia, color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ij_comb_bse_size, my_ij_comb_bse_start, my_ij_comb_bse_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_ij_bse, fm_mat_dummy_1, fm_mat_dummy_2, &
                               nrow_block_mat, ncol_block_mat, &
                               .TRUE., fm_mat_Q%matrix_struct%context, fm_mat_Q%matrix_struct%context)

         CALL cp_fm_release(fm_mat_dummy_1)
         CALL cp_fm_release(fm_mat_dummy_2)

         CALL create_integ_mat(BIb_C_2D_bse_ab, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_virt_square, dimen_ia, color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ab_comb_bse_size, my_ab_comb_bse_start, my_ab_comb_bse_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_ab_bse, fm_mat_dummy_3, fm_mat_dummy_4, &
                               nrow_block_mat, ncol_block_mat, &
                               .TRUE., fm_mat_Q%matrix_struct%context, fm_mat_Q%matrix_struct%context)

         CALL cp_fm_release(fm_mat_dummy_3)
         CALL cp_fm_release(fm_mat_dummy_4)

      END IF

      IF (my_open_shell) THEN

         IF (do_ri_sos_laplace_mp2 .AND. .NOT. do_im_time) THEN
            ! go with laplace MINIMAX MP2
            CALL laplace_minimax_approx(Erpa, para_env, para_env_RPA, unit_nr, homo, virtual, dimen_RI, dimen_ia, Eigenval, &
                                        num_integ_points, num_integ_group, color_rpa_group, &
                                        fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                                        homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta, &
                                        fm_mat_Q_gemm_beta, fm_mat_Q_beta)
         ELSE

            ! go with clenshaw-curtius/minimax quadrature
            CALL rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                             homo, virtual, dimen_RI, dimen_ia, dimen_nm_gw, &
                             Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                             fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                             fm_mat_S_ij_bse, fm_mat_S_ij_bse, &
                             my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                             do_minimax_quad, &
                             do_im_time, do_mao, fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                             fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                             mat_munu, mat_dm_occ_local, mat_dm_virt_local, mat_P_local, mat_P_global, mat_M, &
                             mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, matrix_s, &
                             mao_coeff_occ, mao_coeff_virt, eps_filter, &
                             starts_array, ends_array, sizes_array, color_sub, &
                             fm_mo_coeff_occ_beta=fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta=fm_mo_coeff_virt_beta, &
                             homo_beta=homo_beta, virtual_beta=virtual_beta, &
                             dimen_ia_beta=dimen_ia_beta, Eigenval_beta=Eigenval_beta, fm_mat_S_beta=fm_mat_S_beta, &
                             fm_mat_Q_gemm_beta=fm_mat_Q_gemm_beta, fm_mat_Q_beta=fm_mat_Q_beta, &
                             fm_mat_S_gw_beta=fm_mat_S_gw_beta, gw_corr_lev_occ_beta=gw_corr_lev_occ_beta, &
                             gw_corr_lev_virt_beta=gw_corr_lev_virt_beta, mo_coeff_beta=mo_coeff_beta, &
                             do_ri_sos_laplace_mp2=do_ri_sos_laplace_mp2)
         END IF
      ELSE
         IF (do_ri_sos_laplace_mp2 .AND. .NOT. do_im_time) THEN
            ! go with laplace MINIMAX MP2
            CALL laplace_minimax_approx(Erpa, para_env, para_env_RPA, unit_nr, homo, virtual, dimen_RI, dimen_ia, Eigenval, &
                                        num_integ_points, num_integ_group, color_rpa_group, &
                                        fm_mat_S, fm_mat_Q_gemm, fm_mat_Q)
         ELSE
            ! go with clenshaw-curtius/minimax quadrature
            ! here, we also do the quasi-particle-energy correction for G0W0
            CALL rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                             homo, virtual, dimen_RI, dimen_ia, dimen_nm_gw, &
                             Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                             fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                             fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                             my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                             do_minimax_quad, &
                             do_im_time, do_mao, fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                             fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                             mat_munu, mat_dm_occ_local, mat_dm_virt_local, mat_P_local, mat_P_global, mat_M, &
                             mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, matrix_s, &
                             mao_coeff_occ, mao_coeff_virt, &
                             eps_filter, starts_array, ends_array, sizes_array, color_sub, &
                             do_ri_sos_laplace_mp2=do_ri_sos_laplace_mp2)
         END IF
      END IF

      DEALLOCATE (sub_proc_map)

      DEALLOCATE (sizes_array)
      DEALLOCATE (starts_array)
      DEALLOCATE (ends_array)

      CALL cp_para_env_release(para_env_RPA)

      IF (.NOT. do_im_time) THEN
         CALL cp_fm_release(fm_mat_S)
         CALL cp_fm_release(fm_mat_Q_gemm)
      END IF
      CALL cp_fm_release(fm_mat_Q)
      IF (my_open_shell .AND. .NOT. do_im_time) THEN
         CALL cp_fm_release(fm_mat_S_beta)
         CALL cp_fm_release(fm_mat_Q_gemm_beta)
         CALL cp_fm_release(fm_mat_Q_beta)
      END IF
      IF (my_open_shell .AND. do_ri_sos_laplace_mp2 .AND. do_im_time) THEN
         CALL cp_fm_release(fm_mat_Q_beta)
      END IF

      IF (my_do_gw) THEN
         CALL cp_fm_release(fm_mat_S_gw)
         CALL cp_fm_release(fm_mat_R_gw)
         IF (my_open_shell) THEN
            CALL cp_fm_release(fm_mat_S_gw_beta)
         END IF
      END IF

      IF (do_bse) THEN
         CALL cp_fm_release(fm_mat_S_ij_bse)
         CALL cp_fm_release(fm_mat_S_ab_bse)
      END IF

      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         CALL dbcsr_release(mp2_env%ri_rpa%mo_coeff_o)
         DEALLOCATE (mp2_env%ri_rpa%mo_coeff_o)
         CALL dbcsr_release(mp2_env%ri_rpa%mo_coeff_v)
         DEALLOCATE (mp2_env%ri_rpa%mo_coeff_v)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE rpa_ri_compute_en

! **************************************************************************************************
!> \brief reorder the local data in such a way to help the next stage of matrix creation;
!>        now the data inside the group are divided into a ia x K matrix (BIb_C_2D);
!>        Subroutine created to avoid massive double coding
!> \param BIb_C_2D ...
!> \param BIb_C ...
!> \param para_env_sub ...
!> \param dimen_ia ...
!> \param homo ...
!> \param virtual ...
!> \param sizes_B_virtual ...
!> \param starts_B_virtual ...
!> \param ends_B_virtual ...
!> \param sub_proc_map ...
!> \param my_ia_size ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \author Jan Wilhelm, 03/2015
! **************************************************************************************************
   SUBROUTINE calculate_BIb_C_2D(BIb_C_2D, BIb_C, para_env_sub, dimen_ia, homo, virtual, &
                                 sizes_B_virtual, starts_B_virtual, ends_B_virtual, &
                                 sub_proc_map, my_ia_size, my_ia_start, my_ia_end, my_group_L_size)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: BIb_C_2D
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: BIb_C
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, INTENT(OUT)                               :: dimen_ia
      INTEGER, INTENT(IN)                                :: homo, virtual
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: sizes_B_virtual, starts_B_virtual, &
                                                            ends_B_virtual
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sub_proc_map
      INTEGER                                            :: my_ia_size, my_ia_start, my_ia_end, &
                                                            my_group_L_size

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_BIb_C_2D', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: occ_chunk = 128

      INTEGER :: ia_global, iiB, iproc, itmp(2), jjB, my_B_size, my_B_virtual_start, occ_high, &
         occ_low, proc_receive, proc_send, proc_shift, rec_B_size, rec_B_virtual_end, &
         rec_B_virtual_start
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_ia, sizes_ia, starts_ia
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C_rec

      dimen_ia = homo*virtual

      ALLOCATE (sizes_ia(0:para_env_sub%num_pe-1))
      sizes_ia = 0
      ALLOCATE (starts_ia(0:para_env_sub%num_pe-1))
      starts_ia = 0
      ALLOCATE (ends_ia(0:para_env_sub%num_pe-1))
      ends_ia = 0

      DO iproc = 0, para_env_sub%num_pe-1
         itmp = get_limit(dimen_ia, para_env_sub%num_pe, iproc)
         starts_ia(iproc) = itmp(1)
         ends_ia(iproc) = itmp(2)
         sizes_ia(iproc) = itmp(2)-itmp(1)+1
      END DO

      my_ia_start = starts_ia(para_env_sub%mepos)
      my_ia_end = ends_ia(para_env_sub%mepos)
      my_ia_size = sizes_ia(para_env_sub%mepos)

      my_B_size = sizes_B_virtual(para_env_sub%mepos)
      my_B_virtual_start = starts_B_virtual(para_env_sub%mepos)

      ! reorder data
      ALLOCATE (BIb_C_2D(my_ia_size, my_group_L_size))

!$OMP     PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
!$OMP              SHARED(homo,my_B_size,virtual,my_B_virtual_start,my_ia_start,my_ia_end,BIb_C,BIb_C_2D,&
!$OMP              my_group_L_size)
      DO iiB = 1, homo
         DO jjB = 1, my_B_size
            ia_global = (iiB-1)*virtual+my_B_virtual_start+jjB-1
            IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
               BIb_C_2D(ia_global-my_ia_start+1, 1:my_group_L_size) = BIb_C(1:my_group_L_size, jjB, iiB)
            END IF
         END DO
      END DO

      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_send = sub_proc_map(para_env_sub%mepos+proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos-proc_shift)

         rec_B_size = sizes_B_virtual(proc_receive)
         rec_B_virtual_end = ends_B_virtual(proc_receive)
         rec_B_virtual_start = starts_B_virtual(proc_receive)

         ! do this in chunks to avoid high memory overhead  for both BIb_C_rec and buffers in mp_sendrecv
         ! TODO: fix this more cleanly with a rewrite sending only needed data etc.
         ! TODO: occ_chunk should presumably be precomputed so that messages are limited to e.g. 100MiB.
         ALLOCATE (BIb_C_rec(my_group_L_size, rec_B_size, MIN(homo, occ_chunk)))

         DO occ_low = 1, homo, occ_chunk
            occ_high = MIN(homo, occ_low+occ_chunk-1)
            CALL mp_sendrecv(BIb_C(:, :, occ_low:occ_high), proc_send, &
                             BIb_C_rec(:, :, 1:occ_high-occ_low+1), proc_receive, &
                             para_env_sub%group)
!$OMP          PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
!$OMP                   SHARED(occ_low,occ_high,rec_B_size,virtual,rec_B_virtual_start,my_ia_start,my_ia_end,BIb_C_rec,BIb_C_2D,&
!$OMP                          my_group_L_size)
            DO iiB = occ_low, occ_high
               DO jjB = 1, rec_B_size
                  ia_global = (iiB-1)*virtual+rec_B_virtual_start+jjB-1
                  IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
                     BIb_C_2D(ia_global-my_ia_start+1, 1:my_group_L_size) = BIb_C_rec(1:my_group_L_size, jjB, iiB-occ_low+1)
                  END IF
               END DO
            END DO
         ENDDO

         DEALLOCATE (BIb_C_rec)
      END DO

      DEALLOCATE (starts_ia, ends_ia, sizes_ia)
   END SUBROUTINE calculate_BIb_C_2D

! **************************************************************************************************
!> \brief ...
!> \param BIb_C_2D ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param ngroup ...
!> \param integ_group_size ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param dimen_ia_for_block_size ...
!> \param color_rpa_group ...
!> \param ext_row_block_size ...
!> \param ext_col_block_size ...
!> \param unit_nr ...
!> \param my_ia_size ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param para_env_RPA ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_gemm ...
!> \param fm_mat_Q ...
!> \param nrow_block_mat ...
!> \param ncol_block_mat ...
!> \param beta_case ...
!> \param blacs_env_ext ...
!> \param blacs_env_ext_S ...
!> \param do_gw_open_shell ...
!> \param do_im_time ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param mo_coeff ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param nmo ...
!> \param homo ...
!> \param do_mao ...
!> \param mao_coeff_occ_A ...
!> \param mao_coeff_virt_A ...
!> \param matrix_s ...
!> \param fm_mat_Q_beta ...
! **************************************************************************************************
   SUBROUTINE create_integ_mat(BIb_C_2D, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_ia, dimen_ia_for_block_size, color_rpa_group, &
                               ext_row_block_size, ext_col_block_size, unit_nr, &
                               my_ia_size, my_ia_start, my_ia_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                               nrow_block_mat, ncol_block_mat, beta_case, &
                               blacs_env_ext, blacs_env_ext_S, do_gw_open_shell, &
                               do_im_time, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, mo_coeff, &
                               fm_mo_coeff_occ, fm_mo_coeff_virt, nmo, homo, &
                               do_mao, mao_coeff_occ_A, mao_coeff_virt_A, &
                               matrix_s, fm_mat_Q_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: BIb_C_2D
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN) :: color_sub, ngroup, integ_group_size, dimen_RI, dimen_ia, &
         dimen_ia_for_block_size, color_rpa_group, ext_row_block_size, ext_col_block_size, &
         unit_nr, my_ia_size, my_ia_start, my_ia_end, my_group_L_size, my_group_L_start, &
         my_group_L_end
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S, fm_mat_Q_gemm, fm_mat_Q
      INTEGER, INTENT(INOUT)                             :: nrow_block_mat, ncol_block_mat
      LOGICAL, INTENT(IN), OPTIONAL                      :: beta_case
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext, blacs_env_ext_S
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_gw_open_shell, do_im_time
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau, mo_coeff, &
                                                            fm_mo_coeff_occ, fm_mo_coeff_virt
      INTEGER, INTENT(IN), OPTIONAL                      :: nmo, homo
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_mao
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mao_coeff_occ_A, mao_coeff_virt_A, &
                                                            matrix_s
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_Q_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_integ_mat', &
         routineP = moduleN//':'//routineN

      INTEGER :: comm_exchange, comm_rpa, dummy_proc, end_col_block, end_row_block, grid_2D(2), &
         handle, handle2, handle3, i, i_global, i_local, i_sub, iiB, iii, iproc, iproc_col, &
         iproc_row, itmp(2), j_global, j_local, j_sub, jjB, mepos_in_RPA_group, my_num_col_blocks, &
         my_num_row_blocks, mypcol, myprow, ncol_block, ncol_local, nmao_occ, npcol, nprow, &
         nrow_block, nrow_local, num_rec_cols, number_of_rec, number_of_send, proc_receive, &
         proc_receive_static, proc_send, proc_send_static, proc_shift, rec_counter, rec_ia_end, &
         rec_ia_size, rec_ia_start, rec_L_end, rec_L_size, rec_L_start, rec_pcol
      INTEGER :: rec_prow, ref_send_pcol, ref_send_prow, row_col_proc_ratio, send_counter, &
         send_pcol, send_prow, size_rec_buffer, size_send_buffer, start_col_block, &
         start_row_block, sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii_vet, index_col_rec, map_rec_size, &
                                                            map_send_size, req_send, RPA_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: blocks_ranges_col, blocks_ranges_row, grid_2_mepos, &
         grid_ref_2_send_pos, group_grid_2_mepos, mepos_2_grid, RPA_info_end, RPA_info_size, &
         RPA_info_start
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_beta_case, my_blacs_ext, &
                                                            my_blacs_S_ext, my_do_im_time, &
                                                            my_do_mao, my_gw_open_shell
      REAL(KIND=dp)                                      :: part_ia, part_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_rec
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_Q
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      my_beta_case = .FALSE.
      IF (PRESENT(beta_case)) my_beta_case = beta_case

      my_blacs_ext = .FALSE.
      IF (PRESENT(blacs_env_ext)) my_blacs_ext = .TRUE.

      my_blacs_S_ext = .FALSE.
      IF (PRESENT(blacs_env_ext_S)) my_blacs_S_ext = .TRUE.

      my_gw_open_shell = .FALSE.
      IF (PRESENT(do_gw_open_shell)) my_gw_open_shell = do_gw_open_shell

      my_do_im_time = .FALSE.
      IF (PRESENT(do_im_time)) my_do_im_time = do_im_time

      my_do_mao = .FALSE.
      IF (PRESENT(do_mao)) my_do_mao = do_mao

      ! create the RPA para_env
      IF (.NOT. my_beta_case) THEN
         CALL mp_comm_split_direct(para_env%group, comm_rpa, color_rpa_group)
         NULLIFY (para_env_RPA)
         CALL cp_para_env_create(para_env_RPA, comm_rpa)
      END IF

      ! create the RPA blacs env
      IF (my_blacs_S_ext) THEN
         NULLIFY (blacs_env)
         blacs_env => blacs_env_ext_S
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_ia, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA)
      ELSE
         NULLIFY (blacs_env)
         IF (para_env_RPA%num_pe > 1) THEN
            row_col_proc_ratio = dimen_ia_for_block_size/dimen_RI
            row_col_proc_ratio = MAX(1, row_col_proc_ratio)
            IF (.FALSE.) THEN
               iproc_row = MAX(INT(SQRT(REAL(para_env_RPA%num_pe*row_col_proc_ratio, KIND=dp))), 1)-1
               DO iproc = 1, para_env_RPA%num_pe
                  iproc_row = iproc_row+1
                  IF (MOD(para_env_RPA%num_pe, iproc_row) == 0) EXIT
               END DO
            ELSE
               iproc_row = MIN(MAX(INT(SQRT(REAL(para_env_RPA%num_pe*row_col_proc_ratio, KIND=dp))), 1), para_env_RPA%num_pe)+1
               DO iproc = 1, para_env_RPA%num_pe
                  iproc_row = iproc_row-1
                  IF (MOD(para_env_RPA%num_pe, iproc_row) == 0) EXIT
               END DO
            END IF

            iproc_col = para_env_RPA%num_pe/iproc_row
            grid_2D(1) = iproc_row
            grid_2D(2) = iproc_col
         ELSE
            grid_2D = 1
         END IF
         CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_RPA, grid_2d=grid_2d)

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Number row processes:", grid_2D(1)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Number column processes:", grid_2D(2)
         END IF
         ! CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_RPA)

         ! define the block_size for the row
         IF (ext_row_block_size > 0) THEN
            nrow_block_mat = ext_row_block_size
         ELSE
            nrow_block_mat = dimen_ia_for_block_size/grid_2D(1)/2
            nrow_block_mat = MAX(nrow_block_mat, 1)
         END IF

         ! define the block_size for the column
         IF (ext_col_block_size > 0) THEN
            ncol_block_mat = ext_col_block_size
         ELSE
            ncol_block_mat = dimen_RI/grid_2D(2)/2
            ncol_block_mat = MAX(ncol_block_mat, 1)
         END IF

         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Row block size:", nrow_block_mat
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Column block size:", ncol_block_mat
         END IF

         NULLIFY (fm_struct)
         ! CALL cp_fm_struct_create(fm_struct,context=blacs_env,nrow_global=dimen_ia,&
         !                          ncol_global=dimen_RI,para_env=para_env_RPA)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_ia, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA, &
                                  nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)

      END IF ! external blacs_env

      ! for imaginary time we do not need the fm_mat_S
      IF (.NOT. my_do_im_time) THEN

         ! create the S full matrix, that is the (ia|K) matrix with K colomns
         ! and homo*virtual rows
         NULLIFY (fm_mat_S)
         CALL cp_fm_create(fm_mat_S, fm_struct, name="fm_mat_S")

         CALL cp_fm_set_all(matrix=fm_mat_S, alpha=0.0_dp)

         ! fill the matrix
         CALL cp_fm_get_info(matrix=fm_mat_S, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices, &
                             nrow_block=nrow_block, &
                             ncol_block=ncol_block)
         myprow = fm_mat_S%matrix_struct%context%mepos(1)
         mypcol = fm_mat_S%matrix_struct%context%mepos(2)
         nprow = fm_mat_S%matrix_struct%context%num_pe(1)
         npcol = fm_mat_S%matrix_struct%context%num_pe(2)

         ! create the RPA proc_map
         ALLOCATE (RPA_proc_map(-para_env_RPA%num_pe:2*para_env_RPA%num_pe-1))
         RPA_proc_map = 0
         DO i = 0, para_env_RPA%num_pe-1
            RPA_proc_map(i) = i
            RPA_proc_map(-i-1) = para_env_RPA%num_pe-i-1
            RPA_proc_map(para_env_RPA%num_pe+i) = i
         END DO

         ! create the info array, first index: 1-> L info, 2-> ia info
         ALLOCATE (RPA_info_start(2, 0:para_env_RPA%num_pe-1))
         RPA_info_start = 0
         ALLOCATE (RPA_info_end(2, 0:para_env_RPA%num_pe-1))
         RPA_info_end = 0
         ALLOCATE (RPA_info_size(2, 0:para_env_RPA%num_pe-1))
         RPA_info_size = 0

         RPA_info_start(1, para_env_RPA%mepos) = my_group_L_start
         RPA_info_start(2, para_env_RPA%mepos) = my_ia_start

         RPA_info_end(1, para_env_RPA%mepos) = my_group_L_end
         RPA_info_end(2, para_env_RPA%mepos) = my_ia_end

         RPA_info_size(1, para_env_RPA%mepos) = my_group_L_size
         RPA_info_size(2, para_env_RPA%mepos) = my_ia_size

         CALL mp_sum(RPA_info_start, para_env_RPA%group)
         CALL mp_sum(RPA_info_end, para_env_RPA%group)
         CALL mp_sum(RPA_info_size, para_env_RPA%group)

         IF (.TRUE.) THEN
            CALL timeset(routineN//"_comm", handle2)
            ! new replication scheme
            ! 0) create array with processes positions
            CALL timeset(routineN//"_info", handle3)
            ALLOCATE (grid_2_mepos(0:nprow-1, 0:npcol-1))
            grid_2_mepos = 0
            ALLOCATE (mepos_2_grid(0:para_env_RPA%num_pe-1, 2))
            mepos_2_grid = 0

            grid_2_mepos(myprow, mypcol) = para_env_RPA%mepos
            mepos_2_grid(para_env_RPA%mepos, 1) = myprow
            mepos_2_grid(para_env_RPA%mepos, 2) = mypcol

            CALL mp_sum(grid_2_mepos, para_env_RPA%group)
            CALL mp_sum(mepos_2_grid, para_env_RPA%group)

            ! 1) loop over my local data and define a map for the proc to send data
            ALLOCATE (map_send_size(0:para_env_RPA%num_pe-1))
            map_send_size = 0
            dummy_proc = 0
            DO jjB = my_group_L_start, my_group_L_end
               send_pcol = cp_fm_indxg2p(jjB, ncol_block, dummy_proc, &
                                         fm_mat_S%matrix_struct%first_p_pos(2), npcol)
               DO iiB = my_ia_start, my_ia_end
                  send_prow = cp_fm_indxg2p(iiB, nrow_block, dummy_proc, &
                                            fm_mat_S%matrix_struct%first_p_pos(1), nprow)
                  proc_send = grid_2_mepos(send_prow, send_pcol)
                  map_send_size(proc_send) = map_send_size(proc_send)+1
               END DO
            END DO

            ! 2) loop over my local data of fm_mat_S and define a map for the proc from which rec data
            ALLOCATE (map_rec_size(0:para_env_RPA%num_pe-1))
            map_rec_size = 0
            mepos_in_RPA_group = MOD(color_sub, integ_group_size)
            ALLOCATE (group_grid_2_mepos(0:para_env_sub%num_pe-1, 0:integ_group_size-1))
            group_grid_2_mepos = 0
            group_grid_2_mepos(para_env_sub%mepos, mepos_in_RPA_group) = para_env_RPA%mepos
            CALL mp_sum(group_grid_2_mepos, para_env_RPA%group)
            part_ia = REAL(dimen_ia, KIND=dp)/REAL(para_env_sub%num_pe, KIND=dp)
            part_RI = REAL(dimen_RI, KIND=dp)/REAL(ngroup, KIND=dp)

            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               ! check the group holding this element
               ! dirty way, if someone has a better idea ...
               rec_pcol = INT(REAL(j_global-1, KIND=dp)/part_RI)
               rec_pcol = MAX(0, rec_pcol)
               rec_pcol = MIN(rec_pcol, ngroup-1)
               DO
                  itmp = get_limit(dimen_RI, ngroup, rec_pcol)
                  IF (j_global >= itmp(1) .AND. j_global <= itmp(2)) EXIT
                  IF (j_global < itmp(1)) rec_pcol = rec_pcol-1
                  IF (j_global > itmp(2)) rec_pcol = rec_pcol+1
               END DO
               ! if the group is not in the same RPA group cycle
               IF ((rec_pcol/integ_group_size) .NE. color_rpa_group) CYCLE
               ! convert global position to position into the RPA group
               rec_pcol = MOD(rec_pcol, integ_group_size)

               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  ! check the process in the group holding this element
                  rec_prow = INT(REAL(i_global-1, KIND=dp)/part_ia)
                  rec_prow = MAX(0, rec_prow)
                  rec_prow = MIN(rec_prow, para_env_sub%num_pe-1)
                  DO
                     itmp = get_limit(dimen_ia, para_env_sub%num_pe, rec_prow)
                     IF (i_global >= itmp(1) .AND. i_global <= itmp(2)) EXIT
                     IF (i_global < itmp(1)) rec_prow = rec_prow-1
                     IF (i_global > itmp(2)) rec_prow = rec_prow+1
                  END DO

                  proc_receive = group_grid_2_mepos(rec_prow, rec_pcol)

                  map_rec_size(proc_receive) = map_rec_size(proc_receive)+1

               END DO ! i_global
            END DO ! j_global

            ! 3) check if the local data has to be stored in the new fm_mat_S
            IF (map_rec_size(para_env_RPA%mepos) > 0) THEN
               DO jjB = 1, ncol_local
                  j_global = col_indices(jjB)
                  IF (j_global >= my_group_L_start .AND. j_global <= my_group_L_end) THEN
                     DO iiB = 1, nrow_local
                        i_global = row_indices(iiB)
                        IF (i_global >= my_ia_start .AND. i_global <= my_ia_end) THEN
                           fm_mat_S%local_data(iiB, jjB) = BIb_C_2D(i_global-my_ia_start+1, j_global-my_group_L_start+1)
                        END IF
                     END DO
                  END IF
               END DO
            END IF
            CALL timestop(handle3)

            ! 4) reorder data in the send_buffer
            CALL timeset(routineN//"_buffer_s", handle3)
            number_of_send = 0
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
               IF (map_send_size(proc_send) > 0) THEN
                  number_of_send = number_of_send+1
               END IF
            END DO

            ALLOCATE (buffer_send(number_of_send))

            ! this array given a pair (ref_send_prow,ref_send_pcol) returns
            ! the position in the buffer_send associated to that process
            ALLOCATE (grid_ref_2_send_pos(0:nprow-1, 0:npcol-1))
            grid_ref_2_send_pos = 0

            send_counter = 0
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
               size_send_buffer = map_send_size(proc_send)

               IF (map_send_size(proc_send) > 0) THEN
                  send_counter = send_counter+1
                  ! prepare the sending buffer
                  ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
                  buffer_send(send_counter)%msg = 0.0_dp
                  buffer_send(send_counter)%proc = proc_send

                  ref_send_prow = mepos_2_grid(proc_send, 1)
                  ref_send_pcol = mepos_2_grid(proc_send, 2)

                  grid_ref_2_send_pos(ref_send_prow, ref_send_pcol) = send_counter

                  ! iii=0
                  ! DO iiB=my_ia_start, my_ia_end
                  !   send_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                  !                           fm_mat_S%matrix_struct%first_p_pos(1),nprow)
                  !   IF(send_prow.NE.ref_send_prow) CYCLE
                  !   DO jjB=my_group_L_start, my_group_L_end
                  !     send_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                  !                             fm_mat_S%matrix_struct%first_p_pos(2),npcol)
                  !     IF(send_pcol.NE.ref_send_pcol) CYCLE
                  !     iii=iii+1
                  !     buffer_send(send_counter)%msg(iii)=BIb_C_2D(iiB-my_ia_start+1,jjB-my_group_L_start+1)
                  !   END DO
                  ! END DO

               END IF

            END DO

            ! loop over the locally held data and fill the buffer_send
            ! for doing that we need an array index
            ALLOCATE (iii_vet(number_of_send))

            iii_vet = 0
            DO iiB = my_ia_start, my_ia_end
               send_prow = cp_fm_indxg2p(iiB, nrow_block, dummy_proc, &
                                         fm_mat_S%matrix_struct%first_p_pos(1), nprow)
               DO jjB = my_group_L_start, my_group_L_end
                  send_pcol = cp_fm_indxg2p(jjB, ncol_block, dummy_proc, &
                                            fm_mat_S%matrix_struct%first_p_pos(2), npcol)
                  ! we don't need to send to ourselves
                  IF (grid_2_mepos(send_prow, send_pcol) == para_env_RPA%mepos) CYCLE

                  send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
                  iii_vet(send_counter) = iii_vet(send_counter)+1
                  iii = iii_vet(send_counter)
                  buffer_send(send_counter)%msg(iii) = BIb_C_2D(iiB-my_ia_start+1, jjB-my_group_L_start+1)
               END DO
            END DO

            DEALLOCATE (iii_vet)
            DEALLOCATE (grid_ref_2_send_pos)
            DEALLOCATE (BIb_C_2D)
            CALL timestop(handle3)

            ! 5) create the buffer for receive and post the irecv
            CALL timeset(routineN//"_isendrecv", handle3)
            number_of_rec = 0
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
               IF (map_rec_size(proc_receive) > 0) THEN
                  number_of_rec = number_of_rec+1
               END IF
            END DO

            ALLOCATE (buffer_rec(number_of_rec))

            rec_counter = 0
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
               size_rec_buffer = map_rec_size(proc_receive)

               IF (map_rec_size(proc_receive) > 0) THEN
                  rec_counter = rec_counter+1
                  ! prepare the buffer for receive
                  ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
                  buffer_rec(rec_counter)%msg = 0.0_dp
                  buffer_rec(rec_counter)%proc = proc_receive

                  ! post the receiving message
                  CALL mp_irecv(buffer_rec(rec_counter)%msg, proc_receive, para_env_RPA%group, buffer_rec(rec_counter)%msg_req)
               END IF
            END DO

            ! 6) send data
            ALLOCATE (req_send(number_of_send))
            send_counter = 0
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
               IF (map_send_size(proc_send) > 0) THEN
                  send_counter = send_counter+1
                  CALL mp_isend(buffer_send(send_counter)%msg, proc_send, para_env_RPA%group, buffer_send(send_counter)%msg_req)
                  req_send(send_counter) = buffer_send(send_counter)%msg_req
               END IF
            END DO
            CALL timestop(handle3)

            ! 8) fill the fm_mat_S matrix
            CALL timeset(routineN//"_fill", handle3)
            ! In order to perform this step efficiently first we have to know the
            ! ranges of the blocks over which a given process hold its local data.
            ! Start with the rows ...
            my_num_row_blocks = 1
            DO iiB = 1, nrow_local-1
               IF (ABS(row_indices(iiB+1)-row_indices(iiB)) == 1) CYCLE
               my_num_row_blocks = my_num_row_blocks+1
            END DO
            ALLOCATE (blocks_ranges_row(2, my_num_row_blocks))
            blocks_ranges_row = 0
            blocks_ranges_row(1, 1) = row_indices(1)
            iii = 1
            DO iiB = 1, nrow_local-1
               IF (ABS(row_indices(iiB+1)-row_indices(iiB)) == 1) CYCLE
               iii = iii+1
               blocks_ranges_row(2, iii-1) = row_indices(iiB)
               blocks_ranges_row(1, iii) = row_indices(iiB+1)
            END DO
            blocks_ranges_row(2, my_num_row_blocks) = row_indices(MAX(nrow_local, 1))

            ! and columns
            my_num_col_blocks = 1
            DO jjB = 1, ncol_local-1
               IF (ABS(col_indices(jjB+1)-col_indices(jjB)) == 1) CYCLE
               my_num_col_blocks = my_num_col_blocks+1
            END DO
            ALLOCATE (blocks_ranges_col(2, my_num_col_blocks))
            blocks_ranges_col = 0
            blocks_ranges_col(1, 1) = col_indices(1)
            iii = 1
            DO jjB = 1, ncol_local-1
               IF (ABS(col_indices(jjB+1)-col_indices(jjB)) == 1) CYCLE
               iii = iii+1
               blocks_ranges_col(2, iii-1) = col_indices(jjB)
               blocks_ranges_col(1, iii) = col_indices(jjB+1)
            END DO
            blocks_ranges_col(2, my_num_col_blocks) = col_indices(MAX(ncol_local, 1))

            rec_counter = 0
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)
               size_rec_buffer = map_rec_size(proc_receive)

               IF (map_rec_size(proc_receive) > 0) THEN
                  rec_counter = rec_counter+1

                  rec_L_size = RPA_info_size(1, proc_receive)
                  rec_L_start = RPA_info_start(1, proc_receive)
                  rec_L_end = RPA_info_end(1, proc_receive)
                  ! precompute the number of received columns and relative index
                  num_rec_cols = 0
                  DO jjB = 1, my_num_col_blocks
                     start_col_block = MAX(blocks_ranges_col(1, jjB), rec_L_start)
                     end_col_block = MIN(blocks_ranges_col(2, jjB), rec_L_end)
                     DO j_sub = start_col_block, end_col_block
                        num_rec_cols = num_rec_cols+1
                     END DO
                  END DO
                  ALLOCATE (index_col_rec(num_rec_cols))
                  index_col_rec = 0
                  iii = 0
                  DO jjB = 1, my_num_col_blocks
                     start_col_block = MAX(blocks_ranges_col(1, jjB), rec_L_start)
                     end_col_block = MIN(blocks_ranges_col(2, jjB), rec_L_end)
                     DO j_sub = start_col_block, end_col_block
                        iii = iii+1
                        j_local = cp_fm_indxg2l(j_sub, ncol_block, dummy_proc, &
                                                fm_mat_S%matrix_struct%first_p_pos(2), npcol)
                        index_col_rec(iii) = j_local
                     END DO
                  END DO

                  rec_ia_size = RPA_info_size(2, proc_receive)
                  rec_ia_start = RPA_info_start(2, proc_receive)
                  rec_ia_end = RPA_info_end(2, proc_receive)

                  ! wait for the message
                  CALL mp_wait(buffer_rec(rec_counter)%msg_req)

                  iii = 0
                  DO iiB = 1, my_num_row_blocks
                     start_row_block = MAX(blocks_ranges_row(1, iiB), rec_ia_start)
                     end_row_block = MIN(blocks_ranges_row(2, iiB), rec_ia_end)
                     DO i_sub = start_row_block, end_row_block
                        i_local = cp_fm_indxg2l(i_sub, nrow_block, dummy_proc, &
                                                fm_mat_S%matrix_struct%first_p_pos(1), nprow)
                        ! DO jjB=1, my_num_col_blocks
                        !   start_col_block=MAX(blocks_ranges_col(1,jjB),rec_L_start)
                        !   end_col_block=MIN(blocks_ranges_col(2,jjB),rec_L_end)
                        !   DO j_sub=start_col_block, end_col_block
                        !     j_local=cp_fm_indxg2l(j_sub,ncol_block,dummy_proc,&
                        !                           fm_mat_S%matrix_struct%first_p_pos(2),npcol)
                        !     iii=iii+1
                        !     fm_mat_S%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
                        !   END DO
                        ! END DO
                        DO jjB = 1, num_rec_cols
                           iii = iii+1
                           j_local = index_col_rec(jjB)
                           fm_mat_S%local_data(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
                        END DO
                     END DO
                  END DO

                  ! iii=0
                  ! DO iiB=rec_ia_start, rec_ia_end
                  !   rec_prow=cp_fm_indxg2p(iiB,nrow_block,dummy_proc,&
                  !                    fm_mat_S%matrix_struct%first_p_pos(1),nprow)
                  !   IF(rec_prow.NE.myprow) CYCLE
                  !   i_local=cp_fm_indxg2l(iiB,nrow_block,dummy_proc,&
                  !                         fm_mat_S%matrix_struct%first_p_pos(1),nprow)
                  !   DO jjB=rec_L_start, rec_L_end
                  !     rec_pcol=cp_fm_indxg2p(jjB,ncol_block,dummy_proc,&
                  !                      fm_mat_S%matrix_struct%first_p_pos(2),npcol)
                  !     IF(rec_pcol.NE.mypcol) CYCLE
                  !     j_local=cp_fm_indxg2l(jjB,ncol_block,dummy_proc,&
                  !                     fm_mat_S%matrix_struct%first_p_pos(2),npcol)
                  !     iii=iii+1
                  !     fm_mat_S%local_data(i_local,j_local)=buffer_rec(rec_counter)%msg(iii)
                  !   END DO
                  ! END DO

                  DEALLOCATE (buffer_rec(rec_counter)%msg)
                  DEALLOCATE (index_col_rec)

               END IF
            END DO
            DEALLOCATE (buffer_rec)

            DEALLOCATE (blocks_ranges_row)
            DEALLOCATE (blocks_ranges_col)

            CALL timestop(handle3)

            ! wait for all messeges to be sent
            CALL timeset(routineN//"_waitall", handle3)
            CALL mp_waitall(req_send(:))
            DO send_counter = 1, number_of_send
               DEALLOCATE (buffer_send(send_counter)%msg)
            END DO
            DEALLOCATE (buffer_send)
            CALL timestop(handle3)

            DEALLOCATE (group_grid_2_mepos)
            DEALLOCATE (map_send_size)
            DEALLOCATE (map_rec_size)
            DEALLOCATE (grid_2_mepos)
            DEALLOCATE (mepos_2_grid)

            CALL timestop(handle2)
         ELSE
            ! old inefficient replication scheme
            ! local data
            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               IF (j_global >= my_group_L_start .AND. j_global <= my_group_L_end) THEN
                  DO iiB = 1, nrow_local
                     i_global = row_indices(iiB)
                     IF (i_global >= my_ia_start .AND. i_global <= my_ia_end) THEN
                        fm_mat_S%local_data(iiB, jjB) = BIb_C_2D(i_global-my_ia_start+1, j_global-my_group_L_start+1)
                     END IF
                  END DO
               END IF
            END DO

            proc_send_static = RPA_proc_map(para_env_RPA%mepos+1)
            proc_receive_static = RPA_proc_map(para_env_RPA%mepos-1)

            ! start collect data from other proc in the RPA group
            DO proc_shift = 1, para_env_RPA%num_pe-1
               proc_send = RPA_proc_map(para_env_RPA%mepos+proc_shift)
               proc_receive = RPA_proc_map(para_env_RPA%mepos-proc_shift)

               rec_L_size = RPA_info_size(1, proc_receive)
               rec_L_start = RPA_info_start(1, proc_receive)
               rec_L_end = RPA_info_end(1, proc_receive)

               rec_ia_size = RPA_info_size(2, proc_receive)
               rec_ia_start = RPA_info_start(2, proc_receive)
               rec_ia_end = RPA_info_end(2, proc_receive)

               ALLOCATE (BIb_C_rec(rec_ia_size, rec_L_size))
               BIb_C_rec = 0.0_dp

               CALL mp_sendrecv(BIb_C_2D, proc_send_static, &
                                BIb_C_rec, proc_receive_static, &
                                para_env_RPA%group)

               DO jjB = 1, ncol_local
                  j_global = col_indices(jjB)
                  IF (j_global >= rec_L_start .AND. j_global <= rec_L_end) THEN
                     DO iiB = 1, nrow_local
                        i_global = row_indices(iiB)
                        IF (i_global >= rec_ia_start .AND. i_global <= rec_ia_end) THEN
                           fm_mat_S%local_data(iiB, jjB) = BIb_C_rec(i_global-rec_ia_start+1, j_global-rec_L_start+1)
                        END IF
                     END DO
                  END IF
               END DO

               DEALLOCATE (BIb_C_2D)
               ALLOCATE (BIb_C_2D(rec_ia_size, rec_L_size))
               BIb_C_2D(:, :) = BIb_C_rec

               DEALLOCATE (BIb_C_rec)
            END DO

            DEALLOCATE (BIb_C_2D)
         END IF

         ! deallocaete the info array
         DEALLOCATE (RPA_info_start)
         DEALLOCATE (RPA_info_end)
         DEALLOCATE (RPA_info_size)

         ! mp_sum the local data across processes belonging to different RPA group.
         ! first create the para_env then mp_sum
         sub_sub_color = para_env_RPA%mepos
         CALL mp_comm_split_direct(para_env%group, comm_exchange, sub_sub_color)
         NULLIFY (para_env_exchange)
         CALL cp_para_env_create(para_env_exchange, comm_exchange)

         CALL timeset(routineN//"_sum", handle2)
         CALL mp_sum(fm_mat_S%local_data, para_env_exchange%group)
         CALL timestop(handle2)

         CALL cp_para_env_release(para_env_exchange)

      END IF ! not imag. time

      CALL cp_fm_struct_release(fm_struct)

      IF (.NOT. my_gw_open_shell) THEN

         IF (.NOT. my_do_im_time) THEN
            ! create the Q matrix dimen_RIxdimen_RI where the result of the mat-mat-mult will be stored
            NULLIFY (fm_mat_Q_gemm)
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                                     ncol_global=dimen_RI, para_env=para_env_RPA, &
                                     nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)
            CALL cp_fm_create(fm_mat_Q_gemm, fm_struct, name="fm_mat_Q_gemm")
            CALL cp_fm_struct_release(fm_struct)

            CALL cp_fm_set_all(matrix=fm_mat_Q_gemm, alpha=0.0_dp)
         END IF

         ! create the Q matrix with a different blacs env
         NULLIFY (blacs_env_Q)
         IF (my_blacs_ext) THEN
            blacs_env_Q => blacs_env_ext
         ELSE
            CALL cp_blacs_env_create(blacs_env=blacs_env_Q, para_env=para_env_RPA)
         END IF

         NULLIFY (fm_mat_Q)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_Q, nrow_global=dimen_RI, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA)
         CALL cp_fm_create(fm_mat_Q, fm_struct, name="fm_mat_Q")

         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_set_all(matrix=fm_mat_Q, alpha=0.0_dp)

         IF (PRESENT(fm_mat_Q_beta)) THEN
            NULLIFY (fm_mat_Q_beta)
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, context=blacs_env_Q, nrow_global=dimen_RI, &
                                     ncol_global=dimen_RI, para_env=para_env_RPA)
            CALL cp_fm_create(fm_mat_Q_beta, fm_struct, name="fm_mat_Q_beta")

            CALL cp_fm_struct_release(fm_struct)

            CALL cp_fm_set_all(matrix=fm_mat_Q_beta, alpha=0.0_dp)
         END IF
      END IF

      ! in case we do imaginary time, we allocate fm_scaled_dm_occ and fm_scaled_dm_virt
      IF (my_do_im_time) THEN

         IF (my_do_mao) THEN
            CALL dbcsr_get_info(mao_coeff_occ_A(1)%matrix, &
                                nfullcols_total=nmao_occ)
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, context=mo_coeff%matrix_struct%context, nrow_global=nmao_occ, &
                                     ncol_global=nmao_occ, para_env=mo_coeff%matrix_struct%para_env)
         END IF

         IF (my_do_mao) THEN
            CALL cp_fm_create(fm_scaled_dm_occ_tau, fm_struct)
         ELSE
            ! default: no MAOs
            CALL cp_fm_create(fm_scaled_dm_occ_tau, mo_coeff%matrix_struct)
         END IF
         CALL cp_fm_set_all(fm_scaled_dm_occ_tau, 0.0_dp)

         CALL cp_fm_create(fm_scaled_dm_virt_tau, mo_coeff%matrix_struct)
         CALL cp_fm_set_all(fm_scaled_dm_virt_tau, 0.0_dp)

         CALL create_occ_virt_mo_coeffs(fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, &
                                        nmo, homo, my_do_mao, mao_coeff_occ_A, &
                                        mao_coeff_virt_A, matrix_s, 1)

         IF (my_do_mao) THEN
            CALL cp_fm_struct_release(fm_struct)
         END IF

      END IF ! imag. time

      ! release blacs_env
      IF (.NOT. my_blacs_ext) CALL cp_blacs_env_release(blacs_env_Q)
      IF (.NOT. my_blacs_S_ext) CALL cp_blacs_env_release(blacs_env)

      CALL timestop(handle)

   END SUBROUTINE create_integ_mat

! **************************************************************************************************
!> \brief ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param mo_coeff ...
!> \param nmo ...
!> \param homo ...
!> \param do_mao ...
!> \param mao_coeff_occ_A ...
!> \param mao_coeff_virt_A ...
!> \param matrix_s ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE create_occ_virt_mo_coeffs(fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, &
                                        nmo, homo, do_mao, mao_coeff_occ_A, &
                                        mao_coeff_virt_A, matrix_s, ispin)
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            mo_coeff
      INTEGER, INTENT(IN)                                :: nmo, homo
      LOGICAL, INTENT(IN)                                :: do_mao
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mao_coeff_occ_A, mao_coeff_virt_A, &
                                                            matrix_s
      INTEGER, INTENT(IN)                                :: ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_occ_virt_mo_coeffs', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icol_global, irow_global

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mo_coeff_occ, mo_coeff%matrix_struct)
      CALL cp_fm_set_all(fm_mo_coeff_occ, 0.0_dp)
      CALL cp_fm_to_fm(mo_coeff, fm_mo_coeff_occ)

      ! set all virtual MO coeffs to zero
      DO irow_global = 1, nmo
         DO icol_global = homo+1, nmo
            CALL cp_fm_set_element(fm_mo_coeff_occ, irow_global, icol_global, 0.0_dp)
         END DO
      END DO

      CALL cp_fm_create(fm_mo_coeff_virt, mo_coeff%matrix_struct)
      CALL cp_fm_set_all(fm_mo_coeff_virt, 0.0_dp)
      CALL cp_fm_to_fm(mo_coeff, fm_mo_coeff_virt)

      ! set all occupied MO coeffs to zero
      DO irow_global = 1, nmo
         DO icol_global = 1, homo
            CALL cp_fm_set_element(fm_mo_coeff_virt, irow_global, icol_global, 0.0_dp)
         END DO
      END DO

      ! transform the AO index of fm_mo_coeff_occ to the MAO index
      IF (do_mao) THEN

         CALL transform_MO_coeff_to_MAO_basis(mao_coeff_occ_A, matrix_s, ispin, nmo, fm_mo_coeff_occ)
         CALL transform_MO_coeff_to_MAO_basis(mao_coeff_virt_A, matrix_s, ispin, nmo, fm_mo_coeff_virt)

      END IF

      CALL timestop(handle)

   END SUBROUTINE create_occ_virt_mo_coeffs

! **************************************************************************************************
!> \brief ...
!> \param mao_coeff_A ...
!> \param matrix_s ...
!> \param ispin ...
!> \param nmo ...
!> \param fm_mo_coeff ...
! **************************************************************************************************
   SUBROUTINE transform_MO_coeff_to_MAO_basis(mao_coeff_A, matrix_s, ispin, nmo, fm_mo_coeff)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff_A, matrix_s
      INTEGER, INTENT(IN)                                :: ispin, nmo
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_MO_coeff_to_MAO_basis', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nmao
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff, mo_coeff_mao_basis

      CALL dbcsr_get_info(mao_coeff_A(ispin)%matrix, &
                          row_blk_size=row_blk_sizes, &
                          col_blk_size=col_blk_sizes, &
                          nfullcols_total=nmao)

      NULLIFY (mo_coeff_mao_basis)
      CALL dbcsr_init_p(mo_coeff_mao_basis)
      CALL dbcsr_create(matrix=mo_coeff_mao_basis, &
                        template=mao_coeff_A(ispin)%matrix, &
                        row_blk_size=col_blk_sizes, &
                        col_blk_size=row_blk_sizes, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_mo_coeff)
      CALL dbcsr_init_p(mat_mo_coeff)
      CALL dbcsr_create(matrix=mat_mo_coeff, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mo_coeff, mat_mo_coeff, keep_sparsity=.FALSE.)

      CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coeff_A(ispin)%matrix, mat_mo_coeff, &
                          0.0_dp, mo_coeff_mao_basis)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=fm_mo_coeff%matrix_struct%context, nrow_global=nmao, &
                               ncol_global=nmo, para_env=fm_mo_coeff%matrix_struct%para_env)

      ! get net fm_mo_coeff in the MAO basis
      CALL cp_fm_release(fm_mo_coeff)
      CALL cp_fm_create(fm_mo_coeff, fm_struct, name="mo_coeffs_in_the_mao_basis")

      CALL copy_dbcsr_to_fm(mo_coeff_mao_basis, fm_mo_coeff)

      CALL cp_fm_struct_release(fm_struct)
      CALL dbcsr_release_p(mo_coeff_mao_basis)
      CALL dbcsr_release_p(mat_mo_coeff)

   END SUBROUTINE transform_MO_coeff_to_MAO_basis

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param Erpa ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param para_env_sub ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param dimen_nm_gw ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param fm_mat_S ...
!> \param fm_mat_Q_gemm ...
!> \param fm_mat_Q ...
!> \param fm_mat_S_gw ...
!> \param fm_mat_R_gw ...
!> \param fm_mat_S_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param my_do_gw ...
!> \param do_bse ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param do_minimax_quad ...
!> \param do_im_time ...
!> \param do_mao ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param mo_coeff ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param mat_munu ...
!> \param mat_dm_occ_local ...
!> \param mat_dm_virt_local ...
!> \param mat_P_local ...
!> \param mat_P_global ...
!> \param mat_M ...
!> \param mat_3c_overl_int ...
!> \param mat_3c_overl_int_mao_for_occ ...
!> \param mat_3c_overl_int_mao_for_virt ...
!> \param matrix_s ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param eps_filter ...
!> \param starts_array ...
!> \param ends_array ...
!> \param sizes_array ...
!> \param color_sub ...
!> \param fm_mo_coeff_occ_beta ...
!> \param fm_mo_coeff_virt_beta ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
!> \param fm_mat_Q_gemm_beta ...
!> \param fm_mat_Q_beta ...
!> \param fm_mat_S_gw_beta ...
!> \param gw_corr_lev_occ_beta ...
!> \param gw_corr_lev_virt_beta ...
!> \param mo_coeff_beta ...
!> \param do_ri_sos_laplace_mp2 ...
! **************************************************************************************************
   SUBROUTINE rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                          homo, virtual, dimen_RI, dimen_ia, dimen_nm_gw, &
                          Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                          fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                          fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                          my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                          do_minimax_quad, do_im_time, do_mao, fm_mo_coeff_occ, &
                          fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                          fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, mat_munu, mat_dm_occ_local, &
                          mat_dm_virt_local, mat_P_local, &
                          mat_P_global, mat_M, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                          mat_3c_overl_int_mao_for_virt, matrix_s, mao_coeff_occ, mao_coeff_virt, &
                          eps_filter, starts_array, ends_array, sizes_array, color_sub, &
                          fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, &
                          homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta, &
                          fm_mat_Q_gemm_beta, fm_mat_Q_beta, fm_mat_S_gw_beta, &
                          gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, mo_coeff_beta, do_ri_sos_laplace_mp2)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: Erpa
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA, para_env_sub
      INTEGER, INTENT(IN)                                :: unit_nr, homo, virtual, dimen_RI, &
                                                            dimen_ia, dimen_nm_gw
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                                                            fm_mat_S_gw, fm_mat_R_gw, &
                                                            fm_mat_S_ij_bse, fm_mat_S_ab_bse
      LOGICAL, INTENT(IN)                                :: my_do_gw, do_bse
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, INTENT(IN)                                :: do_minimax_quad, do_im_time, do_mao
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(cp_fm_type), POINTER                          :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu, mat_dm_occ_local, &
                                                            mat_dm_virt_local, mat_P_local, &
                                                            mat_P_global, mat_M
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, mao_coeff_occ, mao_coeff_virt
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: starts_array, ends_array, sizes_array
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mo_coeff_occ_beta, &
                                                            fm_mo_coeff_virt_beta
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta, virtual_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_beta, fm_mat_Q_gemm_beta, &
                                                            fm_mat_Q_beta, fm_mat_S_gw_beta
      INTEGER, OPTIONAL                                  :: gw_corr_lev_occ_beta, &
                                                            gw_corr_lev_virt_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mo_coeff_beta
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_num_int', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vec_Sigma_c_gw, vec_Sigma_c_gw_beta
      INTEGER :: col_start_local, color_sub_col, color_sub_row, count_ev_sc_GW, cut_memory, &
         cut_RI, group_size_P, gw_corr_lev_tot, handle, handle3, handle4, i_cell, i_cut_RI, i_dim, &
         i_global, i_kp, i_mem, i_size, ierr, iiB, ikp, info_chol, iter_ev_sc, j_cell, j_global, &
         j_mem, j_size, jjB, jquad, LLL, max_iter_bse, mm_style, my_group_L_size_im_time, &
         my_num_dgemm_call, n_group_col, n_group_row, n_local_col, n_local_row, nblkrows_total, &
         ncol_local, ngroup_RI_orig, nkp, nmo, nrow_local, num_3c_repl, num_cells_dm, &
         num_fit_points, num_points_corr, num_points_per_magnitude, num_Z_vectors, number_of_rec
      INTEGER :: number_of_rec_axk, number_of_rec_beta, number_of_send, number_of_send_axk, &
         number_of_send_beta, row, row_start_local, size_P
      INTEGER, ALLOCATABLE, DIMENSION(:) :: map_rec_size, map_rec_size_axk, map_rec_size_beta, &
         map_send_size, map_send_size_axk, map_send_size_beta, mepos_P_from_RI_row, &
         my_group_L_sizes_im_time, my_group_L_starts_im_time, row_from_LLL, RPA_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ends_array_prim_col, ends_array_prim_fullcol, &
         ends_array_prim_fullrow, ends_array_prim_row, index_to_cell_3c, local_size_source, &
         local_size_source_axk, local_size_source_beta, sizes_array_prim_col, &
         sizes_array_prim_row, starts_array_prim_col, starts_array_prim_fullcol, &
         starts_array_prim_fullrow, starts_array_prim_row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index_3c, non_zero_blocks_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: non_zero_blocks_3c_cut_col
      INTEGER, DIMENSION(3)                              :: cell_grid_dm
      INTEGER, DIMENSION(:), POINTER :: col_blk_size, col_indices, ends_array_cm, &
         ends_array_cm_mao_occ, ends_array_cm_mao_virt, prim_blk_sizes, RI_blk_sizes, &
         row_blk_offset, row_blk_size, row_indices, starts_array_cm, starts_array_cm_mao_occ, &
         starts_array_cm_mao_virt
      LOGICAL :: do_apply_ic_corr_to_gw, do_gw_im_time, do_ic_model, do_ic_opt_homo_lumo, &
         do_kpoints_cubic_RPA, do_periodic, do_ri_Sigma_x, first_cycle, &
         first_cycle_periodic_correction, my_open_shell
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: do_GW_corr
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: cycle_due_to_sparse_dm, &
                                                            multiply_needed_occ, &
                                                            multiply_needed_virt
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: needed_cutRI_mem_R1vec_R2vec_for_kp
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :, :)     :: has_mat_P_blocks
      REAL(KIND=dp) :: a_scaling, alpha, e_axk, e_axk_corr, e_fermi, e_fermi_beta, E_Range, Emax, &
         Emax_beta, Emin, Emin_beta, eps_filter_im_time, eps_min_trans, FComega, &
         fermi_level_offset, max_error_min, my_flop_rate, omega, omega_old, scaling, &
         stabilize_exp, tau, tau_old, trace_XX
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: delta_corr, Eigenval_last, Eigenval_last_beta, &
         Eigenval_scf, Eigenval_scf_beta, m_value, m_value_beta, Q_log, tau_tj, tau_wj, tj, &
         trace_Qomega, vec_gw_energ, vec_gw_energ_beta, vec_gw_energ_error_fit, &
         vec_gw_energ_error_fit_beta, vec_omega_fit_gw, vec_Sigma_x_gw, vec_Sigma_x_gw_beta, &
         vec_W_gw, vec_W_gw_beta, wj, x_tw, z_value, z_value_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_mat_M, weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_iaQ_bse_local
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ic_corr_list, ic_corr_list_beta
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W_tau
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_type), POINTER :: fm_mat_L_transposed, fm_mat_Q_static_bse, &
         fm_mat_Q_static_bse_gemm, fm_mat_S_gw_work, fm_mat_S_gw_work_beta, fm_mat_work, &
         fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type)                                 :: mat_dm, mat_L, mat_M_P_munu_occ, &
                                                            mat_M_P_munu_virt, mat_P_global_copy, &
                                                            mat_SinvVSinv
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: mat_3c_overl_int_gw, mat_3c_overl_int_gw_beta, &
         mat_3c_overl_nnP_ic, mat_3c_overl_nnP_ic_beta, mat_3c_overl_nnP_ic_reflected, &
         mat_3c_overl_nnP_ic_reflected_beta, mat_greens_fct_occ, mat_greens_fct_occ_beta, &
         mat_greens_fct_virt, mat_greens_fct_virt_beta, mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, &
         mat_W, matrix_berry_im_mo_mo, matrix_berry_re_mo_mo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega, mat_P_omega_beta, &
                                                            mat_P_omega_im_part
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_dm_loc_occ_cut, mat_dm_loc_virt_cut
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut, &
         mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut
      TYPE(dbcsr_type), POINTER                          :: mat_contr_gf_occ, mat_contr_gf_virt, &
                                                            mat_contr_W, mat_dm_loc_occ, &
                                                            mat_dm_loc_virt
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_rec_axk, &
                                                            buffer_rec_beta, buffer_send, &
                                                            buffer_send_axk, buffer_send_beta
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(two_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: offset_combi_block

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(homo_beta) .AND. &
          PRESENT(virtual_beta) .AND. &
          PRESENT(dimen_ia_beta) .AND. &
          PRESENT(Eigenval_beta) .AND. &
          PRESENT(fm_mat_S_beta) .AND. &
          PRESENT(fm_mat_Q_gemm_beta) .AND. &
          PRESENT(fm_mat_Q_beta)) my_open_shell = .TRUE.

      nmo = homo+virtual

      do_gw_im_time = mp2_env%ri_rpa_im_time%do_gw_im_time
      do_ri_Sigma_x = mp2_env%ri_g0w0%do_ri_Sigma_x
      do_ic_model = mp2_env%ri_g0w0%do_ic_model
      do_ic_opt_homo_lumo = mp2_env%ri_g0w0%do_opt_homo_lumo
      do_periodic = mp2_env%ri_g0w0%do_periodic
      ic_corr_list => mp2_env%ri_g0w0%ic_corr_list
      ic_corr_list_beta => mp2_env%ri_g0w0%ic_corr_list_beta
      do_kpoints_cubic_RPA = mp2_env%ri_rpa_im_time%do_im_time_kpoints
      mm_style = mp2_env%ri_rpa%mm_style
      fermi_level_offset = mp2_env%ri_g0w0%fermi_level_offset

      IF (do_kpoints_cubic_RPA .AND. do_ri_sos_laplace_mp2) THEN
         CPABORT("RI-SOS-Laplace-MP2 with k-point-sampling is not implemented.")
      END IF

      do_apply_ic_corr_to_gw = .FALSE.
      IF (ic_corr_list(1) > 0.0_dp) do_apply_ic_corr_to_gw = .TRUE.

      ALLOCATE (tj(num_integ_points))
      tj = 0.0_dp

      ALLOCATE (wj(num_integ_points))
      wj = 0.0_dp

      ALLOCATE (Q_log(dimen_RI))

      ALLOCATE (trace_Qomega(dimen_RI))

      IF (do_im_time) THEN
         ! imag. time RPA only with Minimax
         CPASSERT(do_minimax_quad .OR. do_ri_sos_laplace_mp2)
      END IF

      IF (do_ic_model) THEN
         ! image charge model only implemented for cubic scaling GW
         CPASSERT(do_gw_im_time)
         CPASSERT(.NOT. do_periodic)
      END IF

      ! set up the least-square time grid and other matrices specifically for imag time
      IF (do_im_time) THEN

         CALL timeset(routineN//"_im_t_alloc_mat", handle4)

         ALLOCATE (tau_tj(num_integ_points))
         tau_tj = 0.0_dp

         ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
         weights_cos_tf_t_to_w = 0.0_dp

         ALLOCATE (tau_wj(num_integ_points))
         tau_wj = 0.0_dp

         group_size_P = mp2_env%ri_rpa_im_time%group_size_P
         cut_memory = mp2_env%ri_rpa_im_time%cut_memory
         cut_RI = mp2_env%ri_rpa_im_time_util(1)%cut_RI
         eps_filter_im_time = mp2_env%ri_rpa_im_time%eps_filter_im_time
         ngroup_RI_orig = mp2_env%ri_rpa_im_time_util(1)%n_group_RI_orig
         stabilize_exp = mp2_env%ri_rpa_im_time%stabilize_exp

         ALLOCATE (my_group_L_starts_im_time(cut_RI))
         my_group_L_starts_im_time(:) = mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time
         ALLOCATE (my_group_L_sizes_im_time(cut_RI))
         my_group_L_sizes_im_time(:) = mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time

         num_3c_repl = SIZE(mat_3c_overl_int, 3)

         DO i_cut_RI = 1, cut_RI

            my_group_L_size_im_time = my_group_L_sizes_im_time(i_cut_RI)

            DO i_cell = 1, num_3c_repl
               DO j_cell = 1, num_3c_repl

                  CALL dbcsr_filter(mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix, &
                                    eps_filter*SQRT(REAL(my_group_L_size_im_time, KIND=dp)))

               END DO
            END DO

         END DO

         CALL timestop(handle4)

         CALL timeset(routineN//"_im_t_alloc_mat_1", handle4)

         NULLIFY (mat_dm_loc_occ)
         CALL dbcsr_init_p(mat_dm_loc_occ)
         CALL dbcsr_desymmetrize(mat_dm_occ_local%matrix, mat_dm_loc_occ)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_2", handle4)

         CALL get_non_zero_blocks_3c(mat_3c_overl_int, para_env_sub, cut_RI, non_zero_blocks_3c)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_3", handle4)

         IF (do_kpoints_cubic_RPA) THEN
            CALL get_qs_env(qs_env, &
                            kpoints=kpoints)
            ! we always use an odd number of image cells
            ! CAUTION: also at another point, cell_grid_dm is defined, these definitions have to be identical
            DO i_dim = 1, 3
               cell_grid_dm(i_dim) = (kpoints%nkp_grid(i_dim)/2)*2-1
            END DO
            num_cells_dm = cell_grid_dm(1)*cell_grid_dm(2)*cell_grid_dm(3)
            ALLOCATE (index_to_cell_3c(3, SIZE(kpoints%index_to_cell, 2)))
            CPASSERT(SIZE(kpoints%index_to_cell, 1) == 3)
            index_to_cell_3c(:, :) = kpoints%index_to_cell(:, :)
            ALLOCATE (cell_to_index_3c(LBOUND(kpoints%cell_to_index, 1):UBOUND(kpoints%cell_to_index, 1), &
                                       LBOUND(kpoints%cell_to_index, 2):UBOUND(kpoints%cell_to_index, 2), &
                                       LBOUND(kpoints%cell_to_index, 3):UBOUND(kpoints%cell_to_index, 3)))
            cell_to_index_3c(:, :, :) = kpoints%cell_to_index(:, :, :)

            NULLIFY (cfm_mat_Q)
            CALL cp_cfm_create(cfm_mat_Q, fm_mat_Q%matrix_struct)
            CALL cp_cfm_set_all(cfm_mat_Q, z_zero)

         ELSE
            ALLOCATE (index_to_cell_3c(3, 1))
            index_to_cell_3c(:, 1) = 0
            ALLOCATE (cell_to_index_3c(0:0, 0:0, 0:0))
            cell_to_index_3c(0, 0, 0) = 1
            num_cells_dm = 1
         END IF

         NULLIFY (mat_dm_loc_occ_cut)
         CALL dbcsr_allocate_matrix_set(mat_dm_loc_occ_cut, cut_RI, cut_memory, num_cells_dm)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_4", handle4)

         DO i_mem = 1, cut_memory
            DO i_cut_RI = 1, cut_RI
               DO i_cell = 1, num_cells_dm

                  ALLOCATE (mat_dm_loc_occ_cut(i_cut_RI, i_mem, i_cell)%matrix)
                  CALL dbcsr_create(matrix=mat_dm_loc_occ_cut(i_cut_RI, i_mem, i_cell)%matrix, &
                                    template=mat_dm_loc_occ)

               END DO
            END DO
         END DO

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_5", handle4)

         NULLIFY (mat_dm_loc_virt)
         CALL dbcsr_init_p(mat_dm_loc_virt)
         CALL dbcsr_desymmetrize(mat_dm_virt_local%matrix, mat_dm_loc_virt)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_6", handle4)

         NULLIFY (mat_dm_loc_virt_cut)
         CALL dbcsr_allocate_matrix_set(mat_dm_loc_virt_cut, cut_RI, cut_memory, num_cells_dm)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_7", handle4)

         DO i_mem = 1, cut_memory
            DO i_cut_RI = 1, cut_RI
               DO i_cell = 1, num_cells_dm

                  ALLOCATE (mat_dm_loc_virt_cut(i_cut_RI, i_mem, i_cell)%matrix)
                  CALL dbcsr_create(matrix=mat_dm_loc_virt_cut(i_cut_RI, i_mem, i_cell)%matrix, &
                                    template=mat_dm_loc_virt)

               END DO
            END DO
         END DO

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_8", handle4)

         CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
         CALL dbcsr_filter(mat_munu%matrix, 1.0_dp)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_9", handle4)

         NULLIFY (mat_M_P_munu_occ%matrix)
         ALLOCATE (mat_M_P_munu_occ%matrix)
         CALL dbcsr_create(mat_M_P_munu_occ%matrix, template=mat_M%matrix)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_10", handle4)

         NULLIFY (mat_M_P_munu_virt%matrix)
         ALLOCATE (mat_M_P_munu_virt%matrix)
         CALL dbcsr_create(mat_M_P_munu_virt%matrix, template=mat_M%matrix)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_11", handle4)

         IF (.NOT. do_mao) THEN
            mat_3c_overl_int_mao_for_occ => mat_3c_overl_int
            mat_3c_overl_int_mao_for_virt => mat_3c_overl_int
         END IF

         NULLIFY (mat_M_mu_Pnu_occ)
         CALL dbcsr_allocate_matrix_set(mat_M_mu_Pnu_occ, cut_RI)
         DO i_cut_RI = 1, cut_RI
            ALLOCATE (mat_M_mu_Pnu_occ(i_cut_RI)%matrix)
            CALL dbcsr_create(matrix=mat_M_mu_Pnu_occ(i_cut_RI)%matrix, &
                              template=mat_3c_overl_int_mao_for_occ(i_cut_RI, 1, 1)%matrix)
         END DO

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_12", handle4)

         NULLIFY (mat_M_mu_Pnu_virt)
         CALL dbcsr_allocate_matrix_set(mat_M_mu_Pnu_virt, cut_RI)
         DO i_cut_RI = 1, cut_RI
            ALLOCATE (mat_M_mu_Pnu_virt(i_cut_RI)%matrix)
            CALL dbcsr_create(matrix=mat_M_mu_Pnu_virt(i_cut_RI)%matrix, &
                              template=mat_3c_overl_int_mao_for_virt(i_cut_RI, 1, 1)%matrix)
         END DO

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_13", handle4)

         ! if we do kpoints, mat_P has a kpoint and mat_P_omega has the inted
         ! mat_P(tau, kpoint)
         IF (do_kpoints_cubic_RPA) THEN
            CALL get_kpoint_info(kpoints, nkp=nkp)
            size_P = MAX(num_cells_dm/2+1, nkp)
         ELSE
            size_P = 1
         END IF

         NULLIFY (mat_P_omega)
         CALL dbcsr_allocate_matrix_set(mat_P_omega, num_integ_points, size_P)
         DO jquad = 1, num_integ_points
            DO i_kp = 1, size_P
               ALLOCATE (mat_P_omega(jquad, i_kp)%matrix)
               CALL dbcsr_create(matrix=mat_P_omega(jquad, i_kp)%matrix, &
                                 template=mat_P_global%matrix)
               CALL dbcsr_set(mat_P_omega(jquad, i_kp)%matrix, 0.0_dp)
            END DO
         END DO

         IF (my_open_shell .AND. do_im_time .AND. do_ri_sos_laplace_mp2) THEN
            NULLIFY (mat_P_omega_beta)
            CALL dbcsr_allocate_matrix_set(mat_P_omega_beta, num_integ_points, size_P)
            DO jquad = 1, num_integ_points
               DO i_kp = 1, size_P
                  ALLOCATE (mat_P_omega_beta(jquad, i_kp)%matrix)
                  CALL dbcsr_create(matrix=mat_P_omega_beta(jquad, i_kp)%matrix, &
                                    template=mat_P_global%matrix)
                  CALL dbcsr_set(mat_P_omega_beta(jquad, i_kp)%matrix, 0.0_dp)
               END DO
            END DO
         END IF

         IF (do_kpoints_cubic_RPA) THEN

            NULLIFY (mat_P_omega_im_part)
            CALL dbcsr_allocate_matrix_set(mat_P_omega_im_part, num_integ_points, size_P)
            DO jquad = 1, num_integ_points
               DO i_kp = 1, size_P
                  ALLOCATE (mat_P_omega_im_part(jquad, i_kp)%matrix)
                  CALL dbcsr_create(matrix=mat_P_omega_im_part(jquad, i_kp)%matrix, &
                                    template=mat_P_global%matrix)
                  CALL dbcsr_set(mat_P_omega_im_part(jquad, i_kp)%matrix, 0.0_dp)
               END DO
            END DO

         END IF

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_14", handle4)

         NULLIFY (mat_P_global_copy%matrix)
         ALLOCATE (mat_P_global_copy%matrix)
         CALL dbcsr_create(mat_P_global_copy%matrix, template=mat_P_global%matrix)
         CALL dbcsr_copy(mat_P_global_copy%matrix, mat_P_global%matrix)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_15", handle4)

         n_group_row = mp2_env%ri_rpa_im_time_util(1)%n_group_row
         ALLOCATE (sizes_array_prim_row(0:n_group_row-1, cut_memory))
         DO i_mem = 1, cut_memory
            sizes_array_prim_row(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(:)
         END DO
         ALLOCATE (starts_array_prim_row(0:n_group_row-1, cut_memory))
         DO i_mem = 1, cut_memory
            starts_array_prim_row(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(:)
         END DO
         ALLOCATE (ends_array_prim_row(0:n_group_row-1, cut_memory))
         DO i_mem = 1, cut_memory
            ends_array_prim_row(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(:)
         END DO

         ALLOCATE (starts_array_prim_fullrow(0:n_group_row-1, cut_memory))
         DO i_mem = 1, cut_memory
            starts_array_prim_fullrow(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(:)
         END DO
         ALLOCATE (ends_array_prim_fullrow(0:n_group_row-1, cut_memory))
         DO i_mem = 1, cut_memory
            ends_array_prim_fullrow(:, i_mem) = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(:)
         END DO

         n_group_col = mp2_env%ri_rpa_im_time_util(1)%n_group_col
         ALLOCATE (sizes_array_prim_col(0:n_group_col-1, cut_memory))
         DO j_mem = 1, cut_memory
            sizes_array_prim_col(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(:)
         END DO
         ALLOCATE (starts_array_prim_col(0:n_group_col-1, cut_memory))
         DO j_mem = 1, cut_memory
            starts_array_prim_col(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(:)
         END DO
         ALLOCATE (ends_array_prim_col(0:n_group_col-1, cut_memory))
         DO j_mem = 1, cut_memory
            ends_array_prim_col(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(:)
         END DO

         ALLOCATE (starts_array_prim_fullcol(0:n_group_col-1, cut_memory))
         DO j_mem = 1, cut_memory
            starts_array_prim_fullcol(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(:)
         END DO
         ALLOCATE (ends_array_prim_fullcol(0:n_group_col-1, cut_memory))
         DO j_mem = 1, cut_memory
            ends_array_prim_fullcol(:, j_mem) = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(:)
         END DO

         ALLOCATE (offset_combi_block(cut_memory, cut_memory))

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_16", handle4)

         color_sub_row = mp2_env%ri_rpa_im_time_util(1)%color_sub_row
         color_sub_col = mp2_env%ri_rpa_im_time_util(1)%color_sub_col

         DO i_mem = 1, cut_memory
            DO j_mem = 1, cut_memory

               n_local_row = sizes_array_prim_row(color_sub_row, i_mem)
               row_start_local = starts_array_prim_row(color_sub_row, i_mem)

               n_local_col = sizes_array_prim_col(color_sub_col, j_mem)
               col_start_local = starts_array_prim_col(color_sub_col, j_mem)

               ALLOCATE (offset_combi_block(i_mem, j_mem)%array(row_start_local:row_start_local+n_local_row-1, &
                                                                col_start_local:col_start_local+n_local_col-1))
               offset_combi_block(i_mem, j_mem)%array(:, :) = &
                  mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block(:, :)

            END DO
         END DO

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_17", handle4)

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

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_18", handle4)

         NULLIFY (fm_mat_work)
         CALL cp_fm_create(fm_mat_work, fm_mat_Q%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_Q, fm_mat_work)
         CALL cp_fm_set_all(matrix=fm_mat_work, alpha=0.0_dp)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_19", handle4)

         ALLOCATE (fm_mat_L(SIZE(fm_matrix_L_RI_metric, 1), SIZE(fm_matrix_L_RI_metric, 2)))
         DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
            DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 2)
               CALL cp_fm_create(fm_mat_L(i_size, j_size)%matrix, fm_mat_Q%matrix_struct)
               CALL cp_fm_set_all(fm_mat_L(i_size, j_size)%matrix, 0.0_dp)
            END DO
         END DO

         ALLOCATE (row_from_LLL(dimen_RI))
         row_from_LLL = 0

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_20", handle4)

         CALL dbcsr_get_info(mat_M_P_munu_occ%matrix, &
                             nblkrows_total=nblkrows_total, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_size)

         DO LLL = 1, dimen_RI
            DO row = 1, nblkrows_total
               IF (row_blk_offset(row) <= LLL .AND. LLL < row_blk_offset(row)+row_blk_size(row)) THEN
                  row_from_LLL(LLL) = row
               END IF
            END DO
         END DO

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_20", handle4)

         NULLIFY (fm_mo_coeff_occ_scaled)
         CALL cp_fm_create(fm_mo_coeff_occ_scaled, fm_mo_coeff_occ%matrix_struct)
         CALL cp_fm_to_fm(fm_mo_coeff_occ, fm_mo_coeff_occ_scaled)
         CALL cp_fm_set_all(matrix=fm_mo_coeff_occ_scaled, alpha=0.0_dp)

         NULLIFY (fm_mo_coeff_virt_scaled)
         CALL cp_fm_create(fm_mo_coeff_virt_scaled, fm_mo_coeff_virt%matrix_struct)
         CALL cp_fm_to_fm(fm_mo_coeff_virt, fm_mo_coeff_virt_scaled)
         CALL cp_fm_set_all(matrix=fm_mo_coeff_virt_scaled, alpha=0.0_dp)

         CALL timestop(handle4)

         IF (do_gw_im_time) THEN
            CALL allocate_matrices_gw_im_time(cut_RI, gw_corr_lev_occ, gw_corr_lev_occ_beta, gw_corr_lev_tot, &
                                              gw_corr_lev_virt, gw_corr_lev_virt_beta, homo, homo_beta, nmo, &
                                              num_integ_points, unit_nr, my_group_L_sizes_im_time, my_group_L_starts_im_time, &
                                              row_from_LLL, prim_blk_sizes, RI_blk_sizes, do_ic_model, do_ic_opt_homo_lumo, &
                                              my_open_shell, para_env, para_env_sub, fm_mat_W_tau, fm_mat_Q, &
                                              mo_coeff, mo_coeff_beta, mat_dm_virt_local, mat_P_global, mat_3c_overl_int_gw, &
                                              mat_3c_overl_int_gw_beta, mat_3c_overl_nnP_ic, &
                                              mat_3c_overl_nnP_ic_beta, mat_3c_overl_nnP_ic_reflected, &
                                              mat_3c_overl_nnP_ic_reflected_beta, matrix_s, mat_W, &
                                              mat_3c_overl_int, mat_contr_gf_occ, mat_contr_gf_virt, &
                                              mat_contr_W, mp2_env, qs_env)

         END IF

         CALL timeset(routineN//"_im_t_alloc_mat_21", handle4)

         CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_cut, mat_3c_overl_int, cut_memory, cut_RI, &
                                       starts_array_cm, ends_array_cm, my_group_L_sizes_im_time, eps_filter)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_22", handle4)

         IF (do_mao) THEN

            ! in setup_mat_for_mem_cut_3c, one deallocates mat_3c_overl_int, therefore for the beginning, deallocate
            ! the mao 3c overlap integrals here
            CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_occ, &
                                          cut_memory, cut_RI, starts_array_cm_mao_virt, ends_array_cm_mao_virt, &
                                          my_group_L_sizes_im_time, eps_filter)
            CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_mao_for_virt_cut, mat_3c_overl_int_mao_for_virt, &
                                          cut_memory, cut_RI, starts_array_cm_mao_occ, ends_array_cm_mao_occ, &
                                          my_group_L_sizes_im_time, eps_filter)

         ELSE

            mat_3c_overl_int_mao_for_occ_cut => mat_3c_overl_int_cut
            mat_3c_overl_int_mao_for_virt_cut => mat_3c_overl_int_cut

         END IF

         CALL get_non_zero_blocks_3c_cut_col(mat_3c_overl_int_cut, para_env_sub, cut_RI, &
                                             cut_memory, non_zero_blocks_3c_cut_col)

         CALL check_sparsity_arrays_for_kp(needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                           mat_3c_overl_int_cut, cut_RI, cut_memory, para_env_sub, &
                                           do_kpoints_cubic_RPA)

         ALLOCATE (cycle_due_to_sparse_dm(cut_memory, cut_memory, num_integ_points))
         cycle_due_to_sparse_dm = .FALSE.

         ALLOCATE (multiply_needed_occ(cut_memory, cut_RI, num_cells_dm))
         multiply_needed_occ = .TRUE.

         ALLOCATE (multiply_needed_virt(cut_memory, cut_RI, num_cells_dm))
         multiply_needed_virt = .TRUE.

         ALLOCATE (has_mat_P_blocks(num_cells_dm/2+1, cut_memory, cut_memory, num_3c_repl, num_3c_repl))
         has_mat_P_blocks = .TRUE.

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_23", handle4)

         NULLIFY (fm_mat_L_transposed)
         CALL cp_fm_create(fm_mat_L_transposed, fm_mat_L(1, 1)%matrix%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_L(1, 1)%matrix, fm_mat_L_transposed)
         CALL cp_fm_set_all(matrix=fm_mat_L_transposed, alpha=0.0_dp)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_24", handle4)

         CALL cp_fm_get_info(fm_mat_L_transposed, context=blacs_env)

         ! without kpoints, the size of fm_mat_L is 1x1. with kpoints, the size is N_kpoints x 2 (2 for real/complex)
         DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
         DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 2)

            CALL cp_fm_copy_general(fm_matrix_L_RI_metric(i_size, j_size)%matrix, fm_mat_L_transposed, blacs_env%para_env)
            IF (do_kpoints_cubic_RPA) THEN
               CALL cp_fm_to_fm(fm_mat_L_transposed, fm_mat_L(i_size, j_size)%matrix)
            ELSE
               CALL cp_fm_transpose(fm_mat_L_transposed, fm_mat_L(i_size, j_size)%matrix)
            END IF

         END DO
         END DO

         DO i_size = 1, SIZE(fm_matrix_L_RI_metric, 1)
         DO j_size = 1, SIZE(fm_matrix_L_RI_metric, 2)
            CALL cp_fm_release(fm_matrix_L_RI_metric(i_size, j_size)%matrix)
         END DO
         END DO
         DEALLOCATE (fm_matrix_L_RI_metric)

         CALL cp_fm_release(fm_mat_L_transposed)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_25", handle4)

         NULLIFY (mat_L%matrix)
         ALLOCATE (mat_L%matrix)
         CALL dbcsr_create(mat_L%matrix, template=mat_P_global%matrix)

         CALL timestop(handle4)
         CALL timeset(routineN//"_im_t_alloc_mat_26", handle4)

         IF (.NOT. do_kpoints_cubic_RPA) THEN
            CALL copy_fm_to_dbcsr(fm_mat_L(1, 1)%matrix, mat_L%matrix)
         END IF

         CALL timestop(handle4)

         IF (do_ri_Sigma_x .OR. do_ic_model) THEN

            NULLIFY (mat_SinvVSinv%matrix)
            ALLOCATE (mat_SinvVSinv%matrix)
            CALL dbcsr_create(mat_SinvVSinv%matrix, template=mat_P_global%matrix)

            !  get the Coulomb matrix for Sigma_x = G*V
            CALL dbcsr_multiply("T", "N", 1.0_dp, mat_L%matrix, mat_L%matrix, &
                                0.0_dp, mat_SinvVSinv%matrix, filter_eps=eps_filter_im_time)

         END IF

         IF (do_ri_Sigma_x) THEN

            NULLIFY (mat_dm%matrix)
            ALLOCATE (mat_dm%matrix)
            CALL dbcsr_create(mat_dm%matrix, template=matrix_s(1)%matrix)

         END IF

      END IF

      IF (do_minimax_quad .OR. do_ri_sos_laplace_mp2) THEN
         ! MINIMAX quadrature

         IF (do_kpoints_cubic_RPA) THEN

            CALL gap_and_max_eig_diff_kpoints(qs_env, para_env, Emin, Emax, e_fermi)

            E_Range = Emax/Emin

         ELSE

            Emin = Eigenval(homo+1)-Eigenval(homo)
            Emax = MAXVAL(Eigenval)-MINVAL(Eigenval)
            IF (my_open_shell) THEN
               IF (homo_beta > 0) THEN
                  Emin_beta = Eigenval_beta(homo_beta+1)-Eigenval_beta(homo_beta)
                  Emax_beta = MAXVAL(Eigenval_beta)-MINVAL(Eigenval_beta)
                  Emin = MIN(Emin, Emin_beta)
                  Emax = MAX(Emax, Emax_beta)
               END IF
            END IF
            E_Range = Emax/Emin

         END IF

         IF (.NOT. do_ri_sos_laplace_mp2) THEN
            ALLOCATE (x_tw(2*num_integ_points))
            x_tw = 0.0_dp
            ierr = 0
            CALL get_rpa_minimax_coeff(num_integ_points, E_Range, x_tw, ierr)

            DO jquad = 1, num_integ_points
               tj(jquad) = x_tw(jquad)
               wj(jquad) = x_tw(jquad+num_integ_points)
            END DO

            DEALLOCATE (x_tw)

            IF (unit_nr > 0 .AND. .NOT. do_ic_model) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "INTEG_INFO| Range for the minimax approximation:", E_Range
               WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") "INTEG_INFO| Minimax parameters:", "Weights", "Abscissas"
               DO jquad = 1, num_integ_points
                  WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") wj(jquad), tj(jquad)
               END DO
               CALL m_flush(unit_nr)
            END IF

            ! scale the minimax parameters
            tj(:) = tj(:)*Emin
            wj(:) = wj(:)*Emin
         ELSE
            ! When we perform SOS-MP2, we need an additional factor of 2 for the energies (compare with mp2_laplace.F)
            ! We do not need weights etc. for the cosine transform
            ! We do not scale Emax because it is not needed for SOS-MP2
            Emin = Emin*2.0_dp
         END IF

         ! set up the minimax time grid
         IF (do_im_time) THEN

            ALLOCATE (x_tw(2*num_integ_points))
            x_tw = 0.0_dp

            CALL get_exp_minimax_coeff(num_integ_points, E_Range, x_tw)

            ! For RPA we include already a factor of two (see later steps)
            scaling = 2.0_dp
            IF (do_ri_sos_laplace_mp2) scaling = 1.0_dp

            DO jquad = 1, num_integ_points
               tau_tj(jquad) = x_tw(jquad)/scaling
               tau_wj(jquad) = x_tw(jquad+num_integ_points)/scaling
            END DO

            DEALLOCATE (x_tw)

            IF (unit_nr > 0 .AND. .NOT. do_ic_model) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "INTEG_INFO| Range for the minimax approximation:", E_Range
               ! For testing the gap
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "INTEG_INFO| Gap:", Emin
               WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
                  "INTEG_INFO| Minimax parameters of the time grid:", "Weights", "Abscissas"
               DO jquad = 1, num_integ_points
                  WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") tau_wj(jquad), tau_tj(jquad)
               END DO
               CALL m_flush(unit_nr)
            END IF

            ! scale grid from [1,R] to [Emin,Emax]
            tau_tj(:) = tau_tj(:)/Emin
            tau_wj(:) = tau_wj(:)/Emin

            IF (.NOT. do_ri_sos_laplace_mp2) THEN
               num_points_per_magnitude = mp2_env%ri_rpa_im_time%num_points_per_magnitude
               CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                                 Emin, Emax, max_error_min, num_points_per_magnitude)

               IF (do_gw_im_time) THEN
                  CALL get_weights_gw(num_integ_points, num_points_per_magnitude, unit_nr, Emin, Emax, &
                                      max_error_min, tau_tj, tj, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)
               END IF

            END IF

         END IF

      ELSE
         ! Clenshaw-Curtius quadrature
         DO jquad = 1, num_integ_points-1
            tj(jquad) = jquad*pi/(2.0_dp*num_integ_points)
            wj(jquad) = pi/(num_integ_points*SIN(tj(jquad))**2)
         END DO
         tj(num_integ_points) = pi/2.0_dp
         wj(num_integ_points) = pi/(2.0_dp*num_integ_points*SIN(tj(num_integ_points))**2)

         a_scaling = 1.0_dp
         IF (my_open_shell) THEN
            CALL calc_scaling_factor(a_scaling, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                     num_integ_points, num_integ_group, color_rpa_group, &
                                     tj, wj, fm_mat_S, &
                                     homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta)
         ELSE
            CALL calc_scaling_factor(a_scaling, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                     num_integ_points, num_integ_group, color_rpa_group, &
                                     tj, wj, fm_mat_S)
         END IF

         ! for G0W0, we may set the scaling factor by hand
         IF (my_do_gw) THEN
            a_scaling = 0.2_dp
         END IF

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.5)') 'INTEG_INFO| Scaling parameter:', a_scaling

         wj(:) = wj(:)*a_scaling

      END IF

      IF (.NOT. do_im_time) THEN

         ! initialize buffer for matrix redistribution
         CALL initialize_buffer(fm_mat_Q_gemm, fm_mat_Q, RPA_proc_map, buffer_rec, buffer_send, &
                                number_of_rec, number_of_send, &
                                map_send_size, map_rec_size, local_size_source, para_env_RPA)
         IF (my_open_shell) THEN
            CALL initialize_buffer(fm_mat_Q_gemm_beta, fm_mat_Q_beta, RPA_proc_map, buffer_rec_beta, buffer_send_beta, &
                                   number_of_rec_beta, number_of_send_beta, &
                                   map_send_size_beta, map_rec_size_beta, local_size_source_beta, para_env_RPA)
         END IF
         ! Another buffer for AXK RPA
         IF (mp2_env%ri_rpa%do_ri_axk) THEN
            CALL initialize_buffer(fm_mat_Q, fm_mat_Q_gemm, RPA_proc_map, buffer_rec_axk, buffer_send_axk, &
                                   number_of_rec_axk, number_of_send_axk, &
                                   map_send_size_axk, map_rec_size_axk, local_size_source_axk, para_env_RPA)
         ENDIF

      END IF

      IF (my_open_shell .OR. do_ri_sos_laplace_mp2) THEN
         alpha = 2.0_dp
      ELSE
         alpha = 4.0_dp
      END IF

      IF (my_do_gw .OR. do_gw_im_time) THEN

         CALL allocate_matrices_gw(vec_Sigma_c_gw, vec_Sigma_c_gw_beta, color_rpa_group, dimen_nm_gw, &
                                   gw_corr_lev_occ, gw_corr_lev_occ_beta, gw_corr_lev_virt, homo, homo_beta, &
                                   nmo, num_integ_points, num_integ_group, unit_nr, &
                                   gw_corr_lev_tot, num_fit_points, &
                                   do_minimax_quad, do_periodic, do_ri_Sigma_x, my_do_gw, my_open_shell, &
                                   first_cycle_periodic_correction, do_GW_corr, &
                                   a_scaling, Eigenval, Eigenval_beta, &
                                   tj, vec_omega_fit_gw, vec_Sigma_x_gw, &
                                   vec_Sigma_x_gw_beta, delta_corr, Eigenval_last, Eigenval_last_beta, &
                                   Eigenval_scf, Eigenval_scf_beta, m_value, m_value_beta, &
                                   vec_gw_energ, vec_gw_energ_beta, vec_gw_energ_error_fit, vec_gw_energ_error_fit_beta, &
                                   vec_W_gw, vec_W_gw_beta, z_value, z_value_beta, fm_mat_S, &
                                   fm_mat_S_gw, fm_mat_S_gw_work, fm_mat_S_gw_work_beta, &
                                   fm_mat_S_beta, fm_mat_S_gw_beta, para_env, mp2_env)
         IF (do_bse) THEN

            num_Z_vectors = mp2_env%ri_g0w0%num_z_vectors
            eps_min_trans = mp2_env%ri_g0w0%eps_min_trans
            max_iter_bse = mp2_env%ri_g0w0%max_iter_bse

            CALL cp_fm_create(fm_mat_Q_static_bse_gemm, fm_mat_Q_gemm%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q_gemm, fm_mat_Q_static_bse_gemm)
            CALL cp_fm_set_all(fm_mat_Q_static_bse_gemm, 0.0_dp)

            CALL cp_fm_create(fm_mat_Q_static_bse, fm_mat_Q%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q, fm_mat_Q_static_bse)
            CALL cp_fm_set_all(fm_mat_Q_static_bse, 0.0_dp)

         END IF

      END IF

      Erpa = 0.0_dp
      IF (mp2_env%ri_rpa%do_ri_axk) e_axk = 0.0_dp
      first_cycle = .TRUE.
      omega_old = 0.0_dp
      my_num_dgemm_call = 0
      my_flop_rate = 0.0_dp

      IF (my_do_gw .OR. do_gw_im_time) THEN
         iter_ev_sc = mp2_env%ri_g0w0%iter_ev_sc
      ELSE
         iter_ev_sc = 1
      END IF

!      IF (do_gw_im_time) THEN
!         ! for cubic-scaling GW, we only do G0W0 and no evGW
!         CPASSERT(iter_ev_sc == 1)
!      END IF

      DO count_ev_sc_GW = 1, iter_ev_sc

         IF (do_ic_model) CYCLE

         ! reset some values, important when doing eigenvalue self-consistent GW
         IF (my_do_gw .OR. do_gw_im_time) THEN
            Erpa = 0.0_dp
            vec_Sigma_c_gw = (0.0_dp, 0.0_dp)
            vec_gw_energ = 0.0_dp
            vec_gw_energ_error_fit = 0.0_dp
            z_value = 0.0_dp
            m_value = 0.0_dp
            first_cycle = .TRUE.
            IF (my_open_shell) THEN
               vec_Sigma_c_gw_beta = (0.0_dp, 0.0_dp)
               vec_gw_energ_beta = 0.0_dp
               vec_gw_energ_error_fit_beta = 0.0_dp
               z_value_beta = 0.0_dp
               m_value_beta = 0.0_dp
            END IF
         END IF

         ! calculate Q_PQ(it)
         IF (do_im_time) THEN

            IF (.NOT. do_kpoints_cubic_RPA) THEN
               e_fermi = (Eigenval(homo)+Eigenval(homo+1))*0.5_dp
               IF (my_open_shell) THEN
                  e_fermi_beta = (Eigenval_beta(homo_beta)+Eigenval_beta(homo_beta+1))*0.5_dp
               END IF
            END IF

            tau = 0.0_dp
            tau_old = 0.0_dp

            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,i15)") &
               "MEMORY_INFO| Memory cut:", cut_memory
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,i15)") &
               "MEMORY_INFO| Im. time group size for RI functions:", mp2_env%ri_rpa_im_time%group_size_3c
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,i15)") &
               "MEMORY_INFO| Im. time group size for local P matrix:", mp2_env%ri_rpa_im_time%group_size_P
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
               "SPARSITY_INFO| Eps pgf orb for imaginary time:", mp2_env%mp2_gpw%eps_grid
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
               "SPARSITY_INFO| Eps filter for imaginary time:", eps_filter
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
               "SPARSITY_INFO| Second eps filter for imaginary time:", eps_filter_im_time

            ! for evGW, we have to ensure that mat_P_omega is zero
            CALL zero_mat_P_omega(mat_P_omega, num_integ_points, size_P)

            ! compute the matrix Q(it) and Fourier transform it directly to mat_P_omega(iw)
            CALL compute_mat_P_omega(mat_P_omega, mat_P_omega_im_part, fm_scaled_dm_occ_tau, &
                                     fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                     fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                     mat_P_local, mat_P_global, mat_P_global_copy, &
                                     mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, matrix_s, &
                                     mao_coeff_occ, mao_coeff_virt, 1, &
                                     mat_M_P_munu_occ, mat_M_P_munu_virt, mat_3c_overl_int_cut, &
                                     mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut, &
                                     mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, &
                                     weights_cos_tf_t_to_w, tj, tau_tj, e_fermi, eps_filter, alpha, &
                                     eps_filter_im_time, Eigenval, nmo, n_group_col, &
                                     group_size_P, num_integ_points, jquad, cut_memory, cut_RI, &
                                     unit_nr, mp2_env, para_env, para_env_sub, &
                                     starts_array_prim_col, ends_array_prim_col, &
                                     starts_array_prim_row, ends_array_prim_row, &
                                     starts_array_prim_fullcol, ends_array_prim_fullcol, &
                                     starts_array_prim_fullrow, ends_array_prim_fullrow, &
                                     my_group_L_starts_im_time, my_group_L_sizes_im_time, &
                                     offset_combi_block, starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                     ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                     mepos_P_from_RI_row, row_from_LLL, cycle_due_to_sparse_dm, &
                                     multiply_needed_occ, multiply_needed_virt, &
                                     non_zero_blocks_3c, non_zero_blocks_3c_cut_col, buffer_mat_M, &
                                     do_mao, stabilize_exp, qs_env, index_to_cell_3c, cell_to_index_3c, &
                                     needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                     has_mat_P_blocks, do_ri_sos_laplace_mp2)

            ! the same for open shell, use fm_mo_coeff_occ_beta and fm_mo_coeff_virt_beta
            IF (my_open_shell) THEN
               ! For SOS-MP2 we need the same calculation for alpha and beta spin independently, for RPA the sum of both
               IF (do_ri_sos_laplace_mp2) THEN
                  CALL zero_mat_P_omega(mat_P_omega_beta, num_integ_points, size_P)

                  CALL compute_mat_P_omega(mat_P_omega_beta, mat_P_omega_im_part, fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           mat_P_local, mat_P_global, mat_P_global_copy, &
                                           mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, matrix_s, &
                                           mao_coeff_occ, mao_coeff_virt, 2, &
                                           mat_M_P_munu_occ, mat_M_P_munu_virt, mat_3c_overl_int_cut, &
                                           mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut, &
                                           mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, &
                                           weights_cos_tf_t_to_w, tj, tau_tj, e_fermi_beta, eps_filter, alpha, &
                                           eps_filter_im_time, Eigenval_beta, nmo, n_group_col, &
                                           group_size_P, num_integ_points, jquad, cut_memory, cut_RI, &
                                           unit_nr, mp2_env, para_env, para_env_sub, &
                                           starts_array_prim_col, ends_array_prim_col, &
                                           starts_array_prim_row, ends_array_prim_row, &
                                           starts_array_prim_fullcol, ends_array_prim_fullcol, &
                                           starts_array_prim_fullrow, ends_array_prim_fullrow, &
                                           my_group_L_starts_im_time, my_group_L_sizes_im_time, &
                                           offset_combi_block, starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                           ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                           mepos_P_from_RI_row, row_from_LLL, cycle_due_to_sparse_dm, &
                                           multiply_needed_occ, multiply_needed_virt, &
                                           non_zero_blocks_3c, non_zero_blocks_3c_cut_col, buffer_mat_M, &
                                           do_mao, stabilize_exp, qs_env, index_to_cell_3c, cell_to_index_3c, &
                                           needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                           has_mat_P_blocks, do_ri_sos_laplace_mp2)
               ELSE
                  CALL compute_mat_P_omega(mat_P_omega, mat_P_omega_im_part, fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           mat_P_local, mat_P_global, mat_P_global_copy, &
                                           mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, matrix_s, &
                                           mao_coeff_occ, mao_coeff_virt, 2, &
                                           mat_M_P_munu_occ, mat_M_P_munu_virt, mat_3c_overl_int_cut, &
                                           mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut, &
                                           mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, &
                                           weights_cos_tf_t_to_w, tj, tau_tj, e_fermi_beta, eps_filter, alpha, &
                                           eps_filter_im_time, Eigenval_beta, nmo, n_group_col, &
                                           group_size_P, num_integ_points, jquad, cut_memory, cut_RI, &
                                           unit_nr, mp2_env, para_env, para_env_sub, &
                                           starts_array_prim_col, ends_array_prim_col, &
                                           starts_array_prim_row, ends_array_prim_row, &
                                           starts_array_prim_fullcol, ends_array_prim_fullcol, &
                                           starts_array_prim_fullrow, ends_array_prim_fullrow, &
                                           my_group_L_starts_im_time, my_group_L_sizes_im_time, &
                                           offset_combi_block, starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                           ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                           mepos_P_from_RI_row, row_from_LLL, cycle_due_to_sparse_dm, &
                                           multiply_needed_occ, multiply_needed_virt, &
                                           non_zero_blocks_3c, non_zero_blocks_3c_cut_col, buffer_mat_M, &
                                           do_mao, stabilize_exp, qs_env, index_to_cell_3c, cell_to_index_3c, &
                                           needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                           has_mat_P_blocks, do_ri_sos_laplace_mp2)
               END IF ! do_ri_sos_laplace_mp2

            END IF ! my_open_shell

         END IF ! do im time

         DO jquad = 1, num_integ_points

            IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE

            CALL timeset(routineN//"_RPA_matrix_operations", handle3)

            IF (.NOT. do_ri_sos_laplace_mp2) THEN

               IF (do_minimax_quad) THEN
                  omega = tj(jquad)
               ELSE
                  omega = a_scaling/TAN(tj(jquad))
               END IF

            END IF ! do_ri_sos_laplace_mp2

            IF (do_im_time) THEN
               ! in case we do imag time, we already calculated fm_mat_Q by a Fourier transform from im. time

               IF (.NOT. do_kpoints_cubic_RPA) THEN

                  ! multiplication with RI metric/Coulomb operator
                  CALL dbcsr_multiply("N", "T", 1.0_dp, mat_P_omega(jquad, 1)%matrix, mat_L%matrix, &
                                      0.0_dp, mat_P_global%matrix, filter_eps=eps_filter_im_time)

                  CALL copy_dbcsr_to_fm(mat_P_global%matrix, fm_mat_work)

                  CALL cp_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_mat_L(1, 1)%matrix, fm_mat_work, &
                               0.0_dp, fm_mat_Q)

                  ! For open-shell SOS-MP2 we have two different matrices we have to deal with
                  IF (my_open_shell .AND. do_ri_sos_laplace_mp2) THEN
                     CALL cp_fm_set_all(matrix=fm_mat_work, alpha=0.0_dp)

                     ! multiplication with RI metric/Coulomb operator
                     CALL dbcsr_multiply("N", "T", 1.0_dp, mat_P_omega_beta(jquad, 1)%matrix, mat_L%matrix, &
                                         0.0_dp, mat_P_global%matrix, filter_eps=eps_filter_im_time)

                     CALL copy_dbcsr_to_fm(mat_P_global%matrix, fm_mat_work)

                     CALL cp_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_mat_L(1, 1)%matrix, fm_mat_work, &
                                  0.0_dp, fm_mat_Q_beta)
                  END IF
               END IF

            ELSE

               CALL matrix_operations_canonical(fm_mat_S, first_cycle, count_ev_sc_GW, virtual, Eigenval, Eigenval_last, homo, &
                                                omega, omega_old, jquad, num_integ_points, &
                                                mm_style, dimen_RI, dimen_ia, alpha, fm_mat_Q, fm_mat_Q_gemm, para_env_RPA, &
                                                do_bse, fm_mat_Q_static_bse_gemm, RPA_proc_map, buffer_rec, buffer_send, &
                                                number_of_send, map_send_size, map_rec_size, local_size_source, &
                                                my_num_dgemm_call, my_flop_rate)

               IF (my_open_shell) THEN
                  CALL matrix_operations_canonical(fm_mat_S_beta, first_cycle, count_ev_sc_GW, virtual_beta, Eigenval_beta, &
                                                   Eigenval_last_beta, homo_beta, omega, omega_old, jquad, num_integ_points, &
                                                   mm_style, dimen_RI, dimen_ia_beta, alpha, fm_mat_Q_beta, fm_mat_Q_gemm_beta, &
                                                   para_env_RPA, do_bse, fm_mat_Q_static_bse_gemm, RPA_proc_map, buffer_rec_beta, &
                                                   buffer_send_beta, number_of_send_beta, map_send_size_beta, map_rec_size_beta, &
                                                   local_size_source_beta, my_num_dgemm_call, my_flop_rate)

                  CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_Q, beta=1.0_dp, matrix_b=fm_mat_Q_beta)

               END IF ! open shell

            END IF ! im time

            ! Calculate AXK energy correction
            IF (mp2_env%ri_rpa%do_ri_axk) THEN
               CALL compute_axk_ener(qs_env, fm_mat_Q, fm_mat_Q_gemm, dimen_RI, dimen_ia, &
                                     para_env_sub, &
                                     para_env_RPA, RPA_proc_map, eigenval, fm_mat_S, homo, virtual, omega, &
                                     buffer_send_axk, buffer_rec_axk, &
                                     number_of_send_axk, map_send_size_axk, map_rec_size_axk, &
                                     local_size_source_axk, mp2_env, mat_munu, unit_nr, e_axk_corr)

               ! Evaluate the final AXK energy correction
               e_axk = e_axk+e_axk_corr*wj(jquad)
            ENDIF ! do_ri_axk

            ! get info of fm_mat_Q
            CALL cp_fm_get_info(matrix=fm_mat_Q, &
                                nrow_local=nrow_local, &
                                ncol_local=ncol_local, &
                                row_indices=row_indices, &
                                col_indices=col_indices)

            IF (do_ri_sos_laplace_mp2) THEN
               ! calculate the trace of the product Q*Q
               trace_XX = 0.0_dp
               IF (my_open_shell) THEN
                  DO jjB = 1, ncol_local
                     trace_XX = trace_XX+DOT_PRODUCT(fm_mat_Q%local_data(:, jjB), fm_mat_Q_beta%local_data(:, jjB))
                  END DO
               ELSE
                  trace_XX = NORM2(fm_mat_Q%local_data)
                  trace_XX = trace_XX*trace_XX
               END IF

               Erpa = Erpa-trace_XX*tau_wj(jquad)

            ELSE
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
                        fm_mat_Q%local_data(iiB, jjB) = fm_mat_Q%local_data(iiB, jjB)+1.0_dp
                     END IF
                  END DO
               END DO
               CALL mp_sum(trace_Qomega, para_env_RPA%group)

               IF (.NOT. do_kpoints_cubic_RPA) THEN

                  ! calculate Trace(Log(Matrix)) as Log(DET(Matrix)) via cholesky decomposition
                  CALL cp_fm_cholesky_decompose(matrix=fm_mat_Q, n=dimen_RI, info_out=info_chol)
                  CPASSERT(info_chol == 0)

                  ! get info of cholesky_decomposed(fm_mat_Q)
                  CALL cp_fm_get_info(matrix=fm_mat_Q, &
                                      nrow_local=nrow_local, &
                                      ncol_local=ncol_local, &
                                      row_indices=row_indices, &
                                      col_indices=col_indices)

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
                     ! FComega=FComega+(LOG(Q_log(iiB))-trace_Qomega(iiB))/2.0_dp
                     FComega = FComega+(Q_log(iiB)-trace_Qomega(iiB))/2.0_dp
                  END DO

                  Erpa = Erpa+FComega*wj(jquad)

               END IF

               IF (do_kpoints_cubic_RPA) THEN

                  CALL get_kpoint_info(kpoints, nkp=nkp)

                  DO ikp = 1, nkp

                     ! 1. multiplication Q(iw,k) = K^T(k)P(iw,k)K(k)
                     CALL compute_Q(cfm_mat_Q, &
                                    mat_P_omega(jquad, ikp)%matrix, &
                                    mat_P_omega_im_part(jquad, ikp)%matrix, &
                                    fm_mat_L(ikp, 1)%matrix, &
                                    fm_mat_L(ikp, 2)%matrix, &
                                    dimen_RI)

                     ! 2. Cholesky decomposition of Id + Q(iw,k)
                     CALL cholesky_decomp_Q(cfm_mat_Q, para_env_RPA, trace_Qomega, dimen_RI)

                     ! 3. Computing E_c^RPA = E_c^RPA + a_w/N_k*sum_k ln[det(1+Q(iw,k))-Tr(Q(iw,k))]
                     CALL frequency_and_kpoint_integration(Erpa, cfm_mat_Q, para_env_RPA, trace_Qomega, &
                                                           dimen_RI, wj(jquad), kpoints%wkp(ikp))

                  END DO

               END IF
            END IF ! do_ri_sos_laplace_mp2

            ! save omega and reset the first_cycle flag
            first_cycle = .FALSE.
            omega_old = omega

            CALL timestop(handle3)

            ! the actual G0W0 calculation
            IF (my_do_gw .OR. do_gw_im_time) THEN
               CALL GW_matrix_operations(vec_Sigma_c_gw, vec_Sigma_c_gw_beta, &
                                         dimen_nm_gw, dimen_RI, gw_corr_lev_occ, &
                                         gw_corr_lev_occ_beta, gw_corr_lev_virt, homo, homo_beta, jquad, &
                                         nmo, num_fit_points, num_integ_points, ncol_local, nrow_local, &
                                         col_indices, row_indices, do_bse, &
                                         do_gw_im_time, do_periodic, my_do_gw, my_open_shell, &
                                         first_cycle_periodic_correction, fermi_level_offset, &
                                         omega, Eigenval, Eigenval_beta, delta_corr, tau_tj, tj, vec_omega_fit_gw, &
                                         vec_W_gw, vec_W_gw_beta, wj, weights_cos_tf_w_to_t, fm_mat_W_tau, fm_mat_L, &
                                         fm_mat_Q, fm_mat_Q_static_bse, fm_mat_R_gw, fm_mat_S_gw, fm_mat_S_gw_beta, &
                                         fm_mat_S_gw_work, fm_mat_S_gw_work_beta, fm_mat_work, mo_coeff, para_env, &
                                         para_env_RPA, matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, &
                                         kpoints, qs_env, mp2_env)
            END IF

         END DO ! jquad

         CALL mp_sum(Erpa, para_env%group)

         IF (.NOT. do_ri_sos_laplace_mp2) THEN
            Erpa = Erpa/(pi*2.0_dp)
            IF (do_minimax_quad) Erpa = Erpa/2.0_dp
         END IF

         IF (mp2_env%ri_rpa%do_ri_axk) THEN
            CALL mp_sum(E_axk, para_env%group)
            E_axk = E_axk/(pi*2.0_dp)
            IF (do_minimax_quad) E_axk = E_axk/2.0_dp
            mp2_env%ri_rpa%ener_axk = E_axk
         ENDIF

         IF (.NOT. do_im_time) THEN

            IF (para_env_RPA%mepos == 0) my_flop_rate = my_flop_rate/REAL(MAX(my_num_dgemm_call, 1), KIND=dp)/1.0E9_dp
            CALL mp_sum(my_flop_rate, para_env%group)
            my_flop_rate = my_flop_rate/para_env%num_pe
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.2)") &
               "PERFORMANCE| PDGEMM flop rate (Gflops / MPI rank):", my_flop_rate
         END IF

         ! G0W0 postprocessing: Fitting + correction of MO energies
         IF (my_do_gw .OR. do_gw_im_time) THEN
            CALL gw_postprocessing(vec_Sigma_c_gw, vec_Sigma_c_gw_beta, count_ev_sc_GW, &
                                   gw_corr_lev_occ, gw_corr_lev_occ_beta, &
                                   gw_corr_lev_tot, gw_corr_lev_virt, gw_corr_lev_virt_beta, &
                                   homo, homo_beta, &
                                   nmo, num_fit_points, num_integ_points, &
                                   num_points_corr, unit_nr, &
                                   do_apply_ic_corr_to_gw, do_gw_im_time, &
                                   do_periodic, do_ri_Sigma_x, &
                                   first_cycle_periodic_correction, my_do_gw, my_open_shell, &
                                   do_GW_corr, e_fermi, &
                                   e_fermi_beta, eps_filter, &
                                   fermi_level_offset, stabilize_exp, delta_corr, Eigenval, &
                                   Eigenval_beta, Eigenval_last, Eigenval_last_beta, Eigenval_scf, &
                                   Eigenval_scf_beta, m_value, m_value_beta, &
                                   tau_tj, tj, vec_gw_energ, vec_gw_energ_beta, vec_gw_energ_error_fit, &
                                   vec_gw_energ_error_fit_beta, vec_omega_fit_gw, &
                                   vec_Sigma_x_gw, vec_Sigma_x_gw_beta, &
                                   z_value, z_value_beta, ic_corr_list, ic_corr_list_beta, &
                                   weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, &
                                   fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, fm_mo_coeff_occ, &
                                   fm_mo_coeff_virt, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                                   fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, mo_coeff, fm_mat_W_tau, &
                                   para_env, para_env_RPA, mat_dm, mat_SinvVSinv, &
                                   mat_3c_overl_int_gw, mat_3c_overl_int_gw_beta, matrix_berry_im_mo_mo, &
                                   matrix_berry_re_mo_mo, mat_greens_fct_occ, mat_greens_fct_occ_beta, &
                                   mat_greens_fct_virt, mat_greens_fct_virt_beta, mat_W, matrix_s, &
                                   mat_contr_gf_occ, mat_contr_gf_virt, mat_contr_W, kpoints, mp2_env, qs_env)

            ! if HOMO-LUMO gap differs by less then mp2_env%ri_g0w0%eps_ev_sc_iter, exit ev sc GW loop
            IF (ABS(Eigenval(homo)-Eigenval_last(homo)-Eigenval(homo+1)+Eigenval_last(homo+1)) &
                < mp2_env%ri_g0w0%eps_ev_sc_iter) THEN
               EXIT
            END IF

         END IF ! my_do_gw if

      END DO !ev_sc_gw_loop

      IF (do_ic_model) THEN

         IF (my_open_shell) THEN

            CALL calculate_ic_correction(Eigenval, mat_SinvVSinv, mat_3c_overl_nnP_ic, &
                                         mat_3c_overl_nnP_ic_reflected, mat_contr_gf_occ, matrix_s, gw_corr_lev_tot, &
                                         gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, dimen_RI, unit_nr, &
                                         do_ic_opt_homo_lumo, fm_mat_Q, para_env, mp2_env, do_alpha=.TRUE.)

            CALL calculate_ic_correction(Eigenval_beta, mat_SinvVSinv, mat_3c_overl_nnP_ic_beta, &
                                         mat_3c_overl_nnP_ic_reflected_beta, mat_contr_gf_occ, matrix_s, gw_corr_lev_tot, &
                                         gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, homo_beta, nmo, dimen_RI, unit_nr, &
                                         do_ic_opt_homo_lumo, fm_mat_Q, para_env, mp2_env, do_beta=.TRUE.)

         ELSE

            CALL calculate_ic_correction(Eigenval, mat_SinvVSinv, mat_3c_overl_nnP_ic, &
                                         mat_3c_overl_nnP_ic_reflected, mat_contr_gf_occ, matrix_s, gw_corr_lev_tot, &
                                         gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, dimen_RI, unit_nr, &
                                         do_ic_opt_homo_lumo, fm_mat_Q, para_env, mp2_env)

         END IF

      END IF

      ! postprocessing after GW for Bethe-Salpeter
      IF (do_bse) THEN
         CALL mult_B_with_W_and_fill_local_3c_arrays(fm_mat_S_ij_bse, fm_mat_S_ab_bse, fm_mat_S, fm_mat_Q_static_bse, &
                                                     fm_mat_Q_static_bse_gemm, &
                                                     B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                     B_iaQ_bse_local, dimen_RI, homo, virtual, dimen_ia, &
                                                     starts_array, ends_array, sizes_array, color_sub, para_env)

         CALL do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo, virtual, num_Z_vectors, &
                                     max_iter_bse, eps_min_trans, Eigenval, para_env)

      END IF

      ! release buffer
      CALL release_buffer(RPA_proc_map, buffer_rec, buffer_send, &
                          number_of_rec, number_of_send, &
                          map_send_size, map_rec_size, local_size_source)

      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         CALL release_buffer(RPA_proc_map, buffer_rec_axk, buffer_send_axk, &
                             number_of_rec_axk, number_of_send_axk, &
                             map_send_size_axk, map_rec_size_axk, local_size_source_axk)
      ENDIF

      IF (my_open_shell) THEN
         CALL release_buffer(RPA_proc_map, buffer_rec_beta, buffer_send_beta, &
                             number_of_rec_beta, number_of_send_beta, &
                             map_send_size_beta, map_rec_size_beta, local_size_source_beta)
      END IF

      IF (my_do_gw .OR. do_gw_im_time) THEN
         CALL gw_cleanup(vec_Sigma_c_gw, vec_Sigma_c_gw_beta, num_integ_points, &
                         do_gw_im_time, do_ic_model, do_im_time, do_periodic, do_ri_Sigma_x, &
                         my_do_gw, my_open_shell, do_GW_corr, vec_omega_fit_gw, vec_Sigma_x_gw, &
                         vec_Sigma_x_gw_beta, Eigenval_last, Eigenval_last_beta, Eigenval_scf, &
                         Eigenval_scf_beta, m_value, m_value_beta, vec_gw_energ, &
                         vec_gw_energ_beta, vec_gw_energ_error_fit, vec_gw_energ_error_fit_beta, &
                         vec_W_gw, vec_W_gw_beta, z_value, z_value_beta, weights_cos_tf_w_to_t, &
                         weights_sin_tf_t_to_w, fm_mat_W_tau, fm_mat_S_gw_work, fm_mat_S_gw_work_beta, &
                         mat_3c_overl_int_gw, mat_3c_overl_int_gw_beta, mat_3c_overl_nnP_ic, &
                         mat_3c_overl_nnP_ic_beta, mat_3c_overl_nnP_ic_reflected, &
                         mat_3c_overl_nnP_ic_reflected_beta, mat_greens_fct_occ, mat_greens_fct_occ_beta, &
                         mat_greens_fct_virt, mat_greens_fct_virt_beta, matrix_berry_im_mo_mo, &
                         matrix_berry_re_mo_mo, mat_W, mat_contr_gf_occ, mat_contr_gf_virt, &
                         mat_contr_W, kpoints, mp2_env)
      END IF

      IF (do_im_time) THEN
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
               CALL cp_fm_release(fm_mat_L(i_size, j_size)%matrix)
            END DO
         END DO
         DEALLOCATE (fm_mat_L)
         CALL cp_fm_release(fm_mat_work)

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

         CALL dbcsr_deallocate_matrix_set(mat_M_mu_Pnu_occ)
         CALL dbcsr_deallocate_matrix_set(mat_M_mu_Pnu_virt)

         CALL dbcsr_deallocate_matrix_set(mat_P_omega)
         IF (my_open_shell .AND. do_im_time .AND. do_ri_sos_laplace_mp2) CALL dbcsr_deallocate_matrix_set(mat_P_omega_beta)

         CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_cut)
         IF (do_mao) THEN
            CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_mao_for_occ_cut)
            CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int_mao_for_virt_cut)
         END IF

         DEALLOCATE (buffer_mat_M)
         DEALLOCATE (non_zero_blocks_3c)
         DEALLOCATE (non_zero_blocks_3c_cut_col)
         DEALLOCATE (cycle_due_to_sparse_dm, multiply_needed_occ, multiply_needed_virt)
         DEALLOCATE (row_from_LLL)
         DEALLOCATE (needed_cutRI_mem_R1vec_R2vec_for_kp)
         DEALLOCATE (has_mat_P_blocks)

         IF (do_kpoints_cubic_RPA) THEN
            CALL dbcsr_deallocate_matrix_set(mat_P_omega_im_part)
            CALL cp_cfm_release(cfm_mat_Q)
         END IF

      END IF

      DEALLOCATE (tj)
      DEALLOCATE (wj)
      DEALLOCATE (Q_log)
      DEALLOCATE (trace_Qomega)

      CALL timestop(handle)

   END SUBROUTINE rpa_num_int

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
!> \param jquad ...
!> \param num_integ_points ...
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
   SUBROUTINE matrix_operations_canonical(fm_mat_S, first_cycle, count_ev_sc_GW, virtual, Eigenval, Eigenval_last, homo, &
                                          omega, omega_old, jquad, num_integ_points, &
                                          mm_style, dimen_RI, dimen_ia, alpha, fm_mat_Q, fm_mat_Q_gemm, para_env_RPA, &
                                          do_bse, fm_mat_Q_static_bse_gemm, RPA_proc_map, buffer_rec, buffer_send, &
                                    number_of_send, map_send_size, map_rec_size, local_size_source, my_num_dgemm_call, my_flop_rate)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: first_cycle
      INTEGER, INTENT(IN)                                :: count_ev_sc_GW, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval, Eigenval_last
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), INTENT(IN)                          :: omega, omega_old
      INTEGER, INTENT(IN)                                :: jquad, num_integ_points, mm_style, &
                                                            dimen_RI, dimen_ia
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, fm_mat_Q_gemm
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      LOGICAL, INTENT(IN)                                :: do_bse
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q_static_bse_gemm
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: RPA_proc_map
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: buffer_rec, buffer_send
      INTEGER, INTENT(IN)                                :: number_of_send
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: map_send_size, map_rec_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: local_size_source
      INTEGER, INTENT(INOUT)                             :: my_num_dgemm_call
      REAL(KIND=dp), INTENT(INOUT)                       :: my_flop_rate

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_operations_canonical', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: avirt, handle, i_global, iiB, iocc, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: actual_flop_rate, eigen_diff, t_end, &
                                                            t_start

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

               iocc = MAX(1, i_global-1)/virtual+1
               avirt = i_global-(iocc-1)*virtual
               eigen_diff = Eigenval_last(avirt+homo)-Eigenval_last(iocc)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)/ &
                                               SQRT(eigen_diff/(eigen_diff**2+omega_old**2))

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

               iocc = MAX(1, i_global-1)/virtual+1
               avirt = i_global-(iocc-1)*virtual
               eigen_diff = Eigenval(avirt+homo)-Eigenval(iocc)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)* &
                                               SQRT(eigen_diff/(eigen_diff**2+omega**2))

            END DO
         END DO
      ELSE
         ! In this case the update has to remove the old omega component thus
         ! SQRT(((epsi_a-epsi_i)**2+omega_old**2)/((epsi_a-epsi_i)**2+omega**2))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,&
!$OMP                    fm_mat_S,virtual,homo,omega,omega_old,jquad,num_integ_points)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)

               iocc = MAX(1, i_global-1)/virtual+1
               avirt = i_global-(iocc-1)*virtual
               eigen_diff = Eigenval(avirt+homo)-Eigenval(iocc)

               fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)* &
                                               SQRT((eigen_diff**2+omega_old**2)/(eigen_diff**2+omega**2))

            END DO
         END DO
      END IF

      ! alpha=4.0 is valid only for closed shell systems (alpha is 2.0 in the spin-orbital basis)
      t_start = m_walltime()
      SELECT CASE (mm_style)
      CASE (wfc_mm_style_gemm)
         ! waste-fully computes the full symmetrix matrix, but maybe faster than cp_fm_syrk for optimized cp_fm_gemm
         ! RPA for imaginary frequencies
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

      actual_flop_rate = 2.0_dp*REAL(dimen_ia, KIND=dp)*dimen_RI*REAL(dimen_RI, KIND=dp)
      actual_flop_rate = actual_flop_rate/(MAX(0.01_dp, t_end-t_start))
      IF (para_env_RPA%mepos == 0) my_flop_rate = my_flop_rate+actual_flop_rate
      my_num_dgemm_call = my_num_dgemm_call+1

      IF (do_bse .AND. jquad == 1) THEN
         CALL cp_fm_to_fm(fm_mat_Q_gemm, fm_mat_Q_static_bse_gemm)
      END IF

      ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
      CALL cp_fm_set_all(matrix=fm_mat_Q, alpha=0.0_dp)
      CALL fm_redistribute(fm_mat_Q_gemm, fm_mat_Q, RPA_proc_map, buffer_rec, buffer_send, &
                           number_of_send, &
                           map_send_size, map_rec_size, local_size_source, para_env_RPA)

      CALL timestop(handle)
   END SUBROUTINE matrix_operations_canonical

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param mat_P_omega_re ...
!> \param mat_P_omega_im ...
!> \param fm_mat_L_re ...
!> \param fm_mat_L_im ...
!> \param dimen_RI ...
! **************************************************************************************************
   SUBROUTINE compute_Q(cfm_mat_Q, mat_P_omega_re, mat_P_omega_im, fm_mat_L_re, fm_mat_L_im, dimen_RI)

      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(dbcsr_type), POINTER                          :: mat_P_omega_re, mat_P_omega_im
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L_re, fm_mat_L_im
      INTEGER, INTENT(IN)                                :: dimen_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Q', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_L, cfm_mat_work
      TYPE(cp_fm_type), POINTER                          :: fm_mat_work

      CALL timeset(routineN, handle)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_L)
      CALL cp_cfm_create(cfm_mat_L, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      NULLIFY (fm_mat_work)
      CALL cp_fm_create(fm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      ! 1. Copy dbcsr matrices mat_P_omega_re and mat_P_omega_im to cfm_mat_Q
      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_P_omega_re, fm_mat_work)
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_Q, z_one, fm_mat_work)

      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)
      CALL copy_dbcsr_to_fm(mat_P_omega_im, fm_mat_work)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_Q, gaussi, fm_mat_work)

      ! 2. Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      ! 3. work = P(iw,k)*L(k)
      CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                       z_zero, cfm_mat_work)

      ! 4. Q(iw,k) = L^T(k)*work
      CALL cp_cfm_gemm('C', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                       z_zero, cfm_mat_Q)

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_L)
      CALL cp_fm_release(fm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE compute_Q

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param para_env_RPA ...
!> \param trace_Qomega ...
!> \param dimen_RI ...
! **************************************************************************************************
   SUBROUTINE cholesky_decomp_Q(cfm_mat_Q, para_env_RPA, trace_Qomega, dimen_RI)

      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cholesky_decomp_Q', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, info_chol, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      ! get info of fm_mat_Q
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ! calculate the trace of Q and add 1 on the diagonal
      trace_Qomega = 0.0_dp
!$OMP     PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                 SHARED(ncol_local,nrow_local,col_indices,row_indices,trace_Qomega,cfm_mat_Q,dimen_RI)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               trace_Qomega(i_global) = REAL(cfm_mat_Q%local_data(iiB, jjB))
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB)+z_one
            END IF
         END DO
      END DO
      CALL mp_sum(trace_Qomega, para_env_RPA%group)

      CALL cp_cfm_cholesky_decompose(matrix=cfm_mat_Q, n=dimen_RI, info_out=info_chol)
      CPASSERT(info_chol == 0)

      CALL timestop(handle)

   END SUBROUTINE cholesky_decomp_Q

! **************************************************************************************************
!> \brief ...
!> \param Erpa ...
!> \param cfm_mat_Q ...
!> \param para_env_RPA ...
!> \param trace_Qomega ...
!> \param dimen_RI ...
!> \param freq_weight ...
!> \param kp_weight ...
! **************************************************************************************************
   SUBROUTINE frequency_and_kpoint_integration(Erpa, cfm_mat_Q, para_env_RPA, trace_Qomega, &
                                               dimen_RI, freq_weight, kp_weight)

      REAL(KIND=dp), INTENT(INOUT)                       :: Erpa
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: trace_Qomega
      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), INTENT(IN)                          :: freq_weight, kp_weight

      CHARACTER(LEN=*), PARAMETER :: routineN = 'frequency_and_kpoint_integration', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: FComega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Q_log

      CALL timeset(routineN, handle)

      ! get info of cholesky_decomposed(fm_mat_Q)
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ALLOCATE (Q_log(dimen_RI))
      Q_log = 0.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                SHARED(ncol_local,nrow_local,col_indices,row_indices,Q_log,cfm_mat_Q,dimen_RI)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               Q_log(i_global) = 2.0_dp*LOG(REAL(cfm_mat_Q%local_data(iiB, jjB)))
            END IF
         END DO
      END DO
      CALL mp_sum(Q_log, para_env_RPA%group)

      FComega = 0.0_dp
      DO iiB = 1, dimen_RI
         IF (MODULO(iiB, para_env_RPA%num_pe) /= para_env_RPA%mepos) CYCLE
         ! FComega=FComega+(LOG(Q_log(iiB))-trace_Qomega(iiB))/2.0_dp
         FComega = FComega+(Q_log(iiB)-trace_Qomega(iiB))/2.0_dp
      END DO

      Erpa = Erpa+FComega*freq_weight*kp_weight

      DEALLOCATE (Q_log)

      CALL timestop(handle)

   END SUBROUTINE frequency_and_kpoint_integration

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param mat_SinvVSinv ...
!> \param mat_3c_overl_nnP_ic ...
!> \param mat_3c_overl_nnP_ic_reflected ...
!> \param mat_work ...
!> \param matrix_s ...
!> \param gw_corr_lev_tot ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param nmo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param do_ic_opt_homo_lumo ...
!> \param fm_mat_Q ...
!> \param para_env ...
!> \param mp2_env ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_ic_correction(Eigenval, mat_SinvVSinv, mat_3c_overl_nnP_ic, &
                                      mat_3c_overl_nnP_ic_reflected, mat_work, matrix_s, gw_corr_lev_tot, &
                                      gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, dimen_RI, unit_nr, &
                                      do_ic_opt_homo_lumo, fm_mat_Q, para_env, mp2_env, &
                                      do_alpha, do_beta)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: Eigenval
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_SinvVSinv
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_3c_overl_nnP_ic, &
                                                            mat_3c_overl_nnP_ic_reflected
      TYPE(dbcsr_type), POINTER                          :: mat_work
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: gw_corr_lev_tot, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, homo, nmo, dimen_RI, &
                                                            unit_nr
      LOGICAL, INTENT(IN)                                :: do_ic_opt_homo_lumo
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_ic_correction', &
         routineP = moduleN//':'//routineN

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, n_level_gw, n_level_gw_ref, &
                                                            sc_iter
      LOGICAL                                            :: do_closed_shell, do_neaton, my_do_alpha, &
                                                            my_do_beta
      REAL(KIND=dp)                                      :: alpha, old_energy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff_homo, coeff_lumo, &
                                                            Delta_Sigma_Neaton, Eigenval_M_occ, &
                                                            Eigenval_M_virt, vec_P_occ, vec_P_virt
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_occ, fm_struct_virt
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_occ, fm_mat_M_virt, &
                                                            fm_mat_U_occ, fm_mat_U_virt
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_B_nm_P, mat_N_occ_dbcsr, &
                                                            mat_N_virt_dbcsr
      TYPE(dbcsr_type), POINTER                          :: coeff_homo_dbcsr, coeff_lumo_dbcsr, &
                                                            temp_RI_vector, vec_P_occ_dbcsr, &
                                                            vec_P_virt_dbcsr

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

      IF (PRESENT(do_beta)) THEN
         my_do_beta = do_beta
      ELSE
         my_do_beta = .FALSE.
      END IF

      do_closed_shell = .NOT. (my_do_alpha .OR. my_do_beta)

      do_neaton = .NOT. do_ic_opt_homo_lumo

      IF (do_neaton) THEN

         ALLOCATE (Delta_Sigma_Neaton(gw_corr_lev_tot))
         Delta_Sigma_Neaton = 0.0_dp

         DO n_level_gw = 1, gw_corr_lev_tot

            IF (n_level_gw > gw_corr_lev_occ) THEN
               alpha = -0.5_dp
            ELSE
               alpha = 0.5_dp
            END IF

            CALL dbcsr_multiply("N", "N", alpha, mat_SinvVSinv%matrix, &
                                mat_3c_overl_nnP_ic(n_level_gw)%matrix, 0.0_dp, mat_work)

            CALL dbcsr_trace(mat_work, &
                             mat_3c_overl_nnP_ic_reflected(n_level_gw)%matrix, &
                             Delta_Sigma_Neaton(n_level_gw))

         END DO

         IF (unit_nr > 0) THEN

            WRITE (unit_nr, *) ' '

            IF (do_closed_shell) THEN
               WRITE (unit_nr, '(T3,A)') 'Single-electron energies with image charge (ic) correction'
               WRITE (unit_nr, '(T3,A)') '----------------------------------------------------------'
            ELSE IF (my_do_alpha) THEN
               WRITE (unit_nr, '(T3,A)') 'Single-electron energies of alpha spins with image charge (ic) correction'
               WRITE (unit_nr, '(T3,A)') '-------------------------------------------------------------------------'
            ELSE IF (my_do_beta) THEN
               WRITE (unit_nr, '(T3,A)') 'Single-electron energies of beta spins with image charge (ic) correction'
               WRITE (unit_nr, '(T3,A)') '------------------------------------------------------------------------'
            END IF

            WRITE (unit_nr, *) ' '
            WRITE (unit_nr, '(T3,A)') 'Reference for the ic: Neaton et al., PRL 97, 216405 (2006)'
            WRITE (unit_nr, *) ' '

            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T14,2A)') 'MO     E_n before ic corr           Delta E_ic', &
               '    E_n after ic corr'

            DO n_level_gw = 1, gw_corr_lev_tot
               n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
               IF (n_level_gw <= gw_corr_lev_occ) THEN
                  occ_virt = 'occ'
               ELSE
                  occ_virt = 'vir'
               END IF

               WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')  ', &
                  Eigenval(n_level_gw_ref)*evolt, &
                  Delta_Sigma_Neaton(n_level_gw)*evolt, &
                  (Eigenval(n_level_gw_ref)+Delta_Sigma_Neaton(n_level_gw))*evolt

            END DO

            IF (do_closed_shell) THEN
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A,F57.2)') 'IC HOMO-LUMO gap (eV)', (Eigenval(homo+1)+ &
                                                                         Delta_Sigma_Neaton(gw_corr_lev_occ+1)- &
                                                                         Eigenval(homo)- &
                                                                         Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
            ELSE IF (my_do_alpha) THEN
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A,F51.2)') 'Alpha IC HOMO-LUMO gap (eV)', (Eigenval(homo+1)+ &
                                                                               Delta_Sigma_Neaton(gw_corr_lev_occ+1)- &
                                                                               Eigenval(homo)- &
                                                                               Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
            ELSE IF (my_do_beta) THEN
               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A,F52.2)') 'Beta IC HOMO-LUMO gap (eV)', (Eigenval(homo+1)+ &
                                                                              Delta_Sigma_Neaton(gw_corr_lev_occ+1)- &
                                                                              Eigenval(homo)- &
                                                                              Delta_Sigma_Neaton(gw_corr_lev_occ))*evolt
            END IF

            IF (mp2_env%ri_g0w0%print_ic_values) THEN

               WRITE (unit_nr, '(T3,A)') ' '
               WRITE (unit_nr, '(T3,A)') 'Horizontal list for copying the image charge corrections for use as input:'
               WRITE (unit_nr, '(*(F7.3))') (Delta_Sigma_Neaton(n_level_gw)*evolt, &
                                             n_level_gw=1, gw_corr_lev_tot)

            END IF

         END IF

         Eigenval(homo-gw_corr_lev_occ+1:homo+gw_corr_lev_virt) = Eigenval(homo-gw_corr_lev_occ+1: &
                                                                           homo+gw_corr_lev_virt) &
                                                                  +Delta_Sigma_Neaton(1:gw_corr_lev_tot)

      END IF

      IF (do_ic_opt_homo_lumo) THEN

         IF (SIZE(mp2_env%ri_g0w0%gw_eigenvalues) > 1) THEN

            CPASSERT(SIZE(mp2_env%ri_g0w0%gw_eigenvalues) == gw_corr_lev_tot)

            Eigenval(homo-gw_corr_lev_occ+1:homo+gw_corr_lev_virt) = mp2_env%ri_g0w0%gw_eigenvalues(1:gw_corr_lev_tot)

         END IF

         ALLOCATE (coeff_homo(homo))
         coeff_homo = 0.0_dp
         ! initial guess
         coeff_homo(homo) = 1.0_dp

         ALLOCATE (coeff_lumo(nmo-homo))
         coeff_lumo = 0.0_dp
         ! initial guess
         coeff_lumo(1) = 1.0_dp

         ALLOCATE (vec_P_occ(dimen_RI))
         ALLOCATE (vec_P_virt(dimen_RI))

         NULLIFY (mat_B_nm_P)
         CALL dbcsr_allocate_matrix_set(mat_B_nm_P, gw_corr_lev_tot)

         DO n_level_gw = 1, gw_corr_lev_tot

            ALLOCATE (mat_B_nm_P(n_level_gw)%matrix)
            CALL dbcsr_create(matrix=mat_B_nm_P(n_level_gw)%matrix, &
                              template=mat_3c_overl_nnP_ic(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)

            IF (n_level_gw > gw_corr_lev_occ) THEN
               alpha = -0.5_dp
            ELSE
               alpha = 0.5_dp
            END IF

            CALL dbcsr_multiply("N", "N", alpha, mat_SinvVSinv%matrix, &
                                mat_3c_overl_nnP_ic(n_level_gw)%matrix, 0.0_dp, &
                                mat_B_nm_P(n_level_gw)%matrix)

         END DO

         NULLIFY (coeff_homo_dbcsr)
         CALL dbcsr_init_p(coeff_homo_dbcsr)
         CALL cp_dbcsr_m_by_n_from_row_template(coeff_homo_dbcsr, template=matrix_s(1)%matrix, n=1, &
                                                sym=dbcsr_type_no_symmetry, data_type=dbcsr_type_real_default)
         CALL dbcsr_reserve_all_blocks(coeff_homo_dbcsr)
         CALL dbcsr_set(coeff_homo_dbcsr, 0.0_dp)

         NULLIFY (coeff_lumo_dbcsr)
         CALL dbcsr_init_p(coeff_lumo_dbcsr)
         CALL cp_dbcsr_m_by_n_from_row_template(coeff_lumo_dbcsr, template=matrix_s(1)%matrix, n=1, &
                                                sym=dbcsr_type_no_symmetry, data_type=dbcsr_type_real_default)
         CALL dbcsr_reserve_all_blocks(coeff_lumo_dbcsr)
         CALL dbcsr_set(coeff_lumo_dbcsr, 0.0_dp)

         NULLIFY (vec_P_occ_dbcsr)
         CALL dbcsr_init_p(vec_P_occ_dbcsr)
         CALL cp_dbcsr_m_by_n_from_row_template(vec_P_occ_dbcsr, template=mat_SinvVSinv%matrix, n=1, &
                                                sym=dbcsr_type_no_symmetry, data_type=dbcsr_type_real_default)
         CALL dbcsr_reserve_all_blocks(vec_P_occ_dbcsr)
         CALL dbcsr_set(vec_P_occ_dbcsr, 0.0_dp)

         NULLIFY (vec_P_virt_dbcsr)
         CALL dbcsr_init_p(vec_P_virt_dbcsr)
         CALL cp_dbcsr_m_by_n_from_row_template(vec_P_virt_dbcsr, template=mat_SinvVSinv%matrix, n=1, &
                                                sym=dbcsr_type_no_symmetry, data_type=dbcsr_type_real_default)
         CALL dbcsr_reserve_all_blocks(vec_P_virt_dbcsr)
         CALL dbcsr_set(vec_P_virt_dbcsr, 0.0_dp)

         NULLIFY (temp_RI_vector)
         CALL dbcsr_init_p(temp_RI_vector)
         CALL dbcsr_create(matrix=temp_RI_vector, &
                           template=vec_P_occ_dbcsr, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_N_occ_dbcsr)
         CALL dbcsr_allocate_matrix_set(mat_N_occ_dbcsr, gw_corr_lev_occ)

         DO n_level_gw = 1, gw_corr_lev_occ

            ALLOCATE (mat_N_occ_dbcsr(n_level_gw)%matrix)
            CALL dbcsr_create(matrix=mat_N_occ_dbcsr(n_level_gw)%matrix, &
                              template=coeff_homo_dbcsr, &
                              matrix_type=dbcsr_type_no_symmetry)

         END DO

         NULLIFY (mat_N_virt_dbcsr)
         CALL dbcsr_allocate_matrix_set(mat_N_virt_dbcsr, gw_corr_lev_virt)

         DO n_level_gw = 1, gw_corr_lev_virt

            ALLOCATE (mat_N_virt_dbcsr(n_level_gw)%matrix)
            CALL dbcsr_create(matrix=mat_N_virt_dbcsr(n_level_gw)%matrix, &
                              template=coeff_lumo_dbcsr, &
                              matrix_type=dbcsr_type_no_symmetry)

         END DO

         NULLIFY (fm_struct_occ)
         CALL cp_fm_struct_create(fm_struct_occ, context=fm_mat_Q%matrix_struct%context, nrow_global=gw_corr_lev_occ, &
                                  ncol_global=gw_corr_lev_occ, para_env=fm_mat_Q%matrix_struct%para_env)

         NULLIFY (fm_mat_M_occ)
         CALL cp_fm_create(fm_mat_M_occ, fm_struct_occ)
         CALL cp_fm_set_all(fm_mat_M_occ, 0.0_dp)

         NULLIFY (fm_mat_U_occ)
         CALL cp_fm_create(fm_mat_U_occ, fm_struct_occ)
         CALL cp_fm_set_all(fm_mat_U_occ, 0.0_dp)

         NULLIFY (fm_struct_virt)
         CALL cp_fm_struct_create(fm_struct_virt, context=fm_mat_Q%matrix_struct%context, nrow_global=gw_corr_lev_virt, &
                                  ncol_global=gw_corr_lev_virt, para_env=fm_mat_Q%matrix_struct%para_env)

         NULLIFY (fm_mat_M_virt)
         CALL cp_fm_create(fm_mat_M_virt, fm_struct_virt)
         CALL cp_fm_set_all(fm_mat_M_virt, 0.0_dp)

         NULLIFY (fm_mat_U_virt)
         CALL cp_fm_create(fm_mat_U_virt, fm_struct_virt)
         CALL cp_fm_set_all(fm_mat_U_virt, 0.0_dp)

         ALLOCATE (Eigenval_M_occ(gw_corr_lev_occ))

         ALLOCATE (Eigenval_M_virt(gw_corr_lev_virt))

         old_energy = Eigenval(homo)

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T40,F41.2)') 'Gas phase HOMO (eV): ', Eigenval(homo)*evolt

         DO sc_iter = 1, 20

            CALL fill_coeff_dbcsr_occ(coeff_homo_dbcsr, coeff_homo, homo)

            CALL dbcsr_set(vec_P_occ_dbcsr, 0.0_dp)

            ! build vec_P_occ
            DO n_level_gw = 1, gw_corr_lev_occ

               n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ

               CALL dbcsr_multiply("N", "N", 1.0_dp, mat_B_nm_P(n_level_gw)%matrix, &
                                   coeff_homo_dbcsr, 0.0_dp, temp_RI_vector)

               CALL dbcsr_add(vec_P_occ_dbcsr, temp_RI_vector, 1.0_dp, coeff_homo(n_level_gw_ref))

            END DO

            ! build N matrix
            DO n_level_gw = 1, gw_corr_lev_occ

               CALL dbcsr_multiply("T", "N", 1.0_dp, mat_3c_overl_nnP_ic_reflected(n_level_gw)%matrix, &
                                   vec_P_occ_dbcsr, 0.0_dp, mat_N_occ_dbcsr(n_level_gw)%matrix)

            END DO

            CALL fill_fm_mat_M_occ(fm_mat_M_occ, mat_N_occ_dbcsr, matrix_s, Eigenval, gw_corr_lev_occ, homo, para_env)

            IF (sc_iter == 1) THEN
               CALL print_Neaton_value(fm_mat_M_occ, unit_nr, gw_corr_lev_occ, para_env, do_homo=.TRUE.)
            END IF

            CALL cp_fm_scale(-1.0_dp, fm_mat_M_occ)

            CALL choose_eigv_solver(fm_mat_M_occ, fm_mat_U_occ, Eigenval_M_occ)

            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,I5,A,F35.2)') 'HOMO SCF Cycle: ', sc_iter, &
               '   HOMO with IC (eV): ', &
               -Eigenval_M_occ(1)*evolt

            IF (ABS(old_energy+Eigenval_M_occ(1)) < 1.0E-5) EXIT

            old_energy = -Eigenval_M_occ(1)

            CALL update_coeff_homo(coeff_homo, fm_mat_U_occ, para_env, homo, gw_corr_lev_occ)

         END DO

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T40,F41.2)') 'Gas phase LUMO (eV): ', Eigenval(homo+1)*evolt

         DO sc_iter = 1, 20

            CALL fill_coeff_dbcsr_virt(coeff_lumo_dbcsr, coeff_lumo, homo)

            CALL dbcsr_set(vec_P_virt_dbcsr, 0.0_dp)

            ! build vec_P_virt
            DO n_level_gw = gw_corr_lev_occ+1, gw_corr_lev_tot

               n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ

               CALL dbcsr_multiply("N", "N", 1.0_dp, mat_B_nm_P(n_level_gw)%matrix, &
                                   coeff_lumo_dbcsr, 0.0_dp, temp_RI_vector)

               CALL dbcsr_add(vec_P_virt_dbcsr, temp_RI_vector, 1.0_dp, coeff_lumo(n_level_gw_ref-homo))

            END DO

            ! build N matrix
            DO n_level_gw = 1, gw_corr_lev_virt

               CALL dbcsr_multiply("T", "N", 1.0_dp, mat_3c_overl_nnP_ic_reflected(n_level_gw+gw_corr_lev_occ)%matrix, &
                                   vec_P_virt_dbcsr, 0.0_dp, mat_N_virt_dbcsr(n_level_gw)%matrix)

            END DO

            CALL fill_fm_mat_M_virt(fm_mat_M_virt, mat_N_virt_dbcsr, matrix_s, Eigenval, gw_corr_lev_virt, homo, para_env)

            IF (sc_iter == 1) THEN
               CALL print_Neaton_value(fm_mat_M_virt, unit_nr, 1, para_env, do_lumo=.TRUE.)
            END IF

            CALL choose_eigv_solver(fm_mat_M_virt, fm_mat_U_virt, Eigenval_M_virt)

            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,I5,A,F35.2)') 'LUMO SCF Cycle: ', sc_iter, &
               '   LUMO with IC (eV): ', &
               Eigenval_M_virt(1)*evolt

            IF (ABS(old_energy-Eigenval_M_virt(1)) < 1.0E-5) EXIT

            old_energy = Eigenval_M_virt(1)

            CALL update_coeff_homo(coeff_lumo, fm_mat_U_virt, para_env, nmo-homo, &
                                   gw_corr_lev_virt, do_lumo=.TRUE.)

         END DO

         DEALLOCATE (coeff_homo, coeff_lumo)
         DEALLOCATE (Eigenval_M_occ)
         DEALLOCATE (vec_P_occ, vec_P_virt)
         CALL dbcsr_deallocate_matrix_set(mat_B_nm_P)
         CALL dbcsr_deallocate_matrix_set(mat_N_occ_dbcsr)
         CALL dbcsr_deallocate_matrix_set(mat_N_virt_dbcsr)
         CALL dbcsr_release_p(coeff_homo_dbcsr)
         CALL dbcsr_release_p(coeff_lumo_dbcsr)
         CALL dbcsr_release_p(vec_P_occ_dbcsr)
         CALL dbcsr_release_p(vec_P_virt_dbcsr)
         CALL dbcsr_release_p(temp_RI_vector)
         CALL cp_fm_release(fm_mat_M_occ)
         CALL cp_fm_release(fm_mat_M_virt)
         CALL cp_fm_release(fm_mat_U_occ)
         CALL cp_fm_release(fm_mat_U_virt)
         CALL cp_fm_struct_release(fm_struct_occ)
         CALL cp_fm_struct_release(fm_struct_virt)

      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_ic_correction

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_M_occ ...
!> \param unit_nr ...
!> \param gw_corr_lev_occ ...
!> \param para_env ...
!> \param do_homo ...
!> \param do_lumo ...
! **************************************************************************************************
   SUBROUTINE print_Neaton_value(fm_mat_M_occ, unit_nr, gw_corr_lev_occ, para_env, do_homo, do_lumo)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_occ
      INTEGER, INTENT(IN)                                :: unit_nr, gw_corr_lev_occ
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_homo, do_lumo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_Neaton_value', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_homo, my_do_lumo
      REAL(KIND=dp)                                      :: Neaton_value

      CALL timeset(routineN, handle)

      my_do_homo = .FALSE.
      IF (PRESENT(do_homo)) my_do_homo = do_homo

      my_do_lumo = .FALSE.
      IF (PRESENT(do_lumo)) my_do_lumo = do_lumo

      CPASSERT(my_do_homo .NEQV. my_do_lumo)

      CALL cp_fm_get_info(matrix=fm_mat_M_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      Neaton_value = 0.0_dp

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == gw_corr_lev_occ .AND. i_global == gw_corr_lev_occ) THEN
               Neaton_value = fm_mat_M_occ%local_data(iiB, jjB)
            END IF
         END DO
      END DO

      CALL mp_sum(Neaton_value, para_env%group)

      IF (unit_nr > 0 .AND. my_do_homo) WRITE (unit_nr, '(T3,A,F47.2)') 'Neaton value of the HOMO (eV): ', &
         Neaton_value*evolt
      IF (unit_nr > 0 .AND. my_do_lumo) WRITE (unit_nr, '(T3,A,F47.2)') 'Neaton value of the LUMO (eV): ', &
         Neaton_value*evolt

      CALL timestop(handle)

   END SUBROUTINE print_Neaton_value

! **************************************************************************************************
!> \brief ...
!> \param coeff_homo ...
!> \param fm_mat_U_occ ...
!> \param para_env ...
!> \param homo ...
!> \param gw_corr_lev_occ ...
!> \param do_lumo ...
! **************************************************************************************************
   SUBROUTINE update_coeff_homo(coeff_homo, fm_mat_U_occ, para_env, homo, gw_corr_lev_occ, do_lumo)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: coeff_homo
      TYPE(cp_fm_type), POINTER                          :: fm_mat_U_occ
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: homo, gw_corr_lev_occ
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_lumo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'update_coeff_homo', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_homo
      REAL(KIND=dp)                                      :: norm_coeff_homo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff_homo_update, &
                                                            coeff_homo_update_orthog

      CALL timeset(routineN, handle)

      my_do_homo = .TRUE.
      IF (PRESENT(do_lumo)) my_do_homo = .NOT. do_lumo

      CALL cp_fm_get_info(matrix=fm_mat_U_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ALLOCATE (coeff_homo_update(homo))
      coeff_homo_update = 0.0_dp

      ALLOCATE (coeff_homo_update_orthog(homo))
      coeff_homo_update_orthog = 0.0_dp

      ! take the eigenvector belongning to the largest eigenvalue
      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local

            j_global = col_indices(jjB)

            IF (j_global .NE. 1) CYCLE

            IF (my_do_homo) THEN
               coeff_homo_update(i_global+homo-gw_corr_lev_occ) = fm_mat_U_occ%local_data(iiB, jjB)
            ELSE
               coeff_homo_update(i_global) = fm_mat_U_occ%local_data(iiB, jjB)
            END IF
         END DO
      END DO

      CALL mp_sum(coeff_homo_update, para_env%group)

      norm_coeff_homo = NORM2(coeff_homo_update)

      coeff_homo(:) = coeff_homo_update(:)/norm_coeff_homo

      CALL timestop(handle)

   END SUBROUTINE update_coeff_homo

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_M_virt ...
!> \param mat_N_virt_dbcsr ...
!> \param matrix_s ...
!> \param Eigenval ...
!> \param gw_corr_lev_virt ...
!> \param homo ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fill_fm_mat_M_virt(fm_mat_M_virt, mat_N_virt_dbcsr, matrix_s, Eigenval, gw_corr_lev_virt, homo, para_env)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_virt
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_N_virt_dbcsr, matrix_s
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: gw_corr_lev_virt, homo
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_mat_M_virt', &
         routineP = moduleN//':'//routineN

      INTEGER :: col, col_offset, col_size, handle, i_col, i_global, i_index, i_row, iiB, imepos, &
         j_global, jjB, m_level_gw, n_level_gw, nblkrows_total, ncol_local, nfullrows_total, &
         nrow_local, offset, row, row_block, row_index, row_offset, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_from_indx, entry_counter, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_blk_offset, &
                                                            row_blk_sizes, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_tmp, matrix_tmp_2
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      NULLIFY (matrix_tmp)
      CALL dbcsr_init_p(matrix_tmp)
      CALL cp_dbcsr_m_by_n_from_template(matrix_tmp, template=matrix_s(1)%matrix, m=gw_corr_lev_virt, n=gw_corr_lev_virt, &
                                         sym=dbcsr_type_no_symmetry)

      NULLIFY (matrix_tmp_2)
      CALL dbcsr_init_p(matrix_tmp_2)
      CALL dbcsr_create(matrix=matrix_tmp_2, &
                        template=matrix_tmp, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_reserve_all_blocks(matrix_tmp)
      CALL dbcsr_set(matrix_tmp, 0.0_dp)

      CALL dbcsr_get_info(matrix_tmp, &
                          nblkrows_total=nblkrows_total, &
                          nfullrows_total=nfullrows_total, &
                          row_blk_offset=row_blk_offset, &
                          row_blk_size=row_blk_sizes)

      ALLOCATE (blk_from_indx(nfullrows_total))

      DO row_index = 1, nfullrows_total
         DO row_block = 1, nblkrows_total
            IF (row_index >= row_blk_offset(row_block) .AND. &
                row_index <= row_blk_offset(row_block)+row_blk_sizes(row_block)-1) THEN

               blk_from_indx(row_index) = row_block
            END IF
         END DO
      END DO

      ALLOCATE (num_entries_send(0:para_env%num_pe-1))
      num_entries_send = 0
      ALLOCATE (num_entries_rec(0:para_env%num_pe-1))
      num_entries_rec = 0

      DO n_level_gw = 1, gw_corr_lev_virt
         CALL dbcsr_iterator_start(iter, mat_N_virt_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            DO i_row = 1, row_size
               m_level_gw = row_offset-1+i_row-homo

               IF (m_level_gw < 1) CYCLE
               IF (m_level_gw > gw_corr_lev_virt) CYCLE

               CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                 blk_from_indx(n_level_gw), imepos)

               num_entries_send(imepos) = num_entries_send(imepos)+1
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

      ALLOCATE (buffer_rec(0:para_env%num_pe-1))
      ALLOCATE (buffer_send(0:para_env%num_pe-1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe-1))
      entry_counter(:) = 1

      DO n_level_gw = 1, gw_corr_lev_virt

         CALL dbcsr_iterator_start(iter, mat_N_virt_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            DO i_row = 1, row_size
               m_level_gw = row_offset-1+i_row-homo

               IF (m_level_gw < 1) CYCLE
               IF (m_level_gw > gw_corr_lev_virt) CYCLE

               CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                 blk_from_indx(n_level_gw), imepos)

               offset = entry_counter(imepos)

               buffer_send(imepos)%msg(offset) = data_block(i_row, 1)
               buffer_send(imepos)%indx(offset, 1) = m_level_gw
               buffer_send(imepos)%indx(offset, 2) = n_level_gw

               entry_counter(imepos) = entry_counter(imepos)+1
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array)

      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size, &
                                        col_offset=col_offset, col_size=col_size)

         DO imepos = 0, para_env%num_pe-1
            DO i_index = 1, num_entries_rec(imepos)
               IF (buffer_rec(imepos)%indx(i_index, 1) >= row_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 1) <= row_offset+row_size-1 .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) >= col_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) <= col_offset+col_size-1) THEN

                  i_row = buffer_rec(imepos)%indx(i_index, 1)-row_offset+1
                  i_col = buffer_rec(imepos)%indx(i_index, 2)-col_offset+1

                  data_block(i_row, i_col) = buffer_rec(imepos)%msg(i_index)
               END IF
            END DO
         END DO
      END DO

      CALL dbcsr_iterator_stop(iter)

      ! symmetrize the result
      CALL dbcsr_transposed(matrix_tmp_2, matrix_tmp)
      CALL dbcsr_add(matrix_tmp, matrix_tmp_2, 0.5_dp, 0.5_dp)

      CALL copy_dbcsr_to_fm(matrix_tmp, fm_mat_M_virt)

      ! add the eigenvalue on the diag of M
      CALL cp_fm_get_info(matrix=fm_mat_M_virt, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global) THEN
               fm_mat_M_virt%local_data(iiB, jjB) = fm_mat_M_virt%local_data(iiB, jjB)+ &
                                                    Eigenval(i_global+homo)
            END IF
         END DO
      END DO

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      CALL dbcsr_release_p(matrix_tmp)
      CALL dbcsr_release_p(matrix_tmp_2)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_mat_M_virt

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_M_occ ...
!> \param mat_N_occ_dbcsr ...
!> \param matrix_s ...
!> \param Eigenval ...
!> \param gw_corr_lev_occ ...
!> \param homo ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fill_fm_mat_M_occ(fm_mat_M_occ, mat_N_occ_dbcsr, matrix_s, Eigenval, gw_corr_lev_occ, homo, para_env)
      TYPE(cp_fm_type), POINTER                          :: fm_mat_M_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_N_occ_dbcsr, matrix_s
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, homo
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_mat_M_occ', &
         routineP = moduleN//':'//routineN

      INTEGER :: col, col_offset, col_size, handle, i_col, i_global, i_index, i_row, iiB, imepos, &
         j_global, jjB, m_level_gw, n_level_gw, nblkrows_total, ncol_local, nfullrows_total, &
         nrow_local, offset, row, row_block, row_index, row_offset, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_from_indx, entry_counter, &
                                                            num_entries_rec, num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_blk_offset, &
                                                            row_blk_sizes, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_tmp
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      NULLIFY (matrix_tmp)
      CALL dbcsr_init_p(matrix_tmp)
      CALL cp_dbcsr_m_by_n_from_template(matrix_tmp, template=matrix_s(1)%matrix, m=gw_corr_lev_occ, n=gw_corr_lev_occ, &
                                         sym=dbcsr_type_no_symmetry)

      CALL dbcsr_reserve_all_blocks(matrix_tmp)
      CALL dbcsr_set(matrix_tmp, 0.0_dp)

      CALL dbcsr_get_info(matrix_tmp, &
                          nblkrows_total=nblkrows_total, &
                          nfullrows_total=nfullrows_total, &
                          row_blk_offset=row_blk_offset, &
                          row_blk_size=row_blk_sizes)

      ALLOCATE (blk_from_indx(nfullrows_total))

      DO row_index = 1, nfullrows_total
         DO row_block = 1, nblkrows_total
            IF (row_index >= row_blk_offset(row_block) .AND. &
                row_index <= row_blk_offset(row_block)+row_blk_sizes(row_block)-1) THEN

               blk_from_indx(row_index) = row_block
            END IF
         END DO
      END DO

      ALLOCATE (num_entries_send(0:para_env%num_pe-1))
      num_entries_send = 0
      ALLOCATE (num_entries_rec(0:para_env%num_pe-1))
      num_entries_rec = 0

      DO n_level_gw = 1, gw_corr_lev_occ

         CALL dbcsr_iterator_start(iter, mat_N_occ_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            IF (row_offset+row_size-1 <= homo) THEN
               DO i_row = 1, row_size
                  m_level_gw = row_offset-1+i_row-(homo-gw_corr_lev_occ)

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  num_entries_send(imepos) = num_entries_send(imepos)+1
               END DO
            ELSE IF (row_offset <= homo) THEN
               DO m_level_gw = row_offset-(homo-gw_corr_lev_occ), gw_corr_lev_occ

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  num_entries_send(imepos) = num_entries_send(imepos)+1
               END DO
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      CALL mp_alltoall(num_entries_send, num_entries_rec, 1, para_env%group)

      ALLOCATE (buffer_rec(0:para_env%num_pe-1))
      ALLOCATE (buffer_send(0:para_env%num_pe-1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_rec(imepos), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_send(imepos), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe-1))
      entry_counter(:) = 1

      DO n_level_gw = 1, gw_corr_lev_occ

         CALL dbcsr_iterator_start(iter, mat_N_occ_dbcsr(n_level_gw)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_offset=row_offset, row_size=row_size)

            IF (row_offset+row_size-1 <= homo) THEN
               DO i_row = 1, row_size
                  m_level_gw = row_offset-1+i_row-(homo-gw_corr_lev_occ)

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  offset = entry_counter(imepos)

                  buffer_send(imepos)%msg(offset) = data_block(i_row, 1)
                  buffer_send(imepos)%indx(offset, 1) = m_level_gw
                  buffer_send(imepos)%indx(offset, 2) = n_level_gw

                  entry_counter(imepos) = entry_counter(imepos)+1
               END DO
            ELSE IF (row_offset <= homo) THEN
               DO m_level_gw = row_offset-(homo-gw_corr_lev_occ), gw_corr_lev_occ

                  IF (m_level_gw < 1) CYCLE

                  CALL dbcsr_get_stored_coordinates(matrix_tmp, blk_from_indx(m_level_gw), &
                                                    blk_from_indx(n_level_gw), imepos)

                  offset = entry_counter(imepos)

                  i_row = m_level_gw+(homo-gw_corr_lev_occ)-row_offset+1

                  buffer_send(imepos)%msg(offset) = data_block(i_row, 1)
                  buffer_send(imepos)%indx(offset, 1) = m_level_gw
                  buffer_send(imepos)%indx(offset, 2) = n_level_gw

                  entry_counter(imepos) = entry_counter(imepos)+1
               END DO
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      CALL communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array)

      CALL dbcsr_iterator_start(iter, matrix_tmp)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size, &
                                        col_offset=col_offset, col_size=col_size)

         DO imepos = 0, para_env%num_pe-1
            DO i_index = 1, num_entries_rec(imepos)
               IF (buffer_rec(imepos)%indx(i_index, 1) >= row_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 1) <= row_offset+row_size-1 .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) >= col_offset .AND. &
                   buffer_rec(imepos)%indx(i_index, 2) <= col_offset+col_size-1) THEN

                  i_row = buffer_rec(imepos)%indx(i_index, 1)-row_offset+1
                  i_col = buffer_rec(imepos)%indx(i_index, 2)-col_offset+1

                  data_block(i_row, i_col) = buffer_rec(imepos)%msg(i_index)
               END IF
            END DO
         END DO
      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL copy_dbcsr_to_fm(matrix_tmp, fm_mat_M_occ)

      ! add the eigenvalue on the diag of M
      CALL cp_fm_get_info(matrix=fm_mat_M_occ, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global) THEN
               fm_mat_M_occ%local_data(iiB, jjB) = fm_mat_M_occ%local_data(iiB, jjB)+ &
                                                   Eigenval(i_global+homo-gw_corr_lev_occ)
            END IF
         END DO
      END DO

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      CALL dbcsr_release_p(matrix_tmp)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_mat_M_occ

! **************************************************************************************************
!> \brief ...
!> \param coeff_dbcsr ...
!> \param coeff ...
!> \param homo ...
! **************************************************************************************************
   SUBROUTINE fill_coeff_dbcsr_occ(coeff_dbcsr, coeff, homo)
      TYPE(dbcsr_type), POINTER                          :: coeff_dbcsr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coeff
      INTEGER, INTENT(IN)                                :: homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_coeff_dbcsr_occ', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, end_data_block, handle, row, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, coeff_dbcsr)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size)

         IF (row_offset+row_size-1 <= homo) THEN
            data_block(1:row_size, 1) = coeff(row_offset:row_offset+row_size-1)
         ELSE IF (row_offset <= homo) THEN
            end_data_block = homo-row_offset+1

            data_block(1:end_data_block, 1) = coeff(row_offset:row_offset+end_data_block-1)
         END IF
      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE fill_coeff_dbcsr_occ

! **************************************************************************************************
!> \brief ...
!> \param coeff_dbcsr ...
!> \param coeff ...
!> \param homo ...
! **************************************************************************************************
   SUBROUTINE fill_coeff_dbcsr_virt(coeff_dbcsr, coeff, homo)
      TYPE(dbcsr_type), POINTER                          :: coeff_dbcsr
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coeff
      INTEGER, INTENT(IN)                                :: homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_coeff_dbcsr_virt', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, handle, row, row_offset, row_size, &
                                                            start_data_block
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, coeff_dbcsr)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_offset=row_offset, row_size=row_size)

         IF (row_offset > homo) THEN
            data_block(1:row_size, 1) = coeff(row_offset-homo:row_offset+row_size-homo-1)
         ELSE IF (row_offset+row_size-1 > homo) THEN
            start_data_block = homo-row_offset+2

            data_block(start_data_block:row_size, 1) = coeff(1:row_offset+row_size-homo)
         END IF
      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE fill_coeff_dbcsr_virt

! **************************************************************************************************
!> \brief ...
!> \param a_scaling_ext ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param tj_ext ...
!> \param wj_ext ...
!> \param fm_mat_S ...
!> \param homo_beta ...
!> \param virtual_beta ...
!> \param dimen_ia_beta ...
!> \param Eigenval_beta ...
!> \param fm_mat_S_beta ...
! **************************************************************************************************
   SUBROUTINE calc_scaling_factor(a_scaling_ext, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                  num_integ_points, num_integ_group, color_rpa_group, &
                                  tj_ext, wj_ext, fm_mat_S, &
                                  homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta)
      REAL(KIND=dp), INTENT(INOUT)                       :: a_scaling_ext
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      INTEGER, INTENT(IN)                                :: homo, virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj_ext, wj_ext
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      INTEGER, INTENT(IN), OPTIONAL                      :: homo_beta, virtual_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: Eigenval_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_scaling_factor', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icycle, jquad, nrow_local, &
                                                            nrow_local_beta
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp) :: a_high, a_low, a_scaling, conv_param, eps, first_deriv, left_term, &
         right_term, right_term_ref, right_term_ref_beta, step
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cottj, D_ia, D_ia_beta, iaia_RI, &
                                                            iaia_RI_beta, M_ia, M_ia_beta
      TYPE(cp_para_env_type), POINTER                    :: para_env_row, para_env_row_beta

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(homo_beta) .AND. &
          PRESENT(virtual_beta) .AND. &
          PRESENT(dimen_ia_beta) .AND. &
          PRESENT(Eigenval_beta) .AND. &
          PRESENT(fm_mat_S_beta)) my_open_shell = .TRUE.

      eps = 1.0E-10_dp

      ALLOCATE (cottj(num_integ_points))

      ! calculate the cotangent of the abscissa tj
      DO jquad = 1, num_integ_points
         cottj(jquad) = 1.0_dp/TAN(tj_ext(jquad))
      END DO

      CALL calc_ia_ia_integrals(para_env_RPA, homo, virtual, nrow_local, right_term_ref, Eigenval, D_ia, iaia_RI, M_ia, fm_mat_S, &
                                para_env_row)

      ! In the open shell case do point 1-2-3 for the beta spin
      IF (my_open_shell) THEN
  CALL calc_ia_ia_integrals(para_env_RPA, homo_beta, virtual_beta, nrow_local_beta, right_term_ref_beta, Eigenval_beta, D_ia_beta, &
                                   iaia_RI_beta, M_ia_beta, fm_mat_S_beta, para_env_row_beta)

         right_term_ref = right_term_ref+right_term_ref_beta
      END IF

      ! bcast the result
      IF (para_env%mepos == 0) THEN
         CALL mp_bcast(right_term_ref, 0, para_env%group)
      ELSE
         right_term_ref = 0.0_dp
         CALL mp_bcast(right_term_ref, 0, para_env%group)
      END IF

      ! 5) start iteration for solving the non-linear equation by bisection
      ! find limit, here step=0.5 seems a good compromise
      conv_param = 100.0_dp*EPSILON(right_term_ref)
      step = 0.5_dp
      a_low = 0.0_dp
      a_high = step
      right_term = -right_term_ref
      DO icycle = 1, num_integ_points*2
         a_scaling = a_high

         CALL calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                M_ia, cottj, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_row, para_env_row_beta)
         left_term = left_term/4.0_dp/pi*a_scaling

         IF (ABS(left_term) > ABS(right_term) .OR. ABS(left_term+right_term) <= conv_param) EXIT
         a_low = a_high
         a_high = a_high+step

      END DO

      IF (ABS(left_term+right_term) >= conv_param) THEN
         IF (a_scaling >= 2*num_integ_points*step) THEN
            a_scaling = 1.0_dp
         ELSE

            DO icycle = 1, num_integ_points*2
               a_scaling = (a_low+a_high)/2.0_dp

               CALL calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                      M_ia, cottj, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                      nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                      para_env, para_env_row, para_env_row_beta)
               left_term = left_term/4.0_dp/pi*a_scaling

               IF (ABS(left_term) > ABS(right_term)) THEN
                  a_high = a_scaling
               ELSE
                  a_low = a_scaling
               END IF

               IF (ABS(a_high-a_low) < 1.0e-5_dp) EXIT

            END DO

         END IF
      END IF

      a_scaling_ext = a_scaling
      CALL mp_bcast(a_scaling_ext, 0, para_env%group)

      DEALLOCATE (cottj)
      DEALLOCATE (iaia_RI)
      DEALLOCATE (D_ia)
      DEALLOCATE (M_ia)
      CALL cp_para_env_release(para_env_row)

      IF (my_open_shell) THEN
         DEALLOCATE (iaia_RI_beta)
         DEALLOCATE (D_ia_beta)
         DEALLOCATE (M_ia_beta)
         CALL cp_para_env_release(para_env_row_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_scaling_factor

! **************************************************************************************************
!> \brief ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param nrow_local ...
!> \param right_term_ref ...
!> \param Eigenval ...
!> \param D_ia ...
!> \param iaia_RI ...
!> \param M_ia ...
!> \param fm_mat_S ...
!> \param para_env_row ...
! **************************************************************************************************
   SUBROUTINE calc_ia_ia_integrals(para_env_RPA, homo, virtual, nrow_local, right_term_ref, Eigenval, &
                                   D_ia, iaia_RI, M_ia, fm_mat_S, para_env_row)

      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      INTEGER, INTENT(IN)                                :: homo, virtual
      INTEGER, INTENT(OUT)                               :: nrow_local
      REAL(KIND=dp), INTENT(OUT)                         :: right_term_ref
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: D_ia, iaia_RI, M_ia
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      TYPE(cp_para_env_type), POINTER                    :: para_env_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ia_ia_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: avirt, color_col, color_row, comm_col, &
                                                            comm_row, handle, i_global, iiB, iocc, &
                                                            jjB, ncol_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigen_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: iaia_RI_dp
      TYPE(cp_para_env_type), POINTER                    :: para_env_col

      CALL timeset(routineN, handle)

      ! calculate the (ia|ia) RI integrals
      ! ----------------------------------
      ! 1) get info fm_mat_S
      !XXX CALL cp_fm_to_fm(source=fm_mat_S,destination=fm_mat_G)
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! allocate the local buffer of iaia_RI integrals (dp kind)
      ALLOCATE (iaia_RI_dp(nrow_local))
      iaia_RI_dp = 0.0_dp

      ! 2) perform the local multiplication SUM_K (ia|K)*(ia|K)
      DO jjB = 1, ncol_local
         DO iiB = 1, nrow_local
            iaia_RI_dp(iiB) = iaia_RI_dp(iiB)+fm_mat_S%local_data(iiB, jjB)*fm_mat_S%local_data(iiB, jjB)
         END DO
      END DO

      ! 3) sum the result with the processes of the RPA_group having the same rows
      !          _______K_______               _
      !         |   |   |   |   |             | |
      !     --> | 1 | 5 | 9 | 13|   SUM -->   | |
      !         |___|__ |___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 2 | 6 | 10| 14|   SUM -->   | |
      !      ia |___|___|___|___|             |_|   (ia|ia)_RI
      !         |   |   |   |   |             | |
      !     --> | 3 | 7 | 11| 15|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 4 | 8 | 12| 16|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !

      color_row = fm_mat_S%matrix_struct%context%mepos(1)
      CALL mp_comm_split_direct(para_env_RPA%group, comm_row, color_row)
      NULLIFY (para_env_row)
      CALL cp_para_env_create(para_env_row, comm_row)

      CALL mp_sum(iaia_RI_dp, para_env_row%group)

      ! convert the iaia_RI_dp into double-double precision
      ALLOCATE (iaia_RI(nrow_local))
      DO iiB = 1, nrow_local
         iaia_RI(iiB) = iaia_RI_dp(iiB)
      END DO
      DEALLOCATE (iaia_RI_dp)

      ! 4) calculate the right hand term, D_ia is the matrix containing the
      ! orbital energy differences, M_ia is the diagonal of the full RPA 'excitation'
      ! matrix
      ALLOCATE (D_ia(nrow_local))

      ALLOCATE (M_ia(nrow_local))

      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)

         iocc = MAX(1, i_global-1)/virtual+1
         avirt = i_global-(iocc-1)*virtual
         eigen_diff = Eigenval(avirt+homo)-Eigenval(iocc)

         D_ia(iiB) = eigen_diff
      END DO

      DO iiB = 1, nrow_local
         M_ia(iiB) = D_ia(iiB)*D_ia(iiB)+2.0_dp*D_ia(iiB)*iaia_RI(iiB)
      END DO

      right_term_ref = 0.0_dp
      DO iiB = 1, nrow_local
         right_term_ref = right_term_ref+(SQRT(M_ia(iiB))-D_ia(iiB)-iaia_RI(iiB))
      END DO
      right_term_ref = right_term_ref/2.0_dp

      ! sum the result with the processes of the RPA_group having the same col
      color_col = fm_mat_S%matrix_struct%context%mepos(2)
      CALL mp_comm_split_direct(para_env_RPA%group, comm_col, color_col)
      NULLIFY (para_env_col)
      CALL cp_para_env_create(para_env_col, comm_col)

      ! allocate communication array for columns
      CALL mp_sum(right_term_ref, para_env_col%group)

      CALL cp_para_env_release(para_env_col)

      CALL timestop(handle)

   END SUBROUTINE calc_ia_ia_integrals

! **************************************************************************************************
!> \brief ...
!> \param a_scaling ...
!> \param left_term ...
!> \param first_deriv ...
!> \param num_integ_points ...
!> \param my_open_shell ...
!> \param M_ia ...
!> \param cottj ...
!> \param wj ...
!> \param D_ia ...
!> \param D_ia_beta ...
!> \param M_ia_beta ...
!> \param nrow_local ...
!> \param nrow_local_beta ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param para_env ...
!> \param para_env_row ...
!> \param para_env_row_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                M_ia, cottj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_row, para_env_row_beta)
      REAL(KIND=dp), INTENT(IN)                          :: a_scaling
      REAL(KIND=dp), INTENT(INOUT)                       :: left_term, first_deriv
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: my_open_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: M_ia, cottj, wj, D_ia, D_ia_beta, &
                                                            M_ia_beta
      INTEGER, INTENT(IN)                                :: nrow_local, nrow_local_beta, &
                                                            num_integ_group, color_rpa_group
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_row, para_env_row_beta

      INTEGER                                            :: iiB, jquad
      REAL(KIND=dp)                                      :: first_deriv_beta, left_term_beta, omega

      left_term = 0.0_dp
      first_deriv = 0.0_dp
      left_term_beta = 0.0_dp
      first_deriv_beta = 0.0_dp
      DO jquad = 1, num_integ_points
         ! parallelize over integration points
         IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE
         omega = a_scaling*cottj(jquad)

         DO iiB = 1, nrow_local
            ! parallelize over ia elements in the para_env_row group
            IF (MODULO(iiB, para_env_row%num_pe) /= para_env_row%mepos) CYCLE
            ! calculate left_term
            left_term = left_term+wj(jquad)* &
                        (LOG(1.0_dp+(M_ia(iiB)-D_ia(iiB)**2)/(omega**2+D_ia(iiB)**2))- &
                         (M_ia(iiB)-D_ia(iiB)**2)/(omega**2+D_ia(iiB)**2))
            first_deriv = first_deriv+wj(jquad)*cottj(jquad)**2* &
                          ((-M_ia(iiB)+D_ia(iiB)**2)**2/((omega**2+D_ia(iiB)**2)**2*(omega**2+M_ia(iiB))))
         END DO

         IF (my_open_shell) THEN
            DO iiB = 1, nrow_local_beta
               ! parallelize over ia elements in the para_env_row group
               IF (MODULO(iiB, para_env_row_beta%num_pe) /= para_env_row_beta%mepos) CYCLE
               ! calculate left_term
               left_term_beta = left_term_beta+wj(jquad)* &
                                (LOG(1.0_dp+(M_ia_beta(iiB)-D_ia_beta(iiB)**2)/(omega**2+D_ia_beta(iiB)**2))- &
                                 (M_ia_beta(iiB)-D_ia_beta(iiB)**2)/(omega**2+D_ia_beta(iiB)**2))
               first_deriv_beta = &
                  first_deriv_beta+wj(jquad)*cottj(jquad)**2* &
                  ((-M_ia_beta(iiB)+D_ia_beta(iiB)**2)**2/((omega**2+D_ia_beta(iiB)**2)**2*(omega**2+M_ia_beta(iiB))))
            END DO
         END IF

      END DO

      ! sum the contribution from all proc, starting form the row group
      CALL mp_sum(left_term, para_env%group)
      CALL mp_sum(first_deriv, para_env%group)

      IF (my_open_shell) THEN
         CALL mp_sum(left_term_beta, para_env%group)
         CALL mp_sum(first_deriv_beta, para_env%group)

         left_term = left_term+left_term_beta
         first_deriv = first_deriv+first_deriv_beta
      END IF

   END SUBROUTINE calculate_objfunc

! **************************************************************************************************
!> \brief Calculate integration weights for the tau grid (in dependency of the omega node)
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param weights_cos_tf_t_to_w ...
!> \param omega_tj ...
!> \param E_min ...
!> \param E_max ...
!> \param max_error ...
!> \param num_points_per_magnitude ...
! **************************************************************************************************
   SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_l_sq_wghts_cos_tf_t_to_w', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      CALL timeset(routineN, handle)

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min))+1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes))
      x_values = 0.0_dp
      ALLOCATE (y_values(num_x_nodes))
      y_values = 0.0_dp
      ALLOCATE (mat_A(num_x_nodes, num_integ_points))
      mat_A = 0.0_dp
      ALLOCATE (tau_wj_work(num_integ_points))
      tau_wj_work = 0.0_dp
      ALLOCATE (work_array(2*num_integ_points))
      work_array = 0.0_dp
      ALLOCATE (sing_values(num_integ_points))
      sing_values = 0.0_dp
      ALLOCATE (mat_U(num_x_nodes, num_x_nodes))
      mat_U = 0.0_dp
      ALLOCATE (mat_SinvVSinvT(num_x_nodes, num_integ_points))

      mat_SinvVSinvT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points+12*num_integ_points+2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_SinvVSinvSigma(num_integ_points, num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

      max_error = 0.0_dp

      ! loop over all omega frequency points
      DO jquad = 1, num_integ_points

         chi2_min_jquad = 100.0_dp

         ! set the x-values logarithmically in the interval [Emin,Emax]
         multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp)-1.0_dp))
         DO iii = 1, num_x_nodes
            x_values(iii) = E_min*multiplicator**(iii-1)
         END DO

         omega = omega_tj(jquad)

         ! y=2x/(x^2+omega_k^2)
         DO iii = 1, num_x_nodes
            y_values(iii) = 2.0_dp*x_values(iii)/((x_values(iii))**2+omega**2)
         END DO

         ! calculate mat_A
         DO jjj = 1, num_integ_points
            DO iii = 1, num_x_nodes
               mat_A(iii, jjj) = COS(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
            END DO
         END DO

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)

         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, tau_wj_work, num_integ_points)

         weights_cos_tf_t_to_w(jquad, :) = tau_wj_work(:)

         CALL calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                      y_values, num_integ_points, num_x_nodes)

      END DO ! jquad

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, work_array, sing_values, mat_U, mat_SinvVSinvT, &
                  work, iwork, mat_SinvVSinvSigma, vec_UTy)

      CALL timestop(handle)

   END SUBROUTINE get_l_sq_wghts_cos_tf_t_to_w

! **************************************************************************************************
!> \brief ...
!> \param max_error ...
!> \param omega ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param x_values ...
!> \param y_values ...
!> \param num_integ_points ...
!> \param num_x_nodes ...
! **************************************************************************************************
   PURE SUBROUTINE calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                           y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_tau_grid_with_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, max_error_tmp

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_cosine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk)-func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk)-func_val)
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

   END SUBROUTINE calc_max_error_fit_tau_grid_with_cosine

! **************************************************************************************************
!> \brief ...
!> \param max_error ...
!> \param omega ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param x_values ...
!> \param y_values ...
!> \param num_integ_points ...
!> \param num_x_nodes ...
! **************************************************************************************************
   PURE SUBROUTINE calc_max_error_fit_tau_grid_with_sine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                         y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp), INTENT(INOUT)                       :: max_error
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER, INTENT(IN)                                :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_tau_grid_with_sine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: kkk
      REAL(KIND=dp)                                      :: func_val, max_error_tmp

      max_error_tmp = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_sine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk)-func_val) > max_error_tmp) THEN
            max_error_tmp = ABS(y_values(kkk)-func_val)
         END IF

      END DO

      IF (max_error_tmp > max_error) THEN

         max_error = max_error_tmp

      END IF

   END SUBROUTINE calc_max_error_fit_tau_grid_with_sine

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for cosine transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_tau_grid_cosine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp), INTENT(OUT)                         :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_tau_grid_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val+tau_wj_work(iii)*COS(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE eval_fit_func_tau_grid_cosine

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for sine transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   PURE SUBROUTINE eval_fit_func_tau_grid_sine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp), INTENT(OUT)                         :: func_val
      REAL(KIND=dp), INTENT(IN)                          :: x_value
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_tj, tau_wj_work
      REAL(KIND=dp), INTENT(IN)                          :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_tau_grid_sine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val+tau_wj_work(iii)*SIN(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE eval_fit_func_tau_grid_sine

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

      ALLOCATE (non_zero_blocks_3c_tmp(1:cut_RI, 1:nblkrows_total, 0:(para_env_sub%num_pe-1)))
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

      DO imepos = 0, para_env_sub%num_pe-1
         DO i_cut_RI = 1, cut_RI
            maxlength_tmp = 0
            DO iblk = 1, nblkrows_total
               IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                  maxlength_tmp = maxlength_tmp+1
               END IF
            END DO
            IF (maxlength_tmp > maxlength) THEN
               maxlength = maxlength_tmp
            END IF
         END DO
      END DO

      ! save memory with smaller non_zero_blocks_3c
      ALLOCATE (non_zero_blocks_3c(1:cut_RI, 1:maxlength, 0:(para_env_sub%num_pe-1)))
      non_zero_blocks_3c = 0

      DO imepos = 0, para_env_sub%num_pe-1
         DO i_cut_RI = 1, cut_RI
            block_counter = 0
            DO iblk = 1, nblkrows_total
               IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                  block_counter = block_counter+1
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

      ALLOCATE (non_zero_blocks_3c_tmp(1:cut_RI, 1:nblkrows_total_max, 0:(para_env_sub%num_pe-1)))
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

         DO imepos = 0, para_env_sub%num_pe-1
            DO i_cut_RI = 1, cut_RI
               maxlength_tmp = 0
               DO iblk = 1, nblkrows_total
                  IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                     maxlength_tmp = maxlength_tmp+1
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
      ALLOCATE (non_zero_blocks_3c_cut(1:cut_RI, 1:maxlength, 0:(para_env_sub%num_pe-1), 1:cut_memory))
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

         DO imepos = 0, para_env_sub%num_pe-1
            DO i_cut_RI = 1, cut_RI
               block_counter = 0
               DO iblk = 1, nblkrows_total
                  IF (non_zero_blocks_3c_tmp(i_cut_RI, iblk, imepos) .NE. 0) THEN
                     block_counter = block_counter+1
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

                     occ_local_sum_j_cell = occ_local_sum_j_cell+occ_local

                  END DO ! j_cell

               END DO ! i_cut_RI

            END DO ! i_cell
         END DO ! i_mem

      END IF

      CALL timestop(handle)

   END SUBROUTINE check_sparsity_arrays_for_kp

! **************************************************************************************************
!> \brief test the singular value decomposition for the computation of integration weights for the
!>         Fourier transform between time and frequency grid in cubic-scaling RPA
!> \param nR ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE test_least_square_ft(nR, iw)
      INTEGER, INTENT(IN)                                :: nR, iw

      INTEGER                                            :: ierr, iR, jquad, num_integ_points
      REAL(KIND=dp)                                      :: max_error, multiplicator, Rc, Rc_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj, tj, wj, x_tw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_t_to_w

      Rc_max = 1.0E+7

      multiplicator = Rc_max**(1.0_dp/(REAL(nR, KIND=dp)-1.0_dp))

      DO num_integ_points = 1, 20

         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ALLOCATE (tau_tj(num_integ_points))
         tau_tj = 0.0_dp
         ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
         weights_cos_tf_t_to_w = 0.0_dp
         ALLOCATE (tau_wj(num_integ_points))
         tau_wj = 0.0_dp
         ALLOCATE (tj(num_integ_points))
         tj = 0.0_dp
         ALLOCATE (wj(num_integ_points))
         wj = 0.0_dp

         DO iR = 0, nR-1

            Rc = 2.0_dp*multiplicator**iR

            ierr = 0
            CALL get_rpa_minimax_coeff(num_integ_points, Rc, x_tw, ierr, print_warning=.FALSE.)

            DO jquad = 1, num_integ_points
               tj(jquad) = x_tw(jquad)
               wj(jquad) = x_tw(jquad+num_integ_points)
            END DO

            x_tw = 0.0_dp

            CALL get_exp_minimax_coeff(num_integ_points, Rc, x_tw)

            DO jquad = 1, num_integ_points
               tau_tj(jquad) = x_tw(jquad)/2.0_dp
               tau_wj(jquad) = x_tw(jquad+num_integ_points)/2.0_dp
            END DO

            CALL get_l_sq_wghts_cos_tf_t_to_w(num_integ_points, tau_tj, weights_cos_tf_t_to_w, tj, &
                                              1.0_dp, Rc, max_error, 200)

            IF (iw > 0) THEN
               WRITE (iw, '(T2, I3, F12.1, ES12.3)') num_integ_points, Rc, max_error
            END IF

         END DO

         DEALLOCATE (x_tw, tau_tj, weights_cos_tf_t_to_w, tau_wj, wj, tj)

      END DO

   END SUBROUTINE test_least_square_ft

END MODULE rpa_ri_gpw
