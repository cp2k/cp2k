!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rountines to calculate RI-RPA energy
!> \par History
!>      06.2012 created [Mauro Del Ben]
!>      04.2015 GW routines added [Jan Wilhelm]
!>      10.2015 Cubic-scaling RPA routines added [Jan Wilhelm]
! **************************************************************************************************
MODULE rpa_ri_gpw
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE bibliography,                    ONLY: Bates2013,&
                                              DelBen2013,&
                                              DelBen2015,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_invert,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_syrk,&
                                              cp_fm_transpose,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_syevd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_indxg2l, cp_fm_indxg2p, &
        cp_fm_p_type, cp_fm_release, cp_fm_set_all, cp_fm_set_element, cp_fm_to_fm, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_allocate_matrix_set, dbcsr_copy, dbcsr_create, &
        dbcsr_deallocate_matrix_set, dbcsr_desymmetrize, dbcsr_filter, dbcsr_get_info, &
        dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_release_p, dbcsr_scale, dbcsr_set, dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_types,                       ONLY: hfx_release,&
                                              hfx_type
   USE input_constants,                 ONLY: ri_rpa_g0w0_crossing_newton,&
                                              ri_rpa_g0w0_crossing_none,&
                                              ri_rpa_g0w0_crossing_z_shot,&
                                              wfc_mm_style_gemm,&
                                              wfc_mm_style_syrk
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_sym_create,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_memory,&
                                              m_walltime
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE message_passing,                 ONLY: &
        mp_alltoall, mp_bcast, mp_comm_split_direct, mp_irecv, mp_isend, mp_min, mp_sendrecv, &
        mp_sum, mp_sync, mp_wait, mp_waitall
   USE minimax_exp,                     ONLY: check_exp_minimax_range,&
                                              get_exp_minimax_coeff
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff
   USE mp2_laplace,                     ONLY: laplace_minimax_approx
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type,&
                                              two_dim_int_array
   USE particle_types,                  ONLY: particle_type
   USE qs_band_structure,               ONLY: calculate_kp_orbitals
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE rpa_axk_util,                    ONLY: hfx_create_subgroup
   USE rpa_communication,               ONLY: fm_redistribute,&
                                              gamma_fm_to_dbcsr,&
                                              initialize_buffer,&
                                              release_buffer
   USE rpa_im_time,                     ONLY: communicate_buffer,&
                                              compute_mat_P_omega,&
                                              setup_mat_for_mem_cut_3c
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
   SUBROUTINE rpa_ri_compute_en(qs_env, Erpa, mp2_env, BIb_C, BIb_C_gw, para_env, para_env_sub, color_sub, &
                                ends_array, ends_B_virtual, ends_B_all, sizes_array, sizes_B_virtual, sizes_B_all, &
                                starts_array, starts_B_virtual, starts_B_all, mo_coeff, fm_matrix_L_RI_metric, &
                                Eigenval, nmo, homo, dimen_RI, gw_corr_lev_occ, gw_corr_lev_virt, &
                                unit_nr, do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_mao, matrix_s, &
                                mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                mat_munu, mat_dm_occ_local, mat_dm_virt_local, &
                                mat_P_local, mat_P_global, &
                                mat_M, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                                mat_3c_overl_int_mao_for_virt, eps_filter, BIb_C_beta, homo_beta, Eigenval_beta, &
                                ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, &
                                mo_coeff_beta, BIb_C_gw_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp)                                      :: Erpa
      TYPE(mp2_type), POINTER                            :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C, BIb_C_gw
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER                                            :: color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_virtual, ends_B_all, sizes_array, &
         sizes_B_virtual, sizes_B_all, starts_array, starts_B_virtual, starts_B_all
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, fm_matrix_L_RI_metric
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: nmo, homo, dimen_RI, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt, unit_nr
      LOGICAL                                            :: do_ri_sos_laplace_mp2, my_do_gw, &
                                                            do_im_time, do_mao
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, mao_coeff_occ, mao_coeff_virt, &
                                                            mao_coeff_occ_A, mao_coeff_virt_A
      TYPE(dbcsr_p_type)                                 :: mat_munu, mat_dm_occ_local, &
                                                            mat_dm_virt_local, mat_P_local, &
                                                            mat_P_global, mat_M
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      REAL(KIND=dp)                                      :: eps_filter
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: BIb_C_beta
      INTEGER, OPTIONAL                                  :: homo_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      INTEGER, ALLOCATABLE, DIMENSION(:), OPTIONAL       :: ends_B_virtual_beta, &
                                                            sizes_B_virtual_beta, &
                                                            starts_B_virtual_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mo_coeff_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: BIb_C_gw_beta
      INTEGER, OPTIONAL                                  :: gw_corr_lev_occ_beta, &
                                                            gw_corr_lev_virt_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_ri_compute_en', &
         routineP = moduleN//':'//routineN

      INTEGER :: best_integ_group_size, best_num_integ_point, color_rpa_group, crossing_search, &
         dimen_ia, dimen_ia_beta, dimen_nm_gw, handle, handle2, i, ierr, iiB, &
         input_integ_group_size, integ_group_size, jjB, max_iter_fit, min_integ_group_size, &
         my_B_all_end, my_B_all_size, my_B_all_start, my_B_size, my_B_size_beta, my_B_virtual_end, &
         my_B_virtual_end_beta, my_B_virtual_start, my_B_virtual_start_beta, my_group_L_end, &
         my_group_L_size, my_group_L_start, my_homo_beta, my_ia_end, my_ia_end_beta, my_ia_size, &
         my_ia_size_beta, my_ia_start, my_ia_start_beta, my_nm_gw_end, my_nm_gw_size, &
         my_nm_gw_start
      INTEGER :: ncol_block_mat, ngroup, nrow_block_mat, num_integ_group, num_integ_points, &
         num_poles, pos_integ_group, virtual, virtual_beta
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_ia, ends_ia_beta, ends_nm_gw, sizes_ia, &
         sizes_ia_beta, sizes_nm_gw, starts_ia, starts_ia_beta, starts_nm_gw, sub_proc_map
      LOGICAL                                            :: check_fit, do_minimax_quad, &
                                                            my_open_shell, skip_integ_group_opt
      REAL(KIND=dp) :: allowed_memory, avail_mem, E_Range, Emax, Emax_beta, Emin, Emin_beta, &
         ext_scaling, fermi_level_offset, mem_for_iaK, mem_for_QK, mem_min, mem_per_group, &
         mem_real, needed_mem, omega_max_fit, stop_crit
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D, BIb_C_2D_beta, BIb_C_2D_gw, &
                                                            BIb_C_2D_gw_beta
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_P_omega, fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(cp_fm_type), POINTER :: fm_mat_Q, fm_mat_Q_beta, fm_mat_Q_gemm, fm_mat_Q_gemm_beta, &
         fm_mat_R_gw, fm_mat_R_gw_beta, fm_mat_R_gw_gemm, fm_mat_R_gw_gemm_beta, fm_mat_S, &
         fm_mat_S_beta, fm_mat_S_gw, fm_mat_S_gw_beta, fm_mo_coeff_occ, fm_mo_coeff_occ_beta, &
         fm_mo_coeff_virt, fm_mo_coeff_virt_beta
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

      ! for imaginary time or periodic GW, we use all processors for a single frequency/time point
      IF (do_im_time .OR. mp2_env%ri_g0w0%do_periodic) THEN

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

         ! here we try to find the best rpa/lapleca group size
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

      ! reordering is not necessary for imaginary time
      IF (.NOT. do_im_time) THEN

         my_B_virtual_start = starts_B_virtual(para_env_sub%mepos)
         my_B_virtual_end = ends_B_virtual(para_env_sub%mepos)
         my_B_size = sizes_B_virtual(para_env_sub%mepos)

         IF (my_open_shell) THEN
            my_B_virtual_start_beta = starts_B_virtual_beta(para_env_sub%mepos)
            my_B_virtual_end_beta = ends_B_virtual_beta(para_env_sub%mepos)
            my_B_size_beta = sizes_B_virtual_beta(para_env_sub%mepos)
            my_homo_beta = homo_beta
         ELSE
            my_B_virtual_start_beta = my_B_virtual_start
            my_B_virtual_end_beta = my_B_virtual_end
            my_B_size_beta = my_B_size
            my_homo_beta = homo
         END IF

      END IF ! not imaginary time

      IF (my_do_gw) THEN
         my_B_all_start = starts_B_all(para_env_sub%mepos)
         my_B_all_end = ends_B_all(para_env_sub%mepos)
         my_B_all_size = sizes_B_all(para_env_sub%mepos)

         ext_scaling = mp2_env%ri_g0w0%scaling
         num_poles = mp2_env%ri_g0w0%num_poles
         omega_max_fit = mp2_env%ri_g0w0%omega_max_fit
         stop_crit = mp2_env%ri_g0w0%stop_crit
         max_iter_fit = mp2_env%ri_g0w0%max_iter_fit
         check_fit = mp2_env%ri_g0w0%check_fit
         crossing_search = mp2_env%ri_g0w0%crossing_search
         fermi_level_offset = mp2_env%ri_g0w0%fermi_level_offset
      END IF

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
                                 sizes_ia, starts_ia, ends_ia, &
                                 sizes_B_virtual, starts_B_virtual, ends_B_virtual, &
                                 sub_proc_map, my_ia_size, my_ia_start, my_ia_end, my_group_L_size, &
                                 my_B_size, my_B_virtual_start)
         DEALLOCATE (BIb_C)
         DEALLOCATE (starts_B_virtual)
         DEALLOCATE (ends_B_virtual)
         DEALLOCATE (sizes_B_virtual)

         ! The same for open shell
         IF (my_open_shell) THEN
            CALL calculate_BIb_C_2D(BIb_C_2D_beta, BIb_C_beta, para_env_sub, &
                                    dimen_ia_beta, homo_beta, virtual_beta, &
                                    sizes_ia_beta, starts_ia_beta, ends_ia_beta, &
                                    sizes_B_virtual_beta, starts_B_virtual_beta, ends_B_virtual_beta, &
                                    sub_proc_map, my_ia_size_beta, my_ia_start_beta, my_ia_end_beta, my_group_L_size, &
                                    my_B_size_beta, my_B_virtual_start_beta)

            DEALLOCATE (BIb_C_beta)
            DEALLOCATE (starts_B_virtual_beta)
            DEALLOCATE (ends_B_virtual_beta)
            DEALLOCATE (sizes_B_virtual_beta)

         END IF

      END IF

      ! in the GW case, BIb_C_2D_gw is an nm x K matrix, with n: number of corr GW levels, m=nmo
      IF (my_do_gw) THEN
         CALL calculate_BIb_C_2D(BIb_C_2D_gw, BIb_C_gw, para_env_sub, dimen_nm_gw, &
                                 gw_corr_lev_occ+gw_corr_lev_virt, nmo, &
                                 sizes_nm_gw, starts_nm_gw, ends_nm_gw, &
                                 sizes_B_all, starts_B_all, ends_B_all, &
                                 sub_proc_map, my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, my_group_L_size, &
                                 my_B_all_size, my_B_all_start)

         ! The same for open shell
         IF (my_open_shell) THEN
            DEALLOCATE (sizes_nm_gw, starts_nm_gw, ends_nm_gw)
            CALL calculate_BIb_C_2D(BIb_C_2D_gw_beta, BIb_C_gw_beta, para_env_sub, dimen_nm_gw, &
                                    gw_corr_lev_occ+gw_corr_lev_virt, nmo, &
                                    sizes_nm_gw, starts_nm_gw, ends_nm_gw, &
                                    sizes_B_all, starts_B_all, ends_B_all, &
                                    sub_proc_map, my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, my_group_L_size, &
                                    my_B_all_size, my_B_all_start)
            DEALLOCATE (BIb_C_gw_beta)
         END IF

         DEALLOCATE (BIb_C_gw)
         DEALLOCATE (starts_B_all)
         DEALLOCATE (ends_B_all)
         DEALLOCATE (sizes_B_all)

      END IF

      CALL timestop(handle2)

      ! now create the matrices needed for the calculation, Q, S and G
      ! Q and G will have omega dependence
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
                            nmo=nmo, homo=homo, num_integ_points=num_integ_points, &
                            fm_mat_P_omega=fm_mat_P_omega, do_mao=do_mao, &
                            mao_coeff_occ_A=mao_coeff_occ_A, mao_coeff_virt_A=mao_coeff_virt_A, &
                            matrix_s=matrix_s)

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
                                  .TRUE., fm_mat_Q%matrix_struct%context, fm_mat_Q%matrix_struct%context, .TRUE.)

         END IF

      END IF

      IF (my_open_shell) THEN

         IF (do_ri_sos_laplace_mp2) THEN
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
                             my_do_gw, gw_corr_lev_occ, gw_corr_lev_virt, num_poles, ext_scaling, omega_max_fit, &
                             stop_crit, check_fit, fermi_level_offset, crossing_search, &
                             max_iter_fit, mp2_env%ri_rpa%mm_style, do_minimax_quad, &
                             do_im_time, do_mao, fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                             fm_mat_P_omega, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                             mat_munu, mat_dm_occ_local, mat_dm_virt_local, mat_P_local, mat_P_global, mat_M, &
                             mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, matrix_s, &
                             mao_coeff_occ, mao_coeff_virt, eps_filter, &
                             fm_mo_coeff_occ_beta=fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta=fm_mo_coeff_virt_beta, &
                             homo_beta=homo_beta, virtual_beta=virtual_beta, &
                             dimen_ia_beta=dimen_ia_beta, Eigenval_beta=Eigenval_beta, fm_mat_S_beta=fm_mat_S_beta, &
                             fm_mat_Q_gemm_beta=fm_mat_Q_gemm_beta, fm_mat_Q_beta=fm_mat_Q_beta, &
                             fm_mat_S_gw_beta=fm_mat_S_gw_beta, gw_corr_lev_occ_beta=gw_corr_lev_occ_beta, &
                             gw_corr_lev_virt_beta=gw_corr_lev_virt_beta)
         END IF
      ELSE
         IF (do_ri_sos_laplace_mp2) THEN
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
                             my_do_gw, gw_corr_lev_occ, gw_corr_lev_virt, num_poles, ext_scaling, omega_max_fit, &
                             stop_crit, check_fit, fermi_level_offset, crossing_search, &
                             max_iter_fit, mp2_env%ri_rpa%mm_style, do_minimax_quad, &
                             do_im_time, do_mao, fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                             fm_mat_P_omega, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                             mat_munu, mat_dm_occ_local, mat_dm_virt_local, mat_P_local, mat_P_global, mat_M, &
                             mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, matrix_s, &
                             mao_coeff_occ, mao_coeff_virt, eps_filter)
         END IF
      END IF

      IF (.NOT. do_im_time) THEN

         DEALLOCATE (sizes_ia)
         DEALLOCATE (starts_ia)
         DEALLOCATE (ends_ia)

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
      IF (my_open_shell .AND. (.NOT. do_im_time)) THEN
         DEALLOCATE (sizes_ia_beta)
         DEALLOCATE (starts_ia_beta)
         DEALLOCATE (ends_ia_beta)
         CALL cp_fm_release(fm_mat_S_beta)
         CALL cp_fm_release(fm_mat_Q_gemm_beta)
         CALL cp_fm_release(fm_mat_Q_beta)
      END IF

      IF (my_do_gw) THEN
         DEALLOCATE (sizes_nm_gw)
         DEALLOCATE (starts_nm_gw)
         DEALLOCATE (ends_nm_gw)
         CALL cp_fm_release(fm_mat_S_gw)
         CALL cp_fm_release(fm_mat_R_gw)
         IF (my_open_shell) THEN
            CALL cp_fm_release(fm_mat_S_gw_beta)
         END IF
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
!>        Subroutine created to avoid massive double coding for GW
!> \param BIb_C_2D ...
!> \param BIb_C ...
!> \param para_env_sub ...
!> \param dimen_ia ...
!> \param homo ...
!> \param virtual ...
!> \param sizes_ia ...
!> \param starts_ia ...
!> \param ends_ia ...
!> \param sizes_B_virtual ...
!> \param starts_B_virtual ...
!> \param ends_B_virtual ...
!> \param sub_proc_map ...
!> \param my_ia_size ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \param my_B_size ...
!> \param my_B_virtual_start ...
!> \author Jan Wilhelm, 03/2015
! **************************************************************************************************
   SUBROUTINE calculate_BIb_C_2D(BIb_C_2D, BIb_C, para_env_sub, dimen_ia, homo, virtual, &
                                 sizes_ia, starts_ia, ends_ia, &
                                 sizes_B_virtual, starts_B_virtual, ends_B_virtual, &
                                 sub_proc_map, my_ia_size, my_ia_start, my_ia_end, my_group_L_size, &
                                 my_B_size, my_B_virtual_start)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER                                            :: dimen_ia, homo, virtual
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_ia, starts_ia, ends_ia, &
                                                            sizes_B_virtual, starts_B_virtual, &
                                                            ends_B_virtual, sub_proc_map
      INTEGER                                            :: my_ia_size, my_ia_start, my_ia_end, &
                                                            my_group_L_size, my_B_size, &
                                                            my_B_virtual_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_BIb_C_2D', &
         routineP = moduleN//':'//routineN
      INTEGER, PARAMETER                                 :: occ_chunk = 128

      INTEGER :: ia_global, iiB, iproc, itmp(2), jjB, occ_high, occ_low, proc_receive, proc_send, &
         proc_shift, rec_B_size, rec_B_virtual_end, rec_B_virtual_start
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
!> \param num_integ_points ...
!> \param fm_mat_P_omega ...
!> \param do_mao ...
!> \param mao_coeff_occ_A ...
!> \param mao_coeff_virt_A ...
!> \param matrix_s ...
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
                               fm_mo_coeff_occ, fm_mo_coeff_virt, nmo, homo, num_integ_points, &
                               fm_mat_P_omega, do_mao, mao_coeff_occ_A, mao_coeff_virt_A, &
                               matrix_s)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER :: color_sub, ngroup, integ_group_size, dimen_RI, dimen_ia, dimen_ia_for_block_size, &
         color_rpa_group, ext_row_block_size, ext_col_block_size, unit_nr, my_ia_size, &
         my_ia_start, my_ia_end, my_group_L_size, my_group_L_start, my_group_L_end
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S, fm_mat_Q_gemm, fm_mat_Q
      INTEGER                                            :: nrow_block_mat, ncol_block_mat
      LOGICAL, OPTIONAL                                  :: beta_case
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext, blacs_env_ext_S
      LOGICAL, OPTIONAL                                  :: do_gw_open_shell, do_im_time
      TYPE(cp_fm_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: mo_coeff, fm_mo_coeff_occ, &
                                                            fm_mo_coeff_virt
      INTEGER, OPTIONAL                                  :: nmo, homo, num_integ_points
      TYPE(cp_fm_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: fm_mat_P_omega
      LOGICAL, OPTIONAL                                  :: do_mao
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mao_coeff_occ_A, mao_coeff_virt_A, &
                                                            matrix_s

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

         ALLOCATE (fm_scaled_dm_occ_tau(num_integ_points))
         DO iii = 1, num_integ_points
            IF (my_do_mao) THEN
               CALL cp_fm_create(fm_scaled_dm_occ_tau(iii)%matrix, fm_struct)
            ELSE
               ! default: no MAOs
               CALL cp_fm_create(fm_scaled_dm_occ_tau(iii)%matrix, mo_coeff%matrix_struct)
            END IF
            CALL cp_fm_set_all(fm_scaled_dm_occ_tau(iii)%matrix, 0.0_dp)
         END DO

         ALLOCATE (fm_scaled_dm_virt_tau(num_integ_points))
         DO iii = 1, num_integ_points
            CALL cp_fm_create(fm_scaled_dm_virt_tau(iii)%matrix, mo_coeff%matrix_struct)
            CALL cp_fm_set_all(fm_scaled_dm_virt_tau(iii)%matrix, 0.0_dp)
         END DO

         CALL create_occ_virt_mo_coeffs(fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, &
                                        nmo, homo, my_do_mao, mao_coeff_occ_A, &
                                        mao_coeff_virt_A, matrix_s, 1)

         ALLOCATE (fm_mat_P_omega(num_integ_points))
         DO iii = 1, num_integ_points
            CALL cp_fm_create(fm_mat_P_omega(iii)%matrix, fm_mat_Q%matrix_struct)
            CALL cp_fm_set_all(fm_mat_P_omega(iii)%matrix, 0.0_dp)
         END DO

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
      INTEGER                                            :: nmo, homo
      LOGICAL                                            :: do_mao
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mao_coeff_occ_A, mao_coeff_virt_A, &
                                                            matrix_s
      INTEGER                                            :: ispin

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
      INTEGER                                            :: ispin, nmo
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
                        template=matrix_s(ispin)%matrix, &
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

   END SUBROUTINE

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
!> \param my_do_gw ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param num_poles ...
!> \param ext_scaling ...
!> \param omega_max_fit ...
!> \param stop_crit ...
!> \param check_fit ...
!> \param fermi_level_offset ...
!> \param crossing_search ...
!> \param max_iter_fit ...
!> \param mm_style ...
!> \param do_minimax_quad ...
!> \param do_im_time ...
!> \param do_mao ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param mo_coeff ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_mat_P_omega ...
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
! **************************************************************************************************
   SUBROUTINE rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                          homo, virtual, dimen_RI, dimen_ia, dimen_nm_gw, &
                          Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                          fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                          my_do_gw, gw_corr_lev_occ, gw_corr_lev_virt, num_poles, ext_scaling, &
                          omega_max_fit, stop_crit, check_fit, fermi_level_offset, crossing_search, &
                          max_iter_fit, mm_style, do_minimax_quad, do_im_time, do_mao, fm_mo_coeff_occ, &
                          fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, fm_mat_P_omega, &
                          fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, mat_munu, mat_dm_occ_local, &
                          mat_dm_virt_local, mat_P_local, &
                          mat_P_global, mat_M, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                          mat_3c_overl_int_mao_for_virt, matrix_s, mao_coeff_occ, mao_coeff_virt, &
                          eps_filter, fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, &
                          homo_beta, virtual_beta, dimen_ia_beta, Eigenval_beta, fm_mat_S_beta, &
                          fm_mat_Q_gemm_beta, fm_mat_Q_beta, fm_mat_S_gw_beta, &
                          gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp)                                      :: Erpa
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA, para_env_sub
      INTEGER                                            :: unit_nr, homo, virtual, dimen_RI, &
                                                            dimen_ia, dimen_nm_gw
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                                                            fm_mat_S_gw, fm_mat_R_gw
      LOGICAL                                            :: my_do_gw
      INTEGER                                            :: gw_corr_lev_occ, gw_corr_lev_virt, &
                                                            num_poles
      REAL(KIND=dp)                                      :: ext_scaling, omega_max_fit, stop_crit
      LOGICAL                                            :: check_fit
      REAL(KIND=dp)                                      :: fermi_level_offset
      INTEGER                                            :: crossing_search, max_iter_fit, mm_style
      LOGICAL                                            :: do_minimax_quad, do_im_time, do_mao
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            mo_coeff, fm_matrix_L_RI_metric
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_P_omega, fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(dbcsr_p_type)                                 :: mat_munu, mat_dm_occ_local, &
                                                            mat_dm_virt_local, mat_P_local, &
                                                            mat_P_global, mat_M
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
         mat_3c_overl_int_mao_for_virt, matrix_s, mao_coeff_occ, mao_coeff_virt
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mo_coeff_occ_beta, &
                                                            fm_mo_coeff_virt_beta
      INTEGER, OPTIONAL                                  :: homo_beta, virtual_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_beta, fm_mat_Q_gemm_beta, &
                                                            fm_mat_Q_beta, fm_mat_S_gw_beta
      INTEGER, OPTIONAL                                  :: gw_corr_lev_occ_beta, &
                                                            gw_corr_lev_virt_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_num_int', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit, re_unit
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda, Lambda_without_offset
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vec_Sigma_c_gw, vec_Sigma_c_gw_beta
      INTEGER :: avirt, col_start_local, color_sub_col, color_sub_row, contour_def_end, &
         contour_def_level_tot, contour_def_start, count_ev_sc_GW, cut_memory, cut_RI, &
         group_size_P, gw_corr_lev_tot, handle, handle2, handle3, handle4, i_cut_RI, i_ener, &
         i_global, i_mem, i_mo, ierr, iiB, index_contour_def, info_chol, iocc, iquad, iter_ev_sc, &
         j_global, j_mem, j_mo, jjB, jquad, LLL, m_global, m_global_beta, my_group_L_size_im_time, &
         my_num_dgemm_call, n_global, n_global_beta, n_group_col, n_group_row, n_level_gw, &
         n_level_gw_ref, n_local_col, n_local_row, nblkrows_total, ncol_local, ngroup_RI_orig
      INTEGER :: nm_global, nmo, nrow_local, num_fit_points, num_points_per_magnitude, &
         num_real_freq_cont_def, num_residues, num_var, number_of_rec, number_of_rec_axk, &
         number_of_rec_beta, number_of_send, number_of_send_axk, number_of_send_beta, row, &
         row_start_local
      INTEGER, ALLOCATABLE, DIMENSION(:) :: alpha_beta_from_jquad, map_rec_size, map_rec_size_axk, &
         map_rec_size_beta, map_send_size, map_send_size_axk, map_send_size_beta, &
         mepos_P_from_RI_row, mo_from_jquad, my_group_L_sizes_im_time, my_group_L_starts_im_time, &
         residue_from_jquad, row_from_LLL, RPA_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: ends_array_prim_col, ends_array_prim_fullcol, &
         ends_array_prim_fullrow, ends_array_prim_row, local_size_source, local_size_source_axk, &
         local_size_source_beta, sizes_array_prim_col, sizes_array_prim_row, &
         starts_array_prim_col, starts_array_prim_fullcol, starts_array_prim_fullrow, &
         starts_array_prim_row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: non_zero_blocks_3c_cut_col
      INTEGER, DIMENSION(:), POINTER :: col_blk_size, col_indices, ends_array_cm, &
         ends_array_cm_mao_occ, ends_array_cm_mao_virt, row_blk_offset, row_blk_size, row_indices, &
         starts_array_cm, starts_array_cm_mao_occ, starts_array_cm_mao_virt
      LOGICAL :: do_contour_def, do_mix_exchange, do_periodic, do_ri_Sigma_x, first_cycle, &
         first_cycle_periodic_correction, my_open_shell, remove_neg_virt_energies
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: multiply_needed_occ, multiply_needed_virt
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: cycle_due_to_sparse_dm
      REAL(KIND=dp) :: a_scaling, actual_flop_rate, alpha, contour_def_offset, delta_corr_nn, &
         delta_corr_nn_beta, e_axk, e_axk_corr, e_fermi, e_fermi_beta, E_Range, eigen_diff, Emax, &
         Emax_beta, Emin, Emin_beta, eps_filter_im_time, FComega, frac_rix, max_error_min, &
         my_flop_rate, omega, omega_i, omega_old, sign_occ_virt, t_end, t_start, tau, tau_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: delta_corr, Eigenval_last, Eigenval_last_beta, &
         Eigenval_scf, Eigenval_scf_beta, Lambda_Im, Lambda_Re, m_value, m_value_beta, Q_log, &
         real_freq, tau_tj, tau_wj, tj, trace_Qomega, vec_gw_energ, vec_gw_energ_beta, &
         vec_gw_energ_error_fit, vec_gw_energ_error_fit_beta, vec_omega_fit_gw, &
         vec_omega_fit_gw_sign, vec_omega_gw, vec_Sigma_c_gw_real_freq, &
         vec_Sigma_c_gw_real_freq_beta, vec_Sigma_x_gw, vec_Sigma_x_gw_beta, vec_W_gw, &
         vec_W_gw_beta, wj, x_tw, z_value, z_value_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_mat_M, delta_corr_nm, tau_wj_2d
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), POINTER :: fm_mat_L, fm_mat_L_transposed, fm_mat_S_contour_def, &
         fm_mat_S_contour_def_beta, fm_mat_S_gw_work, fm_mat_S_gw_work_beta, fm_mat_work
      TYPE(dbcsr_p_type)                                 :: mat_L, mat_M_P_munu_occ, &
                                                            mat_M_P_munu_virt, mat_P_global_copy
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, &
                                                            mat_P_omega, matrix_berry_im_mo_mo, &
                                                            matrix_berry_re_mo_mo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: mat_3c_overl_int_cut, &
         mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut, mat_dm_loc_occ_cut, &
         mat_dm_loc_virt_cut
      TYPE(dbcsr_type), POINTER                          :: mat_dm_loc_occ, mat_dm_loc_virt
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

      do_contour_def = .FALSE.

      ALLOCATE (tj(num_integ_points))
      tj = 0.0_dp

      ALLOCATE (wj(num_integ_points))
      wj = 0.0_dp

      ALLOCATE (Q_log(dimen_RI))

      ALLOCATE (trace_Qomega(dimen_RI))

      ! set up the least-square time grid and other matrices specifically for imag time
      IF (do_im_time) THEN

         CALL timeset(routineN//"_im_time_alloc_mat", handle4)

         ALLOCATE (tau_tj(num_integ_points))
         tau_tj = 0.0_dp

         ALLOCATE (tau_wj_2d(num_integ_points, num_integ_points))
         tau_wj_2d = 0.0_dp

         ALLOCATE (tau_wj(num_integ_points))
         tau_wj = 0.0_dp

         group_size_P = mp2_env%ri_rpa_im_time%group_size_P
         cut_memory = mp2_env%ri_rpa_im_time%cut_memory
         cut_RI = mp2_env%ri_rpa_im_time_util(1)%cut_RI
         eps_filter_im_time = mp2_env%ri_rpa_im_time%eps_filter_im_time
         ngroup_RI_orig = mp2_env%ri_rpa_im_time_util(1)%n_group_RI_orig

         ALLOCATE (my_group_L_starts_im_time(cut_RI))
         my_group_L_starts_im_time(:) = mp2_env%ri_rpa_im_time_util(1)%my_group_L_starts_im_time
         ALLOCATE (my_group_L_sizes_im_time(cut_RI))
         my_group_L_sizes_im_time(:) = mp2_env%ri_rpa_im_time_util(1)%my_group_L_sizes_im_time

         DO i_cut_RI = 1, cut_RI

            my_group_L_size_im_time = my_group_L_sizes_im_time(i_cut_RI)

            CALL dbcsr_filter(mat_3c_overl_int(i_cut_RI)%matrix, &
                              eps_filter*SQRT(REAL(my_group_L_size_im_time, KIND=dp)))

         END DO

         NULLIFY (mat_dm_loc_occ)
         CALL dbcsr_init_p(mat_dm_loc_occ)
         CALL dbcsr_desymmetrize(mat_dm_occ_local%matrix, mat_dm_loc_occ)

         CALL get_non_zero_blocks_3c(mat_3c_overl_int, para_env_sub, cut_RI, non_zero_blocks_3c)

         NULLIFY (mat_dm_loc_occ_cut)
         CALL dbcsr_allocate_matrix_set(mat_dm_loc_occ_cut, cut_RI, cut_memory)

         DO i_mem = 1, cut_memory
            DO i_cut_RI = 1, cut_RI

               ALLOCATE (mat_dm_loc_occ_cut(i_cut_RI, i_mem)%matrix)
               CALL dbcsr_create(matrix=mat_dm_loc_occ_cut(i_cut_RI, i_mem)%matrix, &
                                 !                                    template=mat_dm_loc)
                                 template=mat_dm_loc_occ)
            END DO
         END DO

         NULLIFY (mat_dm_loc_virt)
         CALL dbcsr_init_p(mat_dm_loc_virt)
         CALL dbcsr_desymmetrize(mat_dm_virt_local%matrix, mat_dm_loc_virt)

         NULLIFY (mat_dm_loc_virt_cut)
         CALL dbcsr_allocate_matrix_set(mat_dm_loc_virt_cut, cut_RI, cut_memory)

         DO i_mem = 1, cut_memory
            DO i_cut_RI = 1, cut_RI

               ALLOCATE (mat_dm_loc_virt_cut(i_cut_RI, i_mem)%matrix)
               CALL dbcsr_create(matrix=mat_dm_loc_virt_cut(i_cut_RI, i_mem)%matrix, &
                                 !                                    template=mat_dm_loc)
                                 template=mat_dm_loc_virt)

            END DO
         END DO

         CALL dbcsr_set(mat_munu%matrix, 0.0_dp)
         CALL dbcsr_filter(mat_munu%matrix, 1.0_dp)

         NULLIFY (mat_M_P_munu_occ%matrix)
         ALLOCATE (mat_M_P_munu_occ%matrix)
         CALL dbcsr_create(mat_M_P_munu_occ%matrix, template=mat_M%matrix)

         NULLIFY (mat_M_P_munu_virt%matrix)
         ALLOCATE (mat_M_P_munu_virt%matrix)
         CALL dbcsr_create(mat_M_P_munu_virt%matrix, template=mat_M%matrix)

         IF (.NOT. do_mao) THEN
            mat_3c_overl_int_mao_for_occ => mat_3c_overl_int
            mat_3c_overl_int_mao_for_virt => mat_3c_overl_int
         END IF

         NULLIFY (mat_M_mu_Pnu_occ)
         CALL dbcsr_allocate_matrix_set(mat_M_mu_Pnu_occ, cut_RI)
         DO i_cut_RI = 1, cut_RI
            ALLOCATE (mat_M_mu_Pnu_occ(i_cut_RI)%matrix)
            CALL dbcsr_create(matrix=mat_M_mu_Pnu_occ(i_cut_RI)%matrix, &
                              !                                 template=mat_3c_overl_int(i_cut_RI)%matrix)
                              template=mat_3c_overl_int_mao_for_occ(i_cut_RI)%matrix)
         END DO

         NULLIFY (mat_M_mu_Pnu_virt)
         CALL dbcsr_allocate_matrix_set(mat_M_mu_Pnu_virt, cut_RI)
         DO i_cut_RI = 1, cut_RI
            ALLOCATE (mat_M_mu_Pnu_virt(i_cut_RI)%matrix)
            CALL dbcsr_create(matrix=mat_M_mu_Pnu_virt(i_cut_RI)%matrix, &
                              !                                 template=mat_3c_overl_int(i_cut_RI)%matrix)
                              template=mat_3c_overl_int_mao_for_virt(i_cut_RI)%matrix)
         END DO

         NULLIFY (mat_P_omega)
         CALL dbcsr_allocate_matrix_set(mat_P_omega, num_integ_points)
         DO jquad = 1, num_integ_points
            ALLOCATE (mat_P_omega(jquad)%matrix)
            CALL dbcsr_create(matrix=mat_P_omega(jquad)%matrix, &
                              template=mat_P_global%matrix)
            CALL dbcsr_set(mat_P_omega(jquad)%matrix, 0.0_dp)
         END DO

         NULLIFY (mat_P_global_copy%matrix)
         ALLOCATE (mat_P_global_copy%matrix)
         CALL dbcsr_create(mat_P_global_copy%matrix, template=mat_P_global%matrix)
         CALL dbcsr_copy(mat_P_global_copy%matrix, mat_P_global%matrix)

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

!         ALLOCATE (starts_array_cm(cut_memory))
!         starts_array_cm(:) = mp2_env%ri_rpa_im_time%starts_array_cm(:)
!         ALLOCATE (ends_array_cm(cut_memory))
!         ends_array_cm(:) = mp2_env%ri_rpa_im_time%ends_array_cm(:)

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

         NULLIFY (fm_mat_work)
         CALL cp_fm_create(fm_mat_work, fm_mat_Q%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_Q, fm_mat_work)
         CALL cp_fm_set_all(matrix=fm_mat_work, alpha=0.0_dp)

         NULLIFY (fm_mat_L)
         CALL cp_fm_create(fm_mat_L, fm_mat_Q%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_Q, fm_mat_L)
         CALL cp_fm_set_all(matrix=fm_mat_L, alpha=0.0_dp)

         CALL setup_mat_for_mem_cut_3c(mat_3c_overl_int_cut, mat_3c_overl_int, cut_memory, cut_RI, &
                                       starts_array_cm, ends_array_cm, my_group_L_sizes_im_time, eps_filter)

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

         ALLOCATE (cycle_due_to_sparse_dm(cut_memory, cut_memory, num_integ_points))
         cycle_due_to_sparse_dm = .FALSE.

         ALLOCATE (multiply_needed_occ(cut_memory, cut_RI))
         multiply_needed_occ = .TRUE.

         ALLOCATE (multiply_needed_virt(cut_memory, cut_RI))
         multiply_needed_virt = .TRUE.

         ALLOCATE (row_from_LLL(dimen_RI))
         row_from_LLL = 0

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

         NULLIFY (fm_mat_L_transposed)
         CALL cp_fm_create(fm_mat_L_transposed, fm_mat_L%matrix_struct)
         CALL cp_fm_to_fm(fm_mat_L, fm_mat_L_transposed)
         CALL cp_fm_set_all(matrix=fm_mat_L_transposed, alpha=0.0_dp)

         CALL cp_fm_get_info(fm_mat_L_transposed, context=blacs_env)

         CALL cp_fm_copy_general(fm_matrix_L_RI_metric, fm_mat_L_transposed, blacs_env%para_env)

         CALL cp_fm_transpose(fm_mat_L_transposed, fm_mat_L)

         CALL cp_fm_release(fm_mat_L_transposed)

         NULLIFY (mat_L%matrix)
         ALLOCATE (mat_L%matrix)
         CALL dbcsr_create(mat_L%matrix, template=mat_P_global%matrix)

         CALL copy_fm_to_dbcsr(fm_mat_L, mat_L%matrix)

         CALL timestop(handle4)
      END IF

      IF (do_minimax_quad) THEN
         ! MINIMAX quadrature
         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp

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

         ierr = 0
         CALL get_rpa_minimax_coeff(num_integ_points, E_Range, x_tw, ierr)

         DO jquad = 1, num_integ_points
            tj(jquad) = x_tw(jquad)
            wj(jquad) = x_tw(jquad+num_integ_points)
         END DO

         DEALLOCATE (x_tw)

         IF (unit_nr > 0) THEN
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

         ! set up the minimax time grid
         IF (do_im_time) THEN

            ALLOCATE (x_tw(2*num_integ_points))
            x_tw = 0.0_dp

            CALL get_exp_minimax_coeff(num_integ_points, E_Range, x_tw)

            DO jquad = 1, num_integ_points
               tau_tj(jquad) = x_tw(jquad)/2.0_dp
               tau_wj(jquad) = x_tw(jquad+num_integ_points)/2.0_dp
            END DO

            DEALLOCATE (x_tw)

            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "INTEG_INFO| Range for the minimax approximation:", E_Range
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

            num_points_per_magnitude = mp2_env%ri_rpa_im_time%num_points_per_magnitude
            CALL get_rpa_least_sq_weights_tau(num_integ_points, tau_tj, tau_wj_2d, tj, &
                                              Emin, Emax, max_error_min, num_points_per_magnitude)

            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
                  "INTEG_INFO| Maximum deviation of the imag. time fit:", max_error_min
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
         IF (my_do_gw .AND. ext_scaling > 0.0_dp) THEN
            a_scaling = ext_scaling
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

      IF (my_open_shell) THEN
         alpha = 2.0_dp
      ELSE
         alpha = 4.0_dp
      END IF

      nmo = homo+virtual

      IF (my_do_gw) THEN
         gw_corr_lev_tot = gw_corr_lev_occ+gw_corr_lev_virt

         ALLOCATE (Eigenval_scf(nmo))
         Eigenval_scf(:) = Eigenval(:)

         ALLOCATE (Eigenval_last(nmo))
         Eigenval_last(:) = Eigenval(:)

         ! in the case of HF_diag approach of X. Blase (PRB 83, 115103 (2011), Sec. IV), we subtract the
         ! XC potential and add exact exchange
         IF (mp2_env%ri_g0w0%hf_like_ev_start) THEN
            DO n_level_gw = 1, gw_corr_lev_tot
               n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
               Eigenval(n_level_gw_ref) = Eigenval(n_level_gw_ref)+ &
                                          mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 1)
            END DO
         END IF

         ! Eigenval for beta
         IF (my_open_shell) THEN
            ALLOCATE (Eigenval_scf_beta(nmo))
            Eigenval_scf_beta(:) = Eigenval_beta(:)

            ALLOCATE (Eigenval_last_beta(nmo))
            Eigenval_last_beta(:) = Eigenval_beta(:)

            ! in the case of HF_diag approach of X. Blase (PRB 83, 115103 (2011), Sec. IV), we subtract the
            ! XC potential and add exact exchange
            IF (mp2_env%ri_g0w0%hf_like_ev_start) THEN
               DO n_level_gw = 1, gw_corr_lev_tot
                  n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
                  Eigenval_beta(n_level_gw_ref) = Eigenval_beta(n_level_gw_ref)+ &
                                                  mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 2)
               END DO
            END IF
         END IF

         ! create temporary matrix to store B*([1+Q(iw')]^-1-1), has the same size as B
         NULLIFY (fm_mat_S_gw_work)
         CALL cp_fm_create(fm_mat_S_gw_work, fm_mat_S_gw%matrix_struct)
         CALL cp_fm_set_all(matrix=fm_mat_S_gw_work, alpha=0.0_dp)

         IF (my_open_shell) THEN
            NULLIFY (fm_mat_S_gw_work_beta)
            CALL cp_fm_create(fm_mat_S_gw_work_beta, fm_mat_S_gw%matrix_struct)
            CALL cp_fm_set_all(matrix=fm_mat_S_gw_work_beta, alpha=0.0_dp)
         END IF

         ALLOCATE (vec_W_gw(dimen_nm_gw))
         vec_W_gw = 0.0_dp

         IF (my_open_shell) THEN
            ALLOCATE (vec_W_gw_beta(dimen_nm_gw))
            vec_W_gw_beta = 0.0_dp
         END IF

         im_unit = (0.0_dp, 1.0_dp)
         re_unit = (1.0_dp, 0.0_dp)

         ! fill the omega_frequency vector
         ALLOCATE (vec_omega_gw(num_integ_points))
         vec_omega_gw = 0.0_dp

         DO jquad = 1, num_integ_points
            IF (do_minimax_quad) THEN
               omega = tj(jquad)
            ELSE
               omega = a_scaling/TAN(tj(jquad))
            END IF
            vec_omega_gw(jquad) = omega
         END DO

         ! determine number of fit points in the interval [0,w_max] for virt, or [-w_max,0] for occ
         num_fit_points = 0

         DO jquad = 1, num_integ_points
            IF (vec_omega_gw(jquad) < omega_max_fit) THEN
               num_fit_points = num_fit_points+1
            END IF
         END DO

         ! create new arrays containing omega values at which we calculate vec_Sigma_c_gw
         ALLOCATE (vec_omega_fit_gw(num_fit_points))

         ALLOCATE (vec_omega_fit_gw_sign(num_fit_points))

         ALLOCATE (vec_Sigma_c_gw(gw_corr_lev_tot, num_fit_points))

         vec_Sigma_c_gw = (0.0_dp, 0.0_dp)

         IF (my_open_shell) THEN
            ALLOCATE (vec_Sigma_c_gw_beta(gw_corr_lev_tot, num_fit_points))
            vec_Sigma_c_gw_beta = (0.0_dp, 0.0_dp)
         END IF

         do_periodic = mp2_env%ri_g0w0%do_periodic

         IF (do_periodic) THEN

            ALLOCATE (delta_corr(1+homo-gw_corr_lev_occ:homo+gw_corr_lev_virt))
            delta_corr(:) = 0.0_dp

            ALLOCATE (delta_corr_nm(1+homo-gw_corr_lev_occ:homo+gw_corr_lev_virt, 1:nmo))
            delta_corr_nm(:, :) = 0.0_dp

            first_cycle_periodic_correction = .TRUE.

         END IF

         do_ri_Sigma_x = mp2_env%ri_g0w0%do_ri_Sigma_x

         do_mix_exchange = mp2_env%ri_g0w0%mix_exchange

         ! in case we do RI for Sigma_x, we calculate Sigma_x right here
         IF (do_ri_Sigma_x .OR. do_mix_exchange) THEN

            IF (do_mix_exchange) THEN
               frac_rix = 1.0_dp-mp2_env%ri_g0w0%frac_exx
            ELSE
               frac_rix = 1.0_dp
            END IF

            ALLOCATE (vec_Sigma_x_gw(nmo))
            vec_Sigma_x_gw = 0.0_dp

            IF (my_open_shell) THEN
               ALLOCATE (vec_Sigma_x_gw_beta(nmo))
               vec_Sigma_x_gw_beta = 0.0_dp
            END IF

            CALL cp_fm_get_info(matrix=fm_mat_S_gw, &
                                nrow_local=nrow_local, &
                                ncol_local=ncol_local, &
                                row_indices=row_indices, &
                                col_indices=col_indices)

            CALL mp_sync(para_env%group)

            ! loop over (nm) index
            DO iiB = 1, nrow_local

               ! this is needed for correct values within parallelization, see analogue
               ! statement below
               IF (MODULO(1, num_integ_group) /= color_rpa_group) CYCLE

               nm_global = row_indices(iiB)

               ! transform the index nm to n and m, formulae copied from Mauro's code
               n_global = MAX(1, nm_global-1)/nmo+1
               m_global = nm_global-(n_global-1)*nmo
               n_global = n_global+homo-gw_corr_lev_occ

               IF (my_open_shell) THEN
                  n_global_beta = MAX(1, nm_global-1)/nmo+1
                  m_global_beta = nm_global-(n_global_beta-1)*nmo
                  n_global_beta = n_global_beta+homo_beta-gw_corr_lev_occ_beta
               END IF

               IF (m_global <= homo) THEN

                  ! loop over auxiliary basis functions
                  DO jjB = 1, ncol_local

                     ! Sigma_x_n = -sum_m^occ sum_P (B_(nm)^P)^2
                     vec_Sigma_x_gw(n_global) = &
                        vec_Sigma_x_gw(n_global)- &
                        fm_mat_S_gw%local_data(iiB, jjB)**2

                  END DO

               END IF

               ! The same for beta
               IF (my_open_shell) THEN
                  IF (m_global_beta <= homo_beta) THEN
                     ! loop over auxiliary basis functions
                     DO jjB = 1, ncol_local

                        ! Sigma_x_n = -sum_m^occ sum_P (B_(nm)^P)^2
                        vec_Sigma_x_gw_beta(n_global_beta) = &
                           vec_Sigma_x_gw_beta(n_global_beta)- &
                           fm_mat_S_gw_beta%local_data(iiB, jjB)**2

                     END DO
                  END IF
               END IF
            END DO

            CALL mp_sum(vec_Sigma_x_gw, para_env%group)

            mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1) = &
               mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1)+ &
               frac_rix*vec_Sigma_x_gw(:)

            IF (my_open_shell) THEN

               CALL mp_sum(vec_Sigma_x_gw_beta, para_env%group)

               mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2) = &
                  mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2)+ &
                  frac_rix*vec_Sigma_x_gw_beta(:)

            END IF

         END IF

         ! fill the omega vector with frequencies, where we calculate the self-energy
         iquad = 0
         DO jquad = 1, num_integ_points
            IF (vec_omega_gw(jquad) < omega_max_fit) THEN
               iquad = iquad+1
               vec_omega_fit_gw(iquad) = vec_omega_gw(jquad)
            END IF
         END DO

         DEALLOCATE (vec_omega_gw)

         ! check if we want to do the contour deformation technique instead of analytical continuation
         contour_def_start = mp2_env%ri_g0w0%contour_def_start
         contour_def_end = mp2_env%ri_g0w0%contour_def_end

         IF (contour_def_start > 0 .AND. contour_def_end >= contour_def_start) THEN
            do_contour_def = .TRUE.
         END IF
         contour_def_offset = mp2_env%ri_g0w0%contour_def_offset

         IF (do_contour_def) THEN

            ! check range
            IF (contour_def_start < homo-gw_corr_lev_occ+1) THEN
               contour_def_start = homo-gw_corr_lev_occ+1
            END IF

            IF (contour_def_end > homo+gw_corr_lev_virt) THEN
               contour_def_end = homo+gw_corr_lev_virt
            END IF

            contour_def_level_tot = contour_def_end-contour_def_start+1

            num_real_freq_cont_def = 2*contour_def_level_tot

            ! set up array for the self-energy as function of real energy. one energy is slightly
            ! above and one below the DFT MO energy => the space is 2* number of levels which are
            ! treated by contour deformation technique
            ALLOCATE (vec_Sigma_c_gw_real_freq(num_real_freq_cont_def))

            IF (my_open_shell) THEN
               ALLOCATE (vec_Sigma_c_gw_real_freq_beta(num_real_freq_cont_def))
            END IF

            NULLIFY (fm_mat_S_contour_def)
            CALL cp_fm_create(fm_mat_S_contour_def, fm_mat_S%matrix_struct)
            CALL cp_fm_to_fm(source=fm_mat_S, destination=fm_mat_S_contour_def)

            IF (my_open_shell) THEN
               NULLIFY (fm_mat_S_contour_def_beta)
               CALL cp_fm_create(fm_mat_S_contour_def_beta, fm_mat_S_beta%matrix_struct)
               CALL cp_fm_to_fm(source=fm_mat_S_beta, destination=fm_mat_S_contour_def_beta)
            END IF

         END IF

         ! arrays storing the complex fit parameters a_0, a_1, b_1, a_2, b_2, ...
         num_var = 2*num_poles+1
         ALLOCATE (Lambda(num_var))
         Lambda = (0.0_dp, 0.0_dp)
         ALLOCATE (Lambda_without_offset(num_var))
         Lambda_without_offset = (0.0_dp, 0.0_dp)
         ALLOCATE (Lambda_Re(num_var))
         Lambda_Re = 0.0_dp
         ALLOCATE (Lambda_Im(num_var))
         Lambda_Im = 0.0_dp

         ! arrays storing the correlation self-energy, stat. error and z-shot value
         ALLOCATE (vec_gw_energ(gw_corr_lev_tot))
         vec_gw_energ = 0.0_dp
         ALLOCATE (vec_gw_energ_error_fit(gw_corr_lev_tot))
         vec_gw_energ_error_fit = 0.0_dp
         ALLOCATE (z_value(gw_corr_lev_tot))
         z_value = 0.0_dp
         ALLOCATE (m_value(gw_corr_lev_tot))
         m_value = 0.0_dp

         ! the same for beta
         IF (my_open_shell) THEN
            ALLOCATE (vec_gw_energ_beta(gw_corr_lev_tot))
            vec_gw_energ_beta = 0.0_dp
            ALLOCATE (vec_gw_energ_error_fit_beta(gw_corr_lev_tot))
            vec_gw_energ_error_fit_beta = 0.0_dp
            ALLOCATE (z_value_beta(gw_corr_lev_tot))
            z_value_beta = 0.0_dp
            ALLOCATE (m_value_beta(gw_corr_lev_tot))
            m_value_beta = 0.0_dp
         END IF

      END IF

      Erpa = 0.0_dp
      IF (mp2_env%ri_rpa%do_ri_axk) e_axk = 0.0_dp
      first_cycle = .TRUE.
      omega_old = 0.0_dp
      my_num_dgemm_call = 0
      my_flop_rate = 0.0_dp

      IF (my_do_gw) THEN
         iter_ev_sc = mp2_env%ri_g0w0%iter_ev_sc
      ELSE
         iter_ev_sc = 1
      END IF

      DO count_ev_sc_GW = 1, iter_ev_sc

         ! reset some values, important when doing eigenvalue self-consistent GW
         IF (my_do_gw) THEN
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

            e_fermi = (Eigenval(homo)+Eigenval(homo+1))*0.5_dp
            IF (my_open_shell) THEN
               e_fermi_beta = (Eigenval_beta(homo_beta)+Eigenval_beta(homo_beta+1))*0.5_dp
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

            ! imag. time RPA only with Minimax
            CPASSERT(do_minimax_quad)

            ! compute the matrix Q(it) and Fourier transform it directly to fm_mat_P_omega(iw)
            CALL compute_mat_P_omega(mat_P_omega, fm_scaled_dm_occ_tau, &
                                     fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                     mat_P_local, mat_P_global, mat_P_global_copy, &
                                     mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, matrix_s, &
                                     mao_coeff_occ, mao_coeff_virt, 1, &
                                     mat_M_P_munu_occ, mat_M_P_munu_virt, mat_3c_overl_int_cut, &
                                     mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut, &
                                     mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, &
                                     tau_wj_2d, tj, tau_tj, e_fermi, eps_filter, alpha, &
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
                                     do_mao)

            ! the same for open shell, use fm_mo_coeff_occ_beta and fm_mo_coeff_virt_beta
            IF (my_open_shell) THEN

               CALL compute_mat_P_omega(mat_P_omega, fm_scaled_dm_occ_tau, &
                                        fm_scaled_dm_virt_tau, fm_mo_coeff_occ_beta, fm_mo_coeff_virt_beta, &
                                        mat_P_local, mat_P_global, mat_P_global_copy, &
                                        mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, matrix_s, &
                                        mao_coeff_occ, mao_coeff_virt, 2, &
                                        mat_M_P_munu_occ, mat_M_P_munu_virt, mat_3c_overl_int_cut, &
                                        mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut, &
                                        mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, &
                                        tau_wj_2d, tj, tau_tj, e_fermi_beta, eps_filter, alpha, &
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
                                        do_mao)

            END IF

         END IF ! do im time

         num_residues = 0

         ! get the energies for evaluating the residues
         IF (do_contour_def) THEN
            vec_Sigma_c_gw_real_freq = 0.0_dp
            IF (my_open_shell) THEN
               vec_Sigma_c_gw_real_freq_beta = 0.0_dp
            END IF

            ! determine how many residues we need for the contour deformation technique
            DO i_mo = contour_def_start, contour_def_end
               IF (i_mo <= homo) THEN
                  DO j_mo = 1, homo
                     IF (Eigenval(i_mo)-contour_def_offset < Eigenval(j_mo)) THEN
                        num_residues = num_residues+1
                     END IF
                  END DO
               ELSE
                  DO j_mo = homo+1, nmo
                     IF (Eigenval(i_mo)+contour_def_offset > Eigenval(j_mo)) THEN
                        num_residues = num_residues+1
                     END IF
                  END DO
               END IF
            END DO

            ! add residues for beta
            IF (my_open_shell) THEN
               DO i_mo = contour_def_start, contour_def_end
                  IF (i_mo <= homo_beta) THEN
                     DO j_mo = 1, homo_beta
                        IF (Eigenval_beta(i_mo)-contour_def_offset < Eigenval_beta(j_mo)) THEN
                           num_residues = num_residues+1
                        END IF
                     END DO
                  ELSE
                     DO j_mo = homo_beta+1, nmo
                        IF (Eigenval_beta(i_mo)+contour_def_offset > Eigenval_beta(j_mo)) THEN
                           num_residues = num_residues+1
                        END IF
                     END DO
                  END IF
               END DO
            END IF

            ! allocate real frequencies of the residues
            ALLOCATE (real_freq(num_residues))

            ! for getting the MO to which the residue is beloning
            ALLOCATE (mo_from_jquad(num_residues))

            ! for getting the residue (connected with the MO m)
            ALLOCATE (residue_from_jquad(num_residues))

            ! for spin information
            ALLOCATE (alpha_beta_from_jquad(num_residues))

            i_ener = 0
            ! determine which residue number belongs to which MO, residue and energy
            DO i_mo = contour_def_start, contour_def_end
               IF (i_mo <= homo) THEN
                  DO j_mo = 1, homo
                     IF (Eigenval(i_mo)-contour_def_offset < Eigenval(j_mo)) THEN
                        i_ener = i_ener+1
                        mo_from_jquad(i_ener) = i_mo
                        residue_from_jquad(i_ener) = j_mo
                        real_freq(i_ener) = Eigenval(j_mo)-Eigenval(i_mo)
                        alpha_beta_from_jquad(i_ener) = 1
                     END IF
                  END DO
               ELSE
                  DO j_mo = homo+1, nmo
                     IF (Eigenval(i_mo)+contour_def_offset > Eigenval(j_mo)) THEN
                        i_ener = i_ener+1
                        mo_from_jquad(i_ener) = i_mo
                        residue_from_jquad(i_ener) = j_mo
                        real_freq(i_ener) = Eigenval(j_mo)-Eigenval(i_mo)
                        alpha_beta_from_jquad(i_ener) = 1
                     END IF
                  END DO
               END IF
            END DO

            ! residues for beta
            IF (my_open_shell) THEN

               DO i_mo = contour_def_start, contour_def_end
                  IF (i_mo <= homo_beta) THEN
                     DO j_mo = 1, homo_beta
                        IF (Eigenval_beta(i_mo)-contour_def_offset < Eigenval_beta(j_mo)) THEN
                           i_ener = i_ener+1
                           mo_from_jquad(i_ener) = i_mo
                           residue_from_jquad(i_ener) = j_mo
                           real_freq(i_ener) = Eigenval_beta(j_mo)-Eigenval_beta(i_mo)
                           alpha_beta_from_jquad(i_ener) = 2
                        END IF
                     END DO
                  ELSE
                     DO j_mo = homo_beta+1, nmo
                        IF (Eigenval_beta(i_mo)+contour_def_offset > Eigenval_beta(j_mo)) THEN
                           i_ener = i_ener+1
                           mo_from_jquad(i_ener) = i_mo
                           residue_from_jquad(i_ener) = j_mo
                           real_freq(i_ener) = Eigenval_beta(j_mo)-Eigenval_beta(i_mo)
                           alpha_beta_from_jquad(i_ener) = 2
                        END IF
                     END DO
                  END IF
               END DO

            END IF

         END IF

         DO jquad = 1, num_integ_points+num_residues

            IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE

            CALL timeset(routineN//"_RPA_matrix_operations", handle3)

            IF (jquad <= num_integ_points) THEN

               IF (do_minimax_quad) THEN
                  omega = tj(jquad)
               ELSE
                  omega = a_scaling/TAN(tj(jquad))
               END IF

            ELSE
               ! we do contour deformation for GW, where we have a real frequency
               omega = real_freq(jquad-num_integ_points)

            END IF

            IF (do_im_time) THEN

               ! in case we do imag time, we already calculated fm_mat_Q by a Fourier transform from im. time

               ! multiplication with RI metric/Coulomb operator
               CALL dbcsr_multiply("N", "T", 1.0_dp, mat_P_omega(jquad)%matrix, mat_L%matrix, &
                                   0.0_dp, mat_P_global%matrix, filter_eps=eps_filter_im_time)

               CALL copy_dbcsr_to_fm(mat_P_global%matrix, fm_mat_work)

               CALL cp_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, 1.0_dp, fm_mat_L, fm_mat_work, &
                            0.0_dp, fm_mat_Q)

            ELSE

               ! get info of fm_mat_S
               CALL cp_fm_get_info(matrix=fm_mat_S, &
                                   nrow_local=nrow_local, &
                                   ncol_local=ncol_local, &
                                   row_indices=row_indices, &
                                   col_indices=col_indices)

               ! remove eigenvalue part of S matrix from the last eigenvalue self-c. cycle
               IF (first_cycle .AND. count_ev_sc_GW > 1) THEN
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                         SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval_last,fm_mat_S,virtual,homo,omega_old,&
!$OMP                                do_contour_def,fm_mat_S_contour_def)
                  DO jjB = 1, ncol_local
                     j_global = col_indices(jjB)
                     DO iiB = 1, nrow_local
                        i_global = row_indices(iiB)

                        iocc = MAX(1, i_global-1)/virtual+1
                        avirt = i_global-(iocc-1)*virtual
                        eigen_diff = Eigenval_last(avirt+homo)-Eigenval_last(iocc)

                        ! for the contour deformation, we have the original B-matrix
                        IF (do_contour_def) THEN
                           fm_mat_S%local_data(iiB, jjB) = fm_mat_S_contour_def%local_data(iiB, jjB)
                        ELSE
                           fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)/ &
                                                           SQRT(eigen_diff/(eigen_diff**2+omega_old**2))
                        END IF

                     END DO
                  END DO

               END IF

               ! update G matrix with the new value of omega
               IF (first_cycle) THEN
                  ! In this case just update the matrix (symmetric form) with
                  ! SQRT((epsi_a-epsi_i)/((epsi_a-epsi_i)**2+omega**2))
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                         SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,fm_mat_S,virtual,homo,omega)
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
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                         SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval,&
!$OMP                                fm_mat_S,virtual,homo,omega,omega_old,real_freq,jquad,&
!$OMP                                num_integ_points,fm_mat_S_contour_def)
                  DO jjB = 1, ncol_local
                     j_global = col_indices(jjB)
                     DO iiB = 1, nrow_local
                        i_global = row_indices(iiB)

                        iocc = MAX(1, i_global-1)/virtual+1
                        avirt = i_global-(iocc-1)*virtual
                        eigen_diff = Eigenval(avirt+homo)-Eigenval(iocc)

                        IF (jquad <= num_integ_points) THEN

                           fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)* &
                                                           SQRT((eigen_diff**2+omega_old**2)/(eigen_diff**2+omega**2))

                        ELSE

                           ! reset
                           fm_mat_S%local_data(iiB, jjB) = fm_mat_S_contour_def%local_data(iiB, jjB)

                           ! update
                           fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)*eigen_diff/(eigen_diff**2-omega**2)

                        END IF
                     END DO
                  END DO
               END IF

               ! alpha=4.0 is valid only for closed shell systems (alpha is 2.0 in the spin-orbital basis)
               t_start = m_walltime()
               SELECT CASE (mm_style)
               CASE (wfc_mm_style_gemm)
                  ! waste-fully computes the full symmetrix matrix, but maybe faster than cp_fm_syrk for optimized cp_fm_gemm
                  IF (jquad <= num_integ_points) THEN
                     ! RPA for imaginary frequencies
                     CALL cp_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia, alpha=alpha, &
                                  matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                                  matrix_c=fm_mat_Q_gemm)
                  ELSE
                     ! GW contour deformation for real frequencies with fm_mat_S_contour_def
                     CALL cp_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia, alpha=alpha, &
                                  matrix_a=fm_mat_S, matrix_b=fm_mat_S_contour_def, beta=0.0_dp, &
                                  matrix_c=fm_mat_Q_gemm)
                  END IF
               CASE (wfc_mm_style_syrk)
                  ! will only compute the upper half of the matrix, which is fine, since we only use it for cholesky later
                  CALL cp_fm_syrk(uplo='U', trans='T', k=dimen_ia, alpha=alpha, matrix_a=fm_mat_S, &
                                  ia=1, ja=1, beta=0.0_dp, matrix_c=fm_mat_Q_gemm)
                  IF (jquad > num_integ_points) THEN
                     CPABORT("GW contour deformation does not work with wfc_mm_style_syrk.")
                  END IF
               CASE DEFAULT
                  CPABORT("")
               END SELECT
               t_end = m_walltime()

               actual_flop_rate = 2.0_dp*REAL(dimen_ia, KIND=dp)*dimen_RI*REAL(dimen_RI, KIND=dp)
               actual_flop_rate = actual_flop_rate/(MAX(0.01_dp, t_end-t_start))
               IF (para_env_RPA%mepos == 0) my_flop_rate = my_flop_rate+actual_flop_rate
               my_num_dgemm_call = my_num_dgemm_call+1

               ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
               CALL cp_fm_set_all(matrix=fm_mat_Q, alpha=0.0_dp)
               CALL fm_redistribute(fm_mat_Q_gemm, fm_mat_Q, RPA_proc_map, buffer_rec, buffer_send, &
                                    number_of_send, &
                                    map_send_size, map_rec_size, local_size_source, para_env_RPA)

               CALL cp_fm_get_info(matrix=fm_mat_Q, &
                                   nrow_local=nrow_local, &
                                   ncol_local=ncol_local)

               IF (my_open_shell) THEN
                  ! the same for the beta spin
                  !XXX ! copy fm_mat_S into fm_mat_G
                  !XXX CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta)
                  !XXX ! get info of fm_mat_G_beta
                  !XXX CALL cp_fm_get_info(matrix=fm_mat_G_beta,&
                  !XXX                     nrow_local=nrow_local,&
                  !XXX                     ncol_local=ncol_local,&
                  !XXX                     row_indices=row_indices,&
                  !XXX                     col_indices=col_indices)
                  !XXX

                  ! get info of fm_mat_S_beta
                  CALL cp_fm_get_info(matrix=fm_mat_S_beta, &
                                      nrow_local=nrow_local, &
                                      ncol_local=ncol_local, &
                                      row_indices=row_indices, &
                                      col_indices=col_indices)

                  ! remove eigenvalue part of S matrix from the last eigenvalue self-c. cycle
                  IF (first_cycle .AND. count_ev_sc_GW > 1) THEN
!$OMP               PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                           SHARED(ncol_local,nrow_local,col_indices,row_indices,Eigenval_last_beta,&
!$OMP                                  fm_mat_S_beta,virtual_beta,homo_beta,omega_old,do_contour_def,&
!$OMP                                  fm_mat_S_contour_def_beta)
                     DO jjB = 1, ncol_local
                        j_global = col_indices(jjB)
                        DO iiB = 1, nrow_local
                           i_global = row_indices(iiB)

                           iocc = MAX(1, i_global-1)/virtual_beta+1
                           avirt = i_global-(iocc-1)*virtual_beta
                           eigen_diff = Eigenval_last_beta(avirt+homo_beta)-Eigenval_last_beta(iocc)

                           ! for the contour deformation, the last omega was a residue, so negative sign
                           IF (do_contour_def) THEN
                              fm_mat_S_beta%local_data(iiB, jjB) = fm_mat_S_contour_def_beta%local_data(iiB, jjB)
                           ELSE
                              fm_mat_S_beta%local_data(iiB, jjB) = fm_mat_S_beta%local_data(iiB, jjB)/ &
                                                                   SQRT(eigen_diff/(eigen_diff**2+omega_old**2))
                           END IF

                        END DO

                     END DO

                  END IF

                  ! update G matrix with the new value of omega
                  IF (first_cycle) THEN
!$OMP               PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                           SHARED(ncol_local,nrow_local,col_indices,row_indices,&
!$OMP                                  Eigenval_beta,fm_mat_S_beta,virtual_beta,homo_beta,omega)
                     DO jjB = 1, ncol_local
                        j_global = col_indices(jjB)
                        DO iiB = 1, nrow_local
                           i_global = row_indices(iiB)

                           iocc = MAX(1, i_global-1)/virtual_beta+1
                           avirt = i_global-(iocc-1)*virtual_beta
                           eigen_diff = Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

                           fm_mat_S_beta%local_data(iiB, jjB) = fm_mat_S_beta%local_data(iiB, jjB)* &
                                                                SQRT(eigen_diff/(eigen_diff**2+omega**2))

                        END DO
                     END DO
                  ELSE
!$OMP               PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                           SHARED(ncol_local,nrow_local,col_indices,row_indices,jquad,num_integ_points,&
!$OMP                                  Eigenval_beta,fm_mat_S_beta,virtual_beta,homo_beta,omega,omega_old, &
!$OMP                                  fm_mat_S_contour_def_beta)
                     DO jjB = 1, ncol_local
                        j_global = col_indices(jjB)
                        DO iiB = 1, nrow_local
                           i_global = row_indices(iiB)

                           iocc = MAX(1, i_global-1)/virtual_beta+1
                           avirt = i_global-(iocc-1)*virtual_beta
                           eigen_diff = Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

                           IF (jquad <= num_integ_points) THEN

                              fm_mat_S_beta%local_data(iiB, jjB) = fm_mat_S_beta%local_data(iiB, jjB)* &
                                                                   SQRT((eigen_diff**2+omega_old**2)/(eigen_diff**2+omega**2))

                           ELSE

                              ! reset
                              fm_mat_S_beta%local_data(iiB, jjB) = fm_mat_S_contour_def_beta%local_data(iiB, jjB)

                              ! update
                              fm_mat_S_beta%local_data(iiB, jjB) = &
                                 fm_mat_S_beta%local_data(iiB, jjB)*eigen_diff/(eigen_diff**2-omega**2)

                           END IF

                        END DO
                     END DO
                  END IF

                  ! alpha=4.0 is valid only for closed shell systems (alpha is 2.0 in the spin-orbital basis)
                  t_start = m_walltime()
                  SELECT CASE (mm_style)
                  CASE (wfc_mm_style_gemm)
                     IF (jquad <= num_integ_points) THEN
                        CALL cp_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia_beta, alpha=alpha, &
                                     matrix_a=fm_mat_S_beta, matrix_b=fm_mat_S_beta, beta=0.0_dp, &
                                     matrix_c=fm_mat_Q_gemm_beta)
                     ELSE
                        CALL cp_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia_beta, alpha=alpha, &
                                     matrix_a=fm_mat_S_beta, matrix_b=fm_mat_S_contour_def_beta, beta=0.0_dp, &
                                     matrix_c=fm_mat_Q_gemm_beta)
                     END IF
                  CASE (wfc_mm_style_syrk)
                     ! will only compute the upper half of the matrix, which is fine, since we only use it for cholesky later
                     CALL cp_fm_syrk(uplo='U', trans='T', k=dimen_ia, alpha=alpha, matrix_a=fm_mat_S, &
                                     ia=1, ja=1, beta=0.0_dp, matrix_c=fm_mat_Q_gemm)
                  CASE DEFAULT
                     CPABORT("")
                  END SELECT

                  t_end = m_walltime()
                  actual_flop_rate = 2.0_dp*REAL(dimen_ia_beta, KIND=dp)* &
                                     dimen_RI*REAL(dimen_RI, KIND=dp)/(MAX(0.01_dp, t_end-t_start))
                  IF (para_env_RPA%mepos == 0) my_flop_rate = my_flop_rate+actual_flop_rate
                  my_num_dgemm_call = my_num_dgemm_call+1

                  ! copy/redistribute fm_mat_Q_gemm to fm_mat_Q
                  CALL cp_fm_set_all(matrix=fm_mat_Q_beta, alpha=0.0_dp)
                  CALL fm_redistribute(fm_mat_Q_gemm_beta, fm_mat_Q_beta, RPA_proc_map, buffer_rec_beta, buffer_send_beta, &
                                       number_of_send_beta, &
                                       map_send_size_beta, map_rec_size_beta, local_size_source_beta, para_env_RPA)

                  CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_Q, beta=1.0_dp, matrix_b=fm_mat_Q_beta)

               END IF ! open shell

            END IF ! im time

            IF (jquad <= num_integ_points) THEN
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
            ENDIF ! jquad

            ! get info of fm_mat_Q
            CALL cp_fm_get_info(matrix=fm_mat_Q, &
                                nrow_local=nrow_local, &
                                ncol_local=ncol_local, &
                                row_indices=row_indices, &
                                col_indices=col_indices)

            ! calcualte the trace of Q and add 1 on the diagonal
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

            IF (jquad <= num_integ_points) THEN

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

            ! save omega and reset the first_cycle flag
            first_cycle = .FALSE.
            omega_old = omega

            CALL timestop(handle3)

            ! the actual G0W0 calculation
            IF (my_do_gw) THEN
               CALL timeset(routineN//"_G0W0_matrix_operations", handle2)

               IF (jquad <= num_integ_points) THEN

                  ! calculate [1+Q(iw')]^-1
                  CALL cp_fm_cholesky_invert(fm_mat_Q)
                  ! symmetrize the result, fm_mat_R_gw is only temporary work matrix
                  CALL cp_fm_upper_to_full(fm_mat_Q, fm_mat_R_gw)

               ELSE

                  ! inverted matrix is written to fm_mat_R_gw
                  CALL cp_fm_invert(fm_mat_Q, fm_mat_R_gw)

                  CALL cp_fm_to_fm(source=fm_mat_R_gw, destination=fm_mat_Q)

               END IF

               ! periodic correction for GW
               IF (do_periodic) THEN
                  CALL calc_periodic_correction(delta_corr, delta_corr_nm, qs_env, para_env, para_env_RPA, &
                                                mp2_env%ri_g0w0%kp_grid, homo, nmo, gw_corr_lev_occ, &
                                                gw_corr_lev_virt, omega, mo_coeff, Eigenval, &
                                                matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                                first_cycle_periodic_correction, kpoints, &
                                                mp2_env%ri_g0w0%do_mo_coeff_gamma, &
                                                mp2_env%ri_g0w0%num_kp_grids, mp2_env%ri_g0w0%eps_kpoint, &
                                                mp2_env%ri_g0w0%do_extra_kpoints, &
                                                mp2_env%ri_g0w0%do_aux_bas_gw, mp2_env%ri_g0w0%frac_aux_mos)
               END IF

               ! subtract 1 from the diagonal to get rid of exchange self-energy
!$OMP           PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,i_global,j_global) &
!$OMP                       SHARED(ncol_local,nrow_local,col_indices,row_indices,fm_mat_Q,dimen_RI)
               DO jjB = 1, ncol_local
                  j_global = col_indices(jjB)
                  DO iiB = 1, nrow_local
                     i_global = row_indices(iiB)
                     IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
                        fm_mat_Q%local_data(iiB, jjB) = fm_mat_Q%local_data(iiB, jjB)-1.0_dp
                     END IF
                  END DO
               END DO

               ! S_work_(nm)Q = B_(nm)P * ([1+Q]^-1-1)_PQ
               CALL timeset(routineN//"_mult_B_f(Pi)_gw", handle3)
               CALL cp_gemm(transa="N", transb="N", m=dimen_nm_gw, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_mat_S_gw, matrix_b=fm_mat_Q, beta=0.0_dp, &
                            matrix_c=fm_mat_S_gw_work)
               IF (my_open_shell) THEN
                  CALL cp_gemm(transa="N", transb="N", m=dimen_nm_gw, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                               matrix_a=fm_mat_S_gw_beta, matrix_b=fm_mat_Q, beta=0.0_dp, &
                               matrix_c=fm_mat_S_gw_work_beta)
               END IF

               CALL timestop(handle3)

               ! vector W_(nm) = S_work_(nm)Q * [B_(nm)Q]^T

               CALL cp_fm_get_info(matrix=fm_mat_S_gw, &
                                   nrow_local=nrow_local, &
                                   ncol_local=ncol_local, &
                                   row_indices=row_indices, &
                                   col_indices=col_indices)

               vec_W_gw = 0.0_dp
               IF (my_open_shell) THEN
                  vec_W_gw_beta = 0.0_dp
               END IF

               DO iiB = 1, nrow_local
                  nm_global = row_indices(iiB)
                  DO jjB = 1, ncol_local
                     vec_W_gw(nm_global) = vec_W_gw(nm_global)+ &
                                           fm_mat_S_gw_work%local_data(iiB, jjB)*fm_mat_S_gw%local_data(iiB, jjB)
                     IF (my_open_shell) THEN
                        vec_W_gw_beta(nm_global) = vec_W_gw_beta(nm_global)+ &
                                                   fm_mat_S_gw_work_beta%local_data(iiB, jjB)*fm_mat_S_gw_beta%local_data(iiB, jjB)
                     END IF
                  END DO

                  ! transform the index nm of vec_W_gw back to n and m, formulae copied from Mauro's code
                  n_global = MAX(1, nm_global-1)/nmo+1
                  m_global = nm_global-(n_global-1)*nmo
                  n_global = n_global+homo-gw_corr_lev_occ

                  IF (my_open_shell) THEN
                     n_global_beta = MAX(1, nm_global-1)/nmo+1
                     m_global_beta = nm_global-(n_global_beta-1)*nmo
                     n_global_beta = n_global_beta+homo_beta-gw_corr_lev_occ_beta
                  END IF

                  ! compute self-energy for imaginary frequencies
                  IF (jquad <= num_integ_points) THEN

                     DO iquad = 1, num_fit_points

                        ! for occ orbitals, we compute the self-energy for negative frequencies
                        IF (n_global <= homo) THEN
                           sign_occ_virt = -1.0_dp
                        ELSE
                           sign_occ_virt = 1.0_dp
                        END IF

                        omega_i = vec_omega_fit_gw(iquad)*sign_occ_virt

                        ! set the Fermi energy for occ orbitals slightly above the HOMO and
                        ! for virt orbitals slightly below the LUMO
                        IF (n_global <= homo) THEN
                           e_fermi = Eigenval(homo)+fermi_level_offset
                        ELSE
                           e_fermi = Eigenval(homo+1)-fermi_level_offset
                        END IF

                        ! add here the periodic correction
                        IF (do_periodic .AND. col_indices(1) == 1) THEN
                           delta_corr_nn = delta_corr_nm(n_global, m_global)
                        ELSE
                           delta_corr_nn = 0.0_dp
                        END IF

                        ! update the self-energy (use that vec_W_gw(iw) is symmetric), divide the integration
                        ! weight by 2, because the integration is from -infty to +infty and not just 0 to +infty
                        ! as for RPA, also we need for virtual orbitals a complex conjugate
                        vec_Sigma_c_gw(n_global-homo+gw_corr_lev_occ, iquad) = &
                           vec_Sigma_c_gw(n_global-homo+gw_corr_lev_occ, iquad)- &
                           0.5_dp/pi*wj(jquad)/2.0_dp*(vec_W_gw(nm_global)+delta_corr_nn)* &
                           (1.0_dp/(im_unit*(omega+omega_i)+e_fermi-Eigenval(m_global))+ &
                            1.0_dp/(im_unit*(-omega+omega_i)+e_fermi-Eigenval(m_global)))

                        ! the same for beta
                        IF (my_open_shell) THEN
                           ! for occ orbitals, we compute the self-energy for negative frequencies
                           IF (n_global_beta <= homo_beta) THEN
                              sign_occ_virt = -1.0_dp
                           ELSE
                              sign_occ_virt = 1.0_dp
                           END IF

                           omega_i = vec_omega_fit_gw(iquad)*sign_occ_virt

                           ! set the Fermi energy for occ orbitals slightly above the HOMO and
                           ! for virt orbitals slightly below the LUMO
                           IF (n_global_beta <= homo_beta) THEN
                              e_fermi = Eigenval_beta(homo_beta)+fermi_level_offset
                           ELSE
                              e_fermi = Eigenval_beta(homo_beta+1)-fermi_level_offset
                           END IF

                           ! add here the periodic correction
                           IF (do_periodic .AND. col_indices(1) == 1) THEN
                              delta_corr_nn_beta = delta_corr_nm(n_global_beta, m_global)
                           ELSE
                              delta_corr_nn_beta = 0.0_dp
                           END IF

                           ! update the self-energy (use that vec_W_gw(iw) is symmetric), divide the integration
                           ! weight by 2, because the integration is from -infty to +infty and not just 0 to +infty
                           ! as for RPA, also we need for virtual orbitals a complex conjugate
                           vec_Sigma_c_gw_beta(n_global_beta-homo_beta+gw_corr_lev_occ_beta, iquad) = &
                              vec_Sigma_c_gw_beta(n_global_beta-homo_beta+gw_corr_lev_occ_beta, iquad)- &
                              !                              0.5_dp/pi*wj(jquad)/2.0_dp*vec_W_gw_beta(nm_global)/ &
                              !                              (im_unit*(omega+omega_i)+e_fermi-Eigenval_beta(m_global))- &
                              !                              0.5_dp/pi*wj(jquad)/2.0_dp*vec_W_gw_beta(nm_global)/ &
                              !                              (im_unit*(-omega+omega_i)+e_fermi-Eigenval_beta(m_global))
                              0.5_dp/pi*wj(jquad)/2.0_dp*(vec_W_gw_beta(nm_global)+delta_corr_nn_beta)* &
                              (1.0_dp/(im_unit*(omega+omega_i)+e_fermi-Eigenval_beta(m_global))+ &
                               1.0_dp/(im_unit*(-omega+omega_i)+e_fermi-Eigenval_beta(m_global)))

                        END IF
                     END DO ! iquad
                  END IF ! check imaginary frequency

                  IF (do_contour_def) THEN

                     IF (my_open_shell) THEN

                        CALL contour_def_integrate_and_residues(vec_Sigma_c_gw_real_freq, wj, vec_W_gw, Eigenval, &
                                                                omega, contour_def_offset, e_fermi, &
                                                                mo_from_jquad, residue_from_jquad, alpha_beta_from_jquad, &
                                                                contour_def_start, contour_def_end, n_global, jquad, &
                                                                num_integ_points, m_global, nm_global, &
                                                                alpha_beta_case=.TRUE., alpha=.TRUE.)

                        CALL contour_def_integrate_and_residues(vec_Sigma_c_gw_real_freq_beta, wj, vec_W_gw_beta, Eigenval_beta, &
                                                                omega, contour_def_offset, e_fermi, &
                                                                mo_from_jquad, residue_from_jquad, alpha_beta_from_jquad, &
                                                                contour_def_start, contour_def_end, n_global_beta, jquad, &
                                                                num_integ_points, m_global_beta, nm_global, &
                                                                alpha_beta_case=.TRUE., beta=.TRUE.)

                     ELSE

                        CALL contour_def_integrate_and_residues(vec_Sigma_c_gw_real_freq, wj, vec_W_gw, Eigenval, &
                                                                omega, contour_def_offset, e_fermi, &
                                                                mo_from_jquad, residue_from_jquad, alpha_beta_from_jquad, &
                                                                contour_def_start, contour_def_end, n_global, jquad, &
                                                                num_integ_points, m_global, nm_global)
                     END IF

                  END IF ! contour deformation

               END DO ! iiB

               CALL timestop(handle2)

            END IF ! GW

         END DO ! jquad

         CALL mp_sum(Erpa, para_env%group)
         Erpa = Erpa/(pi*2.0_dp)
         IF (do_minimax_quad) Erpa = Erpa/2.0_dp

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
         IF (my_do_gw) THEN

            CALL timeset(routineN//"_G0W0_fit", handle3)

            CALL mp_sum(vec_Sigma_c_gw, para_env%group)

            IF (do_periodic .AND. mp2_env%ri_g0w0%do_average_deg_levels) THEN

               CALL average_degenerate_levels(vec_Sigma_c_gw, Eigenval(1+homo-gw_corr_lev_occ:homo+gw_corr_lev_virt), &
                                              mp2_env%ri_g0w0%eps_eigenval)
               IF (my_open_shell) THEN
                  CALL average_degenerate_levels(vec_Sigma_c_gw_beta, &
                                                 Eigenval_beta(1+homo_beta-gw_corr_lev_occ_beta: &
                                                               homo_beta+gw_corr_lev_virt_beta), &
                                                 mp2_env%ri_g0w0%eps_eigenval)
               END IF
            END IF

            IF (my_open_shell) THEN
               CALL mp_sum(vec_Sigma_c_gw_beta, para_env%group)
            END IF

            IF (do_contour_def) THEN

               CALL mp_sum(vec_Sigma_c_gw_real_freq, para_env%group)

               IF (my_open_shell) THEN
                  CALL mp_sum(vec_Sigma_c_gw_real_freq_beta, para_env%group)
               END IF

            END IF

            CALL mp_sync(para_env%group)

            ! fit the self-energy on imaginary frequency axis and evaluate the fit on the MO energy of the SCF
            DO n_level_gw = 1, gw_corr_lev_tot
               ! processes perform different fits
               IF (MODULO(n_level_gw, para_env%num_pe) /= para_env%mepos) CYCLE

               CALL fit_and_continuation(vec_gw_energ, vec_gw_energ_error_fit, vec_omega_fit_gw, &
                                         z_value, m_value, vec_Sigma_c_gw, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1), &
                                         Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, num_poles, &
                                         num_fit_points, max_iter_fit, crossing_search, homo, check_fit, stop_crit, &
                                         fermi_level_offset)

               IF (my_open_shell) THEN
                  CALL fit_and_continuation( &
                     vec_gw_energ_beta, vec_gw_energ_error_fit_beta, vec_omega_fit_gw, &
                     z_value_beta, m_value_beta, vec_Sigma_c_gw_beta, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2), &
                     Eigenval_beta, Eigenval_scf_beta, n_level_gw, &
                     gw_corr_lev_occ_beta, num_poles, &
                     num_fit_points, max_iter_fit, crossing_search, homo_beta, check_fit, stop_crit, &
                     fermi_level_offset)

               END IF

            END DO ! n_level_gw

            CALL mp_sum(vec_gw_energ_error_fit, para_env%group)
            CALL mp_sum(vec_gw_energ, para_env%group)
            CALL mp_sum(z_value, para_env%group)
            CALL mp_sum(m_value, para_env%group)

            IF (my_open_shell) THEN
               CALL mp_sum(vec_gw_energ_error_fit_beta, para_env%group)
               CALL mp_sum(vec_gw_energ_beta, para_env%group)
               CALL mp_sum(z_value_beta, para_env%group)
               CALL mp_sum(m_value_beta, para_env%group)
            END IF

            remove_neg_virt_energies = mp2_env%ri_g0w0%remove_neg_virt_energies

            ! in case we do the contour deformation technique, we transfer vec_Sigma_c_gw_real_freq to
            ! vec_gw_energ and we compute the z-value, m-value
            ! the results from the analytic continuation are therefore overwritten
            IF (do_contour_def) THEN

               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) ''
                  WRITE (unit_nr, '(T3,A,I4,A,I3)') 'The following MOs have been corrected by contour deformation: MO', &
                     contour_def_start, '  -  MO', contour_def_end
                  WRITE (unit_nr, *) ''
               END IF

               index_contour_def = 1

               ! only correct levels for which contour deformation is enabled
               DO n_level_gw = contour_def_start-(homo-gw_corr_lev_occ), contour_def_end-(homo-gw_corr_lev_occ)

                  ! reset the values from analytic continuation
                  vec_gw_energ(n_level_gw) = 0.0_dp
                  vec_gw_energ_error_fit(n_level_gw) = 0.0_dp
                  z_value(n_level_gw) = 0.0_dp
                  m_value(n_level_gw) = 0.0_dp

                  CALL compute_z_and_m_contour_def(vec_gw_energ, z_value, m_value, &
                                                   vec_Sigma_c_gw_real_freq, contour_def_offset, &
                                                   n_level_gw, index_contour_def)

                  index_contour_def = index_contour_def+1

               END DO

               ! the same for beta
               IF (my_open_shell) THEN

                  index_contour_def = 1

                  ! only correct levels for which contour deformation is enabled
                  DO n_level_gw = contour_def_start-(homo_beta-gw_corr_lev_occ_beta), &
                                  contour_def_end-(homo_beta-gw_corr_lev_occ_beta)

                     ! reset the values from analytic continuation
                     vec_gw_energ_beta(n_level_gw) = 0.0_dp
                     vec_gw_energ_error_fit_beta(n_level_gw) = 0.0_dp
                     z_value_beta(n_level_gw) = 0.0_dp
                     m_value_beta(n_level_gw) = 0.0_dp

                     CALL compute_z_and_m_contour_def(vec_gw_energ_beta, z_value_beta, m_value_beta, &
                                                      vec_Sigma_c_gw_real_freq_beta, contour_def_offset, &
                                                      n_level_gw, index_contour_def)

                     index_contour_def = index_contour_def+1

                  END DO

               END IF

               DEALLOCATE (real_freq)
               DEALLOCATE (mo_from_jquad)
               DEALLOCATE (residue_from_jquad)
               DEALLOCATE (alpha_beta_from_jquad)

            END IF

            ! print the quasiparticle energies and update Eigenval in case you do eigenvalue self-consistent GW
            IF (my_open_shell) THEN

               CALL print_and_update_for_ev_sc( &
                  vec_gw_energ, vec_gw_energ_error_fit, &
                  z_value, m_value, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1), Eigenval, &
                  Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                  count_ev_sc_GW, crossing_search, homo, nmo, unit_nr, mp2_env%ri_g0w0%print_gw_details, &
                  do_alpha=.TRUE.)

               CALL print_and_update_for_ev_sc( &
                  vec_gw_energ_beta, vec_gw_energ_error_fit_beta, &
                  z_value_beta, m_value_beta, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 2), Eigenval_beta, &
                  Eigenval_last_beta, Eigenval_scf_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta, gw_corr_lev_tot, &
                  count_ev_sc_GW, crossing_search, homo_beta, nmo, unit_nr, mp2_env%ri_g0w0%print_gw_details, &
                  do_beta=.TRUE.)

            ELSE

               CALL print_and_update_for_ev_sc( &
                  vec_gw_energ, vec_gw_energ_error_fit, &
                  z_value, m_value, mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, 1), Eigenval, &
                  Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                  count_ev_sc_GW, crossing_search, homo, nmo, unit_nr, mp2_env%ri_g0w0%print_gw_details)

            END IF

            CALL timestop(handle3)

            ! if HOMO-LUMO gap differs by less then mp2_env%ri_g0w0%eps_ev_sc_iter, exit ev sc GW loop
            IF (ABS(Eigenval(homo)-Eigenval_last(homo)-Eigenval(homo+1)+Eigenval_last(homo+1)) &
                < mp2_env%ri_g0w0%eps_ev_sc_iter) THEN
               EXIT
            END IF

         END IF ! my_do_gw if

      END DO !ev_sc_gw_loop

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

      IF (my_do_gw) THEN
         CALL cp_fm_release(fm_mat_S_gw_work)
         DEALLOCATE (vec_Sigma_c_gw)
         DEALLOCATE (vec_omega_fit_gw)
         DEALLOCATE (vec_omega_fit_gw_sign)
         DEALLOCATE (Lambda)
         DEALLOCATE (Lambda_without_offset)
         DEALLOCATE (Lambda_Re)
         DEALLOCATE (Lambda_Im)
         DEALLOCATE (z_value)
         DEALLOCATE (m_value)
         DEALLOCATE (vec_gw_energ)
         DEALLOCATE (vec_gw_energ_error_fit)
         DEALLOCATE (mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw)
         DEALLOCATE (Eigenval_last)
         DEALLOCATE (Eigenval_scf)
         IF (do_ri_Sigma_x) THEN
            DEALLOCATE (vec_Sigma_x_gw)
         END IF
         IF (my_open_shell) THEN
            CALL cp_fm_release(fm_mat_S_gw_work_beta)
            DEALLOCATE (vec_W_gw_beta)
            DEALLOCATE (vec_Sigma_c_gw_beta)
            DEALLOCATE (z_value_beta)
            DEALLOCATE (m_value_beta)
            DEALLOCATE (vec_gw_energ_beta)
            DEALLOCATE (vec_gw_energ_error_fit_beta)
            DEALLOCATE (Eigenval_last_beta)
            DEALLOCATE (Eigenval_scf_beta)
            IF (do_ri_Sigma_x) THEN
               DEALLOCATE (vec_Sigma_x_gw_beta)
            END IF
         END IF
         IF (do_contour_def) THEN
            DEALLOCATE (vec_Sigma_c_gw_real_freq)
            CALL cp_fm_release(fm_mat_S_contour_def)
            IF (my_open_shell) THEN
               DEALLOCATE (vec_Sigma_c_gw_real_freq_beta)
               CALL cp_fm_release(fm_mat_S_contour_def_beta)
            END IF
         END IF
         IF (do_periodic) THEN
            CALL dbcsr_deallocate_matrix_set(matrix_berry_re_mo_mo)
            CALL dbcsr_deallocate_matrix_set(matrix_berry_im_mo_mo)
            CALL kpoint_release(kpoints)
         END IF
      END IF

      IF (do_im_time) THEN

         DO iquad = 1, num_integ_points
            CALL cp_fm_release(fm_scaled_dm_occ_tau(iquad)%matrix)
         END DO
         DEALLOCATE (fm_scaled_dm_occ_tau)

         DO iquad = 1, num_integ_points
            CALL cp_fm_release(fm_scaled_dm_virt_tau(iquad)%matrix)
         END DO
         DEALLOCATE (fm_scaled_dm_virt_tau)

         CALL cp_fm_release(fm_mo_coeff_occ)
         CALL cp_fm_release(fm_mo_coeff_virt)

         IF (my_open_shell) THEN
            CALL cp_fm_release(fm_mo_coeff_occ_beta)
            CALL cp_fm_release(fm_mo_coeff_virt_beta)
         END IF

         CALL cp_fm_release(fm_mat_L)
         CALL cp_fm_release(fm_mat_work)

         CALL dbcsr_release_p(mat_dm_loc_occ)
         CALL dbcsr_release_p(mat_dm_loc_virt)

         CALL dbcsr_deallocate_matrix_set(mat_dm_loc_occ_cut)
         CALL dbcsr_deallocate_matrix_set(mat_dm_loc_virt_cut)

         DO iquad = 1, num_integ_points
            CALL cp_fm_release(fm_mat_P_omega(iquad)%matrix)
         END DO
         DEALLOCATE (fm_mat_P_omega)

         CALL dbcsr_release(mat_M_P_munu_occ%matrix)
         DEALLOCATE (mat_M_P_munu_occ%matrix)
         CALL dbcsr_release(mat_M_P_munu_virt%matrix)
         DEALLOCATE (mat_M_P_munu_virt%matrix)
         CALL dbcsr_release(mat_P_global_copy%matrix)
         DEALLOCATE (mat_P_global_copy%matrix)
         CALL dbcsr_release(mat_L%matrix)
         DEALLOCATE (mat_L%matrix)

!         CALL dbcsr_release(mat_M_mu_Pnu_occ%matrix)
!         DEALLOCATE(mat_M_mu_Pnu_occ%matrix)
!         CALL dbcsr_release(mat_M_mu_Pnu_virt%matrix)
!         DEALLOCATE(mat_M_mu_Pnu_virt%matrix)

         CALL dbcsr_deallocate_matrix_set(mat_M_mu_Pnu_occ)
         CALL dbcsr_deallocate_matrix_set(mat_M_mu_Pnu_virt)

         CALL dbcsr_deallocate_matrix_set(mat_P_omega)

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

      END IF

      DEALLOCATE (tj)
      DEALLOCATE (wj)
      DEALLOCATE (Q_log)
      DEALLOCATE (trace_Qomega)

      CALL timestop(handle)

   END SUBROUTINE rpa_num_int

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
      REAL(KIND=dp)                                      :: a_scaling_ext
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      INTEGER                                            :: homo, virtual
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tj_ext, wj_ext
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      INTEGER, OPTIONAL                                  :: homo_beta, virtual_beta, dimen_ia_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_S_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_scaling_factor', &
         routineP = moduleN//':'//routineN

      INTEGER :: avirt, color_col, color_col_beta, color_row, color_row_beta, comm_col, &
         comm_col_beta, comm_row, comm_row_beta, handle, i_global, icycle, iiB, iocc, jjB, jquad, &
         ncol_local, ncol_local_beta, nrow_local, nrow_local_beta
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, col_indices_beta, &
                                                            row_indices, row_indices_beta
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp) :: a_high, a_low, a_scaling, conv_param, eigen_diff, eps, first_deriv, four, &
         left_term, one, pig, right_term, right_term_ref, right_term_ref_beta, step, two, zero
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cottj, D_ia, D_ia_beta, iaia_RI, &
                                                            iaia_RI_beta, iaia_RI_dp, &
                                                            iaia_RI_dp_beta, M_ia, M_ia_beta, tj, &
                                                            wj
      TYPE(cp_para_env_type), POINTER                    :: para_env_col, para_env_col_beta, &
                                                            para_env_row, para_env_row_beta

      CALL timeset(routineN, handle)

      my_open_shell = .FALSE.
      IF (PRESENT(homo_beta) .AND. &
          PRESENT(virtual_beta) .AND. &
          PRESENT(dimen_ia_beta) .AND. &
          PRESENT(Eigenval_beta) .AND. &
          PRESENT(fm_mat_S_beta)) my_open_shell = .TRUE.

      ZERO = 0.0_dp
      ONE = 1.0_dp
      TWO = 2.0_dp
      FOUR = 4.0_dp
      PIG = pi
      eps = 1.0E-10_dp

      ALLOCATE (cottj(num_integ_points))

      ALLOCATE (tj(num_integ_points))

      ALLOCATE (wj(num_integ_points))

      ! calculate the cotangent of the abscissa tj
      DO jquad = 1, num_integ_points
         tj(jquad) = tj_ext(jquad)
         wj(jquad) = wj_ext(jquad)
         cottj(jquad) = ONE/TAN(tj(jquad))
      END DO

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
         M_ia(iiB) = D_ia(iiB)*D_ia(iiB)+TWO*D_ia(iiB)*iaia_RI(iiB)
      END DO

      right_term_ref = ZERO
      DO iiB = 1, nrow_local
         right_term_ref = right_term_ref+(SQRT(M_ia(iiB))-D_ia(iiB)-iaia_RI(iiB))
      END DO
      right_term_ref = right_term_ref/TWO
      ! right_term_ref=accurate_sum((SQRT(M_ia)-D_ia-iaia_RI))/2.0_dp

      ! sum the result with the processes of the RPA_group having the same col
      color_col = fm_mat_S%matrix_struct%context%mepos(2)
      CALL mp_comm_split_direct(para_env_RPA%group, comm_col, color_col)
      NULLIFY (para_env_col)
      CALL cp_para_env_create(para_env_col, comm_col)

      ! allocate communication array for columns
      CALL mp_sum(right_term_ref, para_env_col%group)

      ! In the open shell case do point 1-2-3 for the beta spin
      IF (my_open_shell) THEN
         !XXX CALL cp_fm_to_fm(source=fm_mat_S_beta,destination=fm_mat_G_beta)
         CALL cp_fm_get_info(matrix=fm_mat_S_beta, &
                             nrow_local=nrow_local_beta, &
                             ncol_local=ncol_local_beta, &
                             row_indices=row_indices_beta, &
                             col_indices=col_indices_beta)

         ALLOCATE (iaia_RI_dp_beta(nrow_local_beta))
         iaia_RI_dp_beta = 0.0_dp

         DO jjB = 1, ncol_local_beta
            DO iiB = 1, nrow_local_beta
               iaia_RI_dp_beta(iiB) = iaia_RI_dp_beta(iiB)+fm_mat_S_beta%local_data(iiB, jjB)*fm_mat_S_beta%local_data(iiB, jjB)
            END DO
         END DO

         color_row_beta = fm_mat_S_beta%matrix_struct%context%mepos(1)
         CALL mp_comm_split_direct(para_env_RPA%group, comm_row_beta, color_row_beta)
         NULLIFY (para_env_row_beta)
         CALL cp_para_env_create(para_env_row_beta, comm_row_beta)

         CALL mp_sum(iaia_RI_dp_beta, para_env_row_beta%group)

         ALLOCATE (iaia_RI_beta(nrow_local_beta))
         DO iiB = 1, nrow_local_beta
            iaia_RI_beta(iiB) = iaia_RI_dp_beta(iiB)
         END DO
         DEALLOCATE (iaia_RI_dp_beta)

         ALLOCATE (D_ia_beta(nrow_local_beta))

         ALLOCATE (M_ia_beta(nrow_local_beta))

         DO iiB = 1, nrow_local_beta
            i_global = row_indices_beta(iiB)

            iocc = MAX(1, i_global-1)/virtual_beta+1
            avirt = i_global-(iocc-1)*virtual_beta
            eigen_diff = Eigenval_beta(avirt+homo_beta)-Eigenval_beta(iocc)

            D_ia_beta(iiB) = eigen_diff
         END DO

         DO iiB = 1, nrow_local_beta
            M_ia_beta(iiB) = D_ia_beta(iiB)*D_ia_beta(iiB)+TWO*D_ia_beta(iiB)*iaia_RI_beta(iiB)
         END DO

         right_term_ref_beta = ZERO
         DO iiB = 1, nrow_local_beta
            right_term_ref_beta = right_term_ref_beta+(SQRT(M_ia_beta(iiB))-D_ia_beta(iiB)-iaia_RI_beta(iiB))
         END DO
         right_term_ref_beta = right_term_ref_beta/TWO

         ! sum the result with the processes of the RPA_group having the same col
         color_col_beta = fm_mat_S_beta%matrix_struct%context%mepos(2)
         CALL mp_comm_split_direct(para_env_RPA%group, comm_col_beta, color_col_beta)
         NULLIFY (para_env_col_beta)
         CALL cp_para_env_create(para_env_col_beta, comm_col_beta)

         CALL mp_sum(right_term_ref_beta, para_env_col_beta%group)

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
                                ZERO, ONE, M_ia, cottj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_row, para_env_row_beta)
         left_term = left_term/FOUR/PIG*a_scaling

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
                                      ZERO, ONE, M_ia, cottj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                      nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                      para_env, para_env_row, para_env_row_beta)
               left_term = left_term/FOUR/PIG*a_scaling

               IF (ABS(left_term) > ABS(right_term)) THEN
                  a_high = a_scaling
               ELSE
                  a_low = a_scaling
               END IF

               IF (ABS(a_high-a_low) < 1.0D-5) EXIT

            END DO

         END IF
      END IF

      a_scaling_ext = a_scaling
      CALL mp_bcast(a_scaling_ext, 0, para_env%group)

      DEALLOCATE (cottj)
      DEALLOCATE (tj)
      DEALLOCATE (wj)
      DEALLOCATE (iaia_RI)
      DEALLOCATE (D_ia)
      DEALLOCATE (M_ia)
      CALL cp_para_env_release(para_env_row)
      CALL cp_para_env_release(para_env_col)

      IF (my_open_shell) THEN
         DEALLOCATE (iaia_RI_beta)
         DEALLOCATE (D_ia_beta)
         DEALLOCATE (M_ia_beta)
         CALL cp_para_env_release(para_env_row_beta)
         CALL cp_para_env_release(para_env_col_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param a_scaling ...
!> \param left_term ...
!> \param first_deriv ...
!> \param num_integ_points ...
!> \param my_open_shell ...
!> \param ZERO ...
!> \param ONE ...
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
                                ZERO, ONE, M_ia, cottj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                nrow_local, nrow_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_row, para_env_row_beta)
      REAL(KIND=dp)                                      :: a_scaling, left_term, first_deriv
      INTEGER                                            :: num_integ_points
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp)                                      :: zero, one
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: M_ia, cottj, wj, D_ia, D_ia_beta, &
                                                            M_ia_beta
      INTEGER                                            :: nrow_local, nrow_local_beta, &
                                                            num_integ_group, color_rpa_group
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_row, para_env_row_beta

      INTEGER                                            :: iiB, jquad
      REAL(KIND=dp)                                      :: first_deriv_beta, left_term_beta, omega

      left_term = ZERO
      first_deriv = ZERO
      left_term_beta = ZERO
      first_deriv_beta = ZERO
      DO jquad = 1, num_integ_points
         ! parallelize over integration points
         IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE
         omega = a_scaling*cottj(jquad)

         DO iiB = 1, nrow_local
            ! parallelize over ia elements in the para_env_row group
            IF (MODULO(iiB, para_env_row%num_pe) /= para_env_row%mepos) CYCLE
            ! calculate left_term
            left_term = left_term+wj(jquad)* &
                        (LOG(ONE+(M_ia(iiB)-D_ia(iiB)**2)/(omega**2+D_ia(iiB)**2))- &
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
                                (LOG(ONE+(M_ia_beta(iiB)-D_ia_beta(iiB)**2)/(omega**2+D_ia_beta(iiB)**2))- &
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
!> \param tau_wj_2d ...
!> \param omega_tj ...
!> \param E_min ...
!> \param E_max ...
!> \param max_error ...
!> \param num_points_per_magnitude ...
! **************************************************************************************************
   SUBROUTINE get_rpa_least_sq_weights_tau(num_integ_points, tau_tj, tau_wj_2d, omega_tj, &
                                           E_min, E_max, max_error, num_points_per_magnitude)

      INTEGER                                            :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tau_wj_2d
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: omega_tj
      REAL(KIND=dp)                                      :: E_min, E_max, max_error
      INTEGER                                            :: num_points_per_magnitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_rpa_least_sq_weights_tau', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      REAL(KIND=dp)                                      :: chi2_min_jquad, multiplicator, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, tau_wj_work, vec_UTy, work, &
                                                            work_array, x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_U, mat_VSigma, mat_VT

      CALL timeset(routineN, handle)

      ! take ten points per 10-interval
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
      ALLOCATE (mat_VT(num_x_nodes, num_integ_points))

      mat_VT = 0.0_dp
      ! double the value nessary for 'A' to achieve good performance
      lwork = 8*num_integ_points*num_integ_points+12*num_integ_points+2*num_x_nodes
      ALLOCATE (work(lwork))
      work = 0.0_dp
      ALLOCATE (iwork(8*num_integ_points))
      iwork = 0
      ALLOCATE (mat_VSigma(num_integ_points, num_x_nodes))
      mat_VSigma = 0.0_dp
      ALLOCATE (vec_UTy(num_x_nodes))
      vec_UTy = 0.0_dp

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
                     mat_VT, num_x_nodes, work, lwork, iwork, info)

         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         ! 1) V*Sigma
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_VSigma(iii, jjj) = mat_VT(jjj, iii)/sing_values(jjj)
            END DO
         END DO

         ! 2) U^T y
         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)

         ! 3) (V*Sigma) * (U^T y)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_VSigma, num_integ_points, vec_UTy, &
                    num_x_nodes, 0.0_dp, tau_wj_work, num_integ_points)

         tau_wj_2d(jquad, :) = tau_wj_work(:)

      END DO ! jquad

      CALL calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                   y_values, num_integ_points, num_x_nodes)

      DEALLOCATE (x_values, y_values, mat_A, tau_wj_work, work_array, sing_values, mat_U, mat_VT, &
                  work, iwork, mat_VSigma, vec_UTy)

      CALL timestop(handle)

   END SUBROUTINE get_rpa_least_sq_weights_tau

! **************************************************************************************************
!> \brief Calculate the matrix mat_N_gw containing the second derivatives
!>        with respect to the fitting parameters. The second derivatives are
!>        calculated numerically by finite differences.
!> \param N_ij matrix element
!> \param Lambda fitting parameters
!> \param Sigma_c ...
!> \param vec_omega_fit_gw ...
!> \param i ...
!> \param j ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param n_level_gw ...
!> \param h  ...
! **************************************************************************************************
   SUBROUTINE calc_mat_N(N_ij, Lambda, Sigma_c, vec_omega_fit_gw, i, j, &
                         num_poles, num_fit_points, n_level_gw, h)
      REAL(KIND=dp)                                      :: N_ij
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sigma_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_fit_gw
      INTEGER                                            :: i, j, num_poles, num_fit_points, &
                                                            n_level_gw
      REAL(KIND=dp)                                      :: h

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mat_N', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit, re_unit
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda_tmp
      INTEGER                                            :: handle, num_var
      REAL(KIND=dp)                                      :: chi2, chi2_sum

      CALL timeset(routineN, handle)

      num_var = 2*num_poles+1
      ALLOCATE (Lambda_tmp(num_var))
      Lambda_tmp = (0.0_dp, 0.0_dp)
      chi2_sum = 0.0_dp
      re_unit = (1.0_dp, 0.0_dp)
      im_unit = (0.0_dp, 1.0_dp)

      !test
      Lambda_tmp(:) = Lambda(:)
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)

      ! Fitting parameters with offset h
      Lambda_tmp(:) = Lambda(:)
      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2)+h*re_unit
      ELSE
         Lambda_tmp((i+1)/2) = Lambda_tmp((i+1)/2)+h*im_unit
      END IF
      IF (MODULO(j, 2) == 0) THEN
         Lambda_tmp(j/2) = Lambda_tmp(j/2)+h*re_unit
      ELSE
         Lambda_tmp((j+1)/2) = Lambda_tmp((j+1)/2)+h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum+chi2

      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2)-2.0_dp*h*re_unit
      ELSE
         Lambda_tmp((i+1)/2) = Lambda_tmp((i+1)/2)-2.0_dp*h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum-chi2

      IF (MODULO(j, 2) == 0) THEN
         Lambda_tmp(j/2) = Lambda_tmp(j/2)-2.0_dp*h*re_unit
      ELSE
         Lambda_tmp((j+1)/2) = Lambda_tmp((j+1)/2)-2.0_dp*h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum+chi2

      IF (MODULO(i, 2) == 0) THEN
         Lambda_tmp(i/2) = Lambda_tmp(i/2)+2.0_dp*h*re_unit
      ELSE
         Lambda_tmp((i+1)/2) = Lambda_tmp((i+1)/2)+2.0_dp*h*im_unit
      END IF
      CALL calc_chi2(chi2, Lambda_tmp, Sigma_c, vec_omega_fit_gw, num_poles, &
                     num_fit_points, n_level_gw)
      chi2_sum = chi2_sum-chi2

      ! Second derivative with symmetric difference quotient
      N_ij = 1.0_dp/2.0_dp*chi2_sum/(4.0_dp*h*h)

      DEALLOCATE (Lambda_tmp)

      CALL timestop(handle)

   END SUBROUTINE calc_mat_N

! **************************************************************************************************
!> \brief Calculate chi2
!> \param chi2 ...
!> \param Lambda fitting parameters
!> \param Sigma_c ...
!> \param vec_omega_fit_gw ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param n_level_gw ...
! **************************************************************************************************
   SUBROUTINE calc_chi2(chi2, Lambda, Sigma_c, vec_omega_fit_gw, num_poles, &
                        num_fit_points, n_level_gw)
      REAL(KIND=dp)                                      :: chi2
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sigma_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_fit_gw
      INTEGER                                            :: num_poles, num_fit_points, n_level_gw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_chi2', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: func_val, im_unit
      INTEGER                                            :: handle, iii, jjj, kkk

      CALL timeset(routineN, handle)

      im_unit = (0.0_dp, 1.0_dp)
      chi2 = 0.0_dp
      DO kkk = 1, num_fit_points
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            ! calculate value of the fit function
            func_val = func_val+Lambda(jjj)/(im_unit*vec_omega_fit_gw(kkk)-Lambda(jjj+1))
         END DO
         chi2 = chi2+(ABS(Sigma_c(n_level_gw, kkk)-func_val))**2
      END DO

      CALL timestop(handle)

   END SUBROUTINE calc_chi2

! **************************************************************************************************
!> \brief Calculate chi2
!> \param chi2 ...
!> \param beta fitting parameters
!> \param Sigma_c ...
!> \param vec_tau_fit_gw ...
!> \param num_poles ...
!> \param num_fit_points ...
! **************************************************************************************************
   SUBROUTINE calc_chi2_im_time(chi2, beta, Sigma_c, vec_tau_fit_gw, num_poles, &
                                num_fit_points)
      REAL(KIND=dp)                                      :: chi2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: beta, Sigma_c, vec_tau_fit_gw
      INTEGER                                            :: num_poles, num_fit_points

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_chi2_im_time', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: func_val
      INTEGER                                            :: handle, iii, jjj, kkk

      CALL timeset(routineN, handle)

      chi2 = 0.0_dp
      DO kkk = 1, num_fit_points
         func_val = 0.0_dp
         DO iii = 1, num_poles
            jjj = iii*3
            ! calculate value of the fit function
            IF (vec_tau_fit_gw(kkk) > 0) THEN
               func_val = func_val+beta(jjj-2)*EXP(-vec_tau_fit_gw(kkk)*beta(jjj))
            ELSE
               func_val = func_val+beta(jjj-1)*EXP(vec_tau_fit_gw(kkk)*beta(jjj))
            END IF
         END DO
         chi2 = chi2+(ABS(Sigma_c(kkk)-func_val))**2

      END DO

      CALL timestop(handle)

   END SUBROUTINE calc_chi2_im_time

! **************************************************************************************************
!> \brief Calculate chi2 when setting up the grid for tau integration with cosine
!> \param chi2 ...
!> \param omega ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param x_values ...
!> \param y_values ...
!> \param num_integ_points ...
!> \param num_x_nodes ...
! **************************************************************************************************
   SUBROUTINE calc_chi2_fit_tau_grid_with_cosine(chi2, omega, tau_tj, tau_wj_work, x_values, &
                                                 y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp)                                      :: chi2, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER                                            :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_chi2_fit_tau_grid_with_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, kkk
      REAL(KIND=dp)                                      :: func_val

      CALL timeset(routineN, handle)

      chi2 = 0.0_dp
      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_cosine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         chi2 = chi2+(ABS(y_values(kkk)-func_val))**2

      END DO

      CALL timestop(handle)

   END SUBROUTINE

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
   SUBROUTINE calc_max_error_fit_tau_grid_with_cosine(max_error, omega, tau_tj, tau_wj_work, x_values, &
                                                      y_values, num_integ_points, num_x_nodes)

      REAL(KIND=dp)                                      :: max_error, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj_work, x_values, y_values
      INTEGER                                            :: num_integ_points, num_x_nodes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_error_fit_tau_grid_with_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, kkk
      REAL(KIND=dp)                                      :: func_val

      CALL timeset(routineN, handle)

      max_error = 0.0_dp

      DO kkk = 1, num_x_nodes

         func_val = 0.0_dp

         CALL eval_fit_func_tau_grid_cosine(func_val, x_values(kkk), num_integ_points, tau_tj, tau_wj_work, omega)

         IF (ABS(y_values(kkk)-func_val) > max_error) THEN
            max_error = ABS(y_values(kkk)-func_val)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate the maximum deviation between the fit and the computed self-ener
!> \param max_dev ...
!> \param Lambda fitting parameters
!> \param Sigma_c ...
!> \param vec_omega_fit_gw ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param n_level_gw ...
! **************************************************************************************************
   SUBROUTINE calc_max_dev(max_dev, Lambda, Sigma_c, vec_omega_fit_gw, num_poles, &
                           num_fit_points, n_level_gw)
      REAL(KIND=dp)                                      :: max_dev
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: Sigma_c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_omega_fit_gw
      INTEGER                                            :: num_poles, num_fit_points, n_level_gw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_max_dev', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: func_val, im_unit
      INTEGER                                            :: handle, iii, jjj, kkk

      CALL timeset(routineN, handle)

      im_unit = (0.0_dp, 1.0_dp)
      max_dev = 0.0_dp
      DO kkk = 1, num_fit_points
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            ! calculate value of the fit function
            func_val = func_val+Lambda(jjj)/(im_unit*vec_omega_fit_gw(kkk)-Lambda(jjj+1))
         END DO
         IF (ABS(Sigma_c(n_level_gw, kkk)-func_val) > max_dev) THEN
            max_dev = ABS(Sigma_c(n_level_gw, kkk)-func_val)
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE calc_max_dev

! **************************************************************************************************
!> \brief Evaluate fit function
!> \param func_val result of fit evaluation
!> \param Lambda fitting parameters
!> \param omega real or complex energy
!> \param e_fermi the Fermi level
!> \param num_poles ...
! **************************************************************************************************
   SUBROUTINE evaluate_fit_function(func_val, Lambda, omega, e_fermi, num_poles)
      COMPLEX(KIND=dp)                                   :: func_val
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: Lambda
      COMPLEX(KIND=dp)                                   :: omega
      REAL(KIND=dp)                                      :: e_fermi
      INTEGER                                            :: num_poles

      CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_fit_function', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii, jjj

      func_val = Lambda(1)
      DO iii = 1, num_poles
         jjj = iii*2
         ! calculate value of the fit function
         func_val = func_val+Lambda(jjj)/(omega-e_fermi-Lambda(jjj+1))
      END DO

   END SUBROUTINE evaluate_fit_function

! **************************************************************************************************
!> \brief Fits the complex self-energy of n_level_gw to a multi-pole model and evaluates the
!>        fit at the energy eigenvalue of the SCF. Real part of the correlation self-energy
!>        is written to vec_gw_energ. Also calculates the statistical error of the correlation
!>        self-energy due to the fit
!> \param vec_gw_energ ...
!> \param vec_gw_energ_error_fit ...
!> \param vec_omega_fit_gw ...
!> \param z_value ...
!> \param m_value ...
!> \param vec_Sigma_c_gw ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval ...
!> \param Eigenval_scf ...
!> \param n_level_gw ...
!> \param gw_corr_lev_occ ...
!> \param num_poles ...
!> \param num_fit_points ...
!> \param max_iter_fit ...
!> \param crossing_search ...
!> \param homo ...
!> \param check_fit ...
!> \param stop_crit ...
!> \param fermi_level_offset ...
! **************************************************************************************************
   SUBROUTINE fit_and_continuation(vec_gw_energ, vec_gw_energ_error_fit, vec_omega_fit_gw, &
                                   z_value, m_value, vec_Sigma_c_gw, vec_Sigma_x_minus_vxc_gw, &
                                   Eigenval, Eigenval_scf, n_level_gw, gw_corr_lev_occ, num_poles, &
                                   num_fit_points, max_iter_fit, crossing_search, homo, check_fit, stop_crit, &
                                   fermi_level_offset)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_gw_energ, vec_gw_energ_error_fit, &
                                                            vec_omega_fit_gw, z_value, m_value
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:)                        :: vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                                            Eigenval_scf
      INTEGER                                            :: n_level_gw, gw_corr_lev_occ, num_poles, &
                                                            num_fit_points, max_iter_fit, &
                                                            crossing_search, homo
      LOGICAL                                            :: check_fit
      REAL(KIND=dp)                                      :: stop_crit, fermi_level_offset

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_and_continuation', &
         routineP = moduleN//':'//routineN

      CHARACTER(5)                                       :: MO_number
      COMPLEX(KIND=dp)                                   :: func_val, im_unit, one, re_unit, rho1, &
                                                            zero
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: dLambda, dLambda_2, Lambda, &
                                                            Lambda_without_offset, vec_b_gw, &
                                                            vec_b_gw_copy
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: mat_A_gw, mat_B_gw
      INTEGER                                            :: handle4, ierr, iii, iiter, info, &
                                                            integ_range, jjj, jquad, kkk, &
                                                            n_level_gw_ref, num_var, output_unit, &
                                                            xpos
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      LOGICAL                                            :: could_exit
      REAL(KIND=dp) :: chi2, chi2_old, delta, deriv_val_real, e_fermi, gw_energ, Ldown, &
         level_energ_GW, Lup, range_step, ScalParam, sign_occ_virt, stat_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Lambda_Im, Lambda_Re, stat_errors, &
                                                            vec_N_gw, vec_omega_fit_gw_sign
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_N_gw

      output_unit = cp_logger_get_default_io_unit()

      im_unit = (0.0_dp, 1.0_dp)
      re_unit = (1.0_dp, 0.0_dp)

      num_var = 2*num_poles+1
      ALLOCATE (Lambda(num_var))
      Lambda = (0.0_dp, 0.0_dp)
      ALLOCATE (Lambda_without_offset(num_var))
      Lambda_without_offset = (0.0_dp, 0.0_dp)
      ALLOCATE (Lambda_Re(num_var))
      Lambda_Re = 0.0_dp
      ALLOCATE (Lambda_Im(num_var))
      Lambda_Im = 0.0_dp

      ALLOCATE (vec_omega_fit_gw_sign(num_fit_points))

      IF (n_level_gw <= gw_corr_lev_occ) THEN
         sign_occ_virt = -1.0_dp
      ELSE
         sign_occ_virt = 1.0_dp
      END IF

      DO jquad = 1, num_fit_points
         vec_omega_fit_gw_sign(jquad) = ABS(vec_omega_fit_gw(jquad))*sign_occ_virt
      END DO

      ! initial guess
      range_step = (vec_omega_fit_gw_sign(num_fit_points)-vec_omega_fit_gw_sign(1))/(num_poles-1)
      DO iii = 1, num_poles
         Lambda_Im(2*iii+1) = vec_omega_fit_gw_sign(1)+(iii-1)*range_step
      END DO
      range_step = (vec_omega_fit_gw_sign(num_fit_points)-vec_omega_fit_gw_sign(1))/num_poles
      DO iii = 1, num_poles
         Lambda_Re(2*iii+1) = ABS(vec_omega_fit_gw_sign(1)+(iii-0.5_dp)*range_step)
      END DO

      DO iii = 1, num_var
         Lambda(iii) = Lambda_Re(iii)+im_unit*Lambda_Im(iii)
      END DO

      CALL calc_chi2(chi2_old, Lambda, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                     num_fit_points, n_level_gw)

      ALLOCATE (mat_A_gw(num_poles+1, num_poles+1))
      ALLOCATE (vec_b_gw(num_poles+1))
      ALLOCATE (ipiv(num_poles+1))
      mat_A_gw = (0.0_dp, 0.0_dp)
      vec_b_gw = 0.0_dp

      DO iii = 1, num_poles+1
         mat_A_gw(iii, 1) = (1.0_dp, 0.0_dp)
      END DO
      integ_range = num_fit_points/num_poles
      DO kkk = 1, num_poles+1
         xpos = (kkk-1)*integ_range+1
         xpos = MIN(xpos, num_fit_points)
         ! calculate coefficient at this point
         DO iii = 1, num_poles
            jjj = iii*2
            func_val = (1.0_dp, 0.0_dp)/(im_unit*vec_omega_fit_gw_sign(xpos)- &
                                         CMPLX(Lambda_Re(jjj+1), Lambda_Im(jjj+1), KIND=dp))
            mat_A_gw(kkk, iii+1) = func_val
         END DO
         vec_b_gw(kkk) = vec_Sigma_c_gw(n_level_gw, xpos)
      END DO

      ! Solve system of linear equations
      CALL ZGETRF(num_poles+1, num_poles+1, mat_A_gw, num_poles+1, ipiv, info)

      CALL ZGETRS('N', num_poles+1, 1, mat_A_gw, num_poles+1, ipiv, vec_b_gw, num_poles+1, info)

      Lambda_Re(1) = REAL(vec_b_gw(1))
      Lambda_Im(1) = AIMAG(vec_b_gw(1))
      DO iii = 1, num_poles
         jjj = iii*2
         Lambda_Re(jjj) = REAL(vec_b_gw(iii+1))
         Lambda_Im(jjj) = AIMAG(vec_b_gw(iii+1))
      END DO

      DEALLOCATE (mat_A_gw)
      DEALLOCATE (vec_b_gw)
      DEALLOCATE (ipiv)

      ALLOCATE (mat_A_gw(num_var*2, num_var*2))
      ALLOCATE (mat_B_gw(num_fit_points, num_var*2))
      ALLOCATE (dLambda(num_fit_points))
      ALLOCATE (dLambda_2(num_fit_points))
      ALLOCATE (vec_b_gw(num_var*2))
      ALLOCATE (vec_b_gw_copy(num_var*2))
      ALLOCATE (ipiv(num_var*2))

      ScalParam = 0.01_dp
      Ldown = 1.5_dp
      Lup = 10.0_dp
      could_exit = .FALSE.

      ! iteration loop for fitting
      DO iiter = 1, max_iter_fit

         CALL timeset(routineN//"_fit_loop_1", handle4)

         ! calc delta lambda
         DO iii = 1, num_var
            Lambda(iii) = Lambda_Re(iii)+im_unit*Lambda_Im(iii)
         END DO
         dLambda = (0.0_dp, 0.0_dp)

         DO kkk = 1, num_fit_points
            func_val = Lambda(1)
            DO iii = 1, num_poles
               jjj = iii*2
               func_val = func_val+Lambda(jjj)/(vec_omega_fit_gw_sign(kkk)*im_unit-Lambda(jjj+1))
            END DO
            dLambda(kkk) = vec_Sigma_c_gw(n_level_gw, kkk)-func_val
         END DO
         rho1 = SUM(dLambda*dLambda)

         ! fill matrix
         mat_B_gw = (0.0_dp, 0.0_dp)
         DO iii = 1, num_fit_points
            mat_B_gw(iii, 1) = 1.0_dp
            mat_B_gw(iii, num_var+1) = im_unit
         END DO
         DO iii = 1, num_poles
            jjj = iii*2
            DO kkk = 1, num_fit_points
               mat_B_gw(kkk, jjj) = 1.0_dp/(im_unit*vec_omega_fit_gw_sign(kkk)-Lambda(jjj+1))
               mat_B_gw(kkk, jjj+num_var) = im_unit/(im_unit*vec_omega_fit_gw_sign(kkk)-Lambda(jjj+1))
               mat_B_gw(kkk, jjj+1) = Lambda(jjj)/(im_unit*vec_omega_fit_gw_sign(kkk)-Lambda(jjj+1))**2
               mat_B_gw(kkk, jjj+1+num_var) = (-Lambda_Im(jjj)+im_unit*Lambda_Re(jjj))/ &
                                              (im_unit*vec_omega_fit_gw_sign(kkk)-Lambda(jjj+1))**2
            END DO
         END DO

         CALL timestop(handle4)

         CALL timeset(routineN//"_fit_matmul_1", handle4)

         one = (1.0_dp, 0.0_dp)
         zero = (0.0_dp, 0.0_dp)
         CALL zgemm('C', 'N', num_var*2, num_var*2, num_fit_points, one, mat_B_gw, num_fit_points, mat_B_gw, num_fit_points, &
                    zero, mat_A_gw, num_var*2)
         CALL timestop(handle4)

         CALL timeset(routineN//"_fit_zgemv_1", handle4)
         CALL zgemv('C', num_fit_points, num_var*2, one, mat_B_gw, num_fit_points, dLambda, 1, &
                    zero, vec_b_gw, 1)

         CALL timestop(handle4)

         ! scale diagonal elements of a_mat
         DO iii = 1, num_var*2
            mat_A_gw(iii, iii) = mat_A_gw(iii, iii)+ScalParam*mat_A_gw(iii, iii)
         END DO

         ! solve linear system
         ierr = 0
         ipiv = 0

         CALL timeset(routineN//"_fit_lin_eq_2", handle4)

         CALL ZGETRF(2*num_var, 2*num_var, mat_A_gw, 2*num_var, ipiv, info)

         CALL ZGETRS('N', 2*num_var, 1, mat_A_gw, 2*num_var, ipiv, vec_b_gw, 2*num_var, info)

         CALL timestop(handle4)

         DO iii = 1, num_var
            Lambda(iii) = Lambda_Re(iii)+im_unit*Lambda_Im(iii)+vec_b_gw(iii)+vec_b_gw(iii+num_var)
         END DO

         ! calculate chi2
         CALL calc_chi2(chi2, Lambda, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                        num_fit_points, n_level_gw)

         IF (chi2 < chi2_old) THEN
            ScalParam = MAX(ScalParam/Ldown, 1E-12_dp)
            DO iii = 1, num_var
               Lambda_Re(iii) = Lambda_Re(iii)+REAL(vec_b_gw(iii)+vec_b_gw(iii+num_var))
               Lambda_Im(iii) = Lambda_Im(iii)+AIMAG(vec_b_gw(iii)+vec_b_gw(iii+num_var))
            END DO
            IF (chi2_old/chi2-1.0_dp < stop_crit) could_exit = .TRUE.
            chi2_old = chi2
         ELSE
            ScalParam = ScalParam*Lup
         END IF
         IF (ScalParam > 100.0_dp .AND. could_exit) EXIT

         IF (ScalParam > 1E+10_dp) ScalParam = 1E-4_dp

         n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
         IF (iiter == max_iter_fit) THEN
            WRITE (MO_number, "(I3)") n_level_gw_ref
            CALL cp_warn(__LOCATION__, &
                         "The fit for corrected level n ="//MO_number//" did not converge. "// &
                         "For levels close to HOMO or LUMO, this is normally no issue. "// &
                         "To avoid this warning, you can (1) increase the "// &
                         "number of fit iterations MAX_ITER_FIT, or you can (2) increase the number "// &
                         "of RPA integration points (then, Sigma_c(i*omega) is more accurate) "// &
                         "or  you can (3) decrease "// &
                         "the fit range by setting the keyword OMEGA_MAX_FIT (in Hartree). ")

         END IF

      END DO

      ! change a_0 [Lambda(1)], so that Sigma(i0) = Fit(i0)
      func_val = Lambda(1)
      DO iii = 1, num_poles
         jjj = iii*2
         ! calculate value of the fit function
         func_val = func_val+Lambda(jjj)/(-Lambda(jjj+1))
      END DO

      Lambda_Re(1) = Lambda_Re(1)-REAL(func_val)+REAL(vec_Sigma_c_gw(n_level_gw, num_fit_points))
      Lambda_Im(1) = Lambda_Im(1)-AIMAG(func_val)+AIMAG(vec_Sigma_c_gw(n_level_gw, num_fit_points))

      Lambda_without_offset(:) = Lambda(:)

      DO iii = 1, num_var
         Lambda(iii) = CMPLX(Lambda_Re(iii), Lambda_Im(iii), KIND=dp)
      END DO

      ! print self_energy and fit on the imaginary frequency axis if required
      IF (check_fit) THEN

         IF (output_unit > 0) THEN

            WRITE (output_unit, *) ' '
            WRITE (output_unit, '(T3,A,I5)') 'Check the GW fit for molecular orbital', n_level_gw_ref
            WRITE (output_unit, '(T3,A)') '-------------------------------------------'
            WRITE (output_unit, *)
            WRITE (output_unit, '(T3,5A)') '  omega (i*eV)    ', 'Re(fit) (eV)    ', &
               'Im(fit) (eV)  ', 'Re(Sig_c) (eV)  ', &
               'Im(Sig_c) (eV)'

         END IF

         DO kkk = 1, num_fit_points
            func_val = Lambda(1)
            DO iii = 1, num_poles
               jjj = iii*2
               ! calculate value of the fit function
               func_val = func_val+Lambda(jjj)/(im_unit*vec_omega_fit_gw_sign(kkk)-Lambda(jjj+1))
            END DO
            WRITE (output_unit, '(1F16.3,4F16.5)') vec_omega_fit_gw_sign(kkk)*27.211_dp, REAL(func_val)*27.211_dp, &
               AIMAG(func_val)*27.211_dp, REAL(vec_Sigma_c_gw(n_level_gw, kkk))*27.211_dp, &
               AIMAG(vec_Sigma_c_gw(n_level_gw, kkk))*27.211_dp
         END DO

         WRITE (output_unit, *) ' '

      END IF

      ! correct the electronic levels
      IF (n_level_gw <= gw_corr_lev_occ) THEN
         e_fermi = Eigenval(homo)+fermi_level_offset
      ELSE
         e_fermi = Eigenval(homo+1)-fermi_level_offset
      END IF

      ! either Z-shot or no crossing search for evaluating Sigma_c
      IF (crossing_search == ri_rpa_g0w0_crossing_none) THEN

         ! calculate func val on the real axis
         ! gw_energ = only correlation part of the self energy
         func_val = Lambda(1)
         DO iii = 1, num_poles
            jjj = iii*2
            func_val = func_val+Lambda(jjj)/(Eigenval(n_level_gw_ref)-e_fermi-Lambda(jjj+1))
         END DO

         gw_energ = REAL(func_val)
         vec_gw_energ(n_level_gw) = gw_energ

      ELSE IF (crossing_search == ri_rpa_g0w0_crossing_z_shot .OR. &
               crossing_search == ri_rpa_g0w0_crossing_newton) THEN

         ! calculate Sigma_c_fit(e_n) and Z
         func_val = Lambda(1)
         z_value(n_level_gw) = 1.0_dp
         DO iii = 1, num_poles
            jjj = iii*2
            z_value(n_level_gw) = z_value(n_level_gw)+REAL(Lambda(jjj)/ &
                                                           (Eigenval(n_level_gw_ref)-e_fermi-Lambda(jjj+1))**2)
            func_val = func_val+Lambda(jjj)/(Eigenval(n_level_gw_ref)-e_fermi-Lambda(jjj+1))
         END DO
         ! m is the slope of the correl self-energy
         m_value(n_level_gw) = 1.0_dp-z_value(n_level_gw)
         z_value(n_level_gw) = 1.0_dp/z_value(n_level_gw)
         gw_energ = REAL(func_val)
         vec_gw_energ(n_level_gw) = gw_energ

         ! in case one wants to do Newton-Raphson on top of the Z-shot
         IF (crossing_search == ri_rpa_g0w0_crossing_newton) THEN

            level_energ_GW = (Eigenval_scf(n_level_gw_ref)- &
                              m_value(n_level_gw)*Eigenval(n_level_gw_ref)+ &
                              vec_gw_energ(n_level_gw)+ &
                              vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))* &
                             z_value(n_level_gw)

            ! Newton-Raphson iteration
            DO kkk = 1, 1000

               ! calculate the value of the fit function for level_energ_GW
               func_val = Lambda(1)
               z_value(n_level_gw) = 1.0_dp
               DO iii = 1, num_poles
                  jjj = iii*2
                  func_val = func_val+Lambda(jjj)/(level_energ_GW-e_fermi-Lambda(jjj+1))
               END DO

               ! calculate the derivative of the fit function for level_energ_GW
               deriv_val_real = -1.0_dp
               DO iii = 1, num_poles
                  jjj = iii*2
                  deriv_val_real = deriv_val_real+REAL(Lambda(jjj))/((ABS(level_energ_GW-e_fermi-Lambda(jjj+1)))**2) &
                                   -(REAL(Lambda(jjj))*(level_energ_GW-e_fermi)-REAL(Lambda(jjj)*CONJG(Lambda(jjj+1))))* &
                                   2.0_dp*(level_energ_GW-e_fermi-REAL(Lambda(jjj+1)))/ &
                                   ((ABS(level_energ_GW-e_fermi-Lambda(jjj+1)))**2)

               END DO

               delta = (Eigenval_scf(n_level_gw_ref)+vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)+REAL(func_val)-level_energ_GW)/ &
                       deriv_val_real

               level_energ_GW = level_energ_GW-delta

               IF (ABS(delta) < 1.0E-08) EXIT

            END DO

            ! update the GW-energy by Newton-Raphson and set the Z-value to 1

            vec_gw_energ(n_level_gw) = REAL(func_val)
            z_value(n_level_gw) = 1.0_dp
            m_value(n_level_gw) = 0.0_dp

         END IF ! Newton-Raphson on top of Z-shot

      END IF ! decision crossing search none, Z-shot

      !   --------------------------------------------
      !  | calculate statistical error due to fitting |
      !   --------------------------------------------

      ! estimate the statistical error of the calculated Sigma_c(i*omega)
      ! by sqrt(chi2/n), where n is the number of fit points

      CALL calc_chi2(chi2, Lambda_without_offset, vec_Sigma_c_gw, vec_omega_fit_gw_sign, num_poles, &
                     num_fit_points, n_level_gw)

      ! Estimate the statistical error of every fit point
      stat_error = SQRT(chi2/num_fit_points)

      ! allocate N array containing the second derivatives of chi^2
      ALLOCATE (vec_N_gw(num_var*2))
      vec_N_gw = 0.0_dp

      ALLOCATE (mat_N_gw(num_var*2, num_var*2))
      mat_N_gw = 0.0_dp

      DO iii = 1, num_var*2
         CALL calc_mat_N(vec_N_gw(iii), Lambda_without_offset, vec_Sigma_c_gw, vec_omega_fit_gw_sign, &
                         iii, iii, num_poles, num_fit_points, n_level_gw, 0.001_dp)
      END DO

      DO iii = 1, num_var*2
         DO jjj = 1, num_var*2
            CALL calc_mat_N(mat_N_gw(iii, jjj), Lambda_without_offset, vec_Sigma_c_gw, vec_omega_fit_gw_sign, &
                            iii, jjj, num_poles, num_fit_points, n_level_gw, 0.001_dp)
         END DO
      END DO

      CALL DGETRF(2*num_var, 2*num_var, mat_N_gw, 2*num_var, ipiv, info)

      ! vec_b_gw is only working array
      CALL DGETRI(2*num_var, mat_N_gw, 2*num_var, ipiv, vec_b_gw, 2*num_var, info)

      ALLOCATE (stat_errors(2*num_var))
      stat_errors = 0.0_dp

      DO iii = 1, 2*num_var
         stat_errors(iii) = SQRT(ABS(mat_N_gw(iii, iii)))*stat_error
      END DO

      ! Compute error of Sigma_c on real axis according to error propagation

      vec_gw_energ_error_fit(n_level_gw) = 0.0_dp

      DO kkk = 1, num_poles
         vec_gw_energ_error_fit(n_level_gw) = vec_gw_energ_error_fit(n_level_gw)+ &
                                              (stat_errors(4*kkk-1)+stat_errors(4*kkk))* &
                                              ABS(1.0_dp/(Eigenval(n_level_gw_ref)-e_fermi-Lambda(2*kkk+1))- &
                                                  1.0_dp/(-Lambda(2*kkk+1)))+ &
                                              (stat_errors(4*kkk+1)+stat_errors(4*kkk+2))*ABS(Lambda(2*kkk))* &
                                              ABS(1.0_dp/(Eigenval(n_level_gw_ref)-e_fermi-Lambda(2*kkk+1))**2- &
                                                  1.0_dp/(-Lambda(2*kkk+1))**2)
      END DO

      DEALLOCATE (mat_N_gw)
      DEALLOCATE (vec_N_gw)
      DEALLOCATE (mat_A_gw)
      DEALLOCATE (mat_B_gw)
      DEALLOCATE (stat_errors)
      DEALLOCATE (dLambda)
      DEALLOCATE (dLambda_2)
      DEALLOCATE (vec_b_gw)
      DEALLOCATE (vec_b_gw_copy)
      DEALLOCATE (ipiv)
      DEALLOCATE (vec_omega_fit_gw_sign)
      DEALLOCATE (Lambda)
      DEALLOCATE (Lambda_without_offset)
      DEALLOCATE (Lambda_Re)
      DEALLOCATE (Lambda_Im)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param vec_gw_energ ...
!> \param z_value ...
!> \param m_value ...
!> \param vec_Sigma_c_gw_real_freq ...
!> \param contour_def_offset ...
!> \param n_level_gw ...
!> \param index_contour_def ...
! **************************************************************************************************
   SUBROUTINE compute_z_and_m_contour_def(vec_gw_energ, z_value, m_value, &
                                          vec_Sigma_c_gw_real_freq, contour_def_offset, &
                                          n_level_gw, index_contour_def)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_gw_energ, z_value, m_value, &
                                                            vec_Sigma_c_gw_real_freq
      REAL(KIND=dp)                                      :: contour_def_offset
      INTEGER                                            :: n_level_gw, index_contour_def

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_z_and_m_contour_def', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: derivative

! average of Sigma_c(e_n+delta) and Sigma_c(e_n-delta)

      vec_gw_energ(n_level_gw) = 0.5_dp*(vec_Sigma_c_gw_real_freq(2*index_contour_def-1) &
                                         +vec_Sigma_c_gw_real_freq(2*index_contour_def))

      derivative = 0.5_dp*(vec_Sigma_c_gw_real_freq(2*index_contour_def) &
                           -vec_Sigma_c_gw_real_freq(2*index_contour_def-1)) &
                   /contour_def_offset

      z_value(n_level_gw) = 1.0_dp/(1.0_dp-derivative)

      m_value(n_level_gw) = derivative

   END SUBROUTINE

! **************************************************************************************************
!> \brief Evaluate the fit function at imaginary time
!> \param func_val ...
!> \param tau ...
!> \param num_poles ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE eval_fit_func_im_time(func_val, tau, num_poles, beta)

      REAL(KIND=dp)                                      :: func_val, tau
      INTEGER                                            :: num_poles
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii, jjj

      func_val = 0.0_dp
      DO iii = 1, num_poles
         jjj = iii*3
         ! calculate value of the fit function
         IF (tau > 0) THEN
            func_val = func_val+beta(jjj-2)*EXP(-tau*beta(jjj))
         ELSE
            func_val = func_val+beta(jjj-1)*EXP(tau*beta(jjj))
         END IF
      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tj ...
!> \param wj ...
! **************************************************************************************************
   SUBROUTINE eval_fit_func_tau_grid(func_val, x_value, num_integ_points, &
                                     tj, wj)

      REAL(KIND=dp)                                      :: func_val, x_value
      INTEGER                                            :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tj, wj

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_tau_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val+wj(iii)*EXP(-2.0_dp*x_value*tj(iii)**2)

      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief Evaluate fit function when calculating tau grid for Fourier transform
!> \param func_val ...
!> \param x_value ...
!> \param num_integ_points ...
!> \param tau_tj ...
!> \param tau_wj_work ...
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE eval_fit_func_tau_grid_cosine(func_val, x_value, num_integ_points, tau_tj, tau_wj_work, omega)

      REAL(KIND=dp)                                      :: func_val, x_value
      INTEGER                                            :: num_integ_points
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tau_wj_work
      REAL(KIND=dp)                                      :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'eval_fit_func_tau_grid_cosine', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iii

      func_val = 0.0_dp

      DO iii = 1, num_integ_points

         ! calculate value of the fit function
         func_val = func_val+tau_wj_work(iii)*COS(omega*tau_tj(iii))*EXP(-x_value*tau_tj(iii))

      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief Prints the GW stuff to the output and optinally to an external file.
!>        Also updates the eigenvalues for eigenvalue-self-consistent GW
!> \param vec_gw_energ ...
!> \param vec_gw_energ_error_fit ...
!> \param z_value ...
!> \param m_value ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param Eigenval ...
!> \param Eigenval_last ...
!> \param Eigenval_scf ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param gw_corr_lev_tot ...
!> \param count_ev_sc_GW ...
!> \param crossing_search ...
!> \param homo ...
!> \param nmo ...
!> \param unit_nr ...
!> \param print_gw_details ...
!> \param do_alpha ...
!> \param do_beta ...
! **************************************************************************************************
   SUBROUTINE print_and_update_for_ev_sc(vec_gw_energ, vec_gw_energ_error_fit, &
                                         z_value, m_value, vec_Sigma_x_minus_vxc_gw, Eigenval, &
                                         Eigenval_last, Eigenval_scf, gw_corr_lev_occ, gw_corr_lev_virt, gw_corr_lev_tot, &
                                         count_ev_sc_GW, crossing_search, homo, nmo, unit_nr, print_gw_details, &
                                         do_alpha, do_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_gw_energ, vec_gw_energ_error_fit, &
                                                            z_value, m_value
      REAL(KIND=dp), DIMENSION(:)                        :: vec_Sigma_x_minus_vxc_gw, Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval_last, Eigenval_scf
      INTEGER                                            :: gw_corr_lev_occ, gw_corr_lev_virt, &
                                                            gw_corr_lev_tot, count_ev_sc_GW, &
                                                            crossing_search, homo, nmo, unit_nr
      LOGICAL                                            :: print_gw_details
      LOGICAL, OPTIONAL                                  :: do_alpha, do_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_and_update_for_ev_sc', &
         routineP = moduleN//':'//routineN

      CHARACTER(4)                                       :: occ_virt
      INTEGER                                            :: handle, n_level_gw, n_level_gw_ref
      LOGICAL                                            :: do_closed_shell, my_do_alpha, my_do_beta
      REAL(KIND=dp)                                      :: eigen_diff

      CALL timeset(routineN, handle)

      IF (PRESENT(do_alpha)) THEN
         my_do_alpha = do_alpha
      ELSE
         my_do_alpha = .FALSE.
      END IF

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

      Eigenval_last(:) = Eigenval(:)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ' '

         IF (do_closed_shell) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies'
            WRITE (unit_nr, '(T3,A)') '-------------------------'
         ELSE IF (my_do_alpha) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of alpha spins'
            WRITE (unit_nr, '(T3,A)') '----------------------------------------'
         ELSE IF (my_do_beta) THEN
            WRITE (unit_nr, '(T3,A)') 'GW quasiparticle energies of beta spins'
            WRITE (unit_nr, '(T3,A)') '---------------------------------------'
         END IF

      END IF

      IF (unit_nr > 0 .AND. (.NOT. print_gw_details)) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T5,A)') 'Molecular orbital        MO energy after SCF (eV)        G0W0 QP energy (eV)'
      END IF

      IF (unit_nr > 0 .AND. print_gw_details) THEN
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') 'The GW quasiparticle energies are calculated according to: '
      END IF

      IF (crossing_search == ri_rpa_g0w0_crossing_none) THEN

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
            Eigenval(n_level_gw_ref) = Eigenval_scf(n_level_gw_ref)+vec_gw_energ(n_level_gw)+ &
                                       vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)
         END DO

         IF (unit_nr > 0 .AND. print_gw_details) THEN
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'E_GW = E_SCF +  Sigc(E_SCF) + Sigx - vxc'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'The energy unit of the following table is eV. Sigc_fit is a very conservative'
            WRITE (unit_nr, '(T3,A)') 'estimate of the statistical error of the correlation self-energy caused by the'
            WRITE (unit_nr, '(T3,A)') 'fitting.'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T14,2A)') 'MO        E_SCF         Sigc', &
               '     Sigc_fit     Sigx-vxc         E_GW'
         END IF

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            IF (unit_nr > 0 .AND. (.NOT. print_gw_details)) THEN
               WRITE (unit_nr, '(T5,I9,3A,2F27.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')     ', &
                  Eigenval_last(n_level_gw_ref)*27.211_dp, &
                  Eigenval(n_level_gw_ref)*27.211_dp
            END IF

            IF (unit_nr > 0 .AND. print_gw_details) THEN
               WRITE (unit_nr, '(T4,I4,3A,5F13.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')', &
                  Eigenval_last(n_level_gw_ref)*27.211_dp, &
                  vec_gw_energ(n_level_gw)*27.211_dp, &
                  vec_gw_energ_error_fit(n_level_gw)*27.211_dp, &
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)*27.211_dp, &
                  Eigenval(n_level_gw_ref)*27.211_dp
            END IF
         END DO

         ! z-shot
      ELSE

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
            Eigenval(n_level_gw_ref) = (Eigenval_scf(n_level_gw_ref)- &
                                        m_value(n_level_gw)*Eigenval(n_level_gw_ref)+ &
                                        vec_gw_energ(n_level_gw)+ &
                                        vec_Sigma_x_minus_vxc_gw(n_level_gw_ref))* &
                                       z_value(n_level_gw)
         END DO

         IF (unit_nr > 0 .AND. print_gw_details) THEN
            WRITE (unit_nr, '(T3,A)') 'E_GW = E_SCF + Z * ( Sigc(E_SCF) + Sigx - vxc )'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T3,A)') 'The energy unit of the following table is eV.  Sigc_fit is a very conservative'
            WRITE (unit_nr, '(T3,2A)') 'estimate of the statistical error of the fitting.'
            WRITE (unit_nr, '(T3,A)') ' '
            WRITE (unit_nr, '(T13,2A)') 'MO      E_SCF       Sigc', &
               '   Sigc_fit   Sigx-vxc          Z       E_GW'
         END IF

         DO n_level_gw = 1, gw_corr_lev_tot
            n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
            IF (n_level_gw <= gw_corr_lev_occ) THEN
               occ_virt = 'occ'
            ELSE
               occ_virt = 'vir'
            END IF

            IF (unit_nr > 0 .AND. (.NOT. print_gw_details)) THEN
               WRITE (unit_nr, '(T5,I9,3A,2F27.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')     ', &
                  Eigenval_last(n_level_gw_ref)*27.211_dp, &
                  Eigenval(n_level_gw_ref)*27.211_dp
            END IF

            IF (unit_nr > 0 .AND. print_gw_details) THEN
               WRITE (unit_nr, '(T3,I4,3A,6F11.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')', &
                  Eigenval_last(n_level_gw_ref)*27.211_dp, &
                  vec_gw_energ(n_level_gw)*27.211_dp, &
                  vec_gw_energ_error_fit(n_level_gw)*27.211_dp, &
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref)*27.211_dp, &
                  z_value(n_level_gw), &
                  Eigenval(n_level_gw_ref)*27.211_dp
            END IF
         END DO

      END IF ! z-shot vs. no crossing

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
      END IF

      ! for eigenvalue self-consistent GW, all eigenvalues have to be corrected
      ! 1) the occupied; check if there are occupied MOs not being corrected by GW
      IF (gw_corr_lev_occ < homo .AND. gw_corr_lev_occ > 0) THEN

         ! calculate average GW correction for occupied orbitals
         eigen_diff = 0.0_dp

         DO n_level_gw = 1, gw_corr_lev_occ
            n_level_gw_ref = n_level_gw+homo-gw_corr_lev_occ
            eigen_diff = eigen_diff+Eigenval(n_level_gw_ref)-Eigenval_last(n_level_gw_ref)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_occ

         ! correct the eigenvalues of the occupied orbitals which have not been corrected by GW
         DO n_level_gw = 1, homo-gw_corr_lev_occ
            Eigenval(n_level_gw) = Eigenval(n_level_gw)+eigen_diff
         END DO

      END IF

      ! 2) the virtual: check if there are virtual orbitals not being corrected by GW
      IF (gw_corr_lev_virt < nmo-homo .AND. gw_corr_lev_virt > 0) THEN

         ! calculate average GW correction for virtual orbitals
         eigen_diff = 0.0_dp
         DO n_level_gw = 1, gw_corr_lev_virt
            n_level_gw_ref = n_level_gw+homo
            eigen_diff = eigen_diff+Eigenval(n_level_gw_ref)-Eigenval_last(n_level_gw_ref)
         END DO
         eigen_diff = eigen_diff/gw_corr_lev_virt

         ! correct the eigenvalues of the virtual orbitals which have not been corrected by GW
         DO n_level_gw = homo+gw_corr_lev_virt+1, nmo
            Eigenval(n_level_gw) = Eigenval(n_level_gw)+eigen_diff
         END DO

      END IF

      IF ((gw_corr_lev_occ == 0 .OR. gw_corr_lev_virt == 0) .AND. count_ev_sc_GW > 1) THEN
         CALL cp_warn(__LOCATION__, &
                      "Please increase for eigenvalue-self-consistent GW, the number of "// &
                      "corrected occupied and/or virtual orbitals above 0.")
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_and_update_for_ev_sc

! **************************************************************************************************
!> \brief Do the numerical integration for the contour deformation technique:
!>        First, numerical integration on the imag. freq. axis, then residues
!> \param vec_Sigma_c_gw_real_freq ...
!> \param wj ...
!> \param vec_W_gw ...
!> \param Eigenval ...
!> \param omega ...
!> \param contour_def_offset ...
!> \param e_fermi ...
!> \param mo_from_jquad ...
!> \param residue_from_jquad ...
!> \param alpha_beta_from_jquad ...
!> \param contour_def_start ...
!> \param contour_def_end ...
!> \param n_global ...
!> \param jquad ...
!> \param num_integ_points ...
!> \param m_global ...
!> \param nm_global ...
!> \param alpha_beta_case ...
!> \param alpha ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE contour_def_integrate_and_residues(vec_Sigma_c_gw_real_freq, wj, vec_W_gw, Eigenval, &
                                                 omega, contour_def_offset, e_fermi, &
                                                 mo_from_jquad, residue_from_jquad, alpha_beta_from_jquad, &
                                                 contour_def_start, contour_def_end, n_global, jquad, &
                                                 num_integ_points, m_global, nm_global, &
                                                 alpha_beta_case, alpha, beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: vec_Sigma_c_gw_real_freq, wj, vec_W_gw
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      REAL(KIND=dp)                                      :: omega, contour_def_offset, e_fermi
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mo_from_jquad, residue_from_jquad, &
                                                            alpha_beta_from_jquad
      INTEGER                                            :: contour_def_start, contour_def_end, &
                                                            n_global, jquad, num_integ_points, &
                                                            m_global, nm_global
      LOGICAL, OPTIONAL                                  :: alpha_beta_case, alpha, beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contour_def_integrate_and_residues', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, index_contour_def, &
                                                            m_level_of_real_energy, &
                                                            n_level_of_real_energy, spin, &
                                                            spin_of_residue
      LOGICAL                                            :: my_alpha, my_alpha_beta_case, my_beta
      REAL(KIND=dp)                                      :: real_energy

      CALL timeset(routineN, handle)

      my_alpha_beta_case = .FALSE.
      IF (PRESENT(alpha_beta_case)) my_alpha_beta_case = alpha_beta_case

      my_alpha = .FALSE.
      IF (PRESENT(alpha)) my_alpha = alpha

      my_beta = .FALSE.
      IF (PRESENT(beta)) my_beta = beta

      spin = 1
      IF (my_alpha_beta_case) THEN
         IF (my_alpha) spin = 1
         IF (my_beta) spin = 2
      END IF

      ! check whether we apply contour deformation to the level n_global
      IF (contour_def_start <= n_global .AND. n_global <= contour_def_end) THEN

         ! check whether we do integration over imag. frequencies or we compute the residues
         IF (jquad <= num_integ_points) THEN

            ! do the integration on the imaginary frequency axis, first for energy slightly below e_n
            real_energy = Eigenval(n_global)-contour_def_offset

            index_contour_def = 2*(n_global-contour_def_start+1)-1

            CALL numerical_integrate_contour_def(vec_Sigma_c_gw_real_freq(index_contour_def), vec_W_gw(nm_global), &
                                                 wj(jquad), omega, real_energy, Eigenval(m_global))

            ! now, for energy slightly above e_n
            real_energy = Eigenval(n_global)+contour_def_offset

            index_contour_def = 2*(n_global-contour_def_start+1)

            CALL numerical_integrate_contour_def(vec_Sigma_c_gw_real_freq(index_contour_def), vec_W_gw(nm_global), &
                                                 wj(jquad), omega, real_energy, Eigenval(m_global))

            ! residues
         ELSE

            ! for jquad > num_integ_points, we add the residue, but only for a single MO n and for a single residue m
            n_level_of_real_energy = mo_from_jquad(jquad-num_integ_points)
            m_level_of_real_energy = residue_from_jquad(jquad-num_integ_points)

            spin_of_residue = 1
            IF (my_alpha_beta_case) THEN
               spin_of_residue = alpha_beta_from_jquad(jquad-num_integ_points)
            END IF

            ! only add residue to a specific level
            IF (n_level_of_real_energy == n_global .AND. &
                m_level_of_real_energy == m_global .AND. &
                spin == spin_of_residue) THEN

               ! add residue, first for energy slightly below e_n
               index_contour_def = 2*(n_global-contour_def_start+1)-1

               IF (e_fermi-Eigenval(n_global) < omega .AND. omega < -contour_def_offset) THEN
                  ! update for residue
                  vec_Sigma_c_gw_real_freq(index_contour_def) = &
                     vec_Sigma_c_gw_real_freq(index_contour_def) &
                     +vec_W_gw(nm_global)
               ELSE IF (-contour_def_offset < omega .AND. omega < e_fermi-Eigenval(n_global)) THEN
                  vec_Sigma_c_gw_real_freq(index_contour_def) = &
                     vec_Sigma_c_gw_real_freq(index_contour_def) &
                     -vec_W_gw(nm_global)
               END IF

               ! add residue, first for energy slightly above e_n
               index_contour_def = 2*(n_global-contour_def_start+1)

               IF (e_fermi-Eigenval(n_global) < omega .AND. omega < contour_def_offset) THEN
                  ! update for residue
                  vec_Sigma_c_gw_real_freq(index_contour_def) = &
                     vec_Sigma_c_gw_real_freq(index_contour_def) &
                     +vec_W_gw(nm_global)
               ELSE IF (contour_def_offset < omega .AND. omega < e_fermi-Eigenval(n_global)) THEN
                  vec_Sigma_c_gw_real_freq(index_contour_def) = &
                     vec_Sigma_c_gw_real_freq(index_contour_def) &
                     -vec_W_gw(nm_global)
               END IF

            END IF ! check whether we have the right residue

         END IF ! imag. frequency / real frequency

      END IF ! check contour deformation for n_global

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int ...
!> \param para_env_sub ...
!> \param cut_RI ...
!> \param non_zero_blocks_3c ...
! **************************************************************************************************
   SUBROUTINE get_non_zero_blocks_3c(mat_3c_overl_int, para_env_sub, cut_RI, non_zero_blocks_3c)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_3c_overl_int
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER                                            :: cut_RI
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_non_zero_blocks_3c', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, block_counter, col, handle, &
                                                            i_cut_RI, iblk, imepos, maxlength, &
                                                            maxlength_tmp, nblkrows_total, row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c_tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_3c_overl_int(1)%matrix, nblkrows_total=nblkrows_total)

      ALLOCATE (non_zero_blocks_3c_tmp(1:cut_RI, 1:nblkrows_total, 0:(para_env_sub%num_pe-1)))
      non_zero_blocks_3c_tmp = 0

      DO i_cut_RI = 1, cut_RI

         CALL dbcsr_iterator_start(iter, mat_3c_overl_int(i_cut_RI)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk)

            non_zero_blocks_3c_tmp(i_cut_RI, row, para_env_sub%mepos) = 1

         ENDDO

         CALL dbcsr_iterator_stop(iter)

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

   END SUBROUTINE

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

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_3c_overl_int_cut_col
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER                                            :: cut_RI, cut_memory
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: non_zero_blocks_3c_cut

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_non_zero_blocks_3c_cut_col', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, block_counter, col, handle, i_cut_RI, i_mem, iblk, imepos, maxlength, &
         maxlength_tmp, nblkrows_total, nblkrows_total_max, row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c_tmp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nblkrows_total_max = 0

      DO i_cut_RI = 1, cut_RI
         DO i_mem = 1, cut_memory
            CALL dbcsr_get_info(mat_3c_overl_int_cut_col(i_cut_RI, i_mem)%matrix, nblkrows_total=nblkrows_total)
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

            CALL dbcsr_iterator_start(iter, mat_3c_overl_int_cut_col(i_cut_RI, i_mem)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk)

               non_zero_blocks_3c_tmp(i_cut_RI, col, para_env_sub%mepos) = 1

            ENDDO

            CALL dbcsr_iterator_stop(iter)

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

            CALL dbcsr_iterator_start(iter, mat_3c_overl_int_cut_col(i_cut_RI, i_mem)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk)

               non_zero_blocks_3c_tmp(i_cut_RI, col, para_env_sub%mepos) = 1

            ENDDO

            CALL dbcsr_iterator_stop(iter)

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

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Sigma_c_to_update ...
!> \param W_gw_element ...
!> \param weight ...
!> \param omega ...
!> \param real_energy ...
!> \param eigenval ...
! **************************************************************************************************
   SUBROUTINE numerical_integrate_contour_def(Sigma_c_to_update, W_gw_element, &
                                              weight, omega, real_energy, eigenval)

      REAL(KIND=dp)                                      :: Sigma_c_to_update, W_gw_element, weight, &
                                                            omega, real_energy, eigenval

      CHARACTER(LEN=*), PARAMETER :: routineN = 'numerical_integrate_contour_def', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: im_unit
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      im_unit = (0.0_dp, 1.0_dp)

      Sigma_c_to_update = Sigma_c_to_update- &
                          REAL( &
                          0.5_dp/pi*weight/2.0_dp*W_gw_element/(im_unit*omega+real_energy-eigenval)+ &
                          0.5_dp/pi*weight/2.0_dp*W_gw_element/(-im_unit*omega+real_energy-eigenval) &
                          )

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Main driver for RPA-AXK energies
!> \param qs_env ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param dimen_RI ...
!> \param dimen_ia ...
!> \param para_env_sub ...
!> \param para_env_RPA ...
!> \param RPA_proc_map ...
!> \param eig ...
!> \param fm_mat_S ...
!> \param homo ...
!> \param virtual ...
!> \param omega ...
!> \param buffer_send_axk ...
!> \param buffer_rec_axk ...
!> \param number_of_send_axk ...
!> \param map_send_size_axk ...
!> \param map_rec_size_axk ...
!> \param local_size_source_axk ...
!> \param mp2_env ...
!> \param mat_munu ...
!> \param unit_nr ...
!> \param e_axk_corr ...          AXK energy correctrion for a quadrature point
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE compute_axk_ener(qs_env, fm_mat_Q, fm_mat_Q_gemm, dimen_RI, dimen_ia, &
                               para_env_sub, para_env_RPA, &
                               RPA_proc_map, eig, fm_mat_S, homo, virtual, omega, &
                               buffer_send_axk, buffer_rec_axk, &
                               number_of_send_axk, map_send_size_axk, map_rec_size_axk, &
                               local_size_source_axk, mp2_env, mat_munu, unit_nr, e_axk_corr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q, fm_mat_Q_gemm
      INTEGER                                            :: dimen_RI, dimen_ia
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub, para_env_RPA
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: RPA_proc_map
      REAL(KIND=dp), DIMENSION(:)                        :: Eig
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      INTEGER                                            :: homo, virtual
      REAL(KIND=dp)                                      :: omega
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_send_axk, buffer_rec_axk
      INTEGER                                            :: number_of_send_axk
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: map_send_size_axk, map_rec_size_axk
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_size_source_axk
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(dbcsr_p_type)                                 :: mat_munu
      INTEGER                                            :: unit_nr
      REAL(KIND=dp)                                      :: e_axk_corr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_axk_ener', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.0000001_dp

      INTEGER :: avirt, handle, i_global, iib, iocc, j_global, jjb, kkb, L_counter, &
         my_group_L_end, my_group_L_size, my_group_L_start, ncol_local, ngroup, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigen_diff, eps_filter, trace_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_mat_Gamma_3, fm_mat_Q_tmp, &
                                                            fm_mat_R_half, fm_mat_R_half_gemm, &
                                                            fm_mat_U
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_3, dbcsr_Gamma_inu_P, &
                                                            dbcsr_Gamma_munu_P
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v

      CALL timeset(routineN, handle)

      ! Eigenvalues
      ALLOCATE (eigenval(dimen_RI))
      eigenval = 0.0_dp
      ! create the R_half and U matrices with a different blacs env similar to Q
      ! and a tmp_Q needed for diagonalization

      NULLIFY (fm_mat_U)
      NULLIFY (fm_mat_Q_tmp)
      NULLIFY (fm_mat_R_half)
      NULLIFY (fm_mat_Q_tmp)
      NULLIFY (fm_struct)

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          matrix_struct=fm_struct)
      CALL cp_fm_create(fm_mat_U, fm_struct, name="fm_mat_U")
      CALL cp_fm_create(fm_mat_R_half, fm_struct, name="fm_mat_R_half")
      CALL cp_fm_create(fm_mat_Q_tmp, fm_struct, name="fm_mat_Q_tmp")
      CALL cp_fm_set_all(matrix=fm_mat_Q_tmp, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=fm_mat_U, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=fm_mat_R_half, alpha=0.0_dp)

      ! Copy Q to Q_tmp
      CALL cp_fm_to_fm(fm_mat_Q, fm_mat_Q_tmp)

      CALL cp_fm_get_info(matrix=fm_mat_Q, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local)

      CALL cp_fm_scale(0.50_dp, fm_mat_Q_tmp)
      ! Diagonalize Q
      CALL choose_eigv_solver(fm_mat_Q_tmp, fm_mat_U, eigenval)

      !Calculate diagonal matrix for R_half

      ! U*diag stored in U, whereas eigenvectors are in fm_mat_Q_tmp
      !CALL cp_fm_to_fm(fm_mat_Q_tmp, fm_mat_U)
      CALL cp_fm_to_fm(fm_mat_U, fm_mat_Q_tmp)

      ! Manipulate eigenvalues to get diagonal matrix
      CALL cp_fm_get_info(matrix=fm_mat_U, &
                          nrow_local=nrow_local, &
                          row_indices=row_indices)

      DO iib = 1, dimen_RI
         IF (ABS(eigenval(iib)) .GE. thresh) THEN
            eigenval(iib) = &
               SQRT((1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp+eigenval(iib)) &
                    -1.0_dp/(eigenval(iib)*(eigenval(iib)+1.0_dp)))
         ELSE
            eigenval(iib) = 0.707_dp
         ENDIF
      ENDDO

      CALL cp_fm_column_scale(fm_mat_U, eigenval)

      ! Release memory
      DEALLOCATE (eigenval)

      ! Get R_half by multiplication
      CALL cp_gemm(transa="N", transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_U, matrix_b=fm_mat_Q_tmp, beta=0.0_dp, &
                   matrix_c=fm_mat_R_half)

      ! get info of fm_mat_S and initialize Gamma_3
      NULLIFY (fm_mat_Gamma_3)
      CALL cp_fm_create(fm_mat_Gamma_3, fm_mat_S%matrix_struct)
      CALL cp_fm_set_all(matrix=fm_mat_Gamma_3, alpha=0.0_dp)
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! Update G with a new value of Omega: in practice, it is G*S

      ! Here eig are orbital energies, don't confuse with eigenval, which are eigenvalues of Q!

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP             SHARED(ncol_local,nrow_local,col_indices,row_indices,Eig,fm_mat_S,virtual,homo,omega)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)

            iocc = MAX(1, i_global-1)/virtual+1
            avirt = i_global-(iocc-1)*virtual
            eigen_diff = eig(avirt+homo)-eig(iocc)

            fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)* &
                                            SQRT(eigen_diff/(eigen_diff**2+omega**2))

         END DO
      END DO
!$OMP END PARALLEL DO

      ! Redistribute fm_mat_R_half for "rectangular" multiplication: ia*P P*P
      NULLIFY (fm_mat_R_half_gemm)
      CALL cp_fm_create(fm_mat_R_half_gemm, fm_mat_Q_gemm%matrix_struct)
      CALL cp_fm_set_all(matrix=fm_mat_R_half_gemm, alpha=0.0_dp)

      CALL fm_redistribute(fm_mat_R_half, fm_mat_R_half_gemm, RPA_proc_map, buffer_rec_axk, buffer_send_axk, &
                           number_of_send_axk, &
                           map_send_size_axk, map_rec_size_axk, local_size_source_axk, para_env_RPA)

      ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2) )
      CALL cp_gemm(transa="N", transb="N", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                   matrix_a=fm_mat_S, matrix_b=fm_mat_R_half_gemm, beta=0.0_dp, &
                   matrix_c=fm_mat_Gamma_3)

      ! Remove extra factor from S after the multiplication
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,iiB,iocc,avirt,eigen_diff,i_global,j_global) &
!$OMP                         SHARED(ncol_local,nrow_local,col_indices,row_indices,eig,fm_mat_S,virtual,homo,omega)
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)

            iocc = MAX(1, i_global-1)/virtual+1
            avirt = i_global-(iocc-1)*virtual
            eigen_diff = eig(avirt+homo)-eig(iocc)

            fm_mat_S%local_data(iiB, jjB) = fm_mat_S%local_data(iiB, jjB)/ &
                                            SQRT(eigen_diff/(eigen_diff**2+omega**2))

         END DO
      END DO
!$OMP              END PARALLEL DO

      ! Release full matrix stuff
      CALL cp_fm_release(fm_mat_Q_tmp)
      CALL cp_fm_release(fm_mat_U)
      CALL cp_fm_release(fm_mat_R_half)
      CALL cp_fm_release(fm_mat_R_half_gemm)

      ! Retrieve mo coefficients in dbcsr format
      NULLIFY (mo_coeff_o, mo_coeff_v)
      mo_coeff_o => mp2_env%ri_rpa%mo_coeff_o
      mo_coeff_v => mp2_env%ri_rpa%mo_coeff_v

      ! Get aux sizes
      ngroup = para_env_RPA%num_pe/para_env_sub%num_pe
      CALL get_aux_sizes(dimen_RI, ngroup, my_group_L_size, my_group_L_start, my_group_L_end, &
                         para_env_sub%num_pe, para_env_RPA%mepos)

      ! Copy Gamma_ia_P^3 to dbcsr matrix set
      CALL gamma_fm_to_dbcsr(fm_mat_Gamma_3, dbcsr_Gamma_3, para_env_RPA, para_env_sub, &
                             homo, virtual, mo_coeff_o, ngroup, my_group_L_start, &
                             my_group_L_end, my_group_L_size, dimen_RI)

      ! Release Gamma_3
      CALL cp_fm_release(fm_mat_Gamma_3)

      ! Create more dbcsr matrices

      NULLIFY (dbcsr_Gamma_inu_P)
      !CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_inu_P, ncol_local)
      CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_inu_P, my_group_L_size)
      NULLIFY (dbcsr_Gamma_munu_P)
      !CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_munu_P, ncol_local)
      CALL dbcsr_allocate_matrix_set(dbcsr_Gamma_munu_P, my_group_L_size)
      eps_filter = mp2_env%mp2_gpw%eps_filter

      L_counter = 0
      DO kkb = my_group_L_start, my_group_L_end
         L_counter = L_counter+1
         ! One-index transformed Gamma_3
         ALLOCATE (dbcsr_Gamma_inu_P(L_counter)%matrix)
         CALL dbcsr_init_p(dbcsr_Gamma_inu_P(L_counter)%matrix)
         CALL dbcsr_create(dbcsr_Gamma_inu_P(L_counter)%matrix, template=mo_coeff_o)
         CALL dbcsr_copy(dbcsr_Gamma_inu_P(L_counter)%matrix, mo_coeff_o)
         CALL dbcsr_set(dbcsr_Gamma_inu_P(L_counter)%matrix, 0.0_dp)
         ! Init Gamma_3 in AO basis
         ALLOCATE (dbcsr_Gamma_munu_P(L_counter)%matrix)
         CALL dbcsr_init_p(dbcsr_Gamma_munu_P(L_counter)%matrix)
         CALL dbcsr_create(dbcsr_Gamma_munu_P(L_counter)%matrix, template=mat_munu%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(dbcsr_Gamma_munu_P(L_counter)%matrix, mat_munu%matrix)
         CALL dbcsr_set(dbcsr_Gamma_munu_P(L_counter)%matrix, 0.0_dp)
      ENDDO

      !! Loup over auxiliary basis functions: multiplication
      L_counter = 0
      DO kkb = my_group_L_start, my_group_L_end
         L_counter = L_counter+1
         ! Do dbcsr multiplication: transform the virtual index
         CALL dbcsr_multiply("N", "T", 1.0_dp, mo_coeff_v, dbcsr_Gamma_3(L_counter)%matrix, &
                             0.0_dp, dbcsr_Gamma_inu_P(L_counter)%matrix, filter_eps=eps_filter)

         !Do dbcsr multiplication: transform the occupied index
         CALL dbcsr_multiply("N", "T", 1.0_dp, dbcsr_Gamma_inu_P(L_counter)%matrix, mo_coeff_o, &
                             0.0_dp, dbcsr_Gamma_munu_P(L_counter)%matrix, filter_eps=eps_filter)
         !
         CALL dbcsr_trace(dbcsr_Gamma_munu_P(L_counter)%matrix, trace_corr)
      ENDDO

      ! Gamma_3 not needed anymore
      L_counter = 0
      DO kkb = my_group_L_start, my_group_L_end
         L_counter = L_counter+1
         CALL dbcsr_release(dbcsr_Gamma_3(L_counter)%matrix)
         DEALLOCATE (dbcsr_Gamma_3(L_counter)%matrix)
      ENDDO
      DEALLOCATE (dbcsr_Gamma_3)

      ! Contract DM with exchange integrals
      !CALL integrate_exchange(qs_env, dbcsr_Gamma_munu_P, mat_munu, para_env_sub, ncol_local, eps_filter, e_axk_corr)
      CALL integrate_exchange(qs_env, dbcsr_Gamma_munu_P, mat_munu, para_env_sub, my_group_L_size, eps_filter, e_axk_corr, &
                              my_group_L_start, my_group_L_end)

      !CALL mp_sum(e_axk_corr,para_env_RPA%group)

      ! Print AXK correlation energy to the file
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F25.14,A4)') 'AXK correlation energy for a quadrature point:', &
         e_axk_corr, ' a.u.'

      L_counter = 0
      DO kkb = my_group_L_start, my_group_L_end
         L_counter = L_counter+1
         CALL dbcsr_release(dbcsr_Gamma_inu_P(L_counter)%matrix)
         CALL dbcsr_release(dbcsr_Gamma_munu_P(L_counter)%matrix)
         DEALLOCATE (dbcsr_Gamma_inu_P(L_counter)%matrix)
         DEALLOCATE (dbcsr_Gamma_munu_P(L_counter)%matrix)
      ENDDO
      DEALLOCATE (dbcsr_Gamma_inu_P)
      DEALLOCATE (dbcsr_Gamma_munu_P)

      CALL timestop(handle)

   END SUBROUTINE compute_axk_ener

! **************************************************************************************************
!> \brief ...
!> \param dimen_RI ...
!> \param ngroup ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param num_proc ...
!> \param mepos ...
! **************************************************************************************************
   SUBROUTINE get_aux_sizes(dimen_RI, ngroup, my_group_L_size, my_group_L_start, my_group_L_end, num_proc, mepos)
      INTEGER                                            :: dimen_RI, ngroup, my_group_L_size, &
                                                            my_group_L_start, my_group_L_end, &
                                                            num_proc, mepos

      INTEGER                                            :: color_sub, igroup, iitmp(2)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array, sizes_array, starts_array

      color_sub = mepos/num_proc

      ALLOCATE (sizes_array(0:ngroup-1))
      sizes_array = 0
      ALLOCATE (starts_array(0:ngroup-1))
      starts_array = 0
      ALLOCATE (ends_array(0:ngroup-1))
      ends_array = 0

      DO igroup = 0, ngroup-1
         iitmp = get_limit(dimen_RI, ngroup, igroup)
         starts_array(igroup) = iitmp(1)
         ends_array(igroup) = iitmp(2)
         sizes_array(igroup) = iitmp(2)-iitmp(1)+1
      ENDDO

      my_group_L_size = sizes_array(color_sub)
      my_group_L_start = starts_array(color_sub)
      my_group_L_end = ends_array(color_sub)

      DEALLOCATE (sizes_array)
      DEALLOCATE (starts_array)
      DEALLOCATE (ends_array)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Contract RPA-AXK density matrix with HF exchange integrals and evaluate the correction
!> \param qs_env ...
!> \param dbcsr_Gamma_munu_P ...  AXK density matrix in AO basis to be contracted
!> \param mat_munu ...
!> \param para_env_sub ...
!> \param P_stack_size ...
!> \param eps_filter ...
!> \param axk_corr ...  The AXK energy correction
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \author Vladimir Rybkin, 08/2016
! **************************************************************************************************
   SUBROUTINE integrate_exchange(qs_env, dbcsr_Gamma_munu_P, mat_munu, para_env_sub, P_stack_size, &
                                 eps_filter, axk_corr, &
                                 my_group_L_start, my_group_L_end)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_Gamma_munu_P
      TYPE(dbcsr_p_type)                                 :: mat_munu
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER                                            :: P_stack_size
      REAL(KIND=dp)                                      :: eps_filter, axk_corr
      INTEGER                                            :: my_group_L_start, my_group_L_end

      CHARACTER(LEN=*), PARAMETER :: routineN = 'integrate_exchange', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: aux, handle, irep, kkb, n_rep_hf, ns
      LOGICAL                                            :: my_recalc_hfx_integrals
      REAL(KIND=dp)                                      :: e_axk_P, ex_energy
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_work_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_2d, rho_ao_2d
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: hfx_sections

      CALL timeset(routineN, handle)

      ! Get qs environment
      NULLIFY (energy)
      CALL get_qs_env(qs_env, &
                      energy=energy)

      ! hfx section
      CALL hfx_create_subgroup(qs_env, para_env_sub, hfx_sections, x_data, n_rep_hf)

      ! create a working rho environment
      NULLIFY (rho_work_ao)
      CALL dbcsr_allocate_matrix_set(rho_work_ao, 1)
      ALLOCATE (rho_work_ao(1)%matrix)
      CALL dbcsr_init_p(rho_work_ao(1)%matrix)
      CALL dbcsr_create(rho_work_ao(1)%matrix, template=mat_munu%matrix)

      ! For the first aux function in the group we recalculate integrals, but only for the first
      my_recalc_hfx_integrals = .TRUE.

      NULLIFY (mat_2d)
      CALL dbcsr_allocate_matrix_set(mat_2d, 1, 1)
      ALLOCATE (mat_2d(1, 1)%matrix)
      CALL dbcsr_init_p(mat_2d(1, 1)%matrix)
      CALL dbcsr_create(mat_2d(1, 1)%matrix, template=mat_munu%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(mat_2d(1, 1)%matrix, mat_munu%matrix)

      ! The loop over auxiliary basis functions
      axk_corr = 0.0_dp
      !DO aux = 1, P_stack_size
      P_stack_size = P_stack_size
      aux = 0
      DO kkb = my_group_L_start, my_group_L_end
         aux = aux+1

         CALL dbcsr_copy(rho_work_ao(1)%matrix, dbcsr_Gamma_munu_P(aux)%matrix)
         ! save old EX energy
         ex_energy = energy%ex

         DO irep = 1, n_rep_hf
            ns = SIZE(rho_work_ao)
            rho_ao_2d(1:ns, 1:1) => rho_work_ao(1:ns)

            CALL dbcsr_set(mat_2d(1, 1)%matrix, 0.0_dp)
            CALL integrate_four_center(qs_env, mat_2d, energy, rho_ao_2d, hfx_sections, &
                                       para_env_sub, my_recalc_hfx_integrals, irep, .TRUE., &
                                       ispin=1, axk_x_data=x_data)
         END DO

         ! restore original EX energy
         energy%ex = ex_energy
         my_recalc_hfx_integrals = .FALSE.
         ! One more dbcsr multiplication and trace
         CALL dbcsr_multiply("T", "N", 1.0_dp, mat_2d(1, 1)%matrix, dbcsr_Gamma_munu_P(aux)%matrix, &
                             0.0_dp, dbcsr_Gamma_munu_P(aux)%matrix, filter_eps=eps_filter)
         CALL dbcsr_trace(dbcsr_Gamma_munu_P(aux)%matrix, e_axk_p)
         axk_corr = axk_corr+e_axk_P
      ENDDO

      CALL dbcsr_release(mat_2d(1, 1)%matrix)
      ! release rho stuff
      CALL dbcsr_release(mat_2d(1, 1)%matrix)
      DEALLOCATE (mat_2d(1, 1)%matrix)
      DEALLOCATE (mat_2d)
      CALL dbcsr_release(rho_work_ao(1)%matrix)
      DEALLOCATE (rho_work_ao(1)%matrix)
      DEALLOCATE (rho_work_ao)
      CALL hfx_release(x_data)

      CALL timestop(handle)

   END SUBROUTINE

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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tau_wj_2d

      Rc_max = 1.0E+7

      multiplicator = Rc_max**(1.0_dp/(REAL(nR, KIND=dp)-1.0_dp))

      DO num_integ_points = 1, 20

         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ALLOCATE (tau_tj(num_integ_points))
         tau_tj = 0.0_dp
         ALLOCATE (tau_wj_2d(num_integ_points, num_integ_points))
         tau_wj_2d = 0.0_dp
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

            CALL get_rpa_least_sq_weights_tau(num_integ_points, tau_tj, tau_wj_2d, tj, &
                                              1.0_dp, Rc, max_error, 200)

            IF (iw > 0) THEN
               WRITE (iw, '(T2, I3, F12.1, ES12.3)') num_integ_points, Rc, max_error
            END IF

         END DO

         DEALLOCATE (x_tw, tau_tj, tau_wj_2d, tau_wj, wj, tj)

      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param delta_corr ...
!> \param delta_corr_nm ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param kp_grid ...
!> \param homo ...
!> \param nmo ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param omega ...
!> \param fm_mo_coeff ...
!> \param Eigenval ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param first_cycle_periodic_correction ...
!> \param kpoints ...
!> \param do_mo_coeff_Gamma_only ...
!> \param num_kp_grids ...
!> \param eps_kpoint ...
!> \param do_extra_kpoints ...
!> \param do_aux_bas ...
!> \param frac_aux_mos ...
! **************************************************************************************************
   SUBROUTINE calc_periodic_correction(delta_corr, delta_corr_nm, qs_env, para_env, para_env_RPA, kp_grid, homo, nmo, &
                                       gw_corr_lev_occ, gw_corr_lev_virt, omega, fm_mo_coeff, Eigenval, &
                                       matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                       first_cycle_periodic_correction, kpoints, do_mo_coeff_Gamma_only, &
                                       num_kp_grids, eps_kpoint, do_extra_kpoints, do_aux_bas, frac_aux_mos)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: delta_corr_nm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA
      INTEGER, DIMENSION(:), POINTER                     :: kp_grid
      INTEGER                                            :: homo, nmo, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt
      REAL(KIND=dp)                                      :: omega
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      LOGICAL :: first_cycle_periodic_correction
      TYPE(kpoint_type), POINTER                         :: kpoints
      LOGICAL                                            :: do_mo_coeff_Gamma_only
      INTEGER                                            :: num_kp_grids
      REAL(KIND=dp)                                      :: eps_kpoint
      LOGICAL                                            :: do_extra_kpoints, do_aux_bas
      REAL(KIND=dp)                                      :: frac_aux_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_periodic_correction', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_head, zeta_head
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv

      CALL timeset(routineN, handle)

      IF (first_cycle_periodic_correction) THEN

         CALL get_kpoints(qs_env, kpoints, kp_grid, num_kp_grids, para_env, h_inv, nmo, do_mo_coeff_Gamma_only, &
                          do_extra_kpoints)

         CALL get_berry_phase(qs_env, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, fm_mo_coeff, &
                              para_env, do_mo_coeff_Gamma_only, homo, nmo, gw_corr_lev_virt, eps_kpoint, do_aux_bas, &
                              frac_aux_mos)

      END IF

      CALL compute_eps_head_Berry(eps_head, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, para_env_RPA, &
                                  qs_env, homo, Eigenval, omega)

      CALL compute_zeta_head(zeta_head, eps_head, kpoints)

      CALL kpoint_sum_for_zeta_head_Berry(delta_corr, delta_corr_nm, zeta_head, kpoints, qs_env, &
                                          matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                          homo, nmo, gw_corr_lev_occ, gw_corr_lev_virt, para_env_RPA, &
                                          do_extra_kpoints)

      DEALLOCATE (eps_head, zeta_head)

      first_cycle_periodic_correction = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param eps_head ...
!> \param kpoints ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param para_env_RPA ...
!> \param qs_env ...
!> \param homo ...
!> \param Eigenval ...
!> \param omega ...
! **************************************************************************************************
   SUBROUTINE compute_eps_head_Berry(eps_head, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, para_env_RPA, &
                                     qs_env, homo, Eigenval, omega)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_head
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: homo
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      REAL(KIND=dp)                                      :: omega

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_eps_head_Berry', &
         routineP = moduleN//':'//routineN

      INTEGER :: col, col_end_in_block, col_offset, col_size, handle, i_col, i_row, ikp, nkp, row, &
         row_offset, row_size, row_start_in_block
      REAL(KIND=dp)                                      :: abs_k_square, cell_volume, &
                                                            correct_kpoint(3), cos_square, &
                                                            eigen_diff, relative_kpoint(3), &
                                                            sin_square
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: P_head
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, cell=cell)
      CALL get_cell(cell=cell, deth=cell_volume)

      NULLIFY (data_block)

      nkp = kpoints%nkp

      ALLOCATE (P_head(nkp))
      P_head(:) = 0.0_dp

      ALLOCATE (eps_head(nkp))
      eps_head(:) = 0.0_dp

      DO ikp = 1, nkp

         relative_kpoint(1:3) = MATMUL(cell%hmat, kpoints%xkp(1:3, ikp))

         correct_kpoint(1:3) = twopi*kpoints%xkp(1:3, ikp)

         abs_k_square = (correct_kpoint(1))**2+(correct_kpoint(2))**2+(correct_kpoint(3))**2

         ! real part of the Berry phase
         CALL dbcsr_iterator_start(iter, matrix_berry_re_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            IF (row_offset+row_size <= homo .OR. col_offset > homo) CYCLE

            IF (row_offset <= homo) THEN
               row_start_in_block = homo-row_offset+2
            ELSE
               row_start_in_block = 1
            END IF

            IF (col_offset+col_size-1 > homo) THEN
               col_end_in_block = homo-col_offset+1
            ELSE
               col_end_in_block = col_size
            END IF

            DO i_row = row_start_in_block, row_size

               DO i_col = 1, col_end_in_block

                  eigen_diff = Eigenval(i_col+col_offset-1)-Eigenval(i_row+row_offset-1)

                  cos_square = (data_block(i_row, i_col))**2

                  P_head(ikp) = P_head(ikp)+2.0_dp*eigen_diff/(omega**2+eigen_diff**2)*cos_square/abs_k_square

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

         ! imaginary part of the Berry phase
         CALL dbcsr_iterator_start(iter, matrix_berry_im_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            IF (row_offset+row_size <= homo .OR. col_offset > homo) CYCLE

            IF (row_offset <= homo) THEN
               row_start_in_block = homo-row_offset+2
            ELSE
               row_start_in_block = 1
            END IF

            IF (col_offset+col_size-1 > homo) THEN
               col_end_in_block = homo-col_offset+1
            ELSE
               col_end_in_block = col_size
            END IF

            DO i_row = row_start_in_block, row_size

               DO i_col = 1, col_end_in_block

                  eigen_diff = Eigenval(i_col+col_offset-1)-Eigenval(i_row+row_offset-1)

                  sin_square = (data_block(i_row, i_col))**2

                  P_head(ikp) = P_head(ikp)+2.0_dp*eigen_diff/(omega**2+eigen_diff**2)*sin_square/abs_k_square

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO

      CALL mp_sum(P_head, para_env_RPA%group)

      ! normalize eps_head
      ! 2.0_dp due to closed shell
      eps_head(:) = 1.0_dp-2.0_dp*P_head(:)/cell_volume*fourpi

      DEALLOCATE (P_head)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param fm_mo_coeff ...
!> \param para_env ...
!> \param do_mo_coeff_Gamma_only ...
!> \param homo ...
!> \param nmo ...
!> \param gw_corr_lev_virt ...
!> \param eps_kpoint ...
!> \param do_aux_bas ...
!> \param frac_aux_mos ...
! **************************************************************************************************
   SUBROUTINE get_berry_phase(qs_env, kpoints, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, fm_mo_coeff, para_env, &
                              do_mo_coeff_Gamma_only, homo, nmo, gw_corr_lev_virt, eps_kpoint, do_aux_bas, &
                              frac_aux_mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: do_mo_coeff_Gamma_only
      INTEGER                                            :: homo, nmo, gw_corr_lev_virt
      REAL(KIND=dp)                                      :: eps_kpoint
      LOGICAL                                            :: do_aux_bas
      REAL(KIND=dp)                                      :: frac_aux_mos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_berry_phase', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col_index, handle, i_col_local, iab, &
                                                            ikind, ikp, nao_aux, ncol_local, &
                                                            nkind, nkp, nmo_for_aux_bas
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(dp)                                           :: abs_kpoint, correct_kpoint(3), &
                                                            scale_kpoint
      REAL(KIND=dp), DIMENSION(:), POINTER               :: evals_P, evals_P_sqrt_inv
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_aux_aux
      TYPE(cp_fm_type), POINTER :: fm_mat_eigv_P, fm_mat_P, fm_mat_P_sqrt_inv, &
         fm_mat_s_aux_aux_inv, fm_mat_scaled_eigv_P, fm_mat_work_aux_aux
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_s_aux_aux, &
                                                            matrix_s_aux_orb
      TYPE(dbcsr_type), POINTER :: cosmat, cosmat_desymm, mat_mo_coeff_aux, mat_mo_coeff_aux_2, &
         mat_mo_coeff_Gamma_all, mat_mo_coeff_Gamma_occ_and_GW, mat_mo_coeff_im, mat_mo_coeff_re, &
         mat_work_aux_orb, mat_work_aux_orb_2, matrix_P, matrix_P_sqrt, matrix_P_sqrt_inv, &
         matrix_s_inv_aux_aux, sinmat, sinmat_desymm, tmp
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: gw_aux_basis_set_list, orb_basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_gw_aux
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sab_orb_mic, sgwgw_list, &
                                                            sgworb_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      nkp = kpoints%nkp

      NULLIFY (matrix_berry_re_mo_mo, matrix_s, cell, matrix_berry_im_mo_mo, sinmat, cosmat, tmp, &
               cosmat_desymm, sinmat_desymm, qs_kind_set, orb_basis_set_list, sab_orb_mic)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      matrix_s=matrix_s, &
                      qs_kind_set=qs_kind_set, &
                      nkind=nkind, &
                      ks_env=ks_env, &
                      sab_orb=sab_orb)

      ALLOCATE (orb_basis_set_list(nkind))
      CALL basis_set_list_setup(orb_basis_set_list, "ORB", qs_kind_set)

      CALL setup_neighbor_list(sab_orb_mic, orb_basis_set_list, qs_env=qs_env, mic=.FALSE.)

      ! create dbcsr matrix of mo_coeff for multiplcation
      NULLIFY (mat_mo_coeff_re)
      CALL dbcsr_init_p(mat_mo_coeff_re)
      CALL dbcsr_create(matrix=mat_mo_coeff_re, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_mo_coeff_im)
      CALL dbcsr_init_p(mat_mo_coeff_im)
      CALL dbcsr_create(matrix=mat_mo_coeff_im, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_mo_coeff_Gamma_all)
      CALL dbcsr_init_p(mat_mo_coeff_Gamma_all)
      CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_all, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mo_coeff, mat_mo_coeff_Gamma_all, keep_sparsity=.FALSE.)

      NULLIFY (mat_mo_coeff_Gamma_occ_and_GW)
      CALL dbcsr_init_p(mat_mo_coeff_Gamma_occ_and_GW)
      CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_occ_and_GW, &
                        template=matrix_s(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL copy_fm_to_dbcsr(fm_mo_coeff, mat_mo_coeff_Gamma_occ_and_GW, keep_sparsity=.FALSE.)

      IF (.NOT. do_aux_bas) THEN

         ! allocate intermediate matrices
         CALL dbcsr_init_p(cosmat)
         CALL dbcsr_init_p(sinmat)
         CALL dbcsr_init_p(tmp)
         CALL dbcsr_init_p(cosmat_desymm)
         CALL dbcsr_init_p(sinmat_desymm)
         CALL dbcsr_create(matrix=cosmat, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=sinmat, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=tmp, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=cosmat_desymm, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=sinmat_desymm, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(cosmat, matrix_s(1)%matrix)
         CALL dbcsr_copy(sinmat, matrix_s(1)%matrix)
         CALL dbcsr_set(cosmat, 0.0_dp)
         CALL dbcsr_set(sinmat, 0.0_dp)

         CALL dbcsr_allocate_matrix_set(matrix_berry_re_mo_mo, nkp)
         CALL dbcsr_allocate_matrix_set(matrix_berry_im_mo_mo, nkp)

      END IF

      IF (do_aux_bas) THEN

         NULLIFY (gw_aux_basis_set_list)
         ALLOCATE (gw_aux_basis_set_list(nkind))

         DO ikind = 1, nkind

            NULLIFY (gw_aux_basis_set_list(ikind)%gto_basis_set)

            NULLIFY (basis_set_gw_aux)

            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_gw_aux, basis_type="AUX_GW")
            CPASSERT(ASSOCIATED(basis_set_gw_aux))

            basis_set_gw_aux%kind_radius = orb_basis_set_list(ikind)%gto_basis_set%kind_radius

            gw_aux_basis_set_list(ikind)%gto_basis_set => basis_set_gw_aux

         END DO

         ! neighbor lists
         NULLIFY (sgwgw_list, sgworb_list)
         CALL setup_neighbor_list(sgwgw_list, gw_aux_basis_set_list, qs_env=qs_env)
         CALL setup_neighbor_list(sgworb_list, gw_aux_basis_set_list, orb_basis_set_list, qs_env=qs_env)

         NULLIFY (matrix_s_aux_aux, matrix_s_aux_orb)

         ! build overlap matrix in gw aux basis and the mixed gw aux basis-orb basis
         CALL build_overlap_matrix_simple(ks_env, matrix_s_aux_aux, &
                                          gw_aux_basis_set_list, gw_aux_basis_set_list, sgwgw_list)

         CALL build_overlap_matrix_simple(ks_env, matrix_s_aux_orb, &
                                          gw_aux_basis_set_list, orb_basis_set_list, sgworb_list)

         CALL dbcsr_get_info(matrix_s_aux_aux(1)%matrix, nfullrows_total=nao_aux)

         nmo_for_aux_bas = FLOOR(frac_aux_mos*REAL(nao_aux, KIND=dp))

         CALL cp_fm_struct_create(fm_struct_aux_aux, &
                                  context=fm_mo_coeff%matrix_struct%context, &
                                  nrow_global=nao_aux, &
                                  ncol_global=nao_aux, &
                                  para_env=para_env)

         NULLIFY (mat_work_aux_orb)
         CALL dbcsr_init_p(mat_work_aux_orb)
         CALL dbcsr_create(matrix=mat_work_aux_orb, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_work_aux_orb_2)
         CALL dbcsr_init_p(mat_work_aux_orb_2)
         CALL dbcsr_create(matrix=mat_work_aux_orb_2, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_mo_coeff_aux)
         CALL dbcsr_init_p(mat_mo_coeff_aux)
         CALL dbcsr_create(matrix=mat_mo_coeff_aux, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_mo_coeff_aux_2)
         CALL dbcsr_init_p(mat_mo_coeff_aux_2)
         CALL dbcsr_create(matrix=mat_mo_coeff_aux_2, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_s_inv_aux_aux)
         CALL dbcsr_init_p(matrix_s_inv_aux_aux)
         CALL dbcsr_create(matrix=matrix_s_inv_aux_aux, &
                           template=matrix_s_aux_aux(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_P)
         CALL dbcsr_init_p(matrix_P)
         CALL dbcsr_create(matrix=matrix_P, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_P_sqrt)
         CALL dbcsr_init_p(matrix_P_sqrt)
         CALL dbcsr_create(matrix=matrix_P_sqrt, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (matrix_P_sqrt_inv)
         CALL dbcsr_init_p(matrix_P_sqrt_inv)
         CALL dbcsr_create(matrix=matrix_P_sqrt_inv, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (fm_mat_s_aux_aux_inv)
         CALL cp_fm_create(fm_mat_s_aux_aux_inv, fm_struct_aux_aux, name="inverse overlap mat")

         NULLIFY (fm_mat_work_aux_aux)
         CALL cp_fm_create(fm_mat_work_aux_aux, fm_struct_aux_aux, name="work mat")

         NULLIFY (fm_mat_P)
         CALL cp_fm_create(fm_mat_P, fm_mo_coeff%matrix_struct)

         NULLIFY (fm_mat_eigv_P)
         CALL cp_fm_create(fm_mat_eigv_P, fm_mo_coeff%matrix_struct)

         NULLIFY (fm_mat_scaled_eigv_P)
         CALL cp_fm_create(fm_mat_scaled_eigv_P, fm_mo_coeff%matrix_struct)

         NULLIFY (fm_mat_P_sqrt_inv)
         CALL cp_fm_create(fm_mat_P_sqrt_inv, fm_mo_coeff%matrix_struct)

         NULLIFY (evals_P)
         ALLOCATE (evals_P(nmo))

         NULLIFY (evals_P_sqrt_inv)
         ALLOCATE (evals_P_sqrt_inv(nmo))

         CALL copy_dbcsr_to_fm(matrix_s_aux_aux(1)%matrix, fm_mat_s_aux_aux_inv)
         ! Calculate S_inverse
         CALL cp_fm_cholesky_decompose(fm_mat_s_aux_aux_inv)
         CALL cp_fm_cholesky_invert(fm_mat_s_aux_aux_inv)
         ! Symmetrize the guy
         CALL cp_fm_upper_to_full(fm_mat_s_aux_aux_inv, fm_mat_work_aux_aux)

         CALL copy_fm_to_dbcsr(fm_mat_s_aux_aux_inv, matrix_s_inv_aux_aux, keep_sparsity=.FALSE.)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_inv_aux_aux, matrix_s_aux_orb(1)%matrix, 0.0_dp, mat_work_aux_orb, &
                             filter_eps=1.0E-15_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, mat_work_aux_orb, mat_mo_coeff_Gamma_all, 0.0_dp, mat_mo_coeff_aux_2, &
                             last_column=nmo_for_aux_bas, filter_eps=1.0E-15_dp)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_s_aux_aux(1)%matrix, mat_mo_coeff_aux_2, 0.0_dp, mat_work_aux_orb, &
                             filter_eps=1.0E-15_dp)

         CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_aux_2, mat_work_aux_orb, 0.0_dp, matrix_P, &
                             filter_eps=1.0E-15_dp)

         CALL copy_dbcsr_to_fm(matrix_P, fm_mat_P)

         CALL cp_fm_syevd(fm_mat_P, fm_mat_eigv_P, evals_P)

         ! only invert the eigenvalues which correspond to the MOs used in the aux. basis
         evals_P_sqrt_inv(1:nmo-nmo_for_aux_bas) = 0.0_dp
         evals_P_sqrt_inv(nmo-nmo_for_aux_bas+1:nmo) = 1.0_dp/SQRT(evals_P(nmo-nmo_for_aux_bas+1:nmo))

         CALL cp_fm_to_fm(fm_mat_eigv_P, fm_mat_scaled_eigv_P)

         CALL cp_fm_get_info(matrix=fm_mat_scaled_eigv_P, &
                             ncol_local=ncol_local, &
                             col_indices=col_indices)

         ! multiply eigenvectors with inverse sqrt of eigenvalues
         DO i_col_local = 1, ncol_local

            col_index = col_indices(i_col_local)

            fm_mat_scaled_eigv_P%local_data(:, i_col_local) = &
               fm_mat_scaled_eigv_P%local_data(:, i_col_local)*evals_P_sqrt_inv(col_index)

         END DO

         CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mat_eigv_P, matrix_b=fm_mat_scaled_eigv_P, beta=0.0_dp, &
                      matrix_c=fm_mat_P_sqrt_inv)

         CALL copy_fm_to_dbcsr(fm_mat_P_sqrt_inv, matrix_P_sqrt_inv, keep_sparsity=.FALSE.)

         CALL dbcsr_multiply('N', 'N', 1.0_dp, mat_mo_coeff_aux_2, matrix_P_sqrt_inv, 0.0_dp, mat_mo_coeff_aux, &
                             filter_eps=1.0E-15_dp)

         ! allocate intermediate matrices
         CALL dbcsr_init_p(cosmat)
         CALL dbcsr_init_p(sinmat)
         CALL dbcsr_init_p(tmp)
         CALL dbcsr_init_p(cosmat_desymm)
         CALL dbcsr_init_p(sinmat_desymm)
         CALL dbcsr_create(matrix=cosmat, template=matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_create(matrix=sinmat, template=matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_create(matrix=tmp, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=cosmat_desymm, &
                           template=matrix_s_aux_aux(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix=sinmat_desymm, &
                           template=matrix_s_aux_aux(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(cosmat, matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_copy(sinmat, matrix_s_aux_aux(1)%matrix)
         CALL dbcsr_set(cosmat, 0.0_dp)
         CALL dbcsr_set(sinmat, 0.0_dp)

         CALL dbcsr_allocate_matrix_set(matrix_berry_re_mo_mo, nkp)
         CALL dbcsr_allocate_matrix_set(matrix_berry_im_mo_mo, nkp)

         ! allocate the new MO coefficients in the aux basis
         CALL dbcsr_release_p(mat_mo_coeff_Gamma_all)
         CALL dbcsr_release_p(mat_mo_coeff_Gamma_occ_and_GW)

         NULLIFY (mat_mo_coeff_Gamma_all)
         CALL dbcsr_init_p(mat_mo_coeff_Gamma_all)
         CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_all, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_copy(mat_mo_coeff_Gamma_all, mat_mo_coeff_aux)

         NULLIFY (mat_mo_coeff_Gamma_occ_and_GW)
         CALL dbcsr_init_p(mat_mo_coeff_Gamma_occ_and_GW)
         CALL dbcsr_create(matrix=mat_mo_coeff_Gamma_occ_and_GW, &
                           template=matrix_s_aux_orb(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_copy(mat_mo_coeff_Gamma_occ_and_GW, mat_mo_coeff_aux)

         DEALLOCATE (evals_P, evals_P_sqrt_inv)

      END IF

      CALL remove_unnecessary_blocks(mat_mo_coeff_Gamma_occ_and_GW, homo, gw_corr_lev_virt)

      DO ikp = 1, nkp

         ALLOCATE (matrix_berry_re_mo_mo(ikp)%matrix)
         CALL dbcsr_init_p(matrix_berry_re_mo_mo(ikp)%matrix)
         CALL dbcsr_create(matrix_berry_re_mo_mo(ikp)%matrix, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, matrix_berry_re_mo_mo(ikp)%matrix)
         CALL dbcsr_set(matrix_berry_re_mo_mo(ikp)%matrix, 0.0_dp)

         ALLOCATE (matrix_berry_im_mo_mo(ikp)%matrix)
         CALL dbcsr_init_p(matrix_berry_im_mo_mo(ikp)%matrix)
         CALL dbcsr_create(matrix_berry_im_mo_mo(ikp)%matrix, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, matrix_berry_im_mo_mo(ikp)%matrix)
         CALL dbcsr_set(matrix_berry_im_mo_mo(ikp)%matrix, 0.0_dp)

         correct_kpoint(1:3) = -twopi*kpoints%xkp(1:3, ikp)

         abs_kpoint = SQRT(correct_kpoint(1)**2+correct_kpoint(2)**2+correct_kpoint(3)**2)

         IF (abs_kpoint < eps_kpoint) THEN

            scale_kpoint = eps_kpoint/abs_kpoint
            correct_kpoint(:) = correct_kpoint(:)*scale_kpoint

         END IF

         ! get the Berry phase
         IF (do_aux_bas) THEN
            CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, correct_kpoint, sab_orb_external=sab_orb_mic, &
                                           basis_type="AUX_GW")
         ELSE
            CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, correct_kpoint, sab_orb_external=sab_orb_mic, &
                                           basis_type="ORB")
         END IF

         IF (do_mo_coeff_Gamma_only) THEN

            CALL dbcsr_desymmetrize(cosmat, cosmat_desymm)

            CALL dbcsr_multiply('N', 'N', 1.0_dp, cosmat_desymm, mat_mo_coeff_Gamma_occ_and_GW, 0.0_dp, tmp, &
                                filter_eps=1.0E-15_dp)

            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_re_mo_mo(ikp)%matrix, filter_eps=1.0E-15_dp)

            CALL dbcsr_desymmetrize(sinmat, sinmat_desymm)

            CALL dbcsr_multiply('N', 'N', 1.0_dp, sinmat_desymm, mat_mo_coeff_Gamma_occ_and_GW, 0.0_dp, tmp, &
                                filter_eps=1.0E-15_dp)

            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_im_mo_mo(ikp)%matrix, filter_eps=1.0E-15_dp)

         ELSE

            ! get mo coeff at the ikp
            CALL copy_fm_to_dbcsr(kpoints%kp_env(ikp)%kpoint_env%mos(1, 1)%mo_set%mo_coeff, &
                                  mat_mo_coeff_re, keep_sparsity=.FALSE.)

            CALL copy_fm_to_dbcsr(kpoints%kp_env(ikp)%kpoint_env%mos(2, 1)%mo_set%mo_coeff, &
                                  mat_mo_coeff_im, keep_sparsity=.FALSE.)

            CALL dbcsr_desymmetrize(cosmat, cosmat_desymm)

            CALL dbcsr_desymmetrize(sinmat, sinmat_desymm)

            ! I.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, cosmat_desymm, mat_mo_coeff_re, 0.0_dp, tmp)

            ! I.1
            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_re_mo_mo(ikp)%matrix)

            ! II.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, sinmat_desymm, mat_mo_coeff_re, 0.0_dp, tmp)

            ! II.5
            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 0.0_dp, &
                                matrix_berry_im_mo_mo(ikp)%matrix)

            ! III.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, cosmat_desymm, mat_mo_coeff_im, 0.0_dp, tmp)

            ! III.7
            CALL dbcsr_multiply('T', 'N', 1.0_dp, mat_mo_coeff_Gamma_all, tmp, 1.0_dp, &
                                matrix_berry_im_mo_mo(ikp)%matrix)

            ! IV.
            CALL dbcsr_multiply('N', 'N', 1.0_dp, sinmat_desymm, mat_mo_coeff_im, 0.0_dp, tmp)

            ! IV.3
            CALL dbcsr_multiply('T', 'N', -1.0_dp, mat_mo_coeff_Gamma_all, tmp, 1.0_dp, &
                                matrix_berry_re_mo_mo(ikp)%matrix)

         END IF

         IF (abs_kpoint < eps_kpoint) THEN

            CALL dbcsr_scale(matrix_berry_im_mo_mo(ikp)%matrix, 1.0_dp/scale_kpoint)
            CALL dbcsr_set(matrix_berry_re_mo_mo(ikp)%matrix, 0.0_dp)
            CALL dbcsr_add_on_diag(matrix_berry_re_mo_mo(ikp)%matrix, 1.0_dp)

         END IF

      END DO

      CALL dbcsr_release_p(cosmat)
      CALL dbcsr_release_p(sinmat)
      CALL dbcsr_release_p(mat_mo_coeff_re)
      CALL dbcsr_release_p(mat_mo_coeff_im)
      CALL dbcsr_release_p(mat_mo_coeff_Gamma_all)
      CALL dbcsr_release_p(mat_mo_coeff_Gamma_occ_and_GW)
      CALL dbcsr_release_p(tmp)
      CALL dbcsr_release_p(cosmat_desymm)
      CALL dbcsr_release_p(sinmat_desymm)
      DEALLOCATE (orb_basis_set_list)

      IF (ASSOCIATED(sab_orb_mic)) THEN
         DO iab = 1, SIZE(sab_orb_mic)
            CALL deallocate_neighbor_list_set(sab_orb_mic(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (sab_orb_mic)
      END IF

      IF (do_aux_bas) THEN

         DEALLOCATE (gw_aux_basis_set_list)
         CALL dbcsr_deallocate_matrix_set(matrix_s_aux_aux)
         CALL dbcsr_deallocate_matrix_set(matrix_s_aux_orb)
         CALL dbcsr_release_p(mat_work_aux_orb)
         CALL dbcsr_release_p(mat_work_aux_orb_2)
         CALL dbcsr_release_p(mat_mo_coeff_aux)
         CALL dbcsr_release_p(mat_mo_coeff_aux_2)
         CALL dbcsr_release_p(matrix_s_inv_aux_aux)
         CALL dbcsr_release_p(matrix_P)
         CALL dbcsr_release_p(matrix_P_sqrt)
         CALL dbcsr_release_p(matrix_P_sqrt_inv)

         CALL cp_fm_struct_release(fm_struct_aux_aux)

         CALL cp_fm_release(fm_mat_s_aux_aux_inv)
         CALL cp_fm_release(fm_mat_work_aux_aux)
         CALL cp_fm_release(fm_mat_P)
         CALL cp_fm_release(fm_mat_eigv_P)
         CALL cp_fm_release(fm_mat_scaled_eigv_P)
         CALL cp_fm_release(fm_mat_P_sqrt_inv)

         ! Deallocate the neighbor list structure
         IF (ASSOCIATED(sgwgw_list)) THEN
            DO iab = 1, SIZE(sgwgw_list)
               CALL deallocate_neighbor_list_set(sgwgw_list(iab)%neighbor_list_set)
            END DO
            DEALLOCATE (sgwgw_list)
         END IF

         IF (ASSOCIATED(sgworb_list)) THEN
            DO iab = 1, SIZE(sgworb_list)
               CALL deallocate_neighbor_list_set(sgworb_list(iab)%neighbor_list_set)
            END DO
            DEALLOCATE (sgworb_list)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_mo_coeff_Gamma_occ_and_GW ...
!> \param homo ...
!> \param gw_corr_lev_virt ...
! **************************************************************************************************
   SUBROUTINE remove_unnecessary_blocks(mat_mo_coeff_Gamma_occ_and_GW, homo, gw_corr_lev_virt)

      TYPE(dbcsr_type), POINTER                          :: mat_mo_coeff_Gamma_occ_and_GW
      INTEGER                                            :: homo, gw_corr_lev_virt

      INTEGER                                            :: col, col_offset, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, mat_mo_coeff_Gamma_occ_and_GW)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        col_offset=col_offset)

         IF (col_offset > homo+gw_corr_lev_virt) THEN

            data_block = 0.0_dp

         END IF

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_filter(mat_mo_coeff_Gamma_occ_and_GW, 1.0E-15_dp)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Berry_nn_local ...
!> \param matrix_berry_mo_mo ...
!> \param para_env ...
!> \param homo ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
! **************************************************************************************************
   SUBROUTINE replicate_Berry_nn(Berry_nn_local, matrix_berry_mo_mo, para_env, homo, gw_corr_lev_occ, gw_corr_lev_virt)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Berry_nn_local
      TYPE(dbcsr_type), POINTER                          :: matrix_berry_mo_mo
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: homo, gw_corr_lev_occ, gw_corr_lev_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_Berry_nn', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, handle, i_row, n_level_gw, row, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      Berry_nn_local = 0.0_dp

      ! fill buffer_send
      CALL dbcsr_iterator_start(iter, matrix_berry_mo_mo)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, row_offset=row_offset)

         IF (row .NE. col) CYCLE

         DO i_row = 1, row_size

            n_level_gw = i_row+row_offset-1

            IF (n_level_gw < 1+homo-gw_corr_lev_occ .OR. n_level_gw > homo+gw_corr_lev_virt) CYCLE

            Berry_nn_local(n_level_gw) = data_block(i_row, i_row)

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL mp_sum(Berry_nn_local, para_env%group)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Berry_ia_local ...
!> \param matrix_berry_mo_mo ...
!> \param para_env ...
!> \param homo ...
! **************************************************************************************************
   SUBROUTINE replicate_Berry_ia(Berry_ia_local, matrix_berry_mo_mo, para_env, homo)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Berry_ia_local
      TYPE(dbcsr_type), POINTER                          :: matrix_berry_mo_mo
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_Berry_ia', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_counter, block_size, col, col_offset, col_size, col_size_modified, &
         col_start_loc_array, entry_counter, handle, iblk, imepos, msg_offset, num_blocks_send, &
         num_entries_send, row, row_offset, row_offset_modified, row_size, row_size_modified, &
         row_start_block, row_start_loc_array
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_rec, num_entries_rec, &
                                                            sizes_rec, sizes_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      num_blocks_send = 0
      num_entries_send = 0

      NULLIFY (data_block)

      CALL dbcsr_iterator_start(iter, matrix_berry_mo_mo)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         IF (row_offset+row_size <= homo .OR. col_offset > homo) CYCLE

         IF (row_offset <= homo) THEN
            row_size_modified = row_size+row_offset-homo-1
         ELSE
            row_size_modified = row_size
         END IF

         IF (col_offset+col_size-1 > homo) THEN
            col_size_modified = homo-col_offset+1
         ELSE
            col_size_modified = col_size
         END IF

         num_blocks_send = num_blocks_send+1
         num_entries_send = num_entries_send+row_size_modified*col_size_modified

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (buffer_rec(0:para_env%num_pe-1))
      ALLOCATE (buffer_send(0:para_env%num_pe-1))

      ALLOCATE (num_entries_rec(0:para_env%num_pe-1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe-1))

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe-1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe-1))

         DO imepos = 0, para_env%num_pe-1
            sizes_send(2*imepos) = num_entries_send
            sizes_send(2*imepos+1) = num_blocks_send
         END DO

         CALL mp_alltoall(sizes_send, sizes_rec, 2, para_env%group)

         DO imepos = 0, para_env%num_pe-1
            num_entries_rec(imepos) = sizes_rec(2*imepos)
            num_blocks_rec(imepos) = sizes_rec(2*imepos+1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send
         num_blocks_rec(0) = num_blocks_send

      END IF

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 5))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send, 5))
         buffer_send(imepos)%indx = 0

      END DO

      entry_counter = 0
      block_counter = 0

      NULLIFY (data_block)

      ! fill buffer_send
      CALL dbcsr_iterator_start(iter, matrix_berry_mo_mo)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         IF (row_offset+row_size <= homo .OR. col_offset > homo) CYCLE

         IF (row_offset <= homo) THEN
            row_size_modified = row_size+row_offset-homo-1
            row_start_block = row_size-row_size_modified+1
            row_offset_modified = homo+1
         ELSE
            row_start_block = 1
            row_size_modified = row_size
            row_offset_modified = row_offset
         END IF

         IF (col_offset+col_size-1 > homo) THEN
            col_size_modified = homo-col_offset+1
         ELSE
            col_size_modified = col_size
         END IF

         block_size = row_size_modified*col_size_modified

         DO imepos = 0, para_env%num_pe-1

            buffer_send(imepos)%msg(entry_counter+1:entry_counter+block_size) = &
               RESHAPE(data_block(row_start_block:row_size, 1:col_size_modified), (/block_size/))

            buffer_send(imepos)%indx(block_counter+1, 1) = row_offset_modified
            buffer_send(imepos)%indx(block_counter+1, 2) = row_size_modified
            buffer_send(imepos)%indx(block_counter+1, 3) = col_offset
            buffer_send(imepos)%indx(block_counter+1, 4) = col_size_modified
            buffer_send(imepos)%indx(block_counter+1, 5) = entry_counter

         END DO

         entry_counter = entry_counter+block_size

         block_counter = block_counter+1

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (sizes_rec(0:para_env%num_pe-1))
      ALLOCATE (sizes_send(0:para_env%num_pe-1))

      DO imepos = 0, para_env%num_pe-1

         sizes_send(imepos) = num_entries_send
         sizes_rec(imepos) = num_entries_rec(imepos)

      END DO

      CALL communicate_buffer(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array, sizes_rec, sizes_send)

      Berry_ia_local = 0.0_dp

      ! fill Berry_ia_local
      DO imepos = 0, para_env%num_pe-1

         DO iblk = 1, num_blocks_rec(imepos)

            row_start_loc_array = buffer_rec(imepos)%indx(iblk, 1)-homo
            row_size = buffer_rec(imepos)%indx(iblk, 2)
            col_start_loc_array = buffer_rec(imepos)%indx(iblk, 3)
            col_size = buffer_rec(imepos)%indx(iblk, 4)
            msg_offset = buffer_rec(imepos)%indx(iblk, 5)

            Berry_ia_local(row_start_loc_array:row_start_loc_array+row_size-1, &
                           col_start_loc_array:col_start_loc_array+col_size-1) = &
               RESHAPE(buffer_rec(imepos)%msg(msg_offset+1:msg_offset+row_size*col_size), &
                       (/row_size, col_size/))

         END DO

      END DO

      ! deallocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%indx)

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
! ...
! **************************************************************************************************
!> \brief ...
!> \param delta_corr ...
!> \param delta_corr_nm ...
!> \param zeta_head ...
!> \param kpoints ...
!> \param qs_env ...
!> \param matrix_berry_re_mo_mo ...
!> \param matrix_berry_im_mo_mo ...
!> \param homo ...
!> \param nmo ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param para_env_RPA ...
!> \param do_extra_kpoints ...
! **************************************************************************************************
   SUBROUTINE kpoint_sum_for_zeta_head_Berry(delta_corr, delta_corr_nm, zeta_head, kpoints, qs_env, matrix_berry_re_mo_mo, &
                                             matrix_berry_im_mo_mo, homo, nmo, gw_corr_lev_occ, gw_corr_lev_virt, &
                                             para_env_RPA, do_extra_kpoints)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta_corr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: delta_corr_nm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zeta_head
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_re_mo_mo, &
                                                            matrix_berry_im_mo_mo
      INTEGER                                            :: homo, nmo, gw_corr_lev_occ, &
                                                            gw_corr_lev_virt
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env_RPA
      LOGICAL                                            :: do_extra_kpoints

      INTEGER                                            :: col, col_offset, col_size, i_col, i_row, &
                                                            ikp, m_level, n_level_gw, nkp, row, &
                                                            row_offset, row_size
      REAL(KIND=dp)                                      :: abs_k_square, cell_volume, &
                                                            check_int_one_over_ksq, contribution, &
                                                            weight
      REAL(KIND=dp), DIMENSION(3)                        :: correct_kpoint
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block, delta_corr_nm_extra
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter, iter_new

      CALL get_qs_env(qs_env=qs_env, cell=cell)

      CALL get_cell(cell=cell, deth=cell_volume)

      nkp = kpoints%nkp

      delta_corr = 0.0_dp
      delta_corr_nm = 0.0_dp

      IF (do_extra_kpoints) THEN
         NULLIFY (delta_corr_nm_extra)
         ALLOCATE (delta_corr_nm_extra(1+homo-gw_corr_lev_occ:homo+gw_corr_lev_virt, 1:nmo))
         delta_corr_nm_extra = 0.0_dp
      END IF

      check_int_one_over_ksq = 0.0_dp

      DO ikp = 1, nkp

         weight = kpoints%wkp(ikp)

         correct_kpoint(1:3) = twopi*kpoints%xkp(1:3, ikp)

         abs_k_square = (correct_kpoint(1))**2+(correct_kpoint(2))**2+(correct_kpoint(3))**2

         ! cos part of the Berry phase
         CALL dbcsr_iterator_start(iter, matrix_berry_re_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO i_col = 1, col_size

               DO n_level_gw = 1+homo-gw_corr_lev_occ, homo+gw_corr_lev_virt

                  IF (n_level_gw == i_col+col_offset-1) THEN

                     DO i_row = 1, row_size

                        contribution = weight*(zeta_head(ikp)-1.0_dp)/abs_k_square*(data_block(i_row, i_col))**2

                        m_level = i_row+row_offset-1

                        IF (.NOT. do_extra_kpoints) THEN

                           delta_corr_nm(n_level_gw, m_level) = delta_corr_nm(n_level_gw, m_level)+contribution

                        ELSE

                           IF (ikp <= nkp*8/9) THEN

                              delta_corr_nm(n_level_gw, m_level) = delta_corr_nm(n_level_gw, m_level)+contribution

                           ELSE

                              delta_corr_nm_extra(n_level_gw, m_level) = delta_corr_nm_extra(n_level_gw, m_level)+contribution

                           END IF

                        END IF

                     END DO

                  END IF

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

         ! the same for the im. part of the Berry phase
         CALL dbcsr_iterator_start(iter_new, matrix_berry_im_mo_mo(ikp)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter_new))

            CALL dbcsr_iterator_next_block(iter_new, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO i_col = 1, col_size

               DO n_level_gw = 1+homo-gw_corr_lev_occ, homo+gw_corr_lev_virt

                  IF (n_level_gw == i_col+col_offset-1) THEN

                     DO i_row = 1, row_size

                        m_level = i_row+row_offset-1

                        contribution = weight*(zeta_head(ikp)-1.0_dp)/abs_k_square*(data_block(i_row, i_col))**2

                        IF (.NOT. do_extra_kpoints) THEN

                           delta_corr_nm(n_level_gw, m_level) = delta_corr_nm(n_level_gw, m_level)+contribution

                        ELSE

                           IF (ikp <= nkp*8/9) THEN

                              delta_corr_nm(n_level_gw, m_level) = delta_corr_nm(n_level_gw, m_level)+contribution

                           ELSE

                              delta_corr_nm_extra(n_level_gw, m_level) = delta_corr_nm_extra(n_level_gw, m_level)+contribution

                           END IF

                        END IF

                     END DO

                  END IF

               END DO

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter_new)

         check_int_one_over_ksq = check_int_one_over_ksq+weight/abs_k_square

      END DO

      ! normalize by the cell volume
      delta_corr_nm = delta_corr_nm/cell_volume*fourpi

      check_int_one_over_ksq = check_int_one_over_ksq/cell_volume

      CALL mp_sum(delta_corr_nm, para_env_RPA%group)

      IF (do_extra_kpoints) THEN

         delta_corr_nm_extra = delta_corr_nm_extra/cell_volume*fourpi

         CALL mp_sum(delta_corr_nm_extra, para_env_RPA%group)

         delta_corr_nm(:, :) = delta_corr_nm(:, :)+(delta_corr_nm(:, :)-delta_corr_nm_extra(:, :))

         DEALLOCATE (delta_corr_nm_extra)

      END IF

      DO n_level_gw = 1+homo-gw_corr_lev_occ, homo+gw_corr_lev_virt

         delta_corr_nm(n_level_gw, 1:n_level_gw-1) = 0.0_dp
         delta_corr_nm(n_level_gw, n_level_gw+1:nmo) = 0.0_dp

      END DO

   END SUBROUTINE

! **************************************************************************************************
! ...
! **************************************************************************************************
!> \brief ...
!> \param zeta_head ...
!> \param eps_head ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE compute_zeta_head(zeta_head, eps_head, kpoints)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zeta_head, eps_head
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_zeta_head', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp, nkp

      CALL timeset(routineN, handle)

      nkp = kpoints%nkp

      ALLOCATE (zeta_head(nkp))

      DO ikp = 1, nkp

         zeta_head(ikp) = 1.0_dp/eps_head(ikp)

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
! ...
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param kp_grid ...
!> \param num_kp_grids ...
!> \param para_env ...
!> \param h_inv ...
!> \param nmo ...
!> \param do_mo_coeff_Gamma_only ...
!> \param do_extra_kpoints ...
! **************************************************************************************************
   SUBROUTINE get_kpoints(qs_env, kpoints, kp_grid, num_kp_grids, para_env, h_inv, nmo, &
                          do_mo_coeff_Gamma_only, do_extra_kpoints)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, DIMENSION(:), POINTER                     :: kp_grid
      INTEGER                                            :: num_kp_grids
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv
      INTEGER                                            :: nmo
      LOGICAL                                            :: do_mo_coeff_Gamma_only, do_extra_kpoints

      INTEGER                                            :: end_kp, i, i_grid_level, ix, iy, iz, &
                                                            nkp_inner_grid, nkp_outer_grid, &
                                                            npoints, start_kp
      INTEGER, DIMENSION(3)                              :: outer_kp_grid
      REAL(KIND=dp)                                      :: kpoint_weight_left, single_weight
      REAL(KIND=dp), DIMENSION(3)                        :: kpt_latt, reducing_factor
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_environment_type), POINTER                 :: qs_env_kp_Gamma_only

      NULLIFY (kpoints, cell, particle_set, qs_env_kp_Gamma_only)

      ! check whether kp_grid includes the Gamma point. If so, abort.
      CPASSERT(MOD(kp_grid(1)*kp_grid(2)*kp_grid(3), 2) == 0)
      IF (do_extra_kpoints) THEN
         CPASSERT(do_mo_coeff_Gamma_only)
      END IF

      IF (do_mo_coeff_Gamma_only) THEN

         outer_kp_grid(1) = kp_grid(1)-1
         outer_kp_grid(2) = kp_grid(2)-1
         outer_kp_grid(3) = kp_grid(3)-1

         CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

         CALL get_cell(cell, h_inv=h_inv)

         CALL kpoint_create(kpoints)

         kpoints%kp_scheme = "GENERAL"
         kpoints%symmetry = .FALSE.
         kpoints%verbose = .FALSE.
         kpoints%full_grid = .FALSE.
         kpoints%use_real_wfn = .FALSE.
         kpoints%eps_geo = 1.e-6_dp
         npoints = kp_grid(1)*kp_grid(2)*kp_grid(3)/2+ &
                   (num_kp_grids-1)*((outer_kp_grid(1)+1)/2*outer_kp_grid(2)*outer_kp_grid(3)-1)

         IF (do_extra_kpoints) THEN

            CPASSERT(num_kp_grids == 1)
            CPASSERT(MOD(kp_grid(1), 4) == 0)
            CPASSERT(MOD(kp_grid(2), 4) == 0)
            CPASSERT(MOD(kp_grid(3), 4) == 0)

         END IF

         IF (do_extra_kpoints) THEN

            npoints = kp_grid(1)*kp_grid(2)*kp_grid(3)/2+kp_grid(1)*kp_grid(2)*kp_grid(3)/2/8

         END IF

         kpoints%full_grid = .TRUE.
         kpoints%nkp = npoints
         ALLOCATE (kpoints%xkp(3, npoints), kpoints%wkp(npoints))
         kpoints%xkp = 0.0_dp
         kpoints%wkp = 0.0_dp

         nkp_outer_grid = outer_kp_grid(1)*outer_kp_grid(2)*outer_kp_grid(3)
         nkp_inner_grid = kp_grid(1)*kp_grid(2)*kp_grid(3)

         i = 0
         reducing_factor(:) = 1.0_dp
         kpoint_weight_left = 1.0_dp

         ! the outer grids
         DO i_grid_level = 1, num_kp_grids-1

            single_weight = kpoint_weight_left/REAL(nkp_outer_grid)

            start_kp = i+1

            DO ix = 1, outer_kp_grid(1)
               DO iy = 1, outer_kp_grid(2)
                  DO iz = 1, outer_kp_grid(3)

                     ! exclude Gamma
                     IF (2*ix-outer_kp_grid(1)-1 == 0 .AND. 2*iy-outer_kp_grid(2)-1 == 0 .AND. &
                         2*iz-outer_kp_grid(3)-1 == 0) CYCLE

                     ! use time reversal symmetry k<->-k
                     IF (2*ix-outer_kp_grid(1)-1 < 0) CYCLE

                     i = i+1
                     kpt_latt(1) = REAL(2*ix-outer_kp_grid(1)-1, KIND=dp)/(2._dp*REAL(outer_kp_grid(1), KIND=dp)) &
                                   *reducing_factor(1)
                     kpt_latt(2) = REAL(2*iy-outer_kp_grid(2)-1, KIND=dp)/(2._dp*REAL(outer_kp_grid(2), KIND=dp)) &
                                   *reducing_factor(2)
                     kpt_latt(3) = REAL(2*iz-outer_kp_grid(3)-1, KIND=dp)/(2._dp*REAL(outer_kp_grid(3), KIND=dp)) &
                                   *reducing_factor(3)
                     kpoints%xkp(1:3, i) = MATMUL(TRANSPOSE(h_inv), kpt_latt(:))

                     IF (2*ix-outer_kp_grid(1)-1 == 0) THEN
                        kpoints%wkp(i) = single_weight
                     ELSE
                        kpoints%wkp(i) = 2._dp*single_weight
                     END IF

                  END DO
               END DO
            END DO

            end_kp = i

            kpoint_weight_left = kpoint_weight_left-SUM(kpoints%wkp(start_kp:end_kp))

            reducing_factor(1) = reducing_factor(1)/REAL(outer_kp_grid(1))
            reducing_factor(2) = reducing_factor(2)/REAL(outer_kp_grid(2))
            reducing_factor(3) = reducing_factor(3)/REAL(outer_kp_grid(3))

         END DO

         single_weight = kpoint_weight_left/REAL(nkp_inner_grid)

         ! the inner grid
         DO ix = 1, kp_grid(1)
            DO iy = 1, kp_grid(2)
               DO iz = 1, kp_grid(3)

                  ! use time reversal symmetry k<->-k
                  IF (2*ix-kp_grid(1)-1 < 0) CYCLE

                  i = i+1
                  kpt_latt(1) = REAL(2*ix-kp_grid(1)-1, KIND=dp)/(2._dp*REAL(kp_grid(1), KIND=dp))*reducing_factor(1)
                  kpt_latt(2) = REAL(2*iy-kp_grid(2)-1, KIND=dp)/(2._dp*REAL(kp_grid(2), KIND=dp))*reducing_factor(2)
                  kpt_latt(3) = REAL(2*iz-kp_grid(3)-1, KIND=dp)/(2._dp*REAL(kp_grid(3), KIND=dp))*reducing_factor(3)

                  kpoints%xkp(1:3, i) = MATMUL(TRANSPOSE(h_inv), kpt_latt(:))

                  kpoints%wkp(i) = 2._dp*single_weight

               END DO
            END DO
         END DO

         IF (do_extra_kpoints) THEN

            single_weight = kpoint_weight_left/REAL(kp_grid(1)*kp_grid(2)*kp_grid(3)/8)

            DO ix = 1, kp_grid(1)/2
               DO iy = 1, kp_grid(2)/2
                  DO iz = 1, kp_grid(3)/2

                     ! use time reversal symmetry k<->-k
                     IF (2*ix-kp_grid(1)/2-1 < 0) CYCLE

                     i = i+1
                     kpt_latt(1) = REAL(2*ix-kp_grid(1)/2-1, KIND=dp)/(REAL(kp_grid(1), KIND=dp))
                     kpt_latt(2) = REAL(2*iy-kp_grid(2)/2-1, KIND=dp)/(REAL(kp_grid(2), KIND=dp))
                     kpt_latt(3) = REAL(2*iz-kp_grid(3)/2-1, KIND=dp)/(REAL(kp_grid(3), KIND=dp))

                     kpoints%xkp(1:3, i) = MATMUL(TRANSPOSE(h_inv), kpt_latt(:))

                     kpoints%wkp(i) = 2._dp*single_weight

                  END DO
               END DO
            END DO

         END IF

         ! default: no symmetry settings
         ALLOCATE (kpoints%kp_sym(kpoints%nkp))
         DO i = 1, kpoints%nkp
            NULLIFY (kpoints%kp_sym(i)%kpoint_sym)
            CALL kpoint_sym_create(kpoints%kp_sym(i)%kpoint_sym)
         END DO

      ELSE

         CALL create_kp_from_gamma(qs_env, qs_env_kp_Gamma_only)

         CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

         CALL calculate_kp_orbitals(qs_env_kp_Gamma_only, kpoints, "MONKHORST-PACK", nadd=nmo, mp_grid=kp_grid(1:3), &
                                    group_size_ext=para_env%num_pe)

         CALL qs_env_release(qs_env_kp_Gamma_only)

      END IF

   END SUBROUTINE

! **************************************************************************************************
! ...
! **************************************************************************************************
!> \brief ...
!> \param vec_Sigma_c_gw ...
!> \param Eigenval_DFT ...
!> \param eps_eigenval ...
! **************************************************************************************************
   SUBROUTINE average_degenerate_levels(vec_Sigma_c_gw, Eigenval_DFT, eps_eigenval)
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vec_Sigma_c_gw
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval_DFT
      REAL(KIND=dp)                                      :: eps_eigenval

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: avg_self_energy
      INTEGER :: degeneracy, first_degenerate_level, i_deg_level, i_level_gw, j_deg_level, jquad, &
         num_deg_levels, num_integ_points, num_levels_gw
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: list_degenerate_levels

      num_levels_gw = SIZE(vec_Sigma_c_gw, 1)

      ALLOCATE (list_degenerate_levels(num_levels_gw))
      list_degenerate_levels = 1

      num_integ_points = SIZE(vec_Sigma_c_gw, 2)

      ALLOCATE (avg_self_energy(num_integ_points))

      DO i_level_gw = 2, num_levels_gw

         IF (ABS(Eigenval_DFT(i_level_gw)-Eigenval_DFT(i_level_gw-1)) < eps_eigenval) THEN

            list_degenerate_levels(i_level_gw) = list_degenerate_levels(i_level_gw-1)

         ELSE

            list_degenerate_levels(i_level_gw) = list_degenerate_levels(i_level_gw-1)+1

         END IF

      END DO

      num_deg_levels = list_degenerate_levels(num_levels_gw)

      DO i_deg_level = 1, num_deg_levels

         degeneracy = 0

         DO i_level_gw = 1, num_levels_gw

            IF (degeneracy == 0 .AND. i_deg_level == list_degenerate_levels(i_level_gw)) THEN

               first_degenerate_level = i_level_gw

            END IF

            IF (i_deg_level == list_degenerate_levels(i_level_gw)) THEN

               degeneracy = degeneracy+1

            END IF

         END DO

         DO jquad = 1, num_integ_points

            avg_self_energy(jquad) = SUM(vec_Sigma_c_gw(first_degenerate_level:first_degenerate_level+degeneracy-1, jquad)) &
                                     /REAL(degeneracy)

         END DO

         DO j_deg_level = 0, degeneracy-1

            vec_Sigma_c_gw(first_degenerate_level+j_deg_level, :) = avg_self_energy(:)

         END DO

      END DO

   END SUBROUTINE

END MODULE rpa_ri_gpw
