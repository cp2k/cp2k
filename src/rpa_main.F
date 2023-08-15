!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
MODULE rpa_main
   USE bibliography,                    ONLY: &
        Bates2013, DelBen2013, DelBen2015, Ren2011, Ren2013, Wilhelm2016a, Wilhelm2016b, &
        Wilhelm2017, Wilhelm2018, cite_reference
   USE bse,                             ONLY: do_subspace_iterations,&
                                              mult_B_with_W_and_fill_local_3c_arrays
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_clear,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dbt_api,                         ONLY: dbt_type
   USE dgemm_counter_types,             ONLY: dgemm_counter_init,&
                                              dgemm_counter_type,&
                                              dgemm_counter_write
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              maxsize,&
                                              release_group_dist
   USE hfx_types,                       ONLY: block_ind_type,&
                                              hfx_compression_type
   USE input_constants,                 ONLY: wfc_mm_style_gemm
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_memory
   USE mathconstants,                   ONLY: pi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE minimax_exp,                     ONLY: check_exp_minimax_range
   USE mp2_grids,                       ONLY: get_clenshaw_grid,&
                                              get_minimax_grid
   USE mp2_laplace,                     ONLY: SOS_MP2_postprocessing
   USE mp2_ri_grad_util,                ONLY: array2fm
   USE mp2_types,                       ONLY: mp2_type,&
                                              three_dim_real_array,&
                                              two_dim_int_array,&
                                              two_dim_real_array
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_axk,                         ONLY: compute_axk_ener
   USE rpa_grad,                        ONLY: rpa_grad_copy_Q,&
                                              rpa_grad_create,&
                                              rpa_grad_finalize,&
                                              rpa_grad_matrix_operations,&
                                              rpa_grad_needed_mem,&
                                              rpa_grad_type
   USE rpa_gw,                          ONLY: allocate_matrices_gw,&
                                              allocate_matrices_gw_im_time,&
                                              compute_GW_self_energy,&
                                              compute_QP_energies,&
                                              compute_W_cubic_GW,&
                                              deallocate_matrices_gw,&
                                              deallocate_matrices_gw_im_time,&
                                              get_fermi_level_offset
   USE rpa_gw_ic,                       ONLY: calculate_ic_correction
   USE rpa_gw_kpoints_util,             ONLY: get_bandstruc_and_k_dependent_MOs,&
                                              invert_eps_compute_W_and_Erpa_kp
   USE rpa_im_time,                     ONLY: compute_mat_P_omega,&
                                              zero_mat_P_omega
   USE rpa_im_time_force_methods,       ONLY: calc_laplace_loop_forces,&
                                              calc_post_loop_forces,&
                                              calc_rpa_loop_forces,&
                                              init_im_time_forces,&
                                              keep_initial_quad
   USE rpa_im_time_force_types,         ONLY: im_time_force_release,&
                                              im_time_force_type
   USE rpa_util,                        ONLY: Q_trace_and_add_unit_matrix,&
                                              alloc_im_time,&
                                              calc_mat_Q,&
                                              compute_Erpa_by_freq_int,&
                                              contract_P_omega_with_mat_L,&
                                              dealloc_im_time,&
                                              remove_scaling_factor_rpa
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_main'

   PUBLIC :: rpa_ri_compute_en

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
!> \param gd_array ...
!> \param gd_B_virtual ...
!> \param gd_B_all ...
!> \param gd_B_occ_bse ...
!> \param gd_B_virt_bse ...
!> \param mo_coeff ...
!> \param fm_matrix_PQ ...
!> \param fm_matrix_L_kpoints ...
!> \param fm_matrix_Minv_L_kpoints ...
!> \param fm_matrix_Minv ...
!> \param fm_matrix_Minv_Vtrunc_Minv ...
!> \param kpoints ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param dimen_RI_red ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param unit_nr ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param my_do_gw ...
!> \param do_im_time ...
!> \param do_bse ...
!> \param matrix_s ...
!> \param mat_munu ...
!> \param mat_P_global ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param starts_array_mc_block ...
!> \param ends_array_mc_block ...
!> \param calc_forces ...
! **************************************************************************************************
   SUBROUTINE rpa_ri_compute_en(qs_env, Erpa, mp2_env, BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, &
                                para_env, para_env_sub, color_sub, &
                                gd_array, gd_B_virtual, gd_B_all, gd_B_occ_bse, gd_B_virt_bse, &
                                mo_coeff, fm_matrix_PQ, fm_matrix_L_kpoints, fm_matrix_Minv_L_kpoints, &
                                fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv, kpoints, &
                                Eigenval, nmo, homo, dimen_RI, dimen_RI_red, gw_corr_lev_occ, gw_corr_lev_virt, &
                                unit_nr, do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_bse, matrix_s, &
                                mat_munu, mat_P_global, &
                                t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                starts_array_mc, ends_array_mc, &
                                starts_array_mc_block, ends_array_mc_block, calc_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: Erpa
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(three_dim_real_array), DIMENSION(:), &
         INTENT(INOUT)                                   :: BIb_C, BIb_C_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: BIb_C_bse_ij, BIb_C_bse_ab
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(INOUT)                             :: color_sub
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_array
      TYPE(group_dist_d1_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: gd_B_virtual
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_B_all, gd_B_occ_bse, gd_B_virt_bse
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_matrix_PQ
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_matrix_L_kpoints, &
                                                            fm_matrix_Minv_L_kpoints, &
                                                            fm_matrix_Minv, &
                                                            fm_matrix_Minv_Vtrunc_Minv
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Eigenval
      INTEGER, INTENT(IN)                                :: nmo
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ, gw_corr_lev_virt
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2, my_do_gw, &
                                                            do_im_time, do_bse
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global
      TYPE(dbt_type)                                     :: t_3c_M
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_O
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_compressed
      TYPE(block_ind_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_ind
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_ri_compute_en'

      INTEGER :: best_integ_group_size, best_num_integ_point, color_rpa_group, dimen_homo_square, &
         dimen_nm_gw, dimen_virt_square, handle, handle2, handle3, ierr, iiB, &
         input_num_integ_groups, integ_group_size, ispin, jjB, min_integ_group_size, &
         my_ab_comb_bse_end, my_ab_comb_bse_size, my_ab_comb_bse_start, my_group_L_end, &
         my_group_L_size, my_group_L_start, my_ij_comb_bse_end, my_ij_comb_bse_size, &
         my_ij_comb_bse_start, my_nm_gw_end, my_nm_gw_size, my_nm_gw_start, ncol_block_mat, &
         ngroup, nrow_block_mat, nspins, num_integ_group, num_integ_points, pos_integ_group
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dimen_ia, my_ia_end, my_ia_size, &
                                                            my_ia_start, virtual
      LOGICAL                                            :: do_kpoints_from_Gamma, do_minimax_quad, &
                                                            my_open_shell, skip_integ_group_opt
      REAL(KIND=dp) :: allowed_memory, avail_mem, E_Range, Emax, Emin, mem_for_iaK, mem_for_QK, &
         mem_min, mem_per_group, mem_per_rank, mem_per_repl, mem_real
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Eigenval_kp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_mat_Q, fm_mat_Q_gemm, fm_mat_S, &
                                                            fm_mat_S_gw
      TYPE(cp_fm_type), DIMENSION(1)                     :: fm_mat_R_gw, fm_mat_S_ab_bse, &
                                                            fm_mat_S_ij_bse
      TYPE(mp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:)                                    :: BIb_C_2D, BIb_C_2D_gw
      TYPE(two_dim_real_array), DIMENSION(1)             :: BIb_C_2D_bse_ab, BIb_C_2D_bse_ij

      CALL timeset(routineN, handle)

      CALL cite_reference(DelBen2013)
      CALL cite_reference(DelBen2015)

      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         CALL cite_reference(Bates2013)
      END IF
      IF (mp2_env%ri_rpa%do_rse) THEN
         CALL cite_reference(Ren2011)
         CALL cite_reference(Ren2013)
      END IF

      IF (my_do_gw) THEN
         CALL cite_reference(Wilhelm2016a)
         CALL cite_reference(Wilhelm2017)
         CALL cite_reference(Wilhelm2018)
      END IF

      IF (do_im_time) THEN
         CALL cite_reference(Wilhelm2016b)
      END IF

      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)
      ALLOCATE (virtual(nspins), dimen_ia(nspins), my_ia_end(nspins), my_ia_start(nspins), my_ia_size(nspins))
      virtual(:) = nmo - homo(:)
      dimen_ia(:) = virtual(:)*homo(:)

      ALLOCATE (Eigenval_kp(nmo, 1, nspins))
      Eigenval_kp(:, 1, :) = Eigenval(:, :)

      IF (do_im_time) mp2_env%ri_rpa%minimax_quad = .TRUE.
      do_minimax_quad = mp2_env%ri_rpa%minimax_quad

      IF (do_ri_sos_laplace_mp2) THEN
         num_integ_points = mp2_env%ri_laplace%n_quadrature
         input_num_integ_groups = mp2_env%ri_laplace%num_integ_groups

         ! check the range for the minimax approximation
         E_Range = mp2_env%e_range
         IF (mp2_env%e_range <= 1.0_dp .OR. mp2_env%e_gap <= 0.0_dp) THEN
            Emin = HUGE(dp)
            Emax = 0.0_dp
            DO ispin = 1, nspins
               IF (homo(ispin) > 0) THEN
                  Emin = MIN(Emin, 2.0_dp*(Eigenval(homo(ispin) + 1, ispin) - Eigenval(homo(ispin), ispin)))
                  Emax = MAX(Emax, 2.0_dp*(MAXVAL(Eigenval(:, ispin)) - MINVAL(Eigenval(:, ispin))))
               END IF
            END DO
            E_Range = Emax/Emin
         END IF
         IF (E_Range < 2.0_dp) E_Range = 2.0_dp
         ierr = 0
         CALL check_exp_minimax_range(num_integ_points, E_Range, ierr)
         IF (ierr /= 0) THEN
            jjB = num_integ_points - 1
            DO iiB = 1, jjB
               num_integ_points = num_integ_points - 1
               ierr = 0
               CALL check_exp_minimax_range(num_integ_points, E_Range, ierr)
               IF (ierr == 0) EXIT
            END DO
         END IF
         CPASSERT(num_integ_points >= 1)
      ELSE
         num_integ_points = mp2_env%ri_rpa%rpa_num_quad_points
         input_num_integ_groups = mp2_env%ri_rpa%rpa_num_integ_groups

         IF (my_do_gw .AND. do_minimax_quad) THEN
            IF (num_integ_points > 34) THEN
               IF (unit_nr > 0) &
                  CALL cp_warn(__LOCATION__, &
                               "The required number of quadrature point exceeds the maximum possible in the "// &
                               "Minimax quadrature scheme. The number of quadrature point has been reset to 30.")
               num_integ_points = 30
            END IF
         ELSE
            IF (do_minimax_quad .AND. num_integ_points > 20) THEN
               IF (unit_nr > 0) &
                  CALL cp_warn(__LOCATION__, &
                               "The required number of quadrature point exceeds the maximum possible in the "// &
                               "Minimax quadrature scheme. The number of quadrature point has been reset to 20.")
               num_integ_points = 20
            END IF
         END IF
      END IF
      allowed_memory = mp2_env%mp2_memory

      CALL get_group_dist(gd_array, color_sub, my_group_L_start, my_group_L_end, my_group_L_size)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! for imaginary time or periodic GW or BSE, we use all processors for a single frequency/time point
      IF (do_im_time .OR. mp2_env%ri_g0w0%do_periodic .OR. do_bse) THEN

         integ_group_size = ngroup
         best_num_integ_point = num_integ_points

      ELSE

         ! Calculate available memory and create integral group according to that
         ! mem_for_iaK is the memory needed for storing the 3 centre integrals
         mem_for_iaK = REAL(SUM(dimen_ia), KIND=dp)*dimen_RI_red*8.0_dp/(1024_dp**2)
         mem_for_QK = REAL(dimen_RI_red, KIND=dp)*nspins*dimen_RI_red*8.0_dp/(1024_dp**2)

         CALL m_memory(mem)
         mem_real = (mem + 1024*1024 - 1)/(1024*1024)
         CALL para_env%min(mem_real)

         mem_per_rank = 0.0_dp

         ! B_ia_P
         mem_per_repl = mem_for_iaK
         ! Q (regular and for dgemm)
         mem_per_repl = mem_per_repl + 2.0_dp*mem_for_QK

         IF (calc_forces) CALL rpa_grad_needed_mem(homo, virtual, dimen_RI_red, mem_per_rank, mem_per_repl, do_ri_sos_laplace_mp2)

         mem_min = mem_per_repl/para_env%num_pe + mem_per_rank

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Minimum required memory per MPI process:', mem_min, ' MiB'
            WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Available memory per MPI process:', mem_real, ' MiB'
         END IF

         ! Use only the allowed amount of memory
         mem_real = MIN(mem_real, allowed_memory)
         ! For the memory estimate, we require the amount of required memory per replication group and the available memory
         mem_real = mem_real - mem_per_rank

         mem_per_group = mem_real*para_env_sub%num_pe

         ! here we try to find the best rpa/laplace group size
         skip_integ_group_opt = .FALSE.

         ! Check the input number of integration groups
         IF (input_num_integ_groups > 0) THEN
            IF (MOD(num_integ_points, input_num_integ_groups) == 0) THEN
               IF (MOD(ngroup, input_num_integ_groups) == 0) THEN
                  best_integ_group_size = ngroup/input_num_integ_groups
                  best_num_integ_point = num_integ_points/input_num_integ_groups
                  skip_integ_group_opt = .TRUE.
               ELSE
                  IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Total number of groups not multiple of NUM_INTEG_GROUPS'
               END IF
            ELSE
               IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'NUM_QUAD_POINTS not multiple of NUM_INTEG_GROUPS'
            END IF
         END IF

         IF (.NOT. skip_integ_group_opt) THEN
            best_integ_group_size = ngroup
            best_num_integ_point = num_integ_points

            min_integ_group_size = MAX(1, ngroup/num_integ_points)

            integ_group_size = min_integ_group_size - 1
            DO iiB = min_integ_group_size + 1, ngroup
               integ_group_size = integ_group_size + 1

               ! check that the ngroup is a multiple of integ_group_size
               IF (MOD(ngroup, integ_group_size) /= 0) CYCLE

               ! check for memory
               avail_mem = integ_group_size*mem_per_group
               IF (avail_mem < mem_per_repl) CYCLE

               ! check that the integration groups have the same size
               num_integ_group = ngroup/integ_group_size
               IF (MOD(num_integ_points, num_integ_group) /= 0) CYCLE

               best_num_integ_point = num_integ_points/num_integ_group
               best_integ_group_size = integ_group_size

               EXIT

            END DO
         END IF

         integ_group_size = best_integ_group_size

      END IF

      IF (unit_nr > 0 .AND. .NOT. do_im_time) THEN
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

      CALL timeset(routineN//"_reorder", handle2)

      ! not necessary for imaginary time

      ALLOCATE (BIb_C_2D(nspins))
      IF (.NOT. do_im_time) THEN

         ! reorder the local data in such a way to help the next stage of matrix creation
         ! now the data inside the group are divided into a ia x K matrix
         DO ispin = 1, nspins
            CALL calculate_BIb_C_2D(BIb_C_2D(ispin)%array, BIb_C(ispin)%array, para_env_sub, dimen_ia(ispin), &
                                    homo(ispin), virtual(ispin), gd_B_virtual(ispin), &
                                    my_ia_size(ispin), my_ia_start(ispin), my_ia_end(ispin), my_group_L_size)

            DEALLOCATE (BIb_C(ispin)%array)
            CALL release_group_dist(gd_B_virtual(ispin))

         END DO

         ! in the GW case, BIb_C_2D_gw is an nm x K matrix, with n: number of corr GW levels, m=nmo
         IF (my_do_gw) THEN
            ALLOCATE (BIb_C_2D_gw(nspins))

            CALL timeset(routineN//"_reorder_gw", handle3)

            dimen_nm_gw = nmo*(gw_corr_lev_occ(1) + gw_corr_lev_virt(1))

            ! The same for open shell
            DO ispin = 1, nspins
               CALL calculate_BIb_C_2D(BIb_C_2D_gw(ispin)%array, BIb_C_gw(ispin)%array, para_env_sub, dimen_nm_gw, &
                                       gw_corr_lev_occ(ispin) + gw_corr_lev_virt(ispin), nmo, gd_B_all, &
                                       my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, my_group_L_size)
               DEALLOCATE (BIb_C_gw(ispin)%array)
            END DO

            CALL release_group_dist(gd_B_all)

            CALL timestop(handle3)

         END IF
      END IF

      IF (do_bse) THEN

         CALL timeset(routineN//"_reorder_bse1", handle3)

         dimen_homo_square = homo(1)**2

         CALL calculate_BIb_C_2D(BIb_C_2D_bse_ij(1)%array, BIb_C_bse_ij, para_env_sub, dimen_homo_square, &
                                 homo(1), homo(1), gd_B_occ_bse, &
                                 my_ij_comb_bse_size, my_ij_comb_bse_start, my_ij_comb_bse_end, my_group_L_size)

         DEALLOCATE (BIb_C_bse_ij)
         CALL release_group_dist(gd_B_occ_bse)

         CALL timestop(handle3)

         CALL timeset(routineN//"_reorder_bse2", handle3)

         dimen_virt_square = virtual(1)**2

         CALL calculate_BIb_C_2D(BIb_C_2D_bse_ab(1)%array, BIb_C_bse_ab, para_env_sub, dimen_virt_square, &
                                 virtual(1), virtual(1), gd_B_virt_bse, &
                                 my_ab_comb_bse_size, my_ab_comb_bse_start, my_ab_comb_bse_end, my_group_L_size)

         DEALLOCATE (BIb_C_bse_ab)
         CALL release_group_dist(gd_B_virt_bse)

         CALL timestop(handle3)

      END IF

      CALL timestop(handle2)

      IF (num_integ_group > 1) THEN
         ALLOCATE (para_env_RPA)
         CALL para_env_RPA%from_split(para_env, color_rpa_group)
      ELSE
         para_env_RPA => para_env
      END IF

      ! now create the matrices needed for the calculation, Q, S and G
      ! Q and G will have omega dependence

      IF (do_im_time) THEN
         ALLOCATE (fm_mat_Q(nspins), fm_mat_Q_gemm(1), fm_mat_S(1))
      ELSE
         ALLOCATE (fm_mat_Q(nspins), fm_mat_Q_gemm(nspins), fm_mat_S(nspins))
      END IF

      CALL create_integ_mat(BIb_C_2D, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                            dimen_RI_red, dimen_ia, color_rpa_group, &
                            mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                            my_ia_size, my_ia_start, my_ia_end, &
                            my_group_L_size, my_group_L_start, my_group_L_end, &
                            para_env_RPA, fm_mat_S, nrow_block_mat, ncol_block_mat, &
                            dimen_ia_for_block_size=dimen_ia(1), &
                            do_im_time=do_im_time, fm_mat_Q_gemm=fm_mat_Q_gemm, fm_mat_Q=fm_mat_Q, qs_env=qs_env)

      DEALLOCATE (BIb_C_2D, my_ia_end, my_ia_size, my_ia_start)

      ! for GW, we need other matrix fm_mat_S, always allocate the container to prevent crying compilers
      ALLOCATE (fm_mat_S_gw(nspins))
      IF (my_do_gw .AND. .NOT. do_im_time) THEN

         CALL create_integ_mat(BIb_C_2D_gw, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, [dimen_nm_gw, dimen_nm_gw], color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               [my_nm_gw_size, my_nm_gw_size], [my_nm_gw_start, my_nm_gw_start], [my_nm_gw_end, my_nm_gw_end], &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_gw, nrow_block_mat, ncol_block_mat, &
                               fm_mat_Q(1)%matrix_struct%context, fm_mat_Q(1)%matrix_struct%context, &
                               fm_mat_Q=fm_mat_R_gw)
         DEALLOCATE (BIb_C_2D_gw)

      END IF

      ! for Bethe-Salpeter, we need other matrix fm_mat_S
      IF (do_bse) THEN

         CALL create_integ_mat(BIb_C_2D_bse_ij, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, [dimen_homo_square], color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               [my_ij_comb_bse_size], [my_ij_comb_bse_start], [my_ij_comb_bse_end], &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_ij_bse, nrow_block_mat, ncol_block_mat, &
                               fm_mat_Q(1)%matrix_struct%context, fm_mat_Q(1)%matrix_struct%context)

         CALL create_integ_mat(BIb_C_2D_bse_ab, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, [dimen_virt_square], color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               [my_ab_comb_bse_size], [my_ab_comb_bse_start], [my_ab_comb_bse_end], &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_ab_bse, nrow_block_mat, ncol_block_mat, &
                               fm_mat_Q(1)%matrix_struct%context, fm_mat_Q(1)%matrix_struct%context)

      END IF

      do_kpoints_from_Gamma = qs_env%mp2_env%ri_rpa_im_time%do_kpoints_from_Gamma
      IF (do_kpoints_from_Gamma) THEN
         CALL get_bandstruc_and_k_dependent_MOs(qs_env, Eigenval_kp)
      END IF

      ! Now start the RPA calculation
      ! fm_mo_coeff_occ, fm_mo_coeff_virt will be deallocated here
      CALL rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                       homo, virtual, dimen_RI, dimen_RI_red, dimen_ia, dimen_nm_gw, &
                       Eigenval_kp, num_integ_points, num_integ_group, color_rpa_group, &
                       fm_matrix_PQ, fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw(1), &
                       fm_mat_S_ij_bse(1), fm_mat_S_ab_bse(1), &
                       my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                       do_minimax_quad, &
                       do_im_time, mo_coeff, &
                       fm_matrix_L_kpoints, fm_matrix_Minv_L_kpoints, &
                       fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv, &
                       mat_munu, mat_P_global, &
                       t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                       starts_array_mc, ends_array_mc, &
                       starts_array_mc_block, ends_array_mc_block, &
                       matrix_s, do_kpoints_from_Gamma, kpoints, gd_array, color_sub, &
                       do_ri_sos_laplace_mp2=do_ri_sos_laplace_mp2, calc_forces=calc_forces)

      CALL release_group_dist(gd_array)

      IF (num_integ_group > 1) CALL mp_para_env_release(para_env_RPA)

      IF (.NOT. do_im_time) THEN
         CALL cp_fm_release(fm_mat_Q_gemm)
         CALL cp_fm_release(fm_mat_S)
      END IF
      CALL cp_fm_release(fm_mat_Q)

      IF (my_do_gw .AND. .NOT. do_im_time) THEN
         CALL cp_fm_release(fm_mat_S_gw)
         CALL cp_fm_release(fm_mat_R_gw(1))
      END IF

      IF (do_bse) THEN
         CALL cp_fm_release(fm_mat_S_ij_bse(1))
         CALL cp_fm_release(fm_mat_S_ab_bse(1))
      END IF

      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         CALL dbcsr_release(mp2_env%ri_rpa%mo_coeff_o)
         DEALLOCATE (mp2_env%ri_rpa%mo_coeff_o)
         CALL dbcsr_release(mp2_env%ri_rpa%mo_coeff_v)
         DEALLOCATE (mp2_env%ri_rpa%mo_coeff_v)
      END IF

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
!> \param gd_B_virtual ...
!> \param my_ia_size ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \author Jan Wilhelm, 03/2015
! **************************************************************************************************
   SUBROUTINE calculate_BIb_C_2D(BIb_C_2D, BIb_C, para_env_sub, dimen_ia, homo, virtual, &
                                 gd_B_virtual, &
                                 my_ia_size, my_ia_start, my_ia_end, my_group_L_size)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: BIb_C_2D
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: BIb_C
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_sub
      INTEGER, INTENT(IN)                                :: dimen_ia, homo, virtual
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_B_virtual
      INTEGER                                            :: my_ia_size, my_ia_start, my_ia_end, &
                                                            my_group_L_size

      INTEGER, PARAMETER                                 :: occ_chunk = 128

      INTEGER :: ia_global, iiB, itmp(2), jjB, my_B_size, my_B_virtual_start, occ_high, occ_low, &
         proc_receive, proc_send, proc_shift, rec_B_size, rec_B_virtual_end, rec_B_virtual_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: BIb_C_rec_1D
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: BIb_C_rec

      itmp = get_limit(dimen_ia, para_env_sub%num_pe, para_env_sub%mepos)
      my_ia_start = itmp(1)
      my_ia_end = itmp(2)
      my_ia_size = my_ia_end - my_ia_start + 1

      CALL get_group_dist(gd_B_virtual, para_env_sub%mepos, sizes=my_B_size, starts=my_B_virtual_start)

      ! reorder data
      ALLOCATE (BIb_C_2D(my_group_L_size, my_ia_size))

!$OMP     PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
!$OMP              SHARED(homo,my_B_size,virtual,my_B_virtual_start,my_ia_start,my_ia_end,BIb_C,BIb_C_2D,&
!$OMP              my_group_L_size)
      DO iiB = 1, homo
         DO jjB = 1, my_B_size
            ia_global = (iiB - 1)*virtual + my_B_virtual_start + jjB - 1
            IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
               BIb_C_2D(1:my_group_L_size, ia_global - my_ia_start + 1) = BIb_C(1:my_group_L_size, jjB, iiB)
            END IF
         END DO
      END DO

      IF (para_env_sub%num_pe > 1) THEN
         ALLOCATE (BIb_C_rec_1D(INT(my_group_L_size, int_8)*maxsize(gd_B_virtual)*MIN(homo, occ_chunk)))
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
            proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

            CALL get_group_dist(gd_B_virtual, proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

            ! do this in chunks to avoid high memory overhead
            DO occ_low = 1, homo, occ_chunk
               occ_high = MIN(homo, occ_low + occ_chunk - 1)
               BIb_C_rec(1:my_group_L_size, 1:rec_B_size, 1:occ_high - occ_low + 1) => &
                  BIb_C_rec_1D(1:INT(my_group_L_size, int_8)*rec_B_size*(occ_high - occ_low + 1))
               CALL para_env_sub%sendrecv(BIb_C(:, :, occ_low:occ_high), proc_send, &
                                          BIb_C_rec(:, :, 1:occ_high - occ_low + 1), proc_receive)
!$OMP          PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
!$OMP                   SHARED(occ_low,occ_high,rec_B_size,virtual,rec_B_virtual_start,my_ia_start,my_ia_end,BIb_C_rec,BIb_C_2D,&
!$OMP                          my_group_L_size)
               DO iiB = occ_low, occ_high
                  DO jjB = 1, rec_B_size
                     ia_global = (iiB - 1)*virtual + rec_B_virtual_start + jjB - 1
                     IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
                     BIb_C_2D(1:my_group_L_size, ia_global - my_ia_start + 1) = BIb_C_rec(1:my_group_L_size, jjB, iiB - occ_low + 1)
                     END IF
                  END DO
               END DO
            END DO

         END DO
         DEALLOCATE (BIb_C_rec_1D)
      END IF

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
!> \param nrow_block_mat ...
!> \param ncol_block_mat ...
!> \param blacs_env_ext ...
!> \param blacs_env_ext_S ...
!> \param dimen_ia_for_block_size ...
!> \param do_im_time ...
!> \param fm_mat_Q_gemm ...
!> \param fm_mat_Q ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE create_integ_mat(BIb_C_2D, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_ia, color_rpa_group, &
                               ext_row_block_size, ext_col_block_size, unit_nr, &
                               my_ia_size, my_ia_start, my_ia_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S, nrow_block_mat, ncol_block_mat, &
                               blacs_env_ext, blacs_env_ext_S, dimen_ia_for_block_size, &
                               do_im_time, fm_mat_Q_gemm, fm_mat_Q, qs_env)

      TYPE(two_dim_real_array), DIMENSION(:), &
         INTENT(INOUT)                                   :: BIb_C_2D
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: color_sub, ngroup, integ_group_size, &
                                                            dimen_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dimen_ia
      INTEGER, INTENT(IN)                                :: color_rpa_group, ext_row_block_size, &
                                                            ext_col_block_size, unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_ia_size, my_ia_start, my_ia_end
      INTEGER, INTENT(IN)                                :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env_RPA
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: fm_mat_S
      INTEGER, INTENT(INOUT)                             :: nrow_block_mat, ncol_block_mat
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext, blacs_env_ext_S
      INTEGER, INTENT(IN), OPTIONAL                      :: dimen_ia_for_block_size
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_im_time
      TYPE(cp_fm_type), DIMENSION(:), OPTIONAL           :: fm_mat_Q_gemm, fm_mat_Q
      TYPE(qs_environment_type), INTENT(IN), OPTIONAL, &
         POINTER                                         :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_integ_mat'

      INTEGER                                            :: col_row_proc_ratio, grid_2D(2), handle, &
                                                            iproc, iproc_col, iproc_row, ispin, &
                                                            mepos_in_RPA_group
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: group_grid_2_mepos
      LOGICAL                                            :: my_blacs_ext, my_blacs_S_ext, &
                                                            my_do_im_time
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_Q
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(group_dist_d1_type)                           :: gd_ia, gd_L

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(blacs_env_ext) .OR. PRESENT(dimen_ia_for_block_size))

      my_blacs_ext = .FALSE.
      IF (PRESENT(blacs_env_ext)) my_blacs_ext = .TRUE.

      my_blacs_S_ext = .FALSE.
      IF (PRESENT(blacs_env_ext_S)) my_blacs_S_ext = .TRUE.

      my_do_im_time = .FALSE.
      IF (PRESENT(do_im_time)) my_do_im_time = do_im_time

      NULLIFY (blacs_env)
      ! create the RPA blacs env
      IF (my_blacs_S_ext) THEN
         blacs_env => blacs_env_ext_S
      ELSE
         IF (para_env_RPA%num_pe > 1) THEN
            col_row_proc_ratio = MAX(1, dimen_ia_for_block_size/dimen_RI)

            iproc_col = MIN(MAX(INT(SQRT(REAL(para_env_RPA%num_pe*col_row_proc_ratio, KIND=dp))), 1), para_env_RPA%num_pe) + 1
            DO iproc = 1, para_env_RPA%num_pe
               iproc_col = iproc_col - 1
               IF (MOD(para_env_RPA%num_pe, iproc_col) == 0) EXIT
            END DO

            iproc_row = para_env_RPA%num_pe/iproc_col
            grid_2D(1) = iproc_row
            grid_2D(2) = iproc_col
         ELSE
            grid_2D = 1
         END IF
         CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_RPA, grid_2d=grid_2d)

         IF (unit_nr > 0 .AND. .NOT. my_do_im_time) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Number row processes:", grid_2D(1)
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Number column processes:", grid_2D(2)
         END IF

         ! define the block_size for the row
         IF (ext_row_block_size > 0) THEN
            nrow_block_mat = ext_row_block_size
         ELSE
            nrow_block_mat = MAX(1, dimen_RI/grid_2D(1)/2)
         END IF

         ! define the block_size for the column
         IF (ext_col_block_size > 0) THEN
            ncol_block_mat = ext_col_block_size
         ELSE
            ncol_block_mat = MAX(1, dimen_ia_for_block_size/grid_2D(2)/2)
         END IF

         IF (unit_nr > 0 .AND. .NOT. my_do_im_time) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Row block size:", nrow_block_mat
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Column block size:", ncol_block_mat
         END IF
      END IF

      IF (.NOT. my_do_im_time) THEN
         DO ispin = 1, SIZE(BIb_C_2D)
            NULLIFY (fm_struct)
            IF (my_blacs_ext) THEN
               CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                                        ncol_global=dimen_ia(ispin), para_env=para_env_RPA)
            ELSE
               CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                                        ncol_global=dimen_ia(ispin), para_env=para_env_RPA, &
                                        nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)

            END IF ! external blacs_env

            CALL create_group_dist(gd_ia, my_ia_start(ispin), my_ia_end(ispin), my_ia_size(ispin), para_env_RPA)
            CALL create_group_dist(gd_L, my_group_L_start, my_group_L_end, my_group_L_size, para_env_RPA)

            ! create the info array

            mepos_in_RPA_group = MOD(color_sub, integ_group_size)
            ALLOCATE (group_grid_2_mepos(0:integ_group_size - 1, 0:para_env_sub%num_pe - 1))
            group_grid_2_mepos = 0
            group_grid_2_mepos(mepos_in_RPA_group, para_env_sub%mepos) = para_env_RPA%mepos
            CALL para_env_RPA%sum(group_grid_2_mepos)

            CALL array2fm(BIb_C_2D(ispin)%array, fm_struct, my_group_L_start, my_group_L_end, &
                          my_ia_start(ispin), my_ia_end(ispin), gd_L, gd_ia, &
                          group_grid_2_mepos, ngroup, para_env_sub%num_pe, fm_mat_S(ispin), &
                          integ_group_size, color_rpa_group)

            DEALLOCATE (group_grid_2_mepos)
            CALL cp_fm_struct_release(fm_struct)

            ! deallocate the info array
            CALL release_group_dist(gd_L)
            CALL release_group_dist(gd_ia)

            ! sum the local data across processes belonging to different RPA group.
            IF (para_env_RPA%num_pe /= para_env%num_pe) THEN
               BLOCK
                  TYPE(mp_comm_type) :: comm_exchange
                  comm_exchange = fm_mat_S(ispin)%matrix_struct%context%interconnect(para_env)
                  CALL comm_exchange%sum(fm_mat_S(ispin)%local_data)
                  CALL comm_exchange%free()
               END BLOCK
            END IF
         END DO
      END IF

      IF (PRESENT(fm_mat_Q_gemm) .AND. .NOT. my_do_im_time) THEN
         ! create the Q matrix dimen_RIxdimen_RI where the result of the mat-mat-mult will be stored
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA, &
                                  nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)
         DO ispin = 1, SIZE(fm_mat_Q_gemm)
            CALL cp_fm_create(fm_mat_Q_gemm(ispin), fm_struct, name="fm_mat_Q_gemm")
         END DO
         CALL cp_fm_struct_release(fm_struct)
      END IF

      IF (PRESENT(fm_mat_Q)) THEN
         NULLIFY (blacs_env_Q)
         IF (my_blacs_ext) THEN
            blacs_env_Q => blacs_env_ext
         ELSE IF (para_env_RPA%num_pe == para_env%num_pe .AND. PRESENT(qs_env)) THEN
            CALL get_qs_env(qs_env, blacs_env=blacs_env_Q)
         ELSE
            CALL cp_blacs_env_create(blacs_env=blacs_env_Q, para_env=para_env_RPA)
         END IF
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_Q, nrow_global=dimen_RI, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA)
         DO ispin = 1, SIZE(fm_mat_Q)
            CALL cp_fm_create(fm_mat_Q(ispin), fm_struct, name="fm_mat_Q")
         END DO

         CALL cp_fm_struct_release(fm_struct)

         IF (.NOT. (my_blacs_ext .OR. (para_env_RPA%num_pe == para_env%num_pe .AND. PRESENT(qs_env)))) &
            CALL cp_blacs_env_release(blacs_env_Q)
      END IF

      ! release blacs_env
      IF (.NOT. my_blacs_S_ext) CALL cp_blacs_env_release(blacs_env)

      CALL timestop(handle)

   END SUBROUTINE create_integ_mat

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
!> \param dimen_RI_red ...
!> \param dimen_ia ...
!> \param dimen_nm_gw ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param fm_matrix_PQ ...
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
!> \param mo_coeff ...
!> \param fm_matrix_L_kpoints ...
!> \param fm_matrix_Minv_L_kpoints ...
!> \param fm_matrix_Minv ...
!> \param fm_matrix_Minv_Vtrunc_Minv ...
!> \param mat_munu ...
!> \param mat_P_global ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param starts_array_mc_block ...
!> \param ends_array_mc_block ...
!> \param matrix_s ...
!> \param do_kpoints_from_Gamma ...
!> \param kpoints ...
!> \param gd_array ...
!> \param color_sub ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param calc_forces ...
! **************************************************************************************************
   SUBROUTINE rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                          homo, virtual, dimen_RI, dimen_RI_red, dimen_ia, dimen_nm_gw, &
                          Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                          fm_matrix_PQ, fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                          fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                          my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                          do_minimax_quad, do_im_time, mo_coeff, &
                          fm_matrix_L_kpoints, fm_matrix_Minv_L_kpoints, &
                          fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv, &
                          mat_munu, mat_P_global, &
                          t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                          starts_array_mc, ends_array_mc, &
                          starts_array_mc_block, ends_array_mc_block, &
                          matrix_s, do_kpoints_from_Gamma, kpoints, gd_array, color_sub, &
                          do_ri_sos_laplace_mp2, calc_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: Erpa
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_RPA, para_env_sub
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dimen_ia
      INTEGER, INTENT(IN)                                :: dimen_nm_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_matrix_PQ
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                                                            fm_mat_S_gw
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_R_gw, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      LOGICAL, INTENT(IN)                                :: my_do_gw, do_bse
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, INTENT(IN)                                :: do_minimax_quad, do_im_time
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_matrix_L_kpoints, &
                                                            fm_matrix_Minv_L_kpoints, &
                                                            fm_matrix_Minv, &
                                                            fm_matrix_Minv_Vtrunc_Minv
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_M
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: t_3c_O
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_compressed
      TYPE(block_ind_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_ind
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      LOGICAL                                            :: do_kpoints_from_Gamma
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2, calc_forces

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_num_int'

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: vec_Sigma_c_gw
      INTEGER :: count_ev_sc_GW, cut_memory, group_size_P, gw_corr_lev_tot, handle, handle3, i, &
         ikp_local, ispin, iter_evGW, iter_sc_GW0, j, jquad, max_iter_bse, min_bsize, mm_style, &
         nkp, nkp_self_energy, nmo, nspins, num_3c_repl, num_cells_dm, num_fit_points, &
         num_Z_vectors, Pspin, Qspin, size_P
      INTEGER(int_8)                                     :: dbcsr_nflop
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index_3c
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, prim_blk_sizes, &
                                                            RI_blk_sizes
      LOGICAL :: do_apply_ic_corr_to_gw, do_gw_im_time, do_ic_model, do_kpoints_cubic_RPA, &
         do_periodic, do_print, do_ri_Sigma_x, exit_ev_gw, first_cycle, &
         first_cycle_periodic_correction, my_open_shell, print_ic_values
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :, :)     :: has_mat_P_blocks
      REAL(KIND=dp) :: a_scaling, alpha, dbcsr_time, e_axk, e_axk_corr, eps_filter, &
         eps_filter_im_time, eps_min_trans, ext_scaling, fermi_level_offset, &
         fermi_level_offset_input, my_flop_rate, omega, omega_max_fit, omega_old, tau, tau_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta_corr, e_fermi, tau_tj, tau_wj, tj, &
                                                            trace_Qomega, vec_omega_fit_gw, wj, &
                                                            wkp_W
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: vec_W_gw, weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: B_abQ_bse_local, B_bar_iaQ_bse_local, &
         B_bar_ijQ_bse_local, B_iaQ_bse_local, Eigenval_last, Eigenval_scf, vec_Sigma_x_gw
      TYPE(cp_cfm_type)                                  :: cfm_mat_Q
      TYPE(cp_fm_type) :: fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, fm_mat_RI_global_work, &
         fm_mat_work, fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, fm_scaled_dm_occ_tau, &
         fm_scaled_dm_virt_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_mat_S_gw_work, fm_mat_W, &
                                                            fm_mo_coeff_occ, fm_mo_coeff_virt
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_mat_L_kpoints, fm_mat_Minv_L_kpoints
      TYPE(dbcsr_p_type)                                 :: mat_dm, mat_L, mat_M_P_munu_occ, &
                                                            mat_M_P_munu_virt, mat_MinvVMinv
      TYPE(dbcsr_p_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: mat_P_omega
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_im_mo_mo, &
                                                            matrix_berry_re_mo_mo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega_kp
      TYPE(dbcsr_type), POINTER                          :: mat_W, mat_work
      TYPE(dbt_type)                                     :: t_3c_overl_int_ao_mo
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(dgemm_counter_type)                           :: dgemm_counter
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: t_3c_O_mo_compressed
      TYPE(im_time_force_type)                           :: force_data
      TYPE(rpa_grad_type)                                :: rpa_grad
      TYPE(two_dim_int_array), ALLOCATABLE, DIMENSION(:) :: t_3c_O_mo_ind

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      nmo = homo(1) + virtual(1)
      my_open_shell = (nspins == 2)

      do_gw_im_time = my_do_gw .AND. do_im_time
      do_ri_Sigma_x = mp2_env%ri_g0w0%do_ri_Sigma_x
      do_ic_model = mp2_env%ri_g0w0%do_ic_model
      print_ic_values = mp2_env%ri_g0w0%print_ic_values
      do_periodic = mp2_env%ri_g0w0%do_periodic
      do_kpoints_cubic_RPA = mp2_env%ri_rpa_im_time%do_im_time_kpoints

      ! For SOS-MP2 only gemm is implemented
      mm_style = wfc_mm_style_gemm
      IF (.NOT. do_ri_sos_laplace_mp2) mm_style = mp2_env%ri_rpa%mm_style

      IF (my_do_gw) THEN
         ext_scaling = 0.2_dp
         omega_max_fit = mp2_env%ri_g0w0%omega_max_fit
         fermi_level_offset_input = mp2_env%ri_g0w0%fermi_level_offset
         iter_evGW = mp2_env%ri_g0w0%iter_evGW
         iter_sc_GW0 = mp2_env%ri_g0w0%iter_sc_GW0
         IF ((.NOT. do_im_time)) THEN
            IF (iter_sc_GW0 .NE. 1 .AND. iter_evGW .NE. 1) CPABORT("Mixed scGW0/evGW not implemented.")
            ! in case of scGW0 with the N^4 algorithm, we use the evGW code but use the DFT eigenvalues for W
            IF (iter_sc_GW0 .NE. 1) iter_evGW = iter_sc_GW0
         END IF
      ELSE
         ext_scaling = 0.0_dp
         iter_evGW = 1
         iter_sc_GW0 = 1
      END IF

      IF (do_kpoints_cubic_RPA .AND. do_ri_sos_laplace_mp2) THEN
         CPABORT("RI-SOS-Laplace-MP2 with k-point-sampling is not implemented.")
      END IF

      do_apply_ic_corr_to_gw = .FALSE.
      IF (mp2_env%ri_g0w0%ic_corr_list(1)%array(1) > 0.0_dp) do_apply_ic_corr_to_gw = .TRUE.

      IF (do_im_time) THEN
         CPASSERT(do_minimax_quad .OR. do_ri_sos_laplace_mp2)

         group_size_P = mp2_env%ri_rpa_im_time%group_size_P
         cut_memory = mp2_env%ri_rpa_im_time%cut_memory
         eps_filter = mp2_env%ri_rpa_im_time%eps_filter
         eps_filter_im_time = mp2_env%ri_rpa_im_time%eps_filter* &
                              mp2_env%ri_rpa_im_time%eps_filter_factor

         min_bsize = mp2_env%ri_rpa_im_time%min_bsize

         CALL alloc_im_time(qs_env, para_env, dimen_RI, dimen_RI_red, &
                            num_integ_points, nspins, fm_mat_Q(1), fm_mo_coeff_occ, fm_mo_coeff_virt, &
                            fm_matrix_Minv_L_kpoints, fm_matrix_L_kpoints, mat_P_global, &
                            t_3c_O, matrix_s, kpoints, eps_filter_im_time, &
                            cut_memory, nkp, num_cells_dm, num_3c_repl, &
                            size_P, ikp_local, &
                            index_to_cell_3c, &
                            cell_to_index_3c, &
                            col_blk_size, &
                            do_ic_model, do_kpoints_cubic_RPA, &
                            do_kpoints_from_Gamma, do_ri_Sigma_x, my_open_shell, &
                            has_mat_P_blocks, wkp_W, &
                            cfm_mat_Q, fm_mat_Minv_L_kpoints, fm_mat_L_kpoints, &
                            fm_mat_RI_global_work, fm_mat_work, fm_mo_coeff_occ_scaled, &
                            fm_mo_coeff_virt_scaled, mat_dm, mat_L, mat_M_P_munu_occ, mat_M_P_munu_virt, &
                            mat_MinvVMinv, mat_P_omega, mat_P_omega_kp, mat_work, mo_coeff, &
                            fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, homo, nmo)

         IF (calc_forces) CALL init_im_time_forces(force_data, fm_matrix_PQ, t_3c_M, unit_nr, mp2_env, qs_env)

         IF (my_do_gw) THEN

            CALL dbcsr_get_info(mat_P_global%matrix, &
                                row_blk_size=RI_blk_sizes)

            CALL dbcsr_get_info(matrix_s(1)%matrix, &
                                row_blk_size=prim_blk_sizes)

            gw_corr_lev_tot = gw_corr_lev_occ(1) + gw_corr_lev_virt(1)

            IF (.NOT. do_kpoints_cubic_RPA) THEN
               CALL allocate_matrices_gw_im_time(gw_corr_lev_occ, gw_corr_lev_virt, homo, nmo, &
                                                 num_integ_points, unit_nr, &
                                                 RI_blk_sizes, do_ic_model, &
                                                 para_env, fm_mat_W, fm_mat_Q(1), &
                                                 mo_coeff, &
                                                 t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                                 t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                                 starts_array_mc, ends_array_mc, &
                                                 t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                                 matrix_s, mat_W, t_3c_O, &
                                                 t_3c_O_compressed, t_3c_O_ind, &
                                                 qs_env)

            END IF
         END IF

      END IF

      IF (do_ic_model) THEN
         ! image charge model only implemented for cubic scaling GW
         CPASSERT(do_gw_im_time)
         CPASSERT(.NOT. do_periodic)
         IF (cut_memory .NE. 1) CPABORT("For IC, use MEMORY_CUT 1 in the LOW_SCALING section.")
      END IF

      ALLOCATE (e_fermi(nspins))
      IF (do_minimax_quad .OR. do_ri_sos_laplace_mp2) THEN
         do_print = .NOT. do_ic_model
         CALL get_minimax_grid(para_env, unit_nr, homo, Eigenval, num_integ_points, do_im_time, &
                               do_ri_sos_laplace_mp2, do_print, &
                               tau_tj, tau_wj, qs_env, do_gw_im_time, &
                               do_kpoints_cubic_RPA, e_fermi(1), tj, wj, &
                               weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, &
                               qs_env%mp2_env%ri_g0w0%regularization_minimax)

         !For sos_laplace_mp2 and low-scaling RPA, potentially need to store/retrieve the initial weights
         IF (qs_env%mp2_env%ri_rpa_im_time%keep_quad) THEN
            CALL keep_initial_quad(tj, wj, tau_tj, tau_wj, weights_cos_tf_t_to_w, &
                                   weights_cos_tf_w_to_t, do_ri_sos_laplace_mp2, do_im_time, &
                                   num_integ_points, unit_nr, qs_env)
         END IF
      ELSE
         IF (calc_forces) CPABORT("Forces with Clenshaw-Curtis grid not implemented.")
         CALL get_clenshaw_grid(para_env, para_env_RPA, unit_nr, homo, virtual, Eigenval, num_integ_points, &
                                num_integ_group, color_rpa_group, fm_mat_S, my_do_gw, &
                                ext_scaling, a_scaling, tj, wj)
      END IF

      ! This array is needed for RPA
      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         ALLOCATE (trace_Qomega(dimen_RI_red))
      END IF

      IF (do_ri_sos_laplace_mp2 .AND. .NOT. do_im_time) THEN
         alpha = 1.0_dp
      ELSE IF (my_open_shell .OR. do_ri_sos_laplace_mp2) THEN
         alpha = 2.0_dp
      ELSE
         alpha = 4.0_dp
      END IF

      IF (my_do_gw) THEN
         CALL allocate_matrices_gw(vec_Sigma_c_gw, color_rpa_group, dimen_nm_gw, &
                                   gw_corr_lev_occ, gw_corr_lev_virt, homo, &
                                   nmo, num_integ_group, num_integ_points, unit_nr, &
                                   gw_corr_lev_tot, num_fit_points, omega_max_fit, &
                                   do_minimax_quad, do_periodic, do_ri_Sigma_x,.NOT. do_im_time, &
                                   first_cycle_periodic_correction, &
                                   a_scaling, Eigenval, tj, vec_omega_fit_gw, vec_Sigma_x_gw, &
                                   delta_corr, Eigenval_last, Eigenval_scf, vec_W_gw, &
                                   fm_mat_S_gw, fm_mat_S_gw_work, &
                                   para_env, mp2_env, kpoints, nkp, nkp_self_energy, &
                                   do_kpoints_cubic_RPA, do_kpoints_from_Gamma)

         IF (do_bse) THEN

            num_Z_vectors = mp2_env%ri_g0w0%num_z_vectors
            eps_min_trans = mp2_env%ri_g0w0%eps_min_trans
            max_iter_bse = mp2_env%ri_g0w0%max_iter_bse

            CALL cp_fm_create(fm_mat_Q_static_bse_gemm, fm_mat_Q_gemm(1)%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q_gemm(1), fm_mat_Q_static_bse_gemm)
            CALL cp_fm_set_all(fm_mat_Q_static_bse_gemm, 0.0_dp)

            CALL cp_fm_create(fm_mat_Q_static_bse, fm_mat_Q(1)%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q_gemm(1), fm_mat_Q_static_bse)
            CALL cp_fm_set_all(fm_mat_Q_static_bse, 0.0_dp)

         END IF

      END IF

      IF (calc_forces .AND. .NOT. do_im_time) CALL rpa_grad_create(rpa_grad, fm_mat_Q(1), &
                                                                   fm_mat_S, homo, virtual, mp2_env, Eigenval(:, 1, :), &
                                                                   unit_nr, do_ri_sos_laplace_mp2)

      Erpa = 0.0_dp
      IF (mp2_env%ri_rpa%do_ri_axk) e_axk = 0.0_dp
      first_cycle = .TRUE.
      omega_old = 0.0_dp
      CALL dgemm_counter_init(dgemm_counter, unit_nr, mp2_env%ri_rpa%print_dgemm_info)

      DO count_ev_sc_GW = 1, iter_evGW

         dbcsr_time = 0.0_dp
         dbcsr_nflop = 0

         IF (do_ic_model) CYCLE

         ! reset some values, important when doing eigenvalue self-consistent GW
         IF (my_do_gw) THEN
            Erpa = 0.0_dp
            vec_Sigma_c_gw = z_zero
            first_cycle = .TRUE.
         END IF

         ! calculate Q_PQ(it)
         IF (do_im_time) THEN

            IF (.NOT. do_kpoints_cubic_RPA) THEN
               DO ispin = 1, nspins
                  e_fermi(ispin) = (Eigenval(homo(ispin), 1, ispin) + Eigenval(homo(ispin) + 1, 1, ispin))*0.5_dp
               END DO
            END IF

            tau = 0.0_dp
            tau_old = 0.0_dp

            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(/T3,A,T66,i15)") &
               "MEMORY_INFO| Memory cut:", cut_memory
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
               "SPARSITY_INFO| Eps filter for M virt/occ tensors:", eps_filter
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
               "SPARSITY_INFO| Eps filter for P matrix:", eps_filter_im_time
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,i15)") &
               "SPARSITY_INFO| Minimum tensor block size:", min_bsize

            ! for evGW, we have to ensure that mat_P_omega is zero
            CALL zero_mat_P_omega(mat_P_omega(:, :, 1))

            ! compute the matrix Q(it) and Fourier transform it directly to mat_P_omega(iw)
            CALL compute_mat_P_omega(mat_P_omega(:, :, 1), fm_scaled_dm_occ_tau, &
                                     fm_scaled_dm_virt_tau, fm_mo_coeff_occ(1), fm_mo_coeff_virt(1), &
                                     fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                     mat_P_global, matrix_s, 1, &
                                     t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                     starts_array_mc, ends_array_mc, &
                                     starts_array_mc_block, ends_array_mc_block, &
                                     weights_cos_tf_t_to_w, tj, tau_tj, e_fermi(1), eps_filter, alpha, &
                                     eps_filter_im_time, Eigenval(:, 1, 1), nmo, &
                                     num_integ_points, cut_memory, &
                                     unit_nr, mp2_env, para_env, &
                                     qs_env, do_kpoints_from_Gamma, &
                                     index_to_cell_3c, cell_to_index_3c, &
                                     has_mat_P_blocks, do_ri_sos_laplace_mp2, &
                                     dbcsr_time, dbcsr_nflop)

            ! the same for open shell, use fm_mo_coeff_occ_beta and fm_mo_coeff_virt_beta
            IF (my_open_shell) THEN
               CALL zero_mat_P_omega(mat_P_omega(:, :, 2))
               CALL compute_mat_P_omega(mat_P_omega(:, :, 2), fm_scaled_dm_occ_tau, &
                                        fm_scaled_dm_virt_tau, fm_mo_coeff_occ(2), &
                                        fm_mo_coeff_virt(2), &
                                        fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                        mat_P_global, matrix_s, 2, &
                                        t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                        starts_array_mc, ends_array_mc, &
                                        starts_array_mc_block, ends_array_mc_block, &
                                        weights_cos_tf_t_to_w, tj, tau_tj, e_fermi(2), eps_filter, alpha, &
                                        eps_filter_im_time, Eigenval(:, 1, 2), nmo, &
                                        num_integ_points, cut_memory, &
                                        unit_nr, mp2_env, para_env, &
                                        qs_env, do_kpoints_from_Gamma, &
                                        index_to_cell_3c, cell_to_index_3c, &
                                        has_mat_P_blocks, do_ri_sos_laplace_mp2, &
                                        dbcsr_time, dbcsr_nflop)

               !For RPA, we sum up the P matrices. If no force needed, can clean-up the beta spin one
               IF (.NOT. do_ri_sos_laplace_mp2) THEN
                  DO j = 1, SIZE(mat_P_omega, 2)
                     DO i = 1, SIZE(mat_P_omega, 1)
                        CALL dbcsr_add(mat_P_omega(i, j, 1)%matrix, mat_P_omega(i, j, 2)%matrix, 1.0_dp, 1.0_dp)
                        IF (.NOT. calc_forces) CALL dbcsr_clear(mat_P_omega(i, j, 2)%matrix)
                     END DO
                  END DO
               END IF
            END IF ! my_open_shell

         END IF ! do im time

         DO jquad = 1, num_integ_points

            IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE

            CALL timeset(routineN//"_RPA_matrix_operations", handle3)

            IF (do_ri_sos_laplace_mp2) THEN
               omega = tau_tj(jquad)
            ELSE
               IF (do_minimax_quad) THEN
                  omega = tj(jquad)
               ELSE
                  omega = a_scaling/TAN(tj(jquad))
               END IF
            END IF ! do_ri_sos_laplace_mp2

            IF (do_im_time) THEN
               ! in case we do imag time, we already calculated fm_mat_Q by a Fourier transform from im. time

               IF (.NOT. (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma)) THEN

                  DO ispin = 1, SIZE(mat_P_omega, 3)
                     CALL contract_P_omega_with_mat_L(mat_P_omega(jquad, 1, ispin)%matrix, mat_L%matrix, mat_work, &
                                                      eps_filter_im_time, fm_mat_work, dimen_RI, dimen_RI_red, &
                                                      fm_mat_Minv_L_kpoints(1, 1), fm_mat_Q(ispin))
                  END DO
               END IF

            ELSE

               IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3, A, 1X, I3, 1X, A, 1X, I3)") &
                  "INTEG_INFO| Started with Integration point", jquad, "of", num_integ_points

               IF (first_cycle .AND. count_ev_sc_gw > 1) THEN
                  IF (iter_sc_gw0 == 1) THEN
                     DO ispin = 1, nspins
                        CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virtual(ispin), &
                                                       Eigenval_last(:, 1, ispin), homo(ispin), omega_old)
                     END DO
                  ELSE
                     DO ispin = 1, nspins
                        CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virtual(ispin), &
                                                       Eigenval_scf(:, 1, ispin), homo(ispin), omega_old)
                     END DO
                  END IF
               END IF

               CALL calc_mat_Q(fm_mat_S(1), do_ri_sos_laplace_mp2, first_cycle, iter_sc_GW0, virtual(1), &
                               Eigenval(:, 1, 1), Eigenval_scf(:, 1, 1), &
                               homo(1), omega, omega_old, jquad, mm_style, dimen_RI_red, dimen_ia(1), alpha, fm_mat_Q(1), &
                               fm_mat_Q_gemm(1), do_bse, fm_mat_Q_static_bse_gemm, dgemm_counter)

               IF (my_open_shell) THEN
                  CALL calc_mat_Q(fm_mat_S(2), do_ri_sos_laplace_mp2, first_cycle, iter_sc_GW0, virtual(2), &
                                  Eigenval(:, 1, 2), Eigenval_scf(:, 1, 2), &
                                  homo(2), omega, omega_old, jquad, mm_style, &
                                  dimen_RI_red, dimen_ia(2), alpha, fm_mat_Q(2), fm_mat_Q_gemm(2), do_bse, &
                                  fm_mat_Q_static_bse_gemm, dgemm_counter)

                  ! For SOS-MP2 we need both matrices separately
                  IF (.NOT. do_ri_sos_laplace_mp2) THEN
                     CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_Q(1), beta=1.0_dp, matrix_b=fm_mat_Q(2))
                  END IF

               END IF ! open shell

            END IF ! im time

            ! Calculate AXK energy correction
            IF (mp2_env%ri_rpa%do_ri_axk) THEN
               CALL compute_axk_ener(qs_env, fm_mat_Q(1), fm_mat_Q_gemm(1), dimen_RI_red, dimen_ia(1), para_env_sub, &
                                     para_env_RPA, Eigenval(:, 1, 1), fm_mat_S(1), homo(1), virtual(1), omega, &
                                     mp2_env, mat_munu, unit_nr, e_axk_corr)

               ! Evaluate the final AXK energy correction
               e_axk = e_axk + e_axk_corr*wj(jquad)
            END IF ! do_ri_axk

            IF (do_ri_sos_laplace_mp2) THEN
               CALL SOS_MP2_postprocessing(fm_mat_Q, Erpa, tau_wj(jquad))

               IF (calc_forces .AND. .NOT. do_im_time) CALL rpa_grad_matrix_operations(mp2_env, rpa_grad, do_ri_sos_laplace_mp2, &
                                                           fm_mat_Q, fm_mat_Q_gemm, dgemm_counter, fm_mat_S, omega, homo, virtual, &
                                                                                       Eigenval(:, 1, :), tau_wj(jquad), unit_nr)
            ELSE
               IF (calc_forces .AND. .NOT. do_im_time) CALL rpa_grad_copy_Q(fm_mat_Q(1), rpa_grad)

               CALL Q_trace_and_add_unit_matrix(dimen_RI_red, trace_Qomega, fm_mat_Q(1), para_env_RPA)

               IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
                  CALL invert_eps_compute_W_and_Erpa_kp(dimen_RI, num_integ_points, jquad, nkp, count_ev_sc_GW, para_env, &
                                                        Erpa, tau_tj, tj, wj, weights_cos_tf_w_to_t, &
                                                        wkp_W, do_gw_im_time, do_ri_Sigma_x, do_kpoints_from_Gamma, &
                                                        cfm_mat_Q, ikp_local, &
                                                        mat_P_omega(:, :, 1), mat_P_omega_kp, qs_env, eps_filter_im_time, unit_nr, &
                                                        kpoints, fm_mat_Minv_L_kpoints, fm_mat_L_kpoints, &
                                                        fm_mat_W, fm_mat_RI_global_work, mat_MinvVMinv, &
                                                        fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv)
               ELSE
                  CALL compute_Erpa_by_freq_int(dimen_RI_red, trace_Qomega, fm_mat_Q(1), para_env_RPA, Erpa, wj(jquad))
               END IF

               IF (calc_forces .AND. .NOT. do_im_time) CALL rpa_grad_matrix_operations(mp2_env, rpa_grad, do_ri_sos_laplace_mp2, &
                                                           fm_mat_Q, fm_mat_Q_gemm, dgemm_counter, fm_mat_S, omega, homo, virtual, &
                                                                                       Eigenval(:, 1, :), wj(jquad), unit_nr)
            END IF ! do_ri_sos_laplace_mp2

            ! save omega and reset the first_cycle flag
            first_cycle = .FALSE.
            omega_old = omega

            CALL timestop(handle3)

            IF (my_do_gw) THEN

               CALL get_fermi_level_offset(fermi_level_offset, fermi_level_offset_input, Eigenval(:, 1, :), homo)

               ! do_im_time = TRUE means low-scaling calculation
               IF (do_im_time) THEN
                  ! only for molecules
                  IF (.NOT. (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma)) THEN
                    CALL compute_W_cubic_GW(fm_mat_W, fm_mat_Q(1), fm_mat_work, dimen_RI, fm_mat_Minv_L_kpoints, num_integ_points, &
                                             tj, tau_tj, weights_cos_tf_w_to_t, jquad, omega)
                  END IF
               ELSE
                  CALL compute_GW_self_energy(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI_red, gw_corr_lev_occ, &
                                              gw_corr_lev_virt, homo, jquad, nmo, num_fit_points, &
                                              do_bse, do_im_time, do_periodic, first_cycle_periodic_correction, &
                                              fermi_level_offset, &
                                              omega, Eigenval(:, 1, :), delta_corr, vec_omega_fit_gw, vec_W_gw, wj, &
                                              fm_mat_Q(1), fm_mat_Q_static_bse, fm_mat_R_gw, fm_mat_S_gw, &
                                              fm_mat_S_gw_work, mo_coeff(1), para_env, &
                                              para_env_RPA, matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, &
                                              kpoints, qs_env, mp2_env)
               END IF
            END IF

            IF (unit_nr > 0) CALL m_flush(unit_nr)
            CALL para_env%sync() ! sync to see output

         END DO ! jquad

         CALL para_env%sum(Erpa)

         IF (.NOT. do_ri_sos_laplace_mp2) THEN
            Erpa = Erpa/(pi*2.0_dp)
            IF (do_minimax_quad) Erpa = Erpa/2.0_dp
         END IF

         IF (mp2_env%ri_rpa%do_ri_axk) THEN
            CALL para_env%sum(E_axk)
            E_axk = E_axk/(pi*2.0_dp)
            IF (do_minimax_quad) E_axk = E_axk/2.0_dp
            mp2_env%ri_rpa%ener_axk = E_axk
         END IF

         IF (calc_forces .AND. do_ri_sos_laplace_mp2 .AND. do_im_time) THEN
            IF (my_open_shell) THEN
               Pspin = 1
               Qspin = 2
               CALL calc_laplace_loop_forces(force_data, mat_P_omega(:, 1, :), t_3c_M, t_3c_O(1, 1), &
                                             t_3c_O_compressed(1, 1, :), t_3c_O_ind(1, 1, :), fm_scaled_dm_occ_tau, &
                                             fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                             fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                             starts_array_mc, ends_array_mc, starts_array_mc_block, &
                                             ends_array_mc_block, num_integ_points, nmo, Eigenval(:, 1, :), &
                                             tau_tj, tau_wj, cut_memory, Pspin, Qspin, my_open_shell, &
                                             unit_nr, dbcsr_time, dbcsr_nflop, mp2_env, qs_env)
               Pspin = 2
               Qspin = 1
               CALL calc_laplace_loop_forces(force_data, mat_P_omega(:, 1, :), t_3c_M, t_3c_O(1, 1), &
                                             t_3c_O_compressed(1, 1, :), t_3c_O_ind(1, 1, :), fm_scaled_dm_occ_tau, &
                                             fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                             fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                             starts_array_mc, ends_array_mc, starts_array_mc_block, &
                                             ends_array_mc_block, num_integ_points, nmo, Eigenval(:, 1, :), &
                                             tau_tj, tau_wj, cut_memory, Pspin, Qspin, my_open_shell, &
                                             unit_nr, dbcsr_time, dbcsr_nflop, mp2_env, qs_env)

            ELSE
               Pspin = 1
               Qspin = 1
               CALL calc_laplace_loop_forces(force_data, mat_P_omega(:, 1, :), t_3c_M, t_3c_O(1, 1), &
                                             t_3c_O_compressed(1, 1, :), t_3c_O_ind(1, 1, :), fm_scaled_dm_occ_tau, &
                                             fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                             fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                             starts_array_mc, ends_array_mc, starts_array_mc_block, &
                                             ends_array_mc_block, num_integ_points, nmo, Eigenval(:, 1, :), &
                                             tau_tj, tau_wj, cut_memory, Pspin, Qspin, my_open_shell, &
                                             unit_nr, dbcsr_time, dbcsr_nflop, mp2_env, qs_env)
            END IF
            CALL calc_post_loop_forces(force_data, unit_nr, qs_env)
         END IF !laplace SOS-MP2

         IF (calc_forces .AND. do_im_time .AND. .NOT. do_ri_sos_laplace_mp2) THEN
            DO ispin = 1, nspins
               CALL calc_rpa_loop_forces(force_data, mat_P_omega(:, 1, :), t_3c_M, t_3c_O(1, 1), &
                                         t_3c_O_compressed(1, 1, :), t_3c_O_ind(1, 1, :), fm_scaled_dm_occ_tau, &
                                         fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                         fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                         starts_array_mc, ends_array_mc, starts_array_mc_block, &
                                         ends_array_mc_block, num_integ_points, nmo, Eigenval(:, 1, :), &
                                         e_fermi(ispin), weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, tj, &
                                         wj, tau_tj, cut_memory, ispin, my_open_shell, unit_nr, dbcsr_time, &
                                         dbcsr_nflop, mp2_env, qs_env)
            END DO
            CALL calc_post_loop_forces(force_data, unit_nr, qs_env)
         END IF

         IF (do_im_time) THEN

            my_flop_rate = REAL(dbcsr_nflop, dp)/(1.0E09_dp*dbcsr_time)
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(/T3,A,T73,ES8.2)") &
               "PERFORMANCE| DBCSR total number of flops:", REAL(dbcsr_nflop*para_env%num_pe, dp)
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.2)") &
               "PERFORMANCE| DBCSR total execution time:", dbcsr_time
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.2)") &
               "PERFORMANCE| DBCSR flop rate (Gflops / MPI rank):", my_flop_rate

         ELSE

            CALL dgemm_counter_write(dgemm_counter, para_env)

         END IF

         ! GW: for low-scaling calculation: Compute self-energy Sigma(i*tau), Sigma(i*omega)
         ! for low-scaling and ordinary-scaling: analytic continuation from Sigma(iw) -> Sigma(w)
         !                                       and correction of quasiparticle energies e_n^GW
         IF (my_do_gw) THEN

            CALL compute_QP_energies(vec_Sigma_c_gw, count_ev_sc_GW, gw_corr_lev_occ, &
                                     gw_corr_lev_tot, gw_corr_lev_virt, homo, &
                                     nmo, num_fit_points, num_integ_points, &
                                     unit_nr, do_apply_ic_corr_to_gw, do_im_time, &
                                     do_periodic, do_ri_Sigma_x, first_cycle_periodic_correction, &
                                     e_fermi, eps_filter, fermi_level_offset, &
                                     delta_corr, Eigenval, &
                                     Eigenval_last, Eigenval_scf, iter_sc_GW0, exit_ev_gw, tau_tj, tj, &
                                     vec_omega_fit_gw, vec_Sigma_x_gw, mp2_env%ri_g0w0%ic_corr_list, &
                                     weights_cos_tf_t_to_w, weights_sin_tf_t_to_w, &
                                     fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, fm_mo_coeff_occ, &
                                     fm_mo_coeff_virt, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                                     mo_coeff(1), fm_mat_W, para_env, para_env_RPA, mat_dm, mat_MinvVMinv, &
                                     t_3c_O, t_3c_M, t_3c_overl_int_ao_mo, t_3c_O_compressed, t_3c_O_mo_compressed, &
                                     t_3c_O_ind, t_3c_O_mo_ind, &
                                     t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                     matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, mat_W, matrix_s, &
                                     kpoints, mp2_env, qs_env, nkp_self_energy, do_kpoints_cubic_RPA, &
                                     starts_array_mc, ends_array_mc)

            ! if HOMO-LUMO gap differs by less than mp2_env%ri_g0w0%eps_ev_sc_iter, exit ev sc GW loop
            IF (exit_ev_gw) EXIT

         END IF ! my_do_gw if

      END DO ! evGW loop

      IF (do_ic_model) THEN

         IF (my_open_shell) THEN

            CALL calculate_ic_correction(Eigenval(:, 1, 1), mat_MinvVMinv%matrix, &
                                         t_3c_overl_nnP_ic(1), t_3c_overl_nnP_ic_reflected(1), &
                                         gw_corr_lev_tot, &
                                         gw_corr_lev_occ(1), gw_corr_lev_virt(1), homo(1), unit_nr, &
                                         print_ic_values, para_env, do_alpha=.TRUE.)

            CALL calculate_ic_correction(Eigenval(:, 1, 2), mat_MinvVMinv%matrix, &
                                         t_3c_overl_nnP_ic(2), t_3c_overl_nnP_ic_reflected(2), &
                                         gw_corr_lev_tot, &
                                         gw_corr_lev_occ(2), gw_corr_lev_virt(2), homo(2), unit_nr, &
                                         print_ic_values, para_env, do_beta=.TRUE.)

         ELSE

            CALL calculate_ic_correction(Eigenval(:, 1, 1), mat_MinvVMinv%matrix, &
                                         t_3c_overl_nnP_ic(1), t_3c_overl_nnP_ic_reflected(1), &
                                         gw_corr_lev_tot, &
                                         gw_corr_lev_occ(1), gw_corr_lev_virt(1), homo(1), unit_nr, &
                                         print_ic_values, para_env)

         END IF

      END IF

      ! postprocessing after GW for Bethe-Salpeter
      IF (do_bse) THEN
         CALL mult_B_with_W_and_fill_local_3c_arrays(fm_mat_S_ij_bse, fm_mat_S_ab_bse, fm_mat_S(1), fm_mat_Q_static_bse, &
                                                     fm_mat_Q_static_bse_gemm, &
                                                     B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                     B_iaQ_bse_local, dimen_RI_red, homo(1), virtual(1), dimen_ia(1), &
                                                     gd_array, color_sub, para_env)

         CALL do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo(1), virtual(1), num_Z_vectors, &
                                     max_iter_bse, eps_min_trans, Eigenval(:, 1, 1), para_env)

      END IF

      IF (my_do_gw) THEN
         CALL deallocate_matrices_gw(fm_mat_S_gw_work, vec_W_gw, vec_Sigma_c_gw, vec_omega_fit_gw, &
                                     mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw, &
                                     Eigenval_last, Eigenval_scf, do_periodic, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                     kpoints, vec_Sigma_x_gw,.NOT. do_im_time)
      END IF

      IF (do_im_time) THEN

         CALL dealloc_im_time(fm_mo_coeff_occ, fm_mo_coeff_virt, &
                              fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, index_to_cell_3c, &
                              cell_to_index_3c, do_ic_model, &
                              do_kpoints_cubic_RPA, do_kpoints_from_Gamma, do_ri_Sigma_x, &
                              has_mat_P_blocks, &
                              wkp_W, cfm_mat_Q, fm_mat_Minv_L_kpoints, fm_mat_L_kpoints, &
                              fm_matrix_Minv, fm_matrix_Minv_Vtrunc_Minv, fm_mat_RI_global_work, fm_mat_work, &
                              fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, mat_dm, mat_L, &
                              mat_MinvVMinv, mat_P_omega, mat_P_omega_kp, &
                              t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, mat_work, qs_env)

         IF (my_do_gw) THEN
            CALL deallocate_matrices_gw_im_time(weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, do_ic_model, &
                                                do_kpoints_cubic_RPA, fm_mat_W, &
                                                t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                                t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                                t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                                mat_W, qs_env)
         END IF

      END IF

      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         DEALLOCATE (tj)
         DEALLOCATE (wj)
         DEALLOCATE (trace_Qomega)
      END IF

      IF (do_im_time .OR. do_ri_sos_laplace_mp2) THEN
         DEALLOCATE (tau_tj)
         DEALLOCATE (tau_wj)
      END IF

      IF (do_im_time .AND. calc_forces) THEN
         CALL im_time_force_release(force_data)
      END IF

      IF (calc_forces .AND. .NOT. do_im_time) CALL rpa_grad_finalize(rpa_grad, mp2_env, para_env_sub, para_env, &
                                                                     qs_env, gd_array, color_sub, do_ri_sos_laplace_mp2, &
                                                                     homo, virtual)

      CALL timestop(handle)

   END SUBROUTINE rpa_num_int

END MODULE rpa_main
