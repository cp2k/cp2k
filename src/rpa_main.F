!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_cfm_types,                    ONLY: cp_cfm_p_type,&
                                              cp_cfm_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_split
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dbcsr_tensor_api,                ONLY: dbcsr_t_type
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
   USE message_passing,                 ONLY: mp_min,&
                                              mp_sendrecv,&
                                              mp_sum
   USE minimax_exp,                     ONLY: check_exp_minimax_range
   USE mp2_grids,                       ONLY: get_clenshaw_grid,&
                                              get_minimax_grid
   USE mp2_laplace,                     ONLY: SOS_MP2_postprocessing
   USE mp2_ri_grad_util,                ONLY: array2fm
   USE mp2_types,                       ONLY: dgemm_counter_init,&
                                              dgemm_counter_type,&
                                              dgemm_counter_write,&
                                              mp2_type,&
                                              three_dim_real_array,&
                                              two_dim_int_array,&
                                              two_dim_real_array
   USE qs_band_structure,               ONLY: calculate_kp_orbitals
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE rpa_axk,                         ONLY: compute_axk_ener
   USE rpa_gw,                          ONLY: GW_matrix_operations,&
                                              allocate_matrices_gw,&
                                              allocate_matrices_gw_im_time,&
                                              compute_QP_energies,&
                                              deallocate_matrices_gw,&
                                              deallocate_matrices_gw_im_time
   USE rpa_gw_ic,                       ONLY: calculate_ic_correction
   USE rpa_im_time,                     ONLY: compute_mat_P_omega,&
                                              zero_mat_P_omega
   USE rpa_kpoints,                     ONLY: RPA_postprocessing_kp
   USE rpa_util,                        ONLY: RPA_postprocessing_nokp,&
                                              RPA_postprocessing_start,&
                                              alloc_im_time,&
                                              calc_mat_Q,&
                                              contract_P_omega_with_mat_L,&
                                              dealloc_im_time
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
!> \param fm_matrix_L_RI_metric ...
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
! **************************************************************************************************
   SUBROUTINE rpa_ri_compute_en(qs_env, Erpa, mp2_env, BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, &
                                para_env, para_env_sub, color_sub, &
                                gd_array, gd_B_virtual, gd_B_all, gd_B_occ_bse, gd_B_virt_bse, &
                                mo_coeff, fm_matrix_L_RI_metric, kpoints, &
                                Eigenval, nmo, homo, dimen_RI, dimen_RI_red, gw_corr_lev_occ, gw_corr_lev_virt, &
                                unit_nr, do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_bse, matrix_s, &
                                mat_munu, mat_P_global, &
                                t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                starts_array_mc, ends_array_mc, &
                                starts_array_mc_block, ends_array_mc_block)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: Erpa
      TYPE(mp2_type), INTENT(INOUT), POINTER             :: mp2_env
      TYPE(three_dim_real_array), DIMENSION(:), &
         INTENT(INOUT)                                   :: BIb_C, BIb_C_gw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: BIb_C_bse_ij, BIb_C_bse_ab
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(INOUT)                             :: color_sub
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_array
      TYPE(group_dist_d1_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: gd_B_virtual
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_B_all, gd_B_occ_bse, gd_B_virt_bse
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Eigenval
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
      TYPE(dbcsr_t_type)                                 :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_O
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_compressed
      TYPE(block_ind_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(INOUT)               :: t_3c_O_ind
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_ri_compute_en'

      INTEGER :: best_integ_group_size, best_num_integ_point, color_rpa_group, dimen_homo_square, &
         dimen_nm_gw, dimen_virt_square, handle, handle2, handle3, i, ierr, iiB, &
         input_integ_group_size, integ_group_size, ispin, jjB, min_integ_group_size, &
         my_ab_comb_bse_end, my_ab_comb_bse_size, my_ab_comb_bse_start, my_group_L_end, &
         my_group_L_size, my_group_L_start, my_ij_comb_bse_end, my_ij_comb_bse_size, &
         my_ij_comb_bse_start, my_nm_gw_end, my_nm_gw_size, my_nm_gw_start, ncol_block_mat, &
         ngroup, nrow_block_mat, nspins, num_integ_group, num_integ_points, pos_integ_group
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dimen_ia, my_ia_end, my_ia_size, &
                                                            my_ia_start, sub_proc_map, virtual
      LOGICAL                                            :: do_kpoints_from_Gamma, do_minimax_quad, &
                                                            my_open_shell, skip_integ_group_opt
      REAL(KIND=dp) :: allowed_memory, avail_mem, E_Range, Emax, Emin, mem_for_iaK, mem_for_QK, &
         mem_min, mem_per_group, mem_real, needed_mem
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D_bse_ab, BIb_C_2D_bse_ij
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: fm_mat_Q, fm_mat_Q_gemm, fm_mat_S, &
                                                            fm_mat_S_gw, fm_mo_coeff_occ, &
                                                            fm_mo_coeff_virt
      TYPE(cp_fm_type), POINTER                          :: fm_mat_R_gw, fm_mat_S_ab_bse, &
                                                            fm_mat_S_ij_bse, fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(two_dim_real_array), ALLOCATABLE, &
         DIMENSION(:)                                    :: BIb_C_2D, BIb_C_2D_gw

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

      IF (do_im_time) mp2_env%ri_rpa%minimax_quad = .TRUE.
      do_minimax_quad = mp2_env%ri_rpa%minimax_quad

      IF (do_ri_sos_laplace_mp2) THEN
         num_integ_points = mp2_env%ri_laplace%n_quadrature
         input_integ_group_size = mp2_env%ri_laplace%integ_group_size

         ! check the range for the minimax approximation
         Emin = HUGE(dp)
         Emax = 0.0_dp
         DO ispin = 1, nspins
            IF (homo(ispin) > 0) THEN
               Emin = MIN(Emin, 2.0_dp*(Eigenval(homo(ispin) + 1, ispin) - Eigenval(homo(ispin), ispin)))
               Emax = MAX(Emax, 2.0_dp*(MAXVAL(Eigenval(:, ispin)) - MINVAL(Eigenval(:, ispin))))
            END IF
         END DO
         E_Range = Emax/Emin
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
         input_integ_group_size = mp2_env%ri_rpa%rpa_integ_group_size

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
         ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
         ! has not been given back to the OS yet.
         CALL mp_min(mem_real, para_env%group)

         mem_min = 2.0_dp*REAL(homo(1), KIND=dp)*maxsize(gd_B_virtual(1))*maxsize(gd_array)*8.0_dp/(1024**2)

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Minimum required memory per MPI process:', &
            mem_min, ' MiB'

         mem_real = allowed_memory - mem_real
         mem_real = MAX(mem_real, mem_min)

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Available memory per MPI process:', &
            mem_real, ' MiB'

         mem_per_group = mem_real*para_env_sub%num_pe

         needed_mem = mem_for_iaK*2.0_dp + mem_for_QK*3.0_dp

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

            integ_group_size = min_integ_group_size - 1
            DO iiB = min_integ_group_size + 1, ngroup
               integ_group_size = integ_group_size + 1

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

      ! create the sub_proc_map
      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe - 1))
      DO i = 0, para_env_sub%num_pe - 1
         sub_proc_map(i) = i
         sub_proc_map(-i - 1) = para_env_sub%num_pe - i - 1
         sub_proc_map(para_env_sub%num_pe + i) = i
      END DO

      ! not necessary for imaginary time

      ALLOCATE (BIb_C_2D(nspins))
      IF (.NOT. do_im_time) THEN

         ! reorder the local data in such a way to help the next stage of matrix creation
         ! now the data inside the group are divided into a ia x K matrix
         DO ispin = 1, nspins
            CALL calculate_BIb_C_2D(BIb_C_2D(ispin)%array, BIb_C(ispin)%array, para_env_sub, dimen_ia(ispin), &
                                    homo(ispin), virtual(ispin), gd_B_virtual(ispin), &
                                    sub_proc_map, my_ia_size(ispin), my_ia_start(ispin), my_ia_end(ispin), my_group_L_size)

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
                                       sub_proc_map, my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, my_group_L_size)
               DEALLOCATE (BIb_C_gw(ispin)%array)
            END DO

            CALL release_group_dist(gd_B_all)

            CALL timestop(handle3)

         END IF
      END IF

      IF (do_bse) THEN

         CALL timeset(routineN//"_reorder_bse1", handle3)

         dimen_homo_square = homo(1)**2

         CALL calculate_BIb_C_2D(BIb_C_2D_bse_ij, BIb_C_bse_ij, para_env_sub, dimen_homo_square, &
                                 homo(1), homo(1), gd_B_occ_bse, &
                                 sub_proc_map, my_ij_comb_bse_size, my_ij_comb_bse_start, my_ij_comb_bse_end, my_group_L_size)

         DEALLOCATE (BIb_C_bse_ij)
         CALL release_group_dist(gd_B_occ_bse)

         CALL timestop(handle3)

         CALL timeset(routineN//"_reorder_bse2", handle3)

         dimen_virt_square = virtual(1)**2

         CALL calculate_BIb_C_2D(BIb_C_2D_bse_ab, BIb_C_bse_ab, para_env_sub, dimen_virt_square, &
                                 virtual(1), virtual(1), gd_B_virt_bse, &
                                 sub_proc_map, my_ab_comb_bse_size, my_ab_comb_bse_start, my_ab_comb_bse_end, my_group_L_size)

         DEALLOCATE (BIb_C_bse_ab)
         CALL release_group_dist(gd_B_virt_bse)

         CALL timestop(handle3)

      END IF

      CALL timestop(handle2)

      ! now create the matrices needed for the calculation, Q, S and G
      ! Q and G will have omega dependence

      ALLOCATE (fm_mat_Q(nspins), fm_mat_Q_gemm(nspins), fm_mat_S(nspins), fm_mo_coeff_occ(nspins), fm_mo_coeff_virt(nspins))
      DO ispin = 1, nspins
         NULLIFY (fm_mat_Q(ispin)%matrix, fm_mat_Q_gemm(ispin)%matrix, fm_mat_S(ispin)%matrix, &
                  fm_mo_coeff_occ(ispin)%matrix, fm_mo_coeff_virt(ispin)%matrix)
      END DO
      IF (my_open_shell .AND. do_im_time .AND. do_ri_sos_laplace_mp2) THEN
         CALL create_integ_mat(BIb_C_2D(1)%array, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, dimen_ia(1), dimen_ia(1), color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ia_size(1), my_ia_start(1), my_ia_end(1), &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S(1)%matrix, nrow_block_mat, ncol_block_mat, &
                               do_im_time=do_im_time, fm_mat_Q_gemm=fm_mat_Q_gemm(1)%matrix, fm_mat_Q=fm_mat_Q(1)%matrix, &
                               fm_scaled_dm_occ_tau=fm_scaled_dm_occ_tau, &
                               fm_scaled_dm_virt_tau=fm_scaled_dm_virt_tau, mo_coeff=mo_coeff(1)%matrix, &
                               fm_mo_coeff_occ=fm_mo_coeff_occ(1)%matrix, fm_mo_coeff_virt=fm_mo_coeff_virt(1)%matrix, &
                               nmo=nmo, homo=homo(1), fm_mat_Q_beta=fm_mat_Q(2)%matrix)
      ELSE
         CALL create_integ_mat(BIb_C_2D(1)%array, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, dimen_ia(1), dimen_ia(1), color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ia_size(1), my_ia_start(1), my_ia_end(1), &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S(1)%matrix, nrow_block_mat, ncol_block_mat, &
                               do_im_time=do_im_time, fm_mat_Q_gemm=fm_mat_Q_gemm(1)%matrix, fm_mat_Q=fm_mat_Q(1)%matrix, &
                               fm_scaled_dm_occ_tau=fm_scaled_dm_occ_tau, &
                               fm_scaled_dm_virt_tau=fm_scaled_dm_virt_tau, mo_coeff=mo_coeff(1)%matrix, &
                               fm_mo_coeff_occ=fm_mo_coeff_occ(1)%matrix, fm_mo_coeff_virt=fm_mo_coeff_virt(1)%matrix, &
                               nmo=nmo, homo=homo(1))
      END IF

      IF (my_open_shell) THEN

         ! for imaginary time, we only have to build beta mo coefficients
         IF (do_im_time) THEN

            CALL create_occ_virt_mo_coeffs(fm_mo_coeff_occ(2)%matrix, fm_mo_coeff_virt(2)%matrix, mo_coeff(2)%matrix, &
                                           nmo, homo(2))

            ! for RPA with imaginary frequency, we have to build the same matrices for beta as for alpha
         ELSE

            CALL create_integ_mat(BIb_C_2D(2)%array, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                                  dimen_RI_red, dimen_ia(2), dimen_ia(2), color_rpa_group, &
                                  mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                                  my_ia_size(2), my_ia_start(2), my_ia_end(2), &
                                  my_group_L_size, my_group_L_start, my_group_L_end, &
                                  para_env_RPA, fm_mat_S(2)%matrix, nrow_block_mat, ncol_block_mat, &
                                  .TRUE., fm_mat_Q(1)%matrix%matrix_struct%context, &
                                  fm_mat_Q_gemm=fm_mat_Q_gemm(2)%matrix, fm_mat_Q=fm_mat_Q(2)%matrix)

         END IF

      END IF
      DEALLOCATE (BIb_C_2D, my_ia_end, my_ia_size, my_ia_start)

      ! for GW, we need other matrix fm_mat_S, always allocate the container to prevent crying compilers
      ALLOCATE (fm_mat_S_gw(nspins))
      IF (my_do_gw .AND. .NOT. do_im_time) THEN

         CALL create_integ_mat(BIb_C_2D_gw(1)%array, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, dimen_nm_gw, dimen_ia(1), color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_gw(1)%matrix, nrow_block_mat, ncol_block_mat, &
                               .TRUE., fm_mat_Q(1)%matrix%matrix_struct%context, fm_mat_Q(1)%matrix%matrix_struct%context, &
                               fm_mat_Q=fm_mat_R_gw)

         IF (my_open_shell) THEN
            CALL create_integ_mat(BIb_C_2D_gw(2)%array, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                                  dimen_RI_red, dimen_nm_gw, dimen_ia(1), color_rpa_group, &
                                  mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                                  my_nm_gw_size, my_nm_gw_start, my_nm_gw_end, &
                                  my_group_L_size, my_group_L_start, my_group_L_end, &
                                  para_env_RPA, fm_mat_S_gw(2)%matrix, nrow_block_mat, ncol_block_mat, &
                                  .TRUE., fm_mat_Q(1)%matrix%matrix_struct%context, fm_mat_Q(1)%matrix%matrix_struct%context)

         END IF
         DEALLOCATE (BIb_C_2D_gw)

      END IF

      ! for Bethe-Salpeter, we need other matrix fm_mat_S
      IF (do_bse) THEN

         CALL create_integ_mat(BIb_C_2D_bse_ij, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, dimen_homo_square, dimen_ia(1), color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ij_comb_bse_size, my_ij_comb_bse_start, my_ij_comb_bse_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_ij_bse, nrow_block_mat, ncol_block_mat, &
                               .TRUE., fm_mat_Q(1)%matrix%matrix_struct%context, fm_mat_Q(1)%matrix%matrix_struct%context)

         CALL create_integ_mat(BIb_C_2D_bse_ab, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI_red, dimen_virt_square, dimen_ia(1), color_rpa_group, &
                               mp2_env%block_size_row, mp2_env%block_size_col, unit_nr, &
                               my_ab_comb_bse_size, my_ab_comb_bse_start, my_ab_comb_bse_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S_ab_bse, &
                               nrow_block_mat, ncol_block_mat, &
                               .TRUE., fm_mat_Q(1)%matrix%matrix_struct%context, fm_mat_Q(1)%matrix%matrix_struct%context)

      END IF

      do_kpoints_from_Gamma = SUM(mp2_env%ri_rpa_im_time%kp_grid) > 0
      IF (do_kpoints_from_Gamma) THEN
         CALL get_bandstruc_and_k_dependent_MOs(qs_env, kpoints, virtual(1))
      END IF

      ! Now start the RPA calculation
      ! fm_mo_coeff_occ, fm_mo_coeff_virt will be deallocated here
      CALL rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                       homo, virtual, dimen_RI, dimen_RI_red, dimen_ia, dimen_nm_gw, &
                       Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                       fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                       fm_mat_S_ij_bse, fm_mat_S_ij_bse, &
                       my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                       do_minimax_quad, &
                       do_im_time, fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                       fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                       mat_munu, mat_P_global, &
                       t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                       starts_array_mc, ends_array_mc, &
                       starts_array_mc_block, ends_array_mc_block, &
                       matrix_s, do_kpoints_from_Gamma, kpoints, gd_array, color_sub, &
                       do_ri_sos_laplace_mp2=do_ri_sos_laplace_mp2)

      DEALLOCATE (sub_proc_map)

      CALL release_group_dist(gd_array)

      CALL cp_para_env_release(para_env_RPA)

      IF (.NOT. do_im_time) THEN
         DO ispin = 1, nspins
            CALL cp_fm_release(fm_mat_Q_gemm(ispin)%matrix)
            CALL cp_fm_release(fm_mat_S(ispin)%matrix)
         END DO
         DEALLOCATE (fm_mat_Q_gemm, fm_mat_S)
      END IF
      DO ispin = 1, nspins
         IF (ASSOCIATED(fm_mat_Q(ispin)%matrix)) CALL cp_fm_release(fm_mat_Q(ispin)%matrix)
      END DO
      DEALLOCATE (fm_mat_Q)

      IF (my_do_gw .AND. .NOT. do_im_time) THEN
         DO ispin = 1, nspins
            CALL cp_fm_release(fm_mat_S_gw(ispin)%matrix)
         END DO
         DEALLOCATE (fm_mat_S_gw)
         CALL cp_fm_release(fm_mat_R_gw)
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
      END IF

      CALL timestop(handle)

   END SUBROUTINE rpa_ri_compute_en

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE get_bandstruc_and_k_dependent_MOs(qs_env, kpoints, virtual)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_bandstruc_and_k_dependent_MOs'

      INTEGER                                            :: handle, i_dim
      INTEGER, DIMENSION(3)                              :: nkp_grid_G, nkp_grid_W
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (qs_env%mp2_env%ri_rpa_im_time%kpoints_G, &
               qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma, &
               qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma_no_xc, &
               para_env)

      nkp_grid_W(1:3) = kpoints%nkp_grid(1:3)
      DO i_dim = 1, 3
         IF (nkp_grid_W(i_dim) == 1) THEN
            nkp_grid_G(i_dim) = 1
         ELSE
            nkp_grid_G(i_dim) = nkp_grid_W(i_dim)*qs_env%mp2_env%ri_rpa_im_time%k_mesh_g_factor
         END IF
      END DO

      CALL get_qs_env(qs_env=qs_env, para_env=para_env)

      CALL create_kp_and_calc_kp_orbitals(qs_env, qs_env%mp2_env%ri_rpa_im_time%kpoints_G, &
                                          "MONKHORST-PACK", virtual, para_env%num_pe, &
                                          mp_grid=nkp_grid_G(1:3))

      CALL timestop(handle)

   END SUBROUTINE get_bandstruc_and_k_dependent_MOs

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpoints ...
!> \param scheme ...
!> \param nadd ...
!> \param group_size_ext ...
!> \param mp_grid ...
!> \param kpgeneral ...
!> \param with_xc_terms ...
! **************************************************************************************************
   SUBROUTINE create_kp_and_calc_kp_orbitals(qs_env, kpoints, scheme, nadd, &
                                             group_size_ext, mp_grid, kpgeneral, with_xc_terms)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*), INTENT(IN)                       :: scheme
      INTEGER                                            :: nadd, group_size_ext
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: mp_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: kpgeneral
      LOGICAL, OPTIONAL                                  :: with_xc_terms

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_kp_and_calc_kp_orbitals'

      INTEGER                                            :: handle
      LOGICAL                                            :: my_with_xc_terms
      TYPE(qs_environment_type), POINTER                 :: qs_env_kp_from_Gamma

      CALL timeset(routineN, handle)

      NULLIFY (qs_env_kp_from_Gamma)

      my_with_xc_terms = .TRUE.
      IF (PRESENT(with_xc_terms)) my_with_xc_terms = with_xc_terms

      CALL create_kp_from_gamma(qs_env, qs_env_kp_from_Gamma, with_xc_terms=my_with_xc_terms)

      ! k-dep. MO coeff. and band structure for Green's function
      CALL calculate_kp_orbitals(qs_env_kp_from_Gamma, kpoints, scheme, nadd=nadd, kpgeneral=kpgeneral, &
                                 mp_grid=mp_grid, group_size_ext=group_size_ext)

      CALL qs_env_release(qs_env_kp_from_Gamma)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param kpgeneral ...
! **************************************************************************************************
   SUBROUTINE get_kpgeneral_for_Sigma_kpoints(qs_env, kpgeneral)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: kpgeneral

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_kpgeneral_for_Sigma_kpoints'

      INTEGER                                            :: handle, i_kp_in_kp_line, i_special_kp, &
                                                            ikk, n_kp_in_kp_line, n_special_kp, &
                                                            nkp_bandstruc

      CALL timeset(routineN, handle)

      n_special_kp = qs_env%mp2_env%ri_g0w0%n_special_kp
      n_kp_in_kp_line = qs_env%mp2_env%ri_g0w0%n_kp_in_kp_line

      nkp_bandstruc = n_kp_in_kp_line*(n_special_kp - 1) + 1

      ALLOCATE (kpgeneral(3, nkp_bandstruc))

      kpgeneral(1:3, 1) = qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, 1)

      ikk = 1

      DO i_special_kp = 2, n_special_kp
         DO i_kp_in_kp_line = 1, n_kp_in_kp_line

            ikk = ikk + 1
            kpgeneral(1:3, ikk) = qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, i_special_kp - 1) + &
                                  REAL(i_kp_in_kp_line, KIND=dp)/REAL(n_kp_in_kp_line, KIND=dp)* &
                                  (qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, i_special_kp) - &
                                   qs_env%mp2_env%ri_g0w0%xkp_special_kp(1:3, i_special_kp - 1))

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_kpgeneral_for_Sigma_kpoints

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
!> \param sub_proc_map ...
!> \param my_ia_size ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \author Jan Wilhelm, 03/2015
! **************************************************************************************************
   SUBROUTINE calculate_BIb_C_2D(BIb_C_2D, BIb_C, para_env_sub, dimen_ia, homo, virtual, &
                                 gd_B_virtual, &
                                 sub_proc_map, my_ia_size, my_ia_start, my_ia_end, my_group_L_size)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: BIb_C_2D
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: BIb_C
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER, INTENT(IN)                                :: dimen_ia, homo, virtual
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_B_virtual
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sub_proc_map
      INTEGER                                            :: my_ia_size, my_ia_start, my_ia_end, &
                                                            my_group_L_size

      INTEGER, PARAMETER                                 :: occ_chunk = 128

      INTEGER :: ia_global, iiB, itmp(2), jjB, my_B_size, my_B_virtual_start, occ_high, occ_low, &
         proc_receive, proc_send, proc_shift, rec_B_size, rec_B_virtual_end, rec_B_virtual_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C_rec

      itmp = get_limit(dimen_ia, para_env_sub%num_pe, para_env_sub%mepos)
      my_ia_start = itmp(1)
      my_ia_end = itmp(2)
      my_ia_size = my_ia_end - my_ia_start + 1

      CALL get_group_dist(gd_B_virtual, para_env_sub%mepos, sizes=my_B_size, starts=my_B_virtual_start)

      ! reorder data
      ALLOCATE (BIb_C_2D(my_ia_size, my_group_L_size))

!$OMP     PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
!$OMP              SHARED(homo,my_B_size,virtual,my_B_virtual_start,my_ia_start,my_ia_end,BIb_C,BIb_C_2D,&
!$OMP              my_group_L_size)
      DO iiB = 1, homo
         DO jjB = 1, my_B_size
            ia_global = (iiB - 1)*virtual + my_B_virtual_start + jjB - 1
            IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
               BIb_C_2D(ia_global - my_ia_start + 1, 1:my_group_L_size) = BIb_C(1:my_group_L_size, jjB, iiB)
            END IF
         END DO
      END DO

      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

         CALL get_group_dist(gd_B_virtual, proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

         ! do this in chunks to avoid high memory overhead  for both BIb_C_rec and buffers in mp_sendrecv
         ! TODO: fix this more cleanly with a rewrite sending only needed data etc.
         ! TODO: occ_chunk should presumably be precomputed so that messages are limited to e.g. 100MiB.
         ALLOCATE (BIb_C_rec(my_group_L_size, rec_B_size, MIN(homo, occ_chunk)))

         DO occ_low = 1, homo, occ_chunk
            occ_high = MIN(homo, occ_low + occ_chunk - 1)
            CALL mp_sendrecv(BIb_C(:, :, occ_low:occ_high), proc_send, &
                             BIb_C_rec(:, :, 1:occ_high - occ_low + 1), proc_receive, &
                             para_env_sub%group)
!$OMP          PARALLEL DO DEFAULT(NONE) PRIVATE(jjB,ia_global) &
!$OMP                   SHARED(occ_low,occ_high,rec_B_size,virtual,rec_B_virtual_start,my_ia_start,my_ia_end,BIb_C_rec,BIb_C_2D,&
!$OMP                          my_group_L_size)
            DO iiB = occ_low, occ_high
               DO jjB = 1, rec_B_size
                  ia_global = (iiB - 1)*virtual + rec_B_virtual_start + jjB - 1
                  IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
                     BIb_C_2D(ia_global - my_ia_start + 1, 1:my_group_L_size) = BIb_C_rec(1:my_group_L_size, jjB, iiB - occ_low + 1)
                  END IF
               END DO
            END DO
         END DO

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
!> \param nrow_block_mat ...
!> \param ncol_block_mat ...
!> \param beta_case ...
!> \param blacs_env_ext ...
!> \param blacs_env_ext_S ...
!> \param do_im_time ...
!> \param fm_mat_Q_gemm ...
!> \param fm_mat_Q ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param mo_coeff ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param nmo ...
!> \param homo ...
!> \param fm_mat_Q_beta ...
! **************************************************************************************************
   SUBROUTINE create_integ_mat(BIb_C_2D, para_env, para_env_sub, color_sub, ngroup, integ_group_size, &
                               dimen_RI, dimen_ia, dimen_ia_for_block_size, color_rpa_group, &
                               ext_row_block_size, ext_col_block_size, unit_nr, &
                               my_ia_size, my_ia_start, my_ia_end, &
                               my_group_L_size, my_group_L_start, my_group_L_end, &
                               para_env_RPA, fm_mat_S, nrow_block_mat, ncol_block_mat, &
                               beta_case, blacs_env_ext, blacs_env_ext_S, &
                               do_im_time, fm_mat_Q_gemm, fm_mat_Q, fm_scaled_dm_occ_tau, &
                               fm_scaled_dm_virt_tau, mo_coeff, fm_mo_coeff_occ, &
                               fm_mo_coeff_virt, nmo, homo, &
                               fm_mat_Q_beta)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: BIb_C_2D
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN) :: color_sub, ngroup, integ_group_size, dimen_RI, dimen_ia, &
         dimen_ia_for_block_size, color_rpa_group, ext_row_block_size, ext_col_block_size, &
         unit_nr, my_ia_size, my_ia_start, my_ia_end, my_group_L_size, my_group_L_start, &
         my_group_L_end
      TYPE(cp_para_env_type), POINTER                    :: para_env_RPA
      TYPE(cp_fm_type), POINTER                          :: fm_mat_S
      INTEGER, INTENT(INOUT)                             :: nrow_block_mat, ncol_block_mat
      LOGICAL, INTENT(IN), OPTIONAL                      :: beta_case
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env_ext, blacs_env_ext_S
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_im_time
      TYPE(cp_fm_type), OPTIONAL, POINTER :: fm_mat_Q_gemm, fm_mat_Q, fm_scaled_dm_occ_tau, &
         fm_scaled_dm_virt_tau, mo_coeff, fm_mo_coeff_occ, fm_mo_coeff_virt
      INTEGER, INTENT(IN), OPTIONAL                      :: nmo, homo
      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mat_Q_beta

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_integ_mat'

      INTEGER                                            :: grid_2D(2), handle, handle2, iproc, &
                                                            iproc_col, iproc_row, &
                                                            mepos_in_RPA_group, &
                                                            row_col_proc_ratio, sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: group_grid_2_mepos
      LOGICAL                                            :: my_beta_case, my_blacs_ext, &
                                                            my_blacs_S_ext, my_do_im_time
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env, blacs_env_Q
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      TYPE(group_dist_d1_type)                           :: gd_ia, gd_L

      CALL timeset(routineN, handle)

      my_beta_case = .FALSE.
      IF (PRESENT(beta_case)) my_beta_case = beta_case

      my_blacs_ext = .FALSE.
      IF (PRESENT(blacs_env_ext)) my_blacs_ext = .TRUE.

      my_blacs_S_ext = .FALSE.
      IF (PRESENT(blacs_env_ext_S)) my_blacs_S_ext = .TRUE.

      my_do_im_time = .FALSE.
      IF (PRESENT(do_im_time)) my_do_im_time = do_im_time

      ! create the RPA para_env
      IF (.NOT. my_beta_case) THEN
         CALL cp_para_env_split(para_env_RPA, para_env, color_rpa_group)
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

            iproc_row = MIN(MAX(INT(SQRT(REAL(para_env_RPA%num_pe*row_col_proc_ratio, KIND=dp))), 1), para_env_RPA%num_pe) + 1
            DO iproc = 1, para_env_RPA%num_pe
               iproc_row = iproc_row - 1
               IF (MOD(para_env_RPA%num_pe, iproc_row) == 0) EXIT
            END DO

            iproc_col = para_env_RPA%num_pe/iproc_row
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

         IF (unit_nr > 0 .AND. .NOT. my_do_im_time) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Row block size:", nrow_block_mat
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MATRIX_INFO| Column block size:", ncol_block_mat
         END IF

         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_ia, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA, &
                                  nrow_block=nrow_block_mat, ncol_block=ncol_block_mat, force_block=.TRUE.)

      END IF ! external blacs_env

      ! for imaginary time we do not need the fm_mat_S
      IF (.NOT. my_do_im_time) THEN

         CALL create_group_dist(gd_ia, my_ia_start, my_ia_end, my_ia_size, para_env_RPA)
         CALL create_group_dist(gd_L, my_group_L_start, my_group_L_end, my_group_L_size, para_env_RPA)

         ! create the info array

         mepos_in_RPA_group = MOD(color_sub, integ_group_size)
         ALLOCATE (group_grid_2_mepos(0:para_env_sub%num_pe - 1, 0:integ_group_size - 1))
         group_grid_2_mepos = 0
         group_grid_2_mepos(para_env_sub%mepos, mepos_in_RPA_group) = para_env_RPA%mepos
         CALL mp_sum(group_grid_2_mepos, para_env_RPA%group)

         CALL array2fm(BIb_C_2D, fm_struct, dimen_ia, dimen_RI, para_env_RPA, &
                       my_ia_start, my_ia_end, my_group_L_start, my_group_L_end, gd_ia, gd_L, &
                       group_grid_2_mepos, para_env_sub%num_pe, ngroup, fm_mat_S, &
                       integ_group_size, color_rpa_group)

         CALL cp_fm_struct_release(fm_struct)

         ! deallocate the info array
         CALL release_group_dist(gd_L)
         CALL release_group_dist(gd_ia)

         ! mp_sum the local data across processes belonging to different RPA group.
         ! first create the para_env then mp_sum
         sub_sub_color = para_env_RPA%mepos
         CALL cp_para_env_split(para_env_exchange, para_env, sub_sub_color)

         CALL timeset(routineN//"_sum", handle2)
         CALL mp_sum(fm_mat_S%local_data, para_env_exchange%group)
         CALL timestop(handle2)

         CALL cp_para_env_release(para_env_exchange)

      END IF ! not imag. time

      CALL cp_fm_struct_release(fm_struct)

      IF (PRESENT(fm_mat_Q_gemm)) THEN
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
      END IF

      IF (PRESENT(fm_mat_Q)) THEN
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

      IF (PRESENT(fm_mat_Q_beta)) THEN
         NULLIFY (fm_mat_Q_beta)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_Q, nrow_global=dimen_RI, &
                                  ncol_global=dimen_RI, para_env=para_env_RPA)
         CALL cp_fm_create(fm_mat_Q_beta, fm_struct, name="fm_mat_Q_beta")

         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_set_all(matrix=fm_mat_Q_beta, alpha=0.0_dp)
      END IF

      ! in case we do imaginary time, we allocate fm_scaled_dm_occ and fm_scaled_dm_virt
      IF (my_do_im_time) THEN

         CALL cp_fm_create(fm_scaled_dm_occ_tau, mo_coeff%matrix_struct)
         CALL cp_fm_set_all(fm_scaled_dm_occ_tau, 0.0_dp)

         CALL cp_fm_create(fm_scaled_dm_virt_tau, mo_coeff%matrix_struct)
         CALL cp_fm_set_all(fm_scaled_dm_virt_tau, 0.0_dp)

         CALL create_occ_virt_mo_coeffs(fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, &
                                        nmo, homo)

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
! **************************************************************************************************
   SUBROUTINE create_occ_virt_mo_coeffs(fm_mo_coeff_occ, fm_mo_coeff_virt, mo_coeff, &
                                        nmo, homo)

      TYPE(cp_fm_type), OPTIONAL, POINTER                :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            mo_coeff
      INTEGER, INTENT(IN)                                :: nmo, homo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_occ_virt_mo_coeffs'

      INTEGER                                            :: handle, icol_global, irow_global

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_mo_coeff_occ, mo_coeff%matrix_struct)
      CALL cp_fm_set_all(fm_mo_coeff_occ, 0.0_dp)
      CALL cp_fm_to_fm(mo_coeff, fm_mo_coeff_occ)

      ! set all virtual MO coeffs to zero
      DO irow_global = 1, nmo
         DO icol_global = homo + 1, nmo
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

      CALL timestop(handle)

   END SUBROUTINE create_occ_virt_mo_coeffs

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
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param mo_coeff ...
!> \param fm_matrix_L_RI_metric ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
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
! **************************************************************************************************
   SUBROUTINE rpa_num_int(qs_env, Erpa, mp2_env, para_env, para_env_RPA, para_env_sub, unit_nr, &
                          homo, virtual, dimen_RI, dimen_RI_red, dimen_ia, dimen_nm_gw, &
                          Eigenval, num_integ_points, num_integ_group, color_rpa_group, &
                          fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, fm_mat_S_gw, fm_mat_R_gw, &
                          fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                          my_do_gw, do_bse, gw_corr_lev_occ, gw_corr_lev_virt, &
                          do_minimax_quad, do_im_time, fm_mo_coeff_occ, &
                          fm_mo_coeff_virt, mo_coeff, fm_matrix_L_RI_metric, &
                          fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, mat_munu, &
                          mat_P_global, &
                          t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                          starts_array_mc, ends_array_mc, &
                          starts_array_mc_block, ends_array_mc_block, &
                          matrix_s, do_kpoints_from_Gamma, kpoints, gd_array, color_sub, &
                          do_ri_sos_laplace_mp2)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(OUT)                         :: Erpa
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_RPA, para_env_sub
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, dimen_RI_red
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dimen_ia
      INTEGER, INTENT(IN)                                :: dimen_nm_gw
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S, fm_mat_Q_gemm, fm_mat_Q, &
                                                            fm_mat_S_gw
      TYPE(cp_fm_type), POINTER                          :: fm_mat_R_gw, fm_mat_S_ij_bse, &
                                                            fm_mat_S_ab_bse
      LOGICAL, INTENT(IN)                                :: my_do_gw, do_bse
      INTEGER, DIMENSION(:), INTENT(IN)                  :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, INTENT(IN)                                :: do_minimax_quad, do_im_time
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            mo_coeff
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(cp_fm_type), POINTER                          :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :), &
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
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_num_int'

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: vec_Sigma_c_gw
      INTEGER :: count_ev_sc_GW, cut_memory, group_size_P, gw_corr_lev_tot, handle, handle3, &
         ispin, iter_evGW, iter_sc_GW0, jquad, max_iter_bse, min_bsize, mm_style, nkp, &
         nkp_self_energy, nmo, nspins, num_3c_repl, num_cells_dm, num_fit_points, num_Z_vectors, &
         size_P
      INTEGER(int_8)                                     :: dbcsr_nflop
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ikp_local
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index_3c
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, prim_blk_sizes, &
                                                            RI_blk_sizes
      LOGICAL :: do_apply_ic_corr_to_gw, do_gw_im_time, do_ic_model, do_kpoints_cubic_RPA, &
         do_periodic, do_ri_Sigma_x, exit_ev_gw, first_cycle, first_cycle_periodic_correction, &
         my_open_shell, print_ic_values
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :, :)     :: has_mat_P_blocks
      REAL(KIND=dp) :: a_scaling, alpha, dbcsr_time, e_axk, e_axk_corr, eps_filter, &
         eps_filter_im_time, eps_min_trans, ext_scaling, fermi_level_offset, &
         fermi_level_offset_input, my_flop_rate, omega, omega_max_fit, omega_old, tau, tau_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: delta_corr, e_fermi, tau_tj, tau_wj, tj, &
                                                            trace_Qomega, vec_omega_fit_gw, wj, &
                                                            wkp_V, wkp_W
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: Eigenval_kp, Eigenval_last, Eigenval_scf, &
         Eigenval_scf_kp, vec_W_gw, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, &
         weights_sin_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                            B_bar_ijQ_bse_local, B_iaQ_bse_local, &
                                                            vec_Sigma_x_gw
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:, :)  :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: fm_mat_S_gw_work, fm_mat_W
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_mat_L
      TYPE(cp_fm_type), POINTER :: fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm, &
         fm_mat_RI_global_work, fm_mat_work, fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type)                                 :: mat_dm, mat_L, mat_M_P_munu_occ, &
                                                            mat_M_P_munu_virt, mat_P_global_copy, &
                                                            mat_SinvVSinv
      TYPE(dbcsr_p_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: mat_P_omega
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_berry_im_mo_mo, &
                                                            matrix_berry_re_mo_mo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega_kp
      TYPE(dbcsr_t_type)                                 :: t_3c_overl_int_ao_mo
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:)      :: t_3c_overl_int_gw_AO, &
                                                            t_3c_overl_int_gw_RI, &
                                                            t_3c_overl_nnP_ic, &
                                                            t_3c_overl_nnP_ic_reflected
      TYPE(dbcsr_type), POINTER                          :: mat_W, mat_work
      TYPE(dgemm_counter_type)                           :: dgemm_counter
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: t_3c_O_mo_compressed
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
                            num_integ_points, fm_mat_Q(1)%matrix, fm_mo_coeff_occ(1)%matrix, fm_mo_coeff_virt(1)%matrix, &
                            fm_matrix_L_RI_metric, mat_P_global, &
                            t_3c_O, matrix_s, kpoints, eps_filter_im_time, &
                            do_ri_sos_laplace_mp2, cut_memory, nkp, num_cells_dm, num_3c_repl, &
                            size_P, ikp_local, &
                            index_to_cell_3c, &
                            cell_to_index_3c, &
                            col_blk_size, &
                            do_ic_model, do_kpoints_cubic_RPA, &
                            do_kpoints_from_Gamma, do_ri_Sigma_x, my_open_shell, &
                            has_mat_P_blocks, wkp_W, wkp_V, &
                            cfm_mat_Q, fm_mat_L, fm_mat_RI_global_work, fm_mat_work, fm_mo_coeff_occ_scaled, &
                            fm_mo_coeff_virt_scaled, mat_dm, mat_L, mat_M_P_munu_occ, mat_M_P_munu_virt, mat_P_global_copy, &
                            mat_SinvVSinv, mat_P_omega, mat_P_omega_kp, mat_work)

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
                                                 para_env, fm_mat_W, fm_mat_Q(1)%matrix, &
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
         CALL get_minimax_grid(para_env, unit_nr, homo, Eigenval, num_integ_points, do_im_time, &
                               do_ri_sos_laplace_mp2,.NOT. do_ic_model, tau_tj, tau_wj, qs_env, do_gw_im_time, &
                               do_kpoints_cubic_RPA, a_scaling, e_fermi(1), tj, wj, mp2_env, &
                               weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w)
      ELSE
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

            CALL cp_fm_create(fm_mat_Q_static_bse_gemm, fm_mat_Q_gemm(1)%matrix%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q_gemm(1)%matrix, fm_mat_Q_static_bse_gemm)
            CALL cp_fm_set_all(fm_mat_Q_static_bse_gemm, 0.0_dp)

            CALL cp_fm_create(fm_mat_Q_static_bse, fm_mat_Q(1)%matrix%matrix_struct)
            CALL cp_fm_to_fm(fm_mat_Q_gemm(1)%matrix, fm_mat_Q_static_bse)
            CALL cp_fm_set_all(fm_mat_Q_static_bse, 0.0_dp)

         END IF

      END IF

      Erpa = 0.0_dp
      IF (mp2_env%ri_rpa%do_ri_axk) e_axk = 0.0_dp
      first_cycle = .TRUE.
      omega_old = 0.0_dp
      CALL dgemm_counter_init(dgemm_counter)

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
                  e_fermi(ispin) = (Eigenval(homo(ispin), ispin) + Eigenval(homo(ispin) + 1, ispin))*0.5_dp
               END DO
            END IF

            tau = 0.0_dp
            tau_old = 0.0_dp

            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(/T3,A,T66,i15)") &
               "MEMORY_INFO| Memory cut:", cut_memory
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
               "SPARSITY_INFO| Eps pgf orb for imaginary time:", mp2_env%mp2_gpw%eps_grid
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
                                     fm_scaled_dm_virt_tau, fm_mo_coeff_occ(1)%matrix, fm_mo_coeff_virt(1)%matrix, &
                                     fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                     mat_P_global, &
                                     matrix_s, &
                                     1, &
                                     t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                     starts_array_mc, ends_array_mc, &
                                     starts_array_mc_block, ends_array_mc_block, &
                                     weights_cos_tf_t_to_w, tj, tau_tj, e_fermi(1), eps_filter, alpha, &
                                     eps_filter_im_time, Eigenval(:, 1), nmo, &
                                     num_integ_points, cut_memory, &
                                     unit_nr, mp2_env, para_env, &
                                     qs_env, do_kpoints_from_Gamma, &
                                     index_to_cell_3c, cell_to_index_3c, &
                                     has_mat_P_blocks, do_ri_sos_laplace_mp2, &
                                     dbcsr_time, dbcsr_nflop)

            ! the same for open shell, use fm_mo_coeff_occ_beta and fm_mo_coeff_virt_beta
            IF (my_open_shell) THEN
               ! For SOS-MP2 we need the same calculation for alpha and beta spin independently, for RPA the sum of both
               IF (do_ri_sos_laplace_mp2) THEN
                  CALL zero_mat_P_omega(mat_P_omega(:, :, 2))

                  CALL compute_mat_P_omega(mat_P_omega(:, :, 2), fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, fm_mo_coeff_occ(2)%matrix, fm_mo_coeff_virt(2)%matrix, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           mat_P_global, &
                                           matrix_s, &
                                           2, &
                                           t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                           starts_array_mc, ends_array_mc, &
                                           starts_array_mc_block, ends_array_mc_block, &
                                           weights_cos_tf_t_to_w, tj, tau_tj, e_fermi(2), eps_filter, alpha, &
                                           eps_filter_im_time, Eigenval(:, 2), nmo, &
                                           num_integ_points, cut_memory, &
                                           unit_nr, mp2_env, para_env, &
                                           qs_env, do_kpoints_from_Gamma, &
                                           index_to_cell_3c, cell_to_index_3c, &
                                           has_mat_P_blocks, do_ri_sos_laplace_mp2, &
                                           dbcsr_time, dbcsr_nflop)
               ELSE
                  CALL compute_mat_P_omega(mat_P_omega(:, :, 1), fm_scaled_dm_occ_tau, &
                                           fm_scaled_dm_virt_tau, fm_mo_coeff_occ(2)%matrix, fm_mo_coeff_virt(2)%matrix, &
                                           fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                           mat_P_global, &
                                           matrix_s, &
                                           2, &
                                           t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                           starts_array_mc, ends_array_mc, &
                                           starts_array_mc_block, ends_array_mc_block, &
                                           weights_cos_tf_t_to_w, tj, tau_tj, e_fermi(2), eps_filter, alpha, &
                                           eps_filter_im_time, Eigenval(:, 2), nmo, &
                                           num_integ_points, cut_memory, &
                                           unit_nr, mp2_env, para_env, &
                                           qs_env, do_kpoints_from_Gamma, &
                                           index_to_cell_3c, cell_to_index_3c, &
                                           has_mat_P_blocks, do_ri_sos_laplace_mp2, &
                                           dbcsr_time, dbcsr_nflop)
               END IF ! do_ri_sos_laplace_mp2

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
                                                      fm_mat_L(1, 1)%matrix, fm_mat_Q(ispin)%matrix)
                  END DO
               END IF

            ELSE

               IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3, A, 1X, I3, 1X, A, 1X, I3)") &
                  "INTEG_INFO| Started with Integration point", jquad, "of", num_integ_points

               CALL calc_mat_Q(fm_mat_S(1)%matrix, do_ri_sos_laplace_mp2, first_cycle, count_ev_sc_GW, iter_sc_GW0, virtual(1), &
                               Eigenval(:, 1), Eigenval_last(:, 1), Eigenval_scf(:, 1), &
                               homo(1), omega, omega_old, jquad, mm_style, dimen_RI_red, dimen_ia(1), alpha, fm_mat_Q(1)%matrix, &
                               fm_mat_Q_gemm(1)%matrix, do_bse, fm_mat_Q_static_bse_gemm, dgemm_counter)

               IF (my_open_shell) THEN
                  CALL calc_mat_Q(fm_mat_S(2)%matrix, do_ri_sos_laplace_mp2, first_cycle, count_ev_sc_GW, iter_sc_GW0, virtual(2), &
                                  Eigenval(:, 2), Eigenval_last(:, 2), Eigenval_scf(:, 2), &
                                  homo(2), omega, omega_old, jquad, mm_style, &
                                  dimen_RI_red, dimen_ia(2), alpha, fm_mat_Q(2)%matrix, fm_mat_Q_gemm(2)%matrix, do_bse, &
                                  fm_mat_Q_static_bse_gemm, dgemm_counter)

                  ! For SOS-MP2 we need both matrices separately
                  IF (.NOT. do_ri_sos_laplace_mp2) THEN
                     CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_Q(1)%matrix, beta=1.0_dp, matrix_b=fm_mat_Q(2)%matrix)
                  END IF

               END IF ! open shell

            END IF ! im time

            ! Calculate AXK energy correction
            IF (mp2_env%ri_rpa%do_ri_axk) THEN
               CALL compute_axk_ener(qs_env, fm_mat_Q(1)%matrix, fm_mat_Q_gemm(1)%matrix, dimen_RI_red, dimen_ia(1), para_env_sub, &
                                     para_env_RPA, Eigenval(:, 1), fm_mat_S(1)%matrix, homo(1), virtual(1), omega, &
                                     mp2_env, mat_munu, unit_nr, e_axk_corr)

               ! Evaluate the final AXK energy correction
               e_axk = e_axk + e_axk_corr*wj(jquad)
            END IF ! do_ri_axk

            IF (do_ri_sos_laplace_mp2) THEN
               CALL SOS_MP2_postprocessing(fm_mat_Q, Erpa, tau_wj(jquad))
            ELSE
               CALL RPA_postprocessing_start(dimen_RI_red, trace_Qomega, fm_mat_Q(1)%matrix, para_env_RPA)

               IF (do_kpoints_cubic_RPA .OR. do_kpoints_from_Gamma) THEN
                  CALL RPA_postprocessing_kp(dimen_RI, num_integ_points, jquad, nkp, count_ev_sc_GW, para_env, &
                                             para_env_RPA, Erpa, tau_tj, tj, wj, weights_cos_tf_w_to_t, &
                                             wkp_W, wkp_V, do_gw_im_time, do_ri_Sigma_x, do_kpoints_from_Gamma, &
                                             do_kpoints_cubic_RPA, cfm_mat_W_kp_tau, cfm_mat_Q, ikp_local, &
                                             mat_P_omega(:, :, 1), mat_P_omega_kp, qs_env, eps_filter_im_time, unit_nr, &
                                             kpoints, fm_mat_L, fm_mat_W, fm_mat_RI_global_work, mat_SinvVSinv)
               ELSE
                  CALL RPA_postprocessing_nokp(dimen_RI_red, trace_Qomega, fm_mat_Q(1)%matrix, para_env_RPA, Erpa, wj(jquad))
               END IF
            END IF ! do_ri_sos_laplace_mp2

            ! save omega and reset the first_cycle flag
            first_cycle = .FALSE.
            omega_old = omega

            CALL timestop(handle3)

            IF (my_do_gw) THEN
               CALL GW_matrix_operations(vec_Sigma_c_gw, dimen_nm_gw, dimen_RI_red, gw_corr_lev_occ, &
                                         gw_corr_lev_virt, homo, jquad, nmo, num_fit_points, num_integ_points, &
                                         do_bse, do_im_time, do_periodic, first_cycle_periodic_correction, &
                                         fermi_level_offset, fermi_level_offset_input, &
                                         omega, Eigenval, delta_corr, tau_tj, tj, vec_omega_fit_gw, &
                                         vec_W_gw, wj, weights_cos_tf_w_to_t, fm_mat_W, fm_mat_L, &
                                         fm_mat_Q(1)%matrix, fm_mat_Q_static_bse, fm_mat_R_gw, fm_mat_S_gw, &
                                         fm_mat_S_gw_work, fm_mat_work, mo_coeff(1)%matrix, para_env, &
                                         para_env_RPA, matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, &
                                         kpoints, qs_env, mp2_env, do_kpoints_cubic_RPA, do_kpoints_from_Gamma)
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
         END IF

         IF (do_im_time) THEN

            my_flop_rate = REAL(dbcsr_nflop, dp)/(1.0E09_dp*dbcsr_time)
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T73,ES8.2)") &
               "PERFORMANCE| DBCSR total number of flops:", REAL(dbcsr_nflop*para_env%num_pe, dp)
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.2)") &
               "PERFORMANCE| DBCSR total execution time:", dbcsr_time
            IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.2)") &
               "PERFORMANCE| DBCSR flop rate (Gflops / MPI rank):", my_flop_rate

         ELSE

            CALL dgemm_counter_write(dgemm_counter, para_env, unit_nr)

         END IF

         ! GW: Compute self-energy Sigma, analytic continuation from iw -> w and correction of quasiparticle energies
         IF (my_do_gw) THEN
            IF (my_open_shell) THEN
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
                                        mo_coeff(1)%matrix, fm_mat_W, para_env, para_env_RPA, mat_dm, mat_SinvVSinv, &
                                        t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                        t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                        matrix_berry_im_mo_mo, matrix_berry_re_mo_mo, mat_W, matrix_s, &
                                        kpoints, mp2_env, qs_env, &
                                        nkp_self_energy, do_kpoints_cubic_RPA, Eigenval_kp, Eigenval_scf_kp)
            ELSE
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
                                        mo_coeff(1)%matrix, fm_mat_W, para_env, para_env_RPA, mat_dm, mat_SinvVSinv, &
                                        t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                        t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, matrix_berry_im_mo_mo, &
                                        matrix_berry_re_mo_mo, mat_W, matrix_s, &
                                        kpoints, mp2_env, qs_env, &
                                        nkp_self_energy, do_kpoints_cubic_RPA, Eigenval_kp, Eigenval_scf_kp)
            END IF

            ! if HOMO-LUMO gap differs by less than mp2_env%ri_g0w0%eps_ev_sc_iter, exit ev sc GW loop
            IF (exit_ev_gw) EXIT

         END IF ! my_do_gw if

      END DO ! evGW loop

      IF (do_ic_model) THEN

         IF (my_open_shell) THEN

            CALL calculate_ic_correction(Eigenval(:, 1), mat_SinvVSinv%matrix, &
                                         t_3c_overl_nnP_ic(1), t_3c_overl_nnP_ic_reflected(1), &
                                         gw_corr_lev_tot, &
                                         gw_corr_lev_occ(1), gw_corr_lev_virt(1), homo(1), unit_nr, &
                                         print_ic_values, para_env, do_alpha=.TRUE.)

            CALL calculate_ic_correction(Eigenval(:, 2), mat_SinvVSinv%matrix, &
                                         t_3c_overl_nnP_ic(2), t_3c_overl_nnP_ic_reflected(2), &
                                         gw_corr_lev_tot, &
                                         gw_corr_lev_occ(2), gw_corr_lev_virt(2), homo(2), unit_nr, &
                                         print_ic_values, para_env, do_beta=.TRUE.)

         ELSE

            CALL calculate_ic_correction(Eigenval(:, 1), mat_SinvVSinv%matrix, &
                                         t_3c_overl_nnP_ic(1), t_3c_overl_nnP_ic_reflected(1), &
                                         gw_corr_lev_tot, &
                                         gw_corr_lev_occ(1), gw_corr_lev_virt(1), homo(1), unit_nr, &
                                         print_ic_values, para_env)

         END IF

      END IF

      ! postprocessing after GW for Bethe-Salpeter
      IF (do_bse) THEN
         CALL mult_B_with_W_and_fill_local_3c_arrays(fm_mat_S_ij_bse, fm_mat_S_ab_bse, fm_mat_S(1)%matrix, fm_mat_Q_static_bse, &
                                                     fm_mat_Q_static_bse_gemm, &
                                                     B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                                     B_iaQ_bse_local, dimen_RI_red, homo(1), virtual(1), dimen_ia(1), &
                                                     gd_array, color_sub, para_env)

         CALL do_subspace_iterations(B_bar_ijQ_bse_local, B_abQ_bse_local, B_bar_iaQ_bse_local, &
                                     B_iaQ_bse_local, homo(1), virtual(1), num_Z_vectors, &
                                     max_iter_bse, eps_min_trans, Eigenval(:, 1), para_env)

      END IF

      IF (my_do_gw) THEN
         IF (my_open_shell) THEN
            CALL deallocate_matrices_gw(fm_mat_S_gw_work, vec_W_gw, vec_Sigma_c_gw, vec_omega_fit_gw, &
                                        mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw, &
                                        Eigenval_last, Eigenval_scf, do_periodic, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                        kpoints, do_ri_Sigma_x, vec_Sigma_x_gw,.NOT. do_im_time)
         ELSE
            CALL deallocate_matrices_gw(fm_mat_S_gw_work, vec_W_gw, vec_Sigma_c_gw, vec_omega_fit_gw, &
                                        mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw, &
                                        Eigenval_last, Eigenval_scf, do_periodic, matrix_berry_re_mo_mo, matrix_berry_im_mo_mo, &
                                        kpoints, do_ri_Sigma_x, vec_Sigma_x_gw,.NOT. do_im_time)
         END IF
      END IF

      IF (do_im_time) THEN
         CALL dealloc_im_time(fm_mo_coeff_occ, fm_mo_coeff_virt, &
                              fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, ikp_local, index_to_cell_3c, &
                              cell_to_index_3c, do_ic_model, &
                              do_kpoints_cubic_RPA, do_kpoints_from_Gamma, do_ri_Sigma_x, &
                              has_mat_P_blocks, &
                              wkp_W, wkp_V, cfm_mat_Q, fm_mat_L, fm_mat_RI_global_work, fm_mat_work, &
                              fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, mat_dm, mat_L, &
                              mat_SinvVSinv, mat_P_omega, mat_P_omega_kp, &
                              t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, mat_work, &
                              qs_env%mp2_env%ri_rpa_im_time%kpoints_G, &
                              qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma, &
                              qs_env%mp2_env%ri_rpa_im_time%kpoints_Sigma_no_xc)

         IF (my_do_gw) THEN
            CALL deallocate_matrices_gw_im_time(weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, do_ic_model, &
                                                do_kpoints_cubic_RPA, fm_mat_W, &
                                                t_3c_overl_int_ao_mo, t_3c_O_mo_compressed, t_3c_O_mo_ind, &
                                                t_3c_overl_int_gw_RI, t_3c_overl_int_gw_AO, &
                                                t_3c_overl_nnP_ic, t_3c_overl_nnP_ic_reflected, &
                                                mat_W, ikp_local, cfm_mat_W_kp_tau)
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

      CALL timestop(handle)

   END SUBROUTINE rpa_num_int

END MODULE rpa_main
