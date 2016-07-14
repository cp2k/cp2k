!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rountines for RPA with imaginary time
!> \par History
!>      10.2015 created [Jan Wilhelm]
! **************************************************************************************************
MODULE rpa_im_time
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_allocate_matrix_set, cp_dbcsr_copy, cp_dbcsr_create, &
        cp_dbcsr_deallocate_matrix_set, cp_dbcsr_filter, cp_dbcsr_finalize, cp_dbcsr_get_info, &
        cp_dbcsr_get_occupation, cp_dbcsr_get_stored_coordinates, cp_dbcsr_init, &
        cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
        cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_p_type, &
        cp_dbcsr_reserve_blocks, cp_dbcsr_set, cp_dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_sendrecv,&
                                              mp_sync,&
                                              mp_waitall
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type,&
                                              two_dim_int_array
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time'

   PUBLIC :: compute_fm_mat_Q_omega, &
             print_occupation_2c, &
             print_occupation_3c, &
             fill_fm_with_2d_array, &
             setup_mat_for_mem_cut_3c, &
             communicate_buffer
CONTAINS

! **************************************************************************************************
!> \brief compute the matrix Q(it) (intermediate) and Fourier transform it
!>        directly to fm_mat_Q_omega(iw) (output)
!> \param fm_mat_Q_omega ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mat_Q ...
!> \param mat_P_local ...
!> \param mat_P_global ...
!> \param mat_M_mu_Pnu_occ ...
!> \param mat_M_mu_Pnu_virt ...
!> \param matrix_s ...
!> \param mat_M_muP_occ ...
!> \param mat_M_muP_virt ...
!> \param mat_3c_overl_int_cut ...
!> \param mat_dm_loc_occ_cut ...
!> \param mat_dm_loc_virt_cut ...
!> \param tau_wj_2d ...
!> \param tj ...
!> \param wj ...
!> \param tau_tj ...
!> \param do_minimax_quad ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param eps_comm ...
!> \param a_scaling ...
!> \param alpha ...
!> \param eps_filter_im_time ...
!> \param Eigenval ...
!> \param nmo ...
!> \param n_group_col ...
!> \param dimen_RI ...
!> \param group_size_P ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param starts_array_prim_col ...
!> \param ends_array_prim_col ...
!> \param starts_array_prim_row ...
!> \param ends_array_prim_row ...
!> \param starts_array_prim_fullcol ...
!> \param ends_array_prim_fullcol ...
!> \param starts_array_prim_fullrow ...
!> \param ends_array_prim_fullrow ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param offset_combi_block ...
!> \param starts_array_cm ...
!> \param ends_array_cm ...
!> \param mepos_P_from_RI_row ...
!> \param non_zero_rows_3c ...
!> \param buffer_mat_M ...
! **************************************************************************************************
   SUBROUTINE compute_fm_mat_Q_omega(fm_mat_Q_omega, fm_scaled_dm_occ_tau, &
                                     fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                     fm_mat_Q, mat_P_local, mat_P_global, &
                                     mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, matrix_s, mat_M_muP_occ, &
                                     mat_M_muP_virt, mat_3c_overl_int_cut, mat_dm_loc_occ_cut, &
                                     mat_dm_loc_virt_cut, tau_wj_2d, tj, wj, tau_tj, &
                                     do_minimax_quad, e_fermi, eps_filter, eps_comm, a_scaling, alpha, &
                                     eps_filter_im_time, Eigenval, nmo, n_group_col, dimen_RI, &
                                     group_size_P, num_integ_points, jquad, cut_memory, cut_RI, unit_nr, &
                                     mp2_env, para_env, para_env_sub, &
                                     starts_array_prim_col, ends_array_prim_col, &
                                     starts_array_prim_row, ends_array_prim_row, &
                                     starts_array_prim_fullcol, ends_array_prim_fullcol, &
                                     starts_array_prim_fullrow, ends_array_prim_fullrow, &
                                     my_group_L_starts_im_time, my_group_L_sizes_im_time, &
                                     offset_combi_block, starts_array_cm, &
                                     ends_array_cm, mepos_P_from_RI_row, non_zero_rows_3c, &
                                     buffer_mat_M)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_Q_omega, fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            fm_mat_Q
      TYPE(cp_dbcsr_p_type)                              :: mat_P_local, mat_P_global
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, &
                                                            matrix_s
      TYPE(cp_dbcsr_p_type)                              :: mat_M_muP_occ, mat_M_muP_virt
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: mat_3c_overl_int_cut, &
                                                            mat_dm_loc_occ_cut, mat_dm_loc_virt_cut
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: tau_wj_2d
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tj, wj, tau_tj
      LOGICAL                                            :: do_minimax_quad
      REAL(KIND=dp)                                      :: e_fermi, eps_filter, eps_comm, &
                                                            a_scaling, alpha, eps_filter_im_time
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: nmo, n_group_col, dimen_RI, &
                                                            group_size_P, num_integ_points, jquad, &
                                                            cut_memory, cut_RI, unit_nr
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: starts_array_prim_col, ends_array_prim_col, &
         starts_array_prim_row, ends_array_prim_row, starts_array_prim_fullcol, &
         ends_array_prim_fullcol, starts_array_prim_fullrow, ends_array_prim_fullrow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_starts_im_time, &
                                                            my_group_L_sizes_im_time
      TYPE(two_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: offset_combi_block
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: starts_array_cm, ends_array_cm, &
                                                            mepos_P_from_RI_row
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_rows_3c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_mat_M

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_mat_Q_omega', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle4, handle5, i_cut_RI, &
                                                            i_mem, iquad, j_mem, my_group_L_size, &
                                                            my_group_L_start
      LOGICAL                                            :: first_cycle_im_time, &
                                                            first_cycle_omega_loop, memory_info
      REAL(KIND=dp)                                      :: omega, omega_old, tau, weight, weight_old
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_dm_occ_global, mat_dm_virt_global

      ! has to be set, otherwise, compiler complains
      tau = 0.0_dp

      memory_info = mp2_env%ri_rpa_im_time%memory_info

      CALL timeset(routineN, handle)

      CALL precompute_fm_scaled_dm_tau(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, nmo, &
                                       fm_mo_coeff_occ, fm_mo_coeff_virt, mat_dm_occ_global, mat_dm_virt_global, &
                                       matrix_s, Eigenval, e_fermi, eps_filter, memory_info, unit_nr, para_env)

      first_cycle_im_time = .TRUE.

      ! iteration over time points
      DO jquad = 1, num_integ_points

         tau = tau_tj(jquad)

         IF (memory_info .AND. unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "RPA_IM_TIME_INFO| Started with time point: ", jquad

         CALL timeset(routineN//"_im_time_repl_subgr", handle4)

         CALL mp_sync(para_env%group)

         CALL timestop(handle4)

         DO j_mem = 1, cut_memory

            CALL replicate_dm_to_subgroup(para_env, para_env_sub, mat_dm_virt_global(jquad)%matrix, nmo, &
                                          mat_dm_loc_virt_cut, starts_array_cm, ends_array_cm, &
                                          j_mem, cut_RI, non_zero_rows_3c)

            DO i_mem = 1, cut_memory

               CALL replicate_dm_to_subgroup(para_env, para_env_sub, mat_dm_occ_global(jquad)%matrix, nmo, &
                                             mat_dm_loc_occ_cut, starts_array_cm, ends_array_cm, &
                                             i_mem, cut_RI, non_zero_rows_3c)

               CALL timeset(routineN//"_calc_M", handle5)

               DO i_cut_RI = 1, cut_RI

                  my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
                  my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

                  ! D^occ*(munuP)
                  CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_dm_loc_occ_cut(i_cut_RI, i_mem)%matrix, &
                                         mat_3c_overl_int_cut(i_cut_RI, j_mem)%matrix, &
                                         0.0_dp, mat_M_mu_Pnu_occ(i_cut_RI)%matrix, &
                                         ! first_row=starts_array_cm(i_mem), &
                                         ! last_row=ends_array_cm(i_mem), &
                                         ! first_column=(starts_array_cm(j_mem)-1)*my_group_L_size+1, &
                                         ! last_column=ends_array_cm(j_mem)*my_group_L_size, &
                                         filter_eps=eps_filter*SQRT(REAL(my_group_L_size)))

                  ! D^virt*(munuP)
                  CALL cp_dbcsr_multiply("N", "N", alpha/2.0_dp, mat_dm_loc_virt_cut(i_cut_RI, j_mem)%matrix, &
                                         mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix, &
                                         0.0_dp, mat_M_mu_Pnu_virt(i_cut_RI)%matrix, &
                                         ! first_row=starts_array_cm(j_mem), &
                                         ! last_row=ends_array_cm(j_mem), &
                                         ! first_column=(starts_array_cm(i_mem)-1)*my_group_L_size+1, &
                                         ! last_column=ends_array_cm(i_mem)*my_group_L_size, &
                                         filter_eps=eps_filter*SQRT(REAL(my_group_L_size)))

               END DO

               ! free the occupied density matrix
               DO i_cut_RI = 1, cut_RI
                  CALL cp_dbcsr_set(mat_dm_loc_occ_cut(i_cut_RI, i_mem)%matrix, 0.0_dp)
                  CALL cp_dbcsr_filter(mat_dm_loc_occ_cut(i_cut_RI, i_mem)%matrix, 1.0_dp)
               END DO

               IF (first_cycle_im_time .AND. memory_info) THEN
                  CALL print_occupation_3c(mat_M_mu_Pnu_occ, unit_nr, "M occ before comm", para_env, cut_RI)
                  CALL print_occupation_3c(mat_M_mu_Pnu_virt, unit_nr, "M virt before comm", para_env, cut_RI)
               END IF

               CALL mp_sync(para_env%group)

               CALL timestop(handle5)

               CALL timeset(routineN//"_fill_M_tot", handle5)

               CALL fill_mat_M_muP_from_mat_M_munu_P(mat_M_muP_occ, mat_M_mu_Pnu_occ, para_env, my_group_L_starts_im_time, &
                                                     my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                                     n_group_col, cut_RI, eps_comm, &
                                                     starts_array_prim_row, ends_array_prim_row, &
                                                     starts_array_prim_col, ends_array_prim_col, &
                                                     offset_combi_block(i_mem, j_mem)%array, starts_array_prim_fullcol, &
                                                     ends_array_prim_fullcol, starts_array_prim_fullrow, &
                                                     ends_array_prim_fullrow, group_size_P, &
                                                     dimen_RI, i_mem, j_mem, buffer_mat_M, eps_filter_im_time, &
                                                     do_occ=.TRUE., do_virt=.FALSE.)

               CALL fill_mat_M_muP_from_mat_M_munu_P(mat_M_muP_virt, mat_M_mu_Pnu_virt, para_env, my_group_L_starts_im_time, &
                                                     my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                                     n_group_col, cut_RI, eps_comm, &
                                                     starts_array_prim_row, ends_array_prim_row, &
                                                     starts_array_prim_col, ends_array_prim_col, &
                                                     offset_combi_block(i_mem, j_mem)%array, starts_array_prim_fullcol, &
                                                     ends_array_prim_fullcol, starts_array_prim_fullrow, &
                                                     ends_array_prim_fullrow, group_size_P, &
                                                     dimen_RI, i_mem, j_mem, buffer_mat_M, eps_filter_im_time, &
                                                     do_occ=.FALSE., do_virt=.TRUE.)

               IF (first_cycle_im_time .AND. memory_info) THEN
                  CALL print_occupation_2c(mat_M_muP_occ%matrix, unit_nr, "M occ after comm", para_env)
                  CALL print_occupation_2c(mat_M_muP_virt%matrix, unit_nr, "M virt after comm", para_env)
               END IF

               CALL mp_sync(para_env%group)

               CALL timestop(handle5)

               CALL timeset(routineN//"_calc_P", handle5)

               ! P_RT = sum_mu sigma M^occ_P_mu_sigma M^virt_R_mu_sigma
               CALL cp_dbcsr_multiply("N", "T", 1.0_dp, mat_M_muP_occ%matrix, mat_M_muP_virt%matrix, &
                                      0.0_dp, mat_P_local%matrix, filter_eps=eps_filter_im_time/REAL(cut_memory))

               IF (first_cycle_im_time .AND. memory_info) THEN
                  CALL print_occupation_2c(mat_P_local%matrix, unit_nr, "local P matrix", para_env)
               END IF

               ! release memory
               CALL cp_dbcsr_set(mat_M_muP_occ%matrix, 0.0_dp)
               CALL cp_dbcsr_filter(mat_M_muP_occ%matrix, 1.0_dp)

               CALL cp_dbcsr_set(mat_M_muP_virt%matrix, 0.0_dp)
               CALL cp_dbcsr_filter(mat_M_muP_virt%matrix, 1.0_dp)

               IF (first_cycle_im_time .AND. memory_info) THEN
                  CALL print_occupation_2c(mat_P_global%matrix, unit_nr, "global P matrix", para_env)
               END IF

               CALL mp_sync(para_env%group)

               CALL timestop(handle5)

               CALL cp_dbcsr_set(mat_P_global%matrix, 0.0_dp)

               CALL fill_mat_P_global_from_mat_P_local(mat_P_global, mat_P_local, para_env)

               CALL cp_dbcsr_set(mat_P_local%matrix, 0.0_dp)
               CALL cp_dbcsr_filter(mat_P_local%matrix, 1.0_dp)

               CALL copy_dbcsr_to_fm(mat_P_global%matrix, fm_mat_Q)

               CALL timeset(routineN//"_Fourier_transform", handle5)

               ! Fourier transform of P(it) to P(iw) (fm_mat_Q_omega)
               first_cycle_omega_loop = .TRUE.

               DO iquad = 1, num_integ_points

                  IF (do_minimax_quad) THEN
                     omega = tj(iquad)
                     weight = tau_wj_2d(iquad, jquad)
                  ELSE
                     omega = a_scaling/TAN(tj(iquad))
                     weight = wj(jquad)
                  END IF

                  IF (first_cycle_omega_loop) THEN
                     ! no multiplication with 2.0 as in Kresses paper (Kaltak, JCTC 10, 2498 (2014), Eq. 12)
                     ! because this factor is already absorbed in the weight w_j
                     CALL cp_fm_scale(COS(omega*tau)*weight, fm_mat_Q)
                  ELSE
                     CALL cp_fm_scale(COS(omega*tau)/COS(omega_old*tau)*weight/weight_old, fm_mat_Q)
                  END IF

                  CALL cp_fm_scale_and_add(1.0_dp, fm_mat_Q_omega(iquad)%matrix, 1.0_dp, fm_mat_Q)

                  first_cycle_omega_loop = .FALSE.

                  omega_old = omega
                  weight_old = weight

               END DO

               CALL timestop(handle5)

               first_cycle_im_time = .FALSE.

            END DO ! cut memory i

            ! free the virtual density matrix
            DO i_cut_RI = 1, cut_RI
               CALL cp_dbcsr_set(mat_dm_loc_virt_cut(i_cut_RI, j_mem)%matrix, 0.0_dp)
               CALL cp_dbcsr_filter(mat_dm_loc_virt_cut(i_cut_RI, j_mem)%matrix, 1.0_dp)
            END DO

         END DO ! cut memory j

      END DO ! time points

      CALL clean_up(mat_dm_occ_global, mat_dm_virt_global)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param mat_dm_global ...
!> \param nmo ...
!> \param mat_dm_loc_cut ...
!> \param starts_array_cm ...
!> \param ends_array_cm ...
!> \param i_mem ...
!> \param cut_RI ...
!> \param non_zero_rows_3c ...
! **************************************************************************************************
   SUBROUTINE replicate_dm_to_subgroup(para_env, para_env_sub, mat_dm_global, nmo, &
                                       mat_dm_loc_cut, starts_array_cm, ends_array_cm, &
                                       i_mem, cut_RI, non_zero_rows_3c)

      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(cp_dbcsr_type), POINTER                       :: mat_dm_global
      INTEGER                                            :: nmo
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: mat_dm_loc_cut
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: starts_array_cm, ends_array_cm
      INTEGER                                            :: i_mem, cut_RI
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_rows_3c

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_dm_to_subgroup', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_counter, block_offset, block_size, col, col_from_buffer, col_offset, &
         col_size, handle, handle1, i_block, i_cut_RI, i_entry, i_mepos, igroup, imepos, &
         imepos_sub, msg_offset, nblkrows_total, ngroup, num_blocks, num_pe_sub, offset, row, &
         row_from_buffer, row_offset, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_counter, cols_to_allocate_all, &
         counter_blk_to_alloc, entry_counter, num_entries_blocks_rec, num_entries_blocks_send, &
         row_block_from_index, rows_to_allocate_all, sizes_rec, sizes_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: cols_to_allocate, rows_to_allocate
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_offset, row_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset("get_sizes_D", handle1)

      NULLIFY (data_block)

      ALLOCATE (num_entries_blocks_send(0:2*para_env%num_pe-1))
      num_entries_blocks_send(:) = 0

      ALLOCATE (num_entries_blocks_rec(0:2*para_env%num_pe-1))
      num_entries_blocks_rec(:) = 0

      ngroup = para_env%num_pe/para_env_sub%num_pe

      num_pe_sub = para_env_sub%num_pe

      CALL cp_dbcsr_iterator_start(iter, mat_dm_global)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

         ! check whether block is in the range of the memory cutoff
         IF (row_offset+row_size-1 >= starts_array_cm(i_mem) .AND. row_offset < ends_array_cm(i_mem)) THEN

            CALL cp_dbcsr_get_stored_coordinates(mat_dm_loc_cut(1, 1)%matrix, row, col, imepos_sub)

            DO igroup = 0, ngroup-1

               imepos = imepos_sub+igroup*num_pe_sub

               num_entries_blocks_send(2*imepos) = num_entries_blocks_send(2*imepos)+row_size*col_size
               num_entries_blocks_send(2*imepos+1) = num_entries_blocks_send(2*imepos+1)+1

            END DO

         END IF

      END DO

      CALL cp_dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("send_sizes_D", handle1)

      IF (para_env%num_pe > 1) THEN
         CALL mp_alltoall(num_entries_blocks_send, num_entries_blocks_rec, 2, para_env%group)
      ELSE
         num_entries_blocks_rec(0:1) = num_entries_blocks_send(0:1)
      END IF

      CALL timestop(handle1)

      CALL timeset("get_data_D", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe-1))
      ALLOCATE (buffer_send(0:para_env%num_pe-1))

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_blocks_rec(2*imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_blocks_send(2*imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_entries_blocks_rec(2*imepos+1), 3))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_entries_blocks_send(2*imepos+1), 3))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe-1))
      entry_counter(:) = 0

      ALLOCATE (blk_counter(0:para_env%num_pe-1))
      blk_counter = 0

      CALL cp_dbcsr_iterator_start(iter, mat_dm_global)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

         ! check whether block is in the range of the memory cutoff
         IF (row_offset+row_size-1 >= starts_array_cm(i_mem) .AND. row_offset < ends_array_cm(i_mem)) THEN

            CALL cp_dbcsr_get_stored_coordinates(mat_dm_loc_cut(1, 1)%matrix, row, col, imepos_sub)

            DO igroup = 0, ngroup-1

               imepos = imepos_sub+igroup*num_pe_sub

               msg_offset = entry_counter(imepos)

               block_size = row_size*col_size

               buffer_send(imepos)%msg(msg_offset+1:msg_offset+block_size) = &
                  RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

               entry_counter(imepos) = entry_counter(imepos)+block_size

               blk_counter(imepos) = blk_counter(imepos)+1

               block_offset = blk_counter(imepos)

               buffer_send(imepos)%indx(block_offset, 1) = row
               buffer_send(imepos)%indx(block_offset, 2) = col
               buffer_send(imepos)%indx(block_offset, 3) = msg_offset

            END DO

         END IF

      END DO

      CALL cp_dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("send_data_D", handle1)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      ALLOCATE (sizes_rec(0:para_env%num_pe-1))
      ALLOCATE (sizes_send(0:para_env%num_pe-1))

      DO imepos = 0, para_env%num_pe-1

         sizes_send(imepos) = num_entries_blocks_send(2*imepos)
         sizes_rec(imepos) = num_entries_blocks_rec(2*imepos)

      END DO

      CALL communicate_buffer(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send, req_array)

      DEALLOCATE (req_array, sizes_rec, sizes_send)

      CALL timestop(handle1)

      CALL timeset("row_block_from_index", handle1)

      CALL cp_dbcsr_get_info(mat_dm_loc_cut(1, 1)%matrix, &
                             nblkrows_total=nblkrows_total, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_size)

      ALLOCATE (row_block_from_index(nmo))
      row_block_from_index = 0

      DO i_entry = 1, nmo
         DO i_block = 1, nblkrows_total

            IF (i_entry >= row_blk_offset(i_block) .AND. i_entry <= row_blk_offset(i_block)+row_blk_size(i_block)-1) THEN

               row_block_from_index(i_entry) = i_block

            END IF

         END DO
      END DO

      CALL timestop(handle1)

      CALL timeset("reserve_blocks_D", handle1)

      num_blocks = 0

      ! get the number of blocks, which have to be allocated
      DO imepos = 0, para_env%num_pe-1
         num_blocks = num_blocks+num_entries_blocks_rec(2*imepos+1)
      END DO

      ALLOCATE (rows_to_allocate_all(num_blocks))
      rows_to_allocate_all = 0

      ALLOCATE (cols_to_allocate_all(num_blocks))
      cols_to_allocate_all = 0

      block_counter = 0

      DO i_mepos = 0, para_env%num_pe-1

         DO i_block = 1, num_entries_blocks_rec(2*i_mepos+1)

            block_counter = block_counter+1

            rows_to_allocate_all(block_counter) = buffer_rec(i_mepos)%indx(i_block, 1)
            cols_to_allocate_all(block_counter) = buffer_rec(i_mepos)%indx(i_block, 2)

         END DO

      END DO

      ALLOCATE (rows_to_allocate(num_blocks, cut_RI))
      rows_to_allocate = 0
      ALLOCATE (cols_to_allocate(num_blocks, cut_RI))
      cols_to_allocate = 0

      ALLOCATE (counter_blk_to_alloc(cut_RI))
      counter_blk_to_alloc = 0

      DO i_cut_RI = 1, cut_RI

         DO i_block = 1, num_blocks

            row = rows_to_allocate_all(i_block)
            col = cols_to_allocate_all(i_block)

            IF (ANY(non_zero_rows_3c(i_cut_RI, :, :) == col)) THEN

               counter_blk_to_alloc(i_cut_RI) = counter_blk_to_alloc(i_cut_RI)+1
               block_counter = counter_blk_to_alloc(i_cut_RI)

               rows_to_allocate(block_counter, i_cut_RI) = row
               cols_to_allocate(block_counter, i_cut_RI) = col

            END IF

         END DO

         CALL cp_dbcsr_set(mat_dm_loc_cut(i_cut_RI, i_mem)%matrix, 0.0_dp)
         CALL cp_dbcsr_filter(mat_dm_loc_cut(i_cut_RI, i_mem)%matrix, 1.0_dp)
         CALL cp_dbcsr_reserve_blocks(mat_dm_loc_cut(i_cut_RI, i_mem)%matrix, &
                                      rows=rows_to_allocate(1:counter_blk_to_alloc(i_cut_RI), i_cut_RI), &
                                      cols=cols_to_allocate(1:counter_blk_to_alloc(i_cut_RI), i_cut_RI))
         CALL cp_dbcsr_finalize(mat_dm_loc_cut(i_cut_RI, i_mem)%matrix)
         CALL cp_dbcsr_set(mat_dm_loc_cut(i_cut_RI, i_mem)%matrix, 0.0_dp)

      END DO

      DEALLOCATE (rows_to_allocate_all, cols_to_allocate_all)

      DO i_cut_RI = 1, cut_RI

         CALL cp_dbcsr_iterator_start(iter, mat_dm_loc_cut(i_cut_RI, i_mem)%matrix)

         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

            CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                              row_size=row_size, col_size=col_size)

            DO imepos = 0, para_env%num_pe-1

               DO i_block = 1, num_entries_blocks_rec(2*imepos+1)

                  row_from_buffer = buffer_rec(imepos)%indx(i_block, 1)
                  col_from_buffer = buffer_rec(imepos)%indx(i_block, 2)
                  offset = buffer_rec(imepos)%indx(i_block, 3)

                  IF (row == row_from_buffer .AND. col == col_from_buffer) THEN

                     data_block(1:row_size, 1:col_size) = &
                        RESHAPE(buffer_rec(imepos)%msg(offset+1:offset+row_size*col_size), &
                                (/row_size, col_size/))

                  END IF

               END DO

            END DO

         END DO ! blocks

         CALL cp_dbcsr_iterator_stop(iter)

      END DO ! cut_RI

      CALL timestop(handle1)

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (num_entries_blocks_send, num_entries_blocks_rec)
      DEALLOCATE (entry_counter)
      DEALLOCATE (blk_counter)
      DEALLOCATE (row_block_from_index)
      DEALLOCATE (rows_to_allocate, cols_to_allocate)
      DEALLOCATE (counter_blk_to_alloc)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_P_global ...
!> \param mat_P_local ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE fill_mat_P_global_from_mat_P_local(mat_P_global, mat_P_local, para_env)

      TYPE(cp_dbcsr_p_type)                              :: mat_P_global, mat_P_local
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_mat_P_global_from_mat_P_local', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: block, block_size, col, col_size, &
                                                            handle, handle1, imepos, offset, &
                                                            rec_counter, row, row_size, &
                                                            send_counter
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, entry_counter, num_blocks_rec, &
         num_blocks_send, num_entries_rec, num_entries_send, sizes_rec, sizes_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset("fill_Q_2_get_coord", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe-1))
      ALLOCATE (buffer_send(0:para_env%num_pe-1))

      ALLOCATE (num_entries_rec(0:para_env%num_pe-1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe-1))
      ALLOCATE (num_entries_send(0:para_env%num_pe-1))
      ALLOCATE (num_blocks_send(0:para_env%num_pe-1))
      num_entries_rec = 0
      num_blocks_rec = 0
      num_entries_send = 0
      num_blocks_send = 0

      CALL cp_dbcsr_iterator_start(iter, mat_P_local%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size)

         CALL cp_dbcsr_get_stored_coordinates(mat_P_global%matrix, row, col, imepos)

         num_entries_send(imepos) = num_entries_send(imepos)+row_size*col_size
         num_blocks_send(imepos) = num_blocks_send(imepos)+1

      END DO

      CALL cp_dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_comm_size", handle1)

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe-1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe-1))

         DO imepos = 0, para_env%num_pe-1

            sizes_send(2*imepos) = num_entries_send(imepos)
            sizes_send(2*imepos+1) = num_blocks_send(imepos)

         END DO

         CALL mp_alltoall(sizes_send, sizes_rec, 2, para_env%group)

         DO imepos = 0, para_env%num_pe-1
            num_entries_rec(imepos) = sizes_rec(2*imepos)
            num_blocks_rec(imepos) = sizes_rec(2*imepos+1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send(0)
         num_blocks_rec(0) = num_blocks_send(0)

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_fill_buffer", handle1)

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 5))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send(imepos), 5))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (block_counter(0:para_env%num_pe-1))
      block_counter(:) = 0

      ALLOCATE (entry_counter(0:para_env%num_pe-1))
      entry_counter(:) = 0

      ! fill buffer_send
      CALL cp_dbcsr_iterator_start(iter, mat_P_local%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size)

         CALL cp_dbcsr_get_stored_coordinates(mat_P_global%matrix, row, col, imepos)

         block_size = row_size*col_size

         offset = entry_counter(imepos)

         buffer_send(imepos)%msg(offset+1:offset+block_size) = &
            RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

         block = block_counter(imepos)+1

         buffer_send(imepos)%indx(block, 1) = row
         buffer_send(imepos)%indx(block, 2) = col
         buffer_send(imepos)%indx(block, 3) = offset

         entry_counter(imepos) = entry_counter(imepos)+block_size

         block_counter(imepos) = block_counter(imepos)+1

      END DO

      CALL cp_dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_comm_data", handle1)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      IF (para_env%num_pe > 1) THEN

!         ! send around the data and indices
!         DO imepos = 0, para_env%num_pe-1
!            CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(imepos, 1), tag=4)
!            CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(imepos, 2), tag=4)
!            CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(imepos, 3), tag=7)
!            CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(imepos, 4), tag=7)
!         END DO
!
!         CALL mp_waitall(req_array(:, 1:4))

         send_counter = 0
         rec_counter = 0

!         ! send around the data and indices
!         DO imepos = 0, para_env%num_pe-1
!            IF (num_entries_send(imepos) > 0) THEN
!               send_counter = send_counter+1
!               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
!            END IF
!            IF (num_entries_rec(imepos) > 0) THEN
!               rec_counter = rec_counter+1
!               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
!            END IF
!            IF (num_entries_send(imepos) > 0) THEN
!               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
!            END IF
!            IF (num_entries_rec(imepos) > 0) THEN
!               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
!            END IF
!         END DO

         DO imepos = 0, para_env%num_pe-1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter+1
               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
            END IF
            IF (num_entries_rec(imepos) > 0) THEN
               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO imepos = 0, para_env%num_pe-1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter+1
               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
            END IF
            IF (num_entries_send(imepos) > 0) THEN
               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         CALL mp_waitall(req_array(1:send_counter, 1:2))
         CALL mp_waitall(req_array(1:rec_counter, 3:4))

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_set_blocks", handle1)

      ! fill mat_P_global
      CALL cp_dbcsr_iterator_start(iter, mat_P_global%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size)

         block_size = row_size*col_size

         DO imepos = 0, para_env%num_pe-1

            DO block = 1, num_blocks_rec(imepos)

               IF (row == buffer_rec(imepos)%indx(block, 1)) THEN

                  IF (col == buffer_rec(imepos)%indx(block, 2)) THEN

                     offset = buffer_rec(imepos)%indx(block, 3)

                     data_block(1:row_size, 1:col_size) = data_block(1:row_size, 1:col_size)+ &
                                                          RESHAPE(buffer_rec(imepos)%msg(offset+1:offset+row_size*col_size), &
                                                                  (/row_size, col_size/))

                  END IF

               END IF

            END DO

         END DO

      END DO

      CALL cp_dbcsr_iterator_stop(iter)

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)

      DEALLOCATE (block_counter, entry_counter)

      DEALLOCATE (req_array)

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_munu_array ...
!> \param unit_nr ...
!> \param matrix_name ...
!> \param para_env ...
!> \param cut_RI ...
! **************************************************************************************************
   SUBROUTINE print_occupation_3c(mat_munu_array, unit_nr, matrix_name, para_env, cut_RI)

      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_munu_array
      INTEGER                                            :: unit_nr
      CHARACTER(len=*)                                   :: matrix_name
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: cut_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_occupation_3c', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cut_RI, imepos
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp)                                      :: local_occupation, max_occupation, &
                                                            min_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupation
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      local_occupation = 0.0_dp

      DO i_cut_RI = 1, cut_RI

         local_occupation = local_occupation+cp_dbcsr_get_occupation(mat_munu_array(i_cut_RI)%matrix)

      END DO

      local_occupation = local_occupation/REAL(cut_RI)

      ALLOCATE (buffer_send(0:para_env%num_pe-1))
      ALLOCATE (buffer_rec(0:para_env%num_pe-1))

      DO imepos = 0, para_env%num_pe-1
         ALLOCATE (buffer_send(imepos)%msg(1))
         buffer_send(imepos)%msg(1) = local_occupation

         ALLOCATE (buffer_rec(imepos)%msg(1))
      END DO

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (req_array(0:para_env%num_pe-1, 2))

         DO imepos = 0, para_env%num_pe-1
            CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(imepos, 1), tag=2)
            CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(imepos, 2), tag=2)
         END DO

         CALL mp_waitall(req_array(:, 1:2))

         DEALLOCATE (req_array)

      ELSE

         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      ALLOCATE (occupation(0:para_env%num_pe-1))

      DO imepos = 0, para_env%num_pe-1

         occupation(imepos) = buffer_rec(imepos)%msg(1)

      END DO

      max_occupation = MAXVAL(occupation)

      min_occupation = MINVAL(occupation)

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%msg)
      END DO

      DEALLOCATE (buffer_send, buffer_rec, occupation)

      ! print fraction of non-zero blocks
      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T64,ES7.1,A,T74,ES7.1)") &
         "MATRIX_INFO| Frac. of non-zero bl. in ", matrix_name, ":", min_occupation, ' -', max_occupation

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_munu ...
!> \param unit_nr ...
!> \param matrix_name ...
!> \param para_env ...
!> \param one_number ...
! **************************************************************************************************
   SUBROUTINE print_occupation_2c(mat_munu, unit_nr, matrix_name, para_env, one_number)
      TYPE(cp_dbcsr_type), POINTER                       :: mat_munu
      INTEGER                                            :: unit_nr
      CHARACTER(len=*)                                   :: matrix_name
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, OPTIONAL                                  :: one_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_occupation_2c', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imepos
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL                                            :: my_one_number
      REAL(KIND=dp)                                      :: local_occupation, max_occupation, &
                                                            min_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupation
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      IF (PRESENT(one_number)) THEN
         my_one_number = one_number
      ELSE
         my_one_number = .FALSE.
      END IF

      local_occupation = cp_dbcsr_get_occupation(mat_munu)

      ALLOCATE (buffer_send(0:para_env%num_pe-1))
      ALLOCATE (buffer_rec(0:para_env%num_pe-1))

      DO imepos = 0, para_env%num_pe-1
         ALLOCATE (buffer_send(imepos)%msg(1))
         buffer_send(imepos)%msg(1) = local_occupation

         ALLOCATE (buffer_rec(imepos)%msg(1))
      END DO

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (req_array(0:para_env%num_pe-1, 2))

         DO imepos = 0, para_env%num_pe-1
            CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(imepos, 1), tag=2)
            CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(imepos, 2), tag=2)
         END DO

         CALL mp_waitall(req_array(:, 1:2))

         DEALLOCATE (req_array)

      ELSE

         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      ALLOCATE (occupation(0:para_env%num_pe-1))

      DO imepos = 0, para_env%num_pe-1

         occupation(imepos) = buffer_rec(imepos)%msg(1)

      END DO

      max_occupation = MAXVAL(occupation)

      min_occupation = MINVAL(occupation)

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%msg)
      END DO

      DEALLOCATE (buffer_send, buffer_rec, occupation)

      IF (my_one_number) THEN

         ! print fraction of non-zero blocks
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T74,ES7.1)") &
            "MATRIX_INFO| Frac. of non-zero bl. in ", matrix_name, ":", max_occupation

      ELSE

         ! print fraction of non-zero blocks
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T64,ES7.1,A,T74,ES7.1)") &
            "MATRIX_INFO| Frac. of non-zero bl. in ", matrix_name, ":", min_occupation, ' -', max_occupation

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief fill_fm_with_2d_array
!> \param fm_mat_Q ...
!> \param trace_2d ...
!> \param para_env ...
!> \param my_group_L_start ...
!> \param my_group_L_size ...
!> \param dimen_RI ...
!> \param ngroup ...
! **************************************************************************************************
   SUBROUTINE fill_fm_with_2d_array(fm_mat_Q, trace_2d, para_env, my_group_L_start, my_group_L_size, &
                                    dimen_RI, ngroup)

      TYPE(cp_fm_type), POINTER                          :: fm_mat_Q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: trace_2d
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: my_group_L_start, my_group_L_size, &
                                                            dimen_RI, ngroup

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_with_2d_array', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igroup, my_group_K_size, &
                                                            my_group_K_start, proc_rec_M, &
                                                            proc_send_M
      INTEGER, DIMENSION(:), POINTER                     :: my_group_K_rec, my_group_K_send
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: trace_2d_rec, trace_2d_send

      CALL timeset(routineN, handle)

      CALL cp_fm_set_submatrix(fm=fm_mat_Q, new_values=trace_2d, start_row=my_group_L_start, start_col=1, &
                               n_rows=my_group_L_size, n_cols=dimen_RI)

      ! Get process, from which we get and to whom we send trace2d

      proc_send_M = para_env%mepos+para_env%num_pe/ngroup

      IF (proc_send_M >= para_env%num_pe) THEN
         proc_send_M = proc_send_M-para_env%num_pe
      END IF

      proc_rec_M = para_env%mepos-para_env%num_pe/ngroup

      IF (proc_rec_M < 0) THEN
         proc_rec_M = proc_rec_M+para_env%num_pe
      END IF

      ALLOCATE (my_group_K_send(2))
      ALLOCATE (my_group_K_rec(2))
      my_group_K_send(1) = my_group_L_start
      my_group_K_send(2) = my_group_L_size

      ALLOCATE (trace_2d_send(my_group_L_size, dimen_RI))

      trace_2d_send = trace_2d

      DO igroup = 1, ngroup

         ! send and receive my_group_K_start and my_group_K_size
         CALL mp_sendrecv(my_group_K_send, proc_send_M, my_group_K_rec, proc_rec_M, para_env%group)

         my_group_K_start = my_group_K_rec(1)
         my_group_K_size = my_group_K_rec(2)

         ALLOCATE (trace_2d_rec(my_group_K_size, dimen_RI))

         CALL mp_sendrecv(trace_2d_send, proc_send_M, trace_2d_rec, proc_rec_M, para_env%group)

         ! send and receive trace_2d
         CALL cp_fm_set_submatrix(fm=fm_mat_Q, new_values=trace_2d_rec, start_row=my_group_K_start, start_col=1, &
                                  n_rows=my_group_K_size, n_cols=dimen_RI)

         DEALLOCATE (trace_2d_send)

         ALLOCATE (trace_2d_send(my_group_K_size, dimen_RI))

         trace_2d_send = trace_2d_rec

         DEALLOCATE (trace_2d_rec)

         my_group_K_send(1) = my_group_K_start
         my_group_K_send(2) = my_group_K_size

      END DO

      DEALLOCATE (trace_2d_send)

      DEALLOCATE (my_group_K_send, my_group_K_rec)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_with_2d_array

! **************************************************************************************************
!> \brief ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param tau_tj ...
!> \param num_integ_points ...
!> \param nmo ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param mat_dm_occ_global ...
!> \param mat_dm_virt_global ...
!> \param matrix_s ...
!> \param Eigenval ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param memory_info ...
!> \param unit_nr ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE precompute_fm_scaled_dm_tau(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, nmo, &
                                          fm_mo_coeff_occ, fm_mo_coeff_virt, mat_dm_occ_global, mat_dm_virt_global, &
                                          matrix_s, Eigenval, e_fermi, eps_filter, memory_info, unit_nr, para_env)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj
      INTEGER                                            :: num_integ_points, nmo
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_dm_occ_global, mat_dm_virt_global, &
                                                            matrix_s
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      REAL(KIND=dp)                                      :: e_fermi, eps_filter
      LOGICAL                                            :: memory_info
      INTEGER                                            :: unit_nr
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precompute_fm_scaled_dm_tau', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, jjB, jquad, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: first_cycle_im_time
      REAL(KIND=dp)                                      :: tau, tau_old

      CALL timeset(routineN, handle)

      ! has to be set, otherwise compiler complains
      tau = 0.0_dp

      first_cycle_im_time = .TRUE.

      ! iteration over time points
      DO jquad = 1, num_integ_points

         tau_old = tau

         tau = tau_tj(jquad)

         ! get info of fm_mo_coeff_occ
         CALL cp_fm_get_info(matrix=fm_mo_coeff_occ, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         ! Multiply the occupied and the virtual MO coefficients with the factor exp((-e_i-e_F)*tau/2).
         ! Then, we simply get the sum over all occ states and virt. states by a simple matrix-matrix
         ! multiplication.

         ! first, the occ
         IF (first_cycle_im_time) THEN

            DO jjB = 1, nrow_local
               DO iiB = 1, ncol_local
                  i_global = col_indices(iiB)

                  fm_mo_coeff_occ%local_data(jjB, iiB) = &
                     fm_mo_coeff_occ%local_data(jjB, iiB)*EXP(tau*0.5_dp*(Eigenval(i_global)-e_fermi))
               END DO
            END DO

         ELSE

            DO jjB = 1, nrow_local
               DO iiB = 1, ncol_local
                  i_global = col_indices(iiB)

                  fm_mo_coeff_occ%local_data(jjB, iiB) = &
                     fm_mo_coeff_occ%local_data(jjB, iiB)*EXP((tau-tau_old)*0.5_dp*(Eigenval(i_global)-e_fermi))
               END DO
            END DO

         END IF

         ! get info of fm_mo_coeff_virt
         CALL cp_fm_get_info(matrix=fm_mo_coeff_virt, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         ! the same for virt
         IF (first_cycle_im_time) THEN

            DO jjB = 1, nrow_local
               DO iiB = 1, ncol_local
                  i_global = col_indices(iiB)

                  fm_mo_coeff_virt%local_data(jjB, iiB) = &
                     fm_mo_coeff_virt%local_data(jjB, iiB)*EXP(-tau*0.5_dp*(Eigenval(i_global)-e_fermi))
               END DO
            END DO

         ELSE

            DO jjB = 1, nrow_local
               DO iiB = 1, ncol_local
                  i_global = col_indices(iiB)

                  fm_mo_coeff_virt%local_data(jjB, iiB) = &
                     fm_mo_coeff_virt%local_data(jjB, iiB)*EXP(-(tau-tau_old)*0.5_dp*(Eigenval(i_global)-e_fermi))
               END DO
            END DO

         END IF

         CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mo_coeff_occ, matrix_b=fm_mo_coeff_occ, beta=0.0_dp, &
                      matrix_c=fm_scaled_dm_occ_tau(jquad)%matrix)

         CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mo_coeff_virt, matrix_b=fm_mo_coeff_virt, beta=0.0_dp, &
                      matrix_c=fm_scaled_dm_virt_tau(jquad)%matrix)

         first_cycle_im_time = .FALSE.

      END DO ! tau points

      ! transfer fm density matrices to dbcsr matrix
      NULLIFY (mat_dm_occ_global)
      CALL cp_dbcsr_allocate_matrix_set(mat_dm_occ_global, num_integ_points)
      DO jquad = 1, num_integ_points
         ALLOCATE (mat_dm_occ_global(jquad)%matrix)
         CALL cp_dbcsr_init(mat_dm_occ_global(jquad)%matrix)
         CALL cp_dbcsr_create(matrix=mat_dm_occ_global(jquad)%matrix, &
                              template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_scaled_dm_occ_tau(jquad)%matrix, &
                               mat_dm_occ_global(jquad)%matrix, &
                               keep_sparsity=.FALSE.)

         CALL cp_dbcsr_filter(mat_dm_occ_global(jquad)%matrix, eps_filter)

      END DO

      NULLIFY (mat_dm_virt_global)
      CALL cp_dbcsr_allocate_matrix_set(mat_dm_virt_global, num_integ_points)
      DO jquad = 1, num_integ_points
         ALLOCATE (mat_dm_virt_global(jquad)%matrix)
         CALL cp_dbcsr_init(mat_dm_virt_global(jquad)%matrix)
         CALL cp_dbcsr_create(matrix=mat_dm_virt_global(jquad)%matrix, &
                              template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_scaled_dm_virt_tau(jquad)%matrix, &
                               mat_dm_virt_global(jquad)%matrix, &
                               keep_sparsity=.FALSE.)

         CALL cp_dbcsr_filter(mat_dm_virt_global(jquad)%matrix, eps_filter)

      END DO

      IF (memory_info) THEN

         DO jquad = 1, num_integ_points
            CALL print_occupation_2c(mat_dm_occ_global(jquad)%matrix, unit_nr, "D_occ(it)", para_env, one_number=.TRUE.)
         END DO

         DO jquad = 1, num_integ_points
            CALL print_occupation_2c(mat_dm_virt_global(jquad)%matrix, unit_nr, "D_virt(it)", para_env, one_number=.TRUE.)
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_M_muP ...
!> \param mat_M_mu_Pnu ...
!> \param para_env ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param mepos_P_from_RI_row ...
!> \param n_group_col ...
!> \param cut_RI ...
!> \param eps_comm ...
!> \param starts_array_prim_row ...
!> \param ends_array_prim_row ...
!> \param starts_array_prim_col ...
!> \param ends_array_prim_col ...
!> \param offset_combi_block ...
!> \param starts_array_prim_fullcol ...
!> \param ends_array_prim_fullcol ...
!> \param starts_array_prim_fullrow ...
!> \param ends_array_prim_fullrow ...
!> \param group_size_P ...
!> \param dimen_RI ...
!> \param i_mem ...
!> \param j_mem ...
!> \param buffer_mat_M ...
!> \param eps_filter_im_time ...
!> \param do_occ ...
!> \param do_virt ...
! **************************************************************************************************
   SUBROUTINE fill_mat_M_muP_from_mat_M_munu_P(mat_M_muP, mat_M_mu_Pnu, para_env, my_group_L_starts_im_time, &
                                               my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                               n_group_col, cut_RI, eps_comm, &
                                               starts_array_prim_row, ends_array_prim_row, &
                                               starts_array_prim_col, ends_array_prim_col, &
                                               offset_combi_block, starts_array_prim_fullcol, &
                                               ends_array_prim_fullcol, starts_array_prim_fullrow, &
                                               ends_array_prim_fullrow, group_size_P, &
                                               dimen_RI, i_mem, j_mem, buffer_mat_M, eps_filter_im_time, &
                                               do_occ, do_virt)

      TYPE(cp_dbcsr_p_type)                              :: mat_M_muP
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_M_mu_Pnu
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_group_L_starts_im_time, &
                                                            my_group_L_sizes_im_time, &
                                                            mepos_P_from_RI_row
      INTEGER                                            :: n_group_col, cut_RI
      REAL(KIND=dp)                                      :: eps_comm
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: starts_array_prim_row, ends_array_prim_row, &
         starts_array_prim_col, ends_array_prim_col, offset_combi_block, &
         starts_array_prim_fullcol, ends_array_prim_fullcol, starts_array_prim_fullrow, &
         ends_array_prim_fullrow
      INTEGER                                            :: group_size_P, dimen_RI, i_mem, j_mem
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_mat_M
      REAL(KIND=dp)                                      :: eps_filter_im_time
      LOGICAL                                            :: do_occ, do_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_mat_M_muP_from_mat_M_munu_P', &
         routineP = moduleN//':'//routineN

      INTEGER :: block, block_counter_int, block_size, col, col_end_in_data_block, col_offset, &
         col_offset_data_block, col_offset_orig, col_offset_prim, col_offset_rec_prim, &
         col_rec_prim, col_rec_prim_rel, col_size, col_size_in_data_block, col_size_orig, &
         col_size_to_send, col_start_in_data_block, color_sub_col, color_sub_P, color_sub_row, &
         end_col, end_row, handle, handle1, i_cut_RI, imepos, LLL, mepos_P, my_group_L_size, &
         my_group_L_start, n_entries_rec, nblkrows_total, num_blocks, offset, offset_rec, &
         old_block, rec_counter, row, row_end_in_data_block, row_offset, row_offset_data_block
      INTEGER :: row_offset_prim, row_offset_rec_prim, row_rec_prim, row_rec_prim_rel, row_size, &
         row_size_in_data_block, row_size_to_send, row_start_in_data_block, send_counter, &
         start_col, start_row
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, cols_to_allocate, entry_counter, &
         num_blocks_rec, num_blocks_send, num_entries_rec, num_entries_send, row_from_LLL, &
         rows_to_allocate, sizes_rec, sizes_send
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_offset, row_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL                                            :: is_large_enough, is_new_block
      REAL(KIND=dp)                                      :: checksum
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CPASSERT(do_occ .NEQV. do_virt)

      CALL timeset(routineN, handle)

      CALL cp_dbcsr_set(mat_M_muP%matrix, 0.0_dp)

      NULLIFY (data_block)

      ALLOCATE (num_entries_send(0:para_env%num_pe-1))
      num_entries_send(:) = 0

      ALLOCATE (num_blocks_send(0:para_env%num_pe-1))
      num_blocks_send(:) = 0

      ALLOCATE (row_from_LLL(dimen_RI))
      row_from_LLL = 0

      CALL cp_dbcsr_get_info(mat_M_muP%matrix, &
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

      CALL timeset("get_sizes_M", handle1)

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         CALL cp_dbcsr_iterator_start(iter, mat_M_mu_Pnu(i_cut_RI)%matrix)

         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

            CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                              row_size=row_size, col_size=col_size, &
                                              row_offset=row_offset, col_offset=col_offset)

            ! process to send to
            DO imepos = 0, para_env%num_pe-1

               color_sub_P = imepos/group_size_P
               mepos_P = MODULO(imepos, group_size_P)
               color_sub_row = color_sub_P/n_group_col
               color_sub_col = MODULO(color_sub_P, n_group_col)

               IF (do_occ) THEN

                  IF (row >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      row <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      col >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      col <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     IF (row == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_row = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_row = ends_array_prim_fullrow(color_sub_row, i_mem)
                        row_size_to_send = end_row-start_row+1

                     ELSE IF (row == starts_array_prim_row(color_sub_row, i_mem)) THEN

                        start_row = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_row = row_offset+row_size-1
                        row_size_to_send = end_row-start_row+1

                     ELSE IF (row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_row = row_offset
                        end_row = ends_array_prim_fullrow(color_sub_row, i_mem)
                        row_size_to_send = end_row-start_row+1

                     ELSE

                        row_size_to_send = row_size

                     END IF

                     IF (col == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_col = (starts_array_prim_fullcol(color_sub_col, j_mem)-1)*my_group_L_size+1
                        end_col = ends_array_prim_fullcol(color_sub_col, j_mem)*my_group_L_size
                        col_size_to_send = end_col-start_col+1

                     ELSE IF (col == starts_array_prim_col(color_sub_col, j_mem)) THEN

                        start_col = (starts_array_prim_fullcol(color_sub_col, j_mem)-1)*my_group_L_size+1
                        end_col = col_offset+col_size-1
                        col_size_to_send = end_col-start_col+1

                     ELSE IF (col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_col = col_offset
                        end_col = ends_array_prim_fullcol(color_sub_col, j_mem)*my_group_L_size
                        col_size_to_send = end_col-start_col+1

                     ELSE

                        col_size_to_send = col_size

                     END IF

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P == mepos_P_from_RI_row(row_from_LLL(LLL+my_group_L_start-1))) THEN

                           num_entries_send(imepos) = num_entries_send(imepos)+ &
                                                      row_size_to_send*col_size_to_send/my_group_L_size

                           num_blocks_send(imepos) = num_blocks_send(imepos)+1

                        END IF

                     END DO

!                 num_entries_send(imepos) = num_entries_send(imepos) + row_size_to_send*col_size_to_send
!                 num_blocks_send(imepos)  = num_blocks_send(imepos) + my_group_L_size

                  END IF

                  ! everything transposed for virtuals
               ELSE IF (do_virt) THEN

                  col_size_orig = col_size/my_group_L_size

                  !             col_offset_orig = col_offset/my_group_L_size

                  col_offset_orig = (col_offset-1)/my_group_L_size+1

                  IF (col >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      col <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      row >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      row <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     IF (col == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_col = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_col = ends_array_prim_fullrow(color_sub_row, i_mem)
                        col_size_to_send = (end_col-start_col+1)*my_group_L_size

                     ELSE IF (col == starts_array_prim_row(color_sub_row, i_mem)) THEN

                        start_col = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_col = col_offset_orig+col_size_orig-1
                        col_size_to_send = (end_col-start_col+1)*my_group_L_size

                     ELSE IF (col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_col = col_offset_orig
                        end_col = ends_array_prim_fullrow(color_sub_row, i_mem)
                        col_size_to_send = (end_col-start_col+1)*my_group_L_size

                     ELSE

                        col_size_to_send = col_size

                     END IF

                     IF (row == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_row = starts_array_prim_fullcol(color_sub_col, j_mem)
                        end_row = ends_array_prim_fullcol(color_sub_col, j_mem)
                        row_size_to_send = end_row-start_row+1

                     ELSE IF (row == starts_array_prim_col(color_sub_col, j_mem)) THEN

                        start_row = starts_array_prim_fullcol(color_sub_col, j_mem)
                        end_row = row_offset+row_size-1
                        row_size_to_send = end_row-start_row+1

                     ELSE IF (row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_row = row_offset
                        end_row = ends_array_prim_fullcol(color_sub_col, j_mem)
                        row_size_to_send = end_row-start_row+1

                     ELSE

                        row_size_to_send = row_size

                     END IF

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P == mepos_P_from_RI_row(row_from_LLL(LLL+my_group_L_start-1))) THEN

                           num_entries_send(imepos) = num_entries_send(imepos)+ &
                                                      row_size_to_send*col_size_to_send/my_group_L_size

                           num_blocks_send(imepos) = num_blocks_send(imepos)+1

                        END IF

                     END DO

!                 num_entries_send(imepos) = num_entries_send(imepos) + row_size_to_send*col_size_to_send
!                 num_blocks_send(imepos)  = num_blocks_send(imepos) + my_group_L_size

                  END IF

               END IF

            END DO

         END DO

         CALL cp_dbcsr_iterator_stop(iter)

      END DO ! i_cut_RI

      CALL timestop(handle1)

      CALL timeset("send_sizes_M", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe-1))
      ALLOCATE (buffer_send(0:para_env%num_pe-1))

      ALLOCATE (num_entries_rec(0:para_env%num_pe-1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe-1))

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe-1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe-1))

         DO imepos = 0, para_env%num_pe-1
            sizes_send(2*imepos) = num_entries_send(imepos)
            sizes_send(2*imepos+1) = num_blocks_send(imepos)
         END DO

         CALL mp_alltoall(sizes_send, sizes_rec, 2, para_env%group)

         DO imepos = 0, para_env%num_pe-1
            num_entries_rec(imepos) = sizes_rec(2*imepos)
            num_blocks_rec(imepos) = sizes_rec(2*imepos+1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send(0)
         num_blocks_rec(0) = num_blocks_send(0)

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_buffer_send_M", handle1)

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe-1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 10))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send(imepos), 10))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe-1))
      entry_counter(:) = 0

      ALLOCATE (block_counter(0:para_env%num_pe-1))
      block_counter(:) = 0

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         ! fill buffer_send
         CALL cp_dbcsr_iterator_start(iter, mat_M_mu_Pnu(i_cut_RI)%matrix)
         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

            CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                              row_size=row_size, col_size=col_size, &
                                              row_offset=row_offset, col_offset=col_offset)

            DO imepos = 0, para_env%num_pe-1

!             color_sub_row = imepos/n_group_col
!             color_sub_col = MODULO(imepos,n_group_col)

               color_sub_P = imepos/group_size_P
               mepos_P = MODULO(imepos, group_size_P)
               color_sub_row = color_sub_P/n_group_col
               color_sub_col = MODULO(color_sub_P, n_group_col)

               IF (do_occ) THEN

                  IF (row >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      row <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      col >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      col <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     col_size_orig = col_size/my_group_L_size

                     ! For terminal blocks, we have to compute the sizes
                     IF (row == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem)-row_offset+1
                        row_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem)-row_offset+1

                     ELSE IF (row == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              row .NE. ends_array_prim_row(color_sub_row, i_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem)-row_offset+1
                        row_end_in_data_block = row_size

                     ELSE IF (row .NE. starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        row_start_in_data_block = 1
                        row_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem)-row_offset+1

                     ELSE

                        row_start_in_data_block = 1
                        row_end_in_data_block = row_size
                        row_size_in_data_block = row_size

                     END IF

                     row_size_in_data_block = row_end_in_data_block-row_start_in_data_block+1

                     col_offset_orig = (col_offset-1)/my_group_L_size+1

                     ! For terminal blocks, we have to compute the sizes
                     IF (col == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem)-col_offset_orig+1
                        col_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem)-col_offset_orig+1

                     ELSE IF (col == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              col .NE. ends_array_prim_col(color_sub_col, j_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem)-col_offset_orig+1
                        col_end_in_data_block = col_size_orig

                     ELSE IF (col .NE. starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        col_start_in_data_block = 1
                        col_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem)-col_offset_orig+1

                     ELSE

                        col_start_in_data_block = 1
                        col_end_in_data_block = col_size_orig

                     END IF

                     col_size_in_data_block = col_end_in_data_block-col_start_in_data_block+1

                     block_size = row_size_in_data_block*col_size_in_data_block

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P .NE. mepos_P_from_RI_row(row_from_LLL(LLL+my_group_L_start-1))) CYCLE

                        offset = entry_counter(imepos)

                        col_offset_data_block = (LLL-1)*col_size_orig+col_start_in_data_block

                        buffer_send(imepos)%msg(offset+1:offset+block_size) = &
                           RESHAPE(data_block(row_start_in_data_block:row_end_in_data_block, &
                                              col_offset_data_block:col_offset_data_block+col_size_in_data_block-1), &
                                   (/block_size/))

                        block = block_counter(imepos)+1

                        buffer_send(imepos)%indx(block, 1) = LLL+my_group_L_start-1
                        buffer_send(imepos)%indx(block, 2) = row
                        buffer_send(imepos)%indx(block, 3) = row_offset
                        buffer_send(imepos)%indx(block, 4) = row_size
                        buffer_send(imepos)%indx(block, 5) = col
                        buffer_send(imepos)%indx(block, 6) = col_offset
                        buffer_send(imepos)%indx(block, 7) = col_size_orig
                        buffer_send(imepos)%indx(block, 8) = offset
                        buffer_send(imepos)%indx(block, 9) = block_size

                        entry_counter(imepos) = entry_counter(imepos)+block_size

                        block_counter(imepos) = block_counter(imepos)+1

                     END DO

                  END IF

               ELSE IF (do_virt) THEN

                  IF (col >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      col <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      row >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      row <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     col_size_orig = col_size/my_group_L_size

                     col_offset_orig = (col_offset-1)/my_group_L_size+1

                     ! For terminal blocks, we have to compute the sizes
                     IF (col == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem)-col_offset_orig+1
                        col_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem)-col_offset_orig+1

                     ELSE IF (col == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              col .NE. ends_array_prim_row(color_sub_row, i_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem)-col_offset_orig+1
                        col_end_in_data_block = col_size_orig

                     ELSE IF (col .NE. starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        col_start_in_data_block = 1
                        col_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem)-col_offset_orig+1

                     ELSE

                        col_start_in_data_block = 1
                        col_end_in_data_block = col_size_orig

                     END IF

                     col_size_in_data_block = col_end_in_data_block-col_start_in_data_block+1

                     ! For terminal blocks, we have to compute the sizes
                     IF (row == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem)-row_offset+1
                        row_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem)-row_offset+1

                     ELSE IF (row == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              row .NE. ends_array_prim_col(color_sub_col, j_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem)-row_offset+1
                        row_end_in_data_block = row_size

                     ELSE IF (row .NE. starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        row_start_in_data_block = 1
                        row_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem)-row_offset+1

                     ELSE

                        row_start_in_data_block = 1
                        row_end_in_data_block = row_size
                        row_size_in_data_block = row_size

                     END IF

                     row_size_in_data_block = row_end_in_data_block-row_start_in_data_block+1

                     block_size = row_size_in_data_block*col_size_in_data_block

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P .NE. mepos_P_from_RI_row(row_from_LLL(LLL+my_group_L_start-1))) CYCLE

                        offset = entry_counter(imepos)

                        col_offset_data_block = (LLL-1)*col_size_orig+col_start_in_data_block

                        buffer_send(imepos)%msg(offset+1:offset+block_size) = &
                           RESHAPE(TRANSPOSE(data_block(row_start_in_data_block:row_end_in_data_block, &
                                                        col_offset_data_block:col_offset_data_block+col_size_in_data_block-1)), &
                                   (/block_size/))

                        block = block_counter(imepos)+1

                        buffer_send(imepos)%indx(block, 1) = LLL+my_group_L_start-1
                        buffer_send(imepos)%indx(block, 2) = col
                        buffer_send(imepos)%indx(block, 3) = col_offset
                        buffer_send(imepos)%indx(block, 4) = col_size_orig
                        buffer_send(imepos)%indx(block, 5) = row
                        buffer_send(imepos)%indx(block, 6) = row_offset
                        buffer_send(imepos)%indx(block, 7) = row_size
                        buffer_send(imepos)%indx(block, 8) = offset
                        buffer_send(imepos)%indx(block, 9) = block_size

                        entry_counter(imepos) = entry_counter(imepos)+block_size

                        block_counter(imepos) = block_counter(imepos)+1

                     END DO

                  END IF

               END IF

            END DO

         END DO

         CALL cp_dbcsr_iterator_stop(iter)

         CALL cp_dbcsr_set(mat_M_mu_Pnu(i_cut_RI)%matrix, 0.0_dp)
         CALL cp_dbcsr_filter(mat_M_mu_Pnu(i_cut_RI)%matrix, 1.0_dp)

      END DO ! i_cut_RI

      CALL timestop(handle1)

      CALL timeset("comm_data_M", handle1)

      DEALLOCATE (entry_counter, block_counter)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

!old
!         ! send around the data and indices
!         DO imepos = 0, para_env%num_pe-1
!            IF (num_entries_send(imepos) > 0) THEN
!               send_counter = send_counter+1
!               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
!            END IF
!            IF (num_entries_rec(imepos) > 0) THEN
!               rec_counter = rec_counter+1
!               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
!            END IF
!            IF (num_entries_send(imepos) > 0) THEN
!               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
!            END IF
!            IF (num_entries_rec(imepos) > 0) THEN
!               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
!            END IF
!         END DO

         DO imepos = 0, para_env%num_pe-1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter+1
               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
            END IF
            IF (num_entries_rec(imepos) > 0) THEN
               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO imepos = 0, para_env%num_pe-1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter+1
               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
            END IF
            IF (num_entries_send(imepos) > 0) THEN
               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         CALL mp_waitall(req_array(1:send_counter, 1:2))
         CALL mp_waitall(req_array(1:rec_counter, 3:4))

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      CALL timestop(handle1)

      CALL timeset("reserve_blocks_M", handle1)

      num_blocks = 0

      ! get the number of blocks, which have to be allocated
      DO imepos = 0, para_env%num_pe-1

         num_blocks = num_blocks+num_blocks_rec(imepos)

      END DO

      ALLOCATE (rows_to_allocate(num_blocks))
      rows_to_allocate = 0
      ALLOCATE (cols_to_allocate(num_blocks))
      cols_to_allocate = 1

      block_counter_int = 0

      DO imepos = 0, para_env%num_pe-1

         DO block = 1, num_blocks_rec(imepos)

            is_new_block = .TRUE.

            ! check whether block is already there
            DO old_block = 1, block_counter_int

               IF (row_from_LLL(buffer_rec(imepos)%indx(block, 1)) == rows_to_allocate(old_block)) THEN

                  is_new_block = .FALSE.

               END IF

            END DO

            ! check whether block has non-vanishing entries
            offset_rec = buffer_rec(imepos)%indx(block, 8)
            block_size = buffer_rec(imepos)%indx(block, 9)

            checksum = SUM(ABS(buffer_rec(imepos)%msg(offset_rec+1:offset_rec+block_size)))

            is_large_enough = (checksum > eps_comm)

            IF (is_new_block .AND. is_large_enough) THEN

               block_counter_int = block_counter_int+1

               rows_to_allocate(block_counter_int) = row_from_LLL(buffer_rec(imepos)%indx(block, 1))

            END IF

         END DO

      END DO

      CALL cp_dbcsr_set(mat_M_muP%matrix, 0.0_dp)

      CALL cp_dbcsr_filter(mat_M_muP%matrix, 1.0_dp)

      CALL cp_dbcsr_reserve_blocks(mat_M_muP%matrix, rows=rows_to_allocate(1:block_counter_int), &
                                   cols=cols_to_allocate(1:block_counter_int))

      CALL cp_dbcsr_finalize(mat_M_muP%matrix)

      CALL timestop(handle1)

      CALL cp_dbcsr_set(mat_M_muP%matrix, 0.0_dp)

      n_entries_rec = SUM(num_entries_rec)

      ! Fill the dbcsr matrix
      CALL timeset("fill_dbcsr_mat_M", handle1)

      color_sub_P = para_env%mepos/group_size_P

      color_sub_row = color_sub_P/n_group_col
      color_sub_col = MODULO(para_env%mepos, n_group_col)

      row_offset_prim = starts_array_prim_row(color_sub_row, i_mem)
      col_offset_prim = starts_array_prim_col(color_sub_col, j_mem)

      CALL cp_dbcsr_iterator_start(iter, mat_M_muP%matrix)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))

         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

         buffer_mat_M(:, :) = 0.0_dp

         DO imepos = 0, para_env%num_pe-1

            DO block = 1, num_blocks_rec(imepos)

               LLL = buffer_rec(imepos)%indx(block, 1)

               IF (LLL >= row_offset .AND. LLL < row_offset+row_size) THEN

                  row_rec_prim = buffer_rec(imepos)%indx(block, 2)
                  row_offset_rec_prim = buffer_rec(imepos)%indx(block, 3)
                  col_rec_prim = buffer_rec(imepos)%indx(block, 5)
                  col_offset_rec_prim = buffer_rec(imepos)%indx(block, 6)
                  offset_rec = buffer_rec(imepos)%indx(block, 8)
                  block_size = buffer_rec(imepos)%indx(block, 9)

                  row_rec_prim_rel = row_rec_prim-row_offset_prim+1
                  col_rec_prim_rel = col_rec_prim-col_offset_prim+1

                  row_offset_data_block = LLL-row_offset

                  col_offset_data_block = offset_combi_block(row_rec_prim, col_rec_prim)

!                 data_block(row_offset_data_block+1,col_offset_data_block+1:col_offset_data_block+block_size) = &
!                   buffer_rec(imepos)%msg(offset_rec+1:offset_rec+block_size)

                  buffer_mat_M(row_offset_data_block+1, col_offset_data_block+1:col_offset_data_block+block_size) = &
                     buffer_rec(imepos)%msg(offset_rec+1:offset_rec+block_size)

               END IF

            END DO

         END DO

         data_block(1:row_size, 1:col_size) = buffer_mat_M(1:row_size, 1:col_size)

      END DO

      CALL cp_dbcsr_iterator_stop(iter)

      CALL cp_dbcsr_filter(mat_M_muP%matrix, eps_filter_im_time)

      DO imepos = 0, para_env%num_pe-1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
      END DO

      DEALLOCATE (req_array)

      DEALLOCATE (buffer_rec, buffer_send)

      DEALLOCATE (num_entries_send, num_entries_rec, num_blocks_send, num_blocks_rec)

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int_cut ...
!> \param mat_3c_overl_int ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param starts_array_cm ...
!> \param ends_array_cm ...
!> \param my_group_L_sizes_im_time ...
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE setup_mat_for_mem_cut_3c(mat_3c_overl_int_cut, mat_3c_overl_int, cut_memory, cut_RI, &
                                       starts_array_cm, ends_array_cm, &
                                       my_group_L_sizes_im_time, eps_filter)

      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: mat_3c_overl_int_cut
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_3c_overl_int
      INTEGER                                            :: cut_memory, cut_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: starts_array_cm, ends_array_cm, &
                                                            my_group_L_sizes_im_time
      REAL(KIND=dp)                                      :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_mat_for_mem_cut_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, col, col_end_in_data_block, col_offset, col_size, col_start_in_data_block, &
         handle, i_cut_RI, i_mem, my_group_L_size, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      NULLIFY (mat_3c_overl_int_cut)
      CALL cp_dbcsr_allocate_matrix_set(mat_3c_overl_int_cut, cut_RI, cut_memory)

      DO i_mem = 1, cut_memory
         DO i_cut_RI = 1, cut_RI
            ALLOCATE (mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix)
            CALL cp_dbcsr_init(mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix)
            CALL cp_dbcsr_create(matrix=mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix, &
                                 template=mat_3c_overl_int(i_cut_RI)%matrix)

            CALL cp_dbcsr_copy(mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix, mat_3c_overl_int(i_cut_RI)%matrix)

            my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

            CALL cp_dbcsr_iterator_start(iter, mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix)
            DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
               CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                                 col_size=col_size, col_offset=col_offset)

               ! set the block to zero if it is outside the range for cutting the mu-sigma combined index
               IF (col_offset > ends_array_cm(i_mem)*my_group_L_size .OR. &
                   col_offset+col_size-1 < (starts_array_cm(i_mem)-1)*my_group_L_size+1) THEN

                  data_block = 0.0_dp

               END IF

               IF (col_offset >= ends_array_cm(i_mem)*my_group_L_size .AND. &
                   col_offset+col_size-1 < ends_array_cm(i_mem)*my_group_L_size) THEN

                  col_end_in_data_block = ends_array_cm(i_mem)*my_group_L_size-col_offset+1

                  data_block(:, col_end_in_data_block+1:col_size) = 0.0_dp

               END IF

               IF (col_offset > (starts_array_cm(i_mem)-1)*my_group_L_size+1 .AND. &
                   col_offset+col_size-1 <= (starts_array_cm(i_mem)-1)*my_group_L_size+1) THEN

                  col_start_in_data_block = (starts_array_cm(i_mem)-1)*my_group_L_size+1

                  data_block(:, 1:col_start_in_data_block-1) = 0.0_dp

               END IF

            END DO

            CALL cp_dbcsr_iterator_stop(iter)

            ! remove the zeroed blocks
            CALL cp_dbcsr_filter(mat_3c_overl_int_cut(i_cut_RI, i_mem)%matrix, eps_filter)

         END DO
      END DO

      CALL cp_dbcsr_deallocate_matrix_set(mat_3c_overl_int)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param num_entries_rec ...
!> \param num_entries_send ...
!> \param buffer_rec ...
!> \param buffer_send ...
!> \param req_array ...
! **************************************************************************************************
   SUBROUTINE communicate_buffer(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, req_array)

      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_entries_rec, num_entries_send
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array

      CHARACTER(LEN=*), PARAMETER :: routineN = 'communicate_buffer', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imepos, rec_counter, send_counter

      CALL timeset(routineN, handle)

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imepos = 0, para_env%num_pe-1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter+1
               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
            END IF
            IF (num_entries_rec(imepos) > 0) THEN
               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO imepos = 0, para_env%num_pe-1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter+1
               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
            END IF
            IF (num_entries_send(imepos) > 0) THEN
               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         CALL mp_waitall(req_array(1:send_counter, 1:2))
         CALL mp_waitall(req_array(1:rec_counter, 3:4))

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_occ_global ...
!> \param mat_dm_virt_global ...
! **************************************************************************************************
   SUBROUTINE clean_up(mat_dm_occ_global, mat_dm_virt_global)
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mat_dm_occ_global, mat_dm_virt_global

      CALL cp_dbcsr_deallocate_matrix_set(mat_dm_occ_global)
      CALL cp_dbcsr_deallocate_matrix_set(mat_dm_virt_global)

   END SUBROUTINE

END MODULE rpa_im_time
