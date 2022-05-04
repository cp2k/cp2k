!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate RI-GPW-MP2 energy using pw
!> \par History
!>      06.2012 created [Mauro Del Ben]
!>      03.2019 Refactored from mp2_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_ri_gpw
   USE ISO_C_BINDING,                   ONLY: C_NULL_PTR,&
                                              c_associated
   USE cp_log_handling,                 ONLY: cp_to_string
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_split
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE group_dist_types,                ONLY: get_group_dist,&
                                              group_dist_d1_type,&
                                              maxsize,&
                                              release_group_dist
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: compare_potential_types
   USE machine,                         ONLY: m_flush,&
                                              m_memory,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_min,&
                                              mp_recv,&
                                              mp_send,&
                                              mp_sendrecv,&
                                              mp_sum
   USE mp2_ri_grad_util,                ONLY: complete_gamma
   USE mp2_types,                       ONLY: dgemm_counter_init,&
                                              dgemm_counter_start,&
                                              dgemm_counter_stop,&
                                              dgemm_counter_type,&
                                              dgemm_counter_write,&
                                              mp2_type,&
                                              three_dim_real_array,&
                                              two_dim_int_array
   USE offload_gemm,                    ONLY: OFFLOAD_GEMM_PU_GPU,&
                                              offload_dgemm,&
                                              offload_gemm_create,&
                                              offload_gemm_destroy,&
                                              offload_gemm_set_op_threshold_gpu

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_gpw'

   PUBLIC :: mp2_ri_gpw_compute_en

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
!> \param Emp2_S ...
!> \param Emp2_T ...
!> \param BIb_C ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param gd_array ...
!> \param gd_B_virtual ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param calc_forces ...
!> \param calc_ex ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_gpw_compute_en(Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                                    gd_array, gd_B_virtual, &
                                    Eigenval, nmo, homo, dimen_RI, unit_nr, calc_forces, calc_ex)
      REAL(KIND=dp), INTENT(INOUT)                       :: Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T
      TYPE(three_dim_real_array), DIMENSION(:), &
         INTENT(INOUT)                                   :: BIb_C
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: color_sub
      TYPE(group_dist_d1_type), INTENT(INOUT)            :: gd_array
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      INTEGER, INTENT(IN)                                :: nmo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      TYPE(group_dist_d1_type), DIMENSION(SIZE(homo)), &
         INTENT(INOUT)                                   :: gd_B_virtual
      INTEGER, INTENT(IN)                                :: dimen_RI, unit_nr
      LOGICAL, INTENT(IN)                                :: calc_forces, calc_ex

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_gpw_compute_en'

      INTEGER :: a, a_global, b, b_global, block_size, decil, end_point, handle, handle2, handle3, &
         iiB, ij_counter, ij_counter_send, ij_index, integ_group_size, ispin, jjB, jspin, kspin, &
         max_ij_pairs, my_block_size, my_group_L_end, my_group_L_size, my_group_L_size_orig, &
         my_group_L_start, my_i, my_ij_pairs, my_j, my_new_group_L_size, ngroup, nspins, &
         num_IJ_blocks, num_integ_group, proc_receive, proc_send, proc_shift, rec_B_size, &
         rec_B_virtual_end, rec_B_virtual_start, rec_L_size, send_B_size, send_B_virtual_end, &
         send_B_virtual_start, send_block_size, send_i, send_ij_index, send_j, start_point
      INTEGER :: sub_P_color, sub_sub_color, tag, total_ij_pairs
      INTEGER, ALLOCATABLE, DIMENSION(:) :: integ_group_pos2color_sub, my_B_size, &
         my_B_virtual_end, my_B_virtual_start, num_ij_pairs, proc_map, proc_map_rep, &
         sizes_array_orig, sub_proc_map, virtual
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: ranges_info_array
      LOGICAL                                            :: my_alpha_alpha_case, my_alpha_beta_case, &
                                                            my_beta_beta_case, my_open_shell_SS
      REAL(KIND=dp)                                      :: amp_fac, my_Emp2_Cou, my_Emp2_EX, &
                                                            sym_fac, t_new, t_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         TARGET                                          :: external_ab, external_i_aL, local_ab, &
                                                            local_ba, t_ab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         TARGET                                          :: BI_C_rec, local_i_aL, local_j_aL, &
                                                            Y_i_aP, Y_j_aP
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange, para_env_P, &
                                                            para_env_rep
      TYPE(dgemm_counter_type)                           :: dgemm_counter
      TYPE(three_dim_real_array), ALLOCATABLE, &
         DIMENSION(:)                                    :: B_ia_Q

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)

      ALLOCATE (virtual(nspins))
      virtual(:) = nmo - homo(:)

      ALLOCATE (my_B_size(nspins), my_B_virtual_start(nspins), my_B_virtual_end(nspins))
      DO ispin = 1, nspins
         CALL get_group_dist(gd_B_virtual(ispin), para_env_sub%mepos, &
                             my_B_virtual_start(ispin), my_B_virtual_end(ispin), my_B_size(ispin))
      END DO

      CALL get_group_dist(gd_array, color_sub, my_group_L_start, my_group_L_end, my_group_L_size)

      CALL dgemm_counter_init(dgemm_counter, unit_nr, mp2_env%ri_mp2%print_dgemm_info)

      ! offload_gemm_ctx has a very footprint the first time this routine is
      ! called.
      IF (.NOT. c_associated(mp2_env%offload_gemm_ctx)) THEN
         CALL offload_gemm_create(mp2_env%offload_gemm_ctx, OFFLOAD_GEMM_PU_GPU)
         CALL offload_gemm_set_op_threshold_gpu(mp2_env%offload_gemm_ctx, 128*128*128*2)
      END IF

      ! We cannot fix the tag because of the recv routine
      tag = 42

      DO jspin = 1, nspins
      DO ispin = jspin, 1, -1

         my_open_shell_SS = (nspins == 2) .AND. (ispin == jspin)

         ! t_ab = amp_fac*(:,a|:,b)-(:,b|:,a)
         ! If we calculate the gradient we need to distinguish
         ! between alpha-alpha and beta-beta cases for UMP2

         my_alpha_alpha_case = .FALSE.
         my_beta_beta_case = .FALSE.
         my_alpha_beta_case = .FALSE.
         IF (ispin /= jspin) THEN
            my_alpha_beta_case = .TRUE.
         ELSE IF (my_open_shell_SS) THEN
            IF (ispin == 1) my_alpha_alpha_case = .TRUE.
            IF (ispin == 2) my_beta_beta_case = .TRUE.
         END IF

         amp_fac = mp2_env%scale_S + mp2_env%scale_T
         IF (my_alpha_beta_case .OR. my_open_shell_SS) amp_fac = mp2_env%scale_T

         CALL mp2_ri_get_integ_group_size( &
            mp2_env, para_env, para_env_sub, gd_array, gd_B_virtual(ispin), &
            homo(ispin), dimen_RI, unit_nr, &
            integ_group_size, ngroup, &
            num_integ_group, virtual(ispin), my_alpha_beta_case, &
            my_open_shell_SS, calc_forces)

         CALL mp2_ri_get_block_size( &
            mp2_env, para_env, para_env_sub, gd_array, gd_B_virtual(ispin), &
            homo(ispin), dimen_RI, unit_nr, block_size, &
            ngroup, num_IJ_blocks, &
            virtual(ispin), my_alpha_beta_case, &
            my_open_shell_SS, calc_forces)

         ! now create a group that contains all the proc that have the same virtual starting point
         ! in the integ group
         ! sub_sub_color=para_env_sub%mepos
         CALL mp2_ri_create_group( &
            BIb_C(ispin)%array, para_env, para_env_sub, homo(ispin), color_sub, &
            gd_array%sizes, calc_forces, &
            integ_group_size, my_B_size(ispin), iiB, my_group_L_end, &
            my_group_L_size, my_group_L_size_orig, my_group_L_start, my_new_group_L_size, &
            sub_sub_color, integ_group_pos2color_sub, proc_map, proc_map_rep, sizes_array_orig, &
            sub_proc_map, ranges_info_array, para_env_exchange, para_env_rep, num_integ_group)

         ! *****************************************************************
         ! **********  REPLICATION-BLOCKED COMMUNICATION SCHEME  ***********
         ! *****************************************************************
         ! introduce block size, the number of occupied orbitals has to be a
         ! multiple of the block size

         ! Calculate the maximum number of ij pairs that have to be computed
         ! among groups
         CALL mp2_ri_communication(my_alpha_beta_case, total_ij_pairs, homo(ispin), homo(jspin), num_IJ_blocks, &
                                   block_size, ngroup, ij_map, color_sub, my_ij_pairs, my_open_shell_SS, unit_nr)

         ALLOCATE (num_ij_pairs(0:para_env_exchange%num_pe - 1))
         CALL mp_allgather(my_ij_pairs, num_ij_pairs, para_env_exchange%group)

         max_ij_pairs = MAXVAL(num_ij_pairs)

         ! start real stuff
         CALL mp2_ri_allocate(local_ab, t_ab, mp2_env, homo, virtual, dimen_RI, my_B_size, &
                              block_size, my_group_L_size, local_i_aL, &
                              local_j_aL, calc_forces, Y_i_aP, Y_j_aP, &
                              ispin, jspin, local_ba)

         IF (unit_nr > 0) THEN
            IF (nspins == 1) THEN
               WRITE (unit_nr, *) "Start loop run"
            ELSE IF (ispin == 1 .AND. jspin == 1) THEN
               WRITE (unit_nr, *) "Start loop run alpha-alpha"
            ELSE IF (ispin == 1 .AND. jspin == 2) THEN
               WRITE (unit_nr, *) "Start loop run alpha-beta"
            ELSE IF (ispin == 2 .AND. jspin == 2) THEN
               WRITE (unit_nr, *) "Start loop run beta-beta"
            END IF
            CALL m_flush(unit_nr)
            t_start = m_walltime()
         END IF

         CALL timeset(routineN//"_RI_loop", handle2)
         my_Emp2_Cou = 0.0_dp
         my_Emp2_EX = 0.0_dp
         DO ij_index = 1, max_ij_pairs

            IF (unit_nr > 0) THEN
               decil = ij_index*10/max_ij_pairs
               IF (decil /= (ij_index - 1)*10/max_ij_pairs) THEN
                  t_new = m_walltime()
                  t_new = (t_new - t_start)/60.0_dp*(max_ij_pairs - ij_index + 1)/MAX(ij_index - 1, 1)
                  WRITE (unit_nr, FMT="(T3,A)") "Percentage of finished loop: "// &
                     cp_to_string(decil*10)//". Minutes left: "//cp_to_string(t_new)
                  CALL m_flush(unit_nr)
               END IF
            END IF

            IF (calc_forces) THEN
               Y_i_aP = 0.0_dp
               Y_j_aP = 0.0_dp
            END IF

            IF (ij_index <= my_ij_pairs) THEN
               ! We have work to do
               ij_counter = (ij_index - MIN(1, color_sub))*ngroup + color_sub
               my_i = ij_map(1, ij_counter)
               my_j = ij_map(2, ij_counter)
               my_block_size = ij_map(3, ij_counter)

               local_i_aL = 0.0_dp
               CALL fill_local_i_aL(local_i_aL(:, :, 1:my_block_size), ranges_info_array(:, :, para_env_exchange%mepos), &
                                    BIb_C(ispin)%array(:, :, my_i:my_i + my_block_size - 1))

               local_j_aL = 0.0_dp
               CALL fill_local_i_aL(local_j_aL(:, :, 1:my_block_size), ranges_info_array(:, :, para_env_exchange%mepos), &
                                    BIb_C(jspin)%array(:, :, my_j:my_j + my_block_size - 1))

               ! collect data from other proc
               CALL timeset(routineN//"_comm", handle3)
               DO proc_shift = 1, para_env_exchange%num_pe - 1
                  proc_send = proc_map(para_env_exchange%mepos + proc_shift)
                  proc_receive = proc_map(para_env_exchange%mepos - proc_shift)

                  send_ij_index = num_ij_pairs(proc_send)

                  CALL get_group_dist(gd_array, proc_receive, sizes=rec_L_size)
                  ALLOCATE (BI_C_rec(rec_L_size, MAXVAL(my_B_size), my_block_size))

                  IF (ij_index <= send_ij_index) THEN
                     ij_counter_send = (ij_index - MIN(1, integ_group_pos2color_sub(proc_send)))*ngroup + &
                                       integ_group_pos2color_sub(proc_send)
                     send_i = ij_map(1, ij_counter_send)
                     send_j = ij_map(2, ij_counter_send)
                     send_block_size = ij_map(3, ij_counter_send)

                     ! occupied i
                     BI_C_rec = 0.0_dp
                     CALL mp_sendrecv(BIb_C(ispin)%array(:, :, send_i:send_i + send_block_size - 1), &
                                      proc_send, BI_C_rec, proc_receive, &
                                      para_env_exchange%group, tag)

                     CALL fill_local_i_aL(local_i_aL(:, :, 1:my_block_size), ranges_info_array(:, :, proc_receive), &
                                          BI_C_rec(:, 1:my_B_size(ispin), :))

                     ! occupied j
                     BI_C_rec = 0.0_dp
                     CALL mp_sendrecv(BIb_C(jspin)%array(:, :, send_j:send_j + send_block_size - 1), &
                                      proc_send, BI_C_rec, proc_receive, &
                                      para_env_exchange%group, tag)

                     CALL fill_local_i_aL(local_j_aL(:, :, 1:my_block_size), ranges_info_array(:, :, proc_receive), &
                                          BI_C_rec(:, 1:my_B_size(jspin), :))

                  ELSE
                     ! we send nothing while we know that we have to receive something

                     ! occupied i
                     BI_C_rec = 0.0_dp
                     CALL mp_recv(BI_C_rec, proc_receive, tag, para_env_exchange%group)

                     CALL fill_local_i_aL(local_i_aL(:, :, 1:my_block_size), ranges_info_array(:, :, proc_receive), &
                                          BI_C_rec(:, 1:my_B_size(ispin), :))

                     ! occupied j
                     BI_C_rec = 0.0_dp
                     CALL mp_recv(BI_C_rec, proc_receive, tag, para_env_exchange%group)

                     CALL fill_local_i_aL(local_j_aL(:, :, 1:my_block_size), ranges_info_array(:, :, proc_receive), &
                                          BI_C_rec(:, 1:my_B_size(jspin), :))

                  END IF

                  DEALLOCATE (BI_C_rec)

               END DO
               CALL timestop(handle3)

               ! loop over the block elements
               DO iiB = 1, my_block_size
                  DO jjB = 1, my_block_size
                     CALL timeset(routineN//"_expansion", handle3)
                     ASSOCIATE (my_local_i_aL => local_i_aL(:, :, iiB), my_local_j_aL => local_j_aL(:, :, jjB))

                        ! calculate the integrals (ia|jb) strating from my local data ...
                        local_ab = 0.0_dp
                        IF ((my_alpha_beta_case) .AND. (calc_forces)) THEN
                           local_ba = 0.0_dp
                        END IF
                        CALL dgemm_counter_start(dgemm_counter)
                        CALL offload_dgemm('T', 'N', my_B_size(ispin), my_B_size(jspin), dimen_RI, 1.0_dp, &
                                           my_local_i_aL, dimen_RI, my_local_j_aL, dimen_RI, &
                                         0.0_dp, local_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :), my_B_size(ispin), &
                                           mp2_env%offload_gemm_ctx)
                        ! Additional integrals only for alpha_beta case and forces
                        IF (my_alpha_beta_case .AND. calc_forces) THEN
                           local_ba(my_B_virtual_start(jspin):my_B_virtual_end(jspin), :) = &
                              TRANSPOSE(local_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :))
                        END IF
                        ! ... and from the other of my subgroup
                        DO proc_shift = 1, para_env_sub%num_pe - 1
                           proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
                           proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

                           CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, &
                                               rec_B_virtual_end, rec_B_size)

                           ALLOCATE (external_i_aL(dimen_RI, rec_B_size))
                           external_i_aL = 0.0_dp

                           CALL mp_sendrecv(my_local_i_aL, proc_send, &
                                            external_i_aL, proc_receive, &
                                            para_env_sub%group, tag)

                           CALL offload_dgemm('T', 'N', rec_B_size, my_B_size(jspin), dimen_RI, 1.0_dp, &
                                              external_i_aL, dimen_RI, my_local_j_aL, dimen_RI, &
                                          0.0_dp, local_ab(rec_B_virtual_start:rec_B_virtual_end, 1:my_B_size(jspin)), rec_B_size, &
                                              mp2_env%offload_gemm_ctx)

                           DEALLOCATE (external_i_aL)
                           ! Additional integrals only for alpha_beta case and forces
                           IF (my_alpha_beta_case .AND. calc_forces) THEN

                              CALL get_group_dist(gd_B_virtual(jspin), proc_receive, rec_B_virtual_start, &
                                                  rec_B_virtual_end, rec_B_size)

                              ALLOCATE (external_i_aL(dimen_RI, rec_B_size))
                              external_i_aL = 0.0_dp

                              CALL mp_sendrecv(my_local_j_aL, proc_send, &
                                               external_i_aL, proc_receive, &
                                               para_env_sub%group, tag)

                              CALL offload_dgemm('T', 'N', rec_B_size, my_B_size(ispin), dimen_RI, 1.0_dp, &
                                                 external_i_aL, dimen_RI, my_local_i_aL, dimen_RI, &
                                          0.0_dp, local_ba(rec_B_virtual_start:rec_B_virtual_end, 1:my_B_size(ispin)), rec_B_size, &
                                                 mp2_env%offload_gemm_ctx)
                              DEALLOCATE (external_i_aL)
                           END IF
                        END DO
                        IF (my_alpha_beta_case .AND. calc_forces) THEN
                           ! Is just an approximation, but the call does not allow it, it ought to be (virtual_i*B_size_j+virtual_j*B_size_i)*dimen_RI
                           CALL dgemm_counter_stop(dgemm_counter, virtual(ispin), my_B_size(ispin) + my_B_size(jspin), dimen_RI)
                        ELSE
                           CALL dgemm_counter_stop(dgemm_counter, virtual(ispin), my_B_size(jspin), dimen_RI)
                        END IF
                        CALL timestop(handle3)

                        !sample peak memory
                        CALL m_memory()

                        CALL timeset(routineN//"_ener", handle3)
                        ! calculate coulomb only MP2
                        sym_fac = 2.0_dp
                        IF (my_i == my_j) sym_fac = 1.0_dp
                        IF (my_alpha_beta_case) sym_fac = 0.5_dp
                        DO b = 1, my_B_size(jspin)
                           b_global = b + my_B_virtual_start(jspin) - 1
                           DO a = 1, virtual(ispin)
                              my_Emp2_Cou = my_Emp2_Cou - sym_fac*2.0_dp*local_ab(a, b)**2/ &
                                            (Eigenval(homo(ispin) + a, ispin) + Eigenval(homo(jspin) + b_global, jspin) - &
                                             Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, jspin))
                           END DO
                        END DO

                        IF (calc_ex) THEN
                           ! contract integrals with orbital energies for exchange MP2 energy
                           ! starting with local ...
                           IF (calc_forces .AND. (.NOT. my_alpha_beta_case)) t_ab = 0.0_dp
                           DO b = 1, my_B_size(ispin)
                              b_global = b + my_B_virtual_start(ispin) - 1
                              DO a = 1, my_B_size(ispin)
                                 a_global = a + my_B_virtual_start(ispin) - 1
                                 my_Emp2_Ex = my_Emp2_Ex + sym_fac*local_ab(a_global, b)*local_ab(b_global, a)/ &
                                              (Eigenval(homo(ispin) + a_global, ispin) + Eigenval(homo(ispin) + b_global, ispin) - &
                                               Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, ispin))
                                 IF (calc_forces .AND. (.NOT. my_alpha_beta_case)) THEN
                                    t_ab(a_global, b) = -(amp_fac*local_ab(a_global, b) - mp2_env%scale_T*local_ab(b_global, a))/ &
                                                        (Eigenval(homo(ispin) + a_global, ispin) + &
                                                         Eigenval(homo(ispin) + b_global, ispin) - &
                                                         Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, ispin))
                                 END IF
                              END DO
                           END DO
                           ! ... and then with external data
                           DO proc_shift = 1, para_env_sub%num_pe - 1
                              proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
                              proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

                          CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)
                          CALL get_group_dist(gd_B_virtual(ispin), proc_send, send_B_virtual_start, send_B_virtual_end, send_B_size)

                              ALLOCATE (external_ab(my_B_size(ispin), rec_B_size))
                              external_ab = 0.0_dp

                              CALL mp_sendrecv(local_ab(send_B_virtual_start:send_B_virtual_end, 1:my_B_size(ispin)), proc_send, &
                                               external_ab(1:my_B_size(ispin), 1:rec_B_size), proc_receive, &
                                               para_env_sub%group, tag)

                              DO b = 1, my_B_size(ispin)
                                 b_global = b + my_B_virtual_start(ispin) - 1
                                 DO a = 1, rec_B_size
                                    a_global = a + rec_B_virtual_start - 1
                                    my_Emp2_Ex = my_Emp2_Ex + sym_fac*local_ab(a_global, b)*external_ab(b, a)/ &
                                              (Eigenval(homo(ispin) + a_global, ispin) + Eigenval(homo(ispin) + b_global, ispin) - &
                                                  Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, ispin))
                                    IF (calc_forces .AND. (.NOT. my_alpha_beta_case)) &
                                       t_ab(a_global, b) = -(amp_fac*local_ab(a_global, b) - mp2_env%scale_T*external_ab(b, a))/ &
                                                           (Eigenval(homo(ispin) + a_global, ispin) + &
                                                            Eigenval(homo(ispin) + b_global, ispin) - &
                                                            Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, ispin))
                                 END DO
                              END DO

                              DEALLOCATE (external_ab)
                           END DO
                        END IF
                        CALL timestop(handle3)

                        IF (calc_forces) THEN
                           ! update P_ab, Gamma_P_ia
                           CALL mp2_update_P_gamma(mp2_env, para_env_sub, gd_B_virtual, &
                                                   Eigenval, homo, dimen_RI, iiB, jjB, my_B_size, &
                                                   my_B_virtual_end, my_B_virtual_start, my_i, my_j, virtual, &
                                                   sub_proc_map, local_ab, t_ab, my_local_i_aL, my_local_j_aL, &
                                                   my_open_shell_ss, Y_i_aP(:, :, iiB), Y_j_aP(:, :, jjB), local_ba, &
                                                   ispin, jspin, dgemm_counter)

                        END IF

                     END ASSOCIATE

                  END DO ! jjB
               END DO ! iiB

            ELSE
               CALL timeset(routineN//"_comm", handle3)
               ! No work to do and we know that we have to receive nothing, but send something
               ! send data to other proc
               DO proc_shift = 1, para_env_exchange%num_pe - 1
                  proc_send = proc_map(para_env_exchange%mepos + proc_shift)
                  proc_receive = proc_map(para_env_exchange%mepos - proc_shift)

                  send_ij_index = num_ij_pairs(proc_send)

                  IF (ij_index <= send_ij_index) THEN
                     ! something to send
                     ij_counter_send = (ij_index - MIN(1, integ_group_pos2color_sub(proc_send)))*ngroup + &
                                       integ_group_pos2color_sub(proc_send)
                     send_i = ij_map(1, ij_counter_send)
                     send_j = ij_map(2, ij_counter_send)
                     send_block_size = ij_map(3, ij_counter_send)

                     ! occupied i
                     CALL mp_send(BIb_C(ispin)%array(:, :, send_i:send_i + send_block_size - 1), &
                                  proc_send, tag, para_env_exchange%group)
                     ! occupied j
                     CALL mp_send(BIb_C(jspin)%array(:, :, send_j:send_j + send_block_size - 1), &
                                  proc_send, tag, para_env_exchange%group)
                  END IF
               END DO
               CALL timestop(handle3)
            END IF

            ! redistribute gamma
            IF (calc_forces) THEN
               CALL mp2_redistribute_gamma(mp2_env%ri_grad%Gamma_P_ia(ispin)%array, ij_index, my_B_size(ispin), &
                                           my_block_size, my_group_L_size, my_i, my_ij_pairs, ngroup, &
                                           num_integ_group, integ_group_pos2color_sub, num_ij_pairs, proc_map, &
                                           ij_map, ranges_info_array, Y_i_aP(:, :, 1:my_block_size), para_env_exchange, &
                                           gd_array%sizes, 1)
               CALL mp2_redistribute_gamma(mp2_env%ri_grad%Gamma_P_ia(jspin)%array, ij_index, my_B_size(jspin), &
                                           my_block_size, my_group_L_size, my_j, my_ij_pairs, ngroup, &
                                           num_integ_group, integ_group_pos2color_sub, num_ij_pairs, proc_map, &
                                           ij_map, ranges_info_array, Y_j_aP(:, :, 1:my_block_size), para_env_exchange, &
                                           gd_array%sizes, 2)
            END IF

         END DO
         CALL timestop(handle2)

         DEALLOCATE (local_i_aL)
         DEALLOCATE (local_j_aL)
         DEALLOCATE (ij_map)
         DEALLOCATE (num_ij_pairs)
         DEALLOCATE (local_ab)

         IF (calc_forces) THEN
            DEALLOCATE (Y_i_aP)
            DEALLOCATE (Y_j_aP)
            IF (ALLOCATED(t_ab)) THEN
               DEALLOCATE (t_ab)
            END IF
            DEALLOCATE (local_ba)

            ! here we check if there are almost degenerate ij
            ! pairs and we update P_ij with these contribution.
            ! If all pairs are degenerate with each other this step will scale O(N^6),
            ! if the number of degenerate pairs scales linearly with the system size
            ! this step will scale O(N^5).
            ! Start counting the number of almost degenerate ij pairs according
            ! to eps_canonical
            CALL quasi_degenerate_P_ij( &
               mp2_env, Eigenval(:, ispin:jspin), homo(ispin:jspin), virtual(ispin:jspin), my_open_shell_ss, &
               my_beta_beta_case, Bib_C(ispin:jspin), unit_nr, dimen_RI, &
               my_B_size(ispin:jspin), ngroup, num_integ_group, my_group_L_size, &
               color_sub, ranges_info_array, para_env_exchange, para_env_sub, proc_map, &
               my_B_virtual_start(ispin:jspin), my_B_virtual_end(ispin:jspin), gd_array%sizes, gd_B_virtual(ispin:jspin), &
               sub_proc_map, integ_group_pos2color_sub, dgemm_counter)

         END IF

         DEALLOCATE (integ_group_pos2color_sub)

         CALL mp_sum(my_Emp2_Cou, para_env%group)
         CALL mp_sum(my_Emp2_Ex, para_env%group)

         IF (calc_forces) THEN
            ! sum P_ab
            IF (.NOT. my_open_shell_ss) THEN
               IF (nspins == 1) mp2_env%ri_grad%P_ab(1)%array(:, :) = mp2_env%ri_grad%P_ab(1)%array(:, :)*2.0_dp
               sub_P_color = para_env_sub%mepos
               CALL cp_para_env_split(para_env_P, para_env, sub_P_color)
               DO kspin = ispin, jspin
                  CALL mp_sum(mp2_env%ri_grad%P_ab(kspin)%array, para_env_P%group)
               END DO
               ! release para_env_P
               CALL cp_para_env_release(para_env_P)
            END IF

            ! recover original information (before replication)
            DEALLOCATE (gd_array%sizes)
            iiB = SIZE(sizes_array_orig)
            ALLOCATE (gd_array%sizes(0:iiB - 1))
            gd_array%sizes(:) = sizes_array_orig
            DEALLOCATE (sizes_array_orig)

            ! make a copy of the original integrals (ia|Q)
            my_group_L_size = my_group_L_size_orig
            ALLOCATE (B_ia_Q(ispin:jspin))
            DO kspin = ispin, jspin
               ALLOCATE (B_ia_Q(kspin)%array(homo(kspin), my_B_size(kspin), my_group_L_size))
               B_ia_Q(kspin)%array = 0.0_dp
               DO jjB = 1, homo(kspin)
                  DO iiB = 1, my_B_size(kspin)
                     B_ia_Q(kspin)%array(jjB, iiB, 1:my_group_L_size) = &
                        BIb_C(kspin)%array(1:my_group_L_size, iiB, jjB)
                  END DO
               END DO
               DEALLOCATE (BIb_C(kspin)%array)
            END DO

            ! sum Gamma and dereplicate
            DO kspin = ispin, jspin
               ALLOCATE (BIb_C(kspin)%array(my_B_size(kspin), homo(kspin), my_group_L_size))
            END DO
            DO proc_shift = 1, para_env_rep%num_pe - 1
               ! invert order
               proc_send = proc_map_rep(para_env_rep%mepos - proc_shift)
               proc_receive = proc_map_rep(para_env_rep%mepos + proc_shift)

               start_point = ranges_info_array(3, proc_shift, para_env_exchange%mepos)
               end_point = ranges_info_array(4, proc_shift, para_env_exchange%mepos)

               DO kspin = ispin, jspin
                  CALL mp_sendrecv(mp2_env%ri_grad%Gamma_P_ia(kspin)%array(:, :, start_point:end_point), &
                                   proc_send, BIb_C(kspin)%array, proc_receive, para_env_rep%group, tag)
!$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
!$OMP          SHARED(mp2_env,BIb_C,kspin,homo,my_B_size,my_group_L_size)
                  mp2_env%ri_grad%Gamma_P_ia(kspin)%array(1:my_B_size(kspin), 1:homo(kspin), 1:my_group_L_size) = &
                     mp2_env%ri_grad%Gamma_P_ia(kspin)%array(1:my_B_size(kspin), 1:homo(kspin), 1:my_group_L_size) &
                     + BIb_C(kspin)%array(:, :, :)
!$OMP END PARALLEL WORKSHARE
               END DO
            END DO

            IF (.NOT. my_open_shell_ss) THEN
               ! Reduce size of Gamma_P_ia to save memory
               DO kspin = ispin, jspin
                  BIb_C(kspin)%array(:, :, :) = mp2_env%ri_grad%Gamma_P_ia(kspin)%array(:, :, 1:my_group_L_size)
                  DEALLOCATE (mp2_env%ri_grad%Gamma_P_ia(kspin)%array)
                  CALL MOVE_ALLOC(BIb_C(kspin)%array, mp2_env%ri_grad%Gamma_P_ia(kspin)%array)
               END DO

               ! B_ia_Q(kspin)%array will be deallocated inside of complete_gamma
               DO kspin = ispin, jspin
                  CALL complete_gamma(mp2_env, B_ia_Q(kspin)%array, dimen_RI, homo(kspin), &
                                      virtual(kspin), para_env, para_env_sub, ngroup, &
                                      my_group_L_size, my_group_L_start, my_group_L_end, &
                                      my_B_size(kspin), my_B_virtual_start(kspin), &
                                      gd_array, gd_B_virtual(kspin), &
                                      sub_proc_map, kspin)
               END DO
            ELSE
               DO kspin = ispin, jspin
                  DEALLOCATE (BIb_C(kspin)%array)
               END DO
               ALLOCATE (BIb_C(ispin)%array(my_group_L_size, my_B_size(ispin), homo(ispin)))
               BIb_C(ispin)%array = 0.0_dp
               ! copy the integrals (ia|Q) back
               DO jjB = 1, homo(ispin)
                  DO iiB = 1, my_B_size(ispin)
                     BIb_C(ispin)%array(1:my_group_L_size, iiB, jjB) = &
                        B_ia_Q(ispin)%array(jjB, iiB, 1:my_group_L_size)
                  END DO
               END DO
               ! There is only one B_ia_Q allocated in this case
               DEALLOCATE (B_ia_Q(ispin)%array)
            END IF
            DEALLOCATE (B_ia_Q)

         END IF

         IF (my_open_shell_SS .OR. my_alpha_beta_case) THEN
            IF (my_alpha_beta_case) THEN
               Emp2_S = Emp2_S + my_Emp2_Cou
               Emp2_Cou = Emp2_Cou + my_Emp2_Cou
            ELSE
               my_Emp2_Cou = my_Emp2_Cou*0.25_dp
               my_Emp2_EX = my_Emp2_EX*0.5_dp
               Emp2_T = Emp2_T + my_Emp2_Cou + my_Emp2_EX
               Emp2_Cou = Emp2_Cou + my_Emp2_Cou
               Emp2_EX = Emp2_EX + my_Emp2_EX
            END IF
         ELSE
            Emp2_Cou = Emp2_Cou + my_Emp2_Cou
            Emp2_EX = Emp2_EX + my_Emp2_EX
         END IF

         DEALLOCATE (proc_map)
         DEALLOCATE (sub_proc_map)
         DEALLOCATE (proc_map_rep)
         DEALLOCATE (ranges_info_array)

         IF (.NOT. my_open_shell_SS) THEN
            ! keep the array for the next calculations
            IF (ALLOCATED(BIb_C(ispin)%array)) DEALLOCATE (BIb_C(ispin)%array)
            CALL release_group_dist(gd_array)
            DO kspin = ispin, jspin
               CALL release_group_dist(gd_B_virtual(kspin))
            END DO
         END IF

         CALL cp_para_env_release(para_env_exchange)
         CALL cp_para_env_release(para_env_rep)

      END DO
      END DO

      ! We do not need this matrix later, so deallocate it here to safe memory
      IF (calc_forces) DEALLOCATE (mp2_env%ri_grad%PQ_half)
      IF (calc_forces .AND. .NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) &
         DEALLOCATE (mp2_env%ri_grad%operator_half)

      CALL dgemm_counter_write(dgemm_counter, para_env)

      ! release memory allocated by offload_dgemm when run on GPU. offload_gemm_ctx is null on cpu only runs
      CALL offload_gemm_destroy(mp2_env%offload_gemm_ctx)
      mp2_env%offload_gemm_ctx = C_NULL_PTR
      CALL timestop(handle)

   END SUBROUTINE mp2_ri_gpw_compute_en

! **************************************************************************************************
!> \brief ...
!> \param local_i_aL ...
!> \param ranges_info_array ...
!> \param BIb_C_rec ...
! **************************************************************************************************
   SUBROUTINE fill_local_i_aL(local_i_aL, ranges_info_array, BIb_C_rec)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: local_i_aL
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: ranges_info_array
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: BIb_C_rec

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fill_local_i_aL'

      INTEGER                                            :: end_point, handle, irep, Lend_pos, &
                                                            Lstart_pos, start_point

      CALL timeset(routineN, handle)

      DO irep = 1, SIZE(ranges_info_array, 2)
         Lstart_pos = ranges_info_array(1, irep)
         Lend_pos = ranges_info_array(2, irep)
         start_point = ranges_info_array(3, irep)
         end_point = ranges_info_array(4, irep)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
!$OMP          SHARED(BIb_C_rec,local_i_aL,Lstart_pos,Lend_pos,start_point,end_point)
         local_i_aL(Lstart_pos:Lend_pos, :, :) = BIb_C_rec(start_point:end_point, :, :)
!$OMP END PARALLEL WORKSHARE
      END DO

      CALL timestop(handle)

   END SUBROUTINE fill_local_i_aL

! **************************************************************************************************
!> \brief ...
!> \param BIb_C ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param para_env_exchange ...
!> \param para_env_rep ...
!> \param homo ...
!> \param proc_map_rep ...
!> \param sizes_array ...
!> \param my_B_size ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_new_group_L_size ...
!> \param new_sizes_array ...
!> \param ranges_info_array ...
! **************************************************************************************************
   SUBROUTINE replicate_iaK_2intgroup(BIb_C, para_env, para_env_sub, para_env_exchange, para_env_rep, &
                                      homo, proc_map_rep, &
                                      sizes_array, &
                                      my_B_size, &
                                      my_group_L_size, my_group_L_start, my_group_L_end, &
                                      my_new_group_L_size, new_sizes_array, ranges_info_array)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: BIb_C
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange, para_env_rep
      INTEGER, INTENT(IN)                                :: homo
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: proc_map_rep
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sizes_array
      INTEGER, INTENT(IN)                                :: my_B_size, my_group_L_size, &
                                                            my_group_L_start, my_group_L_end
      INTEGER, INTENT(INOUT)                             :: my_new_group_L_size
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: new_sizes_array
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: ranges_info_array

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_iaK_2intgroup'

      INTEGER                                            :: end_point, handle, i, max_L_size, &
                                                            proc_receive, proc_shift, start_point, &
                                                            sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: rep_ends_array, rep_sizes_array, &
                                                            rep_starts_array
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: my_info
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C_copy
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: BIb_C_gather

      CALL timeset(routineN, handle)

      ! create the replication group
      sub_sub_color = para_env_sub%mepos*para_env_exchange%num_pe + para_env_exchange%mepos
      CALL cp_para_env_split(para_env_rep, para_env, sub_sub_color)

      ! crate the proc maps
      ALLOCATE (proc_map_rep(-para_env_rep%num_pe:2*para_env_rep%num_pe - 1))
      DO i = 0, para_env_rep%num_pe - 1
         proc_map_rep(i) = i
         proc_map_rep(-i - 1) = para_env_rep%num_pe - i - 1
         proc_map_rep(para_env_rep%num_pe + i) = i
      END DO

      ! create the new limits for K according to the size
      ! of the integral group

      ! info array for replication
      ALLOCATE (rep_ends_array(0:para_env_rep%num_pe - 1))
      ALLOCATE (rep_starts_array(0:para_env_rep%num_pe - 1))
      ALLOCATE (rep_sizes_array(0:para_env_rep%num_pe - 1))

      CALL mp_allgather(my_group_L_size, rep_sizes_array, para_env_rep%group)
      CALL mp_allgather(my_group_L_start, rep_starts_array, para_env_rep%group)
      CALL mp_allgather(my_group_L_end, rep_ends_array, para_env_rep%group)

      ! calculate my_new_group_L_size according to sizes_array
      my_new_group_L_size = my_group_L_size

      ! Info of this process
      ALLOCATE (my_info(4, 0:para_env_rep%num_pe - 1))
      my_info(1, 0) = my_group_L_start
      my_info(2, 0) = my_group_L_end
      my_info(3, 0) = 1
      my_info(4, 0) = my_group_L_size

      DO proc_shift = 1, para_env_rep%num_pe - 1
         proc_receive = proc_map_rep(para_env_rep%mepos - proc_shift)

         my_new_group_L_size = my_new_group_L_size + rep_sizes_array(proc_receive)

         my_info(1, proc_shift) = rep_starts_array(proc_receive)
         my_info(2, proc_shift) = rep_ends_array(proc_receive)
         my_info(3, proc_shift) = my_info(4, proc_shift - 1) + 1
         my_info(4, proc_shift) = my_new_group_L_size

      END DO

      ALLOCATE (new_sizes_array(0:para_env_exchange%num_pe - 1))
      ALLOCATE (ranges_info_array(4, 0:para_env_rep%num_pe - 1, 0:para_env_exchange%num_pe - 1))
      CALL mp_allgather(my_new_group_L_size, new_sizes_array, para_env_exchange%group)
      CALL mp_allgather(my_info, ranges_info_array, para_env_exchange%group)

      DEALLOCATE (rep_sizes_array)
      DEALLOCATE (rep_starts_array)
      DEALLOCATE (rep_ends_array)

      ! replication scheme using mp_allgather
      ! get the max L size of the
      max_L_size = MAXVAL(sizes_array)

      ALLOCATE (BIb_C_copy(max_L_size, my_B_size, homo))
      BIb_C_copy = 0.0_dp
      BIb_C_copy(1:my_group_L_size, 1:my_B_size, 1:homo) = BIb_C

      DEALLOCATE (BIb_C)

      ALLOCATE (BIb_C_gather(max_L_size, my_B_size, homo, 0:para_env_rep%num_pe - 1))
      BIb_C_gather = 0.0_dp

      CALL mp_allgather(BIb_C_copy, BIb_C_gather, para_env_rep%group)

      DEALLOCATE (BIb_C_copy)

      ALLOCATE (BIb_C(my_new_group_L_size, my_B_size, homo))
      BIb_C = 0.0_dp

      ! reorder data
      DO proc_shift = 0, para_env_rep%num_pe - 1
         proc_receive = proc_map_rep(para_env_rep%mepos - proc_shift)

         start_point = my_info(3, proc_shift)
         end_point = my_info(4, proc_shift)

         BIb_C(start_point:end_point, 1:my_B_size, 1:homo) = &
            BIb_C_gather(1:end_point - start_point + 1, 1:my_B_size, 1:homo, proc_receive)

      END DO

      DEALLOCATE (BIb_C_gather)

      CALL timestop(handle)

   END SUBROUTINE replicate_iaK_2intgroup

! **************************************************************************************************
!> \brief ...
!> \param local_ab ...
!> \param t_ab ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param my_B_size ...
!> \param block_size ...
!> \param my_group_L_size ...
!> \param local_i_aL ...
!> \param local_j_aL ...
!> \param calc_forces ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param ispin ...
!> \param jspin ...
!> \param local_ba ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_allocate(local_ab, t_ab, mp2_env, homo, virtual, dimen_RI, my_B_size, &
                              block_size, my_group_L_size, &
                              local_i_aL, local_j_aL, calc_forces, &
                              Y_i_aP, Y_j_aP, ispin, jspin, &
                              local_ba)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: local_ab, t_ab
      TYPE(mp2_type)                                     :: mp2_env
      INTEGER, INTENT(IN)                                :: homo(2), virtual(2), dimen_RI, &
                                                            my_B_size(2), block_size, &
                                                            my_group_L_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: local_i_aL, local_j_aL
      LOGICAL, INTENT(IN)                                :: calc_forces
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: Y_i_aP, Y_j_aP
      INTEGER, INTENT(IN)                                :: ispin, jspin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: local_ba

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mp2_ri_allocate'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ALLOCATE (local_i_aL(dimen_RI, my_B_size(ispin), block_size))
      ALLOCATE (local_j_aL(dimen_RI, my_B_size(jspin), block_size))
      ALLOCATE (local_ab(virtual(ispin), my_B_size(jspin)))

      IF (calc_forces) THEN
         ALLOCATE (Y_i_aP(my_B_size(ispin), dimen_RI, block_size))
         Y_i_aP = 0.0_dp
         ! For  closed-shell, alpha-alpha and beta-beta my_B_size_beta=my_b_size
         ! Not for alpha-beta case: Y_j_aP_beta is sent and received as Y_j_aP
         ALLOCATE (Y_j_aP(my_B_size(jspin), dimen_RI, block_size))
         Y_j_aP = 0.0_dp
         ! Closed shell or alpha-alpha case
         IF (ispin == jspin) THEN
            ALLOCATE (mp2_env%ri_grad%P_ij(ispin)%array(homo(ispin), homo(ispin)))
            ALLOCATE (mp2_env%ri_grad%P_ab(ispin)%array(my_B_size(ispin), virtual(ispin)))
            mp2_env%ri_grad%P_ij(ispin)%array = 0.0_dp
            mp2_env%ri_grad%P_ab(ispin)%array = 0.0_dp
            ALLOCATE (mp2_env%ri_grad%Gamma_P_ia(ispin)%array(my_B_size(ispin), homo(ispin), my_group_L_size))
            mp2_env%ri_grad%Gamma_P_ia(ispin)%array = 0.0_dp

            ! For non-alpha-beta case we need amplitudes
            ALLOCATE (t_ab(virtual(ispin), my_B_size(jspin)))

            ! That is just a dummy. In that way, we can pass it as array to other routines w/o requirement for allocatable array
            ALLOCATE (local_ba(1, 1))
         ELSE
            ! We need more integrals
            ALLOCATE (local_ba(virtual(jspin), my_B_size(ispin)))
         END IF
      END IF
      !

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_allocate

! **************************************************************************************************
!> \brief ...
!> \param my_alpha_beta_case ...
!> \param total_ij_pairs ...
!> \param homo ...
!> \param homo_beta ...
!> \param num_IJ_blocks ...
!> \param block_size ...
!> \param ngroup ...
!> \param ij_map ...
!> \param color_sub ...
!> \param my_ij_pairs ...
!> \param my_open_shell_SS ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_communication(my_alpha_beta_case, total_ij_pairs, homo, homo_beta, num_IJ_blocks, &
                                   block_size, ngroup, ij_map, color_sub, my_ij_pairs, my_open_shell_SS, unit_nr)
      LOGICAL, INTENT(IN)                                :: my_alpha_beta_case
      INTEGER, INTENT(OUT)                               :: total_ij_pairs
      INTEGER, INTENT(IN)                                :: homo, homo_beta
      INTEGER, INTENT(OUT)                               :: num_IJ_blocks
      INTEGER, INTENT(IN)                                :: block_size, ngroup
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: ij_map
      INTEGER, INTENT(IN)                                :: color_sub
      INTEGER, INTENT(OUT)                               :: my_ij_pairs
      LOGICAL, INTENT(IN)                                :: my_open_shell_SS
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_communication'

      INTEGER :: assigned_blocks, first_I_block, first_J_block, handle, iiB, ij_block_counter, &
         ij_counter, jjB, last_i_block, last_J_block, num_block_per_group, total_ij_block, &
         total_ij_pairs_blocks
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_marker

! Calculate the maximum number of ij pairs that have to be computed
! among groups

      CALL timeset(routineN, handle)

      IF (.NOT. my_alpha_beta_case) THEN
         total_ij_pairs = homo*(1 + homo)/2
         num_IJ_blocks = homo/block_size - 1

         first_I_block = 1
         last_i_block = block_size*(num_IJ_blocks - 1)

         first_J_block = block_size + 1
         last_J_block = block_size*(num_IJ_blocks + 1)

         ij_block_counter = 0
         DO iiB = first_I_block, last_i_block, block_size
            DO jjB = iiB + block_size, last_J_block, block_size
               ij_block_counter = ij_block_counter + 1
            END DO
         END DO

         total_ij_block = ij_block_counter
         num_block_per_group = total_ij_block/ngroup
         assigned_blocks = num_block_per_group*ngroup

         total_ij_pairs_blocks = assigned_blocks + (total_ij_pairs - assigned_blocks*(block_size**2))

         ALLOCATE (ij_marker(homo, homo))
         ij_marker = 0
         ALLOCATE (ij_map(3, total_ij_pairs_blocks))
         ij_map = 0
         ij_counter = 0
         my_ij_pairs = 0
         DO iiB = first_I_block, last_i_block, block_size
            DO jjB = iiB + block_size, last_J_block, block_size
               IF (ij_counter + 1 > assigned_blocks) EXIT
               ij_counter = ij_counter + 1
               ij_marker(iiB:iiB + block_size - 1, jjB:jjB + block_size - 1) = 1
               ij_map(1, ij_counter) = iiB
               ij_map(2, ij_counter) = jjB
               ij_map(3, ij_counter) = block_size
               IF (MOD(ij_counter, ngroup) == color_sub) my_ij_pairs = my_ij_pairs + 1
            END DO
         END DO
         DO iiB = 1, homo
            DO jjB = iiB, homo
               IF (ij_marker(iiB, jjB) == 0) THEN
                  ij_counter = ij_counter + 1
                  ij_map(1, ij_counter) = iiB
                  ij_map(2, ij_counter) = jjB
                  ij_map(3, ij_counter) = 1
                  IF (MOD(ij_counter, ngroup) == color_sub) my_ij_pairs = my_ij_pairs + 1
               END IF
            END DO
         END DO
         DEALLOCATE (ij_marker)

         IF ((.NOT. my_open_shell_SS)) THEN
            IF (unit_nr > 0) THEN
               IF (block_size == 1) THEN
                  WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.1)") &
                     "RI_INFO| Percentage of ij pairs communicated with block size 1:", 100.0_dp
               ELSE
                  WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.1)") &
                     "RI_INFO| Percentage of ij pairs communicated with block size 1:", &
                     100.0_dp*REAL((total_ij_pairs - assigned_blocks*(block_size**2)), KIND=dp)/REAL(total_ij_pairs, KIND=dp)
               END IF
               CALL m_flush(unit_nr)
            END IF
         END IF

      ELSE
         ! alpha-beta case no index symmetry
         total_ij_pairs = homo*homo_beta
         ALLOCATE (ij_map(3, total_ij_pairs))
         ij_map = 0
         ij_counter = 0
         my_ij_pairs = 0
         DO iiB = 1, homo
            DO jjB = 1, homo_beta
               ij_counter = ij_counter + 1
               ij_map(1, ij_counter) = iiB
               ij_map(2, ij_counter) = jjB
               ij_map(3, ij_counter) = 1
               IF (MOD(ij_counter, ngroup) == color_sub) my_ij_pairs = my_ij_pairs + 1
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_communication

! **************************************************************************************************
!> \brief ...
!> \param BIb_C ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param homo ...
!> \param color_sub ...
!> \param sizes_array ...
!> \param calc_forces ...
!> \param integ_group_size ...
!> \param my_B_size ...
!> \param iiB ...
!> \param my_group_L_end ...
!> \param my_group_L_size ...
!> \param my_group_L_size_orig ...
!> \param my_group_L_start ...
!> \param my_new_group_L_size ...
!> \param sub_sub_color ...
!> \param integ_group_pos2color_sub ...
!> \param proc_map ...
!> \param proc_map_rep ...
!> \param sizes_array_orig ...
!> \param sub_proc_map ...
!> \param ranges_info_array ...
!> \param para_env_exchange ...
!> \param para_env_rep ...
!> \param num_integ_group ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_create_group(BIb_C, para_env, para_env_sub, homo, color_sub, &
                                  sizes_array, calc_forces, &
                                  integ_group_size, my_B_size, iiB, my_group_L_end, &
                                  my_group_L_size, my_group_L_size_orig, my_group_L_start, my_new_group_L_size, &
                                  sub_sub_color, integ_group_pos2color_sub, &
                                  proc_map, proc_map_rep, sizes_array_orig, &
                                  sub_proc_map, ranges_info_array, para_env_exchange, para_env_rep, num_integ_group)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: BIb_C
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: homo, color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: sizes_array
      LOGICAL, INTENT(IN)                                :: calc_forces
      INTEGER, INTENT(IN)                                :: integ_group_size, my_B_size
      INTEGER, INTENT(OUT)                               :: iiB
      INTEGER, INTENT(IN)                                :: my_group_L_end
      INTEGER, INTENT(INOUT)                             :: my_group_L_size
      INTEGER, INTENT(OUT)                               :: my_group_L_size_orig
      INTEGER, INTENT(IN)                                :: my_group_L_start
      INTEGER, INTENT(INOUT)                             :: my_new_group_L_size
      INTEGER, INTENT(OUT)                               :: sub_sub_color
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: integ_group_pos2color_sub, proc_map, &
                                                            proc_map_rep, sizes_array_orig, &
                                                            sub_proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: ranges_info_array
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange, para_env_rep
      INTEGER, INTENT(IN)                                :: num_integ_group

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_create_group'

      INTEGER                                            :: handle, i
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: new_sizes_array

      CALL timeset(routineN, handle)
      !
      sub_sub_color = para_env_sub%mepos*num_integ_group + color_sub/integ_group_size
      CALL cp_para_env_split(para_env_exchange, para_env, sub_sub_color)

      ! create the proc maps
      ALLOCATE (proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe - 1))
      DO i = 0, para_env_exchange%num_pe - 1
         proc_map(i) = i
         proc_map(-i - 1) = para_env_exchange%num_pe - i - 1
         proc_map(para_env_exchange%num_pe + i) = i
      END DO

      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe - 1))
      DO i = 0, para_env_sub%num_pe - 1
         sub_proc_map(i) = i
         sub_proc_map(-i - 1) = para_env_sub%num_pe - i - 1
         sub_proc_map(para_env_sub%num_pe + i) = i
      END DO

      CALL replicate_iaK_2intgroup(BIb_C, para_env, para_env_sub, para_env_exchange, para_env_rep, &
                                   homo, proc_map_rep, &
                                   sizes_array, &
                                   my_B_size, &
                                   my_group_L_size, my_group_L_start, my_group_L_end, &
                                   my_new_group_L_size, new_sizes_array, ranges_info_array)

      ALLOCATE (integ_group_pos2color_sub(0:para_env_exchange%num_pe - 1))
      CALL mp_allgather(color_sub, integ_group_pos2color_sub, para_env_exchange%group)

      IF (calc_forces) THEN
         iiB = SIZE(sizes_array)
         ALLOCATE (sizes_array_orig(0:iiB - 1))
         sizes_array_orig(:) = sizes_array
      END IF

      my_group_L_size_orig = my_group_L_size
      my_group_L_size = my_new_group_L_size
      DEALLOCATE (sizes_array)

      ALLOCATE (sizes_array(0:integ_group_size - 1))
      sizes_array(:) = new_sizes_array

      DEALLOCATE (new_sizes_array)
      !
      CALL timestop(handle)

   END SUBROUTINE mp2_ri_create_group

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param gd_array ...
!> \param gd_B_virtual ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param integ_group_size ...
!> \param ngroup ...
!> \param num_integ_group ...
!> \param virtual ...
!> \param my_alpha_beta_case ...
!> \param my_open_shell_SS ...
!> \param calc_forces ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_get_integ_group_size(mp2_env, para_env, para_env_sub, gd_array, gd_B_virtual, &
                                          homo, dimen_RI, unit_nr, &
                                          integ_group_size, &
                                          ngroup, num_integ_group, &
                                          virtual, my_alpha_beta_case, &
                                          my_open_shell_SS, calc_forces)
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array, gd_B_virtual
      INTEGER, INTENT(IN)                                :: homo, dimen_RI, unit_nr
      INTEGER, INTENT(OUT)                               :: integ_group_size, ngroup, num_integ_group
      INTEGER, INTENT(IN)                                :: virtual
      LOGICAL, INTENT(IN)                                :: my_alpha_beta_case, my_open_shell_SS, &
                                                            calc_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_get_integ_group_size'

      INTEGER                                            :: best_block_size, best_integ_group_size, &
                                                            handle, iiB, min_integ_group_size
      INTEGER(KIND=int_8)                                :: mem
      REAL(KIND=dp)                                      :: mem_for_aK, mem_for_iaK, mem_for_rep, &
                                                            mem_min, mem_per_group, mem_real

      CALL timeset(routineN, handle)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! Calculate available memory and create integral group according to that
      ! mem_for_iaK is the memory needed for storing the 3 centre integrals
      mem_for_iaK = REAL(homo, KIND=dp)*virtual*dimen_RI*8.0_dp/(1024_dp**2)
      mem_for_aK = REAL(virtual, KIND=dp)*dimen_RI*8.0_dp/(1024_dp**2)

      CALL m_memory(mem)
      mem_real = (mem + 1024*1024 - 1)/(1024*1024)
      ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
      ! has not been given back to the OS yet.
      CALL mp_min(mem_real, para_env%group)

      ! BIB_C_copy/external_i_aL/Gamma_P_ia
      mem_min = MAX(REAL(homo, KIND=dp)*maxsize(gd_array), REAL(dimen_RI, KIND=dp))*maxsize(gd_B_virtual)*8.0_dp/(1024**2)
      IF (calc_forces) THEN
         mem_min = MAX(mem_min, 2.0_dp*maxsize(gd_array)*maxsize(gd_B_virtual)*8.0_dp/(1024**2))
      END IF
      ! BIB_C
      mem_min = REAL(homo, KIND=dp)*maxsize(gd_B_virtual)*maxsize(gd_array)*8.0_dp/(1024**2)
      ! local_i_aL+local_j_aL
      mem_min = mem_min + 2.0_dp*maxsize(gd_B_virtual)*REAL(dimen_RI, KIND=dp)*8.0_dp/(1024**2)
      ! local_ab
      mem_min = mem_min + REAL(virtual, KIND=dp)*maxsize(gd_B_virtual)*8.0_dp/(1024**2)

      IF (calc_forces) THEN
         ! Y_i_aP+Y_j_aP
         mem_min = mem_min + 2.0_dp*maxsize(gd_B_virtual)*dimen_RI*8.0_dp/(1024**2)
         ! local_ba/t_ab
         mem_min = mem_min + REAL(virtual, KIND=dp)*maxsize(gd_B_virtual)*8.0_dp/(1024**2)
         IF (.NOT. my_alpha_beta_case) THEN
            ! P_ij
            mem_min = mem_min + REAL(homo, KIND=dp)*homo*8.0_dp/(1024**2)
            ! P_ab
            mem_min = mem_min + REAL(virtual, KIND=dp)*maxsize(gd_B_virtual)*8.0_dp/(1024**2)
         END IF
      END IF

      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'RI_INFO| Minimum required memory per MPI process:', &
            mem_min, ' MiB'
      END IF

      mem_real = mp2_env%mp2_memory

      mem_per_group = mem_real*para_env_sub%num_pe

      ! here we try to find the best block_size and integ_group_size
      best_integ_group_size = ngroup
      best_block_size = 1

      ! in the open shell case no replication and no block communication is done
      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         ! Here we split the memory half for the communication, half for replication
         IF (mp2_env%ri_mp2%block_size > 0) THEN
            best_block_size = mp2_env%ri_mp2%block_size
            mem_for_rep = MAX(mem_min, mem_per_group - 2.0_dp*mem_for_aK*best_block_size)
         ELSE
            mem_for_rep = mem_per_group/2.0_dp
         END IF
         ! calculate the minimum replication group size according to the available memory
         min_integ_group_size = CEILING(2.0_dp*mem_for_iaK/mem_for_rep)

         integ_group_size = MIN(min_integ_group_size, ngroup) - 1
         DO iiB = min_integ_group_size + 1, ngroup
            integ_group_size = integ_group_size + 1
            ! check that the ngroup is a multiple of  integ_group_size
            IF (MOD(ngroup, integ_group_size) /= 0) CYCLE
            ! check that the integ group size is not too small (10% is empirical for now)
            IF (REAL(integ_group_size, KIND=dp)/REAL(ngroup, KIND=dp) < 0.1_dp) CYCLE

            best_integ_group_size = integ_group_size
            EXIT
         END DO
      END IF

      integ_group_size = best_integ_group_size

      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Group size for integral replication:", integ_group_size*para_env_sub%num_pe
            CALL m_flush(unit_nr)
         END IF
      END IF

      num_integ_group = ngroup/integ_group_size

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_get_integ_group_size

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param gd_array ...
!> \param gd_B_virtual ...
!> \param homo ...
!> \param dimen_RI ...
!> \param unit_nr ...
!> \param block_size ...
!> \param ngroup ...
!> \param num_IJ_blocks ...
!> \param virtual ...
!> \param my_alpha_beta_case ...
!> \param my_open_shell_SS ...
!> \param calc_forces ...
! **************************************************************************************************
   SUBROUTINE mp2_ri_get_block_size(mp2_env, para_env, para_env_sub, gd_array, gd_B_virtual, &
                                    homo, dimen_RI, unit_nr, &
                                    block_size, ngroup, num_IJ_blocks, &
                                    virtual, my_alpha_beta_case, &
                                    my_open_shell_SS, calc_forces)
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array, gd_B_virtual
      INTEGER, INTENT(IN)                                :: homo, dimen_RI, unit_nr
      INTEGER, INTENT(OUT)                               :: block_size, ngroup, num_IJ_blocks
      INTEGER, INTENT(IN)                                :: virtual
      LOGICAL, INTENT(IN)                                :: my_alpha_beta_case, my_open_shell_SS, &
                                                            calc_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_ri_get_block_size'

      INTEGER                                            :: best_block_size, handle
      INTEGER(KIND=int_8)                                :: mem
      REAL(KIND=dp)                                      :: mem_for_aK, mem_for_iaK, mem_min, &
                                                            mem_per_group, mem_real

      CALL timeset(routineN, handle)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! Calculate available memory and create integral group according to that
      ! mem_for_iaK is the memory needed for storing the 3 centre integrals
      mem_for_iaK = REAL(homo, KIND=dp)*virtual*dimen_RI*8.0_dp/(1024_dp**2)
      mem_for_aK = REAL(virtual, KIND=dp)*dimen_RI*8.0_dp/(1024_dp**2)

      CALL m_memory(mem)
      mem_real = (mem + 1024*1024 - 1)/(1024*1024)
      ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
      ! has not been given back to the OS yet.
      CALL mp_min(mem_real, para_env%group)

      ! BIB_C_copy/external_i_aL/Gamma_P_ia
      mem_min = MAX(REAL(homo, KIND=dp)*maxsize(gd_array), REAL(dimen_RI, KIND=dp))*maxsize(gd_B_virtual)*8.0_dp/(1024**2)
      IF (calc_forces) THEN
         mem_min = MAX(mem_min, 2.0_dp*maxsize(gd_array)*maxsize(gd_B_virtual)*8.0_dp/(1024**2))
      END IF
      ! BIB_C
      mem_min = REAL(homo, KIND=dp)*maxsize(gd_B_virtual)*maxsize(gd_array)*8.0_dp/(1024**2)
      ! local_i_aL+local_j_aL
      mem_min = mem_min + 2.0_dp*maxsize(gd_B_virtual)*REAL(dimen_RI, KIND=dp)*8.0_dp/(1024**2)
      ! local_ab
      mem_min = mem_min + REAL(virtual, KIND=dp)*maxsize(gd_B_virtual)*8.0_dp/(1024**2)

      IF (calc_forces) THEN
         ! Y_i_aP+Y_j_aP
         mem_min = mem_min + 2.0_dp*maxsize(gd_B_virtual)*dimen_RI*8.0_dp/(1024**2)
         ! local_ba/t_ab
         mem_min = mem_min + REAL(virtual, KIND=dp)*maxsize(gd_B_virtual)*8.0_dp/(1024**2)
         IF (.NOT. my_alpha_beta_case) THEN
            ! P_ij
            mem_min = mem_min + REAL(homo, KIND=dp)*homo*8.0_dp/(1024**2)
            ! P_ab
            mem_min = mem_min + REAL(virtual, KIND=dp)*maxsize(gd_B_virtual)*8.0_dp/(1024**2)
         END IF
      END IF

      mem_real = mp2_env%mp2_memory

      mem_per_group = mem_real*para_env_sub%num_pe

      best_block_size = 1

      ! in the open shell case no replication and no block communication is done
      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         ! Here we split the memory half for the communication, half for replication
         IF (mp2_env%ri_mp2%block_size > 0) THEN
            best_block_size = mp2_env%ri_mp2%block_size
         END IF

         IF (.NOT. (mp2_env%ri_mp2%block_size > 0)) THEN
            DO
               num_IJ_blocks = (homo/best_block_size)
               num_IJ_blocks = (num_IJ_blocks*num_IJ_blocks - num_IJ_blocks)/2
               IF (num_IJ_blocks > ngroup .OR. best_block_size == 1) THEN
                  EXIT
               ELSE
                  best_block_size = best_block_size - 1
               END IF
            END DO
         END IF

         ! check that best_block_size is not bigger than homo/2-1
         best_block_size = MIN(MAX(homo/2 - 1 + MOD(homo, 2), 1), best_block_size)
      END IF

      block_size = best_block_size

      IF ((.NOT. my_open_shell_SS) .AND. (.NOT. my_alpha_beta_case)) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "RI_INFO| Block size:", block_size
            CALL m_flush(unit_nr)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE mp2_ri_get_block_size

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env_sub ...
!> \param gd_B_virtual ...
!> \param Eigenval ...
!> \param homo ...
!> \param dimen_RI ...
!> \param iiB ...
!> \param jjB ...
!> \param my_B_size ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
!> \param my_i ...
!> \param my_j ...
!> \param virtual ...
!> \param sub_proc_map ...
!> \param local_ab ...
!> \param t_ab ...
!> \param my_local_i_aL ...
!> \param my_local_j_aL ...
!> \param open_ss ...
!> \param Y_i_aP ...
!> \param Y_j_aP ...
!> \param local_ba ...
!> \param ispin ...
!> \param jspin ...
!> \param dgemm_counter ...
! **************************************************************************************************
   SUBROUTINE mp2_update_P_gamma(mp2_env, para_env_sub, gd_B_virtual, &
                                 Eigenval, homo, dimen_RI, iiB, jjB, my_B_size, &
                                 my_B_virtual_end, my_B_virtual_start, my_i, my_j, virtual, sub_proc_map, local_ab, &
                                 t_ab, my_local_i_aL, my_local_j_aL, open_ss, Y_i_aP, Y_j_aP, &
                                 local_ba, ispin, jspin, dgemm_counter)
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_sub
      TYPE(group_dist_d1_type), DIMENSION(:), INTENT(IN) :: gd_B_virtual
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      INTEGER, INTENT(IN)                                :: dimen_RI, iiB, jjB
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_B_size, my_B_virtual_end, &
                                                            my_B_virtual_start
      INTEGER, INTENT(IN)                                :: my_i, my_j
      INTEGER, DIMENSION(:), INTENT(IN)                  :: virtual
      INTEGER, DIMENSION(-para_env_sub%num_pe:), &
         INTENT(IN)                                      :: sub_proc_map
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: local_ab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), TARGET :: t_ab, my_local_i_aL, my_local_j_aL
      LOGICAL, INTENT(IN)                                :: open_ss
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         TARGET                                          :: Y_i_aP, Y_j_aP, local_ba
      INTEGER, INTENT(IN)                                :: ispin, jspin
      TYPE(dgemm_counter_type), INTENT(INOUT)            :: dgemm_counter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_update_P_gamma'

      INTEGER :: a, b, b_global, handle, proc_receive, proc_send, proc_shift, rec_B_size, &
         rec_B_virtual_end, rec_B_virtual_start, send_B_size, send_B_virtual_end, &
         send_B_virtual_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         TARGET                                          :: external_ab, send_ab
      REAL(KIND=dp)                                      :: factor, P_ij_diag
      LOGICAL                                            :: alpha_beta

!
!  In alpha-beta case Y_j_aP_beta is sent and is received as Y_j_aP
!

      CALL timeset(routineN//"_Pia", handle)

! Find out whether we have an alpha-beta case
! update P_ab, Gamma_P_ia
! First, P_ab
      alpha_beta = .NOT. (ispin == jspin)
      IF (open_ss) THEN
         factor = 1.0_dp
      ELSE
         factor = 2.0_dp
      END IF
      ! divide the (ia|jb) integrals by Delta_ij^ab
      DO b = 1, my_B_size(jspin)
         b_global = b + my_B_virtual_start(jspin) - 1
         DO a = 1, virtual(ispin)
            local_ab(a, b) = -local_ab(a, b)/ &
                             (Eigenval(homo(ispin) + a, ispin) + Eigenval(homo(jspin) + b_global, jspin) - &
                              Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, jspin))
         END DO
      END DO
      IF (.NOT. (alpha_beta)) THEN
         P_ij_diag = -SUM(local_ab*t_ab)*factor
      ELSE
         ! update diagonal part of P_ij
         P_ij_diag = -SUM(local_ab*local_ab)*mp2_env%scale_S
         ! More integrals needed only for alpha-beta case: local_ba
         DO b = 1, my_B_size(ispin)
            b_global = b + my_B_virtual_start(ispin) - 1
            DO a = 1, virtual(jspin)
               local_ba(a, b) = -local_ba(a, b)/ &
                                (Eigenval(homo(jspin) + a, jspin) + Eigenval(homo(ispin) + b_global, ispin) - &
                                 Eigenval(my_i + iiB - 1, ispin) - Eigenval(my_j + jjB - 1, jspin))
            END DO
         END DO
      END IF

      ! P_ab and add diagonal part of P_ij

      CALL dgemm_counter_start(dgemm_counter)
      IF (.NOT. (alpha_beta)) THEN
         CALL offload_dgemm('T', 'N', my_B_size(ispin), my_B_size(ispin), virtual(ispin), 1.0_dp, &
                            t_ab, virtual(ispin), local_ab, virtual(ispin), &
                            1.0_dp, mp2_env%ri_grad%P_ab(ispin)%array(:, &
                                     my_B_virtual_start(ispin):my_B_virtual_end(ispin)), my_B_size(ispin), mp2_env%offload_gemm_ctx)
         mp2_env%ri_grad%P_ij(ispin)%array(my_i + iiB - 1, my_i + iiB - 1) = &
            mp2_env%ri_grad%P_ij(ispin)%array(my_i + iiB - 1, my_i + iiB - 1) + P_ij_diag
      ELSE
         CALL offload_dgemm('T', 'N', my_B_size(ispin), my_B_size(ispin), virtual(jspin), mp2_env%scale_S, &
                            local_ba, virtual(jspin), local_ba, virtual(jspin), 1.0_dp, &
                        mp2_env%ri_grad%P_ab(ispin)%array(:, my_B_virtual_start(ispin):my_B_virtual_end(ispin)), my_B_size(ispin), &
                            mp2_env%offload_gemm_ctx)

         mp2_env%ri_grad%P_ij(ispin)%array(my_i + iiB - 1, my_i + iiB - 1) = &
            mp2_env%ri_grad%P_ij(ispin)%array(my_i + iiB - 1, my_i + iiB - 1) + P_ij_diag

         CALL offload_dgemm('T', 'N', my_B_size(jspin), my_B_size(jspin), virtual(ispin), mp2_env%scale_S, &
                            local_ab, virtual(ispin), local_ab, virtual(ispin), 1.0_dp, &
                        mp2_env%ri_grad%P_ab(jspin)%array(:, my_B_virtual_start(jspin):my_B_virtual_end(jspin)), my_B_size(jspin), &
                            mp2_env%offload_gemm_ctx)

         mp2_env%ri_grad%P_ij(jspin)%array(my_j + jjB - 1, my_j + jjB - 1) = &
            mp2_env%ri_grad%P_ij(jspin)%array(my_j + jjB - 1, my_j + jjB - 1) + P_ij_diag
      END IF
      ! The summation is over unique pairs. In alpha-beta case, all pairs are unique: subroutine is called for
      ! both i^alpha,j^beta and i^beta,j^alpha. Formally, my_i can be equal to my_j, but they are different
      ! due to spin in alpha-beta case.
      IF ((my_i /= my_j) .AND. (.NOT. alpha_beta)) THEN

         CALL offload_dgemm('N', 'T', my_B_size(ispin), virtual(ispin), my_B_size(ispin), 1.0_dp, &
                            t_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :), my_B_size(ispin), &
                            local_ab, virtual(ispin), &
                            1.0_dp, mp2_env%ri_grad%P_ab(ispin)%array, my_B_size(ispin), &
                            mp2_env%offload_gemm_ctx)

         mp2_env%ri_grad%P_ij(ispin)%array(my_j + jjB - 1, my_j + jjB - 1) = &
            mp2_env%ri_grad%P_ij(ispin)%array(my_j + jjB - 1, my_j + jjB - 1) + P_ij_diag
      END IF
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

         CALL get_group_dist(gd_B_virtual(jspin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)
         CALL get_group_dist(gd_B_virtual(jspin), proc_send, send_B_virtual_start, send_B_virtual_end, send_B_size)

         ALLOCATE (external_ab(virtual(ispin), rec_B_size))
         external_ab = 0.0_dp

         CALL mp_sendrecv(local_ab, proc_send, &
                          external_ab, proc_receive, &
                          para_env_sub%group)

         IF (.NOT. (alpha_beta)) THEN
            CALL offload_dgemm('T', 'N', my_B_size(ispin), rec_B_size, virtual(ispin), 1.0_dp, &
                               t_ab, virtual(ispin), external_ab, virtual(ispin), &
                            1.0_dp, mp2_env%ri_grad%P_ab(ispin)%array(:, rec_B_virtual_start:rec_B_virtual_end), my_B_size(ispin), &
                               mp2_env%offload_gemm_ctx)
         ELSE
            CALL offload_dgemm('T', 'N', my_B_size(jspin), rec_B_size, virtual(ispin), mp2_env%scale_S, &
                               local_ab, virtual(ispin), external_ab, virtual(ispin), &
                               1.0_dp, mp2_env%ri_grad%P_ab(jspin)%array(:, rec_B_virtual_start:rec_B_virtual_end), &
                               my_B_size(jspin), mp2_env%offload_gemm_ctx)

            ! For alpha-beta part of alpha-density we need a new parallel code
            DEALLOCATE (external_ab)
            ! And new external_ab (of a different size)
            CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)
            CALL get_group_dist(gd_B_virtual(ispin), proc_send, send_B_virtual_start, send_B_virtual_end, send_B_size)
            ALLOCATE (external_ab(virtual(jspin), rec_B_size))
            external_ab = 0.0_dp
            CALL mp_sendrecv(local_ba, proc_send, &
                             external_ab, proc_receive, &
                             para_env_sub%group)
            CALL offload_dgemm('T', 'N', my_B_size(ispin), rec_B_size, virtual(jspin), mp2_env%scale_S, &
                               local_ba, virtual(jspin), external_ab, virtual(jspin), &
                            1.0_dp, mp2_env%ri_grad%P_ab(ispin)%array(:, rec_B_virtual_start:rec_B_virtual_end), my_B_size(ispin), &
                               mp2_env%offload_gemm_ctx)
         END IF

         DEALLOCATE (external_ab)

         IF ((my_i /= my_j) .AND. (.NOT. alpha_beta)) THEN
            ALLOCATE (external_ab(my_B_size(ispin), virtual(ispin)))
            external_ab = 0.0_dp

            ALLOCATE (send_ab(send_B_size, virtual(ispin)))
            send_ab = 0.0_dp

            CALL offload_dgemm('N', 'T', send_B_size, virtual(ispin), my_B_size(ispin), 1.0_dp, &
                               t_ab(send_B_virtual_start:send_B_virtual_end, :), send_B_size, &
                               local_ab, virtual(ispin), &
                               0.0_dp, send_ab, send_B_size, &
                               mp2_env%offload_gemm_ctx)
            CALL mp_sendrecv(send_ab, proc_send, &
                             external_ab, proc_receive, &
                             para_env_sub%group)

            mp2_env%ri_grad%P_ab(ispin)%array(:, :) = mp2_env%ri_grad%P_ab(ispin)%array + external_ab

            DEALLOCATE (external_ab)
            DEALLOCATE (send_ab)
         END IF

      END DO
      IF (.NOT. alpha_beta) THEN
         IF (my_i /= my_j) THEN
            CALL dgemm_counter_stop(dgemm_counter, 2*my_B_size(ispin), virtual(ispin), virtual(ispin))
         ELSE
            CALL dgemm_counter_stop(dgemm_counter, my_B_size(ispin), virtual(ispin), virtual(ispin))
         END IF
      ELSE
         CALL dgemm_counter_stop(dgemm_counter, SUM(my_B_size), virtual(ispin), virtual(jspin))
      END IF
      CALL timestop(handle)

      ! Now, Gamma_P_ia (made of Y_ia_P)

      CALL timeset(routineN//"_Gamma", handle)
      CALL dgemm_counter_start(dgemm_counter)
      IF (.NOT. alpha_beta) THEN
         ! Alpha-alpha, beta-beta and closed shell
         CALL offload_dgemm('N', 'T', my_B_size(ispin), dimen_RI, my_B_size(ispin), 1.0_dp, &
                            t_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :), my_B_size(ispin), &
                            my_local_j_aL, dimen_RI, 1.0_dp, Y_i_aP, my_B_size(ispin), &
                            mp2_env%offload_gemm_ctx)
      ELSE ! Alpha-beta
         CALL offload_dgemm('N', 'T', my_B_size(ispin), dimen_RI, my_B_size(jspin), mp2_env%scale_S, &
                            local_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :), my_B_size(ispin), &
                            my_local_j_aL, dimen_RI, 1.0_dp, Y_i_aP, my_B_size(ispin), mp2_env%offload_gemm_ctx)
         CALL offload_dgemm('T', 'T', my_B_size(jspin), dimen_RI, my_B_size(ispin), mp2_env%scale_S, &
                            local_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :), my_B_size(ispin), &
                            my_local_i_aL, dimen_RI, 1.0_dp, Y_j_aP, my_B_size(jspin), mp2_env%offload_gemm_ctx)
      END IF

      ALLOCATE (external_ab(my_B_size(ispin), dimen_RI))
      external_ab = 0.0_dp
      !
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
         proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

         CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)
         CALL get_group_dist(gd_B_virtual(ispin), proc_send, send_B_virtual_start, send_B_virtual_end, send_B_size)

         IF (.NOT. alpha_beta) THEN
            ALLOCATE (send_ab(send_B_size, dimen_RI))
            send_ab = 0.0_dp
            CALL offload_dgemm('N', 'T', send_B_size, dimen_RI, my_B_size(ispin), 1.0_dp, &
                               t_ab(send_B_virtual_start:send_B_virtual_end, :), send_B_size, &
                               my_local_j_aL, dimen_RI, 0.0_dp, send_ab, send_B_size, &
                               mp2_env%offload_gemm_ctx)
            CALL mp_sendrecv(send_ab, proc_send, external_ab, proc_receive, para_env_sub%group)

            Y_i_aP(:, :) = Y_i_aP + external_ab

            DEALLOCATE (send_ab)
         ELSE ! Alpha-beta case
            ! Alpha-alpha part
            ALLOCATE (send_ab(send_B_size, dimen_RI))
            send_ab = 0.0_dp
            CALL offload_dgemm('N', 'T', send_B_size, dimen_RI, my_B_size(jspin), mp2_env%scale_S, &
                               local_ab(send_B_virtual_start:send_B_virtual_end, :), send_B_size, &
                               my_local_j_aL, dimen_RI, 0.0_dp, send_ab, send_B_size, &
                               mp2_env%offload_gemm_ctx)
            CALL mp_sendrecv(send_ab, proc_send, external_ab, proc_receive, para_env_sub%group)
            Y_i_aP(:, :) = Y_i_aP + external_ab
            DEALLOCATE (send_ab)
         END IF
      END DO
      DEALLOCATE (external_ab)

      IF (alpha_beta) THEN
         ! For beta-beta part (in alpha-beta case) we need a new parallel code
         ALLOCATE (external_ab(my_B_size(jspin), dimen_RI))
         external_ab = 0.0_dp
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
            proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

            CALL get_group_dist(gd_B_virtual(jspin), proc_send, send_B_virtual_start, send_B_virtual_end, send_B_size)
            ALLOCATE (send_ab(send_B_size, dimen_RI))
            send_ab = 0.0_dp
            CALL offload_dgemm('N', 'T', send_B_size, dimen_RI, my_B_size(ispin), mp2_env%scale_S, &
                               local_ba(send_B_virtual_start:send_B_virtual_end, :), send_B_size, &
                               my_local_i_aL, dimen_RI, 0.0_dp, send_ab, send_B_size, &
                               mp2_env%offload_gemm_ctx)
            CALL mp_sendrecv(send_ab, proc_send, external_ab, proc_receive, para_env_sub%group)
            Y_j_aP(:, :) = Y_j_aP + external_ab
            DEALLOCATE (send_ab)

         END DO
         DEALLOCATE (external_ab)

         ! Here, we just use approximate bounds. For large systems virtual(ispin) is approx virtual(jspin), same for B_size
         CALL dgemm_counter_stop(dgemm_counter, 3*virtual(ispin), dimen_RI, my_B_size(jspin))
      ELSE
         CALL dgemm_counter_stop(dgemm_counter, virtual(ispin), dimen_RI, my_B_size(ispin))
      END IF

      IF ((my_i /= my_j) .AND. (.NOT. alpha_beta)) THEN
         ! Alpha-alpha, beta-beta and closed shell
         CALL dgemm_counter_start(dgemm_counter)
         CALL offload_dgemm('T', 'T', my_B_size(ispin), dimen_RI, my_B_size(ispin), 1.0_dp, &
                            t_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), :), my_B_size(ispin), &
                            my_local_i_aL, dimen_RI, 1.0_dp, Y_j_aP, my_B_size(ispin), &
                            mp2_env%offload_gemm_ctx)
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
            proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

            CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

            ALLOCATE (external_ab(dimen_RI, rec_B_size))
            external_ab = 0.0_dp

            CALL mp_sendrecv(my_local_i_aL, proc_send, &
                             external_ab, proc_receive, para_env_sub%group)

            ! Alpha-alpha, beta-beta and closed shell
            CALL offload_dgemm('T', 'T', my_B_size(ispin), dimen_RI, rec_B_size, 1.0_dp, &
                               t_ab(rec_B_virtual_start:rec_B_virtual_end, :), rec_B_size, &
                               external_ab, dimen_RI, 1.0_dp, Y_j_aP, my_B_size(ispin), mp2_env%offload_gemm_ctx)
            DEALLOCATE (external_ab)
         END DO

         CALL dgemm_counter_stop(dgemm_counter, my_B_size(ispin), dimen_RI, virtual(ispin))
      END IF

      CALL timestop(handle)
   END SUBROUTINE mp2_update_P_gamma

! **************************************************************************************************
!> \brief ...
!> \param Gamma_P_ia ...
!> \param ij_index ...
!> \param my_B_size ...
!> \param my_block_size ...
!> \param my_group_L_size ...
!> \param my_i ...
!> \param my_ij_pairs ...
!> \param ngroup ...
!> \param num_integ_group ...
!> \param integ_group_pos2color_sub ...
!> \param num_ij_pairs ...
!> \param proc_map ...
!> \param ij_map ...
!> \param ranges_info_array ...
!> \param Y_i_aP ...
!> \param para_env_exchange ...
!> \param sizes_array ...
!> \param spin ...
! **************************************************************************************************
   SUBROUTINE mp2_redistribute_gamma(Gamma_P_ia, ij_index, my_B_size, &
                                     my_block_size, my_group_L_size, my_i, my_ij_pairs, ngroup, &
                                     num_integ_group, integ_group_pos2color_sub, num_ij_pairs, proc_map, &
                                     ij_map, ranges_info_array, Y_i_aP, para_env_exchange, &
                                     sizes_array, spin)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: Gamma_P_ia
      INTEGER, INTENT(IN)                                :: ij_index, my_B_size, my_block_size, &
                                                            my_group_L_size, my_i, my_ij_pairs, &
                                                            ngroup, num_integ_group
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: integ_group_pos2color_sub, num_ij_pairs, &
                                                            proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: ij_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: ranges_info_array
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Y_i_aP
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_exchange
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sizes_array
      INTEGER, INTENT(IN)                                :: spin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_redistribute_gamma'

      INTEGER :: end_point, handle, handle2, iiB, ij_counter_rec, irep, kkk, lll, Lstart_pos, &
         proc_receive, proc_send, proc_shift, rec_block_size, rec_i, rec_ij_index, send_L_size, &
         start_point, tag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BI_C_rec, BI_C_send

! In alpha-beta case Y_i_aP_beta is sent as Y_j_aP

      CALL timeset(routineN//"_comm2", handle)

      tag = 43

      IF (ij_index <= my_ij_pairs) THEN
         ! somethig to send
         ! start with myself
         CALL timeset(routineN//"_comm2_w", handle2)
         DO irep = 0, num_integ_group - 1
            Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
            start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
            end_point = ranges_info_array(4, irep, para_env_exchange%mepos)
!$OMP       PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll,iiB) &
!$OMP                                 SHARED(start_point,end_point,Lstart_pos,my_block_size,&
!$OMP                                        Gamma_P_ia,my_i,my_B_size,Y_i_aP)
            DO kkk = start_point, end_point
               lll = kkk - start_point + Lstart_pos
               DO iiB = 1, my_block_size
                  Gamma_P_ia(1:my_B_size, my_i + iiB - 1, kkk) = &
                     Gamma_P_ia(1:my_B_size, my_i + iiB - 1, kkk) + &
                     Y_i_aP(1:my_B_size, lll, iiB)
               END DO
            END DO
!$OMP              END PARALLEL DO
         END DO
         CALL timestop(handle2)

         ! Y_i_aP(my_B_size,dimen_RI,block_size)

         DO proc_shift = 1, para_env_exchange%num_pe - 1
            proc_send = proc_map(para_env_exchange%mepos + proc_shift)
            proc_receive = proc_map(para_env_exchange%mepos - proc_shift)

            send_L_size = sizes_array(proc_send)
            ALLOCATE (BI_C_send(my_B_size, my_block_size, send_L_size))
            CALL timeset(routineN//"_comm2_w", handle2)
            BI_C_send = 0.0_dp
            DO irep = 0, num_integ_group - 1
               Lstart_pos = ranges_info_array(1, irep, proc_send)
               start_point = ranges_info_array(3, irep, proc_send)
               end_point = ranges_info_array(4, irep, proc_send)
!$OMP             PARALLEL DO DEFAULT(NONE) PRIVATE(kkk,lll,iiB) &
!$OMP                                       SHARED(start_point,end_point,Lstart_pos,my_block_size,&
!$OMP                                              BI_C_send,my_B_size,Y_i_aP)
               DO kkk = start_point, end_point
                  lll = kkk - start_point + Lstart_pos
                  DO iiB = 1, my_block_size
                     BI_C_send(1:my_B_size, iiB, kkk) = Y_i_aP(1:my_B_size, lll, iiB)
                  END DO
               END DO
!$OMP                END PARALLEL DO
            END DO
            CALL timestop(handle2)

            rec_ij_index = num_ij_pairs(proc_receive)

            IF (ij_index <= rec_ij_index) THEN
               ! we know that proc_receive has something to send for us, let's see what
               ij_counter_rec = &
                  (ij_index - MIN(1, integ_group_pos2color_sub(proc_receive)))*ngroup + integ_group_pos2color_sub(proc_receive)

               rec_i = ij_map(spin, ij_counter_rec)
               rec_block_size = ij_map(3, ij_counter_rec)

               ALLOCATE (BI_C_rec(my_B_size, rec_block_size, my_group_L_size))
               BI_C_rec = 0.0_dp

               CALL mp_sendrecv(BI_C_send, proc_send, &
                                BI_C_rec, proc_receive, &
                                para_env_exchange%group, tag)

               CALL timeset(routineN//"_comm2_w", handle2)
               DO irep = 0, num_integ_group - 1
                  start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
                  end_point = ranges_info_array(4, irep, para_env_exchange%mepos)
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE) SHARED(start_point,end_point,rec_block_size,&
!$OMP                                           Gamma_P_ia,rec_i,iiB,my_B_size,BI_C_rec)
                  Gamma_P_ia(:, rec_i:rec_i + rec_block_size - 1, start_point:end_point) = &
                     Gamma_P_ia(:, rec_i:rec_i + rec_block_size - 1, start_point:end_point) + &
                     BI_C_rec(1:my_B_size, :, start_point:end_point)
!$OMP             END PARALLEL WORKSHARE
               END DO
               CALL timestop(handle2)

               DEALLOCATE (BI_C_rec)

            ELSE
               ! we have something to send but nothing to receive
               CALL mp_send(BI_C_send, proc_send, tag, para_env_exchange%group)

            END IF

            DEALLOCATE (BI_C_send)
         END DO

      ELSE
         ! noting to send check if we have to receive
         DO proc_shift = 1, para_env_exchange%num_pe - 1
            proc_send = proc_map(para_env_exchange%mepos + proc_shift)
            proc_receive = proc_map(para_env_exchange%mepos - proc_shift)
            rec_ij_index = num_ij_pairs(proc_receive)

            IF (ij_index <= rec_ij_index) THEN
               ! we know that proc_receive has something to send for us, let's see what
               ij_counter_rec = &
                  (ij_index - MIN(1, integ_group_pos2color_sub(proc_receive)))*ngroup + integ_group_pos2color_sub(proc_receive)

               rec_i = ij_map(spin, ij_counter_rec)
               rec_block_size = ij_map(3, ij_counter_rec)

               ALLOCATE (BI_C_rec(my_B_size, rec_block_size, my_group_L_size))

               BI_C_rec = 0.0_dp

               CALL mp_recv(BI_C_rec, proc_receive, tag, para_env_exchange%group)

               CALL timeset(routineN//"_comm2_w", handle2)
               DO irep = 0, num_integ_group - 1
                  start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
                  end_point = ranges_info_array(4, irep, para_env_exchange%mepos)
!$OMP             PARALLEL WORKSHARE DEFAULT(NONE) SHARED(start_point,end_point,rec_block_size,&
!$OMP                                           Gamma_P_ia,rec_i,iiB,my_B_size,BI_C_rec)
                  Gamma_P_ia(:, rec_i:rec_i + rec_block_size - 1, start_point:end_point) = &
                     Gamma_P_ia(:, rec_i:rec_i + rec_block_size - 1, start_point:end_point) + &
                     BI_C_rec(1:my_B_size, :, start_point:end_point)
!$OMP             END PARALLEL WORKSHARE
               END DO
               CALL timestop(handle2)

               DEALLOCATE (BI_C_rec)

            END IF
         END DO

      END IF
      CALL timestop(handle)

   END SUBROUTINE mp2_redistribute_gamma

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param open_shell ...
!> \param beta_beta ...
!> \param Bib_C ...
!> \param unit_nr ...
!> \param dimen_RI ...
!> \param my_B_size ...
!> \param ngroup ...
!> \param num_integ_group ...
!> \param my_group_L_size ...
!> \param color_sub ...
!> \param ranges_info_array ...
!> \param para_env_exchange ...
!> \param para_env_sub ...
!> \param proc_map ...
!> \param my_B_virtual_start ...
!> \param my_B_virtual_end ...
!> \param sizes_array ...
!> \param gd_B_virtual ...
!> \param sub_proc_map ...
!> \param integ_group_pos2color_sub ...
!> \param dgemm_counter ...
! **************************************************************************************************
   SUBROUTINE quasi_degenerate_P_ij(mp2_env, Eigenval, homo, virtual, open_shell, &
                                    beta_beta, Bib_C, unit_nr, dimen_RI, &
                                    my_B_size, ngroup, num_integ_group, my_group_L_size, &
                                    color_sub, ranges_info_array, para_env_exchange, para_env_sub, proc_map, &
                                    my_B_virtual_start, my_B_virtual_end, sizes_array, gd_B_virtual, &
                                    sub_proc_map, integ_group_pos2color_sub, dgemm_counter)
      TYPE(mp2_type)                                     :: mp2_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      LOGICAL, INTENT(IN)                                :: open_shell, beta_beta
      TYPE(three_dim_real_array), DIMENSION(:), &
         INTENT(IN)                                      :: BIb_C
      INTEGER, INTENT(IN)                                :: unit_nr, dimen_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_B_size
      INTEGER, INTENT(IN)                                :: ngroup, num_integ_group, &
                                                            my_group_L_size, color_sub
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: ranges_info_array
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_exchange, para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: proc_map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_B_virtual_start, my_B_virtual_end
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sizes_array
      TYPE(group_dist_d1_type), DIMENSION(:), INTENT(IN) :: gd_B_virtual
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: sub_proc_map, integ_group_pos2color_sub
      TYPE(dgemm_counter_type), INTENT(INOUT)            :: dgemm_counter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'quasi_degenerate_P_ij'

      INTEGER :: a, a_global, b, b_global, end_point, handle, handle2, ijk, ijk_counter, &
         ijk_counter_send, ijk_index, irep, ispin, kspin, Lend_pos, Lstart_pos, my_i, my_j, my_k, &
         my_virtual, nspins, proc_receive, proc_send, proc_shift, rec_B_size, rec_B_virtual_end, &
         rec_B_virtual_start, rec_L_size, send_B_size, send_B_virtual_end, send_B_virtual_start, &
         send_i, send_ijk_index, send_j, send_k, size_B_i, size_B_ij, size_B_ijk, size_B_k, &
         start_point, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: max_ijk, my_ijk
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: num_ijk
      LOGICAL                                            :: alpha_beta
      REAL(KIND=dp)                                      :: amp_fac, P_ij_elem
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         TARGET                                          :: BI_C_rec, external_ab, external_aL, &
                                                            local_ab, local_aL, t_ab
      TYPE(two_dim_int_array), ALLOCATABLE, DIMENSION(:) :: ijk_map

      CALL timeset(routineN//"_ij_sing", handle)
! Define the number of loops over orbital triplets

      tag = 44

      nspins = SIZE(BIb_C)
      alpha_beta = (nspins == 2)

      ! Find the number of quasi-degenerate orbitals and orbital triplets

      CALL Find_quasi_degenerate_ij(my_ijk, homo, Eigenval, mp2_env, ijk_map, unit_nr, ngroup, &
                                    beta_beta, para_env_exchange, num_ijk, max_ijk, color_sub)

      ! Set amplitude factor
      amp_fac = mp2_env%scale_S + mp2_env%scale_T
      IF (open_shell) amp_fac = mp2_env%scale_T

      ! Loop(s) over orbital triplets
      DO ispin = 1, nspins
         size_B_i = my_B_size(ispin)
         IF (ispin == 1 .AND. alpha_beta) THEN
            kspin = 2
         ELSE
            kspin = 1
         END IF
         size_B_k = my_B_size(kspin)
         my_virtual = virtual(ispin)
         size_B_ij = 2*size_B_i
         size_B_ijk = size_B_ij + size_B_k

         ALLOCATE (local_aL(dimen_RI, size_B_ijk))
         ALLOCATE (t_ab(my_virtual, size_B_k))

         DO ijk_index = 1, max_ijk(ispin)
            ijk = my_ijk(ispin)
            IF (ijk_index <= ijk) THEN
               ! work to be done
               ijk_counter = (ijk_index - MIN(1, color_sub))*ngroup + color_sub
               my_i = ijk_map(ispin)%array(ijk_counter, 1)
               my_j = ijk_map(ispin)%array(ijk_counter, 2)
               my_k = ijk_map(ispin)%array(ijk_counter, 3)

               local_aL = 0.0_dp
               DO irep = 0, num_integ_group - 1
                  Lstart_pos = ranges_info_array(1, irep, para_env_exchange%mepos)
                  Lend_pos = ranges_info_array(2, irep, para_env_exchange%mepos)
                  start_point = ranges_info_array(3, irep, para_env_exchange%mepos)
                  end_point = ranges_info_array(4, irep, para_env_exchange%mepos)

!$OMP PARALLEL DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_i,size_B_ij,start_point,end_point,&
!$OMP                               my_i,my_j,my_k,BIb_C,ispin,kspin,size_B_k)
!$OMP WORKSHARE
                  local_aL(Lstart_pos:Lend_pos, :size_B_i) = BIb_C(ispin)%array(start_point:end_point, 1:size_B_i, my_i)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
                  local_aL(Lstart_pos:Lend_pos, size_B_i + 1:size_B_ij) = &
                     BIb_C(ispin)%array(start_point:end_point, 1:size_B_i, my_j)
!$OMP END WORKSHARE NOWAIT
!$OMP WORKSHARE
                  local_aL(Lstart_pos:Lend_pos, size_B_ij + 1:) = BIb_C(kspin)%array(start_point:end_point, 1:size_B_k, my_k)
!$OMP END WORKSHARE NOWAIT
!$OMP END PARALLEL
               END DO

               CALL timeset(routineN//"_comm", handle2)
               DO proc_shift = 1, para_env_exchange%num_pe - 1
                  proc_send = proc_map(para_env_exchange%mepos + proc_shift)
                  proc_receive = proc_map(para_env_exchange%mepos - proc_shift)

                  send_ijk_index = num_ijk(proc_send, ispin)

                  rec_L_size = sizes_array(proc_receive)
                  ALLOCATE (BI_C_rec(rec_L_size, size_B_i))

                  IF (ijk_index <= send_ijk_index) THEN
                     ! something to send
                     ijk_counter_send = (ijk_index - MIN(1, integ_group_pos2color_sub(proc_send)))* &
                                        ngroup + integ_group_pos2color_sub(proc_send)
                     send_i = ijk_map(ispin)%array(ijk_counter_send, 1)
                     send_j = ijk_map(ispin)%array(ijk_counter_send, 2)
                     send_k = ijk_map(ispin)%array(ijk_counter_send, 3)
                  END IF

                  ! occupied i
                  BI_C_rec = 0.0_dp
                  IF (ijk_index <= send_ijk_index) THEN
                     CALL mp_sendrecv(BIb_C(ispin)%array(1:my_group_L_size, 1:size_B_i, send_i), proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:size_B_i), proc_receive, &
                                      para_env_exchange%group, tag)
                  ELSE
                     ! nothing to send
                     CALL mp_recv(BI_C_rec(1:rec_L_size, 1:size_B_i), proc_receive, tag, &
                                  para_env_exchange%group)
                  END IF
                  DO irep = 0, num_integ_group - 1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_i,BI_C_rec,start_point,end_point)
                     local_aL(Lstart_pos:Lend_pos, 1:size_B_i) = BI_C_rec(start_point:end_point, 1:size_B_i)
!$OMP END PARALLEL WORKSHARE
                  END DO
                  IF (my_i == my_k .AND. .NOT. alpha_beta) THEN
                     DO irep = 0, num_integ_group - 1
                        Lstart_pos = ranges_info_array(1, irep, proc_receive)
                        Lend_pos = ranges_info_array(2, irep, proc_receive)
                        start_point = ranges_info_array(3, irep, proc_receive)
                        end_point = ranges_info_array(4, irep, proc_receive)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_ij,BI_C_rec,start_point,end_point,size_B_k)
                        local_aL(Lstart_pos:Lend_pos, size_B_ij + 1:) = BI_C_rec(start_point:end_point, 1:size_B_k)
!$OMP END PARALLEL WORKSHARE
                     END DO
                  END IF

                  ! occupied j
                  BI_C_rec = 0.0_dp
                  IF (ijk_index <= send_ijk_index) THEN
                     CALL mp_sendrecv(BIb_C(ispin)%array(1:my_group_L_size, 1:size_B_i, send_j), proc_send, &
                                      BI_C_rec(1:rec_L_size, 1:size_B_i), proc_receive, &
                                      para_env_exchange%group, tag)
                  ELSE
                     ! nothing to send
                     CALL mp_recv(BI_C_rec(1:rec_L_size, 1:size_B_i), proc_receive, tag, &
                                  para_env_exchange%group)
                  END IF
                  DO irep = 0, num_integ_group - 1
                     Lstart_pos = ranges_info_array(1, irep, proc_receive)
                     Lend_pos = ranges_info_array(2, irep, proc_receive)
                     start_point = ranges_info_array(3, irep, proc_receive)
                     end_point = ranges_info_array(4, irep, proc_receive)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_ij,BI_C_rec,start_point,end_point,size_B_i)
                     local_aL(Lstart_pos:Lend_pos, size_B_i + 1:size_B_ij) = BI_C_rec(start_point:end_point, 1:size_B_i)
!$OMP END PARALLEL WORKSHARE
                  END DO
                  IF (my_j == my_k .AND. .NOT. alpha_beta) THEN
                     DO irep = 0, num_integ_group - 1
                        Lstart_pos = ranges_info_array(1, irep, proc_receive)
                        Lend_pos = ranges_info_array(2, irep, proc_receive)
                        start_point = ranges_info_array(3, irep, proc_receive)
                        end_point = ranges_info_array(4, irep, proc_receive)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_ij,BI_C_rec,start_point,end_point,size_B_k)
                        local_aL(Lstart_pos:Lend_pos, size_B_ij + 1:) = BI_C_rec(start_point:end_point, 1:size_B_k)
!$OMP END PARALLEL WORKSHARE
                     END DO
                  END IF

                  IF ((my_i /= my_k .AND. my_j /= my_k) .OR. alpha_beta) THEN
                     IF ((send_i /= send_k .AND. send_j /= send_k) .OR. alpha_beta) THEN
                        ! occupied k
                        BI_C_rec = 0.0_dp
                        DEALLOCATE (BI_C_rec)
                        ALLOCATE (BI_C_rec(rec_L_size, size_B_k))
                        BI_C_rec = 0.0_dp
                        IF (ijk_index <= send_ijk_index) THEN
                           CALL mp_sendrecv(BIb_C(kspin)%array(1:my_group_L_size, 1:size_B_k, send_k), proc_send, &
                                            BI_C_rec(1:rec_L_size, 1:size_B_k), proc_receive, &
                                            para_env_exchange%group, tag)
                        ELSE
                           ! nothing to send
                           BI_C_rec = 0.0_dp
                           CALL mp_recv(BI_C_rec(1:rec_L_size, 1:size_B_k), proc_receive, tag, &
                                        para_env_exchange%group)
                        END IF
                        DO irep = 0, num_integ_group - 1
                           Lstart_pos = ranges_info_array(1, irep, proc_receive)
                           Lend_pos = ranges_info_array(2, irep, proc_receive)
                           start_point = ranges_info_array(3, irep, proc_receive)
                           end_point = ranges_info_array(4, irep, proc_receive)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_ij,BI_C_rec,start_point,end_point,size_B_k)
                           local_aL(Lstart_pos:Lend_pos, size_B_ij + 1:) = BI_C_rec(start_point:end_point, 1:size_B_k)
!$OMP END PARALLEL WORKSHARE
                        END DO
                     ELSE IF (.NOT. alpha_beta) THEN
                        ! occupied k
                        BI_C_rec = 0.0_dp
                        DEALLOCATE (BI_C_rec)
                        ALLOCATE (BI_C_rec(rec_L_size, size_B_k))
                        BI_C_rec = 0.0_dp
                        BI_C_rec = 0.0_dp
                        CALL mp_recv(BI_C_rec(1:rec_L_size, 1:size_B_k), proc_receive, tag, &
                                     para_env_exchange%group)
                        DO irep = 0, num_integ_group - 1
                           Lstart_pos = ranges_info_array(1, irep, proc_receive)
                           Lend_pos = ranges_info_array(2, irep, proc_receive)
                           start_point = ranges_info_array(3, irep, proc_receive)
                           end_point = ranges_info_array(4, irep, proc_receive)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE) SHARED(local_aL,Lstart_pos,Lend_pos,size_B_ij,BI_C_rec,start_point,end_point,size_B_k)
                           local_aL(Lstart_pos:Lend_pos, size_B_ij + 1:) = BI_C_rec(start_point:end_point, 1:size_B_k)
!$OMP END PARALLEL WORKSHARE
                        END DO
                     END IF
                  ELSE IF (send_i /= send_k .AND. send_j /= send_k .AND. .NOT. alpha_beta .AND. ijk_index <= send_ijk_index) THEN
                     CALL mp_send(BIb_C(kspin)%array(1:my_group_L_size, 1:size_B_k, send_k), proc_send, &
                                  tag, para_env_exchange%group)
                  END IF

                  DEALLOCATE (BI_C_rec)
               END DO
               CALL timestop(handle2)

               ! expand integrals
               CALL timeset(routineN//"_exp_ik", handle2)
               CALL dgemm_counter_start(dgemm_counter)
               ALLOCATE (local_ab(my_virtual, size_B_k))
               local_ab = 0.0_dp
               CALL offload_dgemm('T', 'N', size_B_i, size_B_k, dimen_RI, 1.0_dp, &
                                  local_aL(:, :size_B_i), dimen_RI, local_aL(:, size_B_ij + 1:), dimen_RI, &
                                  0.0_dp, local_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), 1:size_B_k), size_B_i, &
                                  mp2_env%offload_gemm_ctx)
               DO proc_shift = 1, para_env_sub%num_pe - 1
                  proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
                  proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

                  CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

                  ALLOCATE (external_aL(dimen_RI, rec_B_size))
                  external_aL = 0.0_dp

                  CALL mp_sendrecv(local_aL(:, :size_B_i), proc_send, &
                                   external_aL, proc_receive, &
                                   para_env_sub%group, tag)

                  CALL offload_dgemm('T', 'N', rec_B_size, size_B_k, dimen_RI, 1.0_dp, &
                                     external_aL, dimen_RI, local_aL(:, size_B_ij + 1:), dimen_RI, &
                                     0.0_dp, local_ab(rec_B_virtual_start:rec_B_virtual_end, 1:size_B_k), rec_B_size, &
                                     mp2_env%offload_gemm_ctx)
                  DEALLOCATE (external_aL)
               END DO
               CALL dgemm_counter_stop(dgemm_counter, my_virtual, size_B_k, dimen_RI)
               CALL timestop(handle2)

               ! Amplitudes
               CALL timeset(routineN//"_tab", handle2)
               t_ab = 0.0_dp
               ! Alpha-alpha, beta-beta and closed shell
               IF (.NOT. alpha_beta) THEN
                  DO b = 1, size_B_k
                     b_global = b + my_B_virtual_start(1) - 1
                     DO a = 1, my_B_size(1)
                        a_global = a + my_B_virtual_start(1) - 1
                        t_ab(a_global, b) = (amp_fac*local_ab(a_global, b) - mp2_env%scale_T*local_ab(b_global, a))/ &
                                            (Eigenval(my_i, 1) + Eigenval(my_k, 1) &
                                             - Eigenval(homo(1) + a_global, 1) - Eigenval(homo(1) + b_global, 1))
                     END DO
                  END DO
               ELSE
                  DO b = 1, size_B_k
                     b_global = b + my_B_virtual_start(kspin) - 1
                     DO a = 1, my_B_size(ispin)
                        a_global = a + my_B_virtual_start(ispin) - 1
                        t_ab(a_global, b) = mp2_env%scale_S*local_ab(a_global, b)/ &
                                            (Eigenval(my_i, ispin) + Eigenval(my_k, kspin) &
                                             - Eigenval(homo(ispin) + a_global, ispin) - Eigenval(homo(kspin) + b_global, kspin))
                     END DO
                  END DO
               END IF

               IF (.NOT. alpha_beta) THEN
                  DO proc_shift = 1, para_env_sub%num_pe - 1
                     proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
                     proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)
                     CALL get_group_dist(gd_B_virtual(1), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)
                     CALL get_group_dist(gd_B_virtual(1), proc_send, send_B_virtual_start, send_B_virtual_end, send_B_size)

                     ALLOCATE (external_ab(size_B_i, rec_B_size))
                     external_ab = 0.0_dp
                     CALL mp_sendrecv(local_ab(send_B_virtual_start:send_B_virtual_end, 1:size_B_k), proc_send, &
                                      external_ab(1:size_B_i, 1:rec_B_size), proc_receive, para_env_sub%group, tag)

                     DO b = 1, my_B_size(1)
                        b_global = b + my_B_virtual_start(1) - 1
                        DO a = 1, rec_B_size
                           a_global = a + rec_B_virtual_start - 1
                           t_ab(a_global, b) = (amp_fac*local_ab(a_global, b) - mp2_env%scale_T*external_ab(b, a))/ &
                                               (Eigenval(my_i, 1) + Eigenval(my_k, 1) &
                                                - Eigenval(homo(1) + a_global, 1) - Eigenval(homo(1) + b_global, 1))
                        END DO
                     END DO

                     DEALLOCATE (external_ab)
                  END DO
               END IF
               CALL timestop(handle2)

               ! Expand the second set of integrals
               CALL timeset(routineN//"_exp_jk", handle2)
               local_ab = 0.0_dp
               CALL dgemm_counter_start(dgemm_counter)
               CALL offload_dgemm('T', 'N', size_B_i, size_B_k, dimen_RI, 1.0_dp, &
                                  local_aL(:, size_B_i + 1:size_B_ij), dimen_RI, local_aL(:, size_B_ij + 1:), dimen_RI, &
                                  0.0_dp, local_ab(my_B_virtual_start(ispin):my_B_virtual_end(ispin), 1:size_B_k), size_B_i, &
                                  mp2_env%offload_gemm_ctx)
               DO proc_shift = 1, para_env_sub%num_pe - 1
                  proc_send = sub_proc_map(para_env_sub%mepos + proc_shift)
                  proc_receive = sub_proc_map(para_env_sub%mepos - proc_shift)

                  CALL get_group_dist(gd_B_virtual(ispin), proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

                  ALLOCATE (external_aL(dimen_RI, rec_B_size))
                  external_aL = 0.0_dp

                  CALL mp_sendrecv(local_aL(:, size_B_i + 1:size_B_ij), proc_send, &
                                   external_aL, proc_receive, &
                                   para_env_sub%group, tag)
                  CALL offload_dgemm('T', 'N', rec_B_size, size_B_k, dimen_RI, 1.0_dp, &
                                     external_aL, dimen_RI, local_aL(:, size_B_ij + 1:), dimen_RI, &
                                     0.0_dp, local_ab(rec_B_virtual_start:rec_B_virtual_end, 1:size_B_k), rec_B_size, &
                                     mp2_env%offload_gemm_ctx)
                  DEALLOCATE (external_aL)
               END DO
               CALL dgemm_counter_stop(dgemm_counter, my_virtual, size_B_k, dimen_RI)
               CALL timestop(handle2)

               CALL timeset(routineN//"_Pij", handle2)
               DO b = 1, size_B_k
                  b_global = b + my_B_virtual_start(kspin) - 1
                  DO a = 1, my_B_size(ispin)
                     a_global = a + my_B_virtual_start(ispin) - 1
                     local_ab(a_global, b) = &
                        local_ab(a_global, b)/(Eigenval(my_j, ispin) + Eigenval(my_k, kspin) &
                                               - Eigenval(homo(ispin) + a_global, ispin) - Eigenval(homo(kspin) + b_global, kspin))
                  END DO
               END DO
               !
               P_ij_elem = SUM(local_ab*t_ab)
               DEALLOCATE (local_ab)
               IF ((.NOT. open_shell) .AND. (.NOT. alpha_beta)) THEN
                  P_ij_elem = P_ij_elem*2.0_dp
               END IF
               IF (beta_beta) THEN
                  mp2_env%ri_grad%P_ij(2)%array(my_i, my_j) = mp2_env%ri_grad%P_ij(2)%array(my_i, my_j) - P_ij_elem
                  mp2_env%ri_grad%P_ij(2)%array(my_j, my_i) = mp2_env%ri_grad%P_ij(2)%array(my_j, my_i) - P_ij_elem
               ELSE
                  mp2_env%ri_grad%P_ij(ispin)%array(my_i, my_j) = mp2_env%ri_grad%P_ij(ispin)%array(my_i, my_j) - P_ij_elem
                  mp2_env%ri_grad%P_ij(ispin)%array(my_j, my_i) = mp2_env%ri_grad%P_ij(ispin)%array(my_j, my_i) - P_ij_elem
               END IF
               CALL timestop(handle2)
            ELSE
               CALL timeset(routineN//"_comm", handle2)
               ! no work to be done, possible messeges to be exchanged
               DO proc_shift = 1, para_env_exchange%num_pe - 1
                  proc_send = proc_map(para_env_exchange%mepos + proc_shift)
                  proc_receive = proc_map(para_env_exchange%mepos - proc_shift)

                  send_ijk_index = num_ijk(proc_send, ispin)

                  IF (ijk_index <= send_ijk_index) THEN
                     ! somethig to send
                     ijk_counter_send = (ijk_index - MIN(1, integ_group_pos2color_sub(proc_send)))*ngroup + &
                                        integ_group_pos2color_sub(proc_send)
                     send_i = ijk_map(ispin)%array(ijk_counter_send, 1)
                     send_j = ijk_map(ispin)%array(ijk_counter_send, 2)
                     send_k = ijk_map(ispin)%array(ijk_counter_send, 3)
                     ! occupied i
                     CALL mp_send(BIb_C(ispin)%array(1:my_group_L_size, 1:size_B_i, send_i), proc_send, tag, &
                                  para_env_exchange%group)
                     ! occupied j
                     CALL mp_send(BIb_C(ispin)%array(1:my_group_L_size, 1:size_B_i, send_j), proc_send, tag, &
                                  para_env_exchange%group)
                     IF ((send_i /= send_k .AND. send_j /= send_k) .OR. alpha_beta) THEN
                        ! occupied k
                        CALL mp_send(BIb_C(kspin)%array(1:my_group_L_size, 1:size_B_k, send_k), proc_send, tag, &
                                     para_env_exchange%group)
                     END IF
                  END IF

               END DO ! proc loop
               CALL timestop(handle2)
            END IF
         END DO ! ijk_index loop
         DEALLOCATE (local_aL)
         DEALLOCATE (t_ab)
      END DO ! over number of loops (iloop)
      !
      DO ispin = 1, SIZE(ijk_map)
         DEALLOCATE (ijk_map(ispin)%array)
      END DO
      DEALLOCATE (ijk_map)
      DEALLOCATE (num_ijk)
      CALL timestop(handle)

   END SUBROUTINE Quasi_degenerate_P_ij

! **************************************************************************************************
!> \brief ...
!> \param my_ijk ...
!> \param homo ...
!> \param Eigenval ...
!> \param mp2_env ...
!> \param ijk_map ...
!> \param unit_nr ...
!> \param ngroup ...
!> \param beta_beta ...
!> \param para_env_exchange ...
!> \param num_ijk ...
!> \param max_ijk ...
!> \param color_sub ...
! **************************************************************************************************
   SUBROUTINE Find_quasi_degenerate_ij(my_ijk, homo, Eigenval, mp2_env, ijk_map, unit_nr, ngroup, &
                                       beta_beta, para_env_exchange, num_ijk, max_ijk, color_sub)

      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: my_ijk
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      TYPE(two_dim_int_array), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: ijk_map
      INTEGER, INTENT(IN)                                :: unit_nr, ngroup
      LOGICAL, INTENT(IN)                                :: beta_beta
      TYPE(cp_para_env_type), INTENT(IN)                 :: para_env_exchange
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: num_ijk
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: max_ijk
      INTEGER, INTENT(IN)                                :: color_sub

      INTEGER                                            :: iib, ijk_counter, ispin, jjb, kkb, &
                                                            my_homo, nspins, num_sing_ij, total_ijk
      LOGICAL                                            :: alpha_beta

      nspins = SIZE(homo)
      alpha_beta = (nspins == 2)

      ALLOCATE (ijk_map(nspins))
      ALLOCATE (my_ijk(nspins))
      ALLOCATE (num_ijk(0:para_env_exchange%num_pe - 1, nspins))
      ALLOCATE (max_ijk(nspins))

      ! General case
      DO ispin = 1, nspins
         my_homo = homo(nspins)
         IF (ispin == 2) my_homo = homo(1)

         num_sing_ij = 0
         DO iiB = 1, homo(ispin)
            ! diagonal elements already updated
            DO jjB = iiB + 1, homo(ispin)
               IF (ABS(Eigenval(jjB, ispin) - Eigenval(iiB, ispin)) < mp2_env%ri_mp2%eps_canonical) &
                  num_sing_ij = num_sing_ij + 1
            END DO
         END DO

         IF (unit_nr > 0) THEN
         IF (.NOT. beta_beta .AND. ispin /= 2) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MO_INFO| Number of ij pairs below EPS_CANONICAL:", num_sing_ij
         ELSE
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MO_INFO| Number of ij pairs (spin beta) below EPS_CANONICAL:", num_sing_ij
         END IF
         CALL m_flush(unit_nr)
         END IF
         total_ijk = my_homo*num_sing_ij
         ALLOCATE (ijk_map(ispin)%array(total_ijk, 3))
         ijk_map(ispin)%array = 0

         my_ijk(ispin) = 0
         ijk_counter = 0
         DO iiB = 1, homo(ispin)
            ! diagonal elements already updated
            DO jjB = iiB + 1, homo(ispin)
               IF (ABS(Eigenval(jjB, ispin) - Eigenval(iiB, ispin)) >= mp2_env%ri_mp2%eps_canonical) CYCLE
               DO kkB = 1, my_homo
                  ijk_counter = ijk_counter + 1
                  ijk_map(ispin)%array(ijk_counter, 1) = iiB
                  ijk_map(ispin)%array(ijk_counter, 2) = jjB
                  ijk_map(ispin)%array(ijk_counter, 3) = kkB
                  IF (MOD(ijk_counter, ngroup) == color_sub) my_ijk(ispin) = my_ijk(ispin) + 1
               END DO
            END DO
         END DO

         CALL mp_allgather(my_ijk(ispin), num_ijk(:, ispin), para_env_exchange%group)
         max_ijk(ispin) = MAXVAL(num_ijk(:, ispin))
      END DO

   END SUBROUTINE Find_quasi_degenerate_ij

END MODULE mp2_ri_gpw
