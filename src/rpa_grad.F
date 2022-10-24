!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate RI-RPA and SOS-MP2 gradients
!> \par History
!>      10.2021 created [Frederick Stein]
! **************************************************************************************************
MODULE rpa_grad
   USE ISO_C_BINDING,                   ONLY: C_NULL_PTR,&
                                              c_associated
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_fm_basic_linalg,              ONLY: cp_fm_geadd,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_indxg2l, cp_fm_indxg2p, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_to_fm, cp_fm_to_fm_submat_general, cp_fm_type
   USE cp_para_env,                     ONLY: cp_para_env_release,&
                                              cp_para_env_split
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dgemm_counter_types,             ONLY: dgemm_counter_start,&
                                              dgemm_counter_stop,&
                                              dgemm_counter_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              group_dist_proc,&
                                              maxsize,&
                                              release_group_dist
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: compare_potential_types
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU,&
                                              local_gemm,&
                                              local_gemm_create,&
                                              local_gemm_destroy,&
                                              local_gemm_set_op_threshold_gpu
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_allgather,&
                                              mp_alltoall,&
                                              mp_recv,&
                                              mp_send,&
                                              mp_sendrecv,&
                                              mp_sum
   USE mp2_laplace,                     ONLY: calc_fm_mat_s_laplace
   USE mp2_ri_grad_util,                ONLY: array2fm,&
                                              create_dbcsr_gamma,&
                                              fm2array,&
                                              prepare_redistribution
   USE mp2_types,                       ONLY: mp2_type,&
                                              one_dim_int_array,&
                                              two_dim_int_array,&
                                              two_dim_real_array
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rpa_util,                        ONLY: calc_fm_mat_S_rpa,&
                                              remove_scaling_factor_rpa
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_grad'

   PUBLIC :: rpa_grad_needed_mem, rpa_grad_type, rpa_grad_create, rpa_grad_finalize, rpa_grad_matrix_operations, rpa_grad_copy_Q

   TYPE sos_mp2_grad_work_type
      PRIVATE
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: pair_list
      TYPE(one_dim_int_array), DIMENSION(:), ALLOCATABLE :: index2send, index2recv
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: P
   END TYPE

   TYPE rpa_grad_work_type
      TYPE(cp_fm_type) :: fm_mat_Q_copy
      TYPE(one_dim_int_array), DIMENSION(:, :), ALLOCATABLE :: index2send
      TYPE(two_dim_int_array), DIMENSION(:, :), ALLOCATABLE :: index2recv
      TYPE(group_dist_d1_type), DIMENSION(:), ALLOCATABLE :: gd_homo, gd_virtual
      INTEGER, DIMENSION(2) :: grid, mepos
      TYPE(two_dim_real_array), DIMENSION(:), ALLOCATABLE :: P_ij, P_ab
   END TYPE

   TYPE rpa_grad_type
      PRIVATE
      TYPE(cp_fm_type) :: fm_Gamma_PQ
      TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE :: fm_Y
      TYPE(sos_mp2_grad_work_type), ALLOCATABLE, DIMENSION(:) :: sos_mp2_work_occ, sos_mp2_work_virt
      TYPE(rpa_grad_work_type) :: rpa_work
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Calculates the necessary minimum memory for the Gradient code ion MiB
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param mem_per_rank ...
!> \param mem_per_repl ...
!> \param do_ri_sos_laplace_mp2 ...
!> \return ...
! **************************************************************************************************
   PURE SUBROUTINE rpa_grad_needed_mem(homo, virtual, dimen_RI, mem_per_rank, mem_per_repl, do_ri_sos_laplace_mp2)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI
      REAL(KIND=dp), INTENT(INOUT)                       :: mem_per_rank, mem_per_repl
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2

      REAL(KIND=dp)                                      :: mem_iaK, mem_KL, mem_pab, mem_pij

      mem_iaK = SUM(REAL(virtual, KIND=dp)*homo)*dimen_RI
      mem_pij = SUM(REAL(homo, KIND=dp)**2)
      mem_pab = SUM(REAL(virtual, KIND=dp)**2)
      mem_KL = REAL(dimen_RI, KIND=dp)*dimen_RI

      ! Required matrices iaK
      ! Ytot_iaP = sum_tau Y_iaP(tau)
      ! Y_iaP(tau) = S_iaP(tau)*Q_PQ(tau) (work array)
      ! Required matrices density matrices
      ! Pij (local)
      ! Pab (local)
      ! Additionally with SOS-MP2
      ! Send and receive buffers for degenerate orbital pairs (rough estimate: everything)
      ! Additionally with RPA
      ! copy of work matrix
      ! receive buffer for calculation of density matrix
      ! copy of matrix Q
      mem_per_rank = mem_per_rank + (mem_pij + mem_pab)*8.0_dp/(1024**2)
      mem_per_repl = mem_per_repl + (mem_iaK + 2.0_dp*mem_iaK/SIZE(homo) + mem_KL)*8.0_dp/(1024**2)
      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         mem_per_repl = mem_per_rank + (mem_iaK/SIZE(homo) + mem_KL)*8.0_dp/(1024**2)
      END IF

   END SUBROUTINE rpa_grad_needed_mem

! **************************************************************************************************
!> \brief Creates the arrays of a rpa_grad_type
!> \param rpa_grad ...
!> \param fm_mat_Q ...
!> \param fm_mat_S ...
!> \param homo ...
!> \param virtual ...
!> \param mp2_env ...
!> \param Eigenval ...
!> \param unit_nr ...
!> \param do_ri_sos_laplace_mp2 ...
! **************************************************************************************************
   SUBROUTINE rpa_grad_create(rpa_grad, fm_mat_Q, fm_mat_S, &
                              homo, virtual, mp2_env, Eigenval, unit_nr, do_ri_sos_laplace_mp2)
      TYPE(rpa_grad_type), INTENT(OUT)                   :: rpa_grad
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_grad_create'

      INTEGER                                            :: handle, ispin, nspins

      CALL timeset(routineN, handle)

      CALL cp_fm_create(rpa_grad%fm_Gamma_PQ, matrix_struct=fm_mat_Q%matrix_struct)
      CALL cp_fm_set_all(rpa_grad%fm_Gamma_PQ, 0.0_dp)

      nspins = SIZE(fm_mat_S)

      ALLOCATE (rpa_grad%fm_Y(nspins))
      DO ispin = 1, nspins
         CALL cp_fm_create(rpa_grad%fm_Y(ispin), fm_mat_S(ispin)%matrix%matrix_struct)

      END DO

      IF (do_ri_sos_laplace_mp2) THEN
         CALL sos_mp2_work_type_create(rpa_grad%sos_mp2_work_occ, rpa_grad%sos_mp2_work_virt, &
                                       unit_nr, Eigenval, homo, virtual, mp2_env%ri_grad%eps_canonical, fm_mat_S)
      ELSE
         CALL rpa_work_type_create(rpa_grad%rpa_work, fm_mat_Q, fm_mat_S, homo, virtual)
      END IF

      CALL timestop(handle)

   END SUBROUTINE rpa_grad_create

! **************************************************************************************************
!> \brief ...
!> \param sos_mp2_work_occ ...
!> \param sos_mp2_work_virt ...
!> \param unit_nr ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param eps_degenerate ...
!> \param fm_mat_S ...
! **************************************************************************************************
   SUBROUTINE sos_mp2_work_type_create(sos_mp2_work_occ, sos_mp2_work_virt, unit_nr, &
                                       Eigenval, homo, virtual, eps_degenerate, fm_mat_S)
      TYPE(sos_mp2_grad_work_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: sos_mp2_work_occ, sos_mp2_work_virt
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), INTENT(IN)                          :: eps_degenerate
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sos_mp2_work_type_create'

      INTEGER                                            :: handle, ispin, nspins

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_mat_S)
      ALLOCATE (sos_mp2_work_occ(nspins), sos_mp2_work_virt(nspins))
      DO ispin = 1, nspins

         CALL create_list_nearly_degen_pairs(Eigenval(1:homo(ispin), ispin), &
                                             eps_degenerate, sos_mp2_work_occ(ispin)%pair_list)
         IF (unit_nr > 0) WRITE (unit_nr, "(T3,A,T75,i6)") &
            "MO_INFO| Number of ij pairs below EPS_CANONICAL:", SIZE(sos_mp2_work_occ(ispin)%pair_list, 2)
         ALLOCATE (sos_mp2_work_occ(ispin)%P(homo(ispin) + SIZE(sos_mp2_work_occ(ispin)%pair_list, 2)))
         sos_mp2_work_occ(ispin)%P = 0.0_dp
         CALL prepare_comm_Pij(sos_mp2_work_occ(ispin), virtual(ispin), fm_mat_S(ispin)%matrix)

         CALL create_list_nearly_degen_pairs(Eigenval(homo(ispin) + 1:, ispin), &
                                             eps_degenerate, sos_mp2_work_virt(ispin)%pair_list)
         IF (unit_nr > 0) WRITE (unit_nr, "(T3,A,T75,i6)") &
            "MO_INFO| Number of ab pairs below EPS_CANONICAL:", SIZE(sos_mp2_work_virt(ispin)%pair_list, 2)
         ALLOCATE (sos_mp2_work_virt(ispin)%P(virtual(ispin) + SIZE(sos_mp2_work_virt(ispin)%pair_list, 2)))
         sos_mp2_work_virt(ispin)%P = 0.0_dp
         CALL prepare_comm_Pab(sos_mp2_work_virt(ispin), virtual(ispin), fm_mat_S(ispin)%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rpa_work ...
!> \param fm_mat_Q ...
!> \param fm_mat_S ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE rpa_work_type_create(rpa_work, fm_mat_Q, fm_mat_S, homo, virtual)
      TYPE(rpa_grad_work_type), INTENT(OUT)              :: rpa_work
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_work_type_create'

      INTEGER :: avirt, col_global, col_local, first_p_pos(2), first_p_pos_col, handle, iocc, &
         ispin, my_a, my_a_end, my_a_size, my_a_start, my_i, my_i_end, my_i_size, my_i_start, &
         my_pcol, ncol_block, ncol_local, nspins, num_pe_col, proc_homo, proc_homo_send, &
         proc_recv, proc_send, proc_virtual, proc_virtual_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data2recv, data2send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices

      CALL timeset(routineN, handle)

      CALL cp_fm_create(rpa_work%fm_mat_Q_copy, matrix_struct=fm_mat_Q%matrix_struct)

      CALL get_blacs_info(fm_mat_S(1)%matrix%matrix_struct%context, &
                          number_of_process_columns=num_pe_col, my_process_column=my_pcol)

      nspins = SIZE(fm_mat_S)

      ALLOCATE (rpa_work%index2send(0:num_pe_col - 1, nspins), &
                rpa_work%index2recv(0:num_pe_col - 1, nspins), &
                rpa_work%gd_homo(nspins), rpa_work%gd_virtual(nspins), &
                data2send(0:num_pe_col - 1), data2recv(0:num_pe_col - 1), &
                rpa_work%P_ij(nspins), rpa_work%P_ab(nspins))

      ! Determine new process grid
      proc_homo = MAX(1, CEILING(SQRT(REAL(num_pe_col, KIND=dp))))
      DO WHILE (MOD(num_pe_col, proc_homo) /= 0)
         proc_homo = proc_homo - 1
      END DO
      proc_virtual = num_pe_col/proc_homo

      rpa_work%grid(1) = proc_virtual
      rpa_work%grid(2) = proc_homo

      rpa_work%mepos(1) = MOD(my_pcol, proc_virtual)
      rpa_work%mepos(2) = my_pcol/proc_virtual

      DO ispin = 1, nspins

         ! Determine distributions of the orbitals
         CALL create_group_dist(rpa_work%gd_homo(ispin), proc_homo, homo(ispin))
         CALL create_group_dist(rpa_work%gd_virtual(ispin), proc_virtual, virtual(ispin))

         CALL cp_fm_struct_get(fm_mat_S(ispin)%matrix%matrix_struct, ncol_local=ncol_local, col_indices=col_indices, &
                               first_p_pos=first_p_pos, ncol_block=ncol_block)
         first_p_pos_col = first_p_pos(2)

         data2send = 0
         ! Count the amount of data2send to each process
         DO col_local = 1, ncol_local
            col_global = col_indices(col_local)

            iocc = (col_global - 1)/virtual(ispin) + 1
            avirt = col_global - (iocc - 1)*virtual(ispin)

            proc_homo_send = group_dist_proc(rpa_work%gd_homo(ispin), iocc)
            proc_virtual_send = group_dist_proc(rpa_work%gd_virtual(ispin), avirt)

            proc_send = proc_homo_send*proc_virtual + proc_virtual_send

            data2send(proc_send) = data2send(proc_send) + 1
         END DO

         DO proc_send = 0, num_pe_col - 1
            ALLOCATE (rpa_work%index2send(proc_send, ispin)%array(data2send(proc_send)))
         END DO

         ! Prepare the indices
         data2send = 0
         DO col_local = 1, ncol_local
            col_global = col_indices(col_local)

            iocc = (col_global - 1)/virtual(ispin) + 1
            avirt = col_global - (iocc - 1)*virtual(ispin)

            proc_homo_send = group_dist_proc(rpa_work%gd_homo(ispin), iocc)
            proc_virtual_send = group_dist_proc(rpa_work%gd_virtual(ispin), avirt)

            proc_send = proc_homo_send*proc_virtual + proc_virtual_send

            data2send(proc_send) = data2send(proc_send) + 1

            rpa_work%index2send(proc_send, ispin)%array(data2send(proc_send)) = col_local
         END DO

         ! Count the amount of data2recv from each process
         CALL get_group_dist(rpa_work%gd_homo(ispin), my_pcol/proc_virtual, my_i_start, my_i_end, my_i_size)
         CALL get_group_dist(rpa_work%gd_virtual(ispin), MOD(my_pcol, proc_virtual), my_a_start, my_a_end, my_a_size)

         data2recv = 0
         DO my_i = my_i_start, my_i_end
         DO my_a = my_a_start, my_a_end
            proc_recv = cp_fm_indxg2p((my_i - 1)*virtual(ispin) + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
            data2recv(proc_recv) = data2recv(proc_recv) + 1
         END DO
         END DO

         DO proc_recv = 0, num_pe_col - 1
            ALLOCATE (rpa_work%index2recv(proc_recv, ispin)%array(2, data2recv(proc_recv)))
         END DO

         data2recv = 0
         DO my_i = my_i_start, my_i_end
         DO my_a = my_a_start, my_a_end
            proc_recv = cp_fm_indxg2p((my_i - 1)*virtual(ispin) + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
            data2recv(proc_recv) = data2recv(proc_recv) + 1

            rpa_work%index2recv(proc_recv, ispin)%array(2, data2recv(proc_recv)) = my_i - my_i_start + 1
            rpa_work%index2recv(proc_recv, ispin)%array(1, data2recv(proc_recv)) = my_a - my_a_start + 1
         END DO
         END DO

         ALLOCATE (rpa_work%P_ij(ispin)%array(my_i_size, homo(ispin)), &
                   rpa_work%P_ab(ispin)%array(my_a_size, virtual(ispin)))
         rpa_work%P_ij(ispin)%array = 0.0_dp
         rpa_work%P_ab(ispin)%array = 0.0_dp

      END DO

      DEALLOCATE (data2send, data2recv)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param Eigenval ...
!> \param eps_degen ...
!> \param pair_list ...
! **************************************************************************************************
   SUBROUTINE create_list_nearly_degen_pairs(Eigenval, eps_degen, pair_list)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: eps_degen
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: pair_list

      INTEGER                                            :: my_i, my_j, num_orbitals, num_pairs, &
                                                            pair_counter

      num_orbitals = SIZE(Eigenval)

! Determine number of nearly degenerate orbital pairs
! Trivial cases: diagonal elements
      num_pairs = 0
      DO my_i = 1, num_orbitals
      DO my_j = my_i + 1, num_orbitals
         IF (ABS(Eigenval(my_i) - Eigenval(my_j)) < eps_degen) num_pairs = num_pairs + 1
      END DO
      END DO
      ALLOCATE (pair_list(2, num_pairs))

! Print the required pairs
      pair_counter = 1
      DO my_i = 1, num_orbitals
      DO my_j = my_i + 1, num_orbitals
      IF (ABS(Eigenval(my_i) - Eigenval(my_j)) < eps_degen) THEN
         pair_list(1, pair_counter) = my_i
         pair_list(2, pair_counter) = my_j
         pair_counter = pair_counter + 1
      END IF
      END DO
      END DO

   END SUBROUTINE create_list_nearly_degen_pairs

! **************************************************************************************************
!> \brief ...
!> \param sos_mp2_work ...
!> \param virtual ...
!> \param fm_mat_S ...
! **************************************************************************************************
   SUBROUTINE prepare_comm_Pij(sos_mp2_work, virtual, fm_mat_S)
      TYPE(sos_mp2_grad_work_type), INTENT(INOUT)        :: sos_mp2_work
      INTEGER, INTENT(IN)                                :: virtual
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prepare_comm_Pij'

      INTEGER :: avirt, col_global, col_local, counter, first_p_pos(2), first_p_pos_col, handle, &
         ij_counter, iocc, my_i, my_j, my_pcol, my_prow, ncol_block, ncol_local, nrow_local, &
         num_ij_pairs, num_pe_col, pcol, pcol_recv, pcol_send, proc_shift, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data2recv, data2send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, ncol_locals
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_exchange

      CALL timeset(routineN, handle)

      tag = 44

      CALL get_blacs_info(fm_mat_S%matrix_struct%context, number_of_process_columns=num_pe_col)
      ALLOCATE (sos_mp2_work%index2send(0:num_pe_col - 1), &
                sos_mp2_work%index2recv(0:num_pe_col - 1))

      ALLOCATE (data2send(0:num_pe_col - 1))
      ALLOCATE (data2recv(0:num_pe_col - 1))

      CALL cp_fm_struct_get(fm_mat_S%matrix_struct, para_env=para_env, ncol_locals=ncol_locals, &
                            ncol_local=ncol_local, col_indices=col_indices, &
                            context=context, ncol_block=ncol_block, first_p_pos=first_p_pos, nrow_local=nrow_local)
      CALL get_blacs_info(context, my_process_row=my_prow, my_process_column=my_pcol, &
                          blacs2mpi=blacs2mpi)
      first_p_pos_col = first_p_pos(2)

      num_ij_pairs = SIZE(sos_mp2_work%pair_list, 2)

      IF (num_ij_pairs > 0) THEN

         NULLIFY (para_env_exchange)
         CALL cp_para_env_split(para_env_exchange, para_env, my_prow)

         data2send = 0
         data2recv = 0

         DO proc_shift = 0, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            counter = 0
            DO col_local = 1, ncol_local
               col_global = col_indices(col_local)

               iocc = MAX(1, col_global - 1)/virtual + 1
               avirt = col_global - (iocc - 1)*virtual

               DO ij_counter = 1, num_ij_pairs

                  my_i = sos_mp2_work%pair_list(1, ij_counter)
                  my_j = sos_mp2_work%pair_list(2, ij_counter)

                  IF (iocc /= my_j) CYCLE
                  pcol = cp_fm_indxg2p((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)
                  IF (pcol /= pcol_send) CYCLE

                  counter = counter + 1

                  EXIT

               END DO
            END DO
            data2send(pcol_send) = counter
         END DO

         CALL mp_alltoall(data2send, data2recv, 1, para_env_exchange%group)

         DO proc_shift = 0, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            ! Collect indices and exchange
            ALLOCATE (sos_mp2_work%index2send(pcol_send)%array(data2send(pcol_send)))

            counter = 0
            DO col_local = 1, ncol_local
               col_global = col_indices(col_local)

               iocc = MAX(1, col_global - 1)/virtual + 1
               avirt = col_global - (iocc - 1)*virtual

               DO ij_counter = 1, num_ij_pairs

                  my_i = sos_mp2_work%pair_list(1, ij_counter)
                  my_j = sos_mp2_work%pair_list(2, ij_counter)

                  IF (iocc /= my_j) CYCLE
                  pcol = cp_fm_indxg2p((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)
                  IF (pcol /= pcol_send) CYCLE

                  counter = counter + 1

                  sos_mp2_work%index2send(pcol_send)%array(counter) = col_global

                  EXIT

               END DO
            END DO

            ALLOCATE (sos_mp2_work%index2recv(pcol_recv)%array(data2recv(pcol_recv)))
            !
            CALL mp_sendrecv(sos_mp2_work%index2send(pcol_send)%array, blacs2mpi(my_prow, pcol_send), &
                             sos_mp2_work%index2recv(pcol_recv)%array, blacs2mpi(my_prow, pcol_recv), &
                             para_env%group, tag)
         END DO

         CALL cp_para_env_release(para_env_exchange)
      END IF

      DEALLOCATE (data2send, data2recv)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param sos_mp2_work ...
!> \param virtual ...
!> \param fm_mat_S ...
! **************************************************************************************************
   SUBROUTINE prepare_comm_Pab(sos_mp2_work, virtual, fm_mat_S)
      TYPE(sos_mp2_grad_work_type), INTENT(INOUT)        :: sos_mp2_work
      INTEGER, INTENT(IN)                                :: virtual
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prepare_comm_Pab'

      INTEGER :: ab_counter, avirt, col_global, col_local, counter, first_p_pos(2), &
         first_p_pos_col, handle, iocc, my_a, my_b, my_pcol, my_prow, ncol_block, ncol_local, &
         nrow_local, num_ab_pairs, num_pe_col, pcol, pcol_recv, pcol_send, proc_shift, tag
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data2recv, data2send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, ncol_locals
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_exchange

      CALL timeset(routineN, handle)

      tag = 44

      CALL get_blacs_info(fm_mat_S%matrix_struct%context, number_of_process_columns=num_pe_col)
      ALLOCATE (sos_mp2_work%index2send(0:num_pe_col - 1), &
                sos_mp2_work%index2recv(0:num_pe_col - 1))

      num_ab_pairs = SIZE(sos_mp2_work%pair_list, 2)
      IF (num_ab_pairs > 0) THEN

         CALL cp_fm_struct_get(fm_mat_S%matrix_struct, para_env=para_env, ncol_locals=ncol_locals, &
                               ncol_local=ncol_local, col_indices=col_indices, &
                               context=context, ncol_block=ncol_block, first_p_pos=first_p_pos, nrow_local=nrow_local)
         CALL get_blacs_info(context, my_process_row=my_prow, my_process_column=my_pcol, &
                             blacs2mpi=blacs2mpi)
         first_p_pos_col = first_p_pos(2)

         NULLIFY (para_env_exchange)
         CALL cp_para_env_split(para_env_exchange, para_env, my_prow)

         ALLOCATE (data2send(0:num_pe_col - 1))
         ALLOCATE (data2recv(0:num_pe_col - 1))

         data2send = 0
         data2recv = 0
         DO proc_shift = 0, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            counter = 0
            DO col_local = 1, ncol_local
               col_global = col_indices(col_local)

               iocc = MAX(1, col_global - 1)/virtual + 1
               avirt = col_global - (iocc - 1)*virtual

               DO ab_counter = 1, num_ab_pairs

                  my_a = sos_mp2_work%pair_list(1, ab_counter)
                  my_b = sos_mp2_work%pair_list(2, ab_counter)

                  IF (avirt /= my_b) CYCLE
                  pcol = cp_fm_indxg2p((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
                  IF (pcol /= pcol_send) CYCLE

                  counter = counter + 1

                  EXIT

               END DO
            END DO
            data2send(pcol_send) = counter
         END DO

         CALL mp_alltoall(data2send, data2recv, 1, para_env_exchange%group)

         DO proc_shift = 0, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            ! Collect indices and exchange
            ALLOCATE (sos_mp2_work%index2send(pcol_send)%array(data2send(pcol_send)))

            counter = 0
            DO col_local = 1, ncol_local
               col_global = col_indices(col_local)

               iocc = MAX(1, col_global - 1)/virtual + 1
               avirt = col_global - (iocc - 1)*virtual

               DO ab_counter = 1, num_ab_pairs

                  my_a = sos_mp2_work%pair_list(1, ab_counter)
                  my_b = sos_mp2_work%pair_list(2, ab_counter)

                  IF (avirt /= my_b) CYCLE
                  pcol = cp_fm_indxg2p((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
                  IF (pcol /= pcol_send) CYCLE

                  counter = counter + 1

                  sos_mp2_work%index2send(pcol_send)%array(counter) = col_global

                  EXIT

               END DO
            END DO

            ALLOCATE (sos_mp2_work%index2recv(pcol_recv)%array(data2recv(pcol_recv)))
            !
            CALL mp_sendrecv(sos_mp2_work%index2send(pcol_send)%array, blacs2mpi(my_prow, pcol_send), &
                             sos_mp2_work%index2recv(pcol_recv)%array, blacs2mpi(my_prow, pcol_recv), &
                             para_env%group, tag)
         END DO

         CALL cp_para_env_release(para_env_exchange)
         DEALLOCATE (data2send, data2recv)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_Q ...
!> \param rpa_grad ...
! **************************************************************************************************
   SUBROUTINE rpa_grad_copy_Q(fm_mat_Q, rpa_grad)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      TYPE(rpa_grad_type), INTENT(INOUT)                 :: rpa_grad

      CALL cp_fm_to_fm(fm_mat_Q, rpa_grad%rpa_work%fm_mat_Q_copy)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param rpa_grad ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param dgemm_counter ...
!> \param fm_mat_S ...
!> \param omega ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param weight ...
! **************************************************************************************************
   SUBROUTINE rpa_grad_matrix_operations(mp2_env, rpa_grad, do_ri_sos_laplace_mp2, fm_mat_Q, fm_mat_Q_gemm, &
                                         dgemm_counter, fm_mat_S, omega, homo, virtual, Eigenval, weight)
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(rpa_grad_type), INTENT(INOUT)                 :: rpa_grad
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_Q, fm_mat_Q_gemm
      TYPE(dgemm_counter_type), INTENT(INOUT)            :: dgemm_counter
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: fm_mat_S
      REAL(KIND=dp), INTENT(IN)                          :: omega
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: weight

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_grad_matrix_operations'

      INTEGER                                            :: col_global, col_local, dimen_ia, &
                                                            dimen_RI, handle, handle2, ispin, &
                                                            jspin, ncol_local, nrow_local, nspins, &
                                                            row_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         TARGET                                          :: mat_S_3D, mat_work_iaP_3D
      TYPE(cp_fm_type)                                   :: fm_work_iaP, fm_work_PQ

      CALL timeset(routineN, handle)

      nspins = SIZE(fm_mat_Q)

      CALL cp_fm_get_info(fm_mat_Q(1)%matrix, nrow_global=dimen_RI, nrow_local=nrow_local, ncol_local=ncol_local, &
                          col_indices=col_indices, row_indices=row_indices)

      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         CALL cp_fm_create(fm_work_PQ, fm_mat_Q(1)%matrix%matrix_struct)

         ! calculate [1+Q(iw')]^-1
         CALL cp_fm_cholesky_invert(fm_mat_Q(1)%matrix)
         ! symmetrize the result, fm_work_PQ is only a work matrix
         CALL cp_fm_upper_to_full(fm_mat_Q(1)%matrix, fm_work_PQ)

         CALL cp_fm_release(fm_work_PQ)

         DO col_local = 1, ncol_local
            col_global = col_indices(col_local)
            DO row_local = 1, nrow_local
            IF (col_global == row_indices(row_local)) THEN
               fm_mat_Q(1)%matrix%local_data(row_local, col_local) = fm_mat_Q(1)%matrix%local_data(row_local, col_local) - 1.0_dp
               EXIT
            END IF
            END DO
         END DO

         CALL timeset(routineN//"_PQ", handle2)
         CALL dgemm_counter_start(dgemm_counter)
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=weight, &
                            matrix_a=rpa_grad%rpa_work%fm_mat_Q_copy, matrix_b=fm_mat_Q(1)%matrix, beta=1.0_dp, &
                            matrix_c=rpa_grad%fm_Gamma_PQ)
         CALL dgemm_counter_stop(dgemm_counter, dimen_RI, dimen_RI, dimen_RI)
         CALL timestop(handle2)

         CALL cp_fm_to_fm_submat_general(fm_mat_Q(1)%matrix, fm_mat_Q_gemm(1)%matrix, dimen_RI, dimen_RI, 1, 1, 1, 1, &
                                         fm_mat_Q_gemm(1)%matrix%matrix_struct%context%group)
      END IF

      DO ispin = 1, nspins
         IF (do_ri_sos_laplace_mp2) THEN
            ! The spin of the other Q matrix is always the other spin
            jspin = nspins - ispin + 1
! or the first matrix in the case of RPA
         ELSE
            jspin = 1
         END IF

         IF (do_ri_sos_laplace_mp2) THEN
            CALL timeset(routineN//"_PQ", handle2)
            CALL dgemm_counter_start(dgemm_counter)
            CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=weight, &
                               matrix_a=fm_mat_Q(ispin)%matrix, matrix_b=fm_mat_Q(jspin)%matrix, beta=1.0_dp, &
                               matrix_c=rpa_grad%fm_Gamma_PQ)
            CALL dgemm_counter_stop(dgemm_counter, dimen_RI, dimen_RI, dimen_RI)
            CALL timestop(handle2)

            CALL cp_fm_to_fm_submat_general(fm_mat_Q(jspin)%matrix, fm_mat_Q_gemm(jspin)%matrix, dimen_RI, dimen_RI, 1, 1, 1, 1, &
                                            fm_mat_Q_gemm(jspin)%matrix%matrix_struct%context%group)
         ELSE
            CALL calc_fm_mat_S_rpa(fm_mat_S(ispin)%matrix, .TRUE., virtual(ispin), Eigenval(:, ispin), &
                                   homo(ispin), omega, 0.0_dp)
         END IF

         CALL timeset(routineN//"_contr_S", handle2)
         CALL cp_fm_create(fm_work_iaP, rpa_grad%fm_Y(ispin)%matrix_struct)

         CALL cp_fm_get_info(fm_mat_S(ispin)%matrix, ncol_global=dimen_ia)

         CALL dgemm_counter_start(dgemm_counter)
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_ia, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_mat_Q_gemm(jspin)%matrix, matrix_b=fm_mat_S(ispin)%matrix, beta=0.0_dp, &
                            matrix_c=fm_work_iaP)
         CALL dgemm_counter_stop(dgemm_counter, dimen_ia, dimen_RI, dimen_RI)
         CALL timestop(handle2)

         IF (do_ri_sos_laplace_mp2) THEN
            CALL calc_P_sos_mp2(homo(ispin), fm_mat_S(ispin)%matrix, fm_work_iaP, &
                                rpa_grad%sos_mp2_work_occ(ispin), rpa_grad%sos_mp2_work_virt(ispin), &
                                omega, weight, virtual(ispin), Eigenval(:, ispin))

            CALL calc_fm_mat_S_laplace(fm_work_iaP, homo(ispin), virtual(ispin), Eigenval(:, ispin), omega)

            CALL cp_fm_scale_and_add(1.0_dp, rpa_grad%fm_Y(ispin), -weight, fm_work_iaP)

            CALL cp_fm_release(fm_work_iaP)
         ELSE
            ! To save memory, we add it now
            CALL cp_fm_scale_and_add(1.0_dp, rpa_grad%fm_Y(ispin), -weight, fm_work_iaP)

            ! Redistribute both matrices and deallocate fm_work_iaP
            CALL redistribute_fm_mat_S(rpa_grad%rpa_work%index2send(:, ispin), rpa_grad%rpa_work%index2recv(:, ispin), &
                                       fm_work_iaP, mat_work_iaP_3D, &
                                       rpa_grad%rpa_work%gd_homo(ispin), rpa_grad%rpa_work%gd_virtual(ispin), &
                                       rpa_grad%rpa_work%mepos)
            CALL cp_fm_release(fm_work_iaP)

            CALL redistribute_fm_mat_S(rpa_grad%rpa_work%index2send(:, ispin), rpa_grad%rpa_work%index2recv(:, ispin), &
                                       fm_mat_S(ispin)%matrix, mat_S_3D, &
                                       rpa_grad%rpa_work%gd_homo(ispin), rpa_grad%rpa_work%gd_virtual(ispin), &
                                       rpa_grad%rpa_work%mepos)

            ! Now collect the density matrix
            CALL calc_P_rpa(mat_S_3D, mat_work_iaP_3D, rpa_grad%rpa_work%gd_homo(ispin), rpa_grad%rpa_work%gd_virtual(ispin), &
                            rpa_grad%rpa_work%grid, rpa_grad%rpa_work%mepos, &
                            fm_mat_S(ispin)%matrix%matrix_struct, &
                            rpa_grad%rpa_work%P_ij(ispin)%array, rpa_grad%rpa_work%P_ab(ispin)%array, &
                            weight, omega, Eigenval(:, ispin), homo(ispin), mp2_env)

            DEALLOCATE (mat_work_iaP_3D, mat_S_3D)

            CALL remove_scaling_factor_rpa(fm_mat_S(ispin)%matrix, virtual(ispin), Eigenval(:, ispin), homo(ispin), omega)

         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param homo ...
!> \param fm_mat_S ...
!> \param fm_work_iaP ...
!> \param sos_mp2_work_occ ...
!> \param sos_mp2_work_virt ...
!> \param omega ...
!> \param weight ...
!> \param virtual ...
!> \param Eigenval ...
! **************************************************************************************************
   SUBROUTINE calc_P_sos_mp2(homo, fm_mat_S, fm_work_iaP, sos_mp2_work_occ, sos_mp2_work_virt, omega, weight, virtual, Eigenval)
      INTEGER, INTENT(IN)                                :: homo
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_work_iaP
      TYPE(sos_mp2_grad_work_type), INTENT(INOUT)        :: sos_mp2_work_occ, sos_mp2_work_virt
      REAL(KIND=dp), INTENT(IN)                          :: omega, weight
      INTEGER, INTENT(IN)                                :: virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_P_sos_mp2'

      INTEGER                                            :: avirt, col_global, col_local, handle, &
                                                            handle2, iocc, my_a, my_i, ncol_local, &
                                                            nrow_local, num_ab_pairs, num_ij_pairs
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp)                                      :: ddot, trace

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_mat_S, col_indices=col_indices, ncol_local=ncol_local, nrow_local=nrow_local)

      CALL timeset(routineN//"_Pij_diag", handle2)
      DO my_i = 1, homo
         ! Collect the contributions of the matrix elements

         trace = 0.0_dp

         DO col_local = 1, ncol_local
            col_global = col_indices(col_local)

            iocc = MAX(1, col_global - 1)/virtual + 1
            avirt = col_global - (iocc - 1)*virtual

            IF (iocc == my_i) trace = trace + &
                                     ddot(nrow_local, fm_mat_S%local_data(:, col_local), 1, fm_work_iaP%local_data(:, col_local), 1)
         END DO

         sos_mp2_work_occ%P(my_i) = sos_mp2_work_occ%P(my_i) - trace*omega*weight

      END DO
      CALL timestop(handle2)

      CALL timeset(routineN//"_Pab_diag", handle2)
      DO my_a = 1, virtual
         ! Collect the contributions of the matrix elements

         trace = 0.0_dp

         DO col_local = 1, ncol_local
            col_global = col_indices(col_local)

            iocc = MAX(1, col_global - 1)/virtual + 1
            avirt = col_global - (iocc - 1)*virtual

            IF (avirt == my_a) trace = trace + &
                                     ddot(nrow_local, fm_mat_S%local_data(:, col_local), 1, fm_work_iaP%local_data(:, col_local), 1)
         END DO

         sos_mp2_work_virt%P(my_a) = sos_mp2_work_virt%P(my_a) + trace*omega*weight

      END DO
      CALL timestop(handle2)

      ! Loop over list and carry out operations
      num_ij_pairs = SIZE(sos_mp2_work_occ%pair_list, 2)
      num_ab_pairs = SIZE(sos_mp2_work_virt%pair_list, 2)
      IF (num_ij_pairs > 0) THEN
         CALL calc_Pij_degen(fm_work_iaP, fm_mat_S, sos_mp2_work_occ%pair_list, &
                             virtual, sos_mp2_work_occ%P(homo + 1:), Eigenval(:homo), omega, weight, &
                             sos_mp2_work_occ%index2send, sos_mp2_work_occ%index2recv)
      END IF
      IF (num_ab_pairs > 0) THEN
         CALL calc_Pab_degen(fm_work_iaP, fm_mat_S, sos_mp2_work_virt%pair_list, &
                             virtual, sos_mp2_work_virt%P(virtual + 1:), Eigenval(homo + 1:), omega, weight, &
                             sos_mp2_work_virt%index2send, sos_mp2_work_virt%index2recv)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_S_1D ...
!> \param mat_work_iaP_3D ...
!> \param gd_homo ...
!> \param gd_virtual ...
!> \param grid ...
!> \param mepos ...
!> \param fm_struct_S ...
!> \param P_ij ...
!> \param P_ab ...
!> \param weight ...
!> \param omega ...
!> \param Eigenval ...
!> \param homo ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE calc_P_rpa(mat_S_1D, mat_work_iaP_3D, gd_homo, gd_virtual, grid, mepos, &
                         fm_struct_S, P_ij, P_ab, weight, omega, Eigenval, homo, mp2_env)
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT), TARGET :: mat_S_1D
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: mat_work_iaP_3D
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_homo, gd_virtual
      INTEGER, DIMENSION(2), INTENT(IN)                  :: grid, mepos
      TYPE(cp_fm_struct_type), INTENT(IN), POINTER       :: fm_struct_S
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: P_ij, P_ab
      REAL(KIND=dp), INTENT(IN)                          :: weight, omega
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: homo
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_P_rpa'

      INTEGER :: handle, handle2, my_a, my_a_size, my_a_start, my_i, my_i_size, my_i_start, &
         my_P_size, my_prow, proc_a_recv, proc_a_send, proc_i_recv, proc_i_send, proc_recv, &
         proc_send, proc_shift, recv_a_end, recv_a_size, recv_a_start, recv_i_end, recv_i_size, &
         recv_i_start, tag
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: buffer_1D
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: buffer_2D, mat_S_2D
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: buffer_3D, mat_S_3D
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ! We allocate it at every step to reduce potential memory conflicts with COSMA
      IF (.NOT. c_associated(mp2_env%local_gemm_ctx)) THEN
         CALL local_gemm_create(mp2_env%local_gemm_ctx, LOCAL_GEMM_PU_GPU)
         CALL local_gemm_set_op_threshold_gpu(mp2_env%local_gemm_ctx, 128*128*128*2)
      END IF

      tag = 47

      my_P_size = SIZE(mat_work_iaP_3D, 1)

      CALL cp_fm_struct_get(fm_struct_S, para_env=para_env)
      CALL get_blacs_info(fm_struct_S%context, my_process_row=my_prow, blacs2mpi=blacs2mpi, para_env=para_env)

      CALL get_group_dist(gd_virtual, mepos(1), my_a_start, sizes=my_a_size)
      CALL get_group_dist(gd_homo, mepos(2), my_i_start, sizes=my_i_size)

      ! We have to remap the indices because mp_sendrecv requires a 3D array (because of mat_work_iaP_3D)
      ! and dgemm requires 2D arrays
      ! Fortran 2008 does allow pointer remapping independently of the ranks but GCC 7 does not properly support it
      mat_S_2D(1:my_P_size*my_a_size, 1:my_i_size) => mat_S_1D(1:INT(my_P_size, int_8)*my_a_size*my_i_size)
      mat_S_3D(1:my_P_size, 1:my_a_size, 1:my_i_size) => mat_S_1D(1:INT(my_P_size, int_8)*my_a_size*my_i_size)

      ALLOCATE (buffer_1D(MAX(INT(maxsize(gd_homo), KIND=int_8)*my_a_size, &
                              INT(maxsize(gd_virtual), KIND=int_8)*my_i_size)*my_P_size))

      DO proc_shift = 0, grid(1) - 1
         proc_a_send = MODULO(mepos(1) + proc_shift, grid(1))
         proc_a_recv = MODULO(mepos(1) - proc_shift, grid(1))

         proc_send = mepos(2)*grid(1) + proc_a_send
         proc_recv = mepos(2)*grid(1) + proc_a_recv

         CALL get_group_dist(gd_virtual, proc_a_recv, recv_a_start, recv_a_end, recv_a_size)

         buffer_3D(1:my_P_size, 1:recv_a_size, 1:my_i_size) => buffer_1D(1:INT(my_P_size, KIND=int_8)*recv_a_size*my_i_size)

         CALL timeset(routineN//"_comm_a", handle2)
         CALL mp_sendrecv(mat_work_iaP_3D, blacs2mpi(my_prow, proc_send), &
                          buffer_3D, blacs2mpi(my_prow, proc_recv), &
                          para_env%group, tag)
         CALL timestop(handle2)

         CALL timeset(routineN//"_dgemm_a", handle2)
         ! This loop imitates the actual tensor contraction
         DO my_i = 1, my_i_size
            CALL local_gemm("T", "N", my_a_size, recv_a_size, my_P_size, &
                            -weight, mat_S_3D(:, :, my_i), my_P_size, &
                            buffer_3D(:, :, my_i), my_P_size, &
                            1.0_dp, P_ab(:, recv_a_start:recv_a_end), my_a_size, mp2_env%local_gemm_ctx)
         END DO
         CALL timestop(handle2)

         CALL timeset(routineN//"_scale_a", handle2)
         DO my_i = 1, my_i_size
            ! Scale both matrices
            DO my_a = 1, my_a_size
               mat_S_3D(:, my_a, my_i) = mat_S_3D(:, my_a, my_i)* &
                                         (omega/(Eigenval(homo + my_a + my_a_start - 1) - Eigenval(my_i + my_i_start - 1)))
            END DO
            DO my_a = 1, recv_a_size
               buffer_3D(:, my_a, my_i) = buffer_3D(:, my_a, my_i)* &
                                          (omega/(Eigenval(homo + my_a + recv_a_start - 1) - Eigenval(my_i + my_i_start - 1)))
            END DO
         END DO
         CALL timestop(handle2)

         CALL timeset(routineN//"_dgemm_a", handle2)
         DO my_i = 1, my_i_size
            ! Repeat contraction
            CALL local_gemm("T", "N", my_a_size, recv_a_size, my_P_size, &
                            weight, mat_S_3D(:, :, my_i), my_P_size, &
                            buffer_3D(:, :, my_i), my_P_size, &
                            1.0_dp, P_ab(:, recv_a_start:recv_a_end), my_a_size, mp2_env%local_gemm_ctx)
         END DO
         CALL timestop(handle2)

         CALL timeset(routineN//"_scale_a", handle2)
         ! Remove scaling of mat_S for later
         DO my_i = 1, my_i_size
            DO my_a = 1, my_a_size
               mat_S_3D(:, my_a, my_i) = mat_S_3D(:, my_a, my_i)* &
                                         ((Eigenval(homo + my_a + my_a_start - 1) - Eigenval(my_i + my_i_start - 1))/omega)
            END DO
            DO my_a = 1, recv_a_size
               buffer_3D(:, my_a, my_i) = buffer_3D(:, my_a, my_i)* &
                                          ((Eigenval(homo + my_a + recv_a_start - 1) - Eigenval(my_i + my_i_start - 1))/omega)
            END DO
         END DO
         CALL timestop(handle2)
      END DO

      DO proc_shift = 0, grid(2) - 1
         proc_i_send = MODULO(mepos(2) + proc_shift, grid(2))
         proc_i_recv = MODULO(mepos(2) - proc_shift, grid(2))

         proc_send = proc_i_send*grid(1) + mepos(1)
         proc_recv = proc_i_recv*grid(1) + mepos(1)

         CALL get_group_dist(gd_homo, proc_i_recv, recv_i_start, recv_i_end, recv_i_size)

         ! Both pointers point to the same memory but we need the first for communication, the second one for dgemm
         buffer_3D(1:my_P_size, 1:my_a_size, 1:recv_i_size) => buffer_1D(1:INT(my_P_size, KIND=int_8)*my_a_size*recv_i_size)
         buffer_2D(1:my_P_size*my_a_size, 1:recv_i_size) => buffer_1D(1:INT(my_P_size, KIND=int_8)*my_a_size*recv_i_size)

         CALL timeset(routineN//"_comm_i", handle2)
         CALL mp_sendrecv(mat_work_iaP_3D, blacs2mpi(my_prow, proc_send), &
                          buffer_3D, blacs2mpi(my_prow, proc_recv), &
                          para_env%group, tag)
         CALL timestop(handle2)

         CALL timeset(routineN//"_dgemm_i", handle2)
         ! This loop imitates the actual tensor contraction
         CALL local_gemm("T", "N", my_i_size, recv_i_size, my_a_size*my_P_size, &
                         weight, mat_S_2D, my_P_size*my_a_size, &
                         buffer_2D, my_P_size*my_a_size, &
                         1.0_dp, P_ij(:, recv_i_start:recv_i_end), my_i_size, mp2_env%local_gemm_ctx)
         CALL timestop(handle2)

         CALL timeset(routineN//"_scale_i", handle2)
         ! Scale both matrices
         DO my_i = 1, my_i_size
         DO my_a = 1, my_a_size
            mat_S_3D(:, my_a, my_i) = mat_S_3D(:, my_a, my_i)* &
                                      (omega/(Eigenval(homo + my_a + my_a_start - 1) - Eigenval(my_i + my_i_start - 1)))
         END DO
         END DO
         DO my_i = 1, recv_i_size
         DO my_a = 1, my_a_size
            buffer_3D(:, my_a, my_i) = buffer_3D(:, my_a, my_i)* &
                                       (omega/(Eigenval(homo + my_a + my_a_start - 1) - Eigenval(my_i + recv_i_start - 1)))
         END DO
         END DO
         CALL timestop(handle2)

         CALL timeset(routineN//"_dgemm_i", handle2)
         ! Repeat contraction
         CALL local_gemm("T", "N", my_i_size, recv_i_size, my_a_size*my_P_size, &
                         -weight, mat_S_2D, my_P_size*my_a_size, &
                         buffer_2D, my_P_size*my_a_size, &
                         1.0_dp, P_ij(:, recv_i_start:recv_i_end), my_i_size, mp2_env%local_gemm_ctx)
         CALL timestop(handle2)

         CALL timeset(routineN//"_scale_i", handle2)
         ! Remove scaling factor of mat_S for later
         DO my_i = 1, my_i_size
         DO my_a = 1, my_a_size
            mat_S_3D(:, my_a, my_i) = mat_S_3D(:, my_a, my_i)* &
                                      ((Eigenval(homo + my_a + my_a_start - 1) - Eigenval(my_i + recv_i_start - 1))/omega)
         END DO
         END DO
         CALL timestop(handle2)
      END DO

      ! release memory allocated by local_gemm when run on GPU. local_gemm_ctx is null on cpu only runs
      CALL local_gemm_destroy(mp2_env%local_gemm_ctx)
      mp2_env%local_gemm_ctx = C_NULL_PTR

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION sinh_over_x(x) RESULT(res)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: res

      ! Calculate sinh(x)/x
      ! Split the intervall to prevent numerical instabilities
      IF (ABS(x) > 3.0e-4_dp) THEN
         res = SINH(x)/x
      ELSE
         res = 1.0_dp + x**2/6.0_dp
      END IF

   END FUNCTION sinh_over_x

! **************************************************************************************************
!> \brief ...
!> \param fm_work_iaP ...
!> \param fm_mat_S ...
!> \param pair_list ...
!> \param virtual ...
!> \param P_ij ...
!> \param Eigenval ...
!> \param omega ...
!> \param weight ...
!> \param index2send ...
!> \param index2recv ...
! **************************************************************************************************
   SUBROUTINE calc_Pij_degen(fm_work_iaP, fm_mat_S, pair_list, virtual, P_ij, Eigenval, omega, weight, index2send, index2recv)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_work_iaP, fm_mat_S
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: pair_list
      INTEGER, INTENT(IN)                                :: virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: P_ij
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: omega, weight
      TYPE(one_dim_int_array), DIMENSION(0:), INTENT(IN) :: index2send, index2recv

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_Pij_degen'

      INTEGER :: avirt, col_global, col_local, counter, first_p_pos(2), first_p_pos_col, handle, &
         handle2, ij_counter, iocc, my_col_local, my_i, my_j, my_pcol, my_prow, ncol_block, &
         ncol_local, nrow_local, num_ij_pairs, num_pe_col, pcol, pcol_recv, pcol_send, proc_shift, &
         recv_size, send_size, size_recv_buffer, size_send_buffer, tag
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, ncol_locals
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      REAL(KIND=dp)                                      :: ddot, trace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_recv, buffer_send
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL cp_fm_struct_get(fm_work_iaP%matrix_struct, para_env=para_env, ncol_locals=ncol_locals, &
                            ncol_local=ncol_local, col_indices=col_indices, &
                            context=context, ncol_block=ncol_block, first_p_pos=first_p_pos, nrow_local=nrow_local)
      CALL get_blacs_info(context, my_process_row=my_prow, my_process_column=my_pcol, &
                          number_of_process_columns=num_pe_col, blacs2mpi=blacs2mpi)
      first_p_pos_col = first_p_pos(2)

      num_ij_pairs = SIZE(pair_list, 2)

      tag = 42

      DO col_local = 1, ncol_local
         col_global = col_indices(col_local)

         iocc = MAX(1, col_global - 1)/virtual + 1
         avirt = col_global - (iocc - 1)*virtual

         DO ij_counter = 1, num_ij_pairs

            my_i = pair_list(1, ij_counter)
            my_j = pair_list(2, ij_counter)

            IF (iocc /= my_j) CYCLE
            pcol = cp_fm_indxg2p((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)
            IF (pcol /= pcol_send) CYCLE

            my_col_local = cp_fm_indxg2l((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)

            trace = ddot(nrow_local, fm_mat_S%local_data(:, my_col_local), 1, fm_work_iaP%local_data(:, col_local), 1)

            P_ij(ij_counter) = P_ij(ij_counter) &
                               - trace*sinh_over_x(0.5_dp*(Eigenval(my_i) - Eigenval(my_j))*omega)*omega*weight

            EXIT

         END DO
      END DO

      IF (num_pe_col > 1) THEN
         size_send_buffer = 0
         size_recv_buffer = 0
         DO proc_shift = 1, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            IF (ALLOCATED(index2send(pcol_send)%array)) &
               size_send_buffer = MAX(size_send_buffer, SIZE(index2send(pcol_send)%array))

            IF (ALLOCATED(index2recv(pcol_recv)%array)) &
               size_recv_buffer = MAX(size_recv_buffer, SIZE(index2recv(pcol_recv)%array))
         END DO

         ALLOCATE (buffer_send(nrow_local, size_send_buffer), buffer_recv(nrow_local, size_recv_buffer))

         DO proc_shift = 1, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            ! Collect data and exchange
            send_size = 0
            IF (ALLOCATED(index2send(pcol_send)%array)) send_size = SIZE(index2send(pcol_send)%array)

            IF (send_size > 0) THEN
               counter = 0
               DO col_local = 1, ncol_local
                  col_global = col_indices(col_local)

                  iocc = MAX(1, col_global - 1)/virtual + 1
                  avirt = col_global - (iocc - 1)*virtual

                  DO ij_counter = 1, num_ij_pairs

                     my_i = pair_list(1, ij_counter)
                     my_j = pair_list(2, ij_counter)

                     IF (iocc /= my_j) CYCLE
                     pcol = cp_fm_indxg2p((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)
                     IF (pcol /= pcol_send) CYCLE

                     counter = counter + 1

                     buffer_send(:, counter) = fm_work_iaP%local_data(:, col_local)

                     EXIT

                  END DO
               END DO
            END IF

            recv_size = 0
            IF (ALLOCATED(index2recv(pcol_recv)%array)) recv_size = SIZE(index2recv(pcol_recv)%array)
            IF (recv_size > 0) THEN
               CALL timeset(routineN//"_send", handle2)
               IF (send_size > 0) THEN
                  CALL mp_sendrecv(buffer_send(:, :send_size), blacs2mpi(my_prow, pcol_send), &
                                   buffer_recv(:, :recv_size), blacs2mpi(my_prow, pcol_recv), &
                                   para_env%group, tag)
               ELSE
                  CALL mp_recv(buffer_recv(:, :recv_size), blacs2mpi(my_prow, pcol_recv), tag, para_env%group)
               END IF
               CALL timestop(handle2)

               DO ij_counter = 1, num_ij_pairs
                  ! Collect the contributions of the matrix elements

                  my_i = pair_list(1, ij_counter)
                  my_j = pair_list(2, ij_counter)

                  trace = 0.0_dp

                  DO col_local = 1, recv_size
                     col_global = index2recv(pcol_recv)%array(col_local)

                     iocc = MAX(1, col_global - 1)/virtual + 1
                     IF (iocc /= my_j) CYCLE
                     avirt = col_global - (iocc - 1)*virtual
                     pcol = cp_fm_indxg2p((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)
                     IF (pcol /= my_pcol) CYCLE

                     my_col_local = cp_fm_indxg2l((my_i - 1)*virtual + avirt, ncol_block, 0, first_p_pos_col, num_pe_col)

                     trace = trace + ddot(nrow_local, fm_mat_S%local_data(:, my_col_local), 1, buffer_recv(:, col_local), 1)
                  END DO

                  P_ij(ij_counter) = P_ij(ij_counter) &
                                     - trace*sinh_over_x(0.5_dp*(Eigenval(my_i) - Eigenval(my_j))*omega)*omega*weight
               END DO
            ELSE IF (send_size > 0) THEN
               CALL timeset(routineN//"_send", handle2)
               CALL mp_send(buffer_send(:, :send_size), blacs2mpi(my_prow, pcol_send), tag, para_env%group)
               CALL timestop(handle2)
            END IF
         END DO
         IF (ALLOCATED(buffer_send)) DEALLOCATE (buffer_send)
         IF (ALLOCATED(buffer_recv)) DEALLOCATE (buffer_recv)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_work_iaP ...
!> \param fm_mat_S ...
!> \param pair_list ...
!> \param virtual ...
!> \param P_ab ...
!> \param Eigenval ...
!> \param omega ...
!> \param weight ...
!> \param index2send ...
!> \param index2recv ...
! **************************************************************************************************
   SUBROUTINE calc_Pab_degen(fm_work_iaP, fm_mat_S, pair_list, virtual, P_ab, Eigenval, omega, weight, index2send, index2recv)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_work_iaP, fm_mat_S
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: pair_list
      INTEGER, INTENT(IN)                                :: virtual
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: P_ab
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: omega, weight
      TYPE(one_dim_int_array), DIMENSION(0:), INTENT(IN) :: index2send, index2recv

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_Pab_degen'

      INTEGER :: ab_counter, avirt, col_global, col_local, counter, first_p_pos(2), &
         first_p_pos_col, handle, handle2, iocc, my_a, my_b, my_col_local, my_pcol, my_prow, &
         ncol_block, ncol_local, nrow_local, num_ab_pairs, num_pe_col, pcol, pcol_recv, pcol_send, &
         proc_shift, recv_size, send_size, size_recv_buffer, size_send_buffer, tag
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, ncol_locals
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      REAL(KIND=dp)                                      :: ddot, trace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_recv, buffer_send
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL cp_fm_struct_get(fm_work_iaP%matrix_struct, para_env=para_env, ncol_locals=ncol_locals, &
                            ncol_local=ncol_local, col_indices=col_indices, &
                            context=context, ncol_block=ncol_block, first_p_pos=first_p_pos, nrow_local=nrow_local)
      CALL get_blacs_info(context, my_process_row=my_prow, my_process_column=my_pcol, &
                          number_of_process_columns=num_pe_col, blacs2mpi=blacs2mpi)
      first_p_pos_col = first_p_pos(2)

      num_ab_pairs = SIZE(pair_list, 2)

      tag = 43

      DO col_local = 1, ncol_local
         col_global = col_indices(col_local)

         iocc = MAX(1, col_global - 1)/virtual + 1
         avirt = col_global - (iocc - 1)*virtual

         DO ab_counter = 1, num_ab_pairs

            my_a = pair_list(1, ab_counter)
            my_b = pair_list(2, ab_counter)

            IF (avirt /= my_b) CYCLE
            pcol = cp_fm_indxg2p((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
            IF (pcol /= my_pcol) CYCLE
            my_col_local = cp_fm_indxg2l((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)

            trace = ddot(nrow_local, fm_mat_S%local_data(:, my_col_local), 1, fm_work_iaP%local_data(:, col_local), 1)

            P_ab(ab_counter) = P_ab(ab_counter) &
                               - trace*sinh_over_x(0.5_dp*(Eigenval(my_a) - Eigenval(my_b))*omega)*omega*weight

            EXIT

         END DO
      END DO

      IF (num_pe_col > 1) THEN
         size_send_buffer = 0
         size_recv_buffer = 0
         DO proc_shift = 1, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            IF (ALLOCATED(index2send(pcol_send)%array)) &
               size_send_buffer = MAX(size_send_buffer, SIZE(index2send(pcol_send)%array))

            IF (ALLOCATED(index2recv(pcol_recv)%array)) &
               size_recv_buffer = MAX(size_recv_buffer, SIZE(index2recv(pcol_recv)%array))
         END DO

         ALLOCATE (buffer_send(nrow_local, size_send_buffer), buffer_recv(nrow_local, size_recv_buffer))

         DO proc_shift = 1, num_pe_col - 1
            pcol_send = MODULO(my_pcol + proc_shift, num_pe_col)
            pcol_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            ! Collect data and exchange
            send_size = 0
            IF (ALLOCATED(index2send(pcol_send)%array)) send_size = SIZE(index2send(pcol_send)%array)

            IF (send_size > 0) THEN
               counter = 0
               DO col_local = 1, ncol_local
                  col_global = col_indices(col_local)

                  iocc = MAX(1, col_global - 1)/virtual + 1
                  avirt = col_global - (iocc - 1)*virtual

                  DO ab_counter = 1, num_ab_pairs

                     my_a = pair_list(1, ab_counter)
                     my_b = pair_list(2, ab_counter)

                     IF (avirt /= my_b) CYCLE
                     pcol = cp_fm_indxg2p((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
                     IF (pcol /= pcol_send) CYCLE

                     counter = counter + 1

                     buffer_send(:, counter) = fm_work_iaP%local_data(:, col_local)

                     EXIT

                  END DO
               END DO
            END IF

            recv_size = 0
            IF (ALLOCATED(index2recv(pcol_recv)%array)) recv_size = SIZE(index2recv(pcol_recv)%array)
            IF (recv_size > 0) THEN
               CALL timeset(routineN//"_send", handle2)
               IF (send_size > 0) THEN
                  CALL mp_sendrecv(buffer_send(:, :send_size), blacs2mpi(my_prow, pcol_send), &
                                   buffer_recv(:, :recv_size), blacs2mpi(my_prow, pcol_recv), &
                                   para_env%group, tag)
               ELSE
                  CALL mp_recv(buffer_recv(:, :recv_size), blacs2mpi(my_prow, pcol_recv), tag, para_env%group)
               END IF
               CALL timestop(handle2)

               DO ab_counter = 1, num_ab_pairs
                  ! Collect the contributions of the matrix elements

                  my_a = pair_list(1, ab_counter)
                  my_b = pair_list(2, ab_counter)

                  trace = 0.0_dp

                  DO col_local = 1, SIZE(index2recv(pcol_recv)%array)
                     col_global = index2recv(pcol_recv)%array(col_local)

                     iocc = MAX(1, col_global - 1)/virtual + 1
                     avirt = col_global - (iocc - 1)*virtual
                     IF (avirt /= my_b) CYCLE
                     pcol = cp_fm_indxg2p((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)
                     IF (pcol /= my_pcol) CYCLE

                     my_col_local = cp_fm_indxg2l((iocc - 1)*virtual + my_a, ncol_block, 0, first_p_pos_col, num_pe_col)

                     trace = trace + ddot(nrow_local, fm_mat_S%local_data(:, my_col_local), 1, buffer_recv(:, col_local), 1)
                  END DO

                  P_ab(ab_counter) = P_ab(ab_counter) &
                                     - trace*sinh_over_x(0.5_dp*(Eigenval(my_a) - Eigenval(my_b))*omega)*omega*weight

               END DO
            ELSE IF (send_size > 0) THEN
               CALL timeset(routineN//"_send", handle2)
               CALL mp_send(buffer_send(:, :send_size), blacs2mpi(my_prow, pcol_send), tag, para_env%group)
               CALL timestop(handle2)
            END IF
         END DO
         IF (ALLOCATED(buffer_send)) DEALLOCATE (buffer_send)
         IF (ALLOCATED(buffer_recv)) DEALLOCATE (buffer_recv)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param index2send ...
!> \param index2recv ...
!> \param fm_mat_S ...
!> \param mat_S_3D ...
!> \param gd_homo ...
!> \param gd_virtual ...
!> \param mepos ...
! **************************************************************************************************
   SUBROUTINE redistribute_fm_mat_S(index2send, index2recv, fm_mat_S, mat_S_3D, gd_homo, gd_virtual, mepos)
      TYPE(one_dim_int_array), DIMENSION(0:), INTENT(IN) :: index2send
      TYPE(two_dim_int_array), DIMENSION(0:), INTENT(IN) :: index2recv
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: mat_S_3D
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_homo, gd_virtual
      INTEGER, DIMENSION(2), INTENT(IN)                  :: mepos

      CHARACTER(LEN=*), PARAMETER :: routineN = 'redistribute_fm_mat_S'

      INTEGER :: col_local, handle, my_a, my_homo, my_i, my_pcol, my_prow, my_virtual, nrow_local, &
         num_pe_col, proc_recv, proc_send, proc_shift, recv_size, send_size, size_recv_buffer, &
         size_send_buffer, tag
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_recv, buffer_send
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      tag = 46

      CALL get_blacs_info(fm_mat_S%matrix_struct%context, my_process_row=my_prow, my_process_column=my_pcol, &
                          number_of_process_columns=num_pe_col, blacs2mpi=blacs2mpi)

      CALL cp_fm_struct_get(fm_mat_S%matrix_struct, nrow_local=nrow_local, para_env=para_env)

      CALL get_group_dist(gd_homo, mepos(2), sizes=my_homo)
      CALL get_group_dist(gd_virtual, mepos(1), sizes=my_virtual)

      ALLOCATE (mat_S_3D(nrow_local, my_virtual, my_homo))

      IF (ALLOCATED(index2send(my_pcol)%array)) THEN
         DO col_local = 1, SIZE(index2send(my_pcol)%array)
            my_a = index2recv(my_pcol)%array(1, col_local)
            my_i = index2recv(my_pcol)%array(2, col_local)
            mat_S_3D(:, my_a, my_i) = fm_mat_S%local_data(:, index2send(my_pcol)%array(col_local))
         END DO
      END IF

      IF (num_pe_col > 1) THEN
         size_send_buffer = 0
         size_recv_buffer = 0
         DO proc_shift = 1, num_pe_col - 1
            proc_send = MODULO(my_pcol + proc_shift, num_pe_col)
            proc_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            send_size = 0
            IF (ALLOCATED(index2send(proc_send)%array)) send_size = SIZE(index2send(proc_send)%array)
            size_send_buffer = MAX(size_send_buffer, send_size)

            recv_size = 0
            IF (ALLOCATED(index2recv(proc_recv)%array)) recv_size = SIZE(index2recv(proc_recv)%array)
            size_recv_buffer = MAX(size_recv_buffer, recv_size)

         END DO

         ALLOCATE (buffer_send(nrow_local, size_send_buffer), buffer_recv(nrow_local, size_recv_buffer))

         DO proc_shift = 1, num_pe_col - 1
            proc_send = MODULO(my_pcol + proc_shift, num_pe_col)
            proc_recv = MODULO(my_pcol - proc_shift, num_pe_col)

            send_size = 0
            IF (ALLOCATED(index2send(proc_send)%array)) send_size = SIZE(index2send(proc_send)%array)
            IF (send_size > 0) THEN
               DO col_local = 1, send_size
                  buffer_send(:, col_local) = fm_mat_S%local_data(:, index2send(proc_send)%array(col_local))
               END DO
            END IF

            recv_size = 0
            IF (ALLOCATED(index2recv(proc_recv)%array)) recv_size = SIZE(index2recv(proc_recv)%array, 2)
            IF (recv_size > 0) THEN
               IF (send_size > 0) THEN
                  CALL mp_sendrecv(buffer_send(:, :send_size), blacs2mpi(my_prow, proc_send), &
                                   buffer_recv(:, :recv_size), blacs2mpi(my_prow, proc_recv), &
                                   para_env%group, tag)
               ELSE
                  CALL mp_recv(buffer_recv(:, :recv_size), blacs2mpi(my_prow, proc_recv), &
                               tag, para_env%group)
               END IF

               DO col_local = 1, recv_size
                  my_a = index2recv(proc_recv)%array(1, col_local)
                  my_i = index2recv(proc_recv)%array(2, col_local)
                  mat_S_3D(:, my_a, my_i) = buffer_recv(:, col_local)
               END DO
            ELSE IF (send_size > 0) THEN
               CALL mp_send(buffer_send(:, :send_size), blacs2mpi(my_prow, proc_send), &
                            tag, para_env%group)
            END IF

         END DO

         IF (ALLOCATED(buffer_send)) DEALLOCATE (buffer_send)
         IF (ALLOCATED(buffer_recv)) DEALLOCATE (buffer_recv)
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rpa_grad ...
!> \param mp2_env ...
!> \param para_env_sub ...
!> \param para_env ...
!> \param qs_env ...
!> \param gd_array ...
!> \param color_sub ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE rpa_grad_finalize(rpa_grad, mp2_env, para_env_sub, para_env, qs_env, gd_array, &
                                color_sub, do_ri_sos_laplace_mp2, homo, virtual)
      TYPE(rpa_grad_type), INTENT(INOUT)                 :: rpa_grad
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env_sub, para_env
      TYPE(qs_environment_type), INTENT(IN), POINTER     :: qs_env
      TYPE(group_dist_d1_type)                           :: gd_array
      INTEGER, INTENT(IN)                                :: color_sub
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_grad_finalize'

      INTEGER :: dimen_ia, dimen_RI, handle, iiB, ispin, my_group_L_end, my_group_L_size, &
         my_group_L_start, my_ia_end, my_ia_size, my_ia_start, my_P_end, my_P_size, my_P_start, &
         ngroup, nspins, pos_group, pos_sub, proc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pos_info
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: group_grid_2_mepos, mepos_2_grid_group
      REAL(KIND=dp)                                      :: my_scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Gamma_2D
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_G_P_ia, fm_PQ, fm_PQ_2, fm_PQ_half, &
                                                            fm_work_PQ, fm_work_PQ_2, fm_Y, &
                                                            operator_half
      TYPE(group_dist_d1_type)                           :: gd_array_new, gd_ia, gd_P, gd_P_new

      CALL timeset(routineN, handle)

      ! Release unnecessary matrices to save memory for next steps

      nspins = SIZE(rpa_grad%fm_Y)

      ! Scaling factor is required to scale the density matrices and the Gamma matrices later
      IF (do_ri_sos_laplace_mp2) THEN
         my_scale = mp2_env%scale_s
      ELSE
         my_scale = -mp2_env%ri_rpa%scale_rpa/(2.0_dp*pi)
         IF (mp2_env%ri_rpa%minimax_quad) my_scale = my_scale/2.0_dp
      END IF

      IF (do_ri_sos_laplace_mp2) THEN
         CALL sos_mp2_grad_finalize(rpa_grad%sos_mp2_work_occ, rpa_grad%sos_mp2_work_virt, &
                                    para_env, para_env_sub, homo, virtual, mp2_env)
      ELSE
         CALL rpa_grad_work_finalize(rpa_grad%rpa_work, mp2_env, homo, &
                                     virtual, para_env, para_env_sub)
      END IF

      CALL get_qs_env(qs_env, blacs_env=blacs_env)

      CALL cp_fm_get_info(rpa_grad%fm_Gamma_PQ, ncol_global=dimen_RI)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env)
      CALL cp_fm_create(fm_PQ, fm_struct)
      CALL cp_fm_create(fm_work_PQ, fm_struct)
      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         CALL cp_fm_create(fm_PQ_2, fm_struct)
      END IF
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_set_all(fm_PQ, 0.0_dp)

      ! We still have to left- and right multiply it with PQhalf
      CALL dereplicate_and_sum_fm(rpa_grad%fm_Gamma_PQ, fm_PQ)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      CALL prepare_redistribution(para_env, para_env_sub, ngroup, &
                                  group_grid_2_mepos, mepos_2_grid_group, pos_info=pos_info)

      ! Create fm_PQ_half
      CALL create_group_dist(gd_P, para_env_sub%num_pe, dimen_RI)
      CALL get_group_dist(gd_P, para_env_sub%mepos, my_P_start, my_P_end, my_P_size)

      CALL get_group_dist(gd_array, color_sub, my_group_L_start, my_group_L_end, my_group_L_size)

      CALL create_group_dist(gd_P_new, para_env%num_pe)
      CALL create_group_dist(gd_array_new, para_env%num_pe)

      DO proc = 0, para_env%num_pe - 1
         ! calculate position of the group
         pos_group = proc/para_env_sub%num_pe
         ! calculate position in the subgroup
         pos_sub = pos_info(proc)
         ! 1 -> rows, 2 -> cols
         CALL get_group_dist(gd_array, pos_group, gd_array_new, proc)
         CALL get_group_dist(gd_P, pos_sub, gd_P_new, proc)
      END DO

      DEALLOCATE (pos_info)
      CALL release_group_dist(gd_P)

      CALL array2fm(mp2_env%ri_grad%PQ_half, fm_PQ%matrix_struct, &
                    my_P_start, my_P_end, &
                    my_group_L_start, my_group_L_end, &
                    gd_P_new, gd_array_new, &
                    group_grid_2_mepos, para_env_sub%num_pe, ngroup, &
                    fm_PQ_half)

      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         CALL array2fm(mp2_env%ri_grad%operator_half, fm_PQ%matrix_struct, my_P_start, my_P_end, &
                       my_group_L_start, my_group_L_end, &
                       gd_P_new, gd_array_new, &
                       group_grid_2_mepos, para_env_sub%num_pe, ngroup, &
                       operator_half)
      END IF

      ! deallocate the info array
      CALL release_group_dist(gd_P_new)
      CALL release_group_dist(gd_array_new)

      IF (compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
! Finish Gamma_PQ
         CALL parallel_gemm(transa="N", transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_PQ, matrix_b=fm_PQ_half, beta=0.0_dp, &
                            matrix_c=fm_work_PQ)

         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=-my_scale, &
                            matrix_a=fm_PQ_half, matrix_b=fm_work_PQ, beta=0.0_dp, &
                            matrix_c=fm_PQ)

         CALL cp_fm_release(fm_work_PQ)
      ELSE
         CALL parallel_gemm(transa="N", transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_PQ, matrix_b=operator_half, beta=0.0_dp, &
                            matrix_c=fm_work_PQ)

         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=my_scale, &
                            matrix_a=operator_half, matrix_b=fm_work_PQ, beta=0.0_dp, &
                            matrix_c=fm_PQ)
         CALL cp_fm_release(operator_half)

         CALL cp_fm_create(fm_work_PQ_2, fm_PQ%matrix_struct, name="fm_Gamma_PQ_2")
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=-my_scale, &
                            matrix_a=fm_PQ_half, matrix_b=fm_work_PQ, beta=0.0_dp, &
                            matrix_c=fm_work_PQ_2)
         CALL cp_fm_to_fm(fm_work_PQ_2, fm_PQ_2)
         CALL cp_fm_geadd(1.0_dp, "T", fm_work_PQ_2, 1.0_dp, fm_PQ_2)
         CALL cp_fm_release(fm_work_PQ_2)
         CALL cp_fm_release(fm_work_PQ)
      END IF

      ALLOCATE (mp2_env%ri_grad%Gamma_PQ(my_P_size, my_group_L_size))
      CALL fm2array(mp2_env%ri_grad%Gamma_PQ, &
                    my_P_size, my_P_start, my_P_end, &
                    my_group_L_size, my_group_L_start, my_group_L_end, &
                    group_grid_2_mepos, mepos_2_grid_group, &
                    para_env_sub%num_pe, ngroup, &
                    fm_PQ)

      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         ALLOCATE (mp2_env%ri_grad%Gamma_PQ_2(my_P_size, my_group_L_size))
         CALL fm2array(mp2_env%ri_grad%Gamma_PQ_2, my_P_size, my_P_start, my_P_end, &
                       my_group_L_size, my_group_L_start, my_group_L_end, &
                       group_grid_2_mepos, mepos_2_grid_group, &
                       para_env_sub%num_pe, ngroup, &
                       fm_PQ_2)
      END IF

! Now, Gamma_Pia
      ALLOCATE (mp2_env%ri_grad%G_P_ia(my_group_L_size, nspins))
      DO ispin = 1, nspins
      DO iiB = 1, my_group_L_size
         NULLIFY (mp2_env%ri_grad%G_P_ia(iiB, ispin)%matrix)
      END DO
      END DO

      ! Redistribute the Y matrix
      DO ispin = 1, nspins
         ! Collect all data of columns for the own sub group locally
         CALL cp_fm_get_info(rpa_grad%fm_Y(ispin), ncol_global=dimen_ia)

         CALL get_qs_env(qs_env, blacs_env=blacs_env)

         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, template_fmstruct=fm_PQ_half%matrix_struct, nrow_global=dimen_ia)
         CALL cp_fm_create(fm_Y, fm_struct)
         CALL cp_fm_struct_release(fm_struct)
         CALL cp_fm_set_all(fm_Y, 0.0_dp)

         CALL dereplicate_and_sum_fm(rpa_grad%fm_Y(ispin), fm_Y)

         CALL cp_fm_create(fm_G_P_ia, fm_Y%matrix_struct)
         CALL cp_fm_set_all(fm_G_P_ia, 0.0_dp)

         CALL parallel_gemm(transa="N", transb="T", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=my_scale, &
                            matrix_a=fm_Y, matrix_b=fm_PQ_half, beta=0.0_dp, &
                            matrix_c=fm_G_P_ia)

         CALL cp_fm_release(fm_Y)

         CALL create_group_dist(gd_ia, para_env_sub%num_pe, dimen_ia)
         CALL get_group_dist(gd_ia, para_env_sub%mepos, my_ia_start, my_ia_end, my_ia_size)

         CALL fm2array(Gamma_2D, my_ia_size, my_ia_start, my_ia_end, &
                       my_group_L_size, my_group_L_start, my_group_L_end, &
                       group_grid_2_mepos, mepos_2_grid_group, &
                       para_env_sub%num_pe, ngroup, &
                       fm_G_P_ia)

         ! create the Gamma_ia_P in DBCSR style
         CALL create_dbcsr_gamma(Gamma_2D, homo(ispin), virtual(ispin), dimen_ia, para_env_sub, &
                                 my_ia_start, my_ia_end, my_group_L_size, gd_ia, &
                                 mp2_env%ri_grad%G_P_ia(:, ispin), mp2_env%ri_grad%mo_coeff_o(ispin)%matrix)

         CALL release_group_dist(gd_ia)

      END DO
      DEALLOCATE (rpa_grad%fm_Y)
      CALL cp_fm_release(fm_PQ_half)

      CALL timestop(handle)

   END SUBROUTINE rpa_grad_finalize

! **************************************************************************************************
!> \brief ...
!> \param sos_mp2_work_occ ...
!> \param sos_mp2_work_virt ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param homo ...
!> \param virtual ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE sos_mp2_grad_finalize(sos_mp2_work_occ, sos_mp2_work_virt, para_env, para_env_sub, homo, virtual, mp2_env)
      TYPE(sos_mp2_grad_work_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: sos_mp2_work_occ, sos_mp2_work_virt
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env, para_env_sub
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sos_mp2_grad_finalize'

      INTEGER                                            :: ab_counter, handle, ij_counter, ispin, &
                                                            itmp(2), my_a, my_b, my_B_size, my_i, &
                                                            my_j, nspins, pcol
      REAL(KIND=dp)                                      :: my_scale

      CALL timeset(routineN, handle)

      nspins = SIZE(sos_mp2_work_occ)
      my_scale = mp2_env%scale_s

      DO ispin = 1, nspins
         DO pcol = 0, SIZE(sos_mp2_work_occ(ispin)%index2send, 1) - 1
            IF (ALLOCATED(sos_mp2_work_occ(ispin)%index2send(pcol)%array)) &
               DEALLOCATE (sos_mp2_work_occ(ispin)%index2send(pcol)%array)
            IF (ALLOCATED(sos_mp2_work_occ(ispin)%index2send(pcol)%array)) &
               DEALLOCATE (sos_mp2_work_occ(ispin)%index2send(pcol)%array)
            IF (ALLOCATED(sos_mp2_work_virt(ispin)%index2recv(pcol)%array)) &
               DEALLOCATE (sos_mp2_work_virt(ispin)%index2recv(pcol)%array)
            IF (ALLOCATED(sos_mp2_work_virt(ispin)%index2recv(pcol)%array)) &
               DEALLOCATE (sos_mp2_work_virt(ispin)%index2recv(pcol)%array)
         END DO
         DEALLOCATE (sos_mp2_work_occ(ispin)%index2send, &
                     sos_mp2_work_occ(ispin)%index2recv, &
                     sos_mp2_work_virt(ispin)%index2send, &
                     sos_mp2_work_virt(ispin)%index2recv)
      END DO

      ! Sum P_ij and P_ab and redistribute them
      DO ispin = 1, nspins
         CALL mp_sum(sos_mp2_work_occ(ispin)%P, para_env%group)

         ALLOCATE (mp2_env%ri_grad%P_ij(ispin)%array(homo(ispin), homo(ispin)))
         mp2_env%ri_grad%P_ij(ispin)%array = 0.0_dp
         DO my_i = 1, homo(ispin)
            mp2_env%ri_grad%P_ij(ispin)%array(my_i, my_i) = my_scale*sos_mp2_work_occ(ispin)%P(my_i)
         END DO
         DO ij_counter = 1, SIZE(sos_mp2_work_occ(ispin)%pair_list, 2)
            my_i = sos_mp2_work_occ(ispin)%pair_list(1, ij_counter)
            my_j = sos_mp2_work_occ(ispin)%pair_list(2, ij_counter)

            mp2_env%ri_grad%P_ij(ispin)%array(my_i, my_j) = my_scale*sos_mp2_work_occ(ispin)%P(homo(ispin) + ij_counter)
            mp2_env%ri_grad%P_ij(ispin)%array(my_j, my_i) = my_scale*sos_mp2_work_occ(ispin)%P(homo(ispin) + ij_counter)
         END DO
         DEALLOCATE (sos_mp2_work_occ(ispin)%P, sos_mp2_work_occ(ispin)%pair_list)

         ! The first index of P_ab has to be distributed within the subgroups, so sum it up first and add the required elements later
         CALL mp_sum(sos_mp2_work_virt(ispin)%P, para_env%group)

         itmp = get_limit(virtual(ispin), para_env_sub%num_pe, para_env_sub%mepos)
         my_B_size = itmp(2) - itmp(1) + 1

         ALLOCATE (mp2_env%ri_grad%P_ab(ispin)%array(my_B_size, virtual(ispin)))
         mp2_env%ri_grad%P_ab(ispin)%array = 0.0_dp
         DO my_a = itmp(1), itmp(2)
            mp2_env%ri_grad%P_ab(ispin)%array(my_a - itmp(1) + 1, my_a) = my_scale*sos_mp2_work_virt(ispin)%P(my_a)
         END DO
         DO ab_counter = 1, SIZE(sos_mp2_work_virt(ispin)%pair_list, 2)
            my_a = sos_mp2_work_virt(ispin)%pair_list(1, ab_counter)
            my_b = sos_mp2_work_virt(ispin)%pair_list(2, ab_counter)

            IF (my_a >= itmp(1) .AND. my_a <= itmp(2)) mp2_env%ri_grad%P_ab(ispin)%array(my_a - itmp(1) + 1, my_b) = &
               my_scale*sos_mp2_work_virt(ispin)%P(virtual(ispin) + ab_counter)
            IF (my_b >= itmp(1) .AND. my_b <= itmp(2)) mp2_env%ri_grad%P_ab(ispin)%array(my_b - itmp(1) + 1, my_a) = &
               my_scale*sos_mp2_work_virt(ispin)%P(virtual(ispin) + ab_counter)
         END DO

         DEALLOCATE (sos_mp2_work_virt(ispin)%P, sos_mp2_work_virt(ispin)%pair_list)

      END DO
      DEALLOCATE (sos_mp2_work_occ, sos_mp2_work_virt)
      IF (nspins == 1) THEN
         mp2_env%ri_grad%P_ij(1)%array(:, :) = 2.0_dp*mp2_env%ri_grad%P_ij(1)%array
         mp2_env%ri_grad%P_ab(1)%array(:, :) = 2.0_dp*mp2_env%ri_grad%P_ab(1)%array
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rpa_work ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param para_env ...
!> \param para_env_sub ...
! **************************************************************************************************
   SUBROUTINE rpa_grad_work_finalize(rpa_work, mp2_env, homo, virtual, para_env, para_env_sub)
      TYPE(rpa_grad_work_type), INTENT(INOUT)            :: rpa_work
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env, para_env_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rpa_grad_work_finalize'

      INTEGER :: handle, ispin, itmp(2), my_a_end, my_a_size, my_a_start, my_B_end, my_B_size, &
         my_B_start, my_i_end, my_i_size, my_i_start, nspins, proc, proc_recv, proc_send, &
         proc_shift, recv_a_end, recv_a_size, recv_a_start, recv_end, recv_start, send_a_end, &
         send_a_size, send_a_start, send_end, send_start, size_recv_buffer, size_send_buffer
      REAL(KIND=dp)                                      :: my_scale
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_recv, buffer_send
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      TYPE(group_dist_d1_type)                           :: gd_a_sub, gd_virtual_sub

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      my_scale = mp2_env%ri_rpa%scale_rpa/(2.0_dp*pi)
      IF (mp2_env%ri_rpa%minimax_quad) my_scale = my_scale/2.0_dp

      CALL cp_fm_release(rpa_work%fm_mat_Q_copy)

      DO ispin = 1, nspins
      DO proc = 0, SIZE(rpa_work%index2send, 1) - 1
         IF (ALLOCATED(rpa_work%index2send(proc, ispin)%array)) DEALLOCATE (rpa_work%index2send(proc, ispin)%array)
         IF (ALLOCATED(rpa_work%index2recv(proc, ispin)%array)) DEALLOCATE (rpa_work%index2recv(proc, ispin)%array)
      END DO
      END DO
      DEALLOCATE (rpa_work%index2send, rpa_work%index2recv)

      DO ispin = 1, nspins
         CALL get_group_dist(rpa_work%gd_homo(ispin), rpa_work%mepos(2), my_i_start, my_i_end, my_i_size)
         CALL release_group_dist(rpa_work%gd_homo(ispin))

         ALLOCATE (mp2_env%ri_grad%P_ij(ispin)%array(homo(ispin), homo(ispin)))
         mp2_env%ri_grad%P_ij(ispin)%array = 0.0_dp
         mp2_env%ri_grad%P_ij(ispin)%array(my_i_start:my_i_end, :) = my_scale*rpa_work%P_ij(ispin)%array
         DEALLOCATE (rpa_work%P_ij(ispin)%array)
         CALL mp_sum(mp2_env%ri_grad%P_ij(ispin)%array, para_env%group)
         mp2_env%ri_grad%P_ij(ispin)%array(:, :) = 0.5_dp*(mp2_env%ri_grad%P_ij(ispin)%array + &
                                                           TRANSPOSE(mp2_env%ri_grad%P_ij(ispin)%array))

         itmp = get_limit(virtual(ispin), para_env_sub%num_pe, para_env_sub%mepos)
         my_B_start = itmp(1)
         my_B_end = itmp(2)
         my_B_size = my_B_end - my_B_start + 1

         ALLOCATE (mp2_env%ri_grad%P_ab(ispin)%array(my_B_size, virtual(ispin)))
         mp2_env%ri_grad%P_ab(ispin)%array = 0.0_dp

         CALL get_group_dist(rpa_work%gd_virtual(ispin), rpa_work%mepos(1), my_a_start, my_a_end, my_a_size)
         CALL release_group_dist(rpa_work%gd_virtual(ispin))
         ! This group dist contains the info which parts of Pab a process currently owns
         CALL create_group_dist(gd_a_sub, my_a_start, my_a_end, my_a_size, para_env_sub)
         ! This group dist contains the info which parts of Pab a process is supposed to own later
         CALL create_group_dist(gd_virtual_sub, para_env_sub%num_pe, virtual(ispin))

         ! Calculate local indices of the common range of own matrix and send process
         send_start = MAX(1, my_B_start - my_a_start + 1)
         send_end = MIN(my_a_size, my_B_end - my_a_start + 1)

         ! Same for recv process but with reverse positions
         recv_start = MAX(1, my_a_start - my_B_start + 1)
         recv_end = MIN(my_B_size, my_a_end - my_B_start + 1)

         mp2_env%ri_grad%P_ab(ispin)%array(recv_start:recv_end, :) = &
            my_scale*rpa_work%P_ab(ispin)%array(send_start:send_end, :)

         IF (para_env_sub%num_pe > 1) THEN
            size_send_buffer = 0
            size_recv_buffer = 0
            DO proc_shift = 1, para_env_sub%num_pe - 1
               proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
               proc_recv = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

               CALL get_group_dist(gd_virtual_sub, proc_send, send_a_start, send_a_end)
               CALL get_group_dist(gd_a_sub, proc_recv, recv_a_start, recv_a_end)

               ! Calculate local indices of the common range of own matrix and send process
               send_start = MAX(1, send_a_start - my_a_start + 1)
               send_end = MIN(my_a_size, send_a_end - my_a_start + 1)

               size_send_buffer = MAX(size_send_buffer, MAX(send_end - send_start + 1, 0))

               ! Same for recv process but with reverse positions
               recv_start = MAX(1, recv_a_start - my_B_start + 1)
               recv_end = MIN(my_B_size, recv_a_end - my_B_start + 1)

               size_recv_buffer = MAX(size_recv_buffer, MAX(recv_end - recv_start + 1, 0))
            END DO
            ALLOCATE (buffer_send(size_send_buffer, virtual(ispin)), buffer_recv(size_recv_buffer, virtual(ispin)))

            DO proc_shift = 1, para_env_sub%num_pe - 1
               proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
               proc_recv = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

               CALL get_group_dist(gd_virtual_sub, proc_send, send_a_start, send_a_end)
               CALL get_group_dist(gd_a_sub, proc_recv, recv_a_start, recv_a_end)

               ! Calculate local indices of the common range of own matrix and send process
               send_start = MAX(1, send_a_start - my_a_start + 1)
               send_end = MIN(my_a_size, send_a_end - my_a_start + 1)
               buffer_send(1:MAX(send_end - send_start + 1, 0), :) = rpa_work%P_ab(ispin)%array(send_start:send_end, :)

               ! Same for recv process but with reverse positions
               recv_start = MAX(1, recv_a_start - my_B_start + 1)
               recv_end = MIN(my_B_size, recv_a_end - my_B_start + 1)

               CALL mp_sendrecv(buffer_send(1:MAX(send_end - send_start + 1, 0), :), proc_send, &
                                buffer_recv(1:MAX(recv_end - recv_start + 1, 0), :), proc_recv, para_env_sub%group)

               mp2_env%ri_grad%P_ab(ispin)%array(recv_start:recv_end, :) = &
                  mp2_env%ri_grad%P_ab(ispin)%array(recv_start:recv_end, :) + &
                  my_scale*buffer_recv(1:MAX(recv_end - recv_start + 1, 0), :)

            END DO

            IF (ALLOCATED(buffer_send)) DEALLOCATE (buffer_send)
            IF (ALLOCATED(buffer_recv)) DEALLOCATE (buffer_recv)
         END IF
         DEALLOCATE (rpa_work%P_ab(ispin)%array)

         CALL release_group_dist(gd_a_sub)

         CALL cp_para_env_split(para_env_exchange, para_env, para_env_sub%mepos)
         CALL mp_sum(mp2_env%ri_grad%P_ab(ispin)%array, para_env_exchange%group)
         CALL cp_para_env_release(para_env_exchange)

         ! Symmetrize P_ab
         IF (para_env_sub%num_pe > 1) THEN

            mp2_env%ri_grad%P_ab(ispin)%array(:, my_B_start:my_B_end) = &
               0.5_dp*(mp2_env%ri_grad%P_ab(ispin)%array(:, my_B_start:my_B_end) &
                       + TRANSPOSE(mp2_env%ri_grad%P_ab(ispin)%array(:, my_B_start:my_B_end)))

            ALLOCATE (buffer_send(my_B_size, maxsize(gd_virtual_sub)))
            ALLOCATE (buffer_recv(my_B_size, maxsize(gd_virtual_sub)))

            DO proc_shift = 1, para_env_sub%num_pe - 1

               proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
               proc_recv = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

               CALL get_group_dist(gd_virtual_sub, proc_send, send_a_start, send_a_end, send_a_size)
               CALL get_group_dist(gd_virtual_sub, proc_recv, recv_a_start, recv_a_end, recv_a_size)

               buffer_send(:, :send_a_size) = mp2_env%ri_grad%P_ab(ispin)%array(:, send_a_start:send_a_end)
               CALL mp_sendrecv(buffer_send(:, :send_a_size), proc_send, &
                                buffer_recv(:, :recv_a_size), proc_recv, para_env_sub%group)

               mp2_env%ri_grad%P_ab(ispin)%array(:, recv_a_start:recv_a_end) = &
                  0.5_dp*(mp2_env%ri_grad%P_ab(ispin)%array(:, recv_a_start:recv_a_end) + buffer_send)

            END DO

            DEALLOCATE (buffer_send, buffer_recv)
         ELSE
            mp2_env%ri_grad%P_ab(ispin)%array(:, :) = 0.5_dp*(mp2_env%ri_grad%P_ab(ispin)%array + &
                                                              TRANSPOSE(mp2_env%ri_grad%P_ab(ispin)%array))
         END IF

         CALL release_group_dist(gd_virtual_sub)

      END DO
      DEALLOCATE (rpa_work%gd_homo, rpa_work%gd_virtual, rpa_work%P_ij, rpa_work%P_ab)
      IF (nspins == 1) THEN
         mp2_env%ri_grad%P_ij(1)%array(:, :) = 2.0_dp*mp2_env%ri_grad%P_ij(1)%array
         mp2_env%ri_grad%P_ab(1)%array(:, :) = 2.0_dp*mp2_env%ri_grad%P_ab(1)%array
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Dereplicate data from fm_sub and collect in fm_global, overlapping data will be added
!> \param fm_sub replicated matrix, all subgroups have the same size, will be release on output
!> \param fm_global global matrix, on output it will contain the sum of the replicated matrices redistributed
! **************************************************************************************************
   SUBROUTINE dereplicate_and_sum_fm(fm_sub, fm_global)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_sub, fm_global

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dereplicate_and_sum_fm'

      INTEGER :: col_local, elements2recv_col, elements2recv_row, elements2send_col, &
         elements2send_row, first_p_pos_global(2), first_p_pos_sub(2), handle, handle2, &
         mypcol_global, myprow_global, ncol_block_global, ncol_block_sub, ncol_local_global, &
         ncol_local_sub, npcol_global, npcol_sub, nprow_global, nprow_sub, nrow_block_global, &
         nrow_block_sub, nrow_local_global, nrow_local_sub, pcol_recv, pcol_send, proc_recv, &
         proc_send, proc_send_global, proc_shift, prow_recv, prow_send, row_local, tag
      INTEGER(int_8)                                     :: size_recv_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data2recv_col, data2recv_row, &
                                                            data2send_col, data2send_row, &
                                                            subgroup2mepos
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_global, col_indices_sub, &
                                                            row_indices_global, row_indices_sub
      INTEGER, DIMENSION(:, :), POINTER                  :: blacs2mpi_global, blacs2mpi_sub, &
                                                            mpi2blacs_global, mpi2blacs_sub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: recv_buffer_1D, send_buffer_1D
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: recv_buffer, send_buffer
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_exchange, para_env_sub
      TYPE(one_dim_int_array), ALLOCATABLE, DIMENSION(:) :: index2recv_col, index2recv_row, &
                                                            index2send_col, index2send_row

      CALL timeset(routineN, handle)

      tag = 1

      nprow_sub = fm_sub%matrix_struct%context%num_pe(1)
      npcol_sub = fm_sub%matrix_struct%context%num_pe(2)

      myprow_global = fm_global%matrix_struct%context%mepos(1)
      mypcol_global = fm_global%matrix_struct%context%mepos(2)
      nprow_global = fm_global%matrix_struct%context%num_pe(1)
      npcol_global = fm_global%matrix_struct%context%num_pe(2)

      CALL cp_fm_get_info(fm_sub, col_indices=col_indices_sub, row_indices=row_indices_sub, &
                          nrow_local=nrow_local_sub, ncol_local=ncol_local_sub)
      CALL cp_fm_struct_get(fm_sub%matrix_struct, para_env=para_env_sub, first_p_pos=first_p_pos_sub, &
                            nrow_block=nrow_block_sub, ncol_block=ncol_block_sub)
      CALL cp_fm_struct_get(fm_global%matrix_struct, first_p_pos=first_p_pos_global, nrow_block=nrow_block_global, &
                            ncol_block=ncol_block_global, para_env=para_env, &
                            col_indices=col_indices_global, row_indices=row_indices_global, &
                            nrow_local=nrow_local_global, ncol_local=ncol_local_global)
      CALL get_blacs_info(fm_sub%matrix_struct%context, blacs2mpi=blacs2mpi_sub, mpi2blacs=mpi2blacs_sub)
      CALL get_blacs_info(fm_global%matrix_struct%context, blacs2mpi=blacs2mpi_global, mpi2blacs=mpi2blacs_global)

      IF (para_env%num_pe /= para_env_sub%num_pe) THEN
         CALL timeset(routineN//"_exchange", handle2)
         NULLIFY (para_env_exchange)
         CALL cp_para_env_split(para_env_exchange, para_env, para_env_sub%mepos)
         CALL mp_sum(fm_sub%local_data, para_env_exchange%group)
         CALL cp_para_env_release(para_env_exchange)
         CALL timestop(handle2)
      END IF

      ALLOCATE (subgroup2mepos(0:para_env_sub%num_pe - 1))
      CALL mp_allgather(para_env%mepos, subgroup2mepos, para_env_sub%group)

      CALL timeset(routineN//"_data2", handle2)
      ! Create a map how much data has to be sent to what process coordinate, interchange rows and columns to transpose the matrices
      CALL get_elements2send(data2send_col, npcol_global, row_indices_sub, ncol_block_global, first_p_pos_global(2), index2send_col)
      CALL get_elements2send(data2send_row, nprow_global, col_indices_sub, nrow_block_global, first_p_pos_global(1), index2send_row)

      ! Create a map how much data has to be sent to what process coordinate, interchange rows and columns to transpose the matrices
      ! Do the reverse for the recieve processes
      CALL get_elements2send(data2recv_col, npcol_sub, row_indices_global, ncol_block_sub, first_p_pos_sub(2), index2recv_col)
      CALL get_elements2send(data2recv_row, nprow_sub, col_indices_global, nrow_block_sub, first_p_pos_sub(1), index2recv_row)
      CALL timestop(handle2)

      CALL timeset(routineN//"_local", handle2)
      ! Loop over local data and transpose
      prow_send = mpi2blacs_global(1, para_env%mepos)
      pcol_send = mpi2blacs_global(2, para_env%mepos)
      prow_recv = mpi2blacs_sub(1, para_env_sub%mepos)
      pcol_recv = mpi2blacs_sub(2, para_env_sub%mepos)
      elements2recv_col = data2recv_col(pcol_recv)
      elements2recv_row = data2recv_row(prow_recv)

!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(row_local,col_local) &
!$OMP                SHARED(elements2recv_col,elements2recv_row,recv_buffer,fm_global,&
!$OMP                       index2recv_col,index2recv_row,pcol_recv,prow_recv, &
!$OMP                       fm_sub,index2send_col,index2send_row,pcol_send,prow_send)
      DO col_local = 1, elements2recv_col
         DO row_local = 1, elements2recv_row
            fm_global%local_data(index2recv_col(pcol_recv)%array(col_local), &
                                 index2recv_row(prow_recv)%array(row_local)) &
               = fm_sub%local_data(index2send_col(pcol_send)%array(row_local), &
                                   index2send_row(prow_send)%array(col_local))
         END DO
      END DO
!$OMP    END PARALLEL DO
      CALL timestop(handle2)

      IF (para_env_sub%num_pe > 1) THEN
         size_send_buffer = 0_int_8
         size_recv_buffer = 0_int_8
         ! Loop over all processes in para_env_sub
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
            proc_recv = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

            proc_send_global = subgroup2mepos(proc_send)
            prow_send = mpi2blacs_global(1, proc_send_global)
            pcol_send = mpi2blacs_global(2, proc_send_global)
            elements2send_col = data2send_col(pcol_send)
            elements2send_row = data2send_row(prow_send)

            size_send_buffer = MAX(size_send_buffer, INT(elements2send_col, int_8)*elements2send_row)

            prow_recv = mpi2blacs_sub(1, proc_recv)
            pcol_recv = mpi2blacs_sub(2, proc_recv)
            elements2recv_col = data2recv_col(pcol_recv)
            elements2recv_row = data2recv_row(prow_recv)

            size_recv_buffer = MAX(size_recv_buffer, INT(elements2recv_col, int_8)*elements2recv_row)
         END DO
         ALLOCATE (send_buffer_1D(size_send_buffer), recv_buffer_1D(size_recv_buffer))

         ! Loop over all processes in para_env_sub
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
            proc_recv = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

            proc_send_global = subgroup2mepos(proc_send)
            prow_send = mpi2blacs_global(1, proc_send_global)
            pcol_send = mpi2blacs_global(2, proc_send_global)
            elements2send_col = data2send_col(pcol_send)
            elements2send_row = data2send_row(prow_send)

            CALL timeset(routineN//"_pack", handle2)
            ! Loop over local data and pack the buffer
            ! Transpose the matrix already
          send_buffer(1:elements2send_row, 1:elements2send_col) => send_buffer_1D(1:INT(elements2send_row, int_8)*elements2send_col)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(row_local,col_local) &
!$OMP                SHARED(elements2send_col,elements2send_row,send_buffer,fm_sub,&
!$OMP                       index2send_col,index2send_row,pcol_send,prow_send)
            DO row_local = 1, elements2send_col
               DO col_local = 1, elements2send_row
                  send_buffer(col_local, row_local) = &
                     fm_sub%local_data(index2send_col(pcol_send)%array(row_local), &
                                       index2send_row(prow_send)%array(col_local))
               END DO
            END DO
!$OMP    END PARALLEL DO
            CALL timestop(handle2)

            prow_recv = mpi2blacs_sub(1, proc_recv)
            pcol_recv = mpi2blacs_sub(2, proc_recv)
            elements2recv_col = data2recv_col(pcol_recv)
            elements2recv_row = data2recv_row(prow_recv)

            ! Send data
          recv_buffer(1:elements2recv_col, 1:elements2recv_row) => recv_buffer_1D(1:INT(elements2recv_row, int_8)*elements2recv_col)
            IF (SIZE(recv_buffer) > 0_int_8) THEN
            IF (SIZE(send_buffer) > 0_int_8) THEN
               CALL mp_sendrecv(send_buffer, proc_send, recv_buffer, proc_recv, para_env_sub%group, tag)
            ELSE
               CALL mp_recv(recv_buffer, proc_recv, tag, para_env_sub%group)
            END IF

            CALL timeset(routineN//"_unpack", handle2)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(row_local,col_local) &
!$OMP                SHARED(elements2recv_col,elements2recv_row,recv_buffer,fm_global,&
!$OMP                       index2recv_col,index2recv_row,pcol_recv,prow_recv)
            DO col_local = 1, elements2recv_col
               DO row_local = 1, elements2recv_row
                  fm_global%local_data(index2recv_col(pcol_recv)%array(col_local), &
                                       index2recv_row(prow_recv)%array(row_local)) &
                     = recv_buffer(col_local, row_local)
               END DO
            END DO
!$OMP    END PARALLEL DO
            CALL timestop(handle2)
            ELSE IF (SIZE(send_buffer) > 0_int_8) THEN
            CALL mp_send(send_buffer, proc_send, tag, para_env_sub%group)
            END IF
         END DO
      END IF

      DEALLOCATE (data2send_col, data2send_row, data2recv_col, data2recv_row)
      DO proc_shift = 0, npcol_global - 1
         DEALLOCATE (index2send_col(proc_shift)%array)
      END DO
      DO proc_shift = 0, npcol_sub - 1
         DEALLOCATE (index2recv_col(proc_shift)%array)
      END DO
      DO proc_shift = 0, nprow_global - 1
         DEALLOCATE (index2send_row(proc_shift)%array)
      END DO
      DO proc_shift = 0, nprow_sub - 1
         DEALLOCATE (index2recv_row(proc_shift)%array)
      END DO
      DEALLOCATE (index2send_col, index2recv_col, index2send_row, index2recv_row)

      CALL cp_fm_release(fm_sub)

      CALL timestop(handle)

   END SUBROUTINE dereplicate_and_sum_fm

! **************************************************************************************************
!> \brief ...
!> \param data2send ...
!> \param np_global ...
!> \param indices_sub ...
!> \param n_block_global ...
!> \param first_p_pos_global ...
!> \param index2send ...
! **************************************************************************************************
   SUBROUTINE get_elements2send(data2send, np_global, indices_sub, n_block_global, first_p_pos_global, index2send)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: data2send
      INTEGER, INTENT(IN)                                :: np_global
      INTEGER, DIMENSION(:), INTENT(IN)                  :: indices_sub
      INTEGER, INTENT(IN)                                :: n_block_global, first_p_pos_global
      TYPE(one_dim_int_array), ALLOCATABLE, &
         DIMENSION(:), INTENT(OUT)                       :: index2send

      INTEGER                                            :: i_global, i_local, proc

      ALLOCATE (data2send(0:np_global - 1))
      data2send = 0
      DO i_local = 1, SIZE(indices_sub)
         i_global = indices_sub(i_local)
         proc = cp_fm_indxg2p(i_global, n_block_global, 0, first_p_pos_global, np_global)
         data2send(proc) = data2send(proc) + 1
      END DO

      ALLOCATE (index2send(0:np_global - 1))
      DO proc = 0, np_global - 1
         ALLOCATE (index2send(proc)%array(data2send(proc)))
         ! We want to crash if there is an error
         index2send(proc)%array = -1
      END DO

      data2send = 0
      DO i_local = 1, SIZE(indices_sub)
         i_global = indices_sub(i_local)
         proc = cp_fm_indxg2p(i_global, n_block_global, 0, first_p_pos_global, np_global)
         data2send(proc) = data2send(proc) + 1
         index2send(proc)%array(data2send(proc)) = i_local
      END DO

   END SUBROUTINE

END MODULE rpa_grad
