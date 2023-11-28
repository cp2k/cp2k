!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      11.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_util

   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_indxg2p,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type
   USE mp2_types,                       ONLY: integ_mat_buffer_type
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_util'

   PUBLIC :: redist_fm

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_out ...
!> \param fm_mat_in ...
!> \param beta ...
!> \param nrow_secidx_in ...
!> \param ncol_secidx_in ...
!> \param nrow_secidx_out ...
!> \param ncol_secidx_out ...
!> \param unit_nr ...
!> \param reordering ...
! **************************************************************************************************
   SUBROUTINE redist_fm(fm_mat_out, fm_mat_in, beta, nrow_secidx_in, ncol_secidx_in, &
                        nrow_secidx_out, ncol_secidx_out, unit_nr, reordering)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_out, fm_mat_in
      REAL(kind=dp)                                      :: beta
      INTEGER, INTENT(IN)                                :: nrow_secidx_in, ncol_secidx_in, &
                                                            nrow_secidx_out, ncol_secidx_out
      INTEGER                                            :: unit_nr
      INTEGER, DIMENSION(4)                              :: reordering

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'redist_fm'

      INTEGER :: col_idx_loc, dummy_proc, handle, handle2, i_entry_rec, i_row_global, idx_col_out, &
         idx_row_out, ii, iproc, j_col_global, jj, ncol_block_in, ncol_block_out, ncol_local_in, &
         ncol_local_out, npcol, nprocs, nprow, nrow_block_in, nrow_block_out, nrow_local_in, &
         nrow_local_out, proc_send, row_idx_loc, send_pcol, send_prow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(4)                              :: indices_in
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_in, col_indices_out, &
                                                            row_indices_in, row_indices_out
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env_out
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req_array

      CALL timeset(routineN, handle)
      CALL timeset(routineN//"_1_setup", handle2)

      para_env_out => fm_mat_out%matrix_struct%para_env
      ! A_iajb
      ! We start by moving data from local parts of W_ijab to the full matrix A_iajb using buffers
      CALL cp_fm_get_info(matrix=fm_mat_out, &
                          nrow_local=nrow_local_out, &
                          ncol_local=ncol_local_out, &
                          row_indices=row_indices_out, &
                          col_indices=col_indices_out, &
                          nrow_block=nrow_block_out, &
                          ncol_block=ncol_block_out)

      nprow = fm_mat_out%matrix_struct%context%num_pe(1)
      npcol = fm_mat_out%matrix_struct%context%num_pe(2)

      ALLOCATE (num_entries_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env_out%num_pe - 1))

      num_entries_rec(:) = 0
      num_entries_send(:) = 0

      dummy_proc = 0

      CALL cp_fm_get_info(matrix=fm_mat_in, &
                          nrow_local=nrow_local_in, &
                          ncol_local=ncol_local_in, &
                          row_indices=row_indices_in, &
                          col_indices=col_indices_in, &
                          nrow_block=nrow_block_in, &
                          ncol_block=ncol_block_in)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A14,A10,T58,I10)') '|', 'Row number of ', fm_mat_out%name, &
            fm_mat_out%matrix_struct%nrow_global
         WRITE (unit_nr, '(T3,A1,T8,A17,A10,T58,I10)') '|', 'Column number of ', fm_mat_out%name, &
            fm_mat_out%matrix_struct%ncol_global

         WRITE (unit_nr, '(T3,A1,T8,A18,A10,T58,I10)') '|', 'Row block size of ', fm_mat_out%name, nrow_block_out
         WRITE (unit_nr, '(T3,A1,T8,A21,A10,T58,I10)') '|', 'Column block size of ', fm_mat_out%name, ncol_block_out

         WRITE (unit_nr, '(T3,A1,T8,A14,A10,T58,I10)') '|', 'Row number of ', fm_mat_in%name, &
            fm_mat_in%matrix_struct%nrow_global
         WRITE (unit_nr, '(T3,A1,T8,A17,A10,T58,I10)') '|', 'Column number of ', fm_mat_in%name, &
            fm_mat_in%matrix_struct%ncol_global

         WRITE (unit_nr, '(T3,A1,T8,A18,A10,T58,I10)') '|', 'Row block size of ', fm_mat_in%name, nrow_block_in
         WRITE (unit_nr, '(T3,A1,T8,A21,A10,T58,I10)') '|', 'Column block size of ', fm_mat_in%name, ncol_block_in
      END IF

      !MG to del: cf. mp2_r_2c.F and rpa_communication.F
      ! Use scalapack wrapper to find process index in A, where ijab -> iajb is needed
      ! To that end, we obtain the global index in A from the level indices
      indices_in(:) = 0
      DO row_idx_loc = 1, nrow_local_in
         indices_in(1) = (row_indices_in(row_idx_loc) - 1)/nrow_secidx_in + 1
         indices_in(2) = MOD(row_indices_in(row_idx_loc) - 1, nrow_secidx_in) + 1
         DO col_idx_loc = 1, ncol_local_in
            indices_in(3) = (col_indices_in(col_idx_loc) - 1)/ncol_secidx_in + 1
            indices_in(4) = MOD(col_indices_in(col_idx_loc) - 1, ncol_secidx_in) + 1

            idx_row_out = indices_in(reordering(2)) + (indices_in(reordering(1)) - 1)*nrow_secidx_out
            idx_col_out = indices_in(reordering(4)) + (indices_in(reordering(3)) - 1)*ncol_secidx_out

            send_prow = cp_fm_indxg2p(idx_row_out, nrow_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_mat_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2_comm_entry_nums", handle2)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A27)') '|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A18)') '|', 'Allocating buffers'
      END IF

      ! Buffers for entries and their indices
      ALLOCATE (buffer_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env_out%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 222.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 333.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%indx(num_entries_rec(iproc), 2))
         buffer_rec(iproc)%indx = 0

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%indx(num_entries_send(iproc), 2))
         buffer_send(iproc)%indx = 0

      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_4_buf_from_fmin_"//fm_mat_out%name, handle2)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A18,A10,A13)') '|', 'Writing data from ', fm_mat_in%name, ' into buffers'
      END IF

      ALLOCATE (entry_counter(0:para_env_out%num_pe - 1))
      entry_counter(:) = 0

      ! Now we can write the actual data and indices to the send-buffer
      DO row_idx_loc = 1, nrow_local_in
         indices_in(1) = (row_indices_in(row_idx_loc) - 1)/nrow_secidx_in + 1
         indices_in(2) = MOD(row_indices_in(row_idx_loc) - 1, nrow_secidx_in) + 1
         DO col_idx_loc = 1, ncol_local_in
            indices_in(3) = (col_indices_in(col_idx_loc) - 1)/ncol_secidx_in + 1
            indices_in(4) = MOD(col_indices_in(col_idx_loc) - 1, ncol_secidx_in) + 1

            idx_row_out = indices_in(reordering(2)) + (indices_in(reordering(1)) - 1)*nrow_secidx_out
            idx_col_out = indices_in(reordering(4)) + (indices_in(reordering(3)) - 1)*ncol_secidx_out

            send_prow = cp_fm_indxg2p(idx_row_out, nrow_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_mat_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)
            entry_counter(proc_send) = entry_counter(proc_send) + 1

            buffer_send(proc_send)%msg(entry_counter(proc_send)) = &
               fm_mat_in%local_data(row_idx_loc, col_idx_loc)

            buffer_send(proc_send)%indx(entry_counter(proc_send), 1) = idx_row_out
            buffer_send(proc_send)%indx(entry_counter(proc_send), 2) = idx_col_out

         END DO
      END DO

      ALLOCATE (req_array(1:para_env_out%num_pe, 4))

      CALL timestop(handle2)

      CALL timeset(routineN//"_5_comm_buffer", handle2)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A21)') '|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_mat_out%name, handle2)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A24,A10)') '|', 'Writing from buffers to ', fm_mat_out%name
      END IF

      ! MG to do: OMP parallel
      ! fill fm_mat_out with the entries from buffer_rec, i.e. buffer_rec are parts of W_ijab
      nprocs = para_env_out%num_pe

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_mat_out, nprocs, nrow_local_out,ncol_local_out, &
!$OMP          num_entries_rec, buffer_rec, row_indices_out, col_indices_out, beta ) &
!$OMP PRIVATE(iproc, i_entry_rec, ii, jj, &
!$OMP          i_row_global, j_col_global )
      DO iproc = 0, nprocs - 1
         DO i_entry_rec = 1, num_entries_rec(iproc)
            DO ii = 1, nrow_local_out

               i_row_global = row_indices_out(ii)

               DO jj = 1, ncol_local_out
                  j_col_global = col_indices_out(jj)

                  IF (i_row_global == buffer_rec(iproc)%indx(i_entry_rec, 1) .AND. &
                      j_col_global == buffer_rec(iproc)%indx(i_entry_rec, 2)) THEN
                     fm_mat_out%local_data(ii, jj) = fm_mat_out%local_data(ii, jj) + beta*buffer_rec(iproc)%msg(i_entry_rec)

                  END IF
               END DO
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7_cleanup", handle2)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A41)') '|', 'Starting cleanup of communication buffers'
      END IF

      !Clean up all the arrays from the communication process
      DO iproc = 0, para_env_out%num_pe - 1
         DEALLOCATE (buffer_rec(iproc)%msg)
         DEALLOCATE (buffer_rec(iproc)%indx)
         DEALLOCATE (buffer_send(iproc)%msg)
         DEALLOCATE (buffer_send(iproc)%indx)
      END DO
      DEALLOCATE (buffer_rec, buffer_send)
      DEALLOCATE (req_array)
      DEALLOCATE (entry_counter)
      DEALLOCATE (num_entries_rec, num_entries_send)

      CALL timestop(handle2)
      CALL timestop(handle)

   END SUBROUTINE redist_fm

END MODULE
