!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE rpa_communication,               ONLY: communicate_buffer
   USE util,                            ONLY: sort,&
                                              sort_unique
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_util'

   PUBLIC :: redist_fm, truncate_fm, fm_write_thresh, print_BSE_start_flag, &
             deallocate_matrices_bse, comp_eigvec_coeff_BSE, sort_excitations

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
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE redist_fm(fm_mat_out, fm_mat_in, beta, nrow_secidx_in, ncol_secidx_in, &
                        nrow_secidx_out, ncol_secidx_out, unit_nr, reordering, mp2_env)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_out, fm_mat_in
      REAL(kind=dp)                                      :: beta
      INTEGER, INTENT(IN)                                :: nrow_secidx_in, ncol_secidx_in, &
                                                            nrow_secidx_out, ncol_secidx_out
      INTEGER                                            :: unit_nr
      INTEGER, DIMENSION(4)                              :: reordering
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env

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

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T7,A14,A10,T71,I10)') 'BSE|', 'Row number of ', fm_mat_out%name, &
            fm_mat_out%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A4,T7,A17,A10,T71,I10)') 'BSE|', 'Column number of ', fm_mat_out%name, &
            fm_mat_out%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A4,T7,A18,A10,T71,I10)') 'BSE|', 'Row block size of ', fm_mat_out%name, nrow_block_out
         WRITE (unit_nr, '(T2,A4,T7,A21,A10,T71,I10)') 'BSE|', 'Column block size of ', fm_mat_out%name, ncol_block_out

         WRITE (unit_nr, '(T2,A4,T7,A14,A10,T71,I10)') 'BSE|', 'Row number of ', fm_mat_in%name, &
            fm_mat_in%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A4,T7,A17,A10,T71,I10)') 'BSE|', 'Column number of ', fm_mat_in%name, &
            fm_mat_in%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A4,T7,A18,A10,T71,I10)') 'BSE|', 'Row block size of ', fm_mat_in%name, nrow_block_in
         WRITE (unit_nr, '(T2,A4,T7,A21,A10,T71,I10)') 'BSE|', 'Column block size of ', fm_mat_in%name, ncol_block_in
      END IF

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
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A27)') 'BSE|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A18)') 'BSE|', 'Allocating buffers'
      END IF

      ! Buffers for entries and their indices
      ALLOCATE (buffer_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env_out%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 0.0_dp

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
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A18,A10,A13)') 'BSE|', 'Writing data from ', fm_mat_in%name, ' into buffers'
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
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A21)') 'BSE|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_mat_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A24,A10)') 'BSE|', 'Writing from buffers to ', fm_mat_out%name
      END IF

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
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A41)') 'BSE|', 'Starting cleanup of communication buffers'
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

! **************************************************************************************************
!> \brief Logic: Matrices have some dimension dimen_RI x nrow_in*ncol_in  for the incoming (untruncated) matrix
!>  and dimen_RI x nrow_out*ncol_out for the truncated matrix. The truncation is done by resorting the indices
!>  via parallel communication.
!> \param fm_mat_out ...
!> \param fm_mat_in ...
!> \param nrow_in ...
!> \param ncol_in ...
!> \param nrow_out ...
!> \param ncol_out ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param nrow_offset ...
!> \param ncol_offset ...
! **************************************************************************************************
   SUBROUTINE truncate_fm(fm_mat_out, fm_mat_in, nrow_in, ncol_in, &
                          nrow_out, ncol_out, unit_nr, mp2_env, &
                          nrow_offset, ncol_offset)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_out
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_in
      INTEGER                                            :: nrow_in, ncol_in, nrow_out, ncol_out, &
                                                            unit_nr
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN), OPTIONAL                      :: nrow_offset, ncol_offset

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'truncate_fm'

      INTEGER :: col_idx_loc, dummy_proc, handle, handle2, i_entry_rec, i_row_global, &
         idx_col_first, idx_col_in, idx_col_out, idx_col_sec, idx_row_in, ii, iproc, j_col_global, &
         jj, ncol_block_in, ncol_block_out, ncol_local_in, ncol_local_out, npcol, nprocs, nprow, &
         nrow_block_in, nrow_block_out, nrow_local_in, nrow_local_out, proc_send, row_idx_loc, &
         send_pcol, send_prow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_in, col_indices_out, &
                                                            row_indices_in, row_indices_out
      LOGICAL                                            :: correct_ncol, correct_nrow
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env_out
      TYPE(mp_request_type), DIMENSION(:, :), POINTER    :: req_array

      CALL timeset(routineN, handle)
      CALL timeset(routineN//"_1_setup", handle2)

      correct_nrow = .FALSE.
      correct_ncol = .FALSE.
      !In case of truncation in the occupied space, we need to correct the interval of indices
      IF (PRESENT(nrow_offset)) THEN
         correct_nrow = .TRUE.
      END IF
      IF (PRESENT(ncol_offset)) THEN
         correct_ncol = .TRUE.
      END IF

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

      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T7,A14,A10,T71,I10)') 'BSE|', 'Row number of ', fm_mat_out%name, &
            fm_mat_out%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A4,T7,A17,A10,T71,I10)') 'BSE|', 'Column number of ', fm_mat_out%name, &
            fm_mat_out%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A4,T7,A18,A10,T71,I10)') 'BSE|', 'Row block size of ', fm_mat_out%name, nrow_block_out
         WRITE (unit_nr, '(T2,A4,T7,A21,A10,T71,I10)') 'BSE|', 'Column block size of ', fm_mat_out%name, ncol_block_out

         WRITE (unit_nr, '(T2,A4,T7,A14,A10,T71,I10)') 'BSE|', 'Row number of ', fm_mat_in%name, &
            fm_mat_in%matrix_struct%nrow_global
         WRITE (unit_nr, '(T2,A4,T7,A17,A10,T71,I10)') 'BSE|', 'Column number of ', fm_mat_in%name, &
            fm_mat_in%matrix_struct%ncol_global

         WRITE (unit_nr, '(T2,A4,T7,A18,A10,T71,I10)') 'BSE|', 'Row block size of ', fm_mat_in%name, nrow_block_in
         WRITE (unit_nr, '(T2,A4,T7,A21,A10,T71,I10)') 'BSE|', 'Column block size of ', fm_mat_in%name, ncol_block_in
      END IF

      ! We find global indices in S with nrow_in and ncol_in for truncation
      DO row_idx_loc = 1, nrow_local_in
         idx_row_in = row_indices_in(row_idx_loc)
         DO col_idx_loc = 1, ncol_local_in
            idx_col_in = col_indices_in(col_idx_loc)

            idx_col_first = (idx_col_in - 1)/ncol_in + 1
            idx_col_sec = MOD(idx_col_in - 1, ncol_in) + 1
            ! If occupied orbitals are included, these have to be handled differently
            ! due to their reversed indexing
            IF (correct_nrow) THEN
               idx_col_first = idx_col_first - nrow_offset + 1
               IF (idx_col_first .LE. 0) CYCLE
            ELSE
               IF (idx_col_first > nrow_out) EXIT
            END IF
            IF (correct_ncol) THEN
               idx_col_sec = idx_col_sec - ncol_offset + 1
               IF (idx_col_sec .LE. 0) CYCLE
            ELSE
               IF (idx_col_sec > ncol_out) CYCLE
            END IF

            idx_col_out = idx_col_sec + (idx_col_first - 1)*ncol_out

            send_prow = cp_fm_indxg2p(idx_row_in, nrow_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_mat_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)

            num_entries_send(proc_send) = num_entries_send(proc_send) + 1

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_2_comm_entry_nums", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A27)') 'BSE|', 'Communicating entry numbers'
      END IF

      CALL para_env_out%alltoall(num_entries_send, num_entries_rec, 1)

      CALL timestop(handle2)

      CALL timeset(routineN//"_3_alloc_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A18)') 'BSE|', 'Allocating buffers'
      END IF

      ! Buffers for entries and their indices
      ALLOCATE (buffer_rec(0:para_env_out%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env_out%num_pe - 1))

      ! allocate data message and corresponding indices
      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_rec(iproc)%msg(num_entries_rec(iproc)))
         buffer_rec(iproc)%msg = 0.0_dp

      END DO

      DO iproc = 0, para_env_out%num_pe - 1

         ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)))
         buffer_send(iproc)%msg = 0.0_dp

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
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T7,A18,A10,A13)') 'BSE|', 'Writing data from ', fm_mat_in%name, ' into buffers'
      END IF

      ALLOCATE (entry_counter(0:para_env_out%num_pe - 1))
      entry_counter(:) = 0

      ! Now we can write the actual data and indices to the send-buffer
      DO row_idx_loc = 1, nrow_local_in
         idx_row_in = row_indices_in(row_idx_loc)
         DO col_idx_loc = 1, ncol_local_in
            idx_col_in = col_indices_in(col_idx_loc)

            idx_col_first = (idx_col_in - 1)/ncol_in + 1
            idx_col_sec = MOD(idx_col_in - 1, ncol_in) + 1
            ! If occupied orbitals are included, these have to be handled differently
            ! due to their reversed indexing
            IF (correct_nrow) THEN
               idx_col_first = idx_col_first - nrow_offset + 1
               IF (idx_col_first .LE. 0) CYCLE
            ELSE
               IF (idx_col_first > nrow_out) EXIT
            END IF
            IF (correct_ncol) THEN
               idx_col_sec = idx_col_sec - ncol_offset + 1
               IF (idx_col_sec .LE. 0) CYCLE
            ELSE
               IF (idx_col_sec > ncol_out) CYCLE
            END IF

            idx_col_out = idx_col_sec + (idx_col_first - 1)*ncol_out

            send_prow = cp_fm_indxg2p(idx_row_in, nrow_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(1), nprow)

            send_pcol = cp_fm_indxg2p(idx_col_out, ncol_block_out, dummy_proc, &
                                      fm_mat_out%matrix_struct%first_p_pos(2), npcol)

            proc_send = fm_mat_out%matrix_struct%context%blacs2mpi(send_prow, send_pcol)
            entry_counter(proc_send) = entry_counter(proc_send) + 1

            buffer_send(proc_send)%msg(entry_counter(proc_send)) = &
               fm_mat_in%local_data(row_idx_loc, col_idx_loc)
            !No need to create row_out, since it is identical to incoming
            !We dont change the RI index for any fm_mat_XX_BSE
            buffer_send(proc_send)%indx(entry_counter(proc_send), 1) = idx_row_in
            buffer_send(proc_send)%indx(entry_counter(proc_send), 2) = idx_col_out

         END DO
      END DO

      ALLOCATE (req_array(1:para_env_out%num_pe, 4))

      CALL timestop(handle2)

      CALL timeset(routineN//"_5_comm_buffer", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A21)') 'BSE|', 'Communicating buffers'
      END IF

      ! communicate the buffer
      CALL communicate_buffer(para_env_out, num_entries_rec, num_entries_send, buffer_rec, &
                              buffer_send, req_array)

      CALL timestop(handle2)

      CALL timeset(routineN//"_6_buffer_to_fmout"//fm_mat_out%name, handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A24,A10)') 'BSE|', 'Writing from buffers to ', fm_mat_out%name
      END IF

      ! fill fm_mat_out with the entries from buffer_rec, i.e. buffer_rec are parts of W_ijab
      nprocs = para_env_out%num_pe

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_mat_out, nprocs, nrow_local_out,ncol_local_out, &
!$OMP          num_entries_rec, buffer_rec, row_indices_out, col_indices_out) &
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
                     fm_mat_out%local_data(ii, jj) = fm_mat_out%local_data(ii, jj) + buffer_rec(iproc)%msg(i_entry_rec)

                  END IF
               END DO
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_7_cleanup", handle2)
      IF (unit_nr > 0 .AND. mp2_env%ri_g0w0%bse_debug_print) THEN
         WRITE (unit_nr, '(T2,A4,T10,A41)') 'BSE|', 'Starting cleanup of communication buffers'
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

   END SUBROUTINE truncate_fm

! **************************************************************************************************
!> \brief Debug function to write elements of a full matrix to file, if they are larger than a given threshold
!> \param fm_mat ...
!> \param thresh ...
!> \param header ...
!> \param unit_nr ...
!> \param abs_vals ...
! **************************************************************************************************
   SUBROUTINE fm_write_thresh(fm_mat, thresh, header, unit_nr, abs_vals)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      REAL(KIND=dp), INTENT(IN)                          :: thresh
      CHARACTER(LEN=*), INTENT(IN)                       :: header
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, OPTIONAL                                  :: abs_vals

      CHARACTER(LEN=*), PARAMETER :: my_footer = " | ENDING WRITING OF MATRIX", &
         routineN = 'fm_write_thresh'

      INTEGER                                            :: handle, i, j, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_abs_vals

      CALL timeset(routineN, handle)

      IF (PRESENT(abs_vals)) THEN
         my_abs_vals = abs_vals
      ELSE
         my_abs_vals = .FALSE.
      END IF

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      !MG to debug
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) header
      END IF
      IF (my_abs_vals) THEN
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               IF (ABS(fm_mat%local_data(i, j)) > thresh) THEN
                  WRITE (*, "(A7,T10,I5,T20,I5,T30,F13.5)") header, row_indices(i), col_indices(j), &
                     ABS(fm_mat%local_data(i, j))
               END IF
            END DO
         END DO
      ELSE
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               IF (ABS(fm_mat%local_data(i, j)) > thresh) THEN
                  WRITE (*, "(A7,T10,I5,T20,I5,T30,F13.5)") header, row_indices(i), col_indices(j), &
                     fm_mat%local_data(i, j)
               END IF
            END DO
         END DO
      END IF
      CALL fm_mat%matrix_struct%para_env%sync()
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) my_footer
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param bse_tda ...
!> \param bse_abba ...
!> \param mp2_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE print_BSE_start_flag(bse_tda, bse_abba, mp2_env, unit_nr)

      LOGICAL, INTENT(IN)                                :: bse_tda, bse_abba
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_BSE_start_flag'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T2,A79)') '*******************************************************************************'
         WRITE (unit_nr, '(T2,A79)') '**                                                                           **'
         WRITE (unit_nr, '(T2,A79)') '**           Bethe Salpeter equation (BSE) for excitation energies           **'
         IF (bse_tda .AND. bse_abba) THEN
            WRITE (unit_nr, '(T2,A79)') '**          solved with and without Tamm-Dancoff approximation (TDA)         **'
         ELSE IF (bse_tda) THEN
            WRITE (unit_nr, '(T2,A79)') '**                solved with Tamm-Dancoff approximation (TDA)               **'
         ELSE
            WRITE (unit_nr, '(T2,A79)') '**               solved without Tamm-Dancoff approximation (TDA)             **'
         END IF

         WRITE (unit_nr, '(T2,A79)') '**                                                                           **'
         WRITE (unit_nr, '(T2,A79)') '*******************************************************************************'
         WRITE (unit_nr, *) ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_trunc ...
!> \param fm_mat_S_ij_trunc ...
!> \param fm_mat_S_ab_trunc ...
!> \param fm_mat_Q_static_bse ...
!> \param fm_mat_Q_static_bse_gemm ...
! **************************************************************************************************
   SUBROUTINE deallocate_matrices_bse(fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, &
                                      fm_mat_S_trunc, fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, &
                                      fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm)

      TYPE(cp_fm_type), INTENT(INOUT) :: fm_mat_S_bar_ia_bse, fm_mat_S_bar_ij_bse, fm_mat_S_trunc, &
         fm_mat_S_ij_trunc, fm_mat_S_ab_trunc, fm_mat_Q_static_bse, fm_mat_Q_static_bse_gemm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_matrices_bse'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_release(fm_mat_S_bar_ia_bse)
      CALL cp_fm_release(fm_mat_S_bar_ij_bse)
      CALL cp_fm_release(fm_mat_S_trunc)
      CALL cp_fm_release(fm_mat_S_ij_trunc)
      CALL cp_fm_release(fm_mat_S_ab_trunc)
      CALL cp_fm_release(fm_mat_Q_static_bse)
      CALL cp_fm_release(fm_mat_Q_static_bse_gemm)

      CALL timestop(handle)

   END SUBROUTINE deallocate_matrices_bse

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_work ...
!> \param eig_vals ...
!> \param beta ...
!> \param gamma ...
!> \param do_transpose ...
! **************************************************************************************************
   SUBROUTINE comp_eigvec_coeff_BSE(fm_mat_work, eig_vals, beta, gamma, do_transpose)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: eig_vals
      REAL(KIND=dp), INTENT(IN)                          :: beta
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: gamma
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_transpose

      CHARACTER(LEN=*), PARAMETER :: routineN = 'comp_eigvec_coeff_BSE'

      INTEGER                                            :: handle, i_row_global, ii, j_col_global, &
                                                            jj, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: my_do_transpose
      REAL(KIND=dp)                                      :: coeff, my_gamma

      CALL timeset(routineN, handle)

      IF (PRESENT(gamma)) THEN
         my_gamma = gamma
      ELSE
         my_gamma = 2.0_dp
      END IF

      IF (PRESENT(do_transpose)) THEN
         my_do_transpose = do_transpose
      ELSE
         my_do_transpose = .FALSE.
      END IF

      CALL cp_fm_get_info(matrix=fm_mat_work, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      IF (my_do_transpose) THEN
         DO jj = 1, ncol_local
            j_col_global = col_indices(jj)
            DO ii = 1, nrow_local
               coeff = (eig_vals(j_col_global)**beta)/my_gamma
               fm_mat_work%local_data(ii, jj) = fm_mat_work%local_data(ii, jj)*coeff
            END DO
         END DO
      ELSE
         DO jj = 1, ncol_local
            DO ii = 1, nrow_local
               i_row_global = row_indices(ii)
               coeff = (eig_vals(i_row_global)**beta)/my_gamma
               fm_mat_work%local_data(ii, jj) = fm_mat_work%local_data(ii, jj)*coeff
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param idx_prim ...
!> \param idx_sec ...
!> \param eigvec_entries ...
! **************************************************************************************************
   SUBROUTINE sort_excitations(idx_prim, idx_sec, eigvec_entries)

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_prim, idx_sec
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sort_excitations'

      INTEGER                                            :: handle, ii, jj, kk, num_entries, num_mults
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: idx_prim_work, idx_sec_work, tmp_index
      LOGICAL                                            :: unique_entries
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigvec_entries_work

      CALL timeset(routineN, handle)

      num_entries = SIZE(idx_prim)

      ALLOCATE (tmp_index(num_entries))

      CALL sort(idx_prim, num_entries, tmp_index)

      ALLOCATE (idx_sec_work(num_entries))
      ALLOCATE (eigvec_entries_work(num_entries))

      DO ii = 1, num_entries
         idx_sec_work(ii) = idx_sec(tmp_index(ii))
         eigvec_entries_work(ii) = eigvec_entries(tmp_index(ii))
      END DO

      DEALLOCATE (tmp_index)
      DEALLOCATE (idx_sec)
      DEALLOCATE (eigvec_entries)

      CALL MOVE_ALLOC(idx_sec_work, idx_sec)
      CALL MOVE_ALLOC(eigvec_entries_work, eigvec_entries)

      !Now check for multiple entries in first idx to check necessity of sorting in second idx
      CALL sort_unique(idx_prim, unique_entries)
      IF (.NOT. unique_entries) THEN
         ALLOCATE (idx_prim_work(num_entries))
         idx_prim_work = idx_prim
         ! Find duplicate entries in idx_prim
         DO ii = 1, num_entries
            IF (idx_prim_work(ii) == 0) CYCLE
            num_mults = COUNT(idx_prim_work == idx_prim_work(ii))
            IF (num_mults > 1) THEN
               !Set all duplicate entries to 0
               idx_prim_work(ii:ii + num_mults - 1) = 0
               !Start sorting in secondary index
               ALLOCATE (idx_sec_work(num_mults))
               ALLOCATE (eigvec_entries_work(num_mults))
               idx_sec_work = idx_sec(ii:ii + num_mults - 1)
               eigvec_entries_work = eigvec_entries(ii:ii + num_mults - 1)
               ALLOCATE (tmp_index(num_mults))
               CALL sort(idx_sec_work, num_mults, tmp_index)

               !Now write newly sorted indices to original arrays
               DO kk = ii, ii + num_mults - 1
                  idx_sec(kk) = idx_sec_work(kk - ii +1 )
                  eigvec_entries(kk) = eigvec_entries_work(tmp_index(kk - ii + 1))
               END DO
               !Deallocate work arrays
               DEALLOCATE (tmp_index)
               DEALLOCATE (idx_sec_work)
               DEALLOCATE (eigvec_entries_work)
            END IF
            idx_prim_work(ii) = idx_prim(ii)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE sort_excitations

END MODULE
