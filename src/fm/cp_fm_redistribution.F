!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Redistribution of full matrices and distributed arrays.
! **************************************************************************************************
MODULE cp_fm_redistribution
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE matrix_communication,            ONLY: exchange_matrix_buffers,&
                                              matrix_buffer_type
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_redistribution'
   PUBLIC :: cp_fm_permute

CONTAINS

! **************************************************************************************************
!> \brief Add a permuted or restricted four-index tensor stored as a full matrix.
!> \param fm_out destination; entries outside the mapped region remain unchanged
!> \param fm_in source; each matrix coordinate encodes two tensor indices
!> \param beta factor multiplying the source contribution
!> \param nrow_secidx_in size of the second tensor index in source rows
!> \param ncol_secidx_in size of the second tensor index in source columns
!> \param nrow_secidx_out stride of the first tensor index in destination rows
!> \param ncol_secidx_out stride of the first tensor index in destination columns
!> \param reordering permutation of the four source indices
!> \param row_offset zero-based destination row offset
!> \param col_offset zero-based destination column offset
!> \param first_index first retained value of each source index (default: one)
!> \param last_index last retained value of each source index (default: full extent)
!> \note The index map must be one-to-one and fit within the destination.
! **************************************************************************************************
   SUBROUTINE cp_fm_permute(fm_out, fm_in, beta, nrow_secidx_in, ncol_secidx_in, &
                            nrow_secidx_out, ncol_secidx_out, reordering, row_offset, col_offset, &
                            first_index, last_index)
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_out
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_in
      REAL(KIND=dp), INTENT(IN)                          :: beta
      INTEGER, INTENT(IN)                                :: nrow_secidx_in, ncol_secidx_in, &
                                                            nrow_secidx_out, ncol_secidx_out
      INTEGER, DIMENSION(4), INTENT(IN)                  :: reordering
      INTEGER, INTENT(IN), OPTIONAL                      :: row_offset, col_offset
      INTEGER, DIMENSION(4), INTENT(IN), OPTIONAL        :: first_index, last_index

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_fm_permute'

      INTEGER                                            :: handle, i, idx_col, idx_row, ientry, ii, &
                                                            iproc, j, jj, my_col_offset, &
                                                            my_row_offset, ncol_local, nprocs, &
                                                            nrow_local, pass
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: entry_counter, num_entries_rec, &
                                                            num_entries_send
      INTEGER, DIMENSION(4)                              :: first, indices, last, mapped
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(matrix_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      para_env => fm_out%matrix_struct%para_env
      nprocs = para_env%num_pe
      my_row_offset = 0
      my_col_offset = 0
      IF (PRESENT(row_offset)) my_row_offset = row_offset
      IF (PRESENT(col_offset)) my_col_offset = col_offset
      first = 1
      last = [fm_in%matrix_struct%nrow_global/nrow_secidx_in, nrow_secidx_in, &
              fm_in%matrix_struct%ncol_global/ncol_secidx_in, ncol_secidx_in]
      IF (PRESENT(first_index)) first = first_index
      IF (PRESENT(last_index)) last = last_index

      CALL cp_fm_get_info(fm_in, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)
      ALLOCATE (num_entries_send(0:nprocs - 1), num_entries_rec(0:nprocs - 1), &
                entry_counter(0:nprocs - 1))
      num_entries_send = 0

      ! Count and then pack using the same coordinate map.
      DO pass = 1, 2
         IF (pass == 2) THEN
            CALL para_env%alltoall(num_entries_send, num_entries_rec, 1)
            ALLOCATE (buffer_send(0:nprocs - 1), buffer_rec(0:nprocs - 1))
            DO iproc = 0, nprocs - 1
               ALLOCATE (buffer_send(iproc)%msg(num_entries_send(iproc)), &
                         buffer_send(iproc)%indx(num_entries_send(iproc), 2), &
                         buffer_rec(iproc)%msg(num_entries_rec(iproc)), &
                         buffer_rec(iproc)%indx(num_entries_rec(iproc), 2))
            END DO
            entry_counter = 0
         END IF

         DO i = 1, nrow_local
            indices(1) = (row_indices(i) - 1)/nrow_secidx_in + 1
            indices(2) = MOD(row_indices(i) - 1, nrow_secidx_in) + 1
            DO j = 1, ncol_local
               indices(3) = (col_indices(j) - 1)/ncol_secidx_in + 1
               indices(4) = MOD(col_indices(j) - 1, ncol_secidx_in) + 1
               IF (ANY(indices < first) .OR. ANY(indices > last)) CYCLE
               mapped = indices - first + 1
               idx_row = my_row_offset + mapped(reordering(2)) + &
                         (mapped(reordering(1)) - 1)*nrow_secidx_out
               idx_col = my_col_offset + mapped(reordering(4)) + &
                         (mapped(reordering(3)) - 1)*ncol_secidx_out
               iproc = fm_out%matrix_struct%context%blacs2mpi( &
                       fm_out%matrix_struct%g2p_row(idx_row), fm_out%matrix_struct%g2p_col(idx_col))
               IF (pass == 1) THEN
                  num_entries_send(iproc) = num_entries_send(iproc) + 1
               ELSE
                  entry_counter(iproc) = entry_counter(iproc) + 1
                  ientry = entry_counter(iproc)
                  buffer_send(iproc)%msg(ientry) = fm_in%local_data(i, j)
                  buffer_send(iproc)%indx(ientry, :) = [idx_row, idx_col]
               END IF
            END DO
         END DO
      END DO

      CALL exchange_matrix_buffers(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(fm_out, nprocs, num_entries_rec, buffer_rec, beta) &
!$OMP PRIVATE(iproc, ientry, ii, jj)
      DO iproc = 0, nprocs - 1
         DO ientry = 1, num_entries_rec(iproc)
            ii = fm_out%matrix_struct%g2l_row(buffer_rec(iproc)%indx(ientry, 1))
            jj = fm_out%matrix_struct%g2l_col(buffer_rec(iproc)%indx(ientry, 2))
            fm_out%local_data(ii, jj) = fm_out%local_data(ii, jj) + beta*buffer_rec(iproc)%msg(ientry)
         END DO
      END DO
!$OMP END PARALLEL DO
      CALL timestop(handle)
   END SUBROUTINE cp_fm_permute

END MODULE cp_fm_redistribution
