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
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE matrix_communication,            ONLY: exchange_matrix_buffers,&
                                              matrix_buffer_type
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_redistribution'
   PUBLIC :: cp_array_to_fm, &
             prepare_subgroup_redistribution, &
             cp_fm_permute

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

! **************************************************************************************************
!> \brief Assemble a full matrix from distributed rectangular array tiles.
!> \param mat2D local tile; tiles on the FM communicator must not overlap
!> \param fm_struct structure defining the destination layout and communicator
!> \param my_start_row first global row of the tile
!> \param my_end_row last global row of the tile
!> \param my_start_col first global column of the tile
!> \param my_end_col last global column of the tile
!> \param fm_mat destination; entries outside the supplied tiles are zero
!> \param do_release_mat whether to deallocate the input tile (default: true)
! **************************************************************************************************
   SUBROUTINE cp_array_to_fm(mat2D, fm_struct, my_start_row, my_end_row, &
                             my_start_col, my_end_col, fm_mat, do_release_mat)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: mat2D
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER, INTENT(IN)                                :: my_start_row, my_end_row, my_start_col, &
                                                            my_end_col
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_mat
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_release_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_array_to_fm'

      INTEGER                                            :: handle, i, ientry, j, pass, peer
      INTEGER, ALLOCATABLE                               :: cursor(:), recv_count(:), send_count(:)
      LOGICAL                                            :: release_mat
      TYPE(matrix_buffer_type), ALLOCATABLE              :: buffer_rec(:), buffer_send(:)
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      para_env => fm_struct%para_env
      release_mat = .TRUE.
      IF (PRESENT(do_release_mat)) release_mat = do_release_mat
      CALL cp_fm_create(fm_mat, fm_struct, name="fm_mat", set_zero=.TRUE.)
      ALLOCATE (send_count(0:para_env%num_pe - 1), recv_count(0:para_env%num_pe - 1), &
                cursor(0:para_env%num_pe - 1))
      send_count = 0

      DO pass = 1, 2
         IF (pass == 2) THEN
            CALL para_env%alltoall(send_count, recv_count, 1)
            ALLOCATE (buffer_send(0:para_env%num_pe - 1), buffer_rec(0:para_env%num_pe - 1))
            DO peer = 0, para_env%num_pe - 1
               ALLOCATE (buffer_send(peer)%msg(send_count(peer)), buffer_send(peer)%indx(send_count(peer), 2), &
                         buffer_rec(peer)%msg(recv_count(peer)), buffer_rec(peer)%indx(recv_count(peer), 2))
            END DO
            cursor = 0
         END IF
         DO j = my_start_col, my_end_col
            DO i = my_start_row, my_end_row
               peer = fm_struct%context%blacs2mpi(fm_struct%g2p_row(i), fm_struct%g2p_col(j))
               IF (pass == 1) THEN
                  send_count(peer) = send_count(peer) + 1
               ELSE
                  cursor(peer) = cursor(peer) + 1
                  ientry = cursor(peer)
                  buffer_send(peer)%indx(ientry, :) = [i, j]
                  buffer_send(peer)%msg(ientry) = mat2D(i - my_start_row + 1, j - my_start_col + 1)
               END IF
            END DO
         END DO
      END DO
      IF (release_mat) DEALLOCATE (mat2D)
      CALL exchange_matrix_buffers(para_env, recv_count, send_count, buffer_rec, buffer_send)
      DO peer = 0, para_env%num_pe - 1
         DO ientry = 1, recv_count(peer)
            i = fm_struct%g2l_row(buffer_rec(peer)%indx(ientry, 1))
            j = fm_struct%g2l_col(buffer_rec(peer)%indx(ientry, 2))
            fm_mat%local_data(i, j) = buffer_rec(peer)%msg(ientry)
         END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE cp_array_to_fm

! **************************************************************************************************
!> \brief prepare array for redistribution
!> \param para_env ...
!> \param para_env_sub ...
!> \param ngroup ...
!> \param group_grid_2_mepos ...
!> \param mepos_2_grid_group ...
!> \param pos_info ...
! **************************************************************************************************
   SUBROUTINE prepare_subgroup_redistribution(para_env, para_env_sub, ngroup, &
                                              group_grid_2_mepos, mepos_2_grid_group, &
                                              pos_info)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: ngroup
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: group_grid_2_mepos, mepos_2_grid_group
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: pos_info

      INTEGER                                            :: i, pos_group, pos_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_pos_info

      ALLOCATE (my_pos_info(0:para_env%num_pe - 1))
      CALL para_env%allgather(para_env_sub%mepos, my_pos_info)

      ALLOCATE (group_grid_2_mepos(0:para_env_sub%num_pe - 1, 0:ngroup - 1))
      group_grid_2_mepos = 0
      ALLOCATE (mepos_2_grid_group(2, 0:para_env%num_pe - 1))
      mepos_2_grid_group = 0

      DO i = 0, para_env%num_pe - 1
         ! calculate position of the group
         pos_group = i/para_env_sub%num_pe
         ! calculate position in the subgroup
         pos_sub = my_pos_info(i)
         ! fill the map from the grid of groups to process
         group_grid_2_mepos(pos_sub, pos_group) = i
         ! and the opposite, from the global pos to the grid pos
         mepos_2_grid_group(1, i) = pos_sub
         mepos_2_grid_group(2, i) = pos_group
      END DO

      IF (PRESENT(pos_info)) CALL move_alloc(my_pos_info, pos_info)

   END SUBROUTINE prepare_subgroup_redistribution

END MODULE cp_fm_redistribution
