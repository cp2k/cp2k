!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief communication routines to reshape / replicate / merge tall-and-skinny matrices.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_reshape_ops
   USE dbm_api,                         ONLY: &
        dbm_clear, dbm_distribution_col_dist, dbm_distribution_obj, dbm_distribution_row_dist, &
        dbm_finalize, dbm_get_col_block_sizes, dbm_get_distribution, dbm_get_name, &
        dbm_get_row_block_sizes, dbm_get_stored_coordinates, dbm_iterator, &
        dbm_iterator_blocks_left, dbm_iterator_next_block, dbm_iterator_start, dbm_iterator_stop, &
        dbm_put_block, dbm_reserve_blocks, dbm_type
   USE dbt_tas_base,                    ONLY: &
        dbt_repl_get_stored_coordinates, dbt_tas_blk_sizes, dbt_tas_clear, dbt_tas_create, &
        dbt_tas_distribution_new, dbt_tas_finalize, dbt_tas_get_stored_coordinates, dbt_tas_info, &
        dbt_tas_iterator_blocks_left, dbt_tas_iterator_next_block, dbt_tas_iterator_start, &
        dbt_tas_iterator_stop, dbt_tas_put_block, dbt_tas_reserve_blocks
   USE dbt_tas_global,                  ONLY: dbt_tas_blk_size_arb,&
                                              dbt_tas_blk_size_repl,&
                                              dbt_tas_dist_arb,&
                                              dbt_tas_dist_repl,&
                                              dbt_tas_distribution,&
                                              dbt_tas_rowcol_data
   USE dbt_tas_split,                   ONLY: colsplit,&
                                              dbt_tas_get_split_info,&
                                              rowsplit
   USE dbt_tas_types,                   ONLY: dbt_tas_distribution_type,&
                                              dbt_tas_iterator,&
                                              dbt_tas_split_info,&
                                              dbt_tas_type
   USE dbt_tas_util,                    ONLY: index_unique,&
                                              swap
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_environ,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_waitall
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_reshape_ops'

   PUBLIC :: &
      dbt_tas_merge, &
      dbt_tas_replicate, &
      dbt_tas_reshape

   TYPE dbt_buffer_type
      INTEGER :: nblock = -1
      INTEGER(KIND=int_8), DIMENSION(:, :), ALLOCATABLE :: indx
      REAL(dp), DIMENSION(:), ALLOCATABLE :: msg
      INTEGER :: endpos = -1
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief copy data (involves reshape)
!> \param matrix_in ...
!> \param matrix_out ...
!> \param summation whether matrix_out = matrix_out + matrix_in
!> \param transposed ...
!> \param move_data memory optimization: move data to matrix_out such that matrix_in is empty on return
!> \author Patrick Seewald
! **************************************************************************************************
   RECURSIVE SUBROUTINE dbt_tas_reshape(matrix_in, matrix_out, summation, transposed, move_data)
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_in, matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation, transposed, move_data

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbt_tas_reshape'

      INTEGER                                            :: bcount, handle, handle2, iproc, mp_comm, &
                                                            mynode, nblk, ndata, numnodes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: blks_to_allocate, index_recv
      INTEGER(KIND=int_8), DIMENSION(2)                  :: blk_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array
      INTEGER, DIMENSION(2)                              :: blk_size
      LOGICAL                                            :: move_prv, tr_in
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(dbt_buffer_type), ALLOCATABLE, DIMENSION(:)   :: buffer_recv, buffer_send
      TYPE(dbt_tas_iterator)                             :: iter
      TYPE(dbt_tas_split_info)                           :: info

      CALL timeset(routineN, handle)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL dbm_clear(matrix_out%matrix)
      ELSE
         CALL dbm_clear(matrix_out%matrix)
      END IF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      END IF

      IF (PRESENT(transposed)) THEN
         tr_in = transposed
      ELSE
         tr_in = .FALSE.
      END IF

      IF (.NOT. matrix_out%valid) THEN
         CPABORT("can not reshape into invalid matrix")
      END IF

      info = dbt_tas_info(matrix_in)
      mp_comm = info%mp_comm
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (buffer_send(0:numnodes - 1))
      ALLOCATE (buffer_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_send(0:numnodes - 1))
      ALLOCATE (num_entries_recv(0:numnodes - 1))
      ALLOCATE (num_entries_send(0:numnodes - 1))
      ALLOCATE (num_rec(0:2*numnodes - 1))
      ALLOCATE (num_send(0:2*numnodes - 1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))
      CALL dbt_tas_iterator_start(iter, matrix_in)

      CALL timeset(routineN//"_get_coord", handle2)
      DO WHILE (dbt_tas_iterator_blocks_left(iter))
         CALL dbt_tas_iterator_next_block(iter, blk_index(1), blk_index(2), &
                                          row_size=blk_size(1), col_size=blk_size(2))

         IF (tr_in) THEN
            CALL dbt_tas_get_stored_coordinates(matrix_out, blk_index(2), blk_index(1), iproc)
         ELSE
            CALL dbt_tas_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         END IF

         num_send(2*iproc) = num_send(2*iproc) + PRODUCT(blk_size)
         num_send(2*iproc + 1) = num_send(2*iproc + 1) + 1
      END DO
      CALL dbt_tas_iterator_stop(iter)
      CALL timestop(handle2)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      CALL timeset(routineN//"_buffer_fill", handle2)
      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL dbt_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc))

         CALL dbt_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc))

      END DO

      CALL dbt_tas_iterator_start(iter, matrix_in)
      DO WHILE (dbt_tas_iterator_blocks_left(iter))
         CALL dbt_tas_iterator_next_block(iter, blk_index(1), blk_index(2), block, &
                                          row_size=blk_size(1), col_size=blk_size(2))
         IF (tr_in) THEN
            CALL dbt_tas_get_stored_coordinates(matrix_out, blk_index(2), blk_index(1), iproc)
         ELSE
            CALL dbt_tas_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         END IF
         CALL dbt_buffer_add_block(buffer_send(iproc), blk_index, block, transposed=tr_in)
      END DO
      CALL dbt_tas_iterator_stop(iter)

      IF (move_prv) CALL dbt_tas_clear(matrix_in)

      CALL timestop(handle2)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL dbt_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

      DO iproc = 0, numnodes - 1
         CALL dbt_buffer_destroy(buffer_send(iproc))
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_buffer_obtain", handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL dbt_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = index_recv(:, :)
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      END DO

      CALL dbt_tas_reserve_blocks(matrix_out, blks_to_allocate(:, 1), blks_to_allocate(:, 2))
      DEALLOCATE (blks_to_allocate)

      DO iproc = 0, numnodes - 1
         ! First, we need to get the index to create block
         DO WHILE (dbt_buffer_blocks_left(buffer_recv(iproc)))
            CALL dbt_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index)
            CALL dbt_tas_blk_sizes(matrix_out, blk_index(1), blk_index(2), blk_size(1), blk_size(2))
            ALLOCATE (block(blk_size(1), blk_size(2)))
            CALL dbt_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index, block)
            CALL dbt_tas_put_block(matrix_out, blk_index(1), blk_index(2), block, summation=summation)
            DEALLOCATE (block)
         END DO
         CALL dbt_buffer_destroy(buffer_recv(iproc))
      END DO

      CALL timestop(handle2)

      CALL dbt_tas_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Replicate matrix_in such that each submatrix of matrix_out is an exact copy of matrix_in
!> \param matrix_in ...
!> \param info ...
!> \param matrix_out ...
!> \param nodata Don't copy data but create matrix_out
!> \param move_data memory optimization: move data to matrix_out such that matrix_in is empty on return
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_replicate(matrix_in, info, matrix_out, nodata, move_data)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix_in
      TYPE(dbt_tas_split_info), INTENT(IN)               :: info
      TYPE(dbt_tas_type), INTENT(OUT)                    :: matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata, move_data

      INTEGER                                            :: nblkcols, nblkrows
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, row_blk_size, &
                                                            row_dist
      TYPE(dbm_distribution_obj)                         :: dbm_dist
      TYPE(dbt_tas_dist_arb), TARGET                     :: dir_dist
      TYPE(dbt_tas_dist_repl), TARGET                    :: repl_dist

      CLASS(dbt_tas_distribution), ALLOCATABLE :: col_dist_obj, row_dist_obj
      CLASS(dbt_tas_rowcol_data), ALLOCATABLE :: row_bsize_obj, col_bsize_obj
      TYPE(dbt_tas_blk_size_repl), TARGET :: repl_blksize
      TYPE(dbt_tas_blk_size_arb), TARGET :: dir_blksize
      TYPE(dbt_tas_distribution_type) :: dist
      INTEGER :: mp_comm, numnodes, mynode
      TYPE(dbt_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: req_array, blks_to_allocate
      INTEGER, DIMENSION(2) :: blk_size
      INTEGER, DIMENSION(2) :: blk_index
      INTEGER(KIND=int_8), DIMENSION(2) :: blk_index_i8
      TYPE(dbm_iterator) :: iter
      INTEGER :: i, iproc, bcount, nblk
      INTEGER, DIMENSION(:), ALLOCATABLE :: iprocs
      LOGICAL :: nodata_prv, move_prv
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :) :: index_recv
      INTEGER :: ndata

      REAL(KIND=dp), DIMENSION(:, :), POINTER :: block

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_replicate'

      INTEGER :: handle, handle2

      NULLIFY (col_blk_size, row_blk_size)

      CALL timeset(routineN, handle)

      IF (PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      END IF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      END IF

      row_blk_size => dbm_get_row_block_sizes(matrix_in)
      col_blk_size => dbm_get_col_block_sizes(matrix_in)
      nblkrows = SIZE(row_blk_size)
      nblkcols = SIZE(col_blk_size)
      dbm_dist = dbm_get_distribution(matrix_in)
      row_dist => dbm_distribution_row_dist(dbm_dist)
      col_dist => dbm_distribution_col_dist(dbm_dist)

      mp_comm = info%mp_comm

      CALL mp_environ(numnodes, mynode, mp_comm)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm)

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         repl_dist = dbt_tas_dist_repl(row_dist, pdims(1), nblkrows, info%ngroup, info%pgrid_split_size)
         dir_dist = dbt_tas_dist_arb(col_dist, pdims(2), INT(nblkcols, KIND=int_8))
         repl_blksize = dbt_tas_blk_size_repl(row_blk_size, info%ngroup)
         dir_blksize = dbt_tas_blk_size_arb(col_blk_size)
         ALLOCATE (row_dist_obj, source=repl_dist)
         ALLOCATE (col_dist_obj, source=dir_dist)
         ALLOCATE (row_bsize_obj, source=repl_blksize)
         ALLOCATE (col_bsize_obj, source=dir_blksize)
      CASE (colsplit)
         dir_dist = dbt_tas_dist_arb(row_dist, pdims(1), INT(nblkrows, KIND=int_8))
         repl_dist = dbt_tas_dist_repl(col_dist, pdims(2), nblkcols, info%ngroup, info%pgrid_split_size)
         dir_blksize = dbt_tas_blk_size_arb(row_blk_size)
         repl_blksize = dbt_tas_blk_size_repl(col_blk_size, info%ngroup)
         ALLOCATE (row_dist_obj, source=dir_dist)
         ALLOCATE (col_dist_obj, source=repl_dist)
         ALLOCATE (row_bsize_obj, source=dir_blksize)
         ALLOCATE (col_bsize_obj, source=repl_blksize)
      END SELECT

      CALL dbt_tas_distribution_new(dist, mp_comm, row_dist_obj, col_dist_obj, split_info=info)
      CALL dbt_tas_create(matrix_out, TRIM(dbm_get_name(matrix_in))//" replicated", &
                          dist, row_bsize_obj, col_bsize_obj, own_dist=.TRUE.)

      IF (nodata_prv) THEN
         CALL dbt_tas_finalize(matrix_out)
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (buffer_send(0:numnodes - 1))
      ALLOCATE (buffer_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_send(0:numnodes - 1))
      ALLOCATE (num_entries_recv(0:numnodes - 1))
      ALLOCATE (num_entries_send(0:numnodes - 1))
      ALLOCATE (num_rec(0:2*numnodes - 1))
      ALLOCATE (num_send(0:2*numnodes - 1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))

      ALLOCATE (iprocs(info%ngroup))
      CALL dbm_iterator_start(iter, matrix_in)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, blk_index(1), blk_index(2), &
                                      row_size=blk_size(1), col_size=blk_size(2))
         CALL dbt_repl_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iprocs)
         DO i = 1, SIZE(iprocs)
            num_send(2*iprocs(i)) = num_send(2*iprocs(i)) + PRODUCT(blk_size)
            num_send(2*iprocs(i) + 1) = num_send(2*iprocs(i) + 1) + 1
         END DO
      END DO
      CALL dbm_iterator_stop(iter)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL dbt_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc))

         CALL dbt_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc))

      END DO

      CALL dbm_iterator_start(iter, matrix_in)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, blk_index(1), blk_index(2), block, &
                                      row_size=blk_size(1), col_size=blk_size(2))
         CALL dbt_repl_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iprocs)
         DO i = 1, SIZE(iprocs)
            CALL dbt_buffer_add_block(buffer_send(iprocs(i)), INT(blk_index, KIND=int_8), block)
         END DO
      END DO
      CALL dbm_iterator_stop(iter)

      IF (move_prv) CALL dbm_clear(matrix_in)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL dbt_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

      DO iproc = 0, numnodes - 1
         CALL dbt_buffer_destroy(buffer_send(iproc))
      END DO

      CALL timestop(handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL dbt_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = INT(index_recv(:, :))
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      END DO

      CALL dbm_reserve_blocks(matrix_out%matrix, blks_to_allocate(:, 1), blks_to_allocate(:, 2))
      DEALLOCATE (blks_to_allocate)

      DO iproc = 0, numnodes - 1
         ! First, we need to get the index to create block
         DO WHILE (dbt_buffer_blocks_left(buffer_recv(iproc)))
            CALL dbt_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8)
            CALL dbt_tas_blk_sizes(matrix_out, blk_index_i8(1), blk_index_i8(2), blk_size(1), blk_size(2))
            ALLOCATE (block(blk_size(1), blk_size(2)))
            CALL dbt_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8, block)
            CALL dbm_put_block(matrix_out%matrix, INT(blk_index_i8(1)), INT(blk_index_i8(2)), block)
            DEALLOCATE (block)
         END DO

         CALL dbt_buffer_destroy(buffer_recv(iproc))
      END DO

      CALL dbt_tas_finalize(matrix_out)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Merge submatrices of matrix_in to matrix_out by sum
!> \param matrix_out ...
!> \param matrix_in ...
!> \param summation ...
!> \param move_data memory optimization: move data to matrix_out such that matrix_in is empty on return
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_merge(matrix_out, matrix_in, summation, move_data)
      TYPE(dbm_type), INTENT(INOUT)                      :: matrix_out
      TYPE(dbt_tas_type), INTENT(INOUT)                  :: matrix_in
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation, move_data

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbt_tas_merge'

      INTEGER                                            :: bcount, handle, handle2, iproc, mp_comm, &
                                                            mynode, nblk, ndata, numnodes
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: index_recv
      INTEGER(KIND=int_8), DIMENSION(2)                  :: blk_index_i8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iprocs, num_blocks_recv, &
                                                            num_blocks_send, num_entries_recv, &
                                                            num_entries_send, num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: blks_to_allocate, blks_to_allocate_u, &
                                                            req_array
      INTEGER, DIMENSION(2)                              :: blk_index, blk_size
      INTEGER, DIMENSION(:), POINTER                     :: col_block_sizes, row_block_sizes
      LOGICAL                                            :: move_prv
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      TYPE(dbm_iterator)                                 :: iter
      TYPE(dbt_buffer_type), ALLOCATABLE, DIMENSION(:)   :: buffer_recv, buffer_send
      TYPE(dbt_tas_split_info)                           :: info

!!

      CALL timeset(routineN, handle)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL dbm_clear(matrix_out)
      ELSE
         CALL dbm_clear(matrix_out)
      END IF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      END IF

      info = dbt_tas_info(matrix_in)
      CALL dbt_tas_get_split_info(info, mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)

      ALLOCATE (buffer_send(0:numnodes - 1))
      ALLOCATE (buffer_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_send(0:numnodes - 1))
      ALLOCATE (num_entries_recv(0:numnodes - 1))
      ALLOCATE (num_entries_send(0:numnodes - 1))
      ALLOCATE (num_rec(0:2*numnodes - 1))
      ALLOCATE (num_send(0:2*numnodes - 1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))

      ALLOCATE (iprocs(info%ngroup))

      CALL dbm_iterator_start(iter, matrix_in%matrix)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, blk_index(1), blk_index(2), &
                                      row_size=blk_size(1), col_size=blk_size(2))
         CALL dbm_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         num_send(2*iproc) = num_send(2*iproc) + PRODUCT(blk_size)
         num_send(2*iproc + 1) = num_send(2*iproc + 1) + 1
      END DO
      CALL dbm_iterator_stop(iter)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL dbt_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc))

         CALL dbt_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc))

      END DO

      CALL dbm_iterator_start(iter, matrix_in%matrix)
      DO WHILE (dbm_iterator_blocks_left(iter))
         CALL dbm_iterator_next_block(iter, blk_index(1), blk_index(2), block, &
                                      row_size=blk_size(1), col_size=blk_size(2))
         CALL dbm_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         CALL dbt_buffer_add_block(buffer_send(iproc), INT(blk_index, KIND=int_8), block)
      END DO

      CALL dbm_iterator_stop(iter)

      IF (move_prv) CALL dbt_tas_clear(matrix_in)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL dbt_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

      DO iproc = 0, numnodes - 1
         CALL dbt_buffer_destroy(buffer_send(iproc))
      END DO

      CALL timestop(handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL dbt_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = INT(index_recv(:, :))
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      END DO

      CALL index_unique(blks_to_allocate, blks_to_allocate_u)

      CALL dbm_reserve_blocks(matrix_out, blks_to_allocate_u(:, 1), blks_to_allocate_u(:, 2))
      DEALLOCATE (blks_to_allocate, blks_to_allocate_u)

      DO iproc = 0, numnodes - 1
         ! First, we need to get the index to create block
         DO WHILE (dbt_buffer_blocks_left(buffer_recv(iproc)))
            CALL dbt_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8)
            row_block_sizes => dbm_get_row_block_sizes(matrix_out)
            col_block_sizes => dbm_get_col_block_sizes(matrix_out)
            blk_size(1) = row_block_sizes(INT(blk_index_i8(1)))
            blk_size(2) = col_block_sizes(INT(blk_index_i8(2)))
            ALLOCATE (block(blk_size(1), blk_size(2)))
            CALL dbt_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8, block)
            CALL dbm_put_block(matrix_out, INT(blk_index_i8(1)), INT(blk_index_i8(2)), block, summation=.TRUE.)
            DEALLOCATE (block)
         END DO
         CALL dbt_buffer_destroy(buffer_recv(iproc))
      END DO

      CALL dbm_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief get all indices from buffer
!> \param buffer ...
!> \param index ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_buffer_get_index(buffer, index)
      TYPE(dbt_buffer_type), INTENT(IN)                  :: buffer
      INTEGER(KIND=int_8), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: index

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_buffer_get_index'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(2)                              :: indx_shape

      CALL timeset(routineN, handle)

      indx_shape = SHAPE(buffer%indx) - [0, 1]
      ALLOCATE (INDEX(indx_shape(1), indx_shape(2)))
      INDEX(:, :) = buffer%indx(1:indx_shape(1), 1:indx_shape(2))
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief how many blocks left in iterator
!> \param buffer ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION dbt_buffer_blocks_left(buffer)
      TYPE(dbt_buffer_type), INTENT(IN)                  :: buffer
      LOGICAL                                            :: dbt_buffer_blocks_left

      dbt_buffer_blocks_left = buffer%endpos .LT. buffer%nblock
   END FUNCTION

! **************************************************************************************************
!> \brief Create block buffer for MPI communication.
!> \param buffer block buffer
!> \param nblock number of blocks
!> \param ndata total number of block entries
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_buffer_create(buffer, nblock, ndata)
      TYPE(dbt_buffer_type), INTENT(OUT)                 :: buffer
      INTEGER, INTENT(IN)                                :: nblock, ndata

      buffer%nblock = nblock
      buffer%endpos = 0
      ALLOCATE (buffer%msg(ndata))
      ALLOCATE (buffer%indx(nblock, 3))
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_buffer_destroy(buffer)
      TYPE(dbt_buffer_type), INTENT(INOUT)               :: buffer

      DEALLOCATE (buffer%msg)
      DEALLOCATE (buffer%indx)
      buffer%nblock = -1
      buffer%endpos = -1
   END SUBROUTINE dbt_buffer_destroy

! **************************************************************************************************
!> \brief insert a block into block buffer (at current iterator position)
!> \param buffer ...
!> \param index index of block
!> \param block ...
!> \param transposed ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_buffer_add_block(buffer, index, block, transposed)
      TYPE(dbt_buffer_type), INTENT(INOUT)               :: buffer
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(IN)      :: index
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed

      INTEGER                                            :: ndata, p, p_data
      INTEGER(KIND=int_8), DIMENSION(2)                  :: index_prv
      LOGICAL                                            :: tr

      IF (PRESENT(transposed)) THEN
         tr = transposed
      ELSE
         tr = .FALSE.
      END IF

      index_prv(:) = INDEX(:)
      IF (tr) THEN
         CALL swap(index_prv)
      END IF
      ndata = PRODUCT(SHAPE(block))

      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = INT(buffer%indx(p, 3))
      END IF

      IF (tr) THEN
         buffer%msg(p_data + 1:p_data + ndata) = RESHAPE(TRANSPOSE(block), [ndata])
      ELSE
         buffer%msg(p_data + 1:p_data + ndata) = RESHAPE(block, [ndata])
      END IF

      buffer%indx(p + 1, 1:2) = index_prv(:)
      IF (p > 0) THEN
         buffer%indx(p + 1, 3) = buffer%indx(p, 3) + INT(ndata, KIND=int_8)
      ELSE
         buffer%indx(p + 1, 3) = INT(ndata, KIND=int_8)
      END IF
      buffer%endpos = buffer%endpos + 1
   END SUBROUTINE

! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_buffer_get_next_block(buffer, ndata, index, block, advance_iter)
      TYPE(dbt_buffer_type), INTENT(INOUT)               :: buffer
      INTEGER, INTENT(OUT)                               :: ndata
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(OUT)     :: index
      REAL(dp), DIMENSION(:, :), INTENT(OUT), OPTIONAL   :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: advance_iter

      INTEGER                                            :: p, p_data
      LOGICAL                                            :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      END IF

      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = INT(buffer%indx(p, 3))
      END IF

      IF (p > 0) THEN
         ndata = INT(buffer%indx(p + 1, 3) - buffer%indx(p, 3))
      ELSE
         ndata = INT(buffer%indx(p + 1, 3))
      END IF
      INDEX(:) = buffer%indx(p + 1, 1:2)

      IF (PRESENT(block)) THEN
         block(:, :) = RESHAPE(buffer%msg(p_data + 1:p_data + ndata), SHAPE(block))
      END IF

      IF (do_advance) buffer%endpos = buffer%endpos + 1
   END SUBROUTINE

! **************************************************************************************************
!> \brief communicate buffer
!> \param mp_comm ...
!> \param buffer_recv ...
!> \param buffer_send ...
!> \param req_array ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      INTEGER, INTENT(IN)                                :: mp_comm
      TYPE(dbt_buffer_type), DIMENSION(0:), &
         INTENT(INOUT)                                   :: buffer_recv, buffer_send
      INTEGER, DIMENSION(:, :), INTENT(OUT)              :: req_array

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_tas_communicate_buffer'

      INTEGER                                            :: handle, iproc, mynode, numnodes, &
                                                            rec_counter, send_counter

      CALL timeset(routineN, handle)
      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (numnodes > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO iproc = 0, numnodes - 1
            IF (buffer_recv(iproc)%nblock > 0) THEN
               rec_counter = rec_counter + 1
               CALL mp_irecv(buffer_recv(iproc)%indx, iproc, mp_comm, req_array(rec_counter, 3), tag=4)
               CALL mp_irecv(buffer_recv(iproc)%msg, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO iproc = 0, numnodes - 1
            IF (buffer_send(iproc)%nblock > 0) THEN
               send_counter = send_counter + 1
               CALL mp_isend(buffer_send(iproc)%indx, iproc, mp_comm, req_array(send_counter, 1), tag=4)
               CALL mp_isend(buffer_send(iproc)%msg, iproc, mp_comm, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         IF (send_counter > 0) THEN
            CALL mp_waitall(req_array(1:send_counter, 1:2))
         END IF
         IF (rec_counter > 0) THEN
            CALL mp_waitall(req_array(1:rec_counter, 3:4))
         END IF

      ELSE
         IF (buffer_recv(0)%nblock > 0) THEN
            buffer_recv(0)%indx(:, :) = buffer_send(0)%indx(:, :)
            buffer_recv(0)%msg(:) = buffer_send(0)%msg(:)
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE

END MODULE
