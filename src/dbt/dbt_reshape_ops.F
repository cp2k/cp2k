!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to reshape / redistribute tensors
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_reshape_ops
   #:include "dbt_macros.fypp"
   #:set maxdim = maxrank
   #:set ndims = range(2,maxdim+1)

   USE dbt_allocate_wrap, ONLY: allocate_any
   USE dbt_tas_base, ONLY: dbt_tas_copy, dbt_tas_get_info, dbt_tas_info
   USE dbt_block, ONLY: &
      block_nd, create_block, destroy_block, dbt_iterator_type, dbt_iterator_next_block, &
      dbt_iterator_blocks_left, dbt_iterator_start, dbt_iterator_stop, dbt_get_block, &
      dbt_reserve_blocks, dbt_put_block
   USE dbt_types, ONLY: dbt_blk_sizes, &
                        dbt_create, &
                        dbt_type, &
                        ndims_tensor, &
                        dbt_get_stored_coordinates, &
                        dbt_clear
   USE kinds, ONLY: default_string_length
   USE kinds, ONLY: dp, dp
   USE message_passing, ONLY: mp_alltoall, &
                              mp_environ, &
                              mp_irecv, &
                              mp_isend, &
                              mp_waitall

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_reshape_ops'

   PUBLIC :: dbt_reshape

   TYPE block_buffer_type
      INTEGER                                    :: ndim = -1
      INTEGER                                    :: nblock = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE      :: indx
      REAL(dp), DIMENSION(:), ALLOCATABLE        :: msg
      INTEGER                                    :: endpos = -1
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief copy data (involves reshape)
!>        tensor_out = tensor_out + tensor_in move_data memory optimization:
!>        transfer data from tensor_in to tensor_out s.t. tensor_in is empty on return
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_reshape(tensor_in, tensor_out, summation, move_data)

      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in, tensor_out
      LOGICAL, INTENT(IN), OPTIONAL                    :: summation
      LOGICAL, INTENT(IN), OPTIONAL                    :: move_data

      INTEGER                                            :: iproc, mp_comm, mynode, ndata, &
                                                            numnodes, bcount, nblk
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array, index_recv, blks_to_allocate
      TYPE(dbt_iterator_type)                        :: iter
      TYPE(block_nd)                                     :: blk_data
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(ndims_tensor(tensor_in))       :: blk_size, ind_nd, index
      LOGICAL :: found, summation_prv, move_prv

      IF (PRESENT(summation)) THEN
         summation_prv = summation
      ELSE
         summation_prv = .FALSE.
      END IF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      END IF

      CPASSERT(tensor_out%valid)

      IF (.NOT. summation_prv) CALL dbt_clear(tensor_out)

      mp_comm = tensor_in%pgrid%mp_comm_2d
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
      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size)
         CALL dbt_get_stored_coordinates(tensor_out, ind_nd, iproc)
         num_send(2*iproc) = num_send(2*iproc) + PRODUCT(blk_size)
         num_send(2*iproc + 1) = num_send(2*iproc + 1) + 1
      END DO
      CALL dbt_iterator_stop(iter)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc), &
                                  ndims_tensor(tensor_in))
         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc), &
                                  ndims_tensor(tensor_in))
      END DO
      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size)
         CALL dbt_get_block(tensor_in, ind_nd, blk_data, found)
         CPASSERT(found)
         CALL dbt_get_stored_coordinates(tensor_out, ind_nd, iproc)
         CALL block_buffer_add_anyd_block(buffer_send(iproc), ind_nd, blk_data)
         CALL destroy_block(blk_data)
      END DO
      CALL dbt_iterator_stop(iter)

      IF (move_prv) CALL dbt_clear(tensor_in)

      CALL dbt_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      DO iproc = 0, numnodes - 1
         CALL block_buffer_destroy(buffer_send(iproc))
      END DO

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, ndims_tensor(tensor_in)))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = index_recv(:, :)
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      END DO

      CALL dbt_reserve_blocks(tensor_out, blks_to_allocate)
      DEALLOCATE (blks_to_allocate)

      DO iproc = 0, numnodes - 1
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_anyd_block(buffer_recv(iproc), ndata, index)
            CALL dbt_blk_sizes(tensor_in, index, blk_size)
            ! create block
            CALL create_block(blk_data, blk_size)
            ! get actual block data
            CALL block_buffer_get_next_anyd_block(buffer_recv(iproc), ndata, index, blk_data)
            CALL dbt_put_block(tensor_out, index, blk_data, summation=summation)
            CALL destroy_block(blk_data)
         END DO
         CALL block_buffer_destroy(buffer_recv(iproc))
      END DO
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create block buffer for MPI communication.
!> \param buffer block buffer
!> \param nblock number of blocks
!> \param ndata total number of block entries
!> \param ndim number of dimensions
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_create(buffer, nblock, ndata, ndim)
      TYPE(block_buffer_type), INTENT(OUT) :: buffer
      INTEGER, INTENT(IN)                  :: nblock, ndata, ndim

      buffer%nblock = nblock
      buffer%endpos = 0
      buffer%ndim = ndim
      ALLOCATE (buffer%msg(ndata))
      ALLOCATE (buffer%indx(nblock, ndim + 1))
   END SUBROUTINE block_buffer_create

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_destroy(buffer)
      TYPE(block_buffer_type), INTENT(INOUT) :: buffer

      DEALLOCATE (buffer%msg)
      DEALLOCATE (buffer%indx)
      buffer%nblock = -1
      buffer%ndim = -1
      buffer%endpos = -1
   END SUBROUTINE block_buffer_destroy

! **************************************************************************************************
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION ndims_buffer(buffer)
      TYPE(block_buffer_type), INTENT(IN) :: buffer
      INTEGER                             :: ndims_buffer

      ndims_buffer = buffer%ndim
   END FUNCTION

! **************************************************************************************************
!> \brief insert a block into block buffer (at current iterator position)
!> \param index index of block
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_add_anyd_block(buffer, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                               :: index
      TYPE(block_nd), INTENT(IN)                  :: block

      CALL block_buffer_add_block(buffer, SIZE(block%blk), index, block%blk)
   END SUBROUTINE

! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_anyd_block(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(OUT)                        :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                              :: index
      TYPE(block_nd), INTENT(INOUT), OPTIONAL     :: block
      LOGICAL, INTENT(IN), OPTIONAL               :: advance_iter

      IF (PRESENT(block)) THEN
         CALL block_buffer_get_next_block(buffer, ndata, index, block%blk, advance_iter=advance_iter)
      ELSE
         CALL block_buffer_get_next_block(buffer, ndata, index, advance_iter=advance_iter)
      END IF
   END SUBROUTINE

! **************************************************************************************************
!> \brief Get all indices from buffer
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_get_index(buffer, index)
      TYPE(block_buffer_type), INTENT(IN)               :: buffer
      INTEGER, INTENT(OUT), DIMENSION(:, :), ALLOCATABLE :: index
      INTEGER, DIMENSION(2)                             :: indx_shape

      indx_shape = SHAPE(buffer%indx) - [0, 1]
      ALLOCATE (index, source=buffer%indx(1:indx_shape(1), 1:indx_shape(2)))
   END SUBROUTINE

! **************************************************************************************************
!> \brief how many blocks left in iterator
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION block_buffer_blocks_left(buffer)
      TYPE(block_buffer_type), INTENT(IN) :: buffer
      LOGICAL                             :: block_buffer_blocks_left

      block_buffer_blocks_left = buffer%endpos .LT. buffer%nblock
   END FUNCTION

! **************************************************************************************************
!> \brief communicate buffer
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      INTEGER, INTENT(IN)                    :: mp_comm
      TYPE(block_buffer_type), DIMENSION(0:), INTENT(INOUT) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(:, :), INTENT(OUT)               :: req_array

      INTEGER                                :: iproc, mynode, numnodes, rec_counter, &
                                                send_counter
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_communicate_buffer'

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

! **************************************************************************************************
!> \brief add block to buffer.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block(buffer, ndata, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)               :: buffer
      INTEGER, INTENT(IN)                                  :: ndata
      REAL(dp), DIMENSION(ndata), INTENT(IN)               :: block
      INTEGER, DIMENSION(ndims_buffer(buffer)), INTENT(IN) :: index
      INTEGER                                              :: p, ndims, p_data
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims + 1)
      END IF

      buffer%msg(p_data + 1:p_data + ndata) = block(:)
      buffer%indx(p + 1, 1:ndims) = index(:)
      IF (p > 0) THEN
         buffer%indx(p + 1, ndims + 1) = buffer%indx(p, ndims + 1) + ndata
      ELSE
         buffer%indx(p + 1, ndims + 1) = ndata
      END IF
      buffer%endpos = buffer%endpos + 1
   END SUBROUTINE

! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block(buffer, ndata, index, block, advance_iter)

      TYPE(block_buffer_type), INTENT(INOUT)                      :: buffer
      INTEGER, INTENT(OUT)                                        :: ndata
      REAL(dp), DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: block
      INTEGER, DIMENSION(ndims_buffer(buffer)), INTENT(OUT)       :: index
      INTEGER                                                     :: p, ndims, p_data
      LOGICAL, INTENT(IN), OPTIONAL                               :: advance_iter
      LOGICAL                                                     :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      END IF

      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims + 1)
      END IF
      IF (p > 0) THEN
         ndata = buffer%indx(p + 1, ndims + 1) - buffer%indx(p, ndims + 1)
      ELSE
         ndata = buffer%indx(p + 1, ndims + 1)
      END IF
      index(:) = buffer%indx(p + 1, 1:ndims)
      IF (PRESENT(block)) THEN
         CALL allocate_any(block, source=buffer%msg(p_data + 1:p_data + ndata))
      END IF

      IF (do_advance) buffer%endpos = buffer%endpos + 1
   END SUBROUTINE

END MODULE dbt_reshape_ops
