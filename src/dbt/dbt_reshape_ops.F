!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
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

   USE OMP_LIB, ONLY: omp_get_num_threads, omp_get_thread_num
   USE dbt_block, ONLY: &
      dbt_iterator_type, dbt_iterator_next_block, dbt_iterator_blocks_left, dbt_iterator_start, &
      dbt_iterator_stop, dbt_get_block_data, dbt_put_block_data, dbt_reserve_blocks
   USE dbt_types, ONLY: dbt_blk_sizes, &
                        dbt_type, &
                        ndims_tensor, &
                        dbt_get_stored_coordinates, &
                        dbt_clear
   USE kinds, ONLY: dp
   USE matrix_communication, ONLY: exchange_matrix_buffers, matrix_buffer_type
   USE message_passing, ONLY: mp_comm_type
   USE util, ONLY: get_limit

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_reshape_ops'

   PUBLIC :: dbt_reshape, dbt_reshape_release

   ! Communication buffers persist across calls such that pages are faulted in once rather than per
   ! call; per block, the index holds the data offset (relative to the peer's data) and nd-index.
   TYPE(matrix_buffer_type), SAVE :: buffer_send, buffer_recv

CONTAINS

! **************************************************************************************************
!> \brief copy data (involves reshape)
!>        tensor_out = tensor_out + tensor_in move_data memory optimization:
!>        transfer data from tensor_in to tensor_out s.t. tensor_in is empty on return
!> \param bounds crop tensor data while copying (as dbt_crop): start and end index per dimension
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbt_reshape(tensor_in, tensor_out, summation, move_data, bounds)

      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in, tensor_out
      LOGICAL, INTENT(IN), OPTIONAL                    :: summation
      LOGICAL, INTENT(IN), OPTIONAL                    :: move_data
      INTEGER, DIMENSION(2, ndims_tensor(tensor_in)), INTENT(IN), OPTIONAL :: bounds

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reshape'

      INTEGER                                            :: iproc, numnodes, nindex, &
                                                            handle, handle2, iblk, offset, ndata
      INTEGER, DIMENSION(2)                              :: block_range
      TYPE(dbt_iterator_type)                            :: iter
      CHARACTER(LEN=32)                                  :: name_pack, name_unpack
      INTEGER, DIMENSION(ndims_tensor(tensor_in))        :: blk_size, blk_offset, ind_nd
      INTEGER, DIMENSION(2, ndims_tensor(tensor_in))     :: bounds_prv
      LOGICAL :: found, summation_prv, move_prv, do_bounds, overlap, inside

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nblks_send_total, ndata_send_total, &
                                                            nblks_recv_total, ndata_recv_total, &
                                                            nblks_send_mythread, ndata_send_mythread
      TYPE(mp_comm_type) :: mp_comm

      CALL timeset(routineN, handle)

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

      do_bounds = PRESENT(bounds)
      IF (do_bounds) THEN
         bounds_prv(:, :) = bounds(:, :)
      ELSE
         bounds_prv(:, :) = 0
      END IF

      CPASSERT(tensor_out%valid)

      IF (.NOT. summation_prv) THEN
         CALL timeset(routineN//"_clear", handle2)
         CALL dbt_clear(tensor_out)
         CALL timestop(handle2)
      END IF

      mp_comm = tensor_in%pgrid%mp_comm_2d
      numnodes = mp_comm%num_pe
      nindex = ndims_tensor(tensor_in) + 1
      ! separate timers for blocks that need a reshape (index order differs from matrix layout)
      name_pack = routineN//"_pack"
      IF (permuted(tensor_in)) name_pack = routineN//"_pack_perm"
      name_unpack = routineN//"_unpack"
      IF (permuted(tensor_out)) name_unpack = routineN//"_unpack_perm"
      ALLOCATE (nblks_send_total(0:numnodes - 1), ndata_send_total(0:numnodes - 1), source=0)
      ALLOCATE (nblks_recv_total(0:numnodes - 1), ndata_recv_total(0:numnodes - 1), source=0)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor_in,tensor_out,summation,mp_comm,numnodes,nindex,name_pack,name_unpack) &
!$OMP SHARED(do_bounds,bounds_prv,buffer_send,buffer_recv) &
!$OMP SHARED(nblks_send_total,ndata_send_total,nblks_recv_total,ndata_recv_total) &
!$OMP PRIVATE(nblks_send_mythread,ndata_send_mythread) &
!$OMP PRIVATE(iter,ind_nd,blk_size,blk_offset,found,overlap,inside,iproc,handle2) &
!$OMP PRIVATE(block_range,offset,ndata,iblk)
      CALL timeset(routineN//"_count", handle2)
      ALLOCATE (nblks_send_mythread(0:numnodes - 1), ndata_send_mythread(0:numnodes - 1), source=0)

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size, blk_offset=blk_offset)
         IF (do_bounds) THEN
            CALL block_bounds(bounds_prv, blk_offset, blk_size, overlap, inside)
            IF (.NOT. overlap) CYCLE
         END IF
         CALL dbt_get_stored_coordinates(tensor_out, ind_nd, iproc)
         nblks_send_mythread(iproc) = nblks_send_mythread(iproc) + 1
         ndata_send_mythread(iproc) = ndata_send_mythread(iproc) + PRODUCT(blk_size)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP CRITICAL(omp_dbt_reshape)
      nblks_send_total(:) = nblks_send_total(:) + nblks_send_mythread(:)
      ndata_send_total(:) = ndata_send_total(:) + ndata_send_mythread(:)
      nblks_send_mythread(:) = nblks_send_total(:) ! current totals indicate slot for this thread
      ndata_send_mythread(:) = ndata_send_total(:)
!$OMP END CRITICAL(omp_dbt_reshape)
!$OMP BARRIER

!$OMP MASTER
      CALL mp_comm%alltoall(nblks_send_total, nblks_recv_total, 1)
      CALL mp_comm%alltoall(ndata_send_total, ndata_recv_total, 1)
!$OMP END MASTER
!$OMP BARRIER
      CALL timestop(handle2)

      CALL timeset(TRIM(name_pack), handle2)
!$OMP MASTER
      CALL grow_buffer(buffer_send, ndata_send_total, nblks_send_total, nindex)
      CALL grow_buffer(buffer_recv, ndata_recv_total, nblks_recv_total, nindex)
!$OMP END MASTER
!$OMP BARRIER

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size, blk_offset=blk_offset)
         inside = .TRUE.
         IF (do_bounds) THEN
            CALL block_bounds(bounds_prv, blk_offset, blk_size, overlap, inside)
            IF (.NOT. overlap) CYCLE
         END IF
         CALL dbt_get_stored_coordinates(tensor_out, ind_nd, iproc)
         ! insert block data
         ndata = PRODUCT(blk_size)
         ndata_send_mythread(iproc) = ndata_send_mythread(iproc) - ndata
         offset = buffer_send%msg_offset(iproc) + ndata_send_mythread(iproc)
         CALL dbt_get_block_data(tensor_in, ind_nd, blk_size, buffer_send%msg(offset + 1:offset + ndata), found)
         CPASSERT(found)
         IF (.NOT. inside) CALL block_crop(bounds_prv, blk_offset, blk_size, buffer_send%msg(offset + 1:offset + ndata))
         ! insert block index
         nblks_send_mythread(iproc) = nblks_send_mythread(iproc) - 1
         iblk = buffer_send%indx_offset(iproc) + nblks_send_mythread(iproc) + 1
         buffer_send%indx(2:, iblk) = ind_nd(:)
         ! Payload offsets in the metadata are relative to the peer's segment.
         buffer_send%indx(1, iblk) = ndata_send_mythread(iproc)
      END DO
      CALL dbt_iterator_stop(iter)
      DEALLOCATE (nblks_send_mythread, ndata_send_mythread)
!$OMP BARRIER
      CALL timestop(handle2)

!$OMP MASTER
      IF (numnodes > 1) THEN
         CALL exchange_matrix_buffers(mp_comm, buffer_recv, buffer_send)
      ELSE
         ! A single rank receives what it sends: the buffers trade places.
         CALL swap_buffers(buffer_send, buffer_recv)
      END IF
!$OMP END MASTER
!$OMP BARRIER

      CALL timeset(routineN//"_reserve", handle2)
      block_range = get_limit(buffer_recv%indx_offset(numnodes), omp_get_num_threads(), omp_get_thread_num())
      CALL dbt_reserve_blocks(tensor_out, TRANSPOSE(buffer_recv%indx(2:, block_range(1):block_range(2))))
!$OMP BARRIER
      CALL timestop(handle2)

      CALL timeset(TRIM(name_unpack), handle2)
      DO iproc = 0, numnodes - 1
!$OMP DO
         DO iblk = buffer_recv%indx_offset(iproc) + 1, buffer_recv%indx_offset(iproc + 1)
            ind_nd(:) = buffer_recv%indx(2:, iblk)
            CALL dbt_blk_sizes(tensor_out, ind_nd, blk_size)
            offset = buffer_recv%msg_offset(iproc) + buffer_recv%indx(1, iblk)
            ndata = PRODUCT(blk_size)
            CALL dbt_put_block_data(tensor_out, ind_nd, blk_size, buffer_recv%msg(offset + 1:offset + ndata), &
                                    summation=summation)
         END DO
!$OMP END DO
      END DO
      CALL timestop(handle2)
!$OMP END PARALLEL

      DEALLOCATE (nblks_recv_total, ndata_recv_total)
      DEALLOCATE (nblks_send_total, ndata_send_total)

      IF (move_prv) THEN
         CALL timeset(routineN//"_clear", handle2)
         CALL dbt_clear(tensor_in)
         CALL timestop(handle2)
      END IF

      CALL timestop(handle)
   END SUBROUTINE dbt_reshape

! **************************************************************************************************
!> \brief Release the communication buffers that persist across calls of dbt_reshape.
!> \author Hans Pabst
! **************************************************************************************************
   SUBROUTINE dbt_reshape_release()
      CALL release_buffer(buffer_send)
      CALL release_buffer(buffer_recv)
   END SUBROUTINE dbt_reshape_release

! **************************************************************************************************
!> \brief Whether the tensor's index order differs from its matrix layout (blocks need a reshape).
!> \author Hans Pabst
! **************************************************************************************************
   FUNCTION permuted(tensor)
      TYPE(dbt_type), INTENT(IN)                            :: tensor
      LOGICAL                                               :: permuted

      INTEGER                                               :: i

      ASSOCIATE (map1_2d => tensor%nd_index_blk%map1_2d, map2_2d => tensor%nd_index_blk%map2_2d)
         permuted = .NOT. ALL([map1_2d, map2_2d] == [(i, i=1, SIZE(map1_2d) + SIZE(map2_2d))])
      END ASSOCIATE
   END FUNCTION permuted

! **************************************************************************************************
!> \brief Whether a block (element offset and size per dimension) overlaps the bounds, and whether
!>        it lies entirely inside them.
!> \author Hans Pabst
! **************************************************************************************************
   PURE SUBROUTINE block_bounds(bounds, blk_offset, blk_size, overlap, inside)
      INTEGER, DIMENSION(:, :), INTENT(IN)                  :: bounds
      INTEGER, DIMENSION(:), INTENT(IN)                     :: blk_offset, blk_size
      LOGICAL, INTENT(OUT)                                  :: overlap, inside

      overlap = ALL(blk_offset <= bounds(2, :) .AND. bounds(1, :) <= blk_offset + blk_size - 1)
      inside = ALL(bounds(1, :) <= blk_offset .AND. blk_offset + blk_size - 1 <= bounds(2, :))
   END SUBROUTINE block_bounds

! **************************************************************************************************
!> \brief Zero the elements of a block (column-major data) outside the bounds (as dbt_crop).
!> \author Hans Pabst
! **************************************************************************************************
   PURE SUBROUTINE block_crop(bounds, blk_offset, blk_size, data)
      INTEGER, DIMENSION(:, :), INTENT(IN)                  :: bounds
      INTEGER, DIMENSION(:), INTENT(IN)                     :: blk_offset, blk_size
      REAL(dp), DIMENSION(:), INTENT(INOUT)                 :: data

      INTEGER                                               :: d, i
      INTEGER, DIMENSION(SIZE(blk_size))                    :: hi, idx, lo

      lo(:) = MAX(bounds(1, :) - blk_offset + 1, 1)
      hi(:) = MIN(bounds(2, :) - blk_offset + 1, blk_size)
      idx(:) = 1
      DO i = 1, SIZE(data)
         IF (ANY(idx < lo .OR. hi < idx)) data(i) = 0.0_dp
         DO d = 1, SIZE(idx) ! next column-major index
            idx(d) = idx(d) + 1
            IF (idx(d) <= blk_size(d)) EXIT
            idx(d) = 1
         END DO
      END DO
   END SUBROUTINE block_crop

! **************************************************************************************************
!> \brief Lay out a persistent buffer for the given per-peer counts (as allocate_matrix_buffers),
!>        growing its arrays only when they are too small.
!> \author Hans Pabst
! **************************************************************************************************
   SUBROUTINE grow_buffer(buffer, ndata, nblks, nindex)
      TYPE(matrix_buffer_type), INTENT(INOUT)               :: buffer
      INTEGER, DIMENSION(0:), INTENT(IN)                    :: ndata, nblks
      INTEGER, INTENT(IN)                                   :: nindex

      INTEGER                                               :: iproc, nprocs

      nprocs = SIZE(ndata)
      IF (ALLOCATED(buffer%msg_offset)) THEN
         IF (SIZE(buffer%msg_offset) /= nprocs + 1) DEALLOCATE (buffer%msg_offset, buffer%indx_offset)
      END IF
      IF (.NOT. ALLOCATED(buffer%msg_offset)) ALLOCATE (buffer%msg_offset(0:nprocs), buffer%indx_offset(0:nprocs))
      buffer%msg_offset(0) = 0
      buffer%indx_offset(0) = 0
      DO iproc = 0, nprocs - 1
         buffer%msg_offset(iproc + 1) = buffer%msg_offset(iproc) + ndata(iproc)
         buffer%indx_offset(iproc + 1) = buffer%indx_offset(iproc) + nblks(iproc)
      END DO

      IF (ALLOCATED(buffer%msg)) THEN
         IF (SIZE(buffer%msg) < buffer%msg_offset(nprocs)) DEALLOCATE (buffer%msg)
      END IF
      IF (.NOT. ALLOCATED(buffer%msg)) THEN
         ALLOCATE (buffer%msg(MAX(buffer%msg_offset(nprocs) + buffer%msg_offset(nprocs)/4, 1)))
      END IF
      ! records are exchanged as whole columns: the leading dimension must match
      IF (ALLOCATED(buffer%indx)) THEN
         IF (SIZE(buffer%indx, 1) /= nindex .OR. SIZE(buffer%indx, 2) < buffer%indx_offset(nprocs)) THEN
            DEALLOCATE (buffer%indx)
         END IF
      END IF
      IF (.NOT. ALLOCATED(buffer%indx)) THEN
         ALLOCATE (buffer%indx(nindex, MAX(buffer%indx_offset(nprocs) + buffer%indx_offset(nprocs)/4, 1)))
      END IF
   END SUBROUTINE grow_buffer

! **************************************************************************************************
!> \brief Exchange the contents of two buffers without copying.
!> \author Hans Pabst
! **************************************************************************************************
   SUBROUTINE swap_buffers(a, b)
      TYPE(matrix_buffer_type), INTENT(INOUT)               :: a, b

      TYPE(matrix_buffer_type)                              :: t

      CALL MOVE_ALLOC(a%msg, t%msg)
      CALL MOVE_ALLOC(b%msg, a%msg)
      CALL MOVE_ALLOC(t%msg, b%msg)
      CALL MOVE_ALLOC(a%indx, t%indx)
      CALL MOVE_ALLOC(b%indx, a%indx)
      CALL MOVE_ALLOC(t%indx, b%indx)
      CALL MOVE_ALLOC(a%msg_offset, t%msg_offset)
      CALL MOVE_ALLOC(b%msg_offset, a%msg_offset)
      CALL MOVE_ALLOC(t%msg_offset, b%msg_offset)
      CALL MOVE_ALLOC(a%indx_offset, t%indx_offset)
      CALL MOVE_ALLOC(b%indx_offset, a%indx_offset)
      CALL MOVE_ALLOC(t%indx_offset, b%indx_offset)
   END SUBROUTINE swap_buffers

! **************************************************************************************************
!> \brief Deallocate the arrays of a buffer.
!> \author Hans Pabst
! **************************************************************************************************
   SUBROUTINE release_buffer(buffer)
      TYPE(matrix_buffer_type), INTENT(INOUT)               :: buffer

      IF (ALLOCATED(buffer%msg)) DEALLOCATE (buffer%msg)
      IF (ALLOCATED(buffer%indx)) DEALLOCATE (buffer%indx)
      IF (ALLOCATED(buffer%msg_offset)) DEALLOCATE (buffer%msg_offset)
      IF (ALLOCATED(buffer%indx_offset)) DEALLOCATE (buffer%indx_offset)
   END SUBROUTINE release_buffer

END MODULE dbt_reshape_ops
