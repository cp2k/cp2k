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
      dbt_iterator_stop, dbt_get_block_data, dbt_put_block_data, dbt_reserve_blocks, dbt_copy_block
   USE dbt_tas_mm, ONLY: dbt_tas_batched_mm_complete
   USE dbt_types, ONLY: dbt_blk_sizes, &
                        dbt_type, &
                        dbt_get_info, &
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

   PUBLIC :: dbt_reshape, dbt_copy_indexed, dbt_reshape_release

   ! Non-owning references allow scalar copies and image gathers to share one transport.
   TYPE tensor_source_type
      TYPE(dbt_type), POINTER :: tensor => NULL()
      INTEGER :: stack_offset = 0
   END TYPE tensor_source_type

   INTERFACE dbt_copy_indexed
      MODULE PROCEDURE dbt_reshape, dbt_gather
   END INTERFACE dbt_copy_indexed

   ! Communication buffers persist across calls such that pages are faulted in once rather than per
   ! call; per block, the index holds the data offset (relative to the peer's data) and nd-index.
   TYPE(matrix_buffer_type), SAVE :: buffer_send, buffer_recv

CONTAINS

! **************************************************************************************************
!> \brief copy data (involves reshape)
!>        tensor_out = tensor_out + tensor_in move_data memory optimization:
!>        transfer data from tensor_in to tensor_out s.t. tensor_in is empty on return
!> \param bounds crop tensor data while copying (as dbt_crop): start and end index per dimension
!> \param tensor_in source tensor on the parent communicator
!> \param tensor_out destination tensor; smaller communicators denote equal, contiguous subgroups
!> \param summation add to destination instead of replacing it
!> \param move_data clear the source after copying
!> \param index_map optional source-to-destination block indices; zero skips a block
!> \param group_mask optional subgroup selection by source block index
!> \param group_dim source dimension indexing group_mask
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbt_reshape(tensor_in, tensor_out, summation, move_data, bounds, index_map, group_mask, group_dim)

      TYPE(dbt_type), TARGET, INTENT(INOUT) :: tensor_in
      TYPE(dbt_type), INTENT(INOUT) :: tensor_out
      LOGICAL, INTENT(IN), OPTIONAL :: summation, move_data
      INTEGER, INTENT(IN), OPTIONAL :: bounds(:, :), index_map(:, :), group_dim
      LOGICAL, INTENT(IN), OPTIONAL :: group_mask(:, :)

      TYPE(tensor_source_type) :: sources(1)

      sources(1)%tensor => tensor_in
      CALL dbt_reshape_sources(sources, tensor_out, tensor_in%pgrid%mp_comm_2d, summation, move_data, bounds, &
                               index_map=index_map, group_mask=group_mask, group_dim=group_dim)
   END SUBROUTINE dbt_reshape

! **************************************************************************************************
!> \brief Gather whole blocks from selected tensors, optionally relabelling each block index.
!>        Zero images or mapped indices are skipped. Blocks retain their sizes. Each image occupies
!>        one source tensor extent along stack_dim, including skipped images. Retained index maps
!>        must be injective; the stack-dimension map stays within one source extent.
!> \param tensor_in source tensors on one communicator, with identical shapes
!> \param tensor_out destination tensor
!> \param images source tensor indices; zero skips an image
!> \param index_map source-to-destination block indices, padded to the largest source dimension
!> \param stack_dim dimension along which images are stacked
! **************************************************************************************************
   SUBROUTINE dbt_gather(tensor_in, tensor_out, images, index_map, stack_dim)
      TYPE(dbt_type), TARGET, INTENT(INOUT) :: tensor_in(:)
      TYPE(dbt_type), INTENT(INOUT) :: tensor_out
      INTEGER, INTENT(IN) :: images(:), index_map(:, :), stack_dim

      INTEGER :: i, nblks(SIZE(index_map, 2))
      TYPE(tensor_source_type) :: sources(SIZE(images))

      CPASSERT(SIZE(tensor_in) > 0)
      CPASSERT(SIZE(index_map, 2) == ndims_tensor(tensor_out))
      CPASSERT(stack_dim >= 1 .AND. stack_dim <= SIZE(nblks))
      CALL dbt_get_info(tensor_in(1), nblks_total=nblks)
      CPASSERT(SIZE(index_map, 1) >= MAXVAL(nblks))
      CPASSERT(ALL(index_map(:nblks(stack_dim), stack_dim) <= nblks(stack_dim)))
      DO i = 1, SIZE(images)
         CPASSERT(images(i) >= 0 .AND. images(i) <= SIZE(tensor_in))
         IF (images(i) == 0) CYCLE
         CALL dbt_tas_batched_mm_complete(tensor_in(images(i))%matrix_rep, warn=.TRUE.)
         sources(i)%tensor => tensor_in(images(i))
         sources(i)%stack_offset = (i - 1)*nblks(stack_dim)
      END DO
      ! As in dbt_copy, invalidate any input replicas kept by the previous contraction.
      CALL dbt_tas_batched_mm_complete(tensor_out%matrix_rep, warn=.TRUE.)
      CALL dbt_reshape_sources(sources, tensor_out, tensor_in(1)%pgrid%mp_comm_2d, &
                               index_map=index_map, stack_dim=stack_dim)
   END SUBROUTINE dbt_gather

! **************************************************************************************************
!> \brief Shared block transport for copies, indexed gathers and subgroup replication.
!>        Subgroups have identical layouts and contiguous rank ranges in comm. group_mask selects
!>        replicas by a source block index; with no mask, every subgroup receives the block.
! **************************************************************************************************
   SUBROUTINE dbt_reshape_sources(sources, tensor_out, comm, summation, move_data, bounds, &
                                  index_map, stack_dim, group_mask, group_dim)

      TYPE(tensor_source_type), INTENT(IN) :: sources(:)
      TYPE(dbt_type), INTENT(INOUT) :: tensor_out
      CLASS(mp_comm_type), INTENT(IN) :: comm
      LOGICAL, INTENT(IN), OPTIONAL                    :: summation
      LOGICAL, INTENT(IN), OPTIONAL                    :: move_data
      INTEGER, INTENT(IN), OPTIONAL :: bounds(:, :), index_map(:, :), stack_dim, group_dim
      LOGICAL, INTENT(IN), OPTIONAL :: group_mask(:, :)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reshape'

      INTEGER                                            :: iproc, numnodes, nindex, &
                                                            handle, handle2, iblk, offset, ndata, &
                                                            isource, d, igroup, ngroups, group_size, peer, first_offset
      INTEGER, DIMENSION(2)                              :: block_range
      TYPE(dbt_iterator_type)                            :: iter
      CHARACTER(LEN=32)                                  :: name_pack, name_unpack
      INTEGER, DIMENSION(ndims_tensor(tensor_out)) :: blk_size, blk_offset, ind_nd, dest, dest_size
      INTEGER, DIMENSION(2, ndims_tensor(tensor_out)) :: bounds_prv
      LOGICAL :: found, summation_prv, move_prv, do_bounds, overlap, inside
      TYPE(dbt_type), POINTER :: tensor_in
      INTEGER :: source_shape(ndims_tensor(tensor_out))

      INTEGER, ALLOCATABLE :: counts_send(:, :), counts_recv(:, :)
      INTEGER, ALLOCATABLE :: nblks_send_mythread(:), ndata_send_mythread(:)

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
      group_size = tensor_out%pgrid%mp_comm_2d%num_pe
      CPASSERT(MOD(comm%num_pe, group_size) == 0)
      CPASSERT(MOD(comm%mepos, group_size) == tensor_out%pgrid%mp_comm_2d%mepos)
      ngroups = comm%num_pe/group_size
      CPASSERT(PRESENT(group_mask) .EQV. PRESENT(group_dim))
      IF (PRESENT(group_mask)) THEN
         CPASSERT(SIZE(group_mask, 1) == ngroups)
         CPASSERT(group_dim >= 1 .AND. group_dim <= SIZE(ind_nd))
      END IF
      IF (PRESENT(index_map)) THEN
         CPASSERT(SIZE(index_map, 2) == SIZE(ind_nd))
         CPASSERT(ALL(index_map >= 0))
      END IF

      IF (.NOT. summation_prv) THEN
         CALL timeset(routineN//"_clear", handle2)
         CALL dbt_clear(tensor_out)
         CALL timestop(handle2)
      END IF

      numnodes = comm%num_pe
      nindex = ndims_tensor(tensor_out) + 1
      ! separate timers for blocks that need a reshape (index order differs from matrix layout)
      name_pack = routineN//"_pack"
      DO isource = 1, SIZE(sources)
         IF (.NOT. ASSOCIATED(sources(isource)%tensor)) CYCLE
         CALL dbt_get_info(sources(isource)%tensor, nblks_total=source_shape)
         IF (PRESENT(index_map)) THEN
            CPASSERT(SIZE(index_map, 1) >= MAXVAL(source_shape))
         END IF
         IF (PRESENT(group_mask)) THEN
            CPASSERT(SIZE(group_mask, 2) == source_shape(group_dim))
         END IF
         IF (permuted(sources(isource)%tensor)) name_pack = routineN//"_pack_perm"
      END DO
      name_unpack = routineN//"_unpack"
      IF (permuted(tensor_out)) name_unpack = routineN//"_unpack_perm"
      IF (numnodes == 1) THEN
         CALL dbt_reshape_local(sources, tensor_out, summation_prv, bounds, index_map, stack_dim, group_mask, group_dim)
      ELSE
         ALLOCATE (counts_send(2, 0:numnodes - 1), counts_recv(2, 0:numnodes - 1), source=0)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(sources,tensor_out,summation,comm,numnodes,nindex,name_pack,name_unpack) &
!$OMP SHARED(index_map,stack_dim,group_mask,group_dim,group_size,ngroups) &
!$OMP SHARED(do_bounds,bounds_prv,buffer_send,buffer_recv) &
!$OMP SHARED(counts_send,counts_recv) &
!$OMP PRIVATE(nblks_send_mythread,ndata_send_mythread) &
!$OMP PRIVATE(iter,ind_nd,blk_size,blk_offset,found,overlap,inside,iproc,handle2) &
!$OMP PRIVATE(block_range,offset,ndata,iblk,tensor_in,isource,d,dest,dest_size,igroup,peer,first_offset)
         CALL timeset(routineN//"_count", handle2)
         ALLOCATE (nblks_send_mythread(0:numnodes - 1), ndata_send_mythread(0:numnodes - 1), source=0)

         DO isource = 1, SIZE(sources)
            tensor_in => sources(isource)%tensor
            IF (.NOT. ASSOCIATED(tensor_in)) CYCLE
            CALL dbt_iterator_start(iter, tensor_in)
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size, blk_offset=blk_offset)
               IF (do_bounds) THEN
                  CALL block_bounds(bounds_prv, blk_offset, blk_size, overlap, inside)
                  IF (.NOT. overlap) CYCLE
               END IF
               dest = ind_nd
               IF (PRESENT(index_map)) THEN
                  DO d = 1, SIZE(dest)
                     dest(d) = index_map(ind_nd(d), d)
                  END DO
                  IF (ANY(dest == 0)) CYCLE
               END IF
               IF (PRESENT(stack_dim)) dest(stack_dim) = dest(stack_dim) + sources(isource)%stack_offset
               IF (PRESENT(index_map)) THEN
                  CALL dbt_blk_sizes(tensor_out, dest, dest_size)
                  CPASSERT(ALL(dest_size == blk_size))
               END IF
               CALL dbt_get_stored_coordinates(tensor_out, dest, iproc)
               DO igroup = 1, ngroups
                  IF (PRESENT(group_mask)) THEN
                     IF (.NOT. group_mask(igroup, ind_nd(group_dim))) CYCLE
                  END IF
                  peer = (igroup - 1)*group_size + iproc
                  nblks_send_mythread(peer) = nblks_send_mythread(peer) + 1
                  ndata_send_mythread(peer) = ndata_send_mythread(peer) + PRODUCT(blk_size)
               END DO
            END DO
            CALL dbt_iterator_stop(iter)
         END DO
!$OMP CRITICAL(omp_dbt_reshape)
         counts_send(1, :) = counts_send(1, :) + ndata_send_mythread(:)
         counts_send(2, :) = counts_send(2, :) + nblks_send_mythread(:)
         ndata_send_mythread(:) = counts_send(1, :) ! current totals indicate slot for this thread
         nblks_send_mythread(:) = counts_send(2, :)
!$OMP END CRITICAL(omp_dbt_reshape)
!$OMP BARRIER

!$OMP MASTER
         CALL comm%alltoall(counts_send, counts_recv, 2)
!$OMP END MASTER
!$OMP BARRIER
         CALL timestop(handle2)

         CALL timeset(TRIM(name_pack), handle2)
!$OMP MASTER
         CALL grow_buffer(buffer_send, counts_send(1, :), counts_send(2, :), nindex)
         CALL grow_buffer(buffer_recv, counts_recv(1, :), counts_recv(2, :), nindex)
!$OMP END MASTER
!$OMP BARRIER

         DO isource = 1, SIZE(sources)
            tensor_in => sources(isource)%tensor
            IF (.NOT. ASSOCIATED(tensor_in)) CYCLE
            CALL dbt_iterator_start(iter, tensor_in)
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size, blk_offset=blk_offset)
               inside = .TRUE.
               IF (do_bounds) THEN
                  CALL block_bounds(bounds_prv, blk_offset, blk_size, overlap, inside)
                  IF (.NOT. overlap) CYCLE
               END IF
               dest = ind_nd
               IF (PRESENT(index_map)) THEN
                  DO d = 1, SIZE(dest)
                     dest(d) = index_map(ind_nd(d), d)
                  END DO
                  IF (ANY(dest == 0)) CYCLE
               END IF
               IF (PRESENT(stack_dim)) dest(stack_dim) = dest(stack_dim) + sources(isource)%stack_offset
               CALL dbt_get_stored_coordinates(tensor_out, dest, iproc)
               first_offset = -1
               DO igroup = 1, ngroups
                  IF (PRESENT(group_mask)) THEN
                     IF (.NOT. group_mask(igroup, ind_nd(group_dim))) CYCLE
                  END IF
                  peer = (igroup - 1)*group_size + iproc
                  ! Insert block data, reading the source only once for all replicas.
                  ndata = PRODUCT(blk_size)
                  ndata_send_mythread(peer) = ndata_send_mythread(peer) - ndata
                  offset = buffer_send%msg_offset(peer) + ndata_send_mythread(peer)
                  IF (first_offset < 0) THEN
                     CALL dbt_get_block_data(tensor_in, ind_nd, blk_size, buffer_send%msg(offset + 1:offset + ndata), found)
                     CPASSERT(found)
                     IF (.NOT. inside) THEN
                        CALL block_crop(bounds_prv, blk_offset, blk_size, buffer_send%msg(offset + 1:offset + ndata))
                     END IF
                     first_offset = offset
                  ELSE
                     buffer_send%msg(offset + 1:offset + ndata) = buffer_send%msg(first_offset + 1:first_offset + ndata)
                  END IF
                  ! Payload offsets in the metadata are relative to the peer's segment.
                  nblks_send_mythread(peer) = nblks_send_mythread(peer) - 1
                  iblk = buffer_send%indx_offset(peer) + nblks_send_mythread(peer) + 1
                  buffer_send%indx(2:, iblk) = dest(:)
                  buffer_send%indx(1, iblk) = ndata_send_mythread(peer)
               END DO
            END DO
            CALL dbt_iterator_stop(iter)
         END DO
         DEALLOCATE (nblks_send_mythread, ndata_send_mythread)
!$OMP BARRIER
         CALL timestop(handle2)

!$OMP MASTER
         CALL exchange_matrix_buffers(comm, buffer_recv, buffer_send)
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

         DEALLOCATE (counts_send, counts_recv)
      END IF

      IF (move_prv) THEN
         CALL timeset(routineN//"_clear", handle2)
         DO isource = 1, SIZE(sources)
            IF (ASSOCIATED(sources(isource)%tensor)) CALL dbt_clear(sources(isource)%tensor)
         END DO
         CALL timestop(handle2)
      END IF

      CALL timestop(handle)
   END SUBROUTINE dbt_reshape_sources

! **************************************************************************************************
!> \brief Copy single-rank tensors without communication buffers or a separate reserve pass.
!>        Matching block layouts are copied directly; other layouts and partial crops use one
!>        reusable buffer per thread.
! **************************************************************************************************
   SUBROUTINE dbt_reshape_local(sources, tensor_out, summation, bounds, index_map, stack_dim, group_mask, group_dim)
      TYPE(tensor_source_type), INTENT(IN) :: sources(:)
      TYPE(dbt_type), INTENT(INOUT) :: tensor_out
      LOGICAL, INTENT(IN) :: summation
      INTEGER, INTENT(IN), OPTIONAL :: bounds(:, :), index_map(:, :), stack_dim, group_dim
      LOGICAL, INTENT(IN), OPTIONAL :: group_mask(:, :)

      INTEGER :: isource, d, ndata
      INTEGER :: ind(ndims_tensor(tensor_out)), dest(ndims_tensor(tensor_out)), &
                 sizes(ndims_tensor(tensor_out)), offsets(ndims_tensor(tensor_out)), &
                 dest_sizes(ndims_tensor(tensor_out))
      LOGICAL :: direct_copy(SIZE(sources)), found, overlap, inside
      REAL(dp), ALLOCATABLE :: data(:)
      TYPE(dbt_type), POINTER :: tensor_in
      TYPE(dbt_iterator_type) :: iter

      DO isource = 1, SIZE(sources)
         tensor_in => sources(isource)%tensor
         IF (.NOT. ASSOCIATED(tensor_in)) CYCLE
         direct_copy(isource) = &
            SIZE(tensor_in%nd_index_blk%map1_2d) == SIZE(tensor_out%nd_index_blk%map1_2d) .AND. &
            ALL(tensor_in%nd_index_blk%map_nd == tensor_out%nd_index_blk%map_nd)
      END DO

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(sources,tensor_out,summation,bounds,index_map,stack_dim,group_mask,group_dim,direct_copy) &
!$OMP PRIVATE(iter,tensor_in,isource,ind,dest,sizes,offsets,dest_sizes,d,ndata,data,found,overlap,inside)
      DO isource = 1, SIZE(sources)
         tensor_in => sources(isource)%tensor
         IF (.NOT. ASSOCIATED(tensor_in)) CYCLE
         CALL dbt_iterator_start(iter, tensor_in)
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind, blk_size=sizes, blk_offset=offsets)
            IF (PRESENT(group_mask)) THEN
               IF (.NOT. group_mask(1, ind(group_dim))) CYCLE
            END IF
            inside = .TRUE.
            IF (PRESENT(bounds)) THEN
               CALL block_bounds(bounds, offsets, sizes, overlap, inside)
               IF (.NOT. overlap) CYCLE
            END IF
            dest(:) = ind(:)
            IF (PRESENT(index_map)) THEN
               DO d = 1, SIZE(dest)
                  dest(d) = index_map(ind(d), d)
               END DO
               IF (ANY(dest == 0)) CYCLE
            END IF
            IF (PRESENT(stack_dim)) dest(stack_dim) = dest(stack_dim) + sources(isource)%stack_offset
            IF (PRESENT(index_map)) THEN
               CALL dbt_blk_sizes(tensor_out, dest, dest_sizes)
               CPASSERT(ALL(dest_sizes == sizes))
            END IF
            IF (direct_copy(isource) .AND. inside) THEN
               CALL dbt_copy_block(tensor_in, tensor_out, ind, summation=summation, ind_out=dest)
            ELSE
               ndata = PRODUCT(sizes)
               IF (.NOT. ALLOCATED(data)) ALLOCATE (data(ndata))
               IF (SIZE(data) < ndata) THEN
                  DEALLOCATE (data)
                  ALLOCATE (data(ndata))
               END IF
               CALL dbt_get_block_data(tensor_in, ind, sizes, data(:ndata), found)
               CPASSERT(found)
               IF (.NOT. inside) CALL block_crop(bounds, offsets, sizes, data(:ndata))
               CALL dbt_put_block_data(tensor_out, dest, sizes, data(:ndata), summation=summation)
            END IF
         END DO
         CALL dbt_iterator_stop(iter)
      END DO
      IF (ALLOCATED(data)) DEALLOCATE (data)
!$OMP END PARALLEL
   END SUBROUTINE dbt_reshape_local

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
