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

   USE OMP_LIB, ONLY: omp_get_num_threads
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
      INTEGER, DIMENSION(:, :), ALLOCATABLE      :: blocks
      REAL(dp), DIMENSION(:), ALLOCATABLE        :: data
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief copy data (involves reshape)
!>        tensor_out = tensor_out + tensor_in move_data memory optimization:
!>        transfer data from tensor_in to tensor_out s.t. tensor_in is empty on return
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbt_reshape(tensor_in, tensor_out, summation, move_data)

      TYPE(dbt_type), INTENT(INOUT)               :: tensor_in, tensor_out
      LOGICAL, INTENT(IN), OPTIONAL                    :: summation
      LOGICAL, INTENT(IN), OPTIONAL                    :: move_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_reshape'

      INTEGER                                            :: iproc, mp_comm, mynode, &
                                                            numnodes, bcount, handle, &
                                                            max_nblk_per_thread, iblk, offset, ndata
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: blks_to_allocate
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(block_nd)                                     :: blk_data
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(ndims_tensor(tensor_in))        :: blk_size, ind_nd
      LOGICAL :: found, summation_prv, move_prv

      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nblks_send_total, ndata_send_total, &
                                                            nblks_recv_total, ndata_recv_total, &
                                                            nblks_send_mythread, ndata_send_mythread, &
                                                            nblks_recv_mythread, ndata_recv_mythread

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

      CPASSERT(tensor_out%valid)

      IF (.NOT. summation_prv) CALL dbt_clear(tensor_out)

      mp_comm = tensor_in%pgrid%mp_comm_2d
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (buffer_send(0:numnodes - 1), buffer_recv(0:numnodes - 1))

      ALLOCATE (nblks_send_total(0:numnodes - 1), ndata_send_total(0:numnodes - 1))
      ALLOCATE (nblks_recv_total(0:numnodes - 1), ndata_recv_total(0:numnodes - 1))
      nblks_send_total(:) = 0
      ndata_send_total(:) = 0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor_in,tensor_out,summation) &
!$OMP SHARED(buffer_send,buffer_recv,mp_comm,numnodes) &
!$OMP SHARED(nblks_send_total,ndata_send_total,nblks_recv_total,ndata_recv_total) &
!$OMP PRIVATE(nblks_send_mythread,ndata_send_mythread,nblks_recv_mythread,ndata_recv_mythread) &
!$OMP PRIVATE(iter,ind_nd,blk_size,blk_data,found,iproc) &
!$OMP PRIVATE(blks_to_allocate,max_nblk_per_thread,bcount,offset,ndata,iblk)
      ALLOCATE (nblks_send_mythread(0:numnodes - 1), ndata_send_mythread(0:numnodes - 1))
      ALLOCATE (nblks_recv_mythread(0:numnodes - 1), ndata_recv_mythread(0:numnodes - 1))
      nblks_send_mythread(:) = 0
      ndata_send_mythread(:) = 0

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size)
         CALL dbt_get_stored_coordinates(tensor_out, ind_nd, iproc)
         nblks_send_mythread(iproc) = nblks_send_mythread(iproc) + 1
         ndata_send_mythread(iproc) = ndata_send_mythread(iproc) + PRODUCT(blk_size)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP CRITICAL
      nblks_send_total(:) = nblks_send_total(:) + nblks_send_mythread(:)
      ndata_send_total(:) = ndata_send_total(:) + ndata_send_mythread(:)
      nblks_send_mythread(:) = nblks_send_total(:) ! current totals indicate slot for this thread
      ndata_send_mythread(:) = ndata_send_total(:)
!$OMP END CRITICAL
!$OMP BARRIER

!$OMP MASTER
      CALL mp_alltoall(nblks_send_total, nblks_recv_total, 1, mp_comm)
      CALL mp_alltoall(ndata_send_total, ndata_recv_total, 1, mp_comm)
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO
      DO iproc = 0, numnodes - 1
         ALLOCATE (buffer_send(iproc)%data(ndata_send_total(iproc)))
         ALLOCATE (buffer_recv(iproc)%data(ndata_recv_total(iproc)))
         ! going to use buffer%blocks(:,0) to store data offsets
         ALLOCATE (buffer_send(iproc)%blocks(nblks_send_total(iproc), 0:ndims_tensor(tensor_in)))
         ALLOCATE (buffer_recv(iproc)%blocks(nblks_recv_total(iproc), 0:ndims_tensor(tensor_in)))
      END DO
!$OMP END DO
!$OMP BARRIER

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size)
         CALL dbt_get_stored_coordinates(tensor_out, ind_nd, iproc)
         CALL dbt_get_block(tensor_in, ind_nd, blk_data, found)
         CPASSERT(found)
         ! insert block data
         ndata = PRODUCT(blk_size)
         ndata_send_mythread(iproc) = ndata_send_mythread(iproc) - ndata
         offset = ndata_send_mythread(iproc)
         buffer_send(iproc)%data(offset + 1:offset + ndata) = blk_data%blk(:)
         ! insert block index
         nblks_send_mythread(iproc) = nblks_send_mythread(iproc) - 1
         iblk = nblks_send_mythread(iproc) + 1
         buffer_send(iproc)%blocks(iblk, 1:) = ind_nd(:)
         buffer_send(iproc)%blocks(iblk, 0) = offset
         CALL destroy_block(blk_data)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP BARRIER

      CALL dbt_communicate_buffer(mp_comm, buffer_recv, buffer_send)
!$OMP BARRIER

!$OMP DO
      DO iproc = 0, numnodes - 1
         DEALLOCATE (buffer_send(iproc)%blocks, buffer_send(iproc)%data)
      END DO
!$OMP END DO

      max_nblk_per_thread = SUM(nblks_recv_total)/omp_get_num_threads() + 1
      ALLOCATE (blks_to_allocate(max_nblk_per_thread, ndims_tensor(tensor_in)))
      bcount = 0
      DO iproc = 0, numnodes - 1
!$OMP DO
         DO iblk = 1, nblks_recv_total(iproc)
            bcount = bcount + 1
            blks_to_allocate(bcount, :) = buffer_recv(iproc)%blocks(iblk, 1:)
         END DO
!$OMP END DO
      END DO
      CALL dbt_reserve_blocks(tensor_out, blks_to_allocate(1:bcount, :))
      DEALLOCATE (blks_to_allocate)
!$OMP BARRIER

      DO iproc = 0, numnodes - 1
!$OMP DO
         DO iblk = 1, nblks_recv_total(iproc)
            ind_nd(:) = buffer_recv(iproc)%blocks(iblk, 1:)
            CALL dbt_blk_sizes(tensor_out, ind_nd, blk_size)
            offset = buffer_recv(iproc)%blocks(iblk, 0)
            ndata = PRODUCT(blk_size)
            CALL create_block(blk_data, blk_size, &
                              array=buffer_recv(iproc)%data(offset + 1:offset + ndata))
            CALL dbt_put_block(tensor_out, ind_nd, blk_data, summation=summation)
            CALL destroy_block(blk_data)
         END DO
!$OMP END DO
      END DO
!$OMP END PARALLEL

      DO iproc = 0, numnodes - 1
         DEALLOCATE (buffer_recv(iproc)%blocks, buffer_recv(iproc)%data)
      END DO

      IF (move_prv) CALL dbt_clear(tensor_in)

      CALL timestop(handle)
   END SUBROUTINE dbt_reshape

! **************************************************************************************************
!> \brief communicate buffer
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE dbt_communicate_buffer(mp_comm, buffer_recv, buffer_send)
      INTEGER, INTENT(IN)                    :: mp_comm
      TYPE(block_buffer_type), DIMENSION(0:), INTENT(INOUT) :: buffer_recv, buffer_send

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_communicate_buffer'

      INTEGER                                               :: iproc, mynode, numnodes, &
                                                               rec_counter, send_counter, i
      INTEGER, ALLOCATABLE, DIMENSION(:, :)                 :: req_array
      INTEGER                                               :: handle

      CALL timeset(routineN, handle)
      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (numnodes > 1) THEN
!$OMP MASTER
         send_counter = 0
         rec_counter = 0

         ALLOCATE (req_array(1:numnodes, 4))

         DO iproc = 0, numnodes - 1
            IF (SIZE(buffer_recv(iproc)%blocks) > 0) THEN
               rec_counter = rec_counter + 1
               CALL mp_irecv(buffer_recv(iproc)%blocks, iproc, mp_comm, req_array(rec_counter, 3), tag=4)
               CALL mp_irecv(buffer_recv(iproc)%data, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO iproc = 0, numnodes - 1
            IF (SIZE(buffer_send(iproc)%blocks) > 0) THEN
               send_counter = send_counter + 1
               CALL mp_isend(buffer_send(iproc)%blocks, iproc, mp_comm, req_array(send_counter, 1), tag=4)
               CALL mp_isend(buffer_send(iproc)%data, iproc, mp_comm, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         IF (send_counter > 0) THEN
            CALL mp_waitall(req_array(1:send_counter, 1:2))
         END IF
         IF (rec_counter > 0) THEN
            CALL mp_waitall(req_array(1:rec_counter, 3:4))
         END IF
!$OMP END MASTER

      ELSE
!$OMP DO SCHEDULE(static, 512)
         DO i = 1, SIZE(buffer_send(0)%blocks, 1)
            buffer_recv(0)%blocks(i, :) = buffer_send(0)%blocks(i, :)
         END DO
!$OMP END DO
!$OMP DO SCHEDULE(static, 512)
         DO i = 1, SIZE(buffer_send(0)%data)
            buffer_recv(0)%data(i) = buffer_send(0)%data(i)
         END DO
!$OMP END DO
      END IF
      CALL timestop(handle)

   END SUBROUTINE dbt_communicate_buffer

END MODULE dbt_reshape_ops
