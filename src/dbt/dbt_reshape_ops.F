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
      block_nd, create_block, destroy_block, dbt_iterator_type, dbt_iterator_next_block, &
      dbt_iterator_blocks_left, dbt_iterator_start, dbt_iterator_stop, dbt_get_block, &
      dbt_reserve_blocks, dbt_put_block
   USE dbt_types, ONLY: dbt_blk_sizes, &
                        dbt_type, &
                        ndims_tensor, &
                        dbt_get_stored_coordinates, &
                        dbt_clear
   USE matrix_communication, ONLY: allocate_matrix_buffers, exchange_matrix_buffers, matrix_buffer_type
   USE message_passing, ONLY: mp_comm_type
   USE util, ONLY: get_limit

#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_reshape_ops'

   PUBLIC :: dbt_reshape

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

      INTEGER                                            :: iproc, numnodes, &
                                                            handle, iblk, offset, ndata
      INTEGER, DIMENSION(2)                              :: block_range
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(block_nd)                                     :: blk_data
      TYPE(matrix_buffer_type)                          :: buffer_recv, buffer_send
      INTEGER, DIMENSION(ndims_tensor(tensor_in))        :: blk_size, ind_nd
      LOGICAL :: found, summation_prv, move_prv

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

      CPASSERT(tensor_out%valid)

      IF (.NOT. summation_prv) CALL dbt_clear(tensor_out)

      mp_comm = tensor_in%pgrid%mp_comm_2d
      numnodes = mp_comm%num_pe
      ALLOCATE (nblks_send_total(0:numnodes - 1), ndata_send_total(0:numnodes - 1), source=0)
      ALLOCATE (nblks_recv_total(0:numnodes - 1), ndata_recv_total(0:numnodes - 1), source=0)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor_in,tensor_out,summation) &
!$OMP SHARED(buffer_send,buffer_recv,mp_comm,numnodes) &
!$OMP SHARED(nblks_send_total,ndata_send_total,nblks_recv_total,ndata_recv_total) &
!$OMP PRIVATE(nblks_send_mythread,ndata_send_mythread) &
!$OMP PRIVATE(iter,ind_nd,blk_size,blk_data,found,iproc) &
!$OMP PRIVATE(block_range,offset,ndata,iblk)
      ALLOCATE (nblks_send_mythread(0:numnodes - 1), ndata_send_mythread(0:numnodes - 1), source=0)

      CALL dbt_iterator_start(iter, tensor_in)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind_nd, blk_size=blk_size)
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
      CALL allocate_matrix_buffers(ndata_send_total, ndata_recv_total, ndims_tensor(tensor_in) + 1, &
                                   buffer_send, buffer_recv, nblks_send_total, nblks_recv_total)
!$OMP END MASTER
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
         offset = buffer_send%msg_offset(iproc) + ndata_send_mythread(iproc)
         buffer_send%msg(offset + 1:offset + ndata) = blk_data%blk(:)
         ! insert block index
         nblks_send_mythread(iproc) = nblks_send_mythread(iproc) - 1
         iblk = buffer_send%indx_offset(iproc) + nblks_send_mythread(iproc) + 1
         buffer_send%indx(2:, iblk) = ind_nd(:)
         ! Payload offsets in the metadata are relative to the peer's segment.
         buffer_send%indx(1, iblk) = ndata_send_mythread(iproc)
         CALL destroy_block(blk_data)
      END DO
      CALL dbt_iterator_stop(iter)
      DEALLOCATE (nblks_send_mythread, ndata_send_mythread)
!$OMP BARRIER

!$OMP MASTER
      IF (numnodes > 1) THEN
         CALL exchange_matrix_buffers(mp_comm, buffer_recv, buffer_send)
         DEALLOCATE (buffer_send%msg, buffer_send%indx)
      ELSE
         ! The single-rank receive buffer takes ownership of the packed arrays.
         CALL MOVE_ALLOC(buffer_send%msg, buffer_recv%msg)
         CALL MOVE_ALLOC(buffer_send%indx, buffer_recv%indx)
      END IF
      DEALLOCATE (buffer_send%msg_offset, buffer_send%indx_offset)
!$OMP END MASTER
!$OMP BARRIER

      block_range = get_limit(SIZE(buffer_recv%indx, 2), omp_get_num_threads(), omp_get_thread_num())
      CALL dbt_reserve_blocks(tensor_out, TRANSPOSE(buffer_recv%indx(2:, block_range(1):block_range(2))))

      DO iproc = 0, numnodes - 1
!$OMP DO
         DO iblk = buffer_recv%indx_offset(iproc) + 1, buffer_recv%indx_offset(iproc + 1)
            ind_nd(:) = buffer_recv%indx(2:, iblk)
            CALL dbt_blk_sizes(tensor_out, ind_nd, blk_size)
            offset = buffer_recv%msg_offset(iproc) + buffer_recv%indx(1, iblk)
            ndata = PRODUCT(blk_size)
            CALL create_block(blk_data, blk_size, &
                              array=buffer_recv%msg(offset + 1:offset + ndata))
            CALL dbt_put_block(tensor_out, ind_nd, blk_data, summation=summation)
            CALL destroy_block(blk_data)
         END DO
!$OMP END DO
      END DO

!$OMP MASTER
      DEALLOCATE (buffer_recv%msg, buffer_recv%indx, buffer_recv%msg_offset, buffer_recv%indx_offset)
!$OMP END MASTER
!$OMP END PARALLEL

      DEALLOCATE (nblks_recv_total, ndata_recv_total)
      DEALLOCATE (nblks_send_total, ndata_send_total)

      IF (move_prv) CALL dbt_clear(tensor_in)

      CALL timestop(handle)
   END SUBROUTINE dbt_reshape

END MODULE dbt_reshape_ops
