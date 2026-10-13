!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exchange numerical matrix payloads and integer metadata between MPI ranks.
! **************************************************************************************************
MODULE matrix_communication
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_request_type,&
                                              mp_waitall
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'matrix_communication'
   PUBLIC :: matrix_buffer_type, allocate_matrix_buffers, exchange_matrix_buffers

   TYPE matrix_buffer_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: msg
      ! Metadata records are columns, contiguous within each rank's segment.
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: indx
      ! Zero-based segment offsets, including the final total at index nprocs.
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: msg_offset, indx_offset
   END TYPE matrix_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief Allocate flat, rank-grouped matrix buffers and their zero-based offsets.
!> \param send_count numerical payload length sent to each rank
!> \param recv_count numerical payload length received from each rank
!> \param n_index number of integer coordinates per metadata record
!> \param buffer_send output send buffer
!> \param buffer_rec output receive buffer
!> \param send_blocks optional metadata-record counts, distinct from payload counts
!> \param recv_blocks optional received metadata-record counts
! **************************************************************************************************
   SUBROUTINE allocate_matrix_buffers(send_count, recv_count, n_index, buffer_send, buffer_rec, &
                                      send_blocks, recv_blocks)
      INTEGER, DIMENSION(0:), INTENT(IN)                 :: send_count, recv_count
      INTEGER, INTENT(IN)                                :: n_index
      TYPE(matrix_buffer_type), INTENT(OUT)              :: buffer_send, buffer_rec
      INTEGER, DIMENSION(0:), INTENT(IN), OPTIONAL       :: send_blocks, recv_blocks

      INTEGER                                            :: nprocs, peer

      CPASSERT(PRESENT(send_blocks) .EQV. PRESENT(recv_blocks))
      nprocs = SIZE(send_count)
      ALLOCATE (buffer_send%msg_offset(0:nprocs), buffer_rec%msg_offset(0:nprocs), &
                buffer_send%indx_offset(0:nprocs), buffer_rec%indx_offset(0:nprocs))
      buffer_send%msg_offset(0) = 0
      buffer_rec%msg_offset(0) = 0
      buffer_send%indx_offset(0) = 0
      buffer_rec%indx_offset(0) = 0
      DO peer = 0, nprocs - 1
         buffer_send%msg_offset(peer + 1) = buffer_send%msg_offset(peer) + send_count(peer)
         buffer_rec%msg_offset(peer + 1) = buffer_rec%msg_offset(peer) + recv_count(peer)
         buffer_send%indx_offset(peer + 1) = buffer_send%indx_offset(peer) + send_count(peer)
         buffer_rec%indx_offset(peer + 1) = buffer_rec%indx_offset(peer) + recv_count(peer)
         IF (PRESENT(send_blocks)) THEN
            buffer_send%indx_offset(peer + 1) = buffer_send%indx_offset(peer) + send_blocks(peer)
            buffer_rec%indx_offset(peer + 1) = buffer_rec%indx_offset(peer) + recv_blocks(peer)
         END IF
      END DO
      ALLOCATE (buffer_send%msg(buffer_send%msg_offset(nprocs)), &
                buffer_rec%msg(buffer_rec%msg_offset(nprocs)), &
                buffer_send%indx(n_index, buffer_send%indx_offset(nprocs)), &
                buffer_rec%indx(n_index, buffer_rec%indx_offset(nprocs)))
   END SUBROUTINE allocate_matrix_buffers

! **************************************************************************************************
!> \brief Synchronously exchange flat matrix buffers with nonempty peers.
!> \param comm communicator covering all sending and receiving ranks
!> \param buffer_rec preallocated receive buffer
!> \param buffer_send populated send buffer
!> \param do_indx whether to exchange integer metadata (default: true)
!> \param do_msg whether to exchange numerical payloads (default: true)
! **************************************************************************************************
   SUBROUTINE exchange_matrix_buffers(comm, buffer_rec, buffer_send, do_indx, do_msg)
      CLASS(mp_comm_type), INTENT(IN)                    :: comm
      TYPE(matrix_buffer_type), INTENT(INOUT)            :: buffer_rec, buffer_send
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_indx, do_msg

      CHARACTER(LEN=*), PARAMETER :: routineN = 'exchange_matrix_buffers'

      INTEGER                                            :: first, handle, last, nprocs, nreq, peer
      LOGICAL                                            :: my_do_indx, my_do_msg
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: requests

      CALL timeset(routineN, handle)
      nprocs = comm%num_pe
      my_do_indx = .TRUE.
      IF (PRESENT(do_indx)) my_do_indx = do_indx
      my_do_msg = .TRUE.
      IF (PRESENT(do_msg)) my_do_msg = do_msg

      ! Copy the self segment directly, including the single-rank case.
      peer = comm%mepos
      IF (my_do_indx) THEN
         buffer_rec%indx(:, buffer_rec%indx_offset(peer) + 1:buffer_rec%indx_offset(peer + 1)) = &
            buffer_send%indx(:, buffer_send%indx_offset(peer) + 1:buffer_send%indx_offset(peer + 1))
      END IF
      IF (my_do_msg) THEN
         buffer_rec%msg(buffer_rec%msg_offset(peer) + 1:buffer_rec%msg_offset(peer + 1)) = &
            buffer_send%msg(buffer_send%msg_offset(peer) + 1:buffer_send%msg_offset(peer + 1))
      END IF

      IF (nprocs > 1) THEN
         nreq = 0
         IF (my_do_indx) THEN
            nreq = COUNT(buffer_rec%indx_offset(1:) > buffer_rec%indx_offset(:nprocs - 1)) + &
                   COUNT(buffer_send%indx_offset(1:) > buffer_send%indx_offset(:nprocs - 1))
         END IF
         IF (my_do_msg) THEN
            nreq = nreq + COUNT(buffer_rec%msg_offset(1:) > buffer_rec%msg_offset(:nprocs - 1)) + &
                   COUNT(buffer_send%msg_offset(1:) > buffer_send%msg_offset(:nprocs - 1))
         END IF
         ALLOCATE (requests(nreq))
         nreq = 0

         ! Post both channels before waiting so metadata and payload can overlap.
         DO peer = 0, nprocs - 1
            IF (peer == comm%mepos) CYCLE
            IF (my_do_indx) THEN
               first = buffer_rec%indx_offset(peer) + 1
               last = buffer_rec%indx_offset(peer + 1)
               IF (last >= first) THEN
                  nreq = nreq + 1
                  CALL comm%irecv(buffer_rec%indx(:, first:last), peer, requests(nreq), tag=4)
               END IF
            END IF
            IF (my_do_msg .AND. buffer_rec%msg_offset(peer + 1) > buffer_rec%msg_offset(peer)) THEN
               first = buffer_rec%msg_offset(peer) + 1
               last = buffer_rec%msg_offset(peer + 1)
               nreq = nreq + 1
               CALL comm%irecv(buffer_rec%msg(first:last), peer, requests(nreq), tag=7)
            END IF
         END DO

         DO peer = 0, nprocs - 1
            IF (peer == comm%mepos) CYCLE
            IF (my_do_indx) THEN
               first = buffer_send%indx_offset(peer) + 1
               last = buffer_send%indx_offset(peer + 1)
               IF (last >= first) THEN
                  nreq = nreq + 1
                  CALL comm%isend(buffer_send%indx(:, first:last), peer, requests(nreq), tag=4)
               END IF
            END IF
            IF (my_do_msg .AND. buffer_send%msg_offset(peer + 1) > buffer_send%msg_offset(peer)) THEN
               first = buffer_send%msg_offset(peer) + 1
               last = buffer_send%msg_offset(peer + 1)
               nreq = nreq + 1
               CALL comm%isend(buffer_send%msg(first:last), peer, requests(nreq), tag=7)
            END IF
         END DO
         IF (nreq > 0) CALL mp_waitall(requests(:nreq))
      END IF
      CALL timestop(handle)
   END SUBROUTINE exchange_matrix_buffers

END MODULE matrix_communication
