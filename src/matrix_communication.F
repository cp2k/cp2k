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
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_type,&
                                              mp_waitall
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'matrix_communication'
   PUBLIC :: matrix_buffer_type, exchange_matrix_buffers

   TYPE matrix_buffer_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: msg
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: sizes
      INTEGER, DIMENSION(:, :), ALLOCATABLE              :: indx
      INTEGER                                            :: proc = -1
      TYPE(mp_request_type)                              :: msg_req = mp_request_type()
   END TYPE matrix_buffer_type

CONTAINS

! **************************************************************************************************
!> \brief Synchronously exchange per-rank buffers.
!> \param para_env communicator covering all sending and receiving ranks
!> \param num_entries_rec payload lengths received from each rank, indexed from zero
!> \param num_entries_send payload lengths sent to each rank, indexed from zero
!> \param buffer_rec preallocated receive buffers
!> \param buffer_send populated send buffers
!> \param do_indx whether to exchange integer metadata (default: true)
!> \param do_msg whether to exchange numerical payloads (default: true)
! **************************************************************************************************
   SUBROUTINE exchange_matrix_buffers(para_env, num_entries_rec, num_entries_send, buffer_rec, buffer_send, &
                                      do_indx, do_msg)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: num_entries_rec, num_entries_send
      TYPE(matrix_buffer_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(INOUT)                     :: buffer_rec, buffer_send
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_indx, do_msg

      CHARACTER(LEN=*), PARAMETER :: routineN = 'exchange_matrix_buffers'

      INTEGER                                            :: handle, imepos, rec_counter, send_counter
      LOGICAL                                            :: my_do_indx, my_do_msg
      TYPE(mp_request_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: req_array

      CALL timeset(routineN, handle)

      ALLOCATE (req_array(para_env%num_pe, 4))

      my_do_indx = .TRUE.
      IF (PRESENT(do_indx)) my_do_indx = do_indx
      my_do_msg = .TRUE.
      IF (PRESENT(do_msg)) my_do_msg = do_msg

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter + 1
               IF (my_do_indx) THEN
                  CALL para_env%irecv(buffer_rec(imepos)%indx, imepos, req_array(rec_counter, 3), tag=4)
               END IF
               IF (my_do_msg) THEN
                  CALL para_env%irecv(buffer_rec(imepos)%msg, imepos, req_array(rec_counter, 4), tag=7)
               END IF
            END IF
         END DO

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter + 1
               IF (my_do_indx) THEN
                  CALL para_env%isend(buffer_send(imepos)%indx, imepos, req_array(send_counter, 1), tag=4)
               END IF
               IF (my_do_msg) THEN
                  CALL para_env%isend(buffer_send(imepos)%msg, imepos, req_array(send_counter, 2), tag=7)
               END IF
            END IF
         END DO

         IF (my_do_indx) THEN
            CALL mp_waitall(req_array(1:send_counter, 1))
            CALL mp_waitall(req_array(1:rec_counter, 3))
         END IF

         IF (my_do_msg) THEN
            CALL mp_waitall(req_array(1:send_counter, 2))
            CALL mp_waitall(req_array(1:rec_counter, 4))
         END IF

      ELSE

         IF (my_do_indx) buffer_rec(0)%indx(:, :) = buffer_send(0)%indx
         IF (my_do_msg) buffer_rec(0)%msg(:) = buffer_send(0)%msg

      END IF

      CALL timestop(handle)

   END SUBROUTINE exchange_matrix_buffers

END MODULE matrix_communication
