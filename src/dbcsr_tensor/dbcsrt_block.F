!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods to operate on n-dimensional tensor blocks.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsrt_block

   USE allocate_wrap,                   ONLY: allocate_any
   USE dbcsr_api,                       ONLY: dbcsr_type_complex_4,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_8
   USE dbcsrt_block_c,                  ONLY: block_nd_c,&
                                              reshape_2d_to_nd_block,&
                                              reshape_nd_to_2d_block
   USE dbcsrt_block_d,                  ONLY: block_nd_d,&
                                              reshape_2d_to_nd_block,&
                                              reshape_nd_to_2d_block
   USE dbcsrt_block_s,                  ONLY: block_nd_s,&
                                              reshape_2d_to_nd_block,&
                                              reshape_nd_to_2d_block
   USE dbcsrt_block_z,                  ONLY: block_nd_z,&
                                              reshape_2d_to_nd_block,&
                                              reshape_nd_to_2d_block
   USE kinds,                           ONLY: real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_environ,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_waitall
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsrt_block'

   PUBLIC :: &
      block_buffer_add_anyd_block, &
      block_buffer_add_block, &
      block_buffer_blocks_left, &
      block_buffer_create, &
      block_buffer_destroy, &
      block_buffer_get_next_anyd_block, &
      block_buffer_iterator_reset, &
      block_buffer_type, &
      block_nd, &
      block_size, &
      communicate_buffer, &
      create_block, &
      ndims_buffer, &
      destroy_block, &
      reshape_2d_to_nd_block, &
      reshape_nd_to_2d_block

   TYPE block_nd
      TYPE(block_nd_d) :: r_dp
      TYPE(block_nd_s) :: r_sp
      TYPE(block_nd_z) :: c_dp
      TYPE(block_nd_c) :: c_sp
      INTEGER          :: data_type = -1
   END TYPE

   INTERFACE create_block
      MODULE PROCEDURE create_block_data_d
      MODULE PROCEDURE create_block_data_s
      MODULE PROCEDURE create_block_data_z
      MODULE PROCEDURE create_block_data_c
      MODULE PROCEDURE create_block_nodata
   END INTERFACE

   TYPE block_buffer_type
      INTEGER                                          :: ndim = -1
      INTEGER                                          :: nblock = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE            :: indx
      REAL(kind=real_8), DIMENSION(:), ALLOCATABLE     :: msg_r_dp
      REAL(kind=real_4), DIMENSION(:), ALLOCATABLE     :: msg_r_sp
      COMPLEX(kind=real_8), DIMENSION(:), ALLOCATABLE  :: msg_c_dp
      COMPLEX(kind=real_4), DIMENSION(:), ALLOCATABLE  :: msg_c_sp
      INTEGER                                          :: data_type = -1
      INTEGER                                          :: endpos = -1
   END TYPE

   INTERFACE block_buffer_add_block
      MODULE PROCEDURE block_buffer_add_block_d
      MODULE PROCEDURE block_buffer_add_block_s
      MODULE PROCEDURE block_buffer_add_block_z
      MODULE PROCEDURE block_buffer_add_block_c
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param block ...
!> \param sizes ...
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata(block, sizes, data_type)
      TYPE(block_nd), INTENT(OUT)                        :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      INTEGER, INTENT(IN)                                :: data_type

      block%data_type = data_type
      SELECT CASE (data_type)
      CASE (dbcsr_type_real_8)
         CALL create_block_nodata_d(block%r_dp, sizes)
      CASE (dbcsr_type_real_4)
         CALL create_block_nodata_s(block%r_sp, sizes)
      CASE (dbcsr_type_complex_8)
         CALL create_block_nodata_z(block%c_dp, sizes)
      CASE (dbcsr_type_complex_4)
         CALL create_block_nodata_c(block%c_sp, sizes)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block(block)
      TYPE(block_nd), INTENT(INOUT)                      :: block

      SELECT CASE (block%data_type)
      CASE (dbcsr_type_real_8)
         CALL destroy_block_d(block%r_dp)
      CASE (dbcsr_type_real_4)
         CALL destroy_block_s(block%r_sp)
      CASE (dbcsr_type_complex_8)
         CALL destroy_block_z(block%c_dp)
      CASE (dbcsr_type_complex_4)
         CALL destroy_block_c(block%c_sp)
      END SELECT

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block ...
!> \retval blk_size ...
! **************************************************************************************************
   PURE FUNCTION block_size(block)
      TYPE(block_nd), INTENT(IN)                         :: block
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: block_size

      SELECT CASE (block%data_type)
      CASE (dbcsr_type_real_8)
         CALL allocate_any(SHAPE(block%r_dp%sizes), block_size, source=block%r_dp%sizes)
      CASE (dbcsr_type_real_4)
         CALL allocate_any(SHAPE(block%r_sp%sizes), block_size, source=block%r_sp%sizes)
      CASE (dbcsr_type_complex_8)
         CALL allocate_any(SHAPE(block%c_dp%sizes), block_size, source=block%c_dp%sizes)
      CASE (dbcsr_type_complex_4)
         CALL allocate_any(SHAPE(block%c_sp%sizes), block_size, source=block%c_sp%sizes)
      END SELECT
   END FUNCTION

! **************************************************************************************************
!> \brief Create block buffer for MPI communication.
!> \param buffer block buffer
!> \param nblock number of blocks
!> \param ndata total number of block entries
!> \param data_type ...
!> \param ndim number of dimensions
! **************************************************************************************************
   SUBROUTINE block_buffer_create(buffer, nblock, ndata, data_type, ndim)
      TYPE(block_buffer_type), INTENT(OUT)               :: buffer
      INTEGER, INTENT(IN)                                :: nblock, ndata, data_type, ndim

      buffer%nblock = nblock
      buffer%data_type = data_type
      buffer%endpos = 0
      buffer%ndim = ndim
      SELECT CASE (data_type)
      CASE (dbcsr_type_real_4)
         ALLOCATE (buffer%msg_r_sp(ndata))
      CASE (dbcsr_type_real_8)
         ALLOCATE (buffer%msg_r_dp(ndata))
      CASE (dbcsr_type_complex_4)
         ALLOCATE (buffer%msg_c_sp(ndata))
      CASE (dbcsr_type_complex_8)
         ALLOCATE (buffer%msg_c_dp(ndata))
      END SELECT
      ALLOCATE (buffer%indx(nblock, ndim+1))
   END SUBROUTINE block_buffer_create

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE block_buffer_destroy(buffer)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer

      SELECT CASE (buffer%data_type)
      CASE (dbcsr_type_real_4)
         DEALLOCATE (buffer%msg_r_sp)
      CASE (dbcsr_type_real_8)
         DEALLOCATE (buffer%msg_r_dp)
      CASE (dbcsr_type_complex_4)
         DEALLOCATE (buffer%msg_c_sp)
      CASE (dbcsr_type_complex_8)
         DEALLOCATE (buffer%msg_c_dp)
      END SELECT
      DEALLOCATE (buffer%indx)
      buffer%nblock = -1
      buffer%data_type = -1
      buffer%ndim = -1
      buffer%endpos = -1
   END SUBROUTINE block_buffer_destroy

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \retval ndims_buffer ...
! **************************************************************************************************
   PURE FUNCTION ndims_buffer(buffer)
      TYPE(block_buffer_type), INTENT(IN)                :: buffer
      INTEGER                                            :: ndims_buffer

      ndims_buffer = buffer%ndim
   END FUNCTION

! **************************************************************************************************
!> \brief insert a block into block buffer (at current iterator position)
!> \param buffer ...
!> \param index index of block
!> \param block block
! **************************************************************************************************
   SUBROUTINE block_buffer_add_anyd_block(buffer, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                                      :: index
      TYPE(block_nd), INTENT(IN)                         :: block

      SELECT CASE (block%data_type)
      CASE (dbcsr_type_real_8)
         CALL block_buffer_add_block_d(buffer, SIZE(block%r_dp%blk), index, block%r_dp%blk)
      CASE (dbcsr_type_real_4)
         CALL block_buffer_add_block_s(buffer, SIZE(block%r_sp%blk), index, block%r_sp%blk)
      CASE (dbcsr_type_complex_8)
         CALL block_buffer_add_block_z(buffer, SIZE(block%c_dp%blk), index, block%c_dp%blk)
      CASE (dbcsr_type_complex_4)
         CALL block_buffer_add_block_c(buffer, SIZE(block%c_sp%blk), index, block%c_sp%blk)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_anyd_block(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(OUT)                               :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                                     :: index
      TYPE(block_nd), INTENT(INOUT), OPTIONAL            :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: advance_iter

      SELECT CASE (buffer%data_type)
      CASE (dbcsr_type_real_8)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_d(buffer, ndata, index, block%r_dp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_d(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      CASE (dbcsr_type_real_4)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_s(buffer, ndata, index, block%r_sp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_s(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      CASE (dbcsr_type_complex_8)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_z(buffer, ndata, index, block%c_dp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_z(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      CASE (dbcsr_type_complex_4)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_c(buffer, ndata, index, block%c_sp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_c(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE block_buffer_iterator_reset(buffer)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer

      buffer%endpos = 0
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \retval block_buffer_blocks_left ...
! **************************************************************************************************
   PURE FUNCTION block_buffer_blocks_left(buffer)
      TYPE(block_buffer_type), INTENT(IN)                :: buffer
      LOGICAL                                            :: block_buffer_blocks_left

      block_buffer_blocks_left = buffer%endpos .LT. buffer%nblock
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param mp_comm ...
!> \param buffer_recv ...
!> \param buffer_send ...
!> \param req_array ...
! **************************************************************************************************
   SUBROUTINE communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      INTEGER, INTENT(IN)                                :: mp_comm
      TYPE(block_buffer_type), DIMENSION(0:)             :: buffer_recv, buffer_send
      INTEGER, DIMENSION(:, :)                           :: req_array

      INTEGER                                            :: iproc, mynode, numnodes, rec_counter, &
                                                            send_counter

      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (numnodes > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO iproc = 0, numnodes-1
            IF (buffer_recv(iproc)%nblock > 0) THEN
               rec_counter = rec_counter+1
               CALL mp_irecv(buffer_recv(iproc)%indx, iproc, mp_comm, req_array(rec_counter, 3), tag=4)
            END IF
            IF (buffer_recv(iproc)%nblock > 0) THEN
               SELECT CASE (buffer_recv (iproc)%data_type)
               CASE (dbcsr_type_real_8)
                  CALL mp_irecv(buffer_recv(iproc)%msg_r_dp, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
               CASE (dbcsr_type_real_4)
                  CALL mp_irecv(buffer_recv(iproc)%msg_r_sp, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
               CASE (dbcsr_type_complex_8)
                  CALL mp_irecv(buffer_recv(iproc)%msg_c_dp, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
               CASE (dbcsr_type_complex_4)
                  CALL mp_irecv(buffer_recv(iproc)%msg_c_sp, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
               END SELECT
            END IF
         END DO

         DO iproc = 0, numnodes-1
            IF (buffer_send(iproc)%nblock > 0) THEN
               send_counter = send_counter+1
               CALL mp_isend(buffer_send(iproc)%indx, iproc, mp_comm, req_array(send_counter, 1), tag=4)
            END IF
            IF (buffer_send(iproc)%nblock > 0) THEN
               SELECT CASE (buffer_recv (iproc)%data_type)
               CASE (dbcsr_type_real_8)
                  CALL mp_isend(buffer_send(iproc)%msg_r_dp, iproc, mp_comm, req_array(send_counter, 2), tag=7)
               CASE (dbcsr_type_real_4)
                  CALL mp_isend(buffer_send(iproc)%msg_r_sp, iproc, mp_comm, req_array(send_counter, 2), tag=7)
               CASE (dbcsr_type_complex_8)
                  CALL mp_isend(buffer_send(iproc)%msg_c_dp, iproc, mp_comm, req_array(send_counter, 2), tag=7)
               CASE (dbcsr_type_complex_4)
                  CALL mp_isend(buffer_send(iproc)%msg_c_sp, iproc, mp_comm, req_array(send_counter, 2), tag=7)
               END SELECT
            END IF
         END DO

         IF (send_counter > 0) THEN
            CALL mp_waitall(req_array(1:send_counter, 1:2))
         ENDIF
         IF (rec_counter > 0) THEN
            CALL mp_waitall(req_array(1:rec_counter, 3:4))
         ENDIF

      ELSE
         IF (buffer_recv(0)%nblock > 0) THEN
            buffer_recv(0)%indx(:, :) = buffer_send(0)%indx(:, :)
            SELECT CASE (buffer_recv (0)%data_type)
            CASE (dbcsr_type_real_8)
               buffer_recv(0)%msg_r_dp(:) = buffer_send(0)%msg_r_dp(:)
            CASE (dbcsr_type_real_4)
               buffer_recv(0)%msg_r_sp(:) = buffer_send(0)%msg_r_sp(:)
            CASE (dbcsr_type_complex_8)
               buffer_recv(0)%msg_c_dp(:) = buffer_send(0)%msg_c_dp(:)
            CASE (dbcsr_type_complex_4)
               buffer_recv(0)%msg_c_sp(:) = buffer_send(0)%msg_c_sp(:)
            END SELECT
         ENDIF
      ENDIF

   END SUBROUTINE

#include "dbcsrt_block_inc_c.f90"
#include "dbcsrt_block_inc_d.f90"
#include "dbcsrt_block_inc_s.f90"
#include "dbcsrt_block_inc_z.f90"
END MODULE dbcsrt_block
