!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Generated from dbcsr_tensor_block.fypp using Fypp.                                               !
! **DO NOT** modify this file, edit dbcsr_tensor_block.fypp instead.                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods to operate on n-dimensional tensor blocks.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_block

   USE allocate_wrap,                   ONLY: allocate_any
   USE dbcsr_api,                       ONLY: dbcsr_type_complex_4,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_8
   USE dbcsr_tensor_index,              ONLY: nd_to_2d_mapping,&
                                              ndims_mapping
   USE kinds,                           ONLY: real_4,&
                                              real_8
   USE message_passing,                 ONLY: mp_environ,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_waitall
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_block'

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

   PUBLIC :: block_nd_r_dp
   PUBLIC :: block_nd_r_sp
   PUBLIC :: block_nd_c_dp
   PUBLIC :: block_nd_c_sp

   TYPE block_nd_r_dp
      INTEGER, DIMENSION(:), ALLOCATABLE   :: sizes
      REAL(kind=real_8), DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

   TYPE block_nd_r_sp
      INTEGER, DIMENSION(:), ALLOCATABLE   :: sizes
      REAL(kind=real_4), DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

   TYPE block_nd_c_dp
      INTEGER, DIMENSION(:), ALLOCATABLE   :: sizes
      COMPLEX(kind=real_8), DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

   TYPE block_nd_c_sp
      INTEGER, DIMENSION(:), ALLOCATABLE   :: sizes
      COMPLEX(kind=real_4), DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

   TYPE block_nd
      TYPE(block_nd_r_dp) :: r_dp
      TYPE(block_nd_r_sp) :: r_sp
      TYPE(block_nd_c_dp) :: c_dp
      TYPE(block_nd_c_sp) :: c_sp
      INTEGER          :: data_type = -1
   END TYPE

   TYPE block_buffer_type
      INTEGER                                    :: ndim = -1
      INTEGER                                    :: nblock = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE      :: indx
      REAL(kind=real_8), DIMENSION(:), ALLOCATABLE       :: msg_r_dp
      REAL(kind=real_4), DIMENSION(:), ALLOCATABLE       :: msg_r_sp
      COMPLEX(kind=real_8), DIMENSION(:), ALLOCATABLE       :: msg_c_dp
      COMPLEX(kind=real_4), DIMENSION(:), ALLOCATABLE       :: msg_c_sp
      INTEGER                                    :: data_type = -1
      INTEGER                                    :: endpos = -1
   END TYPE

   INTERFACE create_block
      MODULE PROCEDURE create_block_data_r_dp
      MODULE PROCEDURE create_block_data_r_sp
      MODULE PROCEDURE create_block_data_c_dp
      MODULE PROCEDURE create_block_data_c_sp
      MODULE PROCEDURE create_block_nodata
   END INTERFACE

   INTERFACE block_buffer_add_block
      MODULE PROCEDURE block_buffer_add_block_r_dp
      MODULE PROCEDURE block_buffer_add_block_r_sp
      MODULE PROCEDURE block_buffer_add_block_c_dp
      MODULE PROCEDURE block_buffer_add_block_c_sp
   END INTERFACE

   INTERFACE reshape_nd_to_2d_block
      MODULE PROCEDURE reshape_2d_to_2d_block_r_dp
      MODULE PROCEDURE reshape_3d_to_2d_block_r_dp
      MODULE PROCEDURE reshape_4d_to_2d_block_r_dp
      MODULE PROCEDURE reshape_2d_to_2d_block_r_sp
      MODULE PROCEDURE reshape_3d_to_2d_block_r_sp
      MODULE PROCEDURE reshape_4d_to_2d_block_r_sp
      MODULE PROCEDURE reshape_2d_to_2d_block_c_dp
      MODULE PROCEDURE reshape_3d_to_2d_block_c_dp
      MODULE PROCEDURE reshape_4d_to_2d_block_c_dp
      MODULE PROCEDURE reshape_2d_to_2d_block_c_sp
      MODULE PROCEDURE reshape_3d_to_2d_block_c_sp
      MODULE PROCEDURE reshape_4d_to_2d_block_c_sp
   END INTERFACE reshape_nd_to_2d_block

   INTERFACE reshape_2d_to_nd_block
      MODULE PROCEDURE reshape_2d_to_2d_block_r_dp
      MODULE PROCEDURE reshape_2d_to_3d_block_r_dp
      MODULE PROCEDURE reshape_2d_to_4d_block_r_dp
      MODULE PROCEDURE reshape_2d_to_2d_block_r_sp
      MODULE PROCEDURE reshape_2d_to_3d_block_r_sp
      MODULE PROCEDURE reshape_2d_to_4d_block_r_sp
      MODULE PROCEDURE reshape_2d_to_2d_block_c_dp
      MODULE PROCEDURE reshape_2d_to_3d_block_c_dp
      MODULE PROCEDURE reshape_2d_to_4d_block_c_dp
      MODULE PROCEDURE reshape_2d_to_2d_block_c_sp
      MODULE PROCEDURE reshape_2d_to_3d_block_c_sp
      MODULE PROCEDURE reshape_2d_to_4d_block_c_sp
   END INTERFACE reshape_2d_to_nd_block

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
         CALL create_block_nodata_r_dp(block%r_dp, sizes)
      CASE (dbcsr_type_real_4)
         CALL create_block_nodata_r_sp(block%r_sp, sizes)
      CASE (dbcsr_type_complex_8)
         CALL create_block_nodata_c_dp(block%c_dp, sizes)
      CASE (dbcsr_type_complex_4)
         CALL create_block_nodata_c_sp(block%c_sp, sizes)
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
         CALL destroy_block_r_dp(block%r_dp)
      CASE (dbcsr_type_real_4)
         CALL destroy_block_r_sp(block%r_sp)
      CASE (dbcsr_type_complex_8)
         CALL destroy_block_c_dp(block%c_dp)
      CASE (dbcsr_type_complex_4)
         CALL destroy_block_c_sp(block%c_sp)
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
      CASE (dbcsr_type_real_8)
         ALLOCATE (buffer%msg_r_dp(ndata))
      CASE (dbcsr_type_real_4)
         ALLOCATE (buffer%msg_r_sp(ndata))
      CASE (dbcsr_type_complex_8)
         ALLOCATE (buffer%msg_c_dp(ndata))
      CASE (dbcsr_type_complex_4)
         ALLOCATE (buffer%msg_c_sp(ndata))
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
      CASE (dbcsr_type_real_8)
         DEALLOCATE (buffer%msg_r_dp)
      CASE (dbcsr_type_real_4)
         DEALLOCATE (buffer%msg_r_sp)
      CASE (dbcsr_type_complex_8)
         DEALLOCATE (buffer%msg_c_dp)
      CASE (dbcsr_type_complex_4)
         DEALLOCATE (buffer%msg_c_sp)
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
         CALL block_buffer_add_block_r_dp(buffer, SIZE(block%r_dp%blk), index, block%r_dp%blk)
      CASE (dbcsr_type_real_4)
         CALL block_buffer_add_block_r_sp(buffer, SIZE(block%r_sp%blk), index, block%r_sp%blk)
      CASE (dbcsr_type_complex_8)
         CALL block_buffer_add_block_c_dp(buffer, SIZE(block%c_dp%blk), index, block%c_dp%blk)
      CASE (dbcsr_type_complex_4)
         CALL block_buffer_add_block_c_sp(buffer, SIZE(block%c_sp%blk), index, block%c_sp%blk)
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
            CALL block_buffer_get_next_block_r_dp(buffer, ndata, index, block%r_dp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_r_dp(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      CASE (dbcsr_type_real_4)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_r_sp(buffer, ndata, index, block%r_sp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_r_sp(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      CASE (dbcsr_type_complex_8)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_c_dp(buffer, ndata, index, block%c_dp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_c_dp(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
      CASE (dbcsr_type_complex_4)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_c_sp(buffer, ndata, index, block%c_sp%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_c_sp(buffer, ndata, index, advance_iter=advance_iter)
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

! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \param block ...
!> \param sizes ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE create_block_data_r_dp(block, sizes, array)
      TYPE(block_nd), INTENT(OUT)                        :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      REAL(kind=real_8), DIMENSION(PRODUCT(sizes)), &
         INTENT(IN)                                      :: array

      ASSOCIATE(blk=>block%r_dp)
      block%data_type = dbcsr_type_real_8
      ALLOCATE (blk%sizes(SIZE(sizes)))
      blk%sizes(:) = sizes(:)
      ALLOCATE (blk%blk(PRODUCT(sizes)))
      blk%blk(:) = array(:)
      END ASSOCIATE
   END SUBROUTINE
! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \param block ...
!> \param sizes ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE create_block_data_r_sp(block, sizes, array)
      TYPE(block_nd), INTENT(OUT)                        :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      REAL(kind=real_4), DIMENSION(PRODUCT(sizes)), &
         INTENT(IN)                                      :: array

      ASSOCIATE(blk=>block%r_sp)
      block%data_type = dbcsr_type_real_4
      ALLOCATE (blk%sizes(SIZE(sizes)))
      blk%sizes(:) = sizes(:)
      ALLOCATE (blk%blk(PRODUCT(sizes)))
      blk%blk(:) = array(:)
      END ASSOCIATE
   END SUBROUTINE
! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \param block ...
!> \param sizes ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE create_block_data_c_dp(block, sizes, array)
      TYPE(block_nd), INTENT(OUT)                        :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_8), DIMENSION(PRODUCT(sizes)), &
         INTENT(IN)                                      :: array

      ASSOCIATE(blk=>block%c_dp)
      block%data_type = dbcsr_type_complex_8
      ALLOCATE (blk%sizes(SIZE(sizes)))
      blk%sizes(:) = sizes(:)
      ALLOCATE (blk%blk(PRODUCT(sizes)))
      blk%blk(:) = array(:)
      END ASSOCIATE
   END SUBROUTINE
! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \param block ...
!> \param sizes ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE create_block_data_c_sp(block, sizes, array)
      TYPE(block_nd), INTENT(OUT)                        :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      COMPLEX(kind=real_4), DIMENSION(PRODUCT(sizes)), &
         INTENT(IN)                                      :: array

      ASSOCIATE(blk=>block%c_sp)
      block%data_type = dbcsr_type_complex_4
      ALLOCATE (blk%sizes(SIZE(sizes)))
      blk%sizes(:) = sizes(:)
      ALLOCATE (blk%blk(PRODUCT(sizes)))
      blk%blk(:) = array(:)
      END ASSOCIATE
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \param block ...
!> \param sizes ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata_r_dp(block, sizes)
      TYPE(block_nd_r_dp), INTENT(OUT)                   :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes

      ALLOCATE (block%sizes(SIZE(sizes)))
      block%sizes(:) = sizes(:)
      ALLOCATE (block%blk(PRODUCT(sizes)))
   END SUBROUTINE
! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \param block ...
!> \param sizes ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata_r_sp(block, sizes)
      TYPE(block_nd_r_sp), INTENT(OUT)                   :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes

      ALLOCATE (block%sizes(SIZE(sizes)))
      block%sizes(:) = sizes(:)
      ALLOCATE (block%blk(PRODUCT(sizes)))
   END SUBROUTINE
! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \param block ...
!> \param sizes ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata_c_dp(block, sizes)
      TYPE(block_nd_c_dp), INTENT(OUT)                   :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes

      ALLOCATE (block%sizes(SIZE(sizes)))
      block%sizes(:) = sizes(:)
      ALLOCATE (block%blk(PRODUCT(sizes)))
   END SUBROUTINE
! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \param block ...
!> \param sizes ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata_c_sp(block, sizes)
      TYPE(block_nd_c_sp), INTENT(OUT)                   :: block
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes

      ALLOCATE (block%sizes(SIZE(sizes)))
      block%sizes(:) = sizes(:)
      ALLOCATE (block%blk(PRODUCT(sizes)))
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block_r_dp(block)
      TYPE(block_nd_r_dp), INTENT(INOUT)                 :: block

      DEALLOCATE (block%blk)
      DEALLOCATE (block%sizes)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block_r_sp(block)
      TYPE(block_nd_r_sp), INTENT(INOUT)                 :: block

      DEALLOCATE (block%blk)
      DEALLOCATE (block%sizes)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block_c_dp(block)
      TYPE(block_nd_c_dp), INTENT(INOUT)                 :: block

      DEALLOCATE (block%blk)
      DEALLOCATE (block%sizes)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block_c_sp(block)
      TYPE(block_nd_c_sp), INTENT(INOUT)                 :: block

      DEALLOCATE (block%blk)
      DEALLOCATE (block%sizes)
   END SUBROUTINE

! **************************************************************************************************
!> \brief add block to buffer.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block_r_dp(buffer, ndata, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(IN)                                :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                                      :: index
      REAL(kind=real_8), DIMENSION(ndata), INTENT(IN)    :: block

      INTEGER                                            :: ndims, p, p_data

      CPASSERT(buffer%data_type .EQ. dbcsr_type_real_8)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF

      buffer%msg_r_dp(p_data+1:p_data+ndata) = block(:)
      buffer%indx(p+1, 1:ndims) = INDEX(:)
      IF (p > 0) THEN
         buffer%indx(p+1, ndims+1) = buffer%indx(p, ndims+1)+ndata
      ELSE
         buffer%indx(p+1, ndims+1) = ndata
      ENDIF
      buffer%endpos = buffer%endpos+1
   END SUBROUTINE
! **************************************************************************************************
!> \brief add block to buffer.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block_r_sp(buffer, ndata, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(IN)                                :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                                      :: index
      REAL(kind=real_4), DIMENSION(ndata), INTENT(IN)    :: block

      INTEGER                                            :: ndims, p, p_data

      CPASSERT(buffer%data_type .EQ. dbcsr_type_real_4)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF

      buffer%msg_r_sp(p_data+1:p_data+ndata) = block(:)
      buffer%indx(p+1, 1:ndims) = INDEX(:)
      IF (p > 0) THEN
         buffer%indx(p+1, ndims+1) = buffer%indx(p, ndims+1)+ndata
      ELSE
         buffer%indx(p+1, ndims+1) = ndata
      ENDIF
      buffer%endpos = buffer%endpos+1
   END SUBROUTINE
! **************************************************************************************************
!> \brief add block to buffer.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block_c_dp(buffer, ndata, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(IN)                                :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                                      :: index
      COMPLEX(kind=real_8), DIMENSION(ndata), INTENT(IN) :: block

      INTEGER                                            :: ndims, p, p_data

      CPASSERT(buffer%data_type .EQ. dbcsr_type_complex_8)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF

      buffer%msg_c_dp(p_data+1:p_data+ndata) = block(:)
      buffer%indx(p+1, 1:ndims) = INDEX(:)
      IF (p > 0) THEN
         buffer%indx(p+1, ndims+1) = buffer%indx(p, ndims+1)+ndata
      ELSE
         buffer%indx(p+1, ndims+1) = ndata
      ENDIF
      buffer%endpos = buffer%endpos+1
   END SUBROUTINE
! **************************************************************************************************
!> \brief add block to buffer.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block_c_sp(buffer, ndata, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(IN)                                :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                                      :: index
      COMPLEX(kind=real_4), DIMENSION(ndata), INTENT(IN) :: block

      INTEGER                                            :: ndims, p, p_data

      CPASSERT(buffer%data_type .EQ. dbcsr_type_complex_4)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF

      buffer%msg_c_sp(p_data+1:p_data+ndata) = block(:)
      buffer%indx(p+1, 1:ndims) = INDEX(:)
      IF (p > 0) THEN
         buffer%indx(p+1, ndims+1) = buffer%indx(p, ndims+1)+ndata
      ELSE
         buffer%indx(p+1, ndims+1) = ndata
      ENDIF
      buffer%endpos = buffer%endpos+1
   END SUBROUTINE

! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block_r_dp(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(OUT)                               :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                                     :: index
      REAL(kind=real_8), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: advance_iter

      INTEGER                                            :: ndims, p, p_data
      LOGICAL                                            :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      CPASSERT(buffer%data_type .EQ. dbcsr_type_real_8)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF
      IF (p > 0) THEN
         ndata = buffer%indx(p+1, ndims+1)-buffer%indx(p, ndims+1)
      ELSE
         ndata = buffer%indx(p+1, ndims+1)
      ENDIF
      INDEX(:) = buffer%indx(p+1, 1:ndims)
      IF (PRESENT(block)) THEN
         ALLOCATE (block(ndata))
         block(:) = buffer%msg_r_dp(p_data+1:p_data+ndata)
      ENDIF

      IF (do_advance) buffer%endpos = buffer%endpos+1
   END SUBROUTINE
! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block_r_sp(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(OUT)                               :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                                     :: index
      REAL(kind=real_4), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: advance_iter

      INTEGER                                            :: ndims, p, p_data
      LOGICAL                                            :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      CPASSERT(buffer%data_type .EQ. dbcsr_type_real_4)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF
      IF (p > 0) THEN
         ndata = buffer%indx(p+1, ndims+1)-buffer%indx(p, ndims+1)
      ELSE
         ndata = buffer%indx(p+1, ndims+1)
      ENDIF
      INDEX(:) = buffer%indx(p+1, 1:ndims)
      IF (PRESENT(block)) THEN
         ALLOCATE (block(ndata))
         block(:) = buffer%msg_r_sp(p_data+1:p_data+ndata)
      ENDIF

      IF (do_advance) buffer%endpos = buffer%endpos+1
   END SUBROUTINE
! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block_c_dp(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(OUT)                               :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                                     :: index
      COMPLEX(kind=real_8), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: advance_iter

      INTEGER                                            :: ndims, p, p_data
      LOGICAL                                            :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      CPASSERT(buffer%data_type .EQ. dbcsr_type_complex_8)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF
      IF (p > 0) THEN
         ndata = buffer%indx(p+1, ndims+1)-buffer%indx(p, ndims+1)
      ELSE
         ndata = buffer%indx(p+1, ndims+1)
      ENDIF
      INDEX(:) = buffer%indx(p+1, 1:ndims)
      IF (PRESENT(block)) THEN
         ALLOCATE (block(ndata))
         block(:) = buffer%msg_c_dp(p_data+1:p_data+ndata)
      ENDIF

      IF (do_advance) buffer%endpos = buffer%endpos+1
   END SUBROUTINE
! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
!> \param advance_iter ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block_c_sp(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)             :: buffer
      INTEGER, INTENT(OUT)                               :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                                     :: index
      COMPLEX(kind=real_4), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT), OPTIONAL                           :: block
      LOGICAL, INTENT(IN), OPTIONAL                      :: advance_iter

      INTEGER                                            :: ndims, p, p_data
      LOGICAL                                            :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      CPASSERT(buffer%data_type .EQ. dbcsr_type_complex_4)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF
      IF (p > 0) THEN
         ndata = buffer%indx(p+1, ndims+1)-buffer%indx(p, ndims+1)
      ELSE
         ndata = buffer%indx(p+1, ndims+1)
      ENDIF
      INDEX(:) = buffer%indx(p+1, 1:ndims)
      IF (PRESENT(block)) THEN
         ALLOCATE (block(ndata))
         block(:) = buffer%msg_c_sp(p_data+1:p_data+ndata)
      ENDIF

      IF (do_advance) buffer%endpos = buffer%endpos+1
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_2d_block_r_dp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_8), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_8), &
         DIMENSION(map%dims_nd(1), map%dims_nd(2)), &
         INTENT(OUT)                                     :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_3d_block_r_dp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_8), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_8), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3)), INTENT(OUT)        :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_4d_block_r_dp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_8), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_8), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4)), &
         INTENT(OUT)                                     :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_2d_block_r_sp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_4), &
         DIMENSION(map%dims_nd(1), map%dims_nd(2)), &
         INTENT(OUT)                                     :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_3d_block_r_sp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3)), INTENT(OUT)        :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_4d_block_r_sp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4)), &
         INTENT(OUT)                                     :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_2d_block_c_dp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_8), DIMENSION(:, :), INTENT(IN)  :: block_2d
      COMPLEX(kind=real_8), &
         DIMENSION(map%dims_nd(1), map%dims_nd(2)), &
         INTENT(OUT)                                     :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_3d_block_c_dp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_8), DIMENSION(:, :), INTENT(IN)  :: block_2d
      COMPLEX(kind=real_8), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3)), INTENT(OUT)    :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_4d_block_c_dp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_8), DIMENSION(:, :), INTENT(IN)  :: block_2d
      COMPLEX(kind=real_8), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3), map%dims_nd(4))&
         , INTENT(OUT)                                   :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_2d_block_c_sp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_4), DIMENSION(:, :), INTENT(IN)  :: block_2d
      COMPLEX(kind=real_4), &
         DIMENSION(map%dims_nd(1), map%dims_nd(2)), &
         INTENT(OUT)                                     :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_3d_block_c_sp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_4), DIMENSION(:, :), INTENT(IN)  :: block_2d
      COMPLEX(kind=real_4), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3)), INTENT(OUT)    :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_4d_block_c_sp(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_4), DIMENSION(:, :), INTENT(IN)  :: block_2d
      COMPLEX(kind=real_4), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3), map%dims_nd(4))&
         , INTENT(OUT)                                   :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_3d_to_2d_block_r_dp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_8), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3)), INTENT(IN)         :: block
      REAL(kind=real_8), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      INTEGER                                            :: i
      REAL(kind=real_8), ALLOCATABLE, DIMENSION(:, :, :) :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_4d_to_2d_block_r_dp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_8), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4)), &
         INTENT(IN)                                      :: block
      REAL(kind=real_8), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      INTEGER                                            :: i
      REAL(kind=real_8), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_3d_to_2d_block_r_sp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3)), INTENT(IN)         :: block
      REAL(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      INTEGER                                            :: i
      REAL(kind=real_4), ALLOCATABLE, DIMENSION(:, :, :) :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_4d_to_2d_block_r_sp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4)), &
         INTENT(IN)                                      :: block
      REAL(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      INTEGER                                            :: i
      REAL(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_3d_to_2d_block_c_dp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_8), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3)), INTENT(IN)     :: block
      COMPLEX(kind=real_8), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      COMPLEX(kind=real_8), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered
      INTEGER                                            :: i

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_4d_to_2d_block_c_dp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_8), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3), map%dims_nd(4))&
         , INTENT(IN)                                    :: block
      COMPLEX(kind=real_8), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      COMPLEX(kind=real_8), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered
      INTEGER                                            :: i

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_3d_to_2d_block_c_sp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_4), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3)), INTENT(IN)     :: block
      COMPLEX(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered
      INTEGER                                            :: i

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block ...
!> \param block_2d ...
! **************************************************************************************************
   SUBROUTINE reshape_4d_to_2d_block_c_sp(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      COMPLEX(kind=real_4), DIMENSION(map%dims_nd(1), &
         map%dims_nd(2), map%dims_nd(3), map%dims_nd(4))&
         , INTENT(IN)                                    :: block
      COMPLEX(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d

      COMPLEX(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered
      INTEGER                                            :: i

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE

END MODULE
