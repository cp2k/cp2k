!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods to operate on n-dimensional tensor blocks.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_block

#:include "dbcsr_tensor.fypp"
#:set maxdim = maxrank
#:set ndims = range(2,maxdim+1)

   USE allocate_wrap,                   ONLY: allocate_any
   USE dbcsr_api,                       ONLY: ${uselist(dtype_float_param)}$
   USE kinds,                           ONLY: ${uselist(dtype_float_prec)}$
   USE message_passing,                 ONLY: mp_environ,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_waitall
   USE dbcsr_tensor_index,              ONLY: nd_to_2d_mapping, ndims_mapping
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

#:for dparam, dtype, dsuffix in dtype_float_list
   PUBLIC :: block_nd_${dsuffix}$
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
   TYPE block_nd_${dsuffix}$
      INTEGER, DIMENSION(:), ALLOCATABLE   :: sizes
      ${dtype}$, DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

#:endfor

   TYPE block_nd
#:for dparam, dtype, dsuffix in dtype_float_list
      TYPE(block_nd_${dsuffix}$) :: ${dsuffix}$
#:endfor
      INTEGER          :: data_type = -1
   END TYPE

   TYPE block_buffer_type
      INTEGER                                    :: ndim = -1
      INTEGER                                    :: nblock = -1
      INTEGER, DIMENSION(:, :), ALLOCATABLE      :: indx
#:for dparam, dtype, dsuffix in dtype_float_list
      ${dtype}$, DIMENSION(:), ALLOCATABLE       :: msg_${dsuffix}$
#:endfor
      INTEGER                                    :: data_type = -1
      INTEGER                                    :: endpos = -1
   END TYPE

   INTERFACE create_block
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE create_block_data_${dsuffix}$
#:endfor
      MODULE PROCEDURE create_block_nodata
   END INTERFACE

   INTERFACE block_buffer_add_block
#:for dparam, dtype, dsuffix in dtype_float_list
      MODULE PROCEDURE block_buffer_add_block_${dsuffix}$
#:endfor
   END INTERFACE

   INTERFACE reshape_nd_to_2d_block
#:for dparam, dtype, dsuffix in dtype_float_list
#:for ndim in ndims
      MODULE PROCEDURE reshape_${ndim}$d_to_2d_block_${dsuffix}$
#:endfor
#:endfor
   END INTERFACE reshape_nd_to_2d_block

   INTERFACE reshape_2d_to_nd_block
#:for dparam, dtype, dsuffix in dtype_float_list
#:for ndim in ndims
      MODULE PROCEDURE reshape_2d_to_${ndim}$d_block_${dsuffix}$
#:endfor
#:endfor
   END INTERFACE reshape_2d_to_nd_block

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param block ...
!> \param sizes ...
!> \param data_type ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata(block, sizes, data_type)
      TYPE(block_nd), INTENT(OUT)       :: block
      INTEGER, DIMENSION(:), INTENT(IN) :: sizes
      INTEGER, INTENT(IN)               :: data_type

      block%data_type = data_type
      SELECT CASE (data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         CALL create_block_nodata_${dsuffix}$(block%${dsuffix}$, sizes)
#:endfor
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block(block)
      TYPE(block_nd), INTENT(INOUT) :: block

      SELECT CASE (block%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         CALL destroy_block_${dsuffix}$(block%${dsuffix}$)
#:endfor
      END SELECT

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block ...
!> \retval blk_size ...
! **************************************************************************************************
   PURE FUNCTION block_size(block)
      TYPE(block_nd), INTENT(IN)         :: block
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_size

      SELECT CASE (block%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         CALL allocate_any(SHAPE(block%${dsuffix}$%sizes), block_size, source=block%${dsuffix}$%sizes)
#:endfor
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
      TYPE(block_buffer_type), INTENT(OUT) :: buffer
      INTEGER, INTENT(IN)                  :: nblock, ndata, data_type, ndim

      buffer%nblock = nblock
      buffer%data_type = data_type
      buffer%endpos = 0
      buffer%ndim = ndim
      SELECT CASE (data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         ALLOCATE (buffer%msg_${dsuffix}$(ndata))
#:endfor
      END SELECT
      ALLOCATE (buffer%indx(nblock, ndim+1))
   END SUBROUTINE block_buffer_create

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE block_buffer_destroy(buffer)
      TYPE(block_buffer_type), INTENT(INOUT) :: buffer

      SELECT CASE (buffer%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         DEALLOCATE (buffer%msg_${dsuffix}$)
#:endfor
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
      TYPE(block_buffer_type), INTENT(IN) :: buffer
      INTEGER                             :: ndims_buffer

      ndims_buffer = buffer%ndim
   END FUNCTION

! **************************************************************************************************
!> \brief insert a block into block buffer (at current iterator position)
!> \param buffer ...
!> \param index index of block
!> \param block block
! **************************************************************************************************
   SUBROUTINE block_buffer_add_anyd_block(buffer, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(IN)                               :: index
      TYPE(block_nd), INTENT(IN)                  :: block

      SELECT CASE (block%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         CALL block_buffer_add_block_${dsuffix}$(buffer, SIZE(block%${dsuffix}$%blk), index, block%${dsuffix}$%blk)
#:endfor
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
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(OUT)                        :: ndata
      INTEGER, DIMENSION(ndims_buffer(buffer)), &
         INTENT(OUT)                              :: index
      TYPE(block_nd), INTENT(INOUT), OPTIONAL     :: block
      LOGICAL, INTENT(IN), OPTIONAL               :: advance_iter

      SELECT CASE (buffer%data_type)
#:for dparam, dtype, dsuffix in dtype_float_list
      CASE (${dparam}$)
         IF (PRESENT(block)) THEN
            CALL block_buffer_get_next_block_${dsuffix}$(buffer, ndata, index, block%${dsuffix}$%blk, advance_iter=advance_iter)
         ELSE
            CALL block_buffer_get_next_block_${dsuffix}$(buffer, ndata, index, advance_iter=advance_iter)
         ENDIF
#:endfor
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
! **************************************************************************************************
   SUBROUTINE block_buffer_iterator_reset(buffer)
      TYPE(block_buffer_type), INTENT(INOUT) :: buffer

      buffer%endpos = 0
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param buffer ...
!> \retval block_buffer_blocks_left ...
! **************************************************************************************************
   PURE FUNCTION block_buffer_blocks_left(buffer)
      TYPE(block_buffer_type), INTENT(IN) :: buffer
      LOGICAL                             :: block_buffer_blocks_left

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
      INTEGER, INTENT(IN)                    :: mp_comm
      TYPE(block_buffer_type), DIMENSION(0:) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(:, :)               :: req_array

      INTEGER                                :: iproc, mynode, numnodes, rec_counter, &
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
#:for dparam, dtype, dsuffix in dtype_float_list
               CASE (${dparam}$)
                  CALL mp_irecv(buffer_recv(iproc)%msg_${dsuffix}$, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
#:endfor
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
#:for dparam, dtype, dsuffix in dtype_float_list
               CASE (${dparam}$)
                  CALL mp_isend(buffer_send(iproc)%msg_${dsuffix}$, iproc, mp_comm, req_array(send_counter, 2), tag=7)
#:endfor
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
#:for dparam, dtype, dsuffix in dtype_float_list
            CASE (${dparam}$)
               buffer_recv(0)%msg_${dsuffix}$(:) = buffer_send(0)%msg_${dsuffix}$(:)
#:endfor
            END SELECT
         ENDIF
      ENDIF

   END SUBROUTINE

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief Create block from array, array can be n-dimensional.
!> \param block ...
!> \param sizes ...
!> \param array ...
! **************************************************************************************************
   SUBROUTINE create_block_data_${dsuffix}$(block, sizes, array)
      TYPE(block_nd), INTENT(OUT)                       :: block
      INTEGER, DIMENSION(:), INTENT(IN)                 :: sizes
      ${dtype}$, DIMENSION(PRODUCT(sizes)),  INTENT(IN) :: array

      ASSOCIATE(blk => block%${dsuffix}$)
      block%data_type = ${dparam}$
      ALLOCATE(blk%sizes(SIZE(sizes)))
      blk%sizes(:) = sizes(:)
      ALLOCATE(blk%blk(PRODUCT(sizes)))
      blk%blk(:) = array(:)
      END ASSOCIATE
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief Create and allocate block, but no data.
!> \param block ...
!> \param sizes ...
! **************************************************************************************************
   SUBROUTINE create_block_nodata_${dsuffix}$(block, sizes)
      INTEGER, INTENT(IN), DIMENSION(:)       :: sizes
      TYPE(block_nd_${dsuffix}$), INTENT(OUT) :: block
      ALLOCATE(block%sizes(SIZE(sizes)))
      block%sizes(:) = sizes(:)
      ALLOCATE(block%blk(PRODUCT(sizes)))
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE destroy_block_${dsuffix}$ (block)
      TYPE(block_nd_${dsuffix}$), INTENT(INOUT) :: block
      DEALLOCATE(block%blk)
      DEALLOCATE(block%sizes)
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief add block to buffer.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE block_buffer_add_block_${dsuffix}$(buffer, ndata, index, block)
      TYPE(block_buffer_type), INTENT(INOUT)               :: buffer
      INTEGER, INTENT(IN)                                  :: ndata
      ${dtype}$, DIMENSION(ndata), INTENT(IN)              :: block
      INTEGER, DIMENSION(ndims_buffer(buffer)), INTENT(IN) :: index
      INTEGER                                              :: p, ndims, p_data
      CPASSERT(buffer%data_type .EQ. ${dparam}$)
      ndims = ndims_buffer(buffer)
      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = buffer%indx(p, ndims+1)
      ENDIF

      buffer%msg_${dsuffix}$(p_data+1:p_data+ndata) = block(:)
      buffer%indx(p+1, 1:ndims) = index(:)
      IF (p > 0) THEN
         buffer%indx(p+1,ndims+1) = buffer%indx(p,ndims+1)+ndata
      ELSE
         buffer%indx(p+1, ndims+1) = ndata
      ENDIF
      buffer%endpos = buffer%endpos+1
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
! **************************************************************************************************
!> \brief get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
!> \param buffer ...
!> \param ndata ...
!> \param index ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE block_buffer_get_next_block_${dsuffix}$(buffer, ndata, index, block, advance_iter)
      TYPE(block_buffer_type), INTENT(INOUT)                      :: buffer
      INTEGER, INTENT(OUT)                                        :: ndata
      ${dtype}$, DIMENSION(:), ALLOCATABLE, OPTIONAL, INTENT(OUT) :: block
      INTEGER, DIMENSION(ndims_buffer(buffer)), INTENT(OUT)       :: index
      INTEGER                                                     :: p, ndims, p_data
      LOGICAL, INTENT(IN), OPTIONAL                               :: advance_iter
      LOGICAL                                                     :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      CPASSERT(buffer%data_type .EQ. ${dparam}$)
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
      index(:) = buffer%indx(p+1,1:ndims)
      IF (PRESENT(block)) THEN
         ALLOCATE (block(ndata))
         block(:) = buffer%msg_${dsuffix}$(p_data+1:p_data+ndata)
      ENDIF

      IF(do_advance) buffer%endpos = buffer%endpos+1
   END SUBROUTINE
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
#:for ndim in ndims
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_${ndim}$d_block_${dsuffix}$(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                        :: map
      ${dtype}$, DIMENSION(:, :), INTENT(IN)                    :: block_2d
      ${dtype}$, DIMENSION(${arrlist("map%dims_nd", nmax=ndim)}$), &
         INTENT(OUT)                                            :: block
      INTEGER, DIMENSION(ndims_mapping(map))                    :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(${shape_colon(ndim)}$) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE
#:endfor
#:endfor

#:for dparam, dtype, dsuffix in dtype_float_list
#:for ndim in ndims
#! need to exclude ndim == 2 because reshape_2d_to_2d_block has already been defined
#:if ndim > 2
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_${ndim}$d_to_2d_block_${dsuffix}$(map, block, block_2d)
      TYPE(nd_to_2d_mapping), INTENT(IN)                                 :: map
      ${dtype}$, DIMENSION(${arrlist("map%dims_nd", nmax=ndim)}$), INTENT(IN) :: block
      ${dtype}$, DIMENSION(map%dims_2d(1), map%dims_2d(2)), INTENT(OUT)  :: block_2d

      INTEGER                                                            :: i
      INTEGER, DIMENSION(ndims_mapping(map))                             :: shape_in, shape_reordered
      ${dtype}$, ALLOCATABLE, DIMENSION(${shape_colon(ndim)}$)           :: block_tmp

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(${shape_colon(ndim)}$) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
#:endif
#:endfor
#:endfor

END MODULE
