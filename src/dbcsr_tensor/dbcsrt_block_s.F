!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods to operate on n-dimensional tensor blocks.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsrt_block_s
   USE allocate_wrap,                   ONLY: allocate_any
   USE dbcsrt_index_types,              ONLY: nd_to_2d_mapping,&
                                              ndims_mapping
   USE kinds,                           ONLY: real_4
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE block_nd_s
      INTEGER, DIMENSION(:), ALLOCATABLE :: sizes
      REAL(kind=real_4), DIMENSION(:), ALLOCATABLE :: blk
   END TYPE

   INTERFACE reshape_nd_to_2d_block
      MODULE PROCEDURE reshape_3d_to_2d_block_s
      MODULE PROCEDURE reshape_4d_to_2d_block_s
      MODULE PROCEDURE reshape_5d_to_2d_block_s
      MODULE PROCEDURE reshape_6d_to_2d_block_s
   END INTERFACE reshape_nd_to_2d_block

   INTERFACE reshape_2d_to_nd_block
      MODULE PROCEDURE reshape_2d_to_3d_block_s
      MODULE PROCEDURE reshape_2d_to_4d_block_s
      MODULE PROCEDURE reshape_2d_to_5d_block_s
      MODULE PROCEDURE reshape_2d_to_6d_block_s
   END INTERFACE reshape_2d_to_nd_block

   PUBLIC :: &
      block_nd_s, &
      reshape_2d_to_nd_block, &
      reshape_nd_to_2d_block

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_3d_block_s(map, block_2d, block)
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
   SUBROUTINE reshape_3d_to_2d_block_s(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3)), INTENT(IN)         :: block

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
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_4d_block_s(map, block_2d, block)
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
   SUBROUTINE reshape_4d_to_2d_block_s(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4)), &
         INTENT(IN)                                      :: block

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
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_5d_block_s(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4), map&
         %dims_nd(5)), INTENT(OUT)                       :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_5d_to_2d_block_s(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4), map&
         %dims_nd(5)), INTENT(IN)                        :: block

      INTEGER                                            :: i
      REAL(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_2d_to_6d_block_s(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), DIMENSION(:, :), INTENT(IN)     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4), map&
         %dims_nd(5), map%dims_nd(6)), INTENT(OUT)       :: block

      INTEGER, DIMENSION(ndims_mapping(map))             :: map_tmp

      map_tmp(:) = [map%map1_2d, map%map2_2d]
      block(:, :, :, :, :, :) = RESHAPE(block_2d, SHAPE=SHAPE(block), ORDER=map_tmp)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \param block_2d ...
!> \param block ...
! **************************************************************************************************
   SUBROUTINE reshape_6d_to_2d_block_s(map, block_2d, block)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      REAL(kind=real_4), &
         DIMENSION(map%dims_2d(1), map%dims_2d(2)), &
         INTENT(OUT)                                     :: block_2d
      REAL(kind=real_4), DIMENSION(map%dims_nd(1), map%&
         dims_nd(2), map%dims_nd(3), map%dims_nd(4), map&
         %dims_nd(5), map%dims_nd(6)), INTENT(IN)        :: block

      INTEGER                                            :: i
      REAL(kind=real_4), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :, :)                     :: block_tmp
      INTEGER, DIMENSION(ndims_mapping(map))             :: shape_in, shape_reordered

      shape_in = SHAPE(block)
      DO i = 1, SIZE(shape_in)
         shape_reordered(map%map_nd(i)) = shape_in(i)
      ENDDO
      CALL allocate_any(shape_reordered, block_tmp)
      block_tmp(:, :, :, :, :, :) = RESHAPE(block, SHAPE=SHAPE(block_tmp), ORDER=map%map_nd)
      block_2d(:, :) = RESHAPE(block_tmp, SHAPE=SHAPE(block_2d))
   END SUBROUTINE
END MODULE
