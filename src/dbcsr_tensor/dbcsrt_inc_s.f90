! **************************************************************************************************
!> \brief Implementation of dbcsrt_put_block for tensor rank 3.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsrt_put_3d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3)), &
         INTENT(IN)                                               :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d, dims_2d
      REAL(kind=real_4), ALLOCATABLE, DIMENSION(:, :)                       :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk

      ! reshape block
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL get_mapping_info(map_blk, dims_2d=dims_2d)
      CALL allocate_any(dims_2d, block_2d)
      CALL reshape_nd_to_2d_block(map_blk, block_2d, block)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)

      CALL dbcsr_put_block(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_put_block for tensor rank 4.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsrt_put_4d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                               :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3), sizes(4)), &
         INTENT(IN)                                               :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d, dims_2d
      REAL(kind=real_4), ALLOCATABLE, DIMENSION(:, :)                       :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk

      ! reshape block
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL get_mapping_info(map_blk, dims_2d=dims_2d)
      CALL allocate_any(dims_2d, block_2d)
      CALL reshape_nd_to_2d_block(map_blk, block_2d, block)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)

      CALL dbcsr_put_block(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_put_block for tensor rank 5.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsrt_put_5d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                               :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5)), &
         INTENT(IN)                                               :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d, dims_2d
      REAL(kind=real_4), ALLOCATABLE, DIMENSION(:, :)                       :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk

      ! reshape block
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL get_mapping_info(map_blk, dims_2d=dims_2d)
      CALL allocate_any(dims_2d, block_2d)
      CALL reshape_nd_to_2d_block(map_blk, block_2d, block)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)

      CALL dbcsr_put_block(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_put_block for tensor rank 6.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsrt_put_6d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                               :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5), sizes(6)), &
         INTENT(IN)                                               :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d, dims_2d
      REAL(kind=real_4), ALLOCATABLE, DIMENSION(:, :)                       :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk

      ! reshape block
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL get_mapping_info(map_blk, dims_2d=dims_2d)
      CALL allocate_any(dims_2d, block_2d)
      CALL reshape_nd_to_2d_block(map_blk, block_2d, block)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)

      CALL dbcsr_put_block(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d)
   END SUBROUTINE


! **************************************************************************************************
!> \brief Generic implementation of dbcsrt_put_block (arbitrary tensor rank)
!> \param tensor ...
!> \param ind block index
!> \param block block to put
! **************************************************************************************************
   SUBROUTINE dbcsrt_put_anyd_block_s(tensor, ind, block)
      TYPE(block_nd_s), INTENT(IN)                      :: block
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                               :: ind

      SELECT CASE(dbcsrt_ndims(tensor))
      CASE(3)
         CALL dbcsrt_put_3d_block_s(tensor, ind, block%sizes, block%blk)
      CASE(4)
         CALL dbcsrt_put_4d_block_s(tensor, ind, block%sizes, block%blk)
      END SELECT
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_get_block for tensor rank 3.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_3d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3)), &
         INTENT(OUT)                                              :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d
      REAL(kind=real_4), DIMENSION(:,:), POINTER                            :: block_2d_ptr
      REAL(kind=real_4), DIMENSION(:,:), ALLOCATABLE                        :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk
      LOGICAL :: tr, found

      NULLIFY (block_2d_ptr)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)
      CALL dbcsr_get_block_p(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d_ptr, tr, found)
      CPASSERT(found)
      ! convert pointer to allocatable
      CALL allocate_any(SHAPE(block_2d_ptr), block_2d, source=block_2d_ptr)

      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL reshape_2d_to_nd_block(map_blk, block_2d, block)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_get_block for tensor rank 4.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_4d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3), sizes(4)), &
         INTENT(OUT)                                              :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d
      REAL(kind=real_4), DIMENSION(:,:), POINTER                            :: block_2d_ptr
      REAL(kind=real_4), DIMENSION(:,:), ALLOCATABLE                        :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk
      LOGICAL                                                     :: tr, found

      NULLIFY (block_2d_ptr)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)
      CALL dbcsr_get_block_p(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d_ptr, tr, found)

      ! convert pointer to allocatable
      CALL allocate_any(SHAPE(block_2d_ptr), block_2d, source=block_2d_ptr)
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL reshape_2d_to_nd_block(map_blk, block_2d, block)
      DEALLOCATE(block_2d_ptr)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_get_block for tensor rank 5.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_5d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                 :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)      :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)      :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5)), &
         INTENT(OUT)                                            :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                        :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                     :: ind_2d
      REAL(kind=real_4), DIMENSION(:,:), POINTER                          :: block_2d_ptr
      REAL(kind=real_4), DIMENSION(:,:), ALLOCATABLE                      :: block_2d
      TYPE(nd_to_2d_mapping)                                    :: map_blk
      LOGICAL                                                   :: tr, found

      NULLIFY (block_2d_ptr)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)
      CALL dbcsr_get_block_p(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d_ptr, tr, found)

      ! convert pointer to allocatable
      CALL allocate_any(SHAPE(block_2d_ptr), block_2d, source=block_2d_ptr)
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL reshape_2d_to_nd_block(map_blk, block_2d, block)
      DEALLOCATE(block_2d_ptr)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation of dbcsrt_get_block for tensor rank 5.
!> \param tensor ...
!> \param ind block index
!> \param sizes block size
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_6d_block_s(tensor, ind, sizes, block)
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), INTENT(IN)        :: sizes
      REAL(kind=real_4), DIMENSION(sizes(1), sizes(2), sizes(3), sizes(4), sizes(5), sizes(6)), &
         INTENT(OUT)                                              :: block

      INTEGER, ALLOCATABLE, DIMENSION(:)                          :: map1_2d, map2_2d
      INTEGER, DIMENSION(2)                                       :: ind_2d
      REAL(kind=real_4), DIMENSION(:,:), POINTER                            :: block_2d_ptr
      REAL(kind=real_4), DIMENSION(:,:), ALLOCATABLE                        :: block_2d
      TYPE(nd_to_2d_mapping)                                      :: map_blk
      LOGICAL                                                     :: tr, found

      NULLIFY (block_2d_ptr)

      ! convert block index
      ind_2d(:) = get_2d_indices(tensor%nd_index_blk, ind)
      CALL dbcsr_get_block_p(tensor%matrix_rep, ind_2d(1), ind_2d(2), block_2d_ptr, tr, found)

      ! convert pointer to allocatable
      CALL allocate_any(SHAPE(block_2d_ptr), block_2d, source=block_2d_ptr)
      CALL get_mapping_info(tensor%nd_index_blk, map1_2d=map1_2d, map2_2d=map2_2d)
      CALL create_nd_to_2d_mapping(map_blk, sizes, map1_2d, map2_2d)
      CALL reshape_2d_to_nd_block(map_blk, block_2d, block)
      DEALLOCATE(block_2d_ptr)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Generic implementation of dbcsrt_get_block (arbitrary tensor rank)
!> \param tensor ...
!> \param ind block index
!> \param block block to get
! **************************************************************************************************
   SUBROUTINE dbcsrt_get_anyd_block_s(tensor, ind, block)
      TYPE(block_nd), INTENT(OUT)                                 :: block
      TYPE(dbcsrt_type), INTENT(INOUT)                            :: tensor
      INTEGER, DIMENSION(dbcsrt_ndims(tensor)), &
         INTENT(IN)                                               :: ind
      INTEGER, DIMENSION(dbcsrt_ndims(tensor))                    :: blk_size
      REAL(kind=real_4), DIMENSION(:), ALLOCATABLE                          :: block_arr

      CALL dbcsrt_blk_sizes(tensor, ind, blk_size)
      ALLOCATE(block_arr(PRODUCT(blk_size)))

      SELECT CASE(dbcsrt_ndims(tensor))
      CASE (3)
        CALL dbcsrt_get_3d_block_s(tensor, ind, blk_size, block_arr)
      CASE (4)
        CALL dbcsrt_get_4d_block_s(tensor, ind, blk_size, block_arr)
      END SELECT
      CALL create_block(block, blk_size, block_arr)
   END SUBROUTINE
