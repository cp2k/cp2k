!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief tensor index and mapping to DBCSR index.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_index
   USE allocate_wrap,                   ONLY: allocate_any
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_index'

   PUBLIC :: &
      combine_index, &
      create_nd_to_2d_mapping, &
      destroy_nd_to_2d_mapping, &
      get_2d_indices, &
      get_mapping_info, &
      get_nd_indices, &
      nd_to_2d_mapping, &
      ndims_mapping, &
      split_index

   TYPE nd_to_2d_mapping
      INTEGER                                      :: ndim_nd
      INTEGER                                      :: ndim1_2d
      INTEGER                                      :: ndim2_2d

      INTEGER, DIMENSION(:), ALLOCATABLE           :: dims_nd
      INTEGER, DIMENSION(2)                        :: dims_2d
      INTEGER, DIMENSION(:), ALLOCATABLE           :: dims1_2d
      INTEGER, DIMENSION(:), ALLOCATABLE           :: dims2_2d

      INTEGER, DIMENSION(:), ALLOCATABLE           :: map1_2d
      INTEGER, DIMENSION(:), ALLOCATABLE           :: map2_2d
      INTEGER, DIMENSION(:), ALLOCATABLE           :: map_nd

      INTEGER                                      :: base
      LOGICAL                                      :: col_major
   END TYPE nd_to_2d_mapping

CONTAINS

! **************************************************************************************************
!> \brief Create all data needed to quickly map between nd index and 2d index.
!> \param map index mapping data
!> \param dims nd sizes
!> \param map1_2d which nd-indices map to first matrix index and in which order
!> \param map2_2d which nd-indices map to second matrix index and in which order
!> \param base base index (1 for Fortran-style, 0 for C-style, default is 1)
!> \param col_major whether index should be column major order (.TRUE. for Fortran-style,
!>        .FALSE. for C-style, default is .TRUE.).
! **************************************************************************************************
   SUBROUTINE create_nd_to_2d_mapping(map, dims, map1_2d, map2_2d, base, col_major)
      TYPE(nd_to_2d_mapping), INTENT(OUT)                :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dims, map1_2d, map2_2d
      INTEGER, INTENT(IN), OPTIONAL                      :: base
      LOGICAL, INTENT(IN), OPTIONAL                      :: col_major

      INTEGER                                            :: i

      IF (PRESENT(col_major)) THEN
         map%col_major = col_major
      ELSE
         map%col_major = .TRUE.
      ENDIF

      IF (PRESENT(base)) THEN
         map%base = base
      ELSE
         map%base = 1
      ENDIF

      map%ndim1_2d = SIZE(map1_2d)
      map%ndim2_2d = SIZE(map2_2d)
      map%ndim_nd = SIZE(dims)

      ALLOCATE (map%map1_2d(map%ndim1_2d))
      ALLOCATE (map%map2_2d(map%ndim2_2d))
      ALLOCATE (map%dims1_2d(map%ndim1_2d))
      ALLOCATE (map%dims2_2d(map%ndim2_2d))
      ALLOCATE (map%dims_nd(map%ndim1_2d+map%ndim2_2d))
      ALLOCATE (map%map_nd(map%ndim_nd))

      map%dims_nd(:) = dims(:)
      map%map1_2d(:) = map1_2d(:)
      map%map2_2d(:) = map2_2d(:)

      DO i = 1, SIZE(map1_2d)
         map%map_nd(map1_2d(i)) = i
         map%dims1_2d(i) = dims(map1_2d(i))
      ENDDO

      DO i = 1, SIZE(map2_2d)
         map%map_nd(map2_2d(i)) = i+SIZE(map1_2d)
         map%dims2_2d(i) = dims(map2_2d(i))
      ENDDO

      map%dims_2d = [PRODUCT(map%dims1_2d), PRODUCT(map%dims2_2d)]

   END SUBROUTINE create_nd_to_2d_mapping

! **************************************************************************************************
!> \brief ...
!> \param map ...
! **************************************************************************************************
   SUBROUTINE destroy_nd_to_2d_mapping(map)
      TYPE(nd_to_2d_mapping), INTENT(INOUT)              :: map

      DEALLOCATE (map%dims1_2d)
      DEALLOCATE (map%dims2_2d)
      DEALLOCATE (map%map1_2d)
      DEALLOCATE (map%map2_2d)
      DEALLOCATE (map%map_nd)
      DEALLOCATE (map%dims_nd)
   END SUBROUTINE destroy_nd_to_2d_mapping

! **************************************************************************************************
!> \brief ...
!> \param map ...
!> \retval ndims_mapping ...
! **************************************************************************************************
   PURE FUNCTION ndims_mapping(map)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER                                            :: ndims_mapping

      ndims_mapping = map%ndim_nd
   END FUNCTION

! **************************************************************************************************
!> \brief get mapping info
!> \param map index mapping data.
!> \param ndim_nd number of dimensions
!> \param ndim1_2d number of dimensions that map to first 2d index
!> \param ndim2_2d number of dimensions that map to first 2d index
!> \param dims_2d 2d dimensions
!> \param dims_nd nd dimensions
!> \param dims1_2d dimensions that map to first 2d index
!> \param dims2_2d dimensions that map to second 2d index
!> \param map1_2d indices that map to first 2d index
!> \param map2_2d indices that map to second 2d index
!> \param map_nd inverse of [map1_2d, map2_2d]
!> \param base base index
!> \param col_major is index in column major order
! **************************************************************************************************
   SUBROUTINE get_mapping_info(map, ndim_nd, ndim1_2d, ndim2_2d, dims_2d, dims_nd, dims1_2d, dims2_2d, &
                               map1_2d, map2_2d, map_nd, base, col_major)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, INTENT(OUT), OPTIONAL                     :: ndim_nd, ndim1_2d, ndim2_2d
      INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL       :: dims_2d
      INTEGER, DIMENSION(ndims_mapping(map)), &
         INTENT(OUT), OPTIONAL                           :: dims_nd
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: dims1_2d, dims2_2d, map1_2d, map2_2d
      INTEGER, DIMENSION(ndims_mapping(map)), &
         INTENT(OUT), OPTIONAL                           :: map_nd
      INTEGER, INTENT(OUT), OPTIONAL                     :: base
      LOGICAL, INTENT(OUT), OPTIONAL                     :: col_major

      IF (PRESENT(ndim_nd)) ndim_nd = map%ndim_nd
      IF (PRESENT(ndim1_2d)) ndim1_2d = map%ndim1_2d
      IF (PRESENT(ndim2_2d)) ndim2_2d = map%ndim2_2d
      IF (PRESENT(dims_2d)) dims_2d(:) = map%dims_2d(:)
      IF (PRESENT(dims_nd)) THEN
         dims_nd(:) = map%dims_nd(:)
      ENDIF
      IF (PRESENT(dims1_2d)) THEN
         CALL allocate_any(SHAPE(map%dims1_2d), dims1_2d, source=map%dims1_2d)
      ENDIF
      IF (PRESENT(dims2_2d)) THEN
         CALL allocate_any(SHAPE(map%dims2_2d), dims2_2d, source=map%dims2_2d)
      ENDIF
      IF (PRESENT(map1_2d)) THEN
         CALL allocate_any(SHAPE(map%map1_2d), map1_2d, source=map%map1_2d)
      ENDIF
      IF (PRESENT(map2_2d)) THEN
         CALL allocate_any(SHAPE(map%map2_2d), map2_2d, map%map2_2d)
      ENDIF
      IF (PRESENT(map_nd)) THEN
         map_nd(:) = map%map_nd(:)
      ENDIF
      IF (PRESENT(base)) THEN
         base = map%base
      ENDIF
      IF (PRESENT(col_major)) THEN
         col_major = map%col_major
      ENDIF

   END SUBROUTINE get_mapping_info

! **************************************************************************************************
!> \brief transform nd index to flat index
!> \param ind_in nd index
!> \param dims nd dimensions
!> \param base base index, default 1 (Fortran style)
!> \param col_major column major ordering, default .TRUE. (Fortran style)
!> \retval ind_out flat index
! **************************************************************************************************
   FUNCTION combine_index(ind_in, dims, base, col_major) RESULT(ind_out)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ind_in, dims
      INTEGER, INTENT(IN), OPTIONAL                      :: base
      LOGICAL, INTENT(IN), OPTIONAL                      :: col_major
      INTEGER                                            :: ind_out

      INTEGER                                            :: i_dim, my_base
      LOGICAL                                            :: my_col_major

      CPASSERT(SIZE(ind_in) .EQ. SIZE(dims))

      IF (PRESENT(base)) THEN
         my_base = base
      ELSE
         my_base = 1
      ENDIF

      IF (PRESENT(col_major)) THEN
         my_col_major = col_major
      ELSE
         my_col_major = .TRUE.
      ENDIF

      IF (my_col_major) THEN
         i_dim = SIZE(dims)-1
         ind_out = ind_in(i_dim+1)
         DO WHILE (i_dim .GE. 1)
            ind_out = (ind_out-my_base)*dims(i_dim)+(ind_in(i_dim)-my_base)+my_base
            i_dim = i_dim-1
         ENDDO
      ELSE
         i_dim = 2
         ind_out = ind_in(i_dim-1)
         DO WHILE (i_dim .LE. SIZE(dims))
            ind_out = (ind_out-my_base)*dims(i_dim)+(ind_in(i_dim)-my_base)+my_base
            i_dim = i_dim+1
         ENDDO
      ENDIF
   END FUNCTION combine_index

! **************************************************************************************************
!> \brief transform flat index to nd index
!> \param ind_in flat index
!> \param dims nd dimensions
!> \param base base index, default 1 (Fortran style)
!> \param col_major column major ordering, default .TRUE. (Fortran style)
!> \retval ind_out nd index
! **************************************************************************************************
   FUNCTION split_index(ind_in, dims, base, col_major) RESULT(ind_out)
      INTEGER, INTENT(IN)                                :: ind_in
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dims
      INTEGER, INTENT(IN), OPTIONAL                      :: base
      LOGICAL, INTENT(IN), OPTIONAL                      :: col_major
      INTEGER, DIMENSION(SIZE(dims))                     :: ind_out

      INTEGER                                            :: i_dim, my_base, tmp
      LOGICAL                                            :: my_col_major

      IF (PRESENT(base)) THEN
         my_base = base
      ELSE
         my_base = 1
      ENDIF

      IF (PRESENT(col_major)) THEN
         my_col_major = col_major
      ELSE
         my_col_major = .TRUE.
      ENDIF

      IF (my_col_major) THEN
         i_dim = 1
         tmp = ind_in
         DO WHILE (i_dim .LE. SIZE(dims))
            ind_out(i_dim) = MOD(tmp-my_base, dims(i_dim))+my_base
            tmp = (tmp-my_base)/dims(i_dim)+my_base
            i_dim = i_dim+1
         END DO
      ELSE
         i_dim = SIZE(dims)
         tmp = ind_in
         DO WHILE (i_dim .GE. 1)
            ind_out(i_dim) = MOD(tmp-my_base, dims(i_dim))+my_base
            tmp = (tmp-my_base)/dims(i_dim)+my_base
            i_dim = i_dim-1
         END DO
      ENDIF
   END FUNCTION split_index

! **************************************************************************************************
!> \brief transform nd index to 2d index, using info from index mapping.
!> \param map index mapping
!> \param ind_in nd index
!> \retval ind_out 2d index
! **************************************************************************************************
   FUNCTION get_2d_indices(map, ind_in) RESULT(ind_out)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, DIMENSION(ndims_mapping(map)), INTENT(IN) :: ind_in
      INTEGER, DIMENSION(2)                              :: ind_out

      INTEGER                                            :: i
      INTEGER, DIMENSION(map%ndim2_2d)                   :: ind_2
      INTEGER, DIMENSION(map%ndim1_2d)                   :: ind_1

      DO i = 1, map%ndim1_2d
         ind_1(i) = ind_in(map%map1_2d(i))
      ENDDO
      DO i = 1, map%ndim2_2d
         ind_2(i) = ind_in(map%map2_2d(i))
      ENDDO
      ind_out(1) = combine_index(ind_1, map%dims1_2d, base=map%base, col_major=map%col_major)
      ind_out(2) = combine_index(ind_2, map%dims2_2d, base=map%base, col_major=map%col_major)
   END FUNCTION get_2d_indices

! **************************************************************************************************
!> \brief transform 2d index to nd index, using info from index mapping.
!> \param map index mapping
!> \param ind_in 2d index
!> \retval ind_out nd index
! **************************************************************************************************
   FUNCTION get_nd_indices(map, ind_in) RESULT(ind_out)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, DIMENSION(2), INTENT(IN)                  :: ind_in
      INTEGER, DIMENSION(ndims_mapping(map))             :: ind_out

      INTEGER                                            :: i
      INTEGER, DIMENSION(map%ndim2_2d)                   :: ind_2
      INTEGER, DIMENSION(map%ndim1_2d)                   :: ind_1

      ind_1(:) = split_index(ind_in(1), map%dims1_2d, base=map%base, col_major=map%col_major)
      ind_2(:) = split_index(ind_in(2), map%dims2_2d, base=map%base, col_major=map%col_major)
      DO i = 1, map%ndim1_2d
         ind_out(map%map1_2d(i)) = ind_1(i)
      ENDDO
      DO i = 1, map%ndim2_2d
         ind_out(map%map2_2d(i)) = ind_2(i)
      ENDDO
   END FUNCTION get_nd_indices

END MODULE
