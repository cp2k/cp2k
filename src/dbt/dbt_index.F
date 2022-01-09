!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief tensor index and mapping to DBM index
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_index
   USE dbt_allocate_wrap, ONLY: allocate_any
   USE kinds, ONLY: int_8
#include "../base/base_uses.f90"
   #:include "dbt_macros.fypp"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_index'

   PUBLIC :: &
      combine_tensor_index, &
      combine_pgrid_index, &
      create_nd_to_2d_mapping, &
      destroy_nd_to_2d_mapping, &
      get_2d_indices_tensor, &
      get_2d_indices_pgrid, &
      dbt_get_mapping_info, &
      get_nd_indices_tensor, &
      get_nd_indices_pgrid, &
      nd_to_2d_mapping, &
      ndims_mapping, &
      split_tensor_index, &
      split_pgrid_index, &
      ndims_mapping_row, &
      ndims_mapping_column, &
      dbt_inverse_order, &
      permute_index

   TYPE nd_to_2d_mapping
      INTEGER                                      :: ndim_nd
      INTEGER                                      :: ndim1_2d
      INTEGER                                      :: ndim2_2d

      INTEGER, DIMENSION(:), ALLOCATABLE           :: dims_nd
      INTEGER(KIND=int_8), DIMENSION(2)            :: dims_2d
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
!> \param col_major whether index should be column major order
!>                  (.TRUE. for Fortran-style, .FALSE. for C-style, default is .TRUE.).
!> \author Patrick Seewald
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
      END IF

      IF (PRESENT(base)) THEN
         map%base = base
      ELSE
         map%base = 1
      END IF

      map%ndim1_2d = SIZE(map1_2d)
      map%ndim2_2d = SIZE(map2_2d)
      map%ndim_nd = SIZE(dims)

      ALLOCATE (map%map1_2d, source=map1_2d)
      ALLOCATE (map%map2_2d, source=map2_2d)
      ALLOCATE (map%dims_nd, source=dims)
      ALLOCATE (map%dims1_2d, source=dims(map1_2d))
      ALLOCATE (map%dims2_2d, source=dims(map2_2d))

      ALLOCATE (map%map_nd(map%ndim_nd))
      map%map_nd(map1_2d) = (/(i, i=1, SIZE(map1_2d))/)
      map%map_nd(map2_2d) = (/(i + SIZE(map1_2d), i=1, SIZE(map2_2d))/)

      map%dims_2d = [PRODUCT(INT(map%dims1_2d, KIND=int_8)), PRODUCT(INT(map%dims2_2d, KIND=int_8))]

   END SUBROUTINE create_nd_to_2d_mapping

! **************************************************************************************************
!> \brief
!> \author Patrick Seewald
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
!> \brief
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION ndims_mapping(map)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER                                            :: ndims_mapping

      ndims_mapping = map%ndim_nd
   END FUNCTION

! **************************************************************************************************
!> \brief how many tensor dimensions are mapped to matrix row
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION ndims_mapping_row(map)
      TYPE(nd_to_2d_mapping), INTENT(IN) :: map
      INTEGER :: ndims_mapping_row
      ndims_mapping_row = map%ndim1_2d
   END FUNCTION

! **************************************************************************************************
!> \brief how many tensor dimensions are mapped to matrix column
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION ndims_mapping_column(map)
      TYPE(nd_to_2d_mapping), INTENT(IN) :: map
      INTEGER :: ndims_mapping_column
      ndims_mapping_column = map%ndim2_2d
   END FUNCTION

! **************************************************************************************************
!> \brief get mapping info
!> \param map index mapping data
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
!> \author Patrick Seewald
! **************************************************************************************************
   PURE SUBROUTINE dbt_get_mapping_info(map, ndim_nd, ndim1_2d, ndim2_2d, dims_2d_i8, &
                                        dims_2d, dims_nd, dims1_2d, dims2_2d, &
                                        map1_2d, map2_2d, map_nd, base, col_major)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, INTENT(OUT), OPTIONAL                     :: ndim_nd, ndim1_2d, ndim2_2d
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(OUT), OPTIONAL       :: dims_2d_i8
      INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL :: dims_2d
      INTEGER, DIMENSION(ndims_mapping(map)), &
         INTENT(OUT), OPTIONAL                           :: dims_nd
      INTEGER, DIMENSION(ndims_mapping_row(map)), INTENT(OUT), &
         OPTIONAL                                        :: dims1_2d
      INTEGER, DIMENSION(ndims_mapping_column(map)), INTENT(OUT), &
         OPTIONAL                                        :: dims2_2d
      INTEGER, DIMENSION(ndims_mapping_row(map)), INTENT(OUT), &
         OPTIONAL                                        :: map1_2d
      INTEGER, DIMENSION(ndims_mapping_column(map)), INTENT(OUT), &
         OPTIONAL                                        :: map2_2d
      INTEGER, DIMENSION(ndims_mapping(map)), &
         INTENT(OUT), OPTIONAL                           :: map_nd
      INTEGER, INTENT(OUT), OPTIONAL                     :: base
      LOGICAL, INTENT(OUT), OPTIONAL                     :: col_major

      IF (PRESENT(ndim_nd)) ndim_nd = map%ndim_nd
      IF (PRESENT(ndim1_2d)) ndim1_2d = map%ndim1_2d
      IF (PRESENT(ndim2_2d)) ndim2_2d = map%ndim2_2d
      IF (PRESENT(dims_2d_i8)) dims_2d_i8(:) = map%dims_2d(:)
      IF (PRESENT(dims_2d)) dims_2d(:) = INT(map%dims_2d(:))
      IF (PRESENT(dims_nd)) THEN
         dims_nd(:) = map%dims_nd(:)
      END IF
      IF (PRESENT(dims1_2d)) THEN
         dims1_2d(:) = map%dims1_2d
      END IF
      IF (PRESENT(dims2_2d)) THEN
         dims2_2d(:) = map%dims2_2d
      END IF
      IF (PRESENT(map1_2d)) THEN
         map1_2d(:) = map%map1_2d
      END IF
      IF (PRESENT(map2_2d)) THEN
         map2_2d(:) = map%map2_2d
      END IF
      IF (PRESENT(map_nd)) THEN
         map_nd(:) = map%map_nd(:)
      END IF
      IF (PRESENT(base)) THEN
         base = map%base
      END IF
      IF (PRESENT(col_major)) THEN
         col_major = map%col_major
      END IF

   END SUBROUTINE dbt_get_mapping_info

! **************************************************************************************************
!> \brief transform nd index to flat index
!> \param ind_in nd index
!> \param dims nd dimensions
!> \param ind_out flat index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION combine_tensor_index(ind_in, dims) RESULT(ind_out)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ind_in, dims
      INTEGER(KIND=int_8)                                :: ind_out
      INTEGER                                            :: i_dim

      ind_out = ind_in(SIZE(dims))
      DO i_dim = SIZE(dims) - 1, 1, -1
         ind_out = (ind_out - 1)*dims(i_dim) + ind_in(i_dim)
      END DO

   END FUNCTION

! **************************************************************************************************
!> \brief transform nd index to flat index
!> \param ind_in nd index
!> \param dims nd dimensions
!> \param ind_out flat index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION combine_pgrid_index(ind_in, dims) RESULT(ind_out)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ind_in, dims
      INTEGER                                            :: ind_out

      INTEGER                                            :: i_dim

      ind_out = ind_in(1)
      DO i_dim = 2, SIZE(dims)
         ind_out = ind_out*dims(i_dim) + ind_in(i_dim)
      END DO
   END FUNCTION

! **************************************************************************************************
!> \brief transform flat index to nd index
!> \param ind_in flat index
!> \param dims nd dimensions
!> \param ind_out nd index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION split_tensor_index(ind_in, dims) RESULT(ind_out)
      INTEGER(KIND=int_8), INTENT(IN)                    :: ind_in
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dims
      INTEGER, DIMENSION(SIZE(dims))                     :: ind_out

      INTEGER(KIND=int_8)                                :: tmp
      INTEGER                                            :: i_dim

      tmp = ind_in
      DO i_dim = 1, SIZE(dims)
         ind_out(i_dim) = INT(MOD(tmp - 1, INT(dims(i_dim), int_8)) + 1)
         tmp = (tmp - 1)/dims(i_dim) + 1
      END DO

   END FUNCTION

! **************************************************************************************************
!> \brief transform flat index to nd index
!> \param ind_in flat index
!> \param dims nd dimensions
!> \param ind_out nd index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION split_pgrid_index(ind_in, dims) RESULT(ind_out)
      INTEGER, INTENT(IN)                                :: ind_in
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dims
      INTEGER, DIMENSION(SIZE(dims))                     :: ind_out

      INTEGER                                            :: tmp
      INTEGER                                            :: i_dim

      tmp = ind_in
      DO i_dim = SIZE(dims), 1, -1
         ind_out(i_dim) = MOD(tmp, dims(i_dim))
         tmp = tmp/dims(i_dim)
      END DO
   END FUNCTION

! **************************************************************************************************
!> \brief transform nd index to 2d index, using info from index mapping.
!> \param map index mapping
!> \param ind_in nd index
!> \param ind_out 2d index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION get_2d_indices_tensor(map, ind_in) RESULT(ind_out)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, DIMENSION(map%ndim_nd), INTENT(IN) :: ind_in
      INTEGER(KIND=int_8), DIMENSION(2)                  :: ind_out
      INTEGER :: i
      INTEGER, DIMENSION(${maxrank}$)                    :: ind_tmp

      DO i = 1, map%ndim1_2d
         ind_tmp(i) = ind_in(map%map1_2d(i))
      END DO
      ind_out(1) = combine_tensor_index(ind_tmp(:map%ndim1_2d), map%dims1_2d)

      DO i = 1, map%ndim2_2d
         ind_tmp(i) = ind_in(map%map2_2d(i))
      END DO
      ind_out(2) = combine_tensor_index(ind_tmp(:map%ndim2_2d), map%dims2_2d)
   END FUNCTION

! **************************************************************************************************
!> \brief transform nd index to 2d index, using info from index mapping.
!> \param map index mapping
!> \param ind_in nd index
!> \param ind_out 2d index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION get_2d_indices_pgrid(map, ind_in) RESULT(ind_out)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, DIMENSION(map%ndim_nd), INTENT(IN) :: ind_in
      INTEGER, DIMENSION(2)                              :: ind_out
      INTEGER :: i
      INTEGER, DIMENSION(${maxrank}$)                    :: ind_tmp

      DO i = 1, map%ndim1_2d
         ind_tmp(i) = ind_in(map%map1_2d(i))
      END DO
      ind_out(1) = combine_pgrid_index(ind_tmp(:map%ndim1_2d), map%dims1_2d)

      DO i = 1, map%ndim2_2d
         ind_tmp(i) = ind_in(map%map2_2d(i))
      END DO
      ind_out(2) = combine_pgrid_index(ind_tmp(:map%ndim2_2d), map%dims2_2d)
   END FUNCTION

! **************************************************************************************************
!> \brief transform 2d index to nd index, using info from index mapping.
!> \param map index mapping
!> \param ind_in 2d index
!> \param ind_out nd index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION get_nd_indices_tensor(map, ind_in) RESULT(ind_out)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(IN)      :: ind_in
      INTEGER, DIMENSION(map%ndim_nd)                    :: ind_out
      INTEGER, DIMENSION(${maxrank}$)                    :: ind_tmp
      INTEGER                                            :: i

      ind_tmp(:map%ndim1_2d) = split_tensor_index(ind_in(1), map%dims1_2d)

      DO i = 1, map%ndim1_2d
         ind_out(map%map1_2d(i)) = ind_tmp(i)
      END DO

      ind_tmp(:map%ndim2_2d) = split_tensor_index(ind_in(2), map%dims2_2d)

      DO i = 1, map%ndim2_2d
         ind_out(map%map2_2d(i)) = ind_tmp(i)
      END DO

   END FUNCTION

! **************************************************************************************************
!> \brief transform 2d index to nd index, using info from index mapping.
!> \param map index mapping
!> \param ind_in 2d index
!> \param ind_out nd index
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION get_nd_indices_pgrid(map, ind_in) RESULT(ind_out)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map
      INTEGER, DIMENSION(2), INTENT(IN)                  :: ind_in
      INTEGER, DIMENSION(map%ndim_nd)                    :: ind_out

      ind_out(map%map1_2d) = split_pgrid_index(ind_in(1), map%dims1_2d)
      ind_out(map%map2_2d) = split_pgrid_index(ind_in(2), map%dims2_2d)

   END FUNCTION

! **************************************************************************************************
!> \brief Invert order
!> \author Patrick Seewald
! **************************************************************************************************
   PURE FUNCTION dbt_inverse_order(order)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: order
      INTEGER, DIMENSION(SIZE(order))                    :: dbt_inverse_order

      INTEGER                                            :: i

      dbt_inverse_order(order) = (/(i, i=1, SIZE(order))/)
   END FUNCTION

! **************************************************************************************************
!> \brief reorder tensor index (no data)
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE permute_index(map_in, map_out, order)
      TYPE(nd_to_2d_mapping), INTENT(IN)                 :: map_in
      TYPE(nd_to_2d_mapping), INTENT(OUT)                :: map_out
      INTEGER, DIMENSION(ndims_mapping(map_in)), &
         INTENT(IN)                                      :: order

      INTEGER                                            :: ndim_nd
      INTEGER, DIMENSION(ndims_mapping_row(map_in))       :: map1_2d, map1_2d_reorder
      INTEGER, DIMENSION(ndims_mapping_column(map_in))    :: map2_2d, map2_2d_reorder
      INTEGER, DIMENSION(ndims_mapping(map_in))          :: dims_nd, dims_reorder

      CALL dbt_get_mapping_info(map_in, ndim_nd, dims_nd=dims_nd, map1_2d=map1_2d, map2_2d=map2_2d)

      dims_reorder(order) = dims_nd

      map1_2d_reorder(:) = order(map1_2d)
      map2_2d_reorder(:) = order(map2_2d)

      CALL create_nd_to_2d_mapping(map_out, dims_reorder, map1_2d_reorder, map2_2d_reorder)
   END SUBROUTINE

END MODULE dbt_index
