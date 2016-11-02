!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods for tensor index and mapping to DBCSR index.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsrt_index
   USE dbcsrt_index_types,              ONLY: nd_to_2d_mapping,&
                                              ndims_mapping
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsrt_index'

   PUBLIC :: &
      get_2d_indices, &
      get_nd_indices, &
      combine_index, &
      split_index

CONTAINS

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

END MODULE dbcsrt_index
