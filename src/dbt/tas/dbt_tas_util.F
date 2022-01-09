!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief often used utilities for tall-and-skinny matrices
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbt_tas_util
   USE kinds, ONLY: int_4, int_8
   USE util, ONLY: sort

#include "../../base/base_uses.f90"

! #if TO_VERSION(1, 11) <= TO_VERSION(LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
!    USE libxsmm, ONLY: libxsmm_diff
! #  define PURE_ARRAY_EQ
! #else
#  define PURE_ARRAY_EQ PURE
! #endif

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_util'

   PUBLIC :: &
      array_eq, &
      index_unique, &
      swap, &
      generate_larnv_seed

   INTERFACE swap
      MODULE PROCEDURE swap_i8
      MODULE PROCEDURE swap_i
   END INTERFACE

   INTERFACE array_eq
      MODULE PROCEDURE array_eq_i8
      MODULE PROCEDURE array_eq_i
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param arr ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE swap_i8(arr)
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(INOUT)   :: arr

      INTEGER(KIND=int_8)                                :: tmp

      tmp = arr(1)
      arr(1) = arr(2)
      arr(2) = tmp
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param arr ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE swap_i(arr)
      INTEGER, DIMENSION(2), INTENT(INOUT)               :: arr

      INTEGER                                            :: tmp

      tmp = arr(1)
      arr(1) = arr(2)
      arr(2) = tmp
   END SUBROUTINE

! **************************************************************************************************
!> \brief Get all unique elements in index_in
!> \param index_in ...
!> \param index_out ...
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE index_unique(index_in, index_out)
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_in
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: index_out

      INTEGER                                            :: blk, count, orig_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_tmp
      INTEGER, DIMENSION(2)                              :: prev_index
      INTEGER, DIMENSION(1:SIZE(index_in, 1), 1:SIZE(&
         index_in, 2))                                   :: index_sorted

      orig_size = SIZE(index_in, 1)
      ALLOCATE (index_tmp(orig_size, 2))
      index_sorted(:, :) = index_in(:, :)
      CALL sort_indices(orig_size, index_sorted(:, 1), index_sorted(:, 2))
      count = 0
      prev_index(:) = [0, 0]
      DO blk = 1, orig_size
         IF (ANY(index_sorted(blk, :) .NE. prev_index(:))) THEN
            count = count + 1
            index_tmp(count, :) = index_sorted(blk, :)
            prev_index(:) = index_sorted(blk, :)
         END IF
      END DO

      ALLOCATE (index_out(count, 2))
      index_out(:, :) = index_tmp(1:count, :)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Sorts the rows & columns of a work matrix
!>        Sorts the row and column indices so that the rows monotonically
!>        increase and the columns monotonically increase within each row.
!>        Passing the blk_p array rearranges the block pointers accordingly.
!>        This must be done if they are pointing to valid data, otherwise
!>        they become invalid.
!>        TODO: Simplify using array_sort.fypp and merge with index_unique above.
!> \param n number of blocks (elements) to sort
!> \param row_i row indices
!> \param col_i column indices
!> \param blk_p block pointers
!> \param blk_d data storage
!> \author Patrick Seewald
! **************************************************************************************************
   SUBROUTINE sort_indices(n, row_i, col_i, blk_p, blk_d)
      INTEGER, INTENT(IN)                                :: n
      INTEGER, DIMENSION(1:), INTENT(INOUT)              :: row_i, col_i
      INTEGER, DIMENSION(1:), INTENT(INOUT), OPTIONAL    :: blk_p, blk_d

      CHARACTER(len=*), PARAMETER                        :: routineN = 'sort_indices'
      INTEGER(KIND=int_8), PARAMETER                     :: lmask8 = 4294967295_int_8

      INTEGER                                            :: error_handle, i
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: sort_keys
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: buf, buf_d

!   ---------------------------------------------------------------------------

      IF (SIZE(row_i) .EQ. 0) RETURN

      CALL timeset(routineN, error_handle)

      IF (SIZE(row_i) < n) CPABORT('row_i too small')
      IF (SIZE(col_i) < n) CPABORT('col_i too small')
      IF (PRESENT(blk_p)) THEN
         IF (SIZE(blk_p) < n) CPABORT('blk_p too small')
         ALLOCATE (buf(n))
         buf(1:n) = blk_p(1:n)
      END IF
      IF (PRESENT(blk_d)) THEN
         ALLOCATE (buf_d(n))
         buf_d(1:n) = blk_d(1:n)
      END IF
      ! Create an ordering for both rows and columns. If the blk_p must
      ! be rearranged, then the col_i array will be used as a
      ! permutation vector.
      ALLOCATE (sort_keys(n))
      sort_keys(:) = IOR(ISHFT(INT(row_i(1:n), int_8), 32), INT(col_i(1:n), int_8))
      IF (PRESENT(blk_p)) col_i(1:n) = (/(i, i=1, n)/)
      ! Now do a nice quicksort.
      CALL sort(sort_keys, n, col_i)
      ! Since blk_d is usually not present we can have two loops that
      ! are essentially the same.
      IF (PRESENT(blk_p)) THEN
         DO i = 1, n
            blk_p(i) = buf(col_i(i))
         END DO
         DEALLOCATE (buf)
      END IF
      IF (PRESENT(blk_d)) THEN
         DO i = 1, n
            blk_d(i) = buf_d(col_i(i))
         END DO
         DEALLOCATE (buf_d)
      END IF
      DO i = 1, n
         col_i(i) = INT(IAND(sort_keys(i), lmask8), int_4)
         row_i(i) = INT(ISHFT(sort_keys(i), -32), int_4)
      END DO
      DEALLOCATE (sort_keys)
      CALL timestop(error_handle)

   END SUBROUTINE sort_indices

! **************************************************************************************************
!> \brief ...
!> \param arr1 ...
!> \param arr2 ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   PURE_ARRAY_EQ FUNCTION array_eq_i(arr1, arr2)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: arr1, arr2
      LOGICAL                                            :: array_eq_i

! #if TO_VERSION(1, 11) <= TO_VERSION(LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
!       array_eq_i = .NOT. libxsmm_diff(arr1, arr2)
! #else
      array_eq_i = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)
! #endif
   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param arr1 ...
!> \param arr2 ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   PURE_ARRAY_EQ FUNCTION array_eq_i8(arr1, arr2)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: arr1, arr2
      LOGICAL                                            :: array_eq_i8

! #if TO_VERSION(1, 11) <= TO_VERSION(LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
!       array_eq_i8 = .NOT. libxsmm_diff(arr1, arr2)
! #else
      array_eq_i8 = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i8 = ALL(arr1 == arr2)
! #endif
   END FUNCTION

! **************************************************************************************************
!> \brief Generate a seed respecting the lapack constraints,
!>        - values between 0..4095 (2**12-1)
!>        - iseed(4) odd
!>        also try to avoid iseed that are zero.
!>        Have but with a unique 2D mapping (irow,icol), and a 'mini-seed' ival
!>
!>           TODO: Move to some utility module.
!>
!> \param irow 1..nrow
!> \param nrow ...
!> \param icol 1..ncol
!> \param ncol ...
!> \param ival mini-seed
!> \return a lapack compatible seed
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION generate_larnv_seed(irow, nrow, icol, ncol, ival) RESULT(iseed)

      INTEGER, INTENT(IN)                                :: irow, nrow, icol, ncol, ival
      INTEGER                                            :: iseed(4)

      INTEGER(KIND=int_8)                                :: map

      map = ((irow - 1 + icol*INT(nrow, int_8))*(1 + MODULO(ival, 2**16)))*2 + 1 + 0*ncol ! ncol used
      iseed(4) = INT(MODULO(map, 2_int_8**12)); map = map/2_int_8**12; ! keep odd
      iseed(3) = INT(MODULO(IEOR(map, 3541_int_8), 2_int_8**12)); map = map/2_int_8**12
      iseed(2) = INT(MODULO(IEOR(map, 1153_int_8), 2_int_8**12)); map = map/2_int_8**12
      iseed(1) = INT(MODULO(IEOR(map, 2029_int_8), 2_int_8**12)); map = map/2_int_8**12
   END FUNCTION generate_larnv_seed

END MODULE
