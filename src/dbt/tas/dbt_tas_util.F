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
