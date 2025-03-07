!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
#if defined(__LIBXSMM)
#include "libxsmm_version.h"
#endif

#if CPVERSION_CHECK(1, 11, <=, LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
   USE libxsmm, ONLY: libxsmm_diff
#  define PURE_ARRAY_EQ
#else
#  define PURE_ARRAY_EQ PURE
#endif

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbt_tas_util'

   PUBLIC :: array_eq, swap

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

#if CPVERSION_CHECK(1, 11, <=, LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
      array_eq_i = .NOT. libxsmm_diff(arr1, arr2)
#else
      array_eq_i = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)
#endif
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

#if CPVERSION_CHECK(1, 11, <=, LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
      array_eq_i8 = .NOT. libxsmm_diff(arr1, arr2)
#else
      array_eq_i8 = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i8 = ALL(arr1 == arr2)
#endif
   END FUNCTION

END MODULE
