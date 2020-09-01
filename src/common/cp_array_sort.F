!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routine for sorting an array
!> \note
!>      CP2K:
!>      Please use the interface defined in util.F for calling sort().
!>
!>      DBCSR:
!>      Please use the interface defined in dbcsr_toollib.F for calling sort().
!> \par History
!>      12.2012 first version [ole]
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_array_sort

   USE kinds, ONLY: sp, dp, int_4, int_8

#include "../base/base_uses.f90"

#:include 'array_sort.fypp'

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cp_1d_s_sort, cp_1d_r_sort, cp_1d_i4_sort, cp_1d_i8_sort

CONTAINS

#:call array_sort(prefix='cp_1d_s', type='REAL(kind=sp)')
#:endcall

   PURE FUNCTION cp_1d_s_less_than(a, b) RESULT(res)
      REAL(kind=sp), INTENT(IN) :: a, b
      LOGICAL                   :: res
      res = a < b
   END FUNCTION cp_1d_s_less_than

#:call array_sort(prefix='cp_1d_r', type='REAL(kind=dp)')
#:endcall

   PURE FUNCTION cp_1d_r_less_than(a, b) RESULT(res)
      REAL(kind=dp), INTENT(IN) :: a, b
      LOGICAL                   :: res
      res = a < b
   END FUNCTION cp_1d_r_less_than

#:call array_sort(prefix='cp_1d_i4', type='INTEGER(kind=int_4)')
#:endcall

   PURE FUNCTION cp_1d_i4_less_than(a, b) RESULT(res)
      INTEGER(kind=int_4), INTENT(IN) :: a, b
      LOGICAL                         :: res
      res = a < b
   END FUNCTION cp_1d_i4_less_than

#:call array_sort(prefix='cp_1d_i8', type='INTEGER(kind=int_8)')
#:endcall

   PURE FUNCTION cp_1d_i8_less_than(a, b) RESULT(res)
      INTEGER(kind=int_8), INTENT(IN) :: a, b
      LOGICAL                         :: res
      res = a < b
   END FUNCTION cp_1d_i8_less_than

END MODULE cp_array_sort
