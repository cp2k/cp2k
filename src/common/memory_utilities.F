!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for the memory handling.
!> \par History
!>      (12.2017) remove stop_memory
!> \author Matthias Krack (25.06.1999)
! **************************************************************************************************
MODULE memory_utilities

   USE kinds, ONLY: dp, int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'memory_utilities'

   PUBLIC :: reallocate

   INTERFACE reallocate
      MODULE PROCEDURE reallocate_c1, reallocate_c2, reallocate_c3, reallocate_c4, &
         reallocate_i1, reallocate_i2, reallocate_i3, reallocate_i4, &
         reallocate_r1, reallocate_r2, reallocate_r3, reallocate_r4, &
         reallocate_r5, reallocate_s1, reallocate_l1, reallocate_8i1, &
         reallocate_8i2
   END INTERFACE

CONTAINS

#! *************************************************************************************************
#!> \brief Fypp macro for common subroutine body
#!> \author Ole Schuett
#!> \author Tiziano Müller
#! *************************************************************************************************
#:def reallocate(suffix, rank, type, zero, worktype=None)
#:set bounds_vars = ','.join("lb{0}_new,ub{0}_new".format(i+1) for i in range(rank))
#:set old_bounds = ','.join(['lb{0}:ub{0}'.format(i+1) for i in range(rank)])
#:set new_bounds = ','.join(['lb{0}_new:ub{0}_new'.format(i+1) for i in range(rank)])
#:set arr_exp = ','.join(':'*rank)
! **************************************************************************************************
!> \brief (Re)Allocate a ${rank}$D vector of type ${type}$ with new dimensions (but same shape)
!> \param p pointer to the existing data, if NULL() calling this is equivalent to an ALLOCATE(...)
#:for i in range(1, rank+1)
!> \param lb${i}$_new new lower bound for dimension ${i}$
!> \param ub${i}$_new new upper bound for dimension ${i}$
#:endfor
! **************************************************************************************************
   SUBROUTINE reallocate_${suffix}$${rank}$ (p, ${bounds_vars}$)
      ${type}$, &
         DIMENSION(${arr_exp}$), &
         POINTER, INTENT(INOUT) :: p

      INTEGER, INTENT(IN) :: &
         ${bounds_vars}$

#:for i in range(1, rank+1)
      INTEGER :: lb${i}$, lb${i}$_old, ub${i}$, ub${i}$_old
#:endfor

#:if worktype
      ${worktype}$, &
         #:else
         ${type}$, &
         #:endif
         DIMENSION(${arr_exp}$), &
         POINTER :: work

      NULLIFY (work)

      IF (ASSOCIATED(p)) THEN
#:for i in range(1, rank+1)
         lb${i}$_old = LBOUND(p, ${i}$)
         ub${i}$_old = UBOUND(p, ${i}$)
         lb${i}$ = MAX(lb${i}$_new, lb${i}$_old)
         ub${i}$ = MIN(ub${i}$_new, ub${i}$_old)
#:endfor
         work => p
      END IF

      ALLOCATE (p(${new_bounds}$))
      p = ${zero}$

      IF (ASSOCIATED(work)) THEN
         p(${old_bounds}$) = work(${old_bounds}$)
         DEALLOCATE (work)
      END IF

   END SUBROUTINE reallocate_${suffix}$${rank}$
#:enddef

$: reallocate(suffix="c",  rank=1, type="COMPLEX(KIND=dp)", zero="(0.0_dp, 0.0_dp)")
$: reallocate(suffix="c",  rank=2, type="COMPLEX(KIND=dp)", zero="(0.0_dp, 0.0_dp)")
$: reallocate(suffix="c",  rank=3, type="COMPLEX(KIND=dp)", zero="(0.0_dp, 0.0_dp)")
$: reallocate(suffix="c",  rank=4, type="COMPLEX(KIND=dp)", zero="(0.0_dp, 0.0_dp)")
$: reallocate(suffix="i",  rank=1, type="INTEGER", zero="0")
$: reallocate(suffix="i",  rank=2, type="INTEGER", zero="0")
$: reallocate(suffix="i",  rank=3, type="INTEGER", zero="0")
$: reallocate(suffix="i",  rank=4, type="INTEGER", zero="0")
$: reallocate(suffix="8i", rank=1, type="INTEGER(KIND=int_8)", zero="0")
$: reallocate(suffix="8i", rank=2, type="INTEGER(KIND=int_8)", zero="0")
$: reallocate(suffix="r",  rank=1, type="REAL(KIND=dp)", zero="0.0_dp")
$: reallocate(suffix="r",  rank=2, type="REAL(KIND=dp)", zero="0.0_dp")
$: reallocate(suffix="r",  rank=3, type="REAL(KIND=dp)", zero="0.0_dp")
$: reallocate(suffix="r",  rank=4, type="REAL(KIND=dp)", zero="0.0_dp")
$: reallocate(suffix="r",  rank=5, type="REAL(KIND=dp)", zero="0.0_dp")
$: reallocate(suffix="l",  rank=1, type="LOGICAL", zero=".FALSE.")
$: reallocate(suffix="s",  rank=1, type="CHARACTER(LEN=*)", zero='""', worktype="CHARACTER(LEN=LEN(p))")

END MODULE memory_utilities
