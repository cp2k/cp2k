!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for evaluating the residual part (1/r^3) of Integrals for
!>        semi-empiric methods
!> \author Teodoro Laino (11.2008) [tlaino]
! *****************************************************************************
MODULE semi_empirical_int3_utils

  USE f77_blas
  USE input_constants,                 ONLY: do_method_pchg
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_arrays,       ONLY: clm_d,&
                                             indexb
  USE semi_empirical_types,            ONLY: semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module   = .FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int3_utils'

  PUBLIC ::   ijkl_low_3, charg_int_3, dcharg_int_3, coeff_int_3

CONTAINS

! *****************************************************************************
!> \brief Low level general driver for computing residual part of semi-empirical
!>        integrals <ij|kl> and their derivatives
!>        The residual part is the leading 1/r^3 term
!>
!> \author Teodoro Laino [tlaino]
!> \date 11.2008 [tlaino]
! *****************************************************************************
  FUNCTION ijkl_low_3(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, itype, eval, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    INTEGER, INTENT(IN)                      :: itype
    REAL(KIND=dp)                            :: eval
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl_low_3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2, lij, lkl
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: add, ccc, chrg, pij, pkl, sum

    failure = .FALSE.
    sum = 0.0_dp
    l1  = ABS (li-lj)
    lij = indexb(li+1, lj+1)
    l2  = ABS (lk-ll)
    lkl = indexb(lk+1, ll+1)

    ! Standard value of the integral
    IF (l1 == 0) THEN
       IF (lij == 1) THEN
          pij = sepi%ko(1)
          IF (ic == 1) THEN
             pij = sepi%ko(9)
          END IF
       ELSE IF (lij == 3) THEN
          pij = sepi%ko(7)
       ELSE IF (lij == 6) THEN
          pij = sepi%ko(8)
       END IF
    END IF
    !
    IF (l2 == 0) THEN
       IF (lkl == 1) THEN
          pkl = sepj%ko(1)
          IF (ic == 2) THEN
             pkl = sepj%ko(9)
          END IF
       ELSE IF (lkl == 3) THEN
          pkl = sepj%ko(7)
       ELSE IF (lkl == 6) THEN
          pkl = sepj%ko(8)
       END IF
    END IF
    IF (l1==0.AND.l2==0) THEN
       IF (itype==do_method_pchg) THEN
          add = 0.0_dp
       ELSE
          add = (pij+pkl) ** 2
       END IF
       ccc = clm_d(ij, l1, 0) * clm_d(kl, l2, 0)
       IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
          chrg = eval(r, l1, l2, add, error)
          sum  = chrg
       END IF
    END IF
    res = sum
  END FUNCTION ijkl_low_3

! *****************************************************************************
!> \brief Evaluates the residual Interaction function between two point-charges
!>        The term evaluated is the 1/r^3 (for short range interactions)
!>        r    -  Distance r12
!>        l1   -  Quantum numbers for multipole of configuration 1
!>        l2   -  Quantum numbers for multipole of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino]
!> \date 11.2008 [tlaino]
! *****************************************************************************
  FUNCTION charg_int_3(r, l1, l2, add, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1, l2
    REAL(KIND=dp), INTENT(in)                :: add
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'charg_int_3', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Computing only residual Integral Values
    charg = 0.0_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg =  -add/(2.0_dp*r**3)
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION charg_int_3

! *****************************************************************************
!> \brief Evaluates the coefficient for the residual Interaction function
!>        between two point-charges
!>        l1   -  Quantum numbers for multipole of configuration 1
!>        l2   -  Quantum numbers for multipole of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino]
!> \date 11.2008 [tlaino]
! *****************************************************************************
  FUNCTION coeff_int_3(r, l1, l2, add, error) RESULT(coeff)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1, l2
    REAL(KIND=dp), INTENT(in)                :: add
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: coeff

    CHARACTER(len=*), PARAMETER :: routineN = 'coeff_int_3', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Computing only residual Integral Values
    coeff = 0.0_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       coeff =  -add/2.0_dp
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION coeff_int_3

! *****************************************************************************
!> \brief Derivatives of residual interaction function between two point-charges
!>
!>        r    -  Distance r12
!>        l1   -  Quantum numbers for multipole of configuration 1
!>        l2   -  Quantum numbers for multipole of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino]
!> \date 11.2008 [tlaino]
! *****************************************************************************
  FUNCTION dcharg_int_3(r, l1, l2, add, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1, l2
    REAL(KIND=dp), INTENT(in)                :: add
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'dcharg_int_3', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Computing only residual Integral Derivatives
    charg = 0.0_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg =  3.0_dp*add/(2.0_dp*r**4)
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION dcharg_int_3

END MODULE semi_empirical_int3_utils
