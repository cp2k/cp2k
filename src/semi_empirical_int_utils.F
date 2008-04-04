!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Integrals for semi-empiric methods
!> \author Teodoro Laino (03.2008) [tlaino] - University of Zurich 
! *****************************************************************************
MODULE semi_empirical_int_utils

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_arrays,       ONLY: clm,&
                                             indexb
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             taper_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_utils'

  PUBLIC :: taper_eval, dtaper_eval
  
CONTAINS

! *****************************************************************************
!> \brief Taper functions
! *****************************************************************************
  FUNCTION taper_eval (taper,rij, error) RESULT(ft)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rij
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: ft

    CHARACTER(len=*), PARAMETER :: routineN = 'taper_eval', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr, ep

    ft = 1._dp
    dr = -taper%al*(rij-taper%r0)
    IF ( dr < 20._dp ) THEN
       ep = EXP(dr)
       ft = ep/(1._dp+ep)
    END IF
  END FUNCTION taper_eval
  
! *****************************************************************************
!> \brief Analytical derivatives for taper function
! *****************************************************************************
  FUNCTION dtaper_eval (taper, rij, error ) RESULT(dft)
    TYPE(taper_type), POINTER                :: taper
    REAL(KIND=dp), INTENT(IN)                :: rij
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: dft

    CHARACTER(len=*), PARAMETER :: routineN = 'dtaper_eval', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dr, ep

    dft = 0.0_dp
    dr = -taper%al*(rij-taper%r0)
    IF ( dr < 20._dp ) THEN
       ep = EXP(dr)
       dft = 1._dp/(1._dp+ep)**2 * ep * (-taper%al)
    END IF

  END FUNCTION dtaper_eval

! *****************************************************************************
!> \brief General driver for computing semi-empirical integrals <ij|kl>
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 03.2008 [tlaino]
! *****************************************************************************
  FUNCTION ijkl (sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, error) RESULT(res)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l1max, l1min, l2, l2max, &
                                                l2min, lij, lkl, lmin, m, mm
    REAL(KIND=dp)                            :: add, ccc, dij, dkl, pij, pkl, &
                                                s1, sum

    l1min = ABS (li-lj)
    l1max = li + lj
    lij   = indexb(li+1, lj+1)
    l2min = ABS (lk-ll)
    l2max = lk + ll
    lkl   = indexb(lk+1, ll+1)

    l1max = MIN (l1max, 2)
    l1min = MIN (l1min, 2)
    l2max = MIN (l2max, 2)
    l2min = MIN (l2min, 2)
    sum = 0.0_dp
    
    DO l1 = l1min, l1max
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
       ELSE
          dij = sepi%cs(lij)
          pij = sepi%ko(lij)
       END IF
       !
       DO l2 = l2min, l2max
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
          ELSE
             dkl = sepj%cs(lkl)
             pkl = sepj%ko(lkl)
          END IF
          !
          add = (pij+pkl) ** 2
          lmin = MIN(l1, l2)
          s1 = 0.d00
          DO m = -lmin, lmin
             ccc = clm(ij, l1, m) * clm(kl, l2, m)
             IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
                mm = ABS (m)
                s1 = s1 + charg_int(r, l1, l2, mm, dij, dkl, add, error) * ccc
             END IF
          END DO
          sum = sum + s1
       END DO
    END DO
    res = sum
  END FUNCTION ijkl

! *****************************************************************************
!> \brief Interaction function between two point-charge configurations (MNDO/d)
!>        r    -  Distance r12
!>        l1,m -  Quantum numbers for multipole of configuration 1
!>        l2,m -  Quantum numbers for multipole of configuration 2
!>        da   -  charge separation of configuration 1
!>        db   -  charge separation of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 03.2008 [tlaino]
! *****************************************************************************
  FUNCTION charg_int(r, l1, l2, m, da, db, add, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1, l2, m
    REAL(KIND=dp), INTENT(in)                :: da, db, add
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'charg_int', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: aa, ab, dxdx, dxqxz, dzdz, &
                                                dzqzz, qqzz, qxzdx, qxzqxz, &
                                                qzzdz, qzzq, xyxy, zzzz

    charg = 0.0_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg = 1.0_dp/SQRT (r**2+add)
       RETURN
    END IF
    ! Z - Q.
    IF (l1 == 1 .AND. l2 == 0) THEN
       charg = -1.0_dp/SQRT((r+da)**2+add) + 1.0_dp/SQRT((r-da)**2+add)
       charg = charg/2.0_dp
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1) THEN
       charg =  1.0_dp/SQRT((r+db)**2+add) - 1.0_dp/SQRT((r-db)**2+add)
       charg = charg/2.0_dp
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 0) THEN
       dzdz  = &
            +1.0_dp/SQRT((r+da-db)**2+add) + 1.0_dp/SQRT((r-da+db)**2+add) &
            -1.0_dp/SQRT((r-da-db)**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzdz/4.0_dp
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 1) THEN
       dxdx  = 2.0_dp/SQRT(r**2+(da-db)**2+add) - 2.0_dp/SQRT(r**2+(da+db)**2+add)
       charg = dxdx*0.25_dp
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2) THEN
       qqzz  = 1.0_dp/SQRT((r-db)**2+add) - 2.0_dp/SQRT(r**2+db**2+add) + 1.0_dp/SQRT((r+db)**2+add)
       charg = qqzz/4.0_dp
       RETURN
    END IF
    ! ZZ -Q
    IF (l1 == 2 .AND. l2 == 0) THEN
       qzzq  = 1.0_dp/SQRT((r-da)**2+add) - 2.0_dp/SQRT(r**2+da**2+add) + 1.0_dp/SQRT((r+da)**2+add)
       charg = qzzq/4.0_dp
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 0) THEN
       dzqzz = &
            +1.0_dp/SQRT((r-da-db)**2+add)    - 2.0_dp/SQRT((r-da)**2+db**2+add) &
            +1.0_dp/SQRT((r+db-da)**2+add)    - 1.0_dp/SQRT((r-db+da)**2+add)    &
            +2.0_dp/SQRT((r+da)**2+db**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzqzz / 8.0_dp
       RETURN
    END IF
    ! ZZ - Z
    IF (l1 == 2 .AND. l2 == 1 .AND. m == 0) THEN
       qzzdz =&
            -1.0_dp/SQRT((r-da-db)**2+add)    + 2.0_dp/SQRT((r-db)**2+da**2+add) &
            -1.0_dp/SQRT((r+da-db)**2+add)    + 1.0_dp/SQRT((r-da+db)**2+add)    &
            -2.0_dp/SQRT((r+db)**2+da**2+add) + 1.0_dp/SQRT((r+da+db)**2+add)
       charg = qzzdz / 8.0_dp
       RETURN
    END IF
    ! ZZ - ZZ
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 0) THEN
       zzzz = &
            +1.0_dp/SQRT((r-da-db)**2+add)    + 1.0_dp/SQRT((r+da+db)**2+add)    &
            +1.0_dp/SQRT((r-da+db)**2+add)    + 1.0_dp/SQRT((r+da-db)**2+add)    &
            -2.0_dp/SQRT((r-da)**2+db**2+add) - 2.0_dp/SQRT((r-db)**2+da**2+add) &
            -2.0_dp/SQRT((r+da)**2+db**2+add) - 2.0_dp/SQRT((r+db)**2+da**2+add) &
            +2.0_dp/SQRT(r**2+(da-db)**2+add) + 2.0_dp/SQRT(r**2+(da+db)**2+add)
       xyxy = &
            +4.0_dp/SQRT(r**2+(da-db)**2+add) + 4.0_dp/SQRT(r**2+(da+db)**2+add) &
            -8.0_dp/SQRT(r**2+da**2+db**2+add)
       charg = (zzzz/16.0_dp - xyxy/64.0_dp)
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 1) THEN
       ab = db/SQRT(2.0_dp)
       dxqxz = &
            -2.0_dp/SQRT((r-ab)**2+(da-ab)**2+add) + 2.0_dp/SQRT((r+ab)**2+(da-ab)**2+add) &
            +2.0_dp/SQRT((r-ab)**2+(da+ab)**2+add) - 2.0_dp/SQRT((r+ab)**2+(da+ab)**2+add)
       charg = dxqxz/8.0_dp
       RETURN
    END IF
    ! ZX - X
    IF (l1 == 2 .AND. l2 == 1 .AND. m == 1) THEN
       aa = da/SQRT(2.0_dp)
       qxzdx = &
            -2.0_dp/SQRT((r+aa)**2+(aa-db)**2+add) + 2.0_dp/SQRT((r-aa)**2+(aa-db)**2+add) &
            +2.0_dp/SQRT((r+aa)**2+(aa+db)**2+add) - 2.0_dp/SQRT((r-aa)**2+(aa+db)**2+add)
       charg = qxzdx/8.0_dp
       RETURN
    END IF
    ! ZX - ZX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 1) THEN
       aa = da/SQRT (2.0_dp)
       ab = db/SQRT (2.0_dp)
       qxzqxz = &
            +2.0_dp/SQRT((r+aa-ab)**2+(aa-ab)**2+add) - 2.0_dp/SQRT((r+aa+ab)**2+(aa-ab)**2+add) &
            -2.0_dp/SQRT((r-aa-ab)**2+(aa-ab)**2+add) + 2.0_dp/SQRT((r-aa+ab)**2+(aa-ab)**2+add) &
            -2.0_dp/SQRT((r+aa-ab)**2+(aa+ab)**2+add) + 2.0_dp/SQRT((r+aa+ab)**2+(aa+ab)**2+add) &
            +2.0_dp/SQRT((r-aa-ab)**2+(aa+ab)**2+add) - 2.0_dp/SQRT((r-aa+ab)**2+(aa+ab)**2+add)
       charg = qxzqxz/16.0_dp
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 2) THEN
       xyxy = 4.0_dp/SQRT(r**2+(da-db)**2+add) + 4.0_dp/SQRT(r**2+(da+db)**2+add) - 8.0_dp/SQRT(r**2+da**2+db**2+add)
       charg = xyxy/16.0_dp
       RETURN
    END IF
  END FUNCTION charg_int

END MODULE semi_empirical_int_utils
