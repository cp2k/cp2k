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
                                             ijkl_ind,&
                                             indexa,&
                                             indexb,&
                                             int2c_type
  USE semi_empirical_types,            ONLY: rotmat_type,&
                                             semi_empirical_type,&
                                             taper_type
#include "cp_common_uses.h"

  IMPLICIT NONE

#include "semi_empirical_int_debug.h"

  PRIVATE
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_utils'

  PUBLIC :: taper_eval, dtaper_eval, ijkl, rotmat, rot_2el_2c_first, store_2el_2c_diag
  
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
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l1max, l1min, l2, l2max, &
                                                l2min, lij, lkl, lmin, m, mm
    REAL(KIND=dp)                            :: add, ccc, chrg, dij, dkl, &
                                                pij, pkl, s1, sum

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
          s1 = 0.0_dp
          DO m = -lmin, lmin
             ccc = clm(ij, l1, m) * clm(kl, l2, m)
             IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
                mm   = ABS (m)
                chrg = charg_int(r, l1, l2, mm, dij, dkl, add, error)
                s1   = s1 + chrg * ccc
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
       charg = 1.0_dp/SQRT(r**2+add)
       RETURN
    END IF
    ! Z - Q.
    IF (l1 == 1 .AND. l2 == 0) THEN
       charg = -1.0_dp/SQRT((r+da)**2+add) + 1.0_dp/SQRT((r-da)**2+add)
       charg = charg*0.5_dp
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1) THEN
       charg =  1.0_dp/SQRT((r+db)**2+add) - 1.0_dp/SQRT((r-db)**2+add)
       charg = charg*0.5_dp
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 0) THEN
       dzdz  = &
            +1.0_dp/SQRT((r+da-db)**2+add) + 1.0_dp/SQRT((r-da+db)**2+add) &
            -1.0_dp/SQRT((r-da-db)**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzdz*0.25_dp
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
       charg = qqzz*0.25_dp
       RETURN
    END IF
    ! ZZ -Q
    IF (l1 == 2 .AND. l2 == 0) THEN
       qzzq  = 1.0_dp/SQRT((r-da)**2+add) - 2.0_dp/SQRT(r**2+da**2+add) + 1.0_dp/SQRT((r+da)**2+add)
       charg = qzzq*0.25_dp
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 0) THEN
       dzqzz = &
            +1.0_dp/SQRT((r-da-db)**2+add)    - 2.0_dp/SQRT((r-da)**2+db**2+add) &
            +1.0_dp/SQRT((r+db-da)**2+add)    - 1.0_dp/SQRT((r-db+da)**2+add)    &
            +2.0_dp/SQRT((r+da)**2+db**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzqzz*0.125_dp
       RETURN
    END IF
    ! ZZ - Z
    IF (l1 == 2 .AND. l2 == 1 .AND. m == 0) THEN
       qzzdz =&
            -1.0_dp/SQRT((r-da-db)**2+add)    + 2.0_dp/SQRT((r-db)**2+da**2+add) &
            -1.0_dp/SQRT((r+da-db)**2+add)    + 1.0_dp/SQRT((r-da+db)**2+add)    &
            -2.0_dp/SQRT((r+db)**2+da**2+add) + 1.0_dp/SQRT((r+da+db)**2+add)
       charg = qzzdz*0.125_dp
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
       charg = (zzzz*0.0625_dp - xyxy*0.015625_dp)
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 1) THEN
       ab = db/SQRT(2.0_dp)
       dxqxz = &
            -2.0_dp/SQRT((r-ab)**2+(da-ab)**2+add) + 2.0_dp/SQRT((r+ab)**2+(da-ab)**2+add) &
            +2.0_dp/SQRT((r-ab)**2+(da+ab)**2+add) - 2.0_dp/SQRT((r+ab)**2+(da+ab)**2+add)
       charg = dxqxz*0.125_dp
       RETURN
    END IF
    ! ZX - X
    IF (l1 == 2 .AND. l2 == 1 .AND. m == 1) THEN
       aa = da/SQRT(2.0_dp)
       qxzdx = &
            -2.0_dp/SQRT((r+aa)**2+(aa-db)**2+add) + 2.0_dp/SQRT((r-aa)**2+(aa-db)**2+add) &
            +2.0_dp/SQRT((r+aa)**2+(aa+db)**2+add) - 2.0_dp/SQRT((r-aa)**2+(aa+db)**2+add)
       charg = qxzdx*0.125_dp
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
       charg = qxzqxz*0.0625_dp
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 2) THEN
       xyxy = 4.0_dp/SQRT(r**2+(da-db)**2+add) + 4.0_dp/SQRT(r**2+(da+db)**2+add) - 8.0_dp/SQRT(r**2+da**2+db**2+add)
       charg = xyxy*0.0625_dp
       RETURN
    END IF
  END FUNCTION charg_int

! *****************************************************************************
!> \brief Computes the general rotation matrices for the integrals
!>        between I and J (J>=I)
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  RECURSIVE SUBROUTINE rotmat (sepi, sepj, rjiv, r, ij_matrix, do_derivatives, do_invert, debug_invert, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rjiv
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(rotmat_type), POINTER               :: ij_matrix
    LOGICAL, INTENT(IN)                      :: do_derivatives
    LOGICAL, INTENT(OUT), OPTIONAL           :: do_invert
    LOGICAL, INTENT(IN), OPTIONAL            :: debug_invert
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imap(3), k, l
    LOGICAL                                  :: dbg_inv, eval, failure
    REAL(KIND=dp)                            :: b, c2a, c2b, ca, ca2, cb, &
                                                cb2, check, pt5sq3, r2i, s2a, &
                                                s2b, sa, sb, sb2, sqb, sqb2i, &
                                                x11, x22, x33
    REAL(KIND=dp), DIMENSION(3) :: b_d, c2a_d, c2b_d, ca2_d, ca_d, cb2_d, &
      cb_d, r_d, s2a_d, s2b_d, sa_d, sb2_d, sb_d, sqb_d, x11_d, x22_d, x33_d
    REAL(KIND=dp), DIMENSION(3, 3)           :: p
    REAL(KIND=dp), DIMENSION(3, 3, 3)        :: p_d
    REAL(KIND=dp), DIMENSION(5, 5)           :: d
    REAL(KIND=dp), DIMENSION(5, 5, 3)        :: d_d

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(ij_matrix),cp_failure_level,routineP,error,failure)
    ! Compute Geomtric data and interatomic distance
    !     CA  = COS(PHI)    , SA  = SIN(PHI)
    !     CB  = COS(THETA)  , SB  = SIN(THETA)
    !     C2A = COS(2*PHI)  , S2A = SIN(2*PHI)
    !     C2B = COS(2*THETA), S2B = SIN(2*PHI)
    eval    = .TRUE.
    dbg_inv = .FALSE.
    pt5sq3  = 0.5_dp*SQRT(3.0_dp)
    imap(1) = 1
    imap(2) = 2
    imap(3) = 3
    check   = rjiv(3)/r
    IF (PRESENT(debug_invert)) dbg_inv = debug_invert
    ! Check for special case: both atoms lie on the Z Axis
    IF (ABS(check) > 0.99999999_dp) THEN
       IF (PRESENT(do_invert)) THEN
          imap(1)   = 1
          imap(2)   = 3
          imap(3)   = 2
          eval      = .TRUE.
          do_invert = .TRUE.
       ELSE
          sa      =  0.0_dp
          sb      =  0.0_dp
          IF      (check < 0.0_dp) THEN
             ca   = -1.0_dp
             cb   = -1.0_dp
          ELSE IF (check > 0.0_dp) THEN
             ca   =  1.0_dp
             cb   =  1.0_dp
          ELSE
             ca   =  0.0_dp
             cb   =  0.0_dp
          END IF
          eval    = .FALSE.
       END IF
    ELSE
       IF (PRESENT(do_invert)) do_invert = .FALSE.
    END IF
    IF (dbg_inv) THEN
       ! When debugging the derivatives of the rotation matrix we must possibly force
       ! the inversion of the reference frame
       CPPostcondition(.NOT.do_derivatives,cp_failure_level,routineP,error,failure)
       imap(1)   = 1
       imap(2)   = 3
       imap(3)   = 2
       eval      = .TRUE.
    END IF
    IF (eval) THEN
       x11   = rjiv(imap(1))
       x22   = rjiv(imap(2))
       x33   = rjiv(imap(3))
       cb    = x33/r
       b     = x11**2 + x22**2
       sqb   = SQRT(b)
       ca    = x11/sqb
       sa    = x22/sqb
       sb    = sqb/r
    END IF
    ! Calculate rotation matrix elements
    p(1, 1) =  ca * sb
    p(2, 1) =  ca * cb
    p(3, 1) = -sa
    p(1, 2) =  sa * sb
    p(2, 2) =  sa * cb
    p(3, 2) =  ca
    p(1, 3) =  cb
    p(2, 3) = -sb
    p(3, 3) =  0.0_dp
    IF (sepi%dorb.OR.sepj%dorb) THEN
       ca2 = ca**2
       cb2 = cb**2
       sb2 = sb**2
       c2a = 2.0_dp * ca2 - 1.0_dp
       c2b = 2.0_dp * cb2 - 1.0_dp
       s2a = 2.0_dp * sa * ca
       s2b = 2.0_dp * sb * cb
       d(1, 1) =  pt5sq3 * c2a * sb2
       d(2, 1) =  0.5_dp * c2a * s2b
       d(3, 1) = -s2a * sb
       d(4, 1) =  c2a * (cb2+0.5_dp*sb2)
       d(5, 1) = -s2a * cb
       d(1, 2) =  pt5sq3 * ca * s2b
       d(2, 2) =  ca * c2b
       d(3, 2) = -sa * cb
       d(4, 2) = -0.5_dp * ca * s2b
       d(5, 2) =  sa * sb
       d(1, 3) =  cb2 - 0.5_dp * sb2
       d(2, 3) = -pt5sq3 * s2b
       d(3, 3) =  0.0_dp
       d(4, 3) =  pt5sq3 * sb2
       d(5, 3) =  0.0_dp
       d(1, 4) =  pt5sq3 * sa * s2b
       d(2, 4) =  sa * c2b
       d(3, 4) =  ca * cb
       d(4, 4) = -0.5_dp * sa * s2b
       d(5, 4) = -ca * sb
       d(1, 5) =  pt5sq3 * s2a * sb2
       d(2, 5) =  0.5_dp * s2a * s2b
       d(3, 5) =  c2a * sb
       d(4, 5) =  s2a * (cb2+0.5_dp*sb2)
       d(5, 5) =  c2a * cb
    END IF
    !  Rotation Elements for : S-P
    DO k = 1, 3
       DO l = 1, 3
          ij_matrix%sp(k, l) = p(k, l)
       END DO
    END DO
    !  Rotation Elements for : P-P
    DO k = 1, 3
       ij_matrix%pp(1, k, k) = p(k, 1) * p(k, 1)
       ij_matrix%pp(2, k, k) = p(k, 2) * p(k, 2)
       ij_matrix%pp(3, k, k) = p(k, 3) * p(k, 3)
       ij_matrix%pp(4, k, k) = p(k, 1) * p(k, 2)
       ij_matrix%pp(5, k, k) = p(k, 1) * p(k, 3)
       ij_matrix%pp(6, k, k) = p(k, 2) * p(k, 3)
       IF (k /= 1) THEN
          DO l = 1, k - 1
             ij_matrix%pp(1, k, l) = 2.0_dp * p(k, 1) * p(l, 1)
             ij_matrix%pp(2, k, l) = 2.0_dp * p(k, 2) * p(l, 2)
             ij_matrix%pp(3, k, l) = 2.0_dp * p(k, 3) * p(l, 3)
             ij_matrix%pp(4, k, l) = p(k, 1) * p(l, 2) + p(k, 2) * p(l, 1)
             ij_matrix%pp(5, k, l) = p(k, 1) * p(l, 3) + p(k, 3) * p(l, 1)
             ij_matrix%pp(6, k, l) = p(k, 2) * p(l, 3) + p(k, 3) * p(l, 2)
          END DO
       END IF
    END DO
    IF (sepi%dorb.OR.sepj%dorb) THEN
       !  Rotation Elements for : S-D
       DO k = 1, 5
          DO l = 1, 5
             ij_matrix%sd(k, l) = d(k, l)
          END DO
       END DO
       !  Rotation Elements for : D-P
       DO k = 1, 5
          DO l = 1, 3
             ij_matrix%pd( 1, k, l) = d(k, 1) * p(l, 1)
             ij_matrix%pd( 2, k, l) = d(k, 1) * p(l, 2)
             ij_matrix%pd( 3, k, l) = d(k, 1) * p(l, 3)
             ij_matrix%pd( 4, k, l) = d(k, 2) * p(l, 1)
             ij_matrix%pd( 5, k, l) = d(k, 2) * p(l, 2)
             ij_matrix%pd( 6, k, l) = d(k, 2) * p(l, 3)
             ij_matrix%pd( 7, k, l) = d(k, 3) * p(l, 1)
             ij_matrix%pd( 8, k, l) = d(k, 3) * p(l, 2)
             ij_matrix%pd( 9, k, l) = d(k, 3) * p(l, 3)
             ij_matrix%pd(10, k, l) = d(k, 4) * p(l, 1)
             ij_matrix%pd(11, k, l) = d(k, 4) * p(l, 2)
             ij_matrix%pd(12, k, l) = d(k, 4) * p(l, 3)
             ij_matrix%pd(13, k, l) = d(k, 5) * p(l, 1)
             ij_matrix%pd(14, k, l) = d(k, 5) * p(l, 2)
             ij_matrix%pd(15, k, l) = d(k, 5) * p(l, 3)
          END DO
       END DO
       !  Rotation Elements for : D-D
       DO k = 1, 5
          ij_matrix%dd( 1, k, k) = d(k, 1) * d(k, 1)
          ij_matrix%dd( 2, k, k) = d(k, 2) * d(k, 2)
          ij_matrix%dd( 3, k, k) = d(k, 3) * d(k, 3)
          ij_matrix%dd( 4, k, k) = d(k, 4) * d(k, 4)
          ij_matrix%dd( 5, k, k) = d(k, 5) * d(k, 5)
          ij_matrix%dd( 6, k, k) = d(k, 1) * d(k, 2)
          ij_matrix%dd( 7, k, k) = d(k, 1) * d(k, 3)
          ij_matrix%dd( 8, k, k) = d(k, 2) * d(k, 3)
          ij_matrix%dd( 9, k, k) = d(k, 1) * d(k, 4)
          ij_matrix%dd(10, k, k) = d(k, 2) * d(k, 4)
          ij_matrix%dd(11, k, k) = d(k, 3) * d(k, 4)
          ij_matrix%dd(12, k, k) = d(k, 1) * d(k, 5)
          ij_matrix%dd(13, k, k) = d(k, 2) * d(k, 5)
          ij_matrix%dd(14, k, k) = d(k, 3) * d(k, 5)
          ij_matrix%dd(15, k, k) = d(k, 4) * d(k, 5)
          IF (k /= 1) THEN
             DO l = 1, k - 1
                ij_matrix%dd( 1, k, l) = 2.0_dp  * d(k, 1) * d(l, 1)
                ij_matrix%dd( 2, k, l) = 2.0_dp  * d(k, 2) * d(l, 2)
                ij_matrix%dd( 3, k, l) = 2.0_dp  * d(k, 3) * d(l, 3)
                ij_matrix%dd( 4, k, l) = 2.0_dp  * d(k, 4) * d(l, 4)
                ij_matrix%dd( 5, k, l) = 2.0_dp  * d(k, 5) * d(l, 5)
                ij_matrix%dd( 6, k, l) = d(k, 1) * d(l, 2) + d(k, 2) * d(l, 1)
                ij_matrix%dd( 7, k, l) = d(k, 1) * d(l, 3) + d(k, 3) * d(l, 1)
                ij_matrix%dd( 8, k, l) = d(k, 2) * d(l, 3) + d(k, 3) * d(l, 2)
                ij_matrix%dd( 9, k, l) = d(k, 1) * d(l, 4) + d(k, 4) * d(l, 1)
                ij_matrix%dd(10, k, l) = d(k, 2) * d(l, 4) + d(k, 4) * d(l, 2)
                ij_matrix%dd(11, k, l) = d(k, 3) * d(l, 4) + d(k, 4) * d(l, 3)
                ij_matrix%dd(12, k, l) = d(k, 1) * d(l, 5) + d(k, 5) * d(l, 1)
                ij_matrix%dd(13, k, l) = d(k, 2) * d(l, 5) + d(k, 5) * d(l, 2)
                ij_matrix%dd(14, k, l) = d(k, 3) * d(l, 5) + d(k, 5) * d(l, 3)
                ij_matrix%dd(15, k, l) = d(k, 4) * d(l, 5) + d(k, 5) * d(l, 4)
             END DO
          END IF
       END DO
    ELSE
       ij_matrix%sd = 0.0_dp
       ij_matrix%pd = 0.0_dp
       ij_matrix%dd = 0.0_dp
    END IF
    ! Evaluate Derivatives if required
    IF (do_derivatives) THEN
       ! This condition is necessary because derivatives uses the invertion of the
       ! axis for treating the divergence point along z-axis
       CPPostcondition(eval,cp_failure_level,routineP,error,failure)
       x11_d = 0.0_dp; x11_d(1) = -1.0_dp
       x22_d = 0.0_dp; x22_d(2) = -1.0_dp
       x33_d = 0.0_dp; x33_d(3) = -1.0_dp
       r_d   = -(/x11,x22,x33/)/r
       b_d   = 2.0_dp*(x11*x11_d + x22*x22_d)
       sqb_d = (0.5_dp/sqb)*b_d
       r2i   = 1.0_dp/r**2
       sqb2i = 1.0_dp/sqb**2
       cb_d  = (x33_d*r-x33*r_d)*r2i
       ca_d  = (x11_d*sqb-x11*sqb_d)*sqb2i
       sa_d  = (x22_d*sqb-x22*sqb_d)*sqb2i
       sb_d  = (sqb_d*r-sqb*r_d)*r2i
       ! Calculate derivatives of rotation matrix elements
       p_d(1, 1, :) =  ca_d * sb + ca * sb_d
       p_d(2, 1, :) =  ca_d * cb + ca * cb_d
       p_d(3, 1, :) = -sa_d
       p_d(1, 2, :) =  sa_d * sb + sa * sb_d
       p_d(2, 2, :) =  sa_d * cb + sa * cb_d
       p_d(3, 2, :) =  ca_d
       p_d(1, 3, :) =  cb_d
       p_d(2, 3, :) = -sb_d
       p_d(3, 3, :) =  0.0_dp
       IF (sepi%dorb.OR.sepj%dorb) THEN
          ca2_d = 2.0_dp * ca * ca_d
          cb2_d = 2.0_dp * cb * cb_d
          sb2_d = 2.0_dp * sb * sb_d
          c2a_d = 2.0_dp * ca2_d
          c2b_d = 2.0_dp * cb2_d
          s2a_d = 2.0_dp * ( sa_d * ca + sa * ca_d )
          s2b_d = 2.0_dp * ( sb_d * cb + sb * cb_d )
          d_d(1, 1, :) =  pt5sq3 * ( c2a_d * sb2 + c2a * sb2_d )
          d_d(2, 1, :) =  0.5_dp * ( c2a_d * s2b + c2a * s2b_d )
          d_d(3, 1, :) = -s2a_d * sb - s2a * sb_d
          d_d(4, 1, :) =  c2a_d * (cb2+0.5_dp*sb2) + c2a * (cb2_d+0.5_dp*sb2_d)
          d_d(5, 1, :) = -s2a_d * cb - s2a * cb_d
          d_d(1, 2, :) =  pt5sq3 * ( ca_d * s2b + ca * s2b_d )
          d_d(2, 2, :) =  ca_d * c2b + ca * c2b_d
          d_d(3, 2, :) = -sa_d * cb - sa * cb_d
          d_d(4, 2, :) = -0.5_dp * ( ca_d * s2b + ca * s2b_d ) 
          d_d(5, 2, :) =  sa_d * sb + sa * sb_d
          d_d(1, 3, :) =  cb2_d - 0.5_dp * sb2_d
          d_d(2, 3, :) = -pt5sq3 * s2b_d
          d_d(3, 3, :) =  0.0_dp
          d_d(4, 3, :) =  pt5sq3 * sb2_d
          d_d(5, 3, :) =  0.0_dp
          d_d(1, 4, :) =  pt5sq3 * ( sa_d * s2b + sa * s2b_d )  
          d_d(2, 4, :) =  sa_d * c2b + sa * c2b_d
          d_d(3, 4, :) =  ca_d * cb + ca * cb_d 
          d_d(4, 4, :) = -0.5_dp * ( sa_d * s2b + sa * s2b_d )
          d_d(5, 4, :) = -ca_d * sb -ca * sb_d
          d_d(1, 5, :) =  pt5sq3 * ( s2a_d * sb2 + s2a * sb2_d )
          d_d(2, 5, :) =  0.5_dp * ( s2a_d * s2b + s2a * s2b_d )
          d_d(3, 5, :) =  c2a_d * sb + c2a * sb_d 
          d_d(4, 5, :) =  s2a_d * (cb2+0.5_dp*sb2) + s2a * (cb2_d+0.5_dp*sb2_d)
          d_d(5, 5, :) =  c2a_d * cb + c2a * cb_d 
       END IF
       !  Derivative for Rotation Elements for : S-P
       DO k = 1, 3
          DO l = 1, 3
             ij_matrix%sp_d(k, l, :) = p_d(k, l, :)
          END DO
       END DO
       !  Derivative for Rotation Elements for : P-P
       DO k = 1, 3
          ij_matrix%pp_d(1, k, k, :) = p_d(k, 1, :) * p(k, 1) +  p(k, 1) * p_d(k, 1, :)
          ij_matrix%pp_d(2, k, k, :) = p_d(k, 2, :) * p(k, 2) +  p(k, 2) * p_d(k, 2, :)
          ij_matrix%pp_d(3, k, k, :) = p_d(k, 3, :) * p(k, 3) +  p(k, 3) * p_d(k, 3, :)
          ij_matrix%pp_d(4, k, k, :) = p_d(k, 1, :) * p(k, 2) +  p(k, 1) * p_d(k, 2, :)
          ij_matrix%pp_d(5, k, k, :) = p_d(k, 1, :) * p(k, 3) +  p(k, 1) * p_d(k, 3, :)
          ij_matrix%pp_d(6, k, k, :) = p_d(k, 2, :) * p(k, 3) +  p(k, 2) * p_d(k, 3, :)
          IF (k /= 1) THEN
             DO l = 1, k - 1
                ij_matrix%pp_d(1, k, l, :) = 2.0_dp * ( p_d(k, 1, :) * p(l, 1) + p(k, 1) * p_d(l, 1, :) ) 
                ij_matrix%pp_d(2, k, l, :) = 2.0_dp * ( p_d(k, 2, :) * p(l, 2) + p(k, 2) * p_d(l, 2, :) )
                ij_matrix%pp_d(3, k, l, :) = 2.0_dp * ( p_d(k, 3, :) * p(l, 3) + p(k, 3) * p_d(l, 3, :) )
                ij_matrix%pp_d(4, k, l, :) = ( p_d(k, 1, :) * p(l, 2) + p(k, 1) * p_d(l, 2, :) ) + &
                                             ( p_d(k, 2, :) * p(l, 1) + p(k, 2) * p_d(l, 1, :) ) 
                ij_matrix%pp_d(5, k, l, :) = ( p_d(k, 1, :) * p(l, 3) + p(k, 1) * p_d(l, 3, :) ) + &
                                             ( p_d(k, 3, :) * p(l, 1) + p(k, 3) * p_d(l, 1, :) ) 
                ij_matrix%pp_d(6, k, l, :) = ( p_d(k, 2, :) * p(l, 3) + p(k, 2) * p_d(l, 3, :) ) + &
                                             ( p_d(k, 3, :) * p(l, 2) + p(k, 3) * p_d(l, 2, :) ) 
             END DO
          END IF
       END DO
       IF (sepi%dorb.OR.sepj%dorb) THEN
          !  Rotation Elements for : S-D
          DO k = 1, 5
             DO l = 1, 5
                ij_matrix%sd_d(k, l, :) = d_d(k, l, :)
             END DO
          END DO
          !  Rotation Elements for : D-P
          DO k = 1, 5
             DO l = 1, 3
                ij_matrix%pd_d( 1, k, l, :) = ( d_d(k, 1, :) * p(l, 1) + d(k, 1) * p_d(l, 1, :) )
                ij_matrix%pd_d( 2, k, l, :) = ( d_d(k, 1, :) * p(l, 2) + d(k, 1) * p_d(l, 2, :) )
                ij_matrix%pd_d( 3, k, l, :) = ( d_d(k, 1, :) * p(l, 3) + d(k, 1) * p_d(l, 3, :) )
                ij_matrix%pd_d( 4, k, l, :) = ( d_d(k, 2, :) * p(l, 1) + d(k, 2) * p_d(l, 1, :) )
                ij_matrix%pd_d( 5, k, l, :) = ( d_d(k, 2, :) * p(l, 2) + d(k, 2) * p_d(l, 2, :) )
                ij_matrix%pd_d( 6, k, l, :) = ( d_d(k, 2, :) * p(l, 3) + d(k, 2) * p_d(l, 3, :) )
                ij_matrix%pd_d( 7, k, l, :) = ( d_d(k, 3, :) * p(l, 1) + d(k, 3) * p_d(l, 1, :) )
                ij_matrix%pd_d( 8, k, l, :) = ( d_d(k, 3, :) * p(l, 2) + d(k, 3) * p_d(l, 2, :) )
                ij_matrix%pd_d( 9, k, l, :) = ( d_d(k, 3, :) * p(l, 3) + d(k, 3) * p_d(l, 3, :) )
                ij_matrix%pd_d(10, k, l, :) = ( d_d(k, 4, :) * p(l, 1) + d(k, 4) * p_d(l, 1, :) )
                ij_matrix%pd_d(11, k, l, :) = ( d_d(k, 4, :) * p(l, 2) + d(k, 4) * p_d(l, 2, :) )
                ij_matrix%pd_d(12, k, l, :) = ( d_d(k, 4, :) * p(l, 3) + d(k, 4) * p_d(l, 3, :) )
                ij_matrix%pd_d(13, k, l, :) = ( d_d(k, 5, :) * p(l, 1) + d(k, 5) * p_d(l, 1, :) )
                ij_matrix%pd_d(14, k, l, :) = ( d_d(k, 5, :) * p(l, 2) + d(k, 5) * p_d(l, 2, :) )
                ij_matrix%pd_d(15, k, l, :) = ( d_d(k, 5, :) * p(l, 3) + d(k, 5) * p_d(l, 3, :) )
             END DO
          END DO
          !  Rotation Elements for : D-D
          DO k = 1, 5
             ij_matrix%dd_d( 1, k, k, :) = ( d_d(k, 1, :) * d(k, 1) + d(k, 1) * d_d(k, 1, :) ) 
             ij_matrix%dd_d( 2, k, k, :) = ( d_d(k, 2, :) * d(k, 2) + d(k, 2) * d_d(k, 2, :) )
             ij_matrix%dd_d( 3, k, k, :) = ( d_d(k, 3, :) * d(k, 3) + d(k, 3) * d_d(k, 3, :) )
             ij_matrix%dd_d( 4, k, k, :) = ( d_d(k, 4, :) * d(k, 4) + d(k, 4) * d_d(k, 4, :) )
             ij_matrix%dd_d( 5, k, k, :) = ( d_d(k, 5, :) * d(k, 5) + d(k, 5) * d_d(k, 5, :) )
             ij_matrix%dd_d( 6, k, k, :) = ( d_d(k, 1, :) * d(k, 2) + d(k, 1) * d_d(k, 2, :) )
             ij_matrix%dd_d( 7, k, k, :) = ( d_d(k, 1, :) * d(k, 3) + d(k, 1) * d_d(k, 3, :) )
             ij_matrix%dd_d( 8, k, k, :) = ( d_d(k, 2, :) * d(k, 3) + d(k, 2) * d_d(k, 3, :) )
             ij_matrix%dd_d( 9, k, k, :) = ( d_d(k, 1, :) * d(k, 4) + d(k, 1) * d_d(k, 4, :) )
             ij_matrix%dd_d(10, k, k, :) = ( d_d(k, 2, :) * d(k, 4) + d(k, 2) * d_d(k, 4, :) )
             ij_matrix%dd_d(11, k, k, :) = ( d_d(k, 3, :) * d(k, 4) + d(k, 3) * d_d(k, 4, :) )
             ij_matrix%dd_d(12, k, k, :) = ( d_d(k, 1, :) * d(k, 5) + d(k, 1) * d_d(k, 5, :) )
             ij_matrix%dd_d(13, k, k, :) = ( d_d(k, 2, :) * d(k, 5) + d(k, 2) * d_d(k, 5, :) )
             ij_matrix%dd_d(14, k, k, :) = ( d_d(k, 3, :) * d(k, 5) + d(k, 3) * d_d(k, 5, :) )
             ij_matrix%dd_d(15, k, k, :) = ( d_d(k, 4, :) * d(k, 5) + d(k, 4) * d_d(k, 5, :) )
             IF (k /= 1) THEN
                DO l = 1, k - 1
                   ij_matrix%dd_d( 1, k, l, :) = 2.0_dp  * ( d_d(k, 1, :) * d(l, 1) + d(k, 1) * d_d(l, 1, :) )
                   ij_matrix%dd_d( 2, k, l, :) = 2.0_dp  * ( d_d(k, 2, :) * d(l, 2) + d(k, 2) * d_d(l, 2, :) )
                   ij_matrix%dd_d( 3, k, l, :) = 2.0_dp  * ( d_d(k, 3, :) * d(l, 3) + d(k, 3) * d_d(l, 3, :) )
                   ij_matrix%dd_d( 4, k, l, :) = 2.0_dp  * ( d_d(k, 4, :) * d(l, 4) + d(k, 4) * d_d(l, 4, :) )
                   ij_matrix%dd_d( 5, k, l, :) = 2.0_dp  * ( d_d(k, 5, :) * d(l, 5) + d(k, 5) * d_d(l, 5, :) )
                   ij_matrix%dd_d( 6, k, l, :) = ( d_d(k, 1, :) * d(l, 2) + d(k, 1) * d_d(l, 2, :) ) + &
                                                 ( d_d(k, 2, :) * d(l, 1) + d(k, 2) * d_d(l, 1, :) )
                   ij_matrix%dd_d( 7, k, l, :) = ( d_d(k, 1, :) * d(l, 3) + d(k, 1) * d_d(l, 3, :) ) + &
                                                 ( d_d(k, 3, :) * d(l, 1) + d(k, 3) * d_d(l, 1, :) )
                   ij_matrix%dd_d( 8, k, l, :) = ( d_d(k, 2, :) * d(l, 3) + d(k, 2) * d_d(l, 3, :) ) + &
                                                 ( d_d(k, 3, :) * d(l, 2) + d(k, 3) * d_d(l, 2, :) )
                   ij_matrix%dd_d( 9, k, l, :) = ( d_d(k, 1, :) * d(l, 4) + d(k, 1) * d_d(l, 4, :) ) + &
                                                 ( d_d(k, 4, :) * d(l, 1) + d(k, 4) * d_d(l, 1, :) )
                   ij_matrix%dd_d(10, k, l, :) = ( d_d(k, 2, :) * d(l, 4) + d(k, 2) * d_d(l, 4, :) ) + &
                                                 ( d_d(k, 4, :) * d(l, 2) + d(k, 4) * d_d(l, 2, :) )
                   ij_matrix%dd_d(11, k, l, :) = ( d_d(k, 3, :) * d(l, 4) + d(k, 3) * d_d(l, 4, :) ) + &
                                                 ( d_d(k, 4, :) * d(l, 3) + d(k, 4) * d_d(l, 3, :) )
                   ij_matrix%dd_d(12, k, l, :) = ( d_d(k, 1, :) * d(l, 5) + d(k, 1) * d_d(l, 5, :) ) + &
                                                 ( d_d(k, 5, :) * d(l, 1) + d(k, 5) * d_d(l, 1, :) )
                   ij_matrix%dd_d(13, k, l, :) = ( d_d(k, 2, :) * d(l, 5) + d(k, 2) * d_d(l, 5, :) ) + &
                                                 ( d_d(k, 5, :) * d(l, 2) + d(k, 5) * d_d(l, 2, :) )
                   ij_matrix%dd_d(14, k, l, :) = ( d_d(k, 3, :) * d(l, 5) + d(k, 3) * d_d(l, 5, :) ) + &
                                                 ( d_d(k, 5, :) * d(l, 3) + d(k, 5) * d_d(l, 3, :) )
                   ij_matrix%dd_d(15, k, l, :) = ( d_d(k, 4, :) * d(l, 5) + d(k, 4) * d_d(l, 5, :) ) + &
                                                 ( d_d(k, 5, :) * d(l, 4) + d(k, 5) * d_d(l, 4, :) )
                END DO
             END IF
          END DO
       ELSE
          ij_matrix%sd_d = 0.0_dp
          ij_matrix%pd_d = 0.0_dp
          ij_matrix%dd_d = 0.0_dp
       END IF
       IF (debug_this_module) THEN 
          CALL check_rotmat_der( sepi, sepj, -rjiv, ij_matrix, do_invert=do_invert, error=error)
       END IF
    END IF

  END SUBROUTINE rotmat

! *****************************************************************************
!> \brief First Step of the rotation of the two-electron two-center integrals in 
!>        SPD basis
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  RECURSIVE SUBROUTINE rot_2el_2c_first (sepi, sepj, rijv, shortrange, taper, invert, ii, kk,&
       rep, logv, ij_matrix, v, lgrad, rep_d, v_d, logv_d, drij, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rijv
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(taper_type), POINTER                :: taper
    LOGICAL, INTENT(IN)                      :: invert
    INTEGER, INTENT(IN)                      :: ii, kk
    REAL(KIND=dp), DIMENSION(491), &
      INTENT(IN)                             :: rep
    LOGICAL, DIMENSION(45, 45), INTENT(OUT)  :: logv
    TYPE(rotmat_type), POINTER               :: ij_matrix
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(OUT)                            :: v
    LOGICAL, INTENT(IN)                      :: lgrad
    REAL(KIND=dp), DIMENSION(491), &
      INTENT(IN), OPTIONAL                   :: rep_d
    REAL(KIND=dp), DIMENSION(45, 45, 3), &
      INTENT(OUT), OPTIONAL                  :: v_d
    LOGICAL, DIMENSION(45, 45), &
      INTENT(OUT), OPTIONAL                  :: logv_d
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(IN), OPTIONAL                   :: drij
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rot_2el_2c_first', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, ij, j1, k, k1, kl, l, &
                                                l1, limkl, ll, mm, nd
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: wrepp, wrepp_d(3)

    failure = .FALSE.
    IF (lgrad) THEN
       CPPostcondition(PRESENT( rep_d),cp_failure_level,routineP,error,failure)
       CPPostcondition(PRESENT(   v_d),cp_failure_level,routineP,error,failure)
       CPPostcondition(PRESENT(logv_d),cp_failure_level,routineP,error,failure)
       CPPostcondition(PRESENT(  drij),cp_failure_level,routineP,error,failure)
    END IF
    limkl = indexb(kk, kk)
    DO i = 1, 45
       DO k = 1, limkl
          logv(i, k) = .FALSE.
          v(i, k) = 0.0_dp
       END DO
    END DO
    !
    DO i1 = 1, ii
       DO j1 = 1, i1
          ij = indexa(i1, j1)
          !
          DO k1 = 1, kk
             !
             DO l1 = 1, k1
                kl = indexa(k1, l1)
                nd = ijkl_ind(ij, kl)
                IF (nd /= 0) THEN
                   !
                   wrepp = rep(nd)
                   ll = indexb(k1, l1)
                   mm = int2c_type(ll)
                   !
                   IF (mm == 1) THEN
                      v(ij, 1) = wrepp
                   ELSE IF (mm == 2) THEN
                      k = k1 - 1
                      v(ij,  2) = v(ij,  2) + ij_matrix%sp(k, 1) * wrepp
                      v(ij,  4) = v(ij,  4) + ij_matrix%sp(k, 2) * wrepp
                      v(ij,  7) = v(ij,  7) + ij_matrix%sp(k, 3) * wrepp
                   ELSE IF (mm == 3) THEN
                      k = k1 - 1
                      l = l1 - 1
                      v(ij,  3) = v(ij,  3) + ij_matrix%pp(1, k, l) * wrepp
                      v(ij,  6) = v(ij,  6) + ij_matrix%pp(2, k, l) * wrepp
                      v(ij, 10) = v(ij, 10) + ij_matrix%pp(3, k, l) * wrepp
                      v(ij,  5) = v(ij,  5) + ij_matrix%pp(4, k, l) * wrepp
                      v(ij,  8) = v(ij,  8) + ij_matrix%pp(5, k, l) * wrepp
                      v(ij,  9) = v(ij,  9) + ij_matrix%pp(6, k, l) * wrepp
                   ELSE IF (mm == 4) THEN
                      k = k1 - 4
                      v(ij, 11) = v(ij, 11) + ij_matrix%sd(k, 1) * wrepp
                      v(ij, 16) = v(ij, 16) + ij_matrix%sd(k, 2) * wrepp
                      v(ij, 22) = v(ij, 22) + ij_matrix%sd(k, 3) * wrepp
                      v(ij, 29) = v(ij, 29) + ij_matrix%sd(k, 4) * wrepp
                      v(ij, 37) = v(ij, 37) + ij_matrix%sd(k, 5) * wrepp
                   ELSE IF (mm == 5) THEN
                      k = k1 - 4
                      l = l1 - 1
                      v(ij, 12) = v(ij, 12) + ij_matrix%pd( 1, k, l) * wrepp
                      v(ij, 13) = v(ij, 13) + ij_matrix%pd( 2, k, l) * wrepp
                      v(ij, 14) = v(ij, 14) + ij_matrix%pd( 3, k, l) * wrepp
                      v(ij, 17) = v(ij, 17) + ij_matrix%pd( 4, k, l) * wrepp
                      v(ij, 18) = v(ij, 18) + ij_matrix%pd( 5, k, l) * wrepp
                      v(ij, 19) = v(ij, 19) + ij_matrix%pd( 6, k, l) * wrepp
                      v(ij, 23) = v(ij, 23) + ij_matrix%pd( 7, k, l) * wrepp
                      v(ij, 24) = v(ij, 24) + ij_matrix%pd( 8, k, l) * wrepp
                      v(ij, 25) = v(ij, 25) + ij_matrix%pd( 9, k, l) * wrepp
                      v(ij, 30) = v(ij, 30) + ij_matrix%pd(10, k, l) * wrepp
                      v(ij, 31) = v(ij, 31) + ij_matrix%pd(11, k, l) * wrepp
                      v(ij, 32) = v(ij, 32) + ij_matrix%pd(12, k, l) * wrepp
                      v(ij, 38) = v(ij, 38) + ij_matrix%pd(13, k, l) * wrepp
                      v(ij, 39) = v(ij, 39) + ij_matrix%pd(14, k, l) * wrepp
                      v(ij, 40) = v(ij, 40) + ij_matrix%pd(15, k, l) * wrepp
                   ELSE IF (mm == 6) THEN
                      k = k1 - 4
                      l = l1 - 4
                      v(ij, 15) = v(ij, 15) + ij_matrix%dd( 1, k, l) * wrepp
                      v(ij, 21) = v(ij, 21) + ij_matrix%dd( 2, k, l) * wrepp
                      v(ij, 28) = v(ij, 28) + ij_matrix%dd( 3, k, l) * wrepp
                      v(ij, 36) = v(ij, 36) + ij_matrix%dd( 4, k, l) * wrepp
                      v(ij, 45) = v(ij, 45) + ij_matrix%dd( 5, k, l) * wrepp
                      v(ij, 20) = v(ij, 20) + ij_matrix%dd( 6, k, l) * wrepp
                      v(ij, 26) = v(ij, 26) + ij_matrix%dd( 7, k, l) * wrepp
                      v(ij, 27) = v(ij, 27) + ij_matrix%dd( 8, k, l) * wrepp
                      v(ij, 33) = v(ij, 33) + ij_matrix%dd( 9, k, l) * wrepp
                      v(ij, 34) = v(ij, 34) + ij_matrix%dd(10, k, l) * wrepp
                      v(ij, 35) = v(ij, 35) + ij_matrix%dd(11, k, l) * wrepp
                      v(ij, 41) = v(ij, 41) + ij_matrix%dd(12, k, l) * wrepp
                      v(ij, 42) = v(ij, 42) + ij_matrix%dd(13, k, l) * wrepp
                      v(ij, 43) = v(ij, 43) + ij_matrix%dd(14, k, l) * wrepp
                      v(ij, 44) = v(ij, 44) + ij_matrix%dd(15, k, l) * wrepp
                   END IF
                END IF
             END DO
          END DO
          DO kl = 1, limkl
             logv(ij, kl) = (ABS (v(ij, kl)) > 0.0_dp)
          END DO
       END DO
    END DO
    ! Gradients
    IF (lgrad) THEN
       DO i = 1, 45
          DO k = 1, limkl
             logv_d(i, k) = .FALSE.
             v_d(i, k, :) = 0.0_dp
          END DO
       END DO
       DO i1 = 1, ii
          DO j1 = 1, i1
             ij = indexa(i1, j1)
             !
             DO k1 = 1, kk
                !
                DO l1 = 1, k1
                   kl = indexa(k1, l1)
                   nd = ijkl_ind(ij, kl)
                   IF (nd /= 0) THEN
                      !
                      wrepp   = rep(nd)
                      wrepp_d = rep_d(nd)*drij
                      ll = indexb(k1, l1)
                      mm = int2c_type(ll)
                      !
                      IF (mm == 1) THEN
                         v_d(ij, 1,:) = wrepp_d
                      ELSE IF (mm == 2) THEN
                         k = k1 - 1
                         v_d(ij,  2, :) = v_d(ij,  2, :) + ij_matrix%sp_d(k, 1, :) * wrepp + ij_matrix%sp(k, 1) * wrepp_d
                         v_d(ij,  4, :) = v_d(ij,  4, :) + ij_matrix%sp_d(k, 2, :) * wrepp + ij_matrix%sp(k, 2) * wrepp_d
                         v_d(ij,  7, :) = v_d(ij,  7, :) + ij_matrix%sp_d(k, 3, :) * wrepp + ij_matrix%sp(k, 3) * wrepp_d
                      ELSE IF (mm == 3) THEN
                         k = k1 - 1
                         l = l1 - 1
                         v_d(ij,  3, :) = v_d(ij,  3, :) + ij_matrix%pp_d(1, k, l, :) * wrepp + ij_matrix%pp(1, k, l) * wrepp_d
                         v_d(ij,  6, :) = v_d(ij,  6, :) + ij_matrix%pp_d(2, k, l, :) * wrepp + ij_matrix%pp(2, k, l) * wrepp_d
                         v_d(ij, 10, :) = v_d(ij, 10, :) + ij_matrix%pp_d(3, k, l, :) * wrepp + ij_matrix%pp(3, k, l) * wrepp_d
                         v_d(ij,  5, :) = v_d(ij,  5, :) + ij_matrix%pp_d(4, k, l, :) * wrepp + ij_matrix%pp(4, k, l) * wrepp_d
                         v_d(ij,  8, :) = v_d(ij,  8, :) + ij_matrix%pp_d(5, k, l, :) * wrepp + ij_matrix%pp(5, k, l) * wrepp_d
                         v_d(ij,  9, :) = v_d(ij,  9, :) + ij_matrix%pp_d(6, k, l, :) * wrepp + ij_matrix%pp(6, k, l) * wrepp_d
                      ELSE IF (mm == 4) THEN
                         k = k1 - 4
                         v_d(ij, 11, :) = v_d(ij, 11, :) + ij_matrix%sd_d(k, 1, :) * wrepp + ij_matrix%sd(k, 1) * wrepp_d
                         v_d(ij, 16, :) = v_d(ij, 16, :) + ij_matrix%sd_d(k, 2, :) * wrepp + ij_matrix%sd(k, 2) * wrepp_d
                         v_d(ij, 22, :) = v_d(ij, 22, :) + ij_matrix%sd_d(k, 3, :) * wrepp + ij_matrix%sd(k, 3) * wrepp_d
                         v_d(ij, 29, :) = v_d(ij, 29, :) + ij_matrix%sd_d(k, 4, :) * wrepp + ij_matrix%sd(k, 4) * wrepp_d
                         v_d(ij, 37, :) = v_d(ij, 37, :) + ij_matrix%sd_d(k, 5, :) * wrepp + ij_matrix%sd(k, 5) * wrepp_d
                      ELSE IF (mm == 5) THEN
                         k = k1 - 4
                         l = l1 - 1
                         v_d(ij, 12, :) = v_d(ij, 12, :) + ij_matrix%pd_d( 1, k, l, :) * wrepp + ij_matrix%pd( 1, k, l) * wrepp_d
                         v_d(ij, 13, :) = v_d(ij, 13, :) + ij_matrix%pd_d( 2, k, l, :) * wrepp + ij_matrix%pd( 2, k, l) * wrepp_d
                         v_d(ij, 14, :) = v_d(ij, 14, :) + ij_matrix%pd_d( 3, k, l, :) * wrepp + ij_matrix%pd( 3, k, l) * wrepp_d
                         v_d(ij, 17, :) = v_d(ij, 17, :) + ij_matrix%pd_d( 4, k, l, :) * wrepp + ij_matrix%pd( 4, k, l) * wrepp_d
                         v_d(ij, 18, :) = v_d(ij, 18, :) + ij_matrix%pd_d( 5, k, l, :) * wrepp + ij_matrix%pd( 5, k, l) * wrepp_d
                         v_d(ij, 19, :) = v_d(ij, 19, :) + ij_matrix%pd_d( 6, k, l, :) * wrepp + ij_matrix%pd( 6, k, l) * wrepp_d
                         v_d(ij, 23, :) = v_d(ij, 23, :) + ij_matrix%pd_d( 7, k, l, :) * wrepp + ij_matrix%pd( 7, k, l) * wrepp_d
                         v_d(ij, 24, :) = v_d(ij, 24, :) + ij_matrix%pd_d( 8, k, l, :) * wrepp + ij_matrix%pd( 8, k, l) * wrepp_d
                         v_d(ij, 25, :) = v_d(ij, 25, :) + ij_matrix%pd_d( 9, k, l, :) * wrepp + ij_matrix%pd( 9, k, l) * wrepp_d
                         v_d(ij, 30, :) = v_d(ij, 30, :) + ij_matrix%pd_d(10, k, l, :) * wrepp + ij_matrix%pd(10, k, l) * wrepp_d
                         v_d(ij, 31, :) = v_d(ij, 31, :) + ij_matrix%pd_d(11, k, l, :) * wrepp + ij_matrix%pd(11, k, l) * wrepp_d
                         v_d(ij, 32, :) = v_d(ij, 32, :) + ij_matrix%pd_d(12, k, l, :) * wrepp + ij_matrix%pd(12, k, l) * wrepp_d
                         v_d(ij, 38, :) = v_d(ij, 38, :) + ij_matrix%pd_d(13, k, l, :) * wrepp + ij_matrix%pd(13, k, l) * wrepp_d
                         v_d(ij, 39, :) = v_d(ij, 39, :) + ij_matrix%pd_d(14, k, l, :) * wrepp + ij_matrix%pd(14, k, l) * wrepp_d
                         v_d(ij, 40, :) = v_d(ij, 40, :) + ij_matrix%pd_d(15, k, l, :) * wrepp + ij_matrix%pd(15, k, l) * wrepp_d
                      ELSE IF (mm == 6) THEN
                         k = k1 - 4
                         l = l1 - 4
                         v_d(ij, 15, :) = v_d(ij, 15, :) + ij_matrix%dd_d( 1, k, l, :) * wrepp + ij_matrix%dd( 1, k, l) * wrepp_d
                         v_d(ij, 21, :) = v_d(ij, 21, :) + ij_matrix%dd_d( 2, k, l, :) * wrepp + ij_matrix%dd( 2, k, l) * wrepp_d
                         v_d(ij, 28, :) = v_d(ij, 28, :) + ij_matrix%dd_d( 3, k, l, :) * wrepp + ij_matrix%dd( 3, k, l) * wrepp_d
                         v_d(ij, 36, :) = v_d(ij, 36, :) + ij_matrix%dd_d( 4, k, l, :) * wrepp + ij_matrix%dd( 4, k, l) * wrepp_d
                         v_d(ij, 45, :) = v_d(ij, 45, :) + ij_matrix%dd_d( 5, k, l, :) * wrepp + ij_matrix%dd( 5, k, l) * wrepp_d
                         v_d(ij, 20, :) = v_d(ij, 20, :) + ij_matrix%dd_d( 6, k, l, :) * wrepp + ij_matrix%dd( 6, k, l) * wrepp_d
                         v_d(ij, 26, :) = v_d(ij, 26, :) + ij_matrix%dd_d( 7, k, l, :) * wrepp + ij_matrix%dd( 7, k, l) * wrepp_d
                         v_d(ij, 27, :) = v_d(ij, 27, :) + ij_matrix%dd_d( 8, k, l, :) * wrepp + ij_matrix%dd( 8, k, l) * wrepp_d
                         v_d(ij, 33, :) = v_d(ij, 33, :) + ij_matrix%dd_d( 9, k, l, :) * wrepp + ij_matrix%dd( 9, k, l) * wrepp_d
                         v_d(ij, 34, :) = v_d(ij, 34, :) + ij_matrix%dd_d(10, k, l, :) * wrepp + ij_matrix%dd(10, k, l) * wrepp_d
                         v_d(ij, 35, :) = v_d(ij, 35, :) + ij_matrix%dd_d(11, k, l, :) * wrepp + ij_matrix%dd(11, k, l) * wrepp_d
                         v_d(ij, 41, :) = v_d(ij, 41, :) + ij_matrix%dd_d(12, k, l, :) * wrepp + ij_matrix%dd(12, k, l) * wrepp_d
                         v_d(ij, 42, :) = v_d(ij, 42, :) + ij_matrix%dd_d(13, k, l, :) * wrepp + ij_matrix%dd(13, k, l) * wrepp_d
                         v_d(ij, 43, :) = v_d(ij, 43, :) + ij_matrix%dd_d(14, k, l, :) * wrepp + ij_matrix%dd(14, k, l) * wrepp_d
                         v_d(ij, 44, :) = v_d(ij, 44, :) + ij_matrix%dd_d(15, k, l, :) * wrepp + ij_matrix%dd(15, k, l) * wrepp_d
                      END IF
                   END IF
                END DO
             END DO
             DO kl = 1, limkl
                logv_d(ij, kl) = (ANY(ABS(v_d(ij, kl, :)) > 0.0_dp))
             END DO
          END DO
       END DO
       IF (debug_this_module) THEN
          CALL rot_2el_2c_first_debug(sepi, sepj, rijv, shortrange, taper, invert, ii, kk, v_d, error)
       END IF
    END IF
  END SUBROUTINE rot_2el_2c_first

! *****************************************************************************
!> \brief Store the two-electron two-center integrals in diagonal form
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE store_2el_2c_diag (limij, limkl, ww, w, ww_dx, ww_dy, ww_dz, dw)

    INTEGER, INTENT(IN)                      :: limij, limkl
    REAL(KIND=dp), DIMENSION(limkl, limij), &
      INTENT(IN), OPTIONAL                   :: ww
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: w
    REAL(KIND=dp), DIMENSION(limkl, limij), &
      INTENT(IN), OPTIONAL                   :: ww_dx, ww_dy, ww_dz
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: dw

    CHARACTER(len=*), PARAMETER :: routineN = 'store_2el_2c_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ij, kl, l

    l = 0
    IF      (PRESENT(ww).AND.PRESENT(w)) THEN
       DO ij = 1, limij
          DO kl = 1, limkl
             l = l + 1
             w(l) = ww(kl, ij)
          END DO
       END DO
    ELSE IF (PRESENT(ww_dx).AND.PRESENT(ww_dy).AND.PRESENT(ww_dz).AND.PRESENT(dw)) THEN
       DO ij = 1, limij
          DO kl = 1, limkl
             l = l + 1
             dw(l,1) = ww_dx(kl, ij)
             dw(l,2) = ww_dy(kl, ij)
             dw(l,3) = ww_dz(kl, ij)
          END DO
       END DO
    ELSE
       STOP "error"
    END IF

  END SUBROUTINE store_2el_2c_diag

END MODULE semi_empirical_int_utils
