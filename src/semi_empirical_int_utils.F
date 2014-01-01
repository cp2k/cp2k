!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for Integrals for semi-empiric methods
!> \author Teodoro Laino (03.2008) [tlaino] - University of Zurich
! *****************************************************************************
MODULE semi_empirical_int_utils

  USE f77_blas
  USE input_constants,                 ONLY: do_method_pchg,&
                                             do_se_IS_kdso_d
  USE kinds,                           ONLY: dp
  USE semi_empirical_int3_utils,       ONLY: charg_int_3,&
                                             dcharg_int_3,&
                                             ijkl_low_3
  USE semi_empirical_int_arrays,       ONLY: &
       CLMp, CLMxx, CLMxy, CLMyy, CLMz, CLMzp, CLMzz, clm_d, clm_sp, &
       ijkl_ind, indexa, indexb, int2c_type
  USE semi_empirical_types,            ONLY: rotmat_type,&
                                             se_int_control_type,&
                                             se_int_screen_type,&
                                             se_taper_type,&
                                             semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE
#include "semi_empirical_int_debug.h"

  PRIVATE
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_utils'

  PUBLIC ::   ijkl_sp,   ijkl_d, rotmat, rot_2el_2c_first, store_2el_2c_diag,&
            d_ijkl_sp, d_ijkl_d

CONTAINS

! *****************************************************************************
!> \brief General driver for computing semi-empirical integrals <ij|kl> with
!>        sp basis set. This code uses the old definitions of quadrupoles and
!>        therefore cannot be used for integrals involving d-orbitals (which
!>        require a definition of quadrupoles based on the rotational invariant
!>        property)
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  FUNCTION ijkl_sp(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_control,&
       se_int_screen, itype, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    INTEGER, INTENT(IN)                      :: itype
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl_sp', &
      routineP = moduleN//':'//routineN

    res = ijkl_sp_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                      se_int_control%integral_screening, se_int_control%shortrange,&
                      se_int_control%pc_coulomb_int, se_int_control%max_multipole,&
                      itype, charg_int_nri, error)

    ! If only the shortrange component is requested we can skip the rest
    IF ((.NOT.se_int_control%pc_coulomb_int).AND.(itype/=do_method_pchg)) THEN
       ! Handle the 1/r^3 term, this term is ALWAYS false for KDSO-D integrals
       IF (se_int_control%shortrange.AND.se_int_control%do_ewald_r3) THEN
          res = res - ijkl_low_3(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r,&
                                 itype, charg_int_3, error)
       END IF
    END IF
  END FUNCTION ijkl_sp

! *****************************************************************************
!> \brief General driver for computing derivatives of semi-empirical integrals
!>        <ij|kl> with sp basis set.
!>        This code uses the old definitions of quadrupoles and therefore
!>        cannot be used for integrals involving d-orbitals (which requires a
!>        definition of quadrupoles based on the rotational invariant property)
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 05.2008 [tlaino]
! *****************************************************************************
  FUNCTION d_ijkl_sp(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_control,&
       se_int_screen, itype, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    INTEGER, INTENT(IN)                      :: itype
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'd_ijkl_sp', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dfs, srd

    IF (se_int_control%integral_screening==do_se_IS_kdso_d) THEN
       res = ijkl_sp_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                         se_int_control%integral_screening, .FALSE.,&
                         se_int_control%pc_coulomb_int, se_int_control%max_multipole,&
                         itype, dcharg_int_nri, error)

       IF (.NOT.se_int_control%pc_coulomb_int) THEN
          dfs = ijkl_sp_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                            se_int_control%integral_screening, .FALSE., .FALSE., &
                            se_int_control%max_multipole, itype, dcharg_int_nri_fs, error)
          res = res + dfs*se_int_screen%dft

          ! In case we need the shortrange part we have to evaluate an additional derivative
          ! to handle the derivative of the Tapering term
          IF (se_int_control%shortrange) THEN
             srd = ijkl_sp_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                               se_int_control%integral_screening, .FALSE., .TRUE., &
                               se_int_control%max_multipole, itype, dcharg_int_nri, error)
             res = res - srd
          END IF
       END IF
    ELSE
       res = ijkl_sp_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                         se_int_control%integral_screening, se_int_control%shortrange,&
                         se_int_control%pc_coulomb_int, se_int_control%max_multipole,&
                         itype, dcharg_int_nri, error)
    END IF

    ! If only the shortrange component is requested we can skip the rest
    IF ((.NOT.se_int_control%pc_coulomb_int).AND.(itype/=do_method_pchg)) THEN
       ! Handle the 1/r^3 term, this term is ALWAYS false for KDSO-D integrals
       IF (se_int_control%shortrange.AND.se_int_control%do_ewald_r3) THEN
          res = res - ijkl_low_3(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r,&
                                 itype, dcharg_int_3, error)
       END IF
    END IF

  END FUNCTION d_ijkl_sp

! *****************************************************************************
!> \brief Low level general driver for computing semi-empirical integrals
!>        <ij|kl> and their derivatives with sp basis set only.
!>        This code uses the old definitions of quadrupoles and
!>        therefore cannot be used for integrals involving d-orbitals (which
!>        require a definition of quadrupoles based on the rotational invariant
!>        property)
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 05.2008 [tlaino]
! *****************************************************************************
  FUNCTION ijkl_sp_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
       iscreen, shortrange, pc_coulomb_int, max_multipole, itype, eval, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    INTEGER, INTENT(IN)                      :: iscreen
    LOGICAL, INTENT(IN)                      :: shortrange, pc_coulomb_int
    INTEGER, INTENT(IN)                      :: max_multipole, itype
    REAL(KIND=dp)                            :: eval
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl_sp_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ccc, l1, l1max, l1min, l2, &
                                                l2max, l2min, lij, lkl, lmin, &
                                                m
    REAL(KIND=dp)                            :: add, chrg, dij, dkl, fact_ij, &
                                                fact_kl, fact_screen, pij, &
                                                pkl, s1, sum

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
    dij = 0.0_dp
    dkl = 0.0_dp
    fact_ij = 1.0_dp
    fact_kl = 1.0_dp
    fact_screen = 1.0_dp
    IF (lij==3) fact_ij = SQRT(2.0_dp)
    IF (lkl==3) fact_kl = SQRT(2.0_dp)
    IF (.NOT.pc_coulomb_int) THEN
       IF (iscreen==do_se_IS_kdso_d) fact_screen = se_int_screen%ft
       ! Standard value of the integral
       DO l1 = l1min, l1max
          IF (l1 == 0) THEN
             IF (lij == 1) THEN
                pij = sepi%ko(1)
                IF (ic == -1 .OR. ic == 1) THEN
                   pij = sepi%ko(9)
                END IF
             ELSE IF (lij == 3) THEN
                pij = sepi%ko(7)
             END IF
          ELSE
             dij = sepi%cs(lij)*fact_ij
             pij = sepi%ko(lij)
          END IF
          !
          DO l2 = l2min, l2max
             IF (l2 == 0) THEN
                IF (lkl == 1) THEN
                   pkl = sepj%ko(1)
                   IF (ic == -1 .OR. ic == 2) THEN
                      pkl = sepj%ko(9)
                   END IF
                ELSE IF (lkl == 3) THEN
                   pkl = sepj%ko(7)
                END IF
             ELSE
                dkl = sepj%cs(lkl)*fact_kl
                pkl = sepj%ko(lkl)
             END IF
             IF (itype==do_method_pchg) THEN
                add = 0.0_dp
             ELSE
                add = (pij+pkl) ** 2
             END IF
             lmin = MAX(l1, l2)
             s1   = 0.0_dp
             DO m = -lmin, lmin
                ccc = clm_sp(ij, l1, m) * clm_sp(kl, l2, m)
                IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
                   chrg = eval(r, l1, l2, clm_sp(ij, l1, m), clm_sp(kl, l2, m), dij, dkl, add, fact_screen, error)
                   s1   = s1 + chrg
                END IF
             END DO
             sum = sum + s1
          END DO
       END DO
       res = sum
    END IF
    ! Shortrange: Possibly computes pure Coulomb and subtract from the original integral value
    IF (shortrange.OR.pc_coulomb_int) THEN
       sum = 0.0_dp
       dij = 0.0_dp
       dkl = 0.0_dp
       add = 0.0_dp
       fact_screen = 0.0_dp
       DO l1 = l1min, l1max
          IF (l1 > max_multipole) CYCLE
          IF (l1 /= 0) THEN
             dij = sepi%cs(lij)*fact_ij
          END IF
          !
          DO l2 = l2min, l2max
             IF (l2 > max_multipole) CYCLE
             IF (l2 /= 0) THEN
                dkl = sepj%cs(lkl)*fact_kl
             END IF
             lmin = MAX(l1, l2)
             s1   = 0.0_dp
             DO m = -lmin, lmin
                ccc = clm_sp(ij, l1, m) * clm_sp(kl, l2, m)
                IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
                   chrg = eval(r, l1, l2, clm_sp(ij, l1, m), clm_sp(kl, l2, m), dij, dkl, add, fact_screen, error)
                   s1   = s1 + chrg
                END IF
             END DO
             sum = sum + s1
          END DO
       END DO
       IF (pc_coulomb_int) res = sum
       IF (shortrange)     res = res - sum
    END IF
  END FUNCTION ijkl_sp_low

! *****************************************************************************
!> \brief Interaction function between two point-charge configurations NDDO sp-code
!>        Non-Rotational Invariant definition of quadrupoles
!>        r    -  Distance r12
!>        l1,m -  Quantum numbers for multipole of configuration 1
!>        l2,m -  Quantum numbers for multipole of configuration 2
!>        da   -  charge separation of configuration 1
!>        db   -  charge separation of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  FUNCTION charg_int_nri(r, l1_i, l2_i, m1_i, m2_i, da_i, db_i, add0, fact_screen, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1_i, l2_i, m1_i, m2_i
    REAL(KIND=dp), INTENT(in)                :: da_i, db_i, add0, fact_screen
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'charg_int_nri', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2, m1, m2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: add, da, db, dxdx, dxqxz, &
                                                dzdz, dzqxx, dzqzz, fact, &
                                                qqxx, qqzz, qxxqxx, qxxqyy, &
                                                qxzqxz, xyxy, zzzz

    failure = .FALSE.
    ! Computing only Integral Values
    IF      ( l1_i  < l2_i) THEN
       l1 = l1_i
       l2 = l2_i
       m1 = m1_i
       m2 = m2_i
       da = da_i
       db = db_i
       fact  = 1.0_dp
    ELSE IF ( l1_i  > l2_i ) THEN
       l1 = l2_i
       l2 = l1_i
       m1 = m2_i
       m2 = m1_i
       da = db_i
       db = da_i
       fact  = (-1.0_dp)**(l1+l2)
    ELSE IF ( l1_i == l2_i ) THEN
       l1 = l1_i
       l2 = l2_i
       IF   ( m1_i <= m2_i ) THEN
          m1 = m1_i
          m2 = m2_i
          da = da_i
          db = db_i
       ELSE
          m1 = m2_i
          m2 = m1_i
          da = db_i
          db = da_i
       END IF
       fact = 1.0_dp
    END IF
    add   = add0*fact_screen
    charg = 0.0_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg =  fact/SQRT(r**2+add)
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1 .AND. m2 == CLMz) THEN
       charg =  1.0_dp/SQRT((r+db)**2+add) - 1.0_dp/SQRT((r-db)**2+add)
       charg = charg*0.5_dp*fact
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m1 == CLMz .AND. m2 == CLMz ) THEN
       dzdz  = &
            +1.0_dp/SQRT((r+da-db)**2+add) + 1.0_dp/SQRT((r-da+db)**2+add) &
            -1.0_dp/SQRT((r-da-db)**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzdz*0.25_dp*fact
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m1 == CLMp .AND. m2 == CLMp ) THEN
       dxdx  = 2.0_dp/SQRT(r**2+(da-db)**2+add) - 2.0_dp/SQRT(r**2+(da+db)**2+add)
       charg = dxdx*0.25_dp*fact
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2 .AND. m2 == CLMzz) THEN
       qqzz  = 1.0_dp/SQRT((r-db)**2+add) - 2.0_dp/SQRT(r**2+add) + 1.0_dp/SQRT((r+db)**2+add)
       charg = qqzz*0.25_dp*fact
       RETURN
    END IF
    ! Q - XX
    IF (l1 == 0 .AND. l2 == 2 .AND. (m2 == CLMyy .OR. m2 == CLMxx)) THEN
       qqxx  = - 1.0_dp/SQRT(r**2+add) + 1.0_dp/SQRT(r**2+add+db**2)
       charg = qqxx*0.5_dp*fact
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMz .AND. m2 == CLMzz ) THEN
       dzqzz = &
            +1.0_dp/SQRT((r-da-db)**2+add)    - 2.0_dp/SQRT((r-da)**2+add) &
            +1.0_dp/SQRT((r-da+db)**2+add)    - 1.0_dp/SQRT((r+da-db)**2+add) &
            +2.0_dp/SQRT((r+da)**2+add)       - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzqzz*0.125_dp*fact
       RETURN
    END IF
    ! Z - XX
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMz .AND. ( m2 == CLMyy .OR. m2 == CLMxx )) THEN
       dzqxx = &
            +1.0_dp/SQRT((r+da)**2+add)    - 1.0_dp/SQRT((r+da)**2+add+db**2) &
            -1.0_dp/SQRT((r-da)**2+add)    + 1.0_dp/SQRT((r-da)**2+add+db**2)
       charg = dzqxx*0.25_dp*fact
       RETURN
    END IF
    ! ZZ - ZZ
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzz .AND. m2 == CLMzz ) THEN
       zzzz = &
            +1.0_dp/SQRT((r-da-db)**2+add)    + 1.0_dp/SQRT((r+da+db)**2+add) &
            +1.0_dp/SQRT((r-da+db)**2+add)    + 1.0_dp/SQRT((r+da-db)**2+add)
       xyxy = &
            +1.0_dp/SQRT((r-da)**2+add)       + 1.0_dp/SQRT((r+da)**2+add) &
            +1.0_dp/SQRT((r-db)**2+add)       + 1.0_dp/SQRT((r+db)**2+add) &
            -2.0_dp/SQRT(r**2+add)
       charg = (zzzz*0.0625_dp - xyxy*0.125_dp)*fact
       RETURN
    END IF
    ! ZZ - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzz .AND. ( m2 == CLMxx .OR. m2 == CLMyy )) THEN
       zzzz = &
            -1.0_dp/SQRT((r+da)**2+add)    + 1.0_dp/SQRT((r+da)**2+db**2+add) &
            -1.0_dp/SQRT((r-da)**2+add)    + 1.0_dp/SQRT((r-da)**2+db**2+add)
       xyxy = &
            +1.0_dp/SQRT(r**2+db**2+add)     - 1.0_dp/SQRT(r**2+add)
       charg = (zzzz*0.125_dp - xyxy*0.25_dp)*fact
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMp .AND. m2 == CLMzp ) THEN
       db = db/2.0_dp
       dxqxz = &
            -1.0_dp/SQRT((r-db)**2+(da-db)**2+add) + 1.0_dp/SQRT((r+db)**2+(da-db)**2+add) &
            +1.0_dp/SQRT((r-db)**2+(da+db)**2+add) - 1.0_dp/SQRT((r+db)**2+(da+db)**2+add)
       charg = dxqxz*0.25_dp*fact
       RETURN
    END IF
    ! ZX - ZX
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzp .AND. m2 == CLMzp ) THEN
       da = da / 2.0_dp
       db = db / 2.0_dp
       qxzqxz = &
            +1.0_dp/SQRT((r+da-db)**2+(da-db)**2+add) - 1.0_dp/SQRT((r+da+db)**2+(da-db)**2+add) &
            -1.0_dp/SQRT((r-da-db)**2+(da-db)**2+add) + 1.0_dp/SQRT((r-da+db)**2+(da-db)**2+add) &
            -1.0_dp/SQRT((r+da-db)**2+(da+db)**2+add) + 1.0_dp/SQRT((r+da+db)**2+(da+db)**2+add) &
            +1.0_dp/SQRT((r-da-db)**2+(da+db)**2+add) - 1.0_dp/SQRT((r-da+db)**2+(da+db)**2+add)
       charg = qxzqxz*0.125_dp*fact
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. (((m1 == CLMyy) .AND. (m2 == CLMyy)).OR.((m1 == CLMxx) .AND. (m2 == CLMxx)))) THEN
       qxxqxx = &
            + 2.0_dp/SQRT(r**2+add)                   + 1.0_dp/SQRT(r**2+(da-db)**2+add) &
            + 1.0_dp/SQRT(r**2+(da+db)**2+add)        - 2.0_dp/SQRT(r**2+da**2+add)      &
            - 2.0_dp/SQRT(r**2+db**2+add)
       charg = qxxqxx*0.125_dp*fact
       RETURN
    END IF
    ! XX - YY
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMyy .AND. m2 == CLMxx) THEN
       qxxqyy = &
            + 1.0_dp/SQRT(r**2+add)                   - 1.0_dp/SQRT(r**2+da**2+add) &
            - 1.0_dp/SQRT(r**2+db**2+add)             + 1.0_dp/SQRT(r**2+da**2+db**2+add)
       charg = qxxqyy*0.25_dp*fact
       RETURN
    END IF
    ! XY - XY
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMxy .AND. m2 == CLMxy) THEN
       qxxqxx = &
            + 2.0_dp/SQRT(r**2+add)                   + 1.0_dp/SQRT(r**2+(da-db)**2+add) &
            + 1.0_dp/SQRT(r**2+(da+db)**2+add)        - 2.0_dp/SQRT(r**2+da**2+add)      &
            - 2.0_dp/SQRT(r**2+db**2+add)
       qxxqyy = &
            + 1.0_dp/SQRT(r**2+add)                   - 1.0_dp/SQRT(r**2+da**2+add) &
            - 1.0_dp/SQRT(r**2+db**2+add)             + 1.0_dp/SQRT(r**2+da**2+db**2+add)
       charg = 0.5_dp * (qxxqxx*0.125_dp - qxxqyy*0.25_dp) * fact
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION charg_int_nri

! *****************************************************************************
!> \brief Derivatives of interaction function between two point-charge
!>        configurations NDDO sp-code.
!>        Non-Rotational Invariant definition of quadrupoles
!>
!>        r    -  Distance r12
!>        l1,m -  Quantum numbers for multipole of configuration 1
!>        l2,m -  Quantum numbers for multipole of configuration 2
!>        da   -  charge separation of configuration 1
!>        db   -  charge separation of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  FUNCTION dcharg_int_nri(r, l1_i, l2_i, m1_i, m2_i, da_i, db_i, add0, fact_screen, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1_i, l2_i, m1_i, m2_i
    REAL(KIND=dp), INTENT(in)                :: da_i, db_i, add0, fact_screen
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'dcharg_int_nri', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2, m1, m2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: add, da, db, dxdx, dxqxz, &
                                                dzdz, dzqxx, dzqzz, fact, &
                                                qqxx, qqzz, qxxqxx, qxxqyy, &
                                                qxzqxz, xyxy, zzzz

    failure = .FALSE.
    ! Computing only Integral Derivatives
    IF      ( l1_i  < l2_i) THEN
       l1 = l1_i
       l2 = l2_i
       m1 = m1_i
       m2 = m2_i
       da = da_i
       db = db_i
       fact  = 1.0_dp
    ELSE IF ( l1_i  > l2_i ) THEN
       l1 = l2_i
       l2 = l1_i
       m1 = m2_i
       m2 = m1_i
       da = db_i
       db = da_i
       fact  = (-1.0_dp)**(l1+l2)
    ELSE IF ( l1_i == l2_i ) THEN
       l1 = l1_i
       l2 = l2_i
       IF   ( m1_i <= m2_i ) THEN
          m1 = m1_i
          m2 = m2_i
          da = da_i
          db = db_i
       ELSE
          m1 = m2_i
          m2 = m1_i
          da = db_i
          db = da_i
       END IF
       fact = 1.0_dp
    END IF
    charg = 0.0_dp
    add   = add0 * fact_screen
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg = r/SQRT(r**2+add)**3
       charg = -charg*fact
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1 .AND. m2 == CLMz) THEN
       charg = (r+db)/SQRT((r+db)**2+add)**3 - (r-db)/SQRT((r-db)**2+add)**3
       charg = -charg*0.5_dp*fact
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m1 == CLMz .AND. m2 == CLMz ) THEN
       dzdz  = &
            +(r+da-db)/SQRT((r+da-db)**2+add)**3 + (r-da+db)/SQRT((r-da+db)**2+add)**3 &
            -(r-da-db)/SQRT((r-da-db)**2+add)**3 - (r+da+db)/SQRT((r+da+db)**2+add)**3
       charg = -dzdz*0.25_dp*fact
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m1 == CLMp .AND. m2 == CLMp ) THEN
       dxdx  = 2.0_dp*r/SQRT(r**2+(da-db)**2+add)**3 - 2.0_dp*r/SQRT(r**2+(da+db)**2+add)**3
       charg = -dxdx*0.25_dp*fact
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2 .AND. m2 == CLMzz) THEN
       qqzz  = (r-db)/SQRT((r-db)**2+add)**3 - 2.0_dp*r/SQRT(r**2+add)**3 + (r+db)/SQRT((r+db)**2+add)**3
       charg = -qqzz*0.25_dp*fact
       RETURN
    END IF
    ! Q - XX
    IF (l1 == 0 .AND. l2 == 2 .AND. (m2 == CLMyy .OR. m2 == CLMxx)) THEN
       qqxx  = -r/SQRT(r**2+add)**3 + r/SQRT(r**2+add+db**2)**3
       charg = -qqxx*0.5_dp*fact
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMz .AND. m2 == CLMzz ) THEN
       dzqzz = &
            +(r-da-db)/SQRT((r-da-db)**2+add)**3    - 2.0_dp*(r-da)/SQRT((r-da)**2+add)**3 &
            +(r-da+db)/SQRT((r-da+db)**2+add)**3    - (r+da-db)/SQRT((r+da-db)**2+add)**3 &
            +2.0_dp*(r+da)/SQRT((r+da)**2+add)**3   - (r+da+db)/SQRT((r+da+db)**2+add)**3
       charg = -dzqzz*0.125_dp*fact
       RETURN
    END IF
    ! Z - XX
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMz .AND. ( m2 == CLMyy .OR. m2 == CLMxx )) THEN
       dzqxx = &
            +(r+da)/SQRT((r+da)**2+add)**3    - (r+da)/SQRT((r+da)**2+add+db**2)**3 &
            -(r-da)/SQRT((r-da)**2+add)**3    + (r-da)/SQRT((r-da)**2+add+db**2)**3
       charg = -dzqxx*0.25_dp*fact
       RETURN
    END IF
    ! ZZ - ZZ
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzz .AND. m2 == CLMzz ) THEN
       zzzz = &
            +(r-da-db)/SQRT((r-da-db)**2+add)**3    + (r+da+db)/SQRT((r+da+db)**2+add)**3 &
            +(r-da+db)/SQRT((r-da+db)**2+add)**3    + (r+da-db)/SQRT((r+da-db)**2+add)**3
       xyxy = &
            +(r-da)/SQRT((r-da)**2+add)**3       + (r+da)/SQRT((r+da)**2+add)**3 &
            +(r-db)/SQRT((r-db)**2+add)**3       + (r+db)/SQRT((r+db)**2+add)**3 &
            -2.0_dp*r/SQRT(r**2+add)**3
       charg = -(zzzz*0.0625_dp - xyxy*0.125_dp)*fact
       RETURN
    END IF
    ! ZZ - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzz .AND. ( m2 == CLMxx .OR. m2 == CLMyy )) THEN
       zzzz = &
            -(r+da)/SQRT((r+da)**2+add)**3    + (r+da)/SQRT((r+da)**2+db**2+add)**3 &
            -(r-da)/SQRT((r-da)**2+add)**3    + (r-da)/SQRT((r-da)**2+db**2+add)**3
       xyxy =   r/SQRT(r**2+db**2+add)**3     - r/SQRT(r**2+add)**3
       charg = -(zzzz*0.125_dp - xyxy*0.25_dp)*fact
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMp .AND. m2 == CLMzp ) THEN
       db = db/2.0_dp
       dxqxz = &
            -(r-db)/SQRT((r-db)**2+(da-db)**2+add)**3 + (r+db)/SQRT((r+db)**2+(da-db)**2+add)**3 &
            +(r-db)/SQRT((r-db)**2+(da+db)**2+add)**3 - (r+db)/SQRT((r+db)**2+(da+db)**2+add)**3
       charg = -dxqxz*0.25_dp*fact
       RETURN
    END IF
    ! ZX - ZX
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzp .AND. m2 == CLMzp ) THEN
       da = da / 2.0_dp
       db = db / 2.0_dp
       qxzqxz = &
            +(r+da-db)/SQRT((r+da-db)**2+(da-db)**2+add)**3 - (r+da+db)/SQRT((r+da+db)**2+(da-db)**2+add)**3 &
            -(r-da-db)/SQRT((r-da-db)**2+(da-db)**2+add)**3 + (r-da+db)/SQRT((r-da+db)**2+(da-db)**2+add)**3 &
            -(r+da-db)/SQRT((r+da-db)**2+(da+db)**2+add)**3 + (r+da+db)/SQRT((r+da+db)**2+(da+db)**2+add)**3 &
            +(r-da-db)/SQRT((r-da-db)**2+(da+db)**2+add)**3 - (r-da+db)/SQRT((r-da+db)**2+(da+db)**2+add)**3
       charg = -qxzqxz*0.125_dp*fact
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. (((m1 == CLMyy) .AND. (m2 == CLMyy)).OR.((m1 == CLMxx) .AND. (m2 == CLMxx)))) THEN
       qxxqxx = &
            + 2.0_dp*r/SQRT(r**2+add)**3                   +        r/SQRT(r**2+(da-db)**2+add)**3 &
            +        r/SQRT(r**2+(da+db)**2+add)**3        - 2.0_dp*r/SQRT(r**2+da**2+add)**3      &
            - 2.0_dp*r/SQRT(r**2+db**2+add)**3
       charg = -qxxqxx*0.125_dp*fact
       RETURN
    END IF
    ! XX - YY
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMyy .AND. m2 == CLMxx) THEN
       qxxqyy = &
            + r/SQRT(r**2+add)**3                   - r/SQRT(r**2+da**2+add)**3 &
            - r/SQRT(r**2+db**2+add)**3             + r/SQRT(r**2+da**2+db**2+add)**3
       charg = -qxxqyy*0.25_dp*fact
       RETURN
    END IF
    ! XY - XY
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMxy .AND. m2 == CLMxy) THEN
       qxxqxx = &
            + 2.0_dp*r/SQRT(r**2+add)**3                    +        r/SQRT(r**2+(da-db)**2+add)**3  &
            +        r/SQRT(r**2+(da+db)**2+add)**3         - 2.0_dp*r/SQRT(r**2+da**2+add)**3       &
            - 2.0_dp*r/SQRT(r**2+db**2+add)**3
       qxxqyy = &
            + r/SQRT(r**2+add)**3                   - r/SQRT(r**2+da**2+add)**3 &
            - r/SQRT(r**2+db**2+add)**3             + r/SQRT(r**2+da**2+db**2+add)**3
       charg = -0.5_dp * (qxxqxx*0.125_dp - qxxqyy*0.25_dp) * fact
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION dcharg_int_nri

! *****************************************************************************
!> \brief Derivatives of interaction function between two point-charge
!>        configurations NDDO sp-code. The derivative takes care of the screening
!>        term only.
!>        Non-Rotational Invariant definition of quadrupoles
!>
!>        r    -  Distance r12
!>        l1,m -  Quantum numbers for multipole of configuration 1
!>        l2,m -  Quantum numbers for multipole of configuration 2
!>        da   -  charge separation of configuration 1
!>        db   -  charge separation of configuration 2
!>        add  -  additive term
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  FUNCTION dcharg_int_nri_fs(r, l1_i, l2_i, m1_i, m2_i, da_i, db_i, add0, fact_screen, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1_i, l2_i, m1_i, m2_i
    REAL(KIND=dp), INTENT(in)                :: da_i, db_i, add0, fact_screen
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'dcharg_int_nri_fs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2, m1, m2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: add, da, db, dxdx, dxqxz, &
                                                dzdz, dzqxx, dzqzz, fact, &
                                                qqxx, qqzz, qxxqxx, qxxqyy, &
                                                qxzqxz, xyxy, zzzz

    failure = .FALSE.
    ! Computing only Integral Derivatives
    IF      ( l1_i  < l2_i) THEN
       l1 = l1_i
       l2 = l2_i
       m1 = m1_i
       m2 = m2_i
       da = da_i
       db = db_i
       fact  = 1.0_dp
    ELSE IF ( l1_i  > l2_i ) THEN
       l1 = l2_i
       l2 = l1_i
       m1 = m2_i
       m2 = m1_i
       da = db_i
       db = da_i
       fact  = (-1.0_dp)**(l1+l2)
    ELSE IF ( l1_i == l2_i ) THEN
       l1 = l1_i
       l2 = l2_i
       IF   ( m1_i <= m2_i ) THEN
          m1 = m1_i
          m2 = m2_i
          da = da_i
          db = db_i
       ELSE
          m1 = m2_i
          m2 = m1_i
          da = db_i
          db = da_i
       END IF
       fact = 1.0_dp
    END IF
    charg = 0.0_dp
    add   = add0 * fact_screen
    ! The 0.5 factor handles the derivative of the SQRT
    fact  = fact * 0.5_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg = add0/SQRT(r**2+add)**3
       charg = -charg*fact
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1 .AND. m2 == CLMz) THEN
       charg = add0/SQRT((r+db)**2+add)**3 - add0/SQRT((r-db)**2+add)**3
       charg = -charg*0.5_dp*fact
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m1 == CLMz .AND. m2 == CLMz ) THEN
       dzdz  = &
            +add0/SQRT((r+da-db)**2+add)**3 + add0/SQRT((r-da+db)**2+add)**3 &
            -add0/SQRT((r-da-db)**2+add)**3 - add0/SQRT((r+da+db)**2+add)**3
       charg = -dzdz*0.25_dp*fact
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m1 == CLMp .AND. m2 == CLMp ) THEN
       dxdx  = 2.0_dp*add0/SQRT(r**2+(da-db)**2+add)**3 - 2.0_dp*add0/SQRT(r**2+(da+db)**2+add)**3
       charg = -dxdx*0.25_dp*fact
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2 .AND. m2 == CLMzz) THEN
       qqzz  = add0/SQRT((r-db)**2+add)**3 - 2.0_dp*add0/SQRT(r**2+add)**3 + add0/SQRT((r+db)**2+add)**3
       charg = -qqzz*0.25_dp*fact
       RETURN
    END IF
    ! Q - XX
    IF (l1 == 0 .AND. l2 == 2 .AND. (m2 == CLMyy .OR. m2 == CLMxx)) THEN
       qqxx  = -add0/SQRT(r**2+add)**3 + add0/SQRT(r**2+add+db**2)**3
       charg = -qqxx*0.5_dp*fact
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMz .AND. m2 == CLMzz ) THEN
       dzqzz = &
            +add0/SQRT((r-da-db)**2+add)**3      - 2.0_dp*add0/SQRT((r-da)**2+add)**3 &
            +add0/SQRT((r-da+db)**2+add)**3      - add0/SQRT((r+da-db)**2+add)**3 &
            +2.0_dp*add0/SQRT((r+da)**2+add)**3  - add0/SQRT((r+da+db)**2+add)**3
       charg = -dzqzz*0.125_dp*fact
       RETURN
    END IF
    ! Z - XX
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMz .AND. ( m2 == CLMyy .OR. m2 == CLMxx )) THEN
       dzqxx = &
            +add0/SQRT((r+da)**2+add)**3    - add0/SQRT((r+da)**2+add+db**2)**3 &
            -add0/SQRT((r-da)**2+add)**3    + add0/SQRT((r-da)**2+add+db**2)**3
       charg = -dzqxx*0.25_dp*fact
       RETURN
    END IF
    ! ZZ - ZZ
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzz .AND. m2 == CLMzz ) THEN
       zzzz = &
            +add0/SQRT((r-da-db)**2+add)**3    + add0/SQRT((r+da+db)**2+add)**3 &
            +add0/SQRT((r-da+db)**2+add)**3    + add0/SQRT((r+da-db)**2+add)**3
       xyxy = &
            +add0/SQRT((r-da)**2+add)**3       + add0/SQRT((r+da)**2+add)**3 &
            +add0/SQRT((r-db)**2+add)**3       + add0/SQRT((r+db)**2+add)**3 &
            -2.0_dp*add0/SQRT(r**2+add)**3
       charg = -(zzzz*0.0625_dp - xyxy*0.125_dp)*fact
       RETURN
    END IF
    ! ZZ - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzz .AND. ( m2 == CLMxx .OR. m2 == CLMyy )) THEN
       zzzz = &
            -add0/SQRT((r+da)**2+add)**3    + add0/SQRT((r+da)**2+db**2+add)**3 &
            -add0/SQRT((r-da)**2+add)**3    + add0/SQRT((r-da)**2+db**2+add)**3
       xyxy =add0/SQRT(r**2+db**2+add)**3   - add0/SQRT(r**2+add)**3
       charg = -(zzzz*0.125_dp - xyxy*0.25_dp)*fact
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m1 == CLMp .AND. m2 == CLMzp ) THEN
       db = db/2.0_dp
       dxqxz = &
            -add0/SQRT((r-db)**2+(da-db)**2+add)**3 + add0/SQRT((r+db)**2+(da-db)**2+add)**3 &
            +add0/SQRT((r-db)**2+(da+db)**2+add)**3 - add0/SQRT((r+db)**2+(da+db)**2+add)**3
       charg = -dxqxz*0.25_dp*fact
       RETURN
    END IF
    ! ZX - ZX
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMzp .AND. m2 == CLMzp ) THEN
       da = da / 2.0_dp
       db = db / 2.0_dp
       qxzqxz = &
            +add0/SQRT((r+da-db)**2+(da-db)**2+add)**3 - add0/SQRT((r+da+db)**2+(da-db)**2+add)**3 &
            -add0/SQRT((r-da-db)**2+(da-db)**2+add)**3 + add0/SQRT((r-da+db)**2+(da-db)**2+add)**3 &
            -add0/SQRT((r+da-db)**2+(da+db)**2+add)**3 + add0/SQRT((r+da+db)**2+(da+db)**2+add)**3 &
            +add0/SQRT((r-da-db)**2+(da+db)**2+add)**3 - add0/SQRT((r-da+db)**2+(da+db)**2+add)**3
       charg = -qxzqxz*0.125_dp*fact
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. (((m1 == CLMyy) .AND. (m2 == CLMyy)).OR.((m1 == CLMxx) .AND. (m2 == CLMxx)))) THEN
       qxxqxx = &
            + 2.0_dp*add0/SQRT(r**2+add)**3                   +        add0/SQRT(r**2+(da-db)**2+add)**3 &
            +        add0/SQRT(r**2+(da+db)**2+add)**3        - 2.0_dp*add0/SQRT(r**2+da**2+add)**3      &
            - 2.0_dp*add0/SQRT(r**2+db**2+add)**3
       charg = -qxxqxx*0.125_dp*fact
       RETURN
    END IF
    ! XX - YY
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMyy .AND. m2 == CLMxx) THEN
       qxxqyy = &
            + add0/SQRT(r**2+add)**3                   - add0/SQRT(r**2+da**2+add)**3 &
            - add0/SQRT(r**2+db**2+add)**3             + add0/SQRT(r**2+da**2+db**2+add)**3
       charg = -qxxqyy*0.25_dp*fact
       RETURN
    END IF
    ! XY - XY
    IF (l1 == 2 .AND. l2 == 2 .AND. m1 == CLMxy .AND. m2 == CLMxy) THEN
       qxxqxx = &
            + 2.0_dp*add0/SQRT(r**2+add)**3                    +        add0/SQRT(r**2+(da-db)**2+add)**3  &
            +        add0/SQRT(r**2+(da+db)**2+add)**3         - 2.0_dp*add0/SQRT(r**2+da**2+add)**3       &
            - 2.0_dp*add0/SQRT(r**2+db**2+add)**3
       qxxqyy = &
            + add0/SQRT(r**2+add)**3                   - add0/SQRT(r**2+da**2+add)**3 &
            - add0/SQRT(r**2+db**2+add)**3             + add0/SQRT(r**2+da**2+db**2+add)**3
       charg = -0.5_dp * (qxxqxx*0.125_dp - qxxqyy*0.25_dp) * fact
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION dcharg_int_nri_fs

! *****************************************************************************
!> \brief General driver for computing semi-empirical integrals <ij|kl>
!>        involving d-orbitals.
!>        The choice of the linear quadrupole was REALLY unhappy
!>        in the first development of the NDDO codes. That choice makes
!>        impossible the merging of the integral code with or without d-orbitals
!>        unless a reparametrization of all NDDO codes for s and p orbitals will
!>        be performed.. more over the choice of the linear quadrupole does not make
!>        calculations rotational invariants (of course the rotational invariant
!>        can be forced). The definitions of quadrupoles for d-orbitals is the
!>        correct one in order to have the rotational invariant property by
!>        construction..
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 03.2008 [tlaino]
! *****************************************************************************
  FUNCTION ijkl_d(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_control, &
       se_int_screen, itype, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    INTEGER, INTENT(IN)                      :: itype
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl_d', &
      routineP = moduleN//':'//routineN

    res = ijkl_d_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                     se_int_control%integral_screening, se_int_control%shortrange,&
                     se_int_control%pc_coulomb_int, se_int_control%max_multipole,&
                     itype, charg_int_ri, error)

    ! If only the shortrange component is requested we can skip the rest
    IF ((.NOT.se_int_control%pc_coulomb_int).AND.(itype/=do_method_pchg)) THEN
       ! Handle the 1/r^3 term, this term is ALWAYS false for KDSO-D integrals
       IF (se_int_control%shortrange.AND.se_int_control%do_ewald_r3) THEN
          res = res - ijkl_low_3(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r,&
                                 itype, charg_int_3, error)
       END IF
    END IF
  END FUNCTION ijkl_d

! *****************************************************************************
!> \brief General driver for computing the derivatives of semi-empirical integrals <ij|kl>
!>        involving d-orbitals.
!>        The choice of the linear quadrupole was REALLY unhappy
!>        in the first development of the NDDO codes. That choice makes
!>        impossible the merging of the integral code with or without d-orbitals
!>        unless a reparametrization of all NDDO codes for s and p orbitals will
!>        be performed.. more over the choice of the linear quadrupole does not make
!>        calculations rotational invariants (of course the rotational invariant
!>        can be forced). The definitions of quadrupoles for d-orbitals is the
!>        correct one in order to have the rotational invariant property by
!>        construction..
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 03.2008 [tlaino]
! *****************************************************************************
  FUNCTION d_ijkl_d(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_control,&
       se_int_screen, itype, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    INTEGER, INTENT(IN)                      :: itype
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'd_ijkl_d', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: dfs, srd

    IF (se_int_control%integral_screening==do_se_IS_kdso_d) THEN
       res = ijkl_d_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                        se_int_control%integral_screening, .FALSE.,&
                        se_int_control%pc_coulomb_int, se_int_control%max_multipole,&
                        itype, dcharg_int_ri, error)

       IF (.NOT.se_int_control%pc_coulomb_int) THEN
          dfs = ijkl_d_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                        se_int_control%integral_screening, .FALSE., .FALSE.,&
                        se_int_control%max_multipole, itype, dcharg_int_ri_fs, error)
          res = res + dfs*se_int_screen%dft

          ! In case we need the shortrange part we have to evaluate an additional derivative
          ! to handle the derivative of the Tapering term
          IF (se_int_control%shortrange) THEN
             srd = ijkl_d_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                              se_int_control%integral_screening, .FALSE., .TRUE.,&
                              se_int_control%max_multipole, itype, dcharg_int_ri, error)
             res = res - srd
          END IF
       END IF
    ELSE
       res = ijkl_d_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
                        se_int_control%integral_screening, se_int_control%shortrange,&
                        se_int_control%pc_coulomb_int, se_int_control%max_multipole,&
                        itype, dcharg_int_ri, error)
    END IF

    ! If only the shortrange component is requested we can skip the rest
    IF ((.NOT.se_int_control%pc_coulomb_int).AND.(itype/=do_method_pchg)) THEN
       ! Handle the 1/r^3 term, this term is ALWAYS false for KDSO-D integrals
       IF (se_int_control%shortrange.AND.se_int_control%do_ewald_r3) THEN
          res = res - ijkl_low_3(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r,&
                                 itype, dcharg_int_3, error)
       END IF
    END IF

  END FUNCTION d_ijkl_d

! *****************************************************************************
!> \brief Low level general driver for computing semi-empirical integrals <ij|kl>
!>        and their derivatives involving d-orbitals.
!>        The choice of the linear quadrupole was REALLY unhappy
!>        in the first development of the NDDO codes. That choice makes
!>        impossible the merging of the integral code with or without d-orbitals
!>        unless a reparametrization of all NDDO codes for s and p orbitals will
!>        be performed.. more over the choice of the linear quadrupole does not make
!>        calculations rotational invariants (of course the rotational invariant
!>        can be forced). The definitions of quadrupoles for d-orbitals is the
!>        correct one in order to have the rotational invariant property by
!>        construction..
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 03.2008 [tlaino]
! *****************************************************************************
  FUNCTION ijkl_d_low(sepi, sepj, ij, kl, li, lj, lk, ll, ic, r, se_int_screen,&
       iscreen, shortrange, pc_coulomb_int, max_multipole, itype, eval, error) RESULT(res)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ij, kl, li, lj, lk, ll, ic
    REAL(KIND=dp), INTENT(IN)                :: r
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    INTEGER, INTENT(IN)                      :: iscreen
    LOGICAL, INTENT(IN)                      :: shortrange, pc_coulomb_int
    INTEGER, INTENT(IN)                      :: max_multipole, itype
    REAL(KIND=dp)                            :: eval
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'ijkl_d_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l1max, l1min, l2, l2max, &
                                                l2min, lij, lkl, lmin, m, mm
    REAL(KIND=dp)                            :: add, ccc, chrg, dij, dkl, &
                                                fact_screen, pij, pkl, s1, sum

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
    dij = 0.0_dp
    dkl = 0.0_dp
    fact_screen = 1.0_dp
    IF (.NOT.pc_coulomb_int) THEN
       IF (iscreen==do_se_IS_kdso_d) fact_screen = se_int_screen%ft
       ! Standard value of the integral
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
             IF (itype==do_method_pchg) THEN
                add = 0.0_dp
             ELSE
                add = (pij+pkl) ** 2
             END IF
             lmin = MIN(l1, l2)
             s1 = 0.0_dp
             DO m = -lmin, lmin
                ccc = clm_d(ij, l1, m) * clm_d(kl, l2, m)
                IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
                   mm   = ABS (m)
                   chrg = eval(r, l1, l2, mm, dij, dkl, add, fact_screen, error)
                   s1   = s1 + chrg * ccc
                END IF
             END DO
             sum = sum + s1
          END DO
       END DO
       res = sum
    END IF
    ! Shortrange: Possibly computes pure Coulomb and subtract from the original integral valeu
    IF (shortrange.OR.pc_coulomb_int) THEN
       sum = 0.0_dp
       dij = 0.0_dp
       dkl = 0.0_dp
       add = 0.0_dp
       fact_screen = 0.0_dp
       DO l1 = l1min, l1max
          IF (l1 > max_multipole) CYCLE
          IF (l1 /= 0) THEN
             dij = sepi%cs(lij)
          END IF
          !
          DO l2 = l2min, l2max
             IF (l2 > max_multipole) CYCLE
             IF (l2 /= 0) THEN
                dkl = sepj%cs(lkl)
             END IF
             lmin = MIN(l1, l2)
             s1   = 0.0_dp
             DO m = -lmin, lmin
                ccc = clm_d(ij, l1, m) * clm_d(kl, l2, m)
                IF (ABS(ccc) > EPSILON(0.0_dp)) THEN
                   mm   = ABS (m)
                   chrg = eval(r, l1, l2, mm, dij, dkl, add, fact_screen, error)
                   s1   = s1 + chrg * ccc
                END IF
             END DO
             sum = sum + s1
          END DO
       END DO
       IF (pc_coulomb_int) res = sum
       IF (shortrange)     res = res - sum
    END IF
  END FUNCTION ijkl_d_low

! *****************************************************************************
!> \brief Interaction function between two point-charge configurations (MNDO/d)
!>        Rotational invariant property built-in in the quadrupole definition
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
  FUNCTION charg_int_ri(r, l1_i, l2_i, m, da_i, db_i, add0, fact_screen, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1_i, l2_i, m
    REAL(KIND=dp), INTENT(in)                :: da_i, db_i, add0, fact_screen
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'charg_int_ri', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, ab, add, da, db, dxdx, &
                                                dxqxz, dzdz, dzqzz, fact, &
                                                qqzz, qxzqxz, xyxy, zzzz

    failure = .FALSE.
    IF      ( l1_i  <= l2_i) THEN
       l1 = l1_i
       l2 = l2_i
       da = da_i
       db = db_i
       fact  = 1.0_dp
    ELSE IF ( l1_i  > l2_i ) THEN
       l1 = l2_i
       l2 = l1_i
       da = db_i
       db = da_i
       fact  = (-1.0_dp)**(l1+l2)
    END IF
    charg = 0.0_dp
    add   = add0*fact_screen
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg = fact/SQRT(r**2+add)
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1) THEN
       charg =  1.0_dp/SQRT((r+db)**2+add) - 1.0_dp/SQRT((r-db)**2+add)
       charg = charg*0.5_dp*fact
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 0) THEN
       dzdz  = &
            +1.0_dp/SQRT((r+da-db)**2+add) + 1.0_dp/SQRT((r-da+db)**2+add) &
            -1.0_dp/SQRT((r-da-db)**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzdz*0.25_dp*fact
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 1) THEN
       dxdx  = 2.0_dp/SQRT(r**2+(da-db)**2+add) - 2.0_dp/SQRT(r**2+(da+db)**2+add)
       charg = dxdx*0.25_dp*fact
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2) THEN
       qqzz  = 1.0_dp/SQRT((r-db)**2+add) - 2.0_dp/SQRT(r**2+db**2+add) + 1.0_dp/SQRT((r+db)**2+add)
       charg = qqzz*0.25_dp*fact
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 0) THEN
       dzqzz = &
            +1.0_dp/SQRT((r-da-db)**2+add)    - 2.0_dp/SQRT((r-da)**2+db**2+add) &
            +1.0_dp/SQRT((r+db-da)**2+add)    - 1.0_dp/SQRT((r-db+da)**2+add)    &
            +2.0_dp/SQRT((r+da)**2+db**2+add) - 1.0_dp/SQRT((r+da+db)**2+add)
       charg = dzqzz*0.125_dp*fact
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
       charg = (zzzz*0.0625_dp - xyxy*0.015625_dp)*fact
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 1) THEN
       ab = db/SQRT(2.0_dp)
       dxqxz = &
            -2.0_dp/SQRT((r-ab)**2+(da-ab)**2+add) + 2.0_dp/SQRT((r+ab)**2+(da-ab)**2+add) &
            +2.0_dp/SQRT((r-ab)**2+(da+ab)**2+add) - 2.0_dp/SQRT((r+ab)**2+(da+ab)**2+add)
       charg = dxqxz*0.125_dp*fact
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
       charg = qxzqxz*0.0625_dp*fact
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 2) THEN
       xyxy = 4.0_dp/SQRT(r**2+(da-db)**2+add) + 4.0_dp/SQRT(r**2+(da+db)**2+add) - 8.0_dp/SQRT(r**2+da**2+db**2+add)
       charg = xyxy*0.0625_dp*fact
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION charg_int_ri

! *****************************************************************************
!> \brief Derivatives of the interaction function between two point-charge
!>        configurations (MNDO/d)
!>        Rotational invariant property built-in in the quadrupole definition
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
  FUNCTION dcharg_int_ri(r, l1_i, l2_i, m, da_i, db_i, add0, fact_screen, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1_i, l2_i, m
    REAL(KIND=dp), INTENT(in)                :: da_i, db_i, add0, fact_screen
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'dcharg_int_ri', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, ab, add, da, db, dxdx, &
                                                dxqxz, dzdz, dzqzz, fact, &
                                                qqzz, qxzqxz, xyxy, zzzz

    failure = .FALSE.
    IF      ( l1_i  <= l2_i) THEN
       l1 = l1_i
       l2 = l2_i
       da = da_i
       db = db_i
       fact  = 1.0_dp
    ELSE IF ( l1_i  > l2_i ) THEN
       l1 = l2_i
       l2 = l1_i
       da = db_i
       db = da_i
       fact  = (-1.0_dp)**(l1+l2)
    END IF
    charg = 0.0_dp
    add   = add0*fact_screen
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg = r/SQRT(r**2+add)**3
       charg = -fact*charg
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1) THEN
       charg =  (r+db)/SQRT((r+db)**2+add)**3 - (r-db)/SQRT((r-db)**2+add)**3
       charg = -charg*0.5_dp*fact
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 0) THEN
       dzdz  = &
            +(r+da-db)/SQRT((r+da-db)**2+add)**3 + (r-da+db)/SQRT((r-da+db)**2+add)**3 &
            -(r-da-db)/SQRT((r-da-db)**2+add)**3 - (r+da+db)/SQRT((r+da+db)**2+add)**3
       charg = -dzdz*0.25_dp*fact
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 1) THEN
       dxdx  = 2.0_dp*r/SQRT(r**2+(da-db)**2+add)**3 - 2.0_dp*r/SQRT(r**2+(da+db)**2+add)**3
       charg = -dxdx*0.25_dp*fact
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2) THEN
       qqzz  = (r-db)/SQRT((r-db)**2+add)**3 - 2.0_dp*r/SQRT(r**2+db**2+add)**3 + (r+db)/SQRT((r+db)**2+add)**3
       charg = -qqzz*0.25_dp*fact
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 0) THEN
       dzqzz = &
            +(r-da-db)/SQRT((r-da-db)**2+add)**3        - 2.0_dp*(r-da)/SQRT((r-da)**2+db**2+add)**3 &
            +(r+db-da)/SQRT((r+db-da)**2+add)**3        -     (r-db+da)/SQRT((r-db+da)**2+add)**3    &
            +2.0_dp*(r+da)/SQRT((r+da)**2+db**2+add)**3 -     (r+da+db)/SQRT((r+da+db)**2+add)**3
       charg = -dzqzz*0.125_dp*fact
       RETURN
    END IF
    ! ZZ - ZZ
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 0) THEN
       zzzz = &
            +(r-da-db)/SQRT((r-da-db)**2+add)**3        + (r+da+db)/SQRT((r+da+db)**2+add)**3    &
            +(r-da+db)/SQRT((r-da+db)**2+add)**3        + (r+da-db)/SQRT((r+da-db)**2+add)**3    &
            -2.0_dp*(r-da)/SQRT((r-da)**2+db**2+add)**3 - 2.0_dp*(r-db)/SQRT((r-db)**2+da**2+add)**3 &
            -2.0_dp*(r+da)/SQRT((r+da)**2+db**2+add)**3 - 2.0_dp*(r+db)/SQRT((r+db)**2+da**2+add)**3 &
            +2.0_dp*r/SQRT(r**2+(da-db)**2+add)**3      + 2.0_dp*r/SQRT(r**2+(da+db)**2+add)**3
       xyxy = &
            +4.0_dp*r/SQRT(r**2+(da-db)**2+add)**3 + 4.0_dp*r/SQRT(r**2+(da+db)**2+add)**3 &
            -8.0_dp*r/SQRT(r**2+da**2+db**2+add)**3
       charg = -(zzzz*0.0625_dp - xyxy*0.015625_dp)*fact
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 1) THEN
       ab = db/SQRT(2.0_dp)
       dxqxz = &
            -2.0_dp*(r-ab)/SQRT((r-ab)**2+(da-ab)**2+add)**3 + 2.0_dp*(r+ab)/SQRT((r+ab)**2+(da-ab)**2+add)**3 &
            +2.0_dp*(r-ab)/SQRT((r-ab)**2+(da+ab)**2+add)**3 - 2.0_dp*(r+ab)/SQRT((r+ab)**2+(da+ab)**2+add)**3
       charg = -dxqxz*0.125_dp*fact
       RETURN
    END IF
    ! ZX - ZX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 1) THEN
       aa = da/SQRT (2.0_dp)
       ab = db/SQRT (2.0_dp)
       qxzqxz = &
            +2.0_dp*(r+aa-ab)/SQRT((r+aa-ab)**2+(aa-ab)**2+add)**3 - 2.0_dp*(r+aa+ab)/SQRT((r+aa+ab)**2+(aa-ab)**2+add)**3 &
            -2.0_dp*(r-aa-ab)/SQRT((r-aa-ab)**2+(aa-ab)**2+add)**3 + 2.0_dp*(r-aa+ab)/SQRT((r-aa+ab)**2+(aa-ab)**2+add)**3 &
            -2.0_dp*(r+aa-ab)/SQRT((r+aa-ab)**2+(aa+ab)**2+add)**3 + 2.0_dp*(r+aa+ab)/SQRT((r+aa+ab)**2+(aa+ab)**2+add)**3 &
            +2.0_dp*(r-aa-ab)/SQRT((r-aa-ab)**2+(aa+ab)**2+add)**3 - 2.0_dp*(r-aa+ab)/SQRT((r-aa+ab)**2+(aa+ab)**2+add)**3
       charg = -qxzqxz*0.0625_dp*fact
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 2) THEN
       xyxy = 4.0_dp*r/SQRT(r**2+(da-db)**2+add)**3 + 4.0_dp*r/SQRT(r**2+(da+db)**2+add)**3 - 8.0_dp*r/SQRT(r**2+da**2+db**2+add)**3
       charg = -xyxy*0.0625_dp*fact
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION dcharg_int_ri

! *****************************************************************************
!> \brief Derivatives of the interaction function between two point-charge
!>        configurations (MNDO/d). This derivative takes into account only the
!>        tapering term
!>        Rotational invariant property built-in in the quadrupole definition
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
  FUNCTION dcharg_int_ri_fs(r, l1_i, l2_i, m, da_i, db_i, add0, fact_screen, error) RESULT(charg)
    REAL(KIND=dp), INTENT(in)                :: r
    INTEGER, INTENT(in)                      :: l1_i, l2_i, m
    REAL(KIND=dp), INTENT(in)                :: da_i, db_i, add0, fact_screen
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: charg

    CHARACTER(len=*), PARAMETER :: routineN = 'dcharg_int_ri_fs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l1, l2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, ab, add, da, db, dxdx, &
                                                dxqxz, dzdz, dzqzz, fact, &
                                                qqzz, qxzqxz, xyxy, zzzz

    failure = .FALSE.
    IF      ( l1_i  <= l2_i) THEN
       l1 = l1_i
       l2 = l2_i
       da = da_i
       db = db_i
       fact  = 1.0_dp
    ELSE IF ( l1_i  > l2_i ) THEN
       l1 = l2_i
       l2 = l1_i
       da = db_i
       db = da_i
       fact  = (-1.0_dp)**(l1+l2)
    END IF
    charg = 0.0_dp
    add   = add0*fact_screen
    ! The 0.5 factor handles the derivative of the SQRT
    fact  = fact * 0.5_dp
    ! Q - Q.
    IF (l1 == 0 .AND. l2 == 0) THEN
       charg = add0/SQRT(r**2+add)**3
       charg = -fact*charg
       RETURN
    END IF
    ! Q - Z.
    IF (l1 == 0 .AND. l2 == 1) THEN
       charg =  add0/SQRT((r+db)**2+add)**3 - add0/SQRT((r-db)**2+add)**3
       charg = -charg*0.5_dp*fact
       RETURN
    END IF
    ! Z - Z.
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 0) THEN
       dzdz  = &
            +add0/SQRT((r+da-db)**2+add)**3 + add0/SQRT((r-da+db)**2+add)**3 &
            -add0/SQRT((r-da-db)**2+add)**3 - add0/SQRT((r+da+db)**2+add)**3
       charg = -dzdz*0.25_dp*fact
       RETURN
    END IF
    ! X - X
    IF (l1 == 1 .AND. l2 == 1 .AND. m == 1) THEN
       dxdx  = 2.0_dp*add0/SQRT(r**2+(da-db)**2+add)**3 - 2.0_dp*add0/SQRT(r**2+(da+db)**2+add)**3
       charg = -dxdx*0.25_dp*fact
       RETURN
    END IF
    ! Q - ZZ
    IF (l1 == 0 .AND. l2 == 2) THEN
       qqzz  = add0/SQRT((r-db)**2+add)**3 - 2.0_dp*add0/SQRT(r**2+db**2+add)**3 + add0/SQRT((r+db)**2+add)**3
       charg = -qqzz*0.25_dp*fact
       RETURN
    END IF
    ! Z - ZZ
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 0) THEN
       dzqzz = &
            +add0/SQRT((r-da-db)**2+add)**3        - 2.0_dp*add0/SQRT((r-da)**2+db**2+add)**3 &
            +add0/SQRT((r+db-da)**2+add)**3        -        add0/SQRT((r-db+da)**2+add)**3    &
            +2.0_dp*add0/SQRT((r+da)**2+db**2+add)**3 -     add0/SQRT((r+da+db)**2+add)**3
       charg = -dzqzz*0.125_dp*fact
       RETURN
    END IF
    ! ZZ - ZZ
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 0) THEN
       zzzz = &
            +add0/SQRT((r-da-db)**2+add)**3        + add0/SQRT((r+da+db)**2+add)**3    &
            +add0/SQRT((r-da+db)**2+add)**3        + add0/SQRT((r+da-db)**2+add)**3    &
            -2.0_dp*add0/SQRT((r-da)**2+db**2+add)**3 - 2.0_dp*add0/SQRT((r-db)**2+da**2+add)**3 &
            -2.0_dp*add0/SQRT((r+da)**2+db**2+add)**3 - 2.0_dp*add0/SQRT((r+db)**2+da**2+add)**3 &
            +2.0_dp*add0/SQRT(r**2+(da-db)**2+add)**3 + 2.0_dp*add0/SQRT(r**2+(da+db)**2+add)**3
       xyxy = &
            +4.0_dp*add0/SQRT(r**2+(da-db)**2+add)**3 + 4.0_dp*add0/SQRT(r**2+(da+db)**2+add)**3 &
            -8.0_dp*add0/SQRT(r**2+da**2+db**2+add)**3
       charg = -(zzzz*0.0625_dp - xyxy*0.015625_dp)*fact
       RETURN
    END IF
    ! X - ZX
    IF (l1 == 1 .AND. l2 == 2 .AND. m == 1) THEN
       ab = db/SQRT(2.0_dp)
       dxqxz = &
            -2.0_dp*add0/SQRT((r-ab)**2+(da-ab)**2+add)**3 + 2.0_dp*add0/SQRT((r+ab)**2+(da-ab)**2+add)**3 &
            +2.0_dp*add0/SQRT((r-ab)**2+(da+ab)**2+add)**3 - 2.0_dp*add0/SQRT((r+ab)**2+(da+ab)**2+add)**3
       charg = -dxqxz*0.125_dp*fact
       RETURN
    END IF
    ! ZX - ZX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 1) THEN
       aa = da/SQRT (2.0_dp)
       ab = db/SQRT (2.0_dp)
       qxzqxz = &
            +2.0_dp*add0/SQRT((r+aa-ab)**2+(aa-ab)**2+add)**3 - 2.0_dp*add0/SQRT((r+aa+ab)**2+(aa-ab)**2+add)**3 &
            -2.0_dp*add0/SQRT((r-aa-ab)**2+(aa-ab)**2+add)**3 + 2.0_dp*add0/SQRT((r-aa+ab)**2+(aa-ab)**2+add)**3 &
            -2.0_dp*add0/SQRT((r+aa-ab)**2+(aa+ab)**2+add)**3 + 2.0_dp*add0/SQRT((r+aa+ab)**2+(aa+ab)**2+add)**3 &
            +2.0_dp*add0/SQRT((r-aa-ab)**2+(aa+ab)**2+add)**3 - 2.0_dp*add0/SQRT((r-aa+ab)**2+(aa+ab)**2+add)**3
       charg = -qxzqxz*0.0625_dp*fact
       RETURN
    END IF
    ! XX - XX
    IF (l1 == 2 .AND. l2 == 2 .AND. m == 2) THEN
       xyxy = 4.0_dp*add0/SQRT(r**2+(da-db)**2+add)**3   + 4.0_dp*add0/SQRT(r**2+(da+db)**2+add)**3 -&
              8.0_dp*add0/SQRT(r**2+da**2+db**2+add)**3
       charg = -xyxy*0.0625_dp*fact
       RETURN
    END IF
    ! We should NEVER reach this point
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END FUNCTION dcharg_int_ri_fs

! *****************************************************************************
!> \brief Computes the general rotation matrices for the integrals
!>        between I and J (J>=I)
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  RECURSIVE SUBROUTINE rotmat (sepi, sepj, rjiv, r, ij_matrix, do_derivatives,&
       do_invert, debug_invert, error)
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
    REAL(KIND=dp), DIMENSION(3, 5, 5)        :: d_d
    REAL(KIND=dp), DIMENSION(5, 5)           :: d

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(ij_matrix),cp_failure_level,routineP,error,failure)
    IF (PRESENT(do_invert)) do_invert = .FALSE.
    IF ((sepi%natorb>1).OR.(sepj%natorb>1)) THEN
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
             imap(1)   = 3
             imap(2)   = 2
             imap(3)   = 1
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
       END IF
       IF (dbg_inv) THEN
          ! When debugging the derivatives of the rotation matrix we must possibly force
          ! the inversion of the reference frame
          CPPostcondition(.NOT.do_derivatives,cp_failure_level,routineP,error,failure)
          imap(1)   = 3
          imap(2)   = 2
          imap(3)   = 1
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
       END IF
       ! Evaluate Derivatives if required
       IF (do_derivatives) THEN
          ! This condition is necessary because derivatives uses the invertion of the
          ! axis for treating the divergence point along z-axis
          CPPostcondition(eval,cp_failure_level,routineP,error,failure)
          x11_d = 0.0_dp; x11_d(1) = 1.0_dp
          x22_d = 0.0_dp; x22_d(2) = 1.0_dp
          x33_d = 0.0_dp; x33_d(3) = 1.0_dp
          r_d   = (/x11,x22,x33/)/r
          b_d   = 2.0_dp*(x11*x11_d + x22*x22_d)
          sqb_d = (0.5_dp/sqb)*b_d
          r2i   = 1.0_dp/r**2
          sqb2i = 1.0_dp/sqb**2
          cb_d  = (x33_d*r-x33*r_d)*r2i
          ca_d  = (x11_d*sqb-x11*sqb_d)*sqb2i
          sa_d  = (x22_d*sqb-x22*sqb_d)*sqb2i
          sb_d  = (sqb_d*r-sqb*r_d)*r2i
          ! Calculate derivatives of rotation matrix elements
          p_d(:, 1, 1) =  ca_d * sb + ca * sb_d
          p_d(:, 2, 1) =  ca_d * cb + ca * cb_d
          p_d(:, 3, 1) = -sa_d
          p_d(:, 1, 2) =  sa_d * sb + sa * sb_d
          p_d(:, 2, 2) =  sa_d * cb + sa * cb_d
          p_d(:, 3, 2) =  ca_d
          p_d(:, 1, 3) =  cb_d
          p_d(:, 2, 3) = -sb_d
          p_d(:, 3, 3) =  0.0_dp
          IF (sepi%dorb.OR.sepj%dorb) THEN
             ca2_d = 2.0_dp * ca * ca_d
             cb2_d = 2.0_dp * cb * cb_d
             sb2_d = 2.0_dp * sb * sb_d
             c2a_d = 2.0_dp * ca2_d
             c2b_d = 2.0_dp * cb2_d
             s2a_d = 2.0_dp * ( sa_d * ca + sa * ca_d )
             s2b_d = 2.0_dp * ( sb_d * cb + sb * cb_d )
             d_d(:, 1, 1) =  pt5sq3 * ( c2a_d * sb2 + c2a * sb2_d )
             d_d(:, 2, 1) =  0.5_dp * ( c2a_d * s2b + c2a * s2b_d )
             d_d(:, 3, 1) = -s2a_d * sb - s2a * sb_d
             d_d(:, 4, 1) =  c2a_d * (cb2+0.5_dp*sb2) + c2a * (cb2_d+0.5_dp*sb2_d)
             d_d(:, 5, 1) = -s2a_d * cb - s2a * cb_d
             d_d(:, 1, 2) =  pt5sq3 * ( ca_d * s2b + ca * s2b_d )
             d_d(:, 2, 2) =  ca_d * c2b + ca * c2b_d
             d_d(:, 3, 2) = -sa_d * cb - sa * cb_d
             d_d(:, 4, 2) = -0.5_dp * ( ca_d * s2b + ca * s2b_d )
             d_d(:, 5, 2) =  sa_d * sb + sa * sb_d
             d_d(:, 1, 3) =  cb2_d - 0.5_dp * sb2_d
             d_d(:, 2, 3) = -pt5sq3 * s2b_d
             d_d(:, 3, 3) =  0.0_dp
             d_d(:, 4, 3) =  pt5sq3 * sb2_d
             d_d(:, 5, 3) =  0.0_dp
             d_d(:, 1, 4) =  pt5sq3 * ( sa_d * s2b + sa * s2b_d )
             d_d(:, 2, 4) =  sa_d * c2b + sa * c2b_d
             d_d(:, 3, 4) =  ca_d * cb + ca * cb_d
             d_d(:, 4, 4) = -0.5_dp * ( sa_d * s2b + sa * s2b_d )
             d_d(:, 5, 4) = -ca_d * sb -ca * sb_d
             d_d(:, 1, 5) =  pt5sq3 * ( s2a_d * sb2 + s2a * sb2_d )
             d_d(:, 2, 5) =  0.5_dp * ( s2a_d * s2b + s2a * s2b_d )
             d_d(:, 3, 5) =  c2a_d * sb + c2a * sb_d
             d_d(:, 4, 5) =  s2a_d * (cb2+0.5_dp*sb2) + s2a * (cb2_d+0.5_dp*sb2_d)
             d_d(:, 5, 5) =  c2a_d * cb + c2a * cb_d
          END IF
          !  Derivative for Rotation Elements for : S-P
          DO k = 1, 3
             DO l = 1, 3
                ij_matrix%sp_d(:, k, l) = p_d(:, k, l)
             END DO
          END DO
          !  Derivative for Rotation Elements for : P-P
          DO k = 1, 3
             ij_matrix%pp_d(:, 1, k, k) = p_d(:, k, 1) * p(k, 1) +  p(k, 1) * p_d(:, k, 1)
             ij_matrix%pp_d(:, 2, k, k) = p_d(:, k, 2) * p(k, 2) +  p(k, 2) * p_d(:, k, 2)
             ij_matrix%pp_d(:, 3, k, k) = p_d(:, k, 3) * p(k, 3) +  p(k, 3) * p_d(:, k, 3)
             ij_matrix%pp_d(:, 4, k, k) = p_d(:, k, 1) * p(k, 2) +  p(k, 1) * p_d(:, k, 2)
             ij_matrix%pp_d(:, 5, k, k) = p_d(:, k, 1) * p(k, 3) +  p(k, 1) * p_d(:, k, 3)
             ij_matrix%pp_d(:, 6, k, k) = p_d(:, k, 2) * p(k, 3) +  p(k, 2) * p_d(:, k, 3)
             IF (k /= 1) THEN
                DO l = 1, k - 1
                   ij_matrix%pp_d(:, 1, k, l) = 2.0_dp * ( p_d(:, k, 1) * p(l, 1) + p(k, 1) * p_d(:, l, 1) )
                   ij_matrix%pp_d(:, 2, k, l) = 2.0_dp * ( p_d(:, k, 2) * p(l, 2) + p(k, 2) * p_d(:, l, 2) )
                   ij_matrix%pp_d(:, 3, k, l) = 2.0_dp * ( p_d(:, k, 3) * p(l, 3) + p(k, 3) * p_d(:, l, 3) )
                   ij_matrix%pp_d(:, 4, k, l) = ( p_d(:, k, 1) * p(l, 2) + p(k, 1) * p_d(:, l, 2) ) + &
                                                ( p_d(:, k, 2) * p(l, 1) + p(k, 2) * p_d(:, l, 1) )
                   ij_matrix%pp_d(:, 5, k, l) = ( p_d(:, k, 1) * p(l, 3) + p(k, 1) * p_d(:, l, 3) ) + &
                                                ( p_d(:, k, 3) * p(l, 1) + p(k, 3) * p_d(:, l, 1) )
                   ij_matrix%pp_d(:, 6, k, l) = ( p_d(:, k, 2) * p(l, 3) + p(k, 2) * p_d(:, l, 3) ) + &
                                                ( p_d(:, k, 3) * p(l, 2) + p(k, 3) * p_d(:, l, 2) )
                END DO
             END IF
          END DO
          IF (sepi%dorb.OR.sepj%dorb) THEN
             !  Rotation Elements for : S-D
             DO k = 1, 5
                DO l = 1, 5
                   ij_matrix%sd_d(:, k, l) = d_d(:, k, l)
                END DO
             END DO
             !  Rotation Elements for : D-P
             DO k = 1, 5
                DO l = 1, 3
                   ij_matrix%pd_d(:,  1, k, l) = ( d_d(:, k, 1) * p(l, 1) + d(k, 1) * p_d(:, l, 1) )
                   ij_matrix%pd_d(:,  2, k, l) = ( d_d(:, k, 1) * p(l, 2) + d(k, 1) * p_d(:, l, 2) )
                   ij_matrix%pd_d(:,  3, k, l) = ( d_d(:, k, 1) * p(l, 3) + d(k, 1) * p_d(:, l, 3) )
                   ij_matrix%pd_d(:,  4, k, l) = ( d_d(:, k, 2) * p(l, 1) + d(k, 2) * p_d(:, l, 1) )
                   ij_matrix%pd_d(:,  5, k, l) = ( d_d(:, k, 2) * p(l, 2) + d(k, 2) * p_d(:, l, 2) )
                   ij_matrix%pd_d(:,  6, k, l) = ( d_d(:, k, 2) * p(l, 3) + d(k, 2) * p_d(:, l, 3) )
                   ij_matrix%pd_d(:,  7, k, l) = ( d_d(:, k, 3) * p(l, 1) + d(k, 3) * p_d(:, l, 1) )
                   ij_matrix%pd_d(:,  8, k, l) = ( d_d(:, k, 3) * p(l, 2) + d(k, 3) * p_d(:, l, 2) )
                   ij_matrix%pd_d(:,  9, k, l) = ( d_d(:, k, 3) * p(l, 3) + d(k, 3) * p_d(:, l, 3) )
                   ij_matrix%pd_d(:, 10, k, l) = ( d_d(:, k, 4) * p(l, 1) + d(k, 4) * p_d(:, l, 1) )
                   ij_matrix%pd_d(:, 11, k, l) = ( d_d(:, k, 4) * p(l, 2) + d(k, 4) * p_d(:, l, 2) )
                   ij_matrix%pd_d(:, 12, k, l) = ( d_d(:, k, 4) * p(l, 3) + d(k, 4) * p_d(:, l, 3) )
                   ij_matrix%pd_d(:, 13, k, l) = ( d_d(:, k, 5) * p(l, 1) + d(k, 5) * p_d(:, l, 1) )
                   ij_matrix%pd_d(:, 14, k, l) = ( d_d(:, k, 5) * p(l, 2) + d(k, 5) * p_d(:, l, 2) )
                   ij_matrix%pd_d(:, 15, k, l) = ( d_d(:, k, 5) * p(l, 3) + d(k, 5) * p_d(:, l, 3) )
                END DO
             END DO
             !  Rotation Elements for : D-D
             DO k = 1, 5
                ij_matrix%dd_d(:,  1, k, k) = ( d_d(:, k, 1) * d(k, 1) + d(k, 1) * d_d(:, k, 1) )
                ij_matrix%dd_d(:,  2, k, k) = ( d_d(:, k, 2) * d(k, 2) + d(k, 2) * d_d(:, k, 2) )
                ij_matrix%dd_d(:,  3, k, k) = ( d_d(:, k, 3) * d(k, 3) + d(k, 3) * d_d(:, k, 3) )
                ij_matrix%dd_d(:,  4, k, k) = ( d_d(:, k, 4) * d(k, 4) + d(k, 4) * d_d(:, k, 4) )
                ij_matrix%dd_d(:,  5, k, k) = ( d_d(:, k, 5) * d(k, 5) + d(k, 5) * d_d(:, k, 5) )
                ij_matrix%dd_d(:,  6, k, k) = ( d_d(:, k, 1) * d(k, 2) + d(k, 1) * d_d(:, k, 2) )
                ij_matrix%dd_d(:,  7, k, k) = ( d_d(:, k, 1) * d(k, 3) + d(k, 1) * d_d(:, k, 3) )
                ij_matrix%dd_d(:,  8, k, k) = ( d_d(:, k, 2) * d(k, 3) + d(k, 2) * d_d(:, k, 3) )
                ij_matrix%dd_d(:,  9, k, k) = ( d_d(:, k, 1) * d(k, 4) + d(k, 1) * d_d(:, k, 4) )
                ij_matrix%dd_d(:, 10, k, k) = ( d_d(:, k, 2) * d(k, 4) + d(k, 2) * d_d(:, k, 4) )
                ij_matrix%dd_d(:, 11, k, k) = ( d_d(:, k, 3) * d(k, 4) + d(k, 3) * d_d(:, k, 4) )
                ij_matrix%dd_d(:, 12, k, k) = ( d_d(:, k, 1) * d(k, 5) + d(k, 1) * d_d(:, k, 5) )
                ij_matrix%dd_d(:, 13, k, k) = ( d_d(:, k, 2) * d(k, 5) + d(k, 2) * d_d(:, k, 5) )
                ij_matrix%dd_d(:, 14, k, k) = ( d_d(:, k, 3) * d(k, 5) + d(k, 3) * d_d(:, k, 5) )
                ij_matrix%dd_d(:, 15, k, k) = ( d_d(:, k, 4) * d(k, 5) + d(k, 4) * d_d(:, k, 5) )
                IF (k /= 1) THEN
                   DO l = 1, k - 1
                      ij_matrix%dd_d(:,  1, k, l) = 2.0_dp  * ( d_d(:, k, 1) * d(l, 1) + d(k, 1) * d_d(:, l, 1) )
                      ij_matrix%dd_d(:,  2, k, l) = 2.0_dp  * ( d_d(:, k, 2) * d(l, 2) + d(k, 2) * d_d(:, l, 2) )
                      ij_matrix%dd_d(:,  3, k, l) = 2.0_dp  * ( d_d(:, k, 3) * d(l, 3) + d(k, 3) * d_d(:, l, 3) )
                      ij_matrix%dd_d(:,  4, k, l) = 2.0_dp  * ( d_d(:, k, 4) * d(l, 4) + d(k, 4) * d_d(:, l, 4) )
                      ij_matrix%dd_d(:,  5, k, l) = 2.0_dp  * ( d_d(:, k, 5) * d(l, 5) + d(k, 5) * d_d(:, l, 5) )
                      ij_matrix%dd_d(:,  6, k, l) = ( d_d(:, k, 1) * d(l, 2) + d(k, 1) * d_d(:, l, 2) ) + &
                                                    ( d_d(:, k, 2) * d(l, 1) + d(k, 2) * d_d(:, l, 1) )
                      ij_matrix%dd_d(:,  7, k, l) = ( d_d(:, k, 1) * d(l, 3) + d(k, 1) * d_d(:, l, 3) ) + &
                                                    ( d_d(:, k, 3) * d(l, 1) + d(k, 3) * d_d(:, l, 1) )
                      ij_matrix%dd_d(:,  8, k, l) = ( d_d(:, k, 2) * d(l, 3) + d(k, 2) * d_d(:, l, 3) ) + &
                                                    ( d_d(:, k, 3) * d(l, 2) + d(k, 3) * d_d(:, l, 2) )
                      ij_matrix%dd_d(:,  9, k, l) = ( d_d(:, k, 1) * d(l, 4) + d(k, 1) * d_d(:, l, 4) ) + &
                                                    ( d_d(:, k, 4) * d(l, 1) + d(k, 4) * d_d(:, l, 1) )
                      ij_matrix%dd_d(:, 10, k, l) = ( d_d(:, k, 2) * d(l, 4) + d(k, 2) * d_d(:, l, 4) ) + &
                                                    ( d_d(:, k, 4) * d(l, 2) + d(k, 4) * d_d(:, l, 2) )
                      ij_matrix%dd_d(:, 11, k, l) = ( d_d(:, k, 3) * d(l, 4) + d(k, 3) * d_d(:, l, 4) ) + &
                                                    ( d_d(:, k, 4) * d(l, 3) + d(k, 4) * d_d(:, l, 3) )
                      ij_matrix%dd_d(:, 12, k, l) = ( d_d(:, k, 1) * d(l, 5) + d(k, 1) * d_d(:, l, 5) ) + &
                                                    ( d_d(:, k, 5) * d(l, 1) + d(k, 5) * d_d(:, l, 1) )
                      ij_matrix%dd_d(:, 13, k, l) = ( d_d(:, k, 2) * d(l, 5) + d(k, 2) * d_d(:, l, 5) ) + &
                                                    ( d_d(:, k, 5) * d(l, 2) + d(k, 5) * d_d(:, l, 2) )
                      ij_matrix%dd_d(:, 14, k, l) = ( d_d(:, k, 3) * d(l, 5) + d(k, 3) * d_d(:, l, 5) ) + &
                                                    ( d_d(:, k, 5) * d(l, 3) + d(k, 5) * d_d(:, l, 3) )
                      ij_matrix%dd_d(:, 15, k, l) = ( d_d(:, k, 4) * d(l, 5) + d(k, 4) * d_d(:, l, 5) ) + &
                                                    ( d_d(:, k, 5) * d(l, 4) + d(k, 5) * d_d(:, l, 4) )
                   END DO
                END IF
             END DO
          END IF
          IF (debug_this_module) THEN
             CALL check_rotmat_der(sepi, sepj, rjiv, ij_matrix, do_invert=do_invert, error=error)
          END IF
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
  RECURSIVE SUBROUTINE rot_2el_2c_first (sepi, sepj, rijv, se_int_control, se_taper,&
       invert, ii, kk, rep, logv, ij_matrix, v, lgrad, rep_d, v_d, logv_d, drij,&
       error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rijv
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
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
    REAL(KIND=dp), DIMENSION(3, 45, 45), &
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
    DO k = 1, limkl
       DO i = 1, 45
          logv(i, k) = .FALSE.
          v(i, k)    = 0.0_dp
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
       DO k = 1, limkl
          DO i = 1, 45
             logv_d(i, k) = .FALSE.
             v_d(1, i, k) = 0.0_dp
             v_d(2, i, k) = 0.0_dp
             v_d(3, i, k) = 0.0_dp
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
                         v_d(1, ij, 1) = wrepp_d(1)
                         v_d(2, ij, 1) = wrepp_d(2)
                         v_d(3, ij, 1) = wrepp_d(3)
                      ELSE IF (mm == 2) THEN
                         k = k1 - 1
                         v_d(1, ij,  2) = v_d(1, ij,  2) + ij_matrix%sp_d(1, k, 1) * wrepp + ij_matrix%sp(k, 1) * wrepp_d(1)
                         v_d(1, ij,  4) = v_d(1, ij,  4) + ij_matrix%sp_d(1, k, 2) * wrepp + ij_matrix%sp(k, 2) * wrepp_d(1)
                         v_d(1, ij,  7) = v_d(1, ij,  7) + ij_matrix%sp_d(1, k, 3) * wrepp + ij_matrix%sp(k, 3) * wrepp_d(1)

                         v_d(2, ij,  2) = v_d(2, ij,  2) + ij_matrix%sp_d(2, k, 1) * wrepp + ij_matrix%sp(k, 1) * wrepp_d(2)
                         v_d(2, ij,  4) = v_d(2, ij,  4) + ij_matrix%sp_d(2, k, 2) * wrepp + ij_matrix%sp(k, 2) * wrepp_d(2)
                         v_d(2, ij,  7) = v_d(2, ij,  7) + ij_matrix%sp_d(2, k, 3) * wrepp + ij_matrix%sp(k, 3) * wrepp_d(2)

                         v_d(3, ij,  2) = v_d(3, ij,  2) + ij_matrix%sp_d(3, k, 1) * wrepp + ij_matrix%sp(k, 1) * wrepp_d(3)
                         v_d(3, ij,  4) = v_d(3, ij,  4) + ij_matrix%sp_d(3, k, 2) * wrepp + ij_matrix%sp(k, 2) * wrepp_d(3)
                         v_d(3, ij,  7) = v_d(3, ij,  7) + ij_matrix%sp_d(3, k, 3) * wrepp + ij_matrix%sp(k, 3) * wrepp_d(3)
                      ELSE IF (mm == 3) THEN
                         k = k1 - 1
                         l = l1 - 1
                         v_d(1, ij,  3) = v_d(1, ij,  3) + ij_matrix%pp_d(1, 1, k, l) * wrepp + ij_matrix%pp(1, k, l) * wrepp_d(1)
                         v_d(1, ij,  6) = v_d(1, ij,  6) + ij_matrix%pp_d(1, 2, k, l) * wrepp + ij_matrix%pp(2, k, l) * wrepp_d(1)
                         v_d(1, ij, 10) = v_d(1, ij, 10) + ij_matrix%pp_d(1, 3, k, l) * wrepp + ij_matrix%pp(3, k, l) * wrepp_d(1)
                         v_d(1, ij,  5) = v_d(1, ij,  5) + ij_matrix%pp_d(1, 4, k, l) * wrepp + ij_matrix%pp(4, k, l) * wrepp_d(1)
                         v_d(1, ij,  8) = v_d(1, ij,  8) + ij_matrix%pp_d(1, 5, k, l) * wrepp + ij_matrix%pp(5, k, l) * wrepp_d(1)
                         v_d(1, ij,  9) = v_d(1, ij,  9) + ij_matrix%pp_d(1, 6, k, l) * wrepp + ij_matrix%pp(6, k, l) * wrepp_d(1)

                         v_d(2, ij,  3) = v_d(2, ij,  3) + ij_matrix%pp_d(2, 1, k, l) * wrepp + ij_matrix%pp(1, k, l) * wrepp_d(2)
                         v_d(2, ij,  6) = v_d(2, ij,  6) + ij_matrix%pp_d(2, 2, k, l) * wrepp + ij_matrix%pp(2, k, l) * wrepp_d(2)
                         v_d(2, ij, 10) = v_d(2, ij, 10) + ij_matrix%pp_d(2, 3, k, l) * wrepp + ij_matrix%pp(3, k, l) * wrepp_d(2)
                         v_d(2, ij,  5) = v_d(2, ij,  5) + ij_matrix%pp_d(2, 4, k, l) * wrepp + ij_matrix%pp(4, k, l) * wrepp_d(2)
                         v_d(2, ij,  8) = v_d(2, ij,  8) + ij_matrix%pp_d(2, 5, k, l) * wrepp + ij_matrix%pp(5, k, l) * wrepp_d(2)
                         v_d(2, ij,  9) = v_d(2, ij,  9) + ij_matrix%pp_d(2, 6, k, l) * wrepp + ij_matrix%pp(6, k, l) * wrepp_d(2)

                         v_d(3, ij,  3) = v_d(3, ij,  3) + ij_matrix%pp_d(3, 1, k, l) * wrepp + ij_matrix%pp(1, k, l) * wrepp_d(3)
                         v_d(3, ij,  6) = v_d(3, ij,  6) + ij_matrix%pp_d(3, 2, k, l) * wrepp + ij_matrix%pp(2, k, l) * wrepp_d(3)
                         v_d(3, ij, 10) = v_d(3, ij, 10) + ij_matrix%pp_d(3, 3, k, l) * wrepp + ij_matrix%pp(3, k, l) * wrepp_d(3)
                         v_d(3, ij,  5) = v_d(3, ij,  5) + ij_matrix%pp_d(3, 4, k, l) * wrepp + ij_matrix%pp(4, k, l) * wrepp_d(3)
                         v_d(3, ij,  8) = v_d(3, ij,  8) + ij_matrix%pp_d(3, 5, k, l) * wrepp + ij_matrix%pp(5, k, l) * wrepp_d(3)
                         v_d(3, ij,  9) = v_d(3, ij,  9) + ij_matrix%pp_d(3, 6, k, l) * wrepp + ij_matrix%pp(6, k, l) * wrepp_d(3)
                      ELSE IF (mm == 4) THEN
                         k = k1 - 4
                         v_d(1, ij, 11) = v_d(1, ij, 11) + ij_matrix%sd_d(1, k, 1) * wrepp + ij_matrix%sd(k, 1) * wrepp_d(1)
                         v_d(1, ij, 16) = v_d(1, ij, 16) + ij_matrix%sd_d(1, k, 2) * wrepp + ij_matrix%sd(k, 2) * wrepp_d(1)
                         v_d(1, ij, 22) = v_d(1, ij, 22) + ij_matrix%sd_d(1, k, 3) * wrepp + ij_matrix%sd(k, 3) * wrepp_d(1)
                         v_d(1, ij, 29) = v_d(1, ij, 29) + ij_matrix%sd_d(1, k, 4) * wrepp + ij_matrix%sd(k, 4) * wrepp_d(1)
                         v_d(1, ij, 37) = v_d(1, ij, 37) + ij_matrix%sd_d(1, k, 5) * wrepp + ij_matrix%sd(k, 5) * wrepp_d(1)

                         v_d(2, ij, 11) = v_d(2, ij, 11) + ij_matrix%sd_d(2, k, 1) * wrepp + ij_matrix%sd(k, 1) * wrepp_d(2)
                         v_d(2, ij, 16) = v_d(2, ij, 16) + ij_matrix%sd_d(2, k, 2) * wrepp + ij_matrix%sd(k, 2) * wrepp_d(2)
                         v_d(2, ij, 22) = v_d(2, ij, 22) + ij_matrix%sd_d(2, k, 3) * wrepp + ij_matrix%sd(k, 3) * wrepp_d(2)
                         v_d(2, ij, 29) = v_d(2, ij, 29) + ij_matrix%sd_d(2, k, 4) * wrepp + ij_matrix%sd(k, 4) * wrepp_d(2)
                         v_d(2, ij, 37) = v_d(2, ij, 37) + ij_matrix%sd_d(2, k, 5) * wrepp + ij_matrix%sd(k, 5) * wrepp_d(2)

                         v_d(3, ij, 11) = v_d(3, ij, 11) + ij_matrix%sd_d(3, k, 1) * wrepp + ij_matrix%sd(k, 1) * wrepp_d(3)
                         v_d(3, ij, 16) = v_d(3, ij, 16) + ij_matrix%sd_d(3, k, 2) * wrepp + ij_matrix%sd(k, 2) * wrepp_d(3)
                         v_d(3, ij, 22) = v_d(3, ij, 22) + ij_matrix%sd_d(3, k, 3) * wrepp + ij_matrix%sd(k, 3) * wrepp_d(3)
                         v_d(3, ij, 29) = v_d(3, ij, 29) + ij_matrix%sd_d(3, k, 4) * wrepp + ij_matrix%sd(k, 4) * wrepp_d(3)
                         v_d(3, ij, 37) = v_d(3, ij, 37) + ij_matrix%sd_d(3, k, 5) * wrepp + ij_matrix%sd(k, 5) * wrepp_d(3)
                      ELSE IF (mm == 5) THEN
                         k = k1 - 4
                         l = l1 - 1
                         v_d(1, ij, 12) = v_d(1, ij, 12) + ij_matrix%pd_d(1,  1, k, l) * wrepp + ij_matrix%pd( 1, k, l) * wrepp_d(1)
                         v_d(1, ij, 13) = v_d(1, ij, 13) + ij_matrix%pd_d(1,  2, k, l) * wrepp + ij_matrix%pd( 2, k, l) * wrepp_d(1)
                         v_d(1, ij, 14) = v_d(1, ij, 14) + ij_matrix%pd_d(1,  3, k, l) * wrepp + ij_matrix%pd( 3, k, l) * wrepp_d(1)
                         v_d(1, ij, 17) = v_d(1, ij, 17) + ij_matrix%pd_d(1,  4, k, l) * wrepp + ij_matrix%pd( 4, k, l) * wrepp_d(1)
                         v_d(1, ij, 18) = v_d(1, ij, 18) + ij_matrix%pd_d(1,  5, k, l) * wrepp + ij_matrix%pd( 5, k, l) * wrepp_d(1)
                         v_d(1, ij, 19) = v_d(1, ij, 19) + ij_matrix%pd_d(1,  6, k, l) * wrepp + ij_matrix%pd( 6, k, l) * wrepp_d(1)
                         v_d(1, ij, 23) = v_d(1, ij, 23) + ij_matrix%pd_d(1,  7, k, l) * wrepp + ij_matrix%pd( 7, k, l) * wrepp_d(1)
                         v_d(1, ij, 24) = v_d(1, ij, 24) + ij_matrix%pd_d(1,  8, k, l) * wrepp + ij_matrix%pd( 8, k, l) * wrepp_d(1)
                         v_d(1, ij, 25) = v_d(1, ij, 25) + ij_matrix%pd_d(1,  9, k, l) * wrepp + ij_matrix%pd( 9, k, l) * wrepp_d(1)
                         v_d(1, ij, 30) = v_d(1, ij, 30) + ij_matrix%pd_d(1, 10, k, l) * wrepp + ij_matrix%pd(10, k, l) * wrepp_d(1)
                         v_d(1, ij, 31) = v_d(1, ij, 31) + ij_matrix%pd_d(1, 11, k, l) * wrepp + ij_matrix%pd(11, k, l) * wrepp_d(1)
                         v_d(1, ij, 32) = v_d(1, ij, 32) + ij_matrix%pd_d(1, 12, k, l) * wrepp + ij_matrix%pd(12, k, l) * wrepp_d(1)
                         v_d(1, ij, 38) = v_d(1, ij, 38) + ij_matrix%pd_d(1, 13, k, l) * wrepp + ij_matrix%pd(13, k, l) * wrepp_d(1)
                         v_d(1, ij, 39) = v_d(1, ij, 39) + ij_matrix%pd_d(1, 14, k, l) * wrepp + ij_matrix%pd(14, k, l) * wrepp_d(1)
                         v_d(1, ij, 40) = v_d(1, ij, 40) + ij_matrix%pd_d(1, 15, k, l) * wrepp + ij_matrix%pd(15, k, l) * wrepp_d(1)

                         v_d(2, ij, 12) = v_d(2, ij, 12) + ij_matrix%pd_d(2,  1, k, l) * wrepp + ij_matrix%pd( 1, k, l) * wrepp_d(2)
                         v_d(2, ij, 13) = v_d(2, ij, 13) + ij_matrix%pd_d(2,  2, k, l) * wrepp + ij_matrix%pd( 2, k, l) * wrepp_d(2)
                         v_d(2, ij, 14) = v_d(2, ij, 14) + ij_matrix%pd_d(2,  3, k, l) * wrepp + ij_matrix%pd( 3, k, l) * wrepp_d(2)
                         v_d(2, ij, 17) = v_d(2, ij, 17) + ij_matrix%pd_d(2,  4, k, l) * wrepp + ij_matrix%pd( 4, k, l) * wrepp_d(2)
                         v_d(2, ij, 18) = v_d(2, ij, 18) + ij_matrix%pd_d(2,  5, k, l) * wrepp + ij_matrix%pd( 5, k, l) * wrepp_d(2)
                         v_d(2, ij, 19) = v_d(2, ij, 19) + ij_matrix%pd_d(2,  6, k, l) * wrepp + ij_matrix%pd( 6, k, l) * wrepp_d(2)
                         v_d(2, ij, 23) = v_d(2, ij, 23) + ij_matrix%pd_d(2,  7, k, l) * wrepp + ij_matrix%pd( 7, k, l) * wrepp_d(2)
                         v_d(2, ij, 24) = v_d(2, ij, 24) + ij_matrix%pd_d(2,  8, k, l) * wrepp + ij_matrix%pd( 8, k, l) * wrepp_d(2)
                         v_d(2, ij, 25) = v_d(2, ij, 25) + ij_matrix%pd_d(2,  9, k, l) * wrepp + ij_matrix%pd( 9, k, l) * wrepp_d(2)
                         v_d(2, ij, 30) = v_d(2, ij, 30) + ij_matrix%pd_d(2, 10, k, l) * wrepp + ij_matrix%pd(10, k, l) * wrepp_d(2)
                         v_d(2, ij, 31) = v_d(2, ij, 31) + ij_matrix%pd_d(2, 11, k, l) * wrepp + ij_matrix%pd(11, k, l) * wrepp_d(2)
                         v_d(2, ij, 32) = v_d(2, ij, 32) + ij_matrix%pd_d(2, 12, k, l) * wrepp + ij_matrix%pd(12, k, l) * wrepp_d(2)
                         v_d(2, ij, 38) = v_d(2, ij, 38) + ij_matrix%pd_d(2, 13, k, l) * wrepp + ij_matrix%pd(13, k, l) * wrepp_d(2)
                         v_d(2, ij, 39) = v_d(2, ij, 39) + ij_matrix%pd_d(2, 14, k, l) * wrepp + ij_matrix%pd(14, k, l) * wrepp_d(2)
                         v_d(2, ij, 40) = v_d(2, ij, 40) + ij_matrix%pd_d(2, 15, k, l) * wrepp + ij_matrix%pd(15, k, l) * wrepp_d(2)

                         v_d(3, ij, 12) = v_d(3, ij, 12) + ij_matrix%pd_d(3,  1, k, l) * wrepp + ij_matrix%pd( 1, k, l) * wrepp_d(3)
                         v_d(3, ij, 13) = v_d(3, ij, 13) + ij_matrix%pd_d(3,  2, k, l) * wrepp + ij_matrix%pd( 2, k, l) * wrepp_d(3)
                         v_d(3, ij, 14) = v_d(3, ij, 14) + ij_matrix%pd_d(3,  3, k, l) * wrepp + ij_matrix%pd( 3, k, l) * wrepp_d(3)
                         v_d(3, ij, 17) = v_d(3, ij, 17) + ij_matrix%pd_d(3,  4, k, l) * wrepp + ij_matrix%pd( 4, k, l) * wrepp_d(3)
                         v_d(3, ij, 18) = v_d(3, ij, 18) + ij_matrix%pd_d(3,  5, k, l) * wrepp + ij_matrix%pd( 5, k, l) * wrepp_d(3)
                         v_d(3, ij, 19) = v_d(3, ij, 19) + ij_matrix%pd_d(3,  6, k, l) * wrepp + ij_matrix%pd( 6, k, l) * wrepp_d(3)
                         v_d(3, ij, 23) = v_d(3, ij, 23) + ij_matrix%pd_d(3,  7, k, l) * wrepp + ij_matrix%pd( 7, k, l) * wrepp_d(3)
                         v_d(3, ij, 24) = v_d(3, ij, 24) + ij_matrix%pd_d(3,  8, k, l) * wrepp + ij_matrix%pd( 8, k, l) * wrepp_d(3)
                         v_d(3, ij, 25) = v_d(3, ij, 25) + ij_matrix%pd_d(3,  9, k, l) * wrepp + ij_matrix%pd( 9, k, l) * wrepp_d(3)
                         v_d(3, ij, 30) = v_d(3, ij, 30) + ij_matrix%pd_d(3, 10, k, l) * wrepp + ij_matrix%pd(10, k, l) * wrepp_d(3)
                         v_d(3, ij, 31) = v_d(3, ij, 31) + ij_matrix%pd_d(3, 11, k, l) * wrepp + ij_matrix%pd(11, k, l) * wrepp_d(3)
                         v_d(3, ij, 32) = v_d(3, ij, 32) + ij_matrix%pd_d(3, 12, k, l) * wrepp + ij_matrix%pd(12, k, l) * wrepp_d(3)
                         v_d(3, ij, 38) = v_d(3, ij, 38) + ij_matrix%pd_d(3, 13, k, l) * wrepp + ij_matrix%pd(13, k, l) * wrepp_d(3)
                         v_d(3, ij, 39) = v_d(3, ij, 39) + ij_matrix%pd_d(3, 14, k, l) * wrepp + ij_matrix%pd(14, k, l) * wrepp_d(3)
                         v_d(3, ij, 40) = v_d(3, ij, 40) + ij_matrix%pd_d(3, 15, k, l) * wrepp + ij_matrix%pd(15, k, l) * wrepp_d(3)
                      ELSE IF (mm == 6) THEN
                         k = k1 - 4
                         l = l1 - 4
                         v_d(1, ij, 15) = v_d(1, ij, 15) + ij_matrix%dd_d(1,  1, k, l) * wrepp + ij_matrix%dd( 1, k, l) * wrepp_d(1)
                         v_d(1, ij, 21) = v_d(1, ij, 21) + ij_matrix%dd_d(1,  2, k, l) * wrepp + ij_matrix%dd( 2, k, l) * wrepp_d(1)
                         v_d(1, ij, 28) = v_d(1, ij, 28) + ij_matrix%dd_d(1,  3, k, l) * wrepp + ij_matrix%dd( 3, k, l) * wrepp_d(1)
                         v_d(1, ij, 36) = v_d(1, ij, 36) + ij_matrix%dd_d(1,  4, k, l) * wrepp + ij_matrix%dd( 4, k, l) * wrepp_d(1)
                         v_d(1, ij, 45) = v_d(1, ij, 45) + ij_matrix%dd_d(1,  5, k, l) * wrepp + ij_matrix%dd( 5, k, l) * wrepp_d(1)
                         v_d(1, ij, 20) = v_d(1, ij, 20) + ij_matrix%dd_d(1,  6, k, l) * wrepp + ij_matrix%dd( 6, k, l) * wrepp_d(1)
                         v_d(1, ij, 26) = v_d(1, ij, 26) + ij_matrix%dd_d(1,  7, k, l) * wrepp + ij_matrix%dd( 7, k, l) * wrepp_d(1)
                         v_d(1, ij, 27) = v_d(1, ij, 27) + ij_matrix%dd_d(1,  8, k, l) * wrepp + ij_matrix%dd( 8, k, l) * wrepp_d(1)
                         v_d(1, ij, 33) = v_d(1, ij, 33) + ij_matrix%dd_d(1,  9, k, l) * wrepp + ij_matrix%dd( 9, k, l) * wrepp_d(1)
                         v_d(1, ij, 34) = v_d(1, ij, 34) + ij_matrix%dd_d(1, 10, k, l) * wrepp + ij_matrix%dd(10, k, l) * wrepp_d(1)
                         v_d(1, ij, 35) = v_d(1, ij, 35) + ij_matrix%dd_d(1, 11, k, l) * wrepp + ij_matrix%dd(11, k, l) * wrepp_d(1)
                         v_d(1, ij, 41) = v_d(1, ij, 41) + ij_matrix%dd_d(1, 12, k, l) * wrepp + ij_matrix%dd(12, k, l) * wrepp_d(1)
                         v_d(1, ij, 42) = v_d(1, ij, 42) + ij_matrix%dd_d(1, 13, k, l) * wrepp + ij_matrix%dd(13, k, l) * wrepp_d(1)
                         v_d(1, ij, 43) = v_d(1, ij, 43) + ij_matrix%dd_d(1, 14, k, l) * wrepp + ij_matrix%dd(14, k, l) * wrepp_d(1)
                         v_d(1, ij, 44) = v_d(1, ij, 44) + ij_matrix%dd_d(1, 15, k, l) * wrepp + ij_matrix%dd(15, k, l) * wrepp_d(1)

                         v_d(2, ij, 15) = v_d(2, ij, 15) + ij_matrix%dd_d(2,  1, k, l) * wrepp + ij_matrix%dd( 1, k, l) * wrepp_d(2)
                         v_d(2, ij, 21) = v_d(2, ij, 21) + ij_matrix%dd_d(2,  2, k, l) * wrepp + ij_matrix%dd( 2, k, l) * wrepp_d(2)
                         v_d(2, ij, 28) = v_d(2, ij, 28) + ij_matrix%dd_d(2,  3, k, l) * wrepp + ij_matrix%dd( 3, k, l) * wrepp_d(2)
                         v_d(2, ij, 36) = v_d(2, ij, 36) + ij_matrix%dd_d(2,  4, k, l) * wrepp + ij_matrix%dd( 4, k, l) * wrepp_d(2)
                         v_d(2, ij, 45) = v_d(2, ij, 45) + ij_matrix%dd_d(2,  5, k, l) * wrepp + ij_matrix%dd( 5, k, l) * wrepp_d(2)
                         v_d(2, ij, 20) = v_d(2, ij, 20) + ij_matrix%dd_d(2,  6, k, l) * wrepp + ij_matrix%dd( 6, k, l) * wrepp_d(2)
                         v_d(2, ij, 26) = v_d(2, ij, 26) + ij_matrix%dd_d(2,  7, k, l) * wrepp + ij_matrix%dd( 7, k, l) * wrepp_d(2)
                         v_d(2, ij, 27) = v_d(2, ij, 27) + ij_matrix%dd_d(2,  8, k, l) * wrepp + ij_matrix%dd( 8, k, l) * wrepp_d(2)
                         v_d(2, ij, 33) = v_d(2, ij, 33) + ij_matrix%dd_d(2,  9, k, l) * wrepp + ij_matrix%dd( 9, k, l) * wrepp_d(2)
                         v_d(2, ij, 34) = v_d(2, ij, 34) + ij_matrix%dd_d(2, 10, k, l) * wrepp + ij_matrix%dd(10, k, l) * wrepp_d(2)
                         v_d(2, ij, 35) = v_d(2, ij, 35) + ij_matrix%dd_d(2, 11, k, l) * wrepp + ij_matrix%dd(11, k, l) * wrepp_d(2)
                         v_d(2, ij, 41) = v_d(2, ij, 41) + ij_matrix%dd_d(2, 12, k, l) * wrepp + ij_matrix%dd(12, k, l) * wrepp_d(2)
                         v_d(2, ij, 42) = v_d(2, ij, 42) + ij_matrix%dd_d(2, 13, k, l) * wrepp + ij_matrix%dd(13, k, l) * wrepp_d(2)
                         v_d(2, ij, 43) = v_d(2, ij, 43) + ij_matrix%dd_d(2, 14, k, l) * wrepp + ij_matrix%dd(14, k, l) * wrepp_d(2)
                         v_d(2, ij, 44) = v_d(2, ij, 44) + ij_matrix%dd_d(2, 15, k, l) * wrepp + ij_matrix%dd(15, k, l) * wrepp_d(2)

                         v_d(3, ij, 15) = v_d(3, ij, 15) + ij_matrix%dd_d(3,  1, k, l) * wrepp + ij_matrix%dd( 1, k, l) * wrepp_d(3)
                         v_d(3, ij, 21) = v_d(3, ij, 21) + ij_matrix%dd_d(3,  2, k, l) * wrepp + ij_matrix%dd( 2, k, l) * wrepp_d(3)
                         v_d(3, ij, 28) = v_d(3, ij, 28) + ij_matrix%dd_d(3,  3, k, l) * wrepp + ij_matrix%dd( 3, k, l) * wrepp_d(3)
                         v_d(3, ij, 36) = v_d(3, ij, 36) + ij_matrix%dd_d(3,  4, k, l) * wrepp + ij_matrix%dd( 4, k, l) * wrepp_d(3)
                         v_d(3, ij, 45) = v_d(3, ij, 45) + ij_matrix%dd_d(3,  5, k, l) * wrepp + ij_matrix%dd( 5, k, l) * wrepp_d(3)
                         v_d(3, ij, 20) = v_d(3, ij, 20) + ij_matrix%dd_d(3,  6, k, l) * wrepp + ij_matrix%dd( 6, k, l) * wrepp_d(3)
                         v_d(3, ij, 26) = v_d(3, ij, 26) + ij_matrix%dd_d(3,  7, k, l) * wrepp + ij_matrix%dd( 7, k, l) * wrepp_d(3)
                         v_d(3, ij, 27) = v_d(3, ij, 27) + ij_matrix%dd_d(3,  8, k, l) * wrepp + ij_matrix%dd( 8, k, l) * wrepp_d(3)
                         v_d(3, ij, 33) = v_d(3, ij, 33) + ij_matrix%dd_d(3,  9, k, l) * wrepp + ij_matrix%dd( 9, k, l) * wrepp_d(3)
                         v_d(3, ij, 34) = v_d(3, ij, 34) + ij_matrix%dd_d(3, 10, k, l) * wrepp + ij_matrix%dd(10, k, l) * wrepp_d(3)
                         v_d(3, ij, 35) = v_d(3, ij, 35) + ij_matrix%dd_d(3, 11, k, l) * wrepp + ij_matrix%dd(11, k, l) * wrepp_d(3)
                         v_d(3, ij, 41) = v_d(3, ij, 41) + ij_matrix%dd_d(3, 12, k, l) * wrepp + ij_matrix%dd(12, k, l) * wrepp_d(3)
                         v_d(3, ij, 42) = v_d(3, ij, 42) + ij_matrix%dd_d(3, 13, k, l) * wrepp + ij_matrix%dd(13, k, l) * wrepp_d(3)
                         v_d(3, ij, 43) = v_d(3, ij, 43) + ij_matrix%dd_d(3, 14, k, l) * wrepp + ij_matrix%dd(14, k, l) * wrepp_d(3)
                         v_d(3, ij, 44) = v_d(3, ij, 44) + ij_matrix%dd_d(3, 15, k, l) * wrepp + ij_matrix%dd(15, k, l) * wrepp_d(3)
                      END IF
                   END IF
                END DO
             END DO
             DO kl = 1, limkl
                logv_d(ij, kl) = (ANY(ABS(v_d(1:3, ij, kl)) > 0.0_dp))
             END DO
          END DO
       END DO
       IF (debug_this_module) THEN
          CALL rot_2el_2c_first_debug(sepi, sepj, rijv, se_int_control, se_taper, invert, ii, kk, v_d, error)
       END IF
    END IF
  END SUBROUTINE rot_2el_2c_first

! *****************************************************************************
!> \brief Store the two-electron two-center integrals in diagonal form
!>
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE store_2el_2c_diag (limij, limkl, ww, w, ww_dx, ww_dy, ww_dz, dw, error)

    INTEGER, INTENT(IN)                      :: limij, limkl
    REAL(KIND=dp), DIMENSION(limkl, limij), &
      INTENT(IN), OPTIONAL                   :: ww
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT), OPTIONAL                :: w
    REAL(KIND=dp), DIMENSION(limkl, limij), &
      INTENT(IN), OPTIONAL                   :: ww_dx, ww_dy, ww_dz
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: dw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'store_2el_2c_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ij, kl, l
    LOGICAL                                  :: failure

    failure = .FALSE.
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
             dw(1,l) = ww_dx(kl, ij)
             dw(2,l) = ww_dy(kl, ij)
             dw(3,l) = ww_dz(kl, ij)
          END DO
       END DO
    ELSE
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE store_2el_2c_diag

END MODULE semi_empirical_int_utils
