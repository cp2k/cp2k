!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Coulomb integrals over Cartesian Gaussian-type functions
!>      (electron repulsion integrals, ERIs).
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author J. Hutter (07.2009)
! *****************************************************************************
MODULE ai_eri_debug
  USE f77_blas
  USE gamma,                           ONLY: fgamma_ref
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_eri_debug'

  INTEGER, PARAMETER            :: lmax = 5

  REAL(dp)                      :: xa, xb, xc, xd
  REAL(dp), DIMENSION(3)        :: A, B, C, D
  REAL(dp), DIMENSION(3)        :: P, Q, W
  REAL(dp)                      :: xsi, eta, rho, T

  REAL(dp), DIMENSION(0:4*lmax) :: fm, I0M

  PRIVATE
  PUBLIC :: init_os, os

CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Cartesian Gaussian-type functions.
!> \author  J. Hutter
!> \date    07.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_os(ya,yb,yc,yd,rA,rB,rC,rD)
    REAL(dp)                                 :: ya, yb, yc, yd
    REAL(dp), DIMENSION(3)                   :: rA, rB, rC, rD

    REAL(dp)                                 :: eab, ecd, kab, kcd

    xa = ya
    xb = yb
    xc = yc
    xd = yd
    A  = rA
    B  = rB
    C  = rC
    D  = rD

    xsi = xa + xb
    eta = xc + xd

    P = (xa*A + xb*B)/xsi
    Q = (xc*C + xd*D)/eta
    W = (xsi*P+eta*Q)/(xsi+eta)

    rho = xsi*eta/(xsi+eta)

    T = rho*SUM( (P-Q)**2 )

    fm = fgamma_ref(4*lmax,T)

    eab = -xa*xb/xsi*SUM( (A-B)**2 )
    kab = SQRT(2._dp)*pi**1.25_dp/xsi * EXP(eab)

    ecd = -xc*xd/eta*SUM( (C-D)**2 )
    kcd = SQRT(2._dp)*pi**1.25_dp/eta * EXP(ecd)

    I0M = kab*kcd/SQRT(xsi+eta) * fm

  END SUBROUTINE init_os

! *****************************************************************************

  RECURSIVE FUNCTION os(an,bn,cn,dn,mi) RESULT(IABCD)
    INTEGER, DIMENSION(3)                    :: an, bn, cn, dn
    INTEGER, OPTIONAL                        :: mi
    REAL(dp)                                 :: IABCD

    INTEGER, DIMENSION(3), PARAMETER         :: i1 = (/1,0,0/), &
                                                i2 = (/0,1,0/), i3 = (/0,0,1/)

    INTEGER                                  :: m

    m = 0
    IF(PRESENT(mi)) m=mi

    IABCD = 0._dp
    IF ( ANY(an < 0) ) RETURN
    IF ( ANY(bn < 0) ) RETURN
    IF ( ANY(cn < 0) ) RETURN
    IF ( ANY(dn < 0) ) RETURN

    IF ( SUM(an+bn+cn+dn) == 0 ) THEN
      IABCD = I0M(m)
      RETURN
    END IF

    IF ( dn(1) > 0 ) THEN
      IABCD = os(an,bn,cn+i1,dn-i1)-(D(1)-C(1))*os(an,bn,cn,dn-i1)
    ELSEIF ( dn(2) > 0 ) THEN
      IABCD = os(an,bn,cn+i2,dn-i2)-(D(2)-C(2))*os(an,bn,cn,dn-i2)
    ELSEIF ( dn(3) > 0 ) THEN
      IABCD = os(an,bn,cn+i3,dn-i3)-(D(3)-C(3))*os(an,bn,cn,dn-i3)
    ELSE
      IF ( bn(1) > 0 ) THEN
        IABCD = os(an+i1,bn-i1,cn,dn)-(B(1)-A(1))*os(an,bn-i1,cn,dn)
      ELSEIF ( bn(2) > 0 ) THEN
        IABCD = os(an+i2,bn-i2,cn,dn)-(B(2)-A(2))*os(an,bn-i2,cn,dn)
      ELSEIF ( bn(3) > 0 ) THEN
        IABCD = os(an+i3,bn-i3,cn,dn)-(B(3)-A(3))*os(an,bn-i3,cn,dn)
      ELSE
        IF ( cn(1) > 0 ) THEN
          IABCD = ((Q(1)-C(1)) + xsi/eta*(P(1)-A(1))) * os(an,bn,cn-i1,dn) + &
                  0.5_dp*an(1)/eta*os(an-i1,bn,cn-i1,dn) + &
                  0.5_dp*(cn(1)-1)/eta*os(an,bn,cn-i1-i1,dn) - &
                  xsi/eta*os(an+i1,bn,cn-i1,dn)
        ELSEIF ( cn(2) > 0 ) THEN
          IABCD = ((Q(2)-C(2)) + xsi/eta*(P(2)-A(2))) * os(an,bn,cn-i2,dn) + &
                  0.5_dp*an(2)/eta*os(an-i2,bn,cn-i2,dn) + &
                  0.5_dp*(cn(2)-1)/eta*os(an,bn,cn-i2-i2,dn) - &
                  xsi/eta*os(an+i2,bn,cn-i2,dn)
        ELSEIF ( cn(3) > 0 ) THEN
          IABCD = ((Q(3)-C(3)) + xsi/eta*(P(3)-A(3))) * os(an,bn,cn-i3,dn) + &
                  0.5_dp*an(3)/eta*os(an-i3,bn,cn-i3,dn) + &
                  0.5_dp*(cn(3)-1)/eta*os(an,bn,cn-i3-i3,dn) - &
                  xsi/eta*os(an+i3,bn,cn-i3,dn)
        ELSE
          IF ( an(1) > 0 ) THEN
            IABCD = (P(1)-A(1))*os(an-i1,bn,cn,dn,m) + &
                    (W(1)-P(1))*os(an-i1,bn,cn,dn,m+1) + &
                    0.5_dp*(an(1)-1)/xsi*os(an-i1-i1,bn,cn,dn,m) - &
                    0.5_dp*(an(1)-1)/xsi*rho/xsi*os(an-i1-i1,bn,cn,dn,m+1)
          ELSEIF ( an(2) > 0 ) THEN
            IABCD = (P(2)-A(2))*os(an-i2,bn,cn,dn,m) + &
                    (W(2)-P(2))*os(an-i2,bn,cn,dn,m+1) + &
                    0.5_dp*(an(2)-1)/xsi*os(an-i2-i2,bn,cn,dn,m) - &
                    0.5_dp*(an(2)-1)/xsi*rho/xsi*os(an-i2-i2,bn,cn,dn,m+1)
          ELSEIF ( an(3) > 0 ) THEN
            IABCD = (P(3)-A(3))*os(an-i3,bn,cn,dn,m) + &
                    (W(3)-P(3))*os(an-i3,bn,cn,dn,m+1) + &
                    0.5_dp*(an(3)-1)/xsi*os(an-i3-i3,bn,cn,dn,m) - &
                    0.5_dp*(an(3)-1)/xsi*rho/xsi*os(an-i3-i3,bn,cn,dn,m+1)
          ELSE
            STOP "I(0000)"
          END IF
        END IF
      END IF
    END IF

  END FUNCTION os


! *****************************************************************************

END MODULE ai_eri_debug
