!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Three-center integrals over Cartesian Gaussian-type functions
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author Dorothea Golze
! *****************************************************************************
MODULE ai_overlap3_debug
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap3_debug'

  INTEGER, PARAMETER            :: lmax = 5

  REAL(dp)                      :: xa, xb, xc
  REAL(dp), DIMENSION(3)        :: A, B, C
  REAL(dp), DIMENSION(3)        :: P, G 
  REAL(dp)                      :: xsi, zeta, sss

  PRIVATE
  PUBLIC :: init_os_overlap3, os_overlap3

CONTAINS

! *****************************************************************************
!> \brief   Calculation of three-center integrals over
!>          Cartesian Gaussian-type functions.
! *****************************************************************************
  SUBROUTINE init_os_overlap3(ya,yb,yc,rA,rB,rC)
    REAL(dp)                                 :: ya, yb, yc
    REAL(dp), DIMENSION(3)                   :: rA, rB, rC

    REAL(dp)                                 :: fpc, ss

    xa = ya
    xb = yb
    xc = yc
    A  = rA
    B  = rB
    C  = rC

    xsi = xa + xb
    zeta = xa*xb/xsi

    P = (xa*A + xb*B)/xsi
    G = (xsi*P + xc*C)/(xsi + xc)

    ss = (pi/xsi)**(3._dp/2._dp)*EXP(-zeta*SUM((A - B)**2))
    
    fpc = EXP(-xsi*xc/(xsi + xc)*SUM((P - C)**2))
    sss = (xsi/(xsi + xc))**(3._dp/2._dp) * ss * fpc 
  
  END SUBROUTINE init_os_overlap3

! *****************************************************************************

  RECURSIVE FUNCTION os_overlap3(an,cn,bn) RESULT(IACB)
    INTEGER, DIMENSION(3)                    :: an, cn, bn
    REAL(dp)                                 :: IACB

    INTEGER, DIMENSION(3), PARAMETER         :: i1 = (/1,0,0/), &
                                                i2 = (/0,1,0/), i3 = (/0,0,1/)

    IACB = 0._dp
    IF ( ANY(an < 0) ) RETURN
    IF ( ANY(bn < 0) ) RETURN
    IF ( ANY(cn < 0) ) RETURN

    IF ( SUM(an+cn+bn) == 0 ) THEN
      IACB = sss
      RETURN
    END IF

    IF( bn(1) > 0 )  THEN
      IACB = os_overlap3(an,cn+i1,bn-i1)+(C(1)-B(1))*os_overlap3(an,cn,bn-i1)
    ELSEIF( bn(2) > 0 )  THEN
      IACB = os_overlap3(an,cn+i2,bn-i2)+(C(2)-B(2))*os_overlap3(an,cn,bn-i2)
    ELSEIF( bn(3) > 0 )  THEN
      IACB = os_overlap3(an,cn+i3,bn-i3)+(C(3)-B(3))*os_overlap3(an,cn,bn-i3)
    ELSE
      IF( cn(1) > 0) THEN
        IACB = os_overlap3(an+i1,cn-i1,bn)+(A(1)-C(1))*os_overlap3(an,cn-i1,bn)
      ELSEIF( cn(2) > 0) THEN
        IACB = os_overlap3(an+i2,cn-i2,bn)+(A(2)-C(2))*os_overlap3(an,cn-i2,bn)
      ELSEIF( cn(3) > 0) THEN
        IACB = os_overlap3(an+i3,cn-i3,bn)+(A(3)-C(3))*os_overlap3(an,cn-i3,bn)
      ELSE
        IF( an(1) > 0) THEN
          IACB =(G(1)-A(1))*os_overlap3(an-i1,cn,bn) + &
                 0.5_dp*(an(1)-1)/(xsi+xc)*os_overlap3(an-i1-i1,cn,bn) 
        ELSEIF( an(2) > 0) THEN
          IACB =(G(2)-A(2))*os_overlap3(an-i2,cn,bn) + &
                 0.5_dp*(an(2)-1)/(xsi+xc)*os_overlap3(an-i2-i2,cn,bn) 
        ELSEIF( an(3) > 0) THEN
          IACB =(G(3)-A(3))*os_overlap3(an-i3,cn,bn) + &
                 0.5_dp*(an(3)-1)/(xsi+xc)*os_overlap3(an-i3-i3,cn,bn)
        ELSE
            STOP "I(0000)"
        ENDIF
      ENDIF
    ENDIF
         
  END FUNCTION os_overlap3


! *****************************************************************************

END MODULE ai_overlap3_debug
