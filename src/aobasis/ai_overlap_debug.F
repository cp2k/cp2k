!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Two-center overlap integrals over Cartesian Gaussian-type functions
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author Dorothea Golze
! *****************************************************************************
MODULE ai_overlap_debug
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap_debug'

  INTEGER, PARAMETER            :: lmax = 5

  REAL(dp)                      :: xa, xb 
  REAL(dp), DIMENSION(3)        :: A, B 
  REAL(dp), DIMENSION(3)        :: P 
  REAL(dp)                      :: xsi, zeta, ss

  PRIVATE
  PUBLIC :: init_os_overlap2, os_overlap2

CONTAINS

! *****************************************************************************
!> \brief   Calculation of overlap integrals over
!>          Cartesian Gaussian-type functions.
! *****************************************************************************
  SUBROUTINE init_os_overlap2(ya,yb,rA,rB)
    REAL(dp)                                 :: ya, yb
    REAL(dp), DIMENSION(3)                   :: rA, rB

    xa = ya
    xb = yb
    A  = rA
    B  = rB

    xsi = xa + xb
    zeta = xa*xb/xsi

    P = (xa*A + xb*B)/xsi

    ss = (pi/xsi)**(3._dp/2._dp)*EXP(-zeta*SUM((A - B)**2))
  
  END SUBROUTINE init_os_overlap2

! *****************************************************************************

  RECURSIVE FUNCTION os_overlap2(an,bn) RESULT(IAB)
    INTEGER, DIMENSION(3)                    :: an, bn
    REAL(dp)                                 :: IAB

    INTEGER, DIMENSION(3), PARAMETER         :: i1 = (/1,0,0/), &
                                                i2 = (/0,1,0/), i3 = (/0,0,1/)

    IAB = 0._dp
    IF ( ANY(an < 0) ) RETURN
    IF ( ANY(bn < 0) ) RETURN

    IF ( SUM(an+bn) == 0 ) THEN
      IAB = ss
      RETURN
    END IF

    IF( bn(1) > 0 )  THEN
      IAB = os_overlap2(an+i1,bn-i1)+(A(1)-B(1))*os_overlap2(an,bn-i1)
    ELSEIF( bn(2) > 0 )  THEN
      IAB = os_overlap2(an+i2,bn-i2)+(A(2)-B(2))*os_overlap2(an,bn-i2)
    ELSEIF( bn(3) > 0 )  THEN
      IAB = os_overlap2(an+i3,bn-i3)+(A(3)-B(3))*os_overlap2(an,bn-i3)
    ELSE
     IF( an(1) > 0) THEN
       IAB =(P(1)-A(1))*os_overlap2(an-i1,bn) + &
              0.5_dp*(an(1)-1)/xsi*os_overlap2(an-i1-i1,bn) 
     ELSEIF( an(2) > 0) THEN
       IAB =(P(2)-A(2))*os_overlap2(an-i2,bn) + &
              0.5_dp*(an(2)-1)/xsi*os_overlap2(an-i2-i2,bn) 
     ELSEIF( an(3) > 0) THEN
       IAB =(P(3)-A(3))*os_overlap2(an-i3,bn) + &
              0.5_dp*(an(3)-1)/xsi*os_overlap2(an-i3-i3,bn)
     ELSE
         STOP "I(0000)"
     ENDIF
    ENDIF
         
  END FUNCTION os_overlap2


! *****************************************************************************

END MODULE ai_overlap_debug
