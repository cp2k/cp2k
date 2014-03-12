!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of Coulomb integrals over Correlated Cartesian Gaussian
!>        Geminals (electron repulsion integrals, ERIs).
!>        This is a slow reference implementation for debugging perposes
!> \par Literature
!>      S.L. Saito and Y. Suzuki, J. Chem. Phys. 114 1115 (2001)
!> \par History
!>      none
!> \author Juerg Hutter (07.2009)
! *****************************************************************************
MODULE ai_geminals_debug
  USE ai_geminals_utils,               ONLY: acoff,&
                                             matv
  USE f77_blas
  USE gamma,                           ONLY: fgamma_ref
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE orbital_pointers,                ONLY: coset
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_geminals_debug'

  PRIVATE

  PUBLIC  init_grs, grs, transform1, transform2

  INTEGER, PARAMETER            :: lmax = 10

  REAL(dp), DIMENSION(2,2)      :: xa, xb, xg, xgm
  REAL(dp), DIMENSION(6)        :: A, B, Q, AB
  REAL(dp)                      :: detg, rho, T
  REAL(dp)                      :: kab

  REAL(dp), DIMENSION(0:2*lmax) :: gm, I0M

  CONTAINS

! *****************************************************************************
!> \brief   Calculation of the primitive two-center Coulomb integrals over
!>          Correlated Cartesian Gaussian Geminals
!> \author  Juerg Hutter
!> \date    06.2009
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_grs(ya,yb,rA,rB)
    REAL(dp), DIMENSION(2, 2)                :: ya, yb
    REAL(dp), DIMENSION(6)                   :: rA, rB

    REAL(dp)                                 :: eab
    REAL(dp), DIMENSION(2, 2)                :: ee
    REAL(dp), DIMENSION(6)                   :: C, D

    xa = ya
    xb = yb
    A  = rA
    B  = rB
    AB = A-B

    xg = xa + xb
    detg = xg(1,1)*xg(2,2)-xg(1,2)*xg(2,1)
    xgm(1,1) = xg(2,2)
    xgm(1,2) = -xg(2,1)
    xgm(2,1) = -xg(1,2)
    xgm(2,2) = xg(1,1)
    xgm = xgm/detg

    C = matv(xa,A)+matv(xb,B)
    Q = matv(xgm,C)

    rho = detg/xg(1,1)
    T = rho*(Q(4)*Q(4) + Q(5)*Q(5) + Q(6)*Q(6))

    gm = fgamma_ref(2*lmax,T) * 2._dp*pi/rho

    C = A-B
    ee = MATMUL(xa,MATMUL(xgm,xb))
    D = matv(ee,C)
    eab = SUM(C*D)
    kab = EXP(-eab)

    I0M = 4._dp*kab*(pi/xg(1,1))**1.5_dp * gm

  END SUBROUTINE init_grs

  RECURSIVE FUNCTION grs(an,bn,mi) RESULT(IAB)
    INTEGER, DIMENSION(6)                    :: an, bn
    INTEGER, OPTIONAL                        :: mi
    REAL(dp)                                 :: IAB

    INTEGER, DIMENSION(6), PARAMETER :: i1 = (/1,0,0,0,0,0/), &
      i2 = (/0,1,0,0,0,0/), i3 = (/0,0,1,0,0,0/), i4 = (/0,0,0,1,0,0/), &
      i5 = (/0,0,0,0,1,0/), i6 = (/0,0,0,0,0,1/)

    INTEGER                                  :: m
    REAL(dp)                                 :: Im0, Im1, Ir0, Ir1, Is0, Is1, &
                                                v(2), w(2)

    m = 0
    IF(PRESENT(mi)) m=mi

    IAB = 0._dp
    IF ( ANY(an < 0) ) RETURN
    IF ( ANY(bn < 0) ) RETURN

    IF ( SUM(an+bn) == 0 ) THEN
      IAB = I0M(m)
      RETURN
    END IF

    IF(bn(1) > 0) THEN
      IAB = grs(an+i1,bn-i1)+AB(1)*grs(an,bn-i1)
    ELSEIF(bn(2) > 0) THEN
      IAB = grs(an+i2,bn-i2)+AB(2)*grs(an,bn-i2)
    ELSEIF(bn(3) > 0) THEN
      IAB = grs(an+i3,bn-i3)+AB(3)*grs(an,bn-i3)
    ELSEIF(bn(4) > 0) THEN
      IAB = grs(an+i4,bn-i4)+AB(4)*grs(an,bn-i4)
    ELSEIF(bn(5) > 0) THEN
      IAB = grs(an+i5,bn-i5)+AB(5)*grs(an,bn-i5)
    ELSEIF(bn(6) > 0) THEN
      IAB = grs(an+i6,bn-i6)+AB(6)*grs(an,bn-i6)
    ELSE
      IF(an(1) > 0) THEN
        Im0 = grs(an-i1,bn,m)
        Im1 = grs(an-i1,bn,m+1)
        Ir0 = grs(an-i1-i1,bn,m)
        Ir1 = grs(an-i1-i1,bn,m+1)
        Is0 = grs(an-i1-i4,bn,m)
        Is1 = grs(an-i1-i4,bn,m+1)
        v(1) = AB(1)*(Im1-Im0)
        v(2) = AB(4)*(Im1-Im0)
        v = MATMUL(xb,v)
        v(1) = v(1)+0.5_dp*(an(1)-1)*(Ir0-Ir1)
        v(2) = v(2)+0.5_dp*an(4)*(Is0-Is1)
        w = MATMUL(xgm,v)
        w(1) = w(1) + 0.5_dp/xg(1,1)*(an(1)-1)*Ir1
        w(1) = w(1) - 1.0_dp/xg(1,1)*(xb(1,1)*AB(1)+xb(1,2)*AB(4)-xg(1,2)*A(4))*Im1
        IAB = w(1)
      ELSEIF(an(2) > 0) THEN
        Im0 = grs(an-i2,bn,m)
        Im1 = grs(an-i2,bn,m+1)
        Ir0 = grs(an-i2-i2,bn,m)
        Ir1 = grs(an-i2-i2,bn,m+1)
        Is0 = grs(an-i2-i5,bn,m)
        Is1 = grs(an-i2-i5,bn,m+1)
        v(1) = AB(2)*(Im1-Im0)
        v(2) = AB(5)*(Im1-Im0)
        v = MATMUL(xb,v)
        v(1) = v(1)+0.5_dp*(an(2)-1)*(Ir0-Ir1)
        v(2) = v(2)+0.5_dp*an(5)*(Is0-Is1)
        w = MATMUL(xgm,v)
        w(1) = w(1) + 0.5_dp/xg(1,1)*(an(2)-1)*Ir1
        w(1) = w(1) - 1.0_dp/xg(1,1)*(xb(1,1)*AB(2)+xb(1,2)*AB(5)-xg(1,2)*A(5))*Im1
        IAB = w(1)
      ELSEIF(an(3) > 0) THEN
        Im0 = grs(an-i3,bn,m)
        Im1 = grs(an-i3,bn,m+1)
        Ir0 = grs(an-i3-i3,bn,m)
        Ir1 = grs(an-i3-i3,bn,m+1)
        Is0 = grs(an-i3-i6,bn,m)
        Is1 = grs(an-i3-i6,bn,m+1)
        v(1) = AB(3)*(Im1-Im0)
        v(2) = AB(6)*(Im1-Im0)
        v = MATMUL(xb,v)
        v(1) = v(1)+0.5_dp*(an(3)-1)*(Ir0-Ir1)
        v(2) = v(2)+0.5_dp*an(6)*(Is0-Is1)
        w = MATMUL(xgm,v)
        w(1) = w(1) + 0.5_dp/xg(1,1)*(an(3)-1)*Ir1
        w(1) = w(1) - 1.0_dp/xg(1,1)*(xb(1,1)*AB(3)+xb(1,2)*AB(6)-xg(1,2)*A(6))*Im1
        IAB = w(1)
      ELSEIF(an(4) > 0) THEN
        Im0 = grs(an-i4,bn,m)
        Im1 = grs(an-i4,bn,m+1)
        Ir0 = grs(an-i4-i1,bn,m)
        Ir1 = grs(an-i4-i1,bn,m+1)
        Is0 = grs(an-i4-i4,bn,m)
        Is1 = grs(an-i4-i4,bn,m+1)
        v(1) = AB(1)*(Im1-Im0)
        v(2) = AB(4)*(Im1-Im0)
        v = MATMUL(xb,v)
        v(1) = v(1)+0.5_dp*(an(1))*(Ir0-Ir1)
        v(2) = v(2)+0.5_dp*(an(4)-1)*(Is0-Is1)
        w = MATMUL(xgm,v)
        w(2) = w(2) - A(4)*Im1
        IAB = w(2)
      ELSEIF(an(5) > 0) THEN
        Im0 = grs(an-i5,bn,m)
        Im1 = grs(an-i5,bn,m+1)
        Ir0 = grs(an-i5-i2,bn,m)
        Ir1 = grs(an-i5-i2,bn,m+1)
        Is0 = grs(an-i5-i5,bn,m)
        Is1 = grs(an-i5-i5,bn,m+1)
        v(1) = AB(2)*(Im1-Im0)
        v(2) = AB(5)*(Im1-Im0)
        v = MATMUL(xb,v)
        v(1) = v(1)+0.5_dp*(an(2))*(Ir0-Ir1)
        v(2) = v(2)+0.5_dp*(an(5)-1)*(Is0-Is1)
        w = MATMUL(xgm,v)
        w(2) = w(2) - A(5)*Im1
        IAB = w(2)
      ELSEIF(an(6) > 0) THEN
        Im0 = grs(an-i6,bn,m)
        Im1 = grs(an-i6,bn,m+1)
        Ir0 = grs(an-i6-i3,bn,m)
        Ir1 = grs(an-i6-i3,bn,m+1)
        Is0 = grs(an-i6-i6,bn,m)
        Is1 = grs(an-i6-i6,bn,m+1)
        v(1) = AB(3)*(Im1-Im0)
        v(2) = AB(6)*(Im1-Im0)
        v = MATMUL(xb,v)
        v(1) = v(1)+0.5_dp*(an(3))*(Ir0-Ir1)
        v(2) = v(2)+0.5_dp*(an(6)-1)*(Is0-Is1)
        w = MATMUL(xgm,v)
        w(2) = w(2) - A(6)*Im1
        IAB = w(2)
      ELSE
        STOP "(an|bn) out of range"
      END IF
    END IF

  END FUNCTION grs

! *****************************************************************************

  SUBROUTINE transform1(in,iout,ac,la,lb)
    REAL(dp), DIMENSION(0:, 0:, 0:, 0:)      :: in, iout
    REAL(dp), DIMENSION(0:, 0:, 0:)          :: ac
    INTEGER, INTENT(IN)                      :: la, lb

    INTEGER                                  :: i1, i2, i3, iu, iv, iw, j1, &
                                                j2, j3, k1, k2, l1, l2, m1, &
                                                m2, nmax

     nmax=MAX(la,lb)
     CALL acoff(ac,nmax)

     iout=0._dp
     DO l1=0,la
       DO l2=0,lb

         DO i1=0,l1
           DO i2=0,l1-i1
             i3=l1-i1-i2
             k1=coset(i1,i2,i3)
             DO j1=0,l2
               DO j2=0,l2-j1
                 j3=l2-j1-j2
                 k2=coset(j1,j2,j3)
                 !transform integrals (ij|.)
                 DO iu=0,i1+j1
                   DO iv=0,i2+j2
                     DO iw=0,i3+j3
                       m1=coset(i1+j1-iu,i2+j2-iv,i3+j3-iw)
                       m2=coset(iu,iv,iw)
                       iout(k1,k2,:,:) = iout(k1,k2,:,:) + &
                            ac(i1,j1,iu)*ac(i2,j2,iv)*ac(i3,j3,iw)*in(m1,m2,:,:)
                     END DO
                   END DO
                 END DO
               END DO
             END DO
           END DO
         END DO

       END DO
     END DO

  END SUBROUTINE transform1

! *****************************************************************************

  SUBROUTINE transform2(in,iout,ac,la,lb)
    REAL(dp), DIMENSION(0:, 0:, 0:, 0:)      :: in, iout
    REAL(dp), DIMENSION(0:, 0:, 0:)          :: ac
    INTEGER, INTENT(IN)                      :: la, lb

    INTEGER                                  :: i1, i2, i3, iu, iv, iw, j1, &
                                                j2, j3, k1, k2, l1, l2, m1, &
                                                m2, nmax

     nmax=MAX(la,lb)
     CALL acoff(ac,nmax)

     iout=0._dp
     DO l1=0,la
       DO l2=0,lb

         DO i1=0,l1
           DO i2=0,l1-i1
             i3=l1-i1-i2
             k1=coset(i1,i2,i3)
             DO j1=0,l2
               DO j2=0,l2-j1
                 j3=l2-j1-j2
                 k2=coset(j1,j2,j3)
                 !transform integrals (.|ij)
                 DO iu=0,i1+j1
                   DO iv=0,i2+j2
                     DO iw=0,i3+j3
                       m1=coset(i1+j1-iu,i2+j2-iv,i3+j3-iw)
                       m2=coset(iu,iv,iw)
                       iout(:,:,k1,k2) = iout(:,:,k1,k2) + &
                            ac(i1,j1,iu)*ac(i2,j2,iv)*ac(i3,j3,iw)*in(:,:,m1,m2)
                     END DO
                   END DO
                 END DO
               END DO
             END DO
           END DO
         END DO

       END DO
     END DO

  END SUBROUTINE transform2

! *****************************************************************************

END MODULE ai_geminals_debug
