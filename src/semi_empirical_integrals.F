!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Integrals for semi-empiric methods
!> \par History
!>      JGH (27.10.2004) : separate routine for nuclear attraction integrals
!>      JGH (13.03.2006) : tapering function
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!> \author JGH (11.10.2004)
! *****************************************************************************
MODULE semi_empirical_integrals

  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             evolt
  USE semi_empirical_int_arrays,       ONLY: pxy,&
                                             rij_threshold
  USE semi_empirical_int_utils,        ONLY: taper
  USE semi_empirical_types,            ONLY: semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_integrals'
  PUBLIC :: rotint, rotnuc, drotint, drotnuc, nucint, terep

CONTAINS

! *****************************************************************************
!> \brief Computes the two particle interactions in the lab frame
!> \param on input  sepi   = Atomic parameters of first atom
!>                  sepj   = Atomic parameters of second atom
!>                  rijv   = Coordinate vector i -> j
!>            shortrange   = compute only the short range part of the integrals
!> 
!>    on output     w      = Array of two-electron repulsion integrals.
!> 
!> \note routine adapted from mopac7 (rotate)
!>       written by Ernest R. Davidson, Indiana University.
! *****************************************************************************
  SUBROUTINE rotint (sepi,sepj,rijv,w,shortrange,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: w
    LOGICAL                                  :: shortrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: si, sj
    REAL(dp) :: a, rij, xx11, xx21, xx22, xx31, xx32, xx33, xy11, xy21, xy22, &
      xy31, xy32, xz11, xz21, xz22, xz31, xz32, xz33, yy11, yy21, yy22, &
      yyzz11, yyzz21, yyzz22, yz11, yz21, yz22, yz31, yz32, zz11, zz21, zz22, &
      zz31, zz32, zz33
    REAL(dp), DIMENSION(22)                  :: ri
    REAL(dp), DIMENSION(3)                   :: x, y, z

    x(:)=-rijv(:)
    rij=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    w(:)=0._dp
    IF (rij > rij_threshold) THEN
       ! The repulsion integrals over molecular frame (w) are stored in the
       ! order in which they will later be used.  ie.  (i,j/k,l) where
       ! j.le.i  and  l.le.k     and l varies most rapidly and i least
       ! rapidly.  (anti-normal computer storage)
       rij = SQRT(rij)

       ! Compute Integrals in Diatomic Frame
       CALL terep(sepi,sepj,rij,ri,shortrange,error=error)
 
       a=1._dp/rij
       x(1) = x(1)*a
       x(2) = x(2)*a
       x(3) = x(3)*a
       IF (ABS(x(3)) > 0.99999999_dp) THEN
          x(3) = SIGN(1._dp,x(3))
          y(1) = 0._dp
          y(2) = 1._dp
          y(3) = 0._dp
          z(1) = 1._dp
          z(2) = 0._dp
          z(3) = 0._dp
       ELSE
          z(3)=SQRT(1._dp-x(3)*x(3))
          a=1._dp/z(3)
          y(1)=-a*x(2)*SIGN(1._dp,x(1))
          y(2)=ABS(a*x(1))
          y(3)=0._dp
          z(1)=-a*x(1)*x(3)
          z(2)=-a*x(2)*x(3)
       ENDIF
       si = (sepi%natorb > 1)
       sJ = (sepj%natorb > 1)
       IF ( si .OR. sj ) THEN
          xx11 = x(1)*x(1)
          xx21 = x(2)*x(1)
          xx22 = x(2)*x(2)
          xx31 = x(3)*x(1)
          xx32 = x(3)*x(2)
          xx33 = x(3)*x(3)
          yy11 = y(1)*y(1)
          YY21 = Y(2)*Y(1)
          yy22 = y(2)*y(2)
          zz11 = z(1)*z(1)
          zz21 = z(2)*z(1)
          zz22 = z(2)*z(2)
          zz31 = z(3)*z(1)
          zz32 = z(3)*z(2)
          zz33 = z(3)*z(3)
          yyzz11 = yy11+zz11
          yyzz21 = yy21+zz21
          yyzz22 = yy22+zz22
          xy11 = 2._dp*x(1)*y(1)
          xy21 =      x(1)*y(2)+x(2)*y(1)
          xy22 = 2._dp*x(2)*y(2)
          xy31 =      x(3)*y(1)
          xy32 =      x(3)*y(2)
          xz11 = 2._dp*x(1)*z(1)
          xz21 =      x(1)*z(2)+x(2)*z(1)
          xz22 = 2._dp*x(2)*z(2)
          xz31 =      x(1)*z(3)+x(3)*z(1)
          xz32 =      x(2)*z(3)+x(3)*z(2)
          xz33 = 2._dp*x(3)*z(3)
          yz11 = 2._dp*y(1)*z(1)
          yz21 =      y(1)*z(2)+y(2)*z(1)
          yz22 = 2._dp*y(2)*z(2)
          yz31 =      y(1)*z(3)
          yz32 =      y(2)*z(3)
       ENDIF
 
!   (s s/s s)
       w(1)=ri(1)
       IF (sj) THEN
!   (s s/px s)
          w(2)=ri(5)*x(1)
!   (s s/px px)
          w(3)=ri(11)*xx11+ri(12)*yyzz11
!   (s s/py s)
          w(4)=ri(5)*x(2)
!   (s s/py px)
          w(5)=ri(11)*xx21+ri(12)*yyzz21
!   (s s/py py)
          w(6)=ri(11)*xx22+ri(12)*yyzz22
!   (s s/pz s)
          w(7)=ri(5)*x(3)
!   (s s/pz px)
          w(8)=ri(11)*xx31+ri(12)*zz31
!   (s s/pz py)
          w(9)=ri(11)*xx32+ri(12)*zz32
!   (s s/pz pz)
          w(10)=ri(11)*xx33+ri(12)*zz33
       END IF
 
       IF (si) THEN
!   (px s/s s)
          w(11)=ri(2)*x(1)
          IF (sj) THEN
!   (px s/px s)
             w(12)=ri(6)*xx11+ri(7)*yyzz11
!   (px s/px px)
             w(13)=x(1)*(ri(13)*xx11+ri(14)*yyzz11) &
                  +ri(15)*(y(1)*xy11+z(1)*xz11)
!   (px s/py s)
             w(14)=ri(6)*xx21+ri(7)*yyzz21
!   (px s/py px)
             w(15)=x(1)*(ri(13)*xx21+ri(14)*yyzz21) &
                  +ri(15)*(y(1)*xy21+z(1)*xz21)
!   (px s/py py)
             w(16)=x(1)*(ri(13)*xx22+ri(14)*yyzz22) &
                  +ri(15)*(y(1)*xy22+z(1)*xz22)
!   (px s/pz s)
             w(17)=ri(6)*xx31+ri(7)*zz31
!   (px s/pz px)
             w(18)=x(1)*(ri(13)*xx31+ri(14)*zz31) &
                  +ri(15)*(y(1)*xy31+z(1)*xz31)
!   (px s/pz py)
             w(19)=x(1)*(ri(13)*xx32+ri(14)*zz32) &
                  +ri(15)*(y(1)*xy32+z(1)*xz32)
!   (px s/pz pz)
             w(20)=x(1)*(ri(13)*xx33+ri(14)*zz33) &
                  +ri(15)*(          z(1)*xz33)
!   (px px/s s)
             w(21)=ri(3)*xx11+ri(4)*yyzz11
!   (px px/px s)
             w(22)=x(1)*(ri(8)*xx11+ri(9)*yyzz11) &
                  +ri(10)*(y(1)*xy11+z(1)*xz11)
!   (px px/px px)
             w(23) = &
                  (ri(16)*xx11+ri(17)*yyzz11)*xx11+ri(18)*xx11*yyzz11 &
                  +ri(19)*(yy11*yy11+zz11*zz11) &
                  +ri(20)*(xy11*xy11+xz11*xz11) &
                  +ri(21)*(yy11*zz11+zz11*yy11) &
                  +ri(22)*yz11*yz11
!   (px px/py s)
             w(24)=x(2)*(ri(8)*xx11+ri(9)*yyzz11) &
                  +ri(10)*(y(2)*xy11+z(2)*xz11)
!   (px px/py px)
             w(25) = &
                  (ri(16)*xx11+ri(17)*yyzz11)*xx21+ri(18)*xx11*yyzz21 &
                  +ri(19)*(yy11*yy21+zz11*zz21) &
                  +ri(20)*(xy11*xy21+xz11*xz21) &
                  +ri(21)*(yy11*zz21+zz11*yy21) &
                  +ri(22)*yz11*yz21
!   (px px/py py)
             w(26) = &
                  (ri(16)*xx11+ri(17)*yyzz11)*xx22+ri(18)*xx11*yyzz22 &
                  +ri(19)*(yy11*yy22+zz11*zz22) &
                  +ri(20)*(xy11*xy22+xz11*xz22) &
                  +ri(21)*(yy11*zz22+zz11*yy22) &
                  +ri(22)*yz11*yz22
!   (px px/pz s)
             w(27)=x(3)*(ri(8)*xx11+ri(9)*yyzz11) &
                  +ri(10)*(         +z(3)*xz11)
!   (px px/pz px)
             w(28) = &
                  (ri(16)*xx11+ri(17)*yyzz11)*xx31 &
                  +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz31 &
                  +ri(20)*(xy11*xy31+xz11*xz31) &
                  +ri(22)*yz11*yz31
!   (px px/pz py)
             w(29) = &
                  (ri(16)*xx11+ri(17)*yyzz11)*xx32 &
                  +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz32 &
                  +ri(20)*(xy11*xy32+xz11*xz32) &
                  +ri(22)*yz11*yz32
!   (px px/pz pz)
             w(30) = &
                  (ri(16)*xx11+ri(17)*yyzz11)*xx33 &
                  +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz33 &
                  +ri(20)*xz11*xz33
!   (py s/s s)
             w(31)=ri(2)*x(2)
!   (py s/px s)
             w(32)=ri(6)*xx21+ri(7)*yyzz21
!   (py s/px px)
             w(33)=x(2)*(ri(13)*xx11+ri(14)*yyzz11) &
                  +ri(15)*(y(2)*xy11+z(2)*xz11)
!   (py s/py s)
             w(34)=ri(6)*xx22+ri(7)*yyzz22
!   (py s/py px)
             w(35)=x(2)*(ri(13)*xx21+ri(14)*yyzz21) &
                  +ri(15)*(y(2)*xy21+z(2)*xz21)
!   (py s/py py)
             w(36)=x(2)*(ri(13)*xx22+ri(14)*yyzz22) &
                  +ri(15)*(y(2)*xy22+z(2)*xz22)
!   (py s/pz s)
             w(37)=ri(6)*xx32+ri(7)*zz32
!   (py s/pz px)
             w(38)=x(2)*(ri(13)*xx31+ri(14)*zz31) &
                  +ri(15)*(y(2)*xy31+z(2)*xz31)
!   (py s/pz py)
             w(39)=x(2)*(ri(13)*xx32+ri(14)*zz32) &
                  +ri(15)*(y(2)*xy32+z(2)*xz32)
!   (py s/pz pz)
             w(40)=x(2)*(ri(13)*xx33+ri(14)*zz33) &
                  +ri(15)*(         +z(2)*xz33)
!   (py px/s s)
             w(41)=ri(3)*xx21+ri(4)*yyzz21
!   (py px/px s)
             w(42)=x(1)*(ri(8)*xx21+ri(9)*yyzz21) &
                  +ri(10)*(y(1)*xy21+z(1)*xz21)
!   (py px/px px)
             w(43) = &
                  (ri(16)*xx21+ri(17)*yyzz21)*xx11+ri(18)*xx21*yyzz11 &
                  +ri(19)*(yy21*yy11+zz21*zz11) &
                  +ri(20)*(xy21*xy11+xz21*xz11) &
                  +ri(21)*(yy21*zz11+zz21*yy11) &
                  +ri(22)*yz21*yz11
!   (py px/py s)
             w(44)=x(2)*(ri(8)*xx21+ri(9)*yyzz21) &
                  +ri(10)*(y(2)*xy21+z(2)*xz21)
!   (py px/py px)
             w(45) = &
                  (ri(16)*xx21+ri(17)*yyzz21)*xx21+ri(18)*xx21*yyzz21 &
                  +ri(19)*(yy21*yy21+zz21*zz21) &
                  +ri(20)*(xy21*xy21+xz21*xz21) &
                  +ri(21)*(yy21*zz21+zz21*yy21) &
                  +ri(22)*yz21*yz21
!   (py px/py py)
             w(46) = &
                  (ri(16)*xx21+ri(17)*yyzz21)*xx22+ri(18)*xx21*yyzz22 &
                  +ri(19)*(yy21*yy22+zz21*zz22) &
                  +ri(20)*(xy21*xy22+xz21*xz22) &
                  +ri(21)*(yy21*zz22+zz21*yy22) &
                  +ri(22)*yz21*yz22
!   (py px/pz s)
             w(47)=x(3)*(ri(8)*xx21+ri(9)*yyzz21) &
                  +ri(10)*(         +z(3)*xz21)
!   (py px/pz px)
             w(48) = &
                  (ri(16)*xx21+ri(17)*yyzz21)*xx31 &
                  +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz31 &
                  +ri(20)*(xy21*xy31+xz21*xz31) &
                  +ri(22)*yz21*yz31
!   (py px/pz py)
             w(49) = &
                  (ri(16)*xx21+ri(17)*yyzz21)*xx32 &
                  +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz32 &
                  +ri(20)*(xy21*xy32+xz21*xz32) &
                  +ri(22)*yz21*yz32
!   (py px/pz pz)
             w(50) = &
                  (ri(16)*xx21+ri(17)*yyzz21)*xx33 &
                  +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz33 &
                  +ri(20)*xz21*xz33
!   (py py/s s)
             w(51)=ri(3)*xx22+ri(4)*yyzz22
!   (py py/px s)
             w(52)=x(1)*(ri(8)*xx22+ri(9)*yyzz22) &
                  +ri(10)*(y(1)*xy22+z(1)*xz22)
!   (py py/px px)
             w(53) = &
                  (ri(16)*xx22+ri(17)*yyzz22)*xx11+ri(18)*xx22*yyzz11 &
                  +ri(19)*(yy22*yy11+zz22*zz11) &
                  +ri(20)*(xy22*xy11+xz22*xz11) &
                  +ri(21)*(yy22*zz11+zz22*yy11) &
                  +ri(22)*yz22*yz11
!   (py py/py s)
             w(54)=x(2)*(ri(8)*xx22+ri(9)*yyzz22) &
                  +ri(10)*(y(2)*xy22+z(2)*xz22)
!   (py py/py px)
             w(55) = &
                  (ri(16)*xx22+ri(17)*yyzz22)*xx21+ri(18)*xx22*yyzz21 &
                  +ri(19)*(yy22*yy21+zz22*zz21) &
                  +ri(20)*(xy22*xy21+xz22*xz21) &
                  +ri(21)*(yy22*zz21+zz22*yy21) &
                  +ri(22)*yz22*yz21
!   (py py/py py)
             w(56) = &
                  (ri(16)*xx22+ri(17)*yyzz22)*xx22+ri(18)*xx22*yyzz22 &
                  +ri(19)*(yy22*yy22+zz22*zz22) &
                  +ri(20)*(xy22*xy22+xz22*xz22) &
                  +ri(21)*(yy22*zz22+zz22*yy22) &
                  +ri(22)*yz22*yz22
!   (py py/pz s)
             w(57)=x(3)*(ri(8)*xx22+ri(9)*yyzz22) &
                  +ri(10)*(         +z(3)*xz22)
!   (py py/pz px)
             w(58) = &
                  (ri(16)*xx22+ri(17)*yyzz22)*xx31 &
                  +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz31 &
                  +ri(20)*(xy22*xy31+xz22*xz31) &
                  +ri(22)*yz22*yz31
!   (py py/pz py)
             w(59) = &
                  (ri(16)*xx22+ri(17)*yyzz22)*xx32 &
                  +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz32 &
                  +ri(20)*(xy22*xy32+xz22*xz32) &
                  +ri(22)*yz22*yz32
!   (py py/pz pz)
             w(60) = &
                  (ri(16)*xx22+ri(17)*yyzz22)*xx33 &
                  +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz33 &
                  +ri(20)*xz22*xz33
!   (pz s/ss)
             w(61)=ri(2)*x(3)
!   (pz s/px s)
             w(62)=ri(6)*xx31+ri(7)*zz31
!   (pz s/px px)
             w(63)=x(3)*(ri(13)*xx11+ri(14)*yyzz11) &
                  +ri(15)*(         +z(3)*xz11)
!   (pz s/py s)
             w(64)=ri(6)*xx32+ri(7)*zz32
!   (pz s/py px)
             w(65)=x(3)*(ri(13)*xx21+ri(14)*yyzz21) &
                  +ri(15)*(         +z(3)*xz21)
!   (pz s/py py)
             w(66)=x(3)*(ri(13)*xx22+ri(14)*yyzz22) &
                  +ri(15)*(         +z(3)*xz22)
!   (pz s/pz s)
             w(67)=ri(6)*xx33+ri(7)*zz33
!   (pz s/pz px)
             w(68)=x(3)*(ri(13)*xx31+ri(14)*zz31) &
                  +ri(15)*(         +z(3)*xz31)
!   (pz s/pz py)
             w(69)=x(3)*(ri(13)*xx32+ri(14)*zz32) &
                  +ri(15)*(         +z(3)*xz32)
!   (pz s/pz pz)
             w(70)=x(3)*(ri(13)*xx33+ri(14)*zz33) &
                  +ri(15)*(         +z(3)*xz33)
!   (pz px/s s)
             w(71)=ri(3)*xx31+ri(4)*zz31
!   (pz px/px s)
             w(72)=x(1)*(ri(8)*xx31+ri(9)*zz31) &
                  +ri(10)*(y(1)*xy31+z(1)*xz31)
!   (pz px/px px)
             w(73) = &
                  (ri(16)*xx31+ri(17)*zz31)*xx11+ri(18)*xx31*yyzz11 &
                  +ri(19)*zz31*zz11 &
                  +ri(20)*(xy31*xy11+xz31*xz11) &
                  +ri(21)*zz31*yy11 &
                  +ri(22)*yz31*yz11
!   (pz px/py s)
             w(74)=x(2)*(ri(8)*xx31+ri(9)*zz31) &
                  +ri(10)*(y(2)*xy31+z(2)*xz31)
!   (pz px/py px)
             w(75) = &
                  (ri(16)*xx31+ri(17)*zz31)*xx21+ri(18)*xx31*yyzz21 &
                  +ri(19)*zz31*zz21 &
                  +ri(20)*(xy31*xy21+xz31*xz21) &
                  +ri(21)*zz31*yy21 &
                  +ri(22)*yz31*yz21
!   (pz px/py py)
             w(76) = &
                  (ri(16)*xx31+ri(17)*zz31)*xx22+ri(18)*xx31*yyzz22 &
                  +ri(19)*zz31*zz22 &
                  +ri(20)*(xy31*xy22+xz31*xz22) &
                  +ri(21)*zz31*yy22 &
                  +ri(22)*yz31*yz22
!   (pz px/pz s)
             w(77)=x(3)*(ri(8)*xx31+ri(9)*zz31) &
                  +ri(10)*(         +z(3)*xz31)
!   (pz px/pz px)
             w(78) = &
                  (ri(16)*xx31+ri(17)*zz31)*xx31 &
                  +(ri(18)*xx31+ri(19)*zz31)*zz31 &
                  +ri(20)*(xy31*xy31+xz31*xz31) &
                  +ri(22)*yz31*yz31
!   (pz px/pz py)
             w(79) = &
                  (ri(16)*xx31+ri(17)*zz31)*xx32 &
                  +(ri(18)*xx31+ri(19)*zz31)*zz32 &
                  +ri(20)*(xy31*xy32+xz31*xz32) &
                  +ri(22)*yz31*yz32
!   (pz px/pz pz)
             w(80) = &
                  (ri(16)*xx31+ri(17)*zz31)*xx33 &
                  +(ri(18)*xx31+ri(19)*zz31)*zz33 &
                  +ri(20)*xz31*xz33
!   (pz py/s s)
             w(81)=ri(3)*xx32+ri(4)*zz32
!   (pz py/px s)
             w(82)=x(1)*(ri(8)*xx32+ri(9)*zz32) &
                  +ri(10)*(y(1)*xy32+z(1)*xz32)
!   (pz py/px px)
             w(83) = &
                  (ri(16)*xx32+ri(17)*zz32)*xx11+ri(18)*xx32*yyzz11 &
                  +ri(19)*zz32*zz11 &
                  +ri(20)*(xy32*xy11+xz32*xz11) &
                  +ri(21)*zz32*yy11 &
                  +ri(22)*yz32*yz11
!   (pz py/py s)
             w(84)=x(2)*(ri(8)*xx32+ri(9)*zz32) &
                  +ri(10)*(y(2)*xy32+z(2)*xz32)
!   (pz py/py px)
             w(85) = &
                  (ri(16)*xx32+ri(17)*zz32)*xx21+ri(18)*xx32*yyzz21 &
                  +ri(19)*zz32*zz21 &
                  +ri(20)*(xy32*xy21+xz32*xz21) &
                  +ri(21)*zz32*yy21 &
                  +ri(22)*yz32*yz21
!   (pz py/py py)
             w(86) = &
                  (ri(16)*xx32+ri(17)*zz32)*xx22+ri(18)*xx32*yyzz22 &
                  +ri(19)*zz32*zz22 &
                  +ri(20)*(xy32*xy22+xz32*xz22) &
                  +ri(21)*zz32*yy22 &
                  +ri(22)*yz32*yz22
!   (pz py/pz s)
             w(87)=x(3)*(ri(8)*xx32+ri(9)*zz32) &
               +ri(10)*(         +z(3)*xz32)
!   (pz py/pz px)
             w(88) = &
                  (ri(16)*xx32+ri(17)*zz32)*xx31 &
                  +(ri(18)*xx32+ri(19)*zz32)*zz31 &
                  +ri(20)*(xy32*xy31+xz32*xz31) &
                  +ri(22)*yz32*yz31
!   (pz py/pz py)
             w(89) = &
                  (ri(16)*xx32+ri(17)*zz32)*xx32 &
                  +(ri(18)*xx32+ri(19)*zz32)*zz32 &
                  +ri(20)*(xy32*xy32+xz32*xz32) &
                  +ri(22)*yz32*yz32
!   (pz py/pz pz)
             w(90) = &
                  (ri(16)*xx32+ri(17)*zz32)*xx33 &
                  +(ri(18)*xx32+ri(19)*zz32)*zz33 &
                  +ri(20)*xz32*xz33
!   (pz pz/s s)
             w(91)=ri(3)*xx33+ri(4)*zz33
!   (pz pz/px s)
             w(92)=x(1)*(ri(8)*xx33+ri(9)*zz33) &
                  +ri(10)*(          z(1)*xz33)
!   (pz pz/px px)
             w(93) = &
                  (ri(16)*xx33+ri(17)*zz33)*xx11+ri(18)*xx33*yyzz11 &
                  +ri(19)*zz33*zz11 &
                  +ri(20)*xz33*xz11 &
                  +ri(21)*zz33*yy11
!   (pz pz/py s)
             w(94)=x(2)*(ri(8)*xx33+ri(9)*zz33) &
                  +ri(10)*(         +z(2)*xz33)
!   (pz pz/py px)
             w(95) = &
                  (ri(16)*xx33+ri(17)*zz33)*xx21+ri(18)*xx33*yyzz21 &
                  +ri(19)*zz33*zz21 &
                  +ri(20)*xz33*xz21 &
                  +ri(21)*zz33*yy21
!   (pz pz/py py)
             w(96) = &
                  (ri(16)*xx33+ri(17)*zz33)*xx22+ri(18)*xx33*yyzz22 &
                  +ri(19)*zz33*zz22 &
                  +ri(20)*xz33*xz22 &
                  +ri(21)*zz33*yy22
!   (pz pz/pz s)
             w(97)=x(3)*(ri(8)*xx33+ri(9)*zz33) &
                  +ri(10)*(         +z(3)*xz33)
!   (pz pz/pz px)
             w(98) = &
                  (ri(16)*xx33+ri(17)*zz33)*xx31 &
                  +(ri(18)*xx33+ri(19)*zz33)*zz31 &
                  +ri(20)*xz33*xz31
!   (pz pz/pz py)
             w(99) = &
                  (ri(16)*xx33+ri(17)*zz33)*xx32 &
                  +(ri(18)*xx33+ri(19)*zz33)*zz32 &
                  +ri(20)*xz33*xz32
!   (pz pz/pz pz)
             w(100) = &
                  (ri(16)*xx33+ri(17)*zz33)*xx33 &
                  +(ri(18)*xx33+ri(19)*zz33)*zz33 &
                  +ri(20)*xz33*xz33
          ELSE
!   (px s/s s)
             w(2)=ri(2)*x(1)
!   (px px/s s)
             w(3)=ri(3)*xx11+ri(4)*yyzz11
!   (py s/s s)
             w(4)=ri(2)*x(2)
!   (py px/s s)
             w(5)=ri(3)*xx21+ri(4)*yyzz21
!   (py py/s s)
             w(6)=ri(3)*xx22+ri(4)*yyzz22
!   (pz s/ss)
             w(7)=ri(2)*x(3)
!   (pz px/s s)
             w(8)=ri(3)*xx31+ri(4)*zz31
!   (pz py/s s)
             w(9)=ri(3)*xx32+ri(4)*zz32
!   (pz pz/s s)
             w(10)=ri(3)*xx33+ri(4)*zz33
          END IF
       END IF
    ENDIF
    
  END SUBROUTINE rotint

! *****************************************************************************
!> \brief Computes the two-particle interactions.
!> \param on input  sepi   = Atomic parameters of first atom
!>                  sepj   = Atomic parameters of second atom
!>                  rijv   = Coordinate vector i -> j
!>
!>    on output e1b,e2a= Array of electron-nuclear attraction integrals,
!>                       e1b = Electron on atom ni attracting nucleus of nj.
!>                       e2a = Electron on atom nj attracting nucleus of ni.
!>                       enuc= nuclear-nuclear repulsion term.
!>
!> \note routine adapted from mopac7 (rotate)
!>       written by Ernest R. Davidson, Indiana University.
! *****************************************************************************
  SUBROUTINE rotnuc (sepi,sepj,rijv,e1b,e2a,enuc,itype,shortrange,error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: e1b, e2a
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotnuc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ig, nt
    LOGICAL                                  :: si, sj
    REAL(dp) :: a, aab, alpi, alpj, apdg, ax, cppp1, cppp2, cpps1, cpps2, &
      csp1, csp2, css1, css2, dai, daj, dbi, dbj, pai, paj, pbi, pbj, poi, &
      poj, qcorr, rij, rija, scale, ssss, xab, xx11, xx21, xx22, xx31, xx32, &
      xx33, yy11, yy21, yy22, yyzz11, yyzz21, yyzz22, zaf, zbf, zz, zz11, &
      zz21, zz22, zz31, zz32, zz33
    REAL(dp), DIMENSION(3)                   :: x, y, z
    REAL(dp), DIMENSION(4)                   :: fni1, fni2, fni3, fnj1, fnj2, &
                                                fnj3
    REAL(dp), DIMENSION(4, 2)                :: ccore

    x(:)=-rijv(:)
    rij=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
    e1b(:)=0._dp
    e2a(:)=0._dp
    enuc=0._dp
    IF (rij > rij_threshold) THEN
       rij = SQRT(rij)
       ! Compute Integrals in Diatamic Frame
       CALL nucint(sepi,sepj,rij,ssss,ccore,itype,shortrange,error=error)

       a=1._dp/rij
       x(1) = x(1)*a
       x(2) = x(2)*a
       x(3) = x(3)*a
       IF (ABS(x(3)) > 0.99999999_dp) THEN
          x(3) = SIGN(1._dp,x(3))
          y(1) = 0._dp
          y(2) = 1._dp
          y(3) = 0._dp
          z(1) = 1._dp
          z(2) = 0._dp
          z(3) = 0._dp
       ELSE
          z(3)=SQRT(1._dp-x(3)*x(3))
          a=1._dp/z(3)
          y(1)=-a*x(2)*SIGN(1._dp,x(1))
          y(2)=ABS(a*x(1))
          y(3)=0._dp
          z(1)=-a*x(1)*x(3)
          z(2)=-a*x(2)*x(3)
       ENDIF
       si = (sepi%natorb > 1)
       sJ = (sepj%natorb > 1)
       IF ( si .OR. sj ) THEN
          xx11 = x(1)*x(1)
          xx21 = x(2)*x(1)
          xx22 = x(2)*x(2)
          xx31 = x(3)*x(1)
          xx32 = x(3)*x(2)
          xx33 = x(3)*x(3)
          yy11 = y(1)*y(1)
          YY21 = Y(2)*Y(1)
          yy22 = y(2)*y(2)
          zz11 = z(1)*z(1)
          zz21 = z(2)*z(1)
          zz22 = z(2)*z(2)
          zz31 = z(3)*z(1)
          zz32 = z(3)*z(2)
          zz33 = z(3)*z(3)
          yyzz11 = yy11+zz11
          yyzz21 = yy21+zz21
          yyzz22 = yy22+zz22
       ENDIF

       css1  = ccore(1,1)
       csp1  = ccore(2,1)
       cpps1 = ccore(3,1)
       cppp1 = ccore(4,1)
       css2  = ccore(1,2)
       csp2  = ccore(2,2)
       cpps2 = ccore(3,2)
       cppp2 = ccore(4,2)

       e1b(1)=-css1
       IF(sepi%natorb == 4) THEN
          e1b(2) = -csp1 *x(1)
          e1b(3) = -cpps1*xx11-cppp1*yyzz11
          e1b(4) = -csp1 *x(2)
          e1b(5) = -cpps1*xx21-cppp1*yyzz21
          e1b(6) = -cpps1*xx22-cppp1*yyzz22
          e1b(7) = -csp1 *x(3)
          e1b(8) = -cpps1*xx31-cppp1*zz31
          e1b(9) = -cpps1*xx32-cppp1*zz32
          e1b(10)= -cpps1*xx33-cppp1*zz33
       END IF
       e2a(1)=-css2
       IF(sepj%natorb == 4) THEN
          e2a(2) = -csp2 *x(1)
          e2a(3) = -cpps2*xx11-cppp2*yyzz11
          e2a(4) = -csp2 *x(2)
          e2a(5) = -cpps2*xx21-cppp2*yyzz21
          e2a(6) = -cpps2*xx22-cppp2*yyzz22
          e2a(7) = -csp2 *x(3)
          e2a(8) = -cpps2*xx31-cppp2*zz31
          e2a(9) = -cpps2*xx32-cppp2*zz32
          e2a(10)= -cpps2*xx33-cppp2*zz33
       END IF

       IF(itype /= do_method_pm6) THEN

          alpi = sepi%alp
          alpj = sepj%alp
          scale = EXP(-alpi*rij)+EXP(-alpj*rij)

          nt=sepi%z+sepj%z
          IF(nt == 8 .OR. nt == 9) THEN
             IF(sepi%z == 7.OR.sepi%z == 8) &
                  scale=scale+(0.5291772083_dp*rij-1._dp)*EXP(-alpi*rij)
             IF(sepj%z == 7.OR.sepj%z == 8) &
                  scale=scale+(0.5291772083_dp*rij-1._dp)*EXP(-alpj*rij)
          ENDIF
          zz = sepi%zeff*sepj%zeff
          enuc = zz*ssss
          scale=ABS(scale*enuc)
          zz=zz/rij
          IF(itype == do_method_am1 .OR. itype == do_method_pm3 .OR. itype == do_method_pdg) THEN
             IF(itype == do_method_am1 .AND. sepi%z == 5) THEN
                !special case AM1 Boron
                SELECT CASE (sepj%z)
                CASE DEFAULT
                   nt=1
                CASE (1)
                   nt=2
                CASE (6)
                   nt=3
                CASE (9,17,35,53)
                   nt=4
                END SELECT
                fni1(:) = sepi%bfn1(:,nt)
                fni2(:) = sepi%bfn2(:,nt)
                fni3(:) = sepi%bfn3(:,nt)
             ELSE
                fni1(:) = sepi%fn1(:)
                fni2(:) = sepi%fn2(:)
                fni3(:) = sepi%fn3(:)
             END IF
             IF(itype == do_method_am1 .AND. sepj%z == 5) THEN
                !special case AM1 Boron
                SELECT CASE (sepi%z)
                CASE DEFAULT
                   nt=1
                CASE (1)
                   nt=2
                CASE (6)
                   nt=3
                CASE (9,17,35,53)
                   nt=4
                END SELECT
                fnj1(:) = sepj%bfn1(:,nt)
                fnj2(:) = sepj%bfn2(:,nt)
                fnj3(:) = sepj%bfn3(:,nt)
             ELSE
                fnj1(:) = sepj%fn1(:)
                fnj2(:) = sepj%fn2(:)
                fnj3(:) = sepj%fn3(:)
             END IF
             ! AM1/PM3/PDG correction to nuclear repulsion
             DO ig=1,SIZE(fni1)
                IF(ABS(fni1(ig)) > 0._dp) THEN
                   AX = fni2(ig)*(rij-fni3(ig))**2
                   IF(ax <= 25._dp) THEN
                      scale=scale +zz*fni1(ig)*EXP(-ax)
                   ENDIF
                ENDIF
                IF(ABS(fnj1(ig)) > 0._dp) THEN
                   ax = fnj2(ig)*(rij-fnj3(ig))**2
                   IF(ax <= 25._dp) THEN
                      scale=scale +zz*fnj1(ig)*EXP(-ax)
                   ENDIF
                ENDIF
             END DO
          ENDIF
          IF(itype == do_method_pdg) THEN
             ! PDDG function 
             zaf=sepi%zeff/nt
             zbf=sepj%zeff/nt
             pai=sepi%pre(1)
             pbi=sepi%pre(2)
             paj=sepj%pre(1)
             pbj=sepj%pre(2)
             dai=sepi%d(1)
             dbi=sepi%d(2)
             daj=sepj%d(1)
             dbj=sepj%d(2)
             apdg=10._dp*0.5291772083_dp*0.5291772083_dp
             qcorr = (zaf*pai+zbf*paj)*EXP(-apdg*(rij-dai-daj)**2)+ &
                  (zaf*pai+zbf*pbj)*EXP(-apdg*(rij-dai-dbj)**2)+ &
                  (zaf*pbi+zbf*paj)*EXP(-apdg*(rij-dbi-daj)**2)+ &
                  (zaf*pbi+zbf*pbj)*EXP(-apdg*(rij-dbi-dbj)**2)
          ELSEIF (itype==do_method_pchg) THEN
             qcorr = 0.0_dp
             scale = 0.0_dp
          ELSE
             qcorr = 0.0_dp
          END IF
          enuc=enuc+scale+qcorr

       ELSE
          enuc = 0._dp
          rija = rij*angstrom
          IF ( sepi%rho > 0._dp ) THEN
             poi = sepi%rho
          ELSE
             poi = 0.5_dp/sepi%am
          END IF
          IF ( sepj%rho > 0._dp ) THEN
             poj = sepj%rho
          ELSE
             poj = 0.5_dp/sepj%am
          END IF
          scale = sepi%zeff*sepj%zeff/SQRT(rij*rij+(poi+poj)**2)
          IF(shortrange) scale = scale -  sepi%zeff*sepj%zeff/rij
          ! PM6 core-core terms
          xab = sepi%xab(sepj%z)
          aab = sepi%aab(sepj%z)
          IF ( (sepi%z == 1 .AND. (sepj%z == 6 .OR. sepj%z == 7 .OR. sepj%z == 8)) .OR. &
               (sepj%z == 1 .AND. (sepi%z == 6 .OR. sepi%z == 7 .OR. sepi%z == 8)) ) THEN
             ! Special Case O-H or N-H or C-H
             enuc=enuc+scale*(1._dp+2._dp*xab*EXP(-aab*rija*rija))
          ELSEIF ( sepi%z == 6 .AND. sepj%z == 6 ) THEN
             ! Special Case C-C
             enuc=enuc+scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6)) &
                  +9.28_dp*EXP(-5.98*rija))
          ELSEIF ( (sepi%z == 8 .AND. sepj%z == 14) .OR. &
               (sepj%z == 8 .AND. sepi%z == 14) ) THEN
             ! Special Case Si-O
             enuc=enuc+scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6)) &
                  -0.0007_dp*EXP(-(rija-2.9_dp)**2) )
          ELSE
             ! General Case
             ! Factor of 2 found by experiment 
             enuc=enuc+scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6)))
          END IF
          ! General correction term a*exp(-b*(rij-c)^2)
          enuc=enuc+(sepi%a*EXP(-sepi%b*(rij-sepi%c)**2))*sepi%zeff*sepj%zeff/rij
          enuc=enuc+(sepj%a*EXP(-sepj%b*(rij-sepj%c)**2))*sepi%zeff*sepj%zeff/rij
          ! Hard core repulsion
          enuc=enuc+1.e-8_dp/evolt &
               * ( (REAL(sepi%z,dp)**(1._dp/3._dp) + REAL(sepj%z,dp)**(1._dp/3._dp))/rija )**12
       END IF
    END IF

  END SUBROUTINE rotnuc

! *****************************************************************************
!> \brief Calculates the derivative pf two-electron repulsion integrals and the
!>      nuclear attraction integrals w.r.t. |r|
!> \param on input rij   = interatomic distance
!>               sepi    = paramters of atom i
!>               sepj    = paramters of atom j
!>        short_range    = requests only the evaluation of the short range part
!>
!>     on output ri      = array of two-electron repulsion integrals
!>
!> \par History
!>      03.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE terep ( sepi, sepj, rij, ri, shortrange, error )
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: ri
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'terep', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: di, dj, failure

    failure = .FALSE.
    di = (sepi%natorb > 4)
    dj = (sepj%natorb > 4)
    ri = 0.0_dp

    CALL terep_sp ( sepi, sepj, rij, ri, shortrange, error=error)
    IF (di.OR.dj) THEN
       ! compute the contribution from d-orbitals
       STOP "D-ORBITALS: terep"
    END IF
  END SUBROUTINE terep

! *****************************************************************************
!> \brief Calculates the  two-electron repulsion integrals and the
!>      nuclear attraction integrals w.r.t. |r| - sp core only
!> \param on input rij   = interatomic distance
!>               sepi    = paramters of atom i
!>               sepj    = paramters of atom j
!>        short_range    = requests only the evaluation of the short range part
!>        do_longrange   = requests only the evaluation of the long range part (recursive use only)
!>
!>     on output ri      = array of two-electron repulsion integrals
!>
!>        The two-centre repulsion integrals (over local coordinates) are
!>        stored as follows (where p-sigma = O,  and p-pi = P and P* )
!>        (SS/SS)=1,     (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!>        (SO/SO)=6,     (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!>        (SS/OO)=11,    (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!>        (OO/OO)=16,    (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!>        (PP/P*P*)=21,  (P*P/P*P)=22.
!> \note routine adapted from mopac7 (repp)
!>       vector version written by Ernest R. Davidson, Indiana University
!> \ar History
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
! *****************************************************************************
  RECURSIVE SUBROUTINE terep_sp ( sepi, sepj, rij, ri, shortrange, do_longrange, error )
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: ri
    LOGICAL, INTENT(IN)                      :: shortrange
    LOGICAL, INTENT(IN), OPTIONAL            :: do_longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'terep_sp', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: longrange, si, sj
    REAL(dp) :: ade, adi, adj, adq, aed, aee, aeq, ami, amj, aqd, aqe, aqi, &
      aqj, aqq, axx, da, db, ddi, ddj, dxdx, dxqxz, dzdz, dze, dzqxx, dzqzz, &
      edz, ee, eqxx, eqzz, ft, qa, qb, qqi, qqj, qxxdz, qxxe, qxxqxx, qxxqyy, &
      qxxqzz, qxzdx, qxzqxz, qzzdz, qzze, qzzqxx, qzzqzz, r, rsq, www, xxx, &
      yyy, zi, zj, zzz
    REAL(dp), DIMENSION(72)                  :: arg, ri_lr, sqr

    longrange = .FALSE.
    IF (PRESENT(do_longrange))  longrange = do_longrange
    r   = rij
    si  = (sepi%natorb > 1)
    sj  = (sepj%natorb > 1)
    zi  = sepi%zeff
    zj  = sepj%zeff
    ddi = sepi%dd
    ddj = sepj%dd
    qqi = sepi%qq
    qqj = sepj%qq
    rsq = r*r

    !  hydrogen - hydrogen  (SS/SS)
    ami = 0.5_dp/sepi%am
    amj = 0.5_dp/sepj%am
    aee = ami + amj
    aee = aee * aee
    IF (longrange) aee = 0.0_dp
    arg(1) = rsq + aee
    sqr(1) = pxy(1)/SQRT(arg(1))
    ee     = sqr(1)
    ri(1)  = ee 

    IF (si.OR.sj) THEN
       ! In case (i) is an heavy atom
       IF (si) THEN
          da  = ddi
          qa  = qqi * 2.0_dp
          IF (longrange) THEN
             ade = 0.0_dp
             aqe = 0.0_dp
          ELSE
             adi = 0.5_dp/sepi%ad
             aqi = 0.5_dp/sepi%aq
             ade = adi + amj
             ade = ade * ade
             aqe = aqi + amj
             aqe = aqe * aqe
          END IF
          xxx     = r+da
          arg(2)  = xxx*xxx + ade
          xxx     = r-da
          arg(3)  = xxx*xxx + ade
          xxx     = r+qa
          arg(4)  = xxx*xxx + aqe
          xxx     = r-qa
          arg(5)  = xxx*xxx + aqe
          arg(6)  = rsq + aqe
          arg(7)  = arg(6) + qa*qa
          sqr(2:7) = pxy(2:7)/SQRT(arg(2:7))
          dze      = sqr(2) + sqr(3)
          qzze     = sqr(4) + sqr(5) + sqr(6)
          qxxe     = sqr(7) + sqr(6)
          ri(2)    = dze
          ri(3)    = ee + qzze
          ri(4)    = ee + qxxe
       END IF

       ! In case (j) is an heavy atom
       IF (sj) THEN
          db  = ddj
          qb  = qqj * 2.0_dp
          IF (longrange) THEN
             aed = 0.0_dp
             aeq = 0.0_dp
          ELSE
             adj = 0.5_dp/sepj%ad
             aqj = 0.5_dp/sepj%aq
             aed = ami + adj
             aed = aed * aed
             aeq = ami + aqj
             aeq = aeq * aeq
          END IF
          xxx      = r-db
          arg(8)   = xxx*xxx + aed
          xxx      = r+db
          arg(9)   = xxx*xxx + aed
          xxx      = r-qb
          arg(10)  = xxx*xxx + aeq
          xxx      = r+qb
          arg(11)  = xxx*xxx + aeq
          arg(12)  = rsq + aeq
          arg(13)  = arg(12) + qb*qb
          sqr(8:13) = pxy(8:13)/SQRT(arg(8:13))
          edz       = sqr(8) + sqr(9)
          eqzz      = sqr(10) + sqr(11) + sqr(12)
          eqxx      = sqr(12) + sqr(13)
          ri(5)     = edz
          ri(11)    = ee + eqzz
          ri(12)    = ee + eqxx
       END IF

       ! In case both (i) and (j) are heavy atoms
       IF (si.AND.sj) THEN
          IF (longrange) THEN
             axx = 0.0_dp
             adq = 0.0_dp
             aqd = 0.0_dp
             aqq = 0.0_dp
          ELSE
             axx = adi + adj
             axx = axx * axx
             adq = adi + aqj
             adq = adq * adq
             aqd = aqi + adj
             aqd = aqd * aqd
             aqq = aqi + aqj
             aqq = aqq * aqq
          END IF
          xxx = da-db
          arg(14) = rsq + axx + xxx*xxx
          xxx = da+db
          arg(15) = rsq + axx + xxx*xxx
          xxx = r + da - db
          arg(16) = xxx*xxx + axx
          xxx = r - da + db
          arg(17) = xxx*xxx + axx
          xxx = r - da - db
          arg(18) = xxx*xxx + axx
          xxx = r + da + db
          arg(19) = xxx*xxx + axx
          xxx = r + da
          arg(20) = xxx*xxx + adq
          arg(21) = arg(20) + qb*qb
          xxx = r - da
          arg(22) = xxx*xxx + adq
          arg(23) = arg(22) + qb*qb
          xxx = r - db
          arg(24) = xxx*xxx + aqd
          arg(25) = arg(24) + qa*qa
          xxx = r + db
          arg(26) = xxx*xxx + aqd
          arg(27) = arg(26) + qa*qa
          xxx = r + da - qb
          arg(28) = xxx*xxx + adq
          xxx = r - da - qb
          arg(29) = xxx*xxx + adq
          xxx = r + da + qb
          arg(30) = xxx*xxx + adq
          xxx = r - da + qb
          arg(31) = xxx*xxx + adq
          xxx = r + qa - db
          arg(32) = xxx*xxx + aqd
          xxx = r + qa + db
          arg(33) = xxx*xxx + aqd
          xxx = r - qa - db
          arg(34) = xxx*xxx + aqd
          xxx = r - qa + db
          arg(35) = xxx*xxx + aqd
          arg(36) = rsq + aqq
          xxx = qa - qb
          arg(37) = arg(36) + xxx*xxx
          xxx = qa + qb
          arg(38) = arg(36) + xxx*xxx
          arg(39) = arg(36) + qa*qa
          arg(40) = arg(36) + qb*qb
          arg(41) = arg(39) + qb*qb
          xxx = r - qb
          arg(42) = xxx*xxx + aqq
          arg(43) = arg(42) + qa*qa
          xxx = r + qb
          arg(44) = xxx*xxx + aqq
          arg(45) = arg(44) + qa*qa
          xxx = r + qa
          arg(46) = xxx*xxx + aqq
          arg(47) = arg(46) + qb*qb
          xxx = r - qa
          arg(48) = xxx*xxx + aqq
          arg(49) = arg(48) + qb*qb
          xxx = r + qa - qb
          arg(50) = xxx*xxx + aqq
          xxx = r + qa + qb
          arg(51) = xxx*xxx + aqq
          xxx = r - qa - qb
          arg(52) = xxx*xxx + aqq
          xxx = r - qa + qb
          arg(53) = xxx*xxx + aqq
          qa=sepi%qq
          qb=sepj%qq
          xxx = da - qb
          xxx = xxx*xxx
          yyy = r - qb
          yyy = yyy*yyy
          zzz = da + qb
          zzz = zzz*zzz
          www = r + qb
          www = www*www
          arg(54) = xxx + yyy + adq
          arg(55) = xxx + www + adq
          arg(56) = zzz + yyy + adq
          arg(57) = zzz + www + adq
          xxx = qa - db
          xxx = xxx*xxx
          yyy = qa + db
          yyy = yyy*yyy
          zzz = r + qa
          zzz = zzz*zzz
          www = r - qa
          www = www*www
          arg(58) = zzz + xxx + aqd
          arg(59) = www + xxx + aqd
          arg(60) = zzz + yyy + aqd
          arg(61) = www + yyy + aqd
          xxx = qa - qb
          xxx = xxx*xxx
          arg(62) = arg(36) + 2._dp*xxx
          yyy = qa + qb
          yyy = yyy*yyy
          arg(63) = arg(36) + 2._dp*yyy
          arg(64) = arg(36) + 2._dp*(qa*qa+qb*qb)
          zzz = r + qa - qb
          zzz = zzz*zzz
          arg(65) = zzz + xxx + aqq
          arg(66) = zzz + yyy + aqq
          zzz = r + qa + qb
          zzz = zzz*zzz
          arg(67) = zzz + xxx + aqq
          arg(68) = zzz + yyy + aqq
          zzz = r - qa - qb
          zzz = zzz*zzz
          arg(69) = zzz + xxx + aqq
          arg(70) = zzz + yyy + aqq
          zzz = r - qa + qb
          zzz = zzz*zzz
          arg(71) = zzz + xxx + aqq
          arg(72) = zzz + yyy + aqq
          sqr(14:72) = pxy(14:72)/SQRT(arg(14:72))

          dxdx   = sqr(14) + sqr(15)
          dzdz   = sqr(16) + sqr(17) + sqr(18) + sqr(19)
          dzqxx  = sqr(20) + sqr(21) + sqr(22) + sqr(23)
          qxxdz  = sqr(24) + sqr(25) + sqr(26) + sqr(27)
          dzqzz  = sqr(28) + sqr(29) + sqr(30) + sqr(31) + sqr(22) + sqr(20)
          qzzdz  = sqr(32) + sqr(33) + sqr(34) + sqr(35) + sqr(24) + sqr(26)
          qxxqxx = sqr(37) + sqr(38) + sqr(39) + sqr(40) + sqr(36)
          qxxqyy = sqr(41) + sqr(39) + sqr(40) + sqr(36)
          qxxqzz = sqr(43) + sqr(45) + sqr(42) + sqr(44) + sqr(39) + sqr(36)
          qzzqxx = sqr(47) + sqr(49) + sqr(46) + sqr(48) + sqr(40) + sqr(36)
          qzzqzz = sqr(50) + sqr(51) + sqr(52) + sqr(53) + sqr(48) + sqr(46) + sqr(42) + sqr(44) + sqr(36)
          dxqxz  = sqr(54) + sqr(55) + sqr(56) + sqr(57)
          qxzdx  = sqr(58) + sqr(59) + sqr(60) + sqr(61)
          qxzqxz = sqr(65) + sqr(67) + sqr(69) + sqr(71) + sqr(66) + sqr(68) + sqr(70) + sqr(72)

          ri(6)  = dzdz
          ri(7)  = dxdx
          ri(8)  =  edz - qzzdz
          ri(9)  =  edz - qxxdz
          ri(10) =      - qxzdx
          ri(13) =  dze - dzqzz
          ri(14) =  dze - dzqxx
          ri(15) =      - dxqxz
          ri(16) = ee + eqzz + qzze + qzzqzz
          ri(17) = ee + eqzz + qxxe + qxxqzz
          ri(18) = ee + eqxx + qzze + qzzqxx
          ri(19) = ee + eqxx + qxxe + qxxqxx
          ri(20) =                    qxzqxz
          ri(21) = ee + eqxx + qxxe + qxxqyy
          ri(22) = 0.5_dp * (qxxqxx - qxxqyy)
       END IF
    END IF

    IF ( shortrange ) THEN
       CALL terep_sp ( sepi, sepj, rij, ri_lr, shortrange=.FALSE., do_longrange=.TRUE.,&
            error=error)
       ri(1)  =  ri(1) - ri_lr(1) 
       ri(2)  =  ri(2) - ri_lr(2)
       ri(3)  =  ri(3) - ri_lr(3)
       ri(4)  =  ri(4) - ri_lr(4)
       ri(5)  =  ri(5) - ri_lr(5)
       ri(6)  =  ri(6) - ri_lr(6)
       ri(7)  =  ri(7) - ri_lr(7)
       ri(8)  =  ri(8) - ri_lr(8)
       ri(9)  =  ri(9) - ri_lr(9)
       ri(10) = ri(10) - ri_lr(10) 
       ri(11) = ri(11) - ri_lr(11)
       ri(12) = ri(12) - ri_lr(12)
       ri(13) = ri(13) - ri_lr(13)
       ri(14) = ri(14) - ri_lr(14)
       ri(15) = ri(15) - ri_lr(15)
       ri(16) = ri(16) - ri_lr(16)
       ri(17) = ri(17) - ri_lr(17)
       ri(18) = ri(18) - ri_lr(18)
       ri(19) = ri(19) - ri_lr(19)
       ri(20) = ri(20) - ri_lr(20)
       ri(21) = ri(21) - ri_lr(21)
       ri(22) = ri(22) - ri_lr(22)
    END IF

    ! Tapering function
    IF (.NOT.longrange) THEN
       ft = taper ( rij, error )
       ri(:) = ft*ri(:)
    END IF

  END SUBROUTINE terep_sp

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals (main driver)
!> \param on input rij     = interatomic distance
!>                sepi    = paramters of atom i
!>                sepj    = paramters of atom j
!>                itype   = type of semi_empirical model
!>                          extension to the original routine to compute qm/mm integrals
!>      on output core    = derivative of 4 X 2 array of electron-core attraction integrals
!>                ssss    = derivative of (ssss) integral
!>                          derivatives are intended w.r.t. rij
!> \par History
!>      03.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE nucint ( sepi, sepj, rij, ssss, core, itype, shortrange, error )
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(OUT)                    :: ssss
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: core
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nucint', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: di, dj, failure

    failure = .FALSE.
    di = (sepi%natorb > 4)
    dj = (sepj%natorb > 4)
    core  = 0.0_dp

    CALL nucint_spcore(sepi, sepj, rij, ssss, core, itype, shortrange, error=error)
    IF (di.OR.dj) THEN
       ! Compute the contribution from d-orbitals
       CALL nucint_d(sepi, sepj, rij, core, itype, shortrange, error=error)
    END IF

  END SUBROUTINE nucint

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals for sp basis
!> \param on input rij   = interatomic distance
!>                sepi   = paramters of atom i
!>                sepj   = paramters of atom j
!>                itype  = type of semi_empirical model
!>                         extension to the original routine to compute qm/mm integrals
!>        short_range    = requests only the evaluation of the short range part
!>        do_longrange   = requests only the evaluation of the long range part (recursive use only)
!>      on output core   = 4 X 2 array of electron-core attraction integrals
!>                ssss   = (ssss) integral
!>                          derivatives are intended w.r.t. rij
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!>         where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!> \note routine adapted from mopac7 (repp)
!>       vector version written by Ernest R. Davidson, Indiana University
!> \par History
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
! *****************************************************************************
  RECURSIVE SUBROUTINE nucint_spcore ( sepi, sepj, rij, ssss, core, itype, shortrange,&
       do_longrange, error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(INOUT)                  :: ssss
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: core
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: shortrange
    LOGICAL, INTENT(IN), OPTIONAL            :: do_longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nucint_spcore', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: longrange, si, sj
    REAL(dp)                                 :: ade, adi, adj, aed, aee, aeq, &
                                                ami, amj, aqe, aqi, aqj, da, &
                                                db, ddi, ddj, ft, qa, qb, &
                                                qqi, qqj, r, rsq, ssss_lr, &
                                                xxx, zi, zj
    REAL(dp), DIMENSION(7)                   :: arg, ri, sqr
    REAL(KIND=dp), DIMENSION(4, 2)           :: core_lr

    longrange = .FALSE.
    IF (PRESENT(do_longrange))  longrange = do_longrange
    ri   = 0.0_dp
    r    = rij
    si   = (sepi%natorb > 1)
    sj   = (sepj%natorb > 1)
    zi   = sepi%zeff
    zj   = sepj%zeff
    ddi  = sepi%dd
    ddj  = sepj%dd
    qqi  = sepi%qq
    qqj  = sepj%qq

    !  hydrogen - hydrogen  (SS/SS)
    aee = 0.0_dp
    IF (itype/=do_method_pchg) THEN
       ami = 0.5_dp/sepi%am
       amj = 0.5_dp/sepj%am
       aee = ami + amj
       aee = aee * aee
    END IF
    IF (longrange)  aee = 0.0_dp
    rsq       = r*r
    arg(1)    = rsq + aee
    sqr(1)    = pxy(1)/SQRT(arg(1))

    ri(1)     = 1.0_dp/SQRT(rsq+aee)
    core(1,1) = zj*ri(1)
    core(1,2) = zi*ri(1)
    IF (si.OR.sj) THEN
       ! In case (i) is an heavy atom
       IF (si) THEN
          da  = ddi
          qa  = qqi * 2.0_dp
          ade = 0.0_dp
          aqe = 0.0_dp
          IF (itype /=do_method_pchg) THEN
             adi = 0.5_dp/sepi%ad
             aqi = 0.5_dp/sepi%aq
             ade = adi + amj
             aqe = aqi + amj
             ade = ade * ade
             aqe = aqe * aqe
          END IF
          IF (longrange) THEN
             ade = 0.0_dp
             aqe = 0.0_dp
          END IF
          xxx    = r+da
          arg(2) = xxx*xxx + ade
          xxx    = r-da
          arg(3) = xxx*xxx + ade
          xxx    = r+qa
          arg(4) = xxx*xxx + aqe
          xxx    = r-qa
          arg(5) = xxx*xxx + aqe
          arg(6) = rsq + aqe
          arg(7) = arg(6) + qa*qa
          sqr(2:7)  = pxy(2:7)/SQRT(arg(2:7))
          ri(2)     = sqr(2) + sqr(3)
          ri(3)     = sqr(1) + sqr(4) + sqr(5) + sqr(6)
          ri(4)     = sqr(1) + sqr(6) + sqr(7)
          core(2,1) = zj*ri(2)
          core(3,1) = zj*ri(3)
          core(4,1) = zj*ri(4)
       END IF

       ! In case (j) is an heavy atom
       IF (sj) THEN
          db  = ddj
          qb  = qqj * 2.0_dp
          ade = 0.0_dp
          aqe = 0.0_dp
          IF (itype /=do_method_pchg) THEN
             adj = 0.5_dp/sepj%ad
             aqj = 0.5_dp/sepj%aq
             aed = ami + adj
             aeq = ami + aqj
             aed = aed * aed
             aeq = aeq * aeq
          END IF
          IF (longrange) THEN
             ade = 0.0_dp
             aqe = 0.0_dp
          END IF
          xxx = r-db
          arg(2) = xxx*xxx + aed
          xxx = r+db
          arg(3) = xxx*xxx + aed
          xxx = r-qb
          arg(4) = xxx*xxx + aeq
          xxx = r+qb
          arg(5) = xxx*xxx + aeq
          arg(6) = rsq + aeq
          arg(7) = arg(6) + qb*qb
          sqr(2:7)  = pxy(2:7)/SQRT(arg(2:7))
          ri(2)     = sqr(2) + sqr(3)
          ri(3)     = sqr(1) + sqr(4) + sqr(5) + sqr(6)
          ri(4)     = sqr(1) + sqr(6) + sqr(7)
          core(2,2) = zi*ri(2)
          core(3,2) = zi*ri(3)
          core(4,2) = zi*ri(4)
       END IF
    END IF

    ! Store the value for (SS|SS)
    ssss = ri(1)

    ! Modify integrals in case of short-range request
    IF ( shortrange ) THEN
       CALL nucint_spcore( sepi, sepj, rij, ssss_lr, core_lr, itype, shortrange=.FALSE.,&
            do_longrange=.TRUE., error=error)
       core = core - core_lr
       ssss = ssss - ssss_lr
    END IF

    ! If no orbitals are present on the atom let's zero the integral contribution
    IF (sepi%natorb == -1) core(:,1) = 0.0_dp
    IF (sepj%natorb == -1) core(:,2) = 0.0_dp

    IF (.NOT.longrange) THEN
       ! Tapering function
       ft = 1.0_dp
       IF (itype /=do_method_pchg) ft = taper ( rij, error )
       core(:,:) = ft*core(:,:)
       ssss = ft*ssss
    END IF

  END SUBROUTINE nucint_spcore

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals involving d orbitals
!> \param on input rij   = interatomic distance
!>                sepi   = paramters of atom i
!>                sepj   = paramters of atom j
!>                itype  = type of semi_empirical model
!>                         extension to the original routine to compute qm/mm integrals
!>        short_range    = requests only the evaluation of the short range part
!>        do_longrange   = requests only the evaluation of the long range part (recursive use only)
!>      on output core   = 4 X 2 array of electron-core attraction integrals
!>                ssss   = (ssss) integral
!>                          derivatives are intended w.r.t. rij
!>
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!>         where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!> \author
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich 
! *****************************************************************************
  RECURSIVE SUBROUTINE nucint_d ( sepi, sepj, rij, core, itype, shortrange, do_longrange,&
       error)
    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: core
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: shortrange
    LOGICAL, INTENT(IN), OPTIONAL            :: do_longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nucint_d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(SIZE(core,1)==10,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(core,2)==2,cp_failure_level,routineP,error,failure)
    STOP "d-orbitals"

  END SUBROUTINE nucint_d
! *****************************************************************************
!> \brief Numerical Derivatives for rotint
! *****************************************************************************
  SUBROUTINE drotint (sepi,sepj,r,dw,delta,shortrange,error)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: r
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: dw
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: od
    REAL(dp), DIMENSION(100)                 :: wm, wp
    REAL(dp), DIMENSION(3)                   :: rr

    od = 0.5_dp/delta
    DO i=1,3
       rr=r
       rr(i)=rr(i)+delta
       CALL rotint (sepi,sepj,rr,wp,shortrange,error=error)
       rr(i)=rr(i)-2._dp*delta
       CALL rotint (sepi,sepj,rr,wm,shortrange,error=error)
       dw(:,i) = od * (wp(:) - wm(:))
    END DO

  END SUBROUTINE drotint

! *****************************************************************************
!> \brief Numerical Derivatives for rotnuc
! *****************************************************************************
  SUBROUTINE drotnuc (sepi,sepj,r,de1b,de2a,denuc,itype,delta,shortrange,error)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: r
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: de1b, de2a
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotnuc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: enucm, enucp, od
    REAL(dp), DIMENSION(10)                  :: e1m, e1p, e2m, e2p
    REAL(dp), DIMENSION(3)                   :: rr

    od = 0.5_dp/delta
    DO i=1,3
       rr=r
       rr(i)=rr(i)+delta
       CALL rotnuc (sepi,sepj,rr,e1p,e2p,enucp,itype,shortrange,error=error)
       rr(i)=rr(i)-2._dp*delta
       CALL rotnuc (sepi,sepj,rr,e1m,e2m,enucm,itype,shortrange,error=error)
       denuc(i) = od * (enucp - enucm)
       de1b(:,i) = od * (e1p(:) - e1m(:))
       de2a(:,i) = od * (e2p(:) - e2m(:))
    END DO

  END SUBROUTINE drotnuc

END MODULE semi_empirical_integrals
