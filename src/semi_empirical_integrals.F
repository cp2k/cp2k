!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004 CP2K developers group                                  !
!-----------------------------------------------------------------------------!
!!****** cp2k/semi_empirical_integrals [1.0] *
!!
!!   NAME
!!     semi_empirical_integrals
!!
!!   FUNCTION
!!     Integrals for semi-empiric methods
!!
!!   AUTHOR
!!     JGH (11.10.2004)
!!
!!   MODIFICATION HISTORY
!!     JGH (27.10.2004) : separate routine for nuclear attraction integrals
!!
!!***
!******************************************************************************

MODULE semi_empirical_integrals

  USE kinds,                           ONLY: dp
  USE semi_empirical_types,            ONLY: semi_empirical_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_integrals'

  PUBLIC :: rotint, rotnuc, drotint, drotnuc

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE rotint (sepi,sepj,rijv,w)
  !
  ! routine adapted from mopac7 (rotate)
  ! written by Ernest R. Davidson, Indiana University.
  !
  !   rotint calculates the two-particle interactions.
  !
  !   on input  sepi   = Atomic parameters of first atom
  !             sepj   = Atomic parameters of second atom
  !             rijv   = Coordinate vector i -> j 
  !
  ! on output w      = Array of two-electron repulsion integrals.
  !
  !***********************************************************************

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: w

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
      IF (rij < 0.00002_dp) THEN
         ! SMALL RIJ CASE
         ! w is zero 
      ELSE
         ! The repulsion integrals over molecular frame (w) are stored in the
         ! order in which they will later be used.  ie.  (i,j/k,l) where
         ! j.le.i  and  l.le.k     and l varies most rapidly and i least
         ! rapidly.  (anti-normal computer storage)
         rij = SQRT(rij)
 
         ! COMPUTE INTEGRALS IN DIATOMIC FRAME
         CALL terep(sepi,sepj,rij,ri)
 
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
 
!     (s s/s s)
         w(1)=ri(1)
         IF (sj) THEN
!     (s s/px s)
            w(2)=ri(5)*x(1)
!     (s s/px px)
            w(3)=ri(11)*xx11+ri(12)*yyzz11
!     (s s/py s)
            w(4)=ri(5)*x(2)
!     (s s/py px)
            w(5)=ri(11)*xx21+ri(12)*yyzz21
!     (s s/py py)
            w(6)=ri(11)*xx22+ri(12)*yyzz22
!     (s s/pz s)
            w(7)=ri(5)*x(3)
!     (s s/pz px)
            w(8)=ri(11)*xx31+ri(12)*zz31
!     (s s/pz py)
            w(9)=ri(11)*xx32+ri(12)*zz32
!     (s s/pz pz)
            w(10)=ri(11)*xx33+ri(12)*zz33
         END IF
 
         IF (si) THEN
!     (px s/s s)
            w(11)=ri(2)*x(1)
            IF (sj) THEN
!     (px s/px s)
               w(12)=ri(6)*xx11+ri(7)*yyzz11
!     (px s/px px)
               w(13)=x(1)*(ri(13)*xx11+ri(14)*yyzz11) &
                 +ri(15)*(y(1)*xy11+z(1)*xz11)
!     (px s/py s)
               w(14)=ri(6)*xx21+ri(7)*yyzz21
!     (px s/py px)
               w(15)=x(1)*(ri(13)*xx21+ri(14)*yyzz21) &
                 +ri(15)*(y(1)*xy21+z(1)*xz21)
!     (px s/py py)
               w(16)=x(1)*(ri(13)*xx22+ri(14)*yyzz22) &
                 +ri(15)*(y(1)*xy22+z(1)*xz22)
!     (px s/pz s)
               w(17)=ri(6)*xx31+ri(7)*zz31
!     (px s/pz px)
               w(18)=x(1)*(ri(13)*xx31+ri(14)*zz31) &
                 +ri(15)*(y(1)*xy31+z(1)*xz31)
!     (px s/pz py)
               w(19)=x(1)*(ri(13)*xx32+ri(14)*zz32) &
                 +ri(15)*(y(1)*xy32+z(1)*xz32)
!     (px s/pz pz)
               w(20)=x(1)*(ri(13)*xx33+ri(14)*zz33) &
                 +ri(15)*(          z(1)*xz33)
!     (px px/s s)
               w(21)=ri(3)*xx11+ri(4)*yyzz11
!     (px px/px s)
               w(22)=x(1)*(ri(8)*xx11+ri(9)*yyzz11) &
                 +ri(10)*(y(1)*xy11+z(1)*xz11)
!     (px px/px px)
               w(23) = &
           (ri(16)*xx11+ri(17)*yyzz11)*xx11+ri(18)*xx11*yyzz11 &
           +ri(19)*(yy11*yy11+zz11*zz11) &
           +ri(20)*(xy11*xy11+xz11*xz11) &
           +ri(21)*(yy11*zz11+zz11*yy11) &
           +ri(22)*yz11*yz11
!     (px px/py s)
               w(24)=x(2)*(ri(8)*xx11+ri(9)*yyzz11) &
                 +ri(10)*(y(2)*xy11+z(2)*xz11)
!     (px px/py px)
               w(25) = &
           (ri(16)*xx11+ri(17)*yyzz11)*xx21+ri(18)*xx11*yyzz21 &
           +ri(19)*(yy11*yy21+zz11*zz21) &
           +ri(20)*(xy11*xy21+xz11*xz21) &
           +ri(21)*(yy11*zz21+zz11*yy21) &
           +ri(22)*yz11*yz21
!     (px px/py py)
               w(26) = &
           (ri(16)*xx11+ri(17)*yyzz11)*xx22+ri(18)*xx11*yyzz22 &
           +ri(19)*(yy11*yy22+zz11*zz22) &
           +ri(20)*(xy11*xy22+xz11*xz22) &
           +ri(21)*(yy11*zz22+zz11*yy22) &
           +ri(22)*yz11*yz22
!     (px px/pz s)
               w(27)=x(3)*(ri(8)*xx11+ri(9)*yyzz11) &
                 +ri(10)*(         +z(3)*xz11)
!     (px px/pz px)
               w(28) = &
            (ri(16)*xx11+ri(17)*yyzz11)*xx31 &
           +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz31 &
           +ri(20)*(xy11*xy31+xz11*xz31) &
           +ri(22)*yz11*yz31
!     (px px/pz py)
               w(29) = &
            (ri(16)*xx11+ri(17)*yyzz11)*xx32 &
           +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz32 &
           +ri(20)*(xy11*xy32+xz11*xz32) &
           +ri(22)*yz11*yz32
!     (px px/pz pz)
               w(30) = &
            (ri(16)*xx11+ri(17)*yyzz11)*xx33 &
           +(ri(18)*xx11+ri(19)*zz11+ri(21)*yy11)*zz33 &
           +ri(20)*xz11*xz33
!     (py s/s s)
               w(31)=ri(2)*x(2)
!     (py s/px s)
               w(32)=ri(6)*xx21+ri(7)*yyzz21
!     (py s/px px)
               w(33)=x(2)*(ri(13)*xx11+ri(14)*yyzz11) &
                 +ri(15)*(y(2)*xy11+z(2)*xz11)
!     (py s/py s)
               w(34)=ri(6)*xx22+ri(7)*yyzz22
!     (py s/py px)
               w(35)=x(2)*(ri(13)*xx21+ri(14)*yyzz21) &
                 +ri(15)*(y(2)*xy21+z(2)*xz21)
!     (py s/py py)
               w(36)=x(2)*(ri(13)*xx22+ri(14)*yyzz22) &
                 +ri(15)*(y(2)*xy22+z(2)*xz22)
!     (py s/pz s)
               w(37)=ri(6)*xx32+ri(7)*zz32
!     (py s/pz px)
               w(38)=x(2)*(ri(13)*xx31+ri(14)*zz31) &
                 +ri(15)*(y(2)*xy31+z(2)*xz31)
!     (py s/pz py)
               w(39)=x(2)*(ri(13)*xx32+ri(14)*zz32) &
                 +ri(15)*(y(2)*xy32+z(2)*xz32)
!     (py s/pz pz)
               w(40)=x(2)*(ri(13)*xx33+ri(14)*zz33) &
                 +ri(15)*(         +z(2)*xz33)
!     (py px/s s)
               w(41)=ri(3)*xx21+ri(4)*yyzz21
!     (py px/px s)
               w(42)=x(1)*(ri(8)*xx21+ri(9)*yyzz21) &
                 +ri(10)*(y(1)*xy21+z(1)*xz21)
!     (py px/px px)
               w(43) = &
           (ri(16)*xx21+ri(17)*yyzz21)*xx11+ri(18)*xx21*yyzz11 &
           +ri(19)*(yy21*yy11+zz21*zz11) &
           +ri(20)*(xy21*xy11+xz21*xz11) &
           +ri(21)*(yy21*zz11+zz21*yy11) &
           +ri(22)*yz21*yz11
!     (py px/py s)
               w(44)=x(2)*(ri(8)*xx21+ri(9)*yyzz21) &
                 +ri(10)*(y(2)*xy21+z(2)*xz21)
!     (py px/py px)
               w(45) = &
           (ri(16)*xx21+ri(17)*yyzz21)*xx21+ri(18)*xx21*yyzz21 &
           +ri(19)*(yy21*yy21+zz21*zz21) &
           +ri(20)*(xy21*xy21+xz21*xz21) &
           +ri(21)*(yy21*zz21+zz21*yy21) &
           +ri(22)*yz21*yz21
!     (py px/py py)
               w(46) = &
           (ri(16)*xx21+ri(17)*yyzz21)*xx22+ri(18)*xx21*yyzz22 &
           +ri(19)*(yy21*yy22+zz21*zz22) &
           +ri(20)*(xy21*xy22+xz21*xz22) &
           +ri(21)*(yy21*zz22+zz21*yy22) &
           +ri(22)*yz21*yz22
!     (py px/pz s)
               w(47)=x(3)*(ri(8)*xx21+ri(9)*yyzz21) &
                 +ri(10)*(         +z(3)*xz21)
!      (py px/pz px)
               w(48) = &
           (ri(16)*xx21+ri(17)*yyzz21)*xx31 &
           +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz31 &
           +ri(20)*(xy21*xy31+xz21*xz31) &
           +ri(22)*yz21*yz31
!      (py px/pz py)
               w(49) = &
           (ri(16)*xx21+ri(17)*yyzz21)*xx32 &
           +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz32 &
           +ri(20)*(xy21*xy32+xz21*xz32) &
           +ri(22)*yz21*yz32
!      (py px/pz pz)
               w(50) = &
           (ri(16)*xx21+ri(17)*yyzz21)*xx33 &
           +(ri(18)*xx21+ri(19)*zz21+ri(21)*yy21)*zz33 &
           +ri(20)*xz21*xz33
!     (py py/s s)
               w(51)=ri(3)*xx22+ri(4)*yyzz22
!     (py py/px s)
               w(52)=x(1)*(ri(8)*xx22+ri(9)*yyzz22) &
                 +ri(10)*(y(1)*xy22+z(1)*xz22)
!      (py py/px px)
               w(53) = &
           (ri(16)*xx22+ri(17)*yyzz22)*xx11+ri(18)*xx22*yyzz11 &
           +ri(19)*(yy22*yy11+zz22*zz11) &
           +ri(20)*(xy22*xy11+xz22*xz11) &
           +ri(21)*(yy22*zz11+zz22*yy11) &
           +ri(22)*yz22*yz11
!     (py py/py s)
               w(54)=x(2)*(ri(8)*xx22+ri(9)*yyzz22) &
                 +ri(10)*(y(2)*xy22+z(2)*xz22)
!      (py py/py px)
               w(55) = &
           (ri(16)*xx22+ri(17)*yyzz22)*xx21+ri(18)*xx22*yyzz21 &
           +ri(19)*(yy22*yy21+zz22*zz21) &
           +ri(20)*(xy22*xy21+xz22*xz21) &
           +ri(21)*(yy22*zz21+zz22*yy21) &
           +ri(22)*yz22*yz21
!      (py py/py py)
               w(56) = &
           (ri(16)*xx22+ri(17)*yyzz22)*xx22+ri(18)*xx22*yyzz22 &
           +ri(19)*(yy22*yy22+zz22*zz22) &
           +ri(20)*(xy22*xy22+xz22*xz22) &
           +ri(21)*(yy22*zz22+zz22*yy22) &
           +ri(22)*yz22*yz22
!     (py py/pz s)
               w(57)=x(3)*(ri(8)*xx22+ri(9)*yyzz22) &
                 +ri(10)*(         +z(3)*xz22)
!      (py py/pz px)
               w(58) = &
           (ri(16)*xx22+ri(17)*yyzz22)*xx31 &
           +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz31 &
           +ri(20)*(xy22*xy31+xz22*xz31) &
           +ri(22)*yz22*yz31
!      (py py/pz py)
               w(59) = &
           (ri(16)*xx22+ri(17)*yyzz22)*xx32 &
           +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz32 &
           +ri(20)*(xy22*xy32+xz22*xz32) &
           +ri(22)*yz22*yz32
!      (py py/pz pz)
               w(60) = &
           (ri(16)*xx22+ri(17)*yyzz22)*xx33 &
           +(ri(18)*xx22+ri(19)*zz22+ri(21)*yy22)*zz33 &
           +ri(20)*xz22*xz33
!     (pz s/ss)
               w(61)=ri(2)*x(3)
!     (pz s/px s)
               w(62)=ri(6)*xx31+ri(7)*zz31
!     (pz s/px px)
               w(63)=x(3)*(ri(13)*xx11+ri(14)*yyzz11) &
                 +ri(15)*(         +z(3)*xz11)
!     (pz s/py s)
               w(64)=ri(6)*xx32+ri(7)*zz32
!     (pz s/py px)
               w(65)=x(3)*(ri(13)*xx21+ri(14)*yyzz21) &
                 +ri(15)*(         +z(3)*xz21)
!     (pz s/py py)
               w(66)=x(3)*(ri(13)*xx22+ri(14)*yyzz22) &
                 +ri(15)*(         +z(3)*xz22)
!     (pz s/pz s)
               w(67)=ri(6)*xx33+ri(7)*zz33
!     (pz s/pz px)
               w(68)=x(3)*(ri(13)*xx31+ri(14)*zz31) &
                 +ri(15)*(         +z(3)*xz31)
!     (pz s/pz py)
               w(69)=x(3)*(ri(13)*xx32+ri(14)*zz32) &
                 +ri(15)*(         +z(3)*xz32)
!     (pz s/pz pz)
               w(70)=x(3)*(ri(13)*xx33+ri(14)*zz33) &
                 +ri(15)*(         +z(3)*xz33)
!     (pz px/s s)
               w(71)=ri(3)*xx31+ri(4)*zz31
!     (pz px/px s)
               w(72)=x(1)*(ri(8)*xx31+ri(9)*zz31) &
                 +ri(10)*(y(1)*xy31+z(1)*xz31)
!      (pz px/px px)
               w(73) = &
           (ri(16)*xx31+ri(17)*zz31)*xx11+ri(18)*xx31*yyzz11 &
           +ri(19)*zz31*zz11 &
           +ri(20)*(xy31*xy11+xz31*xz11) &
           +ri(21)*zz31*yy11 &
           +ri(22)*yz31*yz11
!     (pz px/py s)
               w(74)=x(2)*(ri(8)*xx31+ri(9)*zz31) &
                 +ri(10)*(y(2)*xy31+z(2)*xz31)
!      (pz px/py px)
               w(75) = &
           (ri(16)*xx31+ri(17)*zz31)*xx21+ri(18)*xx31*yyzz21 &
           +ri(19)*zz31*zz21 &
           +ri(20)*(xy31*xy21+xz31*xz21) &
           +ri(21)*zz31*yy21 &
           +ri(22)*yz31*yz21
!      (pz px/py py)
               w(76) = &
           (ri(16)*xx31+ri(17)*zz31)*xx22+ri(18)*xx31*yyzz22 &
           +ri(19)*zz31*zz22 &
           +ri(20)*(xy31*xy22+xz31*xz22) &
           +ri(21)*zz31*yy22 &
           +ri(22)*yz31*yz22
!     (pz px/pz s)
               w(77)=x(3)*(ri(8)*xx31+ri(9)*zz31) &
                 +ri(10)*(         +z(3)*xz31)
!     (pz px/pz px)
               w(78) = &
            (ri(16)*xx31+ri(17)*zz31)*xx31 &
           +(ri(18)*xx31+ri(19)*zz31)*zz31 &
           +ri(20)*(xy31*xy31+xz31*xz31) &
           +ri(22)*yz31*yz31
!      (pz px/pz py)
               w(79) = &
            (ri(16)*xx31+ri(17)*zz31)*xx32 &
           +(ri(18)*xx31+ri(19)*zz31)*zz32 &
           +ri(20)*(xy31*xy32+xz31*xz32) &
           +ri(22)*yz31*yz32
!      (pz px/pz pz)
               w(80) = &
            (ri(16)*xx31+ri(17)*zz31)*xx33 &
           +(ri(18)*xx31+ri(19)*zz31)*zz33 &
           +ri(20)*xz31*xz33
!     (pz py/s s)
               w(81)=ri(3)*xx32+ri(4)*zz32
!     (pz py/px s)
               w(82)=x(1)*(ri(8)*xx32+ri(9)*zz32) &
                 +ri(10)*(y(1)*xy32+z(1)*xz32)
!      (pz py/px px)
               w(83) = &
           (ri(16)*xx32+ri(17)*zz32)*xx11+ri(18)*xx32*yyzz11 &
           +ri(19)*zz32*zz11 &
           +ri(20)*(xy32*xy11+xz32*xz11) &
           +ri(21)*zz32*yy11 &
           +ri(22)*yz32*yz11
!     (pz py/py s)
               w(84)=x(2)*(ri(8)*xx32+ri(9)*zz32) &
                 +ri(10)*(y(2)*xy32+z(2)*xz32)
!      (pz py/py px)
               w(85) = &
           (ri(16)*xx32+ri(17)*zz32)*xx21+ri(18)*xx32*yyzz21 &
           +ri(19)*zz32*zz21 &
           +ri(20)*(xy32*xy21+xz32*xz21) &
           +ri(21)*zz32*yy21 &
           +ri(22)*yz32*yz21
!      (pz py/py py)
               w(86) = &
           (ri(16)*xx32+ri(17)*zz32)*xx22+ri(18)*xx32*yyzz22 &
           +ri(19)*zz32*zz22 &
           +ri(20)*(xy32*xy22+xz32*xz22) &
           +ri(21)*zz32*yy22 &
           +ri(22)*yz32*yz22
!     (pz py/pz s)
               w(87)=x(3)*(ri(8)*xx32+ri(9)*zz32) &
                 +ri(10)*(         +z(3)*xz32)
!      (pz py/pz px)
               w(88) = &
            (ri(16)*xx32+ri(17)*zz32)*xx31 &
           +(ri(18)*xx32+ri(19)*zz32)*zz31 &
           +ri(20)*(xy32*xy31+xz32*xz31) &
           +ri(22)*yz32*yz31
!      (pz py/pz py)
               w(89) = &
            (ri(16)*xx32+ri(17)*zz32)*xx32 &
           +(ri(18)*xx32+ri(19)*zz32)*zz32 &
           +ri(20)*(xy32*xy32+xz32*xz32) &
           +ri(22)*yz32*yz32
!       (pz py/pz pz)
               w(90) = &
            (ri(16)*xx32+ri(17)*zz32)*xx33 &
           +(ri(18)*xx32+ri(19)*zz32)*zz33 &
           +ri(20)*xz32*xz33
!     (pz pz/s s)
               w(91)=ri(3)*xx33+ri(4)*zz33
!     (pz pz/px s)
               w(92)=x(1)*(ri(8)*xx33+ri(9)*zz33) &
                 +ri(10)*(          z(1)*xz33)
!       (pz pz/px px)
               w(93) = &
           (ri(16)*xx33+ri(17)*zz33)*xx11+ri(18)*xx33*yyzz11 &
           +ri(19)*zz33*zz11 &
           +ri(20)*xz33*xz11 &
           +ri(21)*zz33*yy11
!     (pz pz/py s)
               w(94)=x(2)*(ri(8)*xx33+ri(9)*zz33) &
                 +ri(10)*(         +z(2)*xz33)
!       (pz pz/py px)
               w(95) = &
           (ri(16)*xx33+ri(17)*zz33)*xx21+ri(18)*xx33*yyzz21 &
           +ri(19)*zz33*zz21 &
           +ri(20)*xz33*xz21 &
           +ri(21)*zz33*yy21
!       (pz pz/py py)
               w(96) = &
           (ri(16)*xx33+ri(17)*zz33)*xx22+ri(18)*xx33*yyzz22 &
           +ri(19)*zz33*zz22 &
           +ri(20)*xz33*xz22 &
           +ri(21)*zz33*yy22
!     (pz pz/pz s)
               w(97)=x(3)*(ri(8)*xx33+ri(9)*zz33) &
                 +ri(10)*(         +z(3)*xz33)
!       (pz pz/pz px)
               w(98) = &
            (ri(16)*xx33+ri(17)*zz33)*xx31 &
           +(ri(18)*xx33+ri(19)*zz33)*zz31 &
           +ri(20)*xz33*xz31
!       (pz pz/pz py)
               w(99) = &
            (ri(16)*xx33+ri(17)*zz33)*xx32 &
           +(ri(18)*xx33+ri(19)*zz33)*zz32 &
           +ri(20)*xz33*xz32
!       (pz pz/pz pz)
               w(100) = &
            (ri(16)*xx33+ri(17)*zz33)*xx33 &
           +(ri(18)*xx33+ri(19)*zz33)*zz33 &
           +ri(20)*xz33*xz33
      ELSE
!     (px s/s s)
               w(2)=ri(2)*x(1)
!     (px px/s s)
               w(3)=ri(3)*xx11+ri(4)*yyzz11
!     (py s/s s)
               w(4)=ri(2)*x(2)
!     (py px/s s)
               w(5)=ri(3)*xx21+ri(4)*yyzz21
!     (py py/s s)
               w(6)=ri(3)*xx22+ri(4)*yyzz22
!     (pz s/ss)
               w(7)=ri(2)*x(3)
!     (pz px/s s)
               w(8)=ri(3)*xx31+ri(4)*zz31
!     (pz py/s s)
               w(9)=ri(3)*xx32+ri(4)*zz32
!     (pz pz/s s)
               w(10)=ri(3)*xx33+ri(4)*zz33
            END IF
         END IF
 
      ENDIF

  END SUBROUTINE rotint

! *****************************************************************************

  SUBROUTINE rotnuc (sepi,sepj,rijv,e1b,e2a,enuc,itype)
  !
  ! routine adapted from mopac7 (rotate)
  ! written by Ernest R. Davidson, Indiana University.
  !
  !   rotint calculates the two-particle interactions.
  !
  !   on input  sepi   = Atomic parameters of first atom
  !             sepj   = Atomic parameters of second atom
  !             rijv   = Coordinate vector i -> j 
  !             itype  = 1  (MNDO)
  !                    = 2  (AM1)
  !                    = 3  (PM3)
  !
  ! on output e1b,e2a= Array of electron-nuclear attraction integrals,
  !                    e1b = Electron on atom ni attracting nucleus of nj.
  !                    e2a = Electron on atom nj attracting nucleus of ni.
  !           enuc   = nuclear-nuclear repulsion term.
  !
  !***********************************************************************

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: e1b, e2a
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype

    INTEGER                                  :: ig, nt
    LOGICAL                                  :: si, sj
    REAL(dp) :: a, alpi, alpj, ax, cppp1, cppp2, cpps1, cpps2, csp1, csp2, &
      css1, css2, rij, scale, ssss, xx11, xx21, xx22, xx31, xx32, xx33, yy11, &
      yy21, yy22, yyzz11, yyzz21, yyzz22, zz, zz11, zz21, zz22, zz31, zz32, &
      zz33
    REAL(dp), DIMENSION(3)                   :: x, y, z
    REAL(dp), DIMENSION(4)                   :: fni1, fni2, fni3, fnj1, fnj2, &
                                                fnj3
    REAL(dp), DIMENSION(4, 2)                :: ccore

      x(:)=-rijv(:)
      rij=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      e1b(:)=0._dp
      e2a(:)=0._dp
      enuc=0._dp
      IF (rij > 0.00002_dp) THEN
         rij = SQRT(rij)
         ! COMPUTE INTEGRALS IN DIATOMIC FRAME
         CALL nucint(sepi,sepj,rij,ssss,ccore)
 
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
         IF(itype == 2 .OR. itype == 3) THEN
           IF(itype == 2 .AND. sepi%z == 5) THEN
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
           IF(itype == 2 .AND. sepj%z == 5) THEN
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
           ! AM1 or PM3 correction to nuclear repulsion
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
         enuc=enuc+scale
 
      ENDIF

  END SUBROUTINE rotnuc
! *****************************************************************************

  SUBROUTINE terep ( sepi, sepj, rij, ri )
  !
  ! routine adapted from mopac7 (repp)
  ! vector version written by Ernest R. Davidson, Indiana University
  !       Calculates the two-electron repulsion integrals and the
  !       nuclear attraction integrals.
  !
  !   on input rij     = interatomic distance
  !            sepi    = paramters of atom i
  !            sepj    = paramters of atom j
  !  on output ri      = array of two-electron repulsion integrals
  !     The two-centre repulsion integrals (over local coordinates) are
  !     stored as follows (where p-sigma = O,  and p-pi = P and P* )
  !     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
  !     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
  !     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
  !     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
  !     (PP/P*P*)=21,   (P*P/P*P)=22.
  !
! *****************************************************************************

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: ri

    LOGICAL                                  :: si, sj
    REAL(dp) :: ade, adi, adj, adq, aed, aee, aeq, ami, amj, aqd, aqe, aqi, &
      aqj, aqq, axx, da, db, ddi, ddj, dxdx, dxqxz, dzdz, dze, dzqxx, dzqzz, &
      edz, ee, eqxx, eqzz, qa, qb, qqi, qqj, qxxdz, qxxe, qxxqxx, qxxqyy, &
      qxxqzz, qxzdx, qxzqxz, qzzdz, qzze, qzzqxx, qzzqzz, r, rsq, www, xxx, &
      yyy, zi, zj, zzz
    REAL(dp), DIMENSION(72)                  :: arg, sqr

      ri = 0._dp
      r=rij
      si = (sepi%natorb >= 3)
      sj = (sepj%natorb >= 3)
      zi = sepi%zeff
      zj = sepj%zeff
      ddi = sepi%dd
      ddj = sepj%dd
      qqi = sepi%qq
      qqj = sepj%qq
 
      IF ((.NOT.si) .AND. (.NOT.sj)) THEN
         !
         !  hydrogen - hydrogen  (SS/SS)
         !
         ami = sepi%am
         amj = sepj%am
         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         ri(1) = 1._dp/SQRT(r*r+aee)
 
      ELSE IF (si .AND. (.NOT.sj)) THEN
         !
         !  heavy atom - hydrogen
         !
         ami = sepi%am
         adi = sepi%ad
         aqi = sepi%aq
         amj = sepj%am
         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         da=ddi
         qa=qqi * 2._dp
         ade = 0.5_dp/adi + 0.5_dp/amj
         ade = ade * ade
         aqe = 0.5_dp/aqi + 0.5_dp/amj
         aqe = aqe * aqe
         rsq = r*r
         arg(1) = rsq + aee
         xxx = r+da
         arg(2) = xxx*xxx + ade
         xxx = r-da
         arg(3) = xxx*xxx + ade
         xxx = r+qa
         arg(4) = xxx*xxx + aqe
         xxx = r-qa
         arg(5) = xxx*xxx + aqe
         arg(6) = rsq + aqe
         arg(7) = arg(6) + qa*qa
         sqr(1:7) = SQRT(arg(1:7))
         ee = 1._dp/sqr(1)
         ri(1) = ee
         ri(2) = 0.5_dp/sqr(2) - 0.5_dp/sqr(3)
         ri(3) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ri(4) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
 
      ELSE IF ((.NOT.si).AND.sj) THEN
         !
         !  hydrogen - heavy atom
         !
         ami = sepi%am
         amj = sepj%am
         adj = sepj%ad
         aqj = sepj%aq
         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         db=ddj
         qb=qqj * 2._dp
         aed = 0.5_dp/ami + 0.5_dp/adj
         aed = aed * aed
         aeq = 0.5_dp/ami + 0.5_dp/aqj
         aeq = aeq * aeq
         rsq = r*r
         arg(1) = rsq + aee
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
         sqr(1:7) = SQRT(arg(1:7))
         ee = 1._dp/sqr(1)
         ri(1) = ee
         ri(5) = 0.5_dp/sqr(2)  - 0.5_dp/sqr(3)
         ri(11) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ri(12) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
 
      ELSE
         !
         !  heavy atom - heavy atom
         !
         !  define charge separations.
         da=ddi
         db=ddj
         qa=qqi * 2._dp
         qb=qqj * 2._dp
 
         ami = sepi%am
         amj = sepj%am
         adi = sepi%ad
         adj = sepj%ad
         aqi = sepi%aq
         aqj = sepj%aq

         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         ade = 0.5_dp/adi + 0.5_dp/amj
         ade = ade * ade
         aqe = 0.5_dp/aqi + 0.5_dp/amj
         aqe = aqe * aqe
         aed = 0.5_dp/ami + 0.5_dp/adj
         aed = aed * aed
         aeq = 0.5_dp/ami + 0.5_dp/aqj
         aeq = aeq * aeq
         axx = 0.5_dp/adi + 0.5_dp/adj
         axx = axx * axx
         adq = 0.5_dp/adi + 0.5_dp/aqj
         adq = adq * adq
         aqd = 0.5_dp/aqi + 0.5_dp/adj
         aqd = aqd * aqd
         aqq = 0.5_dp/aqi + 0.5_dp/aqj
         aqq = aqq * aqq
         rsq = r * r
         arg(1) = rsq + aee
         xxx = r + da
         arg(2) = xxx * xxx + ade
         xxx = r - da
         arg(3) = xxx*xxx + ade
         xxx = r - qa
         arg(4) = xxx*xxx + aqe
         xxx = r + qa
         arg(5) = xxx*xxx + aqe
         arg(6) = rsq + aqe
         arg(7) = arg(6) + qa*qa
         xxx = r-db
         arg(8) = xxx*xxx + aed
         xxx = r+db
         arg(9) = xxx*xxx + aed
         xxx = r - qB
         arg(10) = xxx*xxx + aeq
         xxx = r + qb
         arg(11) = xxx*xxx + aeq
         arg(12) = rsq + aeq
         arg(13) = arg(12) + qb*qb
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
         yyy = yyy*yyY
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
         sqr(1:72) = SQRT(arg(1:72))
         ee = 1._dp/sqr(1)
         dze = -0.5_dp/sqr(2) + 0.5_dp/sqr(3)
         qzze = 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         qxxe = 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         edz = - 0.5_dp/sqr(8) + 0.5_dp/sqr(9)
         eqzz  = 0.25_dp/sqr(10) + 0.25_dp/sqr(11) - 0.5_dp/sqr(12)
         eqxx  = 0.5_dp/sqr(13) - 0.5_dp/sqr(12)
         dxdx  = 0.5_dp/sqr(14) - 0.5_dp/sqr(15)
         dzdz  = 0.25_dp/sqr(16) + 0.25_dp/sqr(17) - 0.25_dp/sqr(18) &
             - 0.25_dp/sqr(19)
         dzqxx =  0.25_dp/sqr(20) - 0.25_dp/sqr(21) - 0.25_dp/sqr(22) &
             + 0.25_dp/sqr(23)
         qxxdz =  0.25_dp/sqr(24) - 0.25_dp/sqr(25) - 0.25_dp/sqr(26) &
             + 0.25_dp/sqr(27)
         dzqzz = -0.125_dp/sqr(28) + 0.125_dp/sqr(29) - 0.125_dp/sqr(30) &
             + 0.125_dp/sqr(31) - 0.25_dp/sqr(22) + 0.25_dp/sqr(20)
         qzzdz = -0.125_dp/sqr(32) + 0.125_dp/sqr(33) - 0.125_dp/sqr(34) &
             + 0.125_dp/sqr(35) + 0.25_dp/sqr(24) - 0.25_dp/sqr(26)
         qxxqxx = 0.125_dp/sqr(37) + 0.125_dp/sqr(38) - 0.25_dp/sqr(39) &
             - 0.25_dp/sqr(40) + 0.25_dp/sqr(36)
         qxxqyy = 0.25_dp/sqr(41) - 0.25_dp/sqr(39) - 0.25_dp/sqr(40) &
             + 0.25_dp/sqr(36)
         qxxqzz = 0.125_dp/sqr(43) + 0.125_dp/sqr(45) - 0.125_dp/sqr(42) &
             - 0.125_dp/sqr(44) - 0.25_dp/sqr(39) + 0.25_dp/sqr(36)
         qzzqxx = 0.125_dp/sqr(47) + 0.125_dp/sqr(49) - 0.125_dp/sqr(46) &
             - 0.125_dp/sqr(48) - 0.25_dp/sqr(40) + 0.25_dp/sqr(36)
         qzzqzz = 0.0625_dp/sqr(50) + 0.0625_dp/sqr(51) + 0.0625_dp/sqr(52) &
             + 0.0625_dp/sqr(53) - 0.125_dp/sqr(48) - 0.125_dp/sqr(46) &
             - 0.125_dp/sqr(42) - 0.125_dp/sqr(44) + 0.25_dp/sqr(36)
         dxqxz = -0.25_dp/sqr(54) + 0.25_dp/sqr(55) + 0.25_dp/sqr(56) &
             - 0.25_dp/sqr(57)
         qxzdx = -0.25_dp/sqr(58) + 0.25_dp/sqr(59) + 0.25_dp/sqr(60) &
             - 0.25_dp/sqr(61)
         qxzqxz = 0.125_dp/sqr(65) - 0.125_dp/sqr(67) - 0.125_dp/sqr(69) &
             + 0.125_dp/sqr(71) - 0.125_dp/sqr(66) + 0.125_dp/sqr(68) &
             + 0.125_dp/sqr(70) - 0.125_dp/sqr(72)
         ri(1) = ee
         ri(2) = -dze
         ri(3) = ee + qzze
         ri(4) = ee + qxxe
         ri(5) = -edz
         ri(6) = dzdz
         ri(7) = dxdx
         ri(8) = -edz -qzzdz
         ri(9) = -edz -qxxdz
         ri(10) = -qxzdx
         ri(11) =  ee + eqzz
         ri(12) =  ee + eqxx
         ri(13) = -dze -dzqzz
         ri(14) = -dze -dzqxx
         ri(15) = -dxqxz
         ri(16) = ee +eqzz +qzze +qzzqzz
         ri(17) = ee +eqzz +qxxe +qxxqzz
         ri(18) = ee +eqxx +qzze +qzzqxx
         ri(19) = ee +eqxx +qxxe +qxxqxx
         ri(20) = qxzqxz
         ri(21) = ee +eqxx +qxxe +qxxqyy
         ri(22) = 0.5_dp * (qxxqxx -qxxqyy)
 
      END IF

  END SUBROUTINE terep

! *****************************************************************************

  SUBROUTINE nucint ( sepi, sepj, rij, ssss, core )
  !
  ! routine adapted from mopac7 (repp)
  ! vector version written by Ernest R. Davidson, Indiana University
  !       Calculates the nuclear attraction integrals.
  !
  !   on input rij     = interatomic distance
  !            sepi    = paramters of atom i
  !            sepj    = paramters of atom j
  !  on output core    = 4 X 2 array of electron-core attraction integrals
  !            ssss    = (ssss) integral
  !     The storage of the nuclear attraction integrals  core(kl/ij) iS
  !     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
  !     where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
  !
  ! ***************************************************************************

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(OUT)                    :: ssss
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: core

    LOGICAL                                  :: si, sj
    REAL(dp) :: ade, adi, adj, aed, aee, aeq, ami, amj, aqe, aqi, aqj, da, &
      db, ddi, ddj, dze, edz, ee, eqxx, eqzz, qa, qb, qqi, qqj, qxxe, qzze, &
      r, rsq, xxx, zi, zj
    REAL(dp), DIMENSION(14)                  :: arg, ri, sqr

      ri = 0._dp
      core = 0._dp
      r=rij
      si = (sepi%natorb >= 3)
      sj = (sepj%natorb >= 3)
      zi = sepi%zeff
      zj = sepj%zeff
      ddi = sepi%dd
      ddj = sepj%dd
      qqi = sepi%qq
      qqj = sepj%qq
 
      IF ((.NOT.si) .AND. (.NOT.sj)) THEN
         !
         !  hydrogen - hydrogen  (SS/SS)
         !
         ami = sepi%am
         amj = sepj%am
         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         ri(1) = 1._dp/SQRT(r*r+aee)
         core(1,1) = zj*ri(1)
         core(1,2) = zi*ri(1)
 
      ELSE IF (si .AND. (.NOT.sj)) THEN
         !
         !  heavy atom - hydrogen
         !
         ami = sepi%am
         adi = sepi%ad
         aqi = sepi%aq
         amj = sepj%am
         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         da=ddi
         qa=qqi * 2._dp
         ade = 0.5_dp/adi + 0.5_dp/amj
         ade = ade * ade
         aqe = 0.5_dp/aqi + 0.5_dp/amj
         aqe = aqe * aqe
         rsq = r*r
         arg(1) = rsq + aee
         xxx = r+da
         arg(2) = xxx*xxx + ade
         xxx = r-da
         arg(3) = xxx*xxx + ade
         xxx = r+qa
         arg(4) = xxx*xxx + aqe
         xxx = r-qa
         arg(5) = xxx*xxx + aqe
         arg(6) = rsq + aqe
         arg(7) = arg(6) + qa*qa
         sqr(1:7) = SQRT(arg(1:7))
         ee = 1._dp/sqr(1)
         ri(1) = ee
         ri(2) = 0.5_dp/sqr(2) - 0.5_dp/sqr(3)
         ri(3) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ri(4) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         core(1,1) = zj*ri(1)
         core(1,2) = zi*ri(1)
         core(2,1) = zj*ri(2)
         core(3,1) = zj*ri(3)
         core(4,1) = zj*ri(4)
 
      ELSE IF ((.NOT.si).AND.sj) THEN
         !
         !  hydrogen - heavy atom
         !
         ami = sepi%am
         amj = sepj%am
         adj = sepj%ad
         aqj = sepj%aq
         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         db=ddj
         qb=qqj * 2._dp
         aed = 0.5_dp/ami + 0.5_dp/adj
         aed = aed * aed
         aeq = 0.5_dp/ami + 0.5_dp/aqj
         aeq = aeq * aeq
         rsq = r*r
         arg(1) = rsq + aee
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
         sqr(1:7) = SQRT(arg(1:7))
         ee = 1._dp/sqr(1)
         ri(1) = ee
         ri(5) = 0.5_dp/sqr(2)  - 0.5_dp/sqr(3)
         ri(11) = ee + 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         ri(12) = ee + 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         core(1,1) = zj*ri(1)
         core(1,2) = zi*ri(1)
         core(2,2) = zi*ri(5)
         core(3,2) = zi*ri(11)
         core(4,2) = zi*ri(12)
 
      ELSE
         !
         !  heavy atom - heavy atom
         !
         !  define charge separations.
         da=ddi
         db=ddj
         qa=qqi * 2._dp
         qb=qqj * 2._dp
 
         ami = sepi%am
         amj = sepj%am
         adi = sepi%ad
         adj = sepj%ad
         aqi = sepi%aq
         aqj = sepj%aq

         aee = 0.5_dp/ami + 0.5_dp/amj
         aee = aee * aee
         ade = 0.5_dp/adi + 0.5_dp/amj
         ade = ade * ade
         aqe = 0.5_dp/aqi + 0.5_dp/amj
         aqe = aqe * aqe
         aed = 0.5_dp/ami + 0.5_dp/adj
         aed = aed * aed
         aeq = 0.5_dp/ami + 0.5_dp/aqj
         aeq = aeq * aeq
         rsq = r * r
         arg(1) = rsq + aee
         xxx = r + da
         arg(2) = xxx * xxx + ade
         xxx = r - da
         arg(3) = xxx*xxx + ade
         xxx = r - qa
         arg(4) = xxx*xxx + aqe
         xxx = r + qa
         arg(5) = xxx*xxx + aqe
         arg(6) = rsq + aqe
         arg(7) = arg(6) + qa*qa
         xxx = r-db
         arg(8) = xxx*xxx + aed
         xxx = r+db
         arg(9) = xxx*xxx + aed
         xxx = r - qb
         arg(10) = xxx*xxx + aeq
         xxx = r + qb
         arg(11) = xxx*xxx + aeq
         arg(12) = rsq + aeq
         arg(13) = arg(12) + qb*qb
         sqr(1:13) = SQRT(arg(1:13))
         ee = 1._dp/sqr(1)
         dze = -0.5_dp/sqr(2) + 0.5_dp/sqr(3)
         qzze = 0.25_dp/sqr(4) + 0.25_dp/sqr(5) - 0.5_dp/sqr(6)
         qxxe = 0.5_dp/sqr(7) - 0.5_dp/sqr(6)
         edz = - 0.5_dp/sqr(8) + 0.5_dp/sqr(9)
         eqzz  = 0.25_dp/sqr(10) + 0.25_dp/sqr(11) - 0.5_dp/sqr(12)
         eqxx  = 0.5_dp/sqr(13) - 0.5_dp/sqr(12)
         ri(1) = ee
         ri(2) = -dze
         ri(3) = ee + qzze
         ri(4) = ee + qxxe
         ri(5) = -edz
         ri(11) =  ee + eqzz
         ri(12) =  ee + eqxx
 
         ! Calculate core-electron attractions.
 
         core(1,1) = zj*ri(1)
         core(2,1) = zj*ri(2)
         core(3,1) = zj*ri(3)
         core(4,1) = zj*ri(4)
         core(1,2) = zi*ri(1)
         core(2,2) = zi*ri(5)
         core(3,2) = zi*ri(11)
         core(4,2) = zi*ri(12)
 
      END IF

      ssss = ri(1)

  END SUBROUTINE nucint

! *****************************************************************************

  SUBROUTINE drotint (sepi,sepj,r,dw,delta)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: r
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: dw
    REAL(dp), INTENT(IN)                     :: delta

    INTEGER                                  :: i
    REAL(dp)                                 :: od
    REAL(dp), DIMENSION(100)                 :: wm, wp
    REAL(dp), DIMENSION(3)                   :: rr

      od = 0.5_dp/delta
      DO i=1,3
        rr=r
        rr(i)=rr(i)+delta
        CALL rotint (sepi,sepj,rr,wp)
        rr(i)=rr(i)-2._dp*delta
        CALL rotint (sepi,sepj,rr,wm)
        dw(:,i) = od * (wp(:) - wm(:))
      END DO

  END SUBROUTINE drotint


  SUBROUTINE drotnuc (sepi,sepj,r,de1b,de2a,denuc,itype,delta)

    TYPE(semi_empirical_type), INTENT(IN)    :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(IN)       :: r
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: de1b, de2a
    REAL(dp), DIMENSION(:), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta

    INTEGER                                  :: i
    REAL(dp)                                 :: enucm, enucp, od
    REAL(dp), DIMENSION(10)                  :: e1m, e1p, e2m, e2p
    REAL(dp), DIMENSION(3)                   :: rr

      od = 0.5_dp/delta
      DO i=1,3
        rr=r
        rr(i)=rr(i)+delta
        CALL rotnuc (sepi,sepj,rr,e1p,e2p,enucp,itype)
        rr(i)=rr(i)-2._dp*delta
        CALL rotnuc (sepi,sepj,rr,e1m,e2m,enucm,itype)
        denuc(i) = od * (enucp - enucm)
        de1b(:,i) = od * (e1p(:) - e1m(:))
        de2a(:,i) = od * (e2p(:) - e2m(:))
      END DO

  END SUBROUTINE drotnuc

! *****************************************************************************

END MODULE semi_empirical_integrals
