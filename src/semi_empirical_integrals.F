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
!!
!!***
!******************************************************************************

MODULE semi_empirical_integrals

  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE semi_empirical_types,            ONLY: semi_empirical_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "semi_empirical_integrals"

  PUBLIC :: rotint

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE rotint (sepi,sepj,rijv,w,e1b,e2a,enuc,itype)
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
  ! on output w      = Array of two-electron repulsion integrals.
  !           e1b,e2a= Array of electron-nuclear attraction integrals,
  !                    e1b = Electron on atom ni attracting nucleus of nj.
  !                    e2a = Electron on atom nj attracting nucleus of ni.
  !           ENUC   = NUCLEAR-NUCLEAR REPULSION TERM.
  !
  !***********************************************************************

      TYPE ( semi_empirical_type ), INTENT ( IN )      :: sepi, sepj
      REAL ( dp ), DIMENSION (:),INTENT ( IN )         :: rijv
      REAL ( dp ), DIMENSION (:), INTENT ( OUT )       :: w,e1b,e2a
      REAL ( dp ), INTENT ( OUT )                      :: enuc
      INTEGER, INTENT ( IN )                           :: itype

      REAL ( dp ), DIMENSION (22)  :: ri
      REAL ( dp ), DIMENSION (4,2) :: ccore
      REAL ( dp ), DIMENSION (3)   :: x,y,z

      REAL ( dp ) :: rij,gam,a
      REAL ( dp ) :: xx11,xx21,xx22,xx31,xx32,xx33
      REAL ( dp ) :: yy11,yy21,yy22,yy31,yy32,yy33
      REAL ( dp ) :: zz11,zz21,zz22,zz31,zz32,zz33
      REAL ( dp ) :: xy11,xy21,xy22,xy31,xy32,xy33
      REAL ( dp ) :: xz11,xz21,xz22,xz31,xz32,xz33
      REAL ( dp ) :: yz11,yz21,yz22,yz31,yz32,yz33
      REAL ( dp ) :: yyzz11,yyzz21,yyzz22
      REAL ( dp ) :: CSS1,CSP1,CPPS1,CPPP1,CSS2,CSP2,CPPS2,CPPP2
      REAL ( dp ) :: scale,alpi,alpj,zz,ax

      LOGICAL     :: si,sj
      INTEGER     :: ki,nt,ig

      x(:)=-rijv(:)
      rij=x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
      e1b(:)=0._dp
      e2a(:)=0._dp
      IF (rij < 0.00002_dp) THEN
         ! SMALL RIJ CASE
         w(:)=0._dp
         enuc=0._dp
 
      ELSE
         ! The repulsion integrals over molecular frame (w) are stored in the
         ! order in which they will later be used.  ie.  (i,j/k,l) where
         ! j.le.i  and  l.le.k     and l varies most rapidly and i least
         ! rapidly.  (anti-normal computer storage)
         rij = SQRT(rij)
 
         ! COMPUTE INTEGRALS IN DIATOMIC FRAME
         CALL TEREP(sepi,sepj,rij,ri,ccore)
 
         gam = ri(1)
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
 
!     (S S/S S)
         W(1)=RI(1)
         IF (SJ) THEN
!     (S S/PX S)
            W(2)=RI(5)*X(1)
!     (S S/PX PX)
            W(3)=RI(11)*XX11+RI(12)*YYZZ11
!     (S S/PY S)
            W(4)=RI(5)*X(2)
!     (S S/PY PX)
            W(5)=RI(11)*XX21+RI(12)*YYZZ21
!     (S S/PY PY)
            W(6)=RI(11)*XX22+RI(12)*YYZZ22
!     (S S/PZ S)
            W(7)=RI(5)*X(3)
!     (S S/PZ PX)
            W(8)=RI(11)*XX31+RI(12)*ZZ31
!     (S S/PZ PY)
            W(9)=RI(11)*XX32+RI(12)*ZZ32
!     (S S/PZ PZ)
            W(10)=RI(11)*XX33+RI(12)*ZZ33
         ENDIF
 
         IF (SI) THEN
!     (PX S/S S)
            W(11)=RI(2)*X(1)
            IF (SJ) THEN
!     (PX S/PX S)
               W(12)=RI(6)*XX11+RI(7)*YYZZ11
!     (PX S/PX PX)
               W(13)=X(1)*(RI(13)*XX11+RI(14)*YYZZ11) &
                 +RI(15)*(Y(1)*XY11+Z(1)*XZ11)
!     (PX S/PY S)
               W(14)=RI(6)*XX21+RI(7)*YYZZ21
!     (PX S/PY PX)
               W(15)=X(1)*(RI(13)*XX21+RI(14)*YYZZ21) &
                 +RI(15)*(Y(1)*XY21+Z(1)*XZ21)
!     (PX S/PY PY)
               W(16)=X(1)*(RI(13)*XX22+RI(14)*YYZZ22) &
                 +RI(15)*(Y(1)*XY22+Z(1)*XZ22)
!     (PX S/PZ S)
               W(17)=RI(6)*XX31+RI(7)*ZZ31
!     (PX S/PZ PX)
               W(18)=X(1)*(RI(13)*XX31+RI(14)*ZZ31) &
                 +RI(15)*(Y(1)*XY31+Z(1)*XZ31)
!     (PX S/PZ PY)
               W(19)=X(1)*(RI(13)*XX32+RI(14)*ZZ32) &
                 +RI(15)*(Y(1)*XY32+Z(1)*XZ32)
!     (PX S/PZ PZ)
               W(20)=X(1)*(RI(13)*XX33+RI(14)*ZZ33) &
                 +RI(15)*(          Z(1)*XZ33)
!     (PX PX/S S)
               W(21)=RI(3)*XX11+RI(4)*YYZZ11
!     (PX PX/PX S)
               W(22)=X(1)*(RI(8)*XX11+RI(9)*YYZZ11) &
                 +RI(10)*(Y(1)*XY11+Z(1)*XZ11)
!     (PX PX/PX PX)
               W(23) = &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX11+RI(18)*XX11*YYZZ11 &
           +RI(19)*(YY11*YY11+ZZ11*ZZ11) &
           +RI(20)*(XY11*XY11+XZ11*XZ11) &
           +RI(21)*(YY11*ZZ11+ZZ11*YY11) &
           +RI(22)*YZ11*YZ11
!     (PX PX/PY S)
               W(24)=X(2)*(RI(8)*XX11+RI(9)*YYZZ11) &
                 +RI(10)*(Y(2)*XY11+Z(2)*XZ11)
!     (PX PX/PY PX)
               W(25) = &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX21+RI(18)*XX11*YYZZ21 &
           +RI(19)*(YY11*YY21+ZZ11*ZZ21) &
           +RI(20)*(XY11*XY21+XZ11*XZ21) &
           +RI(21)*(YY11*ZZ21+ZZ11*YY21) &
           +RI(22)*YZ11*YZ21
!     (PX PX/PY PY)
               W(26) = &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX22+RI(18)*XX11*YYZZ22 &
           +RI(19)*(YY11*YY22+ZZ11*ZZ22) &
           +RI(20)*(XY11*XY22+XZ11*XZ22) &
           +RI(21)*(YY11*ZZ22+ZZ11*YY22) &
           +RI(22)*YZ11*YZ22
!     (PX PX/PZ S)
               W(27)=X(3)*(RI(8)*XX11+RI(9)*YYZZ11) &
                 +RI(10)*(         +Z(3)*XZ11)
!     (PX PX/PZ PX)
               W(28) = &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX31 &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ31 &
           +RI(20)*(XY11*XY31+XZ11*XZ31) &
           +RI(22)*YZ11*YZ31
!     (PX PX/PZ PY)
               W(29) = &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX32 &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ32 &
           +RI(20)*(XY11*XY32+XZ11*XZ32) &
           +RI(22)*YZ11*YZ32
!     (PX PX/PZ PZ)
               W(30) = &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX33 &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ33 &
           +RI(20)*XZ11*XZ33
!     (PY S/S S)
               W(31)=RI(2)*X(2)
!     (PY S/PX S)
               W(32)=RI(6)*XX21+RI(7)*YYZZ21
!     (PY S/PX PX)
               W(33)=X(2)*(RI(13)*XX11+RI(14)*YYZZ11) &
                 +RI(15)*(Y(2)*XY11+Z(2)*XZ11)
!     (PY S/PY S)
               W(34)=RI(6)*XX22+RI(7)*YYZZ22
!     (PY S/PY PX)
               W(35)=X(2)*(RI(13)*XX21+RI(14)*YYZZ21) &
                 +RI(15)*(Y(2)*XY21+Z(2)*XZ21)
!     (PY S/PY PY)
               W(36)=X(2)*(RI(13)*XX22+RI(14)*YYZZ22) &
                 +RI(15)*(Y(2)*XY22+Z(2)*XZ22)
!     (PY S/PZ S)
               W(37)=RI(6)*XX32+RI(7)*ZZ32
!     (PY S/PZ PX)
               W(38)=X(2)*(RI(13)*XX31+RI(14)*ZZ31) &
                 +RI(15)*(Y(2)*XY31+Z(2)*XZ31)
!     (PY S/PZ PY)
               W(39)=X(2)*(RI(13)*XX32+RI(14)*ZZ32) &
                 +RI(15)*(Y(2)*XY32+Z(2)*XZ32)
!     (PY S/PZ PZ)
               W(40)=X(2)*(RI(13)*XX33+RI(14)*ZZ33) &
                 +RI(15)*(         +Z(2)*XZ33)
!     (PY PX/S S)
               W(41)=RI(3)*XX21+RI(4)*YYZZ21
!     (PY PX/PX S)
               W(42)=X(1)*(RI(8)*XX21+RI(9)*YYZZ21) &
                 +RI(10)*(Y(1)*XY21+Z(1)*XZ21)
!     (PY PX/PX PX)
               W(43) = &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX11+RI(18)*XX21*YYZZ11 &
           +RI(19)*(YY21*YY11+ZZ21*ZZ11) &
           +RI(20)*(XY21*XY11+XZ21*XZ11) &
           +RI(21)*(YY21*ZZ11+ZZ21*YY11) &
           +RI(22)*YZ21*YZ11
!     (PY PX/PY S)
               W(44)=X(2)*(RI(8)*XX21+RI(9)*YYZZ21) &
                 +RI(10)*(Y(2)*XY21+Z(2)*XZ21)
!     (PY PX/PY PX)
               W(45) = &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX21+RI(18)*XX21*YYZZ21 &
           +RI(19)*(YY21*YY21+ZZ21*ZZ21) &
           +RI(20)*(XY21*XY21+XZ21*XZ21) &
           +RI(21)*(YY21*ZZ21+ZZ21*YY21) &
           +RI(22)*YZ21*YZ21
!     (PY PX/PY PY)
               W(46) = &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX22+RI(18)*XX21*YYZZ22 &
           +RI(19)*(YY21*YY22+ZZ21*ZZ22) &
           +RI(20)*(XY21*XY22+XZ21*XZ22) &
           +RI(21)*(YY21*ZZ22+ZZ21*YY22) &
           +RI(22)*YZ21*YZ22
!     (PY PX/PZ S)
               W(47)=X(3)*(RI(8)*XX21+RI(9)*YYZZ21) &
                 +RI(10)*(         +Z(3)*XZ21)
!      (PY PX/PZ PX)
               W(48) = &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX31 &
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ31 &
           +RI(20)*(XY21*XY31+XZ21*XZ31) &
           +RI(22)*YZ21*YZ31
!      (PY PX/PZ PY)
               W(49) = &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX32 &
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ32 &
           +RI(20)*(XY21*XY32+XZ21*XZ32) &
           +RI(22)*YZ21*YZ32
!      (PY PX/PZ PZ)
               W(50) = &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX33 &
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ33 &
           +RI(20)*XZ21*XZ33
!     (PY PY/S S)
               W(51)=RI(3)*XX22+RI(4)*YYZZ22
!     (PY PY/PX S)
               W(52)=X(1)*(RI(8)*XX22+RI(9)*YYZZ22) &
                 +RI(10)*(Y(1)*XY22+Z(1)*XZ22)
!      (PY PY/PX PX)
               W(53) = &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX11+RI(18)*XX22*YYZZ11 &
           +RI(19)*(YY22*YY11+ZZ22*ZZ11) &
           +RI(20)*(XY22*XY11+XZ22*XZ11) &
           +RI(21)*(YY22*ZZ11+ZZ22*YY11) &
           +RI(22)*YZ22*YZ11
!     (PY PY/PY S)
               W(54)=X(2)*(RI(8)*XX22+RI(9)*YYZZ22) &
                 +RI(10)*(Y(2)*XY22+Z(2)*XZ22)
!      (PY PY/PY PX)
               W(55) = &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX21+RI(18)*XX22*YYZZ21 &
           +RI(19)*(YY22*YY21+ZZ22*ZZ21) &
           +RI(20)*(XY22*XY21+XZ22*XZ21) &
           +RI(21)*(YY22*ZZ21+ZZ22*YY21) &
           +RI(22)*YZ22*YZ21
!      (PY PY/PY PY)
               W(56) = &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX22+RI(18)*XX22*YYZZ22 &
           +RI(19)*(YY22*YY22+ZZ22*ZZ22) &
           +RI(20)*(XY22*XY22+XZ22*XZ22) &
           +RI(21)*(YY22*ZZ22+ZZ22*YY22) &
           +RI(22)*YZ22*YZ22
!     (PY PY/PZ S)
               W(57)=X(3)*(RI(8)*XX22+RI(9)*YYZZ22) &
                 +RI(10)*(         +Z(3)*XZ22)
!      (PY PY/PZ PX)
               W(58) = &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX31 &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ31 &
           +RI(20)*(XY22*XY31+XZ22*XZ31) &
           +RI(22)*YZ22*YZ31
!      (PY PY/PZ PY)
               W(59) = &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX32 &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ32 &
           +RI(20)*(XY22*XY32+XZ22*XZ32) &
           +RI(22)*YZ22*YZ32
!      (PY PY/PZ PZ)
               W(60) = &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX33 &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ33 &
           +RI(20)*XZ22*XZ33
!     (PZ S/SS)
               W(61)=RI(2)*X(3)
!     (PZ S/PX S)
               W(62)=RI(6)*XX31+RI(7)*ZZ31
!     (PZ S/PX PX)
               W(63)=X(3)*(RI(13)*XX11+RI(14)*YYZZ11) &
                 +RI(15)*(         +Z(3)*XZ11)
!     (PZ S/PY S)
               W(64)=RI(6)*XX32+RI(7)*ZZ32
!     (PZ S/PY PX)
               W(65)=X(3)*(RI(13)*XX21+RI(14)*YYZZ21) &
                 +RI(15)*(         +Z(3)*XZ21)
!     (PZ S/PY PY)
               W(66)=X(3)*(RI(13)*XX22+RI(14)*YYZZ22) &
                 +RI(15)*(         +Z(3)*XZ22)
!     (PZ S/PZ S)
               W(67)=RI(6)*XX33+RI(7)*ZZ33
!     (PZ S/PZ PX)
               W(68)=X(3)*(RI(13)*XX31+RI(14)*ZZ31) &
                 +RI(15)*(         +Z(3)*XZ31)
!     (PZ S/PZ PY)
               W(69)=X(3)*(RI(13)*XX32+RI(14)*ZZ32) &
                 +RI(15)*(         +Z(3)*XZ32)
!     (PZ S/PZ PZ)
               W(70)=X(3)*(RI(13)*XX33+RI(14)*ZZ33) &
                 +RI(15)*(         +Z(3)*XZ33)
!     (PZ PX/S S)
               W(71)=RI(3)*XX31+RI(4)*ZZ31
!     (PZ PX/PX S)
               W(72)=X(1)*(RI(8)*XX31+RI(9)*ZZ31) &
                 +RI(10)*(Y(1)*XY31+Z(1)*XZ31)
!      (PZ PX/PX PX)
               W(73) = &
           (RI(16)*XX31+RI(17)*ZZ31)*XX11+RI(18)*XX31*YYZZ11 &
           +RI(19)*ZZ31*ZZ11 &
           +RI(20)*(XY31*XY11+XZ31*XZ11) &
           +RI(21)*ZZ31*YY11 &
           +RI(22)*YZ31*YZ11
!     (PZ PX/PY S)
               W(74)=X(2)*(RI(8)*XX31+RI(9)*ZZ31) &
                 +RI(10)*(Y(2)*XY31+Z(2)*XZ31)
!      (PZ PX/PY PX)
               W(75) = &
           (RI(16)*XX31+RI(17)*ZZ31)*XX21+RI(18)*XX31*YYZZ21 &
           +RI(19)*ZZ31*ZZ21 &
           +RI(20)*(XY31*XY21+XZ31*XZ21) &
           +RI(21)*ZZ31*YY21 &
           +RI(22)*YZ31*YZ21
!      (PZ PX/PY PY)
               W(76) = &
           (RI(16)*XX31+RI(17)*ZZ31)*XX22+RI(18)*XX31*YYZZ22 &
           +RI(19)*ZZ31*ZZ22 &
           +RI(20)*(XY31*XY22+XZ31*XZ22) &
           +RI(21)*ZZ31*YY22 &
           +RI(22)*YZ31*YZ22
!     (PZ PX/PZ S)
               W(77)=X(3)*(RI(8)*XX31+RI(9)*ZZ31) &
                 +RI(10)*(         +Z(3)*XZ31)
!     (PZ PX/PZ PX)
               W(78) = &
            (RI(16)*XX31+RI(17)*ZZ31)*XX31 &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ31 &
           +RI(20)*(XY31*XY31+XZ31*XZ31) &
           +RI(22)*YZ31*YZ31
!      (PZ PX/PZ PY)
               W(79) = &
            (RI(16)*XX31+RI(17)*ZZ31)*XX32 &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ32 &
           +RI(20)*(XY31*XY32+XZ31*XZ32) &
           +RI(22)*YZ31*YZ32
!      (PZ PX/PZ PZ)
               W(80) = &
            (RI(16)*XX31+RI(17)*ZZ31)*XX33 &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ33 &
           +RI(20)*XZ31*XZ33
!     (PZ PY/S S)
               W(81)=RI(3)*XX32+RI(4)*ZZ32
!     (PZ PY/PX S)
               W(82)=X(1)*(RI(8)*XX32+RI(9)*ZZ32) &
                 +RI(10)*(Y(1)*XY32+Z(1)*XZ32)
!      (PZ PY/PX PX)
               W(83) = &
           (RI(16)*XX32+RI(17)*ZZ32)*XX11+RI(18)*XX32*YYZZ11 &
           +RI(19)*ZZ32*ZZ11 &
           +RI(20)*(XY32*XY11+XZ32*XZ11) &
           +RI(21)*ZZ32*YY11 &
           +RI(22)*YZ32*YZ11
!     (PZ PY/PY S)
               W(84)=X(2)*(RI(8)*XX32+RI(9)*ZZ32) &
                 +RI(10)*(Y(2)*XY32+Z(2)*XZ32)
!      (PZ PY/PY PX)
               W(85) = &
           (RI(16)*XX32+RI(17)*ZZ32)*XX21+RI(18)*XX32*YYZZ21 &
           +RI(19)*ZZ32*ZZ21 &
           +RI(20)*(XY32*XY21+XZ32*XZ21) &
           +RI(21)*ZZ32*YY21 &
           +RI(22)*YZ32*YZ21
!      (PZ PY/PY PY)
               W(86) = &
           (RI(16)*XX32+RI(17)*ZZ32)*XX22+RI(18)*XX32*YYZZ22 &
           +RI(19)*ZZ32*ZZ22 &
           +RI(20)*(XY32*XY22+XZ32*XZ22) &
           +RI(21)*ZZ32*YY22 &
           +RI(22)*YZ32*YZ22
!     (PZ PY/PZ S)
               W(87)=X(3)*(RI(8)*XX32+RI(9)*ZZ32) &
                 +RI(10)*(         +Z(3)*XZ32)
!      (PZ PY/PZ PX)
               W(88) = &
            (RI(16)*XX32+RI(17)*ZZ32)*XX31 &
           +(RI(18)*XX32+RI(19)*ZZ32)*ZZ31 &
           +RI(20)*(XY32*XY31+XZ32*XZ31) &
           +RI(22)*YZ32*YZ31
!      (PZ PY/PZ PY)
               W(89) = &
            (RI(16)*XX32+RI(17)*ZZ32)*XX32 &
           +(RI(18)*XX32+RI(19)*ZZ32)*ZZ32 &
           +RI(20)*(XY32*XY32+XZ32*XZ32) &
           +RI(22)*YZ32*YZ32
!       (PZ PY/PZ PZ)
               W(90) = &
            (RI(16)*XX32+RI(17)*ZZ32)*XX33 &
           +(RI(18)*XX32+RI(19)*ZZ32)*ZZ33 &
           +RI(20)*XZ32*XZ33
!     (PZ PZ/S S)
               W(91)=RI(3)*XX33+RI(4)*ZZ33
!     (PZ PZ/PX S)
               W(92)=X(1)*(RI(8)*XX33+RI(9)*ZZ33) &
                 +RI(10)*(          Z(1)*XZ33)
!       (PZ PZ/PX PX)
               W(93) = &
           (RI(16)*XX33+RI(17)*ZZ33)*XX11+RI(18)*XX33*YYZZ11 &
           +RI(19)*ZZ33*ZZ11 &
           +RI(20)*XZ33*XZ11 &
           +RI(21)*ZZ33*YY11
!     (PZ PZ/PY S)
               W(94)=X(2)*(RI(8)*XX33+RI(9)*ZZ33) &
                 +RI(10)*(         +Z(2)*XZ33)
!       (PZ PZ/PY PX)
               W(95) = &
           (RI(16)*XX33+RI(17)*ZZ33)*XX21+RI(18)*XX33*YYZZ21 &
           +RI(19)*ZZ33*ZZ21 &
           +RI(20)*XZ33*XZ21 &
           +RI(21)*ZZ33*YY21
!       (PZ PZ/PY PY)
               W(96) = &
           (RI(16)*XX33+RI(17)*ZZ33)*XX22+RI(18)*XX33*YYZZ22 &
           +RI(19)*ZZ33*ZZ22 &
           +RI(20)*XZ33*XZ22 &
           +RI(21)*ZZ33*YY22
!     (PZ PZ/PZ S)
               W(97)=X(3)*(RI(8)*XX33+RI(9)*ZZ33) &
                 +RI(10)*(         +Z(3)*XZ33)
!       (PZ PZ/PZ PX)
               W(98) = &
            (RI(16)*XX33+RI(17)*ZZ33)*XX31 &
           +(RI(18)*XX33+RI(19)*ZZ33)*ZZ31 &
           +RI(20)*XZ33*XZ31
!       (PZ PZ/PZ PY)
               W(99) = &
            (RI(16)*XX33+RI(17)*ZZ33)*XX32 &
           +(RI(18)*XX33+RI(19)*ZZ33)*ZZ32 &
           +RI(20)*XZ33*XZ32
!       (PZ PZ/PZ PZ)
               W(100) = &
            (RI(16)*XX33+RI(17)*ZZ33)*XX33 &
           +(RI(18)*XX33+RI(19)*ZZ33)*ZZ33 &
           +RI(20)*XZ33*XZ33
            ELSE
!     (PX S/S S)
               W(2)=RI(2)*X(1)
!     (PX PX/S S)
               W(3)=RI(3)*XX11+RI(4)*YYZZ11
!     (PY S/S S)
               W(4)=RI(2)*X(2)
!     (PY PX/S S)
               W(5)=RI(3)*XX21+RI(4)*YYZZ21
!     (PY PY/S S)
               W(6)=RI(3)*XX22+RI(4)*YYZZ22
!     (PZ S/SS)
               W(7)=RI(2)*X(3)
!     (PZ PX/S S)
               W(8)=RI(3)*XX31+RI(4)*ZZ31
!     (PZ PY/S S)
               W(9)=RI(3)*XX32+RI(4)*ZZ32
!     (PZ PZ/S S)
               W(10)=RI(3)*XX33+RI(4)*ZZ33
            END IF
         END IF
 
         ! Now rotate the nuclear attraction integrals.
         ! The storage of the nuclear attraction integrals  core(kl/ij) is
         ! (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
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

         rij = rij
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
         zz = REAL(sepi%zeff,dp)*REAL(sepj%zeff,dp)
         enuc = zz*gam
         scale=ABS(scale*enuc)
         zz=zz/rij
         IF(itype == 2 .OR. itype == 3) THEN
           ! AM1 or PM3 correction to nuclear repulsion
            DO ig=1,SIZE(sepi%fn1)
               IF(ABS(sepi%fn1(ig)) > 0._dp) THEN
                  AX = sepi%fn2(ig)*(rij-sepi%fn3(ig))**2
                  IF(ax <= 25._dp) THEN
                     scale=scale +zz*sepi%fn1(ig)*EXP(-ax)
                  ENDIF
               ENDIF
               IF(ABS(sepj%fn1(ig)) > 0._dp) THEN
                  ax = sepj%fn2(ig)*(rij-sepj%fn3(ig))**2
                  IF(ax <= 25._dp) THEN
                     scale=scale +zz*sepj%fn1(ig)*EXP(-ax)
                  ENDIF
               ENDIF
            END DO
         ENDIF
         enuc=enuc+scale
 
      ENDIF

  END SUBROUTINE rotint

! *****************************************************************************

  SUBROUTINE terep ( sepi, sepj, rij, ri, core )
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
  !            core    = 4 X 2 array of electron-core attraction
  !                        integrals
  !     The two-centre repulsion integrals (over local coordinates) are
  !     stored as follows (where p-sigma = O,  and p-pi = P and P* )
  !     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
  !     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
  !     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
  !     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
  !     (PP/P*P*)=21,   (P*P/P*P)=22.
  !     The storage of the nuclear attraction integrals  core(kl/ij) iS
  !     (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
  !     where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
  !
! *****************************************************************************

      TYPE ( semi_empirical_type ), INTENT ( IN )      :: sepi, sepj
      REAL ( dp ), INTENT ( IN )                       :: rij
      REAL ( dp ), DIMENSION (:), INTENT ( OUT )       :: ri
      REAL ( dp ), DIMENSION (:,:),INTENT ( OUT )      :: core

      LOGICAL        :: si,sj
      REAL ( dp )    :: r
      REAL ( dp )    :: ami,amj,adi,adj,aqi,aqj
      REAL ( dp )    :: zi,zj,ddi,ddj,qqi,qqj
      REAL ( dp )    :: aee,da,qa,ade,aqe,rsq,ee,db,qb,aed,aeq,&
                        axx,adq,aqd,aqq,xxx,yyy,zzz,www,&
                        dze,qzze,qxxe,edz,eqzz,eqxx,dxdx,dzdz,dzqxx,&
                        qxxdz,dzqzz,qzzdz,qxxqxx,qxxqyy,qxxqzz,qzzqxx,&
                        qzzqzz,dxqxz,qxzdx,qxzqxz
      REAL ( dp ), DIMENSION ( 72 ) :: arg,sqr

      ri = 0._dp
      core = 0._dp
      r=rij
      si = (sepi%natorb >= 3)
      sj = (sepj%natorb >= 3)
      zi = REAL ( sepi%zeff, dp )
      zj = REAL ( sepj%zeff, dp )
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

  END SUBROUTINE terep

! *****************************************************************************

END MODULE semi_empirical_integrals
