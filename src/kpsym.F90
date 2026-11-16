!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Crystal-symmetry routines originating from the K290/ACMI code.
!>
!>     GROUP1, PGL1, ATFTM1, and ROT1 are based on the Worlton-Warren
!>     Computer Physics Communications package (1971, 1974).
! **************************************************************************************************
MODULE kpsym

   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
   USE string_utilities,                ONLY: xstring
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC  :: GROUP1s

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpsym'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param iout ...
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \param nat ...
!> \param ty ...
!> \param x ...
!> \param b1 ...
!> \param b2 ...
!> \param b3 ...
!> \param ihg ...
!> \param ihc ...
!> \param isy ...
!> \param li ...
!> \param nc ...
!> \param indpg ...
!> \param ib ...
!> \param ntvec ...
!> \param v ...
!> \param f0 ...
!> \param r ...
!> \param tvec ...
!> \param origin ...
!> \param rx ...
!> \param isc ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE group1s(iout, a1, a2, a3, nat, ty, x, b1, b2, b3, &
                      ihg, ihc, isy, li, nc, indpg, ib, ntvec, &
                      v, f0, r, tvec, origin, rx, isc, delta)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN ON SEPTEMBER 10TH - FROM THE ACMI COMPLEX            ==
      ! == (WORLTON AND WARREN, COMPUT.PHYS.COMMUN. 8,71-84 (1974))     ==
      ! == (AND 3,88-117 (1972))                                        ==
      ! == BASIC CRYSTALLOGRAPHIC INFORMATION                           ==
      ! == ABOUT A GIVEN CRYSTAL STRUCTURE.                             ==
      ! == SUBROUTINES NEEDED: PGL1,ATFTM1,ROT1,RLV3                    ==
      ! ==--------------------------------------------------------------==
      ! == INPUT DATA:                                                  ==
      ! == IOUT ... NUMBER OF THE OUTPUT UNIT FOR ON-LINE PRINTING      ==
      ! ==          OF VARIOUS MESSAGES                                 ==
      ! ==          IF IOUT<=0 NO MESSAGE                             ==
      ! == A1,A2,A3 .. ELEMENTARY TRANSLATIONS OF THE LATTICE, IN SOME  ==
      ! ==          UNIT OF LENGTH                                      ==
      ! == NAT .... NUMBER OF ATOMS IN THE UNIT CELL                    ==
      ! ==          ALL THE DIMENSIONS ARE SET FOR NAT <= 20          ==
      ! == TY ..... INTEGERS DISTINGUISHING BETWEEN THE ATOMS OF        ==
      ! ==          DIFFERENT TYPE. TY(I) IS THE TYPE OF THE I-TH ATOM  ==
      ! ==          OF THE BASIS                                        ==
      ! == X ...... CARTESIAN COORDINATES OF THE NAT ATOMS OF THE BASIS ==
      ! == DELTA... REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)           ==
      ! ==--------------------------------------------------------------==
      ! == OUTPUT DATA:                                                 ==
      ! == B1,B2,B3 .. RECIPROCAL LATTICE VECTORS, NOT MULTIPLIED BY    ==
      ! ==          ANY 2PI, IN UNITS RECIPROCAL TO THOSE OF A1,A2,A3   ==
      ! == IHG .... POINT GROUP OF THE PRIMITIVE LATTICE, HOLOHEDRAL    ==
      ! ==          GROUP NUMBER:                                       ==
      ! ==          IHG=1 STANDS FOR TRICLINIC    SYSTEM                ==
      ! ==          IHG=2 STANDS FOR MONOCLINIC   SYSTEM                ==
      ! ==          IHG=3 STANDS FOR ORTHORHOMBIC SYSTEM                ==
      ! ==          IHG=4 STANDS FOR TETRAGONAL   SYSTEM                ==
      ! ==          IHG=5 STANDS FOR CUBIC        SYSTEM                ==
      ! ==          IHG=6 STANDS FOR TRIGONAL     SYSTEM                ==
      ! ==          IHG=7 STANDS FOR HEXAGONAL    SYSTEM                ==
      ! == IHC .... CODE DISTINGUISHING BETWEEN HEXAGONAL AND CUBIC     ==
      ! ==          GROUPS                                              ==
      ! ==          IHC=0 STANDS FOR HEXAGONAL GROUPS                   ==
      ! ==          IHC=1 STANDS FOR CUBIC GROUPS                       ==
      ! == ISY .... CODE INDICATING WHETHER THE SPACE GROUP IS          ==
      ! ==          SYMMORPHIC OR NONSYMMORPHIC                         ==
      ! ==          ISY= 0 NONSYMMORPHIC GROUP                          ==
      ! ==          ISY= 1 SYMMORPHIC GROUP                             ==
      ! ==          ISY=-1 SYMMORPHIC GROUP WITH NON-STANDARD ORIGIN    ==
      ! ==          ISY=-2 UNDETERMINED (NORMALLY NEVER)                ==
      ! ==          THE GROUP IS CONSIDERED SYMMORPHIC IF FOR EACH      ==
      ! ==          OPERATION OF THE POINT GROUP THE SUM OF THE 3       ==
      ! ==          COMPONENTS OF ABS(V(N)) (NONPRIMITIVE TRANSLATION,  ==
      ! ==          SEE BELOW) IS LT. 0.0001                            ==
      ! == ORIGIN   STANDARD ORIGIN IF SYMMORPHIC (CRYSTAL COORDINATES) ==
      ! == LI ..... CODE INDICATING WHETHER THE POINT GROUP             ==
      ! ==          OF THE CRYSTAL CONTAINS INVERSION OR NOT            ==
      ! ==          (OPERATIONS 13 OR 25 IN RESPECTIVELY HEXAGONAL      ==
      ! ==          OR CUBIC GROUPS).                                   ==
      ! ==          LI=0 MEANS: DOES NOT CONTAIN INVERSION              ==
      ! ==          LI>0 MEANS: THERE IS INVERSION IN THE POINT      ==
      ! ==                         GROUP OF THE CRYSTAL                 ==
      ! == NC ..... TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP OF THE  ==
      ! ==          CRYSTAL                                             ==
      ! == INDPG .. POINT GROUP INDEX (DETERMINED IF SYMMORPHIC GROUP)  ==
      ! == IB ..... LIST OF THE ROTATIONS CONSTITUTING THE POINT GROUP  ==
      ! ==          OF THE CRYSTAL. THE NUMBERING IS THAT DEFINED IN    ==
      ! ==          WORLTON AND WARREN, I.E. THE ONE MATERIALIZED IN THE==
      ! ==          ARRAY R (SEE BELOW)                                 ==
      ! ==          ONLY THE FIRST NC ELEMENTS OF THE ARRAY IB ARE      ==
      ! ==          MEANINGFUL                                          ==
      ! == NTVEC .. NUMBER OF TRANSLATIONAL VECTORS                     ==
      ! ==          ASSOCIATED WITH IDENTITY OPERATOR I.E.              ==
      ! ==          GIVES THE NUMBER OF IDENTICAL PRIMITIVE CELLS       ==
      ! == V ...... NONPRIMITIVE TRANSLATIONS (IN THE CASE OF NONSYMMOR-==
      ! ==          PHIC GROUPS). V(I,N) IS THE I-TH COMPONENT          ==
      ! ==          OF THE TRANSLATION CONNECTED WITH THE N-TH ELEMENT  ==
      ! ==          OF THE POINT GROUP (I.E. WITH THE ROTATION          ==
      ! ==          NUMBER IB(N) ).                                     ==
      ! ==          ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS,       ==
      ! ==          THEY REFER TO THE SYSTEM A1,A2,A3.                  ==
      ! == F0 ..... THE FUNCTION DEFINED IN MARADUDIN, IPATOVA BY       ==
      ! ==          EQ. (3.2.12): ATOM TRANSFORMATION TABLE.            ==
      ! ==          THE ELEMENT F0(N,KAPA) MEANS THAT THE N-TH          ==
      ! ==          OPERATION OF THE SPACE GROUP (I.E. OPERATION NUMBER ==
      ! ==          IB(N), TOGETHER WITH AN EVENTUAL NONPRIMITIVE       ==
      ! ==          TRANSLATION  V(N)) TRANSFERS THE ATOM KAPA INTO THE ==
      ! ==          ATOM F0(N,KAPA).                                    ==
      ! ==          THE 49TH LINE GIVES EQUIVALENT ATOMS FOR            ==
      ! ==          FRACTIONAl TRANSLATIONS ASSOCIATED WITH IDENTITY    ==
      ! == R ...... LIST OF THE 3 X 3 ROTATION MATRICES                 ==
      ! ==          (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)    ==
      ! ==          ALL 48 OR 24 MATRICES ARE LISTED.                   ==
      ! ==          FOLLOW NOTATION OF WORLTON-WARREN(1972)             ==
      ! == TVEC  .. LIST OF NTVEC TRANSLATIONAL VECTORS                 ==
      ! ==          ASSOCIATED WITH IDENTITY OPERATOR                   ==
      ! ==          TVEC(1:3,1) = \(0,0,0\)                             ==
      ! ==          (CRYSTAL COORDINATES)                               ==
      ! == RX ..... SCRATCH ARRAY                                       ==
      ! == ISC .... SCRATCH ARRAY                                       ==
      ! ==--------------------------------------------------------------==
      ! == PRINTED OUTPUT:                                              ==
      ! == PROGRAM PRINTS THE TYPE OF THE LATTICE (IHG, IN WORDS),      ==
      ! == LISTS THE OPERATIONS OF THE  POINT GROUP OF THE              ==
      ! == CRYSTAL, INDICATES WHETHER THE SPACE GROUP IS SYMMORPHIC OR  ==
      ! == NONSYMMORPHIC AND WHETHER THE POINT GROUP OF THE CRYSTAL     ==
      ! == CONTAINS INVERSION.                                          ==
      ! ==--------------------------------------------------------------==
      INTEGER                                            :: iout
      REAL(dp)                                           :: a1(3), a2(3), a3(3)
      INTEGER                                            :: nat, ty(nat)
      REAL(dp)                                           :: x(3, nat), b1(3), b2(3), b3(3)
      INTEGER                                            :: ihg, ihc, isy, li, nc, indpg, ib(48), &
                                                            ntvec
      REAL(dp)                                           :: v(3, 48)
      INTEGER                                            :: f0(49, nat)
      REAL(dp)                                           :: r(3, 3, 48), tvec(3, nat), origin(3), &
                                                            rx(3, nat)
      INTEGER                                            :: isc(nat)
      REAL(dp)                                           :: delta

      INTEGER                                            :: i, ncprim
      REAL(dp)                                           :: a(3, 3), ai(3, 3), ap(3, 3), api(3, 3)

      DO i = 1, 3
         a(i, 1) = a1(i)
         a(i, 2) = a2(i)
         a(i, 3) = a3(i)
      END DO
      ! ==--------------------------------------------------------------==
      ! == A(I,J) IS THE I-TH CARTESIAN COMPONENT OF THE J-TH PRIMITIVE ==
      ! == TRANSLATION VECTOR OF THE DIRECT LATTICE                     ==
      ! == TY(I) IS AN INTEGER DISTINGUISHING ATOMS OF DIFFERENT TYPE,  ==
      ! == I.E., DIFFERENT ATOMIC SPECIES                               ==
      ! == X(J,I) IS THE J-TH CARTESIAN COMPONENT OF THE POSITION       ==
      ! == VECTOR FOR THE I-TH ATOM IN THE UNIT CELL.                   ==
      ! ==--------------------------------------------------------------==
      ! ==DETERMINE PRIMITIVE LATTICE VECTORS FOR THE RECIPROCAL LATTICE==
      ! ==--------------------------------------------------------------==
      CALL calbrec(a, ai)
      DO i = 1, 3
         b1(i) = ai(1, i)
         b2(i) = ai(2, i)
         b3(i) = ai(3, i)
      END DO
      ! ==--------------------------------------------------------------==
      ! Determination of the translation vectors associated with
      ! the Identity matrix i.e. if the cell is duplicated
      ! Give also the ``primitive lattice''
      CALL primlatt(a, ai, ap, api, nat, ty, x, ntvec, tvec, f0, isc, delta)
      ! ==--------------------------------------------------------------==
      ! Determination of the holohedral group (and crystal system)
      CALL pgl1(ap, api, ihc, nc, ib, ihg, r, delta)
      IF (ntvec > 1) THEN
         ! All rotations found by PGL1 have axes in x, y or z cart. axis
         ! So we have too check if we do not loose symmetry
         ncprim = nc
         ! The hexagonal system is found if the z axis is the sixfold axis
         CALL pgl1(a, ai, ihc, nc, ib, ihg, r, delta)
         IF (ncprim > nc) THEN
            ! More symmetry with
            CALL pgl1(ap, api, ihc, nc, ib, ihg, r, delta)
         END IF
      END IF

      ! Determination of the space group
      CALL atftm1(iout, r, v, x, f0, origin, ib, ty, nat, ihg, ihc, rx, &
                  nc, indpg, ntvec, a, ai, li, isy, isc, delta)

      IF (iout > 0) THEN
         IF (li > 0) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(1X,A)') &
                  'KPSYM| THE POINT GROUP OF THE CRYSTAL CONTAINS THE INVERSION'
            END IF
         END IF
         IF (iout > 0) THEN
            WRITE (iout, *)
         END IF
      END IF

   END SUBROUTINE group1s
! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param ai ...
! **************************************************************************************************
   SUBROUTINE calbrec(a, ai)
      ! ==--------------------------------------------------------------==
      ! == CALCULATE RECIPROCAL VECTOR BASIS (AI(1:3,1:3))              ==
      ! == INPUT:                                                       ==
      ! ==   A(3,3) A(I,J) IS THE I-TH CARTESIAN COMPONENT              ==
      ! ==          OF THE J-TH PRIMITIVE TRANSLATION VECTOR OF         ==
      ! ==          THE DIRECT LATTICE                                  ==
      ! == OUTPUT:                                                      ==
      ! ==   AI(3,3) RECIPROCAL VECTOR BASIS                            ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: a(3, 3), ai(3, 3)

      INTEGER                                            :: i, il, iu, j, jl, ju
      REAL(dp)                                           :: det

      det = a(1, 1)*a(2, 2)*a(3, 3) + a(2, 1)*a(1, 3)*a(3, 2) + &
            a(3, 1)*a(1, 2)*a(2, 3) - a(1, 1)*a(2, 3)*a(3, 2) - &
            A(2, 1)*A(1, 2)*A(3, 3) - A(3, 1)*A(1, 3)*A(2, 2)
      det = 1._dp/det
      DO i = 1, 3
         il = 1
         iu = 3
         IF (i == 1) il = 2
         IF (i == 3) iu = 2
         DO j = 1, 3
            jl = 1
            ju = 3
            IF (j == 1) jl = 2
            IF (j == 3) ju = 2
            ai(j, i) = (-1._dp)**(i + j)*det* &
                       (A(IL, JL)*A(IU, JU) - A(IL, JU)*A(IU, JL))
         END DO
      END DO
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE calbrec
   ! ==================================================================
! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param ai ...
!> \param ap ...
!> \param api ...
!> \param nat ...
!> \param ty ...
!> \param x ...
!> \param ntvec ...
!> \param tvec ...
!> \param f0 ...
!> \param isc ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE primlatt(a, ai, ap, api, nat, ty, x, ntvec, tvec, f0, isc, delta)
      ! ==--------------------------------------------------------------==
      ! == DETERMINATION OF THE TRANSLATION VECTORS ASSOCIATED WITH     ==
      ! == THE IDENTITY SYMMETRY I.E. IF THE CELL IS DUPLICATED         ==
      ! == GIVE ALSO THE PRIMITIVE DIRECT AND RECIPROCAL LATTICE VECTOR ==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==   A(3,3)   A(I,J) IS THE I-TH CARTESIAN COMPONENT            ==
      ! ==            OF THE J-TH TRANSLATION VECTOR OF                 ==
      ! ==            THE DIRECT LATTICE                                ==
      ! ==   AI(3,3)  RECIPROCAL VECTOR BASIS (CARTESIAN)               ==
      ! ==   NAT      NUMBER OF ATOMS                                   ==
      ! ==   TY(NAT)  TYPE OF ATOMS                                     ==
      ! ==   X(3,NAT) ATOMIC COORDINATES IN CARTESIAN COORDINATES       ==
      ! == DELTA      REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)         ==
      ! == OUTPUT:                                                      ==
      ! ==   AP(3,3)  COMPONENTS OF THE PRIMITIVE TRANSLATION VECTORS   ==
      ! ==   API(3,3) PRIMITIVE RECIPROCAL BASIS VECTORS                ==
      ! ==            BOTH BAISI ARE IN CARTESIAN COORDINATES           ==
      ! ==   NTVEC    NUMBER OF TRANSLATION VECTORS (FRACTIONNAL)       ==
      ! ==   TVEC(3,NTVEC) COMPONENTS OF TRANSLATIONAL VECTORS          ==
      ! ==                 (CRYSTAL COORDINATES)                        ==
      ! ==   F0(49,NAT)    GIVES INEQUIVALENT ATOM FOR EACH ATOM        ==
      ! ==                 THE 49-TH LINE                               ==
      ! ==   ISC(NAT)      SCRATCH ARRAY                                ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: a(3, 3), ai(3, 3), ap(3, 3), api(3, 3)
      INTEGER                                            :: nat, ty(nat)
      REAL(dp)                                           :: x(3, nat)
      INTEGER                                            :: ntvec
      REAL(dp)                                           :: tvec(3, nat)
      INTEGER                                            :: f0(49, nat), isc(nat)
      REAL(dp)                                           :: delta

      INTEGER                                            :: i, il, iv, j, k2
      LOGICAL                                            :: oksym
      REAL(dp)                                           :: vr(3), xb(3)

! Variables
! ==--------------------------------------------------------------==
! First we check if there exist fractional translational vectors
! associated with Identity operation i.e.
! if the cell is duplicated or not.

      ntvec = 1
      tvec(1, 1) = 0._dp
      tvec(2, 1) = 0._dp
      tvec(3, 1) = 0._dp
      DO i = 1, nat
         f0(49, i) = i
      END DO
      DO k2 = 2, nat
         IF (ty(1) /= ty(k2)) CYCLE
         DO i = 1, 3
            xb(i) = x(i, k2) - x(i, 1)
         END DO
         ! A fractional translation vector VR is defined.
         CALL rlv3(ai, xb, vr, il, delta)
         CALL checkrlv3(1, nat, ty, x, x, vr, f0, ai, isc, .TRUE., oksym, delta)
         IF (oksym) THEN
            ! A fractional translational vector is found
            ntvec = ntvec + 1
            ! F0(49,1:NAT) gives number of equivalent atoms
            ! and has atom indexes of inequivalent atoms (for translation)
            DO i = 1, nat
               IF (f0(49, i) > f0(1, i)) f0(49, i) = f0(1, i)
            END DO
            DO i = 1, 3
               tvec(i, ntvec) = vr(i)
            END DO
         END IF
      END DO
      ! ==-------------------------------------------------------------==
      DO i = 1, 3
         ap(1, i) = a(1, i)
         ap(2, i) = a(2, i)
         ap(3, i) = a(3, i)
         api(1, i) = ai(1, i)
         api(2, i) = ai(2, i)
         api(3, i) = ai(3, i)
      END DO
      IF (ntvec == 1) THEN
         ! The current cell is definitely a primitive one
         ! Copy A and AI to AP and API
      ELSE
         ! We are looking for the primitive lattice vector basis set
         ! AP is our current lattice vector basis
         DO iv = 2, ntvec
            ! TVEC in cartesian coordinates
            DO i = 1, 3
               xb(i) = tvec(1, iv)*a(i, 1) &
                       + TVEC(2, IV)*A(I, 2) &
                       + TVEC(3, IV)*A(I, 3)
            END DO
            ! We calculare TVEC in AP basis
            CALL rlv3(api, xb, vr, il, delta)
            DO i = 1, 3
               IF (ABS(vr(i)) > delta) THEN
                  il = NINT(1._dp/ABS(vr(i)))
                  IF (il > 1) THEN
                     ! We replace AP(1:3,I) by TVEC(1:3,IV)
                     DO j = 1, 3
                        ap(j, i) = xb(j)
                     END DO
                     ! Calculate new API
                     CALL calbrec(ap, api)
                     EXIT
                  END IF
               END IF
            END DO
         END DO
      END IF
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE primlatt
   ! ==================================================================
! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param ai ...
!> \param ihc ...
!> \param nc ...
!> \param ib ...
!> \param ihg ...
!> \param r ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE pgl1(a, ai, ihc, nc, ib, ihg, r, delta)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX          ==
      ! == AUXILIARY SUBROUTINE TO GROUP1                               ==
      ! == SUBROUTINE PGL DETERMINES THE POINT GROUP OF THE LATTICE     ==
      ! == AND THE CRYSTAL SYSTEM.                                      ==
      ! == SUBROUTINES NEEDED: ROT1, RLV3                               ==
      ! ==--------------------------------------------------------------==
      ! == WARNING: FOR THE HEXAGONAL SYSTEM, THE 3RD AXIS SUPPOSE      ==
      ! ==          TO BE THE SIX-FOLD AXIS                             ==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! == A ..... DIRECT LATTICE VECTORS                               ==
      ! == AI .... RECIPROCAL LATTICE VECTORS                           ==
      ! == DELTA.. REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)            ==
      ! ==--------------------------------------------------------------==
      ! == OUTPUT:                                                      ==
      ! == IHC .... CODE DISTINGUISHING BETWEEN HEXAGONAL AND CUBIC     ==
      ! ==          GROUPS                                              ==
      ! ==          IHC=0 STANDS FOR HEXAGONAL GROUPS                   ==
      ! ==          IHC=1 STANDS FOR CUBIC GROUPS                       ==
      ! == NC .... NUMBER OF ROTATIONS IN THE POINT GROUP               ==
      ! == IB .... SET OF ROTATION                                      ==
      ! == IHG .... POINT GROUP OF THE PRIMITIVE LATTICE, HOLOHEDRAL    ==
      ! ==          GROUP NUMBER:                                       ==
      ! ==          IHG=1 STANDS FOR TRICLINIC    SYSTEM                ==
      ! ==          IHG=2 STANDS FOR MONOCLINIC   SYSTEM                ==
      ! ==          IHG=3 STANDS FOR ORTHORHOMBIC SYSTEM                ==
      ! ==          IHG=4 STANDS FOR TETRAGONAL   SYSTEM                ==
      ! ==          IHG=5 STANDS FOR CUBIC        SYSTEM                ==
      ! ==          IHG=6 STANDS FOR TRIGONAL     SYSTEM                ==
      ! ==          IHG=7 STANDS FOR HEXAGONAL    SYSTEM                ==
      ! == R ...... LIST OF THE 3 X 3 ROTATION MATRICES                 ==
      ! ==          (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)    ==
      ! ==          ALL 48 OR 24 MATRICES ARE LISTED.                   ==
      ! ==          FOLLOW NOTATION OF WORLTON-WARREN(1972)             ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: a(3, 3), ai(3, 3)
      INTEGER                                            :: ihc, nc, ib(48), ihg
      REAL(dp)                                           :: r(3, 3, 48), delta

      INTEGER                                            :: i, j, k, lx, n, nr
      REAL(dp)                                           :: tr, vr(3), xa(3)

      DO ihc = 0, 1
         ! IHC is 0 for hexagonal groups and 1 for cubic groups.
         IF (ihc == 0) THEN
            nr = 24
         ELSE
            nr = 48
         END IF
         nc = 0
         ! Constructs rotation operations.
         CALL rot1(ihc, r)
         loop_rotation: DO n = 1, nr
            ib(n) = 0
            ! Rotate the A1,2,3 vectors by rotation No. N
            DO k = 1, 3
               DO i = 1, 3
                  xa(i) = 0._dp
                  DO j = 1, 3
                     xa(i) = xa(i) + r(i, j, n)*a(j, k)
                  END DO
               END DO
               CALL rlv3(ai, xa, vr, lx, delta)
               tr = 0._dp
               DO i = 1, 3
                  tr = tr + ABS(vr(i))
               END DO
               ! If VR.ne.0, then XA cannot be a multiple of a lattice vector
               IF (tr > delta) CYCLE loop_rotation
            END DO
            nc = nc + 1
            ib(nc) = n
         END DO loop_rotation
         ! ==------------------------------------------------------------==
         ! IHG stands for holohedral group number.
         IF (ihc == 0) THEN
            ! Hexagonal group:
            IF (nc == 12) ihg = 6
            IF (nc > 12) ihg = 7
            IF (nc >= 12) RETURN
            ! Too few operations, try cubic group: (IHC=1,NR=48)
         ELSE
            ! Cubic group:
            IF (nc < 4) ihg = 1
            IF (nc == 4) ihg = 2
            IF (nc > 4) ihg = 3
            IF (nc == 16) ihg = 4
            IF (nc > 16) ihg = 5
            RETURN
         END IF
      END DO
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE pgl1
   ! ==================================================================
! **************************************************************************************************
!> \brief ...
!> \param ai ...
!> \param xb ...
!> \param vr ...
!> \param il ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE rlv3(ai, xb, vr, il, delta)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX          ==
      ! == AUXILIARY SUBROUTINE TO GROUP1                               ==
      ! == SUBROUTINE RLV REMOVES A DIRECT LATTICE VECTOR               ==
      ! == FROM XB LEAVING THE REMAINDER IN VR.                         ==
      ! == IF A NONZERO LATTICE VECTOR WAS REMOVED, IL IS MADE NONZERO. ==
      ! == VR STANDS FOR V-REFERENCE.                                   ==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==    AI(I,J) ARE THE RECIPROCAL LATTICE VECTORS,               ==
      ! ==            B(I) = AI(I,J),J=1,2,3                            ==
      ! ==    XB(1:3) VECTOR IN CARTESIAN COORDINATES                   ==
      ! == DELTA      REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)         ==
      ! == OUTPUT:                                                      ==
      ! ==    VR IS NOT GIVEN IN CARTESIAN COORDINATES BUT              ==
      ! ==       IN THE SYSTEM A1,A2,A3 (CRYSTAL COORDINATES)           ==
      ! ==       AND BETWEEN -1/2 AND 1/2                               ==
      ! ==    IL ABS OF VR                                              ==
      ! == K.K., 23.10.1979                                             ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: ai(3, 3), xb(3), vr(3)
      INTEGER                                            :: il
      REAL(dp)                                           :: delta

      INTEGER                                            :: i
      REAL(dp)                                           :: ts

      il = 0
      DO i = 1, 3
         vr(i) = 0._dp
      END DO
      ts = ABS(xb(1)) + ABS(xb(2)) + ABS(xb(3))
      IF (ts <= delta) RETURN
      DO i = 1, 3
         vr(i) = vr(i) + ai(i, 1)*xb(1) + ai(i, 2)*xb(2) + ai(i, 3)*xb(3)
         il = il + NINT(ABS(vr(i)))
         ! Change in order to have correct determination of origin and
         ! symmorphic group (T.D 30/03/98)
         ! VR(I) = - MOD(real(VR(I),kind=dp),1._dp)
         vr(i) = NINT(vr(i)) - vr(i)
      END DO
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE rlv3
   ! ==================================================================
! **************************************************************************************************
!> \brief ...
!> \param iout ...
!> \param r ...
!> \param v ...
!> \param x ...
!> \param f0 ...
!> \param origin ...
!> \param ib ...
!> \param ty ...
!> \param nat ...
!> \param ihg ...
!> \param ihc ...
!> \param rx ...
!> \param nc ...
!> \param indpg ...
!> \param ntvec ...
!> \param a ...
!> \param ai ...
!> \param li ...
!> \param isy ...
!> \param isc ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE atftm1(iout, r, v, x, f0, origin, ib, ty, nat, ihg, ihc, &
                     rx, nc, indpg, ntvec, a, ai, li, isy, isc, delta)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX          ==
      ! == AUXILIARY SUBROUTINE TO GROUP1                               ==
      ! == SUBROUTINE ATFTMT DETERMINES                                 ==
      ! ==            THE POINT GROUP OF THE CRYSTAL,                   ==
      ! ==            THE ATOM TRANSFORMATION TABLE,F0,                 ==
      ! ==            THE FRACTIONAL TRANSLATIONS,V,                    ==
      ! ==            ASSOCIATED WITH EACH ROTATION.                    ==
      ! == SUBROUTINES NEEDED: RLV3 CHECKRLV3 SYMMORPHIC XSTRING        ==
      ! == MAY 14TH,1998: A LOT OF CHANGES (ARGUMENTS)                  ==
      ! ==                BETTER DETERMINATION OF V                     ==
      ! == SEP 15TH,1998: DETERMINATION OF FRACTIONAL TRANSLATIONAL VEC.==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==   IOUT Logical file number (output)                          ==
      ! ==        If IOUT<=0 no message                                 ==
      ! ==   IHG  Holohedral group number (determined by PGL1)          ==
      ! ==   IHC  Code distinguishing between hexagonal and cubic groups==
      ! ==          IHC=0 stands for hexagonal groups                   ==
      ! ==          IHC=1 stands for cubic groups                       ==
      ! ==   NC   Number of rotation operations                         ==
      ! ==   NAT  Number of atoms (used in the routine)                 ==
      ! ==   X    Coordinates of atoms (cartesian)                      ==
      ! ==   TY   Type of atoms                                         ==
      ! ==   R    Sets of transformation operations (cartesian)         ==
      ! ==   IB   Index giving NC operations in R                       ==
      ! ==   AI   Reciprocal lattice vectors                            ==
      ! ==   NTVEC Number of translational vectors                      ==
      ! ==        associated with Identity                              ==
      ! ==        if primitive cell NTVEC=1, TVEC=(0,0,0)               ==
      ! == DELTA  REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)          ==
      ! == OUTPUT:                                                      ==
      ! ==   RX(3,NAT) Scratch array                                    ==
      ! ==   ISC(NAT)  Scratch array                                    ==
      ! ==   NC        is modified (number of symmetry operations)      ==
      ! ==   INDPG     Point group index                                ==
      ! ==   V(3,48)   The fractional translations associated           ==
      ! ==             with each rotation (crystal coordinates)         ==
      ! ==   F0(1:48,NAT)                                               ==
      ! ==        The atom transformation table for rotation (48,NAT)   ==
      ! ==   ORIGIN Standard origin if symmorphic (crystal coordinates) ==
      ! ==   ISY  = 1 Isommorphic group                                 ==
      ! ==        =-1 Isommorphic group with non-standard origin        ==
      ! ==        = 0 Non-Isommorphic group                             ==
      ! ==        =-2 Undetermined (normally never)                     ==
      ! == LI ..... Code indicating whether the point group             ==
      ! ==          of the crystal contains inversion or not            ==
      ! ==          (operations 13 or 25 in respectively hexagonal      ==
      ! ==          or cubic groups).                                   ==
      ! ==          LI=0    : does not contain inversion                ==
      ! ==          LI>0 : there is inversion in the point              ==
      ! ==                    group of the crystal                      ==
      ! ==--------------------------------------------------------------==
      ! INDPG group   indpg   group    indpg  group     indpg   group   ==
      ! == 1    1 (c1)    9    3m (c3v)   17  4/mmm(d4h)  25   222(d2)  ==
      ! == 2   <1>(ci)   10   <3>m(d3d)   18     6 (c6)   26   mm2(c2v) ==
      ! == 3    2 (c2)   11     4 (c4)    19    <6>(c3h)  27   mmm(d2h) ==
      ! == 4    m (c1h)  12    <4>(s4)    20    6/m(c6h)  28   23 (t)   ==
      ! == 5   2/m(c2h)  13    4/m(c4h)   21    622(d6)   29   m3 (th)  ==
      ! == 6    3 (c3)   14    422(d4)    22    6mm(c6v)  30   432(o)   ==
      ! == 7   <3>(c3i)  15    4mm(c4v)   23  <6>m2(d3h)  31 <4>3m(td)  ==
      ! == 8   32 (d3)   16  <4>2m(d2d)   24  6/mmm(d6h)  32   m3m(oh)  ==
      ! ==--------------------------------------------------------------==
      ! rname_cubic: Name of 48 rotations (convention Warren-Worlton)
      INTEGER                                            :: iout
      REAL(dp)                                           :: r(3, 3, 48), v(3, 48), origin(3)
      INTEGER                                            :: ib(48), nat, ty(nat), f0(49, nat)
      REAL(dp)                                           :: x(3, nat)
      INTEGER                                            :: ihg, ihc
      REAL(dp)                                           :: rx(3, nat)
      INTEGER                                            :: nc, indpg, ntvec
      REAL(dp)                                           :: a(3, 3), ai(3, 3)
      INTEGER                                            :: li, isy, isc(nat)
      REAL(dp)                                           :: delta

      CHARACTER(len=10), DIMENSION(48), PARAMETER :: rname_cubic = [' 1        ', ' 2[ 10 0] ', &
         ' 2[ 01 0] ', ' 2[ 00 1] ', ' 3[-1-1-1]', ' 3[ 11-1] ', ' 3[-11 1] ', ' 3[ 1-11] ', &
         ' 3[ 11 1] ', ' 3[-11-1] ', ' 3[-1-11] ', ' 3[ 1-1-1]', ' 2[-11 0] ', ' 4[ 00 1] ', &
         ' 4[ 00-1] ', ' 2[ 11 0] ', ' 2[ 0-11] ', ' 2[ 01 1] ', ' 4[ 10 0] ', ' 4[-10 0] ', &
         ' 2[-10 1] ', ' 4[ 0-10] ', ' 2[ 10 1] ', ' 4[ 01 0] ', '-1        ', '-2[ 10 0] ', &
         '-2[ 01 0] ', '-2[ 00 1] ', '-3[-1-1-1]', '-3[ 11-1] ', '-3[-11 1] ', '-3[ 1-11] ', &
         '-3[ 11 1] ', '-3[-11-1] ', '-3[-1-11] ', '-3[ 1-1-1]', '-2[-11 0] ', '-4[ 00 1] ', &
         '-4[ 00-1] ', '-2[ 11 0] ', '-2[ 0-11] ', '-2[ 01 1] ', '-4[ 10 0] ', '-4[-10 0] ', &
         '-2[-10 1] ', '-4[ 0-10] ', '-2[ 10 1] ', '-4[ 01 0] ']
      CHARACTER(len=11), DIMENSION(24), PARAMETER :: rname_hexai = [' 1         ', ' 6[ 00  1] ', &
         ' 3[ 00  1] ', ' 2[ 00  1] ', ' 3[ 00 -1] ', ' 6[ 00 -1] ', ' 2[ 01  0] ', ' 2[-11  0] ', &
         ' 2[ 10  0] ', ' 2[ 21  0] ', ' 2[ 11  0] ', ' 2[ 12  0] ', '-1         ', '-6[ 00  1] ', &
         '-3[ 00  1] ', '-2[ 00  1] ', '-3[ 00 -1] ', '-6[ 00 -1] ', '-2[ 01  0] ', '-2[-11  0] ', &
         '-2[ 10  0] ', '-2[ 21  0] ', '-2[ 11  0] ', '-2[ 12  0] ']
      CHARACTER(len=12), DIMENSION(7), PARAMETER :: icst = ['TRICLINIC   ', 'MONOCLINIC  ', &
         'ORTHORHOMBIC', 'TETRAGONAL  ', 'CUBIC       ', 'TRIGONAL    ', 'HEXAGONAL   ']
      CHARACTER(len=3), DIMENSION(32), PARAMETER :: pgrd = ['c1 ', 'ci ', 'c2 ', 'c1h', 'c2h', &
         'c3 ', 'c3i', 'd3 ', 'c3v', 'd3 ', 'c4 ', 's4 ', 'c4h', 'd4 ', 'c4v', 'd2d', 'd4h', 'c6 ',&
         'c3h', 'c6h', 'd6 ', 'c6v', 'd3h', 'd6h', 'd2 ', 'c2v', 'd2h', 't  ', 'th ', 'o  ', 'td ',&
         'oh ']
      CHARACTER(len=5), DIMENSION(32), PARAMETER :: pgrp = ['    1', '  <1>', '    2', '    m', &
         '  2/m', '    3', '  <3>', '   32', '   3m', ' <3>m', '    4', '  <4>', '  4/m', '  422', &
         '  4mm', '<4>2m', '4/mmm', '    6', '  <6>', '  6/m', '  622', '  6mm', '<6>m2', '6/mmm', &
         '  222', '  mm2', '  mmm', '   23', '   m3', '  432', '<4>3m', '  m3m']

      INTEGER                                            :: i, iis(48), il, info, j, k, k2, l, n, &
                                                            nca, ni
      LOGICAL                                            :: nodupli, oksym
      REAL(dp)                                           :: vc(3, 48), vr(3), vs, xb(3)

      nodupli = ntvec == 1
      nca = 0
      DO n = 1, 48
         iis(n) = 0
      END DO
      ! Calculate translational vector for each operation
      ! and atom transformation table.
      DO n = 1, nc
         l = ib(n)
         iis(l) = 1
         DO k = 1, nat
            DO i = 1, 3
               rx(i, k) = r(i, 1, l)*x(1, k) + r(i, 2, l)*x(2, k) + r(i, 3, l)*x(3, k)
            END DO
         END DO
         DO k = 1, 3
            vr(k) = 0._dp
         END DO
         ! First we determine for VR=(/0,0,0/)
         ! IMPORTANT IF NOT UNIQUE ATOMS FOR DETERMINATION OF SYMMORPHIC
         CALL checkrlv3(n, nat, ty, rx, x, vr, f0, ai, isc, nodupli, oksym, delta)
         IF (.NOT. oksym) THEN
            ! Now we try other possible VR
            ! F0(49,1:NAT) has only inequivalent atom indexes for translation
            DO k2 = 1, nat
               IF (f0(49, k2) < k2) CYCLE
               IF (ty(1) /= ty(k2)) CYCLE
               DO i = 1, 3
                  xb(i) = rx(i, 1) - x(i, k2)
               END DO
               ! A translation vector VR is defined.
               CALL rlv3(ai, xb, vr, il, delta)
               ! ==----------------------------------------------------------==
               ! == SUBROUTINE RLV3 REMOVES A DIRECT LATTICE VECTOR FROM XB  ==
               ! == LEAVING THE REMAINDER IN VR. IF A NONZERO LATTICE        ==
               ! == VECTOR WAS REMOVED, IL IS MADE NONZERO.                  ==
               ! == VR STANDS FOR V-REFERENCE.                               ==
               ! == VR IS NOT GIVEN IN CARTESIAN COORDINATES BUT             ==
               ! == IN THE SYSTEM A1,A2,A3.     K.K., 23.10.1979             ==
               ! ==----------------------------------------------------------==
               CALL checkrlv3(n, nat, ty, rx, x, vr, f0, ai, isc, nodupli, oksym, delta)
               IF (oksym) EXIT
            END DO
            IF (.NOT. oksym) THEN
               iis(l) = 0
               CYCLE
            END IF
         END IF
         nca = nca + 1
         DO i = 1, 3
            v(i, nca) = vr(i)
         END DO
         ! ==------------------------------------------------------------==
         ! == V(I,N) IS THE I-TH COMPONENT OF THE FRACTIONAL             ==
         ! == TRANSLATION ASSOCIATED WITH THE ROTATION N.                ==
         ! == ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS, THEY ARE     ==
         ! == GIVEN IN THE SYSTEM A1,A2,A3.                              ==
         ! == K.K., 23.10. 1979                                          ==
         ! ==------------------------------------------------------------==
      END DO
      ! Remove unused operations
      i = 0
      ni = 13
      IF (ihg < 6) ni = 25
      li = 0
      DO n = 1, nc
         l = ib(n)
         IF (iis(l) == 0) CYCLE
         i = i + 1
         ib(i) = ib(n)
         IF (ib(i) == ni) li = i
         DO k = 1, nat
            f0(i, k) = f0(n, k)
         END DO
      END DO
      ! ==--------------------------------------------------------------==
      nc = i
      vs = 0._dp
      DO n = 1, nc
         vs = vs + ABS(v(1, n)) + ABS(v(2, n)) + ABS(v(3, n))
      END DO
      ! THE ORIGINAL VALUE DELTA=0.0001 WAS MODIFIED
      ! BY K.K. , SEPTEMBER 1979 TO 0.0005
      ! AND RETURNED TO 0.0001 BY RJN OCT 1987
      IF (vs > delta) THEN
         isy = 0
      ELSE
         isy = 1
      END IF
      ! ==--------------------------------------------------------------==
      ! Determination of the point group
      ! (Thierry Deutsch - 1998 [Maybe not complete!!])
      IF (ihg < 6) THEN
         IF (nc == 0) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(" ATFTM1! IHG=",A," NC=",I2)') icst(ihg), nC
            END IF
            CPABORT('ATFTM1: NUMBER OF ROTATION NULL')
            ! Triclinic system
         ELSE IF (nc == 1) THEN
            ! IB=1
            indpg = 1           ! 1 (c1)
         ELSE IF (nc == 2 .AND. ib(2) == 25) THEN
            ! IB=125
            indpg = 2           ! <1>(ci)
         ELSE IF (nc == 2 .AND. ( &
                  ib(2) == 4 .OR. & ! 2[001]
                  ib(2) == 2 .OR. & ! 2[100]
                  ib(2) == 3)) THEN ! 2[010]
            ! Monoclinic system
            ! IB=14 (z-axis) OR
            ! IB=12 (x-axis) OR
            ! IB=13 (y-axis)
            indpg = 3           ! 2 (c2)
         ELSE IF (nc == 2 .AND. ( &
                  ib(2) == 28 .OR. &
                  ib(2) == 26 .OR. &
                  ib(2) == 27)) THEN
            ! IB=128 (z-axis) OR
            ! IB=126 (x-axis) OR
            ! IB=127 (y-axis)
            indpg = 4           ! m (c1h)
         ELSE IF (nc == 4 .AND. ( &
                  ib(4) == 28 .OR. & ! 2[001]
                  ib(4) == 27 .OR. & ! 2[010]
                  ib(4) == 26 .OR. & ! 2[100]
                  ib(4) == 37 .OR. & ! -2[-110]
                  ib(4) == 40)) THEN ! 2[110]
            ! IB=1  425 28 (z-axis)  OR
            ! IB=1  225 26 (x-axis)  OR
            ! IB=1  325 27 (y-axis)  OR
            ! IB=113 2537 (-xy-axis)OR
            ! IB=116 2540 (xy-axis)
            indpg = 5           ! 2/m(c2h)
         ELSE IF (nc == 4 .AND. ( &
                  ib(4) == 15 .OR. &
                  ib(4) == 20 .OR. &
                  ib(4) == 24)) THEN
            ! Tetragonal system
            ! IB=14 1415 (z-axis) OR
            ! IB=12 1920 (x-axis) OR
            ! IB=13 2224 (y-axis)
            indpg = 11          ! 4 (c4)
         ELSE IF (nc == 4 .AND. ( &
                  ib(4) == 39 .OR. &
                  ib(4) == 44 .OR. &
                  ib(4) == 48)) THEN
            ! IB=14 3839 (z-axis) OR
            ! IB=12 4344 (x-axis) OR
            ! IB=13 4648 (y-axis)
            indpg = 12          ! <4>(s4)
         ELSE IF (nc == 8 .AND. ( &
                  (ib(3) == 14 .AND. ib(8) == 39) .OR. &
                  (ib(3) == 19 .AND. ib(8) == 44) .OR. &
                  (ib(3) == 22 .AND. ib(8) == 48))) THEN
            ! IB=14 1415 2825 3839 (z-axis) OR
            ! IB=12 1920 2625 4344 (x-axis) OR
            ! IB=13 2224 2725 4648 (y-axis)
            indpg = 13          ! 422(d4)
         ELSE IF (nc == 8 .AND. ib(4) == 4 .AND. ( &
                  ib(8) == 16 .OR. &
                  ib(8) == 20 .OR. &
                  ib(8) == 24)) THEN
            ! IB=12 3 413 1415 16 (z-axis) OR
            ! IB=12 3 417 1920 18 (x-axis) OR
            ! IB=12 3 421 2224 23 (y-axis)
            indpg = 14          ! 4/m(c4h)
         ELSE IF (nc == 8 .AND. ( &
                  ib(8) == 40 .OR. &
                  ib(8) == 42 .OR. &
                  ib(8) == 47)) THEN
            ! IB=14 1415 2627 3740 (z-axis) OR
            ! IB=12 1920 2827 4142 (x-axis) OR
            ! IB=13 2224 2628 4547 (y-axis)
            indpg = 15          ! 4mm(c4v)
         ELSE IF (nc == 8 .AND. ( &
                  (ib(3) == 13 .AND. ib(8) == 39) .OR. &
                  (ib(3) == 17 .AND. ib(8) == 44) .OR. &
                  (ib(3) == 21 .AND. ib(8) == 48))) THEN
            ! IB=14 1316 2627 3839 (z-axis) OR
            ! IB=12 1718 2827 4344 (x-axis) OR
            ! IB=13 2123 2628 4648 (y-axis)
            indpg = 16          ! <4>2m(d2d)
         ELSE IF (nc == 16 .AND. ( &
                  ib(16) == 40 .OR. &
                  ib(16) == 44 .OR. &
                  ib(16) == 48)) THEN
            ! IB=12 3 413 1415 1625 2627 2837 3839 40 (z-axis) OR
            ! IB=12 3 417 1920 1825 2627 2841 4344 42 (x-axis) OR
            ! IB=12 3 421 2224 2325 2627 2845 4648 47 (y-axis)
            indpg = 17          ! 4/mmm(d4h)
         ELSE IF (nc == 4 .AND. (ib(4) == 4)) THEN
            ! Orthorhombic system
            ! IB=12  3  4
            indpg = 25          ! 222(d2)
         ELSE IF (nc == 4 .AND. ( &
                  ib(4) == 27 .OR. &
                  ib(4) == 28)) THEN
            ! IB=13 2627 (z-axis) OR
            ! IB=12 2728 (x-axis) OR
            ! IB=14 2628 (y-axis) OR
            indpg = 26          ! mm2(c2v)
         ELSE IF (nc == 8) THEN
            ! IB=12 3 425 2627 28
            indpg = 27          ! mmm(d2h)
         ELSE IF (nc == 12 .AND. ( &
                  ib(12) == 12 .OR. &
                  ib(12) == 47 .OR. &
                  ib(12) == 45)) THEN
            ! Cubic system
            ! IB=12  3  4  5  6  7  8  910 1112 OR
            ! IB=15 1113 1823 2530 3537 4247 OR
            ! IB=18 1016 1821 2532 3440 4245
            indpg = 28          ! 23 (t)
         ELSE IF (nc == 24 .AND. ib(24) == 36) THEN
            ! IB= 1  2  3  4  5  6  7  8  910 1112
            ! 2526 2728 2930 3132 3334 3536
            indpg = 29          ! m3 (th)
         ELSE IF (nc == 24 .AND. ib(24) == 24) THEN
            ! IB=12 3 45 6 78 9 1011 12
            ! 1314 1516 1718 1920 2122 2324
            indpg = 30          ! 432 (o)
         ELSE IF (nc == 24 .AND. ib(24) == 48) THEN
            ! IB=12 3 45 6 78 9 1011 12
            ! 3738 3940 4142 4345 4647 48
            indpg = 31          ! <4>3m(td)
         ELSE IF (nc == 48) THEN
            ! IB=1..48
            indpg = 32          ! m3m(oh)
         ELSE
            ! WRITE(6,'(" ATFTM1! IHG=",A," NC=",I2)') ICST(IHG),NC
            ! WRITE(6,'(" ATFTM1!",19I3)') (IB(I),I=1,NC)
            ! WRITE(6,'(" ATFTM1! THIS CASE IS UNKNOWN IN THE DATABASE")')
            ! Probably a sub-group of 32
            indpg = -32
         END IF
      ELSE IF (ihg >= 6) THEN
         IF (nc == 0) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(" ATFTM1! IHG=",A," NC=",I2)') icst(ihg), nC
            END IF
            CPABORT('ATFTM1: NUMBER OF ROTATION NULL')
            ! Triclinic system
         ELSE IF (nc == 1) THEN
            ! IB=1
            indpg = 1           ! 1 (c1)
         ELSE IF (nc == 2 .AND. ib(2) == 13) THEN
            ! IB=113
            indpg = 2           ! <1>(ci)
         ELSE IF (nc == 2 .AND. ( &
                  ib(2) == 4)) THEN ! 2[001]
            ! Monoclinic system
            ! IB=1 4
            indpg = 3           ! 2 (c2)
         ELSE IF (nc == 2 .AND. ( &
                  ib(2) == 16)) THEN
            ! IB=116
            indpg = 4           ! m (c1h)
         ELSE IF (nc == 4 .AND. ( &
                  ib(4) == 24 .OR. &
                  ib(4) == 20)) THEN
            ! IB=112 1324 OR
            ! IB=1  813 20
            indpg = 5           ! 2/m(c2h)
         ELSE IF (nc == 3 .AND. ib(3) == 5) THEN
            ! Trigonal system
            ! IB=13 5
            indpg = 6           ! 3 (c3)
         ELSE IF (nc == 6 .AND. ib(6) == 17) THEN
            ! IB=113 1517 35
            indpg = 7           ! <3>(c3i)
         ELSE IF (nc == 6 .AND. ib(6) == 11) THEN
            ! IB=17 9 1135
            indpg = 8           ! 32 (d3)
         ELSE IF (nc == 6 .AND. ib(6) == 23) THEN
            ! IB=13 5 1921 23
            indpg = 9           ! 3m (c3v)
         ELSE IF (nc == 12 .AND. ib(12) == 23) THEN
            ! IB=13 5 79 1113 1517 1921 23
            indpg = 10          ! <3>m(d3d)
         ELSE IF (nc == 6 .AND. ib(6) == 6) THEN
            ! Hexagonal system
            ! IB=12 3 45 6
            indpg = 18          ! 6 (c6)
         ELSE IF (nc == 6 .AND. ib(6) == 18) THEN
            ! IB=13 5 1416 18
            indpg = 19          ! <6>(c3h)
         ELSE IF (nc == 12 .AND. ib(12) == 18) THEN
            ! IB=12 3 45 6 1314 1516 1718
            indpg = 20          ! 6/m(c6h)
         ELSE IF (nc == 12 .AND. ib(12) == 12) THEN
            ! IB=12 3 45 6 78 9 1011 12
            indpg = 21          ! 622(d6)
         ELSE IF (nc == 12 .AND. ib(2) == 2 .AND. ib(12) == 24) THEN
            ! IB=12 3 45 6 1920 2122 2324
            indpg = 22          ! 6mm(c6v)
         ELSE IF (nc == 12 .AND. ib(2) == 3 .AND. ib(12) == 24) THEN
            ! IB=13 5 79 1114 1618 2022 24
            indpg = 23          ! <6>m2(d3h)
         ELSE IF (nc == 24) THEN
            ! IB=1..24
            indpg = 24          ! 6/mmm(d6h)
         ELSE
            ! Probably a sub-group of 24
            ! WRITE(6,'(" ATFTM1! IHG=",A," NC=",I2)') ICST(IHG),NC
            ! WRITE(6,'(" ATFTM1!",48I3)') (IB(I),I=1,NC)
            ! WRITE(6,'(" ATFTM1! THIS CASE IS UNKNOWN IN THE DATABASE")')
            indpg = -24
         END IF
      END IF
      ! ==--------------------------------------------------------------==
      ! == Determination if the space group is symmorphic or not        ==
      ! ==--------------------------------------------------------------==
      IF (isy /= 1) THEN
         ! Transform V in cartesian coordinates
         DO n = 1, nc
            vc(1, n) = a(1, 1)*v(1, n) + a(1, 2)*v(2, n) + a(1, 3)*v(3, n)
            vc(2, n) = a(2, 1)*v(1, n) + a(2, 2)*v(2, n) + a(2, 3)*v(3, n)
            vc(3, n) = a(3, 1)*v(1, n) + a(3, 2)*v(2, n) + a(3, 3)*v(3, n)
         END DO
         CALL symmorphic(nc, ib, r, vc, ai, info, origin, delta)
         IF (info == 1) THEN
            CALL rlv3(ai, origin, xb, il, delta)
            ! !!!RLV3 determines -XB in crystal coordinates
            ! !!We want between 0.0 and 1.0
            DO i = 1, 3
               IF (-xb(i) >= 0._dp) THEN
                  origin(i) = -xb(i)
               ELSE
                  origin(i) = 1._dp - xb(i)
               END IF
            END DO
            DO i = 1, 3
               xb(i) = a(i, 1)*origin(1) + a(i, 2)*origin(2) + a(i, 3)*origin(3)
            END DO
            isy = -1
         ELSE IF (info == 0) THEN
            isy = 0
         ELSE
            isy = -2
         END IF
      ELSE
         DO i = 1, 3
            origin(i) = 0._dp
         END DO
      END IF
      ! ==--------------------------------------------------------------==
      ! == Output                                                       ==
      ! ==--------------------------------------------------------------==
      IF (iout > 0) THEN
         IF (iout > 0) THEN
            WRITE (iout, *)
         END IF
         CALL xstring(icst(ihg), i, j)
         IF ((ihg == 7 .AND. nc == 24) .OR. &
             (ihg == 5 .AND. nc == 48)) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(A,A,A)') &
                  ' KPSYM| THE POINT GROUP OF THE CRYSTAL IS THE FULL ', &
                  icst(ihg) (i:j), &
                  ' GROUP'
            END IF
         ELSE
            IF (iout > 0) THEN
               WRITE (iout, '(A,A,A,I2,A)') &
                  ' KPSYM| THE CRYSTAL SYSTEM IS ', &
                  icst(ihg) (i:j), &
                  ' WITH ', nc, ' OPERATIONS:'
            END IF
            IF (ihc == 0) THEN
               IF (iout > 0) THEN
                  WRITE (iout, '( 5(5(A13),/))') (rname_hexai(ib(i)), i=1, nc)
               END IF
            ELSE
               IF (iout > 0) THEN
                  WRITE (iout, '(10(5(A13),/))') (rname_cubic(ib(i)), i=1, nc)
               END IF
            END IF
         END IF
         ! ==------------------------------------------------------------==
         IF (isy == 1) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(A)') &
                  ' KPSYM| THE SPACE GROUP OF THE CRYSTAL IS SYMMORPHIC'
            END IF
         ELSE IF (isy == -1) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(A)') &
                  ' KPSYM| THE SPACE GROUP OF THE CRYSTAL IS SYMMORPHIC'
            END IF
            IF (iout > 0) THEN
               WRITE (iout, '(A,A,/,T3,3F10.6,3X,3F10.6)') &
                  ' KPSYM| THE STANDARD ORIGIN OF COORDINATES IS:   ', &
                  '[CARTESIAN]   [CRYSTAL]', xb, origin
            END IF
         ELSE IF (isy == 0) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(A,/,3X,A,F15.6,A)') &
                  ' KPSYM| THE SPACE GROUP IS NON-SYMMORPHIC,', &
                  ' (SUM OF TRANSLATION VECTORS=', vs, ')'
            END IF
         ELSE IF (isy == -2) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(A,A)') &
                  ' KPSYM| CANNOT DETERMINE IF THE SPACE GROUP IS', &
                  ' SYMMORPHIC OR NOT'
            END IF
            IF (iout > 0) THEN
               WRITE (iout, '(A,/,A,/,3X,A,F15.6,A)') &
                  ' KPSYM| THE SPACE GROUP IS NON-SYMMORPHIC,', &
                  ' KPSYM| OR ELSE A NON STANDARD ORIGIN OF COORDINATES WAS USED.', &
                  ' KPSYM| (SUM OF TRANSLATION VECTORS=', vs, ')'
            END IF
         END IF
         IF (indpg > 0) THEN
            CALL xstring(pgrp(indpg), i, j)
            CALL xstring(pgrd(indpg), k, l)
            IF (iout > 0) THEN
               WRITE (iout, '(A,A,"(",A,")",T56,"[INDEX=",I2,"]")') &
                  ' KPSYM| THE POINT GROUP OF THE CRYSTAL IS  ', pgrp(indpg) (i:j), &
                  pgrd(indpg) (k:l), indpg
            END IF
         ELSE
            CALL xstring(pgrp(-indpg), i, j)
            CALL xstring(pgrd(-indpg), k, l)
            IF (iout > 0) THEN
               WRITE (iout, '(A,I2,A,A,"(",A,")",T56,"[INDEX=",I2,"]")') &
                  ' KPSYM| POINT GROUP: GROUP ORDER=', nc, &
                  '    SUBGROUP OF ', pgrp(-indpg) (i:j), &
                  pgrd(-indpg) (k:l), -indpg
            END IF
         END IF
         IF (ntvec == 1) THEN
            IF (iout > 0) THEN
               WRITE (iout, '(A,T60,I6)') &
                  ' KPSYM| NUMBER OF PRIMITIVE CELL:', ntvec
            END IF
         ELSE
            IF (iout > 0) THEN
               WRITE (iout, '(A,T60,I6)') &
                  ' KPSYM| NUMBER OF PRIMITIVE CELLS:', ntvec
            END IF
         END IF
      END IF

   END SUBROUTINE atftm1

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param nat ...
!> \param ty ...
!> \param rx ...
!> \param x ...
!> \param vr ...
!> \param f0 ...
!> \param ai ...
!> \param isc ...
!> \param nodupli ...
!> \param oksym ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE checkrlv3(n, nat, ty, rx, x, vr, f0, ai, isc, &
                        nodupli, oksym, delta)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN IN MAY 14TH, 1998 (T.D.)                             ==
      ! == CHECK IF RX+VR GIVES THE SAME LATTICE AS X                   ==
      ! == BUILD THE ATOM TRANSFORMATION TABLE                          ==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==   N   ROTATION NUMBER (INDEX USED IN F0 BETWEEN 1 AND 48)    ==
      ! ==   NAT           NUMBER OF ATOMS                              ==
      ! ==   TY(1:NAT)     TYPE OF ATOMS                                ==
      ! ==   RX(1:3,1:NAT) ATOMIC COORDINATES FROM Nth ROTATION (CART.) ==
      ! ==   X(1:3,1:NAT)  ATOMIC COORDINATES (CARTESIAN)               ==
      ! ==   VR(1:3)       TRANSLATION VECTOR (CRYSTAL COOR.)           ==
      ! ==   AI(1:3,1:3)   LATTICE RECIPROCAL VECTORS                   ==
      ! ==   NODUPLI       .TRUE., THE CELL IS A PRIMITIVE ONE          ==
      ! ==                 WE CAN SPEED UP                              ==
      ! == DELTA           REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)    ==
      ! == OUTPUT:                                                      ==
      ! ==   F0(1:49,1:NAT) ATOM TRANSFORMATION TABLE                   ==
      ! ==      F0 IS THE FUNCTION DEFINED IN MARADUDIN AND VOSK0       ==
      ! ==      BY EQ.(2.35).                                           ==
      ! ==      IT DEFINES THE ATOM TRANSFORMATION TABLE                ==
      ! ==   OKSYM          TRUE IF RX+VR = X                           ==
      ! ==   ISC(1:NAT)     SCRATCH ARRAY                               ==
      ! ==                  USED TO SPEED UP THE ROUTINE                ==
      ! ==                  EACH ATOM IS ONLY ONCE AN IMAGE             ==
      ! ==                  IF NO DUPLICATION OF THE CELL               ==
      ! ==--------------------------------------------------------------==
      INTEGER                                            :: n, nat, ty(nat)
      REAL(dp)                                           :: rx(3, nat), x(3, nat), vr(3)
      INTEGER                                            :: f0(49, nat)
      REAL(dp)                                           :: ai(3, 3)
      INTEGER                                            :: isc(nat)
      LOGICAL                                            :: nodupli, oksym
      REAL(dp)                                           :: delta

      INTEGER                                            :: ia, ib, il
      REAL(dp)                                           :: tol, vt(3), xb(3)

      ! Fractional residuals at the tolerance boundary accumulate Cartesian
      ! rotation and lattice-conversion roundoff. Account for that roundoff so
      ! equivalent operations are not accepted or rejected by a few ulps.
      tol = delta + 32.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, &
                                                MAXVAL(SUM(ABS(ai), DIM=2))*MAX(MAXVAL(ABS(rx)), MAXVAL(ABS(x))))

      DO ia = 1, nat
         isc(ia) = 0
      END DO
      ! Now we check if ROT(N)+VR gives a correct symmetry.
      atom: DO ia = 1, nat
         DO ib = 1, nat
            IF (ty(ia) == ty(ib) .AND. isc(ib) == 0) THEN
               xb(1) = rx(1, ia) - x(1, ib)
               xb(2) = rx(2, ia) - x(2, ib)
               xb(3) = rx(3, ia) - x(3, ib)
               CALL rlv3(ai, xb, vt, il, delta)
               ! VT STANDS FOR V-TEST
               oksym = ALL(ABS((vr - vt) - ANINT(vr - vt)) <= tol)
               IF (oksym) THEN
                  IF (nodupli) isc(ib) = 1
                  f0(n, ia) = ib
                  ! IR+VR is the good one: another symmetry operation
                  ! Next atom
                  CYCLE atom
               END IF
            END IF
         END DO
         ! VR is not the correct translation vector
         RETURN
      END DO atom
   END SUBROUTINE checkrlv3
   ! ==================================================================
! **************************************************************************************************
!> \brief ...
!> \param nc ...
!> \param ib ...
!> \param r ...
!> \param v ...
!> \param ai ...
!> \param info ...
!> \param origin ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE symmorphic(nc, ib, r, v, ai, info, origin, delta)
      ! ==--------------------------------------------------------------==
      ! == Check if the group is symmorphic with a non-standard origin  ==
      ! == WARNING: If there are equivalent atoms, this routine could   ==
      ! == not determine if the space group is symmorphic               ==
      ! == So you have to check if the solution V=0 works (see ATFTM1)  ==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==   NC Number of operations                                    ==
      ! ==   IB(NC) Index of operation in R                             ==
      ! ==   R(3,3,48) Rotations                                        ==
      ! ==   V(3,NC) Fractional translations related to R(3,3,IB(NC))   ==
      ! ==           R AND V ARE IN CARTESIAN COORDINATES               ==
      ! ==   AI(I,J) ARE THE RECIPROCAL LATTICE VECTORS,                ==
      ! ==           B(I) = AI(I,J),J=1,2,3                             ==
      ! == DELTA     REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)          ==
      ! ==                                                              ==
      ! == OUTPUT:                                                      ==
      ! ==   ORIGIN(1:3) Give standard origin (cartesian coordinates)   ==
      ! ==            Give the standard origin with smallest coordinates==
      ! ==            if NTVEC /= 1                                     ==
      ! ==   INFO = 1 The group is symmorphic                           ==
      ! ==   INFO = 0 The group is not symmorphic                       ==
      ! ==   INFO =-1 The routine cannot determine                      ==
      ! ==--------------------------------------------------------------==
      INTEGER                                            :: nc, ib(nc)
      REAL(dp)                                           :: r(3, 3, 48), v(3, nc), ai(3, 3)
      INTEGER                                            :: info
      REAL(dp)                                           :: origin(3), delta

      INTEGER                                            :: i, i1, ierror, igood(3), il, imissing2, &
                                                            imissing3, iok(3), ionly, ir, j, j1
      REAL(dp)                                           :: diag, dif, r2(2, 2), r3(3, 3), vr(3), &
                                                            xb(3)

! Variables
! ==--------------------------------------------------------------==
! Find a point A / V_R = (1-R).OA

      DO i = 1, 3
         iok(i) = 0
      END DO
      DO i = 1, 3
         origin(i) = 0._dp
      END DO
      DO ir = 1, nc
         dif = v(1, ir)*v(1, ir) + v(2, ir)*v(2, ir) + v(3, ir)*v(3, ir)
         IF (dif > delta*delta) THEN
            DO i = 1, 3
               igood(i) = 1
            END DO
            ! V is non-zero. Construct matrix 1-R
            DO i = 1, 3
               DO j = 1, 3
                  r3(i, j) = -r(i, j, ib(ir))
               END DO
               r3(i, i) = 1 + r3(i, i)
            END DO
            CALL invmat(r3, ierror)
            IF (ierror == 0) THEN
               ! The matrix 3x3 has an inverse.
               DO i = 1, 3
                  vr(i) = r3(i, 1)*v(1, ir) &
                          + r3(i, 2)*v(2, ir) &
                          + r3(i, 3)*v(3, ir)
               END DO
            ELSE
               ! IERROR gives the column which causes some trouble
               ! Construct matrix 1-R with 2x2
               igood(ierror) = 0
               imissing3 = ierror
               i1 = 0
               DO i = 1, 3
                  IF (i /= ierror) THEN
                     i1 = i1 + 1
                     j1 = 0
                     DO j = 1, 3
                        IF (j /= ierror) THEN
                           j1 = j1 + 1
                           r2(i1, j1) = -r(i, j, ib(ir))
                        END IF
                     END DO
                     r2(i1, i1) = 1 + r2(i1, i1)
                  END IF
               END DO
               CALL invmat(r2, ierror)
               IF (ierror == 0) THEN
                  ! The matrix 2X2 has an inverse.
                  ! Solve Vxy = (1-R).OAxy + OAz R3z (z is IMISSING3)
                  i1 = 0
                  DO i = 1, 3
                     IF (igood(i) == 1) THEN
                        i1 = i1 + 1
                        vr(i) = 0._dp
                        j1 = 0
                        DO j = 1, 3
                           IF (igood(j) == 1) THEN
                              j1 = j1 + 1
                              vr(i) = vr(i) + r2(i1, j1)*(v(j, ir) + &
                                                          origin(imissing3)*r(j, imissing3, ib(ir)))
                           END IF
                        END DO
                     ELSE
                        vr(i) = origin(i)
                     END IF
                  END DO
               ELSE
                  ! Construct matrix 1-R with 1x1
                  i1 = 0
                  DO i = 1, 3
                     IF (i /= imissing3) THEN
                        i1 = i1 + 1
                        IF (i1 == ierror) THEN
                           igood(i) = 0
                           imissing2 = i
                        ELSE
                           ionly = i
                        END IF
                     END IF
                  END DO
                  diag = (1 - r(ionly, ionly, ib(ir)))
                  IF (ABS(diag) > delta) THEN
                     vr(ionly) = 1._dp/diag*(v(ionly, ir) + &
                                             origin(imissing3)*r(ionly, imissing3, ib(ir)) + &
                                             origin(imissing2)*r(ionly, imissing2, ib(ir)))
                  ELSE
                     vr(ionly) = origin(ionly)
                     igood(ionly) = 0
                  END IF
                  vr(imissing3) = origin(imissing3)
                  vr(imissing2) = origin(imissing2)
               END IF
            END IF
            ! ==----------------------------------------------------------==
            ! Compare VR with ORIGIN
            dif = 0._dp
            ! If NTVEC /=1 there are NTVEC possible standard origins
            DO i = 1, 3
               IF (iok(i) == 1) THEN
                  dif = dif + ABS(origin(i) - vr(i))
               END IF
            END DO
            IF (dif > delta) THEN
               ! Non-symmorphic
               info = 0
               RETURN
            ELSE
               DO i = 1, 3
                  IF (iok(i) /= 1 .AND. igood(i) == 1) THEN
                     iok(i) = 1
                     origin(i) = vr(i)
                  END IF
               END DO
            END IF
         END IF
      END DO
      ! ==--------------------------------------------------------------==
      IF (iok(1) == 0 .AND. iok(2) == 0 .AND. iok(3) == 0) THEN
         ! Cannot not determine
         info = -1
         RETURN
      END IF
      ! The group is symmorphic
      info = 1
      ! Check
      DO ir = 1, nc
         DO i = 1, 3
            vr(i) = r(i, 1, ib(ir))*origin(1) &
                    + r(i, 2, ib(ir))*origin(2) &
                    + r(i, 3, ib(ir))*origin(3)
            vr(i) = (origin(i) - vr(i)) - v(i, ir)
         END DO
         CALL rlv3(ai, vr, xb, il, delta)
         dif = ABS(xb(1)) + ABS(xb(2)) + ABS(xb(3))
         IF (dif > delta) THEN
            ! Non-symmorphic
            info = 0
            RETURN
         END IF
      END DO
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE symmorphic
   ! ==================================================================
! **************************************************************************************************
!> \brief ...
!> \param ihc ...
!> \param r ...
! **************************************************************************************************
   SUBROUTINE rot1(ihc, r)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN ON FEBRUARY 17TH, 1976                               ==
      ! == GENERATION OF THE X,Y,Z-TRANSFORMATION MATRICES 3X3          ==
      ! == FOR HEXAGONAL AND CUBIC GROUPS                               ==
      ! == SUBROUTINES NEEDED -- NONE                                   ==
      ! ==--------------------------------------------------------------==
      ! == THIS IS IDENTICAL WITH THE SUBROUTINE ROT OF WORLTON-WARREN  ==
      ! == (IN THE AC-COMPLEX), ONLY THE WAY OF TRANSFERRING THE DATA   ==
      ! == WAS CHANGED                                                  ==
      ! ==--------------------------------------------------------------==
      ! == INPUT DATA:                                                  ==
      ! == IHC SWITCH DETERMINING IF WE DESIRE                          ==
      ! ==     THE HEXAGONAL GROUP(IHC=0) OR THE CUBIC GROUP (IHC=1)    ==
      ! == OUTPUT DATA:                                                 ==
      ! == R...THE 3X3 MATRICES OF THE DESIRED COORDINATE REPRESENTATION==
      ! ==     THEIR NUMBERING CORRESPONDS TO THE SYMMETRY ELEMENTS AS  ==
      ! ==     LISTE IN WORLTON-WARREN                                  ==
      ! ==              (COMPUT. PHYS. COMM. 3(1972) 88--117)           ==
      ! == FOR IHC=0 THE FIRST 24 MATRICES OF THE ARRAY R REPRESENT     ==
      ! ==           THE FULL HEXAGONAL GROUP D(6H)                     ==
      ! == FOR IHC=1 THE FIRST 48 MATRICES OF THE ARRAY R REPRESENT     ==
      ! ==           THE FULL CUBIC GROUP O(H)                          ==
      ! ==--------------------------------------------------------------==
      INTEGER                                            :: ihc
      REAL(dp)                                           :: r(3, 3, 48)

      INTEGER                                            :: i, j, k, n, nv
      REAL(dp)                                           :: c, s

      DO j = 1, 3
         DO i = 1, 3
            DO n = 1, 48
               r(i, j, n) = 0._dp
            END DO
         END DO
      END DO
      IF (ihc == 0) THEN
         ! ==------------------------------------------------------------==
         ! DEFINE THE GENERATORS FOR THE ROTATION MATRICES--HEXAGONAL GROUP
         ! ==------------------------------------------------------------==
         c = 0.5_dp
         s = 0.5_dp*SQRT(3.0_dp)
         r(1, 1, 2) = c
         r(1, 2, 2) = -s
         r(2, 1, 2) = s
         r(2, 2, 2) = c
         r(1, 1, 7) = -c
         r(1, 2, 7) = -s
         r(2, 1, 7) = -s
         r(2, 2, 7) = c
         DO n = 1, 6
            r(3, 3, n) = 1._dp
            r(3, 3, n + 18) = 1._dp
            r(3, 3, n + 6) = -1._dp
            r(3, 3, n + 12) = -1._dp
         END DO
         ! ==------------------------------------------------------------==
         ! == GENERATE THE REST OF THE ROTATION MATRICES                 ==
         ! ==------------------------------------------------------------==
         DO i = 1, 2
            r(i, i, 1) = 1._dp
            DO j = 1, 2
               r(i, j, 6) = r(j, i, 2)
               DO k = 1, 2
                  r(i, j, 3) = r(i, j, 3) + r(i, k, 2)*r(k, j, 2)
                  r(i, j, 8) = r(i, j, 8) + r(i, k, 2)*r(k, j, 7)
                  r(i, j, 12) = r(i, j, 12) + r(i, k, 7)*r(k, j, 2)
               END DO
            END DO
         END DO
         DO i = 1, 2
            DO j = 1, 2
               r(i, j, 5) = r(j, i, 3)
               DO k = 1, 2
                  r(i, j, 4) = r(i, j, 4) + r(i, k, 2)*r(k, j, 3)
                  r(i, j, 9) = r(i, j, 9) + r(i, k, 2)*r(k, j, 8)
                  r(i, j, 10) = r(i, j, 10) + r(i, k, 12)*r(k, j, 3)
                  r(i, j, 11) = r(i, j, 11) + r(i, k, 12)*r(k, j, 2)
               END DO
            END DO
         END DO
         DO n = 1, 12
            nv = n + 12
            DO i = 1, 2
               DO j = 1, 2
                  r(i, j, nv) = -r(i, j, n)
               END DO
            END DO
         END DO
      ELSE
         ! ==------------------------------------------------------------==
         ! == DEFINE THE GENERATORS FOR THE ROTATION MATRICES-CUBIC GROUP==
         ! ==------------------------------------------------------------==
         r(1, 3, 9) = 1._dp
         r(2, 1, 9) = 1._dp
         r(3, 2, 9) = 1._dp
         r(1, 1, 19) = 1._dp
         r(2, 3, 19) = -1._dp
         r(3, 2, 19) = 1._dp
         DO i = 1, 3
            r(i, i, 1) = 1._dp
            DO j = 1, 3
               r(i, j, 20) = r(j, i, 19)
               r(i, j, 5) = r(j, i, 9)
               DO k = 1, 3
                  r(i, j, 2) = r(i, j, 2) + r(i, k, 19)*r(k, j, 19)
                  r(i, j, 16) = r(i, j, 16) + r(i, k, 9)*r(k, j, 19)
                  r(i, j, 23) = r(i, j, 23) + r(i, k, 19)*r(k, j, 9)
               END DO
            END DO
         END DO
         DO i = 1, 3
            DO j = 1, 3
               DO k = 1, 3
                  r(i, j, 6) = r(i, j, 6) + r(i, k, 2)*r(k, j, 5)
                  r(i, j, 7) = r(i, j, 7) + r(i, k, 16)*r(k, j, 23)
                  r(i, j, 8) = r(i, j, 8) + r(i, k, 5)*r(k, j, 2)
                  r(i, j, 10) = r(i, j, 10) + r(i, k, 2)*r(k, j, 9)
                  r(i, j, 11) = r(i, j, 11) + r(i, k, 9)*r(k, j, 2)
                  r(i, j, 12) = r(i, j, 12) + r(i, k, 23)*r(k, j, 16)
                  r(i, j, 14) = r(i, j, 14) + r(i, k, 16)*r(k, j, 2)
                  r(i, j, 15) = r(i, j, 15) + r(i, k, 2)*r(k, j, 16)
                  r(i, j, 22) = r(i, j, 22) + r(i, k, 23)*r(k, j, 2)
                  r(i, j, 24) = r(i, j, 24) + r(i, k, 2)*r(k, j, 23)
               END DO
            END DO
         END DO
         DO i = 1, 3
            DO j = 1, 3
               DO k = 1, 3
                  r(i, j, 3) = r(i, j, 3) + r(i, k, 5)*r(k, j, 12)
                  r(i, j, 4) = r(i, j, 4) + r(i, k, 5)*r(k, j, 10)
                  r(i, j, 13) = r(i, j, 13) + r(i, k, 23)*r(k, j, 11)
                  r(i, j, 17) = r(i, j, 17) + r(i, k, 16)*r(k, j, 12)
                  r(i, j, 18) = r(i, j, 18) + r(i, k, 16)*r(k, j, 10)
                  r(i, j, 21) = r(i, j, 21) + r(i, k, 12)*r(k, j, 15)
               END DO
            END DO
         END DO
         DO n = 1, 24
            nv = n + 24
            r(1, 1, nv) = -r(1, 1, n)
            r(1, 2, nv) = -r(1, 2, n)
            r(1, 3, nv) = -r(1, 3, n)
            r(2, 1, nv) = -r(2, 1, n)
            r(2, 2, nv) = -r(2, 2, n)
            r(2, 3, nv) = -r(2, 3, n)
            r(3, 1, nv) = -r(3, 1, n)
            r(3, 2, nv) = -r(3, 2, n)
            r(3, 3, nv) = -r(3, 3, n)
         END DO
      END IF
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE rot1
   ! ==================================================================
END MODULE kpsym
