!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief K-points and crystal symmetry routines based on
!     K290 code:
!     Written on September 12th, 1979.
!     IBM-retouched on October 27th, 1980.
!     Generation of special points modified on 26-May-82 by ohn.
!     Retouched on January 8th, 1997
!     Integration in CPMD-FEMD Program by Thierry Deutsch
!     ==--------------------------------------------------------------==
!     Playing with special points and creation of 'CRYSTALLOGRAPHIC'
!     File for band structure calculations.
!     Generation of special points in k-space for an arbitrary lattice,
!     Following the method Monkhorst,Pack, Phys. Rev. B13 (1976) 5188
!     Modified by Macdonald, Phys. Rev. B18 (1978) 5897
!     Modified also by Ole Holm Nielsen ("SYMMETRIZATION")
!     ==--------------------------------------------------------------==
!     (GROUP1, PGL1, ATFTM1, ROT1 FROM THE
!      "COMPUTER PHYSICS COMMUNICATIONS" PACKAGE "ACMI" - (1971,1974)
!      Worlton-Warren).
! **************************************************************************************************
MODULE kpsym

   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invmat
   USE string_utilities,                ONLY: xstring
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC  :: K290s, GROUP1s

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpsym'

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param iout ...
!> \param nat ...
!> \param nkpoint ...
!> \param nsp ...
!> \param iq1 ...
!> \param iq2 ...
!> \param iq3 ...
!> \param istriz ...
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \param alat ...
!> \param strain ...
!> \param xkapa ...
!> \param rx ...
!> \param tvec ...
!> \param ty ...
!> \param isc ...
!> \param f0 ...
!> \param ntvec ...
!> \param wvk0 ...
!> \param wvkl ...
!> \param lwght ...
!> \param lrot ...
!> \param nhash ...
!> \param includ ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE k290s(iout, nat, nkpoint, nsp, iq1, iq2, iq3, istriz, &
                    a1, a2, a3, alat, strain, xkapa, rx, tvec, &
                    ty, isc, f0, ntvec, wvk0, wvkl, lwght, lrot, &
                    nhash, includ, list, rlist, delta)
      ! ==================================================================
      ! WRITTEN ON SEPTEMBER 12TH, 1979.
      ! IBM-RETOUCHED ON OCTOBER 27TH, 1980.
      ! Tsukuba-retouched on March 19th, 2008.
      ! GENERATION OF SPECIAL POINTS MODIFIED ON 26-MAY-82 BY OHN.
      ! RETOUCHED ON JANUARY 8TH, 1997
      ! INTEGRATION IN CPMD-FEMD PROGRAM BY THIERRY DEUTSCH
      ! ==--------------------------------------------------------------==
      ! PLAYING WITH SPECIAL POINTS AND CREATION OF 'CRYSTALLOGRAPHIC'
      ! FILE FOR BAND STRUCTURE CALCULATIONS.
      ! GENERATION OF SPECIAL POINTS IN K-SPACE FOR AN ARBITRARY LATTICE,
      ! FOLLOWING THE METHOD MONKHORST,PACK, PHYS. REV. B13 (1976) 5188
      ! MODIFIED BY MACDONALD, PHYS. REV. B18 (1978) 5897
      ! MODIFIED ALSO BY OLE HOLM NIELSEN ("SYMMETRIZATION")
      ! ==--------------------------------------------------------------==
      ! TESTING THEIR EFFICIENCY AND PREPARATION OF THE
      ! "STRUCTURAL" FILE FOR RUNNING THE
      ! SELF-CONSISTENT BAND STRUCTURE PROGRAMS.
      ! IN THE CASES WHERE THE POINT GROUP OF THE CRYSTAL DOES NOT
      ! CONTAIN INVERSION, THE LATTER IS ARTIFICIALLY ADDED, IN ORDER
      ! TO MAKE USE OF THE HERMITICITY OF THE HAMILTONIAN
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==   IOUT    LOGIC FILE NUMBER                                  ==
      ! ==   NAT     NUMBER OF ATOMS                                    ==
      ! ==   NKPOINT MAXIMAL NUMBER OF K POINTS                         ==
      ! ==   NSP     NUMBER OF SPECIES                                  ==
      ! ==   IQ1,IQ2,IQ3 THE MONKHORST-PACK MESH PARAMETERS             ==
      ! ==   ISTRIZ  SWITCH FOR SYMMETRIZATION                          ==
      ! ==   A1(3),A2(3),A3(3) LATTICE VECTORS                          ==
      ! ==   ALAT    LATTICE CONSTANT                                   ==
      ! ==   STRAIN(3,3)  STRAIN APPLIED TO LATTICE IN ORDER            ==
      ! ==           TO HAVE K POINTS WITH SYMMETRY OF STRAINED LATTICE ==
      ! ==   XKAPA(3,NAT)   ATOMS COORDINATES                           ==
      ! ==   TY(NAT)      TYPES OF ATOMS                                ==
      ! ==   WVK0(3) SHIFT FOR K POINTS MESh (MACDONALD ARTICLE)        ==
      ! ==   NHASH   SIZE OF THE HASH TABLES (LIST)                     ==
      ! ==   DELTA   REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)       ==
      ! ==           K-VECTOR < DELTA IS CONSIDERED ZERO                ==
      ! == OUTPUT:                                                      ==
      ! ==   RX(3,NAT) SCRATCH ARRAY USED BY GROUP1 ROUTINE             ==
      ! ==   TVEC(1:3,1:NTVEC) TRANSLATION VECTORS (SEE NTVEC)          ==
      ! ==   ISC(NAT)  SCRATCH ARRAY USED BY GROUP1 ROUTINE             ==
      ! ==   F0(49,NAT) ATOM TRANSFORMATION TABLE                       ==
      ! ==              IF NTVEC/=1 THE 49TH GIVES INEQUIVALENT ATOMS   ==
      ! ==   NTVEC NUMBER OF TRANSLATION VECTORS (IF NOT PRIMITIVE CELL)==
      ! ==   WVKL(3,NKPOINT) SPECIAL KPOINTS GENERATED                  ==
      ! ==   LWGHT(NKPOINT)  WEIGHT FOR EACH K POINT                    ==
      ! ==   LROT(48,NKPOINT) SYMMETRY OPERATION FOR EACH K POINTS      ==
      ! ==   INCLUD(NKPOINT)  SCRATCH ARRAY USED BY SPPT2               ==
      ! ==   LIST(NKPOINT+NHASH) HASH TABLE USED BY SPPT2               ==
      ! ==   RLIST(3,NKPOINT) SCRATCH ARRAY USED BY SPPT2               ==
      ! ==--------------------------------------------------------------==
      ! SUBROUTINES NEEDED:
      ! SPPT2, GROUP1, PGL1, ATFTM1, ROT1, STRUCT,
      ! BZRDUC, INBZ, MESH, BZDEFI
      ! (GROUP1, PGL1, ATFTM1, ROT1 FROM THE
      ! "COMPUTER PHYSICS COMMUNICATIONS" PACKAGE "ACMI" - (1971,1974)
      ! WORLTON-WARREN).
      ! ==================================================================
      INTEGER                                            :: iout, nat, nkpoint, nsp, iq1, iq2, iq3, &
                                                            istriz
      REAL(KIND=dp)                                      :: a1(3), a2(3), a3(3), alat, strain(6), &
                                                            xkapa(3, nat), rx(3, nat), tvec(3, nat)
      INTEGER                                            :: ty(nat), isc(nat), f0(49, nat), ntvec
      REAL(KIND=dp)                                      :: wvk0(3), wvkl(3, nkpoint)
      INTEGER                                            :: lwght(nkpoint), lrot(48, nkpoint), &
                                                            nhash, includ(nkpoint), &
                                                            list(nkpoint + nhash)
      REAL(KIND=dp)                                      :: rlist(3, nkpoint), delta

      CHARACTER(len=10), DIMENSION(48), PARAMETER :: rname_cubic = (/' 1        ', ' 2[ 10 0] ', &
         ' 2[ 01 0] ', ' 2[ 00 1] ', ' 3[-1-1-1]', ' 3[ 11-1] ', ' 3[-11 1] ', ' 3[ 1-11] ', &
         ' 3[ 11 1] ', ' 3[-11-1] ', ' 3[-1-11] ', ' 3[ 1-1-1]', ' 2[-11 0] ', ' 4[ 00 1] ', &
         ' 4[ 00-1] ', ' 2[ 11 0] ', ' 2[ 0-11] ', ' 2[ 01 1] ', ' 4[ 10 0] ', ' 4[-10 0] ', &
         ' 2[-10 1] ', ' 4[ 0-10] ', ' 2[ 10 1] ', ' 4[ 01 0] ', '-1        ', '-2[ 10 0] ', &
         '-2[ 01 0] ', '-2[ 00 1] ', '-3[-1-1-1]', '-3[ 11-1] ', '-3[-11 1] ', '-3[ 1-11] ', &
         '-3[ 11 1] ', '-3[-11-1] ', '-3[-1-11] ', '-3[ 1-1-1]', '-2[-11 0] ', '-4[ 00 1] ', &
         '-4[ 00-1] ', '-2[ 11 0] ', '-2[ 0-11] ', '-2[ 01 1] ', '-4[ 10 0] ', '-4[-10 0] ', &
         '-2[-10 1] ', '-4[ 0-10] ', '-2[ 10 1] ', '-4[ 01 0] '/)
      CHARACTER(len=11), DIMENSION(24), PARAMETER :: rname_hexai = (/' 1         ', ' 6[ 00  1] ', &
         ' 3[ 00  1] ', ' 2[ 00  1] ', ' 3[ 00 -1] ', ' 6[ 00 -1] ', ' 2[ 01  0] ', ' 2[-11  0] ', &
         ' 2[ 10  0] ', ' 2[ 21  0] ', ' 2[ 11  0] ', ' 2[ 12  0] ', '-1         ', '-6[ 00  1] ', &
         '-3[ 00  1] ', '-2[ 00  1] ', '-3[ 00 -1] ', '-6[ 00 -1] ', '-2[ 01  0] ', '-2[-11  0] ', &
         '-2[ 10  0] ', '-2[ 21  0] ', '-2[ 11  0] ', '-2[ 12  0] '/)
      CHARACTER(len=12), DIMENSION(7), PARAMETER :: icst = (/'TRICLINIC   ', 'MONOCLINIC  ', &
         'ORTHORHOMBIC', 'TETRAGONAL  ', 'CUBIC       ', 'TRIGONAL    ', 'HEXAGONAL   '/)

      INTEGER :: i, ib(48), ib0(48), ihc, ihc0, ihg, ihg0, indpg, indpg0, invadd, istrin, iswght, &
         isy, isy0, itype, j, k, l, li, li0, lmax, n, nc, nc0, ntot, ntvec0
      INTEGER, DIMENSION(49, 1)                          :: f00
      REAL(KIND=dp) :: a01(3), a02(3), a03(3), b01(3), b02(3), b03(3), b1(3), b2(3), b3(3), &
         dtotstr, origin(3), origin0(3), proj1, proj2, proj3, r(3, 3, 48), r0(3, 3, 48), totstr, &
         tvec0(3, 1), volum, vv0(3)
      REAL(KIND=dp), DIMENSION(3, 1)                     :: x0
      REAL(KIND=dp), DIMENSION(3, 48)                    :: v, v0

      f00 = 0
      x0 = 0._dp
      v = 0._dp
      v0 = 0._dp
! ==--------------------------------------------------------------==
! READ IN LATTICE STRUCTURE
! ==--------------------------------------------------------------==
      DO i = 1, 3
         a01(i) = a1(i)/alat
         a02(i) = a2(i)/alat
         a03(i) = a3(i)/alat
      END DO
      WRITE (iout, '(" KPSYM| NUMBER OF ATOMS (STRUCT):",I6)') nat
      IF (iout > 0) &
         WRITE (iout, '(" KPSYM|",10X,"K   TYPE",14X,"X(K)")')
      itype = 0
      DO i = 1, nat
         ! Assign an atomic type (for internal purposes)
         IF (i .NE. 1) THEN
            DO j = 1, (i - 1)
               IF (ty(j) .EQ. ty(i)) THEN
                  ! Type located
                  GOTO 178
               END IF
            END DO
            ! New type
         END IF
         itype = itype + 1
         IF (itype .GT. nsp) THEN
            IF (iout > 0) &
               WRITE (iout, '(A,I4,")")') &
               ' KPSYM| NUMBER OF ATOMIC TYPES EXCEEDS DIMENSION (NSP=)', &
               nsp
            IF (iout > 0) &
               WRITE (iout, '(" KPSYM| THE ARRAY TY IS:",/,9(1X,10I7,/))') &
               (ty(j), j=1, nat)
            CALL stopgm('K290', 'FATAL ERROR')
         END IF
178      CONTINUE
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",6X,I5,I6,3F10.5)') &
            i, ty(i), (xkapa(j, i), j=1, 3)
      END DO
      ! ==--------------------------------------------------------------==
      ! IS THE STRAIN SIGNIFICANT ?
      ! ==--------------------------------------------------------------==
      dtotstr = delta*delta
      totstr = 0._dp
      istrin = 0
      DO i = 1, 6
         totstr = totstr + ABS(strain(i))
      END DO
      IF (totstr .GT. dtotstr) istrin = 1
      ! ==--------------------------------------------------------------==
      ! Volume of the cell.
      volum = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) + &
              a3(1)*a1(2)*a2(3) - a1(3)*a2(2)*a3(1) - &
              A2(3)*A3(2)*A1(1) - A3(3)*A1(2)*A2(1)
      volum = ABS(volum)
      b1(1) = (a2(2)*a3(3) - a2(3)*a3(2))/volum
      b1(2) = (a2(3)*a3(1) - a2(1)*a3(3))/volum
      b1(3) = (a2(1)*a3(2) - a2(2)*a3(1))/volum
      b2(1) = (a3(2)*a1(3) - a3(3)*a1(2))/volum
      b2(2) = (a3(3)*a1(1) - a3(1)*a1(3))/volum
      b2(3) = (a3(1)*a1(2) - a3(2)*a1(1))/volum
      b3(1) = (a1(2)*a2(3) - a1(3)*a2(2))/volum
      b3(2) = (a1(3)*a2(1) - a1(1)*a2(3))/volum
      b3(3) = (a1(1)*a2(2) - a1(2)*a2(1))/volum
      ! ==--------------------------------------------------------------==
      DO i = 1, 3
         b01(i) = b1(i)*alat
         b02(i) = b2(i)*alat
         b03(i) = b3(i)*alat
      END DO
      ! ==--------------------------------------------------------------==
      ! == GROUP-THEORY ANALYSIS OF LATTICE                             ==
      ! ==--------------------------------------------------------------==
      CALL group1s(iout, a1, a2, a3, nat, ty, xkapa, b1, b2, b3, &
                   ihg, ihc, isy, li, nc, indpg, ib, ntvec, &
                   v, f0, r, tvec, origin, rx, isc, delta)
      ! ==--------------------------------------------------------------==
      DO n = nc + 1, 48
         ib(n) = 0
      END DO
      ! ==--------------------------------------------------------------==
      invadd = 0
      IF (li .EQ. 0) THEN
         IF (iout > 0) &
            WRITE (iout, '(A,/,A,/,A)') &
            ' KPSYM| ALTHOUGH THE POINT GROUP OF THE CRYSTAL DOES NOT', &
            ' KPSYM| CONTAIN INVERSION, THE SPECIAL POINT GENERATION ALGORITHM', &
            ' KPSYM| WILL CONSIDER IT AS A SYMMETRY OPERATION'
         invadd = 1
      END IF
      ! ==--------------------------------------------------------------==
      ! == CRYSTALLOGRAPHIC DATA                                        ==
      ! ==--------------------------------------------------------------==
      IF (iout > 0) THEN
         WRITE (iout, '(/," KPSYM| CRYSTALLOGRAPHIC DATA:")')
         WRITE (iout, '(4X,"A1",3F10.5,10X,"B1",3F10.5)') a1, b1
         WRITE (iout, '(4X,"A2",3F10.5,10X,"B2",3F10.5)') a2, b2
         WRITE (iout, '(4X,"A3",3F10.5,10X,"B3",3F10.5)') a3, b3
      END IF
      ! ==--------------------------------------------------------------==
      ! == GROUP-THEORETICAL INFORMATION                                ==
      ! ==--------------------------------------------------------------==
      IF (iout > 0) &
         WRITE (iout, '(/," KPSYM| GROUP-THEORETICAL INFORMATION:")')
      ! IHG .... Point group of the primitive lattice, holohedral
      IF (iout > 0) &
         WRITE (iout, &
                '(" KPSYM|    POINT GROUP OF THE PRIMITIVE LATTICE: ",A," SYSTEM")') &
         icst(ihg)
      ! IHC .... Code distinguishing between hexagonal and cubic groups
      ! ISY .... Code indicating whether the space group is symmorphic
      IF (isy .EQ. 0) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",4X,"NONSYMMORPHIC GROUP")')
      ELSEIF (isy .EQ. 1) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",4X,"SYMMORPHIC GROUP")')
      ELSEIF (isy .EQ. -1) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",4X,"SYMMORPHIC GROUP WITH NON-STANDARD ORIGIN")')
      ELSEIF (isy .EQ. -2) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",4X,"NONSYMMORPHIC GROUP???")')
      END IF
      ! LI ..... Inversions symmetry
      IF (li .EQ. 0) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",4X,"NO INVERSION SYMMETRY")')
      ELSEIF (li .GT. 0) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM|",4X,"INVERSION SYMMETRY")')
      END IF
      ! NC ..... Total number of elements in the point group
      IF (iout > 0) &
         WRITE (iout, &
                '(" KPSYM|",4X,"TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP:",I3)') nc
      IF (iout > 0) &
         WRITE (iout, '(" KPSYM|",4X,"TO SUM UP: (",I1,5I3,")")') &
         ihg, ihc, isy, li, nc, indpg
      ! IB ..... List of the rotations constituting the point group
      IF (iout > 0) &
         WRITE (iout, '(/," KPSYM|",4X,"LIST OF THE ROTATIONS:")')
      IF (iout > 0) &
         WRITE (iout, '(7X,12I4)') (ib(i), i=1, nc)
      ! V ...... Nonprimitive translations (for nonsymmorphic groups)
      IF (isy .LE. 0) THEN
         IF (iout > 0) &
            WRITE (iout, '(/," KPSYM|",4X,"NONPRIMITIVE TRANSLATIONS:")')
         IF (iout > 0) &
            WRITE (iout, '(A,A)') &
            ' ROT    V  IN THE BASIS A1, A2, A3      ', &
            'V  IN CARTESIAN COORDINATES'
         ! Cartesian components of nonprimitive translation.
         DO i = 1, nc
            DO j = 1, 3
               vv0(j) = v(1, i)*a1(j) + v(2, i)*a2(j) + v(3, i)*a3(j)
            END DO
            IF (iout > 0) &
               WRITE (iout, '(1X,I3,3F10.5,3X,3F10.5)') &
               ib(i), (v(j, i), j=1, 3), vv0
         END DO
      END IF
      ! F0 ..... The function defined in Maradudin, Ipatova by
      ! eq. (3.2.12): atom transformation table.
      IF (iout > 0) &
         WRITE (iout, &
                '(/," KPSYM|",4X,"ATOM TRANSFORMATION TABLE (MARADUDIN,VOSKO):")')
      IF (iout > 0) &
         WRITE (iout, '(5(4X,"R  AT->AT"))')
      IF (iout > 0) &
         WRITE (iout, '(I5," [Identity]")') 1
      DO k = 2, nc
         DO j = 1, nat
            IF (iout > 0) &
               WRITE (iout, '(I5,2I4)', advance="no") ib(k), j, f0(k, j)
            IF ((MOD(j, 5) .EQ. 0) .AND. iout > 0) &
               WRITE (iout, *)
         END DO
         IF ((MOD(j - 1, 5) .NE. 0) .AND. iout > 0) &
            WRITE (iout, *)
      END DO
      ! R ...... List of the 3 x 3 rotation matrices
      IF (iout > 0) &
         WRITE (iout, '(/," KPSYM|",4X,"LIST OF THE 3 X 3 ROTATION MATRICES:")')
      IF (ihc .EQ. 0) THEN
         DO k = 1, nc
            IF (iout > 0) &
               WRITE (iout, &
                      '(4X,I3," (",I2,": ",A11,")",2(3F14.6,/,25X),3F14.6)') &
               k, ib(k), rname_hexai(ib(k)), ((r(i, j, ib(k)), j=1, 3), i=1, 3)
         END DO
      ELSE
         DO k = 1, nc
            IF (iout > 0) &
               WRITE (iout, &
                      '(4X,I3," (",I2,": ",A10,") ",2(3F14.6,/,25X),3F14.6)') &
               k, ib(k), rname_cubic(ib(k)), ((r(i, j, ib(k)), j=1, 3), i=1, 3)
         END DO
      END IF
      ! ==--------------------------------------------------------------==
      ! GENERATE THE BRAVAIS LATTICE
      ! ==--------------------------------------------------------------==
      CALL group1s(iout, a01, a02, a03, 1, ty, x0, b01, b02, b03, &
                   ihg0, ihc0, isy0, li0, nc0, indpg0, ib0, ntvec0, &
                   v0, f00, r0, tvec0, origin0, rx, isc, delta)
      ! ==--------------------------------------------------------------==
      ! It is assumed that the same 'type' of symmetry operations
      ! (cubic/hexagonal) will apply to the crystal as well as the Bravais
      ! lattice.
      ! ==--------------------------------------------------------------==
      IF (iout > 0) &
         WRITE (iout, '(/,1X,19("*"),A,25("*"))') &
         ' GENERATION OF SPECIAL POINTS '
      ! Parameter Q of Monkhorst and Pack, generalized for 3 axes B1,2,3
      IF (iout > 0) &
         WRITE (iout, '(A,/,1X,3I5)') &
         ' KPSYM| MONKHORST-PACK PARAMETERS (GENERALIZED) IQ1,IQ2,IQ3:', &
         iq1, iq2, iq3
      ! WVK0 is the shift of the whole mesh (see Macdonald)
      IF (iout > 0) &
         WRITE (iout, '(A,/,1X,3F10.5)') &
         ' KPSYM| CONSTANT VECTOR SHIFT (MACDONALD) OF THIS MESH:', wvk0
      IF (iabs(iq1) + iabs(iq2) + iabs(iq3) .EQ. 0) GOTO 710
      IF (ABS(istriz) .NE. 1) THEN
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM| INVALID SWITCH FOR SYMMETRIZATION",I10)') istriz
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM| INVALID SWITCH FOR SYMMETRIZATION",I10)') istriz
         CALL stopgm('K290', 'ISTRIZ WRONG ARGUMENT')
      END IF
      IF (iout > 0) &
         WRITE (iout, '(" KPSYM| SYMMETRIZATION SWITCH: ",I3)', advance="no") istriz
      IF (istriz .EQ. 1) THEN
         IF (iout > 0) &
            WRITE (iout, '(" (SYMMETRIZATION OF MONKHORST-PACK MESH)")')
      ELSE
         IF (iout > 0) &
            WRITE (iout, '(" (NO SYMMETRIZATION OF MONKHORST-PACK MESH)")')
      END IF
      ! Set to 0.
      DO i = 1, nkpoint
         lwght(i) = 0
      END DO
      ! ==--------------------------------------------------------------==
      ! == Generation of the points (they are not multiplied            ==
      ! == by 2*Pi because B1,2,3  were not,either)                     ==
      ! ==--------------------------------------------------------------==
      IF (nc .GT. nc0) THEN
         ! Due to non-use of primitive cell, the crystal has more
         ! rotations than Bravais lattice.
         ! We use only the rotations for Bravais lattices
         IF (ntvec .EQ. 1) THEN
            IF (iout > 0) &
               WRITE (iout, *) ' KPSYM| NUMBER OF ROTATIONS FOR BRAVAIS LATTICE', nc0
            IF (iout > 0) &
               WRITE (iout, *) ' KPSYM| NUMBER OF ROTATIONS FOR CRYSTAL LATTICE', nc
            IF (iout > 0) &
               WRITE (iout, *) ' KPSYM| NO DUPLICATION FOUND'
            CALL stopgm('ERROR', &
                        'SOMETHING IS WRONG IN GROUP DETERMINATION')
         END IF
         nc = nc0
         DO i = 1, nc0
            ib(i) = ib0(i)
         END DO
         IF (iout > 0) &
            WRITE (iout, '(/,1X,20("! "),"WARNING",20("!"))')
         IF (iout > 0) &
            WRITE (iout, '(A)') &
            ' KPSYM| THE CRYSTAL HAS MORE SYMMETRY THAN THE BRAVAIS LATTICE'
         IF (iout > 0) &
            WRITE (iout, '(A)') &
            ' KPSYM| BECAUSE THIS IS NOT A PRIMITIVE CELL'
         IF (iout > 0) &
            WRITE (iout, '(A)') &
            ' KPSYM| USE ONLY SYMMETRY FROM BRAVAIS LATTICE'
         IF (iout > 0) &
            WRITE (iout, '(1X,20("! "),"WARNING",20("!"),/)')
      END IF
      CALL sppt2(iout, iq1, iq2, iq3, wvk0, nkpoint, &
                 a01, a02, a03, b01, b02, b03, &
                 invadd, nc, ib, r, ntot, wvkl, lwght, lrot, nc0, ib0, istriz, &
                 nhash, includ, list, rlist, delta)
      ! ==--------------------------------------------------------------==
      ! == Check on error signals                                       ==
      ! ==--------------------------------------------------------------==
      IF (iout > 0) &
         WRITE (iout, '(/," KPSYM|",1X,I5," SPECIAL POINTS GENERATED")') ntot
      IF (ntot .EQ. 0) THEN
         GOTO 710
      ELSE IF (ntot .LT. 0) THEN
         IF (iout > 0) &
            WRITE (iout, '(A,I5,/,A,/,A)') ' KPSYM| DIMENSION NKPOINT =', nkpoint, &
            ' KPSYM| INSUFFICIENT FOR ACCOMMODATING ALL THE SPECIAL POINTS', &
            ' KPSYM| WHAT FOLLOWS IS AN INCOMPLETE LIST'
         ntot = iabs(ntot)
      END IF
      ! Before using the list WVKL as wave vectors, they have to be
      ! multiplied by 2*Pi
      ! The list of weights LWGHT is not normalized
      iswght = 0
      DO i = 1, ntot
         iswght = iswght + lwght(i)
      END DO
      IF (iout > 0) &
         WRITE (iout, '(8X,A,T33,A,4X,A)') &
         'WAVEVECTOR K', 'WEIGHT', 'UNFOLDING ROTATIONS'
      ! Set near-zeroes equal to zero:
      DO l = 1, ntot
         DO i = 1, 3
            IF (ABS(wvkl(i, l)) .LT. delta) wvkl(i, l) = 0._dp
         END DO
         IF (istrin .NE. 0) THEN
            ! Express special points in (unstrained) basis.
            proj1 = 0._dp
            proj2 = 0._dp
            proj3 = 0._dp
            DO i = 1, 3
               proj1 = proj1 + wvkl(i, l)*a01(i)
               proj2 = proj2 + wvkl(i, l)*a02(i)
               proj3 = proj3 + wvkl(i, l)*a03(i)
            END DO
            DO i = 1, 3
               wvkl(i, l) = proj1*b1(i) + proj2*b2(i) + proj3*b3(i)
            END DO
         END IF
         lmax = lwght(l)
         IF (iout > 0) &
            WRITE (iout, fmt='(1X,I5,3F8.4,I8,T42,12I3)') &
            l, (wvkl(i, l), i=1, 3), lwght(l), (lrot(i, l), i=1, MIN(lmax, 12))
         DO j = 13, lmax, 12
            IF (iout > 0) &
               WRITE (iout, fmt='(T42,12I3)') &
               (lrot(i, l), i=j, MIN(lmax, j - 1 + 12))
         END DO
      END DO
      IF (iout > 0) &
         WRITE (iout, '(24X,"TOTAL:",I8)') iswght
      ! ==--------------------------------------------------------------==
710   CONTINUE
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE k290s
! **************************************************************************************************

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
      ! ==          IF IOUT.LE.0 NO MESSAGE                             ==
      ! == A1,A2,A3 .. ELEMENTARY TRANSLATIONS OF THE LATTICE, IN SOME  ==
      ! ==          UNIT OF LENGTH                                      ==
      ! == NAT .... NUMBER OF ATOMS IN THE UNIT CELL                    ==
      ! ==          ALL THE DIMENSIONS ARE SET FOR NAT .LE. 20          ==
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
      ! ==          LI.GT.0 MEANS: THERE IS INVERSION IN THE POINT      ==
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
      IF (ntvec .GT. 1) THEN
         ! All rotations found by PGL1 have axes in x, y or z cart. axis
         ! So we have too check if we do not loose symmetry
         ncprim = nc
         ! The hexagonal system is found if the z axis is the sixfold axis
         CALL pgl1(a, ai, ihc, nc, ib, ihg, r, delta)
         IF (ncprim .GT. nc) THEN
            ! More symmetry with
            CALL pgl1(ap, api, ihc, nc, ib, ihg, r, delta)
         END IF
      END IF

      ! Determination of the space group
      CALL atftm1(iout, r, v, x, f0, origin, ib, ty, nat, ihg, ihc, rx, &
                  nc, indpg, ntvec, a, ai, li, isy, isc, delta)

      IF (iout .GT. 0) THEN
         IF (li .GT. 0) THEN
            IF (iout > 0) &
               WRITE (iout, '(1X,A)') &
               'KPSYM| THE POINT GROUP OF THE CRYSTAL CONTAINS THE INVERSION'
         END IF
         IF (iout > 0) &
            WRITE (iout, *)
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
         IF (i .EQ. 1) il = 2
         IF (i .EQ. 3) iu = 2
         DO j = 1, 3
            jl = 1
            ju = 3
            IF (j .EQ. 1) jl = 2
            IF (j .EQ. 3) ju = 2
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
         IF (ty(1) .NE. ty(k2)) GOTO 100
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
               IF (f0(49, i) .GT. f0(1, i)) f0(49, i) = f0(1, i)
            END DO
            DO i = 1, 3
               tvec(i, ntvec) = vr(i)
            END DO
         END IF
100      CONTINUE
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
      IF (ntvec .EQ. 1) THEN
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
               IF (ABS(vr(i)) .GT. delta) THEN
                  il = NINT(1._dp/ABS(vr(i)))
                  IF (il .GT. 1) THEN
                     ! We replace AP(1:3,I) by TVEC(1:3,IV)
                     DO j = 1, 3
                        ap(j, i) = xb(j)
                     END DO
                     ! Calculate new API
                     CALL calbrec(ap, api)
                     GOTO 200  ! EXIT
                  END IF
               END IF
            END DO
200         CONTINUE
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
         IF (ihc .EQ. 0) THEN
            nr = 24
         ELSE
            nr = 48
         END IF
         nc = 0
         ! Constructs rotation operations.
         CALL rot1(ihc, r)
         DO n = 1, nr
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
               IF (tr .GT. delta) GOTO 140
            END DO
            nc = nc + 1
            ib(nc) = n
140         CONTINUE
         END DO
         ! ==------------------------------------------------------------==
         ! IHG stands for holohedral group number.
         IF (ihc .EQ. 0) THEN
            ! Hexagonal group:
            IF (nc .EQ. 12) ihg = 6
            IF (nc .GT. 12) ihg = 7
            IF (nc .GE. 12) RETURN
            ! Too few operations, try cubic group: (IHC=1,NR=48)
         ELSE
            ! Cubic group:
            IF (nc .LT. 4) ihg = 1
            IF (nc .EQ. 4) ihg = 2
            IF (nc .GT. 4) ihg = 3
            IF (nc .EQ. 16) ihg = 4
            IF (nc .GT. 16) ihg = 5
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
      IF (ts .LE. delta) RETURN
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
      ! == SUBROUTINES NEEDED: RLV3 CHECKRLV3 SYMMORPHIC STOPGM XSTRING ==
      ! == MAY 14TH,1998: A LOT OF CHANGES (ARGUMENTS)                  ==
      ! ==                BETTER DETERMINATION OF V                     ==
      ! == SEP 15TH,1998: DETERMINATION OF FRACTIONAL TRANSLATIONAL VEC.==
      ! ==--------------------------------------------------------------==
      ! == INPUT:                                                       ==
      ! ==   IOUT Logical file number (output)                          ==
      ! ==        If IOUT.LE.0 no message                               ==
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
      ! == DELTA  REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)             ==
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
      ! ==          LI.GT.0 : there is inversion in the point           ==
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

      CHARACTER(len=10), DIMENSION(48), PARAMETER :: rname_cubic = (/' 1        ', ' 2[ 10 0] ', &
         ' 2[ 01 0] ', ' 2[ 00 1] ', ' 3[-1-1-1]', ' 3[ 11-1] ', ' 3[-11 1] ', ' 3[ 1-11] ', &
         ' 3[ 11 1] ', ' 3[-11-1] ', ' 3[-1-11] ', ' 3[ 1-1-1]', ' 2[-11 0] ', ' 4[ 00 1] ', &
         ' 4[ 00-1] ', ' 2[ 11 0] ', ' 2[ 0-11] ', ' 2[ 01 1] ', ' 4[ 10 0] ', ' 4[-10 0] ', &
         ' 2[-10 1] ', ' 4[ 0-10] ', ' 2[ 10 1] ', ' 4[ 01 0] ', '-1        ', '-2[ 10 0] ', &
         '-2[ 01 0] ', '-2[ 00 1] ', '-3[-1-1-1]', '-3[ 11-1] ', '-3[-11 1] ', '-3[ 1-11] ', &
         '-3[ 11 1] ', '-3[-11-1] ', '-3[-1-11] ', '-3[ 1-1-1]', '-2[-11 0] ', '-4[ 00 1] ', &
         '-4[ 00-1] ', '-2[ 11 0] ', '-2[ 0-11] ', '-2[ 01 1] ', '-4[ 10 0] ', '-4[-10 0] ', &
         '-2[-10 1] ', '-4[ 0-10] ', '-2[ 10 1] ', '-4[ 01 0] '/)
      CHARACTER(len=11), DIMENSION(24), PARAMETER :: rname_hexai = (/' 1         ', ' 6[ 00  1] ', &
         ' 3[ 00  1] ', ' 2[ 00  1] ', ' 3[ 00 -1] ', ' 6[ 00 -1] ', ' 2[ 01  0] ', ' 2[-11  0] ', &
         ' 2[ 10  0] ', ' 2[ 21  0] ', ' 2[ 11  0] ', ' 2[ 12  0] ', '-1         ', '-6[ 00  1] ', &
         '-3[ 00  1] ', '-2[ 00  1] ', '-3[ 00 -1] ', '-6[ 00 -1] ', '-2[ 01  0] ', '-2[-11  0] ', &
         '-2[ 10  0] ', '-2[ 21  0] ', '-2[ 11  0] ', '-2[ 12  0] '/)
      CHARACTER(len=12), DIMENSION(7), PARAMETER :: icst = (/'TRICLINIC   ', 'MONOCLINIC  ', &
         'ORTHORHOMBIC', 'TETRAGONAL  ', 'CUBIC       ', 'TRIGONAL    ', 'HEXAGONAL   '/)
      CHARACTER(len=3), DIMENSION(32), PARAMETER :: pgrd = (/'c1 ', 'ci ', 'c2 ', 'c1h', 'c2h', &
         'c3 ', 'c3i', 'd3 ', 'c3v', 'd3 ', 'c4 ', 's4 ', 'c4h', 'd4 ', 'c4v', 'd2d', 'd4h', 'c6 ',&
         'c3h', 'c6h', 'd6 ', 'c6v', 'd3h', 'd6h', 'd2 ', 'c2v', 'd2h', 't  ', 'th ', 'o  ', 'td ',&
         'oh '/)
      CHARACTER(len=5), DIMENSION(32), PARAMETER :: pgrp = (/'    1', '  <1>', '    2', '    m', &
         '  2/m', '    3', '  <3>', '   32', '   3m', ' <3>m', '    4', '  <4>', '  4/m', '  422', &
         '  4mm', '<4>2m', '4/mmm', '    6', '  <6>', '  6/m', '  622', '  6mm', '<6>m2', '6/mmm', &
         '  222', '  mm2', '  mmm', '   23', '   m3', '  432', '<4>3m', '  m3m'/)

      INTEGER                                            :: i, iis(48), il, info, j, k, k2, l, n, &
                                                            nca, ni
      LOGICAL                                            :: nodupli, oksym
      REAL(dp)                                           :: vc(3, 48), vr(3), vs, xb(3)

      nodupli = ntvec .EQ. 1
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
         IF (oksym) THEN
            GOTO 190
         END IF
         ! Now we try other possible VR
         ! F0(49,1:NAT) has only inequivalent atom indexes for translation
         DO k2 = 1, nat
            IF (f0(49, k2) .LT. k2) GOTO 185
            IF (ty(1) .NE. ty(k2)) GOTO 185
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
            IF (oksym) THEN
               GOTO 190
            END IF
185         CONTINUE
         END DO
         iis(l) = 0
         GOTO 210
190      CONTINUE
         nca = nca + 1
         DO i = 1, 3
            v(i, nca) = vr(i)
         END DO
         ! ==------------------------------------------------------------==
         ! == V(I,N) IS THE I-TH COMPONENT OF THE FRACTIONAL             ==
         ! == TRANSLATION ASSOCIATED WITH THE ROTATION N.                ==
         ! == ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS, THEY ARE     ==
         ! == GIVEN IN THE SYSTEM A1,A2,A3.                              ==
         ! == K.K., 23.10. 1979                                         ==
         ! ==------------------------------------------------------------==
210      CONTINUE
      END DO
      ! Remove unused operations
      i = 0
      ni = 13
      IF (ihg .LT. 6) ni = 25
      li = 0
      DO n = 1, nc
         l = ib(n)
         IF (iis(l) .EQ. 0) GOTO 230 ! CYCLE
         i = i + 1
         ib(i) = ib(n)
         IF (ib(i) .EQ. ni) li = i
         DO k = 1, nat
            f0(i, k) = f0(n, k)
         END DO
230      CONTINUE
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
      IF (vs .GT. delta) THEN
         isy = 0
      ELSE
         isy = 1
      END IF
      ! ==--------------------------------------------------------------==
      ! Determination of the point group
      ! (Thierry Deutsch - 1998 [Maybe not complete!!])
      IF (ihg .LT. 6) THEN
         IF (nc .EQ. 0) THEN
            IF (iout > 0) &
               WRITE (iout, '(" ATFTM1! IHG=",A," NC=",I2)') icst(ihg), nC
            CALL stopgm('ATFTM1', 'NUMBER OF ROTATION NULL')
            ! Triclinic system
         ELSEIF (nc .EQ. 1) THEN
            ! IB=1
            indpg = 1           ! 1 (c1)
         ELSEIF (nc .EQ. 2 .AND. ib(2) .EQ. 25) THEN
            ! IB=125
            indpg = 2           ! <1>(ci)
         ELSEIF (nc .EQ. 2 .AND. ( &
                 ib(2) .EQ. 4 .OR. & ! 2[001]
                 ib(2) .EQ. 2 .OR. & ! 2[100]
                 ib(2) .EQ. 3)) THEN ! 2[010]
            ! Monoclinic system
            ! IB=14 (z-axis) OR
            ! IB=12 (x-axis) OR
            ! IB=13 (y-axis)
            indpg = 3           ! 2 (c2)
         ELSEIF (nc .EQ. 2 .AND. ( &
                 ib(2) .EQ. 28 .OR. &
                 ib(2) .EQ. 26 .OR. &
                 ib(2) .EQ. 27)) THEN
            ! IB=128 (z-axis) OR
            ! IB=126 (x-axis) OR
            ! IB=127 (y-axis)
            indpg = 4           ! m (c1h)
         ELSEIF (nc .EQ. 4 .AND. ( &
                 ib(4) .EQ. 28 .OR. & ! 2[001]
                 ib(4) .EQ. 27 .OR. & ! 2[010]
                 ib(4) .EQ. 26 .OR. & ! 2[100]
                 ib(4) .EQ. 37 .OR. & ! -2[-110]
                 ib(4) .EQ. 40)) THEN ! 2[110]
            ! IB=1  425 28 (z-axis)  OR
            ! IB=1  225 26 (x-axis)  OR
            ! IB=1  325 27 (y-axis)  OR
            ! IB=113 2537 (-xy-axis)OR
            ! IB=116 2540 (xy-axis)
            indpg = 5           ! 2/m(c2h)
         ELSEIF (nc .EQ. 4 .AND. ( &
                 ib(4) .EQ. 15 .OR. &
                 ib(4) .EQ. 20 .OR. &
                 ib(4) .EQ. 24)) THEN
            ! Tetragonal system
            ! IB=14 1415 (z-axis) OR
            ! IB=12 1920 (x-axis) OR
            ! IB=13 2224 (y-axis)
            indpg = 11          ! 4 (c4)
         ELSEIF (nc .EQ. 4 .AND. ( &
                 ib(4) .EQ. 39 .OR. &
                 ib(4) .EQ. 44 .OR. &
                 ib(4) .EQ. 48)) THEN
            ! IB=14 3839 (z-axis) OR
            ! IB=12 4344 (x-axis) OR
            ! IB=13 4648 (y-axis)
            indpg = 12          ! <4>(s4)
         ELSEIF (nc .EQ. 8 .AND. ( &
                 (ib(3) .EQ. 14 .AND. ib(8) .EQ. 39) .OR. &
                 (ib(3) .EQ. 19 .AND. ib(8) .EQ. 44) .OR. &
                 (ib(3) .EQ. 22 .AND. ib(8) .EQ. 48))) THEN
            ! IB=14 1415 2825 3839 (z-axis) OR
            ! IB=12 1920 2625 4344 (x-axis) OR
            ! IB=13 2224 2725 4648 (y-axis)
            indpg = 13          ! 422(d4)
         ELSEIF (nc .EQ. 8 .AND. ib(4) .EQ. 4 .AND. ( &
                 ib(8) .EQ. 16 .OR. &
                 ib(8) .EQ. 20 .OR. &
                 ib(8) .EQ. 24)) THEN
            ! IB=12 3 413 1415 16 (z-axis) OR
            ! IB=12 3 417 1920 18 (x-axis) OR
            ! IB=12 3 421 2224 23 (y-axis)
            indpg = 14          ! 4/m(c4h)
         ELSEIF (nc .EQ. 8 .AND. ( &
                 ib(8) .EQ. 40 .OR. &
                 ib(8) .EQ. 42 .OR. &
                 ib(8) .EQ. 47)) THEN
            ! IB=14 1415 2627 3740 (z-axis) OR
            ! IB=12 1920 2827 4142 (x-axis) OR
            ! IB=13 2224 2628 4547 (y-axis)
            indpg = 15          ! 4mm(c4v)
         ELSEIF (nc .EQ. 8 .AND. ( &
                 (ib(3) .EQ. 13 .AND. ib(8) .EQ. 39) .OR. &
                 (ib(3) .EQ. 17 .AND. ib(8) .EQ. 44) .OR. &
                 (ib(3) .EQ. 21 .AND. ib(8) .EQ. 48))) THEN
            ! IB=14 1316 2627 3839 (z-axis) OR
            ! IB=12 1718 2827 4344 (x-axis) OR
            ! IB=13 2123 2628 4648 (y-axis)
            indpg = 16          ! <4>2m(d2d)
         ELSEIF (nc .EQ. 16 .AND. ( &
                 ib(16) .EQ. 40 .OR. &
                 ib(16) .EQ. 44 .OR. &
                 ib(16) .EQ. 48)) THEN
            ! IB=12 3 413 1415 1625 2627 2837 3839 40 (z-axis) OR
            ! IB=12 3 417 1920 1825 2627 2841 4344 42 (x-axis) OR
            ! IB=12 3 421 2224 2325 2627 2845 4648 47 (y-axis)
            indpg = 17          ! 4/mmm(d4h)
         ELSEIF (nc .EQ. 4 .AND. (ib(4) .EQ. 4)) THEN
            ! Orthorhombic system
            ! IB=12  3  4
            indpg = 25          ! 222(d2)
         ELSEIF (nc .EQ. 4 .AND. ( &
                 ib(4) .EQ. 27 .OR. &
                 ib(4) .EQ. 28)) THEN
            ! IB=13 2627 (z-axis) OR
            ! IB=12 2728 (x-axis) OR
            ! IB=14 2628 (y-axis) OR
            indpg = 26          ! mm2(c2v)
         ELSEIF (nc .EQ. 8) THEN
            ! IB=12 3 425 2627 28
            indpg = 27          ! mmm(d2h)
         ELSEIF (nc .EQ. 12 .AND. ( &
                 ib(12) .EQ. 12 .OR. &
                 ib(12) .EQ. 47 .OR. &
                 ib(12) .EQ. 45)) THEN
            ! Cubic system
            ! IB=12  3  4  5  6  7  8  910 1112 OR
            ! IB=15 1113 1823 2530 3537 4247 OR
            ! IB=18 1016 1821 2532 3440 4245
            indpg = 28          ! 23 (t)
         ELSEIF (nc .EQ. 24 .AND. ib(24) .EQ. 36) THEN
            ! IB= 1  2  3  4  5  6  7  8  910 1112
            ! 2526 2728 2930 3132 3334 3536
            indpg = 29          ! m3 (th)
         ELSEIF (nc .EQ. 24 .AND. ib(24) .EQ. 24) THEN
            ! IB=12 3 45 6 78 9 1011 12
            ! 1314 1516 1718 1920 2122 2324
            indpg = 30          ! 432 (o)
         ELSEIF (nc .EQ. 24 .AND. ib(24) .EQ. 48) THEN
            ! IB=12 3 45 6 78 9 1011 12
            ! 3738 3940 4142 4345 4647 48
            indpg = 31          ! <4>3m(td)
         ELSEIF (nc .EQ. 48) THEN
            ! IB=1..48
            indpg = 32          ! m3m(oh)
         ELSE
            ! WRITE(6,'(" ATFTM1! IHG=",A," NC=",I2)') ICST(IHG),NC
            ! WRITE(6,'(" ATFTM1!",19I3)') (IB(I),I=1,NC)
            ! WRITE(6,'(" ATFTM1! THIS CASE IS UNKNOWN IN THE DATABASE")')
            ! Probably a sub-group of 32
            indpg = -32
         END IF
      ELSEIF (ihg .GE. 6) THEN
         IF (nc .EQ. 0) THEN
            IF (iout > 0) &
               WRITE (iout, '(" ATFTM1! IHG=",A," NC=",I2)') icst(ihg), nC
            CALL stopgm('ATFTM1', 'NUMBER OF ROTATION NULL')
            ! Triclinic system
         ELSEIF (nc .EQ. 1) THEN
            ! IB=1
            indpg = 1           ! 1 (c1)
         ELSEIF (nc .EQ. 2 .AND. ib(2) .EQ. 13) THEN
            ! IB=113
            indpg = 2           ! <1>(ci)
         ELSEIF (nc .EQ. 2 .AND. ( &
                 ib(2) .EQ. 4)) THEN ! 2[001]
            ! Monoclinic system
            ! IB=1 4
            indpg = 3           ! 2 (c2)
         ELSEIF (nc .EQ. 2 .AND. ( &
                 ib(2) .EQ. 16)) THEN
            ! IB=116
            indpg = 4           ! m (c1h)
         ELSEIF (nc .EQ. 4 .AND. ( &
                 ib(4) .EQ. 24 .OR. &
                 ib(4) .EQ. 20)) THEN
            ! IB=112 1324 OR
            ! IB=1  813 20
            indpg = 5           ! 2/m(c2h)
         ELSEIF (nc .EQ. 3 .AND. ib(3) .EQ. 5) THEN
            ! Trigonal system
            ! IB=13 5
            indpg = 6           ! 3 (c3)
         ELSEIF (nc .EQ. 6 .AND. ib(6) .EQ. 17) THEN
            ! IB=113 1517 35
            indpg = 7           ! <3>(c3i)
         ELSEIF (nc .EQ. 6 .AND. ib(6) .EQ. 11) THEN
            ! IB=17 9 1135
            indpg = 8           ! 32 (d3)
         ELSEIF (nc .EQ. 6 .AND. ib(6) .EQ. 23) THEN
            ! IB=13 5 1921 23
            indpg = 9           ! 3m (c3v)
         ELSEIF (nc .EQ. 12 .AND. ib(12) .EQ. 23) THEN
            ! IB=13 5 79 1113 1517 1921 23
            indpg = 10          ! <3>m(d3d)
         ELSEIF (nc .EQ. 6 .AND. ib(6) .EQ. 6) THEN
            ! Hexagonal system
            ! IB=12 3 45 6
            indpg = 18          ! 6 (c6)
         ELSEIF (nc .EQ. 6 .AND. ib(6) .EQ. 18) THEN
            ! IB=13 5 1416 18
            indpg = 19          ! <6>(c3h)
         ELSEIF (nc .EQ. 12 .AND. ib(12) .EQ. 18) THEN
            ! IB=12 3 45 6 1314 1516 1718
            indpg = 20          ! 6/m(c6h)
         ELSEIF (nc .EQ. 12 .AND. ib(12) .EQ. 12) THEN
            ! IB=12 3 45 6 78 9 1011 12
            indpg = 21          ! 622(d6)
         ELSEIF (nc .EQ. 12 .AND. ib(2) .EQ. 2 .AND. ib(12) .EQ. 24) THEN
            ! IB=12 3 45 6 1920 2122 2324
            indpg = 22          ! 6mm(c6v)
         ELSEIF (nc .EQ. 12 .AND. ib(2) .EQ. 3 .AND. ib(12) .EQ. 24) THEN
            ! IB=13 5 79 1114 1618 2022 24
            indpg = 23          ! <6>m2(d3h)
         ELSEIF (nc .EQ. 24) THEN
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
      IF (isy .NE. 1) THEN
         ! Transform V in cartesian coordinates
         DO n = 1, nc
            vc(1, n) = a(1, 1)*v(1, n) + a(1, 2)*v(2, n) + a(1, 3)*v(3, n)
            vc(2, n) = a(2, 1)*v(1, n) + a(2, 2)*v(2, n) + a(2, 3)*v(3, n)
            vc(3, n) = a(3, 1)*v(1, n) + a(3, 2)*v(2, n) + a(3, 3)*v(3, n)
         END DO
         CALL symmorphic(nc, ib, r, vc, ai, info, origin, delta)
         IF (info .EQ. 1) THEN
            CALL rlv3(ai, origin, xb, il, delta)
            ! !!!RLV3 determines -XB in crystal coordinates
            ! !!We want between 0.0 and 1.0
            DO i = 1, 3
               IF (-xb(i) .GE. 0._dp) THEN
                  origin(i) = -xb(i)
               ELSE
                  origin(i) = 1._dp - xb(i)
               END IF
            END DO
            DO i = 1, 3
               xb(i) = a(i, 1)*origin(1) + a(i, 2)*origin(2) + a(i, 3)*origin(3)
            END DO
            isy = -1
         ELSEIF (info .EQ. 0) THEN
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
      IF (iout .GT. 0) THEN
         IF (iout > 0) &
            WRITE (iout, *)
         CALL xstring(icst(ihg), i, j)
         IF ((ihg .EQ. 7 .AND. nc .EQ. 24) .OR. &
             (ihg .EQ. 5 .AND. nc .EQ. 48)) THEN
            IF (iout > 0) &
               WRITE (iout, '(A,A,A)') &
               ' KPSYM| THE POINT GROUP OF THE CRYSTAL IS THE FULL ', &
               icst(ihg) (i:j), &
               ' GROUP'
         ELSE
            IF (iout > 0) &
               WRITE (iout, '(A,A,A,I2,A)') &
               ' KPSYM| THE CRYSTAL SYSTEM IS ', &
               icst(ihg) (i:j), &
               ' WITH ', nc, ' OPERATIONS:'
            IF (ihc .EQ. 0) THEN
               IF (iout > 0) &
                  WRITE (iout, '( 5(5(A13),/))') (rname_hexai(ib(i)), i=1, nc)
            ELSE
               IF (iout > 0) &
                  WRITE (iout, '(10(5(A13),/))') (rname_cubic(ib(i)), i=1, nc)
            END IF
         END IF
         ! ==------------------------------------------------------------==
         IF (isy .EQ. 1) THEN
            IF (iout > 0) &
               WRITE (iout, '(A)') &
               ' KPSYM| THE SPACE GROUP OF THE CRYSTAL IS SYMMORPHIC'
         ELSEIF (isy .EQ. -1) THEN
            IF (iout > 0) &
               WRITE (iout, '(A)') &
               ' KPSYM| THE SPACE GROUP OF THE CRYSTAL IS SYMMORPHIC'
            IF (iout > 0) &
               WRITE (iout, '(A,A,/,T3,3F10.6,3X,3F10.6)') &
               ' KPSYM| THE STANDARD ORIGIN OF COORDINATES IS:   ', &
               '[CARTESIAN]   [CRYSTAL]', xb, origin
         ELSEIF (isy .EQ. 0) THEN
            IF (iout > 0) &
               WRITE (iout, '(A,/,3X,A,F15.6,A)') &
               ' KPSYM| THE SPACE GROUP IS NON-SYMMORPHIC,', &
               ' (SUM OF TRANSLATION VECTORS=', vs, ')'
         ELSEIF (isy .EQ. -2) THEN
            IF (iout > 0) &
               WRITE (iout, '(A,A)') &
               ' KPSYM| CANNOT DETERMINE IF THE SPACE GROUP IS', &
               ' SYMMORPHIC OR NOT'
            IF (iout > 0) &
               WRITE (iout, '(A,/,A,/,3X,A,F15.6,A)') &
               ' KPSYM| THE SPACE GROUP IS NON-SYMMORPHIC,', &
               ' KPSYM| OR ELSE A NON STANDARD ORIGIN OF COORDINATES WAS USED.', &
               ' KPSYM| (SUM OF TRANSLATION VECTORS=', vs, ')'
         END IF
         IF (indpg .GT. 0) THEN
            CALL xstring(pgrp(indpg), i, j)
            CALL xstring(pgrd(indpg), k, l)
            IF (iout > 0) &
               WRITE (iout, '(A,A,"(",A,")",T56,"[INDEX=",I2,"]")') &
               ' KPSYM| THE POINT GROUP OF THE CRYSTAL IS  ', pgrp(indpg) (i:j), &
               pgrd(indpg) (k:l), indpg
         ELSE
            CALL xstring(pgrp(-indpg), i, j)
            CALL xstring(pgrd(-indpg), k, l)
            IF (iout > 0) &
               WRITE (iout, '(A,I2,A,A,"(",A,")",T56,"[INDEX=",I2,"]")') &
               ' KPSYM| POINT GROUP: GROUP ORDER=', nc, &
               '    SUBGROUP OF ', pgrp(-indpg) (i:j), &
               pgrd(-indpg) (k:l), -indpg
         END IF
         IF (ntvec .EQ. 1) THEN
            IF (iout > 0) &
               WRITE (iout, '(A,T60,I6)') &
               ' KPSYM| NUMBER OF PRIMITIVE CELL:', ntvec
         ELSE
            IF (iout > 0) &
               WRITE (iout, '(A,T60,I6)') &
               ' KPSYM| NUMBER OF PRIMITIVE CELLS:', ntvec
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
      REAL(dp)                                           :: vt(3), xb(3)

      DO ia = 1, nat
         isc(ia) = 0
      END DO
      ! Now we check if ROT(N)+VR gives a correct symmetry.
      DO ia = 1, nat
         DO ib = 1, nat
            IF (ty(ia) .EQ. ty(ib) .AND. isc(ib) .EQ. 0) THEN
               xb(1) = rx(1, ia) - x(1, ib)
               xb(2) = rx(2, ia) - x(2, ib)
               xb(3) = rx(3, ia) - x(3, ib)
               CALL rlv3(ai, xb, vt, il, delta)
               ! VT STANDS FOR V-TEST
               oksym = (ABS((vr(1) - vt(1)) - NINT(vr(1) - vt(1))) .LT. delta) .AND. &
                       (ABS((vr(2) - vt(2)) - NINT(vr(2) - vt(2))) .LT. delta) .AND. &
                       (ABS((vr(3) - vt(3)) - NINT(vr(3) - vt(3))) .LT. delta)
               IF (oksym) THEN
                  IF (nodupli) isc(ib) = 1
                  f0(n, ia) = ib
                  ! IR+VR is the good one: another symmetry operation
                  ! Next atom
                  GOTO 100
               END IF
            END IF
         END DO
         ! VR is not the correct translation vector
         RETURN
100      CONTINUE
      END DO
      ! ==--------------------------------------------------------------==
      RETURN
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
         IF (dif .GT. delta*delta) THEN
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
            IF (ierror .EQ. 0) THEN
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
                  IF (i .NE. ierror) THEN
                     i1 = i1 + 1
                     j1 = 0
                     DO j = 1, 3
                        IF (j .NE. ierror) THEN
                           j1 = j1 + 1
                           r2(i1, j1) = -r(i, j, ib(ir))
                        END IF
                     END DO
                     r2(i1, i1) = 1 + r2(i1, i1)
                  END IF
               END DO
               CALL invmat(r2, ierror)
               IF (ierror .EQ. 0) THEN
                  ! The matrix 2X2 has an inverse.
                  ! Solve Vxy = (1-R).OAxy + OAz R3z (z is IMISSING3)
                  i1 = 0
                  DO i = 1, 3
                     IF (igood(i) .EQ. 1) THEN
                        i1 = i1 + 1
                        vr(i) = 0._dp
                        j1 = 0
                        DO j = 1, 3
                           IF (igood(j) .EQ. 1) THEN
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
                     IF (i .NE. imissing3) THEN
                        i1 = i1 + 1
                        IF (i1 .EQ. ierror) THEN
                           igood(i) = 0
                           imissing2 = i
                        ELSE
                           ionly = i
                        END IF
                     END IF
                  END DO
                  diag = (1 - r(ionly, ionly, ib(ir)))
                  IF (ABS(diag) .GT. delta) THEN
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
               IF (iok(i) .EQ. 1) THEN
                  dif = dif + ABS(origin(i) - vr(i))
               END IF
            END DO
            IF (dif .GT. delta) THEN
               ! Non-symmorphic
               info = 0
               RETURN
            ELSE
               DO i = 1, 3
                  IF (iok(i) .NE. 1 .AND. igood(i) .EQ. 1) THEN
                     iok(i) = 1
                     origin(i) = vr(i)
                  END IF
               END DO
            END IF
         END IF
      END DO
      ! ==--------------------------------------------------------------==
      IF (iok(1) .EQ. 0 .AND. iok(2) .EQ. 0 .AND. iok(3) .EQ. 0) THEN
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
         IF (dif .GT. delta) THEN
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
      IF (ihc .EQ. 0) THEN
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
! **************************************************************************************************
!> \brief ...
!> \param iout ...
!> \param iq1 ...
!> \param iq2 ...
!> \param iq3 ...
!> \param wvk0 ...
!> \param nkpoint ...
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \param b1 ...
!> \param b2 ...
!> \param b3 ...
!> \param inv ...
!> \param nc ...
!> \param ib ...
!> \param r ...
!> \param ntot ...
!> \param wvkl ...
!> \param lwght ...
!> \param lrot ...
!> \param ncbrav ...
!> \param ibrav ...
!> \param istriz ...
!> \param nhash ...
!> \param includ ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE sppt2(iout, iq1, iq2, iq3, wvk0, nkpoint, &
                    a1, a2, a3, b1, b2, b3, &
                    inv, nc, ib, r, ntot, wvkl, lwght, lrot, &
                    ncbrav, ibrav, istriz, &
                    nhash, includ, list, rlist, delta)
      ! ==--------------------------------------------------------------==
      ! == WRITTEN ON SEPTEMBER 12-20TH, 1979 BY K.K.                   ==
      ! == MODIFIED 26-MAY-82 BY OLE HOLM NIELSEN                       ==
      ! == GENERATION OF SPECIAL POINTS FOR AN ARBITRARY LATTICE,       ==
      ! == FOLLOWING THE METHOD MONKHORST,PACK,                         ==
      ! == PHYS. REV. B13 (1976) 5188                                   ==
      ! == MODIFIED BY MACDONALD, PHYS. REV. B18 (1978) 5897            ==
      ! == THE SUBROUTINE IS WRITTEN ASSUMING THAT THE POINTS ARE       ==
      ! == GENERATED IN THE RECIPROCAL SPACE.                           ==
      ! == IF, HOWEVER, THE B1,B2,B3 ARE REPLACED BY A1,A2,A3, THEN     ==
      ! == SPECIAL POINTS IN THE DIRECT SPACE CAN BE PRODUCED, AS WELL. ==
      ! == (NO MULTIPLICATION BY 2PI IS THEN NECESSARY.)                ==
      ! == IN THE CASE OF NONSYMMORPHIC GROUPS, THE APPLICATION IN THE  ==
      ! == DIRECT SPACE WOULD PROBABLY REQUIRE A CERTAIN CAUTION.       ==
      ! == SUBROUTINES NEEDED: BZDEFI,BZRDUC,INBZ,MESH                  ==
      ! == IN THE CASES WHERE THE POINT GROUP OF THE CRYSTAL DOES NOT   ==
      ! == CONTAIN INVERSION. THE LATTER MAY BE ADDED IF WE WISH        ==
      ! == (SEE COMMENT TO THE SWITCH INV).                             ==
      ! == REDUCTION TO THE 1ST BRILLOUIN ZONE IS DONE                  ==
      ! == BY ADDING G-VECTORS TO FIND THE SHORTEST WAVE-VECTOR.        ==
      ! == THE ROTATIONS OF THE BRAVAIS LATTICE ARE APPLIED TO THE      ==
      ! == MONKHORST/PACK MESH IN ORDER TO FIND ALL K-POINTS            ==
      ! == THAT ARE RELATED BY SYMMETRY. (OLE HOLM NIELSEN)             ==
      ! ==--------------------------------------------------------------==
      ! == INPUT DATA:                                                  ==
      ! == IOUT:    LOGICAL UNIT FOR OUTPUT                             ==
      ! ==          IF (IOUT.LE.0) NO MESSAGE                           ==
      ! == IQ1,IQ2,IQ3 .. PARAMETER Q OF MONKHORST AND PACK,            ==
      ! ==          GENERALIZED AND DIFFERENT FOR THE 3 DIRECTIONS B1,  ==
      ! ==          B2 AND B3                                           ==
      ! == WVK0 ... THE 'ARBITRARY' SHIFT OF THE WHOLE MESH, DENOTED K0 ==
      ! ==          IN MACDONALD. WVK0 = 0 CORRESPONDS TO THE ORIGINAL  ==
      ! ==          SCHEME OF MONKHORST AND PACK.                       ==
      ! ==          UNITS: 2PI/(UNITS OF LENGTH  USED IN A1, A2, A3),   ==
      ! ==          I.E. THE SAME  UNITS AS THE GENERATED SPECIAL POINTS==
      ! == NKPOINT .. VARIABLE DIMENSION OF THE (OUTPUT) ARRAYS WVKL,   ==
      ! ==          LWGHT,LROT, I.E. SPACE RESERVED FOR THE SPECIAL     ==
      ! ==          POINTS AND ACCESSORIES.                             ==
      ! ==          NKPOINT HAS TO BE .GE. NTOT (TOTAL NUMBER OF SPECIAL==
      ! ==          POINTS. THIS IS CHECKED BY THE SUBROUTINE.          ==
      ! == ISTRIZ . INDICATES WHETHER ADDITIONAL MESH POINTS SHOULD BE  ==
      ! ==          GENERATED BY APPLYING GROUP OPERATIONS TO THE MESH. ==
      ! ==          ISTRIZ=+1 MEANS SYMMETRIZE                          ==
      ! ==          ISTRIZ=-1 MEANS DO NOT SYMMETRIZE                   ==
      ! == THE FOLLOWING INPUT DATA MAY BE OBTAINED FROM THE SBRT.      ==
      ! == B1,B2,B3 .. RECIPROCAL LATTICE VECTORS, NOT MULTIPLIED BY    ==
      ! ==          GROUP1: ANY 2PI (IN UNITS RECIPROCAL TO THOSE       ==
      ! ==                  OF A1,A2,A3)                                ==
      ! == INV .... CODE INDICATING WHETHER WE WISH TO ADD THE INVERSION==
      ! ==          TO THE POINT GROUP OF THE CRYSTAL OR NOT (IN THE    ==
      ! ==          CASE THAT THE POINT GROUP DOES NOT CONTAIN ANY).    ==
      ! ==          INV=0 MEANS: DO NOT ADD INVERSION                   ==
      ! ==          INV.NE.0 MEANS: ADD THE INVERSION                   ==
      ! ==          INV.NE.0 SHOULD BE THE STANDARD CHOICE WHEN SPPT2   ==
      ! ==          IS USED IN RECIPROCAL SPACE - IN ORDER TO MAKE      ==
      ! ==          USE OF THE HERMITICITY OF HAMILTONIAN.              ==
      ! ==          WHEN USED IN DIRECT SPACE, THE RIGHT CHOICE OF INV  ==
      ! ==          WILL DEPEND ON THE NATURE OF THE PHYSICAL PROBLEM.  ==
      ! ==          IN THE CASES WHERE THE INVERSION IS ADDED BY THE    ==
      ! ==          SWITCH INV, THE LIST IB WILL NOT BE MODIFIED BUT IN ==
      ! ==          THE OUTPUT LIST LROT SOME OF THE OPERATIONS WILL    ==
      ! ==          APPEAR WITH NEGATIVE SIGN; THIS MEANS THAT THEY HAVE==
      ! ==          TO BE APPLIED MULTIPLIED BY INVERSION.              ==
      ! == NC ..... TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP OF THE  ==
      ! ==          CRYSTAL                                             ==
      ! == IB ..... LIST OF THE ROTATIONS CONSTITUTING THE POINT GROUP  ==
      ! ==          OF THE CRYSTAL. THE NUMBERING IS THAT DEFINED IN    ==
      ! ==          WORLTON AND WARREN, I.E. THE ONE MATERIALIZED IN THE==
      ! ==          ARRAY R (SEE BELOW)                                 ==
      ! ==          ONLY THE FIRST NC ELEMENTS OF THE ARRAY IB ARE      ==
      ! ==          MEANINGFUL                                          ==
      ! == R ...... LIST OF THE 3 X 3 ROTATION MATRICES                 ==
      ! ==          (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)    ==
      ! ==          ALL 48 OR 24 MATRICES ARE LISTED.                   ==
      ! == NCBRAV . TOTAL NUMBER OF ELEMENTS IN RBRAV                   ==
      ! == IBRAV .. LIST OF NCBRAV OPERATIONS OF THE BRAVAIS LATTICE    ==
      ! == DELTA    REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)           ==
      ! ==--------------------------------------------------------------==
      ! == OUTPUT DATA:                                                 ==
      ! == NTOT ... TOTAL NUMBER OF SPECIAL POINTS                      ==
      ! ==          IF NTOT APPEARS NEGATIVE, THIS IS AN ERROR SIGNAL   ==
      ! ==          WHICH MEANS THAT THE DIMENSION NKPOINT WAS CHOSEN   ==
      ! ==          TOO SMALL SO THAT THE ARRAYS WVKL ETC. CANNOT       ==
      ! ==          ACCOMODATE ALL THE GENERATED SPECIAL POINTS.        ==
      ! ==          IN THIS CASE THE ARRAYS WILL BE FILLED UP TO NKPOINT==
      ! ==          AND FURTHER GENERATION OF NEW POINTS WILL BE        ==
      ! ==          INTERRUPTED.                                        ==
      ! == WVKL ... LIST OF SPECIAL POINTS.                             ==
      ! ==          CARTESIAN COORDINATES AND NOT MULTIPLIED BY 2*PI.   ==
      ! ==          ONLY THE FIRST NTOT VECTORS ARE MEANINGFUL          ==
      ! ==          ALTHOUGH NO 2 POINTS FROM THE LIST ARE EQUIVALENT   ==
      ! ==          BY SYMMETRY, THIS SUBROUTINE STILL HAS A KIND OF    ==
      ! ==          'BEAUTY DEFECT': THE POINTS FINALLY                 ==
      ! ==          SELECTED ARE NOT NECESSARILY SITUATED IN A          ==
      ! ==          'COMPACT' IRREDUCIBLE BRILL.ZONE; THEY MIGHT LIE IN ==
      ! ==          DIFFERENT IRREDUCIBLE PARTS OF THE B.Z. - BUT THEY  ==
      ! ==          DO REPRESENT AN IRREDUCIBLE SET FOR INTEGRATION     ==
      ! ==          OVER THE ENTIRE B.Z.                                ==
      ! == LWGHT ... THE LIST OF WEIGHTS OF THE CORRESPONDING POINTS.   ==
      ! ==          THESE WEIGHTS ARE NOT NORMALIZED (JUST INTEGERS)    ==
      ! == LROT ... FOR EACH SPECIAL POINT THE 'UNFOLDING ROTATIONS'    ==
      ! ==          ARE LISTED. IF E.G. THE WEIGHT OF THE I-TH SPECIAL  ==
      ! ==          POINT IS LWGHT(I), THEN THE ROTATIONS WITH NUMBERS  ==
      ! ==          LROT(J,I), J=1,2,...,LWGHT(I) WILL 'SPREAD' THIS    ==
      ! ==          SINGLE POINT FROM THE IRREDUCIBLE PART OF B.Z. INTO ==
      ! ==          SEVERAL POINTS IN AN ELEMENTARY UNIT CELL           ==
      ! ==          (PARALLELOPIPED) OF THE RECIPROCAL SPACE.           ==
      ! ==          SOME OPERATION NUMBERS IN THE LIST LROT MAY APPEAR  ==
      ! ==          NEGATIVE, THIS MEANS THAT THE CORRESPONDING ROTATION==
      ! ==          HAS TO BE APPLIED WITH INVERSION (THE LATTER HAVING ==
      ! ==          BEEN ARTIFICIALLY ADDED AS SYMMETRY OPERATION IN    ==
      ! ==          CASE INV.NE.0).NO OTHER EFFORT WAS TAKEN,TO RENUMBER==
      ! ==          THE ROTATIONS WITH MINUS SIGN OR TO EXTEND THE      ==
      ! ==          LIST OF THE POINT-GROUP OPERATIONS IN THE LIST NB.  ==
      ! == INCLUD ... INTEGER ARRAY USED BY SPPT2 INCLUD(NKPOINT)       ==
      ! ==          THE FIRST BIT (0) IS USED BY THE ROUTINE.           ==
      ! ==          THE OTHER BITS GIVE THE K-POINT INDEX IN            ==
      ! ==          THE SPECIAL K-POINT TABLE.                          ==
      ! ==--------------------------------------------------------------==
      ! == NHASH    USED BY MESH ROUTINE                                ==
      ! == LIST     INTEGER ARRAY USED BY MESH  LIST(NHASH+NKPOINT)     ==
      ! == RLIST    real(8) :: ARRAY USED BY MESH  RLIST(3,NKPOINT)        ==
      ! ==--------------------------------------------------------------==
      ! == Use bit manipulations functions                              ==
      ! ==  IBSET(I,POS) sets the bit POS to 1 in I integer             ==
      ! ==  IBCLR(I,POS) clears the bit POS to 1 in I integer           ==
      ! ==  BTEST(I,POS) .TRUE. if bit POS is 1 in I integer            ==
      ! ==--------------------------------------------------------------==
      INTEGER                                            :: iout, iq1, iq2, iq3
      REAL(dp)                                           :: wvk0(3)
      INTEGER                                            :: nkpoint
      REAL(dp)                                           :: a1(3), a2(3), a3(3), b1(3), b2(3), b3(3)
      INTEGER                                            :: inv, nc, ib(48)
      REAL(dp)                                           :: r(3, 3, 48)
      INTEGER                                            :: ntot
      REAL(dp)                                           :: wvkl(3, nkpoint)
      INTEGER                                            :: lwght(nkpoint), lrot(48, nkpoint), &
                                                            ncbrav, ibrav(48), istriz, nhash, &
                                                            includ(nkpoint), list(nkpoint + nhash)
      REAL(dp)                                           :: rlist(3, nkpoint), delta

      INTEGER, PARAMETER                                 :: no = 0, nrsdir = 100

      INTEGER                                            :: i, i1, i2, i3, ibsign, igarb0, igarbage, &
                                                            igarbg, ii, imesh, iop, iplace, &
                                                            iremov, iwvk, j, jplace, k, n, nplane
      REAL(dp)                                           :: diff, proja(3), projb(3), &
                                                            rsdir(4, nrsdir), ur1, ur2, ur3, &
                                                            wva(3), wvk(3)

! ==--------------------------------------------------------------==

      ntot = 0
      DO i = 1, nkpoint
         lrot(1, i) = 1
         DO j = 2, 48
            lrot(j, i) = 0
         END DO
      END DO
      DO i = 1, nkpoint
         includ(i) = no
      END DO
      DO i = 1, 3
         wva(i) = 0._dp
      END DO
      ! ==--------------------------------------------------------------==
      ! == DEFINE THE 1ST BRILLOUIN ZONE                                ==
      ! ==--------------------------------------------------------------==
      CALL bzdefine(iout, b1, b2, b3, rsdir, nplane, delta)
      ! ==--------------------------------------------------------------==
      ! == Generation of the mesh (they are not multiplied by 2*pi) by  ==
      ! == the Monkhorst/Pack algorithm, supplemented by all rotations  ==
      ! ==--------------------------------------------------------------==
      ! Initialize the list of vectors
      iplace = -2
      CALL mesh(iout, wva, iplace, igarb0, igarbg, nkpoint, nhash, &
                list, rlist, delta)
      imesh = 0
      DO i1 = 1, iq1
         DO i2 = 1, iq2
            DO i3 = 1, iq3
               ur1 = REAL(1 + iq1 - 2*i1, kind=dp)/REAL(2*iq1, kind=dp)
               ur2 = REAL(1 + iq2 - 2*i2, kind=dp)/REAL(2*iq2, kind=dp)
               ur3 = REAL(1 + iq3 - 2*i3, kind=dp)/REAL(2*iq3, kind=dp)
               DO i = 1, 3
                  wvk(i) = ur1*b1(i) + ur2*b2(i) + ur3*b3(i) + wvk0(i)
               END DO
               ! Reduce WVK to the 1st Brillouin zone
               CALL bzrduc(wvk, a1, a2, a3, b1, b2, b3, rsdir, &
                           nrsdir, nplane, delta)
               IF (istriz .EQ. 1) THEN
                  ! Symmetrization of the k-points mesh.
                  ! Apply all the Bravais lattice operations to WVK
                  DO iop = 1, ncbrav
                     DO i = 1, 3
                        wva(i) = 0._dp
                        DO j = 1, 3
                           wva(i) = wva(i) + r(i, j, ibrav(iop))*wvk(j)
                        END DO
                     END DO
                     ! Check that WVA is inside the 1 Bz.
                     IF (.NOT. inside_bz(wva, rsdir, nplane, delta)) GOTO 450
                     ! Place WVA in list
                     iplace = 0
                     CALL mesh(iout, wva, iplace, igarb0, igarbg, &
                               nkpoint, nhash, list, rlist, delta)
                     ! If WVA was new (and therefore inserted),
                     ! IPLACE is the number.
                     IF (iplace .GT. 0) imesh = iplace
                     IF (iplace .GT. nkpoint) GOTO 470
                  END DO
               ELSE
                  ! Place WVK in list
                  iplace = 0
                  CALL mesh(iout, wvk, iplace, igarb0, igarbg, &
                            nkpoint, nhash, list, rlist, delta)
                  imesh = iplace
                  IF (iplace .GT. nkpoint) GOTO 470
               END IF
            END DO
         END DO
      END DO
!deb
!deb get full mesh
!deb
      IF (iout .GT. 0) THEN
         ! IMESH: Number of k points in the mesh.
         IF (iout > 0) &
            WRITE (iout, &
                   '(" KPSYM| THE WAVEVECTOR MESH CONTAINS ",I5," POINTS")') imesh
         IF (iout > 0) &
            WRITE (iout, '(" KPSYM| THE POINTS ARE:")')
         DO ii = 1, imesh
            i = ii
            CALL mesh(iout, wva, i, igarb0, igarbg, nkpoint, nhash, &
                      list, rlist, delta)
            IF (MOD(i, 2) .EQ. 1) THEN
               IF (iout > 0) &
                  WRITE (iout, '(1X,I5,3F10.4)', advance="no") i, wva
            ELSE
               IF (iout > 0) &
                  WRITE (iout, '(1X,I5,3F10.4)') i, wva
            END IF
         END DO
         IF (iout > 0) &
            WRITE (iout, *)
      END IF
      ! ==--------------------------------------------------------------==
      IF (istriz .EQ. 1) THEN
         ! Now figure out if any special point difference (K - K'') is an
         ! integral multiple of a reciprocal-space vector
         iremov = 0
         DO i = 1, (imesh - 1)
            iplace = i
            CALL mesh(iout, wva, iplace, igarb0, igarbg, &
                      nkpoint, nhash, list, rlist, delta)
            ! Project WVA onto B1,2,3:
            proja(1) = 0._dp
            proja(2) = 0._dp
            proja(3) = 0._dp
            DO k = 1, 3
               proja(1) = proja(1) + wva(k)*a1(k)
               proja(2) = proja(2) + wva(k)*a2(k)
               proja(3) = proja(3) + wva(k)*a3(k)
            END DO
            ! Now loop over all the rest of the mesh points
            DO j = (i + 1), imesh
               jplace = j
               CALL mesh(iout, wvk, jplace, igarb0, igarbg, &
                         nkpoint, nhash, list, rlist, delta)
               ! Project WVK onto B1,2,3:
               projb(1) = 0._dp
               projb(2) = 0._dp
               projb(3) = 0._dp
               DO k = 1, 3
                  projb(1) = projb(1) + wvk(k)*a1(k)
                  projb(2) = projb(2) + wvk(k)*a2(k)
                  projb(3) = projb(3) + wvk(k)*a3(k)
               END DO
               ! Check (PROJA - PROJB): Is it integral ?
               DO k = 1, 3
                  diff = proja(k) - projb(k)
                  IF (ABS(REAL(NINT(diff), kind=dp) - diff) .GT. delta) GOTO 280
               END DO
               ! DIFF is integral: remove WVK from mesh:
               CALL remove(wvk, jplace, igarb0, igarbg, &
                           nkpoint, nhash, list, rlist, delta)
               ! If WVK actually removed, increment IREMOV
               IF (jplace .GT. 0) iremov = iremov + 1
280            CONTINUE
            END DO
         END DO
         IF ((iremov .GT. 0 .AND. iout .GT. 0) .AND. iout > 0) &
            WRITE (iout, '(A,A,/,A,1X,I6,A,/)') &
            ' KPSYM| SOME OF THESE MESH POINTS ARE RELATED BY LATTICE ', &
            'TRANSLATION VECTORS', &
            ' KPSYM|', iremov, ' OF THE MESH POINTS REMOVED.'
      END IF
      ! ==--------------------------------------------------------------==
      ! == IN THE MESH OF WAVEVECTORS, NOW SEARCH FOR EQUIVALENT POINTS:==
      ! == THE INVERSION (TIME REVERSAL !) MAY BE USED.                 ==
      ! ==--------------------------------------------------------------==
      DO iwvk = 1, imesh
         ! IF(INCLUD(IWVK) .EQ. YES) GOTO 350
         IF (BTEST(includ(iwvk), 0)) GOTO 350
         ! IWVK has not been encountered previously: new special point,
         ! (only if WVK is not a garbage vector, however.)
         ! INCLUD(IWVK) = YES
         includ(iwvk) = IBSET(includ(iwvk), 0)
         iplace = iwvk
         CALL mesh(iout, wvk, iplace, igarb0, igarbg, &
                   nkpoint, nhash, list, rlist, delta)
         ! Find out whether Wvk is in the garbage list
         CALL garbag(wvk, igarbage, igarb0, &
                     nkpoint, nhash, list, rlist, delta)
         IF (igarbage .GT. 0) GOTO 350
         ntot = ntot + 1
         ! Give the index in the special k points table.
         includ(iwvk) = includ(iwvk) + ntot*2
         DO i = 1, 3
            wvkl(i, ntot) = wvk(i)
         END DO
         lwght(ntot) = 1
         ! ==-----------------------------------------------------------==
         ! Find all the equivalent points (symmetry given by atoms)
         DO n = 1, nc
            ! Rotate:
            DO i = 1, 3
               wva(i) = 0._dp
               DO j = 1, 3
                  wva(i) = wva(i) + r(i, j, ib(n))*wvk(j)
               END DO
            END DO
            ibsign = +1
363         CONTINUE
            ! Find WVA in the list
            iplace = -1
            CALL mesh(iout, wva, iplace, igarb0, igarbg, &
                      nkpoint, nhash, list, rlist, delta)
            IF (iplace .EQ. 0) THEN
               IF (istriz .EQ. -1) THEN
                  ! No symmetrisation -> WVA not in the list
                  GOTO 364
               ELSE
                  ! Find out whether WVA is in the garbage list
                  CALL garbag(wva, igarbage, igarb0, &
                              nkpoint, nhash, list, rlist, delta)
                  IF (igarbage .EQ. 0) THEN
                     ! I think this case is impossible (NC <= NCBRAV)
                     ! Error message
                     GOTO 490
                  END IF
               END IF
            END IF
            ! Find out whether WVA is in the garbage list
            CALL garbag(wva, igarbage, igarb0, &
                        nkpoint, nhash, list, rlist, delta)
            IF (igarbage .GT. 0) GOTO 370
            ! Was WVA encountered before ?
            ! IF(INCLUD(IPLACE) .EQ. YES) GOTO 364
            IF (BTEST(includ(iplace), 0)) GOTO 364
            ! Increment weight.
            lwght(ntot) = lwght(ntot) + 1
            lrot(lwght(ntot), ntot) = ib(n)*ibsign
            ! INCLUD(IPLACE) = YES
            includ(iplace) = IBSET(includ(iplace), 0)
            ! This k-point is an image of a special k-point.
            ! Put the index of the special k-point.
            includ(iplace) = includ(iplace) + ntot*2
364         CONTINUE
            IF (ibsign .EQ. -1 .OR. inv .EQ. 0) GOTO 370
            ! The case where we also apply the inversion to WVA
            ! Repeat the search, but for -WVA
            ibsign = -1
            DO i = 1, 3
               wva(i) = -wva(i)
            END DO
            GOTO 363
370         CONTINUE
         END DO
350      CONTINUE
      END DO
      ! ==--------------------------------------------------------------==
      ! == TOTAL NUMBER OF SPECIAL POINTS: NTOT                         ==
      ! == BEFORE USING THE LIST WVKL AS WAVE VECTORS, THEY HAVE TO BE  ==
      ! == MULTIPLIED BY 2*PI                                           ==
      ! == THE LIST OF WEIGHTS LWGHT IS NOT NORMALIZED                  ==
      ! ==--------------------------------------------------------------==
      IF (ntot .GT. nkpoint) THEN
         IF (iout > 0) &
            WRITE (iout, *) 'IN SPPT2 NUMBER OF SPECIAL POINTS = ', ntot
         IF (iout > 0) &
            WRITE (iout, *) 'BUT NKPOINT = ', nkpoint
         ntot = -1
      END IF
      IF (iout .GT. 0) THEN
         ! Write the index table relating k points in the mesh
         ! with special k points
         IF (iout > 0) &
            WRITE (iout, '(/,A,4X,A)') &
            ' KPSYM|', 'CROSS TABLE RELATING MESH POINTS WITH SPECIAL POINTS:'
         IF (iout > 0) &
            WRITE (iout, '(5(4X,"IK -> SK"))')
         DO i = 1, imesh
            iplace = includ(i)/2
            IF (iout > 0) &
               WRITE (iout, '(1X,I5,1X,I5)', advance="no") i, iplace
            IF ((MOD(i, 5) .EQ. 0) .AND. iout > 0) &
               WRITE (iout, *)
         END DO
         IF ((MOD(j - 1, 5) .NE. 0) .AND. iout > 0) &
            WRITE (iout, *)
      END IF
      RETURN
      ! ==--------------------------------------------------------------==
      ! == ERROR MESSAGES                                               ==
      ! ==--------------------------------------------------------------==
450   CONTINUE
      IF (iout > 0) &
         WRITE (iout, '(A,/)') ' SUBROUTINE SPPT2 *** FATAL ERROR ***'
      IF (iout > 0) &
         WRITE (iout, '(A,3F10.4,/,A,3F10.4,A,/,A,I3,A)') &
         ' THE VECTOR     ', wva, &
         ' GENERATED FROM ', wvk, ' IN THE BASIC MESH', &
         ' BY ROTATION NO. ', ibrav(iop), ' IS OUTSIDE THE 1BZ'
      CALL stopgm('SPPT2', 'VECTOR OUTSIDE THE 1BZ')
      ! ==--------------------------------------------------------------==
470   CONTINUE
      IF (iout > 0) &
         WRITE (iout, '(A,/)') ' SUBROUTINE SPPT2 *** FATAL ERROR ***'
      IF (iout > 0) &
         WRITE (iout, *) 'MESH SIZE EXCEEDS NKPOINT=', nkpoint
      CALL stopgm('SPPT2', 'MESH SIZE EXCEEDED')
      ! ==--------------------------------------------------------------==
490   CONTINUE
      IF (iout > 0) &
         WRITE (iout, '(A,/)') ' SUBROUTINE SPPT2 *** FATAL ERROR ***'
      IF (iout > 0) &
         WRITE (iout, '(A,3F10.4,/,A,3F10.4,A,/,A,I3,A)') &
         ' THE VECTOR     ', wva, &
         ' GENERATED FROM ', wvk, ' IN THE BASIC MESH', &
         ' BY ROTATION NO. ', ib(n), ' IS NOT IN THE LIST'
      CALL stopgm('SPPT2', 'VECTOR NOT IN THE LIST')
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE sppt2
! **************************************************************************************************
!> \brief ...
!> \param iout ...
!> \param wvk ...
!> \param iplace ...
!> \param igarb0 ...
!> \param igarbg ...
!> \param nmesh ...
!> \param nhash ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE mesh(iout, wvk, iplace, igarb0, igarbg, &
                   nmesh, nhash, list, rlist, delta)
      ! ==--------------------------------------------------------------==
      ! == MESH MAINTAINS A LIST OF VECTORS FOR PLACEMENT AND/OR LOOKUP ==
      ! ==                                                              ==
      ! == ADDITIONAL ENTRY POINTS: REMOVE .... REMOVE VECTOR FROM LIST ==
      ! ==                          GARBAG .... WAS VECTOR REMOVED ?    ==
      ! ==                                                              ==
      ! == WVK ....... VECTOR                                           ==
      ! == IPLACE .... ON INPUT: -2 MEANS: INITIALIZE  THE LIST         ==
      ! ==                                 (AND RETURN)                 ==
      ! ==                       -1 MEANS: FIND WVK IN THE LIST         ==
      ! ==                        0 MEANS: ADD  WVK TO THE LIST         ==
      ! ==                       >0 MEANS: RETURN WVK NO. IPLACE        ==
      ! ==            ON OUTPUT: THE POSITION ASSIGNED TO WVK           ==
      ! ==                       (=0 IF WVK IS NOT IN THE LIST)         ==
      ! == DELTA      REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)         ==
      ! ==--------------------------------------------------------------==
      INTEGER                                            :: iout
      REAL(dp)                                           :: wvk(3)
      INTEGER                                            :: iplace, igarb0, igarbg, nmesh, nhash, &
                                                            list(nhash + nmesh)
      REAL(dp)                                           :: rlist(3, nmesh), delta

      INTEGER, PARAMETER                                 :: nil = 0

      INTEGER                                            :: i, ihash, ipoint, j
      INTEGER, SAVE                                      :: istore
      REAL(dp)                                           :: delta1, rhash

! ==--------------------------------------------------------------==
! == Initialization                                               ==
! ==--------------------------------------------------------------==

      delta1 = 10._dp*delta
      IF (iplace .LE. -2) THEN
         DO i = 1, nhash + nmesh
            list(i) = nil
         END DO
         istore = 1
         ! IGARB0 points to a linked list of removed WVKS (the garbage).
         igarb0 = 0
         igarbg = 0
         RETURN
         ! ==--------------------------------------------------------------==
      ELSEIF ((iplace .GT. -2) .AND. (iplace .LE. 0)) THEN
         ! The particular HASH function used in this case:
         rhash = 0.7890_dp*wvk(1) &
                 + 0.6810_dp*wvk(2) &
                 + 0.5811_dp*wvk(3) + delta
         ihash = INT(ABS(rhash)*REAL(nhash, kind=dp))
         ihash = MOD(ihash, nhash) + nmesh + 1
         ! Search for WVK in linked list
         ipoint = list(ihash)
         DO i = 1, 100
            ! List exhausted
            IF (ipoint .EQ. nil) GOTO 130
            ! Compare WVK with this element
            DO j = 1, 3
               IF (ABS(wvk(j) - rlist(j, ipoint)) .GT. delta1) GOTO 115
            END DO
            ! WVK located
            GOTO 160
            ! Next element of list
115         CONTINUE
            ihash = ipoint
            ipoint = list(ihash)
         END DO
         ! List too long
         IF (iout > 0) &
            WRITE (iout, '(2A,/,A)') &
            ' SUBROUTINE MESH *** FATAL ERROR *** LINKED LIST', &
            ' TOO LONG ***', ' CHOOSE A BETTER HASH-FUNCTION'
         CALL stopgm('MESH', 'WARNING')
         ! WVK was not found
130      CONTINUE
         IF (iplace .EQ. -1) THEN
            ! IPLACE=-1 : search for WVK unsuccessful
            iplace = 0
            RETURN
         ELSE
            ! IPLACE=0: add WVK to the list
            list(ihash) = istore
            IF (istore .GT. nmesh) THEN
               IF (iout > 0) &
                  WRITE (iout, '(A)') 'SUBROUTINE MESH *** FATAL ERROR ***'
               IF (iout > 0) &
                  WRITE (iout, '(A,I10,A,/,A,3F10.5)') &
                  ' ISTORE=', istore, ' EXCEEDS DIMENSIONS', &
                  ' WVK = ', wvk
               CALL stopgm('MESH', 'WARNING')
            END IF
            list(istore) = nil
            DO i = 1, 3
               rlist(i, istore) = wvk(i)
            END DO
            istore = istore + 1
            iplace = istore - 1
            RETURN
         END IF
         ! WVK was found
160      CONTINUE
         IF (iplace .EQ. 0) RETURN
         ! IPLACE=-1
         iplace = ipoint
         RETURN
      ELSE
         ! ==--------------------------------------------------------------==
         ! == Return a wavevector (IPLACE > 0)                             ==
         ! ==--------------------------------------------------------------==
         ipoint = iplace
         IF (ipoint .GE. istore) GOTO 190
         DO i = 1, 3
            wvk(i) = rlist(i, ipoint)
         END DO
         RETURN
      END IF
      ! ==--------------------------------------------------------------==
      ! == Error - beyond list                                          ==
      ! ==--------------------------------------------------------------==
190   CONTINUE
      IF (iout > 0) &
         WRITE (iout, '(A,/,A,I5,A,/)') &
         ' SUBROUTINE MESH *** WARNING ***', &
         ' IPLACE = ', iplace, &
         ' IS BEYOND THE LISTS - WVK SET TO 1.0E38'
      DO i = 1, 3
         wvk(i) = 1.0e38_dp
      END DO
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE mesh
! **************************************************************************************************
!> \brief ...
!> \param wvk ...
!> \param iplace ...
!> \param igarb0 ...
!> \param igarbg ...
!> \param nmesh ...
!> \param nhash ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE remove(wvk, iplace, igarb0, igarbg, &
                     nmesh, nhash, list, rlist, delta)
      ! ==--------------------------------------------------------------==
      ! == ENTRY POINT FOR REMOVING A WAVEVECTOR                        ==
      ! ==                                                              ==
      ! == INPUT:                                                       ==
      ! ==   WVK(3)                                                     ==
      ! ==   DELTA   REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)          ==
      ! == OUTPUT:
      ! ==   IPLACE .....1 IF WVK WAS REMOVED                          ==
      ! ==                0 IF WVK WAS NOT REMOVED                      ==
      ! ==                  (WVK NOT IN THE LINKED LISTS)               ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: wvk(3)
      INTEGER                                            :: iplace, igarb0, igarbg, nmesh, nhash, &
                                                            list(nhash + nmesh)
      REAL(dp)                                           :: rlist(3, nmesh), delta

      INTEGER, PARAMETER                                 :: nil = 0

      INTEGER                                            :: i, ihash, ipoint, j
      REAL(dp)                                           :: delta1, rhash

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

      delta1 = 10._dp*delta
      ! The particular hash function used in this case:
      rhash = 0.7890_dp*wvk(1) &
              + 0.6810_dp*wvk(2) &
              + 0.5811_dp*wvk(3) + delta
      ihash = INT(ABS(rhash)*REAL(nhash, kind=dp))
      ihash = MOD(ihash, nhash) + nmesh + 1
      ! Search for WVK in linked list
      ipoint = list(ihash)
      DO i = 1, 100
         ! List exhausted
         IF (ipoint .EQ. nil) THEN
            ! WVK was not found in the mesh:
            iplace = 0
            RETURN
         END IF
         ! Compare WVK with this element
         DO j = 1, 3
            IF (ABS(wvk(j) - rlist(j, ipoint)) .GT. delta1) GOTO 215
         END DO
         ! WVK located, now remove it from the list:
         list(ihash) = list(ipoint)
         ! LIST(IHASH) now points to the next element in the list,
         ! and the present WVK has become garbage.
         ! Add WVK to the list of garbage:
         IF (igarb0 .EQ. 0) THEN
            ! Start up the garbage list:
            igarb0 = ipoint
         ELSE
            list(igarbg) = ipoint
         END IF
         igarbg = ipoint
         list(igarbg) = nil
         iplace = 1
         RETURN
         ! Next element of list
215      CONTINUE
         ihash = ipoint
         ipoint = list(ihash)
      END DO
      ! List too long
      CALL stopgm('MESH', 'LIST TOO LONG')
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE remove
! **************************************************************************************************
!> \brief ...
!> \param wvk ...
!> \param iplace ...
!> \param igarb0 ...
!> \param nmesh ...
!> \param nhash ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE garbag(wvk, iplace, igarb0, &
                     nmesh, nhash, list, rlist, delta)
      ! ==--------------------------------------------------------------==
      ! == ENTRY POINT FOR CHECKING IF A WAVEVECTOR                     ==
      ! ==             IS IN THE GARBAGE LIST                           ==
      ! == INPUT:                                                       ==
      ! ==   WVK(3)                                                     ==
      ! ==   DELTA      REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)       ==
      ! ==                                                              ==
      ! == OUTPUT:                                                      ==
      ! ==   IPLACE  ..... I > 0 IS THE PLACE IN THE GARBAGE LIST       ==
      ! ==                     0 IF WVK NOT AMONG THE GARBAGE           ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: wvk(3)
      INTEGER                                            :: iplace, igarb0, nmesh, nhash, &
                                                            list(nhash + nmesh)
      REAL(dp)                                           :: rlist(3, nmesh), delta

      INTEGER, PARAMETER                                 :: nil = 0

      INTEGER                                            :: i, ihash, ipoint, j
      REAL(dp)                                           :: delta1

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

      delta1 = 10._dp*delta
      ! Search for WVK in linked list
      ! Point to the garbage list
      ipoint = igarb0
      DO i = 1, nmesh
         ! LIST EXHAUSTED
         IF (ipoint .EQ. nil) THEN
            ! WVK was not found in the mesh:
            iplace = 0
            RETURN
         END IF
         ! Compare WVK with this element
         DO j = 1, 3
            IF (ABS(wvk(j) - rlist(j, ipoint)) .GT. delta1) GOTO 315
         END DO
         ! WVK was located in the garbage list
         iplace = i
         RETURN
         ! Next element of list
315      CONTINUE
         ihash = ipoint
         ipoint = list(ihash)
      END DO
      ! List too long
      CALL stopgm('GARBAG', 'LIST TOO LONG')
      ! ==--------------------------------------------------------------==
      RETURN
   END SUBROUTINE garbag

! **************************************************************************************************
!> \brief ...
!> \param wvk ...
!> \param a1 ...
!> \param a2 ...
!> \param a3 ...
!> \param b1 ...
!> \param b2 ...
!> \param b3 ...
!> \param rsdir ...
!> \param nrsdir ...
!> \param nplane ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE bzrduc(wvk, a1, a2, a3, b1, b2, b3, rsdir, nrsdir, nplane, delta)
      ! ==--------------------------------------------------------------==
      ! == REDUCE WVK TO LIE ENTIRELY WITHIN THE 1ST BRILLOUIN ZONE     ==
      ! == BY ADDING B-VECTORS                                          ==
      ! == DELTA      REQUIRED ACCURACY (1.e-6_dp IS A GOOD VALUE)         ==
      ! ==--------------------------------------------------------------==
      REAL(dp)                                           :: wvk(3), a1(3), a2(3), a3(3), b1(3), &
                                                            b2(3), b3(3)
      INTEGER                                            :: nrsdir
      REAL(dp)                                           :: rsdir(4, nrsdir)
      INTEGER                                            :: nplane
      REAL(dp)                                           :: delta

      INTEGER, PARAMETER                                 :: nzones = 4, nnn = 2*nzones + 1, &
                                                            nn = nzones + 1

      INTEGER                                            :: i, i1, i2, i3, n1, n2, n3, nn1, nn2, nn3
      LOGICAL                                            :: inside
      REAL(dp)                                           :: wb(3), wva(3)

! ==--------------------------------------------------------------==
! Variables
! Look around +/- "NZONES" to locate vector
! NZONES may need to be increased for very anisotropic zones
! ==--------------------------------------------------------------==

      IF (.NOT. inside_bz(wvk, rsdir, nplane, delta)) THEN
         inside = .FALSE.
         ! Express WVK in the basis of B1,2,3.
         ! This permits an estimate of how far WVK is from the 1Bz.
         wb(1) = wvk(1)*a1(1) + wvk(2)*a1(2) + wvk(3)*a1(3)
         wb(2) = wvk(1)*a2(1) + wvk(2)*a2(2) + wvk(3)*a2(3)
         wb(3) = wvk(1)*a3(1) + wvk(2)*a3(2) + wvk(3)*a3(3)
         nn1 = NINT(wb(1))
         nn2 = NINT(wb(2))
         nn3 = NINT(wb(3))
         ! Look around the estimated vector for the one truly inside the 1Bz
         n1_loop: DO n1 = 1, nnn
            i1 = nn - n1 - nn1
            DO n2 = 1, nnn
               i2 = nn - n2 - nn2
               DO n3 = 1, nnn
                  i3 = nn - n3 - nn3
                  DO i = 1, 3
                     wva(i) = wvk(i) + REAL(i1, kind=dp)*b1(i) + REAL(i2, kind=dp)*b2(i) + &
                              REAL(i3, kind=dp)*b3(i)
                  END DO
                  inside = inside_bz(wva, rsdir, nplane, delta)
                  IF (inside) EXIT n1_loop
               END DO
            END DO
         END DO n1_loop
         CPASSERT(inside)
         wvk(1:3) = wva(1:3)
      END IF

   END SUBROUTINE bzrduc

! **************************************************************************************************
!> \brief Is wvk in the 1st Brillouin zone ?
!>        Check whether wvk lies inside all the planes that define the 1bz.
!> \param wvk ...
!> \param rsdir ...
!> \param nplane ...
!> \param delta ...
!> \return ...
! **************************************************************************************************
   FUNCTION inside_bz(wvk, rsdir, nplane, delta) RESULT(inbz)
      REAL(KIND=dp), DIMENSION(3)                        :: wvk
      REAL(KIND=dp), DIMENSION(:, :)                     :: rsdir
      INTEGER                                            :: nplane
      REAL(KIND=dp)                                      :: delta
      LOGICAL                                            :: inbz

      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: projct

      inbz = .TRUE.
      DO n = 1, nplane
         projct = (rsdir(1, n)*wvk(1) + rsdir(2, n)*wvk(2) + rsdir(3, n)*wvk(3))/rsdir(4, n)
         IF (ABS(projct) > 0.5_dp + delta) THEN
            inbz = .FALSE.
            EXIT
         END IF
      END DO

   END FUNCTION inside_bz

! **************************************************************************************************
!> \brief Find the vectors whose halves define the 1st Brillouin zone
!>        Output:
!>        nplane -- How many elements of rsdir contain normal vectors defining the planes
!>        Method:
!>        Starting with the parallelopiped spanned by b1,2,3 around the origin,
!>        vectors inside a sufficiently large sphere are tested to see whether
!>        the planes at 1/2*b will further confine the 1bz.
!>        The resulting vectors are not cleaned to avoid redundant planes
!> \param iout ...
!> \param b1 ...
!> \param b2 ...
!> \param b3 ...
!> \param rsdir ...
!> \param nplane ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE bzdefine(iout, b1, b2, b3, rsdir, nplane, delta)
      INTEGER                                            :: iout
      REAL(KIND=dp), DIMENSION(3)                        :: b1, b2, b3
      REAL(KIND=dp), DIMENSION(:, :)                     :: rsdir
      INTEGER                                            :: nplane
      REAL(KIND=dp)                                      :: delta

      INTEGER                                            :: i, i1, i2, i3, n, n1, n2, n3, nb1, nb2, &
                                                            nb3, nnb1, nnb2, nnb3, nrsdir
      REAL(KIND=dp)                                      :: b1len, b2len, b3len, bmax, projct
      REAL(KIND=dp), DIMENSION(3)                        :: bvec

      nrsdir = SIZE(rsdir, 2)

      b1len = b1(1)**2 + b1(2)**2 + b1(3)**2
      b2len = b2(1)**2 + b2(2)**2 + b2(3)**2
      b3len = b3(1)**2 + b3(2)**2 + b3(3)**2
      ! Lattice containing entirely the Brillouin zone
      bmax = b1len + b2len + b3len
      nb1 = INT(SQRT(bmax/b1len) + delta) + 1
      nb2 = INT(SQRT(bmax/b2len) + delta) + 1
      nb3 = INT(SQRT(bmax/b3len) + delta) + 1
      rsdir(:, :) = 0._dp
      ! 1Bz is certainly confined inside the 1/2(B1,B2,B3) parallelopiped
      rsdir(1:3, 1) = b1(1:3)
      rsdir(1:3, 2) = b2(1:3)
      rsdir(1:3, 3) = b3(1:3)
      rsdir(4, 1) = b1len
      rsdir(4, 2) = b2len
      rsdir(4, 3) = b3len
      ! Starting confinement: 3 planes
      nplane = 3
      nnb1 = 2*nb1 + 1
      nnb2 = 2*nb2 + 1
      nnb3 = 2*nb3 + 1

      DO n1 = 1, nnb1
         i1 = nb1 + 1 - n1
         DO n2 = 1, nnb2
            i2 = nb2 + 1 - n2
            DO n3 = 1, nnb3
               i3 = nb3 + 1 - n3
               IF (i1 .EQ. 0 .AND. i2 .EQ. 0 .AND. i3 .EQ. 0) GOTO 150
               DO i = 1, 3
                  bvec(i) = REAL(i1, kind=dp)*b1(i) + REAL(i2, kind=dp)*b2(i) + &
                            REAL(i3, kind=dp)*b3(i)
               END DO
               ! Does the plane of 1/2*BVEC narrow down the 1Bz ?
               DO n = 1, nplane
                  projct = 0.5_dp*(rsdir(1, n)*bvec(1) + rsdir(2, n)*bvec(2) &
                                   + rsdir(3, n)*bvec(3))/rsdir(4, n)
                  ! 1/2*BVEC is outside the Bz - skip this direction
                  ! The 1.e-6_dp takes care of single points touching the Bz,
                  ! and of the -(plane)
                  IF (ABS(projct) .GT. 0.5_dp - delta) GOTO 150
               END DO
               ! 1/2*BVEC further confines the 1Bz - include into RSDIR
               nplane = nplane + 1
               CPASSERT(nplane <= nrsdir)
               DO i = 1, 3
                  rsdir(i, nplane) = bvec(i)
               END DO
               ! Length squared
               rsdir(4, nplane) = bvec(1)**2 + bvec(2)**2 + bvec(3)**2
150            CONTINUE
            END DO
         END DO
      END DO

      IF (iout > 0) THEN
         WRITE (iout, '(A,I3,A,/,A,/,100(" KPSYM|",1X,3F10.4,/))') &
            ' KPSYM| The 1st Brillouin zone is confined by (at most)', &
            nplane, ' planes', &
            ' KPSYM| as defined by the +/- halves of the vectors:', &
            ((rsdir(i, n), i=1, 3), n=1, nplane)
      END IF

   END SUBROUTINE bzdefine

! **************************************************************************************************
!> \brief ...
!> \param a ...
!> \param b ...
! **************************************************************************************************
   SUBROUTINE stopgm(a, b)
      CHARACTER(LEN=*)                                   :: a, b

      CALL cp_warn(a, b)
      CPABORT("stopgm@kpsym")

   END SUBROUTINE stopgm

! **************************************************************************************************

END MODULE kpsym
