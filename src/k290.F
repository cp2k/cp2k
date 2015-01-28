!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief K-points and crystal symmetry routines
!> \author jgh
! *****************************************************************************
MODULE k290
  
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: invmat
  USE string_utilities,                ONLY: xstring
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: csym_type, pgrp, pgrd
  PUBLIC :: kp_sym_gen, release_csym_type, print_crys_symmetry

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'k290'

! *****************************************************************************
!> \brief CSM type
!> \par   Content:
!>       PLEVEL . Output level
!>       PUNIT .. Output unit
!>       NTVEC .. Number of primitive cells?
!>       DELTA .. Accuracy in determining symmetry elements
!>       IHG .... Point group of the primitive lattice, holohedral
!>       IHC .... Code distinguishing between hexagonal and cubic groups
!>       ISY .... Code indicating whether the space group is symmorphic
!>       LI ..... Inversions symmetry
!>       NC ..... Total number of elements in the point group
!>       INDPG .. Point group index
!>       IB ..... List of the rotations constituting the point group
!>       V ...... Nonprimitive translations (for nonsymmorphic groups)
!>       XB ..... Cartesian standard origin or coordinates
!>       ORIGIN . Crystal standard origin or coordinates
!>       R ...... List of the 3 x 3 rotation matrices
!>
! *****************************************************************************
  TYPE csym_type
      INTEGER                                     :: plevel
      INTEGER                                     :: istriz = -1
      INTEGER                                     :: punit
      INTEGER                                     :: nat, nsp, iq1, iq2, iq3, nkpoint, ntvec
      REAL(KIND=dp), DIMENSION (:,:), POINTER     :: xkappa
      INTEGER, DIMENSION (:), POINTER             :: ty
      REAL(KIND=dp)                               :: a1(3), a2(3), a3(3), alat
      REAL(KIND=dp)                               :: strain(3,3) = 0.0_dp
      REAL(KIND=dp)                               :: wvk0(3) = 0.0_dp
      REAL(KIND=dp)                               :: delta, origin(3), xb(3)
      INTEGER                                     :: ihg, ihc, isy, li, nc, indpg
      INTEGER, DIMENSION(48)                      :: ib
      REAL(KIND=dp), DIMENSION(3,48)              :: v
      REAL(KIND=dp), DIMENSION(3,3,48)            :: r
  END TYPE csym_type

!------------------------------------------------------------------------------!
!     Constants
!     Name of 48 rotations (convention Warren-Worlton)
!------------------------------------------------------------------------------!
  CHARACTER(len=10) :: rname_cubic(48) = (/ ' 1        ', ' 2[ 1 0 0]', &
      ' 2[ 0 1 0]', ' 2[ 0 0 1]', ' 3[-1-1-1]', ' 3[ 1 1-1]', ' 3[-1 1 1]', &
      ' 3[ 1-1 1]', ' 3[ 1 1 1]', ' 3[-1 1-1]', ' 3[-1-1 1]', ' 3[ 1-1-1]', &
      ' 2[-1 1 0]', ' 4[ 0 0 1]', ' 4[ 0 0-1]', ' 2[ 1 1 0]', ' 2[ 0-1 1]', &
      ' 2[ 0 1 1]', ' 4[ 1 0 0]', ' 4[-1 0 0]', ' 2[-1 0 1]', ' 4[ 0-1 0]', &
      ' 2[ 1 0 1]', ' 4[ 0 1 0]', '-1        ', '-2[ 1 0 0]', '-2[ 0 1 0]', &
      '-2[ 0 0 1]', '-3[-1-1-1]', '-3[ 1 1-1]', '-3[-1 1 1]', '-3[ 1-1 1]', &
      '-3[ 1 1 1]', '-3[-1 1-1]', '-3[-1-1 1]', '-3[ 1-1-1]', '-2[-1 1 0]', &
      '-4[ 0 0 1]', '-4[ 0 0-1]', '-2[ 1 1 0]', '-2[ 0-1 1]', '-2[ 0 1 1]', &
      '-4[ 1 0 0]', '-4[-1 0 0]', '-2[-1 0 1]', '-4[ 0-1 0]', '-2[ 1 0 1]', &
      '-4[ 0 1 0]'/)
  CHARACTER(len=11) :: rname_hexa(24) = (/ ' 1         ', ' 6[ 0 0. 1]', &
      ' 3[ 0 0. 1]', ' 2[ 0 0. 1]', ' 3[ 0 0.-1]', ' 6[ 0 0.-1]', &
      ' 2[ 0 1. 0]', ' 2[-1 1. 0]', ' 2[ 1 0. 0]', ' 2[ 2 1. 0]', &
      ' 2[ 1 1. 0]', ' 2[ 1 2. 0]', '-1         ', '-6[ 0 0. 1]', &
      '-3[ 0 0. 1]', '-2[ 0 0. 1]', '-3[ 0 0.-1]', '-6[ 0 0.-1]', &
      '-2[ 0 1. 0]', '-2[-1 1. 0]', '-2[ 1 0. 0]', '-2[ 2 1. 0]', &
      '-2[ 1 1. 0]', '-2[ 1 2. 0]'/)
  CHARACTER(len=12) :: icst(7) = (/ '   TRICLINIC', '  MONOCLINIC', &
      'ORTHORHOMBIC', '  TETRAGONAL', '       CUBIC', '    TRIGONAL', &
      '   HEXAGONAL'/)

  CHARACTER(len=5)                         :: pgrp(32)
  DATA pgrp/'    1', '  <1>', '    2', '    m', '  2/m', '    3', &
   '  <3>', '   32', '   3m', ' <3>m', '    4', '  <4>', '  4/m', &
   '  422', '  4mm', '<4>2m', '4/mmm', '    6', '  <6>', '  6/m', &
   '  622', '  6mm', '<6>m2', '6/mmm', '  222', '  mm2', '  mmm', &
   '   23', ' m<3>', '  432', '<4>3m', 'm<3>m'/

  CHARACTER(len=3)                         :: pgrd(32)
  DATA pgrd/'C1 ', 'Ci ', 'C2 ', 'Cs ', 'C2h', 'C3 ', 'S6 ', 'D3 ', &
     'C3v', 'D3d', 'C4 ', 'S4 ', 'C4h', 'D4 ', 'C4v', 'D2d', 'D4h', &
     'C6 ', 'C3h', 'C6h', 'D6 ', 'C6v', 'D3h', 'D6h', 'D2 ', 'C2v', &
     'D2h', 'T  ', 'Th ', 'O  ', 'Td ', 'Oh '/

!------------------------------------------------------------------------------!
!     INDPG group   indpg   group    indpg  group     indpg   group
!        1    1 (C1)    9    3m (C3v)   17  4/mmm(D4h)  25   222(D2)
!        2   <1>(Ci)   10   <3>m(D3d)   18     6 (C6)   26   mm2(C2v)
!        3    2 (C2)   11     4 (C4)    19    <6>(C3h)  27   mmm(D2h)
!        4    m (Cs)   12    <4>(S4)    20    6/m(C6h)  28   23 (T)
!        5   2/m(C2h)  13    4/m(C4h)   21    622(D6)   29  m<3>(Th)
!        6    3 (C3)   14    422(D4)    22    6mm(C6v)  30   432(O)
!        7   <3>(S6)   15    4mm(C4v)   23  <6>m2(D3h)  31 <4>3m(Td)
!        8   32 (D3)   16  <4>2m(D2d)   24  6/mmm(D6h)  32 m<3>m(Oh)
!------------------------------------------------------------------------------!
CONTAINS

! *****************************************************************************
!> \brief Release the CSYM type
!> \param csym  The CSYM type
!> \param error CP2K error handling
! *****************************************************************************
  SUBROUTINE release_csym_type(csym,error)
    TYPE(csym_type), POINTER                 :: csym
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'release_csym_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(csym),cp_failure_level,routineP,error,failure)

    IF ( ASSOCIATED(csym%xkappa) ) THEN
       DEALLOCATE(csym%xkappa,STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF
    IF ( ASSOCIATED(csym%ty) ) THEN
       DEALLOCATE(csym%ty,STAT=ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(csym,STAT=ierr)
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE release_csym_type

! *****************************************************************************
!> \brief ...
!> \param csym ...
!> \param coor ...
!> \param types ...
!> \param hmat ...
!> \param nk ...
!> \param shift ...
!> \param stress ...
!> \param symm ...
!> \param delta ...
!> \param wkpoint ...
!> \param xkpoint ...
!> \param iounit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE kp_sym_gen(csym,coor,types,hmat,nk,shift,stress,symm,delta,&
                        wkpoint,xkpoint,iounit,error)
    TYPE(csym_type), POINTER                 :: csym
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: coor
    INTEGER, DIMENSION(:), INTENT(IN)        :: types
    REAL(KIND=dp), INTENT(IN)                :: hmat(3,3)
    INTEGER, INTENT(IN), OPTIONAL            :: nk(3)
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: shift(3), stress(3,3)
    LOGICAL, INTENT(IN), OPTIONAL            :: symm
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: delta
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: wkpoint
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: xkpoint
    INTEGER, INTENT(IN), OPTIONAL            :: iounit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kp_sym_gen', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ik, isos, isum, n, &
                                                nat, nhash, nkk, nkpts, nsp
    INTEGER, ALLOCATABLE                     :: f0(:,:), includ(:), isc(:), &
                                                list(:), lrot(:,:), lwght(:)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: wsum
    REAL(KIND=dp), ALLOCATABLE               :: rlist(:,:), rx(:,:), &
                                                tvec(:,:), wvkl(:,:)

    failure = .FALSE.
    CALL timeset(routineN,handle)

    CPPrecondition(.NOT.ASSOCIATED(csym),cp_failure_level,routineP,error,failure)
    ALLOCATE(csym,stat=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)

    !..total number of atoms
    nat = SIZE(coor,2)
    !..allocate arrays for coordinates and atom types
    ALLOCATE (csym%xkappa(3,nat),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (csym%ty(nat),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    !..count number of atom types
    csym%ty = types
    nsp = 0
    DO
       n = MAXVAL(csym%ty)
       IF (n/=-100) THEN
          nsp = nsp + 1
          WHERE (csym%ty==n) csym%ty = -100
       ELSE
          EXIT
       END IF
    END DO

    !..copy coordinates and atom types
    csym%xkappa = coor
    csym%ty = types
    csym%nat = nat
    csym%nsp = nsp

    !..set values
    csym%a1(:) = hmat(:,1)
    csym%a2(:) = hmat(:,2)
    csym%a3(:) = hmat(:,3)
    IF ( ABS(csym%a1(1)) > 1.e-10_dp ) THEN
       csym%alat = csym%a1(1)
    ELSE
       csym%alat = SQRT(SUM(csym%a1(:)**2))
    END IF
    IF (PRESENT(nk)) THEN
       csym%iq1 = nk(1)
       csym%iq2 = nk(2)
       csym%iq3 = nk(3)
    ELSE
       csym%iq1 = 0
       csym%iq2 = 0
       csym%iq3 = 0
    END IF
    IF (PRESENT(shift)) csym%wvk0 = shift
    IF (PRESENT(stress)) csym%strain = stress
    IF (PRESENT(symm)) THEN
       csym%istriz = -1
       IF (symm) csym%istriz = 1
    END IF
    csym%nkpoint = csym%iq1*csym%iq2*csym%iq3
    nkk = csym%nkpoint
    !..allocate intermediate arrays
    ALLOCATE (rx(3,nat),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (tvec(3,nat),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (isc(nat),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (f0(49,nat),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (wvkl(3,nkk),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (lwght(nkk),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (lrot(48,nkk),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (includ(nkk),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    nhash = MAX(2000,nkk/10)
    ALLOCATE (list(nkk+nhash),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (rlist(3,nkk),STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)

    IF ( PRESENT(iounit) ) THEN
       csym%punit = iounit
    ELSE
       csym%punit = -1
    END IF
    IF ( PRESENT(delta) ) THEN
       csym%delta = delta
    ELSE
       csym%delta = 1.e-6_dp
    END IF
    !..calculate symmetry operations and generate kpoints
    CALL k290prg(csym%punit,csym%nat,csym%nkpoint,csym%nsp,csym%iq1,csym%iq2,csym%iq3,&
         csym%istriz,csym%a1,csym%a2,csym%a3,csym%alat,csym%strain,csym%xkappa,&
         rx,tvec,csym%ty,isc,f0,csym%ntvec,csym%wvk0,csym%ihg,csym%ihc,csym%isy,&
         csym%li,csym%nc,csym%indpg,csym%ib,csym%v,csym%r,csym%origin,csym%xb,&
         wvkl,lwght,lrot,nhash,includ,list,rlist,csym%delta,error)
    IF ( PRESENT(xkpoint) ) THEN
       !..return kpoints and weights
       CPPrecondition(PRESENT(wkpoint),cp_failure_level,routineP,error,failure)
       isum=0
       nkpts=0
       DO ik=1,nkk
          IF(lwght(ik) == 0) EXIT
          nkpts=ik
          isum = isum + lwght(nkpts)
       END DO
       wsum=REAL(isum,KIND=dp)
       ALLOCATE (xkpoint(3,nkpts),wkpoint(nkpts),STAT=isos)
       CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
       DO ik=1,nkpts
          wkpoint(ik) = REAL(lwght(ik),KIND=dp)/wsum
          xkpoint(1:3,ik) = wvkl(1:3,ik)
       END DO
    END IF
    !..deallocate intermediate arrays
    DEALLOCATE (rx,tvec,isc,f0,wvkl,lwght,lrot,includ,list,rlist,STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE kp_sym_gen

! *****************************************************************************
!> \brief ...
!> \param csym ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE print_crys_symmetry(csym,error)
    TYPE(csym_type), POINTER                 :: csym
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_crys_symmetry', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, indpg, j, k, l, n, plevel, &
                                                unit
    REAL(KIND=dp)                            :: vs, vv0(3)

    unit = csym%punit
    IF(unit < 0) RETURN
    plevel = csym%plevel
    indpg=csym%indpg

    vs = 0.0_dp
    DO n = 1, csym%nc
       DO i = 1, 3
          vs = vs + ABS(csym%v(i,n))
       END DO
    END DO

    ! do be determined
    csym%xb=0._dp

    WRITE (unit,"(/,T2,A)") "Crystal Symmetry Information"

    WRITE (unit, '(A,T62,A,A)') "       Point group of the primitive lattice: ",&
         ADJUSTR(icst(csym%ihg))," system"
    IF (csym%isy==0) THEN
       WRITE (unit,'(T62,"nonsymmorphic group")')
    ELSE IF (csym%isy==1) THEN
       WRITE (unit,'(T65,"symmorphic group")')
    ELSE IF (csym%isy==-1) THEN
       WRITE (unit,'(T40,"symmorphic group with non-standard origin")')
       WRITE (unit,'(T7,A,A,/,T18,3F10.6,3X,3F10.6)') &
            ' The standard origin of coordinates is:   ', '[CARTESIAN]   [CRYSTAL]', csym%xb, csym%origin
    ELSE IF (csym%isy==-2) THEN
       WRITE (unit,'(T24,A)') ' Cannot determine if the space group is symmorphic or not'
       WRITE (unit,'(T7,A,/,T7,A,/,T7,A,F15.6,A)') ' The space group is non-symmorphic,', &
            ' or else a non standard origin of coordinates was used.',' (sum of translation vectors=', vs, ')'
    END IF
    IF (csym%li==0) THEN
       WRITE (unit,'(T60,"no inversion symmetry")')
    ELSE IF (csym%li>0) THEN
       WRITE (unit,'(T63,"inversion symmetry")')
    END IF

    IF (indpg>0) THEN
       CALL xstring(pgrp(indpg),i,j)
       CALL xstring(pgrd(indpg),k,l)
       WRITE (unit,'(T7,A,T60,"[INDEX=",I2,"]",T76,A5)') &
            ' The point group of the crystal is [Hermann-Mauguin] ',indpg,ADJUSTL(pgrp(indpg)(i:j))
       WRITE (unit,'(T41,A,T78,A3)') ' [Schoenflies] ',ADJUSTL(pgrd(indpg)(k:l))
    ELSE
       CALL xstring(pgrp(-indpg),i,j)
       CALL xstring(pgrd(-indpg),k,l)
       WRITE (unit,'(T7,A,I2,T40,A,A,"(",A,")",T71,"[INDEX=",I2,"]")') &
            '  Group order=', csym%nc, '    Subgroup of ', &
            pgrp(-indpg) (i:j), pgrd(-indpg) (k:l), -indpg
    END IF
    IF (csym%ntvec>1) THEN
       WRITE (unit,'(T7,A,T75,I6)') ' Number of primitive cells:',csym%ntvec
    END IF

    WRITE (unit,'(7X,A,T78,I3)') &
         "Total number of elements in the point group:",csym%nc
    WRITE (unit,'(7X,"to sum up: (",T64,I1,5I3,")")') &
         csym%ihg, csym%ihc, csym%isy, csym%li, csym%nc, csym%indpg
    WRITE (unit,'(/,7X,"List of the rotations:")')
    WRITE (unit,'(7X,T33,12I4)') (csym%ib(i),i=1,csym%nc)

    IF (csym%isy<=0) THEN
       WRITE (unit,'(7X,"Nonprimitive translations:")')
       WRITE (unit,'(T7,A,A)') ' ROT    V  in the basis A1, A2, A3      ', &
            'V  in cartesian coordinates'
       ! Cartesian components of nonprimitive translation.
       DO i = 1, csym%nc
          DO j = 1, 3
             vv0(j) = csym%v(1,i)*csym%a1(j) + csym%v(2,i)*csym%a2(j) + csym%v(3,i)*csym%a3(j)
          END DO
          WRITE (unit,'(7X,I6,3F10.5,7X,3F10.5)') csym%ib(i), (csym%v(j,i),j=1,3), vv0
       END DO
    END IF

    IF ( MOD(plevel,10) /= 0 ) THEN
       WRITE (unit,'(7X,"List of the 3 X 3 rotation matrices:")')
       IF (csym%ihc==0) THEN
          DO k = 1, csym%nc
             WRITE (unit, '(17X,I3," (",I2,": ",A11,")",2(3F14.6,/,38X),3F14.6)') k, &
                  csym%ib(k), rname_hexa(csym%ib(k)), ((csym%r(i,j,k),j=1,3),i=1,3)
          END DO
       ELSE
          DO k = 1, csym%nc
             WRITE (unit, '(17X,I3," (",I2,": ",A10,") ",2(3F14.6,/,38X),3F14.6)') k, &
                  csym%ib(k), rname_cubic(csym%ib(k)), ((csym%r(i,j,k),j=1,3),i=1,3)
          END DO
       END IF
    END IF

  END SUBROUTINE print_crys_symmetry

! *****************************************************************************
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
!> \param ihg ...
!> \param ihc ...
!> \param isy ...
!> \param li ...
!> \param nc ...
!> \param indpg ...
!> \param ib ...
!> \param v ...
!> \param r ...
!> \param xb ...
!> \param origin ...
!> \param wvkl ...
!> \param lwght ...
!> \param lrot ...
!> \param nhash ...
!> \param includ ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
!> \param error ...
!> \par Original Header:
!> ---------------------------------------------------------------------------!
!>      WRITTEN ON SEPTEMBER 12TH, 1979.
!>      IBM-RETOUCHED ON OCTOBER 27TH, 1980.
!>      GENERATION OF SPECIAL POINTS MODIFIED ON 26-MAY-82 BY OHN.
!>      RETOUCHED ON JANUARY 8TH, 1997
!>      INTEGRATION IN CPMD-FEMD PROGRAM BY THIERRY DEUTSCH
!> ---------------------------------------------------------------------------!
!>      PLAYING WITH SPECIAL POINTS AND CREATION OF 'CRYSTALLOGRAPHIC'
!>      FILE FOR BAND STRUCTURE CALCULATIONS.
!>      GENERATION OF SPECIAL POINTS IN K-SPACE FOR AN ARBITRARY LATTICE,
!>      FOLLOWING THE METHOD MONKHORST,PACK, PHYS. REV. B13 (1976) 5188
!>      MODIFIED BY MACDONALD, PHYS. REV. B18 (1978) 5897
!>      MODIFIED ALSO BY OLE HOLM NIELSEN ("SYMMETRIZATION")
!> ---------------------------------------------------------------------------!
!>      TESTING THEIR EFFICIENCY AND PREPARATION OF THE
!>      "STRUCTURAL" FILE FOR RUNNING THE
!>      SELF-CONSISTENT BAND STRUCTURE PROGRAMS.
!>      IN THE CASES WHERE THE POINT GROUP OF THE CRYSTAL DOES NOT
!>      CONTAIN INVERSION, THE LATTER IS ARTIFICIALLY ADDED, IN ORDER
!>      TO MAKE USE OF THE HERMITICITY OF THE HAMILTONIAN
!> ---------------------------------------------------------------------------!
!>         INPUT:
!>           IOUT    LOGIC FILE NUMBER
!>           NAT     NUMBER OF ATOMS
!>           NKPOINT MAXIMAL NUMBER OF K POINTS
!>           NSP     NUMBER OF SPECIES
!>           IQ1,IQ2,IQ3 THE MONKHORST-PACK MESH PARAMETERS
!>           ISTRIZ  SWITCH FOR SYMMETRIZATION
!>           A1(3),A2(3),A3(3) LATTICE VECTORS
!>           ALAT    LATTICE CONSTANT
!>           STRAIN(3,3)  STRAIN APPLIED TO LATTICE IN ORDER
!>                   TO HAVE K POINTS WITH SYMMETRY OF STRAINED LATTICE
!>           XKAPA(3,NAT)   ATOMS COORDINATES
!>           TY(NAT)      TYPES OF ATOMS
!>           WVK0(3) SHIFT FOR K POINTS MESh (MACDONALD ARTICLE)
!>           NHASH   SIZE OF THE HASH TABLES (LIST)
!>           DELTA   REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>                   K-VECTOR < DELTA IS CONSIDERED ZERO
!>         OUTPUT:
!>           RX(3,NAT) SCRATCH ARRAY USED BY GROUP1 ROUTINE
!>           TVEC(1:3,1:NTVEC) TRANSLATION VECTORS (SEE NTVEC)
!>           ISC(NAT)  SCRATCH ARRAY USED BY GROUP1 ROUTINE
!>           F0(49,NAT) ATOM TRANSFORMATION TABLE
!>                      IF NTVEC/=1 THE 49TH GIVES INEQUIVALENT ATOMS
!>           NTVEC NUMBER OF TRANSLATION VECTORS (IF NOT PRIMITIVE CELL)
!>           WVKL(3,NKPOINT) SPECIAL KPOINTS GENERATED
!>           LWGHT(NKPOINT)  WEIGHT FOR EACH K POINT
!>           LROT(48,NKPOINT) SYMMETRY OPERATION FOR EACh K POINTS
!>           INCLUD(NKPOINT)  SCRATCH ARRAY USED BY SPPT2
!>           LIST(NKPOINT+NHASH) HASH TABLE USED BY SPPT2
!>           RLIST(3,NKPOINT) SCRATCH ARRAY USED BY SPPT2
!> ---------------------------------------------------------------------------!
!>      SUBROUTINES NEEDED:
!>         SPPT2, GROUP1, PGL1, ATFTM1, ROT1, STRUCT,
!>         BZRDUC, INBZ, MESH, BZDEFI
!>      (GROUP1, PGL1, ATFTM1, ROT1 FROM THE
!>       "COMPUTER PHYSICS COMMUNICATIONS" PACKAGE "ACMI" - (1971,1974)
!>       WORLTON-WARREN).
! *****************************************************************************
  SUBROUTINE k290prg(iout,nat,nkpoint,nsp,iq1,iq2,iq3,istriz,a1,a2,a3, &
       alat,strain,xkapa,rx,tvec,ty,isc,f0,ntvec,wvk0,ihg,ihc,isy,li,nc,indpg,ib,v,r,&
       xb,origin,wvkl,lwght,lrot,nhash,includ,list,rlist,delta,error)
    INTEGER                                  :: iout, nat, nkpoint, nsp, iq1, &
                                                iq2, iq3, istriz
    REAL(KIND=dp)                            :: a1(3), a2(3), a3(3), alat, &
                                                strain(6), xkapa(3,nat), &
                                                rx(3,nat), tvec(3,nat)
    INTEGER                                  :: ty(nat), isc(nat), &
                                                f0(49,nat), ntvec
    REAL(KIND=dp)                            :: wvk0(3)
    INTEGER                                  :: ihg, ihc, isy, li, nc, indpg, &
                                                ib(48)
    REAL(KIND=dp)                            :: v(3,48), r(3,3,48), xb(3), &
                                                origin(3), wvkl(3,nkpoint)
    INTEGER                                  :: lwght(nkpoint), &
                                                lrot(48,nkpoint), nhash, &
                                                includ(nkpoint), &
                                                list(nkpoint+nhash)
    REAL(KIND=dp)                            :: rlist(3,nkpoint), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'k290prg', &
      routineP = moduleN//':'//routineN

    INTEGER :: f00(49,1), i, ib0(48), ihc0, ihg0, indpg0, invadd, istrin, &
      iswght, isy0, j, l, li0, lmax, n, nc0, ntot, ntvec0
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: a01(3), a02(3), a03(3), b01(3), b02(3), b03(3), b1(3), &
      b2(3), b3(3), dtotstr, origin0(3), proj1, proj2, proj3, r0(3,3,48), &
      totstr, tvec0(3,1), v0(3,48), volum, x0(3,1), xb0(3)

    failure = .FALSE.
    v(:,:) = 0.0_dp
    f00 = 0
    x0 = 0.0_dp
    v0 = 0.0_dp
    !----------------------------------------------------------------------------!
    !   read in lattice structure
    !----------------------------------------------------------------------------!
    DO i = 1, 3
       a01(i) = a1(i)/alat
       a02(i) = a2(i)/alat
       a03(i) = a3(i)/alat
    END DO
    !----------------------------------------------------------------------------!
    !   is the strain significant ?
    !----------------------------------------------------------------------------!
    dtotstr = delta*delta
    totstr = 0.0_dp
    istrin = 0
    DO i = 1, 6
       totstr = totstr + ABS(strain(i))
    END DO
    IF (totstr>dtotstr) istrin = 1
    !----------------------------------------------------------------------------!
    !   Volume of the cell.
    volum = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) + a3(1)*a1(2)*a2(3) - &
         a1(3)*a2(2)*a3(1) - a2(3)*a3(2)*a1(1) - a3(3)*a1(2)*a2(1)
    volum = ABS(volum)
    b1(1) = (a2(2)*a3(3)-a2(3)*a3(2))/volum
    b1(2) = (a2(3)*a3(1)-a2(1)*a3(3))/volum
    b1(3) = (a2(1)*a3(2)-a2(2)*a3(1))/volum
    b2(1) = (a3(2)*a1(3)-a3(3)*a1(2))/volum
    b2(2) = (a3(3)*a1(1)-a3(1)*a1(3))/volum
    b2(3) = (a3(1)*a1(2)-a3(2)*a1(1))/volum
    b3(1) = (a1(2)*a2(3)-a1(3)*a2(2))/volum
    b3(2) = (a1(3)*a2(1)-a1(1)*a2(3))/volum
    b3(3) = (a1(1)*a2(2)-a1(2)*a2(1))/volum
    !----------------------------------------------------------------------------!
    !      GROUP-THEORY ANALYSIS OF LATTICE
    !----------------------------------------------------------------------------!
    CALL group1(iout,a1,a2,a3,nat,ty,xkapa,b1,b2,b3,ihg,ihc,isy,li,nc, &
         indpg,ib,ntvec,v,f0,xb,r,tvec,origin,rx,isc,delta,error)
    !----------------------------------------------------------------------------!
    DO n = nc + 1, 48
       ib(n) = 0
    END DO
    !----------------------------------------------------------------------------!
    IF (ABS(iq1)+ABS(iq2)+ABS(iq3)/=0) THEN

       invadd = 0
       IF (li==0) invadd = 1
       IF (li==0 .AND. iout>0) THEN
          WRITE (iout,'(A)') &
               ' K290| Although the point group of the crystal does not'
          WRITE (iout,'(A,A)') ' K290| contain inversion, ', &
               'the special point generation algorithm'
          WRITE (iout,'(A)') ' K290| will consider it as a symmetry operation'
          invadd = 1
       END IF
       !----------------------------------------------------------------------------!
       !   GENERATE THE BRAVAIS LATTICE
       !----------------------------------------------------------------------------!
       IF(iout>0) WRITE (iout,'(/,1X,20("-"),A,20("-"),/,A)') &
                      ' The (unstrained) bravais lattice ', &
                      ' (used for generating the largest possible mesh in the B.Z.)'
       !----------------------------------------------------------------------------!
       CALL group1(iout,a01,a02,a03,1,ty,x0,b01,b02,b03,ihg0,ihc0,isy0,li0, &
            nc0,indpg0,ib0,ntvec0,v0,f00,xb0,r0,tvec0,origin0,rx,isc,delta,error)
       !----------------------------------------------------------------------------!
       !   It is assumed that the same 'type' of symmetry operations
       !   (cubic/hexagonal) will apply to the crystal as well as the Bravais
       !   lattice.
       !----------------------------------------------------------------------------!
       IF(iout>0) THEN
          WRITE (iout,'(/,1X,19("*"),A,25("*"))') &
            ' Generation of special points '
          !   Parameter Q of Monkhorst and Pack, generalized for 3 axes B1,2,3
          WRITE (iout,'(A,/,1X,3I5)') &
            ' Monkhorst-Pack parameters (generalized) IQ1,IQ2,IQ3:', iq1, iq2, &
            iq3
          !   WVK0 is the shift of the whole mesh (see Macdonald)
          WRITE (iout,'(A,/,1X,3F10.5)') &
            ' constant vector shift (MacDonald) of this mesh:', wvk0
          IF (ABS(istriz)/=1) THEN
             WRITE (iout,'(" invalid switch for symmetrization",I10)') istriz
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
          WRITE (iout,'(" symmetrization switch: ",I3)',ADVANCE="NO") istriz
          IF (istriz==1) THEN
             WRITE (iout,'(" (symmetrization of Monkhorst-Pack mesh)")')
          ELSE
             WRITE (iout,'(" (no symmetrization of Monkhorst-Pack mesh)")')
          END IF
       END IF
       !   Set to 0.
       DO i = 1, nkpoint
          lwght(i) = 0
       END DO
       !----------------------------------------------------------------------------!
       !      Generation of the points (they are not multiplied
       !      by 2*Pi because B1,2,3  were not,either)
       !----------------------------------------------------------------------------!
       IF (nc>nc0) THEN
          !     Due to non-use of primitive cell, the crystal has more
          !     rotations than Bravais lattice.
          !     We use only the rotations for Bravais lattices
          IF (ntvec==1 .AND. iout>0) THEN
             WRITE (iout,*) ' K290| Number of rotations for bravais lattice', &
                  nc0
             WRITE (iout,*) ' K290| Number of rotations for crystal lattice', &
                  nc
             WRITE (iout,*) ' K290| No duplication found'
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
          nc = nc0
          DO i = 1, nc0
             ib(i) = ib0(i)
          END DO
          IF(iout>0) THEN
             WRITE (iout,'(/,1X,20("!"),"WARNING",20("!"))')
             WRITE (iout,'(A)') &
                  ' K290| The crystal has more symmetry than the bravais lattice'
             WRITE (iout,'(A)') ' K290| Because this is not a primitive cell'
             WRITE (iout,'(A)') ' K290| Use only symmetry from bravais lattice'
             WRITE (iout,'(1X,20("!"),"WARNING",20("!"),/)')
          END IF
       END IF
       CALL sppt2(iout,iq1,iq2,iq3,wvk0,nkpoint,a01,a02,a03,b01,b02,b03, &
            invadd,nc,ib,r,ntot,wvkl,lwght,lrot,nc0,ib0,istriz,nhash,includ, &
            list,rlist,delta,error)
       !----------------------------------------------------------------------------!
       !      Check on error signals
       !----------------------------------------------------------------------------!
       IF(iout>0) WRITE (iout,'(/,1X,I5," Special points generated")') ntot
       IF (ntot/=0) THEN
          IF (ntot<0) THEN
             IF(iout>0) WRITE (iout,'(A,I5,/,A,/,A)') ' Dimension nkpoint =', nkpoint, &
                  ' insufficient for accommodating all the special points', &
                  ' what follows is an incomplete list'
             ntot = iabs(ntot)
          END IF
          !   Before using the list WVKL as wave vectors, they have to be
          !   multiplied by 2*Pi
          !   The list of weights LWGHT is not normalized
          iswght = 0
          DO i = 1, ntot
             iswght = iswght + lwght(i)
          END DO
          IF(iout>0) WRITE (iout,'(8X,A,T33,A,4X,A)') 'WAVEVECTOR K', 'WEIGHT', &
               'UNFOLDING ROTATIONS'
          !   Set near-zeroes equal to zero:
          DO l = 1, ntot
             DO i = 1, 3
                IF (ABS(wvkl(i,l))<delta) wvkl(i,l) = 0.0_dp
             END DO
             IF (istrin/=0) THEN
                !       Express special points in (unstrained) basis.
                proj1 = 0.0_dp
                proj2 = 0.0_dp
                proj3 = 0.0_dp
                DO i = 1, 3
                   proj1 = proj1 + wvkl(i,l)*a01(i)
                   proj2 = proj2 + wvkl(i,l)*a02(i)
                   proj3 = proj3 + wvkl(i,l)*a03(i)
                END DO
                DO i = 1, 3
                   wvkl(i,l) = proj1*b1(i) + proj2*b2(i) + proj3*b3(i)
                END DO
             END IF
             lmax = lwght(l)
             IF(iout>0) THEN
                WRITE (iout,fmt='(1X,I5,3F8.4,I8,T42,12I3)') l, &
                     (wvkl(i,l),i=1,3), lwght(l), (lrot(i,l),i=1,MIN(lmax,12))
                DO j = 13, lmax, 12
                   WRITE (iout,fmt='(T42,12I3)') (lrot(i,l),i=j,MIN(lmax,j-1+12))
                END DO
             END IF
          END DO
          IF(iout>0) WRITE (iout,'(24X,"TOTAL:",I8)') iswght
       END IF
    END IF
  END SUBROUTINE k290prg

! *****************************************************************************
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
!> \param xb ...
!> \param r ...
!> \param tvec ...
!> \param origin ...
!> \param rx ...
!> \param isc ...
!> \param delta ...
!> \param error ...
!> \par Original Header:
!> ---------------------------------------------------------------------------!
!>       WRITTEN ON SEPTEMBER 10TH - FROM THE ACMI COMPLEX
!>       (WORLTON AND WARREN, COMPUT.PHYS.COMMUN. 8,71-84 (1974))
!>       (AND 3,88-117 (1972))
!>       BASIC CRYSTALLOGRAPHIC INFORMATION
!>       ABOUT A GIVEN CRYSTAL STRUCTURE.
!>       SUBROUTINES NEEDED: PGL1,ATFTM1,ROT1,RLV3
!> ---------------------------------------------------------------------------!
!>       INPUT DATA:
!>       IOUT ... NUMBER OF THE OUTPUT UNIT FOR ON-LINE PRINTING
!>                OF VARIOUS MESSAGES
!>                IF IOUT.LE.0 NO MESSAGE
!>       A1,A2,A3 .. ELEMENTARY TRANSLATIONS OF THE LATTICE, IN SOME
!>                UNIT OF LENGTH
!>       NAT .... NUMBER OF ATOMS IN THE UNIT CELL
!>                ALL THE DIMENSIONS ARE SET FOR NAT .LE. 20
!>       TY ..... INTEGERS DISTINGUISHING BETWEEN THE ATOMS OF
!>                DIFFERENT TYPE. TY(I) IS THE TYPE OF THE I-TH ATOM
!>                OF THE BASIS
!>       X ...... CARTESIAN COORDINATES OF THE NAT ATOMS OF THE BASIS
!>       DELTA... REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!> ---------------------------------------------------------------------------!
!>       OUTPUT DATA:
!>       B1,B2,B3 .. RECIPROCAL LATTICE VECTORS, NOT MULTIPLIED BY
!>                ANY 2PI, IN UNITS RECIPROCAL TO THOSE OF A1,A2,A3
!>       IHG .... POINT GROUP OF THE PRIMITIVE LATTICE, HOLOHEDRAL
!>                GROUP NUMBER:
!>                IHG=1 STANDS FOR TRICLINIC    SYSTEM
!>                IHG=2 STANDS FOR MONOCLINIC   SYSTEM
!>                IHG=3 STANDS FOR ORTHORHOMBIC SYSTEM
!>                IHG=4 STANDS FOR TETRAGONAL   SYSTEM
!>                IHG=5 STANDS FOR CUBIC        SYSTEM
!>                IHG=6 STANDS FOR TRIGONAL     SYSTEM
!>                IHG=7 STANDS FOR HEXAGONAL    SYSTEM
!>       IHC .... CODE DISTINGUISHING BETWEEN HEXAGONAL AND CUBIC
!>                GROUPS
!>                IHC=0 STANDS FOR HEXAGONAL GROUPS
!>                IHC=1 STANDS FOR CUBIC GROUPS
!>       ISY .... CODE INDICATING WHETHER THE SPACE GROUP IS
!>                SYMMORPHIC OR NONSYMMORPHIC
!>                ISY= 0 NONSYMMORPHIC GROUP
!>                ISY= 1 SYMMORPHIC GROUP
!>                ISY=-1 SYMMORPHIC GROUP WITH NON-STANDARD ORIGIN
!>                ISY=-2 UNDETERMINED (NORMALLY NEVER)
!>                THE GROUP IS CONSIDERED SYMMORPHIC IF FOR EACH
!>                OPERATION OF THE POINT GROUP THE SUM OF THE 3
!>                COMPONENTS OF ABS(V(N)) (NONPRIMITIVE TRANSLATION,
!>                SEE BELOW) IS LT. 0.0001
!>       ORIGIN   STANDARD ORIGIN IF SYMMORPHIC (CRYSTAL COORDINATES)
!>       LI ..... CODE INDICATING WHETHER THE POINT GROUP
!>                OF THE CRYSTAL CONTAINS INVERSION OR NOT
!>                (OPERATIONS 13 OR 25 IN RESPECTIVELY HEXAGONAL
!>                OR CUBIC GROUPS).
!>                LI=0 MEANS: DOES NOT CONTAIN INVERSION
!>                LI.GT.0 MEANS: THERE IS INVERSION IN THE POINT
!>                               GROUP OF THE CRYSTAL
!>       NC ..... TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP OF THE
!>                CRYSTAL
!>       INDPG .. POINT GROUP INDEX (DETERMINED IF SYMMORPHIC GROUP)
!>       IB ..... LIST OF THE ROTATIONS CONSTITUTING THE POINT GROUP
!>                OF THE CRYSTAL. THE NUMBERING IS THAT DEFINED IN
!>                WORLTON AND WARREN, I.E. THE ONE MATERIALIZED IN THE
!>                ARRAY R (SEE BELOW)
!>                ONLY THE FIRST NC ELEMENTS OF THE ARRAY IB ARE
!>                MEANINGFUL
!>       NTVEC .. NUMBER OF TRANSLATIONAL VECTORS
!>                ASSOCIATED WITH IDENTITY OPERATOR I.E.
!>                GIVES THE NUMBER OF IDENTICAL PRIMITIVE CELLS
!>       V ...... NONPRIMITIVE TRANSLATIONS (IN THE CASE OF NONSYMMOR-
!>                PHIC GROUPS). V(I,N) IS THE I-TH COMPONENT
!>                OF THE TRANSLATION CONNECTED WITH THE N-TH ELEMENT
!>                OF THE POINT GROUP (I.E. WITH THE ROTATION
!>                NUMBER IB(N) ).
!>                ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS,
!>                THEY REFER TO THE SYSTEM A1,A2,A3.
!>       F0 ..... THE FUNCTION DEFINED IN MARADUDIN, IPATOVA BY
!>                EQ. (3.2.12): ATOM TRANSFORMATION TABLE.
!>                THE ELEMENT F0(N,KAPA) MEANS THAT THE N-TH
!>                OPERATION OF THE SPACE GROUP (I.E. OPERATION NUMBER
!>                IB(N), TOGETHER WITH AN EVENTUAL NONPRIMITIVE
!>                TRANSLATION  V(N)) TRANSFERS THE ATOM KAPA INTO THE
!>                ATOM F0(N,KAPA).
!>                THE 49TH LINE GIVES EQUIVALENT ATOMS FOR
!>                FRACTIONAl TRANSLATIONS ASSOCIATED WITH IDENTITY
!>       R ...... LIST OF THE 3 X 3 ROTATION MATRICES
!>                (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)
!>                ALL 48 OR 24 MATRICES ARE LISTED.
!>                FOLLOW NOTATION OF WORLTON-WARREN(1972)
!>       TVEC  .. LIST OF NTVEC TRANSLATIONAL VECTORS
!>                ASSOCIATED WITH IDENTITY OPERATOR
!>                TVEC(1:3,1) = \(0,0,0\)
!>                (CRYSTAL COORDINATES)
!>       RX ..... SCRATCH ARRAY
!>       ISC .... SCRATCH ARRAY
!> ---------------------------------------------------------------------------!
!>       PRINTED OUTPUT:
!>       PROGRAM PRINTS THE TYPE OF THE LATTICE (IHG, IN WORDS),
!>       LISTS THE OPERATIONS OF THE  POINT GROUP OF THE
!>       CRYSTAL, INDICATES WHETHER THE SPACE GROUP IS SYMMORPHIC OR
!>       NONSYMMORPHIC AND WHETHER THE POINT GROUP OF THE CRYSTAL
!>       CONTAINS INVERSION.
! *****************************************************************************
  SUBROUTINE group1(iout,a1,a2,a3,nat,ty,x,b1,b2,b3,ihg,ihc,isy,li,nc, &
       indpg,ib,ntvec,v,f0,xb,r,tvec,origin,rx,isc,delta,error)
    INTEGER                                  :: iout
    REAL(KIND=dp)                            :: a1(3), a2(3), a3(3)
    INTEGER                                  :: nat, ty(nat)
    REAL(KIND=dp)                            :: x(3,nat), b1(3), b2(3), b3(3)
    INTEGER                                  :: ihg, ihc, isy, li, nc, indpg, &
                                                ib(48), ntvec
    REAL(KIND=dp)                            :: v(3,48)
    INTEGER                                  :: f0(49,nat)
    REAL(KIND=dp)                            :: xb(3), r(3,3,48), &
                                                tvec(3,nat), origin(3), &
                                                rx(3,nat)
    INTEGER                                  :: isc(nat)
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'group1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ncprim
    REAL(KIND=dp)                            :: a(3,3), ai(3,3), ap(3,3), &
                                                api(3,3)

    DO i = 1, 3
       a(i,1) = a1(i)
       a(i,2) = a2(i)
       a(i,3) = a3(i)
    END DO
    !------------------------------------------------------------------------------!
    !        A(I,J) IS THE I-TH CARTESIAN COMPONENT OF THE J-TH PRIMITIVE
    !        TRANSLATION VECTOR OF THE DIRECT LATTICE
    !        TY(I) IS AN INTEGER DISTINGUISHING ATOMS OF DIFFERENT TYPE,
    !        I.E., DIFFERENT ATOMIC SPECIES
    !        X(J,I) IS THE J-TH CARTESIAN COMPONENT OF THE POSITION
    !        VECTOR FOR THE I-TH ATOM IN THE UNIT CELL.
    !------------------------------------------------------------------------------!
    !     DETERMINE PRIMITIVE LATTICE VECTORS FOR THE RECIPROCAL LATTICE
    !------------------------------------------------------------------------------!

    CALL calbrec(a,ai,error)
    DO i = 1, 3
       b1(i) = ai(1,i)
       b2(i) = ai(2,i)
       b3(i) = ai(3,i)
    END DO
    !------------------------------------------------------------------------------!
    !     Determination of the translation vectors associated with
    !     the Identity matrix i.e. if the cell is duplicated
    !     Give also the ``primitive lattice''
    CALL primlatt(a,ai,ap,api,nat,ty,x,ntvec,tvec,f0,isc,delta,error)
    !------------------------------------------------------------------------------!
    !     Determination of the holohedral group (and crystal system)
    CALL pgl1(ap,api,ihc,nc,ib,ihg,r,delta,error)
    IF (ntvec>1) THEN
       !       All rotations found by PGL1 have axes in x, y or z cart. axis
       !       So we have too check if we do not loose symmetry
       ncprim = nc
       !       The hexagonal system is found if the z axis is the sixfold axis
       CALL pgl1(a,ai,ihc,nc,ib,ihg,r,delta,error)
       IF (ncprim>nc) THEN
          !         More symmetry with
          CALL pgl1(ap,api,ihc,nc,ib,ihg,r,delta,error)
       END IF
    END IF
    !------------------------------------------------------------------------------!
    !     Determination of the space group
    CALL atftm1(iout,r,v,x,f0,xb,origin,ib,ty,nat,ihg,rx,nc,indpg,ntvec,a,ai, &
         li,isy,isc,delta,error)
    !------------------------------------------------------------------------------!
  END SUBROUTINE group1

! *****************************************************************************
!> \brief CALCULATE RECIPROCAL VECTOR BASIS (AI(1:3,1:3))
!> \param a ...
!> \param ai ...
!> \param error ...
!> \par   Note
!>             INPUT:
!>              A(3,3) A(I,J) IS THE I-TH CARTESIAN COMPONENT
!>                OF THE J-TH PRIMITIVE TRANSLATION VECTOR OF
!>                THE DIRECT LATTICE
!>             OUTPUT:
!>              AI(3,3) RECIPROCAL VECTOR BASIS
! *****************************************************************************
  SUBROUTINE calbrec(a,ai,error)
    REAL(KIND=dp)                            :: a(3,3), ai(3,3)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calbrec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, il, iu, j, jl, ju
    REAL(KIND=dp)                            :: det

    det = a(1,1)*a(2,2)*a(3,3) + a(2,1)*a(1,3)*a(3,2) + &
         a(3,1)*a(1,2)*a(2,3) - a(1,1)*a(2,3)*a(3,2) - a(2,1)*a(1,2)*a(3,3) - &
         a(3,1)*a(1,3)*a(2,2)
    det = 1.E0_dp/det
    DO i = 1, 3
       il = 1
       iu = 3
       IF (i==1) il = 2
       IF (i==3) iu = 2
       DO j = 1, 3
          jl = 1
          ju = 3
          IF (j==1) jl = 2
          IF (j==3) ju = 2
          ai(j,i) = (-1.0_dp)**(i+j)*det*(a(il,jl)*a(iu,ju)-a(il,ju)*a(iu,jl &
               ))
       END DO
    END DO
  END SUBROUTINE calbrec

! *****************************************************************************
!> \brief DETERMINATION OF THE TRANSLATION VECTORS ASSOCIATED WITH
!>        THE IDENTITY SYMMETRY I.E. IF THE CELL IS DUPLICATED
!>        GIVE ALSO THE PRIMITIVE DIRECT AND RECIPROCAL LATTICE VECTOR
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
!> \param error ...
!> \par Note
!>        INPUT:
!>          A(3,3)   A(I,J) IS THE I-TH CARTESIAN COMPONENT
!>                   OF THE J-TH TRANSLATION VECTOR OF
!>                   THE DIRECT LATTICE
!>          AI(3,3)  RECIPROCAL VECTOR BASIS (CARTESIAN)
!>          NAT      NUMBER OF ATOMS
!>          TY(NAT)  TYPE OF ATOMS
!>          X(3,NAT) ATOMIC COORDINATES IN CARTESIAN COORDINATES
!>          DELTA    REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>        OUTPUT:
!>          AP(3,3)  COMPONENTS OF THE PRIMITIVE TRANSLATION VECTORS
!>          API(3,3) PRIMITIVE RECIPROCAL BASIS VECTORS
!>                   BOTH BAISI ARE IN CARTESIAN COORDINATES
!>          NTVEC    NUMBER OF TRANSLATION VECTORS (FRACTIONNAL)
!>          TVEC(3,NTVEC) COMPONENTS OF TRANSLATIONAL VECTORS
!>                        (CRYSTAL COORDINATES)
!>          F0(49,NAT)    GIVES INEQUIVALENT ATOM FOR EACH ATOM
!>                        THE 49-TH LINE
!>          ISC(NAT)      SCRATCH ARRAY
! *****************************************************************************
  SUBROUTINE primlatt(a,ai,ap,api,nat,ty,x,ntvec,tvec,f0,isc,delta,error)
    REAL(KIND=dp)                            :: a(3,3), ai(3,3), ap(3,3), &
                                                api(3,3)
    INTEGER                                  :: nat, ty(nat)
    REAL(KIND=dp)                            :: x(3,nat)
    INTEGER                                  :: ntvec
    REAL(KIND=dp)                            :: tvec(3,nat)
    INTEGER                                  :: f0(49,nat), isc(nat)
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'primlatt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, il, iv, j, k2
    LOGICAL                                  :: oksym
    REAL(KIND=dp)                            :: vr(3), xb(3)

    ntvec = 1
    tvec(1,1) = 0.0_dp
    tvec(2,1) = 0.0_dp
    tvec(3,1) = 0.0_dp
    !     First we check if there exist fractional translational vectors
    !     associated with Identity operation i.e.
    !     if the cell is duplicated or not.
    DO i = 1, nat
       f0(49,i) = i
    END DO
    DO k2 = 2, nat
       IF (ty(1)/=ty(k2)) GO TO 100
       DO i = 1, 3
          xb(i) = x(i,k2) - x(i,1)
       END DO
       !       A fractional translation vector VR is defined.
       CALL rlv3(ai,xb,vr,il,delta,error)
       CALL checkrlv3(1,nat,ty,x,x,vr,f0,ai,isc,.TRUE.,oksym,delta,error)
       IF (oksym) THEN
          !         A fractional translational vector is found
          ntvec = ntvec + 1
          !         F0(49,1:NAT) gives number of equivalent atoms
          !         and has atom indexes of inequivalent atoms (for translation)
          DO i = 1, nat
             IF (f0(49,i)>f0(1,i)) f0(49,i) = f0(1,i)
          END DO
          DO i = 1, 3
             tvec(i,ntvec) = vr(i)
          END DO
       END IF
100    CONTINUE
    END DO
    !------------------------------------------------------------------------------!
    DO i = 1, 3
       DO il = 1, 3
          ap(il,i) = a(il,i)
       END DO
    END DO
    DO j = 1, 3
       DO i = 1, 3
          api(i,j) = ai(i,j)
       END DO
    END DO
    IF (ntvec==1) THEN
       !       The current cell is definitely a primitive one
       !       Copy A and AI to AP and API
    ELSE
       !       We are looking for the primitive lattice vector basis set
       !       AP is our current lattice vector basis
       DO iv = 2, ntvec
          !         TVEC in cartesian coordinates
          DO i = 1, 3
             xb(i) = tvec(1,iv)*a(i,1) + tvec(2,iv)*a(i,2) + &
                  tvec(3,iv)*a(i,3)
          END DO
          !         We calculare TVEC in AP basis
          CALL rlv3(api,xb,vr,il,delta,error)
          DO i = 1, 3
             IF (ABS(vr(i))>delta) THEN
                il = NINT(1.0_dp/ABS(vr(i)))
                IF (il>1) THEN
                   !               We replace AP(1:3,I) by TVEC(1:3,IV)
                   DO j = 1, 3
                      ap(j,i) = xb(j)
                   END DO
                   !               Calculate new API
                   CALL calbrec(ap,api,error)
                   GO TO 200
                END IF
             END IF
          END DO
200       CONTINUE
       END DO
    END IF
  END SUBROUTINE primlatt

! *****************************************************************************
!> \brief ...
!> \param a ...
!> \param ai ...
!> \param ihc ...
!> \param nc ...
!> \param ib ...
!> \param ihg ...
!> \param r ...
!> \param delta ...
!> \param error ...
!> \par Original Header
!> ---------------------------------------------------------------------------!
!>         WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX
!>         AUXILIARY SUBROUTINE TO GROUP1
!>         SUBROUTINE PGL DETERMINES THE POINT GROUP OF THE LATTICE
!>         AND THE CRYSTAL SYSTEM.
!>         SUBROUTINES NEEDED: ROT1, RLV3
!> ---------------------------------------------------------------------------!
!>         WARNING: FOR THE HEXAGONAL SYSTEM, THE 3RD AXIS SUPPOSE
!>                  TO BE THE SIX-FOLD AXIS
!> ---------------------------------------------------------------------------!
!>         INPUT:
!>         A ..... DIRECT LATTICE VECTORS
!>         AI .... RECIPROCAL LATTICE VECTORS
!>         DELTA.. REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!> ---------------------------------------------------------------------------!
!>         OUTPUT:
!>         IHC .... CODE DISTINGUISHING BETWEEN HEXAGONAL AND CUBIC
!>                  GROUPS
!>                  IHC=0 STANDS FOR HEXAGONAL GROUPS
!>                  IHC=1 STANDS FOR CUBIC GROUPS
!>         NC .... NUMBER OF ROTATIONS IN THE POINT GROUP
!>         IB .... SET OF ROTATION
!>         IHG .... POINT GROUP OF THE PRIMITIVE LATTICE, HOLOHEDRAL
!>                  GROUP NUMBER:
!>                  IHG=1 STANDS FOR TRICLINIC    SYSTEM
!>                  IHG=2 STANDS FOR MONOCLINIC   SYSTEM
!>                  IHG=3 STANDS FOR ORTHORHOMBIC SYSTEM
!>                  IHG=4 STANDS FOR TETRAGONAL   SYSTEM
!>                  IHG=5 STANDS FOR CUBIC        SYSTEM
!>                  IHG=6 STANDS FOR TRIGONAL     SYSTEM
!>                  IHG=7 STANDS FOR HEXAGONAL    SYSTEM
!>         R ...... LIST OF THE 3 X 3 ROTATION MATRICES
!>                  (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)
!>                  ALL 48 OR 24 MATRICES ARE LISTED.
!>                  FOLLOW NOTATION OF WORLTON-WARREN(1972)
! *****************************************************************************
  SUBROUTINE pgl1(a,ai,ihc,nc,ib,ihg,r,delta,error)
    REAL(KIND=dp)                            :: a(3,3), ai(3,3)
    INTEGER                                  :: ihc, nc, ib(48), ihg
    REAL(KIND=dp)                            :: r(3,3,48), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pgl1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, lx, n, nr
    REAL(KIND=dp)                            :: tr, vr(3), xa(3)

    DO ihc = 0, 1
       !       IHC is 0 for hexagonal groups and 1 for cubic groups.
       IF (ihc==0) THEN
          nr = 24
       ELSE
          nr = 48
       END IF
       nc = 0
       !       Constructs rotation operations.
       CALL rot1(ihc,r,error)
       DO n = 1, nr
          ib(n) = 0
          !         Rotate the A1,2,3 vectors by rotation No. N
          DO k = 1, 3
             DO i = 1, 3
                xa(i) = 0.0_dp
                DO j = 1, 3
                   xa(i) = xa(i) + r(i,j,n)*a(j,k)
                END DO
             END DO
             CALL rlv3(ai,xa,vr,lx,delta,error)
             tr = 0.0_dp
             DO i = 1, 3
                tr = tr + ABS(vr(i))
             END DO
             !           If VR.ne.0, then XA cannot be a multiple of a lattice vector
             IF (tr>delta) GO TO 140
          END DO
          nc = nc + 1
          ib(nc) = n
140       CONTINUE
       END DO
       !------------------------------------------------------------------------------!
       !       IHG stands for holohedral group number.
       IF (ihc==0) THEN
          !         Hexagonal group:
          IF (nc==12) ihg = 6
          IF (nc>12) ihg = 7
          IF (nc>=12) RETURN
          !         Too few operations, try cubic group: (IHC=1,NR=48)
       ELSE
          !         Cubic group:
          IF (nc<4) ihg = 1
          IF (nc==4) ihg = 2
          IF (nc>4) ihg = 3
          IF (nc==16) ihg = 4
          IF (nc>16) ihg = 5
          RETURN
       END IF
    END DO
  END SUBROUTINE pgl1

! *****************************************************************************
!> \brief ...
!> \param ai ...
!> \param xb ...
!> \param vr ...
!> \param il ...
!> \param delta ...
!> \param error ...
!> \par Original Header
!> ---------------------------------------------------------------------------!
!>         WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX
!>         AUXILIARY SUBROUTINE TO GROUP1
!>         SUBROUTINE RLV REMOVES A DIRECT LATTICE VECTOR
!>         FROM XB LEAVING THE REMAINDER IN VR.
!>         IF A NONZERO LATTICE VECTOR WAS REMOVED, IL IS MADE NONZERO.
!>         VR STANDS FOR V-REFERENCE.
!> ---------------------------------------------------------------------------!
!>         INPUT:
!>            AI(I,J) ARE THE RECIPROCAL LATTICE VECTORS,
!>                    B(I) = AI(I,J),J=1,2,3
!>            XB(1:3) VECTOR IN CARTESIAN COORDINATES
!>         DELTA      REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>         OUTPUT:
!>            VR IS NOT GIVEN IN CARTESIAN COORDINATES BUT
!>               IN THE SYSTEM A1,A2,A3 (CRYSTAL COORDINATES)
!>               AND BETWEEN -1/2 AND 1/2
!>            IL ABS OF VR
!>         K.K., 23.10.1979
! *****************************************************************************
  SUBROUTINE rlv3(ai,xb,vr,il,delta,error)

    REAL(KIND=dp)                            :: ai(3,3), xb(3), vr(3)
    INTEGER                                  :: il
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rlv3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    REAL(KIND=dp)                            :: ts

    il = 0
    DO i = 1, 3
       vr(i) = 0.0_dp
    END DO
    ts = ABS(xb(1)) + ABS(xb(2)) + ABS(xb(3))
    IF (ts<=delta) RETURN
    DO i = 1, 3
       DO j = 1, 3
          vr(i) = vr(i) + ai(i,j)*xb(j)
       END DO
       il = il + NINT(ABS(vr(i)))
       !       Change in order to have correct determination of origin and
       !       symmorphic group (T.D 30/03/98)
       !        VR(I) = - MOD(REAL(VR(I),KIND=dp),1.0_dp)
       vr(i) = NINT(vr(i)) - vr(i)
    END DO
  END SUBROUTINE rlv3

! *****************************************************************************
!> \brief ...
!> \param iout ...
!> \param r ...
!> \param v ...
!> \param x ...
!> \param f0 ...
!> \param xb ...
!> \param origin ...
!> \param ib ...
!> \param ty ...
!> \param nat ...
!> \param ihg ...
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
!> \param error ...
!> \par Original Header
!> ---------------------------------------------------------------------------!
!>         WRITTEN ON SEPTEMBER 11TH, 1979 - FROM ACMI COMPLEX
!>         AUXILIARY SUBROUTINE TO GROUP1
!>         SUBROUTINE ATFTMT DETERMINES
!>                    THE POINT GROUP OF THE CRYSTAL,
!>                    THE ATOM TRANSFORMATION TABLE,F0,
!>                    THE FRACTIONAL TRANSLATIONS,V,
!>                    ASSOCIATED WITH EACH ROTATION.
!>         SUBROUTINES NEEDED: RLV3 CHECKRLV3 SYMMORPHIC STOPGM XSTRING
!>         MAY 14TH,1998: A LOT OF CHANGES (ARGUMENTS)
!>                        BETTER DETERMINATION OF V
!>         SEP 15TH,1998: DETERMINATION OF FRACTIONAL TRANSLATIONAL VEC.
!> ---------------------------------------------------------------------------!
!>         INPUT:
!>           IOUT Logical file number (output)
!>                If IOUT.LE.0 no message
!>           IHG  Holohedral group number (determined by PGL1)
!>           NC   Number of rotation operations
!>           NAT  Number of atoms (used in the routine)
!>           X    Coordinates of atoms (cartesian)
!>           TY   Type of atoms
!>           R    Sets of transformation operations (cartesian)
!>           IB   Index giving NC operations in R
!>           AI   Reciprocal lattice vectors
!>           NTVEC Number of translational vectors
!>                associated with Identity
!>                if primitive cell NTVEC=1, TVEC=(0,0,0)
!>         DELTA  REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>         OUTPUT:
!>           RX(3,NAT) Scratch array
!>           ISC(NAT)  Scratch array
!>           NC        is modified (number of symmetry operations)
!>           INDPG     Point group index
!>           V(3,48)   The fractional translations associated
!>                     with each rotation (crystal coordinates)
!>           F0(1:48,NAT)
!>                The atom transformation table for rotation (48,NAT)
!>           ORIGIN Standard origin if symmorphic (crystal coordinates)
!>           ISY  = 1 Isommorphic group
!>                =-1 Isommorphic group with non-standard origin
!>                = 0 Non-Isommorphic group
!>                =-2 Undetermined (normally never)
!>         LI ..... Code indicating whether the point group
!>                  of the crystal contains inversion or not
!>                  (operations 13 or 25 in respectively hexagonal
!>                  or cubic groups).
!>                  LI=0    : does not contain inversion
!>                  LI.GT.0 : there is inversion in the point
!>                            group of the crystal
! *****************************************************************************
  SUBROUTINE atftm1(iout,r,v,x,f0,xb,origin,ib,ty,nat,ihg,rx,nc,indpg,ntvec, &
       a,ai,li,isy,isc,delta,error)
    INTEGER                                  :: iout
    REAL(KIND=dp)                            :: r(3,3,48), v(3,48), xb(3), &
                                                origin(3)
    INTEGER                                  :: ib(48), nat, ty(nat), &
                                                f0(49,nat)
    REAL(KIND=dp)                            :: x(3,nat)
    INTEGER                                  :: ihg
    REAL(KIND=dp)                            :: rx(3,nat)
    INTEGER                                  :: nc, indpg, ntvec
    REAL(KIND=dp)                            :: a(3,3), ai(3,3)
    INTEGER                                  :: li, isy, isc(nat)
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'atftm1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iis(48), il, info, j, k, &
                                                k2, l, n, nca, ni
    LOGICAL                                  :: failure, nodupli, oksym
    REAL(KIND=dp)                            :: vc(3,48), vr(3), vs

    failure = .FALSE.
    nodupli = ntvec == 1
    nca = 0
    DO n = 1, 48
       iis(n) = 0
    END DO
    !     Calculate translational vector for each operation
    !     and atom transformation table.
    DO n = 1, nc
       l = ib(n)
       iis(l) = 1
       DO k = 1, nat
          DO i = 1, 3
             rx(i,k) = 0.0_dp
             DO j = 1, 3
                rx(i,k) = rx(i,k) + r(i,j,l)*x(j,k)
             END DO
          END DO
       END DO
       DO k = 1, 3
          vr(k) = 0.0_dp
       END DO
       !       First we determine for VR=(/0,0,0/)
       !       IMPORTANT IF NOT UNIQUE ATOMS FOR DETERMINATION OF SYMMORPHIC
       CALL checkrlv3(n,nat,ty,rx,x,vr,f0,ai,isc,nodupli,oksym,delta,error)
       IF (oksym) THEN
          GO TO 190
       END IF
       !       Now we try other possible VR
       !       F0(49,1:NAT) has only inequivalent atom indexes for translation
       DO k2 = 1, nat
          IF (f0(49,k2)<k2) GO TO 185
          IF (ty(1)/=ty(k2)) GO TO 185
          DO i = 1, 3
             xb(i) = rx(i,1) - x(i,k2)
          END DO
          !         A translation vector VR is defined.
          CALL rlv3(ai,xb,vr,il,delta,error)
          !------------------------------------------------------------------------------!
          !            SUBROUTINE RLV3 REMOVES A DIRECT LATTICE VECTOR FROM XB
          !            LEAVING THE REMAINDER IN VR. IF A NONZERO LATTICE
          !            VECTOR WAS REMOVED, IL IS MADE NONZERO.
          !            VR STANDS FOR V-REFERENCE.
          !            VR IS NOT GIVEN IN CARTESIAN COORDINATES BUT
          !            IN THE SYSTEM A1,A2,A3.     K.K., 23.10.1979
          !------------------------------------------------------------------------------!
          CALL checkrlv3(n,nat,ty,rx,x,vr,f0,ai,isc,nodupli,oksym,delta,error)
          IF (oksym) THEN
             GO TO 190
          END IF
185       CONTINUE
       END DO
       iis(l) = 0
       GO TO 210
190    CONTINUE
       nca = nca + 1
       DO i = 1, 3
          v(i,nca) = vr(i)
       END DO
       !------------------------------------------------------------------------------!
       !          V(I,N) IS THE I-TH COMPONENT OF THE FRACTIONAL
       !          TRANSLATION ASSOCIATED WITH THE ROTATION N.
       !          ATTENTION: V(I) ARE NOT CARTESIAN COMPONENTS, THEY ARE
       !          GIVEN IN THE SYSTEM A1,A2,A3.
       !          K.K., 23. 10. 1979
       !------------------------------------------------------------------------------!
210    CONTINUE
    END DO
    !     Remove unused operations
    i = 0
    ni = 13
    IF (ihg<6) ni = 25
    li = 0
    DO n = 1, nc
       l = ib(n)
       IF (iis(l)==0) GO TO 230
       i = i + 1
       ib(i) = ib(n)
       IF (ib(i)==ni) li = i
       DO k = 1, nat
          f0(i,k) = f0(n,k)
       END DO
230    CONTINUE
    END DO
    !------------------------------------------------------------------------------!
    nc = i
    vs = 0.0_dp
    DO n = 1, nc
       DO i = 1, 3
          vs = vs + ABS(v(i,n))
       END DO
    END DO
    !     THE ORIGINAL VALUE DELTA=0.0001 WAS MODIFIED
    !     BY K.K. , SEPTEMBER 1979 TO 0.0005
    !     AND RETURNED TO 0.0001 BY RJN OCT 1987
    IF (vs>delta) THEN
       isy = 0
    ELSE
       isy = 1
    END IF
    !------------------------------------------------------------------------------!
    !     Determination of the point group
    !     (Thierry Deutsch - 1998 [Maybe not complete!!])
    IF (ihg<6) THEN
       IF (nc==0) THEN
          WRITE (*,'(" ATFTM1! IHG=",A," NC=",I2)') icst(ihg), nc
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          !       Triclinic system
       ELSE IF (nc==1) THEN
          !         IB=1
          indpg = & !1 (c1)
               1
       ELSE IF (nc==2 .AND. ib(2)==25) THEN
          !         IB=1 25
          indpg = & !<1>(ci)
               2
       ELSE IF (nc==2 .AND. (ib(2)==4 .OR. & !         Monoclinic system
            ib(2)==2 .OR. & !         IB=1 4 (z-axis) OR
            ib(2)==3)) & !         IB=1 2 (x-axis) OR
            THEN
          !         IB=1 3 (y-axis)
          !2[001]
          !2[100]
          !2[010]
          indpg = & !2 (c2)
               3
       ELSE IF (nc==2 .AND. (ib(2)==28 .OR. ib(2)==26 .OR. ib(2)==27)) THEN
          !         IB=1 28 (z-axis) OR
          !         IB=1 26 (x-axis) OR
          !         IB=1 27 (y-axis)
          indpg = & !m (c1h)
               4
       ELSE IF (nc==4 .AND. (ib(4)==28 .OR. & !         IB=1  4 25 28 (z-axis)  OR
            ib(4)==27 .OR. & !         IB=1  2 25 26 (x-axis)  OR
            ib(4)==26 .OR. & !         IB=1  3 25 27 (y-axis)  OR
            ib(4)==37 .OR. & !         IB=1 13 25 37 (-xy-axis)OR
            ib(4)==40)) & !         IB=1 16 25 40 (xy-axis)
            THEN
          !2[001]
          !2[010]
          !2[100]
          !-2[-110]
          !2[110]
          indpg = & !2/m(c2h)
               5
       ELSE IF (nc==4 .AND. (ib(4)==15 .OR. ib(4)==20 .OR. ib(4)==24)) THEN
          !         Tetragonal system
          !         IB=1 4 14 15 (z-axis) OR
          !         IB=1 2 19 20 (x-axis) OR
          !         IB=1 3 22 24 (y-axis)
          indpg = & !4 (c4)
               11
       ELSE IF (nc==4 .AND. (ib(4)==39 .OR. ib(4)==44 .OR. ib(4)==48)) THEN
          !         IB=1 4 38 39 (z-axis) OR
          !         IB=1 2 43 44 (x-axis) OR
          !         IB=1 3 46 48 (y-axis)
          indpg = & !<4>(s4)
               12
       ELSE IF (nc==8 .AND. ((ib(3)==14 .AND. ib(8)==39) .OR. (ib(3)==19 &
            .AND. ib(8)==44) .OR. (ib(3)==22 .AND. ib(8)==48))) THEN
          !         IB=1 4 14 15 28 25 38 39 (z-axis) OR
          !         IB=1 2 19 20 26 25 43 44 (x-axis) OR
          !         IB=1 3 22 24 27 25 46 48 (y-axis)
          indpg = & !422(d4)
               13
       ELSE IF (nc==8 .AND. ib(4)==4 .AND. (ib(8)==16 .OR. ib( &
            8)==20 .OR. ib(8)==24)) THEN
          !         IB=1 2 3 4 13 14 15 16 (z-axis) OR
          !         IB=1 2 3 4 17 19 20 18 (x-axis) OR
          !         IB=1 2 3 4 21 22 24 23 (y-axis)
          indpg = & !4/m(c4h)
               14
       ELSE IF (nc==8 .AND. (ib(8)==40 .OR. ib(8)==42 .OR. ib(8)==47)) THEN
          !         IB=1 4 14 15 26 27 37 40 (z-axis) OR
          !         IB=1 2 19 20 28 27 41 42 (x-axis) OR
          !         IB=1 3 22 24 26 28 45 47 (y-axis)
          indpg = & !4mm(c4v)
               15
       ELSE IF (nc==8 .AND. ((ib(3)==13 .AND. ib(8)==39) .OR. (ib(3)==17 &
            .AND. ib(8)==44) .OR. (ib(3)==21 .AND. ib(8)==48))) THEN
          !         IB=1 4 13 16 26 27 38 39 (z-axis) OR
          !         IB=1 2 17 18 28 27 43 44 (x-axis) OR
          !         IB=1 3 21 23 26 28 46 48 (y-axis)
          indpg = & !<4>2m(d2d)
               16
       ELSE IF (nc==16 .AND. (ib(16)==40 .OR. ib(16)==44 .OR. ib(16)==48)) &
            THEN
          !         IB=1 2 3 4 13 14 15 16 25 26 27 28 37 38 39 40 (z-axis) OR
          !         IB=1 2 3 4 17 19 20 18 25 26 27 28 41 43 44 42 (x-axis) OR
          !         IB=1 2 3 4 21 22 24 23 25 26 27 28 45 46 48 47 (y-axis)
          indpg = & !4/mmm(d4h)
               17
       ELSE IF (nc==4 .AND. (ib(4)==4)) THEN
          !         Orthorhombic system
          !         IB=1 2  3  4
          indpg = & !222(d2)
               25
       ELSE IF (nc==4 .AND. (ib(4)==27 .OR. ib(4)==28)) THEN
          !         IB=1 3 26 27 (z-axis) OR
          !         IB=1 2 27 28 (x-axis) OR
          !         IB=1 4 26 28 (y-axis) OR
          indpg = & !mm2(c2v)
               26
       ELSE IF (nc==8) THEN
          !         IB=1 2 3 4 25 26 27 28
          indpg = & !mmm(d2h)
               27
       ELSE IF (nc==12 .AND. (ib(12)==12 .OR. ib(12)==47 .OR. ib(12)==45)) &
            THEN
          !         Cubic system
          !         IB=1 2  3  4  5  6  7  8  9 10 11 12 OR
          !         IB=1 5 11 13 18 23 25 30 35 37 42 47 OR
          !         IB=1 8 10 16 18 21 25 32 34 40 42 45
          indpg = & !23 (t)
               28
       ELSE IF (nc==24 .AND. ib(24)==36) THEN
          !         IB= 1  2  3  4  5  6  7  8  9 10 11 12
          !            25 26 27 28 29 30 31 32 33 34 35 36
          indpg = & !m3 (th)
               29
       ELSE IF (nc==24 .AND. ib(24)==24) THEN
          !         IB=1 2 3 4 5 6 7 8 9 10 11 12
          !           13 14 15 16 17 18 19 20 21 22 23 24
          indpg = & !432 (o)
               30
       ELSE IF (nc==24 .AND. ib(24)==48) THEN
          !         IB=1 2 3 4 5 6 7 8 9 10 11 12
          !           37 38 39 40 41 42 43 45 46 47 48
          indpg = & !<4>3m(td)
               31
       ELSE IF (nc==48) THEN
          !         IB=1..48
          indpg = & !m3m(oh)
               32
       ELSE
          !          WRITE(*,'(" ATFTM1! IHG=",A," NC=",I2)') ICST(IHG),NC
          !          WRITE(*,'(" ATFTM1!",19I3)') (IB(I),I=1,NC)
          !          WRITE(*,'(" ATFTM1! THIS CASE IS UNKNOWN IN THE DATABASE")')
          !         Probably a sub-group of 32
          indpg = -32
       END IF
    ELSE IF (ihg>=6) THEN
       IF (nc==0) THEN
          WRITE (*,'(" ATFTM1! IHG=",A," NC=",I2)') icst(ihg), nc
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          !       Triclinic system
       ELSE IF (nc==1) THEN
          !         IB=1
          indpg = & !1 (c1)
               1
       ELSE IF (nc==2 .AND. ib(2)==13) THEN
          !         IB=1 13
          indpg = & !<1>(ci)
               2
       ELSE IF (nc==2 .AND. (ib(2)==4)) & !         Monoclinic system
            THEN
          !         IB=1 4
          !2[001]
          indpg = & !2 (c2)
               3
       ELSE IF (nc==2 .AND. (ib(2)==16)) THEN
          !         IB=1 16
          indpg = & !m (c1h)
               4
       ELSE IF (nc==4 .AND. (ib(4)==24 .OR. ib(4)==20)) THEN
          !         IB=1 12 13 24 OR
          !         IB=1  8 13 20
          indpg = & !2/m(c2h)
               5
       ELSE IF (nc==3 .AND. ib(3)==5) THEN
          !         Trigonal system
          !         IB=1 3 5
          indpg = & !3 (c3)
               6
       ELSE IF (nc==6 .AND. ib(6)==17) THEN
          !         IB=1 13 15 17 35
          indpg = & !<3>(c3i)
               7
       ELSE IF (nc==6 .AND. ib(6)==11) THEN
          !         IB=1 7 9 11 35
          indpg = & !32 (d3)
               8
       ELSE IF (nc==6 .AND. ib(6)==23) THEN
          !         IB=1 3 5 19 21 23
          indpg = & !3m (c3v)
               9
       ELSE IF (nc==12 .AND. ib(12)==23) THEN
          !         IB=1 3 5 7 9 11 13 15 17 19 21 23
          indpg = & !<3>m(d3d)
               10
       ELSE IF (nc==6 .AND. ib(6)==6) THEN
          !         Hexagonal system
          !         IB=1 2 3 4 5 6
          indpg = & !6 (c6)
               18
       ELSE IF (nc==6 .AND. ib(6)==18) THEN
          !         IB=1 3 5 14 16 18
          indpg = & !<6>(c3h)
               19
       ELSE IF (nc==12 .AND. ib(12)==18) THEN
          !         IB=1 2 3 4 5 6 13 14 15 16 17 18
          indpg = & !6/m(c6h)
               20
       ELSE IF (nc==12 .AND. ib(12)==12) THEN
          !         IB=1 2 3 4 5 6 7 8 9 10 11 12
          indpg = & !622(d6)
               21
       ELSE IF (nc==12 .AND. ib(2)==2 .AND. ib(12)==24) THEN
          !         IB=1 2 3 4 5 6 19 20 21 22 23 24
          indpg = & !6mm(c6v)
               22
       ELSE IF (nc==12 .AND. ib(2)==3 .AND. ib(12)==24) THEN
          !         IB=1 3 5 7 9 11 14 16 18 20 22 24
          indpg = & !<6>m2(d3h)
               23
       ELSE IF (nc==24) THEN
          !         IB=1..24
          indpg = & !6/mmm(d6h)
               24
       ELSE
          !         Probably a sub-group of 24
          !          WRITE(*,'(" ATFTM1! IHG=",A," NC=",I2)') ICST(IHG),NC
          !          WRITE(*,'(" ATFTM1!",48I3)') (IB(I),I=1,NC)
          !          WRITE(*,'(" ATFTM1! THIS CASE IS UNKNOWN IN THE DATABASE")')
          indpg = -24
       END IF
    END IF
    !------------------------------------------------------------------------------!
    !        Determination if the space group is symmorphic or not
    !------------------------------------------------------------------------------!
    IF (isy/=1) THEN
       !       Transform V in cartesian coordinates
       DO n = 1, nc
          DO i = 1, 3
             vc(i,n) = a(i,1)*v(1,n) + a(i,2)*v(2,n) + a(i,3)*v(3,n)
          END DO
       END DO
       CALL symmorphic(nc,ib,r,vc,ai,info,origin,delta,error)
       IF (info==1) THEN
          CALL rlv3(ai,origin,xb,il,delta,error)
          !         !!!RLV3 determines -XB in crystal coordinates
          !         !!We want between 0.0 and 1.0
          DO i = 1, 3
             IF (-xb(i)>=0.0_dp) THEN
                origin(i) = -xb(i)
             ELSE
                origin(i) = 1.0_dp - xb(i)
             END IF
          END DO
          DO i = 1, 3
             xb(i) = a(i,1)*origin(1) + a(i,2)*origin(2) + a(i,3)*origin(3)
          END DO
          isy = -1
       ELSE IF (info==0) THEN
          isy = 0
       ELSE
          isy = -2
       END IF
    ELSE
       DO i = 1, 3
          origin(i) = 0.0_dp
       END DO
    END IF
  END SUBROUTINE atftm1
! *****************************************************************************
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
!> \param error ...
!> \par Original Header
!> ------------------------------------------------------------------------------!
!>         WRITTEN IN MAY 14TH, 1998 (T.D.)
!>         CHECK IF RX+VR GIVES THE SAME LATTICE AS X
!>         BUILD THE ATOM TRANSFORMATION TABLE
!> ------------------------------------------------------------------------------!
!>         INPUT:
!>           N   ROTATION NUMBER (INDEX USED IN F0 BETWEEN 1 AND 48)
!>           NAT           NUMBER OF ATOMS
!>           TY(1:NAT)     TYPE OF ATOMS
!>           RX(1:3,1:NAT) ATOMIC COORDINATES FROM Nth ROTATION (CART.)
!>           X(1:3,1:NAT)  ATOMIC COORDINATES (CARTESIAN)
!>           VR(1:3)       TRANSLATION VECTOR (CRYSTAL COOR.)
!>           AI(1:3,1:3)   LATTICE RECIPROCAL VECTORS
!>           NODUPLI       .TRUE., THE CELL IS A PRIMITIVE ONE
!>                         WE CAN SPEED UP
!>         DELTA           REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>         OUTPUT:
!>           F0(1:49,1:NAT) ATOM TRANSFORMATION TABLE
!>              F0 IS THE FUNCTION DEFINED IN MARADUDIN AND VOSK0
!>              BY EQ.(2.35).
!>              IT DEFINES THE ATOM TRANSFORMATION TABLE
!>           OKSYM          TRUE IF RX+VR = X
!>           ISC(1:NAT)     SCRATCH ARRAY
!>                          USED TO SPEED UP THE ROUTINE
!>                          EACH ATOM IS ONLY ONCE AN IMAGE
!>                          IF NO DUPLICATION OF THE CELL
! *****************************************************************************
  SUBROUTINE checkrlv3(n,nat,ty,rx,x,vr,f0,ai,isc,nodupli,oksym,delta,error)

    INTEGER                                  :: n, nat, ty(nat)
    REAL(KIND=dp)                            :: rx(3,nat), x(3,nat), vr(3)
    INTEGER                                  :: f0(49,nat)
    REAL(KIND=dp)                            :: ai(3,3)
    INTEGER                                  :: isc(nat)
    LOGICAL                                  :: nodupli, oksym
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'checkrlv3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, ib, il
    REAL(KIND=dp)                            :: vt(3), xb(3)

    DO ia = 1, nat
       isc(ia) = 0
    END DO
    !     Now we check if ROT(N)+VR gives a correct symmetry.
    DO ia = 1, nat
       DO ib = 1, nat
          IF (ty(ia)==ty(ib) .AND. isc(ib)==0) THEN
             xb(1) = rx(1,ia) - x(1,ib)
             xb(2) = rx(2,ia) - x(2,ib)
             xb(3) = rx(3,ia) - x(3,ib)
             CALL rlv3(ai,xb,vt,il,delta,error)
             !           VT STANDS FOR V-TEST
             oksym = (ABS((vr(1)-vt(1))-NINT(vr(1)-vt(1)))<delta) .AND. &
                  (ABS((vr(2)-vt(2))-NINT(vr(2)-vt(2)))<delta) .AND. &
                  (ABS((vr(3)-vt(3))-NINT(vr(3)-vt(3)))<delta)
             IF (oksym) THEN
                IF (nodupli) isc(ib) = 1
                f0(n,ia) = ib
                !             IR+VR is the good one: another symmetry operation
                !             Next atom
                GO TO 100
             END IF
          END IF
       END DO
       !       VR is not the correct translation vector
       RETURN
100    CONTINUE
    END DO
  END SUBROUTINE checkrlv3
! *****************************************************************************
!> \brief ...
!> \param nc ...
!> \param ib ...
!> \param r ...
!> \param v ...
!> \param ai ...
!> \param info ...
!> \param origin ...
!> \param delta ...
!> \param error ...
!> \par Original Header
!> ------------------------------------------------------------------------------!
!>         Check if the group is symmorphic with a non-standard origin
!>         WARNING: If there are equivalent atoms, this routine could
!>         not determine if the space group is symmorphic
!>         So you have to check if the solution V=0 works (see ATFTM1)
!> ------------------------------------------------------------------------------!
!>         INPUT:
!>           NC Number of operations
!>           IB(NC) Index of operation in R
!>           R(3,3,48) Rotations
!>           V(3,NC) Fractional translations related to R(3,3,IB(NC))
!>                   R AND V ARE IN CARTESIAN COORDINATES
!>           AI(I,J) ARE THE RECIPROCAL LATTICE VECTORS,
!>                   B(I) = AI(I,J),J=1,2,3
!>         DELTA     REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>
!>         OUTPUT:
!>           ORIGIN(1:3) Give standard origin (cartesian coordinates)
!>                    Give the standard origin with smallest coordinates
!>                    if NTVEC /= 1
!>           INFO = 1 The group is symmorphic
!>           INFO = 0 The group is not symmorphic
!>           INFO =-1 The routine cannot determine
! *****************************************************************************
  SUBROUTINE symmorphic(nc,ib,r,v,ai,info,origin,delta,error)
    INTEGER                                  :: nc, ib(nc)
    REAL(KIND=dp)                            :: r(3,3,48), v(3,nc), ai(3,3)
    INTEGER                                  :: info
    REAL(KIND=dp)                            :: origin(3), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'symmorphic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, ierror, igood(3), il, &
                                                imissing2, imissing3, iok(3), &
                                                ionly, ir, j, j1
    REAL(KIND=dp)                            :: diag, dif, r2(2,2), r3(3,3), &
                                                vr(3), xb(3)

    DO i = 1, 3
       iok(i) = 0
    END DO
    DO i = 1, 3
       origin(i) = 0.0_dp
    END DO
    !     Find a point A / V_R = (1-R).OA
    DO ir = 1, nc
       dif = v(1,ir)*v(1,ir) + v(2,ir)*v(2,ir) + v(3,ir)*v(3,ir)
       IF (dif>delta*delta) THEN
          DO i = 1, 3
             igood(i) = 1
          END DO
          !         V is non-zero. Construct matrix 1-R
          DO i = 1, 3
             DO j = 1, 3
                r3(i,j) = -r(i,j,ib(ir))
             END DO
             r3(i,i) = 1 + r3(i,i)
          END DO
          CALL invmat(r3(1:3,1:3),ierror,error)
          IF (ierror==0) THEN
             !           The matrix 3x3 has an inverse.
             DO i = 1, 3
                vr(i) = r3(i,1)*v(1,ir) + r3(i,2)*v(2,ir) + r3(i,3)*v(3,ir)
             END DO
          ELSE
             !           IERROR gives the column which causes some trouble
             !           Construct matrix 1-R with 2x2
             igood(ierror) = 0
             imissing3 = ierror
             i1 = 0
             DO i = 1, 3
                IF (i/=ierror) THEN
                   i1 = i1 + 1
                   j1 = 0
                   DO j = 1, 3
                      IF (j/=ierror) THEN
                         j1 = j1 + 1
                         r2(i1,j1) = -r(i,j,ib(ir))
                      END IF
                   END DO
                   r2(i1,i1) = 1 + r2(i1,i1)
                END IF
             END DO
             CALL invmat(r2(1:2,1:2),ierror,error)
             IF (ierror==0) THEN
                !             The matrix 2X2 has an inverse.
                !             Solve Vxy = (1-R).OAxy + OAz R3z (z is IMISSING3)
                i1 = 0
                DO i = 1, 3
                   IF (igood(i)==1) THEN
                      i1 = i1 + 1
                      vr(i) = 0.0_dp
                      j1 = 0
                      DO j = 1, 3
                         IF (igood(j)==1) THEN
                            j1 = j1 + 1
                            vr(i) = vr(i) + r2(i1,j1)*(v(j,ir)+origin(imissing3)*r &
                                 (j,imissing3,ib(ir)))
                         END IF
                      END DO
                   ELSE
                      vr(i) = origin(i)
                   END IF
                END DO
             ELSE
                !             Construct matrix 1-R with 1x1
                i1 = 0
                DO i = 1, 3
                   IF (i/=imissing3) THEN
                      i1 = i1 + 1
                      IF (i1==ierror) THEN
                         igood(i) = 0
                         imissing2 = i
                      ELSE
                         ionly = i
                      END IF
                   END IF
                END DO
                diag = (1-r(ionly,ionly,ib(ir)))
                IF (ABS(diag)>delta) THEN
                   vr(ionly) = 1.0_dp/diag*(v(ionly,ir)+origin(imissing3)*r( &
                        ionly,imissing3,ib(ir))+origin(imissing2)*r(ionly, &
                        imissing2,ib(ir)))
                ELSE
                   vr(ionly) = origin(ionly)
                   igood(ionly) = 0
                END IF
                vr(imissing3) = origin(imissing3)
                vr(imissing2) = origin(imissing2)
             END IF
          END IF
          !------------------------------------------------------------------------------!
          !         Compare VR with ORIGIN
          dif = 0.0_dp
          !         If NTVEC /=1 there are NTVEC possible standard origins
          DO i = 1, 3
             IF (iok(i)==1) THEN
                dif = dif + ABS(origin(i)-vr(i))
             END IF
          END DO
          IF (dif>delta) THEN
             !           Non-symmorphic
             info = 0
             RETURN
          ELSE
             DO i = 1, 3
                IF (iok(i)/=1 .AND. igood(i)==1) THEN
                   iok(i) = 1
                   origin(i) = vr(i)
                END IF
             END DO
          END IF
       END IF
    END DO
    !------------------------------------------------------------------------------!
    IF (iok(1)==0 .AND. iok(2)==0 .AND. iok(3)==0) THEN
       !       Cannot not determine
       info = -1
       RETURN
    END IF
    !     The group is symmorphic
    info = 1
    !     Check
    DO ir = 1, nc
       DO i = 1, 3
          vr(i) = r(i,1,ib(ir))*origin(1) + r(i,2,ib(ir))*origin(2) + &
               r(i,3,ib(ir))*origin(3)
          vr(i) = (origin(i)-vr(i)) - v(i,ir)
       END DO
       CALL rlv3(ai,vr,xb,il,delta,error)
       dif = ABS(xb(1)) + ABS(xb(2)) + ABS(xb(3))
       IF (dif>delta) THEN
          !         Non-symmorphic
          info = 0
          RETURN
       END IF
    END DO
  END SUBROUTINE symmorphic

! *****************************************************************************
!> \brief ...
!> \param ihc ...
!> \param r ...
!> \param error ...
!> \par Original Header
!> ------------------------------------------------------------------------------!
!>         WRITTEN ON FEBRUARY 17TH, 1976
!>         GENERATION OF THE X,Y,Z-TRANSFORMATION MATRICES 3X3
!>         FOR HEXAGONAL AND CUBIC GROUPS
!>         SUBROUTINES NEEDED -- NONE
!> ------------------------------------------------------------------------------!
!>         THIS IS IDENTICAL WITH THE SUBROUTINE ROT OF WORLTON-WARREN
!>         (IN THE AC-COMPLEX), ONLY THE WAY OF TRANSFERRING THE DATA
!>         WAS CHANGED
!> ------------------------------------------------------------------------------!
!>         INPUT DATA:
!>         IHC SWITCH DETERMINING IF WE DESIRE
!>             THE HEXAGONAL GROUP(IHC=0) OR THE CUBIC GROUP (IHC=1)
!>         OUTPUT DATA:
!>         R...THE 3X3 MATRICES OF THE DESIRED COORDINATE REPRESENTATION
!>             THEIR NUMBERING CORRESPONDS TO THE SYMMETRY ELEMENTS AS
!>             LISTE IN WORLTON-WARREN
!>                      (COMPUT. PHYS. COMM. 3(1972) 88--117)
!>         FOR IHC=0 THE FIRST 24 MATRICES OF THE ARRAY R REPRESENT
!>                   THE FULL HEXAGONAL GROUP D(6H)
!>         FOR IHC=1 THE FIRST 48 MATRICES OF THE ARRAY R REPRESENT
!>                   THE FULL CUBIC GROUP O(H)
! *****************************************************************************
  SUBROUTINE rot1(ihc,r,error)
    INTEGER                                  :: ihc
    REAL(KIND=dp)                            :: r(3,3,48)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rot1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, n, nv
    REAL(KIND=dp)                            :: c, s

!     Variables
!------------------------------------------------------------------------------!

    DO j = 1, 3
       DO i = 1, 3
          DO n = 1, 48
             r(i,j,n) = 0.0_dp
          END DO
       END DO
    END DO
    IF (ihc==0) THEN
       !------------------------------------------------------------------------------!
       !       DEFINE THE GENERATORS FOR THE ROTATION MATRICES--HEXAGONAL GROUP
       !------------------------------------------------------------------------------!
       c = 0.5E0_dp
       s = 0.5E0_dp*SQRT(3.0E0_dp)
       r(1,1,2) = c
       r(1,2,2) = -s
       r(2,1,2) = s
       r(2,2,2) = c
       r(1,1,7) = -c
       r(1,2,7) = -s
       r(2,1,7) = -s
       r(2,2,7) = c
       DO n = 1, 6
          r(3,3,n) = 1.0_dp
          r(3,3,n+18) = 1.0_dp
          r(3,3,n+6) = -1.0_dp
          r(3,3,n+12) = -1.0_dp
       END DO
       !------------------------------------------------------------------------------!
       !          GENERATE THE REST OF THE ROTATION MATRICES
       !------------------------------------------------------------------------------!
       DO i = 1, 2
          r(i,i,1) = 1.0_dp
          DO j = 1, 2
             r(i,j,6) = r(j,i,2)
             DO k = 1, 2
                r(i,j,3) = r(i,j,3) + r(i,k,2)*r(k,j,2)
                r(i,j,8) = r(i,j,8) + r(i,k,2)*r(k,j,7)
                r(i,j,12) = r(i,j,12) + r(i,k,7)*r(k,j,2)
             END DO
          END DO
       END DO
       DO i = 1, 2
          DO j = 1, 2
             r(i,j,5) = r(j,i,3)
             DO k = 1, 2
                r(i,j,4) = r(i,j,4) + r(i,k,2)*r(k,j,3)
                r(i,j,9) = r(i,j,9) + r(i,k,2)*r(k,j,8)
                r(i,j,10) = r(i,j,10) + r(i,k,12)*r(k,j,3)
                r(i,j,11) = r(i,j,11) + r(i,k,12)*r(k,j,2)
             END DO
          END DO
       END DO
       DO n = 1, 12
          nv = n + 12
          DO i = 1, 2
             DO j = 1, 2
                r(i,j,nv) = -r(i,j,n)
             END DO
          END DO
       END DO
    ELSE
       !------------------------------------------------------------------------------!
       !          DEFINE THE GENERATORS FOR THE ROTATION MATRICES-CUBIC GROUP
       !------------------------------------------------------------------------------!
       r(1,3,9) = 1.0_dp
       r(2,1,9) = 1.0_dp
       r(3,2,9) = 1.0_dp
       r(1,1,19) = 1.0_dp
       r(2,3,19) = -1.0_dp
       r(3,2,19) = 1.0_dp
       DO i = 1, 3
          r(i,i,1) = 1.0_dp
          DO j = 1, 3
             r(i,j,20) = r(j,i,19)
             r(i,j,5) = r(j,i,9)
             DO k = 1, 3
                r(i,j,2) = r(i,j,2) + r(i,k,19)*r(k,j,19)
                r(i,j,16) = r(i,j,16) + r(i,k,9)*r(k,j,19)
                r(i,j,23) = r(i,j,23) + r(i,k,19)*r(k,j,9)
             END DO
          END DO
       END DO
       DO i = 1, 3
          DO j = 1, 3
             DO k = 1, 3
                r(i,j,6) = r(i,j,6) + r(i,k,2)*r(k,j,5)
                r(i,j,7) = r(i,j,7) + r(i,k,16)*r(k,j,23)
                r(i,j,8) = r(i,j,8) + r(i,k,5)*r(k,j,2)
                r(i,j,10) = r(i,j,10) + r(i,k,2)*r(k,j,9)
                r(i,j,11) = r(i,j,11) + r(i,k,9)*r(k,j,2)
                r(i,j,12) = r(i,j,12) + r(i,k,23)*r(k,j,16)
                r(i,j,14) = r(i,j,14) + r(i,k,16)*r(k,j,2)
                r(i,j,15) = r(i,j,15) + r(i,k,2)*r(k,j,16)
                r(i,j,22) = r(i,j,22) + r(i,k,23)*r(k,j,2)
                r(i,j,24) = r(i,j,24) + r(i,k,2)*r(k,j,23)
             END DO
          END DO
       END DO
       DO i = 1, 3
          DO j = 1, 3
             DO k = 1, 3
                r(i,j,3) = r(i,j,3) + r(i,k,5)*r(k,j,12)
                r(i,j,4) = r(i,j,4) + r(i,k,5)*r(k,j,10)
                r(i,j,13) = r(i,j,13) + r(i,k,23)*r(k,j,11)
                r(i,j,17) = r(i,j,17) + r(i,k,16)*r(k,j,12)
                r(i,j,18) = r(i,j,18) + r(i,k,16)*r(k,j,10)
                r(i,j,21) = r(i,j,21) + r(i,k,12)*r(k,j,15)
             END DO
          END DO
       END DO
       DO n = 1, 24
          nv = n + 24
          DO i = 1, 3
             DO j = 1, 3
                r(i,j,nv) = -r(i,j,n)
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE rot1
! *****************************************************************************
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
!> \param error ...
!> \par Original Header
!> ------------------------------------------------------------------------------!
!>         WRITTEN ON SEPTEMBER 12-20TH, 1979 BY K.K.
!>         MODIFIED 26-MAY-82 BY OLE HOLM NIELSEN
!>         GENERATION OF SPECIAL POINTS FOR AN ARBITRARY LATTICE,
!>         FOLLOWING THE METHOD MONKHORST,PACK,
!>         PHYS. REV. B13 (1976) 5188
!>         MODIFIED BY MACDONALD, PHYS. REV. B18 (1978) 5897
!>         THE SUBROUTINE IS WRITTEN ASSUMING THAT THE POINTS ARE
!>         GENERATED IN THE RECIPROCAL SPACE.
!>         IF, HOWEVER, THE B1,B2,B3 ARE REPLACED BY A1,A2,A3, THEN
!>         SPECIAL POINTS IN THE DIRECT SPACE CAN BE PRODUCED, AS WELL.
!>         (NO MULTIPLICATION BY 2PI IS THEN NECESSARY.)
!>         IN THE CASE OF NONSYMMORPHIC GROUPS, THE APPLICATION IN THE
!>         DIRECT SPACE WOULD PROBABLY REQUIRE A CERTAIN CAUTION.
!>         SUBROUTINES NEEDED: BZDEFI,BZRDUC,INBZ,MESH
!>         IN THE CASES WHERE THE POINT GROUP OF THE CRYSTAL DOES NOT
!>         CONTAIN INVERSION. THE LATTER MAY BE ADDED IF WE WISH
!>         (SEE COMMENT TO THE SWITCH INV).
!>         REDUCTION TO THE 1ST BRILLOUIN ZONE IS DONE
!>         BY ADDING G-VECTORS TO FIND THE SHORTEST WAVE-VECTOR.
!>         THE ROTATIONS OF THE BRAVAIS LATTICE ARE APPLIED TO THE
!>         MONKHORST/PACK MESH IN ORDER TO FIND ALL K-POINTS
!>         THAT ARE RELATED BY SYMMETRY. (OLE HOLM NIELSEN)
!> ------------------------------------------------------------------------------!
!>         INPUT DATA:
!>         IOUT:    LOGICAL UNIT FOR OUTPUT
!>                  IF (IOUT.LE.0) NO MESSAGE
!>         IQ1,IQ2,IQ3 .. PARAMETER Q OF MONKHORST AND PACK,
!>                  GENERALIZED AND DIFFERENT FOR THE 3 DIRECTIONS B1,
!>                  B2 AND B3
!>         WVK0 ... THE 'ARBITRARY' SHIFT OF THE WHOLE MESH, DENOTED K0
!>                  IN MACDONALD. WVK0 = 0 CORRESPONDS TO THE ORIGINAL
!>                  SCHEME OF MONKHORST AND PACK.
!>                  UNITS: 2PI/(UNITS OF LENGTH  USED IN A1, A2, A3),
!>                  I.E. THE SAME  UNITS AS THE GENERATED SPECIAL POINTS
!>         NKPOINT .. VARIABLE DIMENSION OF THE (OUTPUT) ARRAYS WVKL,
!>                  LWGHT,LROT, I.E. SPACE RESERVED FOR THE SPECIAL
!>                  POINTS AND ACCESSORIES.
!>                  NKPOINT HAS TO BE .GE. NTOT (TOTAL NUMBER OF SPECIAL
!>                  POINTS. THIS IS CHECKED BY THE SUBROUTINE.
!>         ISTRIZ . INDICATES WHETHER ADDITIONAL MESH POINTS SHOULD BE
!>                  GENERATED BY APPLYING GROUP OPERATIONS TO THE MESH.
!>                  ISTRIZ=+1 MEANS SYMMETRIZE
!>                  ISTRIZ=-1 MEANS DO NOT SYMMETRIZE
!>         THE FOLLOWING INPUT DATA MAY BE OBTAINED FROM THE SBRT.
!>         B1,B2,B3 .. RECIPROCAL LATTICE VECTORS, NOT MULTIPLIED BY
!>                  GROUP1: ANY 2PI (IN UNITS RECIPROCAL TO THOSE
!>                          OF A1,A2,A3)
!>         INV .... CODE INDICATING WHETHER WE WISH TO ADD THE INVERSION
!>                  TO THE POINT GROUP OF THE CRYSTAL OR NOT (IN THE
!>                  CASE THAT THE POINT GROUP DOES NOT CONTAIN ANY).
!>                  INV=0 MEANS: DO NOT ADD INVERSION
!>                  INV.NE.0 MEANS: ADD THE INVERSION
!>                  INV.NE.0 SHOULD BE THE STANDARD CHOICE WHEN SPPT2
!>                  IS USED IN RECIPROCAL SPACE - IN ORDER TO MAKE
!>                  USE OF THE HERMITICITY OF HAMILTONIAN.
!>                  WHEN USED IN DIRECT SPACE, THE RIGHT CHOICE OF INV
!>                  WILL DEPEND ON THE NATURE OF THE PHYSICAL PROBLEM.
!>                  IN THE CASES WHERE THE INVERSION IS ADDED BY THE
!>                  SWITCH INV, THE LIST IB WILL NOT BE MODIFIED BUT IN
!>                  THE OUTPUT LIST LROT SOME OF THE OPERATIONS WILL
!>                  APPEAR WITH NEGATIVE SIGN; THIS MEANS THAT THEY HAVE
!>                  TO BE APPLIED MULTIPLIED BY INVERSION.
!>         NC ..... TOTAL NUMBER OF ELEMENTS IN THE POINT GROUP OF THE
!>                  CRYSTAL
!>         IB ..... LIST OF THE ROTATIONS CONSTITUTING THE POINT GROUP
!>                  OF THE CRYSTAL. THE NUMBERING IS THAT DEFINED IN
!>                  WORLTON AND WARREN, I.E. THE ONE MATERIALIZED IN THE
!>                  ARRAY R (SEE BELOW)
!>                  ONLY THE FIRST NC ELEMENTS OF THE ARRAY IB ARE
!>                  MEANINGFUL
!>         R ...... LIST OF THE 3 X 3 ROTATION MATRICES
!>                  (XYZ REPRESENTATION OF THE O(H) OR D(6)H GROUPS)
!>                  ALL 48 OR 24 MATRICES ARE LISTED.
!>         NCBRAV . TOTAL NUMBER OF ELEMENTS IN RBRAV
!>         IBRAV .. LIST OF NCBRAV OPERATIONS OF THE BRAVAIS LATTICE
!>         DELTA    REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!> ------------------------------------------------------------------------------!
!>         OUTPUT DATA:
!>         NTOT ... TOTAL NUMBER OF SPECIAL POINTS
!>                  IF NTOT APPEARS NEGATIVE, THIS IS AN ERROR SIGNAL
!>                  WHICH MEANS THAT THE DIMENSION NKPOINT WAS CHOSEN
!>                  TOO SMALL SO THAT THE ARRAYS WVKL ETC. CANNOT
!>                  ACCOMODATE ALL THE GENERATED SPECIAL POINTS.
!>                  IN THIS CASE THE ARRAYS WILL BE FILLED UP TO NKPOINT
!>                  AND FURTHER GENERATION OF NEW POINTS WILL BE
!>                  INTERRUPTED.
!>         WVKL ... LIST OF SPECIAL POINTS.
!>                  CARTESIAN COORDINATES AND NOT MULTIPLIED BY 2*PI.
!>                  ONLY THE FIRST NTOT VECTORS ARE MEANINGFUL
!>                  ALTHOUGH NO 2 POINTS FROM THE LIST ARE EQUIVALENT
!>                  BY SYMMETRY, THIS SUBROUTINE STILL HAS A KIND OF
!>                  'BEAUTY DEFECT': THE POINTS FINALLY
!>                  SELECTED ARE NOT NECESSARILY SITUATED IN A
!>                  'COMPACT' IRREDUCIBLE BRILL.ZONE; THEY MIGHT LIE IN
!>                  DIFFERENT IRREDUCIBLE PARTS OF THE B.Z. - BUT THEY
!>                  DO REPRESENT AN IRREDUCIBLE SET FOR INTEGRATION
!>                  OVER THE ENTIRE B.Z.
!>         LWGHT ... THE LIST OF WEIGHTS OF THE CORRESPONDING POINTS.
!>                  THESE WEIGHTS ARE NOT NORMALIZED (JUST INTEGERS)
!>         LROT ... FOR EACH SPECIAL POINT THE 'UNFOLDING ROTATIONS'
!>                  ARE LISTED. IF E.G. THE WEIGHT OF THE I-TH SPECIAL
!>                  POINT IS LWGHT(I), THEN THE ROTATIONS WITH NUMBERS
!>                  LROT(J,I), J=1,2,...,LWGHT(I) WILL 'SPREAD' THIS
!>                  SINGLE POINT FROM THE IRREDUCIBLE PART OF B.Z. INTO
!>                  SEVERAL POINTS IN AN ELEMENTARY UNIT CELL
!>                  (PARALLELOPIPED) OF THE RECIPROCAL SPACE.
!>                  SOME OPERATION NUMBERS IN THE LIST LROT MAY APPEAR
!>                  NEGATIVE, THIS MEANS THAT THE CORRESPONDING ROTATION
!>                  HAS TO BE APPLIED WITH INVERSION (THE LATTER HAVING
!>                  BEEN ARTIFICIALLY ADDED AS SYMMETRY OPERATION IN
!>                  CASE INV.NE.0).NO OTHER EFFORT WAS TAKEN,TO RENUMBER
!>                  THE ROTATIONS WITH MINUS SIGN OR TO EXTEND THE
!>                  LIST OF THE POINT-GROUP OPERATIONS IN THE LIST NB.
!>         INCLUD ... INTEGER ARRAY USED BY SPPT2 INCLUD(NKPOINT)
!>                  THE FIRST BIT (0) IS USED BY THE ROUTINE.
!>                  THE OTHER BITS GIVE THE K-POINT INDEX IN
!>                  THE SPECIAL K-POINT TABLE.
!> ------------------------------------------------------------------------------!
!>         NHASH    USED BY MESH ROUTINE
!>         LIST     INTEGER ARRAY USED BY MESH  LIST(NHASH+NKPOINT)
!>         RLIST    REAL(KIND=dp)  ARRAY USED BY MESH  RLIST(3,NKPOINT)
!> ------------------------------------------------------------------------------!
!>         Use bit manipulations functions
!>          IBSET(I,POS) sets the bit POS to 1 in I integer
!>          IBCLR(I,POS) clears the bit POS to 1 in I integer
!>          BTEST(I,POS) .TRUE. if bit POS is 1 in I integer
! *****************************************************************************
  SUBROUTINE sppt2(iout,iq1,iq2,iq3,wvk0,nkpoint,a1,a2,a3,b1,b2,b3,inv,nc, &
       ib,r,ntot,wvkl,lwght,lrot,ncbrav,ibrav,istriz,nhash,includ,list, &
       rlist,delta,error)
    INTEGER                                  :: iout, iq1, iq2, iq3
    REAL(KIND=dp)                            :: wvk0(3)
    INTEGER                                  :: nkpoint
    REAL(KIND=dp)                            :: a1(3), a2(3), a3(3), b1(3), &
                                                b2(3), b3(3)
    INTEGER                                  :: inv, nc, ib(48)
    REAL(KIND=dp)                            :: r(3,3,48)
    INTEGER                                  :: ntot
    REAL(KIND=dp)                            :: wvkl(3,nkpoint)
    INTEGER :: lwght(nkpoint), lrot(48,nkpoint), ncbrav, ibrav(48), istriz, &
      nhash, includ(nkpoint), list(nkpoint+nhash)
    REAL(KIND=dp)                            :: rlist(3,nkpoint), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sppt2', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: no = 0, nrsdir = 100, yes = 1

    INTEGER                                  :: i, i1, i2, i3, ibsign, &
                                                igarb0, igarbage, igarbg, &
                                                imesh, iop, iplace, iremov, &
                                                iwvk, j, jplace, k, n, nplane
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: diff, proja(3), projb(3), &
                                                rsdir(4,nrsdir), ur1, ur2, &
                                                ur3, wva(3), wvk(3)

    DATA iplace/ -2/
    !------------------------------------------------------------------------------!
    failure = .FALSE.
    !------------------------------------------------------------------------------!

    ntot = 0
    DO i = 1, nkpoint
       lrot(1,i) = 1
       DO j = 2, 48
          lrot(j,i) = 0
       END DO
    END DO
    DO i = 1, nkpoint
       includ(i) = no
    END DO
    DO i = 1, 3
       wva(i) = 0.0_dp
    END DO
    !------------------------------------------------------------------------------!
    !        DEFINE THE 1ST BRILLOUIN ZONE
    !------------------------------------------------------------------------------!
    CALL bzdefi(iout,b1,b2,b3,rsdir,nrsdir,nplane,delta,error)
    !------------------------------------------------------------------------------!
    !        Generation of the mesh (they are not multiplied by 2*pi) by
    !        the Monkhorst/Pack algorithm, supplemented by all rotations
    !------------------------------------------------------------------------------!
    !     Initialize the list of vectors
    iplace = -2
    CALL mesh(iout,wva,iplace,igarb0,igarbg,nkpoint,nhash,list,rlist, &
         delta,error)
    imesh = 0
    DO i1 = 1, iq1
       DO i2 = 1, iq2
          DO i3 = 1, iq3
             ur1 = REAL(1+iq1-2*i1,KIND=dp)/REAL(2*iq1,KIND=dp)
             ur2 = REAL(1+iq2-2*i2,KIND=dp)/REAL(2*iq2,KIND=dp)
             ur3 = REAL(1+iq3-2*i3,KIND=dp)/REAL(2*iq3,KIND=dp)
             DO i = 1, 3
                wvk(i) = ur1*b1(i) + ur2*b2(i) + ur3*b3(i) + wvk0(i)
             END DO
             !           Reduce WVK to the 1st Brillouin zone
             CALL bzrduc(wvk,a1,a2,a3,b1,b2,b3,rsdir,nrsdir,nplane,delta,error)
             IF (istriz==1) THEN
                !             Symmetrization of the k-points mesh.
                !             Apply all the Bravais lattice operations to WVK
                DO iop = 1, ncbrav
                   DO i = 1, 3
                      wva(i) = 0.0_dp
                      DO j = 1, 3
                         wva(i) = wva(i) + r(i,j,ibrav(iop))*wvk(j)
                      END DO
                   END DO
                   !               Check that WVA is inside the 1 Bz.
                   IF (inbz(wva,rsdir,nrsdir,nplane,delta)==no) GO TO 450
                   !               Place WVA in list
                   iplace = 0
                   CALL mesh(iout,wva,iplace,igarb0,igarbg,nkpoint,nhash,list, &
                        rlist,delta,error)
                   !               If WVA was new (and therefore inserted),
                   !               IPLACE is the number.
                   IF (iplace>0) imesh = iplace
                   IF (iplace>nkpoint) GO TO 470
                END DO
             ELSE
                !             Place WVK in list
                iplace = 0
                CALL mesh(iout,wvk,iplace,igarb0,igarbg,nkpoint,nhash,list, &
                     rlist,delta,error)
                imesh = iplace
                IF (iplace>nkpoint) GO TO 470
             END IF
          END DO
       END DO
    END DO
    IF (iout>0) THEN
       !       IMESH: Number of k points in the mesh.
       WRITE (iout,'(" K290| The wavevector mesh contains ",i5," points")') &
            imesh
       WRITE (iout,'(" K290| The points are:")')
       DO i = 1, imesh
          CALL mesh(iout,wva,i,igarb0,igarbg,nkpoint,nhash,list,rlist,delta,error)
          IF (MOD(i,2)==1) THEN
             WRITE (iout,'(1X,I5,3F10.4)',ADVANCE="NO") i, wva
          ELSE
             WRITE (iout,'(1X,I5,3F10.4)') i, wva
          END IF
       END DO
       WRITE (iout,*)
    END IF
    !------------------------------------------------------------------------------!
    IF (istriz==1) THEN
       !       Now figure out if any special point difference (K - K'') is an
       !       integral multiple of a reciprocal-space vector
       iremov = 0
       DO i = 1, (imesh-1)
          iplace = i
          CALL mesh(iout,wva,iplace,igarb0,igarbg,nkpoint,nhash,list,rlist, &
               delta,error)
          !         Project WVA onto B1,2,3:
          proja(1) = 0.0_dp
          proja(2) = 0.0_dp
          proja(3) = 0.0_dp
          DO k = 1, 3
             proja(1) = proja(1) + wva(k)*a1(k)
             proja(2) = proja(2) + wva(k)*a2(k)
             proja(3) = proja(3) + wva(k)*a3(k)
          END DO
          !         Now loop over all the rest of the mesh points
          DO j = (i+1), imesh
             jplace = j
             CALL mesh(iout,wvk,jplace,igarb0,igarbg,nkpoint,nhash,list, &
                  rlist,delta,error)
             !           Project WVK onto B1,2,3:
             projb(1) = 0.0_dp
             projb(2) = 0.0_dp
             projb(3) = 0.0_dp
             DO k = 1, 3
                projb(1) = projb(1) + wvk(k)*a1(k)
                projb(2) = projb(2) + wvk(k)*a2(k)
                projb(3) = projb(3) + wvk(k)*a3(k)
             END DO
             !           Check (PROJA - PROJB): Is it integral ?
             DO k = 1, 3
                diff = proja(k) - projb(k)
                IF (ABS(REAL(NINT(diff),KIND=dp)-diff)>delta) GO TO 280
             END DO
             !           DIFF is integral: remove WVK from mesh:
             CALL remove(wvk,jplace,igarb0,igarbg,nkpoint,nhash,list,rlist, &
                  delta,error)
             !           If WVK actually removed, increment IREMOV
             IF (jplace>0) iremov = iremov + 1
280          CONTINUE
          END DO
       END DO
       IF (iremov>0 .AND. iout>0) WRITE (iout,'(A,A,/,A,I6,A,/)') &
            ' K290| Some of these mesh points are related by lattice ', &
            'translation vectors', 'K K290| ', iremov, &
            ' of the mesh points removed.'
    END IF
    !------------------------------------------------------------------------------!
    !        IN THE MESH OF WAVEVECTORS, NOW SEARCH FOR EQUIVALENT POINTS:
    !        THE INVERSION (TIME REVERSAL !) MAY BE USED.
    !------------------------------------------------------------------------------!
    DO iwvk = 1, imesh
       !        IF(INCLUD(IWVK) .EQ. YES) GOTO 350
       IF (BTEST(includ(iwvk),0)) GO TO 350
       !       IWVK has not been encountered previously: new special point,
       !       (only if WVK is not a garbage vector, however.)
       !        INCLUD(IWVK) = YES
       includ(iwvk) = IBSET(includ(iwvk),0)
       iplace = iwvk
       CALL mesh(iout,wvk,iplace,igarb0,igarbg,nkpoint,nhash,list,rlist, &
            delta,error)
       !       Find out whether Wvk is in the garbage list
       CALL garbag(wvk,igarbage,igarb0,nkpoint,nhash,list,rlist,delta,error)
       IF (igarbage>0) GO TO 350
       ntot = ntot + 1
       !       Give the index in the special k points table.
       includ(iwvk) = includ(iwvk) + ntot*2
       DO i = 1, 3
          wvkl(i,ntot) = wvk(i)
       END DO
       lwght(ntot) = 1
       !------------------------------------------------------------------------------!
       !       Find all the equivalent points (symmetry given by atoms)
       DO n = 1, nc
          !         Rotate:
          DO i = 1, 3
             wva(i) = 0.0_dp
             DO j = 1, 3
                wva(i) = wva(i) + r(i,j,ib(n))*wvk(j)
             END DO
          END DO
          ibsign = + 1
363       CONTINUE
          !         Find WVA in the list
          iplace = -1
          CALL mesh(iout,wva,iplace,igarb0,igarbg,nkpoint,nhash,list,rlist, &
               delta,error)
          IF (iplace==0) THEN
             IF (istriz==-1) THEN
                !             No symmetrisation -> WVA not in the list
                GO TO 364
             ELSE
                !             I think this case is impossible (NC <= NCBRAV)
                !             Error message
                GO TO 490
             END IF
          END IF
          !         Find out whether WVA is in the garbage list
          CALL garbag(wva,igarbage,igarb0,nkpoint,nhash,list,rlist,delta,error)
          IF (igarbage>0) GO TO 370
          !         Was WVA encountered before ?
          !          IF(INCLUD(IPLACE) .EQ. YES) GOTO 364
          IF (BTEST(includ(iplace),0)) GO TO 364
          !         Increment weight.
          lwght(ntot) = lwght(ntot) + 1
          lrot(lwght(ntot),ntot) = ib(n)*ibsign
          !          INCLUD(IPLACE) = YES
          includ(iplace) = IBSET(includ(iplace),0)
          !         This k-point is an image of a special k-point.
          !         Put the index of the special k-point.
          includ(iplace) = includ(iplace) + ntot*2
364       CONTINUE
          IF (ibsign==-1 .OR. inv==0) GO TO 370
          !         The case where we also apply the inversion to WVA
          !         Repeat the search, but for -WVA
          ibsign = -1
          DO i = 1, 3
             wva(i) = -wva(i)
          END DO
          GO TO 363
370       CONTINUE
       END DO
350    CONTINUE
    END DO
    !------------------------------------------------------------------------------!
    !        TOTAL NUMBER OF SPECIAL POINTS: NTOT
    !        BEFORE USING THE LIST WVKL AS WAVE VECTORS, THEY HAVE TO BE
    !        MULTIPLIED BY 2*PI
    !        THE LIST OF WEIGHTS LWGHT IS NOT NORMALIZED
    !------------------------------------------------------------------------------!
    IF (ntot>nkpoint .AND. iout>0) THEN
       WRITE (iout,*) ' K290| In sppt2 number of special points = ', ntot
       WRITE (iout,*) ' K290| but nkpoint = ', nkpoint
       ntot = -1
    END IF
    IF (iout>0) THEN
       !       Write the index table relating k points in the mesh
       !       with special k points
       WRITE (iout,'(/,A)') &
            ' K290| Cross table relating mesh points with special points:'
       WRITE (iout,'(5(4X,"IK -> SK"))')
       DO i = 1, imesh
          iplace = includ(i)/2
          WRITE (iout,'(1X,I5,1X,I5)',ADVANCE="NO") i, iplace
          IF (MOD(i,5)==0) WRITE (iout,*)
       END DO
       IF (MOD(j-1,5)/=0) WRITE (iout,*)
    END IF
    RETURN
    !------------------------------------------------------------------------------!
    !        ERROR MESSAGES
    !------------------------------------------------------------------------------!
450 CONTINUE
    WRITE (*,'(A,/)') ' SUBROUTINE SPPT2 *** FATAL ERROR ***'
    WRITE (*,'(A,3F10.4,/,A,3F10.4,A,/,A,I3,A)') ' THE VECTOR     ', &
         wva, ' GENERATED FROM ', wvk, ' IN THE BASIC MESH', &
         ' BY ROTATION NO. ', ibrav(iop), ' IS OUTSIDE THE 1BZ'
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    !------------------------------------------------------------------------------!
470 CONTINUE
    WRITE (*,'(A,/)') ' SUBROUTINE SPPT2 *** FATAL ERROR ***'
    WRITE (*,*) 'MESH SIZE EXCEEDS NKPOINT=', nkpoint
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    !------------------------------------------------------------------------------!
490 CONTINUE
    WRITE (*,'(A,/)') ' SUBROUTINE SPPT2 *** FATAL ERROR ***'
    WRITE (*,'(A,3F10.4,/,A,3F10.4,A,/,A,I3,A)') ' THE VECTOR     ', &
         wva, ' GENERATED FROM ', wvk, ' IN THE BASIC MESH', &
         ' BY ROTATION NO. ', ib(n), ' IS NOT IN THE LIST'
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END SUBROUTINE sppt2
! *****************************************************************************
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
!> \param error ...
!> \par Original Header
!>         MESH MAINTAINS A LIST OF VECTORS FOR PLACEMENT AND/OR LOOKUP
!>
!>         ADDITIONAL ENTRY POINTS: REMOVE .... REMOVE VECTOR FROM LIST
!>                                  GARBAG .... WAS VECTOR REMOVED ?
!>
!>         WVK ....... VECTOR
!>         IPLACE .... ON INPUT: -2 MEANS: INITIALIZE  THE LIST
!>                                         (AND RETURN)
!>                               -1 MEANS: FIND WVK IN THE LIST
!>                                0 MEANS: ADD  WVK TO THE LIST
!>                               >0 MEANS: RETURN WVK NO. IPLACE
!>                    ON OUTPUT: THE POSITION ASSIGNED TO WVK
!>                               (=0 IF WVK IS NOT IN THE LIST)
!>         DELTA      REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
! *****************************************************************************
  SUBROUTINE mesh(iout,wvk,iplace,igarb0,igarbg,nmesh,nhash,list,rlist, &
       delta,error)

    INTEGER                                  :: iout
    REAL(KIND=dp)                            :: wvk(3)
    INTEGER                                  :: iplace, igarb0, igarbg, &
                                                nmesh, nhash, &
                                                list(nhash+nmesh)
    REAL(KIND=dp)                            :: rlist(3,nmesh), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mesh', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nil = 0

    INTEGER                                  :: i, ihash, ipoint, j
    INTEGER, SAVE                            :: istore = 0
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: delta1, rhash

    failure = .FALSE.
    delta1 = 10.0_dp*delta
    IF (iplace<=-2) THEN
       DO i = 1, nhash + nmesh
          list(i) = nil
       END DO
       istore = 1
       !       IGARB0 points to a linked list of removed WVKS (the garbage).
       igarb0 = 0
       igarbg = 0
       RETURN
       !------------------------------------------------------------------------------!
    ELSE IF ((iplace>-2) .AND. (iplace<=0)) THEN
       !       The particular HASH function used in this case:
       rhash = 0.7890E0_dp*wvk(1) + 0.6810E0_dp*wvk(2) + 0.5811E0_dp*wvk(3) + delta
       ihash = INT(ABS(rhash)*REAL(nhash,KIND=dp))
       ihash = MOD(ihash,nhash) + nmesh + 1
       !       Search for WVK in linked list
       ipoint = list(ihash)
       DO i = 1, 100
          !         List exhausted
          IF (ipoint==nil) GO TO 130
          !         Compare WVK with this element
          DO j = 1, 3
             IF (ABS(wvk(j)-rlist(j,ipoint))>delta1) GO TO 115
          END DO
          !         WVK located
          GO TO 160
          !         Next element of list
115       CONTINUE
          ihash = ipoint
          ipoint = list(ihash)
       END DO
       !       List too long
       WRITE (*,'(2A,/,A)') &
            ' SUBROUTINE MESH *** FATAL ERROR *** LINKED LIST', &
            ' TOO LONG ***', ' CHOOSE A BETTER HASH-FUNCTION'
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !       WVK was not found
130    CONTINUE
       IF (iplace==-1) THEN
          !         IPLACE=-1 : search for WVK unsuccessful
          iplace = 0
          RETURN
       ELSE
          !         IPLACE=0: add WVK to the list
          list(ihash) = istore
          IF (istore>nmesh) THEN
             WRITE (*,'(A)') 'SUBROUTINE MESH *** FATAL ERROR ***'
             WRITE (*,'(A,I10,A,/,A,3F10.5)') ' ISTORE=', istore, &
                  ' EXCEEDS DIMENSIONS', ' WVK = ', wvk
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
          list(istore) = nil
          DO i = 1, 3
             rlist(i,istore) = wvk(i)
          END DO
          istore = istore + 1
          iplace = istore - 1
          RETURN
       END IF
       !       WVK was found
160    CONTINUE
       IF (iplace==0) RETURN
       !       IPLACE=-1
       iplace = ipoint
       RETURN
    ELSE
       !------------------------------------------------------------------------------!
       !        Return a wavevector (IPLACE > 0)
       !------------------------------------------------------------------------------!
       ipoint = iplace
       IF (ipoint>=istore) GO TO 190
       DO i = 1, 3
          wvk(i) = rlist(i,ipoint)
       END DO
       RETURN
    END IF
    !------------------------------------------------------------------------------!
    !        Error - beyond list
    !------------------------------------------------------------------------------!
190 CONTINUE
    IF (iout>0) WRITE (iout,'(A,/,A,I5,A,/)') &
         ' SUBROUTINE MESH *** WARNING ***', ' IPLACE = ', iplace, &
         ' IS BEYOND THE LISTS - WVK SET TO 1.0E38'
    DO i = 1, 3
       wvk(i) = 1.0E38_dp
    END DO
  END SUBROUTINE mesh

! *****************************************************************************
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
!> \param error ...
!> \par Original Header
!>         ENTRY POINT FOR REMOVING A WAVEVECTOR
!>
!>         INPUT:
!>           WVK(3)
!>           DELTA   REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>         OUTPUT:
!>           IPLACE ..... 1 IF WVK WAS REMOVED
!>                        0 IF WVK WAS NOT REMOVED
!>                          (WVK NOT IN THE LINKED LISTS)
! *****************************************************************************
  SUBROUTINE remove(wvk,iplace,igarb0,igarbg,nmesh,nhash,list,rlist,delta,error)

    REAL(KIND=dp)                            :: wvk(3)
    INTEGER                                  :: iplace, igarb0, igarbg, &
                                                nmesh, nhash, &
                                                list(nhash+nmesh)
    REAL(KIND=dp)                            :: rlist(3,nmesh), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'remove', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nil = 0

    INTEGER                                  :: i, ihash, ipoint, j
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: delta1, rhash

    failure = .FALSE.
    delta1 = 10.0_dp*delta
    !     The particular hash function used in this case:
    rhash = 0.7890E0_dp*wvk(1) + 0.6810E0_dp*wvk(2) + 0.5811E0_dp*wvk(3) + delta
    ihash = INT(ABS(rhash)*REAL(nhash,KIND=dp))
    ihash = MOD(ihash,nhash) + nmesh + 1
    !     Search for WVK in linked list
    ipoint = list(ihash)
    DO i = 1, 100
       !       List exhausted
       IF (ipoint==nil) THEN
          !         WVK was not found in the mesh:
          iplace = 0
          RETURN
       END IF
       !       Compare WVK with this element
       DO j = 1, 3
          IF (ABS(wvk(j)-rlist(j,ipoint))>delta1) GO TO 215
       END DO
       !       WVK located, now remove it from the list:
       list(ihash) = list(ipoint)
       !       LIST(IHASH) now points to the next element in the list,
       !       and the present WVK has become garbage.
       !       Add WVK to the list of garbage:
       IF (igarb0==0) THEN
          !         Start up the garbage list:
          igarb0 = ipoint
       ELSE
          list(igarbg) = ipoint
       END IF
       igarbg = ipoint
       list(igarbg) = nil
       iplace = 1
       RETURN
       !       Next element of list
215    CONTINUE
       ihash = ipoint
       ipoint = list(ihash)
    END DO
    !     List too long
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END SUBROUTINE remove

! *****************************************************************************
!> \brief ...
!> \param wvk ...
!> \param iplace ...
!> \param igarb0 ...
!> \param nmesh ...
!> \param nhash ...
!> \param list ...
!> \param rlist ...
!> \param delta ...
!> \param error ...
!> \par Original Header
!>         ENTRY POINT FOR CHECKING IF A WAVEVECTOR
!>                     IS IN THE GARBAGE LIST
!>         INPUT:
!>           WVK(3)
!>           DELTA      REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
!>
!>         OUTPUT:
!>           IPLACE  ..... I > 0 IS THE PLACE IN THE GARBAGE LIST
!>                             0 IF WVK NOT AMONG THE GARBAGE
! *****************************************************************************
  SUBROUTINE garbag(wvk,iplace,igarb0,nmesh,nhash,list,rlist,delta,error)
    REAL(KIND=dp)                            :: wvk(3)
    INTEGER                                  :: iplace, igarb0, nmesh, nhash, &
                                                list(nhash+nmesh)
    REAL(KIND=dp)                            :: rlist(3,nmesh), delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'garbag', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nil = 0

    INTEGER                                  :: i, ihash, ipoint, j
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: delta1

    failure = .FALSE.
    delta1 = 10.0_dp*delta
    !     Search for WVK in linked list
    !     Point to the garbage list
    ipoint = igarb0
    DO i = 1, nmesh
       !       LIST EXHAUSTED
       IF (ipoint==nil) THEN
          !         WVK was not found in the mesh:
          iplace = 0
          RETURN
       END IF
       !       Compare WVK with this element
       DO j = 1, 3
          IF (ABS(wvk(j)-rlist(j,ipoint))>delta1) THEN
             !       Next element of list
             ihash = ipoint
             ipoint = list(ihash)
             EXIT
          ENDIF
       END DO
       !       WVK was located in the garbage list
       iplace = i
       RETURN

    END DO
    !     List too long
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
  END SUBROUTINE garbag

! *****************************************************************************
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
!> \param error ...
!> \par Original Header
!>         REDUCE WVK TO LIE ENTIRELY WITHIN THE 1ST BRILLOUIN ZONE
!>         BY ADDING B-VECTORS
!>         DELTA      REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
! *****************************************************************************
  SUBROUTINE bzrduc(wvk,a1,a2,a3,b1,b2,b3,rsdir,nrsdir,nplane,delta,error)
    REAL(KIND=dp)                            :: wvk(3), a1(3), a2(3), a3(3), &
                                                b1(3), b2(3), b3(3)
    INTEGER                                  :: nrsdir
    REAL(KIND=dp)                            :: rsdir(4,nrsdir)
    INTEGER                                  :: nplane
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bzrduc', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nzones = 4, &
                                                nnn = 2*nzones + 1, &
                                                nn = nzones + 1, yes = 1

    INTEGER                                  :: i, i1, i2, i3, n1, n2, n3, &
                                                nn1, nn2, nn3
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: wb(3), wva(3)

    failure = .FALSE.
    IF (inbz(wvk,rsdir,nrsdir,nplane,delta)==yes) RETURN
    !------------------------------------------------------------------------------!
    !     Look around +/- "NZONES" to locate vector
    !     NZONES may need to be increased for very anisotropic zones
    !------------------------------------------------------------------------------!
    !     WVK already inside 1Bz
    !------------------------------------------------------------------------------!
    !     Express WVK in the basis of B1,2,3.
    !     This permits an estimate of how far WVK is from the 1Bz.
    wb(1) = wvk(1)*a1(1) + wvk(2)*a1(2) + wvk(3)*a1(3)
    wb(2) = wvk(1)*a2(1) + wvk(2)*a2(2) + wvk(3)*a2(3)
    wb(3) = wvk(1)*a3(1) + wvk(2)*a3(2) + wvk(3)*a3(3)
    nn1 = NINT(wb(1))
    nn2 = NINT(wb(2))
    nn3 = NINT(wb(3))
    !     Look around the estimated vector for the one truly inside the 1Bz
    DO n1 = 1, nnn
       i1 = nn - n1 - nn1
       DO n2 = 1, nnn
          i2 = nn - n2 - nn2
          DO n3 = 1, nnn
             i3 = nn - n3 - nn3
             DO i = 1, 3
                wva(i) = wvk(i) + REAL(i1,KIND=dp)*b1(i) + REAL(i2,KIND=dp)*b2(i) + &
                     REAL(i3,kind=dp)*b3(i)
             END DO
             IF (inbz(wva,rsdir,nrsdir,nplane,delta)==yes) THEN
             !------------------------------------------------------------------------------!
             !     The reduced vector
                DO i = 1, 3
                   wvk(i) = wva(i)
                END DO
                RETURN
             ENDIF
          END DO
       END DO
    END DO
    !------------------------------------------------------------------------------!
    !     Fatal error
    WRITE (*,'(A,/,A,3F10.4,A)') ' SUBROUTINE BZRDUC *** FATAL ERROR ***', &
         ' WAVEVECTOR ', wvk, ' COULD NOT BE REDUCED TO THE 1BZ'
    CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)

  END SUBROUTINE bzrduc
! *****************************************************************************
!> \brief ...
!> \param wvk ...
!> \param rsdir ...
!> \param nrsdir ...
!> \param nplane ...
!> \param delta ...
!> \retval inbz ...
!> \par Original Header
!>         IS WVK IN THE 1ST BRILLOUIN ZONE ?
!>         CHECK WHETHER WVK LIES INSIDE ALL THE PLANES
!>         THAT DEFINE THE 1BZ.
!>         DELTA      REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
! *****************************************************************************
  FUNCTION inbz(wvk,rsdir,nrsdir,nplane,delta)
    REAL(KIND=dp)                            :: wvk(3)
    INTEGER                                  :: nrsdir
    REAL(KIND=dp)                            :: rsdir(4,nrsdir)
    INTEGER                                  :: nplane
    REAL(KIND=dp)                            :: delta
    INTEGER                                  :: inbz

    INTEGER, PARAMETER                       :: no = 0, yes = 1

    INTEGER                                  :: n
    REAL(KIND=dp)                            :: projct

!------------------------------------------------------------------------------!

    inbz = no
    DO n = 1, nplane
       projct = (rsdir(1,n)*wvk(1)+rsdir(2,n)*wvk(2)+rsdir(3,n)*wvk(3))/ &
            rsdir(4,n)
       !       WVK is outside the Bz
       IF (ABS(projct)>0.5E0_dp+delta) RETURN
    END DO
    inbz = yes
  END FUNCTION inbz

! *****************************************************************************
!> \brief ...
!> \param iout ...
!> \param b1 ...
!> \param b2 ...
!> \param b3 ...
!> \param rsdir ...
!> \param nrsdir ...
!> \param nplane ...
!> \param delta ...
!> \param error ...
!> \par Original Header
!>         FIND THE VECTORS WHOSE HALVES DEFINE THE 1ST BRILLOUIN ZONE
!>         OUTPUT:
!>         NPLANE TELLS HOW MANY ELEMENTS OF RSDIR CONTAIN
!>                NORMAL VECTORS DEFINING THE PLANES
!>         METHOD: STARTING WITH THE PARALLELOPIPED SPANNED BY B1,2,3
!>         AROUND THE ORIGIN, VECTORS INSIDE A SUFFICIENTLY LARGE
!>         SPHERE ARE TESTED TO SEE WHETHER THE PLANES AT 1/2*B WILL
!>         FURTHER CONFINE THE 1BZ.
!>         THE RESULTING VECTORS ARE NOT CLEANED TO AVOID REDUNDANT
!>         PLANES.
!>         DELTA      REQUIRED ACCURACY (1.0E-6_dp IS A GOOD VALUE)
! *****************************************************************************
  SUBROUTINE bzdefi(iout,b1,b2,b3,rsdir,nrsdir,nplane,delta,error)
    INTEGER                                  :: iout
    REAL(KIND=dp)                            :: b1(3), b2(3), b3(3)
    INTEGER                                  :: nrsdir
    REAL(KIND=dp)                            :: rsdir(4,nrsdir)
    INTEGER                                  :: nplane
    REAL(KIND=dp)                            :: delta
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'bzdefi', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, i2, i3, initlz, j, n, &
                                                n1, n2, n3, nb1, nb2, nb3, &
                                                nnb1, nnb2, nnb3
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: b1len, b2len, b3len, bmax, &
                                                bvec(3), projct

    DATA initlz/0/
    failure = .FALSE.
    !------------------------------------------------------------------------------!
    IF (initlz/=0) RETURN
    !     Once initialized, we do not repeat the calculation
    initlz = 1
    b1len = b1(1)**2 + b1(2)**2 + b1(3)**2
    b2len = b2(1)**2 + b2(2)**2 + b2(3)**2
    b3len = b3(1)**2 + b3(2)**2 + b3(3)**2
    !     Lattice containing entirely the brillouin zone
    bmax = b1len + b2len + b3len
    nb1 = INT(SQRT(bmax/b1len)+delta) + 1
    nb2 = INT(SQRT(bmax/b2len)+delta) + 1
    nb3 = INT(SQRT(bmax/b3len)+delta) + 1
    !     PRINT *,'NB1,2,3 = ',NB1,NB2,NB3
    DO i = 1, nrsdir
       DO j = 1, 4
          rsdir(j,i) = 0.0_dp
       END DO
    END DO
    !     1Bz is certainly confined inside the 1/2(B1,B2,B3) parallelopiped
    DO i = 1, 3
       rsdir(i,1) = b1(i)
       rsdir(i,2) = b2(i)
       rsdir(i,3) = b3(i)
    END DO
    rsdir(4,1) = b1len
    rsdir(4,2) = b2len
    rsdir(4,3) = b3len
    !     Starting confinement: 3 planes
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
             IF (i1==0 .AND. i2==0 .AND. i3==0) CYCLE
             DO i = 1, 3
                bvec(i) = REAL(i1,KIND=dp)*b1(i) + REAL(i2,KIND=dp)*b2(i) + &
                     REAL(i3,KIND=dp)*b3(i)
             END DO
             !           Does the plane of 1/2*BVEC narrow down the 1Bz ?
             DO n = 1, nplane
                projct = 0.5E0_dp*(rsdir(1,n)*bvec(1)+rsdir(2,n)*bvec(2)+rsdir(3, &
                     n)*bvec(3))/rsdir(4,n)
                !             1/2*BVEC is outside the Bz - skip this direction
                !             The 1.0E-6_dp takes care of single points touching the Bz,
                !             and of the -(plane)
                IF (ABS(projct)>0.5E0_dp-delta) EXIT
             END DO
             IF (ABS(projct)>0.5E0_dp-delta) CYCLE
             !           1/2*BVEC further confines the 1Bz - include into RSDIR
             nplane = nplane + 1
             !           PRINT *,NPLANE,' PLANE INCLUDED, I1,2,3 = ',I1,I2,I3
             IF (nplane>nrsdir) THEN
                !------------------------------------------------------------------------------!
                !     Error messages
                IF (iout>0) THEN
                   WRITE (iout,'(A)') ' SUBROUTINE BZDEFI *** FATAL ERROR ***'
                   WRITE (iout,'(" TOO MANY PLANES, NRSDIR = ",I5)') nrsdir
                END IF
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                RETURN
             ENDIF
             DO i = 1, 3
                rsdir(i,nplane) = bvec(i)
             END DO
             !           Length squared
             rsdir(4,nplane) = bvec(1)**2 + bvec(2)**2 + bvec(3)**2
          END DO
       END DO
    END DO
    !------------------------------------------------------------------------------!
    !     Print information
    IF (iout>0) WRITE (iout,'(A,I3,A,/,A,/,100(1X,3F10.4,/))') &
         ' THE 1ST BRILLOUIN ZONE IS CONFINED BY (AT MOST)', nplane, &
         ' PLANES', ' AS DEFINED BY THE +/- HALVES OF THE VECTORS:', &
         ((rsdir(i,n),i=1,3),n=1,nplane)

  END SUBROUTINE bzdefi

END MODULE k290
