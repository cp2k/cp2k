!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief K-points and crystal symmetry routines
!> \author jgh
! **************************************************************************************************
MODULE cryssym

   USE kinds,                           ONLY: dp
   USE spglib_f08,                      ONLY: &
        spg_get_international, spg_get_ir_reciprocal_mesh, spg_get_major_version, &
        spg_get_micro_version, spg_get_minor_version, spg_get_multiplicity, spg_get_pointgroup, &
        spg_get_schoenflies, spg_get_symmetry
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: csym_type, release_csym_type, print_crys_symmetry, print_kp_symmetry
   PUBLIC :: crys_sym_gen, kpoint_gen, apply_rotation_coord

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cryssym'

! **************************************************************************************************
!> \brief CSM type
!> \par   Content:
!>
! **************************************************************************************************
   TYPE csym_type
      LOGICAL                                     :: symlib = .FALSE.
      LOGICAL                                     :: fullgrid = .FALSE.
      INTEGER                                     :: plevel = 0
      INTEGER                                     :: punit = -1
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: delta = 1.0e-8_dp
      REAL(KIND=dp), DIMENSION(3, 3)              :: hmat
      ! KPOINTS
      REAL(KIND=dp), DIMENSION(3)                 :: wvk0 = 0.0_dp
      INTEGER, DIMENSION(3)                       :: mesh
      INTEGER                                     :: nkpoint
      INTEGER                                     :: nat
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: scoord
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: xkpoint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: wkpoint
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: kpmesh
      INTEGER, DIMENSION(:, :), ALLOCATABLE       :: kplink
      INTEGER, DIMENSION(:), ALLOCATABLE          :: kpop
      !SPGLIB
      CHARACTER(len=11)                           :: international_symbol
      CHARACTER(len=6)                            :: pointgroup_symbol
      CHARACTER(len=10)                           :: schoenflies
      INTEGER                                     :: n_operations
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
   END TYPE csym_type

CONTAINS

! **************************************************************************************************
!> \brief Release the CSYM type
!> \param csym  The CSYM type
! **************************************************************************************************
   SUBROUTINE release_csym_type(csym)
      TYPE(csym_type)                                    :: csym

      CHARACTER(LEN=*), PARAMETER :: routineN = 'release_csym_type', &
         routineP = moduleN//':'//routineN

      IF (ALLOCATED(csym%rotations)) THEN
         DEALLOCATE (csym%rotations)
      END IF
      IF (ALLOCATED(csym%translations)) THEN
         DEALLOCATE (csym%translations)
      END IF
      IF (ALLOCATED(csym%atype)) THEN
         DEALLOCATE (csym%atype)
      END IF
      IF (ALLOCATED(csym%scoord)) THEN
         DEALLOCATE (csym%scoord)
      END IF
      IF (ALLOCATED(csym%xkpoint)) THEN
         DEALLOCATE (csym%xkpoint)
      END IF
      IF (ALLOCATED(csym%wkpoint)) THEN
         DEALLOCATE (csym%wkpoint)
      END IF
      IF (ALLOCATED(csym%kpmesh)) THEN
         DEALLOCATE (csym%kpmesh)
      END IF
      IF (ALLOCATED(csym%kplink)) THEN
         DEALLOCATE (csym%kplink)
      END IF
      IF (ALLOCATED(csym%kpop)) THEN
         DEALLOCATE (csym%kpop)
      END IF

   END SUBROUTINE release_csym_type

! **************************************************************************************************
!> \brief ...
!> \param csym ...
!> \param scoor ...
!> \param types ...
!> \param hmat ...
!> \param delta ...
!> \param iounit ...
! **************************************************************************************************
   SUBROUTINE crys_sym_gen(csym, scoor, types, hmat, delta, iounit)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: scoor
      INTEGER, DIMENSION(:), INTENT(IN)                  :: types
      REAL(KIND=dp), INTENT(IN)                          :: hmat(3, 3)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: delta
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'crys_sym_gen', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ierr, major, micro, minor, nat, &
                                                            nop, tra_mat(3, 3)
      LOGICAL                                            :: spglib

      CALL timeset(routineN, handle)

      !..total number of atoms
      nat = SIZE(scoor, 2)
      csym%nat = nat

      ! output unit
      IF (PRESENT(iounit)) THEN
         csym%punit = iounit
      ELSE
         csym%punit = -1
      END IF

      ! accuracy for symmetry
      IF (PRESENT(delta)) THEN
         csym%delta = delta
      ELSE
         csym%delta = 1.e-6_dp
      END IF

      !..set cell values
      csym%hmat = hmat

      ! atom types
      ALLOCATE (csym%atype(nat))
      csym%atype(1:nat) = types(1:nat)

      ! scaled coordinates
      ALLOCATE (csym%scoord(3, nat))
      csym%scoord(1:3, 1:nat) = scoor(1:3, 1:nat)

      !..try spglib
      major = spg_get_major_version()
      minor = spg_get_minor_version()
      micro = spg_get_micro_version()
      IF (major == 0) THEN
         CALL cp_warn(__LOCATION__, "Symmetry library SPGLIB not available")
         spglib = .FALSE.
      ELSE
         IF (major /= 1 .OR. minor /= 9 .OR. micro /= 9) THEN
            CALL cp_warn(__LOCATION__, "Version of Symmetry Library SPGLIB not tested")
         END IF
         spglib = .TRUE.
         ierr = spg_get_international(csym%international_symbol, TRANSPOSE(hmat), scoor, types, nat, delta)
         IF (ierr == 0) THEN
            CALL cp_warn(__LOCATION__, "Symmetry Library SPGLIB failed")
            spglib = .FALSE.
         ELSE
            nop = spg_get_multiplicity(TRANSPOSE(hmat), scoor, types, nat, delta)
            ALLOCATE (csym%rotations(3, 3, nop), csym%translations(3, nop))
            csym%n_operations = nop
            ierr = spg_get_symmetry(csym%rotations, csym%translations, nop, &
                                    TRANSPOSE(hmat), scoor, types, nat, delta)
            ! Schoenflies Symbol
            csym%schoenflies = ' '
            ierr = spg_get_schoenflies(csym%schoenflies, TRANSPOSE(hmat), scoor, types, nat, delta)
            ! Point Group
            csym%pointgroup_symbol = ' '
            tra_mat = 0
            ierr = spg_get_pointgroup(csym%pointgroup_symbol, tra_mat, &
                                      csym%rotations, csym%n_operations)
         END IF
      END IF
      csym%symlib = spglib

      CALL timestop(handle)

   END SUBROUTINE crys_sym_gen

! **************************************************************************************************
!> \brief ...
!> \param csym ...
!> \param nk ...
!> \param symm ...
!> \param shift ...
!> \param full_grid ...
! **************************************************************************************************
   SUBROUTINE kpoint_gen(csym, nk, symm, shift, full_grid)
      TYPE(csym_type)                                    :: csym
      INTEGER, INTENT(IN)                                :: nk(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: symm
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: shift(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_grid

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_gen', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ik, is_shift(3), &
                                                            is_time_reversal, j, nkp, nkpts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwkp, xptr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ixkp
      LOGICAL                                            :: fullmesh
      REAL(KIND=dp)                                      :: rxkp(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp

      CALL timeset(routineN, handle)

      IF (PRESENT(shift)) csym%wvk0 = shift

      IF (PRESENT(symm)) THEN
         csym%istriz = -1
         IF (symm) csym%istriz = 1
      END IF

      IF (PRESENT(full_grid)) THEN
         fullmesh = full_grid
      ELSE
         fullmesh = .FALSE.
      END IF
      csym%fullgrid = fullmesh

      csym%nkpoint = 0
      csym%mesh(1:3) = nk(1:3)

      nkpts = nk(1)*nk(2)*nk(3)
      ALLOCATE (xkp(3, nkpts), wkp(nkpts))
      ! kp: link
      ALLOCATE (csym%kplink(2, nkpts))
      csym%kplink = 0

      ! go through all the options
      IF (csym%symlib) THEN
         ! symmetry library is available
         IF (fullmesh) THEN
            ! full mesh requested
            CALL full_grid_gen(nk, xkp, wkp, shift)
            IF (csym%istriz == 1) THEN
               ! use inversion symmetry
               CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
            ELSE
               ! full kpoint mesh is used
            END IF
         ELSE
            ! use symmtry operations
            IF (SUM(ABS(csym%wvk0)) /= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, "MacDonald shifted k-point meshes are only "// &
                             "possible without symmetrization.")
            END IF
            ALLOCATE (ixkp(3, nkpts), iwkp(nkpts))
            is_shift = 0
            is_time_reversal = 1
            nkp = spg_get_ir_reciprocal_mesh(ixkp, iwkp, nk, is_shift, is_time_reversal, &
                                             TRANSPOSE(csym%hmat), csym%scoord, csym%atype, &
                                             csym%nat, csym%delta)
            wkp = 0.0_dp
            DO i = 1, nkpts
               xkp(:, i) = REAL(is_shift+2*ixkp(:, i), KIND=dp)/REAL(2*nk(:), KIND=dp)
               j = iwkp(i)+1
               wkp(j) = wkp(j)+1.0_dp
            ENDDO
            csym%kplink(1, 1:nkpts) = iwkp(1:nkpts)+1
            DEALLOCATE (ixkp, iwkp)
         END IF
      ELSE
         ! no symmetry library is available
         CALL full_grid_gen(nk, xkp, wkp, shift)
         IF (csym%istriz /= 1 .AND. fullmesh) THEN
            ! full kpoint mesh is used
            DO i = 1, nkpts
               csym%kplink(1, i) = i
            END DO
         ELSE
            ! use inversion symmetry
            CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
         END IF
      END IF
      ! count kpoints
      nkp = 0
      DO i = 1, nkpts
         IF (wkp(i) > 0.0_dp) nkp = nkp+1
      END DO

      ! store reduced kpoint set
      csym%nkpoint = nkp
      ALLOCATE (csym%xkpoint(3, nkp), csym%wkpoint(nkp))
      ALLOCATE (xptr(nkp))
      j = 0
      DO ik = 1, nkpts
         IF (wkp(ik) > 0.0_dp) THEN
            j = j+1
            csym%wkpoint(j) = wkp(ik)
            csym%xkpoint(1:3, j) = xkp(1:3, ik)
            xptr(j) = ik
         END IF
      END DO
      CPASSERT(j == nkp)

      ! kp: mesh
      ALLOCATE (csym%kpmesh(3, nkpts))
      csym%kpmesh(1:3, 1:nkpts) = xkp(1:3, 1:nkpts)

      ! kp: link
      DO j = 1, nkp
         DO ik = 1, nkpts
            i = csym%kplink(1, ik)
            IF (i == xptr(j)) THEN
               csym%kplink(2, ik) = j
            END IF
         END DO
      END DO
      DEALLOCATE (xptr)

      ! kp: operations
      ALLOCATE (csym%kpop(nkpts))
      DO ik = 1, nkpts
         IF (wkp(ik) > 0.0_dp) THEN
            csym%kpop(ik) = 1
         ELSE
            IF (csym%symlib .AND. .NOT. fullmesh) THEN
               j = csym%kplink(2, ik)
               DO i = 1, csym%n_operations
                  rxkp(1:3) = kp_apply_operation(csym%xkpoint(1:3, j), csym%rotations(:, :, i))
                  IF (ALL((ABS(xkp(:, ik)-rxkp(:))) < 1.e-12_dp)) THEN
                     csym%kpop(ik) = i
                     EXIT
                  END IF
               END DO
            ELSE
               ! only time reversal symmetry
               csym%kpop(ik) = 2
            END IF
         END IF
      END DO

      DEALLOCATE (xkp, wkp)

      CALL timestop(handle)

   END SUBROUTINE kpoint_gen

! **************************************************************************************************
!> \brief ...
!> \param nk ...
!> \param xkp ...
!> \param wkp ...
!> \param shift ...
! **************************************************************************************************
   SUBROUTINE full_grid_gen(nk, xkp, wkp, shift)
      INTEGER, INTENT(IN)                                :: nk(3)
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      REAL(KIND=dp), INTENT(IN)                          :: shift(3)

      INTEGER                                            :: i, ix, iy, iz
      REAL(KIND=dp)                                      :: kpt_latt(3)

      wkp = 0.0_dp
      i = 0
      DO ix = 1, nk(1)
         DO iy = 1, nk(2)
            DO iz = 1, nk(3)
               i = i+1
               kpt_latt(1) = REAL(2*ix-nk(1)-1, KIND=dp)/(2._dp*REAL(nk(1), KIND=dp))
               kpt_latt(2) = REAL(2*iy-nk(2)-1, KIND=dp)/(2._dp*REAL(nk(2), KIND=dp))
               kpt_latt(3) = REAL(2*iz-nk(3)-1, KIND=dp)/(2._dp*REAL(nk(3), KIND=dp))
               xkp(1:3, i) = kpt_latt(1:3)
               wkp(i) = 1.0_dp
            END DO
         END DO
      END DO
      DO i = 1, nk(1)*nk(2)*nk(3)
         xkp(1:3, i) = xkp(1:3, i)+shift(1:3)
      END DO

   END SUBROUTINE full_grid_gen

! **************************************************************************************************
!> \brief ...
!> \param xkp ...
!> \param wkp ...
!> \param link ...
! **************************************************************************************************
   SUBROUTINE inversion_symm(xkp, wkp, link)
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: link

      INTEGER                                            :: i, j, nkpts

      nkpts = SIZE(wkp, 1)

      link(:) = 0
      DO i = 1, nkpts
         IF (link(i) == 0) link(i) = i
         DO j = i+1, nkpts
            IF (wkp(j) == 0) CYCLE
            IF (ALL(xkp(:, i) == -xkp(:, j))) THEN
               wkp(i) = wkp(i)+wkp(j)
               wkp(j) = 0.0_dp
               link(j) = i
               EXIT
            END IF
         END DO
      END DO

   END SUBROUTINE inversion_symm

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \param r ...
!> \retval y ...
! **************************************************************************************************
   FUNCTION kp_apply_operation(x, r) RESULT(y)
      REAL(KIND=dp), INTENT(IN)                          :: x(3)
      INTEGER, INTENT(IN)                                :: r(3, 3)
      REAL(KIND=dp)                                      :: y(3)

      y(1) = REAL(r(1, 1), dp)*x(1)+REAL(r(1, 2), dp)*x(2)+REAL(r(1, 3), dp)*x(3)
      y(2) = REAL(r(2, 1), dp)*x(1)+REAL(r(2, 2), dp)*x(2)+REAL(r(2, 3), dp)*x(3)
      y(3) = REAL(r(3, 1), dp)*x(1)+REAL(r(3, 2), dp)*x(2)+REAL(r(3, 3), dp)*x(3)

   END FUNCTION kp_apply_operation

! **************************************************************************************************
!> \brief ...
!> \param f0 ...
!> \param csym ...
!> \param ir ...
! **************************************************************************************************
   SUBROUTINE apply_rotation_coord(f0, csym, ir)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: f0
      TYPE(csym_type)                                    :: csym
      INTEGER, INTENT(IN)                                :: ir

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_rotation_coord', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ia, ib, natom
      REAL(KIND=dp)                                      :: diff
      REAL(KIND=dp), DIMENSION(3)                        :: rb, ri, ro, tr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rot

      natom = csym%nat
      rot(1:3, 1:3) = csym%rotations(1:3, 1:3, ir)
      tr(1:3) = csym%translations(1:3, ir)

      f0 = 0
      DO ia = 1, natom
         ri(1:3) = csym%scoord(1:3, ia)
         ro(1) = REAL(rot(1, 1), dp)*ri(1)+REAL(rot(2, 1), dp)*ri(2)+REAL(rot(3, 1), dp)*ri(3)+tr(1)
         ro(2) = REAL(rot(1, 2), dp)*ri(1)+REAL(rot(2, 2), dp)*ri(2)+REAL(rot(3, 2), dp)*ri(3)+tr(2)
         ro(3) = REAL(rot(1, 3), dp)*ri(1)+REAL(rot(2, 3), dp)*ri(2)+REAL(rot(3, 3), dp)*ri(3)+tr(3)
         DO ib = 1, natom
            rb(1:3) = csym%scoord(1:3, ib)
            diff = SQRT(SUM((ri(:)-rb(:))**2))
            IF (diff < csym%delta) THEN
               f0(ia) = ib
               EXIT
            END IF
         END DO
         CPASSERT(f0(ia) /= 0)
      END DO

   END SUBROUTINE apply_rotation_coord

! **************************************************************************************************
!> \brief ...
!> \param csym ...
! **************************************************************************************************
   SUBROUTINE print_crys_symmetry(csym)
      TYPE(csym_type)                                    :: csym

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_crys_symmetry', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iunit, j, plevel

      iunit = csym%punit
      IF (iunit >= 0) THEN
         plevel = csym%plevel
         WRITE (iunit, "(/,T2,A)") "Crystal Symmetry Information"
         IF (csym%symlib) THEN
            WRITE (iunit, '(A,T71,A10)') "       International Symbol: ", ADJUSTR(TRIM(csym%international_symbol))
            WRITE (iunit, '(A,T71,A10)') "       Point Group Symbol: ", ADJUSTR(TRIM(csym%pointgroup_symbol))
            WRITE (iunit, '(A,T71,A10)') "       Schoenflies Symbol: ", ADJUSTR(TRIM(csym%schoenflies))
            !
            WRITE (iunit, '(A,T71,I10)') "       Number of Symmetry Operations: ", csym%n_operations
            IF (plevel > 0) THEN
               DO i = 1, csym%n_operations
                  WRITE (iunit, '(A,i4,T51,3I10,/,T51,3I10,/,T51,3I10)') &
                     "           Rotation #: ", i, (csym%rotations(j, :, i), j=1, 3)
                  WRITE (iunit, '(T36,3F15.7)') csym%translations(:, i)
               END DO
            END IF
         ELSE
            WRITE (iunit, "(T2,A)") "SPGLIB for Crystal Symmetry Information determination is not availale"
         END IF
      END IF

   END SUBROUTINE print_crys_symmetry

! **************************************************************************************************
!> \brief ...
!> \param csym ...
! **************************************************************************************************
   SUBROUTINE print_kp_symmetry(csym)
      TYPE(csym_type), INTENT(IN)                        :: csym

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_kp_symmetry', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iunit, plevel

      iunit = csym%punit
      IF (iunit >= 0) THEN
         plevel = csym%plevel
         WRITE (iunit, "(/,T2,A)") "K-point Symmetry Information"
         WRITE (iunit, '(A,T67,I14)') "       Number of Special K-points: ", csym%nkpoint
         WRITE (iunit, '(T19,A,T74,A)') " Wavevector Basis ", " Weight"
         DO i = 1, csym%nkpoint
            WRITE (iunit, '(T2,i10,3F10.5,T71,I10)') i, csym%xkpoint(1:3, i), NINT(csym%wkpoint(i))
         END DO
         WRITE (iunit, '(/,A,T63,3I6)') "       K-point Mesh: ", csym%mesh(1), csym%mesh(2), csym%mesh(3)
         WRITE (iunit, '(T19,A,T54,A)') " Wavevector Basis ", " Special Points    Rotation"
         DO i = 1, csym%mesh(1)*csym%mesh(2)*csym%mesh(3)
            WRITE (iunit, '(T2,i10,3F10.5,T45,3I12)') i, csym%kpmesh(1:3, i), &
               csym%kplink(1:2, i), csym%kpop(i)
         END DO
      END IF

   END SUBROUTINE print_kp_symmetry

END MODULE cryssym
