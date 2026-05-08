!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief K-points and crystal symmetry routines
!> \author jgh
! **************************************************************************************************
MODULE cryssym

   USE bibliography,                    ONLY: Togo2018,&
                                              Worlton1972,&
                                              cite_reference
   USE kinds,                           ONLY: dp
   USE kpsym,                           ONLY: group1s,&
                                              k290s
   USE mathlib,                         ONLY: inv_3x3
   USE spglib_f08,                      ONLY: spg_get_international,&
                                              spg_get_major_version,&
                                              spg_get_micro_version,&
                                              spg_get_minor_version,&
                                              spg_get_multiplicity,&
                                              spg_get_pointgroup,&
                                              spg_get_schoenflies,&
                                              spg_get_symmetry
   USE string_utilities,                ONLY: strip_control_codes
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: csym_type, release_csym_type, print_crys_symmetry, print_kp_symmetry
   PUBLIC :: crys_sym_gen, kpoint_gen

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cryssym'

! **************************************************************************************************
!> \brief CSM type
!> \par   Content:
!>
! **************************************************************************************************
   TYPE csym_type
      LOGICAL                                     :: symlib = .FALSE.
      LOGICAL                                     :: fullgrid = .FALSE.
      LOGICAL                                     :: inversion_only = .FALSE.
      LOGICAL                                     :: spglib_reduction = .FALSE.
      LOGICAL                                     :: spglib_backend = .FALSE.
      INTEGER                                     :: plevel = 0
      INTEGER                                     :: punit = -1
      INTEGER                                     :: istriz = -1
      REAL(KIND=dp)                               :: delta = 1.0e-8_dp
      REAL(KIND=dp), DIMENSION(3, 3)              :: hmat = 0.0_dp
      ! KPOINTS
      REAL(KIND=dp), DIMENSION(3)                 :: wvk0 = 0.0_dp
      INTEGER, DIMENSION(3)                       :: mesh = 0
      INTEGER                                     :: nkpoint = 0
      INTEGER                                     :: nat = 0
      INTEGER, DIMENSION(:), ALLOCATABLE          :: atype
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: scoord
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: xkpoint
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: wkpoint
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: kpmesh
      INTEGER, DIMENSION(:, :), ALLOCATABLE       :: kplink
      INTEGER, DIMENSION(:), ALLOCATABLE          :: kpop
      !SPGLIB
      CHARACTER(len=11)                           :: international_symbol = ""
      CHARACTER(len=6)                            :: pointgroup_symbol = ""
      CHARACTER(len=10)                           :: schoenflies = ""
      INTEGER                                     :: n_operations = 0
      INTEGER, DIMENSION(:, :, :), ALLOCATABLE    :: rotations
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: translations
      !K290
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE :: rt
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: vt
      INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: f0
      INTEGER                                     :: nrtot = 0
      INTEGER, DIMENSION(:), ALLOCATABLE          :: ibrot
   END TYPE csym_type

CONTAINS

! **************************************************************************************************
!> \brief Release the CSYM type
!> \param csym  The CSYM type
! **************************************************************************************************
   SUBROUTINE release_csym_type(csym)
      TYPE(csym_type)                                    :: csym

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
      IF (ALLOCATED(csym%rt)) THEN
         DEALLOCATE (csym%rt)
      END IF
      IF (ALLOCATED(csym%vt)) THEN
         DEALLOCATE (csym%vt)
      END IF
      IF (ALLOCATED(csym%f0)) THEN
         DEALLOCATE (csym%f0)
      END IF
      IF (ALLOCATED(csym%ibrot)) THEN
         DEALLOCATE (csym%ibrot)
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'crys_sym_gen'

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

      csym%n_operations = 0

      !..try spglib
      major = spg_get_major_version()
      minor = spg_get_minor_version()
      micro = spg_get_micro_version()
      IF (major == 0) THEN
         CALL cp_warn(__LOCATION__, "Symmetry library SPGLIB not available")
         spglib = .FALSE.
      ELSE
         spglib = .TRUE.
         CALL cite_reference(Togo2018)
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

            CALL strip_control_codes(csym%international_symbol)
            CALL strip_control_codes(csym%schoenflies)
            CALL strip_control_codes(csym%pointgroup_symbol)
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
!> \param gamma_centered ...
!> \param inversion_symmetry_only ...
!> \param use_spglib_reduction ...
!> \param use_spglib_backend ...
! **************************************************************************************************
   SUBROUTINE kpoint_gen(csym, nk, symm, shift, full_grid, gamma_centered, &
                         inversion_symmetry_only, use_spglib_reduction, use_spglib_backend)
      TYPE(csym_type)                                    :: csym
      INTEGER, INTENT(IN)                                :: nk(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: symm
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: shift(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: full_grid, gamma_centered, &
                                                            inversion_symmetry_only, &
                                                            use_spglib_reduction, &
                                                            use_spglib_backend

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'kpoint_gen'

      INTEGER                                            :: handle, i, ik, j, nkp, nkpts
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kpop, xptr
      LOGICAL                                            :: fullmesh, gamma_mesh, inversion_only, &
                                                            orthorhombic_cell, spglib_backend, &
                                                            spglib_reduction
      REAL(KIND=dp)                                      :: cell_eps
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: wkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp

      CALL timeset(routineN, handle)

      IF (PRESENT(shift)) THEN
         csym%wvk0 = shift
      ELSE
         csym%wvk0 = 0.0_dp
      END IF

      csym%istriz = -1
      IF (PRESENT(symm)) THEN
         IF (symm) csym%istriz = 1
      END IF

      IF (PRESENT(full_grid)) THEN
         fullmesh = full_grid
      ELSE
         fullmesh = .FALSE.
      END IF
      csym%fullgrid = fullmesh

      IF (PRESENT(gamma_centered)) THEN
         gamma_mesh = gamma_centered
      ELSE
         gamma_mesh = .FALSE.
      END IF

      IF (PRESENT(inversion_symmetry_only)) THEN
         inversion_only = inversion_symmetry_only
      ELSE
         inversion_only = .FALSE.
      END IF
      csym%inversion_only = inversion_only

      IF (PRESENT(use_spglib_reduction)) THEN
         spglib_reduction = use_spglib_reduction
      ELSE
         spglib_reduction = .FALSE.
      END IF
      csym%spglib_reduction = spglib_reduction

      IF (PRESENT(use_spglib_backend)) THEN
         spglib_backend = use_spglib_backend
      ELSE
         spglib_backend = .FALSE.
      END IF
      csym%spglib_backend = spglib_backend

      cell_eps = 10.0_dp*csym%delta*MAX(1.0_dp, MAXVAL(ABS(csym%hmat)))
      orthorhombic_cell = ABS(csym%hmat(1, 2)) + ABS(csym%hmat(1, 3)) + &
                          ABS(csym%hmat(2, 1)) + ABS(csym%hmat(2, 3)) + &
                          ABS(csym%hmat(3, 1)) + ABS(csym%hmat(3, 2)) <= cell_eps

      IF (spglib_backend .AND. .NOT. spglib_reduction) THEN
         CALL cp_abort(__LOCATION__, &
                       "SYMMETRY_BACKEND SPGLIB requires SYMMETRY_REDUCTION_METHOD SPGLIB")
      END IF

      csym%nkpoint = 0
      csym%mesh(1:3) = nk(1:3)
      csym%nrtot = 0
      IF (ALLOCATED(csym%rt)) DEALLOCATE (csym%rt)
      IF (ALLOCATED(csym%vt)) DEALLOCATE (csym%vt)
      IF (ALLOCATED(csym%ibrot)) DEALLOCATE (csym%ibrot)
      IF (ALLOCATED(csym%f0)) DEALLOCATE (csym%f0)
      ALLOCATE (csym%rt(3, 3, 0), csym%vt(3, 0), csym%ibrot(0), csym%f0(csym%nat, 0))

      nkpts = nk(1)*nk(2)*nk(3)
      ALLOCATE (xkp(3, nkpts), wkp(nkpts), kpop(nkpts))
      ! kp: link
      ALLOCATE (csym%kplink(2, nkpts))
      csym%kplink = 0
      kpop = 0

      ! go through all the options
      IF (csym%symlib) THEN
         ! symmetry library is available
         IF (fullmesh) THEN
            ! full mesh requested
            CALL full_grid_gen(nk, xkp, wkp, shift, gamma_centered=gamma_mesh)
            IF (csym%istriz == 1) THEN
               ! use inversion symmetry
               CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
            ELSE
               ! full kpoint mesh is used
            END IF
         ELSE IF (csym%istriz /= 1 .OR. inversion_only) THEN
            ! use inversion symmetry
            CALL full_grid_gen(nk, xkp, wkp, shift, gamma_centered=gamma_mesh)
            CALL inversion_symm(xkp, wkp, csym%kplink(1, :))
         ELSE
            ! use symmetry library to reduce k-points
            IF (SUM(ABS(csym%wvk0)) /= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, "MacDonald shifted k-point meshes are only "// &
                             "possible without symmetrization.")
            END IF

            CALL full_grid_gen(nk, xkp, wkp, shift, gamma_centered=gamma_mesh)
            IF (spglib_backend) THEN
               CALL kp_symmetry_spglib(csym, xkp, wkp, kpop)
            ELSE
               CALL kp_symmetry(csym, xkp, wkp, kpop, use_spglib_reduction=spglib_reduction)
            END IF

         END IF
      ELSE
         ! no symmetry library is available
         CALL full_grid_gen(nk, xkp, wkp, shift, gamma_centered=gamma_mesh)
         IF (csym%istriz == 1 .AND. .NOT. fullmesh .AND. .NOT. inversion_only .AND. &
             orthorhombic_cell) THEN
            ! fall back to the K290 atom mapping for simple cells when SPGLIB is not linked
            IF (SUM(ABS(csym%wvk0)) /= 0.0_dp) THEN
               CALL cp_abort(__LOCATION__, "MacDonald shifted k-point meshes are only "// &
                             "possible without symmetrization.")
            END IF
            CALL kp_symmetry(csym, xkp, wkp, kpop, use_spglib_reduction=.FALSE.)
         ELSE IF (csym%istriz /= 1 .AND. fullmesh) THEN
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
         IF (wkp(i) > 0.0_dp) nkp = nkp + 1
      END DO

      ! store reduced kpoint set
      csym%nkpoint = nkp
      ALLOCATE (csym%xkpoint(3, nkp), csym%wkpoint(nkp))
      ALLOCATE (xptr(nkp))
      j = 0
      DO ik = 1, nkpts
         IF (wkp(ik) > 0.0_dp) THEN
            j = j + 1
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
      DO ik = 1, nkpts
         i = csym%kplink(1, ik)
         DO j = 1, nkp
            IF (i == xptr(j)) THEN
               csym%kplink(2, ik) = j
               EXIT
            END IF
         END DO
      END DO
      DEALLOCATE (xptr)

      ! kp: operations
      ALLOCATE (csym%kpop(nkpts))
      IF (csym%nrtot > 0 .AND. csym%istriz == 1 .AND. .NOT. fullmesh .AND. &
          .NOT. inversion_only) THEN
         ! atomic symmetry operations possible
         csym%kpop(1:nkpts) = kpop(1:nkpts)
         DO ik = 1, nkpts
            CPASSERT(csym%kpop(ik) /= 0)
         END DO
      ELSE
         ! only time reversal symmetry
         DO ik = 1, nkpts
            IF (wkp(ik) > 0.0_dp) THEN
               csym%kpop(ik) = 1
            ELSE
               csym%kpop(ik) = 2
            END IF
         END DO
      END IF

      DEALLOCATE (xkp, wkp, kpop)

      CALL timestop(handle)

   END SUBROUTINE kpoint_gen

! **************************************************************************************************
!> \brief ...
!> \param csym ...
!> \param xkp ...
!> \param wkp ...
!> \param kpop ...
!> \param use_spglib_reduction ...
! **************************************************************************************************
   SUBROUTINE kp_symmetry(csym, xkp, wkp, kpop, use_spglib_reduction)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: kpop
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_spglib_reduction

      INTEGER                                            :: i, ihc, ihg, indpg, iou, iq1, iq2, iq3, &
                                                            istriz, isy, li, nat, nc, nhash, &
                                                            nkpoint, nrot, nsp, ntvec
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: includ, isc, list, lwght, ty
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: f0, lrot
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: srot
      INTEGER, DIMENSION(48)                             :: ib
      LOGICAL                                            :: spglib_reduction
      REAL(KIND=dp)                                      :: alat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rlist, rx, tvec, wvkl, xkapa
      REAL(KIND=dp), DIMENSION(3)                        :: a1, a2, a3, b1, b2, b3, origin, wvk0
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat, strain
      REAL(KIND=dp), DIMENSION(3, 3, 48)                 :: r
      REAL(KIND=dp), DIMENSION(3, 48)                    :: vt

      iou = csym%punit
      hmat = csym%hmat
      nat = csym%nat
      iq1 = csym%mesh(1)
      iq2 = csym%mesh(2)
      iq3 = csym%mesh(3)
      nkpoint = 10*iq1*iq2*iq3
      wvk0 = csym%wvk0
      istriz = csym%istriz
      IF (PRESENT(use_spglib_reduction)) THEN
         spglib_reduction = use_spglib_reduction
      ELSE
         spglib_reduction = .FALSE.
      END IF
      a1(1:3) = hmat(1:3, 1)
      a2(1:3) = hmat(1:3, 2)
      a3(1:3) = hmat(1:3, 3)
      alat = hmat(1, 1)
      strain = 0.0_dp
      ALLOCATE (xkapa(3, nat), rx(3, nat), tvec(3, 200), ty(nat), isc(nat), f0(49, nat))
      ty(1:nat) = csym%atype(1:nat)
      nsp = MAXVAL(ty)
      DO i = 1, nat
         xkapa(1:3, i) = MATMUL(hmat, csym%scoord(1:3, i))
      END DO
      nhash = 1000
      ALLOCATE (wvkl(3, nkpoint), rlist(3, nkpoint), includ(nkpoint), list(nhash + nkpoint))
      ALLOCATE (lrot(48, nkpoint), lwght(nkpoint))

      IF (iou > 0) THEN
         WRITE (iou, '(/,(T2,A79))') &
            "*******************************************************************************", &
            "**                      Special K-Point Generation by K290                   **", &
            "*******************************************************************************"
      END IF
      CALL cite_reference(Worlton1972)
      IF (spglib_reduction) CALL cite_reference(Togo2018)

      CALL K290s(iou, nat, nkpoint, nsp, iq1, iq2, iq3, istriz, &
                 a1, a2, a3, alat, strain, xkapa, rx, tvec, &
                 ty, isc, f0, ntvec, wvk0, wvkl, lwght, lrot, &
                 nhash, includ, list, rlist, csym%delta)

      CALL GROUP1s(0, a1, a2, a3, nat, ty, xkapa, b1, b2, b3, &
                   ihg, ihc, isy, li, nc, indpg, ib, ntvec, &
                   vt, f0, r, tvec, origin, rx, isc, csym%delta)

      IF (iou > 0) THEN
         WRITE (iou, '((T2,A79))') &
            "*******************************************************************************", &
            "**                              Finished K290                                **", &
            "*******************************************************************************"
      END IF

      csym%nrtot = nc
      IF (ALLOCATED(csym%rt)) DEALLOCATE (csym%rt)
      IF (ALLOCATED(csym%vt)) DEALLOCATE (csym%vt)
      IF (ALLOCATED(csym%ibrot)) DEALLOCATE (csym%ibrot)
      IF (ALLOCATED(csym%f0)) DEALLOCATE (csym%f0)
      ALLOCATE (csym%rt(3, 3, nc), csym%vt(3, nc), csym%ibrot(nc))
      csym%vt(1:3, 1:nc) = vt(1:3, 1:nc)
      ALLOCATE (csym%f0(nat, nc))
      DO i = 1, nc
         csym%rt(1:3, 1:3, i) = r(1:3, 1:3, ib(i))
         csym%f0(1:nat, i) = f0(i, 1:nat)
      END DO
      csym%ibrot(1:nc) = ib(1:nc)

      IF (csym%n_operations > nc) THEN
         IF (ALLOCATED(srot)) DEALLOCATE (srot)
         ALLOCATE (srot(3, 3, csym%n_operations))
         CALL setup_spglib_operations(csym, srot, nrot)
         CALL reduce_spglib_kpoint_mesh(csym, xkp, wkp, kpop, srot, nrot)
         DEALLOCATE (srot)
      ELSE IF (spglib_reduction) THEN
         ALLOCATE (srot(3, 3, csym%n_operations))
         CALL setup_spglib_reduction_rotations(csym, srot, nrot)
         CALL reduce_spglib_kpoint_mesh_k290(csym, xkp, wkp, kpop, srot, nrot, &
                                             a1, a2, a3, b1, b2, b3, alat)
         DEALLOCATE (srot)
      ELSE
         CALL reduce_kpoint_mesh(csym, xkp, wkp, kpop, nc, ib, r, a1, a2, a3, b1, b2, b3, alat)
      END IF
      DEALLOCATE (xkapa, rx, tvec, ty, isc, f0)
      DEALLOCATE (wvkl, rlist, includ, list)
      DEALLOCATE (lrot, lwght)

   END SUBROUTINE kp_symmetry

! **************************************************************************************************
!> \brief Reduce a CP2K Monkhorst-Pack mesh using SPGLIB symmetry operations
!> \param csym ...
!> \param xkp ...
!> \param wkp ...
!> \param kpop ...
! **************************************************************************************************
   SUBROUTINE kp_symmetry_spglib(csym, xkp, wkp, kpop)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: kpop

      INTEGER                                            :: iou, nrot
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: srot

      iou = csym%punit
      IF (iou > 0) THEN
         WRITE (iou, '(/,(T2,A79))') &
            "*******************************************************************************", &
            "**                     Special K-Point Generation by SPGLIB                 **", &
            "*******************************************************************************"
      END IF
      CALL cite_reference(Togo2018)

      ALLOCATE (srot(3, 3, csym%n_operations))
      CALL setup_spglib_operations(csym, srot, nrot)
      CALL reduce_spglib_kpoint_mesh(csym, xkp, wkp, kpop, srot, nrot)
      DEALLOCATE (srot)

      IF (iou > 0) THEN
         WRITE (iou, '((T2,A79))') &
            "*******************************************************************************", &
            "**                              Finished SPGLIB                             **", &
            "*******************************************************************************"
      END IF

   END SUBROUTINE kp_symmetry_spglib

! **************************************************************************************************
!> \brief Store usable SPGLIB space-group operations for k-point symmetry
!> \param csym ...
!> \param srot integer rotations in fractional coordinates
!> \param nrot number of stored rotations
! **************************************************************************************************
   SUBROUTINE setup_spglib_operations(csym, srot, nrot)
      TYPE(csym_type)                                    :: csym
      INTEGER, DIMENSION(:, :, :), INTENT(OUT)           :: srot
      INTEGER, INTENT(OUT)                               :: nrot

      INTEGER                                            :: iop, jop, pass
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: perm
      INTEGER, DIMENSION(3, 3)                           :: eye, irot
      LOGICAL                                            :: duplicate, identity, valid, &
                                                            zero_translation
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_inv, rfrac

      CPASSERT(csym%symlib)

      srot = 0
      csym%nrtot = 0
      IF (ALLOCATED(csym%rt)) DEALLOCATE (csym%rt)
      IF (ALLOCATED(csym%vt)) DEALLOCATE (csym%vt)
      IF (ALLOCATED(csym%ibrot)) DEALLOCATE (csym%ibrot)
      IF (ALLOCATED(csym%f0)) DEALLOCATE (csym%f0)
      ALLOCATE (csym%rt(3, 3, csym%n_operations), csym%vt(3, csym%n_operations))
      ALLOCATE (csym%ibrot(csym%n_operations), csym%f0(csym%nat, csym%n_operations))
      csym%rt = 0.0_dp
      csym%vt = 0.0_dp
      csym%ibrot = 0
      csym%f0 = 0

      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)
      h_inv = inv_3x3(csym%hmat)
      ALLOCATE (perm(csym%nat))

      eye = 0
      eye(1, 1) = 1
      eye(2, 2) = 1
      eye(3, 3) = 1

      nrot = 0
      ! Operation 1 is used as the untransformed representative k-point.
      ! Prefer integer translations before fractional alternatives with the same rotation.
      DO pass = 1, 4
         DO iop = 1, csym%n_operations
            irot(1:3, 1:3) = csym%rotations(1:3, 1:3, iop)
            identity = ALL(irot == eye)
            zero_translation = ALL(ABS(csym%translations(1:3, iop) - &
                                       ANINT(csym%translations(1:3, iop))) <= eps)
            IF (pass == 1 .AND. (.NOT. identity .OR. .NOT. zero_translation)) CYCLE
            IF (pass == 2 .AND. (identity .OR. .NOT. zero_translation)) CYCLE
            IF (pass == 3 .AND. (.NOT. identity .OR. zero_translation)) CYCLE
            IF (pass == 4 .AND. (identity .OR. zero_translation)) CYCLE

            duplicate = .FALSE.
            DO jop = 1, nrot
               IF (ALL(irot == srot(:, :, jop)) .AND. &
                   ALL(ABS(csym%translations(1:3, iop) - csym%vt(1:3, jop) - &
                           ANINT(csym%translations(1:3, iop) - csym%vt(1:3, jop))) <= eps)) THEN
                  duplicate = .TRUE.
                  EXIT
               END IF
            END DO
            IF (duplicate) CYCLE

            CALL spglib_atom_permutation(csym, irot, csym%translations(:, iop), perm, valid)
            IF (.NOT. valid) CYCLE

            nrot = nrot + 1

            srot(1:3, 1:3, nrot) = irot(1:3, 1:3)
            rfrac(1:3, 1:3) = REAL(irot(1:3, 1:3), KIND=dp)
            csym%rt(1:3, 1:3, nrot) = MATMUL(csym%hmat, MATMUL(rfrac, h_inv))
            csym%vt(1:3, nrot) = csym%translations(1:3, iop)
            csym%ibrot(nrot) = nrot
            csym%f0(1:csym%nat, nrot) = perm(1:csym%nat)
         END DO
      END DO

      DEALLOCATE (perm)
      csym%nrtot = nrot
      IF (nrot == 0) CALL cp_abort(__LOCATION__, "SPGLIB did not return usable symmetry operations")

   END SUBROUTINE setup_spglib_operations

! **************************************************************************************************
!> \brief Store unique SPGLIB rotations for K290-backend diagnostic reduction
!> \param csym ...
!> \param srot integer rotations in fractional coordinates
!> \param nrot number of stored rotations
! **************************************************************************************************
   SUBROUTINE setup_spglib_reduction_rotations(csym, srot, nrot)
      TYPE(csym_type)                                    :: csym
      INTEGER, DIMENSION(:, :, :), INTENT(OUT)           :: srot
      INTEGER, INTENT(OUT)                               :: nrot

      INTEGER                                            :: iop, jop, pass
      INTEGER, DIMENSION(3, 3)                           :: eye, irot
      LOGICAL                                            :: duplicate, identity

      CPASSERT(csym%symlib)

      srot = 0
      eye = 0
      eye(1, 1) = 1
      eye(2, 2) = 1
      eye(3, 3) = 1

      nrot = 0
      ! Keep the identity first, matching the representative k-point operation.
      DO pass = 1, 2
         DO iop = 1, csym%n_operations
            irot(1:3, 1:3) = csym%rotations(1:3, 1:3, iop)
            identity = ALL(irot == eye)
            IF (pass == 1 .AND. .NOT. identity) CYCLE
            IF (pass == 2 .AND. identity) CYCLE

            duplicate = .FALSE.
            DO jop = 1, nrot
               IF (ALL(irot == srot(:, :, jop))) THEN
                  duplicate = .TRUE.
                  EXIT
               END IF
            END DO
            IF (duplicate) CYCLE

            nrot = nrot + 1
            srot(1:3, 1:3, nrot) = irot(1:3, 1:3)
         END DO
      END DO

      IF (nrot == 0) CALL cp_abort(__LOCATION__, "SPGLIB did not return usable symmetry rotations")

   END SUBROUTINE setup_spglib_reduction_rotations

! **************************************************************************************************
!> \brief Determine the atom permutation generated by a SPGLIB space-group operation
!> \param csym ...
!> \param rot integer rotation in fractional coordinates
!> \param trans fractional translation
!> \param perm atom permutation
!> \param valid whether all atoms were mapped
! **************************************************************************************************
   SUBROUTINE spglib_atom_permutation(csym, rot, trans, perm, valid)
      TYPE(csym_type)                                    :: csym
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: rot
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: trans
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: perm
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: i, j, nat
      LOGICAL                                            :: found
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: used
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), DIMENSION(3)                        :: diff, spos
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rfrac

      nat = csym%nat
      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)
      rfrac(1:3, 1:3) = REAL(rot(1:3, 1:3), KIND=dp)
      ALLOCATE (used(nat))
      used = .FALSE.
      perm = 0
      valid = .TRUE.

      DO i = 1, nat
         spos(1:3) = MATMUL(rfrac(1:3, 1:3), csym%scoord(1:3, i)) + trans(1:3)
         found = .FALSE.
         DO j = 1, nat
            IF (used(j)) CYCLE
            IF (csym%atype(i) /= csym%atype(j)) CYCLE
            diff(1:3) = spos(1:3) - csym%scoord(1:3, j)
            diff(1:3) = diff(1:3) - ANINT(diff(1:3))
            IF (ALL(ABS(diff(1:3)) < eps)) THEN
               perm(i) = j
               used(j) = .TRUE.
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            valid = .FALSE.
            EXIT
         END IF
      END DO

      DEALLOCATE (used)

   END SUBROUTINE spglib_atom_permutation

! **************************************************************************************************
!> \brief Reduce a k-point mesh with SPGLIB direct-space operations
!> \param csym ...
!> \param xkp full k-point mesh in reciprocal lattice coordinates
!> \param wkp reduced k-point weights
!> \param kpop symmetry operation mapping the representative k-point to a mesh point
!> \param srot integer rotations in fractional coordinates
!> \param nrot number of stored rotations
! **************************************************************************************************
   SUBROUTINE reduce_spglib_kpoint_mesh(csym, xkp, wkp, kpop, srot, nrot)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: kpop
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: srot
      INTEGER, INTENT(IN)                                :: nrot

      INTEGER                                            :: i, iop, isign, j, kr, nkpts, score
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kscore
      INTEGER, DIMENSION(3, 3)                           :: krot
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), DIMENSION(3)                        :: diff, rr

      nkpts = SIZE(wkp)
      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)
      ALLOCATE (kscore(nkpts))

      wkp = 0.0_dp
      kpop = 0
      csym%kplink(1, :) = 0
      kscore = HUGE(0)

      DO i = 1, nkpts
         IF (csym%kplink(1, i) /= 0) CYCLE

         csym%kplink(1, i) = i
         wkp(i) = 1.0_dp
         kpop(i) = 1
         kscore(i) = 0

         DO iop = 1, nrot
            kr = csym%ibrot(iop)
            krot = reciprocal_rotation(srot(:, :, kr))
            score = spglib_operation_score(csym, iop, srot(:, :, kr))
            DO isign = 1, 2
               rr(1:3) = MATMUL(REAL(krot(1:3, 1:3), KIND=dp), xkp(1:3, i))
               IF (isign == 2) THEN
                  rr(1:3) = -rr(1:3)
                  kr = -csym%ibrot(iop)
               ELSE
                  kr = csym%ibrot(iop)
               END IF

               DO j = 1, nkpts
                  diff(1:3) = xkp(1:3, j) - rr(1:3)
                  diff(1:3) = diff(1:3) - ANINT(diff(1:3))
                  IF (ALL(ABS(diff(1:3)) < eps)) THEN
                     IF (csym%kplink(1, j) == 0) THEN
                        csym%kplink(1, j) = i
                        wkp(i) = wkp(i) + 1.0_dp
                        kpop(j) = kr
                        kscore(j) = score
                     ELSE
                        CPASSERT(csym%kplink(1, j) == i)
                        IF (score < kscore(j)) THEN
                           kpop(j) = kr
                           kscore(j) = score
                        END IF
                     END IF
                     EXIT
                  END IF
               END DO
               IF (j > nkpts) CYCLE
            END DO
         END DO
      END DO

      DO i = 1, nkpts
         CPASSERT(csym%kplink(1, i) /= 0)
         CPASSERT(kpop(i) /= 0)
      END DO
      DEALLOCATE (kscore)

   END SUBROUTINE reduce_spglib_kpoint_mesh

! **************************************************************************************************
!> \brief Reduce a k-point mesh with SPGLIB rotations and K290 operations
!> \param csym ...
!> \param xkp full k-point mesh in reciprocal lattice coordinates
!> \param wkp reduced k-point weights
!> \param kpop K290 operation mapping the representative k-point to a mesh point
!> \param srot SPGLIB integer rotations in fractional coordinates
!> \param nrot number of stored rotations
!> \param a1 first lattice vector
!> \param a2 second lattice vector
!> \param a3 third lattice vector
!> \param b1 first reciprocal lattice vector
!> \param b2 second reciprocal lattice vector
!> \param b3 third reciprocal lattice vector
!> \param alat lattice scaling used by K290
! **************************************************************************************************
   SUBROUTINE reduce_spglib_kpoint_mesh_k290(csym, xkp, wkp, kpop, srot, nrot, &
                                             a1, a2, a3, b1, b2, b3, alat)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: kpop
      INTEGER, DIMENSION(:, :, :), INTENT(IN)            :: srot
      INTEGER, INTENT(IN)                                :: nrot
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: a1, a2, a3, b1, b2, b3
      REAL(KIND=dp), INTENT(IN)                          :: alat

      INTEGER                                            :: i, iop, isign, j, k290_op, nkpts
      INTEGER, DIMENSION(3, 3)                           :: krot
      LOGICAL                                            :: valid
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), DIMENSION(3)                        :: diff, rr

      nkpts = SIZE(wkp)
      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)

      wkp = 0.0_dp
      kpop = 0
      csym%kplink(1, :) = 0

      DO i = 1, nkpts
         IF (csym%kplink(1, i) /= 0) CYCLE

         csym%kplink(1, i) = i
         wkp(i) = 1.0_dp
         kpop(i) = 1

         DO iop = 1, nrot
            krot = reciprocal_rotation(srot(:, :, iop))
            DO isign = 1, 2
               rr(1:3) = MATMUL(REAL(krot(1:3, 1:3), KIND=dp), xkp(1:3, i))
               IF (isign == 2) rr(1:3) = -rr(1:3)

               DO j = 1, nkpts
                  diff(1:3) = xkp(1:3, j) - rr(1:3)
                  diff(1:3) = diff(1:3) - ANINT(diff(1:3))
                  IF (ALL(ABS(diff(1:3)) < eps)) THEN
                     IF (j == i) EXIT
                     IF (csym%kplink(1, j) /= 0) THEN
                        CPASSERT(csym%kplink(1, j) == i)
                        EXIT
                     END IF

                     CALL find_k290_kpoint_operation(csym, xkp(1:3, i), xkp(1:3, j), &
                                                     a1, a2, a3, b1, b2, b3, alat, &
                                                     k290_op, valid)
                     IF (.NOT. valid) THEN
                        CALL cp_abort(__LOCATION__, &
                                      "SPGLIB k-point mapping is not represented by K290 backend")
                     END IF
                     csym%kplink(1, j) = i
                     wkp(i) = wkp(i) + 1.0_dp
                     kpop(j) = k290_op
                     EXIT
                  END IF
               END DO
               IF (j > nkpts) CYCLE
            END DO
         END DO
      END DO

      DO i = 1, nkpts
         CPASSERT(csym%kplink(1, i) /= 0)
         CPASSERT(kpop(i) /= 0)
      END DO

   END SUBROUTINE reduce_spglib_kpoint_mesh_k290

! **************************************************************************************************
!> \brief Find a K290 operation that maps one fractional k-point to another
!> \param csym ...
!> \param xref representative k-point
!> \param xtarget target k-point
!> \param a1 first lattice vector
!> \param a2 second lattice vector
!> \param a3 third lattice vector
!> \param b1 first reciprocal lattice vector
!> \param b2 second reciprocal lattice vector
!> \param b3 third reciprocal lattice vector
!> \param alat lattice scaling used by K290
!> \param k290_op K290 operation identifier
!> \param valid whether a matching K290 operation was found
! **************************************************************************************************
   SUBROUTINE find_k290_kpoint_operation(csym, xref, xtarget, a1, a2, a3, b1, b2, b3, alat, &
                                         k290_op, valid)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xref, xtarget, a1, a2, a3, b1, b2, b3
      REAL(KIND=dp), INTENT(IN)                          :: alat
      INTEGER, INTENT(OUT)                               :: k290_op
      LOGICAL, INTENT(OUT)                               :: valid

      INTEGER                                            :: ibsign, iop, kr
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), DIMENSION(3)                        :: diff, rr, wcart

      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)
      k290_op = 0
      valid = .FALSE.

      DO iop = 1, csym%nrtot
         IF (iop > SIZE(csym%rt, 3)) CYCLE
         IF (csym%ibrot(iop) == 0) CYCLE
         DO ibsign = 1, 2
            wcart(1:3) = alat*(xref(1)*b1(1:3) + xref(2)*b2(1:3) + xref(3)*b3(1:3))
            wcart(1:3) = kp_apply_operation(wcart(1:3), csym%rt(1:3, 1:3, iop))
            IF (ibsign == 2) THEN
               wcart(1:3) = -wcart(1:3)
               kr = -csym%ibrot(iop)
            ELSE
               kr = csym%ibrot(iop)
            END IF
            rr(1) = DOT_PRODUCT(a1(1:3), wcart(1:3))/alat
            rr(2) = DOT_PRODUCT(a2(1:3), wcart(1:3))/alat
            rr(3) = DOT_PRODUCT(a3(1:3), wcart(1:3))/alat

            diff(1:3) = xtarget(1:3) - rr(1:3)
            diff(1:3) = diff(1:3) - ANINT(diff(1:3))
            IF (ALL(ABS(diff(1:3)) < eps)) THEN
               k290_op = kr
               valid = .TRUE.
               RETURN
            END IF
         END DO
      END DO

   END SUBROUTINE find_k290_kpoint_operation

! **************************************************************************************************
!> \brief Score SPGLIB operations to choose stable atom transformations
!> \param csym ...
!> \param iop operation index
!> \param srot integer rotation in fractional coordinates
!> \return score, lower values are preferred
! **************************************************************************************************
   FUNCTION spglib_operation_score(csym, iop, srot) RESULT(score)
      TYPE(csym_type), INTENT(IN)                        :: csym
      INTEGER, INTENT(IN)                                :: iop
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: srot
      INTEGER                                            :: score

      INTEGER                                            :: i, nat
      INTEGER, DIMENSION(3, 3)                           :: eye, r2
      REAL(KIND=dp)                                      :: eps

      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)
      nat = SIZE(csym%f0, 1)
      score = 0
      DO i = 1, nat
         IF (csym%f0(i, iop) /= i) score = score + 100
      END DO
      IF (ANY(ABS(csym%vt(1:3, iop) - ANINT(csym%vt(1:3, iop))) > eps)) score = score + 10

      eye = 0
      eye(1, 1) = 1
      eye(2, 2) = 1
      eye(3, 3) = 1
      r2(1:3, 1:3) = MATMUL(srot(1:3, 1:3), srot(1:3, 1:3))
      IF (ANY(r2(1:3, 1:3) /= eye(1:3, 1:3))) score = score + 1

   END FUNCTION spglib_operation_score

! **************************************************************************************************
!> \brief Reciprocal-space rotation corresponding to a fractional direct-space rotation
!> \param rot direct-space rotation
!> \return reciprocal-space rotation
! **************************************************************************************************
   FUNCTION reciprocal_rotation(rot) RESULT(krot)
      INTEGER, DIMENSION(3, 3), INTENT(IN)               :: rot
      INTEGER, DIMENSION(3, 3)                           :: krot

      REAL(KIND=dp), DIMENSION(3, 3)                     :: rinv

      rinv = inv_3x3(REAL(rot(1:3, 1:3), KIND=dp))
      krot(1:3, 1:3) = NINT(TRANSPOSE(rinv(1:3, 1:3)))

   END FUNCTION reciprocal_rotation

! **************************************************************************************************
!> \brief Reduce a CP2K Monkhorst-Pack mesh using K290 symmetry operations
!> \param csym ...
!> \param xkp full k-point mesh in reciprocal lattice coordinates
!> \param wkp reduced k-point weights
!> \param kpop symmetry operation mapping the representative k-point to a mesh point
!> \param nc number of point group operations
!> \param ib K290 operation identifiers
!> \param r K290 rotation matrices
!> \param a1 first lattice vector
!> \param a2 second lattice vector
!> \param a3 third lattice vector
!> \param b1 first reciprocal lattice vector
!> \param b2 second reciprocal lattice vector
!> \param b3 third reciprocal lattice vector
!> \param alat lattice scaling used by K290
! **************************************************************************************************
   SUBROUTINE reduce_kpoint_mesh(csym, xkp, wkp, kpop, nc, ib, r, a1, a2, a3, b1, b2, b3, alat)
      TYPE(csym_type)                                    :: csym
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      INTEGER, DIMENSION(:)                              :: kpop
      INTEGER, INTENT(IN)                                :: nc
      INTEGER, DIMENSION(48), INTENT(IN)                 :: ib
      REAL(KIND=dp), DIMENSION(3, 3, 48), INTENT(IN)     :: r
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: a1, a2, a3, b1, b2, b3
      REAL(KIND=dp), INTENT(IN)                          :: alat

      INTEGER                                            :: i, ibsign, iop, j, kr, nkpts
      REAL(KIND=dp)                                      :: eps
      REAL(KIND=dp), DIMENSION(3)                        :: diff, rr, wcart

      nkpts = SIZE(wkp)
      eps = MAX(1.e-12_dp, 10.0_dp*csym%delta)

      wkp = 0.0_dp
      kpop = 0
      csym%kplink(1, :) = 0

      DO i = 1, nkpts
         IF (csym%kplink(1, i) /= 0) CYCLE

         csym%kplink(1, i) = i
         wkp(i) = 1.0_dp
         kpop(i) = 1

         DO iop = 1, nc
            DO ibsign = 1, 2
               kr = ib(iop)
               wcart(1:3) = alat*(xkp(1, i)*b1(1:3) + xkp(2, i)*b2(1:3) + xkp(3, i)*b3(1:3))
               wcart(1:3) = kp_apply_operation(wcart(1:3), r(1:3, 1:3, kr))
               IF (ibsign == 2) THEN
                  wcart(1:3) = -wcart(1:3)
                  kr = -kr
               END IF
               rr(1) = DOT_PRODUCT(a1(1:3), wcart(1:3))/alat
               rr(2) = DOT_PRODUCT(a2(1:3), wcart(1:3))/alat
               rr(3) = DOT_PRODUCT(a3(1:3), wcart(1:3))/alat

               DO j = 1, nkpts
                  diff(1:3) = xkp(1:3, j) - rr(1:3)
                  diff(1:3) = diff(1:3) - ANINT(diff(1:3))
                  IF (ALL(ABS(diff(1:3)) < eps)) THEN
                     IF (csym%kplink(1, j) == 0) THEN
                        csym%kplink(1, j) = i
                        wkp(i) = wkp(i) + 1.0_dp
                        kpop(j) = kr
                     ELSE
                        CPASSERT(csym%kplink(1, j) == i)
                     END IF
                     EXIT
                  END IF
               END DO
               ! Some point-group operations are incompatible with the requested Monkhorst-Pack mesh.
               IF (j > nkpts) CYCLE
            END DO
         END DO
      END DO

      DO i = 1, nkpts
         CPASSERT(csym%kplink(1, i) /= 0)
         CPASSERT(kpop(i) /= 0)
      END DO

   END SUBROUTINE reduce_kpoint_mesh

!> \brief ...
!> \param nk ...
!> \param xkp ...
!> \param wkp ...
!> \param shift ...
!> \param gamma_centered ...
! **************************************************************************************************
   SUBROUTINE full_grid_gen(nk, xkp, wkp, shift, gamma_centered)
      INTEGER, INTENT(IN)                                :: nk(3)
      REAL(KIND=dp), DIMENSION(:, :)                     :: xkp
      REAL(KIND=dp), DIMENSION(:)                        :: wkp
      REAL(KIND=dp), INTENT(IN)                          :: shift(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: gamma_centered

      INTEGER                                            :: i, idim, ix, iy, iz
      INTEGER, DIMENSION(3)                              :: ik
      LOGICAL                                            :: gamma_mesh
      REAL(KIND=dp)                                      :: kpt_latt(3)

      IF (PRESENT(gamma_centered)) THEN
         gamma_mesh = gamma_centered
      ELSE
         gamma_mesh = .FALSE.
      END IF

      wkp = 0.0_dp
      i = 0
      DO ix = 1, nk(1)
         DO iy = 1, nk(2)
            DO iz = 1, nk(3)
               i = i + 1
               ik(1) = ix
               ik(2) = iy
               ik(3) = iz
               DO idim = 1, 3
                  IF (gamma_mesh .AND. MOD(nk(idim), 2) == 0) THEN
                     kpt_latt(idim) = REAL(2*ik(idim) - nk(idim), KIND=dp)/ &
                                      (2._dp*REAL(nk(idim), KIND=dp))
                  ELSE
                     kpt_latt(idim) = REAL(2*ik(idim) - nk(idim) - 1, KIND=dp)/ &
                                      (2._dp*REAL(nk(idim), KIND=dp))
                  END IF
               END DO
               xkp(1:3, i) = kpt_latt(1:3)
               wkp(i) = 1.0_dp
            END DO
         END DO
      END DO
      DO i = 1, nk(1)*nk(2)*nk(3)
         xkp(1:3, i) = xkp(1:3, i) + shift(1:3)
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
      REAL(KIND=dp), DIMENSION(3)                        :: diff

      nkpts = SIZE(wkp, 1)

      link(:) = 0
      DO i = 1, nkpts
         IF (link(i) == 0) link(i) = i
         DO j = i + 1, nkpts
            IF (wkp(j) == 0) CYCLE
            diff(1:3) = xkp(1:3, i) + xkp(1:3, j)
            diff(1:3) = diff(1:3) - ANINT(diff(1:3))
            IF (ALL(ABS(diff(1:3)) < 1.e-12_dp)) THEN
               wkp(i) = wkp(i) + wkp(j)
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
!> \return ...
! **************************************************************************************************
   FUNCTION kp_apply_operation(x, r) RESULT(y)
      REAL(KIND=dp), INTENT(IN)                          :: x(3), r(3, 3)
      REAL(KIND=dp)                                      :: y(3)

      y(1) = r(1, 1)*x(1) + r(1, 2)*x(2) + r(1, 3)*x(3)
      y(2) = r(2, 1)*x(1) + r(2, 2)*x(2) + r(2, 3)*x(3)
      y(3) = r(3, 1)*x(1) + r(3, 2)*x(2) + r(3, 3)*x(3)

   END FUNCTION kp_apply_operation

! **************************************************************************************************
!> \brief ...
!> \param csym ...
! **************************************************************************************************
   SUBROUTINE print_crys_symmetry(csym)
      TYPE(csym_type)                                    :: csym

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

      INTEGER                                            :: i, iunit, nat, plevel

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
         IF (csym%nrtot > 0) THEN
            WRITE (iunit, '(/,A)') "       Atom Transformation Table"
            nat = SIZE(csym%f0, 1)
            DO i = 1, csym%nrtot
               WRITE (iunit, '(T10,A,I5,(T21,12I5))') " Rot=", csym%ibrot(i), csym%f0(1:nat, i)
            END DO
         END IF
      END IF

   END SUBROUTINE print_kp_symmetry

END MODULE cryssym
