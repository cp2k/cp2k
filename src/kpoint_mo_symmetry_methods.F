!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utilities for transforming k-point MO coefficients under symmetry operations
! **************************************************************************************************
MODULE kpoint_mo_symmetry_methods
   USE cp_cfm_types,                    ONLY: cp_cfm_get_info,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: find_kpoint_rotation_slot,&
                                              kind_rotmat_type,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_mo_symmetry_methods'

   PUBLIC :: kpoint_same_periodic, &
             kpoint_transform_scf_mo

CONTAINS

! **************************************************************************************************
!> \brief Compare two fractional k-points modulo reciprocal lattice vectors.
!> \param xkp_a first k-point
!> \param xkp_b second k-point
!> \return true if the k-points are equivalent
! **************************************************************************************************
   LOGICAL FUNCTION kpoint_same_periodic(xkp_a, xkp_b) RESULT(same)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp_a, xkp_b

      REAL(KIND=dp), DIMENSION(3)                        :: diff

      diff(1:3) = xkp_a(1:3) - xkp_b(1:3)
      diff(1:3) = diff(1:3) - REAL(NINT(diff(1:3)), KIND=dp)
      same = SUM(ABS(diff(1:3))) < 1.0e-8_dp

   END FUNCTION kpoint_same_periodic

! **************************************************************************************************
!> \brief Transform one SCF MO coefficient matrix to an equivalent full-mesh k-point.
!> \param src source MO coefficients
!> \param dst transformed MO coefficients
!> \param qs_kpoint SCF k-point object containing symmetry operations
!> \param ikred representative k-point index
!> \param isym symmetry operation index; 0 direct copy, -1 time reversal
!> \param success true if the transformation was performed
!> \param reason diagnostic message
! **************************************************************************************************
   SUBROUTINE kpoint_transform_scf_mo(src, dst, qs_kpoint, ikred, &
                                      isym, success, reason)
      TYPE(cp_cfm_type), INTENT(IN)                      :: src, dst
      TYPE(kpoint_type), POINTER                         :: qs_kpoint
      INTEGER, INTENT(IN)                                :: ikred, isym
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason

      COMPLEX(KIND=dp)                                   :: coeff, phase
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: dst_coeff, src_coeff
      INTEGER :: iao, iatom, ikind, imo, irow, irow_source, irow_target, irow_work, nao, natom, &
         nmo, rot_slot, source_atom, source_dim, target_atom, target_dim
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_size, ao_start
      LOGICAL                                            :: reverse_phase, time_reversal
      REAL(KIND=dp)                                      :: arg, coskl, sinkl
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_phase
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rotmat
      TYPE(kind_rotmat_type), DIMENSION(:), POINTER      :: kind_rot
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      success = .FALSE.
      reason = ""
      IF (isym == 0) THEN
         CALL cp_cfm_to_cfm(src, dst)
         success = .TRUE.
         RETURN
      END IF

      CALL cp_cfm_get_info(src, nrow_global=nao, ncol_global=nmo)
      ALLOCATE (src_coeff(nao, nmo), dst_coeff(nao, nmo))
      CALL cp_cfm_get_submatrix(src, src_coeff)
      dst_coeff(:, :) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      IF (isym == -1) THEN
         dst_coeff(:, :) = CONJG(src_coeff)
         CALL cp_cfm_set_submatrix(dst, dst_coeff)
         success = .TRUE.
         RETURN
      END IF

      IF (.NOT. ASSOCIATED(qs_kpoint%kp_sym)) THEN
         reason = "SCF symmetry operation data are not available"
         RETURN
      END IF
      kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
      IF (.NOT. ASSOCIATED(kpsym)) THEN
         reason = "SCF k-point symmetry operation is not available"
         RETURN
      END IF
      IF (isym > kpsym%nwred) THEN
         reason = "SCF k-point symmetry operation index is out of range"
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(qs_kpoint%atype) .OR. .NOT. ASSOCIATED(qs_kpoint%kind_rotmat)) THEN
         reason = "SCF atom mappings or basis rotation matrices are not available"
         RETURN
      END IF

      rot_slot = find_kpoint_rotation_slot(qs_kpoint, kpsym%rotp(isym))
      IF (rot_slot == 0) THEN
         reason = "could not match the SCF symmetry operation to a basis rotation"
         RETURN
      END IF
      kind_rot => qs_kpoint%kind_rotmat(rot_slot, :)
      natom = SIZE(qs_kpoint%atype)
      ALLOCATE (ao_start(natom), ao_size(natom))
      irow = 1
      DO iatom = 1, natom
         ikind = qs_kpoint%atype(iatom)
         IF (.NOT. ASSOCIATED(kind_rot(ikind)%rmat)) THEN
            reason = "a required basis rotation matrix is not available"
            RETURN
         END IF
         ao_start(iatom) = irow
         ao_size(iatom) = SIZE(kind_rot(ikind)%rmat, 2)
         irow = irow + ao_size(iatom)
      END DO
      IF (irow - 1 /= nao) THEN
         reason = "atom-resolved AO dimensions do not match the MO coefficient matrix"
         RETURN
      END IF

      time_reversal = kpsym%rotp(isym) < 0
      reverse_phase = qs_kpoint%gamma_centered .AND. ANY(MOD(qs_kpoint%nkp_grid, 2) == 0)
      IF (ASSOCIATED(kpsym%phase_mode)) THEN
         IF (kpsym%phase_mode(isym) > 0) reverse_phase = kpsym%phase_mode(isym) == 2
      END IF
      xkp_phase(1:3) = kpsym%xkp(1:3, isym)
      DO iatom = 1, natom
         source_atom = iatom
         target_atom = kpsym%f0(iatom, isym)
         ikind = qs_kpoint%atype(source_atom)
         rotmat => kind_rot(ikind)%rmat
         source_dim = ao_size(source_atom)
         target_dim = ao_size(target_atom)
         IF (SIZE(rotmat, 1) /= target_dim .OR. SIZE(rotmat, 2) /= source_dim) THEN
            reason = "basis rotation dimensions do not match the atom/AO symmetry transform"
            RETURN
         END IF
         arg = REAL(kpsym%fcell_gauge(1, source_atom, isym), KIND=dp)*xkp_phase(1) + &
               REAL(kpsym%fcell_gauge(2, source_atom, isym), KIND=dp)*xkp_phase(2) + &
               REAL(kpsym%fcell_gauge(3, source_atom, isym), KIND=dp)*xkp_phase(3)
         IF (ASSOCIATED(kpsym%kgphase)) THEN
            arg = arg + kpsym%kgphase(source_atom, isym)
         END IF
         IF (reverse_phase) arg = -arg
         coskl = COS(twopi*arg)
         sinkl = SIN(twopi*arg)
         DO imo = 1, nmo
            DO irow_work = 1, source_dim
               irow_source = ao_start(source_atom) + irow_work - 1
               coeff = src_coeff(irow_source, imo)
               IF (time_reversal) coeff = CONJG(coeff)
               phase = CMPLX(coskl, sinkl, KIND=dp)*coeff
               DO iao = 1, target_dim
                  irow_target = ao_start(target_atom) + iao - 1
                  dst_coeff(irow_target, imo) = dst_coeff(irow_target, imo) + &
                                                rotmat(iao, irow_work)*phase
               END DO
            END DO
         END DO
      END DO

      CALL cp_cfm_set_submatrix(dst, dst_coeff)
      success = .TRUE.

   END SUBROUTINE kpoint_transform_scf_mo

END MODULE kpoint_mo_symmetry_methods
