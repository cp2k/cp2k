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
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kind_rotmat_type,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
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
!> \param src_real real part of source MO coefficients
!> \param src_imag imaginary part of source MO coefficients
!> \param dst_real real part of transformed MO coefficients
!> \param dst_imag imaginary part of transformed MO coefficients
!> \param qs_kpoint SCF k-point object containing symmetry operations
!> \param ikred representative k-point index
!> \param isym symmetry operation index; 0 direct copy, -1 time reversal
!> \param para_env global parallel environment
!> \param success true if the transformation was performed
!> \param reason diagnostic message
! **************************************************************************************************
   SUBROUTINE kpoint_transform_scf_mo(src_real, src_imag, dst_real, dst_imag, qs_kpoint, ikred, &
                                      isym, para_env, success, reason)
      TYPE(cp_fm_type), INTENT(IN)                       :: src_real, src_imag, dst_real, dst_imag
      TYPE(kpoint_type), POINTER                         :: qs_kpoint
      INTEGER, INTENT(IN)                                :: ikred, isym
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(OUT)                               :: success
      CHARACTER(LEN=*), INTENT(OUT)                      :: reason

      INTEGER :: iao, iatom, ikind, imo, irow, irow_source, irow_target, irow_work, nao, natom, &
         nmo, rot_slot, source_atom, source_dim, target_atom, target_dim
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_size, ao_start
      LOGICAL                                            :: reverse_phase, time_reversal
      REAL(KIND=dp)                                      :: arg, coeff_imag, coeff_real, coskl, &
                                                            phase_imag, phase_real, sinkl
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: dst_i, dst_r, src_i, src_r
      REAL(KIND=dp), DIMENSION(3)                        :: xkp_phase
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rotmat
      TYPE(kind_rotmat_type), DIMENSION(:), POINTER      :: kind_rot
      TYPE(kpoint_sym_type), POINTER                     :: kpsym

      success = .FALSE.
      reason = ""
      IF (isym == 0) THEN
         CALL cp_fm_copy_general(src_real, dst_real, para_env)
         CALL cp_fm_copy_general(src_imag, dst_imag, para_env)
         success = .TRUE.
         RETURN
      END IF

      CALL cp_fm_get_info(src_real, nrow_global=nao, ncol_global=nmo)
      ALLOCATE (src_r(nao, nmo), src_i(nao, nmo), dst_r(nao, nmo), dst_i(nao, nmo))
      CALL cp_fm_get_submatrix(src_real, src_r)
      CALL cp_fm_get_submatrix(src_imag, src_i)
      dst_r(:, :) = 0.0_dp
      dst_i(:, :) = 0.0_dp

      IF (isym == -1) THEN
         dst_r(1:nao, 1:nmo) = src_r(1:nao, 1:nmo)
         dst_i(1:nao, 1:nmo) = -src_i(1:nao, 1:nmo)
         CALL cp_fm_set_submatrix(dst_real, dst_r)
         CALL cp_fm_set_submatrix(dst_imag, dst_i)
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         success = .TRUE.
         RETURN
      END IF

      IF (.NOT. ASSOCIATED(qs_kpoint%kp_sym)) THEN
         reason = "SCF symmetry operation data are not available"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      kpsym => qs_kpoint%kp_sym(ikred)%kpoint_sym
      IF (.NOT. ASSOCIATED(kpsym)) THEN
         reason = "SCF k-point symmetry operation is not available"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      IF (isym > kpsym%nwred) THEN
         reason = "SCF k-point symmetry operation index is out of range"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF
      IF (.NOT. ASSOCIATED(qs_kpoint%atype) .OR. .NOT. ASSOCIATED(qs_kpoint%kind_rotmat)) THEN
         reason = "SCF atom mappings or basis rotation matrices are not available"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
         RETURN
      END IF

      rot_slot = find_kpoint_rotation_slot(qs_kpoint, kpsym%rotp(isym))
      IF (rot_slot == 0) THEN
         reason = "could not match the SCF symmetry operation to a basis rotation"
         DEALLOCATE (src_r, src_i, dst_r, dst_i)
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
            DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
            RETURN
         END IF
         ao_start(iatom) = irow
         ao_size(iatom) = SIZE(kind_rot(ikind)%rmat, 2)
         irow = irow + ao_size(iatom)
      END DO
      IF (irow - 1 /= nao) THEN
         reason = "atom-resolved AO dimensions do not match the MO coefficient matrix"
         DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
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
            DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
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
               coeff_real = src_r(irow_source, imo)
               coeff_imag = src_i(irow_source, imo)
               IF (time_reversal) coeff_imag = -coeff_imag
               phase_real = coskl*coeff_real - sinkl*coeff_imag
               phase_imag = sinkl*coeff_real + coskl*coeff_imag
               DO iao = 1, target_dim
                  irow_target = ao_start(target_atom) + iao - 1
                  dst_r(irow_target, imo) = dst_r(irow_target, imo) + &
                                            rotmat(iao, irow_work)*phase_real
                  dst_i(irow_target, imo) = dst_i(irow_target, imo) + &
                                            rotmat(iao, irow_work)*phase_imag
               END DO
            END DO
         END DO
      END DO

      CALL cp_fm_set_submatrix(dst_real, dst_r)
      CALL cp_fm_set_submatrix(dst_imag, dst_i)
      DEALLOCATE (src_r, src_i, dst_r, dst_i, ao_start, ao_size)
      success = .TRUE.

   END SUBROUTINE kpoint_transform_scf_mo

! **************************************************************************************************
!> \brief Locate the basis-rotation slot corresponding to a signed k-point symmetry operation.
!> \param qs_kpoint SCF k-point object
!> \param rotp signed operation identifier
!> \return rotation slot, or zero when no slot matches
! **************************************************************************************************
   INTEGER FUNCTION find_kpoint_rotation_slot(qs_kpoint, rotp) RESULT(rot_slot)
      TYPE(kpoint_type), POINTER                         :: qs_kpoint
      INTEGER, INTENT(IN)                                :: rotp

      INTEGER                                            :: irot, rot_abs

      rot_slot = 0
      rot_abs = ABS(rotp)
      IF (.NOT. ASSOCIATED(qs_kpoint%ibrot)) RETURN
      DO irot = 1, SIZE(qs_kpoint%ibrot)
         IF (rot_abs == qs_kpoint%ibrot(irot)) THEN
            rot_slot = irot
            RETURN
         END IF
      END DO

   END FUNCTION find_kpoint_rotation_slot

END MODULE kpoint_mo_symmetry_methods
