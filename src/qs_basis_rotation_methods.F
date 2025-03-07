!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_basis_rotation_methods
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: do_method_dftb
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_sym_type,&
                                              kpoint_type
   USE orbital_pointers,                ONLY: nso
   USE orbital_transformation_matrices, ONLY: calculate_rotmat,&
                                              orbrotmat_type,&
                                              release_rotmat
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (only in this module)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_basis_rotation_methods'

   ! Public subroutines

   PUBLIC :: qs_basis_rotation

CONTAINS

! **************************************************************************************************
!> \brief   Construct basis set rotation matrices
!> \param qs_env ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE qs_basis_rotation(qs_env, kpoints)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: ik, ikind, ir, ira, irot, jr, lval, &
                                                            nkind, nrot
      REAL(KIND=dp), DIMENSION(3, 3)                     :: rotmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
      TYPE(kpoint_sym_type), POINTER                     :: kpsym
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(kpoints))
      IF (ASSOCIATED(kpoints%kind_rotmat)) THEN
         CALL get_qs_env(qs_env, cell=cell)
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxlgto=lval)
         nrot = SIZE(kpoints%kind_rotmat, 1)
         nkind = SIZE(kpoints%kind_rotmat, 2)
         ! remove possible old rotation matrices
         DO irot = 1, nrot
            DO ikind = 1, nkind
               IF (ASSOCIATED(kpoints%kind_rotmat(irot, ikind)%rmat)) THEN
                  DEALLOCATE (kpoints%kind_rotmat(irot, ikind)%rmat)
               END IF
            END DO
         END DO
         ! check all rotations needed
         NULLIFY (orbrot)
         CALL get_qs_env(qs_env, dft_control=dft_control)
         DO ik = 1, kpoints%nkp
            kpsym => kpoints%kp_sym(ik)%kpoint_sym
            IF (kpsym%apply_symmetry) THEN
               DO irot = 1, SIZE(kpsym%rotp)
                  ir = kpsym%rotp(irot)
                  ira = 0
                  DO jr = 1, SIZE(kpoints%ibrot)
                     IF (ir == kpoints%ibrot(jr)) ira = jr
                  END DO
                  IF (ira > 0) THEN
                     IF (.NOT. ASSOCIATED(kpoints%kind_rotmat(ira, 1)%rmat)) THEN
                        rotmat(1:3, 1:3) = MATMUL(cell%h_inv, &
                                                  MATMUL(kpsym%rot(:, :, irot), cell%hmat))
                        CALL calculate_rotmat(orbrot, rotmat, lval)
                        IF (dft_control%qs_control%method_id == do_method_dftb) THEN
                           CPABORT("ROTMAT")
                        ELSE
                           DO ikind = 1, nkind
                              CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis)
                              NULLIFY (kpoints%kind_rotmat(ira, ikind)%rmat)
                              CALL set_rotmat_basis(kpoints%kind_rotmat(ira, ikind)%rmat, orbrot, orb_basis)
                           END DO
                        END IF
                     END IF
                  END IF
               END DO
            END IF
         END DO
         CALL release_rotmat(orbrot)
      END IF

   END SUBROUTINE qs_basis_rotation

! **************************************************************************************************
!> \brief ...
!> \param rmat ...
!> \param orbrot ...
!> \param basis ...
! **************************************************************************************************
   SUBROUTINE set_rotmat_basis(rmat, orbrot, basis)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rmat
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot
      TYPE(gto_basis_set_type), POINTER                  :: basis

      INTEGER                                            :: fs1, fs2, iset, ishell, l, nset, nsgf
      INTEGER, DIMENSION(:), POINTER                     :: nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf, lshell

      CALL get_gto_basis_set(gto_basis_set=basis, nsgf=nsgf)
      ALLOCATE (rmat(nsgf, nsgf))
      rmat = 0.0_dp

      CALL get_gto_basis_set(gto_basis_set=basis, nset=nset, nshell=nshell, l=lshell, &
                             first_sgf=first_sgf)
      DO iset = 1, nset
         DO ishell = 1, nshell(iset)
            l = lshell(ishell, iset)
            fs1 = first_sgf(ishell, iset)
            fs2 = fs1 + nso(l) - 1
            rmat(fs1:fs2, fs1:fs2) = orbrot(l)%mat(1:nso(l), 1:nso(l))
         END DO
      END DO

   END SUBROUTINE set_rotmat_basis

END MODULE qs_basis_rotation_methods
