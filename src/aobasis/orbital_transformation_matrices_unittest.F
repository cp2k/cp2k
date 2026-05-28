!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM orbital_transformation_matrices_unittest
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: deallocate_orbital_pointers,&
                                              init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: calculate_rotmat,&
                                              orbrotmat_type,&
                                              release_rotmat
#include "../base/base_uses.f90"

   IMPLICIT NONE

   INTEGER, PARAMETER                                 :: lmax = 5
   REAL(KIND=dp), PARAMETER                          :: eps = 1.0E-12_dp
   REAL(KIND=dp), DIMENSION(3, 3)                    :: rotmat
   TYPE(orbrotmat_type), DIMENSION(:), POINTER       :: orbrot => NULL()

   CALL init_orbital_pointers(lmax)

   rotmat = 0.0_dp
   rotmat(1, 1) = 1.0_dp
   rotmat(2, 2) = 1.0_dp
   rotmat(3, 3) = 1.0_dp
   CALL calculate_rotmat(orbrot, rotmat, lmax)
   CALL check_identity(orbrot)
   CALL check_orthogonal(orbrot)

   CALL make_z_rotation(0.7_dp, rotmat)
   CALL calculate_rotmat(orbrot, rotmat, lmax)
   CALL check_orthogonal(orbrot)

   CALL make_axis_rotation([1.0_dp, 2.0_dp, 3.0_dp], 0.37_dp, rotmat)
   CALL calculate_rotmat(orbrot, rotmat, lmax)
   CALL check_orthogonal(orbrot)

   CALL release_rotmat(orbrot)
   CALL deallocate_orbital_pointers()

CONTAINS

! **************************************************************************************************
!> \brief Check that the identity Cartesian rotation produces identity band matrices.
!> \param orbrot ...
! **************************************************************************************************
   SUBROUTINE check_identity(orbrot)
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot

      INTEGER                                            :: i, l, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ref

      DO l = LBOUND(orbrot, 1), UBOUND(orbrot, 1)
         n = SIZE(orbrot(l)%mat, 1)
         ALLOCATE (ref(n, n))
         ref = 0.0_dp
         DO i = 1, n
            ref(i, i) = 1.0_dp
         END DO
         IF (MAXVAL(ABS(orbrot(l)%mat - ref)) > eps) THEN
            CPABORT("Bad identity rotation")
         END IF
         DEALLOCATE (ref)
      END DO

   END SUBROUTINE check_identity

! **************************************************************************************************
!> \brief Check orthogonality of all band rotation matrices.
!> \param orbrot ...
! **************************************************************************************************
   SUBROUTINE check_orthogonal(orbrot)
      TYPE(orbrotmat_type), DIMENSION(:), POINTER        :: orbrot

      INTEGER                                            :: i, j, k, l, n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ident, prod

      DO l = LBOUND(orbrot, 1), UBOUND(orbrot, 1)
         n = SIZE(orbrot(l)%mat, 1)
         ALLOCATE (ident(n, n), prod(n, n))
         ident = 0.0_dp
         DO i = 1, n
            ident(i, i) = 1.0_dp
         END DO
         prod = 0.0_dp
         DO j = 1, n
            DO k = 1, n
               DO i = 1, n
                  prod(i, j) = prod(i, j) + orbrot(l)%mat(i, k)*orbrot(l)%mat(j, k)
               END DO
            END DO
         END DO
         IF (MAXVAL(ABS(prod - ident)) > eps) THEN
            CPABORT("Bad rotation orthogonality")
         END IF
         DEALLOCATE (ident, prod)
      END DO

   END SUBROUTINE check_orthogonal

! **************************************************************************************************
!> \brief Build a Cartesian rotation around the z axis.
!> \param angle ...
!> \param rotmat ...
! **************************************************************************************************
   SUBROUTINE make_z_rotation(angle, rotmat)
      REAL(KIND=dp), INTENT(IN)                          :: angle
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: rotmat

      REAL(KIND=dp)                                      :: c, s

      c = COS(angle)
      s = SIN(angle)
      rotmat = 0.0_dp
      rotmat(1, 1) = c
      rotmat(1, 2) = -s
      rotmat(2, 1) = s
      rotmat(2, 2) = c
      rotmat(3, 3) = 1.0_dp

   END SUBROUTINE make_z_rotation

! **************************************************************************************************
!> \brief Build a Cartesian rotation around an arbitrary axis.
!> \param axis ...
!> \param angle ...
!> \param rotmat ...
! **************************************************************************************************
   SUBROUTINE make_axis_rotation(axis, angle, rotmat)
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: axis
      REAL(KIND=dp), INTENT(IN)                          :: angle
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: rotmat

      REAL(KIND=dp)                                      :: c, norm, one_c, s, x, y, z

      norm = SQRT(SUM(axis**2))
      CPASSERT(norm > 0.0_dp)
      x = axis(1)/norm
      y = axis(2)/norm
      z = axis(3)/norm
      c = COS(angle)
      s = SIN(angle)
      one_c = 1.0_dp - c

      rotmat(1, 1) = c + x*x*one_c
      rotmat(1, 2) = x*y*one_c - z*s
      rotmat(1, 3) = x*z*one_c + y*s
      rotmat(2, 1) = y*x*one_c + z*s
      rotmat(2, 2) = c + y*y*one_c
      rotmat(2, 3) = y*z*one_c - x*s
      rotmat(3, 1) = z*x*one_c - y*s
      rotmat(3, 2) = z*y*one_c + x*s
      rotmat(3, 3) = c + z*z*one_c

   END SUBROUTINE make_axis_rotation

END PROGRAM orbital_transformation_matrices_unittest
