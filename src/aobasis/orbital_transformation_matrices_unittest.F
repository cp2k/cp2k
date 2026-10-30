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

   IMPLICIT NONE

   INTEGER, PARAMETER :: lmax = 5
   REAL(KIND=dp), PARAMETER :: tolerance = 1.0E-12_dp
   INTEGER :: i, l, n
   REAL(KIND=dp) :: rotmat(3, 3), residual(2*lmax + 1, 2*lmax + 1)
   TYPE(orbrotmat_type), DIMENSION(:), POINTER :: orbrot => NULL()

   CALL init_orbital_pointers(lmax)
   ! Rotation about (1,2,2), with cos(theta)=3/5 and sin(theta)=4/5.
   rotmat = RESHAPE([29.0_dp, 28.0_dp, -20.0_dp, -20.0_dp, 35.0_dp, 20.0_dp, &
                     28.0_dp, -4.0_dp, 35.0_dp], [3, 3])/45.0_dp
   CALL calculate_rotmat(orbrot, rotmat, lmax)
   IF (.NOT. (ABS(orbrot(0)%mat(1, 1) - 1.0_dp) <= tolerance)) ERROR STOP "Incorrect s rotation"
   ! The real p harmonics are ordered (y,z,x).
   rotmat = RESHAPE([35.0_dp, 20.0_dp, -20.0_dp, -4.0_dp, 35.0_dp, 28.0_dp, &
                     28.0_dp, -20.0_dp, 29.0_dp], [3, 3])/45.0_dp
   IF (.NOT. ALL(ABS(orbrot(1)%mat - rotmat) <= tolerance)) ERROR STOP "Incorrect p rotation"
   DO l = 0, lmax
      n = SIZE(orbrot(l)%mat, 1)
      residual(1:n, 1:n) = MATMUL(orbrot(l)%mat, TRANSPOSE(orbrot(l)%mat))
      DO i = 1, n
         residual(i, i) = residual(i, i) - 1.0_dp
      END DO
      IF (.NOT. ALL(ABS(residual(1:n, 1:n)) <= tolerance)) THEN
         PRINT *, "Nonorthogonal orbital rotation for l =", l
         ERROR STOP 1
      END IF
   END DO
   CALL release_rotmat(orbrot)
   CALL deallocate_orbital_pointers()

END PROGRAM orbital_transformation_matrices_unittest
