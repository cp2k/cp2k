!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_alignment
   USE trajana_kinds,                   ONLY: dp
   USE trajana_linear_algebra,          ONLY: diagonalize_symmetric

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: center_positions, fit_rotation, identity3

CONTAINS

   SUBROUTINE center_positions(positions, masses, center)
      REAL(dp), INTENT(INOUT)                            :: positions(:, :)
      REAL(dp), INTENT(IN)                               :: masses(:)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: center(3)

      INTEGER                                            :: atom
      REAL(dp)                                           :: current_center(3)

      current_center = 0.0_dp
      DO atom = 1, SIZE(masses)
         current_center = current_center + masses(atom)*positions(:, atom)
      END DO
      current_center = current_center/SUM(masses)
      DO atom = 1, SIZE(masses)
         positions(:, atom) = positions(:, atom) - current_center
      END DO
      IF (PRESENT(center)) center = current_center
   END SUBROUTINE center_positions

   SUBROUTINE fit_rotation(current, reference, masses, rotation, ierr)
      REAL(dp), INTENT(IN)                               :: current(:, :), reference(:, :), masses(:)
      REAL(dp), INTENT(OUT)                              :: rotation(3, 3)
      INTEGER, INTENT(OUT)                               :: ierr

      INTEGER                                            :: atom, column, row
      REAL(dp)                                           :: covariance(3, 3), quaternion(4), quaternion_matrix(4, 4)
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:), eigenvectors(:, :)

      rotation = identity3()
      ierr = 0
      IF (SIZE(masses) < 2) RETURN
      covariance = 0.0_dp
      DO atom = 1, SIZE(masses)
         DO row = 1, 3
            DO column = 1, 3
               covariance(row, column) = covariance(row, column) + &
                                         masses(atom)*current(row, atom)*reference(column, atom)
            END DO
         END DO
      END DO
      IF (MAXVAL(ABS(covariance)) <= TINY(1.0_dp)) RETURN
      CALL quaternion_eigenproblem(covariance, quaternion_matrix)
      CALL diagonalize_symmetric(quaternion_matrix, eigenvalues, eigenvectors, ierr)
      IF (ierr /= 0) RETURN
      quaternion = eigenvectors(:, 1)
      CALL quaternion_rotation(quaternion, rotation)
   END SUBROUTINE fit_rotation

   SUBROUTINE quaternion_eigenproblem(covariance, matrix)
      REAL(dp), INTENT(IN)                               :: covariance(3, 3)
      REAL(dp), INTENT(OUT)                              :: matrix(4, 4)

      REAL(dp)                                           :: trace, z(3)

      trace = covariance(1, 1) + covariance(2, 2) + covariance(3, 3)
      z = [covariance(2, 3) - covariance(3, 2), &
           covariance(3, 1) - covariance(1, 3), &
           covariance(1, 2) - covariance(2, 1)]
      matrix = 0.0_dp
      matrix(1, 1) = trace
      matrix(1, 2:4) = z
      matrix(2:4, 1) = z
      matrix(2:4, 2:4) = covariance + TRANSPOSE(covariance)
      matrix(2, 2) = matrix(2, 2) - trace
      matrix(3, 3) = matrix(3, 3) - trace
      matrix(4, 4) = matrix(4, 4) - trace
   END SUBROUTINE quaternion_eigenproblem

   SUBROUTINE quaternion_rotation(quaternion, rotation)
      REAL(dp), INTENT(IN)                               :: quaternion(4)
      REAL(dp), INTENT(OUT)                              :: rotation(3, 3)

      REAL(dp)                                           :: w, x, y, z

      w = quaternion(1)
      x = quaternion(2)
      y = quaternion(3)
      z = quaternion(4)
      rotation(1, :) = [w*w + x*x - y*y - z*z, 2.0_dp*(x*y - w*z), 2.0_dp*(x*z + w*y)]
      rotation(2, :) = [2.0_dp*(x*y + w*z), w*w - x*x + y*y - z*z, 2.0_dp*(y*z - w*x)]
      rotation(3, :) = [2.0_dp*(x*z - w*y), 2.0_dp*(y*z + w*x), w*w - x*x - y*y + z*z]
   END SUBROUTINE quaternion_rotation

   PURE FUNCTION identity3() RESULT(identity)
      REAL(dp)                                           :: identity(3, 3)

      identity = 0.0_dp
      identity(1, 1) = 1.0_dp
      identity(2, 2) = 1.0_dp
      identity(3, 3) = 1.0_dp
   END FUNCTION identity3

END MODULE trajana_alignment
