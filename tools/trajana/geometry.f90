!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_geometry
   USE trajana_kinds,                   ONLY: dp
   USE trajana_trajectory_types,        ONLY: cell_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cell_inverse, cell_volume, cart_to_fractional, fractional_to_cart
   PUBLIC :: minimum_image, wrap_cartesian, wrap_fractional, shortest_cell_height

CONTAINS

   PURE REAL(dp) FUNCTION determinant3(matrix)
      REAL(dp), INTENT(IN)                               :: matrix(3, 3)

      determinant3 = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
                     matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
                     matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
   END FUNCTION determinant3

   PURE REAL(dp) FUNCTION cell_volume(cell)
      TYPE(cell_type), INTENT(IN)                        :: cell

      cell_volume = ABS(determinant3(cell%h))
   END FUNCTION cell_volume

   SUBROUTINE cell_inverse(cell, inverse, ok)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(OUT)                              :: inverse(3, 3)
      LOGICAL, INTENT(OUT)                               :: ok

      REAL(dp)                                           :: determinant

      determinant = determinant3(cell%h)
      ok = ABS(determinant) > 100.0_dp*TINY(1.0_dp)
      inverse = 0.0_dp
      IF (.NOT. ok) RETURN

      inverse(1, 1) = cell%h(2, 2)*cell%h(3, 3) - cell%h(2, 3)*cell%h(3, 2)
      inverse(1, 2) = cell%h(1, 3)*cell%h(3, 2) - cell%h(1, 2)*cell%h(3, 3)
      inverse(1, 3) = cell%h(1, 2)*cell%h(2, 3) - cell%h(1, 3)*cell%h(2, 2)
      inverse(2, 1) = cell%h(2, 3)*cell%h(3, 1) - cell%h(2, 1)*cell%h(3, 3)
      inverse(2, 2) = cell%h(1, 1)*cell%h(3, 3) - cell%h(1, 3)*cell%h(3, 1)
      inverse(2, 3) = cell%h(1, 3)*cell%h(2, 1) - cell%h(1, 1)*cell%h(2, 3)
      inverse(3, 1) = cell%h(2, 1)*cell%h(3, 2) - cell%h(2, 2)*cell%h(3, 1)
      inverse(3, 2) = cell%h(1, 2)*cell%h(3, 1) - cell%h(1, 1)*cell%h(3, 2)
      inverse(3, 3) = cell%h(1, 1)*cell%h(2, 2) - cell%h(1, 2)*cell%h(2, 1)
      inverse = inverse/determinant
   END SUBROUTINE cell_inverse

   SUBROUTINE cart_to_fractional(cell, cartesian, fractional, ok)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: cartesian(3)
      REAL(dp), INTENT(OUT)                              :: fractional(3)
      LOGICAL, INTENT(OUT)                               :: ok

      REAL(dp)                                           :: inverse(3, 3)

      CALL cell_inverse(cell, inverse, ok)
      IF (ok) fractional = MATMUL(inverse, cartesian)
   END SUBROUTINE cart_to_fractional

   PURE SUBROUTINE fractional_to_cart(cell, fractional, cartesian)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: fractional(3)
      REAL(dp), INTENT(OUT)                              :: cartesian(3)

      cartesian = MATMUL(cell%h, fractional)
   END SUBROUTINE fractional_to_cart

   PURE SUBROUTINE wrap_fractional(fractional, periodic, centered)
      REAL(dp), INTENT(INOUT)                            :: fractional(3)
      LOGICAL, INTENT(IN)                                :: periodic(3)
      LOGICAL, INTENT(IN), OPTIONAL                      :: centered

      INTEGER                                            :: axis
      LOGICAL                                            :: use_centered

      use_centered = .FALSE.
      IF (PRESENT(centered)) use_centered = centered
      DO axis = 1, 3
         IF (.NOT. periodic(axis)) CYCLE
         IF (use_centered) THEN
            fractional(axis) = fractional(axis) - ANINT(fractional(axis))
         ELSE
            fractional(axis) = fractional(axis) - FLOOR(fractional(axis))
         END IF
      END DO
   END SUBROUTINE wrap_fractional

   SUBROUTINE wrap_cartesian(cell, cartesian, wrapped, ok, centered)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: cartesian(3)
      REAL(dp), INTENT(OUT)                              :: wrapped(3)
      LOGICAL, INTENT(OUT)                               :: ok
      LOGICAL, INTENT(IN), OPTIONAL                      :: centered

      REAL(dp)                                           :: fractional(3)

      CALL cart_to_fractional(cell, cartesian, fractional, ok)
      IF (.NOT. ok) RETURN
      CALL wrap_fractional(fractional, cell%periodic, centered)
      CALL fractional_to_cart(cell, fractional, wrapped)
   END SUBROUTINE wrap_cartesian

   SUBROUTINE minimum_image(cell, displacement, nearest, ok)
      TYPE(cell_type), INTENT(IN)                        :: cell
      REAL(dp), INTENT(IN)                               :: displacement(3)
      REAL(dp), INTENT(OUT)                              :: nearest(3)
      LOGICAL, INTENT(OUT)                               :: ok

      INTEGER                                            :: i, j, k
      REAL(dp)                                           :: base(3), best_squared, fractional(3), &
                                                            trial(3), trial_fractional(3), &
                                                            trial_squared

      CALL cart_to_fractional(cell, displacement, fractional, ok)
      IF (.NOT. ok) RETURN
      base = fractional
      WHERE (cell%periodic) base = base - ANINT(base)
      CALL fractional_to_cart(cell, base, nearest)
      best_squared = DOT_PRODUCT(nearest, nearest)

      ! Fractional rounding alone is not sufficient for a skewed triclinic cell.
      ! Search the neighboring images around that initial lattice point.
      DO i = -1, 1
         IF (.NOT. cell%periodic(1) .AND. i /= 0) CYCLE
         DO j = -1, 1
            IF (.NOT. cell%periodic(2) .AND. j /= 0) CYCLE
            DO k = -1, 1
               IF (.NOT. cell%periodic(3) .AND. k /= 0) CYCLE
               trial_fractional = base - REAL([i, j, k], dp)
               CALL fractional_to_cart(cell, trial_fractional, trial)
               trial_squared = DOT_PRODUCT(trial, trial)
               IF (trial_squared < best_squared) THEN
                  best_squared = trial_squared
                  nearest = trial
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE minimum_image

   REAL(dp) FUNCTION shortest_cell_height(cell)
      TYPE(cell_type), INTENT(IN)                        :: cell

      INTEGER                                            :: axis, other1, other2
      LOGICAL                                            :: found
      REAL(dp)                                           :: area, normal(3)

      shortest_cell_height = HUGE(1.0_dp)
      found = .FALSE.
      DO axis = 1, 3
         IF (.NOT. cell%periodic(axis)) CYCLE
         other1 = MOD(axis, 3) + 1
         other2 = MOD(axis + 1, 3) + 1
         normal = cross_product(cell%h(:, other1), cell%h(:, other2))
         area = NORM2(normal)
         IF (area > TINY(1.0_dp)) THEN
            shortest_cell_height = MIN(shortest_cell_height, cell_volume(cell)/area)
            found = .TRUE.
         END IF
      END DO
      IF (.NOT. found) shortest_cell_height = 0.0_dp
   END FUNCTION shortest_cell_height

   PURE FUNCTION cross_product(a, b) RESULT(c)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END FUNCTION cross_product

END MODULE trajana_geometry
