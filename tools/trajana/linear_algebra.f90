!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_linear_algebra
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: diagonalize_symmetric

CONTAINS

   SUBROUTINE diagonalize_symmetric(matrix, eigenvalues, eigenvectors, ierr)
      REAL(dp), INTENT(IN)                               :: matrix(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: eigenvalues(:), eigenvectors(:, :)
      INTEGER, INTENT(OUT)                               :: ierr

      INTEGER                                            :: i, p, q, sweep, n
      REAL(dp)                                           :: c, h, s, scale, sm, t, tau, theta, threshold
      REAL(dp), ALLOCATABLE                              :: b(:), work(:, :), z(:)

      n = SIZE(matrix, 1)
      ierr = 0
      IF (n < 1 .OR. SIZE(matrix, 2) /= n) THEN
         ierr = 1
         ALLOCATE (eigenvalues(0), eigenvectors(0, 0))
         RETURN
      END IF
      ALLOCATE (eigenvalues(n), eigenvectors(n, n), b(n), work(n, n), z(n))
      work = 0.5_dp*(matrix + TRANSPOSE(matrix))
      eigenvectors = 0.0_dp
      DO i = 1, n
         eigenvectors(i, i) = 1.0_dp
      END DO
      eigenvalues = [(work(i, i), i=1, n)]
      b = eigenvalues
      z = 0.0_dp
      DO sweep = 1, 100
         sm = 0.0_dp
         DO p = 1, n - 1
            DO q = p + 1, n
               sm = sm + ABS(work(p, q))
            END DO
         END DO
         scale = MAX(MAXVAL(ABS(eigenvalues)), 1.0_dp)
         IF (sm <= 1.0E-13_dp*scale) EXIT
         IF (sweep < 4) THEN
            threshold = 0.2_dp*sm/REAL(n*n, dp)
         ELSE
            threshold = 0.0_dp
         END IF

         DO p = 1, n - 1
            DO q = p + 1, n
               IF (sweep > 4 .AND. ABS(work(p, q)) <= EPSILON(1.0_dp)* &
                   MAX(ABS(eigenvalues(p)), ABS(eigenvalues(q)), 1.0_dp)) THEN
                  work(p, q) = 0.0_dp
               ELSE IF (ABS(work(p, q)) > threshold) THEN
                  h = eigenvalues(q) - eigenvalues(p)
                  IF (ABS(h) > TINY(1.0_dp) .AND. &
                      ABS(work(p, q)) <= EPSILON(1.0_dp)*ABS(h)) THEN
                     t = work(p, q)/h
                  ELSE
                     theta = 0.5_dp*h/work(p, q)
                     t = SIGN(1.0_dp, theta)/(ABS(theta) + SQRT(1.0_dp + theta*theta))
                  END IF
                  c = 1.0_dp/SQRT(1.0_dp + t*t)
                  s = t*c
                  tau = s/(1.0_dp + c)
                  h = t*work(p, q)
                  z(p) = z(p) - h
                  z(q) = z(q) + h
                  eigenvalues(p) = eigenvalues(p) - h
                  eigenvalues(q) = eigenvalues(q) + h
                  work(p, q) = 0.0_dp
                  DO i = 1, p - 1
                     CALL rotate_pair(work(i, p), work(i, q), s, tau)
                  END DO
                  DO i = p + 1, q - 1
                     CALL rotate_pair(work(p, i), work(i, q), s, tau)
                  END DO
                  DO i = q + 1, n
                     CALL rotate_pair(work(p, i), work(q, i), s, tau)
                  END DO
                  DO i = 1, n
                     CALL rotate_pair(eigenvectors(i, p), eigenvectors(i, q), s, tau)
                  END DO
               END IF
            END DO
         END DO
         b = b + z
         eigenvalues = b
         z = 0.0_dp
      END DO
      IF (sweep > 100) ierr = 2
      CALL sort_descending(eigenvalues, eigenvectors)
      CALL set_vector_signs(eigenvectors)
   END SUBROUTINE diagonalize_symmetric

   SUBROUTINE rotate_pair(first, second, sine, tau)
      REAL(dp), INTENT(INOUT)                            :: first, second
      REAL(dp), INTENT(IN)                               :: sine, tau

      REAL(dp)                                           :: old_first, old_second

      old_first = first
      old_second = second
      first = old_first - sine*(old_second + old_first*tau)
      second = old_second + sine*(old_first - old_second*tau)
   END SUBROUTINE rotate_pair

   SUBROUTINE sort_descending(values, vectors)
      REAL(dp), INTENT(INOUT)                            :: values(:), vectors(:, :)

      INTEGER                                            :: i, j, largest
      REAL(dp)                                           :: temporary
      REAL(dp), ALLOCATABLE                              :: vector(:)

      ALLOCATE (vector(SIZE(values)))
      DO i = 1, SIZE(values) - 1
         largest = i
         DO j = i + 1, SIZE(values)
            IF (values(j) > values(largest)) largest = j
         END DO
         IF (largest == i) CYCLE
         temporary = values(i)
         values(i) = values(largest)
         values(largest) = temporary
         vector = vectors(:, i)
         vectors(:, i) = vectors(:, largest)
         vectors(:, largest) = vector
      END DO
   END SUBROUTINE sort_descending

   SUBROUTINE set_vector_signs(vectors)
      REAL(dp), INTENT(INOUT)                            :: vectors(:, :)

      INTEGER                                            :: column, location(1)

      DO column = 1, SIZE(vectors, 2)
         location = MAXLOC(ABS(vectors(:, column)))
         IF (vectors(location(1), column) < 0.0_dp) vectors(:, column) = -vectors(:, column)
      END DO
   END SUBROUTINE set_vector_signs

END MODULE trajana_linear_algebra
