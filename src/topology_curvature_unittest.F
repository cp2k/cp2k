!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM topology_curvature_unittest
   USE kinds,                           ONLY: dp
   USE topology_curvature,              ONLY: curvature_density,&
                                              link_plaquette,&
                                              polar_link

   IMPLICIT NONE
   REAL(KIND=dp), PARAMETER :: pi = 3.1415926535897932384626433832795_dp, tolerance = 1.e-11_dp
   COMPLEX(KIND=dp) :: p(2, 2, 6), q(2, 2, 6), g(2, 2), link(2, 2)
   REAL(KIND=dp) :: value, maximum, minimum, coarse, fine
   INTEGER :: i, status
   p(:, :, :) = 0.0_dp
   DO i = 1, 6
      p(1, 1, i) = 1.0_dp
      p(2, 2, i) = 1.0_dp
   END DO
   p(1, 1, 1) = EXP(CMPLX(0.0_dp, 0.2_dp, dp))
   p(2, 2, 1) = EXP(CMPLX(0.0_dp, -0.1_dp, dp))
   p(1, 1, 6) = EXP(CMPLX(0.0_dp, 0.3_dp, dp))
   p(2, 2, 6) = EXP(CMPLX(0.0_dp, 0.4_dp, dp))
   CALL curvature_density(p, 4, pi/2.0_dp, value, maximum, status)
   IF (status /= 0 .OR. ABS(value - 0.02_dp/(4.0_dp*pi*pi)) > tolerance) ERROR STOP 'Incorrect C2 product trace'
   g(1, :) = [CMPLX(1.0_dp, 0.0_dp, dp), CMPLX(0.0_dp, 1.0_dp, dp)]/SQRT(2.0_dp)
   g(2, :) = [CMPLX(0.0_dp, 1.0_dp, dp), CMPLX(1.0_dp, 0.0_dp, dp)]/SQRT(2.0_dp)
   DO i = 1, 6
      q(:, :, i) = MATMUL(CONJG(TRANSPOSE(g)), MATMUL(p(:, :, i), g))
   END DO
   CALL curvature_density(q, 4, pi/2.0_dp, value, maximum, status)
   IF (status /= 0 .OR. ABS(value - 0.02_dp/(4.0_dp*pi*pi)) > tolerance) ERROR STOP 'C2 gauge covariance failed'
   CALL curvature_density(p(:, :, 1:1), 2, pi/2.0_dp, value, maximum, status)
   IF (status /= 0 .OR. ABS(value + 0.1_dp/(2.0_dp*pi)) > tolerance) ERROR STOP 'Incorrect determinant C1 phase'
   CALL curvature_density(p, 4, 0.1_dp, value, maximum, status)
   IF (status /= -3) ERROR STOP 'Unresolved C2 phase accepted'
   p(:, :, 1) = 0.0_dp
   CALL polar_link(p(:, :, 1), link, minimum, status, tolerance)
   IF (status == 0) ERROR STOP 'Singular link accepted'
   CALL curvature_density(p, 4, pi/2.0_dp, value, maximum, status)
   IF (status /= -2) ERROR STOP 'Nonunitary plaquette accepted'
   CALL dirac_pairing(6, coarse)
   CALL dirac_pairing(8, fine)
   IF (ABS(fine + 1.0_dp) >= ABS(coarse + 1.0_dp)) ERROR STOP 'C2 Dirac refinement failed'
   IF (ABS(fine + 0.742412846261164_dp) > 1.e-6_dp) ERROR STOP 'Incorrect non-Abelian Dirac C2'
CONTAINS

! **************************************************************************************************
!> \brief Non-Abelian occupied doublet of the 4D Dirac model, exact C2=-1 at mass -3.
!> \param n ...
!> \param RESULT ...
! **************************************************************************************************
   SUBROUTINE dirac_pairing(n, RESULT)
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), INTENT(OUT)                         :: result

      COMPLEX(KIND=dp)                                   :: gamma(4, 4, 5), h(4, 4), id(2, 2), &
                                                            overlap(2, 2), plaquettes(2, 2, 6), &
                                                            sx(2, 2), sy(2, 2), sz(2, 2), work(32)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: links(:, :, :, :), states(:, :, :)
      INTEGER                                            :: at_mu, at_nu, d, INDEX(4), mu, nu, &
                                                            other, p, rem, status, v
      REAL(KIND=dp)                                      :: evals(4), k(4), local_value, maximum, &
                                                            minimum, rwork(12)

      sx = 0.0_dp
      sx(1, 2) = 1.0_dp
      sx(2, 1) = 1.0_dp
      sy = CMPLX(0.0_dp, 1.0_dp, dp)*sx
      sy(1, 2) = -sy(1, 2)
      sz = 0.0_dp
      sz(1, 1) = 1.0_dp
      sz(2, 2) = -1.0_dp
      id = 0.0_dp
      id(1, 1) = 1.0_dp
      id(2, 2) = 1.0_dp
      CALL tensor(sx, sx, gamma(:, :, 1))
      CALL tensor(sx, sy, gamma(:, :, 2))
      CALL tensor(sx, sz, gamma(:, :, 3))
      CALL tensor(sy, id, gamma(:, :, 4))
      CALL tensor(sz, id, gamma(:, :, 5))
      ALLOCATE (states(4, 2, n**4), links(2, 2, 4, n**4))
      DO v = 1, n**4
         rem = v - 1
         DO d = 1, 4
            INDEX(d) = MOD(rem, n)
            rem = rem/n
         END DO
         k = 2.0_dp*pi*REAL(index, dp)/REAL(n, dp)
         h = CMPLX(0.0_dp, 0.0_dp, dp)
         DO d = 1, 4
            h = h + SIN(k(d))*gamma(:, :, d)
         END DO
         h = h + (-3.0_dp + SUM(COS(k)))*gamma(:, :, 5)
         CALL zheev('V', 'U', 4, h, 4, evals, work, SIZE(work), rwork, status)
         IF (status /= 0) ERROR STOP 'Dirac eigensolver failed'
         states(:, :, v) = h(:, 1:2)
      END DO
      DO v = 1, n**4
         DO d = 1, 4
            other = neighbour(v, d, n)
            overlap = MATMUL(CONJG(TRANSPOSE(states(:, :, v))), states(:, :, other))
            CALL polar_link(overlap, links(:, :, d, v), minimum, status, tolerance)
            IF (status /= 0) ERROR STOP 'Singular Dirac link'
         END DO
      END DO
      RESULT = 0.0_dp
      DO v = 1, n**4
         p = 0
         DO mu = 1, 4
            DO nu = mu + 1, 4
               p = p + 1
               at_mu = neighbour(v, mu, n)
               at_nu = neighbour(v, nu, n)
               CALL link_plaquette(links(:, :, mu, v), links(:, :, nu, at_mu), links(:, :, mu, at_nu), &
                                   links(:, :, nu, v), plaquettes(:, :, p))
            END DO
         END DO
         CALL curvature_density(plaquettes, 4, pi/2.0_dp, local_value, maximum, status)
         IF (status /= 0) ERROR STOP 'Invalid Dirac plaquette'
         RESULT = RESULT + local_value
      END DO
   END SUBROUTINE dirac_pairing

! **************************************************************************************************
!> \brief Periodic forward neighbour in a first-axis-fastest mesh.
!> \param v ...
!> \param d ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION neighbour(v, d, n) RESULT(other)
      INTEGER, INTENT(IN)                                :: v, d, n

      INTEGER                                            :: coordinate, stride

      stride = n**(d - 1)
      coordinate = MOD((v - 1)/stride, n)
      other = v + stride
      IF (coordinate == n - 1) other = other - n*stride
   END FUNCTION neighbour

! **************************************************************************************************
!> \brief Kronecker product of two 2x2 matrices.
!> \param a ...
!> \param b ...
!> \param c ...
! **************************************************************************************************
   SUBROUTINE tensor(a, b, c)
      COMPLEX(KIND=dp), INTENT(IN)                       :: a(2, 2), b(2, 2)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: c(4, 4)

      INTEGER                                            :: i, j

      DO j = 1, 2
         DO i = 1, 2
            c(2*i - 1:2*i, 2*j - 1:2*j) = a(i, j)*b
         END DO
      END DO
   END SUBROUTINE tensor
END PROGRAM topology_curvature_unittest
