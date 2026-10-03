!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM uspp_fourier_unittest
   USE atom_upf,                        ONLY: atom_release_upf,&
                                              atom_upfpot_type
   USE kinds,                           ONLY: dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: indco,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_type
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate,&
                                              uspp_augmentation_fourier,&
                                              uspp_augmentation_init,&
                                              uspp_augmentation_integrate,&
                                              uspp_augmentation_type

   IMPLICIT NONE
   INTEGER, PARAMETER :: nb = 5, np = 19, nao = 7
   REAL(KIND=dp), PARAMETER :: alpha = 1.1_dp, volume = 137.0_dp, center(3) = [0.31_dp, -0.43_dp, 1.17_dp]
   TYPE(atom_upfpot_type) :: pot
   TYPE(upf_projector_type) :: projectors
   TYPE(uspp_augmentation_type) :: aug
   INTEGER :: grid, nr, i, j, l, pair, test
   REAL(KIND=dp) :: dx, g(3), error, max_error, vectors(3, 7), d(np, np), matrix(np, np), matrix_half(np, np)
   REAL(KIND=dp) :: weights(7), energy_grid, energy_matrix, b(nao, np), p(nao, nao), delta_s(nao, nao)
   COMPLEX(KIND=dp) :: qg(np, np), qminus(np, np), reference(np, np), rho(7), vg(7)
   COMPLEX(KIND=dp) :: strain(np, np, 3, 3), strain_reference(np, np, 3, 3)

   CALL cp2k_init()
   max_error = 0.0_dp
   DO grid = 1, 3
      pot%number_of_proj = nb
      pot%is_ultrasoft = .TRUE.
      ALLOCATE (pot%lbeta(nb), pot%dion(nb, nb), pot%qij(nb, nb))
      pot%lbeta(:) = [0, 1, 2, 3, 1]
      pot%dion = 0.0_dp
      ! Deliberately inconsistent PP_Q: the metric must use the implemented Q(G=0).
      pot%qij = 42.0_dp
      nr = 2501
      IF (grid == 2) nr = 2500
      ALLOCATE (pot%r(nr), pot%rab(nr), pot%beta(nr, nb), pot%qijl(nr, nb*(nb + 1)/2, 0:6))
      DO i = 1, nr
         IF (grid == 3) THEN
            dx = LOG(1.0_dp + 12.5_dp/0.001_dp)/(nr - 1)
            pot%r(i) = 0.001_dp*(EXP((i - 1)*dx) - 1.0_dp)
            pot%rab(i) = dx*(pot%r(i) + 0.001_dp)
         ELSE
            pot%r(i) = (i - 1)*0.005_dp
            pot%rab(i) = 0.005_dp
         END IF
      END DO
      pot%qijl = 0.0_dp
      DO j = 1, nb
         pot%beta(:, j) = pot%r**pot%lbeta(j)*EXP(-alpha*pot%r**2)
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               pot%qijl(:, pair, l) = (i + j)*0.1_dp*pot%r**(pot%lbeta(i) + pot%lbeta(j) + 2) &
                                      *EXP(-alpha*pot%r**2)
            END DO
         END DO
      END DO
      CALL upf_projector_init(pot, 3, projectors)
      CALL uspp_augmentation_init(pot, projectors, aug)
      IF (aug%nproj /= np) ERROR STOP 'Wrong magnetic-projector count'
      DO test = 1, 5
         SELECT CASE (test)
         CASE (1)
            g = 0.0_dp
         CASE (2)
            g = [1.e-8_dp, -2.e-8_dp, 3.e-8_dp]
         CASE (3)
            g = [0.2_dp, -0.7_dp, 0.9_dp]
         CASE (4)
            g = [3.7_dp, -2.3_dp, 1.5_dp]
         CASE (5)
            g = [0.0_dp, 0.0_dp, 2.8_dp]
         END SELECT
         CALL uspp_augmentation_fourier(aug, g, center, volume, qg, strain)
         CALL gaussian_reference(g, reference, strain_reference)
         error = volume*MAXVAL(ABS(qg - reference))
         max_error = MAX(max_error, error)
         IF (error > 2.e-10_dp) THEN
            WRITE (*, *) 'Gaussian Fourier error:', grid, test, error
            ERROR STOP 'USPP Fourier transform differs from Cartesian Gaussian integral'
         END IF
         error = volume*MAXVAL(ABS(strain - strain_reference))
         max_error = MAX(max_error, error)
         IF (error > 2.e-10_dp) THEN
            WRITE (*, *) 'Gaussian Fourier strain error:', grid, test, error
            ERROR STOP 'USPP strain derivative differs from Cartesian Gaussian integral'
         END IF
         CALL uspp_augmentation_fourier(aug, -g, center, volume, qminus)
         IF (MAXVAL(ABS(qminus - CONJG(qg))) > 2.e-13_dp) ERROR STOP 'Wrong reciprocal conjugation'
         IF (MAXVAL(ABS(qg - TRANSPOSE(qg))) > 2.e-13_dp) ERROR STOP 'Wrong projector symmetry'
      END DO

      vectors(:, 1) = 0.0_dp
      vectors(:, 2) = [0.3_dp, -0.7_dp, 1.1_dp]
      vectors(:, 3) = [2.1_dp, 0.2_dp, -0.4_dp]
      vectors(:, 4) = [0.0_dp, 0.0_dp, 2.8_dp]
      vectors(:, 5:7) = -vectors(:, 2:4)
      DO j = 1, np
         DO i = 1, nao
            b(i, j) = 0.03_dp*SIN(REAL(i*j, dp))
         END DO
      END DO
      p = 0.01_dp
      DO i = 1, nao
         p(i, i) = 0.3_dp
      END DO
      d = MATMUL(TRANSPOSE(b), MATMUL(p, b))
      delta_s = MATMUL(b, MATMUL(aug%q, TRANSPOSE(b)))
      rho = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      CALL uspp_augmentation_collocate(aug, vectors, center, volume, d, rho)
      IF (ABS(volume*REAL(rho(1), dp) - SUM(p*delta_s)) > 1.e-13_dp) ERROR STOP 'Charge/metric mismatch'
      vg(1) = CMPLX(0.37_dp, 0.0_dp, KIND=dp)
      vg(2:4) = [CMPLX(0.2_dp, -0.3_dp, dp), CMPLX(-0.5_dp, 0.8_dp, dp), CMPLX(1.1_dp, 0.4_dp, dp)]
      vg(5:7) = CONJG(vg(2:4))
      weights = 1.0_dp
      CALL uspp_augmentation_integrate(aug, vectors, center, volume, vg, weights, matrix)
      energy_grid = volume*REAL(SUM(CONJG(vg)*rho), dp)
      energy_matrix = SUM(d*matrix)
      IF (ABS(energy_grid - energy_matrix) > 1.e-13_dp) ERROR STOP 'Augmentation operations are not adjoint'
      weights(2:4) = 2.0_dp
      CALL uspp_augmentation_integrate(aug, vectors(:, :4), center, volume, vg(:4), weights(:4), matrix_half)
      IF (MAXVAL(ABS(matrix_half - matrix)) > 1.e-13_dp) ERROR STOP 'Wrong reciprocal multiplicities'
      vg = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      vg(1) = CMPLX(2.3_dp, 0.0_dp, KIND=dp)
      CALL uspp_augmentation_integrate(aug, vectors, center, volume, vg, weights, matrix)
      IF (MAXVAL(ABS(matrix - 2.3_dp*aug%q)) > 1.e-13_dp) ERROR STOP 'Wrong constant-potential shift'
      CALL atom_release_upf(pot)
   END DO
   WRITE (*, *) 'USPP Fourier: 15 analytic s/p/d/f comparisons, maximum error =', max_error
   ! Release allocatables before LeakSanitizer inspects the main program at exit.
   aug = uspp_augmentation_type()
   projectors = upf_projector_type()
   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Independent full Cartesian Gaussian Fourier integrals, without Gaunt or Bessel sums.
!> \param g Reciprocal vector.
!> \param fourier Analytic transform of solid-harmonic products times a Gaussian.
!> \param strain Independent strain derivative from one-higher Cartesian moments.
! **************************************************************************************************
   SUBROUTINE gaussian_reference(g, fourier, strain)
      REAL(KIND=dp), INTENT(IN)                          :: g(3)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: fourier(np, np), strain(np, np, 3, 3)

      COMPLEX(KIND=dp)                                   :: dterm(3), moments(0:7, 3), phase, term
      INTEGER                                            :: a, axis, c, i, j, k, li, lj, mi, mj, &
                                                            powers(3), raised(3), u, v

      DO axis = 1, 3
         moments(0, axis) = SQRT(pi/alpha)*EXP(-g(axis)**2/(4.0_dp*alpha))
         moments(1, axis) = CMPLX(0.0_dp, -g(axis)/(2.0_dp*alpha), dp)*moments(0, axis)
         DO k = 1, 6
            moments(k + 1, axis) = REAL(k, dp)/(2.0_dp*alpha)*moments(k - 1, axis) &
                                   + CMPLX(0.0_dp, -g(axis)/(2.0_dp*alpha), dp)*moments(k, axis)
         END DO
      END DO
      phase = EXP(CMPLX(0.0_dp, -DOT_PRODUCT(g, center), dp))/volume
      fourier = CMPLX(0.0_dp, 0.0_dp, dp)
      strain = CMPLX(0.0_dp, 0.0_dp, dp)
      DO v = 1, np
         j = projectors%radial_index(v)
         lj = pot%lbeta(j)
         mj = projectors%angular_index(v) - lj*lj
         DO u = 1, np
            i = projectors%radial_index(u)
            li = pot%lbeta(i)
            mi = projectors%angular_index(u) - li*li
            DO c = 1, nco(lj)
               DO a = 1, nco(li)
                  powers = indco(:, ncoset(li - 1) + a) + indco(:, ncoset(lj - 1) + c)
                  term = PRODUCT([moments(powers(1), 1), moments(powers(2), 2), moments(powers(3), 3)])
                  fourier(u, v) = fourier(u, v) + orbtramat(li)%slm(mi, a)*orbtramat(lj)%slm(mj, c)*term
                  DO axis = 1, 3
                     raised = powers
                     raised(axis) = raised(axis) + 1
                     dterm(axis) = CMPLX(0.0_dp, -1.0_dp, dp)* &
                                   PRODUCT([moments(raised(1), 1), moments(raised(2), 2), moments(raised(3), 3)])
                     strain(u, v, :, axis) = strain(u, v, :, axis) - &
                                             g*dterm(axis)*orbtramat(li)%slm(mi, a)*orbtramat(lj)%slm(mj, c)
                  END DO
               END DO
            END DO
            fourier(u, v) = fourier(u, v)*phase*(i + j)*0.1_dp
            strain(u, v, :, :) = strain(u, v, :, :)*phase*(i + j)*0.1_dp
            DO axis = 1, 3
               strain(u, v, axis, axis) = strain(u, v, axis, axis) - fourier(u, v)
            END DO
         END DO
      END DO
   END SUBROUTINE gaussian_reference
END PROGRAM uspp_fourier_unittest
