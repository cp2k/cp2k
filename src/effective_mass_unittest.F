!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM effective_mass_unittest
   USE effective_mass_methods,          ONLY: &
        degenerate_range, displace_fractional, emass_carrier_electron, emass_carrier_flat, &
        emass_carrier_hole, emass_carrier_saddle, kato_rellich_reduce, principal_masses, &
        refine_bz_bounds, spectral_newton_step, sylvester_sqrt_tile, wrap_fractional
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi,&
                                              z_zero
   USE mathlib,                         ONLY: diag_complex,&
                                              geeig_right,&
                                              unit_matrix

   IMPLICIT NONE

   CALL check_kato_rellich()
   CALL check_principal_masses()
   CALL check_degenerate_range()
   CALL check_refinement_primitives()
   CALL check_sylvester_sqrt()

CONTAINS

! **************************************************************************************************
!> \brief Validate kato_rellich_reduce on an unsplit triplet, a complex 3-in-4 crossing against
!>        4th-order centered FD geeig_right, and dense U(3) gauge rotation.
! **************************************************************************************************
   SUBROUTINE check_kato_rellich()
      INTEGER, PARAMETER                                 :: ndeg = 3, norb = 4
      REAL(dp), PARAMETER :: e0 = 0.25_dp, e_spectator = -0.35_dp, eps_fd = 1.0E-8_dp, &
         hstep = 2.0E-3_dp, v_pair = 0.030_dp, v_third = -0.045_dp, v_unsplit = 0.12_dp

      COMPLEX(dp) :: aop(ndeg, ndeg), bop(ndeg, ndeg), c(norb, norb), c_tmp(norb, norb), &
         d2h(norb, norb), d2s(norb, norb), d_all(norb, ndeg), dh(norb, norb), ds(norb, norb), &
         h(norb, norb), s(norb, norb), sop(ndeg, ndeg), sos(ndeg, ndeg), umat(ndeg, ndeg)
      INTEGER                                            :: i, j, k
      LOGICAL                                            :: linear_split
      REAL(dp)                                           :: curv_fd(ndeg), evals(norb), &
                                                            expected_curv(ndeg), f0(ndeg), &
                                                            fm(2, ndeg), fp(2, ndeg), max_err, t, &
                                                            vel_fd(ndeg)
      REAL(dp), ALLOCATABLE                              :: inv_mass_eigs(:), inv_mass_rot(:), &
                                                            velocity_eigs(:), velocity_rot(:)

      ! 1. Unsplit triplet with equal velocity and Pulay sop coupling.
      CALL unit_matrix(aop)
      aop = v_unsplit*aop
      sop = RESHAPE([CMPLX(0.10_dp, 0.0_dp, dp), CMPLX(0.02_dp, -0.01_dp, dp), z_zero, &
                     CMPLX(0.02_dp, 0.01_dp, dp), CMPLX(-0.05_dp, 0.0_dp, dp), CMPLX(0.01_dp, 0.0_dp, dp), &
                     z_zero, CMPLX(0.01_dp, 0.0_dp, dp), CMPLX(0.08_dp, 0.0_dp, dp)], [ndeg, ndeg])
      bop = RESHAPE([CMPLX(-4.2_dp, 0.0_dp, dp), z_zero, z_zero, &
                     z_zero, CMPLX(-1.5_dp, 0.0_dp, dp), z_zero, &
                     z_zero, z_zero, CMPLX(2.3_dp, 0.0_dp, dp)], [ndeg, ndeg]) + 2.0_dp*v_unsplit*sop
      expected_curv = [-4.2_dp, -1.5_dp, 2.3_dp]
      CALL kato_rellich_reduce(aop, bop, sop, velocity_eigs, inv_mass_eigs, linear_split)
      IF (linear_split .OR. MAXVAL(ABS(velocity_eigs - v_unsplit)) > 1.0E-13_dp .OR. &
          MAXVAL(ABS(inv_mass_eigs - expected_curv)) > 1.0E-12_dp) THEN
         ERROR STOP 'Unsplit equal-velocity Pulay check failed'
      END IF

      ! 2. Complex 3-in-4 crossing (mirror pair v_pair, split branch v_third, spectator orbital).
      CALL unit_matrix(s)
      s(1, 2) = CMPLX(0.05_dp, 0.02_dp, dp)
      s(2, 1) = CONJG(s(1, 2))
      h = z_zero
      h(1:3, 1:3) = e0*s(1:3, 1:3)
      h(4, 4) = e_spectator

      DO j = 1, norb
         DO i = 1, j
            ds(i, j) = CMPLX(0.02_dp/REAL(ABS(i - j) + 1, dp), 0.005_dp*REAL(i - j, dp), dp)
            d2h(i, j) = CMPLX(0.30_dp/REAL(ABS(i - j) + 1, dp), 0.02_dp*REAL(i - j, dp), dp)
            d2s(i, j) = CMPLX(0.01_dp/REAL(ABS(i - j) + 1, dp), 0.002_dp*REAL(i - j, dp), dp)
            ds(j, i) = CONJG(ds(i, j))
            d2h(j, i) = CONJG(d2h(i, j))
            d2s(j, i) = CONJG(d2s(i, j))
         END DO
      END DO
      dh = e0*ds
      dh(1:2, 1:2) = dh(1:2, 1:2) + v_pair*s(1:2, 1:2)
      dh(3, 3) = dh(3, 3) + v_third
      dh(1:2, 4) = CMPLX(0.05_dp, 0.02_dp, dp)
      dh(4, 1:2) = CONJG(dh(1:2, 4))
      dh(3, 4) = CMPLX(-0.045_dp, 0.015_dp, dp)
      dh(4, 3) = CONJG(dh(3, 4))

      CALL geeig_right(h, s, evals, c)
      d_all(:, :) = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(dh - e0*ds, c(:, 2:norb)))
      aop(:, :) = d_all(2:norb, 1:ndeg)
      sop(:, :) = MATMUL(CONJG(TRANSPOSE(c(:, 2:norb))), MATMUL(ds, c(:, 2:norb)))
      sos(:, :) = 2.0_dp*MATMUL(CONJG(TRANSPOSE(d_all(1:1, 1:ndeg))), d_all(1:1, 1:ndeg))/ &
                  (e0 - evals(1))
      bop(:, :) = MATMUL(CONJG(TRANSPOSE(c(:, 2:norb))), MATMUL(d2h - e0*d2s, c(:, 2:norb))) + sos

      f0(:) = evals(2:norb)
      DO k = 1, 2
         t = REAL(k, dp)*hstep
         CALL geeig_right(h + t*dh + 0.5_dp*t**2*d2h, s + t*ds + 0.5_dp*t**2*d2s, evals, c_tmp)
         fp(k, :) = evals(2:norb)
         CALL geeig_right(h - t*dh + 0.5_dp*t**2*d2h, s - t*ds + 0.5_dp*t**2*d2s, evals, c_tmp)
         ! Negative t reverses v_third (-0.045) above the v_pair (+0.030) doublet while keeping doublet curvature order.
         fm(k, :) = evals([4, 2, 3])
      END DO
      vel_fd(:) = (fm(2, :) - 8.0_dp*fm(1, :) + 8.0_dp*fp(1, :) - fp(2, :))/(12.0_dp*hstep)
      curv_fd(:) = (-fm(2, :) + 16.0_dp*fm(1, :) - 30.0_dp*f0 + 16.0_dp*fp(1, :) - fp(2, :))/ &
                   (12.0_dp*hstep**2)

      CALL kato_rellich_reduce(aop, bop, sop, velocity_eigs, inv_mass_eigs, linear_split)
      IF (.NOT. linear_split) ERROR STOP 'Linear-split crossing failed'
      max_err = MAX(MAXVAL(ABS(vel_fd - velocity_eigs)), MAXVAL(ABS(curv_fd - inv_mass_eigs)))
      IF (max_err > eps_fd) ERROR STOP 'Linear-split FD mismatch'

      ! 3. Dense U(3) subspace gauge invariance.
      CALL diag_complex(bop, umat, vel_fd)
      CALL kato_rellich_reduce(MATMUL(CONJG(TRANSPOSE(umat)), MATMUL(aop, umat)), &
                               MATMUL(CONJG(TRANSPOSE(umat)), MATMUL(bop, umat)), &
                               MATMUL(CONJG(TRANSPOSE(umat)), MATMUL(sop, umat)), &
                               velocity_rot, inv_mass_rot, linear_split)
      IF (.NOT. linear_split .OR. MAXVAL(ABS(velocity_rot - velocity_eigs)) > 1.0E-13_dp .OR. &
          MAXVAL(ABS(inv_mass_rot - inv_mass_eigs)) > 1.0E-13_dp) THEN
         ERROR STOP 'Subspace unitary gauge invariance violated'
      END IF

      PRINT '(A,ES12.3)', 'Kato-Rellich FD max error: ', max_err
      PRINT *, 'Kato-Rellich reduction checks passed.'
   END SUBROUTINE check_kato_rellich

! **************************************************************************************************
!> \brief Validate principal_masses eigenvalues, eigenvectors, and carrier classification.
! **************************************************************************************************
   SUBROUTINE check_principal_masses()
      INTEGER                                            :: carrier
      REAL(dp)                                           :: dirs(3, 3), eigs(3), eye(3, 3), inv(3, 3)

      CALL unit_matrix(eye)
      inv = RESHAPE([1.2_dp, 0.3_dp, 0.0_dp, &
                     0.3_dp, 0.9_dp, 0.1_dp, &
                     0.0_dp, 0.1_dp, 1.5_dp], [3, 3])
      CALL principal_masses(inv, eigs, dirs, carrier)
      IF (carrier /= emass_carrier_electron) ERROR STOP '3D electron tensor misclassified'
      IF (MAXVAL(ABS(MATMUL(TRANSPOSE(dirs), dirs) - eye)) > 1.0E-13_dp) THEN
         ERROR STOP 'Principal directions not orthonormal'
      END IF

      CALL principal_masses(-inv, eigs, dirs, carrier)
      IF (carrier /= emass_carrier_hole) ERROR STOP '3D hole tensor misclassified'

      inv = 0.0_dp
      inv(1:2, 1:2) = RESHAPE([0.32_dp, 0.0_dp, 0.0_dp, 0.55_dp], [2, 2])
      inv(3, 3) = -1.0E-12_dp
      CALL principal_masses(inv, eigs, dirs, carrier)
      IF (carrier /= emass_carrier_electron) ERROR STOP 'Noisy 2D electron tensor misclassified'

      CALL principal_masses(0.0_dp*eye, eigs, dirs, carrier)
      IF (carrier /= emass_carrier_flat) ERROR STOP 'Flat tensor misclassified'

      inv(2, 2) = -0.55_dp
      inv(3, 3) = 0.48_dp
      CALL principal_masses(inv, eigs, dirs, carrier)
      IF (carrier /= emass_carrier_saddle .OR. eigs(1) >= eigs(2) .OR. eigs(2) >= eigs(3)) THEN
         ERROR STOP 'Saddle tensor classification/order failed'
      END IF

      PRINT *, 'Principal mass classification checks passed.'
   END SUBROUTINE check_principal_masses

! **************************************************************************************************
!> \brief Validate transitive adjacent-gap clustering and bounds in degenerate_range.
! **************************************************************************************************
   SUBROUTINE check_degenerate_range()
      REAL(dp), PARAMETER                                :: eps_deg = 1.0E-5_dp

      INTEGER                                            :: ihi, ilo
      REAL(dp)                                           :: energies(6)

      energies = [-0.50_dp, 0.00_dp, 0.8E-5_dp, 1.6E-5_dp, 0.10_dp, 0.20_dp]
      CALL degenerate_range(energies, 2, eps_deg, ilo, ihi)
      IF (ilo /= 2 .OR. ihi /= 4) ERROR STOP 'Chain anchor 2 failed'
      CALL degenerate_range(energies, 4, eps_deg, ilo, ihi)
      IF (ilo /= 2 .OR. ihi /= 4) ERROR STOP 'Chain anchor 4 failed'
      CALL degenerate_range(energies, 1, eps_deg, ilo, ihi)
      IF (ilo /= 1 .OR. ihi /= 1) ERROR STOP 'Lower boundary band failed'
      CALL degenerate_range(energies, 6, eps_deg, ilo, ihi)
      IF (ilo /= 6 .OR. ihi /= 6) ERROR STOP 'Upper boundary band failed'

      PRINT *, 'Degenerate range clustering checks passed.'
   END SUBROUTINE check_degenerate_range

! **************************************************************************************************
!> \brief Validate spectral_newton_step, refine_bz_bounds, displace_fractional, and wrap_fractional.
! **************************************************************************************************
   SUBROUTINE check_refinement_primitives()
      REAL(dp)                                           :: grad(3), hess(3, 3), hmat_2d(3, 3), &
                                                            kf_out(3), p(3), probe_step, step_cap

      grad = [0.02_dp, -0.04_dp, 0.01_dp]
      hess = RESHAPE([2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
      CALL spectral_newton_step(grad, hess, p)
      IF (MAXVAL(ABS(p - [-0.01_dp, 0.01_dp, -0.01_dp])) > 1.0E-14_dp) ERROR STOP 'Newton step wrong'
      CALL spectral_newton_step(10.0_dp*grad, hess, p, step_cap=0.05_dp)
      IF (ABS(NORM2(p) - 0.05_dp) > 1.0E-14_dp) ERROR STOP 'step_cap clamping failed'

      hess = RESHAPE([-2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [3, 3])
      CALL spectral_newton_step([0.02_dp, 0.05_dp, 1.0E-12_dp], hess, p)
      IF (p(1) >= 0.0_dp .OR. ABS(p(3)) > 0.0_dp) ERROR STOP 'Damped/vacuum step wrong'

      hmat_2d = RESHAPE([6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp], [3, 3])
      CALL refine_bz_bounds(hmat_2d, [.TRUE., .TRUE., .FALSE.], step_cap, probe_step)
      IF (ABS(step_cap - 0.20_dp) > 1.0E-14_dp) ERROR STOP '2D slab step_cap wrong'
      CALL refine_bz_bounds(hmat_2d, [.FALSE., .FALSE., .FALSE.], step_cap, probe_step)
      IF (ABS(step_cap - 0.20_dp) > 1.0E-14_dp .OR. ABS(probe_step - 0.05_dp) > 1.0E-14_dp) THEN
         ERROR STOP '0D molecule caps wrong'
      END IF

      kf_out = displace_fractional([0.45_dp, 0.0_dp, 0.30_dp], [twopi/6.0_dp, 0.0_dp, 1.0_dp], &
                                   hmat_2d, [.TRUE., .TRUE., .FALSE.])
      IF (MAXVAL(ABS(kf_out - [0.45_dp, 0.0_dp, 0.30_dp])) > 1.0E-14_dp) ERROR STOP 'displace_fractional wrong'
      IF (MAXVAL(ABS(wrap_fractional([0.5_dp, -1.25_dp, 1.25_dp]) - [-0.5_dp, -0.25_dp, 0.25_dp])) > 1.0E-14_dp) THEN
         ERROR STOP 'wrap_fractional boundary wrong'
      END IF

      PRINT *, 'Refinement primitive checks passed.'
   END SUBROUTINE check_refinement_primitives

! **************************************************************************************************
!> \brief Validate sylvester_sqrt_tile 1st/2nd Frechet derivatives against centered FD and
!>        U(d) eigenvector gauge invariance across a degenerate overlap spectrum.
! **************************************************************************************************
   SUBROUTINE check_sylvester_sqrt()
      INTEGER, PARAMETER                                 :: n = 5
      REAL(dp), PARAMETER                                :: fd_tol = 1.0E-8_dp, &
                                                            gauge_tol = 1.0E-14_dp, h = 6.0E-4_dp

      COMPLEX(dp) :: f0(n, n), f0_r(n, n), f1_a(n, n), f1_a_fd(n, n), f1_a_r(n, n), f1_b(n, n), &
         f1_b_fd(n, n), f1_b_r(n, n), f2_aa(n, n), f2_aa_fd(n, n), f2_ab(n, n), f2_ab_fd(n, n), &
         f2_ab_r(n, n), fm(n, n), fp(n, n), gen(n, n), q0(n, n), s0(n, n), s1_a(n, n), s1_b(n, n), &
         s2_aa(n, n), s2_ab(n, n), s2_bb(n, n), u(n, n), u_rot(n, n), w_gauge(n, n)
      INTEGER                                            :: i, j
      REAL(dp)                                           :: evals(n), max_fd_err, max_gauge_err, &
                                                            sigma(n)

      DO j = 1, n
         DO i = 1, n
            gen(i, j) = CMPLX(0.3_dp*REAL(i + j, dp), 0.09_dp*REAL(j - i, dp), dp)
            s1_a(i, j) = CMPLX(0.07_dp*COS(REAL(i + 2*j, dp)), 0.05_dp*SIN(REAL(j - i, dp)), dp)
            s1_b(i, j) = CMPLX(-0.06_dp*SIN(REAL(2*i + j, dp)), 0.04_dp*SIN(REAL(i - j, dp))*REAL(i + j, dp)/6.0_dp, dp)
            s2_aa(i, j) = CMPLX(0.09_dp*COS(REAL(i - j, dp)), 0.03_dp*SIN(REAL(j - i, dp)), dp)
            s2_ab(i, j) = CMPLX(0.05_dp*COS(REAL(2*i - j, dp)), -0.04_dp*SIN(REAL(j - i, dp)), dp)
            s2_bb(i, j) = CMPLX(0.08_dp*COS(REAL(i + j, dp)), 0.02_dp*SIN(REAL(j - i, dp)), dp)
         END DO
      END DO
      s1_a = 0.5_dp*(s1_a + CONJG(TRANSPOSE(s1_a)))
      s1_b = 0.5_dp*(s1_b + CONJG(TRANSPOSE(s1_b)))
      s2_aa = 0.5_dp*(s2_aa + CONJG(TRANSPOSE(s2_aa)))
      s2_ab = 0.5_dp*(s2_ab + CONJG(TRANSPOSE(s2_ab)))
      s2_bb = 0.5_dp*(s2_bb + CONJG(TRANSPOSE(s2_bb)))

      ! Build S0 with exact 2-fold (0.64) and 3-fold (1.44) degenerate overlap eigenvalues.
      CALL diag_complex(gen + CONJG(TRANSPOSE(gen)), q0, evals)
      sigma = [0.80_dp, 0.80_dp, 1.20_dp, 1.20_dp, 1.20_dp]
      s0 = MATMUL(q0*SPREAD(sigma**2, 1, n), CONJG(TRANSPOSE(q0)))
      CALL diag_complex(s0, u, evals)
      sigma = SQRT(evals)

      CALL sylvester_derivs(u, sigma, s1_a, s1_b, s2_ab, f0, f1_a, f1_b, f2_ab)
      CALL sylvester_derivs(u, sigma, s1_a, s1_a, s2_aa, f0_r, f1_a_r, f1_b_r, f2_aa)

      fp = mat_sqrt(s0 + h*s1_a + 0.5_dp*h*h*s2_aa)
      fm = mat_sqrt(s0 - h*s1_a + 0.5_dp*h*h*s2_aa)
      f1_a_fd = (fp - fm)/(2.0_dp*h)
      f2_aa_fd = (fp - 2.0_dp*mat_sqrt(s0) + fm)/(h*h)
      f1_b_fd = (mat_sqrt(s0 + h*s1_b + 0.5_dp*h*h*s2_bb) - &
                 mat_sqrt(s0 - h*s1_b + 0.5_dp*h*h*s2_bb))/(2.0_dp*h)
      f2_ab_fd = (mat_sqrt(s0 + h*(s1_a + s1_b) + 0.5_dp*h*h*(s2_aa + 2.0_dp*s2_ab + s2_bb)) - &
                  mat_sqrt(s0 + h*(s1_a - s1_b) + 0.5_dp*h*h*(s2_aa - 2.0_dp*s2_ab + s2_bb)) - &
                  mat_sqrt(s0 + h*(-s1_a + s1_b) + 0.5_dp*h*h*(s2_aa - 2.0_dp*s2_ab + s2_bb)) + &
                  mat_sqrt(s0 - h*(s1_a + s1_b) + 0.5_dp*h*h*(s2_aa + 2.0_dp*s2_ab + s2_bb)))/(4.0_dp*h*h)

      max_fd_err = MAX(MAXVAL(ABS(f0 - mat_sqrt(s0))), MAXVAL(ABS(f1_a - f1_a_fd)), &
                       MAXVAL(ABS(f1_b - f1_b_fd)), MAXVAL(ABS(f2_aa - f2_aa_fd)), &
                       MAXVAL(ABS(f2_ab - f2_ab_fd)))
      WRITE (*, '(A,ES12.3)') 'Sylvester FD max error:     ', max_fd_err
      IF (max_fd_err > fd_tol) ERROR STOP 'sylvester_sqrt_tile FD mismatch'

      ! Dense U(2) x U(3) gauge rotation inside the degenerate overlap eigenspaces.
      w_gauge = z_zero
      CALL diag_complex(s1_a(1:2, 1:2), w_gauge(1:2, 1:2), evals(1:2))
      CALL diag_complex(s1_b(3:5, 3:5), w_gauge(3:5, 3:5), evals(3:5))
      u_rot(:, :) = MATMUL(u, w_gauge)
      CALL sylvester_derivs(u_rot, sigma, s1_a, s1_b, s2_ab, f0_r, f1_a_r, f1_b_r, f2_ab_r)
      max_gauge_err = MAX(MAXVAL(ABS(f0_r - f0)), MAXVAL(ABS(f1_a_r - f1_a)), &
                          MAXVAL(ABS(f1_b_r - f1_b)), MAXVAL(ABS(f2_ab_r - f2_ab)))
      IF (max_gauge_err > gauge_tol) ERROR STOP 'sylvester_sqrt_tile U(d) gauge invariance violated'

      PRINT *, 'Sylvester S^(1/2) Frechet and U(d) gauge checks passed.'
   END SUBROUTINE check_sylvester_sqrt

! **************************************************************************************************
!> \brief Principal Hermitian matrix square root S^(1/2) via diag_complex.
!> \param spt Positive-definite Hermitian matrix (n x n)
!> \return Matrix square root S^(1/2)
! **************************************************************************************************
   FUNCTION mat_sqrt(spt) RESULT(fmat)
      COMPLEX(dp), INTENT(IN)                            :: spt(:, :)
      COMPLEX(dp) :: fmat(SIZE(spt, 1), SIZE(spt, 1))

      COMPLEX(dp)                                        :: v(SIZE(spt, 1), SIZE(spt, 1))
      REAL(dp)                                           :: ev(SIZE(spt, 1))

      CALL diag_complex(spt, v, ev)
      IF (ANY(ev <= 0.0_dp)) ERROR STOP 'mat_sqrt non-positive eigenvalue'
      fmat = MATMUL(v*SPREAD(SQRT(ev), 1, SIZE(spt, 1)), CONJG(TRANSPOSE(v)))
   END FUNCTION mat_sqrt

! **************************************************************************************************
!> \brief Evaluate 0th, 1st, and 2nd Sylvester Frechet derivatives of S^(1/2) from (u_in, sig_in).
!> \param u_in Unitary overlap eigenvectors
!> \param sig_in Positive overlap square-root eigenvalues
!> \param sa First derivative along a
!> \param sb First derivative along b
!> \param sab Mixed second derivative along (a, b)
!> \param f0_out Overlap square root S^(1/2)
!> \param f1a_out First Frechet derivative along a
!> \param f1b_out First Frechet derivative along b
!> \param f2ab_out Mixed second Frechet derivative along (a, b)
! **************************************************************************************************
   SUBROUTINE sylvester_derivs(u_in, sig_in, sa, sb, sab, f0_out, f1a_out, f1b_out, f2ab_out)
      COMPLEX(dp), INTENT(IN)                            :: u_in(:, :)
      REAL(dp), INTENT(IN)                               :: sig_in(:)
      COMPLEX(dp), INTENT(IN)                            :: sa(:, :), sb(:, :), sab(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: f0_out(:, :), f1a_out(:, :), &
                                                            f1b_out(:, :), f2ab_out(:, :)

      COMPLEX(dp) :: uh(SIZE(sig_in), SIZE(sig_in)), fb(SIZE(sig_in), SIZE(sig_in)), &
         fab(SIZE(sig_in), SIZE(sig_in)), fa(SIZE(sig_in), SIZE(sig_in))

      uh = CONJG(TRANSPOSE(u_in))
      f0_out = MATMUL(u_in*SPREAD(sig_in, 1, SIZE(sig_in)), uh)
      fa = MATMUL(uh, MATMUL(sa, u_in))
      CALL sylvester_sqrt_tile(fa, sig_in, sig_in)
      f1a_out = MATMUL(u_in, MATMUL(fa, uh))
      fb = MATMUL(uh, MATMUL(sb, u_in))
      CALL sylvester_sqrt_tile(fb, sig_in, sig_in)
      f1b_out = MATMUL(u_in, MATMUL(fb, uh))
      fab = MATMUL(uh, MATMUL(sab, u_in)) - MATMUL(fa, fb) - MATMUL(fb, fa)
      CALL sylvester_sqrt_tile(fab, sig_in, sig_in)
      f2ab_out = MATMUL(u_in, MATMUL(fab, uh))
   END SUBROUTINE sylvester_derivs

END PROGRAM effective_mass_unittest
