!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM effective_mass_unittest
   USE effective_mass_methods,          ONLY: degenerate_range,&
                                              emass_carrier_electron,&
                                              emass_carrier_flat,&
                                              emass_carrier_hole,&
                                              emass_carrier_saddle,&
                                              kato_rellich_reduce,&
                                              principal_masses,&
                                              spectral_newton_step
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_zero
   USE mathlib,                         ONLY: diag_complex,&
                                              geeig_right,&
                                              unit_matrix

   IMPLICIT NONE

   INTEGER, PARAMETER :: ndeg = 3, norb = 4
   REAL(dp), PARAMETER :: e0 = 0.25_dp, e_spectator = -0.35_dp, eps_fd = 1.0E-8_dp, &
                          hstep = 2.0E-3_dp, v_pair = 0.030_dp, v_third = -0.045_dp, v_unsplit = 0.12_dp

   COMPLEX(dp) :: aop(ndeg, ndeg), bop(ndeg, ndeg), c(norb, norb), c_tmp(norb, norb), &
                  d2h(norb, norb), d2s(norb, norb), d_all(norb, ndeg), dh(norb, norb), ds(norb, norb), &
                  h(norb, norb), s(norb, norb), sop(ndeg, ndeg), sos(ndeg, ndeg), umat(ndeg, ndeg)
   INTEGER :: carrier, i, ihi, ilo, j, k
   LOGICAL :: linear_split
   REAL(dp) :: curv_fd(ndeg), dirs(3, 3), eigs(3), evals(norb), f0(ndeg), fm(2, ndeg), &
               fp(2, ndeg), grad(3), hess(3, 3), inv(3, 3), p(3), t, vel_fd(ndeg)
   REAL(dp), ALLOCATABLE :: inv_mass_eigs(:), inv_mass_rot(:), velocity_eigs(:), velocity_rot(:)

   ! 1. Unsplit triplet with equal velocity and Pulay sop coupling (linear_split = .FALSE.).
   CALL unit_matrix(aop)
   aop = v_unsplit*aop
   sop = RESHAPE([CMPLX(0.10_dp, 0.0_dp, dp), CMPLX(0.02_dp, -0.01_dp, dp), z_zero, &
                  CMPLX(0.02_dp, 0.01_dp, dp), CMPLX(-0.05_dp, 0.0_dp, dp), CMPLX(0.01_dp, 0.0_dp, dp), &
                  z_zero, CMPLX(0.01_dp, 0.0_dp, dp), CMPLX(0.08_dp, 0.0_dp, dp)], [ndeg, ndeg])
   bop = RESHAPE([CMPLX(-4.2_dp, 0.0_dp, dp), z_zero, z_zero, z_zero, CMPLX(-1.5_dp, 0.0_dp, dp), &
                  z_zero, z_zero, z_zero, CMPLX(2.3_dp, 0.0_dp, dp)], [ndeg, ndeg]) + 2.0_dp*v_unsplit*sop
   CALL kato_rellich_reduce(aop, bop, sop, velocity_eigs, inv_mass_eigs, linear_split)
   IF (linear_split .OR. MAXVAL(ABS(velocity_eigs - v_unsplit)) > 1.0E-13_dp .OR. &
       MAXVAL(ABS(inv_mass_eigs - [-4.2_dp, -1.5_dp, 2.3_dp])) > 1.0E-12_dp) THEN
      ERROR STOP "Unsplit equal-velocity Pulay check failed"
   END IF

   ! 2. Complex 3-in-4 crossing against 4th-order centered FD geeig_right and U(3) gauge invariance.
   CALL unit_matrix(s)
   s(1, 2) = CMPLX(0.05_dp, 0.02_dp, dp); s(2, 1) = CONJG(s(1, 2))
   h = z_zero; h(1:3, 1:3) = e0*s(1:3, 1:3); h(4, 4) = e_spectator
   DO j = 1, norb
      DO i = 1, j
         ds(i, j) = CMPLX(0.02_dp/REAL(ABS(i - j) + 1, dp), 0.005_dp*REAL(i - j, dp), dp)
         d2h(i, j) = CMPLX(0.30_dp/REAL(ABS(i - j) + 1, dp), 0.02_dp*REAL(i - j, dp), dp)
         d2s(i, j) = CMPLX(0.01_dp/REAL(ABS(i - j) + 1, dp), 0.002_dp*REAL(i - j, dp), dp)
         ds(j, i) = CONJG(ds(i, j)); d2h(j, i) = CONJG(d2h(i, j)); d2s(j, i) = CONJG(d2s(i, j))
      END DO
   END DO
   dh = e0*ds
   dh(1:2, 1:2) = dh(1:2, 1:2) + v_pair*s(1:2, 1:2); dh(3, 3) = dh(3, 3) + v_third
   dh(1:2, 4) = CMPLX(0.05_dp, 0.02_dp, dp); dh(4, 1:2) = CONJG(dh(1:2, 4))
   dh(3, 4) = CMPLX(-0.045_dp, 0.015_dp, dp); dh(4, 3) = CONJG(dh(3, 4))

   CALL geeig_right(h, s, evals, c)
   d_all = MATMUL(CONJG(TRANSPOSE(c)), MATMUL(dh - e0*ds, c(:, 2:norb)))
   aop = d_all(2:norb, 1:ndeg)
   sop = MATMUL(CONJG(TRANSPOSE(c(:, 2:norb))), MATMUL(ds, c(:, 2:norb)))
   sos = 2.0_dp*MATMUL(CONJG(TRANSPOSE(d_all(1:1, 1:ndeg))), d_all(1:1, 1:ndeg))/(e0 - evals(1))
   bop = MATMUL(CONJG(TRANSPOSE(c(:, 2:norb))), MATMUL(d2h - e0*d2s, c(:, 2:norb))) + sos
   f0 = evals(2:norb)
   DO k = 1, 2
      t = REAL(k, dp)*hstep
      CALL geeig_right(h + t*dh + 0.5_dp*t**2*d2h, s + t*ds + 0.5_dp*t**2*d2s, evals, c_tmp)
      fp(k, :) = evals(2:norb)
      CALL geeig_right(h - t*dh + 0.5_dp*t**2*d2h, s - t*ds + 0.5_dp*t**2*d2s, evals, c_tmp)
      ! Negative t reverses v_third (-0.045) above the v_pair (+0.030) doublet while keeping doublet order.
      fm(k, :) = evals([4, 2, 3])
   END DO
   vel_fd = (fm(2, :) - 8.0_dp*fm(1, :) + 8.0_dp*fp(1, :) - fp(2, :))/(12.0_dp*hstep)
   curv_fd = (-fm(2, :) + 16.0_dp*fm(1, :) - 30.0_dp*f0 + 16.0_dp*fp(1, :) - fp(2, :))/(12.0_dp*hstep**2)
   CALL kato_rellich_reduce(aop, bop, sop, velocity_eigs, inv_mass_eigs, linear_split)
   IF (.NOT. linear_split .OR. MAXVAL(ABS(vel_fd - velocity_eigs)) > eps_fd .OR. &
       MAXVAL(ABS(curv_fd - inv_mass_eigs)) > eps_fd) ERROR STOP "Linear-split FD mismatch"

   CALL diag_complex(bop, umat, vel_fd)
   CALL kato_rellich_reduce(MATMUL(CONJG(TRANSPOSE(umat)), MATMUL(aop, umat)), &
                            MATMUL(CONJG(TRANSPOSE(umat)), MATMUL(bop, umat)), &
                            MATMUL(CONJG(TRANSPOSE(umat)), MATMUL(sop, umat)), &
                            velocity_rot, inv_mass_rot, linear_split)
   IF (.NOT. linear_split .OR. MAXVAL(ABS(velocity_rot - velocity_eigs)) > 1.0E-13_dp .OR. &
       MAXVAL(ABS(inv_mass_rot - inv_mass_eigs)) > 1.0E-13_dp) ERROR STOP "Subspace U(3) invariance failed"
   DEALLOCATE (velocity_eigs, inv_mass_eigs, velocity_rot, inv_mass_rot)

   ! 3. Transitive degenerate_range, principal_masses classification, and damped/vacuum Newton step.
   CALL degenerate_range([-0.50_dp, 0.00_dp, 0.8E-5_dp, 1.6E-5_dp, 0.10_dp, 0.20_dp], 4, 1.0E-5_dp, ilo, ihi)
   IF (ilo /= 2 .OR. ihi /= 4) ERROR STOP "Transitive degenerate_range failed"

   inv = RESHAPE([1.2_dp, 0.3_dp, 0.0_dp, 0.3_dp, 0.9_dp, 0.1_dp, 0.0_dp, 0.1_dp, 1.5_dp], [3, 3])
   CALL principal_masses(inv, eigs, dirs, carrier)
   IF (carrier /= emass_carrier_electron) ERROR STOP "3D electron tensor misclassified"
   CALL principal_masses(-inv, eigs, dirs, carrier)
   IF (carrier /= emass_carrier_hole) ERROR STOP "3D hole tensor misclassified"
   CALL principal_masses(RESHAPE([0.32_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.55_dp, 0.0_dp, &
                                  0.0_dp, 0.0_dp, -1.0E-12_dp], [3, 3]), eigs, dirs, carrier)
   IF (carrier /= emass_carrier_electron) ERROR STOP "Noisy 2D electron tensor misclassified"
   CALL principal_masses(RESHAPE([0.32_dp, 0.0_dp, 0.0_dp, 0.0_dp, -0.55_dp, 0.0_dp, &
                                  0.0_dp, 0.0_dp, 0.48_dp], [3, 3]), eigs, dirs, carrier)
   IF (carrier /= emass_carrier_saddle .OR. eigs(1) >= eigs(2) .OR. eigs(2) >= eigs(3)) ERROR STOP "Saddle order failed"
   CALL principal_masses(0.0_dp*inv, eigs, dirs, carrier)
   IF (carrier /= emass_carrier_flat) ERROR STOP "Flat tensor misclassified"

   grad = [0.02_dp, -0.04_dp, 0.01_dp]
   hess = RESHAPE([2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
   CALL spectral_newton_step(10.0_dp*grad, hess, p, step_cap=0.05_dp)
   IF (ABS(NORM2(p) - 0.05_dp) > 1.0E-14_dp) ERROR STOP "step_cap clamping failed"
   hess = RESHAPE([-2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [3, 3])
   CALL spectral_newton_step([0.02_dp, 0.05_dp, 1.0E-12_dp], hess, p)
   IF (p(1) >= 0.0_dp .OR. ABS(p(3)) > 0.0_dp) ERROR STOP "Damped/vacuum Newton step failed"

END PROGRAM effective_mass_unittest
