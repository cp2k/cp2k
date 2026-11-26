!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Analytical k.p effective-mass reduction kernels and extremum-refinement primitives
!>        in Hartree atomic units (hbar = 1), including two-step Kato-Rellich degenerate
!>        subspace reduction, principal-mass tensor classification, and Brillouin-zone-bounded
!>        spectral Newton stepping.
! **************************************************************************************************
MODULE effective_mass_methods
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE mathlib,                         ONLY: diag_complex,&
                                              diamat_all,&
                                              inv_3x3
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'effective_mass_methods'

   PUBLIC :: kato_rellich_reduce, principal_masses, degenerate_range, sylvester_sqrt_tile
   PUBLIC :: spectral_newton_step, refine_bz_bounds, displace_fractional, wrap_fractional

   !> Carrier classification returned by principal_masses.
   INTEGER, PARAMETER, PUBLIC :: emass_carrier_flat = 0, &
                                 emass_carrier_electron = 1, &
                                 emass_carrier_hole = 2, &
                                 emass_carrier_saddle = 3

   !> Relative curvature threshold used by the carrier-type classification and by the
   !> zero-curvature guard of printed principal masses.
   REAL(dp), PARAMETER, PUBLIC :: carrier_tol = 1.0E-8_dp

   !> Velocity scale (a.u.) below which subspace branches are treated as unsplit along
   !> qhat to first order; it also partitions equal-velocity blocks of the two-step
   !> Kato-Rellich reduction. Symmetry-enforced extrema give ~1e-12.
   REAL(dp), PARAMETER, PRIVATE :: linear_split_tol = 1.0E-7_dp

   !> Refinement: trust-region cap of one Cartesian Newton step (1/bohr) and max BZ fraction.
   REAL(dp), PARAMETER, PRIVATE :: refine_step_cap = 0.20_dp, &
                                   refine_bz_step_frac = 0.35_dp

   !> Refinement: curvature floor of the damped spectral Newton step.
   REAL(dp), PARAMETER, PRIVATE :: spectral_damping = 1.0E-2_dp

   !> Refinement: curvature and projected gradient below which an axis counts as a
   !> dispersionless vacuum direction whose step component stays exactly zero;
   !> vacuum_grad_tol < 1e-7/SQRT(3) guarantees nonzero Newton steps above convergence tol.
   REAL(dp), PARAMETER, PRIVATE :: vacuum_curv_tol = 1.0E-8_dp
   REAL(dp), PARAMETER, PRIVATE :: vacuum_grad_tol = 1.0E-8_dp

   !> Refinement: eigenvector probe displacement for saddle/degenerate escape (1/bohr) and max BZ fraction.
   REAL(dp), PARAMETER, PRIVATE :: refine_probe_step = 5.0E-2_dp, &
                                   refine_bz_probe_frac = 0.10_dp

CONTAINS

! **************************************************************************************************
!> \brief Two-step Kato-Rellich subspace reduction from precomputed nd x nd projected matrices.
!> \param aop Projected velocity operator q_a [H^(a) - E_D S^(a)]_DD
!> \param bop Second-order k.p matrix M_eff in the degenerate subspace D
!> \param sop Projected overlap derivative q_a S^(a)_DD for the moving-block Pulay shift
!> \param velocity_eigs Ascending eigenvalues of aop along q
!> \param inv_mass_eigs Directional inverse masses in 1/m_e, paired with velocity_eigs
!> \param linear_split True when the subspace splits at first order along q
! **************************************************************************************************
   SUBROUTINE kato_rellich_reduce(aop, bop, sop, velocity_eigs, inv_mass_eigs, linear_split)
      COMPLEX(dp), INTENT(IN)                            :: aop(:, :), bop(:, :), sop(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: velocity_eigs(:), inv_mass_eigs(:)
      LOGICAL, INTENT(OUT)                               :: linear_split

      COMPLEX(dp), ALLOCATABLE                           :: btilde(:, :), hblk(:, :), hvec(:, :), &
                                                            stilde(:, :), ueig(:, :)
      INTEGER                                            :: jhi, jlo, m, nd
      REAL(dp)                                           :: vbar

      linear_split = .FALSE.
      nd = SIZE(aop, 1)
      CPASSERT(nd >= 1)
      CPASSERT(ALL(SHAPE(aop) == [nd, nd]))
      CPASSERT(ALL(SHAPE(bop) == [nd, nd]))
      CPASSERT(ALL(SHAPE(sop) == [nd, nd]))

      ALLOCATE (velocity_eigs(nd), inv_mass_eigs(nd), ueig(nd, nd), btilde(nd, nd), stilde(nd, nd))
      CALL diag_complex(aop, ueig, velocity_eigs)

      btilde(:, :) = MATMUL(CONJG(TRANSPOSE(ueig)), MATMUL(bop, ueig))
      stilde(:, :) = MATMUL(CONJG(TRANSPOSE(ueig)), MATMUL(sop, ueig))
      linear_split = (velocity_eigs(nd) - velocity_eigs(1)) > linear_split_tol
      jlo = 1
      DO WHILE (jlo <= nd)
         jhi = jlo
         DO
            IF (jhi >= nd) EXIT
            IF (velocity_eigs(jhi + 1) - velocity_eigs(jhi) > linear_split_tol) EXIT
            jhi = jhi + 1
         END DO
         vbar = SUM(velocity_eigs(jlo:jhi))/REAL(jhi - jlo + 1, dp)
         m = jhi - jlo + 1
         ALLOCATE (hblk(m, m), hvec(m, m))
         hblk(:, :) = btilde(jlo:jhi, jlo:jhi) - 2.0_dp*vbar*stilde(jlo:jhi, jlo:jhi)
         CALL diag_complex(hblk, hvec, inv_mass_eigs(jlo:jhi))
         DEALLOCATE (hblk, hvec)
         jlo = jhi + 1
      END DO
   END SUBROUTINE kato_rellich_reduce

! **************************************************************************************************
!> \brief Contiguous degenerate subspace of a band on the ascending eigenvalue list.
!> \param energies Ascending generalized eigenvalues
!> \param iband Band index anchoring the subspace
!> \param eps_degeneracy Adjacent-gap criterion (Ha); gaps <= eps_degeneracy join bands
!> \param ilo Lowest band index of the subspace
!> \param ihi Highest band index of the subspace
!> \note Clustering on adjacent gaps is transitive, unlike an absolute window around
!>       energies(iband), so every member of a quasi-degenerate multiplet finds the same
!>       subspace independent of the anchor member.
! **************************************************************************************************
   SUBROUTINE degenerate_range(energies, iband, eps_degeneracy, ilo, ihi)
      REAL(dp), INTENT(IN)                               :: energies(:)
      INTEGER, INTENT(IN)                                :: iband
      REAL(dp), INTENT(IN)                               :: eps_degeneracy
      INTEGER, INTENT(OUT)                               :: ilo, ihi

      ilo = iband
      DO
         IF (ilo <= 1) EXIT
         IF (energies(ilo) - energies(ilo - 1) > eps_degeneracy) EXIT
         ilo = ilo - 1
      END DO
      ihi = iband
      DO
         IF (ihi >= SIZE(energies)) EXIT
         IF (energies(ihi + 1) - energies(ihi) > eps_degeneracy) EXIT
         ihi = ihi + 1
      END DO
   END SUBROUTINE degenerate_range

! **************************************************************************************************
!> \brief Diagonalize a 3x3 inverse mass tensor into principal masses and carrier type.
!> \param inv_mass Symmetric inverse mass tensor in 1/m_e
!> \param inv_mass_eigs Ascending principal inverse masses in 1/m_e
!> \param principal_dirs Principal directions as eigenvector columns (sign gauge arbitrary)
!> \param carrier One of emass_carrier_electron/hole/saddle/flat; semi-definite tensors of
!>        2D/1D systems (vacuum axes with vanishing inverse mass) follow the sign of the
!>        remaining nonzero components
! **************************************************************************************************
   SUBROUTINE principal_masses(inv_mass, inv_mass_eigs, principal_dirs, carrier)
      REAL(dp), INTENT(IN)                               :: inv_mass(3, 3)
      REAL(dp), INTENT(OUT)                              :: inv_mass_eigs(3), principal_dirs(3, 3)
      INTEGER, INTENT(OUT)                               :: carrier

      INTEGER                                            :: n_neg, n_pos
      REAL(dp)                                           :: a(3, 3), tol

      a = 0.5_dp*(inv_mass + TRANSPOSE(inv_mass))
      CALL diamat_all(a, inv_mass_eigs)
      principal_dirs = a
      tol = carrier_tol*MAX(1.0_dp, MAXVAL(ABS(inv_mass_eigs)))
      n_pos = COUNT(inv_mass_eigs > tol)
      n_neg = COUNT(inv_mass_eigs < -tol)
      IF (n_pos == 0 .AND. n_neg == 0) THEN
         carrier = emass_carrier_flat
      ELSE IF (n_pos > 0 .AND. n_neg == 0) THEN
         carrier = emass_carrier_electron
      ELSE IF (n_neg > 0 .AND. n_pos == 0) THEN
         carrier = emass_carrier_hole
      ELSE
         carrier = emass_carrier_saddle
      END IF
   END SUBROUTINE principal_masses

! **************************************************************************************************
!> \brief Damped spectral Newton step of the band-extremum refinement, in Cartesian
!>        coordinates and free of any Quickstep dependency.
!> \param gradient Cartesian gradient g = s * velocity of the target band, s = +-1
!> \param hessian Cartesian Hessian H = s * inverse mass tensor of the target band
!> \param step Capped Newton step (1/bohr), with exactly zero components along
!>        dispersionless vacuum axes
!> \param step_cap Optional Cartesian trust-region radius (1/bohr); defaults to refine_step_cap
!> \note The Hessian eigenbasis gives the step p = -V diag(1/max(|lambda_i|, 1e-2)) V^T g;
!>       inside a convex basin (|lambda_i| > 1e-2 everywhere) this is the exact Newton
!>       step with quadratic convergence, while the floor regularizes indefinite or flat
!>       curvature and the norm cap keeps the step within the trust-region radius.
! **************************************************************************************************
   SUBROUTINE spectral_newton_step(gradient, hessian, step, step_cap)
      REAL(dp), INTENT(IN)                               :: gradient(3), hessian(3, 3)
      REAL(dp), INTENT(OUT)                              :: step(3)
      REAL(dp), INTENT(IN), OPTIONAL                     :: step_cap

      REAL(dp)                                           :: a(3, 3), cap, evals(3), gtilde(3), &
                                                            ptilde(3)

      cap = refine_step_cap
      IF (PRESENT(step_cap)) cap = step_cap
      a = 0.5_dp*(hessian + TRANSPOSE(hessian))
      CALL diamat_all(a, evals)
      gtilde = MATMUL(gradient, a)
      ptilde = MERGE(0.0_dp, -gtilde/MAX(ABS(evals), spectral_damping), &
                     ABS(evals) <= vacuum_curv_tol .AND. ABS(gtilde) <= vacuum_grad_tol)
      step = MATMUL(a, ptilde)
      IF (NORM2(step) > cap) step = step*(cap/NORM2(step))
   END SUBROUTINE spectral_newton_step

! **************************************************************************************************
!> \brief Wrap fractional k-point components into the home cell [-0.5, 0.5).
!> \param kf Fractional coordinates, arbitrary periodic image
!> \return Wrapped fractional coordinates
! **************************************************************************************************
   PURE FUNCTION wrap_fractional(kf) RESULT(kw)
      REAL(dp), INTENT(IN)                               :: kf(3)
      REAL(dp)                                           :: kw(3)

      kw = kf - ANINT(kf)
      WHERE (kw >= 0.5_dp) kw = kw - 1.0_dp
   END FUNCTION wrap_fractional

! **************************************************************************************************
!> \brief Compute Brillouin-zone-bounded step cap and probe displacement for refinement.
!> \param hmat Cell vectors in bohr, columns a_i
!> \param periodic Periodicity mask per lattice vector
!> \param step_cap Capped Cartesian Newton step radius (1/bohr)
!> \param probe_step Capped Cartesian saddle/degenerate probe step (1/bohr)
! **************************************************************************************************
   SUBROUTINE refine_bz_bounds(hmat, periodic, step_cap, probe_step)
      REAL(dp), INTENT(IN)                               :: hmat(3, 3)
      LOGICAL, INTENT(IN)                                :: periodic(3)
      REAL(dp), INTENT(OUT)                              :: step_cap, probe_step

      INTEGER                                            :: i
      REAL(dp)                                           :: b_min, bmat(3, 3)

      IF (.NOT. ANY(periodic)) THEN
         step_cap = refine_step_cap
         probe_step = refine_probe_step
         RETURN
      END IF
      bmat(:, :) = twopi*inv_3x3(hmat)
      b_min = HUGE(1.0_dp)
      DO i = 1, 3
         IF (periodic(i)) b_min = MIN(b_min, NORM2(bmat(i, 1:3)))
      END DO
      step_cap = MIN(refine_step_cap, refine_bz_step_frac*b_min)
      probe_step = MIN(refine_probe_step, refine_bz_probe_frac*b_min)
   END SUBROUTINE refine_bz_bounds

! **************************************************************************************************
!> \brief Shift fractional wave vector kf0 by Cartesian step step_cart, holding non-periodic
!>        directions fixed, and wrap into [-0.5, 0.5)^3.
!> \param kf0 Starting fractional wave vector
!> \param step_cart Cartesian k-displacement in 1/bohr
!> \param hmat Cell vectors in bohr, columns a_i
!> \param periodic Periodicity mask per lattice vector
!> \return Wrapped displaced fractional wave vector
! **************************************************************************************************
   PURE FUNCTION displace_fractional(kf0, step_cart, hmat, periodic) RESULT(kf1)
      REAL(dp), INTENT(IN)                               :: kf0(3), step_cart(3), hmat(3, 3)
      LOGICAL, INTENT(IN)                                :: periodic(3)
      REAL(dp)                                           :: kf1(3)

      REAL(dp)                                           :: dk(3)

      dk = MATMUL(step_cart, hmat)/twopi
      WHERE (.NOT. periodic) dk = 0.0_dp
      kf1 = wrap_fractional(kf0 + dk)
   END FUNCTION displace_fractional

! **************************************************************************************************
!> \brief Solve X diag(sigma_col) + diag(sigma_row) X = R elementwise in place on a local tile.
!> \param tile Complex matrix block R, overwritten with X_ij = R_ij / (sigma_row(i) + sigma_col(j))
!> \param sigma_row Positive square-root overlap eigenvalues of the block rows
!> \param sigma_col Positive square-root overlap eigenvalues of the block columns
! **************************************************************************************************
   PURE SUBROUTINE sylvester_sqrt_tile(tile, sigma_row, sigma_col)
      COMPLEX(dp), INTENT(INOUT)                         :: tile(:, :)
      REAL(dp), INTENT(IN)                               :: sigma_row(:), sigma_col(:)

      INTEGER                                            :: i, j

      DO j = 1, SIZE(tile, 2)
         DO i = 1, SIZE(tile, 1)
            tile(i, j) = tile(i, j)/(sigma_row(i) + sigma_col(j))
         END DO
      END DO
   END SUBROUTINE sylvester_sqrt_tile

END MODULE effective_mass_methods
