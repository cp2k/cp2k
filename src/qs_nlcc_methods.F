!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common Gaussian representation of nonlinear core densities.
! **************************************************************************************************
MODULE qs_nlcc_methods
   USE cell_types,                      ONLY: cell_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   TYPE nlcc_atom_grid_type
      REAL(dp), ALLOCATABLE :: exponent(:), coefficient(:, :), radius(:)
   END TYPE nlcc_atom_grid_type

   PUBLIC :: get_nlcc_gaussians, nlcc_atom_grid_type, nlcc_atom_grid_init, nlcc_atom_grid_evaluate

CONTAINS

! **************************************************************************************************
!> \brief Return rho(r) = sum_i exp(-r**2/(2*alpha(i)**2)) *
!>        sum_j cval(j,i)*(r/alpha(i))**(2*j-2).
!>        The same representation is used for collocation and force integration.
!> \param qs_kind Atomic kind carrying a GTH or separable Gaussian potential.
!> \param alpha Gaussian widths; an empty array indicates no core correction.
!> \param nct Number of radial polynomial coefficients for each Gaussian.
!> \param cval Radial polynomial coefficients, including the density normalization.
! **************************************************************************************************
   SUBROUTINE get_nlcc_gaussians(qs_kind, alpha, nct, cval)
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: alpha
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: nct
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: cval

      INTEGER                                            :: nexp
      INTEGER, DIMENSION(:), POINTER                     :: gth_nct
      LOGICAL                                            :: nlcc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_nlcc, c_nlcc, gth_alpha
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gth_cval
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL get_qs_kind(qs_kind, gth_potential=gth_potential, sgp_potential=sgp_potential)
      IF (ASSOCIATED(gth_potential)) THEN
         CALL get_potential(gth_potential, nlcc_present=nlcc, nexp_nlcc=nexp, &
                            alpha_nlcc=gth_alpha, nct_nlcc=gth_nct, cval_nlcc=gth_cval)
         IF (nlcc) THEN
            ALLOCATE (alpha(nexp), nct(nexp), cval(SIZE(gth_cval, 1), nexp))
            alpha(:) = gth_alpha(1:nexp)
            nct(:) = gth_nct(1:nexp)
            cval(:, :) = gth_cval(:, 1:nexp)
            RETURN
         END IF
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         CALL get_potential(sgp_potential, has_nlcc=nlcc, n_nlcc=nexp, &
                            a_nlcc=a_nlcc, c_nlcc=c_nlcc)
         IF (nlcc) THEN
            CPASSERT(ALL(a_nlcc(1:nexp) > 0.0_dp))
            ! SGP stores exponents in exp(-a*r**2), whereas GTH stores widths.
            ALLOCATE (alpha(nexp), nct(nexp), cval(1, nexp))
            alpha(:) = 1.0_dp/SQRT(2.0_dp*a_nlcc(1:nexp))
            nct = 1
            cval(1, :) = c_nlcc(1:nexp)
            RETURN
         END IF
      END IF
      ALLOCATE (alpha(0), nct(0), cval(0, 0))

   END SUBROUTINE get_nlcc_gaussians

! **************************************************************************************************
!> \brief Cache core Gaussian polynomials and conservative density/gradient/Hessian tail radii.
!> \param qs_kind Atomic kind supplying the core density.
!> \param DATA Reusable atom-grid core representation.
!> \param eps Absolute tail threshold in atomic units, divided among the Gaussian terms.
! **************************************************************************************************
   SUBROUTINE nlcc_atom_grid_init(qs_kind, DATA, eps)
      TYPE(qs_kind_type), INTENT(IN)                     :: qs_kind
      TYPE(nlcc_atom_grid_type), INTENT(OUT)             :: data
      REAL(dp), INTENT(IN)                               :: eps

      INTEGER                                            :: i, j, n
      INTEGER, ALLOCATABLE                               :: nct(:)
      REAL(dp)                                           :: bound, radius
      REAL(dp), ALLOCATABLE                              :: alpha(:), cval(:, :)

      CPASSERT(eps > 0.0_dp)
      CALL get_nlcc_gaussians(qs_kind, alpha, nct, cval)
      n = SIZE(alpha)
      ALLOCATE (data%exponent(n), data%coefficient(SIZE(cval, 1), n), data%radius(n))
      data%coefficient = 0.0_dp
      data%radius = 0.0_dp
      DO i = 1, n
         data%exponent(i) = 0.5_dp/alpha(i)**2
         DO j = 1, nct(i)
            data%coefficient(j, i) = cval(j, i)/alpha(i)**(2*j - 2)
         END DO
         IF (ALL(data%coefficient(:, i) == 0.0_dp)) CYCLE
         ! Beyond this radius every monomial in the absolute derivative bound
         ! decreases. The support depends only on the supplied core representation.
         radius = SQRT(REAL(nct(i), dp)/data%exponent(i))
         DO
            bound = nlcc_tail_bound(data%coefficient(:, i), data%exponent(i), radius)
            IF (bound <= eps/REAL(n, dp)) EXIT
            radius = 1.2_dp*radius
         END DO
         data%radius(i) = radius
      END DO
   END SUBROUTINE nlcc_atom_grid_init

! **************************************************************************************************
!> \brief Bound Cartesian values and first/second derivatives of one radial Gaussian polynomial.
!> \param coefficient Coefficients of powers of r**2.
!> \param exponent Gaussian exponent.
!> \param radius Positive radius on the decreasing tail.
!> \return Maximum absolute density, gradient-component and Hessian-component bound.
! **************************************************************************************************
   FUNCTION nlcc_tail_bound(coefficient, exponent, radius) RESULT(bound)
      REAL(dp), INTENT(IN)                               :: coefficient(:), exponent, radius
      REAL(dp)                                           :: bound

      INTEGER                                            :: j, p
      REAL(dp)                                           :: ddf, df, f, r2

      r2 = radius**2
      f = 0.0_dp
      df = 0.0_dp
      ddf = 0.0_dp
      DO j = 1, SIZE(coefficient)
         p = j - 1
         f = f + ABS(coefficient(j))*r2**p
         IF (p > 0) df = df + REAL(p, dp)*ABS(coefficient(j))*r2**(p - 1)
         IF (p > 1) ddf = ddf + REAL(p*(p - 1), dp)*ABS(coefficient(j))*r2**(p - 2)
      END DO
      bound = EXP(-exponent*r2)*MAX(f, 2.0_dp*radius*(df + exponent*f), &
                                    2.0_dp*(df + exponent*f) + &
                                    4.0_dp*r2*(ddf + 2.0_dp*exponent*df + exponent**2*f))
   END FUNCTION nlcc_tail_bound

! **************************************************************************************************
!> \brief Sum a core density and Cartesian derivatives over all supported periodic images.
!> \param DATA Prepared core Gaussians.
!> \param cell Cell matrix and inverse; the periodicity argument selects auxiliary-field images.
!> \param periodicity One for each periodic direction, zero otherwise.
!> \param point Cartesian evaluation point, with arbitrary cell representative.
!> \param center Cartesian core center.
!> \param rho Full core density before its equal spin partition.
!> \param gradient Cartesian gradient.
!> \param hessian Cartesian Hessian.
!> \param rho_image_moment Optional sum of gradient(alpha)*image_translation(beta).
!> \param gradient_image_moment Optional sum of Hessian(gamma,alpha)*image_translation(beta).
!>        Contracting these positive moments with field adjoints gives the image virial;
!>        their negatives are cell-strain derivatives at fixed Cartesian point and center.
! **************************************************************************************************
   SUBROUTINE nlcc_atom_grid_evaluate(DATA, cell, periodicity, point, center, rho, gradient, hessian, &
                                      rho_image_moment, gradient_image_moment)
      TYPE(nlcc_atom_grid_type), INTENT(IN)              :: data
      TYPE(cell_type), INTENT(IN)                        :: cell
      INTEGER, INTENT(IN)                                :: periodicity(3)
      REAL(dp), INTENT(IN)                               :: point(3), center(3)
      REAL(dp), INTENT(OUT)                              :: rho, gradient(3), hessian(3, 3)
      REAL(dp), INTENT(OUT), OPTIONAL                    :: rho_image_moment(3, 3), &
                                                            gradient_image_moment(3, 3, 3)

      INTEGER                                            :: axis, axis2, i, ix, iy, iz, j, lower(3), &
                                                            p, upper(3)
      REAL(dp)                                           :: a, ddf, df, f, first, fractional(3), &
                                                            grad(3), hess(3, 3), padding, r2, &
                                                            radius, second, shift(3), value, x(3)

      CPASSERT(ALLOCATED(data%radius))
      rho = 0.0_dp
      gradient = 0.0_dp
      hessian = 0.0_dp
      IF (PRESENT(rho_image_moment)) rho_image_moment = 0.0_dp
      IF (PRESENT(gradient_image_moment)) gradient_image_moment = 0.0_dp
      fractional = MATMUL(cell%h_inv, point - center)
      DO i = 1, SIZE(data%radius)
         IF (data%radius(i) == 0.0_dp) CYCLE
         lower = 0
         upper = 0
         DO axis = 1, 3
            IF (periodicity(axis) == 0) CYCLE
            radius = data%radius(i)*SQRT(SUM(cell%h_inv(axis, :)**2))
            padding = 64.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, radius, ABS(fractional(axis)))
            lower(axis) = CEILING(fractional(axis) - radius - padding)
            upper(axis) = FLOOR(fractional(axis) + radius + padding)
         END DO
         a = data%exponent(i)
         DO iz = lower(3), upper(3)
            DO iy = lower(2), upper(2)
               DO ix = lower(1), upper(1)
                  shift = MATMUL(cell%hmat, REAL([ix, iy, iz], dp))
                  x = point - center - shift
                  r2 = SUM(x**2)
                  IF (r2 > data%radius(i)**2) CYCLE
                  f = 0.0_dp
                  df = 0.0_dp
                  ddf = 0.0_dp
                  DO j = 1, SIZE(data%coefficient, 1)
                     p = j - 1
                     f = f + data%coefficient(j, i)*r2**p
                     IF (p > 0) df = df + REAL(p, dp)*data%coefficient(j, i)*r2**(p - 1)
                     IF (p > 1) ddf = ddf + REAL(p*(p - 1), dp)*data%coefficient(j, i)*r2**(p - 2)
                  END DO
                  value = EXP(-a*r2)
                  second = value*(ddf - 2.0_dp*a*df + a*a*f)
                  first = value*(df - a*f)
                  rho = rho + value*f
                  grad = 2.0_dp*first*x
                  DO axis = 1, 3
                     hess(:, axis) = 4.0_dp*second*x*x(axis)
                     hess(axis, axis) = hess(axis, axis) + 2.0_dp*first
                  END DO
                  gradient = gradient + grad
                  hessian = hessian + hess
                  IF (PRESENT(rho_image_moment)) THEN
                     DO axis = 1, 3
                        rho_image_moment(:, axis) = rho_image_moment(:, axis) + grad*shift(axis)
                     END DO
                  END IF
                  IF (PRESENT(gradient_image_moment)) THEN
                     DO axis2 = 1, 3
                        DO axis = 1, 3
                           gradient_image_moment(:, axis, axis2) = &
                              gradient_image_moment(:, axis, axis2) + hess(:, axis)*shift(axis2)
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE nlcc_atom_grid_evaluate

END MODULE qs_nlcc_methods
