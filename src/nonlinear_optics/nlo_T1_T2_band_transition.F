!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Smooth band transition between population and coherence relaxation times.
! **************************************************************************************************
MODULE nlo_T1_T2_band_transition
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: gaussi,&
                                              pi
   USE nlo_derivatives,                 ONLY: derivatives_index,&
                                              nlo_derivative_count,&
                                              nlo_derivative_indices,&
                                              nlo_n_derivatives
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo_T1_T2_band_transition'

   PUBLIC :: nlo_band_transition_type, nlo_band_transition_prepare, &
             nlo_band_transition_correction, &
             nlo_band_transition_rate

   TYPE :: nlo_band_transition_type
      REAL(dp), ALLOCATABLE     :: energies(:), pair_rate(:, :)
      REAL(dp)                  :: width = 0.0_dp, rate_difference = 0.0_dp
      COMPLEX(dp), ALLOCATABLE  :: integration_points(:), integration_weights(:)
   END TYPE nlo_band_transition_type

CONTAINS

! **************************************************************************************************
!> \brief Prepare Γ_H(X)=γ₂X+(γ₁−γ₂)exp[−(ad_H/w)²]X, ad_H(X)=[H,X].
!> \param energies Band energies in Hartree
!> \param width Energy scale of the smooth transition in Hartree
!> \param rate1 Population damping energy in Hartree
!> \param rate2 Coherence damping energy in Hartree
!> \param relaxation Prepared band transition state
!> \note Central rates use EXP exactly. Positive grades use Cauchy's formula
!>       on two horizontal contours and Gauss-Hermite integration. The contour
!>       shifts and factors are generated from the Hermite Jacobi matrix; no
!>       fitted or tabulated coefficients enter the relaxation model.
!>       Taylor derivatives follow analytically from the resolvent equation.
!>       Contour points share one workspace; lower grades are reconstructed
!>       once per requested grade from the already solved physical density.
! **************************************************************************************************
   SUBROUTINE nlo_band_transition_prepare(energies, width, rate1, rate2, relaxation)
      REAL(dp), INTENT(IN)                               :: energies(:), width, rate1, rate2
      TYPE(nlo_band_transition_type), INTENT(INOUT)      :: relaxation

      INTEGER                                            :: iband, jband, nband
      REAL(dp)                                           :: difference, gaussian

      IF (.NOT. ALL(ieee_is_finite(energies)) .OR. .NOT. ieee_is_finite(width) .OR. &
          .NOT. ieee_is_finite(rate1) .OR. .NOT. ieee_is_finite(rate2) .OR. &
          width <= 0.0_dp .OR. rate1 <= 0.0_dp .OR. rate2 < rate1) THEN
         CPABORT('Invalid NLO band-transition relaxation parameters.')
      END IF
      IF (ALLOCATED(relaxation%pair_rate)) DEALLOCATE (relaxation%pair_rate)
      IF (ALLOCATED(relaxation%energies)) DEALLOCATE (relaxation%energies)
      ALLOCATE (relaxation%energies(SIZE(energies)))
      relaxation%energies(:) = energies
      relaxation%width = width
      relaxation%rate_difference = rate1 - rate2
      CALL nlo_gaussian_quadrature(relaxation%integration_points, relaxation%integration_weights)
      nband = SIZE(energies)
      CPASSERT(nband >= 2)
      ALLOCATE (relaxation%pair_rate(nband, nband))
      DO jband = 1, nband
         DO iband = 1, nband
            difference = energies(iband) - energies(jband)
            gaussian = 0.0_dp
            IF (width >= ABS(difference)/26.0_dp) gaussian = EXP(-(difference/width)**2)
            relaxation%pair_rate(iband, jband) = rate1 + (rate2 - rate1)*(1.0_dp - gaussian)
         END DO
      END DO
   END SUBROUTINE nlo_band_transition_prepare

! **************************************************************************************************
!> \brief Return the central relaxation rate for one band pair.
!> \param relaxation Prepared band transition state
!> \param iband First band index
!> \param jband Second band index
!> \return Central damping rate of the selected band pair
! **************************************************************************************************
   FUNCTION nlo_band_transition_rate(relaxation, iband, jband) RESULT(rate)
      TYPE(nlo_band_transition_type), INTENT(IN)         :: relaxation
      INTEGER, INTENT(IN)                                :: iband, jband
      REAL(dp)                                           :: rate

      CPASSERT(ALLOCATED(relaxation%pair_rate))
      CPASSERT(iband >= 1)
      CPASSERT(iband <= SIZE(relaxation%pair_rate, 1))
      CPASSERT(jband >= 1)
      CPASSERT(jband <= SIZE(relaxation%pair_rate, 2))
      rate = relaxation%pair_rate(iband, jband)
   END FUNCTION nlo_band_transition_rate

! **************************************************************************************************
!> \brief Sum relaxation corrections for one complete Cartesian Taylor grade.
!> \param relaxation Prepared band transition state
!> \param hamiltonian Hamiltonian Taylor coefficients in the band basis
!> \param density Physical Taylor coefficients; only lower grades are read
!> \param grade Requested total Cartesian derivative order, between one and four
!> \param correction Corrections of this grade in nlo_derivative_indices order
!> \note For Y_l=(ad_H−w z_l)⁻¹ρ and R_l,α=−Σ_{0<β≤α}[H_β,Y_l,α−β],
!>       Y_l,α=(ρ_α+R_l,α)/(Δε−w z_l). Only grades below the requested grade
!>       are needed in R_l,α. Reconstruct them for one contour point, accumulate
!>       all corrections of the requested grade, then reuse the workspace.
!>       No band-sized array carries a contour-point dimension. At grade four,
!>       20 lower-grade matrices replace 35 matrices for each of 96 points.
! **************************************************************************************************
   SUBROUTINE nlo_band_transition_correction(relaxation, hamiltonian, density, grade, correction)
      TYPE(nlo_band_transition_type), INTENT(IN)         :: relaxation
      COMPLEX(dp), INTENT(IN)                            :: hamiltonian(:, :, :), density(:, :, :)
      INTEGER, INTENT(IN)                                :: grade
      COMPLEX(dp), INTENT(OUT)                           :: correction(:, :, :)

      COMPLEX(dp)                                        :: weight
      COMPLEX(dp), ALLOCATABLE                           :: inverse(:, :), lower_density(:, :, :), &
                                                            matrix_term(:, :), remainder(:, :)
      INTEGER :: alpha(3, nlo_n_derivatives), ia, iband, ibeta, igamma, ipoint, jband, nband, &
         ncoeff, nlower, partner(nlo_n_derivatives, nlo_n_derivatives)

      CPASSERT(grade >= 1 .AND. grade <= 4)
      CPASSERT(ALLOCATED(relaxation%pair_rate))
      nband = SIZE(relaxation%pair_rate, 1)
      nlower = nlo_derivative_count(grade - 1)
      ncoeff = nlo_derivative_count(grade)
      CPASSERT(SIZE(hamiltonian, 1) == nband)
      CPASSERT(SIZE(hamiltonian, 2) == nband)
      CPASSERT(SIZE(hamiltonian, 3) >= ncoeff)
      CPASSERT(SIZE(density, 1) == nband)
      CPASSERT(SIZE(density, 2) == nband)
      CPASSERT(SIZE(density, 3) >= nlower)
      CPASSERT(SIZE(correction, 1) == nband)
      CPASSERT(SIZE(correction, 2) == nband)
      CPASSERT(SIZE(correction, 3) == ncoeff - nlower)
      ALLOCATE (inverse(nband, nband), lower_density(nband, nband, nlower), &
                remainder(nband, nband), matrix_term(nband, nband))
      CALL nlo_derivative_indices(alpha)
      partner = 0
      DO ia = 2, ncoeff
         DO ibeta = 2, ia
            partner(ibeta, ia) = derivatives_index(alpha, alpha(:, ia) - alpha(:, ibeta))
         END DO
      END DO
      correction = CMPLX(0.0_dp, 0.0_dp, dp)
      DO ipoint = 1, SIZE(relaxation%integration_points)
         DO jband = 1, nband
            DO iband = 1, nband
               inverse(iband, jband) = 1.0_dp/ &
                                       (relaxation%energies(iband) - relaxation%energies(jband) - &
                                        relaxation%width*relaxation%integration_points(ipoint))
            END DO
         END DO
         weight = relaxation%width*relaxation%integration_weights(ipoint)
         DO ia = 1, ncoeff
            remainder = CMPLX(0.0_dp, 0.0_dp, dp)
            DO ibeta = 2, ia
               igamma = partner(ibeta, ia)
               IF (igamma == 0) CYCLE
               matrix_term(:, :) = MATMUL(hamiltonian(:, :, ibeta), lower_density(:, :, igamma))
               remainder(:, :) = remainder - matrix_term
               matrix_term(:, :) = MATMUL(lower_density(:, :, igamma), hamiltonian(:, :, ibeta))
               remainder(:, :) = remainder + matrix_term
            END DO
            remainder(:, :) = remainder*inverse
            IF (ia <= nlower) THEN
               lower_density(:, :, ia) = density(:, :, ia)*inverse + remainder
            ELSE
               correction(:, :, ia - nlower) = correction(:, :, ia - nlower) + weight*remainder
            END IF
         END DO
      END DO
      correction = CMPLX(relaxation%rate_difference, 0.0_dp, dp)*correction
   END SUBROUTINE nlo_band_transition_correction

! **************************************************************************************************
!> \brief Generate integration points and weights for the Gaussian matrix function.
!> \param points Complex contour integration points
!> \param weights Weights for the complex contour integration
!> \note f(x)=exp(-x²) equals the difference of Cauchy integrals on the
!>       lines Im(z)=±h. Along each line, f(x±ih)=exp(h²)exp(-x²)exp(∓2ixh).
!>       Eigenvalues of the Hermite Jacobi matrix are quadrature nodes; the
!>       squared first-row eigenvector entries determine quadrature weights.
!>       The numerical rule uses 48 nodes per line at height h=2.5. These
!>       choices resolve f and its first four derivatives to about 2e-14
!>       scaled absolute error in checks over −10≤x≤10; they are not model
!>       parameters. Points and weights belong to the calling transition state.
! **************************************************************************************************
   SUBROUTINE nlo_gaussian_quadrature(points, weights)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: points(:), weights(:)

      INTEGER, PARAMETER                                 :: n_hermite = 48
      REAL(dp), PARAMETER                                :: height = 2.5_dp

      COMPLEX(dp)                                        :: factor
      INTEGER                                            :: info, ipoint
      REAL(dp)                                           :: jacobi(n_hermite, n_hermite), &
                                                            nodes(n_hermite), quadrature_weight, &
                                                            work(3*n_hermite), x

      ALLOCATE (points(2*n_hermite), weights(2*n_hermite))
      jacobi = 0.0_dp
      DO ipoint = 1, n_hermite - 1
         jacobi(ipoint, ipoint + 1) = SQRT(REAL(ipoint, dp)/2.0_dp)
         jacobi(ipoint + 1, ipoint) = jacobi(ipoint, ipoint + 1)
      END DO
      CALL dsyev('V', 'U', n_hermite, jacobi, n_hermite, nodes, work, SIZE(work), info)
      IF (info /= 0) CPABORT('NLO Gaussian contour quadrature failed.')
      DO ipoint = 1, n_hermite
         x = nodes(ipoint)
         quadrature_weight = SQRT(pi)*jacobi(1, ipoint)**2
         factor = quadrature_weight*EXP(height**2)/(2.0_dp*pi*gaussi)
         points(ipoint) = CMPLX(x, -height, dp)
         points(ipoint + n_hermite) = CMPLX(x, height, dp)
         weights(ipoint) = -factor*EXP(2.0_dp*gaussi*x*height)
         weights(ipoint + n_hermite) = factor*EXP(-2.0_dp*gaussi*x*height)
      END DO
   END SUBROUTINE nlo_gaussian_quadrature

END MODULE nlo_T1_T2_band_transition
