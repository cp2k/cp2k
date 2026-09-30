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
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: get_nlcc_gaussians

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
            alpha = gth_alpha(1:nexp)
            nct = gth_nct(1:nexp)
            cval = gth_cval(:, 1:nexp)
            RETURN
         END IF
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         CALL get_potential(sgp_potential, has_nlcc=nlcc, n_nlcc=nexp, &
                            a_nlcc=a_nlcc, c_nlcc=c_nlcc)
         IF (nlcc) THEN
            CPASSERT(ALL(a_nlcc(1:nexp) > 0.0_dp))
            ! SGP stores exponents in exp(-a*r**2), whereas GTH stores widths.
            alpha = 1.0_dp/SQRT(2.0_dp*a_nlcc(1:nexp))
            ALLOCATE (nct(nexp), cval(1, nexp))
            nct = 1
            cval(1, :) = c_nlcc(1:nexp)
            RETURN
         END IF
      END IF
      ALLOCATE (alpha(0), nct(0), cval(0, 0))

   END SUBROUTINE get_nlcc_gaussians

END MODULE qs_nlcc_methods
