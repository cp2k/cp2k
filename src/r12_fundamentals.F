!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Radial fundamentals for Gaussian two-electron integrals.
!>        Entry m+1 stores order m. Private kernels provide unnormalized coefficients;
!>        the dispatcher returns their common multiplicative normalization separately.
!>        Gaussian-product attenuation and geometry prefactors belong to the caller.
! **************************************************************************************************
MODULE r12_fundamentals
   USE gamma,                           ONLY: fgamma => fgamma_0
   USE input_constants,                 ONLY: &
        do_potential_TShPSC, do_potential_coulomb, do_potential_gaussian, do_potential_long, &
        do_potential_mix_cl, do_potential_mix_cl_trunc, do_potential_mix_lg, do_potential_short, &
        do_potential_truncated
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE t_c_g0,                          ONLY: t_c_g0_n
   USE t_sh_p_s_c,                      ONLY: trunc_CS_poly_n20
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: compute_r12_fundamentals

CONTAINS

! **************************************************************************************************
!> \brief Compute the unnormalized Coulomb Boys coefficients.
!> \param F coefficients of orders zero through m_max
!> \param m_max maximum order, nonnegative
!> \param T reduced squared distance between Gaussian product centers
! **************************************************************************************************
   SUBROUTINE compute_coulomb_fundamentals(F, m_max, T)
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      REAL(dp), INTENT(IN)                               :: T

      CALL fgamma(m_max, T, F)
   END SUBROUTINE compute_coulomb_fundamentals

! **************************************************************************************************
!> \brief Compute radial coefficients for erf(omega*r)/r.
!> \param F coefficients of orders zero through m_max
!> \param m_max maximum order, nonnegative
!> \param rho reduced Gaussian exponent, positive
!> \param T reduced squared distance between Gaussian product centers
!> \param omega range-separation parameter in inverse length units
! **************************************************************************************************
   SUBROUTINE compute_long_range_fundamentals(F, m_max, rho, T, omega)
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      REAL(dp), INTENT(IN)                               :: rho, T, omega

      INTEGER                                            :: m
      REAL(dp)                                           :: omega2, omega_corr2, scaled_T, weight

      omega2 = omega**2
      omega_corr2 = omega2/(omega2 + rho)
      scaled_T = T*omega_corr2
      CALL fgamma(m_max, scaled_T, F)
      weight = SQRT(omega_corr2)
      DO m = 1, m_max + 1
         F(m) = F(m)*weight
         weight = weight*omega_corr2
      END DO
   END SUBROUTINE compute_long_range_fundamentals

! **************************************************************************************************
!> \brief Compute radial coefficients for erfc(omega*r)/r.
!> \param F coefficients of orders zero through m_max
!> \param m_max maximum order, nonnegative
!> \param rho reduced Gaussian exponent, positive
!> \param T reduced squared distance between Gaussian product centers
!> \param omega range-separation parameter in inverse length units
!> \param work scratch space for the long-range coefficients
! **************************************************************************************************
   SUBROUTINE compute_short_range_fundamentals(F, m_max, rho, T, omega, work)
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      REAL(dp), INTENT(IN)                               :: rho, T, omega
      REAL(dp), INTENT(INOUT)                            :: work(m_max + 1)

      CALL compute_coulomb_fundamentals(F, m_max, T)
      CALL compute_long_range_fundamentals(work, m_max, rho, T, omega)
      F = F - work
   END SUBROUTINE compute_short_range_fundamentals

! **************************************************************************************************
!> \brief Compute radial coefficients for Coulomb truncated at cutoff_radius.
!> \param F coefficients of orders zero through m_max
!> \param m_max maximum order, nonnegative
!> \param rho reduced Gaussian exponent, positive
!> \param T reduced squared distance between Gaussian product centers
!> \param cutoff_radius truncation radius
!> \note The truncated-Coulomb tables must be initialized before this call.
! **************************************************************************************************
   SUBROUTINE compute_truncated_fundamentals(F, m_max, rho, T, cutoff_radius)
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      REAL(dp), INTENT(IN)                               :: rho, T, cutoff_radius

      LOGICAL                                            :: use_gamma
      REAL(dp)                                           :: R

      R = cutoff_radius*SQRT(rho)
      CALL t_c_g0_n(F, use_gamma, R, T, m_max)
      IF (use_gamma) CALL compute_coulomb_fundamentals(F, m_max, T)
   END SUBROUTINE compute_truncated_fundamentals

! **************************************************************************************************
!> \brief Compute radial coefficients for exp(-alpha*r*r).
!> \param F coefficients of orders zero through m_max
!> \param m_max maximum order, nonnegative
!> \param rho reduced Gaussian exponent, positive
!> \param T reduced squared distance between Gaussian product centers
!> \param alpha Gaussian operator exponent in inverse squared length units
! **************************************************************************************************
   SUBROUTINE compute_gaussian_fundamentals(F, m_max, rho, T, alpha)
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      REAL(dp), INTENT(IN)                               :: rho, T, alpha

      INTEGER                                            :: m
      REAL(dp)                                           :: scaled_T, weight

      scaled_T = -alpha*T/(rho + alpha)
      weight = 1.0_dp
      DO m = 1, m_max + 1
         F(m) = EXP(scaled_T)*weight
         weight = weight*alpha/(rho + alpha)
      END DO
   END SUBROUTINE compute_gaussian_fundamentals

! **************************************************************************************************
!> \brief Compose radial kernels before applying their common normalization.
!> \param F coefficients of orders zero through m_max, excluding normalization;
!>        MIX_LG includes its distinct component normalizations in F
!> \param normalization common multiplicative normalization; one for MIX_LG
!> \param m_max maximum order, nonnegative
!> \param potential_type operator selector from input_constants
!> \param rho reduced Gaussian exponent, positive
!> \param T reduced squared distance between Gaussian product centers
!> \param work reusable scratch space for composing radial kernels
!> \param omega range-separation parameter; Gaussian operators use alpha=omega**2
!> \param cutoff_radius cutoff parameter for truncated Coulomb and TShPSC
!> \param scale_coulomb Coulomb coefficient in mixed operators
!> \param scale_longrange long-range coefficient in mixed operators
!> \param scale_gaussian Gaussian coefficient in mixed operators
!> \note Spatial screening and table initialization belong to the caller.
! **************************************************************************************************
   SUBROUTINE compute_r12_fundamentals(F, normalization, m_max, potential_type, rho, T, work, &
                                       omega, cutoff_radius, scale_coulomb, scale_longrange, scale_gaussian)
      REAL(dp), INTENT(OUT)                              :: normalization
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      INTEGER, INTENT(IN)                                :: potential_type
      REAL(dp), INTENT(IN)                               :: rho, T
      REAL(dp), INTENT(INOUT)                            :: work(m_max + 1)
      REAL(dp), INTENT(IN)                               :: omega, cutoff_radius, scale_coulomb, &
                                                            scale_longrange, scale_gaussian

      REAL(dp)                                           :: gaussian_normalization, R

      normalization = 2.0_dp*pi/rho
      SELECT CASE (potential_type)
      CASE (do_potential_coulomb)
         CALL compute_coulomb_fundamentals(F, m_max, T)
      CASE (do_potential_truncated)
         CALL compute_truncated_fundamentals(F, m_max, rho, T, cutoff_radius)
      CASE (do_potential_TShPSC)
         R = cutoff_radius*SQRT(rho)
         F = 0.0_dp
         CALL trunc_CS_poly_n20(F, R, T, m_max)
      CASE (do_potential_short)
         CALL compute_short_range_fundamentals(F, m_max, rho, T, omega, work)
      CASE (do_potential_long)
         CALL compute_long_range_fundamentals(F, m_max, rho, T, omega)
      CASE (do_potential_gaussian)
         CALL compute_gaussian_fundamentals(F, m_max, rho, T, omega**2)
         normalization = (pi/(rho + omega**2))**1.5_dp
      CASE (do_potential_mix_cl)
         CALL compute_coulomb_fundamentals(F, m_max, T)
         CALL compute_long_range_fundamentals(work, m_max, rho, T, omega)
         F = F*scale_coulomb + work*scale_longrange
      CASE (do_potential_mix_cl_trunc)
         CALL compute_truncated_fundamentals(F, m_max, rho, T, cutoff_radius)
         CALL compute_coulomb_fundamentals(work, m_max, T)
         F = F*(scale_coulomb + scale_longrange) - work*scale_longrange
         CALL compute_long_range_fundamentals(work, m_max, rho, T, omega)
         F = F + work*scale_longrange
      CASE (do_potential_mix_lg)
         CALL compute_long_range_fundamentals(F, m_max, rho, T, omega)
         CALL compute_gaussian_fundamentals(work, m_max, rho, T, omega**2)
         gaussian_normalization = (pi/(rho + omega**2))**1.5_dp
         F = work*gaussian_normalization*scale_gaussian + F*normalization*scale_longrange
         normalization = 1.0_dp
      CASE DEFAULT
         CPABORT("Unsupported operator for radial r12 fundamentals")
      END SELECT
   END SUBROUTINE compute_r12_fundamentals

END MODULE r12_fundamentals
