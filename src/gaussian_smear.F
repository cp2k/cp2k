!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gaussian smearing: compute occupation numbers, chemical potential, and
!>        the entropic free-energy correction (-TS) using Gaussian broadening.
!>
!>        Occupation function:   f_i = (maxocc/2) * erfc[(e_i - mu) / sigma]
!>        Entropy correction:    kTS = -(sigma / 2*sqrt(pi)) * SUM_i maxocc * exp[-((e_i-mu)/sigma)^2]
!>
!>        The interface mirrors fermi_utils so that the two modules are
!>        interchangeable at every call site.
! **************************************************************************************************
MODULE gaussian_smear

   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC                               :: Gaussian, GaussianFixed, GaussianFixedDeriv, &
                                           GaussianFixedDerivMV, Gaussiankp, Gaussiankp2

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gaussian_smear'
   INTEGER, PARAMETER, PRIVATE          :: BISECT_MAX_ITER = 400
   REAL(KIND=dp), PARAMETER, PRIVATE    :: sqrtpi = 1.7724538509055159_dp  ! sqrt(pi)

CONTAINS

! **************************************************************************************************
!> \brief  Returns occupations according to Gaussian smearing
!>         for a given set of energies and chemical potential.
!>         f_i = (maxocc/2) * erfc( (e_i - mu) / sigma )
!> \param f       occupations (output)
!> \param N       total electron count (output)
!> \param kTS     entropic correction -T*S (output)
!> \param e       eigenvalues (input)
!> \param mu      chemical potential (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation per orbital (input)
!> \param estate  excited-state index for XAS (optional)
!> \param festate occupation of that excited state (optional)
! **************************************************************************************************
   SUBROUTINE Gaussian(f, N, kTS, e, mu, sigma, maxocc, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: f(:), N, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), mu, sigma, maxocc
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      INTEGER                                            :: i, Nstate
      REAL(KIND=dp)                                      :: occupation, x

      Nstate = SIZE(e)
      kTS = 0.0_dp

      DO i = 1, Nstate
         IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
            IF (i == estate) THEN
               occupation = festate
            ELSE
               occupation = maxocc
            END IF
         ELSE
            occupation = maxocc
         END IF

         x = (e(i) - mu)/sigma
         f(i) = occupation*0.5_dp*ERFC(x)
         ! Entropic correction: kTS = -sigma/(2*sqrt(pi)) * SUM occupation * exp(-x^2)
         kTS = kTS - (sigma/(2.0_dp*sqrtpi))*occupation*EXP(-x*x)
      END DO

      N = accurate_sum(f)

   END SUBROUTINE Gaussian

! **************************************************************************************************
!> \brief  Gaussian smearing for k-point cases (internal, 2D arrays)
!> \param f       occupations (nmo x nkp, output)
!> \param nel     total electron count (output)
!> \param kTS     entropic correction (output)
!> \param e       eigenvalues (nmo x nkp, input)
!> \param mu      chemical potential (input)
!> \param wk      k-point weights (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
! **************************************************************************************************
   SUBROUTINE Gaussian2(f, nel, kTS, e, mu, wk, sigma, maxocc)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: f
      REAL(KIND=dp), INTENT(OUT)                         :: nel, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: e
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma, maxocc

      INTEGER                                            :: ik, is, nkp, nmo
      REAL(KIND=dp)                                      :: x

      nmo = SIZE(e, 1)
      nkp = SIZE(e, 2)
      kTS = 0.0_dp

      IF (sigma > 1.0e-14_dp) THEN
         DO ik = 1, nkp
            DO is = 1, nmo
               x = (e(is, ik) - mu)/sigma
               f(is, ik) = maxocc*0.5_dp*ERFC(x)
               kTS = kTS - (sigma/(2.0_dp*sqrtpi))*maxocc*EXP(-x*x)*wk(ik)
            END DO
         END DO
      ELSE
         ! Zero-width limit: step function
         DO ik = 1, nkp
            DO is = 1, nmo
               IF (e(is, ik) <= mu) THEN
                  f(is, ik) = maxocc
               ELSE
                  f(is, ik) = 0.0_dp
               END IF
            END DO
         END DO
      END IF

      nel = 0.0_dp
      DO ik = 1, nkp
         nel = nel + accurate_sum(f(1:nmo, ik))*wk(ik)
      END DO

   END SUBROUTINE Gaussian2

! **************************************************************************************************
!> \brief  Bisection search for mu so that total electron count = N
!>         (Gamma-point version)
!> \param f       occupations (output)
!> \param mu      chemical potential found by bisection (output)
!> \param kTS     entropic correction (output)
!> \param e       eigenvalues (input)
!> \param N       target electron count (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param estate  XAS excited state (optional)
!> \param festate XAS occupation (optional)
! **************************************************************************************************
   SUBROUTINE GaussianFixed(f, mu, kTS, e, N, sigma, maxocc, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: f(:), mu, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), N, sigma, maxocc
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      INTEGER                                            :: iter, my_estate
      REAL(KIND=dp)                                      :: mu_max, mu_min, mu_now, my_festate, &
                                                            N_now, N_tmp

      IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
         my_estate = estate
         my_festate = festate
      ELSE
         my_estate = NINT(maxocc)
         my_festate = my_estate
      END IF

      ! ----- bracket mu -----
      ! erfc decays much faster than Fermi function, so 10*sigma is very generous
      mu_min = MINVAL(e)
      iter = 0
      DO
         iter = iter + 1
         CALL Gaussian(f, N_tmp, kTS, e, mu_min, sigma, maxocc, my_estate, my_festate)
         IF (N_tmp > N .OR. iter > 20) THEN
            mu_min = mu_min - sigma
         ELSE
            EXIT
         END IF
      END DO

      mu_max = MAXVAL(e)
      iter = 0
      DO
         iter = iter + 1
         CALL Gaussian(f, N_tmp, kTS, e, mu_max, sigma, maxocc, my_estate, my_festate)
         IF (N_tmp < N .OR. iter > 20) THEN
            mu_max = mu_max + sigma
         ELSE
            EXIT
         END IF
      END DO

      ! ----- bisection -----
      iter = 0
      DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
         iter = iter + 1
         mu_now = (mu_max + mu_min)/2.0_dp
         CALL Gaussian(f, N_now, kTS, e, mu_now, sigma, maxocc, my_estate, my_festate)

         IF (N_now <= N) THEN
            mu_min = mu_now
         ELSE
            mu_max = mu_now
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Gaussian smearing: maximum bisection iterations reached while finding mu")
            EXIT
         END IF
      END DO

      mu = (mu_max + mu_min)/2.0_dp
      CALL Gaussian(f, N_now, kTS, e, mu, sigma, maxocc, my_estate, my_festate)

   END SUBROUTINE GaussianFixed

! **************************************************************************************************
!> \brief  Bisection search for mu (k-point case, single spin channel or
!>         spin-degenerate)
!> \param f       occupations (nmo x nkp, output)
!> \param mu      chemical potential (output)
!> \param kTS     entropic correction (output)
!> \param e       eigenvalues (nmo x nkp, input)
!> \param nel     target electron count (input)
!> \param wk      k-point weights (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
! **************************************************************************************************
   SUBROUTINE Gaussiankp(f, mu, kTS, e, nel, wk, sigma, maxocc)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: f
      REAL(KIND=dp), INTENT(OUT)                         :: mu, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: e
      REAL(KIND=dp), INTENT(IN)                          :: nel
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma, maxocc

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER                                            :: iter
      REAL(KIND=dp)                                      :: de, mu_max, mu_min, N_now

      de = MAX(10.0_dp*sigma, 0.5_dp)
      mu_min = MINVAL(e) - de
      mu_max = MAXVAL(e) + de

      iter = 0
      DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
         iter = iter + 1
         mu = (mu_max + mu_min)/2.0_dp
         CALL Gaussian2(f, N_now, kTS, e, mu, wk, sigma, maxocc)

         IF (ABS(N_now - nel) < nel*epsocc) EXIT

         IF (N_now <= nel) THEN
            mu_min = mu
         ELSE
            mu_max = mu
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Gaussian smearing (kp): maximum bisection iterations reached")
            EXIT
         END IF
      END DO

      mu = (mu_max + mu_min)/2.0_dp
      CALL Gaussian2(f, N_now, kTS, e, mu, wk, sigma, maxocc)

   END SUBROUTINE Gaussiankp

! **************************************************************************************************
!> \brief  Bisection search for mu (k-point, spin-polarised, shared mu)
!> \param f       occupations (nmo x nkp x 2, output)
!> \param mu      chemical potential (output)
!> \param kTS     entropic correction (output)
!> \param e       eigenvalues (nmo x nkp x 2, input)
!> \param nel     target total electron count (input)
!> \param wk      k-point weights (input)
!> \param sigma   smearing width (input)
! **************************************************************************************************
   SUBROUTINE Gaussiankp2(f, mu, kTS, e, nel, wk, sigma)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: f
      REAL(KIND=dp), INTENT(OUT)                         :: mu, kTS
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: e
      REAL(KIND=dp), INTENT(IN)                          :: nel
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER                                            :: iter
      REAL(KIND=dp)                                      :: de, kTSa, kTSb, mu_max, mu_min, N_now, &
                                                            na, nb

      ! spin-polarised only
      CPASSERT(SIZE(f, 3) == 2 .AND. SIZE(e, 3) == 2)

      de = MAX(10.0_dp*sigma, 0.5_dp)
      mu_min = MINVAL(e) - de
      mu_max = MAXVAL(e) + de

      iter = 0
      DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
         iter = iter + 1
         mu = (mu_max + mu_min)/2.0_dp
         CALL Gaussian2(f(:, :, 1), na, kTSa, e(:, :, 1), mu, wk, sigma, 1.0_dp)
         CALL Gaussian2(f(:, :, 2), nb, kTSb, e(:, :, 2), mu, wk, sigma, 1.0_dp)
         N_now = na + nb

         IF (ABS(N_now - nel) < nel*epsocc) EXIT

         IF (N_now <= nel) THEN
            mu_min = mu
         ELSE
            mu_max = mu
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Gaussian smearing (kp2): maximum bisection iterations reached")
            EXIT
         END IF
      END DO

      mu = (mu_max + mu_min)/2.0_dp
      CALL Gaussian2(f(:, :, 1), na, kTSa, e(:, :, 1), mu, wk, sigma, 1.0_dp)
      CALL Gaussian2(f(:, :, 2), nb, kTSb, e(:, :, 2), mu, wk, sigma, 1.0_dp)
      kTS = kTSa + kTSb

   END SUBROUTINE Gaussiankp2

! **************************************************************************************************
!> \brief  Analytical Jacobian df_i/de_j for Gaussian smearing under the
!>         electron-number constraint  sum(f) = N.
!>
!>         Derivation
!>         ----------
!>         Direct partial derivatives (mu held fixed):
!>            partial f_i / partial e_i = -g_i        (diagonal)
!>            partial f_i / partial mu  = +g_i
!>         where g_i = occupation_i / (sigma * sqrt(pi)) * exp[-((e_i - mu)/sigma)^2]
!>
!>         Differentiating the constraint  sum_i f_i = N  w.r.t. e_j:
!>            -g_j + G * dmu/de_j = 0    =>    dmu/de_j = g_j / G
!>         where G = sum_i g_i.
!>
!>         Full constrained derivative:
!>            df_i/de_j = -delta_{ij} * g_i  +  g_i * g_j / G
!>
!>         This is a diagonal matrix plus a symmetric rank-1 update,
!>         computed in O(N) time (vs. O(N^2 * bisection) for the numerical version).
!>
!> \param dfde    Jacobian matrix dfde(i,j) = df_i/de_j  (Nstate x Nstate, output)
!> \param f       occupations (output)
!> \param mu      chemical potential (output)
!> \param kTS     entropic correction (output)
!> \param e       eigenvalues (input)
!> \param N       target electron count (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param l       unused, kept for interface compatibility with FermiFixedDeriv
!> \param estate  XAS excited state (optional)
!> \param festate XAS occupation (optional)
! **************************************************************************************************
   SUBROUTINE GaussianFixedDeriv(dfde, f, mu, kTS, e, N, sigma, maxocc, l, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: dfde(:, :), f(:), mu, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), N, sigma, maxocc, l
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      CHARACTER(len=*), PARAMETER :: routineN = 'GaussianFixedDeriv'

      INTEGER                                            :: handle, i, j, my_estate, Nstate
      REAL(KIND=dp)                                      :: Gsum, my_festate, occ_i, prefactor, x
      REAL(KIND=dp), ALLOCATABLE                         :: gvec(:)

      MARK_USED(l)

      CALL timeset(routineN, handle)

      IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
         my_estate = estate
         my_festate = festate
      ELSE
         my_estate = NINT(maxocc)
         my_festate = my_estate
      END IF

      ! Step 1: find mu and f by bisection (one call only)
      CALL GaussianFixed(f, mu, kTS, e, N, sigma, maxocc, my_estate, my_festate)

      ! Step 2: build the Gaussian weight vector g_i
      Nstate = SIZE(e)
      ALLOCATE (gvec(Nstate))

      prefactor = 1.0_dp/(sigma*sqrtpi)
      DO i = 1, Nstate
         IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
            IF (i == estate) THEN
               occ_i = festate
            ELSE
               occ_i = maxocc
            END IF
         ELSE
            occ_i = maxocc
         END IF
         x = (e(i) - mu)/sigma
         gvec(i) = occ_i*prefactor*EXP(-x*x)
      END DO

      Gsum = accurate_sum(gvec)

      ! Step 3: assemble dfde(i,j) = -delta_{ij}*g_i + g_i*g_j/G
      ! Guard against Gsum = 0 (zero-temperature limit, all states fully occupied or empty)
      IF (Gsum > EPSILON(Gsum)) THEN
         DO j = 1, Nstate
            DO i = 1, Nstate
               dfde(i, j) = gvec(i)*gvec(j)/Gsum
            END DO
            ! diagonal correction
            dfde(j, j) = dfde(j, j) - gvec(j)
         END DO
      ELSE
         dfde(:, :) = 0.0_dp
      END IF

      DEALLOCATE (gvec)

      CALL timestop(handle)

   END SUBROUTINE GaussianFixedDeriv

! **************************************************************************************************
!> \brief  Apply the Jacobian TRANSPOSE(df/de) to a vector, WITHOUT forming
!>         the full N x N matrix.   O(N) time and O(N) memory.
!>
!>         Exploiting the rank-1 structure:
!>            [TRANSPOSE(dfde) * v]_j = -g_j * v_j  +  g_j * (g . v) / G
!>
!>         This routine replaces the pattern
!>            ALLOCATE(dfde(nmo,nmo))
!>            CALL GaussianFixedDeriv(dfde, ...)
!>            result = MATMUL(TRANSPOSE(dfde), v)
!>            DEALLOCATE(dfde)
!>         used in qs_mo_occupation / set_mo_occupation_1, turning O(N^2)
!>         storage + O(N^2) MATMUL into O(N) throughout.
!>
!> \param RESULT   output vector  = TRANSPOSE(df/de) * v
!> \param f        occupations (output)
!> \param mu       chemical potential (output)
!> \param kTS      entropic correction (output)
!> \param e        eigenvalues (input)
!> \param N_el     target electron count (input)
!> \param sigma    smearing width (input)
!> \param maxocc   maximum occupation (input)
!> \param v        input vector to multiply (Nstate, input)
!> \param estate   XAS (optional)
!> \param festate  XAS (optional)
! **************************************************************************************************
   SUBROUTINE GaussianFixedDerivMV(RESULT, f, mu, kTS, e, N_el, sigma, maxocc, v, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: result(:), f(:), mu, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), N_el, sigma, maxocc, v(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      CHARACTER(len=*), PARAMETER :: routineN = 'GaussianFixedDerivMV'

      INTEGER                                            :: handle, i, my_estate, Nstate
      REAL(KIND=dp)                                      :: gdotv, Gsum, my_festate, occ_i, &
                                                            prefactor, x
      REAL(KIND=dp), ALLOCATABLE                         :: gvec(:)

      CALL timeset(routineN, handle)

      IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
         my_estate = estate
         my_festate = festate
      ELSE
         my_estate = NINT(maxocc)
         my_festate = my_estate
      END IF

      ! Step 1: find mu and f
      CALL GaussianFixed(f, mu, kTS, e, N_el, sigma, maxocc, my_estate, my_festate)

      ! Step 2: Gaussian weight vector
      Nstate = SIZE(e)
      ALLOCATE (gvec(Nstate))

      prefactor = 1.0_dp/(sigma*sqrtpi)
      DO i = 1, Nstate
         IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
            IF (i == estate) THEN
               occ_i = festate
            ELSE
               occ_i = maxocc
            END IF
         ELSE
            occ_i = maxocc
         END IF
         x = (e(i) - mu)/sigma
         gvec(i) = occ_i*prefactor*EXP(-x*x)
      END DO

      Gsum = accurate_sum(gvec)

      ! Step 3: result_j = -g_j * v_j + g_j * (g . v) / G
      IF (Gsum > EPSILON(Gsum)) THEN
         gdotv = 0.0_dp
         DO i = 1, Nstate
            gdotv = gdotv + gvec(i)*v(i)
         END DO
         DO i = 1, Nstate
            RESULT(i) = gvec(i)*(gdotv/Gsum - v(i))
         END DO
      ELSE
         RESULT(:) = 0.0_dp
      END IF

      DEALLOCATE (gvec)

      CALL timestop(handle)

   END SUBROUTINE GaussianFixedDerivMV

END MODULE gaussian_smear
