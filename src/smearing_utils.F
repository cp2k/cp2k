!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unified smearing module supporting four methods:
!>          smear_fermi_dirac  — Fermi-Dirac distribution
!>          smear_gaussian     — Gaussian broadening
!>          smear_mp           — Methfessel-Paxton first order
!>          smear_mv           — Marzari-Vanderbilt (cold smearing)
!>
!>        All methods share the bisection framework, LFOMO/HOMO logic, and the
!>        analytical rank-1 Jacobian.  Only the per-state math (f, kTS, g_i)
!>        differs, selected via a method integer from input_constants.
!>
!> \par History
!>      09.2008: Created (fermi_utils.F)
!>      02.2026: Extended to more smearing method and renamed
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE smearing_utils

   USE bibliography,                    ONLY: FuHo1983,&
                                              Marzari1999,&
                                              Mermin1965,&
                                              MethfesselPaxton1989,&
                                              cite_reference,&
                                              dosSantos2023
   USE input_constants,                 ONLY: smear_fermi_dirac,&
                                              smear_gaussian,&
                                              smear_mp,&
                                              smear_mv
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: rootpi,&
                                              sqrt2,&
                                              sqrthalf
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Unified interface (method as parameter)
   PUBLIC :: SmearOcc, SmearFixed, SmearFixedDeriv, SmearFixedDerivMV
   PUBLIC :: Smearkp, Smearkp2
   PRIVATE :: cite_smearing

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smearing_utils'
   INTEGER, PARAMETER, PRIVATE                           :: BISECT_MAX_ITER = 400
   INTEGER, PARAMETER, PRIVATE                           :: NEWTON_MAX_ITER = 50

CONTAINS
! **************************************************************************************************
!> \brief   Citation of Smearing methods
!> \param method ...
! **************************************************************************************************
   SUBROUTINE cite_smearing(method)
      INTEGER, INTENT(IN)                                :: method

      SELECT CASE (method)
      CASE (smear_fermi_dirac)
         CALL cite_reference(Mermin1965)
      CASE (smear_gaussian)
         CALL cite_reference(FuHo1983)
      CASE (smear_mp)
         CALL cite_reference(FuHo1983)
         CALL cite_reference(MethfesselPaxton1989)
         CALL cite_reference(dosSantos2023)
      CASE (smear_mv)
         CALL cite_reference(FuHo1983)
         CALL cite_reference(Marzari1999)
         CALL cite_reference(dosSantos2023)
      END SELECT
   END SUBROUTINE cite_smearing

! **************************************************************************************************
!> \brief   Returns occupations and smearing correction for a given set of
!>          energies and chemical potential, using one of four smearing methods.
!>
!>          Fermi-Dirac:  f_i = occ / [1 + exp((e_i - mu)/sigma)]
!>          Gaussian:     f_i = (occ/2) * erfc[(e_i - mu)/sigma]
!>          MP-1:         f_i = (occ/2) * erfc(x) - occ*x/(2*sqrt(pi)) * exp(-x^2)
!>          MV:           f_i = (occ/2) * erfc(u) + occ/(sqrt(2*pi)) * exp(-u^2),  u = x + 1/sqrt(2)
!>
!>          kTS is the smearing correction to the free energy (physically -TS
!>          for Fermi-Dirac; a variational correction term for the other methods).
!>          It enters the total energy and the Gillan extrapolation E(0) = E - kTS/2.
!>
!> \param f       occupations (output)
!> \param N       total number of electrons (output)
!> \param kTS     smearing correction to the free energy (output)
!> \param e       eigenvalues (input)
!> \param mu      chemical potential (input)
!> \param sigma   smearing width: kT for Fermi-Dirac, sigma for others (input)
!> \param maxocc  maximum occupation of an orbital (input)
!> \param method  smearing method selector from input_constants (input)
!> \param estate  excited state index for core-level spectroscopy (optional)
!> \param festate occupation of the excited state (optional)
! **************************************************************************************************
   SUBROUTINE SmearOcc(f, N, kTS, e, mu, sigma, maxocc, method, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: f(:), N, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), mu, sigma, maxocc
      INTEGER, INTENT(IN)                                :: method
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      INTEGER                                            :: i, Nstate
      REAL(KIND=dp)                                      :: arg, expu2, expx2, occupation, term1, &
                                                            term2, tmp, tmp2, tmp3, tmp4, tmplog, &
                                                            u, x

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

         SELECT CASE (method)
         CASE (smear_fermi_dirac)
            IF (e(i) > mu) THEN
               arg = -(e(i) - mu)/sigma
               tmp = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = tmp/tmp4
               tmp3 = 1.0_dp/tmp4
               tmplog = -LOG(tmp4)
               term1 = tmp2*(arg + tmplog)
               term2 = tmp3*tmplog
            ELSE
               arg = (e(i) - mu)/sigma
               tmp = EXP(arg)
               tmp4 = tmp + 1.0_dp
               tmp2 = 1.0_dp/tmp4
               tmp3 = tmp/tmp4
               tmplog = -LOG(tmp4)
               term1 = tmp2*tmplog
               term2 = tmp3*(arg + tmplog)
            END IF
            f(i) = occupation*tmp2
            kTS = kTS + sigma*occupation*(term1 + term2)

         CASE (smear_gaussian)
            x = (e(i) - mu)/sigma
            expx2 = EXP(-x*x)
            f(i) = occupation*0.5_dp*ERFC(x)
            kTS = kTS - (sigma/(2.0_dp*rootpi))*occupation*expx2

         CASE (smear_mp)
            x = (e(i) - mu)/sigma
            expx2 = EXP(-x*x)
            f(i) = occupation*(0.5_dp*ERFC(x) - x/(2.0_dp*rootpi)*expx2)
            kTS = kTS + (sigma/(4.0_dp*rootpi))*occupation*(2.0_dp*x*x - 1.0_dp)*expx2

         CASE (smear_mv)
            x = (e(i) - mu)/sigma
            u = x + sqrthalf
            expu2 = EXP(-u*u)
            f(i) = occupation*(0.5_dp*ERFC(u) + expu2/(sqrt2*rootpi))
            kTS = kTS - (sigma/(sqrt2*rootpi))*occupation*u*expu2

         CASE DEFAULT
            CPABORT("SmearOcc: unknown smearing method")
         END SELECT
      END DO

      N = accurate_sum(f)

   END SUBROUTINE SmearOcc

! **************************************************************************************************
!> \brief   k-point version of SmearOcc (module-private).
!>          Computes occupations and kTS for a 2D array of eigenvalues
!>          (nmo x nkp) weighted by k-point weights.
!>          Falls back to a step function when sigma < 1e-14.
!>
!> \param f       occupations (nmo x nkp, output)
!> \param nel     total number of electrons (output)
!> \param kTS     smearing correction (output)
!> \param e       eigenvalues (nmo x nkp, input)
!> \param mu      chemical potential (input)
!> \param wk      k-point weights (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param method  smearing method selector (input)
! **************************************************************************************************
   SUBROUTINE Smear2(f, nel, kTS, e, mu, wk, sigma, maxocc, method)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: f
      REAL(KIND=dp), INTENT(OUT)                         :: nel, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: e
      REAL(KIND=dp), INTENT(IN)                          :: mu
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma, maxocc
      INTEGER, INTENT(IN)                                :: method

      INTEGER                                            :: ik, is, nkp, nmo
      REAL(KIND=dp)                                      :: arg, expu2, expx2, term1, term2, tmp, &
                                                            tmp2, tmp3, tmp4, tmplog, u, x

      nmo = SIZE(e, 1)
      nkp = SIZE(e, 2)
      kTS = 0.0_dp

      IF (sigma > 1.0e-14_dp) THEN
         DO ik = 1, nkp
            DO is = 1, nmo
               SELECT CASE (method)
               CASE (smear_fermi_dirac)
                  IF (e(is, ik) > mu) THEN
                     arg = -(e(is, ik) - mu)/sigma
                     tmp = EXP(arg)
                     tmp4 = tmp + 1.0_dp
                     tmp2 = tmp/tmp4
                     tmp3 = 1.0_dp/tmp4
                     tmplog = -LOG(tmp4)
                     term1 = tmp2*(arg + tmplog)
                     term2 = tmp3*tmplog
                  ELSE
                     arg = (e(is, ik) - mu)/sigma
                     tmp = EXP(arg)
                     tmp4 = tmp + 1.0_dp
                     tmp2 = 1.0_dp/tmp4
                     tmp3 = tmp/tmp4
                     tmplog = -LOG(tmp4)
                     term1 = tmp2*tmplog
                     term2 = tmp3*(arg + tmplog)
                  END IF
                  f(is, ik) = maxocc*tmp2
                  kTS = kTS + sigma*maxocc*(term1 + term2)*wk(ik)

               CASE (smear_gaussian)
                  x = (e(is, ik) - mu)/sigma
                  expx2 = EXP(-x*x)
                  f(is, ik) = maxocc*0.5_dp*ERFC(x)
                  kTS = kTS - (sigma/(2.0_dp*rootpi))*maxocc*expx2*wk(ik)

               CASE (smear_mp)
                  x = (e(is, ik) - mu)/sigma
                  expx2 = EXP(-x*x)
                  f(is, ik) = maxocc*(0.5_dp*ERFC(x) - x/(2.0_dp*rootpi)*expx2)
                  kTS = kTS + (sigma/(4.0_dp*rootpi))*maxocc*(2.0_dp*x*x - 1.0_dp)*expx2*wk(ik)

               CASE (smear_mv)
                  x = (e(is, ik) - mu)/sigma
                  u = x + sqrthalf
                  expu2 = EXP(-u*u)
                  f(is, ik) = maxocc*(0.5_dp*ERFC(u) + expu2/(sqrt2*rootpi))
                  kTS = kTS - (sigma/(sqrt2*rootpi))*maxocc*u*expu2*wk(ik)

               CASE DEFAULT
                  CPABORT("Smear2: unknown smearing method")
               END SELECT
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

   END SUBROUTINE Smear2

! **************************************************************************************************
!> \brief   Bisection search for the chemical potential mu such that the total
!>          electron count equals N, for a given smearing method (Gamma point).
!>          Brackets mu by expanding outward from [min(e), max(e)] in steps
!>          of sigma, then bisects to machine precision.
!>
!>          For MP-1 and MV: the occupation function is non-monotonic, so it's
!>          possible that pure bisection find a spurious root.
!>          We first bisect with Gaussian smearing to get a reliable initial mu,
!>          then refine with Newton's method using the actual method's dN/dmu.
!>          (dos Santos & Marzari, PRB 2023)
!>
!> \param f       occupations (output)
!> \param mu      chemical potential found by bisection (output)
!> \param kTS     smearing correction (output)
!> \param e       eigenvalues (input)
!> \param N       target number of electrons (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param method  smearing method selector (input)
!> \param estate  excited state index for core-level spectroscopy (optional)
!> \param festate occupation of the excited state (optional)
! **************************************************************************************************
   SUBROUTINE SmearFixed(f, mu, kTS, e, N, sigma, maxocc, method, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: f(:), mu, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), N, sigma, maxocc
      INTEGER, INTENT(IN)                                :: method
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      INTEGER                                            :: iter, my_estate, Nstate
      REAL(KIND=dp)                                      :: Gsum, mu_max, mu_min, mu_now, &
                                                            my_festate, N_now, N_tmp
      REAL(KIND=dp), ALLOCATABLE                         :: gvec(:)

      IF (PRESENT(estate) .AND. PRESENT(festate)) THEN
         my_estate = estate
         my_festate = festate
      ELSE
         my_estate = NINT(maxocc)
         my_festate = my_estate
      END IF

      Nstate = SIZE(e)

      CALL cite_smearing(method)

      SELECT CASE (method)

         ! Non-monotonic methods: Gaussian bisection + Newton refinement
      CASE (smear_mp, smear_mv)
         ! Step 1: Gaussian bisection for a reliable initial mu
         mu_min = MINVAL(e)
         iter = 0
         DO
            iter = iter + 1
            CALL SmearOcc(f, N_tmp, kTS, e, mu_min, sigma, maxocc, smear_gaussian, my_estate, my_festate)
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
            CALL SmearOcc(f, N_tmp, kTS, e, mu_max, sigma, maxocc, smear_gaussian, my_estate, my_festate)
            IF (N_tmp < N .OR. iter > 20) THEN
               mu_max = mu_max + sigma
            ELSE
               EXIT
            END IF
         END DO

         iter = 0
         DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
            iter = iter + 1
            mu_now = (mu_max + mu_min)/2.0_dp
            CALL SmearOcc(f, N_now, kTS, e, mu_now, sigma, maxocc, smear_gaussian, my_estate, my_festate)
            IF (N_now <= N) THEN
               mu_min = mu_now
            ELSE
               mu_max = mu_now
            END IF
            IF (iter > BISECT_MAX_ITER) EXIT
         END DO
         mu = (mu_max + mu_min)/2.0_dp

         ! Step 2: Newton refinement with the actual method
         ALLOCATE (gvec(Nstate))
         DO iter = 1, NEWTON_MAX_ITER
            CALL SmearOcc(f, N_now, kTS, e, mu, sigma, maxocc, method, my_estate, my_festate)
            IF (ABS(N_now - N) < N*1.0e-12_dp) EXIT
            CALL compute_gvec(gvec, f, e, mu, sigma, maxocc, Nstate, method, my_estate, my_festate)
            Gsum = accurate_sum(gvec)
            IF (ABS(Gsum) < EPSILON(Gsum)) EXIT
            mu = mu + (N - N_now)/Gsum
         END DO
         DEALLOCATE (gvec)

         ! Final evaluation
         CALL SmearOcc(f, N_now, kTS, e, mu, sigma, maxocc, method, my_estate, my_festate)

         ! Monotonic methods (FD, Gaussian): pure bisection
      CASE DEFAULT
         mu_min = MINVAL(e)
         iter = 0
         DO
            iter = iter + 1
            CALL SmearOcc(f, N_tmp, kTS, e, mu_min, sigma, maxocc, method, my_estate, my_festate)
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
            CALL SmearOcc(f, N_tmp, kTS, e, mu_max, sigma, maxocc, method, my_estate, my_festate)
            IF (N_tmp < N .OR. iter > 20) THEN
               mu_max = mu_max + sigma
            ELSE
               EXIT
            END IF
         END DO

         iter = 0
         DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
            iter = iter + 1
            mu_now = (mu_max + mu_min)/2.0_dp
            CALL SmearOcc(f, N_now, kTS, e, mu_now, sigma, maxocc, method, my_estate, my_festate)
            IF (N_now <= N) THEN
               mu_min = mu_now
            ELSE
               mu_max = mu_now
            END IF
            IF (iter > BISECT_MAX_ITER) THEN
               CPWARN("SmearFixed: maximum bisection iterations reached")
               EXIT
            END IF
         END DO

         mu = (mu_max + mu_min)/2.0_dp
         CALL SmearOcc(f, N_now, kTS, e, mu, sigma, maxocc, method, my_estate, my_festate)

      END SELECT

   END SUBROUTINE SmearFixed

! **************************************************************************************************
!> \brief   Bisection search for mu given a target electron count (k-point case,
!>          single spin channel or spin-degenerate).
!>          Initial bracket width is max(10*sigma, 0.5) for Gaussian/MP/MV,
!>          or sigma*ln[(1-eps)/eps] for Fermi-Dirac, reflecting the different
!>          tail decay rates.
!>
!>          For MP-1 and MV: Gaussian bisection + Newton refinement
!>          (dos Santos & Marzari, PRB 2023).
!>
!> \param f       occupations (nmo x nkp, output)
!> \param mu      chemical potential (output)
!> \param kTS     smearing correction (output)
!> \param e       eigenvalues (nmo x nkp, input)
!> \param nel     target number of electrons (input)
!> \param wk      k-point weights (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param method  smearing method selector (input)
! **************************************************************************************************
   SUBROUTINE Smearkp(f, mu, kTS, e, nel, wk, sigma, maxocc, method)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: f
      REAL(KIND=dp), INTENT(OUT)                         :: mu, kTS
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: e
      REAL(KIND=dp), INTENT(IN)                          :: nel
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma, maxocc
      INTEGER, INTENT(IN)                                :: method

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER                                            :: bisect_method, ik, is, iter, nkp, nmo
      REAL(KIND=dp)                                      :: de, dNdmu, expu2, expx2, mu_max, mu_min, &
                                                            N_now, u, x

      nmo = SIZE(e, 1)
      nkp = SIZE(e, 2)

      CALL cite_smearing(method)

      ! Choose bisection method: Gaussian for MP/MV, actual method for FD/Gaussian
      SELECT CASE (method)
      CASE (smear_mp, smear_mv)
         bisect_method = smear_gaussian
      CASE DEFAULT
         bisect_method = method
      END SELECT

      ! Initial bracket
      SELECT CASE (bisect_method)
      CASE (smear_fermi_dirac)
         de = sigma*LOG((1.0_dp - epsocc)/epsocc)
      CASE DEFAULT
         de = 10.0_dp*sigma
      END SELECT
      de = MAX(de, 0.5_dp)

      ! Bisection with bisect_method
      mu_min = MINVAL(e) - de
      mu_max = MAXVAL(e) + de
      iter = 0
      DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
         iter = iter + 1
         mu = (mu_max + mu_min)/2.0_dp
         CALL Smear2(f, N_now, kTS, e, mu, wk, sigma, maxocc, bisect_method)

         IF (ABS(N_now - nel) < nel*epsocc) EXIT

         IF (N_now <= nel) THEN
            mu_min = mu
         ELSE
            mu_max = mu
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Smearkp: maximum bisection iterations reached")
            EXIT
         END IF
      END DO
      mu = (mu_max + mu_min)/2.0_dp

      ! Newton refinement for non-monotonic methods
      SELECT CASE (method)
      CASE (smear_mp, smear_mv)
         DO iter = 1, NEWTON_MAX_ITER
            CALL Smear2(f, N_now, kTS, e, mu, wk, sigma, maxocc, method)
            IF (ABS(N_now - nel) < nel*epsocc) EXIT

            ! Compute dN/dmu = sum_{ik} wk * g_i(k)  inline
            dNdmu = 0.0_dp
            DO ik = 1, nkp
               DO is = 1, nmo
                  x = (e(is, ik) - mu)/sigma
                  SELECT CASE (method)
                  CASE (smear_mp)
                     expx2 = EXP(-x*x)
                     dNdmu = dNdmu + maxocc*(3.0_dp - 2.0_dp*x*x)/(2.0_dp*sigma*rootpi)*expx2*wk(ik)
                  CASE (smear_mv)
                     u = x + sqrthalf
                     expu2 = EXP(-u*u)
                     dNdmu = dNdmu + maxocc*(2.0_dp + sqrt2*x)/(sigma*rootpi)*expu2*wk(ik)
                  END SELECT
               END DO
            END DO

            IF (ABS(dNdmu) < EPSILON(dNdmu)) EXIT
            mu = mu + (nel - N_now)/dNdmu
         END DO
      END SELECT

      ! Final evaluation with the actual method
      CALL Smear2(f, N_now, kTS, e, mu, wk, sigma, maxocc, method)

   END SUBROUTINE Smearkp

! **************************************************************************************************
!> \brief   Bisection search for mu (k-point, spin-polarised with a shared
!>          chemical potential across both spin channels).
!>          Asserts that the third dimension of f and e is exactly 2.
!>
!>          For MP-1 and MV: Gaussian bisection + Newton refinement.
!>
!> \param f       occupations (nmo x nkp x 2, output)
!> \param mu      chemical potential (output)
!> \param kTS     smearing correction (output)
!> \param e       eigenvalues (nmo x nkp x 2, input)
!> \param nel     target total number of electrons (input)
!> \param wk      k-point weights (input)
!> \param sigma   smearing width (input)
!> \param method  smearing method selector (input)
! **************************************************************************************************
   SUBROUTINE Smearkp2(f, mu, kTS, e, nel, wk, sigma, method)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: f
      REAL(KIND=dp), INTENT(OUT)                         :: mu, kTS
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: e
      REAL(KIND=dp), INTENT(IN)                          :: nel
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: wk
      REAL(KIND=dp), INTENT(IN)                          :: sigma
      INTEGER, INTENT(IN)                                :: method

      REAL(KIND=dp), PARAMETER                           :: epsocc = 1.0e-12_dp

      INTEGER                                            :: bisect_method, ik, is, ispin, iter, nkp, &
                                                            nmo
      REAL(KIND=dp)                                      :: de, dNdmu, expu2, expx2, kTSa, kTSb, &
                                                            mu_max, mu_min, N_now, na, nb, u, x

      CPASSERT(SIZE(f, 3) == 2 .AND. SIZE(e, 3) == 2)

      nmo = SIZE(e, 1)
      nkp = SIZE(e, 2)

      CALL cite_smearing(method)

      SELECT CASE (method)
      CASE (smear_mp, smear_mv)
         bisect_method = smear_gaussian
      CASE DEFAULT
         bisect_method = method
      END SELECT

      SELECT CASE (bisect_method)
      CASE (smear_fermi_dirac)
         de = sigma*LOG((1.0_dp - epsocc)/epsocc)
      CASE DEFAULT
         de = 10.0_dp*sigma
      END SELECT
      de = MAX(de, 0.5_dp)

      ! Bisection with bisect_method
      mu_min = MINVAL(e) - de
      mu_max = MAXVAL(e) + de
      iter = 0
      DO WHILE (mu_max - mu_min > EPSILON(mu)*MAX(1.0_dp, ABS(mu_max), ABS(mu_min)))
         iter = iter + 1
         mu = (mu_max + mu_min)/2.0_dp
         CALL Smear2(f(:, :, 1), na, kTSa, e(:, :, 1), mu, wk, sigma, 1.0_dp, bisect_method)
         CALL Smear2(f(:, :, 2), nb, kTSb, e(:, :, 2), mu, wk, sigma, 1.0_dp, bisect_method)
         N_now = na + nb

         IF (ABS(N_now - nel) < nel*epsocc) EXIT

         IF (N_now <= nel) THEN
            mu_min = mu
         ELSE
            mu_max = mu
         END IF

         IF (iter > BISECT_MAX_ITER) THEN
            CPWARN("Smearkp2: maximum bisection iterations reached")
            EXIT
         END IF
      END DO
      mu = (mu_max + mu_min)/2.0_dp

      ! Newton refinement for non-monotonic methods
      SELECT CASE (method)
      CASE (smear_mp, smear_mv)
         DO iter = 1, NEWTON_MAX_ITER
            CALL Smear2(f(:, :, 1), na, kTSa, e(:, :, 1), mu, wk, sigma, 1.0_dp, method)
            CALL Smear2(f(:, :, 2), nb, kTSb, e(:, :, 2), mu, wk, sigma, 1.0_dp, method)
            N_now = na + nb
            IF (ABS(N_now - nel) < nel*epsocc) EXIT

            ! dN/dmu across both spin channels (maxocc=1 per spin)
            dNdmu = 0.0_dp
            DO ispin = 1, 2
               DO ik = 1, nkp
                  DO is = 1, nmo
                     x = (e(is, ik, ispin) - mu)/sigma
                     SELECT CASE (method)
                     CASE (smear_mp)
                        expx2 = EXP(-x*x)
                        dNdmu = dNdmu + (3.0_dp - 2.0_dp*x*x)/(2.0_dp*sigma*rootpi)*expx2*wk(ik)
                     CASE (smear_mv)
                        u = x + sqrthalf
                        expu2 = EXP(-u*u)
                        dNdmu = dNdmu + (2.0_dp + sqrt2*x)/(sigma*rootpi)*expu2*wk(ik)
                     END SELECT
                  END DO
               END DO
            END DO

            IF (ABS(dNdmu) < EPSILON(dNdmu)) EXIT
            mu = mu + (nel - N_now)/dNdmu
         END DO
      END SELECT

      ! Final evaluation with the actual method
      CALL Smear2(f(:, :, 1), na, kTSa, e(:, :, 1), mu, wk, sigma, 1.0_dp, method)
      CALL Smear2(f(:, :, 2), nb, kTSb, e(:, :, 2), mu, wk, sigma, 1.0_dp, method)
      kTS = kTSa + kTSb

   END SUBROUTINE Smearkp2

! **************************************************************************************************
!> \brief   Computes the smearing weight vector g_i = -df_i/de_i with mu held
!>          fixed (module-private helper for the analytical Jacobian routines).
!>
!>          Fermi-Dirac:  g_i = occ * f_norm * (1 - f_norm) / sigma
!>                        where f_norm = f_i/occ_i  (overflow-safe, uses
!>                        pre-computed f rather than re-evaluating exp)
!>          Gaussian:     g_i = occ / (sigma*sqrt(pi)) * exp(-x^2)
!>          MP-1:         g_i = occ * (3 - 2*x^2) / (2*sigma*sqrt(pi)) * exp(-x^2)
!>          MV:           g_i = occ * (2 + sqrt(2)*x) / (sigma*sqrt(pi)) * exp(-u^2)
!>
!>          Note: g_i can be negative for MP-1 (|x| > sqrt(3/2)) and MV
!>          (x < -sqrt(2)).  This does not break bisection (N(mu) is still
!>          monotone) but the Jacobian routines must guard against G = sum(g_i)
!>          being near zero.
!>
!> \param gvec    weight vector (Nstate, output)
!> \param f       occupations from a prior SmearOcc/SmearFixed call (input)
!> \param e       eigenvalues (input)
!> \param mu      chemical potential (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param Nstate  number of states (input)
!> \param method  smearing method selector (input)
!> \param estate  excited state index (optional)
!> \param festate occupation of the excited state (optional)
! **************************************************************************************************
   SUBROUTINE compute_gvec(gvec, f, e, mu, sigma, maxocc, Nstate, method, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: gvec(:)
      REAL(KIND=dp), INTENT(IN)                          :: f(:), e(:), mu, sigma, maxocc
      INTEGER, INTENT(IN)                                :: Nstate, method
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: expu2, expx2, fi_norm, occ_i, u, x

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

         IF (occ_i < EPSILON(occ_i)) THEN
            gvec(i) = 0.0_dp
            CYCLE
         END IF

         x = (e(i) - mu)/sigma

         SELECT CASE (method)
         CASE (smear_fermi_dirac)
            fi_norm = f(i)/occ_i
            gvec(i) = occ_i*fi_norm*(1.0_dp - fi_norm)/sigma

         CASE (smear_gaussian)
            expx2 = EXP(-x*x)
            gvec(i) = occ_i/(sigma*rootpi)*expx2

         CASE (smear_mp)
            expx2 = EXP(-x*x)
            gvec(i) = occ_i*(3.0_dp - 2.0_dp*x*x)/(2.0_dp*sigma*rootpi)*expx2

         CASE (smear_mv)
            u = x + sqrthalf
            expu2 = EXP(-u*u)
            gvec(i) = occ_i*(2.0_dp + sqrt2*x)/(sigma*rootpi)*expu2

         END SELECT
      END DO

   END SUBROUTINE compute_gvec

! **************************************************************************************************
!> \brief   Analytical Jacobian df_i/de_j for any smearing method under the
!>          electron-number constraint sum(f) = N.
!>
!>          Differentiating f_i(e, mu(e)) where mu is implicitly defined by
!>          the constraint yields:
!>
!>            df_i/de_j  =  -delta_{ij} * g_i  +  g_i * g_j / G
!>
!>          where g_i = -df_i/de_i (mu fixed) and G = sum(g_i).
!>          This is a diagonal matrix plus a symmetric rank-1 update.
!>          Building it costs O(N) for g, plus O(N^2) for the outer product.
!>
!>          Replaces the original numerical finite-difference FermiFixedDeriv
!>          which required 2N bisection solves.  Exact to machine precision
!>          for all four methods.
!>
!> \param dfde    Jacobian matrix dfde(i,j) = df_i/de_j (Nstate x Nstate, output)
!> \param f       occupations (output)
!> \param mu      chemical potential (output)
!> \param kTS     smearing correction (output)
!> \param e       eigenvalues (input)
!> \param N       target number of electrons (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param method  smearing method selector (input)
!> \param estate  excited state index (optional)
!> \param festate occupation of the excited state (optional)
! **************************************************************************************************
   SUBROUTINE SmearFixedDeriv(dfde, f, mu, kTS, e, N, sigma, maxocc, method, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: dfde(:, :), f(:), mu, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), N, sigma, maxocc
      INTEGER, INTENT(IN)                                :: method
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      CHARACTER(len=*), PARAMETER                        :: routineN = 'SmearFixedDeriv'

      INTEGER                                            :: handle, i, j, Nstate
      REAL(KIND=dp)                                      :: Gsum
      REAL(KIND=dp), ALLOCATABLE                         :: gvec(:)

      CALL timeset(routineN, handle)

      ! Step 1: find mu and f
      CALL SmearFixed(f, mu, kTS, e, N, sigma, maxocc, method, estate, festate)

      ! Step 2: build g vector
      Nstate = SIZE(e)
      ALLOCATE (gvec(Nstate))
      CALL compute_gvec(gvec, f, e, mu, sigma, maxocc, Nstate, method, estate, festate)
      Gsum = accurate_sum(gvec)

      ! Step 3: assemble dfde(i,j) = -delta_{ij}*g_i + g_i*g_j/G
      IF (ABS(Gsum) > EPSILON(Gsum)) THEN
         DO j = 1, Nstate
            DO i = 1, Nstate
               dfde(i, j) = gvec(i)*gvec(j)/Gsum
            END DO
            dfde(j, j) = dfde(j, j) - gvec(j)
         END DO
      ELSE
         dfde(:, :) = 0.0_dp
      END IF

      DEALLOCATE (gvec)
      CALL timestop(handle)

   END SUBROUTINE SmearFixedDeriv

! **************************************************************************************************
!> \brief   Apply TRANSPOSE(df/de) to a vector WITHOUT forming the full N x N
!>          Jacobian.  O(N) time and O(N) memory for all four methods.
!>
!>          Exploiting the rank-1 structure of the constrained Jacobian:
!>
!>            [J^T v]_j  =  g_j * (g . v / G  -  v_j)
!>
!>          This replaces the pattern used in qs_mo_occupation:
!>            ALLOCATE(dfde(nmo,nmo))
!>            CALL SmearFixedDeriv(dfde, ...)
!>            RESULT = MATMUL(TRANSPOSE(dfde), v)
!>            DEALLOCATE(dfde)
!>          turning O(N^2) storage + O(N^2) MATMUL into O(N) throughout.
!>
!>          Currently the sole caller (qs_ot_scf do_ener) is dead code, but
!>          this routine is ready for when it is enabled.
!>
!> \param RESULT  output vector = TRANSPOSE(df/de) * v (Nstate, output)
!> \param f       occupations (output)
!> \param mu      chemical potential (output)
!> \param kTS     smearing correction (output)
!> \param e       eigenvalues (input)
!> \param N_el    target number of electrons (input)
!> \param sigma   smearing width (input)
!> \param maxocc  maximum occupation (input)
!> \param method  smearing method selector (input)
!> \param v       input vector to multiply (Nstate, input)
!> \param estate  excited state index (optional)
!> \param festate occupation of the excited state (optional)
! **************************************************************************************************
   SUBROUTINE SmearFixedDerivMV(RESULT, f, mu, kTS, e, N_el, sigma, maxocc, method, v, estate, festate)

      REAL(KIND=dp), INTENT(OUT)                         :: RESULT(:), f(:), mu, kTS
      REAL(KIND=dp), INTENT(IN)                          :: e(:), N_el, sigma, maxocc
      INTEGER, INTENT(IN)                                :: method
      REAL(KIND=dp), INTENT(IN)                          :: v(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: estate
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: festate

      CHARACTER(len=*), PARAMETER                        :: routineN = 'SmearFixedDerivMV'

      INTEGER                                            :: handle, i, Nstate
      REAL(KIND=dp)                                      :: gdotv, Gsum
      REAL(KIND=dp), ALLOCATABLE                         :: gvec(:)

      CALL timeset(routineN, handle)

      ! Step 1: find mu and f
      CALL SmearFixed(f, mu, kTS, e, N_el, sigma, maxocc, method, estate, festate)

      ! Step 2: build g vector
      Nstate = SIZE(e)
      ALLOCATE (gvec(Nstate))
      CALL compute_gvec(gvec, f, e, mu, sigma, maxocc, Nstate, method, estate, festate)
      Gsum = accurate_sum(gvec)

      ! Step 3: RESULT_j = g_j * (g.v / G - v_j)
      IF (ABS(Gsum) > EPSILON(Gsum)) THEN
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

   END SUBROUTINE SmearFixedDerivMV

END MODULE smearing_utils
