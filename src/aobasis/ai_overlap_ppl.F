!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of three-center overlap integrals over Cartesian
!>      Gaussian-type functions for the second term V(ppl) of the local
!>      part of the Goedecker pseudopotential (GTH):
!>
!>      <a|V(local)|b> = <a|V(erf) + V(ppl)|b>
!>                     = <a|V(erf)|b> + <a|V(ppl)|b>
!>                     = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!>                       (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!>                        C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!>      S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!>      C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!> \par History
!>      - Derivatives added (17.05.2002,MK)
!>      - Complete refactoring (05.2011,jhu)
!> \author Matthias Krack (04.10.2000)
! **************************************************************************************************
MODULE ai_overlap_ppl
   USE ai_oneelectron,                  ONLY: os_2center,&
                                              os_3center
   USE ai_overlap_debug,                ONLY: init_os_overlap2,&
                                              os_overlap2
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fac,&
                                              pi
   USE orbital_pointers,                ONLY: indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_overlap_ppl'

! *** Public subroutines ***

   PUBLIC :: ppl_integral, ppl_integral_ri

CONTAINS

! **************************************************************************************************
!> \brief   Calculation of three-center overlap integrals <a|c|b> over
!>           Cartesian Gaussian functions for the local part of the Goedecker
!>           pseudopotential (GTH). c is a primitive Gaussian-type function
!>           with a set of even angular momentum indices.
!>
!>           <a|V(ppl)|b> = <a| (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!>                               C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!>           zetc = alpha**2/2
!>
!> \param la_max_set ...
!> \param la_min_set ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param lb_max_set ...
!> \param lb_min_set ...
!> \param npgfb ...
!> \param rpgfb ...
!> \param zetb ...
!> \param nexp_ppl ...
!> \param alpha_ppl ...
!> \param nct_ppl ...
!> \param cexp_ppl ...
!> \param rpgfc ...
!> \param rab ...
!> \param dab ...
!> \param rac ...
!> \param dac ...
!> \param rbc ...
!> \param dbc ...
!> \param vab ...
!> \param s ...
!> \param pab ...
!> \param force_a ...
!> \param force_b ...
!> \param fs ...
!> \date    May 2011
!> \author  Juerg Hutter
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE ppl_integral(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                           lb_max_set, lb_min_set, npgfb, rpgfb, zetb, nexp_ppl, alpha_ppl, nct_ppl, cexp_ppl, rpgfc, &
                           rab, dab, rac, dac, rbc, dbc, vab, s, pab, force_a, force_b, fs)
      INTEGER, INTENT(IN)                                :: la_max_set, la_min_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: lb_max_set, lb_min_set, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfb, zetb
      INTEGER, INTENT(IN)                                :: nexp_ppl
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: alpha_ppl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nct_ppl
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cexp_ppl
      REAL(KIND=dp), INTENT(IN)                          :: rpgfc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), INTENT(IN)                          :: dab
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: dac
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rbc
      REAL(KIND=dp), INTENT(IN)                          :: dbc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: vab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: s
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: pab
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT), OPTIONAL :: force_a, force_b
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: fs

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ppl_integral'

      INTEGER                                            :: iexp, ij, ipgf, jpgf, mmax, nexp
      REAL(KIND=dp)                                      :: rho, sab, t, zetc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: auxint
      REAL(KIND=dp), DIMENSION(3)                        :: pci

      IF (PRESENT(pab)) THEN
         CPASSERT(PRESENT(force_a))
         CPASSERT(PRESENT(force_b))
         CPASSERT(PRESENT(fs))
         mmax = la_max_set + lb_max_set + 2
         force_a(:) = 0.0_dp
         force_b(:) = 0.0_dp
      ELSE
         mmax = la_max_set + lb_max_set
      END IF

      ALLOCATE (auxint(0:mmax, npgfa*npgfb))
      auxint = 0._dp

      ! *** Calculate auxiliary integrals ***

      DO ipgf = 1, npgfa
         ! *** Screening ***
         IF (rpgfa(ipgf) + rpgfc < dac) CYCLE
         DO jpgf = 1, npgfb
            ! *** Screening ***
            IF ((rpgfb(jpgf) + rpgfc < dbc) .OR. &
                (rpgfa(ipgf) + rpgfb(jpgf) < dab)) CYCLE
            ij = (ipgf - 1)*npgfb + jpgf
            rho = zeta(ipgf) + zetb(jpgf)
            pci(:) = -(zeta(ipgf)*rac(:) + zetb(jpgf)*rbc(:))/rho
            sab = EXP(-(zeta(ipgf)*zetb(jpgf)/rho*dab*dab))
            t = rho*SUM(pci(:)*pci(:))

            DO iexp = 1, nexp_ppl
               nexp = nct_ppl(iexp)
               zetc = alpha_ppl(iexp)
               CALL ppl_aux(auxint(0:mmax, ij), mmax, t, rho, nexp, cexp_ppl(:, iexp), zetc)
            END DO

            auxint(0:mmax, ij) = sab*auxint(0:mmax, ij)

         END DO
      END DO

      CALL os_3center(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                      lb_max_set, lb_min_set, npgfb, rpgfb, zetb, auxint, rpgfc, &
                      rab, dab, rac, dac, rbc, dbc, vab, s, pab, force_a, force_b, fs)

      DEALLOCATE (auxint)

   END SUBROUTINE ppl_integral
! **************************************************************************************************
!> \brief   Calculation of two-center overlap integrals <a|c> over
!>          Cartesian Gaussian functions for the local part of the Goedecker
!>          pseudopotential (GTH). c is a primitive Gaussian-type function
!>          with a set of even angular momentum indices.
!>
!>          <a|V(ppl)|b> = <a| (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!>                               C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!>          zetc = alpha**2/2
!>
!> \param la_max_set ...
!> \param la_min_set ...
!> \param npgfa ...
!> \param rpgfa ...
!> \param zeta ...
!> \param nexp_ppl ...
!> \param alpha_ppl ...
!> \param nct_ppl ...
!> \param cexp_ppl ...
!> \param rpgfc ...
!> \param rac ...
!> \param dac ...
!> \param va ...
!> \param dva ...
!> \date    December 2017
!> \author  Juerg Hutter
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE ppl_integral_ri(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                              nexp_ppl, alpha_ppl, nct_ppl, cexp_ppl, rpgfc, &
                              rac, dac, va, dva)
      INTEGER, INTENT(IN)                                :: la_max_set, la_min_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: rpgfa, zeta
      INTEGER, INTENT(IN)                                :: nexp_ppl
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: alpha_ppl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nct_ppl
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cexp_ppl
      REAL(KIND=dp), INTENT(IN)                          :: rpgfc
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      REAL(KIND=dp), INTENT(IN)                          :: dac
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: va
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: dva

      CHARACTER(len=*), PARAMETER :: routineN = 'ppl_integral_ri', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iexp, ipgf, iw, mmax, na, nexp
      INTEGER, DIMENSION(3)                              :: ani, anm, anp
      LOGICAL                                            :: debug
      REAL(KIND=dp)                                      :: oint, oref, rho, t, zetc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: auxint
      REAL(KIND=dp), DIMENSION(3)                        :: doint, doref

      debug = .FALSE.

      IF (PRESENT(dva)) THEN
         mmax = la_max_set + 1
      ELSE
         mmax = la_max_set
      END IF

      ALLOCATE (auxint(0:mmax, npgfa))
      auxint = 0._dp

      ! *** Calculate auxiliary integrals ***
      DO ipgf = 1, npgfa
         IF (rpgfa(ipgf) + rpgfc < dac) CYCLE
         rho = zeta(ipgf)
         t = rho*dac*dac

         DO iexp = 1, nexp_ppl
            nexp = nct_ppl(iexp)
            zetc = alpha_ppl(iexp)
            CALL ppl_aux(auxint(0:mmax, ipgf), mmax, t, rho, nexp, cexp_ppl(:, iexp), zetc)
         END DO

      END DO

      IF (PRESENT(dva)) THEN
         CALL os_2center(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                         auxint, rpgfc, rac, dac, va, dva)
      ELSE
         CALL os_2center(la_max_set, la_min_set, npgfa, rpgfa, zeta, &
                         auxint, rpgfc, rac, dac, va)
      END IF

      DEALLOCATE (auxint)

      IF (debug) THEN
         iw = 6
         na = 0
         DO ipgf = 1, npgfa
            IF (rpgfa(ipgf) + rpgfc < dac) THEN
               na = na + ncoset(la_max_set)
               CYCLE
            END IF
            rho = zeta(ipgf)
            DO i = ncoset(la_min_set - 1) + 1, ncoset(la_max_set)
               oref = va(na + i)
               ani(1:3) = indco(1:3, i)
               oint = ppl_ri_test(rho, ani, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl)
               ! test
               IF (ABS(oint - oref) > 1.0e-12_dp) THEN
                  WRITE (iw, '(A,3i2,i5,F10.4,2G24.12)') "PPL int error     ", ani, la_max_set, dac, oint, oref
               END IF
               IF (PRESENT(dva)) THEN
                  anp = ani + (/1, 0, 0/)
                  anm = ani - (/1, 0, 0/)
                  doint(1) = 2._dp*rho*ppl_ri_test(rho, anp, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl) &
                             - ani(1)*ppl_ri_test(rho, anm, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl)
                  anp = ani + (/0, 1, 0/)
                  anm = ani - (/0, 1, 0/)
                  doint(2) = 2._dp*rho*ppl_ri_test(rho, anp, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl) &
                             - ani(2)*ppl_ri_test(rho, anm, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl)
                  anp = ani + (/0, 0, 1/)
                  anm = ani - (/0, 0, 1/)
                  doint(3) = 2._dp*rho*ppl_ri_test(rho, anp, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl) &
                             - ani(3)*ppl_ri_test(rho, anm, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl)
                  doref(1:3) = dva(na + i, 1:3)
                  IF (ANY(ABS(doint - doref) > 1.0e-6_dp)) THEN
                     WRITE (iw, '(A,3i2,i5,F10.4,2G24.12)') " PPL dint error   ", &
                        ani, la_max_set, dac, SUM(ABS(doint)), SUM(ABS(doref))
                  END IF
               END IF
            END DO
            na = na + ncoset(la_max_set)
         END DO
      END IF

   END SUBROUTINE ppl_integral_ri

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param ani ...
!> \param rac ...
!> \param nexp_ppl ...
!> \param nct_ppl ...
!> \param alpha_ppl ...
!> \param cexp_ppl ...
!> \return ...
! **************************************************************************************************
   FUNCTION ppl_ri_test(rho, ani, rac, nexp_ppl, nct_ppl, alpha_ppl, cexp_ppl) RESULT(oint)
      REAL(KIND=dp), INTENT(IN)                          :: rho
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ani
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac
      INTEGER, INTENT(IN)                                :: nexp_ppl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nct_ppl
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: alpha_ppl
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cexp_ppl
      REAL(KIND=dp)                                      :: oint

      INTEGER                                            :: iexp, nexp, ni
      REAL(KIND=dp)                                      :: cn, zetc
      REAL(KIND=dp), DIMENSION(3)                        :: ra

      oint = 0.0_dp
      ra = 0.0_dp
      DO iexp = 1, nexp_ppl
         nexp = nct_ppl(iexp)
         zetc = alpha_ppl(iexp)
         CALL init_os_overlap2(rho, zetc, ra, -rac)
         DO ni = 1, nexp
            cn = cexp_ppl(ni, iexp)
            SELECT CASE (ni)
            CASE (1)
               oint = oint + cn*os_overlap2(ani, (/0, 0, 0/))
            CASE (2)
               oint = oint + cn*os_overlap2(ani, (/2, 0, 0/))
               oint = oint + cn*os_overlap2(ani, (/0, 2, 0/))
               oint = oint + cn*os_overlap2(ani, (/0, 0, 2/))
            CASE (3)
               oint = oint + cn*os_overlap2(ani, (/4, 0, 0/))
               oint = oint + cn*os_overlap2(ani, (/0, 4, 0/))
               oint = oint + cn*os_overlap2(ani, (/0, 0, 4/))
               oint = oint + 2.0_dp*cn*os_overlap2(ani, (/2, 2, 0/))
               oint = oint + 2.0_dp*cn*os_overlap2(ani, (/0, 2, 2/))
               oint = oint + 2.0_dp*cn*os_overlap2(ani, (/2, 0, 2/))
            CASE (4)
               oint = oint + cn*os_overlap2(ani, (/6, 0, 0/))
               oint = oint + cn*os_overlap2(ani, (/0, 6, 0/))
               oint = oint + cn*os_overlap2(ani, (/0, 0, 6/))
               oint = oint + 3.0_dp*cn*os_overlap2(ani, (/4, 2, 0/))
               oint = oint + 3.0_dp*cn*os_overlap2(ani, (/4, 0, 2/))
               oint = oint + 3.0_dp*cn*os_overlap2(ani, (/2, 4, 0/))
               oint = oint + 3.0_dp*cn*os_overlap2(ani, (/0, 4, 2/))
               oint = oint + 3.0_dp*cn*os_overlap2(ani, (/2, 0, 4/))
               oint = oint + 3.0_dp*cn*os_overlap2(ani, (/0, 2, 4/))
               oint = oint + 6.0_dp*cn*os_overlap2(ani, (/2, 2, 2/))
            CASE DEFAULT
               CPABORT("OVERLAP_PPL")
            END SELECT
         END DO
      END DO

   END FUNCTION ppl_ri_test

! **************************************************************************************************
!> \brief ...
!> \param auxint ...
!> \param mmax ...
!> \param t ...
!> \param rho ...
!> \param nexp_ppl ...
!> \param cexp_ppl ...
!> \param zetc ...
! **************************************************************************************************
   SUBROUTINE ppl_aux(auxint, mmax, t, rho, nexp_ppl, cexp_ppl, zetc)
      INTEGER, INTENT(IN)                                :: mmax
      REAL(KIND=dp), DIMENSION(0:mmax)                   :: auxint
      REAL(KIND=dp), INTENT(IN)                          :: t, rho
      INTEGER, INTENT(IN)                                :: nexp_ppl
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: cexp_ppl
      REAL(KIND=dp), INTENT(IN)                          :: zetc

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ppl_aux'

      INTEGER                                            :: i, j, ke, kp, pmax
      REAL(KIND=dp)                                      :: a2, a3, a4, cc, f, q, q2, q4, q6, rho2, &
                                                            rho3, t2, t3
      REAL(KIND=dp), DIMENSION(0:6)                      :: polder
      REAL(KIND=dp), DIMENSION(0:mmax)                   :: expder

      CPASSERT(nexp_ppl > 0)
      q = rho + zetc
      polder = 0._dp
      pmax = 0
      IF (nexp_ppl > 0) THEN
         polder(0) = polder(0) + cexp_ppl(1)
         pmax = 0
      END IF
      IF (nexp_ppl > 1) THEN
         q2 = q*q
         a2 = 0.5_dp/q2*cexp_ppl(2)
         polder(0) = polder(0) + a2*(2._dp*rho*t + 3._dp*q)
         polder(1) = polder(1) - a2*2._dp*rho
         pmax = 1
      END IF
      IF (nexp_ppl > 2) THEN
         q4 = q2*q2
         rho2 = rho*rho
         t2 = t*t
         a3 = 0.25_dp/q4*cexp_ppl(3)
         polder(0) = polder(0) + a3*(4._dp*rho2*t2 + 20._dp*rho*t*q + 15._dp*q2)
         polder(1) = polder(1) - a3*(8._dp*rho2*t + 20._dp*rho*q)
         polder(2) = polder(2) + a3*8._dp*rho2
         pmax = 2
      END IF
      IF (nexp_ppl > 3) THEN
         q6 = q4*q2
         rho3 = rho2*rho
         t3 = t2*t
         a4 = 0.125_dp/q6*cexp_ppl(4)
         polder(0) = polder(0) + a4*(8._dp*rho3*t3 + 84._dp*rho2*t2*q + 210._dp*rho*t*q2 + 105._dp*q*q2)
         polder(1) = polder(1) - a4*(24._dp*rho3*t2 + 168._dp*rho2*t*q + 210._dp*rho*q2)
         polder(2) = polder(2) + a4*(48._dp*rho3*t + 168._dp*rho2*q)
         polder(3) = polder(3) - a4*48_dp*rho3
         pmax = 3
      END IF
      IF (nexp_ppl > 4) THEN
         CPABORT("nexp_ppl > 4")
      END IF

      f = zetc/q
      cc = (pi/q)**1.5_dp*EXP(-t*f)

      IF (mmax >= 0) expder(0) = cc
      DO i = 1, mmax
         expder(i) = f*expder(i - 1)
      END DO

      DO i = 0, mmax
         DO j = 0, MIN(i, pmax)
            kp = j
            ke = i - j
            auxint(i) = auxint(i) + expder(ke)*polder(kp)*choose(i, j)
         END DO
      END DO

   END SUBROUTINE ppl_aux
! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param k ...
!> \return ...
! **************************************************************************************************
   FUNCTION choose(n, k)

      INTEGER, INTENT(IN)                                :: n, k
      REAL(KIND=dp)                                      :: choose

      IF (n >= k) THEN
         choose = REAL(NINT(fac(n)/(fac(k)*fac(n - k))), KIND=dp)
      ELSE
         choose = 0.0_dp
      ENDIF

   END FUNCTION choose
! **************************************************************************************************

END MODULE ai_overlap_ppl
