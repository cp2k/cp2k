!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Methods for testing / debugging.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_test

   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_gaussian,                ONLY: create_gaussian_overlap_dist_to_hermite,&
                                              create_hermite_to_cartesian
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate,&
                                              eri_mme_3c_integrate
   USE eri_mme_types,                   ONLY: eri_mme_coulomb,&
                                              eri_mme_longrange,&
                                              eri_mme_param,&
                                              eri_mme_set_potential,&
                                              eri_mme_yukawa
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_test'

   PUBLIC :: eri_mme_2c_perf_acc_test, &
             eri_mme_3c_perf_acc_test

CONTAINS
! **************************************************************************************************
!> \brief Unit test for performance and accuracy
!> \param param ...
!> \param l_max ...
!> \param zet ...
!> \param rabc ...
!> \param nrep ...
!> \param test_accuracy ...
!> \param para_env ...
!> \param iw ...
!> \param potential ...
!> \param pot_par ...
!> \param G_count ...
!> \param R_count ...
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_perf_acc_test(param, l_max, zet, rabc, nrep, test_accuracy, &
                                       para_env, iw, potential, pot_par, G_count, R_count)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: l_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: zet
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rabc
      INTEGER, INTENT(IN)                                :: nrep
      LOGICAL, INTENT(INOUT)                             :: test_accuracy
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: iw
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par
      INTEGER, INTENT(OUT), OPTIONAL                     :: G_count, R_count

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_2c_perf_acc_test', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iab, irep, izet, l, nR, nzet
      LOGICAL                                            :: acc_check
      REAL(KIND=dp)                                      :: acc, t0, t1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: time
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: I_diff, I_ref, I_test
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ht

      IF (PRESENT(G_count)) G_count = 0
      IF (PRESENT(R_count)) R_count = 0

      nzet = SIZE(zet)
      nR = SIZE(rabc, 2)

      IF (PRESENT(potential)) THEN
         CALL eri_mme_set_potential(param, potential, pot_par)
      ENDIF

      ! Calculate reference values (Exact expression in G space converged to high precision)
      IF (test_accuracy) THEN
         ht = twopi*TRANSPOSE(param%h_inv)

         ALLOCATE (I_ref(ncoset(l_max), ncoset(l_max), nR, nzet))
         I_ref(:, :, :, :) = 0.0_dp

         DO izet = 1, nzet
            DO iab = 1, nR
               CALL eri_mme_2c_integrate(param, 0, l_max, 0, l_max, zet(izet), zet(izet), rabc(:, iab), &
                                         I_ref(:, :, iab, izet), 0, 0, &
                                         normalize=.TRUE., potential=potential, pot_par=pot_par)

            ENDDO
         ENDDO
      ENDIF

      ! test performance and accuracy of MME method
      ALLOCATE (I_test(ncoset(l_max), ncoset(l_max), nR, nzet))
      ALLOCATE (I_diff(ncoset(l_max), ncoset(l_max), nR, nzet))

      ALLOCATE (time(0:l_max, nzet))
      DO l = 0, l_max
         DO izet = 1, nzet
            CALL CPU_TIME(t0)
            DO irep = 1, nrep
               DO iab = 1, nR
                  CALL eri_mme_2c_integrate(param, 0, l, 0, l, zet(izet), zet(izet), rabc(:, iab), &
                                            I_test(:, :, iab, izet), 0, 0, &
                                            G_count=G_count, R_count=R_count, &
                                            normalize=.TRUE.)
               ENDDO
            ENDDO
            CALL CPU_TIME(t1)
            time(l, izet) = t1 - t0
         ENDDO
      ENDDO

      CALL mp_sum(time, para_env%group)

      IF (test_accuracy) THEN
         I_diff(:, :, :, :) = ABS(I_test - I_ref)
      ENDIF

      IF (iw > 0) THEN
         WRITE (iw, '(T2, A)') "ERI_MME| Test results for 2c cpu time"
         WRITE (iw, '(T11, A)') "l, zet, cpu time, accuracy"

         DO l = 0, l_max
            DO izet = 1, nzet
               IF (test_accuracy) THEN
                  acc = MAXVAL(I_diff(ncoset(l - 1) + 1:ncoset(l), ncoset(l - 1) + 1:ncoset(l), :, izet))
               ELSE
                  acc = 0.0_dp
               ENDIF

               WRITE (iw, '(T11, I1, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2)') &
                  l, zet(izet), time(l, izet)/nrep, acc
            ENDDO
         ENDDO

         IF (test_accuracy) THEN
            WRITE (iw, '(/T2, A, 47X, ES9.2)') "ERI_MME| Maximum error:", &
               MAXVAL(I_diff)

            IF (param%is_ortho) THEN
               acc_check = param%err_mm + param%err_c .GE. MAXVAL(I_diff)
            ELSE
               acc_check = .TRUE.
            ENDIF

            IF (.NOT. acc_check) &
               CPABORT("Actual error greater than upper bound estimate.")

         ENDIF
      ENDIF

   END SUBROUTINE eri_mme_2c_perf_acc_test

! **************************************************************************************************
!> \brief ...
!> \param param ...
!> \param l_max ...
!> \param zet ...
!> \param rabc ...
!> \param nrep ...
!> \param nsample ...
!> \param para_env ...
!> \param iw ...
!> \param potential ...
!> \param pot_par ...
!> \param GG_count ...
!> \param GR_count ...
!> \param RR_count ...
! **************************************************************************************************
   SUBROUTINE eri_mme_3c_perf_acc_test(param, l_max, zet, rabc, nrep, nsample, &
                                       para_env, iw, potential, pot_par, GG_count, GR_count, RR_count)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: l_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: zet
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: rabc
      INTEGER, INTENT(IN)                                :: nrep
      INTEGER, INTENT(IN), OPTIONAL                      :: nsample
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: iw
      INTEGER, INTENT(IN), OPTIONAL                      :: potential
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: pot_par
      INTEGER, INTENT(OUT), OPTIONAL                     :: GG_count, GR_count, RR_count

      CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_3c_perf_acc_test', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ira, irb, irc, irep, ixyz, izeta, izetb, &
                                                            izetc, la, lb, lc, nintg, nR, ns, nzet
      REAL(KIND=dp)                                      :: t0, t1, time
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: I_test

      IF (PRESENT(GG_count)) GG_count = 0
      IF (PRESENT(GR_count)) GR_count = 0
      IF (PRESENT(RR_count)) RR_count = 0

      ns = 1
      IF (PRESENT(nsample)) ns = nsample

      nzet = SIZE(zet)
      nR = SIZE(rabc, 2)

      IF (PRESENT(potential)) THEN
         CALL eri_mme_set_potential(param, potential, pot_par)
      ENDIF

      IF (param%debug) THEN
         DO izeta = 1, nzet
         DO izetb = 1, nzet
            DO ira = 1, nR
            DO irb = 1, nR
               DO ixyz = 1, 3
                  CALL overlap_dist_expansion_test(l_max, l_max, zet(izeta), zet(izetb), &
                                                   rabc(ixyz, ira), rabc(ixyz, irb), 0.0_dp, param%debug_delta)
               ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (iw > 0) THEN
         IF (PRESENT(potential)) THEN
            SELECT CASE (potential)
            CASE (eri_mme_coulomb)
               WRITE (iw, '(/T2, A)') "ERI_MME| Potential: Coulomb"
            CASE (eri_mme_yukawa)
               WRITE (iw, '(/T2, A, ES9.2)') "ERI_MME| Potential: Yukawa with a=", pot_par
            CASE (eri_mme_longrange)
               WRITE (iw, '(/T2, A, ES9.2)') "ERI_MME| Potential: long-range Coulomb with a=", pot_par
            END SELECT
         ELSE
            WRITE (iw, '(/T2, A)') "ERI_MME| Potential: Coulomb"
         ENDIF
         WRITE (iw, '(T2, A)') "ERI_MME| Test results for 3c cpu time"
         WRITE (iw, '(T11, A)') "la, lb, lc, zeta, zetb, zetc, cpu time"
      ENDIF

      ALLOCATE (I_test(ncoset(l_max), ncoset(l_max), ncoset(l_max)))

      nintg = 0
      DO la = 0, l_max
      DO lb = 0, l_max
      DO lc = 0, l_max
         DO izeta = 1, nzet
         DO izetb = 1, nzet
         DO izetc = 1, nzet
            nintg = nintg + 1
            IF (MOD(nintg, ns) .EQ. 0) THEN
               I_test(:, :, :) = 0.0_dp
               CALL CPU_TIME(t0)
               DO irep = 1, nrep
                  DO ira = 1, nR
                  DO irb = 1, nR
                  DO irc = 1, nR
                     CALL eri_mme_3c_integrate(param, 0, la, 0, lb, 0, lc, zet(izeta), zet(izetb), zet(izetc), &
                                               rabc(:, ira), rabc(:, irb), rabc(:, irc), I_test, 0, 0, 0, &
                                               GG_count, GR_count, RR_count)
                  ENDDO
                  ENDDO
                  ENDDO
               ENDDO
               CALL CPU_TIME(t1)
               time = t1 - t0
               CALL mp_sum(time, para_env%group)
               IF (iw > 0) THEN
                  WRITE (iw, '(T11, I1, 1X, I1, 1X, I1, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2)') &
                     la, lb, lc, zet(izeta), zet(izetb), zet(izetc), time/nrep
               ENDIF
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE eri_mme_3c_perf_acc_test

! **************************************************************************************************
!> \brief check that expanding an overlap distribution of cartesian/hermite Gaussians into a
!>        lin combi of single cartesian/hermite Gaussians is correct.
!> \param l_max ...
!> \param m_max ...
!> \param zeta ...
!> \param zetb ...
!> \param R1 ...
!> \param R2 ...
!> \param r ...
!> \param tolerance ...
!> \note STATUS: tested
! **************************************************************************************************
   SUBROUTINE overlap_dist_expansion_test(l_max, m_max, zeta, zetb, R1, R2, r, tolerance)
      INTEGER, INTENT(IN)                                :: l_max, m_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta, zetb, R1, R2, r, tolerance

      INTEGER                                            :: l, m, t
      REAL(KIND=dp)                                      :: C_prod_err, H_prod_err, Rp, zetp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: C1, C2, C_ol, H1, H2, H_ol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C_prod_ref, C_prod_test, H_prod_ref, &
                                                            H_prod_test, h_to_c_1, h_to_c_2, &
                                                            h_to_c_ol
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: E_C, E_H

      zetp = zeta + zetb
      Rp = (zeta*R1 + zetb*R2)/zetp
      ALLOCATE (C1(0:l_max), H1(0:l_max))
      ALLOCATE (C2(0:m_max), H2(0:m_max))
      ALLOCATE (C_ol(0:l_max + m_max))
      ALLOCATE (H_ol(0:l_max + m_max))
      ALLOCATE (C_prod_ref(0:l_max, 0:m_max))
      ALLOCATE (C_prod_test(0:l_max, 0:m_max))
      ALLOCATE (H_prod_ref(0:l_max, 0:m_max))
      ALLOCATE (H_prod_test(0:l_max, 0:m_max))

      ALLOCATE (E_C(-1:l_max + m_max + 1, -1:l_max, -1:m_max))
      ALLOCATE (E_H(-1:l_max + m_max + 1, -1:l_max, -1:m_max))
      CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, R1, R2, 1, E_C)
      CALL create_gaussian_overlap_dist_to_hermite(l_max, m_max, zeta, zetb, R1, R2, 2, E_H)
      CALL create_hermite_to_cartesian(zetp, l_max + m_max, h_to_c_ol)
      CALL create_hermite_to_cartesian(zeta, l_max, h_to_c_1)
      CALL create_hermite_to_cartesian(zetb, m_max, h_to_c_2)

      DO t = 0, l_max + m_max
         C_ol(t) = (r - Rp)**t*EXP(-zetp*(r - Rp)**2)
      ENDDO

      DO l = 0, l_max
         C1(l) = (r - R1)**l*EXP(-zeta*(r - R1)**2)
      ENDDO
      DO m = 0, m_max
         C2(m) = (r - R2)**m*EXP(-zetb*(r - R2)**2)
      ENDDO

      H1(:) = MATMUL(TRANSPOSE(h_to_c_1(0:, 0:)), C1)
      H2(:) = MATMUL(TRANSPOSE(h_to_c_2(0:, 0:)), C2)
      H_ol(:) = MATMUL(TRANSPOSE(h_to_c_ol(0:, 0:)), C_ol)

      DO m = 0, m_max
         DO l = 0, l_max
            C_prod_ref(l, m) = C1(l)*C2(m)
            H_prod_ref(l, m) = H1(l)*H2(m)
            C_prod_test(l, m) = 0.0_dp
            H_prod_test(l, m) = 0.0_dp
            DO t = 0, l + m
               C_prod_test(l, m) = C_prod_test(l, m) + E_C(t, l, m)*H_ol(t)
               H_prod_test(l, m) = H_prod_test(l, m) + E_H(t, l, m)*H_ol(t)
            ENDDO
         ENDDO
      ENDDO

      C_prod_err = MAXVAL(ABS(C_prod_test - C_prod_ref)/(0.5_dp*(ABS(C_prod_test) + ABS(C_prod_ref)) + 1.0_dp))
      H_prod_err = MAXVAL(ABS(H_prod_test - H_prod_ref)/(0.5_dp*(ABS(C_prod_test) + ABS(C_prod_ref)) + 1.0_dp))

      CPASSERT(C_prod_err .LE. tolerance)
      CPASSERT(H_prod_err .LE. tolerance)
      MARK_USED(tolerance)

   END SUBROUTINE overlap_dist_expansion_test

END MODULE eri_mme_test
