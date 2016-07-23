!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Methods for testing / debugging.
!> \par History
!>       2015 09 created
!> \author Patrick Seewald
! **************************************************************************************************

MODULE eri_mme_test

   USE cp_para_types,                   ONLY: cp_para_env_type
   USE eri_mme_integrate,               ONLY: eri_mme_2c_integrate
   USE eri_mme_types,                   ONLY: eri_mme_param
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: ncoset
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_test'

   PUBLIC :: eri_mme_2c_perf_acc_test

CONTAINS
! **************************************************************************************************
!> \brief Unit test for performance and accuracy (orthorhombic cell only)
!> \param param ...
!> \param l_max ...
!> \param zet ...
!> \param rabc ...
!> \param nrep ...
!> \param test_accuracy ...
!> \param para_env ...
!> \param iw ...
!> \param G_count ...
!> \param R_count ...
! **************************************************************************************************
   SUBROUTINE eri_mme_2c_perf_acc_test(param, l_max, zet, rabc, nrep, test_accuracy, para_env, iw, G_count, R_count)
      TYPE(eri_mme_param), INTENT(INOUT)                 :: param
      INTEGER, INTENT(IN)                                :: l_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: zet
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rabc
      INTEGER, INTENT(IN)                                :: nrep
      LOGICAL, INTENT(IN)                                :: test_accuracy
      TYPE(cp_para_env_type), INTENT(IN), POINTER        :: para_env
      INTEGER, INTENT(IN)                                :: iw
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

      ! Calculate reference values (Exact expression in G space converged to high precision)
      IF (test_accuracy) THEN
         ht = twopi*TRANSPOSE(param%h_inv)

         ALLOCATE (I_ref(ncoset(l_max), ncoset(l_max), nR, nzet))
         I_ref(:, :, :, :) = 0.0_dp

         DO izet = 1, nzet
            DO iab = 1, nR
               CALL eri_mme_2c_integrate(param, 0, l_max, 0, l_max, zet(izet), zet(izet), rabc(:, iab), &
                                         I_ref(:, :, iab, izet), 0, 0, &
                                         normalize=.TRUE., exact_method=.TRUE.)

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
            time(l, izet) = t1-t0
         ENDDO
      ENDDO

      CALL mp_sum(time, para_env%group)

      IF (test_accuracy) THEN
         I_diff(:, :, :, :) = ABS(I_test-I_ref)
      ENDIF

      IF (iw > 0) THEN
         WRITE (iw, '(/T2, A)') "ERI_MME| Test results for 2c cpu time"
         WRITE (iw, '(T11, A)') "l, zet, cpu time, accuracy"

         DO l = 0, l_max
            DO izet = 1, nzet
               IF (test_accuracy) THEN
                  acc = MAXVAL(I_diff(ncoset(l-1)+1:ncoset(l), ncoset(l-1)+1:ncoset(l), :, izet))
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

            acc_check = param%err_mm+param%err_c .GE. MAXVAL(I_diff)
            IF (.NOT. acc_check) &
               CPABORT("Actual error greater than upper bound estimate.")

         ENDIF
      ENDIF

   END SUBROUTINE eri_mme_2c_perf_acc_test

END MODULE eri_mme_test
