!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_charge_mixing

   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: get_pseudo_inverse_svd
   USE message_passing,                 ONLY: mp_sum
   USE qs_density_mixing_types,         ONLY: broyden_mixing_new_nr,&
                                              broyden_mixing_nr,&
                                              gspace_mixing_nr,&
                                              mixing_storage_type,&
                                              multisecant_mixing_nr,&
                                              pulay_mixing_nr
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_charge_mixing'

   PUBLIC :: charge_mixing

CONTAINS

! **************************************************************************************************
!> \brief  Driver for the charge mixing, calls the proper routine given the requested method
!> \param mixing_method ...
!> \param mixing_store ...
!> \param charges ...
!> \param para_env ...
!> \param iter_count ...
!> \par History
!> \author JGH
! **************************************************************************************************
   SUBROUTINE charge_mixing(mixing_method, mixing_store, charges, para_env, iter_count)
      INTEGER, INTENT(IN)                                :: mixing_method
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: iter_count

      CHARACTER(len=*), PARAMETER                        :: routineN = 'charge_mixing'

      INTEGER                                            :: handle, ia, ii, imin, inow, nbuffer, ns, &
                                                            nvec
      REAL(dp)                                           :: alpha

      CALL timeset(routineN, handle)

      IF (mixing_method >= gspace_mixing_nr) THEN
         CPASSERT(ASSOCIATED(mixing_store))
         mixing_store%ncall = mixing_store%ncall + 1
         ns = SIZE(charges, 2)
         ns = MIN(ns, mixing_store%max_shell)
         alpha = mixing_store%alpha
         nbuffer = mixing_store%nbuffer
         inow = MOD(mixing_store%ncall - 1, nbuffer) + 1
         imin = inow - 1
         IF (imin == 0) imin = nbuffer
         IF (mixing_store%ncall > nbuffer) THEN
            nvec = nbuffer
         ELSE
            nvec = mixing_store%ncall - 1
         END IF
         IF (mixing_store%ncall > 1) THEN
            ! store in/out charge difference
            DO ia = 1, mixing_store%nat_local
               ii = mixing_store%atlist(ia)
               mixing_store%dacharge(ia, 1:ns, imin) = mixing_store%acharge(ia, 1:ns, imin) - charges(ii, 1:ns)
            END DO
         END IF
         IF ((iter_count == 1) .OR. (iter_count + 1 <= mixing_store%nskip_mixing)) THEN
            ! skip mixing
            mixing_store%iter_method = "NoMix"
         ELSEIF (((iter_count + 1 - mixing_store%nskip_mixing) <= mixing_store%n_simple_mix) .OR. (nvec == 1)) THEN
            CALL mix_charges_only(mixing_store, charges, alpha, imin, ns, para_env)
            mixing_store%iter_method = "Mixing"
         ELSEIF (mixing_method == gspace_mixing_nr) THEN
            CPABORT("Kerker method not available for Charge Mixing")
         ELSEIF (mixing_method == pulay_mixing_nr) THEN
            CPABORT("Pulay method not available for Charge Mixing")
         ELSEIF (mixing_method == broyden_mixing_nr) THEN
            CALL broyden_mixing(mixing_store, charges, imin, nvec, ns, para_env)
            mixing_store%iter_method = "Broy."
         ELSEIF (mixing_method == broyden_mixing_new_nr) THEN
            CPABORT("Broyden_mixing_new method not available for Charge Mixing")
         ELSEIF (mixing_method == multisecant_mixing_nr) THEN
            CPABORT("Multisecant_mixing method not available for Charge Mixing")
         END IF

         ! store new 'input' charges
         DO ia = 1, mixing_store%nat_local
            ii = mixing_store%atlist(ia)
            mixing_store%acharge(ia, 1:ns, inow) = charges(ii, 1:ns)
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE charge_mixing

! **************************************************************************************************
!> \brief Simple charge mixing
!> \param mixing_store ...
!> \param charges ...
!> \param alpha ...
!> \param imin ...
!> \param ns ...
!> \param para_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE mix_charges_only(mixing_store, charges, alpha, imin, ns, para_env)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      INTEGER, INTENT(IN)                                :: imin, ns
      TYPE(cp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ia, ii

      charges = 0.0_dp

      DO ia = 1, mixing_store%nat_local
         ii = mixing_store%atlist(ia)
         charges(ii, 1:ns) = alpha*mixing_store%dacharge(ia, 1:ns, imin) - mixing_store%acharge(ia, 1:ns, imin)
      END DO

      CALL mp_sum(charges, para_env%group)

   END SUBROUTINE mix_charges_only

! **************************************************************************************************
!> \brief Broyden charge mixing
!> \param mixing_store ...
!> \param charges ...
!> \param inow ...
!> \param nvec ...
!> \param ns ...
!> \param para_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE broyden_mixing(mixing_store, charges, inow, nvec, ns, para_env)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      INTEGER, INTENT(IN)                                :: inow, nvec, ns
      TYPE(cp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, ia, ii, imin, j, nbuffer, nv
      REAL(KIND=dp)                                      :: alpha, maxw, minw, omega0, rskip, wdf, &
                                                            wfac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cvec, gammab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, beta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dq_last, dq_now, q_last, q_now

      CPASSERT(nvec > 1)

      nbuffer = mixing_store%nbuffer
      alpha = mixing_store%alpha
      imin = inow - 1
      IF (imin == 0) imin = nvec
      nv = nvec - 1

      ! charge vectors
      q_now => mixing_store%acharge(:, :, inow)
      q_last => mixing_store%acharge(:, :, imin)
      dq_now => mixing_store%dacharge(:, :, inow)
      dq_last => mixing_store%dacharge(:, :, imin)

      IF (nvec == nbuffer) THEN
         ! reshuffel Broyden storage n->n-1
         DO i = 1, nv - 1
            mixing_store%wbroy(i) = mixing_store%wbroy(i + 1)
            mixing_store%dfbroy(:, :, i) = mixing_store%dfbroy(:, :, i + 1)
            mixing_store%ubroy(:, :, i) = mixing_store%ubroy(:, :, i + 1)
         END DO
         DO i = 1, nv - 1
            DO j = 1, nv - 1
               mixing_store%abroy(i, j) = mixing_store%abroy(i + 1, j + 1)
            END DO
         END DO
      END IF

      omega0 = 0.01_dp
      minw = 1.0_dp
      maxw = 100000.0_dp
      wfac = 0.01_dp

      mixing_store%wbroy(nv) = SUM(dq_now(:, :)**2)
      CALL mp_sum(mixing_store%wbroy(nv), para_env%group)
      mixing_store%wbroy(nv) = SQRT(mixing_store%wbroy(nv))
      IF (mixing_store%wbroy(nv) > (wfac/maxw)) THEN
         mixing_store%wbroy(nv) = wfac/mixing_store%wbroy(nv)
      ELSE
         mixing_store%wbroy(nv) = maxw
      END IF
      IF (mixing_store%wbroy(nv) < minw) mixing_store%wbroy(nv) = minw

      ! dfbroy
      mixing_store%dfbroy(:, :, nv) = dq_now(:, :) - dq_last(:, :)
      wdf = SUM(mixing_store%dfbroy(:, :, nv)**2)
      CALL mp_sum(wdf, para_env%group)
      wdf = 1.0_dp/SQRT(wdf)
      mixing_store%dfbroy(:, :, nv) = wdf*mixing_store%dfbroy(:, :, nv)

      ! abroy matrix
      DO i = 1, nv
         wfac = SUM(mixing_store%dfbroy(:, :, i)*mixing_store%dfbroy(:, :, nv))
         CALL mp_sum(wfac, para_env%group)
         mixing_store%abroy(i, nv) = wfac
         mixing_store%abroy(nv, i) = wfac
      END DO

      ! broyden matrices
      ALLOCATE (amat(nv, nv), beta(nv, nv), cvec(nv), gammab(nv))
      DO i = 1, nv
         wfac = SUM(mixing_store%dfbroy(:, :, i)*dq_now(:, :))
         CALL mp_sum(wfac, para_env%group)
         cvec(i) = mixing_store%wbroy(i)*wfac
      END DO

      DO i = 1, nv
         DO j = 1, nv
            beta(j, i) = mixing_store%wbroy(j)*mixing_store%wbroy(i)*mixing_store%abroy(j, i)
         END DO
         beta(i, i) = beta(i, i) + omega0*omega0
      END DO

      rskip = 1.e-12_dp
      CALL get_pseudo_inverse_svd(beta, amat, rskip)
      gammab(1:nv) = MATMUL(cvec(1:nv), amat(1:nv, 1:nv))

      ! build ubroy
      mixing_store%ubroy(:, :, nv) = alpha*mixing_store%dfbroy(:, :, nv) + wdf*(q_now(:, :) - q_last(:, :))

      charges = 0.0_dp
      DO ia = 1, mixing_store%nat_local
         ii = mixing_store%atlist(ia)
         charges(ii, 1:ns) = q_now(ia, 1:ns) + alpha*dq_now(ia, 1:ns)
      END DO
      DO i = 1, nv
         DO ia = 1, mixing_store%nat_local
            ii = mixing_store%atlist(ia)
            charges(ii, 1:ns) = charges(ii, 1:ns) - mixing_store%wbroy(i)*gammab(i)*mixing_store%ubroy(ia, 1:ns, i)
         END DO
      END DO
      CALL mp_sum(charges, para_env%group)

      DEALLOCATE (amat, beta, cvec, gammab)

   END SUBROUTINE broyden_mixing

END MODULE qs_charge_mixing
