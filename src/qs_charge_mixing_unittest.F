!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_charge_mixing_unittest
   USE ieee_arithmetic, ONLY: ieee_is_finite
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_para_env_type, &
                              mp_world_finalize, &
                              mp_world_init
   USE qs_charge_mixing, ONLY: charge_mixing
   USE qs_density_mixing_types, ONLY: broyden_mixing_nr, &
                                      mixing_storage_release, &
                                      mixing_storage_type
#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE tblite_scf_mixer_type, ONLY: mixer_type
   USE tblite_scc_mixer, ONLY: new_cp2k_tblite_mixer
#endif

   IMPLICIT NONE

   REAL(KIND=dp), DIMENSION(1, 1) :: charges
   TYPE(mp_para_env_type), POINTER :: para_env
   TYPE(mixing_storage_type), POINTER :: mixing_store

   NULLIFY (mixing_store, para_env)
   ALLOCATE (para_env)
   CALL mp_world_init(para_env)

   ALLOCATE (mixing_store)
   mixing_store%ncall = 0
   mixing_store%nbuffer = 4
   mixing_store%n_simple_mix = 0
   mixing_store%nskip_mixing = 0
   mixing_store%alpha = 0.1_dp
   mixing_store%broy_w0 = 0.01_dp
   mixing_store%nat_local = 1
   mixing_store%max_shell = 1
   ALLOCATE (mixing_store%acharge(1, 1, 4), mixing_store%dacharge(1, 1, 4), &
             mixing_store%dfbroy(1, 1, 4), mixing_store%ubroy(1, 1, 4), &
             mixing_store%abroy(4, 4), mixing_store%wbroy(4), mixing_store%atlist(1))
   mixing_store%acharge = 0.0_dp
   mixing_store%dacharge = 0.0_dp
   mixing_store%dfbroy = 0.0_dp
   mixing_store%ubroy = 0.0_dp
   mixing_store%abroy = 0.0_dp
   mixing_store%wbroy = 0.0_dp
   mixing_store%atlist = 1

   charges = 1.0_dp
   CALL charge_mixing(broyden_mixing_nr, mixing_store, charges, para_env, 1)
   charges = 0.5_dp
   CALL charge_mixing(broyden_mixing_nr, mixing_store, charges, para_env, 2)

   ! Reproduce a degenerate Broyden history by supplying the previous residual again.
   charges(1, 1) = mixing_store%acharge(1, 1, 2) - mixing_store%dacharge(1, 1, 1)
   CALL charge_mixing(broyden_mixing_nr, mixing_store, charges, para_env, 3)

   IF (.NOT. ALL(ieee_is_finite(charges))) THEN
      ERROR STOP "Broyden mixing produced non-finite charges for identical residuals"
   END IF
   IF (.NOT. ALL(mixing_store%dfbroy(:, :, 1) == 0.0_dp)) THEN
      ERROR STOP "Degenerate Broyden direction was not discarded"
   END IF

   CALL mixing_storage_release(mixing_store)
   DEALLOCATE (mixing_store)
#if defined(__TBLITE)
   CALL test_tblite_history()
#endif
   DEALLOCATE (para_env)
   CALL mp_world_finalize()

#if defined(__TBLITE)
CONTAINS

! **************************************************************************************************
!> \brief A wrapped Broyden history must match a fresh replay of the retained input/output pairs.
! **************************************************************************************************
   SUBROUTINE test_tblite_history()
      INTEGER, PARAMETER                                 :: nstep = 16, nvar = 4

      INTEGER                                            :: config, i, j, memory, step
      REAL(KIND=dp)                                      :: damping, max_weight, min_weight, omega0, &
                                                            weight_factor
      REAL(KIND=dp), DIMENSION(nvar)                     :: actual, expected
      REAL(KIND=dp), DIMENSION(nvar, nstep)              :: qin, qout

      CLASS(mixer_type), ALLOCATABLE                    :: mixer, replay
      TYPE(error_type), ALLOCATABLE                      :: error

      DO step = 1, nstep
         DO i = 1, nvar
            qin(i, step) = SIN(0.37_dp*REAL(i*step, dp))
            qout(i, step) = qin(i, step) + 0.05_dp*COS(0.21_dp*REAL(i + 2*step, dp))
         END DO
      END DO
      DO config = 1, 2
         damping = 0.4_dp
         omega0 = 0.01_dp
         min_weight = 1.0_dp
         max_weight = 100000.0_dp
         weight_factor = 0.01_dp
         IF (config == 2) THEN
            damping = 0.25_dp
            omega0 = 0.1_dp
            min_weight = 0.5_dp
            max_weight = 3.0_dp
            weight_factor = 0.2_dp
         END IF
         DO memory = 1, 4
            CALL new_cp2k_tblite_mixer(mixer, memory, nvar, damping, omega0, min_weight, max_weight, weight_factor)
            DO step = 1, nstep
               CALL mixer%set(qin(:, step))
               CALL mixer%diff(qout(:, step))
               CALL mixer%next(error)
               IF (ALLOCATED(error)) ERROR STOP "Wrapped tblite mixer failed"
               CALL mixer%get(actual)

               ! Replay at most memory + 1 pairs: this builds exactly the retained
               ! secant directions without ever overwriting a history slot.
               CALL new_cp2k_tblite_mixer(replay, memory, nvar, damping, omega0, min_weight, max_weight, weight_factor)
               DO j = MAX(1, step - memory), step
                  CALL replay%set(qin(:, j))
                  CALL replay%diff(qout(:, j))
                  CALL replay%next(error)
                  IF (ALLOCATED(error)) ERROR STOP "Replayed tblite mixer failed"
               END DO
               CALL replay%get(expected)
               IF (.NOT. ALL(ieee_is_finite(actual))) ERROR STOP "Non-finite wrapped tblite history"
               IF (MAXVAL(ABS(actual - expected)) > 1.0E-10_dp) THEN
                  ERROR STOP "Wrapped tblite history differs from an unwrapped replay"
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE test_tblite_history
#endif

END PROGRAM qs_charge_mixing_unittest
! vim: set ts=3 sw=3 tw=132 :
