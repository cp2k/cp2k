!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_charge_mixing_unittest
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE qs_charge_mixing,                ONLY: charge_mixing
   USE qs_density_mixing_types,         ONLY: broyden_mixing_nr,&
                                              mixing_storage_release,&
                                              mixing_storage_type

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
   DEALLOCATE (para_env)
   CALL mp_world_finalize()

END PROGRAM qs_charge_mixing_unittest
! vim: set ts=3 sw=3 tw=132 :
