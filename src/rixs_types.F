!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> *************************************************************************************************
!> \brief Define Resonant Inelastic XRAY Scattering (RIXS) control type and associated create,
!>        release, etc subroutines
!> \author BSG (02.2025)
!> *************************************************************************************************
MODULE rixs_types
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE xas_tdp_types,                   ONLY: xas_tdp_env_create,&
                                              xas_tdp_env_release,&
                                              xas_tdp_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************
!> \brief
! **************************************************************************************************
   TYPE rixs_env_type

      TYPE(xas_tdp_env_type), POINTER            :: core_state => NULL()
      TYPE(tddfpt2_valence_type), &
         POINTER                                :: valence_state => NULL()

   END TYPE rixs_env_type

! **************************************************************************************************
!> \brief Valence state coming from the qs_tddfpt2 routine
! **************************************************************************************************
   TYPE tddfpt2_valence_type

      INTEGER                                          :: nstates = 0
      TYPE(cp_fm_type), DIMENSION(:, :), &
         POINTER                                       :: evects => NULL() ! eigenvectors
      REAL(dp), DIMENSION(:), ALLOCATABLE              :: evals  ! energies

      ! entities below are coming from tddfpt_ground_state_mos type
      TYPE(cp_fm_type), DIMENSION(:), &
         POINTER                                       :: mos_occ => NULL()

   END TYPE tddfpt2_valence_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rixs_types'

   PUBLIC :: rixs_env_type, tddfpt2_valence_type, tddfpt2_valence_state_create

   PUBLIC :: rixs_env_create, rixs_env_release

CONTAINS

! **************************************************************************************************
!> \brief Creates a rixs environment type
!> \param rixs_env the type to create
! **************************************************************************************************
   SUBROUTINE rixs_env_create(rixs_env)
      TYPE(rixs_env_type), POINTER                       :: rixs_env

      ALLOCATE (rixs_env)
      NULLIFY (rixs_env%core_state)
      CALL xas_tdp_env_create(rixs_env%core_state)
      NULLIFY (rixs_env%valence_state)
      CALL tddfpt2_valence_state_create(rixs_env%valence_state)

   END SUBROUTINE rixs_env_create

! **************************************************************************************************
!> \brief Releases the rixs environment type
!> \param rixs_env the type to release
! **************************************************************************************************
   SUBROUTINE rixs_env_release(rixs_env)
      TYPE(rixs_env_type), POINTER                       :: rixs_env

      IF (ASSOCIATED(rixs_env)) THEN
         IF (ASSOCIATED(rixs_env%core_state)) THEN
            CALL xas_tdp_env_release(rixs_env%core_state)
         END IF
         IF (ASSOCIATED(rixs_env%valence_state)) THEN
            CALL tddfpt2_valence_state_release(rixs_env%valence_state)
         END IF
      END IF

      DEALLOCATE (rixs_env)

   END SUBROUTINE rixs_env_release

! **************************************************************************************************
!> \brief Creates the valence state type
!> \param valence_state ...
! **************************************************************************************************
   SUBROUTINE tddfpt2_valence_state_create(valence_state)
      TYPE(tddfpt2_valence_type), POINTER                :: valence_state

      ALLOCATE (valence_state)

      NULLIFY (valence_state%evects)

      ! entities below come from tddfpt_ground_state_mos type
      NULLIFY (valence_state%mos_occ)

   END SUBROUTINE tddfpt2_valence_state_create

! **************************************************************************************************
!> \brief Releases the valence state type
!> \param valence_state ...
! **************************************************************************************************
   SUBROUTINE tddfpt2_valence_state_release(valence_state)
      TYPE(tddfpt2_valence_type), POINTER                :: valence_state

      IF (ASSOCIATED(valence_state)) THEN
         IF (ASSOCIATED(valence_state%evects)) THEN
            CALL cp_fm_release(valence_state%evects)
         END IF
         IF (ALLOCATED(valence_state%evals)) THEN
            DEALLOCATE (valence_state%evals)
         END IF
         IF (ASSOCIATED(valence_state%mos_occ)) THEN
            CALL cp_fm_release(valence_state%mos_occ)
         END IF
      END IF

      DEALLOCATE (valence_state)

   END SUBROUTINE tddfpt2_valence_state_release

END MODULE rixs_types
