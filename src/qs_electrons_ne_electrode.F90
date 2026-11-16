!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Ne-electrode core-charge correction for constant-potential MD: adjusts the Ne
!>        pseudopotential core charge to follow the fluctuating total electron number.
!> \author Nanchen DONGFANG
! **************************************************************************************************
MODULE qs_electrons_ne_electrode
   USE IEEE_ARITHMETIC,                 ONLY: IEEE_IS_FINITE
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE external_potential_types,        ONLY: all_potential_type,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: apply_ne_core_correction

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electrons_ne_electrode'

CONTAINS

! **************************************************************************************************
!> \brief   Apply Ne core correction update from the current total electron number.
!>
!> It is called directly on force_env once per
!> MD step, right after qs_electrons_drift updates electrons%q, and before the SCF
!> call, since it changes the Ne pseudopotential and must be in effect for the SCF
!> that immediately follows.
!>
!> Receives the target total electron number and converts it to a per-Ne core
!> correction using core_corr_ne = (target_nelect - neutral_electrons_ref)/N_ne.
!> It applies the resulting value to all Ne atomic kinds in the QS subsystem by
!> updating the underlying all-electron or GTH/SGP potential zeff_correction
!> field. This effectively corresponds to adjusting the CORE_CORRECTION keyword
!> for Ne.
!>
!> The implementation traverses the QS environment, obtains the atomic kind
!> and QS kind arrays, and for each kind whose element symbol is 'Ne' it
!> updates zeff_correction := core_corr_ne and keeps the associated effective
!> charge parameters consistent with a statically parsed CORE_CORRECTION input.
!> The resulting change is reflected in the subsequent charge and electron
!> counts through the existing zeff_correction usage in QS routines.
!> \param force_env ...
!> \param target_nelect ...
!> \param neutral_electrons_ref ...
! **************************************************************************************************
   SUBROUTINE apply_ne_core_correction(force_env, target_nelect, neutral_electrons_ref)

      TYPE(force_env_type), POINTER                      :: force_env
      REAL(KIND=dp), INTENT(IN)                          :: target_nelect, neutral_electrons_ref

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_ne_core_correction'
      REAL(KIND=dp), PARAMETER :: core_corr_refresh_tol = 1.0e-12_dp

      CHARACTER(LEN=2)                                   :: elem
      INTEGER                                            :: handle, ikind, iparticle, n_ne_atoms, &
                                                            nkind, output_unit
      LOGICAL                                            :: core_corr_changed, found_ne
      REAL(KIND=dp)                                      :: base_zeff, core_corr_ne, delta_q, &
                                                            new_zeff, old_zeff, zeff_scale
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CALL timeset(routineN, handle)

      IF (.NOT. IEEE_IS_FINITE(target_nelect) .OR. target_nelect <= 0.0_dp) THEN
         CPABORT("target electron count must be positive.")
      END IF
      IF (.NOT. IEEE_IS_FINITE(neutral_electrons_ref) .OR. neutral_electrons_ref <= 0.0_dp) THEN
         CPABORT("neutral_electrons_ref must be positive.")
      END IF

      ! Obtain QS environment and atomic kind / kind sets. If this is not a
      ! QS calculation there is nothing to do.
      NULLIFY (qs_env)
      CALL force_env_get(force_env, qs_env=qs_env)
      IF (.NOT. ASSOCIATED(qs_env)) THEN
         CPABORT("Ne correction requires a QS force environment.")
      END IF

      NULLIFY (cp_subsys)
      CALL force_env_get(force_env, subsys=cp_subsys)
      IF (.NOT. ASSOCIATED(cp_subsys)) THEN
         CPABORT("Ne correction could not access CP2K subsys.")
      END IF

      NULLIFY (atomic_kind_set)
      NULLIFY (qs_kind_set)
      NULLIFY (ks_env)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, ks_env=ks_env)

      IF (.NOT. (ASSOCIATED(atomic_kind_set) .AND. ASSOCIATED(qs_kind_set))) THEN
         CPABORT("Ne correction could not access kind sets.")
      END IF

      nkind = SIZE(atomic_kind_set, 1)
      found_ne = .FALSE.
      core_corr_changed = .FALSE.
      n_ne_atoms = 0

      DO iparticle = 1, cp_subsys%particles%n_els
         atomic_kind => cp_subsys%particles%els(iparticle)%atomic_kind
         IF (ASSOCIATED(atomic_kind)) THEN
            elem = atomic_kind%element_symbol
            elem = ADJUSTL(elem)
            IF (TRIM(elem) == 'Ne' .OR. TRIM(elem) == 'NE') THEN
               n_ne_atoms = n_ne_atoms + 1
            END IF
         END IF
      END DO

      IF (n_ne_atoms <= 0) THEN
         CPABORT("Ne mode requested but no Ne atoms found.")
      END IF

      delta_q = target_nelect - neutral_electrons_ref
      core_corr_ne = delta_q/REAL(n_ne_atoms, KIND=dp)

      DO ikind = 1, nkind
         atomic_kind => atomic_kind_set(ikind)
         qs_kind => qs_kind_set(ikind)

         ! Identify Ne atomic kinds by their element symbol. Compare against
         ! the canonical two-character symbol "Ne" (and its fully uppercase
         ! variant) after left adjustment.
         elem = atomic_kind%element_symbol
         elem = ADJUSTL(elem)

         IF (TRIM(elem) == 'Ne' .OR. TRIM(elem) == 'NE') THEN
            found_ne = .TRUE.
            NULLIFY (all_potential)
            NULLIFY (gth_potential)
            NULLIFY (sgp_potential)

            CALL get_qs_kind(qs_kind, all_potential=all_potential, gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)

            IF (ASSOCIATED(all_potential)) THEN
               IF (ABS(all_potential%zeff_correction - core_corr_ne) > core_corr_refresh_tol) THEN
                  old_zeff = all_potential%zeff
                  base_zeff = old_zeff - all_potential%zeff_correction
                  new_zeff = base_zeff + core_corr_ne
                  IF (ABS(old_zeff) > EPSILON(1.0_dp)) THEN
                     zeff_scale = new_zeff/old_zeff
                     all_potential%ccore_charge = all_potential%ccore_charge*zeff_scale
                  END IF
                  all_potential%zeff = new_zeff
                  all_potential%zeff_correction = core_corr_ne
                  core_corr_changed = .TRUE.
               END IF
            END IF
            IF (ASSOCIATED(gth_potential)) THEN
               IF (ABS(gth_potential%zeff_correction - core_corr_ne) > core_corr_refresh_tol) THEN
                  old_zeff = gth_potential%zeff
                  base_zeff = old_zeff - gth_potential%zeff_correction
                  new_zeff = base_zeff + core_corr_ne
                  IF (ABS(old_zeff) > EPSILON(1.0_dp)) THEN
                     zeff_scale = new_zeff/old_zeff
                     gth_potential%ccore_charge = gth_potential%ccore_charge*zeff_scale
                     gth_potential%cerf_ppl = gth_potential%cerf_ppl*zeff_scale
                  END IF
                  gth_potential%zeff = new_zeff
                  gth_potential%zeff_correction = core_corr_ne
                  core_corr_changed = .TRUE.
               END IF
            END IF
            IF (ASSOCIATED(sgp_potential)) THEN
               IF (ABS(sgp_potential%zeff_correction - core_corr_ne) > core_corr_refresh_tol) THEN
                  old_zeff = sgp_potential%zeff
                  base_zeff = old_zeff - sgp_potential%zeff_correction
                  new_zeff = base_zeff + core_corr_ne
                  IF (ABS(old_zeff) > EPSILON(1.0_dp)) THEN
                     zeff_scale = new_zeff/old_zeff
                     sgp_potential%ccore_charge = sgp_potential%ccore_charge*zeff_scale
                  END IF
                  sgp_potential%zeff = new_zeff
                  sgp_potential%zeff_correction = core_corr_ne
                  core_corr_changed = .TRUE.
               END IF
            END IF
         END IF
      END DO

      IF (.NOT. found_ne) THEN
         CPABORT("no Ne atomic kind was found.")
      END IF

      output_unit = cp_logger_get_default_io_unit()
      IF (output_unit > 0) THEN
         WRITE (output_unit, *) " @ DRIVER MODE: Ne-electrode target_nelect = ", target_nelect, &
            " neutral_electrons = ", neutral_electrons_ref
         WRITE (output_unit, *) " @ DRIVER MODE: Ne atoms = ", n_ne_atoms, &
            " CORE_CORRECTION per Ne = ", core_corr_ne
         IF (.NOT. core_corr_changed) THEN
            WRITE (output_unit, *) " @ DRIVER MODE: Ne CORE_CORRECTION refresh skipped; change <= ", &
               core_corr_refresh_tol
         END IF
      END IF

      IF (core_corr_changed) THEN
         CALL get_qs_kind_set(qs_kind_set, total_zeff_corr=qs_env%total_zeff_corr)
      END IF

      IF (core_corr_changed .AND. ASSOCIATED(ks_env)) THEN
         CALL set_ks_env(ks_env, potential_changed=.TRUE., &
                         forces_up_to_date=.FALSE.)
      END IF

      CALL timestop(handle)

   END SUBROUTINE apply_ne_core_correction

END MODULE qs_electrons_ne_electrode
