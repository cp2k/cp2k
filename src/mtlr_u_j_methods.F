!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief   Driver for self-consistent minimum tracking linear response U and J calculations.
!> \author  Ziwei Chai
!> \date    29.07.2026
!> \version 1.0
! **************************************************************************************************
MODULE mtlr_u_j_methods

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_constants,                 ONLY: atomic_guess,&
                                              restart_guess
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: do_mtlr_u_j

CONTAINS
! **************************************************************************************************
!> \brief         Driver for self-consistent MTLR U/J iteration.
!>                Each outer iteration performs:
!>                  1) one standard ENERGY SCF
!>                  2) one MTLR evaluation of U and J
!>                  3) one update of the Hubbard parameters
!>                until U and J are converged.
!>                using a method based on Lowdin charges
!>                \f[Q = S^{1/2} P S^{1/2}\f]
!>                where \b P and \b S are the density and the
!>                overlap matrix, respectively.
!> \param[in,out] force_env ...
!> \date          29.07.2026
!> \author        Ziwei Chai
!> \version       1.0
! **************************************************************************************************
   SUBROUTINE do_mtlr_u_j(force_env)

      TYPE(force_env_type), INTENT(INOUT), POINTER       :: force_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'do_mtlr_u_j'

      INTEGER                                            :: handle, ikind, k, max_mtlr_iter, n, &
                                                            nkind, output_unit, p_iter, u_iter
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: any_dft_plus_u, any_mtlr_kind, &
                                                            converged, do_reference_scf, &
                                                            wfn_restart_file_explicit
      REAL(KIND=dp) :: delta_j, delta_u, denominator_minus, denominator_plus, eps_u_j_loop, &
         fhxc_minus, fhxc_plus, intercept_minus, intercept_plus, l, std_j, std_minus, std_plus, &
         std_u, sum_trq_minus, sum_trq_minus_x_trq_minus, sum_trq_minus_x_vhxc_minus, &
         sum_trq_plus, sum_trq_plus_x_trq_plus, sum_trq_plus_x_vhxc_plus, sum_vhxc_minus, &
         sum_vhxc_plus
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: j_new, j_old, perturbation_strength, &
                                                            trq_minus, trq_plus, u_new, u_old, &
                                                            vhxc_minus, vhxc_plus
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, force_env_section

      CALL timeset(routineN, handle)

      NULLIFY (atom_list, qs_kind_set, dft_control, logger, dft_section, force_env_section)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(force_env))
      CPASSERT(ASSOCIATED(force_env%qs_env))

      CALL get_qs_env(force_env%qs_env, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, &
                      scf_control=scf_control, &
                      atomic_kind_set=atomic_kind_set)

      nkind = SIZE(atomic_kind_set)
      any_dft_plus_u = .FALSE.
      any_mtlr_kind = .FALSE.
      DO ikind = 1, nkind
         IF (.NOT. ASSOCIATED(qs_kind_set(ikind)%dft_plus_u)) CYCLE
         any_dft_plus_u = .TRUE.
         IF (qs_kind_set(ikind)%dft_plus_u%do_mtlr) THEN
            any_mtlr_kind = .TRUE.
         END IF
      END DO
      IF (.NOT. any_dft_plus_u) THEN
         CALL cp_abort(__LOCATION__, "RUN_TYPE MTLR requires at least one active "// &
                       "DFT_PLUS_U section.")
      END IF
      IF (.NOT. any_mtlr_kind) THEN
         CALL cp_abort(__LOCATION__, "RUN_TYPE MTLR requires an active "// &
                       "MINIMUM_TRACKING_LINEAR_RESPONSE subsection.")
      END IF

      CALL force_env_get(force_env, force_env_section=force_env_section)
      dft_section => section_vals_get_subs_vals(force_env_section, "DFT")
      CALL section_vals_val_get(dft_section, "WFN_RESTART_FILE_NAME", &
                                explicit=wfn_restart_file_explicit)
      IF (wfn_restart_file_explicit) THEN
         CALL cp_abort(__LOCATION__, "MTLR does not allow an explicit WFN_RESTART_FILE_NAME. "// &
                       "Remove this keyword; the reference WFN is managed internally.")
      END IF

      SELECT CASE (scf_control%density_guess)
      CASE (restart_guess)
         do_reference_scf = .TRUE.
      CASE (atomic_guess)
         do_reference_scf = .FALSE.
      CASE DEFAULT
         CPABORT("MTLR requires SCF_GUESS RESTART or SCF_GUESS ATOMIC.")
      END SELECT

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,78('='))")
         WRITE (UNIT=output_unit, FMT="(T2,A)") &
            "MTLR| SCF initialization settings"
         WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
         IF (scf_control%density_guess == restart_guess) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,T68,A12)") &
               "MTLR| SCF initial guess:", "RESTART"
            WRITE (UNIT=output_unit, FMT="(T2,A,T68,A12)") &
               "MTLR| QS extrapolation:", "USE_GUESS"
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "MTLR| A reference SCF will precede each U/J iteration."
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "MTLR| Every perturbation SCF will restart from the reference WFN."
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A,T68,A12)") &
               "MTLR| SCF initial guess:", "ATOMIC"
            WRITE (UNIT=output_unit, FMT="(T2,A,T68,A12)") &
               "MTLR| QS extrapolation:", "USE_GUESS"
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "MTLR| No separate reference SCF will be performed."
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "MTLR| Every perturbation SCF will start from an atomic guess."
         END IF
         WRITE (UNIT=output_unit, FMT="(T2,78('='))")
      END IF

      ALLOCATE (u_new(nkind))
      ALLOCATE (j_new(nkind))
      ALLOCATE (u_old(nkind))
      ALLOCATE (j_old(nkind))
      u_new(:) = 0.0_dp
      j_new(:) = 0.0_dp
      u_old(:) = 0.0_dp
      j_old(:) = 0.0_dp
      converged = .FALSE.
      eps_u_j_loop = dft_control%eps_u_j_loop
      max_mtlr_iter = dft_control%max_mtlr_iter

      ! Ensure that the DFT+U+J machinery remains active during the
      ! initial MTLR calculation, even when a parameter starts from zero.
      DO ikind = 1, nkind
         IF (ASSOCIATED(qs_kind_set(ikind)%dft_plus_u)) THEN
            IF (.NOT. qs_kind_set(ikind)%dft_plus_u%do_mtlr) CYCLE
         END IF
         IF (qs_kind_set(ikind)%dft_plus_u%u_minus_j == 0.0_dp) THEN
            qs_kind_set(ikind)%dft_plus_u%u_minus_j = EPSILON(1.0_dp)
         END IF
         IF (qs_kind_set(ikind)%dft_plus_u%hund_j == 0.0_dp) THEN
            qs_kind_set(ikind)%dft_plus_u%hund_j = EPSILON(1.0_dp)
         END IF
         j_old(ikind) = qs_kind_set(ikind)%dft_plus_u%hund_j
         u_old(ikind) = qs_kind_set(ikind)%dft_plus_u%u_minus_j + &
                        qs_kind_set(ikind)%dft_plus_u%hund_j
      END DO

      DO u_iter = 1, max_mtlr_iter

         IF (dft_control%nspins /= 2) THEN
            CPABORT("Unrestriced KS has to be used (the number of spin channels should be 2).")
         END IF

         dft_control%mtlr_dft_with_perturbation = .FALSE.

         IF (do_reference_scf) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,78('='))")
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| Starting the unperturbed reference SCF."
               WRITE (UNIT=output_unit, FMT="(T2,A,T72,I8)") &
                  "MTLR| U/J iteration:", u_iter
               WRITE (UNIT=output_unit, FMT="(T2,78('='))")
            END IF
            CALL force_env_calc_energy_force(force_env, calc_force=.FALSE.)
         END IF

         DO ikind = 1, nkind

            IF (.NOT. ASSOCIATED(qs_kind_set(ikind)%dft_plus_u)) CYCLE
            IF (.NOT. qs_kind_set(ikind)%dft_plus_u%do_mtlr) CYCLE
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list)
            IF (.NOT. ANY(atom_list == qs_kind_set(ikind)%dft_plus_u%lr_atom)) THEN
               CPABORT("INDEX_PERTURBED_ATOM does not belong to the KIND containing the MTLR section.")
            END IF
            IF (.NOT. ALLOCATED( &
                qs_kind_set(ikind)%dft_plus_u%perturbation_strength)) THEN
               CPABORT("MTLR target does not contain perturbation strengths.")
            END IF

            dft_control%mtlr_ikind = ikind

            n = SIZE(qs_kind_set(ikind)%dft_plus_u%perturbation_strength)
            l = REAL(n, dp)
            ALLOCATE (perturbation_strength(n))
            ALLOCATE (trq_plus(n))
            ALLOCATE (vhxc_plus(n))
            ALLOCATE (trq_minus(n))
            ALLOCATE (vhxc_minus(n))
            trq_plus = 0.0_dp
            trq_minus = 0.0_dp
            vhxc_plus = 0.0_dp
            vhxc_minus = 0.0_dp
            sum_trq_plus = 0.0_dp
            sum_vhxc_plus = 0.0_dp
            sum_trq_plus_x_trq_plus = 0.0_dp
            sum_trq_plus_x_vhxc_plus = 0.0_dp
            sum_trq_minus = 0.0_dp
            sum_vhxc_minus = 0.0_dp
            sum_trq_minus_x_trq_minus = 0.0_dp
            sum_trq_minus_x_vhxc_minus = 0.0_dp
            std_plus = 0.0_dp
            std_minus = 0.0_dp
            perturbation_strength(:) = qs_kind_set(ikind)%dft_plus_u%perturbation_strength(:)

            DO p_iter = 1, n

               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(/,T2,78('='))")
                  WRITE (UNIT=output_unit, FMT="(T2,A,T72,I8)") &
                     "MTLR| U/J iteration:", u_iter
                  WRITE (UNIT=output_unit, FMT="(T2,A,T74,I1,A4,I1)") &
                     "MTLR| Perturbation SCF:", p_iter, " of ", n
                  WRITE (UNIT=output_unit, FMT="(T2,A,T72,I8)") &
                     "MTLR| Target KIND index:", ikind
                  WRITE (UNIT=output_unit, FMT="(T2,A,T72,I8)") &
                     "MTLR| Target atom index:", &
                     qs_kind_set(ikind)%dft_plus_u%lr_atom
                  WRITE (UNIT=output_unit, FMT="(T2,A,T63,F14.8,A3)") &
                     "MTLR| Perturbation strength:", &
                     perturbation_strength(p_iter)*evolt, " eV"
                  WRITE (UNIT=output_unit, FMT="(T2,78('='))")
               END IF

               dft_control%perturbation_strength = perturbation_strength(p_iter)
               dft_control%mtlr_dft_with_perturbation = .TRUE.

               IF (do_reference_scf) THEN
                  CALL force_env_calc_energy_force(force_env, calc_force=.FALSE.)
               END IF

               trq_plus(p_iter) = dft_control%trq(1) + dft_control%trq(2)
               trq_minus(p_iter) = dft_control%trq(1) - dft_control%trq(2)
               vhxc_plus(p_iter) = dft_control%vhxc(1) + dft_control%vhxc(2)
               vhxc_minus(p_iter) = dft_control%vhxc(1) - dft_control%vhxc(2)

            END DO

            !Calculate all the number about vhxc and trq.
            DO k = 1, n
               sum_trq_plus_x_vhxc_plus = sum_trq_plus_x_vhxc_plus + trq_plus(k)*vhxc_plus(k)
               sum_trq_plus = sum_trq_plus + trq_plus(k)
               sum_vhxc_plus = sum_vhxc_plus + vhxc_plus(k)
               sum_trq_plus_x_trq_plus = sum_trq_plus_x_trq_plus + trq_plus(k)*trq_plus(k)
            END DO

            DO k = 1, n
               sum_trq_minus_x_vhxc_minus = sum_trq_minus_x_vhxc_minus + trq_minus(k)*vhxc_minus(k)
               sum_trq_minus = sum_trq_minus + trq_minus(k)
               sum_vhxc_minus = sum_vhxc_minus + vhxc_minus(k)
               sum_trq_minus_x_trq_minus = sum_trq_minus_x_trq_minus + trq_minus(k)*trq_minus(k)
            END DO

            denominator_plus = l*sum_trq_plus_x_trq_plus - sum_trq_plus**2.0_dp
            denominator_minus = l*sum_trq_minus_x_trq_minus - sum_trq_minus**2.0_dp
            IF (ABS(denominator_plus) < 100.0_dp*EPSILON(1.0_dp)) THEN
               CPABORT("MTLR regression is singular.")
            END IF
            IF (ABS(denominator_minus) < 100.0_dp*EPSILON(1.0_dp)) THEN
               CPABORT("MTLR regression is singular.")
            END IF
            fhxc_minus = (l*sum_trq_minus_x_vhxc_minus - sum_trq_minus*sum_vhxc_minus) &
                         /denominator_minus
            fhxc_plus = (l*sum_trq_plus_x_vhxc_plus - sum_trq_plus*sum_vhxc_plus) &
                        /denominator_plus
            intercept_minus = (sum_vhxc_minus - fhxc_minus*sum_trq_minus)/l
            intercept_plus = (sum_vhxc_plus - fhxc_plus*sum_trq_plus)/l
            DO k = 1, n
               std_minus = std_minus + (vhxc_minus(k) - (trq_minus(k)*fhxc_minus + intercept_minus))**2/(l - 2)
            END DO
            DO k = 1, n
               std_plus = std_plus + (vhxc_plus(k) - (trq_plus(k)*fhxc_plus + intercept_plus))**2/(l - 2)
            END DO
            std_minus = SQRT(std_minus)
            std_plus = SQRT(std_plus)

            u_new(ikind) = 0.5_dp*fhxc_plus
            j_new(ikind) = -0.5_dp*fhxc_minus
            delta_u = u_new(ikind) - u_old(ikind)
            delta_j = j_new(ikind) - j_old(ikind)
            std_u = 0.5_dp*std_plus
            std_j = 0.5_dp*std_minus

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,78('='))")
               WRITE (UNIT=output_unit, FMT="(T2,A,I0,A,I0,A)") &
                  "MTLR| U/J iteration ", u_iter, &
                  " results for KIND ", ikind, "."
               WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| Linear-regression data for Hubbard U"
               WRITE (UNIT=output_unit, &
                      FMT="(T2,A3,T8,A10,T22,A12,T37,A13,T52,A13,T67,A13)") &
                  "Pt", &
                  "Pert.[eV]", &
                  "trq(+)", &
                  "vhxc(+)[eV]", &
                  "vhxc(+)-fit", &
                  "vhxc(+)-resid"
               DO k = 1, n
                  WRITE (UNIT=output_unit, &
                         FMT="(T2,I3,T8,F10.4,T22,ES12.4,T37,ES13.5,"// &
                         " T52,ES13.5,T67,ES13.4)") &
                     k, &
                     perturbation_strength(k)*evolt, &
                     trq_plus(k), &
                     vhxc_plus(k)*evolt, &
                     (fhxc_plus*trq_plus(k) + intercept_plus)*evolt, &
                     (vhxc_plus(k) - &
                      (fhxc_plus*trq_plus(k) + intercept_plus))*evolt
               END DO
               WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| Linear-regression data for Hund J"
               WRITE (UNIT=output_unit, &
                      FMT="(T2,A3,T8,A10,T22,A12,T37,A13,T52,A13,T67,A13)") &
                  "Pt", &
                  "Pert.[eV]", &
                  "trq(-)", &
                  "vhxc(-)[eV]", &
                  "vhxc(-)-fit", &
                  "vhxc(-)-resid"
               DO k = 1, n
                  WRITE (UNIT=output_unit, &
                         FMT="(T2,I3,T8,F10.4,T22,ES12.4,T37,ES13.5,"// &
                         " T52,ES13.5,T67,ES13.4)") &
                     k, &
                     perturbation_strength(k)*evolt, &
                     trq_minus(k), &
                     vhxc_minus(k)*evolt, &
                     (fhxc_minus*trq_minus(k) + intercept_minus)*evolt, &
                     (vhxc_minus(k) - &
                      (fhxc_minus*trq_minus(k) + intercept_minus))*evolt
               END DO
               WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
               WRITE (UNIT=output_unit, &
                      FMT="(T2,A,T22,A16,T43,A16,T64,A16)") &
                  "MTLR| Parameter", &
                  "Old [eV]", &
                  "New [eV]", &
                  "Change [eV]"
               WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
               WRITE (UNIT=output_unit, &
                      FMT="(T2,A,T22,ES16.8,T43,ES16.8,T64,ES16.8)") &
                  "MTLR| Hubbard U", &
                  u_old(ikind)*evolt, &
                  u_new(ikind)*evolt, &
                  (u_new(ikind) - u_old(ikind))*evolt
               WRITE (UNIT=output_unit, &
                      FMT="(T2,A,T22,ES16.8,T43,ES16.8,T64,ES16.8)") &
                  "MTLR| Hund J", &
                  j_old(ikind)*evolt, &
                  j_new(ikind)*evolt, &
                  (j_new(ikind) - j_old(ikind))*evolt
               WRITE (UNIT=output_unit, FMT="(T2,78('-'))")
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES16.8,A3)") &
                  "MTLR| Absolute change in U:", &
                  ABS(u_new(ikind) - u_old(ikind))*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES16.8,A3)") &
                  "MTLR| Absolute change in J:", &
                  ABS(j_new(ikind) - j_old(ikind))*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES16.8,A3)") &
                  "MTLR| U fit residual:", &
                  std_u*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,A,T61,ES16.8,A3)") &
                  "MTLR| J fit residual:", &
                  std_j*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,78('='))")
            END IF

            DEALLOCATE (perturbation_strength)
            DEALLOCATE (trq_plus)
            DEALLOCATE (vhxc_plus)
            DEALLOCATE (trq_minus)
            DEALLOCATE (vhxc_minus)

         END DO

         DO ikind = 1, nkind
            qs_kind_set(ikind)%dft_plus_u%u_minus_j = u_new(ikind) - j_new(ikind)
            qs_kind_set(ikind)%dft_plus_u%hund_j = j_new(ikind)
         END DO

         IF (MAXVAL(ABS(u_new(:) - u_old(:))) < eps_u_j_loop .AND. &
             MAXVAL(ABS(j_new(:) - j_old(:))) < eps_u_j_loop) THEN
            converged = .TRUE.
         END IF

         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T2,78('='))")

            WRITE (UNIT=output_unit, FMT="(T2,A,I0,A)") &
               "MTLR| Iteration ", u_iter, &
               " completed for all atomic kinds."

            WRITE (UNIT=output_unit, FMT="(T2,78('-'))")

            WRITE (UNIT=output_unit, FMT="(T2,A,ES16.8,A)") &
               "MTLR| Maximum change in U:          ", &
               MAXVAL(ABS(u_new(:) - u_old(:)))*evolt, " eV"

            WRITE (UNIT=output_unit, FMT="(T2,A,ES16.8,A)") &
               "MTLR| Maximum change in J:          ", &
               MAXVAL(ABS(j_new(:) - j_old(:)))*evolt, " eV"

            WRITE (UNIT=output_unit, FMT="(T2,A,ES16.8,A)") &
               "MTLR| Convergence threshold:        ", &
               eps_u_j_loop*evolt, " eV"

            WRITE (UNIT=output_unit, FMT="(T2,78('-'))")

            IF (converged) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| U and J have converged."

               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| Both maximum changes are below the convergence threshold."

            ELSE IF (u_iter < max_mtlr_iter) THEN
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| U and J have not yet converged."

               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| Proceeding to the next linear-response U/J iteration."

            ELSE
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| U and J have not converged."

               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| The maximum number of MTLR iterations has been reached."
            END IF

            WRITE (UNIT=output_unit, FMT="(T2,78('='))")
         END IF

         IF (converged) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,78('*'))")
               WRITE (UNIT=output_unit, FMT="(T2,A,I0,A)") &
                  "MTLR| U and J converged after ", &
                  u_iter, " iterations."
               WRITE (UNIT=output_unit, FMT="(T2,78('*'),/)")
            END IF
         ELSE IF (u_iter == max_mtlr_iter) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,78('*'))")
               WRITE (UNIT=output_unit, FMT="(T2,A,I0,A)") &
                  "MTLR| U and J did not converge within the maximum of ", &
                  max_mtlr_iter, " iterations."
               WRITE (UNIT=output_unit, FMT="(T2,A)") &
                  "MTLR| Results from the final iteration will be reported."
               WRITE (UNIT=output_unit, FMT="(T2,78('*'),/)")
            END IF
         END IF

         IF (converged .OR. u_iter == max_mtlr_iter) THEN
            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,78('*'))")
               WRITE (UNIT=output_unit, FMT="(T2,A,ES16.8,A)") &
                  "MTLR| Maximum change in U:          ", &
                  MAXVAL(ABS(u_new(:) - u_old(:)))*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,A,ES16.8,A)") &
                  "MTLR| Maximum change in J:          ", &
                  MAXVAL(ABS(j_new(:) - j_old(:)))*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,A,ES16.8,A)") &
                  "MTLR| Maximum change over U and J: ", &
                  MAX(MAXVAL(ABS(u_new(:) - u_old(:))), &
                      MAXVAL(ABS(j_new(:) - j_old(:))))*evolt, " eV"
               WRITE (UNIT=output_unit, FMT="(T2,78('*'))")
               DO ikind = 1, nkind
                  IF (.NOT. qs_kind_set(ikind)%dft_plus_u%do_mtlr) CYCLE
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,I0,A)") &
                     "MTLR| Final parameters for KIND ", ikind, ":"
                  WRITE (UNIT=output_unit, FMT="(T2,A,F16.8,A)") &
                     "MTLR| Calculated Hubbard U:      ", &
                     u_new(ikind)*evolt, " eV"
                  WRITE (UNIT=output_unit, FMT="(T2,A,F16.8,A)") &
                     "MTLR| Calculated Hund J:         ", &
                     j_new(ikind)*evolt, " eV"
                  WRITE (UNIT=output_unit, FMT="(T2,A)") &
                     "MTLR| Recommended CP2K input parameters:"
                  WRITE (UNIT=output_unit, FMT="(T2,A,F16.8)") &
                     "MTLR|   U_MINUS_J [eV]           ", &
                     (u_new(ikind) - j_new(ikind))*evolt
                  WRITE (UNIT=output_unit, FMT="(T2,A,F16.8)") &
                     "MTLR|   J         [eV]           ", &
                     j_new(ikind)*evolt
               END DO
               WRITE (UNIT=output_unit, FMT="(/,T2,78('*'),/)")
            END IF
         END IF

         IF (converged) THEN
            EXIT
         END IF

         u_old(:) = u_new
         j_old(:) = j_new

      END DO

      DO ikind = 1, nkind
         IF (.NOT. qs_kind_set(ikind)%dft_plus_u%do_mtlr) CYCLE
         qs_kind_set(ikind)%dft_plus_u%u_minus_j = u_new(ikind) - j_new(ikind)
         qs_kind_set(ikind)%dft_plus_u%hund_j = j_new(ikind)
      END DO
      dft_control%mtlr_dft_with_perturbation = .FALSE.
      dft_control%perturbation_strength = 0.0_dp

      CALL force_env_calc_energy_force(force_env, calc_force=.FALSE.)

      DEALLOCATE (u_new)
      DEALLOCATE (j_new)
      DEALLOCATE (u_old)
      DEALLOCATE (j_old)

      CALL timestop(handle)

   END SUBROUTINE do_mtlr_u_j

END MODULE mtlr_u_j_methods
