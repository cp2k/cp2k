!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_ot_minimizer_unittest
   USE ieee_arithmetic,                 ONLY: ieee_positive_inf,&
                                              ieee_quiet_nan,&
                                              ieee_value
   USE kinds,                           ONLY: dp
   USE qs_ot_minimizer,                 ONLY: &
        broyden_history_restart_required, lbfgs_curvature_damping_shift, &
        lbfgs_history_restart_required, lbfgs_response_secant_parameters, &
        lbfgs_step_restart_required, ot_mermin_response_assess, &
        ot_mermin_response_candidate_preferred, ot_mermin_response_compare, &
        ot_mermin_response_preparation_needed, ot_mermin_response_probe, &
        ot_mermin_response_shadow_followup, ot_mermin_secant_curvature
   USE qs_scf_loop_utils,               ONLY: qs_scf_kp_ref_refresh_ready,&
                                              qs_scf_kp_search_endpoint

   IMPLICIT NONE

   LOGICAL                                            :: curvature_valid, response_good, valid_response
   REAL(KIND=dp)                                      :: baseline_drop, damping_shift, infinity, &
                                                         mixed_curvature, &
                                                         measured_drop, not_a_number, predicted_drop, &
                                                         response_advantage, response_gain, response_quality, &
                                                         response_scale, response_weight, residual_ratio, &
                                                         secant_curvature, sy_damped, yy_damped

   IF (.NOT. broyden_history_restart_required(.TRUE., .FALSE., .TRUE.)) THEN
      ERROR STOP "Mermin Broyden did not restart a non-descent product history"
   END IF
   IF (.NOT. broyden_history_restart_required(.TRUE., .TRUE., .FALSE.)) THEN
      ERROR STOP "Explicit Broyden history restart was ignored"
   END IF
   IF (broyden_history_restart_required(.FALSE., .TRUE., .TRUE.)) THEN
      ERROR STOP "Broyden history restarted for a descent direction"
   END IF
   IF (broyden_history_restart_required(.TRUE., .FALSE., .FALSE.)) THEN
      ERROR STOP "Orbital-only Broyden default behavior changed"
   END IF

   damping_shift = lbfgs_curvature_damping_shift(0.6_dp, 1.0_dp, 1.0_dp, 0.5_dp)
   IF (damping_shift /= 0.0_dp) ERROR STOP "Well-conditioned L-BFGS curvature was damped"

   damping_shift = lbfgs_curvature_damping_shift(0.1_dp, 1.0_dp, 1.0_dp, 0.5_dp)
   IF (damping_shift <= 0.0_dp) ERROR STOP "Weak L-BFGS curvature was not damped"
   sy_damped = 0.1_dp + damping_shift
   yy_damped = 1.0_dp + 2.0_dp*damping_shift*0.1_dp + damping_shift**2
   IF (sy_damped <= 0.5_dp*SQRT(yy_damped)) THEN
      ERROR STOP "L-BFGS damping did not satisfy the relative-curvature bound"
   END IF

   damping_shift = lbfgs_curvature_damping_shift(-0.1_dp, 1.0_dp, 1.0_dp, 0.2_dp)
   IF (damping_shift <= 0.0_dp) ERROR STOP "Negative L-BFGS curvature was not regularized"

   IF (lbfgs_curvature_damping_shift(0.1_dp, 0.0_dp, 1.0_dp, 0.5_dp) /= 0.0_dp) THEN
      ERROR STOP "L-BFGS damping accepted a zero-length step"
   END IF

   IF (lbfgs_history_restart_required(100.0_dp, 1.0_dp)) THEN
      ERROR STOP "L-BFGS history restarted at the exact gradient-growth threshold"
   END IF
   IF (.NOT. lbfgs_history_restart_required(100.0001_dp, 1.0_dp)) THEN
      ERROR STOP "L-BFGS history did not restart above the gradient-growth threshold"
   END IF
   IF (.NOT. lbfgs_history_restart_required(400.0004_dp, 4.0_dp)) THEN
      ERROR STOP "L-BFGS gradient-growth restart is not scale invariant"
   END IF
   IF (lbfgs_history_restart_required(1.0_dp, 0.0_dp)) THEN
      ERROR STOP "L-BFGS history restarted without a valid preceding gradient"
   END IF

   IF (lbfgs_step_restart_required(1.0E-6_dp*0.1_dp, 0.1_dp)) THEN
      ERROR STOP "L-BFGS step restarted at the exact line-search collapse threshold"
   END IF
   IF (.NOT. lbfgs_step_restart_required(0.999E-6_dp*0.1_dp, 0.1_dp)) THEN
      ERROR STOP "L-BFGS step did not restart below the line-search collapse threshold"
   END IF
   IF (.NOT. lbfgs_step_restart_required(-0.999E-6_dp*2.5_dp, 2.5_dp)) THEN
      ERROR STOP "L-BFGS step restart does not use the accepted-step magnitude"
   END IF
   IF (lbfgs_step_restart_required(0.0_dp, 0.0_dp)) THEN
      ERROR STOP "L-BFGS step restarted without a valid reference step"
   END IF

   CALL lbfgs_response_secant_parameters( &
      1.0_dp, 4.0_dp, 2.0_dp, 9.0_dp, response_scale, response_weight, valid_response)
   IF (.NOT. valid_response) ERROR STOP "Valid L-BFGS response secant was rejected"
   IF (ABS(response_scale - 1.5_dp) > 1.0E-14_dp) THEN
      ERROR STOP "L-BFGS response secant did not match the H0 direction norm"
   END IF
   IF (response_weight /= 1.0_dp) ERROR STOP "Well-aligned L-BFGS response secant was damped"

   CALL lbfgs_response_secant_parameters( &
      0.01_dp, 4.0_dp, 2.0_dp, 9.0_dp, response_scale, response_weight, valid_response)
   IF (.NOT. valid_response) ERROR STOP "Dampable L-BFGS response secant was rejected"
   mixed_curvature = response_weight*response_scale*0.01_dp + &
                     (1.0_dp - response_weight)*2.0_dp
   IF (ABS(mixed_curvature - 0.4_dp) > 1.0E-14_dp) THEN
      ERROR STOP "L-BFGS response secant does not satisfy the Powell curvature bound"
   END IF
   IF (response_weight <= 0.0_dp .OR. response_weight >= 1.0_dp) THEN
      ERROR STOP "Weak L-BFGS response secant was not convexly damped"
   END IF

   CALL lbfgs_response_secant_parameters( &
      -1.0_dp, 4.0_dp, 2.0_dp, 9.0_dp, response_scale, response_weight, valid_response)
   IF (valid_response) ERROR STOP "Non-descent L-BFGS response secant was accepted"
   IF (response_weight /= 0.0_dp) ERROR STOP "Rejected L-BFGS response secant changed H0"

   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 0.8_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.1_dp, &
      predicted_drop, measured_drop, response_quality, residual_ratio, response_good)
   IF (.NOT. response_good .OR. ABS(predicted_drop - 0.375_dp) > 1.0E-14_dp .OR. &
       ABS(measured_drop - 0.4_dp) > 1.0E-14_dp .OR. &
       ABS(response_quality - 0.4_dp/0.375_dp) > 1.0E-14_dp .OR. &
       ABS(residual_ratio - 0.8_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Useful finite Mermin response sample was not accepted"
   END IF
   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 1.5_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.1_dp, &
      predicted_drop, measured_drop, response_quality, residual_ratio, response_good)
   IF (response_good) ERROR STOP "Residual-growing finite response sample was accepted"
   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 1.05_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.1_dp, &
      predicted_drop, measured_drop, response_quality, residual_ratio, response_good)
   IF (response_good) ERROR STOP "Mildly residual-growing finite response sample was accepted"
   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 0.8_dp, 0.1_dp, 10.0_dp, 1.0_dp, 0.1_dp, &
      predicted_drop, measured_drop, response_quality, residual_ratio, response_good)
   IF (response_good .OR. predicted_drop >= 0.0_dp) THEN
      ERROR STOP "Indefinite non-descent response prediction was accepted"
   END IF
   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 0.8_dp, 1.0_dp, 1.0_dp, 1.0E-10_dp, 0.1_dp, &
      predicted_drop, measured_drop, response_quality, residual_ratio, response_good)
   IF (response_good) ERROR STOP "Collapsed response line search was accepted"

   CALL ot_mermin_secant_curvature(10.0_dp, 9.6_dp, 1.0_dp, 0.5_dp, &
                                   secant_curvature, curvature_valid)
   IF (.NOT. curvature_valid .OR. ABS(secant_curvature - 0.8_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Accepted Mermin secant curvature was not recovered"
   END IF
   CALL ot_mermin_secant_curvature(10.0_dp, 9.6_dp, 1.0_dp, 0.0_dp, &
                                   secant_curvature, curvature_valid)
   IF (curvature_valid .OR. secant_curvature /= 0.0_dp) THEN
      ERROR STOP "Collapsed Mermin secant produced a curvature"
   END IF

   CALL ot_mermin_response_compare( &
      10.0_dp, 9.6_dp, 1.0_dp, 0.8_dp, 1.0_dp, 1.0_dp, 0.5_dp, &
      response_advantage, residual_ratio, response_good)
   IF (.NOT. response_good .OR. response_advantage <= 0.0_dp .OR. &
       ABS(residual_ratio - 0.8_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Improved passive finite response prediction was not retained"
   END IF
   CALL ot_mermin_response_compare( &
      10.0_dp, 9.6_dp, 1.0_dp, 0.8_dp, 1.0_dp, -1.0_dp, 0.5_dp, &
      response_advantage, residual_ratio, response_good)
   IF (response_good .OR. response_advantage >= 0.0_dp) THEN
      ERROR STOP "Worse passive finite response prediction was retained"
   END IF

   CALL ot_mermin_response_candidate_preferred( &
      1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 0.5_dp, &
      baseline_drop, predicted_drop, response_gain, response_good)
   IF (.NOT. response_good .OR. ABS(baseline_drop - 0.25_dp) > 1.0E-14_dp .OR. &
       ABS(predicted_drop - 0.375_dp) > 1.0E-14_dp .OR. response_gain <= 0.0_dp) THEN
      ERROR STOP "Better coupled Mermin candidate was not preferred"
   END IF
   CALL ot_mermin_response_candidate_preferred( &
      1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 0.5_dp, &
      baseline_drop, predicted_drop, response_gain, response_good)
   IF (response_good .OR. response_gain >= 0.0_dp) THEN
      ERROR STOP "Worse coupled Mermin candidate was preferred"
   END IF
   CALL ot_mermin_response_candidate_preferred( &
      4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 8.0_dp, 0.5_dp, &
      baseline_drop, predicted_drop, response_gain, response_good)
   IF (.NOT. response_good .OR. ABS(response_gain - 1.0_dp/3.0_dp) > 1.0E-14_dp) THEN
      ERROR STOP "Coupled Mermin candidate comparison is not scale invariant"
   END IF
   CALL ot_mermin_response_candidate_preferred( &
      1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 0.5_dp, &
      baseline_drop, predicted_drop, response_gain, response_good)
   IF (response_good) ERROR STOP "Coupled Mermin candidate ignored unresolved accepted curvature"

   IF (.NOT. ot_mermin_response_probe(.TRUE., 1.0E-3_dp, 8, 3, 0, 0)) THEN
      ERROR STOP "Uncalibrated finite response probe was not scheduled"
   END IF
   IF (ot_mermin_response_probe(.TRUE., 1.0E-3_dp, 8, 2, 0, 0)) THEN
      ERROR STOP "Finite response probe ignored passive accepted-step calibration"
   END IF
   IF (ot_mermin_response_probe(.TRUE., 1.0E-3_dp, 7, 3, 0, 0)) THEN
      ERROR STOP "Uncalibrated finite response probe was scheduled too early"
   END IF
   IF (.NOT. ot_mermin_response_probe(.TRUE., 1.0E-3_dp, 9, 3, 0, 0)) THEN
      ERROR STOP "Calibrated finite response missed its one-step eligibility window"
   END IF
   IF (.NOT. ot_mermin_response_probe(.TRUE., 1.0E-3_dp, 4, 3, 3, 0)) THEN
      ERROR STOP "Calibrated finite response probe did not use the shorter interval"
   END IF
   IF (ot_mermin_response_probe(.TRUE., 1.0E-3_dp, 8, 3, 3, 1)) THEN
      ERROR STOP "Finite response cooldown was ignored"
   END IF
   IF (ot_mermin_response_probe(.TRUE., 3.0E-3_dp, 8, 3, 3, 0)) THEN
      ERROR STOP "Finite response was probed outside its local residual regime"
   END IF

   IF (.NOT. ot_mermin_response_preparation_needed( &
       1.0_dp, 5, 3, 0, 6, .TRUE.)) THEN
      ERROR STOP "Pending finite-response shadow was not completed"
   END IF
   IF (ot_mermin_response_preparation_needed( &
       5.0E-3_dp, 6, 2, 0, 0, .FALSE.)) THEN
      ERROR STOP "Finite response was prepared outside its calibration residual regime"
   END IF
   IF (.NOT. ot_mermin_response_preparation_needed( &
       3.0E-3_dp, 4, 2, 0, 0, .FALSE.)) THEN
      ERROR STOP "Finite response did not prepare its initial shadow calibration"
   END IF
   IF (.NOT. ot_mermin_response_preparation_needed( &
       1.0E-3_dp, 6, 3, 0, 0, .FALSE.)) THEN
      ERROR STOP "Finite response did not seed the next sparse probe window"
   END IF
   IF (ot_mermin_response_preparation_needed( &
       1.0E-3_dp, 5, 3, 0, 0, .FALSE.)) THEN
      ERROR STOP "Finite response was prepared between sparse probe windows"
   END IF
   IF (ot_mermin_response_preparation_needed( &
       1.0E-3_dp, 6, 3, 0, 2, .FALSE.)) THEN
      ERROR STOP "Finite response preparation ignored an active cooldown"
   END IF

   IF (.NOT. ot_mermin_response_shadow_followup(1.0E-3_dp, 7, 3, 0, 0)) THEN
      ERROR STOP "Finite response did not retain the probe-window seed shadow"
   END IF
   IF (.NOT. ot_mermin_response_shadow_followup(1.0E-3_dp, 8, 3, 0, 0)) THEN
      ERROR STOP "Finite response did not retain the paired probe shadow"
   END IF
   IF (ot_mermin_response_shadow_followup(1.0E-3_dp, 9, 3, 0, 0)) THEN
      ERROR STOP "Finite response retained a shadow beyond its probe window"
   END IF
   IF (.NOT. ot_mermin_response_shadow_followup(3.0E-3_dp, 2, 2, 0, 0)) THEN
      ERROR STOP "Finite response stopped before completing initial calibration"
   END IF

   IF (qs_scf_kp_ref_refresh_ready(11, 1.0E-4_dp, 1.0E-6_dp, .TRUE.)) THEN
      ERROR STOP "K-point REF refresh accepted insufficient history"
   END IF
   IF (qs_scf_kp_ref_refresh_ready(12, 1.0E-4_dp, 1.0E-6_dp, .FALSE.)) THEN
      ERROR STOP "K-point REF refresh accepted a line-search trial point"
   END IF
   IF (qs_scf_kp_ref_refresh_ready(12, 0.999E-5_dp, 1.0E-6_dp, .TRUE.)) THEN
      ERROR STOP "K-point REF refresh accepted a converged residual"
   END IF
   IF (qs_scf_kp_ref_refresh_ready(12, 2.0001E-4_dp, 1.0E-6_dp, .TRUE.)) THEN
      ERROR STOP "K-point REF refresh accepted a residual outside the asymptotic window"
   END IF
   IF (.NOT. qs_scf_kp_ref_refresh_ready(12, 1.0E-4_dp, 1.0E-6_dp, .TRUE.)) THEN
      ERROR STOP "K-point REF refresh rejected a mature asymptotic state"
   END IF
   IF (qs_scf_kp_ref_refresh_ready(12, 0.75E-6_dp, 1.0E-6_dp, .TRUE.)) THEN
      ERROR STOP "K-point REF refresh entered the exit window before a physical rebuild"
   END IF
   IF (.NOT. qs_scf_kp_ref_refresh_ready(12, 0.75E-6_dp, 1.0E-6_dp, .TRUE., .TRUE.)) THEN
      ERROR STOP "K-point REF refresh rejected a mature exit-only state"
   END IF
   IF (.NOT. qs_scf_kp_ref_refresh_ready(12, 5.0E-6_dp, 1.0E-6_dp, .TRUE., .TRUE.)) THEN
      ERROR STOP "K-point REF refresh rejected a stagnant post-rebuild tail"
   END IF
   IF (qs_scf_kp_ref_refresh_ready(12, 0.2E-6_dp, 1.0E-6_dp, .TRUE., .TRUE.)) THEN
      ERROR STOP "K-point REF refresh accepted a safely converged exit state"
   END IF
   IF (.NOT. ALL([qs_scf_kp_search_endpoint("OT CG"), &
                  qs_scf_kp_search_endpoint("OT SD"), &
                  qs_scf_kp_search_endpoint("OT DIIS"), &
                  qs_scf_kp_search_endpoint("OT BROY"), &
                  qs_scf_kp_search_endpoint("OT LBFGS")])) THEN
      ERROR STOP "K-point REF refresh rejected a supported minimizer endpoint"
   END IF
   IF (qs_scf_kp_search_endpoint("OT LS") .OR. qs_scf_kp_search_endpoint("OT Sync")) THEN
      ERROR STOP "K-point REF refresh counted a non-endpoint OT iteration"
   END IF

   infinity = ieee_value(0.0_dp, ieee_positive_inf)
   not_a_number = ieee_value(0.0_dp, ieee_quiet_nan)
   IF (lbfgs_history_restart_required(infinity, 1.0_dp)) THEN
      ERROR STOP "L-BFGS history restarted for an infinite current gradient"
   END IF
   IF (lbfgs_history_restart_required(1.0_dp, infinity)) THEN
      ERROR STOP "L-BFGS history restarted for an infinite preceding gradient"
   END IF
   IF (lbfgs_history_restart_required(not_a_number, 1.0_dp)) THEN
      ERROR STOP "L-BFGS history restarted for a NaN current gradient"
   END IF
   IF (lbfgs_history_restart_required(1.0_dp, not_a_number)) THEN
      ERROR STOP "L-BFGS history restarted for a NaN preceding gradient"
   END IF
   IF (lbfgs_step_restart_required(infinity, 1.0_dp)) THEN
      ERROR STOP "L-BFGS step restarted for an infinite accepted step"
   END IF
   IF (lbfgs_step_restart_required(1.0_dp, infinity)) THEN
      ERROR STOP "L-BFGS step restarted for an infinite reference step"
   END IF
   IF (lbfgs_step_restart_required(not_a_number, 1.0_dp)) THEN
      ERROR STOP "L-BFGS step restarted for a NaN accepted step"
   END IF
   IF (lbfgs_curvature_damping_shift(not_a_number, 1.0_dp, 1.0_dp, 0.5_dp) /= 0.0_dp) THEN
      ERROR STOP "L-BFGS damping accepted NaN curvature"
   END IF
   CALL lbfgs_response_secant_parameters( &
      not_a_number, 1.0_dp, 1.0_dp, 1.0_dp, response_scale, response_weight, valid_response)
   IF (valid_response) ERROR STOP "L-BFGS response secant accepted NaN curvature"

END PROGRAM qs_ot_minimizer_unittest
