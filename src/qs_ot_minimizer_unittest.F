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
   USE qs_ot_minimizer,                 ONLY: lbfgs_curvature_damping_shift,&
                                              lbfgs_history_restart_required,&
                                              lbfgs_step_restart_required
   USE qs_scf_loop_utils,               ONLY: qs_scf_kp_ref_refresh_ready

   IMPLICIT NONE

   REAL(KIND=dp)                                      :: damping_shift, infinity, not_a_number, &
                                                         sy_damped, yy_damped

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
   IF (qs_scf_kp_ref_refresh_ready(12, 0.2E-6_dp, 1.0E-6_dp, .TRUE., .TRUE.)) THEN
      ERROR STOP "K-point REF refresh accepted a safely converged exit state"
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

END PROGRAM qs_ot_minimizer_unittest
