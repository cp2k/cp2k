!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM qs_ot_minimizer_unittest
   USE kinds,                           ONLY: dp
   USE qs_ot_minimizer,                 ONLY: broyden_history_restart_required,&
                                              cg_history_restart_required,&
                                              lbfgs_curvature_damping_shift,&
                                              lbfgs_history_restart_required,&
                                              lbfgs_response_secant_parameters,&
                                              lbfgs_step_restart_required,&
                                              ot_mermin_response_assess

   IMPLICIT NONE

   REAL(KIND=dp), PARAMETER                           :: tolerance = 1.0E-14_dp

   INTEGER                                            :: i
   LOGICAL                                            :: accepted
   REAL(KIND=dp)                                      :: damping, measured, predicted, quality, ratio, &
                                                         scale, sy, weight, yy

   ! Histories must restart when they cease to describe a usable descent step.
   IF (.NOT. broyden_history_restart_required(.TRUE., .FALSE., .TRUE.)) THEN
      ERROR STOP "Broyden retained a non-descent history"
   END IF
   IF (broyden_history_restart_required(.FALSE., .TRUE., .TRUE.)) THEN
      ERROR STOP "Broyden discarded a descent history"
   END IF
   IF (.NOT. cg_history_restart_required(.TRUE., 0.5_dp + 32.0_dp*EPSILON(1.0_dp), 0.5_dp)) THEN
      ERROR STOP "CG retained unresolved free-energy history"
   END IF
   IF (cg_history_restart_required(.TRUE., 0.5_dp + 256.0_dp*EPSILON(1.0_dp), 0.5_dp)) THEN
      ERROR STOP "CG discarded resolved free-energy history"
   END IF
   IF (lbfgs_history_restart_required(100.0_dp, 1.0_dp)) THEN
      ERROR STOP "L-BFGS restarted at the gradient-growth threshold"
   END IF
   IF (.NOT. lbfgs_history_restart_required(100.0001_dp, 1.0_dp)) THEN
      ERROR STOP "L-BFGS missed excessive gradient growth"
   END IF
   IF (lbfgs_step_restart_required(1.0E-6_dp*0.1_dp, 0.1_dp)) THEN
      ERROR STOP "L-BFGS restarted at the step-collapse threshold"
   END IF
   IF (.NOT. lbfgs_step_restart_required(0.999E-6_dp*0.1_dp, 0.1_dp)) THEN
      ERROR STOP "L-BFGS retained a collapsed step"
   END IF

   ! Damping must satisfy the curvature bound for weak and negative curvature.
   IF (lbfgs_curvature_damping_shift(0.6_dp, 1.0_dp, 1.0_dp, 0.5_dp) /= 0.0_dp) THEN
      ERROR STOP "Well-conditioned L-BFGS curvature was damped"
   END IF
   DO i = -1, 1, 2
      sy = REAL(i, KIND=dp)*0.1_dp
      damping = lbfgs_curvature_damping_shift(sy, 1.0_dp, 1.0_dp, 0.5_dp)
      yy = 1.0_dp + 2.0_dp*damping*sy + damping**2
      IF (.NOT. (damping > 0.0_dp .AND. sy + damping > 0.5_dp*SQRT(yy))) THEN
         ERROR STOP "L-BFGS damping failed the relative-curvature bound"
      END IF
   END DO

   ! Keep the norm matching, Powell bound, and rejected-secant behavior together.
   CALL lbfgs_response_secant_parameters(1.0_dp, 4.0_dp, 2.0_dp, 9.0_dp, scale, weight, accepted)
   IF (.NOT. (accepted .AND. ABS(scale - 1.5_dp) <= tolerance .AND. weight == 1.0_dp)) THEN
      ERROR STOP "L-BFGS response secant changed the H0 direction norm"
   END IF
   CALL lbfgs_response_secant_parameters(0.01_dp, 4.0_dp, 2.0_dp, 9.0_dp, scale, weight, accepted)
   sy = weight*scale*0.01_dp + (1.0_dp - weight)*2.0_dp
   IF (.NOT. (accepted .AND. ABS(sy - 0.4_dp) <= tolerance .AND. weight > 0.0_dp .AND. weight < 1.0_dp)) THEN
      ERROR STOP "L-BFGS response secant failed the Powell curvature bound"
   END IF
   CALL lbfgs_response_secant_parameters(-1.0_dp, 4.0_dp, 2.0_dp, 9.0_dp, scale, weight, accepted)
   IF (accepted .OR. weight /= 0.0_dp) ERROR STOP "Rejected response secant changed H0"

   ! A useful Mermin response must lower free energy without growing the residual.
   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 0.8_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.1_dp, &
      predicted, measured, quality, ratio, accepted)
   IF (.NOT. (accepted .AND. ABS(predicted - 0.375_dp) <= tolerance .AND. &
              ABS(measured - 0.4_dp) <= tolerance .AND. ABS(quality - 0.4_dp/0.375_dp) <= tolerance .AND. &
              ABS(ratio - 0.8_dp) <= tolerance)) ERROR STOP "Useful Mermin response was rejected"
   CALL ot_mermin_response_assess( &
      10.0_dp, 9.6_dp, 1.0_dp, 1.5_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.1_dp, &
      predicted, measured, quality, ratio, accepted)
   IF (accepted) ERROR STOP "Residual-growing Mermin response was accepted"

END PROGRAM qs_ot_minimizer_unittest
