!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief orbital transformations
!> \par History
!>      None
!> \author Joost VandeVondele (09.2002)
! **************************************************************************************************
MODULE qs_ot_minimizer

   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_scale,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_contrib,                ONLY: dbcsr_dot,&
                                              dbcsr_init_random
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_para_env_type
   USE preconditioner,                  ONLY: apply_preconditioner
   USE qs_ot,                           ONLY: qs_ot_get_derivative,&
                                              qs_ot_get_derivative_complex,&
                                              qs_ot_get_derivative_ref,&
                                              qs_ot_get_derivative_ref_complex
   USE qs_ot_types,                     ONLY: qs_ot_kpoint_preconditioner_scale,&
                                              qs_ot_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC  :: lbfgs_curvature_damping_shift, &
              lbfgs_history_restart_required, &
              lbfgs_response_secant_parameters, &
              lbfgs_step_restart_required, &
              ot_mini

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_minimizer'

CONTAINS
!
! the minimizer interface
! should present all possible modes of minimization
! these include CG SD DIIS
!
!
! IN the case of nspin != 1 we have a gradient that is distributed over different qs_ot_env.
! still things remain basically the same, since there are no constraints between the different qs_ot_env
! we only should take care that the various scalar products are taken over the full vectors.
! all the information needed and collected can be stored in the fist qs_ot_env only
! (indicating that the data type for the gradient/position and minization should be separated)
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param matrix_hc ...
!> \param matrix_hc_im ...
!> \param matrix_hc_physical occupation-weighted derivative used for slopes and convergence
!> \param matrix_hc_physical_im imaginary component of matrix_hc_physical
!> \param para_env_inter_kp communicator between distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_mini(qs_ot_env, matrix_hc, matrix_hc_im, &
                      matrix_hc_physical, matrix_hc_physical_im, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hc
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_hc_im, matrix_hc_physical, &
                                                            matrix_hc_physical_im
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_mini'

      INTEGER                                            :: handle, ispin, nspin
      LOGICAL                                            :: do_ener, do_ks, &
                                                            separate_occupation_gradient
      REAL(KIND=dp)                                      :: tmp

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      separate_occupation_gradient = PRESENT(matrix_hc_physical)
      IF (separate_occupation_gradient) THEN
         CPASSERT(qs_ot_env(1)%settings%occupation_preconditioner)
         CPASSERT(SIZE(matrix_hc_physical) == nspin)
      END IF

      qs_ot_env(1)%OT_METHOD_FULL = ""

      ! compute the gradient for the variables x
      IF (.NOT. qs_ot_env(1)%energy_only) THEN
         qs_ot_env(1)%gradient = 0.0_dp
         DO ispin = 1, nspin
            IF (do_ks) THEN
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  IF (.NOT. PRESENT(matrix_hc_im)) THEN
                     CPABORT("Complex k-point OT derivative requires imaginary H(k)*C(k).")
                  END IF
                  SELECT CASE (qs_ot_env(1)%settings%ot_algorithm)
                  CASE ("TOD")
                     IF (separate_occupation_gradient) THEN
                        CPASSERT(PRESENT(matrix_hc_physical_im))
                        CALL qs_ot_get_derivative_complex(matrix_hc(ispin)%matrix, &
                                                          matrix_hc_im(ispin)%matrix, &
                                                          qs_ot_env(ispin), &
                                                          matrix_hc_physical(ispin)%matrix, &
                                                          matrix_hc_physical_im(ispin)%matrix)
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx))
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx_im))
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                        qs_ot_env(ispin)%matrix_gx)
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx_im, &
                                        qs_ot_env(ispin)%matrix_gx_im)
                        CALL qs_ot_get_derivative_complex(matrix_hc_physical(ispin)%matrix, &
                                                          matrix_hc_physical_im(ispin)%matrix, &
                                                          qs_ot_env(ispin))
                     ELSE
                        CALL qs_ot_get_derivative_complex(matrix_hc(ispin)%matrix, &
                                                          matrix_hc_im(ispin)%matrix, &
                                                          qs_ot_env(ispin))
                     END IF
                  CASE ("REF")
                     IF (separate_occupation_gradient) THEN
                        CPASSERT(PRESENT(matrix_hc_physical_im))
                        CALL qs_ot_get_derivative_ref_complex(matrix_hc(ispin)%matrix, &
                                                              matrix_hc_im(ispin)%matrix, &
                                                              qs_ot_env(ispin), &
                                                              matrix_hc_physical(ispin)%matrix, &
                                                              matrix_hc_physical_im(ispin)%matrix)
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx))
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx_im))
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                        qs_ot_env(ispin)%matrix_gx)
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx_im, &
                                        qs_ot_env(ispin)%matrix_gx_im)
                        CALL qs_ot_get_derivative_ref_complex(matrix_hc_physical(ispin)%matrix, &
                                                              matrix_hc_physical_im(ispin)%matrix, &
                                                              qs_ot_env(ispin))
                     ELSE
                        CALL qs_ot_get_derivative_ref_complex(matrix_hc(ispin)%matrix, &
                                                              matrix_hc_im(ispin)%matrix, &
                                                              qs_ot_env(ispin))
                     END IF
                  CASE DEFAULT
                     CPABORT("Complex k-point OT derivative requires ALGORITHM STRICT or IRAC")
                  END SELECT
               ELSE
                  SELECT CASE (qs_ot_env(1)%settings%ot_algorithm)
                  CASE ("TOD")
                     CALL qs_ot_get_derivative(matrix_hc(ispin)%matrix, qs_ot_env(ispin)%matrix_x, &
                                               qs_ot_env(ispin)%matrix_sx, &
                                               qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin))
                     IF (separate_occupation_gradient) THEN
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx))
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                        qs_ot_env(ispin)%matrix_gx)
                        CALL qs_ot_get_derivative(matrix_hc_physical(ispin)%matrix, &
                                                  qs_ot_env(ispin)%matrix_x, &
                                                  qs_ot_env(ispin)%matrix_sx, &
                                                  qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin))
                     END IF
                  CASE ("REF")
                     CALL qs_ot_get_derivative_ref(matrix_hc(ispin)%matrix, &
                                                   qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_sx, &
                                                   qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin))
                     IF (separate_occupation_gradient) THEN
                        CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx))
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                        qs_ot_env(ispin)%matrix_gx)
                        CALL qs_ot_get_derivative_ref(matrix_hc_physical(ispin)%matrix, &
                                                      qs_ot_env(ispin)%matrix_x, &
                                                      qs_ot_env(ispin)%matrix_sx, &
                                                      qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin))
                     END IF
                  CASE DEFAULT
                     CPABORT("ALGORITHM NYI")
                  END SELECT
               END IF
            END IF
            ! and also the gradient along the direction
            IF (qs_ot_env(1)%use_dx) THEN
               IF (do_ks) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + tmp
                  END IF
                  IF (qs_ot_env(1)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                       qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                        qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + 0.5_dp*tmp
                     END IF
                  END IF
               END IF
               IF (do_ener) THEN
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_dx)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient + tmp
               END IF
            ELSE
               IF (do_ks) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - tmp
                  END IF
                  IF (qs_ot_env(1)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                     qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                        qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - 0.5_dp*tmp
                     END IF
                  END IF
               END IF
               IF (do_ener) THEN
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
                  qs_ot_env(1)%gradient = qs_ot_env(1)%gradient - tmp
               END IF
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gradient)
      END IF

      SELECT CASE (qs_ot_env(1)%settings%OT_METHOD)
      CASE ("CG")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT CG"
            CALL ot_new_cg_direction(qs_ot_env, para_env_inter_kp)
            qs_ot_env(1)%line_search_count = 0
         ELSE
            qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
         END IF
         CALL do_line_search(qs_ot_env)
      CASE ("SD")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT SD"
            CALL ot_new_sd_direction(qs_ot_env, para_env_inter_kp)
            qs_ot_env(1)%line_search_count = 0
         ELSE
            qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
         END IF
         CALL do_line_search(qs_ot_env)
      CASE ("DIIS")
         qs_ot_env(1)%OT_METHOD_FULL = "OT DIIS"
         CALL ot_diis_step(qs_ot_env, para_env_inter_kp)
      CASE ("BROY")
         qs_ot_env(1)%OT_METHOD_FULL = "OT BROY"
         CALL ot_broyden_step(qs_ot_env)
      CASE ("LBFG")
         IF (current_point_is_fine(qs_ot_env)) THEN
            qs_ot_env(1)%OT_METHOD_FULL = "OT LBFGS"
            CALL ot_new_lbfgs_direction(qs_ot_env, para_env_inter_kp)
            qs_ot_env(1)%line_search_count = 0
         ELSE
            qs_ot_env(1)%OT_METHOD_FULL = "OT LS"
         END IF
         CALL do_line_search(qs_ot_env)
      CASE DEFAULT
         CPABORT("OT_METHOD NYI")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE ot_mini

! **************************************************************************************************
!> \brief Sum a minimizer scalar over distributed k-point groups when requested.
!> \param para_env_inter_kp communicator between distributed k-point groups
!> \param value scalar to sum
! **************************************************************************************************
   SUBROUTINE ot_mini_sum_kpoint_scalar(para_env_inter_kp, value)
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      REAL(KIND=dp), INTENT(INOUT)                       :: value

      IF (PRESENT(para_env_inter_kp)) THEN
         IF (ASSOCIATED(para_env_inter_kp)) CALL para_env_inter_kp%sum(value)
      END IF

   END SUBROUTINE ot_mini_sum_kpoint_scalar

!
! checks if the current point is a good point for finding a new direction
! or if we should improve the line_search, if it is used
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION current_point_is_fine(qs_ot_env) RESULT(res)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      LOGICAL                                            :: res

      res = .FALSE.

      ! only if we have a gradient it can be fine
      IF (.NOT. qs_ot_env(1)%energy_only) THEN

         ! we have not yet started with the line search
         IF (qs_ot_env(1)%line_search_count == 0) THEN
            res = .TRUE.
            RETURN
         END IF

         IF (qs_ot_env(1)%line_search_might_be_done) THEN
            ! here we put the more complicated logic later
            res = .TRUE.
            RETURN
         END IF

      END IF

   END FUNCTION current_point_is_fine

!
! performs various kinds of line searches
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      SELECT CASE (qs_ot_env(1)%settings%line_search_method)
      CASE ("GOLD")
         CALL do_line_search_gold(qs_ot_env)
      CASE ("3PNT")
         CALL do_line_search_3pnt(qs_ot_env)
      CASE ("2PNT")
         IF (use_three_point_mermin_search(qs_ot_env)) THEN
            CALL do_line_search_3pnt(qs_ot_env)
         ELSE
            CALL do_line_search_2pnt(qs_ot_env)
         END IF
      CASE ("ADPT")
         CALL do_line_search_adapt(qs_ot_env)
      CASE ("NONE")
         CALL do_line_search_none(qs_ot_env)
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
   END SUBROUTINE do_line_search

! **************************************************************************************************
!> \brief Add an energy-only guard point to strongly preconditioned complex Mermin steps.
!> \param qs_ot_env OT environments
!> \return Whether to use the guarded three-point search
! **************************************************************************************************
   FUNCTION use_three_point_mermin_search(qs_ot_env) RESULT(res)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      LOGICAL                                            :: res

      INTEGER                                            :: ispin

      res = .FALSE.
      IF (.NOT. qs_ot_env(1)%settings%occupation_preconditioner) RETURN

      DO ispin = 1, SIZE(qs_ot_env)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            res = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION use_three_point_mermin_search

! **************************************************************************************************
!> \brief moves x adding the right amount (ds) of the gradient or search direction
!> \param ds ...
!> \param qs_ot_env ...
!> \par History
!>      08.2004 created [ Joost VandeVondele ] copied here from a larger number of subroutines
! **************************************************************************************************
   SUBROUTINE take_step(ds, qs_ot_env)
      REAL(KIND=dp)                                      :: ds
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'take_step'

      INTEGER                                            :: handle, ispin, nspin
      LOGICAL                                            :: do_ener, do_ks

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener

      ! now update x to take into account this new step
      ! either dx or -gx is the direction to use
      IF (qs_ot_env(1)%use_dx) THEN
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_dx, &
                              alpha_scalar=1.0_dp, beta_scalar=ds)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, qs_ot_env(ispin)%matrix_dx_im, &
                                 alpha_scalar=1.0_dp, beta_scalar=ds)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, qs_ot_env(ispin)%rot_mat_dx, &
                                 alpha_scalar=1.0_dp, beta_scalar=ds)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, qs_ot_env(ispin)%rot_mat_dx_im, &
                                    alpha_scalar=1.0_dp, beta_scalar=ds)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x = qs_ot_env(ispin)%ener_x + ds*qs_ot_env(ispin)%ener_dx
            END DO
         END IF
      ELSE
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_gx, &
                              alpha_scalar=1.0_dp, beta_scalar=-ds)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, qs_ot_env(ispin)%matrix_gx_im, &
                                 alpha_scalar=1.0_dp, beta_scalar=-ds)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, qs_ot_env(ispin)%rot_mat_gx, &
                                 alpha_scalar=1.0_dp, beta_scalar=-ds)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, qs_ot_env(ispin)%rot_mat_gx_im, &
                                    alpha_scalar=1.0_dp, beta_scalar=-ds)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x = qs_ot_env(ispin)%ener_x - ds*qs_ot_env(ispin)%ener_gx
            END DO
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE take_step

! implements a golden ratio search as a robust way of minimizing
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_gold(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_gold'
      REAL(KIND=dp), PARAMETER                           :: gold_sec = 0.3819_dp

      INTEGER                                            :: count, handle
      REAL(KIND=dp)                                      :: ds

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      IF (count + 1 > SIZE(qs_ot_env(1)%OT_pos)) THEN
         ! should not happen, we pass with a warning first
         ! you can increase the size of OT_pos and the like in qs_ot_env
         CPABORT("MAX ITER EXCEEDED : FATAL")
      END IF

      IF (qs_ot_env(1)%line_search_count == 1) THEN
         qs_ot_env(1)%line_search_left = 1
         qs_ot_env(1)%line_search_right = 0
         qs_ot_env(1)%line_search_mid = 1
         qs_ot_env(1)%ot_pos(1) = 0.0_dp
         qs_ot_env(1)%ot_energy(1) = qs_ot_env(1)%etotal
         qs_ot_env(1)%ot_pos(2) = qs_ot_env(1)%ds_min/gold_sec
      ELSE
         qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
         ! it's essentially a book keeping game.
         ! keep left on the left, keep (bring) right on the right
         ! and mid in between these two
         IF (qs_ot_env(1)%line_search_right == 0) THEN ! we do not yet have the right bracket
            IF (qs_ot_env(1)%ot_energy(count - 1) < qs_ot_env(1)%ot_energy(count)) THEN
               qs_ot_env(1)%line_search_right = count
               qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) + &
                                                (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) - &
                                                 qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid))*gold_sec
            ELSE
               qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
               qs_ot_env(1)%line_search_mid = count
               qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ot_pos(count)/gold_sec ! expand
            END IF
         ELSE
            ! first determine where we are and construct the new triplet
            IF (qs_ot_env(1)%ot_pos(count) < qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) THEN
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_right = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_left = count
               END IF
            ELSE
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_right = count
               END IF
            END IF
            ! now find the new point in the largest section
            IF ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                 - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) > &
                (qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                 - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left))) THEN
               qs_ot_env(1)%ot_pos(count + 1) = &
                  qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) + &
                  gold_sec*(qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                            - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid))
            ELSE
               qs_ot_env(1)%ot_pos(count + 1) = &
                  qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left) + &
                  gold_sec*(qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                            - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left))
            END IF
            ! check for termination
            IF (((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_right) &
                  - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) < &
                 qs_ot_env(1)%ds_min*qs_ot_env(1)%settings%gold_target) .AND. &
                ((qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid) &
                  - qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_left)) < &
                 qs_ot_env(1)%ds_min*qs_ot_env(1)%settings%gold_target)) THEN
               qs_ot_env(1)%energy_only = .FALSE.
               qs_ot_env(1)%line_search_might_be_done = .TRUE.
            END IF
         END IF
      END IF
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_gold

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_3pnt(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_3pnt'

      INTEGER                                            :: best_index, count, handle, i
      REAL(KIND=dp)                                      :: best_energy, denom, ds, fa, fb, fc, nom, &
                                                            pos, tol, val, xa, xb, xc

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      ! a three point interpolation based on the energy
      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
      SELECT CASE (count)
      CASE (1)
         qs_ot_env(1)%ot_pos(count) = 0.0_dp
         qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ds_min*0.8_dp
      CASE (2)
         IF (qs_ot_env(1)%OT_energy(count) > qs_ot_env(1)%OT_energy(count - 1)) THEN
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%ds_min*0.5_dp
         ELSE
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%ds_min*1.4_dp
         END IF
      CASE (3)
         xa = qs_ot_env(1)%OT_pos(1)
         xb = qs_ot_env(1)%OT_pos(2)
         xc = qs_ot_env(1)%OT_pos(3)
         fa = qs_ot_env(1)%OT_energy(1)
         fb = qs_ot_env(1)%OT_energy(2)
         fc = qs_ot_env(1)%OT_energy(3)
         nom = (xb - xa)**2*(fb - fc) - (xb - xc)**2*(fb - fa)
         denom = (xb - xa)*(fb - fc) - (xb - xc)*(fb - fa)
         IF (ABS(denom) <= 1.0E-18_dp*MAX(ABS(fb - fc), ABS(fb - fa))) THEN
            pos = xb
         ELSE
            pos = xb - 0.5_dp*nom/denom ! position of the stationary point
         END IF
         val = (pos - xa)*(pos - xb)*fc/((xc - xa)*(xc - xb)) + &
               (pos - xb)*(pos - xc)*fa/((xa - xb)*(xa - xc)) + &
               (pos - xc)*(pos - xa)*fb/((xb - xc)*(xb - xa))
         best_index = 1
         DO i = 2, count
            IF (qs_ot_env(1)%OT_energy(i) < qs_ot_env(1)%OT_energy(best_index)) best_index = i
         END DO
         best_energy = qs_ot_env(1)%OT_energy(best_index)
         tol = 10.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, ABS(best_energy))
         IF (use_three_point_mermin_search(qs_ot_env) .AND. best_index /= 1 .AND. &
             val >= best_energy - tol) THEN
            qs_ot_env(1)%OT_pos(count + 1) = qs_ot_env(1)%OT_pos(best_index)
         ELSE IF (val < fa .AND. val <= fb .AND. val <= fc) THEN ! OK, we go to a minimum
            ! we take a guard against too large steps
            qs_ot_env(1)%OT_pos(count + 1) = MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:3))*0.01_dp, &
                                                 MIN(pos, MAXVAL(qs_ot_env(1)%OT_pos(1:3))*4.0_dp))
         ELSE ! just take an extended step
            qs_ot_env(1)%OT_pos(count + 1) = MAXVAL(qs_ot_env(1)%OT_pos(1:3))*2.0_dp
         END IF
         qs_ot_env(1)%energy_only = .FALSE.
         qs_ot_env(1)%line_search_might_be_done = .TRUE.
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_3pnt

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_2pnt(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_2pnt'

      INTEGER                                            :: count, handle
      REAL(KIND=dp)                                      :: a, b, c, ds, pos, val, x0, x1

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      ! a three point interpolation based on the energy
      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
      SELECT CASE (count)
      CASE (1)
         qs_ot_env(1)%ot_pos(count) = 0.0_dp
         qs_ot_env(1)%ot_grad(count) = qs_ot_env(1)%gradient
         qs_ot_env(1)%ot_pos(count + 1) = qs_ot_env(1)%ds_min*1.0_dp
      CASE (2)
         x0 = 0.0_dp
         c = qs_ot_env(1)%ot_energy(1)
         b = qs_ot_env(1)%ot_grad(1)
         x1 = qs_ot_env(1)%ot_pos(2)
         a = (qs_ot_env(1)%ot_energy(2) - b*x1 - c)/(x1**2)
         IF (a <= 0.0_dp) a = 1.0E-15_dp
         pos = -b/(2.0_dp*a)
         val = a*pos**2 + b*pos + c
         qs_ot_env(1)%energy_only = .FALSE.
         qs_ot_env(1)%line_search_might_be_done = .TRUE.
         IF (val < qs_ot_env(1)%ot_energy(1) .AND. val <= qs_ot_env(1)%ot_energy(2)) THEN
            ! we go to a minimum, but ...
            ! we take a guard against too large steps
            qs_ot_env(1)%OT_pos(count + 1) = MAX(MAXVAL(qs_ot_env(1)%OT_pos(1:2))*0.01_dp, &
                                                 MIN(pos, MAXVAL(qs_ot_env(1)%OT_pos(1:2))*4.0_dp))
         ELSE ! just take an extended step
            qs_ot_env(1)%OT_pos(count + 1) = MAXVAL(qs_ot_env(1)%OT_pos(1:2))*2.0_dp
         END IF
      CASE DEFAULT
         CPABORT("NYI")
      END SELECT
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_2pnt

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_adapt(qs_ot_env)

      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CHARACTER(len=*), PARAMETER :: routineN = 'do_line_search_adapt'
      REAL(KIND=dp), PARAMETER                           :: grow_factor = 2.0_dp, &
                                                            shrink_factor = 0.5_dp

      INTEGER                                            :: count, handle, il, im, ir
      REAL(KIND=dp)                                      :: a, b, c, denom, ds, el, em, er, &
                                                            step_size, xl, xm, xr

      CALL timeset(routineN, handle)

      qs_ot_env(1)%line_search_count = qs_ot_env(1)%line_search_count + 1
      count = qs_ot_env(1)%line_search_count
      qs_ot_env(1)%line_search_might_be_done = .FALSE.
      qs_ot_env(1)%energy_only = .TRUE.

      IF (count + 1 > SIZE(qs_ot_env(1)%OT_pos)) THEN
         ! should not happen, we pass with a warning first
         ! you can increase the size of OT_pos and the like in qs_ot_env
         CPABORT("MAX ITER EXCEEDED : FATAL")
      END IF

      ! Perform an adaptive linesearch
      IF (qs_ot_env(1)%line_search_count == 1) THEN
         qs_ot_env(1)%line_search_left = 1
         qs_ot_env(1)%line_search_right = 0
         qs_ot_env(1)%line_search_mid = 1
         qs_ot_env(1)%ot_Pos(1) = 0.0_dp
         qs_ot_env(1)%ot_energy(1) = qs_ot_env(1)%etotal
         qs_ot_env(1)%ot_Pos(2) = qs_ot_env(1)%ds_min*grow_factor
      ELSE
         qs_ot_env(1)%ot_energy(count) = qs_ot_env(1)%etotal
         ! it's essentially a book keeping game.
         ! keep left on the left, keep (bring) right on the right
         ! and mid in between these two
         IF (qs_ot_env(1)%line_search_right == 0) THEN ! we do not yet have the right bracket
            IF (qs_ot_env(1)%ot_energy(count - 1) < qs_ot_env(1)%ot_energy(count)) THEN
               qs_ot_env(1)%line_search_right = count
               qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(qs_ot_env(1)%line_search_mid) + &
                                                (qs_ot_env(1)%ot_Pos(qs_ot_env(1)%line_search_right) - &
                                                 qs_ot_env(1)%ot_Pos(qs_ot_env(1)%line_search_mid))*shrink_factor
            ELSE
               ! expand further
               qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
               qs_ot_env(1)%line_search_mid = count
               qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(count)*grow_factor
            END IF
         ELSE
            ! first determine where we are and construct the new triplet
            IF (qs_ot_env(1)%ot_pos(count) < qs_ot_env(1)%ot_pos(qs_ot_env(1)%line_search_mid)) THEN
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_right = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_left = count
               END IF
            ELSE
               IF (qs_ot_env(1)%ot_energy(count) < qs_ot_env(1)%ot_energy(qs_ot_env(1)%line_search_mid)) THEN
                  qs_ot_env(1)%line_search_left = qs_ot_env(1)%line_search_mid
                  qs_ot_env(1)%line_search_mid = count
               ELSE
                  qs_ot_env(1)%line_search_right = count
               END IF
            END IF
            il = qs_ot_env(1)%line_search_left
            im = qs_ot_env(1)%line_search_mid
            ir = qs_ot_env(1)%line_search_right
            xl = qs_ot_env(1)%OT_pos(il)
            xm = qs_ot_env(1)%OT_pos(im)
            xr = qs_ot_env(1)%OT_pos(ir)
            el = qs_ot_env(1)%ot_energy(il)
            em = qs_ot_env(1)%ot_energy(im)
            er = qs_ot_env(1)%ot_energy(ir)
            IF (em < el) THEN
               IF (er < em) THEN
                  !extend search
                  qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(ir)*grow_factor
               ELSE
                  ! Cramer's rule
                  denom = (xl - xm)*(xl - xr)*(xm - xr)
                  a = (xr*(em - el) + xm*(el - er) + xl*(er - em))/denom
                  b = (xr**2*(el - em) + xm**2*(er - el) + xl**2*(em - er))/denom
                  c = (xm*xr*(xm - xr)*el + xr*xl*(xr - xl)*em + xr*xm*(xr - xm)*er)/denom

                  IF (ABS(a) /= 0.0_dp) THEN
                     step_size = -b/(2.0_dp*a)
                  ELSE
                     step_size = 0.0_dp
                  END IF
                  CPASSERT(step_size >= 0.0_dp)
                  qs_ot_env(1)%ot_Pos(count + 1) = step_size
                  qs_ot_env(1)%line_search_might_be_done = .TRUE.
                  qs_ot_env(1)%energy_only = .FALSE.
               END IF
            ELSE
               ! contract search
               qs_ot_env(1)%ot_Pos(count + 1) = qs_ot_env(1)%ot_Pos(im) + &
                                                (qs_ot_env(1)%ot_Pos(ir) - qs_ot_env(1)%ot_Pos(im))*shrink_factor
            END IF

         END IF
      END IF
      ds = qs_ot_env(1)%OT_pos(count + 1) - qs_ot_env(1)%OT_pos(count)
      qs_ot_env(1)%ds_min = qs_ot_env(1)%OT_pos(count + 1)

      CALL take_step(ds, qs_ot_env)

      CALL timestop(handle)

   END SUBROUTINE do_line_search_adapt

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
! **************************************************************************************************
   SUBROUTINE do_line_search_none(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      CALL take_step(qs_ot_env(1)%ds_min, qs_ot_env)

   END SUBROUTINE do_line_search_none

!
! creates a new SD direction, using the preconditioner if associated
! also updates the gradient for line search
!

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param para_env_inter_kp communicator between distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_new_sd_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_sd_direction'

      INTEGER                                            :: handle, ispin, itmp, k, n, nener, nspin
      LOGICAL                                            :: do_ener, do_ks
      REAL(KIND=dp)                                      :: channel_gnorm, nvariables, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

!***SCP

      nspin = SIZE(qs_ot_env)
      logger => cp_get_default_logger()
      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener

      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
         IF (.NOT. qs_ot_env(1)%use_dx) CPABORT("use dx")
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                               qs_ot_env(ispin)%matrix_preconditioned_gx_im, &
                                               qs_ot_env(ispin)%matrix_dx, &
                                               qs_ot_env(ispin)%matrix_dx_im)
                  ELSE
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx, &
                                               qs_ot_env(ispin)%matrix_gx_im, &
                                               qs_ot_env(ispin)%matrix_dx, &
                                               qs_ot_env(ispin)%matrix_dx_im)
                  END IF
                  CPASSERT(qs_ot_env(ispin)%kpoint_weight > 0.0_dp)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, &
                                   qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, &
                                   qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
               ELSE
                  IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                               qs_ot_env(ispin)%matrix_dx)
                  ELSE
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx)
                  END IF
               END IF
               channel_gnorm = 0.0_dp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
               channel_gnorm = channel_gnorm + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
                  channel_gnorm = channel_gnorm + tmp
               END IF
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner .AND. &
                   channel_gnorm <= 0.0_dp) THEN
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx, &
                                               qs_ot_env(ispin)%matrix_gx_im, &
                                               qs_ot_env(ispin)%matrix_dx, &
                                               qs_ot_env(ispin)%matrix_dx_im)
                     CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, &
                                      qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
                     CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, &
                                      qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
                  ELSE
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx, &
                                               qs_ot_env(ispin)%matrix_dx)
                  END IF
                  channel_gnorm = 0.0_dp
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
                  channel_gnorm = channel_gnorm + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                    qs_ot_env(ispin)%matrix_dx_im, tmp)
                     channel_gnorm = channel_gnorm + tmp
                  END IF
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + channel_gnorm
            END DO
            IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
               logger => cp_get_default_logger()
               WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
            END IF
            DO ispin = 1, nspin
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, -1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, -1.0_dp)
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  ! right now no preconditioner yet
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  END IF
               END DO
               DO ispin = 1, nspin
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, -1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, -1.0_dp)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_preconditioned_gx
               ELSE
                  qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_gx
               END IF
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_dx, qs_ot_env(ispin)%ener_gx)
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner .AND. tmp <= 0.0_dp) THEN
                  qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_gx
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               qs_ot_env(ispin)%ener_dx = -qs_ot_env(ispin)%ener_dx
            END DO
         END IF
      ELSE
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               channel_gnorm = 0.0_dp
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                                 qs_ot_env(ispin)%matrix_preconditioned_gx, tmp)
               ELSE
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
               END IF
               channel_gnorm = channel_gnorm + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                    qs_ot_env(ispin)%matrix_preconditioned_gx_im, tmp)
                  ELSE
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                  END IF
                  channel_gnorm = channel_gnorm + tmp
               END IF
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  IF (channel_gnorm <= 0.0_dp) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                     qs_ot_env(ispin)%matrix_gx)
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_copy(qs_ot_env(ispin)%matrix_preconditioned_gx_im, &
                                        qs_ot_env(ispin)%matrix_gx_im)
                     END IF
                     channel_gnorm = 0.0_dp
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                                    qs_ot_env(ispin)%matrix_preconditioned_gx, tmp)
                     channel_gnorm = channel_gnorm + tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                       qs_ot_env(ispin)%matrix_preconditioned_gx_im, tmp)
                        channel_gnorm = channel_gnorm + tmp
                     END IF
                  END IF
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx, &
                                  qs_ot_env(ispin)%matrix_preconditioned_gx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_im, &
                                     qs_ot_env(ispin)%matrix_preconditioned_gx_im)
                  END IF
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + channel_gnorm
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                    qs_ot_env(ispin)%ener_preconditioned_gx)
                  IF (tmp <= 0.0_dp) THEN
                     qs_ot_env(ispin)%ener_preconditioned_gx = qs_ot_env(ispin)%ener_gx
                     tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
                  END IF
                  qs_ot_env(ispin)%ener_gx = qs_ot_env(ispin)%ener_preconditioned_gx
               ELSE
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
            END DO
         END IF
      END IF

      k = 0
      n = 0
      nener = 0
      IF (do_ks) THEN
         CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
         DO ispin = 1, nspin
            CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
            k = k + itmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) k = k + itmp
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            nener = nener + SIZE(qs_ot_env(ispin)%ener_x)
         END DO
      END IF
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      nvariables = REAL(INT(n, KIND=int_8)*INT(k, KIND=int_8) + nener, KIND=dp)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvariables)
      ! Handling the case of no free variables to optimize
      IF (nvariables > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/nvariables)
         qs_ot_env(1)%gradient = -qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      CALL timestop(handle)

   END SUBROUTINE ot_new_sd_direction

!
! creates a new CG direction. Implements Polak-Ribierre variant
! using the preconditioner if associated
! also updates the gradient for line search
!
! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param para_env_inter_kp communicator between distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_new_cg_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_cg_direction'

      INTEGER                                            :: handle, ispin, itmp, k, n, nener, nspin
      LOGICAL                                            :: do_ener, do_ks
      REAL(KIND=dp)                                      :: beta_pr, gnorm_cross, nvariables, &
                                                            test_down, tmp
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      nspin = SIZE(qs_ot_env)
      logger => cp_get_default_logger()

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      gnorm_cross = 0.0_dp
      IF (do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
            gnorm_cross = gnorm_cross + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
               gnorm_cross = gnorm_cross + tmp
            END IF
         END DO
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
               ! added 0.5, because we have (antisymmetry) only half the number of variables
               gnorm_cross = gnorm_cross + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                 qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                  gnorm_cross = gnorm_cross + 0.5_dp*tmp
               END IF
            END DO
         END IF
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx_old)
            gnorm_cross = gnorm_cross + tmp
         END DO
      END IF
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm_cross)

      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN

         DO ispin = 1, nspin
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                            qs_ot_env(ispin)%matrix_preconditioned_gx_im, &
                                            qs_ot_env(ispin)%matrix_gx_old, &
                                            qs_ot_env(ispin)%matrix_gx_old_im)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, &
                                            qs_ot_env(ispin)%matrix_gx_im, &
                                            qs_ot_env(ispin)%matrix_gx_old, &
                                            qs_ot_env(ispin)%matrix_gx_old_im)
               END IF
               CPASSERT(qs_ot_env(ispin)%kpoint_weight > 0.0_dp)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, &
                                qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, &
                                qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
            ELSE
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                            qs_ot_env(ispin)%matrix_gx_old)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old)
               END IF
            END IF
         END DO
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               END IF
            END DO
            IF (.NOT. qs_ot_env(1)%settings%occupation_preconditioner) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_im, &
                                     qs_ot_env(ispin)%matrix_gx_old_im)
                  END IF
               END DO
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  ! right now no preconditioner yet
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  END IF
               END DO
               DO ispin = 1, nspin
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_old_im)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_preconditioned_gx
               ELSE
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx
               END IF
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx_old)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
            END DO
         END IF
      ELSE
         IF (do_ks) THEN
            qs_ot_env(1)%gnorm = 0.0_dp
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                                 qs_ot_env(ispin)%matrix_preconditioned_gx, tmp)
               ELSE
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, &
                                  qs_ot_env(ispin)%matrix_preconditioned_gx)
               ELSE
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
               END IF
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                    qs_ot_env(ispin)%matrix_preconditioned_gx_im, tmp)
                  ELSE
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_gx_im, tmp)
                  END IF
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
                  IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, &
                                     qs_ot_env(ispin)%matrix_preconditioned_gx_im)
                  ELSE
                     CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
                  END IF
               END IF
            END DO
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
                  ! added 0.5, because we have (antisymmetry) only half the number of variables
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                  END IF
               END DO
            END IF
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%settings%occupation_preconditioner) THEN
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                    qs_ot_env(ispin)%ener_preconditioned_gx)
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_preconditioned_gx
               ELSE
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx
               END IF
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
            END DO
         END IF
      END IF

      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      IF (ASSOCIATED(qs_ot_env(1)%preconditioner) .AND. &
          .NOT. qs_ot_env(1)%settings%occupation_preconditioner .AND. &
          qs_ot_env(1)%gnorm < 0.0_dp) THEN
         WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
      END IF

      ! Occupation-unweighted orbital and energy blocks are useful preconditioner inputs, but only
      ! their coupled direction is relevant for Mermin descent. Fall back as a whole when that
      ! direction is not descending; rejecting individual blocks would destroy the Schur coupling.
      IF (qs_ot_env(1)%settings%occupation_preconditioner .AND. qs_ot_env(1)%gnorm <= 0.0_dp) THEN
         qs_ot_env(1)%gnorm = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx, &
                                               qs_ot_env(ispin)%matrix_gx_im, &
                                               qs_ot_env(ispin)%matrix_gx_old, &
                                               qs_ot_env(ispin)%matrix_gx_old_im)
                     CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, &
                                      qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
                     CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, &
                                      qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight))
                  ELSE
                     CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                               qs_ot_env(ispin)%matrix_gx, &
                                               qs_ot_env(ispin)%matrix_gx_old)
                  END IF
               ELSE
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, &
                                     qs_ot_env(ispin)%matrix_gx_im)
                  END IF
               END IF
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                 qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + 0.5_dp*tmp
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx
               qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + &
                                    DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                                qs_ot_env(ispin)%ener_gx_old)
            END DO
         END IF
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%gnorm)
      END IF

      k = 0
      n = 0
      nener = 0
      IF (do_ks) THEN
         CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
         DO ispin = 1, nspin
            CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
            k = k + itmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) k = k + itmp
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            nener = nener + SIZE(qs_ot_env(ispin)%ener_x)
         END DO
      END IF
      nvariables = REAL(INT(n, KIND=int_8)*INT(k, KIND=int_8) + nener, KIND=dp)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvariables)
      ! Handling the case of no free variables to optimize
      IF (nvariables > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/nvariables)
         beta_pr = (qs_ot_env(1)%gnorm - gnorm_cross)/qs_ot_env(1)%gnorm_old
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         beta_pr = 0.0_dp
      END IF
      beta_pr = MAX(beta_pr, 0.0_dp) ! reset to SD

      test_down = 0.0_dp
      IF (do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old, &
                           alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
            test_down = test_down + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_old_im, &
                              alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
               test_down = test_down + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx, &
                              alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
               test_down = test_down + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_gx_im, &
                                 alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                 qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  test_down = test_down + 0.5_dp*tmp
               END IF
            END IF
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ener_dx = beta_pr*qs_ot_env(ispin)%ener_dx - &
                                       qs_ot_env(ispin)%ener_gx_old
            tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_dx)
            test_down = test_down + tmp
         END DO
      END IF
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, test_down)

      IF (test_down >= 0.0_dp) THEN ! reset to SD
         beta_pr = 0.0_dp
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old, &
                              alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_old_im, &
                                 alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                                 qs_ot_env(ispin)%rot_mat_gx, alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                                    qs_ot_env(ispin)%rot_mat_gx_im, &
                                    alpha_scalar=beta_pr, beta_scalar=-1.0_dp)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_dx = beta_pr*qs_ot_env(ispin)%ener_dx - &
                                          qs_ot_env(ispin)%ener_gx_old
            END DO
         END IF
      END IF
      ! since we change the direction we have to adjust the gradient
      qs_ot_env(1)%gradient = beta_pr*qs_ot_env(1)%gradient - qs_ot_env(1)%gnorm
      qs_ot_env(1)%gnorm_old = qs_ot_env(1)%gnorm

      CALL timestop(handle)

   END SUBROUTINE ot_new_cg_direction

! **************************************************************************************************
!> \brief Returns the smallest shift y <- y + shift*s that meets relative L-BFGS curvature.
!> \param sy scalar product s.y
!> \param ss scalar product s.s
!> \param yy scalar product y.y
!> \param curvature_tol requested lower bound for (s.y)/sqrt((s.s)(y.y))
!> \return non-negative finite damping shift, or zero when no valid shift is needed or available
! **************************************************************************************************
   PURE FUNCTION lbfgs_curvature_damping_shift(sy, ss, yy, curvature_tol) RESULT(damping_shift)
      REAL(KIND=dp), INTENT(IN)                          :: sy, ss, yy, curvature_tol
      REAL(KIND=dp)                                      :: damping_shift

      REAL(KIND=dp)                                      :: curvature_discriminant, target_sy

      damping_shift = 0.0_dp
      IF (.NOT. ieee_is_finite(sy) .OR. .NOT. ieee_is_finite(ss) .OR. &
          .NOT. ieee_is_finite(yy) .OR. ss <= 0.0_dp .OR. yy <= 0.0_dp .OR. &
          curvature_tol < 0.0_dp .OR. curvature_tol >= 1.0_dp) RETURN

      curvature_discriminant = MAX(0.0_dp, ss*yy - sy*sy)
      target_sy = curvature_tol/SQRT(MAX(TINY(1.0_dp), 1.0_dp - curvature_tol**2))* &
                  SQRT(curvature_discriminant)
      IF (sy <= target_sy) THEN
         target_sy = target_sy*(1.0_dp + SQRT(EPSILON(1.0_dp)))
         damping_shift = MAX(0.0_dp, (target_sy - sy)/ss)
         IF (.NOT. ieee_is_finite(damping_shift)) damping_shift = 0.0_dp
      END IF
   END FUNCTION lbfgs_curvature_damping_shift

! **************************************************************************************************
!> \brief Decides whether an L-BFGS history must be discarded after excessive gradient growth.
!> \param current_gradient_norm_sq squared norm at the current accepted point
!> \param previous_gradient_norm_sq squared norm at the preceding accepted point
!> \return true when both norms are finite and the gradient grew by more than a factor of ten
! **************************************************************************************************
   PURE FUNCTION lbfgs_history_restart_required(current_gradient_norm_sq, &
                                                previous_gradient_norm_sq) RESULT(restart_history)
      REAL(KIND=dp), INTENT(IN)                          :: current_gradient_norm_sq, &
                                                            previous_gradient_norm_sq
      LOGICAL                                            :: restart_history

      REAL(KIND=dp), PARAMETER :: gradient_restart_factor = 10.0_dp

      restart_history = ieee_is_finite(current_gradient_norm_sq) .AND. &
                        ieee_is_finite(previous_gradient_norm_sq) .AND. &
                        previous_gradient_norm_sq > TINY(1.0_dp) .AND. &
                        current_gradient_norm_sq > gradient_restart_factor**2* &
                        previous_gradient_norm_sq
   END FUNCTION lbfgs_history_restart_required

! **************************************************************************************************
!> \brief Decides whether L-BFGS must recover from a collapsed accepted line-search step.
!> \param accepted_step step selected by the preceding line search
!> \param reference_step configured initial line-search step
!> \return true for a finite step smaller than one millionth of a valid reference step
! **************************************************************************************************
   PURE FUNCTION lbfgs_step_restart_required(accepted_step, reference_step) RESULT(restart_step)
      REAL(KIND=dp), INTENT(IN)                          :: accepted_step, reference_step
      LOGICAL                                            :: restart_step

      REAL(KIND=dp), PARAMETER :: step_restart_factor = 1.0E-6_dp

      restart_step = ieee_is_finite(accepted_step) .AND. &
                     ieee_is_finite(reference_step) .AND. &
                     reference_step > 0.0_dp .AND. &
                     ABS(accepted_step) < step_restart_factor*reference_step
   END FUNCTION lbfgs_step_restart_required

! **************************************************************************************************
!> \brief Normalize and damp an occupation-response secant against the conventional H0 direction.
!>
!> The response direction is first normalized to the product norm of H0*g. If its directional
!> derivative is smaller than the Powell bound, it is mixed with H0*g until
!> g^T*p >= 0.2*g^T*H0*g. Equal-norm normalization and convex mixing bound the target direction,
!> while the curvature bound keeps the inverse-BFGS completion positive and well conditioned.
!> \param g_dot_response scalar product of g and the raw response direction
!> \param response_norm_sq squared product norm of the raw response direction
!> \param g_dot_h0_g scalar product of g and H0*g
!> \param h0_g_norm_sq squared product norm of H0*g
!> \param response_scale scale applied to the raw response direction
!> \param response_weight convex weight of the normalized response direction
!> \param valid true if finite positive input permits a response secant
! **************************************************************************************************
   PURE SUBROUTINE lbfgs_response_secant_parameters( &
      g_dot_response, response_norm_sq, g_dot_h0_g, h0_g_norm_sq, &
      response_scale, response_weight, valid)
      REAL(KIND=dp), INTENT(IN)                          :: g_dot_response, response_norm_sq, &
                                                            g_dot_h0_g, h0_g_norm_sq
      REAL(KIND=dp), INTENT(OUT)                         :: response_scale, response_weight
      LOGICAL, INTENT(OUT)                               :: valid

      REAL(KIND=dp), PARAMETER                           :: powell_fraction = 0.2_dp

      REAL(KIND=dp)                                      :: scaled_curvature

      response_scale = 1.0_dp
      response_weight = 0.0_dp
      valid = ieee_is_finite(g_dot_response) .AND. &
              ieee_is_finite(response_norm_sq) .AND. &
              ieee_is_finite(g_dot_h0_g) .AND. &
              ieee_is_finite(h0_g_norm_sq) .AND. &
              g_dot_response > 0.0_dp .AND. &
              response_norm_sq > TINY(1.0_dp) .AND. &
              g_dot_h0_g > 0.0_dp .AND. &
              h0_g_norm_sq > TINY(1.0_dp)
      IF (.NOT. valid) RETURN

      response_scale = SQRT(h0_g_norm_sq/response_norm_sq)
      scaled_curvature = response_scale*g_dot_response
      IF (.NOT. ieee_is_finite(response_scale) .OR. &
          .NOT. ieee_is_finite(scaled_curvature) .OR. scaled_curvature <= 0.0_dp) THEN
         response_scale = 1.0_dp
         valid = .FALSE.
         RETURN
      END IF

      response_weight = 1.0_dp
      IF (scaled_curvature < powell_fraction*g_dot_h0_g) THEN
         response_weight = (1.0_dp - powell_fraction)*g_dot_h0_g/ &
                           (g_dot_h0_g - scaled_curvature)
      END IF
   END SUBROUTINE lbfgs_response_secant_parameters

! **************************************************************************************************
!> \brief Apply the L-BFGS initial inverse Hessian in the OT product space.
!>
!> The conventional OT preconditioner defines the base operator H0. If occupation response is
!> enabled, a norm-bounded Powell-damped response direction p defines a positive inverse-BFGS
!> secant completion. It maps the current physical gradient g to p while remaining a linear
!> operator for the arbitrary vector q produced by the first L-BFGS loop:
!>
!> H = (I - rho*p*g^T)*H0*(I - rho*g*p^T) + rho*p*p^T,
!>
!> where rho = 1/(g^T*p). The m+1 history slot is scratch space for p; callers restore the current
!> physical gradient there after this routine returns.
!> \param qs_ot_env OT environments for all local spin/k-point channels
!> \param gamma scalar H0 scale for orbital and energy variables without an OT preconditioner
!> \param gamma_rotation H0 scale for rotation variables
!> \param scratch_index history slot used to store p
!> \param para_env_inter_kp communicator between distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_apply_lbfgs_initial_inverse( &
      qs_ot_env, gamma, gamma_rotation, scratch_index, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      REAL(KIND=dp), INTENT(IN)                          :: gamma, gamma_rotation
      INTEGER, INTENT(IN)                                :: scratch_index
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      INTEGER                                            :: ispin
      LOGICAL                                            :: valid_response
      REAL(KIND=dp) :: coeff_p, coeff_v, g_dot_h0_g, g_dot_h0_q, g_dot_p, h0_g_norm_sq, &
         kpoint_scale, p_dot_q, p_norm_sq, response_scale, response_weight, rho, tmp

      ! Apply the base operator to q, which is stored in the *_gx_old work arrays.
      DO ispin = 1, SIZE(qs_ot_env)
         IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx_old, &
                                         qs_ot_env(ispin)%matrix_gx_old_im, &
                                         qs_ot_env(ispin)%matrix_dx, &
                                         qs_ot_env(ispin)%matrix_dx_im)
               CPASSERT(qs_ot_env(ispin)%kpoint_weight > 0.0_dp)
               kpoint_scale = qs_ot_kpoint_preconditioner_scale( &
                              qs_ot_env(ispin)%kpoint_weight)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, kpoint_scale)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, kpoint_scale)
            ELSE
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx_old, &
                                         qs_ot_env(ispin)%matrix_dx)
            END IF
         ELSE
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old)
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, gamma)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx_im, &
                               qs_ot_env(ispin)%matrix_gx_old_im)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, gamma)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx_old)
            CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, gamma_rotation)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, &
                               qs_ot_env(ispin)%rot_mat_gx_old_im)
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, gamma_rotation)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            qs_ot_env(ispin)%ener_dx = gamma*qs_ot_env(ispin)%ener_gx_old
         END IF
      END DO
      IF (.NOT. qs_ot_env(1)%settings%occupation_preconditioner) RETURN

      ! Build the desired occupation-preconditioned product direction p in the scratch slot.
      DO ispin = 1, SIZE(qs_ot_env)
         CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx))
         IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx_im))
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                         qs_ot_env(ispin)%matrix_preconditioned_gx_im, &
                                         qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                                         qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix)
               kpoint_scale = qs_ot_kpoint_preconditioner_scale( &
                              qs_ot_env(ispin)%kpoint_weight)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                                kpoint_scale)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                                kpoint_scale)
            ELSE
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_preconditioned_gx, &
                                         qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix)
            END IF
         ELSE
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                            qs_ot_env(ispin)%matrix_preconditioned_gx)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx_im))
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                               qs_ot_env(ispin)%matrix_preconditioned_gx_im)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                            qs_ot_env(ispin)%rot_mat_gx)
            CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                             gamma_rotation)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                               qs_ot_env(ispin)%rot_mat_gx_im)
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                                gamma_rotation)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            CPASSERT(ASSOCIATED(qs_ot_env(ispin)%ener_preconditioned_gx))
            qs_ot_env(ispin)%ener_h_e(scratch_index, :) = &
               qs_ot_env(ispin)%ener_preconditioned_gx
         END IF
      END DO

      g_dot_h0_q = 0.0_dp
      g_dot_p = 0.0_dp
      p_dot_q = 0.0_dp
      p_norm_sq = 0.0_dp
      DO ispin = 1, SIZE(qs_ot_env)
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
         g_dot_h0_q = g_dot_h0_q + tmp
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                        qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, tmp)
         g_dot_p = g_dot_p + tmp
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                        qs_ot_env(ispin)%matrix_gx_old, tmp)
         p_dot_q = p_dot_q + tmp
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                        qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, tmp)
         p_norm_sq = p_norm_sq + tmp
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                           qs_ot_env(ispin)%matrix_dx_im, tmp)
            g_dot_h0_q = g_dot_h0_q + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                           qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, tmp)
            g_dot_p = g_dot_p + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                           qs_ot_env(ispin)%matrix_gx_old_im, tmp)
            p_dot_q = p_dot_q + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                           qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, tmp)
            p_norm_sq = p_norm_sq + tmp
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                           qs_ot_env(ispin)%rot_mat_dx, tmp)
            g_dot_h0_q = g_dot_h0_q + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                           qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, tmp)
            g_dot_p = g_dot_p + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                           qs_ot_env(ispin)%rot_mat_gx_old, tmp)
            p_dot_q = p_dot_q + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                           qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, tmp)
            p_norm_sq = p_norm_sq + 0.5_dp*tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                              qs_ot_env(ispin)%rot_mat_dx_im, tmp)
               g_dot_h0_q = g_dot_h0_q + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                              qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, tmp)
               g_dot_p = g_dot_p + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                              qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
               p_dot_q = p_dot_q + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, tmp)
               p_norm_sq = p_norm_sq + 0.5_dp*tmp
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            g_dot_h0_q = g_dot_h0_q + &
                         DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_dx)
            g_dot_p = g_dot_p + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                            qs_ot_env(ispin)%ener_h_e(scratch_index, :))
            p_dot_q = p_dot_q + DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(scratch_index, :), &
                                            qs_ot_env(ispin)%ener_gx_old)
            p_norm_sq = p_norm_sq + DOT_PRODUCT( &
                        qs_ot_env(ispin)%ener_h_e(scratch_index, :), &
                        qs_ot_env(ispin)%ener_h_e(scratch_index, :))
         END IF
      END DO
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, g_dot_h0_q)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, g_dot_p)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, p_dot_q)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, p_norm_sq)

      ! Overwrite q with H0*g; q is no longer needed after p^T*q has been formed.
      g_dot_h0_g = 0.0_dp
      h0_g_norm_sq = 0.0_dp
      DO ispin = 1, SIZE(qs_ot_env)
         IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx, &
                                         qs_ot_env(ispin)%matrix_gx_im, &
                                         qs_ot_env(ispin)%matrix_gx_old, &
                                         qs_ot_env(ispin)%matrix_gx_old_im)
               kpoint_scale = qs_ot_kpoint_preconditioner_scale( &
                              qs_ot_env(ispin)%kpoint_weight)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, kpoint_scale)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, kpoint_scale)
            ELSE
               CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                         qs_ot_env(ispin)%matrix_gx, &
                                         qs_ot_env(ispin)%matrix_gx_old)
            END IF
         ELSE
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old, gamma)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_gx_old_im, gamma)
            END IF
         END IF
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx_old, tmp)
         g_dot_h0_g = g_dot_h0_g + tmp
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old, &
                        qs_ot_env(ispin)%matrix_gx_old, tmp)
         h0_g_norm_sq = h0_g_norm_sq + tmp
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                           qs_ot_env(ispin)%matrix_gx_old_im, tmp)
            g_dot_h0_g = g_dot_h0_g + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old_im, &
                           qs_ot_env(ispin)%matrix_gx_old_im, tmp)
            h0_g_norm_sq = h0_g_norm_sq + tmp
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
            CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_gx_old, gamma_rotation)
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                           qs_ot_env(ispin)%rot_mat_gx_old, tmp)
            g_dot_h0_g = g_dot_h0_g + 0.5_dp*tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_old, &
                           qs_ot_env(ispin)%rot_mat_gx_old, tmp)
            h0_g_norm_sq = h0_g_norm_sq + 0.5_dp*tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                               qs_ot_env(ispin)%rot_mat_gx_im)
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_gx_old_im, gamma_rotation)
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                              qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
               g_dot_h0_g = g_dot_h0_g + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                              qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
               h0_g_norm_sq = h0_g_norm_sq + 0.5_dp*tmp
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            qs_ot_env(ispin)%ener_gx_old = gamma*qs_ot_env(ispin)%ener_gx
            g_dot_h0_g = g_dot_h0_g + &
                         DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                     qs_ot_env(ispin)%ener_gx_old)
            h0_g_norm_sq = h0_g_norm_sq + &
                           DOT_PRODUCT(qs_ot_env(ispin)%ener_gx_old, &
                                       qs_ot_env(ispin)%ener_gx_old)
         END IF
      END DO
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, g_dot_h0_g)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, h0_g_norm_sq)

      CALL lbfgs_response_secant_parameters( &
         g_dot_p, p_norm_sq, g_dot_h0_g, h0_g_norm_sq, &
         response_scale, response_weight, valid_response)
      IF (.NOT. valid_response) RETURN

      ! Replace p by the bounded, damped response target and update its scalar products.
      DO ispin = 1, SIZE(qs_ot_env)
         CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                          response_weight*response_scale)
         CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                        qs_ot_env(ispin)%matrix_gx_old, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp - response_weight)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                             response_weight*response_scale)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                           qs_ot_env(ispin)%matrix_gx_old_im, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp - response_weight)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                             response_weight*response_scale)
            CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                           qs_ot_env(ispin)%rot_mat_gx_old, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp - response_weight)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                                response_weight*response_scale)
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                              qs_ot_env(ispin)%rot_mat_gx_old_im, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp - response_weight)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            qs_ot_env(ispin)%ener_h_e(scratch_index, :) = &
               response_weight*response_scale*qs_ot_env(ispin)%ener_h_e(scratch_index, :) + &
               (1.0_dp - response_weight)*qs_ot_env(ispin)%ener_gx_old
         END IF
      END DO
      g_dot_p = response_weight*response_scale*g_dot_p + &
                (1.0_dp - response_weight)*g_dot_h0_g
      p_dot_q = response_weight*response_scale*p_dot_q + &
                (1.0_dp - response_weight)*g_dot_h0_q

      rho = 1.0_dp/g_dot_p
      coeff_p = -rho*g_dot_h0_q + p_dot_q*(rho + rho*rho*g_dot_h0_g)
      coeff_v = -rho*p_dot_q
      DO ispin = 1, SIZE(qs_ot_env)
         CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, &
                        qs_ot_env(ispin)%matrix_h_e(scratch_index)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=coeff_p)
         CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old, &
                        alpha_scalar=1.0_dp, beta_scalar=coeff_v)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, &
                           qs_ot_env(ispin)%matrix_h_e_im(scratch_index)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=coeff_p)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, &
                           qs_ot_env(ispin)%matrix_gx_old_im, &
                           alpha_scalar=1.0_dp, beta_scalar=coeff_v)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                           qs_ot_env(ispin)%rot_mat_h_e(scratch_index)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=coeff_p)
            CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx_old, &
                           alpha_scalar=1.0_dp, beta_scalar=coeff_v)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                              qs_ot_env(ispin)%rot_mat_h_e_im(scratch_index)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=coeff_p)
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                              qs_ot_env(ispin)%rot_mat_gx_old_im, &
                              alpha_scalar=1.0_dp, beta_scalar=coeff_v)
            END IF
         END IF
         IF (qs_ot_env(1)%settings%do_ener) THEN
            qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_dx + &
                                       coeff_p*qs_ot_env(ispin)%ener_h_e(scratch_index, :) + &
                                       coeff_v*qs_ot_env(ispin)%ener_gx_old
         END IF
      END DO
   END SUBROUTINE ot_apply_lbfgs_initial_inverse

! **************************************************************************************************
!> \brief Builds an L-BFGS direction in the fixed OT product chart.
!>
!> The selected OT preconditioner is used as the initial inverse Hessian. History slot m+1 stores
!> the previous accepted point and gradient; slots 1:m store the circular (s,y) secant history.
!> Orbital rotations are included with the antisymmetric-matrix metric used by the other OT
!> minimizers. Weak curvature can be regularized along the step, while rejected updates leave
!> older valid history intact. Non-descent directions and collapsed line-search steps reset the
!> history; the latter also restore the configured initial step.
!> \param qs_ot_env OT environments for all local spin/k-point channels
!> \param para_env_inter_kp communicator between distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_new_lbfgs_direction(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER :: routineN = 'ot_new_lbfgs_direction'
      REAL(KIND=dp), PARAMETER                           :: rotation_scale_max = 1.0E3_dp, &
                                                            rotation_scale_min = 1.0E-3_dp

      INTEGER                                            :: handle, i, ispin, itmp, j, k, m_history, &
                                                            n, nener, newest, nhistory, nspin
      INTEGER(KIND=int_8)                                :: nrotation, nvariables
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: history_index
      LOGICAL                                            :: do_ener, restart_history, restart_step
      REAL(KIND=dp) :: beta, current_gradient_norm_sq, curvature_tol, damping_shift, g_dot_z, &
         gamma, gamma_rotation, nvariables_global, previous_gradient_norm_sq, ss, ss_rotation, &
         stq, sy, sy_rotation, tmp, yy, yy_newest, yy_rotation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha, rho

      CALL timeset(routineN, handle)

      IF (.NOT. qs_ot_env(1)%settings%ks) THEN
         CPABORT("MINIMIZER LBFGS currently requires OT orbital variables")
      END IF
      nspin = SIZE(qs_ot_env)
      do_ener = qs_ot_env(1)%settings%do_ener
      m_history = qs_ot_env(1)%settings%diis_m
      curvature_tol = qs_ot_env(1)%settings%lbfgs_curvature_tol
      restart_history = .FALSE.
      restart_step = .FALSE.

      ! A local line search can accept a nearly isoenergetic point whose raw gradient is much
      ! larger than at the preceding point. CG naturally resets after such a loss of conjugacy,
      ! whereas retaining an L-BFGS model can amplify it. Discard only the secant model and use
      ! the physical initial inverse Hessian for the current direction. Likewise, restore the
      ! configured initial step if repeated interpolation has collapsed the accepted step by more
      ! than six orders of magnitude without satisfying the SCF convergence criterion.
      IF (qs_ot_env(1)%diis_iter > 0) THEN
         current_gradient_norm_sq = 0.0_dp
         previous_gradient_norm_sq = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
            current_gradient_norm_sq = current_gradient_norm_sq + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(m_history + 1)%matrix, &
                           qs_ot_env(ispin)%matrix_h_e(m_history + 1)%matrix, tmp)
            previous_gradient_norm_sq = previous_gradient_norm_sq + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                              qs_ot_env(ispin)%matrix_gx_im, tmp)
               current_gradient_norm_sq = current_gradient_norm_sq + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(m_history + 1)%matrix, &
                              qs_ot_env(ispin)%matrix_h_e_im(m_history + 1)%matrix, tmp)
               previous_gradient_norm_sq = previous_gradient_norm_sq + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_gx, tmp)
               current_gradient_norm_sq = current_gradient_norm_sq + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(m_history + 1)%matrix, &
                              qs_ot_env(ispin)%rot_mat_h_e(m_history + 1)%matrix, tmp)
               previous_gradient_norm_sq = previous_gradient_norm_sq + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                 qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                  current_gradient_norm_sq = current_gradient_norm_sq + 0.5_dp*tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(m_history + 1)%matrix, &
                                 qs_ot_env(ispin)%rot_mat_h_e_im(m_history + 1)%matrix, tmp)
                  previous_gradient_norm_sq = previous_gradient_norm_sq + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) THEN
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, qs_ot_env(ispin)%ener_gx)
               current_gradient_norm_sq = current_gradient_norm_sq + tmp
               tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(m_history + 1, :), &
                                 qs_ot_env(ispin)%ener_h_e(m_history + 1, :))
               previous_gradient_norm_sq = previous_gradient_norm_sq + tmp
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, current_gradient_norm_sq)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, previous_gradient_norm_sq)
         restart_step = lbfgs_step_restart_required(qs_ot_env(1)%ds_min, &
                                                    qs_ot_env(1)%settings%ds_min)
         restart_history = lbfgs_history_restart_required(current_gradient_norm_sq, &
                                                          previous_gradient_norm_sq) .OR. &
                           restart_step
         IF (restart_history) THEN
            qs_ot_env(1)%diis_iter = 1
            IF (restart_step) qs_ot_env(1)%ds_min = qs_ot_env(1)%settings%ds_min
            qs_ot_env(1)%OT_METHOD_FULL = "OT L-RST"
         END IF
      END IF

      ! Form a candidate secant pair in the direction workspaces relative to the preceding
      ! accepted line-search point. Commit it to the circular history only after the curvature
      ! test, so a rejected candidate cannot overwrite the oldest valid pair in a full buffer.
      IF (qs_ot_env(1)%diis_iter > 0 .AND. .NOT. restart_history) THEN
         j = MOD(qs_ot_env(1)%diis_iter - 1, m_history) + 1
         ss = 0.0_dp
         sy = 0.0_dp
         sy_rotation = 0.0_dp
         ss_rotation = 0.0_dp
         yy = 0.0_dp
         yy_rotation = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_x)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, &
                           qs_ot_env(ispin)%matrix_h_x(m_history + 1)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old, &
                           qs_ot_env(ispin)%matrix_h_e(m_history + 1)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_dx, tmp)
            ss = ss + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx_old, tmp)
            sy = sy + tmp
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx_old, tmp)
            yy = yy + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_x_im)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, &
                              qs_ot_env(ispin)%matrix_h_x_im(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old_im, &
                              qs_ot_env(ispin)%matrix_h_e_im(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_dx_im, &
                              qs_ot_env(ispin)%matrix_dx_im, tmp)
               ss = ss + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_dx_im, &
                              qs_ot_env(ispin)%matrix_gx_old_im, tmp)
               sy = sy + tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_old_im, &
                              qs_ot_env(ispin)%matrix_gx_old_im, tmp)
               yy = yy + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_x)
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                              qs_ot_env(ispin)%rot_mat_h_x(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old, &
                              qs_ot_env(ispin)%rot_mat_h_e(m_history + 1)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_dx, tmp)
               ss = ss + 0.5_dp*tmp
               ss_rotation = ss_rotation + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
               sy = sy + 0.5_dp*tmp
               sy_rotation = sy_rotation + 0.5_dp*tmp
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx_old, tmp)
               yy = yy + 0.5_dp*tmp
               yy_rotation = yy_rotation + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, qs_ot_env(ispin)%rot_mat_x_im)
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                                 qs_ot_env(ispin)%rot_mat_h_x_im(m_history + 1)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                  qs_ot_env(ispin)%rot_mat_gx_im)
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                 qs_ot_env(ispin)%rot_mat_h_e_im(m_history + 1)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_dx_im, &
                                 qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  ss = ss + 0.5_dp*tmp
                  ss_rotation = ss_rotation + 0.5_dp*tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_dx_im, &
                                 qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                  sy = sy + 0.5_dp*tmp
                  sy_rotation = sy_rotation + 0.5_dp*tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                 qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                  yy = yy + 0.5_dp*tmp
                  yy_rotation = yy_rotation + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) THEN
               qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_x - &
                                          qs_ot_env(ispin)%ener_h_x(m_history + 1, :)
               qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx - &
                                              qs_ot_env(ispin)%ener_h_e(m_history + 1, :)
               ss = ss + DOT_PRODUCT(qs_ot_env(ispin)%ener_dx, qs_ot_env(ispin)%ener_dx)
               sy = sy + DOT_PRODUCT(qs_ot_env(ispin)%ener_dx, qs_ot_env(ispin)%ener_gx_old)
               yy = yy + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx_old, qs_ot_env(ispin)%ener_gx_old)
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, ss)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, sy)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, yy)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, ss_rotation)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, sy_rotation)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, yy_rotation)

         ! Shift y along s just enough to satisfy the requested relative curvature. With the
         ! product-chart metric this preserves orbital-rotation covariance. The target below is
         ! obtained by solving (s.y)^2 >= tol^2 (s.s)(y.y) for y <- y + shift*s.
         damping_shift = 0.0_dp
         IF (qs_ot_env(1)%settings%lbfgs_damping) THEN
            damping_shift = lbfgs_curvature_damping_shift(sy, ss, yy, curvature_tol)
         END IF
         IF (damping_shift > 0.0_dp) THEN
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old, &
                              qs_ot_env(ispin)%matrix_dx, &
                              alpha_scalar=1.0_dp, beta_scalar=damping_shift)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old_im, &
                                 qs_ot_env(ispin)%matrix_dx_im, &
                                 alpha_scalar=1.0_dp, beta_scalar=damping_shift)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old, &
                                 qs_ot_env(ispin)%rot_mat_dx, &
                                 alpha_scalar=1.0_dp, beta_scalar=damping_shift)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                    qs_ot_env(ispin)%rot_mat_dx_im, &
                                    alpha_scalar=1.0_dp, beta_scalar=damping_shift)
                  END IF
               END IF
               IF (do_ener) THEN
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx_old + &
                                                 damping_shift*qs_ot_env(ispin)%ener_dx
               END IF
            END DO
            yy = yy + 2.0_dp*damping_shift*sy + damping_shift**2*ss
            yy_rotation = yy_rotation + 2.0_dp*damping_shift*sy_rotation + &
                          damping_shift**2*ss_rotation
            sy = sy + damping_shift*ss
            sy_rotation = sy_rotation + damping_shift*ss_rotation
         END IF
         IF (ieee_is_finite(sy) .AND. ieee_is_finite(ss) .AND. ieee_is_finite(yy) .AND. &
             ss > 0.0_dp .AND. yy > 0.0_dp .AND. sy > curvature_tol*SQRT(ss*yy)) THEN
            DO ispin = 1, nspin
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                               qs_ot_env(ispin)%matrix_dx)
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                               qs_ot_env(ispin)%matrix_gx_old)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                                  qs_ot_env(ispin)%matrix_dx_im)
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                  qs_ot_env(ispin)%matrix_gx_old_im)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_dx)
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_gx_old)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_dx_im)
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                     qs_ot_env(ispin)%rot_mat_gx_old_im)
                  END IF
               END IF
               IF (do_ener) THEN
                  qs_ot_env(ispin)%ener_h_x(j, :) = qs_ot_env(ispin)%ener_dx
                  qs_ot_env(ispin)%ener_h_e(j, :) = qs_ot_env(ispin)%ener_gx_old
               END IF
            END DO
            qs_ot_env(1)%lbfgs_rho(j) = 1.0_dp/sy
            qs_ot_env(1)%lbfgs_yy(j) = yy
            qs_ot_env(1)%lbfgs_sy_rotation(j) = sy_rotation
            qs_ot_env(1)%lbfgs_yy_rotation(j) = yy_rotation
            qs_ot_env(1)%diis_iter = qs_ot_env(1)%diis_iter + 1
         ELSE
            ! Skip this pair. Older pairs remain valid in the fixed OT chart and are retained;
            ! the descent safeguard below still clears all history if their recursion fails.
            qs_ot_env(1)%OT_METHOD_FULL = "OT LSKIP"
         END IF
      ELSE IF (qs_ot_env(1)%diis_iter <= 0) THEN
         qs_ot_env(1)%diis_iter = 1
      END IF

      DO ispin = 1, nspin
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old, qs_ot_env(ispin)%matrix_gx)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_gx_old_im, qs_ot_env(ispin)%matrix_gx_im)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old, qs_ot_env(ispin)%rot_mat_gx)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                               qs_ot_env(ispin)%rot_mat_gx_im)
            END IF
         END IF
         IF (do_ener) qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx
      END DO

      nhistory = MIN(qs_ot_env(1)%diis_iter - 1, m_history)
      ALLOCATE (alpha(MAX(1, nhistory)), rho(MAX(1, nhistory)), &
                history_index(MAX(1, nhistory)))
      gamma_rotation = 1.0_dp
      yy_newest = 0.0_dp

      ! First loop, newest to oldest: q <- V^T g.
      IF (nhistory > 0) THEN
         newest = MOD(qs_ot_env(1)%diis_iter - 2, m_history) + 1
         yy_newest = qs_ot_env(1)%lbfgs_yy(newest)
         IF (qs_ot_env(1)%settings%do_rotation .AND. &
             ieee_is_finite(qs_ot_env(1)%lbfgs_sy_rotation(newest)) .AND. &
             ieee_is_finite(qs_ot_env(1)%lbfgs_yy_rotation(newest)) .AND. &
             qs_ot_env(1)%lbfgs_sy_rotation(newest) > 0.0_dp .AND. &
             qs_ot_env(1)%lbfgs_yy_rotation(newest) > 0.0_dp) THEN
            gamma_rotation = MAX(rotation_scale_min, MIN(rotation_scale_max, &
                                                         qs_ot_env(1)%lbfgs_sy_rotation(newest)/ &
                                                         qs_ot_env(1)%lbfgs_yy_rotation(newest)))
         END IF
         DO i = 1, nhistory
            history_index(i) = MODULO(newest - i, m_history) + 1
            stq = 0.0_dp
            DO ispin = 1, nspin
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x(history_index(i))%matrix, &
                              qs_ot_env(ispin)%matrix_gx_old, tmp)
               stq = stq + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x_im(history_index(i))%matrix, &
                                 qs_ot_env(ispin)%matrix_gx_old_im, tmp)
                  stq = stq + tmp
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x(history_index(i))%matrix, &
                                 qs_ot_env(ispin)%rot_mat_gx_old, tmp)
                  stq = stq + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x_im(history_index(i))%matrix, &
                                    qs_ot_env(ispin)%rot_mat_gx_old_im, tmp)
                     stq = stq + 0.5_dp*tmp
                  END IF
               END IF
               IF (do_ener) THEN
                  stq = stq + DOT_PRODUCT(qs_ot_env(ispin)%ener_h_x(history_index(i), :), &
                                          qs_ot_env(ispin)%ener_gx_old)
               END IF
            END DO
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, stq)
            rho(i) = qs_ot_env(1)%lbfgs_rho(history_index(i))
            alpha(i) = rho(i)*stq
            DO ispin = 1, nspin
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old, &
                              qs_ot_env(ispin)%matrix_h_e(history_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_gx_old_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(history_index(i))%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old, &
                                 qs_ot_env(ispin)%rot_mat_h_e(history_index(i))%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_gx_old_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(history_index(i))%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=-alpha(i))
                  END IF
               END IF
               IF (do_ener) THEN
                  qs_ot_env(ispin)%ener_gx_old = qs_ot_env(ispin)%ener_gx_old - &
                                                 alpha(i)*qs_ot_env(ispin)%ener_h_e(history_index(i), :)
               END IF
            END DO
         END DO
      END IF

      ! Apply the conventional initial inverse Hessian and, when requested, its positive linear
      ! occupation-response completion to the arbitrary vector from the first loop.
      gamma = 1.0_dp
      IF (nhistory > 0 .AND. yy_newest > 0.0_dp) gamma = 1.0_dp/(rho(1)*yy_newest)
      CALL ot_apply_lbfgs_initial_inverse( &
         qs_ot_env, gamma, gamma_rotation, m_history + 1, para_env_inter_kp)

      ! The m+1 slot was scratch space for the occupation secant completion. Preserve the current
      ! accepted point and physical gradient there for the next L-BFGS update.
      DO ispin = 1, nspin
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(m_history + 1)%matrix, &
                         qs_ot_env(ispin)%matrix_x)
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(m_history + 1)%matrix, &
                         qs_ot_env(ispin)%matrix_gx)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%matrix_x_im)
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%matrix_gx_im)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%rot_mat_x)
            CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(m_history + 1)%matrix, &
                            qs_ot_env(ispin)%rot_mat_gx)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(m_history + 1)%matrix, &
                               qs_ot_env(ispin)%rot_mat_x_im)
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(m_history + 1)%matrix, &
                               qs_ot_env(ispin)%rot_mat_gx_im)
            END IF
         END IF
         IF (do_ener) THEN
            qs_ot_env(ispin)%ener_h_x(m_history + 1, :) = qs_ot_env(ispin)%ener_x
            qs_ot_env(ispin)%ener_h_e(m_history + 1, :) = qs_ot_env(ispin)%ener_gx
         END IF
      END DO

      ! Second loop, oldest to newest: z <- V H0 V^T g.
      DO i = nhistory, 1, -1
         beta = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(history_index(i))%matrix, &
                           qs_ot_env(ispin)%matrix_dx, tmp)
            beta = beta + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(history_index(i))%matrix, &
                              qs_ot_env(ispin)%matrix_dx_im, tmp)
               beta = beta + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(history_index(i))%matrix, &
                              qs_ot_env(ispin)%rot_mat_dx, tmp)
               beta = beta + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(history_index(i))%matrix, &
                                 qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  beta = beta + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) THEN
               beta = beta + DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(history_index(i), :), &
                                         qs_ot_env(ispin)%ener_dx)
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, beta)
         beta = rho(i)*beta
         DO ispin = 1, nspin
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx, &
                           qs_ot_env(ispin)%matrix_h_x(history_index(i))%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_dx_im, &
                              qs_ot_env(ispin)%matrix_h_x_im(history_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx, &
                              qs_ot_env(ispin)%rot_mat_h_x(history_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_dx_im, &
                                 qs_ot_env(ispin)%rot_mat_h_x_im(history_index(i))%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=alpha(i) - beta)
               END IF
            END IF
            IF (do_ener) THEN
               qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_dx + &
                                          (alpha(i) - beta)* &
                                          qs_ot_env(ispin)%ener_h_x(history_index(i), :)
            END IF
         END DO
      END DO

      g_dot_z = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
         g_dot_z = g_dot_z + tmp
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
            g_dot_z = g_dot_z + tmp
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
            g_dot_z = g_dot_z + 0.5_dp*tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                              qs_ot_env(ispin)%rot_mat_dx_im, tmp)
               g_dot_z = g_dot_z + 0.5_dp*tmp
            END IF
         END IF
         IF (do_ener) THEN
            g_dot_z = g_dot_z + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                            qs_ot_env(ispin)%ener_dx)
         END IF
      END DO
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, g_dot_z)

      ! Preserve a guaranteed descent direction even if the supplied initial preconditioner is
      ! indefinite or numerical noise corrupts the secant recursion.
      IF (.NOT. ieee_is_finite(g_dot_z) .OR. g_dot_z <= 0.0_dp) THEN
         qs_ot_env(1)%diis_iter = 1
         qs_ot_env(1)%OT_METHOD_FULL = "OT L-SD"
         IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
            DO ispin = 1, nspin
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, &
                                            qs_ot_env(ispin)%matrix_gx_im, &
                                            qs_ot_env(ispin)%matrix_dx, &
                                            qs_ot_env(ispin)%matrix_dx_im)
                  CPASSERT(qs_ot_env(ispin)%kpoint_weight > 0.0_dp)
                  tmp = qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, tmp)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, tmp)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                  END IF
               END IF
               IF (do_ener) qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_gx
            END DO
         ELSE
            DO ispin = 1, nspin
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_im)
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                  END IF
               END IF
               IF (do_ener) qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_gx
            END DO
         END IF
         g_dot_z = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
            g_dot_z = g_dot_z + tmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, qs_ot_env(ispin)%matrix_dx_im, tmp)
               g_dot_z = g_dot_z + tmp
            END IF
            IF (qs_ot_env(1)%settings%do_rotation) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
               g_dot_z = g_dot_z + 0.5_dp*tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                 qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                  g_dot_z = g_dot_z + 0.5_dp*tmp
               END IF
            END IF
            IF (do_ener) THEN
               g_dot_z = g_dot_z + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                               qs_ot_env(ispin)%ener_dx)
            END IF
         END DO
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, g_dot_z)
         IF (.NOT. ieee_is_finite(g_dot_z) .OR. g_dot_z <= 0.0_dp) THEN
            g_dot_z = 0.0_dp
            DO ispin = 1, nspin
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx, qs_ot_env(ispin)%matrix_gx)
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_dx, tmp)
               g_dot_z = g_dot_z + tmp
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_dx_im, qs_ot_env(ispin)%matrix_gx_im)
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                 qs_ot_env(ispin)%matrix_dx_im, tmp)
                  g_dot_z = g_dot_z + tmp
               END IF
               IF (qs_ot_env(1)%settings%do_rotation) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx, qs_ot_env(ispin)%rot_mat_gx)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, qs_ot_env(ispin)%rot_mat_dx, tmp)
                  g_dot_z = g_dot_z + 0.5_dp*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_dx_im, &
                                     qs_ot_env(ispin)%rot_mat_gx_im)
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                    qs_ot_env(ispin)%rot_mat_dx_im, tmp)
                     g_dot_z = g_dot_z + 0.5_dp*tmp
                  END IF
               END IF
               IF (do_ener) THEN
                  qs_ot_env(ispin)%ener_dx = qs_ot_env(ispin)%ener_gx
                  g_dot_z = g_dot_z + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                                  qs_ot_env(ispin)%ener_dx)
               END IF
            END DO
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, g_dot_z)
         END IF
      END IF

      DO ispin = 1, nspin
         CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx, -1.0_dp)
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_dx_im, -1.0_dp)
         END IF
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx, -1.0_dp)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_dx_im, -1.0_dp)
            END IF
         END IF
         IF (do_ener) qs_ot_env(ispin)%ener_dx = -qs_ot_env(ispin)%ener_dx
      END DO
      qs_ot_env(1)%gnorm = g_dot_z
      qs_ot_env(1)%gradient = -g_dot_z

      k = 0
      n = 0
      nener = 0
      nrotation = 0_int_8
      CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
      DO ispin = 1, nspin
         CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
         k = k + itmp
         IF (qs_ot_env(ispin)%has_complex_kpoint_state) k = k + itmp
         IF (qs_ot_env(1)%settings%do_rotation) THEN
            nrotation = nrotation + INT(itmp, KIND=int_8)*INT(itmp - 1, KIND=int_8)/2_int_8
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               nrotation = nrotation + INT(itmp, KIND=int_8)*INT(itmp + 1, KIND=int_8)/2_int_8
            END IF
         END IF
         IF (do_ener) nener = nener + SIZE(qs_ot_env(ispin)%ener_x)
      END DO
      nvariables = INT(n, KIND=int_8)*INT(k, KIND=int_8) + nrotation + &
                   INT(nener, KIND=int_8)
      nvariables_global = REAL(nvariables, KIND=dp)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvariables_global)
      IF (nvariables_global > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(g_dot_z)/nvariables_global)
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      DEALLOCATE (alpha, rho, history_index)
      CALL timestop(handle)

   END SUBROUTINE ot_new_lbfgs_direction

! **************************************************************************************************
!> \brief Build a preconditioned product residual for history-based OT minimizers.
!> \param qs_ot_env OT environments for all local spin/k-point channels
!> \param history_index destination slot in the residual history
!> \param step_scale scale applied after evaluating the descent product
!> \param para_env_inter_kp communicator between distributed k-point groups
!> \param gnorm product of the physical Mermin gradient and preconditioned residual
! **************************************************************************************************
   SUBROUTINE ot_build_history_residual(qs_ot_env, history_index, step_scale, &
                                        para_env_inter_kp, gnorm)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      INTEGER, INTENT(IN)                                :: history_index
      REAL(KIND=dp), INTENT(IN)                          :: step_scale
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp
      REAL(KIND=dp), INTENT(OUT)                         :: gnorm

      LOGICAL                                            :: use_occupation_response
      TYPE(cp_logger_type), POINTER                      :: logger

      use_occupation_response = qs_ot_env(1)%settings%occupation_preconditioner
      CALL ot_build_history_residual_once(qs_ot_env, history_index, step_scale, &
                                          use_occupation_response, gnorm)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm)

      ! The occupation response is a coupled orbital/energy approximation. Falling back one block
      ! at a time would change its Schur balance, so replace the complete product residual instead.
      IF (use_occupation_response .AND. gnorm <= 0.0_dp) THEN
         CALL ot_build_history_residual_once(qs_ot_env, history_index, step_scale, .FALSE., gnorm)
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, gnorm)
      END IF
      IF (gnorm < 0.0_dp) THEN
         logger => cp_get_default_logger()
         WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
      END IF
   END SUBROUTINE ot_build_history_residual

! **************************************************************************************************
!> \brief Build one physical or occupation-response product residual.
!> \param qs_ot_env OT environments for all local spin/k-point channels
!> \param history_index destination slot in the residual history
!> \param step_scale scale applied to the stored residual
!> \param use_occupation_response use the fixed-N occupation-response approximation
!> \param gnorm product of the physical Mermin gradient and unscaled residual
! **************************************************************************************************
   SUBROUTINE ot_build_history_residual_once(qs_ot_env, history_index, step_scale, &
                                             use_occupation_response, gnorm)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      INTEGER, INTENT(IN)                                :: history_index
      REAL(KIND=dp), INTENT(IN)                          :: step_scale
      LOGICAL, INTENT(IN)                                :: use_occupation_response
      REAL(KIND=dp), INTENT(OUT)                         :: gnorm

      INTEGER                                            :: ispin
      LOGICAL                                            :: do_ener, do_ks
      REAL(KIND=dp)                                      :: kpoint_scale, tmp
      TYPE(dbcsr_type), POINTER                          :: orbital_residual, orbital_residual_im

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      gnorm = 0.0_dp

      IF (do_ks) THEN
         DO ispin = 1, SIZE(qs_ot_env)
            IF (use_occupation_response) THEN
               CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx))
               orbital_residual => qs_ot_env(ispin)%matrix_preconditioned_gx
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CPASSERT(ASSOCIATED(qs_ot_env(ispin)%matrix_preconditioned_gx_im))
                  orbital_residual_im => qs_ot_env(ispin)%matrix_preconditioned_gx_im
               END IF
            ELSE
               orbital_residual => qs_ot_env(ispin)%matrix_gx
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  orbital_residual_im => qs_ot_env(ispin)%matrix_gx_im
               END IF
            END IF

            IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                            orbital_residual, orbital_residual_im, &
                                            qs_ot_env(ispin)%matrix_h_e(history_index)%matrix, &
                                            qs_ot_env(ispin)%matrix_h_e_im(history_index)%matrix)
                  CPASSERT(qs_ot_env(ispin)%kpoint_weight > 0.0_dp)
                  kpoint_scale = qs_ot_kpoint_preconditioner_scale(qs_ot_env(ispin)%kpoint_weight)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(history_index)%matrix, kpoint_scale)
                  CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(history_index)%matrix, kpoint_scale)
               ELSE
                  CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, orbital_residual, &
                                            qs_ot_env(ispin)%matrix_h_e(history_index)%matrix)
               END IF
            ELSE
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e(history_index)%matrix, orbital_residual)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_e_im(history_index)%matrix, &
                                  orbital_residual_im)
               END IF
            END IF

            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                           qs_ot_env(ispin)%matrix_h_e(history_index)%matrix, tmp)
            gnorm = gnorm + tmp
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(history_index)%matrix, step_scale)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                              qs_ot_env(ispin)%matrix_h_e_im(history_index)%matrix, tmp)
               gnorm = gnorm + tmp
               CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e_im(history_index)%matrix, step_scale)
            END IF

            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e(history_index)%matrix, &
                               qs_ot_env(ispin)%rot_mat_gx)
               CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                              qs_ot_env(ispin)%rot_mat_h_e(history_index)%matrix, tmp)
               gnorm = gnorm + 0.5_dp*tmp
               CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e(history_index)%matrix, step_scale)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_e_im(history_index)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_gx_im)
                  CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                 qs_ot_env(ispin)%rot_mat_h_e_im(history_index)%matrix, tmp)
                  gnorm = gnorm + 0.5_dp*tmp
                  CALL dbcsr_scale(qs_ot_env(ispin)%rot_mat_h_e_im(history_index)%matrix, step_scale)
               END IF
            END IF
         END DO
      END IF

      IF (do_ener) THEN
         DO ispin = 1, SIZE(qs_ot_env)
            IF (use_occupation_response) THEN
               CPASSERT(ASSOCIATED(qs_ot_env(ispin)%ener_preconditioned_gx))
               qs_ot_env(ispin)%ener_h_e(history_index, :) = &
                  qs_ot_env(ispin)%ener_preconditioned_gx
            ELSE
               qs_ot_env(ispin)%ener_h_e(history_index, :) = qs_ot_env(ispin)%ener_gx
            END IF
            gnorm = gnorm + DOT_PRODUCT(qs_ot_env(ispin)%ener_gx, &
                                        qs_ot_env(ispin)%ener_h_e(history_index, :))
            qs_ot_env(ispin)%ener_h_e(history_index, :) = &
               step_scale*qs_ot_env(ispin)%ener_h_e(history_index, :)
         END DO
      END IF
   END SUBROUTINE ot_build_history_residual_once

! **************************************************************************************************
!> \brief ...
!> \param qs_ot_env ...
!> \param para_env_inter_kp communicator between distributed k-point groups
! **************************************************************************************************
   SUBROUTINE ot_diis_step(qs_ot_env, para_env_inter_kp)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env_inter_kp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ot_diis_step'

      INTEGER                                            :: diis_bound, diis_m, handle, i, info, &
                                                            ispin, itmp, j, k, n, nener, nspin
      LOGICAL                                            :: do_ener, do_ks, do_ot_sd
      REAL(KIND=dp)                                      :: nvariables, overlap, tmp, tr_xnew_gx, &
                                                            tr_xold_gx
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener
      nspin = SIZE(qs_ot_env)

      diis_m = qs_ot_env(1)%settings%diis_m

      IF (qs_ot_env(1)%diis_iter < diis_m) THEN
         diis_bound = qs_ot_env(1)%diis_iter + 1
      ELSE
         diis_bound = diis_m
      END IF

      j = MOD(qs_ot_env(1)%diis_iter, diis_m) + 1 ! index in the circular array

      ! copy the position and the error vector in the diis buffers

      IF (do_ks) THEN
         DO ispin = 1, nspin
            CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, qs_ot_env(ispin)%matrix_x)
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                               qs_ot_env(ispin)%matrix_x_im)
            END IF
            IF (qs_ot_env(ispin)%settings%do_rotation) THEN
               CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, qs_ot_env(ispin)%rot_mat_x)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_copy(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                  qs_ot_env(ispin)%rot_mat_x_im)
               END IF
            END IF
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            qs_ot_env(ispin)%ener_h_x(j, :) = qs_ot_env(ispin)%ener_x(:)
         END DO
      END IF
      CALL ot_build_history_residual(qs_ot_env, j, -qs_ot_env(1)%ds_min, &
                                     para_env_inter_kp, qs_ot_env(1)%gnorm)
      k = 0
      n = 0
      nener = 0
      IF (do_ks) THEN
         CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
         DO ispin = 1, nspin
            CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
            k = k + itmp
            IF (qs_ot_env(ispin)%has_complex_kpoint_state) k = k + itmp
         END DO
      END IF
      IF (do_ener) THEN
         DO ispin = 1, nspin
            nener = nener + SIZE(qs_ot_env(ispin)%ener_x)
         END DO
      END IF
      nvariables = REAL(INT(n, KIND=int_8)*INT(k, KIND=int_8) + nener, KIND=dp)
      CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, nvariables)
      ! Handling the case of no free variables to optimize
      IF (nvariables > 0.0_dp) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/nvariables)
         qs_ot_env(1)%gradient = -qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      ! make the diis matrix and solve it
      DO i = 1, diis_bound
         ! I think there are two possible options, with and without preconditioner
         ! as a metric
         ! the second option seems most logical to me, and it seems marginally faster
         ! in some of the tests
         IF (.FALSE.) THEN
            qs_ot_env(1)%ls_diis(i, j) = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp)
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, tmp)
                        qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + 0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_e(j, :), qs_ot_env(ispin)%ener_h_e(i, :))
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) + tmp
               END DO
            END IF
         ELSE
            qs_ot_env(1)%ls_diis(i, j) = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 tmp)
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx_im, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    tmp)
                     qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_gx_im, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, tmp)
                        qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - &
                                                     qs_ot_env(1)%ds_min*0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_gx(:), qs_ot_env(ispin)%ener_h_e(i, :))
                  qs_ot_env(1)%ls_diis(i, j) = qs_ot_env(1)%ls_diis(i, j) - qs_ot_env(1)%ds_min*tmp
               END DO
            END IF
         END IF
         CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, qs_ot_env(1)%ls_diis(i, j))
         qs_ot_env(1)%ls_diis(j, i) = qs_ot_env(1)%ls_diis(i, j)
         qs_ot_env(1)%ls_diis(i, diis_bound + 1) = 1.0_dp
         qs_ot_env(1)%ls_diis(diis_bound + 1, i) = 1.0_dp
         qs_ot_env(1)%c_diis(i) = 0.0_dp
      END DO
      qs_ot_env(1)%ls_diis(diis_bound + 1, diis_bound + 1) = 0.0_dp
      qs_ot_env(1)%c_diis(diis_bound + 1) = 1.0_dp
      ! put in buffer, dgesv destroys
      qs_ot_env(1)%lss_diis = qs_ot_env(1)%ls_diis

      CALL DGESV(diis_bound + 1, 1, qs_ot_env(1)%lss_diis, diis_m + 1, qs_ot_env(1)%ipivot, &
                 qs_ot_env(1)%c_diis, diis_m + 1, info)

      IF (info /= 0) THEN
         do_ot_sd = .TRUE.
         WRITE (cp_logger_get_default_unit_nr(logger), *) "Singular DIIS matrix"
      ELSE
         do_ot_sd = .FALSE.
         IF (do_ks) THEN
            DO ispin = 1, nspin
               ! OK, add the vectors now
               CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
               DO i = 1, diis_bound
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_h_e(i)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
               END DO
               DO i = 1, diis_bound
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_h_x(i)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
               END DO
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%matrix_x_im, 0.0_dp)
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_h_e_im(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  END DO
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_h_x_im(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                  END DO
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x, 0.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x_im, 0.0_dp)
                  END IF
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_h_e(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_h_e_im(i)%matrix, &
                                       alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     END IF
                  END DO
                  DO i = 1, diis_bound
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_h_x(i)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_h_x_im(i)%matrix, &
                                       alpha_scalar=1.0_dp, beta_scalar=qs_ot_env(1)%c_diis(i))
                     END IF
                  END DO
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x(:) = 0.0_dp
               DO i = 1, diis_bound
                  qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) &
                                               + qs_ot_env(1)%c_diis(i)*qs_ot_env(ispin)%ener_h_e(i, :)
               END DO
               DO i = 1, diis_bound
                  qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) &
                                               + qs_ot_env(1)%c_diis(i)*qs_ot_env(ispin)%ener_h_x(i, :)
               END DO
            END DO
         END IF
         qs_ot_env(1)%diis_iter = qs_ot_env(1)%diis_iter + 1
         IF (qs_ot_env(1)%settings%safer_diis) THEN
            ! now, final check, is the step in fact in the direction of the -gradient ?
            ! if not we're walking towards a sadle point, and should avoid that
            ! the direction of the step is x_new-x_old
            tr_xold_gx = 0.0_dp
            tr_xnew_gx = 0.0_dp
            IF (do_ks) THEN
               DO ispin = 1, nspin
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                                 qs_ot_env(ispin)%matrix_gx, tmp)
                  tr_xold_gx = tr_xold_gx + tmp
                  CALL dbcsr_dot(qs_ot_env(ispin)%matrix_x, &
                                 qs_ot_env(ispin)%matrix_gx, tmp)
                  tr_xnew_gx = tr_xnew_gx + tmp
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                                    qs_ot_env(ispin)%matrix_gx_im, tmp)
                     tr_xold_gx = tr_xold_gx + tmp
                     CALL dbcsr_dot(qs_ot_env(ispin)%matrix_x_im, &
                                    qs_ot_env(ispin)%matrix_gx_im, tmp)
                     tr_xnew_gx = tr_xnew_gx + tmp
                  END IF
                  IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                    qs_ot_env(ispin)%rot_mat_gx, tmp)
                     tr_xold_gx = tr_xold_gx + 0.5_dp*tmp
                     CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_x, &
                                    qs_ot_env(ispin)%rot_mat_gx, tmp)
                     tr_xnew_gx = tr_xnew_gx + 0.5_dp*tmp
                     IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                        tr_xold_gx = tr_xold_gx + 0.5_dp*tmp
                        CALL dbcsr_dot(qs_ot_env(ispin)%rot_mat_x_im, &
                                       qs_ot_env(ispin)%rot_mat_gx_im, tmp)
                        tr_xnew_gx = tr_xnew_gx + 0.5_dp*tmp
                     END IF
                  END IF
               END DO
            END IF
            IF (do_ener) THEN
               DO ispin = 1, nspin
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_h_x(j, :), qs_ot_env(ispin)%ener_gx(:))
                  tr_xold_gx = tr_xold_gx + tmp
                  tmp = DOT_PRODUCT(qs_ot_env(ispin)%ener_x(:), qs_ot_env(ispin)%ener_gx(:))
                  tr_xnew_gx = tr_xnew_gx + tmp
               END DO
            END IF
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, tr_xold_gx)
            CALL ot_mini_sum_kpoint_scalar(para_env_inter_kp, tr_xnew_gx)
            overlap = (tr_xnew_gx - tr_xold_gx)
            ! OK, bad luck, take a SD step along the preconditioned gradient
            IF (overlap > 0.0_dp) THEN
               do_ot_sd = .TRUE.
            END IF
         END IF
      END IF

      IF (do_ot_sd) THEN
         qs_ot_env(1)%OT_METHOD_FULL = "OT SD"
         IF (do_ks) THEN
            DO ispin = 1, nspin
               CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                              1.0_dp, 1.0_dp)
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(j)%matrix, &
                              1.0_dp, 1.0_dp)
               IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%matrix_x_im, 0.0_dp)
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                 qs_ot_env(ispin)%matrix_h_e_im(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  CALL dbcsr_add(qs_ot_env(ispin)%matrix_x_im, &
                                 qs_ot_env(ispin)%matrix_h_x_im(j)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
               IF (qs_ot_env(ispin)%settings%do_rotation) THEN
                  CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x, 0.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_set(qs_ot_env(ispin)%rot_mat_x_im, 0.0_dp)
                  END IF
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                 qs_ot_env(ispin)%rot_mat_h_e(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                    qs_ot_env(ispin)%rot_mat_h_e_im(j)%matrix, &
                                    1.0_dp, 1.0_dp)
                  END IF
                  CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x, &
                                 qs_ot_env(ispin)%rot_mat_h_x(j)%matrix, &
                                 1.0_dp, 1.0_dp)
                  IF (qs_ot_env(ispin)%has_complex_kpoint_state) THEN
                     CALL dbcsr_add(qs_ot_env(ispin)%rot_mat_x_im, &
                                    qs_ot_env(ispin)%rot_mat_h_x_im(j)%matrix, &
                                    1.0_dp, 1.0_dp)
                  END IF
               END IF
            END DO
         END IF
         IF (do_ener) THEN
            DO ispin = 1, nspin
               qs_ot_env(ispin)%ener_x(:) = 0._dp
               qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) + qs_ot_env(ispin)%ener_h_e(j, :)
               qs_ot_env(ispin)%ener_x(:) = qs_ot_env(ispin)%ener_x(:) + qs_ot_env(ispin)%ener_h_x(j, :)
            END DO
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE ot_diis_step

! **************************************************************************************************
!> \brief Energy minimizer by Broyden's method
!> \param qs_ot_env variable to control minimizer behaviour
!> \author Kurt Baarman (09.2010)
! **************************************************************************************************
   SUBROUTINE ot_broyden_step(qs_ot_env)
      TYPE(qs_ot_type), DIMENSION(:), POINTER            :: qs_ot_env

      INTEGER                                            :: diis_bound, diis_m, i, ispin, itmp, j, &
                                                            k, n, nspin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: circ_index
      LOGICAL                                            :: adaptive_sigma, do_ener, do_ks, &
                                                            enable_flip, forget_history
      REAL(KIND=dp)                                      :: beta, eta, gamma, omega, sigma, &
                                                            sigma_dec, sigma_min, tmp, tmp2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: f, x
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: G, S
      TYPE(cp_logger_type), POINTER                      :: logger

      eta = qs_ot_env(1)%settings%broyden_eta
      omega = qs_ot_env(1)%settings%broyden_omega
      sigma_dec = qs_ot_env(1)%settings%broyden_sigma_decrease
      sigma_min = qs_ot_env(1)%settings%broyden_sigma_min
      forget_history = qs_ot_env(1)%settings%broyden_forget_history
      adaptive_sigma = qs_ot_env(1)%settings%broyden_adaptive_sigma
      enable_flip = qs_ot_env(1)%settings%broyden_enable_flip

      do_ks = qs_ot_env(1)%settings%ks
      do_ener = qs_ot_env(1)%settings%do_ener

      beta = qs_ot_env(1)%settings%broyden_beta
      gamma = qs_ot_env(1)%settings%broyden_gamma
      IF (adaptive_sigma) THEN
         IF (qs_ot_env(1)%broyden_adaptive_sigma < 0.0_dp) THEN
            sigma = qs_ot_env(1)%settings%broyden_sigma
         ELSE
            sigma = qs_ot_env(1)%broyden_adaptive_sigma
         END IF
      ELSE
         sigma = qs_ot_env(1)%settings%broyden_sigma
      END IF

      ! simplify our life....
      IF (do_ener .OR. .NOT. do_ks .OR. qs_ot_env(1)%settings%do_rotation .OR. &
          ANY(qs_ot_env(:)%has_complex_kpoint_state)) THEN
         CPABORT("Not yet implemented")
      END IF
      !
      nspin = SIZE(qs_ot_env)

      diis_m = qs_ot_env(1)%settings%diis_m

      IF (qs_ot_env(1)%diis_iter < diis_m) THEN
         diis_bound = qs_ot_env(1)%diis_iter + 1
      ELSE
         diis_bound = diis_m
      END IF

      ! We want x:s, f:s and one random vector
      k = 2*diis_bound + 1
      ALLOCATE (S(k, k))
      ALLOCATE (G(k, k))
      ALLOCATE (f(k))
      ALLOCATE (x(k))
      ALLOCATE (circ_index(diis_bound))
      G = 0.0
      DO i = 1, k
         G(i, i) = sigma
      END DO
      S = 0.0

      j = MOD(qs_ot_env(1)%diis_iter, diis_m) + 1 ! index in the circular array

      DO ispin = 1, nspin
         CALL dbcsr_copy(qs_ot_env(ispin)%matrix_h_x(j)%matrix, qs_ot_env(ispin)%matrix_x)
      END DO

      IF (ASSOCIATED(qs_ot_env(1)%preconditioner)) THEN
         qs_ot_env(1)%gnorm = 0.0_dp
         DO ispin = 1, nspin
            CALL apply_preconditioner(qs_ot_env(ispin)%preconditioner, &
                                      qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix)
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           tmp)
            qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
         END DO
         IF (qs_ot_env(1)%gnorm < 0.0_dp) THEN
            WRITE (cp_logger_get_default_unit_nr(logger), *) "WARNING Preconditioner not positive definite !"
         END IF
         DO ispin = 1, nspin
            CALL dbcsr_scale(qs_ot_env(ispin)%matrix_h_e(j)%matrix, -1.0_dp)
         END DO
      ELSE
         qs_ot_env(1)%gnorm = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, qs_ot_env(ispin)%matrix_gx, tmp)
            qs_ot_env(1)%gnorm = qs_ot_env(1)%gnorm + tmp
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_h_e(j)%matrix, &
                           qs_ot_env(ispin)%matrix_gx, alpha_scalar=0.0_dp, beta_scalar=-1.0_dp)
         END DO
      END IF

      k = 0
      n = 0
      CALL dbcsr_get_info(qs_ot_env(1)%matrix_x, nfullrows_total=n)
      DO ispin = 1, nspin
         CALL dbcsr_get_info(qs_ot_env(ispin)%matrix_x, nfullcols_total=itmp)
         k = k + itmp
      END DO

      ! Handling the case of no free variables to optimize
      IF (INT(n, KIND=int_8)*INT(k, KIND=int_8) /= 0) THEN
         qs_ot_env(1)%delta = SQRT(ABS(qs_ot_env(1)%gnorm)/(INT(n, KIND=int_8)*INT(k, KIND=int_8)))
         qs_ot_env(1)%gradient = -qs_ot_env(1)%gnorm
      ELSE
         qs_ot_env(1)%delta = 0.0_dp
         qs_ot_env(1)%gradient = 0.0_dp
      END IF

      IF (diis_bound == diis_m) THEN
         DO i = 1, diis_bound
            circ_index(i) = MOD(j + i - 1, diis_m) + 1
         END DO
      ELSE
         DO i = 1, diis_bound
            circ_index(i) = i
         END DO
      END IF

      S = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_init_random(qs_ot_env(ispin)%matrix_x)
         DO i = 1, diis_bound
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, tmp)
            S(i, i) = S(i, i) + tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, tmp)
            S(i + diis_bound, i + diis_bound) = S(i + diis_bound, i + diis_bound) + tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_x, tmp)
            S(i, 2*diis_bound + 1) = S(i, 2*diis_bound + 1) + tmp
            S(i, 2*diis_bound + 1) = S(2*diis_bound + 1, i)
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_x, tmp)
            S(i + diis_bound, 2*diis_bound + 1) = S(i + diis_bound, 2*diis_bound + 1) + tmp
            S(i + diis_bound, 2*diis_bound + 1) = S(2*diis_bound + 1, diis_bound + i)
            DO k = (i + 1), diis_bound
               CALL dbcsr_dot( &
                  qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                  qs_ot_env(ispin)%matrix_h_x(circ_index(k))%matrix, &
                  tmp)
               S(i, k) = S(i, k) + tmp
               S(k, i) = S(i, k)
               CALL dbcsr_dot( &
                  qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                  qs_ot_env(ispin)%matrix_h_e(circ_index(k))%matrix, &
                  tmp)
               S(diis_bound + i, diis_bound + k) = S(diis_bound + i, diis_bound + k) + tmp
               S(diis_bound + k, diis_bound + i) = S(diis_bound + i, diis_bound + k)
            END DO
            DO k = 1, diis_bound
               CALL dbcsr_dot( &
                  qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                  qs_ot_env(ispin)%matrix_h_e(circ_index(k))%matrix, tmp)
               S(i, k + diis_bound) = S(i, k + diis_bound) + tmp
               S(k + diis_bound, i) = S(i, k + diis_bound)
            END DO
         END DO
         CALL dbcsr_dot(qs_ot_env(ispin)%matrix_x, qs_ot_env(ispin)%matrix_x, tmp)
         S(2*diis_bound + 1, 2*diis_bound + 1) = S(2*diis_bound + 1, 2*diis_bound + 1) + tmp
      END DO

      ! normalize
      k = 2*diis_bound + 1
      tmp = SQRT(S(k, k))
      S(k, :) = S(k, :)/tmp
      S(:, k) = S(:, k)/tmp

      IF (diis_bound > 1) THEN
         tmp = 0.0_dp
         tmp2 = 0.0_dp
         i = diis_bound
         DO ispin = 1, nspin
            ! dot product of differences
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               tmp)
            tmp2 = tmp2 + tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i - 1))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
               tmp)
            tmp2 = tmp2 - tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i - 1))%matrix, &
               tmp)
            tmp2 = tmp2 - tmp
            CALL dbcsr_dot( &
               qs_ot_env(ispin)%matrix_h_x(circ_index(i - 1))%matrix, &
               qs_ot_env(ispin)%matrix_h_e(circ_index(i - 1))%matrix, &
               tmp)
            tmp2 = tmp2 + tmp
         END DO
         qs_ot_env(1)%c_broy(i - 1) = tmp2
      END IF

      qs_ot_env(1)%energy_h(j) = qs_ot_env(1)%etotal

      ! If we went uphill, do backtracking line search
      i = MINLOC(qs_ot_env(1)%energy_h(1:diis_bound), dim=1)
      IF (i /= j) THEN
         sigma = sigma_dec*sigma
         qs_ot_env(1)%OT_METHOD_FULL = "OT BTRK"
         DO ispin = 1, nspin
            CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                           qs_ot_env(ispin)%matrix_h_x(i)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=(1.0_dp - gamma))
            CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                           qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=gamma)
         END DO
      ELSE
         ! Construct G
         DO i = 2, diis_bound
            f = 0.0
            x = 0.0
            ! f is df_i
            x(i) = 1.0
            x(i - 1) = -1.0
            ! x is dx_i
            f(diis_bound + i) = 1.0
            f(diis_bound + i - 1) = -1.0
            tmp = 1.0_dp
            ! We want a pos def Hessian
            IF (enable_flip) THEN
               IF (qs_ot_env(1)%c_broy(i - 1) > 0) THEN
                  !qs_ot_env(1)%OT_METHOD_FULL="OT FLIP"
                  tmp = -1.0_dp
               END IF
            END IF

            ! get dx-Gdf
            x(:) = tmp*x - MATMUL(G, f)
            ! dfSdf
            ! we calculate matmul(S, f) twice. They're small...
            tmp = DOT_PRODUCT(f, MATMUL(S, f))
            ! NOTE THAT S IS SYMMETRIC !!!
            f(:) = MATMUL(S, f)/tmp
            ! the spread is an outer vector product
            G(:, :) = G + SPREAD(x, dim=2, ncopies=SIZE(f))*SPREAD(f, dim=1, ncopies=SIZE(x))
         END DO
         f = 0.0_dp
         f(2*diis_bound) = 1.0_dp
         x(:) = -beta*MATMUL(G, f)

         ! OK, add the vectors now, this sums up to the proposed step
         DO ispin = 1, nspin
            CALL dbcsr_set(qs_ot_env(ispin)%matrix_x, 0.0_dp)
            DO i = 1, diis_bound
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_e(circ_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-x(i + diis_bound))
            END DO
            DO i = 1, diis_bound
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(circ_index(i))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=x(i))
            END DO
         END DO

         IF (adaptive_sigma) THEN
            tmp = new_sigma(G, S, diis_bound)
            !tmp = tmp * qs_ot_env ( 1 ) % settings % broyden_sigma
            tmp = tmp*eta
            sigma = MIN(omega*sigma, tmp)
         END IF

         ! compute the inner product of direction of the step and gradient
         tmp = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_dot(qs_ot_env(ispin)%matrix_gx, &
                           qs_ot_env(ispin)%matrix_x, &
                           tmp2)
            tmp = tmp + tmp2
         END DO

         DO ispin = 1, nspin
            ! if the direction of the step is not in direction of the gradient,
            ! change step sign
            IF (tmp >= 0.0_dp) THEN
               qs_ot_env(1)%OT_METHOD_FULL = "OT TURN"
               IF (forget_history) THEN
                  qs_ot_env(1)%diis_iter = 0
               END IF
               sigma = sigma*sigma_dec
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix, &
                              alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
            ELSE
               CALL dbcsr_add(qs_ot_env(ispin)%matrix_x, &
                              qs_ot_env(ispin)%matrix_h_x(circ_index(diis_bound))%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END IF
         END DO
      END IF

      ! get rid of S, G, f, x, circ_index for next round
      DEALLOCATE (S, G, f, x, circ_index)

      ! update for next round
      qs_ot_env(1)%diis_iter = qs_ot_env(1)%diis_iter + 1
      qs_ot_env(1)%broyden_adaptive_sigma = MAX(sigma, sigma_min)

   END SUBROUTINE ot_broyden_step

! **************************************************************************************************
!> \brief ...
!> \param G ...
!> \param S ...
!> \param n ...
!> \return ...
! **************************************************************************************************
   FUNCTION new_sigma(G, S, n) RESULT(sigma)
!
! Calculate new sigma from eigenvalues of full size G by Arnoldi.
!
! **************************************************************************************************

      REAL(KIND=dp), DIMENSION(:, :)                     :: G, S
      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: sigma

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: H

      ALLOCATE (H(n, n))
      CALL hess_G(G, S, H, n)
      ALLOCATE (eigv(n))
      CALL diamat_all(H(1:n, 1:n), eigv)

      SELECT CASE (1)
      CASE (1)
         ! This estimator seems to work well. No theory.
         sigma = SUM(ABS(eigv**2))/SUM(ABS(eigv))
      CASE (2)
         ! Estimator based on Frobenius norm minimizer
         sigma = SUM(ABS(eigv))/MAX(1, SIZE(eigv))
      CASE (3)
         ! Estimator based on induced 2-norm
         sigma = (MAXVAL(ABS(eigv)) + MINVAL(ABS(eigv)))*0.5_dp
      END SELECT

      DEALLOCATE (H, eigv)
   END FUNCTION new_sigma

! **************************************************************************************************
!> \brief ...
!> \param G ...
!> \param S ...
!> \param H ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE hess_G(G, S, H, n)
!
! Make a hessenberg out of G into H. Cf Arnoldi.
! Inner product is weighted by S.
! Possible lucky breakdown at n.
!
! **************************************************************************************************
      REAL(KIND=dp), DIMENSION(:, :)                     :: G, S, H
      INTEGER                                            :: n

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: v
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Q

      i = SIZE(G, 1)
      k = SIZE(H, 1)
      ALLOCATE (Q(i, k))
      ALLOCATE (v(i))
      H = 0.0_dp
      Q = 0.0_dp

      Q(:, 1) = 1.0_dp
      tmp = SQRT(DOT_PRODUCT(Q(:, 1), MATMUL(S, Q(:, 1))))
      Q(:, :) = Q(:, :)/tmp

      DO i = 1, k
         v(:) = MATMUL(G, Q(:, i))
         DO j = 1, i
            H(j, i) = DOT_PRODUCT(Q(:, j), MATMUL(S, v))
            v(:) = v - H(j, i)*Q(:, j)
         END DO
         IF (i < k) THEN
            tmp = DOT_PRODUCT(v, MATMUL(S, v))
            IF (tmp <= 0.0_dp) THEN
               n = i
               EXIT
            END IF
            tmp = SQRT(tmp)
            ! Lucky breakdown
            IF (ABS(tmp) < 1e-9_dp) THEN
               n = i
               EXIT
            END IF
            H(i + 1, i) = tmp
            Q(:, i + 1) = v/H(i + 1, i)
         END IF
      END DO

      DEALLOCATE (Q, v)
   END SUBROUTINE hess_G

END MODULE qs_ot_minimizer
