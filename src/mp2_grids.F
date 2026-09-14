!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate frequency and time grids (integration points and weights)
!>        for correlation methods
!> \par History
!>      05.2019 Refactored from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE mp2_grids
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_type
   USE greenx_interface,                ONLY: greenx_get_minimax_grid
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
   USE minimax_exp_gw,                  ONLY: get_exp_minimax_coeff_gw
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff,&
                                              get_rpa_minimax_coeff_larger_grid
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE time_frequency_grid_types,       ONLY: time_frequency_grid_release,&
                                              time_frequency_grid_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_grids'

   INTEGER, PARAMETER :: cosine_time_to_frequency = 1, &
                         sine_time_to_frequency = 2, &
                         cosine_frequency_to_time = 3

   PUBLIC :: build_minimax_time_frequency_grid, get_minimax_grid, get_clenshaw_grid, test_least_square_ft

CONTAINS

! **************************************************************************************************
!> \brief Build a minimax time/frequency grid without caller-specific energy-range logic.
!> \param num_points Number of minimax points.
!> \param energy_min Lower end of the physical energy interval.
!> \param energy_max Upper end of the physical energy interval.
!> \param regularization Regularization used for the fitted transform weights.
!> \param num_points_per_magnitude Number of fitting points per decade.
!> \param grid Grid object to fill.
!> \param build_frequency Whether to construct the frequency components.
!> \param build_time Whether to construct the imaginary-time components.
!> \param build_transforms Whether to construct the time/frequency transform weights.
!> \param build_sine Whether to construct the sine time-to-frequency weights.
!> \param time_scaling Scaling applied to the imaginary-time abscissas.
!> \param time_weight_scaling Scaling applied to the imaginary-time weights.
!> \param max_fit_error Maximum fitting error across all requested transforms.
!> \param print_warning Whether the minimax coefficient routine prints warnings.
! **************************************************************************************************
   SUBROUTINE build_minimax_time_frequency_grid(num_points, energy_min, energy_max, regularization, &
                                                num_points_per_magnitude, grid, build_frequency, build_time, &
                                                build_transforms, build_sine, time_scaling, time_weight_scaling, &
                                                max_fit_error, print_warning)

      INTEGER, INTENT(IN)                                :: num_points
      REAL(KIND=dp), INTENT(IN)                          :: energy_min, energy_max, regularization
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude
      TYPE(time_frequency_grid_type), INTENT(INOUT)      :: grid
      LOGICAL, INTENT(IN)                                :: build_frequency, build_time, &
                                                            build_transforms, build_sine
      REAL(KIND=dp), INTENT(IN)                          :: time_scaling, time_weight_scaling
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: max_fit_error
      LOGICAL, INTENT(IN), OPTIONAL                      :: print_warning

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_minimax_time_frequency_grid'

      INTEGER                                            :: handle, ierr
      LOGICAL                                            :: my_print_warning
      REAL(KIND=dp)                                      :: e_range, max_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: points_and_weights

      CALL timeset(routineN, handle)

      CPASSERT(num_points > 0)
      CPASSERT(energy_min > 0.0_dp)
      CPASSERT(energy_max >= energy_min)
      CPASSERT(time_scaling > 0.0_dp)
      CPASSERT(time_weight_scaling > 0.0_dp)
      CPASSERT(.NOT. build_transforms .OR. (build_frequency .AND. build_time))
      CPASSERT(.NOT. build_sine .OR. build_transforms)
      max_error = 0.0_dp

      CALL time_frequency_grid_release(grid)

      e_range = energy_max/energy_min
      ALLOCATE (points_and_weights(2*num_points))

      my_print_warning = .TRUE.
      IF (PRESENT(print_warning)) my_print_warning = print_warning

      IF (build_frequency) THEN
         IF (num_points <= 20) THEN
            CALL get_rpa_minimax_coeff(num_points, e_range, points_and_weights, ierr, my_print_warning)
         ELSE
            CALL get_rpa_minimax_coeff_larger_grid(num_points, e_range, points_and_weights)
         END IF

         ALLOCATE (grid%frequency(num_points))
         ALLOCATE (grid%frequency_weights(num_points))
         grid%frequency(:) = points_and_weights(1:num_points)*energy_min
         grid%frequency_weights(:) = points_and_weights(num_points + 1:)*energy_min
         IF (num_points >= 26) grid%frequency_weights(:) = grid%frequency_weights(:)*4.0_dp
      END IF

      IF (build_time) THEN
         IF (num_points <= 20) THEN
            CALL get_exp_minimax_coeff(num_points, e_range, points_and_weights)
         ELSE
            CALL get_exp_minimax_coeff_gw(num_points, e_range, points_and_weights)
         END IF

         ALLOCATE (grid%imaginary_time(num_points))
         ALLOCATE (grid%time_weights_at_zero_frequency(num_points))
         grid%imaginary_time(:) = points_and_weights(1:num_points)/time_scaling/energy_min
         grid%time_weights_at_zero_frequency(:) = points_and_weights(num_points + 1:)/time_weight_scaling/energy_min

         IF (build_transforms) THEN
            ALLOCATE (grid%cosine_time_to_frequency_weights(num_points, num_points))
            CALL fit_l_sq_weights(cosine_time_to_frequency, num_points, grid%imaginary_time, &
                                  grid%cosine_time_to_frequency_weights, grid%frequency, energy_min, energy_max, &
                                  max_error, num_points_per_magnitude, regularization)

            ALLOCATE (grid%cosine_frequency_to_time_weights(num_points, num_points))
            CALL fit_l_sq_weights(cosine_frequency_to_time, num_points, grid%imaginary_time, &
                                  grid%cosine_frequency_to_time_weights, grid%frequency, energy_min, energy_max, &
                                  max_error, num_points_per_magnitude, regularization)

            IF (build_sine) THEN
               ALLOCATE (grid%sine_time_to_frequency_weights(num_points, num_points))
               CALL fit_l_sq_weights(sine_time_to_frequency, num_points, grid%imaginary_time, &
                                     grid%sine_time_to_frequency_weights, grid%frequency, energy_min, energy_max, &
                                     max_error, num_points_per_magnitude, regularization)
            END IF
         END IF
      END IF

      IF (PRESENT(max_fit_error)) max_fit_error = max_error

      DEALLOCATE (points_and_weights)

      CALL timestop(handle)

   END SUBROUTINE build_minimax_time_frequency_grid

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param unit_nr ...
!> \param homo ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param do_im_time ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param do_print ...
!> \param tau_tj ...
!> \param tau_wj ...
!> \param qs_env ...
!> \param do_gw_im_time ...
!> \param do_kpoints_cubic_RPA ...
!> \param e_fermi ...
!> \param tj ...
!> \param wj ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_cos_tf_w_to_t ...
!> \param weights_sin_tf_t_to_w ...
!> \param regularization ...
! **************************************************************************************************
   SUBROUTINE get_minimax_grid(para_env, unit_nr, homo, Eigenval, num_integ_points, &
                               do_im_time, do_ri_sos_laplace_mp2, do_print, tau_tj, tau_wj, qs_env, do_gw_im_time, &
                               do_kpoints_cubic_RPA, e_fermi, tj, wj, weights_cos_tf_t_to_w, &
                               weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, regularization)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: do_im_time, do_ri_sos_laplace_mp2, &
                                                            do_print
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: do_gw_im_time, do_kpoints_cubic_RPA
      REAL(KIND=dp), INTENT(OUT)                         :: e_fermi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: regularization

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_minimax_grid'
      INTEGER, PARAMETER                                 :: num_points_per_magnitude = 200

      INTEGER                                            :: handle, ierr, jquad
      REAL(KIND=dp)                                      :: E_Range, Emax, Emin, max_error_min, &
                                                            my_regularization, scaling
      TYPE(time_frequency_grid_type)                     :: grid

      CALL timeset(routineN, handle)

      CALL determine_energy_range(qs_env, para_env, homo, Eigenval, do_ri_sos_laplace_mp2, &
                                  do_kpoints_cubic_RPA, Emin, Emax, e_range, e_fermi)

      ! Open-shell uses ONE minimax grid for the combined [min gap, max span] over both spins (the
      ! superset covers each channel, so it is accurate; per-spin grids would only be more efficient).
      IF (SIZE(homo) > 1) THEN
         CALL cp_hint(__LOCATION__, &
                      "Open-shell RPA/GW uses one minimax grid spanning [min gap, max span] across "// &
                      "both spin channels; raise QUADRATURE_POINTS if QP convergence is marginal for "// &
                      "strongly spin-asymmetric systems.")
      END IF

      my_regularization = qs_env%mp2_env%ri_g0w0%regularization_minimax
      IF (PRESENT(regularization)) my_regularization = regularization

      CALL greenx_get_minimax_grid(unit_nr, num_integ_points, emin, emax, &
                                   tau_tj, tau_wj, my_regularization, &
                                   tj, wj, weights_cos_tf_t_to_w, &
                                   weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, ierr)

! Shortcut if Greenx was available and successful
      IF (ierr == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (PRESENT(regularization)) THEN
         IF (num_integ_points > 20 .AND. e_range < 100.0_dp) THEN
            IF (unit_nr > 0) THEN
               CALL cp_warn(__LOCATION__, &
                            "You requested a large minimax grid (> 20 points) for a small minimax range R (R < 100). "// &
                            "That may lead to numerical "// &
                            "instabilities when computing minimax grid weights. You can prevent small ranges by choosing "// &
                            "a larger basis set with higher angular momenta or alternatively using all-electron calculations.")
            END IF
         END IF

         scaling = 2.0_dp
         IF (do_ri_sos_laplace_mp2) scaling = 1.0_dp
         CALL build_minimax_time_frequency_grid(num_integ_points, Emin, Emax, my_regularization, &
                                                num_points_per_magnitude, grid, &
                                                build_frequency=.NOT. do_ri_sos_laplace_mp2, &
                                                build_time=do_im_time .OR. do_ri_sos_laplace_mp2, &
                                                build_transforms=do_im_time .AND. .NOT. do_ri_sos_laplace_mp2, &
                                                build_sine=do_im_time .AND. (.NOT. do_ri_sos_laplace_mp2) .AND. do_gw_im_time, &
                                                time_scaling=scaling, time_weight_scaling=scaling, &
                                                max_fit_error=max_error_min)

         IF (.NOT. do_ri_sos_laplace_mp2) THEN
            ALLOCATE (tj, SOURCE=grid%frequency)
            ALLOCATE (wj, SOURCE=grid%frequency_weights)
            IF (unit_nr > 0 .AND. do_print) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
                  "MINIMAX_INFO| Number of integration points:", num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "MINIMAX_INFO| Gap for the minimax approximation:", Emin
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "MINIMAX_INFO| Range for the minimax approximation:", e_range
               WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") "MINIMAX_INFO| Minimax parameters:", "Weights", "Abscissas"
               DO jquad = 1, num_integ_points
                  WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") wj(jquad)/Emin, tj(jquad)/Emin
               END DO
               CALL m_flush(unit_nr)
            END IF
         END IF

         IF (do_im_time .OR. do_ri_sos_laplace_mp2) THEN
            ALLOCATE (tau_tj, SOURCE=grid%imaginary_time)
            ALLOCATE (tau_wj, SOURCE=grid%time_weights_at_zero_frequency)
            IF (do_im_time .AND. .NOT. do_ri_sos_laplace_mp2) THEN
               ALLOCATE (weights_cos_tf_t_to_w, SOURCE=grid%cosine_time_to_frequency_weights)
               ALLOCATE (weights_cos_tf_w_to_t, SOURCE=grid%cosine_frequency_to_time_weights)
               IF (ALLOCATED(grid%sine_time_to_frequency_weights)) THEN
                  ALLOCATE (weights_sin_tf_t_to_w, SOURCE=grid%sine_time_to_frequency_weights)
               END IF
            END IF

            IF (unit_nr > 0 .AND. do_print) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "MINIMAX_INFO| Range for the minimax approximation:", e_range
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,F15.4)") &
                  "MINIMAX_INFO| Gap:", Emin
               WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
                  "MINIMAX_INFO| Minimax parameters of the time grid:", "Weights", "Abscissas"
               DO jquad = 1, num_integ_points
                  WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") tau_wj(jquad)*Emin, tau_tj(jquad)*Emin
               END DO
               CALL m_flush(unit_nr)
            END IF

            IF (unit_nr > 0 .AND. do_im_time .AND. do_gw_im_time .AND. .NOT. do_ri_sos_laplace_mp2) THEN
               WRITE (UNIT=unit_nr, FMT="(T3,A,T66,ES15.2)") &
                  "MINIMAX_INFO| Maximum deviation among requested minimax transform fits:", max_error_min
            END IF
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_minimax_grid

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param unit_nr ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param fm_mat_S ...
!> \param my_do_gw ...
!> \param ext_scaling ...
!> \param a_scaling ...
!> \param tj ...
!> \param wj ...
! **************************************************************************************************
   SUBROUTINE get_clenshaw_grid(para_env, para_env_RPA, unit_nr, homo, virtual, Eigenval, num_integ_points, &
                                num_integ_group, color_rpa_group, fm_mat_S, my_do_gw, &
                                ext_scaling, a_scaling, tj, wj)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_RPA
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      LOGICAL, INTENT(IN)                                :: my_do_gw
      REAL(KIND=dp), INTENT(IN)                          :: ext_scaling
      REAL(KIND=dp), INTENT(OUT)                         :: a_scaling
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tj, wj

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_clenshaw_grid'

      INTEGER                                            :: handle, jquad

      CALL timeset(routineN, handle)

      ! Now, start to prepare the different grid
      ALLOCATE (tj(num_integ_points))
      tj = 0.0_dp

      ALLOCATE (wj(num_integ_points))
      wj = 0.0_dp

      DO jquad = 1, num_integ_points - 1
         tj(jquad) = jquad*pi/(2.0_dp*num_integ_points)
         wj(jquad) = pi/(num_integ_points*SIN(tj(jquad))**2)
      END DO
      tj(num_integ_points) = pi/2.0_dp
      wj(num_integ_points) = pi/(2.0_dp*num_integ_points*SIN(tj(num_integ_points))**2)

      IF (my_do_gw .AND. ext_scaling > 0.0_dp) THEN
         a_scaling = ext_scaling
      ELSE
         CALL calc_scaling_factor(a_scaling, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                  num_integ_points, num_integ_group, color_rpa_group, &
                                  tj, wj, fm_mat_S)
      END IF

      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.5)') 'INTEG_INFO| Scaling parameter:', a_scaling

      wj(:) = wj(:)*a_scaling

      CALL timestop(handle)

   END SUBROUTINE get_clenshaw_grid

! **************************************************************************************************
!> \brief ...
!> \param a_scaling_ext ...
!> \param para_env ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param Eigenval ...
!> \param num_integ_points ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param tj_ext ...
!> \param wj_ext ...
!> \param fm_mat_S ...
! **************************************************************************************************
   SUBROUTINE calc_scaling_factor(a_scaling_ext, para_env, para_env_RPA, homo, virtual, Eigenval, &
                                  num_integ_points, num_integ_group, color_rpa_group, &
                                  tj_ext, wj_ext, fm_mat_S)
      REAL(KIND=dp), INTENT(OUT)                         :: a_scaling_ext
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_RPA
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      INTEGER, INTENT(IN)                                :: num_integ_points, num_integ_group, &
                                                            color_rpa_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj_ext, wj_ext
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_scaling_factor'

      INTEGER                                            :: handle, icycle, jquad, ncol_local, &
                                                            ncol_local_beta, nspins
      LOGICAL                                            :: my_open_shell
      REAL(KIND=dp) :: a_high, a_low, a_scaling, conv_param, eps, first_deriv, left_term, &
         right_term, right_term_ref, right_term_ref_beta, step
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cottj, D_ia, D_ia_beta, iaia_RI, &
                                                            iaia_RI_beta, M_ia, M_ia_beta
      TYPE(mp_para_env_type), POINTER                    :: para_env_col, para_env_col_beta

      CALL timeset(routineN, handle)

      nspins = SIZE(homo)
      my_open_shell = (nspins == 2)

      eps = 1.0E-10_dp

      ALLOCATE (cottj(num_integ_points))

      ! calculate the cotangent of the abscissa tj
      DO jquad = 1, num_integ_points
         cottj(jquad) = 1.0_dp/TAN(tj_ext(jquad))
      END DO

      CALL calc_ia_ia_integrals(para_env_RPA, homo(1), virtual(1), ncol_local, right_term_ref, Eigenval(:, 1, 1), &
                                D_ia, iaia_RI, M_ia, fm_mat_S(1), para_env_col)

      ! In the open shell case do point 1-2-3 for the beta spin
      IF (my_open_shell) THEN
         CALL calc_ia_ia_integrals(para_env_RPA, homo(2), virtual(2), ncol_local_beta, right_term_ref_beta, Eigenval(:, 1, 2), &
                                   D_ia_beta, iaia_RI_beta, M_ia_beta, fm_mat_S(2), para_env_col_beta)

         right_term_ref = right_term_ref + right_term_ref_beta
      END IF

      ! bcast the result
      IF (para_env%mepos == 0) THEN
         CALL para_env%bcast(right_term_ref, 0)
      ELSE
         right_term_ref = 0.0_dp
         CALL para_env%bcast(right_term_ref, 0)
      END IF

      ! 5) start iteration for solving the non-linear equation by bisection
      ! find limit, here step=0.5 seems a good compromise
      conv_param = 100.0_dp*EPSILON(right_term_ref)
      step = 0.5_dp
      a_low = 0.0_dp
      a_high = step
      right_term = -right_term_ref
      DO icycle = 1, num_integ_points*2
         a_scaling = a_high

         CALL calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                M_ia, cottj, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                ncol_local, ncol_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_col, para_env_col_beta)
         left_term = left_term/4.0_dp/pi*a_scaling

         IF (ABS(left_term) > ABS(right_term) .OR. ABS(left_term + right_term) <= conv_param) EXIT
         a_low = a_high
         a_high = a_high + step

      END DO

      IF (ABS(left_term + right_term) >= conv_param) THEN
         IF (a_scaling >= 2*num_integ_points*step) THEN
            a_scaling = 1.0_dp
         ELSE

            DO icycle = 1, num_integ_points*2
               a_scaling = (a_low + a_high)/2.0_dp

               CALL calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                      M_ia, cottj, wj_ext, D_ia, D_ia_beta, M_ia_beta, &
                                      ncol_local, ncol_local_beta, num_integ_group, color_rpa_group, &
                                      para_env, para_env_col, para_env_col_beta)
               left_term = left_term/4.0_dp/pi*a_scaling

               IF (ABS(left_term) > ABS(right_term)) THEN
                  a_high = a_scaling
               ELSE
                  a_low = a_scaling
               END IF

               IF (ABS(a_high - a_low) < 1.0e-5_dp) EXIT

            END DO

         END IF
      END IF

      a_scaling_ext = a_scaling
      CALL para_env%bcast(a_scaling_ext, 0)

      DEALLOCATE (cottj)
      DEALLOCATE (iaia_RI)
      DEALLOCATE (D_ia)
      DEALLOCATE (M_ia)
      CALL mp_para_env_release(para_env_col)

      IF (my_open_shell) THEN
         DEALLOCATE (iaia_RI_beta)
         DEALLOCATE (D_ia_beta)
         DEALLOCATE (M_ia_beta)
         CALL mp_para_env_release(para_env_col_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calc_scaling_factor

! **************************************************************************************************
!> \brief ...
!> \param para_env_RPA ...
!> \param homo ...
!> \param virtual ...
!> \param ncol_local ...
!> \param right_term_ref ...
!> \param Eigenval ...
!> \param D_ia ...
!> \param iaia_RI ...
!> \param M_ia ...
!> \param fm_mat_S ...
!> \param para_env_col ...
! **************************************************************************************************
   SUBROUTINE calc_ia_ia_integrals(para_env_RPA, homo, virtual, ncol_local, right_term_ref, Eigenval, &
                                   D_ia, iaia_RI, M_ia, fm_mat_S, para_env_col)

      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_RPA
      INTEGER, INTENT(IN)                                :: homo, virtual
      INTEGER, INTENT(OUT)                               :: ncol_local
      REAL(KIND=dp), INTENT(OUT)                         :: right_term_ref
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: D_ia, iaia_RI, M_ia
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S
      TYPE(mp_para_env_type), POINTER                    :: para_env_col

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_ia_ia_integrals'

      INTEGER                                            :: avirt, color_col, color_row, handle, &
                                                            i_global, iiB, iocc, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: eigen_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: iaia_RI_dp
      TYPE(mp_para_env_type), POINTER                    :: para_env_row

      CALL timeset(routineN, handle)

      ! calculate the (ia|ia) RI integrals
      ! ----------------------------------
      ! 1) get info fm_mat_S
      CALL cp_fm_get_info(matrix=fm_mat_S, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! allocate the local buffer of iaia_RI integrals (dp kind)
      ALLOCATE (iaia_RI_dp(ncol_local))
      iaia_RI_dp = 0.0_dp

      ! 2) perform the local multiplication SUM_K (ia|K)*(ia|K)
      DO iiB = 1, ncol_local
         iaia_RI_dp(iiB) = iaia_RI_dp(iiB) + DOT_PRODUCT(fm_mat_S%local_data(:, iiB), fm_mat_S%local_data(:, iiB))
      END DO

      ! 3) sum the result with the processes of the RPA_group having the same columns
      !          _______ia______               _
      !         |   |   |   |   |             | |
      !     --> | 1 | 5 | 9 | 13|   SUM -->   | |
      !         |___|__ |___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 2 | 6 | 10| 14|   SUM -->   | |
      !       K |___|___|___|___|             |_|   (ia|ia)_RI
      !         |   |   |   |   |             | |
      !     --> | 3 | 7 | 11| 15|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !         |   |   |   |   |             | |
      !     --> | 4 | 8 | 12| 16|   SUM -->   | |
      !         |___|___|___|___|             |_|
      !

      color_col = fm_mat_S%matrix_struct%context%mepos(2)
      ALLOCATE (para_env_col)
      CALL para_env_col%from_split(para_env_RPA, color_col)

      CALL para_env_col%sum(iaia_RI_dp)

      ! convert the iaia_RI_dp into double-double precision
      ALLOCATE (iaia_RI(ncol_local))
      DO iiB = 1, ncol_local
         iaia_RI(iiB) = iaia_RI_dp(iiB)
      END DO
      DEALLOCATE (iaia_RI_dp)

      ! 4) calculate the right hand term, D_ia is the matrix containing the
      ! orbital energy differences, M_ia is the diagonal of the full RPA 'excitation'
      ! matrix
      ALLOCATE (D_ia(ncol_local))

      ALLOCATE (M_ia(ncol_local))

      DO iiB = 1, ncol_local
         i_global = col_indices(iiB)

         iocc = MAX(1, i_global - 1)/virtual + 1
         avirt = i_global - (iocc - 1)*virtual
         eigen_diff = Eigenval(avirt + homo) - Eigenval(iocc)

         D_ia(iiB) = eigen_diff
      END DO

      DO iiB = 1, ncol_local
         M_ia(iiB) = D_ia(iiB)*D_ia(iiB) + 2.0_dp*D_ia(iiB)*iaia_RI(iiB)
      END DO

      right_term_ref = 0.0_dp
      DO iiB = 1, ncol_local
         right_term_ref = right_term_ref + (SQRT(M_ia(iiB)) - D_ia(iiB) - iaia_RI(iiB))
      END DO
      right_term_ref = right_term_ref/2.0_dp

      ! sum the result with the processes of the RPA_group having the same row
      color_row = fm_mat_S%matrix_struct%context%mepos(1)
      ALLOCATE (para_env_row)
      CALL para_env_row%from_split(para_env_RPA, color_row)

      ! allocate communication array for rows
      CALL para_env_row%sum(right_term_ref)

      CALL mp_para_env_release(para_env_row)

      CALL timestop(handle)

   END SUBROUTINE calc_ia_ia_integrals

! **************************************************************************************************
!> \brief ...
!> \param a_scaling ...
!> \param left_term ...
!> \param first_deriv ...
!> \param num_integ_points ...
!> \param my_open_shell ...
!> \param M_ia ...
!> \param cottj ...
!> \param wj ...
!> \param D_ia ...
!> \param D_ia_beta ...
!> \param M_ia_beta ...
!> \param ncol_local ...
!> \param ncol_local_beta ...
!> \param num_integ_group ...
!> \param color_rpa_group ...
!> \param para_env ...
!> \param para_env_col ...
!> \param para_env_col_beta ...
! **************************************************************************************************
   SUBROUTINE calculate_objfunc(a_scaling, left_term, first_deriv, num_integ_points, my_open_shell, &
                                M_ia, cottj, wj, D_ia, D_ia_beta, M_ia_beta, &
                                ncol_local, ncol_local_beta, num_integ_group, color_rpa_group, &
                                para_env, para_env_col, para_env_col_beta)
      REAL(KIND=dp), INTENT(IN)                          :: a_scaling
      REAL(KIND=dp), INTENT(INOUT)                       :: left_term, first_deriv
      INTEGER, INTENT(IN)                                :: num_integ_points
      LOGICAL, INTENT(IN)                                :: my_open_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: M_ia, cottj, wj, D_ia, D_ia_beta, &
                                                            M_ia_beta
      INTEGER, INTENT(IN)                                :: ncol_local, ncol_local_beta, &
                                                            num_integ_group, color_rpa_group
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_col
      TYPE(mp_para_env_type), POINTER                    :: para_env_col_beta

      INTEGER                                            :: iiB, jquad
      REAL(KIND=dp)                                      :: first_deriv_beta, left_term_beta, omega

      left_term = 0.0_dp
      first_deriv = 0.0_dp
      left_term_beta = 0.0_dp
      first_deriv_beta = 0.0_dp
      DO jquad = 1, num_integ_points
         ! parallelize over integration points
         IF (MODULO(jquad, num_integ_group) /= color_rpa_group) CYCLE
         omega = a_scaling*cottj(jquad)

         DO iiB = 1, ncol_local
            ! parallelize over ia elements in the para_env_row group
            IF (MODULO(iiB, para_env_col%num_pe) /= para_env_col%mepos) CYCLE
            ! calculate left_term
            left_term = left_term + wj(jquad)* &
                        (LOG(1.0_dp + (M_ia(iiB) - D_ia(iiB)**2)/(omega**2 + D_ia(iiB)**2)) - &
                         (M_ia(iiB) - D_ia(iiB)**2)/(omega**2 + D_ia(iiB)**2))
            first_deriv = first_deriv + wj(jquad)*cottj(jquad)**2* &
                          ((-M_ia(iiB) + D_ia(iiB)**2)**2/((omega**2 + D_ia(iiB)**2)**2*(omega**2 + M_ia(iiB))))
         END DO

         IF (my_open_shell) THEN
            DO iiB = 1, ncol_local_beta
               ! parallelize over ia elements in the para_env_row group
               IF (MODULO(iiB, para_env_col_beta%num_pe) /= para_env_col_beta%mepos) CYCLE
               ! calculate left_term
               left_term_beta = left_term_beta + wj(jquad)* &
                                (LOG(1.0_dp + (M_ia_beta(iiB) - D_ia_beta(iiB)**2)/(omega**2 + D_ia_beta(iiB)**2)) - &
                                 (M_ia_beta(iiB) - D_ia_beta(iiB)**2)/(omega**2 + D_ia_beta(iiB)**2))
               first_deriv_beta = &
                  first_deriv_beta + wj(jquad)*cottj(jquad)**2* &
                  ((-M_ia_beta(iiB) + D_ia_beta(iiB)**2)**2/((omega**2 + D_ia_beta(iiB)**2)**2*(omega**2 + M_ia_beta(iiB))))
            END DO
         END IF

      END DO

      ! sum the contribution from all proc, starting form the row group
      CALL para_env%sum(left_term)
      CALL para_env%sum(first_deriv)

      IF (my_open_shell) THEN
         CALL para_env%sum(left_term_beta)
         CALL para_env%sum(first_deriv_beta)

         left_term = left_term + left_term_beta
         first_deriv = first_deriv + first_deriv_beta
      END IF

   END SUBROUTINE calculate_objfunc

! **************************************************************************************************
!> \brief Calculate least-squares weights for a time/frequency transform.
!> \param transform_kind Type of transform to fit.
!> \param num_integ_points Number of integration points.
!> \param tau_tj Imaginary-time integration points.
!> \param weights Transform weights to construct.
!> \param omega_tj Frequency integration points.
!> \param E_min Lower end of the fitting interval.
!> \param E_max Upper end of the fitting interval.
!> \param max_error Maximum fitting error.
!> \param num_points_per_magnitude Number of fitting points per decade.
!> \param regularization Regularization parameter for the pseudoinverse.
! **************************************************************************************************
   SUBROUTINE fit_l_sq_weights(transform_kind, num_integ_points, tau_tj, weights, omega_tj, E_min, E_max, &
                               max_error, num_points_per_magnitude, regularization)

      INTEGER, INTENT(IN)                                :: transform_kind, num_integ_points
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tau_tj
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: weights
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: omega_tj
      REAL(KIND=dp), INTENT(IN)                          :: E_min, E_max
      REAL(KIND=dp), INTENT(OUT)                         :: max_error
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude
      REAL(KIND=dp), INTENT(IN)                          :: regularization

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fit_l_sq_weights'

      INTEGER                                            :: handle, iii, info, jjj, jquad, lwork, &
                                                            num_x_nodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iwork
      LOGICAL                                            :: sine_transform, time_to_frequency
      REAL(KIND=dp)                                      :: func_val, max_error_tmp, multiplicator, &
                                                            omega, tau, x_value
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: sing_values, vec_UTy, weight_work, work, &
                                                            x_values, y_values
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mat_A, mat_SinvVSinvSigma, &
                                                            mat_SinvVSinvT, mat_U

      CALL timeset(routineN, handle)

      CPASSERT(transform_kind >= cosine_time_to_frequency)
      CPASSERT(transform_kind <= cosine_frequency_to_time)
      CPASSERT(SIZE(tau_tj) == num_integ_points)
      CPASSERT(SIZE(omega_tj) == num_integ_points)
      CPASSERT(SIZE(weights, 1) == num_integ_points)
      CPASSERT(SIZE(weights, 2) == num_integ_points)

      time_to_frequency = transform_kind /= cosine_frequency_to_time
      sine_transform = transform_kind == sine_time_to_frequency

      ! take num_points_per_magnitude points per 10-interval
      num_x_nodes = (INT(LOG10(E_max/E_min)) + 1)*num_points_per_magnitude

      ! take at least as many x points as integration points to have clear
      ! input for the singular value decomposition
      num_x_nodes = MAX(num_x_nodes, num_integ_points)

      ALLOCATE (x_values(num_x_nodes), y_values(num_x_nodes), mat_A(num_x_nodes, num_integ_points), &
                weight_work(num_integ_points), sing_values(num_integ_points), mat_U(num_x_nodes, num_x_nodes), &
                mat_SinvVSinvT(num_x_nodes, num_integ_points), work(8*num_integ_points*num_integ_points + &
                                                                  12*num_integ_points + 2*num_x_nodes), iwork(8*num_integ_points), &
                mat_SinvVSinvSigma(num_integ_points, num_x_nodes), vec_UTy(num_x_nodes))
      mat_SinvVSinvSigma = 0.0_dp

      ! double the value nessary for 'A' to achieve good performance
      lwork = SIZE(work)

      multiplicator = (E_max/E_min)**(1.0_dp/(REAL(num_x_nodes, KIND=dp) - 1.0_dp))
      DO iii = 1, num_x_nodes
         x_values(iii) = E_min*multiplicator**(iii - 1)
      END DO

      max_error = 0.0_dp

      DO jquad = 1, num_integ_points

         IF (time_to_frequency) THEN
            omega = omega_tj(jquad)

            IF (sine_transform) THEN
               ! y=2*omega/(x^2+omega^2)
               DO iii = 1, num_x_nodes
                  y_values(iii) = 2.0_dp*omega/(x_values(iii)**2 + omega**2)
               END DO
            ELSE
               ! y=2*x/(x^2+omega^2)
               DO iii = 1, num_x_nodes
                  y_values(iii) = 2.0_dp*x_values(iii)/(x_values(iii)**2 + omega**2)
               END DO
            END IF

            DO jjj = 1, num_integ_points
               DO iii = 1, num_x_nodes
                  IF (sine_transform) THEN
                     mat_A(iii, jjj) = SIN(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
                  ELSE
                     mat_A(iii, jjj) = COS(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
                  END IF
               END DO
            END DO
         ELSE
            tau = tau_tj(jquad)

            ! y=exp(-x*|tau|)
            DO iii = 1, num_x_nodes
               y_values(iii) = EXP(-x_values(iii)*tau)
            END DO

            DO jjj = 1, num_integ_points
               omega = omega_tj(jjj)
               DO iii = 1, num_x_nodes
                  x_value = x_values(iii)
                  mat_A(iii, jjj) = COS(tau*omega)*2.0_dp*x_value/(x_value**2 + omega**2)
               END DO
            END DO
         END IF

         ! Singular value decomposition of mat_A
         CALL DGESDD('A', num_x_nodes, num_integ_points, mat_A, num_x_nodes, sing_values, mat_U, num_x_nodes, &
                     mat_SinvVSinvT, num_x_nodes, work, lwork, iwork, info)
         CPASSERT(info == 0)

         ! integration weights = V Sigma U^T y
         DO jjj = 1, num_integ_points
            DO iii = 1, num_integ_points
               mat_SinvVSinvSigma(iii, jjj) = mat_SinvVSinvT(jjj, iii)*sing_values(jjj) &
                                              /(regularization**2 + sing_values(jjj)**2)
            END DO
         END DO

         CALL DGEMM('T', 'N', num_x_nodes, 1, num_x_nodes, 1.0_dp, mat_U, num_x_nodes, y_values, num_x_nodes, &
                    0.0_dp, vec_UTy, num_x_nodes)
         CALL DGEMM('N', 'N', num_integ_points, 1, num_x_nodes, 1.0_dp, mat_SinvVSinvSigma, num_integ_points, &
                    vec_UTy, num_x_nodes, 0.0_dp, weight_work, num_integ_points)

         weights(jquad, :) = weight_work(:)

         max_error_tmp = 0.0_dp
         DO iii = 1, num_x_nodes
            func_val = 0.0_dp
            IF (time_to_frequency) THEN
               DO jjj = 1, num_integ_points
                  IF (sine_transform) THEN
                     func_val = func_val + weight_work(jjj)*SIN(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
                  ELSE
                     func_val = func_val + weight_work(jjj)*COS(omega*tau_tj(jjj))*EXP(-x_values(iii)*tau_tj(jjj))
                  END IF
               END DO
            ELSE
               DO jjj = 1, num_integ_points
                  omega = omega_tj(jjj)
                  func_val = func_val + weight_work(jjj)*COS(tau*omega)*2.0_dp*x_values(iii) &
                             /(x_values(iii)**2 + omega**2)
               END DO
            END IF
            max_error_tmp = MAX(max_error_tmp, ABS(y_values(iii) - func_val))
         END DO
         max_error = MAX(max_error, max_error_tmp)

      END DO

      DEALLOCATE (x_values, y_values, mat_A, weight_work, sing_values, mat_U, mat_SinvVSinvT, work, iwork, &
                  mat_SinvVSinvSigma, vec_UTy)

      CALL timestop(handle)

   END SUBROUTINE fit_l_sq_weights

! **************************************************************************************************
!> \brief test the singular value decomposition for the computation of integration weights for the
!>         Fourier transform between time and frequency grid in cubic-scaling RPA
!> \param nR ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE test_least_square_ft(nR, iw)
      INTEGER, INTENT(IN)                                :: nR, iw

      INTEGER                                            :: ierr, iR, num_integ_points
      REAL(KIND=dp)                                      :: max_error, multiplicator, Rc, Rc_max
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tj, x_tw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_t_to_w

      Rc_max = 1.0E+7

      multiplicator = Rc_max**(1.0_dp/(REAL(nR, KIND=dp) - 1.0_dp))

      DO num_integ_points = 1, 20

         ALLOCATE (x_tw(2*num_integ_points))
         x_tw = 0.0_dp
         ALLOCATE (tau_tj(num_integ_points))
         tau_tj = 0.0_dp
         ALLOCATE (weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
         weights_cos_tf_t_to_w = 0.0_dp
         ALLOCATE (tj(num_integ_points))
         tj = 0.0_dp

         DO iR = 0, nR - 1

            Rc = 2.0_dp*multiplicator**iR

            ierr = 0
            CALL get_rpa_minimax_coeff(num_integ_points, Rc, x_tw, ierr, print_warning=.FALSE.)

            tj(:) = x_tw(1:num_integ_points)

            x_tw = 0.0_dp

            CALL get_exp_minimax_coeff(num_integ_points, Rc, x_tw)

            tau_tj(:) = x_tw(1:num_integ_points)/2.0_dp

            CALL fit_l_sq_weights(cosine_time_to_frequency, num_integ_points, tau_tj, &
                                  weights_cos_tf_t_to_w, tj, 1.0_dp, Rc, max_error, 200, 0.0_dp)

            IF (iw > 0) THEN
               WRITE (iw, '(T2, I3, F12.1, ES12.3)') num_integ_points, Rc, max_error
            END IF

         END DO

         DEALLOCATE (x_tw, tau_tj, weights_cos_tf_t_to_w, tj)

      END DO

   END SUBROUTINE test_least_square_ft
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param gap ...
!> \param max_eig_diff ...
!> \param e_fermi ...
! **************************************************************************************************
   SUBROUTINE gap_and_max_eig_diff_kpoints(qs_env, para_env, gap, max_eig_diff, e_fermi)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), INTENT(OUT)                         :: gap, max_eig_diff, e_fermi

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gap_and_max_eig_diff_kpoints'

      INTEGER                                            :: handle, homo, ikpgr, ispin, kplocal, &
                                                            nmo, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp)                                      :: e_homo, e_homo_temp, e_lumo, e_lumo_temp
      REAL(KIND=dp), DIMENSION(3)                        :: tmp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      kpoints=kpoint)

      mo_set => kpoint%kp_env(1)%kpoint_env%mos(1, 1)
      CALL get_mo_set(mo_set, nmo=nmo)

      CALL get_kpoint_info(kpoint, kp_range=kp_range)
      kplocal = kp_range(2) - kp_range(1) + 1

      gap = 1000.0_dp
      max_eig_diff = 0.0_dp
      e_homo = -1000.0_dp
      e_lumo = 1000.0_dp

      DO ikpgr = 1, kplocal
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         nspin = SIZE(kp%mos, 2)
         DO ispin = 1, nspin
            mo_set => kp%mos(1, ispin)
            CALL get_mo_set(mo_set, eigenvalues=eigenvalues, homo=homo)
            e_homo_temp = eigenvalues(homo)
            e_lumo_temp = eigenvalues(homo + 1)

            IF (e_homo_temp > e_homo) e_homo = e_homo_temp
            IF (e_lumo_temp < e_lumo) e_lumo = e_lumo_temp
            IF (eigenvalues(nmo) - eigenvalues(1) > max_eig_diff) max_eig_diff = eigenvalues(nmo) - eigenvalues(1)

         END DO
      END DO

      ! Collect all three numbers in an array
      ! Reverse sign of lumo to reduce number of MPI calls
      tmp(1) = e_homo
      tmp(2) = -e_lumo
      tmp(3) = max_eig_diff
      CALL para_env%max(tmp)

      gap = -tmp(2) - tmp(1)
      e_fermi = (tmp(1) - tmp(2))*0.5_dp
      max_eig_diff = tmp(3)

      CALL timestop(handle)

   END SUBROUTINE gap_and_max_eig_diff_kpoints

! **************************************************************************************************
!> \brief returns minimal and maximal energy values for the E_range for the minimax grid selection
!> \param qs_env ...
!> \param para_env ...
!> \param homo index of the homo level for the respective spin channel
!> \param Eigenval eigenvalues
!> \param do_ri_sos_laplace_mp2 flag for SOS-MP2
!> \param do_kpoints_cubic_RPA flag for cubic-scaling RPA with k-points
!> \param Emin minimal eigenvalue difference (gap of the system)
!> \param Emax maximal eigenvalue difference
!> \param e_range ...
!> \param e_fermi Fermi level
! **************************************************************************************************
   SUBROUTINE determine_energy_range(qs_env, para_env, homo, Eigenval, do_ri_sos_laplace_mp2, &
                                     do_kpoints_cubic_RPA, Emin, Emax, e_range, e_fermi)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: Eigenval
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2, &
                                                            do_kpoints_cubic_RPA
      REAL(KIND=dp), INTENT(OUT)                         :: Emin, Emax, e_range, e_fermi

      CHARACTER(LEN=*), PARAMETER :: routineN = 'determine_energy_range'

      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: my_do_kpoints
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      ! Test for spin unrestricted
      nspins = SIZE(homo)

      ! Test whether all necessary variables are available
      my_do_kpoints = .FALSE.
      IF (.NOT. do_ri_sos_laplace_mp2) THEN
         my_do_kpoints = do_kpoints_cubic_RPA
      END IF

      IF (my_do_kpoints) THEN
         CALL gap_and_max_eig_diff_kpoints(qs_env, para_env, Emin, Emax, e_fermi)
         E_Range = Emax/Emin
      ELSE
         IF (qs_env%mp2_env%E_range <= 1.0_dp .OR. qs_env%mp2_env%E_gap <= 0.0_dp) THEN
            Emin = HUGE(dp)
            Emax = 0.0_dp
            DO ispin = 1, nspins
               IF (homo(ispin) > 0) THEN
                  Emin = MIN(Emin, Eigenval(homo(ispin) + 1, 1, ispin) - Eigenval(homo(ispin), 1, ispin))
                  Emax = MAX(Emax, MAXVAL(Eigenval(:, :, ispin)) - MINVAL(Eigenval(:, :, ispin)))
               END IF
            END DO
            E_Range = Emax/Emin
            qs_env%mp2_env%e_range = e_range
            qs_env%mp2_env%e_gap = Emin

            CALL get_qs_env(qs_env, input=input)
            CALL section_vals_val_set(input, "DFT%XC%WF_CORRELATION%E_RANGE", r_val=e_range)
            CALL section_vals_val_set(input, "DFT%XC%WF_CORRELATION%E_GAP", r_val=emin)
         ELSE
            E_range = qs_env%mp2_env%E_range
            Emin = qs_env%mp2_env%E_gap
            Emax = Emin*E_range
         END IF
      END IF

      ! When we perform SOS-MP2, we need an additional factor of 2 for the energies (compare with mp2_laplace.F)
      ! We do not need weights etc. for the cosine transform
      ! We do not scale Emax because it is not needed for SOS-MP2
      IF (do_ri_sos_laplace_mp2) THEN
         Emin = Emin*2.0_dp
         Emax = Emax*2.0_dp
      END IF

      CALL timestop(handle)
   END SUBROUTINE determine_energy_range

END MODULE mp2_grids
