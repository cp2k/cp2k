!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition and construction of time/frequency grids for correlation methods.
!> \par History
!>      05.2019 Refactored from rpa_ri_gpw [Frederick Stein]
! **************************************************************************************************
MODULE time_frequency_grids
   USE greenx_interface,                ONLY: greenx_get_minimax_grid
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE minimax_exp,                     ONLY: get_exp_minimax_coeff
   USE minimax_exp_gw,                  ONLY: get_exp_minimax_coeff_gw
   USE minimax_rpa,                     ONLY: get_rpa_minimax_coeff,&
                                              get_rpa_minimax_coeff_larger_grid
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'time_frequency_grids'

   TYPE time_frequency_grid_type
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: frequency
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: frequency_weights
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: imaginary_time
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: time_weights_at_zero_frequency
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: cosine_time_to_frequency_weights
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: cosine_frequency_to_time_weights
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: sine_time_to_frequency_weights
   END TYPE time_frequency_grid_type

   INTEGER, PARAMETER :: cosine_time_to_frequency = 1, &
                         sine_time_to_frequency = 2, &
                         cosine_frequency_to_time = 3

   PUBLIC :: build_clenshaw_grid, build_minimax_time_frequency_grid, test_least_square_ft, &
             time_frequency_grid_release, time_frequency_grid_type

CONTAINS

! **************************************************************************************************
!> \brief Release all data owned by a time_frequency_grid_type object.
!> \param grid The grid object to release.
! **************************************************************************************************
   SUBROUTINE time_frequency_grid_release(grid)

      TYPE(time_frequency_grid_type), INTENT(INOUT)      :: grid

      IF (ALLOCATED(grid%frequency)) DEALLOCATE (grid%frequency)
      IF (ALLOCATED(grid%frequency_weights)) DEALLOCATE (grid%frequency_weights)
      IF (ALLOCATED(grid%imaginary_time)) DEALLOCATE (grid%imaginary_time)
      IF (ALLOCATED(grid%time_weights_at_zero_frequency)) DEALLOCATE (grid%time_weights_at_zero_frequency)
      IF (ALLOCATED(grid%cosine_time_to_frequency_weights)) DEALLOCATE (grid%cosine_time_to_frequency_weights)
      IF (ALLOCATED(grid%cosine_frequency_to_time_weights)) DEALLOCATE (grid%cosine_frequency_to_time_weights)
      IF (ALLOCATED(grid%sine_time_to_frequency_weights)) DEALLOCATE (grid%sine_time_to_frequency_weights)

   END SUBROUTINE time_frequency_grid_release

! **************************************************************************************************
!> \brief Build a Clenshaw-Curtis frequency grid.
!> \param num_points Number of integration points.
!> \param grid Grid object to fill.
! **************************************************************************************************
   SUBROUTINE build_clenshaw_grid(num_points, grid)

      INTEGER, INTENT(IN)                                :: num_points
      TYPE(time_frequency_grid_type), INTENT(OUT)        :: grid

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_clenshaw_grid'

      INTEGER                                            :: handle, jquad

      CALL timeset(routineN, handle)

      CPASSERT(num_points > 0)

      ALLOCATE (grid%frequency(num_points), grid%frequency_weights(num_points))
      grid%frequency = 0.0_dp
      grid%frequency_weights = 0.0_dp

      DO jquad = 1, num_points - 1
         grid%frequency(jquad) = jquad*pi/(2.0_dp*num_points)
         grid%frequency_weights(jquad) = pi/(num_points*SIN(grid%frequency(jquad))**2)
      END DO
      grid%frequency(num_points) = pi/2.0_dp
      grid%frequency_weights(num_points) = &
         pi/(2.0_dp*num_points*SIN(grid%frequency(num_points))**2)

      CALL timestop(handle)

   END SUBROUTINE build_clenshaw_grid

! **************************************************************************************************
!> \brief Build a minimax time/frequency grid through the common backend boundary.
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
!> \param unit_nr Output unit used by an external backend.
!> \param prefer_external_backend Whether to try an available external backend first.
!> \param used_external_backend Whether the external backend supplied the grid.
! **************************************************************************************************
   SUBROUTINE build_minimax_time_frequency_grid(num_points, energy_min, energy_max, regularization, &
                                                num_points_per_magnitude, grid, build_frequency, build_time, &
                                                build_transforms, build_sine, time_scaling, time_weight_scaling, &
                                                max_fit_error, print_warning, unit_nr, prefer_external_backend, &
                                                used_external_backend)

      INTEGER, INTENT(IN)                                :: num_points
      REAL(KIND=dp), INTENT(IN)                          :: energy_min, energy_max, regularization
      INTEGER, INTENT(IN)                                :: num_points_per_magnitude
      TYPE(time_frequency_grid_type), INTENT(OUT)        :: grid
      LOGICAL, INTENT(IN)                                :: build_frequency, build_time, &
                                                            build_transforms, build_sine
      REAL(KIND=dp), INTENT(IN)                          :: time_scaling, time_weight_scaling
      REAL(KIND=dp), INTENT(OUT)                         :: max_fit_error
      LOGICAL, INTENT(IN)                                :: print_warning
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: prefer_external_backend
      LOGICAL, INTENT(OUT), OPTIONAL                     :: used_external_backend

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_minimax_time_frequency_grid'

      INTEGER                                            :: external_ierr, handle, ierr
      REAL(KIND=dp)                                      :: e_range, max_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: points_and_weights

      CALL timeset(routineN, handle)

      IF (PRESENT(used_external_backend)) used_external_backend = .FALSE.
      IF (prefer_external_backend) THEN
         CALL greenx_get_minimax_grid(unit_nr, num_points, energy_min, energy_max, regularization, &
                                      grid%imaginary_time, grid%time_weights_at_zero_frequency, &
                                      grid%frequency, grid%frequency_weights, &
                                      grid%cosine_time_to_frequency_weights, &
                                      grid%cosine_frequency_to_time_weights, &
                                      grid%sine_time_to_frequency_weights, external_ierr)
         IF (external_ierr == 0) THEN
            max_fit_error = 0.0_dp
            IF (PRESENT(used_external_backend)) used_external_backend = .TRUE.
            CALL timestop(handle)
            RETURN
         END IF
         CALL time_frequency_grid_release(grid)
      END IF

      CPASSERT(num_points > 0)
      CPASSERT(energy_min > 0.0_dp)
      CPASSERT(energy_max >= energy_min)
      CPASSERT(time_scaling > 0.0_dp)
      CPASSERT(time_weight_scaling > 0.0_dp)
      CPASSERT(.NOT. build_transforms .OR. (build_frequency .AND. build_time))
      CPASSERT(.NOT. build_sine .OR. build_transforms)
      max_error = 0.0_dp

      e_range = energy_max/energy_min
      ALLOCATE (points_and_weights(2*num_points))

      IF (build_frequency) THEN
         IF (num_points <= 20) THEN
            CALL get_rpa_minimax_coeff(num_points, e_range, points_and_weights, ierr, print_warning)
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

      max_fit_error = max_error

      DEALLOCATE (points_and_weights)

      CALL timestop(handle)

   END SUBROUTINE build_minimax_time_frequency_grid

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

      REAL(KIND=dp), PARAMETER                           :: Rc_max = 1.0E+07_dp

      INTEGER                                            :: ierr, iR, num_integ_points
      REAL(KIND=dp)                                      :: max_error, multiplicator, Rc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tau_tj, tj, x_tw
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: weights_cos_tf_t_to_w

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

END MODULE time_frequency_grids
