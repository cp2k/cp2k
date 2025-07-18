!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the Greenx library
!> \par History
!>      07.2025 Refactored from RPA and BSE modules [Frederick Stein]
! **************************************************************************************************
MODULE greenx_interface
   USE kinds, ONLY: dp
   USE cp_log_handling, ONLY: cp_logger_type, &
                              cp_get_default_logger
   USE cp_output_handling, ONLY: cp_print_key_unit_nr, &
                                 cp_print_key_finished_output, &
                                 cp_print_key_generate_filename, &
                                 low_print_level, &
                                 medium_print_level
   USE input_section_types, ONLY: section_vals_type
   USE machine, ONLY: m_flush
   USE physcon, ONLY: evolt
#if defined (__GREENX)
   USE gx_ac, ONLY: create_thiele_pade, &
                    evaluate_thiele_pade_at, &
                    free_params, &
                    params
   USE gx_minimax, ONLY: gx_minimax_grid
#endif

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'greenx_interface'

   PUBLIC :: greenx_refine_pade, greenx_output_polarizability, greenx_refine_ft, greenx_get_minimax_grid

CONTAINS

! **************************************************************************************************
!> \brief Refines Pade approximants using GreenX, skips this step if GreenX is not available
!> \param e_min ...
!> \param e_max ...
!> \param x_eval ...
!> \param number_of_simulation_steps ...
!> \param number_of_pade_points ...
!> \param logger ...
!> \param ft_section ...
!> \param bse_unit ...
!> \param omega_series ...
!> \param ft_full_series ...
! **************************************************************************************************
   SUBROUTINE greenx_refine_pade(e_min, e_max, x_eval, number_of_simulation_steps, number_of_pade_points, &
                                 logger, ft_section, bse_unit, omega_series, ft_full_series)
      REAL(KIND=dp), INTENT(IN) :: e_min, e_max
      COMPLEX(KIND=dp), DIMENSION(:), POINTER :: x_eval
      INTEGER, INTENT(IN) :: number_of_simulation_steps, number_of_pade_points
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                          :: ft_section
      INTEGER, INTENT(IN) :: bse_unit
      REAL(KIND=dp), DIMENSION(number_of_simulation_steps + 2), INTENT(INOUT) :: omega_series
      REAL(KIND=dp), DIMENSION(6, number_of_simulation_steps + 2), INTENT(INOUT) :: ft_full_series
#if defined (__GREENX)
      INTEGER                                             :: i, ft_unit
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE         :: omega_complex, &
                                                             moments_ft_complex
      COMPLEX(kind=dp), DIMENSION(:, :), ALLOCATABLE      :: moments_eval_complex

      ! Report Padé refinement
      IF (bse_unit > 0) WRITE (bse_unit, '(A10,A27,E23.8E3,E20.8E3)') &
         " PADE_FT| ", "Evaluation grid bounds [eV]", e_min, e_max
      ALLOCATE (omega_complex(number_of_simulation_steps + 2))
      ALLOCATE (moments_ft_complex(number_of_simulation_steps + 2))
      ALLOCATE (moments_eval_complex(3, number_of_pade_points))
      omega_complex(:) = CMPLX(omega_series(:), 0.0, kind=dp)
      DO i = 1, 3
         moments_ft_complex(:) = CMPLX(ft_full_series(2*i - 1, :), &
                                       ft_full_series(2*i, :), &
                                       kind=dp)
         ! Copy the fitting parameters
         ! TODO : Optional direct setting of parameters?
         CALL greenx_refine_ft(e_min, e_max, omega_complex, moments_ft_complex, x_eval, moments_eval_complex(i, :))
      END DO
      ! Write into alternative file
      ft_unit = cp_print_key_unit_nr(logger, ft_section, extension="_PADE.dat", &
                                     file_form="FORMATTED", file_position="REWIND")
      IF (ft_unit > 0) THEN
         DO i = 1, number_of_pade_points
            WRITE (ft_unit, '(E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3,E20.8E3)') &
               REAL(x_eval(i)), REAL(moments_eval_complex(1, i)), AIMAG(moments_eval_complex(1, i)), &
               REAL(moments_eval_complex(2, i)), AIMAG(moments_eval_complex(2, i)), &
               REAL(moments_eval_complex(3, i)), AIMAG(moments_eval_complex(3, i))
         END DO
      END IF
      CALL cp_print_key_finished_output(ft_unit, logger, ft_section)
      DEALLOCATE (omega_complex)
      DEALLOCATE (moments_ft_complex)
      DEALLOCATE (moments_eval_complex)
#else
      IF (bse_unit > 0) WRITE (bse_unit, '(A10,A70)') &
         " PADE_FT| ", "GreenX library is not available. Refinement is skipped"
      MARK_USED(e_min)
      MARK_USED(e_max)
      MARK_USED(x_eval)
      MARK_USED(number_of_simulation_steps)
      MARK_USED(number_of_pade_points)
      MARK_USED(logger)
      MARK_USED(ft_section)
      MARK_USED(omega_series)
      MARK_USED(ft_full_series)
#endif
   END SUBROUTINE greenx_refine_pade
! **************************************************************************************************
!> \brief Outputs the isotropic polarizability tensor element alpha _ ij = mu_i(omega)/E_j(omega),
!>        where i and j are provided by the configuration. The tensor element is energy dependent and
!>        has real and imaginary parts
!> \param logger ...
!> \param pol_section ...
!> \param bse_unit ...
!> \param pol_elements ...
!> \param x_eval ...
!> \param polarizability_refined ...
! **************************************************************************************************
   SUBROUTINE greenx_output_polarizability(logger, pol_section, bse_unit, pol_elements, x_eval, polarizability_refined)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                          :: pol_section
      INTEGER, INTENT(IN) :: bse_unit
      INTEGER, DIMENSION(:, :), POINTER                         :: pol_elements
      COMPLEX(KIND=dp), DIMENSION(:), POINTER :: x_eval
      COMPLEX(kind=dp), DIMENSION(:, :), INTENT(IN)     :: polarizability_refined
#if defined(__GREENX)
      INTEGER                                            :: pol_unit, &
                                                            i, k, n_elems

      n_elems = SIZE(pol_elements, 1)
      ! Print out the refined polarizability to a file
      pol_unit = cp_print_key_unit_nr(logger, pol_section, extension="_PADE.dat", &
                                      file_form="FORMATTED", file_position="REWIND")
      ! Printing for both the stdout and separate file
      IF (pol_unit > 0) THEN
         IF (pol_unit == bse_unit) THEN
            ! Print the stdout preline
            WRITE (pol_unit, '(A21)', advance="no") " POLARIZABILITY_PADE|"
         ELSE
            ! Print also the energy in atomic units
            WRITE (pol_unit, '(A1,A19)', advance="no") "#", "omega [a.u.]"
         END IF
         ! Common - print the energy in eV
         WRITE (pol_unit, '(A20)', advance="no") "Energy [eV]"
         ! Print a header for each polarizability element
         DO k = 1, n_elems - 1
            WRITE (pol_unit, '(A16,I2,I2,A16,I2,I2)', advance="no") &
               "Real pol.", pol_elements(k, 1), pol_elements(k, 2), &
               "Imag pol.", pol_elements(k, 1), pol_elements(k, 2)
         END DO
         WRITE (pol_unit, '(A16,I2,I2,A16,I2,I2)') &
            "Real pol.", pol_elements(n_elems, 1), pol_elements(n_elems, 2), &
            "Imag pol.", pol_elements(n_elems, 1), pol_elements(n_elems, 2)
         DO i = 1, SIZE(x_eval)
            IF (pol_unit == bse_unit) THEN
               ! Print the stdout preline
               WRITE (pol_unit, '(A21)', advance="no") " POLARIZABILITY_PADE|"
            ELSE
               ! omega in a.u.
               WRITE (pol_unit, '(E20.8E3)', advance="no") REAL(x_eval(i), kind=dp)
            END IF
            ! Common values
            WRITE (pol_unit, '(E20.8E3)', advance="no") REAL(x_eval(i), kind=dp)*evolt
            DO k = 1, n_elems - 1
               WRITE (pol_unit, '(E20.8E3,E20.8E3)', advance="no") &
                  REAL(polarizability_refined(k, i)), AIMAG(polarizability_refined(k, i))
            END DO
            ! Print the final value and advance
            WRITE (pol_unit, '(E20.8E3,E20.8E3)') &
               REAL(polarizability_refined(n_elems, i)), AIMAG(polarizability_refined(n_elems, i))
         END DO
         CALL cp_print_key_finished_output(pol_unit, logger, pol_section)
      END IF
#else
      MARK_USED(logger)
      MARK_USED(pol_section)
      MARK_USED(bse_unit)
      MARK_USED(pol_elements)
      MARK_USED(x_eval)
      MARK_USED(polarizability_refined)
#endif
   END SUBROUTINE greenx_output_polarizability
! **************************************************************************************************
!> \brief Refines the FT grid using Padé approximants
!> \param fit_e_min ...
!> \param fit_e_max ...
!> \param x_fit Input x-variables
!> \param y_fit Input y-variables
!> \param x_eval Refined x-variables
!> \param y_eval Refined y-variables
!> \param n_pade_opt ...
! **************************************************************************************************
   SUBROUTINE greenx_refine_ft(fit_e_min, fit_e_max, x_fit, y_fit, x_eval, y_eval, n_pade_opt)
      REAL(kind=dp)                                      :: fit_e_min, &
                                                            fit_e_max
      COMPLEX(kind=dp), DIMENSION(:)                     :: x_fit, &
                                                            y_fit, &
                                                            x_eval, &
                                                            y_eval
      INTEGER, OPTIONAL                                  :: n_pade_opt
#if defined (__GREENX)
      INTEGER                                            :: fit_start, &
                                                            fit_end, &
                                                            max_fit, &
                                                            n_fit, &
                                                            n_pade, &
                                                            n_eval, &
                                                            i
      TYPE(params)                                       :: pade_params

      ! Get the sizes from arrays
      max_fit = SIZE(x_fit)
      n_eval = SIZE(x_eval)

      ! Search for the fit start and end indices
      fit_start = -1
      fit_end = -1
      ! Search for the subset of FT points which is within energy limits given by
      ! the input
      ! Do not search when automatic request of highest energy is made
      IF (fit_e_max < 0) fit_end = max_fit
      DO i = 1, max_fit
         IF (fit_start == -1 .AND. REAL(x_fit(i)) >= fit_e_min) fit_start = i
         IF (fit_end == -1 .AND. REAL(x_fit(i)) > fit_e_max) fit_end = i - 1
         IF (fit_start > 0 .AND. fit_end > 0) EXIT
      END DO
      IF (fit_start == -1) fit_start = 1
      IF (fit_end == -1) fit_end = max_fit
      n_fit = fit_end - fit_start + 1

      n_pade = n_fit/2
      IF (PRESENT(n_pade_opt)) n_pade = n_pade_opt

      ! Warn about a large number of Padé parameters
      IF (n_pade > 1000) THEN
         CPWARN("More then 1000 Padé parameters requested - may reduce with FIT_E_MIN/FIT_E_MAX.")
      END IF
      ! TODO : Symmetry mode settable?
      ! Here, we assume that ft corresponds to transform of real trace
      pade_params = create_thiele_pade(n_pade, x_fit(fit_start:fit_end), y_fit(fit_start:fit_end), &
                                       enforce_symmetry="conjugate")

      ! Check whetner the splice is needed or not
      y_eval(1:n_eval) = evaluate_thiele_pade_at(pade_params, x_eval)

      CALL free_params(pade_params)
#else
      ! Mark used
      MARK_USED(fit_e_min)
      MARK_USED(fit_e_max)
      MARK_USED(x_fit)
      MARK_USED(y_fit)
      MARK_USED(x_eval)
      MARK_USED(y_eval)
      MARK_USED(n_pade_opt)
      CPABORT("Calls to GreenX require CP2K to be compiled with support for GreenX.")
#endif
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param unit_nr ...
!> \param num_integ_points ...
!> \param emin ...
!> \param emax ...
!> \param tau_tj ...
!> \param tau_wj ...
!> \param regularization_minimax ...
!> \param tj ...
!> \param wj ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_cos_tf_w_to_t ...
!> \param weights_sin_tf_t_to_w ...
!> \param ierr ...
! **************************************************************************************************
   SUBROUTINE greenx_get_minimax_grid(unit_nr, num_integ_points, emin, emax, &
                                      tau_tj, tau_wj, regularization_minimax, &
                                      tj, wj, weights_cos_tf_t_to_w, &
                                      weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, ierr)

      INTEGER, INTENT(IN)                                :: unit_nr, num_integ_points
      REAL(KIND=dp), INTENT(IN)                          :: emin, emax
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: tau_tj, tau_wj
      REAL(KIND=dp), INTENT(IN)                 :: regularization_minimax
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: tj, wj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t, &
                                                            weights_sin_tf_t_to_w
      INTEGER, INTENT(OUT) :: ierr
#if defined (__GREENX)
      INTEGER :: gi
      REAL(KIND=dp)                                      :: cosft_duality_error_greenx, &
                                                            max_errors_greenx(3)

      CALL gx_minimax_grid(num_integ_points, Emin, Emax, tau_tj, tau_wj, tj, wj, &
                           weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, weights_sin_tf_t_to_w, &
                           max_errors_greenx, cosft_duality_error_greenx, ierr, &
                           bare_cos_sin_weights=.TRUE., &
                           regularization=regularization_minimax)
      ! Factor 4 is hard-coded in the RPA weights in the internal CP2K minimax routines
      wj(:) = wj(:)*4.0_dp
      IF (ierr == 0) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "GREENX MINIMAX_INFO| Number of integration points:", num_integ_points
            WRITE (UNIT=unit_nr, FMT="(T3,A,T61,3F20.7)") &
               "GREENX MINIMAX_INFO| Gap (Emin):", Emin
            WRITE (UNIT=unit_nr, FMT="(T3,A,T61,3F20.7)") &
               "GREENX MINIMAX_INFO| Maximum eigenvalue difference (Emax):", Emax
            WRITE (UNIT=unit_nr, FMT="(T3,A,T61,3F20.4)") &
               "GREENX MINIMAX_INFO| Energy range (Emax/Emin):", Emax/Emin
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
               "GREENX MINIMAX_INFO| Frequency grid (scaled):", "Weights", "Abscissas"
            DO gi = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") wj(gi), tj(gi)
            END DO
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
               "GREENX MINIMAX_INFO| Time grid (scaled):", "Weights", "Abscissas"
            DO gi = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") tau_wj(gi), tau_tj(gi)
            END DO
            CALL m_flush(unit_nr)
         END IF
      ELSE
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75)") &
               "GREENX MINIMAX_INFO| Grid not available, use internal CP2K grids"
            CALL m_flush(unit_nr)
         END IF
         IF (ALLOCATED(tau_tj)) THEN
            DEALLOCATE (tau_tj)
         END IF
         IF (ALLOCATED(tau_wj)) THEN
            DEALLOCATE (tau_wj)
         END IF
         IF (ALLOCATED(tj)) THEN
            DEALLOCATE (tj)
         END IF
         IF (ALLOCATED(wj)) THEN
            DEALLOCATE (wj)
         END IF
         IF (ALLOCATED(weights_cos_tf_t_to_w)) THEN
            DEALLOCATE (weights_cos_tf_t_to_w)
         END IF
         IF (ALLOCATED(weights_cos_tf_w_to_t)) THEN
            DEALLOCATE (weights_cos_tf_w_to_t)
         END IF
         IF (ALLOCATED(weights_sin_tf_t_to_w)) THEN
            DEALLOCATE (weights_sin_tf_t_to_w)
         END IF
      END IF
#else
      ierr = 1
      MARK_USED(unit_nr)
      MARK_USED(num_integ_points)
      MARK_USED(emin)
      MARK_USED(emax)
      MARK_USED(tau_tj)
      MARK_USED(tau_wj)
      MARK_USED(regularization_minimax)
      MARK_USED(tj)
      MARK_USED(wj)
      MARK_USED(weights_cos_tf_t_to_w)
      MARK_USED(weights_cos_tf_w_to_t)
      MARK_USED(weights_sin_tf_t_to_w)
#endif

   END SUBROUTINE greenx_get_minimax_grid

END MODULE greenx_interface
