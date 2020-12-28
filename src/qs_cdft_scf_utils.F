!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines for performing a constrained DFT SCF run with Quickstep.
!> \par History
!>      - Separated some routines from qs_scf (03.2018) [Nico Holmberg]
!> \author Nico Holmberg (03.2018)
! **************************************************************************************************
MODULE qs_cdft_scf_utils
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_files,                        ONLY: close_file,&
                                              get_unit_number,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_logger_create,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_constants,                 ONLY: &
        broyden_type_1, broyden_type_1_explicit, broyden_type_1_explicit_ls, broyden_type_1_ls, &
        broyden_type_2, broyden_type_2_explicit, broyden_type_2_explicit_ls, broyden_type_2_ls, &
        jacobian_fd1, jacobian_fd1_backward, jacobian_fd1_central, jacobian_fd2, &
        jacobian_fd2_backward, outer_scf_optimizer_broyden, outer_scf_optimizer_newton, &
        outer_scf_optimizer_newton_ls
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathlib,                         ONLY: invert_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_cdft_scf_utils'

   PUBLIC :: prepare_jacobian_stencil, build_diagonal_jacobian, &
             restart_inverse_jacobian, print_inverse_jacobian, &
             create_tmp_logger, initialize_inverse_jacobian

CONTAINS

! **************************************************************************************************
!> \brief Prepares the finite difference stencil for computing the Jacobian. The constraints
!>        are re-evaluated by perturbing each constraint.
!> \param qs_env the qs_env where to build the Jacobian
!> \param output_unit the output unit number
!> \param nwork the number of perturbations to take in the negative direction
!> \param pwork the number of perturbations to take in the positive direction
!> \param coeff list of coefficients that determine how to sum up the various perturbations
!> \param step_multiplier list of values that determine how large steps to take for each perturbatio
!> \param dh total length of the interval to use for computing the finite difference derivatives
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE prepare_jacobian_stencil(qs_env, output_unit, nwork, pwork, &
                                       coeff, step_multiplier, dh)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: output_unit, nwork, pwork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coeff, step_multiplier, dh

      CHARACTER(len=15)                                  :: fmt_code
      INTEGER                                            :: ivar
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (scf_env, scf_control, dft_control)

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, &
                      scf_control=scf_control, &
                      dft_control=dft_control)

      IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) /= 1 .AND. &
          SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) /= SIZE(scf_env%outer_scf%variables, 1)) &
         CALL cp_abort(__LOCATION__, &
                       cp_to_string(SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step))// &
                       " values passed to keyword JACOBIAN_STEP, expected 1 or "// &
                       cp_to_string(SIZE(scf_env%outer_scf%variables, 1)))

      ALLOCATE (dh(SIZE(scf_env%outer_scf%variables, 1)))
      IF (SIZE(dh) /= SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step)) THEN
         DO ivar = 1, SIZE(dh)
            dh(ivar) = scf_control%outer_scf%cdft_opt_control%jacobian_step(1)
         END DO
      ELSE
         dh(:) = scf_control%outer_scf%cdft_opt_control%jacobian_step
      END IF

      SELECT CASE (scf_control%outer_scf%cdft_opt_control%jacobian_type)
      CASE DEFAULT
         CALL cp_abort(__LOCATION__, &
                       "Unknown Jacobian type: "// &
                       cp_to_string(scf_control%outer_scf%cdft_opt_control%jacobian_type))
      CASE (jacobian_fd1)
         ! f'(x0) = [ f(x0+h) - f(x0) ] / h
         nwork = 0
         pwork = 1
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(nwork) = -1.0_dp
         coeff(pwork) = 1.0_dp
         step_multiplier = 1.0_dp
      CASE (jacobian_fd1_backward)
         ! f'(x0) = [ f(x0) - f(x0-h) ] / h
         nwork = -1
         pwork = 0
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(nwork) = -1.0_dp
         coeff(pwork) = 1.0_dp
         step_multiplier = -1.0_dp
      CASE (jacobian_fd2)
         ! f'(x0) = [ -f(x0+2h) + 4f(x0+h) - 3f(x0) ] / 2h
         nwork = 0
         pwork = 2
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(0) = -3.0_dp
         coeff(1) = 4.0_dp
         coeff(2) = -1.0_dp
         step_multiplier(0) = 0.0_dp
         step_multiplier(1) = 1.0_dp
         step_multiplier(2) = 2.0_dp
         dh(:) = 2.0_dp*dh(:)
      CASE (jacobian_fd2_backward)
         ! f'(x0) = [ 3f(x0) - 4f(x0-h) + f(x0-2h) ] / 2h
         nwork = -2
         pwork = 0
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(0) = 3.0_dp
         coeff(-1) = -4.0_dp
         coeff(-2) = 1.0_dp
         step_multiplier(0) = 0.0_dp
         step_multiplier(-1) = -1.0_dp
         step_multiplier(-2) = -2.0_dp
         dh(:) = 2.0_dp*dh(:)
      CASE (jacobian_fd1_central)
         ! f'(x0) = [ f(x0+h) - f(x0-h) ] / 2h
         nwork = -1
         pwork = 1
         ALLOCATE (coeff(nwork:pwork), step_multiplier(nwork:pwork))
         coeff(0) = 0.0_dp
         coeff(nwork) = -1.0_dp
         coeff(pwork) = 1.0_dp
         step_multiplier(0) = 0.0_dp
         step_multiplier(nwork) = -1.0_dp
         step_multiplier(pwork) = 1.0_dp
         dh(:) = 2.0_dp*dh(:)
      END SELECT
      ! Print some info
      IF (output_unit > 0) THEN
         WRITE (output_unit, FMT="(/,A)") &
            " ================================== JACOBIAN CALCULATION ================================="
         WRITE (output_unit, FMT="(A)") &
            " Evaluating inverse Jacobian using finite differences"
         WRITE (output_unit, '(A,I10,A,I10)') &
            " Energy evaluation: ", dft_control%qs_control%cdft_control%ienergy, &
            ", CDFT SCF iteration: ", scf_env%outer_scf%iter_count
         SELECT CASE (scf_control%outer_scf%cdft_opt_control%jacobian_type)
         CASE (jacobian_fd1)
            WRITE (output_unit, '(A)') " Type               : First order forward difference"
         CASE (jacobian_fd1_backward)
            WRITE (output_unit, '(A)') " Type               : First order backward difference"
         CASE (jacobian_fd2)
            WRITE (output_unit, '(A)') " Type               : Second order forward difference"
         CASE (jacobian_fd2_backward)
            WRITE (output_unit, '(A)') " Type               : Second order backward difference"
         CASE (jacobian_fd1_central)
            WRITE (output_unit, '(A)') " Type               : First order central difference"
         CASE DEFAULT
            CALL cp_abort(__LOCATION__, "Unknown Jacobian type: "// &
                          cp_to_string(scf_control%outer_scf%cdft_opt_control%jacobian_type))
         END SELECT
         IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) == 1) THEN
            WRITE (output_unit, '(A,ES12.4)') " Step size          : ", scf_control%outer_scf%cdft_opt_control%jacobian_step
         ELSE
            WRITE (output_unit, '(A,ES12.4,A)') &
               " Step sizes         : ", scf_control%outer_scf%cdft_opt_control%jacobian_step(1), ' (constraint 1)'
            IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step) < 10) THEN
               fmt_code = '(ES34.4,A,I2,A)'
            ELSE
               fmt_code = '(ES34.4,A,I3,A)'
            END IF
            DO ivar = 2, SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_step)
               WRITE (output_unit, fmt_code) scf_control%outer_scf%cdft_opt_control%jacobian_step(ivar), " (constraint", ivar, ")"
            END DO
         END IF
      END IF

   END SUBROUTINE prepare_jacobian_stencil
! **************************************************************************************************
!> \brief Builds a strictly diagonal inverse Jacobian from MD/SCF history.
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \param used_history flag that determines if history was actually used to prepare the Jacobian
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE build_diagonal_jacobian(qs_env, used_history)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL                                            :: used_history

      INTEGER                                            :: i, ihistory, nvar, outer_scf_ihistory
      LOGICAL                                            :: use_md_history
      REAL(KIND=dp)                                      :: inv_error
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: jacobian
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient_history, inv_jacobian, &
                                                            variable_history
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (scf_control, scf_env)

      CALL get_qs_env(qs_env, scf_env=scf_env, &
                      scf_control=scf_control, &
                      outer_scf_ihistory=outer_scf_ihistory)
      ihistory = scf_env%outer_scf%iter_count

      IF (outer_scf_ihistory .GE. 3 .AND. .NOT. used_history) THEN
         ! First, lets try using the history from previous energy evaluations
         CALL get_qs_env(qs_env, gradient_history=gradient_history, &
                         variable_history=variable_history)
         nvar = SIZE(scf_env%outer_scf%variables, 1)
         use_md_history = .TRUE.
         ! Check that none of the history values are identical in which case we should try something different
         DO i = 1, nvar
            IF (ABS(variable_history(i, 2) - variable_history(i, 1)) .LT. 1.0E-12_dp) &
               use_md_history = .FALSE.
         END DO
         IF (use_md_history) THEN
            ALLOCATE (jacobian(nvar, nvar))
            DO i = 1, nvar
               jacobian(i, i) = (gradient_history(i, 2) - gradient_history(i, 1))/ &
                                (variable_history(i, 2) - variable_history(i, 1))
            END DO
            IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
               ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
            inv_jacobian => scf_env%outer_scf%inv_jacobian
            CALL invert_matrix(jacobian, inv_jacobian, inv_error)
            DEALLOCATE (jacobian)
            ! Mark that an inverse Jacobian was just built and the next outer_loop_optimize should not perform
            ! a Broyden update of it
            scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
            ! Mark that the MD history has been used and should not be reused anymore on this energy evaluation
            used_history = .TRUE.
         END IF
      END IF
      IF (ihistory .GE. 2 .AND. .NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
         ! Next, try history from current SCF procedure
         nvar = SIZE(scf_env%outer_scf%variables, 1)
         IF (SIZE(scf_env%outer_scf%gradient, 2) .LT. 3) &
            CALL cp_abort(__LOCATION__, &
                          "Keyword EXTRAPOLATION_ORDER in section OUTER_SCF must be greater than or equal "// &
                          "to 3 for optimizers that build the Jacobian from SCF history.")
         ALLOCATE (jacobian(nvar, nvar))
         DO i = 1, nvar
            jacobian(i, i) = (scf_env%outer_scf%gradient(i, ihistory) - scf_env%outer_scf%gradient(i, ihistory - 1))/ &
                             (scf_env%outer_scf%variables(i, ihistory) - scf_env%outer_scf%variables(i, ihistory - 1))
         END DO
         IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
            ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
         inv_jacobian => scf_env%outer_scf%inv_jacobian
         CALL invert_matrix(jacobian, inv_jacobian, inv_error)
         DEALLOCATE (jacobian)
         scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
      ELSE
         ! No history => will fall back to SD optimizer in outer_loop_optimize
      END IF

   END SUBROUTINE build_diagonal_jacobian

! **************************************************************************************************
!> \brief Restarts the finite difference inverse Jacobian.
!> \param qs_env the qs_environment_type where to compute the Jacobian
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE restart_inverse_jacobian(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: i, iwork, j, nvar
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: inv_jacobian
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (scf_env, scf_control)
      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, scf_env=scf_env, scf_control=scf_control)

      CPASSERT(ASSOCIATED(scf_control%outer_scf%cdft_opt_control%jacobian_vector))
      nvar = SIZE(scf_env%outer_scf%variables, 1)
      IF (SIZE(scf_control%outer_scf%cdft_opt_control%jacobian_vector) /= nvar**2) &
         CALL cp_abort(__LOCATION__, &
                       "Too many or too few values defined for restarting inverse Jacobian.")
      IF (.NOT. ASSOCIATED(scf_env%outer_scf%inv_jacobian)) &
         ALLOCATE (scf_env%outer_scf%inv_jacobian(nvar, nvar))
      inv_jacobian => scf_env%outer_scf%inv_jacobian
      iwork = 1
      DO i = 1, nvar
         DO j = 1, nvar
            inv_jacobian(i, j) = scf_control%outer_scf%cdft_opt_control%jacobian_vector(iwork)
            iwork = iwork + 1
         END DO
      END DO
      DEALLOCATE (scf_control%outer_scf%cdft_opt_control%jacobian_vector)
      scf_control%outer_scf%cdft_opt_control%jacobian_restart = .FALSE.
      scf_control%outer_scf%cdft_opt_control%broyden_update = .FALSE.
      scf_env%outer_scf%deallocate_jacobian = .FALSE.

   END SUBROUTINE restart_inverse_jacobian

! **************************************************************************************************
!> \brief Prints the finite difference inverse Jacobian to file
!> \param logger the default IO logger
!> \param inv_jacobian the inverse Jacobian matrix
!> \param iter_count the iteration number
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE print_inverse_jacobian(logger, inv_jacobian, iter_count)
      TYPE(cp_logger_type), POINTER                      :: logger
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: inv_jacobian
      INTEGER                                            :: iter_count

      CHARACTER(len=default_path_length)                 :: project_name
      INTEGER                                            :: i, j, lp, nvar, output_unit

      nvar = SIZE(inv_jacobian, 1)
      output_unit = get_unit_number()
      project_name = logger%iter_info%project_name
      lp = LEN_TRIM(project_name)
      project_name(lp + 1:LEN(project_name)) = ".inverseJacobian"
      CALL open_file(file_name=project_name, file_status="UNKNOWN", &
                     file_action="WRITE", file_position="APPEND", &
                     unit_number=output_unit)
      WRITE (output_unit, FMT="(/,A)") "Inverse Jacobian matrix in row major order"
      WRITE (output_unit, FMT="(A,I10)") "Iteration: ", iter_count
      DO i = 1, nvar
         DO j = 1, nvar
            WRITE (output_unit, *) inv_jacobian(i, j)
         END DO
      END DO
      CALL close_file(unit_number=output_unit)

   END SUBROUTINE print_inverse_jacobian

! **************************************************************************************************
!> \brief Creates a temporary logger for redirecting output to a new file
!> \param para_env the para_env
!> \param project_name the project basename
!> \param suffix the suffix
!> \param output_unit the default unit number for the newly created temporary logger
!> \param tmp_logger pointer to the newly created temporary logger
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE create_tmp_logger(para_env, project_name, suffix, output_unit, tmp_logger)
      TYPE(cp_para_env_type), POINTER                    :: para_env
      CHARACTER(len=*)                                   :: project_name, suffix
      INTEGER, INTENT(OUT)                               :: output_unit
      TYPE(cp_logger_type), INTENT(OUT), POINTER         :: tmp_logger

      INTEGER                                            :: lp

      IF (para_env%ionode) THEN
         lp = LEN_TRIM(project_name)
         project_name(lp + 1:LEN(project_name)) = suffix
         CALL open_file(file_name=project_name, file_status="UNKNOWN", &
                        file_action="WRITE", file_position="APPEND", &
                        unit_number=output_unit)
      ELSE
         output_unit = -1
      END IF
      CALL cp_logger_create(tmp_logger, &
                            para_env=para_env, &
                            default_global_unit_nr=output_unit, &
                            close_global_unit_on_dealloc=.FALSE.)

   END SUBROUTINE create_tmp_logger

! **************************************************************************************************
!> \brief Checks if the inverse Jacobian should be calculated and initializes the calculation
!> \param scf_control         the scf_control that holds the Jacobian settings
!> \param scf_env             the scf_env that holds the CDFT iteration information
!> \param explicit_jacobian   flag that determines if the finite difference Jacobian is needed
!> \param should_build        flag that determines if the Jacobian should be built
!> \param used_history        flag that determines if SCF history has been used to build a Jacobian
!> \par History
!>      03.2018 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE initialize_inverse_jacobian(scf_control, scf_env, explicit_jacobian, &
                                          should_build, used_history)
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      LOGICAL                                            :: explicit_jacobian, should_build, &
                                                            used_history

      CPASSERT(ASSOCIATED(scf_control))
      CPASSERT(ASSOCIATED(scf_env))

      SELECT CASE (scf_control%outer_scf%optimizer)
      CASE DEFAULT
         CPABORT("Noncompatible optimizer requested.")
      CASE (outer_scf_optimizer_newton, outer_scf_optimizer_newton_ls)
         CPASSERT(ASSOCIATED(scf_control%outer_scf%cdft_opt_control))
         scf_control%outer_scf%cdft_opt_control%build_jacobian = .TRUE.
         explicit_jacobian = .TRUE.
      CASE (outer_scf_optimizer_broyden)
         CPASSERT(ASSOCIATED(scf_control%outer_scf%cdft_opt_control))
         SELECT CASE (scf_control%outer_scf%cdft_opt_control%broyden_type)
         CASE (broyden_type_1, broyden_type_2, broyden_type_1_ls, broyden_type_2_ls)
            scf_control%outer_scf%cdft_opt_control%build_jacobian = .TRUE.
            explicit_jacobian = .FALSE.
         CASE (broyden_type_1_explicit, broyden_type_2_explicit, broyden_type_1_explicit_ls, broyden_type_2_explicit_ls)
            scf_control%outer_scf%cdft_opt_control%build_jacobian = .TRUE.
            explicit_jacobian = .TRUE.
         END SELECT
      END SELECT
      IF (scf_control%outer_scf%cdft_opt_control%build_jacobian) THEN
         ! Reset counter
         IF (scf_env%outer_scf%iter_count == 1) scf_control%outer_scf%cdft_opt_control%ijacobian(1) = 0
         ! Check if an old Jacobian can be reused avoiding a rebuild
         IF (ASSOCIATED(scf_env%outer_scf%inv_jacobian)) THEN
            ! Rebuild if number of previous energy evaluations exceeds limit
            IF (scf_control%outer_scf%cdft_opt_control%ijacobian(2) .GE. &
                scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) .AND. .NOT. used_history .AND. &
                scf_control%outer_scf%cdft_opt_control%jacobian_freq(2) > 0) THEN
               should_build = .TRUE.
               ! Zero the corresponding counters
               scf_control%outer_scf%cdft_opt_control%ijacobian(:) = 0
               ! Rebuild if number of previous SCF iterations exceeds limit (on this energy eval)
            ELSE IF (scf_control%outer_scf%cdft_opt_control%ijacobian(1) .GE. &
                     scf_control%outer_scf%cdft_opt_control%jacobian_freq(1) .AND. &
                     scf_control%outer_scf%cdft_opt_control%jacobian_freq(1) > 0) THEN
               should_build = .TRUE.
               ! Zero the corresponding counter
               scf_control%outer_scf%cdft_opt_control%ijacobian(1) = 0
            END IF
            IF (should_build) DEALLOCATE (scf_env%outer_scf%inv_jacobian)
         ELSE
            should_build = .TRUE.
            ! Zero the counter
            scf_control%outer_scf%cdft_opt_control%ijacobian(:) = 0
         END IF
      END IF

   END SUBROUTINE initialize_inverse_jacobian

END MODULE qs_cdft_scf_utils
