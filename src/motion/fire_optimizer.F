!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief FIRE driver for the shared geometry and cell optimization environment.
! **************************************************************************************************
MODULE fire_optimizer
   USE bibliography,                    ONLY: Bitzek2006,&
                                              Guenole2020,&
                                              cite_reference
   USE cell_manifold_utils,             ONLY: cell_manifold_chart_eval,&
                                              info_ok,&
                                              n_slot
   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iterate,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE fire_minimizer,                  ONLY: fire_accept_step,&
                                              fire_init,&
                                              fire_parameters_type,&
                                              fire_prepare_step,&
                                              fire_state_type,&
                                              fire_step
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE gopt_f_methods,                  ONLY: cp_eval_at,&
                                              gopt_f_ii,&
                                              gopt_f_io,&
                                              gopt_f_io_finalize,&
                                              gopt_f_io_init,&
                                              gopt_f_project_coord,&
                                              print_geo_opt_header,&
                                              print_geo_opt_nc
   USE gopt_f_types,                    ONLY: gopt_f_type
   USE gopt_param_types,                ONLY: gopt_param_type
   USE input_constants,                 ONLY: default_cell_method_id,&
                                              default_fire2_method_id,&
                                              default_minimization_method_id
   USE input_cp2k_fire,                 ONLY: read_fire_section
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fire_optimizer'

   PUBLIC :: geoopt_fire

CONTAINS

! **************************************************************************************************
!> \brief Optimize particle coordinates and, for CELL_OPT, reduced cell coordinates.
!> \param force_env force environment
!> \param gopt_param convergence and iteration controls
!> \param globenv global environment
!> \param geo_section GEO_OPT or CELL_OPT input
!> \param gopt_env optimization environment
!> \param x current optimizer coordinates
! **************************************************************************************************
   SUBROUTINE geoopt_fire(force_env, gopt_param, globenv, geo_section, gopt_env, x)
      TYPE(force_env_type), POINTER                      :: force_env
      TYPE(gopt_param_type), POINTER                     :: gopt_param
      TYPE(global_environment_type), POINTER             :: globenv
      TYPE(section_vals_type), POINTER                   :: geo_section
      TYPE(gopt_f_type), POINTER                         :: gopt_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: x

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'geoopt_fire'

      CHARACTER(LEN=5)                                   :: label
      INTEGER                                            :: handle, iter_start, its, ncart, output_unit
      LOGICAL                                            :: converged, fire2, should_stop, stalled
      REAL(KIND=dp)                                      :: cell_mass, energy, eold, emin, start_time
      REAL(KIND=dp), ALLOCATABLE                         :: old_x(:), base_x(:), step(:)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: gradient
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(fire_parameters_type)                         :: params
      TYPE(fire_state_type)                              :: state
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: fire_section, root_section

      CALL timeset(routineN, handle)
      IF (gopt_env%type_id /= default_minimization_method_id .AND. &
          gopt_env%type_id /= default_cell_method_id) THEN
         CPABORT("FIRE requires a geometry or cell minimization.")
      END IF
      logger => cp_get_default_logger()
      root_section => force_env%root_section
      CALL force_env_get(force_env, para_env=para_env)
      output_unit = cp_print_key_unit_nr(logger, geo_section, "PRINT%PROGRAM_RUN_INFO", extension=".geoLog")

      fire2 = gopt_param%method_id == default_fire2_method_id
      label = MERGE("FIRE2", "FIRE ", fire2)
      CALL cite_reference(Bitzek2006)
      IF (fire2) CALL cite_reference(Guenole2020)
      CALL print_geo_opt_header(gopt_env, output_unit, TRIM(label))
      fire_section => section_vals_get_subs_vals(geo_section, "FIRE")
      CALL read_fire_section(fire_section, params)

      ncart = SIZE(x)
      cell_mass = 1.0_dp
      IF (gopt_env%type_id == default_cell_method_id) THEN
         ncart = ncart - gopt_env%cell_env%manifold%ndof
         CALL section_vals_val_get(fire_section, "CELL_MASS_FACTOR", r_val=cell_mass)
         IF (.NOT. (cell_mass >= 0.0_dp)) THEN
            CPABORT("FIRE: CELL_MASS_FACTOR must be nonnegative.")
         END IF
         IF (cell_mass == 0.0_dp) cell_mass = REAL(MAX(1, ncart/3), dp)**(1.0_dp/3.0_dp)
      END IF
      CALL fire_init(state, params, SIZE(x), ncart, fire2, cell_mass)

      ALLOCATE (gradient(SIZE(x)), old_x(SIZE(x)), base_x(SIZE(x)), step(SIZE(x)))
      CALL cp_iterate(logger%iter_info, increment=0, iter_nr_out=iter_start)
      its = iter_start
      converged = .FALSE.
      start_time = m_walltime()
      CALL evaluate()
      CALL gopt_f_io_init(gopt_env, output_unit, energy, label, its, m_walltime() - start_time)
      eold = energy
      emin = energy

      DO WHILE (its < gopt_param%max_iter)
         start_time = m_walltime()
         old_x(:) = x
         CALL fire_prepare_step(state, -gradient, step, stalled)
         IF (stalled) THEN
            IF (output_unit > 0) WRITE (output_unit, '(T2,A)') &
               "FIRE2 stopped after repeated velocity resets without convergence."
            EXIT
         END IF

         its = its + 1
         CALL cp_iterate(logger%iter_info, last=(its == gopt_param%max_iter), iter_nr=its)
         CALL section_vals_val_set(geo_section, "STEP_START_VAL", i_val=its)
         CALL gopt_f_ii(its, output_unit)
         IF (ANY(step /= 0.0_dp)) THEN
            CALL set_trial(old_x, step)
            IF (ANY(x /= old_x)) CALL evaluate()
         END IF

         base_x(:) = x
         CALL fire_step(state, -gradient, step)
         CALL set_trial(base_x, step)
         IF (ANY(x /= base_x)) CALL evaluate()
         CALL fire_accept_step(state, x - base_x)

         CALL gopt_f_io(gopt_env, force_env, root_section, its, energy, output_unit, eold, emin, &
                        label, gopt_param, SIZE(x), x - old_x, gradient, converged, &
                        used_time=m_walltime() - start_time)
         eold = energy
         emin = MIN(emin, energy)
         CALL external_control(should_stop, "GEO", globenv=globenv)
         IF (converged .OR. should_stop) EXIT
      END DO

      IF (.NOT. converged .AND. its >= gopt_param%max_iter) CALL print_geo_opt_nc(gopt_env, output_unit)
      CALL cp_iterate(logger%iter_info, last=.TRUE., increment=0)
      CALL gopt_f_io_finalize(gopt_env, force_env, x, converged, its, root_section, para_env, 0, output_unit)
      CALL cp_print_key_finished_output(output_unit, logger, geo_section, "PRINT%PROGRAM_RUN_INFO")
      DEALLOCATE (gradient)
      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief Evaluate and distribute the projected gradient at the current coordinates.
! **************************************************************************************************
      SUBROUTINE evaluate()
         CALL gopt_f_project_coord(gopt_env, x)
         CALL cp_eval_at(gopt_env, x, energy, gradient, 0, para_env)
         CALL para_env%bcast(x, 0)
         CALL para_env%bcast(gradient, 0)
      END SUBROUTINE evaluate

! **************************************************************************************************
!> \brief Shorten a step outside the cell chart before evaluating forces.
!> \param origin coordinates before the trial
!> \param displacement proposed displacement
! **************************************************************************************************
      SUBROUTINE set_trial(origin, displacement)
      REAL(KIND=dp), INTENT(IN)                          :: origin(:), displacement(:)

      INTEGER                                            :: attempt, info
      REAL(KIND=dp)                                      :: fraction, q(n_slot)

         x = origin + displacement
         IF (SIZE(x) == ncart) RETURN
         fraction = 1.0_dp
         DO attempt = 1, 20
            CALL cell_manifold_chart_eval(gopt_env%cell_env%manifold, x(ncart + 1:), q, info)
            IF (info == info_ok) RETURN
            fraction = 0.5_dp*fraction
            x = origin + fraction*displacement
         END DO
         CPABORT("FIRE: could not find a valid cell step. Reduce TIMESTEP or MAX_STEP.")
      END SUBROUTINE set_trial
   END SUBROUTINE geoopt_fire

END MODULE fire_optimizer
