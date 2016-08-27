!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Methods for sampling helium variables
!> \author Lukasz Walewski
!> \date   2009-06-10
! **************************************************************************************************
MODULE helium_sampling

   USE cp_external_control,             ONLY: external_control
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                              cp_iterate,&
                                              cp_rm_iter_level
   USE global_types,                    ONLY: global_environment_type
   USE helium_common,                   ONLY: &
        helium_boxmean_3d, helium_calc_plength, helium_calc_rdf, helium_calc_rho, helium_com, &
        helium_eval_expansion, helium_pbc, helium_rotate, helium_set_rdf_coord_system, &
        helium_spline, helium_total_moment_of_inertia, helium_total_projected_area, &
        helium_total_winding_number, helium_update_transition_matrix
   USE helium_interactions,             ONLY: helium_bead_solute_e_f,&
                                              helium_calc_energy,&
                                              helium_solute_e_f
   USE helium_io,                       ONLY: &
        helium_print_accepts, helium_print_action, helium_print_coordinates, helium_print_energy, &
        helium_print_force, helium_print_perm, helium_print_plength, helium_print_rdf, &
        helium_print_rho, helium_print_vector, helium_write_line
   USE helium_types,                    ONLY: e_id_total,&
                                              helium_solvent_p_type,&
                                              helium_solvent_type
   USE input_constants,                 ONLY: &
        helium_forces_average, helium_forces_last, helium_mdist_exponential, &
        helium_mdist_gaussian, helium_mdist_linear, helium_mdist_quadratic, helium_mdist_singlev, &
        helium_mdist_uniform, helium_sampling_ceperley
   USE input_cp2k_restarts,             ONLY: write_restart
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
   USE parallel_rng_types,              ONLY: next_random_number
   USE physcon,                         ONLY: angstrom
   USE pint_public,                     ONLY: pint_com_pos
   USE pint_types,                      ONLY: pint_env_type
   USE splines_types,                   ONLY: spline_data_p_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'helium_sampling'

   PUBLIC :: helium_do_run
   PUBLIC :: helium_sample
   PUBLIC :: helium_step

CONTAINS

! ***************************************************************************
!> \brief  Performs MC simulation for helium (only)
!> \param helium_env ...
!> \param globenv ...
!> \date   2009-07-14
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
!> \note   This routine gets called only when HELIUM_ONLY is set to .TRUE.,
!>         so do not put any property calculations here!
! **************************************************************************************************
   SUBROUTINE helium_do_run(helium_env, globenv)
      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env
      TYPE(global_environment_type), POINTER             :: globenv

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_do_run', routineP = moduleN//':'//routineN

      INTEGER                                            :: k, num_steps, step, tot_steps
      LOGICAL                                            :: should_stop
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(pint_env_type), POINTER                       :: pint_env

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (pint_env)
      num_steps = 0

      IF (ASSOCIATED(helium_env)) THEN
         DO k = 1, SIZE(helium_env)

            NULLIFY (helium_env(k)%helium%logger)
            helium_env(k)%helium%logger => cp_get_default_logger()

            ! create iteration level
            ! The iteration level 'MD' can be misleading since this is a pure MC
            ! simulation, however the helium code is meant to be run mainly from PINT.
            ! We keep it like that for compatibility with PINT+HELIUM simulations.
            IF (k .EQ. 1) THEN
               CALL cp_add_iter_level(helium_env(k)%helium%logger%iter_info, "MD")
               CALL cp_iterate(helium_env(k)%helium%logger%iter_info, &
                               iter_nr=helium_env(k)%helium%first_step)
            END IF
            tot_steps = helium_env(k)%helium%first_step
            num_steps = helium_env(k)%helium%num_steps

            ! set the properties accumulated over the whole MC process to 0
            helium_env(k)%helium%proarea%accu(:) = 0.0_dp
            helium_env(k)%helium%prarea2%accu(:) = 0.0_dp
            helium_env(k)%helium%wnmber2%accu(:) = 0.0_dp
            helium_env(k)%helium%mominer%accu(:) = 0.0_dp
            IF (helium_env(k)%helium%rho_present) THEN
               helium_env(k)%helium%rho_accu(:, :, :, :) = 0.0_dp
            END IF
            IF (helium_env(k)%helium%rdf_present) THEN
               helium_env(k)%helium%rdf_accu(:, :, :) = 0.0_dp
            END IF
         END DO
      END IF

      ! Distribute steps for processors without helium_env
      CALL mp_bcast(num_steps, logger%para_env%source, logger%para_env%group)
      CALL mp_bcast(tot_steps, logger%para_env%source, logger%para_env%group)

      DO step = 1, num_steps

         tot_steps = tot_steps+1

         IF (ASSOCIATED(helium_env)) THEN
            DO k = 1, SIZE(helium_env)
               IF (k .EQ. 1) THEN
                  CALL cp_iterate(helium_env(k)%helium%logger%iter_info, &
                                  last=(step == helium_env(k)%helium%num_steps), &
                                  iter_nr=tot_steps)
               END IF
               helium_env(k)%helium%current_step = tot_steps
            END DO
         END IF

         CALL helium_step(helium_env, pint_env)

         ! call write_restart here to avoid interference with PINT write_restart
         IF (ASSOCIATED(helium_env)) THEN
            CALL write_restart(root_section=helium_env(1)%helium%input, helium_env=helium_env)
         END IF

         ! exit from the main loop if soft exit has been requested
         CALL external_control(should_stop, "MD", globenv=globenv)
         IF (should_stop) EXIT

      END DO

      ! remove iteration level
      IF (ASSOCIATED(helium_env)) THEN
         CALL cp_rm_iter_level(helium_env(1)%helium%logger%iter_info, "MD")
      END IF

      RETURN
   END SUBROUTINE helium_do_run

! ***************************************************************************
!> \brief  Sample the helium phase space
!> \param helium_env ...
!> \param pint_env ...
!> \date   2009-10-27
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
!> \note   Samples helium variable space according to multilevel Metropolis
!>         MC scheme, calculates the forces exerted by helium solvent on the
!>         solute and stores them in helium%force_avrg array. The forces are
!>         averaged over outer MC loop.
!> \note   The implicit assumption is that we simulate solute _with_ solvent
!>         most of the time, so for performance reasons I do not check that
!>         everywhere I should. This leads to some redundancy in the case of
!>         helium-only calculations.
! **************************************************************************************************
   SUBROUTINE helium_sample(helium_env, pint_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env
      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_sample', routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: msg_str
      INTEGER                                            :: i, irot, iweight, k, nslices, nsteps, &
                                                            num_env, offset, sel_mp_source
      REAL(KIND=dp)                                      :: inv_num_env, inv_xn, rnd, rtmp, rweight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work_2d
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      DO k = 1, SIZE(helium_env)

         IF (helium_env(k)%helium%solute_present) THEN
            helium_env(k)%helium%force_avrg(:, :) = 0.0_dp
            helium_env(k)%helium%current_step = pint_env%iter
            helium_env(k)%helium%origin = pint_com_pos(pint_env)
         END IF

         helium_env(k)%helium%energy_avrg(:) = 0.0_dp
         helium_env(k)%helium%plength_avrg(:) = 0.0_dp
         helium_env(k)%helium%num_accepted(:, :) = 0.0_dp

         ! helium parallelization: each processor gets different RN stream and
         ! runs independent helium simulation, the properties and forces are
         ! averaged over parallel helium environments once per step.
         inv_xn = 0.0_dp
         SELECT CASE (helium_env (k)%helium%sampling_method)
         CASE (helium_sampling_ceperley)
            ! helium sampling (outer MC loop)
            DO irot = 1, helium_env(k)%helium%iter_rot

               ! rotate helium beads in imaginary time at random, store current
               ! 'rotation state' in helium%relrot wich is within (0, helium%beads-1)
               ! (this is needed to sample different fragments of the permutation
               ! paths in try_permutations)
               rnd = next_random_number(helium_env(k)%helium%rng_stream_uniform)
               nslices = INT(rnd*helium_env(k)%helium%beads)
               CALL helium_rotate(helium_env(k)%helium, nslices)

               CALL helium_try_permutations(helium_env(k)%helium, pint_env)

               ! calculate & accumulate instantaneous properties for averaging
               IF (helium_env(k)%helium%solute_present) THEN
                  IF (helium_env(k)%helium%get_helium_forces == helium_forces_average) THEN
                     CALL helium_solute_e_f(pint_env, helium_env(k)%helium, rtmp)
                     helium_env(k)%helium%force_avrg(:, :) = helium_env(k)%helium%force_avrg(:, :)+ &
                                                             helium_env(k)%helium%force_inst(:, :)
                  END IF
               END IF
               CALL helium_calc_energy(helium_env(k)%helium, pint_env)

               helium_env(k)%helium%energy_avrg(:) = helium_env(k)%helium%energy_avrg(:)+helium_env(k)%helium%energy_inst(:)
               CALL helium_calc_plength(helium_env(k)%helium)
               helium_env(k)%helium%plength_avrg(:) = helium_env(k)%helium%plength_avrg(:)+helium_env(k)%helium%plength_inst(:)

               ! instantaneous force output according to HELIUM%PRINT%FORCES_INST
               ! Warning: file I/O here may cost A LOT of cpu time!
               ! switched off here to save cpu
               !CALL helium_print_force_inst( helium_env(k)%helium, error )

            END DO ! outer MC loop

            ! If only the last forces are taken, calculate them for the last outer MC loop configuration
            IF ((helium_env(k)%helium%solute_present) .AND. (helium_env(k)%helium%get_helium_forces == helium_forces_last)) THEN
               CALL helium_solute_e_f(pint_env, helium_env(k)%helium, rtmp)
               helium_env(k)%helium%force_avrg(:, :) = helium_env(k)%helium%force_inst(:, :)
            END IF

            ! restore the original alignment of beads in imaginary time
            ! (this is useful to make a single bead's position easy to follow
            ! in the trajectory, otherwise it's index would change every step)
            CALL helium_rotate(helium_env(k)%helium, -helium_env(k)%helium%relrot)

            inv_xn = 1.0_dp/REAL(helium_env(k)%helium%iter_rot, dp)

         CASE DEFAULT
            WRITE (msg_str, *) helium_env(k)%helium%sampling_method
            msg_str = "Sampling method ("//TRIM(ADJUSTL(msg_str))//") not supported."
            CPABORT(msg_str)

         END SELECT

         ! actually average the properties over the outer MC loop
         IF ((helium_env(k)%helium%solute_present) .AND. (helium_env(k)%helium%get_helium_forces == helium_forces_average)) THEN
            helium_env(k)%helium%force_avrg(:, :) = helium_env(k)%helium%force_avrg(:, :)*inv_xn
         END IF
         helium_env(k)%helium%energy_avrg(:) = helium_env(k)%helium%energy_avrg(:)*inv_xn
         helium_env(k)%helium%plength_avrg(:) = helium_env(k)%helium%plength_avrg(:)*inv_xn

         ! instantaneous properties
         helium_env(k)%helium%proarea%inst(:) = helium_total_projected_area(helium_env(k)%helium)
         helium_env(k)%helium%prarea2%inst(:) = helium_env(k)%helium%proarea%inst(:)*helium_env(k)%helium%proarea%inst(:)
         helium_env(k)%helium%wnumber%inst(:) = helium_total_winding_number(helium_env(k)%helium)
         helium_env(k)%helium%wnmber2%inst(:) = helium_env(k)%helium%wnumber%inst(:)*helium_env(k)%helium%wnumber%inst(:)
         helium_env(k)%helium%mominer%inst(:) = helium_total_moment_of_inertia(helium_env(k)%helium)

         ! properties accumulated over the whole MC process
         helium_env(k)%helium%proarea%accu(:) = helium_env(k)%helium%proarea%accu(:)+helium_env(k)%helium%proarea%inst(:)
         helium_env(k)%helium%prarea2%accu(:) = helium_env(k)%helium%prarea2%accu(:)+helium_env(k)%helium%prarea2%inst(:)
         helium_env(k)%helium%wnmber2%accu(:) = helium_env(k)%helium%wnmber2%accu(:)+helium_env(k)%helium%wnmber2%inst(:)
         helium_env(k)%helium%mominer%accu(:) = helium_env(k)%helium%mominer%accu(:)+helium_env(k)%helium%mominer%inst(:)
         IF (helium_env(k)%helium%rho_present) THEN
            CALL helium_calc_rho(helium_env(k)%helium)
            helium_env(k)%helium%rho_accu(:, :, :, :) = helium_env(k)%helium%rho_accu(:, :, :, :)+ &
                                                        helium_env(k)%helium%rho_inst(:, :, :, :)
         END IF
         IF (helium_env(k)%helium%rdf_present) THEN
            CALL helium_set_rdf_coord_system(helium_env(k)%helium, pint_env)
            CALL helium_calc_rdf(helium_env(k)%helium, helium_env(k)%helium%rdf_centers)
            helium_env(k)%helium%rdf_accu(:, :, :) = helium_env(k)%helium%rdf_accu(:, :, :)+helium_env(k)%helium%rdf_inst(:, :, :)
         END IF

         ! running averages (restart-aware)
         nsteps = helium_env(k)%helium%current_step-helium_env(k)%helium%first_step
         iweight = helium_env(k)%helium%averages_iweight
         rweight = REAL(iweight, dp)
         rtmp = 1.0_dp/(REAL(MAX(1, nsteps+iweight), dp))
         helium_env(k)%helium%proarea%ravr(:) = (helium_env(k)%helium%proarea%accu(:)+ &
                                                 rweight*helium_env(k)%helium%proarea%rstr(:))*rtmp
         helium_env(k)%helium%prarea2%ravr(:) = (helium_env(k)%helium%prarea2%accu(:)+ &
                                                 rweight*helium_env(k)%helium%prarea2%rstr(:))*rtmp
         helium_env(k)%helium%wnmber2%ravr(:) = (helium_env(k)%helium%wnmber2%accu(:)+ &
                                                 rweight*helium_env(k)%helium%wnmber2%rstr(:))*rtmp
         helium_env(k)%helium%mominer%ravr(:) = (helium_env(k)%helium%mominer%accu(:)+ &
                                                 rweight*helium_env(k)%helium%mominer%rstr(:))*rtmp

      END DO

      ! average over helium environments sitting at different processors
!TODO the following involves message passing, maybe move it to i/o routines?
      num_env = helium_env(1)%helium%num_env
      inv_num_env = 1.0_dp/REAL(num_env, dp)

      !energy_avrg:
      DO k = 2, SIZE(helium_env)
         helium_env(1)%helium%energy_avrg(:) = helium_env(1)%helium%energy_avrg(:)+ &
                                               helium_env(k)%helium%energy_avrg(:)
      END DO
      CALL mp_sum(helium_env(1)%helium%energy_avrg(:), helium_env(1)%comm)
      helium_env(1)%helium%energy_avrg(:) = helium_env(1)%helium%energy_avrg(:)*inv_num_env
      DO k = 2, SIZE(helium_env)
         helium_env(k)%helium%energy_avrg(:) = helium_env(1)%helium%energy_avrg(:)
      END DO

      !plength_avrg:
      DO k = 2, SIZE(helium_env)
         helium_env(1)%helium%plength_avrg(:) = helium_env(1)%helium%plength_avrg(:)+ &
                                                helium_env(k)%helium%plength_avrg(:)
      END DO
      CALL mp_sum(helium_env(1)%helium%plength_avrg(:), helium_env(1)%comm)
      helium_env(1)%helium%plength_avrg(:) = helium_env(1)%helium%plength_avrg(:)*inv_num_env
      DO k = 2, SIZE(helium_env)
         helium_env(k)%helium%plength_avrg(:) = helium_env(1)%helium%plength_avrg(:)
      END DO

      !num_accepted:
      DO k = 2, SIZE(helium_env)
         helium_env(1)%helium%num_accepted(:, :) = helium_env(1)%helium%num_accepted(:, :)+ &
                                                   helium_env(k)%helium%num_accepted(:, :)
      END DO
      CALL mp_sum(helium_env(1)%helium%num_accepted(:, :), helium_env(1)%comm)
      helium_env(1)%helium%num_accepted(:, :) = helium_env(1)%helium%num_accepted(:, :)*inv_num_env
      DO k = 2, SIZE(helium_env)
         helium_env(k)%helium%num_accepted(:, :) = helium_env(1)%helium%num_accepted(:, :)
      END DO

      !force_avrg:
      IF (helium_env(1)%helium%solute_present) THEN
         IF (helium_env(1)%helium%get_helium_forces == helium_forces_average) THEN
            DO k = 2, SIZE(helium_env)
               helium_env(1)%helium%force_avrg(:, :) = helium_env(1)%helium%force_avrg(:, :)+ &
                                                       helium_env(k)%helium%force_avrg(:, :)
            END DO
            CALL mp_sum(helium_env(1)%helium%force_avrg(:, :), helium_env(1)%comm)
            helium_env(1)%helium%force_avrg(:, :) = helium_env(1)%helium%force_avrg(:, :)*inv_num_env
            DO k = 2, SIZE(helium_env)
               helium_env(k)%helium%force_avrg(:, :) = helium_env(1)%helium%force_avrg(:, :)
            END DO
         ELSE IF (helium_env(1)%helium%get_helium_forces == helium_forces_last) THEN
            IF (logger%para_env%ionode) THEN
               sel_mp_source = INT(next_random_number(helium_env(1)%helium%rng_stream_uniform) &
                                   *REAL(helium_env(1)%helium%num_env, dp))
            END IF
            CALL mp_bcast(sel_mp_source, logger%para_env%source, helium_env(1)%comm)

            offset = 0
            DO i = 1, logger%para_env%mepos
               offset = offset+helium_env(1)%env_all(i)
            END DO
            ALLOCATE (work_2d(SIZE(helium_env(1)%helium%force_avrg, 1), &
                              SIZE(helium_env(1)%helium%force_avrg, 2)))
            work_2d(:, :) = 0.0_dp
            IF (sel_mp_source .GE. offset .AND. sel_mp_source .LT. offset+SIZE(helium_env)) THEN
               work_2d(:, :) = helium_env(sel_mp_source-offset+1)%helium%force_avrg(:, :)
            END IF
            CALL mp_sum(work_2d(:, :), helium_env(1)%comm)
            DO k = 1, SIZE(helium_env)
               helium_env(k)%helium%force_avrg(:, :) = work_2d(:, :)
            END DO
            DEALLOCATE (work_2d)
         END IF
      END IF

      RETURN
   END SUBROUTINE helium_sample

! ***************************************************************************
!> \brief  Perform MC step for helium
!> \param helium_env ...
!> \param pint_env ...
!> \date   2009-06-12
!> \par    History
!>         2016-07-14 Modified to work with independent helium_env [cschran]
!> \author Lukasz Walewski
! **************************************************************************************************
   SUBROUTINE helium_step(helium_env, pint_env)

      TYPE(helium_solvent_p_type), DIMENSION(:), POINTER :: helium_env
      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_step', routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: msgstr, stmp, time_unit
      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: time_start, time_stop, time_used
      REAL(KIND=dp), DIMENSION(:), POINTER               :: DATA

      CALL timeset(routineN, handle)
      time_start = m_walltime()

      IF (ASSOCIATED(helium_env)) THEN
         ! perform the actual phase space sampling
         CALL helium_sample(helium_env, pint_env)

         ! print properties
         CALL helium_print_energy(helium_env)
         CALL helium_print_plength(helium_env)
         CALL helium_print_accepts(helium_env)
         CALL helium_print_force(helium_env)
         CALL helium_print_perm(helium_env)
         CALL helium_print_coordinates(helium_env)
         CALL helium_print_action(pint_env, helium_env)

         IF (helium_env(1)%helium%rdf_present) CALL helium_print_rdf(helium_env)
         IF (helium_env(1)%helium%rho_present) CALL helium_print_rho(helium_env)

         NULLIFY (DATA)
         ALLOCATE (DATA(3*SIZE(helium_env)))

         WRITE (stmp, *) helium_env(1)%helium%apref
         DATA(:) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DATA((k-1)*3+1:k*3) = helium_env(k)%helium%proarea%inst(:)
         END DO
         CALL helium_print_vector(helium_env, &
                                  "MOTION%PINT%HELIUM%PRINT%PROJECTED_AREA", &
                                  DATA, &
                                  angstrom*angstrom, &
                                  (/"A_x", "A_y", "A_z"/), &
                                  "prefactor = "//TRIM(ADJUSTL(stmp))//" [Angstrom^-2]", &
                                  "helium-pa")

         DATA(:) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DATA((k-1)*3+1:k*3) = helium_env(k)%helium%prarea2%ravr(:)
         END DO
         CALL helium_print_vector(helium_env, &
                                  "MOTION%PINT%HELIUM%PRINT%PROJECTED_AREA_2_AVG", &
                                  DATA, &
                                  angstrom*angstrom*angstrom*angstrom, &
                                  (/"<A_x^2>", "<A_y^2>", "<A_z^2>"/), &
                                  "prefactor = "//TRIM(ADJUSTL(stmp))//" [Angstrom^-2]", &
                                  "helium-pa-avg", &
                                  "REWIND", &
                                  .TRUE.)

         DATA(:) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DATA((k-1)*3+1:k*3) = helium_env(k)%helium%mominer%inst(:)
         END DO
         CALL helium_print_vector(helium_env, &
                                  "MOTION%PINT%HELIUM%PRINT%MOMENT_OF_INERTIA", &
                                  DATA, &
                                  angstrom*angstrom, &
                                  (/"I_x/m", "I_y/m", "I_z/m"/), &
                                  "prefactor = "//TRIM(ADJUSTL(stmp))//" [Angstrom^-2]", &
                                  "helium-mi")

         DATA(:) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DATA((k-1)*3+1:k*3) = helium_env(k)%helium%mominer%ravr
         END DO
         CALL helium_print_vector(helium_env, &
                                  "MOTION%PINT%HELIUM%PRINT%MOMENT_OF_INERTIA_AVG", &
                                  DATA, &
                                  angstrom*angstrom, &
                                  (/"<I_x/m>", "<I_y/m>", "<I_z/m>"/), &
                                  "prefactor = "//TRIM(ADJUSTL(stmp))//" [Angstrom^-2]", &
                                  "helium-mi-avg", &
                                  "REWIND", &
                                  .TRUE.)

         DATA(:) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DATA((k-1)*3+1:k*3) = helium_env(k)%helium%wnumber%inst
         END DO
         WRITE (stmp, *) helium_env(1)%helium%wpref
         CALL helium_print_vector(helium_env, &
                                  "MOTION%PINT%HELIUM%PRINT%WINDING_NUMBER", &
                                  DATA, &
                                  angstrom, &
                                  (/"W_x", "W_y", "W_z"/), &
                                  "prefactor = "//TRIM(ADJUSTL(stmp))//" [Angstrom^-2]", &
                                  "helium-wn")

         DATA(:) = 0.0_dp
         DO k = 1, SIZE(helium_env)
            DATA((k-1)*3+1:k*3) = helium_env(k)%helium%wnmber2%ravr
         END DO
         CALL helium_print_vector(helium_env, &
                                  "MOTION%PINT%HELIUM%PRINT%WINDING_NUMBER_2_AVG", &
                                  DATA, &
                                  angstrom*angstrom, &
                                  (/"<W_x^2>", "<W_y^2>", "<W_z^2>"/), &
                                  "prefactor = "//TRIM(ADJUSTL(stmp))//" [Angstrom^-2]", &
                                  "helium-wn-avg", &
                                  "REWIND", &
                                  .TRUE.)
         DEALLOCATE (DATA)

         time_stop = m_walltime()
         time_used = time_stop-time_start
         time_unit = "sec"
         IF (time_used .GE. 60.0_dp) THEN
            time_used = time_used/60.0_dp
            time_unit = "min"
            IF (time_used .GE. 60.0_dp) THEN
               time_used = time_used/60.0_dp
               time_unit = "hours"
            END IF
         END IF
         msgstr = "MC step"
         stmp = ""
         WRITE (stmp, *) helium_env(1)%helium%current_step
         msgstr = TRIM(ADJUSTL(msgstr))//" "//TRIM(ADJUSTL(stmp))//" of"
         stmp = ""
         WRITE (stmp, *) helium_env(1)%helium%last_step
         msgstr = TRIM(ADJUSTL(msgstr))//" "//TRIM(ADJUSTL(stmp))//" in"
         stmp = ""
         WRITE (stmp, '(F20.1)') time_used
         msgstr = TRIM(ADJUSTL(msgstr))//" "//TRIM(ADJUSTL(stmp))
         msgstr = TRIM(ADJUSTL(msgstr))//" "//TRIM(ADJUSTL(time_unit))//"."
         CALL helium_write_line(TRIM(msgstr))

         ! print out the total energy - for regtest evaluation
         stmp = ""
         WRITE (stmp, *) helium_env(1)%helium%energy_avrg(e_id_total)
         msgstr = "Total energy = "//TRIM(ADJUSTL(stmp))
         CALL helium_write_line(TRIM(msgstr))
      END IF

      CALL timestop(handle)

      RETURN
   END SUBROUTINE helium_step

! ***************************************************************************
!> \brief  ...
!> \param helium ...
!> \param  pint_env    path integral environment
!> \par    History
!>        2010-06-17 added different distributions for m-sampling [lwalewski]
!>        2010-06-17 ratio for m-value added (m-sampling related) [lwalewski]
!> \author hforbert
! **************************************************************************************************
   SUBROUTINE helium_try_permutations(helium, pint_env)
      TYPE(helium_solvent_type), POINTER                 :: helium
      TYPE(pint_env_type), POINTER                       :: pint_env

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_try_permutations', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: err_str, stmp
      INTEGER                                            :: cyclen, ni
      LOGICAL                                            :: accepted, selected
      REAL(KIND=dp)                                      :: r, rnd, x, y, z

      IF (helium%maxcycle > 1) CALL helium_update_transition_matrix(helium)

      ! consecutive calls to helium_slice_metro_cyclic require that
      ! ALL(helium%pos.EQ.helium%work) - see a comment there,
      ! besides there is a magic regarding helium%permutation
      ! you have been warned
      helium%work(:, :, :) = helium%pos(:, :, :)

      ! the inner MC loop (without rotation in imaginary time)
      DO ni = 1, helium%iter_norot

         ! set the probability threshold for m_value: 1/(1+(m-1)/helium%m_ratio)
         r = 1.0_dp/(1.0_dp+(helium%maxcycle-1)/helium%m_ratio)

         ! draw permutation length for this trial from the distribution of choice
         !
         SELECT CASE (helium%m_dist_type)

         CASE (helium_mdist_singlev)
            x = next_random_number(helium%rng_stream_uniform)
            IF (x .LT. r) THEN
               cyclen = 1
            ELSE
               cyclen = helium%m_value
            END IF

         CASE (helium_mdist_uniform)
            x = next_random_number(helium%rng_stream_uniform)
            IF (x .LT. r) THEN
               cyclen = helium%m_value
            ELSE
               DO
                  x = next_random_number(helium%rng_stream_uniform)
                  cyclen = INT(helium%maxcycle*x)+1
                  IF (cyclen .NE. helium%m_value) EXIT
               END DO
            END IF

         CASE (helium_mdist_linear)
            x = next_random_number(helium%rng_stream_uniform)
            IF (x .LT. r) THEN
               cyclen = helium%m_value
            ELSE
               DO
                  x = next_random_number(helium%rng_stream_uniform)
                  y = SQRT(2.0_dp*x)
                  cyclen = INT(helium%maxcycle*y/SQRT(2.0_dp))+1
                  IF (cyclen .NE. helium%m_value) EXIT
               END DO
            END IF

         CASE (helium_mdist_quadratic)
            x = next_random_number(helium%rng_stream_uniform)
            IF (x .LT. r) THEN
               cyclen = helium%m_value
            ELSE
               DO
                  x = next_random_number(helium%rng_stream_uniform)
                  y = (3.0_dp*x)**(1.0_dp/3.0_dp)
                  z = 3.0_dp**(1.0_dp/3.0_dp)
                  cyclen = INT(helium%maxcycle*y/z)+1
                  IF (cyclen .NE. helium%m_value) EXIT
               END DO
            END IF

         CASE (helium_mdist_exponential)
            x = next_random_number(helium%rng_stream_uniform)
            IF (x .LT. r) THEN
               cyclen = helium%m_value
            ELSE
               DO
                  DO
                     x = next_random_number(helium%rng_stream_uniform)
                     IF (x .GE. 0.01_dp) EXIT
                  END DO
                  z = -LOG(0.01_dp)
                  y = LOG(x)/z+1.0_dp; 
                  cyclen = INT(helium%maxcycle*y)+1
                  IF (cyclen .NE. helium%m_value) EXIT
               END DO
            END IF

         CASE (helium_mdist_gaussian)
            x = next_random_number(helium%rng_stream_uniform)
            IF (x .LT. r) THEN
               cyclen = 1
            ELSE
               DO
                  x = next_random_number(helium%rng_stream_gaussian)
                  cyclen = INT(x*0.75_dp+helium%m_value-0.5_dp)+1
                  IF (cyclen .NE. 1) EXIT
               END DO
            END IF

         CASE DEFAULT
            WRITE (stmp, *) helium%m_dist_type
            err_str = "M distribution type unknown ("//TRIM(ADJUSTL(stmp))//")"
            CPABORT(err_str)
            cyclen = 1

         END SELECT

         IF (cyclen < 1) cyclen = 1
         IF (cyclen > helium%maxcycle) cyclen = helium%maxcycle
         helium%num_accepted(1, cyclen) = helium%num_accepted(1, cyclen)+1

         ! check, if permutation of this length can be constructed
         IF (cyclen == 1) THEN
            rnd = next_random_number(helium%rng_stream_uniform)
            helium%ptable(1) = 1+INT(rnd*helium%atoms)
            helium%ptable(2) = -1
            helium%pweight = 0.0_dp
            selected = .TRUE.
         ELSE
            selected = helium_select_permutation(helium, cyclen)
         END IF

         IF (selected) THEN
            ! permutation was successfully selected - actually sample this permutation
            accepted = helium_slice_metro_cyclic(helium, pint_env, cyclen)
         ELSE
            accepted = .FALSE.
         END IF

!if (any(helium%pos .ne. helium%work)) then
!  print *, ""
!  print *, "pos and work are different!"
!  print *, ""
!end if

         IF (accepted) THEN
            ! positions changed
            IF (helium%solute_present) THEN
               ! use COM of the solute, but it did not move here so do nothing to save cpu
            ELSE
               IF (helium%periodic) THEN
                  ! use center of the cell, but it did not move here so do nothing to save cpu
               ELSE
                  ! update center of mass
                  helium%center(:) = helium_com(helium)
               END IF
            END IF
         END IF

      END DO

      RETURN
   END SUBROUTINE helium_try_permutations

! ***************************************************************************
!> \brief  Sample path variables for permutation of length <cyclen>
!> \param helium ...
!> \param pint_env ...
!> \param cyclen ...
!> \retval res ...
!> \author hforbert
! **************************************************************************************************
   FUNCTION helium_slice_metro_cyclic(helium, pint_env, cyclen) RESULT(res)
      TYPE(helium_solvent_type), POINTER                 :: helium
      TYPE(pint_env_type), POINTER                       :: pint_env
      INTEGER, INTENT(IN)                                :: cyclen
      LOGICAL                                            :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'helium_slice_metro_cyclic', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: err_str, stmp
      INTEGER                                            :: c, ifix, j, k, l, level, pk1, pk2, stride
      INTEGER, DIMENSION(:), POINTER                     :: p, perm
      LOGICAL                                            :: nperiodic, should_reject
      REAL(KIND=dp)                                      :: cell_size, ds, dtk, e1, e2, pds, &
                                                            prev_ds, r, rtmp, sigma, x
      REAL(KIND=dp), DIMENSION(3)                        :: bis, rm1, rm2, tmp1, tmp2
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: pos, work
      TYPE(spline_data_p_type), DIMENSION(:, :), POINTER :: uij

! trial permutation implicit in p
! since we (momentarily) only use cyclic permutations:
! cyclen = 1 : no permutation, sample p[0] anew
! cyclen = 2 : p[0] -> p[1] -> p[0]
! cyclen = 3 : p[0] -> p[1] -> p[2] -> p[0]
! cyclen = 4 : p[0] -> p[1] -> p[2] -> p[3] -> p[0]

      p => helium%ptable
      prev_ds = helium%pweight

      helium%num_accepted(2, cyclen) = helium%num_accepted(2, cyclen)+1
      level = 1
      res = .FALSE.

      ! nothing to be done
      IF (helium%bisection == 0) RETURN

      ! On entry helium_slice_metro_cyclic ASSUMES that work is a copy of pos
      ! and constructs the trial move there. If the move is accepted, the
      ! changed parts are copied to pos, but if it fails, only the CHANGED parts
      ! of work are overwritten by the corresponding parts of pos. So on exit work
      ! will AGAIN be a copy of pos (either way accept/reject). This is done here
      ! so we do not have to copy around the whole pos array in the caller, and
      ! the caller does not need to know which parts of work might have changed.
      pos => helium%pos
      work => helium%work
      perm => helium%permutation
      uij => helium%uij
      cell_size = (0.5_dp*helium%cell_size)**2
      nperiodic = .NOT. helium%periodic

      pds = prev_ds
      ifix = helium%beads-helium%bisection+1

      ! sanity checks
      !
      IF (ifix < 1) THEN
         WRITE (stmp, *) ifix
         err_str = "ifix<1 test failed (ifix="//TRIM(ADJUSTL(stmp))//")"
         CPABORT(err_str)
      END IF
      !
      j = 1
      k = helium%bisection
      DO
         IF (k < 2) EXIT
         j = j*2
         k = k/2
      END DO
      !
      IF (j /= helium%bisection) THEN
         WRITE (stmp, *) helium%bisection
         err_str = "helium%bisection not a power of 2 (helium%bisection="//TRIM(ADJUSTL(stmp))//")"
         CPABORT(err_str)
      END IF
      !
      IF (helium%bisection < 2) THEN
         WRITE (stmp, *) helium%bisection
         err_str = "helium%bisection less than 2 (helium%bisection="//TRIM(ADJUSTL(stmp))//")"
         CPABORT(err_str)
      END IF

      stride = helium%bisection
      DO
         IF (stride <= 2) EXIT
         ! calc new trial positions
         ! trial action: 0.5*stride*endpointapprox
         sigma = SQRT(0.5_dp*helium%hb2m*(stride/2)*helium%tau)
         dtk = 0.0_dp
         ds = 0.0_dp

         j = ifix+stride/2
         DO
            IF (j > helium%beads-stride/2) EXIT
            pk1 = j-stride/2
            pk2 = j+stride/2
            ! calculate log(T(s)):
            DO k = 1, cyclen
               CALL helium_boxmean_3d(helium, pos(:, p(k), pk1), pos(:, p(k), pk2), bis)
               tmp1(:) = bis(:)-pos(:, p(k), j)
               CALL helium_pbc(helium, tmp1)
               tmp1(:) = tmp1(:)/sigma
               dtk = dtk-0.5_dp*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
            END DO
            ! calculate log(T(sprime)) and sprime itself
            DO k = 1, cyclen
               CALL helium_boxmean_3d(helium, work(:, p(k), pk1), work(:, p(k), pk2), tmp1)
               DO c = 1, 3
                  x = next_random_number(rng_stream=helium%rng_stream_gaussian, variance=1.0_dp)
                  x = sigma*x
                  tmp1(c) = tmp1(c)+x
                  tmp2(c) = x
               END DO
               CALL helium_pbc(helium, tmp1)
               CALL helium_pbc(helium, tmp2)
               work(:, p(k), j) = tmp1(:)
               tmp2(:) = tmp2(:)/sigma
               dtk = dtk+0.5_dp*(tmp2(1)*tmp2(1)+tmp2(2)*tmp2(2)+tmp2(3)*tmp2(3))
            END DO
            j = j+stride
         END DO

         j = helium%beads-stride/2+1
         pk1 = j-stride/2
         DO k = 1, cyclen
            CALL helium_boxmean_3d(helium, pos(:, p(k), pk1), pos(:, perm(p(k)), 1), bis)
            tmp1(:) = bis(:)-pos(:, p(k), j)
            CALL helium_pbc(helium, tmp1)
            tmp1(:) = tmp1(:)/sigma
            dtk = dtk-0.5_dp*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
         END DO
         DO k = 1, cyclen
            CALL helium_boxmean_3d(helium, work(:, p(k), pk1), work(:, perm(p(1+MOD(k, cyclen))), 1), tmp1)
            DO c = 1, 3
               x = next_random_number(rng_stream=helium%rng_stream_gaussian, variance=1.0_dp)
               x = sigma*x
               tmp1(c) = tmp1(c)+x
               tmp2(c) = x
            END DO
            CALL helium_pbc(helium, tmp1)
            CALL helium_pbc(helium, tmp2)
            work(:, p(k), j) = tmp1(:)
            tmp2(:) = tmp2(:)/sigma
            dtk = dtk+0.5_dp*(tmp2(1)*tmp2(1)+tmp2(2)*tmp2(2)+tmp2(3)*tmp2(3))
         END DO
         ! ok we got the new positions
         ! calculate action_k(s)-action_k(sprime)
         x = 1.0_dp/(helium%tau*helium%hb2m*stride)
         j = ifix
         DO
            IF (j > helium%beads-stride/2) EXIT
            pk1 = j+stride/2
            DO k = 1, cyclen
               tmp1(:) = pos(:, p(k), j)-pos(:, p(k), pk1)
               CALL helium_pbc(helium, tmp1)
               ds = ds+x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
               tmp1(:) = work(:, p(k), j)-work(:, p(k), pk1)
               CALL helium_pbc(helium, tmp1)
               ds = ds-x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
               ! interaction change
               IF (helium%solute_present) THEN
                  CALL helium_bead_solute_e_f(pint_env, helium, p(k), pk1, energy=e1)
                  CALL helium_bead_solute_e_f(pint_env, helium, p(k), pk1, work(:, p(k), pk1), e2)
                  ds = ds+(stride/2)*(e1-e2)*helium%tau
               END IF
               DO l = 1, helium%atoms
                  IF (l /= p(k)) THEN
                     tmp1(:) = pos(:, p(k), pk1)-pos(:, l, pk1)
                     CALL helium_pbc(helium, tmp1)
                     r = tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3)
                     IF ((r < cell_size) .OR. nperiodic) THEN
                        r = SQRT(r)
                        ds = ds+REAL(stride/2, dp)*helium_spline(uij(1, 1)%spline_data, r)
                     END IF
                     tmp1(:) = work(:, p(k), pk1)-work(:, l, pk1)
                     CALL helium_pbc(helium, tmp1)
                     r = tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3)
                     IF ((r < cell_size) .OR. nperiodic) THEN
                        r = SQRT(r)
                        ds = ds-REAL(stride/2, dp)*helium_spline(uij(1, 1)%spline_data, r)
                     END IF
                  END IF
               END DO
               ! counted p[k], p[m] twice. subtract those again
               IF (k < cyclen) THEN
                  DO l = k+1, cyclen
                     tmp1(:) = pos(:, p(k), pk1)-pos(:, p(l), pk1)
                     CALL helium_pbc(helium, tmp1)
                     r = tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3)
                     IF ((r < cell_size) .OR. nperiodic) THEN
                        r = SQRT(r)
                        ds = ds-REAL(stride/2, dp)*helium_spline(uij(1, 1)%spline_data, r)
                     END IF
                     tmp1(:) = work(:, p(k), pk1)-work(:, p(l), pk1)
                     CALL helium_pbc(helium, tmp1)
                     r = tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3)
                     IF ((r < cell_size) .OR. nperiodic) THEN
                        r = SQRT(r)
                        ds = ds+REAL(stride/2, dp)*helium_spline(uij(1, 1)%spline_data, r)
                     END IF
                  END DO
               END IF
            END DO
            j = j+stride/2
         END DO
         ! last link
         pk1 = helium%beads-stride/2+1
         DO k = 1, cyclen
            tmp1(:) = pos(:, p(k), pk1)-pos(:, perm(p(k)), 1)
            CALL helium_pbc(helium, tmp1)
            ds = ds+x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
            tmp1(:) = work(:, p(k), pk1)-work(:, perm(p(1+MOD(k, cyclen))), 1)
            CALL helium_pbc(helium, tmp1)
            ds = ds-x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
         END DO
         ! ok now accept or reject:
         rtmp = next_random_number(helium%rng_stream_uniform)
!      IF ((dtk+ds-pds < 0.0_dp).AND.(EXP(dtk+ds-pds)<rtmp)) THEN
         IF (dtk+ds-pds < 0.0_dp) THEN
            IF (EXP(dtk+ds-pds) < rtmp) THEN
               DO c = ifix, helium%beads
                  DO k = 1, cyclen
                     work(:, p(k), c) = pos(:, p(k), c)
                  END DO
               END DO
               RETURN
            END IF
         END IF
         ! accepted. go on to the next level
         helium%num_accepted(level+2, cyclen) = helium%num_accepted(level+2, cyclen)+1
         level = level+1
         pds = ds
         stride = stride/2
      END DO
      ! we are on the lowest level now
      ! calc new trial positions
      ! trial action: endpointapprox for T, full action otherwise
      sigma = SQRT(0.5_dp*helium%hb2m*helium%tau)
      dtk = 0.0_dp
      ds = 0.0_dp
      DO j = ifix+1, helium%beads-1, 2
         pk1 = j-1
         pk2 = j+1
         ! calculate log(T(s)):
         DO k = 1, cyclen
            CALL helium_boxmean_3d(helium, pos(:, p(k), pk1), pos(:, p(k), pk2), bis)
            tmp1(:) = bis(:)-pos(:, p(k), j)
            CALL helium_pbc(helium, tmp1)
            tmp1(:) = tmp1(:)/sigma
            dtk = dtk-0.5_dp*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
         END DO
         ! calculate log(T(sprime)) and sprime itself
         DO k = 1, cyclen
            CALL helium_boxmean_3d(helium, work(:, p(k), pk1), work(:, p(k), pk2), tmp1)
            DO c = 1, 3
               x = next_random_number(rng_stream=helium%rng_stream_gaussian, variance=1.0_dp)
               x = sigma*x
               tmp1(c) = tmp1(c)+x
               tmp2(c) = x
            END DO
            CALL helium_pbc(helium, tmp1)
            CALL helium_pbc(helium, tmp2)
            work(:, p(k), j) = tmp1(:)
            tmp2(:) = tmp2(:)/sigma
            dtk = dtk+0.5_dp*(tmp2(1)*tmp2(1)+tmp2(2)*tmp2(2)+tmp2(3)*tmp2(3))
         END DO
      END DO
      j = helium%beads
      pk1 = j-1
      DO k = 1, cyclen
         CALL helium_boxmean_3d(helium, pos(:, p(k), pk1), pos(:, perm(p(k)), 1), bis)
         tmp1(:) = bis(:)-pos(:, p(k), j)
         CALL helium_pbc(helium, tmp1)
         tmp1(:) = tmp1(:)/sigma
         dtk = dtk-0.5_dp*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
      END DO
      DO k = 1, cyclen
         CALL helium_boxmean_3d(helium, work(:, p(k), pk1), work(:, perm(p(1+MOD(k, cyclen))), 1), tmp1)
         DO c = 1, 3
            x = next_random_number(rng_stream=helium%rng_stream_gaussian, variance=1.0_dp)
            x = sigma*x
            tmp1(c) = tmp1(c)+x
            tmp2(c) = x
         END DO
         CALL helium_pbc(helium, tmp1)
         CALL helium_pbc(helium, tmp2)
         work(:, p(k), j) = tmp1(:)
         tmp2 = tmp2/sigma
         dtk = dtk+0.5_dp*(tmp2(1)*tmp2(1)+tmp2(2)*tmp2(2)+tmp2(3)*tmp2(3))
      END DO
      ! ok we got the new positions.
      ! calculate action_k(s)-action_k(sprime)
      ! interaction change
!TODO interaction ONLY here? or some simplified 12-6 in the upper part?
      IF (helium%solute_present) THEN
         DO j = ifix, helium%beads
            DO k = 1, cyclen
               CALL helium_bead_solute_e_f(pint_env, helium, p(k), j, energy=e1)
               CALL helium_bead_solute_e_f(pint_env, helium, p(k), j, work(:, p(k), j), e2)
               ds = ds+(e1-e2)*helium%tau
            END DO
         END DO
      END IF
      x = 1.0_dp/(helium%tau*helium%hb2m*stride)
      DO j = ifix, helium%beads-1
         pk1 = j+1
         DO k = 1, cyclen
            tmp1(:) = pos(:, p(k), j)-pos(:, p(k), pk1)
            CALL helium_pbc(helium, tmp1)
            ds = ds+x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
            tmp1(:) = work(:, p(k), j)-work(:, p(k), pk1)
            CALL helium_pbc(helium, tmp1)
            ds = ds-x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
            DO l = 1, helium%atoms
               IF (l /= p(k)) THEN
                  rm1(:) = pos(:, p(k), j)-pos(:, l, j)
                  rm2(:) = pos(:, p(k), pk1)-pos(:, l, pk1)
                  ds = ds+helium_eval_expansion(helium, rm1, rm2, uij, 1)
                  rm1(:) = work(:, p(k), j)-work(:, l, j)
                  rm2(:) = work(:, p(k), pk1)-work(:, l, pk1)
                  ds = ds-helium_eval_expansion(helium, rm1, rm2, uij, 1)
               END IF
            END DO
            ! counted p[k], p[m] twice. subtract those again
            IF (k < cyclen) THEN
               DO l = k+1, cyclen
                  rm1(:) = pos(:, p(k), j)-pos(:, p(l), j)
                  rm2(:) = pos(:, p(k), pk1)-pos(:, p(l), pk1)
                  ds = ds-helium_eval_expansion(helium, rm1, rm2, uij, 1)
                  rm1(:) = work(:, p(k), j)-work(:, p(l), j)
                  rm2(:) = work(:, p(k), pk1)-work(:, p(l), pk1)
                  ds = ds+helium_eval_expansion(helium, rm1, rm2, uij, 1)
               END DO
            END IF
         END DO
      END DO
      j = helium%beads
      pk1 = 1
      DO k = 1, cyclen
         tmp1(:) = pos(:, p(k), j)-pos(:, perm(p(k)), 1)
         CALL helium_pbc(helium, tmp1)
         ds = ds+x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
         tmp1(:) = work(:, p(k), j)-work(:, perm(p(1+MOD(k, cyclen))), 1)
         CALL helium_pbc(helium, tmp1)
         ds = ds-x*(tmp1(1)*tmp1(1)+tmp1(2)*tmp1(2)+tmp1(3)*tmp1(3))
         DO l = 1, helium%atoms
            IF (l /= p(k)) THEN
               rm1(:) = pos(:, p(k), j)-pos(:, l, j)
               rm2(:) = pos(:, perm(p(k)), 1)-pos(:, perm(l), 1)
               ds = ds+helium_eval_expansion(helium, rm1, rm2, uij, 1)
            END IF
         END DO
         ! counted p[k], p[m] twice. subtract those again
         IF (k < cyclen) THEN
            DO l = k+1, cyclen
               rm1(:) = pos(:, p(k), j)-pos(:, p(l), j)
               rm2(:) = pos(:, perm(p(k)), pk1)-pos(:, perm(p(l)), pk1)
               ds = ds-helium_eval_expansion(helium, rm1, rm2, uij, 1)
            END DO
         END IF
      END DO
      IF (cyclen > 1) THEN
         !k,c,l
         c = perm(p(1))
         DO k = 1, cyclen-1
            perm(p(k)) = perm(p(k+1))
         END DO
         perm(p(cyclen)) = c
      END IF
      DO k = 1, cyclen
         DO l = 1, helium%atoms
            IF (l /= p(k)) THEN
               rm1(:) = work(:, p(k), j)-work(:, l, j)
               rm2(:) = work(:, perm(p(k)), 1)-work(:, perm(l), 1)
               ds = ds-helium_eval_expansion(helium, rm1, rm2, uij, 1)
            END IF
         END DO
         ! counted p[k], p[m] twice. subtract those again
         IF (k < cyclen) THEN
            DO l = k+1, cyclen
               rm1(:) = work(:, p(k), j)-work(:, p(l), j)
               rm2(:) = work(:, perm(p(k)), 1)-work(:, perm(p(l)), 1)
               ds = ds+helium_eval_expansion(helium, rm1, rm2, uij, 1)
            END DO
         END IF
      END DO
      ! ok now accept or reject:
      rtmp = next_random_number(helium%rng_stream_uniform)
!   IF ((dtk+ds-pds<0.0_dp).AND.(EXP(dtk+ds-pds)<rtmp)) THEN
      IF (dtk+ds-pds < 0.0_dp) THEN
         IF (EXP(dtk+ds-pds) < rtmp) THEN
            DO c = ifix, helium%beads
               DO k = 1, cyclen
                  work(:, p(k), c) = pos(:, p(k), c)
               END DO
            END DO
            IF (cyclen > 1) THEN
               c = perm(p(cyclen))
               DO k = cyclen-1, 1, -1
                  perm(p(k+1)) = perm(p(k))
               END DO
               perm(p(1)) = c
            END IF
            RETURN
         END IF
      END IF
      ! accepted.

      ! rejection of the whole move if any of the new positions is farther away
      ! from the current center of gravity than HELIUM%DROPLET_RADIUS [lwalewski]
      IF (helium%droplet_radius .LT. HUGE(1.0_dp)) THEN
         ! cycle through new positions
         should_reject = .FALSE.
         outer: DO c = ifix, helium%beads
            DO k = 1, cyclen
               bis(:) = work(:, p(k), c)-helium%center(:)
               rtmp = SQRT(bis(1)*bis(1)+bis(2)*bis(2)+bis(3)*bis(3))
               IF (rtmp .GE. helium%droplet_radius) THEN
                  ! found - this one does not fit within R from the center
                  should_reject = .TRUE.
                  EXIT outer
               END IF
            END DO
         END DO outer
         IF (should_reject) THEN
            ! restore work and perm, then return
            DO c = ifix, helium%beads
               DO k = 1, cyclen
                  work(:, p(k), c) = pos(:, p(k), c)
               END DO
            END DO
            IF (cyclen > 1) THEN
               c = perm(p(cyclen))
               DO k = cyclen-1, 1, -1
                  perm(p(k+1)) = perm(p(k))
               END DO
               perm(p(1)) = c
            END IF
            RETURN
         END IF
         ! Note: after CoM position gets updated according to the new positions
         ! some beads may end up outside of the DROPLET_RADIUS range from the new
         ! CoM, they should not be very far off, though.
      END IF
      ! accepted.

      ! copy trial over to the real thing
      DO c = ifix, helium%beads
         DO k = 1, cyclen
            pos(:, p(k), c) = work(:, p(k), c)
         END DO
      END DO
      DO k = 1, cyclen
         helium%iperm(perm(p(k))) = p(k)
      END DO
      helium%num_accepted(level+2, cyclen) = helium%num_accepted(level+2, cyclen)+1
      res = .TRUE.

      RETURN
   END FUNCTION helium_slice_metro_cyclic

! ***************************************************************************
!> \brief  ...
!> \param helium ...
!> \param len ...
!> \retval res ...
!> \author hforbert
! **************************************************************************************************
   FUNCTION helium_select_permutation(helium, len) RESULT(res)
      TYPE(helium_solvent_type), POINTER                 :: helium
      INTEGER, INTENT(IN)                                :: len
      LOGICAL                                            :: res

      INTEGER                                            :: i, j, k, n
      INTEGER, DIMENSION(:), POINTER                     :: iperm, p, perm
      INTEGER, DIMENSION(:, :), POINTER                  :: nmatrix
      REAL(KIND=dp)                                      :: rnd, s1, s2, t
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ipmatrix, pmatrix, tmatrix

      s1 = 0.0_dp
      s2 = 0.0_dp
      res = .FALSE.
      n = helium%atoms
      tmatrix => helium%tmatrix
      pmatrix => helium%pmatrix
      ipmatrix => helium%ipmatrix
      perm => helium%permutation
      iperm => helium%iperm
      p => helium%ptable
      nmatrix => helium%nmatrix

      p(len+1) = -1
      rnd = next_random_number(helium%rng_stream_uniform)
      p(1) = INT(n*rnd)+1
      DO k = 1, len-1
         t = next_random_number(helium%rng_stream_uniform)
         ! find the corresponding path to connect to
         ! using the precalculated optimal decision tree:
         i = n-1
         DO
            IF (tmatrix(p(k), i) > t) THEN
               i = nmatrix(p(k), 2*i-1)
            ELSE
               i = nmatrix(p(k), 2*i)
            END IF
            IF (i < 0) EXIT
         END DO
         i = -i
         ! which particle was it previously connected to?
         p(k+1) = iperm(i)
         ! is it unique? quit if it was already part of the permutation
         DO j = 1, k
            IF (p(j) == p(k+1)) RETURN
         END DO
         ! acummulate the needed values for the final
         ! accept/reject step:
         s1 = s1+ipmatrix(p(k), i)
         s2 = s2+ipmatrix(p(k), perm(p(k)))
      END DO
      ! close the permutation loop:
      s1 = s1+ipmatrix(p(len), perm(p(1)))
      s2 = s2+ipmatrix(p(len), perm(p(len)))
      ! final accept/reject:
      rnd = next_random_number(helium%rng_stream_uniform)
      t = s1*rnd
      IF (t > s2) RETURN
      ! ok, we have accepted the permutation
      ! calculate the action bias for the subsequent resampling
      ! of the paths:
      s1 = pmatrix(p(len), perm(p(1)))-pmatrix(p(len), perm(p(len)))
      DO k = 1, len-1
         s1 = s1+pmatrix(p(k), perm(p(k+1)))-pmatrix(p(k), perm(p(k)))
      END DO
      helium%pweight = s1
      res = .TRUE.
      RETURN
   END FUNCTION helium_select_permutation

END MODULE helium_sampling
