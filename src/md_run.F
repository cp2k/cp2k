!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Perform a molecular dynamics (MD) run using QUICKSTEP
!> \par History
!>   - Added support for Langevin regions (2014/02/05, LT)
!> \author Matthias Krack (07.11.2002)
! *****************************************************************************
MODULE md_run
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE averages_types,                  ONLY: average_quantities_type
  USE barostat_types,                  ONLY: barostat_type,&
                                             create_barostat_type,&
                                             release_barostat_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc_cp2k_plumed_getset_cell
  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE free_energy_methods,             ONLY: free_energy_evaluate
  USE free_energy_types,               ONLY: fe_env_create,&
                                             free_energy_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       ehrenfest, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, reftraj_ensemble
  USE input_cp2k_check,                ONLY: remove_restart_info
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_walltime
  USE mc_types,                        ONLY: mc_ekin_type
  USE md_ener_types,                   ONLY: create_md_ener,&
                                             md_ener_type,&
                                             release_md_ener
  USE md_energies,                     ONLY: initialize_md_ener,&
                                             md_ener_reftraj,&
                                             md_energy,&
                                             md_write_output
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             need_per_atom_wiener_process,&
                                             set_md_env
  USE md_util,                         ONLY: md_output
  USE md_vel_utils,                    ONLY: angvel_control,&
                                             comvel_control,&
                                             setup_velocities,&
                                             temperature_control
  USE mdctrl_methods,                  ONLY: mdctrl_callback
  USE mdctrl_types,                    ONLY: mdctrl_type
  USE metadynamics,                    ONLY: metadyn_forces,&
                                             metadyn_write_colvar
  USE metadynamics_types,              ONLY: set_meta_env
  USE particle_list_types,             ONLY: particle_list_type
  USE reftraj_types,                   ONLY: create_reftraj,&
                                             reftraj_type,&
                                             release_reftraj
  USE reftraj_util,                    ONLY: initialize_reftraj,&
                                             write_output_reftraj
  USE rt_propagation,                  ONLY: rt_prop_setup,&
                                             rt_write_input_restart
  USE rt_propagation_output,           ONLY: rt_prop_output
  USE simpar_methods,                  ONLY: read_md_section
  USE simpar_types,                    ONLY: create_simpar_type,&
                                             release_simpar_type,&
                                             simpar_type
  USE thermal_region_types,            ONLY: release_thermal_regions,&
                                             thermal_regions_type
  USE thermal_region_utils,            ONLY: create_thermal_regions,&
                                             print_thermal_regions_langevin
  USE thermostat_methods,              ONLY: create_thermostats
  USE thermostat_types,                ONLY: release_thermostats,&
                                             thermostats_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE virial_methods,                  ONLY: virial_evaluate
  USE wiener_process,                  ONLY: create_wiener_process,&
                                             create_wiener_process_cv
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_run'

  PUBLIC :: qs_mol_dyn

CONTAINS

! *****************************************************************************
!> \brief Main driver module for Molecular Dynamics
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_mol_dyn ( force_env, globenv, averages, rm_restart_info, hmc_ekin, mdctrl, error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(average_quantities_type), &
      OPTIONAL, POINTER                      :: averages
    LOGICAL, INTENT(IN), OPTIONAL            :: rm_restart_info
    TYPE(mc_ekin_type), OPTIONAL, POINTER    :: hmc_ekin
    TYPE(mdctrl_type), OPTIONAL, POINTER     :: mdctrl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_rm_restart_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: md_section, motion_section

    failure = .FALSE.
    my_rm_restart_info = .TRUE.
    IF (PRESENT(rm_restart_info)) my_rm_restart_info = rm_restart_info
    NULLIFY(md_env, para_env)
    para_env       => force_env%para_env
    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    md_section     => section_vals_get_subs_vals(motion_section,"MD",error=error)

    ! Real call to MD driver - Low Level
    CALL md_env_create(md_env, md_section, para_env, force_env=force_env, error=error)
    CALL set_md_env(md_env, averages=averages, error=error)
    IF(PRESENT(hmc_ekin))THEN
       CALL qs_mol_dyn_low(md_env, md_section, motion_section, force_env, globenv, &
            hmc_ekin=hmc_ekin, error=error)
    ELSE
       CALL qs_mol_dyn_low(md_env, md_section, motion_section, force_env, globenv, mdctrl=mdctrl, error=error)
    ENDIF
    CALL md_env_release(md_env, error=error)

    ! Clean restartable sections..
    IF (my_rm_restart_info) CALL remove_restart_info(force_env%root_section,error=error)
  END SUBROUTINE qs_mol_dyn

! *****************************************************************************
!> \brief Purpose: Driver routine for MD run using QUICKSTEP.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>   - Cleaning (09.2007) Teodoro Laino [tlaino] - University of Zurich
!>   - Added lines to print out langevin regions (2014/02/04, LT)
!> \author Creation (07.11.2002,MK)
! *****************************************************************************
  SUBROUTINE qs_mol_dyn_low(md_env, md_section, motion_section, force_env, globenv, hmc_ekin, mdctrl, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: md_section, motion_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mc_ekin_type), OPTIONAL, POINTER    :: hmc_ekin
    TYPE(mdctrl_type), OPTIONAL, POINTER     :: mdctrl
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn_low', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: my_act, my_pos
    INTEGER                                  :: handle, i, i_kinds, istep, &
                                                md_stride, natom_plumed, &
                                                nsp_plumed, run_type_id
    INTEGER, DIMENSION(90)                   :: na_plumed
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: check, ehrenfest_md, failure, &
                                                save_mem, should_stop, &
                                                write_binary_restart_file
    REAL(KIND=dp)                            :: dummy, time_iter_start, &
                                                time_iter_stop, &
                                                timestep_plumed
    REAL(KIND=dp), DIMENSION(90)             :: pma0_plumed
    REAL(KIND=dp), POINTER                   :: constant, time, used_time
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER :: constraint_section, &
      force_env_section, free_energy_section, global_section, &
      reftraj_section, subsys_section, work_section
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(thermostats_type), POINTER          :: thermostats

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    failure=.FALSE.
    NULLIFY (particles, cell, cell_ref, simpar, itimes, used_time, subsys, &
         md_ener, thermostats, barostat, reftraj, force_env_section, &
         reftraj_section, work_section, atomic_kinds, &
         local_particles, time, fe_env, free_energy_section, &
         constraint_section, thermal_regions)
    logger   => cp_error_get_logger(error)
    para_env => force_env%para_env

    global_section => section_vals_get_subs_vals(force_env%root_section,"GLOBAL",error=error)
    free_energy_section =>section_vals_get_subs_vals(motion_section,"FREE_ENERGY",error=error)
    constraint_section =>section_vals_get_subs_vals(motion_section,"CONSTRAINT",error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)

    CALL section_vals_val_get(global_section,"RUN_TYPE", i_val=run_type_id,error=error)
    IF(run_type_id==ehrenfest) CALL set_md_env(md_env, ehrenfest_md=.TRUE., error=error)

    CALL create_simpar_type(simpar, error)
    force_env_section => force_env%force_env_section
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    CALL cp_add_iter_level(logger%iter_info,"MD",error=error)
    CALL cp_iterate(logger%iter_info,iter_nr=0,error=error)
    ! Read MD section
    CALL read_md_section(simpar, motion_section, md_section, error)
    ! Setup print_keys
    simpar%info_constraint = cp_print_key_unit_nr(logger,constraint_section,&
         "CONSTRAINT_INFO",extension=".shakeLog",log_filename=.FALSE.,error=error)
    simpar%lagrange_multipliers = cp_print_key_unit_nr(logger,constraint_section,&
         "LAGRANGE_MULTIPLIERS",extension=".LagrangeMultLog",log_filename=.FALSE.,error=error)
    simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info,constraint_section,&
         "LAGRANGE_MULTIPLIERS",error=error),cp_p_file)

    ! Create the structure for the md energies
    CALL create_md_ener(md_ener, error=error)
    CALL set_md_env(md_env, md_ener=md_ener, error=error)
    CALL release_md_ener(md_ener, error=error)

    ! If requested setup Thermostats
    CALL create_thermostats(thermostats,  md_section, force_env, simpar, para_env,&
         globenv, global_section, error )

    ! If requested setup Barostat
    CALL create_barostat_type(barostat, md_section, force_env, simpar, globenv, error )

    ! If requested setup different thermal regions
    CALL create_thermal_regions(thermal_regions,  md_section, simpar, force_env, error )

    ! If doing langevin_ensemble, then print out langevin_regions information upon request
    IF (simpar%ensemble == langevin_ensemble) THEN
       my_pos = "REWIND"
       my_act = "WRITE"
       CALL print_thermal_regions_langevin(thermal_regions, simpar, &
                                           pos=my_pos, act=my_act, error=error)
    END IF

    CALL set_md_env(md_env, thermostats=thermostats, barostat=barostat, thermal_regions=thermal_regions,error=error)

    CALL get_md_env(md_env, ehrenfest_md=ehrenfest_md, error=error)

    !If requested set up the REFTRAJ run
    IF(simpar%ensemble == reftraj_ensemble .AND. ehrenfest_md)&
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Ehrenfest MD does not support reftraj ensemble "//&
                     CPSourceFileRef,&
                     error,failure)
    IF(simpar%ensemble == reftraj_ensemble) THEN
        reftraj_section => section_vals_get_subs_vals(md_section,"REFTRAJ",error=error)
        CALL create_reftraj(reftraj, reftraj_section, para_env, error=error)
        CALL set_md_env(md_env, reftraj=reftraj, error=error)
        CALL release_reftraj(reftraj,error=error)
    END IF

    CALL force_env_get(force_env, subsys=subsys, cell=cell, cell_ref=cell_ref,&
                       force_env_section=force_env_section, error=error )

    ! Set V0 if needed
    IF (simpar%ensemble == nph_uniaxial_ensemble.OR.simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
       CPPrecondition(ASSOCIATED(cell_ref),cp_failure_level,routineP,error,failure)
       IF ( simpar%v0 == 0._dp ) simpar%v0 = cell_ref%deth
    ENDIF

    ! Initialize velocities possibly applying constraints at the zeroth MD step
    CALL section_vals_val_get(motion_section,"PRINT%RESTART%SPLIT_RESTART_FILE",&
                              l_val=write_binary_restart_file,error=error)
    CALL setup_velocities(force_env,simpar,globenv,md_env,md_section,constraint_section, &
                          write_binary_restart_file,error)

    ! Setup Free Energy Calculation (if required)
    CALL fe_env_create (fe_env, free_energy_section, error)


    CALL set_md_env(md_env=md_env, simpar=simpar, fe_env=fe_env, cell=cell,&
                    force_env=force_env, error=error)

    ! Possibly initialize Wiener processes
    ![NB] Tested again within create_wiener_process.  Why??
    IF (need_per_atom_wiener_process(md_env, error=error)) CALL create_wiener_process(md_env,error)

    time_iter_start=m_walltime()

    CALL get_md_env(md_env, force_env=force_env, itimes=itimes, constant=constant,&
                    md_ener=md_ener, t=time, used_time=used_time, error=error)

    ! Attach the time counter of the meta_env to the one of the MD
    CALL set_meta_env(force_env%meta_env, time=time, error=error)

    ! Initialize the md_ener structure
    CALL initialize_md_ener(md_ener, force_env, simpar, error=error)

    ! Check for ensembles requiring the stress tensor - takes into account the possibility for
    ! multiple force_evals
    IF ( (simpar%ensemble==npt_i_ensemble).OR.&
         (simpar%ensemble==npt_f_ensemble).OR.&
         (simpar%ensemble==npe_f_ensemble).OR.&
         (simpar%ensemble==npe_i_ensemble).OR.&
         (simpar%ensemble==nph_uniaxial_ensemble).OR.&
         (simpar%ensemble==nph_uniaxial_damped_ensemble)) THEN
       check = force_env%virial%pv_availability
       CALL cp_assert(check,cp_failure_level,cp_assertion_failed,&
            routineP,"Virial evaluation not requested for this run in the input file! "//&
            "You may consider to switch on the virial evaluation with the keyword: STRESS_TENSOR."//&
            "Be sure the method you are using can compute the virial! "//&
 CPSourceFileRef,&
            error,failure)
       IF (ASSOCIATED(force_env%sub_force_env)) THEN
          DO i = 1, SIZE(force_env%sub_force_env)
             IF (ASSOCIATED(force_env%sub_force_env(i)%force_env)) THEN
                check = check .AND. force_env%sub_force_env(i)%force_env%virial%pv_availability
             END IF
          END DO
       END IF
       CALL cp_assert(check,cp_failure_level,cp_assertion_failed,&
            routineP,"Virial evaluation not requested for all the force_eval sections present in"//&
            " the input file! You have to switch on the virial evaluation with the keyword: STRESS_TENSOR "//&
            " in each force_eval section. Be sure the method you are using can compute the virial!"//&
 CPSourceFileRef,&
            error,failure)
    END IF

    ! Computing Forces at zero MD step
    IF (simpar%ensemble /= reftraj_ensemble) THEN
       CALL section_vals_val_get(md_section,"STEP_START_VAL",i_val=itimes,error=error)
       CALL section_vals_val_get(md_section,"TIME_START_VAL",r_val=time,error=error)
       CALL section_vals_val_get(md_section,"ECONS_START_VAL",r_val=constant,error=error)
       CALL cp_iterate(logger%iter_info,iter_nr=itimes,error=error)
       IF(save_mem) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)
          work_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)
          work_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)
       END IF

       IF(ehrenfest_md)THEN
          CALL rt_prop_setup(force_env,error)
          force_env%qs_env%rtp%dt=simpar%dt
       ELSE
          ![NB] Lets let all methods, even ones without consistent energies, succeed here.  They'll fail in actual integrator if needed
          CALL force_env_calc_energy_force (force_env, calc_force=.TRUE., require_consistent_energy_force=.FALSE., error=error)
       END IF

       IF(ASSOCIATED(force_env%qs_env))THEN
          force_env%qs_env%sim_time=time
          force_env%qs_env%sim_step=itimes
       END IF
       ! Warm-up engines for metadynamics
       IF (ASSOCIATED(force_env%meta_env)) THEN
           ! Setup stuff for plumed if needed
           IF (force_env%meta_env%use_plumed .EQV. .TRUE.) THEN
#ifdef PLUMED_CP2K
               natom_plumed=subsys%particles%n_els
               nsp_plumed=subsys%atomic_kinds%n_els
               timestep_plumed=simpar%dt
               DO i_kinds=1,nsp_plumed
                   na_plumed(i_kinds) = subsys%atomic_kinds%els(i_kinds)%natom
                   pma0_plumed(i_kinds) = subsys%atomic_kinds%els(i_kinds)%mass
               END DO
               CALL pbc_cp2k_plumed_getset_cell (cell, set=.TRUE.) !Store the cell pointer for later use.

               !in a parallel run plumed only needs to be initialized on the
               !master node.
               IF (para_env%ionode .EQV. .TRUE.) THEN
                   CALL INIT_METADYN(natom_plumed, nsp_plumed, na_plumed, &
                       timestep_plumed, itimes, pma0_plumed, force_env%meta_env%plumed_input_file)
               END IF
#else
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                routineP,"Requested to use plumed for metadynamics, but cp2k was"//&
                " not compiled with plumed support."//&
                CPSourceFileRef,&
                error,failure)
#endif
           ELSE
            IF(force_env%meta_env%langevin) THEN
                 CALL create_wiener_process_cv(force_env%meta_env, error=error)
                ENDIF
                IF (force_env%meta_env%well_tempered) THEN
                 force_env%meta_env%wttemperature = simpar%temp_ext
                 IF (force_env%meta_env%wtgamma>EPSILON(1._dp)) THEN
                    dummy=force_env%meta_env%wttemperature*(force_env%meta_env%wtgamma-1._dp)
                    IF (force_env%meta_env%delta_t>EPSILON(1._dp)) THEN
                       check=ABS(force_env%meta_env%delta_t-dummy)<1.E+3_dp*EPSILON(1._dp)
                       CALL cp_assert(check,cp_failure_level,cp_assertion_failed,routineP,&
                          "Inconsistency between DELTA_T and WTGAMMA (both specified):"//&
                          " please, verify that DELTA_T=(WTGAMMA-1)*TEMPERATURE",&
                          error,failure)
                    ELSE
                       force_env%meta_env%delta_t = dummy
                    ENDIF
                 ELSE
                    force_env%meta_env%wtgamma    = 1._dp &
                       + force_env%meta_env%delta_t/force_env%meta_env%wttemperature
                 ENDIF
                 force_env%meta_env%invdt         = 1._dp/force_env%meta_env%delta_t
                ENDIF
                CALL metadyn_forces(force_env,error=error)
                CALL metadyn_write_colvar(force_env,error=error)
            END IF
       ENDIF

       IF (simpar%do_respa)THEN
          CALL force_env_calc_energy_force (force_env%sub_force_env(1)%force_env,&
               calc_force=.TRUE.,error=error)
       END IF

       CALL force_env_get( force_env, subsys=subsys, error=error )

       CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
            particles=particles,error=error)


       CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles,&
            force_env%virial, force_env%para_env%group, error=error)

       CALL md_energy(md_env,md_ener,error)
       CALL md_write_output(md_env, error) !inits the print env at itimes == 0 also writes trajectories
       md_stride = 1
    ELSE
       CALL get_md_env(md_env, reftraj=reftraj, error=error)
       CALL initialize_reftraj(reftraj, reftraj_section, md_env, error=error)
       itimes = reftraj%info%first_snapshot -1
       md_stride = reftraj%info%stride
#ifdef PLUMED_CP2K
       IF (ASSOCIATED(force_env%meta_env)) THEN
           IF (force_env%meta_env%use_plumed .EQV. .TRUE.) THEN
               natom_plumed=subsys%particles%n_els
               nsp_plumed=subsys%atomic_kinds%n_els
               timestep_plumed=simpar%dt
               DO i_kinds=1,nsp_plumed
                   na_plumed(i_kinds) = subsys%atomic_kinds%els(i_kinds)%natom
                   pma0_plumed(i_kinds) = subsys%atomic_kinds%els(i_kinds)%mass
               END DO
               CALL pbc_cp2k_plumed_getset_cell (cell, set=.TRUE.) !Store the cell pointer for later use.

               !in a parallel run plumed only needs to be initialized on the
               !master node.
               IF (para_env%ionode .EQV. .TRUE.) THEN
                   CALL INIT_METADYN(natom_plumed, nsp_plumed, na_plumed, &
                       timestep_plumed, itimes, pma0_plumed, force_env%meta_env%plumed_input_file)
               END IF
            END IF
        END IF
#endif
    END IF

    CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
         constraint_section,"CONSTRAINT_INFO",error=error)
    CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
         constraint_section,"LAGRANGE_MULTIPLIERS",error=error)

! if we need the initial kinetic energy for Hybrid Monte Carlo
    IF(PRESENT(hmc_ekin)) hmc_ekin%initial_ekin=md_ener%ekin

    ! Real MD Loop
    DO istep=1, simpar%nsteps, md_stride
       ! Increase counters
       itimes = itimes + 1
       time   = time + simpar%dt
       !needed when electric field fields are applied
       IF(ASSOCIATED(force_env%qs_env))THEN
          force_env%qs_env%sim_time=time
          force_env%qs_env%sim_step=itimes
       END IF
       IF(ehrenfest_md)force_env%qs_env%rtp%istep=istep

       CALL cp_iterate(logger%iter_info,last=(istep==simpar%nsteps),iter_nr=itimes,error=error)
       ! Open possible Shake output units
       simpar%info_constraint = cp_print_key_unit_nr(logger,constraint_section,"CONSTRAINT_INFO",&
            extension=".shakeLog",log_filename=.FALSE.,error=error)
       simpar%lagrange_multipliers = cp_print_key_unit_nr(logger,constraint_section,&
            "LAGRANGE_MULTIPLIERS",extension=".LagrangeMultLog",log_filename=.FALSE.,error=error)
       simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info,constraint_section,&
            "LAGRANGE_MULTIPLIERS",error=error),cp_p_file)

       ! Velocity Verlet Integrator
       CALL velocity_verlet(md_env,globenv,error)

       ! Close Shake output if requested...
       CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
            constraint_section,"CONSTRAINT_INFO",error=error)
       CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
            constraint_section,"LAGRANGE_MULTIPLIERS",error=error)

       ! Free Energy calculation
       CALL free_energy_evaluate(md_env,should_stop,free_energy_section,error)

       IF (should_stop) EXIT

       ! Test for <PROJECT_NAME>.EXIT_MD or for WALL_TIME to exit
       ! Default:
       ! IF so we don't overwrite the restart or append to the trajectory
       ! because the execution could in principle stop inside the SCF where energy
       ! and forces are not converged.
       ! But:
       ! You can force to print the last step (for example if the method used
       ! to compute energy and forces is not SCF based) activating the print_key
       ! MOTION%MD%PRINT%FORCE_LAST.
       CALL external_control(should_stop,"MD",globenv=globenv,error=error)

       ! call external hook e.g. from global optimization
       IF (PRESENT(mdctrl)) &
          CALL mdctrl_callback(mdctrl, md_env, should_stop, error)

       IF (should_stop) THEN
          CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=itimes,error=error)
          CALL md_output(md_env,md_section,force_env%root_section,should_stop,error=error)
          IF(ehrenfest_md)THEN
             CALL rt_prop_output(force_env%qs_env,ehrenfest,error=error)
             CALL rt_write_input_restart(md_env,force_env,error)
          END IF
          EXIT
       END IF

       IF(simpar%ensemble /= reftraj_ensemble) THEN
          CALL md_energy(md_env, md_ener, error)
          CALL temperature_control(simpar, md_env, md_ener, force_env, logger, error)
          CALL comvel_control(md_ener, force_env, md_section, logger, error)
          CALL angvel_control(md_ener, force_env, md_section, logger, error)
       ELSE
          CALL md_ener_reftraj(md_env, md_ener, error)
       END IF

       time_iter_stop=m_walltime()
       used_time = time_iter_stop - time_iter_start
       time_iter_start=time_iter_stop

       CALL md_output(md_env,md_section,force_env%root_section,should_stop,error=error)
       IF(simpar%ensemble == reftraj_ensemble ) THEN
          CALL write_output_reftraj(md_env,error=error)
       END IF
    END DO

! if we need the final kinetic energy for Hybrid Monte Carlo
    IF(PRESENT(hmc_ekin)) hmc_ekin%final_ekin=md_ener%ekin

    ! Remove the iteration level
    CALL cp_rm_iter_level(logger%iter_info,"MD",error=error)

    ! Deallocate Thermostats and Barostats
    CALL release_thermostats(thermostats, error=error)
    CALL release_barostat_type(barostat, error=error)
    CALL release_simpar_type(simpar, error)
    CALL release_thermal_regions(thermal_regions, error)
    CALL timestop(handle)

  END SUBROUTINE qs_mol_dyn_low

END MODULE md_run
