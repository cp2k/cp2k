!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_run [1.0] *
!!
!!   NAME
!!     md_run
!!
!!   FUNCTION
!!     Perform a molecular dynamics (MD) run using QUICKSTEP
!!
!!   AUTHOR
!!     Matthias Krack (07.11.2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE md_run
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE barostat_types,                  ONLY: barostat_type,&
                                             create_barostat_type,&
                                             release_barostat_type
  USE cell_types,                      ONLY: cell_type
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
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_initialize,&
                                             dyn_coeff_set_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE free_energy_methods,             ONLY: free_energy_evaluate
  USE free_energy_types,               ONLY: fe_env_create,&
                                             free_energy_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: langevin_ensemble,&
                                             reftraj_ensemble
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE md,                              ONLY: read_md_section,&
                                             simulation_parameters_type
  USE md_energies,                     ONLY: initialize_md_ener,&
                                             md_energy,&
                                             md_write_output
  USE md_environment_types,            ONLY: create_md_ener,&
                                             get_md_env,&
                                             md_ener_type,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             release_md_ener,&
                                             set_md_env
  USE md_util,                         ONLY: comvel_control,&
                                             md_output,&
                                             setup_velocities,&
                                             temperature_control
  USE particle_list_types,             ONLY: particle_list_type
  USE reftraj_types,                   ONLY: create_reftraj,&
                                             reftraj_type,&
                                             release_reftraj
  USE reftraj_util,                    ONLY: initialize_reftraj
  USE termination,                     ONLY: external_control
  USE thermostat_methods,              ONLY: create_thermostats
  USE thermostat_types,                ONLY: release_thermostats,&
                                             thermostats_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE virial_methods,                  ONLY: virial_evaluate
  USE wiener_process,                  ONLY: create_wiener_process
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_run'

  PUBLIC :: qs_mol_dyn,&
            qs_mol_dyn_low

CONTAINS

! *****************************************************************************
!!****f* md_run/qs_mol_dyn [1.0] *
!!
!!   NAME
!!     qs_mol_dyn
!!
!!   FUNCTION
!!     Main driver module for Molecular Dynamics
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE qs_mol_dyn ( force_env, globenv, error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn', &
      routineP = moduleN//':'//routineN

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(md_environment_type), POINTER       :: md_env

    NULLIFY(md_env, para_env)
    para_env=> force_env % para_env
    ! Real call to MD driver - Low Level
    CALL md_env_create(md_env, para_env, error=error)
    CALL qs_mol_dyn_low(md_env, force_env, globenv, error)
    CALL md_env_release(md_env, error=error)

  END SUBROUTINE qs_mol_dyn

! *****************************************************************************
!!****f* md_run/qs_mol_dyn_low [1.0] *
!!
!!   NAME
!!     qs_mol_dyn_low
!!
!!   FUNCTION
!!     Purpose: Driver routine for MD run using QUICKSTEP.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Creation (07.11.2002,MK)
!!
!!   MODIFICATION HISTORY
!!     Cleaning (09.2007) Teodoro Laino [tlaino] - University of Zurich
!!
!!*** **********************************************************************
  SUBROUTINE qs_mol_dyn_low ( md_env, force_env, globenv, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, istep, iter_level, &
                                                stat
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, save_mem, should_stop
    REAL(KIND=dp)                            :: time_iter_start, &
                                                time_iter_stop
    REAL(KIND=dp), POINTER                   :: t, used_time
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(section_vals_type), POINTER :: force_env_section, global_section, &
      kg_coeff_section, motion_section, reftraj_section, root_section, &
      subsys_section, work_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(thermostats_type), POINTER          :: thermostats

    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    failure=.FALSE.
    NULLIFY (particles, cell, cell_ref, simpar, itimes, used_time, subsys, &
         md_ener, motion_section, thermostats, barostat, reftraj, force_env_section,&
         kg_coeff_section, reftraj_section, work_section, atomic_kinds, &
         local_particles, t, fe_env, dyn_coeff_set)
    logger => cp_error_get_logger(error)
    root_section => force_env % root_section
    para_env     => force_env % para_env

    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)

    ALLOCATE (simpar,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    force_env_section => force_env%force_env_section
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    CALL cp_add_iter_level(logger%iter_info,"MD",n_rlevel_new=iter_level,&
         error=error)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=0,error=error)
    ! Read MD section
    CALL read_md_section ( simpar, root_section, error )
    ! Setup print_keys
    simpar%info_constraint = cp_print_key_unit_nr(logger,root_section,&
         "MOTION%CONSTRAINT%CONSTRAINT_INFO",&
         extension=".shakeLog",log_filename=.FALSE.,error=error)
    simpar%lagrange_multipliers = cp_print_key_unit_nr(logger,root_section,&
         "MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",&
         extension=".LagrangeMultLog",log_filename=.FALSE.,error=error)
    simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info,root_section,&
         "MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",error=error),cp_p_file)


    ! Create the structure for the md energies
    CALL create_md_ener(md_ener, error=error)
    CALL set_md_env(md_env, md_ener=md_ener, error=error)
    CALL release_md_ener(md_ener, error=error)

    ! If requested setup Thermostats
    CALL create_thermostats(thermostats, motion_section, force_env,&
         simpar, para_env, globenv, global_section, error )

    ! If requested setup Barostat
    CALL create_barostat_type(barostat, motion_section, force_env,&
         simpar, para_env, globenv, error )

    CALL set_md_env(md_env, thermostats=thermostats, barostat=barostat, error=error)

    ! Initialize or possibly restart KG coefficients (after thermostats)
    CALL force_env_get(force_env, dyn_coeff_set=dyn_coeff_set, error=error )
    IF (ASSOCIATED(dyn_coeff_set)) THEN
       kg_coeff_section  => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF",error=error)
       CALL dyn_coeff_set_initialize(dyn_coeff_set, motion_section, para_env,&
            globenv, coeff_section=kg_coeff_section, error=error)
    END IF

    !If requested set up the REFTRAJ run
    IF(simpar%ensemble == reftraj_ensemble) THEN
        reftraj_section => section_vals_get_subs_vals(motion_section,"MD%REFTRAJ",error=error)
        CALL create_reftraj(reftraj, reftraj_section, error=error)
        CALL set_md_env(md_env, reftraj=reftraj, error=error)
        CALL release_reftraj(reftraj,error=error)
    END IF

    
    CALL force_env_get(force_env, subsys=subsys, cell=cell, force_env_section=force_env_section,&
         error=error )

    ! Set V0 if needed
    IF ( ASSOCIATED ( cell_ref ) ) THEN
       IF ( simpar % v0 == 0._dp ) simpar % v0 = cell_ref % deth
    ENDIF

    ! Initialize velocities possibly applying constraints at the zeroth MD step
    CALL setup_velocities(force_env, simpar, globenv, motion_section, error)

    ! Setup Free Energy Calculation (if required)
    CALL fe_env_create (fe_env, motion_section, error)

    CALL set_md_env(md_env=md_env, simpar=simpar, fe_env=fe_env, cell=cell,&
                    force_env=force_env, error=error)
   
    ! Possibly initialize Wiener processes
    IF (simpar%ensemble == langevin_ensemble) CALL create_wiener_process(md_env,error)
    time_iter_start=m_walltime()

    CALL get_md_env ( md_env, force_env=force_env, itimes = itimes, &
                     md_ener=md_ener, t = t, used_time = used_time, error=error )

    ! Initialize the md_ener structure
    CALL initialize_md_ener(md_ener, force_env, simpar, error=error)

    ! Computing Forces at zero MD step
    IF (simpar%ensemble /= reftraj_ensemble) THEN
       IF(save_mem) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)
          work_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)
          work_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)
       END IF
       
       CALL force_env_calc_energy_force ( force_env, calc_force=.TRUE.,error=error)
       
       IF (simpar%do_respa)THEN
          CALL force_env_calc_energy_force ( force_env%sub_force_env(1)%force_env,&
               calc_force=.TRUE.,error=error)
       END IF
       
       CALL force_env_get( force_env, subsys=subsys, error=error )
      
       CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
            particles=particles,error=error)
      
       CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles,&
            force_env%virial, force_env%para_env%group, error=error)

       CALL section_vals_val_get(motion_section,"MD%STEP_START_VAL",i_val=itimes,error=error)
       CALL section_vals_val_get(motion_section,"MD%TIME_START_VAL",r_val=t,error=error)
       CALL md_energy(md_env,md_ener,error) 
       CALL md_write_output ( md_env, error ) !inits the print env at itimes == 0 also writes trajectories
    ELSE
       CALL get_md_env(md_env, reftraj=reftraj, error=error)
       CALL initialize_reftraj(reftraj, reftraj_section, md_env, error=error)
       itimes = reftraj%info%first_snapshot -1
    END IF
    
    CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
         root_section,"MOTION%CONSTRAINT%CONSTRAINT_INFO",error=error)
    CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
         root_section,"MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",error=error)

    ! Real MD Loop
    DO istep=1,simpar%nsteps
      itimes = itimes + 1

      CALL cp_iterate(logger%iter_info,last=(istep==simpar%nsteps),&
           iter_nr=itimes,error=error)
      ! Open possible Shake output units
      simpar%info_constraint = cp_print_key_unit_nr(&
           logger,root_section,"MOTION%CONSTRAINT%CONSTRAINT_INFO",&
           extension=".shakeLog",log_filename=.FALSE.,error=error)
      simpar%lagrange_multipliers = cp_print_key_unit_nr(&
           logger,root_section,"MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",&
           extension=".LagrangeMultLog",log_filename=.FALSE.,error=error)
      simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info,root_section,&
           "MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",error=error),cp_p_file)      

      ! Velocity Verlet Integrator
      CALL velocity_verlet(md_env,error)

      ! Close Shake output if requested...
      CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
           root_section,"MOTION%CONSTRAINT%CONSTRAINT_INFO",error=error)
      CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
           root_section,"MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",error=error)

      ! Free Energy calculation
      CALL free_energy_evaluate(md_env,should_stop,motion_section,error)

      IF (should_stop) EXIT
      ! Test for EXIT_MD or for WALL_TIME to exit 
      ! IF so we don't overwrite the restart or append to the trajectory
      ! because the execution could in principle stop inside the SCF where energy
      ! and forces are not converged. So we don't want to dump trajectory
      ! with wrong coordinates (due to the last integration step with fake forces)..
      CALL external_control(should_stop,"MD",globenv,error=error)
      IF (should_stop) EXIT

      CALL md_energy(md_env,md_ener,error) 
      CALL temperature_control(simpar,md_ener,force_env,logger,error)
      CALL comvel_control(simpar,md_ener,force_env,logger,error)

      time_iter_stop=m_walltime()
      used_time = time_iter_stop - time_iter_start
      time_iter_start=time_iter_stop

      CALL md_output(md_env,force_env_section,root_section,error)
    END DO

    CALL cp_rm_iter_level(logger%iter_info,n_rlevel_att=iter_level,&
         error=error)

    ! Deallocate Thermostats and Barostats
    CALL release_thermostats(thermostats, error=error)
    CALL release_barostat_type(barostat, error=error)

    DEALLOCATE(simpar,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_mol_dyn_low

END MODULE md_run
