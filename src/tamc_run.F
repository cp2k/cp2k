!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Perform a temperature accelarated hybrid monte carlo (TAHMC) run using QUICKSTEP
!> \par History
!>      none
!> \author Alin M Elena
! *****************************************************************************
MODULE tamc_run
!!!!!!!!!!
!!!!!!!!!! integrators
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE averages_types,                  ONLY: average_quantities_type
  USE barostat_types,                  ONLY: barostat_type,&
                                             create_barostat_type,&
                                             release_barostat_type
  USE bibliography,                    ONLY: VandenCic2006
  USE cell_types,                      ONLY: cell_type
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE colvar_types,                    ONLY: HBP_colvar_id,&
                                             WC_colvar_id,&
                                             colvar_p_type
  USE constraint_fxd,                  ONLY: fix_atom_control
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
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE free_energy_types,               ONLY: fe_env_create,&
                                             free_energy_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: &
       ehrenfest, langevin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, reftraj_ensemble
  USE input_cp2k_check,                ONLY: remove_restart_info
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE mc_environment_types,            ONLY: get_mc_env,&
                                             mc_env_create,&
                                             mc_env_release,&
                                             mc_environment_type,&
                                             set_mc_env
  USE mc_misc,                         ONLY: mc_averages_create,&
                                             mc_averages_release
  USE mc_move_control,                 ONLY: init_mc_moves,&
                                             mc_moves_release
  USE mc_types,                        ONLY: get_mc_par,&
                                             mc_averages_type,&
                                             mc_ekin_type,&
                                             mc_moves_type,&
                                             mc_simpar_type,&
                                             set_mc_par
  USE md_ener_types,                   ONLY: create_md_ener,&
                                             md_ener_type,&
                                             release_md_ener
  USE md_energies,                     ONLY: initialize_md_ener,&
                                             md_energy
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             set_md_env
  USE md_run,                          ONLY: qs_mol_dyn
  USE metadynamics_types,              ONLY: meta_env_type,&
                                             metavar_type,&
                                             set_meta_env
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE parallel_rng_types,              ONLY: UNIFORM,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: boltzmann,&
                                             femtoseconds,&
                                             joule,&
                                             kelvin
  USE qmmm_util,                       ONLY: apply_qmmm_walls_reflective
  USE qs_environment_types,            ONLY: get_qs_env
  USE qs_scf_post_gpw,                 ONLY: scf_post_calculation_gpw
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE reference_manager,               ONLY: cite_reference
  USE reftraj_types,                   ONLY: create_reftraj,&
                                             reftraj_type,&
                                             release_reftraj
  USE reftraj_util,                    ONLY: initialize_reftraj
  USE rt_propagation,                  ONLY: rt_prop_setup,&
                                             rt_write_input_restart
  USE rt_propagation_output,           ONLY: rt_prop_output
  USE simpar_methods,                  ONLY: read_md_section
  USE simpar_types,                    ONLY: create_simpar_type,&
                                             release_simpar_type,&
                                             simpar_type
  USE string_utilities,                ONLY: str_comp
  USE thermal_region_types,            ONLY: release_thermal_regions,&
                                             thermal_regions_type
  USE thermal_region_utils,            ONLY: create_thermal_regions
  USE thermostat_methods,              ONLY: create_thermostats
  USE thermostat_types,                ONLY: release_thermostats,&
                                             thermostats_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_methods,                  ONLY: virial_evaluate
  USE virial_types,                    ONLY: virial_type
  USE wiener_process,                  ONLY: create_wiener_process,&
                                             create_wiener_process_cv
!!!!! monte carlo part
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tamc_run'

  PUBLIC :: qs_tamc

CONTAINS

! *****************************************************************************
!> \brief Driver routine for TAHMC
!> \details it computes the forces using QuickStep.
!> \author Alin M Elena
! *****************************************************************************
SUBROUTINE qs_tamc(force_env, globenv,input_file_name,averages,error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    CHARACTER(LEN=*), OPTIONAL               :: input_file_name
    TYPE(average_quantities_type), &
      OPTIONAL, POINTER                      :: averages
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_tamc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=20)                        :: ensemble
    INTEGER :: handle, i, initialStep, iprint, isos, istep, j, md_stride, &
      nmccycles, output_unit, rand2skip, run_type_id
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: check, ehrenfest_md, &
                                                explicit, failure, ionode, &
                                                my_rm_restart_info, save_mem, &
                                                should_stop
    REAL(KIND=dp)                            :: auxRandom, inittime, rval, &
                                                temp, time_iter_start, &
                                                time_iter_stop
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: An, fz, xieta, zbuff
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    REAL(KIND=dp), POINTER                   :: constant, time, used_time
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(mc_averages_type), POINTER          :: MCaverages
    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(mc_moves_type), POINTER             :: gmoves, moves
    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(meta_env_type), POINTER             :: meta_env_saved
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(reftraj_type), POINTER              :: reftraj
    TYPE(rng_stream_type), POINTER           :: rng_stream, rng_stream_mc
    TYPE(section_vals_type), POINTER :: constraint_section, &
      force_env_section, free_energy_section, fs_section, global_section, &
      mc_section, md_section, motion_section, reftraj_section, &
      subsys_section, work_section
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(thermostats_type), POINTER          :: thermostats

    initialStep=0
    inittime=0.0_dp

    CALL timeset(routineN,handle)
    failure = .FALSE.
    my_rm_restart_info = .TRUE.
    NULLIFY(md_env, para_env,fs_section)
    para_env       => force_env%para_env
    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    md_section     => section_vals_get_subs_vals(motion_section,"MD",error=error)

    ! Real call to MD driver - Low Level
    CALL md_env_create(md_env, md_section, para_env, force_env=force_env, error=error)
    IF (PRESENT(averages)) CALL set_md_env(md_env, averages=averages, error=error)



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
    CALL cp_iterate(logger%iter_info,iter_nr=initialStep,error=error)
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
! ! !     CALL section_vals_val_get(motion_section,"PRINT%RESTART%SPLIT_RESTART_FILE",&
! ! !                               l_val=write_binary_restart_file,error=error)
!! let us see if this created all the trouble
!      CALL setup_velocities(force_env,simpar,globenv,md_env,md_section,constraint_section, &
!                            write_binary_restart_file,error)

    ! Setup Free Energy Calculation (if required)
    CALL fe_env_create (fe_env, free_energy_section, error)
    CALL set_md_env(md_env=md_env, simpar=simpar, fe_env=fe_env, cell=cell,&
                    force_env=force_env, error=error)

    ! Possibly initialize Wiener processes
    IF (simpar%ensemble == langevin_ensemble) CALL create_wiener_process(md_env,error)
    time_iter_start=m_walltime()

    CALL get_md_env(md_env, force_env=force_env, itimes=itimes, constant=constant,&
                    md_ener=md_ener, t=time, used_time=used_time, error=error)

    ! Attach the time counter of the meta_env to the one of the MD
    CALL set_meta_env(force_env%meta_env, time=time, error=error)
    ! Initialize the md_ener structure

    force_env%meta_env%dt=force_env%meta_env%zdt
    CALL initialize_md_ener(md_ener, force_env, simpar, error=error)
!     force_env%meta_env%dt=force_env%meta_env%zdt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MC setup up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NULLIFY(mc_env,mc_par,rng_stream_mc,MCaverages)

      CALL section_vals_get(force_env_section,n_repetition=isos,error=error)
      CPPostconditionNoFail(isos==1,cp_warning_level,routineP,error)
! set some values...will use get_globenv if that ever comes around

! initialize the random numbers
!       IF (para_env%ionode) THEN
        CALL create_rng_stream(rng_stream=rng_stream_mc,&
                             name="Random numbers for monte carlo acc/rej",&
                             distribution_type=UNIFORM,error=error)
!       ENDIF
!!!!! this shoudl go in a routine hmc_read

      NULLIFY(mc_section)
      ALLOCATE(mc_par, stat=isos) ! do not forget to clean this mess
      CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)

      mc_section => section_vals_get_subs_vals(force_env%root_section,&
         "MOTION%MC",error=error)
      CALL section_vals_val_get(mc_section,"ENSEMBLE",&
           c_val=ensemble,error=error)
      CPPostcondition(str_comp(ensemble,"TRADITIONAL"),cp_failure_level,routineP,error,failure)
      CALL section_vals_val_get(mc_section,"NSTEP",&
           i_val=nmccycles,error=error)
      CPPostcondition(nmccycles>0,cp_failure_level,routineP,error,failure)
      CALL section_vals_val_get(mc_section,"IPRINT",&
           i_val=iprint,error=error)
      CALL section_vals_val_get(mc_section,"RANDOMTOSKIP",i_val=rand2skip,error=error)
      CPPostcondition(rand2skip>=0,cp_failure_level,routineP,error,failure)
      temp=cp_unit_from_cp2k(simpar%temp_ext,"K",error=error)
!

      CALL set_mc_par(mc_par, ensemble=ensemble, nstep=nmccycles, iprint=iprint, temperature=temp, &
        beta=1.0_dp / temp / boltzmann * joule,exp_max_val = 0.9_dp*LOG(HUGE(0.0_dp)),&
        exp_min_val = 0.9_dp*LOG(TINY(0.0_dp)),max_val=HUGE(0.0_dp),min_val= 0.0_dp, &
        source=para_env%source, group=para_env%group, ionode=para_env%ionode, rand2skip=rand2skip)

      output_unit = cp_logger_get_default_io_unit(logger)
      IF (output_unit > 0) THEN
        WRITE(output_unit,'(a,a)')"HMC| Hybrid Monte Carlo Scheme "
        WRITE(output_unit,'(a,a)')"HMC| Ensemble ", ADJUSTL(ensemble)
        WRITE(output_unit,'(a,i0)')"HMC| MC Cycles ", nmccycles
        WRITE(output_unit,'(a,i0,a)')"HMC| Print every ", iprint, " cycles"
        WRITE(output_unit,'(a,i0)')"HMC| Number of random numbers to skip ", rand2skip
        WRITE(output_unit,'(a,f16.8,a)')"HMC| Temperature ",temp, "K"
      ENDIF

       CALL force_env_get( force_env, subsys=subsys, error=error )

       CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
            particles=particles,error=error)

        DO i=1,rand2skip
          auxRandom=next_random_number(rng_stream_mc,error=error)
          DO j=1,3*SIZE(particles%els)
            auxRandom=next_random_number(globenv%gaussian_rng_stream,error=error)
          ENDDO
        ENDDO



      CALL mc_env_create ( mc_env, error = error )
      CALL set_mc_env( mc_env,mc_par=mc_par,force_env=force_env)
!!!!!!!end mc setup

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
       CALL section_vals_val_set(md_section,"STEP_START_VAL",i_val=initialStep,error=error)
       CALL section_vals_val_set(md_section,"TIME_START_VAL",r_val=inittime,error=error)
       initialStep=itimes
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
!         CALL force_env_calc_energy_force (force_env, calc_force=.TRUE., error=error)
         meta_env_saved=> force_env%meta_env
         NULLIFY(force_env%meta_env)
         CALL force_env_calc_energy_force (force_env, calc_force=.FALSE., error=error)
         force_env%meta_env=>meta_env_saved
       END IF

       IF(ASSOCIATED(force_env%qs_env))THEN
!           force_env%qs_env%sim_time=time
!           force_env%qs_env%sim_step=itimes
          force_env%qs_env%sim_time=0.0_dp
          force_env%qs_env%sim_step=0
       END IF
       ! Warm-up engines for metadynamics
       IF (ASSOCIATED(force_env%meta_env)) THEN
          IF(force_env%meta_env%langevin) THEN
             CALL create_wiener_process_cv(force_env%meta_env, error=error)
             NULLIFY(rng_stream)
             DO j=1, (rand2skip-1)/nmccycles
               DO i = 1 , force_env%meta_env%n_colvar
                 rng_stream => force_env%meta_env%rng(i)%stream
                 auxRandom=next_random_number(rng_stream,error=error)
                 auxRandom=next_random_number(rng_stream,error=error)
               ENDDO
             ENDDO
          ENDIF
!           IF (force_env%meta_env%well_tempered) THEN
!              force_env%meta_env%wttemperature = simpar%temp_ext
!              IF (force_env%meta_env%wtgamma>EPSILON(1._dp)) THEN
!                 dummy=force_env%meta_env%wttemperature*(force_env%meta_env%wtgamma-1._dp)
!                 IF (force_env%meta_env%delta_t>EPSILON(1._dp)) THEN
!                    check=ABS(force_env%meta_env%delta_t-dummy)<1.E+3_dp*EPSILON(1._dp)
!                    CALL cp_assert(check,cp_failure_level,cp_assertion_failed,routineP,&
!                       "Inconsistency between DELTA_T and WTGAMMA (both specified):"//&
!                       " please, verify that DELTA_T=(WTGAMMA-1)*TEMPERATURE",&
!                       error,failure)
!                 ELSE
!                    force_env%meta_env%delta_t = dummy
!                 ENDIF
!              ELSE
!                 force_env%meta_env%wtgamma    = 1._dp &
!                    + force_env%meta_env%delta_t/force_env%meta_env%wttemperature
!              ENDIF
!              force_env%meta_env%invdt         = 1._dp/force_env%meta_env%delta_t
!           ENDIF
            CALL tamc_force(force_env,error=error)
!           CALL metadyn_write_colvar(force_env,error=error)
       ENDIF


       IF (simpar%do_respa)THEN
          CALL force_env_calc_energy_force (force_env%sub_force_env(1)%force_env,&
               calc_force=.TRUE.,error=error)
       END IF

!        CALL force_env_get( force_env, subsys=subsys, error=error )
!
!        CALL cp_subsys_get(subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
!             particles=particles,error=error)

       CALL virial_evaluate(atomic_kinds%els, particles%els, local_particles,&
            force_env%virial, force_env%para_env%group, error=error)

       CALL md_energy(md_env,md_ener,error)
!        CALL md_write_output(md_env, error) !inits the print env at itimes == 0 also writes trajectories
       md_stride = 1
    ELSE
       CALL get_md_env(md_env, reftraj=reftraj, error=error)
       CALL initialize_reftraj(reftraj, reftraj_section, md_env, error=error)
       itimes = reftraj%info%first_snapshot -1
       md_stride = reftraj%info%stride
    END IF

    CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
         constraint_section,"CONSTRAINT_INFO",error=error)
    CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
         constraint_section,"LAGRANGE_MULTIPLIERS",error=error)
      CALL init_mc_moves(moves)
      CALL init_mc_moves(gmoves)
      ALLOCATE (r(1:3,SIZE(particles%els)),STAT=isos)
      CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
!       ALLOCATE (r_old(1:3,size(particles%els)),STAT=isos)
!       CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
      CALL mc_averages_create(MCaverages)
    !!!!! some more buffers
 ! Allocate random number for Langevin Thermostat acting on COLVARS
      ALLOCATE (xieta(2*force_env%meta_env%n_colvar),STAT=isos)
      CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
      xieta(:) = 0.0_dp
      ALLOCATE (An(force_env%meta_env%n_colvar),STAT=isos)
      CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
      An(:) = 0.0_dp
      ALLOCATE (fz(force_env%meta_env%n_colvar),STAT=isos)
      CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
      fz(:) = 0.0_dp
      ALLOCATE (zbuff(2*force_env%meta_env%n_colvar),STAT=isos)
      CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
      zbuff(:) = 0.0_dp

    IF (output_unit>0)THEN
         WRITE(output_unit,'(a)')"HMC|==== Initial average forces"
    ENDIF
    CALL metadyn_write_colvar_header(force_env,error=error)
    moves%hmc%attempts=0
    moves%hmc%successes=0
    gmoves%hmc%attempts=0
    gmoves%hmc%successes=0
    IF (initialStep==0) THEN
!!! if we come from a restart we shall properly compute the average force
!!!      read the average force up to now
      DO i=1,force_env%meta_env%n_colvar
          fs_section => section_vals_get_subs_vals(force_env%meta_env%metadyn_section,"EXT_LAGRANGE_FS",error=error)
          CALL section_vals_get(fs_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(fs_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i, r_val=rval, error=error)
             fz(i) = rval*rand2skip
          END IF
      ENDDO

      CALL HMCsampler(globenv,force_env,MCaverages,r,mc_par,simpar,moves,gmoves,rng_stream_mc,output_unit,&
          fz,zbuff,nskip=rand2skip,logger=logger,iter=0,mcsec=mc_section,mdenv=md_env,error=error)
      CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=0,error=error)
      CALL section_vals_val_set(mc_section,"RANDOMTOSKIP",i_val=rand2skip+nmccycles,error=error)
      CALL write_restart(md_env=md_env,root_section=force_env%root_section, error=error)
    ENDIF
    IF (output_unit>0)THEN
         WRITE(output_unit,'(a)')"HMC|==== end initial average forces"
    ENDIF
!     call set_md_env(md_env, init=.FALSE., error=error)

    CALL metadyn_write_colvar(force_env,error=error)

    DO istep=1, force_env%meta_env%TAMCSteps
       ! Increase counters
       itimes = itimes + 1
       time   = time + force_env%meta_env%dt
       IF (output_unit>0)THEN
         WRITE(output_unit,'(a)')"HMC|==================================="
         WRITE(output_unit,'(a,1x,i0)')"HMC| on z step ", istep
       ENDIF
       !needed when electric field fields are applied
       IF(ASSOCIATED(force_env%qs_env))THEN
          force_env%qs_env%sim_time=time
          force_env%qs_env%sim_step=itimes
          force_env%meta_env%time=force_env%qs_env%sim_time
       END IF

       IF(ehrenfest_md)force_env%qs_env%rtp%istep=istep

       CALL cp_iterate(logger%iter_info,last=(istep==force_env%meta_env%TAMCSteps),iter_nr=itimes,error=error)
       ! Open possible Shake output units
       simpar%info_constraint = cp_print_key_unit_nr(logger,constraint_section,"CONSTRAINT_INFO",&
            extension=".shakeLog",log_filename=.FALSE.,error=error)
       simpar%lagrange_multipliers = cp_print_key_unit_nr(logger,constraint_section,&
            "LAGRANGE_MULTIPLIERS",extension=".LagrangeMultLog",log_filename=.FALSE.,error=error)
       simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info,constraint_section,&
            "LAGRANGE_MULTIPLIERS",error=error),cp_p_file)

       ! Velocity Verlet Integrator

       moves%hmc%attempts=0
       moves%hmc%successes=0
       CALL langevinVEC(md_env,globenv,mc_env,moves,gmoves,r,&
          rng_stream_mc,xieta,An,fz,MCaverages,zbuff,error)

       ! Close Shake output if requested...
       CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
            constraint_section,"CONSTRAINT_INFO",error=error)
       CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
            constraint_section,"LAGRANGE_MULTIPLIERS",error=error)
       CALL cp_iterate(logger%iter_info,iter_nr=initialStep,error=error)
       CALL metadyn_write_colvar(force_env,error=error)
       ! Free Energy calculation
!        CALL free_energy_evaluate(md_env,should_stop,free_energy_section,error)

       ![AME:UB] IF (should_stop) EXIT

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
       IF (should_stop) THEN
          CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=itimes,error=error)
!           CALL md_output(md_env,md_section,force_env%root_section,should_stop,error=error)
          IF(ehrenfest_md)THEN
             CALL rt_prop_output(force_env%qs_env,ehrenfest,error=error)
             CALL rt_write_input_restart(md_env,force_env,error)
          END IF
          EXIT
       END IF

!        IF(simpar%ensemble /= reftraj_ensemble) THEN
!           CALL md_energy(md_env, md_ener, error)
!           CALL temperature_control(simpar, md_env, md_ener, force_env, logger, error)
!           CALL comvel_control(md_ener, force_env, md_section, logger, error)
!           CALL angvel_control(md_ener, force_env, md_section, logger, error)
!        ELSE
!           CALL md_ener_reftraj(md_env, md_ener, error)
!        END IF

       time_iter_stop=m_walltime()
       used_time = time_iter_stop - time_iter_start
       time_iter_start=time_iter_stop

!!!!! this writes the restart...
!         CALL md_output(md_env,md_section,force_env%root_section,should_stop,error=error)

!        IF(simpar%ensemble == reftraj_ensemble ) THEN
!           CALL write_output_reftraj(md_env,error=error)
!        END IF

       IF (output_unit>0)THEN
         WRITE(output_unit,'(a,1x,i0)')"HMC| end z step ", istep
         WRITE(output_unit,'(a)')"HMC|==================================="
       ENDIF
    END DO
    CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=itimes,error=error)
    force_env%qs_env%sim_time=0.0_dp
    force_env%qs_env%sim_step=0
    rand2skip=rand2skip+nmccycles*force_env%meta_env%TAMCSteps
    IF (initialStep == 0) rand2skip=rand2skip+nmccycles
    CALL section_vals_val_set(mc_section,"RANDOMTOSKIP",i_val=rand2skip,error=error)

    CALL write_restart(md_env=md_env,root_section=force_env%root_section, error=error)
! if we need the final kinetic energy for Hybrid Monte Carlo
!     hmc_ekin%final_ekin=md_ener%ekin

    ! Remove the iteration level
    CALL cp_rm_iter_level(logger%iter_info,"MD",error=error)

    ! Deallocate Thermostats and Barostats
    CALL release_thermostats(thermostats, error=error)
    CALL release_barostat_type(barostat, error=error)
    CALL release_simpar_type(simpar, error)
    CALL release_thermal_regions(thermal_regions, error)

    CALL md_env_release(md_env, error=error)
    ! Clean restartable sections..
    IF (my_rm_restart_info) CALL remove_restart_info(force_env%root_section,error=error)
!     IF (para_env%ionode) THEN
      CALL delete_rng_stream(rng_stream_mc,error=error)
!     ENDIF
    CALL MC_ENV_RELEASE(mc_env,error)
    DEALLOCATE(mc_par, stat=isos) ! do not forget to clean this mess
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    CALL MC_MOVES_RELEASE(moves)
    CALL MC_MOVES_RELEASE(gmoves)
    DEALLOCATE(r,STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
!     DEALLOCATE(r_old,STAT=isos)
!     CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(xieta,STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(An,STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(fz,STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(zbuff,STAT=isos)
    CPPostcondition(isos==0,cp_failure_level,routineP,error,failure)
    CALL mc_averages_release(MCaverages)
    CALL timestop(handle)

  END SUBROUTINE qs_tamc



! *****************************************************************************
!> \brief Propagates velocities for z half a step
!> \details   Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \author Alin M Elena
! *****************************************************************************
  SUBROUTINE tamc_velocities_colvar(force_env,An,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: An
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tamc_velocities_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dt, fft, sigma
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    ! Add citation
    IF (meta_env%langevin) CALL cite_reference(VandenCic2006)
    dt = meta_env%dt

    ! Evolve Velocities
    meta_env%epot_walls = 0.0_dp
    DO i_c= 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       fft    = cv%ff_s+cv%ff_hills
       sigma  = SQRT((meta_env%temp_wanted*kelvin)*2.0_dp*(boltzmann/joule)*cv%gamma/cv%mass)
       cv%vvp=cv%vvp+0.5_dp*dt*(fft/cv%mass-cv%gamma*cv%vvp)*(1.0_dp-0.25_dp*dt*cv%gamma)+An(i_c)
       meta_env%epot_walls = meta_env%epot_walls + cv%epot_walls
    ENDDO
    CALL timestop(handle)
  END SUBROUTINE tamc_velocities_colvar

! *****************************************************************************
!> \brief propagates z one step
!> \details  Vanden-Eijnden Ciccotti C.Phys.Letter 429 (2006) 310-316
!> \author Alin M Elena
! *****************************************************************************
  SUBROUTINE tamc_position_colvar(force_env,xieta,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: xieta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tamc_position_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dt, sigma
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
!     IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)

    ! Add citation
    IF (meta_env%langevin) CALL cite_reference(VandenCic2006)
    dt = meta_env%dt

    ! Update of ss0
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       sigma  = SQRT((meta_env%temp_wanted*kelvin)*2.0_dp*(boltzmann/joule)*cv%gamma/cv%mass)
!        cv%ss0 =cv%ss0 +dt*cv%vvp
       cv%ss0 =cv%ss0 +dt*cv%vvp+dt*SQRT(dt/12.0_dp)*sigma*xieta(i_c+meta_env%n_colvar)
       IF (cv%periodic) THEN
          ! A periodic COLVAR is always within [-pi,pi]
          cv%ss0 = SIGN(1.0_dp,ASIN(SIN(cv%ss0)))*ACOS(COS(cv%ss0))
       END IF
    ENDDO
    CALL timestop(handle)

  END SUBROUTINE tamc_position_colvar

! *****************************************************************************
!> \brief Computes forces on z
!> #details also can be used to get the potenzial evergy of z
!> \author Alin M Elena
! *****************************************************************************
  SUBROUTINE tamc_force(force_env,zpot,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), INTENT(inout), OPTIONAL   :: zpot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tamc_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, icolvar, ii
    LOGICAL                                  :: explicit, failure
    REAL(kind=dp)                            :: diff_ss, dt, rval
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ss0_section, ss_section, &
                                                vvp_section

    failure=.FALSE.
    NULLIFY(logger,meta_env)
    meta_env => force_env%meta_env
!     IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    NULLIFY(colvar_p,subsys,cv,ss0_section, vvp_section,ss_section)
    CALL force_env_get(force_env, subsys=subsys, error=error)

    dt = meta_env%dt
    IF (.NOT.meta_env%restart) meta_env%n_steps=meta_env%n_steps+1
    ! compute ss and the derivative of ss with respect to the atomic positions
    DO i_c=1,meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       icolvar = cv%icolvar
       CALL colvar_eval_glob_f(icolvar,force_env,error=error)
       cv%ss       = subsys%colvar_p(icolvar)%colvar%ss
       ! Restart for Extended Lagrangian Metadynamics
       IF (meta_env%restart) THEN
          ! Initialize the position of the collective variable in the extended lagrange
          ss0_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS0",error=error)
          CALL section_vals_get(ss0_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(ss0_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%ss0 = rval
          ELSE
             cv%ss0 = cv%ss
          END IF
          vvp_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_VVP",error=error)
          CALL section_vals_get(vvp_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL setup_velocities_z(force_env,error)
             CALL section_vals_val_get(vvp_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%vvp = rval
          ELSE
             CALL setup_velocities_z(force_env,error)
          ENDIF
          ss_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS",error=error)
          CALL section_vals_get(ss_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(ss_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%ss = rval
          END IF
       END IF
       !
    ENDDO
          ! forces on the atoms
       NULLIFY(particles)
       CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
          particles=particles,error=error)

       meta_env%restart = .FALSE.
       meta_env%epot_s = 0.0_dp
       meta_env%epot_walls = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          diff_ss = cv%ss-cv%ss0
          IF (cv%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          cv%epot_s = 0.5_dp*cv%lambda*diff_ss*diff_ss
          cv%ff_s   = cv%lambda*(diff_ss)
          meta_env%epot_s=meta_env%epot_s+cv%epot_s
          icolvar   = cv%icolvar

          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f- cv%ff_s*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO

       ENDDO
    IF(PRESENT(zpot))zpot=meta_env%epot_s
    CALL fix_atom_control(force_env, error=error)

    CALL timestop(handle)
  END SUBROUTINE tamc_force



! *****************************************************************************
!> \brief propagates one time step both z systems and samples the x system
!> \author Alin M Elena
! *****************************************************************************
  SUBROUTINE langevinVEC ( md_env, globenv,mc_env,moves,gmoves,r,&
              rng_stream_mc,xieta,An,fz,averages,zbuff,error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mc_environment_type), POINTER       :: mc_env
    TYPE(mc_moves_type), POINTER             :: moves, gmoves
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r
    TYPE(rng_stream_type), POINTER           :: rng_stream_mc
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: xieta, An, fz
    TYPE(mc_averages_type), INTENT(INOUT), &
      POINTER                                :: averages
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: zbuff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'langevinVEC', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iprint, ivar, nparticle, &
                                                nparticle_kind, nstep, &
                                                output_unit
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dt, gamma, mass, sigma
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(virial_type), POINTER               :: virial

    NULLIFY(logger, mc_par)
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)

    NULLIFY(rng_stream)
    failure = .FALSE.
! quantitites to be nullified for the get_md_env
    NULLIFY(simpar,force_env,para_env)
! quantities to be nullified for the force_env_get environment
    NULLIFY(subsys,cell)
! quantitites to be nullified for the cp_subsys_get
    NULLIFY(atomic_kinds,local_particles,particles,local_molecules,molecules,molecule_kinds,gci)

    CALL get_md_env(md_env=md_env, simpar=simpar, force_env=force_env,&
          para_env=para_env, error=error)
    CALL get_mc_env(mc_env, mc_par=mc_par)
    CALL get_mc_par(mc_par,nstep=nstep,iprint=iprint)

    dt = simpar%dt
    CALL force_env_get(force_env=force_env,subsys=subsys,cell=cell,&
         virial=virial,error=error)

!!!! this bit should vanish once I understand what the hell is with it

!     ! Do some checks on coordinates and box
     CALL apply_qmmm_walls_reflective(force_env, error=error)

     CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

     nparticle_kind  = atomic_kinds%n_els
     atomic_kind_set => atomic_kinds%els
     molecule_kind_set => molecule_kinds%els

     nparticle    = particles%n_els
     particle_set => particles%els
     molecule_set => molecules%els
     CPPrecondition(ASSOCIATED(force_env%meta_env),cp_failure_level,routineP,error,failure)
     CPPrecondition(force_env%meta_env%langevin,cp_failure_level,routineP,error,failure)
   !    *** Velocity Verlet for Langevin *** v(t)--> v(t+1/2)
         !!!!!! noise xi is in the first half, eta in the second half
      DO ivar = 1 , force_env%meta_env%n_colvar
          rng_stream => force_env%meta_env%rng(ivar)%stream
          xieta(ivar)=next_random_number(rng_stream,error=error)
          xieta(ivar+force_env%meta_env%n_colvar)=next_random_number(rng_stream,error=error)
          gamma=force_env%meta_env%metavar(ivar)%gamma
          mass=force_env%meta_env%metavar(ivar)%mass
          sigma  = SQRT((force_env%meta_env%temp_wanted*kelvin)*2.0_dp*(boltzmann/joule)*gamma/mass)
          An(ivar)=0.5_dp*SQRT(dt)*sigma*(xieta(ivar)*(1.0_dp-0.5_dp*dt*gamma)-&
            dt*gamma*xieta(ivar+force_env%meta_env%n_colvar)/SQRT(12.0_dp))
      ENDDO
!    *** Velocity Verlet for Langeving *** v(t)--> v(t+1/2)
      CALL tamc_velocities_colvar(force_env,An,error=error)
!    *** Velocity Verlet for Langevin S(t)->S(t+1)
      CALL tamc_position_colvar(force_env,xieta,error)
!!!!! start zHMC sampler
      CALL HMCsampler(globenv,force_env,averages,r,mc_par,simpar,moves,gmoves,rng_stream_mc,output_unit,fz,zbuff,error=error)

!     CALL final_mc_write(mc_par,tmp_moves,&
!                output_unit,energy_check,&
!                initial_energy,final_energy,&
!                averages)

!!!!!!!!!!!!!!!!!!!! end zHMC sampler
         !    *** Velocity Verlet for Langeving *** v(t+1/2)--> v(t+1)
      CALL tamc_velocities_colvar(force_env,An,error)
!       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
!          local_particles, virial, para_env%group, error=error)

  END SUBROUTINE langevinVEC


! *****************************************************************************
!> \brief Driver routin for the canonical sampler using modified HMC
!> \details at the end of this routine %ff_s will contain mean force
!> \author Alin M Elena
! *****************************************************************************

   SUBROUTINE HMCsampler(globenv,force_env,averages,r,mc_par,simpar,moves,gmoves,rng_stream_mc,output_unit,&
                    fz,zbuff,nskip,logger,mcsec,mdenv,iter,error)
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(mc_averages_type), POINTER          :: averages
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r
    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(mc_moves_type), POINTER             :: moves, gmoves
    TYPE(rng_stream_type), POINTER           :: rng_stream_mc
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: fz, zbuff
    INTEGER, INTENT(IN), OPTIONAL            :: nskip
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: mcsec
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: mdenv
    INTEGER, INTENT(IN), OPTIONAL            :: iter
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, iprint, ishift, it1, j, &
                                                nsamples, nstep
    REAL(KIND=dp)                            :: energy_check, old_epx, &
                                                old_epz, t1
    TYPE(meta_env_type), POINTER             :: meta_env_saved

     IF (PRESENT(nskip)) THEN
        nsamples=nskip
        ishift=nskip
     ELSE
        nsamples=0
        fz=0.0_dp
        ishift=0
     ENDIF
!      lrestart = .false.
!      if (present(logger) .and. present(iter)) THEN
!       lrestart=.true.
!      ENDIF
     CALL get_mc_par(mc_par,nstep=nstep,iprint=iprint)
     meta_env_saved=> force_env%meta_env
     NULLIFY(force_env%meta_env)
     CALL force_env_get(force_env,potential_energy=old_epx,error=error)
     force_env%meta_env=> meta_env_saved

     old_epz=force_env%meta_env%epot_s
!!! average energy will be wrong on restarts
     averages%ave_energy=0.0_dp
     t1=force_env%qs_env%sim_time
     it1=force_env%qs_env%sim_step
     IF (output_unit > 0) THEN
        WRITE(output_unit,'(a,l4)')"HMC|restart? ",force_env%meta_env%restart
        WRITE(output_unit,'(a,3(f16.8,1x))')"HMC|Ep, Epx, Epz ",old_epx+force_env%meta_env%epot_s,old_epx,force_env%meta_env%epot_s
        WRITE(output_unit,'(a)')"#HMC| No | z.. | theta.. | ff_z... | ff_z/n |"
     ENDIF
     DO i = 1,nstep
      IF(MOD(i,iprint) == 0 .AND. (output_unit>0)) THEN
        WRITE(output_unit,'(a,1x,i0)') "HMC|========== On Monte Carlo cycle ",i+ishift
        WRITE(output_unit,'(a)') "HMC| Attempting a minitrajectory move"
        WRITE(output_unit,'(a,1x,i0)')"HMC| start mini-trajectory", i+ishift
        WRITE(output_unit,'(a,1x,i0,1x)',advance="no")"#HMC|0 ",i+ishift
        DO j=1,force_env%meta_env%n_colvar
          WRITE(output_unit,'(f16.8,1x,f16.8,1x,f16.8)',advance="no")force_env%meta_env%metavar(j)%ss0,&
            force_env%meta_env%metavar(j)%ss,&
            force_env%meta_env%metavar(j)%ff_s!,fz(j)/real(i+ishift,dp)
        ENDDO
        WRITE(output_unit,*)
      ENDIF

      CALL mc_hmc_move(mc_par, force_env,globenv, moves,gmoves,old_epx,old_epz,energy_check,&
          r, output_unit, rng_stream_mc,zbuff,error)
     ! check averages...
     ! force average for z needed too...
        averages%ave_energy=averages%ave_energy*REAL(i-1,dp)/REAL(i,dp)+&
                  old_epx/REAL(i,dp)
        DO j=1,force_env%meta_env%n_colvar
          fz(j)=fz(j)+force_env%meta_env%metavar(j)%ff_s
        ENDDO
        IF (output_unit > 0) THEN
          WRITE(output_unit,'(a,1x,i0)')"HMC|end mini-trajectory",i+ishift
!!!!!!!! this prints z and theta(x) --ss0,ss-- needed to determine an acceptable k then the instanteneous force and some instanteneous average for force
          WRITE(output_unit,'(a,1x,i0,1x)',advance="no")"#HMC|1 ",i+ishift
          DO j=1,force_env%meta_env%n_colvar
            WRITE(output_unit,'(f16.8,1x,f16.8,1x,f16.8,1x,f16.8)',advance="no")force_env%meta_env%metavar(j)%ss0,&
              force_env%meta_env%metavar(j)%ss,&
              force_env%meta_env%metavar(j)%ff_s,fz(j)/REAL(i+ishift,dp)
          ENDDO
          WRITE(output_unit,*)
        ENDIF
        nsamples=nsamples+1
        IF(MOD(i,iprint) == 0 .AND. (output_unit>0)) THEN
          WRITE(output_unit,'(a,f16.8)')"HMC| Running average for potential energy ", averages%ave_energy
          WRITE(output_unit,'(a,1x,i0)') "HMC|======== End Monte Carlo cycle ",i+ishift
        ENDIF
!         IF (lrestart) THEN
!           k=nstep/5
!           IF(MOD(i,k) == 0) THEN
!              force_env%qs_env%sim_time=t1
!              force_env%qs_env%sim_step=it1
!              DO j=1,force_env%meta_env%n_colvar
!                force_env%meta_env%metavar(j)%ff_s=fz(j)/real(i+ishift,dp)
!              ENDDO
! !              CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=-1,error=error)
!              CALL section_vals_val_set(mcsec,"RANDOMTOSKIP",i_val=i+ishift,error=error)
!              CALL write_restart(md_env=mdenv,root_section=force_env%root_section, error=error)
! !              CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter,error=error)
!           ENDIF
!         ENDIF
      ENDDO
    force_env%qs_env%sim_time=t1
    force_env%qs_env%sim_step=it1
    IF (output_unit > 0) THEN
      WRITE(output_unit,'(a,i0,a,i0,a,f16.8)')"HMC| local acceptance ratio: ",moves%hmc%successes,"/" ,&
          moves%hmc%attempts, "=",REAL(moves%hmc%successes,dp)/ REAL(moves%hmc%attempts,dp)
      WRITE(output_unit,'(a,i0,a,i0,a,f16.8)')"HMC| global acceptance ratio: ",gmoves%hmc%successes,"/" ,&
          gmoves%hmc%attempts, "=",REAL(gmoves%hmc%successes,dp)/ REAL(gmoves%hmc%attempts,dp)
    ENDIF
    !average force
    DO j=1,force_env%meta_env%n_colvar
      force_env%meta_env%metavar(j)%ff_s=fz(j)/nsamples
    ENDDO
   END SUBROUTINE HMCsampler

! *****************************************************************************
!> \brief performs a hybrid Monte Carlo move
!> \details It runs a NVE trajectory, followed by localisation and accepts rejects
!> using the biased Hamiltonian, rather than the traditional guiding Hamiltonian
!> \author Alin M Elena
! *****************************************************************************
  SUBROUTINE mc_hmc_move ( mc_par,force_env, globenv, moves,gmoves,old_epx,old_epz,&
                        energy_check,r,output_unit,rng_stream,zbuff,error)

    TYPE(mc_simpar_type), POINTER            :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(mc_moves_type), POINTER             :: moves, gmoves
    REAL(KIND=dp), INTENT(INOUT)             :: old_epx, old_epz, energy_check
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: r
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(rng_stream_type), POINTER           :: rng_stream
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: zbuff
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mc_hmc_move', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, iatom, istat, &
                                                j, nAtoms, source
    LOGICAL                                  :: failure, ionode, localise
    REAL(KIND=dp)                            :: BETA, energy_term, &
                                                exp_max_val, exp_min_val, &
                                                new_energy, new_epx, new_epz, &
                                                rand, value, w
    TYPE(cp_subsys_type), POINTER            :: oldsys
    TYPE(mc_ekin_type), POINTER              :: hmc_ekin
    TYPE(meta_env_type), POINTER             :: meta_env_saved
    TYPE(particle_list_type), POINTER        :: particles_set
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(section_vals_type), POINTER         :: dft_section, input

! begin the timing of the subroutine

      CALL timeset(routineN,handle)
      failure=.TRUE.

      NULLIFY(scf_env)
      CALL get_qs_env(force_env%qs_env,scf_env=scf_env,error=error,input=input)
      dft_section => section_vals_get_subs_vals(input,"DFT",error=error)

! get a bunch of stuff from mc_par
      CALL get_mc_par(mc_par,ionode=ionode,&
         BETA=BETA,exp_max_val=exp_max_val,&
         exp_min_val=exp_min_val,source=source,group=group)

! nullify some pointers
!       NULLIFY(particles_set,oldsys,hmc_ekin)
      NULLIFY(particles_set,oldsys,meta_env_saved,hmc_ekin)
      ! now let's grab the particle positions
      CALL force_env_get(force_env,subsys=oldsys,&
         error=error)
      CALL cp_subsys_get(oldsys,particles=particles_set, &
         error=error)
      nAtoms=SIZE(particles_set%els)
! do some allocation


      ALLOCATE (hmc_ekin,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

! record the attempt
      moves%hmc%attempts=moves%hmc%attempts+1
      gmoves%hmc%attempts=gmoves%hmc%attempts+1

! save the old coordinates just in case we need to go back
      DO iatom=1,nAtoms
         r(1:3,iatom)=particles_set%els(iatom)%r(1:3)
      ENDDO
      localise=.TRUE.
! the same for collective variables data should be made,ss first half and ff_s the last half
      DO j=1,force_env%meta_env%n_colvar
        zbuff(j)=force_env%meta_env%metavar(j)%ss
        zbuff(j+force_env%meta_env%n_colvar)=force_env%meta_env%metavar(j)%ff_s
        IF ((oldsys%colvar_p(force_env%meta_env%metavar(j)%icolvar)%colvar%type_id == HBP_colvar_id) .OR.&
            (oldsys%colvar_p(force_env%meta_env%metavar(j)%icolvar)%colvar%type_id == WC_colvar_id)) THEN
            localise=.FALSE.
        ENDIF
      ENDDO

! now run the MD simulation
      meta_env_saved=> force_env%meta_env
      NULLIFY(force_env%meta_env)
      force_env%qs_env%sim_time=0.0_dp
      force_env%qs_env%sim_step=0
      IF (.NOT. localise) THEN
        CALL section_vals_val_set(dft_section,"LOCALIZE%_SECTION_PARAMETERS_",l_val=.FALSE.,error=error)
      ENDIF
      CALL qs_mol_dyn(force_env,globenv,error=error,hmc_ekin=hmc_ekin)
      IF (.NOT. localise) THEN
        CALL section_vals_val_set(dft_section,"LOCALIZE%_SECTION_PARAMETERS_",l_val=.TRUE.,error=error)
        CALL scf_post_calculation_gpw(dft_section, scf_env, force_env%qs_env, error)
      ENDIF

      CALL force_env_get(force_env, potential_energy=new_epx,error=error)

      force_env%meta_env=>meta_env_saved
      CALL tamc_force(force_env, zpot=new_epz,error=error)
      new_energy=new_epx+new_epz
      IF (output_unit>0) THEN
        WRITE(output_unit,'(a,4(f16.8,1x))')"HMC|old Ep, Ekx, Epz, Epx ",old_epx+old_epz,hmc_ekin%initial_ekin, old_epz,old_epx
        WRITE(output_unit,'(a,4(f16.8,1x))')"HMC|new Ep, Ekx, Epz, Epx ",new_energy,hmc_ekin%final_ekin,new_epz,new_epx
      ENDIF
      energy_term=new_energy-old_epx-old_epz+hmc_ekin%final_ekin-hmc_ekin%initial_ekin

      value=-BETA*(energy_term)
! to prevent overflows
      IF    (value > exp_max_val) THEN
         w=10.0_dp
      ELSEIF(value < exp_min_val) THEN
         w=0.0_dp
      ELSE
         w=EXP(value)
      ENDIF

      rand=next_random_number(rng_stream,error=error)
      IF (rand < w) THEN
! accept the move
         moves%hmc%successes=moves%hmc%successes+1
         gmoves%hmc%successes=gmoves%hmc%successes+1
! update energies
         energy_check=energy_check+(new_energy-old_epx-old_epz)
         old_epx=new_epx
         old_epz=new_epz
      ELSE
! reset the cell and particle positions
         DO iatom=1,nAtoms
            particles_set%els(iatom)%r(1:3)=r(1:3,iatom)
         ENDDO
         DO j=1,force_env%meta_env%n_colvar
          force_env%meta_env%metavar(j)%ss=zbuff(j)
          force_env%meta_env%metavar(j)%ff_s=zbuff(j+force_env%meta_env%n_colvar)
        ENDDO

      ENDIF

      DEALLOCATE(hmc_ekin,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

! end the timing
      CALL timestop(handle)

    END SUBROUTINE mc_hmc_move



SUBROUTINE metadyn_write_colvar_header(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_write_colvar_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=100)                       :: aux, fmt
    CHARACTER(len=255)                       :: label1, label2, label3, &
                                                label4, label5, label6
    INTEGER                                  :: handle, i, iw, m
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env

    failure=.FALSE.
    NULLIFY(logger,meta_env)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)

    iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
         "PRINT%COLVAR",extension=".metadynLog",error=error)
    IF (iw>0) THEN
          label1=""
          label2=""
          label3=""
          label4=""
          label5=""
          label6=""
          DO i=1,meta_env%n_colvar
            WRITE(aux,'(a,i0)')"z_",i
            label1=TRIM(label1)//TRIM(aux)
            m=15*i-LEN_TRIM(label1)-1
            label1=TRIM(label1)//REPEAT(" ",m)//"|"
            WRITE(aux,'(a,i0)')"Theta_",i
            label2=TRIM(label2)//TRIM(aux)
            m=15*i-LEN_TRIM(label2)-1
            label2=TRIM(label2)//REPEAT(" ",m)//"|"
            WRITE(aux,'(a,i0)')"F_z",i
            label3=TRIM(label3)//TRIM(aux)
            m=15*i-LEN_TRIM(label3)-1
            label3=TRIM(label3)//REPEAT(" ",m)//"|"
            WRITE(aux,'(a,i0)')"F_h",i
            label4=TRIM(label4)//TRIM(aux)
            m=15*i-LEN_TRIM(label4)-1
            label4=TRIM(label4)//REPEAT(" ",m)//"|"
            WRITE(aux,'(a,i0)')"F_w",i
            label5=TRIM(label5)//TRIM(aux)
            m=15*i-LEN_TRIM(label5)-1
            label5=TRIM(label5)//REPEAT(" ",m)//"|"
            WRITE(aux,'(a,i0)')"v_z",i
            label6=TRIM(label6)//TRIM(aux)
            m=15*i-LEN_TRIM(label6)-1
            label6=TRIM(label6)//REPEAT(" ",m)//"|"
          ENDDO
          WRITE (fmt, '("(a17,6a",i0 ,",4a15)")') meta_env%n_colvar*15
          WRITE(iw,TRIM(fmt))"#Time[fs] |", &
               TRIM(label1), &
               TRIM(label2), &
               TRIM(label3), &
               TRIM(label4), &
               TRIM(label5), &
               TRIM(label6), &
                "Epot_z |", &
                "Ene hills |", &
                "Epot walls |", &
               "Temperature |"

    END IF
    CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
         "PRINT%COLVAR", error=error)

    CALL timestop(handle)

  END SUBROUTINE metadyn_write_colvar_header

  SUBROUTINE metadyn_write_colvar(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_write_colvar', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: temp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    failure=.FALSE.
    NULLIFY(logger,meta_env,cv)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)

    IF (meta_env%langevin) THEN
       meta_env%ekin_s = 0.0_dp
!        meta_env%epot_s = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
       ENDDO
    END IF

    ! write COLVAR file
    iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
         "PRINT%COLVAR",extension=".metadynLog",error=error)
    IF (iw>0) THEN
       IF (meta_env%extended_lagrange) THEN
          WRITE(iw,'(f16.8,70f15.8)')meta_env%time*femtoseconds, &
               (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ss,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_s,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%ff_walls,i=1,meta_env%n_colvar), &
               (meta_env%metavar(i)%vvp,i=1,meta_env%n_colvar), &
                meta_env%epot_s, &
                meta_env%hills_env%energy, &
                meta_env%epot_walls, &
               (meta_env%ekin_s)*2.0_dp/(REAL(meta_env%n_colvar,KIND=dp))*kelvin
       ELSE
          WRITE(iw,'(f16.8,40f13.5)') meta_env%time*femtoseconds,&
               (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar),&
               (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar),&
               (meta_env%metavar(i)%ff_walls,i=1,meta_env%n_colvar),&
                meta_env%hills_env%energy,&
                meta_env%epot_walls
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
         "PRINT%COLVAR", error=error)
    ! Temperature for COLVAR
    IF (meta_env%extended_lagrange) THEN
       temp = meta_env%ekin_s*2.0_dp/(REAL(meta_env%n_colvar,KIND=dp))*kelvin
       meta_env%avg_temp = (meta_env%avg_temp*REAL(meta_env%n_steps,KIND=dp)+&
            temp)/REAL(meta_env%n_steps+1,KIND=dp)
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%TEMPERATURE_COLVAR",extension=".metadynLog",error=error)
       IF (iw > 0) THEN
          WRITE (iw, '(T2,79("-"))')
          WRITE (iw,'( A,T51,f10.2,T71,f10.2)' )' COLVARS INSTANTANEOUS/AVERAGE TEMPERATURE ',&
               temp, meta_env%avg_temp
          WRITE (iw, '(T2,79("-"))')
       ENDIF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%TEMPERATURE_COLVAR", error=error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE metadyn_write_colvar


  SUBROUTINE setup_velocities_z(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i_c
    REAL(kind=dp)                            :: ekin_w, fac_t
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv

    NULLIFY(meta_env)
    meta_env=>force_env%meta_env
    meta_env%ekin_s = 0.0_dp
    DO i_c=1,meta_env%n_colvar
      cv => meta_env%metavar(i_c)
      cv%vvp = next_random_number(force_env%globenv%gaussian_rng_stream,error=error)
      meta_env%ekin_s = meta_env%ekin_s + 0.5_dp*cv%mass*cv%vvp**2
    END DO
    ekin_w = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
    fac_t  = SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
    DO i_c = 1, meta_env%n_colvar
      cv     => meta_env%metavar(i_c)
      cv%vvp = cv%vvp*fac_t
    ENDDO
  END SUBROUTINE setup_velocities_z
END MODULE tamc_run
