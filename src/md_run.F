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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_initialize,&
                                             dyn_coeff_set_p_type
  USE extended_system_types,           ONLY: extended_system_type,&
                                             extended_type_init,&
                                             lnhc_dealloc,&
                                             lnhc_init,&
                                             npt_dealloc
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_rattle,&
                                             force_env_shake
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_type
  USE free_energy_methods,             ONLY: free_energy_evaluate
  USE free_energy_types,               ONLY: fe_env_create,&
                                             free_energy_type
  USE global_types,                    ONLY: global_environment_type
  USE init_extended_system_types,      ONLY: initialize_nhc_baro,&
                                             initialize_nhc_coef,&
                                             initialize_nhc_part,&
                                             initialize_npt
  USE init_extended_system_variables,  ONLY: init_nhc_forces
  USE input_constants,                 ONLY: langevin_ensemble
  USE input_cp2k_restarts,             ONLY: update_subsys,&
                                             write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE md,                              ONLY: initialize_velocities,&
                                             read_md_section,&
                                             simulation_parameters_type
  USE md_energies,                     ONLY: md_energy
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_env_create,&
                                             md_env_release,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_types_new,              ONLY: global_constraint_type
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: kelvin
  USE shell_potential_types,           ONLY: shell_kind_type
  USE termination,                     ONLY: external_control,&
                                             stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE velocity_verlet_control,         ONLY: velocity_verlet
  USE virial_methods,                  ONLY: virial_evaluate
  USE wiener_process,                  ONLY: create_wiener_process
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_run'

  PUBLIC :: qs_mol_dyn, &
            temperature_control


  INTERFACE qs_mol_dyn
     MODULE PROCEDURE qs_mol_dyn1, qs_mol_dyn2
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE qs_mol_dyn1 ( force_env, globenv, error )

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p

    ALLOCATE(force_env_p(1), stat=stat)
    CPPostconditionNoFail(stat==0,cp_failure_level,routineP,error)
    force_env_p(1)%force_env => force_env
    CALL qs_mol_dyn2 ( force_env_p, globenv, error )

  END SUBROUTINE qs_mol_dyn1

! *****************************************************************************

  SUBROUTINE qs_mol_dyn2 ( force_env_p, globenv, error )

!   Purpose: Driver routine for MD run using QUICKSTEP.

!   History: - Creation (07.11.2002,MK)

!   ***************************************************************************

    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_mol_dyn2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ibead, in_use, &
                                                istat, istep, isubsys, &
                                                iter_level, nbeads
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: apply_cns0, failure, &
                                                save_mem, shell_adiabatic, &
                                                shell_present, should_stop
    REAL(KIND=dp)                            :: time_iter_start, &
                                                time_iter_stop
    REAL(KIND=dp), POINTER                   :: t, used_time
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(extended_system_type)               :: extended_type
    TYPE(free_energy_type), POINTER          :: fe_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(section_vals_type), POINTER :: force_env_section, global_section, &
      motion_section, print_section, root_section, subsys_section, &
      work2_section, work3_section, work_section
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env_p),cp_failure_level,routineP,error,failure)

    failure=.FALSE.
    NULLIFY ( particles, shell_particles, core_particles, local_molecules,&
              cell, cell_ref, dyn_coeff_set,&
              simpar,itimes,used_time,subsys,motion_section,gci )
    NULLIFY (md_env,molecule_kinds_new, molecules_new, force_env_section,&
         work_section, work2_section, work3_section, print_section, atomic_kinds, local_particles, t,&
         fe_env)
    logger => cp_error_get_logger(error)
    root_section=>force_env_p ( 1 ) % force_env % root_section
    para_env=>force_env_p ( 1 ) % force_env % para_env

    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)

    NULLIFY(global_section)
    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)


    nbeads = SIZE ( force_env_p )
    DO i = 1, nbeads
      CPPrecondition ( ASSOCIATED ( force_env_p ( i ) % force_env ), cp_failure_level,routineP,error,failure)
    END DO

    ALLOCATE (simpar,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "simpar",0)
    force_env_section => force_env_p(1)%force_env%force_env_section
    subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    CALL cp_add_iter_level(logger%iter_info,"MD",n_rlevel_new=iter_level,&
         error=error)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=0,error=error)
    CALL read_md_section ( simpar, root_section, error )
    simpar%info_constraint = cp_print_key_unit_nr(logger,root_section,&
         "MOTION%CONSTRAINT%CONSTRAINT_INFO",&
         extension=".shakeLog",log_filename=.FALSE.,error=error)
    simpar%lagrange_multipliers = cp_print_key_unit_nr(logger,root_section,&
         "MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",&
         extension=".LagrangeMultLog",log_filename=.FALSE.,error=error)
    simpar%dump_lm = BTEST(cp_print_key_should_output(logger%iter_info,root_section,&
         "MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",error=error),cp_p_file)

    CALL extended_type_init(extended_type)
    ALLOCATE (extended_type%nhc_part(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (extended_type%nhc_baro(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (extended_type%nhc_coef(nbeads),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ibead=1,nbeads
       CALL lnhc_init(extended_type%nhc_part(ibead))
       CALL lnhc_init(extended_type%nhc_baro(ibead))
       CALL lnhc_init(extended_type%nhc_coef(ibead))
    ENDDO

    DO ibead = 1, nbeads
       CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
                          dyn_coeff_set=dyn_coeff_set,in_use=in_use,&
                          cell=cell, cell_ref = cell_ref,&
                          force_env_section=force_env_section, error=error )

       ! Define a reference cell for compatibility with MD restarts
       DO isubsys=1,SIZE(subsys)
          subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",i_rep_section=isubsys,&
               error=error)

          CALL cp_subsys_get(subsys(isubsys)%subsys,&
              local_molecules_new=local_molecules,&
              molecule_kinds_new=molecule_kinds_new,&
              atomic_kinds=atomic_kinds,&
              local_particles=local_particles,&
              molecules_new=molecules_new,&
              gci=gci,&
              particles=particles,&
              error=error)
          !
          ! Possibly restart Nose particles
          !
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(motion_section,"MD%NOSE",error=error)
          print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)
          CALL initialize_nhc_part(cell,simpar,local_molecules, &
                                   molecules_new%els, molecule_kinds_new%els,&
                                   para_env, globenv, extended_type%nhc_part(ibead), &
                                   print_section=print_section,nose_section=work_section,&
                                   gci=gci,error=error)
          !
          ! Possibly restart barostat
          !
          NULLIFY(work_section)
          work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%NOSE",error=error)
          CALL initialize_nhc_baro(simpar,para_env,globenv, extended_type%nhc_baro(ibead),&
                                   work_section=work_section,error=error)

          IF (ASSOCIATED(dyn_coeff_set)) THEN
             !
             ! Initialize coef thermostats and setup coefficient forces, velocities and effective masses
             ! Possibly restart KG coefficients
             !
             NULLIFY(work_section)
             work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%NOSE",error=error)
             CALL initialize_nhc_coef(dyn_coeff_set(isubsys)%dyn_coeff_set,&
                                   root_section,para_env,globenv,extended_type%nhc_coef(ibead), &
                                   work_section=work_section,error=error)

             NULLIFY(work_section)
             work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF",error=error)
             CALL dyn_coeff_set_initialize(dyn_coeff_set(isubsys)%dyn_coeff_set, root_section, para_env, globenv,&
                  work_section=work_section,error=error)
             CALL init_nhc_forces (nhc_coef=extended_type% nhc_coef(ibead),error=error)
          END IF
       END DO
       !
       ! Possibly restart Barostat
       !
       NULLIFY(work_section)
       work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT",error=error)
       CALL initialize_npt (simpar,para_env,globenv,extended_type%npt_info,&
            cell,work_section=work_section,error=error)
       !
       ! Always initialize velocities
       ! and possibly restart them
       !
       DO isubsys=1,SIZE(subsys)
         subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",&
                                                      i_rep_section=isubsys,&
                                                      error=error)
         CALL cp_subsys_get(subsys(isubsys)%subsys,&
                            atomic_kinds=atomic_kinds,&
                            molecule_kinds_new=molecule_kinds_new,&
                            particles=particles,&
                            shell_particles=shell_particles,&
                            core_particles=core_particles,&
                            error=error)

         CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
              shell_present=shell_present, shell_adiabatic=shell_adiabatic)
         NULLIFY (work_section,work2_section, work3_section)
         work_section  => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
         print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)

         IF(shell_present .AND. shell_adiabatic) THEN
            ! shell-model with constraint NYI
           CPPostcondition(.NOT.simpar%constraint,cp_failure_level,routineP,error,failure)
           work2_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
           work3_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
           CALL initialize_velocities(simpar,particles%els,&
                                      molecule_kinds=molecule_kinds_new,&
                                      para_env=para_env,&
                                      globenv=globenv,&
                                      print_section=print_section,&
                                      work_section=work_section,&
                                      work2_section=work2_section,&
                                      work3_section=work3_section,&
                                      shell_part=shell_particles%els,&
                                      core_part=core_particles%els,&
                                      error=error)
         ELSE
           CALL initialize_velocities(simpar,particles%els,&
                                      molecule_kinds=molecule_kinds_new,&
                                      para_env=para_env,&
                                      globenv=globenv,&
                                      print_section=print_section,&
                                      work_section=work_section,&
                                      error=error)
         END IF
       END DO

       ! Set V0 if needed
       IF ( ASSOCIATED ( cell_ref ) ) THEN
          IF ( simpar % v0 == 0._dp ) simpar % v0 = cell_ref % deth
       ENDIF

       CALL init_nhc_forces ( nhc_part=extended_type%nhc_part(ibead),error=error)
       CALL init_nhc_forces ( nhc_baro=extended_type%nhc_baro(ibead),error=error)

    END DO

    CALL md_env_create ( md_env, para_env, error = error )
    ! Setup Free Energy Calculation (if required)
    CALL fe_env_create ( fe_env, motion_section, error)

    CALL set_md_env(md_env=md_env,&
                    simpar=simpar,&
                    fe_env=fe_env,&
                    nhc_part=extended_type%nhc_part,&
                    nhc_baro=extended_type%nhc_baro,&
                    npt=extended_type%npt_info,&
                    cell=cell,&
                    force_env_p=force_env_p,error=error)

    IF (ASSOCIATED(dyn_coeff_set)) THEN
       CALL set_md_env( md_env, nhc_coef=extended_type%nhc_coef,error=error)
    END IF

    ! loop over number of beads
    IF (simpar%ensemble == langevin_ensemble) THEN
      CALL create_wiener_process(md_env,error)
    END IF

    time_iter_start=m_walltime()

    CALL get_md_env ( md_env, force_env_p=force_env_p, itimes = itimes, &
                      t = t, used_time = used_time )
    apply_cns0 = .FALSE.
    IF (simpar%constraint) &
    CALL section_vals_val_get(motion_section,"CONSTRAINT%CONSTRAIN_INIT",l_val=apply_cns0,error=error)
    DO ibead = 1, nbeads
       IF (apply_cns0) THEN
          CALL force_env_calc_energy_force ( force_env_p(ibead)%force_env,&
               calc_force=.TRUE.,error=error)          
          CALL force_env_shake(force_env_p(ibead)%force_env,shake_tol=simpar%shake_tol,&
               log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
               dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
          CALL force_env_rattle(force_env_p(ibead)%force_env,shake_tol=simpar%shake_tol,&
               log_unit=simpar%info_constraint,lagrange_mult=simpar%lagrange_multipliers,&
               dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)
          IF (simpar%do_respa)THEN
             CALL force_env_calc_energy_force ( &
                  force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
                  calc_force=.TRUE.,error=error)
             CALL force_env_shake(force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
                  shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                  lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,compold=.TRUE.,error=error)
             CALL force_env_rattle(force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
                  shake_tol=simpar%shake_tol,log_unit=simpar%info_constraint,&
                  lagrange_mult=simpar%lagrange_multipliers,dump_lm=simpar%dump_lm,reset=.TRUE.,error=error)          
          END IF
          ! Reinitialize velocities rescaling properly after rattle
          CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
               force_env_section=force_env_section, error=error )
          DO isubsys=1,SIZE(subsys)
             subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",&
                                                          i_rep_section=isubsys,&
                                                          error=error)
             CALL update_subsys(subsys_section,force_env_p(ibead)%force_env,.FALSE.,error)
             CALL cp_subsys_get(subsys(isubsys)%subsys,&
                                molecule_kinds_new=molecule_kinds_new,&
                                particles=particles,&
                                error=error)
             NULLIFY (work_section)
             work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
             print_section => section_vals_get_subs_vals(motion_section,"MD%PRINT",error=error)
             CALL initialize_velocities(simpar,particles%els,&
                                        molecule_kinds=molecule_kinds_new,&
                                        para_env=para_env,&
                                        globenv=globenv,&
                                        print_section=print_section,&
                                        work_section=work_section,&
                                        force_rescaling=.TRUE.,&
                                        error=error)
          END DO
       END IF

       IF(save_mem) THEN

          NULLIFY (work_section)
          work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
          CALL section_vals_remove_values(work_section, error)

          IF(shell_present .AND. shell_adiabatic) THEN
             NULLIFY (work_section)
             work_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
             CALL section_vals_remove_values(work_section, error)
             NULLIFY (work_section)
             work_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
            CALL section_vals_remove_values(work_section, error)
          END IF

       END IF

       CALL force_env_calc_energy_force ( force_env_p(ibead)%force_env,&
            calc_force=.TRUE.,error=error)
 
       IF (simpar%do_respa)THEN
          CALL force_env_calc_energy_force ( &
            force_env_p(ibead)%force_env%sub_force_env(1)%force_env,&
            calc_force=.TRUE.,error=error)
       END IF
       CALL force_env_get( force_env_p(ibead)%force_env, subsys=subsys, &
                           error=error )

       CALL cp_subsys_get(subsys(1)%subsys,&
              atomic_kinds=atomic_kinds,&
              local_particles=local_particles,&
              particles=particles,&
              error=error)

       CALL virial_evaluate(atomic_kinds%els, particles%els, &
              local_particles, force_env_p(ibead)%force_env%virial,&
              force_env_p(ibead)%force_env%para_env%group,&
              error=error)

    END DO

    CALL section_vals_val_get(motion_section,"MD%STEP_START_VAL",i_val=itimes,error=error)
    CALL section_vals_val_get(motion_section,"MD%TIME_START_VAL",r_val=t,error=error)
    CALL md_energy ( md_env, error ) ! inits the print env at itimes == 0

    CALL cp_print_key_finished_output(simpar%info_constraint, logger,&
         root_section,"MOTION%CONSTRAINT%CONSTRAINT_INFO",error=error)
    CALL cp_print_key_finished_output(simpar%lagrange_multipliers, logger,&
         root_section,"MOTION%CONSTRAINT%LAGRANGE_MULTIPLIERS",error=error)

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
      DO ibead=1,nbeads
        CALL temperature_control(simpar,force_env_p(ibead)%force_env,logger,error)
        CALL comvel_control(simpar,force_env_p(ibead)%force_env,logger,error)
      END DO
      time_iter_stop=m_walltime()
      used_time = time_iter_stop - time_iter_start
      time_iter_start=time_iter_stop

      CALL md_output(md_env,force_env_section,motion_section,root_section,error)
    END DO

    CALL cp_rm_iter_level(logger%iter_info,n_rlevel_att=iter_level,&
         error=error)
    CALL md_env_release ( md_env, error )

    ! just deallocate this stuff ??
    ! seem to need proper init / destroy subroutines ... to be fixed
    CALL lnhc_dealloc(extended_type%nhc_part)
    CALL lnhc_dealloc(extended_type%nhc_baro)
    CALL lnhc_dealloc(extended_type%nhc_coef)
    CALL npt_dealloc(extended_type%npt_info)

    DEALLOCATE(simpar,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_mol_dyn2

!!****f* md_run/md_output *
!!
!!   NAME
!!     md_output
!!
!!   FUNCTION
!!     collects the part of the MD that, basically, does the output
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE md_output(md_env,force_env_section,motion_section,root_section,error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                motion_section, root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'md_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)


    CALL md_energy ( md_env, error ) ! also writes trajectories

    CALL write_restart(md_env=md_env,root_section=root_section,&
         force_env_section=force_env_section,error=error)

    CALL timestop(handle)

END SUBROUTINE md_output

! *****************************************************************************

  SUBROUTINE temperature_control(simpar,force_env,logger, error)

    ! Perform all temperature manipulations during a QS MD run.

    ! - Creation (15.09.2003,MK)
    ! - adapted to force_env (05.10.2003,fawzi)

    ! *************************************************************************

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'temperature_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, bo(2), handle, iatom, &
                                                ikind, isubsys, iw, mepos, &
                                                natom, nfree_shell, num_pe, &
                                                shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, is_shell, &
                                                shell_adiabatic
    REAL(KIND=dp) :: e_kinetic, ekin_c, ekin_com, ekin_s, ekin_shell, &
      fac_massc, fac_masss, mass, scale, temperature, total_mass, v2
    REAL(KIND=dp), DIMENSION(3)              :: v, vc, vcom, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(atomic_kind,atomic_kinds,particles,atom_list,subsys,para_env)
    NULLIFY( core_particles, shell_particles, shell)

    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    IF (simpar%temp_tol > 0.0_dp ) THEN

      CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)
      iw = cp_print_key_unit_nr(logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
           extension=".mdLog",error=error)


      e_kinetic = 0.0_dp
      total_mass = 0.0_dp
      vcom(:) = 0.0_dp
      DO isubsys=1,SIZE(subsys)

        CALL cp_subsys_get(subsys(isubsys)%subsys,&
                           atomic_kinds=atomic_kinds,&
                           particles=particles,&
                           error=error)

        ! Calculate the current temperature

        DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass,&
                               natom=natom)
          v2 = 0.0_dp
          v(:) = 0.0_dp
          DO iatom=1,natom
            atom = atom_list(iatom)
            v2 = v2 + particles%els(atom)%v(1)*particles%els(atom)%v(1) +&
                      particles%els(atom)%v(2)*particles%els(atom)%v(2) +&
                      particles%els(atom)%v(3)*particles%els(atom)%v(3)
            v(:) = v(:) + particles%els(atom)%v(:)
          END DO
          e_kinetic = e_kinetic + 0.5_dp*mass*v2
          vcom(:) = vcom(:) + mass*v(:)
          total_mass = total_mass + REAL(natom,KIND=dp)*mass
        END DO

      END DO  ! isubsys
      
      temperature = 0.0_dp
      IF (simpar%nfree /=0) temperature = 2.0_dp*e_kinetic/REAL(simpar%nfree,KIND=dp)

      IF (iw>0) THEN
        WRITE (UNIT=iw,FMT="(/,T2,A,(T58,A3,F20.10))")&
          "Centre of mass motion (COM):","x =",vcom(1)/total_mass,&
                                         "y =",vcom(2)/total_mass,&
                                         "z =",vcom(3)/total_mass
      END IF

      ! Scale temperature to reference temperture, if requested

      IF (ABS(simpar%temp_ext - temperature) > simpar%temp_tol) THEN

        scale = SQRT(simpar%temp_ext/temperature)

        DO isubsys=1,SIZE(subsys)

          CALL cp_subsys_get(subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             particles=particles,&
                             shell_particles=shell_particles,&
                             core_particles=core_particles,&
                             error=error)

          CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
               shell_adiabatic=shell_adiabatic)

          IF(shell_adiabatic) THEN
            nfree_shell = 3*SIZE(shell_particles%els)
            ekin_shell = 0.0_dp
            DO ikind=1,atomic_kinds%n_els
              atomic_kind => atomic_kinds%els(ikind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   atom_list=atom_list,&
                                   mass=mass,&
                                   natom=natom,&
                                   shell_active=is_shell,shell=shell)
              IF(is_shell) THEN
                fac_masss = shell%mass_shell/mass 
                fac_massc = shell%mass_core/mass 
                DO iatom=1,natom
                  atom = atom_list(iatom)
                  particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
                  shell_index = particles%els(atom)%shell_index
                  vs(1:3) = shell_particles%els(shell_index)%v(1:3)
                  vc(1:3) = core_particles%els(shell_index)%v(1:3)
                  shell_particles%els(shell_index)%v(1) = &
                        particles%els(atom)%v(1) + fac_massc*(vs(1)-vc(1))
                  shell_particles%els(shell_index)%v(2) = &
                        particles%els(atom)%v(2) + fac_massc*(vs(2)-vc(2))
                  shell_particles%els(shell_index)%v(3) = &
                        particles%els(atom)%v(3) + fac_massc*(vs(3)-vc(3))
                  core_particles%els(shell_index)%v(1) = &
                        particles%els(atom)%v(1) + fac_masss*(vc(1)-vs(1))
                  core_particles%els(shell_index)%v(2) = &
                        particles%els(atom)%v(2) + fac_masss*(vc(2)-vs(2))
                  core_particles%els(shell_index)%v(3) = &
                        particles%els(atom)%v(3) + fac_masss*(vc(3)-vs(3))
                END DO  ! iatom
              ELSE
                DO iatom=1,natom
                  atom = atom_list(iatom)
                  particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
                END DO  ! iatom
              END IF
            END DO  ! ikind
          ELSE
            DO ikind=1,atomic_kinds%n_els
              atomic_kind => atomic_kinds%els(ikind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   atom_list=atom_list,&
                                   mass=mass,&
                                   natom=natom)
              v2 = 0.0_dp
              DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
              END DO
            END DO
          END IF

        END DO

        IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
            "Temperature scaled to requested temperature:",&
            temperature*kelvin," K  ->",simpar%temp_ext*kelvin," K"
        END IF

      END IF

      CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

!   control the internal core-shell motion
    IF(simpar%temp_sh_tol > 0.0_dp) THEN


      ekin_shell = 0.0_dp
      ekin_c = 0.0_dp
      ekin_s = 0.0_dp
      nfree_shell = 0

      CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)
      iw = cp_print_key_unit_nr(logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
           extension=".mdLog",error=error)

      num_pe = para_env%num_pe
      mepos  = para_env%mepos

      DO isubsys=1,SIZE(subsys)

        CALL cp_subsys_get(subsys(isubsys)%subsys,&
                           atomic_kinds=atomic_kinds,&
                           particles=particles,&
                           shell_particles=shell_particles,&
                           core_particles=core_particles,&
                           error=error)

        ! Calculate the current temperature

        CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
             shell_adiabatic=shell_adiabatic)

        IF(shell_adiabatic) THEN
          nfree_shell = nfree_shell + 3*SIZE(shell_particles%els)

          DO ikind=1,atomic_kinds%n_els
            atomic_kind => atomic_kinds%els(ikind)
            CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                 atom_list=atom_list,&
                                 mass=mass,&
                                 natom=natom,&
                                 shell_active=is_shell,shell=shell)

            bo = get_limit( natom, num_pe, mepos )
            IF(is_shell) THEN
               fac_masss = shell%mass_shell/mass 
               fac_massc = shell%mass_core/mass 
               DO iatom=bo(1),bo(2)
                  atom = atom_list(iatom)
                  shell_index = particles%els(atom)%shell_index
                  vs(1:3) = shell_particles%els(shell_index)%v(1:3)
                  vc(1:3) =  core_particles%els(shell_index)%v(1:3)
                  ekin_c = 0.5_dp * shell%mass_core  * (vc(1)*vc(1) + vc(2)*vc(2) + vc(3)*vc(3))
                  ekin_s = 0.5_dp * shell%mass_shell * (vs(1)*vs(1) + vs(2)*vs(2) + vs(3)*vs(3))

                  v(:) =  particles%els(atom)%v(:)
                  ekin_com = 0.5_dp*mass * (v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
                  ekin_shell = ekin_shell + ekin_c + ekin_s - ekin_com
               END DO  ! iatom
            END IF  ! is_shell
          END DO  ! ikind
        END IF  ! shell_adiabatic
      END DO ! isub
      CALL mp_sum(ekin_shell,para_env%group)

      temperature = 0.0_dp
      IF(nfree_shell > 0) THEN
        temperature = 2.0_dp*ekin_shell/REAL(nfree_shell,KIND=dp)

        IF (ABS(simpar%temp_sh_ext - temperature) > simpar%temp_sh_tol) THEN

           scale = 0.0_dp
           IF (temperature>EPSILON(0.0_dp)) scale = SQRT(simpar%temp_sh_ext/temperature)

           total_mass = 0.0_dp
           vcom(:)    = 0.0_dp
           e_kinetic  = 0.0_dp
           DO isubsys = 1,SIZE(subsys)

             CALL cp_subsys_get(subsys(isubsys)%subsys,&
                           atomic_kinds=atomic_kinds,&
                           particles=particles,&
                           shell_particles=shell_particles,&
                           core_particles=core_particles,&
                           error=error)

             total_mass = 0.0_dp
             DO ikind=1,atomic_kinds%n_els
               atomic_kind => atomic_kinds%els(ikind)
               CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                    atom_list=atom_list,&
                                    mass=mass,&
                                    natom=natom,&
                                    shell_active=is_shell,shell=shell)
               IF(is_shell) THEN
                  fac_massc = 1.0_dp/mass 
                  DO iatom= 1, natom
                     atom = atom_list(iatom)
                     shell_index = particles%els(atom)%shell_index
                     vs(1:3) = shell_particles%els(shell_index)%v(1:3)
                     vc(1:3) = core_particles%els(shell_index)%v(1:3)
                     v(:) = particles%els(atom)%v(:)

                     shell_particles%els(shell_index)%v(1:3) = v(1:3) -&
                              shell%mass_core*scale*fac_massc*(vc(1:3)-vs(1:3))  
                     core_particles%els(shell_index)%v(1:3) =  v(1:3) + &
                              shell%mass_shell*scale*fac_massc*(vc(1:3)-vs(1:3))

                     vcom = vcom + v * mass
                     total_mass = total_mass+mass
                     e_kinetic = e_kinetic + &
                                 0.5_dp*mass*(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
                  END DO  ! iatom
               ELSE
                  DO iatom= 1,natom
                     atom = atom_list(iatom)
                     v(:) = particles%els(atom)%v(:)
                     vcom = vcom + v * mass
                     total_mass = total_mass+mass
                  END DO
               END IF  ! is_shell
             END DO  ! ikind
           END DO ! isub

           IF (iw>0) THEN
             WRITE (UNIT=6,FMT="(/,T2,A,F10.2,A,F10.2,A)")&
                  "Temperature  shell internal motion scaled to requested temperature:",&
                  temperature*kelvin," K  ->",simpar%temp_sh_ext*kelvin," K"
             WRITE (UNIT=6,FMT="(/,T2,A,(T58,A3,F20.10))")&
                  "Centre of mass motion (COM):","x =",vcom(1)/total_mass,&
                                                 "y =",vcom(2)/total_mass,&
                                                 "z =",vcom(3)/total_mass
           END IF

        ENDIF
      END IF  ! nfree_shell
      CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF  ! simpar%temp_sh_tol > 0

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE temperature_control

! *****************************************************************************

  SUBROUTINE comvel_control(simpar,force_env,logger, error)

    ! Set to 0 the velocity of the COM along MD runs, if required.

    ! - Creation (29.04.2007,MI)

    ! *************************************************************************

    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'comvel_control', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, bo(2), handle, iatom, &
                                                ikind, isubsys, iw, mepos, &
                                                natom, num_pe, shell_index
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, is_shell, &
                                                shell_adiabatic
    REAL(KIND=dp)                            :: e_kinetic, fac_massc, &
                                                fac_masss, mass, scale, &
                                                temp_new, temperature, &
                                                total_mass, v2, vel_com
    REAL(KIND=dp), DIMENSION(3)              :: v, vc, vcom, vcom_new, vs
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(shell_kind_type), POINTER           :: shell

! -------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    NULLIFY(atomic_kind,atomic_kinds,particles,atom_list,subsys,para_env)
    NULLIFY( core_particles, shell_particles, shell)

    CPPrecondition(ASSOCIATED(simpar),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(force_env),cp_failure_level,routineP,error,failure)

    IF (simpar%comvel_tol > 0.0_dp ) THEN

      CALL force_env_get(force_env,subsys=subsys,para_env=para_env,error=error)
      iw = cp_print_key_unit_nr(logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO",&
           extension=".mdLog",error=error)


      num_pe = para_env%num_pe
      mepos  = para_env%mepos

      e_kinetic = 0.0_dp
      total_mass = 0.0_dp
      vcom(:) = 0.0_dp
      DO isubsys=1,SIZE(subsys)

        CALL cp_subsys_get(subsys(isubsys)%subsys,&
                           atomic_kinds=atomic_kinds,&
                           particles=particles,&
                           error=error)

        ! Calculate the current com velocity and the temperature

        DO ikind=1,atomic_kinds%n_els
          atomic_kind => atomic_kinds%els(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass,&
                               natom=natom)
          v2 = 0.0_dp
          v(:) = 0.0_dp
          bo = get_limit( natom, num_pe, mepos )
          DO iatom=bo(1),bo(2)
            atom = atom_list(iatom)
            v2 = v2 + particles%els(atom)%v(1)*particles%els(atom)%v(1) +&
                      particles%els(atom)%v(2)*particles%els(atom)%v(2) +&
                      particles%els(atom)%v(3)*particles%els(atom)%v(3)
            v(:) = v(:) + particles%els(atom)%v(:)
          END DO
          e_kinetic = e_kinetic + 0.5_dp*mass*v2
          vcom(:) = vcom(:) + mass*v(:)
          total_mass = total_mass + REAL(natom,KIND=dp)*mass
        END DO
      END DO  ! isubsys
      CALL mp_sum(e_kinetic,para_env%group)
      CALL mp_sum(vcom,para_env%group)
      CALL mp_sum(total_mass,para_env%group)

      temperature = 0.0_dp
      IF (simpar%nfree /= 0) temperature = 2.0_dp*e_kinetic/REAL(simpar%nfree,KIND=dp)
      vcom = vcom / total_mass
      vel_com = SQRT(vcom(1)*vcom(1)+vcom(2)*vcom(2)+vcom(3)*vcom(3))

      ! Subtract the velocity of the COM, if requested

      IF (vel_com > simpar%comvel_tol) THEN

        e_kinetic = 0.0_dp
        v(:) = 0.0_dp
        DO isubsys=1,SIZE(subsys)

          CALL cp_subsys_get(subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             particles=particles,&
                             shell_particles=shell_particles,&
                             core_particles=core_particles,&
                             error=error)

          CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
               shell_adiabatic=shell_adiabatic)

          IF(shell_adiabatic) THEN
            DO ikind=1,atomic_kinds%n_els
              atomic_kind => atomic_kinds%els(ikind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   atom_list=atom_list,&
                                   mass=mass,&
                                   natom=natom,&
                                   shell_active=is_shell,shell=shell)
              IF(is_shell) THEN
                v2 = 0.0_dp
                DO iatom=1,natom
                  atom = atom_list(iatom)
                  shell_index = particles%els(atom)%shell_index
                  particles%els(atom)%v(:) = particles%els(atom)%v(:) - vcom(:)
                  v2 = v2 + particles%els(atom)%v(1)*particles%els(atom)%v(1) +&
                            particles%els(atom)%v(2)*particles%els(atom)%v(2) +&
                            particles%els(atom)%v(3)*particles%els(atom)%v(3)
                  shell_particles%els(shell_index)%v(:) = shell_particles%els(shell_index)%v(:) -vcom(:)
                  core_particles%els(shell_index)%v(:) = core_particles%els(shell_index)%v(:) -vcom(:)
                END DO  ! iatom
                e_kinetic = e_kinetic + 0.5_dp*mass*v2
              ELSE
                v2 = 0.0_dp
                DO iatom=1,natom
                  atom = atom_list(iatom)
                  particles%els(atom)%v(:) = particles%els(atom)%v(:) - vcom(:)
                  v2 = v2 + particles%els(atom)%v(1)*particles%els(atom)%v(1) +&
                            particles%els(atom)%v(2)*particles%els(atom)%v(2) +&
                            particles%els(atom)%v(3)*particles%els(atom)%v(3)
                END DO  ! iatom
                e_kinetic = e_kinetic + 0.5_dp*mass*v2
              END IF
            END DO  ! ikind
          ELSE
            DO ikind=1,atomic_kinds%n_els
              atomic_kind => atomic_kinds%els(ikind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   atom_list=atom_list,&
                                   mass=mass,&
                                   natom=natom)
              v2 = 0.0_dp
              DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(:) = particles%els(atom)%v(:)-vcom(:)
                v2 = v2 + particles%els(atom)%v(1)*particles%els(atom)%v(1) +&
                          particles%els(atom)%v(2)*particles%els(atom)%v(2) +&
                          particles%els(atom)%v(3)*particles%els(atom)%v(3)
              END DO
              e_kinetic = e_kinetic + 0.5_dp*mass*v2
            END DO  ! ikind
          END IF  ! shell_adiabatic

        END DO  ! isub
        
        temp_new = 0.0_dp
        IF (simpar%nfree /=0) temp_new = 2.0_dp*e_kinetic/REAL(simpar%nfree,KIND=dp)
        scale = SQRT( temperature / temp_new)
        vcom_new = 0.0_dp
        DO isubsys=1,SIZE(subsys)

          CALL cp_subsys_get(subsys(isubsys)%subsys,&
                             atomic_kinds=atomic_kinds,&
                             particles=particles,&
                             shell_particles=shell_particles,&
                             core_particles=core_particles,&
                             error=error)

          CALL get_atomic_kind_set(atomic_kind_set=atomic_kinds%els,&
               shell_adiabatic=shell_adiabatic)

          IF(shell_adiabatic) THEN
            DO ikind=1,atomic_kinds%n_els
              atomic_kind => atomic_kinds%els(ikind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   atom_list=atom_list,&
                                   mass=mass,&
                                   natom=natom,&
                                   shell_active=is_shell,shell=shell)
              IF(is_shell) THEN
                fac_masss = shell%mass_shell/mass
                fac_massc = shell%mass_core/mass
                v(:) = 0.0_dp
                DO iatom=1,natom
                  atom = atom_list(iatom)
                  particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
                  v(:) = v(:) + particles%els(atom)%v(:) 
                  shell_index = particles%els(atom)%shell_index
                  vs(1:3) = shell_particles%els(shell_index)%v(1:3)
                  vc(1:3) = core_particles%els(shell_index)%v(1:3)
                  shell_particles%els(shell_index)%v(1) = &
                        particles%els(atom)%v(1) + fac_massc*(vs(1)-vc(1))
                  shell_particles%els(shell_index)%v(2) = &
                        particles%els(atom)%v(2) + fac_massc*(vs(2)-vc(2))
                  shell_particles%els(shell_index)%v(3) = &
                        particles%els(atom)%v(3) + fac_massc*(vs(3)-vc(3))
                  core_particles%els(shell_index)%v(1) = &
                        particles%els(atom)%v(1) + fac_masss*(vc(1)-vs(1))
                  core_particles%els(shell_index)%v(2) = &
                        particles%els(atom)%v(2) + fac_masss*(vc(2)-vs(2))
                  core_particles%els(shell_index)%v(3) = &
                        particles%els(atom)%v(3) + fac_masss*(vc(3)-vs(3))
                END DO  ! iatom
                vcom_new(:) = vcom_new(:) + mass*v(:)
              ELSE
                v(:) = 0.0_dp
                DO iatom=1,natom
                  atom = atom_list(iatom)
                  particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
                  v(:) = v(:) + particles%els(atom)%v(:) 
                END DO  ! iatom
                vcom_new(:) = vcom_new(:) + mass*v(:)
              END IF
            END DO  ! ikind
          ELSE
            DO ikind=1,atomic_kinds%n_els
              atomic_kind => atomic_kinds%els(ikind)
              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   atom_list=atom_list,&
                                   mass=mass,&
                                   natom=natom)
              v(:) = 0.0_dp
              DO iatom=1,natom
                atom = atom_list(iatom)
                particles%els(atom)%v(:) = scale*particles%els(atom)%v(:)
                v(:) = v(:) + particles%els(atom)%v(:) 
              END DO
              vcom_new(:) = vcom_new(:) + mass*v(:)
            END DO
          END IF

        END DO  ! isub
        IF (iw>0) THEN
          WRITE (UNIT=iw,FMT="(/,T2,A)") &
            " #########  The velocity of the COM has been subtracted:"
          WRITE (UNIT=iw,FMT="(T2,10X,A,3F16.10,A)") &
            "VCOM = ",vcom(1:3)," a.u.  ->"
          WRITE (UNIT=iw,FMT="(T2,16X,A,3F16.10,A)") &
            "  ->", vcom_new(1:3)," a.u"
        END IF

      END IF  !vel_com > simpar%comvel_tol

      CALL cp_print_key_finished_output(iw,logger,force_env%root_section,&
           "MOTION%MD%PRINT%PROGRAM_RUN_INFO", error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE comvel_control

  ! ***************************************************************************
END MODULE md_run
