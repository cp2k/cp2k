!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of routines to dump the restart file of CP2K
!> \par History
!>      01.2006 [created] Teodoro Laino
! *****************************************************************************
MODULE input_cp2k_restarts
  USE al_system_types,                 ONLY: al_system_type
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE averages_types,                  ONLY: average_quantities_type
  USE cell_methods,                    ONLY: cell_type,&
                                             real_to_scaled
  USE cp2k_info,                       ONLY: write_restart_header
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_linked_list_val,              ONLY: cp_sll_val_create,&
                                             cp_sll_val_get_length,&
                                             cp_sll_val_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE csvr_system_types,               ONLY: csvr_system_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             map_info_type,&
                                             npt_info_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             multiple_fe_list,&
                                             use_qmmm
  USE gle_system_types,                ONLY: gle_type
  USE helium_types,                    ONLY: helium_solvent_type
  USE input_constants,                 ONLY: &
       do_band_collective, do_coord_cp2k, do_thermo_al, do_thermo_csvr, &
       do_thermo_gle, do_thermo_no_communication, do_thermo_nose, ehrenfest, &
       mol_dyn_run, mon_car_run, pint_run, use_rt_restart
  USE input_cp2k_restarts_util,        ONLY: section_velocity_val_set
  USE input_section_types,             ONLY: &
       section_get_keyword_index, section_type, section_vals_add_values, &
       section_vals_get, section_vals_get_subs_vals, &
       section_vals_get_subs_vals3, section_vals_remove_values, &
       section_vals_type, section_vals_val_get, section_vals_val_set, &
       section_vals_val_unset, section_vals_write
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_gather,&
                                             mp_sum
  USE metadynamics_types,              ONLY: meta_env_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE multipole_types,                 ONLY: multipole_type
  USE neb_types,                       ONLY: neb_var_type
  USE parallel_rng_types,              ONLY: dump_rng_stream,&
                                             get_rng_stream,&
                                             rng_record_length
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: get_particle_pos_or_vel,&
                                             particle_type
  USE physcon,                         ONLY: angstrom
  USE pint_transformations,            ONLY: pint_u2x
  USE pint_types,                      ONLY: pint_env_type,&
                                             thermostat_gle,&
                                             thermostat_nose
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: ascii_to_string,&
                                             string_to_ascii
  USE termination,                     ONLY: stop_memory,&
                                             stop_program,&
                                             stop_write
  USE thermostat_types,                ONLY: thermostat_type
  USE thermostat_utils,                ONLY: communication_thermo_low2,&
                                             get_kin_energies
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_restarts'

  PUBLIC :: write_restart,&
            update_subsys,&
            update_input

CONTAINS

! *****************************************************************************
!> \brief checks if a restart needs to be written and does so, updating all necessary fields
!>      in the input file. This is a relatively simple wrapper routine.
!> \param md_env ...
!> \param force_env ...
!> \param root_section ...
!> \param coords ...
!> \param vels ...
!> \param pint_env ...
!> \param helium_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE write_restart(md_env,force_env,root_section,&
       coords, vels, pint_env, helium_env, error)
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(neb_var_type), OPTIONAL, POINTER    :: coords, vels
    TYPE(pint_env_type), OPTIONAL, POINTER   :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_restart', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=30), DIMENSION(2), PARAMETER :: &
      keys = (/"PRINT%RESTART_HISTORY","PRINT%RESTART        "/)

    INTEGER                                  :: handle, ikey, ires, &
                                                nforce_eval
    LOGICAL                                  :: failure, save_mem, &
                                                write_binary_restart_file
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: global_section, &
                                                motion_section, sections

    CALL timeset(routineN,handle)

    failure = .FALSE.
    logger => cp_error_get_logger(error)
    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)

    NULLIFY(global_section)
    global_section => section_vals_get_subs_vals(root_section,"GLOBAL",error=error)
    CALL section_vals_val_get(global_section,"SAVE_MEM",l_val=save_mem,error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              motion_section,keys(1),error=error),cp_p_file).OR. &
        BTEST(cp_print_key_should_output(logger%iter_info,&
              motion_section,keys(2),error=error),cp_p_file)) THEN

       sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
       CALL section_vals_get(sections, n_repetition=nforce_eval, error=error)
       CALL section_vals_val_get(motion_section,"PRINT%RESTART%SPLIT_RESTART_FILE",&
                                 l_val=write_binary_restart_file,error=error)

       IF (write_binary_restart_file) THEN
          CALL update_subsys_release(md_env,force_env,root_section,error)
          CALL update_motion_release(motion_section,error)
          DO ikey=1,SIZE(keys)
             IF (BTEST(cp_print_key_should_output(logger%iter_info,&
                       motion_section,keys(ikey),error=error),cp_p_file)) THEN
                ires = cp_print_key_unit_nr(logger,motion_section,TRIM(keys(ikey)),&
                                            extension=".restart.bin",&
                                            file_action="READWRITE",&
                                            file_form="UNFORMATTED",&
                                            file_position="REWIND",&
                                            file_status="UNKNOWN",&
                                            do_backup=(ikey==2),&
                                            error=error)
                CALL write_binary_restart(ires,root_section,md_env,force_env,error=error)
                CALL cp_print_key_finished_output(ires,logger,motion_section,&
                                                  TRIM(keys(ikey)),error=error)
             END IF
          END DO
       END IF

       CALL update_input(md_env,force_env,root_section,coords,vels,pint_env,helium_env,&
                         save_mem=save_mem,&
                         write_binary_restart_file=write_binary_restart_file,&
                         error=error)

       DO ikey=1,SIZE(keys)
          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
                    motion_section,keys(ikey),error=error),cp_p_file)) THEN
             ires = cp_print_key_unit_nr(logger,motion_section,TRIM(keys(ikey)),&
                                         extension=".restart",&
                                         file_position="REWIND",&
                                         do_backup=(ikey==2),&
                                         error=error)
             IF (ires > 0) THEN
                CALL write_restart_header(ires)
                CALL section_vals_write(root_section,unit_nr=ires,hide_root=.TRUE.,error=error)
             END IF
             CALL cp_print_key_finished_output(ires,logger,motion_section,TRIM(keys(ikey)),error=error)
          END IF
       END DO

       IF (save_mem) THEN
          CALL update_subsys_release(md_env,force_env,root_section,error)
          CALL update_motion_release(motion_section,error)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE write_restart

! *****************************************************************************
!> \brief deallocate some sub_sections of the section subsys to save some memory
!> \param md_env ...
!> \param force_env ...
!> \param root_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE update_subsys_release(md_env,force_env,root_section,error)

    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_subsys_release', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: unit_str
    INTEGER                                  :: handle, iforce_eval, myid, &
                                                nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: i_force_eval
    LOGICAL                                  :: explicit, failure, scale, &
                                                skip_vel_section
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: my_force_b, my_force_env
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(section_vals_type), POINTER         :: force_env_sections, &
                                                subsys_section, work_section

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(core_particles,my_force_env, my_force_b, particles,&
         shell_particles,subsys,work_section)

    IF (PRESENT(md_env)) THEN
       CALL get_md_env ( md_env=md_env,force_env=my_force_env,error=error)
    ELSEIF (PRESENT(force_env)) THEN
       my_force_env => force_env
    END IF

    IF (ASSOCIATED(my_force_env))  THEN
       NULLIFY(subsys_section)
       CALL section_vals_val_get(root_section,"GLOBAL%RUN_TYPE",i_val=myid,error=error)
       skip_vel_section =(&
            (myid /= mol_dyn_run).AND.&
            (myid /= mon_car_run).AND.&
            (myid /= pint_run))

       force_env_sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
       CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, nforce_eval, error)

       DO iforce_eval=1,nforce_eval
          subsys_section => section_vals_get_subs_vals3(force_env_sections,"SUBSYS",&
               i_rep_section=i_force_eval(iforce_eval),error=error)
          CALL section_vals_get(subsys_section, explicit=explicit, error=error)
          IF (.NOT.explicit) CYCLE ! Nothing to update...

          my_force_b => my_force_env
          IF (iforce_eval>1) my_force_b => my_force_env%sub_force_env(iforce_eval-1)%force_env

          CALL force_env_get(my_force_b, subsys=subsys, error=error)

          CALL cp_subsys_get(subsys, particles=particles, shell_particles=shell_particles,&
               core_particles=core_particles, error=error)

          work_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
          CALL section_vals_get(work_section,explicit=explicit,error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,error=error)
             CALL section_vals_val_get(work_section,"SCALED",l_val=scale,error=error)
          END IF
          CALL section_vals_remove_values(work_section, error)
          IF (explicit) THEN
             CALL section_vals_val_set(work_section,"UNIT",c_val=unit_str,error=error)
             CALL section_vals_val_set(work_section,"SCALED",l_val=scale,error=error)
          END IF

          work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
          IF (.NOT.skip_vel_section) THEN
             CALL section_vals_remove_values(work_section, error)
          END IF

          IF (ASSOCIATED(shell_particles)) THEN
             work_section => section_vals_get_subs_vals(subsys_section,"SHELL_COORD",error=error)
             CALL section_vals_get(work_section,explicit=explicit,error=error)
             IF (explicit) THEN
                CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,error=error)
                CALL section_vals_val_get(work_section,"SCALED",l_val=scale,error=error)
             END IF
             CALL section_vals_remove_values(work_section, error=error)
             IF (explicit) THEN
                CALL section_vals_val_set(work_section,"UNIT",c_val=unit_str,error=error)
                CALL section_vals_val_set(work_section,"SCALED",l_val=scale,error=error)
             END IF

             work_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
             IF (.NOT.skip_vel_section) THEN
                CALL section_vals_remove_values(work_section, error=error)
             END IF
          END IF

          IF (ASSOCIATED(core_particles)) THEN
             work_section => section_vals_get_subs_vals(subsys_section,"CORE_COORD",error=error)
             CALL section_vals_get(work_section,explicit=explicit,error=error)
             IF (explicit) THEN
                CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,error=error)
                CALL section_vals_val_get(work_section,"SCALED",l_val=scale,error=error)
             END IF
             CALL section_vals_remove_values(work_section, error=error)
             IF (explicit) THEN
                CALL section_vals_val_set(work_section,"UNIT",c_val=unit_str,error=error)
                CALL section_vals_val_set(work_section,"SCALED",l_val=scale,error=error)
             END IF

             work_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
             IF (.NOT.skip_vel_section) THEN
                CALL section_vals_remove_values(work_section, error=error)
             END IF
          END IF

       END DO

       DEALLOCATE (i_force_eval)

    END IF

    CALL timestop(handle)

  END SUBROUTINE update_subsys_release

! *****************************************************************************
!> \brief deallocate the nose subsections (coord, vel, force, mass) in the md section
!> \param motion_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE update_motion_release(motion_section,error)

    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_motion_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(section_vals_type), POINTER         :: work_section

    CALL timeset(routineN,handle)

    NULLIFY (work_section)

    work_section => section_vals_get_subs_vals(motion_section,"MD%AVERAGES%RESTART_AVERAGES",error=error)
    CALL section_vals_remove_values(work_section,error)

    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section,error)

    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section,error)

    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section,error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section,error)

    CALL timestop(handle)

  END SUBROUTINE update_motion_release

! *****************************************************************************
!> \brief Updates the whole input file for the restart
!> \param md_env ...
!> \param force_env ...
!> \param root_section ...
!> \param coords ...
!> \param vels ...
!> \param pint_env ...
!> \param helium_env ...
!> \param save_mem ...
!> \param write_binary_restart_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_input(md_env,force_env,root_section,coords,vels,pint_env, &
                          helium_env,save_mem,write_binary_restart_file,error)

    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(neb_var_type), OPTIONAL, POINTER    :: coords, vels
    TYPE(pint_env_type), OPTIONAL, POINTER   :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_env
    LOGICAL, INTENT(IN), OPTIONAL            :: save_mem, &
                                                write_binary_restart_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, lcond, my_save_mem, &
                                                my_write_binary_restart_file
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_env_type), POINTER            :: my_force_env
    TYPE(section_vals_type), POINTER         :: motion_section

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(logger,motion_section,my_force_env)

    IF (PRESENT(save_mem)) THEN
       my_save_mem = save_mem
    ELSE
       my_save_mem = .FALSE.
    END IF

    IF (PRESENT(write_binary_restart_file)) THEN
       my_write_binary_restart_file = write_binary_restart_file
    ELSE
       my_write_binary_restart_file = .FALSE.
    END IF

    logger => cp_error_get_logger(error)

    ! Can handle md_env or force_env
    lcond = PRESENT(md_env).OR.PRESENT(force_env).OR.PRESENT(pint_env).OR.PRESENT(helium_env)
    IF (lcond) THEN
       IF (PRESENT(md_env)) THEN
          CALL get_md_env ( md_env=md_env, force_env=my_force_env, error=error)
       ELSE IF (PRESENT(force_env)) THEN
          my_force_env => force_env
       END IF
       ! The real restart setting...
       motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
       CALL update_motion(motion_section,&
                          md_env=md_env,&
                          force_env=my_force_env,&
                          logger=logger,&
                          coords=coords,&
                          vels=vels,&
                          pint_env=pint_env,&
                          helium_env=helium_env,&
                          save_mem=my_save_mem,&
                          write_binary_restart_file=my_write_binary_restart_file,&
                          error=error)
       ! Update one force_env_section per time..
       IF (ASSOCIATED(my_force_env)) THEN
          CALL update_force_eval(force_env=my_force_env,&
                                 logger=logger,&
                                 root_section=root_section,&
                                 write_binary_restart_file=my_write_binary_restart_file,&
                                 error=error)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE update_input

! *****************************************************************************
!> \brief Updates the motion section of the input file
!> \param motion_section ...
!> \param md_env ...
!> \param force_env ...
!> \param logger ...
!> \param coords ...
!> \param vels ...
!> \param pint_env ...
!> \param helium_env ...
!> \param save_mem ...
!> \param write_binary_restart_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_motion(motion_section,md_env,force_env,logger,&
                           coords,vels,pint_env,helium_env,save_mem,&
                           write_binary_restart_file,error)
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(neb_var_type), OPTIONAL, POINTER    :: coords, vels
    TYPE(pint_env_type), OPTIONAL, POINTER   :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_env
    LOGICAL, INTENT(IN), OPTIONAL            :: save_mem, &
                                                write_binary_restart_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_motion', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, handle, handle2, i, &
                                                irep, isec, j, nhc_len, stat, &
                                                tot_nhcneed
    INTEGER, DIMENSION(:), POINTER           :: walkers_status
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, my_save_mem, &
                                                my_write_binary_restart_file
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer, eta, fnhc, mnhc, &
                                                veta, wrk
    REAL(KIND=dp), POINTER                   :: constant, t
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: replica_section, work_section
    TYPE(thermostat_type), POINTER           :: thermostat_baro, &
                                                thermostat_part, &
                                                thermostat_shell

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY (logger, thermostat_part, thermostat_baro, npt, para_env, nhc, &
         work_section, thermostat_shell, t, averages, constant, &
         walkers_status, itimes, meta_env)
    NULLIFY (particles)
    NULLIFY (subsys)
    IF (PRESENT(md_env)) THEN
       CALL get_md_env ( md_env=md_env,&
                         thermostat_part=thermostat_part,&
                         thermostat_baro=thermostat_baro,&
                         thermostat_shell=thermostat_shell,&
                         npt=npt,&
                         t=t,&
                         constant=constant,&
                         itimes=itimes,&
                         averages=averages,&
                         para_env=para_env,&
                         error=error)
    ELSE
       IF (ASSOCIATED(force_env)) THEN
          para_env => force_env%para_env
       ELSEIF (PRESENT(pint_env)) THEN
          para_env => pint_env%logger%para_env
       ELSEIF (PRESENT(helium_env)) THEN
          para_env => helium_env%logger%para_env
       ELSE
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF

    IF(ASSOCIATED(force_env)) THEN
       meta_env => force_env%meta_env
    END IF

    IF (PRESENT(save_mem)) THEN
       my_save_mem = save_mem
    ELSE
      my_save_mem = .FALSE.
    END IF

    IF (PRESENT(write_binary_restart_file)) THEN
       my_write_binary_restart_file = write_binary_restart_file
    ELSE
       my_write_binary_restart_file = .FALSE.
    END IF

    CALL timeset(routineN//"_COUNTERS",handle2)
    IF (ASSOCIATED(itimes)) THEN
       IF (itimes>=0) THEN
          CALL section_vals_val_set(motion_section,"MD%STEP_START_VAL",i_val=itimes,error=error)
          CPAssert(ASSOCIATED(t),cp_failure_level,routineP,error,failure)
          CALL section_vals_val_set(motion_section,"MD%TIME_START_VAL",r_val=t,error=error)
       END IF
    END IF
    IF (ASSOCIATED(constant)) THEN
       CALL section_vals_val_set(motion_section,"MD%ECONS_START_VAL",r_val=constant,error=error)
    END IF
    CALL timestop(handle2)
    ! AVERAGES
    CALL timeset(routineN//"_AVERAGES",handle2)
    IF (ASSOCIATED(averages)) THEN
       IF ((averages%do_averages).AND.(averages%itimes_start/=-1)) THEN
          work_section => section_vals_get_subs_vals(motion_section,"MD%AVERAGES",error=error)
          CALL section_vals_val_set(work_section,"_SECTION_PARAMETERS_",l_val=averages%do_averages,error=error)
          work_section => section_vals_get_subs_vals(motion_section,"MD%AVERAGES%RESTART_AVERAGES",error=error)
          CALL section_vals_val_set(work_section,"ITIMES_START",i_val=averages%itimes_start,error=error)
          CALL section_vals_val_set(work_section,"AVECPU",r_val=averages%avecpu,error=error)
          CALL section_vals_val_set(work_section,"AVEHUGONIOT",r_val=averages%avehugoniot,error=error)
          CALL section_vals_val_set(work_section,"AVETEMP_BARO",r_val=averages%avetemp_baro,error=error)
          CALL section_vals_val_set(work_section,"AVEPOT",r_val=averages%avepot,error=error)
          CALL section_vals_val_set(work_section,"AVEKIN",r_val=averages%avekin,error=error)
          CALL section_vals_val_set(work_section,"AVETEMP",r_val=averages%avetemp,error=error)
          CALL section_vals_val_set(work_section,"AVEKIN_QM",r_val=averages%avekin_qm,error=error)
          CALL section_vals_val_set(work_section,"AVETEMP_QM",r_val=averages%avetemp_qm,error=error)
          CALL section_vals_val_set(work_section,"AVEVOL",r_val=averages%avevol,error=error)
          CALL section_vals_val_set(work_section,"AVECELL_A",r_val=averages%aveca,error=error)
          CALL section_vals_val_set(work_section,"AVECELL_B",r_val=averages%avecb,error=error)
          CALL section_vals_val_set(work_section,"AVECELL_C",r_val=averages%avecc,error=error)
          CALL section_vals_val_set(work_section,"AVEALPHA",r_val=averages%aveal,error=error)
          CALL section_vals_val_set(work_section,"AVEBETA",r_val=averages%avebe,error=error)
          CALL section_vals_val_set(work_section,"AVEGAMMA",r_val=averages%avega,error=error)
          CALL section_vals_val_set(work_section,"AVE_ECONS",r_val=averages%econs,error=error)
          CALL section_vals_val_set(work_section,"AVE_PRESS",r_val=averages%avepress,error=error)
          CALL section_vals_val_set(work_section,"AVE_PXX",r_val=averages%avepxx,error=error)
          ! Virial averages
          IF (ASSOCIATED(averages%virial)) THEN
             ALLOCATE(buffer(9),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = RESHAPE(averages%virial%pv_total,(/9/))
             CALL section_vals_val_set(work_section,"AVE_PV_TOT",r_vals_ptr=buffer,error=error)

             ALLOCATE(buffer(9),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = RESHAPE(averages%virial%pv_virial,(/9/))
             CALL section_vals_val_set(work_section,"AVE_PV_VIR",r_vals_ptr=buffer,error=error)

             ALLOCATE(buffer(9),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = RESHAPE(averages%virial%pv_kinetic,(/9/))
             CALL section_vals_val_set(work_section,"AVE_PV_KIN",r_vals_ptr=buffer,error=error)

             ALLOCATE(buffer(9),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = RESHAPE(averages%virial%pv_constraint,(/9/))
             CALL section_vals_val_set(work_section,"AVE_PV_CNSTR",r_vals_ptr=buffer,error=error)

             ALLOCATE(buffer(9),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = RESHAPE(averages%virial%pv_xc,(/9/))
             CALL section_vals_val_set(work_section,"AVE_PV_XC",r_vals_ptr=buffer,error=error)

             ALLOCATE(buffer(9),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = RESHAPE(averages%virial%pv_fock_4c,(/9/))
             CALL section_vals_val_set(work_section,"AVE_PV_FOCK_4C",r_vals_ptr=buffer,error=error)
          END IF
          ! Colvars averages
          IF (SIZE(averages%avecolvar)>0) THEN
             ALLOCATE(buffer(SIZE(averages%avecolvar)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = averages%avecolvar
             CALL section_vals_val_set(work_section,"AVE_COLVARS",r_vals_ptr=buffer,error=error)
          END IF
          IF (SIZE(averages%aveMmatrix)>0) THEN
             ALLOCATE(buffer(SIZE(averages%aveMmatrix)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = averages%aveMmatrix
             CALL section_vals_val_set(work_section,"AVE_MMATRIX",r_vals_ptr=buffer,error=error)
          END IF
       END IF
    END IF
    CALL timestop(handle2)

    ! PARTICLE THERMOSTAT
    CALL timeset(routineN//"_THERMOSTAT_PARTICLES",handle2)
    IF (ASSOCIATED(thermostat_part)) THEN
       IF (thermostat_part%type_of_thermostat == do_thermo_nose) THEN
          ! Restart of Nose-Hoover Thermostat for Particles
          IF (.NOT.my_write_binary_restart_file) THEN
             nhc => thermostat_part%nhc
             CALL collect_nose_restart_info(nhc,para_env,eta,veta,fnhc,mnhc,error)
             work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE",error=error)
             CALL set_template_restart(work_section,eta,veta,fnhc,mnhc,error)
          END IF
       ELSE IF (thermostat_part%type_of_thermostat == do_thermo_csvr) THEN
          ! Restart of CSVR Thermostat for Particles
          work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%CSVR",error=error)
          CALL dump_csvr_restart_info(thermostat_part%csvr, para_env, work_section, error)
       ELSE IF (thermostat_part%type_of_thermostat == do_thermo_al) THEN
          ! Restart of AD_LANGEVIN Thermostat for Particles
          work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%AD_LANGEVIN",error=error)
          CALL dump_al_restart_info(thermostat_part%al, para_env, work_section, error)
       ELSE IF (thermostat_part%type_of_thermostat == do_thermo_gle) THEN
          ! Restart of GLE Thermostat for Particles
          work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%GLE",error=error)
          CALL dump_gle_restart_info(thermostat_part%gle, para_env, work_section, error)
       END IF
    END IF
    CALL timestop(handle2)

    ! BAROSTAT - THERMOSTAT
    CALL timeset(routineN//"_BAROSTAT",handle2)
    IF (ASSOCIATED(thermostat_baro)) THEN
       IF (thermostat_baro%type_of_thermostat==do_thermo_nose) THEN
          ! Restart of Nose-Hoover Thermostat for Barostat
          nhc         => thermostat_baro%nhc
          nhc_len     = SIZE(nhc%nvt,1)
          tot_nhcneed = nhc%glob_num_nhc
          ALLOCATE(eta(tot_nhcneed*nhc_len), STAT = stat )
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(veta(tot_nhcneed*nhc_len), STAT = stat )
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(fnhc(tot_nhcneed*nhc_len), STAT = stat )
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mnhc(tot_nhcneed*nhc_len), STAT = stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          counter = 0
          DO i=1,SIZE(nhc%nvt,1)
             DO j=1,SIZE(nhc%nvt,2)
                counter = counter + 1
                eta (counter)= nhc%nvt(i,j)%eta
                veta(counter)= nhc%nvt(i,j)%v
                fnhc(counter)= nhc%nvt(i,j)%f
                mnhc(counter)= nhc%nvt(i,j)%mass
             END DO
          END DO
          work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE",error=error)
          CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)
       ELSE IF (thermostat_baro%type_of_thermostat==do_thermo_csvr) THEN
          ! Restart of CSVR Thermostat for Barostat
          work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%CSVR",error=error)
          CALL dump_csvr_restart_info(thermostat_baro%csvr, para_env, work_section, error)
       END IF
    END IF
    CALL timestop(handle2)

    ! BAROSTAT
    CALL timeset(routineN//"_NPT",handle2)
    IF (ASSOCIATED(npt)) THEN
       ALLOCATE(veta(SIZE(npt,1)*SIZE(npt,2)), STAT = stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mnhc(SIZE(npt,1)*SIZE(npt,2)), STAT = stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       counter = 0
       DO i=1,SIZE(npt,1)
          DO j=1,SIZE(npt,2)
             counter = counter + 1
             veta(counter) = npt(i,j)%v
             mnhc(counter) = npt(i,j)%mass
          END DO
       END DO
       work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT",error=error)
       CALL set_template_restart(work_section, veta=veta, mnhc=mnhc, error=error)
    END IF
    CALL timestop(handle2)

    ! SHELL THERMOSTAT
    CALL timeset(routineN//"_THERMOSTAT_SHELL",handle2)
    IF (ASSOCIATED(thermostat_shell)) THEN
       IF (thermostat_shell%type_of_thermostat == do_thermo_nose) THEN
          ! Restart of Nose-Hoover Thermostat for Shell Particles
          IF (.NOT.my_write_binary_restart_file) THEN
             nhc => thermostat_shell%nhc
             CALL collect_nose_restart_info(nhc,para_env,eta,veta,fnhc,mnhc,error)
             work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE",error=error)
             CALL set_template_restart(work_section,eta,veta,fnhc,mnhc,error)
          END IF
       ELSE IF (thermostat_shell%type_of_thermostat == do_thermo_csvr) THEN
          work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%CSVR",error=error)
          ! Restart of CSVR Thermostat for Shell Particles
          CALL dump_csvr_restart_info(thermostat_shell%csvr, para_env, work_section, error)
       END IF
    END IF
    CALL timestop(handle2)

    CALL timeset(routineN//"_META",handle2)
    IF (ASSOCIATED(meta_env)) THEN
       CALL section_vals_val_set(meta_env%metadyn_section,"STEP_START_VAL",&
            i_val=meta_env%n_steps,error=error)
       CALL section_vals_val_set(meta_env%metadyn_section,"NHILLS_START_VAL",&
            i_val=meta_env%hills_env%n_hills,error=error)
       !RG Adaptive hills
       CALL section_vals_val_set(meta_env%metadyn_section,"MIN_DISP",&
            r_val=meta_env%hills_env%min_disp,error=error)
       CALL section_vals_val_set(meta_env%metadyn_section,"OLD_HILL_NUMBER",&
            i_val=meta_env%hills_env%old_hill_number,error=error)
       CALL section_vals_val_set(meta_env%metadyn_section,"OLD_HILL_STEP",&
            i_val=meta_env%hills_env%old_hill_step,error=error)
       !RG Adaptive hills
       IF (meta_env%do_hills.AND.meta_env%hills_env%n_hills/=0) THEN
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_POS",error=error)
          CALL meta_hills_val_set_ss(work_section, meta_env, error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_SCALE",error=error)
          CALL meta_hills_val_set_ds(work_section, meta_env, error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_HEIGHT",error=error)
          CALL meta_hills_val_set_ww(work_section, meta_env, error=error)
          IF (meta_env%well_tempered) THEN
             work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_INVDT",error=error)
             CALL meta_hills_val_set_dt(work_section, meta_env, error=error)
          END IF
       END IF
       IF (meta_env%extended_lagrange) THEN
          CALL section_vals_val_set(meta_env%metadyn_section,"COLVAR_AVG_TEMPERATURE_RESTART",&
               r_val=meta_env%avg_temp,error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS0",error=error)
          DO irep = 1, meta_env%n_colvar
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_val=meta_env%metavar(irep)%ss0,&
                  i_rep_val=irep,error=error)
          END DO
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_VVP",error=error)
          DO irep = 1, meta_env%n_colvar
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_val=meta_env%metavar(irep)%vvp,&
                  i_rep_val=irep,error=error)
          END DO

          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS",error=error)
          DO irep = 1, meta_env%n_colvar
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_val=meta_env%metavar(irep)%ss,&
                  i_rep_val=irep,error=error)
          END DO
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_FS",error=error)
          DO irep = 1, meta_env%n_colvar
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_val=meta_env%metavar(irep)%ff_s,&
                  i_rep_val=irep,error=error)
          END DO

       END IF
       ! Multiple Walkers
       IF (meta_env%do_multiple_walkers) THEN
          ALLOCATE(walkers_status(meta_env%multiple_walkers%walkers_tot_nr),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          walkers_status = meta_env%multiple_walkers%walkers_status
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"MULTIPLE_WALKERS",error=error)
          CALL section_vals_val_set(work_section,"WALKERS_STATUS",i_vals_ptr=walkers_status,error=error)
       END IF
    END IF
    CALL timestop(handle2)
    CALL timeset(routineN//"_NEB",handle2)
    IF (PRESENT(coords).OR.(PRESENT(vels))) THEN
       ! Update NEB section
       replica_section => section_vals_get_subs_vals(motion_section,"BAND%REPLICA",error=error)
       CALL force_env_get(force_env,subsys=subsys,error=error)
       CALL cp_subsys_get(subsys,particles=particles,error=error)
       IF (PRESENT(coords)) THEN
          ! Allocate possible missing sections
          DO
             IF (coords%size_wrk(2)<=SIZE(replica_section%values,2)) EXIT
             CALL section_vals_add_values(replica_section,error=error)
          END DO
          ! Write Values
          DO isec = 1, coords%size_wrk(2)
             CALL section_vals_val_unset(replica_section,"COORD_FILE_NAME",i_rep_section=isec,error=error)
             work_section => section_vals_get_subs_vals3(replica_section,"COORD",i_rep_section=isec,error=error)
             CALL section_neb_coord_val_set(work_section,coords%xyz(:,isec),SIZE(coords%xyz,1),3*SIZE(particles%els),&
                                            3,particles%els,angstrom,error)
             ! Update Collective Variables
             IF (coords%in_use == do_band_collective) THEN
                ALLOCATE (wrk(coords%size_wrk(1)),stat=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                wrk = coords%wrk(:,isec)
                CALL section_vals_val_set(replica_section,"COLLECTIVE",r_vals_ptr=wrk,&
                     i_rep_section=isec,error=error)
             END IF
          END DO
       END IF
       IF (PRESENT(vels)) THEN
          CALL force_env_get(force_env,subsys=subsys,error=error)
          CALL cp_subsys_get(subsys,particles=particles,error=error)
          ! Allocate possible missing sections
          DO
             IF (vels%size_wrk(2)<=SIZE(replica_section%values,2)) EXIT
             CALL section_vals_add_values(replica_section,error=error)
          END DO
          ! Write Values
          DO isec = 1, vels%size_wrk(2)
             work_section => section_vals_get_subs_vals3(replica_section,"VELOCITY",i_rep_section=isec,error=error)
             IF (vels%in_use == do_band_collective) THEN
                CALL section_neb_coord_val_set(work_section,vels%wrk(:,isec),SIZE(vels%wrk,1),SIZE(vels%wrk,1),&
                                               1,particles%els,1.0_dp,error)
             ELSE
                CALL section_neb_coord_val_set(work_section,vels%wrk(:,isec),SIZE(vels%wrk,1),3*SIZE(particles%els),&
                                               3,particles%els,1.0_dp,error)
             END IF
          END DO
       END IF
    END IF
    CALL timestop(handle2)

    IF (PRESENT(pint_env)) THEN
      ! Update PINT section
      CALL update_motion_pint(motion_section, pint_env, error)
    END IF

    IF (PRESENT(helium_env)) THEN
      ! Update HELIUM section
      CALL update_motion_helium(motion_section, helium_env, error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE update_motion


! ***************************************************************************
!> \brief  Update PINT section in the input structure
!> \param motion_section ...
!> \param pint_env ...
!> \param error ...
!> \date   2010-10-13
!> \author Lukasz Walewski <Lukasz.Walewski@ruhr-uni-bochum.de>
! *****************************************************************************
  SUBROUTINE update_motion_pint(motion_section, pint_env, error)

    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_motion_pint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ibead, &
                                                inos, istat
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r_vals
    TYPE(section_vals_type), POINTER         :: pint_section, tmpsec

    CALL timeset(routineN,handle)

    pint_section => section_vals_get_subs_vals(motion_section,"PINT",error=error)
    CALL section_vals_val_set(pint_section,"ITERATION",i_val=pint_env%iter,error=error)

    ! allocate memory for COORDs and VELOCITYs if the BEADS section was not
    ! explicitly given in the input (this is actually done only once since
    ! after section_vals_add_values section becomes explict)
    NULLIFY(tmpsec)
    tmpsec => section_vals_get_subs_vals(pint_section,"BEADS", error=error)
    CALL section_vals_get(tmpsec,explicit=explicit,error=error)
    IF ( .NOT. explicit ) THEN
      CALL section_vals_add_values(tmpsec,error)
    END IF

    ! update bead coordinates in the global input structure
    NULLIFY(r_vals)
    ALLOCATE(r_vals(pint_env%p*pint_env%ndim),STAT=istat)
    CPAssert(istat==0, cp_failure_level, routineP, error, failure)

    i=1
    CALL pint_u2x(pint_env,error=error)
    DO iatom=1, pint_env%ndim
      DO ibead=1, pint_env%p
        r_vals(i)=pint_env%x(ibead,iatom)
        i = i+1
      END DO
    END DO
    CALL section_vals_val_set(pint_section,"BEADS%COORD%_DEFAULT_KEYWORD_", &
         r_vals_ptr=r_vals,error=error)

    ! update bead velocities in the global input structure
    NULLIFY(r_vals)
    ALLOCATE(r_vals(pint_env%p*pint_env%ndim),STAT=istat)
    CPAssert(istat==0, cp_failure_level, routineP, error, failure)
    i=1
    CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
    DO iatom=1, pint_env%ndim
      DO ibead=1, pint_env%p
        r_vals(i)=pint_env%v(ibead,iatom)
        i = i+1
      END DO
    END DO
    CALL section_vals_val_set(pint_section,"BEADS%VELOCITY%_DEFAULT_KEYWORD_",&
         r_vals_ptr=r_vals,error=error)

    IF(pint_env%pimd_thermostat==thermostat_nose) THEN

      ! allocate memory for COORDs and VELOCITYs if the NOSE section was not
      ! explicitly given in the input (this is actually done only once since
      ! after section_vals_add_values section becomes explict)
      NULLIFY(tmpsec)
      tmpsec => section_vals_get_subs_vals(pint_section,"NOSE", error=error)
      CALL section_vals_get(tmpsec,explicit=explicit,error=error)
      IF ( .NOT. explicit ) THEN
        CALL section_vals_add_values(tmpsec,error)
      END IF

      ! update thermostat coordinates in the global input structure
      NULLIFY(r_vals)
      ALLOCATE(r_vals(pint_env%p*pint_env%ndim*pint_env%nnos),STAT=istat)
      CPAssert(istat==0, cp_failure_level, routineP, error, failure)
      i=1
      DO iatom=1, pint_env%ndim
        DO ibead=1, pint_env%p
          DO inos=1, pint_env%nnos
            r_vals(i)=pint_env%tx(inos,ibead,iatom)
            i = i+1
          END DO
        END DO
      END DO
      CALL section_vals_val_set(pint_section,"NOSE%COORD%_DEFAULT_KEYWORD_", &
           r_vals_ptr=r_vals,error=error)

      ! update thermostat velocities in the global input structure
      NULLIFY(r_vals)
      ALLOCATE(r_vals(pint_env%p*pint_env%ndim*pint_env%nnos),STAT=istat)
      CPAssert(istat==0, cp_failure_level, routineP, error, failure)
      i=1
      DO iatom=1, pint_env%ndim
        DO ibead=1, pint_env%p
          DO inos=1, pint_env%nnos
            r_vals(i)=pint_env%tv(inos,ibead,iatom)
            i = i+1
          END DO
        END DO
      END DO
      CALL section_vals_val_set(pint_section,"NOSE%VELOCITY%_DEFAULT_KEYWORD_", &
           r_vals_ptr=r_vals,error=error)

    ELSEIF(pint_env%pimd_thermostat==thermostat_gle) THEN

      NULLIFY(tmpsec)
      tmpsec => section_vals_get_subs_vals(pint_section,"GLE",error=error)
      CALL dump_gle_restart_info(pint_env%gle, pint_env%replicas%para_env, tmpsec, error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE update_motion_pint

  ! ***************************************************************************
  !> \brief  Update HELIUM section in the input structure.
  !> \date   2009-11-12
  !> \author Lukasz Walewski <Lukasz.Walewski@ruhr-uni-bochum.de>
  !> \descr  Transfer the current helium state from the runtime environment
  !>         to the input structure, so that it can be used for I/O, etc.
  !> \note   Moved from the helium_io module directly, might be done better way
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param motion_section ...
!> \param helium ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE update_motion_helium(motion_section, helium, error)

    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_motion_helium', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: err_str, stmp
    INTEGER                                  :: handle, i, itmp, iweight, &
                                                msglen, nsteps, offset, &
                                                reqlen, status
    INTEGER, DIMENSION(:), POINTER           :: int_msg_gather
    LOGICAL                                  :: explicit, failure, lbf
    REAL(kind=dp)                            :: bf, bu, rtmp
    REAL(kind=dp), DIMENSION(3, 2)           :: bg, cg, ig
    REAL(kind=dp), DIMENSION(:), POINTER     :: real_msg, real_msg_gather
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: message
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: tmpsec

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(helium),cp_failure_level,routineP,error,failure)

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF ( .NOT. helium%solute_present ) THEN
      ! update iteration number
      itmp = logger%iter_info%iteration(2)
      CALL section_vals_val_set( &
        helium%input, &
        "MOTION%PINT%ITERATION", &
        i_val=itmp, &
        error=error )
      ! else - PINT will do that
    END IF

    !
    ! save coordinates
    !
    ! allocate the buffer to be passed and fill it with local coords at each proc
    NULLIFY(real_msg)
    msglen = SIZE(helium%pos)
    ALLOCATE(real_msg(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    real_msg(:) = PACK( helium%pos, .TRUE. )

    ! allocate the buffer for message passing
    NULLIFY(real_msg_gather)
    msglen = SIZE(helium%pos) * logger%para_env%num_pe
    ALLOCATE(real_msg_gather(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    ! pass the message from all processors to logger%para_env%source
    real_msg_gather(:) = 0.0_dp
    CALL mp_gather(real_msg,real_msg_gather,&
         logger%para_env%source,logger%para_env%group)

    ! update coordinates in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%COORD%_DEFAULT_KEYWORD_", &
         r_vals_ptr=real_msg_gather,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(real_msg_gather)

    ! DEALLOCATE since this array is only used locally
    DEALLOCATE(real_msg,STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    !
    ! save permutation state
    !
    ! allocate the buffer for message passing
    NULLIFY(int_msg_gather)
    msglen = SIZE(helium%permutation) * logger%para_env%num_pe
    ALLOCATE(int_msg_gather(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    ! pass the message from all processors to logger%para_env%source
    int_msg_gather(:) = 0
    CALL mp_gather(helium%permutation,int_msg_gather,&
         logger%para_env%source,logger%para_env%group)

    ! update permutation state in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%PERM%_DEFAULT_KEYWORD_", &
         i_vals_ptr=int_msg_gather,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(int_msg_gather)

    !
    ! save RNG state
    !
    ! pack RNG state on each processor to the local array
    NULLIFY(real_msg)
    msglen = 40
    ALLOCATE(real_msg(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    CALL get_rng_stream(helium%rng_stream_uniform,bg=bg,cg=cg,ig=ig,&
         buffer=bu,buffer_filled=lbf,error=error)
    offset = 0
    real_msg(offset+1:offset+6)   = PACK( bg, .TRUE. )
    real_msg(offset+7:offset+12)  = PACK( cg, .TRUE. )
    real_msg(offset+13:offset+18) = PACK( ig, .TRUE. )
    IF ( lbf ) THEN
      bf = 1.0_dp
    ELSE
      bf = -1.0_dp
    END IF
    real_msg(offset+19) = bf
    real_msg(offset+20) = bu
    CALL get_rng_stream(helium%rng_stream_gaussian,bg=bg,cg=cg,ig=ig,&
         buffer=bu,buffer_filled=lbf,error=error)
    offset = 20
    real_msg(offset+1:offset+6)   = PACK( bg, .TRUE. )
    real_msg(offset+7:offset+12)  = PACK( cg, .TRUE. )
    real_msg(offset+13:offset+18) = PACK( ig, .TRUE. )
    IF ( lbf ) THEN
      bf = 1.0_dp
    ELSE
      bf = -1.0_dp
    END IF
    real_msg(offset+19) = bf
    real_msg(offset+20) = bu

    ! Gather RNG state (in real_msg_gather vector) from all processors at
    ! logger%para_env%source
    NULLIFY(real_msg_gather)
    msglen = SIZE(real_msg)*logger%para_env%num_pe
    ALLOCATE(real_msg_gather(msglen),STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
    real_msg_gather(:) = 0.0_dp
    CALL mp_gather(real_msg,real_msg_gather,logger%para_env%source,logger%para_env%group)

    ! update the RNG state in the global input structure
    CALL section_vals_val_set(helium%input, &
         "MOTION%PINT%HELIUM%RNG_STATE%_DEFAULT_KEYWORD_", &
         r_vals_ptr=real_msg_gather,error=error)

    ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
    ! assigned in section_vals_val_set - this memeory will be used later on!
    ! "The val becomes the owner of the array" - from section_vals_val_set docu
    NULLIFY(real_msg_gather)

    ! DEALLOCATE since this array is only used locally
    DEALLOCATE(real_msg,STAT=status)
    CPPostcondition(status==0,cp_failure_level,routineP,error,failure)

    IF (helium%solute_present) THEN
      !
      ! save forces on the solute
      !
      ! check that the number of values match the current runtime
      reqlen = helium%solute_atoms * helium%solute_beads * 3
      msglen = SIZE(helium%force_avrg)
      err_str = "Invalid size of HELIUM%FORCE: received '"
      stmp = ""
      WRITE(stmp,*) msglen
      err_str = TRIM(ADJUSTL(err_str)) // &
                TRIM(ADJUSTL(stmp)) // "' but expected '"
      stmp = ""
      WRITE(stmp,*) reqlen
      err_str = TRIM(ADJUSTL(err_str)) // &
                TRIM(ADJUSTL(stmp)) // "'."
      CALL cp_assert(msgLEN==reqlen,cp_failure_level,&
           cp_assertion_failed, routineP, err_str)

      ! allocate the buffer to be saved and fill it with forces
      ! forces should be the same on all processors, but we don't check that here
      NULLIFY(real_msg_gather)
      ALLOCATE(real_msg_gather(msglen),STAT=status)
      CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
      real_msg_gather(:) = PACK( helium%force_avrg, .TRUE. )

      ! update forces in the global input structure
      CALL section_vals_val_set(helium%input, &
           "MOTION%PINT%HELIUM%FORCE%_DEFAULT_KEYWORD_", &
           r_vals_ptr=real_msg_gather,error=error)

      ! NULLIFY, but do not DEALLOCATE! - a new pointer to this array is silently
      ! assigned in section_vals_val_set - this memeory will be used later on!
      ! "The val becomes the owner of the array" - from section_vals_val_set docu
      NULLIFY(real_msg_gather)
    END IF

    !
    ! save the densities
    !
    IF (helium%rho_present) THEN

      ! make sure that the section is explicitly present in the input structure
      NULLIFY(tmpsec)
      tmpsec => section_vals_get_subs_vals(helium%input, &
                "MOTION%PINT%HELIUM%RHO", &
                error=error)
      CALL section_vals_get(tmpsec,explicit=explicit,error=error)
      IF ( .NOT. explicit ) THEN
        CALL section_vals_add_values(tmpsec,error)
      END IF

      ! work on the temporary array so that accumulated data remains intact
      helium%rho_inst(:,:,:,:) = helium%rho_avrg(:,:,:,:)

      ! average all the density estimators over helium environments
      DO i = 1, helium%rho_num
        NULLIFY(message)
        message => helium%rho_inst(i,1::1,1::1,1::1)
        CALL mp_sum(message,logger%para_env%group)
      END DO
      itmp = logger%para_env%num_pe
      CPPostcondition(itmp>0,cp_failure_level,routineP,error,failure)
      rtmp = 1.0_dp / REAL(itmp,dp)
      helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) * rtmp

      ! average over steps performed so far in this run
      nsteps = helium%current_step-helium%first_step
      CPPostcondition(nsteps>0,cp_failure_level,routineP,error,failure)
      rtmp = 1.0_dp / REAL(nsteps,dp)
      helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) * rtmp

      iweight = helium%rho_iweight
      rtmp = 1.0_dp / REAL(nsteps+iweight,dp)
      IF ( helium%rho_restart ) THEN
        ! average over the old and the current density (observe the weights!)
        helium%rho_inst(:,:,:,:) = nsteps * helium%rho_inst(:,:,:,:) + &
                                   iweight * helium%rho_rstr(:,:,:,:)
        helium%rho_inst(:,:,:,:) = helium%rho_inst(:,:,:,:) * rtmp
      END IF

      ! update the densities in the global input structure
      NULLIFY(real_msg)
      msglen = SIZE(helium%rho_inst)
      ALLOCATE(real_msg(msglen),STAT=status)
      CPPostcondition(status==0,cp_failure_level,routineP,error,failure)
      real_msg(:) = PACK( helium%rho_inst, .TRUE. )
      CALL section_vals_val_set(helium%input, &
           "MOTION%PINT%HELIUM%RHO%CUBE_DATA%_DEFAULT_KEYWORD_", &
           r_vals_ptr=real_msg,error=error)
      NULLIFY(real_msg)

      ! update the weighting factor
      itmp = helium%rho_iweight
      IF ( itmp .LT. 0 ) THEN
        itmp = helium%current_step-helium%first_step
      ELSE
        itmp = itmp + helium%current_step-helium%first_step
      END IF
      CALL section_vals_val_set(helium%input, &
           "MOTION%PINT%HELIUM%RHO%IWEIGHT", &
           i_val=itmp,error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE update_motion_helium

! *****************************************************************************
!> \brief routine to dump thermostat CSVR energies
!> \param thermostat_energy ...
!> \param nsize ...
!> \param work_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [teo]
!> \author Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE dump_csvr_energy_info(thermostat_energy, nsize, work_section, error)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: thermostat_energy
    INTEGER, INTENT(IN)                      :: nsize
    TYPE(section_vals_type), POINTER         :: work_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_csvr_energy_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(work_section),cp_failure_level,routineP,error,failure)
    CPPrecondition((work_section%ref_count>0),cp_failure_level,routineP,error,failure)

    NULLIFY (my_val,old_val,section,vals)

    section => work_section%section

    ik = section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)

    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
                   "section "//TRIM(section%name)//" does not contain keyword "//&
                   "_DEFAULT_KEYWORD_",error,failure)

    DO
       IF (SIZE(work_section%values,2)==1) EXIT
       CALL section_vals_add_values(work_section,error=error)
    END DO

    vals => work_section%values(ik,1)%list
    Nlist = 0

    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF

    DO irk=1,nsize
       CALL val_create(val=my_val,r_val=thermostat_energy(irk),error=error)
       IF (Nlist /= 0) THEN
          IF (irk == 1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk == 1) THEN
             NULLIFY (new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY (new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY (my_val)
    END DO
    work_section%values(ik,1)%list => vals

  END SUBROUTINE dump_csvr_energy_info

! *****************************************************************************
!> \brief Collect all information needed to dump the restart for CSVR
!>      thermostat
!> \param csvr ...
!> \param para_env ...
!> \param csvr_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino] - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE dump_csvr_restart_info(csvr, para_env, csvr_section, error)
    TYPE(csvr_system_type), POINTER          :: csvr
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: csvr_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_csvr_restart_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: i, my_index, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: dwork
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: thermo_energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    ! Thermostat Energies
    ALLOCATE(work(csvr%glob_num_csvr),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(thermo_energy(csvr%loc_num_csvr),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1,csvr%loc_num_csvr
       thermo_energy(i) = csvr%nvt(i)%thermostat_energy
    END DO
    CALL get_kin_energies (csvr%map_info,csvr%loc_num_csvr,&
               csvr%glob_num_csvr,thermo_energy,&
               dum, para_env, array_kin=work, error=error)
    DEALLOCATE(thermo_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! If check passes then let's dump the info on the restart file
    work_section => section_vals_get_subs_vals(csvr_section,"THERMOSTAT_ENERGY",error=error)
    CALL dump_csvr_energy_info(work, csvr%glob_num_csvr, work_section, error)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Thermostat Random Number info for restart
    work_section => section_vals_get_subs_vals(csvr_section,"RNG_INIT",error=error)
    ALLOCATE (dwork(rng_record_length,csvr%glob_num_csvr),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dwork = 0
    DO i=1,csvr%loc_num_csvr
       my_index = csvr%map_info%index(i)
       CALL dump_rng_stream(rng_stream=csvr%nvt(i)%gaussian_rng_stream,&
            rng_record=rng_record, error=error)
       CALL string_to_ascii(rng_record,dwork(:,my_index))
    END DO

    !  Collect data if there was no communication in this thermostat
    IF (csvr%map_info%dis_type==do_thermo_no_communication) THEN
       ! Collect data if there was no communication in this thermostat
       CALL mp_sum(dwork,para_env%group)
    ELSE
       ! Perform some check and collect data in case of communicating thermostats
       CALL communication_thermo_low2(dwork, rng_record_length, csvr%glob_num_csvr, para_env, error)
    END IF
    CALL section_rng_val_set(rng_section=work_section,nsize=csvr%glob_num_csvr,ascii=dwork,&
         error=error)
    DEALLOCATE (dwork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE dump_csvr_restart_info

! *****************************************************************************
!> \brief Collect all information needed to dump the restart for AD_LANGEVIN
!>      thermostat
!> \param al ...
!> \param para_env ...
!> \param al_section ...
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      10.2007 created [tlaino] - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE dump_al_restart_info(al, para_env, al_section, error)
    TYPE(al_system_type), POINTER            :: al
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: al_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_al_restart_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dum
    REAL(KIND=dp), DIMENSION(:), POINTER     :: t_array, work
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.    

    ! chi and mass
    ALLOCATE(work(al%glob_num_al),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(t_array(al%loc_num_al),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! copy chi into temporary
    DO i = 1,al%loc_num_al
       t_array(i) = al%nvt(i)%chi
    END DO
    ! consolidate into work
    CALL get_kin_energies (al%map_info,al%loc_num_al,&
               al%glob_num_al,t_array,&
               dum, para_env, array_kin=work, error=error)

    ! If check passes then let's dump the info on the restart file
    work_section => section_vals_get_subs_vals(al_section,"CHI",error=error)
    CALL dump_csvr_energy_info(work, al%glob_num_al, work_section, error)

    ! copy mass into temporary
    DO i = 1,al%loc_num_al
       t_array(i) = al%nvt(i)%mass
    END DO
    ! consolidate into work
    CALL get_kin_energies (al%map_info,al%loc_num_al,&
               al%glob_num_al,t_array,&
               dum, para_env, array_kin=work, error=error)

    ! If check passes then let's dump the info on the restart file
    work_section => section_vals_get_subs_vals(al_section,"MASS",error=error)
    CALL dump_csvr_energy_info(work, al%glob_num_al, work_section, error)

    DEALLOCATE (t_array,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE dump_al_restart_info

! *****************************************************************************
!> \brief Collect all information needed to dump the restart for GLE
!>      thermostat
!> \param gle ...
!> \param para_env ...
!> \param gle_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author MI
! *****************************************************************************
  SUBROUTINE dump_gle_restart_info(gle, para_env, gle_section, error)
    TYPE(gle_type), POINTER                  :: gle
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: gle_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_gle_restart_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: counter, glob_num, i, iproc, &
                                                j, loc_num, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: dwork
    INTEGER, DIMENSION(:), POINTER           :: gle_per_proc, index
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dum
    REAL(dp), DIMENSION(:), POINTER          :: s_tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: thermo_energy
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.

    ! Thermostat Energies
    ALLOCATE(work(gle%glob_num_gle),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(thermo_energy(gle%loc_num_gle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1,gle%loc_num_gle
       thermo_energy(i) = gle%nvt(i)%thermostat_energy
    END DO
    CALL get_kin_energies (gle%map_info,gle%loc_num_gle,&
               gle%glob_num_gle,thermo_energy,&
               dum, para_env, array_kin=work, error=error)
    DEALLOCATE(thermo_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! If check passes then let's dump the info on the restart file
    work_section => section_vals_get_subs_vals(gle_section,"THERMOSTAT_ENERGY",error=error)
    CALL dump_csvr_energy_info(work, gle%glob_num_gle, work_section, error)
    DEALLOCATE (work,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Thermostat Random Number info for restart
    work_section => section_vals_get_subs_vals(gle_section,"RNG_INIT",error=error)
    glob_num = gle%glob_num_gle
    loc_num = gle%loc_num_gle
    ALLOCATE (dwork(rng_record_length,glob_num),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dwork = 0
    DO i=1,loc_num
       j = gle%map_info%index(i)
       CALL dump_rng_stream(rng_stream= gle%nvt(i)%gaussian_rng_stream,&
            rng_record=rng_record, error=error)
       CALL string_to_ascii(rng_record,dwork(:,j))
    END DO

    !  Collect data if there was no communication in this thermostat
    IF (gle%map_info%dis_type==do_thermo_no_communication) THEN
       ! Collect data if there was no communication in this thermostat
       CALL mp_sum(dwork,para_env%group)
    ELSE
       ! Perform some check and collect data in case of communicating thermostats
       CALL communication_thermo_low2(dwork, rng_record_length, glob_num, para_env, error)
    END IF
    CALL section_rng_val_set(rng_section=work_section,nsize=glob_num,ascii=dwork,&
         error=error)
    DEALLOCATE (dwork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE ( gle_per_proc(para_env%num_pe),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    gle_per_proc(:) = 0
    CALL mp_allgather(gle%loc_num_gle,gle_per_proc,para_env%group)

    ! Thermostat S variable info for restart
    NULLIFY(s_tmp)
    ALLOCATE (s_tmp((gle%ndim)*gle%glob_num_gle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    s_tmp=0.0_dp

    NULLIFY(work,index)
    DO iproc=1,para_env%num_pe
       CALL reallocate(work,1,gle_per_proc(iproc)*(gle%ndim))
       CALL reallocate(index,1,gle_per_proc(iproc))
       IF(para_env%mepos == (iproc-1)) THEN
          INDEX(:)=0
          counter=0
          DO i=1,gle%ndim
             DO j=1,gle%loc_num_gle
                counter=counter+1
                work(counter)=gle%nvt(j)%s(i)
                INDEX(j) = gle%map_info%index(j)
             END DO
          END DO
       ELSE
          work(:) = 0.0_dp
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       CALL mp_bcast(index,iproc-1,para_env%group)
       counter=0
       DO i=1,gle%ndim
          DO j=1,gle_per_proc(iproc)
             counter = counter+1
             s_tmp((INDEX(j)-1)*(gle%ndim)+i) = work(counter)
          END DO
       END DO
    END DO

    IF(SIZE(s_tmp)>0) THEN
    work_section => section_vals_get_subs_vals(gle_section,"S",error=error)
    CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_vals_ptr=s_tmp,error=error)
    ELSE
      DEALLOCATE(s_tmp,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(gle_per_proc, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( work, STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE ( index, STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE dump_gle_restart_info

! *****************************************************************************
!> \brief Collect all information needed to dump the restart for Nose-Hoover
!>      thermostat
!> \param nhc ...
!> \param para_env ...
!> \param eta ...
!> \param veta ...
!> \param fnhc ...
!> \param mnhc ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino] - University of Zurich
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE collect_nose_restart_info(nhc, para_env, eta, veta, fnhc, mnhc, error)
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eta, veta, fnhc, mnhc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'collect_nose_restart_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: counter, i, iproc, j, &
                                                nhc_len, num_nhc, numneed, &
                                                stat, tot_nhcneed
    INTEGER, DIMENSION(:), POINTER           :: index, nhc_per_proc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(map_info_type), POINTER             :: map_info

    nhc_len = SIZE(nhc%nvt,1)
    num_nhc = nhc%loc_num_nhc
    numneed = num_nhc
    map_info=> nhc%map_info
    ALLOCATE (nhc_per_proc(para_env%num_pe), STAT = stat )
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "nhc_per_proc",int_size*para_env%num_pe)
    nhc_per_proc(:) = 0

    CALL mp_allgather(numneed,nhc_per_proc,para_env%group)
    tot_nhcneed = nhc%glob_num_nhc

    NULLIFY(work,index)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! nhc%eta
    !-----------------------------------------------------------------------------
    ALLOCATE (eta(tot_nhcneed*nhc_len),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "eta",dp_size*tot_nhcneed*nhc_len)
    DO iproc=1,para_env%num_pe
       CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
       CALL reallocate(index,1,nhc_per_proc(iproc))
       IF(para_env%mepos == (iproc-1)) THEN
          INDEX(:)=0
          counter=0
          DO i=1,nhc_len
             DO j=1,num_nhc
                counter=counter+1
                work(counter)=nhc%nvt(i,j)%eta
                INDEX(j) = map_info%index(j)
             END DO
          END DO
       ELSE
          work(:) = 0.0_dp
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       CALL mp_bcast(index,iproc-1,para_env%group)
       counter=0
       DO i=1,nhc_len
          DO j=1,nhc_per_proc(iproc)
             counter = counter+1
             eta((INDEX(j)-1)*nhc_len+i) = work(counter)
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! nhc%veta
    !-----------------------------------------------------------------------------
    ALLOCATE (veta(tot_nhcneed*nhc_len),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "veta",dp_size*tot_nhcneed*nhc_len)
    DO iproc=1,para_env%num_pe
       CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
       CALL reallocate(index,1,nhc_per_proc(iproc))
       IF(para_env%mepos == (iproc-1)) THEN
          INDEX(:)=0
          counter=0
          DO i=1,nhc_len
             DO j=1,num_nhc
                counter=counter+1
                work(counter)=nhc%nvt(i,j)%v
                INDEX(j) = map_info%index(j)
             END DO
          END DO
       ELSE
          work(:) = 0.0_dp
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       CALL mp_bcast(index,iproc-1,para_env%group)
       counter=0
       DO i=1,nhc_len
          DO j=1,nhc_per_proc(iproc)
             counter = counter+1
             veta((INDEX(j)-1)*nhc_len+i) = work(counter)
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! nhc%force
    !-----------------------------------------------------------------------------
    ALLOCATE (fnhc(tot_nhcneed*nhc_len),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "fnhc",dp_size*tot_nhcneed*nhc_len)
    DO iproc=1,para_env%num_pe
       CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
       CALL reallocate(index,1,nhc_per_proc(iproc))
       IF(para_env%mepos == (iproc-1)) THEN
          INDEX(:)=0
          counter=0
          DO i=1,nhc_len
             DO j=1,num_nhc
                counter=counter+1
                work(counter)=nhc%nvt(i,j)%f
                INDEX(j) = map_info%index(j)
             END DO
          END DO
       ELSE
          work(:) = 0.0_dp
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       CALL mp_bcast(index,iproc-1,para_env%group)
       counter=0
       DO i=1,nhc_len
          DO j=1,nhc_per_proc(iproc)
             counter = counter+1
             fnhc((INDEX(j)-1)*nhc_len+i) = work(counter)
          END DO
       END DO
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! nhc%mass
    !-----------------------------------------------------------------------------
    ALLOCATE (mnhc(tot_nhcneed*nhc_len),STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "mnhc",dp_size*tot_nhcneed*nhc_len)
    DO iproc=1,para_env%num_pe
       CALL reallocate(work,1,nhc_per_proc(iproc)*nhc_len)
       CALL reallocate(index,1,nhc_per_proc(iproc))
       IF(para_env%mepos == (iproc-1)) THEN
          INDEX(:)=0
          counter=0
          DO i=1,nhc_len
             DO j=1,num_nhc
                counter=counter+1
                work(counter)=nhc%nvt(i,j)%mass
                INDEX(j) = map_info%index(j)
             END DO
          END DO
       ELSE
          work(:) = 0.0_dp
       END IF
       CALL mp_bcast(work,iproc-1,para_env%group)
       CALL mp_bcast(index,iproc-1,para_env%group)
       counter=0
       DO i=1,nhc_len
          DO j=1,nhc_per_proc(iproc)
             counter = counter+1
             mnhc((INDEX(j)-1)*nhc_len+i) = work(counter)
          END DO
       END DO
    END DO

    DEALLOCATE (work,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"work")
    DEALLOCATE (index,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"index")
    DEALLOCATE (nhc_per_proc,STAT=stat)
    IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"nhc_per_proc")

  END SUBROUTINE collect_nose_restart_info

! *****************************************************************************
!> \brief routine to dump NEB coordinates and velocities section.. fast implementation
!> \param coord_section ...
!> \param array ...
!> \param narray ...
!> \param nsize ...
!> \param nfield ...
!> \param particle_set ...
!> \param conv_factor ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_neb_coord_val_set(coord_section,array,narray,nsize,nfield,&
                                       particle_set,conv_factor,error)

    TYPE(section_vals_type), POINTER         :: coord_section
    REAL(KIND=dp), DIMENSION(*)              :: array
    INTEGER, INTENT(IN)                      :: narray, nsize, nfield
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'section_neb_coord_val_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_c
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure = .FALSE.
    NULLIFY (my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(coord_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(coord_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => coord_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(coord_section%values,2)==1) EXIT
       CALL section_vals_add_values(coord_section,error=error)
    END DO
    vals  => coord_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,nsize/nfield
       ALLOCATE (my_c(nfield),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (nfield == 3) THEN
          my_c(1:3) = get_particle_pos_or_vel(irk,particle_set,array(1:narray))
          my_c(1:3) = my_c(1:3)*conv_factor
       ELSE
          my_c(1) = array(irk)
       END IF
       CALL val_create(my_val,r_vals_ptr=my_c,error=error)

       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    coord_section%values(ik,1)%list => vals
  END SUBROUTINE section_neb_coord_val_set

! *****************************************************************************
!> \brief Set the nose structure like restart
!> \param work_section ...
!> \param eta ...
!> \param veta ...
!> \param fnhc ...
!> \param mnhc ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE set_template_restart(work_section, eta, veta, fnhc, mnhc, error)
    TYPE(section_vals_type), POINTER         :: work_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eta, veta, fnhc, mnhc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_template_restart', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: coord, force, mass, velocity

    failure = .FALSE.
    NULLIFY (coord, force, velocity, mass)
    IF (PRESENT(eta)) THEN
       IF (SIZE(eta) > 0) THEN
          coord => section_vals_get_subs_vals(work_section,"COORD",error=error)
          CALL section_vals_val_set(coord,"_DEFAULT_KEYWORD_",r_vals_ptr=eta,error=error)
       ELSE
          DEALLOCATE (eta,STAT=stat)
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eta")
       END IF
    END IF
    IF (PRESENT(veta)) THEN
       IF (SIZE(veta) > 0) THEN
          velocity => section_vals_get_subs_vals(work_section,"VELOCITY",error=error)
          CALL section_vals_val_set(velocity,"_DEFAULT_KEYWORD_",r_vals_ptr=veta,error=error)
       ELSE
          DEALLOCATE (veta,STAT=stat)
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"veta")
       END IF
    END IF
    IF (PRESENT(fnhc)) THEN
       IF (SIZE(fnhc) > 0) THEN
          force => section_vals_get_subs_vals(work_section,"FORCE",error=error)
          CALL section_vals_val_set(force,"_DEFAULT_KEYWORD_",r_vals_ptr=fnhc,error=error)
       ELSE
          DEALLOCATE (fnhc,STAT=stat)
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"fnhc")
       END IF
    END IF
    IF (PRESENT(mnhc)) THEN
       IF (SIZE(mnhc) > 0) THEN
          mass => section_vals_get_subs_vals(work_section,"MASS",error=error)
          CALL section_vals_val_set(mass,"_DEFAULT_KEYWORD_",r_vals_ptr=mnhc,error=error)
       ELSE
          DEALLOCATE (mnhc,STAT=stat)
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mnhc")
       END IF
    END IF

  END SUBROUTINE set_template_restart

! *****************************************************************************
!> \brief Updates the force_eval section of the input file
!> \param force_env ...
!> \param logger ...
!> \param root_section ...
!> \param write_binary_restart_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_force_eval(force_env,logger,root_section,&
                               write_binary_restart_file,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: root_section
    LOGICAL, INTENT(IN)                      :: write_binary_restart_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_force_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iforce_eval, myid, nforce_eval
    INTEGER, DIMENSION(:), POINTER           :: i_force_eval
    LOGICAL                                  :: failure, multiple_subsys, &
                                                skip_vel_section
    TYPE(section_vals_type), POINTER         :: force_env_sections, &
                                                qmmm_section, rng_section, &
                                                subsys_section
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(rng_section, subsys_section, virial)
    ! If it's not a dynamical run don't update the velocity section
    CALL section_vals_val_get(root_section,"GLOBAL%RUN_TYPE",i_val=myid,error=error)
    skip_vel_section =((myid /= mol_dyn_run).AND.&
                       (myid /= mon_car_run).AND.&
                       (myid /= pint_run).AND.&
                       (myid /= ehrenfest))

    ! Go on updatig the force_env_section
    force_env_sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
    CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, nforce_eval, error)
    ! The update of the input MUST be realized only on the main force_eval
    ! All the others will be left not updated because there is no real need to update them...
    iforce_eval=1
    subsys_section => section_vals_get_subs_vals3(force_env_sections,"SUBSYS",&
         i_rep_section=i_force_eval(iforce_eval),error=error)
    CALL update_subsys(subsys_section,force_env,skip_vel_section,&
                       write_binary_restart_file,error)

    rng_section    => section_vals_get_subs_vals(subsys_section,"RNG_INIT",error=error)
    CALL update_rng_particle(rng_section,force_env,error)

    qmmm_section   => section_vals_get_subs_vals3(force_env_sections,"QMMM",&
         i_rep_section=i_force_eval(iforce_eval),error=error)
    CALL update_qmmm(qmmm_section,force_env,error)

    ! And now update only the cells of all other force_eval
    ! This is to make consistent for cell variable runs..
    IF (nforce_eval>1) THEN
       CALL force_env_get(force_env=force_env,&
                          virial=virial,&
                          error=error)
       CALL section_vals_val_get(root_section,"MULTIPLE_FORCE_EVALS%MULTIPLE_SUBSYS",l_val=multiple_subsys,error=error)
       IF (virial%pv_availability.AND.multiple_subsys) THEN
          DO iforce_eval = 2, nforce_eval
             subsys_section => section_vals_get_subs_vals3(force_env_sections,"SUBSYS",&
                  i_rep_section=i_force_eval(iforce_eval),error=error)
             CALL update_cell_section(force_env, subsys_section, error)
          END DO
       END IF
    END IF

    IF(myid == ehrenfest)CALL section_vals_val_set(root_section,"FORCE_EVAL%DFT%REAL_TIME_PROPAGATION%INITIAL_WFN",&
                                                   i_val=use_rt_restart,error=error)
    DEALLOCATE(i_force_eval)

  END SUBROUTINE update_force_eval

! *****************************************************************************
!> \brief Updates the qmmm section if needed
!> \param qmmm_section ...
!> \param force_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2007 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_qmmm(qmmm_section,force_env,error)
    TYPE(section_vals_type), POINTER         :: qmmm_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_qmmm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work

    failure = .FALSE.
    SELECT CASE(force_env%in_use)
    CASE(use_qmmm)
       CALL section_vals_get(qmmm_section, explicit=explicit, error=error)
       CPPostcondition(explicit,cp_failure_level,routineP,error,failure)

       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       work = force_env%qmmm_env%transl_v
       CALL section_vals_val_set(qmmm_section,"INITIAL_TRANSLATION_VECTOR",r_vals_ptr=work,error=error)
    END SELECT

  END SUBROUTINE update_qmmm

! *****************************************************************************
!> \brief Updates the rng section of the input file
!>      Write current status of the parallel random number generator (RNG)
!> \param rng_section ...
!> \param force_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_rng_particle(rng_section,force_env,error)

    TYPE(section_vals_type), POINTER         :: rng_section
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_rng_particle', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, nparticle, &
      nparticle_kind, nparticle_local
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: ascii
    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_list_type), POINTER        :: particles

    CALL force_env_get(force_env, subsys=subsys, para_env=para_env, error=error)
    CALL cp_subsys_get(subsys, atomic_kinds=atomic_kinds, local_particles=local_particles,&
         particles=particles, error=error)

    IF (ASSOCIATED(local_particles%local_particle_set)) THEN
       nparticle_kind = atomic_kinds%n_els
       nparticle = particles%n_els

       ALLOCATE (ascii(rng_record_length,nparticle),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ascii = 0

       DO iparticle=1,nparticle
          DO iparticle_kind=1,nparticle_kind
             nparticle_local = local_particles%n_el(iparticle_kind)
             DO iparticle_local=1,nparticle_local
                IF (iparticle == local_particles%list(iparticle_kind)%array(iparticle_local)) THEN
                   CALL dump_rng_stream(rng_stream=local_particles%local_particle_set(iparticle_kind)%&
                         rng(iparticle_local)%stream,&
                        rng_record=rng_record, error=error)
                   CALL string_to_ascii(rng_record,ascii(:,iparticle))
                END IF
             END DO
          END DO
       END DO

       CALL mp_sum(ascii,para_env%group)

       CALL section_rng_val_set(rng_section=rng_section,nsize=nparticle,ascii=ascii,&
            error=error)

       DEALLOCATE (ascii,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE update_rng_particle

! *****************************************************************************
!> \brief Updates the subsys section of the input file
!> \param subsys_section ...
!> \param force_env ...
!> \param skip_vel_section ...
!> \param write_binary_restart_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_subsys(subsys_section,force_env,skip_vel_section,&
                           write_binary_restart_file,error)
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN)                      :: skip_vel_section, &
                                                write_binary_restart_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_subsys', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: coord_file_name, unit_str
    INTEGER                                  :: coord_file_format, handle, &
                                                output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL                                  :: failure, scale
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    NULLIFY (work_section, core_particles, particles, shell_particles, &
             subsys, cell, para_env, multipoles)
    CALL timeset(routineN,handle)
    CALL force_env_get(force_env, subsys=subsys, cell=cell, para_env=para_env, &
                       error=error)

    CALL cp_subsys_get(subsys, particles=particles, molecules_new=molecules,&
         shell_particles=shell_particles, core_particles=core_particles,&
         multipoles=multipoles, error=error)

    ! Remove the multiple_unit_cell from the input structure.. at this point we have
    ! already all the information we need..
    ALLOCATE(multiple_unit_cell(3),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    multiple_unit_cell = 1
    CALL section_vals_val_set(subsys_section, "TOPOLOGY%MULTIPLE_UNIT_CELL",&
         i_vals_ptr=multiple_unit_cell, error=error)

    ! Coordinates and Velocities
    work_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
    CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,ignore_required=.TRUE.,error=error)
    CALL section_vals_val_get(work_section,"SCALED",l_val=scale,ignore_required=.TRUE.,error=error)
    conv_factor = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
    IF (.NOT.write_binary_restart_file) THEN
       CALL section_vals_val_get(subsys_section,"TOPOLOGY%COORD_FILE_FORMAT",&
                                 i_val=coord_file_format,&
                                 error=error)
       IF (coord_file_format == do_coord_cp2k) THEN
          CALL section_vals_val_get(subsys_section,"TOPOLOGY%COORD_FILE_NAME",&
                                    c_val=coord_file_name,&
                                    error=error)
          output_unit = 0
          IF (para_env%ionode) THEN
             CALL open_file(file_name=TRIM(ADJUSTL(coord_file_name)),&
                            file_action="READWRITE",&
                            file_form="FORMATTED",&
                            file_position="REWIND",&
                            file_status="UNKNOWN",&
                            unit_number=output_unit)
            CALL dump_coordinates_cp2k(particles,molecules,cell,conv_factor,&
                                       output_unit=output_unit,&
                                       core_or_shell=.FALSE.,&
                                       scaled_coordinates=scale,&
                                       error=error)
            CALL close_file(unit_number=output_unit)
          END IF
       ELSE
          CALL section_coord_val_set(work_section, particles, molecules, conv_factor, scale,&
                                     cell, error=error)
       END IF
    END IF
    CALL section_vals_val_set(subsys_section,"TOPOLOGY%NUMBER_OF_ATOMS",i_val=particles%n_els,error=error)
    work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
    IF (.NOT.skip_vel_section) THEN
       IF (.NOT.write_binary_restart_file) THEN
          CALL section_velocity_val_set(work_section, particles, conv_factor=1.0_dp, error=error)
       END IF
    ELSE
       CALL section_vals_remove_values(work_section, error)
    END IF

    ! Write restart input for core-shell model
    IF (.NOT.write_binary_restart_file) THEN
       IF (ASSOCIATED(shell_particles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"SHELL_COORD",error=error)
          CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,ignore_required=.TRUE.,error=error)
          CALL section_vals_val_get(work_section,"SCALED",l_val=scale,ignore_required=.TRUE.,error=error)
          conv_factor = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          CALL section_coord_val_set(work_section, shell_particles, molecules, &
                                     conv_factor, scale, cell, shell=.TRUE., error=error)
          IF (.NOT.skip_vel_section) THEN
             work_section => section_vals_get_subs_vals(subsys_section,"SHELL_VELOCITY",error=error)
             CALL section_velocity_val_set(work_section, shell_particles, conv_factor=1.0_dp, error=error)
          END IF
       END IF
       IF (ASSOCIATED(core_particles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"CORE_COORD",error=error)
          CALL section_vals_val_get(work_section,"UNIT",c_val=unit_str,ignore_required=.TRUE.,error=error)
          CALL section_vals_val_get(work_section,"SCALED",l_val=scale,ignore_required=.TRUE.,error=error)
          conv_factor = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)
          CALL section_coord_val_set(work_section, core_particles, molecules, &
                                     conv_factor, scale, cell, shell=.TRUE., error=error)
          IF (.NOT.skip_vel_section) THEN
             work_section => section_vals_get_subs_vals(subsys_section,"CORE_VELOCITY",error=error)
             CALL section_velocity_val_set(work_section, core_particles, conv_factor=1.0_dp, error=error)
          END IF
       END IF
    END IF

    ! Updating cell info
    CALL update_cell_section(force_env, subsys_section, error)
    ! Updating multipoles
    IF (ASSOCIATED(multipoles)) THEN
       work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES",error=error)
       DO
          IF (SIZE(work_section%values,2)==1) EXIT
          CALL section_vals_add_values(work_section,error=error)
       END DO
       IF (ASSOCIATED(multipoles%dipoles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%DIPOLES",error=error)
          CALL update_dipoles_section(multipoles%dipoles, work_section, error)
       END IF
       IF (ASSOCIATED(multipoles%quadrupoles)) THEN
          work_section => section_vals_get_subs_vals(subsys_section,"MULTIPOLES%QUADRUPOLES",error=error)
          CALL update_quadrupoles_section(multipoles%quadrupoles, work_section, error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE update_subsys

! *****************************************************************************
!> \brief routine to update cell section.. splitted from the main update subsys_section
!> \param force_env ...
!> \param subsys_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2007 created [teo] Zurich University
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_cell_section(force_env, subsys_section, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_cell_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:), POINTER           :: iwork
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    NULLIFY(work, work_section, cell, cell_ref, iwork)
    CALL timeset(routineN,handle)
    ! Handle case for parallel runs
    IF (ASSOCIATED(force_env)) THEN
       CALL force_env_get(force_env,&
                          cell=cell,&
                          cell_ref=cell_ref,&
                          use_ref_cell=use_ref_cell,&
                          error=error)
    END IF

    ! CELL
    work_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
    ALLOCATE(work(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CELL VECTORS - A
    work(1:3) = cell%hmat(1:3,1)
    CALL section_vals_val_set(work_section,"A",r_vals_ptr=work,error=error)
    ALLOCATE(work(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CELL VECTORS - B
    work(1:3) = cell%hmat(1:3,2)
    CALL section_vals_val_set(work_section,"B",r_vals_ptr=work,error=error)
    ALLOCATE(work(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ! CELL VECTORS - C
    work(1:3) = cell%hmat(1:3,3)
    CALL section_vals_val_set(work_section,"C",r_vals_ptr=work,error=error)
    ! MULTIPLE_UNIT_CELL
    ALLOCATE(iwork(3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    iwork = 1
    CALL section_vals_val_set(work_section,"MULTIPLE_UNIT_CELL",i_vals_ptr=iwork,error=error)
    ! Unset unused or misleading information
    CALL section_vals_val_unset(work_section,"ABC",error=error)
    CALL section_vals_val_unset(work_section,"ALPHA_BETA_GAMMA",error=error)
    ! CELL_REF
    IF (use_ref_cell.AND.ASSOCIATED(cell_ref)) THEN
       work_section => section_vals_get_subs_vals(subsys_section,"CELL%CELL_REF",error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! CELL VECTORS - A
       work(1:3) = cell_ref%hmat(1:3,1)
       CALL section_vals_val_set(work_section,"A",r_vals_ptr=work,error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! CELL VECTORS - B
       work(1:3) = cell_ref%hmat(1:3,2)
       CALL section_vals_val_set(work_section,"B",r_vals_ptr=work,error=error)
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! CELL VECTORS - C
       work(1:3) = cell_ref%hmat(1:3,3)
       CALL section_vals_val_set(work_section,"C",r_vals_ptr=work,error=error)
       ! MULTIPLE_UNIT_CELL
       ALLOCATE(iwork(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       iwork = 1
       CALL section_vals_val_set(work_section,"MULTIPLE_UNIT_CELL",i_vals_ptr=iwork,error=error)
       ! Unset unused or misleading information
       CALL section_vals_val_unset(work_section,"ABC",error=error)
       CALL section_vals_val_unset(work_section,"ALPHA_BETA_GAMMA",error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE update_cell_section

! *****************************************************************************
!> \brief routine to dump coordinates.. fast implementation
!> \param coord_section ...
!> \param particles ...
!> \param molecules ...
!> \param conv_factor ...
!> \param scale ...
!> \param cell ...
!> \param shell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_coord_val_set(coord_section, particles, molecules, conv_factor,&
       scale, cell, shell, error)
    TYPE(section_vals_type), POINTER         :: coord_section
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_new_list_type), POINTER         :: molecules
    REAL(KIND=dp)                            :: conv_factor
    LOGICAL, INTENT(IN)                      :: scale
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN), OPTIONAL            :: shell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'section_coord_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: line
    CHARACTER(LEN=default_string_length)     :: my_tag, name
    INTEGER                                  :: handle, ik, imol, irk, &
                                                last_atom, Nlist
    LOGICAL                                  :: failure, ldum, &
                                                molname_generated, my_shell
    REAL(KIND=dp), DIMENSION(3)              :: r, s
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(molecule_type), POINTER             :: molecule_now
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    my_shell = .FALSE.
    IF(PRESENT(shell)) my_shell = shell
    CPPrecondition(ASSOCIATED(coord_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(coord_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => coord_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    DO
       IF (SIZE(coord_section%values,2)==1) EXIT
       CALL section_vals_add_values(coord_section,error=error)
    END DO
    vals  => coord_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    imol=0
    last_atom=0
    DO irk=1,particles%n_els
       CALL get_atomic_kind(particles%els(irk)%atomic_kind,name=name)
       IF (my_shell) THEN
         s = particles%els(irk)%r
         IF (scale) THEN
            r = s
            CALL real_to_scaled(s,r,cell)
         ELSE
            s = s*conv_factor
         END IF
         WRITE (UNIT=line,FMT="(A,3(1X,ES25.16),1X,I0)")&
           TRIM(name),s(1:3),particles%els(irk)%atom_index
         CALL val_create(my_val,lc_val=line,error=error)
         IF (Nlist /= 0) THEN
            IF (irk==1) THEN
               new_pos => vals
            ELSE
               new_pos => new_pos%rest
            END IF
            old_val => new_pos%first_el
            CALL val_release(old_val,error=error)
            new_pos%first_el => my_val
         ELSE
            IF (irk==1) THEN
               NULLIFY(new_pos)
               CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
               vals => new_pos
            ELSE
               NULLIFY(new_pos%rest)
               CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
               new_pos => new_pos%rest
            END IF
         END IF
         NULLIFY(my_val)
       ELSE
         IF(last_atom<irk) THEN
           imol=imol+1
           molecule_now => molecules%els(imol)
           CALL get_molecule(molecule_now,last_atom=last_atom)
           CALL get_molecule_kind(molecule_now%molecule_kind,molname_generated=molname_generated,&
                name=my_tag)
           IF (molname_generated) my_tag=""
         END IF
         ldum = qmmm_ff_precond_only_qm(my_tag)
         ldum = qmmm_ff_precond_only_qm(name)
         s = particles%els(irk)%r
         IF (scale) THEN
            r = s
            CALL real_to_scaled(s,r,cell)
         ELSE
            s = s*conv_factor
         END IF
         IF (LEN_TRIM(my_tag) > 0) THEN
            WRITE (UNIT=line,FMT="(A,3(1X,ES25.16),1X,A,1X,I0)")&
              TRIM(name),s(1:3),TRIM(my_tag),imol
         ELSE
            WRITE (UNIT=line,FMT="(A,3(1X,ES25.16))")&
              TRIM(name),s(1:3)
         END IF
         CALL val_create(my_val,lc_val=line,error=error)

         IF (Nlist /= 0) THEN
            IF (irk==1) THEN
               new_pos => vals
            ELSE
               new_pos => new_pos%rest
            END IF
            old_val => new_pos%first_el
            CALL val_release(old_val,error=error)
            new_pos%first_el => my_val
         ELSE
            IF (irk==1) THEN
               NULLIFY(new_pos)
               CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
               vals => new_pos
            ELSE
               NULLIFY(new_pos%rest)
               CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
               new_pos => new_pos%rest
            END IF
         END IF
         NULLIFY(my_val)
       END IF
    END DO
    coord_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE section_coord_val_set

! *****************************************************************************
!> \brief routine to dump dipoles.. fast implementation
!> \param dipoles ...
!> \param dipoles_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_dipoles_section(dipoles, dipoles_section, error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dipoles
    TYPE(section_vals_type), POINTER         :: dipoles_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_dipoles_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ik, irk, Nlist, &
                                                nloop, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(dipoles_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(dipoles_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => dipoles_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    ! At least one of the two arguments must be present..
    nloop = SIZE(dipoles,2)
    DO
       IF (SIZE(dipoles_section%values,2)==1) EXIT
       CALL section_vals_add_values(dipoles_section,error=error)
    END DO
    vals  => dipoles_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,nloop
       ALLOCATE(work(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       work = dipoles(1:3,irk)
       CALL val_create(my_val,r_vals_ptr=work,error=error)

       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    dipoles_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE update_dipoles_section

! *****************************************************************************
!> \brief routine to dump quadrupoles.. fast implementation
!> \param quadrupoles ...
!> \param quadrupoles_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2007 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_quadrupoles_section(quadrupoles, quadrupoles_section, error)

    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: quadrupoles
    TYPE(section_vals_type), POINTER         :: quadrupoles_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_quadrupoles_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ik, ind, irk, j, &
                                                Nlist, nloop, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(quadrupoles_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(quadrupoles_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => quadrupoles_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    ! At least one of the two arguments must be present..
    nloop = SIZE(quadrupoles,2)
    DO
       IF (SIZE(quadrupoles_section%values,2)==1) EXIT
       CALL section_vals_add_values(quadrupoles_section,error=error)
    END DO
    vals  => quadrupoles_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,nloop
       ALLOCATE(work(6),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       ind = 0
       DO i = 1, 3
          DO j = i, 3
             ind  = ind + 1
             work(ind) = quadrupoles(j,i,irk)
          END DO
       END DO
       CALL val_create(my_val,r_vals_ptr=work,error=error)

       IF (Nlist /= 0) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    quadrupoles_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE update_quadrupoles_section

! *****************************************************************************
!> \brief routine to dump rngs.. fast implementation
!> \param rng_section ...
!> \param nsize ...
!> \param ascii ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!>      - string dump (again) instead of integer ASCII code (07.03.06,MK)
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_rng_val_set(rng_section, nsize, ascii, error)

    TYPE(section_vals_type), POINTER         :: rng_section
    INTEGER, INTENT(IN)                      :: nsize
    INTEGER, DIMENSION(:, :)                 :: ascii
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'section_rng_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=rng_record_length)         :: rng_record
    INTEGER                                  :: ik, irk, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(rng_section),cp_failure_level,routineP,error,failure)
    CPPrecondition((rng_section%ref_count > 0),cp_failure_level,routineP,error,failure)

    NULLIFY (my_val,old_val,section,vals)

    section => rng_section%section

    ik = section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)

    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
                   "section "//TRIM(section%name)//" does not contain keyword "//&
                   "_DEFAULT_KEYWORD_",error,failure)

    DO
      IF (SIZE(rng_section%values,2)==1) EXIT
      CALL section_vals_add_values(rng_section,error=error)
    END DO

    vals => rng_section%values(ik,1)%list
    Nlist = 0

    IF (ASSOCIATED(vals)) THEN
      Nlist = cp_sll_val_get_length(vals,error)
    END IF

    DO irk=1,nsize

      CALL ascii_to_string(ascii(:,irk),rng_record)
      CALL val_create(val=my_val,lc_val=rng_record,error=error)

      IF (Nlist /= 0) THEN
        IF (irk == 1) THEN
          new_pos => vals
        ELSE
          new_pos => new_pos%rest
        END IF
        old_val => new_pos%first_el
        CALL val_release(old_val,error=error)
        new_pos%first_el => my_val
      ELSE
        IF (irk == 1) THEN
          NULLIFY (new_pos)
          CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
          vals => new_pos
        ELSE
          NULLIFY (new_pos%rest)
          CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
          new_pos => new_pos%rest
        END IF
      END IF

      NULLIFY (my_val)

    END DO

    rng_section%values(ik,1)%list => vals

  END SUBROUTINE section_rng_val_set

! *****************************************************************************
!> \brief routine to dump hills information during metadynamics run
!> \param ss_section ...
!> \param meta_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE meta_hills_val_set_ss(ss_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: ss_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'meta_hills_val_set_ss', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ss_val
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(ss_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ss_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => ss_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(ss_section%values,2)==1) EXIT
       CALL section_vals_add_values(ss_section,error=error)
    END DO
    vals  => ss_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = SIZE(meta_env%hills_env%ss_history,1)
    DO irk=1,meta_env%hills_env%n_hills
       ALLOCATE(ss_val(lsize),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       ss_val = meta_env%hills_env%ss_history(:,irk)
       CALL val_create(my_val,r_vals_ptr=ss_val,error=error)

       IF (irk <= Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    ss_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_ss

! *****************************************************************************
!> \brief routine to dump hills information during metadynamics run
!> \param ds_section ...
!> \param meta_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE meta_hills_val_set_ds(ds_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: ds_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'meta_hills_val_set_ds', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ds_val
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(ds_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ds_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => ds_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(ds_section%values,2)==1) EXIT
       CALL section_vals_add_values(ds_section,error=error)
    END DO
    vals  => ds_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = SIZE(meta_env%hills_env%delta_s_history,1)
    DO irk=1,meta_env%hills_env%n_hills
       ALLOCATE(ds_val(lsize),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       ds_val = meta_env%hills_env%delta_s_history(:,irk)
       CALL val_create(my_val,r_vals_ptr=ds_val,error=error)

       IF (irk<=Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    ds_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_ds

! *****************************************************************************
!> \brief routine to dump hills information during metadynamics run
!> \param ww_section ...
!> \param meta_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE meta_hills_val_set_ww(ww_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: ww_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'meta_hills_val_set_ww', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(ww_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ww_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => ww_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(ww_section%values,2)==1) EXIT
       CALL section_vals_add_values(ww_section,error=error)
    END DO
    vals  => ww_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = meta_env%hills_env%n_hills
    DO irk=1,lsize
       CALL val_create(my_val,r_val=meta_env%hills_env%ww_history(irk),error=error)

       IF (irk<=Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    ww_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_ww

! *****************************************************************************
!> \brief routine to dump hills information during metadynamics run
!> \param invdt_section ...
!> \param meta_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2009 created [seb]
!> \author SC
! *****************************************************************************
  SUBROUTINE meta_hills_val_set_dt(invdt_section, meta_env, error)

    TYPE(section_vals_type), POINTER         :: invdt_section
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'meta_hills_val_set_dt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, lsize, Nlist
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(invdt_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(invdt_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => invdt_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)
    DO
       IF (SIZE(invdt_section%values,2)==1) EXIT
       CALL section_vals_add_values(invdt_section,error=error)
    END DO
    vals  => invdt_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    lsize = meta_env%hills_env%n_hills
    DO irk=1,lsize
       CALL val_create(my_val,r_val=meta_env%hills_env%invdt_history(irk),error=error)

       IF (irk<=Nlist) THEN
          IF (irk==1) THEN
             new_pos => vals
          ELSE
             new_pos => new_pos%rest
          END IF
          old_val => new_pos%first_el
          CALL val_release(old_val,error=error)
          new_pos%first_el => my_val
       ELSE
          IF (irk==1) THEN
             NULLIFY(new_pos)
             CALL cp_sll_val_create(new_pos,first_el=my_val,error=error)
             vals => new_pos
          ELSE
             NULLIFY(new_pos%rest)
             CALL cp_sll_val_create(new_pos%rest,first_el=my_val,error=error)
             new_pos => new_pos%rest
          END IF
       END IF
       NULLIFY(my_val)
    END DO
    invdt_section%values(ik,1)%list => vals
  END SUBROUTINE meta_hills_val_set_dt

! *****************************************************************************
!> \brief   Dump atomic, core, or shell coordinates to a file in CP2K &COORD
!>          section format
!> \param particles ...
!> \param molecules ...
!> \param cell ...
!> \param conv_factor ...
!> \param output_unit ...
!> \param core_or_shell ...
!> \param scaled_coordinates ...
!> \param error ...
!> \par History
!>      07.02.2011 (Creation, MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dump_coordinates_cp2k(particles,molecules,cell,conv_factor,&
                                   output_unit,core_or_shell,&
                                   scaled_coordinates,error)

    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(IN)                :: conv_factor
    INTEGER, INTENT(IN)                      :: output_unit
    LOGICAL, INTENT(IN)                      :: core_or_shell, &
                                                scaled_coordinates
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dump_coordinates_cp2k', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: kind_name, tag
    INTEGER                                  :: handle, imolecule, iparticle, &
                                                last_atom
    LOGICAL                                  :: failure, ldum, &
                                                molname_generated
    REAL(KIND=dp), DIMENSION(3)              :: r, s
    TYPE(molecule_type), POINTER             :: molecule

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(particles),cp_failure_level,routineP,error,failure)
    IF (.NOT.core_or_shell) THEN
       CPPrecondition(ASSOCIATED(molecules),cp_failure_level,routineP,error,failure)
    END IF
    IF (scaled_coordinates) THEN
       CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
    END IF

    kind_name = ""
    tag = ""
    imolecule = 0
    last_atom = 0
    DO iparticle=1,particles%n_els
       CALL get_atomic_kind(particles%els(iparticle)%atomic_kind,name=kind_name)
       IF (.NOT.core_or_shell) THEN
          IF (iparticle > last_atom) THEN
             imolecule = imolecule + 1
             molecule => molecules%els(imolecule)
             CALL get_molecule(molecule,last_atom=last_atom)
             CALL get_molecule_kind(molecule%molecule_kind,&
                                    molname_generated=molname_generated,&
                                    name=tag)
             IF (molname_generated) tag = ""
          END IF
          ldum = qmmm_ff_precond_only_qm(tag)
          ldum = qmmm_ff_precond_only_qm(kind_name)
       END IF
       IF (scaled_coordinates) THEN
          CALL real_to_scaled(s,particles%els(iparticle)%r,cell)
          r(1:3) = s(1:3)
       ELSE
          r(1:3) = particles%els(iparticle)%r(1:3)*conv_factor
       END IF
       IF (core_or_shell) THEN
          WRITE (UNIT=output_unit,FMT="(A,3(1X,ES25.16),1X,I0)")&
            TRIM(ADJUSTL(kind_name)),r(1:3),particles%els(iparticle)%atom_index
       ELSE
          IF (LEN_TRIM(tag) > 0) THEN
             WRITE (UNIT=output_unit,FMT="(A,3(1X,ES25.16),1X,A,1X,I0)")&
               TRIM(ADJUSTL(kind_name)),r(1:3),TRIM(tag),imolecule
          ELSE
             WRITE (UNIT=output_unit,FMT="(A,3(1X,ES25.16))")&
               TRIM(ADJUSTL(kind_name)),r(1:3)
          END IF
       END IF
    END DO

    CALL timestop(handle)

  END SUBROUTINE dump_coordinates_cp2k

! *****************************************************************************
!> \brief   Write all input sections scaling in size with the number of atoms
!>          in the system to an external file in binary format
!> \param output_unit ...
!> \param root_section ...
!> \param md_env ...
!> \param force_env ...
!> \param error ...
!> \par History
!>      - Creation (10.02.2011,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_binary_restart(output_unit,root_section,md_env,force_env,&
                                  error)

    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_binary_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: binary_restart_file_name
    CHARACTER(LEN=default_string_length)     :: section_label
    INTEGER :: handle, iatom, icore, ikind, imolecule, ishell, istat, &
      n_char_size, n_dp_size, n_int_size, natom, natomkind, ncore, nhc_size, &
      nmolecule, nmoleculekind, nshell, print_level, run_type
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ibuf, imol
    LOGICAL                                  :: print_info, write_velocities
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: my_force_env
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(thermostat_type), POINTER           :: thermostat_part, &
                                                thermostat_shell

    CALL timeset(routineN,handle)

    NULLIFY (atomic_kinds)
    NULLIFY (core_particles)
    NULLIFY (molecule_kinds)
    NULLIFY (molecules)
    NULLIFY (my_force_env)
    NULLIFY (para_env)
    NULLIFY (particles)
    NULLIFY (shell_particles)
    NULLIFY (subsys)
    NULLIFY (thermostat_part)
    NULLIFY (thermostat_shell)

    IF (PRESENT(md_env)) THEN
       CALL get_md_env(md_env=md_env,&
                       force_env=my_force_env,&
                       thermostat_part=thermostat_part,&
                       thermostat_shell=thermostat_shell,&
                       error=error)
    ELSE IF (PRESENT(force_env)) THEN
       my_force_env => force_env
    END IF

    IF (.NOT.ASSOCIATED(my_force_env)) THEN
       CALL timestop(handle)
       RETURN
    END IF

    CALL section_vals_val_get(root_section,"GLOBAL%PRINT_LEVEL",i_val=print_level,error=error)
    IF (print_level > 1) THEN
       print_info = .TRUE.
    ELSE
       print_info = .FALSE.
    END IF

    CALL section_vals_val_get(root_section,"GLOBAL%RUN_TYPE",i_val=run_type,error=error)
    write_velocities = ((run_type == mol_dyn_run).OR.&
                        (run_type == mon_car_run).OR.&
                        (run_type == pint_run))

    CALL force_env_get(force_env=my_force_env,&
                       para_env=para_env,&
                       subsys=subsys,&
                       error=error)
    CALL cp_subsys_get(subsys,&
                       atomic_kinds=atomic_kinds,&
                       particles=particles,&
                       natom=natom,&
                       core_particles=core_particles,&
                       ncore=ncore,&
                       shell_particles=shell_particles,&
                       nshell=nshell,&
                       molecule_kinds_new=molecule_kinds,&
                       molecules_new=molecules,&
                       error=error)

    natomkind = atomic_kinds%n_els
    IF (ASSOCIATED(molecule_kinds)) THEN
       nmoleculekind = molecule_kinds%n_els
    ELSE
       nmoleculekind = 0
    END IF

    IF (ASSOCIATED(molecules)) THEN
       nmolecule = molecules%n_els
    ELSE
       nmolecule = 0
    END IF

    n_char_size = 0 ! init
    n_int_size = 0 ! init
    n_dp_size = 0 ! init

    IF (output_unit > 0) THEN ! only ionode

       IF (print_info) THEN
          INQUIRE (UNIT=output_unit,NAME=binary_restart_file_name,IOSTAT=istat)
          IF (istat /= 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "An error occurred inquiring logical unit <"//&
                               TRIM(ADJUSTL(cp_to_string(output_unit)))//&
                               "> which should be linked to the binary restart file",&
                               para_env)
          END IF
          WRITE (UNIT=*,FMT="(T2,A,/,/,(T3,A,T71,I10))")&
            "Writing binary restart file "//TRIM(ADJUSTL(binary_restart_file_name)),&
            "Number of atomic kinds:",natomkind,&
            "Number of atoms:",natom,&
            "Number of cores (only core-shell model):",ncore,&
            "Number of shells (only core-shell model):",nshell,&
            "Number of molecule kinds:",nmoleculekind,&
            "Number of molecules",nmolecule

          n_int_size = n_int_size + 6
       END IF

       WRITE (UNIT=output_unit,IOSTAT=istat)&
         natomkind,natom,ncore,nshell,nmoleculekind,nmolecule
       IF (istat /= 0) THEN
          CALL stop_write(routineN,moduleN,__LINE__,&
                          "natomkind,natom,ncore,nshell,nmoleculekind,nmolecule "//&
                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                          output_unit,para_env)
       END IF

       ! Write atomic kind names
       DO ikind=1,natomkind
          WRITE (UNIT=output_unit,IOSTAT=istat) atomic_kinds%els(ikind)%name
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "atomic_kinds%els(ikind)%name "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_char_size = n_char_size + LEN(atomic_kinds%els(ikind)%name)
       END DO

       ! Write atomic kind numbers of all atoms
       ALLOCATE (ibuf(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "ibuf",int_size*natom,para_env)
       DO iatom=1,natom
          ibuf(iatom) = particles%els(iatom)%atomic_kind%kind_number
       END DO
       WRITE (UNIT=output_unit,IOSTAT=istat) ibuf(1:natom)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "ibuf(1:natom) -> atomic kind numbers "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_int_size = n_int_size + natom
       ! Write atomic coordinates
       ALLOCATE (rbuf(3,natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "rbuf",dp_size*3*natom,para_env)
       DO iatom=1,natom
          rbuf(1:3,iatom) = particles%els(iatom)%r(1:3)
       END DO
       WRITE (UNIT=output_unit,IOSTAT=istat) rbuf(1:3,1:natom)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "rbuf(1:3,1:natom) -> atomic coordinates "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_dp_size = n_dp_size + 3*natom
       DEALLOCATE (rbuf,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)

       ! Write molecule information if available
       IF (nmolecule > 0) THEN
          ! Write molecule kind names
          DO ikind=1,nmoleculekind
             WRITE (UNIT=output_unit,IOSTAT=istat) molecule_kinds%els(ikind)%name
             IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                             "molecule_kinds%els(ikind)%name "//&
                                             "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                             output_unit,para_env)
             n_char_size = n_char_size + LEN(molecule_kinds%els(ikind)%name)
          END DO
          ! Write molecule (kind) index numbers for all atoms
          ibuf(:) = 0
          ALLOCATE (imol(natom),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "imol",int_size*natom,para_env)
          imol(:) = 0
          DO imolecule=1,nmolecule
             ikind = molecules%els(imolecule)%molecule_kind%kind_number
             DO iatom=molecules%els(imolecule)%first_atom,&
                      molecules%els(imolecule)%last_atom
                ibuf(iatom) = ikind
                imol(iatom) = imolecule
             END DO
          END DO
          ! Write molecule kind index number for each atom
          WRITE (UNIT=output_unit,IOSTAT=istat) ibuf(1:natom)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "ibuf(1:natom) -> molecule kind index numbers "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_int_size = n_int_size + natom
          ! Write molecule index number for each atom
          WRITE (UNIT=output_unit,IOSTAT=istat) imol(1:natom)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "imol(1:natom) -> molecule index numbers "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_int_size = n_int_size + natom
          DEALLOCATE (imol,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"imol",&
                                           para_env)
       END IF ! molecules

       DEALLOCATE (ibuf,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ibuf",para_env)

       ! Core-shell model only
       section_label = "SHELL COORDINATES"
       WRITE (UNIT=output_unit,IOSTAT=istat) section_label,nshell
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "section_label, nshell "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (nshell > 0) THEN
          ! Write shell coordinates
          ALLOCATE (rbuf(3,nshell),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "rbuf",dp_size*3*nshell,para_env)
          DO ishell=1,nshell
             rbuf(1:3,ishell) = shell_particles%els(ishell)%r(1:3)
          END DO
          WRITE (UNIT=output_unit,IOSTAT=istat) rbuf(1:3,1:nshell)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "rbuf(1:3,1:nshell) -> shell coordinates "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_dp_size = n_dp_size + 3*nshell
          DEALLOCATE (rbuf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
          ! Write atomic indices, i.e. number of the atom the shell belongs to
          ALLOCATE (ibuf(nshell),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "ibuf",int_size*nshell,para_env)
          DO ishell=1,nshell
             ibuf(ishell) = shell_particles%els(ishell)%atom_index
          END DO
          WRITE (UNIT=output_unit,IOSTAT=istat) ibuf(1:nshell)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "ibuf(1:nshell) -> atomic indices "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_int_size = n_int_size + nshell
          DEALLOCATE (ibuf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ibuf",para_env)
       END IF

       section_label = "CORE COORDINATES"
       WRITE (UNIT=output_unit,IOSTAT=istat) section_label,ncore
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "section_label, ncore "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (ncore > 0) THEN
          ! Write core coordinates
          ALLOCATE (rbuf(3,ncore),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "rbuf",dp_size*3*ncore,para_env)
          DO icore=1,ncore
             rbuf(1:3,icore) = core_particles%els(icore)%r(1:3)
          END DO
          WRITE (UNIT=output_unit,IOSTAT=istat) rbuf(1:3,1:ncore)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "rbuf(1:3,1:ncore) -> core coordinates "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_dp_size = n_dp_size + 3*ncore
          DEALLOCATE (rbuf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
          ! Write atomic indices, i.e. number of the atom the core belongs to
          ALLOCATE (ibuf(ncore),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "ibuf",int_size*ncore,para_env)
          DO icore=1,ncore
             ibuf(icore) = core_particles%els(icore)%atom_index
          END DO
          WRITE (UNIT=output_unit,IOSTAT=istat) ibuf(1:ncore)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "ibuf(1:ncore) -> atomic indices "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_int_size = n_int_size + ncore
          DEALLOCATE (ibuf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ibuf",para_env)
       END IF
    END IF ! ionode only

    ! Thermostat information

    ! Particle thermostats
    section_label = "PARTICLE THERMOSTATS"
    IF (ASSOCIATED(thermostat_part)) THEN
       ! Nose-Hoover thermostats
       IF (thermostat_part%type_of_thermostat == do_thermo_nose) THEN
          nhc => thermostat_part%nhc
          CALL write_binary_thermostats_nose(nhc,output_unit,section_label,&
                                             n_char_size,n_dp_size,n_int_size,&
                                             print_info,para_env,error)
       END IF
    ELSE
       nhc_size = 0
       IF (output_unit > 0) THEN
          WRITE (UNIT=output_unit,IOSTAT=istat) section_label,nhc_size
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          TRIM(section_label)//", nhc_size "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
       ENDIF
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (output_unit > 0) THEN ! only ionode
          IF (print_info) THEN
             WRITE (UNIT=*,FMT="(T3,A,T71,I10)")&
              "NHC size ("//TRIM(ADJUSTL(section_label))//")",nhc_size
          END IF
       END IF
    END IF

    ! Shell thermostats (only for core-shell models)
    section_label = "SHELL THERMOSTATS"
    IF (ASSOCIATED(thermostat_shell)) THEN
       ! Nose-Hoover thermostats
       IF (thermostat_shell%type_of_thermostat == do_thermo_nose) THEN
          nhc => thermostat_shell%nhc
          CALL write_binary_thermostats_nose(nhc,output_unit,section_label,&
                                             n_char_size,n_dp_size,n_int_size,&
                                             print_info,para_env,error)
       END IF
    ELSE
       nhc_size = 0
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,IOSTAT=istat) section_label,nhc_size
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "nhc_size "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
       ENDIF
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (output_unit > 0) THEN ! only ionode
          IF (print_info) THEN
             WRITE (UNIT=*,FMT="(T3,A,T71,I10)")&
               "NHC size ("//TRIM(ADJUSTL(section_label))//")",nhc_size
          END IF
       END IF
    END IF

    ! Particle velocities

    IF (output_unit > 0) THEN ! only ionode
       ! Write particle velocities if needed
       section_label = "VELOCITIES"
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,IOSTAT=istat) section_label,MERGE(natom,0,write_velocities)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          TRIM(section_label)//", write_velocities "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
       ENDIF
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (print_info) THEN
          WRITE (UNIT=*,FMT="(T3,A,T78,A3)")&
            "Write "//TRIM(ADJUSTL(section_label))//" section",MERGE("YES"," NO",write_velocities)
       END IF
       IF (write_velocities) THEN
          ALLOCATE (rbuf(3,natom),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "rbuf",dp_size*3*natom,para_env)
          ! Write atomic velocities
          DO iatom=1,natom
             rbuf(1:3,iatom) = particles%els(iatom)%v(1:3)
          END DO
          WRITE (UNIT=output_unit,IOSTAT=istat) rbuf(1:3,1:natom)
          IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                          "rbuf(1:3,1:natom) -> atomic velocities "//&
                                          "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                          output_unit,para_env)
          n_dp_size = n_dp_size + 3*natom
          DEALLOCATE (rbuf,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
       END IF
       ! Write shell velocities
       section_label = "SHELL VELOCITIES"
       WRITE (UNIT=output_unit,IOSTAT=istat) section_label,MERGE(nshell,0,write_velocities)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       TRIM(section_label)//", write_velocities "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (print_info) THEN
          WRITE (UNIT=*,FMT="(T3,A,T78,A3)")&
            "Write "//TRIM(ADJUSTL(section_label))//" section",MERGE("YES"," NO",write_velocities)
       END IF
       IF (nshell > 0) THEN
          IF (write_velocities) THEN
             ALLOCATE (rbuf(3,nshell),STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "rbuf",dp_size*3*nshell,para_env)
             DO ishell=1,nshell
                rbuf(1:3,ishell) = shell_particles%els(ishell)%v(1:3)
             END DO
             WRITE (UNIT=output_unit,IOSTAT=istat) rbuf(1:3,1:nshell)
             IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                             "rbuf(1:3,1:nshell) -> shell velocities "//&
                                             "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                             output_unit,para_env)
             n_dp_size = n_dp_size + 3*nshell
             DEALLOCATE (rbuf,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
          END IF
       END IF
       ! Write core velocities
       section_label = "CORE VELOCITIES"
       WRITE (UNIT=output_unit,IOSTAT=istat) section_label,MERGE(ncore,0,write_velocities)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       TRIM(section_label)//", write_velocities "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (print_info) THEN
          WRITE (UNIT=*,FMT="(T3,A,T78,A3)")&
            "Write "//TRIM(ADJUSTL(section_label))//" section",MERGE("YES"," NO",write_velocities)
       END IF
       IF (ncore > 0) THEN
          IF (write_velocities) THEN
             ALLOCATE (rbuf(3,ncore),STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                              "rbuf",dp_size*3*ncore,para_env)
             DO icore=1,ncore
                rbuf(1:3,icore) = core_particles%els(icore)%v(1:3)
             END DO
             WRITE (UNIT=output_unit,IOSTAT=istat) rbuf(1:3,1:ncore)
             IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                             "rbuf(1:3,1:ncore) -> core velocities "//&
                                             "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                             output_unit,para_env)
             n_dp_size = n_dp_size + 3*ncore
             DEALLOCATE (rbuf,STAT=istat)
             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
          END IF
       END IF
    END IF ! ionode only

    ! Optionally, print a small I/O statistics
    IF (output_unit > 0) THEN ! only ionode
       IF (print_info) THEN
          WRITE (UNIT=*,FMT="(/,(T2,I10,1X,I0,A,T68,I10,A))")&
            n_char_size,int_size,"-byte characters written",n_char_size*int_size/1024," KB",&
            n_dp_size,dp_size,"-byte floating point numbers written",n_dp_size*dp_size/1024," KB",&
            n_int_size,int_size,"-byte integer numbers written",n_int_size*int_size/1024," KB"
          WRITE (UNIT=*,FMT="(/,T2,A)")&
            "Binary restart file "//TRIM(ADJUSTL(binary_restart_file_name))//" written"
       END IF
    END IF ! ionode only

    CALL timestop(handle)

  END SUBROUTINE write_binary_restart

! *****************************************************************************
!> \brief   Write an input section for Nose thermostats to an external file in
!>          binary format
!> \param nhc ...
!> \param output_unit ...
!> \param section_label ...
!> \param n_char_size ...
!> \param n_dp_size ...
!> \param n_int_size ...
!> \param print_info ...
!> \param para_env ...
!> \param error ...
!> \par History
!>      - Creation (23.03.2011,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_binary_thermostats_nose(nhc,output_unit,section_label,&
                                           n_char_size,n_dp_size,n_int_size,&
                                           print_info,para_env,error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    INTEGER, INTENT(IN)                      :: output_unit
    CHARACTER(LEN=default_string_length), &
      INTENT(IN)                             :: section_label
    INTEGER, INTENT(INOUT)                   :: n_char_size, n_dp_size, &
                                                n_int_size
    LOGICAL, INTENT(IN)                      :: print_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'write_binary_thermostats_nose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, nhc_size
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eta, fnhc, mnhc, veta

    CALL timeset(routineN,handle)

    NULLIFY (eta)
    NULLIFY (fnhc)
    NULLIFY (mnhc)
    NULLIFY (veta)

    CALL collect_nose_restart_info(nhc,para_env,eta,veta,fnhc,mnhc,error)

    nhc_size = SIZE(eta)

    IF (output_unit > 0) THEN ! only ionode
       WRITE (UNIT=output_unit,IOSTAT=istat) section_label,nhc_size
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "nhc_size "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_char_size = n_char_size + LEN(section_label)
       n_int_size = n_int_size + 1
       IF (print_info) THEN
          WRITE (UNIT=*,FMT="(T3,A,T71,I10)")&
            "NHC size ("//TRIM(ADJUSTL(section_label))//")",nhc_size
       END IF
       ! eta
       WRITE (UNIT=output_unit,IOSTAT=istat) eta(1:nhc_size)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "eta(1:nhc_size) "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_dp_size = n_dp_size + nhc_size
    END IF ! ionode only

    DEALLOCATE (eta,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eta",para_env)

    ! veta
    IF (output_unit > 0) THEN ! only ionode
       WRITE (UNIT=output_unit,IOSTAT=istat) veta(1:nhc_size)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "veta(1:nhc_size) "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_dp_size = n_dp_size + nhc_size
    END IF ! ionode only

    DEALLOCATE (veta,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"veta",para_env)

    ! mnhc
    IF (output_unit > 0) THEN ! only ionode
       WRITE (UNIT=output_unit,IOSTAT=istat) mnhc(1:nhc_size)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "mnhc(1:nhc_size) "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_dp_size = n_dp_size + nhc_size
    END IF ! ionode only

    DEALLOCATE (mnhc,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mnhc",para_env)

    ! fnhc
    IF (output_unit > 0) THEN ! only ionode
       WRITE (UNIT=output_unit,IOSTAT=istat) fnhc(1:nhc_size)
       IF (istat /= 0) CALL stop_write(routineN,moduleN,__LINE__,&
                                       "fnhc(1:nhc_size) "//&
                                       "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                       output_unit,para_env)
       n_dp_size = n_dp_size + nhc_size
    END IF ! ionode only

    DEALLOCATE (fnhc,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"fnhc",para_env)

    CALL timestop(handle)

  END SUBROUTINE write_binary_thermostats_nose

END MODULE input_cp2k_restarts
