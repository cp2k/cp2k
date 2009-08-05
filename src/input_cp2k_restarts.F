!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set of routines to dump the restart file of CP2K
!> \par History
!>      01.2006 [created] Teodoro Laino
! *****************************************************************************
MODULE input_cp2k_restarts
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE averages_types,                  ONLY: average_quantities_type
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cell_types,                      ONLY: cell_type,&
                                             real_to_scaled
  USE cp2k_info,                       ONLY: compile_arch,&
                                             compile_date,&
                                             compile_host,&
                                             compile_lastcvs
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
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             map_info_type,&
                                             npt_info_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type,&
                                             multiple_fe_list,&
                                             use_qmmm
  USE input_constants,                 ONLY: &
       do_band_collective, do_thermo_csvr, do_thermo_no_communication, &
       do_thermo_nose, ehrenfest, mol_dyn_run, mon_car_run, pint_run, &
       use_rt_restart
  USE input_section_types,             ONLY: &
       section_get_keyword_index, section_type, section_vals_add_values, &
       section_vals_get, section_vals_get_subs_vals, &
       section_vals_get_subs_vals3, section_vals_remove_values, &
       section_vals_type, section_vals_val_get, section_vals_val_set, &
       section_vals_val_unset, section_vals_write
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_datum,&
                                             m_getcwd
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_bcast,&
                                             mp_sum
  USE metadynamics_types,              ONLY: meta_env_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: get_molecule_kind
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE multipole_types,                 ONLY: multipole_type
  USE neb_types,                       ONLY: neb_var_type
  USE parallel_rng_types,              ONLY: dump_rng_stream,&
                                             rng_record_length
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_ff_fist,                    ONLY: qmmm_ff_precond_only_qm
  USE string_utilities,                ONLY: ascii_to_string,&
                                             string_to_ascii
  USE thermostat_types,                ONLY: thermostat_type
  USE thermostat_utils,                ONLY: communication_thermo_low2,&
                                             get_csvr_energies
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_restarts'

  PUBLIC :: write_restart,&
            update_subsys,&
            section_velocity_val_set,&
            write_restart_header

CONTAINS

! *****************************************************************************
!> \brief checks if a restart needs to be written and does so, updating all necessary fields
!>      in the input file. This is a relatively simple wrapper routine.
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE write_restart(md_env,force_env,kg_env,meta_env,root_section,&
       coords, vels,  error)
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(neb_var_type), OPTIONAL, POINTER    :: coords, vels
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=30), DIMENSION(2) :: &
      keys = (/"PRINT%RESTART_HISTORY", "PRINT%RESTART        "/)
    INTEGER                                  :: handle, ikey, ires, &
                                                nforce_eval
    LOGICAL                                  :: failure, save_mem
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
         motion_section,keys(1),error=error),cp_p_file) .OR.  &
         BTEST(cp_print_key_should_output(logger%iter_info,&
         motion_section,keys(2),error=error),cp_p_file) ) THEN

       sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
       CALL section_vals_get(sections, n_repetition=nforce_eval, error=error)

       CALL update_input(md_env,force_env,kg_env,meta_env,root_section,coords,vels,&
            save_mem=save_mem,error=error)

       DO ikey=1,SIZE(keys)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               motion_section,keys(ikey),error=error),cp_p_file)) THEN

             ires = cp_print_key_unit_nr(logger,motion_section,TRIM(keys(ikey)), extension=".restart",&
                                         file_position="REWIND", do_backup=(ikey==2), error=error)

             IF (ires>0) THEN
                CALL write_restart_header(ires, error)
                CALL section_vals_write(root_section,unit_nr=ires,hide_root=.TRUE., error=error)
             ENDIF

             CALL cp_print_key_finished_output(ires,logger,motion_section,&
                  TRIM(keys(ikey)), error=error)
          ENDIF
       ENDDO

       IF(save_mem) THEN
          CALL update_subsys_release(md_env,force_env,root_section,error)
          CALL update_motion_release(motion_section,error)
       END IF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE write_restart

! *****************************************************************************
!> \brief Writes the header for the restart file
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      01.2008 [created] - Splitted from write_restart
!> \author Teodoro Laino - University of Zurich - 01.2008
! *****************************************************************************
  SUBROUTINE write_restart_header(ires, error)
    INTEGER, INTENT(IN)                      :: ires
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_restart_header', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=256)                       :: cwd, datx
    INTEGER                                  :: strend, strstart

    CALL m_datum(datx)
    CALL m_getcwd(cwd)

    WRITE (UNIT=ires,FMT="(T2,A)") "# Version information for this restart file "
    WRITE (UNIT=ires,FMT="(T2,A)") "# current date "//TRIM(datx)
    WRITE (UNIT=ires,FMT="(T2,A)") "# current working dir "//TRIM(cwd)

    WRITE (UNIT=ires,FMT="(T2,A,T31,A50)")&
         "# Program compiled at",&
         ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
    WRITE (UNIT=ires,FMT="(T2,A,T31,A50)")&
         "# Program compiled on",&
         ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
    WRITE (UNIT=ires,FMT="(T2,A,T31,A50)")&
         "# Program compiled for",&
         ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
    strend=MAX(1,LEN(compile_lastcvs)-2)
    strstart=MAX(2,strend-50)
    WRITE (UNIT=ires,FMT="(T2,A,T31,A50)")&
         "# Last CVS entry",&
         ADJUSTR(compile_lastcvs(strstart:strend))

  END SUBROUTINE write_restart_header

! *****************************************************************************
!> \brief deallocate some sub_sections of the section subsys to save some memory 
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

          IF(ASSOCIATED(shell_particles)) THEN
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

          IF(ASSOCIATED(core_particles)) THEN
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
       DEALLOCATE(i_force_eval)
    END IF
    CALL timestop(handle)

  END SUBROUTINE update_subsys_release

! *****************************************************************************
!> \brief deallocate the nose subsections (coord, vel, force, mass) in the md section
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
    NULLIFY( work_section)

    work_section => section_vals_get_subs_vals(motion_section,"MD%AVERAGES%RESTART_AVERAGES",error=error)
    CALL section_vals_remove_values(work_section, error)

    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section, error)

    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%BAROSTAT%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section, error)
    
    work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section, error)

    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%COORD",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%VELOCITY",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%FORCE",error=error)
    CALL section_vals_remove_values(work_section, error)
    work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE%MASS",error=error)
    CALL section_vals_remove_values(work_section, error)

    CALL timestop(handle)
  END SUBROUTINE   update_motion_release

! *****************************************************************************
!> \brief Updates the whole input file for the restart
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_input(md_env,force_env,kg_env,meta_env,root_section,coords,vels,save_mem,error)
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), OPTIONAL, POINTER  :: force_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(meta_env_type), OPTIONAL, POINTER   :: meta_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(neb_var_type), OPTIONAL, POINTER    :: coords, vels
    LOGICAL, INTENT(IN), OPTIONAL            :: save_mem
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, lcond, my_save_mem
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(force_env_type), POINTER            :: my_force_env
    TYPE(kg_environment_type), POINTER       :: my_kg_env
    TYPE(meta_env_type), POINTER             :: my_meta_env
    TYPE(section_vals_type), POINTER         :: motion_section

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(logger,motion_section,my_force_env,my_kg_env,my_meta_env)
    my_save_mem = .FALSE.
    IF(PRESENT(save_mem)) my_save_mem=save_mem

    logger => cp_error_get_logger(error)
    !
    ! Can handle md_env or force_env
    !
    lcond = PRESENT(md_env).OR.PRESENT(force_env).OR.PRESENT(kg_env)
    IF (lcond) THEN
       IF (PRESENT(md_env)) THEN
          CALL get_md_env ( md_env=md_env, force_env=my_force_env, error=error)
       ELSE IF (PRESENT(force_env)) THEN 
          my_force_env => force_env
       END IF
       !
       ! Let's set up kg_env
       !
       IF (PRESENT(kg_env)) THEN
          my_kg_env => kg_env
       ELSEIF (ASSOCIATED(my_force_env)) THEN 
          my_kg_env => my_force_env%kg_env
       END IF
       !
       ! meta_env
       !
       IF (PRESENT(meta_env)) THEN
          my_meta_env => meta_env
       ELSEIF (ASSOCIATED(my_force_env)) THEN 
          my_meta_env => my_force_env%meta_env
       END IF
       !
       ! The real restart setting...
       !
       motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
       CALL update_motion(motion_section,md_env,my_force_env,my_kg_env,my_meta_env,&
            logger,coords,vels,save_mem=my_save_mem,error=error)
       !
       ! Update one force_env_section per time..
       !
       IF (ASSOCIATED(my_force_env)) &
            CALL update_force_eval(force_env=my_force_env,logger=logger,&
            root_section=root_section,error=error)
       
    END IF
    CALL timestop(handle)
       
  END SUBROUTINE update_input

! *****************************************************************************
!> \brief Updates the motion section of the input file
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_motion(motion_section,md_env,force_env,kg_env,meta_env,logger,&
       coords,vels,save_mem,error)
    TYPE(section_vals_type), POINTER         :: motion_section
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(neb_var_type), OPTIONAL, POINTER    :: coords, vels
    LOGICAL, INTENT(IN), OPTIONAL            :: save_mem
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_motion', &
      routineP = moduleN//':'//routineN

    INTEGER :: counter, handle, handle2, i, iatom, icoef, iel, ikind, irep, &
      isec, j, n_el, natom, natom_tot, nfield, nhc_len, nkind, nsgf, &
      nsgf_max, stat, tot_nhcneed
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cindex
    INTEGER, DIMENSION(:), POINTER           :: lindex, walkers_status
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, my_save_mem
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pbuffer, pindex, vbuffer, &
                                                vindex
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buffer, eta, fnhc, mnhc, &
                                                veta, wrk
    REAL(KIND=dp), POINTER                   :: constant, t
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: replica_section, work_section
    TYPE(thermostat_type), POINTER           :: thermostat_baro, &
                                                thermostat_coeff, &
                                                thermostat_part, &
                                                thermostat_shell

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY (logger, thermostat_part, thermostat_baro, npt, para_env, nhc, &
         atomic_kind_set, work_section, dyn_coeff_set, atomic_kind, coeffs, &
         aux_basis_set, thermostat_coeff, thermostat_shell, t, averages, &
         constant, walkers_status, itimes)
    IF (PRESENT(md_env)) THEN
       CALL get_md_env ( md_env=md_env,&
                         thermostat_part=thermostat_part,&
                         thermostat_baro=thermostat_baro,&
                         thermostat_coeff=thermostat_coeff,&
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
       ELSEIF (ASSOCIATED(kg_env)) THEN
          para_env => kg_env%para_env
       ELSE
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    
    my_save_mem=.FALSE.
    IF (PRESENT(save_mem)) my_save_mem=save_mem

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
       IF     (thermostat_part%type_of_thermostat==do_thermo_nose) THEN
          ! Restart of Nose-Hoover Thermostat for Particles
          nhc     => thermostat_part%nhc
          CALL collect_nose_restart_info(nhc, para_env, eta, veta, fnhc, mnhc, error)
          
          work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%NOSE",error=error)
          CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)
       ELSE IF (thermostat_part%type_of_thermostat==do_thermo_csvr) THEN
          ! Restart of CSVR Thermostat for Particles
          work_section => section_vals_get_subs_vals(motion_section,"MD%THERMOSTAT%CSVR",error=error)
          CALL dump_csvr_restart_info(thermostat_part%csvr, para_env, work_section, error)
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
    ! COEFFICIENTS 
    CALL timeset(routineN//"_KG",handle2)
    IF (ASSOCIATED(kg_env)) THEN
       CALL get_kg_env(kg_env=kg_env,&
            atomic_kind_set=atomic_kind_set,&
            particle_set=particle_set, &
            dyn_coeff_set=dyn_coeff_set,&
            error=error)
       
       IF (ASSOCIATED(dyn_coeff_set)) THEN          
          natom_tot = SIZE ( particle_set )
          nsgf_max=0
          n_el=0
          DO ikind=1,SIZE(atomic_kind_set)
             atomic_kind=> atomic_kind_set(ikind)
             CALL get_atomic_kind(atomic_kind=atomic_kind,&
                  aux_basis_set=aux_basis_set,&
                  natom=natom)
             IF (.NOT.ASSOCIATED(aux_basis_set)) CYCLE
             CALL get_gto_basis_set(gto_basis_set=aux_basis_set,&
                  nsgf=nsgf)
             nsgf_max=MAX(nsgf_max,nsgf)
             n_el = n_el + natom
          END DO
          
          ALLOCATE (pbuffer(n_el,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (vbuffer(n_el,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (lindex(n_el),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (cindex(natom_tot),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (pindex(natom_tot,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          ALLOCATE (vindex(natom_tot,nsgf_max),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          
          pindex(:,:)=0.0_dp
          vindex(:,:)=0.0_dp
          cindex(:)=0

          nkind=SIZE(dyn_coeff_set%coeffs_of_kind)
          DO ikind=1, nkind
             coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
             atomic_kind=> atomic_kind_set(ikind)
             IF(.NOT.ASSOCIATED(coeffs)) CYCLE
             DO iel= 1,coeffs%n_els
                iatom=dyn_coeff_set%distribution%list(ikind)%array(iel)
                cindex(iatom)=iatom
                DO icoef=1,coeffs%ncoef_atom
                   pindex(iatom,icoef)=coeffs%pos(iel,icoef)
                   vindex(iatom,icoef)=coeffs%vel(iel,icoef)
                END DO
             END DO
          END DO
          CALL mp_sum(pindex,para_env%group)
          CALL mp_sum(vindex,para_env%group)
          CALL mp_sum(cindex, para_env%group)

          iel = 0
          DO iatom = 1, natom_tot
            IF ( cindex ( iatom ) == 0 ) THEN
              CYCLE 
            ELSE
              iel = iel + 1
              lindex( iel ) = cindex ( iatom )
              pbuffer ( iel, : ) = pindex ( iatom, : )
              vbuffer ( iel, : ) = vindex ( iatom, : )
            END IF
          END DO
          
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%RESTART_INDEX",error=error)
          CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",i_vals_ptr=lindex,error=error)
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%COORD",error=error)
          DO irep = 1, SIZE(lindex)
             ALLOCATE(buffer(SIZE(pbuffer,2)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = pbuffer(irep,:)
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_vals_ptr=buffer,&
                  i_rep_val=irep,error=error)    
          END DO
          work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%VELOCITY",error=error)
          DO irep = 1, SIZE(lindex)
             ALLOCATE(buffer(SIZE(vbuffer,2)),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             buffer = vbuffer(irep,:)       
             CALL section_vals_val_set(work_section,"_DEFAULT_KEYWORD_",r_vals_ptr=buffer,&
                  i_rep_val=irep,error=error)    
          END DO
          DEALLOCATE (pbuffer,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          DEALLOCATE (vbuffer,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)   
          DEALLOCATE (cindex,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (pindex,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (vindex,STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(thermostat_coeff)) THEN
          IF     (thermostat_coeff%type_of_thermostat==do_thermo_nose) THEN
             ! Restart of Nose-Hoover Thermostat for KG coefficients
             nhc         => thermostat_coeff%nhc
             CALL collect_nose_restart_info(nhc, para_env, eta, veta, fnhc, mnhc, error)

             work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT%NOSE",error=error)
             CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)    
          ELSE IF (thermostat_coeff%type_of_thermostat==do_thermo_csvr) THEN
             work_section => section_vals_get_subs_vals(motion_section,"MD%KG_COEFF%THERMOSTAT%CSVR",error=error)
             ! Restart of CSVR Thermostat for KG coefficients
             CALL dump_csvr_restart_info(thermostat_coeff%csvr, para_env, work_section, error)
          END IF
       END IF
    END IF
    CALL timestop(handle2)
    ! SHELL THERMOSTAT
    CALL timeset(routineN//"_THERMOSTAT_SHELL",handle2)
    IF (ASSOCIATED(thermostat_shell)) THEN
       IF     (thermostat_shell%type_of_thermostat==do_thermo_nose) THEN
          ! Restart of Nose-Hoover Thermostat for Shell Particles
          nhc     => thermostat_shell%nhc
          CALL collect_nose_restart_info(nhc, para_env, eta, veta, fnhc, mnhc, error)

          work_section => section_vals_get_subs_vals(motion_section,"MD%SHELL%THERMOSTAT%NOSE",error=error)
          CALL set_template_restart(work_section, eta, veta, fnhc, mnhc, error)    
       ELSE IF (thermostat_shell%type_of_thermostat==do_thermo_csvr) THEN
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
       IF (meta_env%do_hills.AND.meta_env%hills_env%n_hills/=0) THEN
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_POS",error=error)
          CALL meta_hills_val_set_ss(work_section, meta_env, error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_AMPL",error=error)
          CALL meta_hills_val_set_ds(work_section, meta_env, error=error)
          work_section => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_HEIGHT",error=error)
          CALL meta_hills_val_set_ww(work_section, meta_env, error=error)
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
             CALL section_neb_coord_val_set(work_section,coords%xyz(:,isec),SIZE(coords%xyz,1),3,error)
             ! Update Collective Variables
             IF (coords%in_use==do_band_collective) THEN
                ALLOCATE(wrk(coords%size_wrk(1)),stat=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
                wrk = coords%wrk(:,isec)
                CALL section_vals_val_set(replica_section,"COLLECTIVE",r_vals_ptr=wrk,&
                     i_rep_section=isec,error=error) 
             END IF
          END DO
       END IF
       IF (PRESENT(vels)) THEN
          ! Allocate possible missing sections
          DO
             IF (vels%size_wrk(2)<=SIZE(replica_section%values,2)) EXIT
             CALL section_vals_add_values(replica_section,error=error)
          END DO
          nfield = 3
          IF (vels%in_use==do_band_collective) nfield = 1
          ! Write Values
          DO isec = 1, vels%size_wrk(2)
             work_section => section_vals_get_subs_vals3(replica_section,"VELOCITY",i_rep_section=isec,error=error)
             CALL section_neb_coord_val_set(work_section,vels%wrk(:,isec),SIZE(vels%wrk,1),nfield,error)
          END DO
       END IF
    END IF
    CALL timestop(handle2)
    CALL timestop(handle)

  END SUBROUTINE update_motion

! *****************************************************************************
!> \brief routine to dump thermostat CSVR energies 
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

    CHARACTER(len=*), PARAMETER :: routineN = 'dump_csvr_energy_info', &
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
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.    
    ! Thermostat Energies
    ALLOCATE(work(csvr%glob_num_csvr),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL get_csvr_energies ( csvr, dum, para_env, array_kin=work, error=error )

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
!> \brief Collect all information needed to dump the restart for Nose-Hoover
!>      thermostat
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
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: work
    TYPE(map_info_type), POINTER             :: map_info

    nhc_len = SIZE(nhc%nvt,1)
    num_nhc = nhc%loc_num_nhc
    numneed = num_nhc
    map_info=> nhc%map_info
    ALLOCATE ( nhc_per_proc(para_env%num_pe), STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    nhc_per_proc(:) = 0

    CALL mp_allgather(numneed,nhc_per_proc,para_env%group)
    tot_nhcneed = nhc%glob_num_nhc

    NULLIFY(work,index)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! nhc%eta  
    !-----------------------------------------------------------------------------
    ALLOCATE(eta(tot_nhcneed*nhc_len), STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
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
    ALLOCATE(veta(tot_nhcneed*nhc_len), STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
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
    ALLOCATE(fnhc(tot_nhcneed*nhc_len), STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
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
    ALLOCATE(mnhc(tot_nhcneed*nhc_len), STAT = stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)    
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

    DEALLOCATE ( work, STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
    DEALLOCATE ( index, STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
    DEALLOCATE ( nhc_per_proc, STAT = stat )
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 

  END SUBROUTINE collect_nose_restart_info

! *****************************************************************************
!> \brief routine to dump NEB coordinates and velocities section.. fast implementation
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      12.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_neb_coord_val_set(coord_section, array, nsize, nfield,error)

    TYPE(section_vals_type), POINTER         :: coord_section
    REAL(KIND=dp), DIMENSION(*)              :: array
    INTEGER, INTENT(IN)                      :: nsize, nfield
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_neb_coord_val_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ik, irk, Nlist, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_c
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
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
       i = (irk-1)*nfield
       ALLOCATE(my_c(nfield),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       my_c = array(i+1:i+nfield)
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
    NULLIFY(coord, force, velocity, mass)
    IF (PRESENT(eta)) THEN
       IF (SIZE(eta)>0) THEN
          coord => section_vals_get_subs_vals(work_section,"COORD",error=error)
          CALL section_vals_val_set(coord,"_DEFAULT_KEYWORD_",r_vals_ptr=eta,error=error)
       ELSE
          DEALLOCATE(eta, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
    IF (PRESENT(veta)) THEN
       IF (SIZE(veta)>0) THEN
          velocity => section_vals_get_subs_vals(work_section,"VELOCITY",error=error)
          CALL section_vals_val_set(velocity,"_DEFAULT_KEYWORD_",r_vals_ptr=veta,error=error)
       ELSE
          DEALLOCATE(veta, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
    IF (PRESENT(fnhc)) THEN
       IF (SIZE(fnhc)>0) THEN
          force => section_vals_get_subs_vals(work_section,"FORCE",error=error)
          CALL section_vals_val_set(force,"_DEFAULT_KEYWORD_",r_vals_ptr=fnhc,error=error)
       ELSE
          DEALLOCATE(fnhc, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
    IF (PRESENT(mnhc)) THEN
       IF (SIZE(mnhc)>0) THEN
          mass => section_vals_get_subs_vals(work_section,"MASS",error=error)
          CALL section_vals_val_set(mass,"_DEFAULT_KEYWORD_",r_vals_ptr=mnhc,error=error)
       ELSE
          DEALLOCATE(mnhc, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure) 
       END IF
    END IF
  END SUBROUTINE set_template_restart

! *****************************************************************************
!> \brief Updates the force_eval section of the input file
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_force_eval(force_env,logger,root_section,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: root_section
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
    CALL update_subsys(subsys_section,force_env,skip_vel_section,error)

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
                        &rng(iparticle_local)%stream,&
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
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \par History
!>      01.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE update_subsys(subsys_section,force_env,skip_vel_section,error)
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(force_env_type), POINTER            :: force_env
    LOGICAL, INTENT(IN)                      :: skip_vel_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'update_subsys', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: unit_str
    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:), POINTER           :: multiple_unit_cell
    LOGICAL                                  :: failure, scale, use_ref_cell
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(section_vals_type), POINTER         :: work_section

    failure = .FALSE.
    NULLIFY(work_section, core_particles, particles, shell_particles, &
         subsys, cell, cell_ref, multipoles)
    CALL timeset(routineN,handle)
    CALL force_env_get(force_env, subsys=subsys, cell=cell, cell_ref=cell_ref,&
         use_ref_cell=use_ref_cell, error=error)

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
    CALL section_coord_val_set(work_section, particles, molecules, conv_factor, scale,&
         cell, error=error)

    work_section => section_vals_get_subs_vals(subsys_section,"VELOCITY",error=error)
    IF (.NOT.skip_vel_section) THEN
       CALL section_velocity_val_set(work_section, particles, conv_factor=1.0_dp, error=error)
    ELSE
       CALL section_vals_remove_values(work_section, error)
    END IF
    ! write restart input for shell-model
    IF(ASSOCIATED(shell_particles)) THEN
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
    ENDIF
    IF(ASSOCIATED(core_particles)) THEN
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
    ENDIF
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

    CHARACTER(len=*), PARAMETER :: routineN = 'section_coord_val_set', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: line
    CHARACTER(LEN=default_string_length)     :: icval, my_tag, name, rcval
    INTEGER                                  :: handle, ik, imol, irk, k, &
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
       IF(my_shell) THEN
         line = TRIM(name)
         s = particles%els(irk)%r
         IF (scale) THEN
            r = s   
            CALL real_to_scaled(s,r,cell)
         ELSE
            s = s*conv_factor
         END IF
         DO k = 1, 3
            WRITE(rcval,"(ES25.16)") s(k)
            line = TRIM(line)//" "//TRIM(rcval)
         END DO
         WRITE(icval,"(I10)") particles%els(irk)%atom_index
         line = TRIM(line)//" "//TRIM(icval)
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
         line = TRIM(name)
         s = particles%els(irk)%r
         IF (scale) THEN
            r = s   
            CALL real_to_scaled(s,r,cell)
         ELSE
            s = s*conv_factor
         END IF
         DO k = 1, 3
            WRITE(rcval,"(ES25.16)") s(k)
            line = TRIM(line)//" "//TRIM(rcval)
         END DO
         line = TRIM(line)//" "//TRIM(my_tag)
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
!> \brief routine to dump velocities.. fast implementation
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2006 created [teo]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE section_velocity_val_set(velocity_section, particles, velocity, conv_factor, error)

    TYPE(section_vals_type), POINTER         :: velocity_section
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: velocity
    REAL(KIND=dp)                            :: conv_factor
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_velocity_val_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ik, irk, Nlist, &
                                                nloop, stat
    LOGICAL                                  :: check, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vel
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: my_val, old_val

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(my_val, old_val, section, vals)
    CPPrecondition(ASSOCIATED(velocity_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(velocity_section%ref_count>0,cp_failure_level,routineP,error,failure)
    section => velocity_section%section
    ik=section_get_keyword_index(section,"_DEFAULT_KEYWORD_",error=error)
    CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
         "section "//TRIM(section%name)//" does not contain keyword "//&
         "_DEFAULT_KEYWORD_",error,failure)

    ! At least one of the two arguments must be present..
    check = PRESENT(particles).NEQV.PRESENT(velocity) 
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    IF(PRESENT(particles)) nloop = particles%n_els
    IF(PRESENT(velocity))  nloop = SIZE(velocity,2)
    
    DO
       IF (SIZE(velocity_section%values,2)==1) EXIT
       CALL section_vals_add_values(velocity_section,error=error)
    END DO
    vals  => velocity_section%values(ik,1)%list
    Nlist = 0
    IF (ASSOCIATED(vals)) THEN
       Nlist = cp_sll_val_get_length(vals,error)
    END IF
    DO irk=1,nloop
       ALLOCATE(vel(3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Always stored in A.U.
       IF(PRESENT(particles)) vel = particles%els(irk)%v(1:3)*conv_factor
       IF(PRESENT(velocity))  vel = velocity(1:3,irk)*conv_factor

       CALL val_create(my_val,r_vals_ptr=vel,error=error)

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
    velocity_section%values(ik,1)%list => vals

    CALL timestop(handle)
  END SUBROUTINE section_velocity_val_set

! *****************************************************************************
!> \brief routine to dump dipoles.. fast implementation
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

    CHARACTER(len=*), PARAMETER :: routineN = 'update_dipoles_section', &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'update_quadrupoles_section', &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'section_rng_val_set', &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_hills_val_set_ss', &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_hills_val_set_ds', &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_hills_val_set_ww', &
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

END MODULE input_cp2k_restarts

