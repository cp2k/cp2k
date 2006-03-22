!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!


MODULE cp2k_runs

  USE bsse,                            ONLY: do_bsse_calculation
  USE cp_files,                        ONLY: close_file,&
                                             get_unit_number,&
                                             open_file
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release,&
                                             parser_search_string
  USE cp_units,                        ONLY: cp_unit_set_create,&
                                             cp_unit_set_release,&
                                             cp_unit_set_type
  USE eip_main,                        ONLY: eip_create_force_env
  USE environment,                     ONLY: compile_date,&
                                             compile_lastcvs,&
                                             cp2k_finalize,&
                                             cp2k_init,&
                                             cp2k_read,&
                                             cp2k_setup
  USE f77_interface,                   ONLY: f77_default_para_env => default_para_env
  USE farming_methods,                 ONLY: farming_parse_input
  USE farming_types,                   ONLY: deallocate_farming_env,&
                                             farming_env_type,&
                                             init_farming_env
  USE fist_debug,                      ONLY: fist_debug_control => debug_control
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_release,&
                                             force_env_set,&
                                             force_env_type
  USE geo_opt,                         ONLY: qs_geo_opt
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create,&
                                             globenv_release
  USE harris_metropolis_mc,            ONLY: harris_mc_ld
  USE input_constants,                 ONLY: &
       bsse_run, do_farming, do_fist, do_kg, do_pint, do_qmmm, do_qs, &
       do_tbmd, do_eip, do_test, electronic_spectra_run, energy_force_run, &
       energy_run, geo_opt_run, hmc_run, ldmc_run, linear_response_run, &
       mol_dyn_run, mon_car_run, none_run, qs_debug_run, KG_DFT_MC_method, &
       QS_MC_method, Harris_MC_method, EIP_MC_method, SE_MC_method, &
       TBMD_MC_method
  USE input_cp2k,                      ONLY: create_cp2k_root_section,&
                                             create_global_section
  USE input_cp2k_check,                ONLY: check_cp2k_input
  USE input_keyword_types,             ONLY: keyword_release
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: &
       section_describe_html, section_describe_index_html, section_release, &
       section_type, section_typo_match, section_vals_create, &
       section_vals_get_subs_vals, section_vals_release, section_vals_type, &
       section_vals_val_get, section_vals_write, write_section_xml
  USE kg_main,                         ONLY: kg_create_force_env
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE library_tests,                   ONLY: lib_test
  USE machine,                         ONLY: m_chdir,&
                                             m_getcwd
  USE mc_run,                          ONLY: do_mon_car
  USE md_run,                          ONLY: qs_mol_dyn
  USE message_passing,                 ONLY: mp_comm_free,&
                                             mp_comm_split,&
                                             mp_environ,&
                                             mp_sync
  USE metadynamics,                    ONLY: metadyn_read
  USE metadynamics_types,              ONLY: meta_env_release,&
                                             meta_env_type
  USE pint_methods,                    ONLY: do_pint_run
  USE qmmm_init,                       ONLY: qmmm_create_force_env
  USE qs_debug,                        ONLY: qs_debug_energy_and_forces
  USE qs_linres_module,                ONLY: linres_calculation
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE qs_tddfpt_module,                ONLY: tddfpt_calculation
  USE tbmd_main,                       ONLY: tbmd_create_force_env
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: cp2k_run, write_cp2k_html_manual, write_xml_file, run_input, &
             farming_run
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp2k_runs'
CONTAINS

!!****f* cp2k_runs/cp2k_run [1.0] *
!!
!!   NAME
!!     cp2k_run
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     performs an instance of a cp2k run
!!     cp2k_run_OK == .TRUE. if everything went fine
!!
!!   NOTES
!!     para_env should be a valid communicator
!!     output_unit should be writeable by at least the lowest rank of the mpi group
!!
!!     recursive because a given run_type might need to be able to perform 
!!     another cp2k_run as part of its job (e.g. farming, classical equilibration, ...)
!!
!!     the idea is that a cp2k instance should be able to run with just three
!!     arguments, i.e. a given input file, output unit, mpi communicator.
!!     giving these three to cp2k_run should produce a valid run.
!!     the only task of the PROGRAM cp2k is to create valid instances of the
!!     above arguments. Ideally, anything that is called afterwards should be
!!     able to run simultaneously / multithreaded / sequential / parallel / ...
!!     and able to fail safe
!!
!!   INPUTS
!!     input_file_name : name of the file to be opened for input 
!!     ouput_unit      : unit to which output should be written
!!     mpi_env         : the mpi environement for the run (the caller is
!!                       responsible to get rid of it)
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  RECURSIVE FUNCTION cp2k_run(input_file_name,output_unit,mpi_comm) RESULT(cp2k_run_OK)

    CHARACTER(LEN=*), INTENT(IN)             :: input_file_name
    INTEGER, INTENT(IN)                      :: output_unit, mpi_comm
    LOGICAL                                  :: cp2k_run_OK

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp2k_run', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=512)                       :: unknown_string
    INTEGER                                  :: unit_nr, &
                                                mc_method
    LOGICAL                                  :: echo_input, failure, &
                                                I_was_ionode, was_present
    TYPE(cp_error_type)                      :: error, suberror
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_parser_type), POINTER            :: my_parser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(section_type), POINTER              :: input_structure
    TYPE(section_vals_type), POINTER         :: glob_section, input_file, &
                                                metadyn_section, &
                                                hmc_section

    cp2k_run_OK = .FALSE.
    failure=.FALSE.
    was_present = .FALSE.
    NULLIFY(para_env, metadyn_section, default_units)
    CALL cp_para_env_create(para_env, group=mpi_comm,&
         owns_group=.FALSE.)

    NULLIFY(globenv, force_env )

    ! init globenv with the data
    CALL globenv_create(globenv, para_env=para_env,  output_unit=output_unit, &
         input_file_name=input_file_name, error=error)

    ! generate all dependent info in globenv
    CALL cp2k_init ( globenv )
    CALL cp_error_init(error, stop_level=cp_fatal_level)

    ! parse the input
    NULLIFY(input_structure)
    CALL create_cp2k_root_section(input_structure,error=error)
    NULLIFY(input_file)
    CALL section_vals_create(input_file,input_structure, error=error)
    NULLIFY(my_parser)
    CALL parser_create(my_parser,file_name=input_file_name,para_env=para_env,&
         error=error)
    CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
    CALL section_vals_parse(input_file,my_parser,root_section=.FALSE.,&
         default_units=default_units,error=error)
    CALL cp_unit_set_release(default_units,error=error)
    CALL parser_release(my_parser,error=error)
    globenv%input_file => input_file
    NULLIFY(input_file)

    ! put failure to true if this is what happened on parsing the input
    CALL cp_error_check(error,failure)

    IF (failure) THEN
       logger => cp_error_get_logger(error)   
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
       IF (unit_nr>0) THEN
          CALL cp_error_get(error,info=unknown_string)
          WRITE(unit_nr,'(T2,A)') ""
          WRITE(unit_nr,'(T2,A)') "Looking for words in the input similar to the unknown "//TRIM(unknown_string)
          WRITE(unit_nr,'(T2,A)') ""
          CALL section_typo_match(input_structure,unknown_string,location_string="", &
                                  unit_number=unit_nr,error=error)
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "CP2K failed to parse the input file."
          WRITE(unit_nr,'(T15,A)') "A full description of the input can be generated using"
          WRITE(unit_nr,'(T15,A)') "cp2k.sopt --html-manual"
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "This is also online for the latest version of CP2K:"
          WRITE(unit_nr,'(T15,A)') "http://cp2k.berlios.de/input/"
          WRITE(unit_nr,'(T15,A)') "or in single file format"
          WRITE(unit_nr,'(T15,A)') "http://cp2k.berlios.de/cp2k_input.htm "
          WRITE(unit_nr,'(T15,A)') ""
          WRITE(unit_nr,'(T15,A)') "If this input was an input of a previous version"
          WRITE(unit_nr,'(T15,A)') "of cp2k, you can try to convert it with --permissive-echo."
          WRITE(unit_nr,'(T15,A)') "However, this will just ignore the unknown keywords..."
          WRITE(unit_nr,'(T15,A)') ""
       ENDIF
    ELSE
       glob_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",error=error)
       CALL section_vals_val_get(glob_section,"ECHO_INPUT",l_val=echo_input,&
            error=error)
       IF (echo_input.AND.para_env%mepos==para_env%source) THEN
          CALL section_vals_write(globenv%input_file,unit_nr=globenv%scr,&
               hide_root=.TRUE., error=error)
       END IF

       IF (.NOT.failure) THEN
          CALL check_cp2k_input(globenv%input_file,para_env=para_env,error=error)
       END IF
       CALL cp_error_check(error,failure)

       IF (failure) THEN
          logger => cp_error_get_logger(error)
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.)
          IF (unit_nr>0) THEN
             WRITE(unit_nr,'(T15,A)') ""
             WRITE(unit_nr,'(T15,A)') "CP2K failed to parse the input."
             WRITE(unit_nr,'(T15,A)') "It was presumably generated by an incompatible CP2K version."
             WRITE(unit_nr,'(T15,A)') "try --permissive-echo to generate a compatible input file (to be checked!!!)"
             WRITE(unit_nr,'(T15,A)') ""
          ENDIF
       ELSE

          CALL cp_error_init(suberror)
          CALL cp2k_read(globenv)
          CALL cp2k_setup(globenv,suberror)

          IF ((globenv%run_type_id == hmc_run) .OR. &
              (globenv%run_type_id == ldmc_run)) THEN
            hmc_section => section_vals_get_subs_vals(globenv%input_file, &
                           "MOTION%MD%HMC", error=error)
            CALL section_vals_val_get(hmc_section, "MC_METHOD", &
                                      i_val=mc_method, error=error)
            SELECT CASE(mc_method)
            CASE(EIP_MC_method)
              globenv%prog_name_id = do_eip
            CASE(SE_MC_method)
              globenv%prog_name_id = do_qs
            CASE(TBMD_MC_method)
              globenv%prog_name_id = do_tbmd
            CASE(KG_DFT_MC_method)
              globenv%prog_name_id = do_kg
            CASE(QS_MC_method)
              globenv%prog_name_id = do_qs
            CASE(Harris_MC_method)
              globenv%prog_name_id = do_qs
            CASE DEFAULT
              PRINT *, "The standard DFT/QS method is used for HMC energy evaluation."
              globenv%prog_name_id = do_qs
            END SELECT
          END IF

          SELECT CASE (globenv%prog_name_id)
          CASE (do_pint)
             CALL do_pint_run(para_env,globenv%input_file,error=suberror)
          CASE (do_farming)
             CALL farming_run ( globenv, suberror )
          CASE (do_fist)
             CALL fist_create_force_env ( force_env, globenv, error = suberror )
          CASE (do_kg)
             CALL kg_create_force_env ( force_env, globenv, error=suberror)
          CASE (do_qmmm)
             CALL qmmm_create_force_env(force_env, globenv, error=suberror)
          CASE (do_qs)
             CALL quickstep_create_force_env(force_env, globenv, error=suberror)
          CASE (do_tbmd)
             CALL tbmd_create_force_env(force_env, globenv, error=suberror)
          CASE (do_eip)
             CALL eip_create_force_env(force_env, globenv, error=suberror)
          CASE (do_test)
             CALL lib_test(globenv,suberror)
          END SELECT
          IF (ASSOCIATED(force_env)) THEN
             NULLIFY(meta_env)
             metadyn_section => section_vals_get_subs_vals(globenv%input_file,"MOTION%METADYN",error=suberror)
             CALL metadyn_read(meta_env, globenv, metadyn_section, error=suberror)
             CALL force_env_set(force_env,meta_env=meta_env,error=suberror)
             CALL meta_env_release(meta_env,error=suberror)
          END IF

          SELECT CASE (globenv%run_type_id)
          CASE (none_run)
             ! do nothing
          CASE (energy_run, energy_force_run)
             IF(  globenv%prog_name_id /= do_qs .AND.&
                  globenv%prog_name_id /= do_kg .AND.&
                  globenv%prog_name_id  /= do_qmmm)  &
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Energy run available only with QS, KG, and QMMM programs "//&
                  CPSourceFileRef,suberror,failure)
             SELECT CASE(globenv%run_type_id)
             CASE (energy_run)
                CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,error=suberror)
             CASE(energy_force_run)
                CALL force_env_calc_energy_force(force_env,calc_force=.TRUE. ,error=suberror)
             END SELECT
          CASE (mol_dyn_run)
             CALL qs_mol_dyn ( force_env, globenv, error=suberror )
          CASE (geo_opt_run)
             CALL qs_geo_opt(force_env,globenv,error=suberror)
          CASE (mon_car_run)
             CALL do_mon_car ( force_env, globenv, error=suberror )
          CASE (hmc_run, ldmc_run)
             CALL harris_mc_ld(force_env, globenv, error=suberror)
          CASE (electronic_spectra_run)
             IF(globenv%prog_name_id /= do_qs) &
                  CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Electron spectra available only with Quickstep program "//&
                  CPSourceFileRef,suberror,failure)
             CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                  error=suberror)
             CALL tddfpt_calculation(force_env%qs_env, globenv)
          CASE (bsse_run)
             CALL do_bsse_calculation(force_env, globenv, error=suberror)
          CASE (linear_response_run)
             IF(globenv%prog_name_id /= do_qs) CALL cp_assert(.FALSE.,cp_failure_level,&
                  cp_assertion_failed,routineP,&
                  "Property calculations by Linear Response only within the QS program "//&
                  CPSourceFileRef,suberror,failure)
             ! The Ground State is needed, it can be read from Restart
             CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
                  error=suberror)
             CALL linres_calculation(force_env%qs_env, globenv, error=suberror)
          CASE (qs_debug_run)
             IF( (globenv%prog_name_id == do_qs) .OR. &
                  (globenv%prog_name_id == do_qmmm) ) THEN
                CALL qs_debug_energy_and_forces(force_env, globenv, error=suberror)
             ELSE IF(globenv%prog_name_id == do_fist) THEN
                CALL fist_debug_control(force_env, globenv, error=suberror)
             ELSE
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Debug run available only with QS, FIST, and QMMM program "//&
                     CPSourceFileRef,suberror,failure)
             END IF
          END SELECT
          CALL force_env_release(force_env, error=suberror)

       END IF
       CALL cp_error_propagate_error(suberror,routineP,&
               "cp2k ending run", error=error,failure=failure)
       CALL cp_error_dealloc_ref(suberror)

    END IF

    CALL section_release(input_structure,error=error)


    I_was_ionode=globenv%ionode
    CALL cp2k_finalize(globenv,error)
    CPPostconditionNoFail(globenv%ref_count==1,cp_failure_level,routineP,error)
    CALL cp_error_dealloc_ref(error)
    CALL globenv_release(globenv)
    CALL cp_para_env_release(para_env)
    cp2k_run_OK = .NOT. failure

  END FUNCTION cp2k_run

!!****f* cp2k_runs/farming_run [1.0] *
!!
!!   NAME
!!     farming_run
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     performs a farming run that performs several independent cp2k_runs
!!
!!   NOTES
!!     needs to be part of this module as the cp2k_run -> farming_run -> cp2k_run
!!     calling style creates a hard circular dependency
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE farming_run(globenv,error)
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_run', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=5)                         :: label
    CHARACTER(LEN=default_path_length)       :: output_file
    CHARACTER(LEN=default_string_length)     :: keyword, str
    INTEGER :: handle, i, i_job_to_restart, ierr, ijob, ijob_end, ijob_start, &
      iunit, n_jobs_to_run, new_group, new_rank, new_size, ngroups, &
      output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: group_distribution
    LOGICAL                                  :: failure, found, run_OK
    TYPE(cp_parser_type), POINTER            :: my_parser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(farming_env_type), POINTER          :: farming_env
    TYPE(section_type), POINTER              :: g_section
    TYPE(section_vals_type), POINTER         :: g_data

    failure=.FALSE.
    CALL timeset(routineN,'I','',handle)
    NULLIFY(my_parser,g_section,g_data,default_units)

    IF (globenv%ionode) WRITE(globenv%scr,FMT="(T2,A)") "FARMING| Hi, welcome on this farm!"

    ALLOCATE(farming_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL init_farming_env(farming_env)
    ! remember where we started
    CALL m_getcwd(farming_env%cwd)
    CALL farming_parse_input(farming_env,globenv,error)

    ! split the current communicator 
    ALLOCATE(group_distribution(0:globenv%num_pe-1),STAT=stat) 
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (farming_env%group_size_wish_set) THEN
       CALL mp_comm_split(globenv%group,new_group,ngroups,group_distribution,&
               subgroup_min_size=farming_env%group_size_wish)
    ELSE IF (farming_env%ngroup_wish_set) THEN
       CALL mp_comm_split(globenv%group,new_group,ngroups,group_distribution,&
               n_subgroups=farming_env%ngroup_wish,&
               group_partition=farming_env%group_partition)
    ELSE
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (globenv%ionode) THEN
        WRITE(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Number of created MPI groups:",ngroups
        WRITE(globenv%scr,FMT="(T2,A)",ADVANCE="NO") "FARMING| Task to group correspondence:"
        DO i=0,globenv%num_pe-1
           IF (MODULO(i,4)==0) WRITE(globenv%scr,*)
           WRITE(globenv%scr,FMT='(A3,I4,A3,I4,A1)',ADVANCE="NO") "  (",i," : ",group_distribution(i),")"
        END DO
        WRITE(globenv%scr,*)
    ENDIF

    ! protect about too many jobs being run in single go. Not more jobs are allowed than the number in the input file

    ! and determine the restart point
    IF (farming_env%cycle) THEN
       n_jobs_to_run=MIN(farming_env%njobs,farming_env%max_steps*ngroups)
       i_job_to_restart=MODULO(farming_env%restart_n+n_jobs_to_run-1,farming_env%njobs)+1
    ELSE
       n_jobs_to_run=MIN(farming_env%njobs,farming_env%max_steps*ngroups)
       n_jobs_to_run=MIN(n_jobs_to_run,farming_env%njobs-farming_env%restart_n+1)
       i_job_to_restart=n_jobs_to_run+farming_env%restart_n
    ENDIF

    ! and write it now, that's the point where the next job starts, even if this one is running
    IF (globenv%ionode) THEN
       iunit=get_unit_number()
       OPEN(UNIT=iunit,FILE="RESTARTFARM") 
       WRITE(iunit,*) i_job_to_restart
       CLOSE(UNIT=iunit)
    ENDIF

    ! this is the job range to be executed. Will be treated cycle
    ijob_start=farming_env%restart_n
    ijob_end=ijob_start+n_jobs_to_run-1

    IF (globenv%ionode) THEN
      IF (ijob_end-ijob_start<0) THEN
        WRITE(globenv%scr,FMT="(T2,A)") "FARMING| --- WARNING --- NO JOBS NEED EXECUTION ? "
        WRITE(globenv%scr,FMT="(T2,A)") "FARMING| is the cycle keyword required ?"
        WRITE(globenv%scr,FMT="(T2,A)") "FARMING| or is a stray RESTARTFARM present ?"
      ELSE
        WRITE(globenv%scr,FMT="(T2,A)") "FARMING| List of jobs : "
        DO ijob=ijob_start,ijob_end
         i=MODULO(ijob-1,farming_env%njobs)+1
         WRITE(globenv%scr,FMT=*) "Job: ",i," Dir: ",TRIM(farming_env%Job_cwd(i))," Input: ", &
           TRIM(farming_env%Job_input(i))," MPI group:", MODULO(i-1,ngroups)
        ENDDO
      ENDIF
    ENDIF

    DO ijob=ijob_start,ijob_end
       i=MODULO(ijob-1,farming_env%njobs)+1
       ! this farms out the jobs
       IF (MODULO(ijob-1,ngroups)==group_distribution(globenv%mepos)) THEN
          ! change to the new working directory
          CALL m_chdir(TRIM(farming_env%Job_cwd(i)),ierr)

          ! generate a fresh call to cp2k_run
          IF (ierr==0) THEN
              CALL mp_environ(new_size,new_rank,new_group)
              IF (new_rank == 0) THEN 

                 ! generate the output file
                  WRITE(output_file,'(A12,I5.5)') "FARMING_OUT_",i
                  CALL parser_create(my_parser,file_name=TRIM(farming_env%Job_input(i)))
                  label="&GLOBAL"
                  CALL parser_search_string(my_parser,label,ignore_case=.TRUE.,found=found,error=error)
                  IF (found) THEN
                     CALL create_global_section(g_section,error=error)
                     CALL section_vals_create(g_data,g_section,error=error)
                     CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
                     CALL section_vals_parse(g_data,my_parser,default_units,&
                          error=error)
                     CALL cp_unit_set_release(default_units,error=error)
                     CALL section_vals_val_get(g_data,"PROJECT",&
                          c_val=str, error=error)
                     IF (str.NE."") output_file=TRIM(str)//".out"
                     CALL section_vals_val_get(g_data,"OUTPUT_FILE_NAME",&
                             c_val=str,error=error)
                     IF (str.NE."") output_file=str
                     CALL section_vals_release(g_data)
                     CALL section_release(g_section)
                  END IF
                  CALL parser_release(my_parser,error=error)
                  CALL open_file(file_name=TRIM(output_file),&
                        file_action="WRITE",&
                        file_status="UNKNOWN",&
                        file_position="APPEND",&
                        unit_number=output_unit)
              ELSE
                ! this should never be used by cp2k_run, but e.g. still happens, incorrectly, 
                ! on eror, writes now to fort.999 so that error messages can be recovered and fixed
                output_unit=999 
             ENDIF
          
              IF (globenv%ionode) WRITE(globenv%scr,FMT="(T2,A,I5.5,A)") " Running Job ",i, &
                        " output to "//TRIM(farming_env%Job_cwd(i))//"/"//TRIM(output_file)
              run_OK=cp2k_run(TRIM(farming_env%Job_input(i)),output_unit,new_group)
              
              IF (new_rank == 0) THEN
                  CALL close_file(unit_number=output_unit)
              ELSE
              ENDIF

          ELSE ! problem if just one node does not change...
              IF (globenv%ionode) WRITE(globenv%scr,FMT="(T2,A,A)") &
                     " ... Failed to change dir to "//TRIM(farming_env%Job_cwd(i))
          ENDIF

          ! change to the original working directory
          CALL m_chdir(TRIM(farming_env%cwd),ierr)
          IF (ierr/=0) THEN
              IF (globenv%ionode) WRITE(globenv%scr,FMT="(T2,A,A)") &
                         " ... Failed to change back home to "//TRIM(farming_env%cwd)
          ENDIF

       ENDIF
    ENDDO

    ! give back the communicators of the split groups
    CALL mp_comm_free(new_group)
    DEALLOCATE(group_distribution,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! clean the farming env
    CALL deallocate_farming_env(farming_env)

    CALL mp_sync(globenv%group)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE farming_run

!!****f* cp2k_runs/write_cp2k_html_manual *
!!
!!   NAME
!!     write_cp2k_html_manual
!!
!!   FUNCTION
!!     writes a small html description of the cp2k input
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     joost [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE write_cp2k_html_manual(error)
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_cp2k_html_manual', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: root_section

    failure=.FALSE.

    IF (.NOT. failure) THEN
       NULLIFY(root_section)
       CALL create_cp2k_root_section(root_section,error=error)
       ! remove the default keyword that ignores things outside the section
       CALL keyword_release(root_section%keywords(0)%keyword,error=error)

       CALL section_describe_html(root_section,"InputReference",0,6,error=error)

       CALL open_file(unit_number=unit_nr,file_name="index.html",&
            file_action="WRITE", file_status="REPLACE")
       WRITE(unit_nr,*) "<HTML><BODY><HEAD><TITLE>The cp2k input structure</TITLE>"
       WRITE(unit_nr,*) "<H1>CP2K input reference</H1>"
       WRITE(unit_nr,*) "<H2>Version information</H2>"
       WRITE(unit_nr,*) "This html manual can be generated automatically from a given cp2k executable "// &
                        "using the --html-manual command line option. The manual describes exactly that version of the code. "//&
                        "This version has been generated using a cp2k version compiled on "//TRIM(compile_date)//". "//&
                        "The latest CVS entry found is "//TRIM(compile_lastcvs)//"."

       ! google search
       WRITE(unit_nr,*) '<H2>Search the input</H2>'
       WRITE(unit_nr,*) '<FORM method=GET action="http://www.google.com/search">'
       WRITE(unit_nr,*) '<input type=hidden name=ie value=UTF-8><input type=hidden name=oe value=UTF-8>'
       WRITE(unit_nr,*) '<TABLE bgcolor="#FFFFFF"><tr><td><A HREF="http://www.google.com/">'
       WRITE(unit_nr,*) '<IMG SRC="http://www.google.com/logos/Logo_25wht.gif" border="0" ALT="Google"></A>'
       WRITE(unit_nr,*) '</td><td><INPUT TYPE=text name=q size=31 maxlength=255 value="">'
       WRITE(unit_nr,*) '<INPUT type=submit name=btnG VALUE="CP2K Input Search">'
       WRITE(unit_nr,*) '<font size=-2><input type=hidden name=domains value="http://cp2k.berlios.de/input/">'
       WRITE(unit_nr,*) '<br><input type=radio name=sitesearch value=""> WWW'
       WRITE(unit_nr,*) '<input type=radio name=sitesearch value="http://cp2k.berlios.de/input/" checked>'
       WRITE(unit_nr,*) 'http://cp2k.berlios.de/input/ <br></font></td></tr></TABLE></FORM>'

       WRITE(unit_nr,*) "<H2>Input structure</H2>"
       WRITE(unit_nr,*) "All sections that can be part of a cp2k input file are shown with their allowed nestings. "
       WRITE(unit_nr,*) "A description of each section, and a list of keywords can be obtained clicking on the links. "
       WRITE(unit_nr,*) "<BR><UL>"

       CALL section_describe_index_html(root_section,"InputReference",unit_nr,error=error)

       WRITE(unit_nr,*) '</UL><BR><hr>Back to the <A HREF="../index.html">CP2K homepage</A>'
       WRITE(unit_nr,*) "</BODY></HTML>"

       CALL close_file(unit_nr)
       CALL section_release(root_section,error=error)
    END IF
  END SUBROUTINE write_cp2k_html_manual

! *****************************************************************************

  SUBROUTINE write_xml_file(error)

    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_xml_file', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, unit_number
    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: root_section

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    IF (.NOT.failure) THEN

      NULLIFY(root_section)

      CALL create_cp2k_root_section(root_section,error=error)

      CALL keyword_release(root_section%keywords(0)%keyword,error)

      CALL open_file(unit_number=unit_number,&
                     file_name="cp2k_input.xml",&
                     file_action="WRITE",&
                     file_status="REPLACE")

      WRITE (UNIT=unit_number,FMT="(A)")&
        "<?xml version=""1.0"" encoding=""ISO-8859-1""?>",&
        "<?xml-stylesheet type=""text/xsl"" href=""cp2k_input.xsl""?>",&
        "<CP2K_INPUT>",&
        " <COMPILE_DATE>"//TRIM(compile_date)//"</COMPILE_DATE>",&
        " <COMPILE_LASTCVS>"//TRIM(compile_lastcvs)//"</COMPILE_LASTCVS>"

      DO i=1,root_section%n_subsections
        CALL write_section_xml(root_section%subsections(i)%section,1,&
                               unit_number,error)
      END DO

      WRITE (UNIT=unit_number,FMT="(A)") "</CP2K_INPUT>"

      CALL close_file(unit_number=unit_number)
      CALL section_release(root_section,error=error)

    END IF

  END SUBROUTINE write_xml_file

! *****************************************************************************

  !!****f* f77_interface/run_input [1.0] *
  !!
  !!   NAME
  !!     run_input
  !!
  !!   FUNCTION
  !!     runs the given input
  !!
  !!   NOTES
  !!     moved here because of circular dependencies
  !!
  !!   INPUTS
  !!     - input_file_path: the path of the input file
  !!     - output_file_path: path of the output file (to which it is appended)
  !!       if it is "__STD_OUT__" the unit 6 is used
  !!     - mpi_comm: the mpi communicator to be used for this environment
  !!       it will not be freed
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE run_input(input_file_path,output_file_path,ierr,mpi_comm)
    CHARACTER(len=*), INTENT(in)             :: input_file_path, &
                                                output_file_path
    INTEGER, INTENT(out)                     :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'run_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure, success
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure=.FALSE.

    IF (PRESENT(mpi_comm)) THEN
       NULLIFY(para_env)
       CALL cp_para_env_create(para_env, group=mpi_comm, owns_group=.FALSE.)
    ELSE
       para_env => f77_default_para_env
       CALL cp_para_env_retain(para_env)
    END IF
    IF (para_env%mepos==para_env%source) THEN
       IF (output_file_path=="__STD_OUT__") THEN
          unit_nr=6
       ELSE
          CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=unit_nr)
       END IF
    ELSE
       unit_nr=-1
    END IF
    success=cp2k_run(input_file_path,unit_nr,para_env%group)
    IF (.NOT.success) THEN
       ierr=cp_failure_level
    ELSE
       ierr=0
    END IF
    CALL cp_para_env_release(para_env)
  END SUBROUTINE run_input
  !****************************************************************************

END MODULE cp2k_runs
!!*****
