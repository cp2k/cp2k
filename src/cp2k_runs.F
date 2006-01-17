!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!


MODULE cp2k_runs

  USE atom_module,                     ONLY: atomic_code
  USE cp_files,                        ONLY: open_file, &
                                             close_file, &
                                             get_unit_number
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_parser,                       ONLY: cp_parser_type, &
                                             parser_create, &
                                             parser_release,&
                                             parser_search_string
  USE cp_units,                        ONLY: cp_unit_set_type, &
                                             cp_unit_set_create, &
                                             cp_unit_set_release
  USE environment,                     ONLY: cp2k_finalize,&
                                             cp2k_init,&
                                             cp2k_read,&
                                             cp2k_setup,&
                                             compile_date,compile_lastcvs
  USE farming_methods,                 ONLY: farming_parse_input
  USE farming_types,                   ONLY: farming_env_type, &
                                             init_farming_env, &
                                             deallocate_farming_env
  USE fist_debug,                      ONLY: fist_debug_control => debug_control
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_set, force_env_get
  USE geo_opt,                         ONLY: qs_geo_opt
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create,&
                                             globenv_retain,&
                                             globenv_release
  USE harris_metropolis_mc,            ONLY: harris_mc_ld
  USE input_constants,                 ONLY: do_qs,&
                                             do_fist,&
                                             do_kg,&
                                             do_qmmm,&
                                             do_atom,&
                                             do_farming,&
                                             do_tbmd,&
                                             do_test,&
                                             none_run,&
                                             energy_run,&
                                             energy_force_run,&
                                             mol_dyn_run,&
                                             geo_opt_run,&
                                             mon_car_run,&
                                             hmc_run,&
                                             ldmc_run,&
                                             electronic_spectra_run,&
                                             bsse_run,&
                                             linear_response_run,&
                                             qs_debug_run
  USE input_cp2k,                      ONLY: create_cp2k_root_section,&
                                             check_cp2k_input, &
                                             create_global_section
  USE input_keyword_types,             ONLY: keyword_release
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: section_type, &
                                             section_release,&
                                             section_vals_type, &
                                             section_vals_release, &
                                             section_vals_create,&
                                             section_vals_write, &
                                             section_vals_get, &
                                             section_vals_get_subs_vals, &
                                             section_describe_html,&
                                             section_vals_val_get, &
                                             section_describe_index_html, &
                                             write_section_xml
  USE kg_main,                         ONLY: kg_create_force_env
  USE kinds,                           ONLY: default_string_length, &
                                             default_path_length,&
                                             dp
  USE library_tests,                   ONLY: lib_test
  USE machine,                         ONLY: print_memory,m_chdir,m_getcwd
  USE mc_run,                          ONLY: do_mon_car
  USE md_run,                          ONLY: qs_mol_dyn 
  USE message_passing,                 ONLY: mp_sync, &
                                             mp_environ, &
                                             mp_comm_split, &
                                             mp_comm_free
  USE metadynamics_types,              ONLY: meta_env_type, &
                                             meta_env_release
  USE metadynamics,                    ONLY: metadyn_read
  USE print_key_types,                 ONLY: print_keys_release
  USE qmmm_init,                       ONLY: qmmm_create_force_env
  USE qs_debug,                        ONLY: qs_debug_energy_and_forces
  USE qs_linres_module,                ONLY: linres_calculation
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE qs_tddfpt_module,                ONLY: tddfpt_calculation
  USE string_utilities,                ONLY: uppercase
  USE tbmd_main,                       ONLY: tbmd_create_force_env
  USE timings,                         ONLY: timeset, &
                                             timestop
  USE bsse,                            ONLY: do_bsse_calculation
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: cp2k_run, write_cp2k_html_manual, write_xml_file
  CHARACTER(LEN=*), PARAMETER :: moduleN = "cp2k_run"
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

    LOGICAL                              :: cp2k_run_OK

    CHARACTER(LEN=*), INTENT(IN) :: input_file_name
    INTEGER         , INTENT(IN) :: output_unit, mpi_comm

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp2k_run",&
         routineP=moduleN//':'//routineN  

    TYPE(global_environment_type), pointer :: globenv
    TYPE(cp_error_type) :: error
    TYPE(force_env_type), pointer :: force_env
    TYPE(meta_env_type), pointer :: meta_env
    TYPE(section_type), POINTER :: input_structure
    TYPE(section_vals_type), POINTER :: input_file, glob_section, metadyn_section
    TYPE(cp_parser_type), pointer :: my_parser
    TYPE(cp_para_env_type), pointer :: para_env
    TYPE(cp_unit_set_type), pointer :: default_units
    LOGICAL :: failure, echo_input, I_was_ionode, was_present

    cp2k_run_OK = .FALSE.
    failure=.false.
    was_present = .FALSE.
    NULLIFY(para_env, metadyn_section, default_units)
    CALL cp_para_env_create(para_env, group=mpi_comm,&
         owns_group=.false.)

    NULLIFY(globenv, force_env )

    ! init globenv with the data
    CALL globenv_create(globenv, para_env=para_env,  output_unit=output_unit, &
         input_file_name=input_file_name, error=error)

    ! generate all dependent info in globenv
    CALL cp2k_init ( globenv )
    CALL cp_error_init(error)

    ! parse the input
    nullify(input_structure)
    CALL create_cp2k_root_section(input_structure,error=error)
    NULLIFY(input_file)
    CALL section_vals_create(input_file,input_structure, error=error)
    CALL section_release(input_structure,error=error)
    nullify(my_parser)
    CALL parser_create(my_parser,file_name=input_file_name,para_env=para_env,&
         error=error)
    CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
    CALL section_vals_parse(input_file,my_parser,root_section=.false.,&
         default_units=default_units,error=error)
    CALL cp_unit_set_release(default_units,error=error)
    CALL parser_release(my_parser,error=error)
    CALL cp_error_check(error,failure)
    glob_section => section_vals_get_subs_vals(input_file,"GLOBAL",error=error)
    CALL section_vals_val_get(glob_section,"ECHO_INPUT",l_val=echo_input,&
         error=error)
    IF (echo_input.AND.para_env%mepos==para_env%source) THEN
       CALL section_vals_write(input_file,unit_nr=globenv%scr,&
            hide_root=.TRUE., error=error)
    END IF

    IF (.NOT.failure) THEN
       CALL check_cp2k_input(input_file,error=error)
    END IF
    globenv%input_file => input_file
    nullify(input_file)
    CALL cp_error_check(error,failure)

    IF (.NOT.failure) THEN
       CALL cp2k_read(globenv)
       CALL cp2k_setup(globenv,error)
       SELECT CASE (globenv%prog_name_id)
       CASE (do_atom)
          CALL atomic_code ( globenv )
       CASE (do_farming)
          CALL farming_run ( globenv, error )
       CASE (do_fist)
          CALL fist_create_force_env ( force_env, globenv,  error = error )
       CASE (do_kg)
          CALL kg_create_force_env ( force_env, globenv, error=error)
       CASE (do_qmmm)
          CALL qmmm_create_force_env(force_env, globenv, error=error)
       CASE (do_qs)
          CALL quickstep_create_force_env(force_env, globenv, error=error)
       CASE (do_tbmd)
          CALL tbmd_create_force_env(force_env, globenv, error=error)
       CASE (do_test)
          CALL lib_test(globenv,error)
       END SELECT
       IF (ASSOCIATED(force_env)) THEN
          NULLIFY(meta_env)
          metadyn_section => section_vals_get_subs_vals(globenv%input_file,"MOTION%METADYN",error=error)
          CALL metadyn_read(meta_env, globenv, metadyn_section, error=error)
          CALL force_env_set(force_env,meta_env=meta_env,error=error)
          CALL meta_env_release(meta_env,error=error)
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
               CPSourceFileRef,error,failure)
          SELECT CASE(globenv%run_type_id)
          CASE (energy_run)
             CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,error=error)
          CASE(energy_force_run)
             CALL force_env_calc_energy_force(force_env,calc_force=.TRUE. ,error=error)
          END SELECT
       CASE (mol_dyn_run)
          CALL qs_mol_dyn ( force_env, globenv, error=error )
       CASE (geo_opt_run)
          CALL qs_geo_opt(force_env,globenv,error=error)
       CASE (mon_car_run)
          CALL do_mon_car ( force_env, globenv, error=error )
       CASE (hmc_run, ldmc_run)
          CALL harris_mc_ld(force_env, globenv, error=error)
       CASE (electronic_spectra_run)
          IF(globenv%prog_name_id /= do_qs) &
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Electron spectra available only with Quickstep program "//&
               CPSourceFileRef,error,failure)
          CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
               error=error)
          CALL tddfpt_calculation(force_env%qs_env, globenv)
       CASE (bsse_run)
          CALL do_bsse_calculation(force_env, globenv, error=error)
       CASE (linear_response_run)
          IF(globenv%prog_name_id /= do_qs) CALL cp_assert(.FALSE.,cp_failure_level,&
               cp_assertion_failed,routineP,&
               "Property calculations by Linear Response only within the QS program "//&
               CPSourceFileRef,error,failure)
          ! The Ground State is needed, it can be read from Restart
          CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
               error=error)
          CALL linres_calculation(force_env%qs_env, globenv, error=error)
       CASE (qs_debug_run)
          IF( (globenv%prog_name_id == do_qs) .OR. &
              (globenv%prog_name_id == do_qmmm) ) THEN
            CALL qs_debug_energy_and_forces(force_env, globenv, error=error)
          ELSE IF(globenv%prog_name_id == do_fist) THEN
            CALL fist_debug_control(force_env, globenv, error=error)
          ELSE
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Debug run available only with QS, FIST, and QMMM program "//&
            CPSourceFileRef,error,failure)
          END IF
       END SELECT
       CALL force_env_release(force_env, error=error)
       
    END IF
    CALL cp_error_dealloc_ref(error)
 
    I_was_ionode=globenv%ionode
    CALL cp2k_finalize(globenv)
    CPPostconditionNoFail(globenv%ref_count==1,cp_failure_level,routineP,error)
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
    TYPE(global_environment_type), &
      pointer                             :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_run', &
      routineP = moduleN//':'//routineN


    LOGICAL                                  :: failure,run_OK,found
    TYPE(farming_env_type), POINTER          :: farming_env
    INTEGER                                  :: stat,ierr,handle,output_unit,&
         length,i,n_rep_val
    CHARACTER(LEN=default_path_length)       :: output_file
    INTEGER, DIMENSION(:), POINTER           :: group_distribution
    INTEGER                                  :: ngroups
    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=default_string_length)     :: keyword, str
    character(len=5)                         :: label
    INTEGER                                  :: new_group,new_rank,new_size
    INTEGER :: n_jobs_to_run, i_job_to_restart, ijob_start, ijob_end,ijob, iunit
    TYPE(cp_parser_type), pointer :: my_parser
    TYPE(section_type), pointer :: g_section
    TYPE(section_vals_type), pointer :: g_data
    TYPE(cp_unit_set_type), pointer :: default_units

    failure=.FALSE.
    CALL timeset(routineN,'I','',handle)
    NULLIFY(my_parser,g_section,g_data,default_units)

    IF (globenv%ionode) write(globenv%scr,FMT="(T2,A)") "FARMING| Hi, welcome on this farm!"

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
        write(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Number of created MPI groups:",ngroups
        write(globenv%scr,FMT="(T2,A)",ADVANCE="NO") "FARMING| Task to group correspondence:"
        DO i=0,globenv%num_pe-1
           IF (MODULO(i,4)==0) write(globenv%scr,*)
           write(globenv%scr,FMT='(A3,I4,A3,I4,A1)',ADVANCE="NO") "  (",i," : ",group_distribution(i),")"
        END DO
        write(globenv%scr,*)
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
       write(iunit,*) i_job_to_restart
       CLOSE(UNIT=iunit)
    ENDIF

    ! this is the job range to be executed. Will be treated cycle
    ijob_start=farming_env%restart_n
    ijob_end=ijob_start+n_jobs_to_run-1

    IF (globenv%ionode) THEN
      IF (ijob_end-ijob_start<0) THEN
        write(globenv%scr,FMT="(T2,A)") "FARMING| --- WARNING --- NO JOBS NEED EXECUTION ? "
        write(globenv%scr,FMT="(T2,A)") "FARMING| is the cycle keyword required ?"
        write(globenv%scr,FMT="(T2,A)") "FARMING| or is a stray RESTARTFARM present ?"
      ELSE
        write(globenv%scr,FMT="(T2,A)") "FARMING| List of jobs : "
        DO ijob=ijob_start,ijob_end
         i=MODULO(ijob-1,farming_env%njobs)+1
         write(globenv%scr,FMT=*) "Job: ",i," Dir: ",TRIM(farming_env%Job_cwd(i))," Input: ", &
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
                  write(output_file,'(A12,I5.5)') "FARMING_OUT_",i
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
          
              IF (globenv%ionode) write(globenv%scr,FMT="(T2,A,I5.5,A)") " Running Job ",i, &
                        " output to "//TRIM(farming_env%Job_cwd(i))//"/"//TRIM(output_file)
              run_OK=cp2k_run(TRIM(farming_env%Job_input(i)),output_unit,new_group)
              
              IF (new_rank == 0) THEN
                  CALL close_file(unit_number=output_unit)
              ELSE
              ENDIF

          ELSE ! problem if just one node does not change...
              IF (globenv%ionode) write(globenv%scr,FMT="(T2,A,A)") &
                     " ... Failed to change dir to "//TRIM(farming_env%Job_cwd(i))
          ENDIF

          ! change to the original working directory
          CALL m_chdir(TRIM(farming_env%cwd),ierr)
          IF (ierr/=0) THEN
              IF (globenv%ionode) write(globenv%scr,FMT="(T2,A,A)") &
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
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='write_cp2k_html_manual',&
         routineP=moduleN//':'//routineN
    integer :: unit_nr
    TYPE(section_type), POINTER :: root_section

    failure=.false.

    IF (.NOT. failure) THEN
       nullify(root_section)
       CALL create_cp2k_root_section(root_section,error=error)
       ! remove the default keyword that ignores things outside the section
       CALL keyword_release(root_section%keywords(0)%keyword,error=error)

       CALL section_describe_html(root_section,"InputReference",6,error=error)

       CALL open_file(unit_number=unit_nr,file_name="index.html",&
            file_action="WRITE", file_status="REPLACE")
       write(unit_nr,*) "<HTML><BODY><HEAD><TITLE>The cp2k input structure</TITLE>"
       write(unit_nr,*) "<H1>CP2K input reference</H1>"
       write(unit_nr,*) "<H2>Version information</H2>"
       write(unit_nr,*) "This html manual can be generated automatically from a given cp2k executable "// &
                        "using the --html-manual command line option. The manual describes exactly that version of the code. "//&
                        "This version has been generated using a cp2k version compiled on "//TRIM(compile_date)//". "//&
                        "The latest CVS entry found is "//TRIM(compile_lastcvs)//"."
       write(unit_nr,*) "<H2>Input structure</H2>"
       write(unit_nr,*) "All sections that can be part of a cp2k input file are shown with their allowed nestings. "
       write(unit_nr,*) "A description of each section, and a list of keywords can be obtained clicking on the links. "
       write(unit_nr,*) "<BR><UL>"

       CALL section_describe_index_html(root_section,"InputReference",unit_nr,error=error)

       write(unit_nr,*) '</UL><BR><hr>Back to the <A HREF="../index.html">CP2K homepage</A>'
       write(unit_nr,*) "</BODY></HTML>"

       call close_file(unit_nr)
       CALL section_release(root_section,error=error)
    END IF
  END SUBROUTINE write_cp2k_html_manual

! *****************************************************************************

  SUBROUTINE write_xml_file(error)

    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "write_xml_file",&
                                   routineP=moduleN//":"//routineN

    TYPE(section_type), POINTER :: root_section
    INTEGER                     :: i,unit_number
    LOGICAL                     :: failure

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

END MODULE cp2k_runs
!!*****
