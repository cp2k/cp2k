!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

MODULE cp2k_runs

  USE atom_module,                     ONLY: atomic_code
  USE environment,                     ONLY: cp2k_finalize,&
                                             cp2k_init
  USE fist_main,                       ONLY: fist_create_force_env
  USE global_types,                    ONLY: global_environment_type
  USE kg_main,                         ONLY: kg_create_force_env
  USE kinds,                           ONLY: default_string_length, &
                                             default_path_length,&
                                             dbl
  USE library_tests,                   ONLY: lib_test
  USE qs_main,                         ONLY: quickstep_create_force_env
  !  USE qs_ep_methods,                   ONLY: ep_create_force_env
  USE timings,                         ONLY: trace_debug
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE farming_methods,                 ONLY: farming_parse_input
  USE farming_types,                   ONLY: farming_env_type, &
                                             init_farming_env, &
                                             deallocate_farming_env
  USE force_env_methods,               ONLY: force_env_calculate_energy
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_release,&
                                             force_env_p_type,&
                                             force_env_set
  USE cp_log_handling,                 ONLY: cp_logger_type, &
                                             cp_failure_level, &
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_error_type, &
                                             cp_a_l, &
                                             cp_assert, &
                                             cp_error_init, &
                                             cp_error_dealloc_ref, &
                                             cp_assertion_failed
  USE geo_opt,                         ONLY: qs_geo_opt
  USE md_run,                          ONLY: qs_mol_dyn 
  USE machine,                         ONLY: print_memory,m_chdir,m_getcwd
  USE message_passing,                 ONLY: mp_sync, &
                                             mp_environ, &
                                             mp_comm_split, &
                                             mp_comm_free
  USE qs_mc,                           ONLY: qs_mon_car
  USE qs_parser,                       ONLY: open_file,close_file
  USE parser,                          ONLY: parser_init, parser_end,&
       get_next, test_next, search_label, read_line
  USE string_utilities,                ONLY: uppercase
  USE trajana_methods,                 ONLY: trajana
  USE timings,                         ONLY: timeset, &
                                             timestop
  USE metadynamics_types, ONLY: meta_env_type, meta_env_release
  USE metadynamics, only: metadyn_read
  USE qmmm_init, only: qmmm_create_force_env
  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: cp2k_run
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
!!     mpi_comm should be a valid communicator
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
!!     mpi_comm        : mpi_comm that can be used for the calling group
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
    INTEGER         , INTENT(IN) :: output_unit
    INTEGER         , INTENT(IN) :: mpi_comm

    CHARACTER(LEN=*), PARAMETER :: routineN = "cp2k_run",&
                                   routineP=moduleN//':'//routineN  

    TYPE(global_environment_type), pointer :: globenv, para_env_p ( : )
    TYPE(cp_error_type) :: error
    TYPE(force_env_p_type), pointer :: force_env_p ( : )
    TYPE(force_env_type), pointer :: force_env
    TYPE(meta_env_type), pointer :: meta_env
    INTEGER :: handle, stat, i,ierr
    logical :: failure

    ! -----------------------------------------------------------------------------
    cp2k_run_OK = .FALSE.
    failure=.false.
    NULLIFY(globenv, force_env, force_env_p, para_env_p )
    ALLOCATE(globenv,stat=stat)
    IF (stat/=0) STOP "cp2k:cp2k, glob_env allocation"

    ! init globenv with the data
    globenv%scr=output_unit
    globenv%group=mpi_comm
    globenv%input_file_name=input_file_name

    ! generate all dependent info in globenv
    CALL cp2k_init ( globenv )
    CALL cp_error_init(error)


    ! initialize force_env
    SELECT CASE (globenv%program_name)
    CASE ("ATOM")
       CALL atomic_code ( globenv )
    CASE ("FARMING")
       CALL farming_run ( globenv, error )
    CASE ("FIST")
      CALL fist_create_force_env ( force_env, globenv,  error = error )
    CASE ("KG")
      CALL kg_create_force_env ( force_env, globenv, error=error)
    CASE ("QS")
      CALL quickstep_create_force_env(force_env, globenv, error=error)
    CASE ("TRAJANA")
       CALL trajana(globenv)
    CASE ("QMMM")
       CALL qmmm_create_force_env(force_env, globenv, error=error)
    CASE ("EP")
       !     CALL ep_create_force_env(force_env,  globenv, error=error)
    CASE ("TEST")
         CALL lib_test(globenv)
    CASE default
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Invalid run TYPE <"//TRIM(globenv%run_type)//"> was specified, "//&
            CPSourceFileRef,&
          error,failure)
    END SELECT
    
    IF (ASSOCIATED(force_env)) THEN
       NULLIFY(meta_env)
       CALL metadyn_read(meta_env, globenv, error=error)
       CALL force_env_set(force_env,meta_env=meta_env,error=error)
       CALL meta_env_release(meta_env,error=error)
    END IF

    ! perform requested calculation
    IF (globenv%run_type/="NONE") THEN 
       !FM this test should be reduced to the first one
  
       IF (globenv%program_name=="QMMM".OR.globenv%program_name=="QS"&
            .OR.globenv%program_name=="EP") THEN
         SELECT CASE (globenv%run_type)
         CASE ("GEOMETRY OPTIMIZATION")
            CALL qs_geo_opt(force_env,globenv,error=error)
         CASE ("MONTE CARLO")
           CALL qs_mon_car ( force_env, globenv, error=error )
         CASE ("MOLECULAR DYNAMICS")
            CALL qs_mol_dyn ( force_env, globenv, error=error )
         CASE ("WAVEFUNCTION OPTIMIZATION","ENERGY CALCULATION")
            CALL force_env_calculate_energy(force_env,error=error)
         CASE DEFAULT
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                   "Invalid run TYPE <"//TRIM(globenv%run_type)//"> was specified, "//&
               CPSourceFileRef,&
               error,failure)
         END SELECT
         CALL force_env_release(force_env, error=error)
       ELSE IF ( globenv % program_name == "FIST" ) THEN
          SELECT CASE (globenv%run_type)
          CASE ("MOLECULAR DYNAMICS")
             CALL qs_mol_dyn ( force_env, globenv, error=error )
          CASE ( "GEOMETRY OPTIMIZATION" )
             CALL qs_geo_opt(force_env,globenv,error=error)
          CASE ( "MONTE CARLO" )
             CALL qs_mon_car ( force_env, globenv, error=error )
!         CASE ( "DEBUG" )
!            CALL fist_debug ( force_env, error = error )
          END SELECT
          CALL force_env_release(force_env, error=error)
#if defined(__HMD)
         IF (globenv%ionode) CALL _dump_allocated(%VAL(0))
#endif
       ELSE IF(globenv%program_name=="KG") THEN
         SELECT CASE (globenv%run_type)
         CASE ("MOLECULAR DYNAMICS")
            CALL qs_mol_dyn ( force_env, globenv, error=error )
         CASE ("WAVEFUNCTION OPTIMIZATION","ENERGY CALCULATION")
          CALL force_env_calculate_energy(force_env,error=error)
!        CASE ("DEBUG")
!           CALL kg_debug(force_env,error=error)
         END SELECT
         CALL force_env_release(force_env, error=error)
       END IF
    END IF

    CALL cp_error_dealloc_ref(error)

    CALL cp2k_finalize(globenv)
    DEALLOCATE(globenv,stat=stat)
    IF (stat/=0) STOP "cp2k:cp2k, glob_env deallocation"
#if defined(__HMD)
    IF (globenv%ionode) CALL _dump_allocated(%VAL(0))
#endif
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
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_run', &
      routineP = moduleN//':'//routineN


    LOGICAL                                  :: failure,run_OK
    TYPE(farming_env_type), POINTER          :: farming_env
    INTEGER                                  :: stat,i,ierr,output_unit,handle,&
         length
    CHARACTER(LEN=default_path_length)       :: output_file
    INTEGER, DIMENSION(:), POINTER           :: group_distribution
    INTEGER                                  :: ngroups
    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=default_string_length)     :: keyword, str
    character(len=5)                         :: label
    INTEGER                                  :: new_group,new_rank,new_size

    failure=.FALSE.
    CALL timeset(routineN,'I','',handle)

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
               n_subgroups=farming_env%ngroup_wish,group_partition=farming_env%group_partition)
    ELSE
       CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (globenv%ionode) THEN
        write(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Number of created MPI groups:",ngroups
        write(globenv%scr,FMT="(T2,A,T71,I10)",ADVANCE="NO") "FARMING| Task to group correspondence:"
        DO i=0,globenv%num_pe-1
           IF (MODULO(i,4)==0) write(globenv%scr,*)
           write(globenv%scr,FMT='(A3,I4,A3,I4,A1)',ADVANCE="NO") "  (",i," : ",group_distribution(i),")"
        END DO
        write(globenv%scr,*)
    ENDIF
    IF (globenv%ionode) THEN
        write(globenv%scr,FMT="(T2,A)") "FARMING| List of jobs : "
        DO i=1,farming_env%njobs
         write(globenv%scr,FMT=*) "Job: ",i," Dir: ",TRIM(farming_env%Job_cwd(i))," Input: ", &
           TRIM(farming_env%Job_input(i))," MPI group:", MODULO(i-1,ngroups)
        ENDDO
    ENDIF

    DO i=1,farming_env%njobs
       ! this farms out the jobs
       IF (MODULO(i-1,ngroups)==group_distribution(globenv%mepos)) THEN
          ! change to the new working directory
          CALL m_chdir(TRIM(farming_env%Job_cwd(i)),ierr)

          ! generate a fresh call to cp2k_run
          IF (ierr==0) THEN
              CALL mp_environ(new_size,new_rank,new_group)
              IF (new_rank == 0) THEN 

                 ! generate the output file
                  write(output_file,'(A12,I5.5)') "FARMING_OUT_",i
                  CALL parser_init(file_name=TRIM(farming_env%Job_input(i)))
                  label="&CP2K"
                  CALL search_label(label,ierr,ignore_case=.true.)
                  IF (ierr==0) THEN
                     DO
                        CALL read_line()
                        test_result = test_next()
                        IF (test_result == "X") exit
                        length=default_string_length
                        CALL get_next(keyword,length)
                        call uppercase(keyword)
                        SELECT CASE(keyword)
                        CASE("PROJECT")
                           length=default_string_length
                           CALL get_next(str,length)
                           output_file=TRIM(str)//".out"
                        CASE("OUTPUT_FILE_NAME","OUTPUT_FILE")
                           length=default_string_length
                           CALL get_next(str,length)
                           output_file=TRIM(str)
                           EXIT
                        END SELECT
                     END DO
                  END IF
                  CALL parser_end()

                  CALL open_file(file_name=TRIM(output_file),&
                        file_action="WRITE",&
                        file_status="UNKNOWN",&
                        file_position="APPEND",&
                        unit_number=output_unit)
              ELSE
                output_unit=-1 ! this should never be used by cp2k_run
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

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE farming_run

END MODULE cp2k_runs
!!*****
