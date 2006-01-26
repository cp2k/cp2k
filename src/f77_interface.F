!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/f77_interface [1.0] *
!!
!!   NAME
!!     f77_interface
!!
!!   FUNCTION
!!     interface to use cp2k as library
!!
!!   NOTES
!!     useful additions for the future would be:
!!     - string(path) based set/get of simple values (to change the new
!!       input during the run and extract more data (energy types for example).
!!     - set/get of a subset of atoms
!!
!!   AUTHOR
!!     fawzi & Johanna
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [fawzi]
!!     11.2004 parallel version [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE f77_interface
  USE cp2k_runs,                       ONLY: cp2k_run
  USE cp_files,                        ONLY: open_file, close_file 
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE environment,                     ONLY: cp2k_finalize,&
                                             cp2k_init,&
                                             cp2k_read,&
                                             cp2k_setup
  USE mathconstants,                   ONLY: init_mathcon
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release,&
                                             force_env_type,&
                                             force_env_retain
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create
  USE input_cp2k_check,                ONLY: check_cp2k_input
  USE input_cp2k,                      ONLY: create_cp2k_root_section
  USE input_constants,                 ONLY: do_qs,&
                                             do_fist,&
                                             do_kg,&
                                             do_qmmm,&
                                             do_atom,&
                                             do_farming,&
                                             do_test,&
                                             run_type_names
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: section_release,&
                                             section_type,&
                                             section_vals_create,&
                                             section_vals_release,&
                                             section_vals_type,&
                                             section_vals_write
  USE kg_main,                         ONLY: kg_create_force_env
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: MPI_COMM_WORLD,&
                                             mp_max,&
                                             mp_world_finalize,&
                                             mp_world_init,&
                                             add_mp_perf_env,&
                                             rm_mp_perf_env,&
                                             get_mp_perf_env,&
                                             mp_perf_env_type,&
                                             mp_perf_env_release,&
                                             mp_perf_env_retain
  USE particle_list_types,             ONLY: particle_list_type
  USE qmmm_init,                       ONLY: qmmm_create_force_env
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE timings,                         ONLY: timeset,&
                                             timestop,&
                                             add_timer_env,&
                                             rm_timer_env,&
                                             timer_env_type,&
                                             get_timer_env,&
                                             timer_env_retain,&
                                             timer_env_release
  USE cp_units, ONLY: cp_unit_set_type, cp_unit_set_create, cp_unit_set_release
  USE physcon,                         ONLY: init_physcon
  USE cp2k_units,                      ONLY: init_cp2k_units,&
                                             rm_cp2k_units
#include "cp_common_uses.h"
  IMPLICIT NONE
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='f77_interface'
  INTEGER, parameter :: PATH_LENGTH=200

  TYPE f_env_p_type
     TYPE(f_env_type), pointer :: f_env
  END TYPE f_env_p_type

  TYPE f_env_type
     integer :: id_nr
     TYPE(force_env_type), pointer :: force_env
     TYPE(cp_error_type) :: error
     TYPE(timer_env_type), pointer :: timer_env
     TYPE(mp_perf_env_type), pointer :: mp_perf_env
  END TYPE f_env_type

  TYPE(f_env_p_type), DIMENSION(:), POINTER, save :: f_envs
  TYPE(cp_para_env_type), POINTER, save :: default_para_env
  LOGICAL, SAVE :: module_initialized=.FALSE.
  INTEGER, SAVE :: last_f_env_id=0, n_f_envs=0

  PUBLIC :: default_para_env
  PUBLIC :: init_cp2k, finalize_cp2k
  PUBLIC :: create_force_env, destroy_force_env, set_pos, get_pos,&
       get_force, calc_energy_force, run_input, get_energy, &
       calc_energy, calc_force, check_input, get_natom
  !,set_force, set_vel, get_vel 
  !!***
  !****************************************************************************
contains

  !!****f* f77_interface/get_pos_of_env [1.0] *
  !!
  !!   NAME
  !!     get_pos_of_env
  !!
  !!   FUNCTION
  !!     returns the position of the force env corresponding to the given id
  !!
  !!   NOTES
  !!     private utility function
  !!
  !!   INPUTS
  !!     - env_id: the id of the requested environment
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  FUNCTION get_pos_of_env(env_id) RESULT(res)
    INTEGER, intent(in)                      :: env_id
    INTEGER                                  :: res

    INTEGER                                  :: env_pos, isub

    env_pos=-1
    DO isub=1,n_f_envs
       IF (f_envs(isub)%f_env%id_nr==env_id) THEN
          env_pos=isub
       END IF
    END DO
    res=env_pos
  END FUNCTION get_pos_of_env
  !****************************************************************************

  !!****f* f77_interface/init_cp2k *
  !!
  !!   NAME
  !!     init_cp2k
  !!
  !!   FUNCTION
  !!     initializes cp2k, needs to be called once before using any of the
  !!     other functions when using cp2k as library
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - init_mpi: if the mpi environment should be initialized
  !!     - ierr: returns a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE init_cp2k(init_mpi,ierr)
    LOGICAL, intent(in)                      :: init_mpi
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "init_cp2k", &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: default_output_unit = 6

    INTEGER                                  :: mpi_comm_default, stat, &
                                                unit_nr
    TYPE(cp_logger_type), pointer            :: logger

    IF (.NOT. module_initialized) THEN
       IF (init_mpi) THEN
          ! get the default system wide communicator
          CALL mp_world_init(mpi_comm_default)
       ELSE
          mpi_comm_default=MPI_COMM_WORLD
       END IF

       call add_mp_perf_env()
       call add_timer_env()
       nullify(default_para_env)
       CALL cp_para_env_create(default_para_env, group=mpi_comm_default, &
            owns_group=.FALSE.)
       IF (default_para_env%source==default_para_env%mepos) THEN
          unit_nr=6
       ELSE
          unit_nr=-1
       END IF
       nullify(logger)
       CALL cp_logger_create(logger,para_env=default_para_env,&
            default_global_unit_nr=unit_nr, &
            close_global_unit_on_dealloc=.FALSE.)
       CALL cp_add_default_logger(logger)
       call cp_logger_release(logger)

       ALLOCATE(f_envs(0),stat=stat)
       CPPostconditionNoErr(stat==0,cp_fatal_level,routineP)
       module_initialized=.true.
       ierr=0

       !   *** Initialize mathematical constants ***
       
       CALL init_mathcon()
       
       !   *** Initialize physical constants ***
       
       CALL init_physcon()
       CALL init_cp2k_units()
    ELSE
       ierr=cp_failure_level
    END IF

  END SUBROUTINE init_cp2k
  !****************************************************************************

  !!****f* f77_interface/finalize_cp2k *
  !!
  !!   NAME
  !!     finalize_cp2k
  !!
  !!   FUNCTION
  !!     cleanup after you have finished using this interface
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - finalize_mpi: if the mpi environment should be finalized
  !!     - ierr: returns a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE finalize_cp2k(finalize_mpi,ierr)
    LOGICAL, intent(in)                      :: finalize_mpi
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "finalize_mpi", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: ienv, stat
    TYPE(cp_error_type)                      :: error

    IF (.NOT.module_initialized) THEN
       ierr=cp_failure_level
    ELSE
       CALL cp_error_init(error, stop_level=cp_fatal_level)
       DO ienv=n_f_envs,1,-1
          CALL destroy_force_env(f_envs(ienv)%f_env%id_nr,ierr=ierr)
          CPAssertNoFail(ierr==0,cp_warning_level,routineP,error)
       END DO
       DEALLOCATE(f_envs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CALL cp_para_env_release(default_para_env,error)
       ierr=cp_error_get_level(error)
       CALL cp_rm_default_logger()
       CALL rm_cp2k_units()
       CALL rm_timer_env()
       CALL rm_mp_perf_env()
       IF (finalize_mpi) THEN
          CALL mp_world_finalize()
       END IF
    END IF
  END SUBROUTINE finalize_cp2k
  !****************************************************************************

!!****f* f77_interface/f_env_dealloc *
!!
!!   NAME
!!     f_env_dealloc
!!
!!   FUNCTION
!!     deallocates a f_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - f_env: the f_env to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE f_env_dealloc(f_env,error)
  TYPE(f_env_type), pointer :: f_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='f_env_dealloc',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(ASSOCIATED(f_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL force_env_release(f_env%force_env,error=error)
     CALL cp_error_dealloc_ref(f_env%error,error=error)
     CALL timer_env_release(f_env%timer_env)
     call mp_perf_env_release(f_env%mp_perf_env)
  END IF
END SUBROUTINE f_env_dealloc
!***************************************************************************

!!****f* f77_interface/f_env_create *
!!
!!   NAME
!!     f_env_create
!!
!!   FUNCTION
!!     createates a f_env
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - f_env: the f_env to createate
!!     - template_error: the error to be used as template for the one stored
!!     - force_env: the force_environment to be stored
!!     - timer_env: the timer env to be stored
!!     - mp_perf_env: the mp performance environement to be stored
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE f_env_create(f_env,template_error,force_env, timer_env,mp_perf_env,&
     id_nr,error)
  TYPE(f_env_type), pointer :: f_env
  TYPE(force_env_type), pointer :: force_env
  TYPE(timer_env_type), pointer :: timer_env
  TYPE(mp_perf_env_type), pointer :: mp_perf_env
  INTEGER, intent(in) :: id_nr
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error, template_error
  
  logical :: failure
  character(len=*), parameter :: routineN='f_env_create',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  
  ALLOCATE(f_env,stat=stat)
  CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     f_env%force_env => force_env
     CALL force_env_retain(f_env%force_env,error=error)
     CALL cp_error_init(f_env%error,template_error=template_error)
     f_env%timer_env => timer_env
     CALL timer_env_retain(f_env%timer_env)
     f_env%mp_perf_env => mp_perf_env
     call mp_perf_env_retain(f_env%mp_perf_env)
     f_env%id_nr=id_nr
  END IF
END SUBROUTINE f_env_create
!***************************************************************************

!!****f* f77_interface/f_env_add_defaults *
!!
!!   NAME
!!     f_env_add_defaults
!!
!!   FUNCTION
!!     adds the default environments of the f_env to the stack of the 
!!     defaults, and returns a new error and sets failure to true if
!!     something went wrong
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - f_env_id: the f_env from where to take the defaults
!!     - f_env: will contain the f_env corresponding to f_env_id
!!     - new_error: an error that can be used for the given f_env
!!     - failure: will be set to true if something went wrong
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE f_env_add_defaults(f_env_id,f_env,new_error, failure)
  INTEGER, intent(in) :: f_env_id
  TYPE(f_env_type), pointer :: f_env
  type(cp_error_type), intent(out) :: new_error
  logical :: failure
  
  logical :: my_failure
  character(len=*), parameter :: routineN='f_env_add_defaults',&
        routineP=moduleN//':'//routineN
  TYPE(cp_logger_type), pointer :: logger
  integer :: f_env_pos

  my_failure=.false.

  NULLIFY(f_env)
  f_env_pos = get_pos_of_env(f_env_id)
  IF (f_env_pos<1) THEN
     call cp_error_init(new_error)
     call cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "invalid env_id "//cp_to_string(f_env_id),error=new_error,&
          failure=failure)
  ELSE
     f_env => f_envs(f_env_pos)%f_env
     CALL cp_error_init(new_error,template_error=f_env%error)
     logger => cp_error_get_logger(f_env%error)
     CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,new_error,my_failure)
     IF (.NOT.my_failure) THEN
        CALL add_mp_perf_env(f_env%mp_perf_env)
        CALL add_timer_env(f_env%timer_env)
        CALL cp_add_default_logger(logger)
     END IF
     CALL cp_error_check(new_error,failure)
  END IF
END SUBROUTINE f_env_add_defaults
!***************************************************************************

!!****f* f77_interface/f_env_rm_defaults *
!!
!!   NAME
!!     f_env_rm_defaults
!!
!!   FUNCTION
!!     removes the default environments of the f_env to the stack of the 
!!     defaults, and sets ierr accordingly to the failuers stored in error
!!     It also releases the error
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - f_env: the f_env from where to take the defaults
!!     - error: the error to be checked (will be released)
!!     - ierr: variable that will be set to a number different from 0 if
!!       error contains an error (otherwise it will be set to 0)
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE f_env_rm_defaults(f_env,error,ierr)
  TYPE(f_env_type), pointer :: f_env
  type(cp_error_type), intent(inout) :: error
  INTEGER, INTENT(out), optional :: ierr

  logical :: failure
  character(len=*), parameter :: routineN='f_env_rm_defaults',&
       routineP=moduleN//':'//routineN
  TYPE(cp_logger_type), POINTER :: logger,d_logger
  TYPE(timer_env_type), pointer :: d_timer_env
  TYPE(mp_perf_env_type), pointer :: d_mp_perf_env

  failure=.false.

  IF (ASSOCIATED(f_env)) THEN
     logger => cp_error_get_logger(f_env%error)
     d_logger => cp_get_default_logger()
     d_timer_env => get_timer_env()
     d_mp_perf_env => get_mp_perf_env()
     CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(d_logger),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(d_timer_env),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(d_mp_perf_env),cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        CPPrecondition(logger%id_nr==d_logger%id_nr,cp_failure_level,routineP,error,failure)
        CPPrecondition(d_timer_env%id_nr==f_env%timer_env%id_nr,cp_failure_level,routineP,error,failure)
        CPPrecondition(d_mp_perf_env%id_nr==f_env%mp_perf_env%id_nr,cp_failure_level,routineP,error,failure)
     END IF
     if (present(ierr)) ierr=cp_error_get_level(error)
     call cp_error_dealloc_ref(error)
     IF (.NOT.failure) THEN
        CALL cp_rm_default_logger()
        CALL rm_timer_env()
        CALL rm_mp_perf_env()
     END IF
  ELSE
     if (present(ierr)) ierr=cp_error_get_level(error)
     call cp_error_dealloc_ref(error)
  END IF
END SUBROUTINE f_env_rm_defaults
!***************************************************************************

  !!****f* f77_interface/run_input [1.0] *
  !!
  !!   NAME
  !!     run_input
  !!
  !!   FUNCTION
  !!     runs the given input
  !!
  !!   NOTES
  !!     -
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
    INTEGER, intent(out)                     :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure, success
    TYPE(cp_para_env_type), pointer          :: para_env

    failure=.false.

    IF (PRESENT(mpi_comm)) THEN
       NULLIFY(para_env)
       CALL cp_para_env_create(para_env, group=mpi_comm, owns_group=.false.)
    ELSE
       para_env => default_para_env
       call cp_para_env_retain(para_env)
    END IF
    IF (para_env%mepos==para_env%source) THEN
       IF (output_file_path=="__STD_OUT__") THEN
          unit_nr=6
       ELSE
          call open_file(file_name=output_file_path,file_status="UNKNOWN",&
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
    call cp_para_env_release(para_env)
  END SUBROUTINE run_input
  !****************************************************************************

  !!****f* f77_interface/create_force_env [1.0] *
  !!
  !!   NAME
  !!     create_force_env
  !!
  !!   FUNCTION
  !!     
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - new_env_id: will contain the id of the newly created environment
  !!     - input_file_path: the path of the input file
  !!     - output_file_path: path of the output file (to which it is appended)
  !!       if it is "__STD_OUT__" the unit 6 is used
  !!     - mpi_comm: the mpi communicator to be used for this environment
  !!       it will not be freed when you get rid of the force_env
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE create_force_env(new_env_id,input_file_path,&
       output_file_path,mpi_comm,ierr)
    INTEGER, INTENT(out)                     :: new_env_id
    CHARACTER(len=*), INTENT(in)             :: input_file_path, &
                                                output_file_path
    INTEGER, INTENT(out), optional                     :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat, unit_nr
    INTEGER, DIMENSION(:), pointer           :: f_env_ids_old
    logical                                  :: failure
    type(cp_error_type)                      :: error
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(cp_para_env_type), pointer          :: para_env
    TYPE(cp_parser_type), pointer            :: cpparser
    TYPE(f_env_p_type), DIMENSION(:), &
      pointer                                :: f_envs_old
    TYPE(force_env_type), pointer :: force_env
    TYPE(timer_env_type), pointer :: timer_env
    TYPE(mp_perf_env_type), pointer :: mp_perf_env
    TYPE(global_environment_type), pointer   :: globenv
    TYPE(section_type), pointer              :: input_structure
    TYPE(section_vals_type), pointer         :: input_file
    TYPE(cp_unit_set_type), pointer :: default_units

    failure=.false.

    NULLIFY(para_env, force_env,timer_env,mp_perf_env, globenv)
    IF (PRESENT(mpi_comm)) THEN
       CALL cp_para_env_create(para_env, group=mpi_comm, owns_group=.false.)
    ELSE
       para_env => default_para_env
       call cp_para_env_retain(para_env)
    END IF
    IF (para_env%mepos==para_env%source) THEN
       IF (output_file_path=="__STD_OUT__") THEN
          unit_nr=6
       ELSE
          call open_file(file_name=output_file_path,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=unit_nr)
       END IF
    ELSE
       unit_nr=-1
    END IF
    CALL globenv_create(globenv, para_env=para_env,&
         output_unit=unit_nr, input_file_name=input_file_path,&
         close_output_unit=(unit_nr>0 .and. unit_nr/=6))
    call cp2k_init(globenv)
    logger => cp_get_default_logger()
    ! warning this is dangerous, I did not check that all the subfunctions
    ! support it, the program might crash upon error
    CALL cp_error_init(error,stop_level=cp_fatal_level,logger=logger)

    IF (n_f_envs+1>SIZE(f_envs)) THEN
       f_envs_old => f_envs
       ALLOCATE(f_envs(n_f_envs+10), stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       DO i=1,n_f_envs
          f_envs(i)%f_env => f_envs_old(i)%f_env
       END DO
       DO i=n_f_envs+1,size(f_envs)
          NULLIFY(f_envs(i)%f_env)
       END DO
       DEALLOCATE(f_envs_old,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF

    ! parse the input
    NULLIFY(input_structure, default_units)
    CALL create_cp2k_root_section(input_structure,error=error)
    NULLIFY(input_file)
    CALL section_vals_create(input_file,input_structure, error=error)
    CALL section_release(input_structure,error=error)
    nullify(cpparser)
    CALL parser_create(cpparser,file_name=input_file_path,para_env=para_env,&
         error=error)
    CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
    CALL section_vals_parse(input_file,cpparser,root_section=.false.,&
         default_units=default_units,error=error)
    CALL cp_unit_set_release(default_units,error=error)
    CALL parser_release(cpparser,error=error)
    IF (globenv%ionode) THEN
       CALL section_vals_write(input_file,unit_nr=globenv%scr,hide_root=.TRUE.,&
            error=error)
    END IF
    globenv%input_file => input_file
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL check_cp2k_input(input_file,para_env=para_env,error=error)
       CALL cp_error_check(error,failure)
    END IF
    NULLIFY(input_file)

    IF (.NOT.failure) THEN

       CALL cp2k_read(globenv)
       CALL cp2k_setup(globenv,error)

       ! initialize force_env
       SELECT CASE (globenv%prog_name_id)
       CASE (do_fist)
          CALL fist_create_force_env ( force_env,&
               globenv,  error = error )
       CASE (do_kg)
          CALL kg_create_force_env ( force_env,&
               globenv, error=error)
       CASE (do_qs)
          CALL quickstep_create_force_env( force_env,&
               globenv, error=error)
       CASE (do_qmmm)
          CALL qmmm_create_force_env( force_env,&
               globenv, error=error)
       CASE default
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Invalid run TYPE <"//TRIM(run_type_names(globenv%run_type_id))//"> was specified, "//&
               CPSourceFileRef,error,failure)
       END SELECT
    END IF
    timer_env => get_timer_env()
    mp_perf_env => get_mp_perf_env()
    last_f_env_id=last_f_env_id+1
    new_env_id=last_f_env_id
    CALL mp_max(last_f_env_id,globenv%group)
    n_f_envs=n_f_envs+1
    CALL f_env_create(f_envs(n_f_envs)%f_env,template_error=error,&
         timer_env=timer_env,mp_perf_env=mp_perf_env,force_env=force_env,&
         id_nr=last_f_env_id,error=error)
    CALL f_env_rm_defaults(f_envs(n_f_envs)%f_env,error=error,ierr=ierr)

  END SUBROUTINE create_force_env
  !****************************************************************************

  !!****f* f77_interface/destroy_force_env *
  !!
  !!   NAME
  !!     destroy_force_env
  !!
  !!   FUNCTION
  !!     deallocates the force_env with the given id
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - env_id: the id of the force_env to remove
  !!     - ierr: will contain a number different from 0 if 
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE destroy_force_env(env_id,ierr)
    INTEGER, intent(in)                      :: env_id
    INTEGER, intent(out)                     :: ierr

    character(len=*), parameter :: routineN = 'destroy_force_env', &
      routineP = moduleN//':'//routineN

    integer                                  :: env_pos, i, stat
    INTEGER, DIMENSION(:), pointer           :: f_env_ids_old
    logical                                  :: failure
    type(cp_error_type)                      :: error
    TYPE(cp_error_type), DIMENSION(:), &
      pointer                                :: errors_old
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(global_environment_type), pointer   :: globenv
    TYPE(f_env_type), pointer :: f_env

    failure=.false.
    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.not.failure) THEN
       env_pos=get_pos_of_env(env_id)

       CALL force_env_get(f_env%force_env,globenv=globenv,&
            error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL cp2k_finalize(globenv,error)
          CALL f_env_dealloc(f_envs(env_pos)%f_env,error=error)
          CALL cp_error_check(error,failure)

          n_f_envs=n_f_envs-1
       END IF
       ierr=cp_error_get_level(error)
       call cp_error_dealloc_ref(error)
       call cp_rm_default_logger()
    ELSE
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE destroy_force_env
  !***************************************************************************

  !!****f* f77_interface/get_natom [1.0] *
  !!
  !!   NAME
  !!     get_natom
  !!
  !!   FUNCTION
  !!     returns the number of atoms in the given force env
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - new_pos: the array with the new positions
  !!     - n_el: number of positions (3*natom) just to check
  !!     - env_id: id of the force_env
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE get_natom(env_id, n_atom, ierr)
    INTEGER, INTENT(in)                      :: env_id
    INTEGER, INTENT(out)                     :: n_atom,ierr

    CHARACTER(len=*), PARAMETER :: routineN = "get_natom", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles
    TYPE(f_env_type), pointer                :: f_env

    failure=.false.

    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    n_atom=0
    IF (.not.failure) THEN
       nullify(subsys)
       CALL force_env_get(f_env%force_env,&
            subsys=subsys,error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             n_atom=n_atom+particles%n_els
          END DO
       END IF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE get_natom
  !****************************************************************************

  !!****f* f77_interface/set_pos [1.0] *
  !!
  !!   NAME
  !!     set_pos
  !!
  !!   FUNCTION
  !!     sets the positions of the atoms
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - new_pos: the array with the new positions
  !!     - n_el: number of positions (3*natom) just to check
  !!     - env_id: id of the force_env
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE set_pos(env_id, new_pos, n_el, ierr)
    INTEGER, INTENT(in)                      :: env_id, n_el
    REAL(kind=dp), DIMENSION(1:n_el)         :: new_pos
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "set_pos", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub, &
                                                my_n_el
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles
    TYPE(f_env_type), pointer                :: f_env

    failure=.false.

    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.not.failure) THEN
       nullify(subsys)
       CALL force_env_get(f_env%force_env,&
            subsys=subsys,error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          my_n_el=0
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             my_n_el=my_n_el+particles%n_els
          END DO
          my_n_el=3*my_n_el
          CALL cp_assert(my_n_el==n_el,cp_failure_level,cp_assertion_failed,&
               routineP,"wrong pos size ("//cp_to_string(n_el)//" vs"//&
               cp_to_string(my_n_el)//") "//&
               CPSourceFileRef,&
               error,failure)
       END IF
       IF (.NOT.failure) THEN
          ii=0
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             DO ip=1,particles%n_els
                DO idir=1,3
                   ii=ii+1
                   particles%els(ip)%r(idir)=new_pos(ii)
                END DO
             END DO
          END DO
       END IF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE set_pos
  !****************************************************************************

  !!****f* f77_interface/get_pos [1.0] *
  !!
  !!   NAME
  !!     get_pos
  !!
  !!   FUNCTION
  !!     gets the positions of the atoms
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - pos: the array where to write the positions
  !!     - n_el: number of positions (3*natom) just to check
  !!     - env_id: id of the force_env
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE get_pos(env_id, pos, n_el, ierr)
    INTEGER, INTENT(in)                      :: env_id, n_el
    REAL(kind=dp), DIMENSION(1:n_el)         :: pos
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "get_pos", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub, &
                                                my_n_el
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles
    TYPE(f_env_type), pointer :: f_env

    failure=.false.

    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.not.failure) THEN
       IF (.NOT.failure) THEN
          CALL force_env_get(f_env%force_env,&
               subsys=subsys,error=error)
          CALL cp_error_check(error,failure)
       END IF
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          my_n_el=0
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             my_n_el=my_n_el+particles%n_els
          END DO
          my_n_el=3*my_n_el
          CALL cp_assert(my_n_el==n_el,cp_failure_level,cp_assertion_failed,&
               routineP,"wrong pos size ("//cp_to_string(n_el)//" vs"//&
               cp_to_string(my_n_el)//&
               CPSourceFileRef,&
               error,failure)
       END IF
       IF (.NOT.failure) THEN
          ii=0
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             DO ip=1,particles%n_els
                DO idir=1,3
                   ii=ii+1
                   pos(ii)=particles%els(ip)%r(idir)
                END DO
             END DO
          END DO
       END IF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE get_pos
  !****************************************************************************

  !!****f* f77_interface/get_force [1.0] *
  !!
  !!   NAME
  !!     get_force
  !!
  !!   FUNCTION
  !!     gets the forces of the atoms
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - force: the array where to write the forces
  !!     - n_el: number of positions (3*natom) just to check
  !!     - env_id: id of the force_env
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE get_force(env_id, force, n_el, ierr)
    INTEGER, INTENT(in)                      :: env_id, n_el
    REAL(kind=dp), DIMENSION(1:n_el)         :: force
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "get_force", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub, &
                                                my_n_el
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles
    TYPE(f_env_type), pointer                :: f_env

    failure=.false.

    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.not.failure) THEN
       IF (.NOT.failure) THEN
          nullify(subsys)
          CALL force_env_get(f_env%force_env,&
               subsys=subsys,error=error)
          CALL cp_error_check(error,failure)
       END IF
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          my_n_el=0
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             my_n_el=my_n_el+particles%n_els
          END DO
          my_n_el=3*my_n_el
          CALL cp_assert(my_n_el==n_el,cp_failure_level,cp_assertion_failed,&
               routineP,"wrong force size ("//cp_to_string(n_el)//" vs"//&
               cp_to_string(my_n_el)//&
               CPSourceFileRef,&
               error,failure)
       END IF
       IF (.NOT.failure) THEN
          ii=0
          DO isub=1,SIZE(subsys)
             particles => subsys(isub)%subsys%particles
             DO ip=1,particles%n_els
                DO idir=1,3
                   ii=ii+1
                   force(ii)=particles%els(ip)%f(idir)
                END DO
             END DO
          END DO
       END IF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE get_force
  !****************************************************************************

  !!****f* f77_interface/calc_energy_force [1.0] *
  !!
  !!   NAME
  !!     calc_energy_force
  !!
  !!   FUNCTION
  !!     updates the energy and the forces of given force_env
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - env_id: id of the force_env that you want to update
  !!     - calc_force: if the forces should be updated, if false the forces
  !!       might be wrong.
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE calc_energy_force(env_id,calc_force,ierr)
    INTEGER, INTENT(in)                      :: env_id
    LOGICAL, intent(in)                      :: calc_force
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "calc_energy_force", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), pointer :: f_env

    failure=.false.
   
    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_calc_energy_force(f_env%force_env,&
            calc_force=calc_force,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE calc_energy_force
  !****************************************************************************

  !!****f* f77_interface/get_energy *
  !!
  !!   NAME
  !!     get_energy
  !!
  !!   FUNCTION
  !!     returns the energy of the last configuration calculated
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - env_id: id of the force_env that you want to update
  !!     - e_pot: the potential energy of the system
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE get_energy(env_id,e_pot,ierr)
    INTEGER, INTENT(in)                      :: env_id
    REAL(kind=dp), intent(out)               :: e_pot
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "get_energy", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), pointer :: f_env

    failure=.false.

    nullify(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_get(f_env%force_env,&
            potential_energy=e_pot, error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE get_energy
  !****************************************************************************

  !!****f* f77_interface/calc_energy *
  !!
  !!   NAME
  !!     calc_energy
  !!
  !!   FUNCTION
  !!     returns the energy of the configuration given by the positions
  !!     passed as argument
  !!
  !!   NOTES
  !!     utility call
  !!
  !!   INPUTS
  !!     - env_id: id of the force_env that you want to update
  !!     - pos: array with the positions
  !!     - n_el: number of elements in pos (3*natom)
  !!     - e_pot: the potential energy of the system
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE calc_energy(env_id,pos,n_el,e_pot,ierr)
    INTEGER, INTENT(in)                      :: env_id
    INTEGER, INTENT(in)                      :: n_el
    REAL(kind=dp), dimension(1:n_el),&
         INTENT(in)                          :: pos
    REAL(kind=dp), intent(out)               :: e_pot
    INTEGER, intent(out)                     :: ierr

    REAL(kind=dp), dimension(1)              :: dummy_f

    CALL calc_force(env_id,pos,n_el,e_pot,dummy_f,0,ierr)
  END SUBROUTINE calc_energy
  !****************************************************************************

  !!****f* f77_interface/calc_force *
  !!
  !!   NAME
  !!     calc_force
  !!
  !!   FUNCTION
  !!     returns the energy of the configuration given by the positions
  !!     passed as argument
  !!
  !!   NOTES
  !!     utility call, but actually it could be a better and more efficient
  !!     interface to connect to other codes if cp2k would be deeply
  !!     refactored
  !!
  !!   INPUTS
  !!     - env_id: id of the force_env that you want to update
  !!     - pos: array with the positions
  !!     - n_el_pos: number of elements in pos (3*natom)
  !!     - e_pot: the potential energy of the system
  !!     - force: array that will contain the forces
  !!     - n_el_force: number of elements in force (3*natom). If 0 the 
  !!       forces are not calculated
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE calc_force(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)
    INTEGER, INTENT(in)                      :: env_id, n_el_pos
    REAL(kind=dp), dimension(1:n_el_pos), &
      INTENT(in)                             :: pos
    REAL(kind=dp), intent(out)               :: e_pot
    INTEGER, INTENT(in)                      :: n_el_force
    REAL(kind=dp), dimension(1:n_el_force), &
      INTENT(inout)                          :: force
    INTEGER, intent(out)                     :: ierr

    logical                                  :: calc_f

    calc_f=(n_el_force/=0)
    CALL set_pos(env_id,pos,n_el_pos,ierr)
    IF (ierr==0) CALL calc_energy_force(env_id,calc_f,ierr)
    IF (ierr==0) CALL get_energy(env_id,e_pot,ierr)
    IF (calc_f.AND.ierr==0) CALL get_force(env_id,force,n_el_force,ierr)
  END SUBROUTINE calc_force
  !****************************************************************************

  !!****f* f77_interface/check_input *
  !!
  !!   NAME
  !!     check_input
  !!
  !!   FUNCTION
  !!     performs a check of the input
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - input_file_path: the path of the input file to check
  !!     - output_file_path: path of the output file (to which it is appended)
  !!       if it is "__STD_OUT__" the unit 6 is used
  !!     - echo_input: if the parsed input should be written out with all the
  !!       defaults made explicit
  !!     - mpi_comm: the mpi communicator (if not given it uses the default
  !!       one)
  !!     - ierr: error control, if different from 0 there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** **********************************************************************
  SUBROUTINE check_input(input_file_path,output_file_path,&
       echo_input,mpi_comm,ierr)
    CHARACTER(len=*), INTENT(in)             :: input_file_path, &
                                                output_file_path
    LOGICAL, OPTIONAL, intent(in)            :: echo_input
    INTEGER, intent(out)                     :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm

    character(len=*), parameter :: routineN = 'check_input', &
      routineP = moduleN//':'//routineN

    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    integer :: unit_nr
    TYPE(section_type), pointer :: input_structure
    TYPE(section_vals_type), pointer :: input_file
    TYPE(cp_parser_type), pointer :: cpparser
    LOGICAL :: my_echo_input
    TYPE(cp_unit_set_type), pointer :: default_units

    failure=.false.
    my_echo_input=.false.
    if (present(echo_input)) my_echo_input=echo_input

    IF (PRESENT(mpi_comm)) THEN
       NULLIFY(para_env)
       CALL cp_para_env_create(para_env, group=mpi_comm)
    ELSE
       para_env => default_para_env
       call cp_para_env_retain(para_env)
    END IF
    IF (para_env%mepos==para_env%source) THEN
       IF (output_file_path=="__STD_OUT__") THEN
          unit_nr=6
       ELSE
          call open_file(file_name=output_file_path,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=unit_nr)
       END IF
    ELSE
       unit_nr=-1
    END IF

    nullify(logger)
    CALL cp_logger_create(logger, para_env=para_env,&
         default_global_unit_nr=unit_nr, &
         close_global_unit_on_dealloc=.FALSE.)
    CALL cp_add_default_logger(logger)
    CALL cp_error_init(error,logger=logger,stop_level=cp_fatal_level)

    NULLIFY(input_structure, default_units)
    CALL create_cp2k_root_section(input_structure,error=error)
    NULLIFY(input_file)
    CALL section_vals_create(input_file,input_structure, error=error)
    CALL section_release(input_structure,error=error)
    nullify(cpparser)
    CALL parser_create(cpparser,file_name=input_file_path,para_env=para_env,&
         error=error)
    CALL cp_unit_set_create(default_units, "OUTPUT",error=error)
    CALL section_vals_parse(input_file,cpparser,root_section=.false.,&
         default_units=default_units,error=error)
    CALL cp_unit_set_release(default_units,error=error)
    CALL parser_release(cpparser,error=error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       IF (my_echo_input.AND.para_env%mepos==para_env%source) THEN
          CALL section_vals_write(input_file,&
               unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.),hide_root=.TRUE.,&
               error=error)
       END IF
       CALL check_cp2k_input(input_file,para_env=para_env,error=error)
    END IF
    CALL section_vals_release(input_file,error=error)

    CALL cp_logger_release(logger)
    CALL cp_para_env_release(para_env)
    ierr=cp_error_get_level(error)
    CALL cp_error_dealloc_ref(error)
    CALL cp_rm_default_logger()
  END SUBROUTINE check_input
  !***************************************************************************

END MODULE f77_interface
