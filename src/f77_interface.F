!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****h* cp2k/f77_interface [1.0] *
!!
!!   NAME
!!     f77_interface
!!
!!   FUNCTION
!!     interface to use cp2k as library
!!
!!   NOTES
!!     -
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_level,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                             cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_get_default_logger,&
                                             cp_logger_create,&
                                             cp_logger_release,&
                                             cp_logger_type,&
                                             cp_rm_default_logger,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             open_file
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
                                             cp2k_read
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create
  USE input_cp2k,                      ONLY: create_cp2k_root_section
  USE input_parsing,                   ONLY: section_vals_parse
  USE input_section_types,             ONLY: section_release,&
                                             section_type,&
                                             section_vals_create,&
                                             section_vals_type,&
                                             section_vals_write
  USE kg_main,                         ONLY: kg_create_force_env
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: MPI_COMM_WORLD,&
                                             mp_max,&
                                             mp_world_finalize,&
                                             mp_world_init
  USE particle_list_types,             ONLY: particle_list_type
  USE qmmm_init,                       ONLY: qmmm_create_force_env
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE timings,                         ONLY: timeset,&
                                             timestop

  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='f77_interface'
  INTEGER, parameter :: PATH_LENGTH=200

  TYPE(force_env_p_type), DIMENSION(:), POINTER, save :: force_envs
  INTEGER, DIMENSION(:), POINTER, SAVE :: f_env_ids
  TYPE(cp_error_type), DIMENSION(:), POINTER, save :: errors
  TYPE(cp_para_env_type), POINTER, save :: default_para_env
  LOGICAL, SAVE :: module_initialized=.FALSE.
  INTEGER, SAVE :: last_f_env_id=0

  PUBLIC :: init_cp2k, finalize_cp2k
  PUBLIC :: create_force_env, destroy_force_env, set_pos, get_pos,&
       get_force, calc_energy_force, run_input, get_energy, &
       calc_energy, calc_force
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
    DO isub=1,SIZE(f_env_ids)
       IF (f_env_ids(isub)==env_id) THEN
          env_pos=isub
       END IF
    END DO
    res=env_pos
  END FUNCTION get_pos_of_env

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
    TYPE(cp_para_env_type), pointer          :: para_env

    IF (.NOT. module_initialized) THEN
       IF (init_mpi) THEN
          ! get the default system wide communicator
          CALL mp_world_init(mpi_comm_default)
       ELSE
          mpi_comm_default=MPI_COMM_WORLD
       END IF

       nullify(default_para_env)
       CALL cp_para_env_create(default_para_env, group=mpi_comm_default)
       IF (para_env%source==para_env%mepos) THEN
          unit_nr=6
       ELSE
          unit_nr=-1
       END IF
       nullify(logger)
       CALL cp_logger_create(logger,para_env=para_env,&
            default_global_unit_nr=unit_nr, &
            close_global_unit_on_dealloc=.FALSE.)
       CALL cp_add_default_logger(logger)
       call cp_logger_release(logger)

       ALLOCATE(force_envs(0),f_env_ids(0),errors(0),stat=stat)
       CPPostconditionNoErr(stat==0,cp_fatal_level,routineP)
       module_initialized=.true.
       ierr=0
    ELSE
       ierr=cp_failure_level
    END IF

  END SUBROUTINE init_cp2k

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
       DO ienv=1,SIZE(force_envs)
          CALL force_env_release(force_envs(ienv)%force_env)
          CALL cp_error_dealloc_ref(errors(ienv))
       END DO
       DEALLOCATE(force_envs,f_env_ids,errors,stat=stat)
       CPPostconditionNoErr(stat==0,cp_warning_level,routineP)
       call cp_para_env_release(default_para_env)
       ierr=cp_error_get_level(error)
       call cp_rm_default_logger()
       IF (finalize_mpi) THEN
          CALL mp_world_finalize()
       END IF
    END IF
  END SUBROUTINE finalize_cp2k
  
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
  !!       it will be freed when you get rid of the force_env
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
    success=cp2k_run(input_file_path,unit_nr,para_env%group)
    IF (.NOT.success) THEN
       ierr=cp_failure_level
    ELSE
       ierr=0
    END IF
    call cp_para_env_release(para_env)
  END SUBROUTINE run_input


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
  !!       it will be freed when you get rid of the force_env
  !!     - ierr: will return a number different from 0 if there was an error
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
  SUBROUTINE create_force_env(new_env_id,input_file_path,&
       output_file_path,ierr,mpi_comm)
    INTEGER, INTENT(out)                     :: new_env_id
    CHARACTER(len=*), INTENT(in)             :: input_file_path, &
                                                output_file_path
    INTEGER, INTENT(out)                     :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm

    CHARACTER(len=*), PARAMETER :: routineN = 'create_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nfenv, stat, unit_nr
    INTEGER, DIMENSION(:), pointer           :: f_env_ids_old
    logical                                  :: failure
    type(cp_error_type)                      :: error
    TYPE(cp_error_type), DIMENSION(:), &
      pointer                                :: errors_old
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(cp_para_env_type), pointer          :: para_env
    TYPE(cp_parser_type), pointer            :: parser
    TYPE(force_env_p_type), DIMENSION(:), &
      pointer                                :: f_envs_old
    TYPE(global_environment_type), pointer   :: globenv
    TYPE(section_type), pointer              :: input_structure
    TYPE(section_vals_type), pointer         :: input_file

    failure=.false.

    nullify(para_env)
    IF (PRESENT(mpi_comm)) THEN
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
    CALL globenv_create(globenv, para_env=para_env,&
         output_unit=unit_nr, input_file_name=input_file_path,&
         close_output_unit=(unit_nr==6))
    call cp2k_init(globenv)
    logger => cp_get_default_logger()
    ! warning this is dangerous, I did not check that all the subfunctions
    ! support it, the program might crash upon error
    CALL cp_error_init(error,stop_level=cp_fatal_level,logger=logger)
    
    nfenv=SIZE(force_envs)
    f_envs_old => force_envs
    f_env_ids_old => f_env_ids
    errors_old => errors
    ALLOCATE(force_envs(nfenv+1),f_env_ids(nfenv+1),errors(nfenv+1),&
         stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    DO i=1,nfenv
       force_envs(i)%force_env => f_envs_old(i)%force_env
       errors(i)=errors_old(i)
       f_env_ids(i)=f_env_ids_old(i)
    END DO
    DEALLOCATE(f_envs_old,errors,f_envs_old,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

    nfenv=nfenv+1
    errors(nfenv)=error

    ! parse the input
    NULLIFY(input_structure)
    CALL create_cp2k_root_section(input_structure,error=error)
    NULLIFY(input_file)
    CALL section_vals_create(input_file,input_structure, error=error)
    CALL section_release(input_structure,error=error)
    nullify(parser)
    CALL parser_create(parser,file_name=input_file_path,error=error)
    CALL section_vals_parse(input_file,parser,root_section=.false.,&
         error=error)
    CALL parser_release(parser,error=error)
    IF (globenv%ionode) THEN
       CALL section_vals_write(input_file,unit_nr=globenv%scr,hide_root=.TRUE.,&
            error=error)
    END IF
    globenv%input_file => input_file
    nullify(input_file)
    
    nfenv=SIZE(force_envs)
    last_f_env_id=last_f_env_id+1
    new_env_id=last_f_env_id
    CALL mp_max(last_f_env_id,globenv%group)
    f_env_ids(nfenv)=new_env_id
    NULLIFY(force_envs(nfenv)%force_env)
    CALL cp2k_read(globenv)

    ! initialize force_env
    SELECT CASE (globenv%program_name)
    CASE ("FIST")
       CALL fist_create_force_env ( force_envs(nfenv)%force_env,&
            globenv,  error = error )
    CASE ("KG")
       CALL kg_create_force_env ( force_envs(nfenv)%force_env,&
            globenv, error=error)
    CASE ("QS")
       CALL quickstep_create_force_env(force_envs(nfenv)%force_env,&
            globenv, error=error)
    CASE ("QMMM")
       CALL qmmm_create_force_env(force_envs(nfenv)%force_env,&
            globenv, error=error)
    CASE default
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "Invalid run TYPE <"//TRIM(globenv%run_type)//"> was specified, "//&
            CPSourceFileRef,error,failure)
    END SELECT
    ierr=cp_error_get_level(error)
    call cp_rm_default_logger()

  END SUBROUTINE create_force_env
  
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

    integer                                  :: env_pos, i, nfenv, stat
    INTEGER, DIMENSION(:), pointer           :: f_env_ids_old
    logical                                  :: failure
    type(cp_error_type)                      :: error
    TYPE(cp_error_type), DIMENSION(:), &
      pointer                                :: errors_old
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(force_env_p_type), DIMENSION(:), &
      pointer                                :: f_envs_old
    TYPE(global_environment_type), pointer   :: globenv

    failure=.false.

    env_pos=get_pos_of_env(env_id)
    IF (env_pos>0) THEN
       call cp_error_init(error,template_error=errors(env_pos))
       logger => cp_error_get_logger(error)
       CALL cp_add_default_logger(logger)

       CALL force_env_get(force_envs(env_pos)%force_env,globenv=globenv,&
            error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL cp2k_finalize(globenv)
          CALL force_env_release(force_envs(env_pos)%force_env,error=error)
          CALL cp_error_check(error,failure)
          CALL cp_error_dealloc_ref(errors(env_pos))

          nfenv=SIZE(force_envs)
          f_envs_old => force_envs
          f_env_ids_old => f_env_ids
          errors_old => errors
          ALLOCATE(force_envs(nfenv-1),f_env_ids(nfenv-1),errors(nfenv-1),&
               stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          DO i=1,env_pos-1
             force_envs(i)%force_env => f_envs_old(i)%force_env
             errors(i)=errors_old(i)
             f_env_ids(i)=f_env_ids_old(i)
          END DO
          DO i=env_pos+1,nfenv
             force_envs(i-1)%force_env => f_envs_old(i)%force_env
             errors(i-1)=errors_old(i)
             f_env_ids(i-1)=f_env_ids_old(i)
          END DO
          DEALLOCATE(f_envs_old,errors,f_envs_old,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       ierr=cp_error_get_level(error)
       call cp_error_dealloc_ref(error)
    ELSE
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE destroy_force_env
!***************************************************************************

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
    INTEGER, INTENT(in)                      :: env_id
    REAL(kind=dp), DIMENSION(1:)             :: new_pos
    INTEGER, INTENT(in)                      :: n_el
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "set_pos", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub, &
                                                my_n_el
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles

    failure=.false.

    env_pos=get_pos_of_env(env_id)
    IF (env_pos>0) THEN
       logger => cp_error_get_logger(errors(env_pos))
       call cp_add_default_logger(logger)
       CALL cp_error_init(error,template_error=errors(env_pos))
       nullify(subsys)
       CALL force_env_get(force_envs(env_id)%force_env,&
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
       ierr=cp_error_get_level(error)
       CALL cp_error_dealloc_ref(error)
       CALL cp_rm_default_logger()
    ELSE
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE set_pos
 
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
    INTEGER, INTENT(in)                      :: env_id
    REAL(kind=dp), DIMENSION(1:)             :: pos
    INTEGER, INTENT(in)                      :: n_el
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "get_pos", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub, &
                                                my_n_el
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles

    failure=.false.

    env_pos=get_pos_of_env(env_id)
    IF (env_pos>0) THEN
       CALL cp_error_init(error,template_error=errors(env_pos))
       logger => cp_error_get_logger(error)
       CALL cp_add_default_logger(logger)
       IF (.NOT.failure) THEN
          CALL force_env_get(force_envs(env_pos)%force_env,&
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
          CALL cp_assert(my_n_el/=n_el,cp_failure_level,cp_assertion_failed,&
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
       ierr=cp_error_get_level(error)
       CALL cp_error_dealloc_ref(error)
       CALL cp_rm_default_logger()
    else
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE get_pos
   
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
    INTEGER, INTENT(in)                      :: env_id
    REAL(kind=dp), DIMENSION(1:)             :: force
    INTEGER, INTENT(in)                      :: n_el
    INTEGER, intent(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = "get_force", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: env_pos, idir, ii, ip, isub, &
                                                my_n_el
    logical                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), pointer            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), pointer                  :: subsys
    TYPE(particle_list_type), pointer        :: particles

    failure=.false.

    env_pos=get_pos_of_env(env_id)
    IF (env_pos>0) THEN
       CALL cp_error_init(error,template_error=errors(env_pos))
       logger => cp_error_get_logger(error)
       CALL cp_add_default_logger(logger)
       IF (.NOT.failure) THEN
          CALL force_env_get(force_envs(env_pos)%force_env,&
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
          CALL cp_assert(my_n_el/=n_el,cp_failure_level,cp_assertion_failed,&
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
       ierr=cp_error_get_level(error)
       CALL cp_error_dealloc_ref(error)
       CALL cp_rm_default_logger()
    else
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE get_force

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
    TYPE(cp_logger_type), pointer            :: logger

    failure=.false.

    env_pos=get_pos_of_env(env_id)
    IF (env_pos>0) THEN
       CALL cp_error_init(error,template_error=errors(env_pos))
       logger => cp_error_get_logger(error)
       CALL cp_add_default_logger(logger)

       CALL force_env_calc_energy_force(force_envs(env_pos)%force_env,&
            calc_force=calc_force,error=error)

       ierr=cp_error_get_level(error)
       CALL cp_error_dealloc_ref(error)
       CALL cp_rm_default_logger()
    else
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE calc_energy_force
  
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
    TYPE(cp_logger_type), pointer            :: logger

    failure=.false.

    env_pos=get_pos_of_env(env_id)
    IF (env_pos>0) THEN
       CALL cp_error_init(error,template_error=errors(env_pos))
       logger => cp_error_get_logger(error)
       CALL cp_add_default_logger(logger)

       CALL force_env_get(force_envs(env_pos)%force_env,&
            potential_energy=e_pot, error=error)

       ierr=cp_error_get_level(error)
       CALL cp_error_dealloc_ref(error)
       CALL cp_rm_default_logger()
    else
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE get_energy

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
    REAL(kind=dp), dimension(:), INTENT(in)  :: pos
    INTEGER, INTENT(in)                      :: n_el
    REAL(kind=dp), intent(out)               :: e_pot
    INTEGER, intent(out)                     :: ierr

    REAL(kind=dp), dimension(1)              :: dummy_f

    CALL calc_force(env_id,pos,n_el,e_pot,dummy_f,0,ierr)
  END SUBROUTINE calc_energy
  
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
    INTEGER, INTENT(in)                      :: env_id
    REAL(kind=dp), dimension(1:), INTENT(in) :: pos
    INTEGER, INTENT(in)                      :: n_el_pos
    REAL(kind=dp), intent(out)               :: e_pot
    REAL(kind=dp), dimension(1:), &
      INTENT(inout)                          :: force
    INTEGER, INTENT(in)                      :: n_el_force
    INTEGER, intent(out)                     :: ierr

    logical                                  :: calc_f

    calc_f=(n_el_force/=0)
    CALL set_pos(env_id,pos,n_el_pos,ierr)
    IF (ierr/=0) CALL calc_energy_force(env_id,calc_f,ierr)
    IF (ierr/=0) CALL get_energy(env_id,e_pot,ierr)
    IF (calc_f.AND.ierr/=0) CALL get_force(env_id,force,n_el_force,ierr)
  END SUBROUTINE calc_force
  
  
END MODULE f77_interface
