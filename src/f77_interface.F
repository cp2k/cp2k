!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief interface to use cp2k as library
!> \note
!>      useful additions for the future would be:
!>      - string(path) based set/get of simple values (to change the new
!>        input during the run and extract more data (energy types for example).
!>      - set/get of a subset of atoms
!> \par History
!>      07.2004 created [fawzi]
!>      11.2004 parallel version [fawzi]
!> \author fawzi & Johanna
! *****************************************************************************
MODULE f77_interface
  USE acc_device,                      ONLY: acc_get_ndevices,&
                                             acc_set_active_device
  USE bibliography,                    ONLY: add_all_references
  USE cell_types,                      ONLY: cell_type,&
                                             init_cell
  USE cp2k_info,                       ONLY: get_runtime_info
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_finalize_lib,&
                                             cp_dbcsr_init_lib
  USE cp_files,                        ONLY: init_preconnection_list,&
                                             open_file
  USE cp_ma_interface,                 ONLY: cp_ma_default_affinity,&
                                             cp_ma_finalize_lib,&
                                             cp_ma_init_lib,&
                                             cp_ma_mpi_reorder_strategy,&
                                             cp_ma_mpi_reordering,&
                                             cp_ma_set_mpi_reordering,&
                                             has_ma,&
                                             has_ma_topology
  USE cp_output_handling,              ONLY: cp_iterate
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: get_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_set,&
                                             cp_subsys_type,&
                                             unpack_subsys_particles
  USE eip_environment,                 ONLY: eip_init
  USE eip_environment_types,           ONLY: eip_env_create,&
                                             eip_env_release,&
                                             eip_environment_type
  USE environment,                     ONLY: cp2k_finalize,&
                                             cp2k_init,&
                                             cp2k_read,&
                                             cp2k_setup
  USE fist_main,                       ONLY: fist_create_force_env
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_create
  USE force_env_types,                 ONLY: &
       force_env_get, force_env_get_frc, force_env_get_natom, &
       force_env_get_nparticle, force_env_get_pos, force_env_get_vel, &
       force_env_release, force_env_retain, force_env_set, force_env_type, &
       multiple_fe_list
  USE force_env_utils,                 ONLY: force_env_shake
  USE fp_types,                        ONLY: fp_env_create,&
                                             fp_env_read,&
                                             fp_env_release,&
                                             fp_env_write,&
                                             fp_type
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_create,&
                                             globenv_release,&
                                             globenv_retain
  USE input_constants,                 ONLY: do_eip,&
                                             do_fist,&
                                             do_mixed,&
                                             do_qmmm,&
                                             do_qs
  USE input_cp2k_check,                ONLY: check_cp2k_input
  USE input_cp2k_force_eval,           ONLY: create_force_eval_section
  USE input_cp2k_read,                 ONLY: empty_initial_variables,&
                                             read_input
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: &
       section_get_keyword, section_release, section_type, &
       section_vals_duplicate, section_vals_get_subs_vals, &
       section_vals_release, section_vals_remove_values, section_vals_retain, &
       section_vals_type, section_vals_val_get, section_vals_write
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: default_output_unit,&
                                             m_chdir,&
                                             m_getcwd
  USE message_passing,                 ONLY: &
       MPI_COMM_WORLD, add_mp_perf_env, get_mp_perf_env, mp_max, &
       mp_perf_env_release, mp_perf_env_retain, mp_perf_env_type, &
       mp_world_finalize, mp_world_init, rm_mp_perf_env
  USE metadynamics_types,              ONLY: meta_env_release,&
                                             meta_env_type
  USE metadynamics_utils,              ONLY: metadyn_read
  USE mixed_environment_types,         ONLY: mixed_env_release,&
                                             mixed_environment_type
  USE mixed_main,                      ONLY: mixed_create_force_env
  USE periodic_table,                  ONLY: init_periodic_table
  USE pw_cuda,                         ONLY: pw_cuda_finalize,&
                                             pw_cuda_init
  USE qmmm_create,                     ONLY: qmmm_create_force_env
  USE qs_main,                         ONLY: quickstep_create_force_env
  USE reference_manager,               ONLY: remove_all_references
  USE string_table,                    ONLY: string_table_allocate,&
                                             string_table_deallocate
  USE timings,                         ONLY: add_timer_env,&
                                             get_timer_env,&
                                             rm_timer_env,&
                                             timer_env_release,&
                                             timer_env_retain,&
                                             timeset,&
                                             timestop
  USE timings_types,                   ONLY: timer_env_type
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'f77_interface'
#ifdef _USE_GA
#include "mafdecls.fh"
#endif

! *****************************************************************************
  TYPE f_env_p_type
     TYPE(f_env_type), POINTER :: f_env
  END TYPE f_env_p_type

! *****************************************************************************
  TYPE f_env_type
     INTEGER :: id_nr
     TYPE(force_env_type), POINTER      :: force_env
     TYPE(cp_error_type)                :: error
     TYPE(timer_env_type), POINTER      :: timer_env
     TYPE(mp_perf_env_type), POINTER    :: mp_perf_env
     CHARACTER(len=default_path_length) :: my_path,old_path
  END TYPE f_env_type

  TYPE(f_env_p_type), DIMENSION(:), POINTER, SAVE :: f_envs
  TYPE(cp_para_env_type), POINTER, SAVE :: default_para_env
  LOGICAL, SAVE :: module_initialized=.FALSE.
  INTEGER, SAVE :: last_f_env_id=0, n_f_envs=0
  INTEGER       :: mpi_mapping_method

  PUBLIC :: default_para_env
  PUBLIC :: init_f77_default
  PUBLIC :: init_cp2k, finalize_cp2k
  PUBLIC :: create_force_env, destroy_force_env, set_pos, get_pos,&
            get_force, calc_energy_force, get_energy, get_stress_tensor,&
            calc_energy, calc_force, check_input, get_natom, get_nparticle,&
            f_env_add_defaults, f_env_rm_defaults, f_env_type,&
            f_env_get_from_id,&
            set_vel, do_shake, set_cell, get_cell, get_result_r1
CONTAINS

! ****************************************************************************
!> \brief set the default mpi reordering method
!> \param mapping_method the selected method
!> \author Christiane Pousa Ribeiro
! *****************************************************************************
 SUBROUTINE init_f77_default(mapping_method)
    INTEGER, INTENT(IN)                      :: mapping_method

   mpi_mapping_method = mapping_method
 END SUBROUTINE init_f77_default

! *****************************************************************************
!> \brief returns the position of the force env corresponding to the given id
!> \param env_id the id of the requested environment
!> \retval res ...
!> \author fawzi
!> \note
!>      private utility function
! *****************************************************************************
  FUNCTION get_pos_of_env(env_id) RESULT(res)
    INTEGER, INTENT(in)                      :: env_id
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

! *****************************************************************************
!> \brief initializes cp2k, needs to be called once before using any of the
!>      other functions when using cp2k as library
!> \param init_mpi if the mpi environment should be initialized
!> \param ierr returns a number different from 0 if there was an error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE init_cp2k(init_mpi,ierr)
    LOGICAL, INTENT(in)                      :: init_mpi
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'init_cp2k', &
      routineP = moduleN//':'//routineN
 
#if defined __GEMINI || __SEASTAR || __BLUEGENE || __NET
    CHARACTER                                :: strategy
#endif

    INTEGER                                  :: mpi_comm_default, stat, &
                                                unit_nr
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger
 
#ifdef _USE_GA
    INTEGER                                  :: heap, stack
#endif

    IF (.NOT. module_initialized) THEN

       ! Initialise preconnection list
       CALL init_preconnection_list()

       ! get runtime information
       CALL get_runtime_info

       ! initialize with a default error (a proper error will be initialized later)
       CALL cp_error_init(error)

       IF (init_mpi) THEN
          ! get the default system wide communicator
          CALL mp_world_init(mpi_comm_default)
#ifdef _USE_GA
          CALL ga_initialize ( )
          heap = 1000000
          stack = 1000000
          IF (.not.ma_init(MT_DBL, stack, heap)) &
            CALL ga_error("ma_init failed in init_cp2k",-1)
#endif
       ELSE
          mpi_comm_default=MPI_COMM_WORLD
       END IF
 
       CALL string_table_allocate()
       CALL add_mp_perf_env()
       CALL add_timer_env()
       
#if defined __GEMINI || __SEASTAR || __BLUEGENE || __NET
       CALL cp_ma_mpi_reorder_strategy(mpi_mapping_method,strategy)
       CALL cp_ma_mpi_reordering(mpi_comm_default,mp_reorder=strategy, error=error)
#endif

       ! re-create the para_env and log with correct (reordered) ranks
       NULLIFY(default_para_env)
       CALL cp_para_env_create(default_para_env, group=mpi_comm_default, &
            owns_group=.FALSE.,error=error)
       IF (default_para_env%source==default_para_env%mepos) THEN
          unit_nr=default_output_unit
       ELSE
          unit_nr=-1
       END IF
       NULLIFY(logger)

       CALL cp_logger_create(logger,para_env=default_para_env,&
            default_global_unit_nr=unit_nr, &
            close_global_unit_on_dealloc=.FALSE.)
       CALL cp_add_default_logger(logger)
       CALL cp_logger_release(logger)

       ALLOCATE(f_envs(0),stat=stat)
       CPPostconditionNoErr(stat==0,cp_fatal_level,routineP)
       module_initialized=.TRUE.
       ierr=0

#if defined __GEMINI || __SEASTAR || __BLUEGENE || __NET
     CALL cp_ma_set_mpi_reordering(strategy)
#endif

       ! Initialize the MA_ARCH configuration
       IF (has_ma) THEN
         CALL cp_ma_init_lib(init_mpi,group=mpi_comm_default,error=error)
         !default binding only if mpi_init - all MPI tasks will be present
         CALL cp_ma_default_affinity(init_mpi)
       ENDIF

       !   *** Initialize mathematical constants ***
       CALL init_periodic_table()

       !   *** init the bibliography ***
       CALL add_all_references()


       IF(acc_get_ndevices() > 0) &
         CALL acc_set_active_device(MOD(default_para_env%mepos, acc_get_ndevices()))

       CALL pw_cuda_init()

       ! Initialize the DBCSR configuration
       CALL cp_dbcsr_init_lib (group=mpi_comm_default, error=error)

    ELSE
       ierr=cp_failure_level
    END IF

  END SUBROUTINE init_cp2k

! *****************************************************************************
!> \brief cleanup after you have finished using this interface
!> \param finalize_mpi if the mpi environment should be finalized
!> \param ierr returns a number different from 0 if there was an error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE finalize_cp2k(finalize_mpi,ierr)
    LOGICAL, INTENT(in)                      :: finalize_mpi
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'finalize_cp2k', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ienv, stat
    TYPE(cp_error_type)                      :: error

    IF (.NOT.module_initialized) THEN
       ierr=cp_failure_level
    ELSE
       ! Finalize the DBCSR configuration
       CALL cp_dbcsr_finalize_lib (group=default_para_env%group, output_unit=0, error=error)
       IF (has_ma .OR. has_ma_topology) THEN
         CALL cp_ma_finalize_lib(error=error)
       ENDIF
       CALL cp_error_init(error, stop_level=cp_fatal_level)
       CALL pw_cuda_finalize()
       DO ienv=n_f_envs,1,-1
          CALL destroy_force_env(f_envs(ienv)%f_env%id_nr,ierr=ierr)
          CPAssertNoFail(ierr==0,cp_warning_level,routineP,error)
       END DO
       DEALLOCATE(f_envs,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CALL cp_para_env_release(default_para_env,error)
       IF (cp_error_get_level(error)>CP_WARNING_LEVEL) THEN
          ierr=cp_error_get_level(error)
       ELSE
          ierr=0
       ENDIF
       CALL cp_rm_default_logger()

       !   *** deallocate the bibliography ***
       CALL remove_all_references()
       CALL rm_timer_env()
       CALL rm_mp_perf_env()
       CALL string_table_deallocate(0)
       IF (finalize_mpi) THEN
#ifdef _USE_GA
          CALL ga_terminate()
#endif
          CALL mp_world_finalize()
       END IF
    END IF
  END SUBROUTINE finalize_cp2k

! *****************************************************************************
!> \brief deallocates a f_env
!> \param f_env the f_env to deallocate
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE f_env_dealloc(f_env,error)
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'f_env_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(f_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL force_env_release(f_env%force_env,error=error)
       CALL cp_error_dealloc_ref(f_env%error,error=error)
       CALL timer_env_release(f_env%timer_env)
       CALL mp_perf_env_release(f_env%mp_perf_env)
       IF (f_env%old_path/=f_env%my_path) THEN
          CALL m_chdir(f_env%old_path,ierr)
          CPAssert(ierr==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE f_env_dealloc

! *****************************************************************************
!> \brief createates a f_env
!> \param f_env the f_env to createate
!> \param template_error the error to be used as template for the one stored
!> \param force_env the force_environment to be stored
!> \param timer_env the timer env to be stored
!> \param mp_perf_env the mp performance environement to be stored
!> \param id_nr ...
!> \param error ...
!> \param old_dir ...
!> \author fawzi
! *****************************************************************************
  SUBROUTINE f_env_create(f_env,template_error,force_env, timer_env,mp_perf_env,&
       id_nr,error,old_dir)
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: template_error
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(timer_env_type), POINTER            :: timer_env
    TYPE(mp_perf_env_type), POINTER          :: mp_perf_env
    INTEGER, INTENT(in)                      :: id_nr
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: old_dir

    CHARACTER(len=*), PARAMETER :: routineN = 'f_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(f_env,stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       f_env%force_env => force_env
       CALL force_env_retain(f_env%force_env,error=error)
       CALL cp_error_init(f_env%error,template_error=template_error)
       f_env%timer_env => timer_env
       CALL timer_env_retain(f_env%timer_env)
       f_env%mp_perf_env => mp_perf_env
       CALL mp_perf_env_retain(f_env%mp_perf_env)
       f_env%id_nr=id_nr
       CALL m_getcwd(f_env%my_path)
       f_env%old_path=f_env%my_path
       IF (PRESENT(old_dir)) f_env%old_path=old_dir
    END IF
  END SUBROUTINE f_env_create

! *****************************************************************************
!> \brief ...
!> \param f_env_id ...
!> \param f_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE f_env_get_from_id(f_env_id,f_env,error)
    INTEGER, INTENT(in)                      :: f_env_id
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'f_env_get_from_id', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: f_env_pos
    LOGICAL                                  :: failure

    failure=.FALSE.
    NULLIFY(f_env)
    f_env_pos = get_pos_of_env(f_env_id)
    IF (f_env_pos<1) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "invalid env_id "//cp_to_string(f_env_id),error=error,&
            failure=failure)
    ELSE
       f_env => f_envs(f_env_pos)%f_env
    END IF

  END SUBROUTINE f_env_get_from_id

! *****************************************************************************
!> \brief adds the default environments of the f_env to the stack of the
!>      defaults, and returns a new error and sets failure to true if
!>      something went wrong
!> \param f_env_id the f_env from where to take the defaults
!> \param f_env will contain the f_env corresponding to f_env_id
!> \param new_error an error that can be used for the given f_env
!> \param failure will be set to true if something went wrong
!> \param handle ...
!> \author fawzi
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
! *****************************************************************************
  SUBROUTINE f_env_add_defaults(f_env_id,f_env,new_error, failure, handle)
    INTEGER, INTENT(in)                      :: f_env_id
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(cp_error_type), INTENT(out)         :: new_error
    LOGICAL                                  :: failure
    INTEGER, INTENT(out), OPTIONAL           :: handle

    CHARACTER(len=*), PARAMETER :: routineN = 'f_env_add_defaults', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: f_env_pos, ierr
    LOGICAL                                  :: my_failure
    TYPE(cp_logger_type), POINTER            :: logger

    my_failure=.FALSE.

    NULLIFY(f_env)
    f_env_pos = get_pos_of_env(f_env_id)
    IF (f_env_pos<1) THEN
       CALL cp_error_init(new_error)
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "invalid env_id "//cp_to_string(f_env_id),error=new_error,&
            failure=failure)
    ELSE
       f_env => f_envs(f_env_pos)%f_env
       CALL cp_error_init(new_error,template_error=f_env%error)
       logger => cp_error_get_logger(f_env%error)
       CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,new_error,my_failure)
       IF (.NOT.my_failure) THEN
          CALL m_getcwd(f_env%old_path)
          IF (f_env%old_path/=f_env%my_path) THEN
             CALL m_chdir(TRIM(f_env%my_path),ierr)
             CPAssert(ierr==0,cp_failure_level,routineP,new_error,failure)
          END IF
          CALL add_mp_perf_env(f_env%mp_perf_env)
          CALL add_timer_env(f_env%timer_env)
          CALL cp_add_default_logger(logger)
          IF (PRESENT(handle)) handle=cp_default_logger_stack_size()
       END IF
       CALL cp_error_check(new_error,failure)
    END IF
  END SUBROUTINE f_env_add_defaults

! *****************************************************************************
!> \brief removes the default environments of the f_env to the stack of the
!>      defaults, and sets ierr accordingly to the failuers stored in error
!>      It also releases the error
!> \param f_env the f_env from where to take the defaults
!> \param error the error to be checked (will be released)
!> \param ierr variable that will be set to a number different from 0 if
!>        error contains an error (otherwise it will be set to 0)
!> \param handle ...
!> \author fawzi
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
! *****************************************************************************
  SUBROUTINE f_env_rm_defaults(f_env,error,ierr,handle)
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, INTENT(out), OPTIONAL           :: ierr
    INTEGER, INTENT(in), OPTIONAL            :: handle

    CHARACTER(len=*), PARAMETER :: routineN = 'f_env_rm_defaults', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr2
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: d_logger, logger
    TYPE(mp_perf_env_type), POINTER          :: d_mp_perf_env
    TYPE(timer_env_type), POINTER            :: d_timer_env

    failure=.FALSE.

    IF (ASSOCIATED(f_env)) THEN
       IF (PRESENT(handle)) THEN
          CPAssert(handle==cp_default_logger_stack_size(),cp_failure_level,routineP,error,failure)
       END IF

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
          ! CPPrecondition(d_timer_env%id_nr==f_env%timer_env%id_nr,cp_failure_level,routineP,error,failure)
          CPPrecondition(d_mp_perf_env%id_nr==f_env%mp_perf_env%id_nr,cp_failure_level,routineP,error,failure)
       END IF
       IF (f_env%old_path/=f_env%my_path) THEN
          CALL m_chdir(TRIM(f_env%old_path),ierr2)
          CPAssert(ierr2==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (PRESENT(ierr)) THEN
          IF (cp_error_get_level(error)>CP_WARNING_LEVEL) THEN
             ierr=cp_error_get_level(error)
          ELSE
             ierr=0
          ENDIF
       ENDIF
       CALL cp_error_dealloc_ref(error)
       IF (.NOT.failure) THEN
          CALL cp_rm_default_logger()
          CALL rm_timer_env()
          CALL rm_mp_perf_env()
       END IF
    ELSE
       IF (PRESENT(ierr)) THEN
          IF (cp_error_get_level(error)>CP_WARNING_LEVEL) THEN
             ierr=cp_error_get_level(error)
          ELSE
             ierr=0
          ENDIF
       ENDIF
       CALL cp_error_dealloc_ref(error)
    END IF
  END SUBROUTINE f_env_rm_defaults

! *****************************************************************************
!> \brief creates a new force environment using the given input, and writing
!>      the output to the given output unit
!> \param new_env_id will contain the id of the newly created environment
!> \param input_declaration ...
!> \param input_path where to read the input (if the input is given it can
!>        a virtual path)
!> \param output_path filename (or name of the unit) for the output
!> \param mpi_comm the mpi communicator to be used for this environment
!>        it will not be freed when you get rid of the force_env
!> \param output_unit if given it should be the unit for the output
!>        and no file is open(should be valid on the processor with rank 0)
!> \param owns_out_unit if the output unit should be closed upon destroing
!>        of the force_env (defaults to true if not default_output_unit)
!> \param input the parsed input, if given and valid it is used
!>        instead of parsing from file
!> \param ierr will return a number different from 0 if there was an error
!> \param work_dir ...
!> \param initial_variables ...
!> \author fawzi
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
! *****************************************************************************
  RECURSIVE SUBROUTINE create_force_env(new_env_id,input_declaration,input_path,&
       output_path,mpi_comm,output_unit,owns_out_unit,input,ierr,work_dir, initial_variables)
    INTEGER, INTENT(out)                     :: new_env_id
    TYPE(section_type), POINTER              :: input_declaration
    CHARACTER(len=*), INTENT(in)             :: input_path
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: output_path
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm, output_unit
    LOGICAL, INTENT(in), OPTIONAL            :: owns_out_unit
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    INTEGER, INTENT(out), OPTIONAL           :: ierr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: work_dir
    CHARACTER(len=*), DIMENSION(:, :), &
      OPTIONAL                               :: initial_variables

    CHARACTER(len=*), PARAMETER :: routineN = 'create_force_env', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: old_dir, wdir
    INTEGER :: handle, i, ierr2, iforce_eval, isubforce_eval, k, &
      method_name_id, my_group, nforce_eval, ngroups, nsubforce_size, stat, &
      unit_nr
    INTEGER, DIMENSION(:), POINTER           :: group_distribution, &
                                                i_force_eval, &
                                                lgroup_distribution
    LOGICAL :: check, failure, multiple_subsys, my_echo, my_owns_out_unit, &
      use_motion_section, use_multiple_para_env
    TYPE(cp_error_type)                      :: error, my_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: my_para_env, para_env
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(f_env_p_type), DIMENSION(:), &
      POINTER                                :: f_envs_old
    TYPE(force_env_type), POINTER            :: force_env, my_force_env
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(mp_perf_env_type), POINTER          :: mp_perf_env
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER :: fe_section, force_env_section, &
      force_env_sections, fp_section, input_file, root_section, &
      subsys_section, wrk_section
    TYPE(timer_env_type), POINTER            :: timer_env

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(input_declaration),cp_failure_level,routineP,error,failure)
    NULLIFY(para_env, force_env,timer_env,mp_perf_env, globenv,meta_env,&
            fp_env,eip_env,mixed_env)
    new_env_id=-1
    IF (PRESENT(mpi_comm)) THEN
       CALL cp_para_env_create(para_env, group=mpi_comm, owns_group=.FALSE.,error=error)
    ELSE
       para_env => default_para_env
       CALL cp_para_env_retain(para_env,error=error)
    END IF

    CALL timeset(routineN,handle)

    CALL m_getcwd(old_dir)
    wdir = old_dir
    IF (PRESENT(work_dir)) THEN
        IF (work_dir/=" ") THEN
            CALL m_chdir(work_dir,ierr2)
            IF (ierr2/=0) THEN
                IF (PRESENT(ierr)) ierr=ierr2
                RETURN
            END IF
            wdir = work_dir
        END IF
    END IF

    IF (PRESENT(output_unit)) THEN
       unit_nr=output_unit
    ELSE
       IF (para_env%mepos==para_env%source) THEN
          IF (output_path=="__STD_OUT__") THEN
             unit_nr=default_output_unit
          ELSE
             CALL open_file(file_name=output_path,file_status="UNKNOWN",&
                  file_action="WRITE", file_position="APPEND",&
                  unit_number=unit_nr)
          END IF
       ELSE
          unit_nr=-1
       END IF
    END IF
    my_owns_out_unit=unit_nr/=default_output_unit
    IF (PRESENT(owns_out_unit)) my_owns_out_unit=owns_out_unit
    CALL globenv_create(globenv,error=error)
    CALL cp2k_init(para_env, output_unit=unit_nr, globenv=globenv, input_file_name=input_path,&
                   wdir=wdir)
    logger => cp_get_default_logger()
    ! warning this is dangerous, I did not check that all the subfunctions
    ! support it, the program might crash upon error
    CALL cp_error_init(error,stop_level=cp_fatal_level,logger=logger)

    NULLIFY(input_file)
    IF (PRESENT(input)) input_file => input
    IF (.NOT.ASSOCIATED(input_file)) THEN
       IF (PRESENT(initial_variables)) THEN
         input_file => read_input(input_declaration, input_path,initial_variables,para_env=para_env,error=error)
       ELSE
         input_file => read_input(input_declaration, input_path,empty_initial_variables,para_env=para_env,error=error)
       ENDIF
    ELSE
       CALL section_vals_retain(input_file,error=error)
    END IF
    CALL section_vals_val_get(input_file,"GLOBAL%ECHO_INPUT",&
         l_val=my_echo,error=error)
    ! echo after check?
    IF (para_env%ionode.and.my_echo) THEN
       CALL section_vals_write(input_file,unit_nr=cp_logger_get_default_unit_nr(logger), &
               hide_root=.TRUE., hide_defaults=.FALSE., error=error)
    END IF
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! root_section => input_file
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXX
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL check_cp2k_input(input_declaration,input_file,para_env=para_env,output_unit=unit_nr,error=error)
       CALL cp_error_check(error,failure)
    END IF
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! NULLIFY(input_file)
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXX
    root_section => input_file
    CALL section_vals_retain(root_section,error)
    CALL cp_error_dealloc_ref(error)
    CALL cp_error_init(error,stop_level=cp_failure_level,logger=logger)

    IF (n_f_envs+1>SIZE(f_envs)) THEN
       f_envs_old => f_envs
       ALLOCATE(f_envs(n_f_envs+10), stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       DO i=1,n_f_envs
          f_envs(i)%f_env => f_envs_old(i)%f_env
       END DO
       DO i=n_f_envs+1,SIZE(f_envs)
          NULLIFY(f_envs(i)%f_env)
       END DO
       DEALLOCATE(f_envs_old,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF

    IF (.NOT.failure) THEN

       CALL cp2k_read(root_section,para_env,globenv,error=error)

       CALL cp2k_setup(root_section,para_env,globenv,error)
       ! Group Distribution
       ALLOCATE(group_distribution(0:para_env%num_pe-1),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       group_distribution  = 0
       lgroup_distribution => group_distribution
       ! Setup all possible force_env
       force_env_sections => section_vals_get_subs_vals(root_section,"FORCE_EVAL",error=error)
       CALL section_vals_val_get(root_section,"MULTIPLE_FORCE_EVALS%MULTIPLE_SUBSYS",&
            l_val=multiple_subsys,ignore_required=.TRUE.,error=error)
       CALL multiple_fe_list(force_env_sections, root_section, i_force_eval, nforce_eval,&
            error)
       ! Enforce the deletion of the subsys (unless not explicitly required)
       IF (.NOT.multiple_subsys) THEN
          DO iforce_eval = 2, nforce_eval
             wrk_section => section_vals_get_subs_vals(force_env_sections,"SUBSYS",&
                  i_rep_section=i_force_eval(iforce_eval),error=error)
             CALL section_vals_remove_values(wrk_section, error=error)
          END DO
       END IF
       nsubforce_size        = nforce_eval-1
       use_multiple_para_env = .FALSE.
       use_motion_section    = .TRUE.
       DO iforce_eval = 1, nforce_eval
          NULLIFY(force_env_section, my_force_env, subsys_section)
          ! Reference subsys from the first ordered force_eval
          IF (.NOT.multiple_subsys) THEN
             subsys_section => section_vals_get_subs_vals(force_env_sections,"SUBSYS",&
                  i_rep_section=i_force_eval(1),error=error)
          END IF
          ! Handling para_env in case of multiple force_eval
          IF (use_multiple_para_env) THEN
             ! Check that the order of the force_eval is the correct one
             CALL section_vals_val_get(force_env_sections, "METHOD", i_val=method_name_id, &
                  i_rep_section=i_force_eval(1), error=error)
             CALL cp_assert(method_name_id==do_mixed,cp_failure_level,cp_assertion_failed,routineP,&
                  "In case of multiple force_eval the MAIN force_eval (the first in the list of FORCE_EVAL_ORDER or "//&
                  "the one omitted from that order list) must be a MIXED_ENV type calculation. Please check your "//&
                  "input file and possibly correct the MULTIPLE_FORCE_EVAL%FORCE_EVAL_ORDER. "//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             check = ASSOCIATED(force_env%mixed_env%sub_para_env)
             CPPrecondition(check,cp_failure_level,routineP,error,failure)
             ngroups = force_env%mixed_env%ngroups
             my_group= lgroup_distribution(para_env%mepos)
             isubforce_eval = iforce_eval-1
             ! If task not allocated on this procs skip setup..
             IF (MODULO(isubforce_eval-1,ngroups)/=my_group) CYCLE
             my_para_env => force_env%mixed_env%sub_para_env(my_group+1)%para_env
             my_error = force_env%mixed_env%sub_error(my_group+1)
          ELSE
             my_para_env => para_env
             my_error = error
          END IF

          ! Initialize force_env_section
          ! No need to allocate one more force_env_section if only 1 force_eval
          ! is provided.. this is in order to save memory..
          IF (nforce_eval>1) THEN
             CALL section_vals_duplicate(force_env_sections,force_env_section,&
                  i_force_eval(iforce_eval),i_force_eval(iforce_eval),my_error)
             IF (iforce_eval/=1) use_motion_section = .FALSE.
          ELSE
             force_env_section => force_env_sections
             use_motion_section = .TRUE.
          END IF
          CALL section_vals_val_get(force_env_section,"METHOD",i_val=method_name_id,error=my_error)

          SELECT CASE (method_name_id)
          CASE (do_fist)
             CALL fist_create_force_env ( my_force_env, root_section, my_para_env, globenv,&
                  force_env_section=force_env_section, subsys_section=subsys_section,&
                  use_motion_section=use_motion_section, error=my_error )
          CASE (do_qs)
             CALL quickstep_create_force_env( my_force_env, root_section, my_para_env, globenv,&
                  force_env_section=force_env_section, subsys_section=subsys_section,&
                  use_motion_section=use_motion_section, error=my_error)
          CASE (do_qmmm)
             CALL qmmm_create_force_env( my_force_env, root_section, my_para_env, globenv,&
                  force_env_section=force_env_section, subsys_section=subsys_section,&
                  use_motion_section=use_motion_section, error=my_error)
          CASE (do_eip)
             CALL eip_env_create(eip_env, error=error)
             CALL eip_init(eip_env, root_section, para_env, force_env_section=force_env_section,&
                  subsys_section=subsys_section, use_motion_section=use_motion_section, error=error)
             CALL force_env_create(my_force_env, root_section, eip_env=eip_env, para_env=para_env,&
                  globenv=globenv, force_env_section=force_env_section, error=error)
             CALL eip_env_release(eip_env, error=error)
          CASE (do_mixed)
             CALL mixed_create_force_env(mixed_env, root_section, my_para_env, globenv,&
                  force_env_section=force_env_section, n_subforce_eval=nsubforce_size,&
                  use_motion_section=use_motion_section, error=my_error)
             CALL force_env_create(my_force_env, root_section, mixed_env=mixed_env, para_env=para_env,&
                  globenv=globenv, force_env_section=force_env_section, error=error)
             CALL mixed_env_release(mixed_env, error=error)
             use_multiple_para_env = .TRUE.
             lgroup_distribution => my_force_env%mixed_env%group_distribution
          CASE default
             CALL create_force_eval_section(section,error)
             keyword => section_get_keyword(section,"METHOD",error=error)
             CALL keyword_get(keyword,enum=enum,error=error)
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Invalid METHOD <"//TRIM(enum_i2c(enum,method_name_id,error=error))//&
               "> was specified, "//&
CPSourceFileRef,&
               my_error,failure)
             CALL section_release(section,error=error)
          END SELECT

          NULLIFY(meta_env, fp_env)
          IF (use_motion_section) THEN
             ! Metadynamics Setup
             fe_section => section_vals_get_subs_vals(root_section,"MOTION%FREE_ENERGY",error=my_error)
             CALL metadyn_read(meta_env, my_force_env, root_section, my_para_env, fe_section,&
                               error=my_error)
             CALL force_env_set(my_force_env,meta_env=meta_env,error=my_error)
             CALL meta_env_release(meta_env,error=my_error)
             ! Flexible Partition Setup
             fp_section => section_vals_get_subs_vals(root_section,"MOTION%FLEXIBLE_PARTITIONING",error=my_error)
             CALL fp_env_create(fp_env,error=my_error)
             CALL fp_env_read(fp_env,fp_section,error=my_error)
             CALL fp_env_write(fp_env,fp_section,error=my_error)
             CALL force_env_set(my_force_env,fp_env=fp_env,error=my_error)
             CALL fp_env_release(fp_env,error=my_error)
          END IF
          ! Handle multiple force_eval
          IF (nforce_eval>1.AND.iforce_eval==1) THEN
             ALLOCATE(my_force_env%sub_force_env(nsubforce_size),stat=stat)
             CPPrecondition(stat==0,cp_failure_level,routineP,my_error,failure)
             ! Nullify subforce_env
             DO k = 1,nsubforce_size
                NULLIFY(my_force_env%sub_force_env(k)%force_env)
             END DO
          END IF
          ! Reference the right force_env
          IF (iforce_eval==1) THEN
             force_env => my_force_env
          ELSE
             force_env%sub_force_env(iforce_eval-1)%force_env => my_force_env
          END IF
          ! Multiple para env for sub_force_eval
          IF (.NOT.use_multiple_para_env) THEN
             lgroup_distribution = iforce_eval
          END IF
          ! Release force_env_section
          IF (nforce_eval>1) CALL section_vals_release(force_env_section,my_error)
       END DO
       DEALLOCATE(group_distribution,stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(i_force_eval,stat=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    timer_env => get_timer_env()
    mp_perf_env => get_mp_perf_env()
    CALL mp_max(last_f_env_id,para_env%group)
    last_f_env_id=last_f_env_id+1
    new_env_id=last_f_env_id
    n_f_envs=n_f_envs+1
    CALL f_env_create(f_envs(n_f_envs)%f_env,template_error=error,&
         timer_env=timer_env,mp_perf_env=mp_perf_env,force_env=force_env,&
         id_nr=last_f_env_id,error=error,old_dir=old_dir)
    CALL force_env_release(force_env,error=error)
    CALL globenv_release(globenv,error=error)
    CALL section_vals_release(root_section,error)
    CALL cp_para_env_release(para_env,error=error)
    CALL f_env_rm_defaults(f_envs(n_f_envs)%f_env,error=error,ierr=ierr)
    CALL timestop(handle)

  END SUBROUTINE create_force_env

! *****************************************************************************
!> \brief deallocates the force_env with the given id
!> \param env_id the id of the force_env to remove
!> \param ierr will contain a number different from 0 if
!> \author fawzi
!> \note
!>      The following routines need to be synchronized wrt. adding/removing
!>      of the default environments (logging, perormance,error):
!>      environment:cp2k_init, environment:cp2k_finalize,
!>      f77_interface:f_env_add_defaults, f77_interface:f_env_rm_defaults,
!>      f77_interface:create_force_env, f77_interface:destroy_force_env
! *****************************************************************************
  RECURSIVE SUBROUTINE destroy_force_env(env_id,ierr)
    INTEGER, INTENT(in)                      :: env_id
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'destroy_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: env_pos, i, stat
    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: root_section

    failure=.FALSE.
    NULLIFY(f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.not.failure) THEN
       env_pos=get_pos_of_env(env_id)
       n_f_envs=n_f_envs-1
       DO i=env_pos,n_f_envs
          f_envs(i)%f_env => f_envs(i+1)%f_env
       END DO
       NULLIFY(f_envs(n_f_envs+1)%f_env)

       CALL force_env_get(f_env%force_env,globenv=globenv,&
            root_section=root_section,para_env=para_env,error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(globenv),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL globenv_retain(globenv,error=error)
          CALL f_env_dealloc(f_env,error=error)
          CALL cp2k_finalize(root_section,para_env,globenv,f_env%old_path,error)
          CALL section_vals_release(root_section,error=error)
          CALL globenv_release(globenv,error=error)
          CALL cp_error_check(error,failure)

          DEALLOCATE(f_env,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (cp_error_get_level(error)>CP_WARNING_LEVEL) THEN
          ierr=cp_error_get_level(error)
       ELSE
          ierr=0
       ENDIF
       CALL cp_error_dealloc_ref(error)
    ELSE
       ! write out something?
       ierr=cp_failure_level
    END IF
  END SUBROUTINE destroy_force_env

! *****************************************************************************
!> \brief returns the number of atoms in the given force env
!> \param env_id id of the force_env
!> \param n_atom ...
!> \param ierr will return a number different from 0 if there was an error
!> \date   22.11.2010 (MK)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE get_natom(env_id, n_atom, ierr)

    INTEGER, INTENT(IN)                      :: env_id
    INTEGER, INTENT(OUT)                     :: n_atom, ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_natom', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    n_atom = 0
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       n_atom = force_env_get_natom(f_env%force_env,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_natom

! *****************************************************************************
!> \brief returns the number of particles in the given force env
!> \param env_id id of the force_env
!> \param n_particle ...
!> \param ierr will return a number different from 0 if there was an error
!> \author Matthias Krack
!>
! *****************************************************************************
  SUBROUTINE get_nparticle(env_id, n_particle, ierr)

    INTEGER, INTENT(IN)                      :: env_id
    INTEGER, INTENT(OUT)                     :: n_particle, ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_nparticle', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    n_particle = 0
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       n_particle = force_env_get_nparticle(f_env%force_env,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_nparticle

! *****************************************************************************
!> \brief gets a cell
!> \param env_id id of the force_env
!> \param cell the array with the cell matrix
!> \param per periodicity
!> \param ierr will return a number different from 0 if there was an error
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE get_cell(env_id, cell, per, ierr)

    INTEGER, INTENT(IN)                      :: env_id
    REAL(KIND=DP), DIMENSION(3, 3)           :: cell
    INTEGER, DIMENSION(3), OPTIONAL          :: per
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell_full
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       NULLIFY (cell_full)
       CALL force_env_get(f_env%force_env,cell=cell_full,error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(cell_full),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          cell = cell_full%hmat
          IF(PRESENT(per)) per(:) = cell_full%perd(:)
       END IF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_cell


! *****************************************************************************
!> \brief gets a result from CP2K that is a real 1D array
!> \param env_id id of the force_env
!> \param description the tag of the result
!> \param N ...
!> \param RESULT ...
!> \param res_exist ...
!> \param ierr will return a number different from 0 if there was an error
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE get_result_r1(env_id, description, N, RESULT, res_exist, ierr)
    INTEGER                                  :: env_id
    CHARACTER(LEN=default_string_length)     :: description
    INTEGER                                  :: N
    REAL(KIND=dp), DIMENSION(1:N)            :: RESULT
    LOGICAL, OPTIONAL                        :: res_exist
    INTEGER                                  :: ierr

    INTEGER                                  :: nres
    LOGICAL                                  :: exist_res, failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_result_type), POINTER            :: results
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env, subsys, results)
    CALL f_env_add_defaults(env_id,f_env,error,failure)

    IF (.NOT. failure) THEN
      CALL force_env_get(f_env%force_env, subsys=subsys, error=error)
      CALL cp_subsys_get(subsys, results=results, error=error)
      ! first test for the result
      IF(PRESENT(res_exist)) THEN
        res_exist = test_for_result(results, description=description, error=error)
        exist_res = res_exist
      ELSE
        exist_res = .TRUE.
      END IF
      ! if existing (or assuming the existance) read the results
      IF(res_exist) THEN
        CALL get_results(results, description=description, n_rep=nres,error=error)
        CALL get_results(results, description=description, values=RESULT, nval=nres,error=error)
      END IF
    ENDIF

    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_result_r1

! *****************************************************************************
!> \brief gets the forces of the particles
!> \param env_id id of the force_env
!> \param frc the array where to write the forces
!> \param n_el number of positions (3*nparticle) just to check
!> \param ierr will return a number different from 0 if there was an error
!> \date   22.11.2010 (MK)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE get_force(env_id, frc, n_el, ierr)

    INTEGER, INTENT(IN)                      :: env_id, n_el
    REAL(KIND=dp), DIMENSION(1:n_el)         :: frc
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_force', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_get_frc(f_env%force_env,frc,n_el,error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_force

! *****************************************************************************
!> \brief gets the stress tensor
!> \param env_id id of the force_env
!> \param stress_tensor the array where to write the stress tensor
!> \param ierr will return a number different from 0 if there was an error
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE get_stress_tensor(env_id, stress_tensor, ierr)

    INTEGER, INTENT(IN)                      :: env_id
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: stress_tensor
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_stress_tensor', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY (f_env, subsys, virial, cell)
    stress_tensor(:,:) = 0.0_dp

    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_get(f_env%force_env, subsys=subsys, cell=cell, error=error)
       CALL cp_subsys_get(subsys, virial=virial, error=error)
       IF(virial%pv_availability) THEN
          stress_tensor(:,:) = virial%pv_virial(:,:)/cell%deth
       ENDIF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_stress_tensor

! *****************************************************************************
!> \brief gets the positions of the particles
!> \param env_id id of the force_env
!> \param pos the array where to write the positions
!> \param n_el number of positions (3*nparticle) just to check
!> \param ierr will return a number different from 0 if there was an error
!> \date   22.11.2010 (MK)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE get_pos(env_id, pos, n_el, ierr)

    INTEGER, INTENT(IN)                      :: env_id, n_el
    REAL(KIND=DP), DIMENSION(1:n_el)         :: pos
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_pos', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_get_pos(f_env%force_env,pos,n_el,error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_pos

! *****************************************************************************
!> \brief gets the velocities of the particles
!> \param env_id id of the force_env
!> \param vel the array where to write the velocities
!> \param n_el number of velocities (3*nparticle) just to check
!> \param ierr will return a number different from 0 if there was an error
!> \author fawzi
!> date    22.11.2010 (MK)
! *****************************************************************************
  SUBROUTINE get_vel(env_id, vel, n_el, ierr)

    INTEGER, INTENT(IN)                      :: env_id, n_el
    REAL(KIND=DP), DIMENSION(1:n_el)         :: vel
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_vel', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_get_vel(f_env%force_env,vel,n_el,error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_vel

! *****************************************************************************
!> \brief sets a new cell
!> \param env_id id of the force_env
!> \param new_cell the array with the cell matrix
!> \param ierr will return a number different from 0 if there was an error
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE set_cell(env_id, new_cell, ierr)

    INTEGER, INTENT(IN)                      :: env_id
    REAL(KIND=DP), DIMENSION(3, 3)           :: new_cell
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_cell', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env, cell, subsys)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       NULLIFY (cell)
       CALL force_env_get(f_env%force_env,cell=cell,error=error)
       CALL cp_error_check(error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(ASSOCIATED(cell),cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          cell%hmat = new_cell
          CALL init_cell(cell)
          CALL force_env_get(f_env%force_env, subsys=subsys, error=error)
          CALL cp_subsys_set(subsys, cell=cell, error=error)
       END IF
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE set_cell

! *****************************************************************************
!> \brief sets the positions of the particles
!> \param env_id id of the force_env
!> \param new_pos the array with the new positions
!> \param n_el number of positions (3*nparticle) just to check
!> \param ierr will return a number different from 0 if there was an error
!> \date   22.11.2010 updated (MK)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE set_pos(env_id, new_pos, n_el, ierr)

    INTEGER, INTENT(IN)                      :: env_id, n_el
    REAL(KIND=dp), DIMENSION(1:n_el)         :: new_pos
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_pos', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       NULLIFY (subsys)
       CALL force_env_get(f_env%force_env,subsys=subsys,error=error)
       CALL cp_error_check(error,failure)
       CALL unpack_subsys_particles(subsys=subsys,r=new_pos,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE set_pos

! *****************************************************************************
!> \brief sets the velocities of the particles
!> \param env_id id of the force_env
!> \param new_vel the array with the new velocities
!> \param n_el number of velocities (3*nparticle) just to check
!> \param ierr will return a number different from 0 if there was an error
!> \date   22.11.2010 updated (MK)
!> \author fawzi
! *****************************************************************************
  SUBROUTINE set_vel(env_id, new_vel, n_el, ierr)

    INTEGER, INTENT(IN)                      :: env_id, n_el
    REAL(kind=dp), DIMENSION(1:n_el)         :: new_vel
    INTEGER, INTENT(OUT)                     :: ierr

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_vel', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       NULLIFY (subsys)
       CALL force_env_get(f_env%force_env,subsys=subsys,error=error)
       CALL cp_error_check(error,failure)
       CALL unpack_subsys_particles(subsys=subsys,v=new_vel,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE set_vel

! *****************************************************************************
!> \brief updates the energy and the forces of given force_env
!> \param env_id id of the force_env that you want to update
!> \param calc_force if the forces should be updated, if false the forces
!>        might be wrong.
!> \param ierr will return a number different from 0 if there was an error
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE calc_energy_force(env_id,calc_force,ierr)

    INTEGER, INTENT(in)                      :: env_id
    LOGICAL, INTENT(in)                      :: calc_force
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_energy_force', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       logger => cp_error_get_logger(error)
       CALL cp_iterate(logger%iter_info,error=error) ! add one to the iteration count
       CALL force_env_calc_energy_force(f_env%force_env,calc_force=calc_force,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE calc_energy_force

! *****************************************************************************
!> \brief returns the energy of the last configuration calculated
!> \param env_id id of the force_env that you want to update
!> \param e_pot the potential energy of the system
!> \param ierr will return a number different from 0 if there was an error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE get_energy(env_id,e_pot,ierr)

    INTEGER, INTENT(in)                      :: env_id
    REAL(kind=dp), INTENT(out)               :: e_pot
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'get_energy', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure = .FALSE.
    NULLIFY (f_env)
    CALL f_env_add_defaults(env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_get(f_env%force_env,potential_energy=e_pot,error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)

  END SUBROUTINE get_energy

! *****************************************************************************
!> \brief returns the energy of the configuration given by the positions
!>      passed as argument
!> \param env_id id of the force_env that you want to update
!> \param pos array with the positions
!> \param n_el number of elements in pos (3*natom)
!> \param e_pot the potential energy of the system
!> \param ierr will return a number different from 0 if there was an error
!> \author fawzi
!> \note
!>      utility call
! *****************************************************************************
  RECURSIVE SUBROUTINE calc_energy(env_id,pos,n_el,e_pot,ierr)

    INTEGER, INTENT(IN)                      :: env_id, n_el
    REAL(KIND=dp), DIMENSION(1:n_el), &
      INTENT(IN)                             :: pos
    REAL(KIND=dp), INTENT(OUT)               :: e_pot
    INTEGER, INTENT(OUT)                     :: ierr

    REAL(KIND=dp), DIMENSION(1)              :: dummy_f

    CALL calc_force(env_id,pos,n_el,e_pot,dummy_f,0,ierr)

  END SUBROUTINE calc_energy

! *****************************************************************************
!> \brief returns the energy of the configuration given by the positions
!>      passed as argument
!> \param env_id id of the force_env that you want to update
!> \param pos array with the positions
!> \param n_el_pos number of elements in pos (3*natom)
!> \param e_pot the potential energy of the system
!> \param force array that will contain the forces
!> \param n_el_force number of elements in force (3*natom). If 0 the
!>        forces are not calculated
!> \param ierr will return a number different from 0 if there was an error
!> \author fawzi
!> \note
!>      utility call, but actually it could be a better and more efficient
!>      interface to connect to other codes if cp2k would be deeply
!>      refactored
! *****************************************************************************
  RECURSIVE SUBROUTINE calc_force(env_id,pos,n_el_pos,e_pot,force,n_el_force,ierr)

    INTEGER, INTENT(in)                      :: env_id, n_el_pos
    REAL(kind=dp), DIMENSION(1:n_el_pos), &
      INTENT(in)                             :: pos
    REAL(kind=dp), INTENT(out)               :: e_pot
    INTEGER, INTENT(in)                      :: n_el_force
    REAL(kind=dp), DIMENSION(1:n_el_force), &
      INTENT(inout)                          :: force
    INTEGER, INTENT(out)                     :: ierr

    LOGICAL                                  :: calc_f

    calc_f = (n_el_force /= 0)
    CALL set_pos(env_id,pos,n_el_pos,ierr)
    IF (ierr == 0) CALL calc_energy_force(env_id,calc_f,ierr)
    IF (ierr == 0) CALL get_energy(env_id,e_pot,ierr)
    IF (calc_f.AND.(ierr == 0)) CALL get_force(env_id,force,n_el_force,ierr)

  END SUBROUTINE calc_force

! *****************************************************************************
!> \brief performs a check of the input
!> \param input_declaration ...
!> \param input_file_path the path of the input file to check
!> \param output_file_path path of the output file (to which it is appended)
!>        if it is "__STD_OUT__" the default_output_unit is used
!> \param echo_input if the parsed input should be written out with all the
!>        defaults made explicit
!> \param mpi_comm the mpi communicator (if not given it uses the default
!>        one)
!> \param permissive ...
!> \param ierr error control, if different from 0 there was an error
!> \author fawzi
! *****************************************************************************
  SUBROUTINE check_input(input_declaration, input_file_path,output_file_path,&
       echo_input,mpi_comm,permissive,ierr)
    TYPE(section_type), POINTER              :: input_declaration
    CHARACTER(len=*), INTENT(in)             :: input_file_path, &
                                                output_file_path
    LOGICAL, INTENT(in), OPTIONAL            :: echo_input
    INTEGER, INTENT(in), OPTIONAL            :: mpi_comm
    LOGICAL, INTENT(in), OPTIONAL            :: permissive
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'check_input', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr
    LOGICAL                                  :: failure, my_echo_input, &
                                                my_permissive
    TYPE(cp_error_type)                      :: error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input_file

    failure=.FALSE.
    my_echo_input=.FALSE.
    my_permissive=.FALSE.
    IF (PRESENT(permissive)) my_permissive=permissive
    IF (PRESENT(echo_input)) my_echo_input=echo_input

    IF (PRESENT(mpi_comm)) THEN
       NULLIFY(para_env)
       CALL cp_para_env_create(para_env, group=mpi_comm,error=error)
    ELSE
       para_env => default_para_env
       CALL cp_para_env_retain(para_env,error=error)
    END IF
    IF (para_env%mepos==para_env%source) THEN
       IF (output_file_path=="__STD_OUT__") THEN
          unit_nr=default_output_unit
       ELSE
          CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
               file_action="WRITE", file_position="APPEND",&
               unit_number=unit_nr)
       END IF
    ELSE
       unit_nr=-1
    END IF

    NULLIFY(logger)
    CALL cp_logger_create(logger, para_env=para_env,&
         default_global_unit_nr=unit_nr, &
         close_global_unit_on_dealloc=.FALSE.)
    CALL cp_add_default_logger(logger)
    CALL cp_logger_release(logger)
    CALL cp_error_init(error,logger=logger,stop_level=cp_fatal_level)

    input_file => read_input(input_declaration, input_file_path,initial_variables=empty_initial_variables, &
                                            para_env=para_env,error=error)
    IF (my_permissive) CALL cp_error_reset(error)
    CALL cp_error_check(error,failure)
    IF (.NOT.failure) THEN
       CALL check_cp2k_input(input_declaration,input_file,para_env=para_env,output_unit=unit_nr,error=error)
       IF (my_permissive) CALL cp_error_reset(error)
       IF (my_echo_input.AND.para_env%mepos==para_env%source) THEN
          CALL section_vals_write(input_file,&
               unit_nr=cp_logger_get_default_unit_nr(logger,local=.FALSE.),hide_root=.TRUE.,&
               hide_defaults=.FALSE.,error=error)
       END IF
    END IF
    CALL section_vals_release(input_file,error=error)

    CALL cp_logger_release(logger)
    CALL cp_para_env_release(para_env,error=error)
    ierr=cp_error_get_level(error)
    CALL cp_error_dealloc_ref(error)
    CALL cp_rm_default_logger()
  END SUBROUTINE check_input

! *****************************************************************************
!> \brief perform the shake procedure (enforce constraints)
!> \param f_env_id the force env on which to apply the shake
!> \param dt timestep
!> \param shake_tol tolerance for shake
!> \param ierr ...
!> \author fawzi
! *****************************************************************************
  SUBROUTINE do_shake(f_env_id,dt,shake_tol,ierr)
    INTEGER, INTENT(in)                      :: f_env_id
    REAL(kind=dp), INTENT(in)                :: dt, shake_tol
    INTEGER, INTENT(out)                     :: ierr

    CHARACTER(len=*), PARAMETER :: routineN = 'do_shake', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_error_type)                      :: error
    TYPE(f_env_type), POINTER                :: f_env

    failure=.FALSE.

    NULLIFY(f_env)
    CALL f_env_add_defaults(f_env_id,f_env,error,failure)
    IF (.NOT.failure) THEN
       CALL force_env_shake(f_env%force_env,&
            dt=dt,shake_tol=shake_tol, error=error)
    END IF
    CALL f_env_rm_defaults(f_env,error,ierr)
  END SUBROUTINE do_shake

END MODULE f77_interface
