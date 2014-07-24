!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief various routines to log and control the output.
!>      The idea is that decisions about where to log should not be done in
!>      the code that generates the log, but should be globally changeable
!>      a central place.
!>      So some care has been taken to have enough information about the
!>      place from where the log comes so that in the future intelligent and
!>      flexible decisions can be taken by the logger, without having to change
!>      other code.
!> \note
!>      contains also routines to convert to a string.
!>      in my idea they should have been with variable length,
!>      (i.e. they should have returned a trim(adjustl(actual_result)))
!>      As a logger should be robust, at the moment I have given up.
!>
!>      At the moment logging and output refer to the same object
!>      (cp_logger_type)
!>      as these are actually different it might be better to separate them
!>      (they have already separate routines in a separate module
!>      @see cp_output_handling).
!>
!>      some practices (use of print *, no cp_error_type,
!>      manual retain release of some objects) are dictated by the need to
!>      have minimal dependency
!> \par History
!>      08.2002 major update: retain, release, printkeys, para_env,
!>              local logging [fawzi]
!>      02.2004 made a stack of default loggers [Joost VandeVondele]
!> \par
!>      @see cp_error_handling
!> \author Fawzi Mohamed
!>      @version 12.2001
! *****************************************************************************
MODULE cp_log_handling
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_iter_types,                   ONLY: cp_iteration_info_create,&
                                             cp_iteration_info_release,&
                                             cp_iteration_info_retain,&
                                             cp_iteration_info_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_datum,&
                                             m_flush,&
                                             m_getpid,&
                                             m_hostnm
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_comm_free
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: print_message
  USE timings,                         ONLY: print_stack

  IMPLICIT NONE
  PRIVATE

  !API types
  PUBLIC :: cp_logger_type
  !API parameter vars
  PUBLIC :: cp_note_level, cp_warning_level, cp_failure_level, cp_fatal_level
  !API default loggers
  PUBLIC :: cp_get_default_logger, cp_add_default_logger, cp_rm_default_logger,&
            cp_default_logger_stack_size
  !API logger routines
  PUBLIC :: cp_logger_create, cp_logger_retain, cp_logger_release,&
            cp_logger_would_log, cp_logger_set, cp_logger_get_default_unit_nr,&
            cp_logger_get_default_io_unit, cp_logger_get_unit_nr, &
            cp_logger_set_log_level, cp_logger_generate_filename, &
            cp_to_string, cp_log

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_log_handling'
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.

  !! level of an error
  INTEGER, PARAMETER          :: cp_fatal_level=3
  !! level of a failure
  INTEGER, PARAMETER          :: cp_failure_level=2
  !! level of a warning
  INTEGER, PARAMETER          :: cp_warning_level=1
  !! level of a note
  INTEGER, PARAMETER          :: cp_note_level=0

  !! a generic function to trasform different types to strings
  INTERFACE cp_to_string
     MODULE PROCEDURE cp_int_to_string, cp_real_dp_to_string, cp_logical_to_string
  END INTERFACE

  INTERFACE cp_log
     MODULE PROCEDURE cp_logger_log
  END INTERFACE

! *****************************************************************************
!> \brief type of a logger, at the moment it contains just a print level
!>      starting at which level it should be logged
!>      (0 note, 1 warning, 2 failure, 3 fatal)
!>      it could be expanded with the ability to focus on one or more
!>      module/object/thread/processor
!> \param id_nr unique number to identify the logger
!> \param ref_count reference count (see cp2k/doc/ReferenceCounting.html)
!> \param print_level the level starting at which something gets printed
!> \param default_local_unit_nr default unit for local logging (-1 if not
!>        yet initialized). Local logging guarantee to each task its own
!>        file.
!> \param default_global_unit_nr default unit for global logging
!>        (-1 if not yet initialized). This unit is valid only on the
!>        processor with %para_env%mepos==%para_env%source.
!> \param para_env the parallel environement for the output.
!>        this might be a super environment of your computation environment
!>        i.e. be very careful not to do global operations like broadcast
!>        with a subset of its processors (use your computation environment
!>        instead).
!> \param close_local_unit_on_dealloc if the local unit should be closed
!>        when this logger is deallocated
!> \param close_global_unit_on_dealloc whether the global unit should be
!>        closed when this logger is deallocated
!> \param suffix a short string that is used as suffix in all the filenames
!>        created by this logger. Can be used to guarantee the unicity of
!>        generated filename
!> \param local_filename the root of the name of the file used for local
!>        logging (can be different from the name of the file corresponding
!>        to default_local_unit_nr, only the one used if the unit needs to
!>        be opened)
!> \param global_filename the root of the name of the file used for
!>        global logging (can be different from the name of the file
!>        corresponding to default_global_unit_nr, only the one used if
!>        the unit needs to be opened)
!> \param print_keys print keys that tell what should be logged/outputted
!> \note
!>      This should be private, but as the output functions have been
!>      moved to another module and there is no "friend" keyword, they
!>      are public.
!>      DO NOT USE THE INTERNAL COMPONENTS DIRECTLY!!!
!> \par History
!>      04.2002 revised [fawzi]
!>      08.2002 major update: retain, release, printkeys, para_env,
!>              local logging [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_logger_type
     INTEGER :: id_nr, ref_count
     INTEGER :: print_level
     INTEGER :: default_local_unit_nr
     INTEGER :: default_global_unit_nr
     LOGICAL :: close_local_unit_on_dealloc, close_global_unit_on_dealloc
     CHARACTER(len=default_string_length)  :: suffix
     CHARACTER(len=default_path_length)    :: local_filename, global_filename
     TYPE(cp_para_env_type), POINTER       :: para_env
     TYPE(cp_iteration_info_type), POINTER :: iter_info
  END TYPE cp_logger_type

! *****************************************************************************
  TYPE default_logger_stack_type
     TYPE(cp_logger_type), POINTER :: cp_default_logger => Null()
  END TYPE default_logger_stack_type

  INTEGER, PRIVATE            :: stack_pointer=0
  INTEGER, PARAMETER, PRIVATE :: max_stack_pointer=10
  TYPE(default_logger_stack_type), SAVE, DIMENSION(max_stack_pointer) ::  default_logger_stack

  INTEGER, SAVE, PRIVATE :: last_logger_id_nr=0

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_default_logger_stack_size() RESULT(res)
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_default_logger_stack_size', &
      routineP = moduleN//':'//routineN

    res=stack_pointer
  END FUNCTION cp_default_logger_stack_size

! *****************************************************************************
!> \brief adds a default logger.
!>      MUST be called before logging occours
!> \param logger ...
!> \author Fawzi Mohamed
!> \note
!>      increments a stack of default loggers the latest one will be
!>      available within the program
! *****************************************************************************
  SUBROUTINE cp_add_default_logger(logger)
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_add_default_logger', &
      routineP = moduleN//':'//routineN

    IF (stack_pointer+1>max_stack_pointer) THEN
       CALL mp_abort(routineP//&
            "too many default loggers, increase max_stack_pointer in "//moduleN)
    ENDIF

    stack_pointer=stack_pointer+1
    NULLIFY(default_logger_stack(stack_pointer)%cp_default_logger)

    default_logger_stack(stack_pointer)%cp_default_logger => logger
    CALL cp_logger_retain(logger)

  END SUBROUTINE cp_add_default_logger

! *****************************************************************************
!> \brief the cousin of cp_add_default_logger, decrements the stack, so that
!>      the default logger is what it has
!>      been
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE cp_rm_default_logger()
    IF (stack_pointer-1<0) THEN
       CALL mp_abort(moduleN//":cp_rm_default_logger"//&
            "can not destroy default logger "//moduleN)
    ENDIF

    CALL cp_logger_release(default_logger_stack(stack_pointer)%cp_default_logger)
    NULLIFY(default_logger_stack(stack_pointer)%cp_default_logger)
    stack_pointer=stack_pointer-1

  END SUBROUTINE cp_rm_default_logger

! *****************************************************************************
!> \brief returns the default logger
!> \retval res ...
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      initializes the default loggers if necessary
! *****************************************************************************
  FUNCTION cp_get_default_logger() RESULT(res)
    TYPE(cp_logger_type), POINTER            :: res

    IF ( .NOT. stack_pointer>0 ) THEN
       CALL mp_abort("cp_log_handling:cp_get_default_logger"//&
            "default logger not yet initialized (CALL cp_init_default_logger)")
    END IF
    res => default_logger_stack(stack_pointer)%cp_default_logger
    IF (.NOT.ASSOCIATED(res)) THEN
       CALL mp_abort("cp_log_handling:cp_get_default_logger"//&
            "default logger is null (released too much ?)")
    END IF
  END FUNCTION cp_get_default_logger

! ================== log ==================

! *****************************************************************************
!> \brief initializes a logger
!> \param logger the logger to initialize
!> \param para_env the parallel environment (this is most likely the global
!>        parallel environment
!> \param print_level the level starting with which something is written
!>        (defaults to cp_note_level)
!> \param default_global_unit_nr the default unit_nr for output
!>        (if not given, and no file is given defaults to the standard output)
!> \param default_local_unit_nr the default unit number for local (i.e. task)
!>        output. If not given defaults to a out.taskid file created upon
!> \param global_filename a new file to open (can be given instread of the
!>        global_unit_nr)
!> \param local_filename a new file to open (with suffix and para_env%mepos
!>        appended). Can be given instread of the default_local_unit_nr).
!>        the file is created only upon the first local logging request
!> \param close_global_unit_on_dealloc if the unit should be closed when the
!>        logger is deallocated (defaults to true if a local_filename is given,
!>        to false otherwise)
!> \param iter_info ...
!> \param close_local_unit_on_dealloc if the unit should be closed when the
!>        logger is deallocated (defaults to true)
!> \param suffix the suffix that should be added to all the generated filenames
!> \param template_logger a logger from where to take the unspecified things
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      the handling of *_filename, default_*_unit_nr, close_*_unit_on_dealloc
!>      tries to take the right decision with different inputs, and thus is a
!>      little complex.
! *****************************************************************************
  SUBROUTINE cp_logger_create(logger, para_env,print_level,&
       default_global_unit_nr,default_local_unit_nr, global_filename,&
       local_filename, close_global_unit_on_dealloc,iter_info,&
       close_local_unit_on_dealloc, suffix, template_logger)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(in), OPTIONAL            :: print_level, &
                                                default_global_unit_nr, &
                                                default_local_unit_nr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: global_filename, &
                                                local_filename
    LOGICAL, INTENT(in), OPTIONAL            :: close_global_unit_on_dealloc
    TYPE(cp_iteration_info_type), OPTIONAL, &
      POINTER                                :: iter_info
    LOGICAL, INTENT(in), OPTIONAL            :: close_local_unit_on_dealloc
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: suffix
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: template_logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    ALLOCATE(logger, stat=stat)
    IF (stat/=0) CALL mp_abort(routineP//" could not ALLOCATE a logger")

    NULLIFY(logger%para_env)
    NULLIFY(logger%iter_info)
    logger%ref_count=1
    last_logger_id_nr=last_logger_id_nr+1
    logger%id_nr=last_logger_id_nr

    IF (PRESENT(template_logger)) THEN
       IF (template_logger%ref_count<1) &
            CALL mp_abort(routineP//" template_logger%ref_count<1")
       logger%print_level=template_logger%print_level
       logger%default_global_unit_nr=template_logger%default_global_unit_nr
       logger%close_local_unit_on_dealloc=template_logger%close_local_unit_on_dealloc
       IF (logger%close_local_unit_on_dealloc) THEN
          logger%default_local_unit_nr=-1
       ELSE
          logger%default_local_unit_nr=template_logger%default_local_unit_nr
       END IF
       logger%close_global_unit_on_dealloc=template_logger%close_global_unit_on_dealloc
       IF (logger%close_global_unit_on_dealloc) THEN
          logger%default_global_unit_nr=-1
       ELSE
          logger%default_global_unit_nr=template_logger%default_global_unit_nr
       END IF
       logger%local_filename=template_logger%local_filename
       logger%global_filename=template_logger%global_filename
       logger%para_env => template_logger%para_env
       logger%suffix=template_logger%suffix
       logger%iter_info => template_logger%iter_info
    ELSE
       ! create a file if nothing is specified, one can also get the unit from the default logger
       ! which should have something reasonable as the argument is required in that case
       logger%default_global_unit_nr=-1
       logger%close_global_unit_on_dealloc=.TRUE.
       logger%local_filename="localLog"
       logger%global_filename="mainLog"
       logger%print_level=cp_note_level
       ! generate a file for default local logger
       ! except the ionode that should write to the default global logger
       logger%default_local_unit_nr=-1
       logger%close_local_unit_on_dealloc=.TRUE.
       logger%suffix=""
    END IF
    IF (PRESENT(para_env)) logger%para_env => para_env
    IF (.NOT.ASSOCIATED(logger%para_env)) CALL mp_abort(routineP//&
         " para env not associated")
    IF (logger%para_env%ref_count<1) CALL mp_abort(routineP//&
         " para_env%ref_count<1")
    logger%para_env%ref_count = logger%para_env%ref_count+1

    IF (PRESENT(print_level)) logger%print_level=print_level

    IF (PRESENT(default_global_unit_nr)) &
         logger%default_global_unit_nr=default_global_unit_nr
    IF (PRESENT(global_filename)) THEN
       logger%global_filename=global_filename
       logger%close_global_unit_on_dealloc=.TRUE.
       logger%default_global_unit_nr=-1
    END IF
    IF (PRESENT(close_global_unit_on_dealloc)) THEN
       logger%close_global_unit_on_dealloc=close_global_unit_on_dealloc
       IF (PRESENT(default_global_unit_nr).AND.PRESENT(global_filename).AND.&
            (.NOT.close_global_unit_on_dealloc)) THEN
          logger%default_global_unit_nr=default_global_unit_nr
       END IF
    END IF

    IF (PRESENT(default_local_unit_nr)) &
         logger%default_local_unit_nr=default_local_unit_nr
    IF (PRESENT(local_filename)) THEN
       logger%local_filename=local_filename
       logger%close_local_unit_on_dealloc=.TRUE.
       logger%default_local_unit_nr=-1
    END IF
    IF (PRESENT(suffix)) logger%suffix=suffix

    IF (PRESENT(close_local_unit_on_dealloc)) THEN
       logger%close_local_unit_on_dealloc=close_local_unit_on_dealloc
       IF (PRESENT(default_local_unit_nr).AND.PRESENT(local_filename).AND.&
            (.NOT.close_local_unit_on_dealloc)) THEN
          logger%default_local_unit_nr=default_local_unit_nr
       END IF
    END IF

    IF (logger%default_local_unit_nr==-1) THEN
       IF (logger%para_env%mepos==logger%para_env%source) THEN
          logger%default_local_unit_nr=logger%default_global_unit_nr
          logger%close_local_unit_on_dealloc=.FALSE.
       ENDIF
    ENDIF
    IF (PRESENT(iter_info)) logger%iter_info => iter_info
    IF (ASSOCIATED(logger%iter_info)) THEN
       CALL cp_iteration_info_retain(logger%iter_info)
    ELSE
       CALL cp_iteration_info_create(logger%iter_info,"")
    END IF
  END SUBROUTINE cp_logger_create

! *****************************************************************************
!> \brief retains the given logger (to be called to keep a shared copy of
!>      the logger)
!> \param logger the logger to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_logger_retain(logger)
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.ASSOCIATED(logger)) CALL mp_abort(routineP//&
         " logger not associated")
    IF (logger%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")
    logger%ref_count=logger%ref_count+1
  END SUBROUTINE cp_logger_retain

! *****************************************************************************
!> \brief releases this logger
!> \param logger the logger to release
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_logger_release(logger)
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(logger)) THEN
       IF (logger%ref_count<1) CALL mp_abort(routineP//" logger%ref_count<1")
       logger%ref_count=logger%ref_count-1
       IF (logger%ref_count==0) THEN
          IF (logger%close_global_unit_on_dealloc .AND. &
               logger%default_global_unit_nr >= 0) THEN
             CALL close_file(logger%default_global_unit_nr)
             logger%close_global_unit_on_dealloc=.FALSE.
             logger%default_global_unit_nr=-1
          END IF
          IF (logger%close_local_unit_on_dealloc .AND. &
               logger%default_local_unit_nr >= 0) THEN
             CALL close_file(logger%default_local_unit_nr)
             logger%close_local_unit_on_dealloc=.FALSE.
             logger%default_local_unit_nr=-1
          END IF
          CALL my_cp_para_env_release(logger%para_env)
          CALL cp_iteration_info_release(logger%iter_info)
          DEALLOCATE(logger,stat=stat)
          IF (stat/=0) PRINT *, routineP," WARNING logger dealloc stat=",stat
       END IF
    END IF
    NULLIFY(logger)
  END SUBROUTINE cp_logger_release

! *****************************************************************************
!> \brief this function can be called to check if the logger would log
!>      a message with the given level from the given source
!>      you should use this function if you do direct logging
!>      (without using cp_logger_log), or if you want to know if the generation
!>      of some costly log info is necessary
!> \param logger the logger you want to log in
!> \param level describes the of the message: cp_fatal_level(3),
!>     cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!> \retval res ...
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_logger_would_log(logger, level) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(in)                      :: level
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_would_log', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: lggr

    lggr => logger
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")

    res = level>=lggr%print_level
  END FUNCTION cp_logger_would_log

! *****************************************************************************
!> \brief returns the unit nr for the requested kind of log.
!> \param logger the logger you want to log in
!> \param local if true returns a local logger (one per task), otherwise
!>     returns a global logger (only the process with para_env%mepos==
!>     para_env%source should write to the global logger). Defaults to
!>     false
!> \retval res ...
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_logger_get_unit_nr(logger, local) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL, INTENT(in), OPTIONAL            :: local
    INTEGER                                  :: res

    res=cp_logger_get_default_unit_nr(logger,local=local)
  END FUNCTION cp_logger_get_unit_nr


! *****************************************************************************
!> \brief returns the unit nr for the ionode (-1 on all other processors)
!>        skips as well checks if the procs calling this function is not the ionode
!> \param logger the logger you want to log in
!> \param local if true returns a local logger (one per task), otherwise
!>     returns a global logger (only the process with para_env%mepos==
!>     para_env%source should write to the global logger). Defaults to
!>     false UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \retval res ...
!> \par History
!>      12.2009 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  FUNCTION cp_logger_get_default_io_unit(logger) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: res

    res=cp_logger_get_default_unit_nr(logger,local=.FALSE.,skip_not_ionode=.TRUE.)
  END FUNCTION cp_logger_get_default_io_unit

! *****************************************************************************
!> \brief Logs in a standard format
!> \param logger is a type that describes where the logging must go.
!> \param level describes the of the message: cp_fatal_level(3),
!>     cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!> \param fromWhere is a string of the form module:function or file:lineNr
!>     that says where the error happend
!> \param message is the message to print
!> \param local if true this goes to the local logger, otherwise (the default)
!>     it goes to the global logger (and only the master process writes it)
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      In multithreaded situations it should be modified to print also
!>      the thread (and maybe the processor), and proper locking should be
!>      performed.
!>      At the moment it just goes to the default unit of the logger
! *****************************************************************************
  SUBROUTINE cp_logger_log(logger, level, fromWhere , message, local)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(in)                      :: level
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    LOGICAL, INTENT(in), OPTIONAL            :: local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_log', &
      routineP = moduleN//':'//routineN

    CHARACTER(7)                             :: label
    CHARACTER(LEN=26)                        :: datx
    INTEGER                                  :: unitNr
    LOGICAL                                  :: loc
    TYPE(cp_logger_type), POINTER            :: lggr

    lggr => logger
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")

    loc=.TRUE.
    IF (PRESENT(local)) loc=local
    IF (cp_logger_would_log(lggr,level)) THEN
       SELECT CASE (level)
       CASE (cp_fatal_level)
          label = 'ERRORL3'
       CASE (cp_failure_level)
          label = 'ERRORL2'
       CASE (cp_warning_level)
          label = 'WARNING'
       CASE (cp_note_level)
          label = 'NOTE'
       CASE default
          label = 'L'//ADJUSTL(cp_int_to_string(level))
       END SELECT
       IF (loc .OR. lggr%para_env%mepos==lggr%para_env%source) THEN
          unitNr = cp_logger_get_unit_nr(lggr,local=loc)
          CALL m_datum(datx)
          CALL print_message(datx(12:19)//' '//TRIM(ADJUSTL(label))//' in '//&
                             fromWhere//' :: '//TRIM(ADJUSTL(message)),&
                             unitNr,level,1,1)
          CALL m_flush(unitNr)
       END IF
    END IF
  END SUBROUTINE cp_logger_log

! *************************** cp_logger_type settings ***************************

! *****************************************************************************
!> \brief changes the logging level. Log messages with a level less than the one
!>      given wo not be printed.
!> \param logger the logger to change
!> \param level the new logging level for the logger
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_logger_set_log_level(logger, level)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(in)                      :: level

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_set_log_level', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(logger)) CALL mp_abort(routineP//&
         " logger not associated")
    IF (logger%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")
    logger%print_level=level
  END SUBROUTINE cp_logger_set_log_level

! *****************************************************************************
!> \brief returns the logging level of the given logger.
!> \param logger the logger you want info from
!> \retval res ...
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_logger_get_log_level(logger) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_get_log_level', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: lggr

    lggr => logger
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")

    res=lggr%print_level
  END FUNCTION cp_logger_get_log_level

! *****************************************************************************
!> \brief asks the default unit number of the given logger.
!>      try to use cp_logger_get_unit_nr
!> \param logger the logger you want info from
!> \param local if you want the local unit nr (defaults to false)
!> \param skip_not_ionode ...
!> \retval res ...
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  RECURSIVE FUNCTION cp_logger_get_default_unit_nr(logger, local, skip_not_ionode) RESULT(res)
    TYPE(cp_logger_type), OPTIONAL, POINTER  :: logger
    LOGICAL, INTENT(in), OPTIONAL            :: local, skip_not_ionode
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_logger_get_default_unit_nr', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename, host_name
    INTEGER                                  :: iostat, pid
    LOGICAL                                  :: loc, skip
    TYPE(cp_logger_type), POINTER            :: lggr

    loc=.TRUE.
    skip=.FALSE.
    IF (PRESENT(logger)) THEN
       lggr => logger
    ELSE
       NULLIFY(lggr)
    ENDIF
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")

    IF (PRESENT(local)) loc=local
    IF (PRESENT(skip_not_ionode)) skip=skip_not_ionode
    IF (.NOT.loc) THEN
       IF (lggr%default_global_unit_nr <= 0) THEN
          IF (lggr%para_env%mepos==lggr%para_env%source) THEN
             CALL cp_logger_generate_filename(lggr,filename,lggr%global_filename,&
                  ".out",local=.FALSE.)
             CALL open_file(TRIM(filename),file_status="unknown",&
                  file_action="WRITE",file_position="APPEND",&
                  unit_number=lggr%default_global_unit_nr)
          ELSE IF (.NOT.skip) THEN
             lggr%default_global_unit_nr=cp_logger_get_default_unit_nr(lggr,.TRUE.)
             lggr%close_global_unit_on_dealloc=.FALSE.
          ELSE
             lggr%default_global_unit_nr=-1
             lggr%close_global_unit_on_dealloc=.FALSE.
          END IF
       END IF
       IF ((lggr%para_env%mepos/=lggr%para_env%source).AND.(.NOT.skip)) THEN
          WRITE (UNIT=lggr%default_global_unit_nr,FMT='(/,T2,A)',IOSTAT=iostat)&
               ' *** WARNING non ionode asked for global logger ***'
          IF (iostat /= 0) THEN
             CALL m_getpid(pid)
             CALL m_hostnm(host_name)
             PRINT*," *** Error trying to WRITE to the local logger ***"
             PRINT*," *** MPI_id           = ",lggr%para_env%mepos
             PRINT*," *** MPI_Communicator = ",lggr%para_env%group
             PRINT*," *** PID              = ",pid
             PRINT*," *** Hostname         = "//TRIM(host_name)
             CALL print_stack(6)
          ELSE
             CALL print_stack(lggr%default_global_unit_nr)
          END IF
       END IF
       res=lggr%default_global_unit_nr
    ELSE
       IF (lggr%default_local_unit_nr<=0) THEN
          CALL cp_logger_generate_filename(lggr,filename,lggr%local_filename,&
               ".out",local=.TRUE.)
          CALL open_file(TRIM(filename),file_status="unknown",&
               file_action="WRITE",&
               file_position="APPEND",&
               unit_number=lggr%default_local_unit_nr)
          WRITE (UNIT=lggr%default_local_unit_nr,FMT='(/,T2,A,I0,A,I0,A)',IOSTAT=iostat)&
               '*** Local logger file of MPI task ',lggr%para_env%mepos,&
               ' in communicator ',lggr%para_env%group,' ***'
          IF (iostat == 0) THEN
             CALL m_getpid(pid)
             CALL m_hostnm(host_name)
             WRITE (UNIT=lggr%default_local_unit_nr,FMT='(T2,A,I0)',IOSTAT=iostat)&
                  '*** PID      = ',pid,&
                  '*** Hostname = '//host_name
             CALL print_stack(lggr%default_local_unit_nr)
          END IF
          IF (iostat /= 0) THEN
             CALL m_getpid(pid)
             CALL m_hostnm(host_name)
             PRINT*," *** Error trying to WRITE to the local logger ***"
             PRINT*," *** MPI_id           = ",lggr%para_env%mepos
             PRINT*," *** MPI_Communicator = ",lggr%para_env%group
             PRINT*," *** PID              = ",pid
             PRINT*," *** Hostname         = "//TRIM(host_name)
             CALL print_stack(6)
          END IF

       END IF
       res=lggr%default_local_unit_nr
    END IF
  END FUNCTION cp_logger_get_default_unit_nr

! *****************************************************************************
!> \brief releases the para object (to be called when you don't want anymore
!>      the shared copy of this object)
!> \param para_env the new group
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      local private version to avoid circular dependencies
!>      (see cp_para_env:cp_para_env_release)!
! *****************************************************************************
  SUBROUTINE my_cp_para_env_release(para_env)
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'my_cp_para_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(para_env)) THEN
       IF (para_env%ref_count<1) THEN
          CALL mp_abort(routineP//" para_env%ref_count<1")
       END IF
       para_env%ref_count=para_env%ref_count-1
       IF (para_env%ref_count<1) THEN
          IF (para_env%owns_group) THEN
             CALL mp_comm_free(para_env%group)
          END IF
          DEALLOCATE(para_env,stat=stat)
          IF(stat/=0)PRINT *,routineP," WARNING para_env DEALLOCATE stat=",stat
       END IF
    END IF
    NULLIFY(para_env)
  END SUBROUTINE my_cp_para_env_release

! *****************************************************************************
!> \brief generates a unique filename (ie adding eventual suffixes and
!>      process ids)
!> \param logger ...
!> \param res ...
!> \param root the start of filename
!> \param postfix the end of the name
!> \param local if the name should be local to this task (defaults to false)
!> \param result the resulting string UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      this should be a function returning a variable length string.
!>      All spaces are moved to the end of the string.
!>      Not fully optimized: result must be a littel longer than the
!>      resulting compressed filename
! *****************************************************************************
  SUBROUTINE cp_logger_generate_filename(logger,res, root, postfix,&
       local)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(inout)          :: res
    CHARACTER(len=*), INTENT(in)             :: root, postfix
    LOGICAL, INTENT(in), OPTIONAL            :: local

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_generate_filename', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, loc
    TYPE(cp_logger_type), POINTER            :: lggr

    failure=.FALSE.
    loc=.FALSE.
    res=' '
    lggr => logger

    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL mp_abort(routineP//&
         " logger%ref_count<1")
    IF (PRESENT(local)) loc=local
    IF (loc) THEN
       res=TRIM(root)//TRIM(lggr%suffix)//'_p'//&
            cp_to_string(lggr%para_env%mepos)//postfix
    ELSE
       res=TRIM(root)//TRIM(lggr%suffix)//postfix
    END IF
    CALL compress(res,full=.TRUE.)
  END SUBROUTINE cp_logger_generate_filename

! *****************************************************************************
!> \brief sets various attributes of the given logger
!> \param logger the logger you want to change
!> \param local_filename the root of the name of the file used for local
!>        logging
!> \param global_filename the root of the name of the file used for
!>        global logging
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_logger_set(logger,local_filename,global_filename)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: local_filename, &
                                                global_filename

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_set', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(logger)) CALL mp_abort(routineP//&
         " unassociated logger")
    IF (PRESENT(local_filename)) logger%local_filename=local_filename
    IF (PRESENT(global_filename)) logger%global_filename=global_filename
  END SUBROUTINE cp_logger_set

! *****************************************************************************
!> \brief converts an int to a string
!>      (should be a variable length string, but that does not work with
!>      all the compilers)
!> \param i the integer to convert
!> \retval res ...
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_int_to_string(i) RESULT(res)
    INTEGER, INTENT(in)                      :: i
    CHARACTER(len=6)                         :: res

    CHARACTER(len=6)                         :: t_res
    INTEGER                                  :: iostat
    REAL(KIND=dp)                            :: tmp_r

    IF (i>999999 .OR. i<-99999) THEN
       tmp_r=i
       WRITE (t_res,fmt='(es6.1)',iostat=iostat) tmp_r
    ELSE
       WRITE (t_res,fmt='(i6)',iostat=iostat) i
    END IF
    res = t_res
    IF (iostat/=0) THEN
       PRINT *,"cp_int_to_string ioerror",iostat
       CALL print_stack(cp_logger_get_default_unit_nr())
    END IF
  END FUNCTION cp_int_to_string

! *****************************************************************************
!> \brief convert a double precision real in a string
!>      (should be a variable length string, but that does not work with
!>      all the compilers)
!> \param val the number to convert
!> \retval res ...
!> \par History
!>      4.2002 revised [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_real_dp_to_string(val) RESULT(res)
    REAL(KIND=dp), INTENT(in)                :: val
    CHARACTER(len=11)                        :: res

    WRITE (res,'(es11.4)') val
  END FUNCTION cp_real_dp_to_string

! *****************************************************************************
!> \brief convert a logical in a string ('T' or 'F')
!> \param val the number to convert
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_logical_to_string(val) RESULT(res)
    LOGICAL, INTENT(in)                      :: val
    CHARACTER(len=1)                         :: res

    IF (val) THEN
       res='T'
    ELSE
       res='F'
    END IF
  END FUNCTION cp_logical_to_string

END MODULE cp_log_handling

