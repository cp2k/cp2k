!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****h* cp2k/cp_log_handling [1.0] *
!!
!!   NAME
!!     cp_log_handling
!!
!!   FUNCTION
!!     various routines to log and control the output.
!!     The idea is that decisions about where to log should not be done in
!!     the code that generates the log, but should be globally changeable
!!     a central place.
!!     So some care has been taken to have enough information about the 
!!     place from where the log comes so that in the future intelligent and
!!     flexible decisions can be taken by the logger, without having to change
!!     other code.
!!
!!   NOTES
!!     contains also routines to convert to a string.
!!     in my idea they should have been with variable length,
!!     (i.e. they should have returned a trim(adjustl(actual_result)))
!!     As a logger should be robust, at the moment I have given up.
!!
!!     At the moment logging and output refer to the same object 
!!     (cp_logger_type)
!!     as these are actually different it might be better to separate them
!!     (they have already separate routines in a separate module 
!!     @see cp_output_handling).
!!
!!     some practices (use of print *, stop_program, no cp_error_type,
!!     manual retain release of some objects) are dictated by the need to
!!     have minimal dependency
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     08.2002 major update: retain, release, printkeys, para_env, 
!!             local logging [fawzi]
!!
!!     @see cp_error_handling
!!
!!   SOURCE
!******************************************************************************

MODULE cp_log_handling
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: init_print_keys,&
                                             print_key_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_abort,&
                                             m_datum,&
                                             m_flush,&
                                             m_getpid,&
                                             m_hostnm
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_free
  USE output_utilities,                ONLY: print_message
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: print_stack

  IMPLICIT NONE
  PRIVATE

! types
  PUBLIC :: cp_logger_type
! parameter vars
  PUBLIC :: cp_note_level, cp_warning_level, cp_failure_level,&
       cp_fatal_level, cp_std_in_unit_nr, cp_std_out_unit_nr
! default loggers
  PUBLIC :: cp_get_default_logger, cp_init_default_logger
! logger routines
  PUBLIC :: cp_logger_create, cp_logger_retain, cp_logger_release,&
       cp_logger_would_log,&
       cp_logger_get_default_unit_nr, cp_logger_get_unit_nr, &
       cp_logger_set_log_level, cp_logger_generate_filename
  PUBLIC :: cp_to_string, cp_log
! underlying routines
!!FM  public :: cp_int_to_string, cp_real_dp_to_string

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_log_handling'
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  INTEGER, PARAMETER :: cp_std_in_unit_nr=5
  INTEGER, PARAMETER :: cp_std_out_unit_nr=6

!! level of an error
  INTEGER, PARAMETER  :: cp_fatal_level=3
!! level of a failure
  INTEGER, PARAMETER  :: cp_failure_level=2
!! level of a warning
  INTEGER, PARAMETER  :: cp_warning_level=1
!! level of a note
  INTEGER, PARAMETER  :: cp_note_level=0

!! a generic function to trasform different types to strings
  INTERFACE cp_to_string
     MODULE PROCEDURE cp_int_to_string, cp_real_dp_to_string
  END INTERFACE
  INTERFACE cp_log
     MODULE PROCEDURE cp_logger_log
  END INTERFACE

  LOGICAL, SAVE, PRIVATE :: default_logger_initialized=.FALSE.
  INTEGER, SAVE, PRIVATE :: last_logger_id_nr=0

!!***
!****************************************************************************

!!****s* cp_log_handling/cp_logger_type [1.0] *
!!
!!   NAME
!!     cp_logger_type
!!
!!   FUNCTION
!!     type of a logger, at the moment it contains just a print level
!!     starting at which level it should be logged 
!!     (0 note, 1 warning, 2 failure, 3 fatal)
!!     it could be expanded with the ability to focus on one or more
!!     module/object/thread/processor
!!
!!   NOTES
!!     This should be private, but as the output functions have been 
!!     moved to another module and there is no "friend" keyword, they
!!     are public.
!!     DO NOT USE THE INTERNAL COMPONENTS DIRECTLY!!!
!!
!!   ATTRIBUTES
!!     - id_nr: unique number to identify the logger
!!     - ref_count: reference count (see cp2k/doc/ReferenceCounting.html)
!!     - print_level: the level starting at which something gets printed
!!     - default_local_unit_nr: default unit for local logging (-1 if not
!!       yet initialized). Local logging guarantee to each task its own
!!       file.
!!     - default_global_unit_nr: default unit for global logging
!!       (-1 if not yet initialized). This unit is valid only on the 
!!       processor with %para_env%mepos==%para_env%source.
!!     - para_env: the parallel environement for the output.
!!       this might be a super environment of your computation environment
!!       i.e. be very careful not to do global operations like broadcast
!!       with a subset of its processors (use your computation environment
!!       instead).
!!     - close_local_unit_on_dealloc: if the local unit should be closed
!!       when this logger is deallocated
!!     - close_global_unit_on_dealloc: whether the global unit should be
!!       closed when this logger is deallocated
!!     - suffix: a short string that is used as suffix in all the filenames
!!       created by this logger. Can be used to guarantee the unicity of
!!       generated filename
!!     - local_filename: the root of the name of the file used for local
!!       logging
!!     - global_filename: the root of the name of the file used for 
!!       global logging
!!     - print_keys: print keys that tell what should be logged/outputted
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2002 revised [fawzi]
!!     08.2002 major update: retain, release, printkeys, para_env, 
!!             local logging [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_logger_type
     INTEGER :: id_nr, ref_count
     INTEGER :: print_level
     INTEGER :: default_local_unit_nr
     INTEGER :: default_global_unit_nr
     LOGICAL :: close_local_unit_on_dealloc, close_global_unit_on_dealloc
     CHARACTER(len=10) :: suffix
     CHARACTER(len=50) :: local_filename, global_filename
     TYPE(print_key_type), POINTER :: print_keys
     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE cp_logger_type
!!***
!****************************************************************************

! the default logger
  TYPE(cp_logger_type), POINTER, PRIVATE ::cp_default_logger ! => null()

CONTAINS

!!****f* cp_log_handling/cp_init_default_logger [1.0] *
!!
!!   NAME
!!     cp_init_default_loggers
!!
!!   FUNCTION
!!     performs the initializations of the logger.
!!     MUST be called before logging occours
!!
!!   NOTES
!!     initializes only if necessary
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_init_default_logger(para_env,print_level,&
       default_global_unit_nr,default_local_unit_nr, global_filename,&
       local_filename, close_global_unit_on_dealloc,&
       close_local_unit_on_dealloc, print_keys, suffix,&
       template_logger)
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(in), OPTIONAL            :: print_level, &
                                                default_global_unit_nr, &
                                                default_local_unit_nr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: global_filename, &
                                                local_filename
    LOGICAL, INTENT(in), OPTIONAL            :: close_global_unit_on_dealloc, &
                                                close_local_unit_on_dealloc
    TYPE(print_key_type), OPTIONAL, TARGET   :: print_keys
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: suffix
    TYPE(cp_logger_type), OPTIONAL, TARGET   :: template_logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_init_default_logger', &
      routineP = moduleN//':'//routineN

    IF (default_logger_initialized) THEN
       PRINT *,routineP," WARNING called more than once"
       RETURN
       CALL cp_logger_release(cp_default_logger)
    ELSE
       NULLIFY(cp_default_logger)
    END IF
    ! cases because of bug with optional pointers arg in pgf90 3.3
    IF (PRESENT(para_env)) THEN
       CALL cp_logger_create(cp_default_logger,&
            para_env=para_env,print_level=print_level,&
            default_global_unit_nr=default_global_unit_nr,&
            default_local_unit_nr=default_local_unit_nr,&
            global_filename=global_filename,&
            local_filename=local_filename,&
            close_global_unit_on_dealloc=close_global_unit_on_dealloc,&
            close_local_unit_on_dealloc=close_local_unit_on_dealloc,&
            suffix=suffix, print_keys=print_keys,&
            template_logger=template_logger)
    ELSE
       CALL cp_logger_create(cp_default_logger,&
            print_level=print_level,&
            default_global_unit_nr=default_global_unit_nr,&
            default_local_unit_nr=default_local_unit_nr,&
            global_filename=global_filename,&
            local_filename=local_filename,&
            close_global_unit_on_dealloc=close_global_unit_on_dealloc,&
            close_local_unit_on_dealloc=close_local_unit_on_dealloc,&
            suffix=suffix, print_keys=print_keys,&
            template_logger=template_logger)
    END IF
    default_logger_initialized=.TRUE.
  END SUBROUTINE cp_init_default_logger
!**************************************************************************

!!****f* cp_log_handling/cp_get_default_logger [1.0] *
!!
!!   NAME
!!     cp_get_default_logger
!!
!!   SYNOPSIS
!!     Function cp_get_default_logger() Result(res)
!!       Type(cp_logger_type), Pointer:: res
!!     End Function cp_get_default_logger
!!
!!   FUNCTION
!!     returns the default logger
!!
!!   NOTES
!!     initializes the default loggers if necessary
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION cp_get_default_logger() RESULT(res)
    TYPE(cp_logger_type), POINTER            :: res

    IF (.NOT. default_logger_initialized) THEN
       CALL stop_program("cp_log_handling:cp_get_default_logger",&
            "default logger not yet initialized (CALL cp_init_default_logger")
    END IF
    IF (.NOT.ASSOCIATED(cp_default_logger)) THEN
       CALL stop_program("cp_log_handling:cp_get_default_logger",&
            "default logger was null (released one time too much??)")
    END IF
    res => cp_default_logger
  END FUNCTION cp_get_default_logger
!***************************************************************************

! ================== log ==================

!!****f* cp_log_handling/cp_logger_create [1.0] *
!!
!!   NAME
!!     cp_logger_create
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_create(logger, para_env, print_level,&
!!         default_global_unit_nr, default_local_unit_nr, global_filename,&
!!         local_filename, close_global_unit_on_dealloc,&
!!         close_local_unit_on_dealloc, print_keys, suffix,&
!!         template_logger)
!!       Type(cp_logger_type), Pointer:: logger
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Integer, Intent (IN), Optional:: print_level,&
!!         default_global_unit_nr, default_local_unit_nr
!!       Character(Len=*), Intent (IN), Optional:: global_filename,&
!!         local_filename
!!       Logical, Intent (IN), Optional:: close_global_unit_on_dealloc,&
!!         close_local_unit_on_dealloc
!!       Type(print_key_type), Target, Optional:: print_keys
!!       Character(Len=*), Intent (IN), Optional:: suffix
!!       Type(cp_logger_type), Target, Optional:: template_logger
!!     End Subroutine cp_logger_create
!!
!!   FUNCTION
!!     initializes a logger
!!
!!   NOTES
!!     the handling of *_filename, default_*_unit_nr, close_*_unit_on_dealloc
!!     tries to take the right decision with different inputs, and thus is a
!!     little complex.
!!
!!   INPUTS
!!     - logger: the logger to initialize
!!     - para_env: the parallel environment (this is most likely the global
!!       parallel environment
!!     - print_level: the level starting with which something is written
!!       (defaults to cp_note_level)
!!     - default_global_unit_nr: the default unit_nr for output
!!       (if not given, and no file is given defaults to the standard output)
!!     - default_local_unit_nr: the default unit number for local (i.e. task)
!!       output. If not given defaults to a out.taskid file created upon 
!!     - global_filename: a new file to open (can be given instread of the 
!!       global_unit_nr)
!!     - local_filename: a new file to open (with suffix and para_env%mepos 
!!       appended). Can be given instread of the default_local_unit_nr).
!!       the file is created only upon the first local logging request
!!     - close_global_unit_on_dealloc: if the unit should be closed when the
!!       logger is deallocated (defaults to true if a local_filename is given,
!!       to false otherwise)
!!     - close_local_unit_on_dealloc: if the unit should be closed when the
!!       logger is deallocated (defaults to true)
!!     - suffix: the suffix that should be added to all the generated filenames
!!     - template_logger: a logger from where to take the unspecified things
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_logger_create(logger, para_env,print_level,&
       default_global_unit_nr,default_local_unit_nr, global_filename,&
       local_filename, close_global_unit_on_dealloc,&
       close_local_unit_on_dealloc, print_keys, suffix, template_logger)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, INTENT(in), OPTIONAL            :: print_level, &
                                                default_global_unit_nr, &
                                                default_local_unit_nr
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: global_filename, &
                                                local_filename
    LOGICAL, INTENT(in), OPTIONAL            :: close_global_unit_on_dealloc, &
                                                close_local_unit_on_dealloc
    TYPE(print_key_type), OPTIONAL, TARGET   :: print_keys
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: suffix
    TYPE(cp_logger_type), OPTIONAL, TARGET   :: template_logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    ALLOCATE(logger, stat=stat)
    IF (stat/=0) CALL stop_program(routineP,"could not ALLOCATE a logger")

    NULLIFY(logger%para_env,logger%print_keys)
    logger%ref_count=1
    last_logger_id_nr=last_logger_id_nr+1
    logger%id_nr=last_logger_id_nr

    IF (PRESENT(template_logger)) THEN
      IF (template_logger%ref_count<1) &
            CALL stop_program(routineP,"template_logger%ref_count<1")
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
      logger%print_keys => template_logger%print_keys
      logger%para_env => template_logger%para_env
      logger%suffix=template_logger%suffix
   ELSE
    logger%default_global_unit_nr=cp_std_out_unit_nr
    logger%local_filename="localLog"
    logger%global_filename="mainLog"
    logger%print_level=cp_note_level
    logger%close_global_unit_on_dealloc=.FALSE.
#ifdef __parallel
    logger%default_local_unit_nr=-1
    logger%close_local_unit_on_dealloc=.TRUE.
#else
    logger%default_local_unit_nr=cp_std_out_unit_nr
    logger%close_local_unit_on_dealloc=.FALSE.
#endif
    logger%suffix=""
   END IF
    IF (PRESENT(para_env)) logger%para_env => para_env
    IF (.NOT.ASSOCIATED(logger%para_env)) CALL stop_program(routineP,&
         "para env not associated")
    IF (logger%para_env%ref_count<1) CALL stop_program(routineP,&
         "para_env%ref_count<1")
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
            .NOT.close_global_unit_on_dealloc) THEN
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
    IF (PRESENT(close_local_unit_on_dealloc)) THEN
       logger%close_local_unit_on_dealloc=close_local_unit_on_dealloc
       IF (PRESENT(default_local_unit_nr).AND.PRESENT(local_filename).AND.&
            .NOT.close_local_unit_on_dealloc) THEN
         logger%default_local_unit_nr=default_local_unit_nr
       END IF
    END IF

    IF (PRESENT(print_keys)) logger%print_keys => print_keys
    IF (ASSOCIATED(logger%print_keys)) THEN
       IF (logger%print_keys%ref_count<1) THEN
          CALL stop_program(routineP,"print_keys%ref_count<1")
       END IF
       logger%print_keys%ref_count=logger%print_keys%ref_count+1
    ELSE
       ALLOCATE(logger%print_keys,stat=stat)
       IF (stat/=0) THEN
          CALL stop_program(routineP,"print_keys allocation stat="//&
               cp_to_string(stat))
          CALL init_print_keys(logger%print_keys,level=3-logger%print_level)
       END IF
    END IF
    ! first global log to safely b-cast global unitnr
    stat=cp_logger_get_default_unit_nr(logger, local=.FALSE.)
    ! it would also be possible (if you don't care about printing a note)
    ! CALL cp_logger_log(logger,level=cp_note_level,&
    !      fromWhere=routineP, message="global logging activated",&
    !      local=.FALSE.)
  END SUBROUTINE cp_logger_create
!***************************************************************************

!!****f* cp_log_handling/cp_logger_retain [1.0] *
!!
!!   NAME
!!     cp_logger_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_retain(logger)
!!       Type(cp_logger_type), Pointer:: logger
!!     End Subroutine cp_logger_retain
!!
!!   FUNCTION
!!     retains the given logger (to be called to keep a shared copy of 
!!     the logger)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_logger_retain(logger)
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.ASSOCIATED(logger)) CALL stop_program(routineP,&
         "logger not associated")
    IF (logger%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")
    logger%ref_count=logger%ref_count+1
  END SUBROUTINE cp_logger_retain
!***************************************************************************

!!****f* cp_log_handling/cp_logger_release *
!!
!!   NAME
!!     cp_logger_release
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_release(logger)
!!       Type(cp_logger_type), Pointer:: logger
!!     End Subroutine cp_logger_release
!!
!!   FUNCTION
!!     releases this logger
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger to release
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_logger_release(logger)
    TYPE(cp_logger_type), POINTER            :: logger

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(logger)) THEN
       IF (logger%ref_count<1) CALL stop_program(routineP,"logger%ref_count<1")
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
          IF (ASSOCIATED(logger%print_keys)) THEN
             IF (logger%print_keys%ref_count<1) THEN
                CALL stop_program(routineP," logger%print_keys%ref_count<1")
             END IF
             logger%print_keys%ref_count=logger%print_keys%ref_count-1
             IF (logger%print_keys%ref_count==0) THEN
                CALL m_abort()
                DEALLOCATE(logger%print_keys,stat=stat)
                IF (stat/=0) PRINT *, routineP,&
                     " WARNING print_keys dealloc stat=",stat
             END IF
          END IF
          CALL my_cp_para_env_release(logger%para_env)
          DEALLOCATE(logger,stat=stat)
          IF (stat/=0) PRINT *, routineP," WARNING logger dealloc stat=",stat
       END IF
    END IF
    NULLIFY(logger)
  END SUBROUTINE cp_logger_release
!***************************************************************************

!!****f* cp_log_handling/cp_logger_would_log [1.0] *
!!
!!   NAME
!!     cp_logger_would_log
!!
!!   SYNOPSIS
!!     Function cp_logger_would_log(logger, level, fromwhere, local)&
!!         Result(res)
!!       Logical:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Logical, Intent (IN), Optional:: local
!!     End Function cp_logger_would_log
!!
!!   FUNCTION
!!     this function can be called to check if the logger would log
!!     a message with the given level from the given source
!!     you should use this function if you do direct logging
!!     (without using cp_logger_log), or if you want to know if the generation 
!!     of some costly log info is necessary
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!  - logger: the logger you want to log in
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - local: if the output is local (defaults to false)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure
  FUNCTION cp_logger_would_log(logger, level, fromWhere, local) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(in)                      :: level
    CHARACTER(len=*), INTENT(in)             :: fromWhere
    LOGICAL, INTENT(in), OPTIONAL            :: local
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_would_log', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: lggr

    lggr => logger
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")

    res = level>=lggr%print_level
  END FUNCTION cp_logger_would_log

!!****f* cp_log_handling/cp_logger_get_unit_nr [1.0] *
!!
!!   NAME
!!     cp_logger_get_unit_nr
!!
!!   SYNOPSIS
!!     Function cp_logger_get_unit_nr(logger, level, fromwhere, local)&
!!         Result(res)
!!       Integer:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Logical, Intent (IN), Optional:: local
!!     End Function cp_logger_get_unit_nr
!!
!!   FUNCTION
!!     returns the unit nr for the requested kind of log.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!  - logger: the logger you want to log in
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - local: if true returns a local logger (one per task), otherwise
!!    returns a global logger (only the process with para_env%mepos==
!!    para_env%source should write to the global logger). Defaults to
!!    false
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure 
  FUNCTION cp_logger_get_unit_nr(logger, level, fromWhere, local) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(in)                      :: level
    CHARACTER(len=*), INTENT(in)             :: fromWhere
    LOGICAL, INTENT(in), OPTIONAL            :: local
    INTEGER                                  :: res

    res=cp_logger_get_default_unit_nr(logger,local=local)
  END FUNCTION cp_logger_get_unit_nr

!!****f* cp_log_handling/cp_logger_log [1.0] *
!!
!!   NAME
!!     cp_logger_log
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_log(logger, level, fromwhere, message)
!!       Type(cp_logger_type):: logger
!!       Integer, Intent (IN):: level
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!     End Subroutine cp_logger_log
!!
!!   FUNCTION
!!     Logs in a standard format
!!
!!   NOTES
!!     In multithreaded situations it should be modified to print also 
!!     the thread (and maybe the processor), and proper locking should be 
!!     performed.
!!     At the moment it just goes to the default unit of the logger
!!
!!   INPUTS
!!  - logger is a type that describes where the logging must go.
!!  - level describes the of the message: cp_fatal_level(3), 
!!    cp_failure_level(2), cp_warning_level(1), cp_note_level(0).
!!  - fromWhere is a string of the form module:function or file:lineNr
!!    that says where the error happend
!!  - message is the message to print
!!  - local: if true this goes to the local logger, otherwise (the default)
!!    it goes to the global logger (and only the master process writes it)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
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
    IF (lggr%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")

    loc=.FALSE.
    IF (PRESENT(local)) loc=local
    IF (cp_logger_would_log(lggr,level,fromWhere,local=loc)) THEN
       SELECT CASE (level)
       CASE (cp_fatal_level)
          label='ERRORL3'
       CASE (cp_failure_level)
          label='ERRORL2'
       CASE (cp_warning_level)
          label='WARNING'
       CASE (cp_note_level)
          label='NOTE'
       CASE default
          label='L'//ADJUSTL(cp_int_to_string(level))
       END SELECT
       unitNr=cp_logger_get_unit_nr(lggr,level,fromWhere,local=loc)
       IF (loc .OR. lggr%para_env%mepos==lggr%para_env%source) THEN
          CALL m_datum(datx)
          CALL print_message(datx(12:19)//' '//label//' in '//&
               fromWhere//' '//message,unitNr,&
               level, level,level/2)
! old (more failsafe, no trim, compress,...) code:
!write (unit=unitNr,fmt='(a)', iostat=oErr) ' '
!write (unit=unitNr,fmt='(5a)', iostat=oErr) label,' in ', fromWhere,',', message
!if (oErr /= 0) then
!   print *, 'ERROR  ',oErr,' while writing log, ', fromWhere, message
!end if
          CALL m_flush(unitNr)
       END IF
    END IF
  END SUBROUTINE cp_logger_log

! *************************** cp_logger_type settings ***************************

!!****f* cp_log_handling/cp_logger_set_log_level [1.0] *
!!
!!   NAME
!!     cp_logger_set_log_level
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_set_log_level(logger, level)
!!       Type(cp_logger_type), Pointer:: logger
!!       Integer, Intent (IN):: level
!!     End Subroutine cp_logger_set_log_level
!!
!!   FUNCTION
!!     changes the logging level. Log messages with a level less than the one
!!     given wo not be printed.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger to change
!!     - level : the new logging level for the logger
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_logger_set_log_level(logger, level)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER, INTENT(in)                      :: level

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_set_log_level', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(logger)) CALL stop_program(routineP,&
         "logger not associated")
    IF (logger%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")
    logger%print_level=level
  END SUBROUTINE cp_logger_set_log_level

!!****f* cp_log_handling/cp_logger_get_log_level [1.0] *
!!
!!   NAME
!!     cp_logger_get_log_level
!!
!!   FUNCTION
!!     returns the logging level of the given logger.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger you want info from
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
!pure
  FUNCTION cp_logger_get_log_level(logger) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_logger_get_log_level', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: lggr

    lggr => logger
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")

    res=lggr%print_level
  END FUNCTION cp_logger_get_log_level

!!****f* cp_log_handling/cp_logger_get_default_unit_nr [1.0] *
!!
!!   NAME
!!     cp_logger_get_default_unit_nr
!!
!!   SYNOPSIS
!!     Function cp_logger_get_default_unit_nr(logger, local) Result(res)
!!       Integer:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Logical, Intent (IN), Optional:: local
!!     End Function cp_logger_get_default_unit_nr
!!
!!   FUNCTION
!!     asks the default unit number of the given logger.
!!     try to use cp_logger_get_unit_nr
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger you want info from
!!     - local: if you want the local unit nr (defaults to false)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
!pure 
  FUNCTION cp_logger_get_default_unit_nr(logger, local) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL, INTENT(in), OPTIONAL            :: local
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_logger_get_default_unit_nr', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=80)                        :: filename, host_name
    INTEGER                                  :: iostat, pid
    LOGICAL                                  :: exists, loc
    TYPE(cp_logger_type), POINTER            :: lggr

    loc=.FALSE.

    lggr => logger
    IF (.NOT.ASSOCIATED(lggr)) lggr => cp_get_default_logger()
    IF (lggr%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")

    IF (PRESENT(local)) loc=local
    IF (.NOT.loc) THEN
       IF (lggr%default_global_unit_nr<=0) THEN
          IF (lggr%para_env%mepos==lggr%para_env%source) THEN
             CALL cp_logger_generate_filename(lggr,filename,lggr%global_filename,&
                  ".out",local=.FALSE.)
             INQUIRE (FILE=TRIM(filename),EXIST=exists)
             IF (exists) THEN
                PRINT *," WARNING ", routineP, " overwriting file ",filename
             END IF
             CALL open_file(TRIM(filename),file_status="unknown",&
                  file_action="WRITE",&
                  unit_number=lggr%default_global_unit_nr)
          END IF
          ! bcast could be dangerous if called only in source processor
          ! thus first global logging upon logger creation
          CALL mp_bcast(lggr%default_global_unit_nr,lggr%para_env%source,&
               lggr%para_env%group)
       END IF
       res=lggr%default_global_unit_nr
    ELSE
       IF (lggr%default_local_unit_nr<0) THEN
          CALL cp_logger_generate_filename(lggr,filename,lggr%local_filename,&
               ".out",local=.TRUE.)
          INQUIRE (FILE=TRIM(filename),EXIST=exists)
          IF (exists) THEN
             PRINT *," WARNING ", routineP, " overwriting file ",filename
          END IF
          CALL open_file(TRIM(filename),file_status="unknown",&
               file_action="WRITE",&
               unit_number=lggr%default_local_unit_nr)
          CALL m_getpid(pid)
          CALL m_hostnm(host_name)
          WRITE (unit=lggr%default_local_unit_nr,fmt='(a,i8,a,i8)',iostat=iostat)&
               ' **** Local logger file of mpi task',lggr%para_env%mepos,&
               ' in communicator',lggr%para_env%group
          IF (iostat==0) THEN
             WRITE (unit=lggr%default_local_unit_nr,fmt='(a,i8,a,a)',iostat=iostat)&
                  ' **** pid=',pid,' host=',host_name
          END IF
          IF (iostat /= 0) THEN
             PRINT *," Error trying to WRITE to the local logger."
             PRINT *," mpi_id=",lggr%para_env%mepos,&
                  " mpi_communicator=",lggr%para_env%group,&
                  " pid=",pid,'host=',host_name
          END IF
          
       END IF
       res=lggr%default_local_unit_nr
    END IF
  END FUNCTION cp_logger_get_default_unit_nr

!!****f* cp_log_handling/my_cp_para_env_release [1.0] *
!!
!!   NAME
!!     my_cp_para_env_release
!!
!!   FUNCTION
!!     releases the para object (to be called when you don't want anymore
!!     the shared copy of this object)
!!
!!   NOTES
!!     local private version to avoid circular dependencies
!!     (see cp_para_env:cp_para_env_release)!
!!
!!   INPUTS
!!     - para_env: the new group
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE my_cp_para_env_release(para_env)
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'my_cp_para_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat

    IF (ASSOCIATED(para_env)) THEN
       IF (para_env%ref_count<1) THEN
          CALL stop_program(routineP,"para_env%ref_count<1")
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
!***************************************************************************

!!****f* cp_log_handling/cp_logger_generate_filename [1.0] *
!!
!!   NAME
!!     cp_logger_generate_filename
!!
!!   SYNOPSIS
!!     Subroutine cp_logger_generate_filename(logger, res, root, postfix,&
!!         local)
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (INOUT):: res
!!       Character(Len=*), Intent (IN):: root, postfix
!!       Logical, Intent (IN), Optional:: local
!!     End Subroutine cp_logger_generate_filename
!!
!!   FUNCTION
!!     generates a unique filename (ie adding eventual suffixes and 
!!     process ids)
!!
!!   NOTES
!!     this should be a function returning a variable length string.
!!     All spaces are moved to the end of the string.
!!     Not fully optimized: result must be a littel longer than the 
!!     resulting compressed filename
!!
!!   INPUTS
!!     - result: the resulting string
!!     - root: the start of filename
!!     - postfix: the end of the name
!!     - local: if the name should be local to this task (defaults to false)
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
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
    IF (lggr%ref_count<1) CALL stop_program(routineP,&
         "logger%ref_count<1")
    IF (PRESENT(local)) loc=local
    IF (loc) THEN
       res=TRIM(root)//TRIM(lggr%suffix)//'_p'//&
            cp_to_string(lggr%para_env%mepos)//postfix
    ELSE
       res=TRIM(root)//TRIM(lggr%suffix)//postfix
    END IF
    CALL compress(res,full=.TRUE.)
  END SUBROUTINE cp_logger_generate_filename
!***************************************************************************

! ============= to string  ============

!!****f* cp_log_handling/cp_int_to_string [1.0] *
!!
!!   NAME
!!     cp_int_to_string
!!
!!   SYNOPSIS
!!     Function cp_int_to_string(i) Result(res)
!!       Character(Len=6):: res
!!       Integer, Intent (IN):: i
!!     End Function cp_int_to_string
!!
!!   FUNCTION
!!     converts an int to a string
!!     (should be a variable length string, but that does not work with 
!!     all the compilers)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     i: the integer to convert
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure
  FUNCTION cp_int_to_string(i) RESULT(res)
    INTEGER, INTENT(in)                      :: i
    CHARACTER(len=6)                         :: res

    INTEGER                                  :: iostat
    REAL(KIND=dp)                            :: tmp_r

    IF (i>999999 .OR. i<-99999) THEN
       tmp_r=i
       WRITE (res,fmt='(es6.1)',iostat=iostat) tmp_r
    ELSE
       WRITE (res,fmt='(i6)',iostat=iostat) i
    END IF
    IF (iostat/=0) THEN
       PRINT *,"cp_int_to_string ioerror",iostat
       CALL print_stack()
    END IF
  END FUNCTION cp_int_to_string

!!****f* cp_log_handling/cp_real_dp_to_string [1.0] *
!!
!!   NAME
!!     cp_real_dp_to_string
!!
!!   SYNOPSIS
!!     Function cp_real_dp_to_string(val) Result(res)
!!       Character(Len=10):: res
!!       Real(Kind=kind(1.0D0)), Intent (IN):: val
!!     End Function cp_real_dp_to_string
!!
!!   FUNCTION
!!     convert a double precision real in a string
!!     (should be a variable length string, but that does not work with 
!!     all the compilers)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     val: the number to convert
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 revised [fawzi]
!!
!!*** *********************************************************************
! pure 
  FUNCTION cp_real_dp_to_string(val) RESULT(res)
    REAL(KIND=dp), INTENT(in)                :: val
    CHARACTER(len=11)                        :: res

    WRITE (res,'(es11.4)') val
  END FUNCTION cp_real_dp_to_string

END MODULE cp_log_handling

!******************************************************************************
