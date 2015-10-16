!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module that contains the routines for error handling
!>      The idea is that in general routines have an optional
!>      parameter that can be used to change the error behaviour.
!> \note
!>      In general 0 means no error, a negative error code means an internal
!>      error (in the routine, postcondition failure), and a positive error
!>      an error on the side of the caller (precondition failure).
!>      Error numbers in the dange -1000..1000 are reserved for generic error
!>      codes.
!>      Inspired from the nag nl90 library.
!>      Some help on this is also available in the high level dev notes on pao.
!>
!>      The following error levels are defined:
!>      0 - note    : CP_NOTE_LEVEL
!>      1 - warning : CP_WARNING_LEVEL
!>      2 - failure : CP_FAILURE_LEVEL
!>      3 - fatal   : CP_FATAL_LEVEL
!>      You shouldn't use the numbers but the constants defined in
!>      cp_log_handling
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
MODULE cp_error_handling
  USE base_hooks,                      ONLY: cp_abort_hook,&
                                             cp_warn_hook
  USE cp_log_handling,                 ONLY: CP_WARNING_LEVEL,&
                                             cp_logger_get_default_io_unit
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_abort
  USE termination,                     ONLY: print_message
  USE timings,                         ONLY: print_stack
#include "../base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_error_handling'

  !API parameters var
  PUBLIC :: cp_assertion_failed, cp_internal_error, cp_wrong_args_error,&
            cp_precondition_failed, cp_caller_error

  !API public routines
  PUBLIC :: cp_assert, cp_error_handling_setup

  !! error number: no error
  INTEGER, PARAMETER :: cp_no_error = 0
  !! error number: generic error on the side of the caller
  INTEGER, PARAMETER :: cp_caller_error = 1
  !! error number: one or more arguments have and invalid value
  INTEGER, PARAMETER :: cp_wrong_args_error = 100
  !! error number: precondition failed
  INTEGER, PARAMETER :: cp_precondition_failed = 200
  !! error number: generic error inside the routine
  INTEGER, PARAMETER :: cp_internal_error = -1
  !! error number: postcondition failed
  INTEGER, PARAMETER :: cp_postcondition_failed = -200
  !! error number: invariant failed
  INTEGER, PARAMETER :: cp_invariant_failed = -100
  !! error number: assertion failure
  INTEGER, PARAMETER :: cp_assertion_failed = -300

  !API (via pointer assignment to hook, PR67982, not meant to be called directly)
  PUBLIC :: cp_abort_handler, cp_warn_handler

CONTAINS

! *****************************************************************************
!> \brief Registers handlers with base_hooks.F
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_error_handling_setup()
    cp_abort_hook => cp_abort_handler
    cp_warn_hook => cp_warn_handler
  END SUBROUTINE cp_error_handling_setup

! *****************************************************************************
!> \brief routine that checks an assertion, if it fails set the optional
!>      argument failure is given sets it to true.
!>      This can be used to do many assertion one fater the other and
!>      check at the end if one failed.
!> \param condition the condition to check, if true does nothing
!> \param level the level of the error (cp_note_level, cp_warning_level,
!> \param error_nr the number of the error
!> \param fromWhere a string that contains the module name and routine name
!>            where this test happened
!> \param message the error message
!> \param failure is set to true if the condition is false otherwise
!>            it is not changed (optional)
!> \param only_ionode ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!>      @version 1.2002
! *****************************************************************************
  SUBROUTINE cp_assert(condition, level, error_nr, fromWhere, &
       message,failure,only_ionode)
    !! the condition that is checked, if false, an error did happen
    LOGICAL, INTENT(in)                      :: condition
    INTEGER, INTENT(in)                      :: level, error_nr
    CHARACTER(len=*), INTENT(in)             :: fromWhere, message
    LOGICAL, INTENT(inout), OPTIONAL         :: failure
    LOGICAL, INTENT(in), OPTIONAL            :: only_ionode

    INTEGER                                  :: dummy

    IF (.NOT. condition) THEN
       IF (PRESENT(failure)) failure=.TRUE.
       IF (PRESENT(only_ionode)) dummy = 42 !pretend only_ionode is used.

       IF(level > cp_warning_level) THEN
          CALL cp_abort(fromWhere, message)
       ELSE
          CALL cp_warn(fromWhere, message)
       ENDIF
    ENDIF
  END SUBROUTINE cp_assert

! *****************************************************************************
!> \brief Abort program with error message
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_abort_handler(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    INTEGER                                  :: unit_nr

    unit_nr = cp_logger_get_default_io_unit()
    IF(unit_nr>0) THEN
       CALL print_message("ABORT in "//TRIM(location)//' :: '//TRIM(ADJUSTL(message)),unit_nr,2,1,1)
       CALL m_flush(unit_nr)
       CALL print_stack(unit_nr)
       CALL m_flush(unit_nr)
    ENDIF

    CALL mp_abort(message)
  END SUBROUTINE cp_abort_handler

! *****************************************************************************
!> \brief Signal a warning
!> \param location ...
!> \param message ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE cp_warn_handler(location, message)
    CHARACTER(len=*), INTENT(in)             :: location, message

    INTEGER                                  :: unit_nr

    unit_nr = cp_logger_get_default_io_unit()
    IF(unit_nr>0) THEN
       CALL print_message("WARNING in "//TRIM(location)//' :: '//TRIM(ADJUSTL(message)),unit_nr,1,1,1)
       CALL m_flush(unit_nr)
    ENDIF
  END SUBROUTINE cp_warn_handler

END MODULE cp_error_handling
