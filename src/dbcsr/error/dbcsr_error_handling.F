!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief
!> \author VW
!>
!> <b>Modification history:</b>
!> - Created Feb 2010
! *****************************************************************************
MODULE dbcsr_error_handling

#include "../../base/base_uses.f90"
  IMPLICIT NONE
  PRIVATE

  ! procedures
  PUBLIC :: dbcsr_error_set, dbcsr_error_stop
  PUBLIC :: dbcsr_assert, dbcsr_abort
  ! types
  PUBLIC :: dbcsr_error_type
  ! parameters
  PUBLIC :: dbcsr_warning_level, dbcsr_failure_level, dbcsr_fatal_level
  PUBLIC :: dbcsr_caller_error, dbcsr_wrong_args_error,&
            dbcsr_internal_error,&
            dbcsr_unimplemented_error_nr

  ! interfaces
  INTERFACE dbcsr_assert
     MODULE PROCEDURE dbcsr_int_assert
     MODULE PROCEDURE dbcsr_logical_assert
     MODULE PROCEDURE dbcsr_not_assert
     MODULE PROCEDURE dbcsr_char_assert
     MODULE PROCEDURE dbcsr_true_assert!bcsr_assert,
  END INTERFACE

  INTEGER, PARAMETER :: dbcsr_error_stack_size = 30
  !! level of an error
  INTEGER, PARAMETER :: dbcsr_fatal_level=3
  !! level of a failure
  INTEGER, PARAMETER :: dbcsr_failure_level=2
  !! level of a warning
  INTEGER, PARAMETER :: dbcsr_warning_level=1
  !! level of a note
  INTEGER, PARAMETER :: dbcsr_note_level=0
  !! error number: no error
  INTEGER, PARAMETER :: dbcsr_no_error = 0
  !! error number: generic error on the side of the caller
  INTEGER, PARAMETER :: dbcsr_caller_error = 1
  !! error number: one or more arguments have and invalid value
  INTEGER, PARAMETER :: dbcsr_wrong_args_error = 100
  !! error number: precondition failed
  INTEGER, PARAMETER :: dbcsr_precondition_failed = 200
  !! error number: generic error inside the routine
  INTEGER, PARAMETER :: dbcsr_internal_error = -1
  !! error number: postcondition failed
  INTEGER, PARAMETER :: dbcsr_postcondition_failed = -200
  !! error number: invariant failed
  INTEGER, PARAMETER :: dbcsr_invariant_failed = -100
  !! error number: assertion failure
  INTEGER, PARAMETER :: dbcsr_assertion_failed = -300
  !! error number: not implemented
  INTEGER, PARAMETER :: dbcsr_unimplemented_error_nr = -1000


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_error_handling'

! *****************************************************************************
!> \brief error type for the DBCSR library
! *****************************************************************************
  TYPE dbcsr_error_type
     PRIVATE
     INTEGER :: dummy = 42 ! empty types confuse some compilers
  END TYPE dbcsr_error_type

CONTAINS

! *****************************************************************************
!> \brief set the error
!> \param[in]    name    name of the routine
!> \param handler ...
!> \param[inout] error   dbcsr error
! *****************************************************************************
  SUBROUTINE dbcsr_error_set(name, handler, error)
    CHARACTER(len=*), INTENT(in)             :: name
    INTEGER, INTENT(out)                     :: handler
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CALL timeset(name, handler)

  END SUBROUTINE dbcsr_error_set

! *****************************************************************************
!> \brief stop the error
!> \param handler ...
!> \param[inout] error   dbcsr error
! *****************************************************************************
  SUBROUTINE dbcsr_error_stop(handler, error)
    INTEGER, INTENT(inout)                   :: handler
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CALL timestop(handler)

  END SUBROUTINE dbcsr_error_stop


! *****************************************************************************
!> \brief ...
!> \param routine ...
!> \param line ...
!> \param msg ...
!> \param error ...
!> \param[inout]
! *****************************************************************************
  SUBROUTINE dbcsr_abort(routine, line,  msg, error)
    CHARACTER(*), INTENT(in)                 :: routine
    INTEGER, INTENT(in)                      :: line
    CHARACTER(*), INTENT(in)                 :: msg
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CALL dbcsr_abort_low(dbcsr_fatal_level, dbcsr_internal_error, routine, msg, line, error)

  END SUBROUTINE dbcsr_abort

! *****************************************************************************
!> \brief ...
!> \param level ...
!> \param etype ...
!> \param routine ...
!> \param msg ...
!> \param line ...
!> \param error ...
!> \param[inout]
! *****************************************************************************
  SUBROUTINE dbcsr_abort_low (level, etype, routine, msg, line, error)
    INTEGER, INTENT(in)                      :: level, etype
    CHARACTER(*), INTENT(in)                 :: routine, msg
    INTEGER, INTENT(in)                      :: line
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    IF(level==dbcsr_warning_level) THEN
       CALL cp_warn(cp__l(routine, line), msg)
    ELSE
       CALL cp_abort(cp__l(routine, line), msg)
    ENDIF

  END SUBROUTINE dbcsr_abort_low

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_int_assert(left, rel, right, level, etype, routine, msg, line, error)
    INTEGER, INTENT(IN)                      :: left
    CHARACTER(len=2), INTENT(IN)             :: rel
    INTEGER, INTENT(IN)                      :: right, level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQ")
       l = left .EQ. right
    CASE ("LT")
       l = left .LT. right
    CASE ("LE")
       l = left .LE. right
    CASE ("GT")
       l = left .GT. right
    CASE ("GE")
       l = left .GE. right
    CASE ("NE")
       l = left .NE. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel, __LINE__, error)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,I9,A4,I9)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line, error)
    ENDIF
  END SUBROUTINE dbcsr_int_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_logical_assert(left, rel, right, level, etype, routine, msg, line, error)
    LOGICAL, INTENT(IN)                      :: left
    CHARACTER(len=*), INTENT(IN)             :: rel
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQV")
       l = left .EQV. right
    CASE ("NEQV")
       l = left .NEQV. right
    CASE ("OR")
       l = left .OR. right
    CASE ("AND")
       l = left .AND. right
    CASE ("IMP")
       l = .NOT. left .OR. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel, __LINE__, error)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,L1,A,L1)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line, error)
    ENDIF
  END SUBROUTINE dbcsr_logical_assert

! *****************************************************************************
!> \brief Assertion
!> \param[in] left            left value
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_char_assert(left, rel, right, level, etype, routine, msg, line, error)
    CHARACTER, INTENT(IN)                    :: left
    CHARACTER(len=2), INTENT(IN)             :: rel
    CHARACTER, INTENT(IN)                    :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("EQ")
       l = left .EQ. right
    CASE ("NE")
       l = left .NE. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_char_assert", "Invalid relation specified: "//rel, __LINE__, error)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,A1,1X,A4,1X,A1)')"ASSERTION FAILED:",&
            left, "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line, error)
    ENDIF
  END SUBROUTINE dbcsr_char_assert


! *****************************************************************************
!> \brief Assertion
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_true_assert(right, level, etype, routine, msg, line, error)
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    l = right
    IF (.NOT.l) THEN
       WRITE(*,'(1X,A,1X,L1)')"ASSERTION FAILED:",&
            right
       CALL dbcsr_abort_low (level, etype, routine, msg, line, error)
    ENDIF
    !CALL cp_assert(l, level, etype, routine, msg)
  END SUBROUTINE dbcsr_true_assert


! *****************************************************************************
!> \brief Assertion
!> \param[in] rel             relation
!> \param[in] right           right value
!> \param[in] level           error level
!> \param[in] etype           error type
!> \param[in] routine         Routine name
!> \param[in] msg   Message to display if the assertion fails
!> \param line ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_not_assert(rel, right, level, etype, routine, msg, line, error)
    CHARACTER(len=3), INTENT(IN)             :: rel
    LOGICAL, INTENT(IN)                      :: right
    INTEGER, INTENT(IN)                      :: level, etype
    CHARACTER(len=*), INTENT(IN)             :: routine, msg
    INTEGER, INTENT(IN)                      :: line
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    LOGICAL                                  :: l

!   ---------------------------------------------------------------------------

    SELECT CASE (rel)
    CASE ("NOT")
       l = .NOT. right
    CASE default
       CALL dbcsr_assert (.FALSE., dbcsr_fatal_level, dbcsr_wrong_args_error,&
            "dbcsr_int_assert", "Invalid relation specified: "//rel, __LINE__, error)
       l = .FALSE.
    END SELECT
    IF (.NOT. l) THEN
       WRITE(*,'(1X,A,1X,A,L1)')"ASSERTION FAILED:",&
            "."//rel//".", right
       CALL dbcsr_abort_low (level, etype, routine, msg, line, error)
    ENDIF
  END SUBROUTINE dbcsr_not_assert

END MODULE dbcsr_error_handling
