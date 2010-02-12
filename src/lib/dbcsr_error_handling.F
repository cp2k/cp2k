!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   
!> \author VW
!>
!> <b>Modification history:</b>
!> - Created Feb 2010
! *****************************************************************************
MODULE dbcsr_error_handling

  USE dbcsr_kinds,                     ONLY: default_string_length

  IMPLICIT NONE
  PRIVATE

  ! procedures
  PUBLIC :: dbcsr_error_set, dbcsr_error_stop, dbcsr_error_print_stack
  ! types
  PUBLIC :: dbcsr_error_stack_type
  PUBLIC :: dbcsr_error_type
  ! parameters
  PUBLIC :: dbcsr_error_stack_size

  INTEGER, PARAMETER :: dbcsr_error_stack_size = 10

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_error_handling'

! *****************************************************************************
!> \brief stack fo the error type (DBCSR library)
! *****************************************************************************
  TYPE dbcsr_error_stack_type
     INTEGER :: level = 0
     CHARACTER(default_string_length) :: pre = ''
     CHARACTER(default_string_length) :: post = ''
     CHARACTER(default_string_length), DIMENSION(dbcsr_error_stack_size) :: s = ''
     INTEGER, DIMENSION(dbcsr_error_stack_size) :: handler = 0
  END TYPE dbcsr_error_stack_type

! *****************************************************************************
!> \brief error type for the DBCSR library
! *****************************************************************************
  TYPE dbcsr_error_type
     INTEGER :: ierror = 0
     INTEGER :: line = 0
     TYPE(dbcsr_error_stack_type) :: stack
  END TYPE dbcsr_error_type

CONTAINS

! *****************************************************************************
!> \brief set the error
!> \param[in]    name    name of the routine
!> \param[inout] error   dbcsr error
! *****************************************************************************
  SUBROUTINE dbcsr_error_set(name, handler, error)
    CHARACTER(len=*), INTENT(in)             :: name
    INTEGER, INTENT(out)                     :: handler
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    error%stack%level = error%stack%level + 1
    IF(error%stack%level.GT.dbcsr_error_stack_size.OR.&
       error%stack%level.LT.1) STOP
    error%stack%s(error%stack%level) = name
    handler = error%stack%level
    error%stack%handler(error%stack%level) = handler
  END SUBROUTINE dbcsr_error_set

! *****************************************************************************
!> \brief stop the error
!> \param[inout] error   dbcsr error
! *****************************************************************************
  SUBROUTINE dbcsr_error_stop(handler, error)
    INTEGER, INTENT(inout)                   :: handler
    TYPE(dbcsr_error_type)                   :: error

    IF(error%stack%level.GT.dbcsr_error_stack_size.OR.&
       error%stack%level.LT.1) THEN
       WRITE(*,*) 'level=',error%stack%level
       WRITE(*,*) 'dbcsr_error_stack_size=',dbcsr_error_stack_size
       STOP
    ENDIF
    IF(handler.NE.error%stack%handler(error%stack%level)) THEN
       WRITE(*,*) 'get handler=',handler
       WRITE(*,*) 'expected handler=',error%stack%handler(error%stack%level)
       WRITE(*,*) 'you may have forget to call dbcsr_error_stop (may be a RETURN)...'
       STOP
    ENDIF
    error%stack%s(error%stack%level) = ''
    error%stack%handler(error%stack%level) = 0
    error%stack%level = error%stack%level - 1
    ! destroy the handler for any case
    handler = -HUGE(0)
  END SUBROUTINE dbcsr_error_stop

! *****************************************************************************
!> \brief stop the error
!> \param[inout] error   dbcsr error
! *****************************************************************************
  SUBROUTINE dbcsr_error_print_stack(error, output_unit)
    TYPE(dbcsr_error_type)                   :: error
    INTEGER, INTENT(in)                      :: output_unit

    INTEGER                                  :: ilevel

    IF(error%stack%level.GT.dbcsr_error_stack_size.OR.&
       error%stack%level.LT.1) THEN
       WRITE(*,*) 'level=',error%stack%level
       WRITE(*,*) 'dbcsr_error_stack_size=',dbcsr_error_stack_size
       STOP
    ENDIF
    IF(LEN_TRIM(error%stack%pre).GT.0) &
         WRITE(output_unit,'(T3,A,A)') 'prestack: ',TRIM(error%stack%pre)
    WRITE(output_unit,'(T3,A)') 'stack: '
    WRITE(output_unit,'(T3,I3,3A,I5)') error%stack%level,' error in ',&
         TRIM(error%stack%s(error%stack%level)),' at line ',error%line
    DO ilevel = error%stack%level-1,1,-1
       WRITE(output_unit,'(T3,I3,A,A)') ilevel,' called from ',TRIM(error%stack%s(ilevel))
    ENDDO
    IF(LEN_TRIM(error%stack%post).GT.0) &
         WRITE(output_unit,'(T3,A,A)') 'poststack: ',TRIM(error%stack%post)
  END SUBROUTINE dbcsr_error_print_stack

END MODULE dbcsr_error_handling
