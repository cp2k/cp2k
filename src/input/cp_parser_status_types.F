!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow the storage of the parser status
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 08.2008
! *****************************************************************************
MODULE cp_parser_status_types
  USE cp_parser_buffer_types,          ONLY: buffer_type,&
                                             create_buffer_type,&
                                             release_buffer_type
  USE kinds,                           ONLY: max_line_length
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  TYPE status_type
     LOGICAL                                        :: in_use
     INTEGER                                        :: old_input_line_number
     INTEGER                                        :: old_icol
     INTEGER                                        :: old_icol1
     INTEGER                                        :: old_icol2
     CHARACTER(LEN=max_line_length)                 :: old_input_line
     ! Store status of the buffer
     TYPE(buffer_type), POINTER                     :: buffer
  END TYPE status_type

  PUBLIC :: status_type, create_status_type, release_status_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_status_types'

CONTAINS

! ****************************************************************************
!> \brief creates the parser status type
!> \param status ...
!> \param error ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_status_type(status, error)
    TYPE(status_type), POINTER               :: status
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_status_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(status),cp_failure_level,routineP,error,failure)
    ALLOCATE(status, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    status%in_use                = .FALSE.
    status%old_input_line        = ""
    status%old_input_line_number = HUGE(0)
    status%old_icol              = HUGE(0)
    status%old_icol1             = HUGE(0)
    status%old_icol2             = HUGE(0)
    NULLIFY(status%buffer)
    CALL create_buffer_type(status%buffer, error)
  END SUBROUTINE create_status_type

! ****************************************************************************
!> \brief releases the parser status type
!> \param status ...
!> \param error ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE release_status_type(status, error)
    TYPE(status_type), POINTER               :: status
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_status_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(status),cp_failure_level,routineP,error,failure)
    CALL release_buffer_type(status%buffer, error)
    DEALLOCATE(status, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_status_type

END MODULE cp_parser_status_types
