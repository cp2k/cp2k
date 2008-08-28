!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
 
! *****************************************************************************
!> \brief a module to allow simple buffering of read lines of a parser
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 08.2008
! *****************************************************************************
MODULE cp_parser_buffer_types
  USE kinds,                           ONLY: max_line_length
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

! ****************************************************************************
!> \brief  Buffer type for speeding-up the parsing in parallel
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  TYPE buffer_type
     INTEGER                              :: size
     INTEGER                              :: present_line_number,&
                                             last_line_number,&
                                             istat
     INTEGER, DIMENSION(:), POINTER       :: input_line_numbers
     CHARACTER(LEN=max_line_length), &
          DIMENSION(:), POINTER           :: input_lines
  END TYPE buffer_type

  PUBLIC :: buffer_type, create_buffer_type, release_buffer_type, copy_buffer_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_buffer_types'
  INTEGER, PARAMETER, PRIVATE          :: buffer_size=1000

CONTAINS

! ****************************************************************************
!> \brief  Creates the parser buffer type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  SUBROUTINE create_buffer_type(buffer, error)
    TYPE(buffer_type), POINTER               :: buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_buffer_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(buffer),cp_failure_level,routineP,error,failure)
    ALLOCATE(buffer, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    buffer%size = buffer_size
    ALLOCATE(buffer%input_lines(buffer%size), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(buffer%input_line_numbers(buffer%size), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    buffer%input_line_numbers  = 0
    buffer%istat               = 0
    buffer%present_line_number = buffer%size
    buffer%last_line_number    = buffer%size
  END SUBROUTINE create_buffer_type

! ****************************************************************************
!> \brief  Releases the parser buffer type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  SUBROUTINE release_buffer_type(buffer, error)
    TYPE(buffer_type), POINTER               :: buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_buffer_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(buffer),cp_failure_level,routineP,error,failure)
    DEALLOCATE(buffer%input_lines, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(buffer%input_line_numbers, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(buffer, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_buffer_type

! ****************************************************************************
!> \brief  copy buffer types
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  SUBROUTINE copy_buffer_type(buffer_in, buffer_out, error)
    TYPE(buffer_type), POINTER               :: buffer_in, buffer_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_buffer_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(buffer_in),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(buffer_out),cp_failure_level,routineP,error,failure)
    CPPostcondition(buffer_in%size==buffer_out%size,cp_failure_level,routineP,error,failure)
    buffer_out%present_line_number = buffer_in%present_line_number
    buffer_out%last_line_number    = buffer_in%last_line_number
    buffer_out%istat               = buffer_in%istat
    buffer_out%input_line_numbers  = buffer_in%input_line_numbers
    buffer_out%input_lines         = buffer_in%input_lines        
  END SUBROUTINE copy_buffer_type

END MODULE cp_parser_buffer_types
