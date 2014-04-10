!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow simple buffering of read lines of a parser
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 08.2008
! *****************************************************************************
MODULE cp_parser_buffer_types
  USE f77_blas
  USE kinds,                           ONLY: max_line_length
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

! ****************************************************************************
!> \brief  Buffer type for speeding-up the parsing in parallel
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  TYPE buffer_type
     INTEGER                              :: size, buffer_id
     INTEGER                              :: present_line_number,&
                                             last_line_number,&
                                             istat
     INTEGER, DIMENSION(:), POINTER       :: input_line_numbers
     CHARACTER(LEN=max_line_length), &
          DIMENSION(:), POINTER           :: input_lines
     TYPE(buffer_type), POINTER           :: sub_buffer
  END TYPE buffer_type

  PUBLIC :: buffer_type, create_buffer_type, release_buffer_type, copy_buffer_type,&
            initialize_sub_buffer, finalize_sub_buffer
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
    buffer%buffer_id           = 0
    buffer%input_line_numbers  = 0
    buffer%istat               = 0
    buffer%present_line_number = buffer%size
    buffer%last_line_number    = buffer%size
    NULLIFY(buffer%sub_buffer)
  END SUBROUTINE create_buffer_type

! ****************************************************************************
!> \brief  Releases the parser buffer type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE release_buffer_type(buffer, error)
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
    IF (ASSOCIATED(buffer%sub_buffer)) THEN
       CALL release_buffer_type(buffer%sub_buffer, error)
    END IF
    DEALLOCATE(buffer, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_buffer_type

! ****************************************************************************
!> \brief  Copies  buffer types
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE copy_buffer_type(buffer_in, buffer_out, force, error)
    TYPE(buffer_type), POINTER               :: buffer_in, buffer_out
    LOGICAL, INTENT(IN), OPTIONAL            :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'copy_buffer_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_force

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(buffer_in),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(buffer_out),cp_failure_level,routineP,error,failure)
    CPPostcondition(buffer_in%size==buffer_out%size,cp_failure_level,routineP,error,failure)
    my_force = .FALSE.
    IF (PRESENT(force)) my_force = force
    ! Copy buffer structure
    buffer_out%present_line_number = buffer_in%present_line_number
    buffer_out%last_line_number    = buffer_in%last_line_number
    buffer_out%istat               = buffer_in%istat
    ! This part can be quite expensive.. we do it only when strictly necessary..
    IF ((buffer_out%buffer_id/=buffer_in%buffer_id).OR.(my_force)) THEN
       buffer_out%buffer_id           = buffer_in%buffer_id
       buffer_out%input_line_numbers  = buffer_in%input_line_numbers
       ! Explicit loop: bypass a NAG bug..
       DO i = 1, SIZE(buffer_in%input_lines)
          buffer_out%input_lines(i)   = buffer_in%input_lines(i)
       END DO
    END IF
    IF (ASSOCIATED(buffer_in%sub_buffer).AND.ASSOCIATED(buffer_out%sub_buffer)) THEN
       CALL copy_buffer_type(buffer_in%sub_buffer, buffer_out%sub_buffer, force, error)
    END IF
  END SUBROUTINE copy_buffer_type

! ****************************************************************************
!> \brief  Initializes sub buffer structure
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  SUBROUTINE initialize_sub_buffer(sub_buffer, buffer, error)
    TYPE(buffer_type), POINTER               :: sub_buffer, buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_sub_buffer', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(buffer),cp_failure_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(sub_buffer),cp_failure_level,routineP,error,failure)
    CALL create_buffer_type(sub_buffer, error)
    CALL copy_buffer_type(buffer, sub_buffer,error=error)
    sub_buffer%present_line_number = 0
  END SUBROUTINE initialize_sub_buffer


! ****************************************************************************
!> \brief  Finalizes sub buffer structure
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   08.2008
! *****************************************************************************
  SUBROUTINE finalize_sub_buffer(sub_buffer, buffer, error)
    TYPE(buffer_type), POINTER               :: sub_buffer, buffer
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'finalize_sub_buffer', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(buffer),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(sub_buffer),cp_failure_level,routineP,error,failure)
    CALL copy_buffer_type(sub_buffer,buffer,error=error)
    CALL release_buffer_type(sub_buffer, error=error)
  END SUBROUTINE finalize_sub_buffer

END MODULE cp_parser_buffer_types
