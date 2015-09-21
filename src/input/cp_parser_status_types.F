!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_status_type(status)
    TYPE(status_type), POINTER               :: status

    CHARACTER(len=*), PARAMETER :: routineN = 'create_status_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(status))
    ALLOCATE(status)
    status%in_use                = .FALSE.
    status%old_input_line        = ""
    status%old_input_line_number = HUGE(0)
    status%old_icol              = HUGE(0)
    status%old_icol1             = HUGE(0)
    status%old_icol2             = HUGE(0)
    NULLIFY(status%buffer)
    CALL create_buffer_type(status%buffer)
  END SUBROUTINE create_status_type

! ****************************************************************************
!> \brief releases the parser status type
!> \param status ...
!> \date  08.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE release_status_type(status)
    TYPE(status_type), POINTER               :: status

    CHARACTER(len=*), PARAMETER :: routineN = 'release_status_type', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(ASSOCIATED(status))
    CALL release_buffer_type(status%buffer)
    DEALLOCATE(status)
  END SUBROUTINE release_status_type

END MODULE cp_parser_status_types
