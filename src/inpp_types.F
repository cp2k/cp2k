!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
 
! *****************************************************************************
!> \brief a module to allow simple internal preprocessing in input files.
!> \par History
!>      - standalone proof-of-concept implemenation (20.02.2008,AK)
!> \author AK
! *****************************************************************************
MODULE inpp_types

  USE kinds,                           ONLY: default_path_length
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  TYPE inpp_type
     ! current stack pointer
     INTEGER                              :: io_stack_level
     INTEGER, POINTER, DIMENSION(:)       :: io_stack_channel,&
                                             io_stack_lineno
     CHARACTER (len=default_path_length),&
        POINTER, DIMENSION(:)             :: io_stack_filename
  END TYPE inpp_type

  PUBLIC :: inpp_type, create_inpp_type, release_inpp_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'inpp_types'
  
CONTAINS

! ****************************************************************************
!> \brief creates the internal preprocessing type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date  22.02.2008
! *****************************************************************************
  SUBROUTINE create_inpp_type(inpp, error)
    TYPE(inpp_type), POINTER                 :: inpp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_inpp_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(.NOT.ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)
    ALLOCATE(inpp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    inpp%io_stack_level = 0
    NULLIFY(inpp%io_stack_channel,&
            inpp%io_stack_lineno,&
            inpp%io_stack_filename)

  END SUBROUTINE create_inpp_type

! ****************************************************************************
!> \brief creates the internal preprocessing type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date  22.02.2008
! *****************************************************************************
  SUBROUTINE release_inpp_type(inpp, error)
    TYPE(inpp_type), POINTER                 :: inpp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_inpp_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(inpp),cp_failure_level,routineP,error,failure)
    IF (ASSOCIATED(inpp%io_stack_channel)) THEN
       DEALLOCATE(inpp%io_stack_channel,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(inpp%io_stack_lineno)) THEN
       DEALLOCATE(inpp%io_stack_lineno,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(inpp%io_stack_filename)) THEN
       DEALLOCATE(inpp%io_stack_filename,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(inpp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_inpp_type

END MODULE inpp_types
