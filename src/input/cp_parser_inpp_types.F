!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a module to allow simple internal preprocessing in input files.
!> \par History
!>      - standalone proof-of-concept implementation (20.02.2008,AK)
!>      - integration into cp2k (22.02.2008,tlaino)
!>      - variables added (25.02.2008,AK)
!> \author Axel Kohlmeyer [AK] - CMM/UPenn Philadelphia
!> \date 25.02.2008
! *****************************************************************************
MODULE cp_parser_inpp_types
  
  USE kinds,                           ONLY: default_path_length
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  TYPE inpp_type
     ! for '@INCLUDE "some_file.inc"'
     ! currently open include file stack pointer
     INTEGER                              :: io_stack_level
     ! include file stack data
     INTEGER, POINTER, DIMENSION(:)       :: io_stack_channel,&
                                             io_stack_lineno
     CHARACTER (len=default_path_length),&
        POINTER, DIMENSION(:)             :: io_stack_filename
     ! for '@SET VAR value' and '${VAR}'
     ! table size
     INTEGER                              :: num_variables
     ! table entries
     CHARACTER (len=default_path_length), &
          POINTER, DIMENSION(:)       :: variable_name
     CHARACTER (len=default_path_length), &
          POINTER, DIMENSION(:)       :: variable_value
  END TYPE inpp_type

  PUBLIC :: inpp_type, create_inpp_type, release_inpp_type
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_parser_inpp_types'

CONTAINS

! ****************************************************************************
!> \brief creates the internal preprocessing type
!> \param inpp ...
!> \param initial_variables ...
!> \param error ...
!> \date  22.02.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE create_inpp_type(inpp, initial_variables, error)
    TYPE(inpp_type), POINTER                 :: inpp
    CHARACTER(len=default_path_length), &
      DIMENSION(:, :), POINTER               :: initial_variables
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

    inpp%num_variables =0
    NULLIFY(inpp%variable_name,&
            inpp%variable_value)

    IF (ASSOCIATED(initial_variables)) THEN
       inpp%num_variables =SIZE(initial_variables,2)
       ALLOCATE(inpp%variable_name(inpp%num_variables),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       inpp%variable_name=initial_variables(1,:)
       ALLOCATE(inpp%variable_value(inpp%num_variables),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       inpp%variable_value=initial_variables(2,:)
    ENDIF

  END SUBROUTINE create_inpp_type

! ****************************************************************************
!> \brief releases the internal preprocessing type
!> \param inpp ...
!> \param error ...
!> \date  22.02.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
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

    IF (ASSOCIATED(inpp%variable_name)) THEN
       DEALLOCATE(inpp%variable_name,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (ASSOCIATED(inpp%variable_value)) THEN
       DEALLOCATE(inpp%variable_value,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(inpp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE release_inpp_type

END MODULE cp_parser_inpp_types
