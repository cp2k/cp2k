!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Holds information on atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! *****************************************************************************
MODULE atprop_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "f77_int.h"

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atprop_types'

  PUBLIC :: atprop_type, atprop_create, atprop_release, atprop_init

! *****************************************************************************
!> \brief type for the atomic properties
!> \author fawzi
! *****************************************************************************
TYPE atprop_type
   LOGICAL                                   :: energy
   LOGICAL                                   :: stress
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atener
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atecc
   REAL(KIND=dp), DIMENSION(:,:,:), POINTER  :: atstress
END TYPE atprop_type

CONTAINS

! *****************************************************************************
  SUBROUTINE atprop_create(atprop_env,error)
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  CALL atprop_release(atprop_env,error)
  ALLOCATE(atprop_env,stat=stat)
  CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  NULLIFY(atprop_env%atener,atprop_env%atstress,atprop_env%atecc)
  atprop_env%energy = .FALSE.
  atprop_env%stress = .FALSE.

  END SUBROUTINE atprop_create

  SUBROUTINE atprop_init(atprop_env,natom,error)
    TYPE(atprop_type), POINTER               :: atprop_env
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(atprop_env),cp_warning_level,routineP,error,failure)

  IF(atprop_env%energy) THEN
    IF(ASSOCIATED(atprop_env%atener)) THEN
      CPPrecondition(SIZE(atprop_env%atener)==natom,cp_warning_level,routineP,error,failure)
    ELSE
      ALLOCATE(atprop_env%atener(natom),stat=stat)
      CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    atprop_env%atener = 0._dp
  END IF
  NULLIFY(atprop_env%atecc)

  IF(atprop_env%stress) THEN
    IF(ASSOCIATED(atprop_env%atstress)) THEN
      CPPrecondition(SIZE(atprop_env%atstress,3)==natom,cp_warning_level,routineP,error,failure)
    ELSE
      ALLOCATE(atprop_env%atstress(3,3,natom),stat=stat)
      CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    atprop_env%atstress = 0._dp
  END IF

  END SUBROUTINE atprop_init

! *****************************************************************************
!> \brief releases the atprop
!> \param atprop_env the object to release
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author fawzi
! *****************************************************************************
SUBROUTINE atprop_release(atprop_env,error)
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(atprop_env)) THEN
     IF (ASSOCIATED(atprop_env%atener)) THEN
        DEALLOCATE(atprop_env%atener,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
     IF (ASSOCIATED(atprop_env%atecc)) THEN
        DEALLOCATE(atprop_env%atecc,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
     IF (ASSOCIATED(atprop_env%atstress)) THEN
        DEALLOCATE(atprop_env%atstress,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
     DEALLOCATE(atprop_env,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  END IF
  NULLIFY(atprop_env)
END SUBROUTINE atprop_release

END MODULE atprop_types
