!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Holds information on atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! *****************************************************************************
MODULE atprop_types
  
  USE kinds,                           ONLY: dp
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atprop_types'

  PUBLIC :: atprop_type, atprop_create, atprop_release, atprop_init
  PUBLIC :: atprop_array_init, atprop_array_add

! *****************************************************************************
!> \brief type for the atomic properties
!> \author fawzi
! *****************************************************************************
TYPE atprop_type
   LOGICAL                                   :: energy
   LOGICAL                                   :: stress
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atener
   REAL(KIND=dp), DIMENSION(:), POINTER      :: ateb
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atexc
   REAL(KIND=dp), DIMENSION(:), POINTER      :: ateself
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atevdw
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atecc
   REAL(KIND=dp), DIMENSION(:), POINTER      :: ate1c
   REAL(KIND=dp), DIMENSION(:), POINTER      :: atecoul
   REAL(KIND=dp), DIMENSION(:,:,:), POINTER  :: atstress
END TYPE atprop_type

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atprop_env ...
!> \param error ...
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
  NULLIFY(atprop_env%atener,atprop_env%atstress)
  NULLIFY(atprop_env%ateb,atprop_env%atevdw,atprop_env%atecc,atprop_env%atecoul)
  NULLIFY(atprop_env%ateself,atprop_env%atexc,atprop_env%ate1c)
  atprop_env%energy = .FALSE.
  atprop_env%stress = .FALSE.

  END SUBROUTINE atprop_create

! *****************************************************************************
!> \brief ...
!> \param atprop_env ...
!> \param natom ...
!> \param error ...
! *****************************************************************************
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
    CALL atprop_array_init(atprop_env%atener,natom,error)
    CALL atprop_array_release(atprop_env%ateb,error)
    CALL atprop_array_release(atprop_env%atevdw,error)
    CALL atprop_array_release(atprop_env%atecc,error)
    CALL atprop_array_release(atprop_env%atecoul,error)
    CALL atprop_array_release(atprop_env%ateself,error)
    CALL atprop_array_release(atprop_env%atexc,error)
    CALL atprop_array_release(atprop_env%ate1c,error)
  END IF

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
!> \brief ...
!> \param atarray ...
!> \param natom ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atprop_array_init(atarray,natom,error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atarray
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_array_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(atarray)) THEN
      CPPrecondition(SIZE(atarray)==natom,cp_warning_level,routineP,error,failure)
    ELSE
      ALLOCATE(atarray(natom),stat=stat)
      CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    atarray = 0._dp

  END SUBROUTINE atprop_array_init

! *****************************************************************************
!> \brief ...
!> \param atarray ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atprop_array_release(atarray,error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: atarray
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_array_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(atarray)) THEN
      DEALLOCATE(atarray,stat=stat)
      CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
      NULLIFY(atarray)
    END IF

  END SUBROUTINE atprop_array_release

! *****************************************************************************
!> \brief ...
!> \param array_a ...
!> \param array_b ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atprop_array_add(array_a,array_b,error)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: array_a, array_b
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atprop_array_add', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF(ASSOCIATED(array_b)) THEN
      CPPrecondition(ASSOCIATED(array_a),cp_warning_level,routineP,error,failure)
      array_a = array_a + array_b
    END IF

  END SUBROUTINE atprop_array_add

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
     ! energy
     CALL atprop_array_release(atprop_env%atener,error)
     CALL atprop_array_release(atprop_env%ateb,error)
     CALL atprop_array_release(atprop_env%ateself,error)
     CALL atprop_array_release(atprop_env%atexc,error)
     CALL atprop_array_release(atprop_env%atevdw,error)
     CALL atprop_array_release(atprop_env%atecc,error)
     CALL atprop_array_release(atprop_env%ate1c,error)
     CALL atprop_array_release(atprop_env%atecoul,error)
     ! stress
     IF (ASSOCIATED(atprop_env%atstress)) THEN
        DEALLOCATE(atprop_env%atstress,stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
     ! atprop type
     DEALLOCATE(atprop_env,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  END IF
  NULLIFY(atprop_env)
END SUBROUTINE atprop_release

END MODULE atprop_types
