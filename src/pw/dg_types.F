!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      none
! *****************************************************************************
MODULE dg_types

  USE dg_rho0_types,                   ONLY: dg_rho0_create,&
                                             dg_rho0_release,&
                                             dg_rho0_retain,&
                                             dg_rho0_type
  USE f77_blas
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

! Information on the assignment function for Ewald
! *****************************************************************************
  TYPE dg_type
    PRIVATE
    INTEGER :: ref_count, id_nr
    INTEGER :: grid_index
    TYPE ( dg_rho0_type ), POINTER :: dg_rho0
  END TYPE dg_type

! *****************************************************************************
  TYPE dg_p_type
    TYPE ( dg_type ), POINTER :: dg
  END TYPE dg_p_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dg_types'
  INTEGER, SAVE, PRIVATE :: last_dg_id=0

  PRIVATE
  PUBLIC :: dg_type, dg_get,&
            dg_set, dg_retain, dg_release,&
            dg_create

CONTAINS

! *****************************************************************************
!> \brief   Get the dg_type
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dg_get ( dg, id_nr, grid_index, dg_rho0 )
    TYPE(dg_type), POINTER                   :: dg
    INTEGER, OPTIONAL                        :: id_nr, grid_index
    TYPE(dg_rho0_type), OPTIONAL, POINTER    :: dg_rho0

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_get', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( id_nr ) ) id_nr = dg % id_nr
    IF ( PRESENT ( grid_index ) ) grid_index = dg % grid_index
    IF ( PRESENT ( dg_rho0 ) ) dg_rho0 => dg % dg_rho0

  END SUBROUTINE dg_get

! *****************************************************************************
!> \brief   create the dg structure
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dg_create ( dg, error )
    TYPE(dg_type), POINTER                   :: dg
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(dg_rho0_type), POINTER              :: dg_rho0

    ALLOCATE ( dg, stat=istat )
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    NULLIFY ( dg_rho0 )
    CALL dg_rho0_create ( dg_rho0, error )
    dg % dg_rho0 => dg_rho0
    last_dg_id=last_dg_id+1
    dg%id_nr=last_dg_id
    dg%ref_count=1

  END SUBROUTINE dg_create

! *****************************************************************************
!> \brief retains the given dg_type
!> \param dg_type the dg_type to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE dg_retain ( dg, error)
    TYPE(dg_type), POINTER                   :: dg
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(dg),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(dg%ref_count>0,cp_failure_level,routineP,error)
       dg%ref_count=dg%ref_count+1
    END IF
  END SUBROUTINE dg_retain

! *****************************************************************************
!> \brief releases the given dg_type
!> \param dg the dg_type to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      see doc/ReferenceCounting.html
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE dg_release(dg, error)
    TYPE(dg_type), POINTER                   :: dg
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(dg)) THEN
       CPPreconditionNoFail(dg%ref_count>0,cp_failure_level,routineP,error)
       dg%ref_count=dg%ref_count-1
       IF (dg%ref_count==0) THEN
          CALL dg_rho0_release ( dg % dg_rho0, error = error )
          DEALLOCATE (  dg, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(dg)
  END SUBROUTINE dg_release

! *****************************************************************************
!> \brief   Set the double grid environment
!> \version 1.0
! *****************************************************************************
  SUBROUTINE dg_set ( dg, dg_rho0, grid_index, error )
    TYPE(dg_type), POINTER                   :: dg
    TYPE(dg_rho0_type), OPTIONAL, POINTER    :: dg_rho0
    INTEGER, OPTIONAL                        :: grid_index
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_set', &
      routineP = moduleN//':'//routineN

    IF ( PRESENT ( dg_rho0 ) ) THEN
       CALL dg_rho0_retain ( dg_rho0 ,error=error)
       CALL dg_rho0_release ( dg % dg_rho0 ,error=error)
       dg % dg_rho0 => dg_rho0
    END IF
    IF ( PRESENT ( grid_index ) ) dg % grid_index = grid_index
  END SUBROUTINE dg_set

END MODULE dg_types
