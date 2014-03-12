!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types for task lists
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele
! *****************************************************************************
MODULE task_list_types
  USE f77_blas
  USE kinds,                           ONLY: dp,&
                                             int_8
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_types'

! *****************************************************************************
  TYPE task_list_type
    INTEGER(kind=int_8), DIMENSION(:, :), POINTER :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER       :: dist_ab
    INTEGER(kind=int_8), DIMENSION(:), POINTER    :: atom_pair_send, atom_pair_recv
    INTEGER                                       :: ntasks
    INTEGER, DIMENSION(:,:),POINTER               :: taskstart,taskstop
    INTEGER, DIMENSION(:),POINTER                 :: npairs
  END TYPE task_list_type

  PUBLIC :: task_list_type

  PUBLIC :: allocate_task_list,&
            deallocate_task_list

CONTAINS

! *****************************************************************************
!> \brief allocates and initialised the components of the task_list_type
!> \par History
!>      01.2008 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE allocate_task_list(task_list,error)
    TYPE(task_list_type), POINTER            :: task_list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_task_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  ALLOCATE(task_list,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  NULLIFY(task_list%tasks)
  NULLIFY(task_list%dist_ab)
  NULLIFY(task_list%atom_pair_send)
  NULLIFY(task_list%atom_pair_recv)
  NULLIFY(task_list%taskstart)
  NULLIFY(task_list%taskstop)
  NULLIFY(task_list%npairs)
  task_list%ntasks=0
END SUBROUTINE allocate_task_list

! *****************************************************************************
!> \brief deallocates the components and the object itself
!> \par History
!>      01.2008 created [Joost VandeVondele]
! *****************************************************************************
SUBROUTINE deallocate_task_list(task_list,error)
    TYPE(task_list_type), POINTER            :: task_list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_task_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  IF (ASSOCIATED(task_list%tasks)) THEN
     DEALLOCATE(task_list%tasks,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF
  IF (ASSOCIATED(task_list%dist_ab)) THEN
     DEALLOCATE(task_list%dist_ab,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF
  IF (ASSOCIATED(task_list%atom_pair_send)) THEN
     DEALLOCATE(task_list%atom_pair_send,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF
  IF (ASSOCIATED(task_list%atom_pair_recv)) THEN
     DEALLOCATE(task_list%atom_pair_recv,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF
  IF (ASSOCIATED(task_list%taskstart)) THEN
     DEALLOCATE(task_list%taskstart,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF
  IF (ASSOCIATED(task_list%taskstop)) THEN
     DEALLOCATE(task_list%taskstop,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF
  IF (ASSOCIATED(task_list%npairs)) THEN
     DEALLOCATE(task_list%npairs,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  ENDIF

  DEALLOCATE(task_list,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
END SUBROUTINE deallocate_task_list
END MODULE task_list_types
