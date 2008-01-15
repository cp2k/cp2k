!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****h* cp2k/task_list_types *
!!
!!   NAME
!!     task_list_typess
!!
!!   FUNCTION
!!     types for task lists
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Joost VandeVondele 
!!
!!   MODIFICATION HISTORY
!!     01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!!
!!   SOURCE
!****************************************************************************
MODULE task_list_types
  USE kinds,                           ONLY: dp,&
                                             int_8
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_types'

  TYPE task_list_type
    INTEGER(kind=int_8), DIMENSION(:, :), POINTER :: tasks
    REAL(KIND=dp), DIMENSION(:, :), POINTER       :: dist_ab
    INTEGER, DIMENSION(:), POINTER                :: atom_pair_send, atom_pair_recv
    INTEGER                                       :: ntasks
  END TYPE task_list_type

  PUBLIC :: task_list_type

  PUBLIC :: allocate_task_list,&
            deallocate_task_list

CONTAINS
!!****f* task_lists/allocate_task_list *
!!
!!   NAME
!!     allocate_task_list
!!
!!   FUNCTION
!!     allocates and initialised the components of the task_list_type
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE allocate_task_list(task_list,error)
  TYPE(task_list_type), POINTER :: task_list
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
  task_list%ntasks=0
END SUBROUTINE allocate_task_list

!!****f* task_lists/deallocate_task_list*
!!
!!   NAME
!!     deallocate_task_list_type
!!
!!   FUNCTION
!!     deallocates the components and the object itself
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE deallocate_task_list(task_list,error)
  TYPE(task_list_type), POINTER :: task_list
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
  DEALLOCATE(task_list,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
END SUBROUTINE deallocate_task_list
END MODULE task_list_types
