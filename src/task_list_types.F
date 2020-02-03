!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types for task lists
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE task_list_types

   USE kinds,                           ONLY: dp,&
                                              int_8
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_types'

   TYPE task_type
      INTEGER(kind=int_8)                           :: destination = 0
      INTEGER(kind=int_8)                           :: source = 0
      INTEGER(kind=int_8)                           :: pair = 0
      INTEGER(kind=int_8)                           :: dist_type = 0
      INTEGER(kind=int_8)                           :: cost = 0
      INTEGER(kind=int_8)                           :: subpatch_pattern = 0
   END TYPE task_type

   TYPE task_list_type
      TYPE(task_type), DIMENSION(:), POINTER        :: tasks => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER       :: dist_ab
      INTEGER(kind=int_8), DIMENSION(:), POINTER    :: atom_pair_send, atom_pair_recv
      INTEGER                                       :: ntasks
      INTEGER, DIMENSION(:, :), POINTER             :: taskstart, taskstop
      INTEGER, DIMENSION(:), POINTER                :: npairs
   END TYPE task_list_type

   INTEGER, PARAMETER                               :: task_size_in_int8 = 6

   PUBLIC :: task_type, task_list_type
   PUBLIC :: task_size_in_int8, serialize_task, deserialize_task
   PUBLIC :: allocate_task_list, deallocate_task_list, reallocate_tasks

CONTAINS

! **************************************************************************************************
!> \brief allocates and initialised the components of the task_list_type
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE allocate_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_task_list', &
         routineP = moduleN//':'//routineN

      ALLOCATE (task_list)

      NULLIFY (task_list%tasks)
      NULLIFY (task_list%dist_ab)
      NULLIFY (task_list%atom_pair_send)
      NULLIFY (task_list%atom_pair_recv)
      NULLIFY (task_list%taskstart)
      NULLIFY (task_list%taskstop)
      NULLIFY (task_list%npairs)
      task_list%ntasks = 0
   END SUBROUTINE allocate_task_list

! **************************************************************************************************
!> \brief Grow an array of tasks while preserving the existing entries.
!> \param tasks ...
!> \param new_size ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE reallocate_tasks(tasks, new_size)
      TYPE(task_type), DIMENSION(:), INTENT(inout), &
         POINTER                                         :: tasks
      INTEGER, INTENT(in)                                :: new_size

      TYPE(task_type), DIMENSION(:), POINTER             :: old_tasks

      IF (.NOT. ASSOCIATED(tasks)) THEN
         ALLOCATE (tasks(new_size))
         RETURN
      ENDIF

      IF (new_size <= SIZE(tasks)) RETURN

      old_tasks => tasks
      NULLIFY (tasks)

      ALLOCATE (tasks(new_size))
      tasks(1:SIZE(old_tasks)) = old_tasks(:)
      DEALLOCATE (old_tasks)
   END SUBROUTINE reallocate_tasks

! **************************************************************************************************
!> \brief deallocates the components and the object itself
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE deallocate_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_task_list', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(task_list%tasks)) THEN
         DEALLOCATE (task_list%tasks)
      ENDIF
      IF (ASSOCIATED(task_list%dist_ab)) THEN
         DEALLOCATE (task_list%dist_ab)
      ENDIF
      IF (ASSOCIATED(task_list%atom_pair_send)) THEN
         DEALLOCATE (task_list%atom_pair_send)
      ENDIF
      IF (ASSOCIATED(task_list%atom_pair_recv)) THEN
         DEALLOCATE (task_list%atom_pair_recv)
      ENDIF
      IF (ASSOCIATED(task_list%taskstart)) THEN
         DEALLOCATE (task_list%taskstart)
      ENDIF
      IF (ASSOCIATED(task_list%taskstop)) THEN
         DEALLOCATE (task_list%taskstop)
      ENDIF
      IF (ASSOCIATED(task_list%npairs)) THEN
         DEALLOCATE (task_list%npairs)
      ENDIF

      DEALLOCATE (task_list)
   END SUBROUTINE deallocate_task_list

! **************************************************************************************************
!> \brief Serialize a task into an integer array. Used for MPI communication.
!> \param task ...
!> \param serialized_task ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE serialize_task(task, serialized_task)
      TYPE(task_type), INTENT(IN)                        :: task
      INTEGER(KIND=int_8), &
         DIMENSION(task_size_in_int8), INTENT(OUT)       :: serialized_task

      serialized_task(1) = task%destination
      serialized_task(2) = task%source
      serialized_task(3) = task%pair
      serialized_task(4) = task%dist_type
      serialized_task(5) = task%cost
      serialized_task(6) = task%subpatch_pattern
   END SUBROUTINE serialize_task

! **************************************************************************************************
!> \brief De-serialize a task from an integer array. Used for MPI communication.
!> \param task ...
!> \param serialized_task ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE deserialize_task(task, serialized_task)
      TYPE(task_type), INTENT(OUT)                       :: task
      INTEGER(KIND=int_8), &
         DIMENSION(task_size_in_int8), INTENT(IN)        :: serialized_task

      task%destination = serialized_task(1)
      task%source = serialized_task(2)
      task%pair = serialized_task(3)
      task%dist_type = serialized_task(4)
      task%cost = serialized_task(5)
      task%subpatch_pattern = serialized_task(6)
   END SUBROUTINE deserialize_task

END MODULE task_list_types
