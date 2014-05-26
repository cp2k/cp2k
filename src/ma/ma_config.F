!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Configuration options for Machine Architecture
!> \author  Christiane Pousa Ribeiro
!> \date    2011-17-08
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-17-08
! *****************************************************************************
MODULE ma_config

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_config'


  PUBLIC :: ma_set_conf_print_branch, ma_get_conf_print_branch
  PUBLIC :: ma_set_conf_print_resume, ma_get_conf_print_resume
  PUBLIC :: ma_set_conf_print_full, ma_get_conf_print_full
  PUBLIC :: ma_set_conf_print_proc, ma_get_conf_print_proc
  PUBLIC :: ma_set_conf_print_thread, ma_get_conf_print_thread
  PUBLIC :: ma_set_conf_print_thread_cur, ma_get_conf_print_thread_cur
  PUBLIC :: ma_set_conf_sched, ma_get_conf_sched
  PUBLIC :: ma_set_conf_mpiSched, ma_get_conf_mpiSched
  PUBLIC :: ma_set_conf_mempol, ma_get_conf_mempol
  PUBLIC :: ma_set_conf_comm_thread, ma_get_conf_comm_thread
  PUBLIC :: ma_set_all_affinty, ma_valid_inputs
  PUBLIC :: ma_set_conf_mpi_reordering, ma_get_conf_mpi_reordering

  ! These are for Machine Architecture internal use.
  !
  PUBLIC :: print_branch,&
            print_resume,&
            print_full
  PUBLIC :: print_proc,&
            print_thread,&
            print_thread_cur
  PUBLIC :: use_sched,&
            use_mempol,&
            mpi_sched,&
            comm_thread,&
            mpi_reorder
  PUBLIC :: isconfigured, isdefault, hasnet

! this is for manual setting of all affinity
  PUBLIC process, memory, thread

  INTEGER, DIMENSION(:), POINTER :: process
  INTEGER, DIMENSION(:), POINTER :: thread
  INTEGER, DIMENSION(:), POINTER :: memory

  ! Possible prints
  LOGICAL            :: print_branch, print_resume, &
                        print_full, print_proc, print_thread, &
                        print_thread_cur, comm_thread

  CHARACTER(len=1)   :: use_sched, mpi_sched, use_mempol, mpi_reorder

  LOGICAL, SAVE :: isconfigured = .FALSE.
  LOGICAL, SAVE :: isdefault = .FALSE.
  LOGICAL, SAVE :: hasnet = .TRUE.

CONTAINS


  SUBROUTINE ma_init_conf ()
    CHARACTER(len=*), PARAMETER :: routineN = 'ma_init_conf', &
      routineP = moduleN//':'//routineN

    print_branch = .FALSE.
    print_resume = .FALSE.
    print_full = .FALSE.
    print_proc = .FALSE.
    print_thread = .FALSE.
    print_thread_cur = .FALSE.
    use_sched = 'N'
    mpi_sched = 'N'
    mpi_reorder = 'N' 
    comm_thread = .FALSE.
  END SUBROUTINE ma_init_conf

  SUBROUTINE ma_set_conf_print_branch (branch)
    LOGICAL, INTENT(IN)                      :: branch

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_print_branch', &
      routineP = moduleN//':'//routineN

    print_branch = branch
  END SUBROUTINE ma_set_conf_print_branch

  FUNCTION ma_get_conf_print_branch () RESULT (branch)
    LOGICAL                                  :: branch

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_print_branch', &
      routineP = moduleN//':'//routineN

    branch = print_branch
  END FUNCTION ma_get_conf_print_branch

  SUBROUTINE ma_set_conf_print_resume (resume)
    LOGICAL, INTENT(IN)                      :: resume

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_print_resume', &
      routineP = moduleN//':'//routineN

    print_resume = resume
  END SUBROUTINE ma_set_conf_print_resume

  FUNCTION ma_get_conf_print_resume () RESULT (resume)
    LOGICAL                                  :: resume

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_print_resume', &
      routineP = moduleN//':'//routineN

    resume = print_resume
  END FUNCTION ma_get_conf_print_resume

  SUBROUTINE ma_set_conf_print_full (full)
    LOGICAL, INTENT(IN)                      :: full

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_print_full', &
      routineP = moduleN//':'//routineN

    print_full = full
  END SUBROUTINE ma_set_conf_print_full

  FUNCTION ma_get_conf_print_full () RESULT (full)
    LOGICAL                                  :: full

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_print_full', &
      routineP = moduleN//':'//routineN

    full = print_full
  END FUNCTION ma_get_conf_print_full

  SUBROUTINE ma_set_conf_print_proc (proc)
    LOGICAL, INTENT(IN)                      :: proc

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_print_proc', &
      routineP = moduleN//':'//routineN

    print_proc = proc
  END SUBROUTINE ma_set_conf_print_proc

  SUBROUTINE ma_set_conf_comm_thread (comm)
    LOGICAL, INTENT(IN)                      :: comm

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_comm_thread', &
      routineP = moduleN//':'//routineN

    comm_thread = comm

  END SUBROUTINE ma_set_conf_comm_thread

  FUNCTION ma_get_conf_comm_thread () RESULT (comm)
    LOGICAL                                  :: comm

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_comm_thread', &
      routineP = moduleN//':'//routineN

    comm = comm_thread
  END FUNCTION ma_get_conf_comm_thread

  FUNCTION ma_get_conf_print_proc () RESULT (proc)
    LOGICAL                                  :: proc

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_print_proc', &
      routineP = moduleN//':'//routineN

    proc = print_proc
  END FUNCTION ma_get_conf_print_proc

  SUBROUTINE ma_set_conf_print_thread (thread)
    LOGICAL, INTENT(IN)                      :: thread

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_print_thread', &
      routineP = moduleN//':'//routineN

    print_thread = thread
  END SUBROUTINE ma_set_conf_print_thread

  FUNCTION ma_get_conf_print_thread () RESULT (thread)
    LOGICAL                                  :: thread

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_print_thread', &
      routineP = moduleN//':'//routineN

    thread = print_thread
  END FUNCTION ma_get_conf_print_thread

  SUBROUTINE ma_set_conf_print_thread_cur (thread_cur)
    LOGICAL, INTENT(IN)                      :: thread_cur

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_print_thread_cur', &
      routineP = moduleN//':'//routineN

    print_thread_cur = thread_cur
  END SUBROUTINE ma_set_conf_print_thread_cur

  FUNCTION ma_get_conf_print_thread_cur () RESULT (thread_cur)
    LOGICAL                                  :: thread_cur

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_print_thread_cur', &
      routineP = moduleN//':'//routineN

    thread_cur = print_thread_cur
  END FUNCTION ma_get_conf_print_thread_cur

  SUBROUTINE ma_set_conf_sched (sched)
    CHARACTER(len=1), INTENT(IN)             :: sched

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_sched', &
      routineP = moduleN//':'//routineN

    use_sched = sched
  END SUBROUTINE ma_set_conf_sched

  SUBROUTINE ma_set_conf_mpiSched (sched)
    CHARACTER(len=1), INTENT(IN)             :: sched

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_mpiSched', &
      routineP = moduleN//':'//routineN

    mpi_sched = sched
  END SUBROUTINE ma_set_conf_mpiSched

  FUNCTION ma_get_conf_mpiSched () RESULT (sched)
    CHARACTER(len=1)                         :: sched

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_mpiSched', &
      routineP = moduleN//':'//routineN

    sched = mpi_sched

  END FUNCTION ma_get_conf_mpiSched


  FUNCTION ma_get_conf_sched () RESULT (sched)
    CHARACTER(len=1)                         :: sched

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_sched', &
      routineP = moduleN//':'//routineN

    sched = use_sched
  END FUNCTION ma_get_conf_sched

  SUBROUTINE ma_set_conf_mempol (mempol)
    CHARACTER(len=1), INTENT(IN)             :: mempol

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_mempol', &
      routineP = moduleN//':'//routineN

    use_mempol = mempol
  END SUBROUTINE ma_set_conf_mempol

  FUNCTION ma_get_conf_mempol () RESULT (mempol)
    CHARACTER(len=1)                         :: mempol

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_mempol', &
      routineP = moduleN//':'//routineN

    mempol = use_mempol
  END FUNCTION ma_get_conf_mempol

  SUBROUTINE ma_set_all_affinty (proc, thr, mem)
    INTEGER, DIMENSION(:), INTENT(IN)        :: proc, thr, mem

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_all_affinty', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: size_mem, size_proc, size_thr

    size_proc = SIZE(proc)
    size_mem = SIZE(mem)
    size_thr = SIZE(thr)
    ALLOCATE(process(size_proc))
    ALLOCATE(memory(size_mem))
    ALLOCATE(thread(size_thr))
    process(:) = proc(:)
    memory(:)= mem(:)
    thread(:)= thr(:)

  END SUBROUTINE ma_set_all_affinty

  LOGICAL FUNCTION ma_valid_inputs(nMPI) RESULT(istat)
    INTEGER :: nMPI

    istat = .FALSE.
    IF (SIZE(process) .GE. nMPI .AND. &
        SIZE(memory) .GE. nMPI) THEN
      istat = .TRUE.
    END IF
  END FUNCTION ma_valid_inputs

  SUBROUTINE ma_set_conf_mpi_reordering(mpi_reordering)
    CHARACTER(len=1), INTENT(IN)             :: mpi_reordering

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_set_conf_mpi_reordering', &
      routineP = moduleN//':'//routineN

    mpi_reorder = mpi_reordering
  END SUBROUTINE ma_set_conf_mpi_reordering

  FUNCTION ma_get_conf_mpi_reordering () RESULT (mpi_reordering)
    CHARACTER(len=1)                         :: mpi_reordering

    CHARACTER(len=*), PARAMETER :: routineN = 'ma_get_conf_mpi_reordering', &
      routineP = moduleN//':'//routineN

    mpi_reordering = mpi_reorder
  END FUNCTION ma_get_conf_mpi_reordering

END MODULE ma_config
