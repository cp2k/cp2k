!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
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


  PUBLIC :: ma_init_conf
  PUBLIC :: ma_set_conf_print_branch, ma_get_conf_print_branch
  PUBLIC :: ma_set_conf_print_resume, ma_get_conf_print_resume
  PUBLIC :: ma_set_conf_print_full, ma_get_conf_print_full
  PUBLIC :: ma_set_conf_print_proc, ma_get_conf_print_proc
  PUBLIC :: ma_set_conf_print_thread, ma_get_conf_print_thread
  PUBLIC :: ma_set_conf_print_thread_cur, ma_get_conf_print_thread_cur
  PUBLIC :: ma_set_conf_sched, ma_get_conf_sched
  PUBLIC :: ma_set_conf_mpiSched, ma_get_conf_mpiSched  
  PUBLIC :: ma_set_conf_mempol, ma_get_conf_mempol


  ! These are for Machine Architecture internal use.
  !
  PUBLIC :: print_branch,&
            print_resume,&
            print_full
  PUBLIC :: print_proc, &
            print_thread, &
            print_thread_cur
  PUBLIC :: use_sched, &
            use_mempol, &
            mpi_sched
  PUBLIC :: isconfigured


  ! Possible prints 
  LOGICAL            :: print_branch, print_resume, &
                        print_full, print_proc, print_thread, & 
                        print_thread_cur

  CHARACTER(len=1)   :: use_sched, mpi_sched, use_mempol

  LOGICAL, SAVE :: isconfigured = .FALSE.

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

END MODULE ma_config
