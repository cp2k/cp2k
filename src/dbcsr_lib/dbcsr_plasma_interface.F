!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Interface to PLASMA
!> \author  V. Weber
!> \date    2010
!> \version 0.1
!>
!> <b>Modification history:</b>
!  - Created 2010
! *****************************************************************************
MODULE dbcsr_plasma_interface

  USE dbcsr_error_handling

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_plasma_interface'

  PUBLIC :: dbcsr_plasma_init, dbcsr_plasma_finalize

CONTAINS

! *****************************************************************************
!> \brief Initialize plasma
!>
! *****************************************************************************
  SUBROUTINE dbcsr_plasma_init(plasma_is_set, error)
    LOGICAL, INTENT(out)                     :: plasma_is_set
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_plasma_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, info, nthreads

    CALL dbcsr_error_set(routineN, error_handler, error)

    plasma_is_set = .FALSE.
    nthreads = 1
!$  nthreads = OMP_GET_MAX_THREADS()

#if defined(__PLASMA)
    CALL plasma_init( nthreads, info )
    CALL dbcsr_assert( info==0, dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "plasma_init failed", __LINE__, error=error)
    plasma_is_set = .TRUE.
#endif

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_plasma_init

! *****************************************************************************
!> \brief Finalize plasma
!>
! *****************************************************************************
  SUBROUTINE dbcsr_plasma_finalize(error)
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_plasma_finalize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handler, info

    CALL dbcsr_error_set(routineN, error_handler, error)

#if defined(__PLASMA)
    CALL plasma_finalize( info )
    CALL dbcsr_assert( info==0, dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "plasma_finalize failed", __LINE__, error=error)
#endif

    CALL dbcsr_error_stop(error_handler, error)
  END SUBROUTINE dbcsr_plasma_finalize

END MODULE dbcsr_plasma_interface
