!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Routines that affect the DBCSR library as a whole
!> \author  Ole Schuett
!> \date    10-2014
! *****************************************************************************

MODULE dbcsr_lib
  USE dbcsr_config,                    ONLY: is_initialized
  USE dbcsr_error_handling,            ONLY: dbcsr_assert,&
                                             dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type,&
                                             dbcsr_fatal_level,&
                                             dbcsr_internal_error
  USE dbcsr_multiply_api,              ONLY: dbcsr_multiply_clear_mempools,&
                                             dbcsr_multiply_lib_finalize,&
                                             dbcsr_multiply_lib_init
  USE kinds,                           ONLY: int_1_size,&
                                             int_2_size,&
                                             int_4_size,&
                                             int_8_size

  !$ USE OMP_LIB

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_lib'

  PUBLIC :: dbcsr_init_lib, dbcsr_finalize_lib, dbcsr_clear_mempools


 CONTAINS


! *****************************************************************************
!> \brief Initialize the DBCSR library
!>
!> Prepares the DBCSR library for use.
!> \param group ...
!> \param[in,out] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_init_lib (group, error)
    INTEGER, INTENT(IN)                      :: group
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!n_stack_buffers,mem_type, n_stack_mem_regions, stack_size, nstacks,&
!INTEGER, DIMENSION(3)                    :: nxstacks
!   ---------------------------------------------------------------------------
!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

    IF (is_initialized) RETURN
    CALL dbcsr_error_set(routineN, error_handle, error)
    !
    CALL dbcsr_assert (int_1_size, "EQ", 1,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of an 8-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_2_size, "EQ", 2,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 16-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_4_size, "EQ", 4,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 32-bit integer size!",&
         __LINE__, error=error)
    CALL dbcsr_assert (int_8_size, "EQ", 8,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Incorrect assumption of a 64-bit integer size!",&
         __LINE__, error=error)

    !$omp parallel default(none)  shared(error)
    CALL dbcsr_multiply_lib_init(error)
    !$omp end parallel

    is_initialized = .TRUE.
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_init_lib


! *****************************************************************************
!> \brief Finalize the DBCSR library
!>
!> Cleans up after the DBCSR library.  Used to deallocate persistent objects.
!> \param group ...
!> \param output_unit ...
!> \param[in,out] error     error
! *****************************************************************************
  SUBROUTINE dbcsr_finalize_lib (group, output_unit, error)
    INTEGER, INTENT(IN)                      :: group, output_unit
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize_lib', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

    IF (.NOT. is_initialized) RETURN
    CALL dbcsr_error_set(routineN, error_handle, error)

    IF(output_unit>0)THEN
       WRITE (UNIT=output_unit,FMT="(/,T2,A)") REPEAT("-",79)
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=output_unit,FMT="(T2,A,T35,A,T80,A)") "-","DBCSR STATISTICS","-"
       WRITE (UNIT=output_unit,FMT="(T2,A,T80,A)") "-","-"
       WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)
    END IF

    !$omp parallel default(none) shared(output_unit, group, error)
    CALL dbcsr_multiply_lib_finalize(group, output_unit, error)
    !$omp end parallel
    IF(output_unit>0) WRITE (UNIT=output_unit,FMT="(T2,A)") REPEAT("-",79)

    is_initialized = .FALSE.
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_finalize_lib


! *****************************************************************************
!> \brief  Deallocate memory contained in mempools
!> \param error ...
! *****************************************************************************
  SUBROUTINE dbcsr_clear_mempools(error)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    !$omp parallel default(none) shared( error)
    CALL dbcsr_multiply_clear_mempools(error)
    !$omp end parallel
  END SUBROUTINE dbcsr_clear_mempools

END MODULE dbcsr_lib
