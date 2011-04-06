!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Configuration options for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-05
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-05
! *****************************************************************************
MODULE dbcsr_config
  USE dbcsr_error_handling

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_config'


  PUBLIC :: dbcsr_set_conf_subcomm, dbcsr_get_conf_subcomm
  PUBLIC :: dbcsr_set_conf_combtypes, dbcsr_get_conf_combtypes
  PUBLIC :: dbcsr_set_conf_mm_driver, dbcsr_get_conf_mm_driver

  ! These are for DBCSR internal use.
  !
  PUBLIC :: use_subcommunicators
  PUBLIC :: use_combined_types
  PUBLIC :: mm_driver,&
       mm_driver_blas,&
       mm_driver_matmul,&
       mm_driver_smm,&
       mm_driver_plasma,&
       has_blas,&
       has_smm_gemm,&
       has_plasma


  ! Allocates subcommunicators for process rows and columns.
  LOGICAL :: use_subcommunicators = .FALSE.

  ! Use combined data types for MPI transfers.
  LOGICAL :: use_combined_types = .FALSE.

  ! Drivers to use for matrix multiplications
  INTEGER, PARAMETER :: mm_driver_blas       = 1
  INTEGER, PARAMETER :: mm_driver_matmul     = 2
  INTEGER, PARAMETER :: mm_driver_smm        = 3
  INTEGER, PARAMETER :: mm_driver_plasma     = 4

  INTEGER :: mm_driver = mm_driver_smm

  LOGICAL, PARAMETER :: has_blas = .TRUE.

#if defined (__HAS_smm_dnn)
  LOGICAL, PARAMETER :: has_smm_gemm = .TRUE.
#else
  LOGICAL, PARAMETER :: has_smm_gemm = .FALSE.
#endif

#if defined(__PLASMA)
  LOGICAL, PARAMETER :: has_plasma = .TRUE.
#else
  LOGICAL, PARAMETER :: has_plasma = .FALSE.
#endif





CONTAINS


  SUBROUTINE dbcsr_set_conf_mm_driver (driver, error)
    INTEGER, INTENT(IN)                      :: driver
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_mm_driver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set(routineN, error_handle, error)
    SELECT CASE (driver)
    CASE (mm_driver_matmul)
       mm_driver = mm_driver_matmul
    CASE (mm_driver_blas)
       CALL dbcsr_assert (has_blas,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Support for BLAS not compiled in.", __LINE__, error=error)
       mm_driver = mm_driver_blas
    CASE (mm_driver_smm)
       mm_driver = mm_driver_smm
    CASE (mm_driver_plasma)
       CALL dbcsr_assert (has_plasma,&
            dbcsr_fatal_level, dbcsr_caller_error, routineN,&
            "Support for PLASMA not compiled in.", __LINE__, error=error)
       mm_driver = mm_driver_plasma
    END SELECT
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_mm_driver

  FUNCTION dbcsr_get_conf_mm_driver (error) RESULT (driver)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: driver

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_mm_driver', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    driver = mm_driver
    CALL dbcsr_error_stop(error_handle, error)
  END FUNCTION dbcsr_get_conf_mm_driver


  SUBROUTINE dbcsr_set_conf_subcomm (use_subcomms, error)
    LOGICAL, INTENT(IN)                      :: use_subcomms
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_subcomm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_subcommunicators = use_subcomms
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_subcomm

  FUNCTION dbcsr_get_conf_subcomm (error) RESULT (use_subcomms)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: use_subcomms

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_subcomm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_subcomms = use_subcommunicators
    CALL dbcsr_error_stop(error_handle, error)
  END FUNCTION dbcsr_get_conf_subcomm

  SUBROUTINE dbcsr_set_conf_combtypes (use_combtypes, error)
    LOGICAL, INTENT(IN)                      :: use_combtypes
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_set_conf_combtypes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_combined_types = use_combtypes
    CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dbcsr_set_conf_combtypes

  FUNCTION dbcsr_get_conf_combtypes (error) RESULT (use_combtypes)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: use_combtypes

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_get_conf_combtypes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

    CALL dbcsr_error_set(routineN, error_handle, error)
    use_combtypes = use_combined_types
    CALL dbcsr_error_stop(error_handle, error)
  END FUNCTION dbcsr_get_conf_combtypes


END MODULE dbcsr_config
