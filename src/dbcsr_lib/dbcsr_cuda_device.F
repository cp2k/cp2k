 !-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA device support for DBCSR
!> \author  Urban Borstnik
!> \date    2011-09-22
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-09-22
! *****************************************************************************
MODULE dbcsr_cuda_device
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dbcsr_error_handling
  USE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_device'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.



  PUBLIC :: dbcsr_cuda_device_sync, dbcsr_cuda_set_shared_mem_config

  PUBLIC :: cudaSharedMemBankSizeDefault, cudaSharedMemBankSizeFourByte, &
            cudaSharedMemBankSizeEightByte

INTEGER, PARAMETER :: cudaSharedMemBankSizeDefault   = 0
INTEGER, PARAMETER :: cudaSharedMemBankSizeFourByte  = 1
INTEGER, PARAMETER :: cudaSharedMemBankSizeEightByte = 2

#if defined (__DBCSR_CUDA)
  INTERFACE
     FUNCTION cuda_device_sync_cu() RESULT (istat) BIND(C, name="dc_device_sync_cu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_device_sync_cu
  END INTERFACE


  INTERFACE
     FUNCTION cuda_set_shared_mem_config(config) RESULT (istat) &
       BIND(C, name="dc_set_shared_mem_config")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: config
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_set_shared_mem_config
  END INTERFACE

#endif

CONTAINS

! *****************************************************************************
! \brief Fortran wrapper for cudaDeviceSynchronize.
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_device_sync(error)
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_device_sync', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    !$OMP CRITICAL(crit_cuda_dev_sync)
    istat = cuda_device_sync_cu();
    !$OMP END CRITICAL(crit_cuda_dev_sync)
#else
    istat = -1
#endif
    IF (istat /= 0) STOP "dbcsr_cuda_device_sync: Could not synchronize all threads"
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_device_sync


! *****************************************************************************
! \brief Fortran wrapper for cudaDeviceSetSharedMemConfig.
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_set_shared_mem_config(config, error)
    INTEGER, INTENT(IN) :: config
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL  :: error
    INTEGER                                  :: istat

#if defined (__DBCSR_CUDA)
    istat = cuda_set_shared_mem_config(config)
#else
    istat = -1
#endif
    IF (istat /= 0) STOP "dbcsr_cuda_set_shared_mem_config: failed"
  END SUBROUTINE dbcsr_cuda_set_shared_mem_config

END MODULE dbcsr_cuda_device
