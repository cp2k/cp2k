 !-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group   !
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


  PUBLIC :: dbcsr_cuda_init, dbcsr_cuda_device_sync,&
            dbcsr_cuda_get_n_devices


#if defined (__DBCSR_CUDA)
  INTERFACE
     FUNCTION cuda_set_device_cu (device_id) RESULT (istat) &
          BIND(C, name="dc_set_device_cu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: device_id
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_set_device_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_get_ndevices_cu (n_devices) RESULT (istat) &
          BIND(C, name="dc_get_ndevices_cu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(OUT)         :: n_devices
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_get_ndevices_cu
  END INTERFACE


  INTERFACE
     FUNCTION cuda_device_sync_cu() RESULT (istat) BIND(C, name="dc_device_sync_cu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_device_sync_cu
  END INTERFACE

#endif

CONTAINS

! *****************************************************************************
!> \brief Initializes for each thread the device. I.e. assigns the calling
!>        thread to a device. We assume all threads have devices of the same
!>        kind
!>
!> \param[in] card_num           number of the card
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_init (card_num, error)
    INTEGER, INTENT(IN), OPTIONAL            :: card_num
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: icard

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    IF (PRESENT (card_num)) THEN
       icard = card_num
    ELSE
       icard = 0
    ENDIF
    istat = cuda_set_device_cu (icard)
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error selecting GPU device.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_init


! *****************************************************************************
! \brief Returns the number of cuda devices. (wraps cudaGetDeviceCount).
! *****************************************************************************
  FUNCTION dbcsr_cuda_get_n_devices (error) RESULT (n_devices)
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    INTEGER                                  :: n_devices

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_get_n_devices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: ndev

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_get_ndevices_cu (ndev)
    n_devices = INT (ndev)
#else
    istat = -1
    n_devices = 0
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error getting device count",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END FUNCTION dbcsr_cuda_get_n_devices


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



END MODULE dbcsr_cuda_device
