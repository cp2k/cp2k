!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE acc_device
#if defined (__ACC)
  USE ISO_C_BINDING
#endif

  IMPLICIT NONE

  PUBLIC :: acc_get_ndevices, acc_set_active_device

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_device'

#if defined (__ACC)
  INTERFACE
    FUNCTION acc_get_ndevices_cu (n_devices) RESULT (istat) &
          BIND(C, name="acc_get_ndevices")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(OUT)         :: n_devices
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION acc_get_ndevices_cu

    FUNCTION acc_set_active_device_cu (dev_id) RESULT (istat) &
          BIND(C, name="acc_set_active_device")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: dev_id
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION acc_set_active_device_cu
  END INTERFACE

#endif

CONTAINS

! *****************************************************************************
!> \brief Get number of accelerator devices
!> \retval n number of accelerator devices
! *****************************************************************************
  FUNCTION acc_get_ndevices() RESULT(n)
    INTEGER                                  :: n

#if defined (__ACC)
    INTEGER                                  :: istat
#endif

     n = 0
#if defined (__ACC)
    istat = acc_get_ndevices_cu(n)
    IF (istat /= 0) STOP 'acc_get_ndevices: failed'
#endif
  END FUNCTION acc_get_ndevices


! *****************************************************************************
!> \brief Set active accelerator device
!> \param dev_id device ID
! *****************************************************************************
  SUBROUTINE acc_set_active_device(dev_id)
    INTEGER :: dev_id

#if defined (__ACC)
    INTEGER :: istat

!$omp parallel default(none) private(istat) shared(dev_id)
    istat = acc_set_active_device_cu (dev_id)
    IF (istat /= 0) STOP 'acc_set_active_device: failed'
!$omp end parallel

#else
    STOP "__ACC not compiled in"
#endif
  END SUBROUTINE acc_set_active_device

END MODULE acc_device
