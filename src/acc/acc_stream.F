!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE acc_stream
  USE ISO_C_BINDING

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_stream'

  PUBLIC :: acc_stream_type
  PUBLIC :: acc_stream_create, acc_stream_destroy
  PUBLIC :: acc_stream_synchronize
  PUBLIC :: acc_stream_priority_range
  PUBLIC :: acc_stream_equal, acc_stream_associated
  PUBLIC :: acc_stream_cptr

  TYPE acc_stream_type
    PRIVATE
    TYPE(C_PTR) :: cptr = C_NULL_PTR
  END TYPE acc_stream_type

#if defined (__ACC)

  INTERFACE
     FUNCTION cuda_stream_create(stream_ptr, name, priority) RESULT (istat) BIND(C, name="acc_stream_create")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: stream_ptr
    CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: name
    INTEGER(KIND=C_INT), VALUE               :: priority
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_create
  END INTERFACE

  INTERFACE
     FUNCTION cuda_stream_priority_range(least, greatest) RESULT (istat) BIND(C, name="acc_stream_priority_range")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: least, greatest, istat

    END FUNCTION cuda_stream_priority_range
  END INTERFACE

  INTERFACE
     FUNCTION cuda_stream_destroy(stream_ptr) RESULT (istat) BIND(C, name="acc_stream_destroy")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_destroy
  END INTERFACE

  INTERFACE
    FUNCTION cuda_stream_sync(stream_ptr) RESULT (istat) BIND(C, name="acc_stream_sync")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_sync
  END INTERFACE

#endif
CONTAINS




! *****************************************************************************
!> \brief Returns C-pointer of given stream.
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_stream_cptr(this) RESULT(res)
    TYPE(acc_stream_type), INTENT(in)        :: this
    TYPE(C_PTR)                              :: res

    res = this%cptr
END FUNCTION acc_stream_cptr


! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamCreate.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_create(this, name, priority)
    TYPE(acc_stream_type),INTENT(OUT) :: this
    CHARACTER(LEN=*), INTENT(IN)             :: name
    INTEGER, INTENT(IN), OPTIONAL            :: priority

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat, my_priority

    my_priority = -1
    IF(PRESENT(priority)) &
       my_priority = priority

    IF(C_ASSOCIATED(this%cptr)) STOP "acc_stream_create: stream already allocated"

    istat = cuda_stream_create(this%cptr, name//c_null_char, my_priority)

    IF(istat /= 0 .OR. .NOT. C_ASSOCIATED(this%cptr)) STOP "acc_stream_create failed"
#endif
END SUBROUTINE acc_stream_create

! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamDestroy.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_destroy(this)
    TYPE(acc_stream_type), &
      INTENT(INOUT)                          :: this

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) STOP "acc_stream_destroy: stream not allocated"
    istat = cuda_stream_destroy(this%cptr)
    IF(istat /= 0) STOP "acc_stream_destroy failed"
    this%cptr = C_NULL_PTR
#endif
END SUBROUTINE acc_stream_destroy


! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamSynchronize.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_synchronize(this)
    TYPE(acc_stream_type), &
      INTENT(IN)                             :: this

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) STOP "acc_stream_synchronize: stream not allocated"
    istat = cuda_stream_sync(this%cptr)
    IF(istat /= 0) STOP "acc_stream_synchronize failed"
#endif
END SUBROUTINE acc_stream_synchronize


! *****************************************************************************
!> \brief Fortran-wrapper for cudaDeviceGetStreamPriorityRange.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_priority_range(least, greatest)
    INTEGER, INTENT(OUT)                     :: least, greatest

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    istat = cuda_stream_priority_range(least, greatest)
    IF(istat /= 0) STOP "acc_stream_priority_range failed"
#endif
END SUBROUTINE acc_stream_priority_range


! *****************************************************************************
!> \brief Checks if two streams are equal
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_stream_equal(this, other) RESULT(res)
    TYPE(acc_stream_type), INTENT(IN) :: this, other
    LOGICAL                                  :: res
#if ! defined (__ACC)
    res = .TRUE.
#else
    res = C_ASSOCIATED(this%cptr, other%cptr)
#endif
END FUNCTION acc_stream_equal


! *****************************************************************************
!> \brief Checks if a streams is associated
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_stream_associated(this) RESULT(res)
    TYPE(acc_stream_type), INTENT(IN) :: this
    LOGICAL                                  :: res
#if ! defined (__ACC)
    res = .FALSE.
#else
    res = C_ASSOCIATED(this%cptr)
#endif
END FUNCTION acc_stream_associated

END MODULE acc_stream
