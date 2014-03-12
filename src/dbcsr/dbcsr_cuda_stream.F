!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE dbcsr_cuda_stream
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dbcsr_error_handling
  USE dummy_c_bindings
  USE dbcsr_message_passing,           ONLY: mp_abort

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_stream'

  PUBLIC :: dbcsr_cuda_stream_type
  PUBLIC :: dbcsr_cuda_stream_create, dbcsr_cuda_stream_destroy
  PUBLIC :: dbcsr_cuda_stream_synchronize
  PUBLIC :: dbcsr_cuda_stream_priority_range
  PUBLIC :: dbcsr_cuda_stream_set_name
  PUBLIC :: dbcsr_cuda_stream_equal, dbcsr_cuda_stream_associated

#if !defined (__DBCSR_CUDA)
  TYPE dbcsr_cuda_stream_type
!MK PRIVATE                  ! workaround for g95
    INTEGER     :: dummy = 0 ! needed to make g95 happy
  END TYPE dbcsr_cuda_stream_type

#else
  PUBLIC :: dbcsr_cuda_stream_cptr

  TYPE dbcsr_cuda_stream_type
    PRIVATE
    TYPE(C_PTR) :: cptr = C_NULL_PTR
  END TYPE dbcsr_cuda_stream_type

  TYPE(dbcsr_cuda_stream_type), SAVE :: dbcsr_cuda_stream_default

  INTERFACE
     FUNCTION cuda_stream_create(stream_ptr) RESULT (istat) BIND(C, name="cuda_stream_create")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_create
  END INTERFACE

#if ! defined (__HAS_NO_CUDA_STREAM_PRIORITIES)
  INTERFACE
     FUNCTION cuda_stream_create_with_priority(stream_ptr, priority) RESULT (istat) BIND(C, name="cuda_stream_create_with_priority")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: stream_ptr
    INTEGER(KIND=C_INT), VALUE               :: priority
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_create_with_priority
  END INTERFACE

  INTERFACE
     FUNCTION cuda_stream_priority_range(least, greatest) RESULT (istat) BIND(C, name="cuda_stream_priority_range")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: least, greatest, istat

    END FUNCTION cuda_stream_priority_range
  END INTERFACE
#endif

  INTERFACE
     FUNCTION cuda_stream_destroy(stream_ptr) RESULT (istat) BIND(C, name="cuda_stream_destroy")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_destroy
  END INTERFACE

  INTERFACE
    FUNCTION cuda_stream_sync(stream_ptr) RESULT (istat) BIND(C, name="cuda_stream_sync")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_sync
  END INTERFACE

#if defined( __CUDA_PROFILING)
  INTERFACE
     SUBROUTINE cuda_stream_set_name(stream_ptr, name) BIND(C, name="cuda_stream_set_name")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr
    CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: name

    END SUBROUTINE cuda_stream_set_name
  END INTERFACE
#endif

#endif
CONTAINS




! *****************************************************************************
!> \brief Returns C-pointer of given stream.
!> \author  Ole Schuett
! *****************************************************************************
#if defined (__DBCSR_CUDA)
FUNCTION dbcsr_cuda_stream_cptr(this) RESULT(res)
    TYPE(dbcsr_cuda_stream_type), INTENT(in) :: this
    TYPE(C_PTR)                              :: res

    res = this%cptr
END FUNCTION dbcsr_cuda_stream_cptr
#endif


! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamCreate.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_stream_create(this, priority, error)
    TYPE(dbcsr_cuda_stream_type), &
      INTENT(OUT)                            :: this
    INTEGER, INTENT(IN), OPTIONAL            :: priority
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat

    IF(C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_stream_create: stream already allocated")

#if defined (__HAS_NO_CUDA_STREAM_PRIORITIES)
    istat = cuda_stream_create(this%cptr)
#else
    IF(PRESENT(priority)) THEN
       istat = cuda_stream_create_with_priority(this%cptr, priority)
    ELSE
       istat = cuda_stream_create(this%cptr)
    ENDIF
#endif

    IF(istat /= 0 .OR. .NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_stream_create failed")
#endif
END SUBROUTINE dbcsr_cuda_stream_create

! *****************************************************************************
!> \brief Fortran-wrapper for nvtxNameCudaStreamA.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_stream_set_name(this, name, error)
    TYPE(dbcsr_cuda_stream_type), &
      INTENT(IN)                            :: this
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if defined( __CUDA_PROFILING)
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_stream_set_name: stream not allocated")
    CALL cuda_stream_set_name(this%cptr, name//c_null_char)
#endif
END SUBROUTINE dbcsr_cuda_stream_set_name

! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamDestroy.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_stream_destroy(this, error)
    TYPE(dbcsr_cuda_stream_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_stream_destroy: stream not allocated")
    istat = cuda_stream_destroy(this%cptr)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_stream_destroy failed")
    this%cptr = C_NULL_PTR
#endif
END SUBROUTINE dbcsr_cuda_stream_destroy


! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamSynchronize.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_stream_synchronize(this, error)
    TYPE(dbcsr_cuda_stream_type), &
      INTENT(IN)                             :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_stream_synchronize: stream not allocated")
    istat = cuda_stream_sync(this%cptr)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_stream_synchronize failed")
#endif
END SUBROUTINE dbcsr_cuda_stream_synchronize


! *****************************************************************************
!> \brief Fortran-wrapper for cudaDeviceGetStreamPriorityRange.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_stream_priority_range(least, greatest, error)
    INTEGER, INTENT(OUT)                     :: least, greatest
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")

#elif defined (__HAS_NO_CUDA_STREAM_PRIORITIES)
    least = -1
    greatest = -1

#else
    INTEGER                                  :: istat
    istat = cuda_stream_priority_range(least, greatest)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_stream_priority_range failed")
#endif
END SUBROUTINE dbcsr_cuda_stream_priority_range


! *****************************************************************************
!> \brief Checks if two streams are equal
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION dbcsr_cuda_stream_equal(this, other) RESULT(res)
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: this, other
    LOGICAL                                  :: res
#if ! defined (__DBCSR_CUDA)
    res = .TRUE.
#else
    res = C_ASSOCIATED(this%cptr, other%cptr)
#endif
END FUNCTION dbcsr_cuda_stream_equal


! *****************************************************************************
!> \brief Checks if a streams is associated
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION dbcsr_cuda_stream_associated(this) RESULT(res)
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: this
    LOGICAL                                  :: res
#if ! defined (__DBCSR_CUDA)
    res = .FALSE.
#else
    res = C_ASSOCIATED(this%cptr)
#endif
END FUNCTION dbcsr_cuda_stream_associated

END MODULE dbcsr_cuda_stream
