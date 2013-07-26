!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group   !
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

#if !defined (__DBCSR_CUDA)
  TYPE dbcsr_cuda_stream_type
    PRIVATE
    INTEGER     :: dummy !needed to make g95 happy
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
SUBROUTINE dbcsr_cuda_stream_create(this, error)
    TYPE(dbcsr_cuda_stream_type), &
      INTENT(OUT)                            :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    istat = cuda_stream_create(this%cptr)
    IF(istat /= 0 .OR. .NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_stream_create failed")
#endif
END SUBROUTINE dbcsr_cuda_stream_create


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
    istat = cuda_stream_sync(this%cptr)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_stream_synchronize failed")
#endif
END SUBROUTINE dbcsr_cuda_stream_synchronize

END MODULE dbcsr_cuda_stream
