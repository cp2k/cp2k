!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE dbcsr_cuda_event
#if defined (__DBCSR_CUDA)
  USE ISO_C_BINDING
  USE dbcsr_cuda_stream, ONLY: dbcsr_cuda_stream_cptr
#endif

  USE dbcsr_error_handling
  USE dummy_c_bindings
  USE dbcsr_cuda_stream, ONLY: dbcsr_cuda_stream_type
  USE dbcsr_message_passing,           ONLY: mp_abort

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_event'

  PUBLIC :: dbcsr_cuda_event_type
  PUBLIC :: dbcsr_cuda_event_create, dbcsr_cuda_event_destroy
  PUBLIC :: dbcsr_cuda_event_record, dbcsr_cuda_event_query
  PUBLIC :: dbcsr_cuda_stream_wait_event, dbcsr_cuda_event_synchronize

  !PUBLIC :: dbcsr_cuda_event_elapsed_time

#if !defined (__DBCSR_CUDA)
  TYPE dbcsr_cuda_event_type
    PRIVATE
    !nothing, except for g95
    INTEGER :: dummy = 0
  END TYPE dbcsr_cuda_event_type

#else

  TYPE dbcsr_cuda_event_type
    PRIVATE
    TYPE(C_PTR) :: cptr = C_NULL_PTR
  END TYPE dbcsr_cuda_event_type

  INTERFACE
     FUNCTION cuda_event_create(event_ptr) RESULT (istat) BIND(C, name="cuda_event_create")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_create
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_destroy(event_ptr) RESULT (istat) BIND(C, name="cuda_event_destroy")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_destroy
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_query(event_ptr) RESULT (istat) BIND(C, name="cuda_event_query")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_query
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_record(event_ptr, stream_ptr) RESULT (istat) BIND(C, name="cuda_event_record")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr, stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_record
  END INTERFACE

  INTERFACE
     FUNCTION cuda_stream_wait_event(stream_ptr, event_ptr) RESULT (istat) BIND(C, name="cuda_stream_wait_event")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr, event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_wait_event
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_synchronize(event_ptr) RESULT (istat) BIND(C, name="cuda_event_synchronize")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_synchronize
  END INTERFACE
#endif

CONTAINS

! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamWaitEvent.
!> \author  Ole Schuett
!> Because of fortran cyrcular dependency restriction this can not go into dbcsr_cuda_stream.F
! *****************************************************************************
SUBROUTINE dbcsr_cuda_stream_wait_event(stream, event, error)
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_cuda_event_type), INTENT(IN)  :: event
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    stream_cptr = dbcsr_cuda_stream_cptr(stream)
    IF(.NOT. C_ASSOCIATED(event%cptr)) CALL mp_abort("dbcsr_cuda_stream_wait_event: event not allocated")
    IF(.NOT. C_ASSOCIATED(stream_cptr)) CALL mp_abort("dbcsr_cuda_stream_wait_event: stream not allocated")
    istat = cuda_stream_wait_event(stream_cptr, event%cptr)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_stream_wait_event failed")
#endif
END SUBROUTINE dbcsr_cuda_stream_wait_event


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventRecord.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_event_record(this, stream, error)
    TYPE(dbcsr_cuda_event_type), INTENT(IN)  :: this
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    stream_cptr = dbcsr_cuda_stream_cptr(stream)
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_event_record: event not allocated")
    IF(.NOT. C_ASSOCIATED(stream_cptr)) CALL mp_abort("dbcsr_cuda_event_record: stream not allocated")
    istat = cuda_event_record(this%cptr, stream_cptr)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_event_record failed")
#endif
END SUBROUTINE dbcsr_cuda_event_record


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventCreate.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_event_create(this, error)
    TYPE(dbcsr_cuda_event_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat

    IF(C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_event_create: already allocated")
    istat = cuda_event_create(this%cptr)
    IF(istat /= 0 .OR. .NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_event_create: failed")
#endif
END SUBROUTINE dbcsr_cuda_event_create


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventDestroy.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_event_destroy(this, error)
    TYPE(dbcsr_cuda_event_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_event_destroy: event not allocated")
    istat = cuda_event_destroy(this%cptr)
    IF(istat /= 0) CALL mp_abort("dbcsr_cuda_event_destroy failed")
    this%cptr = C_NULL_PTR
#endif
END SUBROUTINE dbcsr_cuda_event_destroy


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventQuery.
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION dbcsr_cuda_event_query(this, error) RESULT(res)
    TYPE(dbcsr_cuda_event_type), INTENT(IN)  :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    LOGICAL                                  :: res

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_event_query: event not allocated")
    istat = cuda_event_query(this%cptr)
    IF(istat < 0) CALL mp_abort("dbcsr_cuda_event_query failed")
    res = (istat==0)
#endif
END FUNCTION dbcsr_cuda_event_query


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventSynchronize,.
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE dbcsr_cuda_event_synchronize(this, error)
    TYPE(dbcsr_cuda_event_type), INTENT(IN)  :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) CALL mp_abort("dbcsr_cuda_event_synchronize: event not allocated")
    istat = cuda_event_synchronize(this%cptr)
    IF(istat < 0) CALL mp_abort("dbcsr_cuda_event_synchronize failed")
#endif
END SUBROUTINE dbcsr_cuda_event_synchronize

END MODULE dbcsr_cuda_event
