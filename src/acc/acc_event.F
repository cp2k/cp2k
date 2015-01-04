!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE acc_event
#if defined (__ACC)
  USE ISO_C_BINDING
#endif
  USE acc_stream,                      ONLY: acc_stream_cptr,&
                                             acc_stream_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_event'

  PUBLIC :: acc_event_type
  PUBLIC :: acc_event_create, acc_event_destroy
  PUBLIC :: acc_event_record, acc_event_query
  PUBLIC :: acc_stream_wait_event, acc_event_synchronize


  TYPE acc_event_type
    PRIVATE
#if defined (__ACC)
    TYPE(C_PTR) :: cptr = C_NULL_PTR
#else
    INTEGER :: dummy = 1
#endif
  END TYPE acc_event_type

#if defined (__ACC)

  INTERFACE
     FUNCTION cuda_event_create(event_ptr) RESULT (istat) BIND(C, name="acc_event_create")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_create
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_destroy(event_ptr) RESULT (istat) BIND(C, name="acc_event_destroy")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_destroy
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_query(event_ptr, has_occured) RESULT (istat) BIND(C, name="acc_event_query")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr
    INTEGER(KIND=C_INT)                      :: has_occured, istat

    END FUNCTION cuda_event_query
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_record(event_ptr, stream_ptr) RESULT (istat) BIND(C, name="acc_event_record")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr, stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_record
  END INTERFACE

  INTERFACE
     FUNCTION cuda_stream_wait_event(stream_ptr, event_ptr) RESULT (istat) BIND(C, name="acc_stream_wait_event")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: stream_ptr, event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_stream_wait_event
  END INTERFACE

  INTERFACE
     FUNCTION cuda_event_synchronize(event_ptr) RESULT (istat) BIND(C, name="acc_event_synchronize")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: event_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION cuda_event_synchronize
  END INTERFACE
#endif

CONTAINS

! *****************************************************************************
!> \brief Fortran-wrapper for cudaStreamWaitEvent.
!>        Because of fortran circular dependency restriction this can not go into acc_stream.F
!> \param[in] stream stream
!> \param[in] event event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_stream_wait_event(stream, event)
    TYPE(acc_stream_type), INTENT(IN) :: stream
    TYPE(acc_event_type), INTENT(IN)  :: event

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    stream_cptr = acc_stream_cptr(stream)
    IF(.NOT. C_ASSOCIATED(event%cptr)) STOP "acc_stream_wait_event: event not allocated"
    IF(.NOT. C_ASSOCIATED(stream_cptr)) STOP "acc_stream_wait_event: stream not allocated"
    istat = cuda_stream_wait_event(stream_cptr, event%cptr)
    IF(istat /= 0) STOP "acc_stream_wait_event failed"
#endif
END SUBROUTINE acc_stream_wait_event


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventRecord.
!> \param[in] this event
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_record(this, stream)
    TYPE(acc_event_type), INTENT(IN)  :: this
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    stream_cptr = acc_stream_cptr(stream)
    IF(.NOT. C_ASSOCIATED(this%cptr)) STOP "acc_event_record: event not allocated"
    IF(.NOT. C_ASSOCIATED(stream_cptr)) STOP "acc_event_record: stream not allocated"
    istat = cuda_event_record(this%cptr, stream_cptr)
    IF(istat /= 0) STOP "acc_event_record failed"
#endif
END SUBROUTINE acc_event_record


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventCreate.
!> \param[in,out] this event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_create(this)
    TYPE(acc_event_type), &
      INTENT(INOUT)                          :: this

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat

    IF(C_ASSOCIATED(this%cptr)) STOP "acc_event_create: already allocated"
    istat = cuda_event_create(this%cptr)
    IF(istat /= 0 .OR. .NOT. C_ASSOCIATED(this%cptr)) STOP "acc_event_create: failed"
#endif
END SUBROUTINE acc_event_create


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventDestroy.
!> \param[in,out] this event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_destroy(this)
    TYPE(acc_event_type), &
      INTENT(INOUT)                          :: this

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) STOP "acc_event_destroy: event not allocated"
    istat = cuda_event_destroy(this%cptr)
    IF(istat /= 0) STOP "acc_event_destroy failed"
    this%cptr = C_NULL_PTR
#endif
END SUBROUTINE acc_event_destroy


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventQuery.
!> \param[in] this event
!> \retval res true if event has occured, false otherwise
!> \author  Ole Schuett
! *****************************************************************************
FUNCTION acc_event_query(this) RESULT(res)
    TYPE(acc_event_type), INTENT(IN)         :: this
    LOGICAL                                  :: res

#if ! defined (__ACC)
    res = .FALSE.
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat, has_occured
    IF(.NOT. C_ASSOCIATED(this%cptr)) STOP "acc_event_query: event not allocated"
    istat = cuda_event_query(this%cptr, has_occured)
    IF(istat /= 0) STOP "acc_event_query failed"
    res = (has_occured == 1)
#endif
END FUNCTION acc_event_query


! *****************************************************************************
!> \brief Fortran-wrapper for cudaEventSynchronize.
!> \param[in] this event
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_event_synchronize(this)
    TYPE(acc_event_type), INTENT(IN)  :: this

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    IF(.NOT. C_ASSOCIATED(this%cptr)) STOP "acc_event_synchronize: event not allocated"
    istat = cuda_event_synchronize(this%cptr)
    IF(istat < 0) STOP "acc_event_synchronize failed"
#endif
END SUBROUTINE acc_event_synchronize

END MODULE acc_event
