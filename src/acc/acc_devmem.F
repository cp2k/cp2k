!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014 CP2K developers group                           !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
! *****************************************************************************
MODULE acc_devmem
  USE ISO_C_BINDING
  USE acc_kinds,                       ONLY: int_4,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE acc_stream,                      ONLY: acc_stream_associated,&
                                             acc_stream_cptr,&
                                             acc_stream_synchronize,&
                                             acc_stream_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_devmem'


  PUBLIC :: acc_devmem_type
  PUBLIC :: acc_devmem_allocate_bytes, acc_devmem_deallocate
  PUBLIC :: acc_devmem_setzero_bytes
  PUBLIC :: acc_devmem_allocated, acc_devmem_info
  PUBLIC :: acc_devmem_dev2host, acc_devmem_host2dev
  PUBLIC :: acc_devmem_size_in_bytes
  PUBLIC :: acc_devmem_ensure_size_bytes
  PUBLIC :: acc_devmem_cptr

 INTERFACE acc_devmem_dev2host
  MODULE PROCEDURE dev2host_i4_1D
  MODULE PROCEDURE dev2host_i8_1D
  MODULE PROCEDURE dev2host_r4_1D
  MODULE PROCEDURE dev2host_r8_1D
  MODULE PROCEDURE dev2host_c4_1D
  MODULE PROCEDURE dev2host_c8_1D
 END INTERFACE acc_devmem_dev2host

 INTERFACE acc_devmem_host2dev
  MODULE PROCEDURE host2dev_i4_1D
  MODULE PROCEDURE host2dev_i8_1D
  MODULE PROCEDURE host2dev_r4_1D
  MODULE PROCEDURE host2dev_r8_1D
  MODULE PROCEDURE host2dev_c4_1D
  MODULE PROCEDURE host2dev_c8_1D
  MODULE PROCEDURE host2dev_i4_2D
  MODULE PROCEDURE host2dev_i8_2D
  MODULE PROCEDURE host2dev_r4_2D
  MODULE PROCEDURE host2dev_r8_2D
  MODULE PROCEDURE host2dev_c4_2D
  MODULE PROCEDURE host2dev_c8_2D
 END INTERFACE acc_devmem_host2dev

 TYPE acc_devmem_type
   PRIVATE
   INTEGER                      :: size_in_bytes = -1
   TYPE(C_PTR)                  :: cptr = C_NULL_PTR
 END TYPE acc_devmem_type

#if  defined (__ACC)

  INTERFACE
     FUNCTION dc_dev_mem_info_cu(free, avail) RESULT (istat) BIND(C, name="acc_dev_mem_info")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(OUT)         :: free, avail
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION dc_dev_mem_info_cu
  END INTERFACE

  INTERFACE
     FUNCTION dc_dev_mem_alloc(mem, n) RESULT (istat) BIND(C, name="acc_dev_mem_allocate")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: n
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_dev_mem_alloc
  END INTERFACE

  INTERFACE
     FUNCTION dc_dev_mem_dealloc(mem) RESULT (istat) BIND(C, name="acc_dev_mem_deallocate")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: mem
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_dev_mem_dealloc
  END INTERFACE

  INTERFACE
     FUNCTION dc_memzero_cu(this, offset, length, stream_ptr) RESULT (istat) BIND(C, name="acc_memset_zero")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: this
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: offset, length
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_memzero_cu
  END INTERFACE

  INTERFACE
     FUNCTION dc_memcpy_h2d_cu(host, dev, count, stream_ptr) RESULT (istat) &
          BIND(C, name="acc_memcpy_h2d")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: host
    TYPE(C_PTR), VALUE                       :: dev
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: count
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION dc_memcpy_h2d_cu
  END INTERFACE

  INTERFACE
     FUNCTION dc_memcpy_d2h_cu(dev, host, count, stream_ptr) RESULT (istat) &
          BIND(C, name="acc_memcpy_d2h")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: dev
    TYPE(C_PTR), VALUE                       :: host
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: count
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION dc_memcpy_d2h_cu
  END INTERFACE

  INTERFACE
     FUNCTION dc_memcpy_d2d_cu(dev_src, dev_dst, count, stream_ptr) RESULT (istat) &
          BIND(C, name="acc_memcpy_d2d")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: dev_src
    TYPE(C_PTR), VALUE                       :: dev_dst
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: count
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION dc_memcpy_d2d_cu
  END INTERFACE

#endif

CONTAINS


! *****************************************************************************
!> \brief Ensures that given devmem has at least the requested size.
!> \param stream on which zeroing and memcopying is performed
!> \param requested_size_in_bytes
!> \param nocopy if after growin old content should NOT be copied over. default:false.
!> \param zero_pad if after growing the new memory should be zeroed. default: false.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_ensure_size_bytes(this, stream, requested_size_in_bytes, nocopy, zero_pad)
    TYPE(acc_devmem_type), &
      INTENT(INOUT)                          :: this
    TYPE(acc_stream_type), INTENT(IN) :: stream
    INTEGER, INTENT(IN)                      :: requested_size_in_bytes
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else

    LOGICAL                                  :: my_nocopy, my_zero_pad
    TYPE(C_PTR)                              :: old_cptr, new_cptr, stream_cptr
    INTEGER                                  :: new_size, old_size, istat

    IF(this%size_in_bytes < 0) STOP "acc_devmem_ensure_size_bytes: not allocated"
    IF(.NOT. acc_stream_associated(stream)) &
         STOP "acc_devmem_ensure_size_bytes: stream not associated"

    IF(this%size_in_bytes < requested_size_in_bytes) THEN
       !WRITE (*,*) "acc_devmem_ensure_size_bytes: growing dev_mem to: ", data_size

       new_size = requested_size_in_bytes
       old_size = this%size_in_bytes
       old_cptr = this%cptr

       new_cptr = C_NULL_PTR
       istat = dc_dev_mem_alloc(new_cptr, INT(new_size, KIND=C_SIZE_T))
       IF(istat /= 0) STOP "acc_devmem_ensure_size_bytes: alloc failed"

       this%cptr = new_cptr
       this%size_in_bytes = requested_size_in_bytes

       my_zero_pad = .FALSE.
       IF(PRESENT(zero_pad)) my_zero_pad=zero_pad
       IF(my_zero_pad) &
          CALL acc_devmem_setzero_bytes(this, first_byte=old_size, stream=stream)

       my_nocopy = .FALSE.
       IF(PRESENT(nocopy)) my_nocopy=nocopy
       IF(.NOT. my_nocopy) THEN
          stream_cptr = acc_stream_cptr(stream)
          istat = dc_memcpy_d2d_cu(old_cptr, new_cptr, INT(old_size, KIND=C_SIZE_T), stream_cptr)
          IF(istat /= 0) STOP "acc_devmem_ensure_size_bytes: memcpy failed"
       ENDIF

       CALL acc_stream_synchronize(stream)
       istat = dc_dev_mem_dealloc(old_cptr)
       IF(istat /= 0) STOP "acc_devmem_ensure_size_bytes: dealloc failed"

    ENDIF
#endif
  END SUBROUTINE acc_devmem_ensure_size_bytes

! *****************************************************************************
!> \brief Returns a logical, which indicates if the given devmem is allocated.
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION acc_devmem_allocated(this) RESULT(res)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    LOGICAL                                  :: res

     res = this%size_in_bytes >= 0
  END FUNCTION acc_devmem_allocated


! *****************************************************************************
!> \brief Returns size of given devmem in terms of item count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION acc_devmem_size_in_bytes(this) RESULT(res)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    INTEGER                                  :: res

     IF(this%size_in_bytes < 0) STOP "acc_devmem_len: not allocated"
     res = this%size_in_bytes
  END FUNCTION acc_devmem_size_in_bytes



! *****************************************************************************
!> \brief Returns C-pointer to data of give devmem.
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION acc_devmem_cptr(this) RESULT(res)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    TYPE(C_PTR)                              :: res

     IF(this%size_in_bytes < 0) STOP "acc_devmem_cptr: not allocated"
     res = this%cptr
  END FUNCTION acc_devmem_cptr


! *****************************************************************************
!> \brief Fortran-wrapper for cudaMemGetInfo.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE acc_devmem_info(free, avail)
    INTEGER, INTENT(OUT)                     :: free, avail

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    INTEGER                                  :: istat
    istat = dc_dev_mem_info_cu(free, avail)
    IF (istat /= 0) STOP "acc_devmem_info: failed."
#endif
 END SUBROUTINE acc_devmem_info


! *****************************************************************************
!> \brief Allocates a given devmem.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_allocate_bytes(this, size_in_bytes)
    TYPE(acc_devmem_type), INTENT(INOUT)     :: this
    INTEGER, INTENT(IN)                      :: size_in_bytes

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_devmem_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF(this%size_in_bytes >= 0) STOP "acc_devmem_alloc: allready allocated"
    this%size_in_bytes = size_in_bytes
    IF(size_in_bytes>0) THEN
      istat = dc_dev_mem_alloc(this%cptr, INT(this%size_in_bytes, KIND=C_SIZE_T))
      IF (istat /= 0) STOP "acc_devmem_allocate: failed."
    END IF
#endif
  END SUBROUTINE acc_devmem_allocate_bytes


! *****************************************************************************
!> \brief Deallocates a given devmem.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_deallocate(this)
    TYPE(acc_devmem_type), INTENT(INOUT) :: this

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_devmem_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat

    IF(this%size_in_bytes < 0) STOP "acc_devmem_deallocate: double free"
    IF (this%size_in_bytes > 0) THEN
       istat = dc_dev_mem_dealloc(this%cptr)
       IF (istat /= 0) STOP "acc_devmem_deallocate: failed"
    END IF

    this%size_in_bytes = -1

#endif
  END SUBROUTINE acc_devmem_deallocate


! *****************************************************************************
!> \brief Sets entries in given devmem to zero, asynchronously.
!> \param first begin of region to zero, defaults to 1 if not given.
!> \param last end of region to zero, defaults to size if not given.
!> \param stream on which zeroing is performed.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_devmem_setzero_bytes(this, first_byte, last_byte, stream)
    TYPE(acc_devmem_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN), OPTIONAL               :: first_byte, last_byte
    TYPE(acc_stream_type), INTENT(IN)    :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_devmem_setzero', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    INTEGER(KIND=C_SIZE_T)                   :: length, offset
    TYPE(C_PTR)                              :: stream_cptr

    offset = 0
    IF(PRESENT(first_byte)) offset = first_byte-1
    length = this%size_in_bytes - offset
    IF(PRESENT(last_byte)) length = last_byte - offset

    stream_cptr = acc_stream_cptr(stream)

    IF (length > 0) THEN
       istat = dc_memzero_cu(this%cptr, offset, length, stream_cptr)
       IF (istat /= 0) STOP "acc_devmem_setzero: failed"
    ENDIF
#endif
  END SUBROUTINE acc_devmem_setzero_bytes



! *****************************************************************************
!> \brief Helper-routine performing actuall host2dev transfers.
!> \author  Ole Schuett
! *****************************************************************************
#if defined (__ACC)
 SUBROUTINE host2dev_raw(this, hostmem_cptr, n_bytes, stream)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    TYPE(C_PTR)                              :: hostmem_cptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(acc_stream_type), INTENT(IN)        :: stream

    CHARACTER(len=*), PARAMETER :: routineN = 'host2dev_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    IF(this%size_in_bytes < n_bytes) STOP "acc_devmem_host2dev: devmem too small"

    stream_cptr = acc_stream_cptr(stream)
    IF (n_bytes > 0) THEN
       istat = dc_memcpy_h2d_cu(hostmem_cptr, this%cptr,&
            INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
       IF(istat /= 0) STOP "acc_devmem_host2dev: failed"
    ENDIF
  END SUBROUTINE host2dev_raw
#endif


! *****************************************************************************
!> \brief Helper-routine performing actual dev2host transfers.
!> \author  Ole Schuett
! *****************************************************************************
#if defined (__ACC)
  SUBROUTINE dev2host_raw(this, hostmem_cptr, n_bytes, stream)
    TYPE(acc_devmem_type), INTENT(IN)        :: this
    TYPE(C_PTR)                              :: hostmem_cptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(acc_stream_type), INTENT(IN)        :: stream

    CHARACTER(len=*), PARAMETER :: routineN = 'dev2host_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    IF(this%size_in_bytes < n_bytes) STOP "acc_devmem_dev2host: this too small"

    stream_cptr = acc_stream_cptr(stream)

    IF (n_bytes > 0) THEN
       istat = dc_memcpy_d2h_cu(this%cptr, hostmem_cptr,&
            INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
       IF (istat /= 0) STOP "acc_devmem_dev2host: failed"
    ENDIF
  END SUBROUTINE dev2host_raw
#endif


#include "acc_devmem_i4.f90"
#include "acc_devmem_i8.f90"
#include "acc_devmem_r4.f90"
#include "acc_devmem_r8.f90"
#include "acc_devmem_c4.f90"
#include "acc_devmem_c8.f90"

END MODULE acc_devmem
