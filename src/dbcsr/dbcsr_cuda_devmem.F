!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Urban Borstnik, Ole Schuett
! *****************************************************************************
MODULE dbcsr_cuda_devmem
#if defined (__DBCSR_CUDA)
  USE ISO_C_BINDING
  USE dbcsr_cuda_stream,               ONLY: dbcsr_cuda_stream_cptr
#endif

  USE dbcsr_cuda_stream,               ONLY: dbcsr_cuda_stream_type,&
                                             dbcsr_cuda_stream_synchronize,&
                                             dbcsr_cuda_stream_associated
  USE dbcsr_cuda_device,               ONLY: dbcsr_cuda_device_sync
  USE dbcsr_datatype,                  ONLY: dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_int_4,&
                                             dbcsr_type_int_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_4,int_8,&
                                             int_4_size,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE dbcsr_message_passing,           ONLY: mp_abort
  USE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_devmem'

  LOGICAL, PARAMETER, PRIVATE :: careful_mod=.FALSE.


  PUBLIC :: dbcsr_cuda_devmem_type
  PUBLIC :: dbcsr_cuda_devmem_allocate, dbcsr_cuda_devmem_deallocate
  PUBLIC :: dbcsr_cuda_devmem_setzero
  PUBLIC :: dbcsr_cuda_devmem_allocated, dbcsr_cuda_devmem_info
  PUBLIC :: dbcsr_cuda_devmem_dev2host, dbcsr_cuda_devmem_host2dev
  PUBLIC :: dbcsr_cuda_devmem_size, dbcsr_cuda_devmem_datatype
  PUBLIC :: dbcsr_cuda_devmem_ensure_size

 INTERFACE dbcsr_cuda_devmem_dev2host
  MODULE PROCEDURE dev2host_i4_1D
  MODULE PROCEDURE dev2host_i8_1D
  MODULE PROCEDURE dev2host_r4_1D
  MODULE PROCEDURE dev2host_r8_1D
  MODULE PROCEDURE dev2host_c4_1D
  MODULE PROCEDURE dev2host_c8_1D
 END INTERFACE dbcsr_cuda_devmem_dev2host

 INTERFACE dbcsr_cuda_devmem_host2dev
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
 END INTERFACE dbcsr_cuda_devmem_host2dev


#if ! defined (__DBCSR_CUDA)
  TYPE dbcsr_cuda_devmem_type
     PRIVATE
     INTEGER                      :: datatype = 0
     INTEGER                      :: size = 0 !number of datatype items
  END TYPE dbcsr_cuda_devmem_type

#else

  PUBLIC :: dbcsr_cuda_devmem_cptr

  TYPE dbcsr_cuda_devmem_type
     PRIVATE
     INTEGER                      :: datatype = 0
     INTEGER                      :: size = 0 !number of datatype items
     TYPE(C_PTR)                  :: cptr = C_NULL_PTR
  END TYPE dbcsr_cuda_devmem_type

  INTERFACE
     FUNCTION dc_dev_mem_info_cu(free, avail) RESULT (istat) BIND(C, name="dc_dev_mem_info_cu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(OUT)         :: free, avail
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION dc_dev_mem_info_cu
  END INTERFACE

  INTERFACE
     FUNCTION dc_dev_mem_alloc(mem, n) RESULT (istat) BIND(C, name="dc_dev_mem_alloc")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: n
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_dev_mem_alloc
  END INTERFACE

  INTERFACE
     FUNCTION dc_dev_mem_dealloc(mem) RESULT (istat) BIND(C, name="dc_dev_mem_dealloc")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: mem
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_dev_mem_dealloc
  END INTERFACE

  INTERFACE
     FUNCTION dc_memzero_cu(this, offset, length, stream_ptr) RESULT (istat) BIND(C, name="dc_memzero_cu")
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
          BIND(C, name="dc_memcpy_h2d_cu")
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
          BIND(C, name="dc_memcpy_d2h_cu")
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
          BIND(C, name="dc_memcpy_d2d_cu")
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
!> \param data_size size request in terms of item-count (not bytes!)
!> \param nocopy if after growin old content should NOT be copied over. default:false.
!> \param zero_pad if after growing the new memory should be zeroed. default: false.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_devmem_ensure_size(this, stream, data_size, nocopy, zero_pad, error)
    TYPE(dbcsr_cuda_devmem_type), &
      INTENT(INOUT)                          :: this
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    INTEGER, INTENT(IN)                      :: data_size
    LOGICAL, INTENT(IN), OPTIONAL            :: nocopy, zero_pad
    TYPE(dbcsr_error_type), INTENT(inout)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else

    LOGICAL                                  :: my_nocopy, my_zero_pad
    INTEGER(KIND=C_SIZE_T)                   :: old_nbytes, new_nbytes
    TYPE(C_PTR)                              :: old_cptr, new_cptr, stream_cptr
    INTEGER                                  :: old_size, istat

    IF(this%datatype <= 0) CALL mp_abort("dbcsr_cuda_devmem_ensure_size: not allocated")
    IF(.NOT. dbcsr_cuda_stream_associated(stream)) &
         CALL mp_abort("dbcsr_cuda_devmem_ensure_size: stream not associated")

    IF(this%size < data_size) THEN
       !WRITE (*,*) "dbcsr_cuda_devmem_ensure_size: growing dev_mem to: ", data_size

       old_size  = this%size
       old_cptr  = this%cptr
       old_nbytes = INT(this%size*sizeof(this%datatype), KIND=C_SIZE_T)
       new_nbytes = INT(data_size*sizeof(this%datatype), KIND=C_SIZE_T)

       new_cptr = C_NULL_PTR
       istat = dc_dev_mem_alloc(new_cptr, new_nbytes)
       IF(istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_ensure_size: alloc failed")

       this%cptr = new_cptr
       this%size = data_size

       my_zero_pad = .FALSE.
       IF(PRESENT(zero_pad)) my_zero_pad=zero_pad
       IF(my_zero_pad) &
          CALL dbcsr_cuda_devmem_setzero(this, first=old_size, stream=stream, error=error)

       my_nocopy = .FALSE.
       IF(PRESENT(nocopy)) my_nocopy=nocopy
       IF(.NOT. my_nocopy) THEN
          stream_cptr = dbcsr_cuda_stream_cptr(stream)
          istat = dc_memcpy_d2d_cu(old_cptr, new_cptr, old_nbytes, stream_cptr)
          IF(istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_ensure_size: memcpy failed")
       ENDIF

       CALL dbcsr_cuda_stream_synchronize(stream, error)
       istat = dc_dev_mem_dealloc(old_cptr)
       IF(istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_ensure_size: dealloc failed")

    ENDIF
#endif
  END SUBROUTINE dbcsr_cuda_devmem_ensure_size

! *****************************************************************************
!> \brief Returns a logical, which indicates if the given devmem is allocated.
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_cuda_devmem_allocated(this) RESULT(res)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    LOGICAL                                  :: res

     res = this%datatype > 0
  END FUNCTION dbcsr_cuda_devmem_allocated


! *****************************************************************************
!> \brief Returns size of given devmem in terms of item count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_cuda_devmem_size(this) RESULT(res)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    INTEGER                                  :: res

     IF(this%datatype <= 0) CALL mp_abort("dbcsr_cuda_devmem_len: not allocated")
     res = this%size
  END FUNCTION dbcsr_cuda_devmem_size


! *****************************************************************************
!> \brief Returns datatype of given devmem, constans are defined in dbcsr_datatype.
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_cuda_devmem_datatype(this) RESULT(res)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    INTEGER                                  :: res

     IF(this%datatype <= 0) CALL mp_abort("dbcsr_cuda_devmem_datatype: not allocated")
     res = this%datatype
  END FUNCTION dbcsr_cuda_devmem_datatype


! *****************************************************************************
!> \brief Returns C-pointer to data of give devmem.
!> \author  Ole Schuett
! *****************************************************************************
#if defined (__DBCSR_CUDA)
  FUNCTION dbcsr_cuda_devmem_cptr(this) RESULT(res)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    TYPE(C_PTR)                              :: res

     IF(this%datatype <= 0) CALL mp_abort("dbcsr_cuda_devmem_cptr: not allocated")
     res = this%cptr
  END FUNCTION dbcsr_cuda_devmem_cptr
#endif


! *****************************************************************************
!> \brief Helper-routine, returns size of given datatype in terms of bytes.
!> \author  Ole Schuett
! *****************************************************************************
  FUNCTION sizeof(datatype) RESULT(size)
    INTEGER, INTENT(IN)                      :: datatype
    INTEGER                                  :: size

    SELECT CASE (datatype)
    CASE (dbcsr_type_int_4)
       size = int_4_size
    CASE (dbcsr_type_int_8)
       size = int_8_size
    CASE (dbcsr_type_real_4)
       size = real_4_size
    CASE (dbcsr_type_real_8)
       size = real_8_size
    CASE (dbcsr_type_complex_4)
       size = (2*real_4_size)
    CASE (dbcsr_type_complex_8)
       size = (2*real_8_size)
    CASE default
       CALL mp_abort("dbcsr_cuda_memory: Invalid data type")
    END SELECT
  END FUNCTION sizeof


! *****************************************************************************
!> \brief Fortran-wrapper for cudaMemGetInfo.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_cuda_devmem_info(free, avail, error)
    INTEGER, INTENT(OUT)                     :: free, avail
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    INTEGER                                  :: istat
    istat = dc_dev_mem_info_cu(free, avail)
    IF (istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_info: failed.")
#endif
 END SUBROUTINE dbcsr_cuda_devmem_info


! *****************************************************************************
!> \brief Allocates a given devmem.
!> \param datatype permissible constants are defined in dbcsr_datatype.
!> \param n amount of space to allocate in terms of item-count (not bytes!).
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_devmem_allocate(this, datatype, n, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN)                      :: datatype
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_devmem_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

    IF(this%datatype > 0) CALL mp_abort("dbcsr_cuda_devmem_alloc: allready allocated")
    this%datatype = datatype
    CALL dbcsr_error_set(routineN, error_handle, error)
    this%size = n
    IF(n>0) THEN
      istat = dc_dev_mem_alloc(this%cptr,&
         INT(this%size*sizeof(datatype), KIND=C_SIZE_T))
      IF (istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_allocate: failed.")
    END IF
    CALL dbcsr_error_stop(error_handle, error)
#endif
  END SUBROUTINE dbcsr_cuda_devmem_allocate


! *****************************************************************************
!> \brief Deallocates a given devmem.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_devmem_deallocate(this, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(INOUT) :: this
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_devmem_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

    CALL dbcsr_error_set(routineN, error_handle, error)
    IF(this%datatype <=0) CALL mp_abort("dbcsr_cuda_devmem_deallocate: double free")
    IF (this%size > 0) THEN
       istat = dc_dev_mem_dealloc(this%cptr)
       IF (istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_deallocate: failed")
    END IF

    this%size = 0
    this%datatype = 0

    CALL dbcsr_error_stop(error_handle, error)
#endif
  END SUBROUTINE dbcsr_cuda_devmem_deallocate


! *****************************************************************************
!> \brief Sets entries in given devmem to zero, asynchronously.
!> \param first begin of region to zero, defaults to 1 if not given.
!> \param last end of region to zero, defaults to size if not given.
!> \param stream on which zeroing is performed.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_devmem_setzero(this, first, last, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(INOUT) :: this
    INTEGER, INTENT(IN), OPTIONAL               :: first, last
    TYPE(dbcsr_cuda_stream_type), INTENT(IN)    :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)       :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_devmem_setzero', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, s, istat
    INTEGER(KIND=C_SIZE_T)                   :: len, off
    TYPE(C_PTR)                              :: stream_cptr

    CALL dbcsr_error_set (routineN, error_handle, error)

    s = sizeof(this%datatype)
    off = 0
    IF(PRESENT(first)) off = (first-1) * s
    len = this%size*s - off
    IF(PRESENT(last)) len = last*s - off

    stream_cptr = dbcsr_cuda_stream_cptr(stream)

    IF (len > 0) THEN
       istat = dc_memzero_cu(this%cptr, off, len, stream_cptr)
       IF (istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_setzero: failed")
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
#endif
  END SUBROUTINE dbcsr_cuda_devmem_setzero



! *****************************************************************************
!> \brief Helper-routine performing actuall host2dev transfers.
!> \author  Ole Schuett
! *****************************************************************************
#if defined (__DBCSR_CUDA)
 SUBROUTINE host2dev_raw(this, hostmem_cptr, n_bytes, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    TYPE(C_PTR)                              :: hostmem_cptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'host2dev_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    TYPE(C_PTR)                              :: stream_cptr

    IF(careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)
    IF(this%size*sizeof(this%datatype) < n_bytes) CALL mp_abort("dbcsr_cuda_devmem_host2dev: devmem too small")

    stream_cptr = dbcsr_cuda_stream_cptr(stream)
    IF (n_bytes > 0) THEN
       istat = dc_memcpy_h2d_cu(hostmem_cptr, this%cptr,&
            INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
       IF(istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_host2dev: failed")
    ENDIF
    IF(careful_mod) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host2dev_raw
#endif


! *****************************************************************************
!> \brief Helper-routine performing actual dev2host transfers.
!> \author  Ole Schuett
! *****************************************************************************
#if defined (__DBCSR_CUDA)
  SUBROUTINE dev2host_raw(this, hostmem_cptr, n_bytes, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    TYPE(C_PTR)                              :: hostmem_cptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dev2host_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    TYPE(C_PTR)                              :: stream_cptr

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF(this%size*sizeof(this%datatype) < n_bytes) CALL mp_abort("dbcsr_cuda_devmem_dev2host: this too small")

    stream_cptr = dbcsr_cuda_stream_cptr(stream)

    IF (n_bytes > 0) THEN
       istat = dc_memcpy_d2h_cu(this%cptr, hostmem_cptr,&
            INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
       IF (istat /= 0) CALL mp_abort("dbcsr_cuda_devmem_dev2host: failed")
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev2host_raw
#endif


#include "dbcsr_cuda_devmem_i4.f90"
#include "dbcsr_cuda_devmem_i8.f90"
#include "dbcsr_cuda_devmem_r4.f90"
#include "dbcsr_cuda_devmem_r8.f90"
#include "dbcsr_cuda_devmem_c4.f90"
#include "dbcsr_cuda_devmem_c8.f90"

END MODULE dbcsr_cuda_devmem
