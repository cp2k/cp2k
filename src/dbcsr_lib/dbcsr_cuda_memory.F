!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-06
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-06
! *****************************************************************************
MODULE dbcsr_cuda_memory
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dbcsr_cuda_methods,              ONLY: dbcsr_cuda_dev_mem_get_type,&
       dbcsr_cuda_dev_mem_setup
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type,&
                                             dbcsr_cuda_mem_type_c4,&
                                             dbcsr_cuda_mem_type_c8,&
                                             dbcsr_cuda_mem_type_i4,&
                                             dbcsr_cuda_mem_type_i8,&
                                             dbcsr_cuda_mem_type_r4,&
                                             dbcsr_cuda_mem_type_r8
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_4,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE dbcsr_types,                     ONLY: dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8,&
                                             dbcsr_type_int_4,&
                                             dbcsr_type_int_8
  USE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_memory'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.


  PUBLIC :: dbcsr_cuda_dev_mem_alloc, dbcsr_cuda_dev_mem_dealloc,&
            dbcsr_cuda_dev_mem_realloc

  PUBLIC :: dbcsr_cuda_host_mem_alloc, dbcsr_cuda_host_mem_dealloc

  PUBLIC :: dbcsr_cuda_dev_mem_zero

  PUBLIC :: dbcsr_cuda_dev_mem_init
  PUBLIC :: dbcsr_cuda_dev_mem_new, dbcsr_cuda_dev_mem_release,&
            dbcsr_cuda_dev_mem_hold
  PUBLIC :: dbcsr_cuda_dev_mem_info

  INTERFACE dbcsr_cuda_dev_mem_alloc
     MODULE PROCEDURE dev_mem_alloc_any
     MODULE PROCEDURE dev_mem_alloc_i, dev_mem_alloc_l
     MODULE PROCEDURE dev_mem_alloc_r, dev_mem_alloc_d
     MODULE PROCEDURE dev_mem_alloc_c, dev_mem_alloc_z
  END INTERFACE

  INTERFACE dbcsr_cuda_dev_mem_realloc
     MODULE PROCEDURE dev_mem_realloc_any
     MODULE PROCEDURE dev_mem_realloc_i, dev_mem_realloc_l
     MODULE PROCEDURE dev_mem_realloc_r, dev_mem_realloc_d
     MODULE PROCEDURE dev_mem_realloc_c, dev_mem_realloc_z
  END INTERFACE

  INTERFACE dbcsr_cuda_dev_mem_dealloc
     MODULE PROCEDURE dev_mem_dealloc_any
     MODULE PROCEDURE dev_mem_dealloc_i, dev_mem_dealloc_l
     MODULE PROCEDURE dev_mem_dealloc_r, dev_mem_dealloc_d
     MODULE PROCEDURE dev_mem_dealloc_c, dev_mem_dealloc_z
  END INTERFACE

  INTERFACE dbcsr_cuda_host_mem_alloc
     MODULE PROCEDURE host_mem_alloc_i, host_mem_alloc_l
     MODULE PROCEDURE host_mem_alloc_r, host_mem_alloc_d
     MODULE PROCEDURE host_mem_alloc_c, host_mem_alloc_z
     MODULE PROCEDURE host_mem_alloc_i_2D, host_mem_alloc_l_2D
     MODULE PROCEDURE host_mem_alloc_r_2D, host_mem_alloc_d_2D
     MODULE PROCEDURE host_mem_alloc_c_2D, host_mem_alloc_z_2D
     MODULE PROCEDURE host_mem_alloc_i_4D, host_mem_alloc_l_4D
     MODULE PROCEDURE host_mem_alloc_r_4D, host_mem_alloc_d_4D
     MODULE PROCEDURE host_mem_alloc_c_4D, host_mem_alloc_z_4D
  END INTERFACE

  INTERFACE dbcsr_cuda_host_mem_dealloc
     MODULE PROCEDURE host_mem_dealloc_i, host_mem_dealloc_l
     MODULE PROCEDURE host_mem_dealloc_r, host_mem_dealloc_d
     MODULE PROCEDURE host_mem_dealloc_c, host_mem_dealloc_z
     MODULE PROCEDURE host_mem_dealloc_i_2D, host_mem_dealloc_l_2D
     MODULE PROCEDURE host_mem_dealloc_r_2D, host_mem_dealloc_d_2D
     MODULE PROCEDURE host_mem_dealloc_c_2D, host_mem_dealloc_z_2D
     MODULE PROCEDURE host_mem_dealloc_i_4D, host_mem_dealloc_l_4D
     MODULE PROCEDURE host_mem_dealloc_r_4D, host_mem_dealloc_d_4D
     MODULE PROCEDURE host_mem_dealloc_c_4D, host_mem_dealloc_z_4D
  END INTERFACE

  INTERFACE dbcsr_cuda_dev_mem_zero
     MODULE PROCEDURE dev_mem_zero_any
     MODULE PROCEDURE dev_mem_zero_i, dev_mem_zero_l
     MODULE PROCEDURE dev_mem_zero_r, dev_mem_zero_d
     MODULE PROCEDURE dev_mem_zero_c, dev_mem_zero_z
  END INTERFACE


#if defined (__DBCSR_CUDA)


  INTERFACE
     FUNCTION dbcsr_cuda_dev_mem_info(free, avail) RESULT (istat) BIND(C, name="dc_dev_mem_info_cu")
       USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(OUT)         :: free, avail
    INTEGER(KIND=C_INT)                      :: istat

    END FUNCTION dbcsr_cuda_dev_mem_info
  END INTERFACE


  INTERFACE
     FUNCTION cuda_dev_mem_alloc_cu(mem, n) RESULT (istat) BIND(C, name="dc_dev_mem_alloc")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: n
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_dev_mem_alloc_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_dev_mem_dealloc_cu(mem) RESULT (istat) BIND(C, name="dc_dev_mem_dealloc")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: mem
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_dev_mem_dealloc_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_dev_mem_realloc_cu(mem, n, old_n, memory_crunch)&
          RESULT (istat) BIND(C, name="dc_dev_mem_realloc")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: n, old_n
    INTEGER(KIND=C_INT), INTENT(INOUT)       :: memory_crunch
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_dev_mem_realloc_cu
  END INTERFACE


  INTERFACE
     FUNCTION cuda_host_mem_alloc_cu(mem, n, wc, port) RESULT (istat) BIND(C, name="dc_host_mem_alloc")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: n
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: wc, port
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_host_mem_alloc_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_host_mem_dealloc_cu(mem) RESULT (istat) bind(C, name="dc_host_mem_dealloc")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: mem
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_host_mem_dealloc_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_dev_zero_cu(dev_mem, offset, length) RESULT (istat) BIND(C, name="dc_memzero_cu")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: dev_mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: offset, length
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_dev_zero_cu
  END INTERFACE

#endif


CONTAINS

#if ! defined (__DBCSR_CUDA)
 FUNCTION dbcsr_cuda_dev_mem_info(free, avail) RESULT (istat)
    INTEGER, INTENT(OUT)                     :: free, avail
    INTEGER                                  :: istat

       free = 0
       avail = 0
       istat = -1
 END FUNCTION dbcsr_cuda_dev_mem_info
#endif
        

#include "dbcsr_cuda_memory_i.F"
#include "dbcsr_cuda_memory_l.F"
#include "dbcsr_cuda_memory_r.F"
#include "dbcsr_cuda_memory_d.F"
#include "dbcsr_cuda_memory_c.F"
#include "dbcsr_cuda_memory_z.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Encapsulated memory routines for the device

! *****************************************************************************
!> \brief
! *****************************************************************************
 SUBROUTINE host_mem_alloc_raw (host_mem_c_ptr, n_bytes, error)
    TYPE(C_PTR), INTENT(OUT)                 :: host_mem_c_ptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_alloc_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_alloc_cu(host_mem_c_ptr, &
               INT(n_bytes, KIND=C_SIZE_T), 0_c_int, 0_c_int)
    IF (istat /= 0) STOP "host_mem_alloc_raw: Could not allocate host pinned memory"
#else
    STOP "host_mem_alloc_raw: DBCSR_CUDA not compiled in."
#endif
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_alloc_raw


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_dev_mem_init (dev_mem, error)
    TYPE(dbcsr_cuda_mem_type), POINTER       :: dev_mem
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    NULLIFY (dev_mem)
  END SUBROUTINE dbcsr_cuda_dev_mem_init


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_dev_mem_new (dev_mem, data_type, error)
    TYPE(dbcsr_cuda_mem_type), POINTER       :: dev_mem
    INTEGER, INTENT(IN), OPTIONAL            :: data_type
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    ALLOCATE (dev_mem)
    dev_mem%data_type = 0
    dev_mem%refcount = 1
    IF (PRESENT (data_type)) THEN
       CALL dbcsr_cuda_dev_mem_setup(dev_mem, data_type)
    ENDIF
  END SUBROUTINE dbcsr_cuda_dev_mem_new


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_dev_mem_hold (dev_mem, error)
    TYPE(dbcsr_cuda_mem_type), POINTER       :: dev_mem
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    dev_mem%refcount = dev_mem%refcount +1
  END SUBROUTINE dbcsr_cuda_dev_mem_hold


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_dev_mem_release (dev_mem, error)
    TYPE(dbcsr_cuda_mem_type), POINTER       :: dev_mem
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    dev_mem%refcount = dev_mem%refcount - 1
    IF (dev_mem%refcount == 0) THEN
       CALL dev_mem_dealloc_any (dev_mem, error=error)
       DEALLOCATE (dev_mem)
       NULLIFY (dev_mem)
    ENDIF
  END SUBROUTINE dbcsr_cuda_dev_mem_release


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_alloc_any (dev_mem, n, error, stat)
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: dev_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_alloc_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set(routineN, error_handle, error)

    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_int_4)
       CALL dbcsr_cuda_dev_mem_alloc (dev_mem%d_i, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_int_8)
       CALL dbcsr_cuda_dev_mem_alloc (dev_mem%d_l, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_cuda_dev_mem_alloc (dev_mem%d_r, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_cuda_dev_mem_alloc (dev_mem%d_d, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_cuda_dev_mem_alloc (dev_mem%d_c, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_cuda_dev_mem_alloc (dev_mem%d_z, n,&
            error=error, stat=stat)
    CASE default
       STOP "dev_mem_alloc_any: Invalid data type"
    END SELECT

    IF (PRESENT (error)) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dev_mem_alloc_any


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_realloc_any (dev_mem, n, error, stat)
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: dev_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_realloc_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set(routineN, error_handle, error)

    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_int_4)
       CALL dbcsr_cuda_dev_mem_realloc (dev_mem%d_i, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_int_8)
       CALL dbcsr_cuda_dev_mem_realloc (dev_mem%d_l, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_cuda_dev_mem_realloc (dev_mem%d_r, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_cuda_dev_mem_realloc (dev_mem%d_d, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_cuda_dev_mem_realloc (dev_mem%d_c, n,&
            error=error, stat=stat)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_cuda_dev_mem_realloc (dev_mem%d_z, n,&
            error=error, stat=stat)
    CASE default
       STOP "dev_mem_realloc_any: Invalid data type"
    END SELECT

    IF (PRESENT (error)) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dev_mem_realloc_any


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_dealloc_any (dev_mem, error, stat)
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: dev_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_dealloc_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set(routineN, error_handle, error)

    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_int_4)
       CALL dbcsr_cuda_dev_mem_dealloc (dev_mem%d_i,&
            error=error, stat=stat)
    CASE (dbcsr_type_int_8)
       CALL dbcsr_cuda_dev_mem_dealloc (dev_mem%d_l,&
            error=error, stat=stat)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_cuda_dev_mem_dealloc (dev_mem%d_r,&
            error=error, stat=stat)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_cuda_dev_mem_dealloc (dev_mem%d_d,&
            error=error, stat=stat)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_cuda_dev_mem_dealloc (dev_mem%d_c,&
            error=error, stat=stat)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_cuda_dev_mem_dealloc (dev_mem%d_z,&
            error=error, stat=stat)
    CASE default
       STOP "dev_mem_dealloc_any: Invalid data type"
    END SELECT

    IF (PRESENT (error)) CALL dbcsr_error_stop(error_handle, error)
  END SUBROUTINE dev_mem_dealloc_any


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_zero_any (dev_mem, first, last, error)
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: dev_mem
    INTEGER, INTENT(IN)                      :: first, last
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_zero_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)

    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_int_4)
       CALL dbcsr_cuda_dev_mem_zero (dev_mem%d_i, first, last, error=error)
    CASE (dbcsr_type_int_8)
       CALL dbcsr_cuda_dev_mem_zero (dev_mem%d_l, first, last, error=error)
    CASE (dbcsr_type_real_4)
       CALL dbcsr_cuda_dev_mem_zero (dev_mem%d_r, first, last, error=error)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_cuda_dev_mem_zero (dev_mem%d_d, first, last, error=error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_cuda_dev_mem_zero (dev_mem%d_c, first, last, error=error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_cuda_dev_mem_zero (dev_mem%d_z, first, last, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type", __LINE__, error=error)
    END SELECT

    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_zero_any

END MODULE dbcsr_cuda_memory
