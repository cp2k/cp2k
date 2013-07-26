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
MODULE dbcsr_cuda_operations
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dbcsr_cuda_methods,              ONLY: dbcsr_cuda_dev_mem_get_type
  USE dbcsr_cuda_types,                ONLY: dbcsr_cuda_mem_type,&
                                             dbcsr_cuda_mem_type_c4,&
                                             dbcsr_cuda_mem_type_c8,&
                                             dbcsr_cuda_mem_type_i4,&
                                             dbcsr_cuda_mem_type_i8,&
                                             dbcsr_cuda_mem_type_r4,&
                                             dbcsr_cuda_mem_type_r8
  USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size_referenced,&
                                             dbcsr_data_get_type,&
                                             dbcsr_get_data,&
                                             dbcsr_get_data_p_c,&
                                             dbcsr_get_data_p_d,&
                                             dbcsr_get_data_p_s,&
                                             dbcsr_get_data_p_z
  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: int_4,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_type_complex_4,&
                                             dbcsr_type_complex_8,&
                                             dbcsr_type_int_4,&
                                             dbcsr_type_int_8,&
                                             dbcsr_type_real_4,&
                                             dbcsr_type_real_8
  USE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_operations'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.

  INTEGER, PARAMETER :: default_stream = 0


  PUBLIC :: dbcsr_cuda_cp_host_to_dev, dbcsr_cuda_cp_dev_to_host

  PUBLIC :: dbcsr_cuda_do_mm_stack


  INTERFACE dbcsr_cuda_do_mm_stack
     MODULE PROCEDURE do_mm_stack_any
     MODULE PROCEDURE do_mm_stack_r, do_mm_stack_d,&
                      do_mm_stack_c, do_mm_stack_z
  END INTERFACE

  INTERFACE dbcsr_cuda_cp_host_to_dev
     MODULE PROCEDURE cp_host_to_dev_any
     MODULE PROCEDURE cp_host_to_dev_i0, cp_host_to_dev_l0
     MODULE PROCEDURE cp_host_to_dev_r0, cp_host_to_dev_d0
     MODULE PROCEDURE cp_host_to_dev_c0, cp_host_to_dev_z0
     MODULE PROCEDURE cp_host_to_dev_i, cp_host_to_dev_l
     MODULE PROCEDURE cp_host_to_dev_r, cp_host_to_dev_d
     MODULE PROCEDURE cp_host_to_dev_c, cp_host_to_dev_z
  END INTERFACE

  INTERFACE dbcsr_cuda_cp_dev_to_host
     MODULE PROCEDURE cp_dev_to_host_any
     MODULE PROCEDURE cp_dev_to_host_i0, cp_dev_to_host_l0
     MODULE PROCEDURE cp_dev_to_host_r0, cp_dev_to_host_d0
     MODULE PROCEDURE cp_dev_to_host_c0, cp_dev_to_host_z0
     MODULE PROCEDURE cp_dev_to_host_i, cp_dev_to_host_l
     MODULE PROCEDURE cp_dev_to_host_r, cp_dev_to_host_d
     MODULE PROCEDURE cp_dev_to_host_c, cp_dev_to_host_z
  END INTERFACE


#if defined (__DBCSR_CUDA)
  INTERFACE
     FUNCTION cuda_memcpy_h2d_cu(host, dev, count, async_type, stream_id) RESULT (istat) &
          BIND(C, name="dc_memcpy_h2d_cu")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: host
    TYPE(C_PTR), VALUE                       :: dev
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: count
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: async_type, stream_id
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_memcpy_h2d_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_memcpy_d2h_cu(dev, host, count, async_type, stream_id) RESULT (istat) &
          BIND(C, name="dc_memcpy_d2h_cu")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: dev
    TYPE(C_PTR), VALUE                       :: host
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: count
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: async_type, stream_id
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_memcpy_d2h_cu
  END INTERFACE


  INTERFACE
     FUNCTION cuda_do_stack_cu(param_stack, stack_size, nparams,&
          data_type,&
          a_data, b_data, c_data, c_locks, m_max, n_max, k_max, def_mnk,&
          stream_id)&
          RESULT (istat) &
          BIND(C, name="dc_do_stack_cu")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: param_stack
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: stack_size, nparams, data_type
    TYPE(C_PTR), INTENT(IN), VALUE           :: a_data, b_data
    TYPE(C_PTR), VALUE                       :: c_data, c_locks
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: m_max, n_max, k_max, def_mnk, &
                                                stream_id
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_do_stack_cu
  END INTERFACE
#endif


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111!!!!!!!!!!!!!!!!!!
! Encapsulated data routines

  SUBROUTINE do_mm_stack_any (param_stack, stack_size, nparams,&
       a_data, b_data, c_data, c_locks, m_max, n_max, k_max, def_mnk, stream,&
       error)
    TYPE(dbcsr_cuda_mem_type_i4), INTENT(IN) :: param_stack
    INTEGER, INTENT(IN)                      :: stack_size, nparams
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: a_data, b_data
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: c_data
    TYPE(dbcsr_cuda_mem_type_i4), INTENT(IN) :: c_locks
    INTEGER, INTENT(IN)                      :: m_max, n_max, k_max
    LOGICAL, INTENT(IN)                      :: def_mnk
    INTEGER, INTENT(IN), OPTIONAL            :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_mm_stack_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (careful_mod) THEN
       CALL dbcsr_assert (dbcsr_cuda_dev_mem_get_type (a_data), "EQ",&
            dbcsr_cuda_dev_mem_get_type (c_data),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Different data types.", __LINE__, error=error)
       CALL dbcsr_assert (dbcsr_cuda_dev_mem_get_type (b_data), "EQ",&
            dbcsr_cuda_dev_mem_get_type (c_data),&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Different data types.", __LINE__, error=error)
    ENDIF
    SELECT CASE (dbcsr_cuda_dev_mem_get_type (c_data))
    CASE (dbcsr_type_real_4)
       CALL dbcsr_cuda_do_mm_stack (param_stack, stack_size, nparams,&
             a_data%d_r, b_data%d_r, c_data%d_r,&
             c_locks, m_max, n_max, k_max, def_mnk, stream=stream, error=error)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_cuda_do_mm_stack (param_stack, stack_size, nparams,&
             a_data%d_d, b_data%d_d, c_data%d_d,&
             c_locks, m_max, n_max, k_max, def_mnk, stream=stream, error=error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_cuda_do_mm_stack (param_stack, stack_size, nparams,&
             a_data%d_c, b_data%d_c, c_data%d_c,&
             c_locks, m_max, n_max, k_max, def_mnk, stream=stream, error=error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_cuda_do_mm_stack (param_stack, stack_size, nparams,&
             a_data%d_z, b_data%d_z, c_data%d_z,&
             c_locks, m_max, n_max, k_max, def_mnk, stream=stream, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type", __LINE__, error=error)
    END SELECT
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE do_mm_stack_any


  SUBROUTINE cp_host_to_dev_any (host_mem, dev_mem, async, stream, error)
    TYPE(dbcsr_data_obj), INTENT(IN)         :: host_mem
    TYPE(dbcsr_cuda_mem_type), INTENT(INOUT) :: dev_mem
    LOGICAL, INTENT(IN), OPTIONAL            :: async
    INTEGER, INTENT(IN), OPTIONAL            :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_host_to_dev_any', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, n

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)

    CALL dbcsr_assert (dbcsr_data_get_type (host_mem), "EQ",&
         dbcsr_cuda_dev_mem_get_type (dev_mem),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Host and device data types not compatible.",&
         __LINE__, error=error)

    n = dbcsr_data_get_size_referenced (host_mem)
    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_real_4)
       CALL dbcsr_cuda_cp_host_to_dev (dbcsr_get_data_p_s(host_mem),&
            dev_mem%d_r, n, async=async, stream=stream,&
            error=error)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_cuda_cp_host_to_dev (dbcsr_get_data_p_d(host_mem),&
            dev_mem%d_d, n, async=async, stream=stream,&
            error=error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_cuda_cp_host_to_dev (dbcsr_get_data_p_c(host_mem),&
            dev_mem%d_c, n, async=async, stream=stream,&
            error=error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_cuda_cp_host_to_dev (dbcsr_get_data_p_z(host_mem),&
            dev_mem%d_z, n, async=async, stream=stream,&
            error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type", __LINE__, error=error)
    END SELECT

    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cp_host_to_dev_any


  SUBROUTINE cp_dev_to_host_any (dev_mem, host_mem, async, stream, error)
    TYPE(dbcsr_cuda_mem_type), INTENT(IN)    :: dev_mem
    TYPE(dbcsr_data_obj), INTENT(INOUT)      :: host_mem
    LOGICAL, INTENT(IN), OPTIONAL            :: async
    INTEGER, INTENT(IN), OPTIONAL            :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dev_to_host_any', &
      routineP = moduleN//':'//routineN

    COMPLEX(kind=real_4), DIMENSION(:), &
      POINTER                                :: d_c4
    COMPLEX(kind=real_8), DIMENSION(:), &
      POINTER                                :: d_c8
    INTEGER                                  :: error_handle, n
    REAL(kind=real_4), DIMENSION(:), POINTER :: d_r4
    REAL(kind=real_8), DIMENSION(:), POINTER :: d_r8

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)

    CALL dbcsr_assert (dbcsr_data_get_type (host_mem), "EQ",&
         dbcsr_cuda_dev_mem_get_type (dev_mem),&
         dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
         "Host and device data types not compatible.",&
         __LINE__, error=error)

    n = dbcsr_data_get_size_referenced (host_mem)
    SELECT CASE (dbcsr_cuda_dev_mem_get_type (dev_mem))
    CASE (dbcsr_type_real_4)
       CALL dbcsr_get_data (host_mem, d_r4)
       CALL dbcsr_cuda_cp_dev_to_host (dev_mem%d_r,&
            d_r4, n, async=async, stream=stream, error=error)
    CASE (dbcsr_type_real_8)
       CALL dbcsr_get_data (host_mem, d_r8)
       CALL dbcsr_cuda_cp_dev_to_host (dev_mem%d_d,&
            d_r8, n, async=async, stream=stream, error=error)
    CASE (dbcsr_type_complex_4)
       CALL dbcsr_get_data (host_mem, d_c4)
       CALL dbcsr_cuda_cp_dev_to_host (dev_mem%d_c,&
            d_c4, n, async=async, stream=stream, error=error)
    CASE (dbcsr_type_complex_8)
       CALL dbcsr_get_data (host_mem, d_c8)
       CALL dbcsr_cuda_cp_dev_to_host (dev_mem%d_z,&
            d_c8, n, async=async, stream=stream, error=error)
    CASE default
       CALL dbcsr_assert (.FALSE.,&
            dbcsr_fatal_level, dbcsr_wrong_args_error, routineN,&
            "Invalid data type", __LINE__, error=error)
    END SELECT

    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cp_dev_to_host_any


#include "dbcsr_cuda_operations_i.F"
#include "dbcsr_cuda_operations_l.F"
#include "dbcsr_cuda_operations_r.F"
#include "dbcsr_cuda_operations_d.F"
#include "dbcsr_cuda_operations_c.F"
#include "dbcsr_cuda_operations_z.F"

END MODULE dbcsr_cuda_operations
