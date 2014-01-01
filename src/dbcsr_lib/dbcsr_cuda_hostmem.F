!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE dbcsr_cuda_hostmem
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
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
  USE dbcsr_message_passing,           ONLY: mp_abort

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_hostmem'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.


  PUBLIC :: dbcsr_cuda_hostmem_allocate, dbcsr_cuda_hostmem_deallocate


  INTERFACE dbcsr_cuda_hostmem_allocate
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_i,    dbcsr_cuda_hostmem_alloc_l
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_r,    dbcsr_cuda_hostmem_alloc_d
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_c,    dbcsr_cuda_hostmem_alloc_z
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_i_2D, dbcsr_cuda_hostmem_alloc_l_2D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_r_2D, dbcsr_cuda_hostmem_alloc_d_2D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_c_2D, dbcsr_cuda_hostmem_alloc_z_2D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_i_3D, dbcsr_cuda_hostmem_alloc_l_3D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_r_3D, dbcsr_cuda_hostmem_alloc_d_3D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_c_3D, dbcsr_cuda_hostmem_alloc_z_3D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_i_4D, dbcsr_cuda_hostmem_alloc_l_4D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_r_4D, dbcsr_cuda_hostmem_alloc_d_4D
     MODULE PROCEDURE dbcsr_cuda_hostmem_alloc_c_4D, dbcsr_cuda_hostmem_alloc_z_4D
  END INTERFACE

  INTERFACE dbcsr_cuda_hostmem_deallocate
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_i,    dbcsr_cuda_hostmem_dealloc_l
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_r,    dbcsr_cuda_hostmem_dealloc_d
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_c,    dbcsr_cuda_hostmem_dealloc_z
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_i_2D, dbcsr_cuda_hostmem_dealloc_l_2D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_r_2D, dbcsr_cuda_hostmem_dealloc_d_2D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_c_2D, dbcsr_cuda_hostmem_dealloc_z_2D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_i_3D, dbcsr_cuda_hostmem_dealloc_l_3D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_r_3D, dbcsr_cuda_hostmem_dealloc_d_3D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_c_3D, dbcsr_cuda_hostmem_dealloc_z_3D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_i_4D, dbcsr_cuda_hostmem_dealloc_l_4D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_r_4D, dbcsr_cuda_hostmem_dealloc_d_4D
     MODULE PROCEDURE dbcsr_cuda_hostmem_dealloc_c_4D, dbcsr_cuda_hostmem_dealloc_z_4D
  END INTERFACE



#if defined (__DBCSR_CUDA)

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

#endif


CONTAINS

#include "dbcsr_cuda_hostmem_i.F"
#include "dbcsr_cuda_hostmem_l.F"
#include "dbcsr_cuda_hostmem_r.F"
#include "dbcsr_cuda_hostmem_d.F"
#include "dbcsr_cuda_hostmem_c.F"
#include "dbcsr_cuda_hostmem_z.F"


! *****************************************************************************
!> \brief Helper-routine performing allocation of host-pinned cuda memory.
! *****************************************************************************
 SUBROUTINE dbcsr_cuda_hostmem_alloc_raw (host_mem_c_ptr, n_bytes, error)
    TYPE(C_PTR), INTENT(OUT)                 :: host_mem_c_ptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_hostmem_alloc_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_alloc_cu(host_mem_c_ptr, &
               INT(n_bytes, KIND=C_SIZE_T), 0_c_int, 0_c_int)
    IF (istat /= 0) CALL mp_abort("dbcsr_cuda_hostmem_alloc_raw: Could not allocate host pinned memory")
#else
    CALL mp_abort("dbcsr_cuda_hostmem_alloc_raw: DBCSR_CUDA not compiled in.")
#endif
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dbcsr_cuda_hostmem_alloc_raw


END MODULE dbcsr_cuda_hostmem
