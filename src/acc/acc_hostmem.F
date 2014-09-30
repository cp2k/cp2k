!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014 CP2K developers group                           !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! *****************************************************************************
MODULE acc_hostmem
#if defined (__ACC)
  USE ISO_C_BINDING
#endif
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
                                             acc_stream_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_hostmem'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.


  PUBLIC :: acc_hostmem_allocate, acc_hostmem_deallocate


  INTERFACE acc_hostmem_allocate
     MODULE PROCEDURE acc_hostmem_alloc_i,    acc_hostmem_alloc_l
     MODULE PROCEDURE acc_hostmem_alloc_r,    acc_hostmem_alloc_d
     MODULE PROCEDURE acc_hostmem_alloc_c,    acc_hostmem_alloc_z
     MODULE PROCEDURE acc_hostmem_alloc_i_2D, acc_hostmem_alloc_l_2D
     MODULE PROCEDURE acc_hostmem_alloc_r_2D, acc_hostmem_alloc_d_2D
     MODULE PROCEDURE acc_hostmem_alloc_c_2D, acc_hostmem_alloc_z_2D
     MODULE PROCEDURE acc_hostmem_alloc_i_3D, acc_hostmem_alloc_l_3D
     MODULE PROCEDURE acc_hostmem_alloc_r_3D, acc_hostmem_alloc_d_3D
     MODULE PROCEDURE acc_hostmem_alloc_c_3D, acc_hostmem_alloc_z_3D
     MODULE PROCEDURE acc_hostmem_alloc_i_4D, acc_hostmem_alloc_l_4D
     MODULE PROCEDURE acc_hostmem_alloc_r_4D, acc_hostmem_alloc_d_4D
     MODULE PROCEDURE acc_hostmem_alloc_c_4D, acc_hostmem_alloc_z_4D
  END INTERFACE

  INTERFACE acc_hostmem_deallocate
     MODULE PROCEDURE acc_hostmem_dealloc_i,    acc_hostmem_dealloc_l
     MODULE PROCEDURE acc_hostmem_dealloc_r,    acc_hostmem_dealloc_d
     MODULE PROCEDURE acc_hostmem_dealloc_c,    acc_hostmem_dealloc_z
     MODULE PROCEDURE acc_hostmem_dealloc_i_2D, acc_hostmem_dealloc_l_2D
     MODULE PROCEDURE acc_hostmem_dealloc_r_2D, acc_hostmem_dealloc_d_2D
     MODULE PROCEDURE acc_hostmem_dealloc_c_2D, acc_hostmem_dealloc_z_2D
     MODULE PROCEDURE acc_hostmem_dealloc_i_3D, acc_hostmem_dealloc_l_3D
     MODULE PROCEDURE acc_hostmem_dealloc_r_3D, acc_hostmem_dealloc_d_3D
     MODULE PROCEDURE acc_hostmem_dealloc_c_3D, acc_hostmem_dealloc_z_3D
     MODULE PROCEDURE acc_hostmem_dealloc_i_4D, acc_hostmem_dealloc_l_4D
     MODULE PROCEDURE acc_hostmem_dealloc_r_4D, acc_hostmem_dealloc_d_4D
     MODULE PROCEDURE acc_hostmem_dealloc_c_4D, acc_hostmem_dealloc_z_4D
  END INTERFACE



#if defined (__ACC)

  INTERFACE
     FUNCTION cuda_host_mem_alloc_cu(mem, n, stream_ptr) RESULT (istat) BIND(C, name="acc_host_mem_allocate")
       USE ISO_C_BINDING
    TYPE(C_PTR)                              :: mem
    INTEGER(KIND=C_SIZE_T), INTENT(IN), &
      VALUE                                  :: n
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_host_mem_alloc_cu
  END INTERFACE

  INTERFACE
     FUNCTION cuda_host_mem_dealloc_cu(mem, stream_ptr) RESULT (istat) bind(C, name="acc_host_mem_deallocate")
       USE ISO_C_BINDING
    TYPE(C_PTR), VALUE                       :: mem, stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION cuda_host_mem_dealloc_cu
  END INTERFACE

#endif


CONTAINS

#include "acc_hostmem_i.f90"
#include "acc_hostmem_l.f90"
#include "acc_hostmem_r.f90"
#include "acc_hostmem_d.f90"
#include "acc_hostmem_c.f90"
#include "acc_hostmem_z.f90"


! *****************************************************************************
!> \brief Helper-routine performing allocation of host-pinned cuda memory.
!> \param host_mem_c_ptr pointer to allocated memory
!> \param n_bytes number of bytes to allocate
!> \param stream ...
! *****************************************************************************
#if defined (__ACC)
 SUBROUTINE acc_hostmem_alloc_raw (host_mem_c_ptr, n_bytes, stream)
    TYPE(C_PTR), INTENT(OUT)                 :: host_mem_c_ptr
    INTEGER, INTENT(IN)                      :: n_bytes
    TYPE(acc_stream_type), INTENT(IN)        :: stream

    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_alloc_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    IF(.NOT. acc_stream_associated(stream)) &
         STOP "acc_hostmem_alloc_raw: stream not associated"

    stream_cptr = acc_stream_cptr(stream)

    istat = cuda_host_mem_alloc_cu(host_mem_c_ptr, &
               INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
    IF (istat /= 0) STOP "acc_hostmem_alloc_raw: Could not allocate host pinned memory"
  END SUBROUTINE acc_hostmem_alloc_raw
#endif

#if defined (__ACC)
! *****************************************************************************
!> \brief ...
!> \param host_mem_c_ptr ...
!> \param stream ...
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_raw (host_mem_c_ptr, stream)
    TYPE(C_PTR), INTENT(IN)                  :: host_mem_c_ptr
    TYPE(acc_stream_type), INTENT(IN)        :: stream

    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_raw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    TYPE(C_PTR)                              :: stream_cptr

    IF(.NOT. acc_stream_associated(stream)) &
         STOP "acc_hostmem_dealloc_raw: stream not associated"

    stream_cptr = acc_stream_cptr(stream)

    istat = cuda_host_mem_dealloc_cu(host_mem_c_ptr, stream_cptr)
    IF (istat /= 0 ) STOP "acc_hostmem_dealloc_raw: Could not deallocate host pinned memory"
  END SUBROUTINE acc_hostmem_dealloc_raw
#endif

END MODULE acc_hostmem
