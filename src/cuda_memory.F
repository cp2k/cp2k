!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  routines involved in management of device memory when cuda support
!> is enabled
!> \par History
!>      02.2008 created [bgl]
!> \author Benjamin G Levine
! *****************************************************************************
MODULE cuda_memory
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER , PRIVATE               :: sp=KIND(0.0)
  !the maximimum number of simultaneous allocations allowed
  INTEGER, PARAMETER, PRIVATE               :: max_allocs = 20
  !the offsets of the various allocations (in units of 4 bytes)
  INTEGER, PRIVATE                          :: ioffset(max_allocs)
  !the size of the various allocations (in units of 4 bytes)
  INTEGER, PRIVATE                          :: isize(max_allocs)
  !the current number of allocations
  INTEGER, PRIVATE                          :: n_allocs
  !the total length of the device memory (in units of 4 bytes)
  INTEGER, PRIVATE                          :: length
  !the pointer to the begining of the device memory array
  REAL(kind=sp), POINTER, PRIVATE           :: device_memory
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cuda_memory'
 
  PUBLIC  :: cuda_device_mem_init, cuda_device_mem_release
CONTAINS

! *****************************************************************************
!> \brief allocates GPU device memory to be used later.  Device memory is
!> allocated up front to eliminate the cost of allocating repeatedly later.
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_device_mem_init(root_section)
    TYPE(section_vals_type), POINTER         :: root_section

    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_device_mem_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: memory
    TYPE(cp_error_type)                      :: error
    TYPE(section_vals_type), POINTER         :: cuda_section

    cuda_section => section_vals_get_subs_vals(root_section,&
      "GLOBAL%CUDA",error=error)
    CALL section_vals_val_get(cuda_section,"MEMORY",i_val=memory,error=error)

#if defined ( __CUDA ) || defined ( __FFTCU )

    CALL cuda_device_mem_init_cu(memory)

#endif

  END SUBROUTINE cuda_device_mem_init

! *****************************************************************************
!> \free all CUDA memory, which was previously allocated by cuda_device_mem_init
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_device_mem_release()

#if defined ( __CUDA ) || defined ( __FFTCU )

    CALL cuda_device_mem_release_cu(device_memory)

#endif 

  END SUBROUTINE cuda_device_mem_release

END MODULE cuda_memory
