!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cuda_memory'

  PUBLIC  :: cuda_device_mem_init,&
             cuda_device_mem_release

  INTERFACE
    SUBROUTINE cuda_device_mem_init_cu (memory) &
             BIND(C, name="cuda_device_mem_init_cu_")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT), INTENT(IN)          :: memory

    END SUBROUTINE cuda_device_mem_init_cu
  END INTERFACE

  INTERFACE
    SUBROUTINE cuda_device_mem_release_cu () &
             BIND(C, name="cuda_device_mem_release_cu_")
      USE ISO_C_BINDING
    END SUBROUTINE cuda_device_mem_release_cu
  END INTERFACE

  INTERFACE
    SUBROUTINE cuda_device_fftcu_release_cu () &
             BIND(C, name="fftcu_release_")
      USE ISO_C_BINDING
    END SUBROUTINE cuda_device_fftcu_release_cu
  END INTERFACE
CONTAINS

! *****************************************************************************
!> \brief allocates GPU device memory to be used later.  Device memory is
!> allocated up front to eliminate the cost of allocating repeatedly later.
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_device_mem_init(root_section, error)
    TYPE(section_vals_type), POINTER, INTENT(IN) :: root_section
    TYPE(cp_error_type), INTENT(INOUT)           :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_device_mem_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                      :: memory
    TYPE(section_vals_type), POINTER             :: cuda_section

!    cuda_section => section_vals_get_subs_vals(root_section,&
!      "GLOBAL%CUDA",error=error)
!    CALL section_vals_val_get(cuda_section, "MEMORY", i_val=memory, error=error)

#if defined ( __PW_CUDA ) || defined ( __CUBLASDP )

    memory = 0
    CALL cuda_device_mem_init_cu(memory)
    CALL cublasinit_cu()

#endif

  END SUBROUTINE cuda_device_mem_init

! *****************************************************************************
!> \free all CUDA memory, which was previously allocated by cuda_device_mem_init
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_device_mem_release()
    CHARACTER(len=*), PARAMETER :: routineN = 'cuda_device_mem_release', &
      routineP = moduleN//':'//routineN


#if defined ( __PW_CUDA ) || defined ( __CUBLASDP )

    CALL cuda_device_fftcu_release_cu()
    CALL cuda_device_mem_release_cu()
    CALL cublasshutdown_cu()

#endif

  END SUBROUTINE cuda_device_mem_release

END MODULE cuda_memory
