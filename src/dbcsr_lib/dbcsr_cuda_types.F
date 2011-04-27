!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group   !
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
MODULE dbcsr_cuda_types
#if !defined (__HAS_NO_ISO_C_BINDING)
  USE ISO_C_BINDING
#endif
  USE dbcsr_error_handling
  USE dummy_c_bindings

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_types'

  LOGICAL, PARAMETER :: careful_mod = .TRUE.

  PUBLIC :: dbcsr_cuda_mem_type

  PUBLIC :: dbcsr_cuda_mem_type_i4, dbcsr_cuda_mem_type_i8,&
            dbcsr_cuda_mem_type_r4, dbcsr_cuda_mem_type_r8,&
            dbcsr_cuda_mem_type_c4, dbcsr_cuda_mem_type_c8

  TYPE dbcsr_cuda_mem_type_i4
     TYPE(C_PTR) :: ref
     INTEGER :: len, req_len
  END TYPE dbcsr_cuda_mem_type_i4

  TYPE dbcsr_cuda_mem_type_i8
     TYPE(C_PTR) :: ref
     INTEGER :: len, req_len
  END TYPE dbcsr_cuda_mem_type_i8

  TYPE dbcsr_cuda_mem_type_r4
     TYPE(C_PTR) :: ref
     INTEGER :: len, req_len
  END TYPE dbcsr_cuda_mem_type_r4

  TYPE dbcsr_cuda_mem_type_r8
     TYPE(C_PTR) :: ref
     INTEGER :: len, req_len
  END TYPE dbcsr_cuda_mem_type_r8

  TYPE dbcsr_cuda_mem_type_c4
     TYPE(C_PTR) :: ref
     INTEGER :: len, req_len
  END TYPE dbcsr_cuda_mem_type_c4

  TYPE dbcsr_cuda_mem_type_c8
     TYPE(C_PTR) :: ref
     INTEGER :: len, req_len
  END TYPE dbcsr_cuda_mem_type_c8


  TYPE dbcsr_cuda_mem_type
     INTEGER :: data_type
     TYPE(dbcsr_cuda_mem_type_i4) :: d_i
     TYPE(dbcsr_cuda_mem_type_i8) :: d_l
     TYPE(dbcsr_cuda_mem_type_r4) :: d_r
     TYPE(dbcsr_cuda_mem_type_r8) :: d_d
     TYPE(dbcsr_cuda_mem_type_c4) :: d_c
     TYPE(dbcsr_cuda_mem_type_c8) :: d_z
  END TYPE dbcsr_cuda_mem_type



END MODULE dbcsr_cuda_types
