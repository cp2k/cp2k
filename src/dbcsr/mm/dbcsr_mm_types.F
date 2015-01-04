!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Datatypes used by the dbcsr matrix-matrix multiplication machinery.
!>	
!> \author  Urban Borstnik
!>
!> <b>Modification history:</b>
!>  - 2013-01    reorganized code (Ole Schuett)
! *****************************************************************************
MODULE dbcsr_mm_types

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_types'


  PUBLIC :: stack_descriptor_type

 ! Definitions for the members of a parameter stack.
  INTEGER, PARAMETER, PUBLIC :: dbcsr_ps_width = 7
  INTEGER, PARAMETER, PUBLIC :: p_m = 1
  INTEGER, PARAMETER, PUBLIC :: p_n = 2
  INTEGER, PARAMETER, PUBLIC :: p_k = 3
  INTEGER, PARAMETER, PUBLIC :: p_a_first = 4
  INTEGER, PARAMETER, PUBLIC :: p_b_first = 5
  INTEGER, PARAMETER, PUBLIC :: p_c_first = 6
  INTEGER, PARAMETER, PUBLIC :: p_c_blk = 7


  TYPE stack_descriptor_type
     INTEGER :: m, n, k, max_m ,max_n, max_k
     LOGICAL :: defined_mnk
  END TYPE stack_descriptor_type


END MODULE dbcsr_mm_types
