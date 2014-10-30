!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Higher-level operations on DBCSR matrices.
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - Created 2009-05-12
! *****************************************************************************

MODULE dbcsr_operations
  USE dbcsr_operations_low,            ONLY: &
       dbcsr_add, dbcsr_add_on_diag, dbcsr_btriu, dbcsr_copy, &
       dbcsr_copy_into_existing, dbcsr_filter, dbcsr_frobenius_norm, &
       dbcsr_gershgorin_norm, dbcsr_get_block_diag, dbcsr_get_diag, &
       dbcsr_hadamard_product, dbcsr_init_random, dbcsr_maxabs, dbcsr_norm, &
       dbcsr_scale, dbcsr_scale_by_vector, dbcsr_scale_mat, dbcsr_set, &
       dbcsr_set_diag, dbcsr_sum_replicated, dbcsr_trace, dbcsr_triu

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_operations'

  PUBLIC :: dbcsr_trace, dbcsr_add_on_diag,&
            dbcsr_set, dbcsr_scale, dbcsr_scale_mat, dbcsr_add, dbcsr_copy,&
            dbcsr_copy_into_existing,&
            dbcsr_get_diag, dbcsr_set_diag,&
            dbcsr_get_block_diag, dbcsr_hadamard_product,&
            dbcsr_filter, dbcsr_scale_by_vector,&
            dbcsr_btriu,&
            dbcsr_triu,&
            dbcsr_init_random,&
            dbcsr_sum_replicated, dbcsr_norm,&
            dbcsr_gershgorin_norm, dbcsr_maxabs, dbcsr_frobenius_norm


 CONTAINS

END MODULE dbcsr_operations
