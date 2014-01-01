!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>       2013.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_aux2_types
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_aux2_types'

  INTEGER, PARAMETER, PUBLIC           :: select_row_col = 1
  INTEGER, PARAMETER, PUBLIC           :: select_row = 2

  PUBLIC :: object01_type, object02_type

  TYPE object01_type
    INTEGER                                       :: domain
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE    :: mdata
    INTEGER                                       :: nbrows
    INTEGER                                       :: nbcols
    INTEGER                                       :: nrows
    INTEGER                                       :: ncols
    INTEGER, DIMENSION(:), ALLOCATABLE            :: dbcsr_row
    INTEGER, DIMENSION(:), ALLOCATABLE            :: dbcsr_col
    INTEGER, DIMENSION(:), ALLOCATABLE            :: size_brow
    INTEGER, DIMENSION(:), ALLOCATABLE            :: size_bcol
    INTEGER                                       :: nnodes
    INTEGER                                       :: groupid
  END TYPE object01_type

  TYPE object02_type
    INTEGER, DIMENSION(:), ALLOCATABLE     :: index1
    INTEGER, DIMENSION(:,:), ALLOCATABLE   :: pairs
  END TYPE object02_type

END MODULE almo_scf_aux2_types 


