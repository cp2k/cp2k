!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types to handle submatrices
!> \par History
!>       2013.01 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE domain_submatrix_types
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'domain_submatrix_types'

   INTEGER, PARAMETER, PUBLIC           :: select_row_col = 1
   INTEGER, PARAMETER, PUBLIC           :: select_row = 2

   PUBLIC :: domain_submatrix_type, domain_map_type

   ! submatrix storage with the meta-data necessary to convert
   ! the submatrix into the DBCSR format
   TYPE domain_submatrix_type
      INTEGER                                       :: domain = 0
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE    :: mdata
      INTEGER                                       :: nbrows = 0
      INTEGER                                       :: nbcols = 0
      INTEGER                                       :: nrows = 0
      INTEGER                                       :: ncols = 0
      INTEGER, DIMENSION(:), ALLOCATABLE            :: dbcsr_row
      INTEGER, DIMENSION(:), ALLOCATABLE            :: dbcsr_col
      INTEGER, DIMENSION(:), ALLOCATABLE            :: size_brow
      INTEGER, DIMENSION(:), ALLOCATABLE            :: size_bcol
      INTEGER                                       :: nnodes = 0
      TYPE(mp_comm_type)                            :: group = mp_comm_type()
   END TYPE domain_submatrix_type

   TYPE domain_map_type
      INTEGER, DIMENSION(:), ALLOCATABLE     :: index1
      INTEGER, DIMENSION(:, :), ALLOCATABLE   :: pairs
   END TYPE domain_map_type

END MODULE domain_submatrix_types

