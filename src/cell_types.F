!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles all functions related to the CELL
!> \par History
!>      11.2008 Teodoro Laino [tlaino] - deeply cleaning cell_type from units
!> \author Matthias KracK (16.01.2002, based on a earlier version of CJM, JGH)
! *****************************************************************************
MODULE cell_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cell_types'

  INTEGER, SAVE, PRIVATE :: last_cell_id = 0

! *****************************************************************************
!> \brief   Type defining parameters related to the simulation cell
!> \version 1.0
! *****************************************************************************
  TYPE cell_type
     INTEGER                           :: ref_count,id_nr
     INTEGER                           :: symmetry_id
     LOGICAL                           :: orthorhombic ! actually means a diagonal hmat
     REAL(KIND = dp)                   :: deth
     INTEGER, DIMENSION(3)             :: perd
     REAL(KIND = dp), DIMENSION(3,3)   :: hmat,h_inv
  END TYPE cell_type

  TYPE cell_p_type
     TYPE(cell_type),POINTER :: cell
  END TYPE cell_p_type

  ! Public data types
  PUBLIC :: cell_type, cell_p_type

END MODULE cell_types
