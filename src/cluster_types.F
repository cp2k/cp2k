!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Cluster methods for atoms
!> \par History
!> \author VW (2010)
! *****************************************************************************
MODULE cluster_types

  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cluster_types'

! *** Public subroutines ***

  PUBLIC :: cluster_type, box_type

  TYPE boxes_type
     INTEGER, DIMENSION(3) :: nboxes = 0
     INTEGER, DIMENSION(:,:,:), POINTER :: boxes_ptr => NULL()
     INTEGER, DIMENSION(:), POINTER :: boxes_data => NULL()
  END TYPE boxes_type

  TYPE box_type
     INTEGER :: nparticles = 0
     INTEGER, DIMENSION(:), POINTER :: atom_list => NULL()
     REAL(dp), DIMENSION(:,:), POINTER :: xyz => NULL()
  END TYPE box_type

  TYPE cluster_type
     INTEGER                        :: nclusters=0, nparticles=0
     INTEGER, DIMENSION(:), POINTER :: cluster_ptr => NULL(), &
                                       particle_idx => NULL()
  END TYPE cluster_type

END MODULE cluster_types
