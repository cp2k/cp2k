!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief collect pointers to a block of reals
!> \note
!>      needed to work around NEC compiler bug
!> \author fxcoudert [2006.12]
! *****************************************************************************
MODULE block_p_types
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: block_p_type
! *****************************************************************************
  TYPE block_p_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: block
  END TYPE block_p_type
END MODULE block_p_types
