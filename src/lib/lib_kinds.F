!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author Teodoro Laino - 01.2011
! *****************************************************************************
MODULE lib_kinds
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: wp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: sp = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER, PUBLIC :: int_size = BIT_SIZE(0)/8

END MODULE lib_kinds

