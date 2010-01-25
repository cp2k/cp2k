!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author vw
! *****************************************************************************
MODULE dbcsr_kinds

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: sp, dp
  PUBLIC :: int_4, int_8
  PUBLIC :: real_4, real_8
  PUBLIC :: default_string_length

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )

  INTEGER, PARAMETER :: real_4 = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: real_8 = SELECTED_REAL_KIND ( 14, 200 )

  INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)
  INTEGER, PARAMETER :: int_4 = SELECTED_INT_KIND(5)

  INTEGER, PARAMETER :: default_string_length=80

END MODULE dbcsr_kinds

