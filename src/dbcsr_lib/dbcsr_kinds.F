!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author vw
! *****************************************************************************
MODULE dbcsr_kinds

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: sp, dp
  PUBLIC :: int_1, int_2, int_4, int_8
  PUBLIC :: real_4, real_8
  PUBLIC :: int_1_size, int_2_size, int_4_size, int_8_size
  PUBLIC :: real_4_size, real_8_size
  PUBLIC :: default_string_length

  INTEGER, PARAMETER :: real_4 = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: real_8 = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER :: real_4_size = 4
  INTEGER, PARAMETER :: real_8_size = 8

  INTEGER, PARAMETER :: int_1  = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: int_1_size = BIT_SIZE(INT(0, int_1))/8

  INTEGER, PARAMETER :: int_2  = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: int_2_size = BIT_SIZE(INT(0, int_2))/8

  INTEGER, PARAMETER :: int_4 = SELECTED_INT_KIND(5)
  INTEGER, PARAMETER :: int_4_size = BIT_SIZE(INT(0, int_4))/8

  INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)
  INTEGER, PARAMETER :: int_8_size = BIT_SIZE(INT(0, int_8))/8

  INTEGER, PARAMETER :: sp = real_4
  INTEGER, PARAMETER :: dp = real_8

  INTEGER, PARAMETER :: default_string_length=80

END MODULE dbcsr_kinds

