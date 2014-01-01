!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author Teodoro Laino - 01.2011
! *****************************************************************************

MODULE fft_kinds

  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER, PUBLIC :: wp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: sp = SELECTED_REAL_KIND ( 6, 30 )

#if defined (__NAG)
  INTEGER*8                  :: dummy_64_bit
  INTEGER, PARAMETER, PUBLIC :: integer8_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER, PUBLIC :: integer8_kind = 8
#endif

END MODULE fft_kinds
