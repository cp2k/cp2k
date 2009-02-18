!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type to store data about a (1D or 3D) FFT, including FFTW plan
!> \par History
!>      IAB 09-Jan-2009 : initial version
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************

MODULE fft_plan

  IMPLICIT NONE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER :: dp = KIND(0.0)
#else
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
#endif

#if defined (__NAG)
  INTEGER*8          :: dummy_64_bit
  INTEGER, PARAMETER :: integer8_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER :: integer8_kind = 8
#endif

  PUBLIC :: fft_plan_type, dp

  TYPE fft_plan_type

    INTEGER                             :: fft_type
    INTEGER                             :: fsign
    LOGICAL                             :: trans, fft_in_place, valid
    INTEGER                             :: n, m
    INTEGER, DIMENSION(3)               :: n_3d
    INTEGER (KIND=integer8_kind)        :: fftw_plan

  END TYPE fft_plan_type


END MODULE fft_plan
