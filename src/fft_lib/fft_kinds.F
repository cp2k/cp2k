!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \author Teodoro Laino - 01.2011
! *****************************************************************************
#if defined ( __FFTMKL )
#include "mkl_dfti.f90"
#endif
 
MODULE fft_kinds

#if defined ( __FFTMKL )
  USE MKL_DFTI
#endif

  IMPLICIT NONE
  PRIVATE
#if defined(__SGL) || defined(__FFTSGL)
  INTEGER, PARAMETER, PUBLIC :: wp = KIND(0.0)
  INTEGER, PARAMETER, PUBLIC :: dp = KIND(0.0)
  INTEGER, PARAMETER, PUBLIC :: sp = KIND(0.0)
#if defined ( __FFTMKL )
  INTEGER, PARAMETER, PUBLIC :: dft_prec = DFTi_SINGLE
#else
  INTEGER, PARAMETER, PUBLIC :: dft_prec = HUGE(0)
#endif
#else
  INTEGER, PARAMETER, PUBLIC :: wp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER, PUBLIC :: sp = KIND(0.0)
#if defined ( __FFTMKL )
  INTEGER, PARAMETER, PUBLIC :: dft_prec = DFTi_DOUBLE
#else
  INTEGER, PARAMETER, PUBLIC :: dft_prec = HUGE(0)
#endif
#endif

#if defined (__NAG)
  INTEGER*8                  :: dummy_64_bit
  INTEGER, PARAMETER, PUBLIC :: integer8_kind = KIND(dummy_64_bit)
#else
  INTEGER, PARAMETER, PUBLIC :: integer8_kind = 8
#endif

END MODULE fft_kinds
