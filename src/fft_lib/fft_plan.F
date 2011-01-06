!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Type to store data about a (1D or 3D) FFT, including FFTW plan
!> \par History
!>      IAB 09-Jan-2009 : initial version
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!>      IAB 09-Oct-2009 : Added additional fields needed when using OpenMP
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! *****************************************************************************

MODULE fft_plan
  USE fft_kinds,                       ONLY: integer8_kind

  IMPLICIT NONE

  PUBLIC :: fft_plan_type

  TYPE fft_plan_type
    INTEGER                             :: fft_type
    INTEGER                             :: fsign
    LOGICAL                             :: trans, fft_in_place, valid
    INTEGER                             :: n, m
    INTEGER, DIMENSION(3)               :: n_3d
    INTEGER (KIND=integer8_kind)        :: fftw_plan
!   Needed for OpenMP case where not all threads have the same number of rows
!$  INTEGER (KIND=integer8_kind)        :: alt_fftw_plan
!$  LOGICAL                             :: need_alt_plan
!$  INTEGER                             :: num_threads_needed, num_rows, alt_num_rows
  END TYPE fft_plan_type

END MODULE fft_plan
