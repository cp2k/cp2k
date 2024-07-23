!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Type to store data about a (1D or 3D) FFT, including FFTW plan
!> \par History
!>      IAB 09-Jan-2009 : initial version
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!>      IAB 09-Oct-2009 : Added additional fields needed when using OpenMP
!>                        (c) The Numerical Algorithms Group (NAG) Ltd, 2009 on behalf of the HECToR project
!> \author JGH
! **************************************************************************************************

MODULE fft_plan
   USE ISO_C_BINDING,                   ONLY: C_NULL_PTR,&
                                              C_PTR

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fft_plan_type, FWFFT, BWFFT

   ! We use the signs as they are also used for FFT-SG
   INTEGER, PARAMETER :: FWFFT = +1, BWFFT = -1

   TYPE fft_plan_type
      LOGICAL                             :: trans = .FALSE., fft_in_place = .FALSE., valid = .FALSE., separated_plans = .FALSE.
      INTEGER                             :: n = -1, m = -1
      INTEGER, DIMENSION(3)               :: n_3d = -1

!   Handle for the FFTW plan
      TYPE(C_PTR)                         :: fftw_plan_fw = C_NULL_PTR
      TYPE(C_PTR)                         :: fftw_plan_bw = C_NULL_PTR

!   Plan for the remaining rows for 1D FFT when number of threads does not divide the number of rows exactly
!$    TYPE(C_PTR)                         :: alt_fftw_plan_fw = C_NULL_PTR
!$    TYPE(C_PTR)                         :: alt_fftw_plan_bw = C_NULL_PTR
!$    LOGICAL                             :: need_alt_plan = .FALSE.
!$    INTEGER                             :: num_threads_needed = -1, num_rows = -1, alt_num_rows = -1

!   Individual plans (used by hand-optimised 3D FFT)
      TYPE(C_PTR)                     :: fftw_plan_nx_fw = C_NULL_PTR, fftw_plan_ny_fw = C_NULL_PTR, fftw_plan_nz_fw = C_NULL_PTR
      TYPE(C_PTR)                     :: fftw_plan_nx_bw = C_NULL_PTR, fftw_plan_ny_bw = C_NULL_PTR, fftw_plan_nz_bw = C_NULL_PTR
!   Plans for the remaining rows (when the number of threads does not divide the number of rows exactly)
   TYPE(C_PTR)                     :: fftw_plan_nx_r_fw = C_NULL_PTR, fftw_plan_ny_r_fw = C_NULL_PTR, fftw_plan_nz_r_fw = C_NULL_PTR
   TYPE(C_PTR)                     :: fftw_plan_nx_r_bw = C_NULL_PTR, fftw_plan_ny_r_bw = C_NULL_PTR, fftw_plan_nz_r_bw = C_NULL_PTR

      TYPE(C_PTR) :: spfft_transform = C_NULL_PTR

   END TYPE fft_plan_type

END MODULE fft_plan
