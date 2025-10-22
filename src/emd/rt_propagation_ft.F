!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Separation of Fourier transform utilities into separate file
!> \author Stepan Marek (08.24)
! **************************************************************************************************
MODULE rt_propagation_ft
   USE kinds,           ONLY: dp
   USE mathconstants,   ONLY: twopi
   USE fft_plan,        ONLY: fft_plan_type
   USE fft_lib,         ONLY: fft_library, &
                              fft_alloc, &
                              fft_dealloc, &
                              fft_create_plan_1dm, &
                              fft_destroy_plan, &
                              fft_1dm

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse_io"

   PUBLIC :: multi_fft, &
             fft_freqs, &
             fft_shift

CONTAINS
! **************************************************************************************************
!> \brief Naively calculates the Fourier transform - it is not the bottleneck of this calculation
!> \note Uses physics ordering in frequencies, those can be constructed by fft_freq
!> \param time_series Timestamps in atomic units of time
!> \param value_series Values to be Fourier transformed - moments, field etc.
!> \param result_series FT of the value series
!> \param damping Applied exponential damping
!> \param subtract_value Value to be subtracted from the value_series (for example initial value)
!> \author Stepan Marek
!> \par History
!>      10.2025 Refactored for use with multi_fft routine, moved to separate file [Stepan Marek]
!>      09.2024 Initial version [Stepan Marek]
! **************************************************************************************************
   SUBROUTINE ft_simple(time_series, value_series, result_series, damping, subtract_value)
      REAL(kind=dp), DIMENSION(:)                        :: time_series
      COMPLEX(kind=dp), DIMENSION(:)                     :: value_series
      COMPLEX(kind=dp), DIMENSION(:)                     :: result_series
      REAL(kind=dp)                                      :: damping
      INTEGER                                            :: i0
      COMPLEX(kind=dp)                                   :: subtract_value
      CHARACTER(len=*), PARAMETER                        :: routineN = "ft_simple"
      INTEGER                                            :: N, i, j, start, handle
      REAL(kind=dp)                                      :: t0, delta_t, delta_omega
      LOGICAL                                            :: subtract_initial

      CALL timeset(routineN, handle)

      N = SIZE(time_series)

      delta_t = time_series(2) - time_series(1)

      IF (MOD(N,2) == 0) THEN
         start = -N/2
      ELSE
         start = -(N-1)/2
      END IF

      ! TODO : At least OMP, but ideally even MPI parallelize, or handle this on higher level?
      DO i = 1, N
         result_series(i) = CMPLX(0.0, 0.0, kind=dp)
         DO j = 1, N
            result_series(i) = result_series(i) + EXP(CMPLX(0.0, twopi*(start + i-1)*(j-1)/N)) * &
                    EXP(- damping * delta_t * (j-1)) * (value_series(j) - subtract_value)
         END DO
      END DO
      result_series(:) = delta_t*result_series(:)

      CALL timestop(handle)

   END SUBROUTINE
! **************************************************************************************************
!> \brief Calculates the Fourier transform - couples to FFT libraries in CP2K, if available
!> \param time_series Timestamps in atomic units of time
!> \param value_series Values to be Fourier transformed - moments, field etc. Real only. Many series can be provided.
!> \param result_series FT of the value series - complex numbers
!> \param damping_opt Supply custom exponential damping - default is 4.0/totalTime, i.e. ratio
!>                    of last and first element in windowed value series is reduced by e^(-4)
!> \param t0_opt Carry the FT only starting from certain time - allows for exclusion of trace before
!>               the pulse application etc.
!> \param subtract_initial Subtract the value at t0 from the values to be transformed. Default is true
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE multi_fft(time_series, value_series, result_series, omega_series, &
                   damping_opt, t0_opt, subtract_initial_opt)
      REAL(kind=dp), DIMENSION(:), POINTER              :: time_series
      ! For value and result series: Index 1 - different series, Index 2 - single series entry
      COMPLEX(kind=dp), DIMENSION(:,:), POINTER         :: value_series
      COMPLEX(kind=dp), DIMENSION(:,:), POINTER         :: result_series
      REAL(kind=dp), DIMENSION(:), POINTER              :: omega_series
      REAL(kind=dp), OPTIONAL                           :: damping_opt
      REAL(kind=dp), OPTIONAL                           :: t0_opt
      LOGICAL, OPTIONAL                                 :: subtract_initial_opt
      CHARACTER(len=*), PARAMETER                       :: routineN="multi_fft"
      INTEGER                                           :: handle, &
                                                           nsamples, &
                                                           nseries, &
                                                           i, &
                                                           i0, &
                                                           j, &
                                                           stat
      REAL(kind=dp)                                     :: damping, &
                                                           t0, &
                                                           t_total
      LOGICAL                                           :: subtract_initial, &
                                                           ionode
      TYPE(fft_plan_type)                               :: fft_plan
      COMPLEX(kind=dp), CONTIGUOUS, DIMENSION(:), &
              POINTER                                   :: samples, &
                                                           ft_samples
      ! Evaluate optional arguments
      ! Start with t0
      t0 = 0.0_dp
      IF (PRESENT(t0_opt)) t0 = t0_opt
      ! Determine zero index
      i0 = 1
      DO i=1,SIZE(time_series)
         IF (time_series(i) >= t0) THEN
            i0 = i
            EXIT
         END IF
      END DO
      ! Determine nsamples
      nsamples = SIZE(time_series) - i0 + 1
      ! Determine total time
      t_total = time_series(SIZE(time_series)) - time_series(i0)
      ! Now can determine default damping
      damping = 4.0_dp/(t_total)
      ! Damping option supplied in au units of time
      IF (PRESENT(damping_opt)) THEN
         IF (damping_opt > 0.0_dp) THEN
            damping = 1.0_dp/damping_opt
         ELSE IF (damping_opt == 0.0_dp) THEN
            ! Special case - zero damping
            damping = 0.0_dp
         END IF
      END IF
      ! subtract initial
      subtract_initial = .TRUE.
      IF (PRESENT(subtract_initial_opt)) subtract_initial = subtract_initial_opt

      ! Determine nseries
      nseries = SIZE(value_series, 1)
      ! Reallocate results if nsamples lower than current size
      IF (nsamples /= SIZE(result_series, 2)) THEN
         DEALLOCATE (result_series)
         ALLOCATE (result_series(nseries, nsamples), source=CMPLX(0.0, 0.0, kind=dp))
      ENDIF

      ! Calculate the omega series values, ordered from negative to positive
      CALL fft_freqs(nsamples, t_total, omega_series, fft_ordering_opt=.FALSE.)

#if defined(__FFTW3)
      ! Use FFTW3 library
      ! Allocate the in-out arrays (on every rank)
      CALL timeset(routineN, handle)
      NULLIFY (samples)
      NULLIFY (ft_samples)
      CALL fft_alloc(samples, [nsamples * nseries])
      CALL fft_alloc(ft_samples, [nsamples * nseries])
      CALL fft_create_plan_1dm(fft_plan, fft_library("FFTW3"), -1, .FALSE., nsamples, nseries, samples, ft_samples, 3)
      ! TODO : Fill the samples with actual data
      DO i = 1, nseries
         DO j = 1, nsamples
            ! Subtract initial value if required
            IF (subtract_initial) THEN
               samples(j + (i-1)*nsamples) = value_series(i, i0 + j - 1) - value_series(i, i0)
            ELSE
               samples(j + (i-1)*nsamples) = value_series(i, i0 + j - 1)
            ENDIF
            ! Apply damping
            samples(j + (i-1)*nsamples) = samples(j + (i-1)*nsamples) * &
                    EXP( - damping * (time_series(i0 + j - 1) - time_series(i0)))
         END DO
      END DO
      ! TODO : Carry out the transform
      ! Scale by dt - to transform to an integral
      CALL fft_1dm(fft_plan, samples, ft_samples, time_series(2) - time_series(1), stat)
      IF (stat /= 0) CPABORT("Failed FFT in RTBSE output routines.")
      DO i=1, nseries
         ! TODO : Reorder/fftshift
         CALL fft_shift(ft_samples((i-1)*nsamples+1:i*nsamples))
         result_series(i, :) = ft_samples((i-1)*nsamples+1:i*nsamples)
      END DO
      ! Deallocate
      CALL fft_dealloc(samples)
      CALL fft_dealloc(ft_samples)
      CALL fft_destroy_plan(fft_plan)
      CALL timestop(handle)
#else
      ! FFTSG can only use a specific length - use custom (slow) FT library instead
      ! For each series, do ft_simple call
      DO i=1,nseries
         CALL ft_simple(time_series(i0:SIZE(time_series)), &
            value_series(i,i0:SIZE(value_series, 2)), result_series(i,1:nsamples), &
            damping, value_series(i, i0))
      END DO
#endif
   END SUBROUTINE multi_fft
! **************************************************************************************************
!> \brief Switches the order in result of FT, so that negative frequencies go first
!> \param source Array containing the FT - buffer is used to reorder it
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE fft_shift(source)
      COMPLEX(kind=dp), DIMENSION(:)                    :: source
      COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE       :: buffer
      ! Boundary indices for positive/negative part of the spectrum
      ! Index 1 : 1 = transformed order, 2 = FFT order
      INTEGER, DIMENSION(2)                             :: pos_lower, &
                                                           pos_upper, &
                                                           neg_lower, &
                                                           neg_upper
      INTEGER                                           :: n

      n = SIZE(source)
      IF (MOD(n, 2) == 0) THEN
         ! Even case
         pos_lower(1) = n/2 + 1
         pos_upper(2) = n/2
         neg_lower(2) = n/2 + 1
         neg_upper(1) = n/2
      ELSE
         pos_lower(1) = (n+1)/2
         pos_upper(2) = (n+1)/2
         neg_lower(2) = (n+1)/2 + 1
         neg_upper(1) = (n-1)/2
      END IF
      ! Parity independent positions
      pos_lower(2) = 1
      pos_upper(1) = n
      neg_lower(1) = 1
      neg_upper(2) = n

      ALLOCATE (buffer(n))
      buffer(neg_lower(1):neg_upper(1)) = source(neg_lower(2):neg_upper(2))
      buffer(pos_lower(1):pos_upper(1)) = source(pos_lower(2):pos_upper(2))
      source(:) = buffer(:)
      DEALLOCATE (buffer)

   END SUBROUTINE fft_shift
! **************************************************************************************************
!> \brief Switches the order in result of FT, so that negative frequencies go first
!> \param source Array containing the FT - buffer is used to reorder it
!> \author Stepan Marek
!> \date 10.2025
! **************************************************************************************************
   SUBROUTINE fft_freqs(n, t_total, omegas, fft_ordering_opt)
      ! Number of FT samples
      INTEGER                                           :: n
      ! Total window time, dt = nsamples / t_total
      REAL(kind=dp)                                     :: t_total
      REAL(kind=dp), DIMENSION(:)                       :: omegas
      LOGICAL, OPTIONAL                                 :: fft_ordering_opt
      LOGICAL                                           :: fft_ordering
      INTEGER                                           :: i, &
                                                           start, &
                                                           finish

      ! Determine the order, by default, use physics order,
      ! i.e. negative frequencies before positive ones
      fft_ordering = .FALSE.
      IF (PRESENT(fft_ordering_opt)) fft_ordering = fft_ordering_opt

      IF (.NOT. fft_ordering) THEN
         ! Physics order case
         ! Unit frequencies at
         !  - for even n : -n/2, -n/2 + 1, -n/2 + 2, ..., -1, 0, 1, ..., n/2 - 1
         !  - for odd n : -(n-1)/2, -(n-1)/2 + 1, ..., -1, 0, 1, ..., (n-1)/2
         IF (MOD(n,2) == 0) THEN
            start = -n/2
         ELSE
            start = -(n-1)/2
         END IF
         DO i=1, n
            omegas(i) = start + i-1
         END DO
      ELSE
         ! FFT order case
         ! Unit frequencies at
         !  - for even n : 0, 1, ..., n/2 - 1, -n/2, -n/2 + 1, -n/2 + 2, ..., -1
         !  - for odd n : 0, 1, ..., (n-1)/2, -(n-1)/2, -(n-1)/2 + 1, ..., -1
         IF (MOD(n,2) == 0) THEN
            finish = n/2 - 1
            start = -n/2
         ELSE
            finish = (n-1)/2
            start = -(n-1)/2
         ENDIF
         ! Positive frequencies
         DO i=1,finish+1
            omegas(i) = (i-1)
         END DO
         ! Negative frequencies
         DO i=finish+2,n
            omegas(i) = start + i - finish - 2
         END DO
      ENDIF

      ! Finally, multiply by the factor to translate to angular frequency
      omegas(:) = omegas(:) * twopi / t_total
   END SUBROUTINE fft_freqs

END MODULE rt_propagation_ft
