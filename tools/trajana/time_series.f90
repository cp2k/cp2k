!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_time_series
   USE trajana_fft,                     ONLY: fft_in_place,&
                                              next_power_of_two
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: real_autocorrelation_sum, complex_autocorrelation_sum
   PUBLIC :: real_power_sum, complex_power_sum, remove_real_means, remove_complex_means

CONTAINS

   SUBROUTINE remove_real_means(series)
      REAL(dp), INTENT(INOUT)                            :: series(:, :)

      INTEGER                                            :: item

      DO item = 1, SIZE(series, 1)
         series(item, :) = series(item, :) - SUM(series(item, :))/REAL(SIZE(series, 2), dp)
      END DO
   END SUBROUTINE remove_real_means

   SUBROUTINE remove_complex_means(series)
      COMPLEX(dp), INTENT(INOUT)                         :: series(:, :)

      INTEGER                                            :: item

      DO item = 1, SIZE(series, 1)
         series(item, :) = series(item, :) - SUM(series(item, :))/REAL(SIZE(series, 2), dp)
      END DO
   END SUBROUTINE remove_complex_means

   SUBROUTINE real_autocorrelation_sum(series, maximum_lag, correlation)
      REAL(dp), INTENT(IN)                               :: series(:, :)
      INTEGER, INTENT(IN)                                :: maximum_lag
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: correlation(:)

      COMPLEX(dp), ALLOCATABLE                           :: work(:)
      INTEGER                                            :: frames, item, lag, nfft

      frames = SIZE(series, 2)
      nfft = next_power_of_two(2*frames)
      ALLOCATE (work(nfft), correlation(0:maximum_lag))
      correlation = 0.0_dp
      DO item = 1, SIZE(series, 1)
         work = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         work(:frames) = CMPLX(series(item, :), 0.0_dp, KIND=dp)
         CALL fft_in_place(work, inverse=.FALSE.)
         work = work*CONJG(work)
         CALL fft_in_place(work, inverse=.TRUE.)
         DO lag = 0, maximum_lag
            correlation(lag) = correlation(lag) + REAL(work(lag + 1), dp)
         END DO
      END DO
   END SUBROUTINE real_autocorrelation_sum

   SUBROUTINE complex_autocorrelation_sum(series, maximum_lag, correlation)
      COMPLEX(dp), INTENT(IN)                            :: series(:, :)
      INTEGER, INTENT(IN)                                :: maximum_lag
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: correlation(:)

      COMPLEX(dp), ALLOCATABLE                           :: work(:)
      INTEGER                                            :: frames, item, lag, nfft

      frames = SIZE(series, 2)
      nfft = next_power_of_two(2*frames)
      ALLOCATE (work(nfft), correlation(0:maximum_lag))
      correlation = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO item = 1, SIZE(series, 1)
         work = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         work(:frames) = series(item, :)
         CALL fft_in_place(work, inverse=.FALSE.)
         work = work*CONJG(work)
         CALL fft_in_place(work, inverse=.TRUE.)
         DO lag = 0, maximum_lag
            correlation(lag) = correlation(lag) + work(lag + 1)
         END DO
      END DO
   END SUBROUTINE complex_autocorrelation_sum

   SUBROUTINE real_power_sum(series, window_name, power, nfft, window_norm)
      REAL(dp), INTENT(IN)                               :: series(:, :)
      CHARACTER(LEN=*), INTENT(IN)                       :: window_name
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: power(:)
      INTEGER, INTENT(OUT)                               :: nfft
      REAL(dp), INTENT(OUT)                              :: window_norm

      COMPLEX(dp), ALLOCATABLE                           :: work(:)
      INTEGER                                            :: frames, frequency, item
      REAL(dp), ALLOCATABLE                              :: window(:)

      frames = SIZE(series, 2)
      nfft = next_power_of_two(frames)
      ALLOCATE (work(nfft), window(frames), power(0:nfft/2))
      CALL make_window(window_name, window)
      window_norm = SUM(window**2)
      power = 0.0_dp
      DO item = 1, SIZE(series, 1)
         work = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         work(:frames) = CMPLX(series(item, :)*window, 0.0_dp, KIND=dp)
         CALL fft_in_place(work, inverse=.FALSE.)
         DO frequency = 0, nfft/2
            power(frequency) = power(frequency) + ABS(work(frequency + 1))**2
         END DO
      END DO
   END SUBROUTINE real_power_sum

   SUBROUTINE complex_power_sum(series, window_name, power, nfft, window_norm)
      COMPLEX(dp), INTENT(IN)                            :: series(:, :)
      CHARACTER(LEN=*), INTENT(IN)                       :: window_name
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: power(:)
      INTEGER, INTENT(OUT)                               :: nfft
      REAL(dp), INTENT(OUT)                              :: window_norm

      COMPLEX(dp), ALLOCATABLE                           :: work(:)
      INTEGER                                            :: frames, frequency, item
      REAL(dp), ALLOCATABLE                              :: window(:)

      frames = SIZE(series, 2)
      nfft = next_power_of_two(frames)
      ALLOCATE (work(nfft), window(frames), power(0:nfft/2))
      CALL make_window(window_name, window)
      window_norm = SUM(window**2)
      power = 0.0_dp
      DO item = 1, SIZE(series, 1)
         work = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         work(:frames) = series(item, :)*window
         CALL fft_in_place(work, inverse=.FALSE.)
         DO frequency = 0, nfft/2
            power(frequency) = power(frequency) + ABS(work(frequency + 1))**2
         END DO
      END DO
   END SUBROUTINE complex_power_sum

   SUBROUTINE make_window(name, window)
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      REAL(dp), INTENT(OUT)                              :: window(:)

      INTEGER                                            :: count, index
      REAL(dp)                                           :: pi

      count = SIZE(window)
      SELECT CASE (TRIM(name))
      CASE ("none")
         window = 1.0_dp
      CASE ("hann")
         IF (count == 1) THEN
            window = 1.0_dp
         ELSE
            pi = ACOS(-1.0_dp)
            DO index = 1, count
               window(index) = 0.5_dp*(1.0_dp - COS(2.0_dp*pi*REAL(index - 1, dp)/REAL(count - 1, dp)))
            END DO
         END IF
      CASE DEFAULT
         ERROR STOP "Unknown window"
      END SELECT
   END SUBROUTINE make_window

END MODULE trajana_time_series
