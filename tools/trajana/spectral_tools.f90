!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_spectral_tools
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: frequency_axes, second_frequency_moment, write_peak_header, write_spectral_peaks

CONTAINS

   SUBROUTINE frequency_axes(index, nfft, dt_fs, wavenumber, terahertz, energy_mev)
      INTEGER, INTENT(IN)                                :: index, nfft
      REAL(dp), INTENT(IN)                               :: dt_fs
      REAL(dp), INTENT(OUT)                              :: wavenumber, terahertz, energy_mev

      REAL(dp)                                           :: frequency_hz

      frequency_hz = REAL(index, dp)/(REAL(nfft, dp)*dt_fs*1.0e-15_dp)
      terahertz = frequency_hz*1.0e-12_dp
      wavenumber = frequency_hz/2.99792458e10_dp
      energy_mev = terahertz*4.135667696_dp
   END SUBROUTINE frequency_axes

   REAL(dp) FUNCTION second_frequency_moment(power, nfft, dt_fs) RESULT(moment)
      REAL(dp), INTENT(IN)                               :: power(0:)
      INTEGER, INTENT(IN)                                :: nfft
      REAL(dp), INTENT(IN)                               :: dt_fs

      INTEGER                                            :: frequency
      REAL(dp)                                           :: denominator, omega

      moment = 0.0_dp
      denominator = SUM(power)
      IF (denominator <= TINY(1.0_dp)) RETURN
      DO frequency = 0, UBOUND(power, 1)
         omega = 2.0_dp*ACOS(-1.0_dp)*REAL(frequency, dp)/(REAL(nfft, dp)*dt_fs)
         moment = moment + omega**2*power(frequency)
      END DO
      moment = moment/denominator
   END FUNCTION second_frequency_moment

   SUBROUTINE write_peak_header(unit)
      INTEGER, INTENT(IN)                                :: unit

      WRITE (unit, "(A)") &
         "# channel   q [angstrom^-1]   wavenumber [cm^-1]   frequency [THz]   energy [meV]"//&
         "   phase_velocity [km/s]   intensity   FWHM [meV]   lifetime [fs]"
      WRITE (unit, "(A)") "# lifetime = 1/(pi*FWHM_frequency); unresolved widths are written as zero"
   END SUBROUTINE write_peak_header

   SUBROUTINE write_spectral_peaks(unit, channel, q_value, power, scale, nfft, dt_fs, threshold)
      INTEGER, INTENT(IN)                                :: unit, nfft
      CHARACTER(LEN=*), INTENT(IN)                       :: channel
      REAL(dp), INTENT(IN)                               :: q_value, power(0:), scale, dt_fs, threshold

      INTEGER                                            :: frequency, left, right
      REAL(dp) :: energy_mev, fraction, fwhm_mev, fwhm_thz, half_height, intensity, left_crossing, &
         lifetime_fs, maximum, phase_velocity, right_crossing, terahertz, wavenumber

      IF (UBOUND(power, 1) < 2) RETURN
      maximum = MAXVAL(power(1:UBOUND(power, 1)))
      IF (maximum <= TINY(1.0_dp)) RETURN
      DO frequency = 1, UBOUND(power, 1) - 1
         IF (power(frequency) < threshold*maximum) CYCLE
         IF (power(frequency) < power(frequency - 1) .OR. power(frequency) <= power(frequency + 1)) CYCLE

         half_height = 0.5_dp*power(frequency)
         left = frequency
         DO WHILE (left > 0 .AND. power(left) > half_height)
            left = left - 1
         END DO
         right = frequency
         DO WHILE (right < UBOUND(power, 1) .AND. power(right) > half_height)
            right = right + 1
         END DO
         fwhm_thz = 0.0_dp
         lifetime_fs = 0.0_dp
         IF (left < frequency .AND. right > frequency) THEN
            fraction = (half_height - power(left))/MAX(power(left + 1) - power(left), TINY(1.0_dp))
            left_crossing = REAL(left, dp) + fraction
            fraction = (power(right - 1) - half_height)/MAX(power(right - 1) - power(right), TINY(1.0_dp))
            right_crossing = REAL(right, dp) - fraction
            fwhm_thz = (right_crossing - left_crossing)*1000.0_dp/(REAL(nfft, dp)*dt_fs)
            IF (fwhm_thz > TINY(1.0_dp)) lifetime_fs = 1000.0_dp/(ACOS(-1.0_dp)*fwhm_thz)
         END IF
         CALL frequency_axes(frequency, nfft, dt_fs, wavenumber, terahertz, energy_mev)
         fwhm_mev = fwhm_thz*4.135667696_dp
         intensity = power(frequency)*scale
         phase_velocity = 0.0_dp
         IF (q_value > TINY(1.0_dp)) phase_velocity = 0.2_dp*ACOS(-1.0_dp)*terahertz/q_value
         WRITE (unit, "(A,1X,8ES24.16)") TRIM(channel), q_value, wavenumber, terahertz, energy_mev, &
            phase_velocity, intensity, fwhm_mev, lifetime_fs
      END DO
   END SUBROUTINE write_spectral_peaks

END MODULE trajana_spectral_tools
