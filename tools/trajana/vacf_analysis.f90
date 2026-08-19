!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_vacf_analysis
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option,&
                                              has_flag
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_kinds,                   ONLY: dp
   USE trajana_selections,              ONLY: build_selection
   USE trajana_spectral_tools,          ONLY: frequency_axes,&
                                              write_peak_header,&
                                              write_spectral_peaks
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_time_series,             ONLY: real_autocorrelation_sum,&
                                              real_power_sum,&
                                              remove_real_means
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_vacf

CONTAINS

   SUBROUTINE run_vacf()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: input_path, output_path, peaks_path, &
                                                            selection, spectrum_path, window
      INTEGER :: atom, capacity, component, expected_atoms, first, frames, frequency, ierr, item, &
         lag, last, maximum_lag, nfft, output_unit, peaks_unit, selected_atoms, spectrum_unit, stride
      LOGICAL                                            :: eof, found, peaks_requested, remove_mean, &
                                                            spectrum_requested
      LOGICAL, ALLOCATABLE                               :: selected(:)
      REAL(dp)                                           :: dt_fs, effective_dt, energy_mev, peak_threshold, &
                                                            spectral_scale, spectral_value, terahertz, &
                                                            wavenumber, window_norm
      REAL(dp), ALLOCATABLE                              :: correlation(:), grown(:, :), power(:), &
                                                            series(:, :)
      TYPE(frame_type)                                   :: frame
      TYPE(xyz_reader_type)                              :: trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--spectrum", spectrum_path, spectrum_requested)
      CALL get_option("--peaks", peaks_path, peaks_requested)
      CALL get_option("--select", selection, found, "all")
      CALL get_option("--window", window, found, "none")
      window = lower_case(window)
      IF (window /= "none" .AND. window /= "hann") CALL fail("--window expects none or hann")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_real_option("--peak-threshold", peak_threshold, 0.05_dp)
      CALL get_integer_option("--max-lag", maximum_lag, -1)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      remove_mean = has_flag("--remove-mean")
      IF (dt_fs <= 0.0_dp) CALL fail("vacf requires a positive --dt-fs")
      IF (peak_threshold <= 0.0_dp .OR. peak_threshold > 1.0_dp) &
         CALL fail("--peak-threshold must be in the interval (0,1]")
      IF (peaks_requested .AND. .NOT. spectrum_requested) CALL fail("--peaks requires --spectrum")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")

      CALL trajectory%open_file(input_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      frames = 0
      capacity = 0
      selected_atoms = -1
      expected_atoms = -1

      DO
         CALL trajectory%read_frame(frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (eof .OR. frame%number > last) EXIT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         IF (expected_atoms < 0) expected_atoms = frame%natoms
         IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")
         CALL build_selection(selection, frame%label, selected, ierr, message)
         IF (ierr /= 0) CALL fail("Invalid --select: "//TRIM(message))
         IF (selected_atoms < 0) selected_atoms = COUNT(selected)
         IF (COUNT(selected) /= selected_atoms) CALL fail("The selected atom count changes along the trajectory")

         IF (frames == capacity) THEN
            IF (capacity == 0) THEN
               capacity = 128
               ALLOCATE (series(3*selected_atoms, capacity))
            ELSE
               ALLOCATE (grown(3*selected_atoms, 2*capacity))
               grown(:, :capacity) = series
               CALL MOVE_ALLOC(grown, series)
               capacity = 2*capacity
            END IF
         END IF
         frames = frames + 1
         item = 0
         DO atom = 1, frame%natoms
            IF (.NOT. selected(atom)) CYCLE
            DO component = 1, 3
               item = item + 1
               series(item, frames) = frame%value(component, atom)
            END DO
         END DO
      END DO
      CALL trajectory%close_file()
      IF (frames < 2) CALL fail("vacf requires at least two selected frames")
      IF (frames < capacity) THEN
         ALLOCATE (grown(3*selected_atoms, frames))
         grown = series(:, :frames)
         CALL MOVE_ALLOC(grown, series)
      END IF
      IF (remove_mean) CALL remove_real_means(series)
      IF (maximum_lag < 0) maximum_lag = frames - 1
      IF (maximum_lag >= frames) CALL fail("--max-lag must be smaller than the number of selected frames")

      CALL real_autocorrelation_sum(series, maximum_lag, correlation)
      DO lag = 0, maximum_lag
         correlation(lag) = correlation(lag)/(REAL(selected_atoms, dp)*REAL(frames - lag, dp))
      END DO
      IF (ABS(correlation(0)) <= TINY(1.0_dp)) CALL fail("VACF at zero lag is zero")

      effective_dt = dt_fs*REAL(stride, dp)
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)") "# lag [fs]   VACF   normalized_VACF"
      WRITE (output_unit, "(A,A)") "# selection: ", TRIM(selection)
      DO lag = 0, maximum_lag
         WRITE (output_unit, "(3ES24.16)") REAL(lag, dp)*effective_dt, correlation(lag), correlation(lag)/correlation(0)
      END DO
      CALL close_output(output_unit)

      IF (spectrum_requested) THEN
         CALL real_power_sum(series, window, power, nfft, window_norm)
         CALL open_output(spectrum_path, spectrum_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         WRITE (spectrum_unit, "(A)") &
            "# wavenumber [cm^-1]   frequency [THz]   energy [meV]   one-sided velocity power [velocity^2 fs]"
         WRITE (spectrum_unit, "(A,A)") "# window: ", TRIM(window)
         spectral_scale = effective_dt/(REAL(selected_atoms, dp)*window_norm)
         DO frequency = 0, nfft/2
            CALL frequency_axes(frequency, nfft, effective_dt, wavenumber, terahertz, energy_mev)
            spectral_value = power(frequency)*spectral_scale
            IF (frequency > 0 .AND. frequency < nfft/2) spectral_value = 2.0_dp*spectral_value
            WRITE (spectrum_unit, "(4ES24.16)") wavenumber, terahertz, energy_mev, spectral_value
         END DO
         CALL close_output(spectrum_unit)
         IF (peaks_requested) THEN
            CALL open_output(peaks_path, peaks_unit, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            CALL write_peak_header(peaks_unit)
            CALL write_spectral_peaks(peaks_unit, "vacf", 0.0_dp, power, spectral_scale, nfft, &
                                      effective_dt, peak_threshold)
            CALL close_output(peaks_unit)
         END IF
      END IF
   END SUBROUTINE run_vacf

END MODULE trajana_vacf_analysis
