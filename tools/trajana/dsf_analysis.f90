!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_dsf_analysis
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option,&
                                              has_flag
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_kinds,                   ONLY: dp
   USE trajana_selections,              ONLY: build_selection
   USE trajana_spectral_tools,          ONLY: frequency_axes,&
                                              second_frequency_moment,&
                                              write_peak_header,&
                                              write_spectral_peaks
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_time_series,             ONLY: complex_autocorrelation_sum,&
                                              complex_power_sum,&
                                              remove_complex_means
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: run_dsf

CONTAINS

   SUBROUTINE run_dsf()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE :: current_path, current_spectrum_path, input_path, output_path, &
         peaks_path, q_path, selection, self_output_path, self_spectrum_path, spectrum_path, summary_path, &
         velocity_path, window
      COMPLEX(dp), ALLOCATABLE :: correlation(:), correlation_l(:), correlation_t(:), current_l(:, :), &
         current_t(:, :), density(:, :), phases(:, :), self_correlation(:)
      INTEGER :: atom, capacity, expected_atoms, first, frames, ierr, iq, item, last, maximum_lag, nfft, &
         nq, output_unit, selected_atoms, stride
      LOGICAL :: current_output_requested, current_spectrum_requested, eof, found, have_velocity, &
         peaks_requested, remove_mean, self_output_requested, self_requested, self_spectrum_requested, &
         spectrum_requested, summary_requested, velocity_eof
      LOGICAL, ALLOCATABLE                               :: selected(:)
      REAL(dp) :: current_scale, current_window_norm, density_scale, density_window_norm, dt_fs, effective_dt, &
         mass_density, mean_q, omega0_squared, omega_l_squared, omega_t_squared, peak_threshold, phase, &
         self_scale, self_window_norm, transverse_scale, velocity_scale
      REAL(dp), ALLOCATABLE :: basis_one(:, :), basis_two(:, :), density_power(:), longitudinal_power(:), &
         q(:, :), self_power(:), transverse_power(:)
      TYPE(frame_type)                                   :: frame, velocity_frame
      TYPE(xyz_reader_type)                              :: trajectory, velocity_trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--spectrum", spectrum_path, spectrum_requested)
      CALL get_option("--self-output", self_output_path, self_output_requested)
      CALL get_option("--self-spectrum", self_spectrum_path, self_spectrum_requested)
      CALL get_option("--velocity", velocity_path, have_velocity)
      CALL get_option("--current", current_path, current_output_requested)
      CALL get_option("--current-spectrum", current_spectrum_path, current_spectrum_requested)
      CALL get_option("--peaks", peaks_path, peaks_requested)
      CALL get_option("--summary", summary_path, summary_requested)
      CALL get_option("--q-file", q_path, found)
      IF (.NOT. found) CALL fail("dsf requires --q-file FILE")
      CALL get_option("--select", selection, found, "all")
      CALL get_option("--window", window, found, "none")
      window = lower_case(window)
      IF (window /= "none" .AND. window /= "hann") CALL fail("--window expects none or hann")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_real_option("--velocity-scale", velocity_scale, 1.0_dp)
      CALL get_real_option("--mass-density", mass_density, -1.0_dp)
      CALL get_real_option("--peak-threshold", peak_threshold, 0.05_dp)
      CALL get_integer_option("--max-lag", maximum_lag, -1)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      remove_mean = has_flag("--remove-mean")
      self_requested = self_output_requested .OR. self_spectrum_requested
      IF (dt_fs <= 0.0_dp) CALL fail("dsf requires a positive --dt-fs")
      IF (velocity_scale <= 0.0_dp) CALL fail("--velocity-scale must be positive")
      IF (peak_threshold <= 0.0_dp .OR. peak_threshold > 1.0_dp) &
         CALL fail("--peak-threshold must be in the interval (0,1]")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")
      IF ((current_output_requested .OR. current_spectrum_requested) .AND. .NOT. have_velocity) &
         CALL fail("Current analysis requires --velocity FILE")

      CALL read_q_vectors(q_path, q, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      nq = SIZE(q, 2)
      CALL prepare_q_shell(q, mean_q, basis_one, basis_two, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL trajectory%open_file(input_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (have_velocity) THEN
         IF (TRIM(input_path) == "-" .OR. TRIM(velocity_path) == "-") &
            CALL fail("Position and velocity trajectories cannot both use standard input")
         CALL velocity_trajectory%open_file(velocity_path, ierr, message)
         IF (ierr /= 0) CALL fail(message)
      END IF
      frames = 0
      capacity = 0
      selected_atoms = -1
      expected_atoms = -1

      DO
         CALL trajectory%read_frame(frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         velocity_eof = eof
         IF (have_velocity) THEN
            CALL velocity_trajectory%read_frame(velocity_frame, velocity_eof, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            IF (eof .NEQV. velocity_eof) CALL fail("Position and velocity trajectories have different lengths")
         END IF
         IF (eof) EXIT
         IF (frame%number > last) EXIT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         IF (expected_atoms < 0) expected_atoms = frame%natoms
         IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")
         IF (have_velocity) CALL validate_velocity_frame(frame, velocity_frame)
         CALL build_selection(selection, frame%label, selected, ierr, message)
         IF (ierr /= 0) CALL fail("Invalid --select: "//TRIM(message))
         IF (selected_atoms < 0) selected_atoms = COUNT(selected)
         IF (selected_atoms == 0) CALL fail("The selection contains no atoms")
         IF (COUNT(selected) /= selected_atoms) CALL fail("The selected atom count changes along the trajectory")

         IF (frames == capacity) THEN
            IF (capacity == 0) THEN
               capacity = 128
               ALLOCATE (density(nq, capacity))
               IF (self_requested) ALLOCATE (phases(nq*selected_atoms, capacity))
               IF (have_velocity) ALLOCATE (current_l(nq, capacity), current_t(2*nq, capacity))
            ELSE
               capacity = 2*capacity
               CALL grow_complex_matrix(density, capacity, frames)
               IF (self_requested) CALL grow_complex_matrix(phases, capacity, frames)
               IF (have_velocity) THEN
                  CALL grow_complex_matrix(current_l, capacity, frames)
                  CALL grow_complex_matrix(current_t, capacity, frames)
               END IF
            END IF
         END IF
         frames = frames + 1
         density(:, frames) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         IF (have_velocity) THEN
            current_l(:, frames) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            current_t(:, frames) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         END IF
         item = 0
         DO atom = 1, frame%natoms
            IF (.NOT. selected(atom)) CYCLE
            item = item + 1
            DO iq = 1, nq
               phase = DOT_PRODUCT(q(:, iq), frame%value(:, atom))
               CALL accumulate_modes(iq, item, selected_atoms, phase, velocity_scale, velocity_frame, atom, &
                                     basis_one, basis_two, density, phases, current_l, current_t, frames, &
                                     self_requested, have_velocity)
            END DO
         END DO
      END DO
      CALL trajectory%close_file()
      IF (have_velocity) CALL velocity_trajectory%close_file()
      IF (frames < 2) CALL fail("dsf requires at least two selected frames")
      CALL trim_complex_matrix(density, frames)
      IF (self_requested) CALL trim_complex_matrix(phases, frames)
      IF (have_velocity) THEN
         CALL trim_complex_matrix(current_l, frames)
         CALL trim_complex_matrix(current_t, frames)
      END IF
      IF (remove_mean) THEN
         CALL remove_complex_means(density)
         IF (self_requested) CALL remove_complex_means(phases)
         IF (have_velocity) THEN
            CALL remove_complex_means(current_l)
            CALL remove_complex_means(current_t)
         END IF
      END IF
      IF (maximum_lag < 0) maximum_lag = frames - 1
      IF (maximum_lag >= frames) CALL fail("--max-lag must be smaller than the number of selected frames")
      effective_dt = dt_fs*REAL(stride, dp)

      CALL normalized_correlation(density, maximum_lag, REAL(nq*selected_atoms, dp), correlation)
      IF (ABS(correlation(0)) <= TINY(1.0_dp)) CALL fail("Intermediate scattering function at zero lag is zero")
      CALL write_correlation(output_path, "F(q,t)", selection, mean_q, effective_dt, correlation)

      IF (self_requested) THEN
         CALL normalized_correlation(phases, maximum_lag, REAL(nq*selected_atoms, dp), self_correlation)
         IF (self_output_requested) &
            CALL write_correlation(self_output_path, "F_self(q,t)", selection, mean_q, effective_dt, self_correlation)
      END IF
      IF (have_velocity) THEN
         CALL normalized_correlation(current_l, maximum_lag, REAL(nq*selected_atoms, dp), correlation_l)
         CALL normalized_correlation(current_t, maximum_lag, REAL(2*nq*selected_atoms, dp), correlation_t)
         IF (current_output_requested) &
            CALL write_current_correlation(current_path, selection, mean_q, effective_dt, correlation_l, correlation_t)
      END IF

      IF (spectrum_requested .OR. current_spectrum_requested .OR. peaks_requested) THEN
         CALL complex_power_sum(density, window, density_power, nfft, density_window_norm)
         density_scale = effective_dt/(REAL(nq*selected_atoms, dp)*density_window_norm)
      END IF
      IF (self_spectrum_requested .OR. (peaks_requested .AND. self_requested)) THEN
         CALL complex_power_sum(phases, window, self_power, nfft, self_window_norm)
         self_scale = effective_dt/(REAL(nq*selected_atoms, dp)*self_window_norm)
      END IF
      IF (have_velocity .AND. (current_spectrum_requested .OR. peaks_requested .OR. summary_requested)) THEN
         CALL complex_power_sum(current_l, window, longitudinal_power, nfft, current_window_norm)
         CALL complex_power_sum(current_t, window, transverse_power, nfft, current_window_norm)
         current_scale = effective_dt/(REAL(nq*selected_atoms, dp)*current_window_norm)
         transverse_scale = 0.5_dp*current_scale
      END IF

      IF (spectrum_requested) &
         CALL write_spectrum(spectrum_path, "S_coherent", mean_q, window, density_power, density_scale, nfft, effective_dt)
      IF (self_spectrum_requested) &
         CALL write_spectrum(self_spectrum_path, "S_incoherent", mean_q, window, self_power, self_scale, nfft, effective_dt)
      IF (have_velocity .AND. current_spectrum_requested) &
         CALL write_current_spectrum(current_spectrum_path, mean_q, window, density_power, density_scale, &
                                     longitudinal_power, current_scale, transverse_power, transverse_scale, &
                                     effective_dt, nfft)

      IF (peaks_requested) THEN
         CALL open_output(peaks_path, output_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL write_peak_header(output_unit)
         CALL write_spectral_peaks(output_unit, "coherent", mean_q, density_power, density_scale, nfft, &
                                   effective_dt, peak_threshold)
         IF (self_requested) &
            CALL write_spectral_peaks(output_unit, "incoherent", mean_q, self_power, self_scale, nfft, &
                                      effective_dt, peak_threshold)
         IF (have_velocity) THEN
            CALL write_spectral_peaks(output_unit, "longitudinal", mean_q, longitudinal_power, current_scale, &
                                      nfft, effective_dt, peak_threshold)
            CALL write_spectral_peaks(output_unit, "transverse", mean_q, transverse_power, transverse_scale, &
                                      nfft, effective_dt, peak_threshold)
         END IF
         CALL close_output(output_unit)
      END IF

      IF (summary_requested) THEN
         omega0_squared = 0.0_dp
         omega_l_squared = 0.0_dp
         omega_t_squared = 0.0_dp
         IF (have_velocity) THEN
            omega0_squared = mean_q**2*REAL(correlation_l(0), dp)/REAL(correlation(0), dp)
            omega_l_squared = second_frequency_moment(longitudinal_power, nfft, effective_dt)
            omega_t_squared = second_frequency_moment(transverse_power, nfft, effective_dt)
         END IF
         CALL write_summary(summary_path, mean_q, REAL(correlation(0), dp), correlation_l, correlation_t, &
                            omega0_squared, omega_l_squared, omega_t_squared, mass_density, have_velocity)
      END IF
   END SUBROUTINE run_dsf

   SUBROUTINE accumulate_modes(iq, item, selected_atoms, phase, velocity_scale, velocity_frame, atom, &
                               basis_one, basis_two, density, phases, current_l, current_t, frame_index, &
                               self_requested, have_velocity)
      INTEGER, INTENT(IN)                                :: iq, item, selected_atoms, atom, frame_index
      REAL(dp), INTENT(IN)                               :: phase, velocity_scale
      TYPE(frame_type), INTENT(IN)                       :: velocity_frame
      REAL(dp), INTENT(IN)                               :: basis_one(:, :), basis_two(:, :)
      COMPLEX(dp), INTENT(INOUT)                         :: density(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT)            :: phases(:, :), current_l(:, :), current_t(:, :)
      LOGICAL, INTENT(IN)                                :: self_requested, have_velocity

      COMPLEX(dp)                                        :: exponential
      REAL(dp)                                           :: velocity(3)

      exponential = CMPLX(COS(phase), SIN(phase), KIND=dp)
      density(iq, frame_index) = density(iq, frame_index) + exponential
      IF (self_requested) phases((iq - 1)*selected_atoms + item, frame_index) = exponential
      IF (have_velocity) THEN
         velocity = velocity_scale*velocity_frame%value(:, atom)
         current_l(iq, frame_index) = current_l(iq, frame_index) + &
                                      DOT_PRODUCT(velocity, cross_product(basis_one(:, iq), basis_two(:, iq)))*exponential
         current_t(2*iq - 1, frame_index) = current_t(2*iq - 1, frame_index) + &
                                            DOT_PRODUCT(velocity, basis_one(:, iq))*exponential
         current_t(2*iq, frame_index) = current_t(2*iq, frame_index) + &
                                        DOT_PRODUCT(velocity, basis_two(:, iq))*exponential
      END IF
   END SUBROUTINE accumulate_modes

   SUBROUTINE normalized_correlation(series, maximum_lag, normalization, correlation)
      COMPLEX(dp), INTENT(IN)                            :: series(:, :)
      INTEGER, INTENT(IN)                                :: maximum_lag
      REAL(dp), INTENT(IN)                               :: normalization
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: correlation(:)

      INTEGER                                            :: lag, frames

      frames = SIZE(series, 2)
      CALL complex_autocorrelation_sum(series, maximum_lag, correlation)
      DO lag = 0, maximum_lag
         correlation(lag) = correlation(lag)/(normalization*REAL(frames - lag, dp))
      END DO
   END SUBROUTINE normalized_correlation

   SUBROUTINE write_correlation(path, name, selection, q_value, dt_fs, correlation)
      CHARACTER(LEN=*), INTENT(IN)                       :: path, name, selection
      REAL(dp), INTENT(IN)                               :: q_value, dt_fs
      COMPLEX(dp), INTENT(IN)                            :: correlation(0:)

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: ierr, lag, unit

      CALL open_output(path, unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (unit, "(A)") "# lag [fs]   Re(correlation)   Im(correlation)   normalized_real"
      WRITE (unit, "(A,A)") "# quantity: ", TRIM(name)
      WRITE (unit, "(A,A)") "# selection: ", TRIM(selection)
      WRITE (unit, "(A,ES24.16)") "# q [angstrom^-1]: ", q_value
      DO lag = 0, UBOUND(correlation, 1)
         WRITE (unit, "(4ES24.16)") REAL(lag, dp)*dt_fs, REAL(correlation(lag), dp), &
            AIMAG(correlation(lag)), REAL(correlation(lag)/correlation(0), dp)
      END DO
      CALL close_output(unit)
   END SUBROUTINE write_correlation

   SUBROUTINE write_current_correlation(path, selection, q_value, dt_fs, longitudinal, transverse)
      CHARACTER(LEN=*), INTENT(IN)                       :: path, selection
      REAL(dp), INTENT(IN)                               :: q_value, dt_fs
      COMPLEX(dp), INTENT(IN)                            :: longitudinal(0:), transverse(0:)

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: ierr, lag, unit

      CALL open_output(path, unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (unit, "(A)") &
         "# lag [fs]   Re(C_L)   Im(C_L)   normalized_Re(C_L)   Re(C_T)   Im(C_T)   normalized_Re(C_T)"
      WRITE (unit, "(A,A)") "# selection: ", TRIM(selection)
      WRITE (unit, "(A,ES24.16)") "# q [angstrom^-1]: ", q_value
      DO lag = 0, UBOUND(longitudinal, 1)
         WRITE (unit, "(7ES24.16)") REAL(lag, dp)*dt_fs, REAL(longitudinal(lag), dp), AIMAG(longitudinal(lag)), &
            REAL(longitudinal(lag)/longitudinal(0), dp), REAL(transverse(lag), dp), AIMAG(transverse(lag)), &
            REAL(transverse(lag)/transverse(0), dp)
      END DO
      CALL close_output(unit)
   END SUBROUTINE write_current_correlation

   SUBROUTINE write_spectrum(path, name, q_value, window, power, scale, nfft, dt_fs)
      CHARACTER(LEN=*), INTENT(IN)                       :: path, name, window
      REAL(dp), INTENT(IN)                               :: q_value, power(0:), scale, dt_fs
      INTEGER, INTENT(IN)                                :: nfft

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: frequency, ierr, unit
      REAL(dp)                                           :: energy_mev, omega, spectral_value, terahertz, wavenumber

      CALL open_output(path, unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (unit, "(A)") &
         "# wavenumber [cm^-1]   frequency [THz]   energy [meV]   spectral_density [fs]   omega^2*S/q^2"
      WRITE (unit, "(A,A)") "# quantity: ", TRIM(name)
      WRITE (unit, "(A,A)") "# window: ", TRIM(window)
      WRITE (unit, "(A,ES24.16)") "# q [angstrom^-1]: ", q_value
      DO frequency = 0, UBOUND(power, 1)
         CALL frequency_axes(frequency, nfft, dt_fs, wavenumber, terahertz, energy_mev)
         omega = 2.0_dp*ACOS(-1.0_dp)*REAL(frequency, dp)/(REAL(nfft, dp)*dt_fs)
         spectral_value = power(frequency)*scale
         WRITE (unit, "(5ES24.16)") wavenumber, terahertz, energy_mev, spectral_value, &
            omega**2*spectral_value/q_value**2
      END DO
      CALL close_output(unit)
   END SUBROUTINE write_spectrum

   SUBROUTINE write_current_spectrum(path, q_value, window, density_power, density_scale, &
                                     longitudinal_power, longitudinal_scale, transverse_power, &
                                     transverse_scale, dt_fs, nfft)
      CHARACTER(LEN=*), INTENT(IN)                       :: path, window
      REAL(dp), INTENT(IN) :: q_value, density_power(0:), density_scale, longitudinal_power(0:), &
         longitudinal_scale, transverse_power(0:), transverse_scale, dt_fs
      INTEGER, INTENT(IN)                                :: nfft

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: frequency, ierr, unit
      REAL(dp)                                           :: energy_mev, omega, terahertz, wavenumber

      CALL open_output(path, unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (unit, "(A)") &
         "# wavenumber [cm^-1]   frequency [THz]   energy [meV]   C_L   C_T   omega^2*S/q^2"
      WRITE (unit, "(A,A)") "# window: ", TRIM(window)
      WRITE (unit, "(A,ES24.16)") "# q [angstrom^-1]: ", q_value
      DO frequency = 0, UBOUND(longitudinal_power, 1)
         CALL frequency_axes(frequency, nfft, dt_fs, wavenumber, terahertz, energy_mev)
         omega = 2.0_dp*ACOS(-1.0_dp)*REAL(frequency, dp)/(REAL(nfft, dp)*dt_fs)
         WRITE (unit, "(6ES24.16)") wavenumber, terahertz, energy_mev, &
            longitudinal_power(frequency)*longitudinal_scale, transverse_power(frequency)*transverse_scale, &
            omega**2*density_power(frequency)*density_scale/q_value**2
      END DO
      CALL close_output(unit)
   END SUBROUTINE write_current_spectrum

   SUBROUTINE write_summary(path, q_value, structure_factor, longitudinal, transverse, omega0_squared, &
                            omega_l_squared, omega_t_squared, mass_density, have_velocity)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      REAL(dp), INTENT(IN) :: q_value, structure_factor, omega0_squared, omega_l_squared, omega_t_squared, &
         mass_density
      COMPLEX(dp), ALLOCATABLE, INTENT(IN)               :: longitudinal(:), transverse(:)
      LOGICAL, INTENT(IN)                                :: have_velocity

      CHARACTER(LEN=512)                                 :: message
      INTEGER                                            :: ierr, unit
      REAL(dp)                                           :: c_l, c_t, modulus_l, modulus_t, value_l, value_t

      value_l = 0.0_dp
      value_t = 0.0_dp
      c_l = 0.0_dp
      c_t = 0.0_dp
      modulus_l = 0.0_dp
      modulus_t = 0.0_dp
      IF (have_velocity) THEN
         value_l = REAL(longitudinal(0), dp)
         value_t = REAL(transverse(0), dp)
         c_l = SQRT(MAX(omega_l_squared, 0.0_dp))/q_value
         c_t = SQRT(MAX(omega_t_squared, 0.0_dp))/q_value
         IF (mass_density > 0.0_dp) THEN
            modulus_l = 1.0e4_dp*mass_density*c_l**2
            modulus_t = 1.0e4_dp*mass_density*c_t**2
         END IF
      END IF
      CALL open_output(path, unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (unit, "(A)") &
         "# q   S(q)   C_L(0)   C_T(0)   omega0^2   omegaL^2   omegaT^2   cL [A/fs]   cT [A/fs]"//&
         "   rho*cL^2 [GPa]   rho*cT^2 [GPa]"
      WRITE (unit, "(A)") "# angular-frequency moments are in fs^-2; moduli are zero unless --mass-density is given"
      WRITE (unit, "(11ES24.16)") q_value, structure_factor, value_l, value_t, omega0_squared, &
         omega_l_squared, omega_t_squared, c_l, c_t, modulus_l, modulus_t
      CALL close_output(unit)
   END SUBROUTINE write_summary

   SUBROUTINE prepare_q_shell(q, mean_q, basis_one, basis_two, ierr, message)
      REAL(dp), INTENT(IN)                               :: q(:, :)
      REAL(dp), INTENT(OUT)                              :: mean_q
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: basis_one(:, :), basis_two(:, :)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: iq
      REAL(dp)                                           :: q_hat(3), q_norm, reference(3)

      ierr = 0
      message = ""
      mean_q = 0.0_dp
      DO iq = 1, SIZE(q, 2)
         mean_q = mean_q + NORM2(q(:, iq))
      END DO
      mean_q = mean_q/REAL(SIZE(q, 2), dp)
      DO iq = 1, SIZE(q, 2)
         q_norm = NORM2(q(:, iq))
         IF (ABS(q_norm - mean_q) > 1.0e-6_dp*MAX(1.0_dp, mean_q)) THEN
            ierr = 1
            message = "All q vectors must belong to one magnitude shell; use one run per q shell"
            RETURN
         END IF
      END DO
      ALLOCATE (basis_one(3, SIZE(q, 2)), basis_two(3, SIZE(q, 2)))
      DO iq = 1, SIZE(q, 2)
         q_hat = q(:, iq)/NORM2(q(:, iq))
         IF (ABS(q_hat(1)) < 0.9_dp) THEN
            reference = [1.0_dp, 0.0_dp, 0.0_dp]
         ELSE
            reference = [0.0_dp, 1.0_dp, 0.0_dp]
         END IF
         basis_one(:, iq) = cross_product(q_hat, reference)
         basis_one(:, iq) = basis_one(:, iq)/NORM2(basis_one(:, iq))
         basis_two(:, iq) = cross_product(q_hat, basis_one(:, iq))
      END DO
   END SUBROUTINE prepare_q_shell

   PURE FUNCTION cross_product(first, second) RESULT(product)
      REAL(dp), INTENT(IN)                               :: first(3), second(3)
      REAL(dp)                                           :: product(3)

      product = [first(2)*second(3) - first(3)*second(2), &
                 first(3)*second(1) - first(1)*second(3), &
                 first(1)*second(2) - first(2)*second(1)]
   END FUNCTION cross_product

   SUBROUTINE validate_velocity_frame(frame, velocity_frame)
      TYPE(frame_type), INTENT(IN)                       :: frame, velocity_frame

      INTEGER                                            :: atom

      IF (velocity_frame%natoms /= frame%natoms) &
         CALL fail("Position and velocity frames have different atom counts")
      DO atom = 1, frame%natoms
         IF (TRIM(frame%label(atom)) /= TRIM(velocity_frame%label(atom))) &
            CALL fail("Position and velocity trajectories have different atom ordering")
      END DO
   END SUBROUTINE validate_velocity_frame

   SUBROUTINE grow_complex_matrix(series, new_capacity, used)
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT)            :: series(:, :)
      INTEGER, INTENT(IN)                                :: new_capacity, used

      COMPLEX(dp), ALLOCATABLE                           :: grown(:, :)

      ALLOCATE (grown(SIZE(series, 1), new_capacity))
      grown(:, :used) = series(:, :used)
      CALL MOVE_ALLOC(grown, series)
   END SUBROUTINE grow_complex_matrix

   SUBROUTINE trim_complex_matrix(series, frames)
      COMPLEX(dp), ALLOCATABLE, INTENT(INOUT)            :: series(:, :)
      INTEGER, INTENT(IN)                                :: frames

      COMPLEX(dp), ALLOCATABLE                           :: trimmed(:, :)

      IF (SIZE(series, 2) == frames) RETURN
      ALLOCATE (trimmed(SIZE(series, 1), frames))
      trimmed = series(:, :frames)
      CALL MOVE_ALLOC(trimmed, series)
   END SUBROUTINE trim_complex_matrix

   SUBROUTINE read_q_vectors(path, q, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: q(:, :)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=4096)                                :: line
      INTEGER                                            :: count, ios, unit
      REAL(dp)                                           :: vector(3)
      REAL(dp), ALLOCATABLE                              :: grown(:, :)

      ierr = 0
      message = ""
      count = 0
      ALLOCATE (q(3, 0))
      OPEN (NEWUNIT=unit, FILE=TRIM(path), STATUS="old", ACTION="read", IOSTAT=ios)
      IF (ios /= 0) THEN
         ierr = ios
         message = "Cannot open q-vector file "//TRIM(path)
         RETURN
      END IF
      DO
         READ (unit, "(A)", IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) EXIT
         IF (ios /= 0) THEN
            ierr = ios
            message = "Error while reading q-vector file"
            EXIT
         END IF
         line = ADJUSTL(line)
         IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
         READ (line, *, IOSTAT=ios) vector
         IF (ios /= 0 .OR. NORM2(vector) <= TINY(1.0_dp)) THEN
            ierr = 1
            message = "Each q-vector line must contain three nonzero Cartesian components"
            EXIT
         END IF
         ALLOCATE (grown(3, count + 1))
         IF (count > 0) grown(:, :count) = q
         grown(:, count + 1) = vector
         CALL MOVE_ALLOC(grown, q)
         count = count + 1
      END DO
      CLOSE (unit)
      IF (ierr == 0 .AND. count == 0) THEN
         ierr = 1
         message = "The q-vector file contains no vectors"
      END IF
   END SUBROUTINE read_q_vectors

END MODULE trajana_dsf_analysis
