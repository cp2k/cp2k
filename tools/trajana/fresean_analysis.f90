!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_fresean_analysis
   USE trajana_alignment,               ONLY: center_positions,&
                                              fit_rotation,&
                                              identity3
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option,&
                                              has_flag
   USE trajana_fft,                     ONLY: fft_any_in_place
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_groups,                  ONLY: mass_table_type,&
                                              read_masses
   USE trajana_kinds,                   ONLY: dp
   USE trajana_linear_algebra,          ONLY: diagonalize_symmetric
   USE trajana_selections,              ONLY: build_selection
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_time_series,             ONLY: complex_autocorrelation_sum,&
                                              remove_real_means
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   REAL(dp), PARAMETER :: bohr_per_atomic_time_to_angstrom_per_ps = 21876.91263641118_dp
   REAL(dp), PARAMETER :: speed_of_light_cm_per_ps = 0.0299792458_dp
   REAL(dp), PARAMETER :: mass_velocity_to_j_mol = 10.0_dp
   REAL(dp), PARAMETER :: gas_constant_j_mol_k = 8.31446261815324_dp
   REAL(dp), PARAMETER :: thz_to_wavenumber = 33.3564095198152_dp
   REAL(dp), PARAMETER :: wavenumber_to_kj_mol = 0.0119626565638697_dp

   PUBLIC :: print_fresean_help, run_fresean

CONTAINS

   SUBROUTINE run_fresean()
      CHARACTER(LEN=32), ALLOCATABLE                     :: labels(:)
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE :: mass_path, mode_path, mode_spectrum_path, &
         mode_timeseries_path, output_path, position_path, reference_path, selection, &
         velocity_path, velocity_unit
      INTEGER                                            :: atom, capacity, component, first, &
                                                            frames, ierr, item, last, mode_count, &
                                                            n_correlation, n_dof, n_selected, &
                                                            stride
      INTEGER, ALLOCATABLE                               :: selected_index(:)
      LOGICAL :: eof, found, mode_requested, mode_spectrum_requested, mode_timeseries_requested, &
         position_eof, position_found, reference_eof, reference_found, remove_mean
      LOGICAL, ALLOCATABLE                               :: selected(:)
      REAL(dp) :: aligned_position(3), constraints, dt_fs, effective_dt_ps, frequency_cm, &
         rotation(3, 3), sigma_cm, velocity(3), velocity_scale, velocity_scale_extra
      REAL(dp), ALLOCATABLE :: current_positions(:, :), grown(:, :), grown_positions(:, :), &
         masses(:), position_series(:, :), reference_positions(:, :), series(:, :)
      TYPE(frame_type)                                   :: position_frame, reference_frame, &
                                                            velocity_frame
      TYPE(mass_table_type)                              :: mass_table
      TYPE(xyz_reader_type)                              :: positions, reference, velocities

      CALL get_option("--velocity", velocity_path, found)
      CALL get_option("--position", position_path, position_found)
      CALL get_option("--reference", reference_path, reference_found)
      CALL get_option("--mass-file", mass_path, found)
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--mode-file", mode_path, mode_requested)
      CALL get_option("--mode-spectrum", mode_spectrum_path, mode_spectrum_requested)
      CALL get_option("--mode-timeseries", mode_timeseries_path, mode_timeseries_requested)
      CALL get_option("--select", selection, found, "all")
      CALL get_option("--velocity-unit", velocity_unit, found, "cp2k")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_real_option("--sigma-cm", sigma_cm, 10.0_dp)
      CALL get_real_option("--frequency-cm", frequency_cm, 0.0_dp)
      CALL get_real_option("--velocity-scale", velocity_scale_extra, 1.0_dp)
      CALL get_real_option("--constraints", constraints, 0.0_dp)
      CALL get_integer_option("--correlation-frames", n_correlation, 500)
      CALL get_integer_option("--mode-count", mode_count, 10)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      remove_mean = has_flag("--remove-mean")

      IF (LEN_TRIM(velocity_path) == 0) CALL fail("fresean.x requires --velocity")
      IF (LEN_TRIM(mass_path) == 0) CALL fail("fresean.x requires --mass-file")
      IF (dt_fs <= 0.0_dp) CALL fail("fresean.x requires a positive --dt-fs")
      IF (sigma_cm <= 0.0_dp) CALL fail("--sigma-cm must be positive")
      IF (frequency_cm < 0.0_dp) CALL fail("--frequency-cm cannot be negative")
      IF (velocity_scale_extra <= 0.0_dp) CALL fail("--velocity-scale must be positive")
      IF (n_correlation < 2) CALL fail("--correlation-frames must be at least two")
      IF (mode_count < 1) CALL fail("--mode-count must be positive")
      IF (constraints < 0.0_dp) CALL fail("--constraints cannot be negative")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")
      IF (reference_found .AND. .NOT. position_found) CALL fail("--reference requires --position")
      IF (mode_timeseries_requested .AND. .NOT. position_found) &
         CALL fail("--mode-timeseries requires --position")
      IF (mode_timeseries_requested .AND. frequency_cm <= 0.0_dp) &
         CALL fail("--mode-timeseries requires a positive --frequency-cm")

      velocity_unit = lower_case(velocity_unit)
      SELECT CASE (velocity_unit)
      CASE ("cp2k", "bohr/au_time", "bohr*au_t^-1")
         velocity_scale = bohr_per_atomic_time_to_angstrom_per_ps
      CASE ("angstrom/ps", "a/ps")
         velocity_scale = 1.0_dp
      CASE ("angstrom/fs", "a/fs")
         velocity_scale = 1000.0_dp
      CASE DEFAULT
         CALL fail("--velocity-unit expects cp2k, angstrom/ps, or angstrom/fs")
      END SELECT
      velocity_scale = velocity_scale*velocity_scale_extra

      CALL read_masses(mass_path, mass_table, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (position_found) THEN
         CALL positions%open_file(position_path, ierr, message)
         IF (ierr /= 0) CALL fail(message)
      END IF
      IF (reference_found) THEN
         CALL reference%open_file(reference_path, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL reference%read_frame(reference_frame, reference_eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (reference_eof) CALL fail("The reference trajectory is empty")
         CALL reference%close_file()
      END IF
      CALL velocities%open_file(velocity_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)

      frames = 0
      capacity = 0
      n_selected = -1
      DO
         CALL velocities%read_frame(velocity_frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (position_found) THEN
            CALL positions%read_frame(position_frame, position_eof, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            IF (eof .NEQV. position_eof) CALL fail("Position and velocity trajectories have different lengths")
         END IF
         IF (eof .OR. velocity_frame%number > last) EXIT
         IF (.NOT. frame_selected(velocity_frame%number, first, last, stride)) CYCLE

         IF (frames == 0) THEN
            CALL build_selection(selection, velocity_frame%label, selected, ierr, message)
            IF (ierr /= 0) CALL fail("Invalid --select: "//TRIM(message))
            n_selected = COUNT(selected)
            IF (n_selected == 0) CALL fail("The selection contains no atoms")
            n_dof = 3*n_selected
            IF (constraints >= REAL(n_dof, dp)) CALL fail("--constraints must be smaller than the selected DOF")
            ALLOCATE (selected_index(n_selected), labels(n_selected), masses(n_selected))
            selected_index = PACK([(atom, atom=1, velocity_frame%natoms)], selected)
            DO item = 1, n_selected
               atom = selected_index(item)
               labels(item) = velocity_frame%label(atom)
               CALL mass_table%lookup(labels(item), masses(item), found)
               IF (.NOT. found) CALL fail("No mass for trajectory label "//TRIM(labels(item)))
            END DO
            capacity = 128
            ALLOCATE (series(n_dof, capacity))
            IF (mode_timeseries_requested) THEN
               ALLOCATE (position_series(n_dof, capacity))
            ELSE
               ALLOCATE (position_series(0, 0))
            END IF
            IF (position_found) THEN
               CALL validate_position_frame(velocity_frame, position_frame)
               ALLOCATE (reference_positions(3, n_selected), current_positions(3, n_selected))
               IF (reference_found) THEN
                  CALL select_reference(reference_frame, velocity_frame%natoms, selected_index, labels, &
                                        reference_positions)
               ELSE
                  reference_positions = position_frame%value(:, selected_index)
               END IF
               CALL center_positions(reference_positions, masses)
            ELSE
               ALLOCATE (reference_positions(0, 0), current_positions(0, 0))
            END IF
         ELSE
            IF (velocity_frame%natoms /= SIZE(selected)) CALL fail("The atom count changes along the velocity trajectory")
            IF (ANY(velocity_frame%label(selected_index) /= labels)) &
               CALL fail("The selected atom labels or ordering change along the velocity trajectory")
            IF (position_found) CALL validate_position_frame(velocity_frame, position_frame)
         END IF

         IF (frames == capacity) THEN
            ALLOCATE (grown(n_dof, 2*capacity))
            grown(:, :capacity) = series
            CALL MOVE_ALLOC(grown, series)
            IF (mode_timeseries_requested) THEN
               ALLOCATE (grown_positions(n_dof, 2*capacity))
               grown_positions(:, :capacity) = position_series
               CALL MOVE_ALLOC(grown_positions, position_series)
            END IF
            capacity = 2*capacity
         END IF
         frames = frames + 1
         IF (position_found) THEN
            current_positions = position_frame%value(:, selected_index)
            CALL center_positions(current_positions, masses)
            CALL fit_rotation(current_positions, reference_positions, masses, rotation, ierr)
            IF (ierr /= 0) CALL fail("Rotational alignment did not converge")
         ELSE
            rotation = identity3()
         END IF
         item = 0
         DO atom = 1, n_selected
            velocity = MATMUL(rotation, velocity_frame%value(:, selected_index(atom)))*velocity_scale
            IF (mode_timeseries_requested) &
               aligned_position = MATMUL(rotation, current_positions(:, atom)) - reference_positions(:, atom)
            DO component = 1, 3
               item = item + 1
               series(item, frames) = SQRT(masses(atom))*velocity(component)
               IF (mode_timeseries_requested) &
                  position_series(item, frames) = SQRT(masses(atom))*aligned_position(component)
            END DO
         END DO
      END DO
      CALL velocities%close_file()
      IF (position_found) CALL positions%close_file()
      IF (frames < 2*n_correlation) &
         CALL fail("FRESEAN requires at least 2*correlation-frames selected trajectory frames")
      IF (frames < capacity) THEN
         ALLOCATE (grown(3*n_selected, frames))
         grown = series(:, :frames)
         CALL MOVE_ALLOC(grown, series)
         IF (mode_timeseries_requested) THEN
            ALLOCATE (grown_positions(3*n_selected, frames))
            grown_positions = position_series(:, :frames)
            CALL MOVE_ALLOC(grown_positions, position_series)
         END IF
      END IF
      IF (remove_mean) CALL remove_real_means(series)

      effective_dt_ps = dt_fs*REAL(stride, dp)/1000.0_dp
      CALL analyze_series(series, position_series, masses, labels, selection, effective_dt_ps, sigma_cm, &
                          constraints, frequency_cm, n_correlation, mode_count, position_found, remove_mean, &
                          output_path, mode_path, mode_requested, mode_spectrum_path, &
                          mode_spectrum_requested, mode_timeseries_path, mode_timeseries_requested)
   END SUBROUTINE run_fresean

   SUBROUTINE analyze_series(series, position_series, masses, labels, selection, dt_ps, sigma_cm, constraints, &
                             requested_frequency, n_correlation, requested_modes, aligned, mean_removed, output_path, &
                             mode_path, mode_requested, mode_spectrum_path, mode_spectrum_requested, &
                             mode_timeseries_path, mode_timeseries_requested)
      REAL(dp), INTENT(IN)                               :: series(:, :), position_series(:, :), &
                                                            masses(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: labels(:), selection
      REAL(dp), INTENT(IN)                               :: dt_ps, sigma_cm, constraints, &
                                                            requested_frequency
      INTEGER, INTENT(IN)                                :: n_correlation, requested_modes
      LOGICAL, INTENT(IN)                                :: aligned, mean_removed
      CHARACTER(LEN=*), INTENT(IN)                       :: output_path, mode_path
      LOGICAL, INTENT(IN)                                :: mode_requested
      CHARACTER(LEN=*), INTENT(IN)                       :: mode_spectrum_path
      LOGICAL, INTENT(IN)                                :: mode_spectrum_requested
      CHARACTER(LEN=*), INTENT(IN)                       :: mode_timeseries_path
      LOGICAL, INTENT(IN)                                :: mode_timeseries_requested

      CHARACTER(LEN=512)                                 :: message
      INTEGER :: frequency, ierr, mode, mode_count, mode_timeseries_unit, mode_unit, n_dof, &
         output_unit, selected_frequency, spectrum_unit
      REAL(dp)                                           :: captured, df_cm, dof, &
                                                            kinetic_temperature, normalization, &
                                                            selected_wavenumber, total
      REAL(dp), ALLOCATABLE :: eigenvalue_table(:, :), eigenvalues(:), eigenvectors(:, :), &
         matrices(:, :, :), mode_vectors(:, :), trace(:)

      n_dof = SIZE(series, 1)
      dof = REAL(n_dof, dp) - constraints
      mode_count = MIN(requested_modes, n_dof)
      kinetic_temperature = SUM(series**2)/REAL(SIZE(series, 2), dp)*mass_velocity_to_j_mol/ &
                            (gas_constant_j_mol_k*dof)
      CALL build_frequency_matrices(series, n_correlation, dt_ps, sigma_cm, matrices, df_cm)
      ALLOCATE (trace(n_correlation))
      DO frequency = 1, n_correlation
         trace(frequency) = matrix_trace(matrices(:, :, frequency))
      END DO
      normalization = SUM(trace)/dof
      IF (normalization <= TINY(1.0_dp)) CALL fail("The frequency-dependent velocity matrix has zero total intensity")
      matrices = matrices/normalization
      trace = trace/normalization

      IF (requested_frequency > REAL(n_correlation - 1, dp)*df_cm + 0.5_dp*df_cm) &
         CALL fail("--frequency-cm is outside the available frequency range")
      selected_frequency = NINT(requested_frequency/df_cm) + 1
      selected_frequency = MAX(1, MIN(n_correlation, selected_frequency))
      selected_wavenumber = REAL(selected_frequency - 1, dp)*df_cm
      ALLOCATE (eigenvalue_table(mode_count, n_correlation), mode_vectors(n_dof, mode_count))
      DO frequency = 1, n_correlation
         CALL diagonalize_symmetric(matrices(:, :, frequency), eigenvalues, eigenvectors, ierr)
         IF (ierr /= 0) CALL fail("A FRESEAN matrix eigendecomposition did not converge")
         eigenvalue_table(:, frequency) = eigenvalues(:mode_count)
         IF (frequency == selected_frequency) mode_vectors = eigenvectors(:, :mode_count)
      END DO

      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)") "# FRESEAN frequency-selective anharmonic vibrational analysis"
      WRITE (output_unit, "(A,I0)") "# frames: ", SIZE(series, 2)
      WRITE (output_unit, "(A,I0)") "# selected_atoms: ", SIZE(masses)
      WRITE (output_unit, "(A,F12.6)") "# effective_degrees_of_freedom: ", dof
      WRITE (output_unit, "(A,F16.8)") "# time_step_ps: ", dt_ps
      WRITE (output_unit, "(A,I0)") "# correlation_frames: ", n_correlation
      WRITE (output_unit, "(A,F16.8)") "# frequency_resolution_cm^-1: ", df_cm
      WRITE (output_unit, "(A,F16.8)") "# gaussian_sigma_cm^-1: ", sigma_cm
      WRITE (output_unit, "(A,F16.8)") "# kinetic_temperature_K: ", kinetic_temperature
      WRITE (output_unit, "(A,L1)") "# rotational_alignment: ", aligned
      WRITE (output_unit, "(A,L1)") "# Cartesian_means_removed: ", mean_removed
      WRITE (output_unit, "(A,A)") "# selection: ", TRIM(selection)
      WRITE (output_unit, "(A)", ADVANCE="no") &
         "# wavenumber[cm^-1] frequency[THz] total_VDoS retained_fraction"
      DO mode = 1, mode_count
         WRITE (output_unit, "(A,I0)", ADVANCE="no") " lambda_", mode
      END DO
      WRITE (output_unit, "(A)") ""
      DO frequency = 1, n_correlation
         total = trace(frequency)
         captured = 0.0_dp
         IF (total > TINY(1.0_dp)) captured = SUM(MAX(eigenvalue_table(:, frequency), 0.0_dp))/total
         WRITE (output_unit, "(4ES24.16)", ADVANCE="no") REAL(frequency - 1, dp)*df_cm, &
            REAL(frequency - 1, dp)*df_cm/thz_to_wavenumber, total, captured
         DO mode = 1, mode_count
            WRITE (output_unit, "(ES24.16)", ADVANCE="no") eigenvalue_table(mode, frequency)
         END DO
         WRITE (output_unit, "(A)") ""
      END DO
      CALL close_output(output_unit)

      IF (mode_requested) THEN
         CALL open_output(mode_path, mode_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL write_modes(mode_unit, labels, masses, selected_wavenumber, &
                          eigenvalue_table(:, selected_frequency), mode_vectors)
         CALL close_output(mode_unit)
      END IF
      IF (mode_spectrum_requested) THEN
         CALL open_output(mode_spectrum_path, spectrum_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL write_mode_spectra(spectrum_unit, df_cm, selected_wavenumber, trace, matrices, mode_vectors)
         CALL close_output(spectrum_unit)
      END IF
      IF (mode_timeseries_requested) THEN
         CALL open_output(mode_timeseries_path, mode_timeseries_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         CALL write_mode_timeseries(mode_timeseries_unit, dt_ps, selected_wavenumber, &
                                    position_series, series, mode_vectors)
         CALL close_output(mode_timeseries_unit)
      END IF
   END SUBROUTINE analyze_series

   SUBROUTINE build_frequency_matrices(series, n_correlation, dt_ps, sigma_cm, matrices, df_cm)
      REAL(dp), INTENT(IN)                               :: series(:, :)
      INTEGER, INTENT(IN)                                :: n_correlation
      REAL(dp), INTENT(IN)                               :: dt_ps, sigma_cm
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: matrices(:, :, :)
      REAL(dp), INTENT(OUT)                              :: df_cm

      COMPLEX(dp), ALLOCATABLE                           :: transformed(:, :), work_long(:), &
                                                            work_short(:)
      INTEGER                                            :: first, frames, frequency, lag, length, &
                                                            n_dof, second
      REAL(dp)                                           :: gaussian_norm
      REAL(dp), ALLOCATABLE                              :: window(:)

      n_dof = SIZE(series, 1)
      frames = SIZE(series, 2)
      length = 2*n_correlation - 1
      df_cm = thz_to_wavenumber/(REAL(length, dp)*dt_ps)
      ALLOCATE (transformed(frames, n_dof), work_long(frames), work_short(length), &
                window(length), matrices(n_dof, n_dof, n_correlation))
      DO first = 1, n_dof
         transformed(:, first) = CMPLX(series(first, :), 0.0_dp, KIND=dp)
         CALL fft_any_in_place(transformed(:, first), inverse=.FALSE.)
      END DO

      gaussian_norm = 1.0_dp/SQRT(2.0_dp*ACOS(-1.0_dp)*sigma_cm*sigma_cm)
      DO frequency = 0, n_correlation - 1
         window(frequency + 1) = gaussian_norm* &
                                 EXP(-0.5_dp*(REAL(frequency, dp)*df_cm/sigma_cm)**2)
      END DO
      window(n_correlation + 1:) = window(n_correlation:2:-1)
      work_short = CMPLX(window, 0.0_dp, KIND=dp)
      CALL fft_any_in_place(work_short, inverse=.TRUE.)
      window = REAL(work_short, dp)

      matrices = 0.0_dp
      DO first = 1, n_dof
         DO second = first, n_dof
            work_long = CMPLX(REAL(transformed(:, first)*CONJG(transformed(:, second)), dp), &
                              0.0_dp, KIND=dp)
            CALL fft_any_in_place(work_long, inverse=.TRUE.)
            work_short = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
            work_short(1) = CMPLX(window(1)*REAL(work_long(1), dp), 0.0_dp, KIND=dp)
            DO lag = 1, n_correlation - 1
               work_short(lag + 1) = CMPLX(window(lag + 1)*0.5_dp* &
                                           REAL(work_long(lag + 1) + work_long(frames - lag + 1), dp), 0.0_dp, KIND=dp)
            END DO
            work_short(n_correlation + 1:) = work_short(n_correlation:2:-1)
            CALL fft_any_in_place(work_short, inverse=.FALSE.)
            DO frequency = 1, n_correlation
               matrices(first, second, frequency) = REAL(work_short(frequency), dp)/REAL(frames, dp)
               matrices(second, first, frequency) = matrices(first, second, frequency)
            END DO
         END DO
      END DO
   END SUBROUTINE build_frequency_matrices

   PURE REAL(dp) FUNCTION matrix_trace(matrix)
      REAL(dp), INTENT(IN)                               :: matrix(:, :)

      INTEGER                                            :: index

      matrix_trace = 0.0_dp
      DO index = 1, MIN(SIZE(matrix, 1), SIZE(matrix, 2))
         matrix_trace = matrix_trace + matrix(index, index)
      END DO
   END FUNCTION matrix_trace

   SUBROUTINE write_modes(unit, labels, masses, frequency_cm, eigenvalues, vectors)
      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(LEN=*), INTENT(IN)                       :: labels(:)
      REAL(dp), INTENT(IN)                               :: masses(:), frequency_cm, eigenvalues(:), &
                                                            vectors(:, :)

      INTEGER                                            :: atom, mode
      REAL(dp)                                           :: displacement(3)

      DO mode = 1, SIZE(vectors, 2)
         WRITE (unit, "(I0)") SIZE(masses)
         WRITE (unit, "(A,F12.6,A,I0,A,ES16.8)") "FRESEAN frequency_cm^-1=", frequency_cm, &
            " mode=", mode, " eigenvalue=", eigenvalues(mode)
         DO atom = 1, SIZE(masses)
            displacement = vectors(3*atom - 2:3*atom, mode)/SQRT(masses(atom))
            WRITE (unit, "(A,3(1X,ES24.16))") TRIM(labels(atom)), displacement
         END DO
      END DO
   END SUBROUTINE write_modes

   SUBROUTINE write_mode_spectra(unit, df_cm, selected_frequency, trace, matrices, vectors)
      INTEGER, INTENT(IN)                                :: unit
      REAL(dp), INTENT(IN)                               :: df_cm, selected_frequency, trace(:), &
                                                            matrices(:, :, :), vectors(:, :)

      INTEGER                                            :: frequency, mode
      REAL(dp)                                           :: projection

      WRITE (unit, "(A,F16.8)") "# eigenvectors selected at wavenumber_cm^-1: ", selected_frequency
      WRITE (unit, "(A)", ADVANCE="no") "# wavenumber[cm^-1] frequency[THz] total_VDoS"
      DO mode = 1, SIZE(vectors, 2)
         WRITE (unit, "(A,I0)", ADVANCE="no") " mode_", mode
      END DO
      WRITE (unit, "(A)") ""
      DO frequency = 1, SIZE(trace)
         WRITE (unit, "(3ES24.16)", ADVANCE="no") REAL(frequency - 1, dp)*df_cm, &
            REAL(frequency - 1, dp)*df_cm/thz_to_wavenumber, trace(frequency)
         DO mode = 1, SIZE(vectors, 2)
            projection = DOT_PRODUCT(vectors(:, mode), MATMUL(matrices(:, :, frequency), vectors(:, mode)))
            WRITE (unit, "(ES24.16)", ADVANCE="no") projection
         END DO
         WRITE (unit, "(A)") ""
      END DO
   END SUBROUTINE write_mode_spectra

   SUBROUTINE write_mode_timeseries(unit, dt_ps, frequency_cm, positions, velocities, vectors)
      INTEGER, INTENT(IN)                                :: unit
      REAL(dp), INTENT(IN)                               :: dt_ps, frequency_cm, positions(:, :), &
                                                            velocities(:, :), vectors(:, :)

      INTEGER                                            :: frame, mode, mode_count
      REAL(dp)                                           :: energy_quantum, omega_ps, total_energy
      REAL(dp), ALLOCATABLE                              :: coordinate(:, :), energy(:, :), &
                                                            energy_fraction(:, :), phase_time(:), &
                                                            projected_velocity(:, :)

      mode_count = SIZE(vectors, 2)
      ALLOCATE (coordinate(mode_count, SIZE(positions, 2)), &
                projected_velocity(mode_count, SIZE(velocities, 2)), &
                energy(mode_count, SIZE(velocities, 2)), energy_fraction(mode_count, SIZE(velocities, 2)), &
                phase_time(mode_count))
      coordinate = MATMUL(TRANSPOSE(vectors), positions)
      projected_velocity = MATMUL(TRANSPOSE(vectors), velocities)
      DO mode = 1, mode_count
         coordinate(mode, :) = coordinate(mode, :) - &
                               SUM(coordinate(mode, :))/REAL(SIZE(coordinate, 2), dp)
      END DO
      omega_ps = 2.0_dp*ACOS(-1.0_dp)*speed_of_light_cm_per_ps*frequency_cm
      energy = 0.5_dp*mass_velocity_to_j_mol* &
               (projected_velocity**2 + (omega_ps*coordinate)**2)/1000.0_dp
      energy_quantum = wavenumber_to_kj_mol*frequency_cm
      energy_fraction = 0.0_dp
      DO frame = 1, SIZE(energy, 2)
         total_energy = SUM(energy(:, frame))
         IF (total_energy > TINY(1.0_dp)) energy_fraction(:, frame) = energy(:, frame)/total_energy
      END DO
      DO mode = 1, mode_count
         phase_time(mode) = phase_correlation_time(coordinate(mode, :), projected_velocity(mode, :), &
                                                   omega_ps, dt_ps)
      END DO

      WRITE (unit, "(A)") "# Time-resolved projection onto fixed FRESEAN modes"
      WRITE (unit, "(A,F16.8)") "# eigenvectors_selected_at_wavenumber_cm^-1: ", frequency_cm
      WRITE (unit, "(A,F16.8)") "# harmonic_energy_quantum_kJ_mol^-1: ", energy_quantum
      WRITE (unit, "(A)") "# mass_weighted_coordinate_unit: sqrt(g/mol)*angstrom"
      WRITE (unit, "(A)") "# mass_weighted_velocity_unit: sqrt(g/mol)*angstrom/ps"
      WRITE (unit, "(A)") "# projected_coordinate_means_removed: T"
      WRITE (unit, "(A)") &
         "# Energies are harmonic proxies; fractions are normalized over the extracted modes only."
      WRITE (unit, "(A)") &
         "# Phase-correlation times describe the continuously sampled trajectory, not an intrinsic IVR lifetime."
      WRITE (unit, "(A)") "# A phase-correlation time of -1 means no 1/e crossing within half the trajectory."
      DO mode = 1, mode_count
         WRITE (unit, "(A,I0,A,ES24.16)") "# mode_", mode, "_phase_correlation_1e_ps: ", phase_time(mode)
      END DO
      WRITE (unit, "(A)", ADVANCE="no") "# time[ps]"
      DO mode = 1, mode_count
         WRITE (unit, "(A,I0,A,I0,A,I0,A,I0,A,I0)", ADVANCE="no") &
            " q_", mode, " velocity_", mode, " energy_kJ_mol_", mode, " quanta_", mode, " fraction_", mode
      END DO
      WRITE (unit, "(A)") ""
      DO frame = 1, SIZE(energy, 2)
         WRITE (unit, "(ES24.16)", ADVANCE="no") REAL(frame - 1, dp)*dt_ps
         DO mode = 1, mode_count
            WRITE (unit, "(5ES24.16)", ADVANCE="no") coordinate(mode, frame), projected_velocity(mode, frame), &
               energy(mode, frame), energy(mode, frame)/energy_quantum, energy_fraction(mode, frame)
         END DO
         WRITE (unit, "(A)") ""
      END DO
   END SUBROUTINE write_mode_timeseries

   REAL(dp) FUNCTION phase_correlation_time(coordinate, velocity, omega_ps, dt_ps)
      REAL(dp), INTENT(IN)                               :: coordinate(:), velocity(:), omega_ps, &
                                                            dt_ps

      COMPLEX(dp), ALLOCATABLE                           :: correlation(:), phase(:, :)
      INTEGER                                            :: frame, lag, maximum_lag
      REAL(dp)                                           :: amplitude, amplitude_cutoff, &
                                                            maximum_amplitude

      maximum_lag = SIZE(coordinate)/2
      ALLOCATE (phase(1, SIZE(coordinate)))
      maximum_amplitude = MAXVAL(SQRT((omega_ps*coordinate)**2 + velocity**2))
      amplitude_cutoff = SQRT(EPSILON(1.0_dp))*maximum_amplitude
      phase = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      DO frame = 1, SIZE(coordinate)
         amplitude = SQRT((omega_ps*coordinate(frame))**2 + velocity(frame)**2)
         IF (amplitude > amplitude_cutoff) &
            phase(1, frame) = CMPLX(omega_ps*coordinate(frame), velocity(frame), KIND=dp)/amplitude
      END DO
      CALL complex_autocorrelation_sum(phase, maximum_lag, correlation)
      DO lag = 0, maximum_lag
         correlation(lag) = correlation(lag)/REAL(SIZE(coordinate) - lag, dp)
      END DO
      phase_correlation_time = -1.0_dp
      IF (ABS(correlation(0)) <= TINY(1.0_dp)) RETURN
      DO lag = 1, maximum_lag
         IF (ABS(correlation(lag)/correlation(0)) <= EXP(-1.0_dp)) THEN
            phase_correlation_time = REAL(lag, dp)*dt_ps
            RETURN
         END IF
      END DO
   END FUNCTION phase_correlation_time

   SUBROUTINE validate_position_frame(velocity_frame, position_frame)
      TYPE(frame_type), INTENT(IN)                       :: velocity_frame, position_frame

      IF (position_frame%natoms /= velocity_frame%natoms) &
         CALL fail("Position and velocity trajectories have different atom counts")
      IF (ANY(position_frame%label /= velocity_frame%label)) &
         CALL fail("Position and velocity trajectories have different atom labels or ordering")
   END SUBROUTINE validate_position_frame

   SUBROUTINE select_reference(reference_frame, full_atoms, selected_index, labels, positions)
      TYPE(frame_type), INTENT(IN)                       :: reference_frame
      INTEGER, INTENT(IN)                                :: full_atoms, selected_index(:)
      CHARACTER(LEN=*), INTENT(IN)                       :: labels(:)
      REAL(dp), INTENT(OUT)                              :: positions(:, :)

      IF (reference_frame%natoms == full_atoms) THEN
         IF (ANY(reference_frame%label(selected_index) /= labels)) &
            CALL fail("The reference atom labels do not match the selection")
         positions = reference_frame%value(:, selected_index)
      ELSE IF (reference_frame%natoms == SIZE(selected_index)) THEN
         IF (ANY(reference_frame%label /= labels)) &
            CALL fail("The reference atom labels do not match the selection")
         positions = reference_frame%value
      ELSE
         CALL fail("The reference must contain either all trajectory atoms or exactly the selected atoms")
      END IF
   END SUBROUTINE select_reference

   SUBROUTINE print_fresean_help()
      WRITE (*, "(A)") "Usage: fresean.x --velocity PROJECT-vel-1.xyz --mass-file masses.dat [options]"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Required:"
      WRITE (*, "(A)") "  --velocity FILE            CP2K XYZ velocity trajectory"
      WRITE (*, "(A)") "  --mass-file FILE           LABEL MASS[g/mol] table"
      WRITE (*, "(A)") "  --dt-fs VALUE              time between input frames"
      WRITE (*, "(A)") "  --correlation-frames N     correlation points (default 500)"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Frequency-selective analysis:"
      WRITE (*, "(A)") "  --sigma-cm VALUE           Gaussian spectral width (default 10 cm^-1)"
      WRITE (*, "(A)") "  --frequency-cm VALUE       frequency for extracted modes (default zero)"
      WRITE (*, "(A)") "  --mode-count N             leading modes reported at each frequency (default 10)"
      WRITE (*, "(A)") "  --output FILE              VDoS and leading eigenvalues (default stdout)"
      WRITE (*, "(A)") "  --mode-file FILE           mass-unweighted displacement vectors in XYZ format"
      WRITE (*, "(A)") "  --mode-spectrum FILE       projected VDoS of the extracted modes"
      WRITE (*, "(A)") "  --mode-timeseries FILE     time-resolved fixed-mode coordinates and energy proxies"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Alignment and conventions:"
      WRITE (*, "(A)") "  --position FILE            synchronized positions; enables rotational alignment"
      WRITE (*, "(A)") "  --reference FILE           alignment reference (default first selected frame)"
      WRITE (*, "(A)") "  --select SPEC              all, name:..., or index:..."
      WRITE (*, "(A)") "  --constraints N            constrained selected degrees of freedom"
      WRITE (*, "(A)") "  --velocity-unit UNIT       cp2k (default), angstrom/ps, or angstrom/fs"
      WRITE (*, "(A)") "  --remove-mean              remove each Cartesian time-series mean"
      WRITE (*, "(A)") "  --first/--last/--stride N  frame selection"
   END SUBROUTINE print_fresean_help

END MODULE trajana_fresean_analysis
