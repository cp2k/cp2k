MODULE trajana_dsf_analysis
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option,&
                                              get_real_option,&
                                              has_flag
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_kinds,                   ONLY: dp
   USE trajana_selections,              ONLY: build_selection
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
      CHARACTER(LEN=:), ALLOCATABLE                      :: input_path, output_path, q_path, &
                                                            selection, spectrum_path, window
      COMPLEX(dp), ALLOCATABLE                           :: correlation(:), density(:, :), &
                                                            grown(:, :)
      INTEGER :: atom, capacity, expected_atoms, first, frames, frequency, ierr, iq, lag, last, &
         maximum_lag, nfft, nq, output_unit, selected_atoms, spectrum_unit, stride
      LOGICAL                                            :: eof, found, remove_mean, &
                                                            spectrum_requested
      LOGICAL, ALLOCATABLE                               :: selected(:)
      REAL(dp)                                           :: dt_fs, effective_dt, frequency_hz, &
                                                            phase, spectral_value, wavenumber, &
                                                            window_norm
      REAL(dp), ALLOCATABLE                              :: power(:), q(:, :)
      TYPE(frame_type)                                   :: frame
      TYPE(xyz_reader_type)                              :: trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--spectrum", spectrum_path, spectrum_requested)
      CALL get_option("--q-file", q_path, found)
      IF (.NOT. found) CALL fail("dsf requires --q-file FILE")
      CALL get_option("--select", selection, found, "all")
      CALL get_option("--window", window, found, "none")
      window = lower_case(window)
      IF (window /= "none" .AND. window /= "hann") CALL fail("--window expects none or hann")
      CALL get_real_option("--dt-fs", dt_fs, -1.0_dp)
      CALL get_integer_option("--max-lag", maximum_lag, -1)
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      remove_mean = has_flag("--remove-mean")
      IF (dt_fs <= 0.0_dp) CALL fail("dsf requires a positive --dt-fs")
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")

      CALL read_q_vectors(q_path, q, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      nq = SIZE(q, 2)
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
               ALLOCATE (density(nq, capacity))
            ELSE
               ALLOCATE (grown(nq, 2*capacity))
               grown(:, :capacity) = density
               CALL MOVE_ALLOC(grown, density)
               capacity = 2*capacity
            END IF
         END IF
         frames = frames + 1
         density(:, frames) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO iq = 1, nq
            DO atom = 1, frame%natoms
               IF (.NOT. selected(atom)) CYCLE
               phase = DOT_PRODUCT(q(:, iq), frame%value(:, atom))
               density(iq, frames) = density(iq, frames) + CMPLX(COS(phase), SIN(phase), KIND=dp)
            END DO
         END DO
      END DO
      CALL trajectory%close_file()
      IF (frames < 2) CALL fail("dsf requires at least two selected frames")
      IF (frames < capacity) THEN
         ALLOCATE (grown(nq, frames))
         grown = density(:, :frames)
         CALL MOVE_ALLOC(grown, density)
      END IF
      IF (remove_mean) CALL remove_complex_means(density)
      IF (maximum_lag < 0) maximum_lag = frames - 1
      IF (maximum_lag >= frames) CALL fail("--max-lag must be smaller than the number of selected frames")

      CALL complex_autocorrelation_sum(density, maximum_lag, correlation)
      DO lag = 0, maximum_lag
         correlation(lag) = correlation(lag)/(REAL(nq, dp)*REAL(selected_atoms, dp)*REAL(frames - lag, dp))
      END DO
      IF (ABS(correlation(0)) <= TINY(1.0_dp)) CALL fail("Intermediate scattering function at zero lag is zero")

      effective_dt = dt_fs*REAL(stride, dp)
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)") "# lag [fs]   Re(F(q,t))   Im(F(q,t))   Re(F(q,t)/F(q,0))"
      WRITE (output_unit, "(A,A)") "# selection: ", TRIM(selection)
      DO lag = 0, maximum_lag
         WRITE (output_unit, "(4ES24.16)") REAL(lag, dp)*effective_dt, REAL(correlation(lag), dp), &
            AIMAG(correlation(lag)), REAL(correlation(lag)/correlation(0), dp)
      END DO
      CALL close_output(output_unit)

      IF (spectrum_requested) THEN
         CALL complex_power_sum(density, window, power, nfft, window_norm)
         CALL open_output(spectrum_path, spectrum_unit, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         WRITE (spectrum_unit, "(A)") "# wavenumber [cm^-1]   S(q,omega) [fs]"
         WRITE (spectrum_unit, "(A,A)") "# window: ", TRIM(window)
         DO frequency = 0, nfft/2
            frequency_hz = REAL(frequency, dp)/(REAL(nfft, dp)*effective_dt*1.0e-15_dp)
            wavenumber = frequency_hz/2.99792458e10_dp
            spectral_value = power(frequency)*effective_dt/(REAL(nq, dp)*REAL(selected_atoms, dp)*window_norm)
            WRITE (spectrum_unit, "(2ES24.16)") wavenumber, spectral_value
         END DO
         CALL close_output(spectrum_unit)
      END IF
   END SUBROUTINE run_dsf

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
         IF (ios /= 0 .OR. SQRT(DOT_PRODUCT(vector, vector)) <= TINY(1.0_dp)) THEN
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
