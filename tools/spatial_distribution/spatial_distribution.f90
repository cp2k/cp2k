!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! Generalized reimplementation of the water-specific rot_anal_H2O.c workflow
! originally written by Thomas D. Kuehne.

MODULE spatial_distribution_types
   IMPLICIT NONE

   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)
   INTEGER, PARAMETER :: path_length = 512
   INTEGER, PARAMETER :: label_length = 64
   REAL(KIND=dp), PARAMETER :: angstrom_to_bohr = 1.0_dp/0.529177210903_dp

   TYPE reference_type
      INTEGER :: origin = 0
      INTEGER :: arm1 = 0
      INTEGER :: arm2 = 0
   END TYPE reference_type

   TYPE target_type
      CHARACTER(LEN=label_length) :: name = ""
      CHARACTER(LEN=label_length) :: selector = ""
      REAL(KIND=dp), ALLOCATABLE :: histogram(:, :, :)
      REAL(KIND=dp) :: normalization_weight = 0.0_dp
      INTEGER :: observations = 0
   END TYPE target_type

   TYPE config_type
      CHARACTER(LEN=path_length) :: trajectory = ""
      CHARACTER(LEN=path_length) :: cell_file = ""
      CHARACTER(LEN=path_length) :: output_prefix = "sdf"
      CHARACTER(LEN=16) :: frame_mode = "BISECTOR"
      CHARACTER(LEN=8) :: periodic = "XYZ"
      LOGICAL :: periodic_axes(3) = [.TRUE., .TRUE., .TRUE.]
      LOGICAL :: exclude_reference = .TRUE.
      LOGICAL :: static_cell_set = .FALSE.
      REAL(KIND=dp) :: static_cell(3, 3) = 0.0_dp
      INTEGER :: grid(3) = [81, 81, 81]
      REAL(KIND=dp) :: lower(3) = [-5.0_dp, -5.0_dp, -5.0_dp]
      REAL(KIND=dp) :: upper(3) = [5.0_dp, 5.0_dp, 5.0_dp]
      REAL(KIND=dp) :: cutoff = 0.0_dp
      INTEGER :: first_frame = 1
      INTEGER :: last_frame = HUGE(1)
      INTEGER :: stride = 1
      TYPE(reference_type), ALLOCATABLE :: references(:)
      TYPE(target_type), ALLOCATABLE :: targets(:)
   END TYPE config_type

CONTAINS

   SUBROUTINE fail(message)
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      WRITE (*, '(A)') "Error: "//TRIM(message)
      ERROR STOP 1
   END SUBROUTINE fail

   PURE FUNCTION uppercase(text) RESULT(RESULT)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      CHARACTER(LEN=LEN(text))                           :: result

      INTEGER                                            :: code, i

      RESULT = text
      DO i = 1, LEN(text)
         code = IACHAR(RESULT(i:i))
         IF (code >= IACHAR("a") .AND. code <= IACHAR("z")) RESULT(i:i) = ACHAR(code - 32)
      END DO
   END FUNCTION uppercase

   PURE FUNCTION strip_comment(input) RESULT(output)
      CHARACTER(LEN=*), INTENT(IN)                       :: input
      CHARACTER(LEN=LEN(input))                          :: output

      INTEGER                                            :: exclamation, hash, marker

      output = input
      exclamation = INDEX(output, "!")
      hash = INDEX(output, "#")
      IF (exclamation == 0) THEN
         marker = hash
      ELSE IF (hash == 0) THEN
         marker = exclamation
      ELSE
         marker = MIN(exclamation, hash)
      END IF
      IF (marker > 0) output(marker:) = ""
      output = ADJUSTL(output)
   END FUNCTION strip_comment

   FUNCTION resolve_path(config_file, value) RESULT(path)
      CHARACTER(LEN=*), INTENT(IN)                       :: config_file, value
      CHARACTER(LEN=path_length)                         :: path

      INTEGER                                            :: slash

      path = TRIM(value)
      IF (LEN_TRIM(value) == 0 .OR. value(1:1) == "/") RETURN
      slash = SCAN(TRIM(config_file), "/", BACK=.TRUE.)
      IF (slash > 0) path = config_file(:slash)//TRIM(value)
   END FUNCTION resolve_path

END MODULE spatial_distribution_types

MODULE spatial_distribution_input
   USE spatial_distribution_types, ONLY: config_type, dp, fail, path_length, reference_type, &
                                         resolve_path, strip_comment, uppercase
   IMPLICIT NONE

CONTAINS

   SUBROUTINE read_config(filename, config)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(config_type), INTENT(OUT)                     :: config

      CHARACTER(LEN=1024)                                :: line
      CHARACTER(LEN=64)                                  :: keyword
      INTEGER :: input_unit, ios, nreferences, ntargets, pattern_arm1, pattern_arm1_stride, &
         pattern_arm2, pattern_arm2_stride, pattern_count, pattern_origin, pattern_origin_stride

      nreferences = 0
      ntargets = 0
      OPEN (NEWUNIT=input_unit, FILE=TRIM(filename), STATUS="old", ACTION="read", IOSTAT=ios)
      IF (ios /= 0) CALL fail("Could not open configuration file "//TRIM(filename))

      DO
         READ (input_unit, '(A)', IOSTAT=ios) line
         IF (ios < 0) EXIT
         IF (ios > 0) CALL fail("Could not read configuration file "//TRIM(filename))
         line = strip_comment(line)
         IF (LEN_TRIM(line) == 0) CYCLE
         READ (line, *, IOSTAT=ios) keyword
         IF (ios /= 0) CALL fail("Malformed configuration line: "//TRIM(line))
         keyword = uppercase(keyword)
         SELECT CASE (TRIM(keyword))
         CASE ("REFERENCE")
            nreferences = nreferences + 1
         CASE ("REFERENCE_PATTERN")
            READ (line, *, IOSTAT=ios) keyword, pattern_origin, pattern_arm1, pattern_arm2, &
               pattern_origin_stride, pattern_arm1_stride, pattern_arm2_stride, pattern_count
            IF (ios /= 0 .OR. pattern_count < 1) CALL fail("Malformed REFERENCE_PATTERN: "//TRIM(line))
            nreferences = nreferences + pattern_count
         CASE ("TARGET")
            ntargets = ntargets + 1
         END SELECT
      END DO

      IF (nreferences == 0) CALL fail("At least one REFERENCE or REFERENCE_PATTERN is required")
      IF (ntargets == 0) CALL fail("At least one TARGET is required")
      ALLOCATE (config%references(nreferences), config%targets(ntargets))

      REWIND (input_unit)
      CALL parse_config(input_unit, config)
      CLOSE (input_unit)
      CALL validate_config(filename, config)
   END SUBROUTINE read_config

   SUBROUTINE parse_config(input_unit, config)
      INTEGER, INTENT(IN)                                :: input_unit
      TYPE(config_type), INTENT(INOUT)                   :: config

      CHARACTER(LEN=1024)                                :: line
      CHARACTER(LEN=64)                                  :: keyword
      CHARACTER(LEN=path_length)                         :: path_value
      INTEGER                                            :: arm1, arm1_stride, arm2, arm2_stride, &
                                                            count, i, ios, origin, origin_stride, &
                                                            ref_index, target_index
      REAL(KIND=dp)                                      :: values(9)

      ref_index = 0
      target_index = 0
      DO
         READ (input_unit, '(A)', IOSTAT=ios) line
         IF (ios < 0) EXIT
         IF (ios > 0) CALL fail("Could not read configuration")
         line = strip_comment(line)
         IF (LEN_TRIM(line) == 0) CYCLE
         READ (line, *, IOSTAT=ios) keyword
         IF (ios /= 0) CALL fail("Malformed configuration line: "//TRIM(line))
         keyword = uppercase(keyword)

         SELECT CASE (TRIM(keyword))
         CASE ("TRAJECTORY")
            CALL parse_path_value(line, path_value, ios)
            IF (ios == 0) config%trajectory = path_value
         CASE ("CELL_FILE")
            CALL parse_path_value(line, path_value, ios)
            IF (ios == 0) config%cell_file = path_value
         CASE ("CELL")
            READ (line, *, IOSTAT=ios) keyword, values
            IF (ios == 0) THEN
               config%static_cell(:, 1) = values(1:3)
               config%static_cell(:, 2) = values(4:6)
               config%static_cell(:, 3) = values(7:9)
               config%static_cell_set = .TRUE.
            END IF
         CASE ("FRAME_MODE")
            READ (line, *, IOSTAT=ios) keyword, config%frame_mode
            config%frame_mode = uppercase(config%frame_mode)
         CASE ("PERIODIC")
            READ (line, *, IOSTAT=ios) keyword, config%periodic
            config%periodic = uppercase(config%periodic)
         CASE ("EXCLUDE_REFERENCE")
            READ (line, *, IOSTAT=ios) keyword, config%exclude_reference
         CASE ("REFERENCE")
            ref_index = ref_index + 1
            READ (line, *, IOSTAT=ios) keyword, config%references(ref_index)%origin, &
               config%references(ref_index)%arm1, config%references(ref_index)%arm2
         CASE ("REFERENCE_PATTERN")
            READ (line, *, IOSTAT=ios) keyword, origin, arm1, arm2, origin_stride, arm1_stride, &
               arm2_stride, count
            IF (ios == 0) THEN
               DO i = 0, count - 1
                  ref_index = ref_index + 1
                  config%references(ref_index) = reference_type(origin + i*origin_stride, &
                                                                arm1 + i*arm1_stride, &
                                                                arm2 + i*arm2_stride)
               END DO
            END IF
         CASE ("TARGET")
            target_index = target_index + 1
            READ (line, *, IOSTAT=ios) keyword, config%targets(target_index)%name, &
               config%targets(target_index)%selector
         CASE ("GRID")
            READ (line, *, IOSTAT=ios) keyword, config%grid
         CASE ("BOUNDS")
            READ (line, *, IOSTAT=ios) keyword, config%lower(1), config%upper(1), &
               config%lower(2), config%upper(2), config%lower(3), config%upper(3)
         CASE ("CUTOFF")
            READ (line, *, IOSTAT=ios) keyword, config%cutoff
         CASE ("FIRST_FRAME")
            READ (line, *, IOSTAT=ios) keyword, config%first_frame
         CASE ("LAST_FRAME")
            READ (line, *, IOSTAT=ios) keyword, config%last_frame
            IF (config%last_frame == 0) config%last_frame = HUGE(1)
         CASE ("STRIDE")
            READ (line, *, IOSTAT=ios) keyword, config%stride
         CASE ("OUTPUT_PREFIX")
            CALL parse_path_value(line, path_value, ios)
            IF (ios == 0) config%output_prefix = path_value
         CASE DEFAULT
            CALL fail("Unknown configuration keyword "//TRIM(keyword))
         END SELECT
         IF (ios /= 0) CALL fail("Malformed configuration line: "//TRIM(line))
      END DO
   END SUBROUTINE parse_config

   SUBROUTINE parse_path_value(line, value, ios)
      CHARACTER(LEN=*), INTENT(IN)                       :: line
      CHARACTER(LEN=*), INTENT(OUT)                      :: value
      INTEGER, INTENT(OUT)                               :: ios

      CHARACTER(LEN=1)                                   :: quote
      INTEGER                                            :: separator, value_length

      value = ""
      separator = SCAN(TRIM(line), " "//ACHAR(9))
      IF (separator == 0) THEN
         ios = 1
         RETURN
      END IF
      value = ADJUSTL(line(separator + 1:))
      value_length = LEN_TRIM(value)
      IF (value_length == 0) THEN
         ios = 1
         RETURN
      END IF
      IF (value(1:1) == '"' .OR. value(1:1) == "'") THEN
         quote = value(1:1)
         IF (value_length <= 2 .OR. value(value_length:value_length) /= quote) THEN
            ios = 1
            RETURN
         END IF
         value = value(2:value_length - 1)
      END IF
      ios = 0
   END SUBROUTINE parse_path_value

   SUBROUTINE validate_config(filename, config)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(config_type), INTENT(INOUT)                   :: config

      INTEGER                                            :: i, j

      IF (LEN_TRIM(config%trajectory) == 0) CALL fail("TRAJECTORY is required")
      IF (config%static_cell_set .AND. LEN_TRIM(config%cell_file) > 0) &
         CALL fail("CELL and CELL_FILE are mutually exclusive")
      IF (ANY(config%grid < 1)) CALL fail("All GRID dimensions must be positive")
      IF (ANY(config%upper <= config%lower)) CALL fail("Each upper BOUNDS value must exceed its lower value")
      IF (config%cutoff < 0.0_dp) CALL fail("CUTOFF must not be negative")
      IF (config%first_frame < 1 .OR. config%stride < 1 .OR. config%last_frame < config%first_frame) &
         CALL fail("Invalid FIRST_FRAME, LAST_FRAME, or STRIDE")
      IF (TRIM(config%frame_mode) /= "BISECTOR" .AND. TRIM(config%frame_mode) /= "TRIAD") &
         CALL fail("FRAME_MODE must be BISECTOR or TRIAD")

      config%periodic_axes = .FALSE.
      IF (TRIM(config%periodic) /= "NONE") THEN
         config%periodic_axes(1) = INDEX(config%periodic, "X") > 0
         config%periodic_axes(2) = INDEX(config%periodic, "Y") > 0
         config%periodic_axes(3) = INDEX(config%periodic, "Z") > 0
         IF (.NOT. ANY(config%periodic_axes)) CALL fail("Invalid PERIODIC setting")
      END IF

      DO i = 1, SIZE(config%references)
         IF (MIN(config%references(i)%origin, config%references(i)%arm1, config%references(i)%arm2) < 1) &
            CALL fail("REFERENCE indices are one-based and must be positive")
         IF (config%references(i)%origin == config%references(i)%arm1 .OR. &
             config%references(i)%origin == config%references(i)%arm2 .OR. &
             config%references(i)%arm1 == config%references(i)%arm2) &
            CALL fail("The three indices of a REFERENCE must be distinct")
      END DO
      DO i = 1, SIZE(config%targets)
         IF (LEN_TRIM(config%targets(i)%name) == 0 .OR. LEN_TRIM(config%targets(i)%selector) == 0) &
            CALL fail("TARGET needs a name and a label selector")
         DO j = 1, i - 1
            IF (uppercase(TRIM(config%targets(i)%name)) == uppercase(TRIM(config%targets(j)%name))) &
               CALL fail("TARGET names must be unique")
         END DO
         ALLOCATE (config%targets(i)%histogram(config%grid(1), config%grid(2), config%grid(3)))
         config%targets(i)%histogram = 0.0_dp
      END DO

      config%trajectory = resolve_path(filename, config%trajectory)
      IF (LEN_TRIM(config%cell_file) > 0) config%cell_file = resolve_path(filename, config%cell_file)
      config%output_prefix = resolve_path(filename, config%output_prefix)
   END SUBROUTINE validate_config

END MODULE spatial_distribution_input

MODULE spatial_distribution_geometry
   USE spatial_distribution_types, ONLY: config_type, dp, fail, reference_type
   IMPLICIT NONE

CONTAINS

   PURE FUNCTION cross_product(a, b) RESULT(c)
      REAL(KIND=dp), INTENT(IN)                          :: a(3), b(3)
      REAL(KIND=dp)                                      :: c(3)

      c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
   END FUNCTION cross_product

   SUBROUTINE normalize(vector, description)
      REAL(KIND=dp), INTENT(INOUT)                       :: vector(3)
      CHARACTER(LEN=*), INTENT(IN)                       :: description

      REAL(KIND=dp)                                      :: norm

      norm = SQRT(DOT_PRODUCT(vector, vector))
      IF (norm <= 100.0_dp*EPSILON(norm)) CALL fail("Cannot define local frame: "//TRIM(description))
      vector = vector/norm
   END SUBROUTINE normalize

   PURE FUNCTION determinant(matrix) RESULT(det)
      REAL(KIND=dp), INTENT(IN)                          :: matrix(3, 3)
      REAL(KIND=dp)                                      :: det

      det = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) &
            - matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) &
            + matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
   END FUNCTION determinant

   SUBROUTINE inverse_matrix(matrix, inverse, det)
      REAL(KIND=dp), INTENT(IN)                          :: matrix(3, 3)
      REAL(KIND=dp), INTENT(OUT)                         :: inverse(3, 3), det

      det = determinant(matrix)
      IF (ABS(det) <= 100.0_dp*EPSILON(det)) CALL fail("Simulation cell is singular")

      inverse(1, 1) = matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)
      inverse(1, 2) = matrix(1, 3)*matrix(3, 2) - matrix(1, 2)*matrix(3, 3)
      inverse(1, 3) = matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2)
      inverse(2, 1) = matrix(2, 3)*matrix(3, 1) - matrix(2, 1)*matrix(3, 3)
      inverse(2, 2) = matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1)
      inverse(2, 3) = matrix(1, 3)*matrix(2, 1) - matrix(1, 1)*matrix(2, 3)
      inverse(3, 1) = matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1)
      inverse(3, 2) = matrix(1, 2)*matrix(3, 1) - matrix(1, 1)*matrix(3, 2)
      inverse(3, 3) = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
      inverse = inverse/det
   END SUBROUTINE inverse_matrix

   PURE FUNCTION minimum_image(displacement, cell, inverse, periodic_axes) RESULT(wrapped)
      REAL(KIND=dp), INTENT(IN)                          :: displacement(3), cell(3, 3), &
                                                            inverse(3, 3)
      LOGICAL, INTENT(IN)                                :: periodic_axes(3)
      REAL(KIND=dp)                                      :: wrapped(3)

      INTEGER                                            :: i, lower(3), n1, n2, n3, shift(3), &
                                                            upper(3)
      REAL(KIND=dp)                                      :: bound, candidate(3), candidate_norm, &
                                                            fractional(3), wrapped_norm

      fractional = MATMUL(inverse, displacement)
      DO i = 1, 3
         IF (periodic_axes(i)) THEN
            shift(i) = NINT(fractional(i))
         ELSE
            shift(i) = 0
         END IF
      END DO
      wrapped = MATMUL(cell, fractional - REAL(shift, KIND=dp))
      wrapped_norm = DOT_PRODUCT(wrapped, wrapped)

      ! Component-wise rounding is not always the nearest image in a skewed cell.
      ! Reciprocal-vector bounds make the following small lattice search exact.
      DO i = 1, 3
         IF (periodic_axes(i)) THEN
            bound = SQRT(wrapped_norm*DOT_PRODUCT(inverse(i, :), inverse(i, :))) + &
                    10.0_dp*EPSILON(1.0_dp)
            lower(i) = CEILING(fractional(i) - bound)
            upper(i) = FLOOR(fractional(i) + bound)
         ELSE
            lower(i) = 0
            upper(i) = 0
         END IF
      END DO
      DO n1 = lower(1), upper(1)
         DO n2 = lower(2), upper(2)
            DO n3 = lower(3), upper(3)
               candidate = MATMUL(cell, fractional - REAL([n1, n2, n3], KIND=dp))
               candidate_norm = DOT_PRODUCT(candidate, candidate)
               IF (candidate_norm < wrapped_norm) THEN
                  wrapped = candidate
                  wrapped_norm = candidate_norm
               END IF
            END DO
         END DO
      END DO
   END FUNCTION minimum_image

   SUBROUTINE local_basis(reference, coordinates, cell, inverse, config, basis, arms)
      TYPE(reference_type), INTENT(IN)                   :: reference
      REAL(KIND=dp), INTENT(IN)                          :: coordinates(:, :), cell(3, 3), &
                                                            inverse(3, 3)
      TYPE(config_type), INTENT(IN)                      :: config
      REAL(KIND=dp), INTENT(OUT)                         :: basis(3, 3), arms(3, 2)

      REAL(KIND=dp)                                      :: arm1(3), arm2(3), ex(3), ey(3), ez(3), &
                                                            unit1(3), unit2(3)

      arm1 = minimum_image(coordinates(:, reference%arm1) - coordinates(:, reference%origin), &
                           cell, inverse, config%periodic_axes)
      arm2 = minimum_image(coordinates(:, reference%arm2) - coordinates(:, reference%origin), &
                           cell, inverse, config%periodic_axes)

      IF (TRIM(config%frame_mode) == "BISECTOR") THEN
         unit1 = arm1
         unit2 = arm2
         CALL normalize(unit1, "first orientation arm has zero length")
         CALL normalize(unit2, "second orientation arm has zero length")
         ez = unit1 + unit2
         ex = unit2 - unit1
         CALL normalize(ez, "BISECTOR arms are collinear and point in opposite directions")
         ex = ex - DOT_PRODUCT(ex, ez)*ez
         CALL normalize(ex, "BISECTOR arms do not define an x axis")
      ELSE
         ez = arm1
         CALL normalize(ez, "first TRIAD arm has zero length")
         ex = arm2 - DOT_PRODUCT(arm2, ez)*ez
         CALL normalize(ex, "TRIAD arms are collinear")
      END IF

      ey = cross_product(ez, ex)
      CALL normalize(ey, "orientation vectors do not define a plane")
      ex = cross_product(ey, ez)
      basis(:, 1) = ex
      basis(:, 2) = ey
      basis(:, 3) = ez
      arms(:, 1) = MATMUL(TRANSPOSE(basis), arm1)
      arms(:, 2) = MATMUL(TRANSPOSE(basis), arm2)
   END SUBROUTINE local_basis

END MODULE spatial_distribution_geometry

MODULE spatial_distribution_io
   USE spatial_distribution_types, ONLY: angstrom_to_bohr, config_type, dp, fail, label_length, &
                                         path_length, uppercase
   IMPLICIT NONE

   TYPE cell_series_type
      INTEGER, ALLOCATABLE :: steps(:)
      REAL(KIND=dp), ALLOCATABLE :: cells(:, :, :)
   END TYPE cell_series_type

CONTAINS

   SUBROUTINE read_cell_series(filename, series)
      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      TYPE(cell_series_type), INTENT(OUT)                :: series

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: cell_unit, count, ios, step
      REAL(KIND=dp)                                      :: time, values(9), volume

      count = 0
      OPEN (NEWUNIT=cell_unit, FILE=TRIM(filename), STATUS="old", ACTION="read", IOSTAT=ios)
      IF (ios /= 0) CALL fail("Could not open CELL_FILE "//TRIM(filename))
      DO
         READ (cell_unit, '(A)', IOSTAT=ios) line
         IF (ios < 0) EXIT
         IF (ios > 0) CALL fail("Could not read CELL_FILE "//TRIM(filename))
         READ (line, *, IOSTAT=ios) step, time, values, volume
         IF (ios == 0) count = count + 1
      END DO
      IF (count == 0) CALL fail("CELL_FILE contains no CP2K cell records")

      ALLOCATE (series%steps(count), series%cells(3, 3, count))
      REWIND (cell_unit)
      count = 0
      DO
         READ (cell_unit, '(A)', IOSTAT=ios) line
         IF (ios < 0) EXIT
         READ (line, *, IOSTAT=ios) step, time, values, volume
         IF (ios /= 0) CYCLE
         count = count + 1
         series%steps(count) = step
         series%cells(:, 1, count) = values(1:3)
         series%cells(:, 2, count) = values(4:6)
         series%cells(:, 3, count) = values(7:9)
      END DO
      CLOSE (cell_unit)
   END SUBROUTINE read_cell_series

   SUBROUTINE select_cell(series, frame_index, trajectory_step, cell)
      TYPE(cell_series_type), INTENT(IN)                 :: series
      INTEGER, INTENT(IN)                                :: frame_index, trajectory_step
      REAL(KIND=dp), INTENT(OUT)                         :: cell(3, 3)

      INTEGER                                            :: i

      IF (trajectory_step /= HUGE(1)) THEN
         DO i = 1, SIZE(series%steps)
            IF (series%steps(i) == trajectory_step) THEN
               cell = series%cells(:, :, i)
               RETURN
            END IF
         END DO
         CALL fail("No CELL_FILE record matches trajectory step")
      ELSE
         IF (frame_index > SIZE(series%steps)) CALL fail("CELL_FILE has fewer records than the trajectory")
         cell = series%cells(:, :, frame_index)
      END IF
   END SUBROUTINE select_cell

   SUBROUTINE cell_from_xyz_comment(comment, cell)
      CHARACTER(LEN=*), INTENT(IN)                       :: comment
      REAL(KIND=dp), INTENT(OUT)                         :: cell(3, 3)

      CHARACTER(LEN=LEN(comment))                        :: work
      INTEGER                                            :: first_quote, ios, lattice_position, &
                                                            second_quote
      REAL(KIND=dp)                                      :: values(9)

      work = uppercase(comment)
      lattice_position = INDEX(work, "LATTICE=")
      IF (lattice_position == 0) &
         CALL fail("No cell was configured and the XYZ comment has no Lattice field")
      first_quote = INDEX(comment(lattice_position:), '"')
      IF (first_quote == 0) CALL fail("Malformed Lattice field in extended XYZ comment")
      first_quote = lattice_position + first_quote - 1
      second_quote = INDEX(comment(first_quote + 1:), '"')
      IF (second_quote == 0) CALL fail("Malformed Lattice field in extended XYZ comment")
      second_quote = first_quote + second_quote
      READ (comment(first_quote + 1:second_quote - 1), *, IOSTAT=ios) values
      IF (ios /= 0) CALL fail("Malformed Lattice values in extended XYZ comment")
      cell(:, 1) = values(1:3)
      cell(:, 2) = values(4:6)
      cell(:, 3) = values(7:9)
   END SUBROUTINE cell_from_xyz_comment

   SUBROUTINE read_xyz_frame(trajectory_unit, labels, coordinates, comment, natoms, end_of_file)
      INTEGER, INTENT(IN)                                :: trajectory_unit
      CHARACTER(LEN=label_length), ALLOCATABLE, &
         INTENT(INOUT)                                   :: labels(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: coordinates(:, :)
      CHARACTER(LEN=1024), INTENT(OUT)                   :: comment
      INTEGER, INTENT(OUT)                               :: natoms
      LOGICAL, INTENT(OUT)                               :: end_of_file

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: i, ios

      end_of_file = .FALSE.
      DO
         READ (trajectory_unit, '(A)', IOSTAT=ios) line
         IF (ios < 0) THEN
            end_of_file = .TRUE.
            natoms = 0
            RETURN
         END IF
         IF (ios > 0) CALL fail("Could not read XYZ trajectory")
         IF (LEN_TRIM(line) > 0) EXIT
      END DO
      READ (line, *, IOSTAT=ios) natoms
      IF (ios /= 0 .OR. natoms < 1) CALL fail("Malformed atom count in XYZ trajectory")

      IF (.NOT. ALLOCATED(labels)) THEN
         ALLOCATE (labels(natoms), coordinates(3, natoms))
      ELSE IF (SIZE(labels) /= natoms) THEN
         CALL fail("The number of atoms changes in the XYZ trajectory")
      END IF

      READ (trajectory_unit, '(A)', IOSTAT=ios) comment
      IF (ios /= 0) CALL fail("Missing XYZ comment line")
      DO i = 1, natoms
         READ (trajectory_unit, '(A)', IOSTAT=ios) line
         IF (ios /= 0) CALL fail("Incomplete XYZ frame")
         READ (line, *, IOSTAT=ios) labels(i), coordinates(:, i)
         IF (ios /= 0) CALL fail("Malformed XYZ atom line: "//TRIM(line))
      END DO
   END SUBROUTINE read_xyz_frame

   FUNCTION trajectory_step(comment) RESULT(step)
      CHARACTER(LEN=*), INTENT(IN)                       :: comment
      INTEGER                                            :: step

      CHARACTER(LEN=LEN(comment))                        :: work
      INTEGER                                            :: comma, equals, ios, position

      step = HUGE(1)
      work = uppercase(comment)
      position = INDEX(work, "I =")
      IF (position == 0) position = INDEX(work, "STEP=")
      IF (position == 0) position = INDEX(work, "STEP =")
      IF (position == 0) RETURN
      equals = INDEX(work(position:), "=")
      IF (equals == 0) RETURN
      work = ADJUSTL(work(position + equals:))
      comma = INDEX(work, ",")
      IF (comma > 0) work(comma:) = ""
      READ (work, *, IOSTAT=ios) step
      IF (ios /= 0) step = HUGE(1)
   END FUNCTION trajectory_step

   PURE FUNCTION matches_selector(label, selector) RESULT(matches)
      CHARACTER(LEN=*), INTENT(IN)                       :: label, selector
      LOGICAL                                            :: matches

      CHARACTER(LEN=label_length)                        :: token
      CHARACTER(LEN=label_length+1)                      :: work
      INTEGER                                            :: comma

      matches = .FALSE.
      work = TRIM(uppercase(selector))//","
      IF (TRIM(uppercase(selector)) == "*") THEN
         matches = .TRUE.
         RETURN
      END IF
      DO WHILE (LEN_TRIM(work) > 0)
         comma = INDEX(work, ",")
         IF (comma == 0) RETURN
         token = ADJUSTL(work(:comma - 1))
         IF (TRIM(uppercase(label)) == TRIM(token)) THEN
            matches = .TRUE.
            RETURN
         END IF
         work = ADJUSTL(work(comma + 1:))
      END DO
   END FUNCTION matches_selector

   PURE FUNCTION safe_name(name) RESULT(RESULT)
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      CHARACTER(LEN=LEN(name))                           :: result

      INTEGER                                            :: code, i

      RESULT = name
      DO i = 1, LEN_TRIM(RESULT)
         code = IACHAR(RESULT(i:i))
         IF (.NOT. ((code >= IACHAR("A") .AND. code <= IACHAR("Z")) .OR. &
                    (code >= IACHAR("a") .AND. code <= IACHAR("z")) .OR. &
                    (code >= IACHAR("0") .AND. code <= IACHAR("9")) .OR. &
                    RESULT(i:i) == "-" .OR. RESULT(i:i) == "_")) RESULT(i:i) = "_"
      END DO
   END FUNCTION safe_name

   FUNCTION atomic_number(label) RESULT(number)
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      INTEGER                                            :: number

      CHARACTER(LEN=2)                                   :: candidate
      INTEGER                                            :: i
      CHARACTER(LEN=2), PARAMETER :: symbols(118) = [CHARACTER(LEN=2) :: "H", "He", "Li", "Be", "B"&
         , "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",&
         "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", &
         "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", &
         "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", &
         "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", &
         "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", &
         "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", &
         "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

      candidate = "  "
      IF (LEN_TRIM(label) > 0) candidate(1:1) = uppercase(label(1:1))
      IF (LEN_TRIM(label) > 1) THEN
         candidate(2:2) = label(2:2)
         IF (candidate(2:2) >= "A" .AND. candidate(2:2) <= "Z") &
            candidate(2:2) = ACHAR(IACHAR(candidate(2:2)) + 32)
      END IF
      number = 0
      DO i = 1, SIZE(symbols)
         IF (candidate == symbols(i)) THEN
            number = i
            RETURN
         END IF
      END DO
      candidate(2:2) = " "
      DO i = 1, SIZE(symbols)
         IF (candidate == symbols(i)) THEN
            number = i
            RETURN
         END IF
      END DO
   END FUNCTION atomic_number

   SUBROUTINE write_cube(config, target_index, reference_labels, reference_coordinates, nframes, nreferences)
      TYPE(config_type), INTENT(IN)                      :: config
      INTEGER, INTENT(IN)                                :: target_index
      CHARACTER(LEN=label_length), INTENT(IN)            :: reference_labels(3)
      REAL(KIND=dp), INTENT(IN)                          :: reference_coordinates(3, 3)
      INTEGER, INTENT(IN)                                :: nframes, nreferences

      CHARACTER(LEN=path_length)                         :: filename
      INTEGER                                            :: cube_unit, i, j, k, number
      REAL(KIND=dp)                                      :: origin(3), SPACING(3), value

      filename = TRIM(config%output_prefix)//"-"//TRIM(safe_name(config%targets(target_index)%name))//".cube"
      OPEN (NEWUNIT=cube_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write")
      WRITE (cube_unit, '(A)') "Spatial distribution function generated from a CP2K XYZ trajectory"
      WRITE (cube_unit, '(A,I0,A,I0)') "Target "//TRIM(config%targets(target_index)%name)// &
         "; sampled frames ", nframes, "; reference observations ", nreferences

      spacing = (config%upper - config%lower)/REAL(config%grid, KIND=dp)
      origin = (config%lower + 0.5_dp*spacing)*angstrom_to_bohr
      WRITE (cube_unit, '(I5,3F14.7)') 3, origin
      WRITE (cube_unit, '(I5,3F14.7)') config%grid(1), SPACING(1)*angstrom_to_bohr, 0.0_dp, 0.0_dp
      WRITE (cube_unit, '(I5,3F14.7)') config%grid(2), 0.0_dp, SPACING(2)*angstrom_to_bohr, 0.0_dp
      WRITE (cube_unit, '(I5,3F14.7)') config%grid(3), 0.0_dp, 0.0_dp, SPACING(3)*angstrom_to_bohr
      DO i = 1, 3
         number = atomic_number(reference_labels(i))
         WRITE (cube_unit, '(I5,4F14.7)') number, REAL(number, KIND=dp), &
            reference_coordinates(:, i)*angstrom_to_bohr
      END DO

      DO i = 1, config%grid(1)
         DO j = 1, config%grid(2)
            DO k = 1, config%grid(3)
               value = config%targets(target_index)%histogram(i, j, k)
               WRITE (cube_unit, '(1X,ES13.5E3)', ADVANCE="no") value
               IF (MOD(k, 6) == 0 .OR. k == config%grid(3)) WRITE (cube_unit, '()')
            END DO
         END DO
      END DO
      CLOSE (cube_unit)
      WRITE (*, '(A)') "Wrote "//TRIM(filename)
   END SUBROUTINE write_cube

   SUBROUTINE write_reference_xyz(config, labels, coordinates)
      TYPE(config_type), INTENT(IN)                      :: config
      CHARACTER(LEN=label_length), INTENT(IN)            :: labels(3)
      REAL(KIND=dp), INTENT(IN)                          :: coordinates(3, 3)

      CHARACTER(LEN=path_length)                         :: filename
      INTEGER                                            :: i, xyz_unit

      filename = TRIM(config%output_prefix)//"-reference.xyz"
      OPEN (NEWUNIT=xyz_unit, FILE=TRIM(filename), STATUS="replace", ACTION="write")
      WRITE (xyz_unit, '(I0)') 3
      WRITE (xyz_unit, '(A)') "Average reference geometry in the local SDF frame; coordinates in angstrom"
      DO i = 1, 3
         WRITE (xyz_unit, '(A,3(1X,F16.9))') TRIM(labels(i)), coordinates(:, i)
      END DO
      CLOSE (xyz_unit)
      WRITE (*, '(A)') "Wrote "//TRIM(filename)
   END SUBROUTINE write_reference_xyz

END MODULE spatial_distribution_io

PROGRAM spatial_distribution
   USE spatial_distribution_geometry, ONLY: inverse_matrix, local_basis, minimum_image
   USE spatial_distribution_input, ONLY: read_config
   USE spatial_distribution_io, ONLY: cell_from_xyz_comment, cell_series_type, matches_selector, &
                                      read_cell_series, read_xyz_frame, select_cell, trajectory_step, &
                                      write_cube, write_reference_xyz
   USE spatial_distribution_types, ONLY: config_type, dp, fail, label_length, path_length
   IMPLICIT NONE

   TYPE(config_type) :: config
   TYPE(cell_series_type) :: cell_series
   CHARACTER(LEN=path_length) :: config_filename
   CHARACTER(LEN=1024) :: comment, message
   CHARACTER(LEN=label_length), ALLOCATABLE :: labels(:)
   CHARACTER(LEN=label_length) :: reference_labels(3)
   REAL(KIND=dp), ALLOCATABLE :: coordinates(:, :)
   REAL(KIND=dp) :: arms(3, 2), basis(3, 3), cell(3, 3), cell_inverse(3, 3)
   REAL(KIND=dp) :: determinant_value, displacement(3), grid_spacing(3), local(3)
   REAL(KIND=dp) :: reference_sum(3, 3), voxel_volume, volume
   INTEGER, ALLOCATABLE :: target_counts(:)
   INTEGER :: atom_index, frame, ios, ix, iy, iz, natoms, reference_index, sampled_frames, target_index
   INTEGER :: trajectory_unit
   INTEGER :: total_references
   LOGICAL :: end_of_file, selected

   IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
      WRITE (*, '(A)') "Usage: spatial_distribution.x INPUT_FILE"
      STOP 1
   END IF
   CALL GET_COMMAND_ARGUMENT(1, config_filename)
   CALL read_config(TRIM(config_filename), config)
   IF (LEN_TRIM(config%cell_file) > 0) CALL read_cell_series(config%cell_file, cell_series)

   OPEN (NEWUNIT=trajectory_unit, FILE=TRIM(config%trajectory), STATUS="old", ACTION="read", IOSTAT=ios)
   IF (ios /= 0) CALL fail("Could not open trajectory "//TRIM(config%trajectory))
   ALLOCATE (target_counts(SIZE(config%targets)))
   grid_spacing = (config%upper - config%lower)/REAL(config%grid, KIND=dp)
   voxel_volume = PRODUCT(grid_spacing)
   frame = 0
   sampled_frames = 0
   total_references = 0
   reference_sum = 0.0_dp
   reference_labels = ""

   DO
      CALL read_xyz_frame(trajectory_unit, labels, coordinates, comment, natoms, end_of_file)
      IF (end_of_file) EXIT
      frame = frame + 1
      IF (frame > config%last_frame) EXIT
      selected = frame >= config%first_frame .AND. MOD(frame - config%first_frame, config%stride) == 0
      IF (.NOT. selected) CYCLE

      IF (config%static_cell_set) THEN
         cell = config%static_cell
      ELSE IF (ALLOCATED(cell_series%steps)) THEN
         CALL select_cell(cell_series, frame, trajectory_step(comment), cell)
      ELSE
         CALL cell_from_xyz_comment(comment, cell)
      END IF
      CALL inverse_matrix(cell, cell_inverse, determinant_value)
      volume = ABS(determinant_value)
      CALL validate_indices(config, natoms)

      target_counts = 0
      DO target_index = 1, SIZE(config%targets)
         DO atom_index = 1, natoms
            IF (matches_selector(labels(atom_index), config%targets(target_index)%selector)) &
               target_counts(target_index) = target_counts(target_index) + 1
         END DO
         config%targets(target_index)%normalization_weight = &
            config%targets(target_index)%normalization_weight + &
            REAL(SIZE(config%references)*target_counts(target_index), KIND=dp)/volume
      END DO

      DO reference_index = 1, SIZE(config%references)
         CALL local_basis(config%references(reference_index), coordinates, cell, cell_inverse, config, basis, arms)
         total_references = total_references + 1
         reference_sum(:, 2:3) = reference_sum(:, 2:3) + arms
         IF (total_references == 1) THEN
            reference_labels(1) = labels(config%references(reference_index)%origin)
            reference_labels(2) = labels(config%references(reference_index)%arm1)
            reference_labels(3) = labels(config%references(reference_index)%arm2)
         END IF

         DO target_index = 1, SIZE(config%targets)
            DO atom_index = 1, natoms
               IF (.NOT. matches_selector(labels(atom_index), config%targets(target_index)%selector)) CYCLE
               IF (config%exclude_reference .AND. is_reference_atom(atom_index, config, reference_index)) CYCLE
               displacement = minimum_image(coordinates(:, atom_index) - &
                                            coordinates(:, config%references(reference_index)%origin), &
                                            cell, cell_inverse, config%periodic_axes)
               local = MATMUL(TRANSPOSE(basis), displacement)
               IF (config%cutoff > 0.0_dp) THEN
                  IF (DOT_PRODUCT(local, local) > config%cutoff**2) CYCLE
               END IF
               IF (ANY(local < config%lower) .OR. ANY(local >= config%upper)) CYCLE
               ix = FLOOR((local(1) - config%lower(1))/grid_spacing(1)) + 1
               iy = FLOOR((local(2) - config%lower(2))/grid_spacing(2)) + 1
               iz = FLOOR((local(3) - config%lower(3))/grid_spacing(3)) + 1
               config%targets(target_index)%histogram(ix, iy, iz) = &
                  config%targets(target_index)%histogram(ix, iy, iz) + 1.0_dp
               config%targets(target_index)%observations = config%targets(target_index)%observations + 1
            END DO
         END DO
      END DO
      sampled_frames = sampled_frames + 1
   END DO
   CLOSE (trajectory_unit)

   IF (sampled_frames == 0) CALL fail("No trajectory frames matched the requested range")
   IF (total_references == 0) CALL fail("No reference observations were accumulated")
   reference_sum = reference_sum/REAL(total_references, KIND=dp)

   WRITE (*, '(A,I0)') "Sampled trajectory frames: ", sampled_frames
   WRITE (*, '(A,I0)') "Reference observations: ", total_references
   DO target_index = 1, SIZE(config%targets)
      IF (config%targets(target_index)%normalization_weight <= 0.0_dp) THEN
         CALL fail("TARGET "//TRIM(config%targets(target_index)%name)//" selects no atoms")
      END IF
      config%targets(target_index)%histogram = config%targets(target_index)%histogram/ &
                                               (voxel_volume*config%targets(target_index)%normalization_weight)
      WRITE (message, '(A,A,A,I0)') "Target ", TRIM(config%targets(target_index)%name), &
         " observations: ", config%targets(target_index)%observations
      WRITE (*, '(A)') TRIM(message)
      CALL write_cube(config, target_index, reference_labels, reference_sum, sampled_frames, total_references)
   END DO
   CALL write_reference_xyz(config, reference_labels, reference_sum)
   CALL release_memory(config, cell_series, labels, coordinates, target_counts)

CONTAINS

   SUBROUTINE release_memory(config, cell_series, labels, coordinates, target_counts)
      TYPE(config_type), INTENT(INOUT)                   :: config
      TYPE(cell_series_type), INTENT(INOUT)              :: cell_series
      CHARACTER(LEN=label_length), ALLOCATABLE, &
         INTENT(INOUT)                                   :: labels(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: coordinates(:, :)
      INTEGER, ALLOCATABLE, INTENT(INOUT)                :: target_counts(:)

      INTEGER                                            :: target_index

      IF (ALLOCATED(labels)) DEALLOCATE (labels)
      IF (ALLOCATED(coordinates)) DEALLOCATE (coordinates)
      IF (ALLOCATED(target_counts)) DEALLOCATE (target_counts)
      IF (ALLOCATED(cell_series%steps)) DEALLOCATE (cell_series%steps)
      IF (ALLOCATED(cell_series%cells)) DEALLOCATE (cell_series%cells)
      IF (ALLOCATED(config%targets)) THEN
         DO target_index = 1, SIZE(config%targets)
            IF (ALLOCATED(config%targets(target_index)%histogram)) &
               DEALLOCATE (config%targets(target_index)%histogram)
         END DO
         DEALLOCATE (config%targets)
      END IF
      IF (ALLOCATED(config%references)) DEALLOCATE (config%references)
   END SUBROUTINE release_memory

   SUBROUTINE validate_indices(config, natoms)
      TYPE(config_type), INTENT(IN)                      :: config
      INTEGER, INTENT(IN)                                :: natoms

      INTEGER                                            :: reference_index

      DO reference_index = 1, SIZE(config%references)
         IF (MAX(config%references(reference_index)%origin, config%references(reference_index)%arm1, &
                 config%references(reference_index)%arm2) > natoms) &
            CALL fail("A REFERENCE index exceeds the number of atoms in the trajectory")
      END DO
   END SUBROUTINE validate_indices

   PURE FUNCTION is_reference_atom(atom_index, config, reference_index) RESULT(is_reference)
      INTEGER, INTENT(IN)                                :: atom_index
      TYPE(config_type), INTENT(IN)                      :: config
      INTEGER, INTENT(IN)                                :: reference_index
      LOGICAL                                            :: is_reference

      is_reference = atom_index == config%references(reference_index)%origin .OR. &
                     atom_index == config%references(reference_index)%arm1 .OR. &
                     atom_index == config%references(reference_index)%arm2
   END FUNCTION is_reference_atom

END PROGRAM spatial_distribution
