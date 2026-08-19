!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM trajconvert
   USE trajana_cell_source,             ONLY: cell_source_type
   USE trajana_command_line,            ONLY: command_name,&
                                              fail,&
                                              get_integer_option,&
                                              get_option,&
                                              has_flag
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry,                ONLY: cart_to_fractional,&
                                              fractional_to_cart,&
                                              minimum_image,&
                                              wrap_cartesian,&
                                              wrap_fractional
   USE trajana_groups,                  ONLY: group_type,&
                                              mass_table_type,&
                                              read_groups,&
                                              read_masses,&
                                              validate_groups
   USE trajana_kinds,                   ONLY: dp
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              write_xyz_frame,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE

   CHARACTER(LEN=:), ALLOCATABLE :: operation, input_path, output_path, cell_text, cell_path
   CHARACTER(LEN=:), ALLOCATABLE :: periodic_text, group_path, mass_path, center_mode
   CHARACTER(LEN=512) :: message
   TYPE(xyz_reader_type) :: trajectory
   TYPE(cell_source_type) :: cells
   TYPE(frame_type) :: frame
   TYPE(group_type), ALLOCATABLE :: groups(:)
   TYPE(mass_table_type) :: masses
   REAL(dp), ALLOCATABLE :: output_values(:, :), previous_fractional(:, :), images(:, :)
   CHARACTER(LEN=32), ALLOCATABLE :: output_labels(:)
   INTEGER :: ierr, output_unit, first, last, stride, expected_atoms
   LOGICAL :: eof, found, vectors, wrap_output, initialized

   operation = lower_case(command_name())
   IF (operation == "" .OR. operation == "help" .OR. operation == "--help" .OR. operation == "-h" .OR. &
       has_flag("--help")) THEN
      CALL print_help()
      STOP
   END IF
   IF (operation /= "wrap" .AND. operation /= "unwrap" .AND. operation /= "fractional" .AND. operation /= "center") THEN
      CALL fail("Unknown trajconvert operation: "//operation)
   END IF

   CALL get_option("--input", input_path, found, "-")
   CALL get_option("--output", output_path, found, "-")
   CALL get_option("--cell", cell_text, found, "")
   CALL get_option("--cell-file", cell_path, found, "")
   CALL get_option("--periodic", periodic_text, found, "XYZ")
   CALL get_integer_option("--first", first, 1)
   CALL get_integer_option("--last", last, HUGE(1))
   CALL get_integer_option("--stride", stride, 1)
   IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")

   CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
   IF (ierr /= 0) CALL fail(message)
   CALL trajectory%open_file(input_path, ierr, message)
   IF (ierr /= 0) CALL fail(message)
   CALL open_output(output_path, output_unit, ierr, message)
   IF (ierr /= 0) CALL fail(message)

   vectors = has_flag("--vectors")
   wrap_output = has_flag("--wrap-output")
   initialized = .FALSE.
   expected_atoms = -1

   IF (operation == "center") THEN
      CALL get_option("--groups", group_path, found)
      IF (.NOT. found) CALL fail("center requires --groups FILE")
      CALL get_option("--mode", center_mode, found, "geometry")
      center_mode = lower_case(center_mode)
      IF (center_mode /= "geometry" .AND. center_mode /= "mass") CALL fail("--mode expects geometry or mass")
      CALL read_groups(group_path, groups, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (center_mode == "mass") THEN
         CALL get_option("--mass-file", mass_path, found)
         IF (.NOT. found) CALL fail("--mode mass requires --mass-file FILE")
         CALL read_masses(mass_path, masses, ierr, message)
         IF (ierr /= 0) CALL fail(message)
      END IF
   END IF

   DO
      CALL trajectory%read_frame(frame, eof, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (eof) EXIT
      CALL cells%attach(frame, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      IF (expected_atoms < 0) expected_atoms = frame%natoms
      IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")
      IF (frame%number > last) EXIT

      SELECT CASE (operation)
      CASE ("unwrap")
         IF (.NOT. frame%cell%valid) CALL fail("unwrap requires --cell, --cell-file, or an extended-XYZ Lattice")
         CALL unwrap_frame(frame, previous_fractional, images, output_values, initialized, ierr)
         IF (ierr /= 0) CALL fail("Singular cell in trajectory")
         IF (frame_selected(frame%number, first, last, stride)) &
            CALL write_xyz_frame(output_unit, frame%label, output_values, "Unwrapped; "//frame%comment)
      CASE DEFAULT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         SELECT CASE (operation)
         CASE ("wrap")
            IF (.NOT. frame%cell%valid) CALL fail("wrap requires --cell, --cell-file, or an extended-XYZ Lattice")
            CALL wrap_frame(frame, output_values, ierr)
            IF (ierr /= 0) CALL fail("Singular cell in trajectory")
            CALL write_xyz_frame(output_unit, frame%label, output_values, "Wrapped; "//frame%comment)
         CASE ("fractional")
            IF (.NOT. frame%cell%valid) CALL fail("fractional requires --cell, --cell-file, or an extended-XYZ Lattice")
            CALL fractional_frame(frame, output_values, ierr)
            IF (ierr /= 0) CALL fail("Singular cell in trajectory")
            CALL write_xyz_frame(output_unit, frame%label, output_values, "Fractional coordinates; "//frame%comment)
         CASE ("center")
            IF (.NOT. initialized) THEN
               CALL validate_groups(groups, frame%natoms, 1, ierr, message)
               IF (ierr /= 0) CALL fail(message)
               ALLOCATE (output_labels(SIZE(groups)))
               output_labels = groups%label
               initialized = .TRUE.
            END IF
            IF (wrap_output .AND. .NOT. frame%cell%valid) CALL fail("--wrap-output requires cell information")
            CALL center_frame(frame, groups, masses, center_mode, vectors, wrap_output, output_values, ierr, message)
            IF (ierr /= 0) CALL fail(message)
            CALL write_xyz_frame(output_unit, output_labels, output_values, "Group centers; "//frame%comment)
         END SELECT
      END SELECT
   END DO

   CALL trajectory%close_file()
   CALL cells%close()
   CALL close_output(output_unit)

CONTAINS

   SUBROUTINE wrap_frame(current, values, status)
      TYPE(frame_type), INTENT(IN)                       :: current
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: values(:, :)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: atom
      LOGICAL                                            :: ok

      status = 0
      ALLOCATE (values(3, current%natoms))
      DO atom = 1, current%natoms
         CALL wrap_cartesian(current%cell, current%value(:, atom), values(:, atom), ok)
         IF (.NOT. ok) status = 1
      END DO
   END SUBROUTINE wrap_frame

   SUBROUTINE fractional_frame(current, values, status)
      TYPE(frame_type), INTENT(IN)                       :: current
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: values(:, :)
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: atom
      LOGICAL                                            :: ok

      status = 0
      ALLOCATE (values(3, current%natoms))
      DO atom = 1, current%natoms
         CALL cart_to_fractional(current%cell, current%value(:, atom), values(:, atom), ok)
         IF (.NOT. ok) THEN
            status = 1
            RETURN
         END IF
         CALL wrap_fractional(values(:, atom), current%cell%periodic)
      END DO
   END SUBROUTINE fractional_frame

   SUBROUTINE unwrap_frame(current, previous, image, values, ready, status)
      TYPE(frame_type), INTENT(IN)                       :: current
      REAL(dp), ALLOCATABLE, INTENT(INOUT)               :: previous(:, :), image(:, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: values(:, :)
      LOGICAL, INTENT(INOUT)                             :: ready
      INTEGER, INTENT(OUT)                               :: status

      INTEGER                                            :: atom, axis
      LOGICAL                                            :: ok
      REAL(dp)                                           :: continuous(3), delta, fractional(3)

      status = 0
      ALLOCATE (values(3, current%natoms))
      IF (.NOT. ready) THEN
         ALLOCATE (previous(3, current%natoms), image(3, current%natoms))
         image = 0.0_dp
      END IF
      DO atom = 1, current%natoms
         CALL cart_to_fractional(current%cell, current%value(:, atom), fractional, ok)
         IF (.NOT. ok) THEN
            status = 1
            RETURN
         END IF
         IF (ready) THEN
            DO axis = 1, 3
               IF (.NOT. current%cell%periodic(axis)) CYCLE
               delta = fractional(axis) - previous(axis, atom)
               image(axis, atom) = image(axis, atom) - ANINT(delta)
            END DO
         END IF
         continuous = fractional + image(:, atom)
         CALL fractional_to_cart(current%cell, continuous, values(:, atom))
         previous(:, atom) = fractional
      END DO
      ready = .TRUE.
   END SUBROUTINE unwrap_frame

   SUBROUTINE center_frame(current, definitions, mass_table, mode, vector_data, wrap_result, values, status, text)
      TYPE(frame_type), INTENT(IN)                       :: current
      TYPE(group_type), INTENT(IN)                       :: definitions(:)
      TYPE(mass_table_type), INTENT(IN)                  :: mass_table
      CHARACTER(LEN=*), INTENT(IN)                       :: mode
      LOGICAL, INTENT(IN)                                :: vector_data, wrap_result
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: values(:, :)
      INTEGER, INTENT(OUT)                               :: status
      CHARACTER(LEN=*), INTENT(OUT)                      :: text

      INTEGER                                            :: atom, group, item
      LOGICAL                                            :: mass_found, ok
      REAL(dp)                                           :: displacement(3), position(3), &
                                                            reference(3), total_weight, weight

      status = 0
      text = ""
      ALLOCATE (values(3, SIZE(definitions)))
      DO group = 1, SIZE(definitions)
         reference = current%value(:, definitions(group)%atom(1))
         values(:, group) = 0.0_dp
         total_weight = 0.0_dp
         DO item = 1, SIZE(definitions(group)%atom)
            atom = definitions(group)%atom(item)
            position = current%value(:, atom)
            IF (.NOT. vector_data .AND. current%cell%valid .AND. item > 1) THEN
               CALL minimum_image(current%cell, position - reference, displacement, ok)
               IF (.NOT. ok) THEN
                  status = 1
                  text = "Singular cell in trajectory"
                  RETURN
               END IF
               position = reference + displacement
            END IF
            IF (mode == "mass") THEN
               CALL mass_table%lookup(current%label(atom), weight, mass_found)
               IF (.NOT. mass_found) THEN
                  status = 1
                  text = "No mass provided for atom label "//TRIM(current%label(atom))
                  RETURN
               END IF
            ELSE
               weight = 1.0_dp
            END IF
            values(:, group) = values(:, group) + weight*position
            total_weight = total_weight + weight
         END DO
         values(:, group) = values(:, group)/total_weight
         IF (wrap_result) THEN
            CALL wrap_cartesian(current%cell, values(:, group), position, ok)
            IF (.NOT. ok) THEN
               status = 1
               text = "Singular cell in trajectory"
               RETURN
            END IF
            values(:, group) = position
         END IF
      END DO
   END SUBROUTINE center_frame

   SUBROUTINE print_help()
      WRITE (*, "(A)") "Usage: trajconvert.x OPERATION [options]"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") "Operations:"
      WRITE (*, "(A)") "  wrap        Wrap atoms into the periodic cell"
      WRITE (*, "(A)") "  unwrap      Reconstruct atomwise continuous trajectories"
      WRITE (*, "(A)") "  fractional  Convert Cartesian to fractional coordinates"
      WRITE (*, "(A)") "  center      Write centers for groups from --groups FILE"
      WRITE (*, "(A)") ""
      WRITE (*, "(A)") 'Common options: --input FILE --output FILE --cell "9 values"'
      WRITE (*, "(A)") "                --cell-file FILE --periodic XYZ --first N --last N --stride N"
      WRITE (*, "(A)") "Center options: --mode geometry|mass --mass-file FILE --vectors --wrap-output"
   END SUBROUTINE print_help

END PROGRAM trajconvert
