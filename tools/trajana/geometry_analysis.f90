!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_geometry_analysis
   USE trajana_cell_source,             ONLY: cell_source_type
   USE trajana_command_line,            ONLY: fail,&
                                              get_integer_option,&
                                              get_option
   USE trajana_frame_controls,          ONLY: frame_selected
   USE trajana_geometry,                ONLY: minimum_image
   USE trajana_kinds,                   ONLY: dp
   USE trajana_text_utils,              ONLY: lower_case,&
                                              next_token
   USE trajana_trajectory_io,           ONLY: close_output,&
                                              open_output,&
                                              xyz_reader_type
   USE trajana_trajectory_types,        ONLY: frame_type

   IMPLICIT NONE
   PRIVATE

   TYPE :: action_type
      INTEGER :: kind = 0
      INTEGER :: atom(4) = 0
   END TYPE action_type

   PUBLIC :: action_type, evaluate_action, read_actions, run_geometry, validate_actions

CONTAINS

   SUBROUTINE run_geometry()
      CHARACTER(LEN=512)                                 :: message
      CHARACTER(LEN=:), ALLOCATABLE                      :: action_path, cell_path, cell_text, &
                                                            input_path, output_path, periodic_text
      INTEGER                                            :: expected_atoms, first, ierr, item, last, &
                                                            output_unit, stride
      LOGICAL                                            :: eof, found
      REAL(dp), ALLOCATABLE                              :: result(:)
      TYPE(action_type), ALLOCATABLE                     :: actions(:)
      TYPE(cell_source_type)                             :: cells
      TYPE(frame_type)                                   :: frame
      TYPE(xyz_reader_type)                              :: trajectory

      CALL get_option("--input", input_path, found, "-")
      CALL get_option("--output", output_path, found, "-")
      CALL get_option("--actions", action_path, found, "trajana.in")
      CALL get_option("--cell", cell_text, found, "")
      CALL get_option("--cell-file", cell_path, found, "")
      CALL get_option("--periodic", periodic_text, found, "XYZ")
      CALL get_integer_option("--first", first, 1)
      CALL get_integer_option("--last", last, HUGE(1))
      CALL get_integer_option("--stride", stride, 1)
      IF (first < 1 .OR. last < first .OR. stride < 1) CALL fail("Invalid frame range")

      CALL read_actions(action_path, actions, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL cells%configure(cell_text, cell_path, periodic_text, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL trajectory%open_file(input_path, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      CALL open_output(output_path, output_unit, ierr, message)
      IF (ierr /= 0) CALL fail(message)
      WRITE (output_unit, "(A)", ADVANCE="no") "# frame"
      DO item = 1, SIZE(actions)
         WRITE (output_unit, "(A,I0)", ADVANCE="no") "   value_", item
      END DO
      WRITE (output_unit, "(A)") ""
      ALLOCATE (RESULT(SIZE(actions)))
      expected_atoms = -1

      DO
         CALL trajectory%read_frame(frame, eof, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (eof) EXIT
         CALL cells%attach(frame, ierr, message)
         IF (ierr /= 0) CALL fail(message)
         IF (frame%number > last) EXIT
         IF (.NOT. frame_selected(frame%number, first, last, stride)) CYCLE
         IF (expected_atoms < 0) THEN
            expected_atoms = frame%natoms
            CALL validate_actions(actions, expected_atoms, ierr, message)
            IF (ierr /= 0) CALL fail(message)
         END IF
         IF (frame%natoms /= expected_atoms) CALL fail("The atom count changes along the trajectory")
         DO item = 1, SIZE(actions)
            CALL evaluate_action(actions(item), frame, RESULT(item), ierr)
            IF (ierr /= 0) CALL fail("Degenerate geometry or singular cell in frame")
         END DO
         WRITE (output_unit, "(I12,*(1X,ES24.16))") frame%number, RESULT
      END DO

      CALL close_output(output_unit)
      CALL trajectory%close_file()
      CALL cells%close()
   END SUBROUTINE run_geometry

   SUBROUTINE read_actions(path, actions, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      TYPE(action_type), ALLOCATABLE, INTENT(OUT)        :: actions(:)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=128)                                 :: token
      CHARACTER(LEN=4096)                                :: line
      INTEGER                                            :: count, ios, item, position, required, &
                                                            unit
      LOGICAL                                            :: found
      TYPE(action_type), ALLOCATABLE                     :: grown(:)

      ierr = 0
      message = ""
      count = 0
      ALLOCATE (actions(0))
      OPEN (NEWUNIT=unit, FILE=TRIM(path), STATUS="old", ACTION="read", IOSTAT=ios)
      IF (ios /= 0) THEN
         ierr = ios
         message = "Cannot open action file "//TRIM(path)
         RETURN
      END IF
      DO
         READ (unit, "(A)", IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) EXIT
         IF (ios /= 0) THEN
            ierr = ios
            message = "Error while reading action file"
            EXIT
         END IF
         line = ADJUSTL(line)
         IF (LEN_TRIM(line) == 0 .OR. INDEX("!#", line(1:1)) > 0) CYCLE
         position = 1
         CALL next_token(line, position, token, found)
         IF (.NOT. found) CYCLE
         SELECT CASE (lower_case(token(1:1)))
         CASE ("d")
            required = 2
         CASE ("a")
            required = 3
         CASE ("t")
            required = 4
         CASE ("c")
            CYCLE
         CASE DEFAULT
            ierr = 1
            message = "Unknown geometry action: "//TRIM(token)
            EXIT
         END SELECT
         ALLOCATE (grown(count + 1))
         IF (count > 0) grown(:count) = actions
         grown(count + 1)%kind = required - 1
         DO item = 1, required
            CALL next_token(line, position, token, found)
            IF (.NOT. found) THEN
               ierr = 1
               message = "Too few atom indices in action: "//TRIM(line)
               EXIT
            END IF
            READ (token, *, IOSTAT=ios) grown(count + 1)%atom(item)
            IF (ios /= 0) THEN
               ierr = 1
               message = "Invalid atom index in action: "//TRIM(line)
               EXIT
            END IF
         END DO
         IF (ierr /= 0) EXIT
         CALL MOVE_ALLOC(grown, actions)
         count = count + 1
      END DO
      CLOSE (unit)
      IF (ierr == 0 .AND. count == 0) THEN
         ierr = 1
         message = "The action file contains no actions"
      END IF
   END SUBROUTINE read_actions

   SUBROUTINE validate_actions(actions, natoms, ierr, message)
      TYPE(action_type), INTENT(IN)                      :: actions(:)
      INTEGER, INTENT(IN)                                :: natoms
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: action, item, required

      ierr = 0
      message = ""
      DO action = 1, SIZE(actions)
         required = actions(action)%kind + 1
         DO item = 1, required
            IF (actions(action)%atom(item) < 1 .OR. actions(action)%atom(item) > natoms) THEN
               ierr = 1
               WRITE (message, "(A,I0)") "Invalid atom index in action ", action
               RETURN
            END IF
         END DO
      END DO
   END SUBROUTINE validate_actions

   SUBROUTINE evaluate_action(action, frame, RESULT, ierr)
      TYPE(action_type), INTENT(IN)                      :: action
      TYPE(frame_type), INTENT(IN)                       :: frame
      REAL(dp), INTENT(OUT)                              :: result
      INTEGER, INTENT(OUT)                               :: ierr

      LOGICAL                                            :: ok
      REAL(dp)                                           :: a(3), b(3), c(3), m1(3), n1(3), n2(3), &
                                                            norm_a, norm_b, norm_c, x, y

      ierr = 0
      SELECT CASE (action%kind)
      CASE (1)
         CALL displacement(frame, action%atom(1), action%atom(2), a, ok)
         IF (.NOT. ok) THEN
            ierr = 1
            RETURN
         END IF
         RESULT = NORM2(a)
      CASE (2)
         CALL displacement(frame, action%atom(2), action%atom(1), a, ok)
         IF (.NOT. ok) THEN
            ierr = 1
            RETURN
         END IF
         CALL displacement(frame, action%atom(2), action%atom(3), b, ok)
         norm_a = NORM2(a)
         norm_b = NORM2(b)
         IF (.NOT. ok .OR. norm_a <= TINY(1.0_dp) .OR. norm_b <= TINY(1.0_dp)) THEN
            ierr = 1
            RETURN
         END IF
         x = MAX(-1.0_dp, MIN(1.0_dp, DOT_PRODUCT(a, b)/(norm_a*norm_b)))
         RESULT = ACOS(x)*180.0_dp/ACOS(-1.0_dp)
      CASE (3)
         CALL displacement(frame, action%atom(1), action%atom(2), a, ok)
         IF (.NOT. ok) THEN
            ierr = 1
            RETURN
         END IF
         CALL displacement(frame, action%atom(2), action%atom(3), b, ok)
         IF (.NOT. ok) THEN
            ierr = 1
            RETURN
         END IF
         CALL displacement(frame, action%atom(3), action%atom(4), c, ok)
         IF (.NOT. ok) THEN
            ierr = 1
            RETURN
         END IF
         n1 = cross_product(a, b)
         n2 = cross_product(b, c)
         norm_a = NORM2(n1)
         norm_b = NORM2(n2)
         norm_c = NORM2(b)
         IF (MIN(norm_a, norm_b, norm_c) <= TINY(1.0_dp)) THEN
            ierr = 1
            RETURN
         END IF
         n1 = n1/norm_a
         n2 = n2/norm_b
         m1 = cross_product(n1, b/norm_c)
         x = DOT_PRODUCT(n1, n2)
         y = DOT_PRODUCT(m1, n2)
         RESULT = ATAN2(y, x)*180.0_dp/ACOS(-1.0_dp)
      END SELECT
   END SUBROUTINE evaluate_action

   SUBROUTINE displacement(frame, from_atom, to_atom, vector, ok)
      TYPE(frame_type), INTENT(IN)                       :: frame
      INTEGER, INTENT(IN)                                :: from_atom, to_atom
      REAL(dp), INTENT(OUT)                              :: vector(3)
      LOGICAL, INTENT(OUT)                               :: ok

      REAL(dp)                                           :: raw(3)

      raw = frame%value(:, to_atom) - frame%value(:, from_atom)
      IF (frame%cell%valid) THEN
         CALL minimum_image(frame%cell, raw, vector, ok)
      ELSE
         vector = raw
         ok = .TRUE.
      END IF
   END SUBROUTINE displacement

   PURE FUNCTION cross_product(a, b) RESULT(c)
      REAL(dp), INTENT(IN)                               :: a(3), b(3)
      REAL(dp)                                           :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END FUNCTION cross_product

END MODULE trajana_geometry_analysis
