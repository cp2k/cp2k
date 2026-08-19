!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_trajectory_io
   USE iso_fortran_env,                 ONLY: input_unit,&
                                              output_unit
   USE trajana_kinds,                   ONLY: dp
   USE trajana_text_utils,              ONLY: lower_case
   USE trajana_trajectory_types,        ONLY: cell_type,&
                                              frame_type

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: xyz_reader_type
      INTEGER :: unit = -1
      INTEGER :: frame_number = 0
      CHARACTER(LEN=:), ALLOCATABLE :: path
   CONTAINS
      PROCEDURE :: open_file => xyz_open
      PROCEDURE :: read_frame => xyz_read_frame
      PROCEDURE :: close_file => xyz_close
   END TYPE xyz_reader_type

   TYPE, PUBLIC :: cell_reader_type
      INTEGER :: unit = -1
      INTEGER :: frame_number = 0
      CHARACTER(LEN=:), ALLOCATABLE :: path
   CONTAINS
      PROCEDURE :: open_file => cell_open
      PROCEDURE :: read_cell => cell_read
      PROCEDURE :: close_file => cell_close
   END TYPE cell_reader_type

   PUBLIC :: open_output, close_output, write_xyz_frame, parse_constant_cell

CONTAINS

   SUBROUTINE xyz_open(reader, path, ierr, message)
      CLASS(xyz_reader_type), INTENT(INOUT) :: reader
      CHARACTER(LEN=*), INTENT(IN) :: path
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(OUT) :: message

      message = ""
      reader%path = TRIM(path)
      IF (TRIM(path) == "-") THEN
         reader%unit = input_unit
         ierr = 0
      ELSE
         OPEN (NEWUNIT=reader%unit, FILE=TRIM(path), STATUS="old", ACTION="read", IOSTAT=ierr)
      END IF
      IF (ierr /= 0) WRITE (message, "(A,1X,A)") "Cannot open trajectory", TRIM(path)
   END SUBROUTINE xyz_open

   SUBROUTINE xyz_read_frame(reader, frame, eof, ierr, message)
      CLASS(xyz_reader_type), INTENT(INOUT) :: reader
      TYPE(frame_type), INTENT(OUT) :: frame
      LOGICAL, INTENT(OUT) :: eof
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(OUT) :: message
      CHARACTER(LEN=4096) :: line
      INTEGER :: atom, ios, natoms

      eof = .FALSE.
      ierr = 0
      message = ""

      DO
         READ (reader%unit, "(A)", IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) THEN
            eof = .TRUE.
            RETURN
         ELSE IF (ios /= 0) THEN
            ierr = ios
            message = "Error while reading the trajectory"
            RETURN
         END IF
         IF (LEN_TRIM(line) > 0) EXIT
      END DO

      READ (line, *, IOSTAT=ios) natoms
      IF (ios /= 0 .OR. natoms <= 0) THEN
         ierr = 1
         WRITE (message, "(A,I0)") "Invalid XYZ atom count before frame ", reader%frame_number + 1
         RETURN
      END IF

      READ (reader%unit, "(A)", IOSTAT=ios) line
      IF (ios /= 0) THEN
         ierr = 1
         message = "Missing XYZ comment line"
         RETURN
      END IF

      frame%natoms = natoms
      reader%frame_number = reader%frame_number + 1
      frame%number = reader%frame_number
      frame%comment = TRIM(line)
      ALLOCATE (frame%label(natoms), frame%value(3, natoms))

      DO atom = 1, natoms
         READ (reader%unit, "(A)", IOSTAT=ios) line
         IF (ios /= 0) THEN
            ierr = 1
            WRITE (message, "(A,I0,A,I0)") "Missing atom ", atom, " in frame ", frame%number
            RETURN
         END IF
         READ (line, *, IOSTAT=ios) frame%label(atom), frame%value(:, atom)
         IF (ios /= 0) THEN
            ierr = 1
            WRITE (message, "(A,I0,A,I0)") "Invalid atom ", atom, " in frame ", frame%number
            RETURN
         END IF
      END DO

      CALL parse_embedded_cell(frame%comment, frame%cell)
   END SUBROUTINE xyz_read_frame

   SUBROUTINE xyz_close(reader)
      CLASS(xyz_reader_type), INTENT(INOUT) :: reader

      IF (reader%unit /= -1 .AND. reader%unit /= input_unit) CLOSE (reader%unit)
      reader%unit = -1
   END SUBROUTINE xyz_close

   SUBROUTINE cell_open(reader, path, ierr, message)
      CLASS(cell_reader_type), INTENT(INOUT) :: reader
      CHARACTER(LEN=*), INTENT(IN) :: path
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(OUT) :: message

      message = ""
      reader%path = TRIM(path)
      OPEN (NEWUNIT=reader%unit, FILE=TRIM(path), STATUS="old", ACTION="read", IOSTAT=ierr)
      IF (ierr /= 0) WRITE (message, "(A,1X,A)") "Cannot open cell file", TRIM(path)
   END SUBROUTINE cell_open

   SUBROUTINE cell_read(reader, cell, eof, ierr, message)
      CLASS(cell_reader_type), INTENT(INOUT) :: reader
      TYPE(cell_type), INTENT(OUT) :: cell
      LOGICAL, INTENT(OUT) :: eof
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(OUT) :: message
      CHARACTER(LEN=4096) :: line
      INTEGER :: ios, step
      REAL(dp) :: time, volume

      eof = .FALSE.
      ierr = 0
      message = ""
      cell%valid = .FALSE.

      DO
         READ (reader%unit, "(A)", IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) THEN
            eof = .TRUE.
            RETURN
         ELSE IF (ios /= 0) THEN
            ierr = ios
            message = "Error while reading the cell file"
            RETURN
         END IF
         line = ADJUSTL(line)
         IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
         EXIT
      END DO

      READ (line, *, IOSTAT=ios) step, time, cell%h(:, 1), cell%h(:, 2), cell%h(:, 3), volume
      IF (ios /= 0) THEN
         ierr = 1
         WRITE (message, "(A,I0)") "Invalid CP2K cell record for frame ", reader%frame_number + 1
         RETURN
      END IF
      reader%frame_number = reader%frame_number + 1
      cell%valid = .TRUE.
   END SUBROUTINE cell_read

   SUBROUTINE cell_close(reader)
      CLASS(cell_reader_type), INTENT(INOUT) :: reader

      IF (reader%unit /= -1) CLOSE (reader%unit)
      reader%unit = -1
   END SUBROUTINE cell_close

   SUBROUTINE parse_embedded_cell(comment, cell)
      CHARACTER(LEN=*), INTENT(IN)                       :: comment
      TYPE(cell_type), INTENT(INOUT)                     :: cell

      CHARACTER(LEN=:), ALLOCATABLE                      :: lowered, values
      INTEGER                                            :: first, ios, last
      REAL(dp)                                           :: lattice(9)

      cell%valid = .FALSE.
      lowered = lower_case(comment)
      first = INDEX(lowered, 'lattice="')
      IF (first == 0) RETURN
      first = first + LEN('lattice="')
      last = INDEX(comment(first:), '"')
      IF (last == 0) RETURN
      values = comment(first:first + last - 2)
      READ (values, *, IOSTAT=ios) lattice
      IF (ios /= 0) RETURN
      cell%h(:, 1) = lattice(1:3)
      cell%h(:, 2) = lattice(4:6)
      cell%h(:, 3) = lattice(7:9)
      cell%valid = .TRUE.
   END SUBROUTINE parse_embedded_cell

   SUBROUTINE parse_constant_cell(text, cell, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      TYPE(cell_type), INTENT(OUT)                       :: cell
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      REAL(dp)                                           :: values(9)

      message = ""
      READ (text, *, IOSTAT=ierr) values
      IF (ierr /= 0) THEN
         message = "--cell expects nine numbers: ax ay az bx by bz cx cy cz"
         RETURN
      END IF
      cell%h(:, 1) = values(1:3)
      cell%h(:, 2) = values(4:6)
      cell%h(:, 3) = values(7:9)
      cell%valid = .TRUE.
   END SUBROUTINE parse_constant_cell

   SUBROUTINE open_output(path, unit, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      INTEGER, INTENT(OUT)                               :: unit, ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      message = ""
      IF (TRIM(path) == "-") THEN
         unit = output_unit
         ierr = 0
      ELSE
         OPEN (NEWUNIT=unit, FILE=TRIM(path), STATUS="replace", ACTION="write", IOSTAT=ierr)
      END IF
      IF (ierr /= 0) WRITE (message, "(A,1X,A)") "Cannot open output", TRIM(path)
   END SUBROUTINE open_output

   SUBROUTINE close_output(unit)
      INTEGER, INTENT(INOUT)                             :: unit

      IF (unit /= -1 .AND. unit /= output_unit) CLOSE (unit)
      unit = -1
   END SUBROUTINE close_output

   SUBROUTINE write_xyz_frame(unit, labels, values, comment)
      INTEGER, INTENT(IN)                                :: unit
      CHARACTER(LEN=*), INTENT(IN)                       :: labels(:)
      REAL(dp), INTENT(IN)                               :: values(:, :)
      CHARACTER(LEN=*), INTENT(IN)                       :: comment

      INTEGER                                            :: atom

      WRITE (unit, "(I0)") SIZE(labels)
      WRITE (unit, "(A)") TRIM(comment)
      DO atom = 1, SIZE(labels)
         WRITE (unit, "(A,3(1X,ES24.16))") TRIM(labels(atom)), values(:, atom)
      END DO
   END SUBROUTINE write_xyz_frame

END MODULE trajana_trajectory_io
