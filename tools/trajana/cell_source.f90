!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_cell_source
   USE trajana_text_utils,              ONLY: parse_logical_axes
   USE trajana_trajectory_io,           ONLY: cell_reader_type,&
                                              parse_constant_cell
   USE trajana_trajectory_types,        ONLY: cell_type,&
                                              frame_type

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: cell_source_type
      LOGICAL :: use_constant = .FALSE.
      LOGICAL :: use_file = .FALSE.
      TYPE(cell_type) :: constant_cell = cell_type()
      TYPE(cell_reader_type) :: reader = cell_reader_type()
      LOGICAL :: periodic(3) = .TRUE.
   CONTAINS
      PROCEDURE :: configure => source_configure
      PROCEDURE :: attach => source_attach
      PROCEDURE :: close => source_close
   END TYPE cell_source_type

CONTAINS

   SUBROUTINE source_configure(source, cell_text, cell_path, periodic_text, ierr, message)
      CLASS(cell_source_type), INTENT(INOUT) :: source
      CHARACTER(LEN=*), INTENT(IN) :: cell_text, cell_path, periodic_text
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(OUT) :: message
      LOGICAL :: ok

      ierr = 0
      message = ""
      IF (LEN_TRIM(cell_text) > 0 .AND. LEN_TRIM(cell_path) > 0) THEN
         ierr = 1
         message = "Use either --cell or --cell-file, not both"
         RETURN
      END IF
      CALL parse_logical_axes(periodic_text, source%periodic, ok)
      IF (.NOT. ok) THEN
         ierr = 1
         message = "--periodic expects XYZ, XY, XZ, YZ, X, Y, Z, or NONE"
         RETURN
      END IF

      IF (LEN_TRIM(cell_text) > 0) THEN
         CALL parse_constant_cell(cell_text, source%constant_cell, ierr, message)
         IF (ierr /= 0) RETURN
         source%constant_cell%periodic = source%periodic
         source%use_constant = .TRUE.
      ELSE IF (LEN_TRIM(cell_path) > 0) THEN
         CALL source%reader%open_file(cell_path, ierr, message)
         IF (ierr /= 0) RETURN
         source%use_file = .TRUE.
      END IF
   END SUBROUTINE source_configure

   SUBROUTINE source_attach(source, frame, ierr, message)
      CLASS(cell_source_type), INTENT(INOUT) :: source
      TYPE(frame_type), INTENT(INOUT) :: frame
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(OUT) :: message
      LOGICAL :: eof
      TYPE(cell_type) :: current

      ierr = 0
      message = ""
      IF (source%use_constant) THEN
         frame%cell = source%constant_cell
      ELSE IF (source%use_file) THEN
         CALL source%reader%read_cell(current, eof, ierr, message)
         IF (ierr /= 0) RETURN
         IF (eof) THEN
            ierr = 1
            message = "The cell file has fewer frames than the trajectory"
            RETURN
         END IF
         current%periodic = source%periodic
         frame%cell = current
      ELSE IF (frame%cell%valid) THEN
         frame%cell%periodic = source%periodic
      END IF
   END SUBROUTINE source_attach

   SUBROUTINE source_close(source)
      CLASS(cell_source_type), INTENT(INOUT) :: source

      IF (source%use_file) CALL source%reader%close_file()
   END SUBROUTINE source_close

END MODULE trajana_cell_source
