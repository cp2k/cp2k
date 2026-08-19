!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_groups
   USE trajana_kinds,                   ONLY: dp
   USE trajana_text_utils,              ONLY: lower_case,&
                                              next_token

   IMPLICIT NONE
   PRIVATE

   TYPE, PUBLIC :: group_type
      CHARACTER(LEN=32) :: label = "X"
      INTEGER, ALLOCATABLE :: atom(:)
   END TYPE group_type

   TYPE, PUBLIC :: mass_table_type
      CHARACTER(LEN=32), ALLOCATABLE :: label(:)
      REAL(dp), ALLOCATABLE :: mass(:)
   CONTAINS
      PROCEDURE :: lookup => mass_lookup
   END TYPE mass_table_type

   PUBLIC :: read_groups, read_masses, validate_groups

CONTAINS

   SUBROUTINE read_groups(path, groups, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      TYPE(group_type), ALLOCATABLE, INTENT(OUT)         :: groups(:)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=128)                                 :: token
      CHARACTER(LEN=4096)                                :: line
      INTEGER                                            :: count, ios, item, ntokens, position, unit
      LOGICAL                                            :: found
      TYPE(group_type), ALLOCATABLE                      :: work(:)

      ierr = 0
      message = ""
      count = 0
      ALLOCATE (groups(0))
      OPEN (NEWUNIT=unit, FILE=TRIM(path), STATUS="old", ACTION="read", IOSTAT=ios)
      IF (ios /= 0) THEN
         ierr = ios
         message = "Cannot open group file "//TRIM(path)
         RETURN
      END IF

      DO
         READ (unit, "(A)", IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) EXIT
         IF (ios /= 0) THEN
            ierr = ios
            message = "Error while reading group file"
            EXIT
         END IF
         line = ADJUSTL(line)
         IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE

         position = 1
         CALL next_token(line, position, token, found)
         IF (.NOT. found) CYCLE
         ntokens = 0
         DO
            CALL next_token(line, position, token, found)
            IF (.NOT. found) EXIT
            ntokens = ntokens + 1
         END DO
         IF (ntokens == 0) THEN
            ierr = 1
            message = "Every group needs a label and at least one atom index"
            EXIT
         END IF

         ALLOCATE (work(count + 1))
         IF (count > 0) work(:count) = groups
         position = 1
         CALL next_token(line, position, token, found)
         work(count + 1)%label = TRIM(token)
         ALLOCATE (work(count + 1)%atom(ntokens))
         DO item = 1, ntokens
            CALL next_token(line, position, token, found)
            READ (token, *, IOSTAT=ios) work(count + 1)%atom(item)
            IF (ios /= 0) THEN
               ierr = 1
               message = "Invalid atom index in group file: "//TRIM(token)
               EXIT
            END IF
         END DO
         IF (ierr /= 0) EXIT
         CALL MOVE_ALLOC(work, groups)
         count = count + 1
      END DO
      CLOSE (unit)
      IF (ierr == 0 .AND. count == 0) THEN
         ierr = 1
         message = "The group file contains no groups"
      END IF
   END SUBROUTINE read_groups

   SUBROUTINE validate_groups(groups, natoms, minimum_size, ierr, message)
      TYPE(group_type), INTENT(IN)                       :: groups(:)
      INTEGER, INTENT(IN)                                :: natoms, minimum_size
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      INTEGER                                            :: group, item

      ierr = 0
      message = ""
      DO group = 1, SIZE(groups)
         IF (SIZE(groups(group)%atom) < minimum_size) THEN
            ierr = 1
            WRITE (message, "(A,I0,A,I0,A)") "Group ", group, " needs at least ", minimum_size, " atoms"
            RETURN
         END IF
         DO item = 1, SIZE(groups(group)%atom)
            IF (groups(group)%atom(item) < 1 .OR. groups(group)%atom(item) > natoms) THEN
               ierr = 1
               WRITE (message, "(A,I0,A,I0)") "Invalid atom index in group ", group, ": ", groups(group)%atom(item)
               RETURN
            END IF
         END DO
      END DO
   END SUBROUTINE validate_groups

   SUBROUTINE read_masses(path, table, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: path
      TYPE(mass_table_type), INTENT(OUT)                 :: table
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=32)                                  :: label
      CHARACTER(LEN=32), ALLOCATABLE                     :: labels(:)
      CHARACTER(LEN=4096)                                :: line
      INTEGER                                            :: count, ios, unit
      REAL(dp)                                           :: mass
      REAL(dp), ALLOCATABLE                              :: masses(:)

      ierr = 0
      message = ""
      count = 0
      ALLOCATE (table%label(0), table%mass(0))
      OPEN (NEWUNIT=unit, FILE=TRIM(path), STATUS="old", ACTION="read", IOSTAT=ios)
      IF (ios /= 0) THEN
         ierr = ios
         message = "Cannot open mass file "//TRIM(path)
         RETURN
      END IF
      DO
         READ (unit, "(A)", IOSTAT=ios) line
         IF (IS_IOSTAT_END(ios)) EXIT
         IF (ios /= 0) THEN
            ierr = ios
            message = "Error while reading mass file"
            EXIT
         END IF
         line = ADJUSTL(line)
         IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") CYCLE
         READ (line, *, IOSTAT=ios) label, mass
         IF (ios /= 0 .OR. mass <= 0.0_dp) THEN
            ierr = 1
            message = "Mass-file lines must contain LABEL MASS"
            EXIT
         END IF
         ALLOCATE (labels(count + 1), masses(count + 1))
         IF (count > 0) THEN
            labels(:count) = table%label
            masses(:count) = table%mass
         END IF
         labels(count + 1) = label
         masses(count + 1) = mass
         CALL MOVE_ALLOC(labels, table%label)
         CALL MOVE_ALLOC(masses, table%mass)
         count = count + 1
      END DO
      CLOSE (unit)
      IF (ierr == 0 .AND. count == 0) THEN
         ierr = 1
         message = "The mass file contains no masses"
      END IF
   END SUBROUTINE read_masses

   SUBROUTINE mass_lookup(table, label, mass, found)
      CLASS(mass_table_type), INTENT(IN) :: table
      CHARACTER(LEN=*), INTENT(IN) :: label
      REAL(dp), INTENT(OUT) :: mass
      LOGICAL, INTENT(OUT) :: found
      INTEGER :: item

      found = .FALSE.
      mass = 0.0_dp
      DO item = 1, SIZE(table%label)
         IF (lower_case(TRIM(table%label(item))) == lower_case(TRIM(label))) THEN
            mass = table%mass(item)
            found = .TRUE.
            RETURN
         END IF
      END DO
   END SUBROUTINE mass_lookup

END MODULE trajana_groups
