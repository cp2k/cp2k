!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_selections
   USE trajana_text_utils,              ONLY: lower_case,&
                                              next_token,&
                                              starts_with

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: build_selection

CONTAINS

   SUBROUTINE build_selection(specification, labels, selected, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: specification, labels(:)
      LOGICAL, ALLOCATABLE, INTENT(OUT)                  :: selected(:)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=:), ALLOCATABLE                      :: body, spec

      ALLOCATE (selected(SIZE(labels)))
      selected = .FALSE.
      ierr = 0
      message = ""
      spec = lower_case(TRIM(specification))
      IF (spec == "" .OR. spec == "all") THEN
         selected = .TRUE.
      ELSE IF (starts_with(spec, "name:")) THEN
         body = specification(6:)
         CALL select_names(body, labels, selected)
      ELSE IF (starts_with(spec, "index:")) THEN
         body = specification(7:)
         CALL select_indices(body, selected, ierr, message)
      ELSE
         ierr = 1
         message = "Selection must be all, name:LABEL,..., or index:N,M-P,..."
      END IF
      IF (ierr == 0 .AND. COUNT(selected) == 0) THEN
         ierr = 1
         message = "Selection does not match any atoms"
      END IF
   END SUBROUTINE build_selection

   SUBROUTINE select_names(body, labels, selected)
      CHARACTER(LEN=*), INTENT(IN)                       :: body, labels(:)
      LOGICAL, INTENT(INOUT)                             :: selected(:)

      CHARACTER(LEN=128)                                 :: token
      INTEGER                                            :: atom, position
      LOGICAL                                            :: found

      position = 1
      DO
         CALL next_token(body, position, token, found)
         IF (.NOT. found) EXIT
         DO atom = 1, SIZE(labels)
            IF (lower_case(TRIM(labels(atom))) == lower_case(TRIM(token))) selected(atom) = .TRUE.
         END DO
      END DO
   END SUBROUTINE select_names

   SUBROUTINE select_indices(body, selected, ierr, message)
      CHARACTER(LEN=*), INTENT(IN)                       :: body
      LOGICAL, INTENT(INOUT)                             :: selected(:)
      INTEGER, INTENT(OUT)                               :: ierr
      CHARACTER(LEN=*), INTENT(OUT)                      :: message

      CHARACTER(LEN=128)                                 :: token
      INTEGER                                            :: atom, dash, first, ios, last, position
      LOGICAL                                            :: found

      ierr = 0
      message = ""
      position = 1
      DO
         CALL next_token(body, position, token, found)
         IF (.NOT. found) EXIT
         dash = INDEX(TRIM(token), "-")
         IF (dash == 0) THEN
            READ (token, *, IOSTAT=ios) first
            last = first
         ELSE
            READ (token(:dash - 1), *, IOSTAT=ios) first
            IF (ios == 0) READ (token(dash + 1:), *, IOSTAT=ios) last
         END IF
         IF (ios /= 0 .OR. first < 1 .OR. last < first .OR. last > SIZE(selected)) THEN
            ierr = 1
            message = "Invalid atom-index selection: "//TRIM(token)
            RETURN
         END IF
         DO atom = first, last
            selected(atom) = .TRUE.
         END DO
      END DO
   END SUBROUTINE select_indices

END MODULE trajana_selections
