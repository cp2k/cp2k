!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_text_utils
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: lower_case, next_token, parse_logical_axes, starts_with

CONTAINS

   PURE FUNCTION lower_case(text) RESULT(lower)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      CHARACTER(LEN=LEN(text))                           :: lower

      INTEGER                                            :: code, i

      lower = text
      DO i = 1, LEN(text)
         code = IACHAR(text(i:i))
         IF (code >= IACHAR("A") .AND. code <= IACHAR("Z")) lower(i:i) = ACHAR(code + 32)
      END DO
   END FUNCTION lower_case

   PURE LOGICAL FUNCTION starts_with(text, prefix)
      CHARACTER(LEN=*), INTENT(IN)                       :: text, prefix

      starts_with = LEN_TRIM(text) >= LEN_TRIM(prefix)
      IF (starts_with) starts_with = text(1:LEN_TRIM(prefix)) == prefix(1:LEN_TRIM(prefix))
   END FUNCTION starts_with

   SUBROUTINE next_token(line, position, token, found)
      CHARACTER(LEN=*), INTENT(IN)                       :: line
      INTEGER, INTENT(INOUT)                             :: position
      CHARACTER(LEN=*), INTENT(OUT)                      :: token
      LOGICAL, INTENT(OUT)                               :: found

      INTEGER                                            :: first, last

      token = ""
      first = position
      DO WHILE (first <= LEN_TRIM(line))
         IF (line(first:first) /= " " .AND. line(first:first) /= ACHAR(9) .AND. line(first:first) /= ",") EXIT
         first = first + 1
      END DO
      IF (first > LEN_TRIM(line)) THEN
         found = .FALSE.
         position = LEN_TRIM(line) + 1
         RETURN
      END IF
      IF (line(first:first) == "#") THEN
         found = .FALSE.
         position = LEN_TRIM(line) + 1
         RETURN
      END IF

      last = first
      DO WHILE (last <= LEN_TRIM(line))
         IF (line(last:last) == " " .OR. line(last:last) == ACHAR(9) .OR. line(last:last) == "," .OR. &
             line(last:last) == "#") EXIT
         last = last + 1
      END DO
      token = line(first:last - 1)
      position = last + 1
      found = .TRUE.
   END SUBROUTINE next_token

   SUBROUTINE parse_logical_axes(text, periodic, ok)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      LOGICAL, INTENT(OUT)                               :: periodic(3), ok

      CHARACTER(LEN=:), ALLOCATABLE                      :: value

      value = lower_case(TRIM(text))
      periodic = .FALSE.
      periodic(1) = INDEX(value, "x") > 0
      periodic(2) = INDEX(value, "y") > 0
      periodic(3) = INDEX(value, "z") > 0
      ok = LEN_TRIM(value) > 0 .AND. VERIFY(value, "xyz") == 0
      IF (value == "none") THEN
         periodic = .FALSE.
         ok = .TRUE.
      END IF
   END SUBROUTINE parse_logical_axes

END MODULE trajana_text_utils
