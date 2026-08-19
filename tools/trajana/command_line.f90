!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_command_line
   USE iso_fortran_env,                 ONLY: error_unit
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: command_name, get_option, has_flag, get_integer_option, get_real_option, fail, print_error

CONTAINS

   FUNCTION command_name() RESULT(command)
      CHARACTER(LEN=:), ALLOCATABLE                      :: command

      CHARACTER(LEN=1024)                                :: buffer

      buffer = ""
      IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, buffer)
      command = TRIM(buffer)
   END FUNCTION command_name

   SUBROUTINE get_option(name, value, found, default)
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT)         :: value
      LOGICAL, INTENT(OUT)                               :: found
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: default

      CHARACTER(LEN=4096)                                :: argument, following
      INTEGER                                            :: index

      found = .FALSE.
      value = ""
      IF (PRESENT(default)) value = TRIM(default)
      DO index = 1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(index, argument)
         IF (TRIM(argument) /= TRIM(name)) CYCLE
         IF (index == COMMAND_ARGUMENT_COUNT()) CALL fail("Missing value after "//TRIM(name))
         CALL GET_COMMAND_ARGUMENT(index + 1, following)
         value = TRIM(following)
         found = .TRUE.
      END DO
   END SUBROUTINE get_option

   LOGICAL FUNCTION has_flag(name)
      CHARACTER(LEN=*), INTENT(IN)                       :: name

      CHARACTER(LEN=4096)                                :: argument
      INTEGER                                            :: index

      has_flag = .FALSE.
      DO index = 1, COMMAND_ARGUMENT_COUNT()
         CALL GET_COMMAND_ARGUMENT(index, argument)
         IF (TRIM(argument) == TRIM(name)) THEN
            has_flag = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION has_flag

   SUBROUTINE get_integer_option(name, value, default)
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      INTEGER, INTENT(OUT)                               :: value
      INTEGER, INTENT(IN)                                :: default

      CHARACTER(LEN=:), ALLOCATABLE                      :: text
      INTEGER                                            :: ios
      LOGICAL                                            :: found

      CALL get_option(name, text, found)
      IF (.NOT. found) THEN
         value = default
         RETURN
      END IF
      READ (text, *, IOSTAT=ios) value
      IF (ios /= 0) CALL fail(TRIM(name)//" expects an integer")
   END SUBROUTINE get_integer_option

   SUBROUTINE get_real_option(name, value, default)
      CHARACTER(LEN=*), INTENT(IN)                       :: name
      REAL(dp), INTENT(OUT)                              :: value
      REAL(dp), INTENT(IN)                               :: default

      CHARACTER(LEN=:), ALLOCATABLE                      :: text
      INTEGER                                            :: ios
      LOGICAL                                            :: found

      CALL get_option(name, text, found)
      IF (.NOT. found) THEN
         value = default
         RETURN
      END IF
      READ (text, *, IOSTAT=ios) value
      IF (ios /= 0) CALL fail(TRIM(name)//" expects a number")
   END SUBROUTINE get_real_option

   SUBROUTINE print_error(message)
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      WRITE (error_unit, "(A)") "Error: "//TRIM(message)
   END SUBROUTINE print_error

   SUBROUTINE fail(message)
      CHARACTER(LEN=*), INTENT(IN)                       :: message

      CALL print_error(message)
      ERROR STOP 2
   END SUBROUTINE fail

END MODULE trajana_command_line
