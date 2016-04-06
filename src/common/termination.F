!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Perform an abnormal program termination.
!> \note These routines are low-level and thus provide also an error recovery
!>       when dependencies do not allow the use of the error logger. Only
!>       the master (root) process will dump, if para_env is available and
!>       properly specified. Otherwise (without any information about the
!>       parallel environment) most likely more than one process or even all
!>       processes will send their error dump to the default output unit.
!> \par History
!>      - Routine external_control moved to a separate module
!> \author Matthias Krack (12.02.2001)
! **************************************************************************************************
MODULE termination
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              int_8
   USE machine,                         ONLY: m_memory
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'termination'

   PUBLIC :: print_message
   PUBLIC :: stop_memory

CONTAINS

! **************************************************************************************************
!> \brief Perform a basic blocking of the text in message and print it
!>        optionally decorated with a frame of stars as defined by declev.
!> \param message ...
!> \param output_unit ...
!> \param declev ...
!> \param before ...
!> \param after ...
!> \date 28.08.1996
!> \par History
!>      - Translated to Fortran 90/95 (07.10.1999, Matthias Krack)
!>      - CP2K by JH 21.08.2000
!>      - Bugs in the dynamic format generation removed (09.02.2001, MK)
!>      - Revised (26.01.2011,MK)
!> \author Matthias Krack (MK)
!> \note
!>       after      : Number of empty lines after the message.
!>       before     : Number of empty lines before the message.
!>       declev     : Decoration level (0,1,2, ... star lines).
!>       message    : String with the message text.
!>       output_unit: Logical unit number of output unit.
! **************************************************************************************************
   SUBROUTINE print_message(message, output_unit, declev, before, after)

      CHARACTER(LEN=*), INTENT(IN)                       :: message
      INTEGER, INTENT(IN)                                :: output_unit
      INTEGER, INTENT(IN), OPTIONAL                      :: declev, before, after

      INTEGER                                            :: blank_lines_after, blank_lines_before, &
                                                            decoration_level, i, ibreak, ipos1, &
                                                            ipos2, maxrowlen, msglen, nrow, rowlen

      IF (PRESENT(after)) THEN
         blank_lines_after = MAX(after, 0)
      ELSE
         blank_lines_after = 1
      END IF

      IF (PRESENT(before)) THEN
         blank_lines_before = MAX(before, 0)
      ELSE
         blank_lines_before = 1
      END IF

      IF (PRESENT(declev)) THEN
         decoration_level = MAX(declev, 0)
      ELSE
         decoration_level = 0
      END IF

      IF (decoration_level == 0) THEN
         rowlen = 78
      ELSE
         rowlen = 70
      END IF

      msglen = LEN_TRIM(message)

      ! Calculate number of rows

      nrow = msglen/(rowlen+1)+1

      ! Calculate appropriate row length

      rowlen = MIN(msglen, rowlen)

      ! Generate the blank lines before the message

      DO i = 1, blank_lines_before
         WRITE (UNIT=output_unit, FMT="(A)") ""
      END DO

      ! Scan for the longest row

      ipos1 = 1
      ipos2 = rowlen
      maxrowlen = 0

      DO
         IF (ipos2 < msglen) THEN
            i = INDEX(message(ipos1:ipos2), " ", BACK=.TRUE.)
            IF (i == 0) THEN
               ibreak = ipos2
            ELSE
               ibreak = ipos1+i-2
            END IF
         ELSE
            ibreak = ipos2
         END IF

         maxrowlen = MAX(maxrowlen, ibreak-ipos1+1)

         ipos1 = ibreak+2
         ipos2 = MIN(msglen, ipos1+rowlen-1)

         ! When the last row is processed, exit loop

         IF (ipos1 > msglen) EXIT

      END DO

      ! Generate the first set of star rows

      IF (decoration_level > 1) THEN
         DO i = 1, decoration_level-1
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", maxrowlen+8)
         END DO
      END IF

      ! Break long messages

      ipos1 = 1
      ipos2 = rowlen

      DO
         IF (ipos2 < msglen) THEN
            i = INDEX(message(ipos1:ipos2), " ", BACK=.TRUE.)
            IF (i == 0) THEN
               ibreak = ipos2
            ELSE
               ibreak = ipos1+i-2
            END IF
         ELSE
            ibreak = ipos2
         END IF

         IF (decoration_level == 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") message(ipos1:ibreak)
         ELSE IF (decoration_level > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "*** "//message(ipos1:ibreak)//REPEAT(" ", ipos1+maxrowlen-ibreak)//"***"
         END IF

         ipos1 = ibreak+2
         ipos2 = MIN(msglen, ipos1+rowlen-1)

         ! When the last row is processed, exit loop

         IF (ipos1 > msglen) EXIT
      END DO

      ! Generate the second set star rows

      IF (decoration_level > 1) THEN
         DO i = 1, decoration_level-1
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", maxrowlen+8)
         END DO
      END IF

      ! Generate the blank lines after the message

      DO i = 1, blank_lines_after
         WRITE (UNIT=output_unit, FMT="(A)") ""
      END DO

   END SUBROUTINE print_message

! **************************************************************************************************
!> \brief The memory allocation for a data object failed. Print an error
!>        message and stop the program execution.
!> \param routineN ...
!> \param moduleN ...
!> \param line_number ...
!> \param object ...
!> \param memory ...
!> \par History
!>      - CP2K by JGH 21.08.2000
!>      - Revised (16.02.2011,MK)
!> \author Matthias Krack (12.10.1999,MK)
!> \note
!>      routineN   : Name of the calling routine
!>      moduleN    : Name of the calling module including the calling routine
!>      line_number: Code line number in the module
!>      object     : Name of the data object for which the allocation failed
!>      memory     : Size of the data object in bytes. If zero no size is
!>                   printed.
! **************************************************************************************************
   SUBROUTINE stop_memory(routineN, moduleN, line_number, object, memory)

      CHARACTER(LEN=*), INTENT(IN)                       :: routineN, moduleN
      INTEGER, INTENT(IN)                                :: line_number
      CHARACTER(LEN=*), INTENT(IN)                       :: object
      INTEGER, INTENT(IN)                                :: memory

      CHARACTER(LEN=default_path_length)                 :: message
      CHARACTER(len=default_string_length)               :: location
      INTEGER(KIND=int_8)                                :: total_memory

      CALL m_memory(total_memory)
      total_memory = total_memory/1024

      IF (memory == 0) THEN
         message = "The memory allocation for the data object <"//TRIM(object)// &
                   "> failed."
      ELSE
         WRITE (message, "(A,I0,A)") &
            "The memory allocation for the data object <"//TRIM(object)// &
            "> failed. The requested memory size is ", memory/1024, " KB."
      END IF

      IF (total_memory > 0) THEN
         WRITE (message, "(A,I0,A)") &
            TRIM(message)//" The current total memory allocated by CP2K is ", &
            total_memory, " KB."
      END IF

      WRITE (location, FMT='(A,I6)') TRIM(routineN)// &
         " (MODULE "//TRIM(moduleN)//")"// &
         ":", line_number

      CALL cp_abort(location, message)

   END SUBROUTINE stop_memory

END MODULE termination
