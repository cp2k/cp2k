!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
!>      - Delete stop_memory routine, rename module
!> \author Matthias Krack (12.02.2001)
! **************************************************************************************************
MODULE print_messages
#include "../base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'print_messages'

   PUBLIC :: print_message

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

      nrow = msglen/(rowlen + 1) + 1

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
               ibreak = ipos1 + i - 2
            END IF
         ELSE
            ibreak = ipos2
         END IF

         maxrowlen = MAX(maxrowlen, ibreak - ipos1 + 1)

         ipos1 = ibreak + 2
         ipos2 = MIN(msglen, ipos1 + rowlen - 1)

         ! When the last row is processed, exit loop

         IF (ipos1 > msglen) EXIT

      END DO

      ! Generate the first set of star rows

      IF (decoration_level > 1) THEN
         DO i = 1, decoration_level - 1
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", maxrowlen + 8)
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
               ibreak = ipos1 + i - 2
            END IF
         ELSE
            ibreak = ipos2
         END IF

         IF (decoration_level == 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") message(ipos1:ibreak)
         ELSE IF (decoration_level > 0) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "*** "//message(ipos1:ibreak)//REPEAT(" ", ipos1 + maxrowlen - ibreak)//"***"
         END IF

         ipos1 = ibreak + 2
         ipos2 = MIN(msglen, ipos1 + rowlen - 1)

         ! When the last row is processed, exit loop

         IF (ipos1 > msglen) EXIT
      END DO

      ! Generate the second set star rows

      IF (decoration_level > 1) THEN
         DO i = 1, decoration_level - 1
            WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("*", maxrowlen + 8)
         END DO
      END IF

      ! Generate the blank lines after the message

      DO i = 1, blank_lines_after
         WRITE (UNIT=output_unit, FMT="(A)") ""
      END DO

   END SUBROUTINE print_message

END MODULE print_messages
