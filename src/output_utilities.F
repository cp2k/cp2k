!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/output_utilities [1.0] *
!!
!!   NAME
!!     output_utilities
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   CREATION DATE
!!     06.10.1999
!!
!!   MODIFICATION HISTORY
!!     - CP2K Version by JGH 21.08.2000
!!     - Bugs in the dynamic format generation removed (09.02.2001, MK)
!!     - Default output unit definition removed (14.02.2001, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE output_utilities

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: default_output_unit = 6

  PUBLIC :: print_message,&
            print_warning

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** output_utilities/print_message [1.0] *
!!
!!   NAME
!!     print_message
!!
!!   SYNOPSIS
!!     Subroutine print_message(message, output_unit, declev, before,&
!!         after)
!!       Character(Len=*), Intent (IN):: message
!!       Integer, Intent (IN):: output_unit
!!       Integer, Optional, Intent (IN):: after, before, declev
!!     End Subroutine print_message
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   CREATION DATE
!!     28.08.1996
!!
!!   MODIFICATION HISTORY
!!     - Translated to Fortran 90/95 (07.10.1999, Matthias Krack)
!!     - CP2K by JH 21.08.2000
!!     - Bugs in the dynamic format generation removed (09.02.2001, MK)
!!
!!    NOTES
!!      after      : Number of empty lines after the message.
!!      before     : Number of empty lines before the message.
!!      declev     : Decoration level (0,1,2, ... star lines).
!!      message    : String with the message text.
!!      output_unit: Logical unit number of output unit.
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE print_message(message,output_unit,declev,before,after)

    CHARACTER(LEN=*), INTENT(IN)             :: message
    INTEGER, INTENT(IN)                      :: output_unit
    INTEGER, INTENT(IN), OPTIONAL            :: declev, before, after

    CHARACTER(LEN=12)                        :: fmtstr1
    CHARACTER(LEN=23)                        :: fmtstr2
    INTEGER :: blank_lines_after, blank_lines_before, decoration_level, i, &
      ibreak, ipos1, ipos2, maxrowlen, msglen, nrow, rowlen

!   ---------------------------------------------------------------------------

    IF (PRESENT(after)) THEN
      blank_lines_after = MAX(after,0)
    ELSE
      blank_lines_after = 1
    END IF

    IF (PRESENT(before)) THEN
      blank_lines_before = MAX(before,0)
    ELSE
      blank_lines_before = 1
    END IF

    IF (PRESENT(declev)) THEN
      decoration_level = MAX(declev,0)
    ELSE
      decoration_level = 0
    END IF

    IF (decoration_level == 0) THEN
      rowlen = 78
    ELSE
      rowlen = 70
    END IF

    msglen = LEN_TRIM(message)

!   *** Calculate number of rows ***

    nrow = msglen/(rowlen + 1) + 1

!   *** Calculate appropriate row length ***

    rowlen = MIN(msglen,rowlen)

!   *** Generate the blank lines before the message ***

    DO i=1,blank_lines_before
      WRITE (UNIT=output_unit,FMT="(A)") ""
    END DO

!   *** Scan for the longest row ***

    ipos1 = 1
    ipos2 = rowlen
    maxrowlen = 0

    DO
      IF (ipos2 < msglen) THEN
         i=INDEX(message(ipos1:ipos2)," ",BACK=.TRUE.)
         IF (i == 0) THEN
            ibreak = ipos2
         ELSE
            ibreak = ipos1 + i - 2
         END IF
      ELSE
        ibreak = ipos2
      END IF

      maxrowlen = MAX(maxrowlen,ibreak - ipos1 + 1)

      ipos1 = ibreak + 2
      ipos2 = MIN(msglen,ipos1 + rowlen - 1)

!     *** When the last row is processed, exit loop ***

      IF (ipos1 > msglen) EXIT

    END DO

!   *** Generate the first set of star rows ***

    IF (decoration_level > 1) THEN
      fmtstr1 = "(T2,  (""*""))"
      WRITE (UNIT=fmtstr1(5:6),FMT="(I2)") maxrowlen + 8
      DO i=1,decoration_level-1
        WRITE (UNIT=output_unit,FMT=fmtstr1)
      END DO
    END IF

!   *** Break long messages ***

    ipos1 = 1
    ipos2 = rowlen

    DO

      IF (ipos2 < msglen) THEN
         i=INDEX(message(ipos1:ipos2)," ",BACK=.TRUE.)
         IF (i == 0) THEN
            ibreak = ipos2
         ELSE
            ibreak = ipos1 + i - 2
         END IF
      ELSE
        ibreak = ipos2
      END IF

      IF (decoration_level == 0) THEN
        WRITE (UNIT=output_unit,FMT="(T2,A)") message(ipos1:ibreak)
      ELSE IF (decoration_level > 0) THEN
        fmtstr2 = "(T2,""*** "",A,  X,""***"")"
        WRITE (UNIT=fmtstr2(14:15),FMT="(I2)") ipos1 + maxrowlen - ibreak
        WRITE (UNIT=output_unit,FMT=fmtstr2) message(ipos1:ibreak)
      END IF

      ipos1 = ibreak + 2
      ipos2 = MIN(msglen,ipos1 + rowlen - 1)

!     *** When the last row is processed, exit loop ***

      IF (ipos1 > msglen) EXIT

    END DO

!   *** Generate the second set star rows ***

    IF (decoration_level > 1) THEN
      DO i=1,decoration_level-1
        WRITE (UNIT=output_unit,FMT=fmtstr1)
      END DO
    END IF

!   *** Generate the blank lines after the message ***

    DO i=1,blank_lines_after
      WRITE (UNIT=output_unit,FMT="(A)") ""
    END DO

  END SUBROUTINE print_message

!!*****
!******************************************************************************
!!****** output_utilities/print_warning [1.0] *
!!
!!   NAME
!!     print_warning
!!
!!   SYNOPSIS
!!     Subroutine print_warning(routine, message, output_unit)
!!       Character(Len=*), Intent (IN):: message, routine
!!       Integer, Optional, Intent (IN):: output_unit
!!     End Subroutine print_warning
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   CREATION DATE
!!     28.08.1996
!!
!!   MODIFICATION HISTORY
!!     - Translated to Fortran 90/95 (07.10.1999, Matthias Krack)
!!     - CP2K by JH 21.08.2000
!!
!!    NOTES
!!      message    : Message text.
!!      routine    : Name of the calling routine.
!!      output_unit: Logical unit number of the output unit
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE print_warning(routine,message,output_unit)

    CHARACTER(LEN=*), INTENT(IN)             :: routine, message
    INTEGER, INTENT(IN), OPTIONAL            :: output_unit

    INTEGER                                  :: lunit

!   ---------------------------------------------------------------------------

    IF (PRESENT(output_unit)) THEN
      lunit = output_unit
    ELSE
      lunit = default_output_unit
    END IF

    CALL print_message("WARNING in "//TRIM(routine),lunit,1,2,0)
    CALL print_message(message,lunit,1,1,1)

  END SUBROUTINE print_warning

!!*****
!******************************************************************************

END MODULE output_utilities

!******************************************************************************
