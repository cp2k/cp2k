!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

! Purpose: Output utilities.
! History: - Creation (06.10.1999, Matthias Krack)
!          - CP2K by JH 21.08.2000

MODULE output_utilities
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: print_message, print_warning
  
  INTEGER, PARAMETER :: out = 6
  
CONTAINS

!******************************************************************************

!   Purpose: Write a message string to the logical output unit lunit.
!   History: - Creation (28.08.1996, Matthias Krack)
!            - Translated to Fortran 90/95 (07.10.1999, Matthias Krack)
!            - CP2K by JH 21.08.2000

SUBROUTINE print_message ( message, lunit, declev, before, after )
  
!   ***************************************************************************

!   after  : Number of empty lines after the message.
!   before : Number of empty lines before the message.
!   declev : Decoration level (0,1,2, ... star lines).
!   lunit  : Number of the logical output unit.
!   message: String with the message text.

!   ***************************************************************************

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: message
  INTEGER, INTENT ( IN ), OPTIONAL :: after, before, declev, lunit
  
! Locals
  CHARACTER ( LEN = 12 ) :: fmtstr1
  CHARACTER ( LEN = 23 ) :: fmtstr2
  INTEGER :: blank_lines_after, blank_lines_before, decoration_level
  INTEGER :: i, ibreak, ipos1, ipos2, maxrowlen
  INTEGER :: msglen, nrow, output_unit, rowlen
  
!------------------------------------------------------------------------------
  
!   IF (my_pe /= io_pe) RETURN
  
  IF ( PRESENT ( after ) ) THEN
     blank_lines_after = MAX ( after, 0 )
  ELSE
     blank_lines_after = 1
  END IF
  
  IF ( PRESENT ( before ) ) THEN
     blank_lines_before = MAX ( before, 0 )
  ELSE
     blank_lines_before = 1
  END IF
  
  IF ( PRESENT ( declev ) ) THEN
     decoration_level = MAX ( declev, 0 )
  ELSE
     decoration_level = 0
  END IF
  
  IF ( PRESENT ( lunit ) ) THEN
     output_unit = lunit
     IF ( output_unit < 1 ) output_unit = out
  ELSE
     output_unit = out
  END IF
  
  IF ( decoration_level == 0) THEN
     rowlen = 78
  ELSE
     rowlen = 70
  END IF
  
  msglen = LEN_TRIM ( message )
  
!   *** Calculate number of rows ***
  
  nrow = msglen / ( rowlen + 1 ) + 1
  
!   *** Calculate appropriate row length ***
  
  rowlen = MIN ( msglen, msglen / nrow, rowlen )
  
!   *** Generate the blank lines before the message ***
  
  DO i = 1, blank_lines_before
     WRITE ( output_unit, ' ( )' )
  END DO
  
!   *** Scan for the longest row ***
  
  ipos1 = 1
  ipos2 = rowlen
  maxrowlen = 0
  
  DO
     IF ( ipos2 < msglen) THEN
        ibreak = ipos1 &
             + INDEX ( message ( ipos1:ipos2 ), " ", BACK = .TRUE. ) - 2
     ELSE
        ibreak = ipos2
     END IF
     
     maxrowlen = MAX ( maxrowlen, ibreak - ipos1 + 1)
     
     ipos1 = ibreak + 2
     ipos2 = MIN ( msglen, ipos1 + rowlen - 1)
     
!     *** When the last row is processed, exit loop ***
     
     IF ( ipos1 > msglen) EXIT
     
  END DO
  
!   *** Generate the first set of star rows ***
  
  IF ( decoration_level > 1) THEN
     fmtstr1 = '( T2, ("*") )'
     WRITE ( fmtstr1 ( 5:6 ), '( I2 )' ) maxrowlen + 8
     DO i = 1, decoration_level - 1
        WRITE ( output_unit, fmtstr1 )
     END DO
  END IF
  
!   *** Break long messages ***
  
  ipos1 = 1
  ipos2 = rowlen
  
  DO
     IF ( ipos2 < msglen ) THEN
        ibreak = ipos1 &
             + INDEX ( message ( ipos1:ipos2 ), " ", BACK = .TRUE. ) - 2
     ELSE
        ibreak = ipos2
     END IF
     
     IF ( decoration_level == 0 ) THEN
        WRITE ( output_unit, '( T2, A )' ) message ( ipos1:ibreak )
     ELSE IF ( decoration_level > 0 ) THEN
        fmtstr2 = '( T2, "*** ", A, X, "***" )'
        WRITE ( fmtstr2 ( 14:15 ), '( I2 )' ) ipos1 + maxrowlen - ibreak
        WRITE ( output_unit, fmtstr2 ) message ( ipos1:ibreak )
     END IF
     
     ipos1 = ibreak + 2
     ipos2 = MIN ( msglen, ipos1 + rowlen - 1 )
     
!     *** When the last row is processed, exit loop ***
     
     IF ( ipos1 > msglen ) EXIT
     
  END DO
  
!   *** Generate the second set star rows ***
  
  IF ( decoration_level > 1 ) THEN
     DO i = 1, decoration_level - 1
        WRITE ( output_unit, fmtstr1 )
     END DO
  END IF
  
!   *** Generate the blank lines after the message ***
  
  DO i = 1, blank_lines_after
     WRITE ( output_unit, '( )' )
  END DO
  
END SUBROUTINE print_message

!******************************************************************************

!   Purpose: Print a warning message.
!   History: - Creation (28.08.1996, Matthias Krack)
!            - Translated to Fortran 90/95 (07.10.1999, Matthias Krack)
!            - CP2K by JH 21.08.2000

SUBROUTINE print_warning ( routine, message, lunit )
  
!   ***************************************************************************
  
!   message: Message text.
!   routine: Name of the calling routine.
!   lunit  : Output unit
  
!   ***************************************************************************
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: message, routine
  INTEGER, INTENT ( IN ), OPTIONAL :: lunit
  
  INTEGER :: io
!------------------------------------------------------------------------------
  
  IF(PRESENT(lunit)) THEN
    io = lunit
  ELSE
    io = out
  ENDIF

  CALL print_message ( "WARNING in " // TRIM ( routine ), io, 1, 2, 0 )
  CALL print_message ( message, io, 1, 1, 1 )
  
END SUBROUTINE print_warning

!******************************************************************************

END MODULE output_utilities
