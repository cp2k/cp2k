!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/stop_program [1.0] *
!!
!!   NAME
!!     stop_program
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE stop_program

  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_stop
  USE string_utilities, ONLY : compress
  
  IMPLICIT NONE
  
  INTERFACE stop_memory
    MODULE PROCEDURE stop_allocate,stop_deallocate
  END INTERFACE
  
  PRIVATE
  PUBLIC :: set_stop, stop_prg, stop_memory
  
  INTEGER :: iw
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** stop_program/set_stop [1.0] *
!!
!!   NAME
!!     set_stop
!!
!!   FUNCTION
!!     Initialise the output channel and the processor number
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_stop ( output_channel )
  
  IMPLICIT NONE
  
  INTEGER :: output_channel
  
  iw = output_channel
  
END SUBROUTINE set_stop

!!*****
!******************************************************************************
!!****** stop_program/stop_prg [1.0] *
!!
!!   NAME
!!     stop_prg
!!
!!   FUNCTION
!!     stop the program and write some final comments to the output
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE stop_prg ( routine, &
     comment1, comment2, comment3, comment4, comment5 )
  
  IMPLICIT NONE
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: &
       routine, comment1, comment2, comment3, comment4, comment5
  
! Locals
  INTEGER :: i
  CHARACTER ( LEN = 50 ) :: outline
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( routine ) ) THEN
     i = LEN_TRIM ( routine )
     outline = routine ( 1:i )
     WRITE ( iw, '( A, T31, A )' ) &
          " Last routine executed ", ADJUSTR ( outline )
  END IF
  
  IF ( PRESENT ( comment1 ) ) THEN
     i = LEN_TRIM ( comment1 )
     outline = comment1 ( 1:i )
     WRITE ( iw, '( A, T31, A )' ) " Comment ", ADJUSTR ( outline )
  END IF
  
  IF ( PRESENT ( comment2 ) ) THEN
     i = LEN_TRIM ( comment2 )
     outline = comment2 ( 1:i )
     WRITE ( iw, '( A, T31, A )' ) " Comment ", ADJUSTR ( outline )
  END IF
  
  IF ( PRESENT ( comment3 ) ) THEN
     i = LEN_TRIM ( comment3 )
     outline = comment3 ( 1:i )
     WRITE ( iw, '( A, T31, A )' ) " Comment ", ADJUSTR ( outline )
  END IF
  
  IF ( PRESENT ( comment4 ) ) THEN
     i = LEN_TRIM ( comment4 )
     outline = comment4 ( 1:i )
     WRITE ( iw, '( A, T31, A )' ) " Comment ", ADJUSTR ( outline )
  END IF
  
  IF ( PRESENT ( comment5 ) ) THEN
     i = LEN_TRIM ( comment5 )
     outline = comment5 ( 1:i )
     WRITE ( iw, '( A, T31, A )' ) " Comment ", ADJUSTR ( outline )
  END IF
  
  CALL mp_stop ( 999 )
  
END SUBROUTINE stop_prg

!!*****
!******************************************************************************
!!****** stop_program/stop_allocate [1.0] *
!!
!!   NAME
!!     stop_allocate
!!
!!   FUNCTION
!!     The memory allocation for an array failed. Print an error message
!!     and stop the program execution.
!!
!!   AUTHOR
!!     (12.10.1999, Matthias Krack)
!!
!!   MODIFICATION HISTORY
!!     CP2K by JGH 21.08.2000
!!   
!!   NOTES
!!     array:   Name of the array.
!!     memory:  Size of array in bytes.
!!     routine: Name of the calling routine.
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE stop_allocate ( routine, array, memory )
  
! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: array, routine
  INTEGER, INTENT ( IN ) :: memory
  
! Locals
  CHARACTER ( LEN = 50 ) :: message_1, message_2
  
!------------------------------------------------------------------------------
  
  IF (memory == 0) THEN
     WRITE (message_1,"(A)")&
          "Memory allocation for <"//TRIM(array)//"> failed"
     CALL stop_prg(routine,message_1)
  ELSE
     WRITE (message_1,"(A,I12,A)")&
          "Memory allocation for <"//TRIM(array)//"> failed" 
     WRITE (message_2,"(A,I12,A)")&
          "Requested memory size is ",memory/1024," Kbytes"
     CALL compress(message_2)
     CALL stop_prg(routine,message_1,message_2)
  END IF
  
END SUBROUTINE stop_allocate

!!*****
!******************************************************************************
!!****** stop_program/stop_deallocate [1.0] *
!!
!!   NAME
!!     stop_deallocate
!!
!!   FUNCTION
!!     The memory deallocation for an array failed. Print an error
!!     message and stop the program execution.
!!
!!   AUTHOR
!!     (20.10.1999, Matthias Krack)
!!
!!   MODIFICATION HISTORY
!!     CP2K by JGH 21.08.2000
!!
!!   NOTES
!!     array:   Name of the array.
!!     routine: Name of the calling routine.
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE stop_deallocate ( routine, array )

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: array, routine
  
! Locals
  CHARACTER ( LEN = 50 ) :: message
  
!------------------------------------------------------------------------------
  
  WRITE ( message, '( A )' ) &
       "Memory allocation for <" // TRIM ( array ) // "> failed"
  CALL stop_prg ( routine, message )
  
END SUBROUTINE stop_deallocate

!!*****
!******************************************************************************

END MODULE stop_program

!******************************************************************************
