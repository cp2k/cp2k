!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/termination [1.0] *
!!
!!   NAME
!!     termination
!!
!!   FUNCTION
!!     Perform an abnormal program termination.
!!
!!   AUTHOR
!!     Matthias Krack (12.02.2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   INDEX
!!     SUBROUTINE get_error_unit(lunit)
!!     SUBROUTINE set_error_unit(lunit)
!!     SUBROUTINE stop_allocate(routine,array,memory)
!!     SUBROUTINE stop_deallocate(routine,array)
!!     SUBROUTINE stop_program(routine,message)
!!
!!   SOURCE
!******************************************************************************

MODULE termination

  USE message_passing,  ONLY: mp_stop
  USE output_utilities, ONLY: print_message
  USE string_utilities, ONLY: compress

  IMPLICIT NONE

  PRIVATE

  INTEGER :: output_unit

  PUBLIC :: get_error_unit,&
            set_error_unit,&
            stop_memory,&
            stop_program

  INTERFACE stop_memory
    MODULE PROCEDURE stop_allocate,stop_deallocate
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** termination/get_error_unit [1.0] *
!!
!!   NAME
!!     get_error_unit
!!
!!   FUNCTION
!!     Get the output unit number for error messages.
!!
!!   AUTHOR
!!     MK
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE get_error_unit(lunit)

    INTEGER, INTENT(OUT) :: lunit

!   ---------------------------------------------------------------------------

    lunit = output_unit

  END SUBROUTINE get_error_unit

!!*****
!******************************************************************************
!!****** termination/set_error_unit [1.0] *
!!
!!   NAME
!!     set_error_unit
!!
!!   FUNCTION
!!     Initialise the output unit number for error messages.
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE set_error_unit(lunit)

    INTEGER, INTENT(IN) :: lunit

!   ---------------------------------------------------------------------------

    output_unit = lunit

  END SUBROUTINE set_error_unit

!!*****
!******************************************************************************
!!****** termination/stop_allocate [1.0] *
!!
!!   NAME
!!     stop_allocate
!!
!!   FUNCTION
!!     The memory allocation for an array failed. Print an error message and
!!     stop the program execution.
!!
!!   AUTHOR
!!     Matthias Krack (12.10.1999)
!!
!!   MODIFICATION HISTORY
!!     CP2K by JGH 21.08.2000
!!
!!   NOTES
!!     array  : Name of the array.
!!     memory : Size of array in bytes.
!!     routine: Name of the calling routine.
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE stop_allocate(routine,array,memory)

    CHARACTER(LEN=*), INTENT(IN) :: array,routine
    INTEGER, INTENT(IN)          :: memory

!   *** Local variables ***

    CHARACTER(LEN=240) :: message

!   ---------------------------------------------------------------------------

    IF (memory == 0) THEN
      WRITE (message,"(A)")&
        "The memory allocation for the data object <"//TRIM(array)//"> failed"
    ELSE
      WRITE (message,"(A,I12,A)")&
        "The memory allocation for the data object <"//TRIM(array)//&
        "> failed. The requested memory size is ",memory/1024," Kbytes"
    END IF

    CALL compress(message)

    CALL stop_program(routine,message)

  END SUBROUTINE stop_allocate

!!*****
!******************************************************************************
!!****** termination/stop_deallocate [1.0] *
!!
!!   NAME
!!     stop_deallocate
!!
!!   FUNCTION
!!     The memory deallocation for an array failed. Print an error message and
!!     stop the program execution.
!!
!!   AUTHOR
!!     Matthias Krack (20.10.1999)
!!
!!   MODIFICATION HISTORY
!!     CP2K by JGH 21.08.2000
!!
!!   NOTES
!!     array  : Name of the array.
!!     routine: Name of the calling routine.
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE stop_deallocate(routine,array)

    CHARACTER(LEN=*), INTENT(IN) :: array,routine

!   *** Local variables ***

    CHARACTER(LEN=240) :: message

!   ---------------------------------------------------------------------------

    WRITE (message,"(A)")&
      "The memory deallocation for the data object <"//TRIM(array)//"> failed"

    CALL stop_program(routine,message)

  END SUBROUTINE stop_deallocate

!!*****
!******************************************************************************
!!****** termination/stop_program [1.0] *
!!
!!   NAME
!!     stop_program
!!
!!   FUNCTION
!!     Stop the program run and write an error message.
!!
!!   AUTHOR
!!     Matthias Krack (28.08.1996)
!!
!!   MODIFICATION HISTORY
!!     Translated to Fortran 90 (07.10.99, MK)
!!
!!   NOTES
!!     routine: Name of the calling routine
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE stop_program(routine,error_message)

    CHARACTER(LEN=*), INTENT(IN) :: error_message,routine

!   ---------------------------------------------------------------------------

!   *** Print the error message ***

    CALL print_message("ERROR in "//TRIM(routine),output_unit,2,2,0)
    CALL print_message(error_message,output_unit,1,1,1)

    CALL mp_stop(0,"stop_program")

  END SUBROUTINE stop_program

!!*****
!******************************************************************************

END MODULE termination

!******************************************************************************
