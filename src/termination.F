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

  USE global_types,                    ONLY: global_environment_type
  USE message_passing,                 ONLY: mp_abort
  USE output_utilities,                ONLY: print_message
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: print_stack

  IMPLICIT NONE

  PRIVATE

  INTEGER :: error_unit

  PUBLIC :: get_error_unit,&
            set_error_unit,&
            stop_memory,&
            stop_program

! *****************************************************************************

  INTERFACE stop_memory
    MODULE PROCEDURE stop_allocate,stop_allocate_new,&
                     stop_deallocate,stop_deallocate_new
  END INTERFACE

  INTERFACE stop_program
    MODULE PROCEDURE stop_program_new,stop_program_old
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

    INTEGER, INTENT(OUT)                     :: lunit

!   ---------------------------------------------------------------------------

    lunit = error_unit

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

    INTEGER, INTENT(IN)                      :: lunit

!   ---------------------------------------------------------------------------

    error_unit = lunit

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

    CHARACTER(LEN=*), INTENT(IN)             :: routine, array
    INTEGER, INTENT(IN)                      :: memory

    CHARACTER(LEN=240)                       :: message

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

! *****************************************************************************

  SUBROUTINE stop_allocate_new(routine_name,module_name,line_number,object,&
                               memory)

    CHARACTER(LEN=*), INTENT(IN)             :: routine_name, module_name
    INTEGER, INTENT(IN)                      :: line_number
    CHARACTER(LEN=*), INTENT(IN)             :: object
    INTEGER, INTENT(IN)                      :: memory

    CHARACTER(LEN=240)                       :: message

!   ---------------------------------------------------------------------------

    IF (memory == 0) THEN
      message = "The memory allocation for the data object <"//TRIM(object)//&
                "> failed"
    ELSE
      WRITE (message,"(A,I12,A)")&
        "The memory allocation for the data object <"//TRIM(object)//&
        "> failed. The requested memory size is ",memory/1024," Kbytes"
      CALL compress(message)
    END IF

    CALL stop_program(routine_name,module_name,line_number,message)

  END SUBROUTINE stop_allocate_new

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

    CHARACTER(LEN=*), INTENT(IN)             :: routine, array

    CHARACTER(LEN=240)                       :: message

!   ---------------------------------------------------------------------------

    WRITE (message,"(A)")&
      "The memory deallocation for the data object <"//TRIM(array)//"> failed"

    CALL stop_program(routine,message)

  END SUBROUTINE stop_deallocate

! *****************************************************************************

  SUBROUTINE stop_deallocate_new(routine_name,module_name,line_number,object)

    CHARACTER(LEN=*), INTENT(IN)             :: routine_name, module_name
    INTEGER, INTENT(IN)                      :: line_number
    CHARACTER(LEN=*), INTENT(IN)             :: object

    CHARACTER(LEN=240)                       :: message

!   ---------------------------------------------------------------------------

    message = "The memory deallocation for the data object <"//TRIM(object)//&
              "> failed"

    CALL stop_program(routine_name,module_name,line_number,message)

  END SUBROUTINE stop_deallocate_new

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

  SUBROUTINE stop_program_old(routine,error_message,globenv)

    CHARACTER(LEN=*), INTENT(IN)             :: routine, error_message
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    LOGICAL                                  :: ionode

!   ---------------------------------------------------------------------------

    IF (PRESENT(globenv)) THEN
      ionode = globenv%ionode
    ELSE
      ionode = .TRUE.
    END IF

!   *** Print the error message ***

    IF (ionode) THEN
      CALL print_message("ERROR in "//TRIM(routine),error_unit,2,2,0)
      CALL print_message(error_message,error_unit,1,1,1)
      CALL print_stack(error_unit)
      CALL mp_abort()
    END IF

  END SUBROUTINE stop_program_old

! *****************************************************************************

  SUBROUTINE stop_program_new(routine_name,module_name,line_number,&
                              error_message,globenv)

    CHARACTER(LEN=*), INTENT(IN)             :: routine_name, module_name
    INTEGER, INTENT(IN)                      :: line_number
    CHARACTER(LEN=*), INTENT(IN)             :: error_message
    TYPE(global_environment_type), &
      INTENT(IN), OPTIONAL                   :: globenv

    CHARACTER(LEN=80)                        :: message
    LOGICAL                                  :: ionode

!   ---------------------------------------------------------------------------

    IF (PRESENT(globenv)) THEN
      ionode = globenv%ionode
    ELSE
      ionode = .TRUE.
    END IF

!   *** Print the error message ***

    IF (ionode) THEN
      message = "ERROR in "//TRIM(routine_name)//&
                " (MODULE "//TRIM(module_name)//")"
      CALL print_message(message,error_unit,2,2,0)
      CALL print_message(error_message,error_unit,1,1,0)
      WRITE (UNIT=message,FMT="(A,I6,A)")&
        "Program stopped at line number ",line_number,&
        " of MODULE "//TRIM(module_name)
      CALL compress(message)
      CALL print_message(message,error_unit,1,1,0)
      CALL print_stack(error_unit)
      CALL mp_abort()
    END IF

  END SUBROUTINE stop_program_new

!!*****
!******************************************************************************

END MODULE termination

!******************************************************************************
