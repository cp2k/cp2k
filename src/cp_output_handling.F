!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_output_handling [1.0] *
!!
!!   NAME
!!     cp_output_handling
!!
!!   FUNCTION
!!     routines to handle the output, The idea is to remove the 
!!     decision of wheter to output and what to output from the code
!!     that does the output, and centralize it here.
!!
!!   NOTES
!!     These were originally together with the log handling routines,
!!     but have been spawned off. Some dependencies are still there,
!!     and some of the comments about log handling also applies to output
!!     handling: @see cp_log_handling
!!
!!     A section with output name, output every (-1: last, 0: never, 1: always)
!!     separate File, discriminate fromWhere,... should be implemented.
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_output_handling
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_output_handling'

  PUBLIC :: cp_would_output, cp_unitnr_for_output,&
       cp_write_output, cp_finished_output

!!***
!****************************************************************************
CONTAINS

!!****f* cp_output_handling/cp_would_output [1.0] *
!!
!!   NAME
!!     cp_would_output
!!
!!   SYNOPSIS
!!     Function cp_would_output(logger, outputname, fromwhere, iter, local,&
!!         error)
!!       Logical:: cp_would_output
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Function cp_would_output
!!
!!   FUNCTION
!!    this function can be called to check if the logger would
!!    write the output. This can be used to know if the generation 
!!    of some costly output is necessary
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger that decides where the logging should go
!!     - outputName: the name of the output, for example 'scf_hamiltonian'
!!       this should be unique (i.e. no other routine should have the same
!!       name, because it is used to build an unique identifier)
!!     - fromWhere: string of the form module:function or file:lineNr
!!       that gives the origin of the output
!!     - iter: the iteration number, or some other number that permit to 
!!       differentiate between different instances of the same object
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     08.2002 updated to new logger [fawzi]
!!
!!*** *********************************************************************
! pure 
  FUNCTION cp_would_output(logger, outputName , fromWhere, iter, local, error)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(in)             :: outputName, fromWhere
    INTEGER, INTENT(in)                      :: iter
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: cp_would_output

    cp_would_output=.TRUE.
  END FUNCTION cp_would_output
!****************************************************************************

!!****f* cp_output_handling/cp_unitnr_for_output [1.0] *
!!
!!   NAME
!!     cp_unitnr_for_output
!!
!!   SYNOPSIS
!!     Function cp_unitnr_for_output(logger, outputname, fromwhere, iter,&
!!         local, error) Result(res)
!!       Integer:: res
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (IN):: outputname, fromwhere
!!       Integer, Intent (IN):: iter
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Function cp_unitnr_for_output
!!
!!   FUNCTION
!!    returns the unit nr for the requested kind of output.
!!
!!    you should call cp_finished_output when an iteration output is
!!    finished (to immediately close the file that might have been opened)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     08.2002 updated to new logger [fawzi]
!!
!!*** *********************************************************************
!pure 
  FUNCTION cp_unitnr_for_output(logger, outputName, fromWhere, iter, local,&
       error) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(in)             :: outputName, fromWhere
    INTEGER, INTENT(in)                      :: iter
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    res=cp_logger_get_default_unit_nr(logger,local=local)
  END FUNCTION cp_unitnr_for_output
!****************************************************************************

!!****f* cp_output_handling/cp_finished_output [1.0] *
!!
!!   NAME
!!     cp_finished_output
!!
!!   SYNOPSIS
!!     Subroutine cp_finished_output(logger, outputname, fromwhere, iter,&
!!         unit_nr, local, error)
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (IN):: outputname, fromwhere
!!       Integer, Intent (IN):: iter, unit_nr
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_finished_output
!!
!!   FUNCTION
!!     should be called after you finish working with a unit obtained with
!!     cp_unitnr_for_output, so that the file that might have been opened
!!     can be closed.
!!
!!   NOTES
!!     at the moment does nothing, but this should change in the future
!!
!!   INPUTS
!!     the inputs should be exactly the same of the corresponding
!!     cp_unitnr_for_output
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_finished_output(logger, outputName, fromWhere, iter, &
       unit_nr,local,error)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(in)             :: outputName, fromWhere
    INTEGER, INTENT(in)                      :: iter, unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_finished_output', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

  END SUBROUTINE cp_finished_output
!***************************************************************************

!!****f* cp_output_handling/cp_write_output [1.0] *
!!
!!   NAME
!!     cp_write_output
!!
!!   SYNOPSIS
!!     Subroutine cp_write_output(logger, outputname, fromwhere, iter,&
!!         message, local, error)
!!       Type(cp_logger_type), Pointer:: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere, message
!!       Integer, Intent (IN):: iter
!!       Logical, Intent (IN), Optional:: local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_write_output
!!
!!   FUNCTION
!!     Used to write the result (output) of something
!!     no frills, no label, just the message is written.
!!     In the future some discrimination based on outputName could be 
!!     introduced.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger that decides where the logging should go
!!     - outputName: the name of the output, for example 'scf_hamiltonian'.
!!       this should be unique (i.e. no other routine should have the same
!!       name, because it is used to build an unique identifier)
!!     - fromWhere: string of the form module:function or file:lineNr
!!       that gives the origin of the output
!!     - iter: the iteration number, or some other number that permit to 
!!       differentiate between different instances of the same object
!!     - local: it the output is local to the task or not. defaults to false
!!     
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     08.2002 updated to new logger [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE cp_write_output(logger, outputName , fromWhere , iter, message,&
       local,error)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(in)             :: outputName, fromWhere
    INTEGER, INTENT(in)                      :: iter
    CHARACTER(len=*), INTENT(in)             :: message
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_output_handling:cp_write_output'

    INTEGER                                  :: oErr, unitNr
    LOGICAL                                  :: failure, my_local

    failure=.FALSE.
    my_local=.FALSE.
    IF (PRESENT(local)) my_local=local
    IF (cp_would_output(logger,outputName,fromWhere,iter,&
         local=my_local,error=error)) THEN
       IF (my_local.OR.logger%para_env%mepos==logger%para_env%source) THEN
          unitNr=cp_unitnr_for_output(logger,outputName,fromWhere,iter,&
               local=my_local,error=error) 
          WRITE (unit=unitNr,fmt='(a)', iostat=oErr) message
          CPPostcondition(oErr == 0,cp_failure_level,routineP,error,failure)
          CALL cp_finished_output(logger,outputName,fromWhere,iter,&
               unit_nr=unitNr,local=my_local,error=error)
       END IF
    END IF
  END SUBROUTINE cp_write_output
!****************************************************************************

END MODULE cp_output_handling

