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
                                             cp_iteration_info_type,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             max_iter_r_level
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: compress

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_output_handling'
! cp_should_output, cp_iterate, cp_add_iter_level, cp_rm_iter_level
! cp_get_output_name, create_print_key
  PUBLIC :: cp_would_output, cp_unitnr_for_output,&
       cp_write_output, cp_finished_output

  INTEGER, PARAMETER, PUBLIC :: silent_print_level=0, low_print_level=1,&
       medium_print_level=2, high_print_level=3, debug_print_level=4
       

!!***
!****************************************************************************
CONTAINS

!!****f* cp_output_handling/cp_print_key_section_create [1.0] *
!!
!!   NAME
!!     cp_print_key_section_create
!!
!!   FUNCTION
!!     creates a print_key section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - print_key_section: the print key to create
!!     - name: the name of the print key
!!     - description: the description of the print key
!!     - print_level: print level starting at which the printing takes place
!!       (defaults to debug_print_level)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  SUBROUTINE cp_print_key_section_create(print_key_section, name, description, &
       print_level,error)
    TYPE(section_type), POINTER              :: print_key_section
    CHARACTER(len=*), INTENT(in)             :: name, description
    INTEGER, INTENT(in), OPTIONAL            :: print_level
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_section_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_print_level
    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(print_key_section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_print_level=debug_print_level
       IF (PRESENT(print_level)) my_print_level=print_level

       CALL section_create(print_key_section,name=name,description=description,&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)
       CALL keyword_create(keyword, name="PRINT",&
            description="Level starting at which this proprety is printed",&
            usage="PRINT",&
            default_i_val=my_print_level,lone_keyword_i_val=silent_print_level,&
            enum_c_vals=(/"silent","low   ","medium","high  ","debug "/),&
            enum_i_vals=(/ silent_print_level, low_print_level,&
            medium_print_level,high_print_level,debug_print_level/),&
            error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EACH",&
            description="How often this proprety is printed, this is matched with"//&
            " the actual iteration level from the right replacing non present "//&
            "levels with 1. last_symbolic uses 'l' as iteration number in the "//&
            "name. Not every process is able to identify the last iteration "//&
            "early enough to be able to output.",&
            usage="EACH 2 last_symbolic # each 2 external iterations prints the last internal iteration",&
            default_i_val=my_print_level,lone_keyword_i_val=silent_print_level,&
            enum_c_vals=(/&
            "never         ",&
            "last          ",&
            "last_symbolic "/),&
            enum_i_vals=(/ 0,-1,-2 /),&
            enum_strict=.FALSE., error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE cp_print_key_section_create
!***************************************************************************

!!****f* cp_output_handling/cp_should_output [1.0] *
!!
!!   NAME
!!     cp_should_output
!!
!!   FUNCTION
!!     returns if the given proprety should be printed out
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - iteration_info: information about the actual iteration level
!!     - print_key: the section values of the key to be printed
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_should_output(iteration_info,print_key,error) RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_should_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, ival, print_level
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
  CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
  CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
  res=.FALSE.
  IF (.NOT. failure) THEN
     NULLIFY(i_vals)
     CALL section_vals_val_get(print_key,"PRINT",i_val=print_level,error=error)
     CALL section_vals_val_get(print_key,"EACH",i_vals=i_vals,error=error)
     CPPrecondition(ASSOCIATED(i_vals),cp_failure_level,routineP,error,failure)
     IF (SIZE(i_vals)>=iteration_info%n_rlevel.AND.&
          iteration_info%print_level>=print_level) THEN
        res=.TRUE.
        DO ilevel=1,SIZE(i_vals)
           ival=i_vals(SIZE(i_vals)-ilevel+1)
           SELECT CASE(ival)
           CASE(0)
              res=.FALSE.
           CASE(-1,-2)
              IF (.NOT.iteration_info%last_iter&
                   (iteration_info%n_rlevel-ilevel+1)) THEN
                 res=.FALSE.
              END IF
           CASE default
              CPPrecondition(ival>0,cp_failure_level,routineP,error,failure)
              IF (MODULO(iteration_info%iteration &
                   (iteration_info%n_rlevel-ilevel+1)-1,ival)/=0) THEN
                 res=.FALSE.
              END IF
           END SELECT
        END DO
     END IF
  END IF
END FUNCTION cp_should_output
!***************************************************************************

!!****f* cp_output_handling/cp_iter_string [1.0] *
!!
!!   NAME
!!     cp_iter_string
!!
!!   FUNCTION
!!     returns the iteration string, a string that is useful to create
!!     unique filenames (once you trim it)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - iter_info: the iteration info from where to take the iteration
!!       number
!!     - print_key: the print key to optionally show the last iteration
!!       symbolically
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_iter_string(iter_info,print_key,error) RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iter_info
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: print_key
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=7*max_iter_r_level)        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iter_string', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, shift
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: my_print_key

  failure=.FALSE.

  res=""
  CPPrecondition(ASSOCIATED(iter_info),cp_failure_level,routineP,error,failure)
  CPPrecondition(iter_info%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(my_print_key)
     IF (PRESENT(print_key)) my_print_key => print_key
     IF (ASSOCIATED(my_print_key)) THEN
        NULLIFY(i_vals)
        CALL section_vals_val_get(my_print_key,"EACH",i_vals=i_vals,&
             error=error)
        CPPrecondition(ASSOCIATED(i_vals),cp_failure_level,routineP,error,failure)
        shift=iter_info%n_rlevel-SIZE(i_vals)
        DO ilevel=1,iter_info%n_rlevel
           IF (ilevel>shift) THEN
              IF (i_vals(ilevel-shift)==-2) THEN
                 WRITE(res(7*ilevel-6:7*ilevel),"('l_')")
              ELSE
                 WRITE(res(7*ilevel-6:7*ilevel),"(i6,'_')") iter_info%iteration(ilevel)
              END IF
           ELSE
              WRITE(res(7*ilevel-6:7*ilevel),"(i6,'_')") iter_info%iteration(ilevel)
           END IF
        END DO
     ELSE
        DO ilevel=1,iter_info%n_rlevel
           WRITE(res(7*ilevel-6:7*ilevel),"(i6,'_')") &
             iter_info%iteration(ilevel)
        END DO
     END IF
     CALL compress(res,.TRUE.)
     IF (LEN_TRIM(res)>0) THEN
        res(LEN_TRIM(res):LEN_TRIM(res))=" "
     END IF
  END IF
END FUNCTION cp_iter_string
!***************************************************************************

!!****f* cp_output_handling/cp_iterate *
!!
!!   NAME
!!     cp_iterate
!!
!!   FUNCTION
!!     adds one to the actual iteration
!!
!!   NOTES
!!     this is supposed to be called at the beginning of each iteration
!!
!!   INPUTS
!!     - iteration_info: the iteration info to update
!!     - last: if this iteration is the last one (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_iterate(iteration_info,last,iter_nr,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    LOGICAL, INTENT(in), OPTIONAL            :: last
    INTEGER, INTENT(in), OPTIONAL            :: iter_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterate', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_last

  failure=.FALSE.
  my_last=.FALSE.
  IF (PRESENT(last)) my_last=last
  
  CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
     IF (PRESENT(iter_nr)) THEN
        iteration_info%iteration(iteration_info%n_rlevel)=iter_nr
     ELSE
        iteration_info%iteration(iteration_info%n_rlevel)=&
             iteration_info%iteration(iteration_info%n_rlevel)+1
     END IF
     iteration_info%last_iter(iteration_info%n_rlevel)=my_last
  END IF
END SUBROUTINE cp_iterate
!***************************************************************************

!!****f* cp_output_handling/cp_add_iter_level *
!!
!!   NAME
!!     cp_add_iter_level
!!
!!   FUNCTION
!!     Adds an iteration level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - iteration_info: the iteration info to which an iteration level has
!!       to be added
!!     - n_level_new: number of iteration levels after this call
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_add_iter_level(iteration_info,n_level_new,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    INTEGER, INTENT(out), OPTIONAL           :: n_level_new
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_add_iter_level', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
  CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(iteration_info%n_rlevel<max_iter_r_level,cp_failure_level,routineP,error,failure)
     iteration_info%n_rlevel=iteration_info%n_rlevel+1
     iteration_info%iteration(iteration_info%n_rlevel)=0
     iteration_info%last_iter(iteration_info%n_rlevel)=.FALSE.
     IF (PRESENT(n_level_new)) n_level_new=iteration_info%n_rlevel
  END IF
END SUBROUTINE cp_add_iter_level
!***************************************************************************

!!****f* cp_output_handling/cp_remove_iter_level *
!!
!!   NAME
!!     cp_remove_iter_level
!!
!!   FUNCTION
!!     Removes an iteration level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - iteration_info: the iteration info to which an iteration level has
!!       to be removed
!!     - n_rlevel_att: iteration level before the call (to do some checks)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_remove_iter_level(iteration_info,n_rlevel_att,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    INTEGER, INTENT(in), OPTIONAL            :: n_rlevel_att
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_remove_iter_level', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
  CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(n_rlevel_att)) THEN
        CPPrecondition(n_rlevel_att==iteration_info%n_rlevel,cp_failure_level,routineP,error,failure)
     END IF
     iteration_info%iteration(iteration_info%n_rlevel)=0
     iteration_info%last_iter(iteration_info%n_rlevel)=.FALSE.
     iteration_info%n_rlevel=iteration_info%n_rlevel-1
  END IF
END SUBROUTINE cp_remove_iter_level
!***************************************************************************

!!****f* cp_output_handling/cp_would_output [1.0] *
!!
!!   NAME
!!     cp_would_output
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

