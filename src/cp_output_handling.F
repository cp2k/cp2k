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
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!     08.2002 updated to new logger [fawzi]
!!     10.2004 big rewrite of the output methods, connected to the new
!!             input, and iteration_info [fawzi]
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_output_handling
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_iteration_info_type,&
                                             cp_log,&
                                             cp_logger_generate_filename,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             max_iter_r_level
  USE input_constants,                 ONLY: debug_print_level,&
                                             high_print_level,&
                                             low_print_level,&
                                             medium_print_level,&
                                             silent_print_level
  USE input_enumeration_types,         ONLY: enum_create,&
                                             enum_release,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE string_utilities,                ONLY: compress

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_output_handling'
  PUBLIC :: cp_should_output, cp_iterate, cp_add_iter_level, cp_rm_iter_level
  PUBLIC :: cp_iter_string, cp_print_key_section_create
  PUBLIC :: cp_unitnr_for_output, cp_finished_output, cp_first_output

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
       print_level,each,add_last,filename,common_iter_levels,error)
    TYPE(section_type), POINTER              :: print_key_section
    CHARACTER(len=*), INTENT(in)             :: name, description
    INTEGER, INTENT(in), OPTIONAL            :: print_level
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: each, add_last
    CHARACTER, INTENT(in), OPTIONAL          :: filename
    INTEGER, INTENT(in), OPTIONAL            :: common_iter_levels
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_section_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: my_filename
    INTEGER                                  :: my_comm_iter_levels, &
                                                my_print_level
    LOGICAL                                  :: failure
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(val_type), POINTER                  :: val

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(print_key_section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_print_level=debug_print_level
       IF (PRESENT(print_level)) my_print_level=print_level

       CALL section_create(print_key_section,name=name,description=description,&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, val,enum)
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

       IF (PRESENT(each)) THEN
          CALL val_create(val,i_vals=each,error=error)
       ELSE
          CALL val_create(val,i_val=1,error=error)
       END IF
       CALL keyword_create(keyword, name="EACH",&
            description="How often this proprety is printed, this is matched with"//&
            " the actual iteration level from the right replacing non present "//&
            "levels with 1. How to handle the last iteration is treated "//&
            "separately in ADD_LAST (this mean that EACH 0 might print the last "//&
            "iteration)",&
            usage="EACH 2 5",&
            default_val=val,error=error)
       CALL val_release(val,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL enum_create(enum,c_vals=(/&
            "no      ",&
            "numeric ",&
            "symbolic"/),&
            i_vals=(/ 0,1,2 /),error=error)
       IF (PRESENT(add_last)) THEN
          CALL val_create(val,i_vals=add_last, enum=enum, error=error)
       ELSE
          CALL val_create(val,i_val=0, enum=enum,error=error)
       END IF
       CALL enum_release(enum,error=error)
       CALL keyword_create(keyword, name="add_last",&
               description="If the last iteration should be added, and if it "//&
               "should be marked symbolically (with l) or with the iteration "//&
               "number."//&
               "Not every process is able to identify the last iteration "//&
               "early enough to be able to output.",&
               usage="ADD_LAST no symbolic",&
               default_val=val,&
               error=error)
       CALL val_release(val,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       my_comm_iter_levels=0
       IF (PRESENT(common_iter_levels)) my_comm_iter_levels=common_iter_levels
       CALL keyword_create(keyword, name="COMMON_ITERATION_LEVELS",&
            description="How many iterations levels should be written"//&
            " in the same file (no extra information about the actual"//&
            " iteration level is written to the file)",&
            usage="COMMON_ITERATION_LEVELS 1",&
            default_i_val=my_comm_iter_levels, error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       my_filename=" "
       IF (PRESENT(filename)) my_filename=filename
       CALL keyword_create(keyword, name="FILENAME",&
            description="Root of the filename used for the output."//&
            ' If not given uses the project name + "'//TRIM(name)//&
            '"; "__STD_OUT__" (case sensitive, and spelled exactly as here)'//&
            ' directs the output to the standard log (normally the standard '//&
            'out).',&
            usage="FILENAME somwhere/nice",&
            default_c_val=my_filename, error=error)
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
!!     - control_val: 0 output is always false, if 1 output is
!!       controlled by the print key if 2 output is always done.
!!       Defaults to 1.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  FUNCTION cp_should_output(iteration_info,print_key,control_val,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    INTEGER, INTENT(in), OPTIONAL            :: control_val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_should_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, ival, my_control_val, &
                                                print_level
    INTEGER, DIMENSION(:), POINTER           :: add_last, each
    LOGICAL                                  :: failure, level_passed

    failure=.FALSE.
    my_control_val=1
    IF (PRESENT(control_val)) my_control_val=control_val
    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (my_control_val==2) THEN
       res=.TRUE.
    ELSE IF (my_control_val==0) THEN
       res=.FALSE.
    ELSE
       IF (.NOT.ASSOCIATED(print_key)) THEN
          res=(iteration_info%print_level > debug_print_level)
       ELSE
          CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
          CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
          res=.FALSE.
          IF (.NOT. failure) THEN
             NULLIFY(each,add_last)
             CALL section_vals_val_get(print_key,"PRINT",i_val=print_level,error=error)
             CALL section_vals_val_get(print_key,"EACH",i_vals=each,error=error)
             CALL section_vals_val_get(print_key,"ADD_LAST",i_vals=add_last,error=error)
             CPPrecondition(ASSOCIATED(each),cp_failure_level,routineP,error,failure)
             CPPrecondition(ASSOCIATED(add_last),cp_failure_level,routineP,error,failure)
             IF (SIZE(each)<=iteration_info%n_rlevel.AND.&
                  iteration_info%print_level>=print_level) THEN
                res=.TRUE.
                DO ilevel=1,SIZE(each)
                   level_passed=.FALSE.
                   ival=each(SIZE(each)-ilevel+1)
                   IF (MODULO(iteration_info%iteration &
                        (iteration_info%n_rlevel-ilevel+1)-1,ival)==0) THEN
                      level_passed=.TRUE.
                   END IF
                   IF (SIZE(add_last)>=ilevel) THEN
                      ival=add_last(SIZE(add_last)-ilevel+1)
                      IF (iteration_info%last_iter(iteration_info%n_rlevel-ilevel+1)&
                           .AND.(ival==1.OR.ival==2)) THEN
                         level_passed=.TRUE.
                      END IF
                   END IF
                   IF (.not.level_passed) res=.FALSE.
                END DO
             END IF
          END IF
       END IF
    END IF
END FUNCTION cp_should_output
!***************************************************************************

!!****f* cp_output_handling/cp_first_output [1.0] *
!!
!!   NAME
!!     cp_first_output
!!
!!   FUNCTION
!!     returns true if the given proprety is printed out for the first time
!!     (so that some header or extra information can be printed out)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - iteration_info: information about the actual iteration level
!!     - print_key: the section values of the key to be printed
!!     - control_val: 0 output is always false, if 1 output is
!!       controlled by the print key if 2 output is always done.
!!       Defaults to 1.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  FUNCTION cp_first_output(iteration_info,print_key,control_val,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    INTEGER, INTENT(in), OPTIONAL            :: control_val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_first_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, ival, my_control_val, &
                                                print_level
    INTEGER, DIMENSION(:), POINTER           :: add_last, each
    LOGICAL                                  :: failure, level_passed

    failure=.FALSE.
    my_control_val=1
    IF (PRESENT(control_val)) my_control_val=control_val
    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (my_control_val==2) THEN
       res=ALL(iteration_info%iteration(:iteration_info%n_rlevel)==1)
    ELSE IF (my_control_val==0) THEN
       res=.FALSE.
    ELSE
       IF (.NOT.ASSOCIATED(print_key)) THEN
          res=(iteration_info%print_level > debug_print_level)
       ELSE
          CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
          CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
          res=.FALSE.
          IF (.NOT. failure) THEN
             NULLIFY(each,add_last)
             CALL section_vals_val_get(print_key,"PRINT",i_val=print_level,error=error)
             CALL section_vals_val_get(print_key,"EACH",i_vals=each,error=error)
             CALL section_vals_val_get(print_key,"ADD_LAST",i_vals=add_last,error=error)
             CPPrecondition(ASSOCIATED(each),cp_failure_level,routineP,error,failure)
             CPPrecondition(ASSOCIATED(add_last),cp_failure_level,routineP,error,failure)
             IF (SIZE(each)<=iteration_info%n_rlevel.AND.&
                  iteration_info%print_level>=print_level) THEN
                res=.TRUE.
                DO ilevel=1,SIZE(each)
                   level_passed=.FALSE.
                   ival=each(SIZE(each)-ilevel+1)
                   IF (iteration_info%iteration &
                        (iteration_info%n_rlevel-ilevel+1)==1) THEN
                      level_passed=.TRUE.
                   END IF
                   IF (.not.level_passed) res=.FALSE.
                END DO
             END IF
          END IF
       END IF
    END IF
END FUNCTION cp_first_output
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
!!     - for_file: if the string is to be used for file generation 
!!       (and should consequently ignore some iteration levels depending
!!       on COMMON_ITERATION_LEVELS).
!!       Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_iter_string(iter_info,print_key,for_file,error) RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iter_info
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: print_key
    LOGICAL, INTENT(in), OPTIONAL            :: for_file
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=7*max_iter_r_level)        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iter_string', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: c_i_level, ilevel, n_rlevel, &
                                                shift
    INTEGER, DIMENSION(:), POINTER           :: add_last, each
    LOGICAL                                  :: failure, my_for_file
    TYPE(section_vals_type), POINTER         :: my_print_key

  failure=.FALSE.

  res=""
  my_for_file=.FALSE.
  IF (PRESENT(for_file)) my_for_file=for_file
  CPPrecondition(ASSOCIATED(iter_info),cp_failure_level,routineP,error,failure)
  CPPrecondition(iter_info%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(my_print_key)
     IF (PRESENT(print_key)) my_print_key => print_key
     IF (ASSOCIATED(my_print_key)) THEN
        NULLIFY(each,add_last)
        CALL section_vals_val_get(my_print_key,"EACH",i_vals=each,&
             error=error)
        CALL section_vals_val_get(my_print_key,"add_last",i_vals=add_last,&
             error=error)
        CALL section_vals_val_get(my_print_key,"COMMON_ITERATION_LEVELS",&
             i_val=c_i_level, error=error)
        CPPrecondition(ASSOCIATED(each),cp_failure_level,routineP,error,failure)
        shift=iter_info%n_rlevel-SIZE(add_last)
        n_rlevel=iter_info%n_rlevel
        IF (my_for_file) n_rlevel=MIN(n_rlevel,MAX(0,n_rlevel-c_i_level))
        DO ilevel=1,n_rlevel
           IF (iter_info%last_iter(ilevel).and.ilevel>shift) THEN
              IF (add_last(ilevel-shift)==2) THEN
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
SUBROUTINE cp_add_iter_level(iteration_info,n_rlevel_new,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    INTEGER, INTENT(out), OPTIONAL           :: n_rlevel_new
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
     IF (PRESENT(n_rlevel_new)) n_rlevel_new=iteration_info%n_rlevel
  END IF
END SUBROUTINE cp_add_iter_level
!***************************************************************************

!!****f* cp_output_handling/cp_rm_iter_level *
!!
!!   NAME
!!     cp_rm_iter_level
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
SUBROUTINE cp_rm_iter_level(iteration_info,n_rlevel_att,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    INTEGER, INTENT(in), OPTIONAL            :: n_rlevel_att
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_rm_iter_level', &
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
END SUBROUTINE cp_rm_iter_level
!***************************************************************************

!!****f* cp_output_handling/cp_unitnr_for_output [1.0] *
!!
!!   NAME
!!     cp_unitnr_for_output
!!
!!   FUNCTION
!!    Utility function that opens a file with a unique filename, generated from
!!    the print_key name and iteration info, and returns the unit nr.
!!    (in future the print key might allow for more control on the filename, or
!!    use of the std_out)
!!
!!    you should call cp_finished_output when an iteration output is
!!    finished (to immediately close the file that has likely been opened)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger for the parallel environment, iteration info
!!       and filename generation
!!     - print_key: section with the information, name, label to be used,...
!!       for the output
!!     - local: if the unit should be local to this task, or global to the
!!       program (defaults to false).
!!     - log_filename: if a small log with the filename should be written
!!       to the main log (defaults to false)
!!     - extension: extension to be applied to the filename (including the ".")
!!     - middle_name: name to be added to the generated filename, useful when
!!       print_key activates different distinct outputs, to be able to 
!!       distinguish them
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!
!!*** *********************************************************************
FUNCTION cp_unitnr_for_output(logger, print_key, extension, middle_name,&
     local, log_filename, error) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key
    CHARACTER(len=*), INTENT(in)             :: extension
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: middle_name
    LOGICAL, INTENT(in), OPTIONAL            :: local, log_filename
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_unitnr_for_output', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename
    CHARACTER(len=default_string_length)     :: my_middle_name, outPath
    INTEGER                                  :: i
    LOGICAL                                  :: failure, is_path, my_local

  failure=.FALSE.
  my_local=.FALSE.
  IF (PRESENT(local)) my_local=local
  CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
     CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     IF (my_local.OR.&
          logger%para_env%mepos==logger%para_env%source) THEN
        CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)
        IF (outPath=='__STD_OUT__') THEN
           res=cp_logger_get_default_unit_nr(logger,local=my_local)
        ELSE
           is_path=.FALSE.
           DO i=1,LEN_TRIM(outPath)
              IF (outPath(i:i)=="/") is_path=.TRUE.
           END DO
           IF (.not.is_path) THEN
              if (outPath==' ') my_middle_name=print_key%section%name
              IF (PRESENT(middle_name)) my_middle_name=middle_name
              CALL cp_logger_generate_filename(logger,res=filename,&
                   root=TRIM(logger%iter_info%project_name)//"-"//&
                   TRIM(outPath)//TRIM(my_middle_name),&
                   postfix="-"//TRIM(cp_iter_string(logger%iter_info,&
                   print_key=print_key,for_file=.TRUE.,error=error))//&
                   extension,local=my_local)
           ELSE
              my_middle_name=' '
              IF (PRESENT(middle_name)) my_middle_name=middle_name
              CALL cp_logger_generate_filename(logger,res=filename,&
                   root=TRIM(outPath)//TRIM(my_middle_name),&
                   postfix="-"//TRIM(cp_iter_string(logger%iter_info,&
                   print_key=print_key,for_file=.TRUE.,error=error))//&
                   extension,local=my_local)
           END IF
           CALL open_file(file_name=filename,file_status="UNKNOWN",&
                file_form="FORMATTED",file_action="WRITE",&
                file_position="APPEND",unit_number=res)
           IF (PRESENT(log_filename)) THEN
              IF (log_filename) THEN
                 CALL cp_log(logger=logger, level=cp_note_level, &
                      fromWhere=routineP , &
                      message="writing "//TRIM(print_key%section%name)//" "//&
                      TRIM(cp_iter_string(logger%iter_info))//" to "//&
                      TRIM(filename),local=my_local)
              END IF
           END IF
        END IF
     ELSE
        res=-1
     END IF
  END IF
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
!!     at the moment always closes the unit
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
  SUBROUTINE cp_finished_output(unit_nr, logger, print_key, local,&
       error)
    INTEGER, INTENT(inout)                   :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_finished_output', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_local

    failure=.FALSE.
    my_local=.FALSE.
    IF (PRESENT(local)) my_local=local
    CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
    CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (my_local.OR.&
         logger%para_env%source==logger%para_env%mepos) THEN
       CALL close_file(unit_nr,"KEEP")
    END IF
  END SUBROUTINE cp_finished_output
!***************************************************************************

END MODULE cp_output_handling

