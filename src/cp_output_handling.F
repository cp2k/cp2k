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
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2001 created [fawzi]
!!     08.2002 updated to new logger [fawzi]
!!     10.2004 big rewrite of the output methods, connected to the new
!!             input, and iteration_info [fawzi]
!!     08.2005 property flags [fawzi]
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_output_handling
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_iteration_info_type,&
                                             cp_log,&
                                             cp_logger_generate_filename,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_logger_release,cp_logger_retain,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             max_iter_r_level,&
                                             cp_iter_undefined,&
                                             cp_iteration_info_retain,&
                                             cp_iteration_info_release
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
                                             section_vals_val_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_retain,&
                                             section_vals_release
  USE input_val_types,                 ONLY: val_create,&
                                             val_release,&
                                             val_type,&
                                             val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE string_utilities,                ONLY: compress

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_output_handling'
  PUBLIC :: cp_print_key_should_output, cp_iterate, cp_add_iter_level, cp_rm_iter_level
  PUBLIC :: cp_iter_string, cp_print_key_section_create
  PUBLIC :: cp_print_key_unit_nr, cp_print_key_finished_output
  public :: cp_out_flags_type
  PUBLIC :: cp_out_flags_create, cp_out_flags_retain, cp_out_flags_release,&
       cp_out_flags_set_flag, cp_out_flags_get_flag, cp_out_flags_get,&
       cp_out_flags_should_output, cp_out_flags_copy
  PUBLIC :: cp_out_flags_unit_nr, cp_out_flags_finished_output

!! flags controlling the printing and storing of a property.
!!
!! cp_out_none: do not calculate the property
!! cp_out_file_if  : if the printkey says it calculate and output the property
!! cp_out_store_if : if the printkey says it calculate and store in memory
!!                   the property
!! cp_out_file_each: calculate and output the property with the same periodicity
!!                   as said in the printkey (irrespective of the activation of 
!!                   the printkey)
!! cp_out_store_each: calculate and store the property with the same periodicity
!!                   as said in the printkey (irrespective of the activation of 
!!                   the printkey)
!! cp_out_file: always calculate and output the property
!! cp_out_store: always calculate and store in memory the property
!! cp_out_calc: just calculate the value (indipendently from the fact that there
!!              should be output)
!! cp_out_default: the default value for proprety flags (cp_out_file_if)
!!
!! this flags can be ior-ed together:
!! ior(cp_out_file_if,cp_out_store_if): if the printkey says it both print 
!!                                          and store the property
!!
!! there is no guarantee that a proprety is not stored if it is not necessary
!! not all printkeys have a control flag
INTEGER, PUBLIC, PARAMETER :: cp_p_file_if=3,cp_p_store_if=4,&
     cp_p_store=2,cp_p_file=1,cp_p_file_each=5,cp_p_store_each=6,cp_p_calc=7
INTEGER, PUBLIC, PARAMETER :: cp_out_none=0, cp_out_file_if=IBSET(0,cp_p_file_if),&
     cp_out_store_if=IBSET(0,cp_p_store_if), cp_out_file=IBSET(0,cp_p_file),&
     cp_out_store=IBSET(0,cp_p_store), cp_out_calc=IBSET(0,cp_p_calc),&
     cp_out_file_each=IBSET(0,cp_p_file_each),&
     cp_out_store_each=IBSET(0,cp_p_store_each),&
     cp_out_default=cp_out_file_if

INTEGER, SAVE, private :: last_flags_id=0

!FMINTERFACE cp_should_output
!FM   MODULE procedure cp_out_flags_should_output, cp_print_key_should_output
!FMEND INTERFACE
!FMINTERFACE cp_unitnr_for_output
!FM   MODULE PROCEDURE cp_print_key_unitnr, cp_out_flags_unit_nr
!FMEND INTERFACE
!!***
!****************************************************************************

  !!****s* cp_output_handling/cp_out_flags_type *
  !!
  !!   NAME
  !!     cp_out_flags_type
  !!
  !!   FUNCTION
  !!     stores the flags_env controlling the output of properties
  !!
  !!   NOTES
!!     Two features of this object should be:
!!       1) easy state storage, one should be able to store the state of the
!!          flags, to some changes to them just for one (or few) force evaluations
!!          and then reset the original state. The actual implementation is good
!!          in this respect
!!       2) work well with subsections. This is a problem at the moment, as
!!          if you pass just a subsection of the input the control flags get lost.
!!       A better implementation should be done storing the flags also in the
!!       input itself to be transparent
  !!
  !!   ATTRIBUTES
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - id_nr: identification number (unique to each istance)
!!     - n_flags: number of flags stored in this type
!!     - names: names of the stored flags
!!     - control_val: value of the flag
!!     - input: the input (with all the printkeys)
!!     - logger: logger and iteration information (to know if output is needed)
!!     - strict: if flags that were not stored can be read
!!     - default_val: default value of the flags that are not explicitly
!!       stored
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_out_flags_type
     INTEGER :: ref_count, id_nr, n_flags
     CHARACTER(default_string_length), DIMENSION(:), POINTER :: names
     INTEGER, DIMENSION(:), POINTER :: control_val
     TYPE(section_vals_type), pointer :: input
     TYPE(cp_logger_type), pointer :: logger
     logical :: strict
     integer :: default_val
  END TYPE cp_out_flags_type
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
    CHARACTER(LEN=*), INTENT(in), OPTIONAL          :: filename
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
       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Level starting at which this proprety is printed",&
            usage="silent",&
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
            description=' controls part of the filename for output. '//&
                        ' use __STD_OUT__ (exactly as written here) for the screen or standard logger. '//&
                        ' use filename to obtain projectname-filename. '//&
                        ' use ./filename to get filename ',&
            usage="FILENAME ./filname ",&
            default_c_val=my_filename, error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE cp_print_key_section_create
!***************************************************************************

!!****f* cp_output_handling/cp_print_key_should_output [1.0] *
!!
!!   NAME
!!     cp_print_key_should_output
!!
!!   FUNCTION
!!     returns what should be done with the given proprety
!!     if btest(res,cp_p_store) then the property should be stored in memory
!!     if btest(res,cp_p_file) then the property should be print ed to a file
!!     if res==0 then nothing should be done
!!
!!   NOTES
!!     not all the propreties support can be stored
!!
!!   ARGUMENTS
!!     - iteration_info: information about the actual iteration level
!!     - print_key: the section values of the key to be printed
!!     - control_val: flags that controls printing and storing of the
!!       printkey (see the constants defined in this module).
!!       Defaults to cp_out_default.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  FUNCTION cp_print_key_should_output(iteration_info,print_key,control_val,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    INTEGER, INTENT(in), OPTIONAL            :: control_val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_should_output', &
      routineP = moduleN//':'//routineN

    LOGICAL :: is_on, is_iter
    INTEGER                                  :: my_control_val
    
    my_control_val=cp_out_default
    is_on=cp_printkey_is_on(iteration_info,print_key,error=error)
    is_iter=cp_printkey_is_iter(iteration_info,print_key,error=error)
    IF (PRESENT(control_val)) my_control_val=control_val
    res=0

    IF (BTEST(my_control_val,cp_p_store)) THEN
       res=IBSET(res,cp_p_store)
    ELSE IF (BTEST(my_control_val,cp_p_store_if).and.is_iter.and.is_on) then
       res=IBSET(res,cp_p_store)
    ELSE IF (BTEST(my_control_val,cp_p_store_each).and.is_iter) then
       res=IBSET(res,cp_p_store)
    END IF

    IF (BTEST(my_control_val,cp_p_file)) THEN
       res=IBSET(res,cp_p_file)
    ELSE IF (BTEST(my_control_val,cp_p_file_if).and.is_iter.and.is_on) then
       res=IBSET(res,cp_p_file)
    ELSE IF (BTEST(my_control_val,cp_p_file_each).and.is_iter) then
       res=IBSET(res,cp_p_file)
    END IF
    IF (BTEST(my_control_val,cp_p_calc).OR.res/=0) THEN
       res=IBSET(res,cp_p_calc)
    END IF
  END FUNCTION cp_print_key_should_output
!***************************************************************************

!!****f* cp_output_handling/cp_printkey_is_on *
!!
!!   NAME
!!     cp_printkey_is_on
!!
!!   FUNCTION
!!     returns true if the printlevel activates this printkey
!!     does not look if this iteration it should be printed
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
FUNCTION cp_printkey_is_on(iteration_info,print_key,error) RESULT(res)
  TYPE(cp_iteration_info_type), POINTER    :: iteration_info
  TYPE(section_vals_type), POINTER         :: print_key
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  logical :: res

  logical :: failure
  character(len=*), parameter :: routineN='cp_printkey_is_on',&
        routineP=moduleN//':'//routineN
  integer :: print_level

  failure=.false.
  CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
  CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
    IF (.NOT.ASSOCIATED(print_key)) THEN
       res=(iteration_info%print_level > debug_print_level)
    ELSE
       CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
       CALL section_vals_val_get(print_key,"_SECTION_PARAMETERS_",i_val=print_level,error=error)
       res=iteration_info%print_level>=print_level
    END IF
 END IF
END FUNCTION cp_printkey_is_on
!***************************************************************************

!!****f* cp_output_handling/cp_printkey_is_iter [1.0] *
!!
!!   NAME
!!     cp_printkey_is_iter
!!
!!   FUNCTION
!!     returns if the actual iteration matches those selected by the
!!     given printkey. Does not check it the prinkey is active (at the 
!!     actual print_level)
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
  FUNCTION cp_printkey_is_iter(iteration_info,print_key,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_printkey_is_on', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ilevel, ival
    INTEGER, DIMENSION(:), POINTER           :: add_last, each
    LOGICAL                                  :: failure, level_passed

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.ASSOCIATED(print_key)) THEN
       res=(iteration_info%print_level > debug_print_level)
    ELSE
       CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
       res=.FALSE.
       IF (.NOT. failure) THEN
          NULLIFY(each,add_last)
          CALL section_vals_val_get(print_key,"EACH",i_vals=each,error=error)
          CALL section_vals_val_get(print_key,"ADD_LAST",i_vals=add_last,error=error)
          CPPrecondition(ASSOCIATED(each),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(add_last),cp_failure_level,routineP,error,failure)
          IF (SIZE(each)<=iteration_info%n_rlevel) THEN
             res=.TRUE.
             DO ilevel=1,SIZE(each)
                level_passed=.FALSE.
                ival=each(SIZE(each)-ilevel+1)
                ! output for 0,ival,2*ival,3*ival ...
                IF (ival>0) THEN
                   IF (MODULO(iteration_info%iteration &
                        (iteration_info%n_rlevel-ilevel+1),ival)==0) THEN
                      level_passed=.TRUE.
                   END IF
                END IF
                IF (SIZE(add_last)>=ilevel) THEN
                   ival=add_last(SIZE(add_last)-ilevel+1)
                   IF (iteration_info%last_iter(iteration_info%n_rlevel-ilevel+1)) THEN
                      level_passed=.TRUE.
                   END IF
                END IF
                IF (.NOT.level_passed) res=.FALSE.
             END DO
          END IF
       END IF
    END IF
  END FUNCTION cp_printkey_is_iter
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
!!     - level_name : the name of this level, for pretty printing only, right now
!!     - n_level_new: number of iteration levels after this call
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_add_iter_level(iteration_info,level_name,n_rlevel_new,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    CHARACTER(LEN=*), INTENT(IN)             :: level_name
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
     iteration_info%level_name(iteration_info%n_rlevel)=level_name
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
     iteration_info%iteration(iteration_info%n_rlevel)=cp_iter_undefined
     iteration_info%last_iter(iteration_info%n_rlevel)=.FALSE.
     iteration_info%level_name(iteration_info%n_rlevel)="DELETED"
     iteration_info%n_rlevel=iteration_info%n_rlevel-1
  END IF
END SUBROUTINE cp_rm_iter_level
!***************************************************************************

!!****f* cp_output_handling/cp_print_key_unit_nr [1.0] *
!!
!!   NAME
!!     cp_print_key_unit_nr
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
FUNCTION cp_print_key_unit_nr(logger, print_key, extension, middle_name,&
     local, log_filename, error) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key
    CHARACTER(len=*), INTENT(in)             :: extension
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: middle_name
    LOGICAL, INTENT(in), OPTIONAL            :: local, log_filename
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_unit_nr', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename,root,postfix
    CHARACTER(len=default_string_length)     :: my_middle_name, outPath
    LOGICAL                                  :: failure, my_local

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
           !
           ! complex logic to build filename:
           !   1)  Try to avoid '--' and '-.'
           !   2)  If outPath contains '/' (as in ./filename) do not prepend the project_name
           !
           ! if it is actually a full path, use it as the root
           IF (INDEX(outPath,"/").EQ.0) THEN
              IF (PRESENT(middle_name)) THEN
                   my_middle_name="-"//middle_name
              ELSE
                   IF (outPath.EQ."") THEN
                      my_middle_name="" !"-"//print_key%section%name
                   ELSE
                      my_middle_name="-"//outPath
                   ENDIF
              ENDIF
              IF (my_middle_name.EQ."-") my_middle_name=""
              root=TRIM(logger%iter_info%project_name)//TRIM(my_middle_name)
           ELSE
              root=outPath
           ENDIF
           ! use the cp_iter_string as a postfix
           postfix="-"//TRIM(cp_iter_string(logger%iter_info,&
                              print_key=print_key,for_file=.TRUE.,error=error))
           IF (TRIM(postfix)=="-") postfix=""
           ! and add the extension
           postfix=TRIM(postfix)//extension
           ! and let the logger generate the filename
           CALL cp_logger_generate_filename(logger,res=filename,&
                   root=root, postfix=postfix,local=my_local)
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
END FUNCTION cp_print_key_unit_nr
!****************************************************************************

!!****f* cp_output_handling/cp_print_key_finished_output [1.0] *
!!
!!   NAME
!!     cp_print_key_finished_output
!!
!!   FUNCTION
!!     should be called after you finish working with a unit obtained with
!!     cp_print_key_unit_nr, so that the file that might have been opened
!!     can be closed.
!!
!!   NOTES
!!     closes if the corresponding filename of the printkey is
!!      not __STD_OUT__
!!
!!   INPUTS
!!     the inputs should be exactly the same of the corresponding
!!     cp_print_key_unit_nr
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_print_key_finished_output(unit_nr, logger, print_key, local,&
       error)
    INTEGER, INTENT(inout)                   :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_finished_output', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_local
    CHARACTER(len=default_string_length)     :: outPath

    failure=.FALSE.
    my_local=.FALSE.
    IF (PRESENT(local)) my_local=local
    CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
    CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (my_local.OR.&
         logger%para_env%source==logger%para_env%mepos) THEN
        CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)
        IF (outPath.NE.'__STD_OUT__') THEN
           CALL close_file(unit_nr,"KEEP")
        ENDIF
    END IF
    unit_nr=-1
  END SUBROUTINE cp_print_key_finished_output
!***************************************************************************

!!****f* cp_output_handling/cp_iter_compare *
!!
!!   NAME
!!     cp_iter_compare
!!
!!   FUNCTION
!!     compares two iteration levels
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - iter1: the first iteration level
!!     - iter2: the second iteration level
!!     - equal: returns true is the two level are equal
!!     - common_equal: returns true if the common iteration levels are equal
!!     - ascending: returns true if iter1 < iter2
!!     - descending: retuns true if iter1 > iter2
!!     - offset: returns the difference in the last iteration level if the
!!       the others are equal; +/- HUGE if they are different
!!     - common_offset: like offset, but keeps into account only the common
!!       iteration levels
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_iter_compare(iter1,iter2,equal,common_equal, ascending,&
     descending, offset, common_offset, error)
  INTEGER, DIMENSION(max_iter_r_level), intent(in) :: iter1, iter2
  LOGICAL, OPTIONAL, INTENT(out) :: equal, common_equal, ascending, &
       descending
  INTEGER, OPTIONAL, intent(out) :: offset, common_offset
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, my_equal,my_ascending, my_common_equal,my_descending
  character(len=*), parameter :: routineN='cp_iter_compare',&
        routineP=moduleN//':'//routineN
  INTEGER :: i,my_common_offset,my_offset

  failure=.false.
  my_equal=.true.
  my_common_equal=.TRUE.
  my_ascending=.false.
  my_descending=.false.
  DO i=1,max_iter_r_level
     IF (iter1(i)/=iter2(i)) THEN
        my_equal=.false.
        IF (iter1(i)==cp_iter_undefined.OR.iter2(i)==cp_iter_undefined) THEN
           my_common_equal=.true.
           IF (iter1(i)==cp_iter_undefined) THEN
              my_offset=HUGE(0)
           ELSE
              my_offset=-HUGE(0)
           END IF
        ELSE
           my_common_equal=.false.
           IF (i<max_iter_r_level) THEN
              IF (iter1(i+1)==cp_iter_undefined.OR.iter2(i+1)==cp_iter_undefined) THEN
                 my_common_offset=iter2(i)-iter1(i)
                 IF (iter1(i+1)==cp_iter_undefined.AND.iter2(i+1)==cp_iter_undefined) THEN
                    my_offset=iter2(i)-iter1(i)
                 ELSE IF (iter1(i+1)==cp_iter_undefined) then
                    my_ascending=.true.
                    my_descending=.false.
                    my_offset=HUGE(0)
                 ELSE
                    my_ascending=.false.
                    my_descending=.true.
                    my_offset=-huge(0)
                 END IF
              ELSE
                 IF (iter1(i)<iter2(i)) THEN
                    my_ascending=.true.
                    my_descending=.false.
                    my_offset=HUGE(0)
                    my_common_offset=HUGE(0)
                 ELSE
                    my_ascending=.false.
                    my_descending=.true.
                    my_common_offset=-HUGE(0)
                    my_offset=-HUGE(0)
                 END IF
              END IF
           ELSE
              my_offset=iter2(i)-iter1(i)
              my_common_offset=my_offset
              IF (iter1(i)<iter2(i)) THEN
                 my_ascending=.TRUE.
                 my_descending=.FALSE.
              ELSE
                 my_ascending=.FALSE.
                 my_descending=.TRUE.
              END IF
           END IF
        END IF
        EXIT
     END IF
  END DO
  IF (PRESENT(equal)) equal=my_equal
  IF (PRESENT(common_equal)) common_equal=my_common_equal
  IF (PRESENT(offset)) offset=my_offset
  IF (PRESENT(common_offset)) common_offset=my_common_offset
  IF (PRESENT(ascending)) ascending=my_ascending
  IF (PRESENT(descending)) descending=my_descending
END SUBROUTINE cp_iter_compare
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_create *
!!
!!   NAME
!!     cp_out_flags_create
!!
!!   FUNCTION
!!     allocates and initializes a store for controlling printing flags
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_create(flags_env,input,logger,capacity,strict,&
     default_val,error)
  TYPE(cp_out_flags_type), pointer :: flags_env
  TYPE(section_vals_type), pointer :: input
  TYPE(cp_logger_type), pointer :: logger
  INTEGER, OPTIONAL, intent(in) :: capacity
  LOGICAL, INTENT(in), OPTIONAL :: strict
  INTEGER, INTENT(in), optional :: default_val
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_out_flags_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: my_capacity, stat

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(flags_env,stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     flags_env%ref_count=1
     last_flags_id=last_flags_id+1
     flags_env%id_nr=last_flags_id
     flags_env%n_flags=0
     my_capacity=10
     if (present(capacity)) my_capacity=capacity
     ALLOCATE(flags_env%names(my_capacity),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     ALLOCATE(flags_env%control_val(my_capacity),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     NULLIFY(flags_env%input,flags_env%logger)
     CALL section_vals_retain(input,error=error)
     flags_env%input => input
     CALL cp_logger_retain(logger)
     flags_env%logger => logger
     flags_env%strict=.false.
     if (present(strict)) flags_env%strict=strict
     flags_env%default_val=cp_out_default
     if (present(default_val)) flags_env%default_val=default_val
  END IF
  IF (.NOT.failure) THEN
     flags_env%names=""
     flags_env%control_val=0
  END IF
END SUBROUTINE cp_out_flags_create
!***************************************************************************

!!****f* cp_output_handling/cp_str_binsearch *
!!
!!   NAME
!!     cp_str_binsearch
!!
!!   FUNCTION
!!     makes a binary search in an array of strings and returns the index
!!     where the element is (or the element after which it should be inserted)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - el: the element to search
!!     - list: the list to search (should be ordered in growing order)
!!     - lb: lower bound for the search (defaults to 1)
!!     - ub: upper bound for the search (defaults to size(list))
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_str_binsearch(el,list,lb,ub,error)RESULT(res)
  CHARACTER(len=*), intent(in) :: el
  CHARACTER(len=default_string_length), DIMENSION(:), intent(in) :: list
  INTEGER, INTENT(in), OPTIONAL :: lb,ub
  type(cp_error_type), optional, intent(inout) :: error
  integer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_str_binsearch',&
        routineP=moduleN//':'//routineN
  INTEGER :: ui,li,mi

  failure=.false.
  li=1
  if (present(lb)) li=lb
  ui=size(list)
  if (present(ub)) ui=ub
  CPPrecondition(ui<=SIZE(list),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO
        if (li>ui) exit
        mi=(li+ui)/2
        IF (el<list(mi)) THEN
           ui=mi-1
        ELSE
           li=mi+1
        END IF
     END DO
     res=ui
  END IF
END FUNCTION cp_str_binsearch
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_set_flag *
!!
!!   NAME
!!     cp_out_flags_set_flag
!!
!!   FUNCTION
!!     sets a flag (adds if it is new)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - flags_env: the flags_env structure
!!     - name: the name of the flag (can be a path)
!!     - value: the value to store
!!     - virtual: if the value doe not directly correspond to a print_key,
!!       defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_set_flag(flags_env,name,value,virtual,error)
  TYPE(cp_out_flags_type), pointer :: flags_env
  CHARACTER(len=*), intent(in) :: name
  LOGICAL, INTENT(in), optional :: virtual
  integer :: value
  type(cp_error_type), optional, intent(inout) :: error
  
  character(len=*), parameter :: routineN='cp_out_flags_set_flag',&
        routineP=moduleN//':'//routineN

  LOGICAL :: failure, my_virtual
  INTEGER :: pos_name, path_end, stat, i
  INTEGER, DIMENSION(:), POINTER :: n_val
  CHARACTER(len=default_string_length), DIMENSION(:), pointer :: n_names
  TYPE(section_vals_type), pointer :: printkey
  TYPE(val_type), pointer :: val

  failure=.false.
  
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(LEN_TRIM(name)<=default_string_length,cp_failure_level,routineP,error,failure)
  my_virtual=.false.
  if (present(virtual)) my_virtual=virtual
  IF (.NOT. failure) THEN
     IF (name=="_DEFAULT_") THEN
        flags_env%default_val=value
        RETURN
     END IF
     pos_name=cp_str_binsearch(name,flags_env%names,ub=flags_env%n_flags,error=error)
     IF (pos_name>0) THEN
        IF (flags_env%names(pos_name)==name) THEN
           flags_env%control_val(pos_name)=value
           return
        END IF
     END IF
     IF (SIZE(flags_env%control_val)==flags_env%n_flags) THEN
        ALLOCATE(n_val(flags_env%n_flags+10),n_names(flags_env%n_flags+10),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        n_val(:pos_name)=flags_env%control_val(:pos_name)
        n_val(pos_name+1)=value
        n_val(pos_name+2:flags_env%n_flags+1)=flags_env%control_val(pos_name+1:flags_env%n_flags)
        n_val(flags_env%n_flags+2:)=0
        n_names(:pos_name)=flags_env%names(:pos_name)
        n_names(pos_name+1)=name
        n_names(pos_name+2:flags_env%n_flags+1)=flags_env%names(pos_name+1:flags_env%n_flags)
        n_names(flags_env%n_flags+2:)=""
        DEALLOCATE(flags_env%control_val,flags_env%names,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        flags_env%control_val => n_val
        flags_env%names => n_names
        flags_env%n_flags=flags_env%n_flags+1
     ELSE
        DO i=flags_env%n_flags,pos_name+1,-1
           flags_env%control_val(i+1)=flags_env%control_val(i)
        END DO
        DO i=flags_env%n_flags,pos_name+1,-1
           flags_env%names(i+1)=flags_env%names(i)
        END DO
        flags_env%names(pos_name+1)=name
        flags_env%control_val(pos_name+1)=value
        flags_env%n_flags=flags_env%n_flags+1
     END IF
     IF (.not.my_virtual) THEN
        path_end=INDEX(name,"/",.TRUE.)
        IF (path_end<1.OR.path_end>LEN_TRIM(name)) path_end=LEN_TRIM(name)
        printkey => section_vals_get_subs_vals(flags_env%input,name(1:path_end-1),&
             can_return_null=.TRUE.,error=error)
        CPAssert(ASSOCIATED(printkey),cp_failure_level,routineP,error,failure)
        IF (path_end<LEN_TRIM(name)) THEN
           CALL section_vals_val_get(printkey,name(path_end+1:LEN_TRIM(name)),&
                val=val,error=error)
           CPAssert(ASSOCIATED(val),cp_failure_level,routineP,error,failure)
        END IF
     END IF
  END IF
END SUBROUTINE cp_out_flags_set_flag
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_get_flag *
!!
!!   NAME
!!     cp_out_flags_get_flag
!!
!!   FUNCTION
!!     returns the value of the given flag
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - flags_env: the plage from where to take the flag
!!     - name: the name of the requested flag
!!     - was_default: returns true if the flags was not explictly present
!!       and the default value was returned.
!!       If given does not give an error with strict checking. 
!!       
!!                   
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_out_flags_get_flag(flags_env,name,was_default,error) RESULT(res)
  TYPE(cp_out_flags_type), pointer :: flags_env
  CHARACTER(len=*), intent(in) :: name
  LOGICAL, INTENT(out), optional :: was_default
  type(cp_error_type), optional, intent(inout) :: error
  integer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_out_flags_get_flag',&
        routineP=moduleN//':'//routineN
  INTEGER :: pos_name

  failure=.false.
  
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res=flags_env%default_val
     pos_name=cp_str_binsearch(name,flags_env%names,ub=flags_env%n_flags,error=error)
     IF (pos_name>0) THEN
        IF (flags_env%names(pos_name)==name) THEN
           if (present(was_default)) was_default=.false.
           res=flags_env%control_val(pos_name)
           return
        END IF
     END IF
     CALL cp_assert(.NOT.flags_env%strict.or.present(was_default),&
          cp_failure_level,cp_assertion_failed,&
          routineP,"no control flag for "//TRIM(name),error,failure)
     if (present(was_default)) was_default=.true.
  END IF
END FUNCTION cp_out_flags_get_flag
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_should_output *
!!
!!   NAME
!!     cp_out_flags_should_output
!!
!!   FUNCTION
!!     returns what operations should be done for the given printkey
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - flags_env: the flags_env to check
!!     - name: name of the printkey (can be a path)
!!     - printkey: optional, a place to *return* the printkey
!!     - virtual: if the keyword can be unassociated with a prinkey
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_out_flags_should_output(flags_env,name,print_key,virtual,error)&
     RESULT(res)
  TYPE(cp_out_flags_type), pointer :: flags_env
  CHARACTER(len=*), intent(in) :: name
  LOGICAL, OPTIONAL, intent(in) :: virtual
  TYPE(section_vals_type), POINTER, OPTIONAL :: print_key
  type(cp_error_type), optional, intent(inout) :: error
  integer :: res
  
  LOGICAL :: failure, my_virtual,sub_stat
  character(len=*), parameter :: routineN='cp_out_flags_should_output',&
        routineP=moduleN//':'//routineN
  TYPE(section_vals_type), pointer :: printkey
  TYPE(val_type), pointer :: val
  INTEGER :: control_val,path_end

  failure=.false.
  
  my_virtual=.false.
  if (present(virtual)) my_virtual=virtual
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     control_val=cp_out_flags_get_flag(flags_env,name,error=error)
     path_end=INDEX(name,"/",.TRUE.)-1
     IF (path_end<1.OR.path_end>LEN_TRIM(name)) path_end=LEN_TRIM(name)
     printkey => section_vals_get_subs_vals(flags_env%input,name(1:path_end),&
          can_return_null=.TRUE.,error=error)
     IF (ASSOCIATED(printkey).and.path_end<LEN_TRIM(name)) THEN
        CALL section_vals_val_get(printkey,name(path_end+2:LEN_TRIM(name)),&
             val=val,error=error)
        IF (ASSOCIATED(val)) THEN
           CALL val_get(val,l_val=sub_stat,error=error)
           IF (.NOT.sub_stat) THEN
              control_val=IBCLR(control_val,cp_p_store_if)
              control_val=IBCLR(control_val,cp_p_file_if)
           END IF
        ELSE
           CALL cp_assert(my_virtual,cp_failure_level,cp_assertion_failed,&
                routineP,"non virtual print_key "//TRIM(name)//&
                " has no value",error,failure)
        END IF
     ELSE IF (.NOT.ASSOCIATED(printkey)) THEN
        CALL cp_assert(my_virtual,cp_failure_level,cp_assertion_failed,&
             routineP,"non virtual print_key "//TRIM(name)//&
             " does not exist",error,failure)
     END IF
     res=cp_print_key_should_output(flags_env%logger%iter_info,printkey,control_val,&
          error=error)
     IF (PRESENT(print_key)) print_key => printkey
  END IF
END FUNCTION cp_out_flags_should_output
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_retain *
!!
!!   NAME
!!     cp_out_flags_retain
!!
!!   FUNCTION
!!     retains the given flags_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - flags_env: the flags_env to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_retain(flags_env,error)
  TYPE(cp_out_flags_type), pointer :: flags_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_out_flags_type',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
     flags_env%ref_count=flags_env%ref_count+1
  END IF
END SUBROUTINE cp_out_flags_retain
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_release *
!!
!!   NAME
!!     cp_out_flags_release
!!
!!   FUNCTION
!!     releases the given flags_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - flags_env: the flags_env to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_release(flags_env,error)
  TYPE(cp_out_flags_type), pointer :: flags_env
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_out_flags_release',&
       routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.false.
  
  IF (ASSOCIATED(flags_env)) THEN
     CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
     flags_env%ref_count=flags_env%ref_count-1
     IF (flags_env%ref_count==0) THEN
        DEALLOCATE(flags_env%names,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(flags_env%control_val,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CALL section_vals_release(flags_env%input,error=error)
        CALL cp_logger_release(flags_env%logger)
        DEALLOCATE(flags_env,stat=stat)
     END IF
  END IF
  NULLIFY(flags_env)
END SUBROUTINE cp_out_flags_release
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_set *
!!
!!   NAME
!!     cp_out_flags_set
!!
!!   FUNCTION
!!     changes various attributes in the given flags_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - flags_env: the flags_env to change
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     see cp_out_flags_type for a description of the arguments
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_set(flags_env,input,logger,default_value,&
     strict,error)
  TYPE(cp_out_flags_type), pointer :: flags_env
  TYPE(section_vals_type), pointer, optional :: input
  TYPE(cp_logger_type), POINTER, optional :: logger
  INTEGER, INTENT(in), optional :: default_value
  LOGICAL, INTENT(in), optional :: strict
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_out_flags_set',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(input)) THEN
        CALL section_vals_retain(input,error=error)
        CALL section_vals_release(flags_env%input,error=error)
        flags_env%input => input
     END IF
     IF (PRESENT(logger)) THEN
        CALL cp_logger_retain(logger)
        call cp_logger_release(flags_env%logger)
        flags_env%logger => logger
     END IF
     IF (PRESENT(default_value)) flags_env%default_val=default_value
     IF (PRESENT(strict)) flags_env%strict=strict
  END IF
END SUBROUTINE cp_out_flags_set
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_get *
!!
!!   NAME
!!     cp_out_flags_get
!!
!!   FUNCTION
!!     returns the value of various attributes of flags_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - flags_env: the flags_env to change
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     see cp_out_flags_type for a description of the arguments
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_get(flags_env,input,logger,default_val,&
     strict,error)
  TYPE(cp_out_flags_type), pointer :: flags_env
  TYPE(section_vals_type), pointer, optional :: input
  TYPE(cp_logger_type), POINTER, optional :: logger
  INTEGER, INTENT(out), optional :: default_val
  LOGICAL, INTENT(out), optional :: strict
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_out_flags_get',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(input)) input => flags_env%input
     IF (PRESENT(logger)) logger => flags_env%logger
     IF (PRESENT(default_val)) default_val=flags_env%default_val
     IF (PRESENT(strict)) strict=flags_env%strict
  END IF
END SUBROUTINE cp_out_flags_get
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_copy *
!!
!!   NAME
!!     cp_out_flags_copy
!!
!!   FUNCTION
!!     copies the contents from a flags_env to another
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - flags_env_in: the flags_env to copy
!!     - flags_env_out: the flags_env where to store the copy
!!     - merge: if the flags explicitly present in flags_env_out but not in
!!       flags_env_inshould be preserved
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE cp_out_flags_copy(flags_env_in,flags_env_out,merge,flags_only,error)
  TYPE(cp_out_flags_type), POINTER :: flags_env_in, flags_env_out
  LOGICAL, INTENT(in), OPTIONAL :: merge,flags_only
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure, my_merge, my_flags_only
  character(len=*), parameter :: routineN='cp_out_flags_copy',&
        routineP=moduleN//':'//routineN
  INTEGER :: i, stat

  failure=.false.
  my_merge=.false.
  if (present(merge)) my_merge=merge
  my_flags_only=.FALSE.
  if (present(flags_only)) my_flags_only=flags_only
  CPPrecondition(ASSOCIATED(flags_env_in),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env_in%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(flags_env_out),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env_out%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     if (flags_env_out%id_nr==flags_env_in%id_nr) return
     IF (merge) THEN
        DO i=1,flags_env_in%n_flags
           CALL cp_out_flags_set_flag(flags_env_out,flags_env_in%names(i),&
                value=flags_env_in%control_val(i),error=error)
        END DO
     ELSE
        IF (SIZE(flags_env_out%names)<flags_env_in%n_flags.OR.&
             SIZE(flags_env_out%names)>flags_env_in%n_flags+10) then
           DEALLOCATE(flags_env_out%control_val,flags_env_out%names,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           ALLOCATE(flags_env_out%names(flags_env_in%n_flags),&
                flags_env_out%control_val(flags_env_in%n_flags),stat=stat)
           CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        END IF
        flags_env_out%control_val(:flags_env_in%n_flags)=&
             flags_env_in%control_val(:flags_env_in%n_flags)
        flags_env_out%control_val(flags_env_in%n_flags+1:)=0
        flags_env_out%names(:flags_env_in%n_flags)=&
             flags_env_in%names(:flags_env_in%n_flags)
        flags_env_out%names(flags_env_in%n_flags+1:)=""
     END IF
     flags_env_out%default_val=flags_env_in%default_val
     IF (.NOT.my_flags_only) THEN
        CALL cp_logger_retain(flags_env_in%logger)
        CALL cp_logger_release(flags_env_out%logger)
        flags_env_out%logger => flags_env_in%logger
        CALL section_vals_retain(flags_env_in%input,error=error)
        CALL section_vals_release(flags_env_out%input,error=error)
        flags_env_out%input => flags_env_in%input
        flags_env_out%strict=flags_env_in%strict
     END IF
  END IF
END SUBROUTINE cp_out_flags_copy
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_unit_nr *
!!
!!   NAME
!!     cp_out_flags_unit_nr
!!
!!   FUNCTION
!!     returns a unit number for the given output (or -1 if no input should
!!     be written (if the output unit should be global only one processor
!!     recives a unit_nr > 1)
!!     You should call cp_out_flags_finished_output when you are done with
!!     the unit.
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - flags_env: the flags environment that controls the output
!!     - name: the name of the flag to check (can be a path)
!!     - extension: the extension to use for the filename
!!     - middle_name: the name to add to the project_name (defaults to
!!       the printsection name)
!!     - local: if the output should be local to each mpi task (defaults to
!!       false)
!!     - log_filename: if the generated filename (if it is needed to generate
!!       one) should be written to the standard log.
!!     - virtual: if the name can refer to a non existing printkey (defaults
!!       to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION cp_out_flags_unit_nr(flags_env, name, extension, middle_name,&
     local, log_filename,virtual,error) RESULT(res)
  type(cp_out_flags_type), pointer :: flags_env
  CHARACTER(len=*), intent(in) :: name
  CHARACTER(len=*), INTENT(in) :: extension
  CHARACTER(len=*), INTENT(in), OPTIONAL :: middle_name
  LOGICAL, INTENT(in), OPTIONAL :: local, log_filename, virtual
  type(cp_error_type), optional, intent(inout) :: error
  integer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_out_flags_unit_nr',&
        routineP=moduleN//':'//routineN
  integer :: out_flags
  TYPE(section_vals_type), pointer :: pk

  failure=.false.
  nullify(pk)
  res=-1
  CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     out_flags=cp_out_flags_should_output(flags_env,name,print_key=pk,&
          virtual=virtual,error=error)
     IF (BTEST(out_flags,cp_p_file)) THEN
        res=cp_print_key_unit_nr(flags_env%logger,pk,extension,middle_name,local,&
             log_filename,error=error)
     END IF
  END IF
END FUNCTION cp_out_flags_unit_nr
!***************************************************************************

!!****f* cp_output_handling/cp_out_flags_finished_output [1.0] *
!!
!!   NAME
!!     cp_out_flags_finished_output
!!
!!   FUNCTION
!!     should be called after you finish working with a unit obtained with
!!     cp_out_flags_unit_nr, so that the file that might have been opened
!!     can be closed.
!!
!!   NOTES
!!     closes if the corresponding filename of the printkey is
!!      not __STD_OUT__
!!
!!   INPUTS
!!     the inputs should be exactly the same of the corresponding
!!     cp_out_flags_unit_nr
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_out_flags_finished_output(unit_nr, flags_env,name, local,&
       error)
    INTEGER, INTENT(inout)                   :: unit_nr
    TYPE(cp_out_flags_type), pointer         :: flags_env
    CHARACTER(len=*), intent(in)             :: name
    LOGICAL, INTENT(in), OPTIONAL            :: local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_out_flags_finished_output', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_vals_type), pointer :: pk
    integer :: out_flags

    failure=.FALSE.
    nullify(pk)
    CPPrecondition(ASSOCIATED(flags_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(flags_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       out_flags=cp_out_flags_should_output(flags_env,name,print_key=pk,&
            virtual=.true.,error=error)
       IF (BTEST(out_flags,cp_p_file)) THEN
          CALL cp_print_key_finished_output(unit_nr,flags_env%logger,pk,local,error=error)
       ENDIF
    END IF
  END SUBROUTINE cp_out_flags_finished_output
!***************************************************************************

END MODULE cp_output_handling

