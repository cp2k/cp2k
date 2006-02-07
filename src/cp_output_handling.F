!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


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
  USE cp_files,                        ONLY: close_file,&
                                             open_file,&
                                             move_file
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
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_output_handling'
  PUBLIC :: cp_print_key_should_output, cp_iterate, cp_add_iter_level, cp_rm_iter_level
  PUBLIC :: cp_iter_string, cp_print_key_section_create
  PUBLIC :: cp_print_key_unit_nr, cp_print_key_finished_output
  PUBLIC :: cp_print_key_generate_filename

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
            enum_c_vals=(/"on    ","off   ","silent","low   ","medium",&
            "high  ","debug "/),&
            enum_i_vals=(/ silent_print_level-1,debug_print_level+1,&
            silent_print_level, low_print_level,&
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
            usage="FILENAME ./filename ",&
            default_lc_val=my_filename, error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="__CONTROL_VAL",&
            description=' hidden parameter that controls storage, printing,...'//&
                        ' of the print_key',&
            default_i_val=cp_out_default, error=error)
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
!!     - basis_section: section that contains the printkey
!!     - print_key_path: path to the printkey- "%" between sections, and
!!       optionally a "/" and a logical flag to check). Might be empty.
!!     - used_print_key: here the print_key that was used is returned
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  FUNCTION cp_print_key_should_output(iteration_info,basis_section,&
       print_key_path,used_print_key,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(in), OPTIONAL :: print_key_path
    TYPE(section_vals_type), POINTER,OPTIONAL :: used_print_key
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_should_output', &
      routineP = moduleN//':'//routineN

    LOGICAL :: is_on, is_iter, failure, flags
    INTEGER                                  :: my_control_val,to_path,end_str
    TYPE(section_vals_type), POINTER         :: print_key
    
    failure=.false.
    res=0
    CPPrecondition(ASSOCIATED(basis_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(basis_section%ref_count>0,cp_failure_level,routineP,error,failure)
    if (present(used_print_key)) nullify(used_print_key)
    IF (failure) THEN
       IF (iteration_info%print_level>=debug_print_level) res=cp_out_default
       RETURN
    END IF

    IF (PRESENT(print_key_path)) THEN
       end_str=LEN_TRIM(print_key_path)
       to_path=INDEX(print_key_path,"/")
       IF (to_path<1) THEN
          to_path=end_str+1
       END IF
       
       IF (to_path>1) THEN
          print_key => section_vals_get_subs_vals(basis_section,&
               print_key_path(1:(to_path-1)),error=error)
       ELSE
          print_key => basis_section
       END IF
       CPPrecondition(ASSOCIATED(print_key),cp_failure_level,routineP,error,failure)
       CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (to_path+1<end_str) THEN
          CALL section_vals_val_get(print_key,print_key_path((to_path+1):end_str),&
               l_val=flags,error=error)
       ELSE
          flags=.TRUE.
       END IF
    ELSE
       print_key => basis_section
       flags=.true.
    END IF
    if (present(used_print_key)) used_print_key => print_key

    IF (.NOT.flags) RETURN
    
    CALL section_vals_val_get(print_key,"__CONTROL_VAL",&
         i_val=my_control_val,error=error)
    is_on=cp_printkey_is_on(iteration_info,print_key,error=error)
    is_iter=cp_printkey_is_iter(iteration_info,print_key,error=error)

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
!!     If the root level is 1 removes it
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
                                                shift, s_level
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
     s_level=1
!     if (iter_info%iteration(1)==1) s_level=2
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
        DO ilevel=s_level,n_rlevel
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
        DO ilevel=s_level,iter_info%n_rlevel
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
!!    Utility function that retuns a unit number to write the print key.
!!    Might open a file with a unique filename, generated from
!!    the print_key name and iteration info.
!!
!!    Normally a valid unit (>0) is returned only if cp_print_key_should_output 
!!    says that the print_key should be printed, and if the unit is global
!!    only the io node has a valid unit.
!!    So in many cases you can decide if you should print just checking if
!!    the returned units is bigger than 0.
!!
!!    IMPORTANT you should call cp_finished_output when an iteration output is
!!    finished (to immediately close the file that might have been opened)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - logger: the logger for the parallel environment, iteration info
!!       and filename generation
!!     - basis_section: section that contains the printkey
!!     - print_key_path: path to the printkey- "%" between sections, and
!!       optionally a "/" and a logical flag to check). Might be empty.
!!     - local: if the unit should be local to this task, or global to the
!!       program (defaults to false).
!!     - log_filename: if a small log with the filename should be written
!!       to the main log (defaults to false)
!!     - extension: extension to be applied to the filename (including the ".")
!!     - middle_name: name to be added to the generated filename, useful when
!!       print_key activates different distinct outputs, to be able to 
!!       distinguish them
!!     - ignore_should_output: if true always returns a valid unit (ignoring
!!       cp_print_key_should_output)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!
!!*** *********************************************************************

FUNCTION cp_print_key_generate_filename(logger,print_key,middle_name,extension,&
     my_local, error) RESULT(filename)
  IMPLICIT NONE
  TYPE(cp_logger_type), POINTER            :: logger
  TYPE(section_vals_type), POINTER         :: print_key
  CHARACTER(len=*), INTENT(in)             :: extension
  CHARACTER(len=*), INTENT(in),OPTIONAL    :: middle_name
  LOGICAL, INTENT(IN)                      :: my_local
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
  CHARACTER(len=default_path_length)       :: filename,root,postfix
  CHARACTER(len=default_string_length)     :: my_middle_name, outPath 
  LOGICAL :: build_complete_name
  INTEGER :: my_ind1, my_ind2

  CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)

  ! find out if the extension is present in the input provided name
  my_ind2 = INDEX(outPath,"/")
  IF (my_ind2 /= 0) THEN
     build_complete_name = .FALSE.
     DO WHILE (INDEX(outPath(my_ind2+1:),"/")/=0)
        my_ind2 = INDEX(outPath(my_ind2+1:),"/") + my_ind2 
     END DO
     my_ind1 = INDEX(outPath(my_ind2:),".")
  ELSE
     build_complete_name = .TRUE.
     my_ind1 = INDEX(outPath,".")     
  END IF
  IF (my_ind1 /= 0) THEN
     build_complete_name  = .FALSE.
  END IF

  IF (build_complete_name) THEN
     IF (PRESENT(middle_name)) THEN
        IF (outPath.EQ."") THEN
           my_middle_name="-"//middle_name
        ELSE
           my_middle_name="-"//TRIM(outPath)//"-"//middle_name
        END IF
     ELSE                 
        IF (outPath.EQ."") THEN
           my_middle_name="" 
        ELSE
           my_middle_name="-"//outPath
        ENDIF
     ENDIF
     IF (my_middle_name.EQ."-") my_middle_name=""
     root=TRIM(logger%iter_info%project_name)//TRIM(my_middle_name)

     ! use the cp_iter_string as a postfix
     postfix="-"//TRIM(cp_iter_string(logger%iter_info,print_key=print_key,for_file=.TRUE.,error=error))
     IF (TRIM(postfix)=="-") postfix=""
     
     ! and add the extension
     postfix=TRIM(postfix)//extension
     
  ELSE
     root=outPath
     postfix=""
  ENDIF
  ! and let the logger generate the filename
  CALL cp_logger_generate_filename(logger,res=filename,&
       root=root, postfix=postfix,local=my_local)

END FUNCTION cp_print_key_generate_filename

FUNCTION cp_print_key_unit_nr(logger, basis_section, print_key_path, extension,&
     middle_name, local, log_filename,ignore_should_output, file_form, file_position,&
     file_action, do_backup, error) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(in)             :: print_key_path
    CHARACTER(len=*), INTENT(in)             :: extension
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: middle_name, file_form, file_position,&
         file_action
    LOGICAL, INTENT(in), OPTIONAL            :: local, log_filename, &
         ignore_should_output, do_backup

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_unit_nr', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename,root,postfix
    CHARACTER(len=default_string_length)     :: my_middle_name, outPath, my_file_form,&
         my_file_position, my_file_action
    LOGICAL                                  :: failure, my_local, &
         my_should_output, found, my_do_backup
    TYPE(section_vals_type), POINTER         :: print_key

  failure=.FALSE.
  my_local=.FALSE.
  my_do_backup = .FALSE.
  found = .FALSE.
  res=-1
  my_file_form="FORMATTED"
  my_file_position="APPEND"
  my_file_action="WRITE"
  IF (PRESENT(file_form)) my_file_form = file_form
  IF (PRESENT(file_position)) my_file_position = file_position
  IF (PRESENT(file_action)) my_file_action = file_action
  IF (PRESENT(do_backup)) my_do_backup = do_backup
  NULLIFY(print_key)
  IF (PRESENT(local)) my_local=local
  CPPrecondition(ASSOCIATED(basis_section),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPrecondition(basis_section%ref_count>0,cp_failure_level,routineP,error,failure)
     CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     my_should_output=btest(cp_print_key_should_output(logger%iter_info,&
          basis_section,print_key_path,used_print_key=print_key,error=error),cp_p_file)
     IF (PRESENT(ignore_should_output)) my_should_output=my_should_output.or.ignore_should_output
     IF (.NOT.my_should_output) RETURN
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
           filename = cp_print_key_generate_filename(logger,print_key,middle_name,extension,&
                my_local,error)
           IF (my_do_backup) THEN
              INQUIRE(FILE=filename,EXIST=found)
              IF (found) THEN
                 CALL move_file(TRIM(filename), TRIM(filename)//".bak")
              END IF
           END IF
           CALL open_file(file_name=filename,file_status="UNKNOWN",&
                file_form=my_file_form,file_action=my_file_action,&
                file_position=my_file_position,unit_number=res)
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
!!     not __STD_OUT__
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
  SUBROUTINE cp_print_key_finished_output(unit_nr, logger, basis_section,&
       print_key_path,local,ignore_should_output,error)
    INTEGER, INTENT(inout)                   :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(in), OPTIONAL :: print_key_path
    LOGICAL, INTENT(in), OPTIONAL            :: local, ignore_should_output
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_finished_output', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_local, my_should_output
    CHARACTER(len=default_string_length)     :: outPath
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    my_local=.false.
    NULLIFY(print_key)
    IF (PRESENT(local)) my_local=local
    CPPrecondition(ASSOCIATED(basis_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
    CPPrecondition(basis_section%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
    my_should_output=btest(cp_print_key_should_output(logger%iter_info,basis_section,&
         print_key_path,used_print_key=print_key,error=error),cp_p_file)
    if (present(ignore_should_output)) my_should_output=my_should_output.or.ignore_should_output
    IF (my_should_output.and.(my_local.OR.&
         logger%para_env%source==logger%para_env%mepos)) THEN
        CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)
        IF (outPath.NE.'__STD_OUT__') THEN
           CPPrecondition(unit_nr>0,cp_failure_level,routineP,error,failure)
           CALL close_file(unit_nr,"KEEP")
           unit_nr=-1
        ELSE
           unit_nr=-1
        ENDIF
    END IF
    CPPostcondition(unit_nr==-1,cp_failure_level,routineP,error,failure)
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

END MODULE cp_output_handling

