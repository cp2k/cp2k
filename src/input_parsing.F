!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/input_parsing [1.0] *
!!
!!   NAME
!!     input_parsing
!!
!!   FUNCTION
!!     routines that parse the input
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     06.2004 created
!!
!!   SOURCE
!****************************************************************************
MODULE input_parsing
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_propagate_error,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_linked_list_logical,          ONLY: cp_create,&
                                             cp_dealloc,&
                                             cp_sll_logical_type,&
                                             cp_to_array
  USE cp_linked_list_val,              ONLY: cp_sll_val_create,&
                                             cp_sll_val_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_parser,                       ONLY: cp_parser_type,&
                                             max_line_length,&
                                             parser_get_object,&
                                             parser_location,&
                                             parser_test_next_token
  USE input_keyword_types,             ONLY: keyword_describe,&
                                             keyword_type
  USE input_section_types,             ONLY: section_describe,&
                                             section_get_keyword,&
                                             section_get_keyword_index,&
                                             section_get_subsection_index,&
                                             section_type,&
                                             section_vals_add_values,&
                                             section_vals_type
  USE input_val_types,                 ONLY: logical_t,&
                                             no_t,&
                                             val_create,&
                                             val_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE output_utilities,                ONLY: print_message
!  USE cp_linked_list_int, ONLY: cp_sll_int_type, cp_create,cp_to_array,&

!       cp_dealloc

!  USE cp_linked_list_real, ONLY: cp_sll_real_type, cp_create,cp_to_array,&

!       cp_dealloc

!  USE cp_linked_list_char, ONLY: cp_sll_char_type, cp_create,cp_to_array,&

!       cp_dealloc

  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_parsing'

!***
!****************************************************************************
CONTAINS

!!****f* input_parsing/section_vals_parse *
!!
!!   NAME
!!     section_vals_parse
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE section_vals_parse(section_vals,parser,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_parse', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=max_line_length)           :: token
    INTEGER                                  :: desc_level, ik, irs, is
    LOGICAL                                  :: compatible_end, failure, &
                                                skip_description, skip_help, &
                                                whole_section
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: last_val, new_val
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: el

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
  CALL cp_error_init(suberror,stop_level=cp_fatal_level,template_error=error)
  logger => cp_error_get_logger(error)
  IF (.NOT.failure) THEN
     CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
     CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
     IF (parser%icol1<=parser%icol2) &
          CALL cp_assert(.FALSE.,cp_failure_level,&
          cp_assertion_failed,routineP,&
          "this routine must be called just after having parsed the start of the section "&
          //TRIM(parser_location(parser,error=error)),error,failure)
  END IF
  IF (.NOT.failure) THEN
     section => section_vals%section
     token=parser%input_line(parser%icol1:parser%icol2)
     CALL uppercase(token)
     IF (token==parser%section_character//section%name)&
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "this routine must be called just after having parsed the start of the section "&
          //TRIM(parser_location(parser,error=error)), error,failure)
  END IF
  IF (.NOT. failure) THEN
     skip_description=.FALSE.
     skip_help=.FALSE.
     CALL section_vals_add_values(section_vals,error=error)
     irs=SIZE(section_vals%values,2)

     IF (ASSOCIATED(section%keywords(-1)%keyword)) THEN ! reads keyword params
        keyword => section%keywords(-1)%keyword
        NULLIFY(el)
        CALL val_create_parsing(el,type_of_var=keyword%type_of_var,&
             n_var=keyword%n_var,default_value=keyword%lone_keyword_value,&
             parser=parser,error=suberror)
        CALL cp_error_propagate_error(suberror,routineP,"parsing section "//&
             TRIM(section%name)//" parameters",error=error,failure=failure)
        IF (.NOT.failure) THEN
           NULLIFY(new_val)
           CALL cp_sll_val_create(new_val,el,error=error)
           section_vals%values(-1,irs)%list => new_val
        END IF
     END IF
     DO WHILE (.NOT.failure)
        CALL parser_get_object(parser,token,newline=.TRUE.,&
             lower_to_upper=.TRUE.,error=suberror)
        CALL cp_error_propagate_error(suberror,routineP,&
             "parsing section "// TRIM(section%name),&
             error=error,failure=failure)
        IF (failure) EXIT
        IF (token(1:1)==parser%section_character) THEN
           IF (token=="&END") THEN
              ! end of section
              compatible_end=.TRUE.
              IF (parser_test_next_token(parser)/="EOS") THEN
                 CALL parser_get_object(parser,token,newline=.TRUE.,&
                      lower_to_upper=.TRUE.,error=suberror)
                 IF (token/="SECTION") compatible_end=.FALSE.
              END IF
              IF (parser_test_next_token(parser)/="EOS") THEN
                 CALL parser_get_object(parser,token,newline=.TRUE.,&
                      lower_to_upper=.TRUE.,error=suberror)
                 IF (token/=section%name) compatible_end=.FALSE.
              END IF
              IF (.NOT.compatible_end) THEN
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP,"unmatched end of section "//&
                      TRIM(parser_location(parser,error=error)),error,failure)
                 skip_description=.TRUE.
              END IF
!              call cp_error_dealloc_ref(suberror)
!              RETURN
              EXIT
           END IF
           is=section_get_subsection_index(section,token(2:),error=error)
           IF (is>0) THEN
              CALL section_vals_parse(section_vals%subs_vals(is,irs)%section_vals,&
                   parser=parser,error=suberror)
              CALL cp_error_propagate_error(suberror,routineP,&
                   "parsing subsection "//TRIM(section%subsections(is)%section%name)&
                   //" of section "// TRIM(section%name),&
                   error=error,failure=failure)
           ELSE
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"unknown subsection "//TRIM(token)//" of section "&
                   //TRIM(section%name),error,failure)
           END IF
        ELSE ! token is a keyword
           IF (token=="DESCRIBE") THEN
              skip_description=.TRUE.
              desc_level=3
              IF (parser_test_next_token(parser)=="INT") THEN
                 CALL parser_get_object(parser,desc_level)
              END IF
              whole_section=.TRUE.
              DO WHILE (parser_test_next_token(parser)=="STR")
                 whole_section=.FALSE.
                 CALL parser_get_object(parser,token,newline=.TRUE.,&
                      lower_to_upper=.TRUE.,error=suberror)
                 keyword => section_get_keyword(section,token,error=error)
                 IF (.NOT.ASSOCIATED(keyword)) THEN
                    CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                         routineP,"unknown keyword to describe "//TRIM(token)//&
                         " in section "//TRIM(section%name),error,failure)
                 ELSE
                    CALL keyword_describe(keyword,&
                         cp_logger_get_default_unit_nr(logger),&
                         desc_level,error=error)
                 END IF
              END DO
              IF (whole_section) THEN
                 CALL section_describe(section,&
                      cp_logger_get_default_unit_nr(logger),&
                      desc_level,error=error)
              END IF
           ELSE ! token is a "normal" keyword
              ik=section_get_keyword_index(section,token,error=error)
              IF (ik<1) THEN ! don't accept pseudo keyword names
                 parser%icol=parser%icol1-1 ! re-read also the actual token
                 ik=0
                 IF (.NOT.ASSOCIATED(section%keywords(0)%keyword)) THEN
                    CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                         routineP,"found an unknown keyword "//TRIM(token)//&
                         " in section "//TRIM(section%name),error,failure)
                 END IF
              END IF
              IF (.NOT.failure) THEN
                 keyword => section%keywords(ik)%keyword
                 NULLIFY(el)
                 CALL val_create_parsing(el,type_of_var=keyword%type_of_var,&
                      n_var=keyword%n_var,default_value=keyword%lone_keyword_value,&
                      parser=parser,error=suberror)
                 CALL cp_error_propagate_error(suberror,routineP,"parsing keyword "//&
                      TRIM(token)//"in section "// TRIM(section%name),&
                      error=error,failure=failure)
                 IF (.NOT.failure.AND.ASSOCIATED(el)) THEN
                    NULLIFY(new_val)
                    CALL cp_sll_val_create(new_val,el,error=error)
                    last_val => section_vals%values(-1,irs)%list
                    IF (.NOT.ASSOCIATED(last_val)) THEN
                       section_vals%values(-1,irs)%list => new_val
                    ELSE
                       DO WHILE(ASSOCIATED(last_val%rest))
                          last_val => last_val%rest
                       END DO
                       last_val%rest => new_val
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
     IF (failure) THEN
        IF (.NOT.skip_help) THEN
           CALL print_message(&
                ' You can have a description of the whole input by using the keyword '//&
                'DESCRIBE_SECTIONS in the CP2K section optionally followed by the level of output.',&
                cp_logger_get_default_unit_nr(logger,local=.TRUE.),0,0,0)
           CALL print_message(&
                ' You can describe a single section using the keyword DESCRIBE in the '//&
                'section optionally followed by the level of output and '//&
                '(IF you just want some keywords) by the keywords you want to describe',&
                cp_logger_get_default_unit_nr(logger,local=.TRUE.),0,0,0)
        END IF
        IF (.NOT.skip_description) THEN
           CALL section_describe(section,&
                cp_logger_get_default_unit_nr(logger),&
                2,error=error)
        END IF
     END IF
  END IF
  CALL cp_error_dealloc_ref(suberror)
END SUBROUTINE section_vals_parse
!***************************************************************************

!!****f* input_parsing/val_create_parsing [1.0] *
!!
!!   NAME
!!     val_create_parsing
!!
!!   FUNCTION
!!     creates a val_type object by parsing the values
!!
!!   NOTES
!!     - no_t does not create a value
!!
!!   ARGUMENTS
!!     - val: the value that will be created
!!     - type_of_var: type of the value to be created
!!     - n_var: number of values to be parsed (-1: undefined)
!!     - default_value: a default value if nothing is found (can be null)
!!     - parser: the parser from where the values should be read
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_create_parsing(val,type_of_var, n_var,default_value, parser,error)
    TYPE(val_type), POINTER                  :: val
    INTEGER, INTENT(in)                      :: type_of_var, n_var
    TYPE(val_type), OPTIONAL, POINTER        :: default_value
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_create_parsing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure, l_val
    LOGICAL, DIMENSION(:), POINTER           :: l_val_p
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_sll_logical_type), POINTER       :: l_first, l_last, l_new

!  real(kind=dp) :: r_val
! , i_val
!  character(len=default_string_length) :: c_val
!  TYPE(cp_sll_real_type),POINTER :: r_new,r_last,r_first
!  TYPE(cp_sll_int_type),POINTER :: i_new,i_last,i_first
!  TYPE(cp_sll_char_type),POINTER :: c_new,c_last,c_first
!  REAL(kind=dp),dimension(:), pointer :: r_val_p
!  INTEGER,dimension(:), pointer :: i_val_p
!  CHARACTER(len=default_string_length),dimension(:), pointer :: c_val_p

  failure=.FALSE.
  CALL cp_error_init(suberror,stop_level=cp_fatal_level,template_error=error)
  CPPrecondition(.not.ASSOCIATED(val),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     SELECT CASE(type_of_var)
     CASE(no_t)
     CASE (logical_t)
        NULLIFY(l_val_p)
        IF (parser_test_next_token(parser,error=error)=="EOS") THEN
           IF (.NOT.ASSOCIATED(default_value)) THEN
              IF (n_var<1) THEN
                 ALLOCATE(l_val_p(0),stat=stat)
                 CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                 CALL val_create(val,l_vals_ptr=l_val_p,error=error)
              ELSE
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                      routineP,"no value was given and there is no default value"//&
                      TRIM(parser_location(parser,error=error)),error,failure)
              END IF
           ELSE
              CPPrecondition(ASSOCIATED(default_value%l_val),cp_failure_level,routineP,error,failure)
              IF (.NOT.failure) THEN
                 CALL val_create(val,l_vals=default_value%l_val,error=error)
              END IF
           END IF
        ELSE
           IF (n_var<1) THEN
              NULLIFY(l_last,l_first)
              CALL parser_get_object(parser,l_val,error=error)
              CALL cp_create(l_first,l_val,error=error)
              l_last => l_first
              DO WHILE(parser_test_next_token(parser,error=error)/="EOS"&
                   .AND..NOT.failure)
                 CALL parser_get_object(parser,l_val,error=error)
                 CALL cp_create(l_new,l_val,error=error)
                 l_last%rest => l_new
                 l_last => l_new
              END DO
              l_val_p => cp_to_array(l_first,error=error)
              CALL cp_dealloc(l_first,error=error)
           ELSE
              ALLOCATE(l_val_p(n_var),stat=stat)
              CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
              IF (.NOT.failure) THEN
                 DO i=1,n_var
                    CALL parser_get_object(parser,l_val_p(i),error=suberror)
                    CALL cp_error_propagate_error(suberror,routineP,error=error,&
                         failure=failure)
                    IF (failure) EXIT
                 END DO
              END IF
           END IF
           IF (ASSOCIATED(l_val_p)) THEN
              CALL val_create(val,l_vals_ptr=l_val_p,error=error)
           END IF
        END IF
     CASE default
        CALL cp_unimplemented_error(routineP, "to do", error,cp_failure_level)
     END SELECT
  END IF
  CALL cp_error_dealloc_ref(suberror)
END SUBROUTINE val_create_parsing
!***************************************************************************

END MODULE input_parsing
