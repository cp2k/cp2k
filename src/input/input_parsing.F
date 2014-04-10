!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that parse the input
!> \par History
!>      06.2004 created
!> \author fawzi
! *****************************************************************************
MODULE input_parsing
  USE cp_linked_list_char,             ONLY: cp_create,&
                                             cp_dealloc,&
                                             cp_sll_char_type,&
                                             cp_to_array
  USE cp_linked_list_int,              ONLY: cp_create,&
                                             cp_dealloc,&
                                             cp_sll_int_type,&
                                             cp_to_array
  USE cp_linked_list_logical,          ONLY: cp_create,&
                                             cp_dealloc,&
                                             cp_sll_logical_type,&
                                             cp_to_array
  USE cp_linked_list_real,             ONLY: cp_create,&
                                             cp_dealloc,&
                                             cp_sll_real_type,&
                                             cp_to_array
  USE cp_linked_list_val,              ONLY: cp_sll_val_create,&
                                             cp_sll_val_type
  USE cp_parser_methods,               ONLY: parser_get_object,&
                                             parser_location,&
                                             parser_skip_space,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type
  USE cp_units,                        ONLY: cp_unit_compatible,&
                                             cp_unit_create,&
                                             cp_unit_desc,&
                                             cp_unit_release,&
                                             cp_unit_set_type,&
                                             cp_unit_to_cp2k1,&
                                             cp_unit_type
  USE f77_blas
  USE input_enumeration_types,         ONLY: enum_c2i,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_describe,&
                                             keyword_type
  USE input_section_types,             ONLY: section_describe,&
                                             section_get_keyword,&
                                             section_get_keyword_index,&
                                             section_get_subsection_index,&
                                             section_type,&
                                             section_vals_add_values,&
                                             section_vals_type
  USE input_val_types,                 ONLY: &
       char_t, enum_t, integer_t, lchar_t, logical_t, no_t, real_t, &
       val_create, val_type, val_write
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             max_line_length
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: print_message
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_parsing'

  PUBLIC :: section_vals_parse, val_create_parsing
!***
CONTAINS

! *****************************************************************************
!> \param root_section if the root section should be parsed (defaults to true)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_vals_parse(section_vals,parser,default_units,root_section,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_unit_set_type), POINTER          :: default_units
    LOGICAL, INTENT(in), OPTIONAL            :: root_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_parse', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=max_line_length)           :: token
    INTEGER                                  :: desc_level, handle, ik, irs, &
                                                is, nsub, output_unit
    LOGICAL :: at_end, compatible_end, failure, root_sect, skip_description, &
      skip_help, suberror_happened, whole_section
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: last_val, new_val, &
                                                previous_last, previous_list
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(val_type), POINTER                  :: el

    CALL timeset(routineN,handle)

    NULLIFY(previous_list,previous_last)

    failure=.FALSE.
    root_sect=.TRUE.
    IF (PRESENT(root_section)) root_sect=root_section

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(parser),cp_failure_level,routineP,error,failure)
    CALL cp_error_init(suberror,stop_level=cp_fatal_level,template_error=error)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    IF (.NOT.failure) THEN
       CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(parser%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (root_sect.and.parser%icol1>parser%icol2) &
            CALL cp_assert(.FALSE.,cp_failure_level,&
            cp_assertion_failed,routineP,&
            "Error 1: this routine must be called just after having parsed the start of the section "&
            //TRIM(parser_location(parser,error=error)),error,failure)
    END IF
    IF (.NOT.failure) THEN
       section => section_vals%section
       IF (root_sect) THEN
          token=parser%input_line(parser%icol1:parser%icol2)
          CALL uppercase(token)
          IF (token/=parser%section_character//section%name)&
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Error 2: this routine must be called just after having parsed the start of the section "&
               //TRIM(parser_location(parser,error=error)), error,failure)
       END IF
    END IF
    IF (.NOT.failure) THEN
       CALL cp_assert(section%repeats.OR.SIZE(section_vals%values,2)==0,&
            cp_failure_level,cp_assertion_failed,routineP,&
            "Section "//TRIM(section%name)//&
            " should not repeat "//TRIM(parser_location(parser,error=error)),&
            error,failure)
    END IF
    IF (.NOT. failure) THEN
       skip_description=.TRUE.
       skip_help=.TRUE.
       CALL section_vals_add_values(section_vals,error=error)
       irs=SIZE(section_vals%values,2)

       IF (ASSOCIATED(section%keywords(-1)%keyword)) THEN ! reads section params
          keyword => section%keywords(-1)%keyword
          NULLIFY(el)
          IF (keyword%type_of_var==lchar_t) CALL parser_skip_space(parser,error=error)
          CALL val_create_parsing(el,type_of_var=keyword%type_of_var,&
               n_var=keyword%n_var,default_value=keyword%lone_keyword_value,&
               enum=keyword%enum,unit=keyword%unit,&
               default_units=default_units,&
               parser=parser,error=suberror)
          CALL cp_error_propagate_error(suberror,routineP,"parsing section "//&
               TRIM(section%name)//" parameters",error=error)
          CALL cp_error_reset(suberror)
          IF (.NOT.failure) THEN
             NULLIFY(new_val)
             CALL cp_sll_val_create(new_val,el,error=error)
             section_vals%values(-1,irs)%list => new_val
          END IF
          NULLIFY(el)
       END IF
       DO WHILE (.NOT.failure)
          CALL parser_get_object(parser,token,newline=.TRUE.,&
               lower_to_upper=.TRUE.,at_end=at_end,error=suberror)
          CALL cp_error_propagate_error(suberror,routineP,&
               "parsing section "// TRIM(section%name),&
               error=error,failure=failure)
          CALL cp_error_reset(suberror)
          IF (at_end) THEN
             CALL cp_assert(.not.root_sect,cp_failure_level,&
                  cp_assertion_failed,routineP,&
                  "unexpected end of file while parsing section "//&
                  TRIM(section%name)//" "//TRIM(parser_location(parser,error=error)),&
                  error,failure)
             EXIT
          END IF
          IF (failure) EXIT
          IF (token(1:1)==parser%section_character) THEN
             IF (token=="&END") THEN
                ! end of section
                compatible_end=.TRUE.
                IF (parser_test_next_token(parser,error=error)/="EOL") THEN
                   CALL parser_get_object(parser,token,newline=.FALSE.,&
                        lower_to_upper=.TRUE.,error=suberror)
                   IF (token/="SECTION".and.token/=section%name) THEN
                      compatible_end=.FALSE.
                   END IF
                END IF
                IF (parser_test_next_token(parser,error=error)/="EOL") THEN
                   CALL parser_get_object(parser,token,newline=.FALSE.,&
                        lower_to_upper=.TRUE.,error=suberror)
                   IF (token/=section%name) THEN
                      PRINT *,TRIM(token),"/=",TRIM(section%name)
                      compatible_end=.FALSE.
                   END IF
                END IF
                IF (.NOT.compatible_end) THEN
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                        routineP,"non-compatible end of section "//TRIM(section%name)//" "//&
                        TRIM(parser_location(parser,error=error)),error,failure)
                   skip_description=.TRUE.
                ELSE
                   CALL cp_error_propagate_error(suberror,routineP,&
                        "parsing end of section "// TRIM(section%name),&
                        error=error)
                   CALL cp_error_reset(suberror)
                END IF
                !              call cp_error_dealloc_ref(suberror)
                !              RETURN
                EXIT
             END IF
             is=section_get_subsection_index(section,token(2:),error=error)
             IF (is>0) THEN
                !              PRINT *,"parsing subsection "//TRIM(section%subsections(is)%section%name)&
                !                   //" of section "// TRIM(section%name)
                CALL section_vals_parse(section_vals%subs_vals(is,irs)%section_vals,&
                     default_units=default_units,parser=parser,error=suberror)
                suberror_happened=.FALSE.
                CALL cp_error_propagate_error(suberror,routineP,&
                     "parsing subsection "//TRIM(section%subsections(is)%section%name)&
                     //" of section "// TRIM(section%name),&
                     error=error,failure=suberror_happened)
                CALL cp_error_reset(suberror)
                IF (.NOT. failure .AND. suberror_happened) THEN
                   skip_description=.TRUE.
                   skip_help=.TRUE.
                END IF
             ELSE
                ! this should be an error (failure)
                CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                     routineP,"unknown subsection "//TRIM(token(2:))//" of section "&
                     //TRIM(section%name),error=error,info=TRIM(token(2:)))
                nSub=1
                DO WHILE (nSub>0)
                   CALL parser_get_object(parser,token,newline=.TRUE.,&
                        lower_to_upper=.TRUE.,error=suberror)
                   CALL cp_error_propagate_error(suberror,routineP,&
                        "parsing end of section "// TRIM(section%name),&
                        error=error,failure=failure)
                   CALL cp_error_reset(suberror)
                   IF (failure) EXIT
                   IF (token(1:1)==parser%section_character) THEN
                      IF (token=="&END") THEN
                         nSub=nSub-1
                      ELSE
                         nSub=nSub+1
                      END IF
                   END IF
                END DO
             END IF
          ELSE ! token is a keyword
             IF (token=="DESCRIBE") THEN
                IF (output_unit>0) WRITE(output_unit,"(/,' ****** DESCRIPTION ******',/)")
                skip_description=.TRUE.
                desc_level=3
                IF (parser_test_next_token(parser,error=error)=="INT") THEN
                   CALL parser_get_object(parser,desc_level,error=error)
                END IF
                whole_section=.TRUE.
                DO WHILE (parser_test_next_token(parser,error=error)=="STR")
                   whole_section=.FALSE.
                   CALL parser_get_object(parser,token,newline=.FALSE.,&
                        lower_to_upper=.TRUE.,error=suberror)
                   keyword => section_get_keyword(section,token,error=error)
                   IF (.NOT.ASSOCIATED(keyword)) THEN
                      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                           routineP,"unknown keyword to describe "//TRIM(token)//&
                           " in section "//TRIM(section%name),error,failure)
                   ELSE
                      CALL keyword_describe(keyword,output_unit,desc_level,error=error)
                   END IF
                END DO
                IF (whole_section) THEN
                   CALL section_describe(section, output_unit, desc_level,hide_root=.NOT.root_sect,&
                        error=error)
                END IF
                IF (output_unit>0) WRITE(output_unit,"(/,' ****** =========== ******',/)")

             ELSE ! token is a "normal" keyword
                ik=section_get_keyword_index(section,token,error=error)
                IF (ik<1) THEN ! don't accept pseudo keyword names
                   parser%icol=parser%icol1-1 ! re-read also the actual token
                   ik=0
                   IF (.NOT.ASSOCIATED(section%keywords(0)%keyword)) THEN
                      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                           routineP,"found an unknown keyword "//TRIM(token)//&
                           " in section "//TRIM(section%name),error=error,info=TRIM(token))
                   END IF
                END IF
                IF (.NOT.failure) THEN
                   keyword => section%keywords(ik)%keyword
                   IF (ASSOCIATED(keyword)) THEN
                      NULLIFY(el)
                      IF (ik/=0.and.keyword%type_of_var==lchar_t) &
                           CALL parser_skip_space(parser,error=error)
                      CALL val_create_parsing(el,type_of_var=keyword%type_of_var,&
                           n_var=keyword%n_var,default_value=keyword%lone_keyword_value,&
                           enum=keyword%enum,unit=keyword%unit,&
                           default_units=default_units,parser=parser,error=suberror)
                      CALL cp_error_propagate_error(suberror,routineP,"parsing keyword "//&
                           TRIM(token)//" in section "// TRIM(section%name),&
                           error=error)
                      CALL cp_error_reset(suberror)
                      IF (.NOT.failure.AND.ASSOCIATED(el)) THEN
                         NULLIFY(new_val)
                         CALL cp_sll_val_create(new_val,el,error=error)
                         last_val => section_vals%values(ik,irs)%list
                         IF (.NOT.ASSOCIATED(last_val)) THEN
                            section_vals%values(ik,irs)%list => new_val
                         ELSE
                            IF (.NOT. keyword%repeats) THEN
                               PRINT *, "Keyword "//TRIM(token)//&
                                    " in section "//TRIM(section%name)//" should not repeat."
                               PRINT *, "new_val="
                               CALL val_write(el,6,keyword%unit,error=error)
                               PRINT *,"old_val="
                               DO
                                  CALL val_write(last_val%first_el,6,keyword%unit,&
                                       error=error)
                                  IF (.not.ASSOCIATED(last_val%rest)) EXIT
                                  last_val => last_val%rest
                               END DO
                            END IF
                            CALL cp_assert(keyword%repeats,cp_failure_level,cp_assertion_failed,&
                                 routineP,"Keyword "//TRIM(token)//&
                                 " in section "//TRIM(section%name)//" should not repeat.",error,failure)
                            IF (ASSOCIATED(last_val,previous_list)) THEN
                               last_val=>previous_last
                            ELSE
                               previous_list=>last_val
                            ENDIF
                            DO WHILE(ASSOCIATED(last_val%rest))
                               last_val => last_val%rest
                            END DO
                            last_val%rest => new_val
                            previous_last => new_val
                         END IF
                      END IF
                   END IF
                END IF
             END IF
          END IF
       END DO
       IF (failure.AND.(output_unit>0)) THEN
          IF (.NOT.skip_help) THEN
             CALL print_message(&
                  ' You can have a description of the input by using the keyword '//&
                  'DESCRIBE in the CP2K section optionally followed by the level of output (1-4).',&
                  output_unit,0,0,0)
             CALL print_message(&
                  ' Putting the keyword DESCRIBE in a subsection describes '//&
                  'just that section and all its subsections. '//&
                  'If you just want the description of some keywords you can '//&
                  'just add their name to the DESCRIBE directive.',&
                  output_unit,0,0,0)
          END IF
          IF (.NOT.skip_description) THEN
             CALL section_describe(section,output_unit,3,error=error)
          END IF
       END IF
    END IF
    CALL cp_error_dealloc_ref(suberror)
    CALL timestop(handle)
  END SUBROUTINE section_vals_parse

! *****************************************************************************
!> \brief creates a val_type object by parsing the values
!> \param val the value that will be created
!> \param type_of_var type of the value to be created
!> \param n_var number of values to be parsed (-1: undefined)
!> \param default_value a default value if nothing is found (can be null)
!> \param parser the parser from where the values should be read
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - no_t does not create a value
!> \author fawzi
! *****************************************************************************
  SUBROUTINE val_create_parsing(val,type_of_var, n_var,enum,&
       parser,unit,default_units,default_value,error)
    TYPE(val_type), POINTER                  :: val
    INTEGER, INTENT(in)                      :: type_of_var, n_var
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_unit_set_type), POINTER          :: default_units
    TYPE(val_type), OPTIONAL, POINTER        :: default_value
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_create_parsing', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_val, info, location
    CHARACTER(len=default_string_length), &
      DIMENSION(:), POINTER                  :: c_val_p
    INTEGER                                  :: handle, i, i_val, stat
    INTEGER, DIMENSION(:), POINTER           :: i_val_p
    LOGICAL                                  :: check, failure, l_val
    LOGICAL, DIMENSION(:), POINTER           :: l_val_p
    REAL(kind=dp)                            :: r_val
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_val_p
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_sll_char_type), POINTER          :: c_first, c_last, c_new
    TYPE(cp_sll_int_type), POINTER           :: i_first, i_last, i_new
    TYPE(cp_sll_logical_type), POINTER       :: l_first, l_last, l_new
    TYPE(cp_sll_real_type), POINTER          :: r_first, r_last, r_new

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CALL cp_error_init(suberror,stop_level=cp_fatal_level,template_error=error)
    CPPrecondition(.NOT.ASSOCIATED(val),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       SELECT CASE(type_of_var)
       CASE(no_t)
       CASE (logical_t)
          NULLIFY(l_val_p)
          IF (parser_test_next_token(parser,error=error)=="EOL") THEN
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
                DO WHILE(parser_test_next_token(parser,error=error)/="EOL"&
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
                      CALL cp_error_reset(suberror)
                      IF (failure) EXIT
                   END DO
                END IF
             END IF
             IF (ASSOCIATED(l_val_p)) THEN
                CALL val_create(val,l_vals_ptr=l_val_p,error=error)
             END IF
          END IF
       CASE (integer_t)
          NULLIFY(i_val_p)
          IF (parser_test_next_token(parser,error=error)=="EOL") THEN
             IF (.NOT.ASSOCIATED(default_value)) THEN
                IF (n_var<1) THEN
                   ALLOCATE(i_val_p(0),stat=stat)
                   CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                   CALL val_create(val,i_vals_ptr=i_val_p,error=error)
                ELSE
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                        routineP,"no value was given and there is no default value"//&
                        TRIM(parser_location(parser,error=error)),error,failure)
                END IF
             ELSE
                check = ASSOCIATED(default_value%i_val)
                CPPrecondition(check,cp_failure_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   CALL val_create(val,i_vals=default_value%i_val,error=error)
                END IF
             END IF
          ELSE
             IF (n_var<1) THEN
                NULLIFY(i_last,i_first)
                CALL parser_get_object(parser,i_val,error=error)
                CALL cp_create(i_first,i_val,error=error)
                i_last => i_first
                DO WHILE(parser_test_next_token(parser,error=error)/="EOL".AND..NOT.failure)
                   CALL parser_get_object(parser,i_val,error=error)
                   CALL cp_create(i_new,i_val,error=error)
                   i_last%rest => i_new
                   i_last => i_new
                END DO
                i_val_p => cp_to_array(i_first,error=error)
                CALL cp_dealloc(i_first,error=error)
             ELSE
                ALLOCATE(i_val_p(n_var),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   DO i=1,n_var
                      CALL parser_get_object(parser,i_val_p(i),error=suberror)
                      CALL cp_error_propagate_error(suberror,routineP,error=error,failure=failure)
                      CALL cp_error_reset(suberror)
                      IF (failure) EXIT
                   END DO
                END IF
             END IF
             IF (ASSOCIATED(i_val_p)) THEN
                CALL val_create(val,i_vals_ptr=i_val_p,error=error)
             END IF
          END IF
       CASE (real_t)
          NULLIFY(r_val_p)
          IF (parser_test_next_token(parser,error=error)=="EOL") THEN
             IF (.NOT.ASSOCIATED(default_value)) THEN
                IF (n_var<1) THEN
                   ALLOCATE(r_val_p(0),stat=stat)
                   CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                   CALL val_create(val,r_vals_ptr=r_val_p,error=error)
                ELSE
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                        routineP,"no value was given and there is no default value"//&
                        TRIM(parser_location(parser,error=error)),error,failure)
                END IF
             ELSE
                CPPrecondition(ASSOCIATED(default_value%r_val),cp_failure_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   CALL val_create(val,r_vals=default_value%r_val,error=error)
                END IF
             END IF
          ELSE
             IF (n_var<1) THEN
                NULLIFY(r_last,r_first)
                c_val = ""
                CALL get_r_val(r_val, parser, unit, default_units, c_val, error)
                CALL cp_create(r_first,r_val,error=error)
                r_last => r_first
                DO WHILE(parser_test_next_token(parser,error=error)/="EOL".AND..NOT.failure)
                   CALL get_r_val(r_val, parser, unit, default_units, c_val, error)
                   CALL cp_create(r_new,r_val,error=error)
                   r_last%rest => r_new
                   r_last => r_new
                END DO
                NULLIFY(r_last)
                r_val_p => cp_to_array(r_first,error=error)
                CALL cp_dealloc(r_first,error=error)
             ELSE
                ALLOCATE(r_val_p(n_var),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                c_val = ""
                IF (.NOT.failure) THEN
                   DO i=1,n_var
                      CALL get_r_val(r_val_p(i), parser, unit, default_units, c_val, error)
                      CALL cp_error_propagate_error(suberror,routineP,error=error,failure=failure)
                      CALL cp_error_reset(suberror)
                      IF (failure) EXIT
                   END DO
                END IF
             END IF
             IF (ASSOCIATED(r_val_p)) THEN
                CALL val_create(val,r_vals_ptr=r_val_p,error=error)
             END IF
          END IF
       CASE (char_t)
          NULLIFY(c_val_p)
          IF (parser_test_next_token(parser,error=error)=="EOL") THEN
             IF (n_var<1) THEN
                ALLOCATE(c_val_p(1),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                c_val_p(1)=' '
                CALL val_create(val,c_vals_ptr=c_val_p,error=error)
             ELSE
                IF (.NOT.ASSOCIATED(default_value)) THEN
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                        routineP,"no value was given and there is no default value"//&
                        TRIM(parser_location(parser,error=error)),error,failure)
                ELSE
                   CPPrecondition(ASSOCIATED(default_value%c_val),cp_failure_level,routineP,error,failure)
                   IF (.NOT.failure) THEN
                      CALL val_create(val,c_vals=default_value%c_val,error=error)
                   END IF
                END IF
             END IF
          ELSE
             IF (n_var<1) THEN
                CPAssert(n_var==-1,cp_failure_level,routineP,error,failure)
                NULLIFY(c_last,c_first)
                CALL parser_get_object(parser,c_val,error=error)
                CALL cp_create(c_first,c_val,error=error)
                c_last => c_first
                DO WHILE(parser_test_next_token(parser,error=error)/="EOL"&
                     .AND..NOT.failure)
                   CALL parser_get_object(parser,c_val,error=error)
                   CALL cp_create(c_new,c_val,error=error)
                   c_last%rest => c_new
                   c_last => c_new
                END DO
                c_val_p => cp_to_array(c_first,error=error)
                CALL cp_dealloc(c_first,error=error)
             ELSE
                ALLOCATE(c_val_p(n_var),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   DO i=1,n_var
                      CALL parser_get_object(parser,c_val_p(i),error=suberror)
                      CALL cp_error_propagate_error(suberror,routineP,error=error,&
                           failure=failure)
                      CALL cp_error_reset(suberror)
                      IF (failure) EXIT
                   END DO
                END IF
             END IF
             IF (ASSOCIATED(c_val_p)) THEN
                CALL val_create(val,c_vals_ptr=c_val_p,error=error)
             END IF
          END IF
       CASE (lchar_t)
          CALL cp_assert(.NOT.ASSOCIATED(default_value),cp_failure_level,cp_assertion_failed,&
               routineP,"input variables of type lchar_t cannot have a lone keyword attribute,"//&
               " no value is interpreted as empty string"//&
               TRIM(parser_location(parser,error=error)),error,failure)
          CALL cp_assert(n_var==1,cp_failure_level,cp_assertion_failed,&
               routineP,"input variables of type lchar_t cannot be repeated,"//&
               " one always represent a whole line, till the end"//&
               TRIM(parser_location(parser,error=error)),error,failure)
          IF (parser_test_next_token(parser,error=error)=="EOL") THEN
             ALLOCATE(c_val_p(1),stat=stat)
             CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
             c_val_p(1)=' '
          ELSE
             NULLIFY(c_last,c_first)
             CALL parser_get_object(parser,c_val,string_length=LEN(c_val),&
                  error=error)
             CALL cp_create(c_first,c_val,error=error)
             c_last => c_first
             DO WHILE(parser_test_next_token(parser,error=error)/="EOL"&
                  .AND..NOT.failure)
                CALL parser_get_object(parser,c_val,string_length=LEN(c_val),error=error)
                CALL cp_create(c_new,c_val,error=error)
                c_last%rest => c_new
                c_last => c_new
             END DO
             c_val_p => cp_to_array(c_first,error=error)
             CALL cp_dealloc(c_first,error=error)
          END IF
          CPPostcondition(ASSOCIATED(c_val_p),cp_failure_level,routineP,error,failure)
          CALL val_create(val,lc_vals_ptr=c_val_p,error=error)
       CASE (enum_t)
          CPPrecondition(ASSOCIATED(enum),cp_failure_level,routineP,error,failure)
          NULLIFY(i_val_p)
          IF (parser_test_next_token(parser,error=error)=="EOL") THEN
             IF (.NOT.ASSOCIATED(default_value)) THEN
                IF (n_var<1) THEN
                   ALLOCATE(i_val_p(0),stat=stat)
                   CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                   CALL val_create(val,i_vals_ptr=i_val_p,error=error)
                ELSE
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                        routineP,"no value was given and there is no default value"//&
                        TRIM(parser_location(parser,error=error)),error,failure)
                END IF
             ELSE
                CPPrecondition(ASSOCIATED(default_value%i_val),cp_failure_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   CALL val_create(val,i_vals=default_value%i_val,&
                        enum=default_value%enum,error=error)
                END IF
             END IF
          ELSE
             IF (n_var<1) THEN
                NULLIFY(i_last,i_first)
                CALL parser_get_object(parser,c_val,error=error)
                CALL cp_create(i_first,enum_c2i(enum,c_val,error=error),error=error)
                i_last => i_first
                DO WHILE(parser_test_next_token(parser,error=error)/="EOL".AND..NOT.failure)
                   CALL parser_get_object(parser,c_val,error=error)
                   CALL cp_create(i_new,enum_c2i(enum,c_val,error=error),error=error)
                   i_last%rest => i_new
                   i_last => i_new
                END DO
                i_val_p => cp_to_array(i_first,error=error)
                CALL cp_dealloc(i_first,error=error)
             ELSE
                ALLOCATE(i_val_p(n_var),stat=stat)
                CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
                IF (.NOT.failure) THEN
                   DO i=1,n_var
                      CALL parser_get_object(parser,c_val,error=suberror)
                      i_val_p(i)=enum_c2i(enum,c_val,error=suberror)
                      CALL cp_error_propagate_error(suberror,routineP,error=error,failure=failure)
                      CALL cp_error_reset(suberror)
                      IF (failure) EXIT
                   END DO
                END IF
             END IF
             IF (ASSOCIATED(i_val_p)) THEN
                CALL val_create(val,i_vals_ptr=i_val_p,enum=enum,error=error)
             END IF
          END IF
       CASE default
          CALL cp_unimplemented_error(routineP, "type "//cp_to_string(type_of_var)//&
               "unknown to the parser"//&
CPSourceFileRef,&
               error,cp_failure_level)
       END SELECT
       IF (parser_test_next_token(parser,error=error).NE."EOL") THEN
          location=TRIM(parser_location(parser,error=error))
          CALL parser_get_object(parser,info,error=error)
          CALL cp_assert(.FALSE.,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "found unexpected extra argument "//TRIM(info)//" at "//&
               location,error,failure,info=info)
       ENDIF
    END IF

    CALL cp_error_dealloc_ref(suberror)
    CALL timestop(handle)

  END SUBROUTINE val_create_parsing

! *****************************************************************************
!> \brief Reads and convert a real number from the input file
!> \param parser the parser from where the values should be read
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino - 11.2007 [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE get_r_val(r_val, parser, unit, default_units, c_val, error)
    REAL(kind=dp), INTENT(OUT)               :: r_val
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(cp_unit_type), POINTER              :: unit
    TYPE(cp_unit_set_type), POINTER          :: default_units
    CHARACTER(len=default_string_length), &
      INTENT(INOUT)                          :: c_val
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_r_val', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure
    TYPE(cp_unit_type), POINTER              :: my_unit

    failure = .FALSE.
    NULLIFY(my_unit)
    IF (ASSOCIATED(unit)) THEN
       IF ('STR'==parser_test_next_token(parser,error=error)) THEN
          CALL parser_get_object(parser,c_val,error=error)
          check = c_val(1:1)=="["
          CPPrecondition(check,cp_failure_level,routineP,error,failure)
          check = c_val(LEN_TRIM(c_val):LEN_TRIM(c_val))=="]"
          CPPrecondition(check,cp_failure_level,routineP,error,failure)
          CALL cp_unit_create(my_unit,c_val(2:LEN_TRIM(c_val)-1),error=error)
       ELSE
          IF (c_val /= "") THEN
             CALL cp_unit_create(my_unit,c_val(2:LEN_TRIM(c_val)-1),error=error)
          ELSE
             my_unit => unit
          END IF
       END IF
       CALL cp_assert(cp_unit_compatible(unit,my_unit,error=error),&
            cp_failure_level,cp_assertion_failed,routineP,"Incompatible units. Defined as ("//&
            TRIM(cp_unit_desc(unit,error=error))//") specified in input as ("//&
            TRIM(cp_unit_desc(my_unit,error=error))//"). These units are incompatible!",error,failure)
    END IF
    CALL parser_get_object(parser,r_val,error=error)
    IF (ASSOCIATED(unit)) THEN
       r_val=cp_unit_to_cp2k1(r_val,my_unit,default_units,error=error)
       IF (.NOT.(ASSOCIATED(my_unit,unit))) CALL cp_unit_release(my_unit,error=error)
    END IF

  END SUBROUTINE get_r_val

END MODULE input_parsing
