!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/input_section_types *
!!
!!   NAME
!!     input_section_types
!!
!!   FUNCTION
!!     objects that represent the structure of input sections and the data
!!     contained in an input section
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_section_types
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_a_l, cp_error_check,cp_error_propagate_error
  USE kinds, ONLY: dp, default_string_length
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE input_val_types, ONLY: val_type, val_create, val_retain, val_release,&
       val_write,val_get,no_t
  USE input_keyword_types, ONLY: keyword_type, keyword_create, keyword_retain,&
       keyword_release, description_string_length, keyword_p_type,&
       keyword_describe, keyword_describe_html, write_keyword_xml
  USE cp_linked_list_val, ONLY: cp_sll_val_type, cp_sll_val_create,&
       cp_sll_val_next, cp_sll_val_dealloc, cp_sll_val_get_length,&
       cp_sll_val_get_el_at, cp_sll_val_p_type
  USE string_utilities, only: uppercase
  USE output_utilities, only: print_message

  implicit none
  private

  INTEGER, SAVE, PRIVATE :: last_section_id=0, last_section_vals_id=0
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_section_types'

  PUBLIC :: section_p_type, section_type
  PUBLIC :: section_create,section_retain, section_release, section_describe,&
       section_get_subsection, section_get_keyword, section_get,&
       section_add_keyword, section_add_subsection, section_describe_html, section_describe_index_html
  PUBLIC :: section_get_subsection_index, section_get_keyword_index

  PUBLIC :: section_vals_p_type,section_vals_type
  PUBLIC :: section_vals_create, section_vals_retain, section_vals_release,&
       section_vals_get, section_vals_get_subs_vals, section_vals_val_get,&
       section_vals_write, section_vals_add_values, section_vals_get_subs_vals2
  PUBLIC :: write_section_xml

!!***
!****************************************************************************

!!****s* cp_keywords/section_p_type *
!!
!!   NAME
!!     section_p_type
!!
!!   FUNCTION
!!     represent a pointer to a section (to make arrays of pointers)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - section: the pointer to the section
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!***************************************************************************
  TYPE section_p_type
     TYPE(section_type), pointer :: section
  END TYPE section_p_type
!!***
!****************************************************************************

!!****s* cp_keywords/section_type *
!!
!!   NAME
!!     section_type
!!
!!   FUNCTION
!!     represent a section of the input file
!!
!!   NOTES
!!     - frozen: if the section has been frozen (and no keyword/subsections
!!       can be added)
!!     - required: if the section is required
!!     - repeats: if the section can be repeated more than once in the same
!!       context
!!     - id_nr: identification number (different in each instance)
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - n_keywords: the number of keywords in this section
!!     - name: name of the section
!!     - description: description of the section
!!     - keywords: array with the keywords of this section (might be
!!       oversized)
!!     - subsections: sections contained in this section
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!***************************************************************************
  TYPE section_type
     LOGICAL :: frozen, required, repeats
     INTEGER :: id_nr,ref_count, n_keywords, n_subsections
     character(len=default_string_length) :: name
     CHARACTER(len=description_string_length) :: description
     TYPE(keyword_p_type), DIMENSION(:), POINTER :: keywords
     TYPE(section_p_type), POINTER, DIMENSION(:) :: subsections
  END TYPE section_type
!!***
!****************************************************************************

!!****s* cp_keywords/section_vals_p_type *
!!
!!   NAME
!!     section_vals_p_type
!!
!!   FUNCTION
!!     repesents a pointer to a parsed section (to make arrays of pointers)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - section_vals the pointer to the parsed section
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!***************************************************************************
  TYPE section_vals_p_type
     TYPE(section_vals_type), pointer :: section_vals
  END TYPE section_vals_p_type
!!***
!****************************************************************************

!!****s* cp_keyword/section_vals_type [1.0] *
!!
!!   NAME
!!     section_vals_type
!!
!!   FUNCTION
!!     stores the values of a section
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
!***************************************************************************
  TYPE section_vals_type
     INTEGER :: ref_count, id_nr
     TYPE(section_type),pointer :: section
     TYPE(cp_sll_val_p_type), DIMENSION(:,:), POINTER :: values
     TYPE(section_vals_p_type), DIMENSION(:,:), POINTER :: subs_vals
  END TYPE section_vals_type
!!***
!****************************************************************************
  
contains

  
!!****f* cp_keywords/section_create *
!!
!!   NAME
!!     section_create
!!
!!   FUNCTION
!!     creates a list of keywords
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the list to be created
!!     - n_keywords: hint about the number of keywords, defaults to 10
!!     - n_sections: a hint about how many sections will be added to this
!!       structure, defaults to 0
!!     - repeats: if this section can repeat (defaults to false)
!!     - required: if this section is required (and thus it cannot be 
!!       totally default). Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_create(section,name,description,n_keywords,&
     n_subsections, repeats, required,error)
  TYPE(section_type), pointer :: section
  CHARACTER(len=*), INTENT(in) :: name,description
  INTEGER, INTENT(in), OPTIONAL :: n_keywords, n_subsections
  LOGICAL, INTENT(in), OPTIONAL :: repeats,required
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: my_n_keywords, my_n_subsections,stat,i

  failure=.false.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_n_keywords=10
     IF (PRESENT(n_keywords)) my_n_keywords=n_keywords
     my_n_subsections=0
     if (present(n_subsections)) my_n_subsections=n_subsections

     ALLOCATE(section,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     last_section_id=last_section_id+1
     section%id_nr=last_section_id
     section%ref_count=1

     section%n_keywords=0
     section%n_subsections=0
     section%name=name
     call uppercase(section%name)
     CPPrecondition(LEN_TRIM(description)<=LEN(section%description),cp_warning_level,routineP,error,failure)
     section%description=description
     section%frozen=.FALSE.
     section%required=.FALSE.
     section%repeats=.FALSE.
     IF (PRESENT(required)) section%required=required
     IF (PRESENT(repeats)) section%repeats=repeats

     ALLOCATE(section%keywords(-1:my_n_keywords),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        DO i=-1,my_n_keywords
           NULLIFY(section%keywords(i)%keyword)
        END DO
     END IF

     ALLOCATE(section%subsections(my_n_subsections),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (.NOT.failure) THEN
        DO i=1,my_n_subsections
           NULLIFY(section%subsections(i)%section)
        END DO
     END IF
  END IF
END SUBROUTINE section_create
!***************************************************************************

!!****f* cp_keywords/section_retain *
!!
!!   NAME
!!     section_retain
!!
!!   FUNCTION
!!     retains the given keyword list (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the list to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_retain(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(section%ref_count>0,cp_failure_level,routineP,error)
     section%ref_count=section%ref_count+1
  END IF
END SUBROUTINE section_retain
!***************************************************************************

!!****f* cp_keywords/section_release *
!!
!!   NAME
!!     section_release
!!
!!   FUNCTION
!!     releases the given keyword list (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the list to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE section_release(section,error)
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: i,stat

  failure=.false.
  
  IF (ASSOCIATED(section)) THEN
     CPPreconditionNoFail(section%ref_count>0,cp_failure_level,routineP,error)
     section%ref_count=section%ref_count-1
     IF (section%ref_count==0) THEN
        IF (ASSOCIATED(section%keywords)) THEN
           DO i=-1,UBOUND(section%keywords,1)
              CALL keyword_release(section%keywords(i)%keyword,error=error)
           END DO
           DEALLOCATE(section%keywords,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        section%n_keywords=0
        IF (ASSOCIATED(section%subsections)) THEN
           DO i=1,SIZE(section%subsections)
              CALL section_release(section%subsections(i)%section,error=error)
           END DO
           DEALLOCATE(section%subsections,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        DEALLOCATE(section,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     nullify(section)
  END IF
END SUBROUTINE section_release
!***************************************************************************

!!****f* cp_keywords/section_describe [1.0] *
!!
!!   NAME
!!     section_describe
!!
!!   FUNCTION
!!     prints a description of the given section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to describe
!!     - unit_nr: the unit to write to
!!     - level: the level of output: 0: just section name, 1:keywords,
!!       then see keyword_describe :-)
!!     - hide_root: if the name of the first section should be hidden
!!       (defaults to false).
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE section_describe(section, unit_nr, level, hide_root, error)
  TYPE(section_type), pointer :: section
  INTEGER, INTENT(in) :: unit_nr, level
  LOGICAL, INTENT(in),optional :: hide_root
  type(cp_error_type), optional, intent(inout) :: error

  LOGICAL :: failure, my_hide_root
  character(len=*), parameter :: routineN='section_desc',&
       routineP=moduleN//':'//routineN
  INTEGER :: isub,ikeyword

  failure=.false.
  my_hide_root=.false.
  if (present(hide_root)) my_hide_root=hide_root
  IF (ASSOCIATED(section)) THEN
     CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)

     if (.not.my_hide_root)&
          WRITE(unit_nr,"('*** section &',a,' ***')")TRIM(section%name)
     IF (level>1) THEN
        CALL print_message(TRIM(section%description),unit_nr,0,0,0)
     END IF
     IF (level>0) THEN
        IF (ASSOCIATED(section%keywords(-1)%keyword)) THEN
           CALL keyword_describe(section%keywords(-1)%keyword,unit_nr,&
                level,error=error)
        END IF
        IF (ASSOCIATED(section%keywords(0)%keyword)) THEN
           CALL keyword_describe(section%keywords(0)%keyword,unit_nr,&
                level,error=error)
        END IF
        DO ikeyword=1,section%n_keywords
           CALL keyword_describe(section%keywords(ikeyword)%keyword,unit_nr,&
                level,error=error)
        END DO
     END IF
     IF (section%n_subsections>0) THEN
        IF (.NOT.my_hide_root)&
             WRITE(unit_nr,"('** subsections **')")
        DO isub=1,section%n_subsections
           CALL section_describe(section%subsections(isub)%section,unit_nr,&
                level,error=error)
        END DO
     END IF
     IF (.NOT.my_hide_root)&
          WRITE(unit_nr,"('*** &end section ',a,' ***')")TRIM(section%name)
  ELSE
     WRITE(unit_nr,"('<section *null*>')")
  END IF
END SUBROUTINE section_describe
!***************************************************************************

RECURSIVE SUBROUTINE section_describe_html(section, prefix, unit_nr,  error)
  TYPE(section_type), pointer :: section
  CHARACTER(LEN=*), INTENT(IN) :: prefix
  INTEGER, INTENT(in) :: unit_nr 
  type(cp_error_type), optional, intent(inout) :: error

  LOGICAL :: failure, has_keywords
  CHARACTER(LEN=1000) :: local_prefix
  character(len=*), parameter :: routineN='section_desc',&
       routineP=moduleN//':'//routineN
  INTEGER :: isub,ikeyword

  failure=.false.
  IF (ASSOCIATED(section)) THEN
     local_prefix=TRIM(prefix//"~"//TRIM(section%name))
     CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
     OPEN(unit=179,FILE=TRIM(local_prefix)//".html")
     write(179,*) "<HTML><BODY>"
     write(179,*) "<HEAD><TITLE> The CP2K project : input section "//section%name//"</TITLE></HEAD>"
     write(179,*) "<H1> Section "//section%name//"</H1>"
     write(179,*) TRIM(section%description)
     write(179,*) '<BR>An <A HREF="index.html"> index of all sections </A> is available'
     write(179,*) "<H2> Subsections</H2>"
     IF (section%n_subsections>0) THEN
        write(179,*) "<UL>"
          DO isub=1,section%n_subsections
             write(179,*) '<LI><A HREF="'//TRIM(local_prefix)//"~"//TRIM(section%subsections(isub)%section%name)//".html"//'">'// &
                                           TRIM(section%subsections(isub)%section%name)//'</A>'
          END DO
       write(179,*) "</UL>"
     ELSE
       write(179,*) "None"
     ENDIF
     write(179,*) "<H2> Section keywords </H2>"
     has_keywords=ASSOCIATED(section%keywords(-1)%keyword) .OR. &
                  ASSOCIATED(section%keywords(0)%keyword) .OR. &
                  section%n_keywords>=1
     IF (has_keywords) THEN
        write(179,*) "<UL>"
           IF (ASSOCIATED(section%keywords(-1)%keyword)) THEN
              write(179,*) '<LI><A HREF="#'//TRIM(section%keywords(-1)%keyword%names(1))//'">'//&
                            TRIM(section%keywords(-1)%keyword%names(1))//"</A>"
           END IF
           IF (ASSOCIATED(section%keywords(0)%keyword)) THEN
              write(179,*) '<LI><A HREF="#'//TRIM(section%keywords(0)%keyword%names(1))//'">'//&
                            TRIM(section%keywords(0)%keyword%names(1))//"</A>"
           END IF
           DO ikeyword=1,section%n_keywords
              write(179,*) '<LI><A HREF="#'//TRIM(section%keywords(ikeyword)%keyword%names(1))//'">'//&
                            TRIM(section%keywords(ikeyword)%keyword%names(1))//"</A>"
           END DO
        write(179,*) "</UL>"
        write(179,*) "<H2> Keyword descriptions </H2>"
        write(179,*) "<TABLE>"
           IF (ASSOCIATED(section%keywords(-1)%keyword)) THEN
              CALL keyword_describe_html(section%keywords(-1)%keyword,179,error=error)
           END IF
           IF (ASSOCIATED(section%keywords(0)%keyword)) THEN
              CALL keyword_describe_html(section%keywords(0)%keyword,179,error=error)
           END IF
           DO ikeyword=1,section%n_keywords
              CALL keyword_describe_html(section%keywords(ikeyword)%keyword,179,error=error)
           END DO
        write(179,*) "</TABLE>"
     ELSE
       write(179,*) "None"
     ENDIF
     write(179,*) '<BR><hr>Back to the <A HREF="../index.html">CP2K homepage</A>'
     write(179,*) "</BODY></HTML>"
     CLOSE(unit=179)
     DO isub=1,section%n_subsections
        CALL section_describe_html(section%subsections(isub)%section,TRIM(local_prefix),unit_nr,&
                error=error)
     END DO
  END IF
END SUBROUTINE section_describe_html

RECURSIVE SUBROUTINE section_describe_index_html(section,prefix, unit_nr, error)
  TYPE(section_type), pointer :: section
  CHARACTER(LEN=*), INTENT(IN) :: prefix
  INTEGER, INTENT(in) :: unit_nr
  type(cp_error_type), optional, intent(inout) :: error

  LOGICAL :: failure 
  CHARACTER(LEN=1000) :: local_prefix
  INTEGER :: isub
  character(len=*), parameter :: routineN='section_describe_index_html',&
       routineP=moduleN//':'//routineN

  failure=.false.
  IF (ASSOCIATED(section)) THEN
     CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
     local_prefix=TRIM(prefix//"~"//TRIM(section%name))
     write(unit_nr,*) '<LI><A HREF="'//TRIM(local_prefix)//'.html">'//section%name//"</A>"
     IF (section%n_subsections>0) THEN
        write(unit_nr,*) "<UL>"
        DO isub=1,section%n_subsections
              CALL section_describe_index_html(section%subsections(isub)%section,TRIM(local_prefix),unit_nr,&
                   error=error)
        END DO
        write(unit_nr,*) "</UL>"
     ENDIF
  ENDIF
END SUBROUTINE section_describe_index_html

!***************************************************************************
!!****f* cp_keywords/section_get_subsection_index *
!!
!!   NAME
!!     section_get_subsection_index
!!
!!   FUNCTION
!!     returns the index of requested subsection (-1 if not found)
!!
!!   NOTES
!!     private utility function
!!
!!   INPUTS
!!     - section: the root section
!!     - subsection_name: the name of the subsection you want to get
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION section_get_subsection_index(section,subsection_name,error) RESULT(res)
  TYPE(section_type),pointer :: section
  CHARACTER(len=*), intent(in) :: subsection_name
  type(cp_error_type), optional, intent(inout) :: error
  INTEGER :: res

  logical :: failure
  character(len=*), parameter :: routineN='section_get_subsection_index',&
       routineP=moduleN//':'//routineN
  CHARACTER(len=default_string_length) :: upc_name
  integer :: isub

  failure=.false.

  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  res=-1
  IF (.NOT. failure) THEN
     upc_name=subsection_name
     call uppercase(upc_name)
     DO isub=1,section%n_subsections
        CPInvariant(ASSOCIATED(section%subsections(isub)%section),cp_failure_level,routineP,error,failure)
        IF (section%subsections(isub)%section%name==upc_name) THEN
           res=isub
           EXIT
        END IF
     END DO
  END IF
END FUNCTION section_get_subsection_index
!***************************************************************************

!!****f* cp_keywords/section_get_subsection *
!!
!!   NAME
!!     section_get_subsection
!!
!!   FUNCTION
!!     returns the requested subsection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the root section
!!     - subsection_name: the name of the subsection you want to get
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION section_get_subsection(section,subsection_name,error) RESULT(res)
  TYPE(section_type),pointer :: section
  CHARACTER(len=*), intent(in) :: subsection_name
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(section_type),pointer :: res
  
  character(len=*), parameter :: routineN='section_get_subsection',&
        routineP=moduleN//':'//routineN
  integer :: isub

  isub=section_get_subsection_index(section,subsection_name,error=error)
  IF (isub>0) THEN
     res => section%subsections(isub)%section
  ELSE
     nullify(res)
  END IF
END FUNCTION section_get_subsection
!***************************************************************************

!!****f* cp_keywords/section_get_keyword_index *
!!
!!   NAME
!!     section_get_keyword_index
!!
!!   FUNCTION
!!     returns the index of the requested keyword (or -2 if not found)
!!
!!   NOTES
!!     private utility function
!!
!!   INPUTS
!!     - section: the section the keyword is in
!!     - keyword_name: the keyword you are interested in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION section_get_keyword_index(section,keyword_name,error) RESULT(res)
  TYPE(section_type),pointer :: section
  CHARACTER(len=*), intent(in) :: keyword_name
  type(cp_error_type), optional, intent(inout) :: error
  INTEGER :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='section_get_keyword_index',&
        routineP=moduleN//':'//routineN
  INTEGER :: ik,in
  CHARACTER(len=default_string_length) :: upc_name

  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(section%keywords),cp_failure_level,routineP,error,failure)
  res=-2
  IF (.NOT. failure) THEN
     upc_name=keyword_name
     call uppercase(upc_name)
     DO ik=-1,0
        IF (ASSOCIATED(section%keywords(ik)%keyword)) THEN
           IF (section%keywords(ik)%keyword%names(1)==upc_name) THEN
              res = ik
           END IF
        END IF
     END DO
     IF (res==-2) THEN
        k_search_loop: DO ik=1,section%n_keywords
           CPInvariant(ASSOCIATED(section%keywords(ik)%keyword),cp_failure_level,routineP,error,failure)
           DO in=1,SIZE(section%keywords(ik)%keyword%names)
              IF (section%keywords(ik)%keyword%names(in)==upc_name) THEN
                 res = ik
                 exit k_search_loop
              END IF
           END DO
        END DO k_search_loop
     END IF
  END IF
END FUNCTION section_get_keyword_index

!!****f* cp_keywords/section_get_keyword *
!!
!!   NAME
!!     section_get_keyword
!!
!!   FUNCTION
!!     returns the requested keyword
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section the keyword is in
!!     - keyword_name: the keyword you are interested in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION section_get_keyword(section,keyword_name,error) RESULT(res)
  TYPE(section_type),pointer :: section
  CHARACTER(len=*), intent(in) :: keyword_name
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(keyword_type),pointer :: res
  
  character(len=*), parameter :: routineN='section_get_keyword',&
        routineP=moduleN//':'//routineN
  integer :: ik

  ik=section_get_keyword_index(section,keyword_name,error)
  IF (ik==-2) THEN
     nullify(res)
  ELSE
     res => section%keywords(ik)%keyword
  END IF
END FUNCTION section_get_keyword
!***************************************************************************

!!****f* cp_keywords/section_get *
!!
!!   NAME
!!     section_get
!!
!!   FUNCTION
!!     returns various attibutes of a section object
!!
!!   NOTES
!!      give direct access to keywords and subsections?
!!
!!   INPUTS
!!     - section: the section you want informations about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     For the other attributes see the section type
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_get(section,frozen, required, repeats,id_nr,ref_count, &
     name,description,error)
  TYPE(section_type), pointer :: section
  LOGICAL, INTENT(out),OPTIONAL :: frozen, required, repeats
  INTEGER, INTENT(out),OPTIONAL :: id_nr,ref_count
  CHARACTER(len=*), INTENT(out), OPTIONAL :: name
  CHARACTER(len=*), INTENT(out), OPTIONAL :: description

  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_get',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(frozen)) frozen=section%frozen
     if (present(required)) required=section%required
     if (present(repeats)) repeats=section%repeats
     if (present(id_nr)) id_nr=section%id_nr
     if (present(ref_count)) ref_count=section%ref_count
     if (present(name)) name=section%name
     if (present(description)) description=section%description
  end if
END SUBROUTINE section_get
!***************************************************************************


!!****f* cp_keywords/section_add_keyword *
!!
!!   NAME
!!     section_add_keyword
!!
!!   FUNCTION
!!     adds a keyword to the given section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to which the keyword should be added
!!     - keyword: the keyword to add
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_add_keyword(section,keyword,error)
  TYPE(section_type), pointer :: section
  TYPE(keyword_type), pointer :: keyword
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_add_keyword',&
        routineP=moduleN//':'//routineN
  TYPE(keyword_p_type), DIMENSION(:), pointer :: new_keywords
  INTEGER :: stat,i,j,k

  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.section%frozen,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL keyword_retain(keyword,error=error)
     IF (keyword%names(1)=="_SECTION_PARAMETERS_") THEN
        CALL keyword_release(section%keywords(-1)%keyword,error=error)
        section%keywords(-1)%keyword => keyword
     ELSE IF (keyword%names(1)=="_DEFAULT_KEYWORD_") then
        CALL keyword_release(section%keywords(0)%keyword,error=error)
        section%keywords(0)%keyword => keyword
     ELSE
        DO k=1,SIZE(keyword%names)
           DO i=1,section%n_keywords
              DO j=1,SIZE(section%keywords(i)%keyword%names)
                 CALL cp_assert(keyword%names(k)/=section%keywords(i)%keyword%names(j),&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "trying to add a keyword with a name ("//&
                      TRIM(keyword%names(k))//") that was already used in section "&
                      //TRIM(section%name),error,failure)
              END DO
           END DO
        END DO

        IF (UBOUND(section%keywords,1)==section%n_keywords) THEN
           ALLOCATE(new_keywords(-1:section%n_keywords+10),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO i=-1,section%n_keywords
              new_keywords(i)%keyword => section%keywords(i)%keyword
           END DO
           DO i=section%n_keywords+1,UBOUND(new_keywords,1)
              nullify(new_keywords(i)%keyword)
           END DO
           DEALLOCATE(section%keywords,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           section%keywords => new_keywords
        END IF
        section%n_keywords=section%n_keywords+1
        section%keywords(section%n_keywords)%keyword => keyword
     END IF
  END IF
END SUBROUTINE section_add_keyword
!***************************************************************************

!!****f* cp_keywords/section_add_subsection *
!!
!!   NAME
!!     section_add_subsection
!!
!!   FUNCTION
!!     adds a subsection to the given section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: to section to which you want to add a subsection
!!     - subsection: the subsection to add
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_add_subsection(section,subsection,error)
  TYPE(section_type), POINTER :: section,subsection
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='section_add_subsection',&
       routineP=moduleN//':'//routineN
  TYPE(section_p_type), DIMENSION(:), pointer :: new_subsections
  INTEGER :: stat,i

  failure=.false.

  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(subsection),cp_failure_level,routineP,error,failure)
  CPPrecondition(subsection%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (SIZE(section%subsections)<section%n_subsections+1) THEN
        ALLOCATE(new_subsections(section%n_subsections+10),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        DO i=1,section%n_subsections
           new_subsections(i)%section => section%subsections(i)%section
        END DO
        DO i=section%n_subsections+1,SIZE(new_subsections)
           nullify(new_subsections(i)%section)
        END DO
        DEALLOCATE(section%subsections,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        section%subsections => new_subsections
     END IF
     DO i=1,section%n_subsections
        CALL cp_assert(subsection%name/=section%subsections(i)%section%name,&
             cp_failure_level,cp_assertion_failed,routineP,&
             "trying to add a subsection with a name ("//&
             TRIM(subsection%name)//") that was already used in section "&
             //TRIM(section%name),error,failure)
     END DO
     CALL section_retain(subsection,error=error)
     section%n_subsections=section%n_subsections+1
     section%subsections(section%n_subsections)%section => subsection
  END IF
END SUBROUTINE section_add_subsection
!***************************************************************************

!!****f* cp_keywords/section_vals_create *
!!
!!   NAME
!!     section_vals_create
!!
!!   FUNCTION
!!     creates a object where to store the values of a section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section_vals: the parsed section that will be created
!!     - section: the structure of the section that you want to parse
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE section_vals_create(section_vals,section,error)
  TYPE(section_vals_type), pointer :: section_vals
  TYPE(section_type), pointer :: section
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_vals_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  ALLOCATE(section_vals,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     last_section_vals_id=last_section_vals_id+1
     section_vals%id_nr=last_section_vals_id
     section_vals%ref_count=1
     CALL section_retain(section,error=error)
     section_vals%section => section
     section%frozen=.true.
     ALLOCATE(section_vals%values(-1:section%n_keywords,0),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     ALLOCATE(section_vals%subs_vals(section%n_subsections,1),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO i=1,section%n_subsections
        nullify(section_vals%subs_vals(i,1)%section_vals)
        CALL section_vals_create(section_vals%subs_vals(i,1)%section_vals,&
             section=section%subsections(i)%section,error=error)
     END DO
  END IF
END SUBROUTINE section_vals_create
!***************************************************************************

!!****f* cp_keywords/section_vals_retain [1.0] *
!!
!!   NAME
!!     section_vals_retain
!!
!!   FUNCTION
!!     retains the given section values (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section_vals: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_vals_retain(section_vals,error)
  TYPE(section_vals_type), pointer :: section_vals
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_vals_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(section_vals%ref_count>0,cp_failure_level,routineP,error)
     section_vals%ref_count=section_vals%ref_count+1
  END IF
END SUBROUTINE section_vals_retain
!***************************************************************************

!!****f* cp_keywords/section_vals_release [1.0] *
!!
!!   NAME
!!     section_vals_release
!!
!!   FUNCTION
!!     releases the given object
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section_vals: the section_vals to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE section_vals_release(section_vals, error)
  TYPE(section_vals_type), pointer :: section_vals
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_vals_release',&
        routineP=moduleN//':'//routineN
  TYPE(cp_sll_val_type), pointer :: vals
  TYPE(val_type), pointer :: el
  INTEGER :: i,j, stat

  failure=.false.
  
  IF (ASSOCIATED(section_vals)) THEN
     CPPreconditionNoFail(section_vals%ref_count>0,cp_failure_level,routineP,error)
     section_vals%ref_count=section_vals%ref_count-1
     IF (section_vals%ref_count==0) THEN
        CALL section_release(section_vals%section,error=error)
        DO j=1,SIZE(section_vals%values,2)
           DO i=-1,UBOUND(section_vals%values,1)
              vals => section_vals%values(i,j)%list
              DO WHILE (cp_sll_val_next(vals,el_att=el,error=error))
                 CALL val_release(el,error=error)
              END DO
              CALL cp_sll_val_dealloc(section_vals%values(i,j)%list,error=error)
           END DO
        END DO
        DEALLOCATE(section_vals%values,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DO j=1,SIZE(section_vals%subs_vals,2)
           DO i=1,SIZE(section_vals%subs_vals,1)
              CALL section_vals_release(section_vals%subs_vals(i,j)%section_vals,&
                   error=error)
           END DO
        END DO
        DEALLOCATE(section_vals%subs_vals,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(section_vals,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE section_vals_release
!***************************************************************************

!!****f* cp_keywords/section_vals_get *
!!
!!   NAME
!!     section_vals_get
!!
!!   FUNCTION
!!     returns various attributes about the section_vals
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section_vals: the section vals you want information from
!!     - n_repetition: number of repetitions of the section
!!     - n_sub_vals_rep: number of repetitions of the subsections values
!!       (max(1,n_repetition))
!!     - explicit: if the section was explicitly present in 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     For the other arguments see the attributes of section_vals_type
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_vals_get(section_vals, ref_count, id_nr, n_repetition,&
     n_subs_vals_rep,section,explicit, error)
  TYPE(section_vals_type), pointer :: section_vals
  INTEGER, INTENT(out), OPTIONAL :: ref_count, id_nr, n_repetition, n_subs_vals_rep
  TYPE(section_type), POINTER, optional :: section
  LOGICAL, INTENT(out), optional :: explicit
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_vals_get',&
        routineP=moduleN//':'//routineN

  failure=.false.
  
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     if (present(ref_count)) ref_count=section_vals%ref_count
     if (present(id_nr)) id_nr=section_vals%id_nr
     if (present(section)) section => section_vals%section
     IF (PRESENT(n_repetition)) THEN
        n_repetition=SIZE(section_vals%values,2)
        IF (.NOT.section_vals%section%required) n_repetition=MAX(1,n_repetition)
     END IF
     IF (PRESENT(n_subs_vals_rep)) n_subs_vals_rep=SIZE(section_vals%subs_vals,2)
     IF (PRESENT(explicit)) explicit=(SIZE(section_vals%values,2)>0)
  END IF
END SUBROUTINE section_vals_get
!***************************************************************************

!!****f* cp_keywords/section_vals_get_subs_vals *
!!
!!   NAME
!!     section_vals_get_subs_vals
!!
!!   FUNCTION
!!     returns the values of the requested subsection
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section_vals: the root section
!!     - subsection_name: the name of the requested subsection
!!     - i_rep_section: index of the repetition of section_vals from which
!!       you want to extract the subsection (defaults to 1)
!!     - can_return_null: if the results can be null (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE FUNCTION section_vals_get_subs_vals(section_vals,subsection_name,&
     i_rep_section,can_return_null,error) RESULT(res)
  TYPE(section_vals_type), pointer :: section_vals
  CHARACTER(len=*), intent(in) :: subsection_name
  INTEGER, INTENT(in), optional :: i_rep_section
  LOGICAL, INTENT(in), optional :: can_return_null
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(section_vals_type), pointer :: res
  
  LOGICAL :: failure, my_can_return_null,is_path
  character(len=*), parameter :: routineN='section_vals_get_subs_vals',&
        routineP=moduleN//':'//routineN
  INTEGER :: isection,irep,my_index

  failure=.false.
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)


  my_can_return_null=.false.
  if (present(can_return_null)) my_can_return_null=can_return_null
  nullify(res)
  IF (.NOT. failure) THEN
     irep=1
     if (present(i_rep_section)) irep=i_rep_section

     ! prepare for recursive parsing of subsections. i_rep_section will be used for last section
     my_index=INDEX(subsection_name,"%")
     IF (my_index.EQ.0) THEN
        is_path=.FALSE.
        my_index=LEN(subsection_name)
     ELSE
        is_path=.TRUE.
        irep=1
        my_index=my_index-1
     ENDIF

     CPPrecondition(irep<=SIZE(section_vals%subs_vals,2),cp_failure_level,routineP,error,failure)
  END IF


  IF (.NOT.failure) THEN
     isection=section_get_subsection_index(section_vals%section,subsection_name(1:my_index),&
          error=error)
     IF (isection>0) res => section_vals%subs_vals(isection,irep)%section_vals
  END IF
  call cp_assert(ASSOCIATED(res).OR.my_can_return_null,cp_failure_level,&
       cp_assertion_failed,routineP,&
       "could not find subsection "//TRIM(subsection_name(1:my_index))//" in section "//&
       TRIM(section_vals%section%name)//" at "//&
       CPSourceFileRef,&
       error,failure)
  IF (is_path .AND. ASSOCIATED(res) ) THEN
     res=>section_vals_get_subs_vals(res,subsection_name(my_index+2:LEN(subsection_name)),&
     i_rep_section,can_return_null,error) 
  ENDIF

END FUNCTION section_vals_get_subs_vals
!***************************************************************************

!!****f* cp_keywords/section_vals_get_subs_vals2 *
!!
!!   NAME
!!     section_vals_get_subs_vals2
!!
!!   FUNCTION
!!     returns the values of the n-th non default subsection (null if no
!!     such section exists (not so many non default section))
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section_vals: the root section
!!     - i_section: index of the section
!!     - i_rep_section: index of the repetition of section_vals from which
!!       you want to extract the subsection (defaults to 1)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION section_vals_get_subs_vals2(section_vals,i_section,&
     i_rep_section,error) RESULT(res)
  TYPE(section_vals_type), pointer :: section_vals
  INTEGER, INTENT(in), OPTIONAL :: i_section,i_rep_section
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(section_vals_type), pointer :: res
  
  LOGICAL :: failure
  character(len=*), parameter :: routineN='section_vals_get_subs_vals',&
        routineP=moduleN//':'//routineN
  INTEGER :: irep,isect_att, i

  failure=.false.
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
  nullify(res)
  IF (.NOT. failure) THEN
     irep=1
     if (present(i_rep_section)) irep=i_rep_section
     CPPrecondition(irep<=SIZE(section_vals%subs_vals,2),cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     isect_att=0
     DO i=1,section_vals%section%n_subsections
        IF (SIZE(section_vals%subs_vals(i,irep)%section_vals%values,2)>0) THEN
           isect_att=isect_att+1
           IF (isect_att==i_section) THEN
              res => section_vals%subs_vals(i,irep)%section_vals
              exit
           END IF
        END IF
     END DO
  END IF
END FUNCTION section_vals_get_subs_vals2
!***************************************************************************

!!****f* cp_keywords/section_vals_add_values *
!!
!!   NAME
!!     section_vals_add_values
!!
!!   FUNCTION
!!     adds the place to store the values of a repetition of the section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section_vals: the section you want to extend
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_vals_add_values(section_vals,error)
  TYPE(section_vals_type), pointer :: section_vals
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='section_vals_add_values',&
       routineP=moduleN//':'//routineN
  INTEGER :: i,j,stat
  TYPE(section_vals_p_type), DIMENSION(:,:), pointer :: new_sps
  TYPE(cp_sll_val_p_type), DIMENSION(:,:), pointer :: new_values

  failure=.false.

  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(new_values(-1:UBOUND(section_vals%values,1),SIZE(section_vals%values,2)+1),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO j=1,SIZE(section_vals%values,2)
        DO i=-1,UBOUND(section_vals%values,1)
           new_values(i,j)%list => section_vals%values(i,j)%list
        END DO
     END DO
     DEALLOCATE(section_vals%values,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     section_vals%values => new_values
     DO i=-1,UBOUND(new_values,1)
        j=SIZE ( new_values, 2 )
        NULLIFY(new_values(i,j)%list)
     END DO

     IF (SIZE(new_values,2)>1) THEN
        ALLOCATE(new_sps(SIZE(section_vals%subs_vals,1),&
             SIZE(section_vals%subs_vals,2)+1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO j=1,SIZE(section_vals%subs_vals,2)
           DO i=1,SIZE(section_vals%subs_vals,1)
              new_sps(i,j)%section_vals => section_vals%subs_vals(i,j)%section_vals
           END DO
        END DO
        DEALLOCATE(section_vals%subs_vals,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        section_vals%subs_vals => new_sps
        DO i=1,SIZE(new_sps,1)
           j = SIZE ( new_sps, 2 )
           NULLIFY(new_sps(i,j)%section_vals)
           CALL section_vals_create(new_sps(i,SIZE(new_sps,2))%section_vals,&
                section=section_vals%section%subsections(i)%section,error=error)
        END DO
     END IF
  END IF
END SUBROUTINE section_vals_add_values
!***************************************************************************

!!****f* cp_keywords/section_vals_val_get *
!!
!!   NAME
!!     section_vals_val_get
!!
!!   FUNCTION
!!     returns the requested value
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the root section
!!     - keyword_name: the name of the keyword you want
!!     - isection: which repetition of the section you are interested in
!!       (defaults to 1)
!!     - ival: which repetition of the keyword/val you are interested in
!!       (defaults to 1)
!!     - n_val: returns number of val available
!!     - l_val,i_val,r_val,c_val: returns the logical,integer,real or
!!       character value
!!     - l_vals,i_vals,r_vals,c_vals: returns the logical,integer,real or
!!       character arrays. The val reamins the owner of the array
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_vals_val_get(section_vals,keyword_name,i_rep_section,&
     i_rep_val,n_rep_val,val,&
     l_val,i_val,r_val,c_val,l_vals,i_vals,r_vals,c_vals,error)
  TYPE(section_vals_type),pointer :: section_vals
  CHARACTER(len=*), intent(in) :: keyword_name
  INTEGER, INTENT(in), OPTIONAL :: i_rep_section, i_rep_val
  INTEGER,INTENT(out), optional :: n_rep_val
  TYPE(val_type), POINTER, OPTIONAL :: val
  LOGICAL, OPTIONAL, INTENT(out) :: l_val
  LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals
  INTEGER, OPTIONAL, INTENT(out) :: i_val
  INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals
  REAL(KIND=DP), OPTIONAL, INTENT(out) :: r_val
  REAL(KIND=DP), DIMENSION(:), OPTIONAL, POINTER :: r_vals
  CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: c_val
  CHARACTER(LEN=DEFAULT_STRING_LENGTH), DIMENSION(:), OPTIONAL, POINTER :: c_vals
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure,valRequested
  character(len=*), parameter :: routineN='section_val_get',&
        routineP=moduleN//':'//routineN
  INTEGER :: ik,irs,irk,my_index,tmp_index,len_key
  TYPE(section_vals_type),pointer :: s_vals
  TYPE(section_type),pointer ::section
  TYPE(keyword_type), POINTER :: keyword
  TYPE(cp_sll_val_type), pointer :: vals
  TYPE(val_type), POINTER :: my_val
  type(cp_error_type) :: sub_error
  
  failure=.false.
  
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)

  my_index=INDEX(keyword_name,'%')+1
  len_key=len(keyword_name)
  IF (my_index/=1) THEN
     DO
        tmp_index=INDEX(keyword_name(my_index:len_key),"%")
        IF (tmp_index==0) EXIT
        my_index=tmp_index+1
     END DO
     s_vals=>section_vals_get_subs_vals(section_vals,keyword_name(1:my_index-1),&
          error=error)
  ELSE
     s_vals=> section_vals
  END IF

  irk=1
  irs=1
  if (present(i_rep_section)) irs=i_rep_section
  if (present(i_rep_val)) irk=i_rep_val
  if (present(val)) nullify(val)
  section => s_vals%section
  valRequested=PRESENT(l_val).or.present(i_val).or.present(r_val).or.&
       PRESENT(c_val).OR.PRESENT(l_vals).or.present(i_vals).or.&
       PRESENT(r_vals).OR.PRESENT(c_vals)
  IF (.NOT. failure) THEN
     ik=section_get_keyword_index(s_vals%section,keyword_name(my_index:len_key),error=error)
     CALL cp_assert(ik/=-2,cp_failure_level,cp_assertion_failed,routineP,&
          "section "//TRIM(section%name)//" does not contain keyword "//&
          TRIM(keyword_name(my_index:len_key)),error,failure)
     IF (.NOT.failure) THEN
        keyword => section%keywords(ik)%keyword
     END IF
     CALL cp_assert(irs>0.AND.irs<=SIZE(s_vals%subs_vals,2),&
          cp_failure_level,cp_assertion_failed,routineP,&
          "section repetition requested ("//cp_to_string(irs)//&
          ") out of bounds (1:"//cp_to_string(SIZE(s_vals%subs_vals,2))&
          //")",error,failure)
     nullify(my_val)
  END IF
  IF (.NOT.failure) THEN
     IF (PRESENT(n_rep_val)) n_rep_val=0
     IF (irs<=SIZE(s_vals%values,2)) THEN ! the section was parsed
        vals => s_vals%values(ik,irs)%list
        IF (PRESENT(n_rep_val)) n_rep_val=cp_sll_val_get_length(vals,error=error)
        IF (.NOT.ASSOCIATED(vals)) THEN
! this keyword was not parsed
           IF (.NOT.keyword%required.AND. ASSOCIATED(keyword%default_value)) THEN
              my_val => keyword%default_value
              IF (PRESENT(n_rep_val)) n_rep_val=1
           END IF
        ELSE
           my_val => cp_sll_val_get_el_at(s_vals%values(ik,irs)%list,&
                irk,error=error)
        END IF
     ELSE IF (section%required) THEN
        CALL cp_assert(.FALSE.,&
             cp_failure_level,cp_assertion_failed,routineP,&
             "section "//TRIM(section%name)//&
             " is required, and needed but was not found",&
             error,failure)
     ELSE IF (.NOT.keyword%required.AND.&
             ASSOCIATED(keyword%default_value)) THEN
        IF (PRESENT(n_rep_val)) n_rep_val=1
        my_val => keyword%default_value
     END IF
     IF (PRESENT(val)) val => my_val
  END IF
  IF (.NOT.failure) THEN
     IF (valRequested) THEN
        CPPrecondition(ASSOCIATED(my_val),cp_failure_level,routineP,error,failure)
        CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
        CALL val_get(my_val,l_val=l_val,i_val=i_val,r_val=r_val,&
             c_val=c_val,l_vals=l_vals,i_vals=i_vals,r_vals=r_vals,&
             c_vals=c_vals,error=sub_error)
        CALL cp_error_propagate_error(sub_error,routineP,&
             "getting value from  keyword "//TRIM(keyword_name(my_index:len_key))//" of section "//&
             TRIM(section%name),error=error,failure=failure)
        call cp_error_dealloc_ref(sub_error)
     END IF
  END IF
END SUBROUTINE section_vals_val_get
!***************************************************************************

!!****f* cp_keywords/section_vals_write *
!!
!!   NAME
!!     section_vals_write
!!
!!   FUNCTION
!!     writes the values in the given section in a way that is suitable to
!!     the automatic parsing
!!
!!   NOTES
!!     skips required sections which weren't read
!!
!!   INPUTS
!!     - section_vals: the section to write out
!!     - unit_nr: the unit where to write to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE section_vals_write(section_vals,unit_nr,hide_root,error)
  TYPE(section_vals_type), pointer :: section_vals
  INTEGER, INTENT(in) :: unit_nr
  LOGICAL, INTENT(in), optional :: hide_root
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='section_vals_write',&
       routineP=moduleN//':'//routineN
  INTEGER :: i_rep_s,nr,isec,ik,ival,nval
  TYPE(section_type), pointer :: section
  TYPE(keyword_type), pointer :: keyword
  TYPE(section_vals_type), pointer :: sval
  type(cp_error_type) :: suberror
  TYPE(val_type), pointer :: val
  LOGICAL :: my_hide_root,defaultSection

  failure=.false.
  my_hide_root=.false.
  if (present(hide_root)) my_hide_root=hide_root

  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
  CALL cp_error_init(suberror,template_error=error)
  IF (.NOT. failure) THEN
     CALL section_vals_get(section_vals,n_repetition=nr, section=section,&
          error=error)
     DO i_rep_s=1,nr
        if (.not.my_hide_root) WRITE(unit=unit_nr,fmt="('&',a)") TRIM(section%name)
        defaultSection=(SIZE(section_vals%values,2)==0)
        IF (.NOT.(section%required.AND.defaultSection)) THEN
           DO ik=-1,section%n_keywords
              keyword => section%keywords(ik)%keyword
              IF (ASSOCIATED(keyword)) THEN
                 IF (keyword%type_of_var/=no_t) THEN
                    CALL section_vals_val_get(section_vals,keyword%names(1),&
                         i_rep_s,n_rep_val=nval,error=error)
                    DO ival=1,nval
                       CALL section_vals_val_get(section_vals,keyword%names(1),&
                            i_rep_s,i_rep_val=ival,val=val,error=suberror)
                       CALL cp_error_propagate_error(suberror, fromWhere=routineP,&
                            message="for val "// TRIM(keyword%names(1))//" in section "//&
                            TRIM(section%name),error=error, failure=failure)
                       WRITE(unit=unit_nr,fmt="(' ',a,' ')",advance="NO")&
                            TRIM(keyword%names(1))
                       CALL val_write(val,unit_nr=unit_nr,error=suberror)
                       CALL cp_error_propagate_error(suberror, fromWhere=routineP,&
                            message="for val "// TRIM(keyword%names(1))//" in section "//&
                            TRIM(section%name),error=error, failure=failure)
                       IF (failure) EXIT
                    END DO
                 END IF
              END IF
              IF (failure) EXIT
           END DO
           IF (.NOT.failure) THEN
              IF (ASSOCIATED(section_vals%subs_vals)) THEN
                 DO isec=1,SIZE(section_vals%subs_vals)
                    sval => section_vals%subs_vals(isec,i_rep_s)%section_vals
                    IF (ASSOCIATED(sval)) THEN
                       CALL section_vals_write(sval,unit_nr=unit_nr,error=error)
                    END IF
                 END DO
              END IF
           END IF
        END IF
        if (.not.my_hide_root) &
             WRITE(unit=unit_nr,fmt="('&END ',a)") TRIM(section%name)
     END DO
  END IF
  CALL cp_error_dealloc_ref(suberror)
END SUBROUTINE section_vals_write

! *****************************************************************************

  RECURSIVE SUBROUTINE write_section_xml(section,level,unit_number,error)

    TYPE(section_type), POINTER                  :: section
    INTEGER, INTENT(IN)                          :: level,unit_number
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "write_section_xml",&
                                   routineP=moduleN//':'//routineN

    INTEGER          :: i,l0,l1
    LOGICAL          :: failure

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    IF (ASSOCIATED(section)) THEN

      CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)

!     *** Indentation for current level, next level, etc. ***

      l0 = level
      l1 = level + 1

      WRITE (UNIT=unit_number,FMT="(A)")&
        REPEAT(" ",l0)//"<SECTION>",&
        REPEAT(" ",l1)//"<NAME>"//TRIM(section%name)//"</NAME>",&
        REPEAT(" ",l1)//"<DESCRIPTION>"//TRIM(section%description)//&
        "</DESCRIPTION>"

      IF (section%required) THEN
        WRITE (UNIT=unit_number,FMT="(A)") REPEAT(" ",l1)//"<REQUIRED/>"
      END IF

      IF (section%repeats) THEN
        WRITE (UNIT=unit_number,FMT="(A)") REPEAT(" ",l1)//"<REPEATS/>"
      END IF

      DO i=1,section%n_keywords
        CALL write_keyword_xml(section%keywords(i)%keyword,l1,unit_number,&
                               error)
      END DO

      DO i=1,section%n_subsections
        CALL write_section_xml(section%subsections(i)%section,l1,unit_number,&
                               error)
      END DO

      WRITE (UNIT=unit_number,FMT="(A)") REPEAT(" ",l0)//"</SECTION>"

    END IF

  END SUBROUTINE write_section_xml

! *****************************************************************************

END MODULE input_section_types
