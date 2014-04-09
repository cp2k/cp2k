!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief objects that represent the structure of input sections and the data
!>      contained in an input section
!> \par History
!>      06.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_section_types
  USE cp_linked_list_val,              ONLY: &
       cp_sll_val_create, cp_sll_val_dealloc, cp_sll_val_get_el_at, &
       cp_sll_val_get_length, cp_sll_val_get_rest, cp_sll_val_insert_el_at, &
       cp_sll_val_next, cp_sll_val_p_type, cp_sll_val_rm_el_at, &
       cp_sll_val_set_el_at, cp_sll_val_type
  USE f77_blas
  USE input_keyword_types,             ONLY: &
       enable_color_tags, keyword_describe, keyword_describe_html, &
       keyword_p_type, keyword_release, keyword_retain, keyword_type, &
       keyword_typo_match, write_keyword_xml
  USE input_val_types,                 ONLY: lchar_t,&
                                             no_t,&
                                             val_create,&
                                             val_duplicate,&
                                             val_get,&
                                             val_release,&
                                             val_type,&
                                             val_write
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE reference_manager,               ONLY: get_citation_key,&
                                             print_format_journal,&
                                             print_reference
  USE string_utilities,                ONLY: a2s,&
                                             compress,&
                                             substitute_special_xml_tokens,&
                                             typo_match,&
                                             uppercase
  USE termination,                     ONLY: print_message
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE, PRIVATE :: last_section_id=0, last_section_vals_id=0
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_section_types'

  PUBLIC :: section_p_type, section_type
  PUBLIC :: section_create,section_retain, section_release, section_describe,&
       section_get_subsection, section_get_keyword, section_get,&
       section_add_keyword, section_add_subsection, section_describe_html,&
       section_describe_index_html
  PUBLIC :: section_get_subsection_index, section_get_keyword_index
  PUBLIC :: section_typo_match

  PUBLIC :: section_vals_p_type,section_vals_type
  PUBLIC :: section_vals_create, section_vals_retain, section_vals_release,&
       section_vals_get, section_vals_get_subs_vals, section_vals_val_get, section_vals_list_get, &
       section_vals_write, section_vals_add_values, section_vals_get_subs_vals2,&
       section_vals_val_set, section_vals_val_unset, section_vals_get_subs_vals3,&
       section_vals_set_subs_vals, section_vals_duplicate, section_vals_remove_values,&
       section_vals_check_release
  PUBLIC :: write_section_xml

  PUBLIC :: section_get_ival,&
            section_get_ivals,&
            section_get_rvals,& 
            section_get_rval,&
            section_get_cval,&
            section_get_lval

! *****************************************************************************
!> \brief represent a pointer to a section (to make arrays of pointers)
!> \param section the pointer to the section
!> \author fawzi
! *****************************************************************************
  TYPE section_p_type
     TYPE(section_type), POINTER :: section
  END TYPE section_p_type

! *****************************************************************************
!> \brief represent a section of the input file
!> \note
!>      - frozen: if the section has been frozen (and no keyword/subsections
!>        can be added)
!>      - required: if the section is required
!>      - repeats: if the section can be repeated more than once in the same
!>        context
!>      - id_nr: identification number (different in each instance)
!>      - ref_count: reference count (see doc/ReferenceCounting.html)
!>      - n_keywords: the number of keywords in this section
!>      - name: name of the section
!>      - description: description of the section
!>      - citations: references to literature associated to this section
!>      - keywords: array with the keywords of this section (might be
!>        oversized)
!>      - subsections: sections contained in this section
!> \author fawzi
! *****************************************************************************
  TYPE section_type
     LOGICAL :: frozen, required, repeats, supported_feature
     INTEGER :: id_nr,ref_count, n_keywords, n_subsections
     CHARACTER(len=default_string_length)        :: name
     CHARACTER,DIMENSION(:),POINTER              :: description=>Null()
     INTEGER, POINTER, DIMENSION(:)              :: citations
     TYPE(keyword_p_type), DIMENSION(:), POINTER :: keywords
     TYPE(section_p_type), POINTER, DIMENSION(:) :: subsections
  END TYPE section_type

! *****************************************************************************
!> \brief repesents a pointer to a parsed section (to make arrays of pointers)
!> \param section_vals the pointer to the parsed section
!> \author fawzi
! *****************************************************************************
  TYPE section_vals_p_type
     TYPE(section_vals_type), POINTER :: section_vals
  END TYPE section_vals_p_type

! *****************************************************************************
!> \brief stores the values of a section
!> \author fawzi
! *****************************************************************************
  TYPE section_vals_type
     INTEGER :: ref_count, id_nr
     INTEGER, POINTER, DIMENSION(:)                     :: ibackup
     TYPE(section_type),POINTER                         :: section
     TYPE(cp_sll_val_p_type), DIMENSION(:,:), POINTER   :: values
     TYPE(section_vals_p_type), DIMENSION(:,:), POINTER :: subs_vals
  END TYPE section_vals_type

CONTAINS

! *****************************************************************************
!> \brief creates a list of keywords
!> \param section the list to be created
!> \param n_keywords hint about the number of keywords, defaults to 10
!> \param n_sections a hint about how many sections will be added to this
!>        structure, defaults to 0
!> \param repeats if this section can repeat (defaults to false)
!> \param required if this section is required (and thus it cannot be
!>        totally default). Defaults to false.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_create(section,name,description,n_keywords,&
       n_subsections, repeats, required, citations, supported_feature, error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: name, description
    INTEGER, INTENT(in), OPTIONAL            :: n_keywords, n_subsections
    LOGICAL, INTENT(in), OPTIONAL            :: repeats, required
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: citations
    LOGICAL, INTENT(in), OPTIONAL            :: supported_feature
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, my_n_keywords, &
                                                my_n_subsections, n, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_n_keywords=10
       IF (PRESENT(n_keywords)) my_n_keywords=n_keywords
       my_n_subsections=0
       IF (PRESENT(n_subsections)) my_n_subsections=n_subsections

       ALLOCATE(section,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       last_section_id=last_section_id+1
       section%id_nr=last_section_id
       section%ref_count=1

       section%n_keywords=0
       section%n_subsections=0
       section%name=name
       CALL uppercase(section%name)

       n = LEN_TRIM(description)
       ALLOCATE(section%description(n))
       DO i=1, n
          section%description(i) = description(i:i)
       END DO

       section%frozen=.FALSE.
       section%required=.FALSE.
       section%supported_feature=.FALSE.
       section%repeats=.FALSE.
       IF (PRESENT(required)) section%required=required
       IF (PRESENT(repeats)) section%repeats=repeats
       IF (PRESENT(supported_feature)) section%supported_feature=supported_feature

       NULLIFY(section%citations)
       IF (PRESENT(citations)) THEN
          ALLOCATE(section%citations(SIZE(citations)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          section%citations=citations
       ENDIF

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

! *****************************************************************************
!> \brief retains the given keyword list (see doc/ReferenceCounting.html)
!> \param section the list to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_retain(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(section%ref_count>0,cp_failure_level,routineP,error)
       section%ref_count=section%ref_count+1
    END IF
  END SUBROUTINE section_retain

! *****************************************************************************
!> \brief releases the given keyword list (see doc/ReferenceCounting.html)
!> \param section the list to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_release(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(section)) THEN
       CPPreconditionNoFail(section%ref_count>0,cp_failure_level,routineP,error)
       section%ref_count=section%ref_count-1
       IF (section%ref_count==0) THEN
          IF (ASSOCIATED(section%citations)) THEN
             DEALLOCATE(section%citations,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ENDIF
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
          DEALLOCATE(section%description)
          DEALLOCATE(section,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       NULLIFY(section)
    END IF
  END SUBROUTINE section_release

! *****************************************************************************
!> \brief collects additional information on the section for IO + documentation
!> \author fawzi
! *****************************************************************************
  FUNCTION get_section_info(section, error) RESULT(message)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=default_path_length)       :: message

    CHARACTER(len=*), PARAMETER :: routineN = 'get_section_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: length

    message = " "
    length  = LEN_TRIM(a2s(section%description))
    IF (length>0) THEN
       IF (section%description(length)/=".") THEN
          message = "."
       ENDIF
    ENDIF
    IF (section%repeats) THEN
       message=TRIM(message)//" This section can be repeated and "
    ELSE
       message=TRIM(message)//" This section can not be repeated and "
    END IF
    IF (section%required) THEN
       message=TRIM(message)//" can not be optional."
    ELSE
       message=TRIM(message)//" can be optional."
    END IF

  END FUNCTION get_section_info

! *****************************************************************************
!> \brief prints a description of the given section
!> \param section the section to describe
!> \param unit_nr the unit to write to
!> \param level the level of output: 0: just section name, 1:keywords,
!>        then see keyword_describe :-)
!> \param hide_root if the name of the first section should be hidden
!>        (defaults to false).
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_describe(section, unit_nr, level, hide_root, recurse,error)
    TYPE(section_type), POINTER              :: section
    INTEGER, INTENT(in)                      :: unit_nr, level
    LOGICAL, INTENT(in), OPTIONAL            :: hide_root
    INTEGER, INTENT(in), OPTIONAL            :: recurse
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_describe', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: message
    INTEGER                                  :: ikeyword, isub, my_recurse
    LOGICAL                                  :: failure, my_hide_root

    failure =.FALSE.
    IF (unit_nr>0) THEN
       my_hide_root =.FALSE.
       IF (PRESENT(hide_root)) my_hide_root=hide_root
       my_recurse=0
       IF (PRESENT(recurse)) my_recurse=recurse
       IF (ASSOCIATED(section)) THEN
          CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)

          IF (.not.my_hide_root)&
               WRITE(unit_nr,"('*** section &',a,' ***')")TRIM(section%name)
          IF (level>1) THEN
             message = get_section_info(section, error)
             CALL print_message(TRIM(a2s(section%description))//TRIM(message),unit_nr,0,0,0)
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
          IF (section%n_subsections>0 .and.my_recurse>=0) THEN
             IF (.NOT.my_hide_root)&
                  WRITE(unit_nr,"('** subsections **')")
             DO isub=1,section%n_subsections
                IF (my_recurse>0) THEN
                   CALL section_describe(section%subsections(isub)%section,unit_nr,&
                        level,recurse=my_recurse-1,error=error)
                ELSE
                   WRITE(unit_nr,"(' ',a)") section%subsections(isub)%section%name
                END IF
             END DO
          END IF
          IF (.NOT.my_hide_root)&
               WRITE(unit_nr,"('*** &end section ',a,' ***')")TRIM(section%name)
       ELSE
          WRITE(unit_nr,"(a)") '<section *null*>'
       END IF
    END IF
  END SUBROUTINE section_describe

! *****************************************************************************
  RECURSIVE SUBROUTINE section_describe_html(section, prefix, depth, unit_nr,  error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: prefix
    INTEGER, INTENT(in)                      :: depth, unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_describe_html', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=100)                       :: color_tag
    CHARACTER(LEN=1000)                      :: local_prefix
    CHARACTER(LEN=256), DIMENSION(50), SAVE  :: location, name
    CHARACTER(LEN=default_path_length)       :: message
    INTEGER                                  :: idepth, ikeyword, iref, isub
    LOGICAL                                  :: failure, has_keywords

    failure=.FALSE.
    IF (ASSOCIATED(section)) THEN
       local_prefix=TRIM(prefix//"~"//TRIM(section%name))
       CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
       OPEN(unit=179,FILE=TRIM(local_prefix)//".html")
       WRITE(179,FMT='(A)') "<HTML><BODY>"
       WRITE(179,FMT='(A)') "<HEAD><TITLE> The CP2K project : input section "//section%name//"</TITLE></HEAD>"
       WRITE(179,FMT='(A)') '<META NAME="description" content="CP2K">'
       WRITE(179,FMT='(A)') &
            '<META NAME="keywords" contents="scientific,computing,chemistry,physics,'//&
            'documentation,help,manual,Fortran,parallel">'
       color_tag='<font color="#000000">'
       IF (enable_color_tags.AND.section%supported_feature) THEN
          color_tag='<font color="#00CC00">'
       END IF
       WRITE(179,FMT='(A)') "<H1> Section "//TRIM(color_tag)//section%name//"</font> </H1>"
       WRITE(179,FMT='(A)') ' <A HREF="index.html">Index of all sections</A>. This section is located at '
       DO idepth=1,depth
          WRITE(179,FMT='(A)',ADVANCE="NO") '<A HREF="'//TRIM(location(idepth))//'">'//TRIM(name(idepth))//'</A>%'
       ENDDO
       location(depth+1)=TRIM(local_prefix)//".html"
       name(depth+1)=TRIM(section%name)
       WRITE(179,FMT='(A)') '<A HREF="'//TRIM(location(depth+1))//'">'//TRIM(name(depth+1))//'</A>.'

       message = get_section_info(section, error)
       WRITE(179,FMT='(A)') '<BR><BR>'//TRIM(a2s(section%description))//TRIM(message)

       IF (ASSOCIATED(section%citations)) THEN
          IF (SIZE(section%citations,1)>1) THEN
             WRITE(179,FMT='(A)') '<BR><BR> This section cites following references: '
          ELSE
             WRITE(179,FMT='(A)') '<BR><BR> This section cites following reference: '
          ENDIF
          DO iref=1,SIZE(section%citations,1)
             WRITE(179,FMT='(A,I0,A)') '<A HREF="references.html#reference_',section%citations(iref),'" TITLE="'
             CALL print_reference(section%citations(iref),FORMAT=print_format_journal,unit=179)
             WRITE(179,FMT='(A)') '">['//TRIM(get_citation_key(section%citations(iref)))//']</A>'
          ENDDO
       ENDIF

       WRITE(179,FMT='(A)') "<H2> Subsections</H2>"
       IF (section%n_subsections>0) THEN
          WRITE(179,FMT='(A)') "<UL>"
          DO isub=1,section%n_subsections
             color_tag='<font color="#0000FF">'
             IF (enable_color_tags.AND.section%subsections(isub)%section%supported_feature) THEN
                color_tag='<font color="#00CC00">'
             END IF
             WRITE(179,FMT='(A)') &
                  '<LI><A HREF="'//TRIM(local_prefix)//"~"//TRIM(section%subsections(isub)%section%name)//".html"//'">'// &
                  TRIM(color_tag)//TRIM(section%subsections(isub)%section%name)//'</font> </A>'
          END DO
          WRITE(179,FMT='(A)') "</UL>"
       ELSE
          WRITE(179,FMT='(A)') "None"
       ENDIF
       WRITE(179,FMT='(A)') "<H2> Section keywords </H2>"
       has_keywords=ASSOCIATED(section%keywords(-1)%keyword) .OR. &
            ASSOCIATED(section%keywords(0)%keyword) .OR. &
            section%n_keywords>=1
       IF (has_keywords) THEN
          WRITE(179,FMT='(A)') "<UL>"
          DO ikeyword=-1,section%n_keywords
             IF (ASSOCIATED(section%keywords(ikeyword)%keyword)) THEN
                color_tag='<font color="#0000FF">'
                IF (enable_color_tags.AND.section%keywords(ikeyword)%keyword%supported_feature) THEN
                   color_tag='<font color="#00CC00">'
                END IF
                WRITE(179,FMT='(A)') &
                     '<LI><A HREF="#'//TRIM(section%keywords(ikeyword)%keyword%names(1))//'">'//&
                     TRIM(color_tag)//TRIM(section%keywords(ikeyword)%keyword%names(1))//"</font> </A>"
             END IF
          END DO
          WRITE(179,FMT='(A)') "</UL>"
          WRITE(179,FMT='(A)') "<H2> Keyword descriptions </H2>"
          WRITE(179,FMT='(A)') "<TABLE>"
          DO ikeyword=-1,section%n_keywords
             IF (ASSOCIATED(section%keywords(ikeyword)%keyword)) THEN
                CALL keyword_describe_html(section%keywords(ikeyword)%keyword,179,error=error)
             END IF
          END DO
          WRITE(179,FMT='(A)') "</TABLE>"
       ELSE
          WRITE(179,FMT='(A)') "None"
       ENDIF
       WRITE(179,FMT='(A)') &
            '<BR><hr>Back to the <A HREF="http://www.cp2k.org/">CP2K homepage</A> or '// &
            'the latest version of <A HREF="http://manual.cp2k.org/trunk/">this manual</A>'
       WRITE(179,FMT='(A)') "</BODY></HTML>"
       CLOSE(unit=179)
       DO isub=1,section%n_subsections
          CALL section_describe_html(section%subsections(isub)%section,TRIM(local_prefix),depth+1,unit_nr,&
               error=error)
       END DO
    END IF
  END SUBROUTINE section_describe_html

! *****************************************************************************
  RECURSIVE SUBROUTINE section_describe_index_html(section,prefix, unit_nr, error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: prefix
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_describe_index_html', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=100)                       :: color_tag
    CHARACTER(LEN=1000)                      :: local_prefix
    INTEGER                                  :: isub
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(section)) THEN
       CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
       local_prefix=TRIM(prefix//"~"//TRIM(section%name))
       color_tag='<font color="#0000FF">'
       IF (enable_color_tags.AND.section%supported_feature) THEN
          color_tag='<font color="#00CC00">'
       END IF
       WRITE(unit_nr,FMT='(A)') &
            '<LI><A HREF="'//TRIM(local_prefix)//'.html">'//TRIM(color_tag)//TRIM(section%name)//"</font></A>"
       IF (section%n_subsections>0) THEN
          WRITE(unit_nr,FMT='(A)') "<UL>"
          DO isub=1,section%n_subsections
             CALL section_describe_index_html(section%subsections(isub)%section,TRIM(local_prefix),unit_nr,&
                  error=error)
          END DO
          WRITE(unit_nr,FMT='(A)') "</UL>"
       ENDIF
    ENDIF
  END SUBROUTINE section_describe_index_html

! *****************************************************************************
!> \brief returns the index of requested subsection (-1 if not found)
!> \param section the root section
!> \param subsection_name the name of the subsection you want to get
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      private utility function
!> \author fawzi
! *****************************************************************************
  FUNCTION section_get_subsection_index(section,subsection_name,error) RESULT(res)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: subsection_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_get_subsection_index', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: upc_name
    INTEGER                                  :: isub
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
    res=-1
    IF (.NOT. failure) THEN
       upc_name=subsection_name
       CALL uppercase(upc_name)
       DO isub=1,section%n_subsections
          CPInvariant(ASSOCIATED(section%subsections(isub)%section),cp_failure_level,routineP,error,failure)
          IF (section%subsections(isub)%section%name==upc_name) THEN
             res=isub
             EXIT
          END IF
       END DO
    END IF
  END FUNCTION section_get_subsection_index

! *****************************************************************************
!> \brief returns the requested subsection
!> \param section the root section
!> \param subsection_name the name of the subsection you want to get
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION section_get_subsection(section,subsection_name,error) RESULT(res)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: subsection_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_type), POINTER              :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_get_subsection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isub

    isub=section_get_subsection_index(section,subsection_name,error=error)
    IF (isub>0) THEN
       res => section%subsections(isub)%section
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION section_get_subsection

! *****************************************************************************
!> \brief returns the index of the requested keyword (or -2 if not found)
!> \param section the section the keyword is in
!> \param keyword_name the keyword you are interested in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      private utility function
!> \author fawzi
! *****************************************************************************
  FUNCTION section_get_keyword_index(section,keyword_name,error) RESULT(res)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_get_keyword_index', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, in
    CHARACTER(len=default_string_length)     :: upc_name
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(section%keywords),cp_failure_level,routineP,error,failure)
    res=-2
    IF (.NOT. failure) THEN
       upc_name=keyword_name
       CALL uppercase(upc_name)
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
                   EXIT k_search_loop
                END IF
             END DO
          END DO k_search_loop
       END IF
    END IF
  END FUNCTION section_get_keyword_index

! *****************************************************************************
!> \brief returns the requested keyword
!> \param section the section the keyword is in
!> \param keyword_name the keyword you are interested in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE FUNCTION section_get_keyword(section,keyword_name,error) RESULT(res)
    TYPE(section_type), POINTER              :: section
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(keyword_type), POINTER              :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_get_keyword', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, my_index
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (INDEX(keyword_name,"%")/=0) THEN
       my_index = INDEX(keyword_name,"%") + 1
       CPPrecondition(ASSOCIATED(section%subsections),cp_failure_level,routineP,error,failure)
       DO ik = LBOUND(section%subsections,1), UBOUND(section%subsections,1)
          IF (section%subsections(ik)%section%name==keyword_name(1:my_index-2)) EXIT
       END DO
       CPPrecondition(ik<=UBOUND(section%subsections,1),cp_failure_level,routineP,error,failure)
       res => section_get_keyword(section%subsections(ik)%section,keyword_name(my_index:),error)
    ELSE
       ik=section_get_keyword_index(section,keyword_name,error)
       IF (ik==-2) THEN
          NULLIFY(res)
       ELSE
          res => section%keywords(ik)%keyword
       END IF
    END IF
  END FUNCTION section_get_keyword

! *****************************************************************************
!> \brief returns various attibutes of a section object
!> \param section the section you want informations about
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For the other attributes see the section type
!> \note
!>       give direct access to keywords and subsections?
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_get(section,frozen, required, repeats,id_nr,ref_count, &
       name,description,citations,error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, INTENT(out), OPTIONAL           :: frozen, required, repeats
    INTEGER, INTENT(out), OPTIONAL           :: id_nr, ref_count
    CHARACTER(len=*), INTENT(out), OPTIONAL  :: name, description
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: citations
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(frozen)) frozen=section%frozen
       IF (PRESENT(required)) required=section%required
       IF (PRESENT(repeats)) repeats=section%repeats
       IF (PRESENT(id_nr)) id_nr=section%id_nr
       IF (PRESENT(ref_count)) ref_count=section%ref_count
       IF (PRESENT(name)) name=section%name
       IF (PRESENT(description)) description=a2s(section%description)
       IF (PRESENT(citations)) citations=>section%citations
    END IF
  END SUBROUTINE section_get

! *****************************************************************************
!> \brief adds a keyword to the given section
!> \param section the section to which the keyword should be added
!> \param keyword the keyword to add
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_add_keyword(section,keyword,error)
    TYPE(section_type), POINTER              :: section
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_add_keyword', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, stat
    LOGICAL                                  :: failure
    TYPE(keyword_p_type), DIMENSION(:), &
      POINTER                                :: new_keywords

    failure=.FALSE.

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
       ELSE IF (keyword%names(1)=="_DEFAULT_KEYWORD_") THEN
          CALL keyword_release(section%keywords(0)%keyword,error=error)
          section%keywords(0)%keyword => keyword
       ELSE
          DO k=1,SIZE(keyword%names)
             DO i=1,section%n_keywords
                DO j=1,SIZE(section%keywords(i)%keyword%names)
                   IF (keyword%names(k)==section%keywords(i)%keyword%names(j)) THEN
                      CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
                           "trying to add a keyword with a name ("//&
                           TRIM(keyword%names(k))//") that was already used in section "&
                           //TRIM(section%name),error,failure)
                   ENDIF
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
                NULLIFY(new_keywords(i)%keyword)
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

! *****************************************************************************
!> \brief adds a subsection to the given section
!> \param section to section to which you want to add a subsection
!> \param subsection the subsection to add
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_add_subsection(section,subsection,error)
    TYPE(section_type), POINTER              :: section, subsection
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_add_subsection', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(section_p_type), DIMENSION(:), &
      POINTER                                :: new_subsections

    failure=.FALSE.

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
             NULLIFY(new_subsections(i)%section)
          END DO
          DEALLOCATE(section%subsections,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          section%subsections => new_subsections
       END IF
       DO i=1,section%n_subsections
          IF(subsection%name==section%subsections(i)%section%name)&
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "trying to add a subsection with a name ("//&
               TRIM(subsection%name)//") that was already used in section "&
               //TRIM(section%name),error,failure)
       END DO
       CALL section_retain(subsection,error=error)
       section%n_subsections=section%n_subsections+1
       section%subsections(section%n_subsections)%section => subsection
    END IF
  END SUBROUTINE section_add_subsection

! *****************************************************************************
!> \brief creates a object where to store the values of a section
!> \param section_vals the parsed section that will be created
!> \param section the structure of the section that you want to parse
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_vals_create(section_vals,section,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    ALLOCATE(section_vals,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_section_vals_id=last_section_vals_id+1
       section_vals%id_nr=last_section_vals_id
       section_vals%ref_count=1
       CALL section_retain(section,error=error)
       section_vals%section => section
       section%frozen=.TRUE.
       ALLOCATE(section_vals%values(-1:section%n_keywords,0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(section_vals%subs_vals(section%n_subsections,1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i=1,section%n_subsections
          NULLIFY(section_vals%subs_vals(i,1)%section_vals)
          CALL section_vals_create(section_vals%subs_vals(i,1)%section_vals,&
               section=section%subsections(i)%section,error=error)
       END DO
       NULLIFY(section_vals%ibackup)
    END IF
  END SUBROUTINE section_vals_create

! *****************************************************************************
!> \brief retains the given section values (see doc/ReferenceCounting.html)
!> \param section_vals the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_retain(section_vals,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(section_vals%ref_count>0,cp_failure_level,routineP,error)
       section_vals%ref_count=section_vals%ref_count+1
    END IF
  END SUBROUTINE section_vals_retain

! *****************************************************************************
!> \brief releases the given object
!> \param section_vals the section_vals to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_vals_release(section_vals, error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: vals
    TYPE(val_type), POINTER                  :: el

    failure=.FALSE.

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
          IF (ASSOCIATED(section_vals%ibackup)) THEN
             DEALLOCATE(section_vals%ibackup,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          DEALLOCATE(section_vals,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE section_vals_release

! *****************************************************************************
!> \brief returns various attributes about the section_vals
!> \param section_vals the section vals you want information from
!> \param n_repetition number of repetitions of the section
!> \param n_sub_vals_rep number of repetitions of the subsections values
!>        (max(1,n_repetition))
!> \param explicit if the section was explicitly present in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>
!>      For the other arguments see the attributes of section_vals_type
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_get(section_vals, ref_count, id_nr, n_repetition,&
       n_subs_vals_rep,section,explicit, error)
    TYPE(section_vals_type), POINTER         :: section_vals
    INTEGER, INTENT(out), OPTIONAL           :: ref_count, id_nr, &
                                                n_repetition, n_subs_vals_rep
    TYPE(section_type), OPTIONAL, POINTER    :: section
    LOGICAL, INTENT(out), OPTIONAL           :: explicit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(ref_count)) ref_count=section_vals%ref_count
       IF (PRESENT(id_nr)) id_nr=section_vals%id_nr
       IF (PRESENT(section)) section => section_vals%section
       IF (PRESENT(n_repetition)) THEN
          n_repetition=SIZE(section_vals%values,2)
          IF (.NOT.section_vals%section%required) n_repetition=MAX(1,n_repetition)
       END IF
       IF (PRESENT(n_subs_vals_rep)) n_subs_vals_rep=SIZE(section_vals%subs_vals,2)
       IF (PRESENT(explicit)) explicit=(SIZE(section_vals%values,2)>0)
    END IF
  END SUBROUTINE section_vals_get

! *****************************************************************************
!> \brief returns the values of the requested subsection
!> \param section_vals the root section
!> \param subsection_name the name of the requested subsection
!> \param i_rep_section index of the repetition of section_vals from which
!>        you want to extract the subsection (defaults to 1)
!> \param can_return_null if the results can be null (defaults to false)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  RECURSIVE FUNCTION section_vals_get_subs_vals(section_vals,subsection_name,&
       i_rep_section,can_return_null,error) RESULT(res)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: subsection_name
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section
    LOGICAL, INTENT(in), OPTIONAL            :: can_return_null
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_get_subs_vals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: irep, isection, my_index
    LOGICAL                                  :: failure, is_path, &
                                                my_can_return_null

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)

    my_can_return_null=.FALSE.
    IF (PRESENT(can_return_null)) my_can_return_null=can_return_null
    NULLIFY(res)
    IF (.NOT. failure) THEN
       irep=1
       IF (PRESENT(i_rep_section)) irep=i_rep_section

       ! prepare for recursive parsing of subsections. i_rep_section will be used for last section
       my_index=INDEX(subsection_name,"%")
       IF (my_index.EQ.0) THEN
          is_path=.FALSE.
          my_index=LEN_TRIM(subsection_name)
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
    IF(.NOT.(ASSOCIATED(res).OR.my_can_return_null))&
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
         "could not find subsection "//TRIM(subsection_name(1:my_index))//" in section "//&
         TRIM(section_vals%section%name)//" at "//&
         CPSourceFileRef,&
         error,failure)
    IF (is_path .AND. ASSOCIATED(res) ) THEN
       res=>section_vals_get_subs_vals(res,subsection_name(my_index+2:LEN_TRIM(subsection_name)),&
            i_rep_section,can_return_null,error)
    ENDIF

  END FUNCTION section_vals_get_subs_vals

! *****************************************************************************
!> \brief returns the values of the n-th non default subsection (null if no
!>      such section exists (not so many non default section))
!> \param section_vals the root section
!> \param i_section index of the section
!> \param i_rep_section index of the repetition of section_vals from which
!>        you want to extract the subsection (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION section_vals_get_subs_vals2(section_vals,i_section,i_rep_section,error) RESULT(res)
    TYPE(section_vals_type), POINTER         :: section_vals
    INTEGER, INTENT(in)                      :: i_section
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_get_subs_vals2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, irep, isect_att
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    NULLIFY(res)
    IF (.NOT. failure) THEN
       irep=1
       IF (PRESENT(i_rep_section)) irep=i_rep_section
       CPPrecondition(irep<=SIZE(section_vals%subs_vals,2),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       isect_att=0
       DO i=1,section_vals%section%n_subsections
          IF (SIZE(section_vals%subs_vals(i,irep)%section_vals%values,2)>0) THEN
             isect_att=isect_att+1
             IF (isect_att==i_section) THEN
                res => section_vals%subs_vals(i,irep)%section_vals
                EXIT
             END IF
          END IF
       END DO
    END IF
  END FUNCTION section_vals_get_subs_vals2

! *****************************************************************************
!> \brief returns the values of the n-th non default subsection (null if no
!>      such section exists (not so many non default section))
!> \param section_vals the root section
!> \param i_section index of the section
!> \param i_rep_section index of the repetition of section_vals from which
!>        you want to extract the subsection (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION section_vals_get_subs_vals3(section_vals,subsection_name,&
       i_rep_section,error) RESULT(res)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(LEN=*), INTENT(IN)             :: subsection_name
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(section_vals_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_get_subs_vals3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_section, irep
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    NULLIFY(res)
    IF (.NOT. failure) THEN
       irep=1
       IF (PRESENT(i_rep_section)) irep=i_rep_section
       CPPrecondition(irep<=SIZE(section_vals%subs_vals,2),cp_failure_level,routineP,error,failure)
    END IF
    i_section = section_get_subsection_index(section_vals%section,subsection_name,error)
    IF (.NOT.failure) THEN
       res => section_vals%subs_vals(i_section,irep)%section_vals
    END IF
  END FUNCTION section_vals_get_subs_vals3

! *****************************************************************************
!> \brief adds the place to store the values of a repetition of the section
!> \param section_vals the section you want to extend
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_add_values(section_vals,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_add_values', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_p_type), &
      DIMENSION(:, :), POINTER               :: new_values
    TYPE(section_vals_p_type), &
      DIMENSION(:, :), POINTER               :: new_sps

    failure=.FALSE.

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
       j=SIZE ( new_values, 2 )
       DO i=-1,UBOUND(new_values,1)
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
          j = SIZE ( new_sps, 2 )
          DO i=1,SIZE(new_sps,1)
             NULLIFY(new_sps(i,j)%section_vals)
             CALL section_vals_create(new_sps(i,SIZE(new_sps,2))%section_vals,&
                  section=section_vals%section%subsections(i)%section,error=error)
          END DO
       END IF
    END IF
  END SUBROUTINE section_vals_add_values

! *****************************************************************************
!> \brief removes the values of a repetition of the section
!> \param section_vals the section you want to extend
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_remove_values(section_vals,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_remove_values', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, stat
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_p_type), &
      DIMENSION(:, :), POINTER               :: new_values
    TYPE(cp_sll_val_type), POINTER           :: vals
    TYPE(val_type), POINTER                  :: el

    failure=.FALSE.

    IF (ASSOCIATED(section_vals)) THEN
       CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          NULLIFY(el, vals)
          ! Allocate a null 0 dimension array of values
          ALLOCATE(new_values(-1:section_vals%section%n_keywords,0),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Release old values
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
          section_vals%values => new_values
       END IF
    END IF
  END SUBROUTINE section_vals_remove_values

! These accessor functions can be used instead of passing a variable
! in the parameter list of a subroutine call. This should make the
! code a lot simpler. See xc_rho_set_and_dset_create in xc.F as
! an example.

! *****************************************************************************
  FUNCTION section_get_cval(section_vals,keyword_name,error) RESULT(res)

    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=default_string_length)     :: res

    CALL section_vals_val_get(section_vals, keyword_name, c_val=res, error=error)

  END FUNCTION section_get_cval

! *****************************************************************************
  FUNCTION section_get_rval(section_vals,keyword_name,error) RESULT(res)

    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp)                            :: res

    CALL section_vals_val_get(section_vals, keyword_name, r_val=res, error=error)

  END FUNCTION section_get_rval

! *****************************************************************************
  FUNCTION section_get_rvals(section_vals,keyword_name,error) RESULT(res)

    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(kind=dp), DIMENSION(:), POINTER     :: res

    CALL section_vals_val_get(section_vals, keyword_name, r_vals=res, error=error)

  END FUNCTION section_get_rvals

! *****************************************************************************
  FUNCTION section_get_ival(section_vals,keyword_name,error) RESULT(res)

    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CALL section_vals_val_get(section_vals, keyword_name, i_val=res, error=error)

  END FUNCTION section_get_ival

! *****************************************************************************
  FUNCTION section_get_ivals(section_vals,keyword_name,error) RESULT(res)

    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(:), POINTER           :: res

    CALL section_vals_val_get(section_vals, keyword_name, i_vals=res, error=error)

  END FUNCTION section_get_ivals

! *****************************************************************************
  FUNCTION section_get_lval(section_vals,keyword_name,error) RESULT(res)

    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res

    CALL section_vals_val_get(section_vals, keyword_name, l_val=res, error=error)

  END FUNCTION section_get_lval

! *****************************************************************************
!> \brief returns the requested value
!> \param section the root section
!> \param keyword_name the name of the keyword you want
!> \param isection which repetition of the section you are interested in
!>        (defaults to 1)
!> \param ival which repetition of the keyword/val you are interested in
!>        (defaults to 1)
!> \param n_rep_val returns number of val available
!> \param l_val ,i_val,r_val,c_val: returns the logical,integer,real or
!>        character value
!> \param l_vals ,i_vals,r_vals,c_vals: returns the logical,integer,real or
!>        character arrays. The val reamins the owner of the array
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_val_get(section_vals,keyword_name,i_rep_section,&
       i_rep_val,n_rep_val,val,l_val,i_val,r_val,c_val,l_vals,i_vals,r_vals,&
       c_vals,ignore_required,explicit,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section, i_rep_val
    INTEGER, INTENT(out), OPTIONAL           :: n_rep_val
    TYPE(val_type), OPTIONAL, POINTER        :: val
    LOGICAL, INTENT(out), OPTIONAL           :: l_val
    INTEGER, INTENT(out), OPTIONAL           :: i_val
    REAL(KIND=DP), INTENT(out), OPTIONAL     :: r_val
    CHARACTER(LEN=*), INTENT(out), OPTIONAL  :: c_val
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals
    REAL(KIND=DP), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_vals
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: c_vals
    LOGICAL, INTENT(in), OPTIONAL            :: ignore_required
    LOGICAL, INTENT(out), OPTIONAL           :: explicit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_val_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, irs, len_key, &
                                                my_index, tmp_index
    LOGICAL                                  :: failure, my_ignore_required, &
                                                valRequested
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_sll_val_type), POINTER           :: vals
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: s_vals
    TYPE(val_type), POINTER                  :: my_val

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)

    my_ignore_required=.FALSE.
    IF (PRESENT(ignore_required)) my_ignore_required=ignore_required
    my_index=INDEX(keyword_name,'%')+1
    len_key=LEN_TRIM(keyword_name)
    IF (my_index>1) THEN
       DO
          tmp_index=INDEX(keyword_name(my_index:len_key),"%")
          IF (tmp_index<=0) EXIT
          my_index=my_index+tmp_index
       END DO
       s_vals => section_vals_get_subs_vals(section_vals,keyword_name(1:my_index-2),&
            error=error)
    ELSE
       s_vals=> section_vals
    END IF

    irk=1
    irs=1
    IF (PRESENT(i_rep_section)) irs=i_rep_section
    IF (PRESENT(i_rep_val))     irk=i_rep_val
    IF (PRESENT(val)) NULLIFY(val)
    IF (PRESENT(explicit)) explicit = .FALSE.
    section => s_vals%section
    valRequested=PRESENT(l_val).or.PRESENT(i_val).or.PRESENT(r_val).OR.&
         PRESENT(c_val).OR.PRESENT(l_vals).or.PRESENT(i_vals).OR.&
         PRESENT(r_vals).OR.PRESENT(c_vals)
    IF (.NOT. failure) THEN
       ik=section_get_keyword_index(s_vals%section,keyword_name(my_index:len_key),error=error)
       IF(ik==-2) &
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "section "//TRIM(section%name)//" does not contain keyword "//&
            TRIM(keyword_name(my_index:len_key)),error,failure)
       IF (.NOT.failure) THEN
          keyword => section%keywords(ik)%keyword
       END IF
       IF(.NOT.(irs>0.AND.irs<=SIZE(s_vals%subs_vals,2)))&
             CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
             routineP,"section repetition requested ("//cp_to_string(irs)//&
            ") out of bounds (1:"//cp_to_string(SIZE(s_vals%subs_vals,2))&
            //")",error,failure)
       NULLIFY(my_val)
    END IF
    IF (.NOT.failure) THEN
       IF (PRESENT(n_rep_val)) n_rep_val=0
       IF (irs<=SIZE(s_vals%values,2)) THEN ! the section was parsed
          vals => s_vals%values(ik,irs)%list
          IF (PRESENT(n_rep_val)) n_rep_val=cp_sll_val_get_length(vals,error=error)
          IF (.NOT.ASSOCIATED(vals)) THEN
             ! this keyword was not parsed
             IF ((.NOT.keyword%required.or.my_ignore_required)&
                  .AND. ASSOCIATED(keyword%default_value)) THEN
                my_val => keyword%default_value
                IF (PRESENT(n_rep_val)) n_rep_val=1
             END IF
          ELSE
             my_val => cp_sll_val_get_el_at(s_vals%values(ik,irs)%list,&
                  irk,error=error)
             IF (PRESENT(explicit)) explicit = .TRUE.
          END IF
       ELSE IF (.NOT.my_ignore_required.AND.section%required) THEN
          CALL cp_assert(.FALSE.,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "section "//TRIM(section%name)//&
               " is required, and needed but was not found",&
               error,failure)
       ELSE IF ((my_ignore_required.or..NOT.keyword%required).AND.&
            ASSOCIATED(keyword%default_value)) THEN
          IF (PRESENT(n_rep_val)) n_rep_val=1
          my_val => keyword%default_value
       END IF
       IF (PRESENT(val)) val => my_val
    END IF
    IF (.NOT.failure) THEN
       IF (valRequested) THEN
          IF(.NOT.ASSOCIATED(my_val))&
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Value requested, but no value set getting value from "//&
               "keyword "//TRIM(keyword_name(my_index:len_key))//" of section "//&
               TRIM(section%name),error,failure)
          CALL cp_error_init(sub_error,template_error=error,stop_level=cp_fatal_level)
          CALL val_get(my_val,l_val=l_val,i_val=i_val,r_val=r_val,&
               c_val=c_val,l_vals=l_vals,i_vals=i_vals,r_vals=r_vals,&
               c_vals=c_vals,error=sub_error)
          IF (cp_error_get_level(sub_error)/=0)&
               CALL cp_error_propagate_error(sub_error,routineP,&
               "getting value from  keyword "//TRIM(keyword_name(my_index:len_key))//" of section "//&
               TRIM(section%name),error=error,failure=failure)
          CALL cp_error_dealloc_ref(sub_error)
       END IF

    END IF
  END SUBROUTINE section_vals_val_get

! *****************************************************************************
!> \brief returns the requested list
!> \param section the root section
!> \param keyword_name the name of the keyword you want
!> \param isection which repetition of the section you are interested in
!>        (defaults to 1)
!> \note
!>      - most useful if the full list is needed anyway, so that faster iteration can be used
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE section_vals_list_get(section_vals,keyword_name,i_rep_section,&
       list,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    INTEGER, OPTIONAL                        :: i_rep_section
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_list_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irs, len_key, my_index, &
                                                tmp_index
    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: s_vals

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    NULLIFY(list)
    my_index=INDEX(keyword_name,'%')+1
    len_key=LEN_TRIM(keyword_name)
    IF (my_index>1) THEN
       DO
          tmp_index=INDEX(keyword_name(my_index:len_key),"%")
          IF (tmp_index<=0) EXIT
          my_index=my_index+tmp_index
       END DO
       s_vals => section_vals_get_subs_vals(section_vals,keyword_name(1:my_index-2),&
            error=error)
    ELSE
       s_vals=> section_vals
    END IF

    irs=1
    IF (PRESENT(i_rep_section)) irs=i_rep_section
    section => s_vals%section
    IF (.NOT. failure) THEN
       ik=section_get_keyword_index(s_vals%section,keyword_name(my_index:len_key),error=error)
       IF(ik==-2)&
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "section "//TRIM(section%name)//" does not contain keyword "//&
            TRIM(keyword_name(my_index:len_key)),error,failure)
       IF(.NOT.(irs>0.AND.irs<=SIZE(s_vals%subs_vals,2)))&
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "section repetition requested ("//cp_to_string(irs)//&
            ") out of bounds (1:"//cp_to_string(SIZE(s_vals%subs_vals,2))&
            //")",error,failure)
    END IF
    IF (.NOT. failure) THEN
       list => s_vals%values(ik,irs)%list
    ENDIF

  END SUBROUTINE section_vals_list_get

! *****************************************************************************
!> \brief sets the requested value
!> \param section the root section
!> \param keyword_name the name of the keyword you want (can be a path
!>        separated by '%')
!> \param isection which repetition of the section you are interested in
!>        (defaults to 1)
!> \param ival which repetition of the keyword/val you are interested in
!>        (defaults to 1)
!> \param l_val ,i_val,r_val,c_val: sets the logical,integer,real or
!>        character value
!> \param l_vals_ptr ,i_vals_ptr,r_vals,c_vals: sets the logical,integer,real or
!>        character arrays. The val becomes the owner of the array
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_val_set(section_vals,keyword_name,i_rep_section,i_rep_val,&
       val,l_val,i_val,r_val,c_val,l_vals_ptr,i_vals_ptr,r_vals_ptr,c_vals_ptr,&
       ignore_required,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section, i_rep_val
    TYPE(val_type), OPTIONAL, POINTER        :: val
    LOGICAL, INTENT(in), OPTIONAL            :: l_val
    INTEGER, INTENT(in), OPTIONAL            :: i_val
    REAL(KIND=DP), INTENT(in), OPTIONAL      :: r_val
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: c_val
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals_ptr
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals_ptr
    REAL(KIND=DP), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_vals_ptr
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: c_vals_ptr
    LOGICAL, INTENT(in), OPTIONAL            :: ignore_required
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_val_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, irs, len_key, &
                                                my_index, tmp_index
    LOGICAL                                  :: failure, my_ignore_required, &
                                                valSet
    TYPE(cp_error_type)                      :: sub_error
    TYPE(cp_sll_val_type), POINTER           :: vals
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: s_vals
    TYPE(val_type), POINTER                  :: my_val, old_val

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)

    my_ignore_required=.FALSE.
    IF (PRESENT(ignore_required)) my_ignore_required=ignore_required
    my_index=INDEX(keyword_name,'%')+1
    len_key=LEN_TRIM(keyword_name)
    IF (my_index>1) THEN
       DO
          tmp_index=INDEX(keyword_name(my_index:len_key),"%")
          IF (tmp_index<=0) EXIT
          my_index=my_index+tmp_index
       END DO
       s_vals => section_vals_get_subs_vals(section_vals,keyword_name(1:my_index-2),&
            error=error)
    ELSE
       s_vals=> section_vals
    END IF

    irk=1
    irs=1
    IF (PRESENT(i_rep_section)) irs=i_rep_section
    IF (PRESENT(i_rep_val))     irk=i_rep_val
    section => s_vals%section
    IF (.NOT. failure) THEN
       ik=section_get_keyword_index(s_vals%section,keyword_name(my_index:len_key),error=error)
       IF(ik==-2)&
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "section "//TRIM(section%name)//" does not contain keyword "//&
            TRIM(keyword_name(my_index:len_key)),error,failure)
       ! Add values..
       DO
          IF (irs<=SIZE(s_vals%values,2)) EXIT
          CALL section_vals_add_values(s_vals,error=error)
       END DO
       IF(.NOT.(irs>0.AND.irs<=SIZE(s_vals%subs_vals,2)))&
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "section repetition requested ("//cp_to_string(irs)//&
            ") out of bounds (1:"//cp_to_string(SIZE(s_vals%subs_vals,2))&
            //")",error,failure)
    END IF
    IF (.NOT.failure) THEN
       keyword => s_vals%section%keywords(ik)%keyword
       NULLIFY(my_val)
       IF (PRESENT(val)) my_val => val
       valSet=PRESENT(l_val).OR.PRESENT(i_val).OR.PRESENT(r_val).OR.&
            PRESENT(c_val).OR.PRESENT(l_vals_ptr).OR.PRESENT(i_vals_ptr).OR.&
            PRESENT(r_vals_ptr).OR.PRESENT(c_vals_ptr)
       IF (ASSOCIATED(my_val)) THEN
          ! check better?
          IF(valSet)&
               CALL cp_assert(.FALSE.,&
               cp_failure_level,cp_assertion_failed,routineP,&
               " both val and values present, in setting "//&
               "keyword "//TRIM(keyword_name(my_index:len_key))//" of section "//&
               TRIM(section%name),error,failure)
       ELSE
          ! ignore ?
          IF(.NOT.valSet)&
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               " empty value in setting "//&
               "keyword "//TRIM(keyword_name(my_index:len_key))//" of section "//&
               TRIM(section%name),error,failure)
          CPPrecondition(valSet,cp_failure_level,routineP,error,failure)
          CALL cp_error_init(sub_error,template_error=error,&
               stop_level=cp_fatal_level)
          IF (keyword%type_of_var==lchar_t) THEN
             CALL val_create(my_val, lc_val=c_val, lc_vals_ptr=c_vals_ptr,&
                  error=sub_error)
          ELSE
             CALL val_create(my_val,l_val=l_val,i_val=i_val,r_val=r_val,&
                  c_val=c_val,l_vals_ptr=l_vals_ptr,i_vals_ptr=i_vals_ptr,&
                  r_vals_ptr=r_vals_ptr,&
                  c_vals_ptr=c_vals_ptr,enum=keyword%enum,error=sub_error)
          END IF
          CPPostcondition(ASSOCIATED(my_val),cp_failure_level,routineP,sub_error,failure)
          CPPostcondition(my_val%type_of_var==keyword%type_of_var,cp_failure_level,routineP,error,failure)
          IF(cp_error_get_level(sub_error)/=0)&
               CALL cp_error_propagate_error(sub_error,routineP,&
               "setting value from  keyword "//TRIM(keyword_name(my_index:len_key))//&
               " of section "// TRIM(section%name),error=error,failure=failure)
          CALL cp_error_dealloc_ref(sub_error)
       END IF
    END IF
    IF (.NOT.failure) THEN
       vals => s_vals%values(ik,irs)%list
       IF (irk==-1) THEN
          CALL cp_sll_val_insert_el_at(vals,my_val,index=-1,error=error)
       ELSE IF (irk <= cp_sll_val_get_length(vals,error)) THEN
          IF(irk<=0)&
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"invalid irk "//TRIM(ADJUSTL(cp_to_string(irk)))//&
               " in keyword "//TRIM(keyword_name(my_index:len_key))//" of section "//&
               TRIM(section%name),error,failure)
          old_val => cp_sll_val_get_el_at(vals,index=irk,error=error)
          CALL val_release(old_val,error=error)
          CALL cp_sll_val_set_el_at(vals,value=my_val,index=irk,error=error)
       ELSE IF (irk>cp_sll_val_get_length(vals,error)+1) THEN
          ! change?
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"cannot add extra keyword repetitions to keyword"&
               //TRIM(keyword_name(my_index:len_key))//" of section "//&
               TRIM(section%name),error,failure)
       ELSE
          CALL cp_sll_val_insert_el_at(vals,my_val,index=irk,error=error)
       END IF
       s_vals%values(ik,irs)%list => vals
       NULLIFY(my_val)
    END IF
  END SUBROUTINE section_vals_val_set

! *****************************************************************************
!> \brief unsets (removes) the requested value (if it is a keyword repetitions
!>      removes the repetition, so be careful: the repetition indices bigger
!>      than the actual change.
!> \param section the root section
!> \param keyword_name the name of the keyword you want (can be a path
!>        separated by '%')
!> \param isection which repetition of the section you are interested in
!>        (defaults to 1)
!> \param ival which repetition of the keyword/val you are interested in
!>        (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_val_unset(section_vals,keyword_name,i_rep_section,&
       i_rep_val,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: keyword_name
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section, i_rep_val
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_val_unset', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ik, irk, irs, len_key, &
                                                my_index, tmp_index
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: pos
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: s_vals
    TYPE(val_type), POINTER                  :: old_val

    failure=.FALSE.
    NULLIFY(pos)
    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)

    my_index=INDEX(keyword_name,'%')+1
    len_key=LEN_TRIM(keyword_name)
    IF (my_index>1) THEN
       DO
          tmp_index=INDEX(keyword_name(my_index:len_key),"%")
          IF (tmp_index<=0) EXIT
          my_index=my_index+tmp_index
       END DO
       s_vals => section_vals_get_subs_vals(section_vals,keyword_name(1:my_index-2),&
            error=error)
    ELSE
       s_vals=> section_vals
    END IF

    irk=1
    irs=1
    IF (PRESENT(i_rep_section)) irs=i_rep_section
    IF (PRESENT(i_rep_val)) irk=i_rep_val
    section => s_vals%section
    IF (.NOT. failure) THEN
       ik=section_get_keyword_index(s_vals%section,keyword_name(my_index:len_key),error=error)
       IF(ik==-2)&
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "section "//TRIM(section%name)//" does not contain keyword "//&
            TRIM(keyword_name(my_index:len_key)),error,failure)
    END IF
    IF (.NOT.failure) THEN
       ! ignore unset of non set values
       IF (irs<=SIZE(s_vals%values,2)) THEN
          IF(.NOT.(irs>0.AND.irs<=SIZE(s_vals%subs_vals,2)))&
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "section repetition requested ("//cp_to_string(irs)//&
               ") out of bounds (1:"//cp_to_string(SIZE(s_vals%subs_vals,2))&
               //")",error,failure)
          IF (irk==-1) THEN
             pos => cp_sll_val_get_rest(s_vals%values(ik,irs)%list, iter=-1,error=error)
          ELSE
             pos => cp_sll_val_get_rest(s_vals%values(ik,irs)%list, iter=irk-1,error=error)
          END IF
          IF (ASSOCIATED(pos)) THEN
             old_val => cp_sll_val_get_el_at(s_vals%values(ik,irs)%list,index=irk,&
                  error=error)
             CALL val_release(old_val,error=error)
             CALL cp_sll_val_rm_el_at(s_vals%values(ik,irs)%list,index=irk,&
                  error=error)
          END IF
       END IF
    END IF

  END SUBROUTINE section_vals_val_unset

! *****************************************************************************
!> \brief Checks if all sections or keywords defined have been supported for
!>      release version
!> \param section_vals the section to check for released features
!> \param unit_nr the unit where to write to
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      skips required sections which weren't read
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE section_vals_check_release(section_vals,enable_unsupported_features,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    LOGICAL, INTENT(IN)                      :: enable_unsupported_features
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_check_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_rep_s, ik, isec, nr, nval, &
                                                output_unit
    LOGICAL                                  :: check_key, defaultSection, &
                                                failure
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: sval

    failure=.FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL cp_error_init(suberror,template_error=error)
    IF (.NOT. failure) THEN
       CALL section_vals_get(section_vals, n_repetition=nr, section=section, error=error)
       DO i_rep_s=1,nr
          defaultSection=(SIZE(section_vals%values,2)==0)
          IF (.NOT.(section%required.AND.defaultSection)) THEN
             IF (.NOT.section%supported_feature) THEN
                IF (output_unit>0) THEN
                   WRITE(output_unit,'(/)')
                   WRITE(output_unit,'(T2,"RELEASE_ERROR|",1X,A)')&
                        "Section <"//TRIM(section%name)//"> is not supported for CP2K release version."
                   IF (enable_unsupported_features) THEN
                      WRITE(output_unit,'(T2,"RELEASE_ERROR|",1X,A)')&
                           "Calculation continues but the results could be totally meaningless!"
                   ELSE
                      WRITE(output_unit,'(T2,"RELEASE_ERROR|",1X,A)')&
                           "You may consider to enable the keyword: GLOBAL%ENABLE_UNSUPPORTED_FEATURES or"
                      WRITE(output_unit,'(T2,"RELEASE_ERROR|",1X,A)')"download a Development Version!"
                   END IF
                END IF
                IF (.NOT.enable_unsupported_features)  failure = .TRUE.
             END IF
             DO ik=-1,section%n_keywords
                keyword => section%keywords(ik)%keyword
                IF (ASSOCIATED(keyword)) THEN
                   CALL section_vals_val_get(section_vals,keyword%names(1),i_rep_s,n_rep_val=nval,error=error)
                   check_key = .FALSE.
                   IF (i_rep_s<=SIZE(section_vals%values,2)) THEN
                      vals  => section_vals%values(ik,i_rep_s)%list
                      new_pos => vals
                      IF (ASSOCIATED(new_pos)) THEN
                         ! Keyword was parsed
                         check_key = .TRUE.
                      END IF
                   ELSEIF (.NOT.keyword%required.AND.ASSOCIATED(keyword%default_value)) THEN
                      ! Keyword was not parsed but a default value is present
                      check_key = .TRUE.
                   END IF
                   IF ((check_key).AND.(.NOT.keyword%supported_feature)) THEN
                      IF (output_unit>0) THEN
                         WRITE(output_unit,'(T2,"RELEASE_ERROR|",2X,A)')&
                              "Keyword <"//TRIM(keyword%names(1))//"> is not supported for CP2K release version."
                         IF (enable_unsupported_features) THEN
                            WRITE(output_unit,'(T2,"RELEASE_ERROR|",2X,A)')&
                                 "Calculation continues but the results could be totally meaningless!"
                         ELSE
                            WRITE(output_unit,'(T2,"RELEASE_ERROR|",2X,A)')&
                                 "You may consider to enable the keyword: GLOBAL%ENABLE_UNSUPPORTED_FEATURES or"
                            WRITE(output_unit,'(T2,"RELEASE_ERROR|",2X,A)')"download a Development Version!"
                         END IF
                      END IF
                      IF (.NOT.enable_unsupported_features)  failure = .TRUE.
                   END IF
                END IF
             END DO
             IF (.NOT.failure) THEN
                IF (ASSOCIATED(section_vals%subs_vals)) THEN
                   DO isec=1,SIZE(section_vals%subs_vals,1)
                      sval => section_vals%subs_vals(isec,i_rep_s)%section_vals
                      IF (ASSOCIATED(sval)) THEN
                         ! Section was parsed
                         CALL section_vals_check_release(sval,enable_unsupported_features,error)
                      END IF
                   END DO
                END IF
             END IF
             CPPrecondition(.NOT.failure,cp_failure_level,routineP,error,failure)
          END IF
       END DO
    END IF
    CALL cp_error_dealloc_ref(suberror)

  END SUBROUTINE section_vals_check_release

! *****************************************************************************
!> \brief writes the values in the given section in a way that is suitable to
!>      the automatic parsing
!> \param section_vals the section to write out
!> \param unit_nr the unit where to write to
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      skips required sections which weren't read
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_vals_write(section_vals,unit_nr,hide_root,hide_defaults,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: hide_root, hide_defaults
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: myfmt
    INTEGER                                  :: i_rep_s, ik, isec, ival, nr, &
                                                nval
    INTEGER, SAVE                            :: indent = 1
    LOGICAL                                  :: defaultSection, explicit, &
                                                failure, my_hide_defaults, &
                                                my_hide_root
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_sll_val_type), POINTER           :: new_pos, vals
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section
    TYPE(section_vals_type), POINTER         :: sval
    TYPE(val_type), POINTER                  :: val

    failure=.FALSE.
    my_hide_root=.FALSE.
    my_hide_defaults=.TRUE.
    IF (PRESENT(hide_root)) my_hide_root=hide_root
    IF (PRESENT(hide_defaults)) my_hide_defaults=hide_defaults

    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    CALL cp_error_init(suberror,template_error=error)
    IF ((.NOT.failure).AND.(unit_nr>0)) THEN
       CALL section_vals_get(section_vals, explicit=explicit, n_repetition=nr, section=section,&
            error=error)
       IF (explicit.OR.(.NOT.my_hide_defaults)) THEN
          DO i_rep_s=1,nr
             IF (.NOT.my_hide_root) THEN
                WRITE(myfmt,*)indent,"X"
                CALL compress(myfmt,full=.TRUE.)
                IF (ASSOCIATED(section%keywords(-1)%keyword)) THEN
                   WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//",'&',a,' ')",advance="NO") TRIM(section%name)
                ELSE
                   WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//",'&',a)") TRIM(section%name)
                END IF
             END IF
             defaultSection=(SIZE(section_vals%values,2)==0)
             IF (.NOT.(section%required.AND.defaultSection)) THEN
                IF (.NOT.my_hide_root) indent = indent + 2
                WRITE(myfmt,*)indent,"X"
                CALL compress(myfmt,full=.TRUE.)
                DO ik=-1,section%n_keywords
                   keyword => section%keywords(ik)%keyword
                   IF (ASSOCIATED(keyword)) THEN
                      IF (keyword%type_of_var/=no_t.AND.keyword%names(1)(1:2)/="__") THEN
                         CALL section_vals_val_get(section_vals,keyword%names(1),&
                              i_rep_s,n_rep_val=nval,error=error)
                         IF (i_rep_s<=SIZE(section_vals%values,2)) THEN
                            ! Section was parsed
                            vals  => section_vals%values(ik,i_rep_s)%list
                            DO ival=1,nval
                               IF (ival==1) THEN
                                  new_pos => vals
                               ELSE
                                  new_pos => new_pos%rest
                               END IF
                               IF (.NOT.ASSOCIATED(new_pos)) THEN
                                  ! this keyword was not parsed
                                  IF (.NOT.keyword%required.AND. ASSOCIATED(keyword%default_value)) THEN
                                     val => keyword%default_value
                                     IF (my_hide_defaults) CYCLE
                                  END IF
                               ELSE
                                  val => new_pos%first_el
                               END IF
                               IF (cp_error_get_level(suberror)/=0)&
                                    CALL cp_error_propagate_error(suberror, fromWhere=routineP,&
                                    message="for val "// TRIM(keyword%names(1))//" in section "//&
                                    TRIM(section%name),error=error, failure=failure)
                               IF (keyword%names(1)/='_DEFAULT_KEYWORD_'.AND.&
                                    keyword%names(1)/='_SECTION_PARAMETERS_') THEN
                                  WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//",a,' ')",advance="NO")&
                                       TRIM(keyword%names(1))
                               ELSEIF (keyword%names(1) =='_DEFAULT_KEYWORD_' .AND.&
                                    keyword%type_of_var/=lchar_t) THEN
                                  WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//")",advance="NO")
                               END IF
                               CALL val_write(val,unit_nr=unit_nr,unit=keyword%unit,&
                                    fmt=myfmt,error=suberror)
                               IF (cp_error_get_level(suberror)/=0)&
                                    CALL cp_error_propagate_error(suberror, fromWhere=routineP,&
                                    message="for val "// TRIM(keyword%names(1))//" in section "//&
                                    TRIM(section%name),error=error, failure=failure)
                               IF (failure) EXIT
                            END DO
                         ELSEIF (.NOT.keyword%required.AND. ASSOCIATED(keyword%default_value)) THEN
                            ! Section was not parsed but default for the keywords may exist
                            IF (my_hide_defaults) CYCLE
                            val => keyword%default_value
                            IF (cp_error_get_level(suberror)/=0)&
                                 CALL cp_error_propagate_error(suberror, fromWhere=routineP,&
                                 message="for val "// TRIM(keyword%names(1))//" in section "//&
                                 TRIM(section%name),error=error, failure=failure)
                            IF (keyword%names(1)/='_DEFAULT_KEYWORD_'.AND.&
                                 keyword%names(1)/='_SECTION_PARAMETERS_') THEN
                               WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//",a,' ')",advance="NO")&
                                    TRIM(keyword%names(1))
                            ELSEIF (keyword%names(1) =='_DEFAULT_KEYWORD_' .AND.&
                                 keyword%type_of_var/=lchar_t) THEN
                               WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//")",advance="NO")
                            END IF
                            CALL val_write(val,unit_nr=unit_nr,unit=keyword%unit,&
                                 fmt=myfmt,error=suberror)
                            IF (cp_error_get_level(suberror)/=0)&
                                 CALL cp_error_propagate_error(suberror, fromWhere=routineP,&
                                 message="for val "// TRIM(keyword%names(1))//" in section "//&
                                 TRIM(section%name),error=error, failure=failure)
                            IF (failure) EXIT
                         END IF
                      END IF
                   END IF
                   IF (failure) EXIT
                END DO
                IF (.NOT.failure) THEN
                   IF (ASSOCIATED(section_vals%subs_vals)) THEN
                      DO isec=1,SIZE(section_vals%subs_vals,1)
                         sval => section_vals%subs_vals(isec,i_rep_s)%section_vals
                         IF (ASSOCIATED(sval)) THEN
                            CALL section_vals_write(sval,unit_nr=unit_nr,hide_defaults=hide_defaults,error=error)
                         END IF
                      END DO
                   END IF
                END IF
             END IF
             IF (.NOT.my_hide_root) THEN
                indent = indent - 2
                WRITE (UNIT=unit_nr,FMT="(A)")&
                     REPEAT(" ",indent)//"&END "//TRIM(section%name)
             END IF
          END DO
       END IF
    END IF
    CALL cp_error_dealloc_ref(suberror)

  END SUBROUTINE section_vals_write

! *****************************************************************************
!> \brief writes the values in the given section in xml
! *****************************************************************************
  RECURSIVE SUBROUTINE write_section_xml(section,level,unit_number,error)

    TYPE(section_type), POINTER              :: section
    INTEGER, INTENT(IN)                      :: level, unit_number
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_section_xml', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: repeats, required
    CHARACTER(LEN=8)                         :: short_string
    INTEGER                                  :: i, l0, l1, l2
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(section)) THEN

       CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)

       ! Indentation for current level, next level, etc.

       l0 = level
       l1 = level + 1
       l2 = level + 2

       IF (section%required) THEN
          required = "yes"
       ELSE
          required = "no "
       END IF

       IF (section%repeats) THEN
          repeats = "yes"
       ELSE
          repeats = "no "
       END IF

       WRITE (UNIT=unit_number,FMT="(A)")&
            REPEAT(" ",l0)//"<SECTION required="""//TRIM(required)//&
            """ repeats="""//TRIM(repeats)//""">",&
            REPEAT(" ",l1)//"<NAME>"//TRIM(section%name)//"</NAME>",&
            REPEAT(" ",l1)//"<DESCRIPTION>"//&
            TRIM(substitute_special_xml_tokens(a2s(section%description)))&
            //"</DESCRIPTION>"

       IF (ASSOCIATED(section%citations)) THEN
          DO i=1,SIZE(section%citations,1)
             short_string = ""
             WRITE (UNIT=short_string,FMT="(I8)") section%citations(i)
             WRITE (UNIT=unit_number,FMT="(A)")&
                  REPEAT(" ",l1)//"<REFERENCE>",&
                  REPEAT(" ",l2)//"<NAME>"//TRIM(get_citation_key(section%citations(i)))//"</NAME>",&
                  REPEAT(" ",l2)//"<NUMBER>"//TRIM(ADJUSTL(short_string))//"</NUMBER>",&
                  REPEAT(" ",l1)//"</REFERENCE>"
          END DO
       END IF

       DO i=-1,section%n_keywords
          IF (ASSOCIATED(section%keywords(i)%keyword)) THEN
             CALL write_keyword_xml(section%keywords(i)%keyword,l1,unit_number,&
                  error)
          END IF
       END DO

       DO i=1,section%n_subsections
          CALL write_section_xml(section%subsections(i)%section,l1,unit_number,&
               error)
       END DO

       WRITE (UNIT=unit_number,FMT="(A)") REPEAT(" ",l0)//"</SECTION>"

    END IF

  END SUBROUTINE write_section_xml

! *****************************************************************************
  RECURSIVE SUBROUTINE section_typo_match(section,unknown_string,location_string,&
       matching_rank,matching_string,error)

    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*)                         :: unknown_string, &
                                                location_string
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: matching_rank
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(INOUT)                          :: matching_string
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'section_typo_match', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=LEN(matching_string(1)))   :: line
    INTEGER                                  :: i, imatch, imax, irank
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(section)) THEN
       CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
       imatch=typo_match(TRIM(section%name),TRIM(unknown_string))
       IF (imatch>0) THEN
          WRITE(line,'(T2,A)') " subsection "//TRIM(section%name)//" in section "//TRIM(location_string)
          imax=SIZE(matching_rank,1)
          irank=imax+1
          DO I=imax,1,-1
             IF (imatch>matching_rank(I)) irank=i
          ENDDO
          IF (irank<=imax) THEN
             matching_rank(irank+1:imax)=matching_rank(irank:imax-1)
             matching_string(irank+1:imax)=matching_string(irank:imax-1)
             matching_rank(irank)=imatch
             matching_string(irank)=line
          ENDIF
       END IF

       DO i=-1,section%n_keywords
          IF (ASSOCIATED(section%keywords(i)%keyword)) THEN
             CALL keyword_typo_match(section%keywords(i)%keyword,unknown_string,location_string// &
                  "%"//TRIM(section%name),matching_rank,matching_string,error)
          END IF
       END DO

       DO i=1,section%n_subsections
          CALL section_typo_match(section%subsections(i)%section,unknown_string,&
               location_string//"%"//TRIM(section%name),matching_rank,matching_string,error)
       END DO

    END IF

  END SUBROUTINE section_typo_match

! *****************************************************************************
!> \brief replaces of the requested subsection with the one given
!> \param section_vals the root section
!> \param subsection_name the name of the subsection to replace
!> \param new_section_vals the new section_vals to use
!> \param i_rep_section index of the repetition of section_vals of which
!>        you want to replace the subsection (defaults to 1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_set_subs_vals(section_vals,subsection_name,&
       new_section_vals,i_rep_section,error)
    TYPE(section_vals_type), POINTER         :: section_vals
    CHARACTER(len=*), INTENT(in)             :: subsection_name
    TYPE(section_vals_type), POINTER         :: new_section_vals
    INTEGER, INTENT(in), OPTIONAL            :: i_rep_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_set_subs_vals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: irep, isection, len_key, &
                                                my_index, tmp_index
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: s_vals

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(new_section_vals),cp_failure_level,routineP,error,failure)
    CPPrecondition(new_section_vals%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       irep=1
       IF (PRESENT(i_rep_section)) irep=i_rep_section

       my_index=INDEX(subsection_name,'%')+1
       len_key=LEN_TRIM(subsection_name)
       IF (my_index>1) THEN
          DO
             tmp_index=INDEX(subsection_name(my_index:len_key),"%")
             IF (tmp_index<=0) EXIT
             my_index=my_index+tmp_index
          END DO
          s_vals => section_vals_get_subs_vals(section_vals,subsection_name(1:my_index-2),&
               error=error)
       ELSE
          s_vals=> section_vals
       END IF

       CPPrecondition(irep<=SIZE(s_vals%subs_vals,2),cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT.failure) THEN
       isection=section_get_subsection_index(s_vals%section,subsection_name(my_index:LEN_TRIM(subsection_name)),&
            error=error)
       IF(isection<=0)&
            CALL cp_assert(.FALSE.,cp_failure_level,&
            cp_assertion_failed,routineP,&
            "could not find subsection "//subsection_name(my_index:LEN_TRIM(subsection_name))//" in section "//&
            TRIM(section_vals%section%name)//" at "//&
            CPSourceFileRef,&
            error,failure)
    END IF
    IF (.NOT.failure) THEN
       CALL section_vals_retain(new_section_vals,error=error)
       CALL section_vals_release(s_vals%subs_vals(isection,irep)%section_vals,error=error)
       s_vals%subs_vals(isection,irep)%section_vals => new_section_vals
    END IF

  END SUBROUTINE section_vals_set_subs_vals

! *****************************************************************************
!> \brief creates a deep copy from section_vals_in to section_vals_out
!> \param section_vals_in the section_vals to copy
!> \param section_vals_out the section_vals to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE section_vals_duplicate(section_vals_in,section_vals_out,&
       i_rep_start, i_rep_end, error)
    TYPE(section_vals_type), POINTER         :: section_vals_in, &
                                                section_vals_out
    INTEGER, INTENT(IN), OPTIONAL            :: i_rep_start, i_rep_end
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_duplicate', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(section_vals_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(section_vals_out),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL section_vals_create(section_vals_out,section_vals_in%section,&
            error=error)
       CALL section_vals_copy(section_vals_in,section_vals_out,i_rep_start,i_rep_end,error=error)
    END IF
  END SUBROUTINE section_vals_duplicate

! *****************************************************************************
!> \brief deep copy from section_vals_in to section_vals_out
!> \param section_vals_in the section_vals to copy
!> \param section_vals_out the section_vals where to copy
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      private, only works with a newly initialized section_vals_out
!> \author fawzi
! *****************************************************************************
  RECURSIVE SUBROUTINE section_vals_copy(section_vals_in,section_vals_out,&
       i_rep_low,i_rep_high,error)
    TYPE(section_vals_type), POINTER         :: section_vals_in, &
                                                section_vals_out
    INTEGER, INTENT(IN), OPTIONAL            :: i_rep_low, i_rep_high
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'section_vals_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iend, irep, isec, istart, ival
    LOGICAL                                  :: failure
    TYPE(cp_sll_val_type), POINTER           :: v1, v2
    TYPE(val_type), POINTER                  :: el

    failure=.FALSE.
    NULLIFY(v2,el)

    CPPrecondition(ASSOCIATED(section_vals_in),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(section_vals_out),cp_failure_level,routineP,error,failure)
    !  IF(section_vals_in%section%id_nr/=section_vals_out%section%id_nr)&
    !     CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,routineP,&
    !          CPSourceFileRef,&
    !          error,failure)
    IF (.NOT.failure) THEN
       istart = 1
       iend   = SIZE(section_vals_in%values,2)
       IF (PRESENT(i_rep_low)) istart=i_rep_low
       IF (PRESENT(i_rep_high)) iend=i_rep_high
       DO irep=istart,iend
          CALL section_vals_add_values(section_vals_out,error=error)
          DO ival=LBOUND(section_vals_in%values,1),UBOUND(section_vals_in%values,1)
             v1=>section_vals_in%values(ival,irep)%list
             IF (ASSOCIATED(v1)) THEN
                CALL val_duplicate(v1%first_el,el,error=error)
                CALL cp_sll_val_create(v2,el,error=error)
                NULLIFY(el)
                section_vals_out%values(ival,irep-istart+1)%list => v2
                DO
                   IF (.not.ASSOCIATED(v1%rest)) EXIT
                   v1 => v1%rest
                   CALL val_duplicate(v1%first_el,el,error=error)
                   CALL cp_sll_val_create(v2%rest,first_el=el,error=error)
                   NULLIFY(el)
                   v2 => v2%rest
                END DO
             END IF
          END DO
       END DO
       IF (.NOT.PRESENT(i_rep_low).AND.(.NOT.PRESENT(i_rep_high))) THEN
          CALL cp_assert(SIZE(section_vals_in%values,2)==SIZE(section_vals_out%values,2),&
               cp_failure_level,cp_assertion_failed,routineP,&
               CPSourceFileRef,&
               error,failure)
          CALL cp_assert(SIZE(section_vals_in%subs_vals,2)==SIZE(section_vals_out%subs_vals,2),&
               cp_failure_level,cp_assertion_failed,routineP,&
               CPSourceFileRef,&
               error,failure)
       END IF
       iend   = SIZE(section_vals_in%subs_vals,2)
       IF (PRESENT(i_rep_high)) iend=i_rep_high
       DO irep=istart,iend
          DO isec=1,SIZE(section_vals_in%subs_vals,1)
             CALL section_vals_copy(section_vals_in%subs_vals(isec,irep)%section_vals,&
                  section_vals_out%subs_vals(isec,irep-istart+1)%section_vals,error=error)
          END DO
       END DO
    END IF
  END SUBROUTINE section_vals_copy

END MODULE input_section_types
