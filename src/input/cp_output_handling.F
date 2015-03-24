!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines to handle the output, The idea is to remove the
!>      decision of wheter to output and what to output from the code
!>      that does the output, and centralize it here.
!> \note
!>      These were originally together with the log handling routines,
!>      but have been spawned off. Some dependencies are still there,
!>      and some of the comments about log handling also applies to output
!>      handling: @see cp_log_handling
!> \par History
!>      12.2001 created [fawzi]
!>      08.2002 updated to new logger [fawzi]
!>      10.2004 big rewrite of the output methods, connected to the new
!>              input, and iteration_info [fawzi]
!>      08.2005 property flags [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_output_handling
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_iter_types,                   ONLY: cp_iteration_info_release,&
                                             cp_iteration_info_retain,&
                                             cp_iteration_info_type,&
                                             each_desc_labels,&
                                             each_possible_labels
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length
  USE machine,                         ONLY: m_mov
  USE memory_utilities,                ONLY: reallocate
  USE string_utilities,                ONLY: compress,&
                                             s2a
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_output_handling'
  PUBLIC :: cp_print_key_should_output, cp_iterate, cp_add_iter_level, cp_rm_iter_level
  PUBLIC :: cp_iter_string, cp_print_key_section_create
  PUBLIC :: cp_print_key_unit_nr, cp_print_key_finished_output
  PUBLIC :: cp_print_key_generate_filename, cp_printkey_is_on

  INTEGER, PARAMETER, PUBLIC               :: add_last_no=0,&
                                              add_last_numeric=1,&
                                              add_last_symbolic=2
  INTEGER, PARAMETER, PUBLIC               :: silent_print_level=0,&
                                              low_print_level=1,&
                                              medium_print_level=2,&
                                              high_print_level=3,&
                                              debug_print_level=4

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

! *****************************************************************************
!> \brief stores the flags_env controlling the output of properties
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param id_nr identification number (unique to each istance)
!> \param n_flags number of flags stored in this type
!> \param names names of the stored flags
!> \param control_val value of the flag
!> \param input the input (with all the printkeys)
!> \param logger logger and iteration information (to know if output is needed)
!> \param strict if flags that were not stored can be read
!> \param default_val default value of the flags that are not explicitly
!>        stored
!> \note
!>      Two features of this object should be:
!>        1) easy state storage, one should be able to store the state of the
!>           flags, to some changes to them just for one (or few) force evaluations
!>           and then reset the original state. The actual implementation is good
!>           in this respect
!>        2) work well with subsections. This is a problem at the moment, as
!>           if you pass just a subsection of the input the control flags get lost.
!>        A better implementation should be done storing the flags also in the
!>        input itself to be transparent
!> \author fawzi
! *****************************************************************************
  TYPE cp_out_flags_type
     INTEGER :: ref_count, id_nr, n_flags
     CHARACTER(default_string_length), DIMENSION(:), POINTER :: names
     INTEGER, DIMENSION(:), POINTER :: control_val
     TYPE(section_vals_type), POINTER :: input
     TYPE(cp_logger_type), POINTER :: logger
     LOGICAL :: strict
     INTEGER :: default_val
  END TYPE cp_out_flags_type

CONTAINS

! *****************************************************************************
!> \brief creates a print_key section
!> \param print_key_section the print key to create
!> \param name the name of the print key
!> \param description the description of the print key
!> \param print_level print level starting at which the printing takes place
!>        (defaults to debug_print_level)
!> \param each_iter_names ...
!> \param each_iter_values ...
!> \param add_last ...
!> \param filename ...
!> \param common_iter_levels ...
!> \param citations ...
!> \param supported_feature ...
!> \param unit_str specifies an unit of measure for output quantity. If not
!>        provided the control is totally left to how the output was coded
!>        (i.e. USERS have no possibility to change it)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_print_key_section_create(print_key_section, name, description, &
       print_level,each_iter_names,each_iter_values,add_last,filename,&
       common_iter_levels,citations,supported_feature,unit_str,error)
    TYPE(section_type), POINTER              :: print_key_section
    CHARACTER(len=*), INTENT(IN)             :: name, description
    INTEGER, INTENT(IN), OPTIONAL            :: print_level
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: each_iter_names
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: each_iter_values
    INTEGER, INTENT(IN), OPTIONAL            :: add_last
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: filename
    INTEGER, INTENT(IN), OPTIONAL            :: common_iter_levels
    INTEGER, DIMENSION(:), OPTIONAL          :: citations
    LOGICAL, INTENT(IN), OPTIONAL            :: supported_feature
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: unit_str
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_section_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: my_filename
    INTEGER                                  :: i_each, i_iter, my_add_last, &
                                                my_comm_iter_levels, &
                                                my_print_level, my_value
    LOGICAL                                  :: check, ext_each, failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(print_key_section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       my_print_level=debug_print_level
       IF (PRESENT(print_level)) my_print_level=print_level

       CALL section_create(print_key_section,name=name,description=description,&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=citations, supported_feature=supported_feature, error=error)

       NULLIFY(keyword, subsection)
       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="Level starting at which this property is printed",&
            usage="silent",&
            default_i_val=my_print_level,lone_keyword_i_val=silent_print_level,&
            enum_c_vals=s2a("on","off","silent","low","medium","high","debug"),&
            enum_i_vals=(/ silent_print_level-1,debug_print_level+1,&
            silent_print_level, low_print_level,&
            medium_print_level,high_print_level,debug_print_level/),&
            supported_feature=supported_feature,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="__CONTROL_VAL",&
            description=' hidden parameter that controls storage, printing,...'//&
            ' of the print_key',&
            default_i_val=cp_out_default, supported_feature=supported_feature,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="EACH",&
            description="This section specifies how often this proprety is printed."//&
            "Each keyword inside this section is mapping to a specific iteration level and "//&
            "the value of each of these keywords is matched with the iteration level during "//&
            "the calculation. How to handle the last iteration is treated "//&
            "separately in ADD_LAST (this mean that each iteration level (MD, GEO_OPT, etc..), "//&
            "though equal to 0, might print the last iteration). If an iteration level is specified "//&
            "that is not present in the flow of the calculation it is just ignored.",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            citations=citations, supported_feature=supported_feature, error=error)

       ! Enforce the presence or absence of both.. or give an error
       check = (PRESENT(each_iter_names)).EQV.(PRESENT(each_iter_values))
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       ext_each = (PRESENT(each_iter_names)).AND.(PRESENT(each_iter_values))

       DO i_each = 1, SIZE(each_possible_labels)
          my_value = 1
          IF (ext_each) THEN
             check    = SUM(INDEX(each_iter_names,each_possible_labels(i_each)))<=1
             CPPrecondition(check,cp_failure_level,routineP,error,failure)
             DO i_iter = 1, SIZE(each_iter_names)
                IF (INDEX(TRIM(each_iter_names(i_iter)),TRIM(each_possible_labels(i_each)))/=0) THEN
                   my_value = each_iter_values(i_iter)
                END IF
             END DO
          END IF
          CALL keyword_create(keyword, name=TRIM(each_possible_labels(i_each)),&
               description=TRIM(each_desc_labels(i_each)),&
               usage=TRIM(each_possible_labels(i_each))//" <INTEGER>",&
               default_i_val=my_value, supported_feature=supported_feature, error=error)
          CALL section_add_keyword(subsection,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END DO
       CALL section_add_subsection(print_key_section,subsection,error=error)
       CALL section_release(subsection,error=error)

       my_add_last = add_last_no
       IF (PRESENT(add_last)) THEN
          my_add_last = add_last
       END IF
       CALL keyword_create(keyword, name="ADD_LAST",&
            description="If the last iteration should be added, and if it "//&
            "should be marked symbolically (with lowercase letter l) or with "//&
            "the iteration number. "//&
            "Not every iteration level is able to identify the last iteration "//&
            "early enough to be able to output. When this keyword is activated "//&
            "all iteration levels are checked for the last iteration step.",&
            usage="ADD_LAST (NO|NUMERIC|SYMBOLIC)",&
            enum_c_vals=s2a("no","numeric","symbolic"),&
            enum_i_vals=(/add_last_no, add_last_numeric, add_last_symbolic/),&
            enum_desc=s2a("Do not mark last iteration specifically",&
                            "Mark last iteration with its iteration number",&
                            "Mark last iteration with lowercase letter l"),&
            default_i_val=my_add_last,supported_feature=supported_feature,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       my_comm_iter_levels=0
       IF (PRESENT(common_iter_levels)) my_comm_iter_levels=common_iter_levels
       CALL keyword_create(keyword, name="COMMON_ITERATION_LEVELS",&
            description="How many iterations levels should be written"//&
            " in the same file (no extra information about the actual"//&
            " iteration level is written to the file)",&
            usage="COMMON_ITERATION_LEVELS <INTEGER>",&
            default_i_val=my_comm_iter_levels, supported_feature=supported_feature,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       my_filename=""
       IF (PRESENT(filename)) my_filename=filename
       CALL keyword_create(keyword, name="FILENAME",&
            description=' controls part of the filename for output. '//&
            ' use __STD_OUT__ (exactly as written here) for the screen or standard logger. '//&
            ' use filename to obtain projectname-filename. '//&
            ' use ./filename to get filename.'//&
            ' A middle name (if present), iteration numbers'//&
            ' and extension are always added to the filename.'//&
            ' if you want to avoid it use =filename, in this'//&
            ' case the filename is always exactly as typed.'//&
            ' Please note that this can lead to clashes of'//&
            ' filenames.',&
            usage="FILENAME ./filename ",&
            default_lc_val=my_filename, supported_feature=supported_feature,error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LOG_PRINT_KEY",&
            description="This keywords enables the logger for the print_key (a message is printed on "//&
            "screen everytime data, controlled by this print_key, are written)",&
            usage="LOG_PRINT_KEY <LOGICAL>", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,&
            supported_feature=supported_feature, error=error)
       CALL section_add_keyword(print_key_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       IF (PRESENT(unit_str)) THEN
          CALL keyword_create(keyword, name="UNIT",&
               description='Specify the unit of measurement for the quantity in output. '//&
               "All available CP2K units can be used.",&
               usage="UNIT angstrom",default_c_val=TRIM(unit_str),supported_feature=supported_feature,&
               error=error)
          CALL section_add_keyword(print_key_section,keyword,error=error)
          CALL keyword_release(keyword,error=error)
       END IF
    END IF
  END SUBROUTINE cp_print_key_section_create

! *****************************************************************************
!> \brief returns what should be done with the given proprety
!>      if btest(res,cp_p_store) then the property should be stored in memory
!>      if btest(res,cp_p_file) then the property should be print ed to a file
!>      if res==0 then nothing should be done
!> \param iteration_info information about the actual iteration level
!> \param basis_section section that contains the printkey
!> \param print_key_path path to the printkey- "%" between sections, and
!>        optionally a "/" and a logical flag to check). Might be empty.
!> \param used_print_key here the print_key that was used is returned
!> \param first_time if it ist the first time that an output is written
!>        (not fully correct, but most of the time)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
!> \note
!>      not all the propreties support can be stored
! *****************************************************************************
  FUNCTION cp_print_key_should_output(iteration_info,basis_section,&
       print_key_path,used_print_key,first_time,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: print_key_path
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: used_print_key
    LOGICAL, INTENT(OUT), OPTIONAL           :: first_time
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_should_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: end_str, my_control_val, &
                                                to_path
    LOGICAL                                  :: failure, flags, is_iter, is_on
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    res=0
    IF (PRESENT(first_time)) first_time=.FALSE.
    CPPrecondition(ASSOCIATED(basis_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(basis_section%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(used_print_key)) NULLIFY(used_print_key)
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
       flags=.TRUE.
    END IF
    IF (PRESENT(used_print_key)) used_print_key => print_key

    IF (.NOT.flags) RETURN

    CALL section_vals_val_get(print_key,"__CONTROL_VAL",&
         i_val=my_control_val,error=error)
    is_on=cp_printkey_is_on(iteration_info,print_key,error=error)

    ! a shortcut for most common case
    IF(my_control_val==cp_out_default .AND. .NOT. is_on) RETURN

    is_iter=cp_printkey_is_iter(iteration_info,print_key,first_time=first_time,&
         error=error)

    IF (BTEST(my_control_val,cp_p_store)) THEN
       res=IBSET(res,cp_p_store)
    ELSE IF (BTEST(my_control_val,cp_p_store_if).and.is_iter.and.is_on) THEN
       res=IBSET(res,cp_p_store)
    ELSE IF (BTEST(my_control_val,cp_p_store_each).and.is_iter) THEN
       res=IBSET(res,cp_p_store)
    END IF

    IF (BTEST(my_control_val,cp_p_file)) THEN
       res=IBSET(res,cp_p_file)
    ELSE IF (BTEST(my_control_val,cp_p_file_if).and.is_iter.and.is_on) THEN
       res=IBSET(res,cp_p_file)
    ELSE IF (BTEST(my_control_val,cp_p_file_each).and.is_iter) THEN
       res=IBSET(res,cp_p_file)
    END IF
    IF (BTEST(my_control_val,cp_p_calc).OR.res/=0) THEN
       res=IBSET(res,cp_p_calc)
    END IF
  END FUNCTION cp_print_key_should_output

! *****************************************************************************
!> \brief returns true if the printlevel activates this printkey
!>      does not look if this iteration it should be printed
!> \param iteration_info information about the actual iteration level
!> \param print_key the section values of the key to be printed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_printkey_is_on(iteration_info,print_key,error) RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_printkey_is_on', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: print_level
    LOGICAL                                  :: failure

    failure=.FALSE.
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

! *****************************************************************************
!> \brief returns if the actual iteration matches those selected by the
!>      given printkey. Does not check it the prinkey is active (at the
!>      actual print_level)
!> \param iteration_info information about the actual iteration level
!> \param print_key the section values of the key to be printed
!> \param first_time returns if it is the first time that output is written
!>        (not fully correct, but most of the time)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
  FUNCTION cp_printkey_is_iter(iteration_info,print_key,first_time,error)&
       RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key
    LOGICAL, INTENT(OUT), OPTIONAL           :: first_time
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_printkey_is_iter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: add_last, ilevel, iter_nr, &
                                                ival
    LOGICAL                                  :: failure, first, level_passed

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.ASSOCIATED(print_key)) THEN
       res=(iteration_info%print_level > debug_print_level)
       first=ALL(iteration_info%iteration(1:iteration_info%n_rlevel)==1)
    ELSE
       CPPrecondition(print_key%ref_count>0,cp_failure_level,routineP,error,failure)
       res   = .FALSE.
       first = .FALSE.
       IF (.NOT. failure) THEN
          CALL section_vals_val_get(print_key,"ADD_LAST",i_val=add_last,error=error)
          res  =.TRUE.
          first=.TRUE.
          DO ilevel=1,iteration_info%n_rlevel
             level_passed=.FALSE.
             CALL section_vals_val_get(print_key,"EACH%"//TRIM(iteration_info%level_name(ilevel)),&
                  i_val=ival,error=error)
             IF (ival>0) THEN
                iter_nr=iteration_info%iteration(ilevel)
                IF (iter_nr/ival>1) first=.FALSE.
                IF (MODULO(iter_nr,ival)==0) THEN
                   level_passed=.TRUE.
                END IF
             END IF
             IF (add_last==add_last_numeric .OR. add_last==add_last_symbolic) THEN
                IF (iteration_info%last_iter(ilevel)) THEN
                   level_passed=.TRUE.
                END IF
             END IF
             IF (.NOT.level_passed) res=.FALSE.
          END DO
       END IF
    END IF
    first=first.AND.res
    IF (PRESENT(first_time)) first_time=first
  END FUNCTION cp_printkey_is_iter

! *****************************************************************************
!> \brief returns the iteration string, a string that is useful to create
!>      unique filenames (once you trim it)
!> \param iter_info the iteration info from where to take the iteration
!>        number
!> \param print_key the print key to optionally show the last iteration
!>        symbolically
!> \param for_file if the string is to be used for file generation
!>        (and should consequently ignore some iteration levels depending
!>        on COMMON_ITERATION_LEVELS).
!>        Defaults to false.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
!> \note
!>      If the root level is 1 removes it
! *****************************************************************************
  FUNCTION cp_iter_string(iter_info,print_key,for_file,error) RESULT(res)
    TYPE(cp_iteration_info_type), POINTER    :: iter_info
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: print_key
    LOGICAL, INTENT(IN), OPTIONAL            :: for_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    CHARACTER(len=default_string_length)     :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iter_string', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: add_last, c_i_level, ilevel, &
                                                n_rlevel, s_level
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
       IF (ASSOCIATED(my_print_key)) THEN
          CALL section_vals_val_get(my_print_key,"ADD_LAST",i_val=add_last,error=error)
          CALL section_vals_val_get(my_print_key,"COMMON_ITERATION_LEVELS",i_val=c_i_level, error=error)
          n_rlevel=iter_info%n_rlevel
          IF (my_for_file) n_rlevel=MIN(n_rlevel,MAX(0,n_rlevel-c_i_level))
          DO ilevel=s_level,n_rlevel
             IF (iter_info%last_iter(ilevel)) THEN
                IF (add_last==add_last_symbolic) THEN
                   WRITE(res(9*ilevel-8:9*ilevel),"('l_')")
                ELSE
                   WRITE(res(9*ilevel-8:9*ilevel),"(i8,'_')") iter_info%iteration(ilevel)
                END IF
             ELSE
                WRITE(res(9*ilevel-8:9*ilevel),"(i8,'_')") iter_info%iteration(ilevel)
             END IF
          END DO
       ELSE
          DO ilevel=s_level,iter_info%n_rlevel
             WRITE(res(9*ilevel-8:9*ilevel),"(i8,'_')") iter_info%iteration(ilevel)
          END DO
       END IF
       CALL compress(res,.TRUE.)
       IF (LEN_TRIM(res)>0) THEN
          res(LEN_TRIM(res):LEN_TRIM(res))=" "
       END IF
    END IF
  END FUNCTION cp_iter_string

! *****************************************************************************
!> \brief adds one to the actual iteration
!> \param iteration_info the iteration info to update
!> \param last if this iteration is the last one (defaults to false)
!> \param iter_nr ...
!> \param increment ...
!> \param iter_nr_out ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
!> \note
!>      this is supposed to be called at the beginning of each iteration
! *****************************************************************************
  SUBROUTINE cp_iterate(iteration_info,last,iter_nr,increment,iter_nr_out,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    LOGICAL, INTENT(IN), OPTIONAL            :: last
    INTEGER, INTENT(IN), OPTIONAL            :: iter_nr, increment
    INTEGER, INTENT(OUT), OPTIONAL           :: iter_nr_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_iterate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_increment
    LOGICAL                                  :: failure, my_last

    failure=.FALSE.
    my_last=.FALSE.
    my_increment = 1
    IF (PRESENT(last))            my_last = last
    IF (PRESENT(increment))  my_increment = increment
    IF (PRESENT(iter_nr_out)) iter_nr_out = -1

    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
       IF (PRESENT(iter_nr)) THEN
          iteration_info%iteration(iteration_info%n_rlevel)=iter_nr
       ELSE
          iteration_info%iteration(iteration_info%n_rlevel)=&
               iteration_info%iteration(iteration_info%n_rlevel)+my_increment
       END IF
       ! If requested provide the value of the iteration level
       IF (PRESENT(iter_nr_out)) iter_nr_out = iteration_info%iteration(iteration_info%n_rlevel)

       ! Possibly setup the LAST flag
       iteration_info%last_iter(iteration_info%n_rlevel)=my_last
    END IF
  END SUBROUTINE cp_iterate

! *****************************************************************************
!> \brief Adds an iteration level
!> \param iteration_info the iteration info to which an iteration level has
!>        to be added
!> \param level_name the name of this level, for pretty printing only, right now
!> \param n_rlevel_new number of iteration levels after this call 
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_add_iter_level(iteration_info,level_name,n_rlevel_new,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    CHARACTER(LEN=*), INTENT(IN)             :: level_name
    INTEGER, INTENT(OUT), OPTIONAL           :: n_rlevel_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_add_iter_level', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, found

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       found = .FALSE.
       DO i = 1, SIZE(each_possible_labels)
          IF (TRIM(level_name)==TRIM(each_possible_labels(i))) THEN
             found = .TRUE.
             EXIT
          END IF
       END DO
       IF (found) THEN
          CALL cp_iteration_info_retain(iteration_info)
          iteration_info%n_rlevel=iteration_info%n_rlevel+1
          CALL reallocate(iteration_info%iteration,  1, iteration_info%n_rlevel)
          CALL reallocate(iteration_info%level_name, 1, iteration_info%n_rlevel)
          CALL reallocate(iteration_info%last_iter,  1, iteration_info%n_rlevel)
          iteration_info%iteration(iteration_info%n_rlevel)  = 0
          iteration_info%level_name(iteration_info%n_rlevel) = level_name
          iteration_info%last_iter(iteration_info%n_rlevel)  = .FALSE.
          IF (PRESENT(n_rlevel_new)) n_rlevel_new=iteration_info%n_rlevel
       ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Trying to create an iteration level ("//TRIM(level_name)//") not defined."//&
               "Please update the module: cp_iter_types."//&
               CPSourceFileRef)
       END IF
    END IF

  END SUBROUTINE cp_add_iter_level

! *****************************************************************************
!> \brief Removes an iteration level
!> \param iteration_info the iteration info to which an iteration level has
!>        to be removed
!> \param level_name level_name to be destroyed (if does not match gives an error)
!> \param n_rlevel_att iteration level before the call (to do some checks)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_rm_iter_level(iteration_info,level_name,n_rlevel_att,error)
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    CHARACTER(LEN=*), INTENT(IN)             :: level_name
    INTEGER, INTENT(IN), OPTIONAL            :: n_rlevel_att
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_rm_iter_level', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(iteration_info),cp_failure_level,routineP,error,failure)
    CPPrecondition(iteration_info%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(n_rlevel_att)) THEN
          CPPrecondition(n_rlevel_att==iteration_info%n_rlevel,cp_failure_level,routineP,error,failure)
       END IF
       CALL cp_iteration_info_release(iteration_info)
       ! This check that the iteration levels are consistently created and destroyed..
       ! Never remove this check..
       check = iteration_info%level_name(iteration_info%n_rlevel)==level_name
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       iteration_info%n_rlevel=iteration_info%n_rlevel-1
       CALL reallocate(iteration_info%iteration,  1, iteration_info%n_rlevel)
       CALL reallocate(iteration_info%level_name, 1, iteration_info%n_rlevel)
       CALL reallocate(iteration_info%last_iter,  1, iteration_info%n_rlevel)
    END IF
  END SUBROUTINE cp_rm_iter_level

! *****************************************************************************
!> \brief Utility function that retuns a unit number to write the print key.
!>     Might open a file with a unique filename, generated from
!>     the print_key name and iteration info.
!>
!>     Normally a valid unit (>0) is returned only if cp_print_key_should_output
!>     says that the print_key should be printed, and if the unit is global
!>     only the io node has a valid unit.
!>     So in many cases you can decide if you should print just checking if
!>     the returned units is bigger than 0.
!>
!>     IMPORTANT you should call cp_finished_output when an iteration output is
!>     finished (to immediately close the file that might have been opened)
!> \param logger the logger for the parallel environment, iteration info
!>        and filename generation
!> \param print_key ...
!> \param middle_name name to be added to the generated filename, useful when
!>        print_key activates different distinct outputs, to be able to
!>        distinguish them
!> \param extension extension to be applied to the filename (including the ".")
!> \param my_local if the unit should be local to this task, or global to the
!>        program (defaults to false).
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval filename ...
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_print_key_generate_filename(logger,print_key,middle_name,extension,&
       my_local, error) RESULT(filename)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: print_key
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: middle_name
    CHARACTER(len=*), INTENT(IN)             :: extension
    LOGICAL, INTENT(IN)                      :: my_local
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    CHARACTER(len=default_path_length)       :: filename

    CHARACTER(len=default_path_length)       :: outPath, postfix, root
    CHARACTER(len=default_string_length)     :: my_middle_name, outName
    INTEGER                                  :: my_ind1, my_ind2
    LOGICAL                                  :: has_root

    CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)
    IF (outPath(1:1)=='=') THEN
       filename=outPath(2:LEN_TRIM(outPath))
       RETURN
    END IF
    IF (outPath=="__STD_OUT__") outPath=""
    outName=outPath
    has_root=.FALSE.
    my_ind1=INDEX(outPath,"/")
    my_ind2 = LEN_TRIM(outPath)
    IF (my_ind1 /= 0) THEN
       has_root = .TRUE.
       DO WHILE (INDEX(outPath(my_ind1+1:my_ind2),"/")/=0)
          my_ind1 = INDEX(outPath(my_ind1+1:my_ind2),"/") + my_ind1
       END DO
       IF (my_ind1 == my_ind2) THEN
          outName=""
       ELSE
          outName = outPath(my_ind1+1:my_ind2)
       END IF
    END IF

    IF (PRESENT(middle_name)) THEN
       IF (outName/="") THEN
          my_middle_name="-"//TRIM(outName)//"-"//middle_name
       ELSE
          my_middle_name="-"//middle_name
       END IF
    ELSE
       IF (outName/="") THEN
          my_middle_name="-"//TRIM(outName)
       ELSE
          my_middle_name = ""
       END IF
    ENDIF

    IF (.not.has_root) THEN
       root=TRIM(logger%iter_info%project_name)//TRIM(my_middle_name)
    ELSE IF (outName=="") THEN
       root=outPath(1:my_ind1)//TRIM(logger%iter_info%project_name)//TRIM(my_middle_name)
    ELSE
       root=outPath(1:my_ind1)//my_middle_name(2:LEN_TRIM(my_middle_name))
    END IF

    ! use the cp_iter_string as a postfix
    postfix="-"//TRIM(cp_iter_string(logger%iter_info,print_key=print_key,for_file=.TRUE.,error=error))
    IF (TRIM(postfix)=="-") postfix=""

    ! and add the extension
    postfix=TRIM(postfix)//extension
    ! and let the logger generate the filename
    CALL cp_logger_generate_filename(logger,res=filename,&
         root=root, postfix=postfix,local=my_local)

  END FUNCTION cp_print_key_generate_filename

! *****************************************************************************
!> \brief ...
!> \param logger ...
!> \param basis_section ...
!> \param print_key_path ...
!> \param extension ...
!> \param middle_name ...
!> \param local ...
!> \param log_filename ...
!> \param ignore_should_output ...
!> \param file_form ...
!> \param file_position ...
!> \param file_action ...
!> \param file_status ...
!> \param do_backup ...
!> \param on_file ...
!> \param is_new_file ...
!> \param error ...
!> \retval res ...
! *****************************************************************************
  FUNCTION cp_print_key_unit_nr(logger, basis_section, print_key_path, extension,&
       middle_name, local, log_filename,ignore_should_output, file_form, file_position,&
       file_action, file_status, do_backup, on_file, is_new_file,error) RESULT(res)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: print_key_path
    CHARACTER(len=*), INTENT(IN)             :: extension
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: middle_name
    LOGICAL, INTENT(IN), OPTIONAL            :: local, log_filename, &
                                                ignore_should_output
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: file_form, file_position, &
                                                file_action, file_status
    LOGICAL, INTENT(IN), OPTIONAL            :: do_backup, on_file
    LOGICAL, INTENT(OUT), OPTIONAL           :: is_new_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_unit_nr', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: filename, filename_bak, &
                                                filename_bak_1, filename_bak_2
    CHARACTER(len=default_string_length)     :: my_file_action, my_file_form, &
                                                my_file_position, &
                                                my_file_status, outPath
    INTEGER                                  :: c_i_level, f_backup_level, i, &
                                                my_backup_level, my_nbak, &
                                                nbak, s_backup_level
    LOGICAL                                  :: do_log, failure, found, &
                                                my_do_backup, my_local, &
                                                my_on_file, my_should_output
    TYPE(cp_iteration_info_type), POINTER    :: iteration_info
    TYPE(section_vals_type), POINTER         :: print_key

    failure      = .FALSE.
    my_local     = .FALSE.
    my_do_backup = .FALSE.
    found        = .FALSE.
    res          = -1
    my_file_form     = "FORMATTED"
    my_file_position = "APPEND"
    my_file_action   = "WRITE"
    my_file_status   = "UNKNOWN"
    my_on_file       = .FALSE.
    IF (PRESENT(file_form))     my_file_form     = file_form
    IF (PRESENT(file_position)) my_file_position = file_position
    IF (PRESENT(file_action))   my_file_action   = file_action
    IF (PRESENT(file_status))   my_file_status   = file_status
    IF (PRESENT(do_backup))     my_do_backup     = do_backup
    IF (PRESENT(on_file))       my_on_file       = on_file
    IF (PRESENT(local))         my_local         = local
    NULLIFY(print_key)
    CPPrecondition(ASSOCIATED(basis_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(basis_section%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       my_should_output=BTEST(cp_print_key_should_output(logger%iter_info,&
            basis_section,print_key_path,used_print_key=print_key,error=error),cp_p_file)
       IF (PRESENT(ignore_should_output)) my_should_output=my_should_output.or.ignore_should_output
       IF (.NOT.my_should_output) RETURN
       IF (my_local.OR.&
            logger%para_env%mepos==logger%para_env%source) THEN

          CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)
          IF (outPath=='__STD_OUT__'.AND..NOT.my_on_file) THEN
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
             ! Give back info about a possible existence of the file if required
             IF (PRESENT(is_new_file)) THEN
                INQUIRE(FILE=filename,EXIST=found)
                is_new_file = .NOT.found
                IF (my_file_position=="REWIND") is_new_file = .TRUE.
             END IF
             ! Check is we have to log any operation performed on the file..
             do_log = .FALSE.
             IF (PRESENT(log_filename)) THEN
                do_log = log_filename
             ELSE
                CALL section_vals_val_get(print_key,"LOG_PRINT_KEY",l_val=do_log,error=error)
             END IF
             ! If required do a backup
             IF (my_do_backup) THEN
                INQUIRE(FILE=filename,EXIST=found)
                CALL section_vals_val_get(print_key,"BACKUP_COPIES",i_val=nbak,error=error)
                IF (nbak/=0) THEN
                   iteration_info => logger%iter_info
                   s_backup_level   = 0
                   IF (ASSOCIATED(print_key%ibackup)) s_backup_level = SIZE(print_key%ibackup)
                   CALL section_vals_val_get(print_key,"COMMON_ITERATION_LEVELS",i_val=c_i_level, error=error)
                   my_backup_level = MAX(1,iteration_info%n_rlevel-c_i_level+1)
                   f_backup_level  = MAX(s_backup_level,my_backup_level)
                   IF (f_backup_level>s_backup_level) THEN
                      CALL reallocate(print_key%ibackup,1,f_backup_level)
                      DO i = s_backup_level+1, f_backup_level
                         print_key%ibackup(i) = 0
                      END DO
                   END IF
                   IF (found) THEN
                      print_key%ibackup(my_backup_level)=print_key%ibackup(my_backup_level)+1
                      my_nbak = print_key%ibackup(my_backup_level)
                      ! Recent backup copies correspond to lower backup indexes
                      DO i = MIN(nbak,my_nbak), 2, -1
                         filename_bak_1=TRIM(filename)//".bak-"//ADJUSTL(cp_to_string(i))
                         filename_bak_2=TRIM(filename)//".bak-"//ADJUSTL(cp_to_string(i-1))
                         IF (do_log) THEN
                            CALL cp_log(logger=logger, level=cp_note_level, fromWhere=routineP,&
                                 message="Moving file "//TRIM(filename_bak_2)//" into file "//&
                                 TRIM(filename_bak_1)//".",local=my_local)
                         END IF
                         INQUIRE(FILE=filename_bak_2,EXIST=found)
                         IF (.NOT.found) THEN
                            IF (do_log) THEN
                               CALL cp_log(logger=logger, level=cp_note_level, fromWhere=routineP,&
                                    message="File "//TRIM(filename_bak_2)//" not existing..",local=my_local)
                            END IF
                         ELSE
                            CALL m_mov(TRIM(filename_bak_2), TRIM(filename_bak_1))
                         END IF
                      END DO
                      ! The last backup is always the one with index 1
                      filename_bak=TRIM(filename)//".bak-"//ADJUSTL(cp_to_string(1))
                      IF (do_log) THEN
                         CALL cp_log(logger=logger, level=cp_note_level, fromWhere=routineP,&
                              message="Moving file "//TRIM(filename)//" into file "//&
                              TRIM(filename_bak)//".",local=my_local)
                      END IF
                      CALL m_mov(TRIM(filename), TRIM(filename_bak))
                   ELSE
                      ! Zero the backup history for this new iteration level..
                      print_key%ibackup(my_backup_level)= 0
                   END IF
                END IF
             END IF
             CALL open_file(file_name=filename,file_status=my_file_status,&
                  file_form=my_file_form,file_action=my_file_action,&
                  file_position=my_file_position,unit_number=res)
             IF (do_log) THEN
                CALL cp_log(logger=logger, level=cp_note_level, fromWhere=routineP,&
                     message="Writing "//TRIM(print_key%section%name)//" "//&
                     TRIM(cp_iter_string(logger%iter_info,error=error))//" to "//&
                     TRIM(filename),local=my_local)
             END IF
          END IF
       ELSE
          res=-1
       END IF
    END IF
  END FUNCTION cp_print_key_unit_nr

! *****************************************************************************
!> \brief should be called after you finish working with a unit obtained with
!>      cp_print_key_unit_nr, so that the file that might have been opened
!>      can be closed.
!>
!>      the inputs should be exactly the same of the corresponding
!>      cp_print_key_unit_nr
!> \param unit_nr ...
!> \param logger ...
!> \param basis_section ...
!> \param print_key_path ...
!> \param local ...
!> \param ignore_should_output ...
!> \param on_file ...
!> \param error ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      closes if the corresponding filename of the printkey is
!>      not __STD_OUT__
! *****************************************************************************
  SUBROUTINE cp_print_key_finished_output(unit_nr, logger, basis_section,&
       print_key_path,local,ignore_should_output,on_file,error)
    INTEGER, INTENT(INOUT)                   :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: print_key_path
    LOGICAL, INTENT(IN), OPTIONAL            :: local, ignore_should_output, &
                                                on_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_finished_output', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: outPath
    LOGICAL                                  :: failure, my_local, &
                                                my_on_file, my_should_output
    TYPE(section_vals_type), POINTER         :: print_key

    failure=.FALSE.
    my_local=.FALSE.
    my_on_file=.FALSE.
    NULLIFY(print_key)
    IF (PRESENT(local)) my_local=local
    IF (PRESENT(on_file)) my_on_file=on_file
    CPPrecondition(ASSOCIATED(basis_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(logger),cp_failure_level,routineP,error,failure)
    CPPrecondition(basis_section%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(logger%ref_count>0,cp_failure_level,routineP,error,failure)
    my_should_output=BTEST(cp_print_key_should_output(logger%iter_info,basis_section,&
         print_key_path,used_print_key=print_key,error=error),cp_p_file)
    IF (PRESENT(ignore_should_output)) my_should_output=my_should_output.or.ignore_should_output
    IF (my_should_output.and.(my_local.OR.&
         logger%para_env%source==logger%para_env%mepos)) THEN
       CALL section_vals_val_get(print_key,"FILENAME",c_val=outPath,error=error)
       IF (my_on_file.OR.outPath.NE.'__STD_OUT__') THEN
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

! *****************************************************************************
!> \brief utility method to print out a string
!> \param logger the logger for the parallel environment, iteration info
!>        and filename generation
!> \param basis_section section that contains the printkey
!> \param print_key_path path to the printkey- "%" between sections, and
!>        optionally a "/" and a logical flag to check). Might be empty.
!> \param extension extension to be applied to the filename (including the ".")
!> \param message the message to write
!> \param middle_name name to be added to the generated filename, useful when
!>        print_key activates different distinct outputs, to be able to
!>        distinguish them
!> \param local if the unit should be local to this task, or global to the
!>        program (defaults to false).
!> \param log_filename if a small log with the filename should be written
!>        to the main log (defaults to false)
!> \param ignore_should_output if true always returns a valid unit (ignoring
!>        cp_print_key_should_output)
!> \param on_file ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  SUBROUTINE cp_print_key_log(logger, basis_section, print_key_path, extension,&
       message,middle_name, local, log_filename,ignore_should_output, &
       on_file, error)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: basis_section
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: print_key_path
    CHARACTER(len=*), INTENT(IN)             :: extension, message
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: middle_name
    LOGICAL, INTENT(IN), OPTIONAL            :: local, log_filename, &
                                                ignore_should_output, on_file
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_print_key_log', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit
    LOGICAL                                  :: failure, my_local

    failure=.FALSE.

    IF (.NOT. failure) THEN
       my_local=.FALSE.
       IF (PRESENT(local)) my_local=local
       IF (my_local .OR. logger%para_env%mepos==logger%para_env%source) THEN
          output_unit=cp_print_key_unit_nr(logger, basis_section=basis_section,&
               print_key_path=print_key_path,extension=extension,&
               middle_name=middle_name, local=local, log_filename=log_filename,&
               ignore_should_output=ignore_should_output, &
               on_file=on_file, error=error)
          IF (output_unit>0) THEN
             WRITE(output_unit,"(a)")message
             CALL cp_print_key_finished_output(output_unit,logger,&
                  basis_section=basis_section,print_key_path=print_key_path,&
                  local=local,ignore_should_output=ignore_should_output,&
                  on_file=on_file,error=error)
          END IF
       END IF
    END IF
  END SUBROUTINE cp_print_key_log

END MODULE cp_output_handling

