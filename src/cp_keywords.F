!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****s* cp2k/cp_keywords *
!!
!!   NAME
!!     cp_keywords
!!
!!   FUNCTION
!!     Defines a structure for defining keywords in the input, listing them
!!     and providing basic help.
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   HISTORY
!!     Creation [05.2004]
!!
!!   SOURCE
!******************************************************************************

MODULE cp_keywords
  USE kinds, ONLY: dp,default_string_length
  IMPLICIT NONE

  INTEGER, SAVE, PRIVATE :: last_section_id=0,last_keyword_id=0,&
       last_section_vals_id=0

  INTEGER, PARAMETER, PUBLIC :: no_parse_t=0,logical_t=1, &
       integer_t=2, real_t=3,char_t=4,char_integer_t=1
  INTEGER, PARAMETER, PUBLIC :: description_string_length=10*default_string_length
  
  !!****s* cp_keywords/kval_p_type [1.0] *
  !!
  !!   NAME
  !!     kval_p_type
  !!
  !!   FUNCTION
  !!     pointer to a kval, do create arrays of pointers
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - kval: to pointer to the kval
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE kval_p_type
     TYPE(kval_type), pointer :: kval
  END TYPE kval_p_type
  !!***
  !****************************************************************************

  !!****s* cp_keywords/kval_type *
  !!
  !!   NAME
  !!     kval_type
  !!
  !!   FUNCTION
  !!     low level type that stores the values parsed from a file
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - l_val, i_val, c_val, r_val: arrays with logical,integer,character
  !!       or real values. Might be unassociated if unused
  !!     - next_k_val: pointer to the next kval (if present)
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE kval_type
     LOGICAL :: popped_l,popped_i,popped_c,popped_r
     LOGICAL, DIMENSION(:), POINTER :: l_val
     INTEGER, DIMENSION(:), POINTER :: i_val
     CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: &
          c_val
     REAL(kind=dp), DIMENSION(:), POINTER :: r_val
     TYPE(kval_type), POINTER :: next_kval
  END TYPE kval_type
  !!***
  !****************************************************************************

  !!****s* cp_keywords/keyword_p_type *
  !!
  !!   NAME
  !!     keyword_p_type
  !!
  !!   FUNCTION
  !!     represent a pointer to a keyword (to make arrays of pointers)
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - keyword: the pointer to the keyword
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE keyword_p_type
     TYPE(keyword_type), pointer :: keyword
  END TYPE keyword_p_type
  !!***
  !****************************************************************************

  !!****s* cp_keywords/keyword_type *
  !!
  !!   NAME
  !!     keyword_type
  !!
  !!   FUNCTION
  !!     represent a keyword in the input
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - names: the names of the current keyword (at least one should be
  !!       present) for example "MAXSCF"
  !!     - usage: how to use it "MAXSCF 10"
  !!     - description: what does it do: "MAXSCF : determines the maximum 
  !!       number of steps in an SCF run"
  !!     - type_of_var: the type of keyword (controls how it is parsed)
  !!       it can be one of: no_parse_t,logical_t, integer_t, real_t,
  !!       char_t,char_integer_t
  !!     - n_var: number of values that should be parsed (-1=unknown)
  !!     - repeats: if the keyword can be present more than once in the
  !!       section
  !!     - required: if the keyword is required (leaving it out will give an
  !!       error)
  !!     - default_value: the default value for the keyword
  !!
  !!   AUTHOR
  !!     Joost & fawzi
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE keyword_type
     private
     INTEGER :: ref_count,id_nr
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER:: names
     CHARACTER(LEN=default_string_length) :: usage
     CHARACTER(LEN=description_string_length) :: description
     INTEGER :: type_of_var,n_var
     LOGICAL :: repeats, required
     
     TYPE(kval_type), POINTER :: default_value
  END TYPE keyword_type
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
     TYPE(section_type), pointer :: keyword
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
  !!     - parsed: if the section has been parsed (and no keyword/subsections
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
     PRIVATE
     LOGICAL :: parsed, required, repeats
     INTEGER :: id_nr,ref_count, n_keywords
     character(len=default_string_length) :: name
     CHARACTER(len=description_string_length) :: description
     TYPE(keyword_p_type), DIMENSION(:), POINTER :: keywords
     TYPE(section_p_type), POINTER, DIMENSION(:) :: subsections
  END TYPE section_type
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
     TYPE(kval_p_type), DIMENSION(:,:), POINTER :: values
     TYPE(section_vals_p_type), DIMENSION(:,:), POINTER :: subs_vals
  END TYPE section_vals_type
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

  !!****s* cp_keywords/input_file_type *
  !!
  !!   NAME
  !!     input_file_type
  !!
  !!   FUNCTION
  !!     represent an input file
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
  TYPE input_file_type
     PRIVATE
     INTEGER :: id_nr,ref_count
     character(len=default_path_length) :: path
     TYPE(section_p_type), DIMENSION(:), POINTER :: section
  END TYPE input_file_type
  !!***
  !****************************************************************************


  PUBLIC :: input_file_type,section_type,section_p_type,&
       keyword_type,keyword_p_type, kval_type,kval_p_type
  PUBLIC :: section_create,section_retain,&
       section_release, section_describe,section_parse,section_get_keyword,&
       section_pop_kval,section_get_kval,section_write,&
       section_add_keyword,section_add_
  PUBLIC :: keyword_create,keyword_retain,keyword_release,keyword_get,&
       keyword_set,keyword_pop

CONTAINS
  
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE section_create(section,name,description,n_keywords,repeats,&
     required,error)
  TYPE(section_type), pointer :: section
  CHARACTER(len=*), INTENT(in) :: name,description
  INTEGER, INTENT(in), OPTIONAL :: n_keywords
  LOGICAL, INTENT(in), OPTIONAL :: repeats,required
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: my_n_keywords,stat,i

  failure=.false.

  CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     my_n_keywords=10
     IF (PRESENT(n_keywords)) my_n_keywords=n_keywords

     ALLOCATE(section,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     last_section_id=last_section_id+1
     section%id_nr=last_section_id
     section%ref_count=1

     section%n_keywords=0
     section%name=name
     call uppercase(section%name)
     CPPrecondition(LEN_TRIM(description)<=LEN(section%description),cp_warning_level,routineP,error,failure)
     section%description=description
     section%parsed=.FALSE.
     section%required=.TRUE.
     section%repeats=.FALSE.
     IF (PRESENT(required)) section%required=required
     IF (PRESENT(repeats)) section%repeats=repeats
     NULLIFY(section%subsections)

     ALLOCATE(section%keywords(-1:my_n_keywords),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO i=-1,my_n_keywords
        NULLIFY(section%keywords(i)%keyword)
     END DO
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
     CPPreconditionNoFail(section%ref_count>0,cp_failure_level,routineP,error,failure)
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
     CPPreconditionNoFail(section%ref_count>0,cp_failure_level,routineP,error,failure)
     section%ref_count=section%ref_count-1
     IF (section%ref_count==0) THEN
        IF (ASSOCIATED(section%keywords)) THEN
           DO i=1,SIZE(section%keywords)
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

!!****f* keywords/kval_create *
!!
!!   NAME
!!     kval_create
!!
!!   FUNCTION
!!     creates a keyword value
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kval: the object to be created
!!     - logical_val: a logical to be stored in the kval
!!     - logical_vals: an array of logicals to be stored in kval
!!     - logical_vals_ptr: an array of logicals 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE kval_create(kval,l_val,l_vals,l_vals_ptr,i_val,i_vals,i_vals_ptr,&
     r_val,r_vals,r_vals_ptr,c_val,c_vals,c_vals_ptr,error)
  TYPE(kval_type), POINTER :: kval
  LOGICAL, OPTIONAL INTENT(in) :: l_val
  LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: l_vals
  LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals_ptr
  INTEGER, OPTIONAL INTENT(in) :: i_val
  INTEGER, DIMENSION(:), OPTIONAL, INTENT(in) :: i_vals
  INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals_ptr
  REAL(KIND=DP), OPTIONAL INTENT(in) :: r_val
  REAL(KIND=DP), DIMENSION(:), OPTIONAL, INTENT(in) :: r_vals
  REAL(KIND=DP), DIMENSION(:), OPTIONAL, POINTER :: r_vals_ptr
  CHARACTER(LEN=DEFAULT_STRING_LENGTH), OPTIONAL INTENT(in) :: c_val
  CHARACTER(LEN=DEFAULT_STRING_LENGTH), DIMENSION(:), OPTIONAL, INTENT(in) :: c_vals
  CHARACTER(LEN=DEFAULT_STRING_LENGTH), DIMENSION(:), OPTIONAL, POINTER :: c_vals_ptr
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kval_create',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(kval),cp_failure_level,routineP,error,failure)
  ALLOCATE(kval,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(kval%l_val,kval%i_val,kval%r_val,kval%c_val,kval%next_kval)
     kval%popped_l=.false.
     kval%popped_i=.false.
     kval%popped_r=.false.
     kval%popped_c=.false.

     IF (PRESENT(l_val)) THEN
        CPPrecondition(.NOT.PRESENT(l_vals),cp_failure_level,routineP,error,failure)
        CPPrecondition(.NOT.PRESENT(l_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%l_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%l_val(1)=l_val
     END IF
     IF (PRESENT(l_vals)) THEN
        CPPrecondition(.NOT.PRESENT(l_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%l_val(SIZE(l_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%l_val=l_vals
     END IF
     IF (PRESENT(l_vals_ptr)) THEN
        kval%l_val => l_vals_ptr
     END IF
     
     IF (PRESENT(r_val)) THEN
        CPPrecondition(.NOT.PRESENT(r_vals),cp_failure_level,routineP,error,failure)
        CPPrecondition(.NOT.PRESENT(r_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%r_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%r_val(1)=r_val
     END IF
     IF (PRESENT(r_vals)) THEN
        CPPrecondition(.NOT.PRESENT(r_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%r_val(SIZE(r_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%r_val=r_vals
     END IF
     IF (PRESENT(r_vals_ptr)) THEN
        kval%r_val => r_vals_ptr
     END IF
     
     IF (PRESENT(i_val)) THEN
        CPPrecondition(.NOT.PRESENT(i_vals),cp_failure_level,routineP,error,failure)
        CPPrecondition(.NOT.PRESENT(i_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%i_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%i_val(1)=i_val
     END IF
     IF (PRESENT(i_vals)) THEN
        CPPrecondition(.NOT.PRESENT(i_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%i_val(SIZE(i_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%i_val=i_vals
     END IF
     IF (PRESENT(i_vals_ptr)) THEN
        kval%i_val => i_vals_ptr
     END IF
     
     IF (PRESENT(c_val)) THEN
        CPPrecondition(.NOT.PRESENT(c_vals),cp_failure_level,routineP,error,failure)
        CPPrecondition(.NOT.PRESENT(c_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%c_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%c_val(1)=c_val
     END IF
     IF (PRESENT(c_vals)) THEN
        CPPrecondition(.NOT.PRESENT(c_vals_ptr),cp_failure_level,routineP,error,failure)
        ALLOCATE(kval%c_val(SIZE(c_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        kval%c_val=c_vals
     END IF
     IF (PRESENT(c_vals_ptr)) THEN
        kval%c_val => c_vals_ptr
     END IF
     
  END IF
END SUBROUTINE kval_create
!***************************************************************************

!!****f* keywords/kval_dealloc *
!!
!!   NAME
!!     kval_dealloc
!!
!!   FUNCTION
!!     deallocates the given kval (and all the connected kvals)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kval: the kval to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE kval_dealloc(kval,error)
  TYPE(kval_type), pointer :: kval
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='kval_dealloc',&
       routineP=moduleN//':'//routineN
  INTEGER :: stat
  TYPE(kval_type), pointer :: o_kval

  failure=.false.

  CPPrecondition(ASSOCIATED(kval),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO WHILE(ASSOCIATED(kval))
        IF (ASSOCIATED(kval%l_val)) THEN
           DEALLOCATE(kval%l_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(kval%i_val)) THEN
           DEALLOCATE(kval%i_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(kval%r_val)) THEN
           DEALLOCATE(kval%r_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(kval%c_val)) THEN
           DEALLOCATE(kval%c_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        o_kval => kval
        kval => kval%next_kval
        DEALLOCATE(o_kval,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END DO
  END IF
END SUBROUTINE kval_dealloc
!***************************************************************************

!!****f* keywords/keyword_create *
!!
!!   NAME
!!     keyword_create
!!
!!   FUNCTION
!!     creates a keyword object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - keyword: the keyword object to be created
!!     - name: the name of the keyword
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE keyword_create(keyword, name, type_of_var, description, usage,&
     nvar,repeats,variants,repeats,required,default_kval_ptr,&
     default_l_val, default_r_val, default_c_val, default_i_val,&
     default_l_vals, default_r_vals, default_c_vals, default_i_vals,&
     error)
  TYPE(keyword_type), pointer :: keyword
  CHARACTER(len=*), INTENT(in) :: name
  CHARACTER(len=*), INTENT(in) :: description
  CHARACTER(len=*), INTENT(in) :: usage
  INTEGER, INTENT(in),OPTIONAL :: type_of_var
  INTEGER, INTENT(in), optional :: nvar
  LOGICAL, INTENT(in), optional :: repeats
  CHARACTER(len=*), INTENT(in),DIMENSION(:),optional :: variants
  TYPE(kval_type), POINTER, optional :: default_kval_ptr
  LOGICAL, OPTIONAL INTENT(in) :: default_l_val
  LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: default_l_vals
  INTEGER, OPTIONAL INTENT(in) :: default_i_val
  INTEGER, DIMENSION(:), OPTIONAL, INTENT(in) :: default_i_vals
  REAL(KIND=DP), OPTIONAL INTENT(in) :: default_r_val
  REAL(KIND=DP), DIMENSION(:), OPTIONAL, INTENT(in) :: default_r_vals
  CHARACTER(LEN=*), OPTIONAL INTENT(in) :: default_c_val
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL, INTENT(in) :: default_c_vals
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='keyword_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i

  failure=.false.
  
  CPPrecondition(.NOT.ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  ALLOCATE(keyword,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     keyword%ref_count=1
     last_keyword_id=last_keyword_id+1
     keyword%id_nr=last_keyword_id

     IF (PRESENT(variants)) THEN
        ALLOCATE(keyword%names(SIZE(variants)+1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        keyword%names(1)=name
        DO i=1,SIZE(variants)
           keyword%names(i+1)=variants(i)
        END DO
     ELSE
        ALLOCATE(keyword%names(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        keyword%names(1)=name
     END IF
     DO i=1,SIZE(keyword%names)
        CALL uppercase(keyword%names(i))
     END DO
     
     CPPrecondition(LEN_TRIM(usage)<=LEN(keyword%usage),cp_failure_level,routineP,error,failure)
     keyword%usage=usage
     CPPrecondition(LEN_TRIM(description)<=LEN(keyword%description),cp_failure_level,routineP,error,failure)
     keyword%description=description
     
     keyword%repeats=.false.
     if (present(repeats)) keyword%repeats=repeats
     keyword%required=.false.
     if (present(required)) keyword%required=required

     NULLIFY(keyword%default_value)
     IF (PRESENT(default_kval_ptr)) THEN
        CALL cp_assert(.NOT.(PRESENT(default_l_val).OR.PRESENT(default_l_vals).or.&
             PRESENT(default_i_val).OR.PRESENT(default_i_vals).OR.&
             PRESENT(default_r_val).or.present(default_r_vals).or.&
             PRESENT(default_c_val).OR.PRESENT(default_c_vals)),cp_failure_level,&
             cp_asserion_failed,routineP,&
             "you should pass either default_kval_ptr or a default value, not both",&
             error,failure)
        keyword%default_value => default_kval_ptr
     END IF
     if (.not.associated(keyword%default_value)) then
        CALL kval_create(keyword%default_value,l_val=default_l_val,&
             l_vals=default_l_vals,i_val=default_i_val,i_vals=default_i_vals,&
             r_val=default_r_val,r_vals=default_r_vals,c_val=default_c_val,&
             c_vals=default_c_vals)
     END IF
     keyword%type_of_var=no_parse_t
     keyword%n_var=0
     IF (ASSOCIATED(keyword%default_value%l_val)) THEN
        keyword%n_var=keyword%n_var+size(keyword%default_value%l_val)
        keyword%type_of_var=logical_t
     END IF
     IF (ASSOCIATED(keyword%default_value%i_val)) THEN
        keyword%n_var=keyword%n_var+SIZE(keyword%default_value%i_val)
        keyword%type_of_var=integer_t
     END IF
     IF (ASSOCIATED(keyword%default_value%r_val)) THEN
        keyword%n_var=keyword%n_var+size(keyword%default_value%r_val)
        keyword%type_of_var=real_t
     END IF
     IF (ASSOCIATED(keyword%default_value%c_val)) THEN
        keyword%n_var=keyword%n_var+size(keyword%default_value%c_val)
        IF (ASSOCIATED(keyword%default_value%i_val)) THEN
           keyword%type_of_var=char_integer_t
        ELSE
           keyword%type_of_var=char_t
        END IF
     END IF
     IF (PRESENT(type_of_var)) THEN
        keyword%type_of_var=type_of_var
     END IF
  END IF
END SUBROUTINE keyword_create
!***************************************************************************

!!****f* keywords/keyword_retain *
!!
!!   NAME
!!     keyword_retain
!!
!!   FUNCTION
!!     retains the given keyword (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - keyword: the keyword to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE keyword_retain(keyword, error)
  TYPE(keyword_type), pointer :: keyword
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='keyword_retain',&
        routineP=moduleN//':'//routineN

  failure=.false.
  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
     keyword%ref_count=keyword%ref_count+1
  END IF
END SUBROUTINE keyword_retain
!***************************************************************************

!!****f* cp_keywords/keyword_release *
!!
!!   NAME
!!     keyword_release
!!
!!   FUNCTION
!!     releases the given keyword (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - keyword: the keyword to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE keyword_release(keyword, error)
  TYPE(keyword_type), pointer :: keyword
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='keyword_release',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  IF (associated(keyword) THEN
     CPPreconditionNoFail(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
     keyword%ref_count=keyword%ref_count-1
     IF (keyword%ref_count==0) THEN
        DEALLOCATE(keyword%names,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CALL kval_dealloc(ASSOCIATED(keyword%default_value),error=error)
        DEALLOCATE(keyword,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(keyword)
END SUBROUTINE keyword_release
!***************************************************************************

!!****f* cp_keywords/keyword_get *
!!
!!   NAME
!!     keyword_get
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
SUBROUTINE keyword_get(keyword,names,usage,description,type_of_var,n_var,&
     default_value, error)
  TYPE(keyword_type), POINTER :: keyword
  CHARACTER(len=default_string_length), DIMENSION(:), POINTER, OPTIONAL :: names
  CHARACTER(len=*), INTENT(out), optional :: usage, description
  INTEGER, INTENT(out), OPTIONAL :: type_of_var,n_var
  LOGICAL, INTENT(out), OPTIONAL :: repeats, required
  TYPE(kval_type), POINTER, OPTIONAL :: default_value
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='keyword_get',&
        routineP=moduleN//':'//routineN
  failure=.false.

  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(names)) names => keyword%names
     if (present(usage)) usage=keyword%usage
     if (present(description)) description=keyword%description
     if (present(type_of_var)) type_of_var=keyword%type_of_var
     if (present(n_var)) n_var=keyword%n_var
     if (present(repeats)) repeats=keyword%repeats
     if (present(required)) required=keyword%required
     if (present(default_value)) default_value => keyword%default_value
  end if
END SUBROUTINE keyword_get
!***************************************************************************

!!****f* cp_keywords/keyword_describe [1.0] *
!!
!!   NAME
!!     keyword_describe
!!
!!   FUNCTION
!!     writes out a description of the keyword
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - keyword: the keyword to describe
!!     - unit_nr: the unit to write to
!!     - level: the description level (0 no description, 1 name
!!       2: +usage, 3: +variants+description+default_value+required+repeats
!!       4: +type_of_var)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE keyword_describe(keyword, unit_nr, level,error)
  TYPE(keyword_type), pointer :: keyword
  INTEGER, INTENT(in) :: unit_nr,level
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  character(len=*), parameter :: routineN='keyword_describe',&
       routineP=moduleN//':'//routineN

  failure=.false.

  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure.AND.level>0) THEN
     WRITE(unit_nr,"(a,a,a)") "                                        ---",&
          TRIM(keyword%names(1)),"---"
     IF (level>1) THEN
        WRITE(unit_nr,"(a,a)") "usage         : ",TRIM(keyword%usage)
     END IF
     IF (level>2) THEN
        WRITE(unit_nr,"(a,a)") "description   : ",TRIM(keyword%description)
        SELECT CASE(keyword%type_of_var)
        CASE (logical_t)
           IF (level>3) THEN
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of logicals is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A logical is expected, the presence of the keyword alone is taken as true')")
              ELSE
                 WRITE(unit_nr,"(i6,'  logicals are expected')") keyword%n_var
              END IF
              WRITE(unit_nr,"('  (T,TRUE,YES,ON) and (F,FALSE,NO,OFF) are synonyms')")
           END IF
           IF (ASSOCIATED(keyword%default_value%l_val)) THEN
              WRITE(unit_nr,"('default_value : ')",advance="NO")
              WRITE(unit_nr,"(t1,' ')") keyword%default_value%l_val
              WRITE(unit_nr,"()")
           END IF
        CASE (integer_t)
           IF (level>3) THEN
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of integers is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  An integer is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' integers are expected')") keyword%n_var
              END IF
           END IF
           IF (ASSOCIATED(keyword%default_value%i_val)) THEN
              WRITE(unit_nr,"('default_value : ')",advance="NO")
              WRITE(unit_nr,"(i6,' ')") keyword%default_value%i_val
              WRITE(unit_nr,"()")
           END IF
        CASE (real_t)
           IF (level>3) THEN
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of reals is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A real is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' reals are expected')") keyword%n_var
              END IF
           END IF
           IF (ASSOCIATED(keyword%default_value%r_val)) THEN
              WRITE(unit_nr,"('default_value : ')",advance="NO")
              WRITE(unit_nr,"(e8.3,' ')") keyword%default_value%r_val
              WRITE(unit_nr,"()")
           END IF
        CASE (char_t)
           IF (level>3) THEN
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of words is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A word is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' words are expected')") keyword%n_var
              END IF
           END IF
           IF (ASSOCIATED(keyword%default_value%c_val)) THEN
              WRITE(unit_nr,"('default_value : ')",advance="NO")
              WRITE(unit_nr,"(a,' ')",advance="NO") trim(keyword%default_value%c_val)
              WRITE(unit_nr,"()")
           END IF
        CASE (char_integer_t)
           IF (level>3) THEN
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A word followed by a list of integers is expected')")
              ELSE
                 CPAssert(keyword%n_var/=1,cp_failure_level,routineP,error,failure)
                 WRITE(unit_nr,"('  A word followed by ',i6,' integers are expected')")&
                      keyword%n_var-1
              END IF
           END IF
           IF (ASSOCIATED(keyword%default_value%c_val)) THEN
              WRITE(unit_nr,"('default_value : ',a)",advance="NO")&
                   trim(keyword%default_value%c_val)
              IF (ASSOCIATED(keyword%default_value%i_val)) THEN
                 WRITE(unit_nr,"(i6,' ')",advance="NO") keyword%default_value%i_val
              ELSE
                 CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
              END IF
              WRITE(unit_nr,"()")
           END IF
        CASE (no_parse_t)
        CASE default
           WRITE(unit_nr,"('  Unknown type_of_var')")
           CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)           
        END SELECT
        IF (keyword%required) THEN
           WRITE(unit_nr,"('  This keyword is required')",advance="NO")
        ELSE
           WRITE(unit_nr,"('  This keyword is optional')",advance="NO")
        END IF
        IF (keyword%repeats) THEN
           WRITE(unit_nr,"(' and it can be repeated more than once')",advance="NO")
        END IF
        WRITE(unit_nr,"()")
        IF (SIZE(keyword%names)>1) THEN
           WRITE(unit_nr,"(a)",advance="NO") "variants    : "
           DO i=2,SIZE(keyword%names)
              WRITE(unit_nr,"(a,' ')",advance="NO") keyword%names(i)
           END DO
           WRITE(unit_nr,"()")
        ENDIF
     END IF
  END IF
END SUBROUTINE keyword_describe
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
RECURSIVE SUBROUTINE section_describe(section, unit_nr, level, error)
  TYPE(section_type), pointer :: section
  INTEGER, INTENT(in) :: unit_nr, level
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='section_desc',&
        routineP=moduleN//':'//routineN
  INTEGER :: isub,ikeyword

  failure=.false.
  IF (ASSOCIATED(section)) THEN
     CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
     
     WRITE(unit_nr,"('*** section &',a,' ***')")trim(section%name)
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
     IF (ASSOCIATED(section%subsections)) THEN
        IF (SIZE(section%subsections)>0) THEN
           WRITE(unit_nr,"('** subsections **')")
           DO isub=1,SIZE(section%subsections)
              CALL section_describe(section%subsections(isub),unit_nr,&
                   level,error=error)
           END DO
        END IF
     END IF
     WRITE(unit_nr,"('*** &end section ',a,' ***')")trim(section%name)
  ELSE
     WRITE(unit_nr,"('<section *null*>')")
  END IF
END SUBROUTINE section_describe
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
     IF (ASSOCIATED(section%subsections)) THEN
        DO isub=1,SIZE(section%subsections)
           CPInvariant(ASSOCIATED(section%subsections(isub)%section%name),cp_failure_level,routineP,error,failure)
           IF (section%subsections(isub)%section%name==upc_name) THEN
              res=isub
              EXIT
           END IF
        END DO
     END IF
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
  integer :: ik
  CHARACTER(len=default_string_length) :: upc_name

  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(section%keywords),cp_failure_level,routineP,error,failure)
  res=-2
  IF (.NOT. failure) THEN
     upc_name=keyword_type
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
           CPInvariant(ASSOCIATED(section%n_keywords(ik)%keyword),cp_failure_level,routineP,error,failure)
           DO in=1,SIZE(section%n_keywords(ik)%keyword%names)
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
SUBROUTINE section_get(section,parsed, required, repeats,id_nr,ref_count, &
     name,description,error)
  TYPE(section_type), pointer :: section
  LOGICAL, INTENT(out),OPTIONAL :: parsed, required, repeats
  INTEGER, INTENT(out),OPTIONAL :: id_nr,ref_count, n_keywords,n_subsections
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
     if (present(parsed)) parsed=section%parsed
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
  CPPrecondition(.NOT.section%parsed,cp_failure_level,routineP,error,failure)
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

        IF (SIZE(section%keywords)==section%n_keywords) THEN
           ALLOCATE(new_keywords(section%n_keywords+10),stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           DO i=-1,section%n_keywords
              new_keywords(i)%keyword => section%keywords(i)%keyword
           END DO
           DO i=section%n_keywords+1,SIZE(new_keywords)
              nullify(new_keywords(i)%keyword)
           END DO
           DEALLOCATE(section%keywords,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           section%keywords => new_keywords
        END IF
        section%n_keywords=section%n_keywords+1
        new_keywords(section%n_keywords)%keyword => keyword
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
  TYPE(subsection_p_type), DIMENSION(:), pointer :: new_subsections
  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(subsection),cp_failure_level,routineP,error,failure)
  CPPrecondition(subsection%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (.NOT.ASSOCIATED(section%subsections)) THEN
        ALLOCATE(section%subsections(1),stat=stat)
     ELSE
        DO i=1,SIZE(section%subsections)
           CALL cp_assert(subsection%name/=section%subsections(i)%section%name,&
                      cp_failure_level,cp_assertion_failed,routineP,&
                      "trying to add a subsection with a name ("//&
                      TRIM(subsection%name)//") that was already used in section "&
                      //TRIM(section%name),error,failure)
        END DO
        ALLOCATE(new_subsections(SIZE(section%subsections)+1),stat=stat)
        DO i=1,SIZE(section%subsections)
           new_subsections(i)%section => section%subsections(i)%section
        END DO
        DEALLOCATE(section%subsections,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        section%subsections => new_subsections
     END IF
     CALL section_retain(subsection,error=error)
     section%subsections(size(section%subsections))%section => subsection
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
  integer :: stat

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
     section%parsed=.true.
     ALLOCATE(section_vals%values(section%n_keywords,0),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     IF (ASSOCIATED(section%subsections)) THEN
        ALLOCATE(section_vals%subs_vals(SIZE(section%subsections),1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO i=1,SIZE(section%subsections)
           CALL section_vals_create(section_vals%subs_vals(i)%section_vals,&
                section=section%subsections(i)%section,error=error)
        END DO
     ELSE
        ALLOCATE(section_vals%subs_vals(0,1),stat=stat)
     END IF
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
     CPPreconditionNoFail(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
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

  failure=.false.
  
  IF (ASSOCIATED(section_vals)) THEN
     CPPreconditionNoFail(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
     section_vals%ref_count=section_vals%ref_count-1
     IF (section_vals%ref_count==0) THEN
        CALL section_release(section_vals%section,error=error)
        DO j=1,SIZE(section_vals%values,2)
           DO i=1,SIZE(section_vals%values,1)
              IF (ASSOCIATED(section_vals%values(i,j)%kval)) THEN
                 CALL kval_dealloc(section_vals%values(i,j)%kval,error=error)
              END IF
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
     n_subs_vals_rep,section, error)
  TYPE(section_vals_type), pointer :: section_vals
  INTEGER, INTENT(out), OPTIONAL :: ref_count, id_nr, n_repetition, n_subs_vals_rep
  TYPE(section_type), pointer :: section
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
     IF (PRESENT(n_repetition)) n_repetition=SIZE(section_vals%values,2)
     IF (PRESENT(n_subs_vals_rep)) n_subs_vals_rep=SIZE(section_vals%subs_vals,2)
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
!!     - i_rep_section: index of the repetition of the root section you are
!!       interested in (defaults to 1)
!!     - can_return_null: if the results can be null (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION section_vals_get_subs_vals(section_vals,subsection_name,&
     i_rep_section,can_return_null,error) RESULT(res)
  TYPE(section_vals_type), pointer :: section_vals
  CHARACTER(len=*), intent(in) :: subsection_name
  INTEGER, INTENT(in), optional :: i_rep_section
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(section_vals_type), pointer :: res
  
  LOGICAL :: failure, my_can_return_null
  character(len=*), parameter :: routineN='section_vals_get_subs_vals',&
        routineP=moduleN//':'//routineN
  INTEGER :: isection,irep

  failure=.false.
  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
  my_can_return_null=.false.
  if (can_return_null) my_can_return_null=can_return_null
  nullify(res)
  IF (.NOT. failure) THEN
     irep=1
     if (present(i_rep_section)) irep=i_rep_section
     CPPrecondition(irep<=SIZE(section_vals%subs_vals,2),cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     isection=section_get_subsection_index(section_vals%section,subsection_name,&
          error=error)
     IF (isection>0) res => section_vals%subs_vals(isection,irep)%section_vals
  END IF
  CPPostcondition(ASSOCIATED(res).OR.my_can_return_null,cp_failure_level,routineP,error,failure)
END FUNCTION section_vals_get_subs_vals
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
  TYPE(kval_p_type), DIMENSION(:,:), pointer :: new_values

  failure=.false.

  CPPrecondition(ASSOCIATED(section_vals),cp_failure_level,routineP,error,failure)
  CPPrecondition(section_vals%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(new_values(SIZE(section_vals%values,1),SIZE(section_vals%values,2)+1),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     DO j=1,SIZE(section_vals%values,2)
        DO i=1,SIZE(section_vals%values,1)
           new_values(i,j)%kval => section_vals%values(i,j)%kval
        END DO
     END DO
     DEALLOCATE(section_vals%values,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     section_vals%values => new_values
     DO i=1,SIZE(new_values,1)
        NULLIFY(new_values(i,SIZE(new_values,2))%kval)
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
           NULLIFY(new_sps(i,SIZE(new_sps,2))%section_vals)
           CALL section_vals_create(new_sps(i,SIZE(new_sps,2))%section_vals,&
                section=section_vals%section%subsections(i)%section,error=error)
        END DO
     END IF
  END IF
END SUBROUTINE section_vals_add_values
!***************************************************************************

!!****f* cp_keywords/section_vals_kval_get *
!!
!!   NAME
!!     section_vals_kval_get
!!
!!   FUNCTION
!!     returns the requested value
!!
!!   NOTES
!!     add the possibility to get the kval itself?
!!
!!   INPUTS
!!     - section: the root section
!!     - keyword_name: the name of the keyword you want
!!     - isection: which repetition of the section you are interested in
!!       (defaults to 1)
!!     - ikval: which repetition of the keyword/kval you are interested in
!!       (defaults to 1)
!!     - n_kval: returns number of kval available
!!     - l_val,i_val,r_val,c_val: returns the logical,integer,real or
!!       character value
!!     - l_vals,i_vals,r_vals,c_vals: returns the logical,integer,real or
!!       character arrays. The kval reamins the owner of the array
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
subroutine section_vals_kval_get(section_vals,keyword_name,i_rep_section,&
     i_rep_kval,n_kval,&
     l_val,i_val,r_val,c_val,l_vals,i_vals,r_vals,c_vals,error)
  TYPE(section_type),pointer :: section
  CHARACTER(len=*), intent(in) :: keyword_name
  INTEGER, INTENT(in), optional :: isection,ikval
  type(cp_error_type), optional, intent(inout) :: error
  TYPE(section_type),pointer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='section_kval_get',&
        routineP=moduleN//':'//routineN
  integer :: ik
  CHARACTER(len=default_string_length) :: upc_name

  failure=.false.
  
  CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)
  CPPrecondition(section%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(section%keywords),cp_failure_level,routineP,error,failure)
  NULLIFY(res)
  IF (.NOT. failure) THEN
     upc_name=keyword_type
     call uppercase(upc_name)
     DO ik=-1,0
        IF (ASSOCIATED(section%keywords(ik)%keyword)) THEN
           IF (section%keywords(ik)%keyword%names(1)==upc_name) THEN
              res => section%keywords(ik)%keyword
           END IF
        END IF
     END DO
     IF (.NOT.ASSOCIATED(res)) THEN
        k_search_loop: DO ik=1,section%n_keywords
           CPInvariant(ASSOCIATED(section%n_keywords(ik)%keyword),cp_failure_level,routineP,error,failure)
           DO in=1,SIZE(section%n_keywords(ik)%keyword%names)
              IF (section%keywords(ik)%keyword%names(in)==upc_name) THEN
                 res => section%keywords(ik)%keyword
                 exit k_search_loop
              END IF
           END DO
        END DO k_search_loop
     END IF
  END IF
END FUNCTION section_get_subsection
!***************************************************************************

! return the location of the matching entry
  FUNCTION find_keyword(keywords,key)
     TYPE(keyword_type), DIMENSION(:) :: keywords
     CHARACTER(LEN=*)                 :: key
     INTEGER                          :: find_keyword,I
     find_keyword=-1
     DO I=1,SIZE(keywords)
        IF (TRIM(keywords(I)%keyword).EQ.TRIM(key))  find_keyword=I
        IF (TRIM(keywords(I)%variant1).EQ.TRIM(key)) find_keyword=I
        IF (TRIM(keywords(I)%variant2).EQ.TRIM(key)) find_keyword=I
        IF (TRIM(keywords(I)%variant3).EQ.TRIM(key)) find_keyword=I
     ENDDO
 END FUNCTION find_keyword

END MODULE  cp_keywords
