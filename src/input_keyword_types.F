!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_keyword_types [1.0] *
!!
!!   NAME
!!     input_keyword_types
!!
!!   FUNCTION
!!     represents keywords in an input
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     06.2004 created, based on Joost cp_keywords proposal [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_keyword_types
  USE cp_units,                        ONLY: cp_unit_create,&
                                             cp_unit_desc,&
                                             cp_unit_release,&
                                             cp_unit_type
  USE input_enumeration_types,         ONLY: enum_create,&
                                             enum_release,&
                                             enum_retain,&
                                             enumeration_type
  USE input_val_types,                 ONLY: char_t,&
                                             enum_t,&
                                             integer_t,&
                                             logical_t,&
                                             no_t,&
                                             real_t,&
                                             val_create,&
                                             val_release,&
                                             val_retain,&
                                             val_type,&
                                             val_write,&
                                             val_write_internal
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE output_utilities,                ONLY: print_message
  USE string_utilities,                ONLY: compress,&
                                             substitute_special_xml_tokens,&
                                             uppercase
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_keyword_types'

  INTEGER, PARAMETER, PUBLIC :: description_string_length=10*default_string_length,&
       usage_string_length=default_string_length*2
  INTEGER, SAVE, PRIVATE :: last_keyword_id=0
  
  PUBLIC :: keyword_p_type, keyword_type, keyword_create, keyword_retain,&
       keyword_release, keyword_get, keyword_describe,&
       keyword_describe_html, write_keyword_xml

!!***
!****************************************************************************

!!****s* input_keyword_types/keyword_p_type *
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
     TYPE(keyword_type), POINTER :: keyword
  END TYPE keyword_p_type
!!***
  !****************************************************************************

!!****s* input_keyword_types/keyword_type *
!!
!!   NAME
!!     keyword_type
!!
!!   FUNCTION
!!     represent a keyword in the input
!!
!!   NOTES
!!     I have expressely avoided a format string for the type of keywords:
!!     they should easily map to basic types of fortran, if you need more
!!     information use a subsection. [fawzi]
!!
!!   ATTRIBUTES
!!     - names: the names of the current keyword (at least one should be
!!       present) for example "MAXSCF"
!!     - usage: how to use it "MAXSCF 10"
!!     - description: what does it do: "MAXSCF : determines the maximum 
!!       number of steps in an SCF run"
!!     - type_of_var: the type of keyword (controls how it is parsed)
!!       it can be one of: no_parse_t,logical_t, integer_t, real_t,
!!       char_t
!!     - n_var: number of values that should be parsed (-1=unknown)
!!     - repeats: if the keyword can be present more than once in the
!!       section
!!     - required: if the keyword is required (leaving it out will give an
!!       error)
!!     - enum: enumeration that defines the mapping between integers and
!!       strings
!!     - unit: the default unit this keyword is read in (to automatically
!!       convert to the internal cp2k units during parsing)
!!     - default_value: the default value for the keyword
!!     - lone_keyword_value: value to be used in presence of the keyword
!!       without any parameter
!!
!!   AUTHOR
!!     Joost & fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE keyword_type
     INTEGER :: ref_count,id_nr
     CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER:: names
     CHARACTER(LEN=usage_string_length) :: usage
     CHARACTER(LEN=description_string_length) :: description
     INTEGER :: type_of_var,n_var
     LOGICAL :: repeats, required
     TYPE(enumeration_type), POINTER :: enum
     TYPE(cp_unit_type), POINTER :: unit
     TYPE(val_type), POINTER :: default_value
     TYPE(val_type), POINTER :: lone_keyword_value
  END TYPE keyword_type
!!***
  !****************************************************************************

CONTAINS


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
SUBROUTINE keyword_create(keyword, name, description, usage, type_of_var,&
     n_var,repeats,variants,required,default_val,&
     default_l_val, default_r_val, default_c_val, default_i_val,&
     default_l_vals, default_r_vals, default_c_vals, default_i_vals,&
     lone_keyword_val, lone_keyword_l_val, lone_keyword_r_val, lone_keyword_c_val,&
     lone_keyword_i_val, lone_keyword_l_vals, lone_keyword_r_vals, &
     lone_keyword_c_vals, lone_keyword_i_vals, enum_c_vals, enum_i_vals, &
     enum,enum_strict,enum_desc,unit_str,error)
    TYPE(keyword_type), POINTER              :: keyword
    CHARACTER(len=*), INTENT(in)             :: name, description
    CHARACTER(len=*), INTENT(in), optional   :: usage
    INTEGER, INTENT(in), OPTIONAL            :: type_of_var, n_var
    LOGICAL, INTENT(in), OPTIONAL            :: repeats
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: variants
    LOGICAL, INTENT(in), OPTIONAL            :: required
    TYPE(val_type), OPTIONAL, POINTER        :: default_val
    LOGICAL, INTENT(in), OPTIONAL            :: default_l_val
    REAL(KIND=DP), INTENT(in), OPTIONAL      :: default_r_val
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: default_c_val
    INTEGER, INTENT(in), OPTIONAL            :: default_i_val
    LOGICAL, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: default_l_vals
    REAL(KIND=DP), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: default_r_vals
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: default_c_vals
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: default_i_vals
    TYPE(val_type), OPTIONAL, POINTER        :: lone_keyword_val
    LOGICAL, INTENT(in), OPTIONAL            :: lone_keyword_l_val
    REAL(KIND=DP), INTENT(in), OPTIONAL      :: lone_keyword_r_val
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: lone_keyword_c_val
    INTEGER, INTENT(in), OPTIONAL            :: lone_keyword_i_val
    LOGICAL, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: lone_keyword_l_vals
    REAL(KIND=DP), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: lone_keyword_r_vals
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: lone_keyword_c_vals
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: lone_keyword_i_vals
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: enum_c_vals,enum_desc
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: enum_i_vals
    TYPE(enumeration_type), OPTIONAL, &
      POINTER                                :: enum
    LOGICAL, INTENT(in), OPTIONAL            :: enum_strict
    CHARACTER(len=*), intent(in), OPTIONAL   :: unit_str
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keyword_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: check, failure

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  ALLOCATE(keyword,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     keyword%ref_count=1
     last_keyword_id=last_keyword_id+1
     keyword%id_nr=last_keyword_id
     nullify(keyword%unit)

     IF (PRESENT(variants)) THEN
        ALLOCATE(keyword%names(SIZE(variants)+1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           keyword%names(1)=name
           DO i=1,SIZE(variants)
              keyword%names(i+1)=variants(i)
           END DO
        END IF
     ELSE
        ALLOCATE(keyword%names(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.not.failure) keyword%names(1)=name
     END IF
     DO i=1,SIZE(keyword%names)
        CALL uppercase(keyword%names(i))
     END DO
     
     IF (PRESENT(usage)) THEN
        CPPrecondition(LEN_TRIM(usage)<=LEN(keyword%usage),cp_failure_level,routineP,error,failure)
        keyword%usage=usage
     ELSE
        keyword%usage=""
     END IF
     CPPrecondition(LEN_TRIM(description)<=LEN(keyword%description),cp_failure_level,routineP,error,failure)
     keyword%description=description
     
     keyword%repeats=.FALSE.
     IF (PRESENT(repeats)) keyword%repeats=repeats
     keyword%required=.FALSE.
     IF (PRESENT(required)) keyword%required=required

     NULLIFY(keyword%enum)
     IF (PRESENT(enum)) THEN
        keyword%enum => enum
        IF (ASSOCIATED(enum)) CALL enum_retain(enum,error=error)
     END IF
     IF (PRESENT(enum_i_vals)) THEN
        CPPrecondition(PRESENT(enum_c_vals),cp_failure_level,routineP,error,failure)
        CPPrecondition(.NOT.ASSOCIATED(keyword%enum),cp_failure_level,routineP,error,failure)
        CALL enum_create(keyword%enum,c_vals=enum_c_vals,i_vals=enum_i_vals,&
             desc=enum_desc,strict=enum_strict,error=error)
     ELSE
        CPPrecondition(.NOT.PRESENT(enum_c_vals),cp_failure_level,routineP,error,failure)
     END IF

     NULLIFY(keyword%default_value, keyword%lone_keyword_value)
     IF (PRESENT(default_val)) THEN
        CALL cp_assert(.NOT.(PRESENT(default_l_val).OR.PRESENT(default_l_vals).OR.&
             PRESENT(default_i_val).OR.PRESENT(default_i_vals).OR.&
             PRESENT(default_r_val).or.PRESENT(default_r_vals).OR.&
             PRESENT(default_c_val).OR.PRESENT(default_c_vals)),cp_failure_level,&
             cp_assertion_failed,routineP,&
             "you should pass either default_val or a default value, not both",&
             error,failure)
        keyword%default_value => default_val
        IF (ASSOCIATED(default_val%enum)) THEN
           IF (ASSOCIATED(keyword%enum)) THEN
              CPAssert(keyword%enum%id_nr==default_val%enum%id_nr,cp_failure_level,routineP,error,failure)
           ELSE
              keyword%enum => default_val%enum
              CALL enum_retain(keyword%enum,error=error)
           END IF
        ELSE
           CPAssert(.NOT.ASSOCIATED(keyword%enum),cp_failure_level,routineP,error,failure)
        END IF
        CALL val_retain(default_val,error=error)
     END IF
     IF (.not.ASSOCIATED(keyword%default_value)) THEN
        CALL val_create(keyword%default_value,l_val=default_l_val,&
             l_vals=default_l_vals,i_val=default_i_val,i_vals=default_i_vals,&
             r_val=default_r_val,r_vals=default_r_vals,c_val=default_c_val,&
             c_vals=default_c_vals,enum=keyword%enum,error=error)
     END IF

     keyword%type_of_var=keyword%default_value%type_of_var
     IF (keyword%default_value%type_of_var==no_t) THEN
        CALL val_release(keyword%default_value,error=error)
     END IF

     IF (keyword%type_of_var==no_t) THEN
        IF (PRESENT(type_of_var)) THEN
           keyword%type_of_var=type_of_var
        ELSE
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                "keyword "//TRIM(keyword%names(1))//&
                " assumed undefined type by default",error,failure)
        END IF
     ELSE IF (PRESENT(type_of_var)) THEN
        CALL cp_assert(keyword%type_of_var==type_of_var,cp_failure_level,&
             cp_assertion_failed,routineP, &
             "keyword "//TRIM(keyword%names(1))//&
             " has a type different from the type of the default_value",&
             error,failure)
        keyword%type_of_var=type_of_var
     END IF

     IF (keyword%type_of_var==no_t) THEN
        CALL val_create(keyword%default_value,error=error)
     END IF
  END IF

  IF (.NOT.failure) THEN
     IF (PRESENT(lone_keyword_val)) THEN
        CALL cp_assert(.NOT.(PRESENT(lone_keyword_l_val).OR.PRESENT(lone_keyword_l_vals).OR.&
             PRESENT(lone_keyword_i_val).OR.PRESENT(lone_keyword_i_vals).OR.&
             PRESENT(lone_keyword_r_val).OR.PRESENT(lone_keyword_r_vals).OR.&
             PRESENT(lone_keyword_c_val).OR.PRESENT(lone_keyword_c_vals)),&
             cp_failure_level, cp_assertion_failed,routineP,&
             "you should pass either lone_keyword_val or a lone_keyword value, not both",&
             error,failure)
        keyword%lone_keyword_value => lone_keyword_val
        CALL val_retain(lone_keyword_val,error=error)
        IF (ASSOCIATED(lone_keyword_val%enum)) THEN
           IF (ASSOCIATED(keyword%enum)) THEN
              CALL cp_assert(keyword%enum%id_nr==lone_keyword_val%enum%id_nr, &
                   cp_failure_level,cp_assertion_failed,routineP, &
                   "keyword%enum%id_nr==lone_keyword_val%enum%id_nr", &
                   error,failure)
           ELSE
              IF (ASSOCIATED(keyword%lone_keyword_value)) THEN
                CALL cp_assert(.FALSE., cp_failure_level, cp_precondition_failed, &
                               routineP, ".NOT. ASSOCIATED(keyword%lone_keyword_value)", &
                               error, failure)
              END IF
              keyword%enum => lone_keyword_val%enum
              CALL enum_retain(keyword%enum,error=error)
           END IF
        ELSE
           CPAssert(.NOT.ASSOCIATED(keyword%enum),cp_failure_level,routineP,error,failure)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(keyword%lone_keyword_value)) THEN
        CALL val_create(keyword%lone_keyword_value,l_val=lone_keyword_l_val,&
             l_vals=lone_keyword_l_vals,i_val=lone_keyword_i_val,i_vals=lone_keyword_i_vals,&
             r_val=lone_keyword_r_val,r_vals=lone_keyword_r_vals,c_val=lone_keyword_c_val,&
             c_vals=lone_keyword_c_vals,enum=keyword%enum,error=error)
     END IF
     IF (ASSOCIATED(keyword%lone_keyword_value)) THEN
        IF (keyword%lone_keyword_value%type_of_var==no_t) THEN
           CALL val_release(keyword%lone_keyword_value,error=error)
        ELSE
           CALL cp_assert(keyword%lone_keyword_value%type_of_var==keyword%type_of_var,&
                cp_failure_level,cp_assertion_failed,routineP,&
                "lone_keyword_value type incompatible with "//&
                "keyword type",error,failure)
           IF (keyword%type_of_var==enum_t) THEN
              IF (keyword%enum%strict) THEN
                 check = .FALSE.
                 DO i=1, SIZE(keyword%enum%i_vals)
                    check = check .OR. (keyword%default_value%i_val(1)==keyword%enum%i_vals(i))
                 END DO
                 CALL cp_assert(check,cp_failure_level,&
                         cp_assertion_failed,routineP,"default value not in enumeration : "//keyword%names(1), &
                           error,failure) 
              ENDIF
           ENDIF
        END IF
     END IF

     keyword%n_var=1
     IF (ASSOCIATED(keyword%default_value)) THEN
        SELECT CASE(keyword%default_value%type_of_var)
        CASE(logical_t)
           keyword%n_var=SIZE(keyword%default_value%l_val)
        CASE(integer_t)
           keyword%n_var=SIZE(keyword%default_value%i_val)
        CASE(enum_t)
           IF (keyword%enum%strict) THEN
              check = .FALSE.
              DO i=1, SIZE(keyword%enum%i_vals)
                 check = check .OR. (keyword%default_value%i_val(1)==keyword%enum%i_vals(i))
              END DO
              CALL cp_assert(check,cp_failure_level,&
                      cp_assertion_failed,routineP,"default value not in enumeration : "//keyword%names(1), &
                          error,failure) 
           ENDIF
           keyword%n_var=SIZE(keyword%default_value%i_val)
        CASE(real_t)
           keyword%n_var=SIZE(keyword%default_value%r_val)
        CASE(char_t)
           keyword%n_var=SIZE(keyword%default_value%c_val)
        CASE(no_t)
           keyword%n_var=0
        CASE default
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END SELECT
     END IF
     IF (PRESENT(n_var)) keyword%n_var=n_var

     IF (PRESENT(unit_str)) THEN
        CALL cp_unit_create(keyword%unit,unit_str,error=error)
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
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keyword_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(keyword%ref_count>0,cp_failure_level,routineP,error)
     keyword%ref_count=keyword%ref_count+1
  END IF
END SUBROUTINE keyword_retain
!***************************************************************************

!!****f* input_keyword_types/keyword_release *
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
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keyword_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(keyword)) THEN
     CPPreconditionNoFail(keyword%ref_count>0,cp_failure_level,routineP,error)
     keyword%ref_count=keyword%ref_count-1
     IF (keyword%ref_count==0) THEN
        DEALLOCATE(keyword%names,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CALL val_release(keyword%default_value,error=error)
        CALL val_release(keyword%lone_keyword_value,error=error)
        CALL enum_release(keyword%enum,error=error)
        CALL cp_unit_release(keyword%unit,error=error)
        DEALLOCATE(keyword,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(keyword)
END SUBROUTINE keyword_release
!***************************************************************************

!!****f* input_keyword_types/keyword_get *
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
     default_value, lone_keyword_value,required,repeats,enum,error)
    TYPE(keyword_type), POINTER              :: keyword
    CHARACTER(len=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: names
    CHARACTER(len=*), INTENT(out), OPTIONAL  :: usage, description
    INTEGER, INTENT(out), OPTIONAL           :: type_of_var, n_var
    TYPE(val_type), OPTIONAL, POINTER        :: default_value, &
                                                lone_keyword_value
    LOGICAL, INTENT(out), OPTIONAL           :: required, repeats
    TYPE(enumeration_type), POINTER, OPTIONAL:: enum
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keyword_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(names)) names => keyword%names
     IF (PRESENT(usage)) usage=keyword%usage
     IF (PRESENT(description)) description=keyword%description
     IF (PRESENT(type_of_var)) type_of_var=keyword%type_of_var
     IF (PRESENT(n_var)) n_var=keyword%n_var
     IF (PRESENT(repeats)) repeats=keyword%repeats
     IF (PRESENT(required)) required=keyword%required
     IF (PRESENT(default_value)) default_value => keyword%default_value
     IF (PRESENT(lone_keyword_value)) lone_keyword_value => keyword%lone_keyword_value
     IF (PRESENT(enum)) enum => keyword%enum
  END IF
END SUBROUTINE keyword_get
!***************************************************************************

!!****f* input_keyword_types/keyword_describe [1.0] *
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
    TYPE(keyword_type), POINTER              :: keyword
    INTEGER, INTENT(in)                      :: unit_nr, level
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keyword_describe', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_string
    INTEGER                                  :: i, l
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure.AND.level>0) THEN
     WRITE(unit_nr,"(a,a,a)") "                           ---",&
          TRIM(keyword%names(1)),"---"
     IF (level>1) THEN
        WRITE(unit_nr,"(a,a)") "usage         : ",TRIM(keyword%usage)
     END IF
     IF (level>2) THEN
        WRITE(unit_nr,"(a)")"description   : "
        CALL print_message(TRIM(keyword%description),unit_nr,0,0,0)
        IF (level>3) THEN
           SELECT CASE(keyword%type_of_var)
           CASE (logical_t)
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of logicals is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A logical is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,'  logicals are expected')") keyword%n_var
              END IF
              WRITE(unit_nr,"('  (T,TRUE,YES,ON) and (F,FALSE,NO,OFF) are synonyms')")
           CASE (integer_t)
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of integers is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  An integer is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' integers are expected')") keyword%n_var
              END IF
           CASE (real_t)
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of reals is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A real is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' reals are expected')") keyword%n_var
              END IF
              IF (ASSOCIATED(keyword%unit)) THEN
                 c_string=cp_unit_desc(keyword%unit,accept_undefined=.TRUE.,&
                      error=error)
                 WRITE(unit_nr,"('the default unit of measure is ',a)")&
                      TRIM(c_string)
              END IF
           CASE (char_t)
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of words is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A word is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' words are expected')") keyword%n_var
              END IF
           CASE (enum_t)
              IF (keyword%n_var==-1) THEN
                 WRITE(unit_nr,"('  A list of keywords is expected')")
              ELSE IF (keyword%n_var==1) THEN
                 WRITE(unit_nr,"('  A keyword is expected')")
              ELSE
                 WRITE(unit_nr,"(i6,' keywords are expected')") keyword%n_var
              END IF
           CASE (no_t)
              WRITE(unit_nr,"('  Non-standard type.')")
           CASE default
              CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)           
           END SELECT
        END IF
        IF (keyword%type_of_var==enum_t) THEN
           IF (level>3) THEN
              WRITE(unit_nr,"('  valid keywords:')")
              DO i=1,SIZE(keyword%enum%c_vals)
                 c_string=keyword%enum%c_vals(i)
                 IF (LEN_TRIM(keyword%enum%desc(i))>0) THEN
                    WRITE (unit_nr,"('  - ',a,' : ',a,'.')") &
                         TRIM(c_string),TRIM(keyword%enum%desc(i))
                 ELSE
                    WRITE (unit_nr,"('  - ',a)") TRIM(c_string)
                 END IF
              END DO
           ELSE
              WRITE(unit_nr,"('  valid keywords:')",advance='NO')
              l=17
              DO i=1,SIZE(keyword%enum%c_vals)
                 c_string=keyword%enum%c_vals(i)
                 IF (l+LEN_TRIM(c_string)>72.AND.l>14) THEN
                    WRITE (unit_nr,"(/,'    ')",advance='NO')
                    l=4
                 END IF
                 WRITE (unit_nr,"(' ',a)",advance='NO') TRIM(c_string)
                 l=LEN_TRIM(c_string)+3
              END DO
              WRITE (unit_nr,"()")
           END IF
           IF (.NOT.keyword%enum%strict) THEN
              WRITE (unit_nr,"('     other integer values are also accepted.')")
           END IF
        END IF
        IF (ASSOCIATED(keyword%default_value).AND.keyword%type_of_var/=no_t) THEN
           WRITE(unit_nr,"('default_value : ')",advance="NO")
           CALL val_write(keyword%default_value,unit_nr=unit_nr,error=error)
        END IF
        IF (ASSOCIATED(keyword%lone_keyword_value).AND.keyword%type_of_var/=no_t) THEN
           WRITE(unit_nr,"('lone_keyword  : ')",advance="NO")
           CALL val_write(keyword%lone_keyword_value,unit_nr=unit_nr,error=error)
        END IF
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

!!****f* input_keyword_types/keyword_describe_html [1.0] *
!!
!!   NAME
!!     keyword_describe_html
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Joost VandeVondele [10.2004], based on keyword_describe
!!
!!*** **********************************************************************
SUBROUTINE keyword_describe_html(keyword, unit_nr, error)
    TYPE(keyword_type), POINTER              :: keyword
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keyword_describe_thml', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_string, my_unit
    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
  CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
      WRITE(unit_nr,'(a)') '<TR><TD WITDH="20%">'// &
                           '<A NAME="'//TRIM(keyword%names(1))//'"><u>'// TRIM(keyword%names(1))//'</u></A>'// &
                            '<TD WIDTH="80%">'
      WRITE(unit_nr,'(a)') '<TR><TD WIDTH="20%"><TD WIDTH="80%">'//TRIM(keyword%usage)
      WRITE(unit_nr,'(a)') '<TR><TD WIDTH="20%"><TD WIDTH="80%"><i>'//TRIM(keyword%description)//'</i>'
      IF (keyword%required) THEN
            WRITE(unit_nr,"(a)",advance="NO") '<TR><TD WIDTH="10%"><TD> This required keyword '
      ELSE
            WRITE(unit_nr,"(a)",advance="NO") '<TR><TD WIDTH="10%"><TD> This optional keyword '
      END IF
      SELECT CASE(keyword%type_of_var)
      CASE (logical_t)
               IF (keyword%n_var==-1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a list of logicals'
               ELSE IF (keyword%n_var==1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a logical'
               ELSE
                  WRITE(unit_nr,'(a,i6,a)',ADVANCE="NO") 'expects precisely',keyword%n_var,' logicals'
               END IF
               ! (provide a link to this info) WRITE(unit_nr,"('  (T,TRUE,YES,ON) and (F,FALSE,NO,OFF) are synonyms')")
      CASE (integer_t)
               IF (keyword%n_var==-1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a list of integers'
               ELSE IF (keyword%n_var==1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects an integer'
               ELSE
                  WRITE(unit_nr,'(a,i6,a)',ADVANCE="NO") 'expects precisely',keyword%n_var,' integers'
               END IF
      CASE (real_t)
               IF (keyword%n_var==-1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a list of reals'
               ELSE IF (keyword%n_var==1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a real'
               ELSE
                  WRITE(unit_nr,'(a,i6,a)',ADVANCE="NO") 'expects precisely',keyword%n_var,' reals'
               END IF
      CASE (char_t)
               IF (keyword%n_var==-1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a list of words'
               ELSE IF (keyword%n_var==1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a word'
               ELSE
                  WRITE(unit_nr,'(a,i6,a)',ADVANCE="NO") 'expects precisely',keyword%n_var,' words'
               END IF
      CASE (enum_t)
               IF (keyword%n_var==-1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a list of keywords'
               ELSE IF (keyword%n_var==1) THEN
                  WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a keyword'
               ELSE
                  WRITE(unit_nr,'(a,i6,a)',ADVANCE="NO") 'expects precisely',keyword%n_var,' keywords'
               END IF
      CASE (no_t)
               WRITE(unit_nr,'(a)',ADVANCE="NO") 'expects a non-standard input type'
      CASE DEFAULT
               CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)           
      END SELECT
      IF (keyword%repeats) THEN
            WRITE(unit_nr,"(', and may repeat')",ADVANCE="NO")
      END IF
      IF (ASSOCIATED(keyword%lone_keyword_value).AND.keyword%type_of_var/=no_t) THEN
            WRITE(unit_nr,'(a)',advance="NO") '<TR><TD WIDTH="10%"><TD>This keyword behaves as a switch'
            CALL val_write(keyword%lone_keyword_value,unit_nr=unit_nr,error=error)
      END IF
      IF (ASSOCIATED(keyword%default_value).AND.keyword%type_of_var/=no_t) THEN
            WRITE(unit_nr,'(a)',advance="NO") '. Default  : '
            CALL val_write(keyword%default_value,unit_nr=unit_nr,error=error)
      END IF
      IF (ASSOCIATED(keyword%unit)) THEN
         my_unit=cp_unit_desc(keyword%unit,accept_undefined=.TRUE.,error=error)
         WRITE(unit_nr,'(a)',advance="NO") '. Default unit of measure : '//TRIM(my_unit)
      END IF

      IF (keyword%type_of_var==enum_t) THEN
            IF (ANY(keyword%enum%desc.NE."")) THEN
                WRITE(unit_nr,'(a)',advance='NO') '<TR><TD WIDTH="10%"><TD>  valid keywords:<ul> '
                DO i=1,SIZE(keyword%enum%c_vals)
                   c_string=keyword%enum%c_vals(i)
                   WRITE (unit_nr,'(a,a,a,a,a)',advance='NO') "<li><code>",&
                        TRIM(c_string), "</code>: ",TRIM(keyword%enum%desc(i)),"</li>"
                END DO
                WRITE (unit_nr,'(a)',advance='NO') '</ul>'
            ELSE
                WRITE(unit_nr,'(a)',advance='NO') '<TR><TD WIDTH="10%"><TD>  valid keywords: '
                IF (SIZE(keyword%enum%c_vals)>0) &
                    WRITE (unit_nr,'(a)',advance='NO') TRIM(keyword%enum%c_vals(1))
                DO i=2,SIZE(keyword%enum%c_vals)
                   c_string=keyword%enum%c_vals(i)
                   WRITE (unit_nr,'(a)',advance='NO') ", "//TRIM(c_string)
                END DO
            ENDIF
            IF (.NOT.keyword%enum%strict) THEN
               WRITE (unit_nr,"('     other integer values are also accepted.')")
            END IF
      END IF
      IF (SIZE(keyword%names)>1) THEN
          WRITE(unit_nr,"(a)",advance="NO") '<TR><TD WIDTH="10%"><TD>variants: '
          DO i=2,SIZE(keyword%names)
             WRITE(unit_nr,"(a)",advance="NO") TRIM(keyword%names(i))
             IF (i .NE. SIZE(keyword%names)) THEN
                 WRITE (unit_nr,'(a)',advance='NO') ', '
             ENDIF
          END DO
      ENDIF
    ENDIF
END SUBROUTINE keyword_describe_html

! *****************************************************************************

  SUBROUTINE write_keyword_xml(keyword,level,unit_number,error)

    TYPE(keyword_type), POINTER                  :: keyword
    INTEGER, INTENT(IN)                          :: level,unit_number
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error

    ! Local parameters

    CHARACTER(LEN=*), PARAMETER :: routineN = "write_keyword_xml",&
                                   routineP = moduleN//':'//routineN

    ! Local variables

    CHARACTER(LEN=1000) :: string
    CHARACTER(LEN=5)    :: keyword_n_var
    CHARACTER(LEN=3)    :: repeats,required
    REAL(KIND=dp)       :: value
    INTEGER             :: i,l0,l1,l2,l3,l4
    LOGICAL             :: failure

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(keyword),cp_failure_level,routineP,error,failure)
    CPPrecondition(keyword%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN

!     *** Indentation for current level, next level, etc. ***

      l0 = level
      l1 = level + 1
      l2 = level + 2
      l3 = level + 3
      l4 = level + 4

      IF (keyword%required) THEN
        required = "yes"
      ELSE
        required = "no "
      END IF

      IF (keyword%repeats) THEN
        repeats = "yes"
      ELSE
        repeats = "no "
      END IF

!     *** Write keyword element ***

      WRITE (UNIT=unit_number,FMT="(A)")&
        REPEAT(" ",l0)//"<KEYWORD required="""//TRIM(required)//&
        """ repeats="""//TRIM(repeats)//""">",&
        REPEAT(" ",l1)//"<NAME type=""default"">"//&
        TRIM(keyword%names(1))//"</NAME>"

      DO i=2,SIZE(keyword%names)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<NAME type=""alias"">"//&
          TRIM(keyword%names(i))//"</NAME>"
      END DO

      SELECT CASE(keyword%type_of_var)
      CASE (logical_t)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DATA_TYPE kind=""BOOLEAN"">"
      CASE (integer_t)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DATA_TYPE kind=""INTEGER"">"
      CASE (real_t)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DATA_TYPE kind=""FLOAT"">"
      CASE (char_t)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DATA_TYPE kind=""STRING"">"
      CASE (enum_t)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DATA_TYPE kind=""ENUMERATION"">"
        IF (keyword%enum%strict) THEN
          WRITE (UNIT=unit_number,FMT="(A)")&
            REPEAT(" ",l2)//"<ENUMERATION strict=""yes"">"
        ELSE
          WRITE (UNIT=unit_number,FMT="(A)")&
            REPEAT(" ",l2)//"<ENUMERATION strict=""no"">"
        END IF
        DO i=1,SIZE(keyword%enum%c_vals)
          CALL substitute_special_xml_tokens(keyword%enum%desc(i),string,&
                                             ltu=.FALSE.)
          WRITE (UNIT=unit_number,FMT="(A)")&
            REPEAT(" ",l3)//"<ITEM>",&
            REPEAT(" ",l4)//"<NAME>"//&
            TRIM(ADJUSTL(keyword%enum%c_vals(i)))//"</NAME>",&
            REPEAT(" ",l4)//"<DESCRIPTION>"//&
            TRIM(ADJUSTL(string))//"</DESCRIPTION>",&
            REPEAT(" ",l3)//"</ITEM>"
        END DO
        WRITE (UNIT=unit_number,FMT="(A)") REPEAT(" ",l2)//"</ENUMERATION>"
!MK   CASE (no_t)
!MK     no standard type, no action?
      CASE DEFAULT
        CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
      END SELECT

      keyword_n_var = ""
      WRITE (UNIT=keyword_n_var,FMT="(I5)") keyword%n_var
      WRITE (UNIT=unit_number,FMT="(A)")&
        REPEAT(" ",l2)//"<N_VAR>"//TRIM(ADJUSTL(keyword_n_var))//"</N_VAR>",&
        REPEAT(" ",l1)//"</DATA_TYPE>"

      CALL substitute_special_xml_tokens(keyword%usage,string,ltu=.TRUE.)
      WRITE (UNIT=unit_number,FMT="(A)")&
        REPEAT(" ",l1)//"<USAGE>"//TRIM(string)//"</USAGE>"

      CALL substitute_special_xml_tokens(keyword%description,string,ltu=.FALSE.)
      WRITE (UNIT=unit_number,FMT="(A)")&
        REPEAT(" ",l1)//"<DESCRIPTION>"//TRIM(string)//"</DESCRIPTION>"

      IF (ASSOCIATED(keyword%default_value).AND.&
          (keyword%type_of_var /= no_t)) THEN
        IF (ASSOCIATED(keyword%unit)) THEN
          CALL val_write_internal(val=keyword%default_value,&
                                  string=string,&
                                  unit=keyword%unit,&
                                  error=error)
        ELSE
          CALL val_write_internal(val=keyword%default_value,&
                                  string=string,&
                                  error=error)
        END IF
        CALL compress(string)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DEFAULT_VALUE>"//&
          TRIM(ADJUSTL(string))//"</DEFAULT_VALUE>"
      END IF

      IF (ASSOCIATED(keyword%unit)) THEN
        string = cp_unit_desc(keyword%unit,accept_undefined=.TRUE.,error=error)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<DEFAULT_UNIT>"//&
          TRIM(ADJUSTL(string))//"</DEFAULT_UNIT>"
      END IF

      IF (ASSOCIATED(keyword%lone_keyword_value).AND.&
          (keyword%type_of_var /= no_t)) THEN
        CALL val_write_internal(val=keyword%lone_keyword_value,&
                                string=string,&
                                error=error)
        WRITE (UNIT=unit_number,FMT="(A)")&
          REPEAT(" ",l1)//"<LONE_KEYWORD_VALUE>"//&
          TRIM(ADJUSTL(string))//"</LONE_KEYWORD_VALUE>"
      END IF

      WRITE (UNIT=unit_number,FMT="(A)") REPEAT(" ",l0)//"</KEYWORD>"

    END IF

  END SUBROUTINE write_keyword_xml

! *****************************************************************************

END MODULE input_keyword_types
