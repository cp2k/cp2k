!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/input_enumeration_types [1.0] *
!!
!!   NAME
!!     input_enumeration_types
!!
!!   FUNCTION
!!     represents an enumeration, i.e. a mapping between integers and strings
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE input_enumeration_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE string_utilities,                ONLY: uppercase

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_enumeration_types'
  INTEGER, SAVE, PRIVATE :: last_enumeration_id=0

  INTEGER, PARAMETER, PUBLIC :: keyword_desc_length=default_string_length
  PUBLIC :: enumeration_type
  PUBLIC :: enum_create, enum_retain, enum_release, enum_i2c, enum_c2i
!!***
!****************************************************************************

!!****s* input_enumeration_types/enumeration_type [1.0] *
!!
!!   NAME
!!     enumeration_type
!!
!!   FUNCTION
!!     represents an enumaration, i.e. a mapping between strings and numbers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number (unique)
!!     - ref_count: reference count
!!     - c_vals: string values
!!     - i_vals: integer values
!!     - strict: if integer values not in the list should be accepted
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE enumeration_type
     INTEGER :: id_nr, ref_count
     CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: c_vals
     CHARACTER(len=keyword_desc_length), DIMENSION(:), POINTER :: desc
     INTEGER, DIMENSION(:), POINTER :: i_vals
     LOGICAL :: strict
  END TYPE enumeration_type
!!***
  !****************************************************************************

CONTAINS

!!****f* input_enumeration_types/enum_create [1.0] *
!!
!!   NAME
!!     enum_create
!!
!!   FUNCTION
!!     creates an enumeration
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - enum: the enumeration to be created
!!     - c_vals: string values
!!     - i_vals: integer values
!!     - strict: if integer values not in the list should be accepted, 
!!       defaults defaults to true
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE enum_create(enum,c_vals,i_vals,desc,strict,error)
    TYPE(enumeration_type), POINTER          :: enum
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: c_vals
    INTEGER, DIMENSION(:), INTENT(in)        :: i_vals
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), optional                   :: desc
    LOGICAL, INTENT(in), OPTIONAL            :: strict
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat,i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(enum),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(c_vals)==SIZE(i_vals),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     ALLOCATE(enum,stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     last_enumeration_id=last_enumeration_id+1
     enum%id_nr=last_enumeration_id
     enum%ref_count=1
     ALLOCATE(enum%c_vals(SIZE(c_vals)),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     DO i=1,SIZE(enum%c_vals)
        enum%c_vals(i)=c_vals(i)
        call uppercase(enum%c_vals(i))
     END DO
     ALLOCATE(enum%i_vals(SIZE(i_vals)),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     enum%i_vals=i_vals
     enum%strict=.TRUE.
     IF (PRESENT(strict)) enum%strict=strict
     ALLOCATE(enum%desc(SIZE(c_vals)),stat=stat)
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     IF (PRESENT(desc)) THEN
        CPPrecondition(SIZE(enum%desc)==SIZE(desc),cp_failure_level,routineP,error,failure)
        DO i=1,SIZE(enum%desc)
           CPPrecondition(LEN_TRIM(desc(i))<=LEN(enum%desc(i)),cp_failure_level,routineP,error,failure)
           enum%desc(i)=desc(i)
        END DO
     ELSE
        DO i=1,SIZE(enum%desc)
           enum%desc(i)=' '
        END DO
     END IF
  END IF
END SUBROUTINE enum_create
!***************************************************************************

!!****f* input_enumeration_types/enum_retain [1.0] *
!!
!!   NAME
!!     enum_retain
!!
!!   FUNCTION
!!     retains the given enumeration
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - enum: the obect to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE enum_retain(enum,error)
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(enum),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
       CPPreconditionNoFail(enum%ref_count>0,cp_failure_level,routineP,error)
       enum%ref_count=enum%ref_count+1
    END IF
  END SUBROUTINE enum_retain
!***************************************************************************

!!****f* input_enumeration_types/enum_release [1.0] *
!!
!!   NAME
!!     enum_release
!!
!!   FUNCTION
!!     releases the given enumeration
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - enum: the obect to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE enum_release(enum,error)
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(enum)) THEN
       CPPreconditionNoFail(enum%ref_count>0,cp_failure_level,routineP,error)
       enum%ref_count=enum%ref_count-1
       IF (enum%ref_count==0) THEN
          DEALLOCATE(enum%c_vals,stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          DEALLOCATE(enum%i_vals,stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          DEALLOCATE(enum%desc,stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          DEALLOCATE(enum,stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(enum)
  END SUBROUTINE enum_release
!***************************************************************************

!!****f* input_enumeration_types/enum_i2c [1.0] *
!!
!!   NAME
!!     enum_i2c
!!
!!   FUNCTION
!!     maps an integer to a string
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - enum: the enumeration to use for the mapping
!!     - i: the value to map
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION enum_i2c(enum,i,error) RESULT(res)
    TYPE(enumeration_type), POINTER          :: enum
    INTEGER, INTENT(in)                      :: i
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    CHARACTER(len=default_string_length)     :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_i2c', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j
    LOGICAL                                  :: failure, found

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(enum),cp_failure_level,routineP,error,failure)
  CPPrecondition(enum%ref_count>0,cp_failure_level,routineP,error,failure)
  res=" "
  IF (.NOT. failure) THEN
     found=.FALSE.
     DO j=1,SIZE(enum%i_vals)
        IF (enum%i_vals(j)==i) THEN
           res=enum%c_vals(j)
           found=.TRUE.
           EXIT
        END IF
     END DO
     IF (.NOT.found) THEN
        if (enum%strict) then
         DO j=1,SIZE(enum%desc)
          print *, trim(enum%desc(j))
          print *, trim(enum%c_vals(j))
         ENDDO
          print *, enum%i_vals
        end if
        CALL cp_assert(.NOT.enum%strict,cp_failure_level,cp_assertion_failed,&
             routineP,"invalid value for enumeration:"//cp_to_string(i),&
             error,failure)
        res=ADJUSTL(cp_to_string(i))
     END IF
  END IF
END FUNCTION enum_i2c
!***************************************************************************

!!****f* input_enumeration_types/enum_c2i [1.0] *
!!
!!   NAME
!!     enum_c2i
!!
!!   FUNCTION
!!     maps a string to an integer
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - enum: the enumeration to use for the mapping
!!     - c: the value to map
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION enum_c2i(enum,c,error) RESULT(res)
    TYPE(enumeration_type), POINTER          :: enum
    CHARACTER(len=*), INTENT(in)             :: c
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_c2i', &
      routineP = moduleN//':'//routineN
    character(len=default_string_length) :: upc
    INTEGER                                  :: iostat, j
    LOGICAL                                  :: failure, found

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(enum),cp_failure_level,routineP,error,failure)
  CPPrecondition(enum%ref_count>0,cp_failure_level,routineP,error,failure)
  upc=c
  call uppercase(upc)
  IF (.NOT. failure) THEN
     found=.FALSE.
     DO j=1,SIZE(enum%c_vals)
        IF (enum%c_vals(j)==upc) THEN
           res=enum%i_vals(j)
           found=.TRUE.
           EXIT
        END IF
     END DO
     IF (.NOT.found) THEN
        CALL cp_assert(.NOT.enum%strict,cp_failure_level,cp_assertion_failed,&
             routineP,"invalid value for enumeration:"//TRIM(c),&
             error,failure)
        READ(c,"(i10)",iostat=iostat) res
        CALL cp_assert(iostat==0,cp_failure_level,cp_assertion_failed,&
             routineP,"invalid value for enumeration2:"//TRIM(c),&
             error,failure)
     END IF
  END IF
END FUNCTION enum_c2i
!***************************************************************************

END MODULE input_enumeration_types
