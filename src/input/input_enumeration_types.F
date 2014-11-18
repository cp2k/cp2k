!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief represents an enumeration, i.e. a mapping between integers and strings
!> \par History
!>      08.2004 created [fawzi]
!> \author fawzi
! *****************************************************************************
MODULE input_enumeration_types
  
  USE kinds,                           ONLY: default_string_length
  USE string_utilities,                ONLY: a2s,&
                                             uppercase
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_enumeration_types'
  INTEGER, SAVE, PRIVATE :: last_enumeration_id=0

  PUBLIC :: enumeration_type
  PUBLIC :: enum_create, enum_retain, enum_release, enum_i2c, enum_c2i

! *****************************************************************************
!> \brief represents an enumaration, i.e. a mapping between strings and numbers
!> \param id_nr identification number (unique)
!> \param ref_count reference count
!> \param c_vals string values
!> \param i_vals integer values
!> \param strict if integer values not in the list should be accepted
!> \author fawzi
! *****************************************************************************
  TYPE char_array
     CHARACTER, DIMENSION(:), POINTER :: chars => Null()
  END TYPE char_array

  TYPE enumeration_type
     INTEGER :: id_nr, ref_count
     CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: c_vals
     TYPE(char_array), DIMENSION(:), POINTER :: desc => Null()
     INTEGER, DIMENSION(:), POINTER :: i_vals
     LOGICAL :: strict
  END TYPE enumeration_type

CONTAINS

! *****************************************************************************
!> \brief creates an enumeration
!> \param enum the enumeration to be created
!> \param c_vals string values
!> \param i_vals integer values
!> \param desc ...
!> \param strict if integer values not in the list should be accepted,
!>        defaults defaults to true
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE enum_create(enum,c_vals,i_vals,desc,strict,error)
    TYPE(enumeration_type), POINTER          :: enum
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: c_vals
    INTEGER, DIMENSION(:), INTENT(in)        :: i_vals
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: desc
    LOGICAL, INTENT(in), OPTIONAL            :: strict
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, n, stat
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
        CALL uppercase(enum%c_vals(i))
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
           n = LEN_TRIM(desc(i))
           ALLOCATE(enum%desc(i)%chars(n))
           DO j=1, n
             enum%desc(i)%chars(j) = desc(i)(j:j)
           END DO
       END DO
     ELSE
        DO i=1,SIZE(enum%desc)
           ALLOCATE(enum%desc(i)%chars(1))
           enum%desc(i)%chars(1:1) = ' '
        END DO
     END IF
  END IF
END SUBROUTINE enum_create

! *****************************************************************************
!> \brief retains the given enumeration
!> \param enum the obect to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE enum_retain(enum,error)
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief releases the given enumeration
!> \param enum the obect to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE enum_release(enum,error)
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
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
          DO i=1, SIZE(enum%desc)
             DEALLOCATE(enum%desc(i)%chars)
          END DO
          DEALLOCATE(enum%desc,stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          DEALLOCATE(enum,stat=stat)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(enum)
  END SUBROUTINE enum_release

! *****************************************************************************
!> \brief maps an integer to a string
!> \param enum the enumeration to use for the mapping
!> \param i the value to map
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
FUNCTION enum_i2c(enum,i,error) RESULT(res)
    TYPE(enumeration_type), POINTER          :: enum
    INTEGER, INTENT(in)                      :: i
    TYPE(cp_error_type), INTENT(inout)       :: error
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
        IF (enum%strict) THEN
         DO j=1,SIZE(enum%desc)
          PRINT *, TRIM(a2s(enum%desc(j)%chars))
          PRINT *, TRIM(enum%c_vals(j))
         ENDDO
          PRINT *, enum%i_vals
        END IF
        CALL cp_assert(.NOT.enum%strict,cp_failure_level,cp_assertion_failed,&
             routineP,"invalid value for enumeration:"//cp_to_string(i),&
             error,failure)
        res=ADJUSTL(cp_to_string(i))
     END IF
  END IF
END FUNCTION enum_i2c

! *****************************************************************************
!> \brief maps a string to an integer
!> \param enum the enumeration to use for the mapping
!> \param c the value to map
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \author fawzi
! *****************************************************************************
FUNCTION enum_c2i(enum,c,error) RESULT(res)
    TYPE(enumeration_type), POINTER          :: enum
    CHARACTER(len=*), INTENT(in)             :: c
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'enum_c2i', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: upc
    INTEGER                                  :: iostat, j
    LOGICAL                                  :: failure, found

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(enum),cp_failure_level,routineP,error,failure)
  CPPrecondition(enum%ref_count>0,cp_failure_level,routineP,error,failure)
  upc=c
  CALL uppercase(upc)
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
             error=error,failure=failure,info=TRIM(upc))
        READ(c,"(i10)",iostat=iostat) res
        CALL cp_assert(iostat==0,cp_failure_level,cp_assertion_failed,&
             routineP,"invalid value for enumeration2:"//TRIM(c),&
             error=error,failure=failure,info=TRIM(upc))
     END IF
  END IF
END FUNCTION enum_c2i

END MODULE input_enumeration_types
