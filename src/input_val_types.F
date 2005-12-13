!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_val_types [1.0] *
!!
!!   NAME
!!     input_val_types
!!
!!   FUNCTION
!!     a wrapper for basic fortran types.
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
MODULE input_val_types
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enum_release,&
                                             enum_retain,&
                                             enumeration_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='input_val_types'

  INTEGER, SAVE, PRIVATE :: last_val_id=0

  PUBLIC :: val_p_type, val_type
  PUBLIC :: val_create, val_retain, val_release, val_get, val_write,&
            val_write_internal, val_duplicate
!***
!****************************************************************************

  INTEGER, PARAMETER, PUBLIC :: no_t=0, logical_t=1, &
       integer_t=2, real_t=3, char_t=4, enum_t=5
  
!!****s* input_val_types/val_p_type [1.0] *
!!
!!   NAME
!!     val_p_type
!!
!!   FUNCTION
!!     pointer to a val, to create arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - val: to pointer to the val
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE val_p_type
     TYPE(val_type), POINTER :: val
  END TYPE val_p_type
!!***
  !****************************************************************************

!!****s* input_val_types/val_type *
!!
!!   NAME
!!     val_type
!!
!!   FUNCTION
!!     a type to  have a wrapper that stores any basic fortran type
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - type_of_var: type stored in the val (should be one of no_t,
!!       integer_t, logical_t, real_t, char_t)
!!     - l_val, i_val, c_val, r_val: arrays with logical,integer,character
!!       or real values. Only one should be associated (and namely the one
!!       specified in type_of_var).
!!     - enum: an enumaration to map char to integers
!!
!!   AUTHOR
!!     fawzi
!!
!!   SOURCE
  !***************************************************************************
  TYPE val_type
     INTEGER :: ref_count, id_nr, type_of_var
     LOGICAL, DIMENSION(:), POINTER :: l_val
     INTEGER, DIMENSION(:), POINTER :: i_val
     CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: &
          c_val
     REAL(kind=dp), DIMENSION(:), POINTER :: r_val
     TYPE(enumeration_type), POINTER :: enum
  END TYPE val_type
!!***
  !****************************************************************************
CONTAINS

!!****f* input_val_types/val_create *
!!
!!   NAME
!!     val_create
!!
!!   FUNCTION
!!     creates a keyword value
!!
!!   NOTES
!!     using an enumeration only i_val/i_vals/i_vals_ptr are accepted
!!
!!   INPUTS
!!     - val: the object to be created
!!     - l_val,i_val,r_val,c_val: a logical,integer,real,string
!!       to be stored in the val
!!     - l_vals, i_vals, r_vals, c_vals: an array of logicals, integers, reals,
!!       characters to be stored in val
!!     - l_vals_ptr, i_vals_ptr, r_vals_ptr, c_vals_ptr: an array of logicals,
!!       ... to be stored in val, val will get the ownership of the pointer 
!!     - enum: the enumaration type this value is using
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_create(val,l_val,l_vals,l_vals_ptr,i_val,i_vals,i_vals_ptr,&
     r_val,r_vals,r_vals_ptr,c_val,c_vals,c_vals_ptr,enum,error)
    TYPE(val_type), POINTER                  :: val
    LOGICAL, INTENT(in), OPTIONAL            :: l_val
    LOGICAL, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: l_vals
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals_ptr
    INTEGER, INTENT(in), OPTIONAL            :: i_val
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: i_vals
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals_ptr
    REAL(KIND=DP), INTENT(in), OPTIONAL      :: r_val
    REAL(KIND=DP), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: r_vals
    REAL(KIND=DP), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_vals_ptr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: c_val
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: c_vals
    CHARACTER(LEN=DEFAULT_STRING_LENGTH), &
      DIMENSION(:), OPTIONAL, POINTER        :: c_vals_ptr
    TYPE(enumeration_type), OPTIONAL, &
      POINTER                                :: enum
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: narg, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(val),cp_failure_level,routineP,error,failure)
  ALLOCATE(val,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     NULLIFY(val%l_val,val%i_val,val%r_val,val%c_val,val%enum)
     val%type_of_var=no_t
     last_val_id=last_val_id+1
     val%id_nr=last_val_id
     val%ref_count=1

     narg=0
     val%type_of_var=no_t
     IF (PRESENT(l_val)) THEN
!FM        CPPrecondition(.NOT.PRESENT(l_vals),cp_failure_level,routineP,error,failure)
!FM        CPPrecondition(.NOT.PRESENT(l_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%l_val(1),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        val%l_val(1)=l_val
        val%type_of_var=logical_t
     END IF
     IF (PRESENT(l_vals)) THEN
!FM        CPPrecondition(.NOT.PRESENT(l_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%l_val(SIZE(l_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%l_val=l_vals
        val%type_of_var=logical_t
     END IF
     IF (PRESENT(l_vals_ptr)) THEN
        narg=narg+1
        val%l_val => l_vals_ptr
        val%type_of_var=logical_t
     END IF
     
     IF (PRESENT(r_val)) THEN
!FM        CPPrecondition(.NOT.PRESENT(r_vals),cp_failure_level,routineP,error,failure)
!FM        CPPrecondition(.NOT.PRESENT(r_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%r_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%r_val(1)=r_val
        val%type_of_var=real_t
     END IF
     IF (PRESENT(r_vals)) THEN
!FM        CPPrecondition(.NOT.PRESENT(r_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%r_val(SIZE(r_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%r_val=r_vals
        val%type_of_var=real_t
     END IF
     IF (PRESENT(r_vals_ptr)) THEN
        narg=narg+1
        val%r_val => r_vals_ptr
        val%type_of_var=real_t
     END IF
     
     IF (PRESENT(i_val)) THEN
!FM        CPPrecondition(.NOT.PRESENT(i_vals),cp_failure_level,routineP,error,failure)
!FM        CPPrecondition(.NOT.PRESENT(i_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%i_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%i_val(1)=i_val
        val%type_of_var=integer_t
     END IF
     IF (PRESENT(i_vals)) THEN
!FM        CPPrecondition(.NOT.PRESENT(i_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%i_val(SIZE(i_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%i_val=i_vals
        val%type_of_var=integer_t
     END IF
     IF (PRESENT(i_vals_ptr)) THEN
        narg=narg+1
        val%i_val => i_vals_ptr
        val%type_of_var=integer_t
     END IF
     
     IF (PRESENT(c_val)) THEN
!FM        CPPrecondition(.NOT.PRESENT(c_vals),cp_failure_level,routineP,error,failure)
!FM        CPPrecondition(.NOT.PRESENT(c_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%c_val(1),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%c_val(1)=c_val
        val%type_of_var=char_t
     END IF
     IF (PRESENT(c_vals)) THEN
!FM        CPPrecondition(.NOT.PRESENT(c_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%c_val(SIZE(c_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%c_val=c_vals
        val%type_of_var=char_t
     END IF
     IF (PRESENT(c_vals_ptr)) THEN
        narg=narg+1
        val%c_val => c_vals_ptr
        val%type_of_var=char_t
     END IF
     CPPostcondition(narg<=1,cp_failure_level,routineP,error,failure)
     IF (PRESENT(enum)) THEN
        IF (ASSOCIATED(enum)) THEN
           IF (val%type_of_var/=no_t.AND.val%type_of_var/=integer_t.AND.&
                val%type_of_var/=enum_t) THEN
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
           END IF
           IF (ASSOCIATED(val%i_val)) THEN
              val%type_of_var=enum_t
              val%enum=>enum
              CALL enum_retain(enum,error=error)
           END IF
        END IF
     END IF
     CPPostcondition(ASSOCIATED(val%enum).eqv.val%type_of_var==enum_t,cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE val_create
!***************************************************************************

!!****f* input_val_types/val_release *
!!
!!   NAME
!!     val_release
!!
!!   FUNCTION
!!     releases the given val
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - val: the val to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_release(val,error)
    TYPE(val_type), POINTER                  :: val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(val)) THEN
     CPPreconditionNoFail(val%ref_count>0,cp_failure_level,routineP,error)
     val%ref_count=val%ref_count-1
     IF (val%ref_count==0) THEN
        IF (ASSOCIATED(val%l_val)) THEN
           DEALLOCATE(val%l_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(val%i_val)) THEN
           DEALLOCATE(val%i_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(val%r_val)) THEN
           DEALLOCATE(val%r_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(val%c_val)) THEN
           DEALLOCATE(val%c_val,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL enum_release(val%enum,error=error)
        val%type_of_var=no_t
        DEALLOCATE(val,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(val)
END SUBROUTINE val_release
!***************************************************************************

!!****f* input_val_types/val_retain *
!!
!!   NAME
!!     val_retain
!!
!!   FUNCTION
!!     retains the given val
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - val: the val to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_retain(val,error)
    TYPE(val_type), POINTER                  :: val
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(val),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(val%ref_count>0,cp_failure_level,routineP,error)
     val%ref_count=val%ref_count+1
  END IF
END SUBROUTINE val_retain
!***************************************************************************

!!****f* input_val_types/val_get *
!!
!!   NAME
!!     val_get
!!
!!   FUNCTION
!!     returns the stored values
!!
!!   NOTES
!!     using an enumeration only i_val/i_vals/i_vals_ptr are accepted
!!
!!   INPUTS
!!     - val: the object from which you want to extract the values
!!     - l_val,i_val,r_val,c_val: gets a logical, integer, real, char
!!       from the val
!!     - logical_vals: gets an array of logicals, integers, reals, chars
!!       from the val
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_get(val,has_l,has_i,has_r,has_c,l_val,l_vals,i_val,&
     i_vals,r_val,r_vals,c_val,c_vals,type_of_var,enum,error)
    TYPE(val_type), POINTER                  :: val
    LOGICAL, INTENT(out), OPTIONAL           :: has_l, has_i, has_r, has_c, &
                                                l_val
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals
    INTEGER, INTENT(out), OPTIONAL           :: i_val
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals
    REAL(KIND=DP), INTENT(out), OPTIONAL     :: r_val
    REAL(KIND=DP), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_vals
    CHARACTER(LEN=*), INTENT(out), OPTIONAL  :: c_val
    CHARACTER(LEN=DEFAULT_STRING_LENGTH), &
      DIMENSION(:), OPTIONAL, POINTER        :: c_vals
    INTEGER, INTENT(out), OPTIONAL           :: type_of_var
    TYPE(enumeration_type), OPTIONAL, &
      POINTER                                :: enum
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "val_get", &
      routineP = moduleN//":"//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (PRESENT(has_l)) has_l=ASSOCIATED(val%l_val)
  IF (PRESENT(has_i)) has_i=ASSOCIATED(val%i_val)
  IF (PRESENT(has_r)) has_r=ASSOCIATED(val%r_val)
  IF (PRESENT(has_c)) has_c=ASSOCIATED(val%c_val)
  IF (PRESENT(l_vals)) l_vals => val%l_val
  IF (PRESENT(l_val)) THEN
     IF (ASSOCIATED(val%l_val)) THEN
        IF (SIZE(val%l_val)>0) THEN
           l_val=val%l_val(1)
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  IF (PRESENT(i_vals)) i_vals => val%i_val
  IF (PRESENT(i_val)) THEN
     IF (ASSOCIATED(val%i_val)) THEN
        IF (SIZE(val%i_val)>0) THEN
           i_val=val%i_val(1)
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  IF (PRESENT(r_vals)) r_vals => val%r_val
  IF (PRESENT(r_val)) THEN
     IF (ASSOCIATED(val%r_val)) THEN
        IF (SIZE(val%r_val)>0) THEN
           r_val=val%r_val(1)
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  IF (PRESENT(c_vals)) c_vals => val%c_val
  IF (PRESENT(c_val)) THEN
     IF (ASSOCIATED(val%c_val)) THEN
        IF (SIZE(val%c_val)>0) THEN
           c_val=val%c_val(1)
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE if (associated(val%i_val).and.associated(val%enum)) then
        IF (SIZE(val%i_val)>0) THEN
           c_val=enum_i2c(val%enum,val%i_val(1),error=error)
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  IF (PRESENT(type_of_var)) type_of_var=val%type_of_var

  IF (PRESENT(enum)) enum => val%enum

END SUBROUTINE val_get
!***************************************************************************

!!****f* input_val_types/val_write *
!!
!!   NAME
!!     val_write
!!
!!   FUNCTION
!!     writes out the valuse stored in the val
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - val: the val to write
!!     - unit_nr: the unit to write to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_write(val,unit_nr,error)
    TYPE(val_type), POINTER                  :: val
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_string
    INTEGER                                  :: i, l
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(val)) THEN
     SELECT CASE(val%type_of_var)
     CASE (logical_t)
        IF (ASSOCIATED(val%l_val)) THEN
           DO i=1,SIZE(val%l_val)
              IF (MODULO(i,20)==0) WRITE(unit=unit_nr,fmt="(' \')")
              WRITE(unit=unit_nr,fmt="(' ',l1)",advance="NO")&
                   val%l_val(i)
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (integer_t)
        IF (ASSOCIATED(val%i_val)) THEN
           DO i=1,SIZE(val%i_val)
              IF (MODULO(i,10)==0) WRITE(unit=unit_nr,fmt="(' \')")
              WRITE(unit=unit_nr,fmt="(' ',i14)",advance="NO")&
                   val%i_val(i)
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (real_t)
        IF (ASSOCIATED(val%r_val)) THEN
           DO i=1,SIZE(val%r_val)
              IF (MODULO(i,10)==0) WRITE(unit=unit_nr,fmt="(' \')")
              WRITE(unit=unit_nr,fmt="(' ',es16.8)",advance="NO")&
                   val%r_val(i)
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (char_t)
        IF (ASSOCIATED(val%c_val)) THEN
           l=0
           DO i=1,SIZE(val%c_val)
              IF (l>10.AND.l+LEN_TRIM(val%c_val(i))>76)THEN
                 WRITE(unit=unit_nr,fmt="(' \')")
                 l=0
              ELSE
                 l=l+LEN_TRIM(val%c_val(i))+3
              END IF
              WRITE(unit=unit_nr,fmt='(a)',advance="NO") TRIM(val%c_val(i))
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (enum_t)
        IF (ASSOCIATED(val%i_val)) THEN
           l=0
           DO i=1,SIZE(val%i_val)
              c_string=enum_i2c(val%enum,val%i_val(i),error=error)
              IF (l>10.AND.l+LEN_TRIM(c_string)>76)THEN
                 WRITE(unit=unit_nr,fmt="(' \')")
                 l=0
              ELSE
                 l=l+LEN_TRIM(c_string)+3
              END IF
              WRITE(unit=unit_nr,fmt="(' ',a)",advance="NO") TRIM(c_string)
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
        
     CASE(no_t)
        WRITE(unit=unit_nr,fmt="(' *empty*')",advance="NO")
     CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"unexpected type_of_var for val ",error,failure)
     END SELECT
  ELSE
     WRITE(unit=unit_nr,fmt="(' *null*')",advance="NO")
  END IF
  WRITE(unit=unit_nr,fmt="()")
END SUBROUTINE val_write

! *****************************************************************************

  SUBROUTINE val_write_internal(val,string,error)

!   Purpose: Write values to an internal file, i.e. string variable.

!   History: - Creation (10.03.2005,MK)

!   ***************************************************************************

    TYPE(val_type), POINTER                  :: val
    CHARACTER(LEN=*), INTENT(OUT)            :: string
    TYPE(cp_error_type), INTENT(INOUT),&
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "val_write_internal", &
                                   routineP = moduleN//":"//routineN

    CHARACTER(LEN=default_string_length) :: enum_string
    INTEGER                              :: i,ipos
    LOGICAL                              :: failure

!   ---------------------------------------------------------------------------

    failure = .FALSE.

    string = ""

    IF (ASSOCIATED(val)) THEN

      SELECT CASE(val%type_of_var)
      CASE (logical_t)
        IF (ASSOCIATED(val%l_val)) THEN
          DO i=1,SIZE(val%l_val)
            WRITE (UNIT=string(2*i-1:),FMT="(L2)") val%l_val(i)
          END DO
        ELSE
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
      CASE (integer_t)
        IF (ASSOCIATED(val%i_val)) THEN
          DO i=1,SIZE(val%i_val)
            WRITE (UNIT=string(12*i-11:),FMT="(I12)") val%i_val(i)
          END DO
        ELSE
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
      CASE (real_t)
        IF (ASSOCIATED(val%r_val)) THEN
          DO i=1,SIZE(val%r_val)
            WRITE (UNIT=string(16*i-15:),FMT="(ES16.8)") val%r_val(i)
          END DO
        ELSE
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
      CASE (char_t)
        IF (ASSOCIATED(val%c_val)) THEN
          ipos = 1
          DO i=1,SIZE(val%c_val)
            WRITE (UNIT=string(ipos:),FMT="(A)") TRIM(ADJUSTL(val%c_val(i)))
            ipos = ipos + LEN_TRIM(ADJUSTL(val%c_val(i))) + 1
          END DO
        ELSE
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
      CASE (enum_t)
        IF (ASSOCIATED(val%i_val)) THEN
          DO i=1,SIZE(val%i_val)
            enum_string = enum_i2c(val%enum,val%i_val(i),error)
            WRITE (UNIT=string,FMT="(A)") TRIM(ADJUSTL(enum_string))
          END DO
        ELSE
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
      CASE default
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                       routineP,"unexpected type_of_var for val ",error,&
                       failure)
      END SELECT

    END IF

  END SUBROUTINE val_write_internal

! *****************************************************************************

!!****f* input_val_types/val_duplicate *
!!
!!   NAME
!!     val_duplicate
!!
!!   FUNCTION
!!     creates a copy of the given value
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - val_in: the value to copy
!!     - val_out: the value tha will be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE val_duplicate(val_in,val_out,error)
  TYPE(val_type), POINTER :: val_in,val_out
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='val_duplicate',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.
  CPPrecondition(ASSOCIATED(val_in),cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(val_out),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(val_out,stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     last_val_id=last_val_id+1
     val_out%id_nr=last_val_id
     val_out%type_of_var=val_in%type_of_var
     val_out%ref_count=1
     val_out%enum => val_in%enum
     IF (ASSOCIATED(val_out%enum)) CALL enum_retain(val_out%enum,error=error)

     NULLIFY(val_out%l_val,val_out%i_val,val_out%c_val,val_out%r_val)
     IF (ASSOCIATED(val_in%l_val)) THEN
        ALLOCATE(val_out%l_val(SIZE(val_in%l_val)),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        val_out%l_val=val_in%l_val
     END IF
  !   IF (ASSOCIATED(val_in%l_val)) THEN
  !      ALLOCATE(val_out%l_val(SIZE(val_in%l_val)),stat=stat)
  !      CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  !      val_out%l_val=val_in%l_val
  !   END IF
     IF (ASSOCIATED(val_in%i_val)) THEN
        ALLOCATE(val_out%i_val(SIZE(val_in%i_val)),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        val_out%i_val=val_in%i_val
     END IF
     IF (ASSOCIATED(val_in%r_val)) THEN
        ALLOCATE(val_out%r_val(SIZE(val_in%r_val)),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        val_out%r_val=val_in%r_val
     END IF
     IF (ASSOCIATED(val_in%c_val)) THEN
        ALLOCATE(val_out%c_val(SIZE(val_in%c_val)),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        val_out%c_val=val_in%c_val
     END IF
  END IF
END SUBROUTINE val_duplicate
!***************************************************************************

END MODULE input_val_types
