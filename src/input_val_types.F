!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief a wrapper for basic fortran types.
!> \par History
!>      06.2004 created
!> \author fawzi
! *****************************************************************************
MODULE input_val_types
  USE cp_parser_types,                 ONLY: default_continuation_character
  USE cp_units,                        ONLY: cp_unit_create,&
                                             cp_unit_desc,&
                                             cp_unit_from_cp2k,&
                                             cp_unit_from_cp2k1,&
                                             cp_unit_release,&
                                             cp_unit_retain,&
                                             cp_unit_type
  USE f77_blas
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_val_types'

  INTEGER, SAVE, PRIVATE :: last_val_id=0

  PUBLIC :: val_p_type, val_type
  PUBLIC :: val_create, val_retain, val_release, val_get, val_write,&
            val_write_internal, val_duplicate
!***

  INTEGER, PARAMETER, PUBLIC :: no_t=0, logical_t=1, &
       integer_t=2, real_t=3, char_t=4, enum_t=5, lchar_t=6

! *****************************************************************************
!> \brief pointer to a val, to create arrays of pointers
!> \param val to pointer to the val
!> \author fawzi
! *****************************************************************************
  TYPE val_p_type
     TYPE(val_type), POINTER :: val
  END TYPE val_p_type

! *****************************************************************************
!> \brief a type to  have a wrapper that stores any basic fortran type
!> \param type_of_var type stored in the val (should be one of no_t,
!>        integer_t, logical_t, real_t, char_t)
!> \param l_val , i_val, c_val, r_val: arrays with logical,integer,character
!>        or real values. Only one should be associated (and namely the one
!>        specified in type_of_var).
!> \param enum an enumaration to map char to integers
!> \author fawzi
! *****************************************************************************
  TYPE val_type
     INTEGER :: ref_count, id_nr, type_of_var
     LOGICAL, DIMENSION(:), POINTER :: l_val
     INTEGER, DIMENSION(:), POINTER :: i_val
     CHARACTER(len=default_string_length), DIMENSION(:), POINTER :: &
          c_val
     REAL(kind=dp), DIMENSION(:), POINTER :: r_val
     TYPE(enumeration_type), POINTER :: enum
  END TYPE val_type
CONTAINS

! *****************************************************************************
!> \brief creates a keyword value
!> \param val the object to be created
!> \param l_val ,i_val,r_val,c_val,lc_val: a logical,integer,real,string, long
!>        string to be stored in the val
!> \param l_vals , i_vals, r_vals, c_vals: an array of logicals,
!>        integers, reals, characters, long strings to be stored in val
!> \param l_vals_ptr , i_vals_ptr, r_vals_ptr, c_vals_ptr: an array of logicals,
!>        ... to be stored in val, val will get the ownership of the pointer
!> \param enum the enumaration type this value is using
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      using an enumeration only i_val/i_vals/i_vals_ptr are accepted
!> \author fawzi
! *****************************************************************************
SUBROUTINE val_create(val,l_val,l_vals,l_vals_ptr,i_val,i_vals,i_vals_ptr,&
     r_val,r_vals,r_vals_ptr,c_val,c_vals,c_vals_ptr,lc_val,lc_vals,&
     lc_vals_ptr,enum,error)
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
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: c_vals_ptr
    CHARACTER(LEN=*), INTENT(in), OPTIONAL   :: lc_val
    CHARACTER(LEN=*), DIMENSION(:), &
      INTENT(in), OPTIONAL                   :: lc_vals
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: lc_vals_ptr
    TYPE(enumeration_type), OPTIONAL, &
      POINTER                                :: enum
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, len_c, narg, nVal, stat
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
        CPPrecondition(LEN_TRIM(c_val)<=default_string_length,cp_failure_level,routineP,error,failure)
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
        CPPrecondition(ALL(LEN_TRIM(c_vals)<=default_string_length),cp_failure_level,routineP,error,failure)
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
     IF (PRESENT(lc_val)) THEN
!FM        CPPrecondition(.NOT.PRESENT(lc_vals),cp_failure_level,routineP,error,failure)
!FM        CPPrecondition(.NOT.PRESENT(lc_vals_ptr),cp_failure_level,routineP,error,failure)
        narg=narg+1
        len_c=LEN_TRIM(lc_val)
        nVal=MAX(1,CEILING(REAL(len_c,dp)/80._dp))
        ALLOCATE(val%c_val(nVal),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

        IF (len_c == 0) THEN
           val%c_val(1) = ""
        ELSE
           DO i=1,nVal
              val%c_val(i)=lc_val((i-1)*default_string_length+1:&
                   MIN(len_c,i*default_string_length))
           END DO
        END IF
        val%type_of_var=lchar_t
     END IF
     IF (PRESENT(lc_vals)) THEN
        CPPrecondition(ALL(LEN_TRIM(lc_vals)<=default_string_length),cp_failure_level,routineP,error,failure)
        narg=narg+1
        ALLOCATE(val%c_val(SIZE(lc_vals)),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        val%c_val=lc_vals
        val%type_of_var=lchar_t
     END IF
     IF (PRESENT(lc_vals_ptr)) THEN
        narg=narg+1
        val%c_val => lc_vals_ptr
        val%type_of_var=lchar_t
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

! *****************************************************************************
!> \brief releases the given val
!> \param val the val to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE val_release(val,error)
    TYPE(val_type), POINTER                  :: val
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief retains the given val
!> \param val the val to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE val_retain(val,error)
    TYPE(val_type), POINTER                  :: val
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief returns the stored values
!> \param val the object from which you want to extract the values
!> \param l_val ,i_val,r_val,c_val: gets a logical, integer, real, char
!>        from the val
!> \param logical_vals gets an array of logicals, integers, reals, chars
!>        from the val
!> \param len_c len_trim of c_val (if it was a lc_val, of type lchar_t
!>        it might be longet than default_string_length)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      using an enumeration only i_val/i_vals/i_vals_ptr are accepted
!>      add something like ignore_string_cut that if true does not warn if
!>      the c_val is too short to contain the string
!> \author fawzi
! *****************************************************************************
SUBROUTINE val_get(val,has_l,has_i,has_r,has_lc,has_c,l_val,l_vals,i_val,&
     i_vals,r_val,r_vals,c_val,c_vals,len_c,type_of_var,enum,error)
    TYPE(val_type), POINTER                  :: val
    LOGICAL, INTENT(out), OPTIONAL           :: has_l, has_i, has_r, has_lc, &
                                                has_c, l_val
    LOGICAL, DIMENSION(:), OPTIONAL, POINTER :: l_vals
    INTEGER, INTENT(out), OPTIONAL           :: i_val
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: i_vals
    REAL(KIND=DP), INTENT(out), OPTIONAL     :: r_val
    REAL(KIND=DP), DIMENSION(:), OPTIONAL, &
      POINTER                                :: r_vals
    CHARACTER(LEN=*), INTENT(out), OPTIONAL  :: c_val
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), OPTIONAL, POINTER        :: c_vals
    INTEGER, INTENT(out), OPTIONAL           :: len_c, type_of_var
    TYPE(enumeration_type), OPTIONAL, &
      POINTER                                :: enum
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, l_in, l_out
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (PRESENT(has_l)) has_l=ASSOCIATED(val%l_val)
  IF (PRESENT(has_i)) has_i=ASSOCIATED(val%i_val)
  IF (PRESENT(has_r)) has_r=ASSOCIATED(val%r_val)
  IF (PRESENT(has_c)) has_c=ASSOCIATED(val%c_val) ! use type_of_var?
  IF (PRESENT(has_lc)) has_lc=(val%type_of_var==lchar_t)
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
     l_out=LEN(c_val)
     IF (ASSOCIATED(val%c_val)) THEN
        IF (SIZE(val%c_val)>0) THEN
           IF (val%type_of_var==lchar_t) THEN
              l_in=default_string_length*(SIZE(val%c_val)-1)+&
                   LEN_TRIM(val%c_val(SIZE(val%c_val)))
              CALL cp_assert(l_out>=l_in,cp_warning_level,cp_assertion_failed,&
                   routineP,"val_get will truncate value, value beginning with '"//&
                   TRIM(val%c_val(1))//"' is too long for variable",error,failure)
              DO i=1,SIZE(val%c_val)
                 c_val((i-1)*default_string_length+1:MIN(l_out,i*default_string_length))=&
                      val%c_val(i)(1:MIN(80,l_out-(i-1)*default_string_length))
                 IF (l_out<=i*default_string_length) EXIT
              END DO
              IF (l_out>SIZE(val%c_val)*default_string_length) &
                   c_val(SIZE(val%c_val)*default_string_length+1:l_out)=""
           ELSE
              l_in=LEN_TRIM(val%c_val(1))
              CALL cp_assert(l_out>=l_in,cp_warning_level,cp_assertion_failed,&
                   routineP,"val_get will truncate value, value '"//&
                   TRIM(val%c_val(1))//"' is too long for variable",error,failure)
              c_val=val%c_val(1)
           END IF
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE IF (ASSOCIATED(val%i_val).AND.ASSOCIATED(val%enum)) THEN
        IF (SIZE(val%i_val)>0) THEN
           c_val=enum_i2c(val%enum,val%i_val(1),error=error)
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     ELSE
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
     END IF
  END IF

  IF (PRESENT(len_c)) THEN
     IF (ASSOCIATED(val%c_val)) THEN
        IF (SIZE(val%c_val)>0) THEN
           IF (val%type_of_var==lchar_t) THEN
              len_c=default_string_length*(SIZE(val%c_val)-1)+&
                   LEN_TRIM(val%c_val(SIZE(val%c_val)))
           ELSE
              len_c=LEN_TRIM(val%c_val(1))
           END IF
        ELSE
           len_c=-HUGE(0)
        END IF
     ELSE IF (ASSOCIATED(val%i_val).AND.ASSOCIATED(val%enum)) THEN
        IF (SIZE(val%i_val)>0) THEN
           len_c=LEN_TRIM(enum_i2c(val%enum,val%i_val(1),error=error))
        ELSE
           len_c=-HUGE(0)
        END IF
     ELSE
        len_c=-HUGE(0)
     END IF
  END IF

  IF (PRESENT(type_of_var)) type_of_var=val%type_of_var

  IF (PRESENT(enum)) enum => val%enum

END SUBROUTINE val_get

! *****************************************************************************
!> \brief writes out the valuse stored in the val
!> \param val the val to write
!> \param unit_nr the number of the unit to write to
!> \param unit the unit of mesure in wich the output should be written
!>        (overrides unit_str)
!> \param unit_str the unit of mesure in wich the output should be written
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      unit of mesure used only for reals
!> \author fawzi
! *****************************************************************************
SUBROUTINE val_write(val,unit_nr,unit,unit_str,fmt,error)
    TYPE(val_type), POINTER                  :: val
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_unit_type), OPTIONAL, POINTER    :: unit
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: unit_str, fmt
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: c_string, myfmt, rcval
    INTEGER                                  :: i, iend, item, j, l
    LOGICAL                                  :: failure
    TYPE(cp_unit_type), POINTER              :: my_unit

  failure=.FALSE.

  NULLIFY(my_unit)
  myfmt = ""
  IF (PRESENT(fmt)) myfmt = fmt
  IF (PRESENT(unit)) my_unit => unit
  IF (ASSOCIATED(my_unit)) THEN
     CALL cp_unit_retain(my_unit,error=error)
  ELSE IF (PRESENT(unit_str)) THEN
     CALL cp_unit_create(my_unit,unit_str,error=error)
  END IF
  IF (ASSOCIATED(val)) THEN
     SELECT CASE(val%type_of_var)
     CASE (logical_t)
        IF (ASSOCIATED(val%l_val)) THEN
           DO i=1,SIZE(val%l_val)
              IF (MODULO(i,20)==0) THEN
                 WRITE(unit=unit_nr,fmt="(' ',A)")default_continuation_character
                 WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//")",advance="NO")
              END IF
              WRITE(unit=unit_nr,fmt="(' ',l1)",advance="NO")&
                   val%l_val(i)
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (integer_t)
        IF (ASSOCIATED(val%i_val)) THEN
           item = 0
           i = 1
           loop_i: DO WHILE (i <= SIZE(val%i_val))
             item = item + 1
             IF (MODULO(item,10) == 0) THEN
               WRITE (UNIT=unit_nr,FMT="(1X,A)") default_continuation_character
               WRITE (UNIT=unit_nr,FMT="("//TRIM(myfmt)//")",ADVANCE="NO")
             END IF
             iend = i
             loop_j: DO j=i+1,SIZE(val%i_val)
               IF (val%i_val(j-1) + 1 == val%i_val(j)) THEN
                 iend = iend + 1
               ELSE
                 EXIT loop_j
               END IF
             END DO loop_j
             IF ((iend - i) > 1) THEN
               WRITE (UNIT=unit_nr,FMT="(1X,I0,A2,I0)",ADVANCE="NO")&
                val%i_val(i),"..",val%i_val(iend)
               i = iend
             ELSE
               WRITE (UNIT=unit_nr,FMT="(1X,I0)",ADVANCE="NO")&
                val%i_val(i)
             END IF
             i = i + 1
           END DO loop_i
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (real_t)
        IF (ASSOCIATED(val%r_val)) THEN
           DO i=1,SIZE(val%r_val)
              IF (MODULO(i,5)==0) THEN
                 WRITE(unit=unit_nr,fmt="(' ',A)")default_continuation_character
                 WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//")",advance="NO")
              END IF
              IF (ASSOCIATED(my_unit)) THEN
                 WRITE(rcval,"(ES25.16)")cp_unit_from_cp2k1(val%r_val(i),my_unit,error=error)
              ELSE
                 WRITE(rcval,"(ES25.16)")val%r_val(i)
              END IF
              WRITE(unit=unit_nr,fmt="(' ',A)",advance="NO")TRIM(rcval)
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (char_t)
        IF (ASSOCIATED(val%c_val)) THEN
           l=0
           DO i=1,SIZE(val%c_val)
              IF (i>1) WRITE (unit=unit_nr,fmt="(' ')",advance="NO")
              l=l+1
              IF (l>10.AND.l+LEN_TRIM(val%c_val(i))>76)THEN
                 WRITE(unit=unit_nr,fmt="('\')")
                 WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//")",advance="NO")
                 l=0
                 WRITE(unit=unit_nr,fmt='(a)',advance="NO") TRIM(val%c_val(i))
                 l=l+LEN_TRIM(val%c_val(i))+3
              ELSE IF (LEN_TRIM(val%c_val(i))>0) THEN
                 l=l+LEN_TRIM(val%c_val(i))
                 WRITE(unit=unit_nr,fmt='(a)',advance="NO") TRIM(val%c_val(i))
              ELSE
                 l=l+3
                 WRITE (unit=unit_nr,fmt="(a)",advance="NO") '" "'
              END IF
           END DO
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (lchar_t)
        IF (ASSOCIATED(val%c_val)) THEN
           l=0
           DO i=1,SIZE(val%c_val)-1
              WRITE(unit=unit_nr,fmt='(a)',advance="NO") val%c_val(i)
           END DO
           IF (SIZE(val%c_val)>0) THEN
              WRITE(unit=unit_nr,fmt='(a)',advance="NO") TRIM(val%c_val(SIZE(val%c_val)))
           END IF
        ELSE
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     CASE (enum_t)
        IF (ASSOCIATED(val%i_val)) THEN
           l=0
           DO i=1,SIZE(val%i_val)
              c_string=enum_i2c(val%enum,val%i_val(i),error=error)
              IF (l>10.AND.l+LEN_TRIM(c_string)>76)THEN
                 WRITE(unit=unit_nr,fmt="(' ',A)")default_continuation_character
                 WRITE(unit=unit_nr,fmt="("//TRIM(myfmt)//")",advance="NO")
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
  IF (ASSOCIATED(my_unit)) CALL cp_unit_release(my_unit,error=error)
  WRITE(unit=unit_nr,fmt="()")
END SUBROUTINE val_write

! *****************************************************************************
!> \brief   Write values to an internal file, i.e. string variable.
!> \author  MK
!> \date    10.03.2005
!> \par History
!>          17.01.2006, MK, Optional argument unit for the conversion to the external unit added
!> \version 1.0
! *****************************************************************************
  SUBROUTINE val_write_internal(val,string,unit,error)

    TYPE(val_type), POINTER                  :: val
    CHARACTER(LEN=*), INTENT(OUT)            :: string
    TYPE(cp_unit_type), OPTIONAL, POINTER    :: unit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_write_internal', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: enum_string
    INTEGER                                  :: i, ipos
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: value

! -------------------------------------------------------------------------

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
          IF (PRESENT(unit)) THEN
            DO i=1,SIZE(val%r_val)
              value = cp_unit_from_cp2k(value=val%r_val(i),&
                                        unit_str=cp_unit_desc(unit=unit,&
                                                              error=error),&
                                        error=error)
              WRITE (UNIT=string(16*i-15:),FMT="(ES16.8)") value
            END DO
          ELSE
            DO i=1,SIZE(val%r_val)
              WRITE (UNIT=string(16*i-15:),FMT="(ES16.8)") val%r_val(i)
            END DO
          END IF
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
      CASE (lchar_t)
        IF (ASSOCIATED(val%c_val)) THEN
          CALL val_get(val,c_val=string,error=error)
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
!> \brief creates a copy of the given value
!> \param val_in the value to copy
!> \param val_out the value tha will be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
SUBROUTINE val_duplicate(val_in,val_out,error)
    TYPE(val_type), POINTER                  :: val_in, val_out
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'val_duplicate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
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

END MODULE input_val_types
