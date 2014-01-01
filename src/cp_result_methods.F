!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  set of type/routines to handle the storage of results in force_envs
!> \author fschiff (12.2007)
!> \par    History
!>         - 10.2008 Teodoro Laino [tlaino] - University of Zurich
!>                   major rewriting:
!>                   - information stored in a proper type (not in a character!)
!>                   - module more lean
! *****************************************************************************
MODULE cp_result_methods
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: &
       cp_result_clean, cp_result_copy, cp_result_create, cp_result_release, &
       cp_result_type, cp_result_value_copy, cp_result_value_create, &
       cp_result_value_init, cp_result_value_p_reallocate
  USE f77_blas
  USE input_val_types,                 ONLY: integer_t,&
                                             logical_t,&
                                             real_t
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_result_methods'

  PUBLIC :: put_results,&
            test_for_result,&
            get_results,&
            cp_results_erase,&
            cp_results_mp_bcast

  INTERFACE put_results
     MODULE PROCEDURE put_result_r1, put_result_r2
  END INTERFACE

  INTERFACE get_results
     MODULE PROCEDURE get_result_r1, get_result_r2, get_nreps
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Store a 1D array of reals in result_list
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE put_result_r1(results,description,values,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length)     :: description
    REAL(KIND=dp), DIMENSION(:)              :: values
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_result_r1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isize, jsize
    LOGICAL                                  :: check, failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    check  = SIZE(results%result_label)==SIZE(results%result_value)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    isize=SIZE(results%result_label)
    jsize=SIZE(values)

    CALL reallocate(results%result_label,1,isize+1)
    CALL cp_result_value_p_reallocate(results%result_value,1,isize+1,error)

    results%result_label(isize+1) = description
    CALL cp_result_value_init(results%result_value(isize+1)%value, real_t, jsize, error)
    results%result_value(isize+1)%value%real_type = values

  END SUBROUTINE put_result_r1

! *****************************************************************************
!> \brief Store a 2D array of reals in result_list
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE put_result_r2(results,description,values,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length)     :: description
    REAL(KIND=dp), DIMENSION(:, :)           :: values
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_result_r2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isize, jsize
    LOGICAL                                  :: check, failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    check  = SIZE(results%result_label)==SIZE(results%result_value)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    isize=SIZE(results%result_label)
    jsize=SIZE(values,1)*SIZE(values,2)

    CALL reallocate(results%result_label,1,isize+1)
    CALL cp_result_value_p_reallocate(results%result_value,1,isize+1,error)

    results%result_label(isize+1) = description
    CALL cp_result_value_init(results%result_value(isize+1)%value, real_t, jsize, error)
    results%result_value(isize+1)%value%real_type = RESHAPE(values,(/jsize/))

  END SUBROUTINE put_result_r2

! *****************************************************************************
!> \brief test for a certain result in the result_list
!> \par History
!>      10.2013
!> \author Mandes
! *****************************************************************************
  FUNCTION test_for_result(results,description, error) RESULT(res_exist)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length)     :: description
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: res_exist

    CHARACTER(len=*), PARAMETER :: routineN = 'test_for_result', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nlist
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    nlist=SIZE(results%result_value)
    res_exist = .FALSE.
    DO i = 1, nlist
       IF(TRIM(results%result_label(i))==TRIM(description)) THEN
         res_exist = .TRUE.
         EXIT
       END IF
    END DO

  END FUNCTION test_for_result

! *****************************************************************************
!> \brief gets the required part out of the result_list
!> \param nval      : if more than one entry for a given description is given you may choose
!>                    which entry you want
!> \param n_rep     : integer indicating how many times the section exists in result_list
!> \param n_entries : gets the number of lines used for a given description
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE get_result_r1(results,description,values,nval,n_rep,n_entries,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length)     :: description
    REAL(KIND=dp), DIMENSION(:)              :: values
    INTEGER, OPTIONAL                        :: nval, n_rep, n_entries
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_result_r1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, nlist, nrep, size_res, &
                                                size_values
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    nlist=SIZE(results%result_value)
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    CPPostcondition(SIZE(results%result_label)==nlist,cp_failure_level,routineP,error,failure)
    nrep = 0
    DO i = 1, nlist
       IF(TRIM(results%result_label(i))==TRIM(description)) nrep = nrep + 1
    END DO

    IF(PRESENT(n_rep))THEN
       n_rep=nrep
    END IF

    CALL cp_assert(nrep.GT.0,cp_failure_level,cp_assertion_failed,routineP,&
         " Trying to access result ("//TRIM(description)//") which was never stored! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    DO i=1,nlist
       IF(TRIM(results%result_label(i))==TRIM(description))THEN
          CALL cp_assert(results%result_value(i)%value%type_in_use==real_t,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "Attempt to retrieve a RESULT which is not a REAL! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)

          size_res=SIZE(results%result_value(i)%value%real_type)
          EXIT
       END IF
    END DO
    IF(PRESENT(n_entries)) n_entries = size_res
    size_values = SIZE(values,1)
    IF(PRESENT(nval))THEN
       CPPostcondition(size_res==size_values,cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(nrep*size_res==size_values,cp_failure_level,routineP,error,failure)
    END IF
    k=0
    DO i = 1,nlist
       IF(TRIM(results%result_label(i))==TRIM(description))THEN
          k = k + 1
          IF(PRESENT(nval))THEN
             IF(k==nval)THEN
                values = results%result_value(i)%value%real_type
                EXIT
             END IF
          ELSE
             values((k-1)*size_res+1:k*size_res) =  results%result_value(i)%value%real_type
          END IF
       END IF
    END DO

  END SUBROUTINE get_result_r1

! *****************************************************************************
!> \brief gets the required part out of the result_list
!> \param nval      : if more than one entry for a given description is given you may choose
!>                    which entry you want
!> \param n_rep     : integer indicating how many times the section exists in result_list
!> \param n_entries : gets the number of lines used for a given description
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE get_result_r2(results,description,values,nval,n_rep,n_entries,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length)     :: description
    REAL(KIND=dp), DIMENSION(:, :)           :: values
    INTEGER, OPTIONAL                        :: nval, n_rep, n_entries
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_result_r2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, k, nlist, nrep, size_res, &
                                                size_values
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    nlist=SIZE(results%result_value)
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    CPPostcondition(SIZE(results%result_label)==nlist,cp_failure_level,routineP,error,failure)
    nrep = 0
    DO i = 1, nlist
       IF(TRIM(results%result_label(i))==TRIM(description)) nrep = nrep + 1
    END DO

    IF(PRESENT(n_rep))THEN
       n_rep=nrep
    END IF

    CALL cp_assert(nrep.GT.0,cp_failure_level,cp_assertion_failed,routineP,&
         " Trying to access result ("//TRIM(description)//") which was never stored! "//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    DO i=1,nlist
       IF(TRIM(results%result_label(i))==TRIM(description))THEN
          CALL cp_assert(results%result_value(i)%value%type_in_use==real_t,&
               cp_failure_level,cp_assertion_failed,routineP,&
               "Attempt to retrieve a RESULT which is not a REAL! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)

          size_res=SIZE(results%result_value(i)%value%real_type)
          EXIT
       END IF
    END DO
    IF(PRESENT(n_entries)) n_entries = size_res
    size_values = SIZE(values,1)*SIZE(values,2)
    IF(PRESENT(nval))THEN
       CPPostcondition(size_res==size_values,cp_failure_level,routineP,error,failure)
    ELSE
       CPPostcondition(nrep*size_res==size_values,cp_failure_level,routineP,error,failure)
    END IF
    k=0
    DO i = 1,nlist
       IF(TRIM(results%result_label(i))==TRIM(description))THEN
          k = k + 1
          IF(PRESENT(nval))THEN
             IF(k==nval)THEN
                values = RESHAPE(results%result_value(i)%value%real_type,(/SIZE(values,1),SIZE(values,2)/))
                EXIT
             END IF
          ELSE
             values((k-1)*size_res+1:k*size_res,:) =   RESHAPE(results%result_value(i)%value%real_type,&
                                                               (/SIZE(values,1),SIZE(values,2)/))
          END IF
       END IF
    END DO

  END SUBROUTINE get_result_r2

! *****************************************************************************
!> \brief gets the required part out of the result_list
!> \param n_rep     : integer indicating how many times the section exists in result_list
!> \param n_entries : gets the number of lines used for a given description
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE get_nreps(results,description,n_rep,n_entries,type_in_use,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length)     :: description
    INTEGER, OPTIONAL                        :: n_rep, n_entries, type_in_use
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_nreps', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: I, nlist
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    nlist=SIZE(results%result_value)
    CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
    CPPostcondition(SIZE(results%result_label)==nlist,cp_failure_level,routineP,error,failure)
    IF(PRESENT(n_rep))THEN
       n_rep=0
       DO i=1,nlist
          IF(TRIM(results%result_label(i))==TRIM(description)) n_rep = n_rep + 1
       END DO
    END IF
    IF(PRESENT(n_entries))THEN
       n_entries = 0
       DO i=1,nlist
          IF(TRIM(results%result_label(i))==TRIM(description))THEN
             SELECT CASE(results%result_value(i)%value%type_in_use)
             CASE(real_t)
                n_entries = n_entries + SIZE(results%result_value(i)%value%real_type)
             CASE(integer_t)
                n_entries = n_entries + SIZE(results%result_value(i)%value%integer_type)
             CASE(logical_t)
                n_entries = n_entries + SIZE(results%result_value(i)%value%logical_type)
             CASE DEFAULT
                ! Type not implemented in cp_result_type
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             EXIT
          END IF
       END DO
    END IF
    IF (PRESENT(type_in_use)) THEN
       DO i=1,nlist
          IF(TRIM(results%result_label(i))==TRIM(description))THEN
             type_in_use = results%result_value(i)%value%type_in_use
             EXIT
          END IF
       END DO
    END IF
  END SUBROUTINE get_nreps

! *****************************************************************************
!> \brief erase a part of  result_list
!> \param nval : if more than one entry for a given description is given you may choose
!>               which entry you want to delete
!> \par History
!>      12.2007 created
!>      10.2008 Teodoro Laino [tlaino] - major rewriting
!> \author fschiff
! *****************************************************************************
  SUBROUTINE cp_results_erase(results,description,nval,error)
    TYPE(cp_result_type), POINTER            :: results
    CHARACTER(LEN=default_string_length), &
      OPTIONAL                               :: description
    INTEGER, OPTIONAL                        :: nval
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_results_erase', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: entry_deleted, i, k, &
                                                new_size, nlist, nrep, stat
    LOGICAL                                  :: failure
    TYPE(cp_result_type), POINTER            :: clean_results

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    new_size = 0
    IF(PRESENT(description))THEN
       CPPostcondition(description(1:1)=='[',cp_failure_level,routineP,error,failure)
       nlist=SIZE(results%result_value)
       nrep=0
       DO i=1,nlist
          IF(TRIM(results%result_label(i))==TRIM(description)) nrep = nrep + 1
       END DO
       IF(nrep.NE.0)THEN
          k             = 0
          entry_deleted = 0
          DO i=1,nlist
             IF(TRIM(results%result_label(i))==TRIM(description))THEN
                k = k + 1
                IF (PRESENT(nval)) THEN
                   IF (nval==k) THEN
                      entry_deleted = entry_deleted + 1
                      EXIT
                   END IF
                ELSE
                   entry_deleted = entry_deleted + 1
                END IF
             END IF
          END DO
          CPPostcondition(nlist-entry_deleted>=0,cp_failure_level,routineP,error,failure)
          new_size = nlist-entry_deleted
          NULLIFY(clean_results)
          CALL cp_result_create(clean_results, error)
          CALL cp_result_clean(clean_results, error)
          ALLOCATE(clean_results%result_label(new_size), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(clean_results%result_value(new_size), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, new_size
             NULLIFY(clean_results%result_value(i)%value)
             CALL cp_result_value_create(clean_results%result_value(i)%value, error)
          END DO
          k = 0
          DO i =1,nlist
             IF(TRIM(results%result_label(i))/=TRIM(description)) THEN
                k = k + 1
                clean_results%result_label(k) = results%result_label(i)
                CALL cp_result_value_copy(clean_results%result_value(k)%value,&
                     results%result_value(i)%value, error)
             END IF
          END DO
          CALL cp_result_copy(clean_results, results, error)
          CALL cp_result_release(clean_results, error)
       END IF
    ELSE
       CALL cp_result_clean(results, error)
       ALLOCATE(results%result_label(new_size), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(results%result_value(new_size), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE cp_results_erase

! *****************************************************************************
!> \brief broadcast results type
!> \author  10.2008 Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE cp_results_mp_bcast(results, source, para_env, error)
    TYPE(cp_result_type), POINTER            :: results
    INTEGER, INTENT(IN)                      :: source
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_results_mp_bcast', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nlist, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: size_value, type_in_use
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(results),cp_failure_level,routineP,error,failure)
    nlist = 0
    IF (para_env%mepos==source) nlist = SIZE(results%result_value)
    CALL mp_bcast(nlist, source, para_env%group)

    ALLOCATE(size_value(nlist), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(type_in_use(nlist), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (para_env%mepos==source) THEN
       DO i = 1, nlist
          CALL get_nreps(results,description=results%result_label(i),&
               n_entries=size_value(i),type_in_use=type_in_use(i),&
               error=error)
       END DO
    END IF
    CALL mp_bcast(size_value, source, para_env%group)
    CALL mp_bcast(type_in_use, source, para_env%group)

    IF (para_env%mepos/=source) THEN
       CALL cp_result_clean(results, error)
       ALLOCATE(results%result_value(nlist),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(results%result_label(nlist),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, nlist
          results%result_label(i) = ""
          NULLIFY(results%result_value(i)%value)
          CALL cp_result_value_create(results%result_value(i)%value,error)
          CALL cp_result_value_init(results%result_value(i)%value,&
               type_in_use=type_in_use(i),size_value=size_value(i), error=error)
       END DO
    END IF
    DO i = 1, nlist
       CALL mp_bcast(results%result_label(i), source, para_env%group)
       SELECT CASE(results%result_value(i)%value%type_in_use)
       CASE(real_t)
          CALL mp_bcast(results%result_value(i)%value%real_type, source, para_env%group)
       CASE(integer_t)
          CALL mp_bcast(results%result_value(i)%value%integer_type, source, para_env%group)
       CASE(logical_t)
          CALL mp_bcast(results%result_value(i)%value%logical_type, source, para_env%group)
       CASE DEFAULT
          ! Type not implemented in cp_result_type
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END DO
    DEALLOCATE(type_in_use, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(size_value, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE cp_results_mp_bcast

END MODULE cp_result_methods
