!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_fm_pool [1.0] *
!!
!!   NAME
!!     cp_fm_pool
!!
!!   FUNCTION
!!     pool for full matrixes
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_pool
  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  USE kinds, ONLY: wp=>dp
  USE cp_para_types, ONLY: cp_para_env_type
  USE timings, ONLY: timeset, timestop
  USE cp_full_matrix, ONLY: cp_full_matrix_type, cp_fm_create, cp_fm_release
  USE cp_f_matrix_struct, ONLY: cp_f_matrix_struct_type, cp_fmstruct_retain,&
       cp_fmstruct_release
  USE cp_linked_list_fm, ONLY: cp_sll_fm_type, cp_sll_fm_dealloc,&
       cp_sll_fm_next, cp_sll_fm_get_first_el, cp_sll_fm_rm_first_el,&
       cp_sll_fm_insert_el
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_pool'
  INTEGER, SAVE, PRIVATE :: last_fm_pool_id_nr=0

  PUBLIC :: cp_fm_pool_type, cp_fm_pool_p_type
  PUBLIC :: fm_pool_create, fm_pool_retain, fm_pool_release,&
       fm_pool_create_matrix, fm_pool_give_back_matrix
!***
!****************************************************************************

!!****s* cp_fm_pool/cp_fm_pool_type [1.0] *
!!
!!   NAME
!!     cp_fm_pool_type
!!
!!   FUNCTION
!!     represent a pool of full matrixes
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (see /cp2k/doc/ReferenceCounting.html)
!!     - matrix_struct: matrix structure of the matrixes in this pool
!!     - cache: linked list with the matrixes in the pool
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_fm_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE(cp_f_matrix_struct_type), POINTER :: matrix_struct
     TYPE(cp_sll_fm_type), POINTER :: cache
  END TYPE cp_fm_pool_type
!!***
!****************************************************************************


!!****s* cp_fm_pool/cp_fm_pool_p_type [1.0] *
!!
!!   NAME
!!     cp_fm_pool_p_type
!!
!!   FUNCTION
!!     to create arrays of pools
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - pool: the pool
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_fm_pool_p_type
     TYPE(cp_fm_pool_type), POINTER :: pool
  END TYPE cp_fm_pool_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* cp_fm_pool/fm_pool_create [1.0] *
!!
!!   NAME
!!     fm_pool_create
!!
!!   FUNCTION
!!     creates a pool of full matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to create
!!     - matrix_struct: the structure of the matrixes that are stored in
!!       this pool
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fm_pool_create(pool, matrix_struct, error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_f_matrix_struct_type), POINTER :: matrix_struct
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    failure=.FALSE.

    ALLOCATE(pool, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pool%matrix_struct=> matrix_struct
       CALL cp_fmstruct_retain(pool%matrix_struct,error=error)
       last_fm_pool_id_nr=last_fm_pool_id_nr+1
       pool%id_nr=last_fm_pool_id_nr
       pool%ref_count=1
       NULLIFY(pool%cache)
    END IF
  END SUBROUTINE fm_pool_create
!***************************************************************************

!!****f* cp_fm_pool/cp_fm_pool_retain [1.0] *
!!
!!   NAME
!!     cp_fm_pool_retain
!!
!!   FUNCTION
!!     retains the pool (see cp2k/doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fm_pool_retain(pool,error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

       pool%ref_count=pool%ref_count+1
    END IF
  END SUBROUTINE fm_pool_retain
!***************************************************************************

!!****f* cp_fm_pool/fm_pool_flush_cache [1.0] *
!!
!!   NAME
!!     fm_pool_flush_cache
!!
!!   FUNCTION
!!     deallocates all the cached matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to flush
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fm_pool_flush_cache(pool,error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_flush_cache',&
         routineP=moduleN//':'//routineN
    TYPE(cp_sll_fm_type), POINTER :: iterator
    TYPE(cp_full_matrix_type), POINTER :: m_att
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       iterator => pool%cache
       DO
          IF (.NOT.cp_sll_fm_next(iterator,el_att=m_att)) EXIT
          CALL cp_fm_release(m_att,error=error)
       END DO
       CALL cp_sll_fm_dealloc(pool%cache,error=error)
    END IF
  END SUBROUTINE fm_pool_flush_cache
!***************************************************************************

!!****f* cp_fm_pool/cp_fm_pool_release [1.0] *
!!
!!   NAME
!!     cp_fm_pool_release
!!
!!   FUNCTION
!!     releases the given pool (see cp2k/doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fm_pool_release(pool,error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_pool_release',&
         routineP=moduleN//':'//routineN
    TYPE(cp_sll_fm_type), POINTER :: iterator
    TYPE(cp_full_matrix_type), POINTER :: matrix_att
    INTEGER :: stat
    failure=.FALSE.

    IF (ASSOCIATED(pool)) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
       pool%ref_count=pool%ref_count-1
       IF (pool%ref_count<1) THEN
          CALL fm_pool_flush_cache(pool,error=error)
          CALL cp_fmstruct_release(pool%matrix_struct,error=error)

          DEALLOCATE(pool,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE fm_pool_release
!***************************************************************************

!!****f* cp_fm_pool/fm_pool_create_matrix [1.0] *
!!
!!   NAME
!!     fm_pool_create_matrix
!!
!!   FUNCTION
!!     returns a matrix, allocating it if none is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool from where you get the matrix
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION fm_pool_create_matrix(pool, error) RESULT(res)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    TYPE(cp_full_matrix_type), POINTER :: res

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_create_matrix',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(pool%cache)) THEN
          res => cp_sll_fm_get_first_el(pool%cache)
          CALL cp_sll_fm_rm_first_el(pool%cache)
       ELSE
          NULLIFY(res)
          CALL cp_fm_create(res,matrix_struct=pool%matrix_struct,&
               error=error)
       END IF
    END IF
    CPPostcondition(ASSOCIATED(res),cp_failure_level,routineP,error,failure)
    CPPostcondition(res%ref_count==1,cp_warning_level,routineP,error,failure)
  END FUNCTION fm_pool_create_matrix
!***************************************************************************

!!****f* cp_fm_pool/fm_pool_give_back_matrix [1.0] *
!!
!!   NAME
!!     fm_pool_give_back_matrix
!!
!!   FUNCTION
!!     returns the matrix to the pool
!!
!!   NOTES
!!     transfere the ownership of the matrix to the pool
!!     (it is as if you had called cp_fm_release)
!!
!!   INPUTS
!!     - pool: the pool where to cache the matrix
!!     - matrix: the matrix to give back
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE fm_pool_give_back_matrix(pool, matrix, error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_full_matrix_type), POINTER :: matrix

    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_give_back_matrix',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_assert(pool%matrix_struct%id_nr==matrix%matrix_struct%id_nr,&
            cp_failure_level, cp_assertion_failed, routineP,&
            "pool cannot reuse matrixes with another structure "//&
            CPSourceFileRef,&
            error=error,failure=failure)
    END IF

    IF (.NOT. failure) THEN
       CPPreconditionNoFail(matrix%ref_count==1,cp_warning_level,routineP,error)
       CALL cp_sll_fm_insert_el(pool%cache, element=matrix)
    END IF
    NULLIFY(matrix)
  END SUBROUTINE fm_pool_give_back_matrix
!***************************************************************************

END MODULE cp_fm_pool
