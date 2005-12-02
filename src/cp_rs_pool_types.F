!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_rs_pool_types [1.0] *
!!
!!   NAME
!!     cp_rs_pool_types
!!
!!   FUNCTION
!!     pool for for elements that are retained and released
!!
!!   NOTES
!!     This is a template
!!
!!     **** DO NOT MODIFY THE .F FILES ****
!!     modify the .template and .instantition
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_rs_pool_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_linked_list_rs,               ONLY: cp_sll_rs_dealloc,&
                                             cp_sll_rs_get_first_el,&
                                             cp_sll_rs_insert_el,&
                                             cp_sll_rs_insert_ordered,&
                                             cp_sll_rs_next,&
                                             cp_sll_rs_rm_first_el,&
                                             cp_sll_rs_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release,&
                                             pw_grid_retain
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type,&
                                             rs_grid_create,&
                                             rs_grid_release
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_rs_pool_types'
  INTEGER, SAVE, PRIVATE :: last_rs_pool_id_nr=0

  PUBLIC :: cp_rs_pool_type, cp_rs_pool_p_type
  PUBLIC :: rs_pool_create, rs_pool_retain,&
       rs_pool_release,&
       rs_pool_create_rs, rs_pool_give_back_rs,&
       rs_pool_get_el_struct, rs_pool_create_rs_vect,&
       rs_pool_give_back_rs_vect
  PUBLIC :: rs_pools_copy, rs_pools_dealloc,&
       rs_pools_flush_cache, rs_pools_create_rs_vect,&
       rs_pools_give_back_rs_vect
!***
!****************************************************************************

!!****s* cp_rs_pool_types/cp_rs_pool_type [1.0] *
!!
!!   NAME
!!     cp_rs_pool_type
!!
!!   FUNCTION
!!     represent a pool of elements with the same structure
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (see /cp2k/doc/ReferenceCounting.html)
!!     - el_struct: the structure of the elements stored in this pool
!!     - cache: linked list with the elements in the pool
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_rs_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE(pw_grid_type), POINTER :: el_struct
     INTEGER :: nsmax
     TYPE(cp_sll_rs_type), POINTER :: cache
  END TYPE cp_rs_pool_type
!!***
!****************************************************************************


!!****s* cp_rs_pool_types/cp_rs_pool_p_type [1.0] *
!!
!!   NAME
!!     cp_rs_pool_p_type
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
  TYPE cp_rs_pool_p_type
     TYPE(cp_rs_pool_type), POINTER :: pool
  END TYPE cp_rs_pool_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* cp_rs_pool_types/rs_pool_create [1.0] *
!!
!!   NAME
!!     rs_pool_create
!!
!!   FUNCTION
!!     creates a pool of elements
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to create
!!     - el_struct: the structure of the elements that are stored in
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
  SUBROUTINE rs_pool_create(pool, el_struct,&
       nsmax, error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(pw_grid_type), POINTER              :: el_struct
    INTEGER, INTENT(IN)                      :: nsmax
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(pool, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pool%el_struct=> el_struct
       CALL pw_grid_retain(pool%el_struct,error=error)
       last_rs_pool_id_nr=last_rs_pool_id_nr+1
       pool%id_nr=last_rs_pool_id_nr
       pool%ref_count=1
       NULLIFY(pool%cache)
       pool%nsmax=nsmax
    END IF
  END SUBROUTINE rs_pool_create
!***************************************************************************

!!****f* cp_rs_pool_types/cp_rs_pool_retain [1.0] *
!!
!!   NAME
!!     cp_rs_pool_retain
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
  SUBROUTINE rs_pool_retain(pool,error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

       pool%ref_count=pool%ref_count+1
    END IF
  END SUBROUTINE rs_pool_retain
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pool_flush_cache [1.0] *
!!
!!   NAME
!!     rs_pool_flush_cache
!!
!!   FUNCTION
!!     deallocates all the cached elements
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
  SUBROUTINE rs_pool_flush_cache(pool,error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_flush_cache', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_sll_rs_type), POINTER            :: iterator
    TYPE(realspace_grid_type), POINTER       :: el_att

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       iterator => pool%cache
       DO
          IF (.NOT.cp_sll_rs_next(iterator,el_att=el_att)) EXIT
          CALL rs_grid_release(el_att,error=error)
       END DO
       CALL cp_sll_rs_dealloc(pool%cache,error=error)
    END IF
  END SUBROUTINE rs_pool_flush_cache
!***************************************************************************

!!****f* cp_rs_pool_types/cp_rs_pool_release [1.0] *
!!
!!   NAME
!!     cp_rs_pool_release
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
  SUBROUTINE rs_pool_release(pool,error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_rs_pool_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(pool)) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
       pool%ref_count=pool%ref_count-1
       IF (pool%ref_count==0) THEN
          pool%ref_count=1
          CALL rs_pool_flush_cache(pool,error=error)
          CALL pw_grid_release(pool%el_struct,error=error)
          pool%ref_count=0

          DEALLOCATE(pool,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE rs_pool_release
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pool_create_rs [1.0] *
!!
!!   NAME
!!     rs_pool_create_rs
!!
!!   FUNCTION
!!     returns an element, allocating it if none is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool from where you get the element
!!     - element: will contain the new element
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
  SUBROUTINE rs_pool_create_rs(pool, element,&
        error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(realspace_grid_type), POINTER       :: element
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_create_rs', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(pool%cache)) THEN
          element => cp_sll_rs_get_first_el(pool%cache)
          CALL cp_sll_rs_rm_first_el(pool%cache)
          
       ELSE
          NULLIFY(element)
          CALL rs_grid_create(element,pw_grid=pool%el_struct,&
         nsmax=pool%nsmax,error=error)
       END IF
       
    END IF
    CPPostcondition(ASSOCIATED(element),cp_failure_level,routineP,error,failure)
    CPPostcondition(element%ref_count==1,cp_warning_level,routineP,error,failure)
  END SUBROUTINE rs_pool_create_rs
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pool_create_rs_vect [1.0] *
!!
!!   NAME
!!     rs_pool_create_rs_vect
!!
!!   FUNCTION
!!     returns a vector of elements, allocating them, if none
!!     is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool from where you get the elements
!!     - n_elements: the number of elements in the vector
!!     - elements: will contain the new vector
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
  SUBROUTINE rs_pool_create_rs_vect(pool, n_elements,&
       elements,  error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    INTEGER, INTENT(in)                      :: n_elements
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: elements
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_create_rs_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    ALLOCATE(elements(n_elements), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO i=1,n_elements
          NULLIFY(elements(i)%rs_grid)
          CALL rs_pool_create_rs(pool,elements(i)%rs_grid,error=error)
!          CALL rs_pool_create_rs(pool,elements(i)%rs_grid,&
!                error=error)
       END DO
       
    END IF
  END SUBROUTINE rs_pool_create_rs_vect
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pool_give_back_rs [1.0] *
!!
!!   NAME
!!     rs_pool_give_back_rs
!!
!!   FUNCTION
!!     returns the element to the pool
!!
!!   NOTES
!!     transfers the ownership of the element to the pool
!!     (it is as if you had called cp_rs_release)
!!     Accept give_backs of non associated elements?
!!
!!   INPUTS
!!     - pool: the pool where to cache the element
!!     - element: the element to give back
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
  SUBROUTINE rs_pool_give_back_rs(pool, element, error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(realspace_grid_type), POINTER       :: element
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_give_back_rs', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: did_insert, failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(element),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
          CALL cp_assert(pool%el_struct%id_nr==element%grid_id,&
            cp_failure_level, cp_assertion_failed, routineP,&
            "pool cannot reuse rs grids of another pw_grid "//&
            CPSourceFileRef,&
            error=error,failure=failure)
    END IF

    IF (.NOT. failure) THEN
       CPPreconditionNoFail(element%ref_count==1,cp_warning_level,routineP,error)
       IF (cp_debug) THEN
          CALL cp_sll_rs_insert_ordered(pool%cache, el=element,&
               insert_equals=.FALSE., did_insert=did_insert, error=error)
          CPPostcondition(did_insert,cp_failure_level,routineP,error,failure)
       ELSE
          CALL cp_sll_rs_insert_el(pool%cache, el=element)
       END IF
    END IF
    NULLIFY(element)
  END SUBROUTINE rs_pool_give_back_rs
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pool_give_back_rs_vect [1.0] *
!!
!!   NAME
!!     rs_pool_give_back_rs_vect
!!
!!   FUNCTION
!!     returns the elements to the pool
!!
!!   NOTES
!!     transfers the ownership of the elements to the pool
!!     and deallocates the array of elements
!!     (it is as if you had called cp_rs_vect_dealloc)
!!     Accept give_backs of non associated arrays?
!!
!!   INPUTS
!!     - pool: the pool where to cache the elements
!!     - elements: the elements to give back
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
  SUBROUTINE rs_pool_give_back_rs_vect(pool, elements, error)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: elements
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'rs_pool_give_back_elements_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elements),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i=1,SIZE(elements)
          IF (ASSOCIATED(elements(i)%rs_grid)) THEN
             CALL rs_pool_give_back_rs(pool,&
                  elements(i)%rs_grid,&
                  error=error)
          END IF
       END DO
       DEALLOCATE(elements,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(elements)
  END SUBROUTINE rs_pool_give_back_rs_vect
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pool_get_el_struct [1.0] *
!!
!!   NAME
!!     rs_pool_get_el_struct
!!
!!   FUNCTION
!!     returns the structure of the elements in this pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool you are interested in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
FUNCTION rs_pool_get_el_struct(pool,error) RESULT(res)
    TYPE(cp_rs_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(pw_grid_type), POINTER              :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pool_get_el_struct', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res => pool%el_struct
  ELSE
     NULLIFY(res)
  END IF
END FUNCTION rs_pool_get_el_struct
!***************************************************************************

!================== pools ================

!!****f* cp_rs_pool_types/rs_pools_copy [1.0] *
!!
!!   NAME
!!     rs_pools_copy
!!
!!   FUNCTION
!!     shallow copy of an array of pools (retains each pool)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - source_pools: the pools to copy
!!     - target_pools: will contains the new pools
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_pools_copy(source_pools, target_pools, error)
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pools_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(source_pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(target_pools(SIZE(source_pools)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     DO i=1,SIZE(source_pools)
        target_pools(i)%pool => source_pools(i)%pool
        CALL rs_pool_retain(source_pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE rs_pools_copy
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pools_dealloc [1.0] *
!!
!!   NAME
!!     rs_pools_dealloc
!!
!!   FUNCTION
!!     deallocate an array of pools (releasing each pool)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_pools_dealloc(pools, error)
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pools_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL rs_pool_release(pools(i)%pool,error=error)
     END DO
     DEALLOCATE(pools,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE rs_pools_dealloc
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pools_flush_cache [1.0] *
!!
!!   NAME
!!     rs_pools_flush_cache
!!
!!   FUNCTION
!!     deallocates all the cached elements of an array of pools
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools to flush
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_pools_flush_cache(pools, error)
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pools_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        CALL rs_pool_flush_cache(pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE rs_pools_flush_cache
!***************************************************************************

!!****f* cp_rs_pool_types/rs_pools_create_rs_vect [1.0] *
!!
!!   NAME
!!     rs_pools_create_rs_vect
!!
!!   FUNCTION
!!     Returns a vector with an element from each pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools to create the elements from
!!     - elements: will contain the vector of elements
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_pools_create_rs_vect(pools,elements,&
     error)
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: elements
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pools_create_rs_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(cp_rs_pool_type), POINTER           :: pool

  failure=.FALSE.
  NULLIFY(pool)
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(elements(SIZE(pools)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        NULLIFY(elements(i)%rs_grid)
        pool => pools(i)%pool
        CALL rs_pool_create_rs(pool,elements(i)%rs_grid,error=error)
!           CALL rs_pool_create_type(realspace_grid_type)(pools(i)%pool,&
!               elements(i)%rs_grid,error=error)
     END DO
     
  END IF
END SUBROUTINE rs_pools_create_rs_vect
!***************************************************************************

!!****f* cp_rs_pool_p_types/rs_pools_give_back_rs_vect [1.0] *
!!
!!   NAME
!!     rs_pools_give_back_rs_vect
!!
!!   FUNCTION
!!     returns a vector to the pools. The vector is deallocated
!!     (like cp_rs_vect_dealloc)
!!
!!   NOTES
!!     accept unassociated vect?
!!
!!   INPUTS
!!     - pools: the pool where to give back the vector
!!     - elements: the vector of elements to give back
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_pools_give_back_rs_vect(pools,elements,error)
    TYPE(cp_rs_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: elements
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pools_give_back_rs_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(elements),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(pools)==SIZE(elements),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL rs_pool_give_back_rs(pools(i)%pool,&
             elements(i)%rs_grid,error=error)
     END DO
     DEALLOCATE(elements, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
  NULLIFY(elements)
END SUBROUTINE rs_pools_give_back_rs_vect
!***************************************************************************

END MODULE cp_rs_pool_types
