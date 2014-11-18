!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief pool for for elements that are retained and released
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_fm_pool_types
  USE cp_fm_struct,                    ONLY: cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_linked_list_fm,               ONLY: cp_sll_fm_dealloc,&
                                             cp_sll_fm_get_first_el,&
                                             cp_sll_fm_insert_el,&
                                             cp_sll_fm_insert_ordered,&
                                             cp_sll_fm_next,&
                                             cp_sll_fm_rm_first_el,&
                                             cp_sll_fm_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_pool_types'
  INTEGER, SAVE, PRIVATE :: last_fm_pool_id_nr=0

  PUBLIC :: cp_fm_pool_type, cp_fm_pool_p_type
  PUBLIC :: fm_pool_create, fm_pool_retain,&
            fm_pool_release,&
            fm_pool_create_fm, fm_pool_give_back_fm,&
            fm_pool_get_el_struct
  PUBLIC :: fm_pools_copy, fm_pools_dealloc,&
            fm_pools_create_fm_vect,&
            fm_pools_give_back_fm_vect
!***

! *****************************************************************************
!> \brief represent a pool of elements with the same structure
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param el_struct the structure of the elements stored in this pool
!> \param cache linked list with the elements in the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_fm_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE(cp_fm_struct_type), POINTER :: el_struct
     
     TYPE(cp_sll_fm_type), POINTER :: cache
  END TYPE cp_fm_pool_type

! *****************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_fm_pool_p_type
     TYPE(cp_fm_pool_type), POINTER :: pool
  END TYPE cp_fm_pool_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a pool of elements
!> \param pool the pool to create
!> \param el_struct the structure of the elements that are stored in
!>        this pool
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_create(pool, el_struct,&
        error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_struct_type), POINTER         :: el_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(pool, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pool%el_struct=> el_struct
       CALL cp_fm_struct_retain(pool%el_struct,error=error)
       last_fm_pool_id_nr=last_fm_pool_id_nr+1
       pool%id_nr=last_fm_pool_id_nr
       pool%ref_count=1
       NULLIFY(pool%cache)
       
    END IF
  END SUBROUTINE fm_pool_create

! *****************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_retain(pool,error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

       pool%ref_count=pool%ref_count+1
    END IF
  END SUBROUTINE fm_pool_retain

! *****************************************************************************
!> \brief deallocates all the cached elements
!> \param pool the pool to flush
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_flush_cache(pool,error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_flush_cache', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_fm_type), POINTER                :: el_att
    TYPE(cp_sll_fm_type), POINTER            :: iterator

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       iterator => pool%cache
       DO
          IF (.NOT.cp_sll_fm_next(iterator,el_att=el_att,error=error)) EXIT
          CALL cp_fm_release(el_att,error=error)
       END DO
       CALL cp_sll_fm_dealloc(pool%cache,error=error)
    END IF
  END SUBROUTINE fm_pool_flush_cache

! *****************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_release(pool,error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(pool)) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
       pool%ref_count=pool%ref_count-1
       IF (pool%ref_count==0) THEN
          pool%ref_count=1
          CALL fm_pool_flush_cache(pool,error=error)
          CALL cp_fm_struct_release(pool%el_struct,error=error)
          pool%ref_count=0

          DEALLOCATE(pool,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE fm_pool_release

! *****************************************************************************
!> \brief returns an element, allocating it if none is in the pool
!> \param pool the pool from where you get the element
!> \param element will contain the new element
!>\param name the name for the new matrix (optional)
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_create_fm(pool, element,&
       name, error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_type), POINTER                :: element
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_create_fm', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(pool%cache)) THEN
          element => cp_sll_fm_get_first_el(pool%cache,error=error)
          CALL cp_sll_fm_rm_first_el(pool%cache,error=error)
          
       ELSE
          NULLIFY(element)
          CALL cp_fm_create(element,matrix_struct=pool%el_struct,error=error)
       END IF
       
    IF (PRESENT(name)) THEN
      element%name=name
      element%print_count=0
    ELSE
      element%name="tmp-"//TRIM(ADJUSTL(cp_to_string(element%id_nr)))
      element%print_count=0
      ! guarantee output unicity?
    END IF
    
    END IF
    CPPostcondition(ASSOCIATED(element),cp_failure_level,routineP,error,failure)
    CPPostcondition(element%ref_count==1,cp_warning_level,routineP,error,failure)
  END SUBROUTINE fm_pool_create_fm

! *****************************************************************************
!> \brief returns a vector of elements, allocating them, if none
!>      is in the pool
!> \param pool the pool from where you get the elements
!> \param n_elements the number of elements in the vector
!> \param elements will contain the new vector
!> \param name the name for the new matrixes (optional)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE fm_pool_create_fm_vect(pool, n_elements,&
       elements, name, error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    INTEGER, INTENT(in)                      :: n_elements
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: elements
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_create_fm_vect', &
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
          NULLIFY(elements(i)%matrix)
          IF (PRESENT(name)) THEN
         CALL fm_pool_create_fm(pool,elements(i)%matrix,&
           name=name//"-"//ADJUSTL(cp_to_string(i)),&
           error=error)
      ELSE
         CALL fm_pool_create_fm(pool,elements(i)%matrix,&
           error=error)
      END IF

       END DO
       
    END IF
  END SUBROUTINE fm_pool_create_fm_vect

! *****************************************************************************
!> \brief returns the element to the pool
!> \param pool the pool where to cache the element
!> \param element the element to give back
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      transfers the ownership of the element to the pool
!>      (it is as if you had called cp_fm_release)
!>      Accept give_backs of non associated elements?
! *****************************************************************************
  SUBROUTINE fm_pool_give_back_fm(pool, element, error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_type), POINTER                :: element
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_give_back_fm', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: did_insert, failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(element),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
          CALL cp_assert(pool%el_struct%id_nr==element%matrix_struct%id_nr,&
            cp_failure_level, cp_assertion_failed, routineP,&
            "pool cannot reuse matrixes with another structure "//&
CPSourceFileRef,&
            error=error,failure=failure)
    END IF

    IF (.NOT. failure) THEN
       CPPreconditionNoFail(element%ref_count==1,cp_warning_level,routineP,error)
       IF (cp_debug) THEN
          CALL cp_sll_fm_insert_ordered(pool%cache, el=element,&
               insert_equals=.FALSE., did_insert=did_insert, error=error)
          CPPostcondition(did_insert,cp_failure_level,routineP,error,failure)
       ELSE
          CALL cp_sll_fm_insert_el(pool%cache, el=element,error=error)
       END IF
    END IF
    NULLIFY(element)
  END SUBROUTINE fm_pool_give_back_fm

! *****************************************************************************
!> \brief returns the elements to the pool
!> \param pool the pool where to cache the elements
!> \param elements the elements to give back
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      transfers the ownership of the elements to the pool
!>      and deallocates the array of elements
!>      (it is as if you had called cp_fm_vect_dealloc)
!>      Accept give_backs of non associated arrays?
! *****************************************************************************
  SUBROUTINE fm_pool_give_back_fm_vect(pool, elements, error)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: elements
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_give_back_fm_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(elements),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i=1,SIZE(elements)
          IF (ASSOCIATED(elements(i)%matrix)) THEN
             CALL fm_pool_give_back_fm(pool,&
                  elements(i)%matrix,&
                  error=error)
          END IF
       END DO
       DEALLOCATE(elements,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(elements)
  END SUBROUTINE fm_pool_give_back_fm_vect

! *****************************************************************************
!> \brief returns the structure of the elements in this pool
!> \param pool the pool you are interested in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \retval res ...
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION fm_pool_get_el_struct(pool,error) RESULT(res)
    TYPE(cp_fm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(cp_fm_struct_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pool_get_el_struct', &
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
END FUNCTION fm_pool_get_el_struct

!================== pools ================

! *****************************************************************************
!> \brief shallow copy of an array of pools (retains each pool)
!> \param source_pools the pools to copy
!> \param target_pools will contains the new pools
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_copy(source_pools, target_pools, error)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_copy', &
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
        CALL fm_pool_retain(source_pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE fm_pools_copy

! *****************************************************************************
!> \brief deallocate an array of pools (releasing each pool)
!> \param pools the pools to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_dealloc(pools, error)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL fm_pool_release(pools(i)%pool,error=error)
     END DO
     DEALLOCATE(pools,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE fm_pools_dealloc

! *****************************************************************************
!> \brief deallocates all the cached elements of an array of pools
!> \param pools the pools to flush
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_flush_cache(pools, error)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        CALL fm_pool_flush_cache(pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE fm_pools_flush_cache

! *****************************************************************************
!> \brief Returns a vector with an element from each pool
!> \param pools the pools to create the elements from
!> \param elements will contain the vector of elements
!> \param name the name for the new matrixes (optional)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE fm_pools_create_fm_vect(pools,elements,&
     name,error)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: elements
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_create_fm_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(cp_fm_pool_type), POINTER           :: pool

  failure=.FALSE.
  NULLIFY(pool)

  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(elements(SIZE(pools)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        NULLIFY(elements(i)%matrix)
        pool => pools(i)%pool
        IF (PRESENT(name)) THEN
         CALL fm_pool_create_fm(pool,elements(i)%matrix,&
           name=name//"-"//ADJUSTL(cp_to_string(i)),&
           error=error)
      ELSE
         CALL fm_pool_create_fm(pool,elements(i)%matrix,&
           error=error)
      END IF

     END DO
     
  END IF
END SUBROUTINE fm_pools_create_fm_vect

! *****************************************************************************
!> \brief returns a vector to the pools. The vector is deallocated
!>      (like cp_fm_vect_dealloc)
!> \param pools the pool where to give back the vector
!> \param elements the vector of elements to give back
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      accept unassociated vect?
! *****************************************************************************
SUBROUTINE fm_pools_give_back_fm_vect(pools,elements,error)
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: elements
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fm_pools_give_back_fm_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(elements),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(pools)==SIZE(elements),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL fm_pool_give_back_fm(pools(i)%pool,&
             elements(i)%matrix,error=error)
     END DO
     DEALLOCATE(elements, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
  NULLIFY(elements)
END SUBROUTINE fm_pools_give_back_fm_vect

END MODULE cp_fm_pool_types
