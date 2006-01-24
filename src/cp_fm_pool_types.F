!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/cp_fm_pool_types [1.0] *
!!
!!   NAME
!!     cp_fm_pool_types
!!
!!   FUNCTION
!!     pool for for elements that are retained and released
!!
!!   NOTES
!!    ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _ 
!!   |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!!   | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!!   | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!!   |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!!                                                       |___/                                        
!!     ____ _                  ___                              _ _       _       _       
!!    / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _ 
!!   | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!!   | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!!    \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!!                                                                                  |___/ 
!!    _____ _     _       _____ _ _      _ 
!!   |_   _| |__ (_)___  |  ___(_) | ___| |
!!     | | | '_ \| / __| | |_  | | |/ _ \ |
!!     | | | | | | \__ \ |  _| | | |  __/_|
!!     |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!!   
!!     This is a template
!!
!!     **** DO NOT MODIFY THE .F FILES ****
!!     modify the .template and .instantition instead
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_pool_types
  USE cp_linked_list_fm,        ONLY: cp_sll_fm_dealloc,&
                                             cp_sll_fm_get_first_el,&
                                             cp_sll_fm_insert_el,&
                                             cp_sll_fm_next,&
                                             cp_sll_fm_rm_first_el,&
                                             cp_sll_fm_type,&
                                             cp_sll_fm_insert_ordered
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"
  USE cp_fm_struct, ONLY: cp_fm_struct_retain, cp_fm_struct_release,&
        cp_fm_struct_type
    USE cp_fm_types,  ONLY: cp_fm_type, cp_fm_p_type,&
        cp_fm_create, cp_fm_release
        
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_pool_types'
  INTEGER, SAVE, PRIVATE :: last_fm_pool_id_nr=0

  PUBLIC :: cp_fm_pool_type, cp_fm_pool_p_type
  PUBLIC :: fm_pool_create, fm_pool_retain,&
       fm_pool_release,&
       fm_pool_create_fm, fm_pool_give_back_fm,&
       fm_pool_get_el_struct, fm_pool_create_fm_vect,&
       fm_pool_give_back_fm_vect
  PUBLIC :: fm_pools_copy, fm_pools_dealloc,&
       fm_pools_flush_cache, fm_pools_create_fm_vect,&
       fm_pools_give_back_fm_vect
!***
!****************************************************************************

!!****s* cp_fm_pool_types/cp_fm_pool_type [1.0] *
!!
!!   NAME
!!     cp_fm_pool_type
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
  TYPE cp_fm_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     type(cp_fm_struct_type), POINTER :: el_struct
     
     TYPE(cp_sll_fm_type), POINTER :: cache
  END TYPE cp_fm_pool_type
!!***
!****************************************************************************


!!****s* cp_fm_pool_types/cp_fm_pool_p_type [1.0] *
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

!!****f* cp_fm_pool_types/fm_pool_create [1.0] *
!!
!!   NAME
!!     fm_pool_create
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
  SUBROUTINE fm_pool_create(pool, el_struct,&
        error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    type(cp_fm_struct_type), POINTER :: el_struct
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
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
!***************************************************************************

!!****f* cp_fm_pool_types/cp_fm_pool_retain [1.0] *
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

!!****f* cp_fm_pool_types/fm_pool_flush_cache [1.0] *
!!
!!   NAME
!!     fm_pool_flush_cache
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
  SUBROUTINE fm_pool_flush_cache(pool,error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_flush_cache',&
         routineP=moduleN//':'//routineN
    TYPE(cp_sll_fm_type), POINTER :: iterator
    type(cp_fm_type), POINTER :: el_att
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       iterator => pool%cache
       DO
          IF (.NOT.cp_sll_fm_next(iterator,el_att=el_att)) EXIT
          CALL cp_fm_release(el_att,error=error)
       END DO
       CALL cp_sll_fm_dealloc(pool%cache,error=error)
    END IF
  END SUBROUTINE fm_pool_flush_cache
!***************************************************************************

!!****f* cp_fm_pool_types/cp_fm_pool_release [1.0] *
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
    type(cp_fm_type), POINTER :: el_att
    INTEGER :: stat
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
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pool_create_fm [1.0] *
!!
!!   NAME
!!     fm_pool_create_fm
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
!!     - name: the name for the new matrix (optional)
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
  SUBROUTINE fm_pool_create_fm(pool, element,&
       name, error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    type(cp_fm_type), POINTER :: element
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    character(len=*), intent(in), optional :: name

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_create_fm',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(pool%cache)) THEN
          element => cp_sll_fm_get_first_el(pool%cache)
          CALL cp_sll_fm_rm_first_el(pool%cache)
          
       ELSE
          NULLIFY(element)
              call cp_fm_create(element,matrix_struct=pool%el_struct,error=error)
       END IF
       
    IF (present(name)) THEN
      element%name=name
      element%print_count=0
    ELSE
      element%name="tmp-"//trim(adjustl(cp_to_string(element%id_nr)))
      element%print_count=0
      ! guarantee output unicity?
    END IF
    
    END IF
    CPPostcondition(ASSOCIATED(element),cp_failure_level,routineP,error,failure)
    CPPostcondition(element%ref_count==1,cp_warning_level,routineP,error,failure)
  END SUBROUTINE fm_pool_create_fm
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pool_create_fm_vect [1.0] *
!!
!!   NAME
!!     fm_pool_create_fm_vect
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
!!     - name: the name for the new matrixes (optional)
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
  SUBROUTINE fm_pool_create_fm_vect(pool, n_elements,&
       elements, name, error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    INTEGER, INTENT(in) :: n_elements
    type(cp_fm_p_type), DIMENSION(:), POINTER :: elements
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    character(len=*),intent(in),optional :: name

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_create_fm_vect',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat,i

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    ALLOCATE(elements(n_elements), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO i=1,n_elements
          NULLIFY(elements(i)%matrix)
          IF (PRESENT(name)) THEN
         call fm_pool_create_fm(pool,elements(i)%matrix,&
           name=name//"-"//ADJUSTL(cp_to_string(i)),&
           error=error)
      ELSE
         call fm_pool_create_fm(pool,elements(i)%matrix,&
           error=error)
      END IF

       END DO
       
    END IF
  END SUBROUTINE fm_pool_create_fm_vect
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pool_give_back_fm [1.0] *
!!
!!   NAME
!!     fm_pool_give_back_fm
!!
!!   FUNCTION
!!     returns the element to the pool
!!
!!   NOTES
!!     transfers the ownership of the element to the pool
!!     (it is as if you had called cp_fm_release)
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
  SUBROUTINE fm_pool_give_back_fm(pool, element, error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    type(cp_fm_type), POINTER :: element
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, did_insert
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_give_back_fm',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
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
          CALL cp_sll_fm_insert_el(pool%cache, el=element)
       END IF
    END IF
    NULLIFY(element)
  END SUBROUTINE fm_pool_give_back_fm
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pool_give_back_fm_vect [1.0] *
!!
!!   NAME
!!     fm_pool_give_back_fm_vect
!!
!!   FUNCTION
!!     returns the elements to the pool
!!
!!   NOTES
!!     transfers the ownership of the elements to the pool
!!     and deallocates the array of elements
!!     (it is as if you had called cp_fm_vect_dealloc)
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
  SUBROUTINE fm_pool_give_back_fm_vect(pool, elements, error)
    TYPE(cp_fm_pool_type), POINTER :: pool
    type(cp_fm_p_type), DIMENSION(:), POINTER :: elements
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='fm_pool_give_back_elements_vect',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat,i
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
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pool_get_el_struct [1.0] *
!!
!!   NAME
!!     fm_pool_get_el_struct
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
FUNCTION fm_pool_get_el_struct(pool,error) RESULT(res)
  type(cp_fm_struct_type), POINTER :: res
  TYPE(cp_fm_pool_type), POINTER :: pool
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='fm_pool_get_el_struct',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res => pool%el_struct
  ELSE
     NULLIFY(res)
  END IF
END FUNCTION fm_pool_get_el_struct
!***************************************************************************

!================== pools ================

!!****f* cp_fm_pool_types/fm_pools_copy [1.0] *
!!
!!   NAME
!!     fm_pools_copy
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
SUBROUTINE fm_pools_copy(source_pools, target_pools, error)
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: source_pools,&
       target_pools
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: stat, i
  CHARACTER(len=*), PARAMETER :: routineN='fm_pools_copy',&
        routineP=moduleN//':'//routineN
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
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pools_dealloc [1.0] *
!!
!!   NAME
!!     fm_pools_dealloc
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
SUBROUTINE fm_pools_dealloc(pools, error)
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: pools
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)   :: error
  
  LOGICAL :: failure
  INTEGER :: stat, i
  CHARACTER(len=*), PARAMETER :: routineN='fm_pools_dealloc',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL fm_pool_release(pools(i)%pool,error=error)
     END DO
     DEALLOCATE(pools,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE fm_pools_dealloc
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pools_flush_cache [1.0] *
!!
!!   NAME
!!     fm_pools_flush_cache
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
SUBROUTINE fm_pools_flush_cache(pools, error)
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: pools
  TYPE(cp_error_type), OPTIONAL, INTENT(inout)   :: error
  
  LOGICAL :: failure
  INTEGER :: stat,i
  CHARACTER(len=*), PARAMETER :: routineN='fm_pools_flush_cache',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        CALL fm_pool_flush_cache(pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE fm_pools_flush_cache
!***************************************************************************

!!****f* cp_fm_pool_types/fm_pools_create_fm_vect [1.0] *
!!
!!   NAME
!!     fm_pools_create_fm_vect
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
!!     - name: the name for the new matrixes (optional)
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
SUBROUTINE fm_pools_create_fm_vect(pools,elements,&
     name,error)
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: pools
  type(cp_fm_p_type), DIMENSION(:), POINTER :: elements
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  character(len=*),intent(in),optional :: name
  
  LOGICAL :: failure
  INTEGER :: stat,i
  CHARACTER(len=*), PARAMETER :: routineN='fm_pools_create_fm_vect',&
        routineP=moduleN//':'//routineN
  TYPE(cp_fm_pool_type), POINTER :: pool

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
         call fm_pool_create_fm(pool,elements(i)%matrix,&
           name=name//"-"//ADJUSTL(cp_to_string(i)),&
           error=error)
      ELSE
         call fm_pool_create_fm(pool,elements(i)%matrix,&
           error=error)
      END IF

     END DO
     
  END IF
END SUBROUTINE fm_pools_create_fm_vect
!***************************************************************************

!!****f* cp_fm_pool_p_types/fm_pools_give_back_fm_vect [1.0] *
!!
!!   NAME
!!     fm_pools_give_back_fm_vect
!!
!!   FUNCTION
!!     returns a vector to the pools. The vector is deallocated
!!     (like cp_fm_vect_dealloc)
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
SUBROUTINE fm_pools_give_back_fm_vect(pools,elements,error)
  TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: pools
  type(cp_fm_p_type), DIMENSION(:), POINTER :: elements
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  INTEGER :: i, stat
  CHARACTER(len=*), PARAMETER :: routineN='fm_pools_give_back_fm_vect',&
        routineP=moduleN//':'//routineN
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
  nullify(elements)
END SUBROUTINE fm_pools_give_back_fm_vect
!***************************************************************************

END MODULE cp_fm_pool_types
