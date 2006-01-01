
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_sm_pool_types [1.0] *
!!
!!   NAME
!!     cp_sm_pool_types
!!
!!   FUNCTION
!!     pool for sparse matrixes, quick port of cp_fm_pool_types, will be
!!     cleaned with the new structure of sm matrixes
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_sm_pool_types
  USE cp_linked_list_sm,               ONLY: cp_sll_sm_dealloc,&
                                             cp_sll_sm_get_first_el,&
                                             cp_sll_sm_insert_el,&
                                             cp_sll_sm_insert_ordered,&
                                             cp_sll_sm_next,&
                                             cp_sll_sm_rm_first_el,&
                                             cp_sll_sm_type
  USE cp_sm_struct,                    ONLY: cp_sm_node_type,&
                                             cp_sm_struct_type,&
                                             sm_struct_get,&
                                             sm_struct_release,&
                                             sm_struct_retain
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE kinds,                           ONLY: dp
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             cp_sm_set,&
                                             deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_sm_pool_types'
  INTEGER, SAVE, PRIVATE :: last_sm_pool_id_nr=0

  PUBLIC :: cp_sm_pool_type, cp_sm_pool_p_type
  PUBLIC :: sm_pool_create, sm_pool_retain, sm_pool_release,&
       sm_pool_create_matrix, sm_pool_give_back_matrix,&
       sm_pool_get_mstruct, sm_pool_create_matrix_vect,&
       sm_pool_give_back_matrix_vect
  PUBLIC :: sm_pools_copy, sm_pools_dealloc, sm_pools_flush_cache,&
       sm_pools_create_matrix_vect, sm_pools_give_back_matrix_vect
!***
!****************************************************************************

!!****s* cp_sm_pool_types/cp_sm_pool_type [1.0] *
!!
!!   NAME
!!     cp_sm_pool_type
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
  TYPE cp_sm_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE(cp_sm_struct_type), POINTER :: matrix_struct
     TYPE(cp_sll_sm_type), POINTER :: cache
  END TYPE cp_sm_pool_type
!!***
!****************************************************************************


!!****s* cp_sm_pool_types/cp_sm_pool_p_type [1.0] *
!!
!!   NAME
!!     cp_sm_pool_p_type
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
  TYPE cp_sm_pool_p_type
     TYPE(cp_sm_pool_type), POINTER :: pool
  END TYPE cp_sm_pool_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* cp_sm_pool_types/sm_pool_create [1.0] *
!!
!!   NAME
!!     sm_pool_create
!!
!!   SYNOPSIS
!!     Subroutine sm_pool_create(pool, matrix_struct, error)
!!       Type(cp_sm_pool_type), Pointer:: pool
!!       Type(cp_sm_struct_type), Pointer:: matrix_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pool_create
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
  SUBROUTINE sm_pool_create(pool, matrix_struct, error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_sm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(pool, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pool%matrix_struct=> matrix_struct
       CALL sm_struct_retain(pool%matrix_struct,error=error)
       last_sm_pool_id_nr=last_sm_pool_id_nr+1
       pool%id_nr=last_sm_pool_id_nr
       pool%ref_count=1
       NULLIFY(pool%cache)
    END IF
  END SUBROUTINE sm_pool_create
!***************************************************************************

!!****f* cp_sm_pool_types/cp_sm_pool_retain [1.0] *
!!
!!   NAME
!!     cp_sm_pool_retain
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
  SUBROUTINE sm_pool_retain(pool,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

       pool%ref_count=pool%ref_count+1
    END IF
  END SUBROUTINE sm_pool_retain
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pool_flush_cache [1.0] *
!!
!!   NAME
!!     sm_pool_flush_cache
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
  SUBROUTINE sm_pool_flush_cache(pool,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_flush_cache', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_sll_sm_type), POINTER            :: iterator
    TYPE(real_matrix_type), POINTER          :: m_att

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       iterator => pool%cache
       DO
          IF (.NOT.cp_sll_sm_next(iterator,el_att=m_att)) EXIT
!FM          CALL cp_sm_release(m_att,error=error)
          CALL deallocate_matrix(m_att)
       END DO
       CALL cp_sll_sm_dealloc(pool%cache,error=error)
    END IF
  END SUBROUTINE sm_pool_flush_cache
!***************************************************************************

!!****f* cp_sm_pool_types/cp_sm_pool_release [1.0] *
!!
!!   NAME
!!     cp_sm_pool_release
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
  SUBROUTINE sm_pool_release(pool,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_pool_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(pool)) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
       pool%ref_count=pool%ref_count-1
       IF (pool%ref_count==0) THEN

          pool%ref_count=1
          CALL sm_pool_flush_cache(pool,error=error)
          CALL sm_struct_release(pool%matrix_struct,error=error)
          pool%ref_count=0

          DEALLOCATE(pool,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE sm_pool_release
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pool_create_matrix [1.0] *
!!
!!   NAME
!!     sm_pool_create_matrix
!!
!!   SYNOPSIS
!!     Subroutine sm_pool_create_matrix(pool, matrix, name, error)
!!       Type(cp_sm_pool_type), Pointer:: pool
!!       Type(real_matrix_type), Pointer:: matrix
!!       Character(Len=*), Optional, Intent (IN):: name
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pool_create_matrix
!!
!!   FUNCTION
!!     returns a matrix, allocating it if none is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool from where you get the matrix
!!     - matrix: will contain the new matrix
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
  SUBROUTINE sm_pool_create_matrix(pool, matrix,name,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_create_matrix', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=40)                        :: symmetry
    INTEGER                                  :: n_blocks_col, n_blocks_row, &
                                                n_cols, n_rows
    INTEGER, DIMENSION(:), POINTER           :: first_col, first_row, &
                                                last_col, last_row
    INTEGER, SAVE                            :: last_m_id = 0
    LOGICAL                                  :: failure
    TYPE(cp_sm_node_type), POINTER           :: next_row, pos_att
    TYPE(distribution_2d_type), POINTER      :: distribution_2d

    failure=.FALSE.
    NULLIFY(first_row, first_col, last_row, last_col, next_row, pos_att,&
         distribution_2d)

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (ASSOCIATED(pool%cache)) THEN
          matrix => cp_sll_sm_get_first_el(pool%cache)
          CALL cp_sll_sm_rm_first_el(pool%cache)
       ELSE
          NULLIFY(matrix)
          CALL sm_struct_get(pool%matrix_struct,n_rows=n_rows,n_cols=n_cols,&
               n_blocks_row=n_blocks_row, n_blocks_col=n_blocks_col,&
               first_row=first_row,last_row=last_row,&
               first_col=first_col, last_col=last_col, symmetry=symmetry,&
               distribution_2d=distribution_2d,error=error)
          CALL allocate_matrix(matrix,nrow=n_rows,&
               ncol=n_cols,nblock_row=n_blocks_row,&
               nblock_col=n_blocks_col,&
               first_row=first_row,last_row=last_row,&
               first_col=first_col,last_col=last_col,&
               distribution_2d=distribution_2d,&
               sparsity_id=pool%matrix_struct%id_nr,&
               matrix_name="tmpName",matrix_symmetry=symmetry)

          ! allocate blocks
          pos_att => pool%matrix_struct%local_blocks
          DO WHILE (ASSOCIATED(pos_att))
             next_row => pos_att%next_row
             DO WHILE (ASSOCIATED(pos_att))
                IF (pos_att%row>0 .AND. pos_att%col>0) THEN
                   CALL add_block_node(matrix,block_row=pos_att%row,&
                        block_col=pos_att%col)
                END IF
                pos_att => pos_att%next_col
             END DO
             pos_att => next_row
          END DO

       END IF
       IF (PRESENT(name)) THEN
          CALL cp_sm_set(matrix,name=name,error=error)
!          matrix%print_count=0
       ELSE
          last_m_id=last_m_id+1
          CALL cp_sm_set(matrix,&
               name="tmp_matrix "//ADJUSTL(cp_to_string(last_m_id)),error=error)
       END IF
    END IF
    CPPostcondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    !CPPostcondition(matrix%ref_count==1,cp_warning_level,routineP,error,failure)
  END SUBROUTINE sm_pool_create_matrix
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pool_create_matrix_vect [1.0] *
!!
!!   NAME
!!     sm_pool_create_matrix_vect
!!
!!   SYNOPSIS
!!     Subroutine sm_pool_create_matrix_vect(pool, n_matrixes, matrixes,&
!!         name, error)
!!       Type(cp_sm_pool_type), Pointer:: pool
!!       Integer, Intent (IN):: n_matrixes
!!       Type(real_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Character(Len=*), Intent (IN), Optional:: name
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pool_create_matrix_vect
!!
!!   FUNCTION
!!     returns a vector of identical matrixes, allocating it if none
!!     is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool from where you get the matrix
!!     - n_matrixes: the number of matrixes in the vector
!!     - matrixes: will contain the new matrix vector
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
  SUBROUTINE sm_pool_create_matrix_vect(pool, n_matrixes, matrixes, name,&
       error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    INTEGER, INTENT(in)                      :: n_matrixes
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_create_matrix_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    ALLOCATE(matrixes(n_matrixes), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO i=1,n_matrixes
          NULLIFY(matrixes(i)%matrix)
          IF (PRESENT(name)) THEN
             CALL sm_pool_create_matrix(pool,matrixes(i)%matrix,&
                  name=name//"-"//ADJUSTL(cp_to_string(i)),&
                  error=error)
          ELSE
             CALL sm_pool_create_matrix(pool,matrixes(i)%matrix,error=error)
          END IF
       END DO
    END IF
  END SUBROUTINE sm_pool_create_matrix_vect
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pool_give_back_matrix [1.0] *
!!
!!   NAME
!!     sm_pool_give_back_matrix
!!
!!   SYNOPSIS
!!     Subroutine sm_pool_give_back_matrix(pool, matrix, error)
!!       Type(cp_sm_pool_type), Pointer:: pool
!!       Type(real_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pool_give_back_matrix
!!
!!   FUNCTION
!!     returns the matrix to the pool
!!
!!   NOTES
!!     transfers the ownership of the matrix to the pool
!!     (it is as if you had called cp_sm_release)
!!     Accept give_backs of non associated matrixes?
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
  SUBROUTINE sm_pool_give_back_matrix(pool, matrix, error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_give_back_matrix', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: did_insert, failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrix),cp_warning_level,routineP,error,failure)
!FM    IF (.NOT. failure) THEN
!FM       CALL cp_assert(pool%matrix_struct%id_nr==matrix%matrix_struct%id_nr,&
!FM            cp_failure_level, cp_assertion_failed, routineP,&
!FM            "pool cannot reuse matrixes with another structure "//&
!FM            CPSourceFileRef,&
!FM            error=error,failure=failure)
!FM    END IF

    IF (.NOT. failure) THEN
       CPPreconditionNoFail(matrix%ref_count==1,cp_warning_level,routineP,error)
       IF (cp_debug) THEN
          CALL cp_sll_sm_insert_ordered(pool%cache, el=matrix,&
               insert_equals=.FALSE., did_insert=did_insert)
          CPAssert(did_insert,cp_failure_level,routineP,error,failure)
       ELSE
          CALL cp_sll_sm_insert_el(pool%cache, el=matrix)
       END IF
    END IF
    NULLIFY(matrix)
  END SUBROUTINE sm_pool_give_back_matrix
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pool_give_back_matrix_vect [1.0] *
!!
!!   NAME
!!     sm_pool_give_back_matrix_vect
!!
!!   SYNOPSIS
!!     Subroutine sm_pool_give_back_matrix_vect(pool, matrixes, error)
!!       Type(cp_sm_pool_type), Pointer:: pool
!!       Type(real_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pool_give_back_matrix_vect
!!
!!   FUNCTION
!!     returns the matrix to the pool
!!
!!   NOTES
!!     transfers the ownership of the matrixes to the pool
!!     and deallocates the array of matrixes
!!     (it is as if you had called cp_sm_vect_dealloc)
!!     Accept give_backs of non associated matrixes?
!!
!!   INPUTS
!!     - pool: the pool where to cache the matrix
!!     - matrixes: the matrixes to give back
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
  SUBROUTINE sm_pool_give_back_matrix_vect(pool, matrixes, error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'sm_pool_give_back_matrix_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(matrixes),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i=1,SIZE(matrixes)
          IF (ASSOCIATED(matrixes(i)%matrix)) THEN
             CALL sm_pool_give_back_matrix(pool, matrixes(i)%matrix,&
                  error=error)
          END IF
       END DO
       DEALLOCATE(matrixes,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
    NULLIFY(matrixes)
  END SUBROUTINE sm_pool_give_back_matrix_vect
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pool_get_mstruct [1.0] *
!!
!!   NAME
!!     sm_pool_get_mstruct
!!
!!   SYNOPSIS
!!     Function sm_pool_get_mstruct(pool, error) Result(res)
!!       Type(cp_sm_struct_type), Pointer:: res
!!       Type(cp_sm_pool_type), Pointer:: pool
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function sm_pool_get_mstruct
!!
!!   FUNCTION
!!     returns the matrix structure of this pool
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
FUNCTION sm_pool_get_mstruct(pool,error) RESULT(res)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_sm_struct_type), POINTER         :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_get_mstruct', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     res => pool%matrix_struct
  ELSE
     NULLIFY(res)
  END IF
END FUNCTION sm_pool_get_mstruct
!***************************************************************************

!================== pools ================

!!****f* cp_sm_pool_types/sm_pools_copy [1.0] *
!!
!!   NAME
!!     sm_pools_copy
!!
!!   SYNOPSIS
!!     Subroutine sm_pools_copy(source_pools, target_pools, error)
!!       Type(cp_sm_pool_p_type), Dimension(:), Pointer:: source_pools,&
!!         target_pools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pools_copy
!!
!!   FUNCTION
!!     shallow copy of an array of pools (retains each pool)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - source_pools: the pools to copy
!!     - target_pools: allocated, will contains the new pools
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
SUBROUTINE sm_pools_copy(source_pools, target_pools, error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pools_copy', &
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
        CALL sm_pool_retain(source_pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE sm_pools_copy
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pools_dealloc [1.0] *
!!
!!   NAME
!!     sm_pools_dealloc
!!
!!   SYNOPSIS
!!     Subroutine sm_pools_dealloc(pools, error)
!!       Type(cp_sm_pool_p_type), Dimension(:), Pointer:: pools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pools_dealloc
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
SUBROUTINE sm_pools_dealloc(pools, error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pools_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL sm_pool_release(pools(i)%pool,error=error)
     END DO
     DEALLOCATE(pools,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE sm_pools_dealloc
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pools_flush_cache [1.0] *
!!
!!   NAME
!!     sm_pools_flush_cache
!!
!!   SYNOPSIS
!!     Subroutine sm_pools_flush_cache(pools, error)
!!       Type(cp_sm_pool_p_type), Dimension(:), Pointer:: pools
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pools_flush_cache
!!
!!   FUNCTION
!!     deallocates all the cached matrixes of an array of pools
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
SUBROUTINE sm_pools_flush_cache(pools, error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pools_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        CALL sm_pool_flush_cache(pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE sm_pools_flush_cache
!***************************************************************************

!!****f* cp_sm_pool_types/sm_pools_create_matrix_vect [1.0] *
!!
!!   NAME
!!     sm_pools_create_matrix_vect
!!
!!   SYNOPSIS
!!     Subroutine sm_pools_create_matrix_vect(pools, matrixes, name, error)
!!       Type(cp_sm_pool_p_type), Dimension(:), Pointer:: pools
!!       Type(real_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Character(Len=*), Optional, Intent (IN):: name
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pools_create_matrix_vect
!!
!!   FUNCTION
!!     Returns a matrix vector with a matrix from each pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools to create the matrix from
!!     - matrixes: will contain the vector of matrixes
!!     - name: optionally a name for the matrixes
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
SUBROUTINE sm_pools_create_matrix_vect(pools,matrixes,name,error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pools_create_matrix_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(matrixes(SIZE(pools)),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     DO i=1,SIZE(pools)
        NULLIFY(matrixes(i)%matrix)
        IF (PRESENT(name)) THEN
           CALL sm_pool_create_matrix(pools(i)%pool,matrixes(i)%matrix,&
                name=name//'-'//ADJUSTL(cp_to_string(i)),error=error)
        ELSE
           CALL sm_pool_create_matrix(pools(i)%pool,matrixes(i)%matrix,&
                error=error)
        END IF
     END DO
  END IF
END SUBROUTINE sm_pools_create_matrix_vect
!***************************************************************************

!!****f* cp_sm_pool_p_types/sm_pools_give_back_matrix_vect [1.0] *
!!
!!   NAME
!!     sm_pools_give_back_matrix_vect
!!
!!   SYNOPSIS
!!     Subroutine sm_pools_give_back_matrix_vect(pools, matrixes, error)
!!       Type(cp_sm_pool_p_type), Dimension(:), Pointer:: pools
!!       Type(real_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_pools_give_back_matrix_vect
!!
!!   FUNCTION
!!     returns a matrix vector to the pools. The vector is deallocated
!!     (like cp_sm_vect_dealloc)
!!
!!   NOTES
!!     accept unassociated vect?
!!
!!   INPUTS
!!     - pools: the pool where to give back the matrix vector
!!     - matrixes: the matrixes to give back
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
SUBROUTINE sm_pools_give_back_matrix_vect(pools,matrixes,error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'sm_pools_give_back_matrix_vect', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrixes),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(matrixes)==SIZE(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL sm_pool_give_back_matrix(pools(i)%pool,&
             matrixes(i)%matrix,error=error)
     END DO
  END IF
  DEALLOCATE(matrixes,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  NULLIFY(matrixes)
END SUBROUTINE sm_pools_give_back_matrix_vect
!***************************************************************************

END MODULE cp_sm_pool_types
