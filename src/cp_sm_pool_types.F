!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief pool for sparse matrixes, quick port of cp_fm_pool_types, will be
!>      cleaned with the new structure of sm matrixes
!> \par History
!>      03.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
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
  USE f77_blas
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             allocate_matrix,&
                                             cp_sm_set,&
                                             deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_sm_pool_types'
  INTEGER, SAVE, PRIVATE :: last_sm_pool_id_nr=0

  PUBLIC :: cp_sm_pool_type, cp_sm_pool_p_type
  PUBLIC :: sm_pool_create, sm_pool_retain, sm_pool_release,&
       sm_pool_create_matrix, sm_pool_give_back_matrix,&
       sm_pool_get_mstruct, sm_pool_create_matrix_vect,&
       sm_pool_give_back_matrix_vect
  PUBLIC :: sm_pools_copy, sm_pools_dealloc, sm_pools_flush_cache,&
       sm_pools_create_matrix_vect, sm_pools_give_back_matrix_vect
!***

! *****************************************************************************
!> \brief represent a pool of full matrixes
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param matrix_struct matrix structure of the matrixes in this pool
!> \param cache linked list with the matrixes in the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_sm_pool_type
     PRIVATE
     INTEGER :: ref_count, id_nr
     TYPE(cp_sm_struct_type), POINTER :: matrix_struct
     TYPE(cp_sll_sm_type), POINTER :: cache
  END TYPE cp_sm_pool_type

! *****************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_sm_pool_p_type
     TYPE(cp_sm_pool_type), POINTER :: pool
  END TYPE cp_sm_pool_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a pool of full matrixes
!> \param pool the pool to create
!> \param matrix_struct the structure of the matrixes that are stored in
!>        this pool
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_create(pool, matrix_struct, error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_sm_struct_type), POINTER         :: matrix_struct
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_retain(pool,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief deallocates all the cached matrixes
!> \param pool the pool to flush
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_flush_cache(pool,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error

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
          IF (.NOT.cp_sll_sm_next(iterator,el_att=m_att,error=error)) EXIT
!FM          CALL cp_sm_release(m_att,error=error)
          CALL deallocate_matrix(m_att,error=error)
       END DO
       CALL cp_sll_sm_dealloc(pool%cache,error=error)
    END IF
  END SUBROUTINE sm_pool_flush_cache

! *****************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_release(pool,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_pool_release', &
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

! *****************************************************************************
!> \brief returns a matrix, allocating it if none is in the pool
!> \param pool the pool from where you get the matrix
!> \param matrix will contain the new matrix
!> \param name the name for the new matrix (optional)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_create_matrix(pool, matrix,name,error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(real_matrix_type), POINTER          :: matrix
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

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
          matrix => cp_sll_sm_get_first_el(pool%cache,error=error)
          CALL cp_sll_sm_rm_first_el(pool%cache,error=error)
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
               matrix_name="tmpName",matrix_symmetry=symmetry,error=error)

          ! allocate blocks
          pos_att => pool%matrix_struct%local_blocks
          DO WHILE (ASSOCIATED(pos_att))
             next_row => pos_att%next_row
             DO WHILE (ASSOCIATED(pos_att))
                IF (pos_att%row>0 .AND. pos_att%col>0) THEN
                   CALL add_block_node(matrix,block_row=pos_att%row,&
                        block_col=pos_att%col,error=error)
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

! *****************************************************************************
!> \brief returns a vector of identical matrixes, allocating it if none
!>      is in the pool
!> \param pool the pool from where you get the matrix
!> \param n_matrixes the number of matrixes in the vector
!> \param matrixes will contain the new matrix vector
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_create_matrix_vect(pool, n_matrixes, matrixes, name,&
       error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    INTEGER, INTENT(in)                      :: n_matrixes
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief returns the matrix to the pool
!> \param pool the pool where to cache the matrix
!> \param matrix the matrix to give back
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      transfers the ownership of the matrix to the pool
!>      (it is as if you had called cp_sm_release)
!>      Accept give_backs of non associated matrixes?
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_give_back_matrix(pool, matrix, error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout)       :: error

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
!FMCPSourceFileRef,&
!FM            error=error,failure=failure)
!FM    END IF

    IF (.NOT. failure) THEN
       CPPreconditionNoFail(matrix%ref_count==1,cp_warning_level,routineP,error)
       IF (cp_debug) THEN
          CALL cp_sll_sm_insert_ordered(pool%cache, el=matrix,&
               insert_equals=.FALSE., did_insert=did_insert,error=error)
          CPAssert(did_insert,cp_failure_level,routineP,error,failure)
       ELSE
          CALL cp_sll_sm_insert_el(pool%cache, el=matrix,error=error)
       END IF
    END IF
    NULLIFY(matrix)
  END SUBROUTINE sm_pool_give_back_matrix

! *****************************************************************************
!> \brief returns the matrix to the pool
!> \param pool the pool where to cache the matrix
!> \param matrixes the matrixes to give back
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      transfers the ownership of the matrixes to the pool
!>      and deallocates the array of matrixes
!>      (it is as if you had called cp_sm_vect_dealloc)
!>      Accept give_backs of non associated matrixes?
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE sm_pool_give_back_matrix_vect(pool, matrixes, error)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief returns the matrix structure of this pool
!> \param pool the pool you are interested in
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION sm_pool_get_mstruct(pool,error) RESULT(res)
    TYPE(cp_sm_pool_type), POINTER           :: pool
    TYPE(cp_error_type), INTENT(inout)       :: error
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

!================== pools ================

! *****************************************************************************
!> \brief shallow copy of an array of pools (retains each pool)
!> \param source_pools the pools to copy
!> \param target_pools allocated, will contains the new pools
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE sm_pools_copy(source_pools, target_pools, error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief deallocate an array of pools (releasing each pool)
!> \param pools the pools to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE sm_pools_dealloc(pools, error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief deallocates all the cached matrixes of an array of pools
!> \param pools the pools to flush
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE sm_pools_flush_cache(pools, error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief Returns a matrix vector with a matrix from each pool
!> \param pools the pools to create the matrix from
!> \param matrixes will contain the vector of matrixes
!> \param name optionally a name for the matrixes
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE sm_pools_create_matrix_vect(pools,matrixes,name,error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

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

! *****************************************************************************
!> \brief returns a matrix vector to the pools. The vector is deallocated
!>      (like cp_sm_vect_dealloc)
!> \param pools the pool where to give back the matrix vector
!> \param matrixes the matrixes to give back
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      accept unassociated vect?
!> \par History
!>      09.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE sm_pools_give_back_matrix_vect(pools,matrixes,error)
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrixes
    TYPE(cp_error_type), INTENT(inout)       :: error

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

END MODULE cp_sm_pool_types
