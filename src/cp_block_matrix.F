!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_block_matrix [1.0] *
!!
!!   NAME
!!     cp_block_matrix
!!
!!   FUNCTION
!!     represent a sparse matrix type
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     3.2002 created [fawzi]
!!     4.2002 adapted to replace pao_matrix, moved iterators away [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_block_matrix
  USE cp_b_matrix_structure
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_precondition_failed,&
                                             cp_unimplemented_error
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: cp_matrix_block_iterator,&
                                             cp_sm_alloc_nonsparse_blocks,&
                                             cp_sm_b_i_dealloc_ref,&
                                             cp_sm_b_i_init,&
                                             cp_sm_b_i_next,&
                                             cp_sm_set_to
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sync
  USE sparse_matrix_types
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_block_matrix'

  PUBLIC :: cp_block_matrix_type, cp_sp_block_handle_type,&
       cp_block_matrix_p_type
  PUBLIC :: cp_create, cp_dealloc, cp_retain, cp_release, &
       cp_get_local_block,&
       cp_release_block, cp_set_local_block, cp_remove_block,&
       cp_get_matrix, cp_set_to, cp_get
  PUBLIC :: cp_sp_create, cp_sp_retain, cp_sp_release,&
       cp_sp_get_local_block,&
       cp_sp_block_handle_release, cp_sp_set_local_block, cp_sp_remove_block,&
       cp_sp_array_dealloc, cp_sp_get_matrix, cp_sp_array_release,&
       cp_sp_alloc_nonsparse_blocks, cp_sp_set_to, cp_sp_block_handle_init,&
       cp_sp_get, cp_sp_multiplies_blacs, cp_sp_array_create

  INTERFACE cp_init
     MODULE PROCEDURE cp_sp_block_handle_init
  END INTERFACE
  INTERFACE cp_create
     MODULE PROCEDURE cp_sp_create, cp_sp_array_create
  END INTERFACE
!FM  interface cp_dealloc
!FM     module procedure cp_sp_dealloc
!FM  end interface
  INTERFACE cp_retain
     MODULE PROCEDURE cp_sp_retain
  END INTERFACE
  INTERFACE cp_release
     MODULE PROCEDURE cp_sp_release, cp_sp_array_release
  END INTERFACE
  INTERFACE cp_get_matrix
     MODULE PROCEDURE cp_sp_get_matrix
  END INTERFACE
  INTERFACE cp_flush_cache
     MODULE PROCEDURE cp_sp_flush_cache
  END INTERFACE
  INTERFACE cp_get_local_block
     MODULE PROCEDURE cp_sp_get_local_block
  END INTERFACE
  INTERFACE cp_release_block
     MODULE PROCEDURE cp_sp_block_handle_release
  END INTERFACE
  INTERFACE cp_set_local_block
     MODULE PROCEDURE cp_sp_set_local_block
  END INTERFACE
  INTERFACE cp_remove_block
     MODULE PROCEDURE cp_sp_remove_block
  END INTERFACE
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_sp_array_dealloc
  END INTERFACE
  INTERFACE cp_set_to
     MODULE PROCEDURE cp_sp_set_to
  END INTERFACE
  INTERFACE cp_get
     MODULE PROCEDURE cp_sp_get
  END INTERFACE
!***
!****************************************************************************

!!****s* cp_block_matrix/cp_block_matrix_type [1.0] *
!!
!!   NAME
!!     cp_block_matrix_type
!!
!!   FUNCTION
!!     represent a sparse matrix
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - initialized: true if the structure was initialized (bug catcher) 
!!       (use ref_count?)
!!     - frozen: true if the matrix is frozen (no change to the values
!!       allowed).
!!     - should_dealloc_matrix: true if the matrix should be 
!!       deallocated when the structure is deallocated
!!     - b_rep_is_primary: if true the block representation is the primary
!!       representation (get works), otherwise the internal is.
!!     - b_i_rep_in_sync: block and internal representation are syncronized
!!     - cache_is_valid: if true assumes that the local non owned blocks
!!       contain an up-to-date snapshot of the corresponding block.
!!     - symmetric: if the matrix is symmetric
!!     - matrix_struct: the matrix structure of this matrix
!!     - ref_count: reference counter, to know when this shared matrix
!!       can be safely deallocated
!!     - id_nr: an identification number to ease equality tests between
!!       matrices and for debugging
!!     - sm: the underlying real matrix (make more efficient?
!!       direct implementation?)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_block_matrix_type
     LOGICAL :: initialized
     LOGICAL :: frozen
     LOGICAL :: should_dealloc_matrix
     LOGICAL :: b_rep_is_primary
     LOGICAL :: b_i_rep_in_sync
     LOGICAL :: cache_is_valid
     TYPE(cp_b_matrix_struct_type), POINTER :: matrix_struct
     INTEGER :: ref_count, id_nr
     TYPE(real_matrix_type), POINTER :: sm
  END TYPE cp_block_matrix_type
!!***
!****************************************************************************

!!****s* cp_block_matrix/cp_block_matrix_p_type [1.0] *
!!
!!   NAME
!!     cp_block_matrix_p_type
!!
!!   FUNCTION
!!     just to have arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - matrix: the pointer to the matrix
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_block_matrix_p_type
     TYPE(cp_block_matrix_type), POINTER :: matrix
  END TYPE cp_block_matrix_p_type
!!***
!****************************************************************************

!!****s* cp_block_matrix_type/cp_sp_block_handle_type [1.0] *
!!
!!   NAME
!!     cp_sp_block_handle_type
!!
!!   FUNCTION
!!     represent a block (to release the data if it is cached, and made
!!     get statement different from matrix type to matrix type)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     matrix: the matrix from where comes the block
!!     ... to do
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_sp_block_handle_type
     TYPE(cp_block_matrix_type), POINTER :: matrix
     INTEGER :: block_row, block_col
     CHARACTER :: trans,trans_to_do
     REAL(KIND = dp), DIMENSION(:,:), POINTER :: block_val
     LOGICAL :: should_dealloc_data, initialized
  END TYPE cp_sp_block_handle_type
!!***
!****************************************************************************

CONTAINS

!!****f* cp_block_matrix/cp_sp_create [1.0] *
!!
!!   NAME
!!     cp_sp_create
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_create(matrix, matrix_struct, u_matrix,&
!!         should_dealloc_matrix, b_rep_is_primary,&
!!         allocate_nonsparse_blocks, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Type(cp_b_matrix_struct_type), Intent (IN), Target:: matrix_struct
!!       Type(real_matrix_type), Optional, Intent (IN), Target:: u_matrix
!!       Logical, Intent (IN), Optional:: should_dealloc_matrix,&
!!         b_rep_is_primary, allocate_nonsparse_blocks
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_create
!!
!!   FUNCTION
!!     allocates a new sparse matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to initialize
!!     - matrix_struct: the structure of this matrix
!!     - u_matrix: the underlying sparse matrix SHARED
!!     - should_dealloc_matrix:if the underlying matrix should be deallocated
!!       (defaults to true)
!!     - b_rep_is_primary: if the block representation is primary
!!       (defaults to true)
!!     - allocate_nonsparse_blocks: if true allocates all the blocks outside
!!       the sparsity. Defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_create(matrix, matrix_struct, u_matrix, &
       should_dealloc_matrix, b_rep_is_primary, allocate_nonsparse_blocks, error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_b_matrix_struct_type), &
      INTENT(in), TARGET                     :: matrix_struct
    TYPE(real_matrix_type), INTENT(in), &
      OPTIONAL, TARGET                       :: u_matrix
    LOGICAL, INTENT(in), OPTIONAL            :: should_dealloc_matrix, &
                                                b_rep_is_primary, &
                                                allocate_nonsparse_blocks
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=10)                        :: m_symm
    INTEGER                                  :: block_cols, block_rows, stat
    INTEGER, SAVE                            :: last_id_nr = 0
    LOGICAL                                  :: failure
    TYPE(cp_b_dims_type), POINTER            :: b_dims

    failure=.FALSE.

    ALLOCATE(matrix,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_id_nr=last_id_nr+1
       matrix%id_nr=last_id_nr
       matrix%frozen=.FALSE.
       matrix%cache_is_valid=.FALSE.
       matrix%should_dealloc_matrix=.TRUE.
       IF (PRESENT(should_dealloc_matrix)) &
            matrix%should_dealloc_matrix=should_dealloc_matrix
       matrix%b_rep_is_primary = .TRUE.
       IF (PRESENT(b_rep_is_primary)) matrix%b_rep_is_primary=b_rep_is_primary
       matrix%b_i_rep_in_sync = .FALSE.
       matrix%matrix_struct => matrix_struct
       CALL cp_retain(matrix%matrix_struct,error=error)
       matrix%ref_count=1
       NULLIFY(matrix%sm)
       matrix%initialized=.TRUE.
       IF (PRESENT(u_matrix)) THEN
          matrix%sm => u_matrix
       ELSE
          IF (matrix%matrix_struct%symmetric) THEN
             m_symm='symmetric'
          ELSE
             m_symm='NONE'
          END IF
          b_dims => matrix%matrix_struct%b_dims
          block_rows=SIZE(b_dims%block_begins_at_row)-1
          block_cols=SIZE(b_dims%block_begins_at_col)-1
          CALL allocate_matrix(matrix%sm,&
               nrow=b_dims%block_begins_at_row(block_rows+1)-1,&
               ncol=b_dims%block_begins_at_col(block_cols+1)-1,&
               nblock_row=block_rows,&
               nblock_col=block_cols,&
               first_row=b_dims%block_begins_at_row(1:block_rows),&
               last_row=b_dims%block_begins_at_row(2:(block_rows+1))-1,&
               first_col=b_dims%block_begins_at_col(1:block_cols),&
               last_col=b_dims%block_begins_at_col(2:(block_cols+1))-1,&
               matrix_name='matrix'//cp_to_string(matrix%id_nr),&
               matrix_symmetry=m_symm)
       END IF
       IF (PRESENT(allocate_nonsparse_blocks)) THEN
          IF (allocate_nonsparse_blocks) THEN
             CALL cp_sp_alloc_nonsparse_blocks(matrix,error=error)
          END IF
       END IF
    END IF
  END SUBROUTINE cp_sp_create
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_dealloc [1.0] *
!!
!!   NAME
!!     cp_sp_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_dealloc(matrix, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_dealloc
!!
!!   FUNCTION
!!     deallocates a sparse matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to be deallocated
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_dealloc(matrix, error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(matrix)) THEN
       CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
       CPPreconditionNoFail(matrix%ref_count==0,cp_warning_level,routineP,error)
       IF (.NOT. failure) THEN
          matrix%initialized=.FALSE.
          CALL cp_release(matrix%matrix_struct, error=error)
          NULLIFY(matrix%matrix_struct)
          IF (matrix%should_dealloc_matrix) THEN
             IF (ASSOCIATED(matrix%sm)) THEN
                CALL deallocate_matrix(matrix%sm)
             END IF
          ELSE
             NULLIFY(matrix%sm)
          END IF
          DEALLOCATE(matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE cp_sp_dealloc
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_retain [1.0] *
!!
!!   NAME
!!     cp_sp_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_retain(matrix, error)
!!       Type(cp_block_matrix_type), Intent (INOUT), Target:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_retain
!!
!!   FUNCTION
!!     augments the retain count by one (to be called when you need  and 
!!     want to keep arround this shared object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix you want to retain
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_retain(matrix, error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_block_matrix_type), POINTER      :: m_ptr

    failure=.FALSE.
    m_ptr => matrix

    CPPrecondition(ASSOCIATED(m_ptr),cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CPPreconditionNoFail(matrix%ref_count>0,cp_warning_level,routineP,error)
       matrix%ref_count=matrix%ref_count+1
    END IF
  END SUBROUTINE cp_sp_retain
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_release [1.0] *
!!
!!   NAME
!!     cp_sp_release
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_release(matrix, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_release
!!
!!   FUNCTION
!!     decreases the retain count by one (to be called when you don't need
!!     a retained object anymore)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix you want to release
!!     - error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_release(matrix, error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(matrix)) THEN
       CPPreconditionNoFail(matrix%ref_count>0,cp_warning_level,routineP,error)
       matrix%ref_count=matrix%ref_count-1
       IF (matrix%ref_count==0) THEN
          CALL cp_sp_dealloc(matrix,error=error)
       END IF
    END IF
    NULLIFY(matrix)
  END SUBROUTINE cp_sp_release
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_get_matrix [1.0] *
!!
!!   NAME
!!     cp_sp_get_matrix
!!
!!   SYNOPSIS
!!     Function cp_sp_get_matrix(matrix, error) Result(res)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Type(real_matrix_type), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sp_get_matrix
!!
!!   FUNCTION
!!     returns the underlying sparse matrix (might go away)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix you want info from
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  FUNCTION cp_sp_get_matrix(matrix, error) RESULT(res)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(real_matrix_type), POINTER          :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_get_matrix', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    
    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    res => matrix%sm
  END FUNCTION cp_sp_get_matrix
!***************************************************************************

!================ block matrix op =============
! operations related with the fact that we have a block matrix and we want
! to work with the blocks

!!****f* cp_block_matrix/cp_sp_set_b_rep [1.0] *
!!
!!   NAME
!!     cp_sp_set_b_rep
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_consolidate_matrix(matrix, error)
!!       Type(cp_block_matrix_type), Intent (INOUT):: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_consolidate_matrix
!!
!!   FUNCTION
!!     Switches to the block representation. Must be called before any get
!!     (even the get of owned blocks).
!!
!!   NOTES
!!     Does nothing, is not really nedeed for this kind of distribution,
!!     but must be called to make the transition to other matrix types easy. 
!!
!!   INPUTS
!!     matrix: the matrix to consolidate
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_set_b_rep(matrix,error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_set_b_rep', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       IF (.not.matrix%b_rep_is_primary) THEN
          IF (cp_debug .AND. debug_this_module) THEN
             CALL mp_sync(matrix%matrix_struct%distribution%para_env%group)
          END IF
          matrix%b_rep_is_primary=.TRUE.
       END IF
    END IF
  END SUBROUTINE cp_sp_set_b_rep
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_set_i_rep [1.0] *
!!
!!   NAME
!!     cp_sp_set_i_rep
!!
!!   FUNCTION
!!     Switches to the internal representation. Must be called before
!!     some operations like matrix multiplication.
!!
!!   NOTES
!!     Does nothing, is not really nedeed for this kind of distribution,
!!     but must be called to make the transition to other matrix types easy. 
!!
!!   INPUTS
!!     matrix: the matrix to consolidate
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_set_i_rep(matrix,error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_set_b_rep', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       IF (matrix%b_rep_is_primary) THEN
          IF (cp_debug .AND. debug_this_module) THEN
             CALL mp_sync(matrix%matrix_struct%distribution%para_env%group)
          END IF
          matrix%b_rep_is_primary=.FALSE.
       END IF
    END IF
  END SUBROUTINE cp_sp_set_i_rep
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_flush_cache [1.0] *
!!
!!   NAME
!!     cp_sp_flush_cache
!!
!!   FUNCTION
!!     removes the cached data (all non owned blocks)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix that should be purged
!!     - flush_secondary_rep: if true releases the secondary representation
!!       (to really free all the possible memory). Defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_flush_cache(matrix,flush_secondary_rep,error)
    TYPE(cp_block_matrix_type), INTENT(in)   :: matrix
    LOGICAL, INTENT(in), OPTIONAL            :: flush_secondary_rep
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bcol, brow, mepos
    LOGICAL                                  :: failure
    TYPE(cp_b_distribution_type), POINTER    :: distribution
    TYPE(cp_matrix_block_iterator)           :: b_iterator

    NULLIFY(distribution)
    failure=.FALSE.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_sm_b_i_init(b_iterator,matrix=matrix%sm,error=error)
       distribution => matrix%matrix_struct%distribution
       mepos=matrix%matrix_struct%distribution%para_env%mepos
       CPPrecondition(ASSOCIATED(distribution),cp_failure_level,routineP,error,failure)
       DO
          IF (.not.cp_sm_b_i_next(b_iterator,block_row=brow,&
               block_col=bcol)) EXIT
          IF (distribution%block_owner(brow,bcol) /= mepos) THEN
             CALL remove_block_node(matrix=matrix%sm, block_row=brow,&
                  block_col=bcol,error=error)
          END IF
       END DO
       CALL cp_sm_b_i_dealloc_ref(b_iterator,error=error)
    END IF
  END SUBROUTINE cp_sp_flush_cache
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_block_handle_init [1.0] *
!!
!!   NAME
!!     cp_sp_block_handle_init
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_block_handle_init(block_handle, error)
!!       Type(cp_sp_block_handle_type), Intent (OUT):: block_handle
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_block_handle_init
!!
!!   FUNCTION
!!     initializes a block handle, to be called before getting a block in it
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_block_handle_init(block_handle,error)
    TYPE(cp_sp_block_handle_type), &
      INTENT(out)                            :: block_handle
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_block_handle_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  block_handle%block_col=-1
  block_handle%block_row=-1
  block_handle%trans='N'
  block_handle%trans_to_do='N'
  block_handle%should_dealloc_data=.FALSE.
  NULLIFY(block_handle%block_val, block_handle%matrix)
  block_handle%initialized=.TRUE.
END SUBROUTINE cp_sp_block_handle_init
!***************************************************************************

!!****f* cp_block_matrix_type/cp_sp_get_local_block [1.0] *
!!
!!   NAME
!!     cp_sp_get_local_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_get_local_block(matrix, block_row, block_col,&
!!         block_handle, trans, trans_to_do, block_val, local_owner, owner,&
!!         alloc_block, error)
!!       Type(cp_block_matrix_type), Intent (INOUT), Target:: matrix
!!       Integer, Intent (IN):: block_row, block_col
!!       Type(cp_sp_block_handle_type), Intent (INOUT):: block_handle
!!       Character(Len=1), Intent (IN), Optional:: trans
!!       Character(Len=1), Intent (OUT), Optional:: trans_to_do
!!       Real(KIND = dp), Dimension(:,:), Optional, Pointer:: block_val
!!       Logical, Intent (OUT), Optional:: local_owner
!!       Integer, Intent (OUT), Optional:: owner
!!       Logical, Intent (IN), Optional:: alloc_block
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_get_local_block
!!
!!   FUNCTION
!!     returns a block of the matrix, that is stored locally
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: matrix to get the info from
!!     - block_row: the block row to get
!!     - block_col: the block col to get
!!     - trans: the transposition wanted (N: None, T: Transpose,
!!       C: Conjugate transpose, B: Bar (just conjugate, no transpose))
!!     - block_handle: a handle that is initialized with the info 
!!       regarding the block
!!     - block_val: the values of the block
!!     - trans_to_do: the transition still to perform (defined like trans)
!!       if block_val is asked but not trans_to_do then all the
!!       transpositions are always performed
!!     - alloc_block: if true allocates a block that is not stored locally,
!!       otherwise Defaults to false.
!!     - local_owner: true if the actual processor is the owner of the
!!       block
!!     - owner: mpi id of the owner of the block
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_get_local_block(matrix,block_row,block_col,&
       block_handle,trans,trans_to_do,block_val,&
       local_owner,owner,alloc_block,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    INTEGER, INTENT(in)                      :: block_row, block_col
    TYPE(cp_sp_block_handle_type), &
      INTENT(inout)                          :: block_handle
    CHARACTER(len=1), INTENT(in), OPTIONAL   :: trans
    CHARACTER(len=1), INTENT(out), OPTIONAL  :: trans_to_do
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: block_val
    LOGICAL, INTENT(out), OPTIONAL           :: local_owner
    INTEGER, INTENT(out), OPTIONAL           :: owner
    LOGICAL, INTENT(in), OPTIONAL            :: alloc_block
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_get_local_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tmp_block

    failure=.FALSE.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(block_handle%initialized,cp_warning_level,routineP,error)
    CPPreconditionNoFail(.not.ASSOCIATED(block_handle%block_val),cp_warning_level,routineP,error)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(.not.matrix%frozen,cp_warning_level,routineP,error)
       block_handle%trans='N'
       IF (PRESENT(trans)) block_handle%trans=trans
       CALL uppercase(block_handle%trans)
       block_handle%trans_to_do=block_handle%trans
       CALL cp_assert(block_handle%trans=='N'.or.block_handle%trans=='T',&
            cp_warning_level,cp_precondition_failed,routineP,&
            "unknown requested transposition '"//block_handle%trans//"' in "//&
            CPSourceFileRef,&
            error=error)
       block_handle%block_row=block_row
       block_handle%block_col=block_col
       IF (matrix%matrix_struct%symmetric .AND. block_row > block_col) THEN
          IF (block_handle%trans=='N') THEN
             block_handle%trans_to_do='T'
          ELSE
             block_handle%trans_to_do='N'
          END IF
          CALL get_block_node(matrix%sm, block_row=block_handle%block_col, &
               block_col=block_handle%block_row, BLOCK=tmp_block)
          IF (PRESENT(alloc_block)) THEN
             IF (alloc_block .AND. .not.ASSOCIATED(tmp_block)) THEN
                CALL add_block_node(matrix%sm,&
                     block_row=block_handle%block_col,&
                     block_col=block_handle%block_row, BLOCK=tmp_block)
             END IF
          END IF
       ELSE
          CALL get_block_node(matrix%sm, block_row=block_handle%block_row, &
               block_col=block_handle%block_col, BLOCK=tmp_block)
          IF (PRESENT(alloc_block)) THEN
             IF (alloc_block .AND. .not.ASSOCIATED(tmp_block)) THEN
                CALL add_block_node(matrix%sm, block_row=block_handle%block_row,&
                     block_col=block_handle%block_col, BLOCK=tmp_block)
             END IF
          END IF
       END IF
       IF (PRESENT(trans_to_do).OR.(.not.PRESENT(block_val)).OR.&
            block_handle%trans_to_do=='N') THEN
          block_handle%block_val => tmp_block
          block_handle%should_dealloc_data=.FALSE.
       ELSE
          block_handle%trans_to_do='N'
          IF (ASSOCIATED(tmp_block)) THEN
             ALLOCATE(block_handle%block_val(SIZE(tmp_block,2),&
                  SIZE(tmp_block,1)), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             IF (.not.failure) THEN
                DO i=1,SIZE(tmp_block,2)
                   block_handle%block_val(i,:)=tmp_block(:,i)
                END DO
             END IF
          ELSE
             NULLIFY(block_handle%block_val)
          END IF
          block_handle%should_dealloc_data=.TRUE.
       END IF
       IF (PRESENT(trans_to_do)) trans_to_do = block_handle%trans_to_do
       IF (PRESENT(block_val)) block_val => block_handle%block_val
       IF (PRESENT(local_owner)) THEN
          local_owner=matrix%matrix_struct%distribution%block_owner&
               (block_handle%block_row,block_handle%block_col)==&
               matrix%matrix_struct%distribution%para_env%mepos
       END IF
       IF (PRESENT(owner)) THEN
          owner=matrix%matrix_struct%distribution%block_owner&
               (block_handle%block_row,block_handle%block_col)
       END IF
       block_handle%initialized=.TRUE.
    END IF
  END SUBROUTINE cp_sp_get_local_block
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_block_handle_release [1.0] *
!!
!!   NAME
!!     cp_sp_block_handle_release
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_block_handle_release(block_handle, changed_values,&
!!         error)
!!       Type(cp_sp_block_handle_type), Intent (INOUT):: block_handle
!!       Logical, Intent (IN):: changed_values
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_block_handle_release
!!
!!   FUNCTION
!!     Releases the memory that might have been allocated to respond
!!     to the get. MUST be called after a block that was obtained with get
!!     is no longer used. If some values were changed in the array then
!!     changed_values should be true.
!!
!!   NOTES
!!     slightly nonstandard use of the release name (rename routine?)
!!
!!   INPUTS
!!     - block_handle: the block handle to be released
!!     - changed_values: if the values of the block have been changed
!!       defaults to false.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_block_handle_release(block_handle,changed_values,error)
    TYPE(cp_sp_block_handle_type), &
      INTENT(inout)                          :: block_handle
    LOGICAL, INTENT(in)                      :: changed_values
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_block_handle_release', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: tt
    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_handle%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (changed_values.and.block_handle%matrix%frozen) THEN
          CALL cp_assert(.FALSE.,&
               cp_warning_level,cp_assertion_failed,routineP,&
               "changed values in frozen matrix",error)
          block_handle%matrix%frozen=.FALSE.
       END IF
       IF (changed_values .AND. block_handle%should_dealloc_data) THEN
          IF (block_handle%trans_to_do==block_handle%trans) THEN
             tt='N'
          ELSE
             tt='T'
          END IF
          CALL cp_sp_set_local_block(block_handle%matrix, &
               block_row=block_handle%block_row,&
               block_col=block_handle%block_col,&
               trans=tt,local_scratch=.TRUE.,&
               block_val=block_handle%block_val,&
               error=error)
       END IF
       IF (block_handle%should_dealloc_data) THEN
          IF (ASSOCIATED(block_handle%block_val)) THEN
             DEALLOCATE(block_handle%block_val,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END IF
       NULLIFY(block_handle%block_val)
    END IF
  END SUBROUTINE cp_sp_block_handle_release
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_set_local_block [1.0] *
!!
!!   NAME
!!     cp_sp_set_local_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_set_local_block(matrix, block_row, block_col,&
!!         block_val, trans, local_scratch, error)
!!       Type(cp_block_matrix_type), Intent (INOUT), Target:: matrix
!!       Integer, Intent (IN):: block_row, block_col
!!       Real(KIND = dp), Intent (IN), Dimension(:,:), Target:: block_val
!!       Character, Intent (IN), Optional:: trans
!!       Logical, Intent (IN), Optional:: local_scratch
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_set_local_block
!!
!!   FUNCTION
!!     sets the value in a local block
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix that should be changed
!!     block_row: the row that should be changed
!!     block_col: the column that should be updated
!!     block_val: the new values of the block
!!     trans: if the data is transposed (N: non transposed, T: transposed)
!!     local_scratch: if the actual block might be not owned, and is just
!!               a scratch block to be stored in the cache
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_set_local_block(matrix,block_row,block_col,block_val,&
       trans,local_scratch,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    INTEGER, INTENT(in)                      :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in), TARGET                     :: block_val
    CHARACTER, INTENT(in), OPTIONAL          :: trans
    LOGICAL, INTENT(in), OPTIONAL            :: local_scratch
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_set_local_block', &
      routineP = moduleN//':'//routineN

    CHARACTER                                :: t_to_do
    INTEGER                                  :: b_col, b_row, i, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tmp_block

    failure=.FALSE.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       CALL cp_assert(.not.matrix%frozen,&
            cp_warning_level,cp_assertion_failed,routineP,&
            "changed values in frozen matrix",error)
       matrix%frozen=.FALSE.
    END IF
    IF (.NOT. failure) THEN
       t_to_do='N'
       IF (PRESENT(trans)) t_to_do=trans
       CALL uppercase(t_to_do)
       b_row=block_col
       b_row=block_row
       IF (matrix%matrix_struct%symmetric.and.block_row > block_col) THEN
          IF (t_to_do=='N') THEN
             t_to_do='T'
          ELSE
             t_to_do='N'
          END IF
          b_row=block_row
          b_row=block_col
       END IF
       IF (t_to_do=='T') THEN
          ALLOCATE(tmp_block(SIZE(block_val,2),SIZE(block_val,1)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,SIZE(block_val,2)
             tmp_block(i,:)=block_val(:,i)
          END DO
          CALL put_block_node(matrix%sm, block_row=b_row,&
               block_col=b_col, BLOCK=tmp_block)
          DEALLOCATE(tmp_block,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ELSE
          tmp_block => block_val
          CALL put_block_node(matrix%sm, block_row=b_row,&
               block_col=b_col, BLOCK=tmp_block)     
       END IF
    ENDIF
  END SUBROUTINE cp_sp_set_local_block
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_remove_block [1.0] *
!!
!!   NAME
!!     cp_sp_remove_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_remove_block(matrix, block_row, block_col, error)
!!       Type(cp_block_matrix_type), Intent (INOUT), Target:: matrix
!!       Integer, Intent (IN):: block_col, block_row
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_remove_block
!!
!!   FUNCTION
!!     Removes the local block at the given position.
!!     Works both for owned and cached blocks.
!!     Does not chenge the distribution structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_remove_block(matrix,block_row, block_col,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    INTEGER, INTENT(in)                      :: block_row, block_col
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_remove_block', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       IF (matrix%matrix_struct%symmetric.and.block_row > block_col) THEN
          CALL remove_block_node(matrix%sm,block_row=block_col,&
               block_col=block_row, error=error)
       ELSE
          CALL remove_block_node(matrix%sm,block_row=block_row,&
               block_col=block_col, error=error)
       END IF
    END IF
  END SUBROUTINE cp_sp_remove_block
!***************************************************************************
!!****f* cp_block_matrix/cp_sp_array_create [1.0] *
!!
!!   NAME
!!     cp_sp_array_create
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_array_create(matrixes, n_matrixes, matrix_struct,&
!!         b_rep_is_primary, allocate_nonsparse_blocks, error)
!!       Type(cp_block_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Type(cp_b_matrix_struct_type), Intent (IN), Target:: matrix_struct
!!       Integer, Intent (IN):: n_matrixes
!!       Logical, Intent (IN), Optional:: b_rep_is_primary,&
!!         allocate_nonsparse_blocks
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_array_create
!!
!!   FUNCTION
!!     allocates a one dimensional array of sparse matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: the array of matrixes to allocate and initialize
!!     - n_matrixes: the number of matrixes to allocate in the vector
!!     - matrix_struct: the structure of this matrix
!!     - b_rep_is_primary: if the block representation is primary
!!       (defaults to true)
!!     - allocate_nonsparse_blocks: if the blocks outside the sparsity 
!!!      should be allocated (defaults to false)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_array_create(matrixes, n_matrixes, matrix_struct, &
      b_rep_is_primary, allocate_nonsparse_blocks, error)
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), POINTER                  :: matrixes
    INTEGER, INTENT(in)                      :: n_matrixes
    TYPE(cp_b_matrix_struct_type), &
      INTENT(in), TARGET                     :: matrix_struct
    LOGICAL, INTENT(in), OPTIONAL            :: b_rep_is_primary, &
                                                allocate_nonsparse_blocks
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_array_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    
    ALLOCATE(matrixes(n_matrixes),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       DO i=1,n_matrixes
          NULLIFY(matrixes(i)%matrix)
          CALL cp_sp_create(matrixes(i)%matrix,&
              matrix_struct=matrix_struct,&
              b_rep_is_primary=b_rep_is_primary,&
              allocate_nonsparse_blocks=allocate_nonsparse_blocks,&
              error=error)
       END DO
    END IF
  END SUBROUTINE cp_sp_array_create

!!****f* cp_block_matrix/cp_sp_dealloc_m_array [1.0] *
!!
!!   NAME
!!     cp_sp_dealloc_m_array
!!
!!   FUNCTION
!!     deallocates an array of matrices (releasing each matrix)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - m_array: the array of matrices
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_array_dealloc(m_array, error)
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), POINTER                  :: m_array
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_array_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(m_array)) THEN
       DO i=1,SIZE(m_array)
          CALL cp_release(m_array(i)%matrix,error=error)
       END DO
       DEALLOCATE(m_array,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
    NULLIFY(m_array)
  END SUBROUTINE cp_sp_array_dealloc
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_array_release [1.0] *
!!
!!   NAME
!!     cp_sp_array_release
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_array_release(array, error)
!!       Type(cp_block_matrix_p_type), Dimension(:), Intent (INOUT):: array
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_array_release
!!
!!   FUNCTION
!!     releases all the elements of an array of cp_block_matrix_p_type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - array: the array to be released
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_array_release(array,error)
    TYPE(cp_block_matrix_p_type), &
      DIMENSION(:), INTENT(inout)            :: array
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_array_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  DO i=1,SIZE(array)
     CALL cp_release(array(i)%matrix,error=error)
     NULLIFY(array(i)%matrix)
  END DO
END SUBROUTINE cp_sp_array_release
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_alloc_nonsparse_blocks *
!!
!!   NAME
!!     cp_sp_alloc_nonsparse_blocks
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_alloc_nonsparse_blocks(matrix, error)
!!       Type(cp_block_matrix_type), Intent (INOUT), Target:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_alloc_nonsparse_blocks
!!
!!   FUNCTION
!!     guarantees that all the blocks that are declared nonzero in the
!!     sparsity of the matrix are allocated
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to be modified
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     5.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_alloc_nonsparse_blocks(matrix,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_alloc_nonsparse_blocks', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL cp_sm_alloc_nonsparse_blocks(matrix%sm,&
          m_struct=matrix%matrix_struct,&
          para_env=matrix%matrix_struct%distribution%para_env,&
          error=error)
  END IF
END SUBROUTINE cp_sp_alloc_nonsparse_blocks
!***************************************************************************

!============= sparsity ===========

!!****f* cp_block_matrix/cp_sp_update_sparsity [1.0] *
!!
!!   NAME
!!     cp_sp_update_sparsity
!!
!!   FUNCTION
!!     updates the sparsity of the matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to update
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_update_sparsity(matrix,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_update_sparsity', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="update sparsity", error=error)
     ! to do
  END IF
END SUBROUTINE cp_sp_update_sparsity
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_check_sparsity [1.0] *
!!
!!   NAME
!!     cp_sp_check_sparsity
!!
!!   FUNCTION
!!    checks if the sparsity of the matrix is respected
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to check
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_check_sparsity(matrix,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_check_sparsity', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="check sparsity", error=error)
     ! to do
  END IF
END SUBROUTINE cp_sp_check_sparsity
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_enforce_sparsity [1.0] *
!!
!!   NAME
!!     cp_sp_enforce_sparsity
!!
!!   FUNCTION
!!    enforces the sparsity of the matrix (removes all the blocks 
!!    outside the sparsity)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to check
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_enforce_sparsity(matrix,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_enforce_sparsity', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="enforce sparsity", error=error)
     ! to do
  END IF
END SUBROUTINE cp_sp_enforce_sparsity
!***************************************************************************

!============= distribution ===========

!!****f* cp_block_matrix/cp_sp_check_distribution [1.0] *
!!
!!   NAME
!!     cp_sp_check_distribution
!!
!!   FUNCTION
!!    checks if the distribution of the matrix is respected
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to check
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_check_distribution(matrix,error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_check_distribution', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="check distribution", error=error)
     ! to do
  END IF
END SUBROUTINE cp_sp_check_distribution
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_get [1.0] *
!!
!!   NAME
!!     cp_sp_get
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_get(matrix, matrix_struct, distribution, dimensions,&
!!         sparsity, n_block_rows, n_block_cols, n_rows, n_cols, frozen,&
!!         b_rep_is_primary, id_nr, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Type(cp_b_matrix_struct_type), Pointer, Optional:: matrix_struct
!!       Type(cp_b_distribution_type), Pointer, Optional:: distribution
!!       Type(cp_b_dims_type), Pointer, Optional:: dimensions
!!       Logical, Dimension(:,:), Pointer, Optional:: sparsity
!!       Integer, Intent (OUT), Optional:: n_block_cols, n_block_rows,&
!!         n_rows, n_cols, id_nr
!!       Logical, Intent (OUT), Optional:: frozen, b_rep_is_primary
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_get
!!
!!   FUNCTION
!!     returns various attributes of the sparese matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the sparse matrix you want info about
!!     - matrix_struct: the structure of this matrix
!!     - distribution: the distribution of this matrix
!!     - dimensions: the size of this matrix
!!     - sparsity: the sparsity of this matrix
!!     - n_block_rows: the number of blocks in a row
!!     - n_block_cols: the number of blocks in a column
!!     - n_rows: the total number of rows
!!     - n_cols: the total number of columns
!!     - frozen: if the matrix is forzen
!!     - b_rep_is_primary: if the block representation is the primary
!!       representation
!!     - id_nr: identification number (to compare matrices)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_get(matrix,matrix_struct, distribution,dimensions,&
     sparsity, n_block_rows, n_block_cols, n_rows, n_cols, frozen,&
     b_rep_is_primary, id_nr, error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    TYPE(cp_b_matrix_struct_type), &
      OPTIONAL, POINTER                      :: matrix_struct
    TYPE(cp_b_distribution_type), OPTIONAL, &
      POINTER                                :: distribution
    TYPE(cp_b_dims_type), OPTIONAL, POINTER  :: dimensions
    LOGICAL, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: sparsity
    INTEGER, INTENT(out), OPTIONAL           :: n_block_rows, n_block_cols, &
                                                n_rows, n_cols
    LOGICAL, INTENT(out), OPTIONAL           :: frozen, b_rep_is_primary
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  IF (.not.failure) THEN
     CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT. failure) THEN
     IF (PRESENT(matrix_struct)) matrix_struct => matrix%matrix_struct
     IF (PRESENT(distribution)) distribution => &
          matrix%matrix_struct%distribution
     IF (PRESENT(dimensions)) dimensions => matrix%matrix_struct%b_dims
     IF (PRESENT(sparsity)) sparsity => matrix%matrix_struct%sparsity
     IF (PRESENT(n_block_rows)) n_block_rows = &
          SIZE(matrix%matrix_struct%b_dims%block_begins_at_row)-1
     IF (PRESENT(n_block_cols)) n_block_cols = &
          SIZE(matrix%matrix_struct%b_dims%block_begins_at_col)-1
     IF (PRESENT(n_rows)) n_rows = &
          matrix%matrix_struct%b_dims%block_begins_at_row &
          (SIZE(matrix%matrix_struct%b_dims%block_begins_at_row))-1
     IF (PRESENT(n_cols)) n_cols = &
          matrix%matrix_struct%b_dims%block_begins_at_col &
          (SIZE(matrix%matrix_struct%b_dims%block_begins_at_col))-1
     IF (PRESENT(frozen)) frozen=matrix%frozen
     IF (PRESENT(b_rep_is_primary)) b_rep_is_primary=matrix%b_rep_is_primary
     IF (PRESENT(id_nr)) id_nr=matrix%id_nr
  END IF
END SUBROUTINE cp_sp_get
!***************************************************************************


!!****f* cp_block_matrix/cp_sp_redistribute [1.0] *
!!
!!   NAME
!!     cp_sp_redistribute
!!
!!   FUNCTION
!!    redistributes the blocks of the matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to redistibute
!!     - new_distribution: the new distribution
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!! 
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_redistribute(matrix, new_distribution, error)
    TYPE(cp_block_matrix_type), &
      INTENT(inout), TARGET                  :: matrix
    TYPE(cp_b_distribution_type), POINTER    :: new_distribution
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_redistribute', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (.NOT. failure) THEN
     CALL cp_unimplemented_error(fromWhere=routineP, &
          message="redistribute", error=error)
     ! to do
  END IF
END SUBROUTINE cp_sp_redistribute
!***************************************************************************\

!!****f* cp_block_matrix/cp_sp_set_to [1.0] *
!!
!!   NAME
!!     cp_sp_set_to
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_set_to(matrix, value, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Real(KIND = dp), Intent (IN), Optional:: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_set_to
!!
!!   FUNCTION
!!     sets the matrix to some value (normally 0)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrix: the matrix to initialize
!!     - value: the new value (defaults to 0)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_sp_set_to(matrix, value, error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sp_set_to', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL cp_sm_set_to(matrix%sm, value=value, error=error)
  END IF
END SUBROUTINE cp_sp_set_to
!***************************************************************************

!!****f* cp_block_matrix/cp_sp_multiplies_blacs [1.0] *
!!
!!   NAME
!!     cp_sp_multiplies_blacs
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_multiplies_blacs(matrix_a, matrix_b, matrix_c,&
!!         alpha, beta, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix_a
!!       Type(cp_fm_type), Pointer:: matrix_b, matrix_c
!!       Real(KIND = dp), Intent (IN), Optional:: alpha, beta
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine cp_sp_multiplies_blacs
!!
!!   FUNCTION
!!     multiplies a sparse matrix with a blacs matrix
!!     matrix_c = alpha op(matrix_a) op(matrix_b) + beta matrix_c
!!     or (if side='R')
!!     matrix_c = alpha op(matrix_b) op(matrix_a) + beta matrix_c
!!
!!   NOTES
!!     at the moment the doc is a "wishlist", I just call
!!     Joost's sparse_time_blacs
!!
!!   INPUTS
!!     - matrix_a: the sparse matrix to multiply
!!     - matrix_b: the blacs matrix to multiply
!!     - matrix_c: the blacs matrix to hold the result
!!     - side: on which side of matrix_b is applied matrix_a ('L' (the default):
!!       matrix_a * matrix_b, 'R' : matrix_b * matrix_a
!!     - alpha: the scaling of the product (defaults to 1.0)
!!     - beta: the scaling of the values in matrix_c (defaults to 0.0)
!!     - m: the number of rows of the resulting matrix
!!     - n: the number of columns of the resulting matrix (defaults to all
!!       the columns of matrix_b)
!!     - k: the of the columns of matrix_a (rows of matrix_b) that gets multiplied
!!     - tr_a: if matrix a should be transposed: if false (the default) 
!!       op(matrix_a)=matrix_a, else op(matrix_a) is the transpose of matrix_a.
!!     - tr_b: if matrix b should be transposed: if tr_b is false (the default)
!!       op(matrix_b)=matrix_b, else op(matrix_b) is the transpose of matrix_b.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_sp_multiplies_blacs(matrix_a,matrix_b,matrix_c,alpha,beta,error)
    TYPE(cp_block_matrix_type), POINTER      :: matrix_a
    TYPE(cp_fm_type), POINTER                :: matrix_b, matrix_c
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha, beta
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: ncol

    CALL cp_fm_get_info(matrix_b,ncol_global=ncol)
    ! JVDV block owner commented away, this should be useless right now  
    CALL cp_sm_fm_multiply(matrix_a%sm, matrix_b, matrix_c,ncol,&
         alpha=alpha, beta=beta)
  END SUBROUTINE cp_sp_multiplies_blacs

END MODULE cp_block_matrix
