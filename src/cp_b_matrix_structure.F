!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_b_matrix_structure [1.0] *
!!
!!   NAME
!!     cp_b_matrix_structure
!!
!!   FUNCTION
!!     represent the structure of a block sparse matrix.
!!     These types should be enough generic so that different types
!!     of matrices should be able to use them.
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!****************************************************************************
MODULE cp_b_matrix_structure
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_b_matrix_structure'

  PUBLIC :: cp_b_distribution_type,cp_b_dims_type,&
       cp_b_matrix_struct_type, cp_bs_iterator_type

  PUBLIC :: cp_create, cp_retain, cp_release, cp_get, cp_next,&
       cp_init, cp_dealloc_ref
  PUBLIC :: cp_freeze_sparsity, cp_get_sparsity_frozen, cp_get_block_info

  PUBLIC :: cp_dstruct_create, cp_dstruct_retain,&
       cp_dstruct_release, cp_dstruct_get, &
       cp_b_dims_create, cp_b_dims_retain,&
       cp_b_dims_release,&
       cp_bmstruct_create, cp_bmstruct_freeze_sparsity,&
       cp_bmstruct_get_sparsity_frozen, cp_bmstruct_retain,&
       cp_bmstruct_release, cp_bmstruct_get_block_info,&
       cp_bs_iter_init, cp_bs_iter_dealloc_ref, cp_bs_iter_next,&
       cp_bmstruct_get

  INTERFACE cp_create
     MODULE PROCEDURE cp_dstruct_create
  END INTERFACE
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_dstruct_dealloc
  END INTERFACE
  INTERFACE cp_retain
     MODULE PROCEDURE cp_dstruct_retain
  END INTERFACE
  INTERFACE cp_release
     MODULE PROCEDURE cp_dstruct_release
  END INTERFACE
  INTERFACE cp_get
     MODULE PROCEDURE cp_dstruct_get
  END INTERFACE

  INTERFACE cp_create
     MODULE PROCEDURE cp_b_dims_create
  END INTERFACE
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_b_dims_dealloc
  END INTERFACE
  INTERFACE cp_retain
     MODULE PROCEDURE cp_b_dims_retain
  END INTERFACE
  INTERFACE cp_release
     MODULE PROCEDURE cp_b_dims_release
  END INTERFACE

  INTERFACE cp_create
     MODULE PROCEDURE cp_bmstruct_create
  END INTERFACE
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_bmstruct_dealloc
  END INTERFACE
  INTERFACE cp_freeze_sparsity
     MODULE PROCEDURE cp_bmstruct_freeze_sparsity
  END INTERFACE
  INTERFACE cp_get
     MODULE PROCEDURE cp_bmstruct_get
  END INTERFACE
  INTERFACE cp_get_sparsity_frozen
     MODULE PROCEDURE cp_bmstruct_get_sparsity_frozen
  END INTERFACE
  INTERFACE cp_retain
     MODULE PROCEDURE cp_bmstruct_retain
  END INTERFACE
  INTERFACE cp_release
     MODULE PROCEDURE cp_bmstruct_release
  END INTERFACE
  INTERFACE cp_get_block_info
     MODULE PROCEDURE cp_bmstruct_get_block_info
  END INTERFACE

  INTERFACE cp_init
     MODULE PROCEDURE cp_bs_iter_init
  END INTERFACE
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_bs_iter_dealloc_ref
  END INTERFACE
  INTERFACE cp_next
     MODULE PROCEDURE cp_bs_iter_next
  END INTERFACE
!!***
!****************************************************************************

!!****s* cp_b_matrix_structure/cp_b_distribution_type [1.0] *
!!
!!   NAME
!!     cp_b_distribution_type
!!
!!   FUNCTION
!!     represent how the blocks of a sparse matrix are distribued among 
!!     the processors.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     initialized: true if the structure is initialized (bug catcher)
!!     block_owner: a matrix with the mpi id of the processor that has the
!!                  corresponding block.
!!     para_env: the mpi communicator, owner ids refer to this environment
!!     ref_count: an integer that keeps track of how many matrices are
!!                sharing this structure
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
  TYPE cp_b_distribution_type
     LOGICAL :: initialized
     TYPE(cp_para_env_type), POINTER :: para_env
     INTEGER, DIMENSION(:,:), POINTER :: block_owner
     INTEGER :: ref_count, id_nr
  END TYPE cp_b_distribution_type
!!***
  !****************************************************************************

!!****s* cp_b_matrix_structure/cp_b_dims_type [1.0] *
!!
!!   NAME
!!     cp_b_dims_type
!!
!!   FUNCTION
!!     Sores the size of the blocks of the matrix.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     initialized: true if the structure is initialized (bug catcher)
!!     block_begins_at_row: a vector of length n_block_rows+1 that 
!!                contains the index of the row at which the corresponding
!!                block begins (the last is one past the end of the last
!!                block)
!!     block_begins_at_col: a vector of length n_block_cols+1 that contains 
!!                the index of the column at which the corresponding block
!!                begins (the last is one past the end of the last block)  
!!     ref_count: an integer that keeps track of how many matrices are
!!                sharing this structure
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
  TYPE cp_b_dims_type
     LOGICAL :: initialized
     INTEGER, DIMENSION(:), POINTER :: block_begins_at_row, block_begins_at_col
     INTEGER :: ref_count, id_nr
  END TYPE cp_b_dims_type
!!***
  !****************************************************************************

!!****s* cp_b_matrix_structure/cp_b_matrix_struct_type [1.0] *
!!
!!   NAME
!!     cp_b_matrix_struct_type
!!
!!   FUNCTION
!!     Describes a matrix type.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     initialized: true if the structure is initialized (bug catcher)
!!     symmetric: true if the matrix is symmetric
!!     distribution: the distribution structure of this matrix
!!     b_dims: the dimensions of the matrix and its blocks
!!     sparsity: true if the corresponding block is non zero, if 
!!               not associated it means always true (no sparsity)
!!     sparsity_frozen: if the sparsity structure is frozen (no new non zero
!!            block can be added)
!!     ref_count: an integer that keeps track of how many matrices are
!!                sharing this structure
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
  TYPE cp_b_matrix_struct_type
     LOGICAL :: initialized, sparsity_frozen, symmetric
     TYPE(cp_b_dims_type), POINTER :: b_dims
     TYPE(cp_b_distribution_type), POINTER :: distribution
     LOGICAL, DIMENSION(:,:), POINTER :: sparsity
     INTEGER :: ref_count, id_nr
  END TYPE cp_b_matrix_struct_type
!!***
  !****************************************************************************

!!****s* cp_b_matrix_structure/cp_bs_iterator_type [1.0] *
!!
!!   NAME
!!     cp_bs_iterator_type
!!
!!   FUNCTION
!!     iterates on the distribution structure.
!!
!!   NOTES
!!     Putting different iterator types into this makes it inefficent 
!!     (many ifs), but this should not be the bottleneck
!!     I future it could iterates along the diagonal (to hopefully
!!     minimize dependencies between the processors)
!!
!!   ATTRIBUTES
!!     initialized: true if the structure has been initialized (bug catcher)
!!     do_global: if all the blocks should be iterated (otherwise only 
!!                the local)
!!     do_zero: if the blocks outside the sparsity structure should 
!!                be iterated (must be true if the sparsity structure is not
!!                frozen)
!!     block_row: the actual block row (-1 if invalid)
!!     block_col: the actual block column
!!     m_struct: the matrix structure you iterate on
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
  TYPE cp_bs_iterator_type
     LOGICAL :: initialized, do_global, do_zero, do_lower_half
     INTEGER :: block_row,block_col,n_block_rows, n_block_cols
     TYPE(cp_b_matrix_struct_type), POINTER :: m_struct
  END TYPE cp_bs_iterator_type
!!***
  !****************************************************************************

CONTAINS

  !======= block distribution structure =======

!!****f* cp_b_matrix_structure/cp_dstruct_create [1.0] *
!!
!!   NAME
!!     cp_dstruct_create
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_create(d_struct, para_env, n_block_rows,&
!!         n_block_cols, error)
!!       Type(cp_b_distribution_type), Pointer:: d_struct
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Integer, Intent (IN):: n_block_rows, n_block_cols
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_dstruct_create
!!
!!   FUNCTION
!!     allocates a distribution structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: the distribution structure
!!     para_env: the mpi environment for communication
!!     n_block_rows: the number of block rows
!!     n_block_cols: the number of block columns
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
  SUBROUTINE cp_dstruct_create(d_struct,para_env,n_block_rows,n_block_cols,&
       error)
    TYPE(cp_b_distribution_type), POINTER    :: d_struct
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(in)                      :: n_block_rows, n_block_cols
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_b_distribution_type', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    INTEGER, SAVE                            :: last_id_nr = 0
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(n_block_cols>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(n_block_rows>0,cp_failure_level,routineP,error,failure)
    IF (.not.failure) THEN
       ALLOCATE(d_struct,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       last_id_nr=last_id_nr+1
       d_struct%para_env => para_env
       CALL cp_para_env_retain(para_env,error=error)
       ALLOCATE(d_struct%block_owner(n_block_rows,n_block_cols),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       d_struct%block_owner=para_env%source
       d_struct%ref_count=1
       d_struct%id_nr=last_id_nr
       d_struct%initialized=.TRUE.
    END IF
  END SUBROUTINE cp_dstruct_create
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_dstruct_dealloc [1.0] *
!!
!!   NAME
!!     cp_dstruct_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_dealloc(d_struct, error)
!!       Type(cp_b_distribution_type), Pointer:: d_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_dstruct_dealloc
!!
!!   FUNCTION
!!     deallocates the memory allocated by the object AND the object itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: the distribution structure to deallocate
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
  SUBROUTINE cp_dstruct_dealloc(d_struct,error)
    TYPE(cp_b_distribution_type), POINTER    :: d_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dstruct_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(d_struct)) THEN
       CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          CPPreconditionNoFail(d_struct%ref_count==0,cp_warning_level,routineP,error)
          d_struct%ref_count=-1
          CALL cp_para_env_release(d_struct%para_env,error=error)
          IF (ASSOCIATED(d_struct%block_owner)) THEN
             DEALLOCATE(d_struct%block_owner,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          d_struct%initialized=.FALSE.
       END IF
       DEALLOCATE(d_struct,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF
  END SUBROUTINE cp_dstruct_dealloc
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_dstruct_retain [1.0] *
!!
!!   NAME
!!     cp_dstruct_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_retain(d_struct, error)
!!       Type(cp_b_distribution_type), Pointer:: d_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_dstruct_retain
!!
!!   FUNCTION
!!     augments the retain count by one (to be called if you want a shared 
!!     copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: the structure to retain
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
  SUBROUTINE cp_dstruct_retain(d_struct,error)
    TYPE(cp_b_distribution_type), POINTER    :: d_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dstruct_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (cp_debug) THEN
       CPPrecondition(ASSOCIATED(d_struct),cp_failure_level,routineP,error,failure)
       CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(d_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       d_struct%ref_count=d_struct%ref_count+1
    END IF
  END SUBROUTINE cp_dstruct_retain
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_dstruct_release [1.0] *
!!
!!   NAME
!!     cp_dstruct_release
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_release(d_struct, error)
!!       Type(cp_b_distribution_type), Pointer:: d_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_dstruct_release
!!
!!   FUNCTION
!!     decreases the ref_count by one (to be called if you do not want anymore
!!     your shared copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: the structure to release
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
  SUBROUTINE cp_dstruct_release(d_struct,error)
    TYPE(cp_b_distribution_type), POINTER    :: d_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dstruct_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(d_struct)) THEN
       IF (cp_debug) THEN
          CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
          CPPrecondition(d_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.not.failure) THEN
          d_struct%ref_count=d_struct%ref_count-1
          IF (d_struct%ref_count<1) THEN
             CALL cp_dstruct_dealloc(d_struct)
          END IF
       END IF
    END IF
    NULLIFY(d_struct)
  END SUBROUTINE cp_dstruct_release
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_dstruct_get [1.0] *
!!
!!   NAME
!!     cp_dstruct_get
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_get(d_struct, block_owner, para_env, error)
!!       Type(cp_b_distribution_type), Intent (IN):: d_struct
!!       Integer, Dimension(:,:), Pointer, Optional:: block_owner
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_dstruct_get
!!
!!   FUNCTION
!!     returns attributes about the given object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: the object you waht info about
!!     block_owner: the array with the block owners
!!     para_env: the mpi environement
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
  SUBROUTINE cp_dstruct_get(d_struct, block_owner, para_env, error)
    TYPE(cp_b_distribution_type), INTENT(in) :: d_struct
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: block_owner
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dstruct_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(block_owner)) block_owner => d_struct%block_owner
       IF (PRESENT(para_env)) para_env = d_struct%para_env
    END IF
  END SUBROUTINE cp_dstruct_get
  !***************************************************************************

  !====== dimensions structure ======
!!****f* cp_b_matrix_structure/cp_b_dims_create [1.0] *
!!
!!   NAME
!!     cp_b_dims_creates
!!
!!   FUNCTION
!!     allocates and initializes the dimensions of the blocks
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     b_dims: the structure that has to be initialized
!!     b_row_sizes: the row sizes of the blocks
!!     b_col_sizes: the column sizes of the blocks
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
  SUBROUTINE cp_b_dims_create(b_dims, b_row_sizes, b_col_sizes, error)
    TYPE(cp_b_dims_type), POINTER            :: b_dims
    INTEGER, DIMENSION(:), INTENT(in)        :: b_row_sizes, b_col_sizes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_b_dims_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    INTEGER, SAVE                            :: last_id_nr = 0
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(SIZE(b_row_sizes)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(b_col_sizes)>0,cp_failure_level,routineP,error,failure)
    IF(cp_debug) THEN
       CPPrecondition(ALL(b_row_sizes>=0),cp_failure_level,routineP,error,failure)
       CPPrecondition(ALL(b_col_sizes>=0),cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       ALLOCATE (b_dims,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       last_id_nr=last_id_nr+1
       b_dims%id_nr=last_id_nr
       b_dims%ref_count=1
       ALLOCATE (b_dims%block_begins_at_row(SIZE(b_row_sizes)+1),&
            b_dims%block_begins_at_col(SIZE(b_col_sizes)+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       b_dims%block_begins_at_row(1)=1
       DO i=1,SIZE(b_row_sizes)
          b_dims%block_begins_at_row(i+1)=b_dims%block_begins_at_row(i)+&
               b_row_sizes(i)
       END DO
       b_dims%block_begins_at_col(1)=1
       DO i=1,SIZE(b_col_sizes)
          b_dims%block_begins_at_col(i+1)=b_dims%block_begins_at_col(i)+&
               b_col_sizes(i)
       END DO
       b_dims%initialized=.TRUE.
    END IF
  END SUBROUTINE cp_b_dims_create
  !***************************************************************************

!!****f* cp_matrix_structure/cp_b_dims_dealloc [1.0] *
!!
!!   NAME
!!     cp_b_dims_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_b_dims_dealloc(b_dims, error)
!!       Type(cp_b_dims_type), Pointer:: b_dims
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_b_dims_dealloc
!!
!!   FUNCTION
!!     deallocates the memory allocated by the object and the object itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     b_dims: the structure 
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
  SUBROUTINE cp_b_dims_dealloc(b_dims,error)
    TYPE(cp_b_dims_type), POINTER            :: b_dims
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_b_dims_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(b_dims%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(b_dims%ref_count==0,cp_warning_level,routineP,error)
       IF (ASSOCIATED(b_dims)) THEN
          IF (ASSOCIATED(b_dims%block_begins_at_row)) THEN
             DEALLOCATE(b_dims%block_begins_at_row, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(b_dims%block_begins_at_col)) THEN
             DEALLOCATE(b_dims%block_begins_at_col,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
          b_dims%initialized=.FALSE.
          DEALLOCATE(b_dims,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE cp_b_dims_dealloc
  !***************************************************************************

!!****f* cp_matrix_structure/cp_b_dims_retain [1.0] *
!!
!!   NAME
!!     cp_b_dims_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_b_dims_retain(b_dims, error)
!!       Type(cp_b_dims_type), Intent (INOUT):: b_dims
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_b_dims_retain
!!
!!   FUNCTION
!!     augments the retain count by one (to be called if you want a shared 
!!     copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     b_dims: the object to retain
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
  SUBROUTINE cp_b_dims_retain(b_dims,error)
    TYPE(cp_b_dims_type), INTENT(inout)      :: b_dims
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_b_dims_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(b_dims%initialized,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(b_dims%ref_count>0,cp_warning_level,routineP,error)
    IF (.NOT. failure) THEN
       b_dims%ref_count=b_dims%ref_count+1
    END IF
  END SUBROUTINE cp_b_dims_retain
  !***************************************************************************

!!****f* cp_matrix_structure/cp_b_dims_release [1.0] *
!!
!!   NAME
!!     cp_b_dims_release
!!
!!   SYNOPSIS
!!     Subroutine cp_b_dims_release(b_dims, error)
!!       Type(cp_b_dims_type), Pointer:: b_dims
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_b_dims_release
!!
!!   FUNCTION
!!     decreases the ref_count by one (to be called if you do not want anymore
!!     your shared copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     b_dims: the object to release
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
  SUBROUTINE cp_b_dims_release(b_dims,error)
    TYPE(cp_b_dims_type), POINTER            :: b_dims
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_b_dims_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(b_dims)) THEN
       IF (cp_debug) THEN
          CPPrecondition(b_dims%initialized,cp_failure_level,routineP,error,failure)
          CPPrecondition(b_dims%ref_count>0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.not.failure) THEN
          b_dims%ref_count=b_dims%ref_count-1
          IF (b_dims%ref_count<1) THEN
             CALL cp_b_dims_dealloc(b_dims,error)
          END IF
       END IF
    END IF
    NULLIFY(b_dims)
  END SUBROUTINE cp_b_dims_release
  !***************************************************************************

  !======== matrix structure ========
!!****f* cp_b_matrix_structure/cp_bmstruct_create [1.0] *
!!
!!   NAME
!!     cp_bmstruct_create
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_create(m_struct, symmetric, d_struct, b_dims,&
!!         sparsity, no_sparsity, para_env, b_row_sizes, b_col_sizes,&
!!         error)
!!       Type(cp_b_matrix_struct_type), Pointer:: m_struct
!!       Type(cp_b_distribution_type), Target, Intent (INOUT), Optional::&
!!         d_struct
!!       Type(cp_b_dims_type), Target, Intent (INOUT), Optional:: b_dims
!!       Logical, Dimension(:,:), Optional:: sparsity
!!       Logical, Intent (IN), Optional:: no_sparsity, symmetric
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Integer, Dimension(:), Intent (IN), Optional:: b_row_sizes,&
!!         b_col_sizes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_create
!!
!!   FUNCTION
!!     allocated and initializes a block matrix structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     m_struct: the matrix structure
!!     d_struct: the distribution structure
!!     b_dims: the block dimensions
!!     sparsity: the spasity structure
!!     no_sparsity: if true there is no sparsity stucture (the same as
!!            passing a non associated sparsity pointer), defaults to false
!!     para_env: the mpi environement
!!     symmetric: if the matrix is symmetric (defaults to true)
!!     b_row_sizes: sizes of the rows of the blocks
!!     b_col_sizes: sizes of the cols of the blocks
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
  SUBROUTINE cp_bmstruct_create(m_struct,symmetric,d_struct,b_dims, &
       sparsity, no_sparsity, para_env, b_row_sizes, b_col_sizes,error)
    TYPE(cp_b_matrix_struct_type), POINTER   :: m_struct
    LOGICAL, INTENT(in), OPTIONAL            :: symmetric
    TYPE(cp_b_distribution_type), &
      INTENT(inout), OPTIONAL, TARGET        :: d_struct
    TYPE(cp_b_dims_type), INTENT(inout), &
      OPTIONAL, TARGET                       :: b_dims
    LOGICAL, DIMENSION(:, :), OPTIONAL       :: sparsity
    LOGICAL, INTENT(in), OPTIONAL            :: no_sparsity
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: b_row_sizes, b_col_sizes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    INTEGER, SAVE                            :: last_id_nr = 0
    LOGICAL                                  :: failure, no_spars

!, intent(in)

    failure=.FALSE.

    no_spars=.FALSE.
    IF (PRESENT(no_sparsity)) no_spars=no_sparsity
    ALLOCATE(m_struct, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       last_id_nr=last_id_nr+1
       m_struct%id_nr=last_id_nr
       m_struct%sparsity_frozen= PRESENT(sparsity).or.no_spars
       m_struct%symmetric=.TRUE. ! make it default to false
       IF (PRESENT(symmetric)) m_struct%symmetric=symmetric
       m_struct%ref_count=1

       IF (PRESENT(d_struct)) THEN
          m_struct%distribution => d_struct
          CALL cp_retain(m_struct%distribution, error=error)
       ELSE
          CPPrecondition(PRESENT(para_env),cp_failure_level,routineP,error,failure)
          IF (.not.PRESENT(b_dims)) THEN
             CPPrecondition(PRESENT(b_col_sizes),cp_failure_level,routineP,error,failure)
             CPPrecondition(PRESENT(b_row_sizes),cp_failure_level,routineP,error,failure)
             IF (.not.failure) THEN
                CALL cp_dstruct_create(m_struct%distribution,para_env=para_env,&
                     n_block_rows=SIZE(b_row_sizes),&
                     n_block_cols=SIZE(b_col_sizes),&
                     error=error)
             END IF
          ELSE
             CALL cp_dstruct_create(m_struct%distribution,para_env=para_env,&
                  n_block_rows=SIZE(b_dims%block_begins_at_row)-1,&
                  n_block_cols=SIZE(b_dims%block_begins_at_col)-1,&
                  error=error)
          END IF
       END IF

       IF (PRESENT(b_dims)) THEN
          m_struct%b_dims => b_dims
          CALL cp_retain(b_dims,error)
       ELSE
          CPPrecondition(PRESENT(b_col_sizes),cp_failure_level,routineP,error,failure)
          CPPrecondition(PRESENT(b_row_sizes),cp_failure_level,routineP,error,failure)
          IF (.not.failure) THEN
             CALL cp_b_dims_create(m_struct%b_dims,&
                  b_row_sizes=b_row_sizes,&
                  b_col_sizes=b_col_sizes,&
                  error=error)
          END IF
       END IF
    END IF
    IF (.not.failure) THEN
       IF (no_spars) THEN
          IF (PRESENT(sparsity)) THEN
             !CPPreconditionNoFail(.not.associated(sparsity),cp_warning_level,routineP,error)
          END IF
          NULLIFY(m_struct%sparsity)
       ELSE IF (PRESENT(sparsity)) THEN
          !if (associated(sparsity)) then
             ALLOCATE(m_struct%sparsity(&
                  SIZE(m_struct%b_dims%block_begins_at_row)-1,&
                  SIZE(m_struct%b_dims%block_begins_at_col)-1),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(m_struct%sparsity,1)==SIZE(sparsity,1),cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(m_struct%sparsity,2)==SIZE(sparsity,2),cp_failure_level,routineP,error,failure)
             IF (.not.failure) m_struct%sparsity=sparsity
          !end if
       ELSE
          ALLOCATE(m_struct%sparsity(&
               SIZE(m_struct%b_dims%block_begins_at_row)-1,&
               SIZE(m_struct%b_dims%block_begins_at_col)-1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.not.failure) THEN
             m_struct%sparsity=.FALSE.
          END IF
       END IF
    END IF
    m_struct%initialized=.not.failure
    CPPostcondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
  END SUBROUTINE cp_bmstruct_create
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bmstruct_dealloc [1.0] *
!!
!!   NAME
!!     cp_bmstruct_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_dealloc(m_struct, error)
!!       Type(cp_b_matrix_struct_type), Pointer:: m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_dealloc
!!
!!   FUNCTION
!!     deallocates the memory allocated by the object and the object
!!     itself
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     m_struct: the object to deallocate
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
  SUBROUTINE cp_bmstruct_dealloc(m_struct,error)
    TYPE(cp_b_matrix_struct_type), POINTER   :: m_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(m_struct)) THEN
       CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
       IF (.not.failure) THEN
          CPPreconditionNoFail(m_struct%ref_count==0,cp_warning_level,routineP,error)
          m_struct%initialized=.FALSE.
          IF (ASSOCIATED(m_struct%b_dims)) THEN
             CALL cp_release(m_struct%b_dims, error)
          END IF
          IF (ASSOCIATED(m_struct%distribution)) THEN
             CALL cp_release(m_struct%distribution,error)
          END IF
          IF (ASSOCIATED(m_struct%sparsity)) THEN
             DEALLOCATE(m_struct%sparsity, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF
       END IF
    END IF
  END SUBROUTINE cp_bmstruct_dealloc
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bmstruct_get [1.0] *
!!
!!   NAME
!!     cp_bmstruct_get
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_get(m_struct, para_env, sparsity,&
!!         sparsity_frozen, symmetric, b_dims, distribution, id_nr,&
!!         ref_count, error)
!!       Type(cp_b_matrix_struct_type), Intent (IN):: m_struct
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Type(cp_b_dims_type), Pointer, Optional:: b_dims
!!       Type(cp_b_distribution_type), Pointer, Optional:: distribution
!!       Logical, Dimension(:,:), Pointer, Optional:: sparsity
!!       Integer, Intent (OUT), Optional:: ref_count, id_nr
!!       Logical, Intent (OUT), Optional:: sparsity_frozen, symmetric
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_get
!!
!!   FUNCTION
!!     returns various attributes of the matrix structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
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
SUBROUTINE cp_bmstruct_get(m_struct, para_env, sparsity, sparsity_frozen,&
     symmetric,&
     b_dims, distribution, id_nr, ref_count, error)
    TYPE(cp_b_matrix_struct_type), &
      INTENT(in)                             :: m_struct
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    LOGICAL, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: sparsity
    LOGICAL, INTENT(out), OPTIONAL           :: sparsity_frozen, symmetric
    TYPE(cp_b_dims_type), OPTIONAL, POINTER  :: b_dims
    TYPE(cp_b_distribution_type), OPTIONAL, &
      POINTER                                :: distribution
    INTEGER, INTENT(out), OPTIONAL           :: id_nr, ref_count
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(sparsity)) sparsity => m_struct%sparsity
     IF (PRESENT(sparsity_frozen)) &
          sparsity_frozen = cp_bmstruct_get_sparsity_frozen(m_struct,error)
     IF (PRESENT(symmetric)) symmetric=m_struct%symmetric
     IF (PRESENT(b_dims)) b_dims => m_struct%b_dims
     IF (PRESENT(distribution)) distribution => m_struct%distribution
     IF (PRESENT(ref_count)) ref_count=m_struct%ref_count
     IF (PRESENT(id_nr)) id_nr=m_struct%id_nr
  END IF
END SUBROUTINE cp_bmstruct_get
!***************************************************************************

!!****f* cp_b_matrix_structureure/cp_bmstruct_freeze_sparsity [1.0] *
!!
!!   NAME
!!     cp_bmstruct_freeze_sparsity
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_freeze_sparsity(matrix_struct, no_sync,&
!!         error)
!!       Type(cp_b_matrix_struct_type), Intent (INOUT):: matrix_struct
!!       Logical, Optional, Intent (IN):: no_sync
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_freeze_sparsity
!!
!!   FUNCTION
!!     freezes and syncs the matrix structure between all the processors.
!!
!!   NOTES
!!     sync should be optimized
!!
!!   INPUTS
!!     matrix_struct: the matrix structure whose sparsity has to be
!!              frozen
!!     no_sync: if true does not sync between the processors
!!              (defaults to false)
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
  SUBROUTINE cp_bmstruct_freeze_sparsity(matrix_struct,no_sync,error)
    TYPE(cp_b_matrix_struct_type), &
      INTENT(inout)                          :: matrix_struct
    LOGICAL, INTENT(in), OPTIONAL            :: no_sync
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_freeze_sparsity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: tmp_distr
    LOGICAL                                  :: failure, sync

    failure=.FALSE.
    sync=.TRUE.

    CALL timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(matrix_struct%initialized,cp_failure_level,routineP,error,failure)
    IF ((.NOT. failure) .AND. (.NOT. matrix_struct%sparsity_frozen)) THEN
       IF (PRESENT(no_sync)) sync=.not.no_sync
       IF (sync) THEN
          WHERE (matrix_struct%distribution%block_owner /= &
               matrix_struct%distribution%para_env%mepos) matrix_struct%sparsity = .FALSE.
          ! would be nice to do this directly with logical (or at least use
          ! the bits of the integer
          ALLOCATE (tmp_distr(SIZE(matrix_struct%sparsity,1),&
               SIZE(matrix_struct%sparsity,2)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          WHERE (matrix_struct%sparsity)
             tmp_distr = 0
          ELSEWHERE
             tmp_distr = 1
          END WHERE
          CALL mp_sum(tmp_distr,matrix_struct%distribution%para_env%group)
          WHERE (tmp_distr == 0)
             matrix_struct%sparsity = .FALSE.
          ELSEWHERE
             matrix_struct%sparsity = .TRUE.
          END WHERE
          DEALLOCATE(tmp_distr,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       matrix_struct%sparsity_frozen=.TRUE.
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_bmstruct_freeze_sparsity
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bmstruct_get_sparsity_frozen [1.0] *
!!
!!   NAME
!!     cp_bmstruct_get_sparsity_frozen
!!
!!   SYNOPSIS
!!     Function cp_bmstruct_get_sparsity_frozen(matrix_struct, error)&
!!         Result(res)
!!       Type(cp_b_matrix_struct_type), Intent (IN):: matrix_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_bmstruct_get_sparsity_frozen
!!
!!   FUNCTION
!!     returns true if the sparsity structure is frozen
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix_struct: the matrix struct you want info about
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
  FUNCTION cp_bmstruct_get_sparsity_frozen(matrix_struct,error) RESULT(res)
    TYPE(cp_b_matrix_struct_type), &
      INTENT(in)                             :: matrix_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_bmstruct_get_sparsity_frozen', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(matrix_struct%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(.FALSE.,cp_warning_level,routineP,error,failure)
       ! to do
    END IF
    res=.FALSE.
    CALL timestop(0.0_dp,handle)
  END FUNCTION cp_bmstruct_get_sparsity_frozen
  !***************************************************************************

!!****f* cp_matrix_structure/cp_bmstruct_retain [1.0] *
!!
!!   NAME
!!     cp_bmstruct_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_retain(m_struct, error)
!!       Type(cp_b_matrix_struct_type), Intent (INOUT):: m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_retain
!!
!!   FUNCTION
!!     augments the retain count by one (to be called if you want a shared 
!!     copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     m_struct: the object to retain
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
  SUBROUTINE cp_bmstruct_retain(m_struct,error)
    TYPE(cp_b_matrix_struct_type), &
      INTENT(inout)                          :: m_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(m_struct%ref_count>0,cp_warning_level,routineP,error)
    IF (.NOT. failure) THEN
       m_struct%ref_count=m_struct%ref_count+1
    END IF
  END SUBROUTINE cp_bmstruct_retain
  !***************************************************************************

!!****f* cp_matrix_structure/cp_bmstruct_release [1.0] *
!!
!!   NAME
!!     cp_bmstruct_release
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_release(m_struct, error)
!!       Type(cp_b_matrix_struct_type), Pointer:: m_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_release
!!
!!   FUNCTION
!!     decreases the ref_count by one (to be called if you do not want anymore
!!     your shared copy of this object)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     m_struct: the object to release
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
  SUBROUTINE cp_bmstruct_release(m_struct,error)
    TYPE(cp_b_matrix_struct_type), POINTER   :: m_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(m_struct)) THEN
       IF (cp_debug) THEN
          CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
          CPPrecondition(m_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.not.failure) THEN
          m_struct%ref_count=m_struct%ref_count-1
          IF (m_struct%ref_count<1) THEN
             CALL cp_dealloc(m_struct,error)
          END IF
       END IF
    END IF
    NULLIFY(m_struct)
  END SUBROUTINE cp_bmstruct_release
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bmstruct_get_block_info [1.0] *
!!
!!   NAME
!!     cp_bmstruct_get_block_info
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_get_block_info(m_struct, block_row, block_col,&
!!         owner, n_local_rows, n_local_cols, begins_at_row, ends_at_row,&
!!         begins_at_col, ends_at_col, out_of_sparsity, me_owner, error)
!!       Type(cp_b_matrix_struct_type), Intent (IN):: m_struct
!!       Integer, Intent (IN):: block_col, block_row
!!       Integer, Intent (OUT), Optional:: owner, n_local_rows,&
!!         n_local_cols, begins_at_row, begins_at_col, ends_at_row,&
!!         ends_at_col
!!       Logical, Intent (OUT), Optional:: me_owner, out_of_sparsity
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bmstruct_get_block_info
!!
!!   FUNCTION
!!     returns information about the given block
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     m_struct: the matrix structure object you get the info from
!!     block_row: the row of the block you are intrested in
!!     block_col: the column of the block you are intrested in
!!     owner: the id of the processor that owns the block
!!     n_local_rows: the number of rows in the block
!!     n_local_cols: the number of columns in the block
!!     begins_at_row: the global row index at which the block begins
!!     ends_at_row: the global row index at which the block ends
!!     begins_at_col: the global column index at which the block begins
!!     ends_at_col: the global column index at which the block ends
!!     out_of_sparsity: if the actual block is outside the sparsity structure
!!     me_owner: true if the actual processor is the owner
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
  SUBROUTINE cp_bmstruct_get_block_info(m_struct,block_row, block_col,&
       owner,n_local_rows,n_local_cols,begins_at_row,ends_at_row,&
       begins_at_col,ends_at_col,out_of_sparsity, me_owner, error)
    TYPE(cp_b_matrix_struct_type), &
      INTENT(in)                             :: m_struct
    INTEGER, INTENT(in)                      :: block_row, block_col
    INTEGER, INTENT(out), OPTIONAL :: owner, n_local_rows, n_local_cols, &
      begins_at_row, ends_at_row, begins_at_col, ends_at_col
    LOGICAL, INTENT(out), OPTIONAL           :: out_of_sparsity, me_owner
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bmstruct_get_block_info', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
    IF (cp_debug) THEN
       CPPrecondition(block_row>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(block_col>0,cp_failure_level,routineP,error,failure)
       ! check also the upper bound?
       CPPrecondition(ASSOCIATED(m_struct%b_dims),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(m_struct%distribution),cp_failure_level,routineP,error,failure)
       IF (.not.failure) THEN
          CPPrecondition(ASSOCIATED(m_struct%b_dims%block_begins_at_col),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(m_struct%b_dims%block_begins_at_row),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(m_struct%distribution%block_owner),cp_failure_level,routineP,error,failure)
       END IF
    END IF
    IF (.NOT. failure) THEN
       IF (PRESENT(owner)) owner=m_struct%distribution%block_owner(block_row,block_col)
       IF (PRESENT(n_local_rows)) THEN
          n_local_rows=m_struct%b_dims%block_begins_at_row(block_row+1)-&
               m_struct%b_dims%block_begins_at_row(block_row)
       END IF
       IF (PRESENT(n_local_cols)) THEN
          n_local_cols=m_struct%b_dims%block_begins_at_col(block_col+1)-&
               m_struct%b_dims%block_begins_at_col(block_col)
       END IF
       IF (PRESENT(begins_at_row)) begins_at_row=m_struct%b_dims%block_begins_at_row(block_row)
       IF (PRESENT(begins_at_col)) begins_at_col=m_struct%b_dims%block_begins_at_col(block_col)
       IF (PRESENT(ends_at_row)) ends_at_row=m_struct%b_dims%block_begins_at_row(block_row+1)-1
       IF (PRESENT(ends_at_col)) ends_at_col=m_struct%b_dims%block_begins_at_col(block_col+1)-1
       IF (PRESENT(out_of_sparsity)) THEN
          IF (.not.m_struct%sparsity_frozen .OR. &
               .not.ASSOCIATED(m_struct%sparsity)) THEN
             out_of_sparsity=.FALSE.
          ELSE
             out_of_sparsity=m_struct%sparsity(block_row,block_col)
          END IF
       END IF
    END IF
  END SUBROUTINE cp_bmstruct_get_block_info
  !***************************************************************************

  !====== block iterator ======

!!****f* cp_b_matrix_structure/cp_bs_iter_init [1.0] *
!!
!!   NAME
!!     cp_bs_iter_init
!!
!!   SYNOPSIS
!!     Subroutine cp_bs_iter_init(bs_iterator, m_struct, do_global, do_zero,&
!!         do_lower_half, error)
!!       Type(cp_bs_iterator_type), Intent (OUT):: bs_iterator
!!       Type(cp_b_matrix_struct_type), Intent (IN), Target:: m_struct
!!       Logical, Intent (IN):: do_global
!!       Logical, Intent (IN), Optional:: do_zero, do_lower_half
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bs_iter_init
!!
!!   FUNCTION
!!     initializes a block structure iterator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     bs_iterator: the iterator that should be initialized
!!     m_struct: the matrix structure it should be iterated on
!!     do_global: if true iterates trough all the blocks, otherwise
!!                does only the owned ones
!!     do_zero: if true does even the blocks outside the sparsity pattern
!!              (defaults to true)
!!     do_lower_half: if true does even the blocks on the lower half of
!!              the matrix, otherwise does only the diagonal and upper 
!!              half (if the mattrix is symmetric defaults to false,
!!              otherwies to true)
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
  SUBROUTINE cp_bs_iter_init(bs_iterator, m_struct, do_global,&
       do_zero, do_lower_half, error)
    TYPE(cp_bs_iterator_type), INTENT(out)   :: bs_iterator
    TYPE(cp_b_matrix_struct_type), &
      INTENT(in), TARGET                     :: m_struct
    LOGICAL, INTENT(in)                      :: do_global
    LOGICAL, INTENT(in), OPTIONAL            :: do_zero, do_lower_half
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bs_iter_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       bs_iterator%do_global=do_global
       bs_iterator%do_zero=.FALSE.
       bs_iterator%do_lower_half=.not.m_struct%symmetric
       IF (PRESENT(do_zero)) bs_iterator%do_zero=do_zero
       IF (PRESENT(do_lower_half)) bs_iterator%do_lower_half=do_lower_half
       bs_iterator%block_row=1
       bs_iterator%block_col=1
       bs_iterator%m_struct => m_struct
       bs_iterator%n_block_rows=SIZE(m_struct%b_dims%block_begins_at_row)
       bs_iterator%n_block_cols=SIZE(m_struct%b_dims%block_begins_at_col)
       bs_iterator%initialized=.TRUE.
    END IF
  END SUBROUTINE cp_bs_iter_init
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bs_iter_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_bs_iter_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_bs_iter_dealloc_ref(bs_iterator, error)
!!       Type(cp_bs_iterator_type), Intent (INOUT):: bs_iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_bs_iter_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by the iterator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     bs_iterator: the iterator that should free its memory
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
  SUBROUTINE cp_bs_iter_dealloc_ref(bs_iterator,error)
    TYPE(cp_bs_iterator_type), INTENT(inout) :: bs_iterator
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bs_iter_dealloc_ref', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(bs_iterator%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(bs_iterator%m_struct)
       bs_iterator%initialized=.FALSE.
    END IF
  END SUBROUTINE cp_bs_iter_dealloc_ref
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bs_iter_next [1.0] *
!!
!!   NAME
!!     cp_bs_iter_next
!!
!!   SYNOPSIS
!!     Function cp_bs_iter_next(bs_iterator, block_row, block_col, owner,&
!!         n_local_rows, n_local_cols, begins_at_row, ends_at_row,&
!!         begins_at_col, ends_at_col, out_of_sparsity, me_owner, error)&
!!         Result(res)
!!       Type(cp_bs_iterator_type), Intent (INOUT):: bs_iterator
!!       Integer, Intent (OUT), Optional:: block_col, block_row, owner,&
!!         n_local_rows, n_local_cols, begins_at_row, begins_at_col,&
!!         ends_at_row, ends_at_col
!!       Logical, Intent (OUT), Optional:: me_owner, out_of_sparsity
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_bs_iter_next
!!
!!   FUNCTION
!!     goes to the next matrix block, returns false if there is no next block.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     bs_iterator: the iterator
!!     block_row: the row of the actual block
!!     block_col: the column of the actual block
!!     owner: the id of the processor that owns the block
!!     n_local_rows: the number of rows in the block
!!     n_local_cols: the number of columns in the block
!!     begins_at_row: the global row index at which the block begins
!!     ends_at_row: the global row index at which the block ends
!!     begins_at_col: the global column index at which the block begins
!!     ends_at_col: the global column index at which the block ends
!!     out_of_sparsity: if the actual block is outside the sparsity structure
!!     me_owner: true if the actual processor is the owner
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
  FUNCTION cp_bs_iter_next(bs_iterator,block_row, block_col,&
       owner,n_local_rows,n_local_cols,begins_at_row,ends_at_row,&
       begins_at_col,ends_at_col,out_of_sparsity, me_owner, error)&
       RESULT(res)
    TYPE(cp_bs_iterator_type), INTENT(inout) :: bs_iterator
    INTEGER, INTENT(out), OPTIONAL :: block_row, block_col, owner, &
      n_local_rows, n_local_cols, begins_at_row, ends_at_row, begins_at_col, &
      ends_at_col
    LOGICAL, INTENT(out), OPTIONAL           :: out_of_sparsity, me_owner
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_bs_iter_next', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (cp_debug) THEN
       CPPrecondition(bs_iterator%initialized,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       bs_iterator%block_col=bs_iterator%block_col+1
       IF (bs_iterator%block_col> bs_iterator%n_block_cols) THEN
          bs_iterator%block_row=bs_iterator%block_row+1
          IF (bs_iterator%block_row<1.or.bs_iterator%block_row> bs_iterator%n_block_rows) THEN
             bs_iterator%block_row=-1
             bs_iterator%block_col= bs_iterator%n_block_cols+1
             IF (PRESENT(block_row)) block_row=-1
             IF (PRESENT(block_col)) block_col=-1
             res=.FALSE.
          ELSE
             IF (bs_iterator%do_lower_half) THEN
                bs_iterator%block_col=1
             ELSE
                bs_iterator%block_col=bs_iterator%block_row
             END IF
          END IF
       END IF
       IF (PRESENT(block_row)) block_row=bs_iterator%block_row
       IF (PRESENT(block_col)) block_col=bs_iterator%block_col
       CALL cp_bmstruct_get_block_info(bs_iterator%m_struct,&
            bs_iterator%block_row, bs_iterator%block_col,&
            owner,n_local_rows,n_local_cols,begins_at_row,ends_at_row,&
            begins_at_col,ends_at_col,out_of_sparsity, me_owner, error)
    END IF
  END FUNCTION cp_bs_iter_next
  !***************************************************************************

END MODULE cp_b_matrix_structure
