!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sm_struct [1.0] *
!!
!!   NAME
!!     cp_sm_struct
!!
!!   FUNCTION
!!     represent the structure of a block matrix (sparse matrix)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_sm_struct
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
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
                                             cp_para_env_retain,&
                                             cp_para_env_write
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_sm_struct'
  INTEGER, SAVE, PRIVATE :: last_sm_struct_id=0

  PUBLIC :: cp_sm_node_type, cp_sm_struct_type, cp_sm_struct_p_type
  PUBLIC :: sm_struct_create, sm_struct_retain, sm_struct_release,&
       sm_struct_add_local_block, sm_struct_clear_local_blocks,&
       sm_struct_get, sm_struct_write
!!***
!****************************************************************************

!!****s* cp_sm_struct/cp_sm_node_type [1.0] *
!!
!!   NAME
!!     cp_sm_node_type
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - pos: position of the actual block (in 1d array)
!!     - row: row of the block
!!     - col: column of the block
!!     - next_row: next block in the same col
!!     - next_col: next block in the same row
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
  !**************************************************************************
  TYPE cp_sm_node_type
     INTEGER :: pos, row, col
     TYPE(cp_sm_node_type), POINTER :: next_row, next_col
  END TYPE cp_sm_node_type
!!***
  !**************************************************************************

!!****s* cp_sm_struct/cp_sm_struct_type [1.0] *
!!
!!   NAME
!!     cp_sm_struct_type
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     split dimension in another structure?
!!
!!   ATTRIBUTES
!!     - n_local_blocks: number of local blocks
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
  !**************************************************************************
  TYPE cp_sm_struct_type
     INTEGER :: ref_count, id_nr, n_local_blocks
     INTEGER, DIMENSION(:), POINTER :: block_begins_at_row, &
          block_begins_at_col, first_col, &
          first_row, last_col, last_row
     TYPE(cp_para_env_type), POINTER :: para_env
     CHARACTER(len=40) :: symmetry
     TYPE(cp_sm_node_type), POINTER :: local_blocks
     TYPE(distribution_2d_type), POINTER :: distribution_2d
  END TYPE cp_sm_struct_type
!!***
  !**************************************************************************

!!****s* cp_sm_struct_p_type [1.0] *
!!
!!   NAME
!!     cp_sm_struct_p_type
!!
!!   FUNCTION
!!     just to build array of sm_struct pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
  !**************************************************************************
  TYPE cp_sm_struct_p_type
     TYPE(cp_sm_struct_type), POINTER :: sm_struct
  END TYPE cp_sm_struct_p_type
!!***
  !**************************************************************************

CONTAINS

!!****f* cp_sm_struct/sm_struct_create [1.0] *
!!
!!   NAME
!!     sm_struct_create
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_create(sm_struct, block_begins_at_row,&
!!         block_begins_at_col, para_env, symmetry, distribution_2d,&
!!         local_blocks_ptr, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Integer, Dimension(:), Intent (IN):: block_begins_at_row,&
!!         block_begins_at_col
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Character(Len=*):: symmetry
!!       Type(distribution_2d_type), Pointer:: distribution_2d
!!       Type(cp_sm_node_type), Optional, Pointer:: local_blocks_ptr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_create(sm_struct, block_begins_at_row,&
     block_begins_at_col, para_env, symmetry, distribution_2d,&
     local_blocks_ptr,error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    INTEGER, DIMENSION(:), INTENT(in)        :: block_begins_at_row, &
                                                block_begins_at_col
    TYPE(cp_para_env_type), POINTER          :: para_env
    CHARACTER(len=*)                         :: symmetry
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_sm_node_type), OPTIONAL, POINTER :: local_blocks_ptr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(sm_struct, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     sm_struct%ref_count=1
     last_sm_struct_id=last_sm_struct_id+1
     sm_struct%id_nr=last_sm_struct_id
     ALLOCATE(sm_struct%block_begins_at_row(SIZE(block_begins_at_row)), &
          sm_struct%block_begins_at_col(SIZE(block_begins_at_col)),&
          stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     sm_struct%block_begins_at_row=block_begins_at_row
     sm_struct%block_begins_at_col=block_begins_at_col
     NULLIFY(sm_struct%first_col, sm_struct%first_row,&
          sm_struct%last_col, sm_struct%last_row)
     sm_struct%para_env => para_env
     CALL cp_para_env_retain(para_env,error=error)
     sm_struct%symmetry=symmetry
     sm_struct%distribution_2d => distribution_2d
     IF (ASSOCIATED(distribution_2d)) THEN
        CALL distribution_2d_retain(distribution_2d,error=error)
     END IF
     IF (PRESENT(local_blocks_ptr)) THEN
        sm_struct%local_blocks => local_blocks_ptr
     ELSE
        NULLIFY(sm_struct%local_blocks)
     END IF
     CALL sm_nodes_renumber(sm_struct%local_blocks,&
          nblocks=sm_struct%n_local_blocks, error=error)
  END IF
END SUBROUTINE sm_struct_create
!***************************************************************************

!!****f* cp_sm_struct/sm_struct_retain [1.0] *
!!
!!   NAME
!!     sm_struct_retain
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_retain(sm_struct, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_retain
!!
!!   FUNCTION
!!     retains the given sm_struct (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_retain(sm_struct,error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(sm_struct),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(sm_struct%ref_count>0,cp_failure_level,routineP,error,failure)
     sm_struct%ref_count=sm_struct%ref_count+1
  END IF
END SUBROUTINE sm_struct_retain
!***************************************************************************

!!****f* cp_sm_struct/sm_nodes_dealloc [1.0] *
!!
!!   NAME
!!     sm_nodes_dealloc
!!
!!   FUNCTION
!!     deallocates the whole node structure staring at sm_nodes
!!
!!   NOTES
!!     assumes that the first column has pointers to each row
!!     private routine
!!
!!   ARGUMENTS
!!     - sm_nodes: the nodes to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_nodes_dealloc(sm_nodes,error)
    TYPE(cp_sm_node_type), POINTER           :: sm_nodes
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_nodes_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_node_type), POINTER           :: next_pos, next_row, pos_att

  failure=.FALSE.
  
  IF (ASSOCIATED(sm_nodes)) THEN
     NULLIFY(next_row, pos_att)
     pos_att => sm_nodes
     DO WHILE (ASSOCIATED(pos_att))
        next_row => pos_att%next_row
        DO WHILE (ASSOCIATED(pos_att))
           next_pos => pos_att%next_col
           DEALLOCATE(pos_att, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           pos_att => next_pos
        END DO
        pos_att => next_row
     END DO
     NULLIFY(sm_nodes)
  END IF
END SUBROUTINE sm_nodes_dealloc
!***************************************************************************

!!****f* cp_sm_struct/sm_struct_release [1.0] *
!!
!!   NAME
!!     sm_struct_release
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_release(sm_struct, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_release
!!
!!   FUNCTION
!!     releases the given sm_struct (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_release(sm_struct, error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(sm_struct)) THEN
     CPPreconditionNoFail(sm_struct%ref_count>0,cp_failure_level,routineP,error)
     sm_struct%ref_count=sm_struct%ref_count-1
     IF (sm_struct%ref_count==0) THEN
        IF (ASSOCIATED(sm_struct%block_begins_at_row)) THEN
           DEALLOCATE(sm_struct%block_begins_at_row,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(sm_struct%block_begins_at_col)) THEN
           DEALLOCATE(sm_struct%block_begins_at_col,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(sm_struct%first_col)) THEN
           DEALLOCATE(sm_struct%first_col,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(sm_struct%first_row)) THEN
           DEALLOCATE(sm_struct%first_row,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(sm_struct%last_row)) THEN
           DEALLOCATE(sm_struct%last_row,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(sm_struct%last_col)) THEN
           DEALLOCATE(sm_struct%last_col,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        CALL cp_para_env_release(sm_struct%para_env,error=error)
        CALL sm_nodes_dealloc(sm_struct%local_blocks, error=error)
        CALL distribution_2d_release(sm_struct%distribution_2d,error=error)

        DEALLOCATE(sm_struct, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE sm_struct_release
!***************************************************************************

!!****f* cp_sm_struct/sm_struct_add_local_block [1.0] *
!!
!!   NAME
!!     sm_struct_add_local_block
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_add_local_block(sm_struct, block_row, block_col,&
!!         error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Integer, Intent (IN):: block_row, block_col
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_add_local_block
!!
!!   FUNCTION
!!     adds a local block to the given row and col
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the structure where to add the block
!!     - block_row: the row of the block to add
!!     - block_col: the column of the block to add
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_add_local_block(sm_struct, block_row, block_col, error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    INTEGER, INTENT(in)                      :: block_row, block_col
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_add_local_block', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_node_type), POINTER           :: col_att, new_node, row_att

  failure=.FALSE.
  NULLIFY(row_att, col_att, new_node)
  
  CPPrecondition(ASSOCIATED(sm_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(sm_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (.not.ASSOCIATED(sm_struct%local_blocks)) THEN
        ALLOCATE(sm_struct%local_blocks,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        sm_struct%local_blocks%row=0
        sm_struct%local_blocks%col=0
        sm_struct%local_blocks%pos=0
        NULLIFY(sm_struct%local_blocks%next_row,&
             sm_struct%local_blocks%next_col)
     END IF
  END IF
  IF (.not.failure) THEN
     row_att => sm_struct%local_blocks
     CPAssert(ASSOCIATED(row_att),cp_failure_level,routineP,error,failure)
     DO WHILE (.TRUE.)
        IF (.not.ASSOCIATED(row_att%next_row)) EXIT
        IF (row_att%next_row%row > block_row) EXIT
        row_att => row_att%next_row
     END DO
     CPPostcondition(ASSOCIATED(row_att),cp_failure_level,routineP,error,failure)
     CPPostcondition(row_att%row<=block_row,cp_failure_level,routineP,error,failure)
     IF (row_att%row /= block_row) THEN
        ALLOCATE(new_node,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        new_node%row=block_row
        new_node%col=0
        new_node%pos=0
        new_node%next_row => row_att%next_row
        row_att%next_row => new_node
        NULLIFY(new_node%next_col)
        row_att => new_node
     END IF
  END IF
  IF (.NOT. failure) THEN
     DO WHILE (.TRUE.)
        IF (.not.ASSOCIATED(row_att%next_col)) EXIT
        IF (row_att%next_col%col > block_col) EXIT
        row_att => row_att%next_col
     END DO
  END IF

  IF (.NOT. failure) THEN
     col_att => sm_struct%local_blocks
     CPAssert(ASSOCIATED(col_att),cp_failure_level,routineP,error,failure)
     DO WHILE (.TRUE.)
        IF (.not.ASSOCIATED(col_att%next_col)) EXIT
        IF (col_att%next_col%col > block_col) EXIT
        col_att => col_att%next_col
     END DO
     CPPostcondition(ASSOCIATED(col_att),cp_failure_level,routineP,error,failure)
     CPPostcondition(col_att%col<=block_col,cp_failure_level,routineP,error,failure)
     IF (col_att%col /= block_col) THEN
        ALLOCATE(new_node,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        new_node%col=block_col
        new_node%row=0
        new_node%pos=0
        new_node%next_col => col_att%next_col
        col_att%next_col => new_node
        NULLIFY(new_node%next_row)
        col_att => new_node
     END IF
  END IF
  IF (.NOT. failure) THEN
     DO WHILE (.TRUE.)
        IF (.not.ASSOCIATED(col_att%next_row)) EXIT
        IF (col_att%next_row%row > block_row) EXIT
        col_att => col_att%next_row
     END DO
  END IF

  IF (.not.failure) THEN
     IF (col_att%row == block_row) THEN
        CPAssert(row_att%col == block_col,cp_failure_level,routineP,error,failure)
     ELSE
        CPAssert(row_att%col /= block_col,cp_failure_level,routineP,error,failure)
        ALLOCATE(new_node,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        new_node%col=block_col
        new_node%row=block_row
        sm_struct%n_local_blocks=sm_struct%n_local_blocks+1
        new_node%pos=sm_struct%n_local_blocks
        new_node%next_col => row_att%next_col
        row_att%next_col => new_node
        new_node%next_row => col_att%next_row
        col_att%next_row => new_node
     END IF
  END IF

END SUBROUTINE sm_struct_add_local_block
!***************************************************************************

!!****f* cp_sm_struct/sm_nodes_renumber [1.0] *
!!
!!   NAME
!!     sm_nodes_renumber
!!
!!   FUNCTION
!!     renumbers the nodes of the structure
!!
!!   NOTES
!!     private routine
!!
!!   ARGUMENTS
!!     - sm_nodes: the nodes to renumber
!!     - nblocks: will contain the number of block (ignores 0 col and row)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_nodes_renumber(sm_nodes, nblocks, error)
    TYPE(cp_sm_node_type), POINTER           :: sm_nodes
    INTEGER, INTENT(out), OPTIONAL           :: nblocks
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_nodes_renumber', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nr
    LOGICAL                                  :: failure
    TYPE(cp_sm_node_type), POINTER           :: next_row, pos_att

  failure=.FALSE.
  nr=0
  
  IF (ASSOCIATED(sm_nodes)) THEN
     NULLIFY(next_row, pos_att)
     pos_att => sm_nodes
     DO WHILE (ASSOCIATED(pos_att))
        next_row => pos_att%next_row
        DO WHILE (ASSOCIATED(pos_att))
           IF (pos_att%row>0 .AND. pos_att%col>0) THEN
              nr=nr+1
              pos_att%pos=nr
           ELSE
              pos_att%pos=0
           END IF
           pos_att => pos_att%next_col
        END DO
        pos_att => next_row
     END DO
  END IF
  
  IF (PRESENT(nblocks)) nblocks=nr

END SUBROUTINE sm_nodes_renumber
!***************************************************************************

!!****f* cp_sm_struct/sm_struct_clear_local_blocks [1.0] *
!!
!!   NAME
!!     sm_struct_clear_local_blocks
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_clear_local_blocks(sm_struct, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_clear_local_blocks
!!
!!   FUNCTION
!!     clears the local blocks
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_clear_local_blocks(sm_struct, error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_clear_local_blocks', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(sm_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(sm_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL sm_nodes_dealloc(sm_struct%local_blocks,error=error)
     sm_struct%n_local_blocks=0
  END IF
END SUBROUTINE sm_struct_clear_local_blocks
!***************************************************************************

!!****f* cp_sm_struct/sm_struct_get [1.0] *
!!
!!   NAME
!!     sm_struct_get
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_get(sm_struct, ref_count, id_nr, n_local_blocks,&
!!         block_begins_at_row, block_begins_at_col, first_col, first_row,&
!!         last_col, last_row, para_env, symmetry, local_blocks,&
!!         n_blocks_row, n_blocks_col, n_rows, n_cols, distribution_2d,&
!!         error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Integer, Intent (OUT), Optional:: ref_count, id_nr, n_local_blocks
!!       Integer, Dimension(:), Pointer, Optional:: block_begins_at_row,&
!!         block_begins_at_col, first_col, first_row, last_col, last_row
!!       Type(cp_para_env_type), Pointer, Optional:: para_env
!!       Character(Len=*), Intent (OUT), Optional:: symmetry
!!       Type(cp_sm_node_type), Pointer, Optional:: local_blocks
!!       Integer, Intent (OUT), Optional:: n_blocks_row, n_blocks_col,&
!!         n_rows, n_cols
!!       Type(distribution_2d_type), Pointer, Optional:: distribution_2d
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_get
!!
!!   FUNCTION
!!     returns various attributes about the given sm_struct
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the structure you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_get(sm_struct, ref_count, id_nr, n_local_blocks,&
     block_begins_at_row, block_begins_at_col, first_col, &
     first_row, last_col, last_row, para_env, symmetry, local_blocks,&
     n_blocks_row,n_blocks_col,n_rows,n_cols,distribution_2d, error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    INTEGER, INTENT(out), OPTIONAL           :: ref_count, id_nr, &
                                                n_local_blocks
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: block_begins_at_row, &
                                                block_begins_at_col, &
                                                first_col, first_row, &
                                                last_col, last_row
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    CHARACTER(len=*), INTENT(out), OPTIONAL  :: symmetry
    TYPE(cp_sm_node_type), OPTIONAL, POINTER :: local_blocks
    INTEGER, INTENT(out), OPTIONAL           :: n_blocks_row, n_blocks_col, &
                                                n_rows, n_cols
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(sm_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(sm_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(ref_count)) ref_count=sm_struct%ref_count
     IF (PRESENT(id_nr)) id_nr=sm_struct%id_nr
     IF (PRESENT(n_local_blocks)) n_local_blocks=sm_struct%n_local_blocks
     IF (PRESENT(block_begins_at_col)) &
          block_begins_at_col => sm_struct%block_begins_at_col
     IF (PRESENT(block_begins_at_row)) &
          block_begins_at_row => sm_struct%block_begins_at_row
     IF (PRESENT(first_col)) THEN !FM avoid allocation? (use subset)
        IF (.not.ASSOCIATED(sm_struct%first_col)) THEN
           ALLOCATE(sm_struct%first_col(SIZE(sm_struct%block_begins_at_col)-1),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           sm_struct%first_col = sm_struct%block_begins_at_col(1:&
                SIZE(sm_struct%first_col))
        END IF
        first_col => sm_struct%first_col
     END IF
     IF (PRESENT(first_row)) THEN !FM avoid allocation? (use subset)
        IF (.not.ASSOCIATED(sm_struct%first_row)) THEN
           ALLOCATE(sm_struct%first_row(SIZE(sm_struct%block_begins_at_row)-1),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           sm_struct%first_row = sm_struct%block_begins_at_row(1:&
                SIZE(sm_struct%first_row))
        END IF
        first_row => sm_struct%first_row
     END IF
     IF (PRESENT(last_row)) THEN
        IF (.not.ASSOCIATED(sm_struct%last_row)) THEN
           ALLOCATE(sm_struct%last_row(SIZE(sm_struct%block_begins_at_row)-1),&
                stat=stat)
           sm_struct%last_row=sm_struct%block_begins_at_row(2:)-1
        END IF
        last_row => sm_struct%last_row
     END IF
     IF (PRESENT(last_col)) THEN
        IF (.not.ASSOCIATED(sm_struct%last_col)) THEN
           ALLOCATE(sm_struct%last_col(SIZE(sm_struct%block_begins_at_col)-1),&
                stat=stat)
           sm_struct%last_col=sm_struct%block_begins_at_col(2:)-1
        END IF
        last_col => sm_struct%last_col
     END IF
     IF (PRESENT(para_env)) THEN
        para_env => sm_struct%para_env
     END IF
     IF (PRESENT(symmetry)) THEN
        symmetry = sm_struct%symmetry
     END IF
     IF (PRESENT(local_blocks)) local_blocks => sm_struct%local_blocks
     IF (PRESENT(n_blocks_row))n_blocks_row=SIZE(sm_struct%block_begins_at_row)-1
     IF (PRESENT(n_blocks_col))n_blocks_col=SIZE(sm_struct%block_begins_at_col)-1
     IF (PRESENT(n_rows)) THEN
        n_rows=sm_struct%block_begins_at_row(&
             SIZE(sm_struct%block_begins_at_row))-1
     END IF
     IF (PRESENT(n_cols)) THEN
        n_cols=sm_struct%block_begins_at_col(&
             SIZE(sm_struct%block_begins_at_col))-1
     END IF
     IF (PRESENT(distribution_2d)) distribution_2d => sm_struct%distribution_2d
  END IF
END SUBROUTINE sm_struct_get
!***************************************************************************

!!****f* cp_sm_struct/sm_struct_write [1.0] *
!!
!!   NAME
!!     sm_struct_write
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_write(sm_struct, unit_nr, local,&
!!         short_description, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: local, short_description
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_write
!!
!!   FUNCTION
!!     writes out information about the given sm_struct
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the sm struct to output
!!     - unit_nr: the unit where to output
!!     - local: if the unit is local to this processor (defaults to false)
!!     - short_description: if only a short description should be written 
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE sm_struct_write(sm_struct, unit_nr, local, short_description, error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: local, short_description
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'sm_struct_write', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_local, &
                                                my_short_description, &
                                                should_write
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sm_node_type), POINTER           :: next_row, pos_att

  failure=.FALSE.
  NULLIFY(pos_att,next_row)
  logger => cp_error_get_logger(error)
  my_local=.FALSE.
  IF (PRESENT(local)) my_local=local
  my_short_description=.TRUE.
  IF (PRESENT(short_description)) my_short_description=short_description
  should_write=local.or. logger%para_env%mepos == logger%para_env%source

  IF (should_write) THEN
     IF (ASSOCIATED(sm_struct)) THEN
        WRITE(unit=unit_nr,fmt="('<sm_struct>={ id_nr=',i6,', ref_count=',i6,"//&
             "', n_local_blocks=',i6,',')")&
             sm_struct%id_nr, sm_struct%ref_count, sm_struct%n_local_blocks
        WRITE(unit=unit_nr,fmt="(' symmetry=',a,',')")&
             sm_struct%symmetry
        IF (my_short_description) THEN

           IF (ASSOCIATED(sm_struct%block_begins_at_row)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' block_begins_at_row=(',i6,' elements),')")&
                   SIZE(sm_struct%block_begins_at_row)
           ELSE
              WRITE(unit=unit_nr,fmt="(' block_begins_at_row=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%block_begins_at_col)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' block_begins_at_col=(',i6,' elements),')")&
                   SIZE(sm_struct%block_begins_at_col)
           ELSE
              WRITE(unit=unit_nr,fmt="(' block_begins_at_col=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%first_col)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' first_col=(',i6,' elements),')")&
                   SIZE(sm_struct%first_col)
           ELSE
              WRITE(unit=unit_nr,fmt="(' first_col=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%first_row)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' first_row=(',i6,' elements),')")&
                   SIZE(sm_struct%first_row)
           ELSE
              WRITE(unit=unit_nr,fmt="(' first_row=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%last_col)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' last_col=(',i6,' elements),')")&
                   SIZE(sm_struct%last_col)
           ELSE
              WRITE(unit=unit_nr,fmt="(' last_col=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%last_row)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' last_row=(',i6,' elements),')")&
                   SIZE(sm_struct%last_row)
           ELSE
              WRITE(unit=unit_nr,fmt="(' last_row=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%para_env)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' para_env=*associated*,')")
!!              &
!!                   SIZE(sm_struct%last_row)
           ELSE
              WRITE(unit=unit_nr,fmt="(' para_env=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%local_blocks)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' local_blocks=*associated*,')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' local_blocks=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%distribution_2d)) THEN
              WRITE(unit=unit_nr,fmt="(' distribution_2d=*associated*,')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' distribution_2d=*null*,')")
           END IF

        ELSE ! long description

           IF (ASSOCIATED(sm_struct%block_begins_at_row)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' block_begins_at_row=(')", advance='no')
              WRITE(unit=unit_nr, fmt="(i6,',')", advance='no')&
                   sm_struct%block_begins_at_row
              WRITE(unit=unit_nr,fmt="(' ),')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' block_begins_at_row=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%block_begins_at_col)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' block_begins_at_col=(')", advance='no')
              WRITE(unit=unit_nr, fmt="(i6,',')", advance='no')&
                   sm_struct%block_begins_at_col
              WRITE(unit=unit_nr,fmt="(' ),')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' block_begins_at_col=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%first_col)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' first_col=(')", advance='no')
              WRITE(unit=unit_nr, fmt="(i6,',')", advance='no')&
                   sm_struct%first_col
              WRITE(unit=unit_nr,fmt="(' ),')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' first_col=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%first_row)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' first_row=(')", advance='no')
              WRITE(unit=unit_nr, fmt="(i6,',')", advance='no')&
                   sm_struct%first_row
              WRITE(unit=unit_nr,fmt="(' ),')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' first_row=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%last_col)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' last_col=(')", advance='no')
              WRITE(unit=unit_nr, fmt="(i6,',')", advance='no')&
                   sm_struct%last_col
              WRITE(unit=unit_nr,fmt="(' ),')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' last_col=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%last_row)) THEN
              WRITE(unit=unit_nr,&
                   fmt="(' last_row=(')", advance='no')
              WRITE(unit=unit_nr, fmt="(i6,',')", advance='no')&
                   sm_struct%last_row
              WRITE(unit=unit_nr,fmt="(' ),')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' last_row=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%para_env)) THEN
              WRITE(unit=unit_nr,fmt="(' para_env=')", advance='no')
              CALL cp_para_env_write(sm_struct%para_env,unit_nr=unit_nr,&
                   error=error)
           ELSE
              WRITE(unit=unit_nr,fmt="(' para_env=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%local_blocks)) THEN

              WRITE(unit=unit_nr, fmt="(' local_blocks(rows{i,j,pos})=(',/,' ')",&
                   advance='no')
              pos_att => sm_struct%local_blocks
              DO WHILE (ASSOCIATED(pos_att))
                 next_row => pos_att%next_row
                 DO WHILE (ASSOCIATED(pos_att))
                    WRITE(unit=unit_nr,fmt="(' {',i6,',',i6,',',i6,'}')",&
                         advance='no')&
                         pos_att%row, pos_att%col, pos_att%pos
                    pos_att => pos_att%next_col
                 END DO
                 WRITE (unit=unit_nr,fmt="(/,' ')",advance='no')
                 pos_att => next_row
              END DO
              WRITE (unit=unit_nr,fmt="(')')")

              IF (debug_this_module) THEN
                 WRITE(unit=unit_nr, fmt="(' local_blocks(cols{i,j,pos})=(',/,' ')",&
                      advance='no')
                 pos_att => sm_struct%local_blocks
                 DO WHILE (ASSOCIATED(pos_att))
                    next_row => pos_att%next_col
                    DO WHILE (ASSOCIATED(pos_att))
                       WRITE(unit=unit_nr,fmt="(' {',i6,',',i6,',',i6,'}')",&
                            advance='no')&
                            pos_att%row, pos_att%col, pos_att%pos
                       pos_att => pos_att%next_row
                    END DO
                    WRITE (unit=unit_nr,fmt="(/,' ')", advance='no')
                    pos_att => next_row
                 END DO
                 WRITE (unit=unit_nr,fmt="(')')")
              END IF
           ELSE
              WRITE(unit=unit_nr,fmt="(' local_blocks=*null*,')")
           END IF

           IF (ASSOCIATED(sm_struct%distribution_2d)) THEN
              WRITE(unit=unit_nr,fmt="(' distribution_2d=*associated*,')")
           ELSE
              WRITE(unit=unit_nr,fmt="(' distribution_2d=*null*,')")
           END IF

        END IF
        WRITE(unit=unit_nr,fmt="('}')")
     ELSE
        WRITE(unit=unit_nr,fmt="(' <cp_sm_struct>=*null*,')")
     END IF
  END IF
END SUBROUTINE sm_struct_write
!***************************************************************************

END MODULE cp_sm_struct
