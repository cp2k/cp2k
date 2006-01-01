!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/cp_matrix_utils [1.0] *
!!
!!   NAME
!!     cp_matrix_utils
!!
!!   FUNCTION
!!     various utilities to do operations with sparse matrices
!!     and 2d array seen as dense matrices
!!
!!   NOTES
!!     this should depend only on sparse_matrix_types
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_matrix_utils
  USE cp_array_utils,                  ONLY: cp_2d_i_guarantee_size,&
                                             cp_2d_logical_guarantee_size,&
                                             cp_2d_r_write
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_struct,                    ONLY: cp_sm_struct_type,&
                                             sm_struct_add_local_block
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sum
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             find_block_node,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_matrix_utils'
! types
  PUBLIC :: cp_matrix_block_iterator

  ! underliying methods
  PUBLIC :: cp_sm_b_i_init,&
            cp_sm_b_i_get,&
            cp_sm_b_i_next,&
            cp_sm_write,&
            sm_struct_add_sm_local_blocks


!!****s* cp_matrix_utils/cp_matrix_block_iterator [1.0] *
!!
!!   NAME
!!     cp_matrix_block_iterator
!!
!!   FUNCTION
!!     goes through all the local blocks of the given matrix
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     initialized: true if the structure was initialized (bug catching)
!!     matrix: the actual matrix
!!     block_row: the actual block row (-1 if past end, -2 in case
!!            of error, 0 after a dealloc)
!!     block_node: the actual block
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  TYPE cp_matrix_block_iterator
     PRIVATE
     LOGICAL :: initialized
     TYPE(real_matrix_type), POINTER :: matrix
     INTEGER :: block_row
     TYPE(real_block_node_type), POINTER :: block_node
  END TYPE cp_matrix_block_iterator
!!***
!****************************************************************************

CONTAINS

!!****f* cp_matrix_utils/cp_sm_b_i_init [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_init
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_b_i_init(block_iterator, matrix, block_row,&
!!         block_node, error)
!!       Type(cp_matrix_block_iterator), Intent (OUT):: block_iterator
!!       Type(real_matrix_type), Target, Intent (IN):: matrix
!!       Integer, Intent (IN), Optional:: block_row
!!       Type(real_block_node_type), Target, Optional:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_b_i_init
!!
!!   FUNCTION
!!     initializes the block iterator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the block iterator to be initialized
!!     matrix: the matrix you iterate on
!!     block_row: the initial block row (defaults to 0, before the first)
!!     block_node: the initial block (defaults to the first of the row)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sm_b_i_init(block_iterator, matrix, block_row,&
       block_node, error)
    TYPE(cp_matrix_block_iterator), &
      INTENT(out)                            :: block_iterator
    TYPE(real_matrix_type), pointer          :: matrix
    INTEGER, INTENT(in), OPTIONAL            :: block_row
    TYPE(real_block_node_type), OPTIONAL, &
         POINTER                                 :: block_node
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_b_i_init'

    LOGICAL                                  :: failure

! well subsequent operations on the iterator could change tha data in the matrix

    failure=.FALSE.

    block_iterator%matrix => matrix
    block_iterator%block_row = 1
    NULLIFY(block_iterator%block_node)
    block_iterator%initialized=.TRUE.
    IF (PRESENT(block_node)) &
         block_iterator%block_node => block_node
    IF (PRESENT(block_row)) &
         block_iterator%block_row = block_row
  END SUBROUTINE cp_sm_b_i_init
!***************************************************************************


!!****f* cp_matrix_utils/cp_sm_b_i_get [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_get
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_b_i_get(block_iterator, block_row, block_col,&
!!         block_val, finished, matrix, block_node, error)
!!       Type(cp_matrix_block_iterator), Intent (IN):: block_iterator
!!       Integer, Optional, Intent (OUT):: block_row, block_col
!!       Real(KIND = dp), Dimension(:,:), Pointer, Optional:: block_val
!!       Logical, Optional, Intent (OUT):: finished
!!       Type(real_matrix_type), Optional, Pointer:: matrix
!!       Type(real_block_node_type), Optional, Pointer:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_b_i_get
!!
!!   FUNCTION
!!     returns the vaious attrubutes of the iterator
!!
!!   NOTES
!!     also att global indexig (firstRow, firstCol)??
!!
!!   INPUTS
!!     block_iterator: the block iterator you get the info from
!!     block_row: the actual block row
!!     block_col: the actual block column
!!     block_val: the values contained in the block (can be modified, 
!!           and the matrix values are modified)
!!     finished: true if the iterator is at end
!!     matrix: the matrix you are iterating on
!!     block_node: the block node (try not to use this)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE cp_sm_b_i_get(block_iterator, block_row, block_col,&
       block_val, finished, matrix, block_node, error)
    TYPE(cp_matrix_block_iterator), &
      INTENT(in)                             :: block_iterator
    INTEGER, INTENT(out), OPTIONAL           :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: block_val
    LOGICAL, INTENT(out), OPTIONAL           :: finished
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix
    TYPE(real_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_b_i_get'

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(block_iterator%matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(block_row)) &
            block_row=block_iterator%block_row
       IF (PRESENT(block_col)) THEN
          IF (ASSOCIATED(block_iterator%block_node)) THEN
             CALL get_block_node(block_iterator%block_node, block_col=block_col)
          ELSE
             CPAssert(.FALSE.,cp_warning_level,routineP,error,failure)
             block_col=-1
          END IF
       END IF
       IF (PRESENT(block_val)) THEN
          CPPrecondition(block_iterator%block_row>0,cp_failure_level,routineP,error,failure)
          IF (.not.failure) THEN
             IF (ASSOCIATED(block_iterator%block_node)) THEN
                CALL get_block_node(block_iterator%block_node,BLOCK=block_val)
             ELSE
                NULLIFY(block_val)
                CPErrorMessage(cp_failure_level,routineP,"invalid BLOCK",error)
             END IF
          END IF
       END IF
       IF (PRESENT(finished)) &
            finished = block_iterator%block_row <= 0
       IF (PRESENT(matrix)) &
            matrix => block_iterator%matrix
       IF (PRESENT(block_node)) &
            block_node => block_iterator%block_node
    END IF
  END SUBROUTINE cp_sm_b_i_get
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_b_i_next [1.0] *
!!
!!   NAME
!!     cp_sm_b_i_next
!!
!!   SYNOPSIS
!!     Function cp_sm_b_i_next(block_iterator, block_row, block_col,&
!!         block_val, finished, matrix, block_node, error) Result(res)
!!       Type(cp_matrix_block_iterator), Intent (INOUT):: block_iterator
!!       Integer, Optional, Intent (OUT):: block_row, block_col
!!       Real(KIND = dp), Dimension(:,:), Pointer, Optional:: block_val
!!       Logical, Optional, Intent (OUT):: finished
!!       Type(real_matrix_type), Optional, Pointer:: matrix
!!       Type(real_block_node_type), Optional, Pointer:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Logical:: res
!!     End Function cp_sm_b_i_next
!!
!!   FUNCTION
!!     Moves the iterator to the next matrix block.
!!     Returns true if the iterator is valid (not gone past end)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - block_iterator: the block iterator that you want to move
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!     the other arguments are the same as the get function, and
!!     are valid only if the result is true
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 1.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  FUNCTION cp_sm_b_i_next(block_iterator, block_row, block_col,&
       block_val, finished, matrix, block_node, error)RESULT(res)
    TYPE(cp_matrix_block_iterator), &
      INTENT(inout)                          :: block_iterator
    INTEGER, INTENT(out), OPTIONAL           :: block_row, block_col
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: block_val
    LOGICAL, INTENT(out), OPTIONAL           :: finished
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: matrix
    TYPE(real_block_node_type), OPTIONAL, &
      POINTER                                :: block_node
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineP = 'cp_matrix_utils:cp_sm_b_i_next'

    INTEGER                                  :: i, nblock_row
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(block_iterator%matrix),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (block_iterator%block_row > 0) THEN
          IF (ASSOCIATED(block_iterator%block_node)) THEN
             block_iterator%block_node => next_block_node(block_iterator%block_node)
          ELSE
             block_iterator%block_node => first_block_node(block_iterator%matrix,&
                  block_iterator%block_row)
          END IF
          IF (.not.ASSOCIATED(block_iterator%block_node)) THEN
             CALL get_matrix_info(block_iterator%matrix,nblock_row=nblock_row)
             DO i=block_iterator%block_row+1,nblock_row
                block_iterator%block_node => &
                     first_block_node(block_iterator%matrix,i)
                IF (ASSOCIATED(block_iterator%block_node)) THEN
                   block_iterator%block_row=i
                   EXIT
                END IF
             END DO
          END IF
          IF (.not.ASSOCIATED(block_iterator%block_node)) THEN
             block_iterator%block_row=-1
          ELSE
             res=.TRUE.
             CALL cp_sm_b_i_get(block_iterator, block_row=block_row,&
                  block_col=block_col, block_val=block_val,&
                  finished=finished, matrix=matrix, block_node=block_node,&
                  error=error)
             RETURN
          END IF
       END IF
       CPPostcondition(.not.ASSOCIATED(block_iterator%block_node),cp_warning_level,routineP,error,failure)
       IF (.NOT. block_iterator%block_row <= 0) THEN
          CPErrorMessage(cp_warning_level,routineP,"inconsistent state in iterator",error)
          block_iterator%block_row=-2
       END IF
    ELSE
       block_iterator%block_row=-2
       NULLIFY(block_iterator%block_node)
    END IF
    res=.FALSE.
  END FUNCTION cp_sm_b_i_next
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_write [1.0] *
!!
!!   NAME
!!     cp_sm_write
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_write(matrix, unit_nr, long_description, local,&
!!         error)
!!       Type(real_matrix_type), Pointer:: matrix
!!       Integer, Intent (IN):: unit_nr
!!       Logical, Intent (IN), Optional:: long_description, local
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_write
!!
!!   FUNCTION
!!     writes a matrix to the given output unit
!!
!!   NOTES
!!     to do
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
SUBROUTINE cp_sm_write(matrix,unit_nr,long_description,local,error)
    TYPE(real_matrix_type), POINTER          :: matrix
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: long_description, local
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sm_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    INTEGER                                  :: iblock_col, iblock_row
    LOGICAL                                  :: failure, my_local, &
                                                my_long_description
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix_block
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(real_block_node_type), POINTER      :: block_node

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  my_local=.FALSE.
  my_long_description=.FALSE.
  IF (PRESENT(local)) my_local=local
  IF (PRESENT(long_description)) my_long_description=long_description
  IF (.not.my_local) my_local=(logger%para_env%mepos==logger%para_env%source)
  
  IF (ASSOCIATED(matrix)) THEN
     IF (my_local) THEN
        CALL get_matrix_info(matrix, matrix_name=matrix_name, &
             matrix_symmetry=matrix_symmetry)
        WRITE(unit=unit_nr,fmt="(' <real_matrix_type>{ name=',a)")&
             TRIM(matrix_name)
        WRITE(unit=unit_nr,fmt="(' symmetry=',a)")&
             TRIM(matrix_symmetry)
        DO iblock_row=1,matrix%nblock_row
           block_node => first_block_node(matrix,iblock_row)
           DO WHILE (ASSOCIATED(block_node))

              CALL get_block_node(block_node=block_node,&
                BLOCK=matrix_block,block_col=iblock_col)
              IF (ASSOCIATED(matrix_block)) THEN
                 IF (my_long_description) THEN
                    WRITE (unit=unit_nr,fmt="(' block(',i6,',',i6,')=')")&
                         iblock_row,iblock_col
                    CALL cp_2d_r_write(matrix_block,unit_nr=unit_nr,error=error)
                 ELSE
                    WRITE (unit=unit_nr,&
                         fmt="(' block(',i6,',',i6,')=array(',i6,',',i6,')')")&
                         iblock_row,iblock_col,&
                         SIZE(matrix_block,1),SIZE(matrix_block,2)
                 END IF
              ELSE
                 WRITE (unit=unit_nr,fmt="(' block(',i6,',',i6,')=*null*')")&
                      iblock_row,iblock_col
              END IF
              
              block_node => next_block_node(block_node)
           END DO
        END DO
        WRITE(unit=unit_nr,fmt="(' }')")
        CALL m_flush(unit_nr)
     END IF
  ELSE
     IF (my_local) THEN
        WRITE(unit=unit_nr,fmt="(' <real_matrix_type *null*>')")
        CALL m_flush(unit_nr)
     END IF
  END IF
END SUBROUTINE cp_sm_write
!***************************************************************************

!!****f* cp_matrix_utils/sm_struct_add_sm_local_blocks [1.0] *
!!
!!   NAME
!!     sm_struct_add_sm_local_blocks
!!
!!   SYNOPSIS
!!     Subroutine sm_struct_add_sm_local_blocks(sm_struct, matrix, error)
!!       Type(cp_sm_struct_type), Pointer:: sm_struct
!!       Type(real_matrix_type), Pointer:: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine sm_struct_add_sm_local_blocks
!!
!!   FUNCTION
!!     Adds the blocks of the given matrix to the local blocks of the
!!     given matrix structure
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - sm_struct: the structure to update
!!     - matrix: the matrix from where to get the blocks
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
SUBROUTINE sm_struct_add_sm_local_blocks(sm_struct,matrix,error)
    TYPE(cp_sm_struct_type), POINTER         :: sm_struct
    TYPE(real_matrix_type), POINTER          :: matrix
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'sm_struct_add_sm_local_blocks', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: block_col, block_row, handle
    LOGICAL                                  :: failure
    TYPE(cp_matrix_block_iterator)           :: m_iter

  failure=.FALSE.
  
  CALL timeset(routineN,"I","",handle)
  CPPrecondition(ASSOCIATED(sm_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(sm_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(matrix),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL cp_sm_b_i_init(m_iter, matrix=matrix, error=error)
     
     !
     ! appears to be quadratic scaling 
     !  
     DO WHILE( cp_sm_b_i_next(m_iter, block_row=block_row, &
          block_col=block_col,error=error))
        CALL sm_struct_add_local_block(sm_struct, block_row=block_row,&
             block_col=block_col, error=error)
     END DO
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE sm_struct_add_sm_local_blocks
!***************************************************************************

END MODULE cp_matrix_utils
