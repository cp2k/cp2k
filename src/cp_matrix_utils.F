!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
module cp_matrix_utils
  use cp_log_handling, only : cp_to_string, cp_note_level, cp_warning_level,&
       cp_failure_level, cp_fatal_level
  use cp_error_handling, only : cp_assert, cp_error_type, cp_init,&
       cp_dealloc_ref, cp_error_propagate_error, cp_debug, cp_assertion_failed,&
       cp_error_message
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use sparse_matrix_types
  use message_passing, only: mp_sum
  !use qs_blacs, only:
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_matrix_utils'
! types
  public :: cp_matrix_block_iterator

! core methods
  public :: cp_init, cp_dealloc_ref, cp_get, cp_next, cp_valid

  public :: cp_trace

  public :: cp_p_norm_matrix, cp_sup_norm_matrix,&
       cp_is_block_diagonal, cp_init_diagonal_blocks,&
       cp_d_struct_create, cp_set_to
  
  ! underliying methods (try to use the generic interface)
  public :: cp_sm_diag_trace, cp_sm_set_to, cp_matrix_b_i_init,&
       cp_matrix_b_i_dealloc_ref, cp_matrix_b_i_get, &
       cp_matrix_b_i_next_block, cp_matrix_b_i_valid, &
       cp_sm_is_block_diagonal, cp_sm_init_diagonal_blocks

  ! return the trace of the given matrices
  interface cp_trace
     module procedure cp_sm_diag_trace
  end interface

!! set all the values to a given value
  interface cp_set_to
     module procedure cp_sm_set_to
  end interface

! = core methods =

!! initializes the object
  interface cp_init
     module procedure cp_matrix_b_i_init
  end interface

!! deallocates the memory allocated by the object
  interface cp_dealloc_ref
     module procedure cp_matrix_b_i_dealloc_ref
  end interface

!! gets various attributes from the object
  interface cp_get
     module procedure cp_matrix_b_i_get
  end interface

!! goes to the next element, returns true if valid
  interface cp_next
     module procedure cp_matrix_b_i_next_block
  end interface

!! return true if it is valid (minimal testing)
  interface cp_valid
     module procedure cp_matrix_b_i_valid
  end interface

!! returns true if the matrix is block diagonal
  interface cp_is_block_diagonal
     module procedure cp_sm_is_block_diagonal
  end interface

!! initializes the diagonal matrix blocks
  interface cp_init_diagonal_blocks
     module procedure cp_sm_init_diagonal_blocks
  end interface

!!***
!****************************************************************************

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
  type cp_matrix_block_iterator
     private
     logical :: initialized
     type(real_matrix_type), pointer :: matrix
     integer :: block_row
     type(real_block_node_type), pointer :: block_node
  end type cp_matrix_block_iterator
!!***
!****************************************************************************

contains
! ============== real full matrix utilities ============

!!****f* cp_matrix_utils/cp_p_norm_matrix [1.0] *
!!
!!   NAME
!!     cp_p_norm_matrix
!!
!!   SYNOPSIS
!!     Function cp_p_norm_matrix(matrix, p) Result(res)
!!       Real(Kind=wp):: res
!!       Real(Kind=wp), Intent (IN), Dimension(:,:):: matrix
!!       Integer, Optional, Intent (IN):: p
!!     End Function cp_p_norm_matrix
!!
!!   FUNCTION
!!     calculates the p-te power of the lp norm on the given matrix:
!!     ||a_{i j}||_{p}^p=\sum_{i,j} a_{i j}^p
!!     The frobenious norm is the one with p=2
!!
!!   NOTES
!!     could be much more efficient
!!
!!   INPUTS
!!     matrix: the matrix whose norm is wanted
!!     p: exponent whose power is used to calculate the norm,
!!        defaults to 2
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  !pure
  function cp_p_norm_matrix(matrix, p) result(res)
    real(kind=wp) :: res
    real(kind=wp), intent(in), dimension(:,:) :: matrix
    integer,optional, intent(in) :: p
!type(cp_error_type), optional, intent(inout) :: error ! removed to make it pure

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_froebnius_norm2'
    integer :: i,j, pVal
    failure=.FALSE.

    res=0
    if (.NOT. failure) then
       pVal=2
       if (.not. present(p)) pVal=p
       do i=1,size(matrix,1)
          do j=1,size(matrix,2)
             res=res+matrix(i,j)**pVal
          end do
       end do
    end if
  end function cp_p_norm_matrix
!***************************************************************************

!!****f* cp_matrix_utils/cp_sup_norm_matrix [1.0] *
!!
!!   NAME
!!     cp_sup_norm_matrix
!!
!!   SYNOPSIS
!!     Function cp_sup_norm_matrix(matrix) Result(res)
!!       Real(Kind=wp):: res
!!       Real(Kind=wp), Dimension(:,:), Intent (IN):: matrix
!!     End Function cp_sup_norm_matrix
!!
!!   FUNCTION
!!     returns the supremum (maximum) norm of the matrix
!!     ||a_{i j}||_{sup}=max_{i,j}(|a_{i j}|)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix whose norm you want to calculate
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  !pure 
  function cp_sup_norm_matrix(matrix) result(res)
    real(kind=wp) :: res
    real(kind=wp), dimension(:,:), intent(in) :: matrix
!type(cp_error_type), optional, intent(inout) :: error ! removed to make it pure

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_sup_norm_matrix'
    failure=.FALSE.

    if (.NOT. failure) then
       res=maxval(abs(matrix))
    else
       res=-1.0_wp
    end if
  end function cp_sup_norm_matrix
!***************************************************************************

! ================== cp_matrix_block_iterator ================

!!****f* cp_matrix_utils/cp_matrix_b_i_init [1.0] *
!!
!!   NAME
!!     cp_matrix_b_i_init
!!
!!   SYNOPSIS
!!     Subroutine cp_matrix_b_i_init(block_iterator, matrix, block_row,&
!!         block_node, error)
!!       Type(cp_matrix_block_iterator), Intent (OUT):: block_iterator
!!       Type(real_matrix_type), Target, Intent (IN):: matrix
!!       Integer, Intent (IN), Optional:: block_row
!!       Type(real_block_node_type), Target, Optional:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_matrix_b_i_init
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
  subroutine cp_matrix_b_i_init(block_iterator, matrix, block_row,&
       block_node, error)
    type(cp_matrix_block_iterator), intent(out) :: block_iterator
    type(real_matrix_type), target, intent(in) :: matrix
! well subsequent operations on the iterator could change tha data in the matrix
    integer, intent(in), optional :: block_row
    type(real_block_node_type), target, optional :: block_node
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_matrix_b_i_init'
    failure=.false.

    block_iterator%matrix => matrix
    block_iterator%block_row = 1
    nullify(block_iterator%block_node)
    block_iterator%initialized=.true.
    if (present(block_node)) &
         block_iterator%block_node => block_node
    if (present(block_row)) &
         block_iterator%block_row = block_row
  end subroutine cp_matrix_b_i_init
!***************************************************************************

!!****f* cp_matrix_utils/cp_matrix_b_i_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_matrix_b_i_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_matrix_b_i_dealloc_ref(block_iterator, error)
!!       Type(cp_matrix_block_iterator), Intent (INOUT):: block_iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_matrix_b_i_dealloc_ref
!!
!!   FUNCTION
!!     deallocates the memory allocated by the matrix block iterator
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the block iterator you want to free
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
  subroutine cp_matrix_b_i_dealloc_ref(block_iterator, error)
    type(cp_matrix_block_iterator), intent(inout) :: block_iterator
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_matrix_b_i_dealloc_ref'
    failure=.false.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       nullify(block_iterator%matrix, block_iterator%block_node)
       block_iterator%block_row=0
    end if
  end subroutine cp_matrix_b_i_dealloc_ref
!***************************************************************************

!!****f* cp_matrix_utils/cp_matrix_b_i_get [1.0] *
!!
!!   NAME
!!     cp_matrix_b_i_get
!!
!!   SYNOPSIS
!!     Subroutine cp_matrix_b_i_get(block_iterator, block_row, block_col,&
!!         block_val, finished, matrix, block_node, error)
!!       Type(cp_matrix_block_iterator), Intent (IN):: block_iterator
!!       Integer, Optional, Intent (OUT):: block_row, block_col
!!       Real(Kind=wp), Dimension(:,:), Pointer, Optional:: block_val
!!       Logical, Optional, Intent (OUT):: finished
!!       Type(real_matrix_type), Optional, Pointer:: matrix
!!       Type(real_block_node_type), Optional, Pointer:: block_node
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_matrix_b_i_get
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
  subroutine cp_matrix_b_i_get(block_iterator, block_row, block_col,&
       block_val, finished, matrix, block_node, error)
    type(cp_matrix_block_iterator), intent(in) :: block_iterator
    integer, optional, intent(out) :: block_row, block_col
    real(kind=wp), dimension(:,:), pointer, optional :: block_val
    logical, optional, intent(out) :: finished
    type(real_matrix_type), optional, pointer :: matrix
    type(real_block_node_type), optional, pointer :: block_node
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_matrix_b_i_get'
    failure=.false.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(block_iterator%matrix),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(block_row)) &
            block_row=block_iterator%block_row
       if (present(block_col)) then
          if (associated(block_iterator%block_node)) then
             call get_block_node(block_iterator%block_node, block_col=block_col)
          else
             CPAssert(.false.,cp_warning_level,routineP,error,failure)
             block_col=-1
          end if
       end if
       if (present(block_val)) then
          CPPrecondition(block_iterator%block_row>0,cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             if (associated(block_iterator%block_node)) then
                call get_block_node(block_iterator%block_node,block=block_val)
             else
                nullify(block_val)
                CPErrorMessage(cp_failure_level,routineP,"invalid block",error)
             end if
          end if
       end if
       if (present(finished)) &
            finished = block_iterator%block_row <= 0
       if (present(matrix)) &
            matrix => block_iterator%matrix
       if (present(block_node)) &
            block_node => block_iterator%block_node
    end if
  end subroutine cp_matrix_b_i_get
!***************************************************************************

!!****f* cp_matrix_utils/cp_matrix_b_i_valid [1.0] *
!!
!!   NAME
!!     cp_matrix_b_i_valid
!!
!!   SYNOPSIS
!!     Function cp_matrix_b_i_valid(block_iterator, error) Result(res)
!!       Logical:: res
!!       Type(cp_matrix_block_iterator), Intent (IN):: block_iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_matrix_b_i_valid
!!
!!   FUNCTION
!!     returns true if the iterator is valid (not past end or invalid).
!!     a fast minimal testing is done, and the errors are written as warnings
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the iterator you get the info from
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
  function cp_matrix_b_i_valid(block_iterator,error) result(res)
    logical ::res
    type(cp_matrix_block_iterator), intent(in) :: block_iterator
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_matrix_b_i_valid'
    failure=.false.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (.not.(associated(block_iterator%matrix).and. block_iterator%block_row>0 ))then
          CPAssert(block_iterator%block_row > 0,cp_warning_level,routineP,error,failure)
          CPAssert(associated(block_iterator%matrix),cp_warning_level,routineP,error,failure)
       end if
    end if
    res=.not.failure
  end function cp_matrix_b_i_valid
!***************************************************************************

!!****f* cp_matrix_utils/cp_matrix_b_i_next_block [1.0] *
!!
!!   NAME
!!     cp_matrix_b_i_next_block
!!
!!   SYNOPSIS
!!     Function cp_matrix_b_i_next_block(block_iterator, error) Result(res)
!!       Logical:: res
!!       Type(cp_matrix_block_iterator), Intent (INOUT):: block_iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_matrix_b_i_next_block
!!
!!   FUNCTION
!!     Moves the iterator to the next matrix block.
!!     Returns true if the iterator is valid (not gone past end)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_iterator: the block iterator that you want to move
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
  function cp_matrix_b_i_next_block(block_iterator, error) result(res)
    logical :: res
    type(cp_matrix_block_iterator), intent(inout) :: block_iterator
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_matrix_b_i_next_block'
    integer :: nblock_row,i
    failure=.false.

    CPPrecondition(block_iterator%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(block_iterator%matrix),cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (block_iterator%block_row > 0) then
          if (associated(block_iterator%block_node)) then
             block_iterator%block_node => next_block_node(block_iterator%block_node)
             if (.not.associated(block_iterator%block_node)) then
                block_iterator%block_row= block_iterator%block_row+1
             end if
          end if
          if (.not.associated(block_iterator%block_node)) then
             call get_matrix_info(block_iterator%matrix,nblock_row=nblock_row)
             do i=block_iterator%block_row,nblock_row
                block_iterator%block_node => &
                     first_block_node(block_iterator%matrix,&
                     block_iterator%block_row)
                if (associated(block_iterator%block_node)) then
                   block_iterator%block_row=i
                   exit
                end if
             end do
          end if
          if (.not.associated(block_iterator%block_node)) then
             block_iterator%block_row=-1
          else
             res=.true.
             return
          end if
       end if
       CPPostcondition(.not.associated(block_iterator%block_node),cp_warning_level,routineP,error,failure)
       if (.not. block_iterator%block_row <= 0) then
          CPErrorMessage(cp_warning_level,routineP,"inconsistent state in iterator",error)
          block_iterator%block_row=-2
       end if
    else
       block_iterator%block_row=-2
       nullify(block_iterator%block_node)
    end if
    res=.false.
  end function cp_matrix_b_i_next_block
!***************************************************************************

! ============ real block (sparse) matrix utils ===========

!!****f* cp_matrix_utils/cp_sm_diag_trace [1.0] *
!!
!!   NAME
!!     cp_sm_diag_trace
!!
!!   SYNOPSIS
!!     Function cp_sm_diag_trace(diagonal, matrix, global_env, error)&
!!         Result(res)
!!       Real(Kind=wp):: res
!!       Real(Kind=wp), Dimension(:), Intent (IN):: diagonal
!!       Type(real_matrix_type), Intent (IN), Target:: matrix
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sm_diag_trace
!!
!!   FUNCTION
!!     calculates the trace of the given matriz with the diagonal matrix
!!     given through his diagonal elements. Missing giagonal elements are
!!     assumed to be 0.
!!     Tr(diag(diagonal) matrix)= sum(diagonal(i)*matrix(i,i))
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
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
function cp_sm_diag_trace(diagonal, matrix, global_env, error) result(res)
  real(kind=wp) :: res
  real(kind=wp), dimension(:), intent(in) :: diagonal
  type(real_matrix_type), intent(in), target :: matrix
  type(global_environment_type), intent(in),target :: global_env
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle, i, j
  character(len=*), parameter :: routineN='cp_sm_diag_trace',&
       routineP=moduleN//':'//routineN
  type(cp_matrix_block_iterator) :: iterator
  real(kind=wp), dimension(:,:), pointer :: block_val
  integer, dimension(:), pointer :: first_row,first_col
  type(real_matrix_type), pointer :: matrix_ptr
  failure=.false.
  matrix_ptr => matrix
  
  call timeset(routineN//','//moduleN,'I',"",handle)
  call cp_init(iterator, matrix=matrix)
  call get_matrix_info(matrix_ptr, first_row=first_row, first_col=first_col)
  res=0.0_wp
  do
     if (.not.cp_next(iterator,error=error)) exit
     call cp_get(iterator, block_row=i,block_col=j, block_val=block_val,&
          error=error)
     if (i==j) then
        CPPrecondition(size(block_val,1)==size(block_val,2),cp_failure_level,routineP,error,failure)
        CPPrecondition(first_col(i)==first_row(i),cp_failure_level,routineP,error,failure)
        do i=1,min(size(block_val,1),size(block_val,2))
           res=res+block_val(i,i)+diagonal(first_col(i)+i)
        end do
     end if
  end do
  call mp_sum(res,global_env%group)
  call timestop(0.0_wp,handle)
end function cp_sm_diag_trace
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_set_to [1.0] *
!!
!!   NAME
!!     cp_sm_set_to
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_set_to(matrix, value, error)
!!       Type(real_matrix_type), Intent (INOUT), Target:: matrix
!!       Real(Kind=wp), Optional, Intent (IN):: value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_set_to
!!
!!   FUNCTION
!!     initializes the allocated blocks of the real matrix with the given
!!     value (defaults to 0)
!!
!!   NOTES
!!     the iteration looks ugly because this function has much error
!!     handling to test the iterator, it should be:
!!     do; if (.not.cp_next(iterator)) exit
!!       call cp_get(..);..
!!     end do
!!
!!   INPUTS
!!     matrix: the matrix to initialize
!!     value: the new value
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
  subroutine cp_sm_set_to(matrix, value, error)
    type(real_matrix_type), intent(inout), target :: matrix
    real(kind=wp), optional, intent(in) :: value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_matrix_utils:cp_sm_set_to'
    type(cp_matrix_block_iterator) :: iterator
    real(kind=wp), dimension(:,:), pointer :: block_val
    real(kind=wp) :: val
    type(cp_error_type) :: iError
    failure=.false.
    call cp_init(iError,template_error=error)
    val=0.0_wp

    if (.not. failure) then
       call cp_init(iterator,matrix,error=iError)
       call cp_error_propagate_error(iError,routineP,&
            CPSourceFileRef,&
            error,failure)
    end if
    if (.not.failure) then
       if (present(value)) val=value
       do
          if (cp_next(iterator,error=iError)) then
             call cp_get(iterator,block_val=block_val,error=iError)
             call cp_error_propagate_error(iError,routineP,&
                  CPSourceFileRef,&
                  error,failure)
             CPPrecondition(associated(block_val),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                block_val=val
             else
                exit
             end if
          else
             call cp_error_propagate_error(iError,routineP,&
                  CPSourceFileRef,&
                  error,failure)
             exit
          end if
       end do
       call cp_dealloc_ref(iterator,error=error)
    end if
    call cp_dealloc_ref(iError,error=error)
  end subroutine cp_sm_set_to
!***************************************************************************

!!****f* cp_matrix_utils/cp_build_d_struct [1.0] *
!!
!!   NAME
!!     cp_build_d_struct
!!
!!   FUNCTION
!!     builds the distribution structure of the diagonal blocks
!!     of the given matrix (the array is allocated and you are responsible
!!     for it)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     d_struct: a pointer that will point to an array with the
!!            distribution structure of the diagonal elements
!!     matrix: the matrix you get the structure from (if not given
!!            distributed in blocks of equal size between the processors)      
!!     n_blocks: the number of diagonal blocks, if no matrix is given 
!!            this must be given
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
subroutine cp_d_struct_create(d_struct, matrix, n_blocks, global_env, error)
  integer, dimension(:), pointer :: d_struct
  type(global_environment_type), intent(in), target :: global_env
  type(real_matrix_type), intent(in), target, optional :: matrix
  integer, optional, intent(in) :: n_blocks
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: handle, i, j, nrows, ncols, stat, rest, nmin
  character(len=*), parameter :: routineN='cp_build_d_struct',&
       routineP=moduleN//':'//routineN
  type(cp_matrix_block_iterator) :: iterator
  type(real_matrix_type), pointer :: matrix_ptr
  failure=.false.
  
  !call timeset(routineN//','//moduleN,'I',"",handle)
  if (present(matrix)) then
     matrix_ptr => matrix
     call get_matrix_info(matrix_ptr, nblock_row=nrows, nblock_col=ncols)
     CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
     ncols=min(nrows,ncols)
     if (present(n_blocks)) then
        CPPrecondition(ncols==n_blocks,cp_warning_level,routineP,error,failure)
     end if
     if (.not.failure) then
        ! assume every diagonal block is on one (and just one) proc
        allocate(d_struct(ncols),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        d_struct=0
        call cp_init(iterator, matrix=matrix_ptr, error=error)
        do
           if (.not.cp_next(iterator)) exit
           call cp_get(iterator,block_row=i,block_col=j,error=error)
           if (i==j) then
              d_struct(i)=global_env%mepos
           end if
        end do
        call cp_dealloc_ref(iterator,error=error)
        call mp_sum(d_struct,global_env%group)
     end if
  else
     CPPrecondition(present(n_blocks),cp_failure_level,routineP,error,failure)
     if (.not.failure) then
        allocate(d_struct(n_blocks),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     end if
     if (.not.failure) then
        nmin= n_blocks/global_env%num_pe
        rest=mod(n_blocks,global_env%num_pe)
        do j=0,(rest-1)
           do i=1,nmin+1
              d_struct(j*(nmin+1)+i)=j
           end do
        end do
        do j=rest,(global_env%num_pe-1)
           do i=1,nmin+1
              d_struct(rest+j*nmin+i)=j
           end do
        end do
     end if
  end if
  !call timestop(0.0_wp,handle)
end subroutine cp_d_struct_create
!***************************************************************************


!!****f* cp_matrix_utils/cp_sm_init_diagonal_blocks [1.0] *
!!
!!   NAME
!!     cp_sm_init_diagonal_blocks
!!
!!   SYNOPSIS
!!     Subroutine cp_sm_init_diagonal_blocks(matrix, global_env, d_struct,&
!!         default_value, diagonal_value, error)
!!       Type(real_matrix_type), Intent (INOUT), Target:: matrix
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Integer, Optional, Dimension(:), Intent (IN), Target:: d_struct
!!       Real(Kind=wp), Intent (IN), Optional:: default_value,&
!!         diagonal_value
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sm_init_diagonal_blocks
!!
!!   FUNCTION
!!     initialize the block diagonal blocks of the given sparse matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix to be initialized
!!     d_struct: the distribution structure of the diagonal blocks
!!     default_value: the default value of the elements (if not given
!!             no initialization is performed)
!!     diagonal_value: the value of the elements on the diagonal
!!             (defaults to default_value)
!!     error: variable to control error logging, stopping,...
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sm_init_diagonal_blocks(matrix, global_env, d_struct, &
       default_value, diagonal_value, error)
    type(real_matrix_type), intent(inout), target :: matrix
    type(global_environment_type), intent(in), target :: global_env
    integer, optional, dimension(:), intent(in), target :: d_struct
    real(kind=wp), intent(in), optional :: default_value, diagonal_value
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, i, j, nrows, ncols, rest, nblocks, stat
    integer, pointer, dimension(:) :: my_d_struct
    type(real_block_node_type), pointer :: m_block
    real(kind=wp), dimension(:,:), pointer :: block_val
    character(len=*), parameter :: routineN='cp_sm_init_diagonal_blocks',&
         routineP=moduleN//':'//routineN
    type(real_matrix_type), pointer :: matrix_ptr
    matrix_ptr => matrix
    failure=.false.

    CPPrecondition(associated(matrix_ptr),cp_failure_level,routineP,error,failure)
    call timeset(routineN//','//moduleN,'I',"",handle)
    if (.not.failure) then
       call get_matrix_info(matrix_ptr, nblock_row=nrows, nblock_col=ncols)
       CPPreconditionNoFail(nrows==ncols,cp_warning_level,routineP,error)
       ncols=min(ncols,nrows)
       if (present(d_struct)) then
          CPPrecondition(ncols==size(d_struct),cp_failure_level,routineP,error,failure)
          my_d_struct => d_struct
       else
          call cp_d_struct_create(my_d_struct, n_blocks=ncols, &
               global_env=global_env, error=error)
       end if
       CPPostcondition(associated(my_d_struct),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          do i=1,size(my_d_struct)
             if (global_env%mepos==my_d_struct(i)) then
                m_block => find_block_node(matrix_ptr,i,i)
                if (.not.associated(m_block)) then
                   call add_block_node(matrix_ptr,i,i,block_val)
                else
                   call get_block_node(m_block,block=block_val)
                end if
                if (present(default_value)) block_val=default_value
                if (present(diagonal_value)) then
                   forall (j=1:min(size(block_val,1),size(block_val,2)))
                      block_val(j,j)=diagonal_value
                   end forall
                end if
             else
                m_block => find_block_node(matrix_ptr,i,i)
                if (associated(m_block)) then
                   CPErrorMessage(cp_failure_level,routineP,"diag block on wrong processor",error)
                   ! call remove_real_block_node(matrix_ptr,i,i) ! to do
                end if
             end if
          end do
       end if
       if (cp_debug.and.debug_this_module) then ! test if diag
          call cp_assert(cp_is_block_diagonal(matrix,global_env),&
               cp_failure_level,&
               cp_assertion_failed,routineP, "matrix is not diagonal in "//&
               CPSourceFileRef,&
               error,failure)
       end if
    end if
  call timestop(0.0_wp,handle)
end subroutine cp_sm_init_diagonal_blocks
!***************************************************************************

!!****f* cp_matrix_utils/cp_sm_is_block_diagonal [1.0] *
!!
!!   NAME
!!     cp_sm_is_block_diagonal
!!
!!   SYNOPSIS
!!     Function cp_sm_is_block_diagonal(matrix, global_env, error)&
!!         Result(res)
!!       Logical:: res
!!       Type(real_matrix_type), Intent (INOUT), Target:: matrix
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sm_is_block_diagonal
!!
!!   FUNCTION
!!     returns true if the given matrix is block diagonal
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     matrix: the matrix to check
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 2.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
function cp_sm_is_block_diagonal(matrix, global_env, error) result(res)
  logical :: res
  type(real_matrix_type), intent(inout), target :: matrix
  type(global_environment_type), intent(in), target :: global_env
  type(cp_error_type), optional, intent(inout) :: error

  logical :: failure
  integer, dimension(:), allocatable :: myDiag
  integer :: handle, i,j, nrows,ncols, stat
  character(len=*), parameter :: routineN='cp_sm_is_block_diagonal',&
       routineP=moduleN//':'//routineN
  type(cp_matrix_block_iterator) :: iterator
  type(real_matrix_type), pointer :: matrix_ptr
  failure=.false.
  matrix_ptr => matrix

  !call timeset(routineN//','//moduleN,'I',"",handle)
  CPPrecondition(associated(matrix_ptr),cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call get_matrix_info(matrix_ptr, nrow=nrows, ncol=ncols)
     CPAssert(ncols==nrows,cp_warning_level,routineP,error,failure)
     ncols=min(nrows,ncols)
     allocate(myDiag(ncols),stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     myDiag=0
     call cp_init(iterator, matrix=matrix, error=error)
     do
        if (.not.cp_next(iterator)) exit
        call cp_get(iterator,block_row=i,block_col=j,error=error)
        if (i==j) then
           myDiag(i)=myDiag(i)+1
        else
           call cp_assert(.false.,cp_warning_level,cp_assertion_failed,&
                routineP,"there is an off diagonal block! in "//&
                CPSourceFileRef,error)
           myDiag(i)=myDiag(i)+100
           myDiag(j)=myDiag(j)+100
        end if
     end do
     call cp_dealloc_ref(iterator,error=error)
     call mp_sum(myDiag,global_env%group)
     ! there are off diagonal blocks
     CPAssert(all(myDiag<100),cp_warning_level,routineP,error,failure)
     do i=1,size(myDiag)
        myDiag(i)= mod(myDiag(i),100)
     end do
     CPAssert(all(myDiag==1),cp_failure_level,routineP,error,failure)
     deallocate(myDiag, stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)     
  end if
  res=.not.failure
  !call timestop(0.0_wp,handle)
end function cp_sm_is_block_diagonal
!***************************************************************************

end module cp_matrix_utils
