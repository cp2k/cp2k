!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sparse_matrix [1.0] *
!!
!!   NAME
!!     cp_sparse_matrix
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
module cp_sparse_matrix
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use string_utilities, only: uppercase
  use timings, only: timeset, timestop
  use cp_b_matrix_structure
  use sparse_matrix_types
  use message_passing, only: mp_sync
  use cp_matrix_utils, only: cp_sm_alloc_nonsparse_blocks
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_sparse_matrix'

  public :: cp_sparse_matrix_type, cp_sp_block_handle_type,&
       cp_sparse_matrix_p_type
  public :: cp_create, cp_dealloc, cp_retain, cp_release, &
       cp_consolidate_matrix, cp_get_local_block,&
       cp_release_block, cp_set_local_block, cp_remove_block,&
       cp_get_matrix
  public :: cp_sp_create, cp_sp_dealloc, cp_sp_retain, cp_sp_release,&
       cp_sp_consolidate_matrix, cp_sp_get_local_block,&
       cp_sp_release_block, cp_sp_set_local_block, cp_sp_remove_block,&
       cp_sp_array_dealloc, cp_sp_get_matrix, cp_sp_array_release,&
       cp_sp_alloc_nonsparse_blocks

  interface cp_create
     module procedure cp_sp_create
  end interface
  interface cp_dealloc
     module procedure cp_sp_dealloc
  end interface
  interface cp_retain
     module procedure cp_sp_retain
  end interface
  interface cp_release
     module procedure cp_sp_release, cp_sp_array_release
  end interface
  interface cp_get_matrix
     module procedure cp_sp_get_matrix
  end interface
  interface cp_consolidate_matrix
     module procedure cp_sp_consolidate_matrix
  end interface
!!FM  interface cp_flush_cache
!!FM     module procedure cp_sp_flush_cache
!!FM  end interface
  interface cp_get_local_block
     module procedure cp_sp_get_local_block
  end interface
  interface cp_release_block
     module procedure cp_sp_release_block
  end interface
  interface cp_set_local_block
     module procedure cp_sp_set_local_block
  end interface
  interface cp_remove_block
     module procedure cp_sp_remove_block
  end interface
  interface cp_dealloc
     module procedure cp_sp_array_dealloc
  end interface

  !***
  !****************************************************************************

!!****s* cp_sparse_matrix/cp_sparse_matrix_type [1.0] *
!!
!!   NAME
!!     cp_sparse_matrix_type
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
!!     - allowed should_dealloc_matrix: true if the matrix should be 
!!       deallocated when the structure is deallocated
!!     - matrix_consolidated: if the matrix has been consolidated (get works)
!!     - symmetric: if the matrix is symmetric
!!     - matrix_struct: the matrix structure of this matrix
!!     - ref_count: reference counter, to know when this shared matrix
!!       can be safely deallocated
!!     - matrix: the underlying real matrix (make more efficient?
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
  type cp_sparse_matrix_type
     logical :: initialized
     logical :: frozen
     logical :: should_dealloc_matrix
     logical :: cache_consolidated
     type(cp_b_matrix_struct_type), pointer :: matrix_struct
     integer :: ref_count, id_nr
     type(real_matrix_type), pointer :: matrix
     type(global_environment_type), pointer :: global_env
  end type cp_sparse_matrix_type
!!***
  !****************************************************************************

!!****s* cp_sparse_matrix/cp_sparse_matrix_p_type [1.0] *
!!
!!   NAME
!!     cp_sparse_matrix_p_type
!!
!!   FUNCTION
!!     just to have arrays of pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - matrix: the pointerto the matrix
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  type cp_sparse_matrix_p_type
     type(cp_sparse_matrix_type), pointer :: matrix
  end type cp_sparse_matrix_p_type
!!***
  !****************************************************************************

!!****s* cp_sparse_matrix_type/cp_sp_block_handle_type [1.0] *
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
!!     matrix: the matrix that 
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
  type cp_sp_block_handle_type
     type(cp_sparse_matrix_type), pointer :: matrix
     integer :: block_row, block_col
     character :: trans,trans_to_do
     real(kind=wp), dimension(:,:), pointer :: block_val
     logical :: should_dealloc_data, initialized
  end type cp_sp_block_handle_type
!!***
  !****************************************************************************

contains

!!****f* cp_sparse_matrix/cp_sp_create [1.0] *
!!
!!   NAME
!!     cp_sp_create
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_create(matrix, global_env, matrix_structure,&
!!         u_matrix, should_dealloc_matrix, error)
!!       Type(cp_sparse_matrix_type), Pointer:: matrix
!!       Type(global_environment_type), Intent (IN), Target:: global_env
!!       Type(cp_b_matrix_struct_type), Intent (IN), Target::&
!!         matrix_structure
!!       Type(real_matrix_type), Optional, Intent (IN), Target:: u_matrix
!!       Logical, Intent (IN), Optional:: should_dealloc_matrix
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
!!     - global_env: the mpi environement that this matrix will use (SHARED)
!!     - matrix_structure: the structure of this matrix
!!     - u_matrix: the underlying sparse matrix SHARED
!!     - should_dealloc_matrix:if the underlying matrix should be deallocated
!!       (defaults to true)
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
  subroutine cp_sp_create(matrix, global_env, matrix_structure, u_matrix, &
       should_dealloc_matrix, error)
    type(cp_sparse_matrix_type), pointer :: matrix
    type(global_environment_type), intent(in), target :: global_env
    type(cp_b_matrix_struct_type), intent(in), target :: matrix_structure
    type(real_matrix_type), optional, intent(in), target :: u_matrix
    logical, intent(in), optional :: should_dealloc_matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_create',&
         routineP=moduleN//':'//routineN
    character(len=10) :: m_symm
    type(cp_b_dims_type), pointer :: b_dims
    integer, save :: last_id_nr=0
    integer :: stat, block_cols, block_rows
    failure=.false.

    allocate(matrix,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       last_id_nr=last_id_nr+1
       matrix%id_nr=last_id_nr
       matrix%frozen=.false.
       matrix%should_dealloc_matrix=.true.
       if (present(should_dealloc_matrix)) &
            matrix%should_dealloc_matrix=should_dealloc_matrix
       matrix%cache_consolidated=.false.
       matrix%matrix_struct => matrix_structure
       call cp_retain(matrix%matrix_struct,error=error)
       matrix%ref_count=1
       matrix%global_env => global_env
       nullify(matrix%matrix)
       matrix%initialized=.true.
       if (present(u_matrix)) then
          matrix%matrix => u_matrix
       else
          if (matrix%matrix_struct%symmetric) then
             m_symm='symmetric'
          else
             m_symm='none'
          end if
          b_dims => matrix%matrix_struct%b_dims
          block_rows=size(b_dims%block_begins_at_row)-1
          block_cols=size(b_dims%block_begins_at_col)-1
          call allocate_matrix(matrix%matrix,&
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
       end if
    end if
  end subroutine cp_sp_create
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_dealloc [1.0] *
!!
!!   NAME
!!     cp_sp_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_dealloc(matrix, error)
!!       Type(cp_sparse_matrix_type), Pointer:: matrix
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
  subroutine cp_sp_dealloc(matrix, error)
    type(cp_sparse_matrix_type), pointer :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_dealloc',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    if (associated(matrix)) then
       CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
       CPPreconditionNoFail(matrix%ref_count==0,cp_warning_level,routineP,error)
       if (.not. failure) then
          matrix%initialized=.false.
          call cp_release(matrix%matrix_struct, error=error)
          nullify(matrix%matrix_struct)
          nullify(matrix%global_env)
          if (matrix%should_dealloc_matrix) then
             if (associated(matrix%matrix)) then
                call deallocate_matrix(matrix%matrix)
             end if
          else
             nullify(matrix%matrix)
          end if
          deallocate(matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       end if
    end if
  end subroutine cp_sp_dealloc
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_retain [1.0] *
!!
!!   NAME
!!     cp_sp_retain
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_retain(matrix, error)
!!       Type(cp_sparse_matrix_type), Intent (INOUT), Target:: matrix
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
  subroutine cp_sp_retain(matrix, error)
    type(cp_sparse_matrix_type), intent(inout), target :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_retain',&
         routineP=moduleN//':'//routineN
    type(cp_sparse_matrix_type), pointer :: m_ptr
    failure=.false.
    m_ptr => matrix

    CPPrecondition(associated(m_ptr),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPreconditionNoFail(matrix%ref_count>0,cp_warning_level,routineP,error)
       matrix%ref_count=matrix%ref_count+1
    end if
  end subroutine cp_sp_retain
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_release [1.0] *
!!
!!   NAME
!!     cp_sp_release
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_release(matrix, error)
!!       Type(cp_sparse_matrix_type), Pointer:: matrix
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
  subroutine cp_sp_release(matrix, error)
    type(cp_sparse_matrix_type), pointer :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_release',&
         routineP=moduleN//':'//routineN
    integer :: stat
    failure=.false.

    if (associated(matrix)) then
       CPPreconditionNoFail(matrix%ref_count>0,cp_warning_level,routineP,error)
       matrix%ref_count=matrix%ref_count-1
       if (matrix%ref_count==0) then
          call cp_sp_dealloc(matrix,error=error)
       end if
    end if
  end subroutine cp_sp_release
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_get_matrix [1.0] *
!!
!!   NAME
!!     cp_sp_get_matrix
!!
!!   SYNOPSIS
!!     Function cp_sp_get_matrix(matrix, error) Result(res)
!!       Type(cp_sparse_matrix_type), Intent (IN), Target:: matrix
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
  function cp_sp_get_matrix(matrix, error) result(res)
    type(cp_sparse_matrix_type), intent(in), target :: matrix
    type(real_matrix_type), pointer ::res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_get_matrix',&
         routineP=moduleN//':'//routineN
    type(cp_sparse_matrix_type), pointer :: m_ptr
    failure=.false.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    res => matrix%matrix
  end function cp_sp_get_matrix
  !***************************************************************************

  !================ block matrix op =============
  ! operations related with the fact that we have a block matrix and we want
  ! to work with the blocks

!!****f* cp_sparse_matrix/cp_sp_consolidate_matrix [1.0] *
!!
!!   NAME
!!     cp_sp_consolidate_matrix
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_consolidate_matrix(matrix, error)
!!       Type(cp_sparse_matrix_type), Intent (INOUT):: matrix
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_consolidate_matrix
!!
!!   FUNCTION
!!     consolidates the cache of this matrix. Must be called before any get
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
!!     @author Fawzi Mohamed
!!     @version 3.2002
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sp_consolidate_matrix(matrix,error)
    type(cp_sparse_matrix_type), intent(inout) :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_consolidate_matrix',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug .and. debug_this_module) then
       CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
       call mp_sync(matrix%global_env%group)
    end if
    matrix%cache_consolidated=.true.
  end subroutine cp_sp_consolidate_matrix
  !***************************************************************************

!!FM  !!****f* cp_sparse_matrix/cp_sp_flush_cache [1.0] *
!!FM  !!
!!FM  !!   NAME
!!FM  !!     cp_sp_flush_cache
!!FM  !!
!!FM  !!   FUNCTION
!!FM  !!     removes all, or part of the cached data
!!FM  !!
!!FM  !!   NOTES
!!FM  !!     -
!!FM  !!
!!FM  !!   INPUTS
!!FM  !!     matrix: the matrix that should be purged
!!FM  !!     keep_local_scratch: if true keeps the local scratch blocks 
!!FM  !!		(defaults to false)
!!FM  !!     keep_valid_local_copy: if true keeps the valid local copies
!!FM  !!		(defaults to false)
!!FM  !!     keep_local_copy: if true keeps all the local copies
!!FM  !!		(defaults to false)
!!FM  !!     error: variable to control error logging, stopping,... 
!!FM  !!            see module cp_error_handling 
!!FM  !!
!!FM  !!   AUTHOR
!!FM  !!     @author Fawzi Mohamed
!!FM  !!     @version 3.2002
!!FM  !!
!!FM  !!   MODIFICATION HISTORY
!!FM  !!     none
!!FM  !!
!!FM  !!*** **********************************************************************
!!FM  subroutine cp_sp_flush_cache(matrix, keep_valid_local_copy, keep_scratch, &
!!FM       keep_local_copy,error)
!!FM    type(cp_sparse_matrix_type), intent(in) :: matrix
!!FM    type(cp_error_type), optional, intent(inout) :: error
!!FM
!!FM    logical :: failure
!!FM    integer :: handle
!!FM    character(len=*), parameter :: routineN='cp_sp_flush_cache',&
!!FM         routineP=moduleN//':'//routineN
!!FM    failure=.false.
!!FM
!!FM    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
!!FM    if (.not. failure) then
!!FM       ! to do
!!FM    end if
!!FM  end subroutine cp_sp_flush_cache
!!FM  !***************************************************************************

!!****f* cp_sparse_matrix_type/cp_sp_get_local_block [1.0] *
!!
!!   NAME
!!     cp_sp_get_local_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_get_local_block(matrix, block_row, block_col,&
!!         block_handle, trans, trans_to_do, block_val, local_scratch,&
!!         local_copy, local_owner, owner, error)
!!       Type(cp_sparse_matrix_type), Intent (INOUT), Target:: matrix
!!       Integer, Intent (IN):: block_row, block_col
!!       Type(cp_sp_block_handle_type), Intent (OUT):: block_handle
!!       Character(Len=1), Intent (IN), Optional:: trans
!!       Character(Len=1), Intent (OUT), Optional:: trans_to_do
!!       Real(Kind=wp), Dimension(:,:), Optional, Pointer:: block_val
!!       Logical, Intent (OUT), Optional:: local_scratch, local_copy,&
!!         local_owner
!!       Integer, Intent (OUT), Optional:: owner
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
!!     matrix: matrix to get the info from
!!     block_row: the block row to get
!!     block_col: the block col to get
!!     trans: the transposition wanted (N: None, T: Transpose,
!!            C: Conjugate transpose, B: Bar (just conjugate, no transpose))
!!     block_handle: a handle that is initialized with the info 
!!            regarding the block
!!     block_val: the values of the block
!!     trans_to_do: the transition still to perfomr (defined like trans)
!!            if not asked all the transpositions are always performed
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
  subroutine cp_sp_get_local_block(matrix,block_row,block_col,&
       block_handle,trans,trans_to_do,block_val,&
       local_scratch,local_copy,local_owner,owner,error)
    type(cp_sparse_matrix_type), intent(inout), target :: matrix
    integer, intent(in) :: block_row,block_col
    type(cp_sp_block_handle_type), intent(out) :: block_handle
    character(len=1), intent(in),optional :: trans
    character(len=1), intent(out),optional :: trans_to_do
    real(kind=wp), dimension(:,:), optional, pointer :: block_val
    logical, intent(out), optional :: local_scratch, local_copy, local_owner
    integer, intent(out), optional :: owner
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_get_local_block',&
         routineP=moduleN//':'//routineN
    integer :: i, stat
    real(kind=wp), dimension(:,:), pointer :: tmp_block
    failure=.false.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)   
    CPPreconditionNoFail(.not.associated(block_handle%block_val),cp_warning_level,routineP,error)
    if (.not. failure) then
       CPPreconditionNoFail(.not.matrix%frozen,cp_warning_level,routineP,error)
       block_handle%trans='N'
       if (present(trans)) block_handle%trans=trans
       call uppercase(block_handle%trans)
       block_handle%trans_to_do=block_handle%trans
       call cp_assert(block_handle%trans=='N'.or.block_handle%trans=='T',&
            cp_warning_level,cp_precondition_failed,routineP,&
            "unknown requested transposition '"//block_handle%trans//"' in "//&
            CPSourceFileRef,&
            error=error)
       block_handle%block_row=block_row
       block_handle%block_col=block_col
       if (matrix%matrix_struct%symmetric .and. block_row > block_col) then
          if (block_handle%trans=='N') then
             block_handle%trans_to_do='T'
          else
             block_handle%trans_to_do='N'
          end if
          call get_block_node(matrix%matrix, block_row=block_handle%block_col, &
               block_col=block_handle%block_row, block=tmp_block)
       else
          call get_block_node(matrix%matrix, block_row=block_handle%block_row, &
               block_col=block_handle%block_col, block=tmp_block)
       end if
       if (present(trans_to_do).or.block_handle%trans_to_do=='N') then
          block_handle%block_val => tmp_block
          block_handle%should_dealloc_data=.false.
       else
          block_handle%trans_to_do='N'
          if (associated(tmp_block)) then
             allocate(block_handle%block_val(size(tmp_block,2),&
                  size(tmp_block,1)), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                do i=1,size(tmp_block,2)
                   block_handle%block_val(i,:)=tmp_block(:,i)
                end do
             end if
          else
             nullify(block_handle%block_val)
          end if
          block_handle%should_dealloc_data=.true.
       end if
       if (present(trans_to_do)) trans_to_do = block_handle%trans_to_do
       if (present(block_val)) block_val => block_handle%block_val
       block_handle%initialized=.true.
    end if
  end subroutine cp_sp_get_local_block
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_release_block [1.0] *
!!
!!   NAME
!!     cp_sp_release_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_release_block(block_handle, changed_values, error)
!!       Type(cp_sp_block_handle_type), Intent (INOUT):: block_handle
!!       Logical, Intent (IN):: changed_values
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_release_block
!!
!!   FUNCTION
!!     Releases the memory that might have been allocated to respond
!!     to the get. MUST be called after a block that was get is no longer
!!     used. If some values were changed in the array then changed_values
!!     should be true.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     block_handle: the block handle to be released
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
  subroutine cp_sp_release_block(block_handle,changed_values,error)
    type(cp_sp_block_handle_type), intent(inout) :: block_handle
    logical, intent(in) :: changed_values
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_sp_release_block',&
         routineP=moduleN//':'//routineN
    character :: tt
    failure=.false.

    CPPrecondition(block_handle%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (changed_values .and. block_handle%should_dealloc_data) then
          if (block_handle%trans_to_do==block_handle%trans) then
             tt='N'
          else
             tt='T'
          end if
          call cp_sp_set_local_block(block_handle%matrix, &
               block_row=block_handle%block_row,&
               block_col=block_handle%block_col,&
               trans=tt,local_scratch=.true.,&
               block_val=block_handle%block_val,&
               error=error)
       end if
       if (block_handle%should_dealloc_data) then
          if (associated(block_handle%block_val)) then
             deallocate(block_handle%block_val,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          end if
       end if
    end if
  end subroutine cp_sp_release_block
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_set_local_block [1.0] *
!!
!!   NAME
!!     cp_sp_set_local_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_set_local_block(matrix, block_row, block_col,&
!!         block_val, trans, local_scratch, error)
!!       Type(cp_sparse_matrix_type), Intent (INOUT), Target:: matrix
!!       Integer, Intent (IN):: block_row, block_col
!!       Real(Kind=wp), Intent (IN), Dimension(:,:), Target:: block_val
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
  subroutine cp_sp_set_local_block(matrix,block_row,block_col,block_val,&
       trans,local_scratch,error)
    type(cp_sparse_matrix_type), intent(inout), target :: matrix
    integer, intent(in) :: block_row, block_col
    real(kind=wp), intent(in), dimension(:,:), target :: block_val
    character, intent(in), optional :: trans
    logical, intent(in), optional :: local_scratch
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_set_local_block',&
         routineP=moduleN//':'//routineN
    character :: t_to_do
    integer :: b_row,b_col, i, stat
    real(kind=wp), dimension(:,:), pointer :: tmp_block
    failure=.false.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       t_to_do='N'
       if (present(trans)) t_to_do=trans
       call uppercase(t_to_do)
       b_row=block_col
       b_row=block_row
       if (matrix%matrix_struct%symmetric.and.block_row > block_col) then
          if (t_to_do=='N') then
             t_to_do='T'
          else
             t_to_do='N'
          end if
          b_row=block_row
          b_row=block_col
       end if
       if (t_to_do=='T') then
          allocate(tmp_block(size(block_val,2),size(block_val,1)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          do i=1,size(block_val,2)
             tmp_block(i,:)=block_val(:,i)
          end do
          call put_block_node(matrix%matrix, block_row=b_row,&
               block_col=b_col, block=tmp_block)
          deallocate(tmp_block,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       else
          tmp_block => block_val
          call put_block_node(matrix%matrix, block_row=b_row,&
               block_col=b_col, block=tmp_block)     
       end if
    endif
  end subroutine cp_sp_set_local_block
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_remove_block [1.0] *
!!
!!   NAME
!!     cp_sp_remove_block
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_remove_block(matrix, block_row, block_col, error)
!!       Type(cp_sparse_matrix_type), Intent (INOUT), Target:: matrix
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
  subroutine cp_sp_remove_block(matrix,block_row, block_col,error)
    type(cp_sparse_matrix_type), intent(inout), target :: matrix
    integer, intent(in):: block_col, block_row
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_remove_block',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CPAssert(.false.,cp_failure_level,routineP,error,failure)
       if (matrix%matrix_struct%symmetric.and.block_row > block_col) then
          ! to do
       else
       end if
    end if
  end subroutine cp_sp_remove_block
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_dealloc_m_array [1.0] *
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
  subroutine cp_sp_array_dealloc(m_array, error)
    type(cp_sparse_matrix_p_type), dimension(:), pointer :: m_array
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_array_dealloc',&
         routineP=moduleN//':'//routineN
    integer :: i, stat
    failure=.false.

    if (associated(m_array)) then
       do i=1,size(m_array)
          call cp_release(m_array(i)%matrix,error=error)
       end do
       deallocate(m_array,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    end if
  end subroutine cp_sp_array_dealloc
  !***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_array_release [1.0] *
!!
!!   NAME
!!     cp_sp_array_release
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_array_release(array, error)
!!       Type(cp_sparse_matrix_p_type), Dimension(:), Intent (INOUT)::&
!!         array
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sp_array_release
!!
!!   FUNCTION
!!     releases all the elements of an array of cp_sparse_matrix_p_type
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
subroutine cp_sp_array_release(array,error)
  type(cp_sparse_matrix_p_type), dimension(:), intent(inout) :: array
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_array_release',&
        routineP=moduleN//':'//routineN
  integer :: i
  failure=.false.
  
  do i=1,size(array)
     call cp_release(array(i)%matrix,error=error)
     nullify(array(i)%matrix)
  end do
end subroutine cp_sp_array_release
!***************************************************************************

!!****f* cp_sparse_matrix/cp_sp_alloc_nonsparse_blocks *
!!
!!   NAME
!!     cp_sp_alloc_nonsparse_blocks
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_alloc_nonsparse_blocks(matrix, error)
!!       Type(cp_sparse_matrix_type), Intent (INOUT), Target:: matrix
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
subroutine cp_sp_alloc_nonsparse_blocks(matrix,error)
  type(cp_sparse_matrix_type), intent(inout), target :: matrix
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_alloc_nonsparse_blocks',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call cp_sm_alloc_nonsparse_blocks(matrix%matrix,&
          m_struct=matrix%matrix_struct,&
          global_env=matrix%global_env,&
          error=error)
  end if
end subroutine cp_sp_alloc_nonsparse_blocks
!***************************************************************************

end module cp_sparse_matrix
