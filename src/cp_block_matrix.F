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
module cp_block_matrix
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger, cp_precondition_failed
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use string_utilities, only: uppercase
  use timings, only: timeset, timestop
  use cp_b_matrix_structure
  use sparse_matrix_types
  use message_passing, only: mp_sync
  use cp_matrix_utils, only: cp_sm_alloc_nonsparse_blocks, &
       cp_matrix_block_iterator, cp_sm_b_i_init, cp_sm_b_i_next, &
       cp_sm_b_i_dealloc_ref, cp_sm_set_to
  use qs_blacs, only:sparse_times_blacs, blacs_matrix_type, get_blacs_matrix_info
  use cp_para_env, only: cp_para_env_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_block_matrix'

  public :: cp_block_matrix_type, cp_sp_block_handle_type,&
       cp_block_matrix_p_type
  public :: cp_create, cp_dealloc, cp_retain, cp_release, &
       cp_get_local_block,&
       cp_release_block, cp_set_local_block, cp_remove_block,&
       cp_get_matrix, cp_set_to, cp_get
  public :: cp_sp_create, cp_sp_retain, cp_sp_release,&
       cp_sp_get_local_block,&
       cp_sp_block_handle_release, cp_sp_set_local_block, cp_sp_remove_block,&
       cp_sp_array_dealloc, cp_sp_get_matrix, cp_sp_array_release,&
       cp_sp_alloc_nonsparse_blocks, cp_sp_set_to, cp_sp_block_handle_init,&
       cp_sp_get, cp_sp_multiplies_blacs, cp_sp_array_create

  interface cp_init
     module procedure cp_sp_block_handle_init
  end interface
  interface cp_create
     module procedure cp_sp_create, cp_sp_array_create
  end interface
!FM  interface cp_dealloc
!FM     module procedure cp_sp_dealloc
!FM  end interface
  interface cp_retain
     module procedure cp_sp_retain
  end interface
  interface cp_release
     module procedure cp_sp_release, cp_sp_array_release
  end interface
  interface cp_get_matrix
     module procedure cp_sp_get_matrix
  end interface
  interface cp_flush_cache
     module procedure cp_sp_flush_cache
  end interface
  interface cp_get_local_block
     module procedure cp_sp_get_local_block
  end interface
  interface cp_release_block
     module procedure cp_sp_block_handle_release
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
  interface cp_set_to
     module procedure cp_sp_set_to
  end interface
  interface cp_get
     module procedure cp_sp_get
  end interface
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
  type cp_block_matrix_type
     logical :: initialized
     logical :: frozen
     logical :: should_dealloc_matrix
     logical :: b_rep_is_primary
     logical :: b_i_rep_in_sync
     logical :: cache_is_valid
     type(cp_b_matrix_struct_type), pointer :: matrix_struct
     integer :: ref_count, id_nr
     type(real_matrix_type), pointer :: sm
  end type cp_block_matrix_type
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
  type cp_block_matrix_p_type
     type(cp_block_matrix_type), pointer :: matrix
  end type cp_block_matrix_p_type
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
  type cp_sp_block_handle_type
     type(cp_block_matrix_type), pointer :: matrix
     integer :: block_row, block_col
     character :: trans,trans_to_do
     real(kind=wp), dimension(:,:), pointer :: block_val
     logical :: should_dealloc_data, initialized
  end type cp_sp_block_handle_type
!!***
!****************************************************************************

contains

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
  subroutine cp_sp_create(matrix, matrix_struct, u_matrix, &
       should_dealloc_matrix, b_rep_is_primary, allocate_nonsparse_blocks, error)
    type(cp_block_matrix_type), pointer :: matrix
    type(cp_b_matrix_struct_type), intent(in), target :: matrix_struct
    type(real_matrix_type), optional, intent(in), target :: u_matrix
    logical, intent(in), optional :: should_dealloc_matrix, b_rep_is_primary,&
        allocate_nonsparse_blocks
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
       matrix%cache_is_valid=.false.
       matrix%should_dealloc_matrix=.true.
       if (present(should_dealloc_matrix)) &
            matrix%should_dealloc_matrix=should_dealloc_matrix
       matrix%b_rep_is_primary = .true.
       if (present(b_rep_is_primary)) matrix%b_rep_is_primary=b_rep_is_primary
       matrix%b_i_rep_in_sync = .false.
       matrix%matrix_struct => matrix_struct
       call cp_retain(matrix%matrix_struct,error=error)
       matrix%ref_count=1
       nullify(matrix%sm)
       matrix%initialized=.true.
       if (present(u_matrix)) then
          matrix%sm => u_matrix
       else
          if (matrix%matrix_struct%symmetric) then
             m_symm='symmetric'
          else
             m_symm='none'
          end if
          b_dims => matrix%matrix_struct%b_dims
          block_rows=size(b_dims%block_begins_at_row)-1
          block_cols=size(b_dims%block_begins_at_col)-1
          call allocate_matrix(matrix%sm,&
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
       if (present(allocate_nonsparse_blocks)) then
          if (allocate_nonsparse_blocks) then
             call cp_sp_alloc_nonsparse_blocks(matrix,error=error)
          end if
       end if
    end if
  end subroutine cp_sp_create
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
  subroutine cp_sp_dealloc(matrix, error)
    type(cp_block_matrix_type), pointer :: matrix
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
          if (matrix%should_dealloc_matrix) then
             if (associated(matrix%sm)) then
                call deallocate_matrix(matrix%sm)
             end if
          else
             nullify(matrix%sm)
          end if
          deallocate(matrix,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       end if
    end if
  end subroutine cp_sp_dealloc
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
  subroutine cp_sp_retain(matrix, error)
    type(cp_block_matrix_type), intent(inout), target :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_retain',&
         routineP=moduleN//':'//routineN
    type(cp_block_matrix_type), pointer :: m_ptr
    failure=.false.
    m_ptr => matrix

    CPPrecondition(associated(m_ptr),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       CPPreconditionNoFail(matrix%ref_count>0,cp_warning_level,routineP,error)
       matrix%ref_count=matrix%ref_count+1
    end if
  end subroutine cp_sp_retain
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
  subroutine cp_sp_release(matrix, error)
    type(cp_block_matrix_type), pointer :: matrix
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
    nullify(matrix)
  end subroutine cp_sp_release
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
  function cp_sp_get_matrix(matrix, error) result(res)
    type(cp_block_matrix_type), pointer :: matrix
    type(real_matrix_type), pointer ::res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_get_matrix',&
         routineP=moduleN//':'//routineN
    type(cp_block_matrix_type), pointer :: m_ptr
    failure=.false.
    
    CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    res => matrix%sm
  end function cp_sp_get_matrix
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
  subroutine cp_sp_set_b_rep(matrix,error)
    type(cp_block_matrix_type), pointer :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_set_b_rep',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)

    if (.not.failure) then
       if (.not.matrix%b_rep_is_primary) then
          if (cp_debug .and. debug_this_module) then
             call mp_sync(matrix%matrix_struct%distribution%para_env%group)
          end if
          matrix%b_rep_is_primary=.true.
       end if
    end if
  end subroutine cp_sp_set_b_rep
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
  subroutine cp_sp_set_i_rep(matrix,error)
    type(cp_block_matrix_type), pointer :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_set_b_rep',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)

    if (.not.failure) then
       if (matrix%b_rep_is_primary) then
          if (cp_debug .and. debug_this_module) then
             call mp_sync(matrix%matrix_struct%distribution%para_env%group)
          end if
          matrix%b_rep_is_primary=.false.
       end if
    end if
  end subroutine cp_sp_set_i_rep
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
  subroutine cp_sp_flush_cache(matrix,flush_secondary_rep,error)
    type(cp_block_matrix_type), intent(in) :: matrix
    logical, optional, intent(in) :: flush_secondary_rep
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: brow,bcol, mepos
    character(len=*), parameter :: routineN='cp_sp_flush_cache',&
         routineP=moduleN//':'//routineN
    type(cp_matrix_block_iterator) :: b_iterator
    type(cp_b_distribution_type), pointer :: distribution
    nullify(distribution)
    failure=.false.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_sm_b_i_init(b_iterator,matrix=matrix%sm,error=error)
       distribution => matrix%matrix_struct%distribution
       mepos=matrix%matrix_struct%distribution%para_env%mepos
       CPPrecondition(associated(distribution),cp_failure_level,routineP,error,failure)
       do
          if (.not.cp_sm_b_i_next(b_iterator,block_row=brow,&
               block_col=bcol)) exit
          if (distribution%block_owner(brow,bcol) /= mepos) then
             call remove_block_node(matrix=matrix%sm, block_row=brow,&
                  block_col=bcol,error=error)
          end if
       end do
       call cp_sm_b_i_dealloc_ref(b_iterator,error=error)
    end if
  end subroutine cp_sp_flush_cache
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
subroutine cp_sp_block_handle_init(block_handle,error)
  type(cp_sp_block_handle_type), intent(out) :: block_handle
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_block_handle_init',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  block_handle%block_col=-1
  block_handle%block_row=-1
  block_handle%trans='N'
  block_handle%trans_to_do='N'
  block_handle%should_dealloc_data=.false.
  nullify(block_handle%block_val, block_handle%matrix)
  block_handle%initialized=.true.
end subroutine cp_sp_block_handle_init
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
!!       Real(Kind=wp), Dimension(:,:), Optional, Pointer:: block_val
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
  subroutine cp_sp_get_local_block(matrix,block_row,block_col,&
       block_handle,trans,trans_to_do,block_val,&
       local_owner,owner,alloc_block,error)
    type(cp_block_matrix_type), intent(inout), target :: matrix
    integer, intent(in) :: block_row,block_col
    type(cp_sp_block_handle_type), intent(inout) :: block_handle
    character(len=1), intent(in),optional :: trans
    character(len=1), intent(out),optional :: trans_to_do
    real(kind=wp), dimension(:,:), optional, pointer :: block_val
    logical, intent(out), optional :: local_owner
    integer, intent(out), optional :: owner
    logical, intent(in), optional :: alloc_block
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_get_local_block',&
         routineP=moduleN//':'//routineN
    integer :: i, stat
    real(kind=wp), dimension(:,:), pointer :: tmp_block
    failure=.false.

    CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(block_handle%initialized,cp_warning_level,routineP,error)
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
          call get_block_node(matrix%sm, block_row=block_handle%block_col, &
               block_col=block_handle%block_row, block=tmp_block)
          if (present(alloc_block)) then
             if (alloc_block .and. .not.associated(tmp_block)) then
                call add_block_node(matrix%sm,&
                     block_row=block_handle%block_col,&
                     block_col=block_handle%block_row, block=tmp_block)
             end if
          end if
       else
          call get_block_node(matrix%sm, block_row=block_handle%block_row, &
               block_col=block_handle%block_col, block=tmp_block)
          if (present(alloc_block)) then
             if (alloc_block .and. .not.associated(tmp_block)) then
                call add_block_node(matrix%sm, block_row=block_handle%block_row,&
                     block_col=block_handle%block_col, block=tmp_block)
             end if
          end if
       end if
       if (present(trans_to_do).or.(.not.present(block_val)).or.&
            block_handle%trans_to_do=='N') then
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
       if (present(local_owner)) then
          local_owner=matrix%matrix_struct%distribution%block_owner&
               (block_handle%block_row,block_handle%block_col)==&
               matrix%matrix_struct%distribution%para_env%mepos
       end if
       if (present(owner)) then
          owner=matrix%matrix_struct%distribution%block_owner&
               (block_handle%block_row,block_handle%block_col)
       end if
       block_handle%initialized=.true.
    end if
  end subroutine cp_sp_get_local_block
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
  subroutine cp_sp_block_handle_release(block_handle,changed_values,error)
    type(cp_sp_block_handle_type), intent(inout) :: block_handle
    logical, intent(in) :: changed_values
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_sp_block_handle_release',&
         routineP=moduleN//':'//routineN
    character :: tt
    failure=.false.

    CPPrecondition(block_handle%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (changed_values.and.block_handle%matrix%frozen) then
          call cp_assert(.false.,&
               cp_warning_level,cp_assertion_failed,routineP,&
               "changed values in frozen matrix",error)
          block_handle%matrix%frozen=.false.
       end if
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
       nullify(block_handle%block_val)
    end if
  end subroutine cp_sp_block_handle_release
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
    type(cp_block_matrix_type), intent(inout), target :: matrix
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
    if (.not.failure) then
       call cp_assert(.not.matrix%frozen,&
            cp_warning_level,cp_assertion_failed,routineP,&
            "changed values in frozen matrix",error)
       matrix%frozen=.false.
    end if
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
          call put_block_node(matrix%sm, block_row=b_row,&
               block_col=b_col, block=tmp_block)
          deallocate(tmp_block,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       else
          tmp_block => block_val
          call put_block_node(matrix%sm, block_row=b_row,&
               block_col=b_col, block=tmp_block)     
       end if
    endif
  end subroutine cp_sp_set_local_block
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
  subroutine cp_sp_remove_block(matrix,block_row, block_col,error)
    type(cp_block_matrix_type), intent(inout), target :: matrix
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
          call remove_block_node(matrix%sm,block_row=block_col,&
               block_col=block_row, error=error)
       else
          call remove_block_node(matrix%sm,block_row=block_row,&
               block_col=block_col, error=error)
       end if
    end if
  end subroutine cp_sp_remove_block
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
  subroutine cp_sp_array_create(matrixes, n_matrixes, matrix_struct, &
      b_rep_is_primary, allocate_nonsparse_blocks, error)
    type(cp_block_matrix_p_type), dimension(:), pointer :: matrixes
    type(cp_b_matrix_struct_type), intent(in), target :: matrix_struct
    integer, intent(in) :: n_matrixes
    logical, intent(in), optional :: b_rep_is_primary,&
        allocate_nonsparse_blocks
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_array_create',&
         routineP=moduleN//':'//routineN
    integer :: stat,i
    failure=.false.
    
    allocate(matrixes(n_matrixes),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       do i=1,n_matrixes
          nullify(matrixes(i)%matrix)
          call cp_sp_create(matrixes(i)%matrix,&
              matrix_struct=matrix_struct,&
              b_rep_is_primary=b_rep_is_primary,&
              allocate_nonsparse_blocks=allocate_nonsparse_blocks,&
              error=error)
       end do
    end if
  end subroutine cp_sp_array_create

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
  subroutine cp_sp_array_dealloc(m_array, error)
    type(cp_block_matrix_p_type), dimension(:), pointer :: m_array
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
    nullify(m_array)
  end subroutine cp_sp_array_dealloc
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
subroutine cp_sp_array_release(array,error)
  type(cp_block_matrix_p_type), dimension(:), intent(inout) :: array
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
subroutine cp_sp_alloc_nonsparse_blocks(matrix,error)
  type(cp_block_matrix_type), intent(inout), target :: matrix
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_alloc_nonsparse_blocks',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call cp_sm_alloc_nonsparse_blocks(matrix%sm,&
          m_struct=matrix%matrix_struct,&
          para_env=matrix%matrix_struct%distribution%para_env,&
          error=error)
  end if
end subroutine cp_sp_alloc_nonsparse_blocks
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
subroutine cp_sp_update_sparsity(matrix,error)
  type(cp_block_matrix_type), intent(inout), target :: matrix
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_update_sparsity',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(fromWhere=routineP, &
          message="update sparsity", error=error)
     ! to do
  end if
end subroutine cp_sp_update_sparsity
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
subroutine cp_sp_check_sparsity(matrix,error)
  type(cp_block_matrix_type), intent(inout), target :: matrix
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_check_sparsity',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(fromWhere=routineP, &
          message="check sparsity", error=error)
     ! to do
  end if
end subroutine cp_sp_check_sparsity
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
subroutine cp_sp_enforce_sparsity(matrix,error)
  type(cp_block_matrix_type), intent(inout), target :: matrix
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_enforce_sparsity',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(fromWhere=routineP, &
          message="enforce sparsity", error=error)
     ! to do
  end if
end subroutine cp_sp_enforce_sparsity
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
subroutine cp_sp_check_distribution(matrix,error)
  type(cp_block_matrix_type), intent(inout), target :: matrix
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_check_distribution',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(fromWhere=routineP, &
          message="check distribution", error=error)
     ! to do
  end if
end subroutine cp_sp_check_distribution
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
subroutine cp_sp_get(matrix,matrix_struct, distribution,dimensions,&
     sparsity, n_block_rows, n_block_cols, n_rows, n_cols, frozen,&
     b_rep_is_primary, id_nr, error)
  type(cp_block_matrix_type), pointer :: matrix
  type(cp_b_matrix_struct_type), pointer, optional :: matrix_struct
  type(cp_b_distribution_type), pointer, optional :: distribution
  type(cp_b_dims_type), pointer, optional :: dimensions
  logical, dimension(:,:), pointer, optional :: sparsity
  integer, intent(out), optional :: n_block_cols, n_block_rows,&
       n_rows, n_cols, id_nr
  logical, intent(out), optional :: frozen, b_rep_is_primary
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_get',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     CPPrecondition(matrix%ref_count>0,cp_failure_level,routineP,error,failure)
  end if
  if (.not. failure) then
     if (present(matrix_struct)) matrix_struct => matrix%matrix_struct
     if (present(distribution)) distribution => &
          matrix%matrix_struct%distribution
     if (present(dimensions)) dimensions => matrix%matrix_struct%b_dims
     if (present(sparsity)) sparsity => matrix%matrix_struct%sparsity
     if (present(n_block_rows)) n_block_rows = &
          size(matrix%matrix_struct%b_dims%block_begins_at_row)-1
     if (present(n_block_cols)) n_block_cols = &
          size(matrix%matrix_struct%b_dims%block_begins_at_col)-1
     if (present(n_rows)) n_rows = &
          matrix%matrix_struct%b_dims%block_begins_at_row &
          (size(matrix%matrix_struct%b_dims%block_begins_at_row))-1
     if (present(n_cols)) n_cols = &
          matrix%matrix_struct%b_dims%block_begins_at_col &
          (size(matrix%matrix_struct%b_dims%block_begins_at_col))-1
     if (present(frozen)) frozen=matrix%frozen
     if (present(b_rep_is_primary)) b_rep_is_primary=matrix%b_rep_is_primary
     if (present(id_nr)) id_nr=matrix%id_nr
  end if
end subroutine cp_sp_get
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
subroutine cp_sp_redistribute(matrix, new_distribution, error)
  type(cp_block_matrix_type), intent(inout), target :: matrix
  type(cp_b_distribution_type), pointer :: new_distribution
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_redistribute',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (.not. failure) then
     call cp_unimplemented_error(fromWhere=routineP, &
          message="redistribute", error=error)
     ! to do
  end if
end subroutine cp_sp_redistribute
!***************************************************************************\

!!****f* cp_block_matrix/cp_sp_set_to [1.0] *
!!
!!   NAME
!!     cp_sp_set_to
!!
!!   SYNOPSIS
!!     Subroutine cp_sp_set_to(matrix, value, error)
!!       Type(cp_block_matrix_type), Pointer:: matrix
!!       Real(Kind=wp), Intent (IN), Optional:: value
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
subroutine cp_sp_set_to(matrix, value, error)
  type(cp_block_matrix_type), pointer :: matrix
  real(kind=wp), intent(in), optional :: value
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sp_set_to',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(associated(matrix),cp_failure_level,routineP,error,failure)
  CPPrecondition(matrix%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     call cp_sm_set_to(matrix%sm, value=value, error=error)
  end if
end subroutine cp_sp_set_to
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
!!       Type(blacs_matrix_type), Pointer:: matrix_b, matrix_c
!!       Real(Kind=wp), Intent (IN), Optional:: alpha, beta
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
  subroutine cp_sp_multiplies_blacs(matrix_a,matrix_b,matrix_c,alpha,beta,error)
    type(cp_block_matrix_type), pointer :: matrix_a
    type(blacs_matrix_type), pointer :: matrix_b, matrix_c
    real(kind=wp), intent(in), optional :: alpha, beta
    type(cp_error_type), intent(inout), optional :: error

    integer :: ncol
    
    call get_blacs_matrix_info(matrix_b,ncol_global=ncol)
    call sparse_times_blacs(matrix_a%sm, matrix_b, matrix_c,ncol,&
         para_env=matrix_a%matrix_struct%distribution%para_env,&
         alpha=alpha, beta=beta,&
         block_owner=matrix_a%matrix_struct%distribution%block_owner)
  end subroutine cp_sp_multiplies_blacs

end module cp_block_matrix
