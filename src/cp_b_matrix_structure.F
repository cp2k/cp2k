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
module cp_b_matrix_structure
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use message_passing, only: mp_sum
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_b_matrix_structure'

  public :: cp_b_distribution_type,cp_b_dims_type,&
       cp_b_matrix_struct_type, cp_bs_iterator_type

  public :: cp_create, cp_dealloc, cp_retain, cp_release, cp_get, cp_next,&
       cp_init, cp_dealloc_ref
  public :: cp_freeze_sparsity, cp_get_sparsity_frozen, cp_get_block_info

  public :: cp_dstruct_create, cp_dstruct_dealloc, cp_dstruct_retain,&
       cp_dstruct_release, cp_dstruct_get, &
       cp_b_dims_create, cp_b_dims_dealloc, cp_b_dims_retain,&
       cp_b_dims_release,&
       cp_bmstruct_create, cp_bmstruct_dealloc, cp_bmstruct_freeze_sparsity,&
       cp_bmstruct_get_sparsity_frozen, cp_bmstruct_retain,&
       cp_bmstruct_release, cp_bmstruct_get_block_info,&
       cp_bs_iter_init, cp_bs_iter_dealloc_ref, cp_bs_iter_next,&
       cp_bmstruct_get

  interface cp_create
     module procedure cp_dstruct_create
  end interface
  interface cp_dealloc
     module procedure cp_dstruct_dealloc
  end interface
  interface cp_retain
     module procedure cp_dstruct_retain
  end interface
  interface cp_release
     module procedure cp_dstruct_release
  end interface
  interface cp_get
     module procedure cp_dstruct_get
  end interface

  interface cp_create
     module procedure cp_b_dims_create
  end interface
  interface cp_dealloc
     module procedure cp_b_dims_dealloc
  end interface
  interface cp_retain
     module procedure cp_b_dims_retain
  end interface
  interface cp_release
     module procedure cp_b_dims_release
  end interface

  interface cp_create
     module procedure cp_bmstruct_create
  end interface
  interface cp_dealloc
     module procedure cp_bmstruct_dealloc
  end interface
  interface cp_freeze_sparsity
     module procedure cp_bmstruct_freeze_sparsity
  end interface
  interface cp_get
     module procedure cp_bmstruct_get
  end interface
  interface cp_get_sparsity_frozen
     module procedure cp_bmstruct_get_sparsity_frozen
  end interface
  interface cp_retain
     module procedure cp_bmstruct_retain
  end interface
  interface cp_release
     module procedure cp_bmstruct_release
  end interface
  interface cp_get_block_info
     module procedure cp_bmstruct_get_block_info
  end interface

  interface cp_init
     module procedure cp_bs_iter_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_bs_iter_dealloc_ref
  end interface
  interface cp_next
     module procedure cp_bs_iter_next
  end interface
  !***
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
!!     global_env: the global environement with the mpi communicator,...
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
  type cp_b_distribution_type
     logical :: initialized
     type(global_environment_type), pointer :: global_env
     integer, dimension(:,:), pointer :: block_owner
     integer :: ref_count, id_nr
  end type cp_b_distribution_type
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
  type cp_b_dims_type
     logical :: initialized
     integer, dimension(:), pointer :: block_begins_at_row, block_begins_at_col
     integer :: ref_count, id_nr
  end type cp_b_dims_type
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
  type cp_b_matrix_struct_type
     logical :: initialized, sparsity_frozen, symmetric
     type(cp_b_dims_type), pointer :: b_dims
     type(cp_b_distribution_type), pointer :: distribution
     logical, dimension(:,:), pointer :: sparsity
     integer :: ref_count, id_nr
  end type cp_b_matrix_struct_type
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
  type cp_bs_iterator_type
     logical :: initialized, do_global, do_zero, do_lower_half
     integer :: block_row,block_col,n_block_rows, n_block_cols
     type(cp_b_matrix_struct_type), pointer :: m_struct
  end type cp_bs_iterator_type
!!***
  !****************************************************************************

contains

  !======= block distribution structure =======

!!****f* cp_b_matrix_structure/cp_dstruct_create [1.0] *
!!
!!   NAME
!!     cp_dstruct_create
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_create(d_struct, global_env, n_block_rows,&
!!         n_block_cols, error)
!!       Type(cp_b_distribution_type), Pointer:: d_struct
!!       Type(global_environment_type), Intent (IN), Target:: global_env
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
!!     global_env: the global environment
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
  subroutine cp_dstruct_create(d_struct,global_env,n_block_rows,n_block_cols,&
       error)
    type(cp_b_distribution_type), pointer :: d_struct
    type(global_environment_type), intent(in), target :: global_env
    integer, intent(in) :: n_block_rows, n_block_cols
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    integer, save :: last_id_nr=0
    character(len=*), parameter :: routineN='cp_b_distribution_type',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(n_block_cols>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(n_block_rows>0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       allocate(d_struct,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       last_id_nr=last_id_nr+1
       d_struct%global_env => global_env
       allocate(d_struct%block_owner(n_block_rows,n_block_cols),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       d_struct%block_owner=global_env%source
       d_struct%ref_count=1
       d_struct%id_nr=last_id_nr
       d_struct%initialized=.true.
    end if
  end subroutine cp_dstruct_create
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
  subroutine cp_dstruct_dealloc(d_struct,error)
    type(cp_b_distribution_type), pointer :: d_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_dstruct_dealloc',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (associated(d_struct)) then
       CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          CPPreconditionNoFail(d_struct%ref_count==0,cp_warning_level,routineP,error)
          d_struct%ref_count=-1
          if (associated(d_struct%block_owner)) then
             deallocate(d_struct%block_owner,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
          nullify(d_struct%global_env)
          d_struct%initialized=.false.
       end if
       deallocate(d_struct,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    end if
  end subroutine cp_dstruct_dealloc
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
  subroutine cp_dstruct_retain(d_struct,error)
    type(cp_b_distribution_type), pointer :: d_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_dstruct_retain',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug) then
       CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(d_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       d_struct%ref_count=d_struct%ref_count+1
    end if
  end subroutine cp_dstruct_retain
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
  subroutine cp_dstruct_release(d_struct,error)
    type(cp_b_distribution_type), pointer :: d_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_dstruct_release',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug) then
       CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(d_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       d_struct%ref_count=d_struct%ref_count-1
       if (d_struct%ref_count<1) then
          call cp_dstruct_dealloc(d_struct)
       end if
    end if
  end subroutine cp_dstruct_release
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_dstruct_get [1.0] *
!!
!!   NAME
!!     cp_dstruct_get
!!
!!   SYNOPSIS
!!     Subroutine cp_dstruct_get(d_struct, block_owner, global_env, error)
!!       Type(cp_b_distribution_type), Intent (IN):: d_struct
!!       Integer, Dimension(:,:), Pointer, Optional:: block_owner
!!       Type(global_environment_type), Pointer, Optional:: global_env
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
!!     global_env: the global environement
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
  subroutine cp_dstruct_get(d_struct, block_owner, global_env, error)
    type(cp_b_distribution_type), intent(in) :: d_struct
    integer, dimension(:,:), pointer, optional :: block_owner
    type(global_environment_type), pointer, optional :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_dstruct_get',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(d_struct%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(block_owner)) block_owner => d_struct%block_owner
       if (present(global_env)) global_env => d_struct%global_env
    end if
  end subroutine cp_dstruct_get
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
  subroutine cp_b_dims_create(b_dims, b_row_sizes, b_col_sizes, error)
    type(cp_b_dims_type), pointer :: b_dims
    integer, dimension(:), intent(in) :: b_row_sizes, b_col_sizes
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat,i
    integer, save :: last_id_nr=0
    character(len=*), parameter :: routineN='cp_b_dims_create',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(size(b_row_sizes)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(size(b_col_sizes)>0,cp_failure_level,routineP,error,failure)
    if(cp_debug) then
       CPPrecondition(all(b_row_sizes>=0),cp_failure_level,routineP,error,failure)
       CPPrecondition(all(b_col_sizes>=0),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       allocate (b_dims,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       last_id_nr=last_id_nr+1
       b_dims%id_nr=last_id_nr
       b_dims%ref_count=1
       allocate (b_dims%block_begins_at_row(size(b_row_sizes)+1),&
            b_dims%block_begins_at_col(size(b_col_sizes)+1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       b_dims%block_begins_at_row(1)=1
       do i=1,size(b_row_sizes)
          b_dims%block_begins_at_row(i+1)=b_dims%block_begins_at_row(i)+&
               b_row_sizes(i)
       end do
       b_dims%block_begins_at_col(1)=1
       do i=1,size(b_col_sizes)
          b_dims%block_begins_at_col(i+1)=b_dims%block_begins_at_col(i)+&
               b_col_sizes(i)
       end do
       b_dims%initialized=.true.
    end if
  end subroutine cp_b_dims_create
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
  subroutine cp_b_dims_dealloc(b_dims,error)
    type(cp_b_dims_type), pointer :: b_dims
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_b_dims_dealloc',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(b_dims%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CPPreconditionNoFail(b_dims%ref_count==0,cp_warning_level,routineP,error)
       if (associated(b_dims)) then
          if (associated(b_dims%block_begins_at_row)) then
             deallocate(b_dims%block_begins_at_row, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
          if (associated(b_dims%block_begins_at_col)) then
             deallocate(b_dims%block_begins_at_col,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
          b_dims%initialized=.false.
          deallocate(b_dims,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
    end if
  end subroutine cp_b_dims_dealloc
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
  subroutine cp_b_dims_retain(b_dims,error)
    type(cp_b_dims_type), intent(inout) :: b_dims
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_b_dims_retain',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(b_dims%initialized,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(b_dims%ref_count>0,cp_warning_level,routineP,error)
    if (.not. failure) then
       b_dims%ref_count=b_dims%ref_count+1
    end if
  end subroutine cp_b_dims_retain
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
  subroutine cp_b_dims_release(b_dims,error)
    type(cp_b_dims_type), pointer :: b_dims
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_b_dims_release',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug) then
       CPPrecondition(b_dims%initialized,cp_failure_level,routineP,error,failure)
       CPPrecondition(b_dims%ref_count>0,cp_failure_level,routineP,error,failure)
    end if
    if (.not.failure) then
       b_dims%ref_count=b_dims%ref_count-1
       if (b_dims%ref_count<1) then
          call cp_b_dims_dealloc(b_dims,error)
       end if
    end if
  end subroutine cp_b_dims_release
  !***************************************************************************

  !======== matrix structure ========
!!****f* cp_b_matrix_structure/cp_bmstruct_create [1.0] *
!!
!!   NAME
!!     cp_bmstruct_create
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_create(m_struct, symmetric, d_struct, b_dims,&
!!         sparsity, no_sparsity, global_env, b_row_sizes, b_col_sizes,&
!!         error)
!!       Type(cp_b_matrix_struct_type), Pointer:: m_struct
!!       Type(cp_b_distribution_type), Target, Intent (INOUT), Optional::&
!!         d_struct
!!       Type(cp_b_dims_type), Target, Intent (INOUT), Optional:: b_dims
!!       Logical, Dimension(:,:), Optional:: sparsity
!!       Logical, Intent (IN), Optional:: no_sparsity, symmetric
!!       Type(global_environment_type), Target, Intent (IN), Optional::&
!!         global_env
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
!!     global_env: the global environement
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
  subroutine cp_bmstruct_create(m_struct,symmetric,d_struct,b_dims, &
       sparsity, no_sparsity,global_env, b_row_sizes, b_col_sizes,error)
    type(cp_b_matrix_struct_type), pointer :: m_struct
    type(cp_b_distribution_type), target, intent(inout), optional :: d_struct
    type(cp_b_dims_type), target, intent(inout), optional :: b_dims
    logical, dimension(:,:), optional :: sparsity !, intent(in)
    logical, intent(in), optional :: no_sparsity, symmetric
    type(global_environment_type), target, intent(in), optional :: global_env
    integer, dimension(:), intent(in), optional :: b_row_sizes, b_col_sizes
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, no_spars
    integer :: stat
    integer, save :: last_id_nr=0
    character(len=*), parameter :: routineN='cp_bmstruct_create',&
         routineP=moduleN//':'//routineN
    failure=.false.

    no_spars=.false.
    if (present(no_sparsity)) no_spars=no_sparsity
    allocate(m_struct, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       last_id_nr=last_id_nr+1
       m_struct%id_nr=last_id_nr
       m_struct%sparsity_frozen= present(sparsity).or.no_spars
       m_struct%symmetric=.true. ! make it default to false
       if (present(symmetric)) m_struct%symmetric=symmetric
       m_struct%ref_count=1

       if (present(d_struct)) then
          m_struct%distribution => d_struct
          call cp_retain(m_struct%distribution, error)
       else
          CPPrecondition(present(global_env),cp_failure_level,routineP,error,failure)
          if (.not.present(b_dims)) then
             CPPrecondition(present(b_col_sizes),cp_failure_level,routineP,error,failure)
             CPPrecondition(present(b_row_sizes),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                call cp_dstruct_create(m_struct%distribution,global_env=global_env,&
                     n_block_rows=size(b_row_sizes),&
                     n_block_cols=size(b_col_sizes),&
                     error=error)
             end if
          else
             call cp_dstruct_create(m_struct%distribution,global_env=global_env,&
                  n_block_rows=size(b_dims%block_begins_at_row)-1,&
                  n_block_cols=size(b_dims%block_begins_at_col)-1,&
                  error=error)
          end if
       end if

       if (present(b_dims)) then
          m_struct%b_dims => b_dims
          call cp_retain(b_dims,error)
       else
          CPPrecondition(present(b_col_sizes),cp_failure_level,routineP,error,failure)
          CPPrecondition(present(b_row_sizes),cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             call cp_b_dims_create(m_struct%b_dims,&
                  b_row_sizes=b_row_sizes,&
                  b_col_sizes=b_col_sizes,&
                  error=error)
          end if
       end if
    end if
    if (.not.failure) then
       if (no_spars) then
          if (present(sparsity)) then
             !CPPreconditionNoFail(.not.associated(sparsity),cp_warning_level,routineP,error)
          end if
          nullify(m_struct%sparsity)
       else if (present(sparsity)) then
          !if (associated(sparsity)) then
             allocate(m_struct%sparsity(&
                  size(m_struct%b_dims%block_begins_at_row)-1,&
                  size(m_struct%b_dims%block_begins_at_col)-1),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             CPPostcondition(size(m_struct%sparsity,1)==size(sparsity,1),cp_failure_level,routineP,error,failure)
             CPPostcondition(size(m_struct%sparsity,2)==size(sparsity,2),cp_failure_level,routineP,error,failure)
             if (.not.failure) m_struct%sparsity=sparsity
          !end if
       else
          allocate(m_struct%sparsity(&
               size(m_struct%b_dims%block_begins_at_row)-1,&
               size(m_struct%b_dims%block_begins_at_col)-1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          if (.not.failure) then
             sparsity=.false.
          end if
       end if
    end if
    m_struct%initialized=.not.failure
    CPPostcondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
  end subroutine cp_bmstruct_create
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
  subroutine cp_bmstruct_dealloc(m_struct,error)
    type(cp_b_matrix_struct_type), pointer :: m_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='cp_bmstruct_dealloc',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (associated(m_struct)) then
       CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPreconditionNoFail(m_struct%ref_count==0,cp_warning_level,routineP,error)
          m_struct%initialized=.false.
          if (associated(m_struct%b_dims)) then
             call cp_release(m_struct%b_dims, error)
          end if
          if (associated(m_struct%distribution)) then
             call cp_release(m_struct%distribution,error)
          end if
          if (associated(m_struct%sparsity)) then
             deallocate(m_struct%sparsity, stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          end if
       end if
    end if
  end subroutine cp_bmstruct_dealloc
  !***************************************************************************

!!****f* cp_b_matrix_structure/cp_bmstruct_get [1.0] *
!!
!!   NAME
!!     cp_bmstruct_get
!!
!!   SYNOPSIS
!!     Subroutine cp_bmstruct_get(m_struct, global_env, sparsity,&
!!         sparsity_frozen, symmetric, b_dims, distribution, id_nr,&
!!         ref_count, error)
!!       Type(cp_b_matrix_struct_type), Intent (IN):: m_struct
!!       Type(global_environment_type), Intent (IN), Target:: global_env
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
subroutine cp_bmstruct_get(m_struct, global_env, sparsity, sparsity_frozen,&
     symmetric,&
     b_dims, distribution, id_nr, ref_count, error)
  type(cp_b_matrix_struct_type), intent(in) :: m_struct
  type(global_environment_type), intent(in), target :: global_env
  type(cp_b_dims_type), pointer, optional :: b_dims
  type(cp_b_distribution_type), pointer, optional :: distribution
  logical, dimension(:,:), pointer, optional :: sparsity
  integer, intent(out), optional:: ref_count, id_nr
  logical, intent(out), optional:: sparsity_frozen, symmetric
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_bmstruct_get',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(sparsity)) sparsity => m_struct%sparsity
     if (present(sparsity_frozen)) &
          sparsity_frozen = cp_bmstruct_get_sparsity_frozen(m_struct,error)
     if (present(symmetric)) symmetric=m_struct%symmetric
     if (present(b_dims)) b_dims => m_struct%b_dims
     if (present(distribution)) distribution => m_struct%distribution
     if (present(ref_count)) ref_count=m_struct%ref_count
     if (present(id_nr)) id_nr=m_struct%id_nr
  end if
end subroutine cp_bmstruct_get
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
  subroutine cp_bmstruct_freeze_sparsity(matrix_struct,no_sync,error)
    type(cp_b_matrix_struct_type), intent(inout) :: matrix_struct
    logical, optional, intent(in) :: no_sync
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, sync
    integer :: handle, stat
    integer, allocatable, dimension(:,:) :: tmp_distr
    character(len=*), parameter :: routineN='cp_bmstruct_freeze_sparsity',&
         routineP=moduleN//':'//routineN
    failure=.false.
    sync=.true.

    call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(matrix_struct%initialized,cp_failure_level,routineP,error,failure)
    if ((.not. failure) .and. (.not. matrix_struct%sparsity_frozen)) then
       if (present(no_sync)) sync=.not.no_sync
       if (sync) then
          where (matrix_struct%distribution%block_owner /= &
               matrix_struct%distribution%global_env%mepos) matrix_struct%sparsity = .false.
          ! would be nice to do this directly with logical (or at least use
          ! the bits of the integer
          allocate (tmp_distr(size(matrix_struct%sparsity,1),&
               size(matrix_struct%sparsity,2)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          where (matrix_struct%sparsity)
             tmp_distr = 0
          elsewhere
             tmp_distr = 1
          end where
          call mp_sum(tmp_distr,matrix_struct%distribution%global_env%group)
          where (tmp_distr == 0)
             matrix_struct%sparsity = .false.
          elsewhere
             matrix_struct%sparsity = .true.
          end where
          deallocate(tmp_distr,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       end if
       matrix_struct%sparsity_frozen=.true.
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_bmstruct_freeze_sparsity
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
  function cp_bmstruct_get_sparsity_frozen(matrix_struct,error) result(res)
    type(cp_b_matrix_struct_type), intent(in) :: matrix_struct
    type(cp_error_type), optional, intent(inout) :: error
    logical :: res

    logical :: failure
    integer :: handle
    character(len=*), parameter :: routineN='cp_bmstruct_get_sparsity_frozen',&
         routineP=moduleN//':'//routineN
    failure=.false.

    call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(matrix_struct%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       CPPrecondition(.false.,cp_warning_level,routineP,error,failure)
       ! to do
    end if
    res=.false.
    call timestop(0.0_wp,handle)
  end function cp_bmstruct_get_sparsity_frozen
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
  subroutine cp_bmstruct_retain(m_struct,error)
    type(cp_b_matrix_struct_type), intent(inout) :: m_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_bmstruct_retain',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(m_struct%ref_count>0,cp_warning_level,routineP,error)
    if (.not. failure) then
       m_struct%ref_count=m_struct%ref_count+1
    end if
  end subroutine cp_bmstruct_retain
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
  subroutine cp_bmstruct_release(m_struct,error)
    type(cp_b_matrix_struct_type), pointer :: m_struct
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_bmstruct_release',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (associated(m_struct)) then
       if (cp_debug) then
          CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
          CPPrecondition(m_struct%ref_count>0,cp_failure_level,routineP,error,failure)
       end if
       if (.not.failure) then
          m_struct%ref_count=m_struct%ref_count-1
          if (m_struct%ref_count<1) then
             call cp_dealloc(m_struct,error)
          end if
       end if
    end if
  end subroutine cp_bmstruct_release
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
  subroutine cp_bmstruct_get_block_info(m_struct,block_row, block_col,&
       owner,n_local_rows,n_local_cols,begins_at_row,ends_at_row,&
       begins_at_col,ends_at_col,out_of_sparsity, me_owner, error)
    type(cp_b_matrix_struct_type), intent(in) :: m_struct
    integer, intent(in) :: block_col, block_row
    integer, intent(out), optional :: owner, n_local_rows,n_local_cols,&
         begins_at_row,begins_at_col, ends_at_row, ends_at_col
    logical, intent(out), optional :: me_owner, out_of_sparsity
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_bmstruct_get_block_info',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
    if (cp_debug) then
       CPPrecondition(block_row>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(block_col>0,cp_failure_level,routineP,error,failure)
       ! check also the upper bound?
       CPPrecondition(associated(m_struct%b_dims),cp_failure_level,routineP,error,failure)
       CPPrecondition(associated(m_struct%distribution),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(associated(m_struct%b_dims%block_begins_at_col),cp_failure_level,routineP,error,failure)
          CPPrecondition(associated(m_struct%b_dims%block_begins_at_row),cp_failure_level,routineP,error,failure)
          CPPrecondition(associated(m_struct%distribution%block_owner),cp_failure_level,routineP,error,failure)
       end if
    end if
    if (.not. failure) then
       if (present(owner)) owner=m_struct%distribution%block_owner(block_row,block_col)
       if (present(n_local_rows)) then
          n_local_rows=m_struct%b_dims%block_begins_at_row(block_row+1)-&
               m_struct%b_dims%block_begins_at_row(block_row)
       end if
       if (present(n_local_cols)) then
          n_local_cols=m_struct%b_dims%block_begins_at_col(block_col+1)-&
               m_struct%b_dims%block_begins_at_col(block_col)
       end if
       if (present(begins_at_row)) begins_at_row=m_struct%b_dims%block_begins_at_row(block_row)
       if (present(begins_at_col)) begins_at_col=m_struct%b_dims%block_begins_at_col(block_col)
       if (present(ends_at_row)) ends_at_row=m_struct%b_dims%block_begins_at_row(block_row+1)-1
       if (present(ends_at_col)) ends_at_col=m_struct%b_dims%block_begins_at_col(block_col+1)-1
       if (present(out_of_sparsity)) then
          if (.not.m_struct%sparsity_frozen .or. &
               .not.associated(m_struct%sparsity)) then
             out_of_sparsity=.false.
          else
             out_of_sparsity=m_struct%sparsity(block_row,block_col)
          end if
       end if
    end if
  end subroutine cp_bmstruct_get_block_info
  !***************************************************************************

  !====== block iterator ======

!!****f* cp_bs_iter_init [1.0] *
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
  subroutine cp_bs_iter_init(bs_iterator, m_struct, do_global,&
       do_zero, do_lower_half, error)
    type(cp_bs_iterator_type), intent(out) :: bs_iterator
    type(cp_b_matrix_struct_type), intent(in), target :: m_struct
    logical, intent(in) :: do_global
    logical, intent(in), optional :: do_zero, do_lower_half
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_bs_iter_init',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(m_struct%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       bs_iterator%do_global=do_global
       bs_iterator%do_zero=.false.
       bs_iterator%do_lower_half=.not.m_struct%symmetric
       if (present(do_zero)) bs_iterator%do_zero=do_zero
       if (present(do_lower_half)) bs_iterator%do_lower_half=do_lower_half
       bs_iterator%block_row=1
       bs_iterator%block_col=1
       bs_iterator%m_struct => m_struct
       bs_iterator%n_block_rows=size(m_struct%b_dims%block_begins_at_row)
       bs_iterator%n_block_cols=size(m_struct%b_dims%block_begins_at_col)
       bs_iterator%initialized=.true.
    end if
  end subroutine cp_bs_iter_init
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
  subroutine cp_bs_iter_dealloc_ref(bs_iterator,error)
    type(cp_bs_iterator_type), intent(inout) :: bs_iterator
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_bs_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(bs_iterator%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       nullify(bs_iterator%m_struct)
       bs_iterator%initialized=.false.
    end if
  end subroutine cp_bs_iter_dealloc_ref
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
  function cp_bs_iter_next(bs_iterator,block_row, block_col,&
       owner,n_local_rows,n_local_cols,begins_at_row,ends_at_row,&
       begins_at_col,ends_at_col,out_of_sparsity, me_owner, error)&
       result(res)
    type(cp_bs_iterator_type), intent(inout) :: bs_iterator
    integer, intent(out), optional :: block_col, block_row,owner,&
         n_local_rows,n_local_cols,begins_at_row,begins_at_col,&
         ends_at_row, ends_at_col
    logical, intent(out), optional :: me_owner, out_of_sparsity
    type(cp_error_type), optional, intent(inout) :: error
    logical :: res

    logical :: failure
    character(len=*), parameter :: routineN='cp_bs_iter_next',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug) then
       CPPrecondition(bs_iterator%initialized,cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       bs_iterator%block_col=bs_iterator%block_col+1
       if (bs_iterator%block_col> bs_iterator%n_block_cols) then
          bs_iterator%block_row=bs_iterator%block_row+1
          if (bs_iterator%block_row<1.or.bs_iterator%block_row> bs_iterator%n_block_rows) then
             bs_iterator%block_row=-1
             bs_iterator%block_col= bs_iterator%n_block_cols+1
             if (present(block_row)) block_row=-1
             if (present(block_col)) block_col=-1
             res=.false.
          else
             if (bs_iterator%do_lower_half) then
                bs_iterator%block_col=1
             else
                bs_iterator%block_col=bs_iterator%block_row
             end if
          end if
       end if
       if (present(block_row)) block_row=bs_iterator%block_row
       if (present(block_col)) block_col=bs_iterator%block_col
       call cp_bmstruct_get_block_info(bs_iterator%m_struct,&
            bs_iterator%block_row, bs_iterator%block_col,&
            owner,n_local_rows,n_local_cols,begins_at_row,ends_at_row,&
            begins_at_col,ends_at_col,out_of_sparsity, me_owner, error)
    end if
  end function cp_bs_iter_next
  !***************************************************************************

end module cp_b_matrix_structure
