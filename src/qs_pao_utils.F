!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "qs_prep_globals.h"

!!****h* cp2k/qs_pao_utils [1.0] *
!!
!!   NAME
!!     qs_pao_utils
!!
!!   FUNCTION
!!     all kinds of utilities
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
module qs_pao_utils
  use qs_log_handling
  use qs_error_handling
  use kinds, only: wp => dp
  use global_types, only: global_environment_type
  use matrix_types
  implicit none
  private
  logical, private, parameter :: debug_this_module=.true.

! =========== public objects ==========
! types
  public :: qs_1d_i_array, qs_matrix_block_iterator

! procedures
  public :: qs_guarantee_size, qs_set_to_identity, qs_set_to, &
       qs_get_random_integer, qs_p_norm_matrix, qs_sup_norm_matrix,&
       qs_set_to_identity

! core methods
  public :: qs_init, qs_dealloc_ref, qs_set, qs_get, qs_next, qs_valid

! ========== util types =========
  type qs_1d_i_array
     integer, dimension(:), pointer :: array
  end type qs_1d_i_array

! ============== interfaces ===========

!! if the pointer does not point to an array of the requested size 
!! allocates it anew
  interface qs_guarantee_size
     module procedure qs_guarantee_size_1d_int, qs_guarantee_size_1d_real, &
          qs_guarantee_size_1d_logical
  end interface

!! set the matrix to the identity matrix
  interface qs_set_to_identity
     module procedure qs_set_real_matrix_to_identity
  end interface

!! set all the values to a given value
  interface qs_set_to
     module procedure qs_set_real_matrix_to
  end interface

! = core methods =

!! initializes the object
  interface qs_init
     module procedure qs_init_matrix_b_i
  end interface

!! deallocates the memory allocated by the object
  interface qs_dealloc_ref
     module procedure qs_dealloc_matrix_b_i_ref
  end interface

!! gets various attributes from the object
  interface qs_get
     module procedure qs_matrix_b_i_get
  end interface

!! goes to the next element, returns true if valid
  interface qs_next
     module procedure qs_next_matrix_block
  end interface

!! return true if it is valid (minimal testing)
  interface qs_valid
     module procedure qs_matrix_b_i_valid
  end interface
!!***
!****************************************************************************

!!****s* qs_pao_utils/qs_matrix_block_iterator [1.0] *
!!
!!   NAME
!!     qs_matrix_block_iterator
!!
!!   FUNCTION
!!     goes through all the blocks of the given matrix
!!
!!   NOTES
!!     should be in matrix_types
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
  type qs_matrix_block_iterator
     private
     logical :: initialized
     type(real_matrix_type), pointer :: matrix
     integer :: block_row
     type(real_block_node_type), pointer :: block_node
  end type qs_matrix_block_iterator
!!***
!****************************************************************************

contains

!!****f* qs_pao_utils/qs_get_random_integer [1.0] *
!!
!!   NAME
!!     qs_get_random_integer
!!
!!   FUNCTION
!!     This function returns a random integer in the requested range
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     from: minumum random number defaults to 1
!!     to: maximum random number defaults to 100
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling, but does NOT report
!!            the error in this variable (to make it pure)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!****************************************************************************
  pure function qs_get_random_integer(from,to,error) result(res)
    integer :: res
    integer, optional, intent(in) :: from, to
    type(qs_error_type), optional, intent(in) :: error

    logical :: failure
    character(len=*), parameter :: routineN='util:qs_get_random_integer'
    real :: rand
    integer :: fromi, toi
    type(qs_error_type) :: iError
    failure=.FALSE.
!call qs_init(iError,template_error=error)
    if (.not. present(from)) then
       fromi=1
    else
       fromi=from
    end if
    if (.not. present(to)) then
       toi=100
    else
       toi=to
    end if
!QSPrecondition(fromi<=toi,qs_failure_level,routineN,iError,failure)
    if (.not.fromi<=toi) failure=.true.
    if (.NOT. failure) then
       call random_number(rand)
       res=fromi+floor(rand*real(toi-fromi+1))
!QSPostcondition(res<=toi.and.res>=fromi,qs_warning_level,routineN,iError,failure)
       if (res>toi) res=toi
       if (res<fromi) res=fromi
    else
       res=-1
    end if
!call qs_dealloc_ref(iError,error=error)
  end function qs_get_random_integer
!!***
!****************************************************************************

! ============== real full matrix utilities ============

!!****f* qs_pao_utils/qs_p_norm_matrix [1.0] *
!!
!!   NAME
!!     qs_p_norm_matrix
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  pure function qs_p_norm_matrix(matrix, p) result(res)
    real(kind=wp) :: res
    real(kind=wp), intent(in), dimension(:,:) :: matrix
    integer,optional :: p
!type(qs_error_type), optional, intent(inout) :: error ! removed to make it pure

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_froebnius_norm2'
    integer :: i,j
    failure=.FALSE.

    res=0
    if (.NOT. failure) then
       if (.not. present(p)) p=2
       do i=1,size(matrix,1)
          do j=1,size(matrix,2)
             res=res+matrix(i,j)**p
          end do
       end do
    end if
  end function qs_p_norm_matrix
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_sup_norm_matrix [1.0] *
!!
!!   NAME
!!     qs_sup_norm_matrix
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
!!            see module qs_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!     @version 12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  pure function qs_sup_norm_matrix(matrix) result(res)
    real(kind=wp) :: res
    real(kind=wp), dimension(:,:), intent(in) :: matrix
!type(qs_error_type), optional, intent(inout) :: error ! removed to make it pure

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_sup_norm_matrix'
    failure=.FALSE.

    if (.NOT. failure) then
       res=maxval(abs(matrix))
    else
       res=-1.0_wp
    end if
  end function qs_sup_norm_matrix
!!***
!***************************************************************************

! ========== memory utilities =========

!!****f* qs_pao_utils/qs_guarantee_size_1d_int [1.0] *
!!
!!   NAME
!!     qs_guarantee_size_1d_int
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     Issues a warning when the size changes (but not on allocation
!!     and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_guarantee_size_1d_int(array, n, init, error)
    integer, dimension(:), pointer :: array
    integer, intent(in) :: n
    integer, intent(in), optional :: init
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='qs_pao_utils:qs_guarantee_size_1d_int'
    failure=.false.

    QSPrecondition(n>=0,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             QSErrorMessage(qs_warning_level,routineN,'size has changed',error)
             deallocate(array, stat=stat)
             QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          QSPostcondition(stat==0, qs_failure_level,routineN,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine qs_guarantee_size_1d_int
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_guarantee_size_1d_logical [1.0] *
!!
!!   NAME
!!     qs_guarantee_size_1d_logical
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     Issues a warning when the size changes (but not on allocation
!!     and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_guarantee_size_1d_logical(array, n, init, error)
    logical, dimension(:), pointer :: array
    integer, intent(in) :: n
    logical, intent(in), optional :: init
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='qs_pao_utils:qs_guarantee_size_1d_logical'
    failure=.false.

    QSPrecondition(n>=0,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             QSErrorMessage(qs_warning_level,routineN,'size has changed',error)
             deallocate(array, stat=stat)
             QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          QSPostcondition(stat==0, qs_failure_level,routineN,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine qs_guarantee_size_1d_logical
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_guarantee_size_1d_real [1.0] *
!!
!!   NAME
!!     qs_guarantee_size_1d_real
!!
!!   FUNCTION
!!     If the size of the array is changes reallocate it.
!!     If warning_on_change is present and true issues a warning 
!!     whe the size changes (but not on allocation and deallocation).
!!     
!!     The data is NOT preserved (if you want to preserve the data see
!!     the realloc in the module memory_utilities)
!!
!!   NOTES
!!     this is a different behaviour than the realloc in the module
!!     memory_utilities. It is quite low level
!!
!!   INPUTS
!!     array: the array to reallocate if necessary
!!     n: the wanted size
!!     init: the initial value of the elements when resized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_guarantee_size_1d_real(array, n, init, error)
    real(kind=wp), dimension(:), pointer :: array
    integer, intent(in) :: n
    real(kind=wp), optional, intent(in) :: init
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: stat
    character(len=*), parameter :: routineN='qs_pao_utils:qs_guarantee_size_1d_real'
    failure=.false.

    QSPrecondition(n>=0,qs_failure_level,routineN,error,failure)
    failureIf: if (.not. failure) then
       if (associated(array)) then
          if (size(array) /= n) then
             QSErrorMessage(qs_warning_level,routineN,'size has changed',error)
             deallocate(array, stat=stat)
             QSPostcondition(stat==0,qs_warning_level,routineN,error,failure)
          end if
       end if
       if (.not.associated(array)) then
          allocate(array(n), stat=stat)
          QSPostcondition(stat==0, qs_failure_level,routineN,error,failure)
          if (present(init).and..not.failure) array=init
       end if
    end if failureIf
  end subroutine qs_guarantee_size_1d_real
!!***
!***************************************************************************

! =========== real block matrix utilities =========

!!****f* qs_pao_utils/qs_set_real_matrix_to_identity [1.0] *
!!
!!   NAME
!!     qs_set_real_matrix_to_identity
!!
!!   FUNCTION
!!     initialized the real matrix to the identity matrix
!!
!!   NOTES
!!     does not remove off diagonal blocks, should be in matrix_types, 
!!     assumes a square matrix form.
!!
!!   INPUTS
!!     matrix: the matrix to be initialized
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_set_real_matrix_to_identity(matrix,error)
    type(real_matrix_type), target, intent(inout) :: matrix
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure, did_ii
    character(len=*), parameter :: routineN='qs_pao_utils:qs_set_real_matrix_to_identity'
    type(real_matrix_type), pointer :: m_ptr
    integer :: block_row_i, block_col_i, i, nblock_row, nblock_col
    type(real_block_node_type), pointer :: blockAtt
    real(kind=wp), dimension(:,:), pointer :: blockVal
    failure=.false.
    nullify(blockAtt)
    m_ptr => matrix

    QSPrecondition(associated(m_ptr),qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       call get_matrix_info(m_ptr,nblock_col=nblock_col, nblock_row=nblock_row)
       do block_row_i=1,nblock_row
          blockAtt => first_block_node(m_ptr,block_row_i)
          did_ii=.false.
          do
             if (.not.associated(blockAtt)) exit
             call get_block_node(blockAtt, block_col=block_col_i,block=blockVal)
             blockVal=0.0_wp
             if (block_col_i==block_row_i) then
                QSPrecondition(size(blockVal,1)==size(blockVal,2),qs_failure_level,routineN,error,failure)
                forall (i=1:size(blockVal,1)) blockVal(i,i)=1.0_wp
                did_ii=.true.
             end if
             blockAtt => next_block_node(blockAtt)
          end do
          if (.not. did_ii) then
             nullify(blockVal)
             call add_block_node(m_ptr,block_row_i,block_row_i,block=blockVal)
             blockVal=0.0_wp
             QSPrecondition(size(blockVal,1)==size(blockVal,2),qs_failure_level,routineN,error,failure)
             forall (i=1:size(blockVal,1)) blockVal(i,i)=1.0_wp
          end if
       end do
    end if
  end subroutine qs_set_real_matrix_to_identity
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_set_real_matrix_to [1.0] *
!!
!!   NAME
!!     qs_set_real_matrix_to
!!
!!   FUNCTION
!!     initializes the allocated blocks of the real matrix with the given
!!     value (defaults to 0)
!!
!!   NOTES
!!     the iteration looks ugly because this function has much error
!!     handling to test the iterator (it should be:
!!     do; if (.not.qs_next(iterator)) exit
!!       call qs_get(..);..
!!     end do
!!
!!   INPUTS
!!     matrix: the matrix to initialize
!!     value: the new value
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
  subroutine qs_set_real_matrix_to(matrix, value, error)
    type(real_matrix_type), intent(inout), target :: matrix
    real(kind=wp), optional, intent(in) :: value
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_set_real_matrix_to'
    type(qs_matrix_block_iterator) :: iterator
    real(kind=wp), dimension(:,:), pointer :: block_val
    real(kind=wp) :: val
    type(qs_error_type) :: iError
    failure=.false.
    call qs_init(iError,template_error=error)
    val=0.0_wp

    if (.not. failure) then
       call qs_init(iterator,matrix,error=iError)
       call qs_propagate_error(iError,routineN,QSSourceFileRef(),error,failure)
    end if
    if (.not.failure) then
       if (present(value)) val=value
       do
          if (qs_next(iterator,error=iError)) then
             call qs_get(iterator,block_val=block_val,error=iError)
             call qs_propagate_error(iError,routineN,QSSourceFileRef(),error,failure)
             QSPrecondition(associated(block_val),qs_failure_level,routineN,error,failure)
             if (.not.failure) then
                block_val=val
             else
                exit
             end if
          else
             call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
                  error,failure)
             exit
          end if
       end do
       call qs_dealloc_ref(iterator,error=error)
    end if
    call qs_dealloc_ref(iError,error=error)
  end subroutine qs_set_real_matrix_to
!!***
!***************************************************************************

! ================== qs_matrix_block_iterator ================

!!****f* qs_pao_utils/qs_init_matrix_b_i [1.0] *
!!
!!   NAME
!!     qs_init_matrix_b_i
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
!!            see module qs_error_handling 
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
  subroutine qs_init_matrix_b_i(block_iterator, matrix, block_row,&
       block_node, error)
    type(qs_matrix_block_iterator), intent(out) :: block_iterator
    type(real_matrix_type), target, intent(in) :: matrix
! well subsequent operations on the iterator could change tha data in the matrix
    integer, intent(in), optional :: block_row
    type(real_block_node_type), target, optional :: block_node
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_init_matrix_b_i'
    failure=.false.

    block_iterator%matrix => matrix
    block_iterator%block_row = 1
    nullify(block_iterator%block_node)
    block_iterator%initialized=.true.
    if (present(block_node)) &
         block_iterator%block_node => block_node
    if (present(block_row)) &
         block_iterator%block_row = block_row
  end subroutine qs_init_matrix_b_i
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_dealloc_matrix_b_i_ref [1.0] *
!!
!!   NAME
!!     qs_dealloc_matrix_b_i_ref
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
!!            see module qs_error_handling 
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
  subroutine qs_dealloc_matrix_b_i_ref(block_iterator, error)
    type(qs_matrix_block_iterator), intent(inout) :: block_iterator
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_dealloc_matrix_b_i_ref'
    failure=.false.

    QSPrecondition(block_iterator%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       nullify(block_iterator%matrix, block_iterator%block_node)
       block_iterator%block_row=0
    end if
  end subroutine qs_dealloc_matrix_b_i_ref
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_matrix_b_i_get [1.0] *
!!
!!   NAME
!!     qs_matrix_b_i_get
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
!!            see module qs_error_handling 
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
  subroutine qs_matrix_b_i_get(block_iterator, block_row, block_col,&
       block_val, finished, matrix, block_node, error)
    type(qs_matrix_block_iterator), intent(in) :: block_iterator
    integer, optional, intent(out) :: block_row, block_col
    real(kind=wp), dimension(:,:), pointer, optional :: block_val
    logical, optional, intent(out) :: finished
    type(real_matrix_type), optional, pointer :: matrix
    type(real_block_node_type), optional, pointer :: block_node
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_matrix_b_i_get'
    failure=.false.

    QSPrecondition(block_iterator%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(block_iterator%matrix),qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (present(block_row)) &
            block_row=block_iterator%block_row
       if (present(block_col)) then
          if (associated(block_iterator%block_node)) then
             call get_block_node(block_iterator%block_node, block_col=block_col)
          else
             QSAssert(.false.,qs_warning_level,routineN,error,failure)
             block_col=-1
          end if
       end if
       if (present(block_val)) then
          QSPrecondition(block_iterator%block_row>0,qs_failure_level,routineN,error,failure)
          if (.not.failure) then
             if (associated(block_iterator%block_node)) then
                call get_block_node(block_iterator%block_node,block=block_val)
             else
                nullify(block_val)
                QSErrorMessage(qs_failure_level,routineN,"invalid block",error)
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
  end subroutine qs_matrix_b_i_get
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_matrix_b_i_valid [1.0] *
!!
!!   NAME
!!     qs_matrix_b_i_valid
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
!!            see module qs_error_handling 
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
  function qs_matrix_b_i_valid(block_iterator,error) result(res)
    logical ::res
    type(qs_matrix_block_iterator), intent(in) :: block_iterator
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_matrix_b_i_valid'
    failure=.false.

    QSPrecondition(block_iterator%initialized,qs_failure_level,routineN,error,failure)
    if (.not. failure) then
       if (.not.(associated(block_iterator%matrix).and. block_iterator%block_row>0 ))then
          QSAssert(block_iterator%block_row > 0,qs_warning_level,routineN,error,failure)
          QSAssert(associated(block_iterator%matrix),qs_warning_level,routineN,error,failure)
       end if
    end if
    res=.not.failure
  end function qs_matrix_b_i_valid
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_next_matrix_block [1.0] *
!!
!!   NAME
!!     qs_next_matrix_block
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
!!            see module qs_error_handling 
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
  function qs_next_matrix_block(block_iterator, error) result(res)
    logical :: res
    type(qs_matrix_block_iterator), intent(inout) :: block_iterator
    type(qs_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_pao_utils:qs_next_matrix_block'
    integer :: nblock_row,i
    failure=.false.

    QSPrecondition(block_iterator%initialized,qs_failure_level,routineN,error,failure)
    QSPrecondition(associated(block_iterator%matrix),qs_failure_level,routineN,error,failure)
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
       QSPostcondition(.not.associated(block_iterator%block_node),qs_warning_level,routineN,error,failure)
       if (.not. block_iterator%block_row <= 0) then
          QSErrorMessage(qs_warning_level,routineN,"inconsistent state in iterator",error)
          block_iterator%block_row=-2
       end if
    else
       block_iterator%block_row=-2
       nullify(block_iterator%block_node)
    end if
    res=.false.
  end function qs_next_matrix_block
!!***
!***************************************************************************

!!****f* qs_pao_utils/qs_sparse_m_sym_mult [1.0] *
!!
!!   NAME
!!     qs_sparse_m_sym_mult
!!
!!   FUNCTION
!!     multiplies two sparsesymmetric matrices 
!!
!!   NOTES
!!     not parallelized, but easy to do: first build cols of a then
!!     do different rows on different processors (balancing  it well 
!!     will give more work)
!!
!!   INPUTS
!!     a: the first matrix
!!     transpose_a: if a should be transposed
!!     b: the second matrix
!!     transpose_b: if b should be transposed
!!     c: the resulting matrix (c = scale_c c + scale_ab a * b 
!!          (maybe with a^T or b^T in place of a and b))
!!     scale_ab: how a*b should be scaled (defaults to 1)
!!     scale_c: how c should be scaled (defaults to 0)
!!     error: variable to control error logging, stopping,... 
!!            see module qs_error_handling 
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
!!$subroutine qs_sparse_m_sym_mult(a,b,c,transpose_a,transpose_b,&
!!$     scale_c,scale_ab,error)
!!$  type(real_matrix_type), intent(in), target :: a,b
!!$  type(real_matrix_type), intent(in), target :: c
!!$  logical, optional, intent(in) :: transpose_b, transpose_a
!!$  real(kind=wp), optional, intent(in) :: scale_ab, scale_c
!!$  type(qs_error_type), optional, intent(inout) :: error
!!$  
!!$  logical :: failure
!!$  character(len=*), parameter :: routineN='qs_pao_utils:qs_sparse_m_sym_mult'
!!$  type(real_matrix_type), pointer :: a_m,b_m,c_m
!!$  real(kind=wp) :: s_c,s_ab
!!$  real(kind=wp), dimension(:,:), pointer :: block_val
!!$  type(qs_error_type) :: iError
!!$  type(qs_sll_sm_node_type), dimension(:), allocatable :: col_a
!!$  
!!$  type(qs_matrix_block_iterator) :: m_iter
!!$  character(len=1), dimension(2) :: trans_a,trans_b
!!$  integer :: block_row, block_col, i, j, old_block_row, ncol_a,&
!!$       stat
!!$  failure=.false.
!!$  nullify(block_val)
!!$  call qs_init(iError, template_error=error)
!!$  a_m => a; b_m => b; c_m => c
!!$  s_ab=1.0_wp; s_c=0.0_wp
!!$
!!$  if (present(scale_c)) s_c=scale_c
!!$  if (present(scale_ab)) s_ab=scale_ab
!!$
!!$  QSPrecondition(associated(a_m),qs_failure_level,routineN,error,failure)
!!$  QSPrecondition(associated(b_m),qs_failure_level,routineN,error,failure)
!!$  QSPrecondition(associated(c_m),qs_failure_level,routineN,error,failure)
!!$  if (failure) return
!!$
!!$  ! scale c
!!$  call qs_init(m_iter,matrix=c,error=iError)
!!$  call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
!!$                error,failure,qs_failure_level)
!!$  if (.not.failure) then
!!$     do
!!$        if (.not.next(m_iter,error=error)) exit
!!$        call qs_get(m_iter, block_val=block_val, error=iError)
!!$        if (s_c==0.0_wp) then
!!$           block_val=0
!!$        else
!!$           forall (i=1:size(block_val,1))
!!$              forall (j=1:size(block_val,2))
!!$                 block_val(i,j) = s_c*block_val(i,j)
!!$              end forall
!!$           end forall
!!$        end if
!!$     end do
!!$     call qs_dealloc_ref(m_iter,error=iError)
!!$     call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
!!$                error,failure,qs_failure_level)
!!$  end if
!!$  
!!$  ! set up transposition arrays
!!$  trans_a = (/ 'T','N' /)
!!$  trans_b = (/ 'T','N' /)
!!$  if (present(transpose_a)) then
!!$     if (transpose_a) then
!!$        trans_a = (/ 'N','T' /)
!!$     end if
!!$  end if
!!$  if (present(transpose_b)) then
!!$     if (transpose_b) then
!!$        trans_b = (/ 'N','T' /)
!!$     end if
!!$  end if
!!$
!!$  ! initialize col linked lists
!!$  call get_matrix_info(a_m,ncol=ncol_a)
!!$  allocate(col_a(ncol_a),stat=stat)
!!$  QSPostcondition(stat==0,qs_failure_level,routineN,error,failure)
!!$
!!$  ! iterates on a
!!$  if (.not. failure) then
!!$     call qs_init(m_iter,matrix=a,error=iError)
!!$     call qs_propagate_error(iError,routineN,QSSourceFileRef(),&
!!$          error,failure,qs_failure_level)
!!$  end if
!!$  if (.not.failure) then
!!$     old_block_row=1
!!$ !    do 
!!$ !       qs
!!$  end if
!!$end subroutine qs_sparse_m_sym_mult
!!***
!***************************************************************************


end module qs_pao_utils

