!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sp_iterators [1.0] *
!!
!!   NAME
!!     cp_sp_iterators
!!
!!   FUNCTION
!!     iterators on the sparse block matrix
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
module cp_sp_iterators
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_sp_iterators'

  public :: cp_sp_local_b_iter_type, cp_sp_owned_b_iter_type,&
       cp_sp_cached_b_iter_type, cp_sp_g_block_iterator
  public :: cp_init, cp_dealloc_ref, cp_next
  public :: cp_sp_local_b_iter_init, cp_sp_local_b_iter_dealloc_ref,&
       cp_sp_local_b_iter_next, cp_sp_owned_b_iter_init,&
       cp_sp_owned_b_iter_dealloc_ref, cp_sp_cached_b_iter_init,&
       cp_sp_owned_b_iter_dealloc_ref

  interface cp_init
     module procedure cp_sp_local_b_iter_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_sp_local_b_iter_dealloc_ref
  end interface
  interface cp_next
     module procedure cp_sp_local_b_iter_next
  end interface
  interface cp_init
     module procedure cp_sp_owned_b_iter_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_sp_owned_b_iter_dealloc_ref
  end interface
  interface cp_init
     module procedure cp_sp_cached_b_iter_init
  end interface
  interface cp_dealloc_ref
     module procedure cp_sp_owned_b_iter_dealloc_ref
  end interface

  !***
  !****************************************************************************
  !!****s* cp_block_matrix_type/cp_sp_local_b_iter_type [1.0] *
  !!
  !!   NAME
  !!     cp_sp_local_b_iter_type
  !!
  !!   FUNCTION
  !!     iterates on the local blocks
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
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
  !***************************************************************************
  type cp_sp_local_b_iter_type
     logical :: initialized, all_blocks
     type(cp_block_matrix_type), pointer :: matrix
     type(cp_matrix_block_iterator) :: iterator
  end type cp_sp_local_b_iter_type
  !!***
  !****************************************************************************

  !!****s* cp_block_matrix_type/cp_sp_owned_b_iter_type [1.0] *
  !!
  !!   NAME
  !!     cp_sp_owned_b_iter_type
  !!
  !!   FUNCTION
  !!     iterates on the owned blocks
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
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
  !***************************************************************************
  type cp_sp_owned_b_iter_type
     type(cp_sp_local_b_iter_type), pointer :: local_iterator
  end type cp_sp_owned_b_iter_type
  !!***
  !****************************************************************************

  !!****s* cp_block_matrix_type/cp_sp_cached_b_iter_type [1.0] *
  !!
  !!   NAME
  !!     cp_sp_cached_b_iter_type
  !!
  !!   FUNCTION
  !!     iterates on the cached blocks
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
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
  !***************************************************************************
  type cp_sp_cached_b_iter_type
     type(cp_sp_local_b_iter_type), pointer :: local_iterator
  end type cp_sp_cached_b_iter_type
  !!***
  !****************************************************************************

  !!****s* cp_block_matrix_type/cp_sp_g_block_iterator [1.0] *
  !!
  !!   NAME
  !!     cp_sp_g_block_iterator
  !!
  !!   FUNCTION
  !!     iterates on all the non zero blocks of the matrix
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
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
  !***************************************************************************
  type cp_sp_g_block_iterator
     logical :: initialized
     type(cp_block_matrix_type) :: matrix
     integer :: i,j
  end type cp_sp_g_block_iterator
  !!***
  !****************************************************************************

contains
  !=========== block iterators ============

  !==== local block iterator ====

  !!****f* cp_sp_iterators/cp_sp_local_b_iter_init [1.0] *
  !!
  !!   NAME
  !!     cp_sp_local_b_iter_init
  !!
  !!   FUNCTION
  !!     an iteratotor that goes trough all the non zero blocks local
  !!     blocks, both owned and cached.
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be initialized
  !!     matrix: the matrix to iterate on
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
  subroutine cp_sp_local_b_iter_init(iterator,matrix,error)
    type(cp_sp_owned_b_iter_type), intent(out) :: iterator
    type(cp_block_matrix_type), intent(in), target :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_local_b_iter_init',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (.not. failure) then

    end if
  end subroutine cp_sp_local_b_iter_init
  !***************************************************************************

  !!****f* cp_sp_local_b_iter_dealloc_ref [1.0] *
  !!
  !!   NAME
  !!      cp_sp_local_b_iter_dealloc_ref
  !!
  !!   FUNCTION
  !!     dealloc the memory used by the iterator
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be released
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
  subroutine cp_sp_local_b_iter_dealloc_ref(iterator, error)
    type(cp_sp_owned_b_iter_dealloc_ref), intent(inout) :: iterator
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_local_b_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (.not. failure) then

    end if
  end subroutine cp_sp_local_b_iter_dealloc_ref
  !***************************************************************************

  !!****f* cp_sp_iterators/cp_sp_local_b_iter_next [1.0] *
  !!
  !!   NAME
  !!     cp_sp_local_b_iter_next
  !!
  !!   FUNCTION
  !!     goes to the next local block, returns true if not at end
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator that sould move to the next block
  !!     block_handler: the handler to the actual block (to do some changes 
  !!            to the values you MUST release the handler with
  !!            call cp_release_block (block_handler,changed_values=.true.)
  !!     block_row: the actual block row
  !!     block_col: the actual block col
  !!     block_val: the actual block value (becomes invalid after
  !!                another call to next)
  !!     trans_to_do: the transposition still to do [remove ?]
  !!     local_scratch: if it is a local scratch block
  !!     owned: if it is an owned block
  !!     local_copy: if it is a local copy
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
  function cp_sp_local_b_iter_next(iterator, block_row, block_col, block_val,&
       trans_to_do, local_scatch, owned, local_copy,error)result(res)
    type(cp_sp_local_b_iter_next), intent(inout) :: iterator
    integer, intent(out), optional :: block_row, block_col
    real(kind=wp), dimension(:,:), pointer, optional :: block_val
    character(len=1), intent(out), optional :: trans_to_do
    logical, intent(out), optional :: local_scratch, owned, local_copy
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_sp_local_b_iter_next',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (.not. failure) then

    end if
  end function cp_sp_local_b_iter_next
  !***************************************************************************

  !=== owned blocks iterator ===

  !!****f* cp_sp_iterators/cp_sp_owned_b_iter_init [1.0] *
  !!
  !!   NAME
  !!     cp_sp_owned_b_iter_init
  !!
  !!   FUNCTION
  !!     an iteratotor that goes trough the non zero owned blocks.
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be initialized
  !!     matrix: the matrix to iterate on
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
  subroutine cp_sp_owned_b_iter_init(iterator,matrix,error)
    type(cp_sp_owned_b_iter_type), intent(out) :: iterator
    type(cp_block_matrix_type), intent(in), target :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineN='cp_sp_owned_b_iter_init',&
         routineP=moduleN//':'//routineN

    call cp_sp_local_b_iter_init(iterator%local_iterator,&
         matrix=matrix,error=error)
  end subroutine cp_sp_owned_b_iter_init
  !***************************************************************************

  !!****f* cp_sp_owned_b_iter_dealloc_ref [1.0] *
  !!
  !!   NAME
  !!      cp_sp_owned_b_iter_dealloc_ref
  !!
  !!   FUNCTION
  !!     dealloc the memory used by the iterator
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be released
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
  subroutine cp_sp_owned_b_iter_dealloc_ref(iterator, error)
    type(cp_sp_owned_b_iter_dealloc_ref), intent(inout) :: iterator
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineN='cp_sp_owned_b_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN

    call cp_sp_local_b_iter_dealloc_ref(iterator%local_iterator,&
         error=error)
  end subroutine cp_sp_owned_b_iter_dealloc_ref
  !***************************************************************************


  !=== cached blocks iterator ===

  !!****f* cp_sp_iterators/cp_sp_cached_b_iter_init [1.0] *
  !!
  !!   NAME
  !!     cp_sp_cached_b_iter_init
  !!
  !!   FUNCTION
  !!     an iteratotor that goes trough the non zero cached blocks.
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be initialized
  !!     matrix: the matrix to iterate on
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
  subroutine cp_sp_cached_b_iter_init(iterator,matrix,error)
    type(cp_sp_cached_b_iter_type), intent(out) :: iterator
    type(cp_block_matrix_type), intent(in), target :: matrix
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineN='cp_sp_cached_b_iter_init',&
         routineP=moduleN//':'//routineN

    call cp_sp_local_b_iter_init(iterator%local_iterator,&
         matrix=matrix,error=error)
  end subroutine cp_sp_cached_b_iter_init
  !***************************************************************************


  !!****f* cp_sp_owned_b_iter_dealloc_ref [1.0] *
  !!
  !!   NAME
  !!      cp_sp_owned_b_iter_dealloc_ref
  !!
  !!   FUNCTION
  !!     dealloc the memory used by the iterator
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be released
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
  subroutine cp_sp_owned_b_iter_dealloc_ref(iterator, error)
    type(cp_sp_owned_b_iter_dealloc_ref), intent(inout) :: iterator
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineN='cp_sp_owned_b_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN

    call cp_sp_local_b_iter_dealloc_ref(iterator%local_iterator,&
         error=error)
  end subroutine cp_sp_owned_b_iter_dealloc_ref
  !***************************************************************************

end module cp_sp_iterators
