!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_projection [1.0] *
!!
!!   NAME
!!     pao_projection
!!
!!   FUNCTION
!!     represent a simple projection (simply keeps some indices)
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
module pao_projection
  use cp_log_handling
  use cp_error_handling
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_projection'

  !***
  !****************************************************************************

!!****s* pao_types/pao_projection_type [1.0] *
!!
!!   NAME
!!     pao_projection_type
!!
!!   FUNCTION
!!     represent a simple projection (ie a subset of indices)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     initialized: true if initialized (bug catcher)
!!     proj_indices: the local indices of the full basis that are kept
!!            by the projection
!!     non_proj_indices: the local indices of the full basis that are
!!            not in the min_bas (i.e. the indices of the excluded basis)
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
  type pao_projection_type
     logical :: initialized
     integer, dimension(:), pointer :: proj_indices
     integer, dimension(:), pointer :: non_proj_indices
  end type pao_projection_type
!!***
  !****************************************************************************

contains


!!****f* pao_projection/cp_proj_init [1.0] *
!!
!!   NAME
!!     cp_proj_init
!!
!!   FUNCTION
!!     initializes the given projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to be initialized
!!     proj_indices: the indices of the full basis that are kept
!!     size_full_bas: the size of the full basis
!!     non_proj_indices: the index of the excluded basis
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
  subroutine cp_proj_init(proj, proj_indices, size_full_bas, &
       error)
    type(pao_projection_type), intent(out) :: proj
    integer, dimension(:), intent(in) ::proj_indices
    integer, intent(in) :: size_full_bas
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure,found
    integer :: handle, stat, i, j, exc_pos
    character(len=*), parameter :: routineN='cp_proj_init',&
         routineP=moduleN//':'//routineN
    failure=.false.

    !call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(all(proj_indices>0),cp_failure_level,routineP,error,failure)
    CPPrecondition(size(proj_indices)<=size_full_bas,cp_failure_level,routineP,error,failure)
    allocate(proj%proj_indices(size(proj_indices)),&
         proj%non_proj_indices(size_full_bas-size(proj_indices)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !do i=2,size(proj_indices)
    !   CPPrecondition(proj_indices(i-1)<proj_indices(i),cp_failure_level,routineP,error,failure)
    !end do
    if (.not.failure) then
       proj%proj_indices=proj_indices
       exc_pos=0
       do i=1,size_full_bas
          found=.false.
          do j=1,size(proj%proj_indices)
             if (i==proj%proj_indices(j)) then
                found=.true.
                exit
             end if
          end do
          if (.not.found) then
             CPInvariant(exc_pos<size(proj%non_proj_indices),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                exc_pos=exc_pos+1
                proj%non_proj_indices(exc_pos)=i
             end if
          end if
       end do
       proj%initialized=.true.
    end if
    !call timestop(0.0_wp,handle)
  end subroutine cp_proj_init
  !***************************************************************************

!!****f* pao_projection/cp_proj_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_proj_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory allocated by this object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to be deallocated
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
  subroutine cp_proj_dealloc_ref(proj, error)
    type(pao_projection_type), intent(inout) :: proj
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, stat
    character(len=*), parameter :: routineN='cp_proj_dealloc_ref',&
         routineP=moduleN//':'//routineN
    failure=.false.

    !call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       deallocate(proj%proj_indices, proj%non_proj_indices, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       proj%initialized=.false.
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_proj_dealloc_ref
  !***************************************************************************

!!****f* pao_projection/cp_proj_valid [1.0] *
!!
!!   NAME
!!     cp_proj_valid
!!
!!   FUNCTION
!!     returns true if the projection is valid (minimal tests)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to test
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
  function cp_proj_valid(proj,error) result(res)
    type(pao_projection_type), intent(in) :: proj
    logical ::res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_proj_valid',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(proj%proj_indices),cp_warning_level,routineP,error,failure)
    CPPrecondition(associated(proj%non_proj_indices),cp_warning_level,routineP,error,failure)
    res=.not.failure
  end function cp_proj_valid
  !***************************************************************************

!!****f* pao_projection/cp_proj_get [1.0] *
!!
!!   NAME
!!     cp_proj_get
!!
!!   FUNCTION
!!     returns various attributes about the given projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection object
!!     size_min_bas: the size of the minimal basis
!!     size_full_bas: the size of the full basis
!!     size_excl_bas: the size of the excluded basis
!!     proj_indices: the local indices of the full basis that are kept
!!            by the projection
!!     non_proj_indices: the local indices of the full basis that are
!!            not in the min_bas (i.e. the indices of the excluded basis)
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
  subroutine cp_proj_get(proj, size_min_bas, size_full_bas, &
       size_excl_bas, proj_indices, non_proj_indices, error)
    type(pao_projection_type), optional, intent(in) :: proj
    integer, optional, intent(out) :: size_excl_bas,size_full_bas,&
         size_min_bas
    integer, optional, dimension(:), pointer :: proj_indices,&
         non_proj_indices
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle
    character(len=*), parameter :: routineN='',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    !call timeset(routineN//','//moduleN,'I',"",handle)
    if (.not. failure) then
       if (present(size_min_bas)) size_min_bas=size(proj%proj_indices)
       if (present(size_full_bas)) &
            size_full_bas=size(proj%proj_indices)+size(proj%non_proj_indices)
       if (present(size_excl_bas)) size_excl_bas=size(proj%non_proj_indices)
       if (present(proj_indices)) proj_indices => proj%proj_indices
       if (present(non_proj_indices)) non_proj_indices => proj%non_proj_indices
    end if
    !call timestop(0.0_wp,handle)
  end subroutine cp_proj_get
  !***************************************************************************

!!****f* pao_projection/cp_proj_get_proj_index [1.0] *
!!
!!   NAME
!!     cp_proj_get_proj_index
!!
!!   FUNCTION
!!     returns the projection index
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to get the info from
!!     i: the local index in the min bas
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
  function cp_proj_get_proj_index(proj, i, error) result(res)
    type(pao_projection_type), intent(in) :: proj
    integer, intent(in):: i
    integer :: res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle
    character(len=*), parameter :: routineN='cp_proj_get_proj_index',&
         routineP=moduleN//':'//routineN
    failure=.false.

    !call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_debug.and.debug_this_module) then
       CPPrecondition(cp_proj_valid(proj,error=error),cp_failure_level,routineP,error,failure)
       CPPrecondition(i>0,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(i<=size(proj%proj_indices),cp_failure_level,routineP,error,failure)
       end if
    end if
    if (.not. failure) then
       res=proj%proj_indices(i)
    else
       res=-1
    end if
    !call timestop(0.0_wp,handle)
  end function cp_proj_get_proj_index
  !***************************************************************************

!!****f* pao_projection/cp_proj_get_proj_index [1.0] *
!!
!!   NAME
!!     cp_proj_get_proj_index
!!
!!   FUNCTION
!!     returns the projection index
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to get the info from
!!     i: the local index in the min bas
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
  function cp_proj_get_non_proj_index(proj, i, error) result(res)
    type(pao_projection_type), intent(in) :: proj
    integer, intent(in):: i
    integer :: res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle
    character(len=*), parameter :: routineN='cp_proj_get_non_proj_index',&
         routineP=moduleN//':'//routineN
    failure=.false.

    !call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_debug.and.debug_this_module) then
       CPPrecondition(cp_proj_valid(proj,error=error),cp_failure_level,routineP,error,failure)
       CPPrecondition(i>0,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(i<=size(proj%non_proj_indices),cp_failure_level,routineP,error,failure)
       end if
    end if
    if (.not. failure) then
       res=proj%non_proj_indices(i)
    else
       res=-1
    end if
    !call timestop(0.0_wp,handle)
  end function cp_proj_get_non_proj_index
  !***************************************************************************

end module pao_projection
