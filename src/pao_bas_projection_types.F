!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pao_bas_projection_types [1.0] *
!!
!!   NAME
!!     pao_bas_projection_types
!!
!!   FUNCTION
!!     represent a simple projection (simply keeps some indexes)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 3.2002 created [fawzi]
!!     - 4.2002 added create, dealloc, write
!!
!!   SOURCE
!****************************************************************************
module pao_bas_projection_types
  use cp_log_handling, only: cp_failure_level, cp_warning_level, cp_to_string
  use cp_error_handling, only: cp_error_type, cp_assert, cp_debug
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type
  use timings, only: timeset, timestop
  use cp_array_utils, only: cp_write
  use cp_linked_list_int, only: cp_create, cp_insert, cp_dealloc,&
       cp_sll_int_type, cp_next
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='pao_bas_projection'
  
  public :: pao_bas_projection_type
  public :: cp_create, cp_dealloc, cp_init, cp_dealloc_ref,&
       cp_valid, cp_get, cp_get_proj_indexes, cp_get_non_proj_indexes,&
       cp_write
  public :: cp_proj_create, cp_proj_dealloc, cp_proj_init,&
       cp_proj_dealloc_ref, cp_proj_valid, cp_proj_get,&
       cp_proj_get_proj_indexes, cp_proj_get_non_proj_indexes,&
       cp_proj_write , cp_proj_init2, cp_proj_create2

  interface cp_create
     module procedure cp_proj_create , cp_proj_create2
  end interface
  interface cp_dealloc
     module procedure cp_proj_dealloc
  end interface
  interface  cp_init
     module procedure cp_proj_init , cp_proj_init2
  end interface
  interface cp_dealloc_ref
     module procedure cp_proj_dealloc_ref
  end interface
  interface cp_valid
     module procedure cp_proj_valid
  end interface
  interface cp_get
     module procedure cp_proj_get
  end interface
  interface cp_get_proj_indexes
     module procedure cp_proj_get_proj_indexes
  end interface
  interface cp_get_non_proj_indexes
     module procedure cp_proj_get_non_proj_indexes
  end interface
  interface cp_write
     module procedure cp_proj_write
  end interface

!!***
!****************************************************************************

!!****s* pao_types/pao_bas_projection_type [1.0] *
!!
!!   NAME
!!     pao_bas_projection_type
!!
!!   FUNCTION
!!     represent a simple projection (ie a subset of indexes)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - initialized: true if initialized (bug catcher)
!!     - proj_indexes: the local indexes of the full basis that are kept
!!       by the projection
!!     - non_proj_indexes: the local indexes of the full basis that are
!!       not in the min_bas (i.e. the indexes of the excluded basis)
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
  type pao_bas_projection_type
     private
     logical :: initialized
     integer, dimension(:), pointer :: proj_indexes
     integer, dimension(:), pointer :: non_proj_indexes
  end type pao_bas_projection_type
!!***
  !****************************************************************************

contains

!!****f* pao_bas_projection/cp_proj_create [1.0] *
!!
!!   NAME
!!     cp_proj_create
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_create(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Pointer:: proj
!!       Integer, Dimension(:), Intent (IN):: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_create
!!
!!   FUNCTION
!!     allocates and initianizes a pao_bas_projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     see cp_proj_init
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine cp_proj_create(proj, proj_indexes, size_full_bas, &
     error)
  type(pao_bas_projection_type), pointer :: proj
  integer, dimension(:), intent(in) ::proj_indexes
  integer, intent(in) :: size_full_bas
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: stat
  character(len=*), parameter :: routineN='cp_proj_create',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  allocate(proj,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     call cp_proj_init(proj, proj_indexes, size_full_bas, error)
  end if
end subroutine cp_proj_create
!***************************************************************************

!!****f* pao_bas_projection/cp_proj_create2 [1.0] *
!!
!!   NAME
!!     cp_proj_create2
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_create2(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Pointer:: proj
!!       Type(cp_sll_int_type), Pointer:: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_create2
!!
!!   FUNCTION
!!     allocates and initianizes a pao_bas_projection (from a linked list)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     see cp_proj_init2
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
subroutine cp_proj_create2(proj, proj_indexes, size_full_bas, &
     error)
  type(pao_bas_projection_type), pointer :: proj
  type(cp_sll_int_type), pointer ::proj_indexes
  integer, intent(in) :: size_full_bas
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: stat
  character(len=*), parameter :: routineN='cp_proj_create2',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  allocate(proj,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  if (.not.failure) then
     call cp_proj_init2(proj, proj_indexes, size_full_bas, error)
  end if
end subroutine cp_proj_create2
!***************************************************************************

!!****f* pao_bas_projection/cp_proj_dealloc [1.0] *
!!
!!   NAME
!!     cp_proj_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_dealloc(proj, error)
!!       Type(pao_bas_projection_type), Pointer:: proj
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_dealloc
!!
!!   FUNCTION
!!     deallocates a projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     proj: the projection to dealloc
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
subroutine cp_proj_dealloc(proj,error)
  type(pao_bas_projection_type), pointer :: proj
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: stat
  character(len=*), parameter :: routineN='cp_proj_dealloc',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  if (associated(proj)) then
     call cp_proj_dealloc_ref(proj,error=error)
     deallocate(proj,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  end if
end subroutine cp_proj_dealloc
!***************************************************************************

!!****f* pao_bas_projection/cp_proj_init [1.0] *
!!
!!   NAME
!!     cp_proj_init
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_init(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Intent (OUT):: proj
!!       Integer, Dimension(:), Intent (IN):: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_init
!!
!!   FUNCTION
!!     initializes the given projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to be initialized
!!     - proj_indexes: the indexes of the full basis that are kept
!!     - size_full_bas: the size of the full basis
!!     - non_proj_indexes: the index of the excluded basis
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
  subroutine cp_proj_init(proj, proj_indexes, size_full_bas, &
       error)
    type(pao_bas_projection_type), intent(out) :: proj
    integer, dimension(:), intent(in) ::proj_indexes
    integer, intent(in) :: size_full_bas
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure,found
    integer :: stat, i, j, exc_pos
    character(len=*), parameter :: routineN='cp_proj_init',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(all(proj_indexes>0),cp_failure_level,routineP,error,failure)
    CPPrecondition(size(proj_indexes)<=size_full_bas,cp_failure_level,routineP,error,failure)
    allocate(proj%proj_indexes(size(proj_indexes)),&
         proj%non_proj_indexes(size_full_bas-size(proj_indexes)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !do i=2,size(proj_indexes)
    !   CPPrecondition(proj_indexes(i-1)<proj_indexes(i),cp_failure_level,routineP,error,failure)
    !end do
    if (.not.failure) then
       proj%proj_indexes=proj_indexes
       exc_pos=0
       do i=1,size_full_bas
          found=.false.
          do j=1,size(proj%proj_indexes)
             if (i==proj%proj_indexes(j)) then
                found=.true.
                exit
             end if
          end do
          if (.not.found) then
             CPInvariant(exc_pos<size(proj%non_proj_indexes),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                exc_pos=exc_pos+1
                proj%non_proj_indexes(exc_pos)=i
             end if
          end if
       end do
       proj%initialized=.true.
    end if
  end subroutine cp_proj_init
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_init2 [1.0] *
!!
!!   NAME
!!     cp_proj_init2
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_init2(proj, proj_indexes, size_full_bas, error)
!!       Type(pao_bas_projection_type), Intent (OUT):: proj
!!       Type(cp_sll_int_type), Pointer:: proj_indexes
!!       Integer, Intent (IN):: size_full_bas
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_init2
!!
!!   FUNCTION
!!     initializes the given projection using a linked list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to be initialized
!!     - proj_indexes: the indexes of the full basis that are kept
!!     - size_full_bas: the size of the full basis
!!     - non_proj_indexes: the index of the excluded basis
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
  subroutine cp_proj_init2(proj, proj_indexes, size_full_bas, &
       error)
    type(pao_bas_projection_type), intent(out) :: proj
    type(cp_sll_int_type), pointer ::proj_indexes
    integer, intent(in) :: size_full_bas
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure,found
    integer :: stat, size_min_bas, i, j, exc_pos, el_att
    character(len=*), parameter :: routineN='cp_proj_init2',&
         routineP=moduleN//':'//routineN
    type(cp_sll_int_type), pointer :: iterator
    failure=.false.

    iterator => proj_indexes
    size_min_bas=0
    do
       if (.not.cp_next(iterator, el_att=el_att, error=error)) exit
       size_min_bas=size_min_bas+1
       CPPrecondition(el_att>0,cp_failure_level,routineP,error,failure)
    end do
    CPPrecondition(size_min_bas<=size_full_bas,cp_failure_level,routineP,error,failure)
    CPPreconditionNoFail(size_min_bas>0,cp_warning_level,routineP,error)
    allocate(proj%proj_indexes(size_min_bas),&
         proj%non_proj_indexes(size_full_bas-size_min_bas), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       iterator => proj_indexes
       do i=1,size_min_bas
          if (.not.cp_next(iterator, el_att=el_att, error=error)) exit
          proj%proj_indexes(i)=el_att
       end do
       exc_pos=0
       do i=1,size_full_bas
          found=.false.
          do j=1,size(proj%proj_indexes)
             if (i==proj%proj_indexes(j)) then
                found=.true.
                exit
             end if
          end do
          if (.not.found) then
             CPInvariant(exc_pos<size(proj%non_proj_indexes),cp_failure_level,routineP,error,failure)
             if (.not.failure) then
                exc_pos=exc_pos+1
                proj%non_proj_indexes(exc_pos)=i
             end if
          end if
       end do
       proj%initialized=.true.
    end if
  end subroutine cp_proj_init2
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_dealloc_ref [1.0] *
!!
!!   NAME
!!     cp_proj_dealloc_ref
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_dealloc_ref(proj, error)
!!       Type(pao_bas_projection_type), Intent (INOUT):: proj
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory allocated by this object
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to be deallocated
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
  subroutine cp_proj_dealloc_ref(proj, error)
    type(pao_bas_projection_type), intent(inout) :: proj
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, stat
    character(len=*), parameter :: routineN='cp_proj_dealloc_ref',&
         routineP=moduleN//':'//routineN
    failure=.false.

    !call timeset(routineN//','//moduleN,'I',"",handle)
    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       deallocate(proj%proj_indexes, proj%non_proj_indexes, stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       proj%initialized=.false.
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_proj_dealloc_ref
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_valid [1.0] *
!!
!!   NAME
!!     cp_proj_valid
!!
!!   SYNOPSIS
!!     Function cp_proj_valid(proj, error) Result(res)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Logical:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_proj_valid
!!
!!   FUNCTION
!!     returns true if the projection is valid (minimal tests)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to test
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
  function cp_proj_valid(proj,error) result(res)
    type(pao_bas_projection_type), intent(in) :: proj
    logical ::res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_proj_valid',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(proj%proj_indexes),cp_warning_level,routineP,error,failure)
    CPPrecondition(associated(proj%non_proj_indexes),cp_warning_level,routineP,error,failure)
    res=.not.failure
  end function cp_proj_valid
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_get [1.0] *
!!
!!   NAME
!!     cp_proj_get
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_get(proj, size_min_bas, size_full_bas,&
!!         size_excl_bas, proj_indexes, non_proj_indexes, error)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Integer, Optional, Intent (OUT):: size_excl_bas, size_full_bas,&
!!         size_min_bas
!!       Integer, Optional, Dimension(:), Pointer:: proj_indexes,&
!!         non_proj_indexes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_get
!!
!!   FUNCTION
!!     returns various attributes about the given projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection object
!!     - size_min_bas: the size of the minimal basis
!!     - size_full_bas: the size of the full basis
!!     - size_excl_bas: the size of the excluded basis
!!     - proj_indexes: the local indexes of the full basis that are kept
!!       by the projection
!!     - non_proj_indexes: the local indexes of the full basis that are
!!       not in the min_bas (i.e. the indexes of the excluded basis)
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
  subroutine cp_proj_get(proj, size_min_bas, size_full_bas, &
       size_excl_bas, proj_indexes, non_proj_indexes, error)
    type(pao_bas_projection_type), intent(in) :: proj
    integer, optional, intent(out) :: size_excl_bas,size_full_bas,&
         size_min_bas
    integer, optional, dimension(:), pointer :: proj_indexes,&
         non_proj_indexes
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='',&
         routineP=moduleN//':'//routineN
    failure=.false.

    CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       if (present(size_min_bas)) size_min_bas=size(proj%proj_indexes)
       if (present(size_full_bas)) &
            size_full_bas=size(proj%proj_indexes)+size(proj%non_proj_indexes)
       if (present(size_excl_bas)) size_excl_bas=size(proj%non_proj_indexes)
       if (present(proj_indexes)) proj_indexes => proj%proj_indexes
       if (present(non_proj_indexes)) non_proj_indexes => proj%non_proj_indexes
    end if
  end subroutine cp_proj_get
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_get_proj_indexes [1.0] *
!!
!!   NAME
!!     cp_proj_get_proj_indexes
!!
!!   SYNOPSIS
!!     Function cp_proj_get_proj_indexes(proj, error) Result(res)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Integer, Dimension(:), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_proj_get_proj_indexes
!!
!!   FUNCTION
!!     returns the projection indexes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to get the info from
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
  function cp_proj_get_proj_indexes(proj, error) result(res)
    type(pao_bas_projection_type), intent(in) :: proj
    integer, dimension(:), pointer :: res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_proj_get_proj_indexes',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug.and.debug_this_module) then
       CPPrecondition(cp_proj_valid(proj,error=error),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       res => proj%proj_indexes
    else
       nullify(res)
    end if
  end function cp_proj_get_proj_indexes
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_get_non_proj_indexes [1.0] *
!!
!!   NAME
!!     cp_proj_get_non_proj_indexes
!!
!!   SYNOPSIS
!!     Function cp_proj_get_non_proj_indexes(proj, error) Result(res)
!!       Type(pao_bas_projection_type), Intent (IN):: proj
!!       Integer, Dimension(:), Pointer:: res
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_proj_get_non_proj_indexes
!!
!!   FUNCTION
!!     returns the indexes excluded by th projection
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to get the info from
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
  function cp_proj_get_non_proj_indexes(proj, error) result(res)
    type(pao_bas_projection_type), intent(in) :: proj
    integer, dimension(:), pointer :: res
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='cp_proj_get_non_proj_indexes',&
         routineP=moduleN//':'//routineN
    failure=.false.

    if (cp_debug.and.debug_this_module) then
       CPPrecondition(cp_proj_valid(proj,error=error),cp_failure_level,routineP,error,failure)
    end if
    if (.not. failure) then
       res => proj%non_proj_indexes
    else
       nullify(res)
    end if
  end function cp_proj_get_non_proj_indexes
  !***************************************************************************

!!****f* pao_bas_projection/cp_proj_write [1.0] *
!!
!!   NAME
!!     cp_proj_write
!!
!!   SYNOPSIS
!!     Subroutine cp_proj_write(proj, unit_nr, error)
!!       Type(pao_bas_projection_type), Intent (IN), Target:: proj
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_proj_write
!!
!!   FUNCTION
!!     writes a description of the projection to the given unit
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - proj: the projection to output
!!     - unit_nr: the unit number to write to
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
subroutine cp_proj_write(proj,unit_nr,error)
  type(pao_bas_projection_type), intent(in), target :: proj
  integer, intent(in) :: unit_nr
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_proj_write',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(proj%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     write(unit=unit_nr,fmt="(t2,a)")'# projection indexes'
     call cp_write(proj%proj_indexes,unit_nr=unit_nr,error=error)
     write(unit=unit_nr,fmt="(t2,a)")'# non projection indexes'
     call cp_write(proj%non_proj_indexes,unit_nr=unit_nr,error=error)
  end if
end subroutine cp_proj_write
!***************************************************************************

end module pao_bas_projection_types
