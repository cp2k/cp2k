!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_linked_list_int_ptr [1.0] *
!!
!!   NAME
!!     cp_linked_list_int_ptr
!!
!!   FUNCTION
!!     describes a generic linked list template.
!!     Linked list are supposed to always use pointers to the nodes for
!!     basically everything: a pointer to the node is a list, an element of
!!     the list, an iterator between the elment of the list.
!!     An empty list is represented by an unassociated pointer.
!!
!!   NOTES
!!     template with one parameter. nametype1 should have less than 10 chars.
!!     DO NOT modify the .F files, modify the .template file
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     1.2002 created
!!     4.2002 changed most of the methods, by making access to the list 
!!            always through pointers (identifying pointer, list and iterators)
!!
!!   SOURCE
!*****************************************************************************
module cp_linked_list_int_ptr
  use cp_log_handling
  use cp_error_handling
  use cp_output_handling
  use kinds, only: wp=>dp
  use timings, only: timeset, timestop
  use cp_array_i_utils, only: cp_i_write
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_linked_list_int_ptr'

! type
  public :: cp_sll_int_ptr_type
! common methods
  public :: cp_create, cp_dealloc, cp_get, cp_set, cp_next, cp_output
! special get
  public :: cp_get_first_el, cp_get_rest, cp_get_empty, cp_get_length,&
       cp_get_element_at
! special set
  public :: cp_set_element_at
! structure manipulation
  public :: cp_insert, cp_remove_first_el, cp_remove_el, cp_remove_all
! low level
  public :: cp_dealloc_node, cp_write

  ! underlying routines
  public :: cp_sll_int_ptr_create, cp_sll_int_ptr_dealloc,&
       cp_sll_int_ptr_dealloc_node,cp_sll_int_ptr_set,&
       cp_sll_int_ptr_get,cp_sll_int_ptr_next,&
       cp_sll_int_ptr_get_first_el, cp_sll_int_ptr_get_rest,&
       cp_sll_int_ptr_get_empty, cp_sll_int_ptr_get_length,&
       cp_sll_int_ptr_get_el_at, cp_sll_int_ptr_set_el_at,&
       cp_sll_int_ptr_insert_el, cp_sll_int_ptr_insert_el_at,&
       cp_sll_int_ptr_rm_first_el, cp_sll_int_ptr_rm_el_at,&
       cp_sll_int_ptr_rm_all_el, cp_sll_int_ptr_output,&
       cp_sll_int_ptr_write

! creation of an object (from a pointer)
  interface cp_create
     module procedure cp_sll_int_ptr_create
  end interface
! destruction of an object (from a pointer)
  interface cp_dealloc
     module procedure cp_sll_int_ptr_dealloc
  end interface
! destruction only of the node (low level)
  interface cp_dealloc_node
     module procedure cp_sll_int_ptr_dealloc_node
  end interface
! modifies attributes of an object
  interface cp_set
     module procedure cp_sll_int_ptr_set
  end interface
! returns attributes of an object
  interface cp_get
     module procedure cp_sll_int_ptr_get
  end interface
! iterates to the next element
  interface cp_next
     module procedure cp_sll_int_ptr_next
  end interface
! returns the first element
  interface cp_get_first_el
     module procedure cp_sll_int_ptr_get_first_el
  end interface
! returns the rest of the list
  interface cp_get_rest
     module procedure cp_sll_int_ptr_get_rest
  end interface
! returns if the list is empty
  interface cp_get_empty
     module procedure cp_sll_int_ptr_get_empty
  end interface
! returns the length of the list
  interface cp_get_length
     module procedure cp_sll_int_ptr_get_length
  end interface
! returns the element at the given position
  interface cp_get_element_at
     module procedure cp_sll_int_ptr_get_el_at
  end interface
! sets the element at the given position
  interface cp_set_element_at
     module procedure cp_sll_int_ptr_set_el_at
  end interface
! inserts one element call cp_insert(list,element,...)
  interface cp_insert
     module procedure cp_sll_int_ptr_insert_el, &
          cp_sll_int_ptr_insert_el_at
  end interface
! removes an element
  interface cp_remove_el
     module procedure cp_sll_int_ptr_rm_first_el, &
          cp_sll_int_ptr_rm_el_at
  end interface
! removes the first el
  interface cp_remove_first_el
     module procedure cp_sll_int_ptr_rm_first_el
  end interface
! remove all the elments
  interface cp_remove_all
     module procedure cp_sll_int_ptr_rm_all_el
  end interface
! outputs the list
  interface cp_output
     module procedure cp_sll_int_ptr_output
  end interface
! writes the list (low-level)
  interface cp_write
     module procedure cp_sll_int_ptr_write
  end interface
!!***
!****************************************************************************

!!****s* cp2k/cp_sll_int_ptr_type [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_type
!!
!!   FUNCTION
!!     represent a single linked list that stores pointers to the elements
!!
!!   NOTES
!!     List are alway accessed through pointers, so every node of the 
!!     linked list can be seen as a list, its first element
!!     a pointer to the position before itself, in a very natural way:
!!     all the insertions take place before the actual element, and 
!!     you still can insert an element at the end.
!!     This way I could identify nodes, lists and pointers between the
!!     elements of the list.
!!     Indexing is 1 based.
!!
!!   ATTRIBUTES
!!     initialized: true if the node was initialized (bugcatcher).
!!                in this structure it can eat up quite some space, remove?
!!     first_el: the element that is stored in this node. 
!!     rest: the rest of the list
!!
!!     empty: true if the list pointer is not associated, if it points to
!!            to a not it is always false (as there is at least the 
!!            first_el in the list)
!!     length: the number of elements in the list
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!***************************************************************************
  type cp_sll_int_ptr_type
     logical :: initialized
     integer, pointer :: first_el
     type(cp_sll_int_ptr_type), pointer :: rest
  end type cp_sll_int_ptr_type
!!***
!****************************************************************************

contains

! =========== creation / distruction ========

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_create [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_create
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_create(sll, first_el, rest, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Intent (IN), Target, Optional:: first_el
!!       Type(cp_sll_int_ptr_type), Intent (IN), Target, Optional:: rest
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_create
!!
!!   FUNCTION
!!    allocates and initializes a single linked list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     sll: the single linked list to initialize
!!     first_el: the first element of this list
!!     rest: the following elements (if not given: empty)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_create(sll,first_el,rest,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, intent(in),target, optional :: first_el
    type(cp_sll_int_ptr_type), intent(in), target , optional :: rest
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_create'
    integer :: stat
    type(cp_sll_int_ptr_type), pointer :: empty_list
    failure=.false.

    if (.not.present(first_el)) then
       nullify(sll)
       if (present(rest)) sll => rest
    else
       allocate(sll, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          sll%first_el => first_el
          nullify(sll%rest)
          if (present(rest)) sll%rest => rest
          sll%initialized=.true.
       end if
    end if
    if (failure) nullify(sll)
  end subroutine cp_sll_int_ptr_create
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_dealloc(sll, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_dealloc
!!
!!   FUNCTION
!!     deallocates the singly linked list starting at sll.
!!     Does not work if loops are present!
!!
!!   NOTES
!!     does not deallocate the elments that are stored in the list
!!     check more?
!!
!!   INPUTS
!!     sll: the list to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_dealloc(sll,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc'

    call cp_sll_int_ptr_rm_all_el(sll,error)
  end subroutine cp_sll_int_ptr_dealloc
!***************************************************************************

! * low-level *

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc_node [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_dealloc_node
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_dealloc_node(sll, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_dealloc_node
!!
!!   FUNCTION
!!     deallocates a node of a singly linked list (low level)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     sll: the node to be deallocated
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_dealloc_node(sll,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc_node'
    integer:: stat
    failure=.false.

    CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
    deallocate(sll, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  end subroutine cp_sll_int_ptr_dealloc_node
!***************************************************************************

! ============= get/set ============

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_set [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_set
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_set(sll, first_el, rest, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Intent (IN), Target, Optional:: first_el
!!       Type(cp_sll_int_ptr_type), Intent (IN), Optional, Target:: rest
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_set
!!
!!   FUNCTION
!!     sets tha various attributes of a single linked list
!!     It is an error to set the rest of an empty list without setting 
!!     the first_el (change and make it simply a sll=>rest ?)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     sll: the single linked list to change
!!     first_el: the element to replace the first element of this list
!!     rest: the rest of the list (can be unassociated)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_set(sll,first_el,rest,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, intent(in),target, optional :: first_el
    type(cp_sll_int_ptr_type), intent(in), optional, target :: rest
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_set_sll_int_ptr'
    failure=.false.

    if (.not.associated(sll)) then
       if (present(first_el)) then
          call cp_sll_int_ptr_create(sll,first_el,rest,error)
       else
          CPAssert(.not.present(rest),cp_failure_level,routineP,error,failure)
       end if
    else
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          if (present(first_el)) sll%first_el => first_el
          if (present(rest)) sll%rest => rest
       end if
    end if
  end subroutine cp_sll_int_ptr_set
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_get(sll, first_el, rest, empty, length,&
!!         error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Pointer, Optional:: first_el
!!       Type(cp_sll_int_ptr_type), Optional, Pointer:: rest
!!       Logical, Optional, Intent (OUT):: empty
!!       Integer, Intent (OUT), Optional:: length
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_get
!!
!!   FUNCTION
!!     returns various attributes from the linked list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     sll: the single linked list to change
!!     first_el: the first element of this list (can be unassociated)
!!     rest: the rest of the list (can be unassociated)
!!     length: the length of the list
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_get(sll,first_el,rest,empty,length,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, pointer, optional :: first_el
    type(cp_sll_int_ptr_type), optional, pointer :: rest
    logical, optional, intent(out) :: empty
    integer, intent(out), optional :: length
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_get'
    failure=.false.

    if (.not.associated(sll)) then
       CPPrecondition(.not.present(first_el),cp_failure_level,routineP,error,failure)
       if (present(rest)) nullify(rest)
       if (present(empty)) empty=.true.
       if (present(length)) length=0
    else
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          if (present(first_el)) first_el => sll%first_el
          if (present(rest)) rest => sll%rest
          if (present(empty)) empty = .false.
          if (present(length)) &
               length = cp_sll_int_ptr_get_length(sll,error=error)
       end if
    end if
  end subroutine cp_sll_int_ptr_get
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_first_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_first_el
!!
!!   SYNOPSIS
!!     Function cp_sll_int_ptr_get_first_el(sll, error) Result(res)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Integer, Pointer:: res
!!     End Function cp_sll_int_ptr_get_first_el
!!
!!   FUNCTION
!!     returns the first element stored in the list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     sll: the single linked list to get the element from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_sll_int_ptr_get_first_el(sll,error) result(res)
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error
    integer, pointer :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_get_first'
    failure=.false.

    if (cp_debug) then
       CPPrecondition(associated(sll),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       end if
    end if
    if (.not. failure) then
       res => sll%first_el
    end if
  end function cp_sll_int_ptr_get_first_el
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_rest [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_rest
!!
!!   SYNOPSIS
!!     Function cp_sll_int_ptr_get_rest(sll, iter, error) Result(res)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Optional:: iter
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Type(cp_sll_int_ptr_type), Pointer:: res
!!     End Function cp_sll_int_ptr_get_rest
!!
!!   FUNCTION
!!     returns the rest of the list
!!
!!   NOTES
!!     split the case iter=1 to make it more optimized?
!!
!!   INPUTS
!!     sll: the single linked list to get the rest from
!!     iter: how many times the call to rest should be iterated,
!!           defaults to 1; -1 means till end of the list.
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_sll_int_ptr_get_rest(sll, iter, error) result(res)
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, optional :: iter
    type(cp_error_type), optional, intent(inout) :: error
    type(cp_sll_int_ptr_type), pointer :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_get_first'
    integer :: i
    failure=.false.

    if (.not.associated(sll)) then
       nullify(res)
    else
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          if (present(iter)) then
             res => sll
             do i=1,iter
                if (associated(res%rest)) then
                   res => res%rest
                else
                   call cp_assert(.false.,cp_warning_level,cp_assertion_failed,&
                        routineP, "tried to go past end in "// &
                        CPSourceFileRef,&
                        error,failure)
                end if
             end do
             if (iter==-1) then
                res => sll
                do i=1,iter
                   if (associated(res%rest)) then
                      res => res%rest
                   else
                      call cp_assert(.false.,cp_warning_level,cp_assertion_failed,&
                           routineP, "tried to go past end in "// &
                           CPSourceFileRef,&
                           error,failure)
                   end if
                end do
             end if
          else
             res => sll%rest ! make the common case fast...
          end if
       else
          nullify(res)
       end if
    end if
  end function cp_sll_int_ptr_get_rest
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_empty [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_empty
!!
!!   SYNOPSIS
!!     Function cp_sll_int_ptr_get_empty(sll, error) Result(res)
!!       Logical:: res
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sll_int_ptr_get_empty
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     sll: the single linked list to get the rest from
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_sll_int_ptr_get_empty(sll,error) result(res)
    logical :: res
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_get_empty'
    failure=.false.

    res = .not.associated(sll)
  end function cp_sll_int_ptr_get_empty
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_length [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_length
!!
!!   SYNOPSIS
!!     Function cp_sll_int_ptr_get_length(sll, error) Result(res)
!!       Type(cp_sll_int_ptr_type), Intent (IN), Target:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Integer:: res
!!     End Function cp_sll_int_ptr_get_length
!!
!!   FUNCTION
!!     returns the length of the list
!!
!!   NOTES
!!     slow (O(n))
!!
!!   INPUTS
!!     sll: the list you want to know the length of
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_sll_int_ptr_get_length(sll,error) result(res)
    type(cp_sll_int_ptr_type), intent(in), target :: sll
    type(cp_error_type), optional, intent(inout) :: error
    integer ::res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_get_length'
    type(cp_sll_int_ptr_type), pointer :: iterator
    failure=.false.

    res=0
    iterator => sll
    do
       if (associated(iterator)) then
          res=res+1
          iterator => iterator%rest
       else
          exit
       end if
    end do
  end function cp_sll_int_ptr_get_length
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_el_at
!!
!!   SYNOPSIS
!!     Function cp_sll_int_ptr_get_el_at(sll, index, error) Result(res)
!!       Integer, Pointer:: res
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Intent (IN):: index
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Function cp_sll_int_ptr_get_el_at
!!
!!   FUNCTION
!!     returns the element at the given index
!!
!!   NOTES
!!     slow (O(index))
!!
!!   INPUTS
!!     sll: the list you get the element from
!!     index: the position of the element (stating at 1)
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_sll_int_ptr_get_el_at(sll,index,error) result(res)
    integer, pointer :: res
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, intent(in) :: index
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_get_el_at'
    type(cp_sll_int_ptr_type), pointer :: pos
    failure=.false.

    if (cp_debug) then
       CPPrecondition(index>0.or.index==-1,cp_failure_level,routineP,error,failure)
    end if
    if (index==-1) then
       pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
    else
       pos => cp_sll_int_ptr_get_rest(sll, iter=index-1,error=error)
    end if
    CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       res => pos%first_el
    end if
  end function cp_sll_int_ptr_get_el_at
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_set_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_set_el_at
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_set_el_at(sll, index, value, error)
!!       Integer, Intent (IN), Target:: value
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Intent (IN):: index
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_set_el_at
!!
!!   FUNCTION
!!     sets the element at the given index
!!
!!   NOTES
!!     slow (O(index))
!!
!!   INPUTS
!!     sll: the list you get the element from
!!     index: the position of the element (stating at 1)
!!            -1 means at the end
!!     value: the new element
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_set_el_at(sll,index,value,error)
    integer, intent(in),target :: value
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, intent(in) :: index
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_set_el_at'
    type(cp_sll_int_ptr_type), pointer :: pos
    failure=.false.

    if (index==-1) then
       pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
    else
       pos => cp_sll_int_ptr_get_rest(sll, iter=index-1,error=error)
    end if
    CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       pos%first_el => value
    end if
  end subroutine cp_sll_int_ptr_set_el_at
!***************************************************************************

! * iteration *

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_next [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_next
!!
!!   SYNOPSIS
!!     Function cp_sll_int_ptr_next(iterator, el_att, error) Result(res)
!!       Type(cp_sll_int_ptr_type), Pointer:: iterator
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Integer, Pointer, Optional:: el_att
!!       Logical:: res
!!     End Function cp_sll_int_ptr_next
!!
!!   FUNCTION
!!     returns true if the actual element is valid (i.e. iterator ont at end)
!!     moves the iterator to the next element
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     iterator: iterator that moves along the list
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!     el_att: the actual element (valid only if the function returns true)
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  function cp_sll_int_ptr_next(iterator,el_att,error) result(res)
    type(cp_sll_int_ptr_type), pointer :: iterator
    type(cp_error_type), optional, intent(inout) :: error
    integer, pointer, optional :: el_att
    logical :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_next'
    failure=.false.

    if (associated(iterator)) then
       if (cp_debug) then
          CPPrecondition(iterator%initialized,cp_failure_level,routineP,error,failure)
       end if
       res=.not.failure
       el_att => iterator%first_el
       iterator => iterator%rest
    else
       res=.false.
    end if
  end function cp_sll_int_ptr_next
!***************************************************************************

! ============ structure modifications ============

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_insert_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_insert_el
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_insert_el(sll, element, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Intent (IN), Target:: element
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_insert_el
!!
!!   FUNCTION
!!     insert an element at the beginning of the list
!!
!!   NOTES
!!     fast (O(1))
!!
!!   INPUTS
!!     sll: the single linked list point at the beginning of which 
!!          you want to add the element
!!     element: the element to add
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_insert_el(sll,element,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, intent(in),target:: element
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr:cp_sll_int_ptr_add_el'
    type(cp_sll_int_ptr_type), pointer :: newSlot
    integer :: stat
    failure=.false.
    nullify(newSlot)

    call cp_sll_int_ptr_create(newSlot,first_el=element,&
         rest=sll,error=error)
    sll => newSlot
  end subroutine cp_sll_int_ptr_insert_el
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_rm_first_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_rm_first_el
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_rm_first_el(sll, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_rm_first_el
!!
!!   FUNCTION
!!     remove the first element of the linked list
!!
!!   NOTES
!!     fast (O(1))
!!
!!   INPUTS
!!     sll: the list whose first element has to be removed
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_rm_first_el(sll,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_rm_first_el'
    type(cp_sll_int_ptr_type), pointer :: node_to_rm
    failure=.false.
    node_to_rm => sll

    if (associated(sll)) then
       sll => sll%rest
       call cp_sll_int_ptr_dealloc_node(node_to_rm,error=error)
    else
       call cp_assert(.false.,cp_warning_level,cp_assertion_failed,&
            routineP,"tried to remove first el of an empty list in "//&
            CPSourceFileRef,&
            error,failure)
    end if
  end subroutine cp_sll_int_ptr_rm_first_el
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_insert_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_insert_el_at
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_insert_el_at(sll, index, element, error)
!!       Integer, Intent (IN), Target:: element
!!       Integer, Intent (IN):: index
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_insert_el_at
!!
!!   FUNCTION
!!     inserts the element at the given index
!!
!!   NOTES
!!     slow (O(index))
!!
!!   INPUTS
!!     sll: the list you get the element from
!!     index: the position of the element (stating at 1). 
!!            If it is -1, it means at end
!!     element: the new element
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_insert_el_at(sll,index,element,error)
    integer, intent(in),target :: element
    integer, intent(in) :: index
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_insert_el_at'
    type(cp_sll_int_ptr_type), pointer :: pos
    failure=.false.

    if (index==0) then
       call cp_sll_int_ptr_insert_el(sll,element,error=error)
    else
       if (index==-1) then
          pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
       else
          pos => cp_sll_int_ptr_get_rest(sll, iter=index-1,error=error)
       end if
       CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          call cp_sll_int_ptr_insert_el(pos%rest,element,error=error)
       end if
    end if
  end subroutine cp_sll_int_ptr_insert_el_at
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_rm_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_rm_el_at
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_rm_el_at(sll, index, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Integer, Intent (IN):: index
!!     End Subroutine cp_sll_int_ptr_rm_el_at
!!
!!   FUNCTION
!!     removes the element at the given index
!!
!!   NOTES
!!     slow (O(index))
!!
!!   INPUTS
!!     sll: the list you get the element from
!!     index: the position of the element (stating at 1)
!!     element: the new element
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_rm_el_at(sll,index,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error
    integer, intent(in)::index

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_remove_el_at'
    type(cp_sll_int_ptr_type), pointer :: pos
    failure=.false.
    
    if (cp_debug) then
       CPPrecondition(index>0.or.index==-1,cp_failure_level,routineP,error,failure)
    end if
    if (index==1) then
       call cp_sll_int_ptr_rm_first_el(sll,error=error)
    else
       if (index==-1) then
          pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
       else
          pos => cp_sll_int_ptr_get_rest(sll, iter=index-2,error=error)
       end if
       CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          call cp_sll_int_ptr_rm_first_el(pos%rest,error=error)
       end if
    end if
  end subroutine cp_sll_int_ptr_rm_el_at
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_rm_all_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_rm_all_el
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_rm_all_el(sll, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_rm_all_el
!!
!!   FUNCTION
!!     removes all the elements from the list
!!
!!   NOTES
!!     check more?
!!
!!   INPUTS
!!     sll: the list that should be removed
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     @author Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  subroutine cp_sll_int_ptr_rm_all_el(sll,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int_ptr/cp_sll_int_ptr_rm_all_el'
    type(cp_sll_int_ptr_type), pointer :: next_node, actual_node
    integer:: stat
    failure=.false.

    actual_node => sll
    do
       if (.not.associated(actual_node)) exit
       next_node => cp_get_rest(actual_node,error=error)
       call cp_sll_int_ptr_dealloc_node(actual_node,error=error)
       actual_node => next_node
    end do
    nullify(sll)
  end subroutine cp_sll_int_ptr_rm_all_el
!***************************************************************************

! ========= output =========

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_output [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_output
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_output(logger, outputname, fromwhere, iter,&
!!         sll, comment, error)
!!       Type(cp_logger_type):: logger
!!       Character(Len=*), Intent (IN):: outputname
!!       Character(Len=*), Intent (IN):: fromwhere
!!       Integer, Intent (IN):: iter
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Character(Len=*), Intent (IN), Optional:: comment
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_output
!!
!!   FUNCTION
!!     outputs the contents of the list as a comma separated parentized list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     logger: the logger that decides where the logging should go
!!     outputName: the name of the output, for example 'hamiltonian'
!!     fromWhere: string of the form module:function or file:lineNr
!!                that says where the error happend
!!     iter: the iteration number, or some other number that permit to 
!!           differentiate between different instances of the same object
!!     sll: the list to output
!!     comment: an (optional) comment that is written before the list
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
  subroutine cp_sll_int_ptr_output(logger, outputName,&
       fromWhere, iter, sll,comment, error)
    type(cp_logger_type) ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    type(cp_sll_int_ptr_type), pointer :: sll
    character(len=*), intent(in),optional :: comment
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    integer :: handle, unitNr, iostat
    character(len=*), parameter :: routineN='cp_sll_int_ptr_output',&
         routineP=moduleN//':'//routineN
    failure=.false.
  
    call timeset(routineN//','//moduleN,'I',"",handle)
    if (cp_would_output(logger, outputName , fromWhere, iter)) then
       unitNr=cp_unitnr_for_output(logger, outputName, fromWhere, iter)
       if (present(comment)) then
          write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
          CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
       end if
       call cp_sll_int_ptr_write(sll,unit_nr=unitNr,error=error)
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_sll_int_ptr_output
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_write [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_write
!!
!!   SYNOPSIS
!!     Subroutine cp_sll_int_ptr_write(sll, unit_nr, error)
!!       Type(cp_sll_int_ptr_type), Pointer:: sll
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_sll_int_ptr_write
!!
!!   FUNCTION
!!     writes the list as a comma separated parentized list to the given
!!     unit
!!
!!   NOTES
!!     low level
!!
!!   INPUTS
!!     sll: the list to output
!!     unit_nr: the unit to write to (defaults to the standard out)
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
  subroutine cp_sll_int_ptr_write(sll,unit_nr,error)
    type(cp_sll_int_ptr_type), pointer :: sll
    integer, intent(in) :: unit_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    type(cp_sll_int_ptr_type), pointer :: iterator
    character(len=*), parameter :: routineN='cp_sll_int_ptr_write',&
         routineP=moduleN//':'//routineN
    integer, pointer :: el
    failure=.false.
    iterator => sll
  
    do
       if (.not.cp_sll_int_ptr_next(iterator,el_att=el,error=error)) exit
       call cp_i_write(el, unit_nr=unit_nr, error=error)
    end do
  end subroutine cp_sll_int_ptr_write
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same (multi-line use change it a bit)
! [template(type1,type1in,type1out,nametype1,use,=,write_el)]
! ARGS:
!  type1in = "integer, intent(in),target"
!  type1 = "integer, pointer"
!  write_el = "cp_i_write"
!  type1out = "integer, pointer"
!  nametype1 = "int_ptr"
!  use = "use cp_array_i_utils, only: cp_i_write"
!  = = "=>"


end module cp_linked_list_int_ptr

