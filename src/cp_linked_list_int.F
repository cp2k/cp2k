!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"
#define CP_SLL_I_LESS_Q(el1,el2) ( el1 < el2 )
#define CP_SLL_I_EQUAL_Q(el1,el2) ( el1 == el2 )

!!****h* cp2k/cp_linked_list_int [1.0] *
!!
!!   NAME
!!     cp_linked_list_int
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
module cp_linked_list_int
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log,&
       cp_logger_get_default_unit_nr, cp_logger_type
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use cp_output_handling
  use kinds, only: wp=>dp
  use timings, only: timeset, timestop
  use cp_array_i_utils, only: cp_i_write
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='cp_linked_list_int'

! type
  public :: cp_sll_int_type
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
  public :: cp_sll_int_create, cp_sll_int_dealloc,&
       cp_sll_int_dealloc_node,cp_sll_int_set,&
       cp_sll_int_get,cp_sll_int_next,&
       cp_sll_int_get_first_el, cp_sll_int_get_rest,&
       cp_sll_int_get_empty, cp_sll_int_get_length,&
       cp_sll_int_get_el_at, cp_sll_int_set_el_at,&
       cp_sll_int_insert_el, cp_sll_int_insert_el_at,&
       cp_sll_int_rm_first_el, cp_sll_int_rm_el_at,&
       cp_sll_int_rm_all_el, cp_sll_int_output,&
       cp_sll_int_write, cp_sll_int_to_array,&
       cp_sll_int_from_array, cp_sll_int_insert_ordered,&
       cp_sll_int_insert_ordered2

! creation of an object (from a pointer)
  interface cp_create
     module procedure cp_sll_int_create
  end interface
! destruction of an object (from a pointer)
  interface cp_dealloc
     module procedure cp_sll_int_dealloc
  end interface
! destruction only of the node (low level)
  interface cp_dealloc_node
     module procedure cp_sll_int_dealloc_node
  end interface
! modifies attributes of an object
  interface cp_set
     module procedure cp_sll_int_set
  end interface
! returns attributes of an object
  interface cp_get
     module procedure cp_sll_int_get
  end interface
! iterates to the next element
  interface cp_next
     module procedure cp_sll_int_next
  end interface
! returns the first element
  interface cp_get_first_el
     module procedure cp_sll_int_get_first_el
  end interface
! returns the rest of the list
  interface cp_get_rest
     module procedure cp_sll_int_get_rest
  end interface
! returns if the list is empty
  interface cp_get_empty
     module procedure cp_sll_int_get_empty
  end interface
! returns the length of the list
  interface cp_get_length
     module procedure cp_sll_int_get_length
  end interface
! returns the element at the given position
  interface cp_get_element_at
     module procedure cp_sll_int_get_el_at
  end interface
! sets the element at the given position
  interface cp_set_element_at
     module procedure cp_sll_int_set_el_at
  end interface
! inserts one element call cp_insert(list,element,...)
  interface cp_insert
     module procedure cp_sll_int_insert_el
  end interface
  interface cp_insert_ordered
     module procedure cp_sll_int_insert_ordered,&
          cp_sll_int_insert_ordered2
  end interface
  interface cp_insert_at
     module procedure cp_sll_int_insert_el_at
  end interface
! removes an element
  interface cp_remove_el
     module procedure cp_sll_int_rm_first_el, &
          cp_sll_int_rm_el_at
  end interface
! removes the first el
  interface cp_remove_first_el
     module procedure cp_sll_int_rm_first_el
  end interface
! remove all the elments
  interface cp_remove_all
     module procedure cp_sll_int_rm_all_el
  end interface
! outputs the list
  interface cp_output
     module procedure cp_sll_int_output
  end interface
! writes the list (low-level)
  interface cp_write
     module procedure cp_sll_int_write
  end interface
!!***
!****************************************************************************

!!****s* cp2k/cp_sll_int_type [1.0] *
!!
!!   NAME
!!     cp_sll_int_type
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
  type cp_sll_int_type
     logical :: initialized
     integer :: first_el
     type(cp_sll_int_type), pointer :: rest
  end type cp_sll_int_type
!!***
!****************************************************************************

contains

! =========== creation / distruction ========

!!****f* cp_linked_list_int/cp_sll_int_create [1.0] *
!!
!!   NAME
!!     cp_sll_int_create
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
  subroutine cp_sll_int_create(sll,first_el,rest,error)
    type(cp_sll_int_type), pointer :: sll
    integer, intent(in), optional :: first_el
    type(cp_sll_int_type), pointer, optional :: rest
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_create'
    integer :: stat
    type(cp_sll_int_type), pointer :: empty_list
    failure=.false.

    if (.not.present(first_el)) then
       nullify(sll)
       if (present(rest)) sll => rest
    else
       allocate(sll, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          sll%first_el = first_el
          nullify(sll%rest)
          if (present(rest)) sll%rest => rest
          sll%initialized=.true.
       end if
    end if
    if (failure) nullify(sll)
  end subroutine cp_sll_int_create
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_dealloc [1.0] *
!!
!!   NAME
!!     cp_sll_int_dealloc
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
  subroutine cp_sll_int_dealloc(sll,error)
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_dealloc'

    call cp_sll_int_rm_all_el(sll,error)
  end subroutine cp_sll_int_dealloc
!***************************************************************************

! * low-level *

!!****f* cp_linked_list_int/cp_sll_int_dealloc_node [1.0] *
!!
!!   NAME
!!     cp_sll_int_dealloc_node
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
  subroutine cp_sll_int_dealloc_node(sll,error)
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_dealloc_node'
    integer:: stat
    failure=.false.

    CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
    deallocate(sll, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  end subroutine cp_sll_int_dealloc_node
!***************************************************************************

! ============= get/set ============

!!****f* cp_linked_list_int/cp_sll_int_set [1.0] *
!!
!!   NAME
!!     cp_sll_int_set
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
  subroutine cp_sll_int_set(sll,first_el,rest,error)
    type(cp_sll_int_type), pointer :: sll
    integer, intent(in), optional :: first_el
    type(cp_sll_int_type), pointer, optional :: rest
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_set_sll_int'
    failure=.false.

    if (.not.associated(sll)) then
       if (present(first_el)) then
          call cp_sll_int_create(sll,first_el,rest,error)
       else
          CPAssert(.not.present(rest),cp_failure_level,routineP,error,failure)
       end if
    else
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          if (present(first_el)) sll%first_el = first_el
          if (present(rest)) sll%rest => rest
       end if
    end if
  end subroutine cp_sll_int_set
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_get [1.0] *
!!
!!   NAME
!!     cp_sll_int_get
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
  subroutine cp_sll_int_get(sll,first_el,rest,empty,length,error)
    type(cp_sll_int_type), pointer :: sll
    integer, intent(out), optional :: first_el
    type(cp_sll_int_type), optional, pointer :: rest
    logical, optional, intent(out) :: empty
    integer, intent(out), optional :: length
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_get'
    failure=.false.

    if (.not.associated(sll)) then
       CPPrecondition(.not.present(first_el),cp_failure_level,routineP,error,failure)
       if (present(rest)) nullify(rest)
       if (present(empty)) empty=.true.
       if (present(length)) length=0
    else
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       if (.not. failure) then
          if (present(first_el)) first_el = sll%first_el
          if (present(rest)) rest => sll%rest
          if (present(empty)) empty = .false.
          if (present(length)) &
               length = cp_sll_int_get_length(sll,error=error)
       end if
    end if
  end subroutine cp_sll_int_get
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_get_first_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_get_first_el
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
  function cp_sll_int_get_first_el(sll,error) result(res)
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error
    integer :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_get_first'
    failure=.false.

    if (cp_debug) then
       CPPrecondition(associated(sll),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       end if
    end if
    if (.not. failure) then
       res = sll%first_el
    end if
  end function cp_sll_int_get_first_el
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_get_rest [1.0] *
!!
!!   NAME
!!     cp_sll_int_get_rest
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
  function cp_sll_int_get_rest(sll, iter, error) result(res)
    type(cp_sll_int_type), pointer :: sll
    integer, optional :: iter
    type(cp_error_type), optional, intent(inout) :: error
    type(cp_sll_int_type), pointer :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_get_first'
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
  end function cp_sll_int_get_rest
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_get_empty [1.0] *
!!
!!   NAME
!!     cp_sll_int_get_empty
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
  function cp_sll_int_get_empty(sll,error) result(res)
    logical :: res
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_get_empty'
    failure=.false.

    res = .not.associated(sll)
  end function cp_sll_int_get_empty
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_get_length [1.0] *
!!
!!   NAME
!!     cp_sll_int_get_length
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
  function cp_sll_int_get_length(sll,error) result(res)
    type(cp_sll_int_type), intent(in), target :: sll
    type(cp_error_type), optional, intent(inout) :: error
    integer ::res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_get_length'
    type(cp_sll_int_type), pointer :: iterator
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
  end function cp_sll_int_get_length
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_get_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_get_el_at
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
  function cp_sll_int_get_el_at(sll,index,error) result(res)
    integer :: res
    type(cp_sll_int_type), pointer :: sll
    integer, intent(in) :: index
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_get_el_at'
    type(cp_sll_int_type), pointer :: pos
    failure=.false.

    if (cp_debug) then
       CPPrecondition(index>0.or.index==-1,cp_failure_level,routineP,error,failure)
    end if
    if (index==-1) then
       pos => cp_sll_int_get_rest(sll, iter=-1,error=error)
    else
       pos => cp_sll_int_get_rest(sll, iter=index-1,error=error)
    end if
    CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       res = pos%first_el
    end if
  end function cp_sll_int_get_el_at
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_set_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_set_el_at
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
  subroutine cp_sll_int_set_el_at(sll,index,value,error)
    integer, intent(in) :: value
    type(cp_sll_int_type), pointer :: sll
    integer, intent(in) :: index
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_set_el_at'
    type(cp_sll_int_type), pointer :: pos
    failure=.false.

    if (index==-1) then
       pos => cp_sll_int_get_rest(sll, iter=-1,error=error)
    else
       pos => cp_sll_int_get_rest(sll, iter=index-1,error=error)
    end if
    CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
    if (.not.failure) then
       pos%first_el = value
    end if
  end subroutine cp_sll_int_set_el_at
!***************************************************************************

! * iteration *

!!****f* cp_linked_list_int/cp_sll_int_next [1.0] *
!!
!!   NAME
!!     cp_sll_int_next
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
  function cp_sll_int_next(iterator,el_att,error) result(res)
    type(cp_sll_int_type), pointer :: iterator
    type(cp_error_type), optional, intent(inout) :: error
    integer, intent(out), optional :: el_att
    logical :: res

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_next'
    failure=.false.

    if (associated(iterator)) then
       if (cp_debug) then
          CPPrecondition(iterator%initialized,cp_failure_level,routineP,error,failure)
       end if
       res=.not.failure
       el_att = iterator%first_el
       iterator => iterator%rest
    else
       res=.false.
    end if
  end function cp_sll_int_next
!***************************************************************************

! ============ structure modifications ============

!!****f* cp_linked_list_int/cp_sll_int_insert_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_insert_el
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
  subroutine cp_sll_int_insert_el(sll,element,error)
    type(cp_sll_int_type), pointer :: sll
    integer, intent(in):: element
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int:cp_sll_int_add_el'
    type(cp_sll_int_type), pointer :: newSlot
    integer :: stat
    failure=.false.
    nullify(newSlot)

    call cp_sll_int_create(newSlot,first_el=element,&
         rest=sll,error=error)
    sll => newSlot
  end subroutine cp_sll_int_insert_el
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_rm_first_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_rm_first_el
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
  subroutine cp_sll_int_rm_first_el(sll,error)
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_rm_first_el'
    type(cp_sll_int_type), pointer :: node_to_rm
    failure=.false.
    node_to_rm => sll

    if (associated(sll)) then
       sll => sll%rest
       call cp_sll_int_dealloc_node(node_to_rm,error=error)
    else
       call cp_assert(.false.,cp_warning_level,cp_assertion_failed,&
            routineP,"tried to remove first el of an empty list in "//&
            CPSourceFileRef,&
            error,failure)
    end if
  end subroutine cp_sll_int_rm_first_el
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_insert_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_insert_el_at
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
  subroutine cp_sll_int_insert_el_at(sll,element,index,error)
    integer, intent(in) :: element
    integer, intent(in) :: index
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_insert_el_at'
    type(cp_sll_int_type), pointer :: pos
    failure=.false.

    if (index==1) then
       call cp_sll_int_insert_el(sll,element,error=error)
    else
       if (index==-1) then
          pos => cp_sll_int_get_rest(sll, iter=-1,error=error)
       else
          pos => cp_sll_int_get_rest(sll, iter=index-1,error=error)
       end if
       CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          call cp_sll_int_insert_el(pos%rest,element,error=error)
       end if
    end if
  end subroutine cp_sll_int_insert_el_at
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_rm_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_rm_el_at
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
  subroutine cp_sll_int_rm_el_at(sll,index,error)
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error
    integer, intent(in)::index

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_remove_el_at'
    type(cp_sll_int_type), pointer :: pos
    failure=.false.
    
    if (cp_debug) then
       CPPrecondition(index>0.or.index==-1,cp_failure_level,routineP,error,failure)
    end if
    if (index==1) then
       call cp_sll_int_rm_first_el(sll,error=error)
    else
       if (index==-1) then
          pos => cp_sll_int_get_rest(sll, iter=-1,error=error)
       else
          pos => cp_sll_int_get_rest(sll, iter=index-2,error=error)
       end if
       CPPrecondition(associated(pos),cp_failure_level,routineP,error,failure)
       if (.not.failure) then
          call cp_sll_int_rm_first_el(pos%rest,error=error)
       end if
    end if
  end subroutine cp_sll_int_rm_el_at
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_rm_all_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_rm_all_el
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
  subroutine cp_sll_int_rm_all_el(sll,error)
    type(cp_sll_int_type), pointer :: sll
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineP='cp_linked_list_int/cp_sll_int_rm_all_el'
    type(cp_sll_int_type), pointer :: next_node, actual_node
    integer:: stat
    failure=.false.

    actual_node => sll
    do
       if (.not.associated(actual_node)) exit
       next_node => cp_get_rest(actual_node,error=error)
       call cp_sll_int_dealloc_node(actual_node,error=error)
       actual_node => next_node
    end do
    nullify(sll)
  end subroutine cp_sll_int_rm_all_el
!***************************************************************************

! ========= output =========

!!****f* cp_linked_list_int/cp_sll_int_output [1.0] *
!!
!!   NAME
!!     cp_sll_int_output
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
!!     local: it the output is local to this task. Defaults to false
!!     finish: if true (the default) calls cp_finished_output
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2002 created [fawzi]
!!     08.2002 added local/global output [fawzi]
!!
!!*** **********************************************************************
  subroutine cp_sll_int_output(logger, outputName,&
       fromWhere, iter, sll,comment,local, finish, error)
    type(cp_logger_type),pointer ::logger
    character(len=*) , intent(in) :: outputName
    character(len=*), intent(in):: fromWhere
    integer, intent(in) :: iter
    type(cp_sll_int_type), pointer :: sll
    character(len=*), intent(in),optional :: comment
    logical, intent(in), optional :: local, finish
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure, loc, fin
    integer :: handle, unitNr, iostat
    character(len=*), parameter :: routineN='cp_sll_int_output',&
         routineP=moduleN//':'//routineN
    failure=.false.; loc=.false.; fin=.true.
  
    call timeset(routineN//','//moduleN,'I',"",handle)
    if (present(local)) loc=local
    if (present(finish)) fin=finish

    if (loc .or. logger%para_env%mepos==logger%para_env%source) then
       if (cp_would_output(logger, outputName , fromWhere, iter,&
            local=loc,error=error)) then

          unitNr=cp_unitnr_for_output(logger, outputName, fromWhere, iter,&
               local=loc,error=error)
          
          if (present(comment)) then
             write(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
             CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
          end if
          call cp_sll_int_write(sll,unit_nr=unitNr,error=error)

          if (fin) then
             call cp_finished_output(logger, outputName, fromWhere, iter,&
                  local=local,error=error)
          end if
       end if
    end if
    call timestop(0.0_wp,handle)
  end subroutine cp_sll_int_output
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_write [1.0] *
!!
!!   NAME
!!     cp_sll_int_write
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
  subroutine cp_sll_int_write(sll,unit_nr,error)
    type(cp_sll_int_type), pointer :: sll
    integer, intent(in) :: unit_nr
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    type(cp_sll_int_type), pointer :: iterator
    character(len=*), parameter :: routineN='cp_sll_int_write',&
         routineP=moduleN//':'//routineN
    integer :: el
    failure=.false.
    iterator => sll
  
    do
       if (.not.cp_sll_int_next(iterator,el_att=el,error=error)) exit
       call cp_i_write(el, unit_nr=unit_nr, error=error)
    end do
  end subroutine cp_sll_int_write
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_to_array [1.0] *
!!
!!   NAME
!!     cp_sll_int_to_array
!!
!!   FUNCTION
!!     returns a newly allocated array with the same contents as
!!     the linked list
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - sll: the list to trasform in array
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
function cp_sll_int_to_array(sll,error) result(res)
  type(cp_sll_int_type), pointer :: sll
  type(cp_error_type), optional, intent(inout) :: error
  integer, dimension(:), pointer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sll_int_to_array',&
        routineP=moduleN//':'//routineN
  integer :: len,stat,i
  type(cp_sll_int_type), pointer :: iter
  failure=.false.
  
  len=cp_sll_int_get_length(sll,error)
  allocate(res(len),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     iter => sll
     do i=1,len
        res(len) = iter%first_el
        if (.not.(cp_sll_int_next(iter,error=error).or.i==len)) then
           CPAssert(.false.,cp_failure_level,routineP,error,failure)
        end if
     end do
  end if
end function cp_sll_int_to_array
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_from_array [1.0] *
!!
!!   NAME
!!     cp_sll_int_from_array
!!
!!   FUNCTION
!!     returns a linked list with the same contents as the given array
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - array: the array you want to copy
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
function cp_sll_int_from_array(array,error) result(res)
  integer, dimension(:), intent(in) :: array
  type(cp_error_type), optional, intent(inout) :: error
  type(cp_sll_int_type), pointer :: res
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sll_int_from_array',&
        routineP=moduleN//':'//routineN
  integer :: i
  type(cp_sll_int_type), pointer :: last_el
  failure=.false.
  
  nullify(res,last_el)
  if (size(array)>0) then
     call cp_sll_int_create(res,&
          first_el=array(1),&
          error=error)
     last_el => res
  end if
  do i=2,size(array)
     call cp_sll_int_create(last_el%rest,&
          first_el=array(i),&
          error=error)
     last_el => last_el%rest
  end do
end function cp_sll_int_from_array
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_insert_ordered [1.0] *
!!
!!   NAME
!!     cp_sll_int_insert_ordered
!!
!!   FUNCTION
!!     insert an element mantaining the order
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - sll: the list you want to insert in
!!     - el: the element that you want to insert
!!     - insert_equals: if equal elements should be inserted
!!       (defaults to false)
!!     - did_insert: true if it did insert the element
!!     - pos: node where the element has been inserted (or of the same
!!       element that was already in the list)
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
subroutine cp_sll_int_insert_ordered(sll,el,insert_equals,&
     did_insert,pos,error)
  type(cp_sll_int_type), pointer :: sll
  integer, intent(in) :: el
  logical, optional, intent(in) :: insert_equals
  logical, optional, intent(out) :: did_insert
  type(cp_sll_int_type), pointer, optional :: pos
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure,i_eq
  character(len=*), parameter :: routineN='cp_sll_int_insert_ordered',&
        routineP=moduleN//':'//routineN
  type(cp_sll_int_type), pointer :: iter
  failure=.false.
  i_eq=.false.
  
  if (present(did_insert)) did_insert=.false.
  if (present(pos)) nullify(pos)
  if (present(insert_equals)) i_eq=insert_equals
  if (.not.associated(sll)) then
     call cp_sll_int_create(sll,first_el=el,error=error)
     if (present(did_insert)) did_insert=.true.
     if (present(pos)) pos=>sll
  else if (.not.CP_SLL_I_LESS_Q(sll%first_el,el)) then
     if (present(pos)) pos=>sll
     if (i_eq.or..not.CP_SLL_I_EQUAL_Q(el,sll%first_el)) then
        call cp_sll_int_insert_el(sll,el,error=error)
        if (present(did_insert)) did_insert=.true.
        if (present(pos)) pos=>sll
     end if
  else
     iter => sll
     do
        if (.not.associated(iter%rest)) then
           call cp_sll_int_insert_el(iter%rest,el,error=error)
           exit
        else if (.not.CP_SLL_I_LESS_Q(iter%rest%first_el,el)) then
           if (present(pos)) pos=>iter
           if (i_eq.or. .not. CP_SLL_I_EQUAL_Q(el,sll%first_el)) then
              call cp_sll_int_insert_el(iter%rest,el,error=error)
              if (present(did_insert)) did_insert=.true.
              if (present(pos)) pos=>iter%rest
           end if
           exit
        end if
        CPInvariant(cp_sll_int_next(iter),cp_failure_level,routineP,error,failure)
     end do
     CPAssert(associated(iter),cp_failure_level,routineP,error,failure)
  end if
end subroutine cp_sll_int_insert_ordered
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_insert_ordered2 [1.0] *
!!
!!   NAME
!!     cp_sll_int_insert_ordered2
!!
!!   FUNCTION
!!     insert an element mantaining the order defined by a user given function
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - sll: the list you want to insert in
!!     - el: the element that you want to insert
!!     - compare_function: a function used to compare two elements,
!!       it should return a negative number if the first argument is less
!!       than the second, 0 if they are equal, otherwise a positive number
!!     - insert_equals: if equal elements should be inserted
!!       (defaults to false)
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
subroutine cp_sll_int_insert_ordered2(sll,el,compare_function,&
     insert_equals,error)
  type(cp_sll_int_type), pointer :: sll
  integer, intent(in) :: el
  interface
     function compare_function(el1,el2)
       use cp_array_i_utils, only: cp_i_write
       integer :: compare_function
       integer, intent(in) :: el1,el2
     end function compare_function
  end interface
  logical, optional, intent(in) :: insert_equals
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure,i_eq
  character(len=*), parameter :: routineN='cp_sll_int_insert_ordered2',&
        routineP=moduleN//':'//routineN
  integer :: comp
  type(cp_sll_int_type), pointer :: iter
  failure=.false.
  i_eq=.false.
  
  if (present(insert_equals)) i_eq=insert_equals
  if (.not.associated(sll)) then
     call cp_sll_int_create(sll,first_el=el,error=error)
     return
  end if
  comp=compare_function(sll%first_el,el)
  if (comp>=0) then
     if (i_eq.or.comp/=0) then
        call cp_sll_int_insert_el(sll,el,error=error)
     end if
  else
     iter => sll
     do
        if (.not.associated(iter%rest)) then
           call cp_sll_int_insert_el(iter%rest,el,error=error)
           exit
        end if
        comp=compare_function(iter%rest%first_el,el)
        if (comp>=0) then
           if (i_eq.or. comp/=0) then
              call cp_sll_int_insert_el(iter%rest,el,error=error)
           end if
           exit
        end if
        CPInvariant(cp_sll_int_next(iter),cp_failure_level,routineP,error,failure)
     end do
     CPAssert(associated(iter),cp_failure_level,routineP,error,failure)
  end if
end subroutine cp_sll_int_insert_ordered2
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_contains [1.0] *
!!
!!   NAME
!!     cp_sll_int_contains
!!
!!   FUNCTION
!!     returns true if the list contains the given element
!!
!!   NOTES
!!     marginally faster if ordered (n/2 instread of n)
!!
!!   INPUTS
!!     - sll: the linked list you want to control
!!     - el: the element to control
!!     - ordered: if true assumes that the list is ordered in growing order.
!!       (defaults to false)
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
function cp_sll_int_contains(sll,el,ordered,error) result(res)
  logical ::res
  type(cp_sll_int_type), pointer :: sll
  integer, intent(in) :: el
  logical, intent(in), optional :: ordered
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='cp_sll_int_contains',&
        routineP=moduleN//':'//routineN
  type(cp_sll_int_type), pointer :: iter
  integer :: el_att
  failure=.false.
  res=.false.
  
  iter => sll
  if (ordered) then
     do
        if (.not.cp_sll_int_next(iter,el_att=el_att,error=error)) exit
        if (.not.CP_SLL_I_LESS_Q(el_att,el)) then
           res=CP_SLL_I_EQUAL_Q(el_att,el)
           return
        end if
     end do
  else
     do
        if (.not.cp_sll_int_next(iter,el_att=el_att,error=error)) exit
        if (CP_SLL_I_EQUAL_Q(el_att,el)) then
           res=.true.
           return
        end if
     end do
  end if
end function cp_sll_int_contains
!***************************************************************************

!!****f* cp_linked_list_int/cp_sll_int_contains [1.0] *
!!
!!   NAME
!!     cp_sll_int_contains
!!
!!   FUNCTION
!!     returns true if the list contains the given element
!!
!!   NOTES
!!     marginally faster if ordered (n/2 instread of n)
!!
!!   INPUTS
!!     - sll: the linked list you want to control
!!     - el: the element to control
!!     - compare_function: a function used to compare two elements,
!!       it should return a negative number if the first argument is less
!!       than the second, 0 if they are equal, otherwise a positive number
!!     - ordered: if true assumes that the list is ordered in growing order.
!!       (defaults to false)
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
function cp_sll_int_contains2(sll,el,compare_function,ordered,error)&
     result(res)
  logical ::res
  type(cp_sll_int_type), pointer :: sll
  integer, intent(in) :: el
  interface
     function compare_function(el1,el2)
       use cp_array_i_utils, only: cp_i_write
       integer :: compare_function
       integer, intent(in) :: el1,el2
     end function compare_function
  end interface
  logical, intent(in), optional :: ordered
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  integer :: comp
  character(len=*), parameter :: routineN='cp_sll_int_contains2',&
        routineP=moduleN//':'//routineN
  type(cp_sll_int_type), pointer :: iter
  integer :: el_att
  failure=.false.
  res=.false.
  
  iter => sll
  if (ordered) then
     do
        if (.not.cp_sll_int_next(iter,el_att=el_att,error=error)) exit
        comp=compare_function(el_att,el)
        if (comp>=0) then
           res= comp==0
           return
        end if
     end do
  else
     do
        if (.not.cp_sll_int_next(iter,el_att=el_att,error=error)) exit
        if (compare_function(el_att,el)==0) then
           res=.true.
           return
        end if
     end do
  end if
end function cp_sll_int_contains2
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same (multi-line use change it a bit)
! [template(defines,nametype1,type1,type1in,type1out,type1arrayEl,arrayEl,array=,=,use,write_el,lessQ,equalQ)]
! ARGS:
!  type1arrayEl = "integer"
!  lessQ = "CP_SLL_I_LESS_Q"
!  type1 = "integer"
!  nametype1 = "int"
!  type1out = "integer, intent(out)"
!  defines = 
!    "#define CP_SLL_I_LESS_Q(el1,el2) ( el1 < el2 )
!     #define CP_SLL_I_EQUAL_Q(el1,el2) ( el1 == el2 )
!     "
!  = = "="
!  type1in = "integer, intent(in)"
!  equalQ = "CP_SLL_I_EQUAL_Q"
!  write_el = "cp_i_write"
!  array= = "="
!  arrayEl = ""
!  use = "use cp_array_i_utils, only: cp_i_write"


end module cp_linked_list_int

