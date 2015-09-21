!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#define CP_SLL_L_LESS_Q(el1,el2) ( .not.el1.and.el2 )
#define CP_SLL_L_EQUAL_Q(el1,el2) ( el1.EQV.el2 )


! *****************************************************************************
!> \brief describes a generic linked list template.
!>      Linked list are supposed to always use pointers to the nodes for
!>      basically everything: a pointer to the node is a list, an element of
!>      the list, an iterator between the elment of the list.
!>      An empty list is represented by an unassociated pointer.
!> \note
!>     ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _
!>    |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!>    | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!>    | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!>    |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>                                                        |___/
!>      ____ _                  ___                              _ _       _       _
!>     / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _
!>    | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!>    | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!>     \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!>                                                                                   |___/
!>     _____ _     _       _____ _ _      _
!>    |_   _| |__ (_)___  |  ___(_) | ___| |
!>      | | | '_ \| / __| | |_  | | |/ _ \ |
!>      | | | | | | \__ \ |  _| | | |  __/_|
!>      |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>
!>      This is a template
!>
!>      **** DO NOT MODIFY THE .F FILES ****
!>      modify the .template and .instantition instead
!> \par History
!>      1.2002 created
!>      4.2002 changed most of the methods, by making access to the list
!>             always through pointers (identifying pointer, list and iterators)
!>      6.2004 removed %initialized from list elements
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_linked_list_logical
   
#include "../common/cp_common_uses.f90"
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_linked_list_logical'

!API type
  PUBLIC :: cp_sll_logical_type, cp_sll_logical_p_type
!API common methods
  PUBLIC :: cp_create, cp_dealloc, cp_get, cp_set, cp_next
!API special get
  PUBLIC :: cp_get_first_el, cp_get_rest, cp_get_empty, cp_get_length,&
       cp_get_element_at, cp_to_array
!API special set
  PUBLIC :: cp_set_element_at
!API structure manipulation
  PUBLIC :: cp_insert, cp_remove_first_el, cp_remove_el, cp_remove_all
!API low level
  PUBLIC :: cp_dealloc_node

  !API underlying routines
  PUBLIC :: cp_sll_logical_create, cp_sll_logical_dealloc,&
       cp_sll_logical_dealloc_node,cp_sll_logical_set,&
       cp_sll_logical_get,cp_sll_logical_next,&
       cp_sll_logical_get_first_el, cp_sll_logical_get_rest,&
       cp_sll_logical_get_empty, cp_sll_logical_get_length,&
       cp_sll_logical_get_el_at, cp_sll_logical_set_el_at,&
       cp_sll_logical_insert_el, cp_sll_logical_insert_el_at,&
       cp_sll_logical_rm_first_el, cp_sll_logical_rm_el_at,&
       cp_sll_logical_rm_all_el, &
       cp_sll_logical_to_array,&
       cp_sll_logical_from_array, cp_sll_logical_insert_ordered,&
       cp_sll_logical_insert_ordered2

! creation of an object (from a pointer)
  INTERFACE cp_create
     MODULE PROCEDURE cp_sll_logical_create
  END INTERFACE
! destruction of an object (from a pointer)
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_sll_logical_dealloc
  END INTERFACE
! destruction only of the node (low level)
  INTERFACE cp_dealloc_node
     MODULE PROCEDURE cp_sll_logical_dealloc_node
  END INTERFACE
! modifies attributes of an object
  INTERFACE cp_set
     MODULE PROCEDURE cp_sll_logical_set
  END INTERFACE
! returns attributes of an object
  INTERFACE cp_get
     MODULE PROCEDURE cp_sll_logical_get
  END INTERFACE
! iterates to the next element
  INTERFACE cp_next
     MODULE PROCEDURE cp_sll_logical_next
  END INTERFACE
! returns the first element
  INTERFACE cp_get_first_el
     MODULE PROCEDURE cp_sll_logical_get_first_el
  END INTERFACE
! returns the rest of the list
  INTERFACE cp_get_rest
     MODULE PROCEDURE cp_sll_logical_get_rest
  END INTERFACE
! returns if the list is empty
  INTERFACE cp_get_empty
     MODULE PROCEDURE cp_sll_logical_get_empty
  END INTERFACE
! returns the length of the list
  INTERFACE cp_get_length
     MODULE PROCEDURE cp_sll_logical_get_length
  END INTERFACE
! returns the element at the given position
  INTERFACE cp_get_element_at
     MODULE PROCEDURE cp_sll_logical_get_el_at
  END INTERFACE
! sets the element at the given position
  INTERFACE cp_set_element_at
     MODULE PROCEDURE cp_sll_logical_set_el_at
  END INTERFACE
! inserts one element call cp_insert(list,element,...)
  INTERFACE cp_insert
     MODULE PROCEDURE cp_sll_logical_insert_el
  END INTERFACE
!MK  INTERFACE cp_insert_ordered
!MK     MODULE PROCEDURE cp_sll_logical_insert_ordered,&
!MK          cp_sll_logical_insert_ordered2
!MK  END INTERFACE
  INTERFACE cp_insert_at
     MODULE PROCEDURE cp_sll_logical_insert_el_at
  END INTERFACE
! removes an element
  INTERFACE cp_remove_el
     MODULE PROCEDURE cp_sll_logical_rm_first_el, &
          cp_sll_logical_rm_el_at
  END INTERFACE
! removes the first el
  INTERFACE cp_remove_first_el
     MODULE PROCEDURE cp_sll_logical_rm_first_el
  END INTERFACE
! remove all the elments
  INTERFACE cp_remove_all
     MODULE PROCEDURE cp_sll_logical_rm_all_el
  END INTERFACE
! transorms the list in array
  INTERFACE cp_to_array
     MODULE PROCEDURE cp_sll_logical_to_array
  END INTERFACE

! *****************************************************************************
!> \brief represent a single linked list that stores pointers to the elements
!> \param first_el the element that is stored in this node.
!> \param rest the rest of the list
!>
!> \param empty true if the list pointer is not associated, if it points to
!>             to a not it is always false (as there is at least the
!>             first_el in the list)
!> \param length the number of elements in the list
!> \note
!>      List are alway accessed through pointers, so every node of the
!>      linked list can be seen as a list, its first element
!>      a pointer to the position before itself, in a very natural way:
!>      all the insertions take place before the actual element, and
!>      you still can insert an element at the end.
!>      This way I could identify nodes, lists and pointers between the
!>      elements of the list.
!>      Indexing is 1 based.
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_sll_logical_type
     LOGICAL :: first_el
     TYPE(cp_sll_logical_type), POINTER :: rest
  END TYPE cp_sll_logical_type

! *****************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_sll_logical_p_type
     TYPE(cp_sll_logical_type), POINTER :: list
  END TYPE cp_sll_logical_p_type

CONTAINS


! =========== creation / distruction ========

! *****************************************************************************
!> \brief allocates and initializes a single linked list
!> \param sll the single linked list to initialize
!> \param first_el the first element of this list
!> \param rest the following elements (if not given: empty)
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_sll_logical_create(sll,first_el,rest)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(in), OPTIONAL            :: first_el
    TYPE(cp_sll_logical_type), OPTIONAL, &
      POINTER                                :: rest

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_create', &
      routineP = moduleN//':'//routineN

    IF (.NOT.PRESENT(first_el)) THEN
       NULLIFY(sll)
       IF (PRESENT(rest)) sll => rest
    ELSE
       ALLOCATE(sll)
       sll%first_el = first_el
       NULLIFY(sll%rest)
       IF (PRESENT(rest)) sll%rest => rest
    END IF
  END SUBROUTINE cp_sll_logical_create

! *****************************************************************************
!> \brief deallocates the singly linked list starting at sll.
!>      Does not work if loops are present!
!> \param sll the list to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      does not deallocate the elments that are stored in the list
!>      check more?
! *****************************************************************************
  SUBROUTINE cp_sll_logical_dealloc(sll)
    TYPE(cp_sll_logical_type), POINTER       :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_dealloc', &
      routineP = moduleN//':'//routineN

    CALL cp_sll_logical_rm_all_el(sll)
  END SUBROUTINE cp_sll_logical_dealloc

! * low-level *

! *****************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_sll_logical_dealloc_node(sll)
    TYPE(cp_sll_logical_type), POINTER       :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_dealloc_node', &
      routineP = moduleN//':'//routineN

    DEALLOCATE(sll)
  END SUBROUTINE cp_sll_logical_dealloc_node

! ============= get/set ============

! *****************************************************************************
!> \brief sets tha various attributes of a single linked list
!>      It is an error to set the rest of an empty list without setting
!>      the first_el (change and make it simply a sll=>rest ?)
!> \param sll the single linked list to change
!> \param first_el the element to replace the first element of this list
!> \param rest the rest of the list (can be unassociated)
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_sll_logical_set(sll,first_el,rest)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(in), OPTIONAL            :: first_el
    TYPE(cp_sll_logical_type), OPTIONAL, &
      POINTER                                :: rest

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_set', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(sll)) THEN
       IF (PRESENT(first_el)) THEN
          CALL cp_sll_logical_create(sll,first_el,rest)
       ELSE
          CPASSERT(.NOT.PRESENT(rest))
       END IF
    ELSE
       IF (PRESENT(first_el)) sll%first_el = first_el
       IF (PRESENT(rest)) sll%rest => rest
    END IF
  END SUBROUTINE cp_sll_logical_set

! *****************************************************************************
!> \brief returns various attributes from the linked list
!> \param sll the single linked list to change
!> \param first_el the first element of this list (can be unassociated)
!> \param rest the rest of the list (can be unassociated)
!> \param empty ...
!> \param length the length of the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_sll_logical_get(sll,first_el,rest,empty,length)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(out), OPTIONAL           :: first_el
    TYPE(cp_sll_logical_type), OPTIONAL, &
      POINTER                                :: rest
    LOGICAL, INTENT(out), OPTIONAL           :: empty
    INTEGER, INTENT(out), OPTIONAL           :: length

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_get', &
      routineP = moduleN//':'//routineN

    IF (.NOT.ASSOCIATED(sll)) THEN
       CPASSERT(.NOT.PRESENT(first_el))
       IF (PRESENT(rest)) NULLIFY(rest)
       IF (PRESENT(empty)) empty=.TRUE.
       IF (PRESENT(length)) length=0
    ELSE
       IF (PRESENT(first_el)) first_el = sll%first_el
       IF (PRESENT(rest)) rest => sll%rest
       IF (PRESENT(empty)) empty = .FALSE.
       IF (PRESENT(length)) &
            length = cp_sll_logical_get_length(sll)
    END IF
  END SUBROUTINE cp_sll_logical_get

! *****************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_sll_logical_get_first_el(sll) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_get_first_el', &
      routineP = moduleN//':'//routineN

    IF (cp_debug) THEN
       CPASSERT(ASSOCIATED(sll))
    END IF

    res = sll%first_el
  END FUNCTION cp_sll_logical_get_first_el

! *****************************************************************************
!> \brief returns the rest of the list
!> \param sll the single linked list to get the rest from
!> \param iter how many times the call to rest should be iterated,
!>            defaults to 1; -1 means till end of the list.
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      split the case iter=1 to make it more optimized?
! *****************************************************************************
  FUNCTION cp_sll_logical_get_rest(sll, iter) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    INTEGER, OPTIONAL                        :: iter
    TYPE(cp_sll_logical_type), POINTER       :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_get_rest', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

    IF (.NOT.ASSOCIATED(sll)) THEN
       NULLIFY(res)
    ELSE
       IF (PRESENT(iter)) THEN
          res => sll
          DO i=1,iter
             IF (ASSOCIATED(res%rest)) THEN
                res => res%rest
             ELSE
                CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                     routineP, "tried to go past end in "// &
                     CPSourceFileRef)
             END IF
          END DO
          IF (iter==-1) THEN
             DO
                IF (.NOT.ASSOCIATED(res%rest)) EXIT
                res => res%rest
             END DO
          END IF
       ELSE
          res => sll%rest ! make the common case fast...
       END IF
    END IF
  END FUNCTION cp_sll_logical_get_rest

! *****************************************************************************
!> \brief ...
!> \param sll the single linked list to get the rest from
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_sll_logical_get_empty(sll) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_get_empty', &
      routineP = moduleN//':'//routineN

    res = .NOT.ASSOCIATED(sll)
  END FUNCTION cp_sll_logical_get_empty

! *****************************************************************************
!> \brief returns the length of the list
!> \param sll the list you want to know the length of
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(n))
! *****************************************************************************
  FUNCTION cp_sll_logical_get_length(sll) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_get_length', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: iterator

    res=0
    iterator => sll
    DO
       IF (ASSOCIATED(iterator)) THEN
          res=res+1
          iterator => iterator%rest
       ELSE
          EXIT
       END IF
    END DO
  END FUNCTION cp_sll_logical_get_length

! *****************************************************************************
!> \brief returns the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! *****************************************************************************
  FUNCTION cp_sll_logical_get_el_at(sll,index) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: index
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_get_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: pos

    IF (cp_debug) THEN
       CPASSERT(index>0.OR.index==-1)
    END IF
    IF (index==-1) THEN
       pos => cp_sll_logical_get_rest(sll, iter=-1)
    ELSE
       pos => cp_sll_logical_get_rest(sll, iter=index-1)
    END IF
    CPASSERT(ASSOCIATED(pos))

    res = pos%first_el
  END FUNCTION cp_sll_logical_get_el_at

! *****************************************************************************
!> \brief sets the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!>             -1 means at the end
!> \param value the new element
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! *****************************************************************************
  SUBROUTINE cp_sll_logical_set_el_at(sll,index,value)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: index
    LOGICAL, INTENT(in)                      :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_set_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: pos

    IF (index==-1) THEN
       pos => cp_sll_logical_get_rest(sll, iter=-1)
    ELSE
       pos => cp_sll_logical_get_rest(sll, iter=index-1)
    END IF
    CPASSERT(ASSOCIATED(pos))

    pos%first_el = value
  END SUBROUTINE cp_sll_logical_set_el_at

! * iteration *

! *****************************************************************************
!> \brief returns true if the actual element is valid (i.e. iterator ont at end)
!>      moves the iterator to the next element
!> \param iterator iterator that moves along the list
!> \param el_att the actual element (valid only if the function returns true)
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_sll_logical_next(iterator,el_att) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: iterator
    LOGICAL, INTENT(out), OPTIONAL           :: el_att
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_next', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(iterator)) THEN
       res=.TRUE.
       IF (PRESENT(el_att)) el_att = iterator%first_el
       iterator => iterator%rest
    ELSE
       res=.FALSE.
    END IF
  END FUNCTION cp_sll_logical_next

! ============ structure modifications ============

! *****************************************************************************
!> \brief insert an element at the beginning of the list
!> \param sll the single linked list point at the beginning of which
!>           you want to add the element
!> \param el the element to add
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! *****************************************************************************
  SUBROUTINE cp_sll_logical_insert_el(sll,el)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(in)                      :: el

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_insert_el', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: newSlot

    NULLIFY(newSlot)

    CALL cp_sll_logical_create(newSlot,first_el=el,&
         rest=sll)
    sll => newSlot
  END SUBROUTINE cp_sll_logical_insert_el

! *****************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! *****************************************************************************
  SUBROUTINE cp_sll_logical_rm_first_el(sll)
    TYPE(cp_sll_logical_type), POINTER       :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_rm_first_el', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: node_to_rm

    node_to_rm => sll

    IF (ASSOCIATED(sll)) THEN
       sll => sll%rest
       CALL cp_sll_logical_dealloc_node(node_to_rm)
    ELSE
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
            routineP,"tried to remove first el of an empty list in "//&
            CPSourceFileRef)
    END IF
  END SUBROUTINE cp_sll_logical_rm_first_el

! *****************************************************************************
!> \brief inserts the element at the given index
!> \param sll the list you get the element from
!> \param el the new element
!> \param index the position of the element (stating at 1).
!>             If it is -1, it means at end
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! *****************************************************************************
  SUBROUTINE cp_sll_logical_insert_el_at(sll,el,index)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(in)                      :: el
    INTEGER, INTENT(in)                      :: index

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_insert_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: pos

    IF (index==1) THEN
       CALL cp_sll_logical_insert_el(sll,el)
    ELSE
       IF (index==-1) THEN
          pos => cp_sll_logical_get_rest(sll, iter=-1)
       ELSE
          pos => cp_sll_logical_get_rest(sll, iter=index-2)
       END IF
       CPASSERT(ASSOCIATED(pos))
       CALL cp_sll_logical_insert_el(pos%rest,el)
    END IF
  END SUBROUTINE cp_sll_logical_insert_el_at

! *****************************************************************************
!> \brief removes the element at the given index
!> \param sll the list you get the element from
!> \param index the position of the element (stating at 1)
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      slow (O(index))
! *****************************************************************************
  SUBROUTINE cp_sll_logical_rm_el_at(sll,index)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: index

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_rm_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: pos

    IF (cp_debug) THEN
       CPASSERT(index>0.OR.index==-1)
    END IF
    IF (index==1) THEN
       CALL cp_sll_logical_rm_first_el(sll)
    ELSE
       IF (index==-1) THEN
          pos => cp_sll_logical_get_rest(sll, iter=-1)
       ELSE
          pos => cp_sll_logical_get_rest(sll, iter=index-2)
       END IF
       CPASSERT(ASSOCIATED(pos))
       CALL cp_sll_logical_rm_first_el(pos%rest)
    END IF
  END SUBROUTINE cp_sll_logical_rm_el_at

! *****************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! *****************************************************************************
  SUBROUTINE cp_sll_logical_rm_all_el(sll)
    TYPE(cp_sll_logical_type), POINTER       :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_rm_all_el', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: actual_node, next_node

    actual_node => sll
    DO
       IF (.NOT.ASSOCIATED(actual_node)) EXIT
       next_node => actual_node%rest
       CALL cp_sll_logical_dealloc_node(actual_node)
       actual_node => next_node
    END DO
    NULLIFY(sll)
  END SUBROUTINE cp_sll_logical_rm_all_el

! *****************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to trasform in array
!> \retval res ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION cp_sll_logical_to_array(sll) RESULT(res)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, DIMENSION(:), POINTER           :: res

    INTEGER                                  :: i, len
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_to_array', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_logical_type), POINTER       :: iter

  len=cp_sll_logical_get_length(sll)
  ALLOCATE(res(len))
  iter => sll
  DO i=1,len
     res(i) = iter%first_el
     CPASSERT(cp_sll_logical_next(iter).OR.i==len)
  END DO
END FUNCTION cp_sll_logical_to_array

! *****************************************************************************
!> \brief returns a linked list with the same contents as the given array
!> \param array the array you want to copy
!> \retval res ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION cp_sll_logical_from_array(array) RESULT(res)
    LOGICAL, DIMENSION(:), INTENT(in)        :: array
    TYPE(cp_sll_logical_type), POINTER       :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_logical_from_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    TYPE(cp_sll_logical_type), POINTER       :: last_el

  NULLIFY(res,last_el)
  IF (SIZE(array)>0) THEN
     CALL cp_sll_logical_create(res,&
          first_el=array(1))
     last_el => res
  END IF
  DO i=2,SIZE(array)
     CALL cp_sll_logical_create(last_el%rest,&
          first_el=array(i))
     last_el => last_el%rest
  END DO
END FUNCTION cp_sll_logical_from_array

! *****************************************************************************
!> \brief insert an element mantaining the order
!> \param sll the list you want to insert in
!> \param el the element that you want to insert
!> \param insert_equals if equal elements should be inserted
!>        (defaults to false)
!> \param did_insert true if it did insert the element
!> \param pos node where the element has been inserted (or of the same
!>        element that was already in the list)
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_sll_logical_insert_ordered(sll,el,insert_equals,&
     did_insert,pos)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(in)                      :: el
    LOGICAL, INTENT(in), OPTIONAL            :: insert_equals
    LOGICAL, INTENT(out), OPTIONAL           :: did_insert
    TYPE(cp_sll_logical_type), OPTIONAL, &
      POINTER                                :: pos

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_sll_logical_insert_ordered', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: i_eq
    TYPE(cp_sll_logical_type), POINTER       :: iter

  i_eq=.FALSE.

  IF (PRESENT(did_insert)) did_insert=.FALSE.
  IF (PRESENT(pos)) NULLIFY(pos)

  IF (PRESENT(insert_equals)) i_eq=insert_equals
  IF (.NOT.ASSOCIATED(sll)) THEN
     CALL cp_sll_logical_create(sll,first_el=el)
     IF (PRESENT(did_insert)) did_insert=.TRUE.
     IF (PRESENT(pos)) pos=>sll
  ELSE IF (.NOT.CP_SLL_L_LESS_Q(sll%first_el,el)) THEN
     IF (PRESENT(pos)) pos=>sll
     IF (i_eq.OR.CP_SLL_L_LESS_Q(el,sll%first_el)) THEN
        CALL cp_sll_logical_insert_el(sll,el)
        IF (PRESENT(did_insert)) did_insert=.TRUE.
        IF (PRESENT(pos)) pos=>sll
     END IF
  ELSE
     iter => sll
     DO
        IF (.NOT.ASSOCIATED(iter%rest)) THEN
           CALL cp_sll_logical_insert_el(iter%rest,el)
           IF (PRESENT(did_insert)) did_insert=.TRUE.
           IF (PRESENT(pos)) pos=>iter%rest
           EXIT
        ELSE IF (.NOT.CP_SLL_L_LESS_Q(iter%rest%first_el,el)) THEN
           IF (PRESENT(pos)) pos=>iter
           IF (i_eq.OR. CP_SLL_L_LESS_Q(el,iter%rest%first_el)) THEN
              CALL cp_sll_logical_insert_el(iter%rest,el)
              IF (PRESENT(did_insert)) did_insert=.TRUE.
              IF (PRESENT(pos)) pos=>iter%rest
           END IF
           EXIT
        END IF
        CPASSERT(cp_sll_logical_next(iter))
     END DO
     CPASSERT(ASSOCIATED(iter))
  END IF
END SUBROUTINE cp_sll_logical_insert_ordered

! *****************************************************************************
!> \brief insert an element mantaining the order defined by a user given function
!> \param sll the list you want to insert in
!> \param el the element that you want to insert
!> \param compare_function a function used to compare two elements,
!>        it should return a negative number if the first argument is less
!>        than the second, 0 if they are equal, otherwise a positive number
!> \param insert_equals if equal elements should be inserted
!>        (defaults to false)
!> \param did_insert ...
!> \param pos ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE cp_sll_logical_insert_ordered2(sll,el,compare_function,&
     insert_equals,did_insert,pos)
    TYPE(cp_sll_logical_type), POINTER       :: sll
    LOGICAL, INTENT(in)                      :: el
  INTERFACE
! *****************************************************************************
     FUNCTION compare_function(el1,el2)
        
       INTEGER :: compare_function
       LOGICAL, INTENT(in) :: el1,el2
     END FUNCTION compare_function
  END INTERFACE
    LOGICAL, INTENT(in), OPTIONAL            :: insert_equals
    LOGICAL, INTENT(out), OPTIONAL           :: did_insert
    TYPE(cp_sll_logical_type), OPTIONAL, &
      POINTER                                :: pos

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_sll_logical_insert_ordered2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comp
    LOGICAL                                  :: i_eq
    TYPE(cp_sll_logical_type), POINTER       :: iter

  i_eq=.FALSE.

  IF (PRESENT(did_insert)) did_insert=.FALSE.
  IF (PRESENT(pos)) NULLIFY(pos)

  IF (PRESENT(insert_equals)) i_eq=insert_equals
  IF (.NOT.ASSOCIATED(sll)) THEN
     CALL cp_sll_logical_create(sll,first_el=el)
     IF (PRESENT(did_insert)) did_insert=.TRUE.
     IF (PRESENT(pos)) pos=>sll%rest
     RETURN
  END IF
  comp=compare_function(sll%first_el,el)
  IF (comp>=0) THEN
     IF (i_eq.OR.comp/=0) THEN
        CALL cp_sll_logical_insert_el(sll,el)
        IF (PRESENT(did_insert)) did_insert=.TRUE.
        IF (PRESENT(pos)) pos=>sll%rest
     END IF
  ELSE
     iter => sll
     DO
        IF (.NOT.ASSOCIATED(iter%rest)) THEN
           CALL cp_sll_logical_insert_el(iter%rest,el)
           IF (PRESENT(did_insert)) did_insert=.TRUE.
           IF (PRESENT(pos)) pos=>iter%rest
           EXIT
        END IF
        comp=compare_function(iter%rest%first_el,el)
        IF (comp>=0) THEN
           IF (i_eq.OR. comp/=0) THEN
              CALL cp_sll_logical_insert_el(iter%rest,el)
              IF (PRESENT(did_insert)) did_insert=.TRUE.
              IF (PRESENT(pos)) pos=>iter%rest
           END IF
           EXIT
        END IF
        CPASSERT(cp_sll_logical_next(iter))
     END DO
     CPASSERT(ASSOCIATED(iter))
  END IF
END SUBROUTINE cp_sll_logical_insert_ordered2

! template def put here so that line numbers in template and derived
! files are almost the same (multi-line use change it a bit)
! [template(defines,nametype1,type1,type1in,type1out,type1arrayEl,arrayEl,array=,=,USE,write_el,lessQ,equalQ,private_routines)]
! ARGS:
!  = = "="
!  USE = " "
!  array= = "="
!  arrayEl = ""
!  common_dir = "../common"
!  defines = 
!    "#define CP_SLL_L_LESS_Q(el1,el2) ( .not.el1.and.el2 )
!     #define CP_SLL_L_EQUAL_Q(el1,el2) ( el1.EQV.el2 )
!     "
!  equalQ = "CP_SLL_L_EQUAL_Q"
!  lessQ = "CP_SLL_L_LESS_Q"
!  nametype1 = "logical"
!  private_routines = ""
!  template = "../common/cp_linked_list__nametype1_.template"
!  type1 = "logical"
!  type1arrayEl = "logical"
!  type1in = "logical, intent(in)"
!  type1out = "logical, intent(out)"
!  write_el = ""


END MODULE cp_linked_list_logical

