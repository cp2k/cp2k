!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#define CP_SLL_DERIV_LESS_Q(el1,el2) cp_sll_deriv_less_q(el1,el2)
#define CP_SLL_DERIV_EQUAL_Q(el1,el2) ( el1%desc == el2%desc )


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
MODULE cp_linked_list_xc_deriv
  USE xc_derivative_types,             ONLY: xc_derivative_p_type,&
                                             xc_derivative_type
#include "../common/../base/base_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_linked_list_xc_deriv'

!API type
  PUBLIC :: cp_sll_xc_deriv_type, cp_sll_xc_deriv_p_type
!API common methods
  PUBLIC :: cp_create, cp_dealloc, cp_next
!API special get
  PUBLIC :: cp_get_first_el, cp_get_rest, cp_get_length, cp_get_element_at, cp_to_array
!API special set
  PUBLIC :: cp_set_element_at
!API structure manipulation
  PUBLIC :: cp_insert, cp_remove_first_el, cp_remove_el, cp_remove_all

  !API underlying routines
  PUBLIC :: cp_sll_xc_deriv_create, cp_sll_xc_deriv_dealloc,&
       cp_sll_xc_deriv_next,&
       cp_sll_xc_deriv_get_length,&
       cp_sll_xc_deriv_get_el_at, cp_sll_xc_deriv_set_el_at,&
       cp_sll_xc_deriv_get_first_el, cp_sll_xc_deriv_get_rest,&
       cp_sll_xc_deriv_insert_el, cp_sll_xc_deriv_insert_el_at,&
       cp_sll_xc_deriv_rm_first_el, cp_sll_xc_deriv_rm_el_at,&
       cp_sll_xc_deriv_rm_all_el, &
       cp_sll_xc_deriv_to_array

! creation of an object (from a pointer)
  INTERFACE cp_create
     MODULE PROCEDURE cp_sll_xc_deriv_create
  END INTERFACE
! destruction of an object (from a pointer)
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_sll_xc_deriv_dealloc
  END INTERFACE
! iterates to the next element
  INTERFACE cp_next
     MODULE PROCEDURE cp_sll_xc_deriv_next
  END INTERFACE
! returns the first element
  INTERFACE cp_get_first_el
     MODULE PROCEDURE cp_sll_xc_deriv_get_first_el
  END INTERFACE
! returns the rest of the list
  INTERFACE cp_get_rest
     MODULE PROCEDURE cp_sll_xc_deriv_get_rest
  END INTERFACE
! returns the length of the list
  INTERFACE cp_get_length
     MODULE PROCEDURE cp_sll_xc_deriv_get_length
  END INTERFACE
! returns the element at the given position
  INTERFACE cp_get_element_at
     MODULE PROCEDURE cp_sll_xc_deriv_get_el_at
  END INTERFACE
! sets the element at the given position
  INTERFACE cp_set_element_at
     MODULE PROCEDURE cp_sll_xc_deriv_set_el_at
  END INTERFACE
! inserts one element call cp_insert(list,element,...)
  INTERFACE cp_insert
     MODULE PROCEDURE cp_sll_xc_deriv_insert_el
  END INTERFACE
  INTERFACE cp_insert_at
     MODULE PROCEDURE cp_sll_xc_deriv_insert_el_at
  END INTERFACE
! removes an element
  INTERFACE cp_remove_el
     MODULE PROCEDURE cp_sll_xc_deriv_rm_first_el, &
          cp_sll_xc_deriv_rm_el_at
  END INTERFACE
! removes the first el
  INTERFACE cp_remove_first_el
     MODULE PROCEDURE cp_sll_xc_deriv_rm_first_el
  END INTERFACE
! remove all the elments
  INTERFACE cp_remove_all
     MODULE PROCEDURE cp_sll_xc_deriv_rm_all_el
  END INTERFACE
! transorms the list in array
  INTERFACE cp_to_array
     MODULE PROCEDURE cp_sll_xc_deriv_to_array
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
  TYPE cp_sll_xc_deriv_type
     TYPE(xc_derivative_type),POINTER :: first_el
     TYPE(cp_sll_xc_deriv_type), POINTER :: rest
  END TYPE cp_sll_xc_deriv_type

! *****************************************************************************
!> \brief pointer to a linked list (to make arrays of pointers)
!> \param list the pointer to the list
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_sll_xc_deriv_p_type
     TYPE(cp_sll_xc_deriv_type), POINTER :: list
  END TYPE cp_sll_xc_deriv_p_type

CONTAINS


! *****************************************************************************
!> \brief private compare function
!> \param el1 ...
!> \param el2 ...
!> \retval res ...
! *****************************************************************************
FUNCTION cp_sll_deriv_less_q(el1,el2) RESULT(res)
    TYPE(xc_derivative_type), POINTER        :: el1, el2
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_deriv_less_q', &
      routineP = moduleN//':'//routineN

  CPASSERT(ASSOCIATED(el1))
  CPASSERT(ASSOCIATED(el1%split_desc))
  CPASSERT(ASSOCIATED(el2))
  CPASSERT(ASSOCIATED(el2%split_desc))
  res=SIZE(el1%split_desc)<SIZE(el2%split_desc).OR.&
      (SIZE(el1%split_desc)==SIZE(el2%split_desc).and.el1%desc<el2%desc)
END FUNCTION

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
  SUBROUTINE cp_sll_xc_deriv_create(sll,first_el,rest)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    TYPE(xc_derivative_type), OPTIONAL, &
      POINTER                                :: first_el
    TYPE(cp_sll_xc_deriv_type), OPTIONAL, &
      POINTER                                :: rest

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_create', &
      routineP = moduleN//':'//routineN

    IF (.NOT.PRESENT(first_el)) THEN
       NULLIFY(sll)
       IF (PRESENT(rest)) sll => rest
    ELSE
       ALLOCATE(sll)
       sll%first_el => first_el
       NULLIFY(sll%rest)
       IF (PRESENT(rest)) sll%rest => rest
    END IF
  END SUBROUTINE cp_sll_xc_deriv_create

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
  SUBROUTINE cp_sll_xc_deriv_dealloc(sll)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_dealloc', &
      routineP = moduleN//':'//routineN

    CALL cp_sll_xc_deriv_rm_all_el(sll)
  END SUBROUTINE cp_sll_xc_deriv_dealloc

! * low-level *

! *****************************************************************************
!> \brief deallocates a node of a singly linked list (low level)
!> \param sll the node to be deallocated
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_sll_xc_deriv_dealloc_node(sll)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_dealloc_node', &
      routineP = moduleN//':'//routineN

    DEALLOCATE(sll)
  END SUBROUTINE cp_sll_xc_deriv_dealloc_node

! ============= get/set ============

! *****************************************************************************
!> \brief returns the first element stored in the list
!> \param sll the single linked list to get the element from
!> \retval res ...
!> \par History
!>      none
!> \author Fawzi Mohamed
! *****************************************************************************
  FUNCTION cp_sll_xc_deriv_get_first_el(sll) RESULT(res)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    TYPE(xc_derivative_type), POINTER        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_get_first_el', &
      routineP = moduleN//':'//routineN

    res => sll%first_el
  END FUNCTION cp_sll_xc_deriv_get_first_el

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
  FUNCTION cp_sll_xc_deriv_get_rest(sll, iter) RESULT(res)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    INTEGER, OPTIONAL                        :: iter
    TYPE(cp_sll_xc_deriv_type), POINTER      :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_get_rest', &
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
                CPABORT("tried to go past end")
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
  END FUNCTION cp_sll_xc_deriv_get_rest

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
  FUNCTION cp_sll_xc_deriv_get_length(sll) RESULT(res)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_get_length', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: iterator

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
  END FUNCTION cp_sll_xc_deriv_get_length

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
  FUNCTION cp_sll_xc_deriv_get_el_at(sll,index) RESULT(res)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    INTEGER, INTENT(in)                      :: index
    TYPE(xc_derivative_type), POINTER        :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_get_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos

    IF (index==-1) THEN
       pos => cp_sll_xc_deriv_get_rest(sll, iter=-1)
    ELSE
       pos => cp_sll_xc_deriv_get_rest(sll, iter=index-1)
    END IF
    CPASSERT(ASSOCIATED(pos))

    res => pos%first_el
  END FUNCTION cp_sll_xc_deriv_get_el_at

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
  SUBROUTINE cp_sll_xc_deriv_set_el_at(sll,index,value)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    INTEGER, INTENT(in)                      :: index
    TYPE(xc_derivative_type), POINTER        :: value

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_set_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos

    IF (index==-1) THEN
       pos => cp_sll_xc_deriv_get_rest(sll, iter=-1)
    ELSE
       pos => cp_sll_xc_deriv_get_rest(sll, iter=index-1)
    END IF
    CPASSERT(ASSOCIATED(pos))

    pos%first_el => value
  END SUBROUTINE cp_sll_xc_deriv_set_el_at

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
  FUNCTION cp_sll_xc_deriv_next(iterator,el_att) RESULT(res)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: iterator
    TYPE(xc_derivative_type), OPTIONAL, &
      POINTER                                :: el_att
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_next', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(iterator)) THEN
       res=.TRUE.
       IF (PRESENT(el_att)) el_att => iterator%first_el
       iterator => iterator%rest
    ELSE
       res=.FALSE.
    END IF
  END FUNCTION cp_sll_xc_deriv_next

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
  SUBROUTINE cp_sll_xc_deriv_insert_el(sll,el)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    TYPE(xc_derivative_type), POINTER        :: el

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_insert_el', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: newSlot

    NULLIFY(newSlot)

    CALL cp_sll_xc_deriv_create(newSlot,first_el=el,&
         rest=sll)
    sll => newSlot
  END SUBROUTINE cp_sll_xc_deriv_insert_el

! *****************************************************************************
!> \brief remove the first element of the linked list
!> \param sll the list whose first element has to be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      fast (O(1))
! *****************************************************************************
  SUBROUTINE cp_sll_xc_deriv_rm_first_el(sll)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_rm_first_el', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: node_to_rm

    node_to_rm => sll

    IF (ASSOCIATED(sll)) THEN
       sll => sll%rest
       CALL cp_sll_xc_deriv_dealloc_node(node_to_rm)
    ELSE
       CPABORT("tried to remove first el of an empty list")
    END IF
  END SUBROUTINE cp_sll_xc_deriv_rm_first_el

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
  SUBROUTINE cp_sll_xc_deriv_insert_el_at(sll,el,index)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    TYPE(xc_derivative_type), POINTER        :: el
    INTEGER, INTENT(in)                      :: index

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_insert_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos

    IF (index==1) THEN
       CALL cp_sll_xc_deriv_insert_el(sll,el)
    ELSE
       IF (index==-1) THEN
          pos => cp_sll_xc_deriv_get_rest(sll, iter=-1)
       ELSE
          pos => cp_sll_xc_deriv_get_rest(sll, iter=index-2)
       END IF
       CPASSERT(ASSOCIATED(pos))
       CALL cp_sll_xc_deriv_insert_el(pos%rest,el)
    END IF
  END SUBROUTINE cp_sll_xc_deriv_insert_el_at

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
  SUBROUTINE cp_sll_xc_deriv_rm_el_at(sll,index)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    INTEGER, INTENT(in)                      :: index

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_rm_el_at', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: pos

    IF (index==1) THEN
       CALL cp_sll_xc_deriv_rm_first_el(sll)
    ELSE
       IF (index==-1) THEN
          pos => cp_sll_xc_deriv_get_rest(sll, iter=-1)
       ELSE
          pos => cp_sll_xc_deriv_get_rest(sll, iter=index-2)
       END IF
       CPASSERT(ASSOCIATED(pos))
       CALL cp_sll_xc_deriv_rm_first_el(pos%rest)
    END IF
  END SUBROUTINE cp_sll_xc_deriv_rm_el_at

! *****************************************************************************
!> \brief removes all the elements from the list
!> \param sll the list that should be removed
!> \par History
!>      none
!> \author Fawzi Mohamed
!> \note
!>      check more?
! *****************************************************************************
  SUBROUTINE cp_sll_xc_deriv_rm_all_el(sll)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_rm_all_el', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: actual_node, next_node

    actual_node => sll
    DO
       IF (.NOT.ASSOCIATED(actual_node)) EXIT
       next_node => actual_node%rest
       CALL cp_sll_xc_deriv_dealloc_node(actual_node)
       actual_node => next_node
    END DO
    NULLIFY(sll)
  END SUBROUTINE cp_sll_xc_deriv_rm_all_el

! *****************************************************************************
!> \brief returns a newly allocated array with the same contents as
!>      the linked list
!> \param sll the list to trasform in array
!> \retval res ...
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
FUNCTION cp_sll_xc_deriv_to_array(sll) RESULT(res)
    TYPE(cp_sll_xc_deriv_type), POINTER      :: sll
    TYPE(xc_derivative_p_type), &
      DIMENSION(:), POINTER                  :: res

    INTEGER                                  :: i, len
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_xc_deriv_to_array', &
      routineP = moduleN//':'//routineN

    TYPE(cp_sll_xc_deriv_type), POINTER      :: iter

  len=cp_sll_xc_deriv_get_length(sll)
  ALLOCATE(res(len))
  iter => sll
  DO i=1,len
     res(i)%deriv => iter%first_el
     CPASSERT(cp_sll_xc_deriv_next(iter).OR.i==len)
  END DO
END FUNCTION cp_sll_xc_deriv_to_array

! template def put here so that line numbers in template and derived
! files are almost the same (multi-line use change it a bit)
! [template(defines,nametype1,type1,type1in,type1out,type1arrayEl,arrayEl,array=,=,USE,write_el,lessQ,equalQ,private_routines)]
! ARGS:
!  = = "=>"
!  USE = 
!    "USE kinds, only: dp
!     USE xc_derivative_types, only: xc_derivative_type,xc_derivative_p_type"
!  array= = "=>"
!  arrayEl = "%deriv"
!  common_dir = "../common"
!  defines = 
!    "#define CP_SLL_DERIV_LESS_Q(el1,el2) cp_sll_deriv_less_q(el1,el2)
!     #define CP_SLL_DERIV_EQUAL_Q(el1,el2) ( el1%desc == el2%desc )
!     "
!  equalQ = "CP_SLL_DERIV_EQUAL_Q"
!  lessQ = "cp_sll_deriv_less_q"
!  nametype1 = "xc_deriv"
!  private_routines = 
!    "
!     ! *****************************************************************************
!     !> \brief private compare function
!     !> \param el1 ...
!     !> \param el2 ...
!     !> \retval res ...
!     ! *****************************************************************************
!     function cp_sll_deriv_less_q(el1,el2) result(res)
!       type(xc_derivative_type), pointer :: el1,el2
!       logical :: res
!     
!       character(len=*),parameter :: routineN='cp_sll_deriv_less_q',&
!           routineP=moduleN//':'//routineN
!     
!       CPASSERT(associated(el1))
!       CPASSERT(associated(el1%split_desc))
!       CPASSERT(associated(el2))
!       CPASSERT(associated(el2%split_desc))
!       res=size(el1%split_desc)<size(el2%split_desc).or.&
!           (size(el1%split_desc)==size(el2%split_desc).and.el1%desc<el2%desc)
!     end function
!     "
!  template = "../common/cp_linked_list__nametype1_.template"
!  type1 = "type(xc_derivative_type),pointer"
!  type1arrayEl = "type(xc_derivative_p_type)"
!  type1in = "type(xc_derivative_type),pointer"
!  type1out = "type(xc_derivative_type),pointer"
!  write_el = ""


END MODULE cp_linked_list_xc_deriv

