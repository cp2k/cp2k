!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"
#define CP_SLL_I_LESS_Q(el1,el2) ( el1 < el2 )
#define CP_SLL_I_EQUAL_Q(el1,el2) ( el1 == el2 )
#define PW_WRITE_UNIMPLEMENTED(el, unit, err)

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
!!     This is a template
!!
!!     **** DO NOT MODIFY THE .F FILES ****
!!     modify the .template and .instantition files
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
MODULE cp_linked_list_int_ptr
  USE cp_array_i_utils,                ONLY: cp_i_write
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_output_handling,              ONLY: cp_finished_output,&
                                             cp_unitnr_for_output,&
                                             cp_would_output
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_linked_list_int_ptr'

! type
  PUBLIC :: cp_sll_int_ptr_type
! common methods
  PUBLIC :: cp_create, cp_dealloc, cp_get, cp_set, cp_next, cp_output
! special get
  PUBLIC :: cp_get_first_el, cp_get_rest, cp_get_empty, cp_get_length,&
       cp_get_element_at
! special set
  PUBLIC :: cp_set_element_at
! structure manipulation
  PUBLIC :: cp_insert, cp_remove_first_el, cp_remove_el, cp_remove_all
! low level
  PUBLIC :: cp_dealloc_node, cp_write

  ! underlying routines
  PUBLIC :: cp_sll_int_ptr_create, cp_sll_int_ptr_dealloc,&
       cp_sll_int_ptr_dealloc_node,cp_sll_int_ptr_set,&
       cp_sll_int_ptr_get,cp_sll_int_ptr_next,&
       cp_sll_int_ptr_get_first_el, cp_sll_int_ptr_get_rest,&
       cp_sll_int_ptr_get_empty, cp_sll_int_ptr_get_length,&
       cp_sll_int_ptr_get_el_at, cp_sll_int_ptr_set_el_at,&
       cp_sll_int_ptr_insert_el, cp_sll_int_ptr_insert_el_at,&
       cp_sll_int_ptr_rm_first_el, cp_sll_int_ptr_rm_el_at,&
       cp_sll_int_ptr_rm_all_el, cp_sll_int_ptr_output,&
       cp_sll_int_ptr_write, cp_sll_int_ptr_to_array,&
       cp_sll_int_ptr_from_array, cp_sll_int_ptr_insert_ordered,&
       cp_sll_int_ptr_insert_ordered2

! creation of an object (from a pointer)
  INTERFACE cp_create
     MODULE PROCEDURE cp_sll_int_ptr_create
  END INTERFACE
! destruction of an object (from a pointer)
  INTERFACE cp_dealloc
     MODULE PROCEDURE cp_sll_int_ptr_dealloc
  END INTERFACE
! destruction only of the node (low level)
  INTERFACE cp_dealloc_node
     MODULE PROCEDURE cp_sll_int_ptr_dealloc_node
  END INTERFACE
! modifies attributes of an object
  INTERFACE cp_set
     MODULE PROCEDURE cp_sll_int_ptr_set
  END INTERFACE
! returns attributes of an object
  INTERFACE cp_get
     MODULE PROCEDURE cp_sll_int_ptr_get
  END INTERFACE
! iterates to the next element
  INTERFACE cp_next
     MODULE PROCEDURE cp_sll_int_ptr_next
  END INTERFACE
! returns the first element
  INTERFACE cp_get_first_el
     MODULE PROCEDURE cp_sll_int_ptr_get_first_el
  END INTERFACE
! returns the rest of the list
  INTERFACE cp_get_rest
     MODULE PROCEDURE cp_sll_int_ptr_get_rest
  END INTERFACE
! returns if the list is empty
  INTERFACE cp_get_empty
     MODULE PROCEDURE cp_sll_int_ptr_get_empty
  END INTERFACE
! returns the length of the list
  INTERFACE cp_get_length
     MODULE PROCEDURE cp_sll_int_ptr_get_length
  END INTERFACE
! returns the element at the given position
  INTERFACE cp_get_element_at
     MODULE PROCEDURE cp_sll_int_ptr_get_el_at
  END INTERFACE
! sets the element at the given position
  INTERFACE cp_set_element_at
     MODULE PROCEDURE cp_sll_int_ptr_set_el_at
  END INTERFACE
! inserts one element call cp_insert(list,element,...)
  INTERFACE cp_insert
     MODULE PROCEDURE cp_sll_int_ptr_insert_el
  END INTERFACE
  INTERFACE cp_insert_ordered
     MODULE PROCEDURE cp_sll_int_ptr_insert_ordered,&
          cp_sll_int_ptr_insert_ordered2
  END INTERFACE
  INTERFACE cp_insert_at
     MODULE PROCEDURE cp_sll_int_ptr_insert_el_at
  END INTERFACE
! removes an element
  INTERFACE cp_remove_el
     MODULE PROCEDURE cp_sll_int_ptr_rm_first_el, &
          cp_sll_int_ptr_rm_el_at
  END INTERFACE
! removes the first el
  INTERFACE cp_remove_first_el
     MODULE PROCEDURE cp_sll_int_ptr_rm_first_el
  END INTERFACE
! remove all the elments
  INTERFACE cp_remove_all
     MODULE PROCEDURE cp_sll_int_ptr_rm_all_el
  END INTERFACE
! outputs the list
  INTERFACE cp_output
     MODULE PROCEDURE cp_sll_int_ptr_output
  END INTERFACE
! writes the list (low-level)
  INTERFACE cp_write
     MODULE PROCEDURE cp_sll_int_ptr_write
  END INTERFACE
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
  TYPE cp_sll_int_ptr_type
     LOGICAL :: initialized
     INTEGER, POINTER :: first_el
     TYPE(cp_sll_int_ptr_type), POINTER :: rest
  END TYPE cp_sll_int_ptr_type
!!***
!****************************************************************************

CONTAINS


! =========== creation / distruction ========

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_create [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_create
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
  SUBROUTINE cp_sll_int_ptr_create(sll,first_el,rest,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), OPTIONAL, TARGET    :: first_el
    TYPE(cp_sll_int_ptr_type), OPTIONAL, &
      POINTER                                :: rest
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_create'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.PRESENT(first_el)) THEN
       NULLIFY(sll)
       IF (PRESENT(rest)) sll => rest
    ELSE
       ALLOCATE(sll, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          sll%first_el => first_el
          NULLIFY(sll%rest)
          IF (PRESENT(rest)) sll%rest => rest
          sll%initialized=.TRUE.
       END IF
    END IF
    IF (failure) NULLIFY(sll)
  END SUBROUTINE cp_sll_int_ptr_create
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_dealloc
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
  SUBROUTINE cp_sll_int_ptr_dealloc(sll,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc'

    CALL cp_sll_int_ptr_rm_all_el(sll,error)
  END SUBROUTINE cp_sll_int_ptr_dealloc
!***************************************************************************

! * low-level *

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc_node [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_dealloc_node
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
  SUBROUTINE cp_sll_int_ptr_dealloc_node(sll,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_dealloc_node'

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
    sll%initialized=.FALSE.
    DEALLOCATE(sll, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  END SUBROUTINE cp_sll_int_ptr_dealloc_node
!***************************************************************************

! ============= get/set ============

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_set [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_set
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
  SUBROUTINE cp_sll_int_ptr_set(sll,first_el,rest,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), OPTIONAL, TARGET    :: first_el
    TYPE(cp_sll_int_ptr_type), OPTIONAL, &
      POINTER                                :: rest
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_set_sll_int_ptr'

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.ASSOCIATED(sll)) THEN
       IF (PRESENT(first_el)) THEN
          CALL cp_sll_int_ptr_create(sll,first_el,rest,error)
       ELSE
          CPAssert(.NOT.PRESENT(rest),cp_failure_level,routineP,error,failure)
       END IF
    ELSE
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          IF (PRESENT(first_el)) sll%first_el => first_el
          IF (PRESENT(rest)) sll%rest => rest
       END IF
    END IF
  END SUBROUTINE cp_sll_int_ptr_set
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get
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
  SUBROUTINE cp_sll_int_ptr_get(sll,first_el,rest,empty,length,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, OPTIONAL, POINTER               :: first_el
    TYPE(cp_sll_int_ptr_type), OPTIONAL, &
      POINTER                                :: rest
    LOGICAL, INTENT(out), OPTIONAL           :: empty
    INTEGER, INTENT(out), OPTIONAL           :: length
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_get'

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.ASSOCIATED(sll)) THEN
       CPPrecondition(.NOT.PRESENT(first_el),cp_failure_level,routineP,error,failure)
       IF (PRESENT(rest)) NULLIFY(rest)
       IF (PRESENT(empty)) empty=.TRUE.
       IF (PRESENT(length)) length=0
    ELSE
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          IF (PRESENT(first_el)) first_el => sll%first_el
          IF (PRESENT(rest)) rest => sll%rest
          IF (PRESENT(empty)) empty = .FALSE.
          IF (PRESENT(length)) &
               length = cp_sll_int_ptr_get_length(sll,error=error)
       END IF
    END IF
  END SUBROUTINE cp_sll_int_ptr_get
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_first_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_first_el
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
  FUNCTION cp_sll_int_ptr_get_first_el(sll,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, POINTER                         :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_get_first'

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (cp_debug) THEN
       CPPrecondition(ASSOCIATED(sll),cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    IF (.NOT. failure) THEN
       res => sll%first_el
    END IF
  END FUNCTION cp_sll_int_ptr_get_first_el
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_rest [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_rest
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
  FUNCTION cp_sll_int_ptr_get_rest(sll, iter, error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, OPTIONAL                        :: iter
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_sll_int_ptr_type), POINTER       :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_get_first'

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.ASSOCIATED(sll)) THEN
       NULLIFY(res)
    ELSE
       CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
       IF (.NOT. failure) THEN
          IF (PRESENT(iter)) THEN
             res => sll
             DO i=1,iter
                IF (ASSOCIATED(res%rest)) THEN
                   res => res%rest
                ELSE
                   CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                        routineP, "tried to go past end in "// &
CPSourceFileRef,&
                        error,failure)
                END IF
             END DO
             IF (iter==-1) THEN
                res => sll
                DO i=1,iter
                   IF (ASSOCIATED(res%rest)) THEN
                      res => res%rest
                   ELSE
                      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
                           routineP, "tried to go past end in "// &
CPSourceFileRef,&
                           error,failure)
                   END IF
                END DO
             END IF
          ELSE
             res => sll%rest ! make the common case fast...
          END IF
       ELSE
          NULLIFY(res)
       END IF
    END IF
  END FUNCTION cp_sll_int_ptr_get_rest
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_empty [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_empty
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
  FUNCTION cp_sll_int_ptr_get_empty(sll,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_get_empty'

    LOGICAL                                  :: failure

    failure=.FALSE.

    res = .NOT.ASSOCIATED(sll)
  END FUNCTION cp_sll_int_ptr_get_empty
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_length [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_length
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
  FUNCTION cp_sll_int_ptr_get_length(sll,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_get_length'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: iterator

    failure=.FALSE.

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
  END FUNCTION cp_sll_int_ptr_get_length
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_get_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_get_el_at
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
  FUNCTION cp_sll_int_ptr_get_el_at(sll,index,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: index
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, POINTER                         :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_get_el_at'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: pos

    failure=.FALSE.

    IF (cp_debug) THEN
       CPPrecondition(index>0.OR.index==-1,cp_failure_level,routineP,error,failure)
    END IF
    IF (index==-1) THEN
       pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
    ELSE
       pos => cp_sll_int_ptr_get_rest(sll, iter=index-1,error=error)
    END IF
    CPPrecondition(ASSOCIATED(pos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       res => pos%first_el
    END IF
  END FUNCTION cp_sll_int_ptr_get_el_at
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_set_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_set_el_at
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
  SUBROUTINE cp_sll_int_ptr_set_el_at(sll,index,value,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: index
    INTEGER, INTENT(in), TARGET              :: value
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_set_el_at'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: pos

    failure=.FALSE.

    IF (index==-1) THEN
       pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
    ELSE
       pos => cp_sll_int_ptr_get_rest(sll, iter=index-1,error=error)
    END IF
    CPPrecondition(ASSOCIATED(pos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pos%first_el => value
    END IF
  END SUBROUTINE cp_sll_int_ptr_set_el_at
!***************************************************************************

! * iteration *

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_next [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_next
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
  FUNCTION cp_sll_int_ptr_next(iterator,el_att,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: iterator
    INTEGER, OPTIONAL, POINTER               :: el_att
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_next'

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(iterator)) THEN
       IF (cp_debug) THEN
          CPPrecondition(iterator%initialized,cp_failure_level,routineP,error,failure)
       END IF
       res=.NOT.failure
       IF (PRESENT(el_att)) el_att => iterator%first_el
       iterator => iterator%rest
    ELSE
       res=.FALSE.
    END IF
  END FUNCTION cp_sll_int_ptr_next
!***************************************************************************

! ============ structure modifications ============

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_insert_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_insert_el
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
!!     el: the element to add
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
  SUBROUTINE cp_sll_int_ptr_insert_el(sll,el,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), TARGET              :: el
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr:cp_sll_int_ptr_add_el'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: newSlot

    failure=.FALSE.
    NULLIFY(newSlot)

    CALL cp_sll_int_ptr_create(newSlot,first_el=el,&
         rest=sll,error=error)
    sll => newSlot
  END SUBROUTINE cp_sll_int_ptr_insert_el
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_rm_first_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_rm_first_el
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
  SUBROUTINE cp_sll_int_ptr_rm_first_el(sll,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_rm_first_el'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: node_to_rm

    failure=.FALSE.
    node_to_rm => sll

    IF (ASSOCIATED(sll)) THEN
       sll => sll%rest
       CALL cp_sll_int_ptr_dealloc_node(node_to_rm,error=error)
    ELSE
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
            routineP,"tried to remove first el of an empty list in "//&
CPSourceFileRef,&
            error,failure)
    END IF
  END SUBROUTINE cp_sll_int_ptr_rm_first_el
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_insert_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_insert_el_at
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
!!     el: the new element
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
  SUBROUTINE cp_sll_int_ptr_insert_el_at(sll,el,index,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), TARGET              :: el
    INTEGER, INTENT(in)                      :: index
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_insert_el_at'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: pos

    failure=.FALSE.

    IF (index==1) THEN
       CALL cp_sll_int_ptr_insert_el(sll,el,error=error)
    ELSE
       IF (index==-1) THEN
          pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
       ELSE
          pos => cp_sll_int_ptr_get_rest(sll, iter=index-1,error=error)
       END IF
       CPPrecondition(ASSOCIATED(pos),cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          CALL cp_sll_int_ptr_insert_el(pos%rest,el,error=error)
       END IF
    END IF
  END SUBROUTINE cp_sll_int_ptr_insert_el_at
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_rm_el_at [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_rm_el_at
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
!!     el: the new element
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
  SUBROUTINE cp_sll_int_ptr_rm_el_at(sll,index,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: index
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_remove_el_at'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: pos

    failure=.FALSE.
    
    IF (cp_debug) THEN
       CPPrecondition(index>0.OR.index==-1,cp_failure_level,routineP,error,failure)
    END IF
    IF (index==1) THEN
       CALL cp_sll_int_ptr_rm_first_el(sll,error=error)
    ELSE
       IF (index==-1) THEN
          pos => cp_sll_int_ptr_get_rest(sll, iter=-1,error=error)
       ELSE
          pos => cp_sll_int_ptr_get_rest(sll, iter=index-2,error=error)
       END IF
       CPPrecondition(ASSOCIATED(pos),cp_failure_level,routineP,error,failure)
       IF (.NOT.failure) THEN
          CALL cp_sll_int_ptr_rm_first_el(pos%rest,error=error)
       END IF
    END IF
  END SUBROUTINE cp_sll_int_ptr_rm_el_at
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_rm_all_el [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_rm_all_el
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
  SUBROUTINE cp_sll_int_ptr_rm_all_el(sll,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineP = 'cp_linked_list_int_ptr/cp_sll_int_ptr_rm_all_el'

    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: actual_node, next_node

    failure=.FALSE.

    actual_node => sll
    DO
       IF (.NOT.ASSOCIATED(actual_node)) EXIT
       next_node => cp_get_rest(actual_node,error=error)
       CALL cp_sll_int_ptr_dealloc_node(actual_node,error=error)
       actual_node => next_node
    END DO
    NULLIFY(sll)
  END SUBROUTINE cp_sll_int_ptr_rm_all_el
!***************************************************************************

! ========= output =========

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_output [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_output
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
  SUBROUTINE cp_sll_int_ptr_output(logger, outputName,&
       fromWhere, iter, sll,comment,local, finish, error)
    TYPE(cp_logger_type), POINTER            :: logger
    CHARACTER(len=*), INTENT(in)             :: outputName, fromWhere
    INTEGER, INTENT(in)                      :: iter
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    CHARACTER(len=*), INTENT(in), OPTIONAL   :: comment
    LOGICAL, INTENT(in), OPTIONAL            :: local, finish
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_int_ptr_output', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iostat, unitNr
    LOGICAL                                  :: failure, fin, loc

    failure=.FALSE.; loc=.FALSE.; fin=.TRUE.
  
    CALL timeset(routineN//','//moduleN,'I',"",handle)
    IF (PRESENT(local)) loc=local
    IF (PRESENT(finish)) fin=finish

    IF (loc .OR. logger%para_env%mepos==logger%para_env%source) THEN
       IF (cp_would_output(logger, outputName , fromWhere, iter,&
            local=loc,error=error)) THEN

          unitNr=cp_unitnr_for_output(logger, outputName, fromWhere, iter,&
               local=loc,error=error)
          
          IF (PRESENT(comment)) THEN
             WRITE(unit=unitNr,fmt='(" ",a,a)',iostat=iostat) '# ',comment
             CPInvariant(iostat==0,cp_failure_level,routineP,error,failure)
          END IF
          CALL cp_sll_int_ptr_write(sll,unit_nr=unitNr,error=error)

          IF (fin) THEN
             CALL cp_finished_output(logger, outputName, fromWhere, iter,&
                  unit_nr=unitNr,local=local,error=error)
          END IF
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_sll_int_ptr_output
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_write [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_write
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
  SUBROUTINE cp_sll_int_ptr_write(sll,unit_nr,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_int_ptr_write', &
      routineP = moduleN//':'//routineN

    INTEGER, POINTER                         :: el
    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: iterator

    failure=.FALSE.
    iterator => sll
  
    DO
       IF (.NOT.cp_sll_int_ptr_next(iterator,el_att=el,error=error)) EXIT
       CALL cp_i_write(el, unit_nr=unit_nr, error=error)
    END DO
  END SUBROUTINE cp_sll_int_ptr_write
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_to_array [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_to_array
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
FUNCTION cp_sll_int_ptr_to_array(sll,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER, DIMENSION(:), POINTER           :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_int_ptr_to_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, len, stat
    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: iter

  failure=.FALSE.
  
  len=cp_sll_int_ptr_get_length(sll,error)
  ALLOCATE(res(len),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     iter => sll
     DO i=1,len
        res(len) = iter%first_el
        IF (.NOT.(cp_sll_int_ptr_next(iter,error=error).OR.i==len)) THEN
           CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
     END DO
  END IF
END FUNCTION cp_sll_int_ptr_to_array
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_from_array [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_from_array
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
FUNCTION cp_sll_int_ptr_from_array(array,error) RESULT(res)
    INTEGER, DIMENSION(:), INTENT(in)        :: array
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cp_sll_int_ptr_type), POINTER       :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_int_ptr_from_array', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: last_el

  failure=.FALSE.
  
  NULLIFY(res,last_el)
  IF (SIZE(array)>0) THEN
     CALL cp_sll_int_ptr_create(res,&
          first_el=array(1),&
          error=error)
     last_el => res
  END IF
  DO i=2,SIZE(array)
     CALL cp_sll_int_ptr_create(last_el%rest,&
          first_el=array(i),&
          error=error)
     last_el => last_el%rest
  END DO
END FUNCTION cp_sll_int_ptr_from_array
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_insert_ordered [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_insert_ordered
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
SUBROUTINE cp_sll_int_ptr_insert_ordered(sll,el,insert_equals,&
     did_insert,pos,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), TARGET              :: el
    LOGICAL, INTENT(in), OPTIONAL            :: insert_equals
    LOGICAL, INTENT(out), OPTIONAL           :: did_insert
    TYPE(cp_sll_int_ptr_type), OPTIONAL, &
      POINTER                                :: pos
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_sll_int_ptr_insert_ordered', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, i_eq
    TYPE(cp_sll_int_ptr_type), POINTER       :: iter

  failure=.FALSE.
  i_eq=.FALSE.
  
  IF (PRESENT(did_insert)) did_insert=.FALSE.
  IF (PRESENT(pos)) NULLIFY(pos)

  IF (PRESENT(insert_equals)) i_eq=insert_equals
  IF (.NOT.ASSOCIATED(sll)) THEN
     CALL cp_sll_int_ptr_create(sll,first_el=el,error=error)
     IF (PRESENT(did_insert)) did_insert=.TRUE.
     IF (PRESENT(pos)) pos=>sll
  ELSE IF (.NOT.CP_SLL_I_LESS_Q(sll%first_el,el)) THEN
     CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
     IF (PRESENT(pos)) pos=>sll
     IF (i_eq.OR.CP_SLL_I_LESS_Q(el,sll%first_el)) THEN
        CALL cp_sll_int_ptr_insert_el(sll,el,error=error)
        IF (PRESENT(did_insert)) did_insert=.TRUE.
        IF (PRESENT(pos)) pos=>sll
     END IF
  ELSE
     CPPrecondition(sll%initialized,cp_failure_level,routineP,error,failure)
     iter => sll
     DO
        IF (.NOT.ASSOCIATED(iter%rest)) THEN
           CALL cp_sll_int_ptr_insert_el(iter%rest,el,error=error)
           IF (PRESENT(did_insert)) did_insert=.TRUE.
           IF (PRESENT(pos)) pos=>iter%rest
           EXIT
        ELSE IF (.NOT.CP_SLL_I_LESS_Q(iter%rest%first_el,el)) THEN
           IF (PRESENT(pos)) pos=>iter
           IF (i_eq.OR. CP_SLL_I_LESS_Q(el,iter%rest%first_el)) THEN
              CALL cp_sll_int_ptr_insert_el(iter%rest,el,error=error)
              IF (PRESENT(did_insert)) did_insert=.TRUE.
              IF (PRESENT(pos)) pos=>iter%rest
           END IF
           EXIT
        END IF
        CPInvariant(cp_sll_int_ptr_next(iter),cp_failure_level,routineP,error,failure)
     END DO
     CPAssert(ASSOCIATED(iter),cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE cp_sll_int_ptr_insert_ordered
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_insert_ordered2 [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_insert_ordered2
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
SUBROUTINE cp_sll_int_ptr_insert_ordered2(sll,el,compare_function,&
     insert_equals,did_insert,pos,error)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), TARGET              :: el
  INTERFACE
     FUNCTION compare_function(el1,el2)
       USE kinds, ONLY: dp
USE cp_array_i_utils, ONLY: cp_i_write
       INTEGER :: compare_function
       INTEGER, INTENT(in),TARGET :: el1,el2
     END FUNCTION compare_function
  END INTERFACE
    LOGICAL, INTENT(in), OPTIONAL            :: insert_equals
    LOGICAL, INTENT(out), OPTIONAL           :: did_insert
    TYPE(cp_sll_int_ptr_type), OPTIONAL, &
      POINTER                                :: pos
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'cp_sll_int_ptr_insert_ordered2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comp
    LOGICAL                                  :: failure, i_eq
    TYPE(cp_sll_int_ptr_type), POINTER       :: iter

  failure=.FALSE.
  i_eq=.FALSE.
  
  IF (PRESENT(did_insert)) did_insert=.FALSE.
  IF (PRESENT(pos)) NULLIFY(pos)

  IF (PRESENT(insert_equals)) i_eq=insert_equals
  IF (.NOT.ASSOCIATED(sll)) THEN
     CALL cp_sll_int_ptr_create(sll,first_el=el,error=error)
     IF (PRESENT(did_insert)) did_insert=.TRUE.
     IF (PRESENT(pos)) pos=>iter%rest
     RETURN
  END IF
  comp=compare_function(sll%first_el,el)
  IF (comp>=0) THEN
     IF (i_eq.OR.comp/=0) THEN
        CALL cp_sll_int_ptr_insert_el(sll,el,error=error)
        IF (PRESENT(did_insert)) did_insert=.TRUE.
        IF (PRESENT(pos)) pos=>iter%rest
     END IF
  ELSE
     iter => sll
     DO
        IF (.NOT.ASSOCIATED(iter%rest)) THEN
           CALL cp_sll_int_ptr_insert_el(iter%rest,el,error=error)
           IF (PRESENT(did_insert)) did_insert=.TRUE.
           IF (PRESENT(pos)) pos=>iter%rest
           EXIT
        END IF
        comp=compare_function(iter%rest%first_el,el)
        IF (comp>=0) THEN
           IF (i_eq.OR. comp/=0) THEN
              CALL cp_sll_int_ptr_insert_el(iter%rest,el,error=error)
              IF (PRESENT(did_insert)) did_insert=.TRUE.
              IF (PRESENT(pos)) pos=>iter%rest
           END IF
           EXIT
        END IF
        CPInvariant(cp_sll_int_ptr_next(iter),cp_failure_level,routineP,error,failure)
     END DO
     CPAssert(ASSOCIATED(iter),cp_failure_level,routineP,error,failure)
  END IF
END SUBROUTINE cp_sll_int_ptr_insert_ordered2
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_contains [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_contains
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
FUNCTION cp_sll_int_ptr_contains(sll,el,ordered,error) RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), TARGET              :: el
    LOGICAL, INTENT(in), OPTIONAL            :: ordered
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_int_ptr_contains', &
      routineP = moduleN//':'//routineN

    INTEGER, POINTER                         :: el_att
    LOGICAL                                  :: failure, ord
    TYPE(cp_sll_int_ptr_type), POINTER       :: iter

  failure=.FALSE.; ord=.FALSE.
  res=.FALSE.
  
  iter => sll
  IF (PRESENT(ordered)) ord=ordered
  IF (ord) THEN
     DO
        IF (.NOT.cp_sll_int_ptr_next(iter,el_att=el_att,error=error)) EXIT
        IF (.NOT.CP_SLL_I_LESS_Q(el_att,el)) THEN
           res=.NOT.CP_SLL_I_LESS_Q(el,el_att)
           RETURN
        END IF
     END DO
  ELSE
     DO
        IF (.NOT.cp_sll_int_ptr_next(iter,el_att=el_att,error=error)) EXIT
        IF (.NOT.CP_SLL_I_EQUAL_Q(el_att,el)) THEN
           res=.TRUE.
           RETURN
        END IF
     END DO
  END IF
END FUNCTION cp_sll_int_ptr_contains
!***************************************************************************

!!****f* cp_linked_list_int_ptr/cp_sll_int_ptr_contains [1.0] *
!!
!!   NAME
!!     cp_sll_int_ptr_contains
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
FUNCTION cp_sll_int_ptr_contains2(sll,el,compare_function,ordered,error)&
     RESULT(res)
    TYPE(cp_sll_int_ptr_type), POINTER       :: sll
    INTEGER, INTENT(in), TARGET              :: el
  INTERFACE
     FUNCTION compare_function(el1,el2)
       USE kinds, ONLY: dp
USE cp_array_i_utils, ONLY: cp_i_write
       INTEGER :: compare_function
       INTEGER, INTENT(in),TARGET :: el1,el2
     END FUNCTION compare_function
  END INTERFACE
    LOGICAL, INTENT(in), OPTIONAL            :: ordered
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL                                  :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_sll_int_ptr_contains2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: comp
    INTEGER, POINTER                         :: el_att
    LOGICAL                                  :: failure
    TYPE(cp_sll_int_ptr_type), POINTER       :: iter

  failure=.FALSE.
  res=.FALSE.
  
  iter => sll
  IF (ordered) THEN
     DO
        IF (.NOT.cp_sll_int_ptr_next(iter,el_att=el_att,error=error)) EXIT
        comp=compare_function(el_att,el)
        IF (comp>=0) THEN
           res= comp==0
           RETURN
        END IF
     END DO
  ELSE
     DO
        IF (.NOT.cp_sll_int_ptr_next(iter,el_att=el_att,error=error)) EXIT
        IF (compare_function(el_att,el)==0) THEN
           res=.TRUE.
           RETURN
        END IF
     END DO
  END IF
END FUNCTION cp_sll_int_ptr_contains2
!***************************************************************************

! template def put here so that line numbers in template and derived 
! files are almost the same (multi-line use change it a bit)
! [template(defines,nametype1,type1,type1in,type1out,type1arrayEl,arrayEl,array=,=,USE,write_el,lessQ,equalQ,private_routines)]
! ARGS:
!  = = "=>"
!  USE = 
!    "USE kinds, only: dp
!     USE cp_array_i_utils, ONLY: cp_i_write"
!  array= = "="
!  arrayEl = ""
!  defines = 
!    "#define CP_SLL_I_LESS_Q(el1,el2) ( el1 < el2 )
!     #define CP_SLL_I_EQUAL_Q(el1,el2) ( el1 == el2 )
!     #define PW_WRITE_UNIMPLEMENTED(el, unit, err)
!     "
!  equalQ = "CP_SLL_I_EQUAL_Q"
!  lessQ = "CP_SLL_I_LESS_Q"
!  nametype1 = "int_ptr"
!  private_routines = ""
!  type1 = "integer, pointer"
!  type1arrayEl = "integer"
!  type1in = "integer, intent(in),target"
!  type1out = "integer, pointer"
!  write_el = "cp_i_write"


END MODULE cp_linked_list_int_ptr

