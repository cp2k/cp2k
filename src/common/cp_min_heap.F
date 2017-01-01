!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE cp_min_heap
   USE kinds,                           ONLY: int_4,&
                                              int_8
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: cp_heap_type, keyt, valt
   PUBLIC :: cp_heap_pop, cp_heap_reset_node, cp_heap_fill
   PUBLIC :: cp_heap_new, cp_heap_release
   PUBLIC :: cp_heap_get_first, cp_heap_reset_first

   ! Sets the types
   INTEGER, PARAMETER         :: keyt = int_4
   INTEGER, PARAMETER         :: valt = int_8

   TYPE cp_heap_node
      INTEGER(KIND=keyt) :: key
      INTEGER(KIND=valt) :: value
   END TYPE cp_heap_node

   TYPE cp_heap_node_e
      TYPE(cp_heap_node) :: node
   END TYPE cp_heap_node_e

   TYPE cp_heap_type
      INTEGER :: n
      INTEGER, DIMENSION(:), POINTER           :: index
      TYPE(cp_heap_node_e), DIMENSION(:), POINTER :: nodes
   END TYPE cp_heap_type

CONTAINS

   ! Lookup functions

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \retval parent ...
! **************************************************************************************************
   ELEMENTAL FUNCTION get_parent(n) RESULT(parent)
      INTEGER, INTENT(IN)                                :: n
      INTEGER                                            :: parent

      parent = INT(n/2)
   END FUNCTION get_parent

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \retval child ...
! **************************************************************************************************
   ELEMENTAL FUNCTION get_left_child(n) RESULT(child)
      INTEGER, INTENT(IN)                                :: n
      INTEGER                                            :: child

      child = 2*n
   END FUNCTION get_left_child

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \retval child ...
! **************************************************************************************************
   ELEMENTAL FUNCTION get_right_child(n) RESULT(child)
      INTEGER, INTENT(IN)                                :: n
      INTEGER                                            :: child

      child = 2*n+1
   END FUNCTION get_right_child

! **************************************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
!> \retval value ...
! **************************************************************************************************
   ELEMENTAL FUNCTION get_value(heap, n) RESULT(value)
      TYPE(cp_heap_type), INTENT(IN)                     :: heap
      INTEGER, INTENT(IN)                                :: n
      INTEGER(KIND=valt)                                 :: value

      value = heap%nodes(n)%node%value
   END FUNCTION get_value

! **************************************************************************************************
!> \brief ...
!> \param heap ...
!> \param key ...
!> \retval value ...
! **************************************************************************************************
   ELEMENTAL FUNCTION get_value_by_key(heap, key) RESULT(value)
      TYPE(cp_heap_type), INTENT(IN)                     :: heap
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      INTEGER(KIND=valt)                                 :: value

      INTEGER                                            :: n

      n = heap%index(key)
      value = get_value(heap, n)
   END FUNCTION get_value_by_key

   ! Initialization functions

! **************************************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE cp_heap_new(heap, n)
      TYPE(cp_heap_type), INTENT(OUT)                    :: heap
      INTEGER, INTENT(IN)                                :: n

      heap%n = n
      ALLOCATE (heap%index(n))
      ALLOCATE (heap%nodes(n))
   END SUBROUTINE cp_heap_new

! **************************************************************************************************
!> \brief ...
!> \param heap ...
! **************************************************************************************************
   SUBROUTINE cp_heap_release(heap)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap

      DEALLOCATE (heap%index)
      DEALLOCATE (heap%nodes)
      heap%n = 0
   END SUBROUTINE cp_heap_release

! **************************************************************************************************
!> \brief Fill heap with given values
!> \param heap ...
!> \param values ...
! **************************************************************************************************
   SUBROUTINE cp_heap_fill(heap, values)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER(KIND=valt), DIMENSION(:), INTENT(IN)       :: values

      INTEGER                                            :: first, i, n

      n = SIZE(values)
      CPASSERT(heap%n >= n)

      DO i = 1, n
         heap%index(i) = i
         heap%nodes(i)%node%key = i
         heap%nodes(i)%node%value = values(i)
      END DO
      ! Sort from the last full subtree
      first = get_parent(n)
      DO i = first, 1, -1
         CALL bubble_down(heap, i)
      ENDDO

   END SUBROUTINE cp_heap_fill

! **************************************************************************************************
!> \brief Returns the first heap element without removing it.
!> \param heap ...
!> \param key ...
!> \param value ...
!> \param found ...
! **************************************************************************************************
   SUBROUTINE cp_heap_get_first(heap, key, value, found)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER(KIND=keyt), INTENT(OUT)                    :: key
      INTEGER(KIND=valt), INTENT(OUT)                    :: value
      LOGICAL, INTENT(OUT)                               :: found

      IF (heap%n .LT. 1) THEN
         found = .FALSE.
      ELSE
         found = .TRUE.
         key = heap%nodes(1)%node%key
         value = heap%nodes(1)%node%value
      ENDIF
   END SUBROUTINE cp_heap_get_first

! **************************************************************************************************
!> \brief Returns and removes the first heap element and rebalances
!>        the heap.
!> \param heap ...
!> \param key ...
!> \param value ...
!> \param found ...
! **************************************************************************************************
   SUBROUTINE cp_heap_pop(heap, key, value, found)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER(KIND=keyt), INTENT(OUT)                    :: key
      INTEGER(KIND=valt), INTENT(OUT)                    :: value
      LOGICAL, INTENT(OUT)                               :: found

!

      CALL cp_heap_get_first(heap, key, value, found)
      IF (found) THEN
         IF (heap%n .GT. 1) THEN
            CALL cp_heap_copy_node(heap, 1, heap%n)
            heap%n = heap%n-1
            CALL bubble_down(heap, 1)
         ELSE
            heap%n = heap%n-1
         ENDIF
      ENDIF
   END SUBROUTINE cp_heap_pop

! **************************************************************************************************
!> \brief Changes the value of the heap element with given key and
!>        rebalances the heap.
!> \param heap ...
!> \param key ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE cp_heap_reset_node(heap, key, value)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER(KIND=keyt), INTENT(IN)                     :: key
      INTEGER(KIND=valt), INTENT(IN)                     :: value

      INTEGER                                            :: n, new_pos

      CPASSERT(heap%n > 0)

      n = heap%index(key)
      CPASSERT(heap%nodes(n)%node%key == key)
      heap%nodes(n)%node%value = value
      CALL bubble_up(heap, n, new_pos)
      CALL bubble_down(heap, new_pos)
   END SUBROUTINE cp_heap_reset_node

! **************************************************************************************************
!> \brief Changes the value of the minimum heap element and rebalances the heap.
!> \param heap ...
!> \param value ...
! **************************************************************************************************
   SUBROUTINE cp_heap_reset_first(heap, value)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER(KIND=valt), INTENT(IN)                     :: value

      CPASSERT(heap%n > 0)
      heap%nodes(1)%node%value = value
      CALL bubble_down(heap, 1)
   END SUBROUTINE cp_heap_reset_first

! **************************************************************************************************
!> \brief ...
!> \param heap ...
!> \param e1 ...
!> \param e2 ...
! **************************************************************************************************
   PURE SUBROUTINE cp_heap_swap(heap, e1, e2)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER, INTENT(IN)                                :: e1, e2

      INTEGER(KIND=keyt)                                 :: key1, key2
      TYPE(cp_heap_node)                                 :: tmp_node

      key1 = heap%nodes(e1)%node%key
      key2 = heap%nodes(e2)%node%key

      tmp_node = heap%nodes(e1)%node
      heap%nodes(e1)%node = heap%nodes(e2)%node
      heap%nodes(e2)%node = tmp_node

      heap%index(key1) = e2
      heap%index(key2) = e1
   END SUBROUTINE cp_heap_swap

! **************************************************************************************************
!> \brief Sets node e1 to e2
!> \param heap ...
!> \param e1 ...
!> \param e2 ...
! **************************************************************************************************
   PURE SUBROUTINE cp_heap_copy_node(heap, e1, e2)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER, INTENT(IN)                                :: e1, e2

      INTEGER(KIND=keyt)                                 :: key1, key2

      key1 = heap%nodes(e1)%node%key
      key2 = heap%nodes(e2)%node%key

      heap%nodes(e1)%node = heap%nodes(e2)%node

      heap%index(key1) = 0
      heap%index(key2) = e1
   END SUBROUTINE cp_heap_copy_node

! **************************************************************************************************
!> \brief Balances a heap by bubbling down from the given element.
!> \param heap ...
!> \param first ...
! **************************************************************************************************
   SUBROUTINE bubble_down(heap, first)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER, INTENT(IN)                                :: first

      INTEGER                                            :: e, left_child, right_child, smallest
      INTEGER(kind=valt)                                 :: left_child_value, min_value, &
                                                            right_child_value
      LOGICAL                                            :: all_done

!
      CPASSERT(0 < first .AND. first <= heap%n)

      e = first
      all_done = .FALSE.
      ! Check whether we are finished, i.e,. whether the element to
      ! bubble down is childless.
      DO WHILE (e .LE. get_parent(heap%n) .AND. .NOT. all_done)
         ! Determines which node (current, left, or right child) has the
         ! smallest value.
         smallest = e
         min_value = get_value(heap, e)
         left_child = get_left_child(e)
         IF (left_child .LE. heap%n) THEN
            left_child_value = get_value(heap, left_child)
            IF (left_child_value .LT. min_value) THEN
               min_value = left_child_value
               smallest = left_child
            ENDIF
         ENDIF
         right_child = left_child+1
         IF (right_child .LE. heap%n) THEN
            right_child_value = get_value(heap, right_child)
            IF (right_child_value .LT. min_value) THEN
               min_value = right_child_value
               smallest = right_child
            ENDIF
         ENDIF
         !
         CALL cp_heap_swap(heap, e, smallest)
         IF (smallest .EQ. e) THEN
            all_done = .TRUE.
         ELSE
            e = smallest
         ENDIF
      ENDDO
   END SUBROUTINE bubble_down

! **************************************************************************************************
!> \brief Balances a heap by bubbling up from the given element.
!> \param heap ...
!> \param first ...
!> \param new_pos ...
! **************************************************************************************************
   SUBROUTINE bubble_up(heap, first, new_pos)
      TYPE(cp_heap_type), INTENT(INOUT)                  :: heap
      INTEGER, INTENT(IN)                                :: first
      INTEGER, INTENT(OUT)                               :: new_pos

      INTEGER                                            :: e, parent
      INTEGER(kind=valt)                                 :: my_value, parent_value
      LOGICAL                                            :: all_done

      CPASSERT(0 < first .AND. first <= heap%n)

      e = first
      all_done = .FALSE.
      IF (e .GT. 1) THEN
         my_value = get_value(heap, e)
      ENDIF
      ! Check whether we are finished, i.e,. whether the element to
      ! bubble up is an orphan.
      new_pos = e
      DO WHILE (e .GT. 1 .AND. .NOT. all_done)
         ! Switches the parent and the current element if the current
         ! element's value is greater than the parent's value.
         parent = get_parent(e)
         parent_value = get_value(heap, parent)
         IF (my_value .LT. parent_value) THEN
            CALL cp_heap_swap(heap, e, parent)
            e = parent
         ELSE
            all_done = .TRUE.
         ENDIF
      ENDDO
      new_pos = e
   END SUBROUTINE bubble_up

END MODULE cp_min_heap
