!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE min_heap

  IMPLICIT NONE

  PUBLIC :: heap_t
  PUBLIC :: heap_get_first, heap_reset_first, heap_fill
  PUBLIC :: heap_new, heap_release

  ! Sets the types to the default integer types.
  INTEGER, PARAMETER :: keyt = KIND(0)
  INTEGER, PARAMETER :: valt = KIND(0)

  TYPE heap_node
     INTEGER(KIND=keyt) :: key
     INTEGER(KIND=valt) :: value
  END TYPE heap_node

  TYPE heap_node_e
     TYPE(heap_node) :: node
  END TYPE heap_node_e

  TYPE heap_t
     INTEGER :: n
     TYPE(heap_node_e), DIMENSION(:), POINTER :: nodes
  END TYPE heap_t

  INTERFACE node_set
     MODULE PROCEDURE node_set_node, node_set_heap
  END INTERFACE

  INTERFACE node_get
     MODULE PROCEDURE node_get_node, node_get_heap
  END INTERFACE

CONTAINS

  ! Lookup functions

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \retval parent ...
! *****************************************************************************
  ELEMENTAL FUNCTION get_parent(n) RESULT (parent)
    INTEGER, INTENT(IN)                      :: n
    INTEGER                                  :: parent

    parent = INT(n/2)
  END FUNCTION get_parent

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \retval child ...
! *****************************************************************************
  ELEMENTAL FUNCTION get_left_child(n) RESULT (child)
    INTEGER, INTENT(IN)                      :: n
    INTEGER                                  :: child

    child = 2*n
  END FUNCTION get_left_child

! *****************************************************************************
!> \brief ...
!> \param n ...
!> \retval child ...
! *****************************************************************************
  ELEMENTAL FUNCTION get_right_child(n) RESULT (child)
    INTEGER, INTENT(IN)                      :: n
    INTEGER                                  :: child

    child = 2*n+1
  END FUNCTION get_right_child

! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
!> \retval value ...
! *****************************************************************************
  ELEMENTAL FUNCTION get_value(heap, n) RESULT (value)
    TYPE(heap_t), INTENT(IN)                 :: heap
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=valt)                       :: value

    value = heap%nodes(n)%node%value
  END FUNCTION get_value


  ! Initialization functions

! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_new (heap, n, error)
    TYPE(heap_t), INTENT(OUT)                :: heap
    INTEGER, INTENT(IN)                      :: n
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER                                  :: stat

!

    heap%n = n
    ALLOCATE (heap%nodes(n), stat=stat)
    error = (stat /= 0)
  END SUBROUTINE heap_new

! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_release (heap, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER                                  :: stat

!

    DEALLOCATE (heap%nodes, stat=stat)
    error = (stat /= 0)
    NULLIFY (heap%nodes)
    heap%n = 0
  END SUBROUTINE heap_release

  !> \brief Fill heap with given keys and values
  !> \param[in,out] heap  Heap to fill.
  !> \param[in] keys      Keys to store in heap
  !> \param[in] values    Values (i.e., costs) of the keys.
  !> \param[out] error    Error occurred.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param keys ...
!> \param values ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_fill (heap, keys, values, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER(KIND=keyt), DIMENSION(:), &
      INTENT(IN)                             :: keys
    INTEGER(KIND=valt), DIMENSION(:), &
      INTENT(IN)                             :: values
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER                                  :: first, i, n

!

    error = .FALSE.
    n = SIZE (keys)
    IF (SIZE (values) .NE. n .OR. heap%n .LT. n) THEN
       error = .TRUE.
    ELSE
       DO i = 1, n
       !   call node_set(heap, i, keys(i), values(i))
          heap%nodes(i)%node%key = keys(i)
          heap%nodes(i)%node%value = values(i)
       END DO
       ! Sort from the last full subtree
       first = get_parent (n)
       DO i = first, 1, -1
          CALL bubble_down (heap, i, error)
       ENDDO
    ENDIF
  END SUBROUTINE heap_fill

  ! \brief Returns the first heap element without removing it.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param key ...
!> \param value ...
!> \param found ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_get_first (heap, key, value, found, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER(KIND=keyt), INTENT(OUT)          :: key
    INTEGER(KIND=valt), INTENT(OUT)          :: value
    LOGICAL, INTENT(OUT)                     :: found, error

!

    error = .FALSE.
    IF (heap%n .LT. 1) THEN
       found = .FALSE.
    ELSE
       found = .TRUE.
       CALL node_get (heap, 1, key, value)
    ENDIF
  END SUBROUTINE heap_get_first

  ! \brief Changes the value of the minimum heap element and
  !        rebalances the heap.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param value ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_reset_first (heap, value, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER(KIND=valt), INTENT(IN)           :: value
    LOGICAL, INTENT(OUT)                     :: error

!

    error = .FALSE.
    IF (heap%n .LT. 1) THEN
       error = .TRUE.
    ELSE
       heap%nodes(1)%node%value = value
       CALL bubble_down (heap, 1, error)
    ENDIF
  END SUBROUTINE heap_reset_first

  ! \brief Returns and removes the first heap element and rebalances
  !        the heap.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param key ...
!> \param value ...
!> \param found ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_pop (heap, key, value, found, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER(KIND=keyt), INTENT(OUT)          :: key
    INTEGER(KIND=valt), INTENT(OUT)          :: value
    LOGICAL, INTENT(OUT)                     :: found, error

!

    CALL heap_get_first (heap, key, value, found, error)
    IF (.NOT. error .AND. found) THEN
       IF (heap%n .GT. 1) THEN
          CALL heap_copy_node (heap, 1, heap%n)
          heap%n = heap%n-1
          CALL bubble_down (heap, 1, error)
       ELSE
          heap%n = heap%n-1
       ENDIF
    ENDIF
  END SUBROUTINE heap_pop

  ! \brief Adds an element to the heap and rebalances it.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param key ...
!> \param value ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_add (heap, key, value, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER(KIND=keyt), INTENT(IN)           :: key
    INTEGER(KIND=valt), INTENT(IN)           :: value
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER                                  :: new_n, new_pos

!

    error = .FALSE.
    new_n = heap%n+1
    CALL heap_reallocate (heap, new_n, error)
    heap%n = new_n
    IF (error) RETURN
    CALL node_set (heap, new_n, key, value)
    CALL bubble_up (heap, new_n, new_pos, error)
    CALL bubble_down (heap, new_pos, error)
  END SUBROUTINE heap_add


! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param e1 ...
!> \param e2 ...
! *****************************************************************************
  PURE SUBROUTINE heap_swap (heap, e1, e2)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER, INTENT(IN)                      :: e1, e2

    TYPE(heap_node)                          :: tmp_node

!

    tmp_node = heap%nodes(e1)%node
    heap%nodes(e1)%node = heap%nodes(e2)%node
    heap%nodes(e2)%node = tmp_node
  END SUBROUTINE heap_swap

  ! \brief Sets node e1 to e2
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param e1 ...
!> \param e2 ...
! *****************************************************************************
  PURE SUBROUTINE heap_copy_node (heap, e1, e2)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER, INTENT(IN)                      :: e1, e2

!

    heap%nodes(e1)%node = heap%nodes(e2)%node
  END SUBROUTINE heap_copy_node


! *****************************************************************************
!> \brief ...
!> \param node ...
!> \param key ...
!> \param value ...
! *****************************************************************************
  PURE SUBROUTINE node_get_node (node, key, value)
    TYPE(heap_node), INTENT(IN)              :: node
    INTEGER(KIND=keyt), INTENT(OUT)          :: key
    INTEGER(KIND=valt), INTENT(OUT)          :: value

!

    key = node%key
    value = node%value
  END SUBROUTINE node_get_node
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
!> \param key ...
!> \param value ...
! *****************************************************************************
  PURE SUBROUTINE node_get_heap (heap, n, key, value)
    TYPE(heap_t), INTENT(IN)                 :: heap
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=keyt), INTENT(OUT)          :: key
    INTEGER(KIND=valt), INTENT(OUT)          :: value

!

    CALL node_get_node (heap%nodes(n)%node, key, value)
  END SUBROUTINE node_get_heap
! *****************************************************************************
!> \brief ...
!> \param node ...
!> \param key ...
!> \param value ...
! *****************************************************************************
  PURE SUBROUTINE node_set_node (node, key, value)
    TYPE(heap_node), INTENT(OUT)             :: node
    INTEGER(KIND=keyt), INTENT(IN)           :: key
    INTEGER(KIND=valt), INTENT(IN)           :: value

!

    node%key = key
    node%value = value
  END SUBROUTINE node_set_node
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
!> \param key ...
!> \param value ...
! *****************************************************************************
  PURE SUBROUTINE node_set_heap (heap, n, key, value)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=keyt), INTENT(IN)           :: key
    INTEGER(KIND=valt), INTENT(IN)           :: value

!

    CALL node_set_node (heap%nodes(n)%node, key, value)
  END SUBROUTINE node_set_heap


  ! \brief Balances a heap by bubbling down from the given element.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param first ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE bubble_down (heap, first, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER, INTENT(IN)                      :: first
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER                                  :: e, left_child, right_child, &
                                                smallest
    INTEGER(kind=valt)                       :: left_child_value, min_value, &
                                                right_child_value
    LOGICAL                                  :: all_done

!

    IF (first .GT. heap%n .OR. first .LT. 1) THEN
       error = .TRUE.
       RETURN
    ENDIF
    error = .FALSE.
    e = first
    all_done = .FALSE.
    ! Check whether we are finished, i.e,. whether the element to
    ! bubble down is childless.
    DO WHILE (e .LE. get_parent (heap%n) .AND. .NOT. all_done)
       ! Determines which node (current, left, or right child) has the
       ! smallest value.
       smallest = e
       min_value = get_value (heap, e)
       left_child = get_left_child (e)
       IF (left_child .LE. heap%n) THEN
          left_child_value = get_value (heap, left_child)
          IF (left_child_value .LT. min_value) THEN
             min_value = left_child_value
             smallest = left_child
          ENDIF
       ENDIF
       right_child = left_child + 1
       IF (right_child .LE. heap%n) THEN
          right_child_value = get_value (heap, right_child)
          IF (right_child_value .LT. min_value) THEN
             min_value = right_child_value
             smallest = right_child
          ENDIF
       ENDIF
       !
       CALL heap_swap (heap, e, smallest)
       IF (smallest .EQ. e) THEN
          all_done = .TRUE.
       ELSE
          e = smallest
       ENDIF
    ENDDO
  END SUBROUTINE bubble_down


  ! \brief Balances a heap by bubbling up from the given element.
! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param first ...
!> \param new_pos ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE bubble_up (heap, first, new_pos, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER, INTENT(IN)                      :: first
    INTEGER, INTENT(OUT)                     :: new_pos
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER                                  :: e, parent
    INTEGER(kind=valt)                       :: my_value, parent_value
    LOGICAL                                  :: all_done

!

    IF (first .GT. heap%n .OR. first .LT. 1) THEN
       error = .TRUE.
       RETURN
    ENDIF
    error = .FALSE.
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
       parent = get_parent (e)
       parent_value = get_value (heap, parent)
       IF (my_value .LT. parent_value) THEN
          CALL heap_swap (heap, e, parent)
          e = parent
       ELSE
          all_done = .TRUE.
       ENDIF
    ENDDO
    new_pos = e
  END SUBROUTINE bubble_up

! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param n ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_reallocate (heap, n, error)
    TYPE(heap_t), INTENT(INOUT)              :: heap
    INTEGER, INTENT(IN)                      :: n
    LOGICAL, INTENT(OUT)                     :: error

    INTEGER, PARAMETER                       :: min_size = 10
    REAL, PARAMETER                          :: oversize = 1.4

    INTEGER                                  :: n_allocate, stat
    TYPE(heap_node_e), DIMENSION(:), POINTER :: new_nodes, tmp_nodes

!

    error = .FALSE.
    IF (n .GT. SIZE (heap%nodes)) THEN
       n_allocate = MAX(MAX(n, min_size), INT(n*oversize))
       ALLOCATE (new_nodes(n_allocate), stat=stat)
       IF (stat /= 0) THEN
          ALLOCATE (new_nodes(n), stat=stat)
       ENDIF
       IF (stat /= 0) THEN
          error = .TRUE.
          RETURN
       ENDIF
       new_nodes(1:heap%n) = heap%nodes(1:heap%n)
       tmp_nodes => heap%nodes
       heap%nodes => new_nodes
       DEALLOCATE (tmp_nodes, stat=stat)
       IF (stat /= 0) error = .TRUE.
    ENDIF
  END SUBROUTINE heap_reallocate

! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_check (heap, error)
    TYPE(heap_t), INTENT(in)                 :: heap
    LOGICAL, INTENT(out)                     :: error

    INTEGER                                  :: i, n
    INTEGER(kind=valt)                       :: left_child_value, my_value, &
                                                right_child_value

!

    error = .FALSE.
    n = heap%n
    DO i = 1, n
       my_value = get_value (heap, i)
       IF (get_left_child (i) .LE. n) THEN
          left_child_value = get_value (heap, get_left_child(i))
          IF (left_child_value .LT. my_value) THEN
             error = .TRUE.
             WRITE(*,*)"Error at",i
          ENDIF
       ENDIF
       IF (get_right_child (i) .LE. n) THEN
          right_child_value = get_value (heap, get_right_child(i))
          IF (right_child_value .LT. my_value) THEN
             error = .TRUE.
             WRITE(*,*)"Error at",i
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE heap_check


! *****************************************************************************
!> \brief ...
!> \param heap ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE heap_print (heap, error)
    TYPE(heap_t), INTENT(in)                 :: heap
    LOGICAL, INTENT(out)                     :: error

    INTEGER                                  :: first, i, last, n
    INTEGER(kind=keyt), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER(kind=valt), ALLOCATABLE, &
      DIMENSION(:)                           :: values

    error = .FALSE.
    CALL heap_check (heap, error)
    IF (error) THEN
       WRITE(*,*)"Error in heap"
    ENDIF
    n = heap%n
    first = 1
    last = first*2-1
    DO WHILE (first .LE. heap%n)
       last = MIN(last, heap%n)
       ALLOCATE (keys(first:last))
       ALLOCATE (values(first:last))
       DO i = first, last
          CALL node_get (heap, i, keys(i), values(i))
       END DO
       WRITE(*,'(8(("[",I3,",",I3,"]",1X)))')&
            (/ ((/keys(i),values(i)/), i=first,last) /)
       !write(*,'("[",I3,",",I3,"]",1X)') (/ (i, i=first,last) /)
       DEALLOCATE (keys)
       DEALLOCATE (values)
       first = first*2
       last = first*2-1
    ENDDO
  END SUBROUTINE heap_print

  !

! *****************************************************************************
!> \brief ...
!> \param error ...
!> \param message ...
! *****************************************************************************
  SUBROUTINE error_test (error, message)
    LOGICAL, INTENT(in)                      :: error
    CHARACTER(len=*), INTENT(in)             :: message

    IF (error) THEN
       WRITE(*,'(1X,A,1X,A)')"Error:", message
    ENDIF
  END SUBROUTINE error_test

! *****************************************************************************
!> \brief ...
!> \param test_size ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE test_heap (test_size, error)
    INTEGER, INTENT(in), OPTIONAL            :: test_size
    LOGICAL, INTENT(out)                     :: error

    INTEGER                                  :: i, n
    INTEGER(kind=keyt)                       :: key
    INTEGER(kind=keyt), ALLOCATABLE, &
      DIMENSION(:)                           :: keys
    INTEGER(kind=valt)                       :: value
    INTEGER(kind=valt), ALLOCATABLE, &
      DIMENSION(:)                           :: values
    LOGICAL                                  :: found
    TYPE(heap_t)                             :: heap

!
!

    IF (PRESENT(test_size)) THEN
       n = test_size
    ELSE
       n = 10
    ENDIF
    CALL heap_new (heap, n, error=error)
    CALL error_test (error, "new")
    !
    ALLOCATE (keys(n))
    ALLOCATE (values(n))
    FORALL (i = 1:n)
       keys(i) = n-i
       values(i) = 10*(n-i)
    END FORALL
    CALL heap_fill (heap, keys, values, error)
    CALL error_test (error, "fill")
    !
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Popping"
    CALL heap_pop (heap, key, value, found, error)
    CALL error_test (error, "pop")
    WRITE(*,*)key, value, found
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Popping"
    CALL heap_pop (heap, key, value, found, error)
    CALL error_test (error, "pop")
    WRITE(*,*)key, value, found
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Resetting first to 55"
    CALL heap_reset_first (heap, 55, error)
    CALL error_test (error, "reset_first")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Resetting first to 1"
    CALL heap_reset_first (heap, 1, error)
    CALL error_test (error, "reset_first")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Resetting first to 333"
    CALL heap_reset_first (heap, 333, error)
    CALL error_test (error, "reset_first")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Adding -2, -20"
    CALL heap_add (heap, -2, -20, error)
    CALL error_test (error, "add")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Adding 66, 66"
    CALL heap_add (heap, 66, 66, error)
    CALL error_test (error, "add")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    found = .TRUE.
    DO WHILE (found)
       key = 0 ; value = 0
       CALL heap_pop (heap, key, value, found, error)
       CALL error_test (error, "pop")
       WRITE(*,*)key, value, found
    END DO
    !
    WRITE(*,*)" Adding -2, -20"
    CALL heap_add (heap, -2, -20, error)
    CALL error_test (error, "add")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    WRITE(*,*)" Adding 66, 66"
    CALL heap_add (heap, 66, 66, error)
    CALL error_test (error, "add")
    CALL heap_print (heap, error)
    CALL error_test (error, "print")
    !
    DEALLOCATE (keys)
    DEALLOCATE (values)
    CALL heap_release (heap, error)
    CALL error_test (error, "release")
  END SUBROUTINE test_heap

END MODULE min_heap
