!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_neighbor_list_types [1.0] *
!!
!!   NAME
!!     qs_neighbor_list_types
!!
!!   FUNCTION
!!     Define the neighbor list data types for Quickstep.
!!
!!   AUTHOR
!!     Matthias Krack (21.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_neighbor_list_types

! *****************************************************************************

  USE kinds, ONLY: int_size,&
                   wp => dp,&
                   wp_size => dp_size

  USE termination, ONLY: stop_memory,&
                         stop_program

  IMPLICIT NONE

  PRIVATE

! *** Definition of the data types for a linked list of neighbors ***

  TYPE neighbor_node_type
    PRIVATE
    TYPE(neighbor_node_type), POINTER :: next_neighbor_node
    INTEGER                           :: neighbor
    REAL(wp), DIMENSION(3)            :: r
    INTEGER, DIMENSION(3)             :: cell
  END TYPE neighbor_node_type

  TYPE neighbor_list_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: next_neighbor_list
    TYPE(neighbor_node_type), POINTER :: first_neighbor_node
    INTEGER                           :: atom,nnode
    INTEGER, DIMENSION(3)             :: cell
  END TYPE neighbor_list_type

  TYPE neighbor_list_set_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: first_neighbor_list
    REAL(wp)                          :: r_max
    INTEGER                           :: nlist
  END TYPE neighbor_list_set_type

  TYPE neighbor_node_p_type
    TYPE(neighbor_node_type), POINTER :: neighbor_node
  END TYPE neighbor_node_p_type

  TYPE neighbor_list_p_type
    TYPE(neighbor_list_type), POINTER :: neighbor_list
  END TYPE neighbor_list_p_type

  TYPE neighbor_list_set_p_type
    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
  END TYPE neighbor_list_set_p_type

! *** Public data types ***

  PUBLIC :: neighbor_list_p_type,&
            neighbor_list_set_type,&
            neighbor_list_set_p_type,&
            neighbor_list_type,&
            neighbor_node_p_type,&
            neighbor_node_type

! *** Public subroutines ***

  PUBLIC :: add_neighbor_list,&
            add_neighbor_node,&
            add_neighbor_nodes,&
            allocate_neighbor_list_set,&
            deallocate_neighbor_list_set,&
            deallocate_neighbor_nodes,&
            extract_neighbor_list,&
            get_neighbor_list,&
            get_neighbor_list_set,&
            get_neighbor_node,&
            set_neighbor_list,&
            set_neighbor_node

! *** Public functions ***

  PUBLIC :: find_neighbor_list,&
            first_list,&
            first_node,&
            next

! *****************************************************************************

  INTERFACE find_neighbor_list
    MODULE PROCEDURE find_neighbor_list_1,&
                     find_neighbor_list_2
  END INTERFACE

  INTERFACE next
    MODULE PROCEDURE next_neighbor_list_1,&
                     next_neighbor_list_2,&
                     next_neighbor_node
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE add_neighbor_list(neighbor_list_set,atom,cell,neighbor_list)

!   Purpose: Allocate a new neighbor list.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER     :: neighbor_list
    INTEGER, INTENT(IN)                   :: atom
    INTEGER, DIMENSION(3), INTENT(IN)     :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_neighbor_list (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: new_neighbor_list

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Allocate a new neighbor list ***

    ALLOCATE (new_neighbor_list,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_list",0)

!   *** Initialize the new neighbor node ***

    NULLIFY (new_neighbor_list%first_neighbor_node)
    new_neighbor_list%atom = atom
    new_neighbor_list%nnode = 0
    new_neighbor_list%cell(:) = cell(:)

!   *** Link the new neighbor list to the neighbor list set ***

    new_neighbor_list%next_neighbor_list =>&
      neighbor_list_set%first_neighbor_list

    neighbor_list_set%first_neighbor_list => new_neighbor_list

    neighbor_list => new_neighbor_list

!   *** Increment the neighbor list counter ***

    neighbor_list_set%nlist = neighbor_list_set%nlist + 1

  END SUBROUTINE add_neighbor_list

! *****************************************************************************

  SUBROUTINE add_neighbor_node(neighbor_list,neighbor,r,cell)

!   Purpose: Allocate a new neighbor list node and store the neighbor data.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER  :: neighbor_list
    INTEGER, INTENT(IN)                :: neighbor
    REAL(wp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, DIMENSION(3), INTENT(IN)  :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_neighbor_node (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: new_neighbor_node

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Allocate a new neighbor data set ***

    ALLOCATE (new_neighbor_node,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_node",0)

!   *** Store the data set of the new neighbor ***

    new_neighbor_node%neighbor = neighbor
    new_neighbor_node%r(:) = r(:)
    new_neighbor_node%cell(:) = cell(:)

!   *** Link the new neighbor node to the list ***

    new_neighbor_node%next_neighbor_node => neighbor_list%first_neighbor_node

    neighbor_list%first_neighbor_node => new_neighbor_node

!   *** Increment the neighbor node counter ***

    neighbor_list%nnode = neighbor_list%nnode + 1

  END SUBROUTINE add_neighbor_node

! *****************************************************************************

  SUBROUTINE add_neighbor_nodes(neighbor_list_set,atom,cell,nnode,neighbors,r)

!   Purpose: Add a set of new neighbor nodes to a neighbor list and allocate
!            a new neighbor list, if necessary.

!   History: - Creation (13.03.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    INTEGER, INTENT(IN)                   :: atom,nnode
    INTEGER, DIMENSION(3), INTENT(IN)     :: cell
    INTEGER, DIMENSION(:), INTENT(IN)     :: neighbors
    REAL(wp), DIMENSION(:,:), INTENT(IN)  :: r

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_neighbor_nodes (MODULE qs_neighbor_list_types)"
    INTEGER, DIMENSION(3), PARAMETER :: cell0 = (/0,0,0/)

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    INTEGER :: inode,istat

!   ---------------------------------------------------------------------------

!   *** Find neighbor list ***

    neighbor_list => find_neighbor_list(neighbor_list_set,atom,cell)

!   *** Add a new neighbor list to the neighbor list set, ***
!   *** if the requested neighbor list was not found      ***

    IF (.NOT.ASSOCIATED(neighbor_list)) THEN
      CALL add_neighbor_list(neighbor_list_set,atom,cell,neighbor_list)
    END IF

!   *** Add the new neighbor nodes to the neighbor list ***

    DO inode=1,nnode
      CALL add_neighbor_node(neighbor_list,neighbors(inode),r(:,inode),cell0)
    END DO

  END SUBROUTINE add_neighbor_nodes

! *****************************************************************************

  SUBROUTINE allocate_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Allocate and initialize a set of neighbor lists.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    REAL(wp), OPTIONAL                    :: r_max

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_neighbor_list_set (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old neighbor list set ***

    IF (ASSOCIATED(neighbor_list_set)) THEN
      CALL deallocate_neighbor_list_set(neighbor_list_set)
    END IF

!   *** Allocate a set of neighbor lists ***

    ALLOCATE (neighbor_list_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"neighbor_list_set",0)

!   *** Initialize the pointer to the first neighbor list ***

    NULLIFY (neighbor_list_set%first_neighbor_list)

!   *** Initialize the neighbor list counter ***

    neighbor_list_set%nlist = 0

!   *** Initialize the maximum interaction radius (optional) ***

    IF (PRESENT(r_max)) THEN
      neighbor_list_set%r_max = r_max
    ELSE
      neighbor_list_set%r_max = 0.0_wp
    END IF

  END SUBROUTINE allocate_neighbor_list_set

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_list_set(neighbor_list_set)

!   Purpose: Deallocate a neighbor list set.

!   History: - Creation (03.11.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_neighbor_list_set (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list,next_neighbor_list
    TYPE(neighbor_node_type), POINTER :: neighbor_node,next_neighbor_node
    INTEGER                           :: istat

!   ---------------------------------------------------------------------------

    neighbor_list => neighbor_list_set%first_neighbor_list

    DO WHILE (ASSOCIATED(neighbor_list))
      next_neighbor_list => neighbor_list%next_neighbor_list
      neighbor_node => neighbor_list%first_neighbor_node
      DO WHILE (ASSOCIATED(neighbor_node))
        next_neighbor_node => neighbor_node%next_neighbor_node
        DEALLOCATE (neighbor_node,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbor_node")
        neighbor_node => next_neighbor_node
      END DO
      DEALLOCATE (neighbor_list,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"neighbor_list")
      neighbor_list => next_neighbor_list
    END DO

    DEALLOCATE (neighbor_list_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"neighbor_list_set")

  END SUBROUTINE deallocate_neighbor_list_set

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_nodes(neighbor_list,neighbor_node)

!   Purpose: Deallocate all neighbor nodes of a neighbor list starting from
!            the supplied neighbor node.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list
    TYPE(neighbor_node_type), POINTER :: neighbor_node

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_neighbor_nodes (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: next_neighbor_node
    INTEGER                           :: istat

!   ---------------------------------------------------------------------------

    DO WHILE (ASSOCIATED(neighbor_node))
      next_neighbor_node => neighbor_node%next_neighbor_node
      DEALLOCATE (neighbor_node,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"neighbor_node")
      neighbor_list%nnode = neighbor_list%nnode - 1
      neighbor_node => next_neighbor_node
    END DO

  END SUBROUTINE deallocate_neighbor_nodes

! *****************************************************************************

  SUBROUTINE extract_neighbor_list(neighbor_list,neighbors,r)

!   Purpose: Return the entire data set of the reqested neighbor list.

!   History: - Creation (02.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER     :: neighbor_list
    INTEGER, DIMENSION(:), INTENT(OUT)    :: neighbors
    REAL(wp), DIMENSION(:,:), INTENT(OUT) :: r

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE extract_neighbor_list (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: neighbor_node

    INTEGER :: inode,istat,nnode

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      nnode = neighbor_list%nnode

      IF (nnode > 0) THEN

        IF (nnode > SIZE(neighbors)) THEN
          CALL stop_program(routine,"Size of object <neighbors> too small")
        END IF

        IF (3*nnode > SIZE(r)) THEN
          CALL stop_program(routine,"Size of object <r> too small")
        END IF

        inode = 0

        neighbor_node => neighbor_list%first_neighbor_node

        DO WHILE (ASSOCIATED(neighbor_node))
          inode = inode + 1
          neighbors(inode) = neighbor_node%neighbor
          r(:,inode) = neighbor_node%r(:)
          neighbor_node => neighbor_node%next_neighbor_node
        END DO

      END IF

    ELSE

      CALL stop_program(routine,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE extract_neighbor_list

! *****************************************************************************

  FUNCTION find_neighbor_list_1(neighbor_list_set,atom) RESULT(neighbor_list)

!   Purpose: Return a pointer to the neighbor list of atom "atom" in a
!            neighbor list set.

!   History: - Creation (14.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    INTEGER, INTENT(IN)                   :: atom

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   ---------------------------------------------------------------------------

    neighbor_list => neighbor_list_set%first_neighbor_list

    DO WHILE (ASSOCIATED(neighbor_list))
      IF (neighbor_list%atom == atom) EXIT
      neighbor_list => neighbor_list%next_neighbor_list
    END DO

  END FUNCTION find_neighbor_list_1

! *****************************************************************************

  FUNCTION find_neighbor_list_2(neighbor_list_set,atom,cell)&
    RESULT(neighbor_list)

!   Purpose: Return a pointer to the neighbor list of atom "atom" in cell
!            "cell" in a neighbor list set.

!   History: - Creation (14.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    INTEGER, INTENT(IN)                   :: atom
    INTEGER, DIMENSION(3), INTENT(IN)     :: cell

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   ---------------------------------------------------------------------------

    neighbor_list => neighbor_list_set%first_neighbor_list

    DO WHILE (ASSOCIATED(neighbor_list))
      IF (neighbor_list%atom == atom) THEN
        IF ((neighbor_list%cell(1) == cell(1)).AND.&
            (neighbor_list%cell(2) == cell(2)).AND.&
            (neighbor_list%cell(3) == cell(3))) EXIT
      END IF
      neighbor_list => neighbor_list%next_neighbor_list
    END DO

  END FUNCTION find_neighbor_list_2

! *****************************************************************************

  FUNCTION first_list(neighbor_list_set) RESULT(first_neighbor_list)

!   Purpose: Return a pointer to the first neighbor list of a neighbor list
!            set.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

    TYPE(neighbor_list_type), POINTER :: first_neighbor_list

!   ---------------------------------------------------------------------------

    first_neighbor_list => neighbor_list_set%first_neighbor_list

  END FUNCTION first_list

! *****************************************************************************

  FUNCTION first_node(neighbor_list) RESULT(first_neighbor_node)

!   Purpose: Return a pointer to the first neighbor node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    TYPE(neighbor_node_type), POINTER :: first_neighbor_node

!   ---------------------------------------------------------------------------

    first_neighbor_node => neighbor_list%first_neighbor_node

  END FUNCTION first_node

! *****************************************************************************

  SUBROUTINE get_neighbor_list(neighbor_list,atom,cell,nnode)

!   Purpose: Return the reqested data of a neighbor list.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER            :: neighbor_list
    INTEGER, OPTIONAL, INTENT(OUT)               :: atom,nnode
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(OUT) :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE get_neighbor_list (MODULE qs_neighbor_list_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN
      IF (PRESENT(atom)) atom = neighbor_list%atom
      IF (PRESENT(nnode)) nnode = neighbor_list%nnode
      IF (PRESENT(cell)) cell(:) = neighbor_list%cell(:)
    ELSE
      CALL stop_program(routine,&
                        "The requested neighbor list is not associated")
    END IF

  END SUBROUTINE get_neighbor_list

! *****************************************************************************

  SUBROUTINE get_neighbor_list_set(neighbor_list_set,r_max,nlist)

!   Purpose: Return the components of a neighbor list set.

!   History: - Creation (10.11.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    REAL(wp), OPTIONAL, INTENT(OUT)       :: r_max
    INTEGER, OPTIONAL, INTENT(OUT)        :: nlist

!   ---------------------------------------------------------------------------

    IF (PRESENT(r_max)) r_max = neighbor_list_set%r_max
    IF (PRESENT(nlist)) nlist = neighbor_list_set%nlist

  END SUBROUTINE get_neighbor_list_set

! *****************************************************************************

  SUBROUTINE get_neighbor_node(neighbor_node,neighbor,r,cell)

!   Purpose: Return the reqested data of a neighbor node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER             :: neighbor_node
    INTEGER, OPTIONAL, INTENT(OUT)                :: neighbor
    REAL(wp), DIMENSION(3), OPTIONAL, INTENT(OUT) :: r
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(OUT)  :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE get_neighbor_node (MODULE qs_neighbor_list_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_node)) THEN
      IF (PRESENT(neighbor)) neighbor = neighbor_node%neighbor
      IF (PRESENT(r)) r(:) = neighbor_node%r(:)
      IF (PRESENT(cell)) cell(:) = neighbor_node%cell(:)
    ELSE
      CALL stop_program(routine,&
                        "The requested neighbor node is not associated")
    END IF

  END SUBROUTINE get_neighbor_node

! *****************************************************************************

  FUNCTION next_neighbor_list_1(neighbor_list) RESULT(next_list)

!   Purpose: Return the pointer to the next neighbor list.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    TYPE(neighbor_list_type), POINTER :: next_list

!   ---------------------------------------------------------------------------

    next_list => neighbor_list%next_neighbor_list

  END FUNCTION next_neighbor_list_1

! *****************************************************************************

  FUNCTION next_neighbor_list_2(neighbor_list,cell) RESULT(next_list)

!   Purpose: Return the pointer to the next neighbor list of an atom in the
!            cell "cell".

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list
    INTEGER, DIMENSION(3), INTENT(IN) :: cell

    TYPE(neighbor_list_type), POINTER :: next_list

!   ---------------------------------------------------------------------------

    next_list => neighbor_list%next_neighbor_list

    DO WHILE (ASSOCIATED(next_list))
      IF ((next_list%cell(1) == cell(1)).AND.&
          (next_list%cell(2) == cell(2)).AND.&
          (next_list%cell(3) == cell(3))) EXIT
      next_list => next_list%next_neighbor_list
    END DO

  END FUNCTION next_neighbor_list_2

! *****************************************************************************

  FUNCTION next_neighbor_node(neighbor_node) RESULT(next_node)

!   Purpose: Return the pointer to the next neighbor node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER :: neighbor_node

    TYPE(neighbor_node_type), POINTER :: next_node

!   ---------------------------------------------------------------------------

    next_node => neighbor_node%next_neighbor_node

  END FUNCTION next_neighbor_node

! *****************************************************************************

  SUBROUTINE set_neighbor_list(neighbor_list,atom,cell)

!   Purpose: Set the reqested data of a neighbor list.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER           :: neighbor_list
    INTEGER, OPTIONAL, INTENT(IN)               :: atom
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(IN) :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE set_neighbor_list (MODULE qs_neighbor_list_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN
      IF (PRESENT(atom)) neighbor_list%atom = atom
      IF (PRESENT(cell)) neighbor_list%cell(:) = cell(:)
    ELSE
      CALL stop_program(routine,&
                        "The requested neighbor list is not associated")
    END IF

  END SUBROUTINE set_neighbor_list

! *****************************************************************************

  SUBROUTINE set_neighbor_node(neighbor_node,neighbor,r,cell)

!   Purpose: Set the reqested data of a neighbor node.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER            :: neighbor_node
    INTEGER, OPTIONAL, INTENT(IN)                :: neighbor
    REAL(wp), DIMENSION(3), OPTIONAL, INTENT(IN) :: r
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(IN)  :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE set_neighbor_node (MODULE qs_neighbor_list_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_node)) THEN
      IF (PRESENT(neighbor)) neighbor_node%neighbor = neighbor
      IF (PRESENT(r)) neighbor_node%r(:) = r(:)
      IF (PRESENT(cell)) neighbor_node%cell(:) = cell(:)
    ELSE
      CALL stop_program(routine,&
                        "The requested neighbor node is not associated")
    END IF

  END SUBROUTINE set_neighbor_node

! *****************************************************************************

END MODULE qs_neighbor_list_types
