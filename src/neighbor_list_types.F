!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/neighbor_list_types [1.0] *
!!
!!   NAME
!!     neighbor_list_types
!!
!!   FUNCTION
!!     Define the neighbor list data types.
!!
!!   AUTHOR
!!     Matthias Krack (21.06.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE neighbor_list_types

! *****************************************************************************

! Index:

! SUBROUTINE add_neighbor_list(neighbor_list_set,atom,cell)
! SUBROUTINE add_neighbor_node(neighbor_list,r2,neighbor,r,cell)
! SUBROUTINE allocate_neighbor_list_set(neighbor_list_set,&
!                                       neighbor_list_set_name)
! SUBROUTINE deallocate_neighbor_list_set(neighbor_list_set)
! SUBROUTINE extract_neighbor_list(neighbor_list,r2,neighbors,r,cell)
! SUBROUTINE get_neighbor_list(neighbor_list,atom,cell)
! SUBROUTINE get_neighbor_list_set_info(neighbor_list_set,name,nlist)
! SUBROUTINE get_neighbor_node(neighbor_node,r2,neighbor,r,cell)

! FUNCTION find_neighbor_list(neighbor_list_set,atom,cell)&
!   RESULT(neighbor_list)
! FUNCTION first_neighbor_list(neighbor_list_set)
! FUNCTION first_neighbor_node(neighbor_list)
! FUNCTION next_neighbor_list(neighbor_list)
! FUNCTION next_neighbor_node(neighbor_node)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  IMPLICIT NONE

  PRIVATE

! *** Definition of the data types for a linked list of neighbors ***

  TYPE neighbor_node_type
    PRIVATE
    TYPE(neighbor_node_type), POINTER :: next_neighbor_node
    REAL(wp)                          :: r2
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
    CHARACTER(LEN=60)                 :: name
    INTEGER                           :: nlist
  END TYPE neighbor_list_set_type

! *** Public data types ***

  PUBLIC :: neighbor_list_set_type,&
            neighbor_list_type,&
            neighbor_node_type

! *** Public subroutines ***

  PUBLIC :: add_neighbor_list,&
            add_neighbor_node,&
            allocate_neighbor_list_set,&
            deallocate_neighbor_list_set,&
            extract_neighbor_list,&
            get_neighbor_list,&
            get_neighbor_list_set_info,&
            get_neighbor_node

! *** Public functions ***

  PUBLIC :: find_neighbor_list,&
            first_neighbor_list,&
            first_neighbor_node,&
            next_neighbor_list,&
            next_neighbor_node

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE add_neighbor_list(neighbor_list_set,atom,cell,neighbor_list)

!   Purpose: Allocate a new neighbor list.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER     :: neighbor_list
    INTEGER, INTENT(IN)                   :: atom
    INTEGER, DIMENSION(3), INTENT(IN)     :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_neighbor_list (MODULE neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: new_neighbor_list

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Allocate a new neighbor list ***

    ALLOCATE (new_neighbor_list,STAT=istat)

    IF (istat /= 0) CALL stop_memory(routine,neighbor_list_set%name,0)

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

  SUBROUTINE add_neighbor_node(neighbor_list,r2,neighbor,r,cell)

!   Purpose: Allocate a new neighbor list node and store the neighbor data.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(neighbor_list_type), POINTER  :: neighbor_list
    REAL(wp), INTENT(IN)               :: r2
    INTEGER, INTENT(IN)                :: neighbor
    REAL(wp), DIMENSION(3), INTENT(IN) :: r
    INTEGER, DIMENSION(3), INTENT(IN)  :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_neighbor_node (MODULE neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: new_neighbor_node

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Allocate a new neighbor data set ***

    ALLOCATE (new_neighbor_node,STAT=istat)

    IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_node",0)

!   *** Store the data set of the new neighbor ***

    new_neighbor_node%neighbor = neighbor
    new_neighbor_node%r2 = r2
    new_neighbor_node%r(:) = r(:)
    new_neighbor_node%cell(:) = cell(:)

!   *** Link the new neighbor node to the list ***

    new_neighbor_node%next_neighbor_node => neighbor_list%first_neighbor_node

    neighbor_list%first_neighbor_node => new_neighbor_node

!   *** Increment the neighbor node counter ***

    neighbor_list%nnode = neighbor_list%nnode + 1

  END SUBROUTINE add_neighbor_node

! *****************************************************************************

  SUBROUTINE allocate_neighbor_list_set(neighbor_list_set,&
                                        neighbor_list_set_name)

!   Purpose: Allocate and initialize a set of neighbor lists.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    CHARACTER(LEN=*), INTENT(IN)          :: neighbor_list_set_name

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_neighbor_list_set (MODULE neighbor_list_types)"

!   *** Local variables ***

    INTEGER :: istat

!   ---------------------------------------------------------------------------

!   *** Deallocate the old neighbor list set ***

    IF (ASSOCIATED(neighbor_list_set)) THEN
      CALL deallocate_neighbor_list_set(neighbor_list_set)
    END IF

!   *** Allocate a set of neighbor lists ***

    ALLOCATE (neighbor_list_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,neighbor_list_set_name,0)

!   *** Store the name of the neighbor list set ***

    neighbor_list_set%name = neighbor_list_set_name

!   *** Initialize the pointer to the first neighbor list ***

    NULLIFY (neighbor_list_set%first_neighbor_list)

!   *** Initialize the neighbor list counter ***

    neighbor_list_set%nlist = 0

  END SUBROUTINE allocate_neighbor_list_set

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_list_set(neighbor_list_set)

!   Purpose: Deallocate a neighbor list set.

!   History: - Creation (03.11.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_memory

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_neighbor_list_set (MODULE neighbor_list_types)"

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
        IF (istat /= 0) THEN
          CALL stop_memory(routine,neighbor_list_set%name//" (neighbor_node)")
        END IF
        neighbor_node => next_neighbor_node
      END DO
      DEALLOCATE (neighbor_list,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine,neighbor_list_set%name//" (neighbor_list)")
      END IF
      neighbor_list => next_neighbor_list
    END DO

    DEALLOCATE (neighbor_list_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,neighbor_list_set%name)

  END SUBROUTINE deallocate_neighbor_list_set

! *****************************************************************************

  SUBROUTINE extract_neighbor_list(neighbor_list,r2,neighbors,r,cell)

!   Purpose: Return the entire data set of the reqested neighbor list.

!   History: - Creation (02.10.2000, Matthias Krack)

!   ***************************************************************************

    USE memory_utilities, ONLY: reallocate
    USE termination,      ONLY: stop_program

    TYPE(neighbor_list_type), POINTER          :: neighbor_list
    REAL(wp), DIMENSION(:), POINTER            :: r2
    INTEGER, DIMENSION(:), POINTER             :: neighbors
    REAL(wp), DIMENSION(:,:), POINTER          :: r
    INTEGER, DIMENSION(:,:), OPTIONAL, POINTER :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE extract_neighbor_list (MODULE neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: neighbor_node

    INTEGER :: inode

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      IF (neighbor_list%nnode > 0) THEN

        r2 => reallocate(r2,1,neighbor_list%nnode)
        neighbors => reallocate(neighbors,1,neighbor_list%nnode)
        r => reallocate(r,1,3,1,neighbor_list%nnode)
        IF (PRESENT(cell)) cell => reallocate(cell,1,3,1,neighbor_list%nnode)

        inode = 0

        neighbor_node => neighbor_list%first_neighbor_node

        DO WHILE (ASSOCIATED(neighbor_node))
          inode = inode + 1
          r2(inode) = neighbor_node%r2
          neighbors(inode) = neighbor_node%neighbor
          r(:,inode) = neighbor_node%r(:)
          IF (PRESENT(cell)) cell(:,inode) = neighbor_node%cell(:)
          neighbor_node => neighbor_node%next_neighbor_node
        END DO

      END IF

    ELSE

      CALL stop_program(routine,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE extract_neighbor_list

! *****************************************************************************

  FUNCTION find_neighbor_list(neighbor_list_set,atom,cell)&
    RESULT(neighbor_list)

!   Purpose: Return a pointer to the requested neighbor list in a neighbor list
!            set.

!   History: - Creation (14.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER       :: neighbor_list_set
    INTEGER, INTENT(IN)                         :: atom
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(IN) :: cell

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   ---------------------------------------------------------------------------

    neighbor_list => neighbor_list_set%first_neighbor_list

    IF (PRESENT(cell)) THEN
      DO WHILE (ASSOCIATED(neighbor_list))
        IF ((neighbor_list%atom == atom).AND.&
            (neighbor_list%cell(1) == cell(1)).AND.&
            (neighbor_list%cell(2) == cell(2)).AND.&
            (neighbor_list%cell(3) == cell(3))) RETURN
        neighbor_list => neighbor_list%next_neighbor_list
      END DO
    ELSE
      DO WHILE (ASSOCIATED(neighbor_list))
        IF (neighbor_list%atom == atom) RETURN
        neighbor_list => neighbor_list%next_neighbor_list
      END DO
    END IF

  END FUNCTION find_neighbor_list

! *****************************************************************************

  FUNCTION first_neighbor_list(neighbor_list_set)

!   Purpose: Return a pointer to the first neighbor list of a neighbor list
!            set.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

    TYPE(neighbor_list_type), POINTER :: first_neighbor_list

!   ---------------------------------------------------------------------------

    first_neighbor_list => neighbor_list_set%first_neighbor_list

  END FUNCTION first_neighbor_list

! *****************************************************************************

  FUNCTION first_neighbor_node(neighbor_list)

!   Purpose: Return a pointer to the first neighbor node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    TYPE(neighbor_node_type), POINTER :: first_neighbor_node

!   ---------------------------------------------------------------------------

    first_neighbor_node => neighbor_list%first_neighbor_node

  END FUNCTION first_neighbor_node

! *****************************************************************************

  SUBROUTINE get_neighbor_list(neighbor_list,atom,cell,nnode)

!   Purpose: Return the reqested data of a neighbor list.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_program

    TYPE(neighbor_list_type), POINTER            :: neighbor_list
    INTEGER, OPTIONAL, INTENT(OUT)               :: atom,nnode
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(OUT) :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE get_neighbor_list (MODULE neighbor_list_types)"

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

  SUBROUTINE get_neighbor_list_set_info(neighbor_list_set,name,nlist)

!   Purpose: Return the name of a neighbor list set.

!   History: - Creation (10.11.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER   :: neighbor_list_set
    CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT)          :: nlist

!   ---------------------------------------------------------------------------

    IF (PRESENT(name)) name = TRIM(neighbor_list_set%name)
    IF (PRESENT(nlist)) nlist = neighbor_list_set%nlist

  END SUBROUTINE get_neighbor_list_set_info

! *****************************************************************************

  SUBROUTINE get_neighbor_node(neighbor_node,r2,neighbor,r,cell)

!   Purpose: Return the reqested data of a neighbor node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    USE termination, ONLY: stop_program

    TYPE(neighbor_node_type), POINTER             :: neighbor_node
    REAL(wp), OPTIONAL, INTENT(OUT)               :: r2
    INTEGER, OPTIONAL, INTENT(OUT)                :: neighbor
    REAL(wp), DIMENSION(3), OPTIONAL, INTENT(OUT) :: r
    INTEGER, DIMENSION(3), OPTIONAL, INTENT(OUT)  :: cell

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE get_neighbor_node (MODULE neighbor_list_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_node)) THEN
      IF (PRESENT(r2)) r2 = neighbor_node%r2
      IF (PRESENT(neighbor)) neighbor = neighbor_node%neighbor
      IF (PRESENT(r)) r(:) = neighbor_node%r(:)
      IF (PRESENT(cell)) cell(:) = neighbor_node%cell(:)
    ELSE
      CALL stop_program(routine,&
                        "The requested neighbor node is not associated")
    END IF

  END SUBROUTINE get_neighbor_node

! *****************************************************************************

  FUNCTION next_neighbor_list(neighbor_list)

!   Purpose: Return the pointer to the subsequent neighbor list.

!   History: - Creation (13.09.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    TYPE(neighbor_list_type), POINTER :: next_neighbor_list

!   ---------------------------------------------------------------------------

    next_neighbor_list => neighbor_list%next_neighbor_list

  END FUNCTION next_neighbor_list

! *****************************************************************************

  FUNCTION next_neighbor_node(neighbor_node)

!   Purpose: Return the pointer to the next neighbor node.

!   History: - Creation (23.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER :: neighbor_node

    TYPE(neighbor_node_type), POINTER :: next_neighbor_node

!   ---------------------------------------------------------------------------

    next_neighbor_node => neighbor_node%next_neighbor_node

  END FUNCTION next_neighbor_node

! *****************************************************************************

END MODULE neighbor_list_types
