!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_neighbor_list_types [1.0] *
!!
!!   NAME
!!     qs_neighbor_list_types
!!
!!   FUNCTION
!!     Define the neighbor list data types and the corresponding functionality
!!
!!   AUTHOR
!!     Matthias Krack (21.06.2000)
!!
!!   MODIFICATION HISTORY
!!     - cleaned (23.07.2003,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_neighbor_list_types

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_neighbor_list_types"

! *** Definition of the data types for a linked list of neighbors ***

  TYPE block_p_type
    REAL(KIND = dp), DIMENSION(:,:), POINTER :: block
  END TYPE block_p_type

  TYPE qlist_type
    TYPE(block_p_type), DIMENSION(:), POINTER :: sac
    REAL(KIND = dp), DIMENSION(:), POINTER    :: r2
    REAL(KIND = dp), DIMENSION(:,:), POINTER  :: r
    INTEGER                                   :: n
    INTEGER, DIMENSION(:), POINTER            :: index_list,list
  END TYPE qlist_type

  TYPE reduced_3c_list_type
    INTEGER                           :: nnode
    INTEGER, DIMENSION(:),  POINTER   :: index_atom
    REAL(dp), DIMENSION(:), POINTER   :: rac2, rbc2
    REAL(dp), DIMENSION(:,:), POINTER :: rac, rbc
  END TYPE reduced_3c_list_type

  TYPE neighbor_node_type
    PRIVATE
    TYPE(neighbor_node_type), POINTER :: next_neighbor_node
    INTEGER                           :: neighbor
    REAL(dp), DIMENSION(3)            :: r,s
    INTEGER, DIMENSION(3)             :: cell
    TYPE(reduced_3c_list_type),&
        DIMENSION(:),POINTER          :: reduced_3c_oce,&
                                         reduced_3c_rho0
  END TYPE neighbor_node_type

  TYPE neighbor_list_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: next_neighbor_list
    TYPE(neighbor_node_type), POINTER :: first_neighbor_node,&
                                         last_neighbor_node
    INTEGER                           :: atom,nnode
    INTEGER, DIMENSION(3)             :: cell
  END TYPE neighbor_list_type

  TYPE neighbor_list_set_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: first_neighbor_list,&
                                         last_neighbor_list
    REAL(dp)                          :: r_max
    INTEGER                           :: nlist
  END TYPE neighbor_list_set_type

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
            neighbor_node_type,&
            qlist_type,&
            reduced_3c_list_type,&
            block_p_type

! *** Public subroutines ***

  PUBLIC :: add_neighbor_list,&
            add_neighbor_node,&
            allocate_neighbor_list_set,&
            clean_neighbor_list_set,&
            deallocate_neighbor_list,&
            deallocate_neighbor_list_set,&
            deallocate_reduced_3c_list,&
            extract_neighbor_list,&
            get_neighbor_list,&
            get_neighbor_list_set,&
            get_neighbor_node,&
            init_neighbor_list,&
            init_neighbor_list_set,&
            set_neighbor_list,&
            set_neighbor_list_set,&
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
    MODULE PROCEDURE next_neighbor_list,&
                     next_neighbor_node
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE add_neighbor_list(neighbor_list_set,atom,cell,neighbor_list)

!   Purpose: Add a new neighbor list to a neighbor list set.

!   History: - Creation (13.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    INTEGER, INTENT(IN)                      :: atom
    INTEGER, DIMENSION(3), INTENT(IN)        :: cell
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(LEN=*), PARAMETER :: routine_name = "add_neighbor_list"

    INTEGER                                  :: istat
    TYPE(neighbor_list_type), POINTER        :: new_neighbor_list

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      IF (ASSOCIATED(neighbor_list_set%last_neighbor_list)) THEN

        new_neighbor_list =>&
          neighbor_list_set%last_neighbor_list%next_neighbor_list

        IF (.NOT.ASSOCIATED(new_neighbor_list)) THEN

!         *** Allocate a new neighbor list ***

          ALLOCATE (new_neighbor_list,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "new_neighbor_list",0)
          END IF

          NULLIFY (new_neighbor_list%next_neighbor_list)
          NULLIFY (new_neighbor_list%first_neighbor_node)

!         *** Link the new neighbor list to the neighbor list set ***

          neighbor_list_set%last_neighbor_list%next_neighbor_list => new_neighbor_list

        END IF

      ELSE

        new_neighbor_list => neighbor_list_set%first_neighbor_list

        IF (.NOT.ASSOCIATED(new_neighbor_list)) THEN

!         *** Allocate a new first neighbor list ***

          ALLOCATE (new_neighbor_list,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "new_neighbor_list",0)
          END IF

          NULLIFY (new_neighbor_list%next_neighbor_list)
          NULLIFY (new_neighbor_list%first_neighbor_node)

!         *** Link the new first neighbor list to the neighbor list set ***

          neighbor_list_set%first_neighbor_list => new_neighbor_list

        END IF

      END IF

!     *** Store the data set of the new neighbor list ***

      NULLIFY (new_neighbor_list%last_neighbor_node)
      new_neighbor_list%atom = atom
      new_neighbor_list%nnode = 0
      new_neighbor_list%cell(:) = cell(:)

!     *** Update the pointer to the last neighbor ***
!     *** list of the neighbor list set           ***

      neighbor_list_set%last_neighbor_list => new_neighbor_list

!     *** Increment the neighbor list counter ***

      neighbor_list_set%nlist = neighbor_list_set%nlist + 1

!     *** Return a pointer to the new neighbor list ***

      neighbor_list => new_neighbor_list

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE add_neighbor_list

! *****************************************************************************

  SUBROUTINE add_neighbor_node(neighbor_list,neighbor,cell,r,s,exclusion_list,&
                               l3c_oce,l3c_rho0,nkind)

!   Purpose: Add a new neighbor list node to a neighbor list.

!   History: - Creation (23.06.2000,MK)
!            - Exclusion list added for classic code (MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    INTEGER, INTENT(IN)                      :: neighbor
    INTEGER, DIMENSION(3), INTENT(IN)        :: cell
    REAL(dp), DIMENSION(3), INTENT(IN)       :: r
    REAL(dp), DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: s
    INTEGER, DIMENSION(:), POINTER, &
      OPTIONAL                               :: exclusion_list
    LOGICAL, INTENT(IN), OPTIONAL            :: l3c_oce, l3c_rho0
    INTEGER, INTENT(IN), OPTIONAL            :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "add_neighbor_node"

    INTEGER                                  :: iatom, istat,my_nkind
    LOGICAL                                  :: my_l3c_oce,my_l3c_rho0
    TYPE(neighbor_node_type), POINTER        :: new_neighbor_node

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

!     *** Check for exclusions ***

      IF (PRESENT(exclusion_list)) THEN
        IF ( ASSOCIATED ( exclusion_list ) ) THEN
          DO iatom=1,SIZE(exclusion_list)
            IF (exclusion_list(iatom) == 0) EXIT
            IF (exclusion_list(iatom) == neighbor) RETURN
          END DO
        END IF
      END IF

      my_l3c_oce = .FALSE.
      IF(PRESENT(l3c_oce)) my_l3c_oce = l3c_oce
      my_l3c_rho0 = .FALSE.
      IF(PRESENT(l3c_rho0)) my_l3c_rho0 = l3c_rho0
      my_nkind = 0
      IF(PRESENT(nkind)) my_nkind = nkind

      IF (ASSOCIATED(neighbor_list%last_neighbor_node)) THEN

        new_neighbor_node => neighbor_list%last_neighbor_node%next_neighbor_node

        IF (.NOT.ASSOCIATED(new_neighbor_node)) THEN

!         *** Allocate a new neighbor node ***

          ALLOCATE (new_neighbor_node,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "new_neighbor_node",0)
          END IF

          NULLIFY (new_neighbor_node%next_neighbor_node)

!         *** Link the new neighbor node to the neighbor list ***

          neighbor_list%last_neighbor_node%next_neighbor_node => new_neighbor_node

        END IF

      ELSE

        new_neighbor_node => neighbor_list%first_neighbor_node

        IF (.NOT.ASSOCIATED(new_neighbor_node)) THEN

!         *** Allocate a new first neighbor node ***

          ALLOCATE (new_neighbor_node,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "new_neighbor_node",0)
          END IF

          NULLIFY (new_neighbor_node%next_neighbor_node)

!         *** Link the new first neighbor node to the neighbor list ***

          neighbor_list%first_neighbor_node => new_neighbor_node

        END IF

      END IF

!     *** Store the data set of the new neighbor ***

      new_neighbor_node%neighbor = neighbor
      new_neighbor_node%cell(:) = cell(:)
      new_neighbor_node%r(:) = r(:)
      IF (PRESENT(s)) new_neighbor_node%s(:) = s(:)
      IF (my_l3c_oce) THEN
        CALL create_reduced_3c_list(new_neighbor_node%reduced_3c_oce,my_nkind)
      ELSE
          NULLIFY(new_neighbor_node%reduced_3c_oce)
      END IF
      IF (my_l3c_rho0) THEN
        CALL create_reduced_3c_list(new_neighbor_node%reduced_3c_rho0,my_nkind)
      ELSE
          NULLIFY(new_neighbor_node%reduced_3c_rho0)
      END IF

!     *** Update the pointer to the last neighbor node of the neighbor list ***

      neighbor_list%last_neighbor_node => new_neighbor_node

!     *** Increment the neighbor node counter ***

      neighbor_list%nnode = neighbor_list%nnode + 1

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE add_neighbor_node

! *****************************************************************************

  SUBROUTINE allocate_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Allocate and initialize a set of neighbor lists.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp)                                 :: r_max

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_neighbor_list_set"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------
!   *** Deallocate the old neighbor list set ***

    IF (ASSOCIATED(neighbor_list_set)) THEN
      CALL deallocate_neighbor_list_set(neighbor_list_set)
    END IF

!   *** Allocate a set of neighbor lists ***

    ALLOCATE (neighbor_list_set,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "neighbor_list_set",0)
    END IF

    NULLIFY (neighbor_list_set%first_neighbor_list)

!   *** Initialize the pointers to the first neighbor list ***

    CALL init_neighbor_list_set(neighbor_list_set,r_max)

  END SUBROUTINE allocate_neighbor_list_set

! *****************************************************************************

  SUBROUTINE clean_neighbor_list_set(neighbor_list_set)

!   Purpose: Deallocate all unused neighbor lists and neighbor nodes in a
!            neighbor list set.

!   History: - Creation (20.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "clean_neighbor_list_set"

    INTEGER                                  :: istat
    TYPE(neighbor_list_type), POINTER        :: neighbor_list, &
                                                next_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node, &
                                                next_neighbor_node

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

!     *** Deallocate all unused neighbor lists ***

      IF (ASSOCIATED(neighbor_list_set%last_neighbor_list)) THEN
        neighbor_list => neighbor_list_set%last_neighbor_list%next_neighbor_list
        NULLIFY (neighbor_list_set%last_neighbor_list%next_neighbor_list)
      ELSE
        neighbor_list => neighbor_list_set%first_neighbor_list
        NULLIFY (neighbor_list_set%first_neighbor_list)
      END IF

      DO WHILE (ASSOCIATED(neighbor_list))
        next_neighbor_list => neighbor_list%next_neighbor_list
        CALL  deallocate_neighbor_list(neighbor_list)
        neighbor_list => next_neighbor_list
      END DO

!     *** Deallocate all unused neighbor nodes in the used neighbor lists ***

      neighbor_list => neighbor_list_set%first_neighbor_list

      DO WHILE (ASSOCIATED(neighbor_list))

        next_neighbor_list => neighbor_list%next_neighbor_list

        IF (ASSOCIATED(neighbor_list%last_neighbor_node)) THEN
          neighbor_node => neighbor_list%last_neighbor_node%next_neighbor_node
          NULLIFY (neighbor_list%last_neighbor_node%next_neighbor_node)
        ELSE
          neighbor_node => neighbor_list%first_neighbor_node
          NULLIFY (neighbor_list%first_neighbor_node)
        END IF

        DO WHILE (ASSOCIATED(neighbor_node))
          next_neighbor_node => neighbor_node%next_neighbor_node
          DEALLOCATE (neighbor_node,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,"neighbor_node")
          END IF
          neighbor_node => next_neighbor_node
        END DO

        neighbor_list => next_neighbor_list

      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE clean_neighbor_list_set

! *****************************************************************************

  SUBROUTINE create_reduced_3c_list(reduced_3c_list,nkind)

    TYPE(reduced_3c_list_type), DIMENSION(:),&
                               POINTER       :: reduced_3c_list
    INTEGER, INTENT(IN)                      :: nkind

    CHARACTER(LEN=*), PARAMETER :: routine_name = "create_reduced_3c_list"

    INTEGER                                  :: istat, ikind
 
!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(reduced_3c_list)) &
           CALL deallocate_reduced_3c_list(reduced_3c_list)
    ALLOCATE(reduced_3c_list(nkind),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routine_name,module_name,__LINE__,"reduced_3c_list")
    END IF
    DO ikind= 1,nkind
      NULLIFY(reduced_3c_list(ikind)%index_atom)
      NULLIFY(reduced_3c_list(ikind)%rac)
      NULLIFY(reduced_3c_list(ikind)%rac2)
      NULLIFY(reduced_3c_list(ikind)%rbc)
      NULLIFY(reduced_3c_list(ikind)%rbc2)
    END DO

  END SUBROUTINE create_reduced_3c_list

! *****************************************************************************

  SUBROUTINE deallocate_reduced_3c_list(reduced_3c_list)

    TYPE(reduced_3c_list_type), DIMENSION(:),&
                               POINTER       :: reduced_3c_list

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_reduced_3c_list"

    INTEGER                                  :: istat, ikind, nkind

!   ---------------------------------------------------------------------------

    IF(ASSOCIATED(reduced_3c_list)) THEN

      nkind = SIZE(reduced_3c_list)

      DO ikind = 1, nkind
        IF(ASSOCIATED(reduced_3c_list(ikind)%index_atom)) THEN
          DEALLOCATE(reduced_3c_list(ikind)%index_atom,STAT=istat)
          DEALLOCATE(reduced_3c_list(ikind)%rac,STAT=istat)
          DEALLOCATE(reduced_3c_list(ikind)%rac2,STAT=istat)
          DEALLOCATE(reduced_3c_list(ikind)%rbc,STAT=istat)
          DEALLOCATE(reduced_3c_list(ikind)%rbc2,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
               "index_atom,rac,rac2,rbc,rbc2")
          END IF
        END IF
      END DO

      DEALLOCATE (reduced_3c_list,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"reduced_3c_list")
      END IF

    END IF 

  END SUBROUTINE deallocate_reduced_3c_list

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_list(neighbor_list)

!   Purpose: Deallocate a neighbor list.

!   History: - Creation (20.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_neighbor_list"

    INTEGER                                  :: istat
    TYPE(neighbor_node_type), POINTER        :: neighbor_node, &
                                                next_neighbor_node

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      neighbor_node => neighbor_list%first_neighbor_node

      DO WHILE (ASSOCIATED(neighbor_node))
        next_neighbor_node => neighbor_node%next_neighbor_node
        IF(ASSOCIATED(neighbor_node%reduced_3c_oce)) &
               CALL deallocate_reduced_3c_list(neighbor_node%reduced_3c_oce)
        IF(ASSOCIATED(neighbor_node%reduced_3c_rho0)) &
               CALL deallocate_reduced_3c_list(neighbor_node%reduced_3c_rho0)
        DEALLOCATE (neighbor_node,STAT=istat)
        IF (istat /= 0) THEN
          CALL stop_memory(routine_name,module_name,__LINE__,"neighbor_node")
        END IF
        neighbor_node => next_neighbor_node
      END DO

      DEALLOCATE (neighbor_list,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"neighbor_list")
      END IF

    END IF

  END SUBROUTINE deallocate_neighbor_list

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_list_set(neighbor_list_set)

!   Purpose: Deallocate a neighbor list set.

!   History: - Creation (03.11.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set

    CHARACTER(LEN=*), PARAMETER :: &
      routine_name = "deallocate_neighbor_list_set"

    INTEGER                                  :: istat
    TYPE(neighbor_list_type), POINTER        :: neighbor_list, &
                                                next_neighbor_list

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      neighbor_list => neighbor_list_set%first_neighbor_list

      DO WHILE (ASSOCIATED(neighbor_list))
        next_neighbor_list => neighbor_list%next_neighbor_list
        CALL  deallocate_neighbor_list(neighbor_list)
        neighbor_list => next_neighbor_list
      END DO

      DEALLOCATE (neighbor_list_set,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"neighbor_list_set")
      END IF

    END IF

  END SUBROUTINE deallocate_neighbor_list_set

! *****************************************************************************

  SUBROUTINE extract_neighbor_list(neighbor_list,neighbors,r)

!   Purpose: Return the entire data set of the reqested neighbor list.

!   History: - Creation (02.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    INTEGER, DIMENSION(:), INTENT(OUT)       :: neighbors
    REAL(dp), DIMENSION(:, :), INTENT(OUT)   :: r

    CHARACTER(LEN=*), PARAMETER :: routine_name = "extract_neighbor_list"

    INTEGER                                  :: inode, nnode
    TYPE(neighbor_node_type), POINTER        :: neighbor_node

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      nnode = neighbor_list%nnode

      IF (nnode > 0) THEN

        IF (nnode > SIZE(neighbors)) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "The size of the data object <neighbors> is too small")
        END IF

        IF (3*nnode > SIZE(r)) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "The size of the data object <r> is too small")
        END IF

        neighbor_node => neighbor_list%first_neighbor_node

        DO inode=1,nnode
          neighbors(inode) = neighbor_node%neighbor
          r(:,inode) = neighbor_node%r(:)
          neighbor_node => neighbor_node%next_neighbor_node
        END DO

      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE extract_neighbor_list

! *****************************************************************************

  FUNCTION find_neighbor_list_1(neighbor_list_set,atom)&
    RESULT(neighbor_list)

!   Purpose: Return a pointer to the neighbor list of atom "atom" in a
!            neighbor list set. The whole list is traversed.

!   History: - Creation (14.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    INTEGER, INTENT(IN)                      :: atom
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(LEN=*), PARAMETER :: routine_name = "find_neighbor_list_1"

    INTEGER                                  :: ilist
    TYPE(neighbor_list_type), POINTER        :: current_neighbor_list

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      NULLIFY (neighbor_list)

      current_neighbor_list => neighbor_list_set%first_neighbor_list

      DO ilist=1,neighbor_list_set%nlist
        IF (current_neighbor_list%atom == atom) THEN
          neighbor_list => current_neighbor_list
          EXIT
        END IF
        current_neighbor_list => current_neighbor_list%next_neighbor_list
      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END FUNCTION find_neighbor_list_1

! *****************************************************************************

  FUNCTION find_neighbor_list_2(neighbor_list_set,atom,cell)&
    RESULT(neighbor_list)

!   Purpose: Return a pointer to the neighbor list of atom "atom" in cell
!            "cell" in a neighbor list set. The whole list is traversed.

!   History: - Creation (14.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    INTEGER, INTENT(IN)                      :: atom
    INTEGER, DIMENSION(3), INTENT(IN)        :: cell
    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(LEN=*), PARAMETER :: routine_name = "find_neighbor_list_2"

    INTEGER                                  :: ilist
    TYPE(neighbor_list_type), POINTER        :: current_neighbor_list

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      NULLIFY (neighbor_list)

      current_neighbor_list => neighbor_list_set%first_neighbor_list

      DO ilist=1,neighbor_list_set%nlist
        IF (current_neighbor_list%atom == atom) THEN
          IF ((current_neighbor_list%cell(1) == cell(1)).AND.&
              (current_neighbor_list%cell(2) == cell(2)).AND.&
              (current_neighbor_list%cell(3) == cell(3))) THEN
            neighbor_list => current_neighbor_list
            EXIT
          END IF
        END IF
        current_neighbor_list => current_neighbor_list%next_neighbor_list
      END DO

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END FUNCTION find_neighbor_list_2

! *****************************************************************************

  FUNCTION first_list(neighbor_list_set)&
    RESULT(first_neighbor_list)

!   Purpose: Return a pointer to the first neighbor list of a neighbor list
!            set.

!   History: - Creation (13.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: first_neighbor_list

!   ---------------------------------------------------------------------------

    first_neighbor_list => neighbor_list_set%first_neighbor_list

  END FUNCTION first_list

! *****************************************************************************

  FUNCTION first_node(neighbor_list)&
    RESULT(first_neighbor_node)

!   Purpose: Return a pointer to the first neighbor node of a neighbor list.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: first_neighbor_node

!   ---------------------------------------------------------------------------

    first_neighbor_node => neighbor_list%first_neighbor_node

  END FUNCTION first_node

! *****************************************************************************

  SUBROUTINE get_neighbor_list(neighbor_list,atom,cell,nnode)

!   Purpose: Return the reqested data of a neighbor list.

!   History: - Creation (13.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    INTEGER, INTENT(OUT), OPTIONAL           :: atom
    INTEGER, DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: cell
    INTEGER, INTENT(OUT), OPTIONAL           :: nnode

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_neighbor_list"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      IF (PRESENT(atom)) atom = neighbor_list%atom
      IF (PRESENT(nnode)) nnode = neighbor_list%nnode
      IF (PRESENT(cell)) cell(:) = neighbor_list%cell(:)

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE get_neighbor_list

! *****************************************************************************

  SUBROUTINE get_neighbor_list_set(neighbor_list_set,r_max,nlist)

!   Purpose: Return the components of a neighbor list set.

!   History: - Creation (10.11.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp), INTENT(OUT), OPTIONAL          :: r_max
    INTEGER, INTENT(OUT), OPTIONAL           :: nlist

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_neighbor_list_set"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      IF (PRESENT(r_max)) r_max = neighbor_list_set%r_max
      IF (PRESENT(nlist)) nlist = neighbor_list_set%nlist

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE get_neighbor_list_set

! *****************************************************************************

  SUBROUTINE get_neighbor_node(neighbor_node,neighbor,cell,r,s,&
                               reduced_3c_oce,reduced_3c_rho0)

!   Purpose: Return the reqested data of a neighbor node.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    INTEGER, INTENT(OUT), OPTIONAL           :: neighbor
    INTEGER, DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: cell
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: r, s
    TYPE(reduced_3c_list_type), &
        DIMENSION(:), POINTER, OPTIONAL      :: reduced_3c_oce, &
                                                reduced_3c_rho0

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_neighbor_node"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_node)) THEN

      IF (PRESENT(neighbor)) neighbor = neighbor_node%neighbor
      IF (PRESENT(r)) r(:) = neighbor_node%r(:)
      IF (PRESENT(s)) s(:) = neighbor_node%s(:)
      IF (PRESENT(cell)) cell(:) = neighbor_node%cell(:)
      IF (PRESENT(reduced_3c_oce)) &
                  reduced_3c_oce => neighbor_node%reduced_3c_oce
      IF (PRESENT(reduced_3c_rho0)) &
                  reduced_3c_rho0 => neighbor_node%reduced_3c_rho0

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor node is not associated")

    END IF

  END SUBROUTINE get_neighbor_node

! *****************************************************************************

  SUBROUTINE init_neighbor_list(neighbor_list)

!   Purpose: Initialize a neighbor list. Nothing is (de)allocated here.
!            This routine is also used to prepare a neighbor list for
!            overwriting.

!   History: - Creation (21.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_neighbor_list"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

!     *** Initialize the pointers to the last neighbor node ***

      NULLIFY (neighbor_list%last_neighbor_node)

!     *** Initialize the neighbor list counter ***

      neighbor_list%nnode = 0

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE init_neighbor_list

! *****************************************************************************

  SUBROUTINE init_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Initialize a neighbor list set. Nothing is (de)allocated here.
!            This routine is also used to prepare a neighbor list set for
!            overwriting.

!   History: - Creation (20.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp), INTENT(IN), OPTIONAL           :: r_max

    CHARACTER(LEN=*), PARAMETER :: routine_name = "init_neighbor_list_set"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

!     *** Initialize the pointers to the last neighbor list ***

      NULLIFY (neighbor_list_set%last_neighbor_list)

!     *** Initialize the neighbor list counter ***

      neighbor_list_set%nlist = 0

!     *** Initialize the maximum interaction radius (optional) ***

      IF (PRESENT(r_max)) neighbor_list_set%r_max = r_max

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE init_neighbor_list_set

! *****************************************************************************

  FUNCTION next_neighbor_list(neighbor_list)&
    RESULT (next_list)

!   Purpose: Return the pointer to the next neighbor list.

!   History: - Creation (13.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list, next_list

!   ---------------------------------------------------------------------------

    next_list => neighbor_list%next_neighbor_list

  END FUNCTION next_neighbor_list

! *****************************************************************************

  FUNCTION next_neighbor_node(neighbor_node) RESULT(next_node)

!   Purpose: Return the pointer to the next neighbor node.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER        :: neighbor_node, next_node

!   ---------------------------------------------------------------------------

    next_node => neighbor_node%next_neighbor_node

  END FUNCTION next_neighbor_node

! *****************************************************************************

  SUBROUTINE set_neighbor_list(neighbor_list,atom,cell)

!   Purpose: Set the reqested data of a neighbor list.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    INTEGER, INTENT(IN), OPTIONAL            :: atom
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: cell

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_neighbor_list"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      IF (PRESENT(atom)) neighbor_list%atom = atom
      IF (PRESENT(cell)) neighbor_list%cell(:) = cell(:)

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list is not associated")

    END IF

  END SUBROUTINE set_neighbor_list

! *****************************************************************************

  SUBROUTINE set_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Set the components of a neighbor list set.

!   History: - Creation (18.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    REAL(dp), INTENT(IN), OPTIONAL           :: r_max

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_neighbor_list_set"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      IF (PRESENT(r_max)) neighbor_list_set%r_max = r_max

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor list set is not associated")

    END IF

  END SUBROUTINE set_neighbor_list_set

! *****************************************************************************

  SUBROUTINE set_neighbor_node(neighbor_node,neighbor,cell,r,s)

!   Purpose: Set the reqested data of a neighbor node.

!   History: - Creation (10.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    INTEGER, INTENT(IN), OPTIONAL            :: neighbor
    INTEGER, DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: cell
    REAL(dp), DIMENSION(3), INTENT(IN), &
      OPTIONAL                               :: r, s

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_neighbor_node"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_node)) THEN

      IF (PRESENT(neighbor)) neighbor_node%neighbor = neighbor
      IF (PRESENT(r)) neighbor_node%r(:) = r(:)
      IF (PRESENT(s)) neighbor_node%s(:) = s(:)
      IF (PRESENT(cell)) neighbor_node%cell(:) = cell(:)

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The requested neighbor node is not associated")

    END IF

  END SUBROUTINE set_neighbor_node

! *****************************************************************************

END MODULE qs_neighbor_list_types
