!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_neighbor_list_types [1.0] *
!!
!!   NAME
!!     qs_neighbor_list_types
!!
!!   FUNCTION
!!     Define the neighbor list data types and the corresponding functionality
!!     for Quickstep.
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_check,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='qs_neighbor_list_types'
  INTEGER, SAVE, PRIVATE :: last_neigh_group_id=0

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
    TYPE(neighbor_node_type), POINTER :: first_neighbor_node,&
                                         last_neighbor_node
    INTEGER                           :: atom,nnode
    INTEGER, DIMENSION(3)             :: cell
  END TYPE neighbor_list_type

  TYPE neighbor_list_set_type
    PRIVATE
    TYPE(neighbor_list_type), POINTER :: first_neighbor_list,&
                                         last_neighbor_list
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
            neighbor_node_type,&
            neighbor_group_type,&
            neighbor_group_p_type

! *** Public subroutines ***

  PUBLIC :: add_neighbor_list,&
            add_neighbor_node,&
            allocate_neighbor_list_set,&
            clean_neighbor_list_set,&
            deallocate_neighbor_list,&
            deallocate_neighbor_list_set,&
            extract_neighbor_list,&
            get_neighbor_list,&
            get_neighbor_list_set,&
            get_neighbor_node,&
            init_neighbor_list,&
            init_neighbor_list_set,&
            set_neighbor_list,&
            set_neighbor_list_set,&
            set_neighbor_node,&
            neigh_group_create,&
            neigh_group_retain,&
            neigh_group_release

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
!!***

!!****s* qs_neighbor_list_types/neighbor_group_type [1.0] *
!!
!!   NAME
!!     neighbor_group_type
!!
!!   FUNCTION
!!     represents the upper triangualr part of a group of kind related
!!     neighbor lists: the neighbors of atoms of kind 1 that are of
!!     kind 1, the neighbors of atoms of kind 1 that are of
!!     kind 2,...,the neighbors of atoms of kind 2 that are of
!!     kind 2,...
!!     
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - n_kinds: number of kinds
!!     - neigh_sets: the neighbor lists
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE neighbor_group_type
     INTEGER :: id_nr, ref_count, n_kinds
     LOGICAL :: should_dealloc_neigh_sets
     TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: neigh_sets
  END TYPE neighbor_group_type
!!***
  !****************************************************************************

!!****s* qs_neighbor_list_p_types/neighbor_group_p_type [1.0] *
!!
!!   NAME
!!     neighbor_group_p_type
!!
!!   FUNCTION
!!     to build arrays of neighbor_group_type pointers
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - neigh_group: the neighbor_group_type
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE neighbor_group_p_type
     TYPE(neighbor_group_type), POINTER :: neigh_group
  END TYPE neighbor_group_p_type
!!***
  !****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE add_neighbor_list(neighbor_list_set,atom,cell,neighbor_list)

!   Purpose: Add a new neighbor list to neighbor list set.

!   History: - Creation (13.09.2000,MK)

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
    INTEGER                           :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set%last_neighbor_list)) THEN

      new_neighbor_list =>&
        neighbor_list_set%last_neighbor_list%next_neighbor_list

      IF (.NOT.ASSOCIATED(new_neighbor_list)) THEN

!       *** Allocate a new neighbor list ***

        ALLOCATE (new_neighbor_list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_list (1)",0)

        NULLIFY (new_neighbor_list%next_neighbor_list)
        NULLIFY (new_neighbor_list%first_neighbor_node)

!       *** Link the new neighbor list to the neighbor list set ***

        neighbor_list_set%last_neighbor_list%next_neighbor_list =>&
          new_neighbor_list

      END IF

    ELSE

      new_neighbor_list => neighbor_list_set%first_neighbor_list

      IF (.NOT.ASSOCIATED(new_neighbor_list)) THEN

!       *** Allocate a new first neighbor list ***

        ALLOCATE (new_neighbor_list,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_list (2)",0)

        NULLIFY (new_neighbor_list%next_neighbor_list)
        NULLIFY (new_neighbor_list%first_neighbor_node)

!       *** Link the new first neighbor list to the neighbor list set ***

        neighbor_list_set%first_neighbor_list => new_neighbor_list

      END IF

    END IF

!   *** Store the data set of the new neighbor list ***

    NULLIFY (new_neighbor_list%last_neighbor_node)
    new_neighbor_list%atom = atom
    new_neighbor_list%nnode = 0
    new_neighbor_list%cell(:) = cell(:)

!   *** Update the pointer to the last neighbor ***
!   *** list of the neighbor list set           ***

    neighbor_list_set%last_neighbor_list => new_neighbor_list

!   *** Increment the neighbor list counter ***

    neighbor_list_set%nlist = neighbor_list_set%nlist + 1

!   *** Return a pointer to the new neighbor list ***

    neighbor_list => new_neighbor_list

  END SUBROUTINE add_neighbor_list

! *****************************************************************************

  SUBROUTINE add_neighbor_node(neighbor_list,neighbor,r,cell,exclusion_list)

!   Purpose: Add a new neighbor list node to neighbor list.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER           :: neighbor_list
    INTEGER, INTENT(IN)                         :: neighbor
    REAL(wp), DIMENSION(3), INTENT(IN)          :: r
    INTEGER, DIMENSION(3), INTENT(IN)           :: cell
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: exclusion_list

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE add_neighbor_node (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: new_neighbor_node
    INTEGER                           :: iatom,istat

!   ---------------------------------------------------------------------------

!   *** Check for exclusions ***

    IF (PRESENT(exclusion_list)) THEN
      DO iatom=1,SIZE(exclusion_list)
        IF (exclusion_list(iatom) == 0) EXIT
        IF (exclusion_list(iatom) == neighbor) RETURN
      END DO
    END IF

    IF (ASSOCIATED(neighbor_list%last_neighbor_node)) THEN

      new_neighbor_node => neighbor_list%last_neighbor_node%next_neighbor_node

      IF (.NOT.ASSOCIATED(new_neighbor_node)) THEN

!       *** Allocate a new neighbor node ***

        ALLOCATE (new_neighbor_node,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_node (1)",0)
        NULLIFY (new_neighbor_node%next_neighbor_node)

!       *** Link the new neighbor node to the neighbor list ***

        neighbor_list%last_neighbor_node%next_neighbor_node =>&
          new_neighbor_node

      END IF

    ELSE

      new_neighbor_node => neighbor_list%first_neighbor_node

      IF (.NOT.ASSOCIATED(new_neighbor_node)) THEN

!       *** Allocate a new first neighbor node ***

        ALLOCATE (new_neighbor_node,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"new_neighbor_node (2)",0)
        NULLIFY (new_neighbor_node%next_neighbor_node)

!       *** Link the new first neighbor node to the neighbor list ***

        neighbor_list%first_neighbor_node => new_neighbor_node

      END IF

    END IF

!   *** Store the data set of the new neighbor ***

    new_neighbor_node%neighbor = neighbor
    new_neighbor_node%r(:) = r(:)
    new_neighbor_node%cell(:) = cell(:)

!   *** Update the pointer to the last neighbor node of the neighbor list ***

    neighbor_list%last_neighbor_node => new_neighbor_node

!   *** Increment the neighbor node counter ***

    neighbor_list%nnode = neighbor_list%nnode + 1

  END SUBROUTINE add_neighbor_node

! *****************************************************************************

  SUBROUTINE allocate_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Allocate and initialize a set of neighbor lists.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    REAL(wp)                              :: r_max

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

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE clean_neighbor_list_set (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list,next_neighbor_list
    TYPE(neighbor_node_type), POINTER :: neighbor_node,next_neighbor_node
    INTEGER                           :: istat

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
          IF (istat /= 0) CALL stop_memory(routine,"neighbor_node")
          neighbor_node => next_neighbor_node
        END DO

        neighbor_list => next_neighbor_list

      END DO

    END IF

  END SUBROUTINE clean_neighbor_list_set

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_list(neighbor_list)

!   Purpose: Deallocate a neighbor list.

!   History: - Creation (20.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_neighbor_list (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_node_type), POINTER :: neighbor_node,next_neighbor_node
    INTEGER                           :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

      neighbor_node => neighbor_list%first_neighbor_node

      DO WHILE (ASSOCIATED(neighbor_node))
        next_neighbor_node => neighbor_node%next_neighbor_node
        DEALLOCATE (neighbor_node,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine,"neighbor_node")
        neighbor_node => next_neighbor_node
      END DO

      DEALLOCATE (neighbor_list,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"neighbor_list")

    END IF

  END SUBROUTINE deallocate_neighbor_list

! *****************************************************************************

  SUBROUTINE deallocate_neighbor_list_set(neighbor_list_set)

!   Purpose: Deallocate a neighbor list set.

!   History: - Creation (03.11.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_neighbor_list_set (MODULE qs_neighbor_list_types)"

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list,next_neighbor_list
    INTEGER                           :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

      neighbor_list => neighbor_list_set%first_neighbor_list

      DO WHILE (ASSOCIATED(neighbor_list))
        next_neighbor_list => neighbor_list%next_neighbor_list
        CALL  deallocate_neighbor_list(neighbor_list)
        neighbor_list => next_neighbor_list
      END DO

      DEALLOCATE (neighbor_list_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"neighbor_list_set")

    END IF

  END SUBROUTINE deallocate_neighbor_list_set

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

        neighbor_node => neighbor_list%first_neighbor_node

        DO inode=1,nnode
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

  FUNCTION find_neighbor_list_1(neighbor_list_set,atom)&
    RESULT(neighbor_list)

!   Purpose: Return a pointer to the neighbor list of atom "atom" in a
!            neighbor list set. The whole list is traversed.

!   History: - Creation (14.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    INTEGER, INTENT(IN)                   :: atom

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: current_neighbor_list
    INTEGER                           :: ilist

!   ---------------------------------------------------------------------------

    NULLIFY (neighbor_list)

    current_neighbor_list => neighbor_list_set%first_neighbor_list

    DO ilist=1,neighbor_list_set%nlist
      IF (current_neighbor_list%atom == atom) THEN
        neighbor_list => current_neighbor_list
        EXIT
      END IF
      current_neighbor_list => current_neighbor_list%next_neighbor_list
    END DO

  END FUNCTION find_neighbor_list_1

! *****************************************************************************

  FUNCTION find_neighbor_list_2(neighbor_list_set,atom,cell)&
    RESULT(neighbor_list)

!   Purpose: Return a pointer to the neighbor list of atom "atom" in cell
!            "cell" in a neighbor list set. The whole list is traversed.

!   History: - Creation (14.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    INTEGER, INTENT(IN)                   :: atom
    INTEGER, DIMENSION(3), INTENT(IN)     :: cell

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: current_neighbor_list
    INTEGER                           :: ilist

!   ---------------------------------------------------------------------------

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

  END FUNCTION find_neighbor_list_2

! *****************************************************************************

  FUNCTION first_list(neighbor_list_set)&
    RESULT(first_neighbor_list)

!   Purpose: Return a pointer to the first neighbor list of a neighbor list
!            set.

!   History: - Creation (13.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set

    TYPE(neighbor_list_type), POINTER :: first_neighbor_list

!   ---------------------------------------------------------------------------

    first_neighbor_list => neighbor_list_set%first_neighbor_list

  END FUNCTION first_list

! *****************************************************************************

  FUNCTION first_node(neighbor_list)&
    RESULT(first_neighbor_node)

!   Purpose: Return a pointer to the first neighbor node of a neighbor list.

!   History: - Creation (23.06.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    TYPE(neighbor_node_type), POINTER :: first_neighbor_node

!   ---------------------------------------------------------------------------

    first_neighbor_node => neighbor_list%first_neighbor_node

  END FUNCTION first_node

! *****************************************************************************

  SUBROUTINE get_neighbor_list(neighbor_list,atom,cell,nnode)

!   Purpose: Return the reqested data of a neighbor list.

!   History: - Creation (13.09.2000,MK)

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

!   History: - Creation (10.11.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    REAL(wp), OPTIONAL, INTENT(OUT)       :: r_max
    INTEGER, OPTIONAL, INTENT(OUT)        :: nlist

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE get_neighbor_list_set (MODULE qs_neighbor_list_types)"

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN
      IF (PRESENT(r_max)) r_max = neighbor_list_set%r_max
      IF (PRESENT(nlist)) nlist = neighbor_list_set%nlist
    ELSE
      CALL stop_program(routine,&
                        "The requested neighbor list set is not associated")
    END IF

  END SUBROUTINE get_neighbor_list_set

! *****************************************************************************

  SUBROUTINE get_neighbor_node(neighbor_node,neighbor,r,cell)

!   Purpose: Return the reqested data of a neighbor node.

!   History: - Creation (23.06.2000,MK)

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

  SUBROUTINE init_neighbor_list(neighbor_list)

!   Purpose: Initialize a neighbor list. Nothing is (de)allocated here.
!            This routine is also used to prepare a neighbor list for
!            overwriting.

!   History: - Creation (21.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list)) THEN

!     *** Initialize the pointers to the last neighbor node ***

      NULLIFY (neighbor_list%last_neighbor_node)

!     *** Initialize the neighbor list counter ***

      neighbor_list%nnode = 0

    END IF

  END SUBROUTINE init_neighbor_list

! *****************************************************************************

  SUBROUTINE init_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Initialize a neighbor list set. Nothing is (de)allocated here.
!            This routine is also used to prepare a neighbor list set for
!            overwriting.

!   History: - Creation (20.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    REAL(wp), OPTIONAL                    :: r_max

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(neighbor_list_set)) THEN

!     *** Initialize the pointers to the last neighbor list ***

      NULLIFY (neighbor_list_set%last_neighbor_list)

!     *** Initialize the neighbor list counter ***

      neighbor_list_set%nlist = 0

!     *** Initialize the maximum interaction radius (optional) ***

      IF (PRESENT(r_max)) neighbor_list_set%r_max = r_max

    END IF

  END SUBROUTINE init_neighbor_list_set

! *****************************************************************************

  FUNCTION next_neighbor_list(neighbor_list)&
    RESULT (next_list)

!   Purpose: Return the pointer to the next neighbor list.

!   History: - Creation (13.09.2000,MK)

!   ***************************************************************************

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    TYPE(neighbor_list_type), POINTER :: next_list

!   ---------------------------------------------------------------------------

    next_list => neighbor_list%next_neighbor_list

  END FUNCTION next_neighbor_list

! *****************************************************************************

  FUNCTION next_neighbor_node(neighbor_node) RESULT(next_node)

!   Purpose: Return the pointer to the next neighbor node.

!   History: - Creation (23.06.2000,MK)

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

  SUBROUTINE set_neighbor_list_set(neighbor_list_set,r_max)

!   Purpose: Set the components of a neighbor list set.

!   History: - Creation (18.09.2002,MK)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    REAL(wp), OPTIONAL, INTENT(IN)        :: r_max

!   ---------------------------------------------------------------------------

    IF (PRESENT(r_max)) neighbor_list_set%r_max = r_max

  END SUBROUTINE set_neighbor_list_set

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

!!****f* qs_neighbor_list_types/neigh_group_create [1.0] *
!!
!!   NAME
!!     neigh_group_create
!!
!!   SYNOPSIS
!!     Subroutine neigh_group_create(neigh_group, n_kinds, neigh_sets_ptr,&
!!         error)
!!       Type(neighbor_group_type), Pointer:: neigh_group
!!       Integer, Intent (IN):: n_kinds
!!       Type(neighbor_list_set_p_type), Dimension(:), Pointer, Optional::&
!!         neigh_sets_ptr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine neigh_group_create
!!
!!   FUNCTION
!!     creates a neighbor group
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - neigh_group: the neigh_group to create
!!     - n_kinds: number of kinds
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE neigh_group_create(neigh_group, n_kinds, neigh_sets_ptr, error)
  TYPE(neighbor_group_type), POINTER :: neigh_group
  INTEGER, INTENT(in) :: n_kinds
  TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER, OPTIONAL :: neigh_sets_ptr
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='neigh_group_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i

  failure=.FALSE.
  
  ALLOCATE(neigh_group,stat=stat)
  CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

  IF (.NOT. failure) THEN
     last_neigh_group_id=last_neigh_group_id+1
     neigh_group%id_nr=last_neigh_group_id
     neigh_group%ref_count=1
     neigh_group%n_kinds=n_kinds
     IF (PRESENT(neigh_sets_ptr)) THEN
        neigh_group%should_dealloc_neigh_sets=.FALSE.
        neigh_group%neigh_sets => neigh_sets_ptr
     ELSE
        neigh_group%should_dealloc_neigh_sets=.TRUE.
        ALLOCATE(neigh_group%neigh_sets(n_kinds*(n_kinds+1)/2), stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DO i=1,SIZE(neigh_group%neigh_sets)
           NULLIFY(neigh_group%neigh_sets(i)%neighbor_list_set)
        END DO
     END IF
  END IF
END SUBROUTINE neigh_group_create
!***************************************************************************

!!****f* qs_neighbor_list_types/neigh_group_retain [1.0] *
!!
!!   NAME
!!     neigh_group_retain
!!
!!   SYNOPSIS
!!     Subroutine neigh_group_retain(neigh_group, error)
!!       Type(neighbor_group_type), Pointer:: neigh_group
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine neigh_group_retain
!!
!!   FUNCTION
!!     retains a neighbor_group
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - neigh_group: the neighbor group to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE neigh_group_retain(neigh_group,error)
  TYPE(neighbor_group_type), POINTER :: neigh_group
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='neigh_group_retain',&
        routineP=moduleN//':'//routineN
  
  failure=.FALSE.
  CPPrecondition(ASSOCIATED(neigh_group),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(neigh_group%ref_count>0,cp_failure_level,routineP,error)
     neigh_group%ref_count=neigh_group%ref_count+1
  END IF
END SUBROUTINE neigh_group_retain
!***************************************************************************

!!****f* qs_neighbor_list_types/neigh_group_release [1.0] *
!!
!!   NAME
!!     neigh_group_release
!!
!!   SYNOPSIS
!!     Subroutine neigh_group_release(neigh_group, error)
!!       Type(neighbor_group_type), Pointer:: neigh_group
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine neigh_group_release
!!
!!   FUNCTION
!!     releases a neighbor_group
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - neigh_group: the neighbor group to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE neigh_group_release(neigh_group,error)
  TYPE(neighbor_group_type), POINTER :: neigh_group
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='neigh_group_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i
  
  failure=.FALSE.
  IF (ASSOCIATED(neigh_group)) THEN
     CPPreconditionNoFail(neigh_group%ref_count>0,cp_failure_level,routineP,error)
     neigh_group%ref_count=neigh_group%ref_count-1
     IF (neigh_group%ref_count==0) THEN
        neigh_group%ref_count=1
        IF (neigh_group%should_dealloc_neigh_sets.AND.&
             ASSOCIATED(neigh_group%neigh_sets)) THEN
           DO i=1,SIZE(neigh_group%neigh_sets)
              CALL deallocate_neighbor_list_set(&
                   neigh_group%neigh_sets(i)%neighbor_list_set)
           END DO
           DEALLOCATE(neigh_group%neigh_sets, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        neigh_group%ref_count=0
        DEALLOCATE(neigh_group, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(neigh_group)
END SUBROUTINE neigh_group_release
!***************************************************************************

END MODULE qs_neighbor_list_types
