!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/nl_utilities [1.0] *
!!
!!   NAME
!!     nl_utilities
!!
!!   FUNCTION
!!     Utility routines for the neighbor list generation.
!!
!!   AUTHOR
!!     Matthias Krack (28.08.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE nl_utilities

! *****************************************************************************

! Index:

! SUBROUTINE build_neighbor_list(neighbor_list_set,atom_a,cell_a,radius_a,&
!                                kind_b,radius_b)
! SUBROUTINE write_neighbor_list_set(neighbor_list_set,lunit)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE atoms,                 ONLY: atom_info,natom
  USE cell_parameters,       ONLY: abc,unit_of_length,unit_of_length_name
  USE atomic_kinds,          ONLY: kind_info_type
  USE method_specifications, ONLY: periodic
  USE neighbor_list_types,   ONLY: add_neighbor_list,&
                                   add_neighbor_node,&
                                   find_neighbor_list,&
                                   first_neighbor_list,&
                                   first_neighbor_node,&
                                   get_neighbor_list,&
                                   get_neighbor_list_set_info,&
                                   get_neighbor_node,&
                                   neighbor_list_set_type,&
                                   neighbor_list_type,&
                                   neighbor_node_type,&
                                   next_neighbor_list,&
                                   next_neighbor_node

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: build_neighbor_list,&
            write_neighbor_list_set

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_neighbor_list(neighbor_list_set,atom_a,cell_a,radius_a,&
                                 kind_b,radius_b)

!   Purpose: Build the neighbor list of atom "atom_a" in cell "cell_a" with all
!            the potential neighbors of the atomic kind "kind_b" based on the
!            specified radii.

!   History: - Creation (10.11.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    TYPE(kind_info_type)                  :: kind_b
    REAL(wp), INTENT(IN)                  :: radius_a,radius_b
    INTEGER, INTENT(IN)                   :: atom_a
    INTEGER, DIMENSION(3), INTENT(IN)     :: cell_a

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list

    REAL(wp) :: rab_max,rab2,rab2_max
    INTEGER  :: atom_b,iatom,icell,jcell,kcell

    REAL(wp), DIMENSION(3) :: ra,rab,rb
    INTEGER, DIMENSION(3)  :: ncell

!   ---------------------------------------------------------------------------

!   *** Find the neighbor list of atom "atom_a" in cell "cell_a" ***

    neighbor_list => find_neighbor_list(neighbor_list_set,atom_a,cell_a)

!   *** Create a new neighbor list, if the neighbor list was not found ***

    IF (.NOT.ASSOCIATED(neighbor_list)) THEN
      CALL add_neighbor_list(neighbor_list_set=neighbor_list_set,&
                             atom=atom_a,&
                             cell=cell_a,&
                             neighbor_list=neighbor_list)
    END IF

!   *** Calculate the position vector of atom "atom_a" in cell "cell_a" ***

    ra(:) = atom_info(atom_a)%r_pbc(:) + REAL(cell_a(:),wp)*abc(:)

!   *** Find all neighbor atoms of atom "atom_a" ***

    IF (periodic) THEN

      DO iatom=1,kind_b%natom

        atom_b = kind_b%atom_list(iatom)

!       *** Calculate the position vector of atom "atom_b" in cell "cell_a" ***

        rb(:) = atom_info(atom_b)%r_pbc(:)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius_a + radius_b
        rab2_max = rab_max*rab_max

!       *** Calculate the number of cells in each direction ***
!       *** scanned for potential neighbor atoms            ***

        ncell(:) = CEILING(rab_max/abc(:))

!       *** Loop over all neighbor cells ***

        DO icell=cell_a(1)-ncell(1),cell_a(1)+ncell(1)
          DO jcell=cell_a(2)-ncell(2),cell_a(2)+ncell(2)
            DO kcell=cell_a(3)-ncell(3),cell_a(3)+ncell(3)

!             *** Calculate the distance vector between two atoms ***

              rab(:) = rb(:) - ra(:) + REAL((/icell,jcell,kcell/),wp)*abc(:)

!             *** Calculate the square of the distance vector ***

              rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

!             *** Check, if the two atoms are neighbors ***

              IF (rab2 < rab2_max) THEN

!               *** Create a new neighbor list entry, if the    ***
!               *** two atoms fulfil the neighborhood criterion ***

                CALL add_neighbor_node(neighbor_list=neighbor_list,&
                                       neighbor=atom_b,&
                                       r2=rab2,&
                                       r=rab,&
                                       cell=(/icell,jcell,kcell/))

              END IF

            END DO
          END DO
        END DO

      END DO

    ELSE

      DO iatom=1,kind_b%natom

        atom_b = kind_b%atom_list(iatom)

!       *** Calculate the position vector of atom "atom_b" ***

        rb(:) = atom_info(atom_b)%r_pbc(:)

!       *** Calculate the square of the maximum interaction distance ***

        rab_max = radius_a + radius_b
        rab2_max = rab_max*rab_max

!       *** Calculate the distance vector between the two atoms ***

        rab(:) = rb(:) - ra(:)

!       *** Calculate the square of the interatomic distance ***

        rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

!       *** Check, if the two atoms are neighbors ***

        IF (rab2 < rab2_max) THEN

!         *** Create a new neighbor list entry, if the    ***
!         *** two atoms fulfil the neighborhood criterion ***

          CALL add_neighbor_node(neighbor_list=neighbor_list,&
                                 neighbor=atom_b,&
                                 r2=rab2,&
                                 r=rab,&
                                 cell=(/0,0,0/))

        END IF

      END DO

    END IF

  END SUBROUTINE build_neighbor_list

! *****************************************************************************

  SUBROUTINE write_neighbor_list_set(neighbor_list_set,lunit)

!   Purpose: Write a set of neighbor lists to the logical unit number lunit.

!   History: - Creation (28.08.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(neighbor_list_set_type), POINTER :: neighbor_list_set
    INTEGER, INTENT(IN)                   :: lunit

!   *** Local variables ***

    TYPE(neighbor_list_type), POINTER :: neighbor_list
    TYPE(neighbor_node_type), POINTER :: neighbor_node

    CHARACTER(LEN=60) :: name
    REAL(wp)          :: rab2
    INTEGER           :: iatom,neighbor,nneighbor

    REAL(wp), DIMENSION(3) :: ra,rab,rb
    INTEGER, DIMENSION(3)  :: cell

!   ---------------------------------------------------------------------------

    CALL get_neighbor_list_set_info(neighbor_list_set=neighbor_list_set,&
                                    name=name)

!   *** Write headline ***

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/,/,T3,A,7X,A,2(11X,A),10X,A)")&
      TRIM(name)//" NEIGHBOR LIST IN "//TRIM(unit_of_length_name),&
      "Atom  Neighbor  Cell(i,j,k)","X","Y","Z","Distance"

!   *** Loop over all atoms and their corresponding neighbor lists ***

    neighbor_list => first_neighbor_list(neighbor_list_set)

    DO WHILE (ASSOCIATED(neighbor_list))

      CALL get_neighbor_list(neighbor_list=neighbor_list,&
                             atom=iatom,&
                             cell=cell,&
                             nnode=nneighbor)

      ra(:) = atom_info(iatom)%r_pbc(:) + REAL(cell(:),wp)*abc(:)

      WRITE (UNIT=lunit,FMT="(/,T2,I5,3X,I6,2X,3I4,3F12.6)")&
        iatom,nneighbor,cell(:),ra(:)/unit_of_length

!     *** Direct the work pointer to the start point of the current list ***

      neighbor_node => first_neighbor_node(neighbor_list)

!     *** Traverse the neighbor list of the current ***
!     *** atom and print the stored information     ***

      DO WHILE (ASSOCIATED(neighbor_node))

        CALL get_neighbor_node(neighbor_node=neighbor_node,&
                               neighbor=neighbor,&
                               r=rab,&
                               r2=rab2,&
                               cell=cell)

        rb(:) = ra(:) + rab(:)

        WRITE (UNIT=lunit,FMT="(T10,I6,2X,3I4,3F12.6,2X,F12.6)")&
          neighbor,cell(:),rb(:)/unit_of_length,SQRT(rab2)/unit_of_length

        neighbor_node => next_neighbor_node(neighbor_node)

      END DO

      neighbor_list => next_neighbor_list(neighbor_list)

    END DO

  END SUBROUTINE write_neighbor_list_set

! *****************************************************************************

END MODULE nl_utilities
