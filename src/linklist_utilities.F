!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE linklist_utilities
  
  USE kinds, ONLY : dbl
  USE molecule_types, ONLY : particle_node_type, linklist_neighbor, &
       linklist_images
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : get_cell_param, cell_type
  USE stop_program, ONLY : stop_memory
  
  IMPLICIT NONE

  PUBLIC :: get_match, update_verlet_list, &
            number_of_image_cells, startlist_images_type
  PRIVATE

  TYPE startlist_images_type
    TYPE ( linklist_images), POINTER :: image
  END TYPE startlist_images_type


  CONTAINS

!------------------------------------------------------------------------------!
  FUNCTION get_match(jatom,elist,nexcl) RESULT (match)
    IMPLICIT NONE
    INTEGER, INTENT ( IN ) :: jatom
    INTEGER, INTENT ( IN ) :: elist ( : )
    INTEGER, INTENT ( IN ) :: nexcl
    LOGICAL :: match
    INTEGER :: iexclude

    match = .FALSE.
    DO iexclude = 1, nexcl
      IF (elist(iexclude)>jatom) EXIT
      IF (elist(iexclude)==jatom) THEN
        match = .TRUE.
        EXIT
      END IF
    END DO
  END FUNCTION get_match
!------------------------------------------------------------------------------!
  SUBROUTINE update_verlet_list(n_images,part,j,pnode,hmat,h_inv,perd, &
                          distsq,current_neighbor)
    IMPLICIT NONE
    TYPE (particle_type ), INTENT ( IN ), TARGET :: part
    TYPE (particle_node_type ), INTENT ( INOUT ) :: pnode
    INTEGER, INTENT ( IN ) :: j
    INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n_images
    REAL ( dbl ), INTENT ( IN ) :: hmat(3,3), h_inv(3,3), perd(3)
    REAL ( dbl ), INTENT ( IN ) :: distsq
    TYPE (linklist_neighbor), POINTER :: current_neighbor
! locals
    REAL ( dbl ), DIMENSION (3) :: rij, s
    REAL ( dbl ) :: rijsq
    INTEGER :: isos

    rij(1) = part % r(1) - pnode % p % r(1)
    rij(2) = part % r(2) - pnode % p % r(2)
    rij(3) = part % r(3) - pnode % p % r(3)
! inline code for general function pbc
    s(1) = h_inv(1,1)*rij(1) + h_inv(1,2)*rij(2) + h_inv(1,3)*rij(3)
    s(2) = h_inv(2,1)*rij(1) + h_inv(2,2)*rij(2) + h_inv(2,3)*rij(3)
    s(3) = h_inv(3,1)*rij(1) + h_inv(3,2)*rij(2) + h_inv(3,3)*rij(3)
    IF (MAXVAL(n_images)==0) THEN
      IF(ABS(s(1)) > 0.5_dbl) THEN
        s(1) = s(1) - perd(1)*INT(s(1)+SIGN(0.5_dbl,s(1)))
      END IF
      IF(ABS(s(2)) > 0.5_dbl) THEN
        s(2) = s(2) - perd(2)*INT(s(2)+SIGN(0.5_dbl,s(2)))
      END IF
      IF(abs(s(3)) > 0.5_dbl) THEN
        s(3) = s(3) - perd(3)*INT(s(3)+sign(0.5_dbl,s(3)))
      END IF
      rij(1) = hmat(1,1)*s(1) + hmat(1,2)*s(2) + hmat(1,3)*s(3)
      rij(2) = hmat(2,1)*s(1) + hmat(2,2)*s(2) + hmat(2,3)*s(3)
      rij(3) = hmat(3,1)*s(1) + hmat(3,2)*s(2) + hmat(3,3)*s(3)
      rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
    END IF
!
!
! if we are within the cutoff, let the current_neighor
! position in the link_list point to the jth atom
!
    IF (rijsq<=distsq.OR.MAXVAL(n_images) /= 0 ) THEN
      pnode%nneighbor = pnode%nneighbor + 1
      current_neighbor%p => part
      current_neighbor%index = j
      current_neighbor%nimages = 0
!
      IF (MAXVAL(n_images) /= 0 ) CALL find_images(n_images,s,hmat, &
                                 perd,distsq,current_neighbor)
!
! advance the current position in the link_link by one.
! however, since we are constantly reusing this list
! and the number of neighbors may fluctuate, we have to check
! to see that the next member of the list exists. If not, create it
! and then advance by one.
      IF ( .NOT. associated(current_neighbor%next)) THEN
        ALLOCATE (current_neighbor%next,STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'update_verlet_list', &
             'current_neighbor%next', 0 )
        current_neighbor => current_neighbor%next
        NULLIFY (current_neighbor%next)
      ELSE
        current_neighbor => current_neighbor%next
      END IF
    END IF

  END SUBROUTINE update_verlet_list
!------------------------------------------------------------------------------!
  SUBROUTINE update_image_list(sij,n,vec,hmat,perd,distsq,current_image)
    IMPLICIT NONE
    INTEGER, INTENT ( INOUT ) :: n
    REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: sij, perd, vec
    REAL ( dbl ), INTENT ( IN ), DIMENSION ( :, : ) :: hmat
    REAL ( dbl ), INTENT ( IN ) :: distsq
    TYPE (linklist_images), POINTER :: current_image
! locals
    REAL ( dbl ), DIMENSION (3) :: rij, s
    REAL ( dbl ) :: rijsq
    INTEGER :: isos
! add the lattice tranlation
    s(1) = sij(1) + perd(1)*vec(1)
    s(2) = sij(2) + perd(2)*vec(2)
    s(3) = sij(3) + perd(3)*vec(3)
! compute distance
    rij(1) = hmat(1,1)*s(1) + hmat(1,2)*s(2) + hmat(1,3)*s(3)
    rij(2) = hmat(2,1)*s(1) + hmat(2,2)*s(2) + hmat(2,3)*s(3)
    rij(3) = hmat(3,1)*s(1) + hmat(3,2)*s(2) + hmat(3,3)*s(3)
!
    rijsq = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
!
! if we are within the cutoff, let the current_image
! position in the link_list point to the jth atom
!
    IF (rijsq<=distsq) THEN
      n = n + 1
      current_image%vec = vec
!
! advance the current position in the link_link by one.
! however, since we are constantly reusing this list
! and the number of images may fluctuate, we have to check
! to see that the next member of the list exists. If not, create it
! and then advance by one.
      IF ( .NOT. ASSOCIATED ( current_image % next ) ) THEN
        ALLOCATE ( current_image % next, STAT = isos )
        IF ( isos /= 0 ) CALL stop_memory ( 'update_verlet_list', &
             'current_image%next', 0 )
        current_image => current_image % next
        NULLIFY ( current_image % next )
      ELSE
         current_image => current_image % next
      END IF
    END IF

  END SUBROUTINE update_image_list
!------------------------------------------------------------------------------!
  SUBROUTINE find_images(n_images,s,hmat,perd,distsq,current_neighbor)
    IMPLICIT NONE
    REAL ( dbl ), INTENT ( IN ), DIMENSION ( : ) :: s
    REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: hmat
    REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: perd
    REAL ( dbl ), INTENT ( IN ) :: distsq
    INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n_images
    TYPE ( linklist_neighbor), POINTER :: current_neighbor
! locals
    TYPE ( linklist_images), POINTER :: current_image
    INTEGER :: icell,jcell,kcell
    REAL ( dbl ), DIMENSION(3) :: vec

    current_image=>current_neighbor%image

    DO kcell=-n_images(3), n_images(3)
      DO jcell=-n_images(2), n_images(2)
        DO icell=-n_images(1), n_images(1)
! exclude the center cell
          IF (icell==0.AND.jcell==0.AND.kcell==0) CYCLE
 
          vec ( : )=(/icell,jcell,kcell/)

          CALL update_image_list(s,current_neighbor%nimages,vec, &
                                 hmat,perd,distsq,current_image)

        END DO
      END DO
    END DO
  END SUBROUTINE find_images
!------------------------------------------------------------------------------!
  SUBROUTINE number_of_image_cells(box,n_images,verlet_skin,list_cut)
    IMPLICIT NONE
    TYPE ( cell_type ), INTENT ( IN ) :: box
    REAL ( dbl ),DIMENSION ( :, : ), INTENT ( IN ) :: list_cut
    INTEGER, DIMENSION(:,:,:), INTENT ( OUT ) :: n_images
    REAL ( dbl ), INTENT ( IN ) :: verlet_skin
! locals
    INTEGER :: itype, jtype, id, ntypes
    REAL ( dbl ), DIMENSION(3) :: quotient, cell_lengths
! get cell parameters
    CALL get_cell_param(box,cell_lengths)
! compute the maximum translation vector for each type
    ntypes=SIZE(list_cut,1)
    DO itype=1, ntypes
      DO jtype=1, ntypes
        quotient ( : )=(list_cut(itype,jtype )-verlet_skin)/cell_lengths ( : )
        DO id=1,3
          IF (quotient(id)<=.5_dbl) THEN
            n_images(itype,jtype,id)=0
          ELSE
            n_images(itype,jtype,id) =  &
            CEILING(list_cut(itype,jtype )/cell_lengths(id))
          END IF
        END DO 
      END DO
    END DO
  END SUBROUTINE number_of_image_cells
!------------------------------------------------------------------------------!
END MODULE linklist_utilities
