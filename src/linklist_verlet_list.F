!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linklist_verlet_list [1.0] *
!!
!!   NAME
!!     linklist_verlet_list
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_verlet_list

  USE kinds, ONLY : dbl

  USE global_types, ONLY : global_environment_type
  USE linklist_cell_types, ONLY : cell_neighborlist_type, &
       neighbor_cells_type, cell_atoms_linklist_type
  USE linklist_utilities, ONLY : get_match, update_verlet_list,  &
       number_of_image_cells, startlist_images_type
  USE molecule_types, ONLY : particle_node_type, linklist_neighbor, &
       linklist_exclusion
  USE message_passing, ONLY : mp_max, mp_min, mp_sum
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory
  USE timings, ONLY : timeset, timestop
  USE util, ONLY : include_list, sort

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: verlet_list

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

! to get the Verlet neighbor list for nonbonded interactions

SUBROUTINE verlet_list(natom_types,n_images,n_cell,pnode,part,box, &
     verlet_skin,rlist_cut,rlist_cutsq,neighbor, &
     cell_ll,start_list_nl,start_list_im, globenv)

  IMPLICIT NONE

!
! neighbor_list lists all the neighbors of all atoms
!

! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ), TARGET :: &
       pnode
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( IN ), TARGET :: part
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( linklist_neighbor), POINTER :: start_list_nl
  TYPE ( startlist_images_type ), POINTER :: start_list_im ( : )
!*apsi*  TYPE ( cell_neighborlist_type ), INTENT ( IN ) :: cell_ll(:,:,:)
!*apsi*  TYPE ( neighbor_cells_type ), INTENT ( IN ) :: neighbor ( :, : )
!*apsi* PGI 3.2 did not like these as arrays
  TYPE ( cell_neighborlist_type ), DIMENSION ( :, :, : ), POINTER :: cell_ll
  TYPE ( neighbor_cells_type ), DIMENSION ( :, : ), POINTER :: neighbor
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  REAL ( dbl ), DIMENSION ( :, : ), INTENT ( IN ) :: rlist_cut, rlist_cutsq
  REAL ( dbl ), INTENT ( IN ) :: verlet_skin
  INTEGER, DIMENSION ( :, :, : ), POINTER :: n_images
  INTEGER, DIMENSION ( : ), INTENT ( IN ) :: n_cell
  INTEGER, INTENT ( IN ) :: natom_types

! Locals
  TYPE (cell_atoms_linklist_type ), POINTER :: ll
  INTEGER :: i, ii, j, natoms, nnodes, handle, ipair
  INTEGER :: iatom, jatom, iexclude, isos
  INTEGER :: nloc, nsum, nmax, nmin, nexcl, ncells
  INTEGER :: k, ipart, jpart, itype, jtype, npairs
  INTEGER, DIMENSION (3) :: index
  INTEGER, PARAMETER :: els = 50
  INTEGER, DIMENSION ( : ), POINTER :: elist, ind
  REAL ( dbl ) :: rijsq, fis
  REAL ( dbl ), DIMENSION (3) :: perd
  REAL ( dbl ), DIMENSION (3,3) :: hmat, h_inv
  LOGICAL :: first_time

!
! the following lines are used for the link_list of
! neighbors
!
  TYPE (linklist_neighbor), POINTER :: current_neighbor
  TYPE (linklist_exclusion), POINTER :: ll_exclude

!------------------------------------------------------------------------------

!
! if we are entering this subroutine for the first time
! create the start of the link_list of neighbors.
! start_list and first_time are both of type save
! so this if statement will not be run or needed
! after the first call to the subroutine
!
  CALL timeset ( 'VERLET_LIST', 'I', 'KIs', handle )
  fis = 0.0_dbl

  ALLOCATE (elist(els),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', 'elist', els )

  ALLOCATE (ind(els),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', 'ind', els )

! local copies of box parameters
  hmat = box%hmat
  h_inv = box%h_inv
  perd = box%perd
  nnodes = size(pnode)
  natoms = size(part)

  first_time = .NOT. ASSOCIATED ( n_images )

  IF ( first_time ) THEN

! compute the number of images for a given cut-off
     ALLOCATE ( n_images ( natom_types, natom_types, 3 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'verlet list', &
          'n_images', 3 * natom_types ** 2 )

     CALL number_of_image_cells ( box, n_images, verlet_skin, rlist_cut )

! allocate the start of the lists
     ALLOCATE (start_list_nl,STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', 'start_list', 0 )
     NULLIFY (start_list_nl%next)

     IF (MAXVAL(n_images) /= 0 ) THEN
        npairs=natoms*(natoms+1)/2
        ALLOCATE(start_list_im(npairs),STAT=isos)
        IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', 'start_list_im', 0 )
        DO ipair=1,npairs
           ALLOCATE(start_list_im(ipair)%image, STAT=isos)
           IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', &
                'start_list_im%image', 0 )
           NULLIFY (start_list_im(ipair)%image%next)
        END DO
     END IF

  END IF

  IF (first_time .AND. globenv % print_level > 4 &
       .OR. globenv % print_level > 9 ) THEN
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T71,i10 )' ) &
          ' verlet_list| # of atoms : ', natoms
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T71,i10 )' ) &
          ' verlet_list| # of atom nodes : ', nnodes
  END IF

!
! let the current_neighbor point to the start of the list
!
  current_neighbor => start_list_nl

!
! loop over all distinct i,j pairs of atoms
!
  ipair=0
  DO i = 1, nnodes
     ii = pnode(i) %p%iatom
     itype = pnode(i) %p%prop%ptype
     pnode(i) %nneighbor = 0
     nexcl = pnode(i) %nexcl
     IF (nexcl>size(elist)) STOP 'els'
     ll_exclude => pnode(i) %ex

! list number of exclusion
     DO iexclude = 1, nexcl
        elist(iexclude) = ll_exclude%p%iatom
        ll_exclude => ll_exclude%next
     END DO

! sort elist
     CALL sort(elist,nexcl,ind)

!
! pnode(i)%nl is the start of the neighbor link_list for
! atom i. Let it point to the current position in
! the link_list
!
     pnode(i) %nl => current_neighbor

! first: loop over all particles in the same cell
     DO jtype = 1, natom_types
        ll => part(ii) % cell % pTYPE ( jtype ) % ll
        DO jpart = 1, part(ii) % cell % natoms(jtype )
           j = ll%atom
! use black/white scheme to avoid double counting
           IF (include_list(ii,j)) THEN
              ipair=ipair+1
              IF (MAXVAL(n_images(itype,jtype,:)) /= 0 )  &
                   current_neighbor%image=>start_list_im(ipair)%image
! look for a match in the exclusion list (j in excl of pnode(i))
              SELECT CASE (nexcl)
              CASE DEFAULT
                 IF ( .NOT. get_match(j,elist,nexcl)) THEN
                    CALL update_verlet_list(n_images ( itype, jtype, : ), &
                         part(j), j, &
                         pnode(i),hmat,h_inv,perd,rlist_cutsq(itype,jtype ), &
                         current_neighbor)
                 END IF
              CASE (0)
                 CALL update_verlet_list(n_images(itype,jtype,:),part(j),j, &
                      pnode(i),hmat,h_inv,perd,rlist_cutsq(itype,jtype ), &
                      current_neighbor)
              CASE (1)
                 IF ( j/=elist(1)) THEN
                    CALL update_verlet_list(n_images(itype,jtype,:), &
                         part(j),j, &
                         pnode(i),hmat,h_inv,perd,rlist_cutsq(itype,jtype ), &
                         current_neighbor)
                 END IF
              END SELECT
           END IF
           ll => ll%next
        END DO
     END DO

! second: loop over all neighbor cells
     DO jtype = 1, natom_types
        IF (ALL(n_cell==1)) EXIT
        ncells=SIZE (neighbor(itype,jtype ) % index, 2)
        DO k = 1, ncells
           index(1) = part(ii) % cell % neighbor(itype,jtype ) %index(1,k)
           index(2) = part(ii) % cell % neighbor(itype,jtype ) %index(2,k)
           index(3) = part(ii) % cell % neighbor(itype,jtype ) %index(3,k)
           ll => cell_ll(index(1),index(2),index(3)) % pTYPE ( jtype ) %ll
           DO jpart = 1, cell_ll(index(1),index(2),index(3)) %natoms( &
                jtype )
              ipair=ipair+1
              j = ll%atom
              SELECT CASE (nexcl)
              CASE DEFAULT
                 IF ( .NOT. get_match(j,elist,nexcl)) THEN
                    CALL update_verlet_list(n_images(itype,jtype,:), &
                         part(j),j, &
                         pnode(i),hmat,h_inv,perd,rlist_cutsq(itype,jtype ), &
                         current_neighbor)
                 END IF
              CASE (0)
                 CALL update_verlet_list(n_images(itype,jtype,:),part(j),j, &
                      pnode(i),hmat,h_inv,perd,rlist_cutsq(itype,jtype ), &
                      current_neighbor)
              CASE (1)
                 IF ( j/=elist(1)) THEN
                    CALL update_verlet_list(n_images(itype,jtype,:), &
                         part(j),j, &
                         pnode(i),hmat,h_inv,perd,rlist_cutsq(itype,jtype ), &
                         current_neighbor)
                 END IF
              END SELECT
              ll => ll%next
           END DO
        END DO
     END DO

     IF ( first_time .AND. globenv % print_level>4 &
          .OR. globenv % print_level>9) THEN
        IF (globenv % ionode) WRITE (globenv % scr,'(A,i10,T54,A,T71,I10 )' ) &
             ' verlet_list| Atom :', i, ' # of neighbors:', &
             pnode(i) %nneighbor
     END IF

  END DO

!
! write some info on the neighbor lists to output
!
  nloc = sum(pnode ( : )%nneighbor)
  IF ( first_time .AND. globenv % print_level > 4 &
       .OR. globenv % print_level > 9 ) THEN

#if defined(__parallel)
     nsum = nloc
     CALL mp_sum(nsum,globenv%group)
     nmax = nloc
     CALL mp_max(nmax,globenv%group)
     nmin = nloc
     CALL mp_min(nmin,globenv%group)
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T69,I12 )' ) &
          ' verlet_list| Total number of interactions ', nsum
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T69,I12 )' ) &
          ' verlet_list| Minimum number of interactions per node', nmin
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T69,I12 )' ) &
          ' verlet_list| Maximum number of interactions per node', nmax
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T69,F10.2,A )' ) &
          ' verlet_list| Load balancing error', 100.0_dbl &
          * REAL ( nmax - nmin, dbl) &
          / REAL ( nsum, dbl ) * globenv % num_pe, ' %'
#else
     IF (globenv % ionode) WRITE (globenv % scr,'(A,T69,I12 )' ) &
          ' verlet_list| Total number of interactions ', nloc
#endif
  END IF
  fis = nloc * 1.0E-3_dbl

  CALL timestop(fis,handle)

  DEALLOCATE (elist,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', 'elist' )

  DEALLOCATE (ind,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'verlet_list', 'ind' )

END SUBROUTINE verlet_list

!******************************************************************************

END MODULE linklist_verlet_list
