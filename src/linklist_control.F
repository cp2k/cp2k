!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2k developers group                                 !
!-----------------------------------------------------------------------------!

MODULE linklist_control
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE linklist_cell_list, ONLY : cell_list
  USE linklist_cell_types, ONLY : cell_neighborlist_type,  &
       cell_atoms_linklist_type, &
       neighbor_cells_type
  USE linklist_utilities, ONLY : startlist_images_type
  USE linklist_verlet_list, ONLY : verlet_list
  USE molecule_types, ONLY : particle_node_type,linklist_neighbor, &
       linklist_exclusion, linklist_images
  USE message_passing, ONLY : mp_max
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type
  USE stop_program, ONLY : stop_prg, stop_memory
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: set_ll_parm, get_ll_parm, list_control
  
! globals specific to the particle-particle verlet list
  TYPE ( cell_neighborlist_type ), DIMENSION (:,:,:), POINTER :: pp_cell_ll
  TYPE ( cell_atoms_linklist_type ), DIMENSION (:,:,:,:), POINTER  :: &
       pp_startlist_cell
  TYPE ( neighbor_cells_type ), DIMENSION ( :, : ), POINTER :: pp_neighbor
  TYPE ( linklist_neighbor ), POINTER :: pp_startlist_nl
  TYPE ( startlist_images_type ), POINTER :: pp_startlist_im ( : )
  INTEGER, DIMENSION ( 3 ) :: pp_cell
  INTEGER, POINTER :: pp_images(:,:,:)
  INTEGER, POINTER :: pp_ncell ( : )
  
  TYPE ( global_environment_type ) :: globenv_ll
  INTEGER :: natom_types=0
  INTEGER :: list_type = 0
  INTEGER :: counter = 0, last_update = 0, num_update = 0
  INTEGER :: print_level = 1
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rlist_cut, rlist_cutsq
  REAL ( dbl ) :: verlet_skin
  REAL ( dbl ) :: lup, aup
  
!print_level >  0 : no output
!print_level <  5 : normal output
!print_level < 10 : initial neighbor list info
!print_level > 10 : full neighbor list info at each update
  
CONTAINS

!******************************************************************************
  
SUBROUTINE set_ll_parm ( globenv_in, vskin, natype, rcut, &
     ppcells, count, printlevel, ltype )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv_in
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: vskin
  REAL ( dbl ), INTENT ( IN ), OPTIONAL :: rcut ( :, : )
  INTEGER, INTENT ( IN ), DIMENSION ( : ), OPTIONAL :: ppcells
  INTEGER, INTENT ( IN ), OPTIONAL :: natype
  INTEGER, INTENT ( IN ), OPTIONAL :: count
  INTEGER, INTENT ( IN ), OPTIONAL :: printlevel
  CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: ltype
  
! Locals
  REAL ( dbl ) :: cut_off
  INTEGER :: iat, jat
  INTEGER :: isos
  
!------------------------------------------------------------------------------
  
  globenv_ll = globenv_in
  
! Defines the number of cells for the particle-particle lists
  IF ( PRESENT ( ppcells)) pp_cell = ppcells
  
! Defines the verlet skin
  IF ( PRESENT ( vskin)) verlet_skin = vskin
  
! Defines the number of atom types
  IF ( PRESENT ( natype ))  natom_types = natype
  
! Defines cut-off for particle-particle lists
  IF ( PRESENT ( rcut)) THEN
     
! allocating rlist and rlist_cut to upper triangular part
!    of interaction matrix
     ALLOCATE (rlist_cut(natom_types,natom_types),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'set_ll_param', &
          'rlist_cut', natom_types ** 2 )
     
     ALLOCATE ( rlist_cutsq ( natom_types, natom_types ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'set_ll_param', &
          'rlist_cutsq', natom_types ** 2 )
     DO iat = 1, natom_types
        DO jat = iat, natom_types
           rlist_cut(iat,jat) = sqrt(rcut(iat,jat)) + verlet_skin
           rlist_cutsq(iat,jat) = rlist_cut(iat,jat)*rlist_cut(iat,jat)
           rlist_cut(jat,iat) = rlist_cut(iat,jat)
           rlist_cutsq(jat,iat) = rlist_cutsq(iat,jat)
        END DO
     END DO
  END IF
  
! Defines the counting variable
  IF ( PRESENT ( count ) ) counter = count
  
! Defines the print level
  IF ( PRESENT ( printlevel ) ) print_level = printlevel
  
! Defines the type of list to be constructed
  IF ( PRESENT ( ltype ) ) THEN
     SELECT CASE ( ltype )
     CASE DEFAULT
        list_type = 0
     CASE ( 'NONBOND')
        list_type = 1
     END SELECT
  END IF
  
END SUBROUTINE set_ll_parm

!******************************************************************************

SUBROUTINE get_ll_parm ( last_up, av_up )
  
  IMPLICIT NONE
  
  REAL ( dbl ), INTENT ( OUT ) :: last_up, av_up
  
  last_up = lup
  av_up = aup
  
END SUBROUTINE get_ll_parm

!******************************************************************************

!
! to decide whether the neighbor list is to be updated or not
! based on a displacement criterion; 
! if any particle has moved by 0.5*skin_length from the previous 
! list_update, then the list routine is called. 
!
SUBROUTINE list_control ( pnode, part, box )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ) :: box
  
! Locals
  REAL ( dbl ), DIMENSION (3) :: disp_vec
  REAL ( dbl ), ALLOCATABLE, SAVE :: r_last_update ( :, : )
  REAL ( dbl ) :: displace, max_displace
  LOGICAL :: list_update_flag, ionode
  INTEGER :: i, nnodes, isos, iw
  CHARACTER ( LEN = 36 ) :: string
  
!------------------------------------------------------------------------------
  
  nnodes = SIZE ( pnode )
  IF ( .NOT. ALLOCATED ( r_last_update ) ) THEN
     ALLOCATE ( r_last_update ( 3, nnodes ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'list_control', &
          'r_last_update', 3 * nnodes )
     
     r_last_update = 1.0E+10_dbl
  END IF
  
! check to see if neighbor list has to be updated
  max_displace = -1.0E+15_dbl
  DO i = 1, nnodes
     disp_vec ( : ) = r_last_update(:,i) - pnode(i) %p%r ( : )
     displace = DOT_PRODUCT ( disp_vec, disp_vec )
     max_displace = MAX ( displace, max_displace )
  END DO
#if defined(__parallel)
  CALL mp_max ( max_displace, globenv_ll%group )
#endif
  
  IF ( max_displace >= 0.25_dbl * verlet_skin ** 2 ) THEN
     list_update_flag = .TRUE.
     DO i = 1, nnodes
        r_last_update(1,i) = pnode(i) % p % r(1)
        r_last_update(2,i) = pnode(i) % p % r(2)
        r_last_update(3,i) = pnode(i) % p % r(3)
     END DO
  ELSE
     list_update_flag = .FALSE.
  END IF
  
  IF ( list_update_flag ) THEN
     SELECT CASE ( list_type )
     CASE ( 0 )
        
        CALL stop_prg ( 'list_control', 'list type not specified' )
        
     CASE ( 1 )
        string = ' cell_index: nonbonded neighbor list'
        
        CALL cell_list(natom_types,pp_cell,rlist_cut,string,pnode,part, &
             box,pp_cell_ll,pp_startlist_cell,pp_neighbor, globenv_ll)
        
        CALL verlet_list(natom_types,pp_images,pp_cell,pnode,part,box, &
             verlet_skin,rlist_cut,rlist_cutsq,pp_neighbor,pp_cell_ll, &
             pp_startlist_nl,pp_startlist_im, globenv_ll)
     END SELECT
     
     iw = globenv_ll % scr
     ionode = globenv_ll % ionode
     IF ( ionode .AND. print_level > 5 ) WRITE (iw,'(a,T71,I10 )' ) &
          ' list_control| Neighbor list updated at step: ', counter
     IF ( counter > 0 ) THEN
        num_update = num_update + 1
        lup = counter + 1 - last_update
        last_update = counter + 1
        aup = aup + (lup-aup)/float(num_update)
     ELSE
        num_update = 0
        lup = 0
        last_update = 1
        aup = 0.0_dbl
     END IF
     
  END IF
  
  counter = counter + 1
  
END SUBROUTINE list_control

!******************************************************************************

END MODULE linklist_control
