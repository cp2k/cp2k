!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2k developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linklist_control [1.0] *
!!
!!   NAME
!!     linklist_control
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Harald Forbert (Dec-2000): Changes for multiple linked lists
!!                                linklist_internal_data_type
!!
!!   SOURCE
!******************************************************************************

MODULE linklist_control

  USE atomic_kinds, ONLY : kind_info_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE linklist_cell_list, ONLY : cell_list
  USE linklist_cell_types, ONLY : cell_neighborlist_type,  &
       cell_atoms_linklist_type, &
       neighbor_cells_type
  USE linklist_utilities, ONLY : startlist_images_type
  USE linklist_verlet_list, ONLY : verlet_list
  USE linklist_overlap_list, ONLY : overlap_list
  USE molecule_types, ONLY : particle_node_type,linklist_neighbor, &
       linklist_exclusion, linklist_images
  USE message_passing, ONLY : mp_max
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop
  USE linklist_types, ONLY : linklist_internal_data_type
  USE structure_types, ONLY : structure_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: set_ll_parm, get_ll_parm, list_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE set_ll_parm ( internal_data, globenv_in, vskin, natype, rcut, &
     ppcells, count, printlevel, ltype )

  IMPLICIT NONE

! Arguments
  TYPE ( linklist_internal_data_type ), INTENT ( INOUT ) :: internal_data
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

  internal_data%globenv_ll = globenv_in

! Defines the number of cells for the particle-particle lists
  IF ( PRESENT ( ppcells)) internal_data%pp_cell = ppcells

! Defines the verlet skin
  IF ( PRESENT ( vskin)) internal_data%verlet_skin = vskin

! Defines the number of atom types
  IF ( PRESENT ( natype ))  internal_data%natom_types = natype

! Defines cut-off for particle-particle lists
  IF ( PRESENT ( rcut)) THEN

! allocating rlist and rlist_cut to upper triangular part
!    of interaction matrix
     ALLOCATE (internal_data%rlist_cut(internal_data%natom_types, &
               internal_data%natom_types),STAT=isos)
     IF ( isos /= 0 ) CALL stop_memory ( 'set_ll_param', &
          'rlist_cut', internal_data%natom_types ** 2 )

     ALLOCATE ( internal_data%rlist_cutsq ( internal_data%natom_types,  &
                internal_data%natom_types ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'set_ll_param', &
          'rlist_cutsq', internal_data%natom_types ** 2 )
     DO iat = 1, internal_data%natom_types
        DO jat = iat, internal_data%natom_types
           internal_data%rlist_cut(iat,jat) = sqrt(rcut(iat,jat)) + &
              internal_data%verlet_skin
           internal_data%rlist_cutsq(iat,jat) = &
              internal_data%rlist_cut(iat,jat)*internal_data%rlist_cut(iat,jat)
           internal_data%rlist_cut(jat,iat) = internal_data%rlist_cut(iat,jat)
           internal_data%rlist_cutsq(jat,iat) = internal_data%rlist_cutsq(iat,jat)
        END DO
     END DO
  END IF

! Defines the counting variable
  IF ( PRESENT ( count ) ) internal_data%counter = count

! Defines the print level
  IF ( PRESENT ( printlevel ) ) internal_data%print_level = printlevel

! Defines the type of list to be constructed
  IF ( PRESENT ( ltype ) ) THEN
     SELECT CASE ( ltype )
     CASE DEFAULT
        internal_data%list_type = 0
     CASE ( 'NONBOND')
        internal_data%list_type = 1
     CASE ( 'OVERLAP')
        internal_data%list_type = 2
     END SELECT
  END IF

END SUBROUTINE set_ll_parm

!******************************************************************************

SUBROUTINE get_ll_parm ( internal_data, last_up, av_up )

  IMPLICIT NONE

  TYPE ( linklist_internal_data_type ), INTENT ( IN ) :: internal_data
  REAL ( dbl ), INTENT ( OUT ) :: last_up, av_up

  last_up = internal_data%lup
  av_up = internal_data%aup

END SUBROUTINE get_ll_parm

!******************************************************************************

!
! to decide whether the neighbor list is to be updated or not
! based on a displacement criterion;
! if any particle has moved by 0.5*skin_length from the previous
! list_update, then the list routine is called.
!
SUBROUTINE list_control ( internal_data, pnode, part, box, ki )

  IMPLICIT NONE

! Arguments
  TYPE ( linklist_internal_data_type ), INTENT ( INOUT ) :: internal_data
  TYPE ( particle_type ), DIMENSION ( : ), INTENT ( INOUT ) :: part
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( INOUT ) :: pnode
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE (kind_info_type), INTENT (IN), DIMENSION (:), OPTIONAL :: ki

! Locals
  REAL ( dbl ), DIMENSION (3) :: disp_vec
  REAL ( dbl ) :: displace, max_displace
  LOGICAL :: list_update_flag, ionode
  INTEGER :: i, nnodes, isos, iw, handle
  CHARACTER ( LEN = 36 ) :: string

!------------------------------------------------------------------------------

  CALL timeset ( 'LIST_CONTROL','I',' ',handle )

  nnodes = SIZE ( pnode )
  IF ( .NOT. ASSOCIATED ( internal_data%r_last_update ) ) THEN
     ALLOCATE ( internal_data%r_last_update ( 3, nnodes ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'list_control', &
          'r_last_update', 3 * nnodes )

     internal_data%r_last_update = 1.0E+10_dbl
  END IF

! check to see if neighbor list has to be updated
  max_displace = -1.0E+15_dbl
  DO i = 1, nnodes
     disp_vec ( : ) = internal_data%r_last_update(:,i) - pnode(i) %p%r ( : )
     displace = DOT_PRODUCT ( disp_vec, disp_vec )
     max_displace = MAX ( displace, max_displace )
  END DO
#if defined(__parallel)
  CALL mp_max ( max_displace, internal_data%globenv_ll%group )
#endif

  IF ( max_displace >= 0.25_dbl * internal_data%verlet_skin ** 2 ) THEN
     list_update_flag = .TRUE.
     DO i = 1, nnodes
        internal_data%r_last_update(1,i) = pnode(i) % p % r(1)
        internal_data%r_last_update(2,i) = pnode(i) % p % r(2)
        internal_data%r_last_update(3,i) = pnode(i) % p % r(3)
     END DO
  ELSE
     list_update_flag = .FALSE.
  END IF

  IF ( list_update_flag ) THEN
     SELECT CASE ( internal_data%list_type )
     CASE ( 0 )

        CALL stop_program ( 'list_control', 'list type not specified' )

     CASE ( 1 )
        string = ' cell_index: nonbonded neighbor list'

        CALL cell_list(internal_data%natom_types,internal_data%pp_cell,&
             internal_data%rlist_cut,string,pnode,part, &
             box,internal_data%pp_cell_ll,internal_data%pp_startlist_cell,&
             internal_data%pp_neighbor, internal_data%globenv_ll, &
	     internal_data%list_type)

        CALL verlet_list(internal_data%natom_types,internal_data%pp_images,&
             internal_data%pp_cell,pnode,part,box, &
             internal_data%verlet_skin,internal_data%rlist_cut,&
             internal_data%rlist_cutsq,internal_data%pp_neighbor,&
             internal_data%pp_cell_ll,internal_data%pp_startlist_nl,&
             internal_data%pp_startlist_im, internal_data%globenv_ll, &
	     internal_data%list_type)

     CASE ( 2 )

        string = ' cell_index: overlap neighbor list'

        CALL cell_list(internal_data%natom_types,internal_data%pp_cell,&
             internal_data%rlist_cut,string,pnode,part, &
             box,internal_data%pp_cell_ll,internal_data%pp_startlist_cell,&
             internal_data%pp_neighbor, internal_data%globenv_ll, & 
	     internal_data%list_type)

        CALL overlap_list(internal_data%pp_images,&
             internal_data%pp_cell,pnode,part,box, ki, &
             internal_data%verlet_skin,internal_data%rlist_cut,&
             internal_data%rlist_cutsq,internal_data%pp_neighbor,&
             internal_data%pp_cell_ll,internal_data%pp_startlist_nl,&
             internal_data%pp_startlist_im, internal_data%globenv_ll, &
	     internal_data%list_type)

     END SELECT

     iw = internal_data%globenv_ll % scr
     ionode = internal_data%globenv_ll % ionode
     IF ( ionode .AND. internal_data%print_level > 5 ) WRITE (iw,'(a,T71,I10 )' ) &
          ' list_control| Neighbor list updated at step: ', internal_data%counter
     IF ( internal_data%counter > 0 ) THEN
        internal_data%num_update = internal_data%num_update + 1
        internal_data%lup = internal_data%counter + 1 - internal_data%last_update
        internal_data%last_update = internal_data%counter + 1
        internal_data%aup = internal_data%aup + &
            (internal_data%lup-internal_data%aup)/float(internal_data%num_update)
     ELSE
        internal_data%num_update = 0
        internal_data%lup = 0
        internal_data%last_update = 1
        internal_data%aup = 0.0_dbl
     END IF

  END IF

  internal_data%counter = internal_data%counter + 1

  CALL timestop ( 0._dbl, handle )

END SUBROUTINE list_control

!******************************************************************************

END MODULE linklist_control

!******************************************************************************
