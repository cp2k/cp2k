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
  USE atomic_kind_types, ONLY : atomic_kind_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE fist_environment_types, ONLY : replica_environment_type
  USE linklist_verlet_list, ONLY : build_verlet_lists, rebuild_verlet_lists
  USE linklist_types, ONLY : linklist_internal_data_type
  USE molecule_types, ONLY : particle_node_type
  USE message_passing, ONLY : mp_max
  USE particle_types, ONLY : particle_type
  USE simulation_cell, ONLY : cell_type, real_to_scaled, scaled_to_real
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop

  USE fist_neighbor_lists, ONLY: build_fist_neighbor_lists,&
                                 rebuild_fist_neighbor_lists

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: list_control

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!
! to decide whether the neighbor list is to be updated or not
! based on a displacement criterion;
! if any particle has moved by 0.5*skin_length from the previous
! list_update, then the list routine is called.
!
SUBROUTINE list_control ( rep_env )

  IMPLICIT NONE

! Arguments
  TYPE ( replica_environment_type ), INTENT ( INOUT ) :: rep_env

! Locals
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( linklist_internal_data_type ), POINTER :: internal_data
  TYPE ( particle_type ), POINTER :: part ( : )
  TYPE ( particle_node_type ), POINTER :: pnode ( : )
  TYPE ( cell_type ), POINTER :: box
  REAL ( dbl ), DIMENSION (3) :: disp_vec
  REAL ( dbl ) :: displace, max_displace
  LOGICAL :: list_update_flag, ionode
  LOGICAL :: first_time
  INTEGER :: i, nnodes, isos, iw, handle
  CHARACTER ( LEN = 36 ) :: string

!------------------------------------------------------------------------------

  CALL timeset ( 'LIST_CONTROL','I',' ',handle )

! Assigning local pointers

  part => rep_env % part 
  pnode => rep_env % pnode 
  box =>  rep_env % box
  internal_data => rep_env % ll_data ( 1 )
  atomic_kind_set => rep_env % atomic_kind_set 
  
!MK  first_time = .NOT.ASSOCIATED ( internal_data % sab_ppnl )
  first_time = .NOT.ASSOCIATED ( internal_data % nonbonded )

  nnodes = SIZE ( pnode )
  IF ( .NOT. ASSOCIATED ( internal_data % r_last_update ) ) THEN
     ALLOCATE ( internal_data%r_last_update ( 3, nnodes ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'list_control', &
          'r_last_update', 3 * nnodes )

     internal_data % r_last_update = 1.0E+10_dbl
  END IF

! check to see if neighbor list has to be updated
  max_displace = -1.0E+15_dbl
  DO i = 1, nnodes
     disp_vec ( : ) = scaled_to_real ( internal_data % r_last_update ( :, i ), box ) -  &
                      pnode ( i ) % p % r ( : )
     displace = DOT_PRODUCT ( disp_vec, disp_vec )
     max_displace = MAX ( displace, max_displace )
  END DO
#if defined(__parallel)
  CALL mp_max ( max_displace, internal_data % globenv % group )
#endif

  IF ( max_displace >= 0.25_dbl * internal_data%verlet_skin ** 2 ) THEN
     list_update_flag = .TRUE.
     DO i = 1, nnodes
        internal_data % r_last_update ( :, i ) = real_to_scaled ( pnode ( i ) % p % r ( : ), box )
      !  internal_data % r_last_update ( 1, i ) =  pnode ( i ) % p % r ( 1 ), box
      !  internal_data % r_last_update ( 2, i ) =  pnode ( i ) % p % r ( 2 ), box
      !  internal_data % r_last_update ( 3, i ) =  pnode ( i ) % p % r ( 3 ), box
     END DO
  ELSE
     list_update_flag = .FALSE.
  END IF

  IF ( list_update_flag ) THEN
!  IF ( .TRUE. ) THEN

     IF ( first_time ) THEN
       CALL build_fist_neighbor_lists(atomic_kind_set,part,pnode,box,&
                                      internal_data%rlist_cut,&
                                      internal_data%nonbonded,&
                                      internal_data%globenv)
!MK       CALL build_verlet_lists (atomic_kind_set, part, box, internal_data, pnode ) 
     ELSE
       CALL rebuild_fist_neighbor_lists(atomic_kind_set,part,pnode,box,&
                                        internal_data%rlist_cut,&
                                        internal_data%nonbonded,&
                                        internal_data%globenv)
!MK       CALL rebuild_verlet_lists ( atomic_kind_set, part, box, internal_data, pnode )
     ENDIF

     iw = internal_data % globenv % scr
     ionode = internal_data % globenv % ionode
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
