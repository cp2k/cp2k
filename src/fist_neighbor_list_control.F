!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2k developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_neighbor_list_control [1.0] *
!!
!!   NAME
!!     fist_neighbor_list_control
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

MODULE fist_neighbor_list_control
  USE atomic_kind_types, ONLY : atomic_kind_type
  USE distribution_1d_types, ONLY: distribution_1d_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl, default_string_length
  USE fist_environment_types, ONLY : fist_environment_type
  USE fist_neighbor_lists, ONLY : build_fist_neighbor_lists
  USE fist_nonbond_env_types, ONLY : fist_nonbond_env_type, &
                                     fist_nonbond_env_get, &
                                     fist_nonbond_env_set
  USE message_passing, ONLY : mp_max
  USE particle_types, ONLY : particle_type
  USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type
  USE simulation_cell, ONLY : cell_type, real_to_scaled, scaled_to_real
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: list_control

!***
!******************************************************************************

CONTAINS

!******************************************************************************
!
! to decide whether the neighbor list is to be updated or not
! based on a displacement criterion;
! if any particle has moved by 0.5*skin_length from the previous
! list_update, then the list routine is called.
!
SUBROUTINE list_control ( atomic_kind_set, particle_set, local_particles, &
                          cell, fist_nonbond_env, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles
  TYPE ( cell_type ), POINTER :: cell
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv
  TYPE ( fist_nonbond_env_type ), POINTER :: fist_nonbond_env

! Locals
  REAL ( dbl ), DIMENSION (3) :: disp_vec
  REAL ( dbl ) :: displace, max_displace, verlet_skin
  LOGICAL :: list_update_flag, ionode
  INTEGER :: counter, num_update, last_update, lup, aup
  INTEGER :: i, nnodes, node, isos, iw, handle, nparticle_local
  INTEGER :: nparticle_kind, iparticle_kind, iparticle, iparticle_local
  CHARACTER ( LEN = default_string_length) :: string
  REAL ( dbl ), POINTER, DIMENSION ( :, : ) :: r_last_update, &
                                               rlist_cutsq
  TYPE ( neighbor_list_set_p_type ),  &
                POINTER, DIMENSION ( : ) :: nonbonded
  

!------------------------------------------------------------------------------

  CALL timeset ( 'LIST_CONTROL','I',' ',handle )

! Assigning local pointers
  CALL fist_nonbond_env_get ( fist_nonbond_env,  nonbonded = nonbonded, &
                              rlist_cutsq = rlist_cutsq, aup = aup, lup = lup, &
                              counter = counter, r_last_update = r_last_update, &
                              num_update = num_update, verlet_skin = verlet_skin ) 

  nparticle_kind = SIZE ( atomic_kind_set )
  nnodes = 0
  DO iparticle_kind = 1, nparticle_kind
     nnodes = nnodes + local_particles%n_el(iparticle_kind) 
  ENDDO 

  IF ( .NOT. ASSOCIATED (  r_last_update ) ) THEN
     ALLOCATE ( r_last_update ( 3, nnodes ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ( 'fist_neighbor_list_control', &
          'r_last_update', 3 * nnodes )
      r_last_update = 1.0E+10_dbl
  END IF

! check to see if neighbor list has to be updated
  max_displace = -1.0E+15_dbl
  node = 0
  DO iparticle_kind = 1, nparticle_kind
     nparticle_local = local_particles%n_el(iparticle_kind)
     DO iparticle_local=1,nparticle_local
        node = node + 1
        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
        disp_vec ( : ) = scaled_to_real ( r_last_update ( :, node ), cell ) -  &
                         particle_set ( iparticle ) % r ( : )
        max_displace = MAX ( displace, max_displace )
     END DO
  END DO

#if defined(__parallel)
  CALL mp_max ( max_displace,  globenv % group )
#endif

  IF ( max_displace >= 0.25_dbl * verlet_skin ** 2 ) THEN
     list_update_flag = .TRUE.
     node = 0
     DO iparticle_kind = 1, nparticle_kind
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
           node = node + 1
           iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
           r_last_update ( :, node ) =  &
             real_to_scaled ( particle_set ( iparticle ) % r ( : ), cell )
        END DO
     END DO
  ELSE
     list_update_flag = .FALSE.
  END IF

  IF ( list_update_flag ) THEN

       CALL build_fist_neighbor_lists ( atomic_kind_set, particle_set, local_particles, &
                                        cell, rlist_cutsq, nonbonded, globenv )

       iw = globenv % scr
       ionode = globenv % ionode

       IF ( counter > 0 ) THEN
        num_update = num_update + 1
        lup = counter + 1 - last_update
        last_update = counter + 1
        aup = aup + &
            (lup-aup)/float(num_update)
       ELSE
          num_update = 0
          lup = 0
          last_update = 1
          aup = 0.0_dbl
       END IF

       CALL fist_nonbond_env_set ( fist_nonbond_env, lup = lup, aup = aup, &
                                  r_last_update = r_last_update, &
                                   num_update=num_update, last_update=last_update )

       IF ( ionode ) THEN
          WRITE ( iw, '( )' )
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( A,A,T31,A,T49,2(1X,F15.2) )' ) &
               ' LIST UPDATES', '[steps]', '= ', counter, aup
          WRITE ( iw, '( A,A )' ) ' **************************************', &
               '*****************************************'
          WRITE ( iw, '( )' )
       ENDIF
  END IF

  counter = counter + 1
  CALL fist_nonbond_env_set ( fist_nonbond_env, counter = counter )

  CALL timestop ( 0._dbl, handle )

END SUBROUTINE list_control

!******************************************************************************

END MODULE fist_neighbor_list_control

!******************************************************************************
