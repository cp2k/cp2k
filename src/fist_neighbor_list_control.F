!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005 CP2K developers group                           !
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
!!     07.02.2005: using real coordinates for r_last_update; cleaned (MK)
!!
!!   SOURCE
!******************************************************************************

MODULE fist_neighbor_list_control

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE fist_neighbor_list_types,        ONLY: neighbor_list_set_p_type
  USE fist_neighbor_lists,             ONLY: build_fist_neighbor_lists
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_set,&
                                             fist_nonbond_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE message_passing,                 ONLY: mp_max
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "fist_neighbor_list_control"

  PUBLIC :: list_control

!***
!******************************************************************************

CONTAINS

!******************************************************************************

! to decide whether the neighbor list is to be updated or not
! based on a displacement criterion;
! if any particle has moved by 0.5*skin_length from the previous
! list_update, then the list routine is called.

SUBROUTINE list_control ( atomic_kind_set, particle_set, local_particles, &
                          cell, fist_nonbond_env, globenv )

    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "list_control"

    INTEGER :: counter, handle, i, iparticle, iparticle_kind, &
      iparticle_local, istat, iw, last_update, nparticle, nparticle_kind, &
      nparticle_local, num_update
    LOGICAL                                  :: ionode, list_update_flag
    REAL(KIND=dp)                            :: aup, dr2, dr2_max, lup, &
                                                verlet_skin
    REAL(KIND=dp), DIMENSION(3)              :: dr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: r_last_update, rlist_cutsq
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nonbonded

!------------------------------------------------------------------------------

  CALL timeset (routine_name,"I"," ",handle )

! Assigning local pointers

  CALL fist_nonbond_env_get(fist_nonbond_env,&
                            nonbonded=nonbonded,&
                            rlist_cutsq=rlist_cutsq,&
                            aup=aup,&
                            lup=lup, &
                            counter=counter,&
                            r_last_update=r_last_update,&
                            num_update=num_update,&
                            verlet_skin=verlet_skin,&
                            last_update=last_update)

  nparticle_kind = SIZE(atomic_kind_set)
  nparticle = SIZE(particle_set)

  IF (.NOT.ASSOCIATED(r_last_update)) THEN
    ALLOCATE (r_last_update(3,nparticle),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,"r_last_update",&
                       3*nparticle*dp_size)
    END IF
    r_last_update = 1.0E+10_dp
  END IF

! check to see if neighbor list has to be updated

  dr2_max = 0.0_dp

  DO iparticle_kind=1,nparticle_kind
    nparticle_local = local_particles%n_el(iparticle_kind)
    DO iparticle_local=1,nparticle_local
      iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
      dr(:) = r_last_update(:,iparticle) - particle_set(iparticle)%r(:)
      dr2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
      dr2_max = MAX(dr2_max,dr2)
    END DO
  END DO

#if defined(__parallel)
  CALL mp_max(dr2_max,globenv%group)
#endif

  IF (dr2_max > 0.25_dp*verlet_skin**2) THEN
    list_update_flag = .TRUE.
    DO iparticle=1,nparticle
      r_last_update(:,iparticle ) = particle_set(iparticle)%r(:)
    END DO
  ELSE
    list_update_flag = .FALSE.
  END IF

  IF (list_update_flag) THEN

       CALL build_fist_neighbor_lists(atomic_kind_set,particle_set,&
                                      local_particles,cell,rlist_cutsq,&
                                      nonbonded,globenv)

       iw = globenv % scr
       ionode = globenv % ionode

       IF ( counter > 0 ) THEN
        num_update = num_update + 1
        lup = counter + 1 - last_update
        last_update = counter + 1
        aup = aup + &
            (lup-aup)/REAL ( num_update,KIND=dp)
       ELSE
          num_update = 0
          lup = 0
          last_update = 1
          aup = 0.0_dp
       END IF

       CALL fist_nonbond_env_set ( fist_nonbond_env,  &
                                   lup = lup,         &
                                   aup = aup,         &
                                   r_last_update = r_last_update, &
                                   nonbonded = nonbonded, &
                                   num_update=num_update,  &
                                   last_update=last_update )

       IF ( ( globenv % print_level .GE. 2 ) .AND. ( ionode ) ) THEN
          WRITE ( iw, "( )" )
          WRITE ( iw, "( A,A )" ) " **************************************", &
               "*****************************************"
          WRITE ( iw, "( T52,A )" ) "INSTANTANEOUS        AVERAGES"
          WRITE ( iw, "( A,A,T31,A,T49,2(1X,F15.2) )" ) &
               " LIST UPDATES", "[steps]", "= ", lup, aup
          WRITE ( iw, "( A,A )" ) " **************************************", &
               "*****************************************"
          WRITE ( iw, "( )" )
       ENDIF
  END IF

  counter = counter + 1
  CALL fist_nonbond_env_set(fist_nonbond_env,counter=counter)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE list_control

!******************************************************************************

END MODULE fist_neighbor_list_control
