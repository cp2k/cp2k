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
!!     CJM,MK
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
  USE cell_types,                      ONLY: cell_type,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE input_section_types, ONLY: section_vals_type, section_vals_get_subs_vals,&
       section_vals_get, section_vals_val_get, section_vals_release
#include "cp_common_uses.h"

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
! if any particle has moved by 0.5*verlet_skin from the previous
! list update, then the list routine is called.

SUBROUTINE list_control ( atomic_kind_set, particle_set, local_particles, &
                          cell, fist_nonbond_env, globenv, force_env_section, error)

    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(global_environment_type), &
      pointer                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    TYPE(section_vals_type), POINTER :: force_env_section

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "list_control"

    INTEGER :: counter, handle, i, iparticle, iparticle_kind, &
      iparticle_local, istat, last_update, nparticle, nparticle_kind, &
      nparticle_local, num_update
    LOGICAL                                  :: build_from_scratch,&
                                                update_neighbor_lists
    REAL(KIND=dp)                            :: aup, dr2, dr2_max, lup, &
                                                verlet_skin
    REAL(KIND=dp), DIMENSION(3)              :: dr
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: r_last_update, rlist_cutsq
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nonbonded

!------------------------------------------------------------------------------

  CALL timeset(routine_name,"I"," ",handle )

! *** Assigning local pointers ***

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

! *** Check, if the neighbor lists have to be built or updated ***

  update_neighbor_lists = .FALSE.
  build_from_scratch = .FALSE.

  IF (ASSOCIATED(r_last_update)) THEN

    dr2_max = 0.0_dp

    DO iparticle_kind=1,nparticle_kind
      nparticle_local = local_particles%n_el(iparticle_kind)
      DO iparticle_local=1,nparticle_local
        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
        dr(:) = scaled_to_real(r_last_update(:,iparticle),cell) -&
                particle_set(iparticle)%r(:)
        dr2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
        dr2_max = MAX(dr2_max,dr2)
      END DO
    END DO

#if defined(__parallel)
    CALL mp_max(dr2_max,globenv%group)
#endif

    IF (dr2_max > 0.25_dp*verlet_skin**2) THEN

      DO iparticle=1,nparticle
        r_last_update(:,iparticle ) = real_to_scaled(particle_set(iparticle)%r(:),cell)
      END DO

      update_neighbor_lists = .TRUE.

    END IF

  ELSE

    ALLOCATE (r_last_update(3,nparticle),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,"r_last_update",&
                       3*nparticle*dp_size)
    END IF

    DO iparticle=1,nparticle
      r_last_update(:,iparticle ) = real_to_scaled(particle_set(iparticle)%r(:),cell)
    END DO

    update_neighbor_lists = .TRUE.
    build_from_scratch = .TRUE.

  END IF

  IF (update_neighbor_lists) THEN

    CALL build_fist_neighbor_lists(atomic_kind_set,particle_set,&
                                   local_particles,cell,rlist_cutsq,&
                                   nonbonded,globenv,&
                                   build_from_scratch=build_from_scratch,&
                                   force_env_section=force_env_section,error=error)

    IF ( counter > 0 ) THEN
      num_update = num_update + 1
      lup = counter + 1 - last_update
      last_update = counter + 1
      aup = aup + (lup - aup)/REAL(num_update,KIND=dp)
    ELSE
      num_update = 0
      lup = 0
      last_update = 1
      aup = 0.0_dp
    END IF

    CALL fist_nonbond_env_set(fist_nonbond_env,&
                              lup=lup,&
                              aup=aup,&
                              r_last_update=r_last_update,&
                              nonbonded=nonbonded,&
                              num_update=num_update,&
                              last_update=last_update)

    IF ((globenv%print_level >= 2).AND.(globenv%ionode)) THEN
      WRITE (UNIT=globenv%scr,&
             FMT="(/,T2,A,/,T52,A,/,A,T31,A,T49,2(1X,F15.2),/,T2,A,/)")&
        REPEAT("*",79),"INSTANTANEOUS        AVERAGES",&
        " LIST UPDATES[steps]","= ",lup,aup,REPEAT("*",79)
    END IF

  END IF

  counter = counter + 1
  CALL fist_nonbond_env_set(fist_nonbond_env,counter=counter)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE list_control

!******************************************************************************

END MODULE fist_neighbor_list_control
