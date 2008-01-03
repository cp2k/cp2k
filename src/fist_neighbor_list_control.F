!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE fist_neighbor_list_types,        ONLY: fist_neighbor_type
  USE fist_neighbor_lists,             ONLY: build_fist_neighbor_lists
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_set,&
                                             fist_nonbond_env_type,&
                                             pos_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_max
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             tersoff_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_neighbor_list_control'

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
       cell, fist_nonbond_env, para_env, force_env_section,shell_particle_set,&
       core_particle_set, force_update, error)

    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(particle_type), POINTER             :: shell_particle_set( : ), &
                                                core_particle_set( : )
    LOGICAL, INTENT(IN), OPTIONAL            :: force_update
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'list_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: counter, handle, i, ikind, iparticle, iparticle_kind, &
      iparticle_local, ishell, jkind, last_update, nparticle, nparticle_kind, &
      nparticle_local, nshell, num_update, output_unit, stat
    LOGICAL :: build_from_scratch, failure, geo_check, shell_adiabatic, &
      shell_present, update_neighbor_lists
    LOGICAL, DIMENSION(:, :), POINTER        :: pot_type
    REAL(KIND=dp)                            :: aup, dr2, dr2_max, lup, &
                                                verlet_skin
    REAL(KIND=dp), DIMENSION(3)              :: dr, rab, s, s2r
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rlist_cutsq, rlist_lowsq
    TYPE(cell_type), POINTER                 :: cell_last_update
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fist_neighbor_type), POINTER        :: nonbonded
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(pos_type), DIMENSION(:), POINTER    :: r_last_update, &
                                                r_last_update_pbc, &
                                                rcore_last_update_pbc, &
                                                rshell_last_update_pbc

    failure = .FALSE.
    CALL timeset(routineN,"I"," ",handle )
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! *** Assigning local pointers ***
    CALL fist_nonbond_env_get(fist_nonbond_env,&
                              nonbonded=nonbonded,&
                              rlist_cutsq=rlist_cutsq,&
                              rlist_lowsq=rlist_lowsq,&
                              aup=aup,&
                              lup=lup, &
                              counter=counter,&
                              r_last_update=r_last_update,&
                              r_last_update_pbc=r_last_update_pbc,&
                              rshell_last_update_pbc=rshell_last_update_pbc,&
                              rcore_last_update_pbc=rcore_last_update_pbc,&
                              cell_last_update=cell_last_update,&
                              num_update=num_update,&
                              potparm=potparm,&
                              last_update=last_update,error=error)
    
    nparticle_kind = SIZE(atomic_kind_set)
    nparticle = SIZE(particle_set)
    nshell = 0
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present, shell_adiabatic=shell_adiabatic)
    IF(shell_present) THEN
       nshell = SIZE(shell_particle_set) 
    END IF

    ! *** Check, if the neighbor lists have to be built or updated ***
    update_neighbor_lists = .FALSE.
    CALL section_vals_val_get(force_env_section,"MM%FORCEFIELD%NEIGHBOR_LISTS%NEIGHBOR_LISTS_FROM_SCRATCH",&
         l_val=build_from_scratch,error=error)
    CALL section_vals_val_get(force_env_section,"MM%FORCEFIELD%NEIGHBOR_LISTS%GEO_CHECK",&
         l_val=geo_check,error=error)
    IF (ASSOCIATED(r_last_update)) THEN
       CALL section_vals_val_get(force_env_section,"MM%FORCEFIELD%NEIGHBOR_LISTS%VERLET_SKIN",&
            r_val=verlet_skin,error=error)
       dr2_max = 0.0_dp
       DO iparticle_kind=1,nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             s2r = r_last_update(iparticle)%r
             s=particle_set(iparticle)%r(:)
             dr(:) = s2r - s
             dr2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
             dr2_max = MAX(dr2_max,dr2)
          END DO
       END DO
       
       CALL mp_max(dr2_max,para_env%group)
       
       IF (dr2_max > 0.25_dp*verlet_skin**2.OR.build_from_scratch) THEN
          DO iparticle=1,nparticle
             r_last_update(iparticle )%r=particle_set(iparticle)%r(:)
          END DO
          update_neighbor_lists = .TRUE.
       END IF
    ELSE
       ALLOCATE (r_last_update(nparticle),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iparticle=1,nparticle
          r_last_update(iparticle )%r=particle_set(iparticle)%r(:)
       END DO
       
       update_neighbor_lists = .TRUE.
       build_from_scratch = .TRUE.
    END IF
    ! Force Update
    IF (PRESENT(force_update)) THEN
       IF (force_update) update_neighbor_lists = .TRUE.
    END IF
    
    ! Update the rlast_update_pbc
    IF (.NOT.ASSOCIATED(r_last_update_pbc)) THEN
       ALLOCATE (r_last_update_pbc(nparticle),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(shell_present .AND. .NOT.ASSOCIATED(rshell_last_update_pbc)) THEN
       ALLOCATE (rshell_last_update_pbc(nshell),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(shell_present .AND. .NOT.ASSOCIATED(rcore_last_update_pbc)) THEN
       ALLOCATE (rcore_last_update_pbc(nshell),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    ! update the neighbor lists
    IF (update_neighbor_lists) THEN
       ALLOCATE(pot_type(nparticle_kind,nparticle_kind),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO ikind = 1, nparticle_kind
          DO jkind = ikind, nparticle_kind
             pot_type(ikind, jkind) = .FALSE.
             IF (ANY(potparm%pot(ikind, jkind)%pot%type==tersoff_type)) THEN
                pot_type(ikind, jkind) = .TRUE.
             END IF
             pot_type(jkind, ikind) = pot_type(ikind, jkind)
          END DO
       END DO
       CALL build_fist_neighbor_lists(atomic_kind_set,particle_set,&
                                      local_particles,cell,rlist_cutsq,rlist_lowsq,&
                                      nonbonded,para_env,&
                                      build_from_scratch=build_from_scratch,&
                                      geo_check=geo_check,&
                                      force_env_section=force_env_section,&
                                      pot_type=pot_type,&
                                      error=error)
       
       CALL cell_release(cell_last_update,error=error)
       CALL cell_create(cell_last_update,error=error)
       CALL cell_clone(cell,cell_last_update,error)
       
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
                                 r_last_update_pbc=r_last_update_pbc,&
                                 rshell_last_update_pbc=rshell_last_update_pbc,&
                                 rcore_last_update_pbc=rcore_last_update_pbc,&
                                 nonbonded=nonbonded,&
                                 num_update=num_update,&
                                 last_update=last_update,&
                                 cell_last_update=cell_last_update,error=error)

       output_unit = cp_print_key_unit_nr(logger,force_env_section,"MM%PRINT%NEIGHBOR_LISTS",&
            extension=".mmLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,&
               FMT="(/,T2,A,/,T52,A,/,A,T31,A,T49,2(1X,F15.2),/,T2,A,/)")&
               REPEAT("*",79),"INSTANTANEOUS        AVERAGES",&
               " LIST UPDATES[steps]","= ",lup,aup,REPEAT("*",79)
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
            "MM%PRINT%NEIGHBOR_LISTS",error=error)
       DEALLOCATE(pot_type,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    
    DO iparticle=1,nparticle
       rab = r_last_update(iparticle)%r
       IF (cell%orthorhombic) THEN
          rab(1) = - cell%hmat(1,1)*cell%perd(1)*ANINT(cell_last_update%h_inv(1,1)*rab(1))
          rab(2) = - cell%hmat(2,2)*cell%perd(2)*ANINT(cell_last_update%h_inv(2,2)*rab(2))
          rab(3) = - cell%hmat(3,3)*cell%perd(3)*ANINT(cell_last_update%h_inv(3,3)*rab(3))
       ELSE
          s(1) = cell_last_update%h_inv(1,1)*rab(1) + cell_last_update%h_inv(1,2)*rab(2) +&
               cell_last_update%h_inv(1,3)*rab(3)
          s(2) = cell_last_update%h_inv(2,1)*rab(1) + cell_last_update%h_inv(2,2)*rab(2) +&
               cell_last_update%h_inv(2,3)*rab(3)
          s(3) = cell_last_update%h_inv(3,1)*rab(1) + cell_last_update%h_inv(3,2)*rab(2) +&
               cell_last_update%h_inv(3,3)*rab(3)
          s(1) = - cell%perd(1)*ANINT(s(1))
          s(2) = - cell%perd(2)*ANINT(s(2))
          s(3) = - cell%perd(3)*ANINT(s(3))
          rab(1) = + cell%hmat(1,1)*s(1) + cell%hmat(1,2)*s(2) + cell%hmat(1,3)*s(3)
          rab(2) = + cell%hmat(2,1)*s(1) + cell%hmat(2,2)*s(2) + cell%hmat(2,3)*s(3)
          rab(3) = + cell%hmat(3,1)*s(1) + cell%hmat(3,2)*s(2) + cell%hmat(3,3)*s(3)
       END IF
       r_last_update_pbc(iparticle )%r=particle_set(iparticle)%r+rab
       ishell = particle_set(iparticle)%shell_index
       IF(ishell/=0) THEN
          rshell_last_update_pbc(ishell )%r =  rab + shell_particle_set(ishell)%r(:)
          IF(shell_adiabatic) THEN
             rcore_last_update_pbc(ishell )%r = rab + core_particle_set(ishell)%r(:)
          ELSE
             rcore_last_update_pbc(ishell )%r = r_last_update_pbc(iparticle )%r(:)
          END IF
       END IF
    END DO
    
    counter = counter + 1
    CALL fist_nonbond_env_set(fist_nonbond_env,counter=counter,error=error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE list_control

!******************************************************************************

END MODULE fist_neighbor_list_control
