!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
MODULE al_system_dynamics
 
  USE al_system_types,                 ONLY: al_system_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: map_info_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_types,                  ONLY: particle_type
  USE thermostat_utils,                ONLY: ke_region_particles,&
                                             vel_rescale_particles
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PARAMETER :: debug_this_module=.FALSE.
  PUBLIC :: al_particles

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'al_system_dynamics'

CONTAINS

! *****************************************************************************
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_particles( al, molecule_kind_set, molecule_set,  &
       particle_set, local_molecules, local_particles, group, vel, error)

    TYPE(al_system_type), POINTER            :: al
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    INTEGER, INTENT(IN)                      :: group
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'al_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: my_shell_adiabatic
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineN,handle)
    my_shell_adiabatic = .FALSE.
    map_info => al%map_info

   IF (debug_this_module) &
      CALL dump_vel(molecule_kind_set, molecule_set, local_molecules, particle_set, vel, "INIT")

    IF (al%tau_nh <= 0.0_dp) THEN
       CALL al_OU_step(0.5_dp, al, map_info, molecule_kind_set, molecule_set, &
         particle_set, local_molecules, local_particles, vel, error)
      IF (debug_this_module) &
         CALL dump_vel(molecule_kind_set, molecule_set, local_molecules, particle_set, vel, "post OU")
    ELSE
       ! quarter step of Langevin using Ornstein-Uhlenbeck
       CALL al_OU_step(0.25_dp, al, map_info, molecule_kind_set, molecule_set, &
         particle_set, local_molecules, local_particles, vel, error)
      IF (debug_this_module) &
         CALL dump_vel(molecule_kind_set, molecule_set, local_molecules, particle_set, vel, "post 1st OU")

       ! Compute the kinetic energy for the region to thermostat for the (T dependent chi step)
       CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
            local_molecules, molecule_set, group, vel=vel, error=error)
       ! quarter step of chi, and set vel drag factors for a half step
       CALL al_NH_quarter_step(al, map_info, set_half_step_vel_factors=.TRUE., error=error)

       ! Now scale the particle velocities for a NH half step
       CALL vel_rescale_particles(map_info, molecule_kind_set, molecule_set, particle_set, &
            local_molecules, my_shell_adiabatic, vel=vel, error=error)
       ! Recompute the kinetic energy for the region to thermostat (for the T dependent chi step)
       CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
            local_molecules, molecule_set, group, vel=vel, error=error)
      IF (debug_this_module) &
         CALL dump_vel(molecule_kind_set, molecule_set, local_molecules, particle_set, vel, "post rescale_vel")

       ! quarter step of chi
       CALL al_NH_quarter_step(al, map_info, set_half_step_vel_factors=.FALSE., error=error)

       ! quarter step of Langevin using Ornstein-Uhlenbeck
       CALL al_OU_step(0.25_dp, al, map_info, molecule_kind_set, molecule_set, &
         particle_set, local_molecules, local_particles, vel, error)
      IF (debug_this_module) &
         CALL dump_vel(molecule_kind_set, molecule_set, local_molecules, particle_set, vel, "post 2nd OU")
   ENDIF

    ! Recompute the final kinetic energy for the region to thermostat
    CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel=vel, error=error)

    CALL timestop(handle)
  END SUBROUTINE al_particles

   SUBROUTINE dump_vel(molecule_kind_set, molecule_set, local_molecules, particle_set, vel, label)
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(dp), OPTIONAL                       :: vel(:,:)
    CHARACTER(len=*)                         :: label

    INTEGER                                  :: first_atom, ikind, imol, &
                                                imol_local, ipart, last_atom, &
                                                nmol_local
    TYPE(molecule_type), POINTER             :: molecule

       DO ikind=1, SIZE(molecule_kind_set)
         nmol_local = local_molecules % n_el (ikind)
         DO imol_local = 1, nmol_local
            imol = local_molecules % list(ikind)  % array (imol_local)
            molecule => molecule_set(imol)
            CALL get_molecule (molecule, first_atom=first_atom, last_atom=last_atom)
            DO ipart = first_atom, last_atom
               IF (PRESENT(vel)) THEN
                  WRITE(unit=*,fmt='("VEL ",A20," IPART ",I6," V ",3F20.10)') TRIM(label), ipart, vel(:,ipart)
               ELSE
                  WRITE(unit=*,fmt='("VEL ",A20," IPART ",I6," V ",3F20.10)') TRIM(label), ipart, particle_set(ipart)%v(:)
               ENDIF
            END DO
         END DO
      END DO
   END SUBROUTINE dump_vel

  SUBROUTINE al_OU_step(step, al, map_info, molecule_kind_set, molecule_set, &
    particle_set, local_molecules, local_particles, vel, error)
    REAL(dp), INTENT(in)                     :: step
    TYPE(al_system_type), POINTER            :: al
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'al_OU_step', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, i, ii, ikind, imap, imol, imol_local, ipart, &
      iparticle_kind, iparticle_local, last_atom, nmol_local, nparticle, &
      nparticle_kind, nparticle_local, stat
    LOGICAL                                  :: check, failure, present_vel
    REAL(KIND=dp)                            :: mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: w(:,:)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(rng_stream_type), POINTER           :: rng_stream

    present_vel = PRESENT(vel)

    ![NB] not a big deal, but could this be done once at init time?
    DO i = 1, al%loc_num_al
       imap = map_info%map_index(i)
       ! drag on velocities
       IF (al%tau_langevin > 0.0_dp) THEN
          map_info%v_scale(imap) = EXP(-step*al%dt/al%tau_langevin)
          map_info%s_kin(imap) = SQRT((al%nvt(i)%nkt/al%nvt(i)%degrees_of_freedom)*(1.0_dp-map_info%v_scale(imap)**2))
       ELSE
          map_info%v_scale(imap) = 1.0_dp
          map_info%s_kin(imap) = 0.0_dp
       ENDIF
       ! magnitude of random force, not including 1/sqrt(mass) part
    END DO

    nparticle = SIZE(particle_set)
    nparticle_kind = SIZE(local_particles%n_el)
    ALLOCATE (w(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    w(:,:) = 0.0_dp
    check=(nparticle_kind <= SIZE(local_particles%n_el) .AND. nparticle_kind <= SIZE(local_particles%list))
    CPPrecondition(check, cp_failure_level, routineP, error, failure)
    check=ASSOCIATED(local_particles%local_particle_set)
    CPPrecondition(check, cp_failure_level, routineP, error, failure)
    DO iparticle_kind=1,nparticle_kind
      nparticle_local = local_particles%n_el(iparticle_kind)
      check=(nparticle_local <= SIZE(local_particles%list(iparticle_kind)%array))
      CPPrecondition( check, cp_failure_level, routineP, error, failure)
      DO iparticle_local = 1, nparticle_local
          ipart = local_particles%list(iparticle_kind)%array(iparticle_local)
          rng_stream => local_particles%local_particle_set(iparticle_kind)%rng(iparticle_local)%stream
          w(1,ipart) = next_random_number(rng_stream,variance=1.0_dp,error=error)
          w(2,ipart) = next_random_number(rng_stream,variance=1.0_dp,error=error)
          w(3,ipart) = next_random_number(rng_stream,variance=1.0_dp,error=error)
       END DO
    END DO

    ii = 0
    DO ikind=1, SIZE(molecule_kind_set)
      nmol_local = local_molecules % n_el (ikind)
      DO imol_local = 1, nmol_local
         imol = local_molecules % list(ikind)  % array (imol_local)
         molecule => molecule_set(imol)
         CALL get_molecule (molecule, first_atom=first_atom, last_atom=last_atom)
         DO ipart = first_atom, last_atom
            ii = ii + 1
            atomic_kind => particle_set(ipart) % atomic_kind
            CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
            IF (present_vel) THEN
               vel(:,ipart) = vel(:,ipart) * map_info%v_scale(ii) + &
                              map_info%s_kin(ii)/SQRT(mass)*w(:,ipart)
            ELSE
               particle_set(ipart)%v(:) = particle_set(ipart)%v(:) * map_info%v_scale(ii) + &
                                          map_info%s_kin(ii)/SQRT(mass)*w(:,ipart)
            ENDIF
         END DO
      END DO
    END DO

    DEALLOCATE (w,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE al_OU_step

! *****************************************************************************
!> \author Noam Bernstein [noamb] 02.2012
! *****************************************************************************
  SUBROUTINE al_NH_quarter_step(al, map_info, set_half_step_vel_factors, error)
    TYPE(al_system_type), POINTER            :: al
    TYPE(map_info_type), POINTER             :: map_info
    LOGICAL, INTENT(in)                      :: set_half_step_vel_factors
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'al_NH_quarter_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imap
    REAL(KIND=dp)                            :: decay, delta_K

![NB] how to deal with dt_fact?

    DO i = 1, al%loc_num_al
       IF (al%nvt(i)%mass > 0.0_dp) THEN
          imap = map_info%map_index(i)
          delta_K = 0.5_dp*(map_info%s_kin(imap)-al%nvt(i)%nkt)
          al%nvt(i)%chi = al%nvt(i)%chi + 0.5_dp*al%dt*delta_K/al%nvt(i)%mass
          IF (set_half_step_vel_factors) THEN
            decay = EXP(-0.5_dp*al%dt*al%nvt(i)%chi)
            map_info%v_scale(imap) = decay
          ENDIF
       ELSE
          al%nvt(i)%chi = 0.0_dp
          IF (set_half_step_vel_factors) THEN
             map_info%v_scale(imap) = 1.0_dp
          ENDIF
       ENDIF
    END DO

  END SUBROUTINE al_NH_quarter_step

END MODULE al_system_dynamics
