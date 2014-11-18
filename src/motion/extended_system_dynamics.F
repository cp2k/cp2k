!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      CJM 20-Feb-2001: Now npt_ifo is allocated to zero when not used
!>      CJM 11-apr-2001: adding routines to thermostat ao_type
!>      CJM 02-Aug-2003: renamed
!> \author CJM
! *****************************************************************************
MODULE extended_system_dynamics

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: lnhc_parameters_type,&
                                             map_info_type,&
                                             npt_info_type
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE thermostat_utils,                ONLY: ke_region_baro,&
                                             ke_region_particles,&
                                             ke_region_shells,&
                                             vel_rescale_baro,&
                                             vel_rescale_particles,&
                                             vel_rescale_shells
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  LOGICAL, PARAMETER :: debug_this_module=.FALSE.
  PUBLIC :: shell_scale_comv,&
            lnhc_particles,&
            lnhc_barostat,&
            lnhc_shells

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'extended_system_dynamics'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param nhc ...
!> \param npt ...
!> \param group ...
!> \param error ...
!> \date 13-DEC-2000
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE lnhc_barostat ( nhc, npt, group, error )

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    INTEGER, INTENT(IN)                      :: group
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lnhc_barostat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineN,handle)
    map_info => nhc%map_info

    ! Compute the kinetic energy of the barostat
    CALL ke_region_baro(map_info, npt, group, error)

    ! Calculate forces on the Nose-Hoover Thermostat and apply chains
    CALL do_nhc(nhc, map_info, error)

    ! Now scale the particle velocities
    CALL vel_rescale_baro(map_info, npt, error)

    CALL timestop(handle)
  END SUBROUTINE lnhc_barostat

! *****************************************************************************
!> \brief ...
!> \param nhc ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param group ...
!> \param shell_adiabatic ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \param error ...
!> \date 14-NOV-2000
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE lnhc_particles( nhc, molecule_kind_set, molecule_set,  &
       particle_set, local_molecules, group, shell_adiabatic,&
       shell_particle_set, core_particle_set, vel, shell_vel, core_vel,&
       error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: shell_adiabatic
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set( : ), &
                                                core_particle_set( : )
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lnhc_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: my_shell_adiabatic
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineN,handle)
    my_shell_adiabatic = .FALSE.
    IF(PRESENT(shell_adiabatic)) my_shell_adiabatic = shell_adiabatic
    map_info => nhc%map_info

    ! Compute the kinetic energy for the region to thermostat
    CALL ke_region_particles(map_info, particle_set, molecule_kind_set,&
         local_molecules, molecule_set, group, vel, error)

    ! Calculate forces on the Nose-Hoover Thermostat and apply chains
    CALL do_nhc(nhc, map_info, error)

    ! Now scale the particle velocities
    CALL vel_rescale_particles(map_info, molecule_kind_set, molecule_set, particle_set,&
         local_molecules, my_shell_adiabatic, shell_particle_set,core_particle_set,&
         vel, shell_vel, core_vel, error)

    CALL timestop(handle)
  END SUBROUTINE lnhc_particles

! *****************************************************************************
!> \brief ...
!> \param nhc ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param group ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \param error ...
!> \date 14-NOV-2000
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE lnhc_shells(nhc, atomic_kind_set, particle_set, local_particles, &
       group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel,&
       error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    INTEGER, INTENT(IN)                      :: group
    TYPE(particle_type), OPTIONAL, POINTER   :: shell_particle_set(:), &
                                                core_particle_set(:)
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: vel(:,:), shell_vel(:,:), &
                                                core_vel(:,:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'lnhc_shells', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(map_info_type), POINTER             :: map_info

    CALL timeset(routineN,handle)
    map_info => nhc%map_info

    ! Compute the kinetic energy of the region to thermostat
    CALL ke_region_shells(map_info, particle_set, atomic_kind_set, local_particles,&
    group, core_particle_set, shell_particle_set, core_vel, shell_vel, error)

    ! Calculate forces on the Nose-Hoover Thermostat and apply chains
    CALL do_nhc(nhc, map_info, error)

    ! Now scale the particle velocities
    CALL vel_rescale_shells(map_info, atomic_kind_set, particle_set, local_particles, &
         shell_particle_set, core_particle_set, shell_vel, core_vel, vel, error)

    CALL timestop(handle)
  END SUBROUTINE lnhc_shells

! *****************************************************************************
!> \brief ...
!> \param nhc ...
!> \param map_info ...
!> \param error ...
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! *****************************************************************************
  SUBROUTINE do_nhc(nhc, map_info, error)
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(map_info_type), POINTER             :: map_info
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_nhc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imap, n
    LOGICAL                                  :: failure

    failure = .FALSE.

    ! Force on the first bead in every thermostat chain
    DO n = 1, nhc%loc_num_nhc
       imap = nhc%map_info%map_index(n)
       IF (nhc%nvt(1,n)%nkt==0.0_dp) CYCLE
       nhc%nvt(1,n)%f = (map_info%s_kin(imap)-nhc%nvt(1,n)%nkt)/nhc%nvt(1,n)%mass
    END DO

    ! Perform multiple time stepping using Yoshida
    CALL multiple_step_yoshida(nhc,error)

  END SUBROUTINE do_nhc

! *****************************************************************************
!> \brief ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param com_vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \date 14-NOV-2000
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE shell_scale_comv(atomic_kind_set, local_particles, particle_set,&
       com_vel, shell_vel, core_vel)

    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(IN)                :: com_vel( :, : )
    REAL(KIND=dp), INTENT(INOUT)             :: shell_vel( :, : ), &
                                                core_vel( :, : )

    CHARACTER(len=*), PARAMETER :: routineN = 'shell_scale_comv', &
      routineP = moduleN//':'//routineN

    INTEGER :: iparticle, iparticle_kind, iparticle_local, nparticle_kind, &
      nparticle_local, shell_index
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: fac_massc, fac_masss, mass, &
                                                vc(3), vs(3)
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    nparticle_kind = SIZE(atomic_kind_set)

    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass,&
            shell_active=is_shell,shell=shell)
       IF(is_shell) THEN
          fac_masss = shell%mass_shell/mass
          fac_massc = shell%mass_core/mass
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             shell_index = particle_set(iparticle)%shell_index
             vs(1:3) = shell_vel( 1:3,shell_index )
             vc(1:3) = core_vel( 1:3,shell_index )
             shell_vel( 1,shell_index ) = com_vel( 1, iparticle ) + fac_massc*(vs(1)-vc(1))
             shell_vel( 2,shell_index ) = com_vel( 2, iparticle ) + fac_massc*(vs(2)-vc(2))
             shell_vel( 3,shell_index ) = com_vel( 3, iparticle ) + fac_massc*(vs(3)-vc(3))
             core_vel( 1,shell_index ) = com_vel( 1, iparticle ) + fac_masss*(vc(1)-vs(1))
             core_vel( 2,shell_index ) = com_vel( 2, iparticle ) + fac_masss*(vc(2)-vs(2))
             core_vel( 3,shell_index ) = com_vel( 3, iparticle ) + fac_masss*(vc(3)-vs(3))
          END DO
       END IF  ! is_shell
    END DO  ! iparticle_kind
  END SUBROUTINE shell_scale_comv

! *****************************************************************************
!> \brief ...
!> \param nhc ...
!> \param error ...
!> \date 14-NOV-2000
!> \par History
!>      none
! *****************************************************************************
  SUBROUTINE multiple_step_yoshida(nhc,error)

    TYPE(lnhc_parameters_type), POINTER      :: nhc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'multiple_step_yoshida', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imap, inc, inhc, iyosh, n, &
                                                nx1, nx2
    REAL(KIND=dp)                            :: scale
    TYPE(map_info_type), POINTER             :: map_info

    nx1 = SIZE ( nhc % nvt ,1 )
    nx2 = SIZE ( nhc % nvt ,2 )
    map_info => nhc%map_info
    ! perform multiple time stepping using Yoshida
    NCLOOP: DO inc = 1, nhc % nc
       YOSH: DO iyosh = 1, nhc % nyosh

          ! update velocity on the last thermostat in the chain    ! O1
          nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
               nhc % nvt(nhc % nhc_len,:) % f * 0.25_dp*nhc % dt_yosh(iyosh)*nhc%dt_fact

          ! update velocity of other thermostats on chain (from nhc_len-1 to 1)  ! O2
          DO n = 1, nhc%loc_num_nhc
             IF (nhc % nvt(1,n) % nkt == 0.0_dp) CYCLE
             DO inhc = nhc % nhc_len - 1, 1, -1
                scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh)*nhc%dt_fact)
                nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
                nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                     nhc % nvt(inhc,n) % f* 0.25_dp*nhc % dt_yosh(iyosh)*nhc%dt_fact ! shift
                nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
             END DO
          END DO

          ! the core of the operator ----- START------
          ! update nhc positions
          nhc % nvt ( :, : ) % eta = nhc % nvt ( :, : ) % eta + &
               0.5_dp*nhc % nvt ( :, : ) % v* nhc % dt_yosh(iyosh)*nhc%dt_fact

          ! now accumulate the scale factor for particle velocities
          DO n = 1, nhc%loc_num_nhc
             imap = nhc%map_info%map_index(n)
             IF (nhc % nvt(1,n) % nkt == 0.0_dp) CYCLE
             map_info%v_scale(imap) = map_info%v_scale(imap)*EXP(-0.5_dp* nhc%dt_yosh(iyosh)*nhc%dt_fact*nhc % nvt(1,n)%v)
          END DO
          ! the core of the operator ------ END ------

          ! update the force on first thermostat again (since particle velocities changed)
          DO n = 1, nhc%loc_num_nhc
             imap = nhc%map_info%map_index(n)
             IF (nhc % nvt(1,n) % nkt == 0.0_dp) CYCLE
             nhc % nvt(1,n) % f = (map_info%s_kin(imap)*map_info%v_scale(imap)*&
                  map_info%v_scale(imap)-nhc % nvt(1,n) % nkt) / nhc % nvt(1,n) % mass
          END DO

          ! update velocity of other thermostats on chain (from 1 to nhc_len-1)  ! O2
          DO inhc = 1, nhc % nhc_len - 1
             DO n = 1, nhc%loc_num_nhc
                IF (nhc % nvt(1,n) % nkt == 0.0_dp) CYCLE
                scale = EXP(-0.125_dp*nhc % nvt(inhc+1,n) % v*nhc % dt_yosh(iyosh)*nhc%dt_fact)
                nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
                nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v + &
                     nhc % nvt(inhc,n) % f*0.25_dp*nhc % dt_yosh(iyosh)*nhc%dt_fact ! shift
                nhc % nvt(inhc,n) % v = nhc % nvt(inhc,n) % v*scale ! scale
             END DO

             ! updating the forces on all the thermostats
             DO n = 1, nhc%loc_num_nhc
                IF (nhc % nvt(1,n) % nkt == 0.0_dp) CYCLE
                nhc % nvt(inhc+1,n) % f = (nhc % nvt(inhc,n) % mass*nhc % nvt(inhc,n) % v &
                     *nhc % nvt(inhc,n) % v-nhc % nvt(inhc+1,n) % nkt)/  nhc % nvt(inhc+1,n) % mass
             END DO
          END DO
          ! update velocity on last thermostat                             ! O1
          nhc % nvt(nhc % nhc_len,:) % v = nhc % nvt(nhc % nhc_len,:) % v + &
               nhc % nvt(nhc % nhc_len,:) % f * 0.25_dp*nhc % dt_yosh(iyosh)*nhc%dt_fact
       END DO YOSH
    END DO NCLOOP
  END SUBROUTINE multiple_step_yoshida

END MODULE extended_system_dynamics
