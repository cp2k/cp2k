!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides integrator utility routines for the integrators
!> \par History
!>      Teodoro Laino [tlaino] 02.2008 - Splitting from integrator and creation
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
! *****************************************************************************
MODULE integrator_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE barostat_types,                  ONLY: do_clv_x,&
                                             do_clv_xy,&
                                             do_clv_xyz,&
                                             do_clv_xz,&
                                             do_clv_y,&
                                             do_clv_yz,&
                                             do_clv_z
  USE cell_types,                      ONLY: cell_type
  USE constraint,                      ONLY: rattle_roll_control
  USE constraint_util,                 ONLY: pv_constraint
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE extended_system_types,           ONLY: debug_isotropic_limit,&
                                             debug_uniaxial_limit,&
                                             npt_info_type
  USE input_constants,                 ONLY: &
       isokin_ensemble, npe_f_ensemble, npe_i_ensemble, &
       nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, npt_f_ensemble, &
       npt_i_ensemble, nve_ensemble, nvt_ensemble
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matvec_3x3,&
                                             transpose_3d
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             global_constraint_type,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type,&
                                             update_particle_set
  USE shell_potential_types,           ONLY: shell_kind_type
  USE simpar_types,                    ONLY: simpar_type
  USE thermostat_types,                ONLY: set_thermostats,&
                                             thermostats_type
  USE virial_types,                    ONLY: virial_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'integrator_utils'

! *****************************************************************************
  TYPE old_variables_type
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: v
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: r
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: eps
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: veps
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: h
  END TYPE old_variables_type

! *****************************************************************************
  TYPE tmp_variables_type
    INTEGER, POINTER :: itimes
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: pos
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: vel
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: shell_pos
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: shell_vel
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: core_pos
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: core_vel
    REAL(KIND=dp) :: max_vel, max_vel_sc
    REAL(KIND=dp) :: max_dvel, max_dvel_sc
    REAL(KIND=dp) :: max_dr, max_dsc
    REAL(KIND=dp) :: arg_r(3), arg_v(3), u(3,3), e_val(3), s, ds
    REAL(KIND=dp), DIMENSION(3) :: poly_r, &
                                   poly_v, scale_r, scale_v
  END TYPE tmp_variables_type

  INTERFACE set
     MODULE PROCEDURE set_particle_set, set_vel
  END INTERFACE

  INTERFACE update_pv
     MODULE PROCEDURE update_pv_particle_set, update_pv_velocity
  END INTERFACE

  INTERFACE damp_v
     MODULE PROCEDURE damp_v_particle_set, damp_v_velocity
  END INTERFACE

  PUBLIC :: old_variables_type, tmp_variables_type, allocate_old, deallocate_old,&
       allocate_tmp, update_dealloc_tmp, damp_v, set, update_pv, get_s_ds,&
       rattle_roll_setup, damp_veps, update_veps, vv_first, vv_second,&
       variable_timestep

CONTAINS


! *****************************************************************************
!> \brief ...
!> \param old ...
!> \param particle_set ...
!> \param npt ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE allocate_old ( old, particle_set, npt)
    TYPE(old_variables_type), POINTER        :: old
    TYPE(particle_type), POINTER             :: particle_set(:)
    TYPE(npt_info_type), POINTER             :: npt(:,:)

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_old', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, jdim, natoms, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    natoms = SIZE ( particle_set)
    idim = SIZE ( npt, 1)
    jdim = SIZE ( npt, 2)
    CPASSERT(.NOT.ASSOCIATED(old))
    ALLOCATE ( old, STAT = stat )
    CPASSERT(stat==0)

    ALLOCATE ( old%v ( natoms, 3 ), STAT = stat )
    CPASSERT(stat==0)
    old%v = 0.0_dp
    ALLOCATE ( old%r ( natoms, 3 ), STAT = stat )
    CPASSERT(stat==0)
    old%r = 0.0_dp
    ALLOCATE ( old%eps ( idim, jdim ),STAT = stat )
    CPASSERT(stat==0)
    old%eps = 0.0_dp
    ALLOCATE ( old%veps ( idim, jdim ),STAT = stat )
    CPASSERT(stat==0)
    old%veps = 0.0_dp
    ALLOCATE ( old%h(3, 3),STAT=stat)
    CPASSERT(stat==0)
    old%h = 0.0_dp

  END SUBROUTINE allocate_old

! *****************************************************************************
!> \brief ...
!> \param old ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE deallocate_old (old)
    TYPE(old_variables_type), POINTER        :: old

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_old', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(old)) THEN
       IF (ASSOCIATED(old%v)) THEN
          DEALLOCATE ( old%v, STAT = stat )
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(old%r)) THEN
          DEALLOCATE ( old%r, STAT = stat )
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(old%eps)) THEN
          DEALLOCATE ( old%eps, STAT = stat )
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(old%veps)) THEN
          DEALLOCATE ( old%veps, STAT = stat )
          CPASSERT(stat==0)
       END IF
       IF (ASSOCIATED(old%h)) THEN
          DEALLOCATE ( old%h, STAT=stat)
          CPASSERT(stat==0)
       END IF
       DEALLOCATE ( old, STAT=stat)
       CPASSERT(stat==0)
    END IF

  END SUBROUTINE deallocate_old

! *****************************************************************************
!> \brief allocate temporary variables to store positions and velocities
!>        used by the velocity-verlet integrator
!> \param md_env ...
!> \param tmp ...
!> \param nparticle ...
!> \param nshell ...
!> \param shell_adiabatic ...
!> \par History
!>      none
!> \author MI (February 2008)
! *****************************************************************************
  SUBROUTINE allocate_tmp ( md_env, tmp, nparticle, nshell, shell_adiabatic)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(tmp_variables_type), POINTER        :: tmp
    INTEGER, INTENT(IN)                      :: nparticle, nshell
    LOGICAL, INTENT(IN)                      :: shell_adiabatic

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_tmp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPASSERT(.NOT.ASSOCIATED(tmp))
    ALLOCATE (tmp, STAT=stat)
    CPASSERT(stat==0)

    ! Nullify Components
    NULLIFY(tmp%pos)
    NULLIFY(tmp%vel)
    NULLIFY(tmp%shell_pos)
    NULLIFY(tmp%shell_vel)
    NULLIFY(tmp%core_pos)
    NULLIFY(tmp%core_vel)
    NULLIFY(tmp%itimes)

    !     *** Allocate work storage for positions and velocities ***
    ALLOCATE (tmp%pos(3,nparticle),STAT=stat)
    CPASSERT(stat==0)

    ALLOCATE (tmp%vel(3,nparticle),STAT=stat)
    CPASSERT(stat==0)

    tmp%pos(:,:) = 0.0_dp
    tmp%vel(:,:) = 0.0_dp

    IF(shell_adiabatic) THEN
       !     *** Allocate work storage for positions and velocities ***
       ALLOCATE (tmp%shell_pos(3,nshell),STAT=stat)
       CPASSERT(stat==0)

       ALLOCATE (tmp%core_pos(3,nshell),STAT=stat)
       CPASSERT(stat==0)

       ALLOCATE (tmp%shell_vel(3,nshell),STAT=stat)
       CPASSERT(stat==0)

       ALLOCATE (tmp%core_vel(3,nshell),STAT=stat)
       CPASSERT(stat==0)

       tmp%shell_pos(:,:) = 0.0_dp
       tmp%shell_vel(:,:) = 0.0_dp
       tmp%core_pos(:,:) = 0.0_dp
       tmp%core_vel(:,:) = 0.0_dp
   END IF

   tmp%arg_r = 0.0_dp
   tmp%arg_v = 0.0_dp
   tmp%u = 0.0_dp
   tmp%e_val = 0.0_dp
   tmp%max_vel = 0.0_dp
   tmp%max_vel_sc = 0.0_dp
   tmp%max_dvel = 0.0_dp
   tmp%max_dvel_sc = 0.0_dp
   tmp%scale_r = 1.0_dp
   tmp%scale_v = 1.0_dp
   tmp%poly_r = 1.0_dp
   tmp%poly_v = 1.0_dp

   CALL get_md_env(md_env=md_env, itimes=tmp%itimes)

  END SUBROUTINE allocate_tmp

! *****************************************************************************
!> \brief update positions and deallocate temporary variable
!> \param tmp ...
!> \param particle_set ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param para_env ...
!> \param shell_adiabatic ...
!> \param pos ...
!> \param vel ...
!> \param should_deall_vel ...
!> \par History
!>      none
!> \author MI (February 2008)
! *****************************************************************************
  SUBROUTINE update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
       core_particle_set, para_env, shell_adiabatic, pos, vel,&
       should_deall_vel)

    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                shell_particle_set, &
                                                core_particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    LOGICAL, INTENT(IN), OPTIONAL            :: pos, vel, should_deall_vel

    CHARACTER(len=*), PARAMETER :: routineN = 'update_dealloc_tmp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure, my_deall, my_pos, &
                                                my_vel

    failure = .FALSE.

    CPASSERT(ASSOCIATED(tmp))
    my_pos =.FALSE.
    my_vel =.FALSE.
    my_deall=.TRUE.
    IF(PRESENT(pos)) my_pos=pos
    IF(PRESENT(vel)) my_vel=vel
    IF(PRESENT(should_deall_vel))my_deall=should_deall_vel

    !      *** Broadcast the new particle positions ***
    IF(my_pos) THEN
       CALL update_particle_set ( particle_set, para_env % group, pos = tmp%pos)
       DEALLOCATE (tmp%pos,STAT=stat)
       CPASSERT(stat==0)

       IF(shell_adiabatic)THEN
          CALL update_particle_set ( shell_particle_set, para_env%group, pos=tmp%shell_pos)
          CALL update_particle_set ( core_particle_set, para_env%group, pos=tmp%core_pos)
          DEALLOCATE (tmp%shell_pos,STAT=stat)
          CPASSERT(stat==0)
          DEALLOCATE (tmp%core_pos,STAT=stat)
          CPASSERT(stat==0)
       END IF
    END IF

    !     *** Broadcast the new particle velocities ***
    IF(my_vel) THEN
       CALL update_particle_set ( particle_set, para_env % group, vel = tmp%vel)
       IF (shell_adiabatic) THEN
          CALL update_particle_set ( shell_particle_set, para_env%group, vel=tmp%shell_vel)
          CALL update_particle_set ( core_particle_set, para_env%group, vel=tmp%core_vel)
       END IF
       IF(my_deall)THEN
          DEALLOCATE (tmp%vel,STAT=stat)
          CPASSERT(stat==0)
          IF (shell_adiabatic) THEN
             DEALLOCATE (tmp%shell_vel,STAT=stat)
             CPASSERT(stat==0)
             DEALLOCATE (tmp%core_vel,STAT=stat)
             CPASSERT(stat==0)
          END IF
          CPASSERT(.NOT.ASSOCIATED(tmp%pos))
          CPASSERT(.NOT.ASSOCIATED(tmp%shell_pos))
          CPASSERT(.NOT.ASSOCIATED(tmp%core_pos))
          DEALLOCATE (tmp,STAT=stat)
          CPASSERT(stat==0)
       END IF
    END IF

  END SUBROUTINE  update_dealloc_tmp

! *****************************************************************************
!> \brief ...
!> \param tmp ...
!> \param nparticle_kind ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param dt ...
!> \param para_env ...
!> \param tmpv ...
!> \par History
!>      - Created [2004-07]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE get_s_ds(tmp, nparticle_kind, atomic_kind_set,  local_particles, particle_set,&
       dt, para_env, tmpv)
    TYPE(tmp_variables_type), POINTER        :: tmp
    INTEGER                                  :: nparticle_kind
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp)                            :: dt
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN), OPTIONAL            :: tmpv

    CHARACTER(len=*), PARAMETER :: routineN = 'get_s_ds', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_local
    LOGICAL                                  :: my_tmpv
    REAL(KIND=dp)                            :: a, b, K, mass, rb
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    my_tmpv=.FALSE.
    IF(PRESENT(tmpv)) my_tmpv = tmpv

    K=0.0_dp
    a=0.0_dp
    b=0.0_dp
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       IF (mass/=0.0_dp) THEN
          nparticle_local = local_particles%n_el(iparticle_kind)
          IF(my_tmpv) THEN
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                K=K+ 0.5_dp*mass *DOT_PRODUCT(tmp%vel(:,iparticle),tmp%vel(:,iparticle))
                a=a+              DOT_PRODUCT(tmp%vel(:,iparticle),particle_set(iparticle)%f(:))
                b=b+(1.0_dp/mass)*DOT_PRODUCT(particle_set(iparticle)%f(:),particle_set(iparticle)%f(:))
             END DO
          ELSE
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                K=K+ 0.5_dp*mass *DOT_PRODUCT(particle_set(iparticle)%v(:),particle_set(iparticle)%v(:))
                a=a+              DOT_PRODUCT(particle_set(iparticle)%v(:),particle_set(iparticle)%f(:))
                b=b+(1.0_dp/mass)*DOT_PRODUCT(particle_set(iparticle)%f(:),particle_set(iparticle)%f(:))
             END DO
          END IF
       END IF
    END DO
    CALL mp_sum(K,para_env % group)
    CALL mp_sum(a,para_env % group)
    CALL mp_sum(b,para_env % group)
    a=a/(2.0_dp*K)
    b=b/(2.0_dp*K)
    rb=SQRT(b)
    tmp%s =(a/b)*(COSH(dt*rb/2.0_dp)-1) +SINH(dt*rb/2.0_dp)/rb
    tmp%ds=(a/b)*(SINH(dt*rb/2.0_dp)*rb)+COSH(dt*rb/2.0_dp)

  END SUBROUTINE get_s_ds

! *****************************************************************************
!> \brief ...
!> \param old ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param cell ...
!> \param npt ...
!> \param char ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE set_particle_set ( old, atomic_kind_set, particle_set, local_particles, cell , npt, char )
    TYPE(old_variables_type), POINTER        :: old
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt
    CHARACTER(LEN=*), INTENT(IN)             :: char

    CHARACTER(len=*), PARAMETER :: routineN = 'set_particle_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, iparticle, &
                                                iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_kind, &
                                                nparticle_local

    nparticle_kind = SIZE ( atomic_kind_set )
    SELECT CASE (char)
    CASE ( 'F') ! forward assigning the old
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el ( iparticle_kind )
          DO iparticle_local = 1, nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             DO idim = 1, 3
                old % v ( iparticle, idim ) =  particle_set ( iparticle ) % v ( idim )
                old % r ( iparticle, idim ) =  particle_set ( iparticle ) % r ( idim )
             END DO
          END DO
       END DO
       old % eps ( :, : ) = npt ( :, : ) % eps
       old % veps ( :, : ) = npt ( :, : ) % v
       old % h ( :, : ) = cell%hmat ( :, : )
    CASE ( 'B') ! back assigning the original variables
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el ( iparticle_kind )
          DO iparticle_local = 1, nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             DO idim = 1, 3
                particle_set ( iparticle ) % v ( idim ) = old % v ( iparticle, idim )
                particle_set ( iparticle ) % r ( idim ) = old % r ( iparticle, idim )
             END DO
          END DO
       END DO
       npt ( :, : ) % eps = old % eps ( :, : )
       npt ( :, : ) % v = old % veps ( :, : )
       cell % hmat ( :, : )=old % h ( :, : )
    END SELECT

  END SUBROUTINE set_particle_set

! *****************************************************************************
!> \brief ...
!> \param old ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param vel ...
!> \param local_particles ...
!> \param cell ...
!> \param npt ...
!> \param char ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE set_vel ( old, atomic_kind_set, particle_set, vel, local_particles, cell , npt, char )
    TYPE(old_variables_type), POINTER        :: old
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cell_type), POINTER                 :: cell
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt
    CHARACTER(LEN=*), INTENT(IN)             :: char

    CHARACTER(len=*), PARAMETER :: routineN = 'set_vel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: idim, iparticle, &
                                                iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_kind, &
                                                nparticle_local

    nparticle_kind = SIZE ( atomic_kind_set )
    SELECT CASE (char)
    CASE ( 'F') ! forward assigning the old
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el ( iparticle_kind )
          DO iparticle_local = 1, nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             DO idim = 1, 3
                old % v ( iparticle, idim ) =  vel ( idim, iparticle )
                old % r ( iparticle, idim ) =  particle_set ( iparticle ) % r ( idim )
             END DO
          END DO
       END DO
       old % eps ( :, : ) = npt ( :, : ) % eps
       old % veps ( :, : ) = npt ( :, : ) % v
       old % h ( :, : ) = cell%hmat ( :, : )
    CASE ( 'B') ! back assigning the original variables
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el ( iparticle_kind )
          DO iparticle_local = 1, nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             DO idim = 1, 3
                vel ( idim, iparticle ) = old % v ( iparticle, idim )
                particle_set ( iparticle ) % r ( idim ) = old % r ( iparticle, idim )
             END DO
          END DO
       END DO
       npt ( :, : ) % eps = old % eps ( :, : )
       npt ( :, : ) % v = old % veps ( :, : )
       cell % hmat ( :, : )=old % h ( :, : )
    END SELECT

  END SUBROUTINE set_vel


! *****************************************************************************
!> \brief overloaded routine provides damping for particles via nph_uniaxial_damped dynamics
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param gamma1 ...
!> \param npt ...
!> \param dt ...
!> \param group ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE damp_v_particle_set (molecule_kind_set, molecule_set,&
       particle_set, local_molecules, gamma1, npt, dt, group )

    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    REAL(KIND=dp), INTENT(IN)                :: gamma1
    TYPE(npt_info_type), INTENT(IN)          :: npt
    REAL(KIND=dp), INTENT(IN)                :: dt
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: first_atom, ikind, imol, &
                                                imol_local, ipart, last_atom, &
                                                nmol_local
    REAL(KIND=dp)                            :: alpha, ikin, kin, mass, scale
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule

    kin = 0.0_dp
    DO ikind = 1, SIZE ( molecule_kind_set )
       ! Compute the total kinetic energy
       nmol_local = local_molecules % n_el ( ikind )
       DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom = first_atom, &
               last_atom = last_atom )
          DO ipart = first_atom, last_atom
             atomic_kind => particle_set ( ipart ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             kin = kin + mass * particle_set ( ipart ) % v ( 1 ) *            &
                  particle_set ( ipart ) % v ( 1 )
             kin = kin + mass * particle_set ( ipart ) % v ( 2 ) *            &
                  particle_set ( ipart ) % v ( 2 )
             kin = kin + mass * particle_set ( ipart ) % v ( 3 ) *            &
                  particle_set ( ipart ) % v ( 3 )
          END DO
       END DO
    END DO
    !
    CALL mp_sum ( kin, group )
    !
    ikin = 1.0_dp/ kin
    scale = 1.0_dp
    alpha = 2.0_dp * npt % mass * npt % v * npt % v * gamma1 * ikin
    scale = scale * SQRT ( 1.0_dp + alpha* 0.5_dp * dt )
    ! Scale
    DO ikind = 1, SIZE ( molecule_kind_set )
       nmol_local = local_molecules % n_el ( ikind )
       DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom = first_atom, &
               last_atom = last_atom )
          DO ipart = first_atom, last_atom
             particle_set ( ipart ) % v ( 1 ) =  particle_set ( ipart ) % v ( 1 ) * scale
             particle_set ( ipart ) % v ( 2 ) =  particle_set ( ipart ) % v ( 2 ) * scale
             particle_set ( ipart ) % v ( 3 ) =  particle_set ( ipart ) % v ( 3 ) * scale
          END DO
       END DO
    END DO

  END SUBROUTINE damp_v_particle_set

! *****************************************************************************
!> \brief overloaded provides damping for particles via nph_uniaxial_damped dynamics
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param vel ...
!> \param gamma1 ...
!> \param npt ...
!> \param dt ...
!> \param group ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE damp_v_velocity (molecule_kind_set, molecule_set, &
       particle_set, local_molecules, vel, gamma1, npt, dt, group )

    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: gamma1
    TYPE(npt_info_type), INTENT(IN)          :: npt
    REAL(KIND=dp), INTENT(IN)                :: dt
    INTEGER, INTENT(IN)                      :: group

    INTEGER                                  :: first_atom, ikind, imol, &
                                                imol_local, ipart, last_atom, &
                                                nmol_local
    REAL(KIND=dp)                            :: alpha, ikin, kin, mass, scale
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_type), POINTER             :: molecule

    kin = 0.0_dp
    ! Compute the total kinetic energy
    DO ikind = 1, SIZE ( molecule_kind_set )
       nmol_local = local_molecules % n_el ( ikind )
       DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom = first_atom, &
               last_atom = last_atom )
          DO ipart = first_atom, last_atom
             atomic_kind => particle_set ( ipart ) % atomic_kind
             CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
             kin = kin + mass * vel ( 1, ipart ) * vel ( 1, ipart )
             kin = kin + mass * vel ( 2, ipart ) * vel ( 2, ipart )
             kin = kin + mass * vel ( 3, ipart ) * vel ( 3, ipart )
          END DO
       END DO
    END DO
    !
    CALL mp_sum ( kin, group )
    !
    ikin = 1.0_dp/ kin
    scale = 1.0_dp
    alpha = 2.0_dp * npt % mass * npt % v * npt % v * gamma1  * ikin
    scale = scale * SQRT ( 1.0_dp + alpha* 0.5_dp * dt )
    ! Scale
    DO ikind = 1, SIZE ( molecule_kind_set )
       nmol_local = local_molecules % n_el ( ikind )
       DO imol_local = 1, nmol_local
          imol = local_molecules % list ( ikind ) % array ( imol_local )
          molecule => molecule_set ( imol )
          CALL get_molecule ( molecule, first_atom = first_atom, &
               last_atom = last_atom )
          DO ipart = first_atom, last_atom
             vel ( 1, ipart ) = vel ( 1, ipart ) * scale
             vel ( 2, ipart ) = vel ( 2, ipart ) * scale
             vel ( 3, ipart ) = vel ( 3, ipart ) * scale
          END DO
       END DO
    END DO

  END SUBROUTINE damp_v_velocity

! *****************************************************************************
!> \brief provides damping for barostat via nph_uniaxial_damped dynamics
!> \param npt ...
!> \param gamma1 ...
!> \param dt ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE damp_veps ( npt, gamma1, dt )

    TYPE(npt_info_type), INTENT(INOUT)       :: npt
    REAL(KIND=dp), INTENT(IN)                :: gamma1, dt

    REAL(KIND=dp)                            :: scale

    scale = 1.0_dp
    scale = scale * EXP ( -gamma1 * 0.25_dp * dt )
    ! Scale
    npt % v = npt % v * scale

  END SUBROUTINE damp_veps

! *****************************************************************************
!> \brief update veps using multiplier obtained from SHAKE
!> \param old ...
!> \param gci ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param local_molecules ...
!> \param vel ...
!> \param dt ...
!> \param cell ...
!> \param npt ...
!> \param simpar ...
!> \param virial ...
!> \param vector_v ...
!> \param roll_tol ...
!> \param iroll ...
!> \param infree ...
!> \param first ...
!> \param para_env ...
!> \param u ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE rattle_roll_setup ( old, gci, atomic_kind_set, particle_set, local_particles, &
       molecule_kind_set, molecule_set, local_molecules, vel, dt, cell, npt, simpar, virial,&
       vector_v, roll_tol, iroll, infree, first, para_env, u)

    TYPE(old_variables_type), POINTER        :: old
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    REAL(KIND=dp), INTENT(IN)                :: dt
    TYPE(cell_type), POINTER                 :: cell
    TYPE(npt_info_type), DIMENSION(:, :), &
      POINTER                                :: npt
    TYPE(simpar_type), INTENT(IN)            :: simpar
    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vector_v
    REAL(KIND=dp), INTENT(OUT)               :: roll_tol
    INTEGER, INTENT(INOUT)                   :: iroll
    REAL(KIND=dp), INTENT(IN)                :: infree
    LOGICAL, INTENT(INOUT)                   :: first
    TYPE(cp_para_env_type), INTENT(IN)       :: para_env
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: u( : , : )

    REAL(KIND=dp)                            :: kin
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin
    TYPE(npt_info_type), DIMENSION(3, 3)     :: npt_loc

    IF ( first ) THEN
       CALL update_pv ( gci, simpar, atomic_kind_set, vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

    ENDIF
    first = .FALSE.

    ! assigning local variable
    SELECT CASE ( simpar % ensemble )
    CASE ( npt_i_ensemble )
       npt_loc ( :, : ) % v = 0.0_dp
       npt_loc ( :, : ) % mass = 0.0_dp
       npt_loc ( 1, 1 ) % v = npt ( 1, 1 ) % v
       npt_loc ( 2, 2 ) % v = npt ( 1, 1 ) % v
       npt_loc ( 3, 3 ) % v = npt ( 1, 1 ) % v
       npt_loc ( 1, 1 ) % mass = npt ( 1, 1 ) % mass
       npt_loc ( 2, 2 ) % mass = npt ( 1, 1 ) % mass
       npt_loc ( 3, 3 ) % mass = npt ( 1, 1 ) % mass
    CASE ( npt_f_ensemble )
       npt_loc = npt
    END SELECT

    ! resetting
    CALL set ( old, atomic_kind_set, particle_set, vel, local_particles, cell , npt, 'B' )
    CALL rattle_roll_control( gci, local_molecules, molecule_set, molecule_kind_set, &
         particle_set, vel, dt, simpar, vector_v, npt_loc%v, roll_tol, iroll,&
         para_env, u, cell, local_particles)

  END SUBROUTINE rattle_roll_setup

! *****************************************************************************
!> \brief Overloaded routine to compute veps given the particles
!>      structure or a local copy of the velocity array
!> \param gci ...
!> \param simpar ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param local_particles ...
!> \param kin ...
!> \param pv_kin ...
!> \param virial ...
!> \param int_group ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE update_pv_particle_set ( gci, simpar, atomic_kind_set, particle_set,  &
       local_molecules, molecule_set, molecule_kind_set, &
       local_particles, kin, pv_kin, virial, int_group )
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(simpar_type), INTENT(IN)            :: simpar
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: kin
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: pv_kin
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: int_group

    INTEGER                                  :: i, iparticle, iparticle_kind, &
                                                iparticle_local, j, &
                                                nparticle_local
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    pv_kin = 0.0_dp
    kin = 0.0_dp
    DO iparticle_kind = 1, SIZE ( atomic_kind_set )
       atomic_kind => atomic_kind_set ( iparticle_kind )
       CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
       nparticle_local = local_particles % n_el ( iparticle_kind )
       DO iparticle_local = 1, nparticle_local
          iparticle = local_particles % list ( iparticle_kind ) % array ( iparticle_local )
          DO i = 1, 3
             DO j = 1, 3
                pv_kin ( i, j ) = pv_kin ( i, j ) + &
                     mass* particle_set ( iparticle ) % v ( i )* &
                     particle_set ( iparticle ) % v ( j )
             END DO
             kin = kin + mass * particle_set ( iparticle ) % v ( i ) * &
                  particle_set ( iparticle ) % v ( i )
          END DO
       END DO
    END DO

    CALL mp_sum( pv_kin, int_group )
    CALL mp_sum( kin, int_group )

    ! updating the constraint virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,    &
         molecule_kind_set, &
         particle_set, virial, &
         int_group )

  END SUBROUTINE update_pv_particle_set

! *****************************************************************************
!> \brief Overloaded routine to compute kinetic virials given the particles
!>      structure or a local copy of the velocity array
!> \param gci ...
!> \param simpar ...
!> \param atomic_kind_set ...
!> \param vel ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param local_particles ...
!> \param kin ...
!> \param pv_kin ...
!> \param virial ...
!> \param int_group ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE update_pv_velocity ( gci, simpar, atomic_kind_set, vel, particle_set,  &
       local_molecules, molecule_set, molecule_kind_set, &
       local_particles, kin, pv_kin, virial, int_group )
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(simpar_type), INTENT(IN)            :: simpar
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: vel( :, : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: kin
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_kin
    TYPE(virial_type), INTENT(INOUT)         :: virial
    INTEGER, INTENT(IN)                      :: int_group

    INTEGER                                  :: i, iparticle, iparticle_kind, &
                                                iparticle_local, j, &
                                                nparticle_local
    REAL(KIND=dp)                            :: mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

    pv_kin = 0.0_dp
    kin = 0._dp
    DO iparticle_kind = 1, SIZE ( atomic_kind_set )
       atomic_kind => atomic_kind_set ( iparticle_kind )
       CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
       nparticle_local = local_particles % n_el ( iparticle_kind )
       DO iparticle_local = 1, nparticle_local
          iparticle = local_particles % list ( iparticle_kind ) % array ( iparticle_local )
          DO i = 1, 3
             DO j = 1, 3
                pv_kin ( i, j ) = pv_kin ( i, j ) + &
                     mass* vel ( i, iparticle ) * vel ( j, iparticle )
             END DO
             kin = kin + mass * vel ( i, iparticle ) * vel ( i, iparticle )
          END DO
       END DO
    END DO

    CALL mp_sum( pv_kin, int_group )
    CALL mp_sum( kin, int_group )

    ! updating the constraint virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,    &
         molecule_kind_set, &
         particle_set, virial, &
         int_group )

  END SUBROUTINE update_pv_velocity

! *****************************************************************************
!> \brief Routine to compute veps
!> \param box ...
!> \param npt ...
!> \param simpar ...
!> \param pv_kin ...
!> \param kin ...
!> \param virial ...
!> \param infree ...
!> \param virial_components ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE update_veps ( box, npt, simpar, pv_kin, kin, virial, infree, virial_components )

    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(INOUT)                          :: npt
    TYPE(simpar_type), INTENT(IN)            :: simpar
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: pv_kin
    REAL(KIND=dp), INTENT(IN)                :: kin
    TYPE(virial_type), INTENT(INOUT)         :: virial
    REAL(KIND=dp), INTENT(IN)                :: infree
    INTEGER, INTENT(IN), OPTIONAL            :: virial_components

    INTEGER                                  :: ii
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fdotr, trace, v, v0, v0i, vi
    REAL(KIND=dp), DIMENSION(3, 3)           :: unit

    failure = .FALSE.
    ! initialize unit
    unit = 0.0_dp
    unit(1,1) = 1.0_dp
    unit(2,2) = 1.0_dp
    unit(3,3) = 1.0_dp

    IF ( simpar % ensemble == npt_i_ensemble .OR.&
         simpar % ensemble == npe_i_ensemble ) THEN
       ! get force on barostat
       fdotr = 0.0_dp
       DO ii = 1, 3
          fdotr = fdotr + virial % pv_virial ( ii, ii ) +  &
               virial % pv_constraint ( ii, ii )
       END DO

       npt ( :, : ) % f = ( 1.0_dp + ( 3.0_dp * infree ) ) * kin + fdotr - &
            3.0_dp * simpar % p_ext * box % deth
    ELSEIF ( simpar % ensemble == npt_f_ensemble .OR.&
             simpar % ensemble == npe_f_ensemble ) THEN
       npt ( :, : ) % f = virial % pv_virial ( :, : ) +                   &
            pv_kin ( :, : ) +  virial % pv_constraint ( :, : ) -  &
            unit ( :, : ) * simpar % p_ext * box % deth +       &
            infree * kin * unit ( :, : )
       IF (debug_isotropic_limit) THEN
          trace = npt ( 1, 1 ) % f + npt ( 2, 2 ) % f + npt ( 3, 3 ) % f
          trace = trace / 3.0_dp
          npt ( :, : ) % f = trace * unit ( :, : )
       END IF
    ELSEIF ( simpar % ensemble == nph_uniaxial_ensemble .OR. &
             simpar % ensemble == nph_uniaxial_damped_ensemble ) THEN
       v = box % deth
       vi = 1._dp/v
       v0 = simpar % v0
       v0i = 1._dp/v0
       ! orthorhombic box ONLY
       ! Chooses only the compressive solution
       IF ( v < v0 ) THEN
          npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +            &
               pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -      &
               simpar % p0 * v - simpar % v_shock * simpar % v_shock *   &
               v * v0i * ( 1._dp -v * v0i ) + infree * kin
       ELSE
          npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +            &
               pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -      &
               simpar % p0 * v + infree * kin
       ENDIF
       IF (debug_uniaxial_limit) THEN
          ! orthorhombic box ONLY
          npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +                   &
               pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -  &
               simpar % p0 * box % deth + infree * kin
       END IF
    ENDIF

    ! update barostat velocities
    npt ( :, : ) % v = npt ( :, : ) % v + &
         0.5_dp * simpar % dt * npt ( :, : ) % f / npt ( :, : ) % mass

    ! Screen the dynamics of the barostat according user request
    IF (PRESENT(virial_components)) THEN
       SELECT CASE(virial_components)
       CASE(do_clv_xyz)
          ! Do nothing..
       CASE(do_clv_x)
          npt(1,2)%v=0.0_dp
          npt(1,3)%v=0.0_dp
          npt(1,2)%f=0.0_dp
          npt(1,3)%f=0.0_dp
          npt(2,1)%v=0.0_dp
          npt(2,2)%v=0.0_dp
          npt(2,3)%v=0.0_dp
          npt(2,1)%f=0.0_dp
          npt(2,2)%f=0.0_dp
          npt(2,3)%f=0.0_dp
          npt(3,1)%v=0.0_dp
          npt(3,2)%v=0.0_dp
          npt(3,3)%v=0.0_dp
          npt(3,1)%f=0.0_dp
          npt(3,2)%f=0.0_dp
          npt(3,3)%f=0.0_dp
       CASE(do_clv_y)
          npt(1,1)%v=0.0_dp
          npt(1,2)%v=0.0_dp
          npt(1,3)%v=0.0_dp
          npt(1,1)%f=0.0_dp
          npt(1,2)%f=0.0_dp
          npt(1,3)%f=0.0_dp
          npt(2,1)%v=0.0_dp
          npt(2,3)%v=0.0_dp
          npt(2,1)%f=0.0_dp
          npt(2,3)%f=0.0_dp
          npt(3,1)%v=0.0_dp
          npt(3,2)%v=0.0_dp
          npt(3,3)%v=0.0_dp
          npt(3,1)%f=0.0_dp
          npt(3,2)%f=0.0_dp
          npt(3,3)%f=0.0_dp
       CASE(do_clv_z)
          npt(1,1)%v=0.0_dp
          npt(1,2)%v=0.0_dp
          npt(1,3)%v=0.0_dp
          npt(1,1)%f=0.0_dp
          npt(1,2)%f=0.0_dp
          npt(1,3)%f=0.0_dp
          npt(2,1)%v=0.0_dp
          npt(2,2)%v=0.0_dp
          npt(2,3)%v=0.0_dp
          npt(2,1)%f=0.0_dp
          npt(2,2)%f=0.0_dp
          npt(2,3)%f=0.0_dp
          npt(3,1)%v=0.0_dp
          npt(3,2)%v=0.0_dp
          npt(3,1)%f=0.0_dp
          npt(3,2)%f=0.0_dp
       CASE(do_clv_xy)
          npt(1,3)%v=0.0_dp
          npt(1,3)%f=0.0_dp
          npt(2,3)%v=0.0_dp
          npt(2,3)%f=0.0_dp
          npt(3,1)%v=0.0_dp
          npt(3,2)%v=0.0_dp
          npt(3,3)%v=0.0_dp
          npt(3,1)%f=0.0_dp
          npt(3,2)%f=0.0_dp
          npt(3,3)%f=0.0_dp
       CASE(do_clv_xz)
          npt(1,2)%v=0.0_dp
          npt(1,2)%f=0.0_dp
          npt(2,1)%v=0.0_dp
          npt(2,2)%v=0.0_dp
          npt(2,3)%v=0.0_dp
          npt(2,1)%f=0.0_dp
          npt(2,2)%f=0.0_dp
          npt(2,3)%f=0.0_dp
          npt(3,2)%v=0.0_dp
          npt(3,2)%f=0.0_dp
       CASE(do_clv_yz)
          npt(1,1)%v=0.0_dp
          npt(1,2)%v=0.0_dp
          npt(1,3)%v=0.0_dp
          npt(1,1)%f=0.0_dp
          npt(1,2)%f=0.0_dp
          npt(1,3)%f=0.0_dp
          npt(2,1)%v=0.0_dp
          npt(2,1)%f=0.0_dp
          npt(3,1)%v=0.0_dp
          npt(3,1)%f=0.0_dp
       END SELECT
    END IF

  END SUBROUTINE update_veps

! *****************************************************************************
!> \brief First half of the velocity-verlet algorithm : update velocity by half
!>        step and positions by full step
!> \param tmp ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param core_particle_set ...
!> \param shell_particle_set ...
!> \param nparticle_kind ...
!> \param shell_adiabatic ...
!> \param dt ...
!> \param u ...
!> \par History
!>      none
!> \author MI (February 2008)
! *****************************************************************************
  SUBROUTINE vv_first ( tmp, atomic_kind_set, local_particles, particle_set, &
       core_particle_set, shell_particle_set, nparticle_kind,&
       shell_adiabatic, dt, u)

    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                core_particle_set, &
                                                shell_particle_set
    INTEGER, INTENT(IN)                      :: nparticle_kind
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp)                            :: dt
    REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: u

    CHARACTER(len=*), PARAMETER :: routineN = 'vv_first', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_local, shell_index
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: dm, dmc, dms, dsc(3), &
                                                dvsc(3), fac_massc, &
                                                fac_masss, mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(atomic_kind,shell)
    tmp%max_vel = 0.0_dp
    tmp%max_vel_sc = 0.0_dp
    tmp%max_dvel = 0.0_dp
    tmp%max_dvel_sc = 0.0_dp
    tmp%max_dr = 0.0_dp
    tmp%max_dsc = 0.0_dp
    dsc = 0.0_dp
    dvsc = 0.0_dp

    !     *** Velocity Verlet (first part) ***
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass,shell_active=is_shell)
       IF (mass/=0.0_dp) THEN
          dm = 0.5_dp*dt/mass
          IF(is_shell .AND. shell_adiabatic) THEN
             CALL get_atomic_kind(atomic_kind=atomic_kind, shell=shell)
             dms = 0.5_dp*dt/shell%mass_shell
             dmc = 0.5_dp*dt/shell%mass_core
             fac_masss = shell%mass_shell/mass
             fac_massc = shell%mass_core/mass
             nparticle_local = local_particles%n_el(iparticle_kind)
             IF(PRESENT(u)) THEN
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   shell_index = particle_set(iparticle)%shell_index
                   ! Transform positions and velocities and forces of the shell
                   CALL transform_first(shell_particle_set,tmp%shell_pos,tmp%shell_vel,&
                        shell_index,u,dms,dt,tmp%poly_v,tmp%poly_r,tmp%scale_v,tmp%scale_r)

                   ! Transform positions and velocities and forces of the core
                   CALL transform_first(core_particle_set,tmp%core_pos,tmp%core_vel,&
                        shell_index,u,dmc,dt,tmp%poly_v,tmp%poly_r,tmp%scale_v,tmp%scale_r)

                   ! Derive velocities and positions of the COM
                   tmp%vel (:,iparticle) = fac_masss*tmp%shell_vel (:,shell_index)+&
                        fac_massc*tmp%core_vel (:,shell_index)
                   tmp%pos (:,iparticle) = fac_masss*tmp%shell_pos (:,shell_index)+&
                        fac_massc*tmp%core_pos (:,shell_index)

                   tmp%max_vel = MAX(tmp%max_vel,ABS(tmp%vel(1,iparticle)),&
                        ABS(tmp%vel (2,iparticle)),ABS(tmp%vel (3,iparticle)))
                   tmp%max_vel_sc = MAX(tmp%max_vel_sc,&
                        ABS(tmp%shell_vel(1,shell_index)-tmp%core_vel(1,shell_index)),&
                        ABS(tmp%shell_vel(2,shell_index)-tmp%core_vel(2,shell_index)),&
                        ABS(tmp%shell_vel(3,shell_index)-tmp%core_vel(3,shell_index)))
                   tmp%max_dr = MAX(tmp%max_dr,&
                        ABS(particle_set(iparticle)%r(1)-tmp%pos (1,iparticle)),&
                        ABS(particle_set(iparticle)%r(2)-tmp%pos (2,iparticle)),&
                        ABS(particle_set(iparticle)%r(3)-tmp%pos (3,iparticle)))
                   tmp%max_dvel = MAX(tmp%max_dvel,&
                        ABS(particle_set(iparticle)%v(1)-tmp%vel (1,iparticle)),&
                        ABS(particle_set(iparticle)%v(2)-tmp%vel (2,iparticle)),&
                        ABS(particle_set(iparticle)%v(3)-tmp%vel (3,iparticle)))

                   dsc(:) =  tmp%shell_pos(:,shell_index)-tmp%core_pos(:,shell_index)- &
                        shell_particle_set(shell_index)%r(:)+core_particle_set(shell_index)%r(:)
                   tmp%max_dsc = MAX(tmp%max_dsc,ABS(dsc(1)),ABS(dsc(2)),ABS(dsc(3)))
                   dvsc(:) =  tmp%shell_vel(:,shell_index)-tmp%core_vel(:,shell_index)- &
                        shell_particle_set(shell_index)%v(:)+core_particle_set(shell_index)%v(:)
                   tmp%max_dvel_sc = MAX(tmp%max_dvel_sc,ABS(dvsc(1)),ABS(dvsc(2)),ABS(dvsc(3)))
                END DO  ! iparticle_local
             ELSE
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   shell_index = particle_set(iparticle)%shell_index
                   tmp%shell_vel (:,shell_index) = &
                        shell_particle_set(shell_index)%v(:) * tmp%scale_v(:) * tmp%scale_v(:) +&
                        tmp%scale_v(:) * tmp%poly_v(:) * dms * shell_particle_set(shell_index)%f(:)
                   tmp%shell_pos (:,shell_index) = &
                        shell_particle_set(shell_index)%r(:) * tmp%scale_r(:) * tmp%scale_r(:)  +&
                        tmp%scale_r(:) * tmp%poly_r(:) *dt * tmp%shell_vel ( :, shell_index )
                   tmp%core_vel (:,shell_index) = &
                        core_particle_set(shell_index)%v(:) * tmp%scale_v(:) * tmp%scale_v(:) +&
                        tmp%scale_v(:) * tmp%poly_v(:) *dmc * core_particle_set(shell_index)%f(:)
                   tmp%core_pos (:,shell_index) = &
                        core_particle_set(shell_index)%r(:) * tmp%scale_r(:) * tmp%scale_r(:) +&
                        tmp%scale_r(:) * tmp%poly_r(:) *dt * tmp%core_vel ( :, shell_index )

                   tmp%vel (:,iparticle) = fac_masss*tmp%shell_vel (:,shell_index)+&
                        fac_massc*tmp%core_vel (:,shell_index)
                   tmp%pos (:,iparticle) = fac_masss*tmp%shell_pos (:,shell_index)+&
                        fac_massc*tmp%core_pos (:,shell_index)
                   tmp%max_vel = MAX(tmp%max_vel,&
                        ABS(tmp%vel(1,iparticle)),ABS(tmp%vel (2,iparticle)),ABS(tmp%vel (3,iparticle)))
                   tmp%max_vel_sc = MAX(tmp%max_vel_sc,&
                        ABS(tmp%shell_vel(1,shell_index)-tmp%core_vel(1,shell_index)),&
                        ABS(tmp%shell_vel(2,shell_index)-tmp%core_vel(2,shell_index)),&
                        ABS(tmp%shell_vel(3,shell_index)-tmp%core_vel(3,shell_index)))
                   tmp%max_dr = MAX(tmp%max_dr,&
                        ABS(particle_set(iparticle)%r(1)-tmp%pos (1,iparticle)),&
                        ABS(particle_set(iparticle)%r(2)-tmp%pos (2,iparticle)),&
                        ABS(particle_set(iparticle)%r(3)-tmp%pos (3,iparticle)))
                   tmp%max_dvel = MAX(tmp%max_dvel,&
                        ABS(particle_set(iparticle)%v(1)-tmp%vel (1,iparticle)),&
                        ABS(particle_set(iparticle)%v(2)-tmp%vel (2,iparticle)),&
                        ABS(particle_set(iparticle)%v(3)-tmp%vel (3,iparticle)))
                   dsc(:) = tmp%shell_pos (:,shell_index)-tmp%core_pos (:,shell_index)-&
                        shell_particle_set(shell_index)%r(:)+core_particle_set(shell_index)%r(:)
                   tmp%max_dsc = MAX(tmp%max_dsc,ABS(dsc(1)),ABS(dsc(2)),ABS(dsc(3)))
                   dvsc(:) =  tmp%shell_vel(:,shell_index)-tmp%core_vel(:,shell_index)- &
                        shell_particle_set(shell_index)%v(:)+core_particle_set(shell_index)%v(:)
                   tmp%max_dvel_sc = MAX(tmp%max_dvel_sc,ABS(dvsc(1)),ABS(dvsc(2)),ABS(dvsc(3)))
                END DO  ! iparticle_local
             END IF
          ELSE
             nparticle_local = local_particles%n_el(iparticle_kind)
             IF(PRESENT(u)) THEN
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   ! Transform positions and velocities and forces
                   CALL transform_first(particle_set,tmp%pos,tmp%vel,&
                        iparticle,u,dm,dt,tmp%poly_v,tmp%poly_r,tmp%scale_v,tmp%scale_r)

                   tmp%max_vel = MAX(tmp%max_vel,&
                        ABS(tmp%vel (1,iparticle)),ABS(tmp%vel (2,iparticle)),ABS(tmp%vel (3,iparticle)))
                   tmp%max_dr = MAX(tmp%max_dr,ABS(particle_set(iparticle)%r(1)-tmp%pos (1,iparticle)),&
                        ABS(particle_set(iparticle)%r(2)-tmp%pos (2,iparticle)),&
                        ABS(particle_set(iparticle)%r(3)-tmp%pos (3,iparticle)))
                   tmp%max_dvel = MAX(tmp%max_dvel,&
                        ABS(particle_set(iparticle)%v(1)-tmp%vel (1,iparticle)),&
                        ABS(particle_set(iparticle)%v(2)-tmp%vel (2,iparticle)),&
                        ABS(particle_set(iparticle)%v(3)-tmp%vel (3,iparticle)))
                END DO  ! iparticle_local
             ELSE
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   tmp%vel (1,iparticle) = &
                        particle_set(iparticle)%v(1) * tmp%scale_v(1) * tmp%scale_v(1)+&
                        tmp%scale_v(1) * tmp%poly_v(1) * dm * particle_set(iparticle)%f(1)
                   tmp%vel (2,iparticle) = &
                        particle_set(iparticle)%v(2) * tmp%scale_v(2) * tmp%scale_v(2)+&
                        tmp%scale_v(2) * tmp%poly_v(2) * dm * particle_set(iparticle)%f(2)
                   tmp%vel (3,iparticle) = &
                        particle_set(iparticle)%v(3) * tmp%scale_v(3) * tmp%scale_v(3)+&
                        tmp%scale_v(3) * tmp%poly_v(3) * dm * particle_set(iparticle)%f(3)
                   tmp%pos (1,iparticle) = &
                        particle_set(iparticle)%r(1) * tmp%scale_r(1) * tmp%scale_r(1)+&
                        tmp%scale_r(1) * tmp%poly_r(1) * dt * tmp%vel ( 1, iparticle )
                   tmp%pos (2,iparticle) = &
                        particle_set(iparticle)%r(2) * tmp%scale_r(2) * tmp%scale_r(2)+&
                        tmp%scale_r(2) * tmp%poly_r(2) * dt * tmp%vel ( 2, iparticle )
                   tmp%pos (3,iparticle) = &
                        particle_set(iparticle)%r(3) * tmp%scale_r(3) * tmp%scale_r(3)+&
                        tmp%scale_r(3) * tmp%poly_r(3) * dt * tmp%vel ( 3, iparticle )

                   tmp%max_vel = MAX(tmp%max_vel,&
                        ABS(tmp%vel (1,iparticle)),ABS(tmp%vel (2,iparticle)),ABS(tmp%vel (3,iparticle)))
                   tmp%max_dr = MAX(tmp%max_dr,&
                        ABS(particle_set(iparticle)%r(1)-tmp%pos (1,iparticle)),&
                        ABS(particle_set(iparticle)%r(2)-tmp%pos (2,iparticle)),&
                        ABS(particle_set(iparticle)%r(3)-tmp%pos (3,iparticle)))
                   tmp%max_dvel = MAX(tmp%max_dvel,&
                        ABS(particle_set(iparticle)%v(1)-tmp%vel (1,iparticle)),&
                        ABS(particle_set(iparticle)%v(2)-tmp%vel (2,iparticle)),&
                        ABS(particle_set(iparticle)%v(3)-tmp%vel (3,iparticle)))
                END DO
             END IF
          END IF
       END IF
    END DO
  END SUBROUTINE vv_first

! *****************************************************************************
!> \brief Second half of the velocity-verlet algorithm : update velocity by half
!>        step  using the new forces
!> \param tmp ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param core_particle_set ...
!> \param shell_particle_set ...
!> \param nparticle_kind ...
!> \param shell_adiabatic ...
!> \param dt ...
!> \param u ...
!> \par History
!>      none
!> \author MI (February 2008)
! *****************************************************************************
  SUBROUTINE vv_second ( tmp, atomic_kind_set, local_particles, particle_set,&
       core_particle_set, shell_particle_set, nparticle_kind,&
       shell_adiabatic, dt, u)

    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                core_particle_set, &
                                                shell_particle_set
    INTEGER, INTENT(IN)                      :: nparticle_kind
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    REAL(KIND=dp)                            :: dt
    REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL :: u

    CHARACTER(len=*), PARAMETER :: routineN = 'vv_second', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iparticle, iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_local, shell_index
    LOGICAL                                  :: is_shell
    REAL(KIND=dp)                            :: dm, dmc, dms, fac_massc, &
                                                fac_masss, mass
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(shell_kind_type), POINTER           :: shell

    NULLIFY(atomic_kind,shell)

    !     *** Velocity Verlet (second part) ***
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass,&
            shell_active=is_shell)
       IF (mass/=0.0_dp) THEN
          dm = 0.5_dp*dt/mass
          IF(is_shell .AND. shell_adiabatic) THEN
             CALL get_atomic_kind(atomic_kind=atomic_kind, shell=shell)
             dms = 0.5_dp*dt/shell%mass_shell
             dmc = 0.5_dp*dt/shell%mass_core
             fac_masss = shell%mass_shell/mass
             fac_massc = shell%mass_core/mass
             nparticle_local = local_particles%n_el(iparticle_kind)
             IF(PRESENT(u)) THEN
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   shell_index = particle_set(iparticle)%shell_index
                   ! Transform velocities and forces shell
                   CALL transform_second(shell_particle_set,tmp%shell_vel,shell_index,&
                        u,dms,tmp%poly_v,tmp%scale_v)

                   ! Transform velocities and forces core
                   CALL transform_second(core_particle_set,tmp%core_vel,shell_index,&
                        u,dmc,tmp%poly_v,tmp%scale_v)

                   ! Derive velocties of the COM
                   tmp%vel(1,iparticle) = fac_masss*tmp%shell_vel(1,shell_index) + &
                        fac_massc*tmp%core_vel(1,shell_index)
                   tmp%vel(2,iparticle) = fac_masss*tmp%shell_vel(2,shell_index)+ &
                        fac_massc*tmp%core_vel(2,shell_index)
                   tmp%vel(3,iparticle) = fac_masss*tmp%shell_vel(3,shell_index) + &
                        fac_massc*tmp%core_vel(3,shell_index)
                END DO  ! iparticle_local
             ELSE
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   shell_index = particle_set(iparticle)%shell_index
                   tmp%shell_vel (1,shell_index) = &
                        tmp%shell_vel (1,shell_index) * tmp%scale_v(1) * tmp%scale_v(1) + &
                        tmp%scale_v(1) * tmp%poly_v(1) * dms*shell_particle_set(shell_index)%f(1)
                   tmp%shell_vel (2,shell_index) =&
                        tmp%shell_vel (2,shell_index) * tmp%scale_v(2) * tmp%scale_v(2) + &
                        tmp%scale_v(2) * tmp%poly_v(2) * dms*shell_particle_set(shell_index)%f(2)
                   tmp%shell_vel (3,shell_index) =&
                        tmp%shell_vel (3,shell_index) * tmp%scale_v(3) * tmp%scale_v(3) + &
                        tmp%scale_v(3) * tmp%poly_v(3) * dms*shell_particle_set(shell_index)%f(3)
                   tmp%core_vel (1,shell_index) =&
                        tmp%core_vel (1,shell_index) * tmp%scale_v(1) * tmp%scale_v(1) + &
                        tmp%scale_v(1) * tmp%poly_v(1) * dmc*core_particle_set(shell_index)%f(1)
                   tmp%core_vel (2,shell_index) =&
                        tmp%core_vel (2,shell_index) * tmp%scale_v(2) * tmp%scale_v(2)+ &
                        tmp%scale_v(2) * tmp%poly_v(2) * dmc*core_particle_set(shell_index)%f(2)
                   tmp%core_vel (3,shell_index) =&
                        tmp%core_vel (3,shell_index) * tmp%scale_v(3) * tmp%scale_v(3)+ &
                        tmp%scale_v(3) * tmp%poly_v(3) * dmc*core_particle_set(shell_index)%f(3)

                   tmp%vel ( 1, iparticle) = fac_masss*tmp%shell_vel (1,shell_index)+ &
                        fac_massc*tmp%core_vel (1,shell_index)
                   tmp%vel ( 2, iparticle) = fac_masss*tmp%shell_vel (2,shell_index)+ &
                        fac_massc*tmp%core_vel (2,shell_index)
                   tmp%vel ( 3, iparticle) = fac_masss*tmp%shell_vel (3,shell_index)+ &
                        fac_massc*tmp%core_vel (3,shell_index)
                END DO  ! iparticle_local
             END IF
          ELSE
             nparticle_local = local_particles%n_el(iparticle_kind)
             IF(PRESENT(u)) THEN
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   CALL transform_second(particle_set,tmp%vel,iparticle,&
                        u,dm,tmp%poly_v,tmp%scale_v)
                END DO  ! iparticle_local
             ELSE
                DO iparticle_local=1,nparticle_local
                   iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                   tmp%vel ( 1, iparticle) = &
                        tmp%vel ( 1, iparticle) * tmp%scale_v(1) * tmp%scale_v(1)  + &
                        tmp%scale_v(1) * tmp%poly_v(1) * dm*particle_set(iparticle)%f(1)
                   tmp%vel ( 2, iparticle) = &
                        tmp%vel ( 2, iparticle) * tmp%scale_v(2) * tmp%scale_v(2)  + &
                        tmp%scale_v(2) * tmp%poly_v(2) * dm*particle_set(iparticle)%f(2)
                   tmp%vel ( 3, iparticle) = &
                        tmp%vel ( 3, iparticle) * tmp%scale_v(3) * tmp%scale_v(3)  + &
                        tmp%scale_v(3) * tmp%poly_v(3) * dm*particle_set(iparticle)%f(3)
                END DO
             END IF
          END IF
       END IF
    END DO

  END SUBROUTINE vv_second

! *****************************************************************************
!> \brief Transform position and velocities for the first half of the
!>        Velocity-Verlet integration in the npt_f ensemble
!> \param particle_set ...
!> \param pos ...
!> \param vel ...
!> \param index ...
!> \param u ...
!> \param dm ...
!> \param dt ...
!> \param poly_v ...
!> \param poly_r ...
!> \param scale_v ...
!> \param scale_r ...
!> \par History
!>      none
!> \author MI (February 2008)
! *****************************************************************************
  SUBROUTINE transform_first ( particle_set,pos,vel,index,u,dm,dt,poly_v, &
       poly_r, scale_v, scale_r )

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pos, vel
    INTEGER, INTENT(IN)                      :: index
    REAL(KIND=dp), INTENT(IN)                :: u(3,3), dm, dt, poly_v(3), &
                                                poly_r(3), scale_v(3), &
                                                scale_r(3)

    CHARACTER(len=*), PARAMETER :: routineN = 'transform_first', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(3)              :: uf, ur, uv

    CALL matvec_3x3 (ur, TRANSPOSE_3D ( u ), particle_set ( index ) % r ( : ) )
    CALL matvec_3x3 (uv, TRANSPOSE_3D ( u ), particle_set ( index ) % v ( : ) )
    CALL matvec_3x3 (uf, TRANSPOSE_3D ( u ), particle_set ( index ) % f ( : ) )
    !
    uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dm
    uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dm
    uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dm

    ur(1) = ur(1)*scale_r(1)*scale_r(1) + uv(1)*scale_r(1)*poly_r(1)*dt
    ur(2) = ur(2)*scale_r(2)*scale_r(2) + uv(2)*scale_r(2)*poly_r(2)*dt
    ur(3) = ur(3)*scale_r(3)*scale_r(3) + uv(3)*scale_r(3)*poly_r(3)*dt
    !
    CALL MATVEC_3x3 (pos ( :, index ), u, ur )
    CALL MATVEC_3x3 (vel ( :, index ), u, uv )

  END SUBROUTINE transform_first

! *****************************************************************************
!> \brief Transform position and velocities for the second half of the
!>        Velocity-Verlet integration in the npt_f ensemble
!> \param particle_set ...
!> \param vel ...
!> \param index ...
!> \param u ...
!> \param dm ...
!> \param poly_v ...
!> \param scale_v ...
!> \par History
!>      none
!> \author MI (February 2008)
! *****************************************************************************
  SUBROUTINE transform_second ( particle_set,vel,index,u,dm,poly_v,scale_v)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vel
    INTEGER, INTENT(IN)                      :: index
    REAL(KIND=dp), INTENT(IN)                :: u(3,3), dm, poly_v(3), &
                                                scale_v(3)

    CHARACTER(len=*), PARAMETER :: routineN = 'transform_second', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp), DIMENSION(3)              :: uf, uv

    CALL matvec_3x3 (uv, TRANSPOSE_3D ( u ), vel(:, index ) )
    CALL matvec_3x3 (uf, TRANSPOSE_3D ( u ), particle_set ( index ) % f ( : ) )
    !
    uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dm
    uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dm
    uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dm

    CALL MATVEC_3x3 (vel ( :, index ), u, uv )

  END SUBROUTINE transform_second

! *****************************************************************************
!> \brief  Compute the timestep rescaling factor
!> \param md_env ...
!> \param tmp ...
!> \param dt ...
!> \param simpar ...
!> \param para_env ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param core_particle_set ...
!> \param shell_particle_set ...
!> \param nparticle_kind ...
!> \param shell_adiabatic ...
!> \param npt ...
!> \par History
!>      none
!> \author MI (October 2008)
! *****************************************************************************
  SUBROUTINE  variable_timestep(md_env,tmp,dt,simpar,para_env,atomic_kind_set, local_particles,&
       particle_set, core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic,npt)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(tmp_variables_type), POINTER        :: tmp
    REAL(KIND=dp), INTENT(INOUT)             :: dt
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                core_particle_set, &
                                                shell_particle_set
    INTEGER, INTENT(IN)                      :: nparticle_kind
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt(:,:)

    CHARACTER(len=*), PARAMETER :: routineN = 'variable_timestep', &
      routineP = moduleN//':'//routineN

    INTEGER, SAVE                            :: itime = 0
    REAL(KIND=dp)                            :: dt_fact, dt_fact_f, &
                                                dt_fact_old, dt_fact_sc, &
                                                dt_fact_v, dt_old
    REAL(KIND=dp), POINTER                   :: time
    TYPE(thermostats_type), POINTER          :: thermostats

    dt_fact    = 1.0_dp
    dt_fact_sc = 1.0_dp
    dt_fact_f = 1.0_dp
    dt_fact_v = 1.0_dp
    dt_old = dt
    dt_fact_old = simpar%dt_fact
    simpar%dt_fact = 1.0_dp
    NULLIFY(thermostats)

    itime = itime + 1
    CALL mp_max(tmp%max_dr,para_env%group)
    IF(tmp%max_dr > simpar%dr_tol) THEN
       CALL mp_max(tmp%max_dvel,para_env%group)
       dt_fact = SQRT(simpar%dr_tol*dt/tmp%max_dvel)/dt
    END IF

    IF(shell_adiabatic) THEN
       CALL mp_max(tmp%max_dsc,para_env%group)
       IF(tmp%max_dsc > simpar%dsc_tol) THEN
          CALL mp_max(tmp%max_dvel_sc,para_env%group)
          dt_fact_sc = SQRT(simpar%dsc_tol*dt/tmp%max_dvel_sc)/dt
       END IF
    END IF

    dt_fact_f = MIN(dt_fact, dt_fact_sc,1.0_dp)
    IF(dt_fact_f < 1.0_dp) THEN
       IF(dt_fact_f < 0.1_dp) dt_fact_f = 0.1_dp
       ! repeat the first vv half-step with the rescaled time step
       dt = dt*dt_fact_f
       simpar%dt_fact = dt_fact_f
       CALL rescaled_vv_first(tmp,dt,simpar,atomic_kind_set, local_particles,&
            particle_set, core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic,npt)
    END IF

    dt_fact = 1.0_dp
    dt_fact_sc = 1.0_dp
    CALL mp_max(tmp%max_dr,para_env%group)
    IF(tmp%max_dr  > simpar%dr_tol) THEN
       CALL mp_max(tmp%max_vel,para_env%group)
       dt_fact = simpar%dr_tol/tmp%max_vel/dt
    END IF

    IF(shell_adiabatic) THEN
       CALL mp_max(tmp%max_dsc,para_env%group)
       IF(tmp%max_dsc > simpar%dsc_tol) THEN
          CALL mp_max(tmp%max_vel_sc,para_env%group)
          dt_fact_sc = simpar%dsc_tol/tmp%max_vel_sc/dt
       END IF
    END IF
    dt_fact_v = MIN(dt_fact, dt_fact_sc,1.0_dp)

    IF(dt_fact_v < 1.0_dp) THEN
       IF(dt_fact_v < 0.1_dp) dt_fact_v = 0.1_dp
       ! repeat the first vv half-step with the rescaled time step
       dt = dt*dt_fact_v
       simpar%dt_fact = dt_fact_f*dt_fact_v
       CALL rescaled_vv_first(tmp,dt,simpar,atomic_kind_set, local_particles,&
            particle_set, core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic,npt)
    END IF

    simpar%dt_fact = dt_fact_f*dt_fact_v
    IF(simpar%dt_fact  < 1.0_dp) THEN
      CALL get_md_env(md_env, t=time, thermostats=thermostats)
      time = time - dt_old + dt_old*simpar%dt_fact
      IF(ASSOCIATED(thermostats)) THEN
         CALL set_thermostats(thermostats, dt_fact=simpar%dt_fact)
      END IF
    END IF

  END SUBROUTINE variable_timestep

! *****************************************************************************
!> \brief Repeat the first step of velocity verlet witht the rescaled time step
!> \param tmp ...
!> \param dt ...
!> \param simpar ...
!> \param atomic_kind_set ...
!> \param local_particles ...
!> \param particle_set ...
!> \param core_particle_set ...
!> \param shell_particle_set ...
!> \param nparticle_kind ...
!> \param shell_adiabatic ...
!> \param npt ...
!> \par History
!>      none
!> \author MI (October 2008)
!>     I soliti ignoti
! *****************************************************************************
  SUBROUTINE rescaled_vv_first(tmp,dt,simpar,atomic_kind_set, local_particles,&
           particle_set, core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic,npt)

    TYPE(tmp_variables_type), POINTER        :: tmp
    REAL(KIND=dp), INTENT(IN)                :: dt
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                core_particle_set, &
                                                shell_particle_set
    INTEGER, INTENT(IN)                      :: nparticle_kind
    LOGICAL, INTENT(IN)                      :: shell_adiabatic
    TYPE(npt_info_type), OPTIONAL, POINTER   :: npt(:,:)

    CHARACTER(len=*), PARAMETER :: routineN = 'rescaled_vv_first', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: e2 = 1.0_dp/6.0_dp, &
                                                e4 = e2/20.0_dp, &
                                                e6 = e4/42.0_dp, &
                                                e8 = e6/72.0_dp

    REAL(KIND=dp)                            :: arg_r(3), arg_v(3), e_val(3), &
                                                infree, trvg, u(3,3)

    infree = 1.0_dp / REAL ( simpar % nfree, dp )
    arg_r = tmp%arg_r
    arg_v = tmp%arg_v
    u     = tmp%u
    e_val = tmp%e_val

    SELECT CASE(simpar%ensemble)
    CASE (nve_ensemble, nvt_ensemble)
      CALL vv_first(tmp,atomic_kind_set, local_particles, particle_set, &
           core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt)
    CASE (isokin_ensemble)
      tmp%scale_v(1:3) = SQRT(1.0_dp/tmp%ds)
      tmp%poly_v(1:3) = 2.0_dp*tmp%s/SQRT(tmp%ds)/dt
      CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
           core_particle_set, shell_particle_set, nparticle_kind,&
           shell_adiabatic, dt)

    CASE (npt_i_ensemble, npe_i_ensemble)
       arg_r = arg_r * simpar%dt_fact*simpar%dt_fact
       tmp%poly_r(1:3) = 1.0_dp + e2*arg_r(1) + e4*arg_r(1)*arg_r(1) + e6*arg_r(1)**3 + e8*arg_r(1)**4
       arg_v = arg_v * simpar%dt_fact*simpar%dt_fact
       tmp%poly_v(1:3) = 1.0_dp + e2*arg_v(1) + e4*arg_v(1)*arg_v(1) + e6*arg_v(1)**3 + e8*arg_v(1)**4
       tmp%scale_r(1:3) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
       tmp%scale_v(1:3) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
                     (1.0_dp + 3.0_dp * infree ) )
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
           core_particle_set, shell_particle_set, nparticle_kind,&
           shell_adiabatic, dt)

    CASE (npt_f_ensemble,npe_f_ensemble)
       trvg = npt ( 1, 1 ) % v  + npt ( 2, 2 ) % v + npt ( 3, 3 ) % v
       arg_r(:) = arg_r(:) * simpar%dt_fact*simpar%dt_fact
       tmp%poly_r = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
       tmp%scale_r ( : ) = EXP ( 0.5_dp * dt * e_val ( : ) )
       arg_v(:) = arg_v(:) * simpar%dt_fact*simpar%dt_fact
       tmp%poly_v = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
       tmp%scale_v ( : ) = EXP ( -0.25_dp * dt * ( &
             e_val ( : ) + trvg * infree ) )

       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
           core_particle_set, shell_particle_set, nparticle_kind,&
           shell_adiabatic, dt, u)

    CASE (nph_uniaxial_ensemble,nph_uniaxial_damped_ensemble)
        arg_r = arg_r * simpar%dt_fact*simpar%dt_fact
        tmp%poly_r(1) = 1._dp + e2*arg_r(1) + e4*arg_r(1)*arg_r(1) + e6*arg_r(1)**3 + e8*arg_r(1)**4
        arg_v(2) = arg_v(2) * simpar%dt_fact*simpar%dt_fact
        arg_v(1) = arg_v(1) * simpar%dt_fact*simpar%dt_fact
        tmp%poly_v(1) = 1._dp + e2*arg_v(1) + e4*arg_v(1)*arg_v(1) + e6*arg_v(1)**3 + e8*arg_v(1)**4
        tmp%poly_v(2) = 1._dp + e2*arg_v(2) + e4*arg_v(2)*arg_v(2) + e6*arg_v(2)**3 + e8*arg_v(2)**4
        tmp%poly_v(3) = 1._dp + e2*arg_v(2) + e4*arg_v(2)*arg_v(2) + e6*arg_v(2)**3 + e8*arg_v(2)**4
        tmp%scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
        tmp%scale_v(1) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
             (1._dp + infree ) )
        tmp%scale_v(2) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * infree )
        tmp%scale_v(3) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * infree )
        CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt)

    END SELECT

  END SUBROUTINE rescaled_vv_first

END MODULE integrator_utils
