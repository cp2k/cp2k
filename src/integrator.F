!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/integrator [1.0] *
!!
!!   NAME
!!     integrator
!!
!!   FUNCTION
!!     Provides integrator routines (velocity verlet) for all the
!!     ensemble types
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (15-Mar-2001) : Pass logical for box change to force routine
!!     Harald Forbert (Apr-2001): added path integral routine nvt_pimd
!!     CJM (15-Apr-2001) : added coef integrators and energy routines
!!     Joost VandeVondele (Juli-2003): simple version of isokinetic ensemble
!!
!!   SOURCE
!******************************************************************************

MODULE integrator

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE constraint,                      ONLY: getold,&
                                             rattle_control,&
                                             shake_control, &
                                             shake_roll_control, &
                                             rattle_roll_control, &
                                             pv_constraint
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_p_type,&
                                             dyn_coeff_set_p_type,&
                                             dyn_coeff_type,&
                                             get_dyn_coeff_set
  USE eigenvalueproblems,              ONLY: diagonalise
  USE extended_system_dynamics,        ONLY: lnhc
  USE extended_system_types,           ONLY: lnhc_parameters_type, npt_info_type
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set,&
                                             force_env_p_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: dyn_coeff_mass_atom,&
                                             dyn_coeff_mass_coeff,&
                                             dyn_coeff_mass_global,&
                                             dyn_coeff_mass_kind
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type
  USE message_passing,                 ONLY: mp_comm_dup,&
                                             mp_sum
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type, get_hinv
  USE termination,                     ONLY: stop_memory, stop_program
!  USE transformations,                 ONLY: fr2fu,&
!                                             u2r

  USE util,                            ONLY: matmul_3x3, matvec_3x3, transpose_3d
  USE virial_methods,                  ONLY: virial_evaluate
  USE virial_types,                    ONLY: virial_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: module_name = "integrator"

  TYPE old_variables_type
     REAL (KIND=dp), POINTER, DIMENSION ( :, : ) :: v
     REAL (KIND=dp), POINTER, DIMENSION ( :, : ) :: r
     REAL (KIND=dp), POINTER, DIMENSION ( :, : ) :: eps
     REAL (KIND=dp), POINTER, DIMENSION ( :, : ) :: veps
     REAL (KIND=dp), POINTER, DIMENSION ( :, : ) :: h
  END TYPE old_variables_type

  TYPE ( old_variables_type ) :: old

  INTEGER :: int_group

  PUBLIC :: set_integrator, nve, nvt, isokin, npt_i, npt_f, nph_uniaxial
  PUBLIC :: nph_uniaxial_damped

  INTERFACE set
     MODULE PROCEDURE set_particle_set, set_vel
  END INTERFACE

  INTERFACE get_veps
     MODULE PROCEDURE get_veps_particle_set, get_veps_velocity
  END INTERFACE
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** integrator/set_integrator [1.0] *
!!
!!   NAME
!!     set_integrator
!!
!!   FUNCTION
!!     setup of parameters for integrator
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (29-Mar-2001) : duplicate communicator, instead of copy
!!
!!   NOTES
!!     Call this subroutine before the first call to energy or velocity_verlet
!!     or if you want to change ionode and/or output files
!!
!!   SOURCE
!******************************************************************************


SUBROUTINE set_integrator ( globenv )


    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

!------------------------------------------------------------------------------

  ! XXXXX
  ! is this one ever freed ???
  ! we should anyway use md_env % para_env % group in these routines

  CALL mp_comm_dup ( globenv % group, int_group )

END SUBROUTINE set_integrator

!!*****
!******************************************************************************
!!****** integrator/nve [1.0] *
!!
!!   NAME
!!     nve
!!
!!   FUNCTION
!!     nve integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     - the local particle lists are used instead of pnode (Sep. 2003,MK)
!!     - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE nve ( md_env )

    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "nve_new"

    INTEGER :: icoef, iparticle, iparticle_kind, iparticle_local, istat, &
      isubsys, mass_storage, nmolecule, nmolecule_kind, nparticle, &
      nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: dm, dmcoef, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      POINTER                                :: coeffs_of_kind
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: nhc_coef
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(virial_type), POINTER               :: virial

!   ---------------------------------------------------------------------------

    NULLIFY (dyn_coeff_set,coeffs_of_kind,coeffs,nhc_coef)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_coef=nhc_coef)
    dt = simpar%dt
    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("nve","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set, virial = virial )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els

!     *** Allocate work storage for positions and velocities ***

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

      IF (simpar%constraint) CALL getold(local_molecules, molecule_set,  &
                                         molecule_kind_set, particle_set)

!     *** Velocity Verlet (first part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel (:,iparticle) = particle_set(iparticle)%v(:) +&
                              dm * particle_set(iparticle)%f(:)
          pos (:,iparticle) = particle_set(iparticle)%r(:) +&
                              dt * vel ( :, iparticle )
        END DO
      END DO


      IF (simpar%constraint) CALL shake_control( local_molecules, molecule_set,   &
                                                molecule_kind_set, particle_set, &
                                                pos, vel, simpar )

!     *** Broadcast the new particle positions ***

      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"pos")

!     *** coefficients
    
      IF(ASSOCIATED(dyn_coeff_set)) THEN 
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(ASSOCIATED(nhc_coef))THEN
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                      int_group)
        END IF
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)     
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                        coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
                coeffs%pos(iparticle_local,icoef)=coeffs%pos(iparticle_local,icoef)+ &
                       dt*coeffs%vel(iparticle_local,icoef)
              END DO
            END DO
          END IF
        END DO
      END IF

!     *** Update forces ***

!      CALL force_env_set(force_env_p(1)%force_env,box_changed=.false.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle) = vel ( :, iparticle) +&
                                dm*particle_set(iparticle)%f(:)
        END DO
      END DO

      IF (simpar%constraint) CALL rattle_control( local_molecules, molecule_set, &
                                                 molecule_kind_set, particle_set, &
                                                 vel, simpar )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

!     *** coefficients

      IF(ASSOCIATED(dyn_coeff_set)) THEN
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind) 
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                         coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
              END DO
            END DO
          END IF
        END DO
        IF(ASSOCIATED(nhc_coef))THEN 
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                     int_group)
        END IF
      END IF


!     *** Update virial ***

      IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                      molecule_set,    & 
                                                      molecule_kind_set, &
                                                      particle_set, virial, &
                                                      "FULL", int_group )

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

    END DO ! next subsys or sub-system

  END SUBROUTINE nve
!!*****
!******************************************************************************
!!****** integrator/isokin [1.0] *
!!
!!   NAME
!!     isokin
!!
!!   FUNCTION
!!     simplest version of the isokinetic gaussian thermostat
!!
!!   NOTES
!!     - time reversible, and conserves the kinetic energy to machine precision
!!     - is not yet supposed to work for e.g. constraints, our the extended version
!!       of this thermostat
!!       see:
!!        - Zhang F. , JCP 106, 6102 (1997)
!!        - Minary P. et al, JCP 118, 2510 (2003)
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     - Created [2004-07]
!!
!!   SOURCE
!******************************************************************************
  SUBROUTINE isokin ( md_env )

    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "isokin"

    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, isubsys, &
      nparticle, nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: a, b, dm, ds, dt, K, mass, &
                                                rb, s
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: work
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

!   ---------------------------------------------------------------------------

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p)
                    
    dt = simpar%dt

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("isokin","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys)

    IF (SIZE(subsys)/=1) THEN
       ! presumably one just needs to get e.g. kinetic energy of the full system by looping 
       ! over the subsystems
       CALL stop_program("isokin","integrator",__LINE__,&
       "force_env with more than (or less) than one subsys not supported atm!")
    END IF

    IF (simpar%constraint) THEN
        CALL stop_program("isokin","integrator",__LINE__,&
             "Constraints not yet implemented")
    ENDIF

    do_subsys: DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nparticle = particles%n_els
      particle_set => particles%els

!     *** Allocate work storage for positions and velocities ***
      ALLOCATE (work(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "work",3*nparticle*dp_size)
      work(:,:) = 0.0_dp

      ! compute K,a,b,s,ds
      CALL get_s_ds()

!     *** 'Velocity Verlet' (first part) ***
      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = s/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          particle_set(iparticle)%v(:) = (1.0_dp/ds)*( particle_set(iparticle)%v(:) +&
                                                       dm*particle_set(iparticle)%f(:) )
          work(:,iparticle) = particle_set(iparticle)%r(:) +&
                              dt*particle_set(iparticle)%v(:)
        END DO
      END DO

!     *** Broadcast the new particle positions ***

      CALL mp_sum(work,int_group)

      DO iparticle=1,nparticle
        particle_set(iparticle)%r(:) = work(:,iparticle)
        work(:,iparticle) = 0.0_dp
      END DO

!      CALL force_env_set(force_env_p(1)%force_env,box_changed=.false.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

      ! compute K,a,b,s,ds
      CALL get_s_ds()

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = s/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          work ( :, iparticle ) = (1.0_dp/ds)*( particle_set(iparticle)%v(:) +&
                                                dm*particle_set(iparticle)%f(:))
        END DO
      END DO

!     *** Broadcast the new particle velocities ***

      CALL mp_sum(work,int_group)

      DO iparticle=1,nparticle
        particle_set(iparticle)%v(:) = work(:,iparticle)
      END DO

      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")

    END DO do_subsys
  CONTAINS
      SUBROUTINE get_s_ds()
         K=0.0_dp
         a=0.0_dp
         b=0.0_dp
         DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             K=K+ 0.5_dp*mass *DOT_PRODUCT(particle_set(iparticle)%v(:),particle_set(iparticle)%v(:))
             a=a+              DOT_PRODUCT(particle_set(iparticle)%v(:),particle_set(iparticle)%f(:))
             b=b+(1.0_dp/mass)*DOT_PRODUCT(particle_set(iparticle)%f(:),particle_set(iparticle)%f(:))
           END DO
         END DO
         CALL mp_sum(K,int_group)
         CALL mp_sum(a,int_group)
         CALL mp_sum(b,int_group)
         a=a/(2.0_dp*K)
         b=b/(2.0_dp*K)
         rb=SQRT(b)
         s =(a/b)*(COSH(dt*rb/2.0_dp)-1) +SINH(dt*rb/2.0_dp)/rb
         ds=(a/b)*(SINH(dt*rb/2.0_dp)*rb)+COSH(dt*rb/2.0_dp)
      END SUBROUTINE get_s_ds
  END SUBROUTINE isokin

!******************************************************************************
!!****** integrator/nvt [1.0] *
!!
!!   NAME
!!     nvt
!!
!!   FUNCTION
!!     nvt integrator for particle positions & momenta
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     - the local particle lists are used instead of pnode (Sep. 2003,MK)
!!     - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE nvt ( md_env )
    IMPLICIT NONE

    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "nvt"

    INTEGER :: icoef, iparticle, iparticle_kind, iparticle_local, istat, &
      isubsys, mass_storage, nmolecule, nmolecule_kind, nparticle, &
      nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: dm, dmcoef, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_p_type), DIMENSION(:), &
      POINTER                                :: coeffs_of_kind
    TYPE(dyn_coeff_set_p_type), &
      DIMENSION(:), POINTER                  :: dyn_coeff_set
    TYPE(dyn_coeff_type), POINTER            :: coeffs
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: force_env_p
    TYPE(lnhc_parameters_type), &
      DIMENSION(:), POINTER                  :: nhc_coef
    TYPE(lnhc_parameters_type), POINTER      :: nhc( : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE ( virial_type ), POINTER            :: virial

!   ---------------------------------------------------------------------------

    NULLIFY (dyn_coeff_set,coeffs_of_kind,coeffs,nhc_coef)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_part=nhc,nhc_coef=nhc_coef)
    dt = simpar%dt
    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("nvt","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       dyn_coeff_set=dyn_coeff_set, virial=virial)

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els


!     *** Allocate work storage for positions and velocities ***

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

! perform Nose'-Hoover chain dynamimcs
      CALL lnhc ( nhc ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, int_group )

      IF (simpar%constraint) CALL getold(local_molecules, molecule_set,  &
                                         molecule_kind_set, particle_set)

!     *** Velocity Verlet (first part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle ) = particle_set(iparticle)%v(:) + dm*particle_set(iparticle)%f(:)
          pos ( :, iparticle ) = particle_set(iparticle)%r(:) + dt*vel(:,iparticle)
        END DO
      END DO


      IF (simpar%constraint) CALL shake_control( local_molecules, molecule_set,   &
                                                 molecule_kind_set, particle_set, &
                                                 pos, vel, simpar )

!     *** Broadcast the new particle positions ***
      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")

!     *** coefficients

      IF(ASSOCIATED(dyn_coeff_set)) THEN
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(ASSOCIATED(nhc_coef))THEN
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                      int_group)
        END IF
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                        coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
                coeffs%pos(iparticle_local,icoef)=coeffs%pos(iparticle_local,icoef)+ &
                       dt*coeffs%vel(iparticle_local,icoef)
              END DO
            END DO
          END IF
        END DO
      END IF

!     *** Update forces ***

!      CALL force_env_set(force_env_p(1)%force_env,box_changed=.false.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle ) = vel ( :, iparticle ) + &
                                 dm*particle_set(iparticle)%f(:)
        END DO
      END DO

      IF (simpar%constraint) CALL rattle_control( local_molecules, molecule_set, &
                                                  molecule_kind_set, particle_set, &  
                                                  vel, simpar )

! perform Nose'-Hoover chain dynamimcs
      CALL lnhc ( nhc ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, vel, int_group )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
      IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                      molecule_set,    & 
                                                      molecule_kind_set, &
                                                      particle_set, virial, &
                                                      "FULL", int_group )

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

!     *** coefficients

      IF(ASSOCIATED(dyn_coeff_set)) THEN
        CALL get_dyn_coeff_set(dyn_coeff_set=dyn_coeff_set(isubsys)%dyn_coeff_set,&
                               coeffs_of_kind=coeffs_of_kind, global_mass=global_mass,&
                               mass_storage=mass_storage)
        IF(mass_storage==dyn_coeff_mass_global) dmcoef=0.5_dp*dt/global_mass
        DO iparticle_kind=1,nparticle_kind
          IF (ASSOCIATED(coeffs_of_kind(iparticle_kind)%coeffs)) THEN
            coeffs=>coeffs_of_kind(iparticle_kind)%coeffs
            IF(mass_storage==dyn_coeff_mass_kind) dmcoef=0.5_dp*dt/coeffs%masses(1,1)
            DO iparticle_local=1,coeffs%n_els
              DO icoef=1,coeffs%ncoef_atom
                SELECT CASE (mass_storage)
                CASE(dyn_coeff_mass_global,dyn_coeff_mass_kind)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+&
                         dmcoef*coeffs%forces(iparticle_local,icoef)
                CASE(dyn_coeff_mass_atom,dyn_coeff_mass_coeff)
                  coeffs%vel(iparticle_local,icoef)=coeffs%vel(iparticle_local,icoef)+0.5_dp*dt/ &
                         coeffs%masses(iparticle_local,icoef)*coeffs%forces(iparticle_local,icoef)
                END SELECT
              END DO
            END DO
          END IF
        END DO
        IF(ASSOCIATED(nhc_coef))THEN
          CALL lnhc (nhc_coef(isubsys),dyn_coeff_set(isubsys)%dyn_coeff_set,&
                     int_group)
        END IF
      END IF

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

    END DO ! next subsys or sub-system

  END SUBROUTINE nvt

!!*****
!******************************************************************************
!!****** integrator/npt_i [1.0] *
!!
!!   NAME
!!     npt_i
!!
!!   SYNOPSIS
!!     Subroutine npt_i(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine npt_i
!!
!!   FUNCTION
!!     npt_i integrator for particle positions & momenta
!!     isotropic box changes
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE npt_i ( md_env )
  
    IMPLICIT NONE
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "npt_i"

    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, &
               isubsys,  nmolecule, nmolecule_kind, nparticle, &
               nparticle_kind, nparticle_local
    INTEGER, POINTER :: itimes
    REAL(KIND=dp)                            :: dm, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part( : ), nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER          :: npt ( :, : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cell_type), POINTER                 :: cell
    TYPE ( virial_type ), POINTER            :: virial
    REAL (KIND=dp), PARAMETER :: e2 = 1.0_dp/6.0_dp, e4 = e2/20.0_dp, &
       e6 = e4/42.0_dp, e8 = e6/72.0_dp
    REAL (KIND=dp), SAVE :: eps_0
    REAL (KIND=dp) :: arg_r, poly_r, scale_r, roll_tol
    REAL (KIND=dp) :: arg_v, poly_v, scale_v, kin, infree
    REAL (KIND=dp), DIMENSION (3) :: vector_r, vector_v
    LOGICAL, PARAMETER :: box_change = .TRUE.
    LOGICAL :: first
    INTEGER :: iroll
! dbg
    INTEGER :: ii,jj


    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_part=nhc_part, nhc_baro=nhc_baro, npt = npt, &
                    itimes = itimes )
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree,KIND=dp) 

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("npt_i","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       virial = virial, cell = cell )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els


!     *** Allocate work storage for positions and velocities ***

      CALL allocate_old ( particle_set, npt )

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

!------------------------------------------------------------------------------

! initialize eps_0 the first time through
      IF ( itimes == 1 ) eps_0 = npt (1,1) % eps

! perform Nose'-Hoover chains for barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

! perform Nose'-Hoover chains for particles
      CALL lnhc ( nhc_part ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, int_group )

! setting up for ROLL: saving old variables
      roll_tol = -1.0E+10_dp
      IF (simpar%constraint) THEN
        iroll = 1
        CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
        CALL getold ( local_molecules, molecule_set,  &
                      molecule_kind_set, particle_set)
      ENDIF


      SR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! SHAKE-ROLL LOOP

          IF (simpar%constraint) THEN
             CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
          END IF
          CALL get_veps ( cell, npt, simpar, atomic_kind_set, particle_set, &
                          local_molecules, molecule_set, molecule_kind_set, &
                          local_particles, virial, infree, int_group )

!
! first half of velocity verlet
!
          arg_r = ( 0.5_dp * npt ( 1, 1 ) % v * dt ) * &
                  ( 0.5_dp * npt ( 1, 1 ) % v * dt )
          poly_r = 1.0_dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

          arg_v = ( 0.25_dp * npt ( 1, 1 ) % v * dt *  &
                  (1.0_dp + 3.0_dp * infree ) )*( 0.25_dp * npt ( 1, 1 ) % v * &
                   dt * ( 1.0_dp + 3.0_dp * infree ) )
          poly_v = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

          scale_r = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
          scale_v = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
                       (1.0_dp + 3.0_dp * infree ) )
         DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           dm = 0.5_dp*dt/mass
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             vel ( :, iparticle ) = particle_set(iparticle)%v(:) * scale_v * scale_v + &
                                    scale_v * poly_v * dm * particle_set(iparticle)%f(:)
             pos ( :, iparticle ) = particle_set(iparticle)%r(:) * scale_r * scale_r + &
                                    scale_r * poly_r * dt *vel(:,iparticle)
           END DO
         END DO

        roll_tol = 0.0_dp
        vector_r ( : ) = scale_r * poly_r
        vector_v ( : ) = scale_v * poly_v

        IF ( simpar % constraint ) CALL shake_roll_control( local_molecules, &
                                        molecule_set, molecule_kind_set, &
                                        particle_set, pos, vel, simpar,  &
                                        roll_tol, iroll, vector_r, vector_v,  &
                                        int_group )
  END DO SR

! Update eps:
      npt ( :, : ) % eps =  npt ( :, : ) % eps + dt * npt ( :, : ) % v

! Update h_mat
      cell % hmat ( :, : ) = cell % hmat ( :, : ) &
      * EXP ( npt ( 1, 1 ) % eps - eps_0 )

      eps_0 =  npt ( 1, 1 ) % eps

! Update the inverse
      CALL get_hinv ( cell )

!     *** Broadcast the new particle positions ***
      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")


!     *** Update forces ***

      CALL force_env_set(force_env_p(1)%force_env,box_changed=.TRUE.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( :, iparticle ) = vel ( :, iparticle ) * scale_v * scale_v + &
                                 scale_v * poly_v * dm * particle_set(iparticle)%f(:)
        END DO
      END DO

      roll_tol = -1.E10_dp
      IF (simpar%constraint) THEN
         first = .TRUE.
         iroll = 1
         CALL set ( atomic_kind_set, particle_set, vel, local_particles, cell, npt, 'F' )
      END IF

      RR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! RATTLE-ROLL LOOP
         roll_tol = 0.0_dp
         IF ( simpar % constraint ) CALL rattle_roll_setup ( atomic_kind_set, &
                                      particle_set, local_particles, &
                                      molecule_kind_set, molecule_set,  &
                                      local_molecules, vel,  &
                                      cell, npt, simpar, virial, vector_v,  &
                                      roll_tol, iroll, infree, first, int_group )

         CALL get_veps ( cell, npt, simpar, atomic_kind_set, vel, particle_set,&
                         local_molecules, molecule_set, molecule_kind_set, &
                         local_particles, virial, infree, int_group )

      END DO RR


! Nose-Hoover for the Particles
      CALL lnhc ( nhc_part ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, vel, int_group )

! Nose-Hoover for the Barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
      IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                      molecule_set,    & 
                                                      molecule_kind_set, &
                                                      particle_set, virial, &
                                                      "FULL", int_group )

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

! Deallocate old variables
    CALL deallocate_old ( )
  END DO

  END SUBROUTINE npt_i
!!*****
!******************************************************************************
!!****** integrator/nph_uniaxial [1.0] *
!!
!!   NAME
!!     nph_uniaxial
!!
!!   SYNOPSIS
!!     Subroutine npt_uniaxial(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine npt_i
!!
!!   FUNCTION
!!     nph_uniaxial integrator (non-Hamiltonian version) 
!!     for particle positions & momenta undergoing
!!     uniaxial stress ( in x-direction of orthorhombic cell) 
!!     due to a shock compression:
!!     Reed et. al. Physical Review Letters 90, 235503 (2003).
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE nph_uniaxial ( md_env )
  
    IMPLICIT NONE
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "npt_uniaxial"

    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, &
               isubsys,  nmolecule, nmolecule_kind, nparticle, &
               nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: dm, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE ( npt_info_type ), POINTER          :: npt ( :, : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cell_type), POINTER                 :: cell
    TYPE ( virial_type ), POINTER            :: virial
    REAL ( dp ), PARAMETER :: e2 = 1._dp/6._dp, e4 = e2/20._dp, &
       e6 = e4/42._dp, e8 = e6/72._dp
    REAL ( dp ) :: arg_r, poly_r, scale_r, roll_tol
    REAL ( dp ) :: arg_v, poly_v, scale_v, kin, infree
    REAL ( dp ), DIMENSION (3) :: vector_r, vector_v
    LOGICAL, PARAMETER :: box_change = .TRUE.
    LOGICAL :: first
    INTEGER :: iroll


    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    npt = npt )
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree, dp ) 

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("npt_i","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       virial = virial, cell = cell )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els


!     *** Allocate work storage for positions and velocities ***

      CALL allocate_old ( particle_set, npt )

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

!------------------------------------------------------------------------------

! setting up for ROLL: saving old variables
      roll_tol = -1.0E+10_dp
      IF (simpar%constraint) THEN
        iroll = 1
        CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
        CALL getold ( local_molecules, molecule_set,  &
                      molecule_kind_set, particle_set)
      ENDIF


      SR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! SHAKE-ROLL LOOP

          IF (simpar%constraint) THEN
             CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
          END IF
          CALL get_veps ( cell, npt, simpar, atomic_kind_set, particle_set, &
                          local_molecules, molecule_set, molecule_kind_set, &
                          local_particles, virial, infree, int_group )

!
! first half of velocity verlet
!
          arg_r = ( 0.5_dp * npt ( 1, 1 ) % v * dt ) * &
                  ( 0.5_dp * npt ( 1, 1 ) % v * dt )
          poly_r = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

          arg_v = ( 0.25_dp * npt ( 1, 1 ) % v * dt *  &
                  (1._dp + infree ) )*( 0.25_dp * npt ( 1, 1 ) % v * &
                   dt * ( 1._dp + infree ) )
          poly_v = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

          scale_r = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
          scale_v = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
                       (1._dp + infree ) )
         DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           dm = 0.5_dp*dt/mass
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             vel ( 1, iparticle ) = particle_set(iparticle)%v(1) * scale_v * scale_v + &
                                    scale_v * poly_v * dm * particle_set(iparticle)%f(1)
             vel ( 2, iparticle ) = particle_set(iparticle)%v(2) +  dm * particle_set(iparticle)%f(2)
             vel ( 3, iparticle ) = particle_set(iparticle)%v(3) +  dm * particle_set(iparticle)%f(3)
             pos ( 1, iparticle ) = particle_set(iparticle)%r(1) * scale_r * scale_r + &
                                    scale_r * poly_r * dt * vel ( 1, iparticle )
             pos ( 2, iparticle ) = particle_set(iparticle)%r(2) + dt*vel(2,iparticle)
             pos ( 3, iparticle ) = particle_set(iparticle)%r(3) + dt*vel(3,iparticle)
           END DO
         END DO

        roll_tol = 0._dp
        vector_r ( : ) = 0._dp
        vector_v ( : ) = 0._dp
        vector_r ( 1 ) = scale_r * poly_r
        vector_v ( 1 ) = scale_v * poly_v

        IF ( simpar % constraint ) CALL shake_roll_control( local_molecules, &
                                        molecule_set, molecule_kind_set, &
                                        particle_set, pos, vel, simpar,  &
                                        roll_tol, iroll, vector_r, vector_v,  &
                                        int_group )
  END DO SR

! Update h_mat

       cell % hmat ( 1, 1 ) = cell % hmat ( 1, 1 ) * scale_r * scale_r 

! Update the inverse
      CALL get_hinv ( cell )

!     *** Broadcast the new particle positions ***
      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")


!     *** Update forces ***

      CALL force_env_set(force_env_p(1)%force_env,box_changed=.TRUE.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( 1, iparticle ) = vel ( 1, iparticle ) * scale_v * scale_v + &
                                 scale_v * poly_v * dm * particle_set(iparticle)%f(1)
          vel ( 2, iparticle ) = vel ( 2, iparticle ) +  dm * particle_set(iparticle)%f(2)
          vel ( 3, iparticle ) = vel ( 3, iparticle ) +  dm * particle_set(iparticle)%f(3)
        END DO
      END DO

      roll_tol = -1.E10_dp
      IF (simpar%constraint) THEN
         first = .TRUE.
         iroll = 1
         CALL set ( atomic_kind_set, particle_set, vel, local_particles, cell, npt, 'F' )
      END IF

      RR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! RATTLE-ROLL LOOP
         roll_tol = 0._dp
         IF ( simpar % constraint ) CALL rattle_roll_setup ( atomic_kind_set, &
                                      particle_set, local_particles, &
                                      molecule_kind_set, molecule_set,  &
                                      local_molecules, vel,  &
                                      cell, npt, simpar, virial, vector_v,  &
                                      roll_tol, iroll, infree, first, int_group )

         CALL get_veps ( cell, npt, simpar, atomic_kind_set, vel, particle_set,&
                         local_molecules, molecule_set, molecule_kind_set, &
                         local_particles, virial, infree, int_group )

      END DO RR

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
      IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                      molecule_set,    & 
                                                      molecule_kind_set, &
                                                      particle_set, virial, &
                                                      "FULL", int_group )

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

! Deallocate old variables
    CALL deallocate_old ( )
  END DO

  END SUBROUTINE nph_uniaxial
!!*****
!******************************************************************************
!!****** integrator/nph_uniaxial_damped [1.0] *
!!
!!   NAME
!!     nph_uniaxial_damped
!!
!!   SYNOPSIS
!!     Subroutine nph_uniaxial(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine npt_i
!!
!!   FUNCTION
!!     nph_uniaxial integrator (non-Hamiltonian version) 
!!     for particle positions & momenta undergoing
!!     uniaxial stress ( in x-direction of orthorhombic cell) 
!!     due to a shock compression:
!!     Reed et. al. Physical Review Letters 90, 235503 (2003).
!!     Added damping (e.g. thermostat to barostat)
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE nph_uniaxial_damped ( md_env )
  
    IMPLICIT NONE
    TYPE(md_environment_type), POINTER       :: md_env

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "nph_uniaxial_damped"

    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, &
               isubsys,  nmolecule, nmolecule_kind, nparticle, &
               nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: dm, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE ( npt_info_type ), POINTER          :: npt ( :, : )
    TYPE(lnhc_parameters_type), POINTER      :: nhc_baro ( : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cell_type), POINTER                 :: cell
    TYPE ( virial_type ), POINTER            :: virial
    REAL ( dp ), PARAMETER :: e2 = 1._dp/6._dp, e4 = e2/20._dp, &
       e6 = e4/42._dp, e8 = e6/72._dp
    REAL ( dp ) :: arg_r, poly_r, scale_r, roll_tol
    REAL ( dp ) :: arg_v, poly_v, scale_v, kin, infree
    REAL ( dp ), DIMENSION (3) :: vector_r, vector_v
    LOGICAL, PARAMETER :: box_change = .TRUE.
    LOGICAL :: first
    INTEGER :: iroll


    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    npt = npt, nhc_baro=nhc_baro )
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree, dp ) 

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("npt_i","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       virial = virial, cell = cell )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els


!     *** Allocate work storage for positions and velocities ***

      CALL allocate_old ( particle_set, npt )

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

!------------------------------------------------------------------------------

! perform Nose'-Hoover chains for barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

! setting up for ROLL: saving old variables
      roll_tol = -1.0E+10_dp
      IF (simpar%constraint) THEN
        iroll = 1
        CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
        CALL getold ( local_molecules, molecule_set,  &
                      molecule_kind_set, particle_set)
      ENDIF


      SR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! SHAKE-ROLL LOOP

          IF (simpar%constraint) THEN
             CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
          END IF
          CALL get_veps ( cell, npt, simpar, atomic_kind_set, particle_set, &
                          local_molecules, molecule_set, molecule_kind_set, &
                          local_particles, virial, infree, int_group )

!
! first half of velocity verlet
!
          arg_r = ( 0.5_dp * npt ( 1, 1 ) % v * dt ) * &
                  ( 0.5_dp * npt ( 1, 1 ) % v * dt )
          poly_r = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

          arg_v = ( 0.25_dp * npt ( 1, 1 ) % v * dt *  &
                  (1._dp + infree ) )*( 0.25_dp * npt ( 1, 1 ) % v * &
                   dt * ( 1._dp + infree ) )
          poly_v = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

          scale_r = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
          scale_v = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
                       (1._dp + infree ) )
         DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           dm = 0.5_dp*dt/mass
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             vel ( 1, iparticle ) = particle_set(iparticle)%v(1) * scale_v * scale_v + &
                                    scale_v * poly_v * dm * particle_set(iparticle)%f(1)
             vel ( 2, iparticle ) = particle_set(iparticle)%v(2) +  dm * particle_set(iparticle)%f(2)
             vel ( 3, iparticle ) = particle_set(iparticle)%v(3) +  dm * particle_set(iparticle)%f(3)
             pos ( 1, iparticle ) = particle_set(iparticle)%r(1) * scale_r * scale_r + &
                                    scale_r * poly_r * dt * vel ( 1, iparticle )
             pos ( 2, iparticle ) = particle_set(iparticle)%r(2) + dt*vel(2,iparticle)
             pos ( 3, iparticle ) = particle_set(iparticle)%r(3) + dt*vel(3,iparticle)
           END DO
         END DO

        roll_tol = 0._dp
        vector_r ( : ) = 0._dp
        vector_v ( : ) = 0._dp
        vector_r ( 1 ) = scale_r * poly_r
        vector_v ( 1 ) = scale_v * poly_v

        IF ( simpar % constraint ) CALL shake_roll_control( local_molecules, &
                                        molecule_set, molecule_kind_set, &
                                        particle_set, pos, vel, simpar,  &
                                        roll_tol, iroll, vector_r, vector_v,  &
                                        int_group )
  END DO SR

! Update h_mat

       cell % hmat ( 1, 1 ) = cell % hmat ( 1, 1 ) * scale_r * scale_r 


! Update the inverse
      CALL get_hinv ( cell )

!     *** Broadcast the new particle positions ***
      CALL update_particle_set ( particle_set, int_group, pos = pos )

      DEALLOCATE (pos,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")


!     *** Update forces ***

      CALL force_env_set(force_env_p(1)%force_env,box_changed=.TRUE.)
      CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

      DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
        dm = 0.5_dp*dt/mass
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel ( 1, iparticle ) = vel ( 1, iparticle ) * scale_v * scale_v + &
                                 scale_v * poly_v * dm * particle_set(iparticle)%f(1)
          vel ( 2, iparticle ) = vel ( 2, iparticle ) +  dm * particle_set(iparticle)%f(2)
          vel ( 3, iparticle ) = vel ( 3, iparticle ) +  dm * particle_set(iparticle)%f(3)
        END DO
      END DO

      roll_tol = -1.E10_dp
      IF (simpar%constraint) THEN
         first = .TRUE.
         iroll = 1
         CALL set ( atomic_kind_set, particle_set, vel, local_particles, cell, npt, 'F' )
      END IF

      RR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! RATTLE-ROLL LOOP
         roll_tol = 0._dp
         IF ( simpar % constraint ) CALL rattle_roll_setup ( atomic_kind_set, &
                                      particle_set, local_particles, &
                                      molecule_kind_set, molecule_set,  &
                                      local_molecules, vel,  &
                                      cell, npt, simpar, virial, vector_v,  &
                                      roll_tol, iroll, infree, first, int_group )

         CALL get_veps ( cell, npt, simpar, atomic_kind_set, vel, particle_set,&
                         local_molecules, molecule_set, molecule_kind_set, &
                         local_particles, virial, infree, int_group )

      END DO RR

! Nose-Hoover for the Barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

!     *** Broadcast the new particle velocities ***
      CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
      IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                      molecule_set,    & 
                                                      molecule_kind_set, &
                                                      particle_set, virial, &
                                                      "FULL", int_group )

      DEALLOCATE (vel,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

      CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                             local_particles, virial, int_group  )

! Deallocate old variables
    CALL deallocate_old ( )
  END DO

  END SUBROUTINE nph_uniaxial_damped
!!*****
!******************************************************************************
!!****** integrator/set [1.0] *
!!
!!   NAME
!!     set_particle_set
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_particle_set ( atomic_kind_set, particle_set, local_particles, cell , npt, char )


! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles
  TYPE ( cell_type ), POINTER :: cell
  TYPE (npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char

! Locals
  INTEGER :: idim, iparticle_kind, nparticle_kind, nparticle_local, iparticle
  INTEGER :: iparticle_local

!------------------------------------------------------------------------------
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
!!*****
!******************************************************************************
!!****** integrator/set [1.0] *
!!
!!   NAME
!!     set_vel
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_vel ( atomic_kind_set, particle_set, vel, local_particles, cell , npt, char )


! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  REAL (KIND=dp), INTENT ( INOUT ) :: vel ( :, : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles
  TYPE ( cell_type ), POINTER :: cell
  TYPE (npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
  CHARACTER ( LEN = * ), INTENT ( IN ) :: char

! Locals
  INTEGER :: idim, iparticle_kind, nparticle_kind, nparticle_local, iparticle
  INTEGER :: iparticle_local

!------------------------------------------------------------------------------
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
!!*****
!******************************************************************************
!!****** integrator/allocate_old [1.0] *
!!
!!   NAME
!!     set
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE allocate_old ( particle_set, npt )

  IMPLICIT NONE

! Arguments
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( npt_info_type ), POINTER :: npt ( :, : )

! Locals
  INTEGER :: isos, idim, jdim, natoms

!------------------------------------------------------------------------------

  natoms = SIZE ( particle_set)
  idim = SIZE ( npt, 1)
  jdim = SIZE ( npt, 2)

  isos = 0
  IF ( .NOT. ASSOCIATED ( old % v ) ) &
    ALLOCATE ( old % v ( natoms, 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  old % v = 0.0_dp
  IF ( .NOT. ASSOCIATED ( old % r ) ) &
    ALLOCATE ( old % r ( natoms, 3 ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  old % r = 0.0_dp
  IF ( .NOT. ASSOCIATED ( old % eps ) )  &
    ALLOCATE ( old % eps ( idim, jdim ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  old % eps = 0.0_dp
  IF ( .NOT. ASSOCIATED ( old % veps ) ) &
    ALLOCATE ( old % veps ( idim, jdim ),STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  old % veps = 0.0_dp
  IF ( .NOT. ASSOCIATED ( old % h ) ) &
    ALLOCATE ( old%h(3, 3),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )
  old % h = 0.0_dp


END SUBROUTINE allocate_old

!!*****
!******************************************************************************
!!****** integrator/deallocate_old [1.0] *
!!
!!   NAME
!!     set
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE deallocate_old 

  IMPLICIT NONE

! Arguments

! Locals
  INTEGER :: isos 

!------------------------------------------------------------------------------


  isos = 0
  DEALLOCATE ( old % v, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%v', 0 )
  DEALLOCATE ( old % r, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%r', 0 )
  DEALLOCATE ( old % eps, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%eps', 0 )
  DEALLOCATE ( old % veps, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%veps', 0 )
  DEALLOCATE ( old % h, STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'set', 'old%h', 0 )
  NULLIFY ( old % v )
  NULLIFY ( old % r )
  NULLIFY ( old % eps )
  NULLIFY ( old % veps )
  NULLIFY ( old % h )

END SUBROUTINE deallocate_old

!!*****
!******************************************************************************
!!****** integrator/npt_f [1.0] *
!!
!!   NAME
!!     npt_f
!!
!!   SYNOPSIS
!!     Subroutine npt_f(md_env)
!!       Implicit None
!!       Type(md_environment_type), Intent (INOUT):: md_env
!!     End Subroutine npt_f
!!
!!   FUNCTION
!!     Velocity Verlet integrator for the NPT ensemble with fully flexible cell
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

 SUBROUTINE npt_f ( md_env )
  
    IMPLICIT NONE

! Arguments
    TYPE ( md_environment_type ), POINTER :: md_env

! Locals
    CHARACTER(LEN=*), PARAMETER              :: routine_name = "npt_f"
  
    INTEGER :: iparticle, iparticle_kind, iparticle_local, istat, &
               isubsys,  nmolecule, nmolecule_kind, nparticle, &
               nparticle_kind, nparticle_local
    INTEGER, POINTER :: itimes
    INTEGER :: i, j, iroll
    REAL(KIND=dp)                            :: dm, dt, global_mass, &
                                                mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_p_type), DIMENSION(:), POINTER :: force_env_p
    TYPE(lnhc_parameters_type), POINTER      :: nhc_part( : ), nhc_baro ( : )
    TYPE ( npt_info_type ), POINTER          :: npt ( :, : )
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(cell_type), POINTER                 :: cell
    TYPE ( virial_type ), POINTER            :: virial
    REAL (KIND=dp), PARAMETER :: e2 = 1.0_dp/6.0_dp, e4 = e2/20.0_dp, &
       e6 = e4/42.0_dp, e8 = e6/72.0_dp
    REAL (KIND=dp), DIMENSION (3,3) ::  u, uh
    REAL (KIND=dp), DIMENSION (3) :: arg_r, poly_r, scale_r, e_val, ur, uv, uf
    REAL (KIND=dp), DIMENSION (3) :: arg_v, poly_v, scale_v
    REAL (KIND=dp) :: roll_tol, infree, trvg
    REAL (KIND=dp), DIMENSION (3) :: vector_v
    LOGICAL, PARAMETER :: box_change = .TRUE.
    LOGICAL :: first

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env_p=force_env_p,&
                    nhc_part=nhc_part, nhc_baro=nhc_baro, npt = npt, &
                    itimes = itimes )
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree,KIND=dp) 

    IF (SIZE(force_env_p)/=1) THEN
       CALL stop_program("npt_f","integrator",__LINE__,&
       "MD_ENV with more than (or less) than one force_env not supported atm!")
    END IF
    CALL force_env_get(force_env=force_env_p(1)%force_env,subsys=subsys,&
                       virial = virial, cell = cell )

    DO isubsys=1,SIZE(subsys)

      CALL cp_subsys_get(subsys=subsys(isubsys)%subsys,&
                        atomic_kinds=atomic_kinds,&
                        local_particles=local_particles,&
                        particles=particles, &
                        local_molecules_new=local_molecules,&
                        molecules_new=molecules, &
                        molecule_kinds_new=molecule_kinds )

      nparticle_kind = atomic_kinds%n_els
      atomic_kind_set => atomic_kinds%els
      nmolecule_kind = molecule_kinds % n_els
      molecule_kind_set => molecule_kinds % els

      nparticle = particles%n_els
      particle_set => particles%els
      nmolecule = molecules % n_els
      molecule_set => molecules % els

!     *** Allocate work storage for positions and velocities ***

      CALL allocate_old ( particle_set, npt )

      ALLOCATE (pos(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "pos",3*nparticle*dp_size)
      ALLOCATE (vel(3,nparticle),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "vel",3*nparticle*dp_size)
      pos(:,:) = 0.0_dp
      vel(:,:) = 0.0_dp

!------------------------------------------------------------------------------

! perform Nose'-Hoover chains for barostat
      CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )

! perform Nose'-Hoover chains for particles
      CALL lnhc ( nhc_part ( isubsys ), molecule_kind_set, molecule_set, &
                  atomic_kind_set, particle_set, local_molecules, int_group )

! setting up for ROLL: saving old variables
      roll_tol = -1.0E+10_dp
      IF (simpar%constraint) THEN
        iroll = 1
        CALL set ( atomic_kind_set, particle_set,  &
                   local_particles, cell, npt, 'F' )
        CALL getold ( local_molecules, molecule_set,  &
                      molecule_kind_set, particle_set)
      ENDIF

      SR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! SHAKE-ROLL LOOP

          IF (simpar%constraint) THEN
             CALL set ( atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
          END IF
          CALL get_veps ( cell, npt, simpar, atomic_kind_set, particle_set, &
                         local_molecules, molecule_set, molecule_kind_set, &
                         local_particles, virial, infree, int_group )


          trvg = npt ( 1, 1 ) % v  + npt ( 2, 2 ) % v + npt ( 3, 3 ) % v 
!     
! find eigenvalues and eigenvectors of npt ( :, : ) % v
!

          CALL diagonalise ( matrix = npt(:,:)%v, mysize = 3, &
          storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )
!
! first half of velocity verlet
!
          arg_r ( : ) =  0.5_dp * e_val ( : ) * dt * &
                         0.5_dp * e_val ( : ) * dt
          poly_r = 1.0_dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
          scale_r ( : ) = EXP ( 0.5_dp * dt * e_val ( : ) )

          arg_v ( : ) =  0.25_dp * dt * ( e_val ( : ) + trvg * infree ) *&
                         0.25_dp * dt * ( e_val ( : ) + trvg * infree )
          poly_v = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
          scale_v ( : ) = EXP ( -0.25_dp * dt * ( &
                                e_val ( : ) + trvg * infree ) )

          DO iparticle_kind=1,nparticle_kind
            atomic_kind => atomic_kind_set(iparticle_kind)
            CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
            dm = 0.5_dp*dt/mass
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
! Transform positions and velocities and forces
              ur = MATVEC_3X3 ( TRANSPOSE_3D ( u ), particle_set ( iparticle ) % r ( : ) )
              uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), particle_set ( iparticle ) % v ( : ) )
              uf = MATVEC_3X3 ( TRANSPOSE_3D ( u ), particle_set ( iparticle ) % f ( : ) )
!
              uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dm
              uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dm
              uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dm
!
              ur(1) = ur(1)*scale_r(1)*scale_r(1) +  &
                      uv(1)*scale_r(1)*poly_r(1)*dt
              ur(2) = ur(2)*scale_r(2)*scale_r(2) +  &
                      uv(2)*scale_r(2)*poly_r(2)*dt
              ur(3) = ur(3)*scale_r(3)*scale_r(3) +  &
                      uv(3)*scale_r(3)*poly_r(3)*dt
!
              pos ( :, iparticle ) = MATVEC_3x3 ( u, ur )
              vel ( :, iparticle ) = MATVEC_3x3 ( u, uv )
            END DO
          END DO

          roll_tol = 0.0_dp
          IF ( simpar % constraint ) CALL shake_roll_control( local_molecules, &
                                          molecule_set, molecule_kind_set, &
                                          particle_set, pos, vel, simpar,  &
                                          roll_tol, iroll, scale_r*poly_r, &
                                          scale_v*poly_v, int_group, u )
  END DO SR

! Update h_mat

       uh = MATMUL_3X3 ( TRANSPOSE_3D ( u ), cell % hmat )

       DO i = 1, 3
          DO j = 1, 3
            uh ( i, j ) = uh ( i, j ) * scale_r ( i ) * scale_r ( i )
          END DO
       END DO

       cell % hmat = MATMUL_3x3 ( u, uh )
! Update the inverse
       CALL get_hinv ( cell )

!     *** Broadcast the new particle positions ***
       CALL update_particle_set ( particle_set, int_group, pos = pos )

       DEALLOCATE (pos,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"work")

!     *** Update forces ***

       CALL force_env_set(force_env_p(1)%force_env,box_changed=.TRUE.)
       CALL force_env_calc_energy_force(force_env_p(1)%force_env)

!     *** Velocity Verlet (second part) ***

       DO iparticle_kind=1,nparticle_kind
         atomic_kind => atomic_kind_set(iparticle_kind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
         dm = 0.5_dp*dt/mass
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local=1,nparticle_local
           iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
           uv = MATVEC_3X3 ( TRANSPOSE_3D ( u ), vel ( :, iparticle ) )
           uf = MATVEC_3x3 ( TRANSPOSE_3D ( u ), particle_set ( iparticle ) % f ( : ) )
           uv(1) = uv(1)*scale_v(1)*scale_v(1) + uf(1)*scale_v(1)*poly_v(1)*dm
           uv(2) = uv(2)*scale_v(2)*scale_v(2) + uf(2)*scale_v(2)*poly_v(2)*dm
           uv(3) = uv(3)*scale_v(3)*scale_v(3) + uf(3)*scale_v(3)*poly_v(3)*dm
           vel ( :, iparticle ) = MATVEC_3x3(u,uv)
         END DO
       END DO
       roll_tol = -1.E10_dp
       IF (simpar%constraint) THEN
          first = .TRUE.
          iroll = 1
          CALL set ( atomic_kind_set, particle_set, vel, local_particles, cell, npt, 'F' )
       END IF
     RR: DO WHILE (ABS(roll_tol)>=1.E-10_dp) ! RATTLE-ROLL LOOP
        roll_tol = 0.0_dp
        IF ( simpar % constraint ) CALL rattle_roll_setup ( atomic_kind_set, &
                                     particle_set, local_particles, &
                                     molecule_kind_set, molecule_set,  &
                                     local_molecules, vel, &
                                     cell, npt, simpar, virial, vector_v,  &
                                     roll_tol, iroll, infree, first, int_group, u )

        CALL get_veps ( cell, npt, simpar, atomic_kind_set, vel, particle_set, &
                        local_molecules, molecule_set, molecule_kind_set, &
                        local_particles, virial, infree, int_group )

     END DO RR

! Nose-Hoover for the Particles
     CALL lnhc ( nhc_part ( isubsys ), molecule_kind_set, molecule_set, &
                 atomic_kind_set, particle_set, local_molecules, vel, int_group )

! Nose-Hoover for the Barostat
     CALL lnhc ( nhc_baro ( isubsys ), npt, int_group )


!     *** Broadcast the new particle velocities ***
     CALL update_particle_set ( particle_set, int_group, vel = vel )

!     *** Update constraint virial ***
     IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                     molecule_set,    & 
                                                     molecule_kind_set, &
                                                     particle_set, virial, &
                                                     "FULL", int_group )
     DEALLOCATE (vel,STAT=istat)
     IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"vel")

     CALL virial_evaluate ( atomic_kind_set, particle_set,  &
                            local_particles, virial, int_group  )

! Deallocate old variables
    CALL deallocate_old ( )
  END DO


 END SUBROUTINE npt_f

!!*****
!******************************************************************************
!!****** integrator/rattle_roll_setup [1.0] *
!!
!!   NAME
!!     rattle_roll_setup
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE rattle_roll_setup ( atomic_kind_set, particle_set, local_particles, &
                               molecule_kind_set, molecule_set, local_molecules, &
                               vel, cell, npt, simpar, virial, vector_v,  &
                               roll_tol, iroll, infree, first, int_group, u )

  IMPLICIT NONE

! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( molecule_type ), POINTER  :: molecule_set ( : )
  TYPE ( molecule_kind_type ), POINTER  :: molecule_kind_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles, local_molecules
  TYPE ( cell_type ), POINTER :: cell
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (virial_type ), POINTER :: virial
  TYPE (npt_info_type ), DIMENSION ( :, : ), POINTER :: npt
  REAL (KIND=dp), INTENT ( IN ) :: infree
  REAL (KIND=dp), INTENT ( INOUT ) :: vel ( :, : )
  LOGICAL, INTENT ( INOUT ) :: first
  REAL (KIND=dp), INTENT ( OUT ) :: roll_tol
  REAL (KIND=dp), INTENT ( IN ), DIMENSION ( : ) :: vector_v
  INTEGER, INTENT ( INOUT ) :: iroll
  INTEGER, INTENT ( IN ) :: int_group
  REAL (KIND=dp), INTENT ( IN ), OPTIONAL :: u ( : , : )

! Locals
  TYPE ( npt_info_type ), DIMENSION ( 3, 3 ) :: npt_loc

!------------------------------------------------------------------------------

  IF ( first ) THEN
     CALL get_veps ( cell, npt, simpar, atomic_kind_set, vel, particle_set, &
                     local_molecules, molecule_set, molecule_kind_set, &
                     local_particles, virial, infree, int_group )

  END IF
  first = .FALSE.

! assigning local variable
  SELECT CASE ( simpar % ensemble )
  CASE ( 'NPT_I' )
    npt_loc ( :, : ) % v = 0.0_dp
    npt_loc ( :, : ) % mass = 0.0_dp
    npt_loc ( 1, 1 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 2, 2 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 3, 3 ) % v = npt ( 1, 1 ) % v
    npt_loc ( 1, 1 ) % mass = npt ( 1, 1 ) % mass
    npt_loc ( 2, 2 ) % mass = npt ( 1, 1 ) % mass
    npt_loc ( 3, 3 ) % mass = npt ( 1, 1 ) % mass
  CASE ( 'NPT_F' )
    npt_loc = npt
  END SELECT

! resetting

   CALL set ( atomic_kind_set, particle_set, vel, local_particles, cell , npt, 'B' )
   CALL rattle_roll_control( local_molecules, molecule_set, molecule_kind_set, &
                             particle_set, vel, simpar, virial, vector_v,   &
                             npt_loc % v, npt_loc %  mass, roll_tol, iroll,    &
                             int_group, u )

END SUBROUTINE rattle_roll_setup

!!*****
!******************************************************************************
!!****** integrator/get_veps_particle_set [1.0] *
!!
!!   NAME
!!     get_veps_particle_set, get_veps_velocity
!!
!!   FUNCTION
!!     Overloaded routine to compute veps given the particles
!!     structure or a local copy of the velocity array
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_veps_particle_set ( box, npt, simpar, atomic_kind_set, particle_set,  &
                                   local_molecules, molecule_set, molecule_kind_set, &
                                   local_particles, virial, infree, int_group )

! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE ( molecule_type ), POINTER  :: molecule_set ( : )
  TYPE ( molecule_kind_type ), POINTER  :: molecule_kind_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles, local_molecules
  TYPE (cell_type ), INTENT ( IN ) :: box
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (virial_type ), INTENT ( INOUT ) :: virial
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL (KIND=dp), INTENT ( IN ) :: infree
  INTEGER, INTENT ( IN ) :: int_group

! locals
  INTEGER :: i, j, iatom, ii, jj
  INTEGER :: iparticle_kind, nparticle_local, iparticle_local, iparticle
  REAL ( KIND=dp ) :: kin, fdotr,mass,v0,v0i,v,vi
  REAL ( KIND=dp), DIMENSION ( 3, 3 ) :: pv_kin, unit
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
! dbg ISOTROPIC LIMIT
!  REAL (KIND=dp) :: trace
! dbg

! initializing locals
  unit = 0.0_dp
  pv_kin = 0.0_dp
  unit(1,1) = 1.0_dp
  unit(2,2) = 1.0_dp
  unit(3,3) = 1.0_dp
   
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
      END DO
    END DO 
  END DO

  CALL mp_sum( pv_kin, int_group )

  kin = 0.0_dp
  DO ii = 1, 3
     kin = kin +  pv_kin ( ii, ii )
  END DO



! updating the constraint virial
  IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                  molecule_set,    & 
                                                  molecule_kind_set, &
                                                  particle_set, virial, &
                                                  "FULL", int_group )
  IF ( simpar % ensemble == 'NPT_I' ) THEN
! get force on barostat
     fdotr = 0.0_dp
     DO ii = 1, 3
       fdotr = fdotr + virial % pv_virial ( ii, ii ) +  &
                       virial % pv_constraint ( ii, ii )
     END DO

     npt ( :, : ) % f = ( 1.0_dp + ( 3.0_dp * infree ) ) * kin + fdotr - &
       3.0_dp * simpar % p_ext * box % deth
  ELSEIF ( simpar % ensemble == 'NPT_F' ) THEN
     npt ( :, : ) % f = virial % pv_virial ( :, : ) +                   &
       pv_kin ( :, : ) +  virial % pv_constraint ( :, : ) -  &
       unit ( :, : ) * simpar % p_ext * box % deth +       &
       infree * kin * unit ( :, : )
! dbg ISOTROPIC LIMIT
!    trace = npt ( 1, 1 ) % f + npt ( 2, 2 ) % f + npt ( 3, 3 ) % f
!    trace = trace / 3.0_dp
!    npt ( :, : ) % f = trace * unit ( :, : )
! dbg
  ELSEIF ( simpar % ensemble (1:3) == 'NPH' ) THEN
     v = box % deth
     vi = 1._dp/v
     v0 = simpar % v0 
     v0i = 1._dp/v0
!!! orthorhombic box only  !!!
     npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +            &
       pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -      &
       simpar % p0 * v - simpar % v_shock * simpar % v_shock *   &
       v * v0i * ( 1._dp -v * v0i ) + infree * kin 
! dbg UNIAXIAL LIMIT
!     npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +                   &
!       pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -  &
!       simpar % p0 * box % deth + infree * kin 
! dbg UNIAXIAL LIMIT
!!! orthorhombic box only  !!!
  ENDIF

! update barostat velocities
     npt ( :, : ) % v = npt ( :, : ) % v + &
     0.5_dp * simpar % dt * npt ( :, : ) % f / npt ( :, : ) % mass


END SUBROUTINE get_veps_particle_set

!!*****
!******************************************************************************
!!****** integrator/get_veps_velocity[1.0] *
!!
!!   NAME
!!     get_veps_particle_set, get_veps_velocity
!!
!!   FUNCTION
!!     Overloaded routine to compute veps given the particles
!!     structure or a local copy of the velocity array
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_veps_velocity ( box, npt, simpar, atomic_kind_set, vel, particle_set, &
                               local_molecules, molecule_set, molecule_kind_set, &
                               local_particles, virial, infree, int_group )

! Arguments
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind_set ( : )
  REAL (KIND=dp), INTENT ( INOUT ) :: vel ( :, : )
  TYPE ( particle_type ), POINTER ::  particle_set ( : )
  TYPE ( molecule_type ), POINTER  :: molecule_set ( : )
  TYPE ( molecule_kind_type ), POINTER  :: molecule_kind_set ( : )
  TYPE ( distribution_1d_type ), POINTER :: local_particles, local_molecules
  TYPE (cell_type ), INTENT ( IN ) :: box
  TYPE (simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE (virial_type ), INTENT ( INOUT ) :: virial
  TYPE (npt_info_type ), INTENT ( INOUT ), DIMENSION ( :, : ) :: npt
  REAL (KIND=dp), INTENT ( IN ) :: infree
  INTEGER, INTENT ( IN ) :: int_group

! locals
  INTEGER :: i, j, iatom, ii, jj
  INTEGER :: iparticle_kind, nparticle_local, iparticle_local, iparticle
  REAL ( KIND=dp ) :: kin, fdotr,mass,v0,v0i,v,vi
  REAL ( KIND=dp ), DIMENSION ( 3, 3 ) :: pv_kin, unit
  TYPE ( atomic_kind_type ), POINTER :: atomic_kind
! dbg ISOTROPIC LIMIT
!  REAL (KIND=dp) :: trace
! dbg

! initializing locals
  unit = 0.0_dp
  pv_kin = 0.0_dp
  unit(1,1) = 1.0_dp
  unit(2,2) = 1.0_dp
  unit(3,3) = 1.0_dp
   
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
      END DO
    END DO 
  END DO

  CALL mp_sum( pv_kin, int_group )

  kin = 0.0_dp
  DO ii = 1, 3
     kin = kin +  pv_kin ( ii, ii )
  END DO

! updating the constraint virial
  IF ( simpar % constraint ) CALL pv_constraint ( local_molecules, &
                                                  molecule_set,    & 
                                                  molecule_kind_set, &
                                                  particle_set, virial, &
                                                  "FULL", int_group )
 
  IF ( simpar % ensemble == 'NPT_I' ) THEN
! get force on barostat
     fdotr = 0.0_dp
     DO ii = 1, 3
       fdotr = fdotr + virial % pv_virial ( ii, ii ) +  &
                       virial % pv_constraint ( ii, ii )
     END DO

     npt ( :, : ) % f = ( 1.0_dp + ( 3.0_dp * infree ) ) * kin + fdotr - &
       3.0_dp * simpar % p_ext * box % deth
  ELSEIF ( simpar % ensemble == 'NPT_F' ) THEN
     npt ( :, : ) % f = virial % pv_virial ( :, : ) +                   &
       pv_kin ( :, : ) +  virial % pv_constraint ( :, : ) -  &
       unit ( :, : ) * simpar % p_ext * box % deth +       &
       infree * kin * unit ( :, : )
! dbg ISOTROPIC LIMIT
!    trace = npt ( 1, 1 ) % f + npt ( 2, 2 ) % f + npt ( 3, 3 ) % f
!    trace = trace / 3.0_dp
!    npt ( :, : ) % f = trace * unit ( :, : )
! dbg
  ELSEIF ( simpar % ensemble (1:3) == 'NPH' ) THEN
     v = box % deth
     vi = 1._dp/v
     v0 = simpar % v0 
     v0i = 1._dp/v0
     npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +            &
       pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -      &
       simpar % p0 * v - simpar % v_shock * simpar % v_shock *   &
       v * v0i * ( 1._dp -v * v0i ) + infree * kin 
! dbg UNIAXIAL LIMIT
!     npt ( 1, 1 ) % f = virial % pv_virial ( 1, 1 ) +                   &
!       pv_kin ( 1, 1 ) +  virial % pv_constraint ( 1, 1 ) -  &
!       simpar % p0 * box % deth + infree * kin 
! dbg UNIAXIAL LIMIT
  ENDIF

! update barostat velocities
     npt ( :, : ) % v = npt ( :, : ) % v + &
     0.5_dp * simpar % dt * npt ( :, : ) % f / npt ( :, : ) % mass

END SUBROUTINE get_veps_velocity

!******************************************************************************
SUBROUTINE update_particle_set ( particle_set, int_group,  pos, vel )
                                
    TYPE(particle_type), POINTER             :: particle_set( : )
    INTEGER, INTENT(IN)                      :: int_group
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL       :: pos( :, : ), vel( :, : )

    INTEGER                                  :: iparticle, nparticle

  nparticle = SIZE ( particle_set )

  IF ( ( PRESENT ( pos ) ).AND. ( .NOT. PRESENT ( vel ) ) ) THEN
    CALL mp_sum( pos,int_group)
    DO iparticle=1,nparticle
      particle_set(iparticle)%r(:) = pos(:,iparticle)
    END DO
  ELSEIF ( ( PRESENT ( vel ) ) .AND. ( .NOT. PRESENT ( pos ) ) ) THEN
    CALL mp_sum( vel,int_group)
    DO iparticle=1,nparticle
      particle_set(iparticle)%v(:) = vel(:,iparticle)
    END DO
  ELSEIF ( ( PRESENT ( vel ) ) .AND. ( PRESENT ( pos ) ) ) THEN
    CALL mp_sum( pos,int_group)
    CALL mp_sum( vel,int_group)
    DO iparticle=1,nparticle
      particle_set(iparticle)%v(:) = vel(:,iparticle)
      particle_set(iparticle)%r(:) = pos(:,iparticle)
    END DO
  ENDIF

END SUBROUTINE update_particle_set 
!!*****
!******************************************************************************

END MODULE integrator

!******************************************************************************
