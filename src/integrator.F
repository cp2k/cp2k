!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides integrator routines (velocity verlet) for all the
!>      ensemble types
!> \par History
!>      JGH (15-Mar-2001) : Pass logical for box change to force routine
!>      Harald Forbert (Apr-2001): added path integral routine nvt_pimd
!>      CJM (15-Apr-2001) : added coef integrators and energy routines
!>      Joost VandeVondele (Juli-2003): simple version of isokinetic ensemble
!>      Teodoro Laino [tlaino] 10.2007 - University of Zurich: Generalization to
!>                                       different kind of thermostats
!>      Teodoro Laino [tlaino] 11.2007 - Metadynamics: now part of the MD modules
!>      Marcella Iannuzzi      02.2008 - Collecting common code (VV and creation of
!>                                       a temporary type)
!>      Teodoro Laino [tlaino] 02.2008 - Splitting integrator module and keeping in
!>                                       integrator only the INTEGRATORS
!> \author CJM
! *****************************************************************************
MODULE integrator
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type,&
                                             get_hinv,&
                                             parse_cell_line
  USE constraint,                      ONLY: rattle_control,&
                                             shake_control,&
                                             shake_roll_control,&
                                             shake_update_targets
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE constraint_util,                 ONLY: getold,&
                                             pv_constraint
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: broadcast_input_information,&
                                             parser_get_next_line,&
                                             parser_read_line
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_methods,         ONLY: dynamical_coeff_verlet1,&
                                             dynamical_coeff_verlet2
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type
  USE eigenvalueproblems,              ONLY: diagonalise
  USE extended_system_dynamics,        ONLY: shell_scale_comv
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force,&
                                             force_env_set_cell
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: npe_f_ensemble,&
                                             npe_i_ensemble
  USE integrator_utils,                ONLY: &
       allocate_old, allocate_tmp, damp_v, damp_veps, deallocate_old, &
       get_s_ds, old_variables_type, rattle_roll_setup, set, &
       tmp_variables_type, update_dealloc_tmp, update_pv, update_veps, &
       vv_first, vv_second
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matmul_3x3,&
                                             transpose_3d
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max
  USE metadynamics,                    ONLY: metadyn_forces
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: global_constraint_type,&
                                             molecule_type
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             update_particle_set
  USE physcon,                         ONLY: femtoseconds
  USE qmmm_util,                       ONLY: apply_qmmm_walls_reflective
  USE reftraj_types,                   ONLY: reftraj_type
  USE reftraj_util,                    ONLY: compute_msd_reftraj
  USE termination,                     ONLY: stop_program
  USE thermostat_methods,              ONLY: apply_thermostat_baro,&
                                             apply_thermostat_particles,&
                                             apply_thermostat_shells
  USE thermostat_types,                ONLY: thermostat_type
  USE virial_methods,                  ONLY: virial_evaluate
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'integrator'

  PUBLIC :: isokin, langevin, nve, nvt, npt_i, npt_f, nve_respa
  PUBLIC :: nph_uniaxial_damped, nph_uniaxial, reftraj


CONTAINS

! *****************************************************************************
!> \brief Langevin integrator for particle positions & momenta (Brownian dynamics)
!> \par Literature
!>      - A. Ricci and G. Ciccotti, Mol. Phys. 101, 1927-1931 (2003)
!> \par History
!>      - Created (01.07.2005,MK)
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE langevin(md_env,error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'langevin', &
      routineP = moduleN//':'//routineN

    INTEGER :: iparticle, iparticle_kind, iparticle_local, nmolecule, &
      nmolecule_kind, nparticle, nparticle_kind, nparticle_local, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: c, c1, c2, c3, c4, dm, dt, &
                                                gam, mass, sigma
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel, w
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY (cell, para_env, gci, force_env)
    CALL get_md_env(md_env=md_env, simpar=simpar, force_env=force_env,&
         para_env=para_env, error=error)
    
    dt = simpar%dt
    gam = simpar%gamma+simpar%shadow_gamma

    CALL force_env_get(force_env=force_env,subsys=subsys,cell=cell,&
         virial=virial,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds%n_els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules%n_els
    molecule_set => molecules%els

    ! Allocate work storage
    ALLOCATE (pos(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    pos(:,:) = 0.0_dp

    ALLOCATE (vel(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    vel(:,:) = 0.0_dp

    ALLOCATE (w(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    w(:,:) = 0.0_dp

    IF (simpar%constraint) CALL getold(gci,local_molecules,molecule_set,&
         molecule_kind_set, particle_set,cell,error)

    ! Generate random variables
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       nparticle_local = local_particles%n_el(iparticle_kind)
       sigma = simpar%var_w*mass
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          rng_stream => local_particles%local_particle_set(iparticle_kind)%&
               &rng(iparticle_local)%stream
          w(1,iparticle) = next_random_number(rng_stream,variance=sigma,error=error)
          w(2,iparticle) = next_random_number(rng_stream,variance=sigma,error=error)
          w(3,iparticle) = next_random_number(rng_stream,variance=sigma,error=error)
       END DO
    END DO

    ! Apply fix atom constraint
    CALL fix_atom_control(force_env,error,w)

    ! Velocity Verlet (first part)
    c = EXP(-0.25_dp*dt*gam)
    c2 = c*c
    c4 = c2*c2
    c1 = dt*c2

    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       nparticle_local = local_particles%n_el(iparticle_kind)
       dm = 0.5_dp*dt/mass
       c3 = dm/c2
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel(:,iparticle) = particle_set(iparticle)%v(:) +&
               c3*particle_set(iparticle)%f(:)
          pos(:,iparticle) = particle_set(iparticle)%r(:) +&
               c1*particle_set(iparticle)%v(:) +&
               c*dm*(dt*particle_set(iparticle)%f(:) +&
               w(:,iparticle))
       END DO
    END DO

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)

       CALL shake_control(gci,local_molecules,molecule_set,molecule_kind_set,&
            particle_set,pos,vel,dt,simpar%shake_tol,&
            simpar%info_constraint,simpar%lagrange_multipliers,&
            simpar%dump_lm,cell,para_env%group,local_particles, error)
    END IF

    ! Broadcast the new particle positions
    CALL update_particle_set(particle_set,para_env%group,pos=pos)

    DEALLOCATE (pos,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Update forces
    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, vel, error=error)

    ! Update Verlet (second part)
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       dm = 0.5_dp*dt/mass
       c3 = dm/c2
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel(1,iparticle) = vel(1,iparticle) + c3*particle_set(iparticle)%f(1)
          vel(2,iparticle) = vel(2,iparticle) + c3*particle_set(iparticle)%f(2)
          vel(3,iparticle) = vel(3,iparticle) + c3*particle_set(iparticle)%f(3)
          vel(1,iparticle) = c4*vel(1,iparticle) + c2*w(1,iparticle)/mass
          vel(2,iparticle) = c4*vel(2,iparticle) + c2*w(2,iparticle)/mass
          vel(3,iparticle) = c4*vel(3,iparticle) + c2*w(3,iparticle)/mass
       END DO
    END DO

    IF (simpar%constraint) THEN
       CALL rattle_control(gci,local_molecules,molecule_set,molecule_kind_set,&
            particle_set,vel,dt,simpar%shake_tol,&
            simpar%info_constraint,simpar%lagrange_multipliers,&
            simpar%dump_lm,cell,para_env%group,local_particles, error)
    END IF

    ! Broadcast the new particle velocities
    CALL update_particle_set(particle_set,para_env%group,vel=vel)

    DEALLOCATE (vel,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE (w,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Update virial
    IF (simpar%constraint) CALL pv_constraint(gci,local_molecules,molecule_set,&
       molecule_kind_set, particle_set,virial,para_env%group)
  
    CALL virial_evaluate(atomic_kind_set,particle_set,local_particles,&
         virial,para_env%group, error=error)

  END SUBROUTINE langevin

! *****************************************************************************
!> \brief nve integrator for particle positions & momenta
!> \par History
!>      - the local particle lists are used instead of pnode (Sep. 2003,MK)
!>      - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!> \author CJM
! *****************************************************************************
  SUBROUTINE nve ( md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nmolecule, nmolecule_kind, &
                                                nparticle, nparticle_kind, &
                                                nshell
    LOGICAL                                  :: failure, shell_adiabatic, &
                                                shell_present
    REAL(KIND=dp)                            :: dt, max_vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_coeff, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE. 

    NULLIFY (dyn_coeff_set, thermostat_coeff, tmp)
    NULLIFY (subsys, simpar, para_env, cell, gci, force_env, virial)
    NULLIFY (atomic_kinds, local_particles, molecules, molecule_kind_set, molecule_set, particle_set)
    NULLIFY (shell_particles, shell_particle_set, core_particles, &
            core_particle_set, thermostat_shell)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_coeff=thermostat_coeff,thermostat_shell=thermostat_shell,&
         para_env=para_env, error=error)
    dt = simpar%dt
    CALL force_env_get(force_env=force_env, subsys=subsys,&
         dyn_coeff_set=dyn_coeff_set, virial=virial, cell=cell, error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules, &
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules % n_els
    molecule_set => molecules % els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,shell_particles=shell_particles,&
            core_particles=core_particles, error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)

       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(tmp, nparticle,nshell,shell_adiabatic,error=error)

    ! Apply thermostat over the full set of shells if required
    IF(shell_adiabatic) THEN
       CALL apply_thermostat_shells(thermostat_shell, atomic_kind_set, particle_set,&
            local_particles, para_env%group, shell_particle_set=shell_particle_set,&
            core_particle_set=core_particle_set, error=error)
    END IF
    
    IF (simpar%constraint) CALL getold(gci,local_molecules, molecule_set,  &
         molecule_kind_set, particle_set, cell, error)

    ! Velocity Verlet (first part)
    CALL vv_first(tmp,atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, max_vel)

    IF(simpar%variable_dt) THEN
        CALL mp_max(max_vel,para_env%group)
        IF(max_vel * dt  > simpar%dr_tol) THEN
           simpar%dt_fact = simpar%dr_tol/max_vel/dt
           ! repeat the first vv half-step with the rescaled time step 
           dt = dt*simpar%dt_fact
           CALL vv_first(tmp,atomic_kind_set, local_particles, particle_set, &
                core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, max_vel)
        ELSE
           simpar%dt_fact = 1.0_dp
        END IF
    END IF

    IF (simpar%constraint) THEN 
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
           molecule_kind_set, dt, force_env%root_section, error)

       CALL shake_control(gci,local_molecules, molecule_set,&
           molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt,simpar%shake_tol, &
           simpar%info_constraint,simpar%lagrange_multipliers,simpar%dump_lm,&
           cell,para_env%group,local_particles, error )
    END IF

    ! Broadcast the new particle positions and deallocate pos part of temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! coefficients
    CALL dynamical_coeff_verlet1(dyn_coeff_set, thermostat_coeff, para_env, dt, error)
    
    ! Update forces
    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,shell_adiabatic, dt )

    IF (simpar%constraint) CALL rattle_control(gci,local_molecules, molecule_set, &
         molecule_kind_set, particle_set, tmp%vel, dt,simpar%shake_tol,&
         simpar%info_constraint,simpar%lagrange_multipliers,simpar%dump_lm,&
         cell,para_env%group,local_particles, error )

    ! Apply thermostat over the full set of shell if required
    IF (shell_adiabatic) THEN
       CALL apply_thermostat_shells(thermostat_shell,atomic_kind_set, particle_set,&
            local_particles, para_env%group, vel=tmp%vel, &
            shell_vel=tmp%shell_vel, core_vel=tmp%core_vel,&
            error=error)
    END IF

    IF (simpar%annealing) THEN
       tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing
       IF (shell_adiabatic) THEN
          CALL shell_scale_comv(atomic_kind_set,local_particles,particle_set,&
               tmp%vel,tmp%shell_vel,tmp%core_vel)
       END IF
    END IF

    ! Broadcast the new particle velocities and deallocate the full temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    ! Coefficients
    CALL dynamical_coeff_verlet2(dyn_coeff_set, thermostat_coeff, para_env, dt, error)

    ! Update virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial,para_env % group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env % group  ,error=error)

  END SUBROUTINE nve

! *****************************************************************************
!> \brief simplest version of the isokinetic gaussian thermostat
!> \note
!>      - time reversible, and conserves the kinetic energy to machine precision
!>      - is not yet supposed to work for e.g. constraints, our the extended version
!>        of this thermostat
!>        see:
!>         - Zhang F. , JCP 106, 6102 (1997)
!>         - Minary P. et al, JCP 118, 2510 (2003)
!> \par History
!>      - Created [2004-07]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE isokin ( md_env, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'isokin', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nparticle, nparticle_kind, &
                                                nshell
    LOGICAL                                  :: failure, shell_adiabatic, &
                                                shell_present
    REAL(KIND=dp)                            :: ds, dt, max_vel, poly_v(3), &
                                                s, scale_v(3)
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(tmp_variables_type), POINTER        :: tmp

    failure = .FALSE.

    NULLIFY(force_env, tmp)
    NULLIFY(atomic_kinds,para_env,subsys,local_particles)
    NULLIFY(core_particles,particles,shell_particles)
    NULLIFY(core_particle_set,particle_set,shell_particle_set)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         para_env=para_env, error=error)

    dt = simpar%dt

    CALL force_env_get(force_env=force_env,subsys=subsys,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    IF (simpar%constraint) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Constraints not yet implemented")
    ENDIF

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,&
         local_particles=local_particles,&
         particles=particles ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nparticle = particles%n_els
    particle_set => particles%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present, shell_adiabatic=shell_adiabatic)

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,shell_particles=shell_particles,&
            core_particles=core_particles, error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)

       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp (tmp, nparticle, nshell, shell_adiabatic, error=error)

    ! compute s,ds
    CALL get_s_ds(tmp, s, ds, nparticle_kind, atomic_kind_set, local_particles, particle_set,&
         dt, para_env)

    ! Velocity Verlet (first part)
    scale_v(1:3) = SQRT(1.0_dp/ds)
    poly_v(1:3) = 2.0_dp*s/SQRT(ds)/dt
    CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, max_vel,  poly_v=poly_v, scale_v=scale_v)


    IF(simpar%variable_dt) THEN
        CALL mp_max(max_vel,para_env%group)
        IF(max_vel * dt  > simpar%dr_tol) THEN
           simpar%dt_fact = simpar%dr_tol/max_vel/dt
           ! repeat the first vv half-step with the rescaled time step 
           dt = simpar%dt*simpar%dt_fact
           scale_v(1:3) = SQRT(1.0_dp/ds)
           poly_v(1:3) = 2.0_dp*s/SQRT(ds)/dt
           CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
                core_particle_set, shell_particle_set, nparticle_kind,&
                shell_adiabatic, dt, max_vel,  poly_v=poly_v, scale_v=scale_v)
        ELSE
           simpar%dt_fact = 1.0_dp
        END IF
    END IF


    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! compute s,ds
    CALL get_s_ds(tmp, s, ds, nparticle_kind, atomic_kind_set,  local_particles, particle_set,&
       dt, para_env, tmpv=.TRUE.)

    ! Velocity Verlet (second part)
    scale_v(1:3) = SQRT(1.0_dp/ds)
    poly_v(1:3) = 2.0_dp*s/SQRT(ds)/dt
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt,  poly_v=poly_v, scale_v=scale_v)

    IF (simpar%annealing)  tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing

    !  Broadcast the new particle velocities and deallocate the temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

  END SUBROUTINE isokin

! *****************************************************************************
!> \brief nvt integrator for particle positions & momenta
!> \par History
!>      - the local particle lists are used instead of pnode (Sep. 2003,MK)
!>      - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!> \author CJM
! *****************************************************************************
  SUBROUTINE nvt ( md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nvt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nmolecule, nmolecule_kind, &
                                                nparticle, nparticle_kind, &
                                                nshell
    LOGICAL                                  :: failure, shell_adiabatic, &
                                                shell_present
    REAL(KIND=dp)                            :: dt, max_vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_coeff, &
                                                thermostat_part, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY(gci, force_env, dyn_coeff_set, thermostat_coeff, tmp, &
         thermostat_part, thermostat_shell, cell, shell_particles, &
         shell_particle_set, core_particles, core_particle_set)
    NULLIFY(para_env, subsys, local_molecules, local_particles, molecule_kinds,&
         molecules, molecule_kind_set, molecule_set)
    NULLIFY(simpar, thermostat_coeff, thermostat_part, thermostat_shell)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_part=thermostat_part, thermostat_coeff=thermostat_coeff,&
         thermostat_shell=thermostat_shell,para_env=para_env,&
         error=error)
    dt = simpar%dt

    CALL force_env_get(force_env=force_env, subsys=subsys, dyn_coeff_set=dyn_coeff_set,&
         virial=virial, cell=cell,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind  = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind  = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle    = particles%n_els
    particle_set => particles%els
    nmolecule    = molecules % n_els
    molecule_set => molecules % els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)

    !  Allocate work storage for positions and velocities
    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,shell_particles=shell_particles,&
            core_particles=core_particles, error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)

       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(tmp, nparticle,nshell,shell_adiabatic,error=error)

    ! Apply Thermostat over the full set of particles
    IF(shell_adiabatic) THEN
       CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
            particle_set, local_molecules, para_env % group, shell_adiabatic=shell_adiabatic,&
            shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
            error=error)

       CALL apply_thermostat_shells(thermostat_shell, atomic_kind_set, particle_set,&
            local_particles, para_env%group, shell_particle_set=shell_particle_set,&
            core_particle_set=core_particle_set, error=error)
    ELSE
       CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
            particle_set, local_molecules, para_env%group ,error=error)
    END IF

    IF (simpar%constraint) CALL getold( gci, local_molecules, molecule_set,  &
         molecule_kind_set, particle_set, cell, error)

   ! Velocity Verlet (first part)
    CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, max_vel)

    IF(simpar%variable_dt) THEN
        CALL mp_max(max_vel,para_env%group)
        IF(max_vel * dt  > simpar%dr_tol) THEN
           simpar%dt_fact = simpar%dr_tol/max_vel/dt
           ! repeat the first vv half-step with the rescaled time step 
           dt = simpar%dt*simpar%dt_fact
           CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, max_vel)

        ELSE
           simpar%dt_fact = 1.0_dp
        END IF
    END IF

    IF (simpar%constraint) THEN 
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
           molecule_kind_set, dt, force_env%root_section, error)

       CALL shake_control( gci, local_molecules, molecule_set,   &
           molecule_kind_set, particle_set,tmp%pos, tmp%vel, dt,simpar%shake_tol,&
           simpar%info_constraint,simpar%lagrange_multipliers,simpar%dump_lm,&
           cell, para_env%group,local_particles,error )
    END IF

    ! Broadcast the new particle positions and deallocate pos components of temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! coefficients
    CALL dynamical_coeff_verlet1 ( dyn_coeff_set, thermostat_coeff, para_env, dt, error)

    ! Update forces
    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt)

    IF (simpar%constraint) CALL rattle_control( gci, local_molecules, molecule_set, &
         molecule_kind_set, particle_set, tmp%vel, dt,simpar%shake_tol,&
         simpar%info_constraint,simpar%lagrange_multipliers,simpar%dump_lm,&
         cell,para_env%group,local_particles, error )

    ! Apply Thermostat over the full set of particles
    IF (shell_adiabatic) THEN
       CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
            particle_set, local_molecules, para_env%group, shell_adiabatic=shell_adiabatic,&
            vel= tmp%vel, shell_vel= tmp%shell_vel, core_vel= tmp%core_vel, error=error)

       CALL apply_thermostat_shells(thermostat_shell,atomic_kind_set, particle_set,&
            local_particles, para_env%group, vel= tmp%vel, shell_vel= tmp%shell_vel, &
            core_vel= tmp%core_vel, error=error)
    ELSE
       CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
            particle_set, local_molecules, para_env%group, vel= tmp%vel, error=error)
    END IF

    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )


    ! Update constraint virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial,para_env % group )

    ! coefficients
    CALL dynamical_coeff_verlet2 ( dyn_coeff_set, thermostat_coeff, para_env, dt, error)

    !     **  Evaluate Virial
    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env % group, error=error)

  END SUBROUTINE nvt

! *****************************************************************************
!> \brief npt_i integrator for particle positions & momenta
!>      isotropic box changes
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE npt_i ( md_env, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'npt_i', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: box_change = .TRUE.
    REAL(KIND=dp), PARAMETER                 :: e2 = 1.0_dp/6.0_dp, &
                                                e4 = e2/20.0_dp, &
                                                e6 = e4/42.0_dp, &
                                                e8 = e6/72.0_dp

    INTEGER                                  :: iroll, nmolecule, &
                                                nmolecule_kind, nparticle, &
                                                nparticle_kind, nshell
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, shell_present
    REAL(KIND=dp)                            :: arg_r, arg_v, dt, infree, &
                                                kin, max_vel, poly_r(3), &
                                                poly_v(3), roll_tol, &
                                                scale_r(3), scale_v(3)
    REAL(KIND=dp), DIMENSION(3)              :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin
    REAL(KIND=dp), SAVE                      :: eps_0
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(npt_info_type), POINTER             :: npt(:,:)
    TYPE(old_variables_type), POINTER        :: old
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_baro, &
                                                thermostat_part, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE. 

    NULLIFY(gci,thermostat_baro,thermostat_part,thermostat_shell,force_env)
    NULLIFY(atomic_kinds,cell,para_env,subsys,local_molecules,local_particles)
    NULLIFY(molecule_kinds,molecules,molecule_kind_set,npt)
    NULLIFY(core_particles,particles,shell_particles, tmp, old)
    NULLIFY(core_particle_set,particle_set,shell_particle_set)
    NULLIFY(simpar,virial)
   
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_part=thermostat_part,&
         thermostat_baro=thermostat_baro, thermostat_shell=thermostat_shell, npt=npt,&
         itimes=itimes, first_time=first_time, para_env=para_env, error=error)
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree,KIND=dp)

    CALL force_env_get(force_env=force_env, subsys=subsys, virial=virial,&
         cell=cell, error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules, &
         gci=gci,molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules % n_els
    molecule_set => molecules % els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env % group, error=error)
    END IF

    ! Allocate work storage for positions and velocities
    CALL allocate_old ( old, particle_set, npt, error=error )

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,&
            shell_particles=shell_particles, core_particles=core_particles,&
            error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)
       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(tmp, nparticle, nshell, shell_adiabatic, error=error)

    ! Initialize eps_0 the first time through
    IF ( first_time ) eps_0 = npt (1,1) % eps

    ! Apply thermostat to  barostat
    CALL apply_thermostat_baro( thermostat_baro, npt, para_env%group, error)

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_i_ensemble) THEN
      IF(shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, para_env%group,shell_adiabatic=shell_adiabatic,&
             shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
             error=error)

      ELSE
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, para_env%group ,error=error)
      END IF
    END IF

    ! Apply Thermostat over the core-shell motion 
    CALL apply_thermostat_shells(thermostat_shell, atomic_kind_set, particle_set,&
         local_particles, para_env%group, shell_particle_set=shell_particle_set,&
         core_particle_set=core_particle_set, error=error)

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)
    END IF

    ! setting up for ROLL: saving old variables
    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ENDIF

    SR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! SHAKE-ROLL LOOP

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF

       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       arg_r = ( 0.5_dp * npt ( 1, 1 ) % v * dt ) * &
            ( 0.5_dp * npt ( 1, 1 ) % v * dt )
       poly_r(1:3) = 1.0_dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4

       arg_v = ( 0.25_dp * npt ( 1, 1 ) % v * dt *  &
            (1.0_dp + 3.0_dp * infree ) )*( 0.25_dp * npt ( 1, 1 ) % v * &
            dt * ( 1.0_dp + 3.0_dp * infree ) )
       poly_v(1:3) = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

       scale_r(1:3) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
       scale_v(1:3) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
            (1.0_dp + 3.0_dp * infree ) )

       ! first half of velocity verlet
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, max_vel, poly_r=poly_r, poly_v=poly_v, &
            scale_r=scale_r, scale_v=scale_v )

       IF(simpar%variable_dt) THEN
           CALL mp_max(max_vel, para_env%group)
           IF(max_vel*dt > simpar%dr_tol) THEN
              simpar%dt_fact = simpar%dr_tol/max_vel/dt
              dt = dt*simpar%dt_fact
              arg_r = arg_r * simpar%dt_fact*simpar%dt_fact
              poly_r(1:3) = 1.0_dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
              arg_v = arg_v * simpar%dt_fact*simpar%dt_fact
              poly_v(1:3) = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
              scale_r(1:3) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
              scale_v(1:3) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
                            (1.0_dp + 3.0_dp * infree ) )
              CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
                  core_particle_set, shell_particle_set, nparticle_kind,&
                  shell_adiabatic, dt, max_vel, poly_r=poly_r, poly_v=poly_v, &
                  scale_r=scale_r, scale_v=scale_v )
           ELSE
              simpar%dt_fact = 1.0_dp
           END IF
       END IF

       roll_tol = 0.0_dp
       vector_r ( : ) = scale_r(:) * poly_r(:)
       vector_v ( : ) = scale_v(:) * poly_v(:)

       IF ( simpar % constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt, simpar,  &
            roll_tol, iroll, vector_r, vector_v, para_env % group, cell=cell, &
            local_particles=local_particles, error=error )
    END DO SR

    ! Update eps:
    npt ( :, : ) % eps =  npt ( :, : ) % eps + dt * npt ( :, : ) % v

    ! Update h_mat
    cell % hmat ( :, : ) = cell % hmat ( :, : ) &
         * EXP ( npt ( 1, 1 ) % eps - eps_0 )

    eps_0 =  npt ( 1, 1 ) % eps

    ! Update the inverse
    CALL get_hinv ( cell )

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! Update forces
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, poly_v=poly_v, scale_v=scale_v )

    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       first = .TRUE.
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    END IF

    RR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! RATTLE-ROLL LOOP
       roll_tol = 0.0_dp
       IF ( simpar % constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set,  &
            local_molecules, tmp%vel, dt,  cell, npt, simpar, virial, vector_v,  &
            roll_tol, iroll, infree, first, para_env,error=error)

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

    END DO RR

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_i_ensemble) THEN
      IF (shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, para_env%group, shell_adiabatic=shell_adiabatic,&
             vel=tmp%vel, shell_vel=tmp%shell_vel, core_vel=tmp%core_vel, error=error)
      ELSE
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, para_env%group, vel=tmp%vel, error=error)
      END IF
    END IF


    ! Apply Thermostat over the core-shell motion 
    IF(ASSOCIATED(thermostat_shell)) THEN
      CALL apply_thermostat_shells(thermostat_shell,atomic_kind_set, particle_set,&
           local_particles, para_env%group, vel=tmp%vel, shell_vel=tmp%shell_vel, &
           core_vel=tmp%core_vel, error=error)
    END IF

    ! Apply Thermostat to Barostat
    CALL apply_thermostat_baro( thermostat_baro, npt, para_env%group, error)

    ! Annealing of particle velocities is only possible when no thermostat is active
    IF (simpar% ensemble == npe_i_ensemble .AND. simpar%annealing) THEN
       tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing
       IF (shell_adiabatic) THEN
          CALL shell_scale_comv(atomic_kind_set,local_particles,particle_set,&
               tmp%vel,tmp%shell_vel,tmp%core_vel)
       END IF
    END IF
    ! Annealing of CELL velocities is only possible when no thermostat is active
    IF (simpar% ensemble == npe_i_ensemble .AND. simpar%annealing_cell) THEN
       npt(1,1)%v = npt(1,1)%v * simpar%f_annealing_cell
    END IF


    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    ! Update constraint virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set, molecule_kind_set, particle_set, virial, para_env % group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env % group , error=error)

    ! Deallocate old variables
    CALL deallocate_old (old, error=error )

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE npt_i

! *****************************************************************************
!> \brief uses coordinates in a file and generates frame after frame of these
!> \note
!>     it can be used to compute some properties on already available trajectories
!> \par History
!>      - 04.2005 created [Joost VandeVondele]
!>      - modified to make it more general [MI] 
! *****************************************************************************
  SUBROUTINE reftraj ( md_env, error )
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: AA
    INTEGER                                  :: cell_itimes, I, nparticle, &
                                                Nread, stat, trj_itimes
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, init, my_end
    REAL(KIND=dp)                            :: cell_time, h(3,3), trj_epot, &
                                                trj_time, vol
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(reftraj_type), POINTER              :: reftraj_env
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    failure = .FALSE.

    NULLIFY(reftraj_env, particle_set, particles, force_env, subsys, simpar, para_env, cell)
    CALL get_md_env(md_env=md_env, init=init, reftraj=reftraj_env, force_env=force_env,&
         para_env=para_env, simpar=simpar, error=error)

    CALL force_env_get(force_env=force_env, cell=cell, subsys=subsys, &
         error=error)
    reftraj_env%isnap = reftraj_env%isnap + reftraj_env%info%stride

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)
    CALL cp_subsys_get(subsys=subsys, particles=particles, error=error)
    nparticle    = particles%n_els
    particle_set => particles%els
    ALLOCATE (pos(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    pos = 0.0_dp
    
    ! SnapShots read from an external file
    IF(para_env%mepos==para_env%source) THEN
       CALL parser_read_line(reftraj_env%info%traj_parser,1,error=error)
       READ(reftraj_env%info%traj_parser%input_line,FMT="(I8)") nread
       CPPostcondition(nread==nparticle,cp_failure_level,routineP,error,failure)
       CALL parser_read_line(reftraj_env%info%traj_parser,1,error=error)
       trj_itimes = HUGE(0)
       READ(reftraj_env%info%traj_parser%input_line,FMT="(T6,I8,T23,F12.3,T41,F20.10)",&
            ERR=999)trj_itimes, trj_time, trj_epot
999    IF (trj_itimes==HUGE(0)) THEN
          CALL get_md_env(md_env, itimes=itimes, error=error)
          trj_itimes = itimes
          trj_time   = 0.0_dp
          trj_epot   = 0.0_dp
       END IF       
       DO i = 1,nread-1
          CALL parser_read_line(reftraj_env%info%traj_parser,1,error=error)
          READ(reftraj_env%info%traj_parser%input_line,*) AA,pos(1:3,i)
          pos(1,i) = cp_unit_to_cp2k(pos(1,i),"angstrom",error=error)
          pos(2,i) = cp_unit_to_cp2k(pos(2,i),"angstrom",error=error)
          pos(3,i) = cp_unit_to_cp2k(pos(3,i),"angstrom",error=error)
       END DO
       ! End of file is properly addressed in the previous call..
       ! Let's check directly (providing some info) also for the last
       ! line of this frame..
       CALL parser_read_line(reftraj_env%info%traj_parser,1,at_end=my_end,error=error)
       READ(reftraj_env%info%traj_parser%input_line,*) AA,pos(1:3,i)
       pos(1,i) = cp_unit_to_cp2k(pos(1,i),"angstrom",error=error)
       pos(2,i) = cp_unit_to_cp2k(pos(2,i),"angstrom",error=error)
       pos(3,i) = cp_unit_to_cp2k(pos(3,i),"angstrom",error=error)
    END IF
    CALL mp_bcast(my_end,para_env%source,para_env%group)
    ! Check if we reached the end of the file and provide some info..
    IF (my_end) THEN
       CALL cp_assert(reftraj_env%isnap==(simpar%nsteps-1),cp_fatal_level,cp_assertion_failed,routineP,&
            "Reached the end of the Trajectory  frames in the TRAJECTORY file. Number of "//&
            "missing frames ("//cp_to_string((simpar%nsteps-1)-reftraj_env%isnap)//")."//& 
CPSourceFileRef,&
            only_ionode=.TRUE.)
    END IF
    CALL broadcast_input_information(reftraj_env%info%traj_parser,error=error)
    CALL mp_bcast (trj_itimes,para_env%source,para_env%group)
    CALL mp_bcast (trj_time,para_env%source,para_env%group)
    CALL mp_bcast (trj_epot,para_env%source,para_env%group)
    CALL update_particle_set (particle_set, para_env%group, pos=pos)
    DEALLOCATE (pos, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    IF(reftraj_env%info%variable_volume)THEN
       CALL parser_get_next_line(reftraj_env%info%cell_parser,1,at_end=my_end,error=error)
       CALL parse_cell_line(reftraj_env%info%cell_parser%input_line, cell_itimes, cell_time, h, vol, error)
       CPPostcondition(trj_itimes==cell_itimes,cp_failure_level,routineP,error,failure)
       ! Check if we reached the end of the file and provide some info..
       IF (my_end) THEN
          CALL cp_assert(reftraj_env%isnap==(simpar%nsteps-1),cp_fatal_level,cp_assertion_failed,routineP,&
               "Reached the end of the cell info frames in the CELL file. Number of "//&
               "missing frames ("//cp_to_string((simpar%nsteps-1)-reftraj_env%isnap)//")."//& 
CPSourceFileRef,&
               only_ionode=.TRUE.)
       END IF
    END IF
   
    IF (init) THEN
       reftraj_env%time0   = trj_time
       reftraj_env%epot0   = trj_epot
       reftraj_env%itimes0 = trj_itimes 
    END IF
    
    IF (trj_itimes/=0.0_dp.AND.trj_time/=0.0_dp) simpar%dt = (trj_time/femtoseconds)/REAL(trj_itimes,KIND=dp)
    
    reftraj_env%epot   = trj_epot
    reftraj_env%itimes = trj_itimes
    reftraj_env%time   = trj_time/femtoseconds
    
    IF(reftraj_env%info%variable_volume)THEN
       cell%hmat(1:3,1) = h(1:3,1)
       cell%hmat(1:3,2) = h(1:3,2)
       cell%hmat(1:3,3) = h(1:3,3)
       CALL get_hinv(cell)
       CALL force_env_set_cell(force_env, cell, error)
    END IF
    
    ! Task to perform on the reference trajectory
    ! Compute energy and forces
    CALL force_env_calc_energy_force(force_env,eval_energy_forces=reftraj_env%info%eval_EF,error=error)
    CALL metadyn_forces(force_env, error=error)
    
    ! Compute MSD with respect to a reference configuration
    IF(reftraj_env%info%msd) THEN
       CALL compute_msd_reftraj(reftraj_env,md_env, particle_set,error=error)
    END IF

    ! Skip according the stride both Trajectory and Cell (if possible)
    CALL parser_get_next_line(reftraj_env%info%traj_parser,(reftraj_env%info%stride-1)*(nparticle+2),error=error)
    IF(reftraj_env%info%variable_volume)THEN
       CALL parser_get_next_line(reftraj_env%info%cell_parser,(reftraj_env%info%stride-1),error=error)
    END IF
    
  END SUBROUTINE reftraj

! *****************************************************************************
!> \brief nph_uniaxial integrator (non-Hamiltonian version)
!>      for particle positions & momenta undergoing
!>      uniaxial stress ( in x-direction of orthorhombic cell)
!>      due to a shock compression:
!>      Reed et. al. Physical Review Letters 90, 235503 (2003).
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE nph_uniaxial ( md_env, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nph_uniaxial', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: box_change = .TRUE.
    REAL(dp), PARAMETER                      :: e2 = 1._dp/6._dp, &
                                                e4 = e2/20._dp, &
                                                e6 = e4/42._dp, e8 = e6/72._dp

    INTEGER                                  :: iroll, nmolecule, &
                                                nmolecule_kind, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, shell_present
    REAL(KIND=dp) :: arg_r, arg_v, arg_vx, dt, infree, kin, max_vel, &
      poly_r(3), poly_v(3), roll_tol, scale_r(3), scale_v(3)
    REAL(KIND=dp), DIMENSION(3)              :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(old_variables_type), POINTER        :: old
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY(gci, force_env)
    NULLIFY(atomic_kinds,cell,para_env,subsys,local_molecules,local_particles)
    NULLIFY(molecule_kinds,molecules,molecule_kind_set,npt)
    NULLIFY(core_particles,particles,shell_particles, tmp, old)
    NULLIFY(core_particle_set,particle_set,shell_particle_set)
    NULLIFY(simpar,virial)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,npt = npt,&
         first_time=first_time, para_env=para_env, error=error)
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree, dp )

    CALL force_env_get(force_env, subsys=subsys, virial = virial, cell=cell ,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules, molecules_new=molecules, gci=gci,&
         molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules % n_els
    molecule_set => molecules % els

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env % group  ,error=error)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)

    ! Allocate work storage for positions and velocities
    CALL allocate_old ( old, particle_set, npt, error=error )

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,&
            shell_particles=shell_particles, core_particles=core_particles,&
            error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)
       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(tmp,nparticle, nshell, shell_adiabatic, error=error)

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)
    END IF

    ! setting up for ROLL: saving old variables
    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ENDIF

    SR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! SHAKE-ROLL LOOP

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF
       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       arg_r = ( 0.5_dp * npt ( 1, 1 ) % v * dt ) * &
            ( 0.5_dp * npt ( 1, 1 ) % v * dt )
       poly_r(1) = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
       poly_r(2) = 1.0_dp
       poly_r(3) = 1.0_dp

       arg_vx = ( 0.25_dp * npt ( 1, 1 ) % v * dt *  &
            (1._dp + infree ) )*( 0.25_dp * npt ( 1, 1 ) % v * &
            dt * ( 1._dp + infree ) )
       arg_v = ( 0.25_dp * npt ( 1, 1 ) % v * dt * infree) * &
            ( 0.25_dp * npt ( 1, 1 ) % v * dt * infree )
       poly_v(1) = 1._dp + e2*arg_vx + e4*arg_vx*arg_vx + e6*arg_vx**3 + e8*arg_vx**4
       poly_v(2) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
       poly_v(3) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

       scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
       scale_r(2) = 1.0_dp
       scale_r(3) = 1.0_dp

       scale_v(1) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
            (1._dp + infree ) )
       scale_v(2) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * infree )
       scale_v(3) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * infree )

       ! first half of velocity verlet
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, max_vel, poly_r=poly_r, poly_v=poly_v, &
            scale_r=scale_r, scale_v=scale_v )

       IF(simpar%variable_dt) THEN
           CALL mp_max(max_vel, para_env%group)
           IF(max_vel*dt > simpar%dr_tol) THEN
              simpar%dt_fact = simpar%dr_tol/max_vel/dt
              dt = dt*simpar%dt_fact
              arg_r = arg_r * simpar%dt_fact*simpar%dt_fact
              poly_r(1) = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
              arg_v = arg_v * simpar%dt_fact*simpar%dt_fact
              poly_v(1) = 1._dp + e2*arg_vx + e4*arg_vx*arg_vx + e6*arg_vx**3 + e8*arg_vx**4
              poly_v(2) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
              poly_v(3) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
              scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
              scale_v(1) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * &
                   (1._dp + infree ) )
              scale_v(2) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * infree )
              scale_v(3) = EXP( -0.25_dp * dt * npt ( 1, 1 ) % v * infree )
              CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
                  core_particle_set, shell_particle_set, nparticle_kind,&
                  shell_adiabatic, dt, max_vel, poly_r=poly_r, poly_v=poly_v, &
                  scale_r=scale_r, scale_v=scale_v )
           ELSE
              simpar%dt_fact = 1.0_dp
           END IF
       END IF

       roll_tol = 0._dp
       vector_r ( : ) = 0._dp
       vector_v ( : ) = scale_v(:) * poly_v(:)
       vector_r ( 1 ) = scale_r(1) * poly_r(1)

       IF ( simpar % constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt, simpar,&
            roll_tol, iroll, vector_r, vector_v, para_env % group, cell=cell, &
            local_particles=local_particles, error=error )
    END DO SR

    ! Update h_mat
    cell % hmat ( 1, 1 ) = cell % hmat ( 1, 1 ) * scale_r(1) * scale_r(1)

    ! Update the inverse
    CALL get_hinv ( cell )

    ! Broadcast the new particle positions and deallocate the pos component of temporary
    CALL update_dealloc_tmp (tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! Update forces (and stress)
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env, error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
     CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
          core_particle_set, shell_particle_set, nparticle_kind,&
          shell_adiabatic, dt, poly_v=poly_v, scale_v=scale_v )

    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       first = .TRUE.
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    END IF

    RR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! RATTLE-ROLL LOOP
       roll_tol = 0._dp
       IF ( simpar % constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set,  &
            local_molecules, tmp%vel, dt, cell, npt, simpar, virial, vector_v,  &
            roll_tol, iroll, infree, first, para_env, error=error )

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )
    END DO RR

    IF (simpar%annealing)  tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing

    ! Broadcast the new particle velocities and deallocate the temporary
    CALL update_dealloc_tmp (tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    ! Update constraint virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial, para_env % group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env % group, error=error)

    ! Deallocate old variables
    CALL deallocate_old ( old, error=error )

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE nph_uniaxial

! *****************************************************************************
!> \brief nph_uniaxial integrator (non-Hamiltonian version)
!>      for particle positions & momenta undergoing
!>      uniaxial stress ( in x-direction of orthorhombic cell)
!>      due to a shock compression:
!>      Reed et. al. Physical Review Letters 90, 235503 (2003).
!>      Added damping (e.g. thermostat to barostat)
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE nph_uniaxial_damped ( md_env, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nph_uniaxial_damped', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: box_change = .TRUE.
    REAL(dp), PARAMETER                      :: e2 = 1._dp/6._dp, &
                                                e4 = e2/20._dp, &
                                                e6 = e4/42._dp, e8 = e6/72._dp

    INTEGER                                  :: iroll, nmolecule, &
                                                nmolecule_kind, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, shell_present
    REAL(KIND=dp) :: aa, aax, arg_r, arg_v, arg_vx, dt, gamma1, infree, kin, &
      max_vel, poly_r(3), poly_v(3), roll_tol, scale_r(3), scale_v(3)
    REAL(KIND=dp), DIMENSION(3)              :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(old_variables_type), POINTER        :: old
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY(gci, force_env)
    NULLIFY(atomic_kinds,cell,para_env,subsys,local_molecules,local_particles)
    NULLIFY(molecule_kinds,molecules,molecule_kind_set,npt)
    NULLIFY(core_particles,particles,shell_particles,tmp,old)
    NULLIFY(core_particle_set,particle_set,shell_particle_set)
    NULLIFY(simpar,virial)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,npt = npt,&
         first_time=first_time, para_env=para_env, error=error)
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree, dp )
    gamma1 = simpar % gamma_nph

    CALL force_env_get(force_env, subsys=subsys, virial = virial, cell = cell ,error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,gci=gci,&
         molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules % n_els
    molecule_set => molecules % els

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env % group , error=error)
    END IF

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)

    ! Allocate work storage for positions and velocities 
    CALL allocate_old ( old, particle_set, npt, error=error )

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,&
            shell_particles=shell_particles, core_particles=core_particles,&
            error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)
       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(tmp,nparticle, nshell, shell_adiabatic, error=error)

    ! perform damping on velocities
    CALL damp_v ( molecule_kind_set, molecule_set, particle_set, local_molecules,&
         gamma1, npt(1,1), dt, para_env % group )

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)
    END IF

    ! setting up for ROLL: saving old variables
    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ENDIF

    SR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! SHAKE-ROLL LOOP

       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF
       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

       arg_r = ( 0.5_dp * npt ( 1, 1 ) % v * dt ) * &
            ( 0.5_dp * npt ( 1, 1 ) % v * dt )
       poly_r(1) = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
       poly_r(2) = 1.0_dp
       poly_r(3) = 1.0_dp

       aax = npt ( 1, 1 ) % v * ( 1.0_dp + infree )
       arg_vx = ( 0.25_dp * dt * aax ) * ( 0.25_dp * dt * aax )
       poly_v(1) = 1._dp + e2*arg_vx + e4*arg_vx*arg_vx + e6*arg_vx**3 + e8*arg_vx**4

       aa = npt ( 1, 1 ) % v *  infree
       arg_v = ( 0.25_dp * dt * aa ) * ( 0.25_dp * dt * aa )
       poly_v(2) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
       poly_v(3) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4

       scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
       scale_r(2) = 1.0_dp
       scale_r(3) = 1.0_dp
       scale_v(1) = EXP( -0.25_dp * dt * aax )
       scale_v(2) = EXP( -0.25_dp * dt * aa )
       scale_v(3) = EXP( -0.25_dp * dt * aa )

       ! first half of velocity verlet
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, max_vel, poly_r=poly_r, poly_v=poly_v, &
            scale_r=scale_r, scale_v=scale_v )


       IF(simpar%variable_dt) THEN
           CALL mp_max(max_vel, para_env%group)
           IF(max_vel*dt > simpar%dr_tol) THEN
              simpar%dt_fact = simpar%dr_tol/max_vel/dt
              dt = dt*simpar%dt_fact
              arg_r = arg_r * simpar%dt_fact*simpar%dt_fact
              poly_r(1) = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
              arg_vx = ( 0.25_dp * dt * aax ) * ( 0.25_dp * dt * aax )
              poly_v(1) = 1._dp + e2*arg_vx + e4*arg_vx*arg_vx + e6*arg_vx**3 + e8*arg_vx**4
              arg_v = ( 0.25_dp * dt * aa ) * ( 0.25_dp * dt * aa )
              poly_v(2) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
              poly_v(3) = 1._dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
       
              scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 ) % v )
              scale_v(1) = EXP( -0.25_dp * dt * aax )
              scale_v(2) = EXP( -0.25_dp * dt * aa )
              scale_v(3) = EXP( -0.25_dp * dt * aa )
              CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
                  core_particle_set, shell_particle_set, nparticle_kind,&
                  shell_adiabatic, dt, max_vel, poly_r=poly_r, poly_v=poly_v, &
                  scale_r=scale_r, scale_v=scale_v )
           ELSE
              simpar%dt_fact = 1.0_dp
           END IF
        END IF

       roll_tol = 0._dp
       vector_r ( : ) = 0._dp
       vector_v ( : ) = scale_v(:) * poly_v(:)
       vector_r ( 1 ) = scale_r(1) * poly_r(1)

       IF ( simpar % constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos,tmp%vel, dt, simpar,&
            roll_tol, iroll, vector_r, vector_v, para_env % group, cell=cell,&
            local_particles=local_particles, error=error )
    END DO SR

    ! Update h_mat
    cell % hmat ( 1, 1 ) = cell % hmat ( 1, 1 ) * scale_r(1) * scale_r(1)

    ! Update the inverse
    CALL get_hinv ( cell )

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! Update forces
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env, error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
     CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
          core_particle_set, shell_particle_set, nparticle_kind,&
          shell_adiabatic, dt, poly_v=poly_v, scale_v=scale_v )

    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       first = .TRUE.
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    END IF

    RR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! RATTLE-ROLL LOOP
       roll_tol = 0._dp
       IF ( simpar % constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set, local_molecules,&
            tmp%vel, dt, cell, npt, simpar, virial, vector_v, roll_tol, iroll, infree, first,&
            para_env, error=error)
       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

    END DO RR

    ! perform damping on velocities
    CALL damp_v ( molecule_kind_set, molecule_set, particle_set, local_molecules,&
         tmp%vel, gamma1,  npt ( 1, 1 ), dt, para_env % group )

    IF (simpar%annealing)  tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing

    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    ! Update constraint virial
    IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set, molecule_kind_set, particle_set, virial, para_env % group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env % group, error=error)

    ! Deallocate old variables
    CALL deallocate_old ( old, error=error)

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE nph_uniaxial_damped

! *****************************************************************************
!> \brief Velocity Verlet integrator for the NPT ensemble with fully flexible cell
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE npt_f ( md_env, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'npt_f', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: box_change = .TRUE.
    REAL(KIND=dp), PARAMETER                 :: e2 = 1.0_dp/6.0_dp, &
                                                e4 = e2/20.0_dp, &
                                                e6 = e4/42.0_dp, &
                                                e8 = e6/72.0_dp

    INTEGER                                  :: i, iroll, j, nmolecule, &
                                                nmolecule_kind, nparticle, &
                                                nparticle_kind, nshell
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, shell_present
    REAL(KIND=dp)                            :: dt, infree, kin, max_vel, &
                                                roll_tol, trvg
    REAL(KIND=dp), DIMENSION(3)              :: arg_r, arg_v, e_val, poly_r, &
                                                poly_v, scale_r, scale_v, &
                                                vector_v
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin, u, uh
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(old_variables_type), POINTER        :: old
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_baro, &
                                                thermostat_part, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(gci,thermostat_baro,thermostat_part,thermostat_shell,force_env)
    NULLIFY(atomic_kinds,cell,para_env,subsys,local_molecules,local_particles)
    NULLIFY(molecule_kinds,molecules,molecule_kind_set,npt)
    NULLIFY(core_particles,particles,shell_particles,tmp,old)
    NULLIFY(core_particle_set,particle_set,shell_particle_set)
    NULLIFY(simpar,virial)
    
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_part=thermostat_part,&
         thermostat_baro=thermostat_baro, thermostat_shell=thermostat_shell,&
         npt=npt, itimes=itimes, first_time=first_time, para_env=para_env, error=error)
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar % nfree,KIND=dp)

    CALL force_env_get(force_env, subsys=subsys,virial = virial, cell = cell ,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,&
         local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules, &
         gci=gci,molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules % n_els
    molecule_set => molecules % els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic)

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env % group  ,error=error)
    END IF

    ! Allocate work storage for positions and velocities
    CALL allocate_old ( old, particle_set, npt, error=error )

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,&
            shell_particles=shell_particles, core_particles=core_particles,&
            error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)
       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(tmp,nparticle, nshell, shell_adiabatic, error=error)

    ! Apply Thermostat to Barostat
    CALL apply_thermostat_baro( thermostat_baro, npt, para_env%group, error)

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_f_ensemble) THEN
      IF(shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, para_env%group,shell_adiabatic=shell_adiabatic,&
             shell_particle_set=shell_particle_set,core_particle_set=core_particle_set,&
             error=error)
      ELSE
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, para_env%group ,error=error)
      END IF
    END IF

    ! Apply Thermostat over the core-shell motion 
    CALL apply_thermostat_shells(thermostat_shell, atomic_kind_set, particle_set,&
         local_particles, para_env%group, shell_particle_set=shell_particle_set,&
         core_particle_set=core_particle_set, error=error)

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)
    END IF

    ! setting up for ROLL: saving old variables
    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ENDIF

    SR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! SHAKE-ROLL LOOP

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF
       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       trvg = npt ( 1, 1 ) % v  + npt ( 2, 2 ) % v + npt ( 3, 3 ) % v
       !
       ! find eigenvalues and eigenvectors of npt ( :, : ) % v
       !

       CALL diagonalise ( matrix = npt(:,:)%v, mysize = 3, &
            storageform = "UPPER", eigenvalues = e_val, eigenvectors = u )

       arg_r ( : ) =  0.5_dp * e_val ( : ) * dt * &
            0.5_dp * e_val ( : ) * dt
       poly_r = 1.0_dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
       scale_r ( : ) = EXP ( 0.5_dp * dt * e_val ( : ) )

       arg_v ( : ) =  0.25_dp * dt * ( e_val ( : ) + trvg * infree ) *&
            0.25_dp * dt * ( e_val ( : ) + trvg * infree )
       poly_v = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
       scale_v ( : ) = EXP ( -0.25_dp * dt * ( &
            e_val ( : ) + trvg * infree ) )

       CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, max_vel, poly_r, poly_v, scale_r, scale_v, u)

       IF(simpar%variable_dt) THEN
           CALL mp_max(max_vel, para_env%group)
           IF(max_vel*dt > simpar%dr_tol) THEN
              simpar%dt_fact = simpar%dr_tol/max_vel/dt
              dt = dt*simpar%dt_fact
              arg_r(:) = arg_r(:) * simpar%dt_fact*simpar%dt_fact
              poly_r = 1._dp + e2*arg_r + e4*arg_r*arg_r + e6*arg_r**3 + e8*arg_r**4
              scale_r ( : ) = EXP ( 0.5_dp * dt * e_val ( : ) )
              arg_v(:) = arg_v(:) * simpar%dt_fact*simpar%dt_fact
              poly_v = 1.0_dp + e2*arg_v + e4*arg_v*arg_v + e6*arg_v**3 + e8*arg_v**4
              scale_v ( : ) = EXP ( -0.25_dp * dt * ( &
                    e_val ( : ) + trvg * infree ) )

              CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
                  core_particle_set, shell_particle_set, nparticle_kind,&
                  shell_adiabatic, dt, max_vel, poly_v, scale_r, scale_v, u)
           ELSE
              simpar%dt_fact = 1.0_dp
           END IF
       END IF


       roll_tol = 0.0_dp
       IF ( simpar % constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt, simpar,  &
            roll_tol, iroll, scale_r*poly_r, scale_v*poly_v, para_env % group, u, cell,&
            local_particles=local_particles, error=error )
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

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE., error=error )

    ! Update forces
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second (tmp, atomic_kind_set, local_particles, particle_set,&
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, poly_v, scale_v, u )

    roll_tol = -simpar%roll_tol
    IF (simpar%constraint) THEN
       first = .TRUE.
       iroll = 1
       CALL set (old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    END IF
    RR: DO WHILE (ABS(roll_tol)>=simpar%roll_tol) ! RATTLE-ROLL LOOP
       roll_tol = 0.0_dp
       IF ( simpar % constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set,  &
            local_molecules, tmp%vel, dt, cell, npt, simpar, virial, vector_v,  &
            roll_tol, iroll, infree, first, para_env, u, error=error)

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env % group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

    END DO RR

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_f_ensemble) THEN
      IF (shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, para_env%group,shell_adiabatic=shell_adiabatic,&
             vel=tmp%vel, shell_vel=tmp%shell_vel, core_vel=tmp%core_vel, error=error)

      ELSE
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, para_env%group, vel=tmp%vel, error=error)
      END IF
    END IF

    ! Apply Thermostat over the core-shell motion 
    IF(ASSOCIATED(thermostat_shell)) THEN
       CALL apply_thermostat_shells(thermostat_shell,atomic_kind_set, particle_set,&
            local_particles, para_env%group, vel=tmp%vel, shell_vel=tmp%shell_vel, &
            core_vel=tmp%core_vel, error=error)
    END IF
         
    ! Apply Thermostat to Barostat
    CALL apply_thermostat_baro( thermostat_baro, npt, para_env%group, error)

    ! Annealing of particle velocities is only possible when no thermostat is active
    IF (simpar% ensemble == npe_f_ensemble .AND. simpar%annealing) THEN
       tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing
       IF (shell_adiabatic) THEN
          CALL shell_scale_comv(atomic_kind_set,local_particles,particle_set,&
               tmp%vel,tmp%shell_vel,tmp%core_vel)
       END IF
    END IF
    ! Annealing of CELL velocities is only possible when no thermostat is active
    IF (simpar% ensemble == npe_f_ensemble .AND. simpar%annealing_cell) THEN
       npt(:,:)%v = npt(:,:)%v * simpar%f_annealing_cell
    END IF

    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE., error=error )

    ! Update constraint virial
    IF ( simpar % constraint ) &
         CALL pv_constraint ( gci, local_molecules, molecule_set, &
         molecule_kind_set, particle_set, virial, para_env % group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env % group, error=error  )

    ! Deallocate old variables
    CALL deallocate_old ( old, error=error )

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE npt_f

! *****************************************************************************
!> \brief RESPA integrator for nve ensemble for particle positions & momenta
!> \author FS
! *****************************************************************************
  SUBROUTINE nve_respa ( md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nve_respa', &
      routineP = moduleN//':'//routineN

    INTEGER :: i_step, iparticle, iparticle_kind, iparticle_local, &
      n_time_steps, nmolecule, nmolecule_kind, nparticle, nparticle_kind, &
      nparticle_local, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dm, dt, mass
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: pos, vel
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys, subsys_respa
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_constraint_type), POINTER    :: gci
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles, particles_respa
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                particle_set_respa
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY (para_env,cell,subsys_respa,particles_respa, particle_set_respa,gci, force_env)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env, para_env=para_env, error=error)
    dt = simpar%dt

    n_time_steps=simpar%n_time_steps

    CALL force_env_get(force_env,subsys=subsys, virial=virial, cell=cell,error=error)
    CALL force_env_get(force_env%sub_force_env(1)%force_env,subsys=subsys_respa,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles, local_molecules_new=local_molecules, molecules_new=molecules,&
         gci=gci, molecule_kinds_new=molecule_kinds ,error=error)

    CALL cp_subsys_get (subsys=subsys_respa, particles=particles_respa,error=error)
    particle_set_respa => particles_respa%els

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    nmolecule_kind = molecule_kinds % n_els
    molecule_kind_set => molecule_kinds % els

    nparticle = particles%n_els
    particle_set => particles%els
    nmolecule = molecules % n_els
    molecule_set => molecules % els

    ! Allocate work storage for positions and velocities
    ALLOCATE (pos(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (vel(3,nparticle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    vel(:,:) = 0.0_dp

    IF (simpar%constraint) CALL getold(gci, local_molecules, molecule_set,  &
         molecule_kind_set, particle_set, cell, error)

    ! Multiple time step (first part)
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       dm = 0.5_dp * dt / mass
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel (:,iparticle) = particle_set(iparticle)%v(:) +&
               dm * (particle_set(iparticle)%f(:)-&
               particle_set_respa(iparticle)%f(:))
       END DO
    END DO

    ! Velocity Verlet (first part)
    DO i_step=1,n_time_steps
       pos(:,:) = 0.0_dp
       DO iparticle_kind=1,nparticle_kind
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          dm = 0.5_dp*dt/(n_time_steps*mass)
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             vel (:,iparticle) = vel(:,iparticle) +&
                  dm * particle_set_respa(iparticle)%f(:)
             pos (:,iparticle) = particle_set(iparticle)%r(:) +&
                  (dt/n_time_steps) * vel ( :, iparticle )
          END DO
       END DO

       IF (simpar%constraint) THEN 
          ! Possibly update the target values
          CALL shake_update_targets(gci, local_molecules, molecule_set, &
              molecule_kind_set, dt, force_env%root_section, error)

          CALL shake_control( gci, local_molecules, molecule_set,&
              molecule_kind_set, particle_set,pos, vel, dt, simpar%shake_tol,&
              simpar%info_constraint, simpar%lagrange_multipliers,simpar%dump_lm, cell,&
              para_env%group, local_particles, error )
       END IF

       ! Broadcast the new particle positions
       CALL update_particle_set ( particle_set, para_env % group, pos = pos )
       DO iparticle = 1, SIZE(particle_set)
          particle_set_respa(iparticle)%r = particle_set(iparticle)%r
       END DO

       ! Update forces
       CALL force_env_calc_energy_force(force_env%sub_force_env(1)%force_env,error=error)
       CALL metadyn_forces(force_env, vel, error=error)

       ! Velocity Verlet (second part)
       DO iparticle_kind=1,nparticle_kind
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
          dm = 0.5_dp*dt /(n_time_steps * mass)
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local=1,nparticle_local
             iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
             vel ( 1, iparticle) = vel ( 1, iparticle) + dm*particle_set_respa(iparticle)%f(1)
             vel ( 2, iparticle) = vel ( 2, iparticle) + dm*particle_set_respa(iparticle)%f(2)
             vel ( 3, iparticle) = vel ( 3, iparticle) + dm*particle_set_respa(iparticle)%f(3)
          END DO
       END DO

       IF (simpar%constraint) CALL rattle_control( gci,local_molecules, molecule_set, &
            molecule_kind_set, particle_set, vel, dt,simpar%shake_tol,&
            simpar%info_constraint, simpar%lagrange_multipliers,&
            simpar%dump_lm, cell, para_env%group,local_particles, error )

       IF (simpar%annealing)  vel(:,:)=vel(:,:)*simpar%f_annealing

       !     *** Update virial ***
       IF ( simpar % constraint ) CALL pv_constraint ( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, virial, para_env%group )

       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env % group , error=error )
    END DO
    DEALLOCATE (pos,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Multiple time step (second part)
    ! Compute forces for respa force_env
    CALL force_env_calc_energy_force(force_env,error=error)
    CALL metadyn_forces(force_env, vel, error=error)

    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       dm = 0.5_dp * dt  / mass
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          vel (1,iparticle) = vel(1, iparticle)+dm*(particle_set(iparticle)%f(1)-particle_set_respa(iparticle)%f(1))
          vel (2,iparticle) = vel(2, iparticle)+dm*(particle_set(iparticle)%f(2)-particle_set_respa(iparticle)%f(2))
          vel (3,iparticle) = vel(3, iparticle)+dm*(particle_set(iparticle)%f(3)-particle_set_respa(iparticle)%f(3))
       END DO
    END DO

    ! Broadcast the new particle velocities
    CALL update_particle_set ( particle_set, para_env % group, vel = vel )

    DEALLOCATE (vel,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE nve_respa

END MODULE integrator

