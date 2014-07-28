!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
!>      Lianheng Tong [LT]     12.2013 - Added regions to Langevin MD
!> \author CJM
! *****************************************************************************
MODULE integrator
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE barostat_types,                  ONLY: barostat_type
  USE cell_methods,                    ONLY: cell_type,&
                                             init_cell,&
                                             parse_cell_line
  USE constraint,                      ONLY: rattle_control,&
                                             shake_control,&
                                             shake_roll_control,&
                                             shake_update_targets
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE constraint_util,                 ONLY: getold,&
                                             pv_constraint
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_read_line
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE eigenvalueproblems,              ONLY: diagonalise
  USE extended_system_dynamics,        ONLY: shell_scale_comv
  USE extended_system_types,           ONLY: npt_info_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_set_cell,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: npe_f_ensemble,&
                                             npe_i_ensemble
  USE integrator_utils,                ONLY: &
       allocate_old, allocate_tmp, damp_v, damp_veps, deallocate_old, &
       get_s_ds, old_variables_type, rattle_roll_setup, set, &
       tmp_variables_type, update_dealloc_tmp, update_pv, update_veps, &
       variable_timestep, vv_first, vv_second
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: matmul_3x3,&
                                             transpose_3d
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE metadynamics,                    ONLY: metadyn_integrator,&
                                             metadyn_velocities_colvar
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
  USE physcon,                         ONLY: bohr,&
                                             femtoseconds
  USE qmmm_update,                     ONLY: qmmm_update_force_mixing_env
  USE qmmm_util,                       ONLY: apply_qmmm_walls_reflective
  USE reftraj_types,                   ONLY: reftraj_type
  USE reftraj_util,                    ONLY: compute_msd_reftraj
  USE shell_opt,                       ONLY: optimize_shell_core
  USE simpar_types,                    ONLY: simpar_type
  USE termination,                     ONLY: stop_program
  USE thermal_region_types,            ONLY: thermal_region_type,&
                                             thermal_regions_type
  USE thermostat_methods,              ONLY: apply_thermostat_baro,&
                                             apply_thermostat_particles,&
                                             apply_thermostat_shells
  USE thermostat_types,                ONLY: thermostat_type
  USE virial_methods,                  ONLY: virial_evaluate
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'integrator'

  PUBLIC :: isokin, langevin, nve, nvt, npt_i, npt_f, nve_respa
  PUBLIC :: nph_uniaxial_damped, nph_uniaxial, nvt_adiabatic, reftraj


CONTAINS

! *****************************************************************************
!> \brief Langevin integrator for particle positions & momenta (Brownian dynamics)
!> \param md_env ...
!> \param error ...
!> \par Literature
!>      - A. Ricci and G. Ciccotti, Mol. Phys. 101, 1927-1931 (2003)
!>      - For langevin regions:
!>        - L. Kantorovich, Phys. Rev. B 78, 094304 (2008)
!>        - L. Kantorovich and N. Rompotis, Phys. Rev. B 78, 094305 (2008)
!> \par History
!>   - Created (01.07.2005,MK)
!>   - Added support for only performing Langevin MD on a region of atoms
!>     (01.12.2013, LT)
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE langevin(md_env,error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'langevin', &
      routineP = moduleN//':'//routineN

    INTEGER :: iparticle, iparticle_kind, iparticle_local, iparticle_reg, &
      ireg, nparticle, nparticle_kind, nparticle_local, nshell, stat
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: do_langevin
    REAL(KIND=dp)                            :: c, c1, c2, c3, c4, dm, dt, &
                                                gam, mass, reg_temp, sigma
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: var_w
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
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermal_region_type), POINTER       :: thermal_region
    TYPE(thermal_regions_type), POINTER      :: thermal_regions
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY(cell, para_env, gci, force_env)
    NULLIFY(atomic_kinds, local_particles,subsys,local_molecules,molecule_kinds,molecules)
    NULLIFY(molecule_kind_set,molecule_set,particles,particle_set,rng_stream,simpar,virial)
    NULLIFY(thermal_region, thermal_regions)

    CALL get_md_env(md_env=md_env, simpar=simpar, force_env=force_env,&
         para_env=para_env, thermal_regions=thermal_regions, error=error)

    dt = simpar%dt
    gam = simpar%gamma+simpar%shadow_gamma
    nshell = 0

    CALL force_env_get(force_env=force_env,subsys=subsys,cell=cell,&
         virial=virial,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,&
                       atomic_kinds=atomic_kinds,&
                       gci=gci,&
                       local_particles=local_particles,&
                       local_molecules_new=local_molecules,&
                       molecules_new=molecules,&
                       molecule_kinds_new=molecule_kinds,&
                       nshell=nshell,&
                       particles=particles,&
                       error=error)
    CALL cp_assert((nshell == 0),cp_failure_level,cp_assertion_failed,routineP,&
                   "Langevin dynamics is not yet implemented for core-shell models",&
                   error,failure)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    ! Setup the langevin regions information
    ALLOCATE(do_langevin(nparticle), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (simpar%do_thermal_region) THEN
       DO iparticle = 1, nparticle
          do_langevin(iparticle) = thermal_regions%do_langevin(iparticle)
       END DO
    ELSE
       do_langevin(1:nparticle) = .TRUE.
    END IF

    ! Allocate the temperature dependent variance (var_w) of the
    ! random variable for each atom. It may be different for different
    ! atoms because of the possiblity of Langevin regions, and var_w
    ! for each region should depend on the temperature defined in the
    ! region
    ALLOCATE(var_w(nparticle), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    var_w(1:nparticle) = simpar%var_w
    IF (simpar%do_thermal_region) THEN
       DO ireg = 1, thermal_regions%nregions
          thermal_region => thermal_regions%thermal_region(ireg)
          DO iparticle_reg = 1, thermal_region%npart
             iparticle = thermal_region%part_index(iparticle_reg)
             reg_temp = thermal_region%temp_expected
             var_w(iparticle) = 2.0_dp*reg_temp*simpar%dt*(simpar%gamma + simpar%noisy_gamma)
          END DO
       END DO
    END IF

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
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (do_langevin(iparticle)) THEN
             sigma = var_w(iparticle)*mass
             rng_stream => local_particles%local_particle_set(iparticle_kind)%&
                   rng(iparticle_local)%stream
             w(1,iparticle) = next_random_number(rng_stream,variance=sigma,error=error)
             w(2,iparticle) = next_random_number(rng_stream,variance=sigma,error=error)
             w(3,iparticle) = next_random_number(rng_stream,variance=sigma,error=error)
          END IF
       END DO
    END DO

    DEALLOCATE(var_w,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
          IF (do_langevin(iparticle)) THEN
             vel(:,iparticle) = particle_set(iparticle)%v(:) +&
                  c3*particle_set(iparticle)%f(:)
             pos(:,iparticle) = particle_set(iparticle)%r(:) +&
                  c1*particle_set(iparticle)%v(:) +&
                  c*dm*(dt*particle_set(iparticle)%f(:) +&
                  w(:,iparticle))
          ELSE
             vel(:,iparticle) = particle_set(iparticle)%v(:) +&
                  dm*particle_set(iparticle)%f(:)
             pos(:,iparticle) = particle_set(iparticle)%r(:) +&
                  dt*particle_set(iparticle)%v(:) +&
                  dm*dt*particle_set(iparticle)%f(:)
          END IF
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
    CALL update_particle_set(particle_set,para_env%group,pos=pos,error=error)

    DEALLOCATE(pos,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Update forces
    CALL force_env_calc_energy_force(force_env,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, vel, error=error)

    ! Update Verlet (second part)
    DO iparticle_kind=1,nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
       dm = 0.5_dp*dt/mass
       c3 = dm/c2
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          IF (do_langevin(iparticle)) THEN
             vel(1,iparticle) = vel(1,iparticle) + c3*particle_set(iparticle)%f(1)
             vel(2,iparticle) = vel(2,iparticle) + c3*particle_set(iparticle)%f(2)
             vel(3,iparticle) = vel(3,iparticle) + c3*particle_set(iparticle)%f(3)
             vel(1,iparticle) = c4*vel(1,iparticle) + c2*w(1,iparticle)/mass
             vel(2,iparticle) = c4*vel(2,iparticle) + c2*w(2,iparticle)/mass
             vel(3,iparticle) = c4*vel(3,iparticle) + c2*w(3,iparticle)/mass
          ELSE
             vel(1,iparticle) = vel(1,iparticle) + dm*particle_set(iparticle)%f(1)
             vel(2,iparticle) = vel(2,iparticle) + dm*particle_set(iparticle)%f(2)
             vel(3,iparticle) = vel(3,iparticle) + dm*particle_set(iparticle)%f(3)
          END IF
       END DO
    END DO

    IF (simpar%constraint) THEN
       CALL rattle_control(gci,local_molecules,molecule_set,molecule_kind_set,&
            particle_set,vel,dt,simpar%shake_tol,&
            simpar%info_constraint,simpar%lagrange_multipliers,&
            simpar%dump_lm,cell,para_env%group,local_particles, error)
    END IF

    ! Broadcast the new particle velocities
    CALL update_particle_set(particle_set,para_env%group,vel=vel,error=error)

    DEALLOCATE(vel,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(w,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DEALLOCATE(do_langevin,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Update virial
    IF (simpar%constraint) CALL pv_constraint(gci,local_molecules,molecule_set,&
       molecule_kind_set, particle_set,virial,para_env%group)

    CALL virial_evaluate(atomic_kind_set,particle_set,local_particles,&
         virial,para_env%group, error=error)

  END SUBROUTINE langevin

! *****************************************************************************
!> \brief nve integrator for particle positions & momenta
!> \param md_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>   - the local particle lists are used instead of pnode (Sep. 2003,MK)
!>   - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!> \author CJM
! *****************************************************************************
  SUBROUTINE nve ( md_env, globenv, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i_iter, n_iter, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL :: deallocate_vel, ehrenfest_md, failure, shell_adiabatic, &
      shell_check_distance, shell_present
    REAL(KIND=dp)                            :: dt
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: v_old
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
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_coeff, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY (thermostat_coeff, tmp)
    NULLIFY (subsys, simpar, para_env, cell, gci, force_env, virial)
    NULLIFY (atomic_kinds, local_particles, molecules, molecule_kind_set, molecule_set, particle_set)
    NULLIFY (shell_particles, shell_particle_set, core_particles, &
            core_particle_set, thermostat_shell)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_coeff=thermostat_coeff,thermostat_shell=thermostat_shell,&
         para_env=para_env, ehrenfest_md=ehrenfest_md, error=error)
    dt = simpar%dt
    CALL force_env_get(force_env=force_env, subsys=subsys,&
         virial=virial, cell=cell, error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules, &
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic,&
         shell_check_distance=shell_check_distance)

    IF(shell_present) THEN
       CALL cp_subsys_get(subsys=subsys,shell_particles=shell_particles,&
            core_particles=core_particles, error=error)
       shell_particle_set =>  shell_particles%els
       nshell = SIZE(shell_particles%els)

       IF(shell_adiabatic) THEN
          core_particle_set =>  core_particles%els
       END IF
    END IF

    CALL allocate_tmp(md_env, tmp, nparticle,nshell,shell_adiabatic,error=error)
    IF(ehrenfest_md) ALLOCATE(v_old(3,SIZE(tmp%vel,2)))

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
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error)

    IF(simpar%variable_dt) CALL variable_timestep(md_env,tmp, dt, simpar, para_env, atomic_kind_set,&
         local_particles, particle_set, core_particle_set, shell_particle_set,&
         nparticle_kind, shell_adiabatic, error=error)

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

    IF(shell_adiabatic .AND. shell_check_distance) THEN
        CALL optimize_shell_core(force_env,  particle_set,&
             shell_particle_set, core_particle_set, globenv, tmp=tmp, check=.TRUE., error=error)
    END IF

    ! Update forces
    ! In case of ehrenfest dynamics, velocities need to be iterated
    IF(ehrenfest_md)THEN
       v_old(:,:)=tmp%vel
       CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error )
       CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
            core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,&
            should_deall_vel=.FALSE.,error=error )
       tmp%vel=v_old
       n_iter=force_env%qs_env%dft_control%rtp_control%max_iter
    ELSE
       n_iter=1
    END IF
    DO i_iter=1,n_iter

       deallocate_vel=.TRUE.
       IF(ehrenfest_md)THEN
          force_env%qs_env%rtp%iter=i_iter
          tmp%vel=v_old
          force_env%qs_env%rtp%converged=.FALSE.
          deallocate_vel=.FALSE.
       END IF

       ![NB] let nve work with force mixing which does not have consistent energies and forces
       CALL force_env_calc_energy_force(force_env,require_consistent_energy_force=.FALSE.,error=error)

       ! Metadynamics
       CALL metadyn_integrator(force_env, tmp%vel, error=error)

       ! Velocity Verlet (second part)
       CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error )

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

       IF(ehrenfest_md)deallocate_vel=force_env%qs_env%rtp%converged
       IF(i_iter.EQ.n_iter) deallocate_vel=.TRUE.
       ! Broadcast the new particle velocities and deallocate the full temporary
       CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
            core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,&
            should_deall_vel=deallocate_vel,error=error )
       IF(ehrenfest_md)THEN
          IF(force_env%qs_env%rtp%converged)EXIT
          !  will be handled better in the next step of the implementation
          IF(i_iter==n_iter)&
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                 routineP,"EMD did not converge, either increase MAX_ITER "//&
                 "or use a smaller TIMESTEP",&
                 error,failure)
       END IF

    END DO

    ! Update virial
    IF ( simpar%constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial,para_env%group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group  ,error=error)

  END SUBROUTINE nve

! *****************************************************************************
!> \brief simplest version of the isokinetic gaussian thermostat
!> \param md_env ...
!> \param error ...
!> \par History
!>   - Created [2004-07]
!> \author Joost VandeVondele
!> \note
!>      - time reversible, and conserves the kinetic energy to machine precision
!>      - is not yet supposed to work for e.g. constraints, our the extended version
!>        of this thermostat
!>        see:
!>         - Zhang F. , JCP 106, 6102 (1997)
!>         - Minary P. et al, JCP 118, 2510 (2003)
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
    REAL(KIND=dp)                            :: dt
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
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(tmp_variables_type), POINTER        :: tmp

    failure = .FALSE.

    NULLIFY(force_env, tmp, simpar)
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
                         "Constraints not yet implemented",para_env)
    END IF

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

    CALL allocate_tmp (md_env,tmp, nparticle, nshell, shell_adiabatic, error=error)

    ! compute s,ds
    CALL get_s_ds(tmp, nparticle_kind, atomic_kind_set, local_particles, particle_set,&
         dt, para_env)

    ! Velocity Verlet (first part)
    tmp%scale_v(1:3) = SQRT(1.0_dp/tmp%ds)
    tmp%poly_v(1:3) = 2.0_dp*tmp%s/SQRT(tmp%ds)/dt
    CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, error=error)


    IF(simpar%variable_dt) CALL variable_timestep(md_env,tmp, dt, simpar,para_env,atomic_kind_set, &
      local_particles, particle_set, core_particle_set, shell_particle_set,&
      nparticle_kind, shell_adiabatic,error=error)

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    CALL force_env_calc_energy_force(force_env,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, error=error)

    ! compute s,ds
    CALL get_s_ds(tmp, nparticle_kind, atomic_kind_set,  local_particles, particle_set,&
       dt, para_env, tmpv=.TRUE.)

    ! Velocity Verlet (second part)
    tmp%scale_v(1:3) = SQRT(1.0_dp/tmp%ds)
    tmp%poly_v(1:3) = 2.0_dp*tmp%s/SQRT(tmp%ds)/dt
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, error=error)

    IF (simpar%annealing)  tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing

    !  Broadcast the new particle velocities and deallocate the temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

  END SUBROUTINE isokin
! *****************************************************************************
!> \brief nvt adiabatic integrator for particle positions & momenta
!> \param md_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>   - the local particle lists are used instead of pnode (Sep. 2003,MK)
!>   - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!> \author CJM
! *****************************************************************************
  SUBROUTINE nvt_adiabatic ( md_env, globenv, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nvt_adiabatic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, ivar, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, shell_adiabatic, &
                                                shell_check_distance, &
                                                shell_present
    REAL(KIND=dp)                            :: dt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rand
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
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_coeff, &
                                                thermostat_fast, &
                                                thermostat_shell, &
                                                thermostat_slow
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.


    NULLIFY(gci, force_env, thermostat_coeff, tmp, &
         thermostat_fast, thermostat_slow, thermostat_shell, cell, shell_particles, &
         shell_particle_set, core_particles, core_particle_set, rand)
    NULLIFY(para_env, subsys, local_molecules, local_particles, molecule_kinds,&
         molecules, molecule_kind_set, molecule_set, atomic_kinds,particles)
    NULLIFY(simpar)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_fast=thermostat_fast, thermostat_slow=thermostat_slow,&
         thermostat_coeff=thermostat_coeff, thermostat_shell=thermostat_shell, &
         para_env=para_env, error=error)
    dt = simpar%dt

    CALL force_env_get(force_env=force_env, subsys=subsys, &
         virial=virial, cell=cell,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind  = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle    = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic,&
         shell_check_distance=shell_check_distance)

    IF (ASSOCIATED(force_env%meta_env)) THEN
       ! Allocate random number for Langevin Thermostat acting on COLVARS
       IF(force_env%meta_env%langevin) THEN
          ALLOCATE (rand(force_env%meta_env%n_colvar),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          rand(:) = 0.0_dp
       ENDIF
    ENDIF

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

    CALL allocate_tmp(md_env, tmp, nparticle,nshell,shell_adiabatic,error=error)

    ! Apply Thermostat over the full set of particles
    IF(shell_adiabatic) THEN
!       CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
!            particle_set, local_molecules, para_env%group, shell_adiabatic=shell_adiabatic,&
!            shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
!            error=error)

       CALL apply_thermostat_shells(thermostat_shell, atomic_kind_set, particle_set,&
            local_particles, para_env%group, shell_particle_set=shell_particle_set,&
            core_particle_set=core_particle_set, error=error)
    ELSE
       CALL apply_thermostat_particles(thermostat_fast, molecule_kind_set, molecule_set,&
            particle_set, local_molecules, local_particles, para_env%group ,error=error)

       CALL apply_thermostat_particles(thermostat_slow, molecule_kind_set, molecule_set,&
            particle_set, local_molecules, local_particles, para_env%group ,error=error)
    END IF

    IF (simpar%constraint) CALL getold( gci, local_molecules, molecule_set,  &
         molecule_kind_set, particle_set, cell, error)

    !    *** Velocity Verlet for Langeving *** v(t)--> v(t+1/2)
    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          DO ivar = 1 , force_env%meta_env%n_colvar
             rng_stream => force_env%meta_env%rng(ivar)%stream
             rand(ivar)=next_random_number(rng_stream,error=error)
          ENDDO
          CALL metadyn_velocities_colvar(force_env,rand,error=error)
       ENDIF
    ENDIF

   ! Velocity Verlet (first part)
    CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error)

    IF(simpar%variable_dt) CALL variable_timestep(md_env,tmp, dt, simpar,para_env,atomic_kind_set,&
      local_particles, particle_set, core_particle_set, shell_particle_set,&
      nparticle_kind, shell_adiabatic,error=error)

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

    IF(shell_adiabatic .AND. shell_check_distance) THEN
        CALL optimize_shell_core(force_env,  particle_set,&
             shell_particle_set, core_particle_set, globenv, tmp=tmp, check=.TRUE., error=error)
    END IF

    ! Update forces
    CALL force_env_calc_energy_force(force_env,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, rand=rand, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error)

    IF (simpar%constraint) CALL rattle_control( gci, local_molecules, molecule_set, &
         molecule_kind_set, particle_set, tmp%vel, dt,simpar%shake_tol,&
         simpar%info_constraint,simpar%lagrange_multipliers,simpar%dump_lm,&
         cell,para_env%group,local_particles, error )

    ! Apply Thermostat over the full set of particles
    IF (shell_adiabatic) THEN
     !  CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
     !       particle_set, local_molecules, para_env%group, shell_adiabatic=shell_adiabatic,&
     !       vel= tmp%vel, shell_vel= tmp%shell_vel, core_vel= tmp%core_vel, error=error)

       CALL apply_thermostat_shells(thermostat_shell,atomic_kind_set, particle_set,&
            local_particles, para_env%group, vel= tmp%vel, shell_vel= tmp%shell_vel, &
            core_vel= tmp%core_vel, error=error)
    ELSE
       CALL apply_thermostat_particles(thermostat_slow,molecule_kind_set, molecule_set, &
            particle_set, local_molecules, local_particles, para_env%group, vel= tmp%vel, error=error)

       CALL apply_thermostat_particles(thermostat_fast,molecule_kind_set, molecule_set, &
            particle_set, local_molecules, local_particles, para_env%group, vel= tmp%vel, error=error)
    END IF

    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          DEALLOCATE(rand,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF

    ! Update constraint virial
    IF ( simpar%constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial,para_env%group )

    !     **  Evaluate Virial
    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group, error=error)

  END SUBROUTINE nvt_adiabatic


! *****************************************************************************
!> \brief nvt integrator for particle positions & momenta
!> \param md_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>   - the local particle lists are used instead of pnode (Sep. 2003,MK)
!>   - usage of fragments retrieved from the force environment (Oct. 2003,MK)
!> \author CJM
! *****************************************************************************
  SUBROUTINE nvt ( md_env, globenv, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nvt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, ivar, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, shell_adiabatic, &
                                                shell_check_distance, &
                                                shell_present
    REAL(KIND=dp)                            :: dt
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rand
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
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_coeff, &
                                                thermostat_part, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.

    NULLIFY(gci, force_env, thermostat_coeff, tmp, &
         thermostat_part, thermostat_shell, cell, shell_particles, &
         shell_particle_set, core_particles, core_particle_set, rand)
    NULLIFY(para_env, subsys, local_molecules, local_particles, molecule_kinds,&
         molecules, molecule_kind_set, molecule_set, atomic_kinds,particles)
    NULLIFY(simpar, thermostat_coeff, thermostat_part, thermostat_shell)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_part=thermostat_part, thermostat_coeff=thermostat_coeff,&
         thermostat_shell=thermostat_shell,para_env=para_env,&
         error=error)
    dt = simpar%dt

    CALL force_env_get(force_env=force_env, subsys=subsys, &
         virial=virial, cell=cell,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,&
         molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind  = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle    = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic,&
         shell_check_distance=shell_check_distance)

    IF (ASSOCIATED(force_env%meta_env)) THEN
       ! Allocate random number for Langevin Thermostat acting on COLVARS
       IF(force_env%meta_env%langevin) THEN
          ALLOCATE (rand(force_env%meta_env%n_colvar),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          rand(:) = 0.0_dp
       ENDIF
    ENDIF

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

    CALL allocate_tmp(md_env, tmp, nparticle,nshell,shell_adiabatic,error=error)

    ! Apply Thermostat over the full set of particles
    IF(shell_adiabatic) THEN
       CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
            particle_set, local_molecules, local_particles, para_env%group, shell_adiabatic=shell_adiabatic,&
            shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
            error=error)

       CALL apply_thermostat_shells(thermostat_shell, atomic_kind_set, particle_set,&
            local_particles, para_env%group, shell_particle_set=shell_particle_set,&
            core_particle_set=core_particle_set, error=error)
    ELSE
       CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
            particle_set, local_molecules, local_particles, para_env%group ,error=error)
    END IF

    IF (simpar%constraint) CALL getold( gci, local_molecules, molecule_set,  &
         molecule_kind_set, particle_set, cell, error)

    !    *** Velocity Verlet for Langeving *** v(t)--> v(t+1/2)
    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          DO ivar = 1 , force_env%meta_env%n_colvar
             rng_stream => force_env%meta_env%rng(ivar)%stream
             rand(ivar)=next_random_number(rng_stream,error=error)
          ENDDO
          CALL metadyn_velocities_colvar(force_env,rand,error=error)
       ENDIF
    ENDIF

   ! Velocity Verlet (first part)
    CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error)

    IF(simpar%variable_dt) CALL variable_timestep(md_env,tmp, dt, simpar,para_env,atomic_kind_set,&
      local_particles, particle_set, core_particle_set, shell_particle_set,&
      nparticle_kind, shell_adiabatic,error=error)

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

    IF(shell_adiabatic .AND. shell_check_distance) THEN
        CALL optimize_shell_core(force_env,  particle_set,&
             shell_particle_set, core_particle_set, globenv, tmp=tmp, check=.TRUE., error=error)
    END IF

    ![ADAPT] update input structure with new coordinates, make new labels
    CALL qmmm_update_force_mixing_env(force_env, force_env%root_section, error=error)

    ![NB] recreate pointers changed by creation of new subsys in qmmm_update_force_mixing_env
    ![NB] ugly hack, which is why adaptivity isn't implemented in most other ensembles
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL force_env_get(force_env=force_env, subsys=subsys, &
       virial=virial, cell=cell,error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
       particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,&
       molecule_kinds_new=molecule_kinds,gci=gci,error=error)

    nparticle_kind  = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle    = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
       shell_present=shell_present,shell_adiabatic=shell_adiabatic,&
       shell_check_distance=shell_check_distance)

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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Update forces
    ![NB] let nvt work with force mixing which does not have consistent energies and forces
    CALL force_env_calc_energy_force(force_env,require_consistent_energy_force=.FALSE.,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, rand=rand, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind, shell_adiabatic, dt, error=error)

    IF (simpar%constraint) CALL rattle_control( gci, local_molecules, molecule_set, &
         molecule_kind_set, particle_set, tmp%vel, dt,simpar%shake_tol,&
         simpar%info_constraint,simpar%lagrange_multipliers,simpar%dump_lm,&
         cell,para_env%group,local_particles, error )

    ! Apply Thermostat over the full set of particles
    IF (shell_adiabatic) THEN
       CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
            particle_set, local_molecules, local_particles, para_env%group, shell_adiabatic=shell_adiabatic,&
            vel= tmp%vel, shell_vel= tmp%shell_vel, core_vel= tmp%core_vel, error=error)

       CALL apply_thermostat_shells(thermostat_shell,atomic_kind_set, particle_set,&
            local_particles, para_env%group, vel= tmp%vel, shell_vel= tmp%shell_vel, &
            core_vel= tmp%core_vel, error=error)
    ELSE
       CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
            particle_set, local_molecules, local_particles, para_env%group, vel= tmp%vel, error=error)
    END IF

    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          DEALLOCATE(rand,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF

    ! Update constraint virial
    IF ( simpar%constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial,para_env%group )

    !     **  Evaluate Virial
    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group, error=error)

  END SUBROUTINE nvt

! *****************************************************************************
!> \brief npt_i integrator for particle positions & momenta
!>      isotropic box changes
!> \param md_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE npt_i ( md_env, globenv, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'npt_i', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: box_change = .TRUE.
    REAL(KIND=dp), PARAMETER                 :: e2 = 1.0_dp/6.0_dp, &
                                                e4 = e2/20.0_dp, &
                                                e6 = e4/42.0_dp, &
                                                e8 = e6/72.0_dp

    INTEGER                                  :: iroll, istat, ivar, &
                                                nparticle, nparticle_kind, &
                                                nshell
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, &
                                                shell_check_distance, &
                                                shell_present
    REAL(KIND=dp)                            :: dt, infree, kin, roll_tol, &
                                                roll_tol_thrs
    REAL(KIND=dp), DIMENSION(3)              :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rand
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
    TYPE(rng_stream_type), POINTER           :: rng_stream
    TYPE(simpar_type), POINTER               :: simpar
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
    NULLIFY(simpar,virial,rand)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_part=thermostat_part, thermostat_baro=thermostat_baro,&
         thermostat_shell=thermostat_shell, npt=npt, first_time=first_time,&
         para_env=para_env, error=error)
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar%nfree,KIND=dp)

    CALL force_env_get(force_env=force_env, subsys=subsys, virial=virial,&
         cell=cell, error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules, &
         gci=gci,molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic,&
         shell_check_distance=shell_check_distance)

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env%group, error=error)
    END IF

    ! Allocate work storage for positions and velocities
    CALL allocate_old ( old, particle_set, npt, error=error )

    IF (ASSOCIATED(force_env%meta_env)) THEN
       ! Allocate random number for Langevin Thermostat acting on COLVARS
       IF(force_env%meta_env%langevin) THEN
          ALLOCATE (rand(force_env%meta_env%n_colvar),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          rand(:) = 0.0_dp
       ENDIF
    ENDIF

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

    CALL allocate_tmp(md_env, tmp, nparticle, nshell, shell_adiabatic, error=error)

    ! Initialize eps_0 the first time through
    IF ( first_time ) eps_0 = npt (1,1)%eps

    ! Apply thermostat to  barostat
    CALL apply_thermostat_baro( thermostat_baro, npt, para_env%group, error)

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_i_ensemble) THEN
      IF(shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, local_particles, para_env%group,shell_adiabatic=shell_adiabatic,&
             shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
             error=error)

      ELSE
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, local_particles, para_env%group ,error=error)
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
    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    !    *** Velocity Verlet for Langeving *** v(t)--> v(t+1/2)
    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          DO ivar = 1 , force_env%meta_env%n_colvar
             rng_stream => force_env%meta_env%rng(ivar)%stream
             rand(ivar)=next_random_number(rng_stream,error=error)
          ENDDO
          CALL metadyn_velocities_colvar(force_env,rand,error=error)
       ENDIF
    ENDIF

    SR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! SHAKE-ROLL LOOP

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF

       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       tmp%arg_r(1) = ( 0.5_dp * npt ( 1, 1 )%v * dt ) * &
            ( 0.5_dp * npt ( 1, 1 )%v * dt )
       tmp%poly_r(1:3) = 1.0_dp + e2*tmp%arg_r(1) + e4*tmp%arg_r(1)*tmp%arg_r(1) + &
            e6*tmp%arg_r(1)**3 + e8*tmp%arg_r(1)**4

       tmp%arg_v(1) = ( 0.25_dp * npt ( 1, 1 )%v * dt *  &
            (1.0_dp + 3.0_dp * infree ) )*( 0.25_dp * npt ( 1, 1 )%v * &
            dt * ( 1.0_dp + 3.0_dp * infree ) )
       tmp%poly_v(1:3) = 1.0_dp + e2*tmp%arg_v(1) + e4*tmp%arg_v(1)*tmp%arg_v(1) + &
            e6*tmp%arg_v(1)**3 + e8*tmp%arg_v(1)**4

       tmp%scale_r(1:3) = EXP( 0.5_dp * dt * npt ( 1, 1 )%v )
       tmp%scale_v(1:3) = EXP( -0.25_dp * dt * npt ( 1, 1 )%v * &
            (1.0_dp + 3.0_dp * infree ) )

       ! first half of velocity verlet
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, error=error )

       IF(simpar%variable_dt) CALL  variable_timestep(md_env,tmp, dt, simpar,para_env,&
            atomic_kind_set, local_particles, particle_set, core_particle_set,&
            shell_particle_set, nparticle_kind, shell_adiabatic, npt=npt,error=error)


       roll_tol = 0.0_dp
       vector_r ( : ) = tmp%scale_r(:) * tmp%poly_r(:)
       vector_v ( : ) = tmp%scale_v(:) * tmp%poly_v(:)

       IF ( simpar%constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt, simpar,  &
            roll_tol, iroll, vector_r, vector_v, para_env%group, cell=cell, &
            local_particles=local_particles, error=error )
    END DO SR

    ! Update eps:
    npt ( :, : )%eps =  npt ( :, : )%eps + dt * npt ( :, : )%v

    ! Update h_mat
    cell%hmat ( :, : ) = cell%hmat ( :, : ) * EXP ( npt ( 1, 1 )%eps - eps_0 )

    eps_0 =  npt ( 1, 1 )%eps

    ! Update the inverse
    CALL init_cell ( cell )

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp ( tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    IF(shell_adiabatic .AND. shell_check_distance) THEN
        CALL optimize_shell_core(force_env,  particle_set,&
             shell_particle_set, core_particle_set, globenv, tmp=tmp, check=.TRUE., error=error)
    END IF

    ! Update forces
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, rand=rand, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, error=error )

    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       first = .TRUE.
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    RR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! RATTLE-ROLL LOOP
       roll_tol = 0.0_dp
       IF ( simpar%constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set,  &
            local_molecules, tmp%vel, dt,  cell, npt, simpar, virial, vector_v,  &
            roll_tol, iroll, infree, first, para_env,error=error)

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )
    END DO RR

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_i_ensemble) THEN
      IF (shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, local_particles, para_env%group, shell_adiabatic=shell_adiabatic,&
             vel=tmp%vel, shell_vel=tmp%shell_vel, core_vel=tmp%core_vel, error=error)
      ELSE
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, local_particles, para_env%group, vel=tmp%vel, error=error)
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
    IF ( simpar%constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set, molecule_kind_set, particle_set, virial, para_env%group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group , error=error)

    ! Deallocate old variables
    CALL deallocate_old (old, error=error )

    IF (ASSOCIATED(force_env%meta_env)) THEN
       IF(force_env%meta_env%langevin) THEN
          DEALLOCATE(rand,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE npt_i

! *****************************************************************************
!> \brief uses coordinates in a file and generates frame after frame of these
!> \param md_env ...
!> \param error ...
!> \par History
!>   - 04.2005 created [Joost VandeVondele]
!>   - modified to make it more general [MI]
!> \note
!>     it can be used to compute some properties on already available trajectories
! *****************************************************************************
  SUBROUTINE reftraj ( md_env, error )
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'reftraj', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: AA
    INTEGER                                  :: cell_itimes, I, nparticle, &
                                                Nread, trj_itimes
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, init, my_end, test_ok
    REAL(KIND=dp)                            :: cell_time, h(3,3), trj_epot, &
                                                trj_time, vol
    REAL(KIND=dp), POINTER                   :: time
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(reftraj_type), POINTER              :: reftraj_env
    TYPE(simpar_type), POINTER               :: simpar

    failure = .FALSE.
    NULLIFY(reftraj_env, particle_set, particles, force_env, subsys, simpar, para_env, cell, itimes, time)
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

    ! SnapShots read from an external file (parsers calls are buffered! please
    ! don't put any additional MPI call!) [tlaino]
    CALL parser_read_line(reftraj_env%info%traj_parser,1,error=error)
    READ(reftraj_env%info%traj_parser%input_line,FMT="(I8)") nread
    CPPostcondition(nread==nparticle,cp_failure_level,routineP,error,failure)
    CALL parser_read_line(reftraj_env%info%traj_parser,1,error=error)
    test_ok = .FALSE.
    READ(reftraj_env%info%traj_parser%input_line,FMT="(T6,I8,T23,F12.3,T41,F20.10)",ERR=999)&
         trj_itimes, trj_time, trj_epot
    test_ok = .TRUE.
999 IF (.NOT.test_ok) THEN
       ! Handling properly the error when reading the title of an XYZ
       CALL get_md_env(md_env, itimes=itimes, error=error)
       trj_itimes = itimes
       trj_time   = 0.0_dp
       trj_epot   = 0.0_dp
    END IF
    DO i = 1,nread-1
       CALL parser_read_line(reftraj_env%info%traj_parser,1,error=error)
       READ(reftraj_env%info%traj_parser%input_line(1:LEN_TRIM(reftraj_env%info%traj_parser%input_line)),*) AA,particle_set(i)%r
       particle_set(i)%r = particle_set(i)%r*bohr
    END DO
    ! End of file is properly addressed in the previous call..
    ! Let's check directly (providing some info) also for the last
    ! line of this frame..
    CALL parser_read_line(reftraj_env%info%traj_parser,1,at_end=my_end,error=error)
    READ(unit=reftraj_env%info%traj_parser%input_line,fmt=*) AA,particle_set(i)%r
    particle_set(i)%r(1) = cp_unit_to_cp2k(particle_set(i)%r(1),"angstrom",error=error)
    particle_set(i)%r(2) = cp_unit_to_cp2k(particle_set(i)%r(2),"angstrom",error=error)
    particle_set(i)%r(3) = cp_unit_to_cp2k(particle_set(i)%r(3),"angstrom",error=error)

    ! Check if we reached the end of the file and provide some info..
    IF (my_end) THEN
       CALL cp_assert(reftraj_env%isnap==(simpar%nsteps-1),cp_fatal_level,cp_assertion_failed,routineP,&
            "Reached the end of the Trajectory  frames in the TRAJECTORY file. Number of "//&
            "missing frames ("//cp_to_string((simpar%nsteps-1)-reftraj_env%isnap)//")."//&
CPSourceFileRef,&
            only_ionode=.TRUE.)
    END IF

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
    CALL get_md_env(md_env, t=time, error=error)
    time = reftraj_env%time


    IF(reftraj_env%info%variable_volume)THEN
       cell%hmat = h
       CALL init_cell(cell)
       CALL force_env_set_cell(force_env, cell, error)
    END IF

    ![ADAPT] update input structure with new coordinates, make new labels
    CALL qmmm_update_force_mixing_env(force_env, force_env%root_section, error=error)
    ! no pointers into force_env%subsys to update

    ! Task to perform on the reference trajectory
    ! Compute energy and forces
    ![NB] let reftraj work with force mixing which does not have consistent energies and forces
    CALL force_env_calc_energy_force(force_env,eval_energy_forces=reftraj_env%info%eval_EF,&
         require_consistent_energy_force=.FALSE., error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, error=error)

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
!> \param md_env ...
!> \param error ...
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

    INTEGER                                  :: iroll, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, shell_present
    REAL(KIND=dp)                            :: dt, infree, kin, roll_tol, &
                                                roll_tol_thrs
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
    TYPE(simpar_type), POINTER               :: simpar
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
    infree = 1.0_dp / REAL ( simpar%nfree, dp )

    CALL force_env_get(force_env, subsys=subsys, virial = virial, cell=cell ,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules, molecules_new=molecules, gci=gci,&
         molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env%group  ,error=error)
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

    CALL allocate_tmp(md_env, tmp,nparticle, nshell, shell_adiabatic, error=error)

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)
    END IF

    ! setting up for ROLL: saving old variables
    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    SR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! SHAKE-ROLL LOOP

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF
       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       tmp%arg_r(1) = ( 0.5_dp * npt ( 1, 1 )%v * dt ) * &
            ( 0.5_dp * npt ( 1, 1 )%v * dt )
       tmp%poly_r(1) = 1._dp + e2*tmp%arg_r(1) + e4*tmp%arg_r(1)*tmp%arg_r(1) + &
           e6*tmp%arg_r(1)**3 + e8*tmp%arg_r(1)**4
       tmp%poly_r(2) = 1.0_dp
       tmp%poly_r(3) = 1.0_dp

       tmp%arg_v(1) = ( 0.25_dp * npt ( 1, 1 )%v * dt *  &
            (1._dp + infree ) )*( 0.25_dp * npt ( 1, 1 )%v * &
            dt * ( 1._dp + infree ) )
       tmp%arg_v(2) = ( 0.25_dp * npt ( 1, 1 )%v * dt * infree) * &
            ( 0.25_dp * npt ( 1, 1 )%v * dt * infree )
       tmp%poly_v(1) = 1._dp + e2*tmp%arg_v(1) + e4*tmp%arg_v(1)*tmp%arg_v(1) + &
            e6*tmp%arg_v(1)**3 + e8*tmp%arg_v(1)**4
       tmp%poly_v(2) = 1._dp + e2*tmp%arg_v(2) + e4*tmp%arg_v(2)*tmp%arg_v(2) + &
            e6*tmp%arg_v(2)**3 + e8*tmp%arg_v(2)**4
       tmp%poly_v(3) = 1._dp + e2*tmp%arg_v(2) + e4*tmp%arg_v(2)*tmp%arg_v(2) + &
            e6*tmp%arg_v(2)**3 + e8*tmp%arg_v(2)**4

       tmp%scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 )%v )
       tmp%scale_r(2) = 1.0_dp
       tmp%scale_r(3) = 1.0_dp

       tmp%scale_v(1) = EXP( -0.25_dp * dt * npt ( 1, 1 )%v * &
            (1._dp + infree ) )
       tmp%scale_v(2) = EXP( -0.25_dp * dt * npt ( 1, 1 )%v * infree )
       tmp%scale_v(3) = EXP( -0.25_dp * dt * npt ( 1, 1 )%v * infree )

       ! first half of velocity verlet
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, error=error )

       IF(simpar%variable_dt) CALL variable_timestep(md_env,tmp, dt, simpar,para_env,&
          atomic_kind_set, local_particles, particle_set, core_particle_set,&
          shell_particle_set, nparticle_kind, shell_adiabatic, npt=npt, error=error)


       roll_tol = 0._dp
       vector_r ( : ) = 0._dp
       vector_v ( : ) = tmp%scale_v(:) * tmp%poly_v(:)
       vector_r ( 1 ) = tmp%scale_r(1) * tmp%poly_r(1)

       IF ( simpar%constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt, simpar,&
            roll_tol, iroll, vector_r, vector_v, para_env%group, cell=cell, &
            local_particles=local_particles, error=error )
    END DO SR

    ! Update h_mat
    cell%hmat ( 1, 1 ) = cell%hmat ( 1, 1 ) * tmp%scale_r(1) * tmp%scale_r(1)

    ! Update the cell
    CALL init_cell ( cell )

    ! Broadcast the new particle positions and deallocate the pos component of temporary
    CALL update_dealloc_tmp (tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! Update forces (and stress)
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env, error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
     CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
          core_particle_set, shell_particle_set, nparticle_kind,&
          shell_adiabatic, dt, error=error )

    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       first = .TRUE.
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    RR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! RATTLE-ROLL LOOP
       roll_tol = 0._dp
       IF ( simpar%constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set,  &
            local_molecules, tmp%vel, dt, cell, npt, simpar, virial, vector_v,  &
            roll_tol, iroll, infree, first, para_env, error=error )

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )
    END DO RR

    IF (simpar%annealing)  tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing

    ! Broadcast the new particle velocities and deallocate the temporary
    CALL update_dealloc_tmp (tmp, particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    ! Update constraint virial
    IF ( simpar%constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set,molecule_kind_set,particle_set, virial, para_env%group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group, error=error)

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
!> \param md_env ...
!> \param error ...
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

    INTEGER                                  :: iroll, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, shell_present
    REAL(KIND=dp)                            :: aa, aax, dt, gamma1, infree, &
                                                kin, roll_tol, roll_tol_thrs
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
    TYPE(simpar_type), POINTER               :: simpar
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
    infree = 1.0_dp / REAL ( simpar%nfree, dp )
    gamma1 = simpar%gamma_nph

    CALL force_env_get(force_env, subsys=subsys, virial = virial, cell = cell ,error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds,local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules,gci=gci,&
         molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env%group , error=error)
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

    CALL allocate_tmp(md_env, tmp,nparticle, nshell, shell_adiabatic, error=error)

    ! perform damping on velocities
    CALL damp_v ( molecule_kind_set, molecule_set, particle_set, local_molecules,&
         gamma1, npt(1,1), dt, para_env%group )

    IF (simpar%constraint) THEN
       ! Possibly update the target values
       CALL shake_update_targets(gci, local_molecules, molecule_set, &
            molecule_kind_set, dt, force_env%root_section, error)
    END IF

    ! setting up for ROLL: saving old variables
    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    SR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! SHAKE-ROLL LOOP

       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF
       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

       tmp%arg_r(1) = ( 0.5_dp * npt ( 1, 1 )%v * dt ) * &
            ( 0.5_dp * npt ( 1, 1 )%v * dt )
       tmp%poly_r(1) = 1._dp + e2*tmp%arg_r(1) + e4*tmp%arg_r(1)*tmp%arg_r(1) +&
           e6*tmp%arg_r(1)**3 + e8*tmp%arg_r(1)**4

       aax = npt ( 1, 1 )%v * ( 1.0_dp + infree )
       tmp%arg_v(1) = ( 0.25_dp * dt * aax ) * ( 0.25_dp * dt * aax )
       tmp%poly_v(1) = 1._dp + e2*tmp%arg_v(1) + e4*tmp%arg_v(1)*tmp%arg_v(1) +&
          e6*tmp%arg_v(1)**3 + e8*tmp%arg_v(1)**4

       aa = npt ( 1, 1 )%v *  infree
       tmp%arg_v(2) = ( 0.25_dp * dt * aa ) * ( 0.25_dp * dt * aa )
       tmp%poly_v(2) = 1._dp + e2*tmp%arg_v(2) + e4*tmp%arg_v(2)*tmp%arg_v(2) +&
           e6*tmp%arg_v(2)**3 + e8*tmp%arg_v(2)**4
       tmp%poly_v(3) = 1._dp + e2*tmp%arg_v(2) + e4*tmp%arg_v(2)*tmp%arg_v(2) +&
           e6*tmp%arg_v(2)**3 + e8*tmp%arg_v(2)**4

       tmp%scale_r(1) = EXP( 0.5_dp * dt * npt ( 1, 1 )%v )
       tmp%scale_v(1) = EXP( -0.25_dp * dt * aax )
       tmp%scale_v(2) = EXP( -0.25_dp * dt * aa )
       tmp%scale_v(3) = EXP( -0.25_dp * dt * aa )

       ! first half of velocity verlet
       CALL vv_first( tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, error=error )


       IF(simpar%variable_dt)  CALL variable_timestep(md_env,tmp, dt, simpar,para_env,&
           atomic_kind_set, local_particles, particle_set, core_particle_set,&
           shell_particle_set, nparticle_kind, shell_adiabatic, npt=npt, error=error)


       roll_tol = 0._dp
       vector_r ( : ) = 0._dp
       vector_v ( : ) = tmp%scale_v(:) * tmp%poly_v(:)
       vector_r ( 1 ) = tmp%scale_r(1) * tmp%poly_r(1)

       IF ( simpar%constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos,tmp%vel, dt, simpar,&
            roll_tol, iroll, vector_r, vector_v, para_env%group, cell=cell,&
            local_particles=local_particles, error=error )
    END DO SR

    ! Update h_mat
    cell%hmat ( 1, 1 ) = cell%hmat ( 1, 1 ) * tmp%scale_r(1) * tmp%scale_r(1)

    ! Update the inverse
    CALL init_cell ( cell )

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE.,  error=error )

    ! Update forces
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env, error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second(tmp, atomic_kind_set, local_particles, particle_set, &
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, error=error )

    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       first = .TRUE.
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    RR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! RATTLE-ROLL LOOP
       roll_tol = 0._dp
       IF ( simpar%constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set, local_molecules,&
            tmp%vel, dt, cell, npt, simpar, virial, vector_v, roll_tol, iroll, infree, first,&
            para_env, error=error)
       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree )

       ! perform damping on the barostat momentum
       CALL damp_veps ( npt ( 1, 1 ), gamma1, dt )

    END DO RR

    ! perform damping on velocities
    CALL damp_v ( molecule_kind_set, molecule_set, particle_set, local_molecules,&
         tmp%vel, gamma1,  npt ( 1, 1 ), dt, para_env%group )

    IF (simpar%annealing)  tmp%vel(:,:)=tmp%vel(:,:)*simpar%f_annealing

    ! Broadcast the new particle velocities and deallocate temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, vel=.TRUE.,  error=error )

    ! Update constraint virial
    IF ( simpar%constraint ) CALL pv_constraint ( gci, local_molecules, &
         molecule_set, molecule_kind_set, particle_set, virial, para_env%group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group, error=error)

    ! Deallocate old variables
    CALL deallocate_old ( old, error=error)

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE nph_uniaxial_damped

! *****************************************************************************
!> \brief Velocity Verlet integrator for the NPT ensemble with fully flexible cell
!> \param md_env ...
!> \param globenv ...
!> \param error ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE npt_f ( md_env, globenv, error )

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'npt_f', &
      routineP = moduleN//':'//routineN
    LOGICAL, PARAMETER                       :: box_change = .TRUE.
    REAL(KIND=dp), PARAMETER                 :: e2 = 1.0_dp/6.0_dp, &
                                                e4 = e2/20.0_dp, &
                                                e6 = e4/42.0_dp, &
                                                e8 = e6/72.0_dp

    INTEGER                                  :: i, iroll, j, nparticle, &
                                                nparticle_kind, nshell
    LOGICAL                                  :: failure, first, first_time, &
                                                shell_adiabatic, &
                                                shell_check_distance, &
                                                shell_present
    REAL(KIND=dp)                            :: dt, infree, kin, roll_tol, &
                                                roll_tol_thrs, trvg
    REAL(KIND=dp), DIMENSION(3)              :: vector_r, vector_v
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_kin, uh
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(barostat_type), POINTER             :: barostat
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
    TYPE(simpar_type), POINTER               :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_baro, &
                                                thermostat_part, &
                                                thermostat_shell
    TYPE(tmp_variables_type), POINTER        :: tmp
    TYPE(virial_type), POINTER               :: virial

    failure = .FALSE.
    NULLIFY(gci,thermostat_baro,thermostat_part,thermostat_shell,force_env)
    NULLIFY(atomic_kinds,cell,para_env,subsys,local_molecules,local_particles)
    NULLIFY(molecule_kinds,molecules,molecule_kind_set,npt,barostat)
    NULLIFY(core_particles,particles,shell_particles,tmp,old)
    NULLIFY(core_particle_set,particle_set,shell_particle_set)
    NULLIFY(simpar,virial)

    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env,&
         thermostat_part=thermostat_part, thermostat_baro=thermostat_baro,&
         thermostat_shell=thermostat_shell, npt=npt, first_time=first_time,&
         para_env=para_env, barostat=barostat, error=error)
    dt = simpar%dt
    infree = 1.0_dp / REAL ( simpar%nfree,KIND=dp)

    CALL force_env_get(force_env, subsys=subsys, virial=virial, cell=cell ,error=error)

    ! Do some checks on coordinates and box
    CALL apply_qmmm_walls_reflective(force_env, error=error)

    CALL cp_subsys_get(subsys=subsys,atomic_kinds=atomic_kinds, local_particles=local_particles,&
         particles=particles,local_molecules_new=local_molecules,molecules_new=molecules, &
         gci=gci,molecule_kinds_new=molecule_kinds ,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         shell_present=shell_present,shell_adiabatic=shell_adiabatic,&
         shell_check_distance=shell_check_distance)

    IF ( first_time ) THEN
       CALL virial_evaluate ( atomic_kind_set, particle_set,  &
            local_particles, virial, para_env%group  ,error=error)
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

    CALL allocate_tmp(md_env, tmp,nparticle, nshell, shell_adiabatic, error=error)

    ! Apply Thermostat to Barostat
    CALL apply_thermostat_baro(thermostat_baro, npt, para_env%group, error)

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_f_ensemble) THEN
      IF(shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, local_particles, para_env%group,shell_adiabatic=shell_adiabatic,&
             shell_particle_set=shell_particle_set,core_particle_set=core_particle_set,&
             error=error)
      ELSE
        CALL apply_thermostat_particles(thermostat_part, molecule_kind_set, molecule_set,&
             particle_set, local_molecules, local_particles, para_env%group ,error=error)
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
    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       iroll = 1
       CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'F' )
       CALL getold ( gci, local_molecules, molecule_set,  &
            molecule_kind_set, particle_set, cell, error)
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    SR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! SHAKE-ROLL LOOP

       IF (simpar%constraint) THEN
          CALL set ( old, atomic_kind_set, particle_set, local_particles, cell, npt, 'B' )
       END IF
       CALL update_pv ( gci, simpar, atomic_kind_set, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree,&
                          virial_components=barostat%virial_components)

       trvg = npt ( 1, 1 )%v  + npt ( 2, 2 )%v + npt ( 3, 3 )%v
       !
       ! find eigenvalues and eigenvectors of npt ( :, : )%v
       !

       CALL diagonalise ( matrix = npt(:,:)%v, mysize = 3, &
            storageform = "UPPER", eigenvalues = tmp%e_val, eigenvectors = tmp%u )

       tmp%arg_r ( : ) =  0.5_dp * tmp%e_val ( : ) * dt * &
            0.5_dp * tmp%e_val ( : ) * dt
       tmp%poly_r = 1.0_dp + e2*tmp%arg_r + e4*tmp%arg_r*tmp%arg_r + &
                e6*tmp%arg_r**3 + e8*tmp%arg_r**4
       tmp%scale_r ( : ) = EXP ( 0.5_dp * dt * tmp%e_val ( : ) )

       tmp%arg_v ( : ) =  0.25_dp * dt * ( tmp%e_val ( : ) + trvg * infree ) *&
            0.25_dp * dt * ( tmp%e_val ( : ) + trvg * infree )
       tmp%poly_v = 1.0_dp + e2*tmp%arg_v + e4*tmp%arg_v*tmp%arg_v +&
            e6*tmp%arg_v**3 + e8*tmp%arg_v**4
       tmp%scale_v ( : ) = EXP ( -0.25_dp * dt * ( tmp%e_val( : ) + trvg * infree ) )

       CALL vv_first(tmp, atomic_kind_set, local_particles, particle_set, &
            core_particle_set, shell_particle_set, nparticle_kind,&
            shell_adiabatic, dt, u=tmp%u, error=error)

       IF(simpar%variable_dt) CALL  variable_timestep(md_env,tmp, dt, simpar,para_env,&
           atomic_kind_set, local_particles, particle_set, core_particle_set,&
           shell_particle_set, nparticle_kind, shell_adiabatic,npt=npt,error=error)

       roll_tol = 0.0_dp
       vector_r = tmp%scale_r*tmp%poly_r
       vector_v = tmp%scale_v*tmp%poly_v

       IF ( simpar%constraint ) CALL shake_roll_control( gci, local_molecules, &
            molecule_set, molecule_kind_set, particle_set, tmp%pos, tmp%vel, dt,&
            simpar, roll_tol, iroll, vector_r, vector_v, &
            para_env%group, u=tmp%u, cell=cell,&
            local_particles=local_particles, error=error )
    END DO SR

    ! Update h_mat
    uh = MATMUL_3X3 ( TRANSPOSE_3D ( tmp%u ), cell%hmat )

    DO i = 1, 3
       DO j = 1, 3
          uh ( i, j ) = uh ( i, j ) * tmp%scale_r ( i ) * tmp%scale_r ( i )
       END DO
    END DO

    cell%hmat = MATMUL_3x3 ( tmp%u, uh )
    ! Update the inverse
    CALL init_cell ( cell )

    ! Broadcast the new particle positions and deallocate the pos components of temporary
    CALL update_dealloc_tmp (tmp,particle_set, shell_particle_set, &
         core_particle_set, para_env, shell_adiabatic, pos=.TRUE., error=error )

    IF(shell_adiabatic .AND. shell_check_distance) THEN
        CALL optimize_shell_core(force_env,  particle_set,&
             shell_particle_set, core_particle_set, globenv, tmp=tmp, check=.TRUE., error=error)
    END IF

    ! Update forces
    CALL force_env_set_cell(force_env, cell, error)
    CALL force_env_calc_energy_force(force_env,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, tmp%vel, error=error)

    ! Velocity Verlet (second part)
    CALL vv_second (tmp, atomic_kind_set, local_particles, particle_set,&
         core_particle_set, shell_particle_set, nparticle_kind,&
         shell_adiabatic, dt, tmp%u, error=error )

    IF (simpar%constraint) THEN
       roll_tol_thrs = simpar%roll_tol
       first = .TRUE.
       iroll = 1
       CALL set (old, atomic_kind_set, particle_set, tmp%vel, local_particles, cell, npt, 'F' )
    ELSE
       roll_tol_thrs = EPSILON(0.0_dp)
    ENDIF
    roll_tol = -roll_tol_thrs

    RR: DO WHILE (ABS(roll_tol)>=roll_tol_thrs) ! RATTLE-ROLL LOOP
       roll_tol = 0.0_dp
       IF ( simpar%constraint ) CALL rattle_roll_setup ( old, gci, atomic_kind_set, &
            particle_set, local_particles, molecule_kind_set, molecule_set,  &
            local_molecules, tmp%vel, dt, cell, npt, simpar, virial, vector_v,  &
            roll_tol, iroll, infree, first, para_env, u=tmp%u, error=error)

       CALL update_pv ( gci, simpar, atomic_kind_set, tmp%vel, particle_set, &
            local_molecules, molecule_set, molecule_kind_set, &
            local_particles, kin, pv_kin, virial, para_env%group )
       CALL update_veps ( cell, npt, simpar, pv_kin, kin, virial, infree,&
                          virial_components=barostat%virial_components)
    END DO RR

    ! Apply Thermostat over the full set of particles
    IF(simpar% ensemble /= npe_f_ensemble) THEN
      IF (shell_adiabatic) THEN
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, local_particles, para_env%group,shell_adiabatic=shell_adiabatic,&
             vel=tmp%vel, shell_vel=tmp%shell_vel, core_vel=tmp%core_vel, error=error)

      ELSE
        CALL apply_thermostat_particles(thermostat_part,molecule_kind_set, molecule_set, &
             particle_set, local_molecules, local_particles, para_env%group, vel=tmp%vel, error=error)
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
    IF ( simpar%constraint ) &
         CALL pv_constraint ( gci, local_molecules, molecule_set, &
         molecule_kind_set, particle_set, virial, para_env%group )

    CALL virial_evaluate ( atomic_kind_set, particle_set,  &
         local_particles, virial, para_env%group, error=error  )

    ! Deallocate old variables
    CALL deallocate_old ( old, error=error )

    IF (first_time) THEN
       first_time = .FALSE.
       CALL set_md_env(md_env, first_time=first_time,error=error)
    END IF

  END SUBROUTINE npt_f

! *****************************************************************************
!> \brief RESPA integrator for nve ensemble for particle positions & momenta
!> \param md_env ...
!> \param error ...
!> \author FS
! *****************************************************************************
  SUBROUTINE nve_respa ( md_env, error)

    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nve_respa', &
      routineP = moduleN//':'//routineN

    INTEGER :: i_step, iparticle, iparticle_kind, iparticle_local, &
      n_time_steps, nparticle, nparticle_kind, nparticle_local, stat
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
    TYPE(simpar_type), POINTER               :: simpar

    failure = .FALSE.

    NULLIFY (para_env,cell,subsys_respa,particles_respa, particle_set_respa, gci, force_env, atomic_kinds)
    NULLIFY (atomic_kind_set,simpar,subsys,particles,particle_set)
    NULLIFY (local_molecules,molecule_kinds, molecules,molecule_kind_set, local_particles)
    CALL get_md_env(md_env=md_env,simpar=simpar,force_env=force_env, para_env=para_env, error=error)
    dt = simpar%dt

    n_time_steps=simpar%n_time_steps

    CALL force_env_get(force_env,subsys=subsys,cell=cell,error=error)
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
    molecule_kind_set => molecule_kinds%els

    nparticle = particles%n_els
    particle_set => particles%els
    molecule_set => molecules%els

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
       CALL update_particle_set ( particle_set, para_env%group, pos = pos, error=error)
       DO iparticle = 1, SIZE(particle_set)
          particle_set_respa(iparticle)%r = particle_set(iparticle)%r
       END DO

       ! Update forces
       CALL force_env_calc_energy_force(force_env%sub_force_env(1)%force_env,error=error)

       ! Metadynamics
       CALL metadyn_integrator(force_env, vel, error=error)

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
    END DO
    DEALLOCATE (pos,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Multiple time step (second part)
    ! Compute forces for respa force_env
    CALL force_env_calc_energy_force(force_env,error=error)

    ! Metadynamics
    CALL metadyn_integrator(force_env, vel, error=error)

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
    CALL update_particle_set ( particle_set, para_env%group, vel = vel, error=error)

    DEALLOCATE (vel,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE nve_respa

END MODULE integrator
