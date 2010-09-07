!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2010  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief evaluete the potential energy and its gradients using an array
!>      with same dimension as the particle_set
!> \param gopt_env the geometry optimization environment
!> \param x the position where the function should be evaluated
!> \param f the function value
!> \param gradient the value of its gradient
!> \param error variable to control error logging, stopping,...
!>             see module cp_error_handling
!> \par History
!>       CELL OPTIMIZATION:  Teodoro Laino [tlaino] - University of Zurich - 03.2008
!>
!> \author Teodoro Laino [tlaino] - University of Zurich - 01.2008
! *****************************************************************************
  RECURSIVE SUBROUTINE cp_eval_at(gopt_env, x, f, gradient, master, para_env, error)
    USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                               get_atomic_kind
    USE averages_types,                  ONLY: average_quantities_type,&
                                               create_averages,&
                                               release_averages
    USE cell_opt_utils,                  ONLY: get_dg_dh,&
                                               gopt_new_error_create,&
                                               gopt_new_error_release
    USE cell_types,                      ONLY: cell_type,&
                                               init_cell,&
                                               real_to_scaled,&
                                               scaled_to_real,&
                                               write_cell
    USE cp_para_types,                   ONLY: cp_para_env_type
    USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                               cp_subsys_type
    USE dimer_methods,                   ONLY: cp_eval_at_ts
    USE f77_blas
    USE force_env_methods,               ONLY: force_env_calc_energy_force
    USE force_env_types,                 ONLY: force_env_get,&
                                               force_env_set_cell
    USE geo_opt,                         ONLY: cp_geo_opt
    USE gopt_f_types,                    ONLY: gopt_f_type
    USE input_constants,                 ONLY: default_minimization_method_id,&
                                               default_ts_method_id,&
                                               default_cell_method_id,&
                                               default_cell_geo_opt_id,&
                                               default_cell_md_id,&
                                               default_shellcore_method_id,&
                                               nvt_ensemble,&
                                               mol_dyn_run,&
                                               geo_opt_run,&
                                               cell_opt_run
    USE input_section_types,             ONLY: section_vals_get,&
                                               section_vals_get_subs_vals,&
                                               section_vals_type,&
                                               section_vals_val_get
    USE md_run,                          ONLY: qs_mol_dyn
    USE kinds,                           ONLY: dp,&
                                               default_string_length
    USE particle_list_types,             ONLY: particle_list_type
    USE particle_types,                  ONLY: write_structure_data
    USE shell_potential_types,           ONLY: shell_kind_type
    USE timings,                         ONLY: timeset,&
                                               timestop
    USE virial_methods,                  ONLY: virial_update
    USE virial_types,                    ONLY: cp_virial,&
                                               virial_create,&
                                               virial_release,&
                                               virial_type
                                               USE physcon
#include "cp_common_uses.h"
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    INTEGER, INTENT(in)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at', moduleN='gopt_f77_methods',&
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ensemble, handle, i, idg, idir, ip, shell_index, stat
    LOGICAL                                  :: failure, explicit, is_shell
    REAL(KIND=dp)                            :: f_ts, mass, potential_energy,&
    mass_core, mass_shell
    REAL(KIND=dp), DIMENSION(3)              :: r, s
    REAL(KIND=dp), DIMENSION(3,3)            :: av_ptens
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient_ts
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_list_type), POINTER        :: shell_particles
    TYPE(particle_list_type), POINTER        :: core_particles
    TYPE(virial_type), POINTER               :: virial, virial_avg
    TYPE(cp_error_type)                      :: new_error
    CHARACTER(LEN=default_string_length)     :: project_name
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(section_vals_type), POINTER         :: work, avgs_section
    TYPE(shell_kind_type), POINTER           :: shell

    failure = .FALSE.

    NULLIFY (averages)
    NULLIFY (cell)
    NULLIFY (core_particles)
    NULLIFY (gradient_ts)
    NULLIFY (particles)
    NULLIFY (shell_particles)
    NULLIFY (subsys)
    NULLIFY (virial)

    CALL timeset(routineN,handle)
    CALL force_env_get(gopt_env%force_env,subsys=subsys,cell=cell,error=error)
    CALL cp_subsys_get(subsys,&
                       core_particles=core_particles,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       error=error)
    SELECT CASE (gopt_env%type_id)
    CASE (default_minimization_method_id,default_ts_method_id)
      idg = 0
      DO ip=1,particles%n_els
        shell_index = particles%els(ip)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            idg = idg + 1
            particles%els(ip)%r(i) = x(idg)
          END DO
        END IF
      END DO
      DO ip=1,particles%n_els
        shell_index = particles%els(ip)%shell_index
        IF (shell_index /= 0) THEN
          DO i=1,3
            idg = idg + 1
            core_particles%els(shell_index)%r(i) = x(idg)
          END DO
          DO i=1,3
            idg = idg + 1
            shell_particles%els(shell_index)%r(i) = x(idg)
          END DO
          mass_core = core_particles%els(shell_index)%atomic_kind%shell%mass_core
          mass_shell = shell_particles%els(shell_index)%atomic_kind%shell%mass_shell
          mass = particles%els(ip)%atomic_kind%mass
          particles%els(ip)%r(:) = (mass_core*core_particles%els(shell_index)%r(:) +&
                                    mass_shell*shell_particles%els(shell_index)%r(:))/mass
        END IF
      END DO
      CALL write_structure_data(particles%els,cell,gopt_env%motion_section,error)
      SELECT CASE (gopt_env%type_id)
      CASE (default_minimization_method_id)
        ! Geometry Minimization
        CALL force_env_calc_energy_force(gopt_env%force_env,&
                                         calc_force=PRESENT(gradient),&
                                         error=error)
        ! Possibly take the potential energy
        IF (PRESENT(f)) THEN
          CALL force_env_get(gopt_env%force_env,potential_energy=f,error=error)
        END IF
        ! Possibly take the gradients
        IF (PRESENT(gradient)) THEN
          IF (master == para_env%mepos) THEN ! we are on the master
            CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
            CALL force_env_get(gopt_env%force_env,subsys=subsys,error=error)
            CALL cp_subsys_get(subsys,&
                               core_particles=core_particles,&
                               particles=particles,&
                               shell_particles=shell_particles,&
                               error=error)
            idg = 0
            DO ip=1,particles%n_els
              shell_index = particles%els(ip)%shell_index
              IF (shell_index == 0) THEN
                DO i=1,3
                  idg = idg + 1
                  CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
                  gradient(idg) = -particles%els(ip)%f(i)
                END DO
              END IF
            END DO
            DO ip=1,particles%n_els
              shell_index = particles%els(ip)%shell_index
              IF (shell_index /= 0) THEN
                DO i=1,3
                  idg = idg + 1
                  CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
                  gradient(idg) = -core_particles%els(shell_index)%f(i)
                END DO
                DO i=1,3
                  idg = idg + 1
                  CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
                  gradient(idg) = -shell_particles%els(shell_index)%f(i)
                END DO
              END IF
            END DO
          END IF
        END IF
      CASE (default_ts_method_id)
          ! Transition State Optimization
          ALLOCATE(gradient_ts(particles%n_els*3),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ! Real calculation of energy and forces for transition state optimization:
          ! When doing dimer methods forces have to be always computed since the function
          ! to minimize is not the energy but the effective force
          CALL cp_eval_at_ts (gopt_env, x, f_ts, gradient_ts, calc_force=.TRUE., error=error)
          ! Possibly take the potential energy
          IF (PRESENT(f)) f = f_ts

          ! Possibly take the gradients
          IF (PRESENT(gradient)) THEN
             IF (master==para_env%mepos) THEN ! we are on the master
                CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
                gradient = gradient_ts
             END IF
          END IF
          DEALLOCATE(gradient_ts, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
      END SELECT
      ! This call is necessary for QM/MM if a Translation is applied
      ! this makes the geometry optimizer consistent
      idg = 0
      DO ip=1,particles%n_els
        shell_index = particles%els(ip)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            idg = idg + 1
            particles%els(ip)%r(i) = x(idg)
          END DO
        END IF
      END DO
      DO ip=1,particles%n_els
        shell_index = particles%els(ip)%shell_index
        IF (shell_index /= 0) THEN
          DO i=1,3
            idg = idg + 1
            core_particles%els(shell_index)%r(i) = x(idg)
          END DO
          DO i=1,3
            idg = idg + 1
            shell_particles%els(shell_index)%r(i) = x(idg)
          END DO
          mass_core = core_particles%els(shell_index)%atomic_kind%shell%mass_core
          mass_shell = shell_particles%els(shell_index)%atomic_kind%shell%mass_shell
          mass = particles%els(ip)%atomic_kind%mass
          particles%els(ip)%r(:) = (mass_core*core_particles%els(shell_index)%r(:) +&
                                    mass_shell*shell_particles%els(shell_index)%r(:))/mass
        END IF
      END DO
    CASE(default_cell_method_id)
       ! Check for VIRIAL
       CALL force_env_get(gopt_env%force_env,virial=virial,error=error)
       CALL cp_assert(virial%pv_availability,cp_failure_level,cp_assertion_failed,routineP,&
            "Cell optimization requested but FORCE_EVAL%STRESS_TENSOR was not defined! "//&
            "Activate the evaluation of the stress tensor for cell optimization! "//&
CPSourceFileRef,&
            only_ionode=.TRUE.)

       ! Getting Scaled coordinates
       DO ip=1,particles%n_els
          CALL real_to_scaled(s,particles%els(ip)%r,cell)
          particles%els(ip)%r = s
       END DO

       ! Set the new cell (according the present optimization step)
       cell%hmat = 0.0_dp
       idg = 0
       DO i = 1, 3
          DO j = 1, i
             idg = idg + 1
             cell%hmat(j,i) = x(idg)
          END DO
       END DO
       CALL init_cell(cell)
       CALL force_env_set_cell(gopt_env%force_env, cell, error)

       ! Converting into reals
       DO ip=1,particles%n_els
          CALL scaled_to_real(r,particles%els(ip)%r,cell)
          particles%els(ip)%r = r
       END DO

       ! Possibly output the new cell used for the next calculation
       CALL write_cell(cell,gopt_env%geo_section,error=error)

       !Compute the pressure tensor
       CALL virial_create(virial_avg, error)
       SELECT CASE(gopt_env%cell_method_id)
       CASE(default_cell_geo_opt_id)
          work => section_vals_get_subs_vals(gopt_env%motion_section,"GEO_OPT",error=error)
          CALL section_vals_get(work,explicit=explicit,error=error)
          CALL cp_assert(explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "Cell optimization at 0K was requested. GEO_OPT section MUST be provided in the "//&
               "input file! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ! Perform a geometry optimization
          CALL gopt_new_error_create(new_error, gopt_env%force_env%root_section, para_env,&
               project_name, id_run=geo_opt_run, error=error)
          CALL cp_geo_opt(gopt_env%force_env, gopt_env%globenv, eval_opt_geo=.FALSE., error=new_error)
          CALL force_env_get(gopt_env%force_env,potential_energy=potential_energy,virial=virial,error=error)
          CALL cp_virial(virial, virial_avg)
       CASE(default_cell_md_id)
          work        => section_vals_get_subs_vals(gopt_env%motion_section,"MD",error=error)
          avgs_section => section_vals_get_subs_vals(work,"AVERAGES",error=error)
          CALL section_vals_get(work,explicit=explicit,error=error)
          CALL cp_assert(explicit,cp_failure_level,cp_assertion_failed,routineP,&
               "Cell optimization at finite temperature was requested. MD section MUST be provided in the "//&
               "input file! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ! Only NVT ensemble is allowed..
          CALL section_vals_val_get(gopt_env%motion_section,"MD%ENSEMBLE",i_val=ensemble,error=error)
          CALL cp_assert(ensemble==nvt_ensemble,cp_failure_level,cp_assertion_failed,routineP,&
               "Cell optimization at finite temperature requires the NVT MD ensemble! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ! Perform a molecular dynamics
          CALL gopt_new_error_create(new_error, gopt_env%force_env%root_section, para_env,&
               project_name, id_run=mol_dyn_run, error=error)
          CALL create_averages(averages, avgs_section, virial_avg=.TRUE., force_env=gopt_env%force_env,&
               error=error)
          CALL qs_mol_dyn(gopt_env%force_env, gopt_env%globenv, averages, rm_restart_info=.FALSE.,&
                          error=new_error)
          ! Retrieve the average of the stress tensor and the average of the potential energy
          potential_energy=averages%avepot
          CALL cp_virial(averages%virial, virial_avg)
          CALL release_averages(averages, error)
       CASE DEFAULT
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       CALL gopt_new_error_release(new_error, gopt_env%force_env%root_section, para_env, project_name,&
            cell_opt_run, error)

       ! Update the virial
       CALL virial_update(virial_avg, subsys, para_env, error)

       ! Possibly take give back the potential energy
       IF (PRESENT(f)) THEN
          f = potential_energy
       END IF

       ! Possibly give back the gradients
       IF (PRESENT(gradient)) THEN
          CPPrecondition(ANY(virial_avg%pv_total/=0),cp_failure_level,routineP,error,failure)
          ! Convert the average ptens
          av_ptens(:,:) = virial_avg%pv_total(:,:)/cell%deth
          IF (master==para_env%mepos) THEN ! we are on the master
             CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
             ! Compute the gradients on the cell
             CALL get_dg_dh(gradient,av_ptens,gopt_env%cell_env%pres_ext,cell,gopt_env%cell_env%mtrx,&
                  keep_angles=gopt_env%cell_env%keep_angles,pres_int=gopt_env%cell_env%pres_int,error=error)
          END IF
       END IF
       CALL virial_release(virial_avg, error)
    CASE(default_shellcore_method_id)
       CALL cp_subsys_get(subsys, shell_particles=shell_particles,&
            core_particles=core_particles,  error=error)
       idg=0
       DO ip=1,particles%n_els
         NULLIFY (atomic_kind)
         atomic_kind => particles%els(ip)%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind,shell_active=is_shell,mass=mass,shell=shell)
         IF (is_shell) THEN
           shell_index = particles%els(ip)%shell_index
           DO idir=1,3
              idg = 3*(shell_index-1) + idir
              shell_particles%els(shell_index)%r(idir) = core_particles%els(ip)%r(idir) - x(idg)
           END DO
         END IF
       END DO
       CALL write_structure_data(particles%els,cell,gopt_env%motion_section,error)

       ! Shell-core optimization
       CALL force_env_calc_energy_force(gopt_env%force_env,&
                                        calc_force=PRESENT(gradient),&
                                        error=error)

       ! Possibly take the potential energy
       IF (PRESENT(f)) THEN
          CALL force_env_get(gopt_env%force_env,potential_energy=f,error=error)
       END IF

       ! Possibly take the gradients
       IF (PRESENT(gradient)) THEN
         IF (master==para_env%mepos) THEN ! we are on the master
            CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
            CALL force_env_get(gopt_env%force_env, subsys=subsys,error=error)
            idg = 0
            CALL cp_subsys_get(subsys, shell_particles=shell_particles, core_particles=core_particles, error=error)
            DO ip=1,shell_particles%n_els
               shell_index = core_particles%els(ip)%shell_index
               DO idir=1,3
                 idg = idg + 1
                 gradient(idg) = -(core_particles%els(ip)%f(idir) - shell_particles%els(ip)%f(idir))
               END DO
            END DO
         END IF
       END IF
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE cp_eval_at
