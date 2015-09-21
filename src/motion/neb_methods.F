!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Module performing a Nudged Elastic Band Calculation
!> \note
!>      Numerical accuracy for parallel runs:
!>       Each replica starts the SCF run from the one optimized
!>       in a previous run. It may happen then energies and derivatives
!>       of a serial run and a parallel run could be slightly different
!>       'cause of a different starting density matrix.
!>       Exact results are obtained using:
!>          EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!> \author Teodoro Laino 09.2006
!> \par  History
!>       - Teodoro Laino 10.2008 [tlaino] - University of Zurich
!>         Extension to a subspace of collective variables
! *****************************************************************************
MODULE neb_methods
  USE colvar_utils,                    ONLY: number_of_colvar
  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: band_header
  USE input_constants,                 ONLY: band_diis_opt,&
                                             band_md_opt,&
                                             do_rep_blocked,&
                                             do_sm
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE neb_io,                          ONLY: dump_neb_info,&
                                             neb_rep_env_map_info,&
                                             read_neb_section
  USE neb_md_utils,                    ONLY: control_vels_a,&
                                             control_vels_b
  USE neb_opt_utils,                   ONLY: accept_diis_step,&
                                             neb_ls
  USE neb_types,                       ONLY: neb_type,&
                                             neb_var_create,&
                                             neb_var_release,&
                                             neb_var_type
  USE neb_utils,                       ONLY: build_replica_coords,&
                                             check_convergence,&
                                             neb_calc_energy_forces,&
                                             reorient_images,&
                                             reparametrize_images
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: massunit
  USE replica_methods,                 ONLY: rep_env_create
  USE replica_types,                   ONLY: rep_env_release,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_methods'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
  PUBLIC :: neb

CONTAINS

! *****************************************************************************
!> \brief Real subroutine for NEB calculations
!> \param input ...
!> \param input_declaration ...
!> \param para_env ...
!> \param globenv ...
!> \author Teodoro Laino 09.2006
!> \note
!>      Based on the use of replica_env
! *****************************************************************************
  SUBROUTINE neb(input, input_declaration, para_env, globenv)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'neb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, iw, iw2, nrep, &
                                                output_unit, prep, &
                                                proc_dist_type
    LOGICAL                                  :: check, failure, row_force
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(neb_var_type), POINTER              :: coords, forces, vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER :: diis_section, force_env_section, &
      md_section, motion_section, neb_section, print_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger, subsys, f_env, rep_env)
    NULLIFY(forces, coords, vels, neb_env)
    logger => cp_get_default_logger()
    CALL cp_add_iter_level(logger%iter_info,"BAND")
    motion_section => section_vals_get_subs_vals(input,"MOTION")
    print_section  => section_vals_get_subs_vals(motion_section,"PRINT")
    neb_section    => section_vals_get_subs_vals(motion_section,"BAND")
    output_unit=cp_print_key_unit_nr(logger,neb_section,"PROGRAM_RUN_INFO",&
         extension=".nebLog")
    CALL section_vals_val_get(neb_section,"NPROC_REP",i_val=prep)
    CALL section_vals_val_get(neb_section,"PROC_DIST_TYPE",i_val=proc_dist_type)
    row_force = (proc_dist_type==do_rep_blocked)
    nrep = MAX(1,para_env%num_pe/prep)
    CALL cp_assert(.NOT.(nrep*prep/=para_env%num_pe.AND.output_unit>0),&
         cp_warning_level,cp_assertion_failed,routineP,&
         "Number of totally requested processors ("//TRIM(ADJUSTL(cp_to_string(para_env%num_pe)))//") "//&
         "is not compatible with the number of processors requested per replica ("//&
         TRIM(ADJUSTL(cp_to_string(prep)))//") and the number of replicas ("//&
         TRIM(ADJUSTL(cp_to_string(nrep)))//") . ["//&
         TRIM(ADJUSTL(cp_to_string(para_env%num_pe-nrep*prep)))//"] processors will be wasted! "//&
         CPSourceFileRef, only_ionode=.TRUE.)
    force_env_section => section_vals_get_subs_vals(input,"FORCE_EVAL")
    ! Create Replica Environments
    IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Replica_env Setup. START"
    CALL rep_env_create(rep_env, para_env=para_env, input=input,&
         input_declaration=input_declaration,nrep=nrep,prep=prep,row_force=row_force)
    CPASSERT(SIZE(rep_env%local_rep_indices)==1)
    IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Replica_env Setup. END"
    IF (ASSOCIATED(rep_env)) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            failure=failure)
       CALL force_env_get(f_env%force_env,subsys=subsys)
       particle_set => subsys%particles%els
       ! Read NEB controlling parameters
       ALLOCATE (neb_env)
       neb_env%force_env            => f_env%force_env
       neb_env%root_section         => input
       neb_env%force_env_section    => force_env_section
       neb_env%motion_print_section => print_section
       neb_env%neb_section          => neb_section
       neb_env%nsize_xyz            = rep_env%ndim
       neb_env%nsize_int            = number_of_colvar(f_env%force_env)
       check = (neb_env%nsize_xyz >= neb_env%nsize_int)
       CPASSERT(check)
       ! Check that teh used colvar are uniquely determined
       check = (number_of_colvar(f_env%force_env) == &
                number_of_colvar(f_env%force_env, unique=.TRUE.))
       CPASSERT(check)
       CALL read_neb_section(neb_env,neb_section)
       ! Print BAND header
       iw2=cp_print_key_unit_nr(logger,neb_section,"BANNER",extension=".nebLog")
       CALL band_header(iw2, neb_env%number_of_replica, nrep, prep)
       CALL cp_print_key_finished_output(iw2,logger,neb_section,"BANNER")
       ! Allocate the principal vectors used in the BAND calculation
       CALL neb_var_create(coords, neb_env, full_allocation=.TRUE.)
       CALL neb_var_create(forces, neb_env)
       CALL neb_var_create(vels,   neb_env)
       ! Collecting the coordinates of the starting replicas of the BAND calculation
       IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Building initial set of coordinates. START"
       iw = cp_print_key_unit_nr(logger,neb_section,"PROGRAM_RUN_INFO/INITIAL_CONFIGURATION_INFO",&
         extension=".nebLog")
       CALL build_replica_coords(neb_section,particle_set,coords,vels,neb_env,iw,globenv,&
                                 rep_env%para_env)
       CALL cp_print_key_finished_output(iw,logger,neb_section,&
            "PROGRAM_RUN_INFO/INITIAL_CONFIGURATION_INFO")
       IF (output_unit > 0) WRITE (output_unit,'(T2,"NEB|",A)')" Building initial set of coordinates. END"
       ! Print some additional info in the replica_env initialization file
       CALL neb_rep_env_map_info(rep_env, neb_env)
       ! Perform NEB optimization
       SELECT CASE(neb_env%opt_type)
       CASE (band_md_opt)
          neb_env%opt_type_label = "MOLECULAR DYNAMICS"
          md_section => section_vals_get_subs_vals(neb_section,"OPTIMIZE_BAND%MD")
          CALL neb_md(rep_env, neb_env, coords, vels, forces, particle_set, output_unit,&
                      md_section, logger, globenv)
       CASE (band_diis_opt)
          neb_env%opt_type_label = "DIIS"
          diis_section => section_vals_get_subs_vals(neb_section,"OPTIMIZE_BAND%DIIS")
          CALL neb_diis(rep_env, neb_env, coords, vels, forces, particle_set, output_unit,&
                        diis_section, logger, globenv)
       END SELECT
       ! Release force_eval
       CALL f_env_rm_defaults(f_env,ierr)
       ! Release coords, vels and forces
       CALL neb_var_release(coords)
       CALL neb_var_release(forces)
       CALL neb_var_release(vels)
       ! At the end let's destroy the environment of the BAND calculation
       DEALLOCATE(neb_env)
    END IF
    CALL rep_env_release(rep_env)
    CALL cp_print_key_finished_output(output_unit,logger,neb_section,&
         "PROGRAM_RUN_INFO")
    CALL cp_rm_iter_level(logger%iter_info,"BAND")
    CALL timestop(handle)
  END SUBROUTINE neb

! *****************************************************************************
!> \brief MD type optimization NEB
!> \param rep_env ...
!> \param neb_env ...
!> \param coords ...
!> \param vels ...
!> \param forces ...
!> \param particle_set ...
!> \param output_unit ...
!> \param md_section ...
!> \param logger ...
!> \param globenv ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_md(rep_env, neb_env, coords, vels, forces, particle_set, output_unit,&
       md_section, logger, globenv)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), OPTIONAL, POINTER        :: neb_env
    TYPE(neb_var_type), POINTER              :: coords, vels, forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_md', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ic, is, istep, &
                                                iw, max_steps, natom, &
                                                shell_index
    LOGICAL                                  :: converged, failure, &
                                                should_stop
    REAL(KIND=dp)                            :: dt
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distances, energies
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mass
    TYPE(neb_var_type), POINTER              :: Dcoords
    TYPE(section_vals_type), POINTER         :: tc_section, vc_section

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY( Dcoords, tc_section, vc_section)
    CPASSERT(ASSOCIATED(coords))
    CPASSERT(ASSOCIATED(vels))
    ! MD band for string methods type does not make anywa sense. Stop calculation.
    IF (neb_env%id_type==do_sm) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="MD band optimization and String Method incompatible.", &
            error_level=cp_warning_level)
    END IF
    ! Output unit
    iw=cp_print_key_unit_nr(logger,neb_env%neb_section,"REPLICA_INFO",&
         extension=".replicaLog")
    tc_section => section_vals_get_subs_vals(md_section,"TEMP_CONTROL")
    vc_section => section_vals_get_subs_vals(md_section,"VEL_CONTROL")
    CALL section_vals_val_get(md_section,"TIMESTEP",r_val=dt)
    CALL section_vals_val_get(md_section,"MAX_STEPS",i_val=max_steps)
    ! Initial setup for MD
    CALL neb_var_create(Dcoords, neb_env)
    ALLOCATE (mass(SIZE(coords%wrk,1),neb_env%number_of_replica))
    ALLOCATE (energies(neb_env%number_of_replica))
    ALLOCATE (distances(neb_env%number_of_replica-1))
    ! Setting up the mass array
    IF (neb_env%use_colvar) THEN
       mass(:,:) = 0.5_dp*dt/massunit
    ELSE
       natom = SIZE(particle_set)
       DO iatom=1,natom
          ic = 3*(iatom - 1)
          shell_index = particle_set(iatom)%shell_index
          IF (shell_index == 0) THEN
             mass(ic+1:ic+3,:) = 0.5_dp*dt/particle_set(iatom)%atomic_kind%mass
          ELSE
             is = 3*(natom + shell_index - 1)
             mass(ic+1:ic+3,:) = 0.5_dp*dt/particle_set(iatom)%atomic_kind%shell%mass_core
             mass(is+1:is+3,:) = 0.5_dp*dt/particle_set(iatom)%atomic_kind%shell%mass_shell
          END IF
       END DO
    END IF
    ! Initializing forces array
    CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
         output_unit, distances, neb_env%number_of_replica)
    neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
    CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
         particle_set, iw)

    CALL dump_neb_info(neb_env=neb_env,&
                       coords=coords,&
                       vels=vels,&
                       forces=forces,&
                       particle_set=particle_set,&
                       logger=logger,&
                       istep=0,&
                       energies=energies,&
                       distances=distances,&
                       output_unit=output_unit)
    md_opt_loop: DO istep = 1, max_steps
       CALL cp_iterate(logger%iter_info,iter_nr=istep)
       ! Save the optimization step counter
       neb_env%istep = istep
       ! Velocity Verlet (first part)
       vels%wrk(:,:) = vels%wrk(:,:) + mass(:,:)*forces%wrk(:,:)
       ! Control on velocity - I part [rescale, annealing]
       CALL control_vels_a(vels, particle_set, tc_section, vc_section, output_unit,&
            istep)
       ! Coordinate step
       Dcoords%wrk(:,:) = dt*vels%wrk(:,:)
       coords%wrk(:,:) = coords%wrk(:,:) + Dcoords%wrk(:,:)

       CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
            output_unit, distances, neb_env%number_of_replica)
       neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
       CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
            particle_set, iw)
       ! Check for an external exit command
       CALL external_control(should_stop,"NEB",globenv=globenv)
       IF(should_stop) EXIT
       ! Control on velocity - II part [check vels VS forces, Steepest Descent like]
       CALL control_vels_b(vels, forces, vc_section)
       ! Velocity Verlet (second part)
       vels%wrk(:,:) = vels%wrk(:,:) + mass(:,:)*forces%wrk(:,:)
       ! Dump Infos
       CALL dump_neb_info(neb_env=neb_env,&
                          coords=coords,&
                          vels=vels,&
                          forces=forces,&
                          particle_set=particle_set,&
                          logger=logger,&
                          istep=istep,&
                          energies=energies,&
                          distances=distances,&
                          output_unit=output_unit)
       converged = check_convergence(neb_env, Dcoords, forces)
       IF (converged) EXIT
    END DO md_opt_loop

    DEALLOCATE (mass)
    DEALLOCATE (energies)
    DEALLOCATE (distances)
    CALL neb_var_release(Dcoords)
    CALL cp_print_key_finished_output(iw,logger,neb_env%neb_section,&
         "REPLICA_INFO")
    CALL timestop(handle)

  END SUBROUTINE neb_md

! *****************************************************************************
!> \brief DIIS type optimization NEB
!> \param rep_env ...
!> \param neb_env ...
!> \param coords ...
!> \param vels ...
!> \param forces ...
!> \param particle_set ...
!> \param output_unit ...
!> \param diis_section ...
!> \param logger ...
!> \param globenv ...
!> \author Teodoro Laino 09.2006
! *****************************************************************************
  SUBROUTINE neb_diis(rep_env, neb_env, coords, vels, forces, particle_set, output_unit,&
       diis_section, logger, globenv)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), OPTIONAL, POINTER        :: neb_env
    TYPE(neb_var_type), POINTER              :: coords, vels, forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(section_vals_type), POINTER         :: diis_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_diis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istep, iw, iw2, &
                                                max_sd_steps, max_steps, &
                                                n_diis
    INTEGER, DIMENSION(:), POINTER           :: set_err
    LOGICAL                                  :: check_diis, converged, &
                                                diis_on, do_ls, failure, &
                                                should_stop, skip_ls
    REAL(KIND=dp)                            :: max_stepsize, norm, stepsize, &
                                                stepsize0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distances, energies
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: crr, err
    TYPE(neb_var_type), POINTER              :: sline

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(sline, crr, err)
    neb_env%opt_type_label = "SD"
    do_ls = .TRUE.
    CPASSERT(ASSOCIATED(coords))
    CPASSERT(ASSOCIATED(vels))
    CPASSERT(ASSOCIATED(forces))
    iw=cp_print_key_unit_nr(logger,neb_env%neb_section,"REPLICA_INFO",&
         extension=".replicaLog")
    CALL section_vals_val_get(diis_section,"MAX_STEPS",i_val=max_steps)
    CALL section_vals_val_get(diis_section,"N_DIIS",i_val=n_diis)
    CALL section_vals_val_get(diis_section,"STEPSIZE",r_val=stepsize0)
    CALL section_vals_val_get(diis_section,"MAX_STEPSIZE",r_val=max_stepsize)
    CALL section_vals_val_get(diis_section,"NO_LS",l_val=skip_ls)
    CALL section_vals_val_get(diis_section,"MAX_SD_STEPS",i_val=max_sd_steps)
    CALL section_vals_val_get(diis_section,"CHECK_DIIS",l_val=check_diis)
    iw2 = cp_print_key_unit_nr(logger,diis_section,"DIIS_INFO",&
         extension=".diisLog")
    ! Initial setup for DIIS
    stepsize = stepsize0
    ! Allocate type for Line Search direction
    CALL neb_var_create(sline, neb_env, full_allocation=.TRUE.)
    ! Array of error vectors
    ALLOCATE (err(PRODUCT(coords%size_wrk),n_diis))
    ALLOCATE (crr(PRODUCT(coords%size_wrk),n_diis))
    ALLOCATE (set_err(n_diis))
    ALLOCATE (energies(neb_env%number_of_replica))
    ALLOCATE (distances(neb_env%number_of_replica-1))
    ! Initializing forces array
    CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
         output_unit, distances, neb_env%number_of_replica)
    CALL reparametrize_images(neb_env%reparametrize_frames, neb_env%spline_order,&
         neb_env%smoothing, coords%wrk, sline%wrk, distances)
    neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
    CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
         particle_set, iw)
    ! Dump Infos
    CALL dump_neb_info(neb_env=neb_env,&
                       coords=coords,&
                       forces=forces,&
                       particle_set=particle_set,&
                       logger=logger,&
                       istep=0,&
                       energies=energies,&
                       distances=distances,&
                       vels=vels,&
                       output_unit=output_unit)
    ! If rotation is requested let's apply it at the beginning of the
    ! Geometry optimization and then let's disable it
    neb_env%rotate_frames= .FALSE.
    ! Main SD/DIIS loop
    set_err = -1
    DO istep = 1, max_steps
       CALL cp_iterate(logger%iter_info,iter_nr=istep)
       neb_env%opt_type_label="SD"
       ! Save the optimization step counter
       neb_env%istep = istep
       ! Perform one step of SD with line search
       norm = SQRT(SUM(forces%wrk*forces%wrk))
       IF (norm<EPSILON(0.0_dp)) THEN
          ! Let's handle the case in which the band is already fully optimized
          converged = .TRUE.
          EXIT
       END IF
       sline%wrk = forces%wrk / norm
       IF (do_ls.AND.(.NOT.skip_ls)) THEN
          CALL neb_ls(stepsize, sline, rep_env, neb_env, coords, energies, forces,&
               vels, particle_set, iw, output_unit, distances, diis_section, iw2)
          IF (iw2>0) &
               WRITE(iw2,'(T2,A,T69,F12.6)')"SD| Stepsize in SD after linesearch",&
               stepsize
       ELSE
          stepsize = MIN(norm*stepsize0,max_stepsize)
          IF (iw2>0) &
               WRITE(iw2,'(T2,A,T69,F12.6)')"SD| Stepsize in SD no linesearch performed",&
               stepsize
       END IF
       sline%wrk = stepsize*sline%wrk
       diis_on = accept_diis_step(istep>max_sd_steps,n_diis,err,crr,set_err,sline,coords,&
                                  check_diis,iw2)
       IF (diis_on) THEN
          neb_env%opt_type_label= "DIIS"
       END IF
       do_ls = .TRUE.
       IF (COUNT(set_err==-1)==1) do_ls = .FALSE.
       coords%wrk = coords%wrk + sline%wrk
       ! Compute forces
       CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
            output_unit, distances, neb_env%number_of_replica)
       CALL reparametrize_images(neb_env%reparametrize_frames, neb_env%spline_order,&
            neb_env%smoothing, coords%wrk, sline%wrk, distances)
       neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
       CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
            particle_set, iw)
       ! Check for an external exit command
       CALL external_control(should_stop,"NEB",globenv=globenv)
       IF(should_stop) EXIT
       ! Dump Infos
       CALL dump_neb_info(neb_env=neb_env,&
                          coords=coords,&
                          forces=forces,&
                          particle_set=particle_set,&
                          logger=logger,&
                          istep=istep,&
                          energies=energies,&
                          distances=distances,&
                          vels=vels,&
                          output_unit=output_unit)

       converged = check_convergence(neb_env, sline, forces)
       IF (converged) EXIT
    END DO
    DEALLOCATE (energies)
    DEALLOCATE (distances)
    DEALLOCATE (err)
    DEALLOCATE (crr)
    DEALLOCATE (set_err)
    CALL neb_var_release(sline)
    CALL timestop(handle)
  END SUBROUTINE neb_diis

END MODULE neb_methods
