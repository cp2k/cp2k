!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/neb_methods [1.0] *
!!
!!   NAME
!!     neb_methods
!!
!!   FUNCTION
!!     Module performing a Nudged Elastic Band Calculation
!!
!!   NOTES
!!     Numerical accuracy for parallel runs:
!!      Each replica starts the SCF run from the one optimized
!!      in a previous run. It may happen then energies and derivatives
!!      of a serial run and a parallel run could be slightly different
!!      'cause of a different starting density matrix.
!!      Exact results are obtained using:
!!         EXTRAPOLATION USE_GUESS in QS section (Teo 09.2006)
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
MODULE neb_methods
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: neb_diis_opt,&
                                             neb_md_opt
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE neb_IO,                          ONLY: dump_neb_info,&
                                             read_neb_section
  USE neb_md_utils,                    ONLY: control_vels_a,&
                                             control_vels_b
  USE neb_opt_utils,                   ONLY: accept_diis_step,&
                                             neb_ls
  USE neb_types,                       ONLY: neb_type
  USE neb_utils,                       ONLY: build_replica_coords,&
                                             check_convergence,&
                                             neb_calc_energy_forces,&
                                             reorient_images
  USE particle_types,                  ONLY: particle_type
  USE replica_types,                   ONLY: rep_env_create,&
                                             rep_env_release,&
                                             rep_env_write,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'neb_methods'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
  PUBLIC :: neb

CONTAINS

!!****f* neb_methods/neb [1.0] *
!!
!!   NAME
!!     neb
!!
!!   FUNCTION
!!     Real subroutine for NEB calculations
!!
!!   NOTES
!!     Based on the use of replica_env
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE neb(input, para_env, globenv, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ierr, ihandle, &
                                                irep, nrep, output_unit, &
                                                prep, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords, vels
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(neb_type), POINTER                  :: neb_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER :: diis_section, force_env_section, &
      md_section, motion_section, neb_section, print_section

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(logger, subsys, f_env, rep_env)
       NULLIFY(coords, vels, neb_env)
       logger => cp_error_get_logger(error)
       CALL cp_add_iter_level(logger%iter_info,"NEB",ihandle,error)
       motion_section => section_vals_get_subs_vals(input,"MOTION",error=error)
       print_section  => section_vals_get_subs_vals(motion_section,"PRINT",error=error)
       neb_section    => section_vals_get_subs_vals(motion_section,"NEB",error=error)
       output_unit=cp_print_key_unit_nr(logger,motion_section,"NEB%PROGRAM_RUN_INFO",&
            extension=".nebLog",error=error)

       CALL section_vals_val_get(neb_section,"NPROC_REP",i_val=prep, error=error)
       nrep = MAX(1,para_env%num_pe/prep)
       IF (output_unit>0) THEN
          WRITE(output_unit,'(T2,"NEB|",3(A,I5))')" NUMBER OF REPLICAS:",nrep,&
               " NUMBER OF PROCS:",para_env%num_pe," NUMBER OF PROCS PER REPLICA:",prep
       END IF
       IF (nrep*prep/=para_env%num_pe.AND.output_unit>0) THEN
          WRITE(output_unit,'(T2,"NEB|",A,I5)')" NUMBER OF TOTALLY REQUESTED PROCESSORS:",&
               para_env%num_pe," IS NOT COMPATIBLE WITH THE NUMBER OF PROCESSORS REQUESTED PER REPLICA:",&
               prep," AND THE NUMBER OF REPLICAS:",nrep
          WRITE(output_unit,'(T2,"NEB| ",I0,A)')para_env%num_pe-nrep*prep," PROCESSORS WILL BE WASTED!"
       END IF
       ! Just one force_env allowed
       force_env_section => section_vals_get_subs_vals(input,"FORCE_EVAL",error=error)
       CALL section_vals_get(force_env_section,n_repetition=irep,error=error)
       CPPostcondition(irep==1,cp_failure_level,routineP,error,failure)
       ! Create Replica Environments
       IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Replica_env Setup. START"
       CALL rep_env_create(rep_env, para_env=para_env, input=input,&
            nrep=nrep,prep=prep,error=error)
       IF (debug_this_module) THEN
          CALL rep_env_write(rep_env,cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               error=error)
       END IF
       IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Replica_env Setup. END"
       IF (ASSOCIATED(rep_env)) THEN
          CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
               new_error=new_error, failure=failure)
          CALL force_env_get(f_env%force_env,subsys=subsys,error=error)
          CPPostcondition(SIZE(subsys)==1,cp_failure_level,routineP,error,failure)
          particle_set => subsys(1)%subsys%particles%els
          ! Read NEB controlling parameters
          ALLOCATE(neb_env, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          neb_env%force_env            => f_env%force_env
          neb_env%root_section         => input
          neb_env%force_env_section    => force_env_section
          neb_env%motion_print_section => print_section
          neb_env%neb_section          => neb_section
          CALL read_neb_section(neb_env,neb_section,error=error)
          ! Collecting the coordinates of the starting replicas of the NEB
          IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Building initial set of coordinates. START"
          CALL build_replica_coords(neb_section, particle_set,&
               coords, vels, neb_env, rep_env%ndim, output_unit, globenv, &
               rep_env%para_env, logger, error)
          IF (output_unit>0) WRITE(output_unit,'(T2,"NEB|",A)')" Building initial set of coordinates. END"
          ! Perform NEB optimization
          SELECT CASE(neb_env%opt_type)
          CASE(neb_md_opt)
             neb_env%opt_type_label = "MOLECULAR DYNAMICS"
             md_section => section_vals_get_subs_vals(neb_section,"MD",error=error)
             CALL neb_md(rep_env, neb_env, coords, vels, particle_set, output_unit,&
                  md_section, logger, error)
          CASE(neb_diis_opt)
             neb_env%opt_type_label = "DIIS"
             diis_section => section_vals_get_subs_vals(neb_section,"DIIS",error=error)
             CALL neb_diis(rep_env, neb_env, coords, vels, particle_set, output_unit,&
                  diis_section, logger, error)
          END SELECT
          CALL f_env_rm_defaults(f_env,new_error,ierr)
          DEALLOCATE(neb_env, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (ASSOCIATED(coords)) THEN
             DEALLOCATE(coords, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(vels)) THEN
             DEALLOCATE(vels, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END IF
       CALL rep_env_release(rep_env,error=error)
       CALL cp_print_key_finished_output(output_unit,logger,motion_section,&
            "NEB%PROGRAM_RUN_INFO", error=error)
       CALL cp_rm_iter_level(logger%iter_info,ihandle,error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE neb

!!****f* neb_methods/neb_md [1.0] *
!!
!!   NAME
!!     neb_md
!!
!!   FUNCTION
!!     MD type optimization NEB
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE neb_md(rep_env, neb_env, coords, vels, particle_set, output_unit,&
       md_section, logger, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), OPTIONAL, POINTER        :: neb_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords, vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_md', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ind, iparticle, &
                                                istep, iw, k, max_steps, stat
    LOGICAL                                  :: converged, failure
    REAL(KIND=dp)                            :: dm, dt, my_mass
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distances, energies
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mass
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Dcoords, forces
    TYPE(section_vals_type), POINTER         :: tc_section, vc_section

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(forces, tc_section, vc_section)
       CPPostcondition(ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(vels),cp_failure_level,routineP,error,failure)
       iw=cp_print_key_unit_nr(logger,neb_env%neb_section,"REPLICA_INFO",&
            extension=".replicaLog",error=error)
       tc_section => section_vals_get_subs_vals(md_section,"TEMP_CONTROL",error=error)
       vc_section => section_vals_get_subs_vals(md_section,"VEL_CONTROL",error=error)
       CALL section_vals_val_get(md_section,"TIMESTEP",r_val=dt, error=error)
       CALL section_vals_val_get(md_section,"MAX_STEPS",i_val=max_steps, error=error)
       ! Initial setup for MD
       ALLOCATE (forces(3*SIZE(particle_set),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (Dcoords(3*SIZE(particle_set),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (mass(3*SIZE(particle_set),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (energies(SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (distances(SIZE(coords,2)-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! setting up the mass array
       DO iparticle = 1, SIZE(particle_set)
          my_mass = particle_set(iparticle)%atomic_kind%mass
          dm = 0.5_dp*dt/my_mass
          DO k = 1, 3
             ind = (iparticle-1)*3 + k
             mass(ind,1) = dm
          END DO
       END DO
       DO k = 2, neb_env%number_of_replica
          mass(:,k) = mass(:,k-1)
       END DO
       ! Initializing forces array
       CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
            output_unit, distances, rep_env%para_env, error=error)
       neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
       CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
            particle_set, iw, error)

       CALL dump_neb_info(neb_env=neb_env,&
                          coords=coords,&
                          vels=vels,&
                          forces=forces,&
                          particle_set=particle_set,&
                          logger=logger,&
                          istep=0,&
                          energies=energies,&
                          distances=distances,&
                          output_unit=output_unit,&
                          error=error)
       md_opt_loop: DO istep = 1, max_steps
          CALL cp_iterate(logger%iter_info,iter_nr=istep,error=error)
          ! Save the optimization step counter
          neb_env%istep = istep
          ! Velocity Verlet (first part)
          vels(:,:)   = vels(:,:)   + mass(:,:) * forces(:,:)
          ! Control on velocity - I part [rescale, annealing]
          CALL control_vels_a(vels, particle_set, tc_section, vc_section, output_unit,&
               istep, error)
          ! Coordinate step
          Dcoords(:,:) = dt * vels(:,:)
          coords (:,:) = coords(:,:) + Dcoords(:,:)

          CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
               output_unit, distances, rep_env%para_env, error=error)
          neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
          CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
               particle_set, iw, error)
          ! Control on velocity - II part [check vels VS forces, Steepest Descent like]
          CALL control_vels_b(vels, forces, vc_section, error)
          ! Velocity Verlet (second part)
          vels(:,:)   = vels(:,:) + mass(:,:) * forces(:,:)
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
                             output_unit=output_unit,&
                             error=error)
          converged = check_convergence(neb_env, Dcoords, forces, error)
          IF (converged) EXIT
       END DO md_opt_loop
       DEALLOCATE (Dcoords,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (forces,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (mass,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (energies,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (distances,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL cp_print_key_finished_output(iw,logger,neb_env%neb_section,&
            "REPLICA_INFO", error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE neb_md

!!****f* neb_methods/neb_diis [1.0] *
!!
!!   NAME
!!     neb_diis
!!
!!   FUNCTION
!!     DIIS type optimization NEB
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     Teodoro Laino 09.2006
!!
!!   MODIFICATION HISTORY
!!             
!!
!!   SOURCE
!****************************************************************************
  SUBROUTINE neb_diis(rep_env, neb_env, coords, vels, particle_set, output_unit,&
       diis_section, logger, error)
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(neb_type), OPTIONAL, POINTER        :: neb_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: coords, vels
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(section_vals_type), POINTER         :: diis_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'neb_diis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istep, iw, &
                                                max_sd_steps, max_steps, &
                                                n_diis, stat
    INTEGER, DIMENSION(:), POINTER           :: set_err
    LOGICAL                                  :: converged, diis_on, do_ls, &
                                                failure, skip_ls, check_diis
    REAL(KIND=dp)                            :: norm, stepsize, stepsize0, max_stepsize
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: distances, energies
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: crr, err, forces, sline

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       neb_env%opt_type_label = "SD"
       do_ls = .TRUE.
       vels  = 0.0_dp
       NULLIFY(forces)
       CPPostcondition(ASSOCIATED(coords),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(vels),cp_failure_level,routineP,error,failure)
       iw=cp_print_key_unit_nr(logger,neb_env%neb_section,"REPLICA_INFO",&
            extension=".replicaLog",error=error)
       CALL section_vals_val_get(diis_section,"MAX_STEPS",i_val=max_steps, error=error)
       CALL section_vals_val_get(diis_section,"N_DIIS",i_val=n_diis, error=error)
       CALL section_vals_val_get(diis_section,"STEPSIZE",r_val=stepsize0, error=error)
       CALL section_vals_val_get(diis_section,"MAX_STEPSIZE",r_val=max_stepsize, error=error)
       CALL section_vals_val_get(diis_section,"NO_LS",l_val=skip_ls, error=error)
       CALL section_vals_val_get(diis_section,"MAX_SD_STEPS",i_val=max_sd_steps, error=error)
       CALL section_vals_val_get(diis_section,"CHECK_DIIS",l_val=check_diis, error=error)
       stepsize = stepsize0
       ! Initial setup for DIIS
       ALLOCATE (forces(3*SIZE(particle_set),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! line search direction
       ALLOCATE (sline(3*SIZE(particle_set),SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Array of error vectors
       ALLOCATE (err(3*SIZE(particle_set)*SIZE(coords,2),n_diis),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (crr(3*SIZE(particle_set)*SIZE(coords,2),n_diis),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (set_err(n_diis),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (energies(SIZE(coords,2)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (distances(SIZE(coords,2)-1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Initializing forces array
       CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
            output_unit, distances, rep_env%para_env, error=error)
       neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
       CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
            particle_set, iw, error)
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
                          output_unit=output_unit,&
                          error=error)
       ! If rotation is requested let's apply it at the beginning of the 
       ! Geometry optimization and then let's disable it
       neb_env%rotate_frames= .FALSE.
       ! Main SD/DIIS loop
       set_err = -1
       DO istep = 1, max_steps
          CALL cp_iterate(logger%iter_info,iter_nr=istep,error=error)
          neb_env%opt_type_label="SD"
          ! Perform one step of SD with line search
          norm = SQRT(SUM(forces*forces))
          IF (norm<EPSILON(0.0_dp)) THEN
             ! Let's handle the case in which the band is already fully optimized
             converged = .TRUE.
             EXIT
          END IF
          sline = forces / norm
          IF (do_ls.AND.(.NOT.skip_ls)) THEN
             CALL neb_ls(stepsize, sline, rep_env, neb_env, coords, energies, forces,&
                  vels, particle_set, iw, output_unit, distances, diis_section, error)
             IF (output_unit>0) &
                  WRITE(output_unit,'(T2,A,T69,F12.6)')"SD| Stepsize in SD after linesearch",&
                  stepsize
          ELSE
             stepsize = min(norm*stepsize0,max_stepsize)
             IF (output_unit>0) &
                  WRITE(output_unit,'(T2,A,T69,F12.6)')"SD| Stepsize in SD no linesearch performed",&
                  stepsize
          END IF
          sline = stepsize*sline
          diis_on = accept_diis_step(istep>max_sd_steps,n_diis,err,crr,set_err,sline,coords,&
                                     output_unit,check_diis,error)
          IF (diis_on) THEN 
             neb_env%opt_type_label= "DIIS"
          END IF
          do_ls = .TRUE.
          IF (count(set_err==-1)==1) do_ls = .FALSE.
          coords = coords + sline
          ! Compute forces
          CALL reorient_images(neb_env%rotate_frames, particle_set, coords, vels,&
               output_unit, distances, rep_env%para_env, error)
          neb_env%avg_distance = SQRT(SUM(distances*distances)/REAL(SIZE(distances),KIND=dp))
          CALL neb_calc_energy_forces(rep_env, neb_env, coords, energies, forces,&
               particle_set, iw, error)
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
                             output_unit=output_unit,&
                             error=error)

          converged = check_convergence(neb_env, sline, forces, error)
          IF (converged) EXIT
       END DO
       DEALLOCATE (energies, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (distances, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (forces, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (sline, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (err, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (crr, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (set_err, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE neb_diis

END MODULE neb_methods
