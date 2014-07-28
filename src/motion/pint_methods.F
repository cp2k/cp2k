!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods to performs a path integral run
!> \author fawzi
!> \par History
!>      02.2005 created [fawzi]
!>           11.2006 modified so it might actually work [hforbert]
!> \note   quick & dirty rewrite of my python program
! *****************************************************************************
MODULE pint_methods

  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_external_control,             ONLY: external_control
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k,&
                                             cp_unit_to_cp2k
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get
  USE gle_system_dynamics,             ONLY: gle_cholesky_stab,&
                                             gle_matrix_exp,&
                                             restart_gle
  USE gle_system_types,                ONLY: gle_dealloc,&
                                             gle_init,&
                                             gle_thermo_create
  USE global_types,                    ONLY: global_environment_type
  USE helium_methods,                  ONLY: helium_create,&
                                             helium_init,&
                                             helium_release
  USE helium_sampling,                 ONLY: helium_do_run,&
                                             helium_step
  USE helium_types,                    ONLY: helium_solvent_type
  USE input_constants,                 ONLY: transformation_normal,&
                                             transformation_stage
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: &
       section_type, section_vals_get, section_vals_get_subs_vals, &
       section_vals_release, section_vals_retain, section_vals_type, &
       section_vals_val_get, section_vals_val_set, section_vals_val_unset
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE machine,                         ONLY: m_walltime
  USE mathconstants,                   ONLY: twopi
  USE mathlib,                         ONLY: gcd
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             next_rng_seed,&
                                             rng_stream_type
  USE particle_list_types,             ONLY: particle_list_type
  USE pint_gle,                        ONLY: pint_calc_gle_energy,&
                                             pint_gle_init,&
                                             pint_gle_step
  USE pint_io,                         ONLY: pint_write_centroids,&
                                             pint_write_com,&
                                             pint_write_ener,&
                                             pint_write_line,&
                                             pint_write_rgyr,&
                                             pint_write_step_info,&
                                             pint_write_trajectory
  USE pint_normalmode,                 ONLY: normalmode_calc_uf_h,&
                                             normalmode_env_create,&
                                             normalmode_init_masses,&
                                             normalmode_release
  USE pint_public,                     ONLY: pint_com_pos,&
                                             pint_levy_walk
  USE pint_staging,                    ONLY: staging_calc_uf_h,&
                                             staging_env_create,&
                                             staging_init_masses,&
                                             staging_release
  USE pint_transformations,            ONLY: pint_f2uf,&
                                             pint_u2x,&
                                             pint_x2u
  USE pint_types,                      ONLY: e_conserved_id,&
                                             e_kin_thermo_id,&
                                             e_kin_virial_id,&
                                             e_potential_id,&
                                             pint_env_type,&
                                             thermostat_gle,&
                                             thermostat_none,&
                                             thermostat_nose
  USE replica_methods,                 ONLY: rep_env_calc_e_f,&
                                             rep_env_create
  USE replica_types,                   ONLY: rep_env_release,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_methods'
  INTEGER, SAVE, PRIVATE :: last_pint_id=0

  PUBLIC :: do_pint_run

CONTAINS

  ! ***************************************************************************
  !> \brief  Create a path integral environment
  !> \author fawzi
  !> \note   Might return an unassociated pointer in parallel on the processors
  !>         that are not needed.
  !> \par    History
  !>           Fixed some bugs [hforbert]
  !>           Added normal mode transformation [hforbert]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param input ...
!> \param input_declaration ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_create(pint_env,input,input_declaration,para_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: output_file_name, project_name
    INTEGER                                  :: handle, iat, idim, idir, &
                                                ierr, itmp, nrep, prep, stat
    LOGICAL                                  :: explicit, failure, ltmp, &
                                                wrong_input
    REAL(kind=dp)                            :: mass
    REAL(kind=dp), DIMENSION(3, 2)           :: seed
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER         :: gle_section, nose_section, &
                                                pint_section, &
                                                transform_section

    CALL timeset(routineN,handle)

    failure=.FALSE.
    NULLIFY(f_env,subsys,particles, nose_section, gle_section)

    CPPrecondition(.NOT.ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
    CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(rep_env)
       pint_section => section_vals_get_subs_vals(input,"MOTION%PINT",&
            error=error)
       CALL section_vals_val_get(pint_section,"p",i_val=nrep,error=error)
       CALL section_vals_val_get(pint_section,"proc_per_replica",&
            i_val=prep,error=error)
       ! Maybe let the user have his/her way as long as prep is
       ! within the bounds of number of CPUs??
       IF ( (prep < 1) .OR. (prep > para_env%num_pe) .OR. &
            (MOD(prep*nrep,para_env%num_pe) /= 0) ) THEN
          prep = para_env%num_pe/gcd(para_env%num_pe,nrep)
          IF (para_env%ionode) THEN
             WRITE (*,*) "PINT WARNING: Adjusting number of processors per replica to ",prep
          END IF
       END IF

       ! replica_env modifies the global input structure - which is wrong - one
       ! of the side effects is the inifite adding of the -r-N string to the
       ! project name and the output file name, which corrupts restart files.
       ! For now: save the project name and output file name and restore them
       ! after the rep_env_create has executed - the initialization of the
       ! replicas will run correctly anyways.
       ! TODO: modify rep_env so that it behaves better
       CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",c_val=project_name,error=error)
       CALL section_vals_val_get(input,"GLOBAL%OUTPUT_FILE_NAME",c_val=output_file_name,error=error)
       CALL rep_env_create(rep_env, para_env=para_env, input=input,&
            input_declaration=input_declaration,nrep=nrep,prep=prep, row_force=.TRUE.,error=error)
       CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",c_val=TRIM(project_name),error=error)
       IF ( LEN_TRIM(output_file_name) .GT. 0 ) THEN
         CALL section_vals_val_set(input,"GLOBAL%OUTPUT_FILE_NAME",c_val=TRIM(output_file_name),error=error)
       ELSE
         CALL section_vals_val_unset(input,"GLOBAL%OUTPUT_FILE_NAME",error=error)
       END IF
       IF (.NOT. ASSOCIATED(rep_env)) RETURN

       ALLOCATE(pint_env,STAT=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
    END IF
    IF (.NOT. failure) THEN
       last_pint_id=last_pint_id+1
       pint_env%id_nr=last_pint_id
       pint_env%ref_count=1
       NULLIFY(pint_env%replicas,pint_env%input,pint_env%staging_env,&
               pint_env%normalmode_env)
       pint_env%p=nrep
       pint_env%replicas => rep_env
       pint_env%ndim=rep_env%ndim
       pint_env%input => input
       CALL section_vals_retain(pint_env%input,error=error)

       ! get first step, last step, number of steps, etc
       CALL section_vals_val_get(input,"MOTION%PINT%ITERATION",&
            i_val=itmp, error=error)
       pint_env%first_step = itmp
       CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
            explicit=explicit, error=error)
       IF ( explicit ) THEN
         CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
              i_val=itmp, error=error)
         pint_env%last_step = itmp
         pint_env%num_steps = pint_env%last_step - pint_env%first_step
       ELSE
         CALL section_vals_val_get(input,"MOTION%PINT%NUM_STEPS",&
              i_val=itmp, error=error)
         pint_env%num_steps = itmp
         pint_env%last_step = pint_env%first_step + pint_env%num_steps
       END IF

       CALL section_vals_val_get(pint_section,"DT",&
            r_val=pint_env%dt,error=error)
       pint_env%t = pint_env%first_step * pint_env%dt

       CALL section_vals_val_get(pint_section,"nrespa",i_val=pint_env%nrespa,&
            error=error)
       CALL section_vals_val_get(pint_section,"Temp",r_val=pint_env%kT,&
            error=error)
       CALL section_vals_val_get(pint_section,"T_TOL",&
            r_val=pint_env%t_tol,error=error)
       CALL section_vals_val_get(pint_section,"transformation",&
            i_val=pint_env%transform, error=error)

       NULLIFY(pint_env%tx,pint_env%tv,pint_env%tv_t,pint_env%tv_old,pint_env%tv_new,pint_env%tf)

       pint_env%nnos = 0
       pint_env%pimd_thermostat = thermostat_none
       nose_section => section_vals_get_subs_vals(input,"MOTION%PINT%NOSE",&
            error=error)
       CALL section_vals_get(nose_section, explicit=explicit, error=error)
       IF(explicit) THEN
         CALL section_vals_val_get(nose_section,"nnos",i_val=pint_env%nnos,&
              error=error)
         IF(pint_env%nnos>0)THEN
           pint_env%pimd_thermostat = thermostat_nose
           ALLOCATE(&
              pint_env%tx(pint_env%nnos,pint_env%p,pint_env%ndim),    &
              pint_env%tv(pint_env%nnos,pint_env%p,pint_env%ndim),    &
              pint_env%tv_t(pint_env%nnos,pint_env%p,pint_env%ndim),  &
              pint_env%tv_old(pint_env%nnos,pint_env%p,pint_env%ndim),&
              pint_env%tv_new(pint_env%nnos,pint_env%p,pint_env%ndim),&
              pint_env%tf(pint_env%nnos,pint_env%p,pint_env%ndim), STAT=stat)
           CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
           pint_env%tx       = 0._dp
           pint_env%tv       = 0._dp
           pint_env%tv_t     = 0._dp
           pint_env%tv_old   = 0._dp
           pint_env%tv_new   = 0._dp
           pint_env%tf       = 0._dp
         END IF
       END IF

       pint_env%beta=1._dp/pint_env%kT
!TODO
! v_tol not in current input structure
! should also probably be part of nose_section
!       CALL section_vals_val_get(transform_section,"v_tol_nose",r_val=pint_env%v_tol,&
!            error=error)
!MK ... but we have to initialise v_tol
       pint_env%v_tol = 0.0_dp ! to be fixed

       NULLIFY (pint_env%randomG)

       seed(:,:) = next_rng_seed(error=error)
       CALL create_rng_stream(pint_env%randomG,&
                              name="pint_randomG",&
                              distribution_type=GAUSSIAN,&
                              extended_precision=.TRUE.,&
                              seed=seed,error=error)

       ALLOCATE(pint_env%e_pot_bead(pint_env%p),STAT=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       pint_env%e_pot_bead=0._dp
       pint_env%e_pot_h=0._dp
       pint_env%e_kin_beads=0._dp
       pint_env%e_pot_t=0._dp
       pint_env%e_gle=0._dp
       pint_env%e_kin_t=0._dp
       pint_env%energy(:) = 0.0_dp

!TODO: rearrange to use standard nose hoover chain functions/data types

       ALLOCATE(&
            pint_env%x(pint_env%p,pint_env%ndim),          &
            pint_env%v(pint_env%p,pint_env%ndim),          &
            pint_env%f(pint_env%p,pint_env%ndim),          &
            pint_env%external_f(pint_env%p,pint_env%ndim), &
            pint_env%ux(pint_env%p,pint_env%ndim),         &
            pint_env%uv(pint_env%p,pint_env%ndim),         &
            pint_env%uv_t(pint_env%p,pint_env%ndim),       &
            pint_env%uv_new(pint_env%p,pint_env%ndim),     &
            pint_env%uf(pint_env%p,pint_env%ndim),         &
            pint_env%uf_h(pint_env%p,pint_env%ndim),       &
            pint_env%rtmp_ndim(pint_env%ndim),             &
            pint_env%rtmp_natom(pint_env%ndim/3),          &
            STAT=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       pint_env%x        = 0._dp
       pint_env%v        = 0._dp
       pint_env%f        = 0._dp
       pint_env%external_f= 0._dp
       pint_env%ux       = 0._dp
       pint_env%uv       = 0._dp
       pint_env%uv_t     = 0._dp
       pint_env%uv_new   = 0._dp
       pint_env%uf       = 0._dp
       pint_env%uf_h     = 0._dp
       pint_env%rtmp_ndim= 0._dp
       pint_env%rtmp_natom= 0._dp
       pint_env%time_per_step = 0.0_dp

       IF (pint_env%transform == transformation_stage) THEN
          transform_section => section_vals_get_subs_vals(input,&
             "MOTION%PINT%STAGING",error=error)
          CALL staging_env_create(pint_env%staging_env,transform_section,&
             p=pint_env%p,kT=pint_env%kT, error=error)
       ELSE
          transform_section => section_vals_get_subs_vals(input,&
             "MOTION%PINT%NORMALMODE",error=error)
          CALL normalmode_env_create(pint_env%normalmode_env,&
             transform_section,p=pint_env%p,kT=pint_env%kT,error=error)
          wrong_input=pint_env%nrespa*twopi/(SQRT(pint_env%p/pint_env%normalmode_env%modefactor)*pint_env%kT)/pint_env%dt>10
          CPPostcondition(wrong_input,cp_warning_level,routineP,error,failure)
       END IF
       ALLOCATE(pint_env%mass(pint_env%ndim),STAT=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
            f_env=f_env,new_error=new_error, failure=failure)
       CALL force_env_get(force_env=f_env%force_env,subsys=subsys,&
            error=new_error)
       CALL cp_subsys_get(subsys,particles=particles,error=new_error)

!TODO length of pint_env%mass is redundant
       idim=0
       DO iat=1,pint_env%ndim/3
          CALL get_atomic_kind(particles%els(iat)%atomic_kind,mass=mass)
          DO idir=1,3
             idim=idim+1
             pint_env%mass(idim)=mass
          END DO
       END DO
       CALL f_env_rm_defaults(f_env,new_error,ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

       ALLOCATE(pint_env%Q(pint_env%p),&
            pint_env%mass_beads(pint_env%p,pint_env%ndim),&
            pint_env%mass_fict(pint_env%p,pint_env%ndim),STAT=stat)
       CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
       IF (pint_env%transform == transformation_stage) THEN
          CALL staging_init_masses(pint_env%staging_env,mass=pint_env%mass,&
               mass_beads=pint_env%mass_beads,mass_fict=pint_env%mass_fict,&
               Q=pint_env%Q,error=error)
       ELSE
          CALL normalmode_init_masses(pint_env%normalmode_env, &
               mass=pint_env%mass, mass_beads=pint_env%mass_beads, &
               mass_fict=pint_env%mass_fict, Q=pint_env%Q, error=error)
       END IF

       NULLIFY(pint_env%gle)
       gle_section => section_vals_get_subs_vals(input,"MOTION%PINT%GLE",&
            error=error)
       CALL section_vals_get(gle_section, explicit=explicit, error=error)
       IF(explicit) THEN
         ALLOCATE(pint_env%gle,stat=stat)
         CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
         CALL gle_init(pint_env%gle, dt=pint_env%dt/pint_env%nrespa, temp=pint_env%kT,&
              section=gle_section, error=error)
         IF (pint_env%pimd_thermostat==thermostat_none .AND. pint_env%gle%ndim .GT. 0) THEN
            pint_env%pimd_thermostat=thermostat_gle

            ! initialize a GLE with ALL degrees of freedom on node 0,
            ! as it seems to me that here everything but force eval is replicated
            pint_env%gle%loc_num_gle=pint_env%p*pint_env%ndim
            pint_env%gle%glob_num_gle=pint_env%gle%loc_num_gle
            ALLOCATE(pint_env%gle%map_info%index(pint_env%gle%loc_num_gle))
            CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
            DO itmp = 1, pint_env%gle%loc_num_gle
                pint_env%gle%map_info%index(itmp)=itmp
            ENDDO
            CALL gle_thermo_create(pint_env%gle,pint_env%gle%loc_num_gle,error)

            ! here we should have read a_mat and c_mat;
            !we can therefore compute the matrices needed for the propagator
            ! deterministic part of the propagator
            CALL gle_matrix_exp((-pint_env%dt/pint_env%nrespa*0.5_dp)*pint_env%gle%a_mat, &
                  pint_env%gle%ndim,15,15,pint_env%gle%gle_t)
            ! stochastic part
            CALL gle_cholesky_stab(pint_env%gle%c_mat-MATMUL(pint_env%gle%gle_t,&
                  MATMUL(pint_env%gle%c_mat,TRANSPOSE(pint_env%gle%gle_t))), &
                                pint_env%gle%gle_s, pint_env%gle%ndim)
           ! and initialize the additional momenta
            CALL pint_gle_init(pint_env, error)
         END IF
       END IF

       CALL section_vals_val_get(pint_section,"FIX_CENTROID_POS",&
            l_val=ltmp,error=error)
       IF ( ltmp .AND. (pint_env%transform .EQ. transformation_normal) ) THEN
         pint_env%first_propagated_mode = 2
       ELSE
         pint_env%first_propagated_mode = 1
       END IF

    END IF

    CALL timestop(handle)

    RETURN
  END SUBROUTINE pint_create

  ! ***************************************************************************
  !> \brief Retain a path integral environment
  !> \param pint_env the pint_env to retain
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author Fawzi Mohamed
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_retain(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
      pint_env%ref_count=pint_env%ref_count+1
    END IF
    RETURN
  END SUBROUTINE pint_retain

  ! ***************************************************************************
  !> \brief Release a path integral environment
  !> \param pint_env the pint_env to release
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author Fawzi Mohamed
  !> \par History
  !>      Added normal mode transformation [hforbert]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_release(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(pint_env)) THEN
      CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
      pint_env%ref_count=pint_env%ref_count-1
      IF (pint_env%ref_count==0) THEN
        CALL rep_env_release(pint_env%replicas,error=error)
        CALL section_vals_release(pint_env%input,error=error)
        IF (ASSOCIATED(pint_env%staging_env)) THEN
           CALL staging_release(pint_env%staging_env,error=error)
        END IF
        IF (ASSOCIATED(pint_env%normalmode_env)) THEN
           CALL normalmode_release(pint_env%normalmode_env,error=error)
        END IF
        CALL delete_rng_stream(pint_env%randomG,error=error)

        DEALLOCATE(pint_env%mass,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%e_pot_bead,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DEALLOCATE(pint_env%x,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%v,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%f,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%external_f,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%mass_beads,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%mass_fict,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%ux,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uv,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uv_t,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uv_new,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uf,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%uf_h,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%rtmp_ndim,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        DEALLOCATE(pint_env%rtmp_natom,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        IF(pint_env%pimd_thermostat==thermostat_nose) THEN
          DEALLOCATE(pint_env%tx,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pint_env%tv,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pint_env%tv_t,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pint_env%tv_old,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pint_env%tv_new,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DEALLOCATE(pint_env%tf,STAT=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ELSEIF(pint_env%pimd_thermostat==thermostat_gle) THEN
          CALL gle_dealloc(pint_env%gle, error=error)
        END IF

        DEALLOCATE(pint_env%Q,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

        DEALLOCATE(pint_env,STAT=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

      END IF
    END IF

    NULLIFY(pint_env)

    RETURN
  END SUBROUTINE pint_release

  ! ***************************************************************************
  !> \brief Tests the path integral methods
  !> \param para_env parallel environment
  !> \param input the input to test
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param para_env ...
!> \param input ...
!> \param input_declaration ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_test(para_env,input,input_declaration,error)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ib, idim, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: c, e_h, err
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: x1
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pint_env_type), POINTER             :: pint_env

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
  CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
  NULLIFY(pint_env)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
      CALL pint_create(pint_env,input,input_declaration,para_env,error=error)
     IF (ASSOCIATED(pint_env)) THEN
        ALLOCATE(x1(pint_env%ndim,pint_env%p),stat=stat)
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
        x1(:,:)=pint_env%x
        CALL pint_x2u(pint_env,error=error)
        pint_env%x=0._dp
        CALL pint_u2x(pint_env,error=error)
        err=0._dp
        DO i=1,pint_env%ndim
           err=MAX(err,ABS(x1(1,i)-pint_env%x(1,i)))
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_r1="//cp_to_string(err),&
             local=.FALSE.)

        CALL pint_calc_uf_h(pint_env,e_h=e_h,error=error)
        c=-pint_env%staging_env%w_p**2
        pint_env%f=0._dp
        DO idim=1,pint_env%ndim
           DO ib=1,pint_env%p
              pint_env%f(ib,idim)=pint_env%f(ib,idim)+&
                   c*(2._dp*pint_env%x(ib,idim)&
                   -pint_env%x(MODULO(ib-2,pint_env%p)+1,idim)&
                   -pint_env%x(MODULO(ib,pint_env%p)+1,idim))
           END DO
        END DO
        CALL pint_f2uf(pint_env,error=error)
        err=0._dp
        DO idim=1,pint_env%ndim
           DO ib=1,pint_env%p
              err=MAX(err,ABS(pint_env%uf(ib,idim)-pint_env%uf_h(ib,idim)))
           END DO
        END DO
        CALL cp_log(logger,cp_note_level+1,routineP,"diff_f_h="//cp_to_string(err),&
             local=.FALSE.)
     END IF
    END IF
    RETURN
  END SUBROUTINE pint_test

  ! ***************************************************************************
  !> \brief  Perform a path integral simulation
  !> \param  para_env parallel environment
  !> \param  input the input to test
  !> \param  error variable to control error logging, stopping,...
  !>         see module cp_error_handling
  !> \author Fawzi Mohamed
  !> \par    History
  !>         2003-11 created [fawzi]
  !>         2009-12-14 globenv parameter added to handle soft exit
  !>           requests [lwalewski]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param para_env ...
!> \param input ...
!> \param input_declaration ...
!> \param globenv ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE do_pint_run( para_env, input, input_declaration, globenv, error )
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_pint_run', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: helium_only_mid = 1, &
                                                solute_only_mid = 2, &
                                                solute_with_helium_mid = 3

    INTEGER                                  :: handle, mode
    LOGICAL                                  :: explicit, failure, &
                                                helium_only, solvent_present
    TYPE(helium_solvent_type), POINTER       :: helium
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(section_vals_type), POINTER         :: helium_section

    CALL timeset(routineN,handle)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(input),cp_failure_level,routineP,error,failure)
    CPPrecondition(para_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(input%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN

      ! check if helium solvent is present
      NULLIFY(helium_section)
      helium_section => section_vals_get_subs_vals(input,&
                        "MOTION%PINT%HELIUM",error=error)
      CALL section_vals_get(helium_section,explicit=explicit,error=error)
      IF ( explicit ) THEN
        CALL section_vals_val_get(helium_section,"_SECTION_PARAMETERS_",&
             l_val=solvent_present,error=error)
      ELSE
        solvent_present = .FALSE.
      END IF

      ! check if there is anything but helium
      IF (solvent_present) THEN
        CALL section_vals_val_get(helium_section,"HELIUM_ONLY",&
          l_val=helium_only,error=error)
      ELSE
        helium_only = .FALSE.
      END IF

      ! pick the mode of operation
      mode = 0
      IF (solvent_present) THEN
        IF (helium_only) THEN
          mode = helium_only_mid
        ELSE
         mode = solute_with_helium_mid
        END IF
      ELSE
        mode = solute_only_mid
      END IF

      NULLIFY(pint_env,helium)

      ! perform the simulation according to the chosen mode
      SELECT CASE (mode)

      CASE (helium_only_mid)
        CALL helium_create(helium,input,error=error)
        CALL helium_init(helium,pint_env,error)
        CALL helium_do_run(helium,globenv,error)
        CALL helium_release(helium,error=error)

      CASE (solute_only_mid)
        CALL pint_create(pint_env,input,input_declaration,para_env,error=error)
        CALL pint_init(pint_env,error)
        CALL pint_do_run(pint_env,globenv,error=error)
        CALL pint_release(pint_env,error=error)

      CASE (solute_with_helium_mid)
        CALL pint_create(pint_env,input,input_declaration,para_env,error=error)
        ! init pint wihtout helium forces (they are not yet initialized)
        CALL pint_init(pint_env,error)
        ! init helium with solute's positions (they are already initialized)
        CALL helium_create(helium,input,solute=pint_env,error=error)
        CALL helium_init(helium,pint_env,error)
        ! reinit pint forces with helium forces (they are now initialized)
        CALL pint_init_f( pint_env, helium_solvent=helium, error=error )

        CALL pint_do_run(pint_env,globenv,helium=helium,error=error)
        CALL helium_release(helium,error=error)
        CALL pint_release(pint_env,error=error)

      CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
          routineP,"Unknown mode ("//TRIM(ADJUSTL(cp_to_string(mode)))//")",&
          error,failure)
      END SELECT

    END IF

    CALL timestop(handle)

    RETURN
  END SUBROUTINE do_pint_run

  ! ***************************************************************************
  !> \brief  Reads the restart, initializes the beads, etc.
  !> \author Fawzi Mohamed
  !> \par    History
  !>           11.2003 created [fawzi]
  !>           actually ASSIGN input pointer [hforbert]
  !>           2010-12-16 turned into a wrapper routine [lwalewski]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_init(pint_env, error)

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( failure ) THEN
      RETURN
    END IF

    CALL pint_init_x( pint_env, error )
    CALL pint_init_v( pint_env, error=error )
    CALL pint_init_t( pint_env, error=error )
    CALL pint_init_f( pint_env, error=error )

    RETURN
  END SUBROUTINE pint_init


  ! ***************************************************************************
  !> \brief  Assign initial postions to the beads.
  !> \author Lukasz Walewski
  !> \date   2010-12-15
  !> \descr  Initialization is done in the following way:
  !>           1. assign all beads with the same classical positions from
  !>              FORCE_EVAL (hot start)
  !>           2. spread the beads around classical positions as if they were
  !>              free particles (if requested)
  !>           3. replace positions generated in steps 1-2 with the explicit
  !>              ones if they are explicitly given in the input structure
  !>           4. apply Gaussian noise to the positions generated so far (if
  !>              requested)
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_init_x( pint_env, error )

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_x', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msg, tmp
    INTEGER                                  :: ia, ib, ic, idim, input_seed, &
                                                istat, n_rep_val
    LOGICAL                                  :: done_init, done_levy, &
                                                done_rand, explicit, failure, &
                                                levycorr, ltmp
    REAL(kind=dp)                            :: tcorr, var
    REAL(kind=dp), DIMENSION(3)              :: x0
    REAL(kind=dp), DIMENSION(3, 2)           :: seed
    REAL(kind=dp), DIMENSION(:), POINTER     :: bx, r_vals
    TYPE(rng_stream_type), POINTER           :: rng_gaussian
    TYPE(section_vals_type), POINTER         :: input_section

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( failure ) THEN
      RETURN
    END IF

    DO idim = 1, pint_env%ndim
      DO ib=1,pint_env%p
        pint_env%x(ib,idim)=pint_env%replicas%r(idim,ib)
      END DO
    END DO

    done_levy = .FALSE.
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%LEVY_POS_SAMPLE",&
         l_val=ltmp,error=error)
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%LEVY_TEMP_FACTOR",&
         r_val=tcorr,error=error)
    IF ( ltmp ) THEN

      NULLIFY(bx)
      ALLOCATE(bx(3*pint_env%p), STAT=istat)
      CPPostcondition(istat==0,cp_fatal_level,routineP,error,failure)
      NULLIFY (rng_gaussian)
      CALL section_vals_val_get(pint_env%input,&
           "MOTION%PINT%INIT%LEVY_SEED",i_val=input_seed,error=error)
      seed(:,:) = REAL(input_seed,KIND=dp)
!      seed(:,:) = next_rng_seed(error=error)
      CALL create_rng_stream(rng_gaussian,&
                             name="tmp_rng_gaussian",&
                             distribution_type=GAUSSIAN,&
                             extended_precision=.TRUE.,&
                             seed=seed,error=error)

      CALL section_vals_val_get(pint_env%input,&
           "MOTION%PINT%INIT%LEVY_CORRELATED",&
           l_val=levycorr,error=error)

      IF ( levycorr ) THEN

        ! correlated Levy walk - the same path for all atoms
        x0 = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        CALL pint_levy_walk( x0, pint_env%p, 1.0_dp, bx, rng_gaussian, error )
        idim = 0
        DO ia = 1, pint_env%ndim/3
          var = SQRT(1.0_dp/(pint_env%kT*tcorr*pint_env%mass(3*ia)))
          DO ic = 1, 3
            idim = idim + 1
            DO ib=1,pint_env%p
              pint_env%x(ib,idim) = pint_env%x(ib,idim) + bx(3*(ib-1)+ic)*var
            END DO
          END DO
        END DO

      ELSE

        ! uncorrelated bead initialization - distinct Levy walk for each atom
        idim = 0
        DO ia = 1, pint_env%ndim/3
          x0(1) = pint_env%x(1,3*(ia-1)+1)
          x0(2) = pint_env%x(1,3*(ia-1)+2)
          x0(3) = pint_env%x(1,3*(ia-1)+3)
          var = SQRT(1.0_dp/(pint_env%kT*tcorr*pint_env%mass(3*ia)))
          CALL pint_levy_walk( x0, pint_env%p, var, bx, rng_gaussian, error )
          DO ic = 1, 3
            idim = idim + 1
            DO ib=1,pint_env%p
              pint_env%x(ib,idim) = pint_env%x(ib,idim) + bx(3*(ib-1)+ic)
            END DO
          END DO
        END DO

      END IF

      CALL delete_rng_stream(rng_gaussian,error=error)
      DEALLOCATE(bx, STAT=istat)
      CPPostcondition(istat==0,cp_fatal_level,routineP,error,failure)
      done_levy = .TRUE.
    END IF

    done_init = .FALSE.
    NULLIFY(input_section)
    input_section => section_vals_get_subs_vals(pint_env%input,&
                     "MOTION%PINT%BEADS%COORD",&
                     error=error)
    CALL section_vals_get(input_section,explicit=explicit,error=error)
    IF (explicit) THEN
      CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
           n_rep_val=n_rep_val,error=error)
      IF (n_rep_val>0) THEN
        CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
             r_vals=r_vals,error=error)
        CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim,&
             cp_failure_level,cp_assertion_failed,&
             "Invalid size of MOTION%PINT%BEADS%COORD "//&
             CPSourceFileRef,&
             routineP,error,failure)
        ic=0
        DO idim=1,pint_env%ndim
          DO ib=1,pint_env%p
            ic=ic+1
            pint_env%x(ib,idim)=r_vals(ic)
          END DO
        END DO
        done_init = .TRUE.
      END IF
    END IF

    done_rand = .FALSE.
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%RANDOMIZE_POS",&
         l_val=ltmp,error=error)
    IF ( ltmp ) THEN
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          pint_env%x(ib,idim) = pint_env%x(ib,idim) + &
                                next_random_number(rng_stream=pint_env%randomG,&
                                variance=pint_env%beta/&
                                SQRT(12.0_dp*pint_env%mass(idim)),&
                                error=error)
        END DO
      END DO
      done_rand = .TRUE.
    END IF

    WRITE(tmp,'(A)') "Bead positions initialization:"
    IF ( done_init ) THEN
      WRITE(msg,'(A,A)') TRIM(tmp), " input structure"
    ELSE IF ( done_levy ) THEN
      WRITE(msg,'(A,A)') TRIM(tmp), " Levy random walk"
    ELSE
      WRITE(msg,'(A,A)') TRIM(tmp), " hot start"
    END IF
    CALL pint_write_line(msg, error)

    IF ( done_levy ) THEN
      WRITE(msg,'(A,F6.3)') "Levy walk at effective temperature: ", tcorr
    END IF

    IF ( done_rand ) THEN
      WRITE(msg,'(A)') "Added gaussian noise to the positions of the beads."
      CALL pint_write_line(msg, error)
    END IF

    RETURN
  END SUBROUTINE pint_init_x


  ! ***************************************************************************
  !> \brief  Initialize velocities
  !> \param  pint_env the pint env in which you should initialize the
  !>         velocity
  !> \param  error variable to control error logging, stopping,...
  !>         see module cp_error_handling
  !> \par    History
  !>         2010-12-16 gathered all velocity-init code here [lwalewski]
  !>         2011-04-05 added centroid velocity initialization [lwalewski]
  !>         2011-12-19 removed optional parameter kT, target temperature is
  !>                    now determined from the input directly [lwalewski]
  !> \note   Initialization is done according to the following protocol:
  !>         1. set all the velocities to FORCE_EVAL%SUBSYS%VELOCITY if present
  !>         2. scale the velocities according to the actual temperature
  !>            (has no effect if vels not present in 1.)
  !>         3. draw vels for the remaining dof from MB distribution
  !>            (all or non-centroid modes only depending on 1.)
  !>         4. add random noise to the centroid vels if CENTROID_SPEED == T
  !>         5. set the vels for all dof to 0.0 if VELOCITY_QUENCH == T
  !>         6. set the vels according to the explicit values from the input
  !>            if present
  !> \author fawzi
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_init_v( pint_env, error )
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_v', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: msg, stmp, stmp1, stmp2, &
                                                unit_str
    INTEGER                                  :: first_mode, ia, ib, ic, idim, &
                                                itmp, n_rep_val
    LOGICAL                                  :: done_init, done_quench, &
                                                done_scale, done_sped, &
                                                explicit, failure, ltmp, &
                                                vels_present
    REAL(kind=dp)                            :: actual_t, ek, rtmp, target_t, &
                                                unit_conv
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_vals
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: input_section

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( failure ) THEN
      RETURN
    END IF

    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    ! read the velocities from the input file if they are given explicitly
    vels_present = .FALSE.
    NULLIFY(input_section)
    input_section => section_vals_get_subs_vals(pint_env%input,&
                     "FORCE_EVAL%SUBSYS%VELOCITY",&
                     error=error)
    CALL section_vals_get(input_section,explicit=explicit,error=error)
    IF ( explicit ) THEN

      CALL section_vals_val_get(input_section,"PINT_UNIT",&
                                c_val=unit_str,error=error)
      unit_conv = cp_unit_to_cp2k(1.0_dp,TRIM(unit_str),error=error)

      ! assign all the beads with the same velocities from FORCE_EVAL%SUBSYS%VELOCITY
      NULLIFY(r_vals)
      CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
           n_rep_val=n_rep_val,error=error)
      stmp = ""
      WRITE(stmp,*) n_rep_val
      msg = "Invalid number of atoms in FORCE_EVAL%SUBSYS%VELOCITY ("//&
            TRIM(ADJUSTL(stmp))//")."
      CALL cp_assert(3*n_rep_val==pint_env%ndim,&
                     cp_failure_level, cp_assertion_failed,&
                     routineP//" from "//&
CPSourceFileRef,&
                     msg, error, failure)
      DO ia = 1, pint_env%ndim/3
        CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_",&
             i_rep_val=ia, r_vals=r_vals, error=error)
        itmp = SIZE(r_vals)
        stmp = ""
        WRITE(stmp,*) itmp
        msg = "Number of coordinates != 3 in FORCE_EVAL%SUBSYS%VELOCITY ("//&
              TRIM(ADJUSTL(stmp))//")."
        CALL cp_assert(itmp==3,&
                       cp_failure_level, cp_assertion_failed,&
                       routineP//" from "//&
CPSourceFileRef,&
                       msg, error, failure)
        DO ib = 1, pint_env%p
          DO ic = 1, 3
          idim = 3*(ia-1)+ic
          pint_env%v(ib,idim) = r_vals(ic) * unit_conv
          END DO
        END DO
      END DO

      vels_present = .TRUE.
    END IF

    ! set the actual temperature...
    IF (vels_present) THEN
      ! ...from the initial velocities
      ek = 0.0_dp
      DO ia = 1, pint_env%ndim/3
        rtmp = 0.0_dp
        DO ic = 1, 3
          idim = 3*(ia-1)+ic
          rtmp = rtmp + pint_env%v(1,idim) * pint_env%v(1,idim)
        END DO
        ek = ek + 0.5_dp * pint_env%mass(idim) * rtmp
      END DO
      actual_t = 2.0_dp * ek / pint_env%ndim
    ELSE
      ! ...using the temperature value from the input
      actual_t = pint_env%kT
    END IF

    ! set the target temperature
    target_t = pint_env%kT
    CALL section_vals_val_get(pint_env%input, &
         "MOTION%PINT%INIT%VELOCITY_SCALE",&
         l_val=done_scale,error=error)
    IF (vels_present) THEN
      IF (done_scale) THEN
        ! rescale the velocities to match the target temperature
        rtmp = SQRT( target_t / actual_t )
        DO ia = 1, pint_env%ndim/3
          DO ib = 1, pint_env%p
            DO ic = 1, 3
              idim = 3*(ia-1)+ic
              pint_env%v(ib,idim) = rtmp * pint_env%v(ib,idim)
            END DO
          END DO
        END DO
      ELSE
        target_t = actual_t
      END IF
    END IF

    ! draw velocities from the M-B distribution...
    IF (vels_present) THEN
      ! ...for non-centroid modes only
      CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
      first_mode = 2
    ELSE
      ! ...for all the modes
      first_mode = 1
    END IF
    DO idim = 1, SIZE(pint_env%uv,2)
      DO ib = first_mode, SIZE(pint_env%uv,1)
        pint_env%uv(ib,idim) = &
          next_random_number(rng_stream=pint_env%randomG,&
                             variance=target_t/pint_env%mass_fict(ib,idim),&
                             error=error)
      END DO
    END DO

    ! add random component to the centroid velocity if requsted
    done_sped = .FALSE.
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%CENTROID_SPEED",&
         l_val=ltmp,error=error)
    IF (ltmp) THEN
      CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
      DO idim=1,pint_env%ndim
        rtmp = next_random_number(rng_stream=pint_env%randomG,&
                                variance=pint_env%mass(idim)*pint_env%kT,&
                                error=error)/pint_env%mass(idim)
        DO ib=1,pint_env%p
          pint_env%v(ib,idim)=pint_env%v(ib,idim)+rtmp
        END DO
      END DO
      CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
      done_sped = .TRUE.
    END IF

    ! quench (set to zero) velocities for all the modes if requested
    ! (disregard all the initialization done so far)
    done_quench = .FALSE.
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%VELOCITY_QUENCH",&
         l_val=ltmp,error=error)
    IF (ltmp) THEN
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          pint_env%v(ib,idim) = 0.0_dp
        END DO
      END DO
      CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
      done_quench = .TRUE.
    END IF

    ! set the velocities to the values from the input if they are explicit
    ! (disregard all the initialization done so far)
    done_init = .FALSE.
    NULLIFY(input_section)
    input_section => section_vals_get_subs_vals(pint_env%input,&
                     "MOTION%PINT%BEADS%VELOCITY",&
                     error=error)
    CALL section_vals_get(input_section,explicit=explicit,error=error)
    IF (explicit) THEN
      CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
           n_rep_val=n_rep_val,error=error)
      IF (n_rep_val>0) THEN
        CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
             r_vals=r_vals,error=error)
        CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim,&
                       cp_failure_level,cp_assertion_failed,&
                       "Invalid size of MOTION%PINT%BEAD%VELOCITY "//&
                       CPSourceFileRef,&
                       routineP,error,failure)
        itmp=0
        DO idim=1,pint_env%ndim
          DO ib=1,pint_env%p
            itmp=itmp+1
            pint_env%v(ib,idim)=r_vals(itmp)
          END DO
        END DO
        CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v,error=error)
        done_init = .TRUE.
      END IF
    END IF

    unit_conv = cp_unit_from_cp2k(1.0_dp,"K",error=error)
    WRITE(stmp1,'(F10.2)') target_t * unit_conv
    msg = "Bead velocities initialization:"
    IF ( done_init ) THEN
      msg = TRIM(msg)//" input structure"
    ELSE IF ( done_quench ) THEN
      msg = TRIM(msg)//" quenching (set to 0.0)"
    ELSE
      IF ( vels_present ) THEN
        msg = TRIM(ADJUSTL(msg))//" centroid +"
      END IF
      msg = TRIM(ADJUSTL(msg))//" Maxwell-Boltzmann at "//TRIM(ADJUSTL(stmp1))//" K."
    END IF
    CALL pint_write_line(msg, error)

    IF ( done_init .AND. done_quench ) THEN
      msg = "WARNING: exclusive options requested (velocity restart and quenching)"
      CALL pint_write_line(msg, error)
      msg = "WARNING: velocity restart took precedence"
      CALL pint_write_line(msg, error)
    END IF

    IF ( (.NOT. done_init) .AND. (.NOT. done_quench) ) THEN
      IF (vels_present .AND. done_scale) THEN
        WRITE(stmp1,'(F10.2)') actual_t * unit_conv
        WRITE(stmp2,'(F10.2)') target_t * unit_conv
        msg = "Scaled initial velocities from "//TRIM(ADJUSTL(stmp1))//&
              " to "//TRIM(ADJUSTL(stmp2))//" K as requested."
        CALL pint_write_line(msg, error)
      END IF
      IF ( done_sped ) THEN
        msg = "Added random component to the initial centroid velocities."
        CALL pint_write_line(msg, error)
      END IF
    END IF

    RETURN
  END SUBROUTINE pint_init_v


  ! ***************************************************************************
  !> \brief  Assign initial postions and velocities to the thermostats.
  !> \author Lukasz Walewski
  !> \date   2010-12-15
  !> \note   Extracted from pint_init
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param kT ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_init_t( pint_env, kT, error )

    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), INTENT(in), OPTIONAL      :: kT
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_t', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, ii, inos, n_rep_val
    LOGICAL                                  :: explicit, failure, gle_restart
    REAL(kind=dp)                            :: mykt
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_vals
    TYPE(section_vals_type), POINTER         :: input_section

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( failure ) THEN
      RETURN
    END IF

    IF(pint_env%pimd_thermostat==thermostat_nose) THEN

      mykt=pint_env%kT
      IF (PRESENT(kT)) mykt=kT
      DO idim=1,SIZE(pint_env%tv,3)
        DO ib=1,SIZE(pint_env%tv,2)
          DO inos=1,SIZE(pint_env%tv,1)
             pint_env%tv(inos,ib,idim) = &
               next_random_number(rng_stream=pint_env%randomG,&
                                  variance=mykt/pint_env%Q(ib),&
                                  error=error)
          END DO
        END DO
      END DO

      NULLIFY(input_section)
      input_section => section_vals_get_subs_vals(pint_env%input,&
                       "MOTION%PINT%NOSE%COORD",&
                       error=error)
      CALL section_vals_get(input_section,explicit=explicit,error=error)
      IF (explicit) THEN
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                  n_rep_val=n_rep_val,error=error)
        IF (n_rep_val>0) THEN
          CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
          CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                    r_vals=r_vals,error=error)
          CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim*pint_env%nnos,&
                         cp_failure_level,cp_assertion_failed,&
                         "Invalid size of MOTION%PINT%NOSE%COORD "//&
                         CPSourceFileRef,&
                         routineP,error,failure)
          ii=0
          DO idim=1,pint_env%ndim
            DO ib=1,pint_env%p
              DO inos=1,pint_env%nnos
                ii=ii+1
                pint_env%tx(inos,ib,idim)=r_vals(ii)
              END DO
            END DO
          END DO
        END IF
      END IF

      NULLIFY(input_section)
      input_section => section_vals_get_subs_vals(pint_env%input,&
                       "MOTION%PINT%NOSE%VELOCITY",&
                       error=error)
      CALL section_vals_get(input_section,explicit=explicit,error=error)
      IF (explicit) THEN
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                  n_rep_val=n_rep_val,error=error)
        IF (n_rep_val>0) THEN
          CPPrecondition(n_rep_val==1,cp_failure_level,routineP,error,failure)
          CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                    r_vals=r_vals,error=error)
          CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim*pint_env%nnos,&
                         cp_failure_level,cp_assertion_failed,&
                         "Invalid size of MOTION%PINT%NOSE%VELOCITY "//&
                         CPSourceFileRef,&
                         routineP,error,failure)
          ii=0
          DO idim=1,pint_env%ndim
            DO ib=1,pint_env%p
              DO inos=1,pint_env%nnos
                ii=ii+1
                pint_env%tv(inos,ib,idim)=r_vals(ii)
              END DO
            END DO
          END DO
        END IF
      END IF

    ELSEIF(pint_env%pimd_thermostat==thermostat_gle) THEN
      NULLIFY(input_section)
      input_section => section_vals_get_subs_vals(pint_env%input,&
                       "MOTION%PINT%GLE",&
                       error=error)
      CALL section_vals_get(input_section,explicit=explicit,error=error)
      IF (explicit) THEN
        CALL restart_gle(pint_env%gle,input_section,save_mem=.FALSE.,&
                         restart=gle_restart,error=error)
      END IF
    END IF

    RETURN
  END SUBROUTINE pint_init_t


  ! ***************************************************************************
  !> \brief  Prepares the forces, etc. to perform an MD step
  !> \author fawzi
  !> \par    History
  !>           Added nh_energy calculation [hforbert]
  !>           Bug fixes for no thermostats [hforbert]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium_solvent ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_init_f( pint_env, helium_solvent, error )
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_solvent
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: e_h

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL pint_x2u(pint_env,error=error)
      CALL pint_calc_uf_h(pint_env=pint_env,e_h=e_h,error=error)
      CALL pint_calc_f(pint_env,error=error)

      ! add helium forces to the solute's internal ones
      ! Assume that helium has been already initialized and helium_solvent
      ! contains proper forces in force_avrg array.
      IF (PRESENT(helium_solvent)) THEN
        IF (ASSOCIATED(helium_solvent)) THEN
          pint_env%f(:,:) = pint_env%f(:,:) + helium_solvent%force_avrg(:,:)
        END IF
      END IF
      CALL pint_f2uf(pint_env,error=error)

      ! set the centroid forces to 0 if FIX_CENTROID_POS
      IF ( pint_env%first_propagated_mode .EQ. 2 ) THEN
        pint_env%uf(1,:) = 0.0_dp
      END IF

      CALL pint_calc_e_kin_beads_u(pint_env,error=error)
      CALL pint_calc_e_vir(pint_env,error=error)
      DO idim=1,SIZE(pint_env%uf_h,2)
        DO ib=pint_env%first_propagated_mode,SIZE(pint_env%uf_h,1)
          pint_env%uf_h(ib,idim)=pint_env%uf_h(ib,idim)&
               +REAL(pint_env%nrespa,dp)*pint_env%uf(ib,idim)
        END DO
      END DO

      IF (pint_env%nnos > 0) THEN
        DO idim=1,SIZE(pint_env%uf_h,2)
          DO ib=1,SIZE(pint_env%uf_h,1)
            pint_env%tf(1,ib,idim)=(pint_env%mass_fict(ib,idim)*&
                 pint_env%uv(ib,idim)**2-pint_env%kT)/pint_env%Q(ib)
          END DO
        END DO

        DO idim=1,pint_env%ndim
          DO ib=1,pint_env%p
            DO inos=1,pint_env%nnos-1
              pint_env%tf(inos+1,ib,idim)=pint_env%tv(inos,ib,idim)**2-&
                pint_env%kT/pint_env%Q(ib)
            END DO
            DO inos=1,pint_env%nnos-1
              pint_env%tf(inos,ib,idim)=pint_env%tf(inos,ib,idim)&
                -pint_env%tv(inos,ib,idim)*pint_env%tv(inos+1,ib,idim)
            END DO
          END DO
        END DO
        CALL pint_calc_nh_energy(pint_env,error=error)
      END IF
    END IF
    RETURN
  END SUBROUTINE pint_init_f


  ! ***************************************************************************
  !> \brief  Perform the PIMD simulation (main MD loop)
  !> \author Fawzi Mohamed
  !> \note   Everything should be read for an md step.
  !> \par    History
  !>         2003-11 created [fawzi]
  !>         renamed from pint_run to pint_do_run because of conflicting name
  !>           of pint_run in input_constants [hforbert]
  !>         2009-12-14 globenv parameter added to handle soft exit
  !>           requests [lwalewski]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param globenv ...
!> \param helium ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_do_run(pint_env, globenv, helium, error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_do_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: step
    LOGICAL                                  :: failure, should_stop
    REAL(kind=dp)                            :: scal

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

    NULLIFY(pint_env%logger)
    pint_env%logger => cp_error_get_logger(error)

    ! create iteration level and initialize iteration info
    CALL cp_add_iter_level(pint_env%logger%iter_info,"MD",error=error)
    CALL cp_iterate(pint_env%logger%iter_info,iter_nr=pint_env%first_step,error=error)
    pint_env%iter = pint_env%first_step

    ! write the properties at 0-th step
    CALL pint_calc_energy(pint_env,error)
    CALL pint_write_ener(pint_env,error=error)
    CALL pint_write_centroids(pint_env,error=error)
    CALL pint_write_trajectory(pint_env,error=error)
    CALL pint_write_com(pint_env,error=error)
    CALL pint_write_rgyr(pint_env,error=error)

    ! main PIMD loop
    DO step = 1, pint_env%num_steps

      pint_env%iter = pint_env%iter + 1
      CALL cp_iterate(pint_env%logger%iter_info,&
                      last=(step==pint_env%num_steps),&
                      iter_nr=pint_env%iter,error=error)
      pint_env%t = pint_env%t + pint_env%dt

      IF (pint_env%t_tol > 0.0_dp) THEN
        IF (ABS(2._dp*pint_env%e_kin_beads/(pint_env%p*pint_env%ndim)&
          -pint_env%kT)>pint_env%t_tol) THEN
          scal=SQRT(pint_env%kT*(pint_env%p*pint_env%ndim)/(2.0_dp*pint_env%e_kin_beads))
          pint_env%uv=scal*pint_env%uv
          CALL pint_init_f(pint_env,helium_solvent=helium,error=error)
        END IF
      END IF
      CALL pint_step(pint_env,helium_solvent=helium,error=error)

      CALL pint_write_ener(pint_env,error=error)
      CALL pint_write_centroids(pint_env,error=error)
      CALL pint_write_trajectory(pint_env,error=error)
      CALL pint_write_com(pint_env,error=error)
      CALL pint_write_rgyr(pint_env,error=error)

      CALL write_restart(root_section=pint_env%input,&
           pint_env=pint_env, helium_env=helium, error=error)

      ! exit from the main loop if soft exit has been requested
      CALL external_control(should_stop,"MD",globenv=globenv,error=error)
      IF (should_stop) EXIT

    END DO

    ! remove iteration level
    CALL cp_rm_iter_level(pint_env%logger%iter_info,"MD",error=error)

    RETURN
  END SUBROUTINE pint_do_run


  ! ***************************************************************************
  !> \brief  Does an MD step (and nrespa harmonic evaluations)
  !> \author fawzi
  !> \par    History
  !>           various bug fixes [hforbert]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param helium_solvent ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_step(pint_env,helium_solvent,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_solvent
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ib, idim, inos, &
                                                iresp
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dti, dti2, dti22, e_h, rn, &
                                                tdti, time_start, time_stop, &
                                                tol
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: tmp

    CALL timeset(routineN,handle)
    time_start = m_walltime()

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     rn=REAL(pint_env%nrespa,dp)
     dti=pint_env%dt/rn
     dti2=dti/2._dp
     tdti=2.*dti
     dti22=dti**2/2._dp

     DO iresp=1,pint_env%nrespa
        ! integrate bead positions, first_propagated_mode = { 1, 2 }
        IF (pint_env%pimd_thermostat==thermostat_nose) THEN
          DO i = pint_env%first_propagated_mode, pint_env%p
            pint_env%ux(i,:) = pint_env%ux(i,:) + &
                               dti*pint_env%uv(i,:) + &
                               dti22*(pint_env%uf_h(i,:) - &
                               pint_env%uv(i,:)*pint_env%tv(1,i,:))
          END DO
          pint_env%tx=pint_env%tx+dti*pint_env%tv+dti22*pint_env%tf
        ELSE
           DO i = pint_env%first_propagated_mode, pint_env%p
             pint_env%ux(i,:) = pint_env%ux(i,:) + &
                                dti*pint_env%uv(i,:) + &
                                dti22*pint_env%uf_h(i,:)
           END DO
        END IF
        ! integrate nh vars
        ! integrate v at half step
        IF (pint_env%pimd_thermostat==thermostat_nose) THEN
           pint_env%uv_t=pint_env%uv+dti2* &
                         (pint_env%uf_h-pint_env%uv*pint_env%tv(1,:,:))
           tmp => pint_env%tv_t
           pint_env%tv_t => pint_env%tv
           pint_env%tv => tmp

           pint_env%tv=pint_env%tv_old+tdti*pint_env%tf
           pint_env%tv_old=pint_env%tv_t
           pint_env%tv_t=pint_env%tv_t+dti2*pint_env%tf
        ELSE
           pint_env%uv_t=pint_env%uv+dti2*pint_env%uf_h
        END IF


        ! calc forces at new pos
        CALL pint_calc_uf_h(pint_env=pint_env,e_h=e_h,error=error)
        IF (iresp==pint_env%nrespa) THEN
           CALL pint_u2x(pint_env,error=error)
           CALL pint_calc_f(pint_env,error=error)

           ! perform helium step and add helium forces
           IF (PRESENT(helium_solvent)) THEN
             IF (ASSOCIATED(helium_solvent)) THEN
               helium_solvent%current_step = pint_env%iter
               helium_solvent%origin = pint_com_pos(pint_env,error)
               CALL helium_step( helium_solvent, pint_env, error )
               pint_env%f(:,:)=pint_env%f(:,:)+helium_solvent%force_avrg(:,:)
             END IF
           END IF

           CALL pint_f2uf(pint_env,error=error)
           pint_env%uf_h=pint_env%uf_h+rn*pint_env%uf
        END IF

        ! set the centroid forces to 0 if FIX_CENTROID_POS
        IF ( pint_env%first_propagated_mode .EQ. 2 ) THEN
          pint_env%uf_h(1,:) = 0.0_dp
        END IF

        ! add the new forces to v
        pint_env%uv_t=pint_env%uv_t+dti2*pint_env%uf_h
        ! iterate v to convergence
        IF (pint_env%pimd_thermostat==thermostat_nose) THEN
           DO i=1,6
              tol=0._dp
              pint_env%uv_new=pint_env%uv_t/(1.+dti2*pint_env%tv(1,:,:))
              DO idim=1,pint_env%ndim
                 DO ib=1,pint_env%p
                    pint_env%tf(1,ib,idim)=(pint_env%mass_fict(ib,idim)*&
                         pint_env%uv_new(ib,idim)**2-pint_env%kT)/&
                         pint_env%Q(ib)
                 END DO
              END DO

              DO idim=1,pint_env%ndim
                 DO ib=1,pint_env%p
                    DO inos=1,pint_env%nnos-1
                       pint_env%tv_new(inos,ib,idim)=&
                            (pint_env%tv_t(inos,ib,idim)+dti2*pint_env%tf(inos,ib,idim))/&
                            (1._dp+dti2*pint_env%tv(inos+1,ib,idim))
                       pint_env%tf(inos+1,ib,idim)=&
                            (pint_env%tv_new(inos,ib,idim)**2-&
                            pint_env%kT/pint_env%Q(ib))
                       tol=MAX(tol,ABS(pint_env%tv(inos,ib,idim)&
                            -pint_env%tv_new(inos,ib,idim)))
                    END DO
                    pint_env%tv_new(pint_env%nnos,ib,idim)=&
                         pint_env%tv_t(pint_env%nnos,ib,idim)+&
                         dti2*pint_env%tf(pint_env%nnos,ib,idim)
                    tol=MAX(tol,ABS(pint_env%tv(pint_env%nnos,ib,idim)&
                         -pint_env%tv_new(pint_env%nnos,ib,idim)))
                    tol=MAX(tol,ABS(pint_env%uv(ib,idim)&
                         -pint_env%uv_new(ib,idim)))
                 END DO
              END DO

              pint_env%uv=pint_env%uv_new
              pint_env%tv=pint_env%tv_new
              IF (tol <= pint_env%v_tol) EXIT
           END DO
           DO inos=1,pint_env%nnos-1
              pint_env%tf(inos,:,:)=pint_env%tf(inos,:,:)&
                   -pint_env%tv(inos,:,:)*pint_env%tv(inos+1,:,:)
           END DO
        ELSEIF(pint_env%pimd_thermostat==thermostat_gle) THEN
           CALL pint_gle_step(pint_env, error=error)
           pint_env%uv=pint_env%uv_t
        ELSE
           pint_env%uv=pint_env%uv_t
        ENDIF
     END DO

     ! calculate the energy components
     CALL pint_calc_energy( pint_env, error )

     ! check that the number of MD steps matches
     ! the number of force evaluations done so far
!TODO make this check valid if we start from ITERATION != 0
!     CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
!          f_env=f_env,new_error=new_error, failure=failure)
!     NULLIFY(logger)
!     logger => cp_error_get_logger(new_error)
!     CALL cp_assert(logger%iter_info%iteration(2)==pint_env%iter+1,&
!          cp_failure_level,cp_assertion_failed,routineP,&
!          "md & force_eval lost sychro "//&
!          CPSourceFileRef,&
!          error,failure)
!     CALL f_env_rm_defaults(f_env,new_error,ierr)

  END IF

    time_stop = m_walltime()
    pint_env%time_per_step = time_stop - time_start
    CALL pint_write_step_info(pint_env,error)
    CALL timestop(handle)

    RETURN
  END SUBROUTINE pint_step


  ! ***************************************************************************
  !> \brief  Calculate the energy components (private wrapper function)
  !> \author Lukasz Walewski
  !> \date   2011-01-07
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_calc_energy( pint_env, error )

    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_energy', &
      routineP = moduleN//':'//routineN

    CALL pint_calc_e_kin_beads_u(pint_env,error=error)
    CALL pint_calc_e_vir(pint_env,error=error)

    IF (pint_env%pimd_thermostat==thermostat_nose) THEN
      CALL pint_calc_nh_energy(pint_env,error=error)
    ELSEIF (pint_env%pimd_thermostat==thermostat_gle) THEN
      CALL pint_calc_gle_energy(pint_env,error=error)
    END IF

    pint_env%energy(e_kin_thermo_id) = &
      0.5_dp * REAL(pint_env%p,dp) * REAL(pint_env%ndim,dp) * pint_env%kT - &
      pint_env%e_pot_h

    pint_env%energy(e_potential_id) = &
      SUM(pint_env%e_pot_bead) / REAL(pint_env%p,dp)

    pint_env%energy(e_conserved_id) = &
      pint_env%energy(e_potential_id) + &
      pint_env%e_pot_h + &
      pint_env%e_kin_beads + &
      pint_env%e_pot_t + &
      pint_env%e_kin_t + pint_env%e_gle

    RETURN
  END SUBROUTINE pint_calc_energy

  ! ***************************************************************************
  !> \brief  Calculate the harmonic force in the u basis
  !> \param  pint_env the path integral environment in which the harmonic
  !>         forces should be calculated
  !> \param  error variable to control error logging, stopping,...
  !>         see module cp_error_handling
  !> \author fawzi
  !> \par    History
  !>           Added normal mode transformation [hforbert]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param e_h ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_calc_uf_h(pint_env,e_h,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(KIND=dp), INTENT(OUT)               :: e_h
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_uf_h', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT. failure) THEN
      IF (pint_env%transform == transformation_stage) THEN
        CALL staging_calc_uf_h(pint_env%staging_env,&
                               pint_env%mass_beads,&
                               pint_env%ux,&
                               pint_env%uf_h,&
                               pint_env%e_pot_h,&
                               error=error)
      ELSE
        CALL normalmode_calc_uf_h(pint_env%normalmode_env,&
                                  pint_env%mass_beads,&
                                  pint_env%ux,&
                                  pint_env%uf_h,&
                                  pint_env%e_pot_h,&
                                  error=error)
      END IF
      e_h=pint_env%e_pot_h
      pint_env%uf_h=pint_env%uf_h/pint_env%mass_fict
    END IF
    RETURN
  END SUBROUTINE pint_calc_uf_h

  ! ***************************************************************************
  !> \brief calculates the force (and energy) in each bead, returns the sum
  !>      of the potential energy
  !> \param pint_env path integral environment on which you want to calculate
  !>        the forces
  !> \param x positions at which you want to evaluate the forces
  !> \param f the forces
  !> \param e potential energy on each bead
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  !> \par    History
  !>           2009-06-15 moved helium calls out from here [lwalewski]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param x ...
!> \param f ...
!> \param e ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_calc_f(pint_env,x,f,e,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: x
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: f
    REAL(kind=dp), DIMENSION(:), &
      INTENT(out), OPTIONAL, TARGET          :: e
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:), POINTER     :: my_e
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_f, my_x

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      my_x => pint_env%x
      IF (PRESENT(x)) my_x => x
      my_f => pint_env%f
      IF (PRESENT(f)) my_f => f
      my_e => pint_env%e_pot_bead
      IF (PRESENT(e)) my_e => e
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          pint_env%replicas%r(idim,ib)=my_x(ib,idim)
        END DO
      END DO
      CALL rep_env_calc_e_f(pint_env%replicas,calc_f=.TRUE.,error=error)
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          !ljw: is that fine ? - idim <-> ib
          my_f(ib,idim) = pint_env%replicas%f(idim,ib)
        END DO
      END DO
      my_e=pint_env%replicas%f(SIZE(pint_env%replicas%f,1),:)

    END IF
    RETURN
  END SUBROUTINE pint_calc_f

  ! ***************************************************************************
  !> \brief  Calculate the kinetic energy of the beads (in the u variables)
  !> \param  error variable to control error logging, stopping,...
  !>         see module cp_error_handling
  !> \author fawzi
  !> \par    History
  !>         Bug fix to give my_uv a default location if not given in call [hforbert]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param uv ...
!> \param e_k ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_calc_e_kin_beads_u(pint_env,uv,e_k,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: uv
    REAL(kind=dp), INTENT(out), OPTIONAL     :: e_k
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_e_kin_beads_u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: res
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_uv

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    res=-1.0_dp
    IF (.NOT. failure) THEN
      my_uv => pint_env%uv
      IF (PRESENT(uv)) my_uv => uv
      res=0._dp
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          res=res+pint_env%mass_fict(ib,idim)*my_uv(ib,idim)**2
        END DO
      END DO
      res=res*0.5
      IF (.not.PRESENT(uv)) pint_env%e_kin_beads=res
    END IF
    IF (PRESENT(e_k)) e_k=res
    RETURN
  END SUBROUTINE pint_calc_e_kin_beads_u

  ! ***************************************************************************
  !> \brief  Calculate the virial estimator of the real (quantum) kinetic energy
  !> \author hforbert
  !> \note   This subroutine modifies pint_env%energy(e_kin_virial_id) global
  !>         variable [lwalewski]
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param e_vir ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_calc_e_vir(pint_env,e_vir,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), INTENT(out), OPTIONAL     :: e_vir
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_e_vir', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: res, xcentroid

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    res=-1.0_dp
    IF (.NOT. failure) THEN
      res=0._dp
      DO idim=1,pint_env%ndim
        ! calculate the centroid
        xcentroid = 0._dp
        DO ib=1,pint_env%p
          xcentroid = xcentroid + pint_env%x(ib,idim)
        END DO
        xcentroid = xcentroid/REAL(pint_env%p,dp)
        DO ib=1,pint_env%p
          res=res+(pint_env%x(ib,idim)-xcentroid)*pint_env%f(ib,idim)
        END DO
      END DO
      res=0.5_dp*(REAL(pint_env%ndim,dp)*pint_env%kT-res/REAL(pint_env%p,dp))
      pint_env%energy(e_kin_virial_id) = res
    END IF
    IF (PRESENT(e_vir)) e_vir=res
    RETURN
  END SUBROUTINE pint_calc_e_vir

  ! ***************************************************************************
  !> \brief calculates the energy (potential and kinetic) of the Nose-Hoover
  !>      chain thermostats
  !> \param pint_env the path integral environment
  !> \param error variable to control error logging, stopping,...
  !>        see module cp_error_handling
  !> \author fawzi
  ! ***************************************************************************
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pint_calc_nh_energy(pint_env,error)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_nh_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ekin, epot

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pint_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(pint_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      ekin=0._dp
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          DO inos=1,pint_env%nnos
            ekin=ekin+pint_env%Q(ib)*pint_env%tv(inos,ib,idim)**2
          END DO
        END DO
      END DO
      pint_env%e_kin_t=0.5_dp*ekin
      epot=0._dp
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          DO inos=1,pint_env%nnos
            epot=epot+pint_env%tx(inos,ib,idim)
          END DO
        END DO
      END DO
      pint_env%e_pot_t=pint_env%kT*epot
  END IF
    RETURN
  END SUBROUTINE pint_calc_nh_energy

END MODULE pint_methods
