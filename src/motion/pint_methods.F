!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
!> \param pint_env ...
!> \param input ...
!> \param input_declaration ...
!> \param para_env ...
!> \par    History
!>           Fixed some bugs [hforbert]
!>           Added normal mode transformation [hforbert]
!> \author fawzi
!> \note   Might return an unassociated pointer in parallel on the processors
!>         that are not needed.
! *****************************************************************************
  SUBROUTINE pint_create(pint_env,input,input_declaration,para_env)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: output_file_name, project_name
    INTEGER                                  :: handle, iat, idim, idir, &
                                                ierr, itmp, nrep, prep, stat
    LOGICAL                                  :: explicit, failure, ltmp, &
                                                wrong_input
    REAL(kind=dp)                            :: mass
    REAL(kind=dp), DIMENSION(3, 2)           :: seed
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

    CPASSERT(.NOT.ASSOCIATED(pint_env))
    CPASSERT(ASSOCIATED(input))
    CPASSERT(input%ref_count>0)
    NULLIFY(rep_env)
    pint_section => section_vals_get_subs_vals(input,"MOTION%PINT")
    CALL section_vals_val_get(pint_section,"p",i_val=nrep)
    CALL section_vals_val_get(pint_section,"proc_per_replica",&
         i_val=prep)
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
    CALL section_vals_val_get(input,"GLOBAL%PROJECT_NAME",c_val=project_name)
    CALL section_vals_val_get(input,"GLOBAL%OUTPUT_FILE_NAME",c_val=output_file_name)
    CALL rep_env_create(rep_env, para_env=para_env, input=input,&
         input_declaration=input_declaration,nrep=nrep,prep=prep, row_force=.TRUE.)
    CALL section_vals_val_set(input,"GLOBAL%PROJECT_NAME",c_val=TRIM(project_name))
    IF ( LEN_TRIM(output_file_name) .GT. 0 ) THEN
      CALL section_vals_val_set(input,"GLOBAL%OUTPUT_FILE_NAME",c_val=TRIM(output_file_name))
    ELSE
      CALL section_vals_val_unset(input,"GLOBAL%OUTPUT_FILE_NAME")
    END IF
    IF (.NOT. ASSOCIATED(rep_env)) RETURN

    ALLOCATE(pint_env,STAT=stat)
    CPASSERT(stat==0)
    last_pint_id=last_pint_id+1
    pint_env%id_nr=last_pint_id
    pint_env%ref_count=1
    NULLIFY(pint_env%replicas,pint_env%input,pint_env%staging_env,&
            pint_env%normalmode_env)
    pint_env%p=nrep
    pint_env%replicas => rep_env
    pint_env%ndim=rep_env%ndim
    pint_env%input => input
    CALL section_vals_retain(pint_env%input)

    ! get first step, last step, number of steps, etc
    CALL section_vals_val_get(input,"MOTION%PINT%ITERATION",&
         i_val=itmp)
    pint_env%first_step = itmp
    CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
         explicit=explicit)
    IF ( explicit ) THEN
      CALL section_vals_val_get(input,"MOTION%PINT%MAX_STEP",&
           i_val=itmp)
      pint_env%last_step = itmp
      pint_env%num_steps = pint_env%last_step - pint_env%first_step
    ELSE
      CALL section_vals_val_get(input,"MOTION%PINT%NUM_STEPS",&
           i_val=itmp)
      pint_env%num_steps = itmp
      pint_env%last_step = pint_env%first_step + pint_env%num_steps
    END IF

    CALL section_vals_val_get(pint_section,"DT",&
         r_val=pint_env%dt)
    pint_env%t = pint_env%first_step * pint_env%dt

    CALL section_vals_val_get(pint_section,"nrespa",i_val=pint_env%nrespa)
    CALL section_vals_val_get(pint_section,"Temp",r_val=pint_env%kT)
    CALL section_vals_val_get(pint_section,"T_TOL",&
         r_val=pint_env%t_tol)
    CALL section_vals_val_get(pint_section,"transformation",&
         i_val=pint_env%transform)

    NULLIFY(pint_env%tx,pint_env%tv,pint_env%tv_t,pint_env%tv_old,pint_env%tv_new,pint_env%tf)

    pint_env%nnos = 0
    pint_env%pimd_thermostat = thermostat_none
    nose_section => section_vals_get_subs_vals(input,"MOTION%PINT%NOSE")
    CALL section_vals_get(nose_section, explicit=explicit)
    IF(explicit) THEN
      CALL section_vals_val_get(nose_section,"nnos",i_val=pint_env%nnos)
      IF(pint_env%nnos>0)THEN
        pint_env%pimd_thermostat = thermostat_nose
        ALLOCATE(&
           pint_env%tx(pint_env%nnos,pint_env%p,pint_env%ndim),    &
           pint_env%tv(pint_env%nnos,pint_env%p,pint_env%ndim),    &
           pint_env%tv_t(pint_env%nnos,pint_env%p,pint_env%ndim),  &
           pint_env%tv_old(pint_env%nnos,pint_env%p,pint_env%ndim),&
           pint_env%tv_new(pint_env%nnos,pint_env%p,pint_env%ndim),&
           pint_env%tf(pint_env%nnos,pint_env%p,pint_env%ndim), STAT=stat)
        CPASSERT(stat==0)
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
!       CALL section_vals_val_get(transform_section,"v_tol_nose",r_val=pint_env%v_tol)
!MK ... but we have to initialise v_tol
    pint_env%v_tol = 0.0_dp ! to be fixed

    NULLIFY (pint_env%randomG)

    seed(:,:) = next_rng_seed()
    CALL create_rng_stream(pint_env%randomG,&
                           name="pint_randomG",&
                           distribution_type=GAUSSIAN,&
                           extended_precision=.TRUE.,&
                           seed=seed)

    ALLOCATE(pint_env%e_pot_bead(pint_env%p),STAT=stat)
    CPASSERT(stat==0)
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
    CPASSERT(stat==0)
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
          "MOTION%PINT%STAGING")
       CALL staging_env_create(pint_env%staging_env,transform_section,&
          p=pint_env%p,kT=pint_env%kT)
    ELSE
       transform_section => section_vals_get_subs_vals(input,&
          "MOTION%PINT%NORMALMODE")
       CALL normalmode_env_create(pint_env%normalmode_env,&
          transform_section,p=pint_env%p,kT=pint_env%kT)
       wrong_input=pint_env%nrespa*twopi/(SQRT(pint_env%p/pint_env%normalmode_env%modefactor)*pint_env%kT)/pint_env%dt>10
       CPASSERT(wrong_input)
    END IF
    ALLOCATE(pint_env%mass(pint_env%ndim),STAT=stat)
    CPASSERT(stat==0)
    CALL f_env_add_defaults(f_env_id=pint_env%replicas%f_env_id,&
         f_env=f_env,failure=failure)
    CALL force_env_get(force_env=f_env%force_env,subsys=subsys)
    CALL cp_subsys_get(subsys,particles=particles)

!TODO length of pint_env%mass is redundant
    idim=0
    DO iat=1,pint_env%ndim/3
       CALL get_atomic_kind(particles%els(iat)%atomic_kind,mass=mass)
       DO idir=1,3
          idim=idim+1
          pint_env%mass(idim)=mass
       END DO
    END DO
    CALL f_env_rm_defaults(f_env,ierr)
    CPASSERT(ierr==0)

    ALLOCATE(pint_env%Q(pint_env%p),&
         pint_env%mass_beads(pint_env%p,pint_env%ndim),&
         pint_env%mass_fict(pint_env%p,pint_env%ndim),STAT=stat)
    CPASSERT(stat==0)
    IF (pint_env%transform == transformation_stage) THEN
       CALL staging_init_masses(pint_env%staging_env,mass=pint_env%mass,&
            mass_beads=pint_env%mass_beads,mass_fict=pint_env%mass_fict,&
            Q=pint_env%Q)
    ELSE
       CALL normalmode_init_masses(pint_env%normalmode_env, &
            mass=pint_env%mass, mass_beads=pint_env%mass_beads, &
            mass_fict=pint_env%mass_fict, Q=pint_env%Q)
    END IF

    NULLIFY(pint_env%gle)
    gle_section => section_vals_get_subs_vals(input,"MOTION%PINT%GLE")
    CALL section_vals_get(gle_section, explicit=explicit)
    IF(explicit) THEN
      ALLOCATE(pint_env%gle,stat=stat)
      CPASSERT(stat==0)
      CALL gle_init(pint_env%gle, dt=pint_env%dt/pint_env%nrespa, temp=pint_env%kT,&
           section=gle_section)
      IF (pint_env%pimd_thermostat==thermostat_none .AND. pint_env%gle%ndim .GT. 0) THEN
         pint_env%pimd_thermostat=thermostat_gle

         ! initialize a GLE with ALL degrees of freedom on node 0,
         ! as it seems to me that here everything but force eval is replicated
         pint_env%gle%loc_num_gle=pint_env%p*pint_env%ndim
         pint_env%gle%glob_num_gle=pint_env%gle%loc_num_gle
         ALLOCATE(pint_env%gle%map_info%index(pint_env%gle%loc_num_gle))
         CPASSERT(stat==0)
         DO itmp = 1, pint_env%gle%loc_num_gle
             pint_env%gle%map_info%index(itmp)=itmp
         ENDDO
         CALL gle_thermo_create(pint_env%gle,pint_env%gle%loc_num_gle)

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
         CALL pint_gle_init(pint_env)
      END IF
    END IF

    CALL section_vals_val_get(pint_section,"FIX_CENTROID_POS",&
         l_val=ltmp)
    IF ( ltmp .AND. (pint_env%transform .EQ. transformation_normal) ) THEN
      pint_env%first_propagated_mode = 2
    ELSE
      pint_env%first_propagated_mode = 1
    END IF


    CALL timestop(handle)

    RETURN
  END SUBROUTINE pint_create

! ***************************************************************************
!> \brief Retain a path integral environment
!> \param pint_env the pint_env to retain
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pint_retain(pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    pint_env%ref_count=pint_env%ref_count+1
    RETURN
  END SUBROUTINE pint_retain

! ***************************************************************************
!> \brief Release a path integral environment
!> \param pint_env the pint_env to release
!> \par History
!>      Added normal mode transformation [hforbert]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pint_release(pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(pint_env)) THEN
      CPASSERT(pint_env%ref_count>0)
      pint_env%ref_count=pint_env%ref_count-1
      IF (pint_env%ref_count==0) THEN
        CALL rep_env_release(pint_env%replicas)
        CALL section_vals_release(pint_env%input)
        IF (ASSOCIATED(pint_env%staging_env)) THEN
           CALL staging_release(pint_env%staging_env)
        END IF
        IF (ASSOCIATED(pint_env%normalmode_env)) THEN
           CALL normalmode_release(pint_env%normalmode_env)
        END IF
        CALL delete_rng_stream(pint_env%randomG)

        DEALLOCATE(pint_env%mass,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%e_pot_bead,STAT=stat)
        CPASSERT(stat==0)

        DEALLOCATE(pint_env%x,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%v,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%f,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%external_f,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%mass_beads,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%mass_fict,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%ux,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%uv,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%uv_t,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%uv_new,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%uf,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%uf_h,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%rtmp_ndim,STAT=stat)
        CPASSERT(stat==0)
        DEALLOCATE(pint_env%rtmp_natom,STAT=stat)
        CPASSERT(stat==0)

        IF(pint_env%pimd_thermostat==thermostat_nose) THEN
          DEALLOCATE(pint_env%tx,STAT=stat)
          CPASSERT(stat==0)
          DEALLOCATE(pint_env%tv,STAT=stat)
          CPASSERT(stat==0)
          DEALLOCATE(pint_env%tv_t,STAT=stat)
          CPASSERT(stat==0)
          DEALLOCATE(pint_env%tv_old,STAT=stat)
          CPASSERT(stat==0)
          DEALLOCATE(pint_env%tv_new,STAT=stat)
          CPASSERT(stat==0)
          DEALLOCATE(pint_env%tf,STAT=stat)
          CPASSERT(stat==0)
        ELSEIF(pint_env%pimd_thermostat==thermostat_gle) THEN
          CALL gle_dealloc(pint_env%gle)
        END IF

        DEALLOCATE(pint_env%Q,STAT=stat)
        CPASSERT(stat==0)

        DEALLOCATE(pint_env,STAT=stat)
        CPASSERT(stat==0)

      END IF
    END IF

    NULLIFY(pint_env)

    RETURN
  END SUBROUTINE pint_release

! ***************************************************************************
!> \brief Tests the path integral methods
!> \param para_env parallel environment
!> \param input the input to test
!> \param input_declaration ...
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_test(para_env,input,input_declaration)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration

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

  CPASSERT(ASSOCIATED(para_env))
  CPASSERT(ASSOCIATED(input))
  CPASSERT(para_env%ref_count>0)
  CPASSERT(input%ref_count>0)
  NULLIFY(pint_env)
  logger => cp_get_default_logger()
   CALL pint_create(pint_env,input,input_declaration,para_env)
  IF (ASSOCIATED(pint_env)) THEN
     ALLOCATE(x1(pint_env%ndim,pint_env%p),stat=stat)
     CPASSERT(stat==0)
     x1(:,:)=pint_env%x
     CALL pint_x2u(pint_env)
     pint_env%x=0._dp
     CALL pint_u2x(pint_env)
     err=0._dp
     DO i=1,pint_env%ndim
        err=MAX(err,ABS(x1(1,i)-pint_env%x(1,i)))
     END DO
     CALL cp_log(logger,cp_note_level+1,routineP,"diff_r1="//cp_to_string(err),&
          local=.FALSE.)

     CALL pint_calc_uf_h(pint_env,e_h=e_h)
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
     CALL pint_f2uf(pint_env)
     err=0._dp
     DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
           err=MAX(err,ABS(pint_env%uf(ib,idim)-pint_env%uf_h(ib,idim)))
        END DO
     END DO
     CALL cp_log(logger,cp_note_level+1,routineP,"diff_f_h="//cp_to_string(err),&
          local=.FALSE.)
  END IF
    RETURN
  END SUBROUTINE pint_test

! ***************************************************************************
!> \brief  Perform a path integral simulation
!> \param  para_env parallel environment
!> \param  input the input to test
!> \param input_declaration ...
!> \param globenv ...
!> \par    History
!>         2003-11 created [fawzi]
!>         2009-12-14 globenv parameter added to handle soft exit
!>           requests [lwalewski]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE do_pint_run( para_env, input, input_declaration, globenv)
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: input
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(global_environment_type), POINTER   :: globenv

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

    CPASSERT(ASSOCIATED(para_env))
    CPASSERT(ASSOCIATED(input))
    CPASSERT(para_env%ref_count>0)
    CPASSERT(input%ref_count>0)

    ! check if helium solvent is present
    NULLIFY(helium_section)
    helium_section => section_vals_get_subs_vals(input,&
                      "MOTION%PINT%HELIUM")
    CALL section_vals_get(helium_section,explicit=explicit)
    IF ( explicit ) THEN
      CALL section_vals_val_get(helium_section,"_SECTION_PARAMETERS_",&
           l_val=solvent_present)
    ELSE
      solvent_present = .FALSE.
    END IF

    ! check if there is anything but helium
    IF (solvent_present) THEN
      CALL section_vals_val_get(helium_section,"HELIUM_ONLY",&
        l_val=helium_only)
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
      CALL helium_create(helium,input)
      CALL helium_init(helium,pint_env)
      CALL helium_do_run(helium,globenv)
      CALL helium_release(helium)

    CASE (solute_only_mid)
      CALL pint_create(pint_env,input,input_declaration,para_env)
      CALL pint_init(pint_env)
      CALL pint_do_run(pint_env,globenv)
      CALL pint_release(pint_env)

    CASE (solute_with_helium_mid)
      CALL pint_create(pint_env,input,input_declaration,para_env)
      ! init pint wihtout helium forces (they are not yet initialized)
      CALL pint_init(pint_env)
      ! init helium with solute's positions (they are already initialized)
      CALL helium_create(helium,input,solute=pint_env)
      CALL helium_init(helium,pint_env)
      ! reinit pint forces with helium forces (they are now initialized)
      CALL pint_init_f( pint_env, helium_solvent=helium)

      CALL pint_do_run(pint_env,globenv,helium=helium)
      CALL helium_release(helium)
      CALL pint_release(pint_env)

    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
        routineP,"Unknown mode ("//TRIM(ADJUSTL(cp_to_string(mode)))//")",&
        failure)
    END SELECT


    CALL timestop(handle)

    RETURN
  END SUBROUTINE do_pint_run

! ***************************************************************************
!> \brief  Reads the restart, initializes the beads, etc.
!> \param pint_env ...
!> \par    History
!>           11.2003 created [fawzi]
!>           actually ASSIGN input pointer [hforbert]
!>           2010-12-16 turned into a wrapper routine [lwalewski]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pint_init(pint_env)

    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    IF ( failure ) THEN
      RETURN
    END IF

    CALL pint_init_x( pint_env)
    CALL pint_init_v( pint_env)
    CALL pint_init_t( pint_env)
    CALL pint_init_f( pint_env)

    RETURN
  END SUBROUTINE pint_init


! ***************************************************************************
!> \brief  Assign initial postions to the beads.
!> \param pint_env ...
!> \date   2010-12-15
!> \author Lukasz Walewski
!> \note  Initialization is done in the following way:
!>           1. assign all beads with the same classical positions from
!>              FORCE_EVAL (hot start)
!>           2. spread the beads around classical positions as if they were
!>              free particles (if requested)
!>           3. replace positions generated in steps 1-2 with the explicit
!>              ones if they are explicitly given in the input structure
!>           4. apply Gaussian noise to the positions generated so far (if
!>              requested)
! *****************************************************************************
  SUBROUTINE pint_init_x( pint_env)

    TYPE(pint_env_type), POINTER             :: pint_env

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
    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
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
         l_val=ltmp)
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%LEVY_TEMP_FACTOR",&
         r_val=tcorr)
    IF ( ltmp ) THEN

      NULLIFY(bx)
      ALLOCATE(bx(3*pint_env%p), STAT=istat)
      CPASSERT(istat==0)
      NULLIFY (rng_gaussian)
      CALL section_vals_val_get(pint_env%input,&
           "MOTION%PINT%INIT%LEVY_SEED",i_val=input_seed)
      seed(:,:) = REAL(input_seed,KIND=dp)
!      seed(:,:) = next_rng_seed()
      CALL create_rng_stream(rng_gaussian,&
                             name="tmp_rng_gaussian",&
                             distribution_type=GAUSSIAN,&
                             extended_precision=.TRUE.,&
                             seed=seed)

      CALL section_vals_val_get(pint_env%input,&
           "MOTION%PINT%INIT%LEVY_CORRELATED",&
           l_val=levycorr)

      IF ( levycorr ) THEN

        ! correlated Levy walk - the same path for all atoms
        x0 = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        CALL pint_levy_walk( x0, pint_env%p, 1.0_dp, bx, rng_gaussian)
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
          CALL pint_levy_walk( x0, pint_env%p, var, bx, rng_gaussian)
          DO ic = 1, 3
            idim = idim + 1
            DO ib=1,pint_env%p
              pint_env%x(ib,idim) = pint_env%x(ib,idim) + bx(3*(ib-1)+ic)
            END DO
          END DO
        END DO

      END IF

      CALL delete_rng_stream(rng_gaussian)
      DEALLOCATE(bx, STAT=istat)
      CPASSERT(istat==0)
      done_levy = .TRUE.
    END IF

    done_init = .FALSE.
    NULLIFY(input_section)
    input_section => section_vals_get_subs_vals(pint_env%input,&
                     "MOTION%PINT%BEADS%COORD")
    CALL section_vals_get(input_section,explicit=explicit)
    IF (explicit) THEN
      CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
           n_rep_val=n_rep_val)
      IF (n_rep_val>0) THEN
        CPASSERT(n_rep_val==1)
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
             r_vals=r_vals)
        CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim,&
             cp_failure_level,cp_assertion_failed,&
             "Invalid size of MOTION%PINT%BEADS%COORD "//&
             CPSourceFileRef,&
             routineP,failure)
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
         l_val=ltmp)
    IF ( ltmp ) THEN
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          pint_env%x(ib,idim) = pint_env%x(ib,idim) + &
                                next_random_number(rng_stream=pint_env%randomG,&
                                variance=pint_env%beta/&
                                SQRT(12.0_dp*pint_env%mass(idim)))
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
    CALL pint_write_line(msg)

    IF ( done_levy ) THEN
      WRITE(msg,'(A,F6.3)') "Levy walk at effective temperature: ", tcorr
    END IF

    IF ( done_rand ) THEN
      WRITE(msg,'(A)') "Added gaussian noise to the positions of the beads."
      CALL pint_write_line(msg)
    END IF

    RETURN
  END SUBROUTINE pint_init_x


! ***************************************************************************
!> \brief  Initialize velocities
!> \param  pint_env the pint env in which you should initialize the
!>         velocity
!> \par    History
!>         2010-12-16 gathered all velocity-init code here [lwalewski]
!>         2011-04-05 added centroid velocity initialization [lwalewski]
!>         2011-12-19 removed optional parameter kT, target temperature is
!>                    now determined from the input directly [lwalewski]
!> \author fawzi
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
! *****************************************************************************
  SUBROUTINE pint_init_v( pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

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
    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    IF ( failure ) THEN
      RETURN
    END IF

    NULLIFY(logger)
    logger => cp_get_default_logger()

    ! read the velocities from the input file if they are given explicitly
    vels_present = .FALSE.
    NULLIFY(input_section)
    input_section => section_vals_get_subs_vals(pint_env%input,&
                     "FORCE_EVAL%SUBSYS%VELOCITY")
    CALL section_vals_get(input_section,explicit=explicit)
    IF ( explicit ) THEN

      CALL section_vals_val_get(input_section,"PINT_UNIT",&
                                c_val=unit_str)
      unit_conv = cp_unit_to_cp2k(1.0_dp,TRIM(unit_str))

      ! assign all the beads with the same velocities from FORCE_EVAL%SUBSYS%VELOCITY
      NULLIFY(r_vals)
      CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
           n_rep_val=n_rep_val)
      stmp = ""
      WRITE(stmp,*) n_rep_val
      msg = "Invalid number of atoms in FORCE_EVAL%SUBSYS%VELOCITY ("//&
            TRIM(ADJUSTL(stmp))//")."
      CALL cp_assert(3*n_rep_val==pint_env%ndim,&
                     cp_failure_level, cp_assertion_failed,&
                     routineP//" from "//&
CPSourceFileRef,&
                     msg,failure)
      DO ia = 1, pint_env%ndim/3
        CALL section_vals_val_get(input_section, "_DEFAULT_KEYWORD_",&
             i_rep_val=ia, r_vals=r_vals)
        itmp = SIZE(r_vals)
        stmp = ""
        WRITE(stmp,*) itmp
        msg = "Number of coordinates != 3 in FORCE_EVAL%SUBSYS%VELOCITY ("//&
              TRIM(ADJUSTL(stmp))//")."
        CALL cp_assert(itmp==3,&
                       cp_failure_level, cp_assertion_failed,&
                       routineP//" from "//&
CPSourceFileRef,&
                       msg,failure)
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
         l_val=done_scale)
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
      CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v)
      first_mode = 2
    ELSE
      ! ...for all the modes
      first_mode = 1
    END IF
    DO idim = 1, SIZE(pint_env%uv,2)
      DO ib = first_mode, SIZE(pint_env%uv,1)
        pint_env%uv(ib,idim) = &
          next_random_number(rng_stream=pint_env%randomG,&
                             variance=target_t/pint_env%mass_fict(ib,idim))
      END DO
    END DO

    ! add random component to the centroid velocity if requsted
    done_sped = .FALSE.
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%CENTROID_SPEED",&
         l_val=ltmp)
    IF (ltmp) THEN
      CALL pint_u2x(pint_env,ux=pint_env%uv,x=pint_env%v)
      DO idim=1,pint_env%ndim
        rtmp = next_random_number(rng_stream=pint_env%randomG,&
                                variance=pint_env%mass(idim)*pint_env%kT)&
                                /pint_env%mass(idim)
        DO ib=1,pint_env%p
          pint_env%v(ib,idim)=pint_env%v(ib,idim)+rtmp
        END DO
      END DO
      CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v)
      done_sped = .TRUE.
    END IF

    ! quench (set to zero) velocities for all the modes if requested
    ! (disregard all the initialization done so far)
    done_quench = .FALSE.
    CALL section_vals_val_get(pint_env%input,&
         "MOTION%PINT%INIT%VELOCITY_QUENCH",&
         l_val=ltmp)
    IF (ltmp) THEN
      DO idim=1,pint_env%ndim
        DO ib=1,pint_env%p
          pint_env%v(ib,idim) = 0.0_dp
        END DO
      END DO
      CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v)
      done_quench = .TRUE.
    END IF

    ! set the velocities to the values from the input if they are explicit
    ! (disregard all the initialization done so far)
    done_init = .FALSE.
    NULLIFY(input_section)
    input_section => section_vals_get_subs_vals(pint_env%input,&
                     "MOTION%PINT%BEADS%VELOCITY")
    CALL section_vals_get(input_section,explicit=explicit)
    IF (explicit) THEN
      CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
           n_rep_val=n_rep_val)
      IF (n_rep_val>0) THEN
        CPASSERT(n_rep_val==1)
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
             r_vals=r_vals)
        CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim,&
                       cp_failure_level,cp_assertion_failed,&
                       "Invalid size of MOTION%PINT%BEAD%VELOCITY "//&
                       CPSourceFileRef,&
                       routineP,failure)
        itmp=0
        DO idim=1,pint_env%ndim
          DO ib=1,pint_env%p
            itmp=itmp+1
            pint_env%v(ib,idim)=r_vals(itmp)
          END DO
        END DO
        CALL pint_x2u(pint_env,ux=pint_env%uv,x=pint_env%v)
        done_init = .TRUE.
      END IF
    END IF

    unit_conv = cp_unit_from_cp2k(1.0_dp,"K")
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
    CALL pint_write_line(msg)

    IF ( done_init .AND. done_quench ) THEN
      msg = "WARNING: exclusive options requested (velocity restart and quenching)"
      CALL pint_write_line(msg)
      msg = "WARNING: velocity restart took precedence"
      CALL pint_write_line(msg)
    END IF

    IF ( (.NOT. done_init) .AND. (.NOT. done_quench) ) THEN
      IF (vels_present .AND. done_scale) THEN
        WRITE(stmp1,'(F10.2)') actual_t * unit_conv
        WRITE(stmp2,'(F10.2)') target_t * unit_conv
        msg = "Scaled initial velocities from "//TRIM(ADJUSTL(stmp1))//&
              " to "//TRIM(ADJUSTL(stmp2))//" K as requested."
        CALL pint_write_line(msg)
      END IF
      IF ( done_sped ) THEN
        msg = "Added random component to the initial centroid velocities."
        CALL pint_write_line(msg)
      END IF
    END IF

    RETURN
  END SUBROUTINE pint_init_v


! ***************************************************************************
!> \brief  Assign initial postions and velocities to the thermostats.
!> \param pint_env ...
!> \param kT ...
!> \date   2010-12-15
!> \author Lukasz Walewski
!> \note   Extracted from pint_init
! *****************************************************************************
  SUBROUTINE pint_init_t( pint_env, kT)

    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), INTENT(in), OPTIONAL      :: kT

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_t', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, ii, inos, n_rep_val
    LOGICAL                                  :: explicit, failure, gle_restart
    REAL(kind=dp)                            :: mykt
    REAL(kind=dp), DIMENSION(:), POINTER     :: r_vals
    TYPE(section_vals_type), POINTER         :: input_section

    failure=.FALSE.
    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
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
                                  variance=mykt/pint_env%Q(ib))
          END DO
        END DO
      END DO

      NULLIFY(input_section)
      input_section => section_vals_get_subs_vals(pint_env%input,&
                       "MOTION%PINT%NOSE%COORD")
      CALL section_vals_get(input_section,explicit=explicit)
      IF (explicit) THEN
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                  n_rep_val=n_rep_val)
        IF (n_rep_val>0) THEN
          CPASSERT(n_rep_val==1)
          CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                    r_vals=r_vals)
          CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim*pint_env%nnos,&
                         cp_failure_level,cp_assertion_failed,&
                         "Invalid size of MOTION%PINT%NOSE%COORD "//&
                         CPSourceFileRef,&
                         routineP,failure)
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
                       "MOTION%PINT%NOSE%VELOCITY")
      CALL section_vals_get(input_section,explicit=explicit)
      IF (explicit) THEN
        CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                  n_rep_val=n_rep_val)
        IF (n_rep_val>0) THEN
          CPASSERT(n_rep_val==1)
          CALL section_vals_val_get(input_section,"_DEFAULT_KEYWORD_",&
                                    r_vals=r_vals)
          CALL cp_assert(SIZE(r_vals)==pint_env%p*pint_env%ndim*pint_env%nnos,&
                         cp_failure_level,cp_assertion_failed,&
                         "Invalid size of MOTION%PINT%NOSE%VELOCITY "//&
                         CPSourceFileRef,&
                         routineP,failure)
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
                       "MOTION%PINT%GLE")
      CALL section_vals_get(input_section,explicit=explicit)
      IF (explicit) THEN
        CALL restart_gle(pint_env%gle,input_section,save_mem=.FALSE.,&
                         restart=gle_restart)
      END IF
    END IF

    RETURN
  END SUBROUTINE pint_init_t


! ***************************************************************************
!> \brief  Prepares the forces, etc. to perform an MD step
!> \param pint_env ...
!> \param helium_solvent ...
!> \par    History
!>           Added nh_energy calculation [hforbert]
!>           Bug fixes for no thermostats [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_init_f( pint_env, helium_solvent)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_solvent

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_init_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: e_h

    failure=.FALSE.
    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    CALL pint_x2u(pint_env)
    CALL pint_calc_uf_h(pint_env=pint_env,e_h=e_h)
    CALL pint_calc_f(pint_env)

    ! add helium forces to the solute's internal ones
    ! Assume that helium has been already initialized and helium_solvent
    ! contains proper forces in force_avrg array.
    IF (PRESENT(helium_solvent)) THEN
      IF (ASSOCIATED(helium_solvent)) THEN
        pint_env%f(:,:) = pint_env%f(:,:) + helium_solvent%force_avrg(:,:)
      END IF
    END IF
    CALL pint_f2uf(pint_env)

    ! set the centroid forces to 0 if FIX_CENTROID_POS
    IF ( pint_env%first_propagated_mode .EQ. 2 ) THEN
      pint_env%uf(1,:) = 0.0_dp
    END IF

    CALL pint_calc_e_kin_beads_u(pint_env)
    CALL pint_calc_e_vir(pint_env)
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
      CALL pint_calc_nh_energy(pint_env)
    END IF
    RETURN
  END SUBROUTINE pint_init_f


! ***************************************************************************
!> \brief  Perform the PIMD simulation (main MD loop)
!> \param pint_env ...
!> \param globenv ...
!> \param helium ...
!> \par    History
!>         2003-11 created [fawzi]
!>         renamed from pint_run to pint_do_run because of conflicting name
!>           of pint_run in input_constants [hforbert]
!>         2009-12-14 globenv parameter added to handle soft exit
!>           requests [lwalewski]
!> \author Fawzi Mohamed
!> \note   Everything should be read for an md step.
! *****************************************************************************
  SUBROUTINE pint_do_run(pint_env, globenv, helium)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_do_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: step
    LOGICAL                                  :: failure, should_stop
    REAL(kind=dp)                            :: scal

    failure=.FALSE.
    CPASSERT(ASSOCIATED(pint_env))

    NULLIFY(pint_env%logger)
    pint_env%logger => cp_get_default_logger()

    ! create iteration level and initialize iteration info
    CALL cp_add_iter_level(pint_env%logger%iter_info,"MD")
    CALL cp_iterate(pint_env%logger%iter_info,iter_nr=pint_env%first_step)
    pint_env%iter = pint_env%first_step

    ! write the properties at 0-th step
    CALL pint_calc_energy(pint_env)
    CALL pint_write_ener(pint_env)
    CALL pint_write_centroids(pint_env)
    CALL pint_write_trajectory(pint_env)
    CALL pint_write_com(pint_env)
    CALL pint_write_rgyr(pint_env)

    ! main PIMD loop
    DO step = 1, pint_env%num_steps

      pint_env%iter = pint_env%iter + 1
      CALL cp_iterate(pint_env%logger%iter_info,&
                      last=(step==pint_env%num_steps),&
                      iter_nr=pint_env%iter)
      pint_env%t = pint_env%t + pint_env%dt

      IF (pint_env%t_tol > 0.0_dp) THEN
        IF (ABS(2._dp*pint_env%e_kin_beads/(pint_env%p*pint_env%ndim)&
          -pint_env%kT)>pint_env%t_tol) THEN
          scal=SQRT(pint_env%kT*(pint_env%p*pint_env%ndim)/(2.0_dp*pint_env%e_kin_beads))
          pint_env%uv=scal*pint_env%uv
          CALL pint_init_f(pint_env,helium_solvent=helium)
        END IF
      END IF
      CALL pint_step(pint_env,helium_solvent=helium)

      CALL pint_write_ener(pint_env)
      CALL pint_write_centroids(pint_env)
      CALL pint_write_trajectory(pint_env)
      CALL pint_write_com(pint_env)
      CALL pint_write_rgyr(pint_env)

      CALL write_restart(root_section=pint_env%input,&
           pint_env=pint_env, helium_env=helium)

      ! exit from the main loop if soft exit has been requested
      CALL external_control(should_stop,"MD",globenv=globenv)
      IF (should_stop) EXIT

    END DO

    ! remove iteration level
    CALL cp_rm_iter_level(pint_env%logger%iter_info,"MD")

    RETURN
  END SUBROUTINE pint_do_run


! ***************************************************************************
!> \brief  Does an MD step (and nrespa harmonic evaluations)
!> \param pint_env ...
!> \param helium_solvent ...
!> \par    History
!>           various bug fixes [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_step(pint_env,helium_solvent)
    TYPE(pint_env_type), POINTER             :: pint_env
    TYPE(helium_solvent_type), OPTIONAL, &
      POINTER                                :: helium_solvent

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
    CPASSERT(ASSOCIATED(pint_env))

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
     CALL pint_calc_uf_h(pint_env=pint_env,e_h=e_h)
     IF (iresp==pint_env%nrespa) THEN
        CALL pint_u2x(pint_env)
        CALL pint_calc_f(pint_env)

        ! perform helium step and add helium forces
        IF (PRESENT(helium_solvent)) THEN
          IF (ASSOCIATED(helium_solvent)) THEN
            helium_solvent%current_step = pint_env%iter
            helium_solvent%origin = pint_com_pos(pint_env)
            CALL helium_step( helium_solvent, pint_env)
            pint_env%f(:,:)=pint_env%f(:,:)+helium_solvent%force_avrg(:,:)
          END IF
        END IF

        CALL pint_f2uf(pint_env)
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
        CALL pint_gle_step(pint_env)
        pint_env%uv=pint_env%uv_t
     ELSE
        pint_env%uv=pint_env%uv_t
     ENDIF
  END DO

  ! calculate the energy components
  CALL pint_calc_energy( pint_env)

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
!          failure)
!     CALL f_env_rm_defaults(f_env,new_error,ierr)


    time_stop = m_walltime()
    pint_env%time_per_step = time_stop - time_start
    CALL pint_write_step_info(pint_env)
    CALL timestop(handle)

    RETURN
  END SUBROUTINE pint_step


! ***************************************************************************
!> \brief  Calculate the energy components (private wrapper function)
!> \param pint_env ...
!> \date   2011-01-07
!> \author Lukasz Walewski
! *****************************************************************************
  SUBROUTINE pint_calc_energy( pint_env)

    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_energy', &
      routineP = moduleN//':'//routineN

    CALL pint_calc_e_kin_beads_u(pint_env)
    CALL pint_calc_e_vir(pint_env)

    IF (pint_env%pimd_thermostat==thermostat_nose) THEN
      CALL pint_calc_nh_energy(pint_env)
    ELSEIF (pint_env%pimd_thermostat==thermostat_gle) THEN
      CALL pint_calc_gle_energy(pint_env)
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
!> \param e_h ...
!> \par    History
!>           Added normal mode transformation [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_calc_uf_h(pint_env,e_h)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(KIND=dp), INTENT(OUT)               :: e_h

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_uf_h', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (pint_env%transform == transformation_stage) THEN
      CALL staging_calc_uf_h(pint_env%staging_env,&
                             pint_env%mass_beads,&
                             pint_env%ux,&
                             pint_env%uf_h,&
                             pint_env%e_pot_h)
    ELSE
      CALL normalmode_calc_uf_h(pint_env%normalmode_env,&
                                pint_env%mass_beads,&
                                pint_env%ux,&
                                pint_env%uf_h,&
                                pint_env%e_pot_h)
    END IF
    e_h=pint_env%e_pot_h
    pint_env%uf_h=pint_env%uf_h/pint_env%mass_fict
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
!> \par    History
!>           2009-06-15 moved helium calls out from here [lwalewski]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_calc_f(pint_env,x,f,e)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: x
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(out), OPTIONAL, TARGET          :: f
    REAL(kind=dp), DIMENSION(:), &
      INTENT(out), OPTIONAL, TARGET          :: e

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_f', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:), POINTER     :: my_e
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_f, my_x

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
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
    CALL rep_env_calc_e_f(pint_env%replicas,calc_f=.TRUE.)
    DO idim=1,pint_env%ndim
      DO ib=1,pint_env%p
        !ljw: is that fine ? - idim <-> ib
        my_f(ib,idim) = pint_env%replicas%f(idim,ib)
      END DO
    END DO
    my_e=pint_env%replicas%f(SIZE(pint_env%replicas%f,1),:)

    RETURN
  END SUBROUTINE pint_calc_f

! ***************************************************************************
!> \brief  Calculate the kinetic energy of the beads (in the u variables)
!> \param pint_env ...
!> \param uv ...
!> \param e_k ...
!> \par    History
!>         Bug fix to give my_uv a default location if not given in call [hforbert]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_calc_e_kin_beads_u(pint_env,uv,e_k)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), DIMENSION(:, :), &
      INTENT(in), OPTIONAL, TARGET           :: uv
    REAL(kind=dp), INTENT(out), OPTIONAL     :: e_k

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_e_kin_beads_u', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: res
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: my_uv

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    res=-1.0_dp
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
    IF (PRESENT(e_k)) e_k=res
    RETURN
  END SUBROUTINE pint_calc_e_kin_beads_u

! ***************************************************************************
!> \brief  Calculate the virial estimator of the real (quantum) kinetic energy
!> \param pint_env ...
!> \param e_vir ...
!> \author hforbert
!> \note   This subroutine modifies pint_env%energy(e_kin_virial_id) global
!>         variable [lwalewski]
! *****************************************************************************
  SUBROUTINE pint_calc_e_vir(pint_env,e_vir)
    TYPE(pint_env_type), POINTER             :: pint_env
    REAL(kind=dp), INTENT(out), OPTIONAL     :: e_vir

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_e_vir', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: res, xcentroid

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
    res=-1.0_dp
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
    IF (PRESENT(e_vir)) e_vir=res
    RETURN
  END SUBROUTINE pint_calc_e_vir

! ***************************************************************************
!> \brief calculates the energy (potential and kinetic) of the Nose-Hoover
!>      chain thermostats
!> \param pint_env the path integral environment
!> \author fawzi
! *****************************************************************************
  SUBROUTINE pint_calc_nh_energy(pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_calc_nh_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ib, idim, inos
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: ekin, epot

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pint_env))
    CPASSERT(pint_env%ref_count>0)
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
    RETURN
  END SUBROUTINE pint_calc_nh_energy

END MODULE pint_methods
