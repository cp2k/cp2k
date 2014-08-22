!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief qs_environement methods that use many other modules
!> \par History
!>      09.2002 created [fawzi]
!>      - local atom distribution (25.06.2003,MK)
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_environment_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_deallocate_matrix_set,&
                                             cp_dbcsr_p_type,&
                                             dbcsr_distribution_obj
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
  USE cp_ddapc_types,                  ONLY: cp_ddapc_release
  USE cp_ddapc_util,                   ONLY: cp_ddapc_init
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE f77_blas
  USE input_constants,                 ONLY: do_ppl_analytic,&
                                             do_ppl_grid,&
                                             kg_tnadd_embed,&
                                             use_aux_fit_basis_set
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_methods,                  ONLY: pw_env_create,&
                                             pw_env_rebuild
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_release,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_transfer
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_ppl_grid,&
                                             calculate_rho_core,&
                                             calculate_rho_nlcc
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools
  USE qs_outer_scf,                    ONLY: outer_loop_variables_count
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_scf_types,                    ONLY: scf_env_did_change
  USE scf_control_types,               ONLY: scf_control_type
  USE scp_environment,                 ONLY: scp_nddo_init
  USE scp_environment_types,           ONLY: scp_env_create,&
                                             scp_env_release,&
                                             scp_environment_type
  USE task_list_methods,               ONLY: generate_qs_task_list
  USE task_list_types,                 ONLY: allocate_task_list,&
                                             deallocate_task_list,&
                                             task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_environment_methods'

  PUBLIC :: qs_env_update_s_mstruct,&
            qs_env_rebuild_pw_env,&
            qs_env_setup
!***
CONTAINS

! *****************************************************************************
!> \brief initializes various components of the qs_env, that need only
!>      atomic_kind_set, cell, dft_control, scf_control, c(i)%nmo,
!>      c(i)%nao, and particle_set to be initialized.
!>      The previous components of qs_env must be valid.
!>      Initializes pools, charges and pw_env.
!> \param qs_env the qs_env to set up
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_env_setup(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_ao_aux_fit, &
                                                nhistory, nvariables, stat
    LOGICAL                                  :: failure, scp_nddo
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: outer_scf_history
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scp_environment_type), POINTER      :: scp_env

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(atomic_kind_set, dft_control, scf_control, qs_charges, para_env,&
         distribution_2d, molecule_kind_set, molecule_set, particle_set, scp_env)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set,&
                    scf_control=scf_control,&
                    para_env=para_env,&
                    error=error)
    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
    scp_nddo = dft_control%qs_control%se_control%scp
    IF (.NOT. failure) THEN
       ! allocate qs_charges
       CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
       CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
       CALL qs_charges_release(qs_charges,error=error)

       ! outer scf setup
       IF (scf_control%outer_scf%have_scf) THEN
          nvariables=outer_loop_variables_count(scf_control,error=error)
          nhistory=scf_control%outer_scf%extrapolation_order
          ALLOCATE(outer_scf_history(nvariables,nhistory),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL set_qs_env(qs_env,outer_scf_history=outer_scf_history,error=error)
          CALL set_qs_env(qs_env,outer_scf_ihistory=0,error=error)
       ENDIF

       ! set up pw_env
       CALL qs_env_rebuild_pw_env(qs_env, error=error)

       ! rebuilds fm_pools

       ! XXXX should get rid of the mpools
       IF (ASSOCIATED(qs_env%mos)) THEN
           CALL mpools_rebuild_fm_pools(qs_env%mpools,mos=qs_env%mos,&
                  blacs_env=qs_env%blacs_env, para_env=para_env,&
                  error=error)
       ENDIF

       ! If we use auxiliary density matrix methods rebuild fm_pools
       IF( qs_env%dft_control%do_admm ) THEN
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                  nsgf=n_ao_aux_fit,&
                                  basis_set_id=use_aux_fit_basis_set)
         CALL mpools_rebuild_fm_pools(qs_env%mpools_aux_fit,mos=qs_env%mos_aux_fit,&
              blacs_env=qs_env%blacs_env, para_env=para_env,&
              error=error)
        END IF

       ! create 2d distribution

       CALL distribute_molecules_2d(cell=qs_env%cell,&
                                    particle_kind_set=atomic_kind_set,&
                                    particle_set=particle_set,&
                                    molecule_kind_set=molecule_kind_set,&
                                    molecule_set=molecule_set,&
                                    distribution_2d=distribution_2d,&
                                    blacs_env=qs_env%blacs_env,&
                                    force_env_section=qs_env%input, error=error)

       ! and use it to create the dbcsr_dist, which should be the sole user of distribution_2d by now.
       ALLOCATE(dbcsr_dist)
       CALL cp_dbcsr_dist2d_to_dist (distribution_2d, dbcsr_dist, error)
       CALL set_qs_env(qs_env,dbcsr_dist=dbcsr_dist,error=error)

       ! also keep distribution_2d in qs_env
       CALL set_qs_env(qs_env,distribution_2d=distribution_2d,error=error)
       CALL distribution_2d_release(distribution_2d, error=error)

       ! SCP
       IF ( scp_nddo ) THEN
       ! create the SCP env
          CALL scp_env_create ( scp_env, error )
          CALL scp_nddo_init ( scp_env, qs_env, error )
          CALL set_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error )
          ! keeping only one copy
          CALL scp_env_release ( scp_env, error )
       END IF
    END IF
    CALL timestop(handle)

END SUBROUTINE qs_env_setup

! *****************************************************************************
!> \brief updates the s_mstruct to reflect the new overlap structure,
!>      and also updates rho_core distribution.
!>      Should be called after the atoms have moved and the new overlap
!>      has been calculated.
!> \param qs_env the environment to update
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_env_update_s_mstruct(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_update_s_mstruct', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: do_ppl, failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(para_env, dft_control,particle_set,distribution_2d)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control,&
         particle_set=particle_set, matrix_ks=matrix_ks, matrix_ks_aux_fit=matrix_ks_aux_fit, &
         distribution_2d=distribution_2d,&
         error=error)

    ! *** deallocate matrices that will have the wrong structure ***
    CALL cp_dbcsr_deallocate_matrix_set(matrix_ks,error=error)
    CALL set_qs_env(qs_env,matrix_ks=matrix_ks,error=error)


    IF( dft_control%do_admm) THEN
      ! *** deallocate matrices that will have the wrong structure ***
      CALL cp_dbcsr_deallocate_matrix_set(matrix_ks_aux_fit,error=error)
      CALL set_qs_env(qs_env,matrix_ks_aux_fit=matrix_ks_aux_fit,error=error)
    END IF

    IF(dft_control%qs_control%gapw) THEN
      qs_env%qs_charges%total_rho_core_rspace=qs_env%local_rho_set%rhoz_tot
      IF(dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
        ! *** updates rho core ***
        CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
        CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
        CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
        IF (.NOT. failure) THEN
           CALL calculate_rho_core(qs_env%rho_core, &
                qs_env%qs_charges%total_rho_core_rspace,qs_env,only_nopaw=.TRUE.,error=error)
        END IF
      ELSE
        IF (ASSOCIATED(qs_env%rho_core)) THEN
          CALL pw_release(qs_env%rho_core%pw,error=error)
          DEALLOCATE(qs_env%rho_core)
        ENDIF
      ENDIF
      ! force analytic ppl calculation
      dft_control%qs_control%do_ppl_method=do_ppl_analytic
    ELSE IF(dft_control%qs_control%semi_empirical) THEN
      !??
    ELSE IF(dft_control%qs_control%dftb) THEN
      !??
    ELSE IF(dft_control%qs_control%scptb) THEN
      !??
    ELSE
      ! *** updates rho core ***
      CPPrecondition(ASSOCIATED(qs_env%rho_core),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
      CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure) THEN
         CALL calculate_rho_core(qs_env%rho_core, &
              qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
      END IF
    END IF

    ! calculate local pseudopotential on grid
    do_ppl = dft_control%qs_control%do_ppl_method==do_ppl_grid
    IF ( do_ppl ) THEN
      CPPrecondition(ASSOCIATED(qs_env%vppl),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
      CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure) THEN
         CALL calculate_ppl_grid(qs_env%vppl,qs_env,error=error)
      END IF
    END IF

    ! compute the rho_nlcc
    IF (ASSOCIATED(qs_env%rho_nlcc)) THEN
      CPPrecondition(ASSOCIATED(qs_env%rho_nlcc),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(qs_env%rho_nlcc_g),cp_failure_level,routineP,error,failure)
      CPPrecondition(ASSOCIATED(qs_env%pw_env),cp_failure_level,routineP,error,failure)
      CPPrecondition(qs_env%pw_env%ref_count>0,cp_failure_level,routineP,error,failure)
      CALL calculate_rho_nlcc(qs_env%rho_nlcc, qs_env,error=error)
      CALL pw_transfer(qs_env%rho_nlcc%pw,qs_env%rho_nlcc_g%pw,error=error)
    ENDIF

    ! allocates and creates the task_list
    CALL qs_create_task_list(qs_env, error)

    ! *** environment for ddapc ***
    IF (ASSOCIATED(qs_env%cp_ddapc_env)) THEN
       CALL cp_ddapc_release(qs_env%cp_ddapc_env,error)
    END IF
    CALL cp_ddapc_init(qs_env, error)

    ! *** tell ks_env ***
    IF (ASSOCIATED(qs_env%ks_env)) THEN
       CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE.,&
            error=error)
    END IF
    !   *** Updates rho structure ***
    CALL qs_env_rebuild_rho(qs_env=qs_env,error=error)

    ! *** tell scf_env ***
    IF (ASSOCIATED(qs_env%scf_env)) THEN
       CALL scf_env_did_change(qs_env%scf_env,&
            error=error)
    END IF

    CALL timestop(handle)

END SUBROUTINE qs_env_update_s_mstruct

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
SUBROUTINE qs_create_task_list(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_create_task_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, isub
    LOGICAL                                  :: failure, &
                                                skip_load_balance_distributed,&
                                                soft_valid
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(task_list_type), POINTER            :: task_list

    CALL timeset(routineN,handle)
    failure=.FALSE.

    dft_control => qs_env%dft_control

    skip_load_balance_distributed=dft_control%qs_control%skip_load_balance_distributed
    IF (.NOT. (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb) ) THEN
       ! generate task lists (non-soft)
       IF (.NOT. dft_control%qs_control%gapw) THEN
          CALL get_qs_env(qs_env=qs_env,task_list=task_list,error=error)
          IF (.NOT. ASSOCIATED(task_list)) THEN
             CALL allocate_task_list(task_list,error)
             CALL set_qs_env(qs_env=qs_env,task_list=task_list,error=error)
          ENDIF
          CALL generate_qs_task_list(qs_env, task_list, &
                 reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                 skip_load_balance_distributed=skip_load_balance_distributed, error=error)
       ENDIF
       ! generate the soft task list
       IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
          CALL get_qs_env(qs_env=qs_env,task_list_soft=task_list,error=error)  ! task_list == soft_task_list
          IF (.NOT. ASSOCIATED(task_list)) THEN
             CALL allocate_task_list(task_list,error)
             CALL set_qs_env(qs_env=qs_env,task_list_soft=task_list,error=error)
          ENDIF
          CALL generate_qs_task_list(qs_env, task_list, &
                 reorder_rs_grid_ranks=.TRUE., soft_valid = .TRUE., &
                 skip_load_balance_distributed=skip_load_balance_distributed, error=error)
       ENDIF
    ENDIF

    IF( dft_control%do_admm) THEN
      ! generate the aux_fit task list
      CALL get_qs_env(qs_env=qs_env,task_list_aux_fit=task_list,error=error)
      IF (.NOT. ASSOCIATED(task_list)) THEN
        CALL allocate_task_list(task_list,error)
        CALL set_qs_env(qs_env=qs_env,task_list_aux_fit=task_list,error=error)
      ENDIF
      CALL generate_qs_task_list(qs_env, task_list,&
           reorder_rs_grid_ranks=.FALSE., soft_valid=.FALSE.,basis_set_id=use_aux_fit_basis_set, &
                 skip_load_balance_distributed=skip_load_balance_distributed, error=error)
    END IF

    IF (dft_control%qs_control%do_kg) THEN
      soft_valid = (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc)

      IF(qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN

         IF (ASSOCIATED(qs_env%kg_env%subset)) THEN
           DO isub=1,qs_env%kg_env%nsubsets
             IF (ASSOCIATED(qs_env%kg_env%subset(isub)%task_list)) &
               CALL deallocate_task_list(qs_env%kg_env%subset(isub)%task_list, error)
           END DO
         ELSE
           ALLOCATE(qs_env%kg_env%subset(qs_env%kg_env%nsubsets))
         END IF

         DO isub=1, qs_env%kg_env%nsubsets

           CALL allocate_task_list(qs_env%kg_env%subset(isub)%task_list, error)

           ! generate the subset task list from the neighborlist
           CALL generate_qs_task_list(qs_env, qs_env%kg_env%subset(isub)%task_list, &
                reorder_rs_grid_ranks=.FALSE., soft_valid = soft_valid, &
                skip_load_balance_distributed=skip_load_balance_distributed, &
                sab_orb_external=qs_env%kg_env%subset(isub)%sab_orb, &
                error=error)

         END DO

       END IF

    END IF

    CALL timestop(handle)

END SUBROUTINE qs_create_task_list

! *****************************************************************************
!> \brief rebuilds the pw_env in the given qs_env, allocating it if necessary
!> \param qs_env the qs_env whose pw_env has to be rebuilt
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_env_rebuild_pw_env(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_pw_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure, nlcc
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: new_pw_env
    TYPE(pw_p_type), POINTER                 :: external_vxc, rho_core, &
                                                rho_nlcc, rho_nlcc_g, vee, &
                                                vppl
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ! rebuild pw_env
    NULLIFY(dft_control)
    CALL get_qs_env(qs_env, pw_env=new_pw_env, dft_control=dft_control,error=error)

    IF (.NOT.ASSOCIATED(new_pw_env)) THEN
       CALL pw_env_create(new_pw_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env,pw_env=new_pw_env,error=error)
       CALL pw_env_release(new_pw_env,error=error)
    ELSE
       CALL pw_env_rebuild(new_pw_env,qs_env=qs_env,error=error)
    END IF

    ! reallocate rho_core
    CALL get_qs_env(qs_env, pw_env=new_pw_env, rho_core=rho_core,error=error)
    CPPrecondition(ASSOCIATED(new_pw_env),cp_failure_level,routineP,error,failure)
    IF(dft_control%qs_control%gapw) THEN
      IF (ASSOCIATED(rho_core)) THEN
         CALL pw_release(rho_core%pw,error=error)
         DEALLOCATE(rho_core,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF(dft_control%qs_control%gapw_control%nopaw_as_gpw ) THEN
        ALLOCATE(rho_core,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (.NOT.failure) THEN
           CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
           CALL pw_pool_create_pw(auxbas_pw_pool, rho_core%pw, &
                use_data=COMPLEXDATA1D, error=error)
           rho_core%pw%in_space=RECIPROCALSPACE
           CALL set_qs_env(qs_env, rho_core=rho_core,error=error)
        END IF
      END IF
    ELSE IF (dft_control%qs_control%semi_empirical) THEN
    ELSE IF (dft_control%qs_control%dftb) THEN
    ELSE IF (dft_control%qs_control%scptb) THEN
    ELSE
      IF (ASSOCIATED(rho_core)) THEN
         CALL pw_release(rho_core%pw,error=error)
         DEALLOCATE(rho_core,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(rho_core,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_core%pw, &
              use_data=COMPLEXDATA1D, error=error)
         rho_core%pw%in_space=RECIPROCALSPACE
         CALL set_qs_env(qs_env, rho_core=rho_core,error=error)
      END IF
    END IF

    ! reallocate vppl (realspace grid of local pseudopotential
    IF (dft_control%qs_control%do_ppl_method==do_ppl_grid) THEN
      NULLIFY(vppl)
      CALL get_qs_env(qs_env,pw_env=new_pw_env,vppl=vppl,error=error)
      IF (ASSOCIATED(vppl)) THEN
         CALL pw_release(vppl%pw,error=error)
         DEALLOCATE(vppl,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(vppl,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, vppl%pw, use_data=REALDATA3D, error=error)
         vppl%pw%in_space=REALSPACE
         CALL set_qs_env(qs_env,vppl=vppl,error=error)
      END IF
    END IF

    ! reallocate rho_nlcc
    CALL has_nlcc(nlcc,qs_env,error)
    IF (nlcc) THEN
      ! right now, not working with gapw/gapw_xc, needs implementation in the GAPW XC routines
      IF (dft_control%qs_control%gapw_xc .OR. dft_control%qs_control%gapw) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Non linear core correction for GAPW not implemented",&
               error=error, error_level=cp_failure_level)
      ENDIF
      ! the realspace version
      NULLIFY(rho_nlcc)
      CALL get_qs_env(qs_env,pw_env=new_pw_env,rho_nlcc=rho_nlcc,error=error)
      IF (ASSOCIATED(rho_nlcc)) THEN
         CALL pw_release(rho_nlcc%pw,error=error)
         DEALLOCATE(rho_nlcc,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(rho_nlcc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_nlcc%pw, use_data=REALDATA3D, error=error)
         rho_nlcc%pw%in_space=REALSPACE
         CALL set_qs_env(qs_env,rho_nlcc=rho_nlcc,error=error)
      END IF
      ! the g-space version
      NULLIFY(rho_nlcc_g)
      CALL get_qs_env(qs_env,pw_env=new_pw_env,rho_nlcc_g=rho_nlcc_g,error=error)
      IF (ASSOCIATED(rho_nlcc_g)) THEN
         CALL pw_release(rho_nlcc_g%pw,error=error)
         DEALLOCATE(rho_nlcc_g,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(rho_nlcc_g,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_nlcc_g%pw, use_data=COMPLEXDATA1D, error=error)
         rho_nlcc_g%pw%in_space=RECIPROCALSPACE
         CALL set_qs_env(qs_env,rho_nlcc_g=rho_nlcc_g,error=error)
      END IF
    END IF


    ! reallocate vee: external electrostatic potential
    IF (dft_control%apply_external_potential) THEN
      NULLIFY(vee)
      CALL get_qs_env(qs_env,pw_env=new_pw_env,vee=vee,error=error)
      IF (ASSOCIATED(vee)) THEN
         CALL pw_release(vee%pw,error=error)
         DEALLOCATE(vee,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(vee,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, vee%pw, use_data=REALDATA3D, error=error)
         vee%pw%in_space=REALSPACE
         CALL set_qs_env(qs_env,vee=vee,error=error)
         dft_control%eval_external_potential=.TRUE.
      END IF
    END IF

! ZMP Reallocate external_vxc: external vxc potential
    IF (dft_control%apply_external_vxc) THEN
      NULLIFY(external_vxc)
      CALL get_qs_env(qs_env,pw_env=new_pw_env,external_vxc=external_vxc,error=error)
      IF (ASSOCIATED(external_vxc)) THEN
        CALL pw_release(external_vxc%pw,error=error)
        DEALLOCATE(external_vxc,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      ENDIF
      ALLOCATE(external_vxc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (.NOT.failure) THEN
        CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
        CALL pw_pool_create_pw(auxbas_pw_pool, external_vxc%pw, use_data=REALDATA3D, error=error)
        external_vxc%pw%in_space=REALSPACE
        CALL set_qs_env(qs_env,external_vxc=external_vxc,error=error)
        dft_control%read_external_vxc=.TRUE.
      END IF
    END IF

    CALL timestop(handle)

END SUBROUTINE qs_env_rebuild_pw_env

! *****************************************************************************
!> \brief finds if a given qs run needs to use nlcc
!> \param nlcc ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE has_nlcc(nlcc,qs_env,error)

    LOGICAL                                  :: nlcc
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'has_nlcc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind
    LOGICAL                                  :: nlcc_present
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gth_potential_type), POINTER        :: gth_potential

    nlcc=.FALSE.

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)
    DO ikind=1,SIZE(atomic_kind_set)
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind, gth_potential=gth_potential)
      IF (.NOT.ASSOCIATED(gth_potential)) CYCLE
      CALL get_potential(potential=gth_potential,nlcc_present=nlcc_present)
      nlcc=nlcc.OR.nlcc_present
    ENDDO

  END SUBROUTINE has_nlcc

! *****************************************************************************
!> \brief rebuilds the rho structure, making sure that everything is allocated
!>      and has the right size
!> \param qs_env the environment in which rho should be rebuilt
!> \param rebuild_ao if it is necessary to rebuild rho_ao. Defaults to true.
!> \param rebuild_grids if it in necessary to rebuild rho_r and rho_g.
!>        Defaults to false.
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs updated  pw pools, s_mstruct and h.
!>      The use of p to keep the structure of h (needed for the forces)
!>      is ugly and should be removed.
!>      If necessary rho is created from scratch.
! *****************************************************************************
SUBROUTINE qs_env_rebuild_rho(qs_env, rebuild_ao, rebuild_grids, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rebuild_ao, rebuild_grids
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_rho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: do_admm, gapw_xc
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env, dft_control=dft_control,error=error)
    gapw_xc=dft_control%qs_control%gapw_xc
    do_admm = dft_control%do_admm
    CALL qs_rho_rebuild(qs_env%rho,qs_env=qs_env,&
         rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, error=error)
    IF(gapw_xc) THEN
      CALL qs_rho_rebuild(qs_env%rho_xc,qs_env=qs_env,&
         rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, gapw_xc=gapw_xc, error=error)
      qs_env%rho_xc%rho_ao => qs_env%rho%rho_ao
    END IF
    IF(do_admm) THEN
      CALL qs_rho_rebuild(qs_env%rho_aux_fit,qs_env=qs_env,&
                          rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, &
                          basis_set_id=use_aux_fit_basis_set, error=error)
    END IF
! ZMP rebuilding external density
    IF (dft_control%apply_external_density) THEN
      CALL qs_rho_rebuild(qs_env%rho_external,qs_env=qs_env,&
                           rebuild_grids=rebuild_grids, &
                           error=error)
      dft_control%read_external_density=.TRUE.
    ENDIF


    CALL timestop(handle)

END SUBROUTINE qs_env_rebuild_rho

END MODULE qs_environment_methods
