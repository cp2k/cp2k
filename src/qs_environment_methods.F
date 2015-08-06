!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief qs_environement methods that use many other modules
!> \par History
!>      09.2002 created [fawzi]
!>      - local atom distribution (25.06.2003,MK)
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE qs_environment_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: dbcsr_distribution_obj
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
  USE cp_ddapc_types,                  ONLY: cp_ddapc_release
  USE cp_ddapc_util,                   ONLY: cp_ddapc_init
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_2d
  USE ewald_environment_types,         ONLY: ewald_environment_type
  USE ewald_pw_methods,                ONLY: ewald_pw_grid_update
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE external_potential_types,        ONLY: get_potential,&
                                             gth_potential_type
  USE input_constants,                 ONLY: do_ppl_analytic,&
                                             do_ppl_grid,&
                                             kg_tnadd_embed
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
                                             pw_release,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_ppl_grid,&
                                             calculate_rho_core,&
                                             calculate_rho_nlcc
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_did_change,&
                                             qs_ks_env_type,&
                                             set_ks_env
  USE qs_matrix_pools,                 ONLY: mpools_rebuild_fm_pools
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_outer_scf,                    ONLY: outer_loop_variables_count
  USE qs_rho0_types,                   ONLY: rho0_mpole_type
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: scf_env_did_change
  USE scf_control_types,               ONLY: scf_control_type
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
!> \param unit_nr ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_env_setup(qs_env,unit_nr,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_ao_aux_fit, &
                                                nhistory, nvariables, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: outer_scf_history
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
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
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(qs_kind_set, atomic_kind_set, dft_control, scf_control, qs_charges, para_env,&
            distribution_2d, molecule_kind_set, molecule_set, particle_set, cell, &
            ks_env, blacs_env)

    CALL get_qs_env(qs_env=qs_env,&
                    qs_kind_set=qs_kind_set,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set,&
                    scf_control=scf_control,&
                    para_env=para_env,&
                    blacs_env=blacs_env,&
                    cell=cell,&
                    ks_env=ks_env,&
                    error=error)

    CPPrecondition(ASSOCIATED(qs_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(atomic_kind_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_control),cp_failure_level,routineP,error,failure)
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
               blacs_env=blacs_env, para_env=para_env,&
               error=error)
    ENDIF

    ! If we use auxiliary density matrix methods rebuild fm_pools
    IF(dft_control%do_admm) THEN
      CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao_aux_fit,&
                               basis_type="AUX_FIT",error=error)
      CALL mpools_rebuild_fm_pools(qs_env%mpools_aux_fit,mos=qs_env%mos_aux_fit,&
           blacs_env=blacs_env, para_env=para_env,&
           error=error)
     END IF

    ! create 2d distribution

    CALL distribute_molecules_2d(cell=cell,&
                                 atomic_kind_set=atomic_kind_set,&
                                 qs_kind_set=qs_kind_set,&
                                 particle_set=particle_set,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 distribution_2d=distribution_2d,&
                                 blacs_env=blacs_env,&
                                 force_env_section=qs_env%input, error=error)

    ! and use it to create the dbcsr_dist, which should be the sole user of distribution_2d by now.
    ALLOCATE(dbcsr_dist)
    CALL cp_dbcsr_dist2d_to_dist(distribution_2d, dbcsr_dist, unit_nr, error)
    CALL set_ks_env(ks_env, dbcsr_dist=dbcsr_dist, error=error)

    ! also keep distribution_2d in qs_env
    CALL set_ks_env(ks_env, distribution_2d=distribution_2d, error=error)
    CALL distribution_2d_release(distribution_2d, error=error)

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

    INTEGER                                  :: handle, istat
    LOGICAL                                  :: do_ppl, failure
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), POINTER                 :: rho_core, rho_nlcc, &
                                                rho_nlcc_g, vppl

    CALL timeset(routineN,handle)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    NULLIFY(dft_control)
    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    error=error)

    ! *** updates rho core ***
    NULLIFY(rho_core)
    CALL get_qs_env(qs_env,rho_core=rho_core,error=error)
    IF(dft_control%qs_control%gapw) THEN
       qs_env%qs_charges%total_rho_core_rspace=qs_env%local_rho_set%rhoz_tot
       IF(dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
          CPPrecondition(ASSOCIATED(rho_core),cp_failure_level,routineP,error,failure)
          CALL calculate_rho_core(rho_core, &
               qs_env%qs_charges%total_rho_core_rspace,qs_env,only_nopaw=.TRUE.,error=error)
       ELSE
          IF (ASSOCIATED(rho_core)) THEN
             CALL pw_release(rho_core%pw,error=error)
             DEALLOCATE(rho_core,STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
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
       CPPrecondition(ASSOCIATED(rho_core),cp_failure_level,routineP,error,failure)
       CALL calculate_rho_core(rho_core, &
            qs_env%qs_charges%total_rho_core_rspace,qs_env,error=error)
    END IF

    ! calculate local pseudopotential on grid
    do_ppl = dft_control%qs_control%do_ppl_method==do_ppl_grid
    IF ( do_ppl ) THEN
       NULLIFY(vppl)
       CALL get_qs_env(qs_env,vppl=vppl,error=error)
       CPPrecondition(ASSOCIATED(vppl),cp_failure_level,routineP,error,failure)
       CALL calculate_ppl_grid(vppl,qs_env,error=error)
    END IF

    ! compute the rho_nlcc
    NULLIFY(rho_nlcc,rho_nlcc_g)
    CALL get_qs_env(qs_env,rho_nlcc=rho_nlcc,rho_nlcc_g=rho_nlcc_g,error=error)
    IF (ASSOCIATED(rho_nlcc)) THEN
       CALL calculate_rho_nlcc(rho_nlcc,qs_env,error=error)
       CALL pw_transfer(rho_nlcc%pw,rho_nlcc_g%pw,error=error)
    ENDIF

    ! allocates and creates the task_list
    CALL qs_create_task_list(qs_env,error)

    ! *** environment for ddapc ***
    IF (ASSOCIATED(qs_env%cp_ddapc_env)) THEN
       CALL cp_ddapc_release(qs_env%cp_ddapc_env,error)
    END IF
    CALL cp_ddapc_init(qs_env,error)

    ! *** tell ks_env ***
    CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE.,error=error)

    !   *** Updates rho structure ***
    CALL qs_env_rebuild_rho(qs_env=qs_env,error=error)

    ! *** tell scf_env ***
    IF (ASSOCIATED(qs_env%scf_env)) THEN
       CALL scf_env_did_change(qs_env%scf_env,error=error)
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
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(task_list_type), POINTER            :: task_list

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(ks_env, dft_control)
    CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, error=error)

    skip_load_balance_distributed=dft_control%qs_control%skip_load_balance_distributed
    IF (.NOT. (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb) ) THEN
       ! generate task lists (non-soft)
       IF (.NOT. dft_control%qs_control%gapw) THEN
          CALL get_ks_env(ks_env, task_list=task_list, error=error)
          IF (.NOT. ASSOCIATED(task_list)) THEN
             CALL allocate_task_list(task_list,error)
             CALL set_ks_env(ks_env, task_list=task_list, error=error)
          ENDIF
          CALL generate_qs_task_list(ks_env, task_list, &
                 reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                 skip_load_balance_distributed=skip_load_balance_distributed, error=error)
       ENDIF
       ! generate the soft task list
       IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
          CALL get_ks_env(ks_env, task_list_soft=task_list, error=error)
          IF (.NOT. ASSOCIATED(task_list)) THEN
             CALL allocate_task_list(task_list,error)
             CALL set_ks_env(ks_env, task_list_soft=task_list, error=error)
          ENDIF
          CALL generate_qs_task_list(ks_env, task_list, &
                 reorder_rs_grid_ranks=.TRUE., soft_valid = .TRUE., &
                 skip_load_balance_distributed=skip_load_balance_distributed, error=error)
       ENDIF
    ENDIF

    IF( dft_control%do_admm) THEN
      ! generate the aux_fit task list
      CALL get_ks_env(ks_env, task_list_aux_fit=task_list, error=error)
      IF (.NOT. ASSOCIATED(task_list)) THEN
        CALL allocate_task_list(task_list,error)
        CALL set_ks_env(ks_env, task_list_aux_fit=task_list, error=error)
      ENDIF
      CALL get_ks_env(ks_env,sab_aux_fit=sab_orb,error=error)
      CALL generate_qs_task_list(ks_env, task_list,&
             reorder_rs_grid_ranks=.FALSE.,soft_valid=.FALSE.,basis_type="AUX_FIT", &
             skip_load_balance_distributed=skip_load_balance_distributed, &
             sab_orb_external=sab_orb,error=error)
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
           CALL generate_qs_task_list(ks_env, qs_env%kg_env%subset(isub)%task_list, &
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
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(pw_env_type), POINTER               :: new_pw_env
    TYPE(pw_p_type), POINTER                 :: external_vxc, rho0_gs, &
                                                rho0_rs, rho_core, rho_nlcc, &
                                                rho_nlcc_g, vee, vppl
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER                   :: v_hartree_rspace
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole

    CALL timeset(routineN,handle)
    failure=.FALSE.
    ! rebuild pw_env
    NULLIFY(dft_control, cell, ks_env, v_hartree_rspace, auxbas_pw_pool)
    NULLIFY(rho0_mpole,rho0_gs,rho0_rs)
    NULLIFY(ewald_env,ewald_pw,new_pw_env,external_vxc,rho_core,rho_nlcc,rho_nlcc_g,vee,vppl)

    CALL get_qs_env(qs_env, ks_env=ks_env, pw_env=new_pw_env, error=error)
    IF (.NOT.ASSOCIATED(new_pw_env)) THEN
       CALL pw_env_create(new_pw_env,error=error)
       CALL set_ks_env(ks_env, pw_env=new_pw_env, error=error)
       CALL pw_env_release(new_pw_env,error=error)
    ENDIF

    CALL get_qs_env(qs_env, pw_env=new_pw_env, dft_control=dft_control,&
                    cell=cell, error=error)

    IF(ANY(new_pw_env%cell_hmat/=cell%hmat)) THEN
       ! only rebuild if nessecary
       new_pw_env%cell_hmat = cell%hmat
       CALL pw_env_rebuild(new_pw_env,qs_env=qs_env,error=error)

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
           CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
           CALL pw_pool_create_pw(auxbas_pw_pool, rho_core%pw, &
                use_data=COMPLEXDATA1D, error=error)
           rho_core%pw%in_space=RECIPROCALSPACE
           CALL set_ks_env(ks_env, rho_core=rho_core,error=error)
         END IF
         CALL get_qs_env(qs_env=qs_env,rho0_mpole=rho0_mpole,error=error)
         IF(ASSOCIATED(rho0_mpole)) THEN
           rho0_rs => rho0_mpole%rho0_s_rs
           rho0_gs => rho0_mpole%rho0_s_gs
           IF(ASSOCIATED(rho0_rs)) THEN
              CALL pw_release(rho0_rs%pw,error=error)
              DEALLOCATE(rho0_rs,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           ALLOCATE(rho0_rs,stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
           CALL pw_pool_create_pw(auxbas_pw_pool, rho0_rs%pw,&
           use_data=REALDATA3D,in_space=REALSPACE,error=error)
           rho0_mpole%rho0_s_rs => rho0_rs
   
           IF(ASSOCIATED(rho0_gs)) THEN
              CALL pw_release(rho0_gs%pw,error=error)
              DEALLOCATE(rho0_gs,stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
           ALLOCATE(rho0_gs,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
           CALL pw_pool_create_pw(auxbas_pw_pool, rho0_gs%pw,&
           use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
           rho0_mpole%rho0_s_gs => rho0_gs
         END IF
       ELSE IF (dft_control%qs_control%semi_empirical) THEN
          IF (dft_control%qs_control%se_control%do_ewald .OR. &
              dft_control%qs_control%se_control%do_ewald_gks) THEN
             ! rebuild Ewald environment
             CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
             CALL ewald_pw_grid_update(ewald_pw,ewald_env,cell%hmat,error)
          END IF
       ELSE IF (dft_control%qs_control%dftb) THEN
          IF (dft_control%qs_control%dftb_control%do_ewald) THEN
             ! rebuild Ewald environment
             CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
             CALL ewald_pw_grid_update(ewald_pw,ewald_env,cell%hmat,error)
          END IF
       ELSE IF (dft_control%qs_control%scptb) THEN
          IF (dft_control%qs_control%scptb_control%do_ewald) THEN
             ! rebuild Ewald environment
             CALL get_qs_env(qs_env=qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
             CALL ewald_pw_grid_update(ewald_pw,ewald_env,cell%hmat,error)
          END IF
       ELSE
         IF (ASSOCIATED(rho_core)) THEN
            CALL pw_release(rho_core%pw,error=error)
            DEALLOCATE(rho_core,stat=stat)
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
         ENDIF
         ALLOCATE(rho_core,stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         CALL pw_env_get(new_pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_core%pw, &
              use_data=COMPLEXDATA1D, error=error)
         rho_core%pw%in_space=RECIPROCALSPACE
         CALL set_ks_env(ks_env, rho_core=rho_core,error=error)
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
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, vppl%pw, use_data=REALDATA3D, error=error)
         vppl%pw%in_space=REALSPACE
         CALL set_ks_env(ks_env,vppl=vppl,error=error)
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
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_nlcc%pw, use_data=REALDATA3D, error=error)
         rho_nlcc%pw%in_space=REALSPACE
         CALL set_ks_env(ks_env,rho_nlcc=rho_nlcc,error=error)
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
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_nlcc_g%pw, use_data=COMPLEXDATA1D, error=error)
         rho_nlcc_g%pw%in_space=RECIPROCALSPACE
         CALL set_ks_env(ks_env,rho_nlcc_g=rho_nlcc_g,error=error)
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
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, vee%pw, use_data=REALDATA3D, error=error)
         vee%pw%in_space=REALSPACE
         CALL set_ks_env(ks_env,vee=vee,error=error)
         dft_control%eval_external_potential=.TRUE.
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
         CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool, external_vxc%pw, use_data=REALDATA3D, error=error)
         external_vxc%pw%in_space=REALSPACE
         CALL set_qs_env(qs_env,external_vxc=external_vxc,error=error)
         dft_control%read_external_vxc=.TRUE.
       END IF

       CALL get_ks_env(ks_env, v_hartree_rspace=v_hartree_rspace, error=error)
       IF (ASSOCIATED(v_hartree_rspace)) &
          CALL pw_release(v_hartree_rspace, error=error)
       CALL get_qs_env(qs_env,pw_env=new_pw_env,error=error)
       CALL pw_env_get(new_pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace,&
            use_data=REALDATA3D, in_space=REALSPACE, error=error)
       CALL set_ks_env(ks_env, v_hartree_rspace=v_hartree_rspace, error=error)
    ENDIF

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
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    nlcc=.FALSE.

    CALL get_qs_env(qs_env=qs_env,qs_kind_set=qs_kind_set,error=error)
    DO ikind=1,SIZE(qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential,error=error)
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
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit, &
                                                rho_aux_fit_buffer, &
                                                rho_external, rho_xc

    NULLIFY(rho)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    rho=rho,&
                    rho_xc=rho_xc,&
                    rho_aux_fit=rho_aux_fit,&
                    rho_aux_fit_buffer=rho_aux_fit_buffer,&
                    rho_external=rho_external,&
                    error=error)

    gapw_xc=dft_control%qs_control%gapw_xc
    do_admm = dft_control%do_admm
    CALL qs_rho_rebuild(rho,qs_env=qs_env,&
         rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, error=error)

    IF(gapw_xc) THEN
      CALL qs_rho_rebuild(rho_xc,qs_env=qs_env,&
         rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, error=error)
    END IF
    IF(do_admm) THEN
      CALL qs_rho_rebuild(rho_aux_fit,qs_env=qs_env,&
                          rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, &
                          admm=.TRUE., error=error)
      CALL qs_rho_rebuild(rho_aux_fit_buffer,qs_env=qs_env,&
                          rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, &
                          admm=.TRUE., error=error)
    END IF
! ZMP rebuilding external density
    IF (dft_control%apply_external_density) THEN
      CALL qs_rho_rebuild(rho_external,qs_env=qs_env,&
                           rebuild_grids=rebuild_grids, error=error)
      dft_control%read_external_density=.TRUE.
    ENDIF


    CALL timestop(handle)

END SUBROUTINE qs_env_rebuild_rho

END MODULE qs_environment_methods
