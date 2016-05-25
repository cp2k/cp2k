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
MODULE qs_update_s_mstruct
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_ddapc_types,                  ONLY: cp_ddapc_release
   USE cp_ddapc_util,                   ONLY: cp_ddapc_init
   USE input_constants,                 ONLY: do_ppl_analytic,&
                                              do_ppl_grid,&
                                              kg_tnadd_embed
   USE pw_methods,                      ONLY: pw_transfer
   USE pw_types,                        ONLY: pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: calculate_ppl_grid,&
                                              calculate_rho_core,&
                                              calculate_rho_nlcc
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_type
   USE qs_scf_types,                    ONLY: scf_env_did_change
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_update_s_mstruct'

   PUBLIC :: qs_env_update_s_mstruct
!***
CONTAINS

! *****************************************************************************
!> \brief updates the s_mstruct to reflect the new overlap structure,
!>      and also updates rho_core distribution.
!>      Should be called after the atoms have moved and the new overlap
!>      has been calculated.
!> \param qs_env the environment to update
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE qs_env_update_s_mstruct(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_update_s_mstruct', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_ppl
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_p_type), POINTER                           :: rho_core, rho_nlcc, rho_nlcc_g, vppl

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control)

      ! *** updates rho core ***
      NULLIFY (rho_core)
      CALL get_qs_env(qs_env, rho_core=rho_core)
      IF (dft_control%qs_control%gapw) THEN
         qs_env%qs_charges%total_rho_core_rspace = qs_env%local_rho_set%rhoz_tot
         IF (dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
            CPASSERT(ASSOCIATED(rho_core))
            CALL calculate_rho_core(rho_core, &
                                    qs_env%qs_charges%total_rho_core_rspace, qs_env, only_nopaw=.TRUE.)
         ELSE
            IF (ASSOCIATED(rho_core)) THEN
               CALL pw_release(rho_core%pw)
               DEALLOCATE (rho_core)
            ENDIF
         ENDIF
         ! force analytic ppl calculation
         dft_control%qs_control%do_ppl_method = do_ppl_analytic
      ELSE IF (dft_control%qs_control%semi_empirical) THEN
         !??
      ELSE IF (dft_control%qs_control%dftb) THEN
         !??
      ELSE IF (dft_control%qs_control%scptb) THEN
         !??
      ELSE
         CPASSERT(ASSOCIATED(rho_core))
         CALL calculate_rho_core(rho_core, &
                                 qs_env%qs_charges%total_rho_core_rspace, qs_env)
      END IF

      ! calculate local pseudopotential on grid
      do_ppl = dft_control%qs_control%do_ppl_method == do_ppl_grid
      IF (do_ppl) THEN
         NULLIFY (vppl)
         CALL get_qs_env(qs_env, vppl=vppl)
         CPASSERT(ASSOCIATED(vppl))
         CALL calculate_ppl_grid(vppl, qs_env)
      END IF

      ! compute the rho_nlcc
      NULLIFY (rho_nlcc, rho_nlcc_g)
      CALL get_qs_env(qs_env, rho_nlcc=rho_nlcc, rho_nlcc_g=rho_nlcc_g)
      IF (ASSOCIATED(rho_nlcc)) THEN
         CALL calculate_rho_nlcc(rho_nlcc, qs_env)
         CALL pw_transfer(rho_nlcc%pw, rho_nlcc_g%pw)
      ENDIF

      ! allocates and creates the task_list
      CALL qs_create_task_list(qs_env)

      ! *** environment for ddapc ***
      IF (ASSOCIATED(qs_env%cp_ddapc_env)) THEN
         CALL cp_ddapc_release(qs_env%cp_ddapc_env)
      END IF
      CALL cp_ddapc_init(qs_env)

      ! *** tell ks_env ***
      CALL qs_ks_did_change(qs_env%ks_env, s_mstruct_changed=.TRUE.)

      !   *** Updates rho structure ***
      CALL qs_env_rebuild_rho(qs_env=qs_env)

      ! *** tell scf_env ***
      IF (ASSOCIATED(qs_env%scf_env)) THEN
         CALL scf_env_did_change(qs_env%scf_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_env_update_s_mstruct

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE qs_create_task_list(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_create_task_list', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, isub
      LOGICAL                                            :: skip_load_balance_distributed, soft_valid
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list

      CALL timeset(routineN, handle)
      NULLIFY (ks_env, dft_control)
      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control)

      skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
      IF (.NOT. (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb)) THEN
         ! generate task lists (non-soft)
         IF (.NOT. dft_control%qs_control%gapw) THEN
            CALL get_ks_env(ks_env, task_list=task_list)
            IF (.NOT. ASSOCIATED(task_list)) THEN
               CALL allocate_task_list(task_list)
               CALL set_ks_env(ks_env, task_list=task_list)
            ENDIF
            CALL generate_qs_task_list(ks_env, task_list, &
                                       reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                                       skip_load_balance_distributed=skip_load_balance_distributed)
         ENDIF
         ! generate the soft task list
         IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
            CALL get_ks_env(ks_env, task_list_soft=task_list)
            IF (.NOT. ASSOCIATED(task_list)) THEN
               CALL allocate_task_list(task_list)
               CALL set_ks_env(ks_env, task_list_soft=task_list)
            ENDIF
            CALL generate_qs_task_list(ks_env, task_list, &
                                       reorder_rs_grid_ranks=.TRUE., soft_valid=.TRUE., &
                                       skip_load_balance_distributed=skip_load_balance_distributed)
         ENDIF
      ENDIF

      IF (dft_control%do_admm) THEN
         ! generate the aux_fit task list
         CALL get_ks_env(ks_env, task_list_aux_fit=task_list)
         IF (.NOT. ASSOCIATED(task_list)) THEN
            CALL allocate_task_list(task_list)
            CALL set_ks_env(ks_env, task_list_aux_fit=task_list)
         ENDIF
         CALL get_ks_env(ks_env, sab_aux_fit=sab_orb)
         CALL generate_qs_task_list(ks_env, task_list, &
                                    reorder_rs_grid_ranks=.FALSE., soft_valid=.FALSE., basis_type="AUX_FIT", &
                                    skip_load_balance_distributed=skip_load_balance_distributed, &
                                    sab_orb_external=sab_orb)
      END IF

      IF (dft_control%qs_control%do_kg) THEN
         soft_valid = (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc)

         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN

            IF (ASSOCIATED(qs_env%kg_env%subset)) THEN
               DO isub = 1, qs_env%kg_env%nsubsets
                  IF (ASSOCIATED(qs_env%kg_env%subset(isub)%task_list)) &
                     CALL deallocate_task_list(qs_env%kg_env%subset(isub)%task_list)
               END DO
            ELSE
               ALLOCATE (qs_env%kg_env%subset(qs_env%kg_env%nsubsets))
            END IF

            DO isub = 1, qs_env%kg_env%nsubsets

               CALL allocate_task_list(qs_env%kg_env%subset(isub)%task_list)

               ! generate the subset task list from the neighborlist
               CALL generate_qs_task_list(ks_env, qs_env%kg_env%subset(isub)%task_list, &
                                          reorder_rs_grid_ranks=.FALSE., soft_valid=soft_valid, &
                                          skip_load_balance_distributed=skip_load_balance_distributed, &
                                          sab_orb_external=qs_env%kg_env%subset(isub)%sab_orb)

            END DO

         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_create_task_list

! *****************************************************************************
!> \brief rebuilds the rho structure, making sure that everything is allocated
!>      and has the right size
!> \param qs_env the environment in which rho should be rebuilt
!> \param rebuild_ao if it is necessary to rebuild rho_ao. Defaults to true.
!> \param rebuild_grids if it in necessary to rebuild rho_r and rho_g.
!>        Defaults to false.
!> \par History
!>      10.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      needs updated  pw pools, s_mstruct and h.
!>      The use of p to keep the structure of h (needed for the forces)
!>      is ugly and should be removed.
!>      If necessary rho is created from scratch.
! **************************************************************************************************
   SUBROUTINE qs_env_rebuild_rho(qs_env, rebuild_ao, rebuild_grids)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in), OPTIONAL                      :: rebuild_ao, rebuild_grids

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_rebuild_rho', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_admm, gapw_xc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit, rho_aux_fit_buffer, &
                                                            rho_external, rho_xc

      NULLIFY (rho)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      rho_xc=rho_xc, &
                      rho_aux_fit=rho_aux_fit, &
                      rho_aux_fit_buffer=rho_aux_fit_buffer, &
                      rho_external=rho_external)

      gapw_xc = dft_control%qs_control%gapw_xc
      do_admm = dft_control%do_admm
      CALL qs_rho_rebuild(rho, qs_env=qs_env, &
                          rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids)

      IF (gapw_xc) THEN
         CALL qs_rho_rebuild(rho_xc, qs_env=qs_env, &
                             rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids)
      END IF
      IF (do_admm) THEN
         CALL qs_rho_rebuild(rho_aux_fit, qs_env=qs_env, &
                             rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, &
                             admm=.TRUE.)
         CALL qs_rho_rebuild(rho_aux_fit_buffer, qs_env=qs_env, &
                             rebuild_ao=rebuild_ao, rebuild_grids=rebuild_grids, &
                             admm=.TRUE.)
      END IF
! ZMP rebuilding external density
      IF (dft_control%apply_external_density) THEN
         CALL qs_rho_rebuild(rho_external, qs_env=qs_env, &
                             rebuild_grids=rebuild_grids)
         dft_control%read_external_density = .TRUE.
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE qs_env_rebuild_rho

END MODULE qs_update_s_mstruct
