!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility routines for qs_scf
! **************************************************************************************************
MODULE qs_scf_initialization
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_row_template,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_dbcsr_output,                 ONLY: write_fm_with_basis_info
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_get_el_struct
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_triangular,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_real_default
   USE input_constants,                 ONLY: &
        broy_mix, broy_mix_new, cholesky_dbcsr, cholesky_inverse, cholesky_off, &
        diag_block_davidson, diag_block_krylov, diag_filter_matrix, diag_ot, diag_standard, &
        direct_p_mix, kerker_mix, multisec_mix, no_mix, ot2cdft, outer_scf_none, plus_u_lowdin, &
        pulay_mix, wfi_frozen_method_nr, wfi_ps_method_nr, wfi_use_guess_method_nr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE pw_types,                        ONLY: pw_p_type
   USE qmmm_image_charge,               ONLY: conditional_calc_image_matrix
   USE qs_block_davidson_types,         ONLY: block_davidson_allocate,&
                                              block_davidson_env_create
   USE qs_cdft_opt_types,               ONLY: cdft_opt_type_copy
   USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                              mixing_storage_create,&
                                              mixing_storage_release,&
                                              no_mixing_nr
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_fb_distribution_methods,      ONLY: fb_distribution_build
   USE qs_fb_env_methods,               ONLY: fb_env_build_atomic_halos,&
                                              fb_env_build_rcut_auto,&
                                              fb_env_read_input,&
                                              fb_env_write_info
   USE qs_fb_env_types,                 ONLY: fb_env_create,&
                                              fb_env_has_data
   USE qs_initial_guess,                ONLY: calculate_first_density_matrix
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type,&
                                              set_qs_kind
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_mixing_utils,                 ONLY: charge_mixing_init,&
                                              mixing_allocate,&
                                              mixing_init
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type
   USE qs_outer_scf,                    ONLY: outer_loop_extrapolate,&
                                              outer_loop_switch,&
                                              outer_loop_variables_count
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_diagonalization,          ONLY: diag_subspace_allocate
   USE qs_scf_lanczos,                  ONLY: krylov_space_allocate
   USE qs_scf_output,                   ONLY: qs_scf_initial_info
   USE qs_scf_types,                    ONLY: &
        block_davidson_diag_method_nr, block_krylov_diag_method_nr, diag_subspace_env_create, &
        filter_matrix_diag_method_nr, general_diag_method_nr, krylov_space_create, &
        ot_diag_method_nr, ot_method_nr, qs_scf_env_type, scf_env_create, special_diag_method_nr
   USE qs_wf_history_methods,           ONLY: reorthogonalize_vectors,&
                                              wfi_extrapolate,&
                                              wfi_get_method_label,&
                                              wfi_update
   USE scf_control_types,               ONLY: scf_control_type
   USE xas_env_types,                   ONLY: xas_environment_type
   USE xas_restart,                     ONLY: xas_initialize_rho
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_initialization'

   PUBLIC:: qs_scf_env_initialize, qs_scf_env_init_basic

CONTAINS

! **************************************************************************************************
!> \brief initializes input parameters if needed or restores values from
!>        previous runs to fill scf_env with the values required for scf
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param scf_env ...
!> \param scf_control ...
!> \param scf_section ...
! **************************************************************************************************
   SUBROUTINE qs_scf_env_initialize(qs_env, scf_env, scf_control, scf_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), OPTIONAL, POINTER          :: scf_control
      TYPE(section_vals_type), OPTIONAL, POINTER         :: scf_section

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_env_initialize', &
         routineP = moduleN//':'//routineN

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(scf_control_type), POINTER                    :: my_scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, my_scf_section

      CALL get_qs_env(qs_env, input=input, dft_control=dft_control)

      IF (PRESENT(scf_control)) THEN
         my_scf_control => scf_control
      ELSE
         CALL get_qs_env(qs_env, scf_control=my_scf_control)
      ENDIF

      dft_section => section_vals_get_subs_vals(input, "DFT")
      IF (PRESENT(scf_section)) THEN
         my_scf_section => scf_section
      ELSE
         my_scf_section => section_vals_get_subs_vals(dft_section, "SCF")
      END IF

      CALL qs_scf_ensure_scf_env(qs_env, scf_env)

      CALL section_vals_val_get(my_scf_section, "CHOLESKY", i_val=scf_env%cholesky_method)

      CALL qs_scf_ensure_mos(qs_env, scf_env)

      ! set flags for diagonalization
      CALL qs_scf_ensure_diagonalization(scf_env, my_scf_section, qs_env, &
                                         my_scf_control, qs_env%has_unit_metric)
      ! set parameters for mixing/DIIS during scf
      CALL qs_scf_ensure_mixing(my_scf_control, my_scf_section, scf_env, dft_control)

      CALL qs_scf_ensure_work_matrices(qs_env, scf_env)

      CALL qs_scf_ensure_mixing_store(qs_env, scf_env)

      ! Initialize outer loop variables: handle CDFT and regular outer loop separately
      IF (dft_control%qs_control%cdft) THEN
         CALL qs_scf_ensure_cdft_loop_vars(qs_env, scf_env, dft_control, &
                                           scf_control=my_scf_control)
      ELSE
         CALL qs_scf_ensure_outer_loop_vars(scf_env, my_scf_control)
      END IF

      CALL init_scf_run(scf_env, qs_env, my_scf_section, my_scf_control)

   END SUBROUTINE qs_scf_env_initialize

! **************************************************************************************************
!> \brief initializes input parameters if needed for non-scf calclulations using diagonalization
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_env_init_basic(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_env_init_basic', &
         routineP = moduleN//':'//routineN

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      CALL get_qs_env(qs_env, input=input, dft_control=dft_control)

      CALL get_qs_env(qs_env, scf_control=scf_control)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      CALL qs_scf_ensure_scf_env(qs_env, scf_env)

      CALL section_vals_val_get(scf_section, "CHOLESKY", i_val=scf_env%cholesky_method)
      scf_control%use_diag = .TRUE.
      scf_control%diagonalization%method = diag_standard

      CALL qs_scf_ensure_mos(qs_env, scf_env)

      ! set flags for diagonalization
      CALL qs_scf_ensure_diagonalization(scf_env, scf_section, qs_env, &
                                         scf_control, qs_env%has_unit_metric)
      CALL qs_scf_ensure_work_matrices(qs_env, scf_env)

      CALL init_scf_run(scf_env, qs_env, scf_section, scf_control)

   END SUBROUTINE qs_scf_env_init_basic

! **************************************************************************************************
!> \brief makes sure scf_env is allocated (might already be from before)
!>        in case it is present the g-space mixing storage is reset
!> \param qs_env ...
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_scf_env(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_scf_env', &
         routineP = moduleN//':'//routineN

      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g
      TYPE(qs_rho_type), POINTER                         :: rho

      NULLIFY (rho_g)

      IF (.NOT. ASSOCIATED(scf_env)) THEN ! i.e. for MD this is associated on the second step (it so seems)
         CALL scf_env_create(scf_env)
      ELSE
         ! Reallocate mixing store, if the g space grid (cell) has changed
         SELECT CASE (scf_env%mixing_method)
         CASE (kerker_mix, pulay_mix, broy_mix, broy_mix_new, multisec_mix)
            IF (ASSOCIATED(scf_env%mixing_store)) THEN
               ! The current mixing_store data structure does not allow for an unique
               ! grid comparison, but the probability that the 1d lengths of the old and
               ! the new grid are accidentily equal is rather low
               CALL get_qs_env(qs_env, rho=rho)
               CALL qs_rho_get(rho, rho_g=rho_g)
               IF (ASSOCIATED(scf_env%mixing_store%rhoin)) THEN
                  IF (SIZE(rho_g(1)%pw%pw_grid%gsq) /= SIZE(scf_env%mixing_store%rhoin(1)%cc)) &
                     CALL mixing_storage_release(scf_env%mixing_store)
               END IF
            END IF
         END SELECT
      END IF

   END SUBROUTINE qs_scf_ensure_scf_env

! **************************************************************************************************
!> \brief performs allocation of outer SCF variables
!> \param scf_env the SCF environment which contains the outer SCF variables
!> \param scf_control control settings for the outer SCF loop
!> \param nvar (optional) set number of outer SCF variables externally if CDFT SCF is active
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_outer_loop_vars(scf_env, scf_control, nvar)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      INTEGER, OPTIONAL                                  :: nvar

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_outer_loop_vars', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nhistory, nvariables

      IF (scf_control%outer_scf%have_scf) THEN
         nhistory = scf_control%outer_scf%max_scf + 1
         IF (PRESENT(nvar)) THEN
            IF (nvar > 0) THEN
               nvariables = nvar
            ELSE
               nvariables = outer_loop_variables_count(scf_control)
            END IF
         ELSE
            nvariables = outer_loop_variables_count(scf_control)
         END IF
         ALLOCATE (scf_env%outer_scf%variables(nvariables, nhistory))
         ALLOCATE (scf_env%outer_scf%count(nhistory))
         scf_env%outer_scf%count = 0
         ALLOCATE (scf_env%outer_scf%gradient(nvariables, nhistory))
         ALLOCATE (scf_env%outer_scf%energy(nhistory))
      ENDIF

   END SUBROUTINE qs_scf_ensure_outer_loop_vars

! **************************************************************************************************
!> \brief performs allocation of CDFT SCF variables
!> \param qs_env the qs_env where to perform the allocation
!> \param scf_env the currently active scf_env
!> \param dft_control the dft_control that holds the cdft_control type
!> \param scf_control the currently active scf_control
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_cdft_loop_vars(qs_env, scf_env, dft_control, scf_control)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(scf_control_type), POINTER                    :: scf_control

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_cdft_loop_vars', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nhistory, nvariables
      LOGICAL                                            :: do_kpoints
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gradient_history, outer_scf_history, &
                                                            variable_history

      NULLIFY (outer_scf_history, gradient_history, variable_history)
      CALL get_qs_env(qs_env=qs_env, do_kpoints=do_kpoints)
      ! Test kpoints
      IF (do_kpoints) &
         CPABORT("CDFT calculation not possible with kpoints")
      ! Check that OUTER_SCF section in DFT&SCF is active
      ! This section must always be active to facilitate
      ! switching of the CDFT and SCF control parameters in outer_loop_switch
      IF (.NOT. scf_control%outer_scf%have_scf) &
         CPABORT("Section SCF&OUTER_SCF must be active for CDFT calculations.")
      ! Initialize CDFT and outer_loop variables (constraint settings active in scf_control)
      IF (dft_control%qs_control%cdft_control%constraint_control%have_scf) THEN
         nhistory = dft_control%qs_control%cdft_control%constraint_control%max_scf + 1
         IF (scf_control%outer_scf%type /= outer_scf_none) THEN
            nvariables = outer_loop_variables_count(scf_control, &
                                                    dft_control%qs_control%cdft_control)
         ELSE
            ! First iteration: scf_control has not yet been updated
            nvariables = SIZE(dft_control%qs_control%cdft_control%target)
         END IF
         ALLOCATE (dft_control%qs_control%cdft_control%constraint%variables(nvariables, nhistory))
         ALLOCATE (dft_control%qs_control%cdft_control%constraint%count(nhistory))
         dft_control%qs_control%cdft_control%constraint%count = 0
         ALLOCATE (dft_control%qs_control%cdft_control%constraint%gradient(nvariables, nhistory))
         ALLOCATE (dft_control%qs_control%cdft_control%constraint%energy(nhistory))
         CALL qs_scf_ensure_outer_loop_vars(scf_env, scf_control, nvariables)
      ENDIF
      ! Executed only on first call (OT settings active in scf_control)
      ! Save OT settings and constraint initial values in CDFT control
      ! Then switch to constraint outer_scf settings for proper initialization of history
      IF (scf_control%outer_scf%have_scf) THEN
         IF (scf_control%outer_scf%type == outer_scf_none) THEN
            dft_control%qs_control%cdft_control%ot_control%have_scf = .TRUE.
            dft_control%qs_control%cdft_control%ot_control%max_scf = scf_control%outer_scf%max_scf
            dft_control%qs_control%cdft_control%ot_control%eps_scf = scf_control%outer_scf%eps_scf
            dft_control%qs_control%cdft_control%ot_control%step_size = scf_control%outer_scf%step_size
            dft_control%qs_control%cdft_control%ot_control%type = scf_control%outer_scf%type
            dft_control%qs_control%cdft_control%ot_control%optimizer = scf_control%outer_scf%optimizer
            dft_control%qs_control%cdft_control%ot_control%diis_buffer_length = scf_control%outer_scf%diis_buffer_length
            dft_control%qs_control%cdft_control%ot_control%bisect_trust_count = scf_control%outer_scf%bisect_trust_count
            CALL cdft_opt_type_copy(dft_control%qs_control%cdft_control%ot_control%cdft_opt_control, &
                                    scf_control%outer_scf%cdft_opt_control)
            ! In case constraint and OT extrapolation orders are different, make sure to use former
            nvariables = SIZE(dft_control%qs_control%cdft_control%target)
            IF (scf_control%outer_scf%extrapolation_order /= &
                dft_control%qs_control%cdft_control%constraint_control%extrapolation_order &
                .OR. nvariables /= 1) THEN
               DEALLOCATE (qs_env%outer_scf_history)
               DEALLOCATE (qs_env%gradient_history)
               DEALLOCATE (qs_env%variable_history)
               nhistory = dft_control%qs_control%cdft_control%constraint_control%extrapolation_order
               ALLOCATE (outer_scf_history(nvariables, nhistory))
               ALLOCATE (gradient_history(nvariables, 2))
               gradient_history = 0.0_dp
               ALLOCATE (variable_history(nvariables, 2))
               variable_history = 0.0_dp
               CALL set_qs_env(qs_env, outer_scf_history=outer_scf_history, &
                               gradient_history=gradient_history, variable_history=variable_history)
            END IF
            CALL outer_loop_switch(scf_env, scf_control, dft_control%qs_control%cdft_control, ot2cdft)
         END IF
      END IF

   END SUBROUTINE qs_scf_ensure_cdft_loop_vars

! **************************************************************************************************
!> \brief performs allocation of the mixing storage
!> \param qs_env ...
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_mixing_store(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_mixing_store', &
         routineP = moduleN//':'//routineN

      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)

      IF (scf_env%mixing_method > 0) THEN
         CALL mixing_allocate(qs_env, scf_env%mixing_method, scf_env%p_mix_new, &
                              scf_env%p_delta, dft_control%nspins, &
                              scf_env%mixing_store)
      ELSE
         NULLIFY (scf_env%p_mix_new)
      END IF

   END SUBROUTINE qs_scf_ensure_mixing_store

! **************************************************************************************************
!> \brief Performs allocation of the SCF work matrices
!>        In case of kpoints we probably don't need most of these matrices,
!>        maybe we have to initialize some matrices in the fm_pool in kpoints
!> \param qs_env ...
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_work_matrices(qs_env, scf_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_ensure_work_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, is, nao, nrow_block, nw
      LOGICAL                                            :: do_kpoints
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct, ao_mo_fmstruct
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (ao_mo_fm_pools, ao_mo_fmstruct, ao_ao_fmstruct, dft_control, matrix_s, mos)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      matrix_s_kp=matrix_s, &
                      mos=mos, &
                      scf_control=scf_control, &
                      do_kpoints=do_kpoints)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)

      ! create an ao_ao parallel matrix structure
      ao_mo_fmstruct => fm_pool_get_el_struct(ao_mo_fm_pools(1)%pool)
      CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block)
      CALL get_mo_set(mos(1)%mo_set, nao=nao)
      CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, &
                               nrow_block=nrow_block, &
                               ncol_block=nrow_block, &
                               nrow_global=nao, &
                               ncol_global=nao, &
                               template_fmstruct=ao_mo_fmstruct)

      IF ((scf_env%method /= ot_method_nr) .AND. &
          (scf_env%method /= block_davidson_diag_method_nr)) THEN
         IF (.NOT. ASSOCIATED(scf_env%scf_work1)) THEN
            nw = dft_control%nspins
            IF (do_kpoints) nw = 4
            ALLOCATE (scf_env%scf_work1(nw))
            DO is = 1, SIZE(scf_env%scf_work1)
               NULLIFY (scf_env%scf_work1(is)%matrix)
               CALL cp_fm_create(scf_env%scf_work1(is)%matrix, &
                                 matrix_struct=ao_ao_fmstruct, &
                                 name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))// &
                                 "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(is))))
            END DO
         END IF
         IF ((.NOT. ASSOCIATED(scf_env%ortho)) .AND. &
             (scf_env%method /= ot_diag_method_nr) .AND. &
             (scf_env%method /= special_diag_method_nr)) THEN
            ! Initialize fm matrix to store the Cholesky decomposition
            CALL cp_fm_create(scf_env%ortho, &
                              matrix_struct=ao_ao_fmstruct, &
                              name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))// &
                              "ORTHO_MATRIX")
            ! Initialize dbcsr matrix to store the Cholesky decomposition
            IF (scf_env%cholesky_method == cholesky_dbcsr) THEN
               ref_matrix => matrix_s(1, 1)%matrix
               CALL dbcsr_init_p(scf_env%ortho_dbcsr)
               CALL dbcsr_create(scf_env%ortho_dbcsr, template=ref_matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_init_p(scf_env%buf1_dbcsr)
               CALL dbcsr_create(scf_env%buf1_dbcsr, template=ref_matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_init_p(scf_env%buf2_dbcsr)
               CALL dbcsr_create(scf_env%buf2_dbcsr, template=ref_matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
            ELSE IF (scf_env%cholesky_method == cholesky_inverse .OR. &
                     (scf_control%level_shift /= 0.0_dp .AND. &
                      scf_env%cholesky_method == cholesky_off)) THEN
               CALL cp_fm_create(scf_env%ortho_m1, &
                                 matrix_struct=ao_ao_fmstruct, &
                                 name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))// &
                                 "ORTHO_MATRIX-1")
            END IF
         END IF
         IF (.NOT. ASSOCIATED(scf_env%scf_work2)) THEN
            CALL cp_fm_create(scf_env%scf_work2, &
                              matrix_struct=ao_ao_fmstruct, &
                              name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))// &
                              "WORK_MATRIX-2")
         END IF
      END IF

      IF (dft_control%dft_plus_u) THEN
         IF (dft_control%plus_u_method_id == plus_u_lowdin) THEN
            IF (.NOT. ASSOCIATED(scf_env%scf_work2)) THEN
               CALL cp_fm_create(scf_env%scf_work2, &
                                 matrix_struct=ao_ao_fmstruct, &
                                 name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))// &
                                 "WORK_MATRIX-2")
            END IF
            IF (.NOT. ASSOCIATED(scf_env%s_half)) THEN
               CALL cp_fm_create(scf_env%s_half, &
                                 matrix_struct=ao_ao_fmstruct, &
                                 name="S**(1/2) MATRIX")
            END IF
         END IF
      END IF

      IF (do_kpoints) THEN
         IF (.NOT. ASSOCIATED(scf_env%scf_work1)) THEN
            nw = 4
            ALLOCATE (scf_env%scf_work1(nw))
            DO is = 1, SIZE(scf_env%scf_work1)
               NULLIFY (scf_env%scf_work1(is)%matrix)
               CALL cp_fm_create(scf_env%scf_work1(is)%matrix, &
                                 matrix_struct=ao_ao_fmstruct, &
                                 name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))// &
                                 "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(is))))
            END DO
         END IF
      END IF

      CALL cp_fm_struct_release(ao_ao_fmstruct)

      CALL timestop(handle)

   END SUBROUTINE qs_scf_ensure_work_matrices

! **************************************************************************************************
!> \brief performs allocation of the MO matrices
!> \param qs_env ...
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_mos(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_mos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ic, ik, ikk, ispin, nmo, nmo_mat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_derivs_aux_fit
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools, ao_mo_fm_pools_aux_fit
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_aux_fit, mo_coeff_last
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mos_aux_fit, mos_last_converged
      TYPE(mo_set_p_type), DIMENSION(:, :), POINTER      :: mos_k
      TYPE(xas_environment_type), POINTER                :: xas_env

      CALL timeset(routineN, handle)

      NULLIFY (ao_mo_fm_pools, dft_control, mos, ao_mo_fm_pools_aux_fit)
      NULLIFY (mo_derivs, mo_coeff, mo_coeff_aux_fit, mo_derivs_aux_fit, xas_env, matrix_s)
      NULLIFY (mos_last_converged, mo_coeff_last)

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)
      MARK_USED(scf_env)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      mos=mos, &
                      matrix_s_kp=matrix_s, &
                      mos_aux_fit=mos_aux_fit, &
                      xas_env=xas_env)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
      IF (dft_control%switch_surf_dip) THEN
         CALL get_qs_env(qs_env, mos_last_converged=mos_last_converged)
      ENDIF

      nmo_mat = dft_control%nspins
      IF (dft_control%restricted) nmo_mat = 1 ! right now, there might be more mos than needed derivs

!   *** finish initialization of the MOs ***
      CPASSERT(ASSOCIATED(mos))
      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
         IF (.NOT. ASSOCIATED(mo_coeff)) THEN
            CALL init_mo_set(mos(ispin)%mo_set, &
                             fm_pool=ao_mo_fm_pools(ispin)%pool, &
                             name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))// &
                             "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
         END IF
         IF (.NOT. ASSOCIATED(mo_coeff_b)) THEN
            CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, ncol_global=nmo)
            CALL dbcsr_init_p(mos(ispin)%mo_set%mo_coeff_b)
            CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_set%mo_coeff_b, template=matrix_s(1, 1)%matrix, n=nmo, &
                                                   sym=dbcsr_type_no_symmetry)
         ENDIF
      END DO
!   *** get the mo_derivs OK if needed ***
      IF (qs_env%requires_mo_derivs) THEN
         CALL get_qs_env(qs_env, mo_derivs=mo_derivs)
         IF (.NOT. ASSOCIATED(mo_derivs)) THEN
            ALLOCATE (mo_derivs(nmo_mat))
            DO ispin = 1, nmo_mat
               CALL get_mo_set(mos(ispin)%mo_set, mo_coeff_b=mo_coeff_b)
               NULLIFY (mo_derivs(ispin)%matrix)
               CALL dbcsr_init_p(mo_derivs(ispin)%matrix)
               CALL dbcsr_create(mo_derivs(ispin)%matrix, template=mo_coeff_b, &
                                 name="mo_derivs", matrix_type=dbcsr_type_no_symmetry, &
                                 nze=0, data_type=dbcsr_type_real_default)
            ENDDO
            CALL set_qs_env(qs_env, mo_derivs=mo_derivs)
         ENDIF

      ELSE
         ! nothing should be done
      ENDIF

!   *** finish initialization of the MOs for ADMM and derivs if needed ***
      IF (dft_control%do_admm) THEN
         IF (dft_control%restricted) &
            CPABORT("ROKS with ADMM is not implemented")
         CPASSERT(ASSOCIATED(mos_aux_fit))
         CALL mpools_get(qs_env%mpools_aux_fit, ao_mo_fm_pools=ao_mo_fm_pools_aux_fit)
         DO ispin = 1, SIZE(mos_aux_fit)
            CALL get_mo_set(mos_aux_fit(ispin)%mo_set, mo_coeff=mo_coeff_aux_fit, mo_coeff_b=mo_coeff_b)
            IF (.NOT. ASSOCIATED(mo_coeff_aux_fit)) THEN
               CALL init_mo_set(mos_aux_fit(ispin)%mo_set, &
                                fm_pool=ao_mo_fm_pools_aux_fit(ispin)%pool, &
                                name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))// &
                                "%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))))
            END IF
            IF (.NOT. ASSOCIATED(mo_coeff_b)) THEN
               CALL cp_fm_get_info(mos_aux_fit(ispin)%mo_set%mo_coeff, ncol_global=nmo)
               CALL dbcsr_init_p(mos_aux_fit(ispin)%mo_set%mo_coeff_b)
               CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=matrix_s_aux_fit)
               CALL cp_dbcsr_m_by_n_from_row_template(mos_aux_fit(ispin)%mo_set%mo_coeff_b, &
                                                      template=matrix_s_aux_fit(1)%matrix, n=nmo, &
                                                      sym=dbcsr_type_no_symmetry)
            ENDIF
         END DO
         IF (qs_env%requires_mo_derivs) THEN
            CALL get_qs_env(qs_env, mo_derivs_aux_fit=mo_derivs_aux_fit)
            IF (.NOT. ASSOCIATED(mo_derivs_aux_fit)) THEN
               ALLOCATE (mo_derivs_aux_fit(nmo_mat))
               DO ispin = 1, nmo_mat
                  CALL get_mo_set(mos_aux_fit(ispin)%mo_set, mo_coeff=mo_coeff_aux_fit)
                  CALL cp_fm_create(mo_derivs_aux_fit(ispin)%matrix, mo_coeff_aux_fit%matrix_struct)
               ENDDO
               CALL set_qs_env(qs_env, mo_derivs_aux_fit=mo_derivs_aux_fit)
            ENDIF
         ELSE
            ! nothing should be done
         ENDIF
      END IF

! *** finish initialization of mos_last_converged *** [SGh]
      IF (dft_control%switch_surf_dip) THEN
         CPASSERT(ASSOCIATED(mos_last_converged))
         DO ispin = 1, SIZE(mos_last_converged)
            CALL get_mo_set(mos_last_converged(ispin)%mo_set, mo_coeff=mo_coeff_last)
            IF (.NOT. ASSOCIATED(mo_coeff_last)) THEN
               CALL init_mo_set(mos_last_converged(ispin)%mo_set, &
                                fm_ref=mos(ispin)%mo_set%mo_coeff, &
                                name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))// &
                                "%mos_last_converged"//TRIM(ADJUSTL(cp_to_string(ispin))))
            END IF
         END DO
      END IF
      ! kpoints: we have to initialize all the k-point MOs
      CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
      IF (kpoints%nkp /= 0) THEN
         ! check for some incompatible options
         IF (qs_env%requires_mo_derivs) THEN
            CPWARN("MO derivative methods flag has been switched off for kpoint calculation")
            ! we switch it off to make band structure calculations
            ! possible for OT gamma point calculations
            qs_env%requires_mo_derivs = .FALSE.
         END IF
         IF (dft_control%do_admm) &
            CPABORT("No ADMM implemented with kpoints")
         IF (dft_control%do_xas_calculation) &
            CPABORT("No XAS implemented with kpoints")
         DO ik = 1, SIZE(kpoints%kp_env)
            CALL mpools_get(kpoints%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
            mos_k => kpoints%kp_env(ik)%kpoint_env%mos
            ikk = kpoints%kp_range(1) + ik - 1
            CPASSERT(ASSOCIATED(mos_k))
            DO ispin = 1, SIZE(mos_k, 2)
               DO ic = 1, SIZE(mos_k, 1)
                  CALL get_mo_set(mos_k(ic, ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
                  IF (.NOT. ASSOCIATED(mo_coeff)) THEN
                     CALL init_mo_set(mos_k(ic, ispin)%mo_set, &
                                      fm_pool=ao_mo_fm_pools(ispin)%pool, &
                                      name="kpoints_"//TRIM(ADJUSTL(cp_to_string(ikk)))// &
                                      "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
                  END IF
                  ! no sparse matrix representation of kpoint MO vectors
                  CPASSERT(.NOT. ASSOCIATED(mo_coeff_b))
               END DO
            END DO
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_ensure_mos

! **************************************************************************************************
!> \brief sets flag for mixing/DIIS during scf
!> \param scf_control ...
!> \param scf_section ...
!> \param scf_env ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_mixing(scf_control, scf_section, scf_env, dft_control)
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_ensure_mixing', &
         routineP = moduleN//':'//routineN

      TYPE(section_vals_type), POINTER                   :: mixing_section

      SELECT CASE (scf_control%mixing_method)
      CASE (no_mix)
         scf_env%mixing_method = no_mixing_nr
         scf_env%p_mix_alpha = 1.0_dp
      CASE (direct_p_mix, kerker_mix, pulay_mix, broy_mix, broy_mix_new, multisec_mix)
         scf_env%mixing_method = scf_control%mixing_method
         mixing_section => section_vals_get_subs_vals(scf_section, "MIXING")
         IF (.NOT. ASSOCIATED(scf_env%mixing_store)) &
            CALL mixing_storage_create(scf_env%mixing_store, mixing_section, scf_env%mixing_method, &
                                       dft_control%qs_control%cutoff)
      CASE DEFAULT
         CPABORT("Unknown mixing method")
      END SELECT

      ! Disable DIIS for OT and g-space density mixing methods
      IF (scf_env%method == ot_method_nr) THEN
         ! No mixing is used with OT
         scf_env%mixing_method = no_mixing_nr
         scf_env%p_mix_alpha = 1.0_dp
         scf_env%skip_diis = .TRUE.
      END IF

      IF (scf_control%use_diag .AND. scf_env%mixing_method == no_mixing_nr) THEN
         CPABORT("Diagonalization procedures without mixing are not recommendable")
      END IF

      IF (scf_env%mixing_method > direct_mixing_nr) THEN
         scf_env%skip_diis = .TRUE.
         scf_env%p_mix_alpha = scf_env%mixing_store%alpha
         IF (scf_env%mixing_store%beta == 0.0_dp) THEN
            CPABORT("Mixing employing the Kerker damping factor needs BETA /= 0.0")
         END IF
      END IF

      IF (scf_env%mixing_method == direct_mixing_nr) THEN
         scf_env%p_mix_alpha = scf_env%mixing_store%alpha
         IF (scf_control%eps_diis < scf_control%eps_scf) THEN
            scf_env%skip_diis = .TRUE.
            CPWARN("the DIIS scheme is disabled, since EPS_DIIS < EPS_SCF")
         END IF
      END IF

   END SUBROUTINE qs_scf_ensure_mixing

! **************************************************************************************************
!> \brief sets flags for diagonalization and ensure that everything is
!>        allocated
!> \param scf_env ...
!> \param scf_section ...
!> \param qs_env ...
!> \param scf_control ...
!> \param has_unit_metric ...
! **************************************************************************************************
   SUBROUTINE qs_scf_ensure_diagonalization(scf_env, scf_section, qs_env, &
                                            scf_control, has_unit_metric)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      LOGICAL                                            :: has_unit_metric

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_ensure_diagonalization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nao, nmo
      LOGICAL                                            :: do_kpoints, need_coeff_b, not_se_or_tb
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL get_qs_env(qs_env=qs_env, do_kpoints=do_kpoints, dft_control=dft_control, mos=mos)
      not_se_or_tb = .NOT. (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb .OR. &
                            dft_control%qs_control%semi_empirical)
      need_coeff_b = .FALSE.
      scf_env%needs_ortho = .FALSE.

      IF (scf_control%use_diag) THEN
         ! sanity check whether combinations are allowed
         IF (dft_control%restricted) &
            CPABORT("OT only for restricted (ROKS)")
         SELECT CASE (scf_control%diagonalization%method)
         CASE (diag_ot, diag_block_krylov, diag_block_davidson)
            IF (.NOT. not_se_or_tb) &
               CPABORT("TB and SE not possible with OT diagonalization")
         END SELECT
         SELECT CASE (scf_control%diagonalization%method)
            ! Diagonalization: additional check whether we are in an orthonormal basis
         CASE (diag_standard)
            scf_env%method = general_diag_method_nr
            scf_env%needs_ortho = (.NOT. has_unit_metric) .AND. (.NOT. do_kpoints)
            IF (has_unit_metric) THEN
               scf_env%method = special_diag_method_nr
            END IF
            ! OT Diagonalization: not possible with ROKS
         CASE (diag_ot)
            IF (dft_control%roks) &
               CPABORT("ROKS with OT diagonalization not possible")
            IF (do_kpoints) &
               CPABORT("OT diagonalization not possible with kpoint calculations")
            scf_env%method = ot_diag_method_nr
            need_coeff_b = .TRUE.
            ! Block Krylov diagonlization: not possible with ROKS,
            ! allocation of additional matrices is needed
         CASE (diag_block_krylov)
            IF (dft_control%roks) &
               CPABORT("ROKS with block PF diagonalization not possible")
            IF (do_kpoints) &
               CPABORT("Block Krylov diagonalization not possible with kpoint calculations")
            scf_env%method = block_krylov_diag_method_nr
            scf_env%needs_ortho = .TRUE.
            IF (.NOT. ASSOCIATED(scf_env%krylov_space)) &
               CALL krylov_space_create(scf_env%krylov_space, scf_section)
            CALL krylov_space_allocate(scf_env%krylov_space, scf_control, mos)
            ! Block davidson diagonlization: allocation of additional matrices is needed
         CASE (diag_block_davidson)
            IF (do_kpoints) &
               CPABORT("Block Davidson diagonalization not possible with kpoint calculations")
            scf_env%method = block_davidson_diag_method_nr
            IF (.NOT. ASSOCIATED(scf_env%block_davidson_env)) &
               CALL block_davidson_env_create(scf_env%block_davidson_env, dft_control%nspins, &
                                              scf_section)
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
               CALL block_davidson_allocate(scf_env%block_davidson_env(ispin), mo_coeff, nao, nmo)
            END DO
            need_coeff_b = .TRUE.
            ! Filter matrix diagonalisation method
         CASE (diag_filter_matrix)
            scf_env%method = filter_matrix_diag_method_nr
            IF (.NOT. fb_env_has_data(scf_env%filter_matrix_env)) THEN
               CALL fb_env_create(scf_env%filter_matrix_env)
            END IF
            CALL fb_env_read_input(scf_env%filter_matrix_env, scf_section)
            CALL fb_env_build_rcut_auto(scf_env%filter_matrix_env, qs_env)
            CALL fb_env_write_info(scf_env%filter_matrix_env, qs_env, scf_section)
            CALL fb_distribution_build(scf_env%filter_matrix_env, qs_env, scf_section)
            CALL fb_env_build_atomic_halos(scf_env%filter_matrix_env, qs_env, scf_section)
         CASE DEFAULT
            CPABORT("Unknown diagonalization method")
         END SELECT
         ! Check if subspace diagonlization is requested: allocation of additional matrices is needed
         IF (scf_control%do_diag_sub) THEN
            scf_env%needs_ortho = .TRUE.
            IF (.NOT. ASSOCIATED(scf_env%subspace_env)) &
               CALL diag_subspace_env_create(scf_env%subspace_env, scf_section, &
                                             dft_control%qs_control%cutoff)
            CALL diag_subspace_allocate(scf_env%subspace_env, qs_env, mos)
            IF (do_kpoints) &
               CPABORT("No subspace diagonlization with kpoint calculation")
         END IF
         ! OT: check if OT is used instead of diagonlization. Not possible with added MOS at the moment
      ELSEIF (scf_control%use_ot) THEN
         scf_env%method = ot_method_nr
         need_coeff_b = .TRUE.
         IF (scf_control%added_mos(1) + scf_control%added_mos(2) > 0) &
            CPABORT("OT with ADDED_MOS/=0 not implemented")
         IF (dft_control%restricted .AND. dft_control%nspins .NE. 2) &
            CPABORT("nspin must be 2 for restricted (ROKS)")
         IF (do_kpoints) &
            CPABORT("OT not possible with kpoint calculations")
      ELSE
         CPABORT("OT or DIAGONALIZATION have to be set")
      END IF
      DO ispin = 1, dft_control%nspins
         mos(ispin)%mo_set%use_mo_coeff_b = need_coeff_b
      END DO

   END SUBROUTINE qs_scf_ensure_diagonalization

! **************************************************************************************************
!> \brief performs those initialisations that need to be done only once
!>       (e.g. that only depend on the atomic positions)
!>       this will be called in scf
!> \param scf_env ...
!> \param qs_env ...
!> \param scf_section ...
!> \param scf_control ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE init_scf_run(scf_env, qs_env, scf_section, scf_control)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(scf_control_type), POINTER                    :: scf_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_run', routineP = moduleN//':'//routineN

      INTEGER                                            :: after, handle, ikind, iw, nao, ndep, &
                                                            output_unit
      LOGICAL                                            :: dft_plus_u_atom, do_kpoints, &
                                                            init_u_ramping_each_scf, omit_headers, &
                                                            s_minus_half_available
      REAL(KIND=dp)                                      :: u_ramping
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xas_environment_type), POINTER                :: xas_env

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, matrix_s, dft_control, mos, qs_kind, rho, xas_env)

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(qs_env%ref_count > 0)
      NULLIFY (para_env)

      s_minus_half_available = .FALSE.
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      qs_kind_set=qs_kind_set, &
                      mos=mos, &
                      rho=rho, &
                      nelectron_total=scf_env%nelectron, &
                      do_kpoints=do_kpoints, &
                      para_env=para_env, &
                      xas_env=xas_env)

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      CALL qs_scf_initial_info(output_unit, mos, dft_control)
      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      ! calc ortho matrix
      ndep = 0
      IF (scf_env%needs_ortho) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, scf_env%ortho)
         IF (scf_env%cholesky_method > cholesky_off) THEN
            CALL cp_fm_cholesky_decompose(scf_env%ortho)
            IF (scf_env%cholesky_method == cholesky_dbcsr) THEN
               CALL cp_fm_triangular_invert(scf_env%ortho)
               CALL cp_fm_set_all(scf_env%scf_work2, 0.0_dp)
               CALL cp_fm_to_fm_triangular(scf_env%ortho, scf_env%scf_work2, "U")
               CALL copy_fm_to_dbcsr(scf_env%scf_work2, scf_env%ortho_dbcsr)
            ELSE IF (scf_env%cholesky_method == cholesky_inverse) THEN
               CALL cp_fm_to_fm(scf_env%ortho, scf_env%ortho_m1)
               CALL cp_fm_triangular_invert(scf_env%ortho_m1)
            END IF
         ELSE
            CALL cp_fm_power(scf_env%ortho, scf_env%scf_work2, -0.5_dp, &
                             scf_control%eps_eigval, ndep)
            IF (scf_control%level_shift /= 0.0_dp) THEN
               CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, scf_env%ortho_m1)
               CALL cp_fm_power(scf_env%ortho_m1, scf_env%scf_work2, 0.5_dp, &
                                scf_control%eps_eigval, ndep)
            END IF
            s_minus_half_available = .TRUE.
         END IF

         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              qs_env%input, "DFT%PRINT%AO_MATRICES/ORTHO"), cp_p_file)) THEN
            iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/ORTHO", &
                                      extension=".Log")
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
            after = MIN(MAX(after, 1), 16)
            CALL write_fm_with_basis_info(scf_env%ortho, 4, after, qs_env, &
                                          para_env, output_unit=iw, omit_headers=omit_headers)
            CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                              "DFT%PRINT%AO_MATRICES/ORTHO")
         END IF
      END IF

      CALL get_mo_set(mo_set=mos(1)%mo_set, nao=nao)

      ! DFT+U methods based on Lowdin charges need S^(1/2)
      IF (dft_control%dft_plus_u) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         IF (dft_control%plus_u_method_id == plus_u_lowdin) THEN
            IF (s_minus_half_available) THEN
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, scf_env%ortho, scf_env%s_half, &
                                            nao)
            ELSE
               CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, scf_env%s_half)
               CALL cp_fm_power(scf_env%s_half, scf_env%scf_work2, 0.5_dp, &
                                scf_control%eps_eigval, ndep)
            END IF
         END IF
         DO ikind = 1, SIZE(qs_kind_set)
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, &
                             dft_plus_u_atom=dft_plus_u_atom, &
                             u_ramping=u_ramping, &
                             init_u_ramping_each_scf=init_u_ramping_each_scf)
            IF (dft_plus_u_atom .AND. (u_ramping /= 0.0_dp)) THEN
               IF (init_u_ramping_each_scf) THEN
                  CALL set_qs_kind(qs_kind=qs_kind, u_minus_j=0.0_dp)
               END IF
            END IF
         END DO
      END IF

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(T2,A,T71,I10)") &
            "Number of independent orbital functions:", nao - ndep
      END IF
      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      ! extrapolate outer loop variables
      IF (scf_control%outer_scf%have_scf) THEN
         CALL outer_loop_extrapolate(qs_env)
      ENDIF

      ! initializes rho and the mos
      IF (ASSOCIATED(qs_env%xas_env)) THEN
         ! if just optimized wfn, e.g. ground state
         ! changes come from a perturbation, e.g., the occupation numbers
         ! it could be generalized for other cases, at the moment used only for core level spectroscopy
         ! initialize the density with the localized mos
         CALL xas_initialize_rho(qs_env, scf_env, scf_control)
      ELSE
         CALL scf_env_initial_rho_setup(scf_env, qs_env=qs_env, &
                                        scf_section=scf_section, scf_control=scf_control)
      END IF

      ! Frozen density approximation
      IF (ASSOCIATED(qs_env%wf_history)) THEN
         IF (qs_env%wf_history%interpolation_method_nr == wfi_frozen_method_nr) THEN
            IF (.NOT. ASSOCIATED(qs_env%wf_history%past_states(1)%snapshot)) THEN
               CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
               CALL qs_rho_create(qs_env%wf_history%past_states(1)%snapshot%rho_frozen)
               CALL duplicate_rho_type(rho_input=rho, &
                                       rho_output=qs_env%wf_history%past_states(1)%snapshot%rho_frozen, &
                                       qs_env=qs_env)
            END IF
         END IF
      END IF

      !image charge method, calculate image_matrix if required
      IF (qs_env%qmmm) THEN
         IF (qs_env%qmmm .AND. qs_env%qmmm_env_qm%image_charge) THEN
            CALL conditional_calc_image_matrix(qs_env=qs_env, &
                                               qmmm_env=qs_env%qmmm_env_qm)
         ENDIF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE init_scf_run

! **************************************************************************************************
!> \brief Initializes rho and the mos, so that an scf cycle can start
!> \param scf_env the scf env in which to do the scf
!> \param qs_env the qs env the scf_env lives in
!> \param scf_section ...
!> \param scf_control ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE scf_env_initial_rho_setup(scf_env, qs_env, scf_section, scf_control)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: scf_section
      TYPE(scf_control_type), POINTER                    :: scf_control

      CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_initial_rho_setup', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: extrapolation_method_nr, handle, ispin, &
                                                            nmo, output_unit
      LOGICAL                                            :: orthogonal_wf
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom

      CALL timeset(routineN, handle)
      NULLIFY (mo_coeff, rho, dft_control, para_env, mos)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env, &
                      rho=rho, &
                      mos=mos, &
                      dft_control=dft_control, &
                      para_env=para_env)

      extrapolation_method_nr = wfi_use_guess_method_nr
      IF (ASSOCIATED(qs_env%wf_history)) THEN
         CALL wfi_extrapolate(qs_env%wf_history, &
                              qs_env=qs_env, dt=1.0_dp, &
                              extrapolation_method_nr=extrapolation_method_nr, &
                              orthogonal_wf=orthogonal_wf)
         ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
         IF ((.NOT. orthogonal_wf) .AND. &
             (scf_env%method == ot_method_nr) .AND. &
             (.NOT. (extrapolation_method_nr == wfi_use_guess_method_nr))) THEN
            DO ispin = 1, SIZE(mos)
               CALL get_mo_set(mos(ispin)%mo_set, &
                               mo_coeff=mo_coeff, nmo=nmo)
               CALL reorthogonalize_vectors(qs_env, &
                                            v_matrix=mo_coeff, n_col=nmo)
               CALL set_mo_occupation(mo_set=mos(ispin)%mo_set, &
                                      smear=scf_control%smear)
            END DO
         END IF
      END IF
      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A,I0)") &
            "Extrapolation method: "// &
            TRIM(wfi_get_method_label(extrapolation_method_nr))
         IF (extrapolation_method_nr == wfi_ps_method_nr) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,I0,A)") &
               "Extrapolation order:  ", &
               MAX((MIN(qs_env%wf_history%memory_depth, qs_env%wf_history%snapshot_count) - 1), 0)
         END IF
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%PROGRAM_RUN_INFO")
      IF (extrapolation_method_nr == wfi_use_guess_method_nr) THEN
         CALL calculate_first_density_matrix(scf_env=scf_env, qs_env=qs_env)
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
      END IF

      ! Some preparation for the mixing
      IF (scf_env%mixing_method > 1) THEN
         IF (dft_control%qs_control%gapw) THEN
            CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom)
            CALL mixing_init(scf_env%mixing_method, rho, scf_env%mixing_store, &
                             para_env, rho_atom=rho_atom)
         ELSEIF (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb) THEN
            CALL charge_mixing_init(scf_env%mixing_store)
         ELSEIF (dft_control%qs_control%semi_empirical) THEN
            CPABORT('SE Code not possible')
         ELSE
            CALL mixing_init(scf_env%mixing_method, rho, scf_env%mixing_store, &
                             para_env)
         END IF
      END IF

      DO ispin = 1, SIZE(mos) !fm->dbcsr
         IF (mos(ispin)%mo_set%use_mo_coeff_b) THEN
            CALL copy_fm_to_dbcsr(mos(ispin)%mo_set%mo_coeff, &
                                  mos(ispin)%mo_set%mo_coeff_b) !fm->dbcsr
         ENDIF
      ENDDO !fm->dbcsr

      CALL timestop(handle)

   END SUBROUTINE scf_env_initial_rho_setup

END MODULE qs_scf_initialization
