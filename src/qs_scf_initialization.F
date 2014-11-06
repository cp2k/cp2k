!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utility routines for qs_scf
! *****************************************************************************
MODULE qs_scf_initialization
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       array_i1d_obj, array_release, cp_create_bl_distribution, &
       cp_dbcsr_col_block_sizes, cp_dbcsr_create, cp_dbcsr_distribution, &
       cp_dbcsr_distribution_release, cp_dbcsr_init_p, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_type, dbcsr_distribution_mp, &
       dbcsr_distribution_new, dbcsr_distribution_obj, &
       dbcsr_distribution_row_dist, dbcsr_mp_npcols, dbcsr_type_no_symmetry, &
       dbcsr_type_real_default
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
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
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       broy_mix, broy_mix_new, cholesky_dbcsr, cholesky_inverse, &
       cholesky_off, densities_guess, diag_block_davidson, diag_block_krylov, &
       diag_ot, diag_standard, direct_p_mix, kerker_mix, multisec_mix, &
       no_mix, plus_u_lowdin, pulay_mix, wfi_frozen_method_nr, &
       wfi_use_guess_method_nr
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE kpoint_types,                    ONLY: kpoint_type
  USE lri_optimize_ri_basis,           ONLY: optimize_lri_basis
  USE qmmm_image_charge,               ONLY: calculate_image_matrix
  USE qs_block_davidson_types,         ONLY: block_davidson_allocate,&
                                             block_davidson_env_create
  USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                             mixing_storage_create,&
                                             mixing_storage_release,&
                                             mixing_storage_type,&
                                             no_mixing_nr
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_gspace_mixing,                ONLY: mixing_allocate,&
                                             mixing_init
  USE qs_initial_guess,                ONLY: calculate_first_density_matrix
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type,&
                                             set_qs_kind
  USE qs_ks_methods,                   ONLY: qs_ks_did_change
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation,&
                                             set_mo_set
  USE qs_outer_scf,                    ONLY: outer_loop_extrapolate,&
                                             outer_loop_variables_count
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                             qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_create,&
                                             qs_rho_type
  USE qs_scf_diagonalization,          ONLY: diag_subspace_allocate
  USE qs_scf_lanczos,                  ONLY: krylov_space_allocate
  USE qs_scf_output,                   ONLY: qs_scf_initial_info
  USE qs_scf_post_dftb,                ONLY: scf_post_calculation_dftb
  USE qs_scf_post_gpw,                 ONLY: scf_post_calculation_gpw
  USE qs_scf_post_scptb,               ONLY: scf_post_calculation_scptb
  USE qs_scf_post_se,                  ONLY: scf_post_calculation_se
  USE qs_scf_types,                    ONLY: &
       block_davidson_diag_method_nr, block_krylov_diag_method_nr, &
       diag_subspace_env_create, general_diag_method_nr, krylov_space_create, &
       ot_diag_method_nr, ot_method_nr, qs_scf_env_type, scf_env_create, &
       scf_env_release, special_diag_method_nr
  USE qs_wf_history_methods,           ONLY: reorthogonalize_vectors,&
                                             wfi_extrapolate,&
                                             wfi_get_method_label,&
                                             wfi_update
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_env_types,                   ONLY: xas_environment_type
  USE xas_restart,                     ONLY: xas_initialize_rho
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_initialization'

  PUBLIC:: qs_scf_env_initialize, qs_scf_compute_properties

CONTAINS

! *****************************************************************************
!> \brief initializes input parameters if needed or restores values from 
!>        previous runs to fill scf_env with the values required for scf
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param scf_env ...
!> \param do_xas_tp ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************

  SUBROUTINE qs_scf_env_initialize(qs_env,scf_env,do_xas_tp,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    LOGICAL, INTENT(IN), OPTIONAL            :: do_xas_tp
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_env_initialize', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: do_my_xas_tp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                scf_section

    do_my_xas_tp=.FALSE.
    IF(PRESENT(do_xas_tp))do_my_xas_tp=do_xas_tp

    CALL get_qs_env(qs_env,error=error,input=input, rho=rho,&
                       dft_control=dft_control,scf_control=scf_control)
    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
    scf_section => section_vals_get_subs_vals(dft_section,"SCF",error=error)

    CALL qs_scf_ensure_scf_env(qs_env,scf_env,error)

    CALL section_vals_val_get(scf_section,"CHOLESKY",i_val=scf_env%cholesky_method,error=error)

    CALL qs_scf_ensure_mos(qs_env,scf_env,do_my_xas_tp,error)

    ! set flags for diagonalization
    CALL qs_scf_ensure_diagonalization(scf_env,scf_section,qs_env,&
                                     scf_control,qs_env%has_unit_metric,error)
    ! set parameters for mixing/DIIS during scf
    CALL qs_scf_ensure_mixing(scf_control,scf_section,scf_env,dft_control,error)

    CALL qs_scf_ensure_work_matrices(qs_env,scf_env,error)    

    CALL qs_scf_ensure_mixing_store(qs_env,do_my_xas_tp,error)

    CALL qs_scf_ensure_outer_loop_vars(scf_env,scf_control,dft_control,error)

    CALL init_scf_run(scf_env, qs_env, scf_section, do_my_xas_tp, error)

  END SUBROUTINE qs_scf_env_initialize

! *****************************************************************************
!> \brief makes sure scf_env is allocated (might already be from before)
!>        in case it is present the g-space mixing storage is reset
!> \param qs_env ...
!> \param scf_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scf_ensure_scf_env(qs_env,scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_scf_env', &
      routineP = moduleN//':'//routineN

    TYPE(qs_rho_type), POINTER               :: rho

    IF (.NOT.ASSOCIATED(scf_env)) THEN ! i.e. for MD this is associated on the second step (it so seems)
       CALL scf_env_create(scf_env, error=error)
       CALL set_qs_env(qs_env,scf_env=scf_env,error=error)
       CALL scf_env_release(scf_env,error=error)
       CALL get_qs_env(qs_env=qs_env,scf_env=scf_env,error=error)
    ELSE
       ! Reallocate mixing store, if the g space grid (cell) has changed
       SELECT CASE (scf_env%mixing_method)
       CASE (kerker_mix,pulay_mix,broy_mix,broy_mix_new,multisec_mix)
         IF (ASSOCIATED(scf_env%mixing_store)) THEN
           ! The current mixing_store data structure does not allow for an unique
           ! grid comparison, but the probability that the 1d lengths of the old and
           ! the new grid are accidentily equal is rather low
           CALL get_qs_env(qs_env,error=error, rho=rho)
           IF (SIZE(rho%rho_g(1)%pw%pw_grid%gsq) /= SIZE(scf_env%mixing_store%rhoin(1)%cc)) &
             CALL mixing_storage_release(scf_env%mixing_store,error=error)
         END IF
       END SELECT
    END IF

   END SUBROUTINE qs_scf_ensure_scf_env
! *****************************************************************************
!> \brief performs allocation of outer SCF variables
!> \param scf_env ...
!> \param scf_control ...
!> \param dft_control ...
!> \param error ...
! *****************************************************************************
   SUBROUTINE qs_scf_ensure_outer_loop_vars(scf_env,scf_control,dft_control,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'qs_scf_ensure_outer_loop_vars', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: nhistory, nvariables, stat
    LOGICAL                                  :: failure

    IF (scf_control%outer_scf%have_scf) THEN
       nhistory = scf_control%outer_scf%max_scf + 1
       nvariables = outer_loop_variables_count(scf_control,error)
       ALLOCATE(scf_env%outer_scf%variables(nvariables,nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(scf_env%outer_scf%count(nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       scf_env%outer_scf%count=0
       ALLOCATE(scf_env%outer_scf%gradient(nvariables,nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(scf_env%outer_scf%energy(nhistory),STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
 
  END SUBROUTINE qs_scf_ensure_outer_loop_vars

! *****************************************************************************
!> \brief performs allocation of the mixing storage
!> \param qs_env ...
!> \param do_xas_tp ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scf_ensure_mixing_store(qs_env,do_xas_tp,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL                                  :: do_xas_tp
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_mixing_store', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: mixing_method
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(xas_environment_type), POINTER      :: xas_env

    NULLIFY(dft_control,scf_env,xas_env)
    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    scf_env=scf_env,xas_env=xas_env,error=error)
 
    IF(do_xas_tp) THEN
      mixing_method = xas_env%mixing_method
      mixing_store => xas_env%mixing_store
    ELSE
      mixing_method = scf_env%mixing_method
      mixing_store => scf_env%mixing_store
    END IF
    IF (mixing_method>0) THEN
      CALL mixing_allocate(qs_env,mixing_method,mixing_store,dft_control%nspins,error)
    ELSE
       NULLIFY(scf_env%p_mix_new)
    END IF

  END SUBROUTINE qs_scf_ensure_mixing_store
! *****************************************************************************
!> \brief Performs allocation of the SCF work matrices
!>        In case of kpoints we probably don't need most of these matrices,
!>        maybe we have to initialize some matrices in the fm_pool in kpoints
!> \param qs_env ...
!> \param scf_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scf_ensure_work_matrices(qs_env,scf_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_ensure_work_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, is, nao, nrow_block, &
                                                nw, stat
    LOGICAL                                  :: do_kpoints, failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: ref_matrix
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct, ao_mo_fmstruct
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)

    NULLIFY(ao_mo_fm_pools, ao_mo_fmstruct, ao_ao_fmstruct, dft_control, matrix_s, mos)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    mos=mos,&
                    do_kpoints=do_kpoints,&
                    error=error)
    CALL mpools_get(qs_env%mpools,ao_mo_fm_pools=ao_mo_fm_pools,error=error)

    ! create an ao_ao parallel matrix structure
    ao_mo_fmstruct => fm_pool_get_el_struct(ao_mo_fm_pools(1)%pool,error=error)
    CALL cp_fm_struct_get(ao_mo_fmstruct, nrow_block=nrow_block,error=error)
    CALL get_mo_set(mos(1)%mo_set,nao=nao)
    CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
                             nrow_block=nrow_block,&
                             ncol_block=nrow_block,&
                             nrow_global=nao,&
                             ncol_global=nao,&
                             template_fmstruct=ao_mo_fmstruct,&
                             error=error)

    IF ((scf_env%method /= ot_method_nr) .AND. &
        (scf_env%method /= block_davidson_diag_method_nr) ) THEN
       IF (.NOT.ASSOCIATED(scf_env%scf_work1)) THEN
          nw = dft_control%nspins
          IF (do_kpoints) nw = 4
          ALLOCATE(scf_env%scf_work1(nw), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO is=1,SIZE(scf_env%scf_work1)
             NULLIFY(scf_env%scf_work1(is)%matrix)
             CALL cp_fm_create(scf_env%scf_work1(is)%matrix,&
                               matrix_struct=ao_ao_fmstruct,&
                               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                               "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(is))),&
                               error=error)
          END DO
       END IF
       IF ((.NOT.ASSOCIATED(scf_env%ortho)).AND.&
           (scf_env%method /= ot_diag_method_nr).AND.&
           (scf_env%method /= special_diag_method_nr)) THEN
          ! Initialize fm matrix to store the Cholesky decomposition
          CALL cp_fm_create(scf_env%ortho,&
                            matrix_struct=ao_ao_fmstruct,&
                            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                            "ORTHO_MATRIX",&
                            error=error)
          ! Initialize dbcsr matrix to store the Cholesky decomposition
          IF (scf_env%cholesky_method == cholesky_dbcsr) THEN
             ref_matrix => matrix_s(1)%matrix
             CALL cp_dbcsr_init_p(scf_env%ortho_dbcsr,error=error)
             CALL cp_dbcsr_create(scf_env%ortho_dbcsr,template=ref_matrix,&
                                  matrix_type=dbcsr_type_no_symmetry,error=error)
             CALL cp_dbcsr_init_p(scf_env%buf1_dbcsr,error=error)
             CALL cp_dbcsr_create(scf_env%buf1_dbcsr,template=ref_matrix,&
                                  matrix_type=dbcsr_type_no_symmetry,error=error)
             CALL cp_dbcsr_init_p(scf_env%buf2_dbcsr,error=error)
             CALL cp_dbcsr_create(scf_env%buf2_dbcsr,template=ref_matrix,&
                                  matrix_type=dbcsr_type_no_symmetry,error=error)
          ELSE IF (scf_env%cholesky_method == cholesky_inverse) THEN
             CALL cp_fm_create(scf_env%ortho_m1,&
                               matrix_struct=ao_ao_fmstruct,&
                               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                               "ORTHO_MATRIX-1",&
                               error=error)
          END IF
       END IF
       IF (.NOT.ASSOCIATED(scf_env%scf_work2)) THEN
          CALL cp_fm_create(scf_env%scf_work2,&
                            matrix_struct=ao_ao_fmstruct,&
                            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                            "WORK_MATRIX-2",&
                            error=error)
       END IF
    END IF

    IF (dft_control%dft_plus_u) THEN
      IF (dft_control%plus_u_method_id == plus_u_lowdin) THEN
        IF (.NOT.ASSOCIATED(scf_env%scf_work2)) THEN
          CALL cp_fm_create(scf_env%scf_work2,&
                            matrix_struct=ao_ao_fmstruct,&
                            name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                            "WORK_MATRIX-2",&
                            error=error)
        END IF
        IF (.NOT.ASSOCIATED(scf_env%s_half)) THEN
          CALL cp_fm_create(scf_env%s_half,&
                            matrix_struct=ao_ao_fmstruct,&
                            name="S**(1/2) MATRIX",&
                            error=error)
        END IF
      END IF
    END IF

    IF (do_kpoints) THEN
       IF (.NOT.ASSOCIATED(scf_env%scf_work1)) THEN
          nw = 4
          ALLOCATE(scf_env%scf_work1(nw), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO is=1,SIZE(scf_env%scf_work1)
             NULLIFY(scf_env%scf_work1(is)%matrix)
             CALL cp_fm_create(scf_env%scf_work1(is)%matrix,&
                               matrix_struct=ao_ao_fmstruct,&
                               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//&
                               "WORK_MATRIX-1-"//TRIM(ADJUSTL(cp_to_string(is))),&
                               error=error)
          END DO
       END IF
    END IF

    CALL cp_fm_struct_release(ao_ao_fmstruct,error=error)

    CALL timestop(handle)

  END SUBROUTINE qs_scf_ensure_work_matrices

! *****************************************************************************
!> \brief performs allocation of the MO matrices  
!> \param qs_env ...
!> \param scf_env ...
!> \param do_xas_tp ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scf_ensure_mos(qs_env,scf_env,do_xas_tp,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    LOGICAL, INTENT(IN)                      :: do_xas_tp
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_ensure_mos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ic, ik, ikk, ispin, &
                                                nmo, nmo_mat
    LOGICAL                                  :: failure
    TYPE(array_i1d_obj)                      :: col_blk_size, col_dist
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, mo_derivs
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs_aux_fit
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools, &
                                                ao_mo_fm_pools_aux_fit
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dbcsr_distribution_obj)             :: dist
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(mo_set_p_type), DIMENSION(:, :), &
      POINTER                                :: mos_k
    TYPE(xas_environment_type), POINTER      :: xas_env

    CALL timeset(routineN,handle)

    NULLIFY(ao_mo_fm_pools, dft_control, mos, ao_mo_fm_pools_aux_fit)
    NULLIFY( mo_derivs, mo_coeff, mo_coeff_aux_fit, mo_derivs_aux_fit, xas_env, matrix_s)

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         mos=mos,&
         matrix_s=matrix_s,&
         mos_aux_fit=mos_aux_fit,&
         xas_env=xas_env,error=error)
    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools,&
         error=error)

    nmo_mat=dft_control%nspins
    IF(dft_control%restricted)nmo_mat=1  ! right now, there might be more mos than needed derivs

!   *** finish initialization of the MOs ***
    CPPrecondition(ASSOCIATED(mos),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       DO ispin=1,SIZE(mos)
          CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,mo_coeff_b=mo_coeff_b)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(mos(ispin)%mo_set,&
                  ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
          IF(.not.ASSOCIATED(mo_coeff_b)) THEN
             CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff,ncol_global=nmo,error=error)
             CALL cp_create_bl_distribution (col_dist, col_blk_size, nmo, &! fm->dbcsr
                  dbcsr_mp_npcols(dbcsr_distribution_mp(cp_dbcsr_distribution(matrix_s(1)%matrix))))! fm->dbcsr
             CALL dbcsr_distribution_new (dist, dbcsr_distribution_mp (cp_dbcsr_distribution(matrix_s(1)%matrix)),&! fm->dbcsr
                  dbcsr_distribution_row_dist(cp_dbcsr_distribution(matrix_s(1)%matrix)), col_dist)! fm->dbcsr
             CALL cp_dbcsr_init_p(mos(ispin)%mo_set%mo_coeff_b,error=error)
             CALL cp_dbcsr_create(mos(ispin)%mo_set%mo_coeff_b, "mo_coeff_b", dist, dbcsr_type_no_symmetry,&! fm->dbcsr
                  cp_dbcsr_row_block_sizes(matrix_s(1)%matrix), col_blk_size, 0, 0, dbcsr_type_real_default,&!fm->dbcsr
                  error=error)! fm->dbcsr
             CALL cp_dbcsr_distribution_release (dist)! fm->dbcsr
             CALL array_release (col_blk_size)! fm->dbcsr
             CALL array_release (col_dist)! fm->dbcsr
          ENDIF
    ! very first tests for xas
          IF(do_xas_tp .AND. ispin==1) THEN
            CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,smear=xas_env%smear,&
                 xas_env=xas_env,error=error)
            CALL set_mo_set(mos(ispin)%mo_set,&
                 uniform_occupation=.FALSE.,error=error)
          END IF
       END DO
    END IF
!   *** get the mo_derivs OK if needed ***

    IF (qs_env%requires_mo_derivs) THEN
       CALL get_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
       IF (.NOT.ASSOCIATED(mo_derivs)) THEN
          ALLOCATE(mo_derivs(nmo_mat))
          DO ispin=1,nmo_mat
             CALL get_mo_set(mos(ispin)%mo_set,mo_coeff_b=mo_coeff_b)
             NULLIFY(mo_derivs(ispin)%matrix)
             CALL cp_dbcsr_init_p(mo_derivs(ispin)%matrix,error=error)
             CALL cp_dbcsr_create(mo_derivs(ispin)%matrix, "mo_derivs",&
                  cp_dbcsr_distribution(mo_coeff_b), dbcsr_type_no_symmetry,&
                  cp_dbcsr_row_block_sizes(mo_coeff_b), cp_dbcsr_col_block_sizes(mo_coeff_b),&
                  0, 0, dbcsr_type_real_default, error=error)
          ENDDO
          CALL set_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
       ENDIF

    ELSE
       ! nothing should be done
    ENDIF

!   *** finish initialization of the MOs for ADMM and derivs if needed ***
    IF(dft_control%do_admm) THEN
      CPPrecondition(ASSOCIATED(mos_aux_fit),cp_failure_level,routineP,error,failure)
      CALL mpools_get(qs_env%mpools_aux_fit, ao_mo_fm_pools=ao_mo_fm_pools_aux_fit,&
           error=error)
      IF (.NOT.failure) THEN
         DO ispin=1,SIZE(mos_aux_fit)
            CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
            IF (.NOT.ASSOCIATED(mo_coeff_aux_fit)) THEN
               CALL init_mo_set(mos_aux_fit(ispin)%mo_set,&
                    ao_mo_fm_pools_aux_fit(ispin)%pool,&
                    name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                    "%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                    error=error)
            END IF
            IF(do_xas_tp .AND. ispin==1) THEN
              CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
            END IF
         END DO
      END IF
      IF (qs_env%requires_mo_derivs) THEN
       CALL get_qs_env(qs_env,mo_derivs_aux_fit=mo_derivs_aux_fit,error=error)
       IF (.NOT.ASSOCIATED(mo_derivs_aux_fit)) THEN
         ALLOCATE(mo_derivs_aux_fit(nmo_mat))
         DO ispin=1,nmo_mat
            CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
            CALL cp_fm_create(mo_derivs_aux_fit(ispin)%matrix,mo_coeff_aux_fit%matrix_struct,error=error)
         ENDDO
         CALL set_qs_env(qs_env,mo_derivs_aux_fit=mo_derivs_aux_fit,error=error)
       ENDIF
      ELSE
         ! nothing should be done
      ENDIF
    END IF

    ! kpoints: we have to initialize all the k-point MOs
    CALL get_qs_env(qs_env=qs_env,kpoints=kpoints,error=error)
    IF(kpoints%nkp /= 0) THEN
       ! check for some incompatible options
       CALL cp_assert(.NOT.(qs_env%requires_mo_derivs),cp_failure_level,cp_assertion_failed,&
            routineP,"No MO derivative methods allowed with kpoints",error,failure)
       CALL cp_assert(.NOT.(dft_control%do_admm),cp_failure_level,cp_assertion_failed,&
            routineP,"No ADMM implemented with kpoints",error,failure)
       CALL cp_assert(.NOT.(do_xas_tp),cp_failure_level,cp_assertion_failed,&
            routineP,"No XAS implemented with kpoints",error,failure)
       DO ik=1,SIZE(kpoints%kp_env)
          CALL mpools_get(kpoints%mpools, ao_mo_fm_pools=ao_mo_fm_pools,error=error)
          mos_k => kpoints%kp_env(ik)%kpoint_env%mos
          ikk = kpoints%kp_range(1) + ik - 1
          CPPrecondition(ASSOCIATED(mos_k),cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO ispin=1,SIZE(mos_k,2)
                DO ic=1,SIZE(mos_k,1)
                   CALL get_mo_set(mos_k(ic,ispin)%mo_set,mo_coeff=mo_coeff,mo_coeff_b=mo_coeff_b)
                   IF (.NOT.ASSOCIATED(mo_coeff)) THEN
                      CALL init_mo_set(mos_k(ic,ispin)%mo_set,ao_mo_fm_pools(ispin)%pool,&
                           name="kpoints_"//TRIM(ADJUSTL(cp_to_string(ikk)))//&
                           "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),error=error)
                   END IF
                   ! no sparse matrix representation of kpoint MO vectors
                   CPPrecondition(.NOT.ASSOCIATED(mo_coeff_b),cp_failure_level,routineP,error,failure)
                END DO
             END DO
          END IF
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_ensure_mos

! *****************************************************************************
!> \brief sets flag for mixing/DIIS during scf 
!> \param scf_control ...
!> \param scf_section ...
!> \param scf_env ...
!> \param dft_control ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_scf_ensure_mixing(scf_control,scf_section,scf_env,dft_control,error)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_ensure_mixing', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: mixing_section

    SELECT CASE(scf_control%mixing_method)
    CASE(no_mix)
       scf_env%mixing_method=no_mixing_nr
       scf_env%p_mix_alpha = 1.0_dp
    CASE(direct_p_mix,kerker_mix,pulay_mix,broy_mix,broy_mix_new,multisec_mix)
      scf_env%mixing_method =  scf_control%mixing_method
      mixing_section => section_vals_get_subs_vals(scf_section,"MIXING",error=error)
      IF(.NOT.ASSOCIATED(scf_env%mixing_store))&
      CALL mixing_storage_create(scf_env%mixing_store, mixing_section, scf_env%mixing_method, &
          dft_control%qs_control%cutoff, error=error)
    CASE DEFAULT
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"Unknown mixing method",error,failure)
    END SELECT

    ! Disable DIIS for OT and g-space density mixing methods
    IF( scf_env%method==ot_method_nr ) THEN
      ! No mixing is used with OT
       scf_env%mixing_method=no_mixing_nr
       scf_env%p_mix_alpha = 1.0_dp
       scf_env%skip_diis = .TRUE.
    END IF

    IF (scf_control%use_diag .AND. scf_env%mixing_method==no_mixing_nr) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Diagonalization procedures without mixing are not recommendable",error,failure)
    END IF

    IF(scf_env%mixing_method>direct_mixing_nr) THEN
      scf_env%skip_diis = .TRUE.
      scf_env%p_mix_alpha = scf_env%mixing_store%alpha
      IF(scf_env%mixing_store%beta==0.0_dp) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Mixing employing the Kerker damping factor needs BETA /= 0.0",error,failure)
      END IF
    END IF

    IF (scf_env%mixing_method==direct_mixing_nr) THEN
      scf_env%p_mix_alpha = scf_env%mixing_store%alpha
      IF(scf_control%eps_diis<scf_control%eps_scf) THEN
        scf_env%skip_diis = .TRUE.
        CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                       "the DIIS scheme is disabled, since EPS_DIIS < EPS_SCF",&
                       error,failure)
      END IF
    END IF

  END SUBROUTINE qs_scf_ensure_mixing

! *****************************************************************************
!> \brief sets flags for diagonalization and ensure that everything is
!>        allocated  
!> \param scf_env ...
!> \param scf_section ...
!> \param qs_env ...
!> \param scf_control ...
!> \param has_unit_metric ...
!> \param error ...
! *****************************************************************************

  SUBROUTINE qs_scf_ensure_diagonalization(scf_env,scf_section,qs_env,&
                                          scf_control,has_unit_metric,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(scf_control_type), POINTER          :: scf_control
    LOGICAL                                  :: has_unit_metric
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'qs_scf_ensure_diagonalization', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin
    LOGICAL                                  :: do_kpoints, failure, &
                                                need_coeff_b, not_se_or_dftb
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL get_qs_env(qs_env=qs_env,do_kpoints=do_kpoints,dft_control=dft_control,mos=mos,error=error)
    not_se_or_dftb=.NOT.(dft_control%qs_control%dftb .OR. dft_control%qs_control%semi_empirical)
    need_coeff_b=.FALSE.
    scf_env%needs_ortho=.FALSE.

    IF (scf_control%use_diag) THEN
! sanity check whether combinations are allowed
      CALL cp_assert(.NOT.dft_control%restricted,cp_failure_level,cp_assertion_failed,&
            routineP,"OT only for restricted (ROKS)",error,failure)
      SELECT CASE (scf_control%diagonalization%method)
      CASE(diag_ot,diag_block_krylov,diag_block_davidson)
         CALL cp_assert(not_se_or_dftb,cp_failure_level,cp_assertion_failed,&
            routineP,"DFTB and SE not possible with OT diagonalization",error,failure)
      END SELECT
      SELECT CASE (scf_control%diagonalization%method)
! Diagonalization: additional check whether we are in an orthonormal basis
      CASE(diag_standard)
        scf_env%method=general_diag_method_nr
        scf_env%needs_ortho=.NOT.has_unit_metric
        IF (has_unit_metric) THEN
            scf_env%method=special_diag_method_nr
        END IF
! OT Diagonalization: not possible with ROKS
      CASE ( diag_ot ) 
        CALL cp_assert(.NOT.(dft_control%roks),cp_failure_level,cp_assertion_failed,&
              routineP,"ROKS with OT diagonalization not possible",error,failure)
        CALL cp_assert(.NOT.do_kpoints,cp_failure_level,cp_assertion_failed,&
              routineP,"OT diagonalization not possible with kpoint calculations",error,failure)
        scf_env%method=ot_diag_method_nr
        need_coeff_b=.TRUE.
! Block Krylov diagonlization: not possible with ROKS, allocation of additional matrices is needed
      CASE ( diag_block_krylov ) 
        CALL cp_assert(.NOT.(dft_control%roks),cp_failure_level,cp_assertion_failed,&
              routineP,"ROKS with block PF diagonalization not possible",error,failure)
        CALL cp_assert(.NOT.do_kpoints,cp_failure_level,cp_assertion_failed,&
              routineP,"Block Krylov diagonalization not possible with kpoint calculations",error,failure)
        scf_env%method=block_krylov_diag_method_nr
        scf_env%needs_ortho=.TRUE.
        IF(.NOT.ASSOCIATED(scf_env%krylov_space))&
           CALL krylov_space_create(scf_env%krylov_space, scf_section,error=error)
        CALL krylov_space_allocate(scf_env%krylov_space, scf_control, mos, error=error)
! Block davidson diagonlization: allocation of additional matrices is needed
      CASE ( diag_block_davidson )
        CALL cp_assert(.NOT.do_kpoints,cp_failure_level,cp_assertion_failed,&
              routineP,"Block Davidson diagonalization not possible with kpoint calculations",error,failure)
        scf_env%method=block_davidson_diag_method_nr
        IF(.NOT.ASSOCIATED(scf_env%block_davidson_env))&
           CALL block_davidson_env_create(scf_env%block_davidson_env,dft_control%nspins,&
                                          scf_section,error=error)
        DO ispin=1,dft_control%nspins
          CALL block_davidson_allocate(scf_env%block_davidson_env(ispin), mos(ispin)%mo_set, error=error)
        END DO
        need_coeff_b=.TRUE.
      CASE DEFAULT
        CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
             routineP,"Unknown diagonalization method",error,failure)
      END SELECT
! Check if subspace diagonlization is requested: allocation of additional matrices is needed
      IF(scf_control%do_diag_sub) THEN
        scf_env%needs_ortho=.TRUE.
        IF(.NOT.ASSOCIATED(scf_env%subspace_env))&
         CALL diag_subspace_env_create(scf_env%subspace_env,scf_section,&
              dft_control%qs_control%cutoff,error=error)
         CALL diag_subspace_allocate(scf_env%subspace_env, qs_env, mos, error=error)
        CALL cp_assert(.NOT.do_kpoints,cp_failure_level,cp_assertion_failed,&
              routineP,"No subspace diagonlization with kpoint calculation",error,failure)
      END IF
! OT: check if OT is used instead of diagonlization. Not possible with added MOS at the moment
    ELSEIF (scf_control%use_ot) THEN
      scf_env%method=ot_method_nr
      need_coeff_b=.TRUE.
      IF(scf_control%added_mos(1)+scf_control%added_mos(2) > 0) &
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"OT with ADDED_MOS/=0 not implemented", error,failure)
      IF(dft_control%restricted)&
         CALL cp_assert(dft_control%nspins.EQ.2,cp_failure_level,cp_assertion_failed,&
               routineP,"nspin must be 2 for restricted (ROKS)", error,failure)
      CALL cp_assert(.NOT.do_kpoints,cp_failure_level,cp_assertion_failed,&
            routineP,"OT not possible with kpoint calculations",error,failure)
    ELSE
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"OT or DIAGONALIZATION have to be set",error,failure)
    END IF
    DO ispin=1,dft_control%nspins
       mos(ispin)%mo_set%use_mo_coeff_b=need_coeff_b
    END DO

 END SUBROUTINE qs_scf_ensure_diagonalization

! *****************************************************************************
!> \brief computes properties for a given hamilonian using the current wfn 
!> \param qs_env ...
!> \param dft_section ...
!> \param error ...
! *****************************************************************************

 SUBROUTINE qs_scf_compute_properties(qs_env,dft_section,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_compute_properties', &
      routineP = moduleN//':'//routineN

    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_scf_env_type), POINTER           :: scf_env

    CALL get_qs_env(qs_env, scf_env=scf_env, dft_control=dft_control, error=error )

    IF (dft_control%qs_control%semi_empirical) THEN
       CALL scf_post_calculation_se (dft_section, qs_env, error)
    ELSEIF (dft_control%qs_control%dftb) THEN
       CALL scf_post_calculation_dftb (dft_section, scf_env, qs_env, error)
    ELSEIF (dft_control%qs_control%scptb) THEN
       CALL scf_post_calculation_scptb (qs_env, error)
    ELSEIF (dft_control%qs_control%do_kg) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="No properties from PRINT section available for KG methods", &
            error=error, error_level=cp_warning_level)
    ELSEIF (dft_control%qs_control%ofgpw) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="No properties from PRINT section available for OFGPW methods", &
            error=error, error_level=cp_warning_level)
    ELSEIF (dft_control%qs_control%lri_optbas.AND.dft_control%qs_control%gpw) THEN
       CALL optimize_lri_basis(qs_env,error)
    ELSE
       CALL scf_post_calculation_gpw (dft_section, scf_env, qs_env, error)
    END IF

  END SUBROUTINE qs_scf_compute_properties

! *****************************************************************************
!> \brief performs those initialisations that need to be done only once
!>       (e.g. that only depend on the atomic positions)
!>       this will be called in scf
!> \param scf_env ...
!> \param qs_env ...
!> \param scf_section ...
!> \param do_xas_tp ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE init_scf_run(scf_env, qs_env, scf_section, do_xas_tp, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    LOGICAL, INTENT(IN), OPTIONAL            :: do_xas_tp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_scf_run', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, iw, nao, ndep, &
                                                output_unit
    LOGICAL :: dft_plus_u_atom, do_kpoints, failure, init_u_ramping_each_scf, &
      my_transition_potential, s_minus_half_available
    REAL(KIND=dp)                            :: u_ramping
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_h, matrix_ks, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    NULLIFY(scf_control, qs_kind_set, matrix_h, matrix_s, matrix_ks, &
         dft_control, mos, rho)
    failure=.FALSE.
    my_transition_potential = .FALSE.
    IF(PRESENT(do_xas_tp)) my_transition_potential = do_xas_tp
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    NULLIFY(para_env)

    s_minus_half_available = .FALSE.
    CALL get_qs_env(qs_env,&
                    scf_control=scf_control,&
                    dft_control=dft_control,&
                    qs_kind_set=qs_kind_set,&
                    mos=mos,matrix_ks=matrix_ks,&
                    rho=rho,&
                    matrix_h=matrix_h,&
                    matrix_s=matrix_s,&
                    nelectron_total=scf_env%nelectron,&
                    do_kpoints=do_kpoints,&
                    para_env=para_env,&
                    error=error)

    CALL get_qs_kind_set(qs_kind_set=qs_kind_set)

    output_unit = cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
                                       extension=".scfLog",error=error)
    CALL qs_scf_initial_info(output_unit,mos,dft_control)
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

    ! calc ortho matrix
    ndep = 0
    IF (scf_env%needs_ortho)THEN
       CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,scf_env%ortho,error=error)
       IF (scf_env%cholesky_method>cholesky_off) THEN
         CALL cp_fm_cholesky_decompose(scf_env%ortho,error=error)
         IF(scf_env%cholesky_method == cholesky_dbcsr) THEN
            CALL cp_fm_triangular_invert(scf_env%ortho,error=error)
            CALL cp_fm_set_all(scf_env%scf_work2,0.0_dp,error=error)
            CALL cp_fm_to_fm_triangular(scf_env%ortho,scf_env%scf_work2,"U",error=error)
            CALL copy_fm_to_dbcsr(scf_env%scf_work2, scf_env%ortho_dbcsr, error=error)
         ELSE IF(scf_env%cholesky_method==cholesky_inverse) THEN
            CALL cp_fm_to_fm(scf_env%ortho,scf_env%ortho_m1,error=error)
            CALL cp_fm_triangular_invert(scf_env%ortho_m1,error=error)
         END IF
       ELSE
            CALL cp_fm_power(scf_env%ortho,scf_env%scf_work2,-0.5_dp,&
                            scf_control%eps_eigval,ndep,error=error)
            s_minus_half_available = .TRUE.
       END IF
 
       IF (BTEST(cp_print_key_should_output(logger%iter_info,&
              qs_env%input,"DFT%PRINT%AO_MATRICES/ORTHO",error=error),cp_p_file)) THEN
            iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/ORTHO",&
                 extension=".Log",error=error)
            CALL write_fm_with_basis_info(scf_env%ortho,4,6,qs_env,para_env,output_unit=iw,error=error)
            CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
                 "DFT%PRINT%AO_MATRICES/ORTHO", error=error)
       END IF
    END IF

    CALL get_mo_set(mo_set=mos(1)%mo_set,nao=nao)

    ! DFT+U methods based on Lowdin charges need S^(1/2)
    IF (dft_control%dft_plus_u) THEN
      CALL cp_assert(.NOT.do_kpoints,cp_failure_level,cp_assertion_failed,&
            routineP,"DFT+U method not allowed with kpoints",error,failure)
      IF (dft_control%plus_u_method_id == plus_u_lowdin) THEN
        IF (s_minus_half_available) THEN
          CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,scf_env%ortho,scf_env%s_half,&
                                 nao,error=error)
        ELSE
          CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,scf_env%s_half,error=error)
          CALL cp_fm_power(scf_env%s_half,scf_env%scf_work2,0.5_dp,&
                           scf_control%eps_eigval,ndep,error=error)
        END IF
      END IF
      DO ikind=1,SIZE(qs_kind_set)
        qs_kind => qs_kind_set(ikind)
        CALL get_qs_kind(qs_kind=qs_kind,&
                             dft_plus_u_atom=dft_plus_u_atom,&
                             u_ramping=u_ramping,&
                             init_u_ramping_each_scf=init_u_ramping_each_scf)
        IF (dft_plus_u_atom.AND.(u_ramping /= 0.0_dp)) THEN
           IF (init_u_ramping_each_scf) THEN
              CALL set_qs_kind(qs_kind=qs_kind,&
                                   u_minus_j=0.0_dp)
           END IF
        END IF
      END DO
    END IF

    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
               extension=".scfLog",error=error)
    IF (output_unit > 0) THEN
       WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
         "Number of independent orbital functions:",nao - ndep
    END IF
    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
          "PRINT%PROGRAM_RUN_INFO", error=error)

    ! extrapolate outer loop variables
    IF (scf_control%outer_scf%have_scf) THEN
       CALL outer_loop_extrapolate(qs_env,error)
    ENDIF

    ! initializes rho and the mos
    IF( my_transition_potential) THEN
      ! initialize the density with the localized mos
        CALL xas_initialize_rho(qs_env,error=error)
    ELSE
      CALL scf_env_initial_rho_setup(scf_env,qs_env=qs_env,&
           scf_section=scf_section, error=error)
    END IF

    ! Frozen density approximation
    IF (ASSOCIATED(qs_env%wf_history)) THEN
      IF (qs_env%wf_history%interpolation_method_nr==wfi_frozen_method_nr) THEN
        IF (.NOT. ASSOCIATED(qs_env%wf_history%past_states(1)%snapshot)) THEN
          CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp, &
                error=error)
          CALL qs_rho_create(qs_env%wf_history%past_states(1)%snapshot%rho_frozen, error)
          CALL duplicate_rho_type(rho_input=rho, &
                rho_output=qs_env%wf_history%past_states(1)%snapshot%rho_frozen, &
                qs_env=qs_env, error=error)
        END IF
      END IF
    END IF

   !image charge method, calculate image_matrix
   IF(qs_env%qmmm) THEN
     IF(qs_env%qmmm_env_qm%image_charge) THEN
       IF(.NOT.qs_env%qmmm_env_qm%image_charge_pot%coeff_iterative) &
          CALL calculate_image_matrix(image_matrix=qs_env%image_matrix,&
               ipiv=qs_env%ipiv,qs_env=qs_env,qmmm_env=qs_env%qmmm_env_qm,&
               error=error)
     ENDIF
   ENDIF

    CALL timestop(handle)

  END SUBROUTINE init_scf_run
! *****************************************************************************
!> \brief Initializes rho and the mos, so that an scf cycle can start
!> \param scf_env the scf env in which to do the scf
!> \param qs_env the qs env the scf_env lives in
!> \param scf_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_env_initial_rho_setup(scf_env, qs_env, scf_section, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_initial_rho_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extrapolation_method_nr, &
                                                handle, ispin, nmo, &
                                                output_unit
    LOGICAL                                  :: failure, orthogonal_wf
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(mo_coeff, rho, dft_control, para_env)
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)

    CALL get_qs_env(qs_env,&
                    rho=rho,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    error=error)

    IF (.NOT.failure) THEN
       extrapolation_method_nr=wfi_use_guess_method_nr
       IF (ASSOCIATED(qs_env%wf_history)) THEN
          CALL wfi_extrapolate(qs_env%wf_history, &
               qs_env=qs_env, dt=1.0_dp, &
               extrapolation_method_nr=extrapolation_method_nr,&
               orthogonal_wf=orthogonal_wf, error=error)
          ! wfi_use_guess_method_nr the wavefunctions are not yet initialized
          IF ((.NOT.orthogonal_wf).AND.&
              (scf_env%method == ot_method_nr).AND.&
              (.NOT.(extrapolation_method_nr == wfi_use_guess_method_nr))) THEN
             DO ispin=1,SIZE(qs_env%mos)
                CALL get_mo_set(qs_env%mos(ispin)%mo_set, &
                     mo_coeff=mo_coeff, nmo=nmo)
                CALL reorthogonalize_vectors(qs_env, &
                     v_matrix=mo_coeff, n_col=nmo,&
                     error=error)
                CALL set_mo_occupation(mo_set=qs_env%mos(ispin)%mo_set, &
                     smear=qs_env%scf_control%smear, error=error)
             END DO
          END IF
       END IF
       output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
            extension=".scfLog",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A)")&
               "Extrapolation method: "//&
               TRIM(wfi_get_method_label(extrapolation_method_nr,error=error))
       END IF

       CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
       IF (extrapolation_method_nr==wfi_use_guess_method_nr) THEN
          CALL calculate_first_density_matrix(scf_env=scf_env,qs_env=qs_env,error=error)
          IF (.NOT.(qs_env%scf_control%density_guess==densities_guess)) THEN
            CALL qs_rho_update_rho(rho,qs_env=qs_env, error=error)
            CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,&
                 error=error)
          END IF
       END IF

       ! Some preparation for the mixing
       IF(qs_env%scf_env%mixing_method>1) THEN
          IF(dft_control % qs_control%gapw) THEN
            CALL get_qs_env(qs_env=qs_env,&
                 rho_atom_set=rho_atom,error=error)
            CALL mixing_init(qs_env%scf_env%mixing_method,rho,qs_env%scf_env%mixing_store,&
                             para_env,rho_atom=rho_atom,error=error)
          ELSE
            CALL mixing_init(qs_env%scf_env%mixing_method,rho,qs_env%scf_env%mixing_store,&
                             para_env,error=error)
          END IF
       END IF

       DO ispin=1,SIZE(qs_env%mos)!fm->dbcsr
          IF(qs_env%mos(ispin)%mo_set%use_mo_coeff_b) THEN
             CALL copy_fm_to_dbcsr(qs_env%mos(ispin)%mo_set%mo_coeff,&
                  qs_env%mos(ispin)%mo_set%mo_coeff_b,error=error)!fm->dbcsr
          ENDIF
       ENDDO!fm->dbcsr

    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_env_initial_rho_setup

END MODULE qs_scf_initialization 
