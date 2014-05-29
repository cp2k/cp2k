!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains methods used in the context of density fitting
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE admm_methods
  USE admm_types,                      ONLY: admm_create_block_list,&
                                             admm_env_create,&
                                             admm_type
  USE cp_control_types,                ONLY: admm_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_distribution, &
       cp_dbcsr_get_block_p, cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, &
       cp_dbcsr_get_num_blocks, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_type, &
       dbcsr_type_no_symmetry, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: &
       do_admm_basis_set_projection, do_admm_block_aux_basis_off, &
       do_admm_block_aux_basis_on, do_admm_block_density_matrix, &
       do_admm_block_purify_blocked, do_admm_block_purify_full, &
       do_admm_block_purify_off, do_admm_purify_cauchy, &
       do_admm_purify_cauchy_subspace, do_admm_purify_mo_diag, &
       do_admm_purify_mo_no_diag, do_admm_purify_none, &
       do_hfx_potential_coulomb, do_hfx_potential_short, &
       do_hfx_potential_truncated, use_aux_fit_basis_set, use_orb_basis_set, &
       xc_funct_no_shortcut
  USE input_section_types,             ONLY: section_vals_duplicate,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_overlap,                      ONLY: build_overlap_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: admm_fit_mo_coeffs,&
            admm_merge_ks_matrix,&
            admm_merge_mo_derivs,&
            calc_aux_mo_derivs_none,&
            calc_mixed_overlap_force,&
            create_admm_xc_section,&
            admm_calculate_density_matrix,&
            fit_mo_coeffs_diag

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_methods'

!***

  CONTAINS

  SUBROUTINE admm_calculate_density_matrix(admm_env,mo_set,density_matrix, density_matrix_aux,&
                                           ispin, nspins, use_dbcsr,error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix, &
                                                density_matrix_aux
    INTEGER                                  :: ispin, nspins
    LOGICAL, INTENT(IN), OPTIONAL            :: use_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'admm_calculate_density_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%method_id)
    CASE (do_admm_basis_set_projection)
      SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_none)
        CALL calculate_density_matrix(mo_set, density_matrix_aux, use_dbcsr, error=error)
      CASE(do_admm_purify_cauchy)
        CALL calculate_dm_mo_no_diag(admm_env, mo_set, density_matrix_aux, ispin, error)
        CALL purify_density_matrix_cauchy(admm_env, mo_set, density_matrix_aux, ispin, error)
      CASE(do_admm_purify_cauchy_subspace)
        CALL calculate_dm_mo_no_diag(admm_env, mo_set, density_matrix_aux, ispin, error)
      CASE(do_admm_purify_mo_diag)
        CALL calculate_density_matrix(mo_set, density_matrix_aux, use_dbcsr, error=error)
      CASE(do_admm_purify_mo_no_diag)
        CALL calculate_dm_mo_no_diag(admm_env, mo_set, density_matrix_aux, ispin, error)
      END SELECT

    CASE(do_admm_block_density_matrix)
      SELECT CASE(admm_env%block_purification_method)
      CASE(do_admm_block_purify_full)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL blockify_density_matrix(admm_env, density_matrix, density_matrix_aux, ispin, nspins, error)
          CALL purify_dm_cauchy_blocked(admm_env, mo_set, density_matrix_aux, ispin, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      CASE(do_admm_block_purify_off)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL blockify_density_matrix(admm_env, density_matrix, density_matrix_aux, ispin, nspins, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      CASE(do_admm_block_purify_blocked)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL blockify_density_matrix(admm_env, density_matrix, density_matrix_aux, ispin, nspins, error)
          CALL purify_dm_cauchy_blocked(admm_env, mo_set, density_matrix_aux, ispin, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      END SELECT
    END SELECT


    CALL timestop(handle)

  END SUBROUTINE admm_calculate_density_matrix

  SUBROUTINE admm_merge_mo_derivs(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                  mo_derivs_aux_fit, matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_merge_mo_derivs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%method_id)
    CASE (do_admm_basis_set_projection)
      SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_none)
      CASE(do_admm_purify_cauchy)
      CASE(do_admm_purify_cauchy_subspace)
      CASE(do_admm_purify_mo_diag)
        CALL merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit,&
                                  mo_derivs,mo_derivs_aux_fit, matrix_ks_aux_fit,&
                                  error)
      CASE(do_admm_purify_mo_no_diag)
        CALL merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit,&
                                     mo_derivs,mo_derivs_aux_fit, matrix_ks_aux_fit,&
                                     error)
    END SELECT

    CASE(do_admm_block_density_matrix)

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_merge_mo_derivs

  SUBROUTINE admm_merge_ks_matrix(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                             matrix_ks, matrix_ks_aux_fit, matrix_s, &
                                             matrix_p_aux_fit, matrix_p, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit, &
                                                matrix_p
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, jatom
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%method_id)
    CASE (do_admm_basis_set_projection)
      SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_none)
        CALL merge_ks_matrix_none(ispin, admm_env, mo_set, mo_coeff,&
                                  mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
                                  matrix_s, matrix_p_aux_fit, error)
      CASE(do_admm_purify_cauchy)
        CALL merge_ks_matrix_cauchy(ispin, admm_env, mo_set, mo_coeff,&
                                    mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
                                    matrix_s, matrix_p_aux_fit, matrix_p, error)
      CASE(do_admm_purify_cauchy_subspace)
        CALL merge_ks_matrix_cauchy_subspace(ispin, admm_env, mo_set, mo_coeff,&
                                             mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
                                             matrix_s, matrix_p_aux_fit, error)
      CASE(do_admm_purify_mo_diag)
      CASE(do_admm_purify_mo_no_diag)
    END SELECT

    CASE(do_admm_block_density_matrix)
      SELECT CASE(admm_env%block_purification_method)
      CASE(do_admm_block_purify_full)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL merge_ks_matrix_cauchy_blocked(ispin, admm_env, mo_set, mo_coeff,&
                                              mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
                                              matrix_s, matrix_p_aux_fit, matrix_p, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      CASE(do_admm_block_purify_off)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL cp_dbcsr_iterator_start(iter, matrix_ks_aux_fit(ispin)%matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
            CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
            IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
              sparse_block = 0.0_dp
            END IF
          END DO
          CALL cp_dbcsr_iterator_stop(iter)
          CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, 1.0_dp, 1.0_dp, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      CASE(do_admm_block_purify_blocked)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL merge_ks_matrix_cauchy_blocked(ispin, admm_env, mo_set, mo_coeff,&
                                              mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
                                              matrix_s, matrix_p_aux_fit, matrix_p, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      END SELECT
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_merge_ks_matrix

  SUBROUTINE admm_fit_mo_coeffs(qs_env, admm_env, admm_control, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(admm_control_type), POINTER         :: admm_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_fit_mo_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, natom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: admm_block_section, input, &
                                                xc_section

    CALL timeset(routineN,handle)

    NULLIFY(xc_section, admm_block_section)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
      admm_env%method_id = admm_control%method_id
      admm_env%purification_method = admm_control%purification_method
      admm_env%block_purification_method = admm_control%block_purification_method
      admm_env%block_projection_method = admm_control%block_projection_method
      IF( admm_env%method_id == do_admm_block_density_matrix) THEN
        admm_block_section => section_vals_get_subs_vals(input,&
                              "DFT%AUXILIARY_DENSITY_MATRIX_METHOD%BLOCK_DENSITY_MATRIX_METHOD",error=error)
        CALL get_qs_env(qs_env=qs_env,&
                        particle_set=particle_set,&
                        error=error)
        natom = SIZE(particle_set,1)

        CALL admm_create_block_list(admm_block_section, admm_env, natom, error)
      END IF
    END IF

    SELECT CASE(admm_env%method_id)
    CASE (do_admm_basis_set_projection)
      SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_none)
        CALL fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)
      CASE(do_admm_purify_cauchy)
        CALL fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)
      CASE(do_admm_purify_cauchy_subspace)
        CALL fit_mo_coeffs_no_diag(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                   mos, mos_aux_fit, geometry_did_change, error)
      CASE(do_admm_purify_mo_diag)
        CALL fit_mo_coeffs_diag(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)
      CASE(do_admm_purify_mo_no_diag)
        CALL fit_mo_coeffs_no_diag(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                   mos, mos_aux_fit, geometry_did_change, error)
    END SELECT

    CASE(do_admm_block_density_matrix)
      SELECT CASE(admm_env%block_purification_method)
      CASE(do_admm_block_purify_full)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                  mos, mos_aux_fit, geometry_did_change, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      CASE(do_admm_block_purify_off)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                  mos, mos_aux_fit, geometry_did_change, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT
      CASE(do_admm_block_purify_blocked)
        SELECT CASE(admm_env%block_projection_method)
        CASE (do_admm_block_aux_basis_off)
          CALL fit_mo_coeffs_blocked(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                     mos, mos_aux_fit, geometry_did_change, error)
        CASE (do_admm_block_aux_basis_on)
        END SELECT

      END SELECT

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_fit_mo_coeffs

! *****************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>
!> \param admm_env The ADMM env
!> \param para_env The parallel env
!> \param matrix_s_aux_fit the overlap matrix of the auxiliary fitting basis set
!> \param matrix_s_mixed the mixed overlap matrix of the auxiliary fitting basis
!>        set and the orbital basis set
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param geometry_did_change flag to indicate if the geomtry changed
!> \param error
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE fit_mo_coeffs_diag(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, istat, &
                                                nao_aux_fit, nao_orb, nmo, &
                                                nspins
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig_work
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
    END IF

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed
    IF( geometry_did_change ) THEN
      CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv,error)
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*P
      CALL cp_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                    admm_env%A,error)

      !! B=Q^(T)*A
      CALL cp_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                      1.0_dp,admm_env%Q,admm_env%A,0.0_dp,&
                      admm_env%B,error)

   END IF

    ! *** Calculate the mo_coeffs for the fitting basis
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      IF( nmo == 0 ) CYCLE
      !! Lambda = C^(T)*B*C
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
      CALL cp_gemm('N','N',nao_orb,nmo,nao_orb,&
                      1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_gemm('T','N',nmo,nmo,nao_orb,&
                      1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                      admm_env%lambda(ispin)%matrix,error)
      CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error)
      CALL cp_fm_syevd(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%R(ispin)%matrix,&
                       admm_env%eigvals_lambda(ispin)%eigvals%data,error=error)
      ALLOCATE(eig_work(nmo), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO i=1,nmo
        eig_work(i) = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
      END DO
      CALL cp_fm_to_fm(admm_env%R(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error)
      CALL cp_fm_column_scale(admm_env%work_nmo_nmo1(ispin)%matrix,eig_work)
      CALL cp_gemm('N','T',nmo,nmo,nmo,&
                      1.0_dp,admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                      admm_env%lambda_inv_sqrt(ispin)%matrix,error)
      CALL cp_gemm('N','N',nao_orb,nmo,nmo,&
                      1.0_dp,mo_coeff,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,admm_env%work_orb_nmo(ispin)%matrix, 0.0_dp,&
                      mo_coeff_aux_fit,error)
      DEALLOCATE(eig_work, STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END DO

    CALL timestop(handle)

  END SUBROUTINE fit_mo_coeffs_diag

  SUBROUTINE fit_mo_coeffs_no_diag(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                   mos, mos_aux_fit, geometry_did_change, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs_no_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
    END IF

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv,error)
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*Q
      CALL cp_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                    admm_env%A,error)

      !! B=Q^(T)*A
      CALL cp_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                      1.0_dp,admm_env%Q,admm_env%A,0.0_dp,&
                      admm_env%B,error)
    END IF

    DO ispin = 1,nspins
      nmo = admm_env%nmo(ispin)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)

      ! Calculate Lambda inverse for later usage
      CALL cp_gemm('N','N',nao_orb,nmo,nao_orb,&
                      1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix,error)
      CALL cp_gemm('T','N',nmo,nmo,nao_orb,&
                      1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                      admm_env%lambda(ispin)%matrix,error)
      CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error)

      CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
      CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,error=error)
      CALL cp_fm_to_fm(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,error=error)


      !! ** C_hat = AC
      CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_orb,&
                      1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                      admm_env%C_hat(ispin)%matrix,error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE fit_mo_coeffs_no_diag

  SUBROUTINE fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo, nmo_mos, nspins
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occ_num, occ_num_aux
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
    END IF

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv,error)
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*Q
      CALL cp_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                    admm_env%A,error)
    END IF

    DO ispin = 1,nspins
      nmo = admm_env%nmo(ispin)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, occupation_numbers=occ_num,nmo=nmo_mos)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,&
           occupation_numbers=occ_num_aux)

      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                      mo_coeff_aux_fit,error)

      occ_num_aux(1:nmo) = occ_num(1:nmo)
      ! XXXX should only be done first time XXXX
      CALL cp_fm_set_all(admm_env%lambda(ispin)%matrix,0.0_dp,1.0_dp,error)
      CALL cp_fm_set_all(admm_env%lambda_inv(ispin)%matrix,0.0_dp,1.0_dp,error)
      CALL cp_fm_set_all(admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,1.0_dp,error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE fit_mo_coeffs_none

  SUBROUTINE fit_mo_coeffs_blocked(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                   mos, mos_aux_fit, geometry_did_change, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs_blocked', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, jatom, &
                                                nao_aux_fit, nao_orb, nspins
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s_tilde
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
    END IF

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      NULLIFY(matrix_s_tilde)
      ALLOCATE(matrix_s_tilde)
      CALL cp_dbcsr_init (matrix_s_tilde, error)
      CALL cp_dbcsr_create(matrix_s_tilde, 'MATRIX s_tilde', &
           cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_symmetric, &
           cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
           cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
           cp_dbcsr_get_num_blocks(matrix_s_aux_fit(1)%matrix), &
           cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
           cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix), &
           error=error)

      CALL cp_dbcsr_copy(matrix_s_tilde, matrix_s_aux_fit(1)%matrix, error=error)

      CALL cp_dbcsr_iterator_start(iter, matrix_s_tilde)
      DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
        CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
        IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
            sparse_block = 0.0_dp
        END IF
      END DO
      CALL cp_dbcsr_iterator_stop(iter)


      CALL copy_dbcsr_to_fm(matrix_s_tilde,admm_env%S_inv,error)

      CALL cp_dbcsr_deallocate_matrix(matrix_s_tilde,error)

      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*Q
      CALL cp_fm_set_all(admm_env%A, 0.0_dp, 1.0_dp, error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE fit_mo_coeffs_blocked


! *****************************************************************************
!> \brief Calculates the product Kohn-Sham-Matrix x mo_coeff for the auxiliary
!>        basis set and transforms it into the orbital basis. This is needed
!>        in order to use OT
!>
!> \param ispin which spin to transform
!> \param admm_env The ADMM env
!> \param mo_coeff the MO coefficients from the orbital basis set
!> \param mo_coeff_aux_fit the MO coefficients from the auxiliary fitting basis set
!> \param mo_derivs KS x mo_coeff from the orbital basis set to which we add the
!>        auxiliary basis set part
!> \param matrix_ks_aux_fit the Kohn-Sham matrix from the auxiliary fitting basis set
!> \param error
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                  mo_derivs_aux_fit, matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nao_aux_fit, &
                                                nao_orb, nmo, stat
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eig_diff, pole, tmp32, tmp52, &
                                                tmp72, tmp92
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL cp_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%H(ispin)%matrix,error)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    scaling_factor = 2.0_dp*occupation_numbers

    CALL cp_fm_column_scale(admm_env%H(ispin)%matrix,scaling_factor)

    CALL cp_fm_to_fm(admm_env%H(ispin)%matrix, mo_derivs_aux_fit(ispin)%matrix, error=error)

    ! *** Add first term
    CALL cp_gemm('N','T', nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%H(ispin)%matrix,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    CALL cp_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)


    ! *** Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nmo
      DO j=i,nmo
        eig_diff = ( admm_env%eigvals_lambda(ispin)%eigvals%data(i) -&
                     admm_env%eigvals_lambda(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( ABS(eig_diff) < 0.0001_dp ) THEN
          tmp32 = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))**3
          tmp52 = tmp32/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
          tmp72 = tmp52/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
          tmp92 = tmp72/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff

          pole = -0.5_dp*tmp32 + 3.0_dp/8.0_dp*tmp52 - 5.0_dp/16.0_dp*tmp72 + 35.0_dp/128.0_dp*tmp92
          CALL cp_fm_set_element(admm_env%M(ispin)%matrix,i,j,pole,error)
        ELSE
          pole = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
          pole = pole - 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))
          pole = pole/(admm_env%eigvals_lambda(ispin)%eigvals%data(i)-&
                       admm_env%eigvals_lambda(ispin)%eigvals%data(j))
          CALL cp_fm_set_element(admm_env%M(ispin)%matrix,i,j,pole,error)
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(admm_env%M(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,error=error)

    ! *** 2nd term to be added to fm_H

    !! Part 1: B^(T)*C* R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T)
    !! Part 2: B*C*(R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T))^(T)

    ! *** H'*R
    CALL cp_gemm('N','N', nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%H(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    ! *** A^(T)*H'*R
    CALL cp_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)
    ! *** c^(T)*A^(T)*H'*R
    CALL cp_gemm('T','N', nmo, nmo, nao_orb,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo1(ispin)%matrix,error)
    ! *** R^(T)*c^(T)*A^(T)*H'*R
    CALL cp_gemm('T','N', nmo, nmo, nmo,&
                    1.0_dp,admm_env%R(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo2(ispin)%matrix,error)
    ! *** R^(T)*c^(T)*A^(T)*H'*R x M
    CALL cp_fm_schur_product(admm_env%work_nmo_nmo2(ispin)%matrix,&
                             admm_env%M(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,error)
    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M)
    CALL cp_gemm('N','N', nmo, nmo, nmo,&
                    1.0_dp,admm_env%R(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo2(ispin)%matrix,error)

    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_gemm('N','T', nmo, nmo, nmo,&
                    1.0_dp,admm_env%work_nmo_nmo2(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                    admm_env%R_schur_R_t(ispin)%matrix,error)

    ! *** B^(T)*c
    CALL cp_gemm('T','N', nao_orb, nmo, nao_orb,&
                    1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)

    ! *** Add first term to fm_H
    ! *** B^(T)*c* R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_gemm('N','N', nao_orb, nmo, nmo,&
                    1.0_dp,admm_env%work_orb_nmo(ispin)%matrix,admm_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)

    ! *** Add second term to fm_H
    ! *** B*C *[ R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)]^(T)
    CALL cp_gemm('N','T', nao_orb, nmo, nmo,&
                    1.0_dp,admm_env%work_orb_nmo(ispin)%matrix,admm_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)

    DO i = 1,SIZE(scaling_factor)
      scaling_factor(i) = 1.0_dp/scaling_factor(i)
    END DO

    CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin)%matrix,scaling_factor)

    CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(ispin)%matrix,1.0_dp,admm_env%mo_derivs_tmp(ispin)%matrix,error)

    DEALLOCATE(scaling_factor, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs_diag

  SUBROUTINE calc_aux_mo_derivs_none(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit,  &
                                  matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_aux_mo_derivs_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo, stat
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    ! just calculate the mo derivs in the aux basis  
    ! only needs to be done on the converged ks matrix for the force calc
    ! Note with OT and purification NONE, the merging of the derivs
    ! happens implicitly because the KS matrices have been already been merged
    ! and adding them here would be double counting.

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL cp_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%H(ispin)%matrix,error)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    scaling_factor = 2.0_dp*occupation_numbers
    CALL cp_fm_column_scale(admm_env%H(ispin)%matrix,scaling_factor)

    DEALLOCATE(scaling_factor, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE calc_aux_mo_derivs_none

  SUBROUTINE merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                  mo_derivs_aux_fit, matrix_ks_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs_no_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo, stat
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    scaling_factor = 0.5_dp


    !! ** calculate first part
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix,error)
    CALL cp_gemm('T', 'N',  nao_orb, nmo, nao_aux_fit,&
                    2.0_dp,admm_env%A,admm_env%work_aux_nmo2(ispin)%matrix,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)
    !! ** calculate second part
    CALL cp_gemm('T', 'N',  nmo, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%work_aux_nmo2(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_orb,error)
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%work_orb_orb,0.0_dp,&
                    admm_env%work_aux_orb,error)
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%S,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    CALL cp_gemm('T', 'N',  nao_orb, nmo, nao_aux_fit,&
                    -2.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)

    CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin)%matrix,scaling_factor)

    CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(ispin)%matrix,1.0_dp,admm_env%mo_derivs_tmp(ispin)%matrix,error)

    DEALLOCATE(scaling_factor, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs_no_diag


! *****************************************************************************
!> \brief Calculates contribution of forces due to basis transformation
!>
!>        dE/dR = dE/dC'*dC'/dR
!>        dE/dC = Ks'*c'*occ = H'
!>
!>        dC'/dR = - tr(A*lambda^(-1/2)*H'^(T)*S^(-1) * dS'/dR)
!>                 - tr(A*C*Y^(T)*C^(T)*Q^(T)*A^(T) * dS'/dR)
!>                 + tr(C*lambda^(-1/2)*H'^(T)*S^(-1) * dQ/dR)
!>                 + tr(A*C*Y^(T)*c^(T) * dQ/dR)
!>                 + tr(C*Y^(T)*C^(T)*A^(T) * dQ/dR)
!>
!>        where
!>
!>        A = S'^(-1)*Q
!>        lambda = C^(T)*B*C
!>        B = Q^(T)*A
!>        Y = R*[ (R^(T)*C^(T)*A^(T)*H'*R) xx M ]*R^(T)
!>        lambda = R*D*R^(T)
!>        Mij = Poles-Matrix (see above)
!>        xx = schur product
!>
!> \param qs_env the QS environment
!> \param ispin which spin to transform
!> \param admm_env The ADMM env
!> \param mo_coeff the MO coefficients from the orbital basis set
!> \param matrix_s_aux_fit overlap matrix from auxiliary fitting basis
!> \param matrix_s_aux_fit_vs_orb mixed orbital/aux fit overlap matrix
!> \param error
!>
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_mixed_overlap_force(qs_env, para_env, ispin, admm_env, mo_coeff,  &
                                      matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, logger,&
                                      iw, use_virial, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_logger_type), POINTER            :: logger
    INTEGER                                  :: iw
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mixed_overlap_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                neighbor_list_id, nmo
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type), POINTER             :: matrix_w_q, matrix_w_s

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(matrix_w_q, matrix_w_s)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    ! *** forces are only implemented for mo_diag or none and basis_projection ***
    IF (.NOT.(admm_env%method_id==do_admm_basis_set_projection .AND. &
              (admm_env%purification_method==do_admm_purify_mo_diag .OR. &
               admm_env%purification_method==do_admm_purify_none))) THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ! *** Create sparse work matrices
    CALL get_qs_env(qs_env=qs_env,&
                    neighbor_list_id=neighbor_list_id,&
                    error=error)

    ALLOCATE(matrix_w_s)
    CALL cp_dbcsr_init (matrix_w_s, error)
    CALL cp_dbcsr_create(matrix_w_s, 'W MATRIX AUX S', &
         cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_no_symmetry, &
         cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix), &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s,qs_env%sab_aux_fit_asymm,error=error)

    ALLOCATE(matrix_w_q)
    CALL cp_dbcsr_init(matrix_w_q, error=error)
    CALL cp_dbcsr_copy(matrix_w_q,matrix_s_aux_fit_vs_orb(1)%matrix,&
                    "W MATRIX AUX Q",error=error)

    ! *** S'^(-T)*H'
    IF (.NOT. admm_env%purification_method==do_admm_purify_none) THEN
    CALL cp_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,qs_env%mo_derivs_aux_fit(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    ELSE

    CALL cp_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%H(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)
    END IF

    ! *** S'^(-T)*H'*Lambda^(-T/2)
    CALL cp_gemm('N', 'T',  nao_aux_fit,nmo, nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix,error)

    ! *** C*Lambda^(-1/2)*H'^(T)*S'^(-1) minus sign due to force = -dE/dR
    CALL cp_gemm('N', 'T',  nao_aux_fit, nao_orb, nmo,&
                    -1.0_dp,admm_env%work_aux_nmo2(ispin)%matrix,mo_coeff,0.0_dp,&
                    admm_env%work_aux_orb,error)

    ! *** A*C*Lambda^(-1/2)*H'^(T)*S'^(-1), minus sign to recover from above
    CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    -1.0_dp,admm_env%work_aux_orb,admm_env%A,0.0_dp,&
                    admm_env%work_aux_aux,error)


    IF (.NOT. (admm_env%purification_method==do_admm_purify_none)) THEN
       ! *** C*Y
       CALL cp_gemm('N', 'N',  nao_orb, nmo, nmo,&
                    1.0_dp,mo_coeff,admm_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)
       ! *** C*Y^(T)*C^(T)
       CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_orb,error)
       ! *** A*C*Y^(T)*C^(T) Add to work aux_orb, minus sign due to force = -dE/dR
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    -1.0_dp,admm_env%A,admm_env%work_orb_orb,1.0_dp,&
                    admm_env%work_aux_orb,error)

       ! *** C*Y^(T)
       CALL cp_gemm('N', 'T',  nao_orb, nmo, nmo,&
                    1.0_dp,mo_coeff,admm_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix,error)
       ! *** C*Y*C^(T)
       CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_orb,error)
       ! *** A*C*Y*C^(T) Add to work aux_orb, minus sign due to -dE/dR
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    -1.0_dp,admm_env%A,admm_env%work_orb_orb,1.0_dp,&
                    admm_env%work_aux_orb,error)
    END IF

    ! *** copy to sparse matrix
    CALL copy_fm_to_dbcsr(admm_env%work_aux_orb, matrix_w_q,keep_sparsity=.TRUE.,&
         error=error)

    IF (.NOT. (admm_env%purification_method==do_admm_purify_none)) THEN
       ! *** A*C*Y^(T)*C^(T)
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    1.0_dp,admm_env%A,admm_env%work_orb_orb,0.0_dp,&
                    admm_env%work_aux_orb,error)
       ! *** A*C*Y^(T)*C^(T)*A^(T) add to aux_aux, minus sign cancels
       CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    1.0_dp,admm_env%work_aux_orb,admm_env%A,1.0_dp,&
                    admm_env%work_aux_aux,error)
    END IF

    ! *** copy to sparse matrix
    CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s,keep_sparsity=.TRUE.,&
         error=error)


    ! *** This can be done in one call w_total = w_alpha + w_beta
    CALL build_overlap_matrix(qs_env, nderivative=-1, &
         basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_aux_fit_basis_set, &
         sab_nl=qs_env%sab_aux_fit_asymm, calculate_forces=.TRUE., admm=.TRUE., &
         matrix_p=matrix_w_s, error=error)
    CALL build_overlap_matrix(qs_env, nderivative=-1, &
         basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_orb_basis_set, &
         sab_nl=qs_env%sab_aux_fit_vs_orb, calculate_forces=.TRUE., admm=.TRUE., &
         matrix_p=matrix_w_q, error=error)


    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",&
            extension=".Log",error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_w_s,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT", error=error)
    END IF
    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",error=error),cp_p_file)) THEN
       iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",&
            extension=".Log",error=error)
       CALL cp_dbcsr_write_sparse_matrix(matrix_w_q,4,6,qs_env,para_env,output_unit=iw,error=error)
       CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
            "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT", error=error)
    END IF

    ! *** Deallocated weighted density matrices
    CALL cp_dbcsr_deallocate_matrix(matrix_w_s,error)
    CALL cp_dbcsr_deallocate_matrix(matrix_w_q,error)

    CALL timestop(handle)

  END SUBROUTINE calc_mixed_overlap_force

! *****************************************************************************
!> \brief This routine modifies the xc section depending on the potential type
!>        used for the HF exchange and the resulting correction term. Currently
!>        three types of corrections are implemented:
!>
!>        coulomb:     Ex,hf = Ex,hf' + (PBEx-PBEx')
!>        shortrange:  Ex,hf = Ex,hf' + (XWPBEX-XWPBEX')
!>        truncated:   Ex,hf = Ex,hf' + ( (XWPBEX0-PBE_HOLE_TC_LR) -(XWPBEX0-PBE_HOLE_TC_LR)' )
!>
!>        with ' denoting the auxiliary basis set and
!>
!>          PBEx:           PBE exchange functional
!>          XWPBEX:         PBE exchange hole for short-range potential (erfc(omega*r)/r)
!>          XWPBEX0:        PBE exchange hole for standard coulomb potential
!>          PBE_HOLE_TC_LR: PBE exchange hole for longrange truncated coulomb potential
!>
!>
!> \param qs_env the qs environment
!> \param xc_section the original xc_section
!> \param admm_env the ADMM environment
!> \param error
!>
!> \par History
!>      12.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_admm_xc_section(qs_env, xc_section, admm_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_admm_xc_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: hfx_potential_type, ifun, nfun
    LOGICAL                                  :: failure, funct_found
    REAL(dp)                                 :: cutoff_radius, hfx_fraction, &
                                                omega, scale_x
    TYPE(section_vals_type), POINTER         :: xc_fun, xc_fun_section

    NULLIFY(admm_env%xc_section_aux, admm_env%xc_section_primary)
    failure=.FALSE.

    !! ** Duplicate existing xc-section

    CALL section_vals_duplicate(xc_section,admm_env%xc_section_aux,error=error)
    CALL section_vals_duplicate(xc_section,admm_env%xc_section_primary,error=error)
    !** Now modify the auxiliary basis
    !** First remove all functionals
    xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux,"XC_FUNCTIONAL",error=error)

    !* Overwrite possible shortcut
    CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                              i_val=xc_funct_no_shortcut,error=error)

    !** Get number of Functionals in the list
    ifun = 0
    nfun = 0
    DO
      ifun = ifun+1
      xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
      IF (.NOT.ASSOCIATED(xc_fun)) EXIT
      nfun = nfun + 1
    END DO

    ifun = 0
    DO ifun = 1,nfun
      xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=1,error=error)
      IF (.NOT.ASSOCIATED(xc_fun)) EXIT
      CALL section_vals_remove_values(xc_fun, error=error)
    END DO

    hfx_potential_type = qs_env%x_data(1,1)%potential_parameter%potential_type
    hfx_fraction = qs_env%x_data(1,1)%general_parameter%fraction

    !! ** Add functionals evaluated with auxiliary basis
    SELECT CASE (hfx_potential_type)
    CASE (do_hfx_potential_coulomb)
      CALL section_vals_val_set(xc_fun_section,"PBE%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                r_val=-hfx_fraction,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_C",&
                                r_val=0.0_dp,error=error)
    CASE (do_hfx_potential_short)
      omega =  qs_env%x_data(1,1)%potential_parameter%omega
      CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                r_val=-hfx_fraction,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                r_val=0.0_dp,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%OMEGA",&
                                r_val=omega,error=error)
    CASE (do_hfx_potential_truncated)
      cutoff_radius = qs_env%x_data(1,1)%potential_parameter%cutoff_radius
      CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                r_val=hfx_fraction,error=error)
      CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                r_val=cutoff_radius,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                l_val=.TRUE.,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                r_val=0.0_dp,error=error)
      CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                r_val=-hfx_fraction,error=error)
    CASE DEFAULT
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT


    !** Now modify the functionals for the primary basis
    xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
    !* Overwrite possible shortcut
    CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                              i_val=xc_funct_no_shortcut,error=error)


    SELECT CASE (hfx_potential_type)
    CASE (do_hfx_potential_coulomb)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "PBE" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"PBE%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                  r_val=hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_C",&
                                  r_val=0.0_dp,error=error)
      ELSE
        CALL section_vals_val_get(xc_fun_section,"PBE%SCALE_X",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x + hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                  r_val=scale_x,error=error)
      END IF
    CASE (do_hfx_potential_short)
      omega =  qs_env%x_data(1,1)%potential_parameter%omega
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "XWPBE" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=0.0_dp,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%OMEGA",&
                                  r_val=omega,error=error)
      ELSE
        CALL section_vals_val_get(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x + hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=scale_x,error=error)
      END IF
    CASE (do_hfx_potential_truncated)
      cutoff_radius =  qs_env%x_data(1,1)%potential_parameter%cutoff_radius
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "PBE_HOLE_T_C_LR" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                  r_val=-hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                  r_val=cutoff_radius,error=error)

      ELSE
        CALL section_vals_val_get(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x - hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                  r_val=scale_x,error=error)
        CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                  r_val=cutoff_radius,error=error)
      END IF
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        IF( xc_fun%section%name == "XWPBE" ) THEN
          funct_found = .TRUE.
        END IF
      END DO
      IF( .NOT. funct_found ) THEN
        CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                  l_val=.TRUE.,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=hfx_fraction,error=error)
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                  r_val=0.0_dp,error=error)

      ELSE
        CALL section_vals_val_get(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=scale_x,error=error)
        scale_x = scale_x + hfx_fraction
        CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                  r_val=scale_x,error=error)
      END IF

    END SELECT


    IF( 1==0 ) THEN
      WRITE(*,*) "primary"
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT

        scale_x=-1000.0_dp
        IF(xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X",&
                                      r_val=scale_x,error=error)
        END IF
        IF( xc_fun%section%name == "XWPBE" ) THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X0",&
                                    r_val=hfx_fraction,error=error)

           WRITE(*,*) xc_fun%section%name, scale_x, hfx_fraction
        ELSE
          WRITE(*,*) xc_fun%section%name, scale_x
        END IF
      END DO

      WRITE(*,*) "auxiliary"
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux,"XC_FUNCTIONAL",error=error)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        scale_x=-1000.0_dp
        IF(xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X",&
                                      r_val=scale_x,error=error)
        END IF
        IF( xc_fun%section%name == "XWPBE" ) THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X0",&
                                    r_val=hfx_fraction,error=error)

           WRITE(*,*) xc_fun%section%name, scale_x, hfx_fraction
        ELSE
          WRITE(*,*) xc_fun%section%name, scale_x
        END IF
      END DO
    END IF

  END SUBROUTINE create_admm_xc_section

  SUBROUTINE merge_ks_matrix_cauchy_subspace(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                             matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'merge_ks_matrix_cauchy_subspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    !! Calculate Lambda^{-2}
    CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix, error=error)
    CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
    CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
    !! Symmetrize the guy
    CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv2(ispin)%matrix,error=error)
    !! Take square
    CALL cp_gemm('N', 'T',  nmo, nmo, nmo,&
                    1.0_dp,admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%lambda_inv2(ispin)%matrix,error)

    !! ** C_hat = AC
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_orb,&
                    1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                    admm_env%C_hat(ispin)%matrix,error)

    !! calc P_tilde from C_hat
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)

    CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,  admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%P_tilde(ispin)%matrix,error)


    !! ** C_hat*Lambda^{-2}
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv2(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix,error)

    !! ** C_hat*Lambda^{-2}*C_hat^T
    CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%C_hat(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)


    !! ** S*C_hat*Lambda^{-2}*C_hat^T
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%S,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux2,error)


    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    !! ** S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux2,admm_env%K(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)

    !! ** P_tilde*S
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%P_tilde(ispin)%matrix,admm_env%S,0.0_dp,&
                    admm_env%work_aux_aux2,error)


    !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    -1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                    admm_env%work_aux_aux3,error)


    !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S+S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
    CALL cp_fm_scale_and_add(1.0_dp,admm_env%work_aux_aux3,1.0_dp,admm_env%work_aux_aux,error)


    !! first_part*A
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux3,admm_env%A,0.0_dp,&
                    admm_env%work_aux_orb,error)

    !! + first_part^T*A
    CALL cp_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux3,admm_env%A,1.0_dp,&
                    admm_env%work_aux_orb,error)



    !! A^T*(first+seccond)=H
    CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix), &
         cp_dbcsr_get_data_size(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)


    CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.,&
         error=error)

    CALL cp_gemm('N', 'N',  nao_orb, nmo, nao_orb,&
                    1.0_dp,admm_env%work_orb_orb,mo_coeff,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix,error)


    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy_subspace

  SUBROUTINE merge_ks_matrix_cauchy_blocked(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                            matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit,&
                                            matrix_p, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit, &
                                                matrix_p
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'merge_ks_matrix_cauchy_blocked', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iatom, j, &
                                                jatom, nao_aux_fit, nao_orb, &
                                                nmo, nspins
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eig_diff, pole, tmp
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)
    nspins = SIZE(admm_env%P_to_be_purified)

    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux, error=error)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2, error=error)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux,error=error)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3, error=error)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,error=error)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T", error=error)

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix, error=error)

    ! *** Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      DO j=i,nao_aux_fit
        eig_diff = ( admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) -&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( ABS(eig_diff) == 0.0_dp ) THEN
          pole = delta(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
          CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole,error)
        ELSE
          pole = 1.0_dp/(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-&
                         admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j))
          tmp = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
          tmp = tmp - Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j)-0.5_dp)
          pole = tmp*pole
          CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole,error)
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)


    !! S^(-1)*R
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)
    !! K*S^(-1)*R
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux2,error)
    !! R^T*S^(-1)*K*S^(-1)*R
    CALL cp_gemm('T', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                    admm_env%work_aux_aux3,error)
    !! R^T*S^(-1)*K*S^(-1)*R x M
    CALL cp_fm_schur_product(admm_env%work_aux_aux3, admm_env%M_purify(ispin)%matrix,&
                             admm_env%work_aux_aux,error)

    !! R^T*A
    CALL cp_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%R_purify(ispin)%matrix, admm_env%A, 0.0_dp,&
                    admm_env%work_aux_orb,error)

    !! (R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_orb, 0.0_dp,&
                    admm_env%work_aux_orb2,error)
    !! A^T*R*(R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
    CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%work_aux_orb, admm_env%work_aux_orb2, 0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric,&
 cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix), cp_dbcsr_get_data_size(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), &
         error=error)

    CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE., &
         error=error)

    ! ** now loop through the list and nullify blocks
    CALL cp_dbcsr_iterator_start(iter, matrix_k_tilde)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
      CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
      IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
         sparse_block = 0.0_dp
      END IF
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy_blocked


  SUBROUTINE merge_ks_matrix_none(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                  matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    !! K*A
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%A,0.0_dp,&
                    admm_env%work_aux_orb,error)
    !! A^T*K*A
    CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix), cp_dbcsr_get_data_size(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.,&
         error=error)

    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_none

  SUBROUTINE calculate_dm_mo_no_diag(admm_env,mo_set,density_matrix,ispin,error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dm_mo_no_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nmo
    REAL(KIND=dp)                            :: alpha

    CALL timeset(routineN,handle)


    CALL cp_dbcsr_set(density_matrix,0.0_dp,error=error)
    nao_aux_fit = admm_env%nao_aux_fit
    nmo = admm_env%nmo(ispin)
    CALL cp_fm_to_fm(admm_env%C_hat(ispin)%matrix, admm_env%work_aux_nmo(ispin)%matrix, error=error)
    CALL cp_fm_column_scale(admm_env%work_aux_nmo(ispin)%matrix,mo_set%occupation_numbers(1:mo_set%homo))

    CALL cp_gemm('N','N',nao_aux_fit,nmo,nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix,error)


    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
      alpha=1.0_dp
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha,error=error)
    ELSE
      alpha=1.0_dp
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE calculate_dm_mo_no_diag

  SUBROUTINE purify_dm_cauchy_blocked(admm_env,mo_set,density_matrix,ispin,error)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'purify_dm_cauchy_blocked', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    REAL(KIND=dp)                            :: pole
    TYPE(cp_fm_type), POINTER                :: mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    nspins = SIZE(admm_env%P_to_be_purified)

    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff_aux_fit)

    !! * For the time beeing, get the P to be purified from the mo_coeffs
    !! * This needs to be replaced with the a block modified P

!    CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nmo,&
!                     1.0_dp,mo_coeff_aux_fit,mo_coeff_aux_fit,0.0_dp,&
!                     admm_env%P_to_be_purified(ispin)%matrix,error)

    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux, error=error)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2, error=error)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux,error=error)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3, error=error)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,error=error)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T", error=error)

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix, error=error)

    ! *** Construct Matrix M for Hadamard Product
    CALL cp_fm_set_all(admm_env%M_purify(ispin)%matrix,0.0_dp,error=error)
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      pole = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
      CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,i,pole,error)
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL copy_dbcsr_to_fm(density_matrix,admm_env%work_aux_aux3,error)
    CALL cp_fm_upper_to_full(admm_env%work_aux_aux3,admm_env%work_aux_aux,error=error)

    ! ** S^(-1)*R
    CALL cp_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)
    ! ** S^(-1)*R*M
    CALL cp_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%M_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux2,error)
    ! ** S^(-1)*R*M*R^T*S^(-1)
    CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux2,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux3,error)

    CALL copy_fm_to_dbcsr(admm_env%work_aux_aux3, density_matrix,keep_sparsity=.TRUE., error=error)

    IF( nspins == 1 ) THEN
      CALL cp_dbcsr_scale(density_matrix, 2.0_dp, error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE purify_dm_cauchy_blocked

  SUBROUTINE blockify_density_matrix(admm_env,density_matrix, density_matrix_aux,&
                                     ispin, nspins, error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix, &
                                                density_matrix_aux
    INTEGER                                  :: ispin, nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'blockify_density_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, jatom
    LOGICAL                                  :: failure, found
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block, sparse_block_aux
    TYPE(cp_dbcsr_iterator)                  :: iter

    failure = .FALSE.

    CALL timeset(routineN,handle)

    ! ** set blocked density matrix to 0
    CALL cp_dbcsr_set(density_matrix_aux, 0.0_dp, error)

    ! ** now loop through the list and copy corresponding blocks
    CALL cp_dbcsr_iterator_start(iter, density_matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
      CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
      IF( admm_env%block_map(iatom,jatom) == 1 ) THEN
        CALL cp_dbcsr_get_block_p(density_matrix_aux,&
                               row=iatom,col=jatom,BLOCK=sparse_block_aux,found=found)
        IF( found ) THEN
          sparse_block_aux = sparse_block
        END IF

      END IF
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL copy_dbcsr_to_fm(density_matrix_aux,admm_env%P_to_be_purified(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_orb_orb2,error=error)

    IF( nspins == 1 ) THEN
      CALL cp_fm_scale(0.5_dp, admm_env%P_to_be_purified(ispin)%matrix, error)
    END IF

    CALL timestop(handle)
  END SUBROUTINE blockify_density_matrix

  SUBROUTINE merge_ks_matrix_cauchy(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                    matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit, matrix_p, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit, &
                                                matrix_p
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_cauchy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nao_aux_fit, &
                                                nao_orb, nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eig_diff, pole, tmp
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    !** Get P from mo_coeffs, otherwise we have troubles with occupation numbers ...
    CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp, mo_coeff, mo_coeff, 0.0_dp,&
                    admm_env%work_orb_orb,error)

    !! A*P
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp,&
                    admm_env%work_aux_orb2,error)
    !! A*P*A^T
    CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    1.0_dp, admm_env%work_aux_orb2, admm_env%A, 0.0_dp,&
                    admm_env%P_to_be_purified(ispin)%matrix,error)


    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux, error=error)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2, error=error)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux,error=error)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3, error=error)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,error=error)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T", error=error)

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix, error=error)

    ! *** Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      DO j=i,nao_aux_fit
        eig_diff = ( admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) -&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( ABS(eig_diff) == 0.0_dp ) THEN
          pole = delta(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
          CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole,error)
        ELSE
          pole = 1.0_dp/(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-&
                         admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j))
          tmp = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
          tmp = tmp - Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j)-0.5_dp)
          pole = tmp*pole
          CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole,error)
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)


    !! S^(-1)*R
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)
    !! K*S^(-1)*R
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux2,error)
    !! R^T*S^(-1)*K*S^(-1)*R
    CALL cp_gemm('T', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                    admm_env%work_aux_aux3,error)
    !! R^T*S^(-1)*K*S^(-1)*R x M
    CALL cp_fm_schur_product(admm_env%work_aux_aux3, admm_env%M_purify(ispin)%matrix,&
                             admm_env%work_aux_aux,error)

    !! R^T*A
    CALL cp_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%R_purify(ispin)%matrix, admm_env%A, 0.0_dp,&
                    admm_env%work_aux_orb,error)

    !! (R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_orb, 0.0_dp,&
                    admm_env%work_aux_orb2,error)
    !! A^T*R*(R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
    CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%work_aux_orb, admm_env%work_aux_orb2, 0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_size( matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)

    CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE., error=error)

    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy

 SUBROUTINE purify_density_matrix_cauchy(admm_env,mo_set,density_matrix,ispin,error)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'purify_density_matrix_cauchy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    REAL(KIND=dp)                            :: pole
    TYPE(cp_fm_type), POINTER                :: mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    nspins = SIZE(admm_env%P_to_be_purified)

    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff_aux_fit)

    !! * For the time beeing, get the P to be purified from the mo_coeffs
    !! * This needs to be replaced with the a block modified P

    CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nmo,&
                     1.0_dp,mo_coeff_aux_fit,mo_coeff_aux_fit,0.0_dp,&
                     admm_env%P_to_be_purified(ispin)%matrix,error)

    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux, error=error)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2, error=error)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux,error=error)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3, error=error)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,error=error)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T", error=error)

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix, error=error)

    ! *** Construct Matrix M for Hadamard Product
    CALL cp_fm_set_all(admm_env%M_purify(ispin)%matrix,0.0_dp,error=error)
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      pole = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
      CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,i,pole,error)
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL copy_dbcsr_to_fm(density_matrix,admm_env%work_aux_aux3,error)
    CALL cp_fm_upper_to_full(admm_env%work_aux_aux3,admm_env%work_aux_aux,error=error)

    ! ** S^(-1)*R
    CALL cp_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)
    ! ** S^(-1)*R*M
    CALL cp_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%M_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux2,error)
    ! ** S^(-1)*R*M*R^T*S^(-1)
    CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux2,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux3,error)

    CALL copy_fm_to_dbcsr(admm_env%work_aux_aux3, density_matrix,keep_sparsity=.TRUE.,error=error)

    IF( nspins == 1 ) THEN
      CALL cp_dbcsr_scale(density_matrix, 2.0_dp, error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE purify_density_matrix_cauchy


  FUNCTION delta(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: delta

    IF( x == 0.0_dp) THEN
      delta = 1.0_dp
    ELSE
      delta = 0.0_dp
    END IF

  END FUNCTION delta

  FUNCTION Heaviside(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: Heaviside

    IF( x < 0.0_dp ) THEN
      Heaviside = 0.0_dp
    ELSE
      Heaviside = 1.0_dp
    END IF
  END FUNCTION Heaviside

END MODULE admm_methods
