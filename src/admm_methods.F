!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains ADMM methods which require molecular orbitals
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE admm_methods
  USE admm_types,                      ONLY: admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
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
  USE input_constants,                 ONLY: do_admm_purify_cauchy,&
                                             do_admm_purify_cauchy_subspace,&
                                             do_admm_purify_mo_diag,&
                                             do_admm_purify_mo_no_diag,&
                                             do_admm_purify_none,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: add_qs_force,&
                                             qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_force
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: admm_mo_calc_rho_aux,&
            admm_mo_merge_ks_matrix,&
            admm_mo_merge_derivs,&
            calc_mixed_overlap_force,&
            calc_aux_mo_derivs_none

  PUBLIC :: admm_fit_mo_coeffs

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_methods'

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE admm_mo_calc_rho_aux(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'admm_mo_calc_rho_aux', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: s_mstruct_changed
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r_aux
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb, &
                                                rho_ao, rho_ao_aux
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g_aux, rho_r_aux
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit

    CALL timeset(routineN,handle)

    NULLIFY(ks_env, admm_env, mos, mos_aux_fit, matrix_s_aux_fit,&
            matrix_s_aux_fit_vs_orb, rho, rho_aux_fit, para_env)
    NULLIFY(rho_g_aux, rho_r_aux, rho_ao, rho_ao_aux, tot_rho_r_aux)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    mos_aux_fit=mos_aux_fit,&
                    mos=mos,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    para_env=para_env,&
                    s_mstruct_changed=s_mstruct_changed,&
                    rho=rho,&
                    rho_aux_fit=rho_aux_fit,&
                    error=error)

    CALL qs_rho_get(rho, rho_ao=rho_ao, error=error)
    CALL qs_rho_get(rho_aux_fit,&
                    rho_ao=rho_ao_aux,&
                    rho_g=rho_g_aux,&
                    rho_r=rho_r_aux,&
                    tot_rho_r=tot_rho_r_aux,&
                    error=error)

    ! convert mos from full to dbcsr matrices
    DO ispin=1, dft_control%nspins
       IF(mos(ispin)%mo_set%use_mo_coeff_b) THEN
          CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b,mos(ispin)%mo_set%mo_coeff,error=error)
       ENDIF
    ENDDO

    ! fit mo coeffcients
    CALL admm_fit_mo_coeffs(admm_env, para_env, &
                            matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,&
                            mos, mos_aux_fit, s_mstruct_changed, error=error)

    DO ispin=1,dft_control%nspins
       IF(admm_env%block_dm) THEN
          CALL blockify_density_matrix(admm_env,&
                                       density_matrix=rho_ao(ispin)%matrix,&
                                       density_matrix_aux=rho_ao_aux(ispin)%matrix,&
                                       ispin=ispin,&
                                       nspins=dft_control%nspins,&
                                       error=error)
       ELSE
          CALL calculate_dm_mo_no_diag(admm_env,&
                                       mo_set=mos(ispin)%mo_set,&
                                       density_matrix=rho_ao_aux(ispin)%matrix,&
                                       ispin=ispin,&
                                       error=error)
       ENDIF

       IF(admm_env%purification_method == do_admm_purify_cauchy)&
         CALL purify_dm_cauchy(admm_env,&
                               mo_set=mos_aux_fit(ispin)%mo_set,&
                               density_matrix=rho_ao_aux(ispin)%matrix,&
                               ispin=ispin,&
                               blocked=admm_env%block_dm,&
                               error=error)

       CALL calculate_rho_elec(ks_env=ks_env,&
                               matrix_p=rho_ao_aux(ispin)%matrix,&
                               rho=rho_r_aux(ispin),&
                               rho_gspace=rho_g_aux(ispin),&
                               total_rho=tot_rho_r_aux(ispin),&
                               soft_valid=.FALSE.,&
                               basis_set_id=use_aux_fit_basis_set,&
                               error=error)
    END DO

    CALL qs_rho_set(rho_aux_fit, rho_r_valid=.TRUE., rho_g_valid=.TRUE., error=error)

    CALL timestop(handle)

  END SUBROUTINE admm_mo_calc_rho_aux

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE admm_mo_merge_ks_matrix(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_mo_merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(admm_type), POINTER                 :: admm_env

    CALL timeset(routineN,handle)
    NULLIFY(admm_env)

    CALL get_qs_env(qs_env, admm_env=admm_env, error=error)

    SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_cauchy)
        CALL merge_ks_matrix_cauchy(qs_env, error)

      CASE(do_admm_purify_cauchy_subspace)
        CALL merge_ks_matrix_cauchy_subspace(qs_env, error)

      CASE(do_admm_purify_none)
         CALL merge_ks_matrix_none(qs_env, error)

      CASE(do_admm_purify_mo_diag, do_admm_purify_mo_no_diag)
         !do nothing
      CASE DEFAULT
         STOP "admm_mo_merge_ks_matrix: unkown purification method"
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_mo_merge_ks_matrix

! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff ...
!> \param mo_coeff_aux_fit ...
!> \param mo_derivs ...
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE admm_mo_merge_derivs(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_mo_merge_derivs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_mo_diag)
        CALL merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit,&
                                  mo_derivs,mo_derivs_aux_fit, matrix_ks_aux_fit, error)

      CASE(do_admm_purify_mo_no_diag)
        CALL merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit,&
                                     mo_derivs,mo_derivs_aux_fit, matrix_ks_aux_fit, error)

      CASE(do_admm_purify_none, do_admm_purify_cauchy,do_admm_purify_cauchy_subspace)
          !do nothing
      CASE DEFAULT
          STOP "admm_mo_merge_derivs: unkown purification method"
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_mo_merge_derivs


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param para_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param geometry_did_change ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE admm_fit_mo_coeffs(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change, error)

    TYPE(admm_type), POINTER                 :: admm_env
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

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL fit_mo_coeffs(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                 mos, mos_aux_fit, geometry_did_change,&
                                 blocked=admm_env%block_fit, error=error)

    SELECT CASE(admm_env%purification_method)
       CASE(do_admm_purify_mo_no_diag, do_admm_purify_cauchy_subspace)
         CALL purify_mo_cholesky(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                 mos, mos_aux_fit, error)

       CASE(do_admm_purify_mo_diag)
         CALL purify_mo_diag(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                             mos, mos_aux_fit, error)

       CASE DEFAULT
         CALL purify_mo_none(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                             mos, mos_aux_fit, error)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_fit_mo_coeffs


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param para_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param geometry_did_change ...
!> \param blocked ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fit_mo_coeffs(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                   mos, mos_aux_fit, geometry_did_change, blocked, error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change, blocked
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, jatom, &
                                                nao_aux_fit, nao_orb, nspins
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s_tilde

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      IF(.NOT. blocked) THEN
         CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv,error)
      ELSE
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
      ENDIF

      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)

      !! Calculate A=S'^(-1)*Q
      IF(blocked) THEN
        CALL cp_fm_set_all(admm_env%A, 0.0_dp, 1.0_dp, error)
      ELSE
        CALL cp_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                      1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                      admm_env%A,error)

         ! this multiplication is apparent not need for purify_none
         !! B=Q^(T)*A
         CALL cp_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                         1.0_dp,admm_env%Q,admm_env%A,0.0_dp,&
                         admm_env%B,error)
      ENDIF
    END IF

    CALL timestop(handle)

  END SUBROUTINE fit_mo_coeffs


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
!> \param error ...
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE purify_mo_cholesky(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, error)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_cholesky', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

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



      CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
      CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin)%matrix,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,error=error)
      CALL cp_fm_to_fm(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,error=error)


      !! ** C_hat = AC
      CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_orb,&
                      1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                      admm_env%C_hat(ispin)%matrix,error)
      CALL cp_fm_to_fm(admm_env%C_hat(ispin)%matrix, mo_coeff_aux_fit, error=error)

    END DO

    CALL timestop(handle)

  END SUBROUTINE purify_mo_cholesky


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
!> \param error ...
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE purify_mo_diag(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, error)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, &
                                                nao_aux_fit, nao_orb, nmo, &
                                                nspins
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig_work
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

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
      ALLOCATE(eig_work(nmo))
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

      CALL cp_fm_to_fm(mo_coeff_aux_fit, admm_env%C_hat(ispin)%matrix, error=error)
      CALL cp_fm_set_all(admm_env%lambda_inv(ispin)%matrix,0.0_dp,1.0_dp,error)
      DEALLOCATE(eig_work)
    END DO

    CALL timestop(handle)

  END SUBROUTINE purify_mo_diag


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param para_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE purify_mo_none(admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo, nmo_mos, nspins
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occ_num, occ_num_aux
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

    DO ispin = 1,nspins
      nmo = admm_env%nmo(ispin)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, occupation_numbers=occ_num,nmo=nmo_mos)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,&
           occupation_numbers=occ_num_aux)

      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                      mo_coeff_aux_fit,error)
      CALL cp_fm_to_fm(mo_coeff_aux_fit, admm_env%C_hat(ispin)%matrix, error=error)

      occ_num_aux(1:nmo) = occ_num(1:nmo)
      ! XXXX should only be done first time XXXX
      CALL cp_fm_set_all(admm_env%lambda(ispin)%matrix,0.0_dp,1.0_dp,error)
      CALL cp_fm_set_all(admm_env%lambda_inv(ispin)%matrix,0.0_dp,1.0_dp,error)
      CALL cp_fm_set_all(admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,1.0_dp,error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE purify_mo_none


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param mo_set ...
!> \param density_matrix ...
!> \param ispin ...
!> \param blocked ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE purify_dm_cauchy(admm_env,mo_set,density_matrix,ispin,blocked,error)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    LOGICAL, INTENT(IN)                      :: blocked
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'purify_dm_cauchy', &
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

    IF(.NOT. blocked) THEN
       CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nmo,&
                     1.0_dp,mo_coeff_aux_fit,mo_coeff_aux_fit,0.0_dp,&
                     admm_env%P_to_be_purified(ispin)%matrix,error)
    ENDIF

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

  END SUBROUTINE purify_dm_cauchy


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE merge_ks_matrix_cauchy(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_cauchy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iatom, ispin, &
                                                j, jatom, nao_aux_fit, &
                                                nao_orb, nmo
    REAL(dp)                                 :: eig_diff, pole, tmp
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)
    NULLIFY(admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, mos, mo_coeff)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    mos=mos,&
                    error=error)

    DO ispin=1, dft_control%nspins
       nao_aux_fit = admm_env%nao_aux_fit
       nao_orb = admm_env%nao_orb
       nmo = admm_env%nmo(ispin)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)

       IF(.NOT. admm_env%block_dm) THEN
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

       ENDIF

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
            cp_dbcsr_distribution(matrix_ks(ispin)%matrix), &
            dbcsr_type_symmetric,&
            cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
            cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix),&
            cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix),&
            cp_dbcsr_get_data_size( matrix_ks(ispin)%matrix),&
            cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)

       CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error)

       CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
       CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
       CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE., error=error)

       IF(admm_env%block_dm) THEN
          ! ** now loop through the list and nullify blocks
          CALL cp_dbcsr_iterator_start(iter, matrix_k_tilde)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
            CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
            IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
               sparse_block = 0.0_dp
            END IF
          END DO
          CALL cp_dbcsr_iterator_stop(iter)
       ENDIF

       CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

       CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    END DO !spin-loop

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE merge_ks_matrix_cauchy_subspace(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'merge_ks_matrix_cauchy_subspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit

    CALL timeset(routineN,handle)
    NULLIFY(admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, mos, mos_aux_fit,&
            mo_coeff, mo_coeff_aux_fit)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    mos=mos,&
                    mos_aux_fit=mos_aux_fit,&
                    error=error)

    DO ispin=1, dft_control%nspins
       nao_aux_fit = admm_env%nao_aux_fit
       nao_orb = admm_env%nao_orb
       nmo = admm_env%nmo(ispin)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL get_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)

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
            cp_dbcsr_distribution(matrix_ks(ispin)%matrix),&
            dbcsr_type_symmetric,&
            cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
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

    ENDDO !spin loop
    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy_subspace


! *****************************************************************************
!> \brief Calculates the product Kohn-Sham-Matrix x mo_coeff for the auxiliary
!>        basis set and transforms it into the orbital basis. This is needed
!>        in order to use OT
!>
!> \param ispin which spin to transform
!> \param admm_env The ADMM env
!> \param mo_set ...
!> \param mo_coeff the MO coefficients from the orbital basis set
!> \param mo_coeff_aux_fit the MO coefficients from the auxiliary fitting basis set
!> \param mo_derivs KS x mo_coeff from the orbital basis set to which we add the
!>        auxiliary basis set part
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit the Kohn-Sham matrix from the auxiliary fitting basis set
!> \param error ...
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
                                                nao_orb, nmo
    REAL(dp)                                 :: eig_diff, pole, tmp32, tmp52, &
                                                tmp72, tmp92
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL cp_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%H(ispin)%matrix,error)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
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

    DEALLOCATE(scaling_factor)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs_diag


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE merge_ks_matrix_none(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, ispin, &
                                                jatom, nao_aux_fit, nao_orb, &
                                                nmo
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)
    NULLIFY(admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, matrix_k_tilde, sparse_block)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    !mos=mos,&
                    !mos_aux_fit=mos_aux_fit,&
                   error=error)

    DO ispin=1, dft_control%nspins
       IF(admm_env%block_dm) THEN
          CALL cp_dbcsr_iterator_start(iter, matrix_ks_aux_fit(ispin)%matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
            CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
            IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
              sparse_block = 0.0_dp
            END IF
          END DO
          CALL cp_dbcsr_iterator_stop(iter)
          CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, 1.0_dp, 1.0_dp, error)

       ELSE

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
          CALL cp_dbcsr_create(matrix_k_tilde,template=matrix_ks(ispin)%matrix,&
                               name='MATRIX K_tilde', matrix_type='S', error=error)
          CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
          CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
          CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.,&
               error=error)

          CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

          CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)
       ENDIF
    ENDDO !spin loop
    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_none


! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff ...
!> \param mo_coeff_aux_fit ...
!> \param matrix_ks_aux_fit ...
!> \param error ...
! *****************************************************************************
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
                                                nmo
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    CALL timeset(routineN,handle)

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
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
    scaling_factor = 2.0_dp*occupation_numbers
    CALL cp_fm_column_scale(admm_env%H(ispin)%matrix,scaling_factor)

    DEALLOCATE(scaling_factor)

    CALL timestop(handle)

  END SUBROUTINE calc_aux_mo_derivs_none


! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff ...
!> \param mo_coeff_aux_fit ...
!> \param mo_derivs ...
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit ...
!> \param error ...
! *****************************************************************************
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
                                                nmo
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
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

    DEALLOCATE(scaling_factor)

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
!> \param error ...
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_mixed_overlap_force(qs_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mixed_overlap_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw, &
                                                nao_aux_fit, nao_orb, natom, &
                                                neighbor_list_id, nmo
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: admm_force
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_dbcsr_type), POINTER             :: matrix_w_q, matrix_w_s
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux_fit, &
                                                sab_aux_fit_asymm, &
                                                sab_aux_fit_vs_orb
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(admm_env, logger, dft_control, para_env, mos, mo_coeff, matrix_w_q, matrix_w_s,&
            sab_aux_fit, sab_aux_fit_asymm, &
            sab_aux_fit_vs_orb, ks_env, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    neighbor_list_id=neighbor_list_id,&
                    sab_aux_fit=sab_aux_fit,&
                    sab_aux_fit_asymm=sab_aux_fit_asymm,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    mos=mos,&
                    para_env=para_env,&
                    error=error)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb

    logger => cp_error_get_logger(error)

    ! *** forces are only implemented for mo_diag or none and basis_projection ***
    IF (admm_env%block_dm) THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    IF (.NOT.(admm_env%purification_method==do_admm_purify_mo_diag .OR. &
               admm_env%purification_method==do_admm_purify_none)) THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ! *** Create sparse work matrices

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
    CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s,sab_aux_fit_asymm,error=error)

    ALLOCATE(matrix_w_q)
    CALL cp_dbcsr_init(matrix_w_q, error=error)
    CALL cp_dbcsr_copy(matrix_w_q,matrix_s_aux_fit_vs_orb(1)%matrix,&
                    "W MATRIX AUX Q",error=error)

    DO ispin=1, dft_control%nspins
        nmo = admm_env%nmo(ispin)
        CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)

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

        ! allocate force vector
        CALL get_qs_env(qs_env=qs_env,natom=natom,error=error)
        ALLOCATE(admm_force(3,natom))
        admm_force = 0.0_dp
        CALL build_overlap_force(ks_env, admm_force,&
             basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_aux_fit_basis_set, &
             sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s, error=error)
        CALL build_overlap_force(ks_env, admm_force,&
             basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_orb_basis_set, &
             sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q, error=error)
        ! add forces
        CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                        force=force,error=error)
        CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set, error)
        DEALLOCATE(admm_force)

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

    ENDDO !spin loop

    ! *** Deallocated weighted density matrices
    CALL cp_dbcsr_deallocate_matrix(matrix_w_s,error)
    CALL cp_dbcsr_deallocate_matrix(matrix_w_q,error)

    CALL timestop(handle)

  END SUBROUTINE calc_mixed_overlap_force


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param mo_set ...
!> \param density_matrix ...
!> \param ispin ...
!> \param error ...
! *****************************************************************************
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
      !alpha=mo_set%maxocc
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE calculate_dm_mo_no_diag

! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param density_matrix ...
!> \param density_matrix_aux ...
!> \param ispin ...
!> \param nspins ...
!> \param error ...
! *****************************************************************************
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
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block, sparse_block_aux
    TYPE(cp_dbcsr_iterator)                  :: iter

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

! *****************************************************************************
!> \brief ...
!> \param x ...
! *****************************************************************************
  FUNCTION delta(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: delta

    IF( x == 0.0_dp) THEN !TODO: exact comparison of reals?
      delta = 1.0_dp
    ELSE
      delta = 0.0_dp
    END IF

  END FUNCTION delta

! *****************************************************************************
!> \brief ...
!> \param x ...
! *****************************************************************************
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
