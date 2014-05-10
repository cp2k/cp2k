!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains methods used in the context of density fitting
!> \par History
!>      04.2008 created [Manuel Guidon]
!>      02.2013 moved from admm_methods
!> \author Manuel Guidon
! *****************************************************************************
MODULE admm_utils
  USE admm_types,                      ONLY: admm_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_distribution, &
       cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, &
       cp_dbcsr_get_num_blocks, cp_dbcsr_init, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_type, dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE input_constants,                 ONLY: do_admm_basis_set_projection,&
                                             do_admm_block_density_matrix,&
                                             do_admm_purify_cauchy,&
                                             do_admm_purify_cauchy_subspace,&
                                             do_admm_purify_mo_diag,&
                                             do_admm_purify_mo_no_diag,&
                                             do_admm_purify_none
  USE kinds,                           ONLY: dp
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC admm_correct_for_eigenvalues,&
         admm_uncorrect_for_eigenvalues

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_utils'

!***

  CONTAINS

  SUBROUTINE admm_correct_for_eigenvalues(ispin, admm_env, ks_matrix, ks_matrix_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_type), POINTER             :: ks_matrix, ks_matrix_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: nao_aux_fit, nao_orb
    TYPE(cp_dbcsr_type), POINTER             :: work

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb


    SELECT CASE(admm_env%method_id)
    CASE (do_admm_basis_set_projection)
      SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_none)
      CASE(do_admm_purify_cauchy)
      CASE(do_admm_purify_cauchy_subspace)
        !* remove what has been added and add the correction
        NULLIFY(work)
        ALLOCATE(work)
        CALL cp_dbcsr_init (work, error)
        CALL cp_dbcsr_create(work, 'work', &
             cp_dbcsr_distribution(ks_matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(ks_matrix),&
             cp_dbcsr_col_block_sizes(ks_matrix), cp_dbcsr_get_num_blocks(ks_matrix), &
             cp_dbcsr_get_data_size(ks_matrix),&
             cp_dbcsr_get_data_type(ks_matrix), error=error)

        CALL cp_dbcsr_copy(work, ks_matrix, error=error)
        CALL cp_dbcsr_set(work, 0.0_dp, error)
        CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp, error)

        ! ** calculate A^T*H_tilde*A
        CALL cp_fm_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                        1.0_dp, admm_env%K(ispin)%matrix,admm_env%A, 0.0_dp,&
                        admm_env%work_aux_orb,error)
        CALL cp_fm_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                        1.0_dp, admm_env%A, admm_env%work_aux_orb ,0.0_dp,&
                        admm_env%H_corr(ispin)%matrix,error)

        CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp, error)
        CALL cp_dbcsr_deallocate_matrix(work,error)

      CASE(do_admm_purify_mo_diag)
        !* remove what has been added and add the correction
        NULLIFY(work)
        ALLOCATE(work)
        CALL cp_dbcsr_init (work, error)
        CALL cp_dbcsr_create(work, 'work', &
             cp_dbcsr_distribution(ks_matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(ks_matrix),&
             cp_dbcsr_col_block_sizes(ks_matrix), cp_dbcsr_get_num_blocks(ks_matrix), cp_dbcsr_get_data_size(ks_matrix),&
             cp_dbcsr_get_data_type(ks_matrix), error=error)

        CALL cp_dbcsr_copy(work, ks_matrix, error=error)
        CALL cp_dbcsr_set(work, 0.0_dp, error)
        CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        ! ** calculate A^T*H_tilde*A
        CALL cp_fm_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                        1.0_dp, admm_env%K(ispin)%matrix,admm_env%A, 0.0_dp,&
                        admm_env%work_aux_orb,error)
        CALL cp_fm_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                        1.0_dp, admm_env%A, admm_env%work_aux_orb ,0.0_dp,&
                        admm_env%H_corr(ispin)%matrix,error)

        CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp, error)
        CALL cp_dbcsr_deallocate_matrix(work,error)

      CASE(do_admm_purify_mo_no_diag)
    END SELECT

    CASE(do_admm_block_density_matrix)

    END SELECT


  END SUBROUTINE admm_correct_for_eigenvalues

  SUBROUTINE admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_matrix, ks_matrix_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_type), POINTER             :: ks_matrix, ks_matrix_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: nao_aux_fit, nao_orb
    TYPE(cp_dbcsr_type), POINTER             :: work

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb


    SELECT CASE(admm_env%method_id)
    CASE (do_admm_basis_set_projection)
      SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_none)
      CASE(do_admm_purify_cauchy)
      CASE(do_admm_purify_cauchy_subspace)
        !* remove what has been added and add the correction
        NULLIFY(work)
        ALLOCATE(work)
        CALL cp_dbcsr_init (work, error)
        CALL cp_dbcsr_create(work, 'work', &
             cp_dbcsr_distribution(ks_matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(ks_matrix),&
             cp_dbcsr_col_block_sizes(ks_matrix), cp_dbcsr_get_num_blocks(ks_matrix), cp_dbcsr_get_data_size(ks_matrix),&
             cp_dbcsr_get_data_type(ks_matrix), error=error)

        CALL cp_dbcsr_copy(work, ks_matrix, error=error)
        CALL cp_dbcsr_set(work, 0.0_dp, error)
        CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp, error)

        CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_set(work, 0.0_dp, error)
        CALL copy_fm_to_dbcsr(admm_env%ks_to_be_merged(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_add(ks_matrix, work, 1.0_dp, 1.0_dp, error)
        CALL cp_dbcsr_deallocate_matrix(work,error)

      CASE(do_admm_purify_mo_diag)
        NULLIFY(work)
        ALLOCATE(work)
        CALL cp_dbcsr_init (work, error)
        CALL cp_dbcsr_create(work, 'work', &
             cp_dbcsr_distribution(ks_matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(ks_matrix),&
             cp_dbcsr_col_block_sizes(ks_matrix), cp_dbcsr_get_num_blocks(ks_matrix), cp_dbcsr_get_data_size(ks_matrix),&
             cp_dbcsr_get_data_type(ks_matrix), error=error)

        CALL cp_dbcsr_copy(work, ks_matrix, error=error)
        CALL cp_dbcsr_set(work, 0.0_dp, error)

        CALL copy_fm_to_dbcsr(admm_env%H_corr(ispin)%matrix, work, keep_sparsity=.TRUE.,&
             error=error)

        CALL cp_dbcsr_add(ks_matrix, work, 1.0_dp, -1.0_dp, error)
        CALL cp_dbcsr_deallocate_matrix(work,error)

      CASE(do_admm_purify_mo_no_diag)
    END SELECT

    CASE(do_admm_block_density_matrix)

    END SELECT

  END SUBROUTINE admm_uncorrect_for_eigenvalues

END MODULE admm_utils
