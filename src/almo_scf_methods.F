!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Subroutines for ALMO SCF
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_methods
  USE almo_scf_types,                  ONLY: almo_scf_env_type
  USE cp_dbcsr_interface,              ONLY: &
       array_data, cp_dbcsr_add, cp_dbcsr_add_on_diag, &
       cp_dbcsr_col_block_sizes, cp_dbcsr_copy, cp_dbcsr_create, &
       cp_dbcsr_desymmetrize, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_frobenius_norm, cp_dbcsr_get_block_p, &
       cp_dbcsr_get_diag, cp_dbcsr_get_info, cp_dbcsr_init, &
       cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
       cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, &
       cp_dbcsr_iterator_stop, cp_dbcsr_multiply, cp_dbcsr_nblkcols_total, &
       cp_dbcsr_nblkrows_total, cp_dbcsr_norm, cp_dbcsr_print, &
       cp_dbcsr_release, cp_dbcsr_reserve_block2d, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_set_diag, cp_dbcsr_transposed, &
       cp_dbcsr_type, cp_dbcsr_work_create, dbcsr_distribution_mp, &
       dbcsr_error_type, dbcsr_init_random, dbcsr_mp_group, &
       dbcsr_mp_numnodes, dbcsr_norm_maxabsnorm, dbcsr_type_no_symmetry, &
       dbcsr_type_symmetric
  USE domain_submatrix_methods,        ONLY: &
       add_submatrices, construct_dbcsr_from_submatrices, &
       construct_submatrices, copy_submatrices, copy_submatrix_data, &
       init_submatrices, multiply_submatrices, print_submatrices, &
       release_submatrices
  USE domain_submatrix_types,          ONLY: domain_map_type,&
                                             domain_submatrix_type,&
                                             select_row,&
                                             select_row_col
  USE f77_blas
  USE input_constants,                 ONLY: almo_domain_layout_molecular,&
                                             almo_mat_distr_atomic,&
                                             almo_scf_diag
  USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                             matrix_sqrt_Newton_Schulz
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_methods'

  PUBLIC almo_scf_ks_to_ks_blk, almo_scf_p_blk_to_t_blk,&
         almo_scf_t_blk_to_p, almo_scf_t_blk_to_t_blk_orthonormal,&
         almo_scf_t_to_p, almo_scf_ks_blk_to_tv_blk,&
         almo_scf_ks_xx_to_tv_xx,&
         apply_projector, get_overlap,&
         generator_to_unitary,&
         newton_grad_to_step, orthogonalize_mos, &
         pseudo_invert_diagonal_blk, construct_test,&
         construct_domain_preconditioner,&
         apply_domain_operators,&
         construct_domain_s_inv,&
         construct_domain_s_sqrt,&
         distribute_domains,&
         almo_scf_ks_to_ks_xx,&
         construct_domain_r_down

CONTAINS

! *****************************************************************************
!> \brief builds projected KS matrices for the overlapping domains
!>        also computes the DIIS error vector as a by-product
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_ks_to_ks_xx(almo_scf_env,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_to_ks_xx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, ndomains
    REAL(KIND=dp)                            :: eps_multiply
    TYPE(cp_dbcsr_type) :: matrix_tmp1, matrix_tmp2, matrix_tmp3, &
      matrix_tmp4, matrix_tmp5, matrix_tmp6, matrix_tmp7, matrix_tmp8, &
      matrix_tmp9
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_tmp1, subm_tmp2, &
                                                subm_tmp3

    CALL timeset(routineN,handle)

    eps_multiply=almo_scf_env%eps_filter

    DO ispin=1,almo_scf_env%nspins

       ndomains = cp_dbcsr_nblkcols_total(almo_scf_env%quench_t(ispin))

       ! 0. Create KS_xx 
       CALL construct_submatrices(&
               almo_scf_env%matrix_ks(ispin),&
               almo_scf_env%domain_ks_xx(:,ispin),&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row_col,&
               error)

       !!!!! RZK-warning MAKE SURE THAT YOU NEED BLOCKS OUTSIDE QUENCH_T
       !!!!! FOR ALL NO-MATRICES NOT COMPUTING THEM CAN SAVE LOTS OF TIME

       ! 1. TMP1=KS.T
       !    Cost: NOn
       !matrix_tmp1 = create NxO, full
       CALL cp_dbcsr_init(matrix_tmp1, error=error)
       CALL cp_dbcsr_create(matrix_tmp1,&
                            template=almo_scf_env%matrix_t(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_ks(ispin),&
                                     almo_scf_env%matrix_t(ispin),&
                                     0.0_dp, matrix_tmp1,&
                                     filter_eps=eps_multiply,&
                                     error=error)

       ! 2. TMP2=TMP1.SigInv=KS.T.SigInv
       !    Cost: NOO
       !matrix_tmp2 = create NxO, full
       CALL cp_dbcsr_init(matrix_tmp2, error=error)
       CALL cp_dbcsr_create(matrix_tmp2,&
                            template=almo_scf_env%matrix_t(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1,&
                                     almo_scf_env%matrix_sigma_inv(ispin),&
                                     0.0_dp, matrix_tmp2,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       ! 3. TMP1=S.T
       !    Cost: NOn
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1),&
               almo_scf_env%matrix_t(ispin),&
               0.0_dp, matrix_tmp1,&
               filter_eps=eps_multiply,&
               error=error)

       ! 4. TMP4=TMP2.tr(TMP1)=KS.T.SigInv.tr(T).S
       !    Cost: NNO
       !matrix_tmp4 = create NxN
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
               template=almo_scf_env%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp2,&
               matrix_tmp1,&
               0.0_dp, matrix_tmp4,&
               filter_eps=eps_multiply,&
               error=error)
       
       ! 5. KS_xx=KS_xx-TMP4_xx-tr(TMP4_xx)
       ALLOCATE(subm_tmp1(ndomains))
       CALL init_submatrices(subm_tmp1,error)
       CALL construct_submatrices(&
               matrix_tmp4,&
               subm_tmp1,&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row_col,&
               error)
       CALL add_submatrices(1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               -1.0_dp,subm_tmp1,'N',error)
       CALL add_submatrices(1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               -1.0_dp,subm_tmp1,'T',error)

       ! 6. TMP3=tr(TMP4).T=S.T.SigInv.tr(T).KS.T 
       !    Cost: NOn
       !matrix_tmp3 = create NxO, full
       CALL cp_dbcsr_init(matrix_tmp3, error=error)
       CALL cp_dbcsr_create(matrix_tmp3,&
                            template=almo_scf_env%matrix_t(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                                     matrix_tmp4,&
                                     almo_scf_env%matrix_t(ispin),&
                                     0.0_dp, matrix_tmp3,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_release(matrix_tmp4,error=error)
       
       ! 8. TMP6=TMP3.SigInv=S.T.SigInv.tr(T).KS.T.SigInv 
       !    Cost: NOO
       !matrix_tmp6 = create NxO, full
       CALL cp_dbcsr_init(matrix_tmp6, error=error)
       CALL cp_dbcsr_create(matrix_tmp6,&
                            template=almo_scf_env%matrix_t(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              matrix_tmp3,&
                              almo_scf_env%matrix_sigma_inv(ispin),&
                              0.0_dp, matrix_tmp6,&
                              filter_eps=eps_multiply,&
                              error=error)
       
       ! 8A. Use intermediate matrices to evaluate the gradient/error 
       !     Err=(TMP2-TMP6)_q=(KS.T.SigInv-S.T.SigInv.tr(T).KS.T.SigInv)_q
       ! error vector in AO-MO basis
       CALL cp_dbcsr_copy(almo_scf_env%matrix_err_xx(ispin),&
               almo_scf_env%quench_t(ispin),error=error)
       CALL cp_dbcsr_copy(almo_scf_env%matrix_err_xx(ispin),&
               matrix_tmp2,keep_sparsity=.TRUE.,error=error)
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
               template=almo_scf_env%matrix_t(ispin),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_copy(matrix_tmp4,&
               almo_scf_env%quench_t(ispin),error=error)
       CALL cp_dbcsr_copy(matrix_tmp4,&
               matrix_tmp6,keep_sparsity=.TRUE.,error=error)
       CALL cp_dbcsr_add(almo_scf_env%matrix_err_xx(ispin),&
               matrix_tmp4,1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_release(matrix_tmp4, error=error)
       !
       ! error vector in AO-AO basis
       ! RZK-warning tmp4 can be created using the sparsity pattern, 
       ! then retain_sparsity can be used to perform the multiply
       ! this will save some time
       CALL cp_dbcsr_copy(matrix_tmp3,&
               matrix_tmp2,error=error)
       CALL cp_dbcsr_add(matrix_tmp3,&
               matrix_tmp6,1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
               template=almo_scf_env%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp,&
               matrix_tmp3,&
               almo_scf_env%matrix_t(ispin),&
               0.0_dp, matrix_tmp4,&
               filter_eps=eps_multiply,&
               error=error)
       CALL construct_submatrices(&
               matrix_tmp4,&
               almo_scf_env%domain_err(:,ispin),&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row_col,&
               error)
       CALL cp_dbcsr_release(matrix_tmp4, error=error)
       ! domain_err submatrices are in down-up representation
       ! bring them into the orthogonalized basis
       ALLOCATE(subm_tmp2(ndomains))
       CALL init_submatrices(subm_tmp2,error)
       CALL multiply_submatrices('N','N',1.0_dp,&
               almo_scf_env%domain_err(:,ispin),&
               almo_scf_env%domain_s_sqrt(:,ispin),0.0_dp,subm_tmp2,error)
       CALL multiply_submatrices('N','N',1.0_dp,&
               almo_scf_env%domain_s_sqrt_inv(:,ispin),&
               subm_tmp2,0.0_dp,almo_scf_env%domain_err(:,ispin),error)
       
       ! 9. TMP5=TMP6.tr(TMP1)=S.T.SigInv.tr(T).KS.T.SigInv.tr(T).S 
       !    Cost: NNO
       !    matrix_tmp5 = create NxN, full
       ! RZK-warning tmp5 can be created using the sparsity pattern, 
       ! then retain_sparsity can be used to perform the multiply
       ! this will save some time
       CALL cp_dbcsr_init(matrix_tmp5, error=error)
       CALL cp_dbcsr_create(matrix_tmp5,&
                            template=almo_scf_env%matrix_s(1),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp,&
                              matrix_tmp6,&
                              matrix_tmp1,&
                              0.0_dp, matrix_tmp5,&
                              filter_eps=eps_multiply,&
                              error=error)

       ! 10. KS_xx=KS_xx+TMP5_xx
       CALL construct_submatrices(&
               matrix_tmp5,&
               subm_tmp1,&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row_col,&
               error)
       CALL cp_dbcsr_release(matrix_tmp5,error=error)
       CALL add_submatrices(1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               1.0_dp,subm_tmp1,'N',error)

       ! 11. KS_xx=KS_xx + [S.T]_xx.[SigInv.tr(T).KS.(1-T.SigInv.tr(T).S)]_xx + transposed 
       ALLOCATE(subm_tmp3(ndomains))
       CALL init_submatrices(subm_tmp3,error)
       CALL construct_submatrices(&
               matrix_tmp2,&
               subm_tmp2,&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row,&
               error)
       CALL construct_submatrices(&
               matrix_tmp6,&
               subm_tmp3,&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row,&
               error)
       CALL cp_dbcsr_release(matrix_tmp6,error=error)
       CALL add_submatrices(1.0_dp,subm_tmp2,&
               -1.0_dp,subm_tmp3,'N',error)
       CALL construct_submatrices(&
               matrix_tmp1,&
               subm_tmp3,&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row,&
               error)
       CALL multiply_submatrices('N','T',1.0_dp,subm_tmp2,&
               subm_tmp3,0.0_dp,subm_tmp1,error)
       CALL add_submatrices(1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               1.0_dp,subm_tmp1,'N',error)
       CALL add_submatrices(1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               1.0_dp,subm_tmp1,'T',error)
      
       ! 12. TMP7=tr(T).KS.T.SigInv 
       CALL cp_dbcsr_init(matrix_tmp7, error=error)
       CALL cp_dbcsr_create(matrix_tmp7,&
                            template=almo_scf_env%matrix_sigma_blk(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                              almo_scf_env%matrix_t(ispin),&
                              matrix_tmp2,&
                              0.0_dp, matrix_tmp7,&
                              filter_eps=eps_multiply,&
                              error=error)

       ! 13. TMP8=[SigInv.tr(T).KS.T.SigInv]_xx 
       CALL cp_dbcsr_init(matrix_tmp8, error=error)
       CALL cp_dbcsr_create(matrix_tmp8,&
                            template=almo_scf_env%matrix_sigma_blk(ispin),&
                            matrix_type=dbcsr_type_symmetric,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp8,almo_scf_env%matrix_sigma_blk(ispin),&
                          error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              almo_scf_env%matrix_sigma_inv(ispin),&
                              matrix_tmp7,&
                              0.0_dp, matrix_tmp8,&
                              retain_sparsity=.TRUE.,&
                              filter_eps=eps_multiply,&
                              error=error)
       CALL cp_dbcsr_release(matrix_tmp7,error=error)

       ! 13. TMP9=[S.T]_xx 
       CALL cp_dbcsr_init(matrix_tmp9, error=error)
       CALL cp_dbcsr_create(matrix_tmp9,&
                            template=almo_scf_env%matrix_t(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,almo_scf_env%quench_t(ispin),error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,matrix_tmp1,keep_sparsity=.TRUE.,&
                          error=error)

       ! 14. TMP3=TMP9.TMP8=[S.T]_xx.[SigInv.tr(T).KS.T.SigInv]_xx
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              matrix_tmp9,&
                              matrix_tmp8,&
                              0.0_dp, matrix_tmp3,&
                              filter_eps=eps_multiply,&
                              error=error)
       CALL cp_dbcsr_release(matrix_tmp8,error=error)
       CALL cp_dbcsr_release(matrix_tmp9,error=error)

       ! 15. KS_xx=KS_xx+[S.T]_xx.[SigInv.tr(T).KS.T.SigInv]_xx.[tr(T).S]_xx
       CALL construct_submatrices(&
               matrix_tmp3,&
               subm_tmp2,&
               almo_scf_env%quench_t(ispin),&
               almo_scf_env%domain_map(ispin),&
               almo_scf_env%cpu_of_domain,&
               select_row,&
               error)
       CALL multiply_submatrices('N','T',1.0_dp,subm_tmp2,&
               subm_tmp3,0.0_dp,subm_tmp1,error)
       CALL add_submatrices(1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               1.0_dp,subm_tmp1,'N',error)
       
       !!!!!!! use intermediate matrices to get the error vector !!!!!!!
       !!!!!!! make sure s_blk_sqrt and its inverse exist (i.e. we use diag algorithm)
       !CPPrecondition(almo_scf_env%almo_update_algorithm.eq.almo_scf_diag,cp_failure_level,routineP,error,failure)
       !! tmp_err = (1-S.T_blk.SigInv.tr(T_blk)).F.T_blk.SigInv
       !CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       !CALL cp_dbcsr_create(matrix_tmp_err,&
       !        template=almo_scf_env%matrix_t(ispin),&
       !        error=error)
       !CALL cp_dbcsr_copy(matrix_tmp_err,&
       !        matrix_tmp2,&
       !        error=error)
       !CALL cp_dbcsr_add(matrix_tmp_err,matrix_tmp3,&
       !        1.0_dp,-1.0_dp,error=error)
       !! err_blk = tmp_err.tr(T_blk)
       !CALL cp_dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
       !        almo_scf_env%matrix_s_blk_sqrt(1),&
       !        error=error)
       !CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp_err,&
       !        almo_scf_env%matrix_t(ispin),&
       !        0.0_dp, almo_scf_env%matrix_err_blk(ispin),&
       !        retain_sparsity=.TRUE.,&
       !        filter_eps=eps_multiply,&
       !        error=error)
       !CALL cp_dbcsr_release(matrix_tmp_err,error=error)
       !! bring to the orthogonal basis
       !! err_blk = (S_blk^-1/2).err_blk.(S_blk^1/2)
       !CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       !CALL cp_dbcsr_create(matrix_tmp_err,&
       !        template=almo_scf_env%matrix_err_blk(ispin),&
       !        error=error)
       !CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
       !        almo_scf_env%matrix_err_blk(ispin),&
       !        almo_scf_env%matrix_s_blk_sqrt(1),&
       !        0.0_dp, matrix_tmp_err,&
       !        filter_eps=eps_multiply,&
       !        error=error)
       !CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
       !        almo_scf_env%matrix_s_blk_sqrt_inv(1),&
       !        matrix_tmp_err,&
       !        0.0_dp, almo_scf_env%matrix_err_blk(ispin),&
       !        filter_eps=eps_multiply,&
       !        error=error)
       !! subtract transpose
       !CALL cp_dbcsr_transposed(matrix_tmp_err,&
       !        almo_scf_env%matrix_err_blk(ispin),error=error)
       !CALL cp_dbcsr_add(almo_scf_env%matrix_err_blk(ispin),&
       !        matrix_tmp_err,&
       !        1.0_dp,-1.0_dp,error=error)
       !CALL cp_dbcsr_release(matrix_tmp_err,error=error)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       CALL release_submatrices(subm_tmp3,error)
       CALL release_submatrices(subm_tmp2,error)
       CALL release_submatrices(subm_tmp1,error)
       DEALLOCATE(subm_tmp3)
       DEALLOCATE(subm_tmp2)
       DEALLOCATE(subm_tmp1)
       CALL cp_dbcsr_release(matrix_tmp3,error=error)
       CALL cp_dbcsr_release(matrix_tmp2,error=error)
       CALL cp_dbcsr_release(matrix_tmp1,error=error)

    ENDDO ! spins

    CALL timestop(handle)

  END SUBROUTINE almo_scf_ks_to_ks_xx

! *****************************************************************************
!> \brief computes the projected KS from the total KS matrix
!>        also computes the DIIS error vector as a by-product
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_ks_to_ks_blk(almo_scf_env,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_to_ks_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps_multiply
    TYPE(cp_dbcsr_type) :: matrix_tmp1, matrix_tmp2, matrix_tmp3, &
      matrix_tmp4, matrix_tmp5, matrix_tmp6, matrix_tmp7, matrix_tmp8, &
      matrix_tmp9, matrix_tmp_err

    CALL timeset(routineN,handle)

    eps_multiply=almo_scf_env%eps_filter

    DO ispin=1,almo_scf_env%nspins

       ! 1. TMP1=KS.T_blk
       !    Cost: NOn
       !matrix_tmp1 = create NxO, full
       CALL cp_dbcsr_init(matrix_tmp1, error=error)
       CALL cp_dbcsr_create(matrix_tmp1,&
                            template=almo_scf_env%matrix_t(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_ks(ispin),&
                                     almo_scf_env%matrix_t_blk(ispin),&
                                     0.0_dp, matrix_tmp1,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       ! 2. TMP2=TMP1.SigInv=KS.T_blk.SigInv
       !    Cost: NOO
       !matrix_tmp2 = create NxO, full
       CALL cp_dbcsr_init(matrix_tmp2, error=error)
       CALL cp_dbcsr_create(matrix_tmp2,&
                            template=almo_scf_env%matrix_t(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1,&
                                     almo_scf_env%matrix_sigma_inv(ispin),&
                                     0.0_dp, matrix_tmp2,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       !!!!!! use intermediate matrices to get the error vector !!!!!!!
       !CALL cp_dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
       !        almo_scf_env%matrix_t_blk(ispin),&
       !        error=error)
       !CALL cp_dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
       !        matrix_tmp2,&
       !        keep_sparsity=.TRUE.,&
       !        error=error)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! 3. TMP1=S.T_blk
       !    Cost: NOn
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1),&
                                     almo_scf_env%matrix_t_blk(ispin),&
                                     0.0_dp, matrix_tmp1,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       ! 4. TMP4_blk=TMP2.tr(TMP1)=KS.T_blk.SigInv.tr(T_blk).S
       !    Cost: NnO
       !matrix_tmp4 = create NxN, blk
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
                            template=almo_scf_env%matrix_s_blk(1),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp4,almo_scf_env%matrix_s_blk(1),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp2,&
                                     matrix_tmp1,&
                                     0.0_dp, matrix_tmp4,&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       ! 5. KS_blk=KS_blk-TMP4_blk
       CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),&
                                     matrix_tmp4,&
                                     1.0_dp,-1.0_dp,error=error)
       
       ! 6. TMP5_blk=tr(TMP4_blk)
       !    KS_blk=KS_blk-tr(TMP4_blk)
       !matrix_tmp5 = create NxN, blk
       CALL cp_dbcsr_init(matrix_tmp5, error=error)
       CALL cp_dbcsr_create(matrix_tmp5,&
                            template=almo_scf_env%matrix_s_blk(1),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_transposed(matrix_tmp5,matrix_tmp4,error=error)
       CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),matrix_tmp5,&
                                     1.0_dp,-1.0_dp,error=error)
       
       ! 7. TMP3=tr(T_blk).TMP2=tr(T_blk).KS.T_blk.SigInv
       !    Cost: OOn
       !matrix_tmp3 = create OxO, full
       CALL cp_dbcsr_init(matrix_tmp3, error=error)
       CALL cp_dbcsr_create(matrix_tmp3,&
                            template=almo_scf_env%matrix_sigma_inv(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                                     almo_scf_env%matrix_t_blk(ispin),&
                                     matrix_tmp2,&
                                     0.0_dp, matrix_tmp3,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       ! 8. TMP6=SigInv.TMP3=SigInv.tr(T_blk).KS.T_blk.SigInv
       !    Cost: OOO
       !matrix_tmp6 = create OxO, full
       CALL cp_dbcsr_init(matrix_tmp6, error=error)
       CALL cp_dbcsr_create(matrix_tmp6,&
                            template=almo_scf_env%matrix_sigma_inv(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              almo_scf_env%matrix_sigma_inv(ispin),&
                              matrix_tmp3,&
                              0.0_dp, matrix_tmp6,&
                              filter_eps=eps_multiply,&
                              error=error)
       
       ! 9. TMP3=TMP1.TMP6=S.T_blk.SigInv.tr(T_blk).KS.T_blk.SigInv
       !    Cost: NOO
       !matrix_tmp3 = re-create NxO, full
       CALL cp_dbcsr_release(matrix_tmp3,error=error)
       CALL cp_dbcsr_init(matrix_tmp3, error=error)
       CALL cp_dbcsr_create(matrix_tmp3,&
                            template=almo_scf_env%matrix_t(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1,&
                                     matrix_tmp6,&
                                     0.0_dp, matrix_tmp3,&
                                     filter_eps=eps_multiply,&
                                     error=error)

       !!!!!! use intermediate matrices to get the error vector !!!!!!!
       !CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       !CALL cp_dbcsr_create(matrix_tmp_err,&
       !        template=almo_scf_env%matrix_t_blk(ispin),&
       !        error=error)
       !CALL cp_dbcsr_copy(matrix_tmp_err,&
       !        almo_scf_env%matrix_t_blk(ispin),&
       !        error=error)
       !CALL cp_dbcsr_copy(matrix_tmp_err,matrix_tmp3,&
       !        keep_sparsity=.TRUE.,&
       !        error=error)
       !CALL cp_dbcsr_add(almo_scf_env%matrix_err_blk(ispin),matrix_tmp_err,&
       !        1.0_dp,-1.0_dp,error=error)
       !CALL cp_dbcsr_release(matrix_tmp_err,error=error)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       !!!!!! use intermediate matrices to get the error vector !!!!!!!
       !!!!!! make sure s_blk_sqrt and its inverse exist (i.e. we use diag algorithm)
       CPPrecondition(almo_scf_env%almo_update_algorithm.eq.almo_scf_diag,cp_failure_level,routineP,error,failure)
       ! tmp_err = (1-S.T_blk.SigInv.tr(T_blk)).F.T_blk.SigInv
       CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       CALL cp_dbcsr_create(matrix_tmp_err,&
               template=almo_scf_env%matrix_t_blk(ispin),&
               error=error)
       CALL cp_dbcsr_copy(matrix_tmp_err,&
               matrix_tmp2,&
               error=error)
       CALL cp_dbcsr_add(matrix_tmp_err,matrix_tmp3,&
               1.0_dp,-1.0_dp,error=error)
       ! err_blk = tmp_err.tr(T_blk)
       CALL cp_dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
               almo_scf_env%matrix_s_blk_sqrt(1),&
               error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp_err,&
               almo_scf_env%matrix_t_blk(ispin),&
               0.0_dp, almo_scf_env%matrix_err_blk(ispin),&
               retain_sparsity=.TRUE.,&
               filter_eps=eps_multiply,&
               error=error)
       CALL cp_dbcsr_release(matrix_tmp_err,error=error)
       ! bring to the orthogonal basis
       ! err_blk = (S_blk^-1/2).err_blk.(S_blk^1/2)
       CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       CALL cp_dbcsr_create(matrix_tmp_err,&
               template=almo_scf_env%matrix_err_blk(ispin),&
               error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
               almo_scf_env%matrix_err_blk(ispin),&
               almo_scf_env%matrix_s_blk_sqrt(1),&
               0.0_dp, matrix_tmp_err,&
               filter_eps=eps_multiply,&
               error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
               almo_scf_env%matrix_s_blk_sqrt_inv(1),&
               matrix_tmp_err,&
               0.0_dp, almo_scf_env%matrix_err_blk(ispin),&
               filter_eps=eps_multiply,&
               error=error)
       ! subtract transpose
       CALL cp_dbcsr_transposed(matrix_tmp_err,&
               almo_scf_env%matrix_err_blk(ispin),error=error)
       CALL cp_dbcsr_add(almo_scf_env%matrix_err_blk(ispin),&
               matrix_tmp_err,&
               1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_release(matrix_tmp_err,error=error)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! later we will need only the blk version of TMP6
       ! create it here and release TMP6
       !matrix_tmp9 = create OxO, blk
       !matrix_tmp9 = copy data from matrix_tmp6, retain sparsity
       !matrix_tmp6 = release
       CALL cp_dbcsr_init(matrix_tmp9, error=error)
       CALL cp_dbcsr_create(matrix_tmp9,&
                            template=almo_scf_env%matrix_sigma_blk(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,almo_scf_env%matrix_sigma_blk(ispin),&
                          error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,matrix_tmp6,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp6,error=error)
       
       !10. KS_blk=KS_blk+TMP3.tr(TMP1)=
       !          =KS_blk+S.T_blk.SigInv.tr.(T_blk).KS.T_blk.SigInv.tr(T_blk).S
       !    Cost: NnO
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp3,&
                                     matrix_tmp1,&
                                     1.0_dp, almo_scf_env%matrix_ks_blk(ispin),&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       ! 11. TMP4_blk=TMP7_blk.tr(TMP8_blk)
       !    Cost: Nnn
       !matrix_tmp7 = create NxO, blk
       !matrix_tmp7 = copy data from matrix_tmp3, retain sparsity
       !matrix_tmp3 = release
       !matrix_tmp8 = create NxO, blk
       !matrix_tmp8 = copy data from matrix_tmp1, retain sparsity
       !matrix_tmp1 = release
       CALL cp_dbcsr_init(matrix_tmp7, error=error)
       CALL cp_dbcsr_create(matrix_tmp7,&
                            template=almo_scf_env%matrix_t_blk(ispin),&
                            error=error)
       ! transfer only the ALMO blocks from tmp3 into tmp7:
       ! first, copy t_blk into tmp7 to transfer the blk structure,
       ! then copy tmp3 into tmp7 with retain_sparsity 
       CALL cp_dbcsr_copy(matrix_tmp7,almo_scf_env%matrix_t_blk(ispin),&
                          error=error)
       CALL cp_dbcsr_copy(matrix_tmp7,matrix_tmp3,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp3,error=error)
       ! do the same for tmp1->tmp8
       CALL cp_dbcsr_init(matrix_tmp8, error=error)
       CALL cp_dbcsr_create(matrix_tmp8,&
                            template=almo_scf_env%matrix_t_blk(ispin),&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp8,almo_scf_env%matrix_t_blk(ispin),&
                          error=error)
       CALL cp_dbcsr_copy(matrix_tmp8,matrix_tmp1,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp1,error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7,&
                              matrix_tmp8,&
                              0.0_dp, matrix_tmp4,&
                              filter_eps=eps_multiply,&
                              retain_sparsity=.TRUE.,&
                              error=error)
       
       ! 12. KS_blk=KS_blk-TMP4_blk
       CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),matrix_tmp4,&
                                     1.0_dp,-1.0_dp,error=error)
       
       ! 13. TMP5_blk=tr(TMP5_blk)
       !     KS_blk=KS_blk-tr(TMP4_blk)
       CALL cp_dbcsr_transposed(matrix_tmp5,matrix_tmp4,error=error)
       CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),matrix_tmp5,&
                                     1.0_dp,-1.0_dp,error=error)
       
       ! 14. TMP4_blk=TMP7_blk.tr(TMP8_blk)
       !     Cost: Nnn
       CALL cp_dbcsr_copy(matrix_tmp7,matrix_tmp2,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp2,error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7,&
                                     matrix_tmp8,&
                                     0.0_dp, matrix_tmp4,&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       ! 15. KS_blk=KS_blk+TMP4_blk
       CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),matrix_tmp4,&
                                     1.0_dp,1.0_dp,error=error)
       
       ! 16. KS_blk=KS_blk+tr(TMP4_blk)
       CALL cp_dbcsr_transposed(matrix_tmp5,matrix_tmp4,error=error)
       CALL cp_dbcsr_release(matrix_tmp4,error=error)
       CALL cp_dbcsr_add(almo_scf_env%matrix_ks_blk(ispin),matrix_tmp5,&
                                     1.0_dp,1.0_dp,error=error)
       CALL cp_dbcsr_release(matrix_tmp5,error=error)
                                     
       ! 17. TMP10_blk=TMP8_blk.TMP9_blk
       !    Cost: Noo
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp8,&
                                     matrix_tmp9,&
                                     0.0_dp, matrix_tmp7,&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_release(matrix_tmp9,error=error)
       
       ! 18. KS_blk=TMP7_blk.tr(TMP8_blk)
       !    Cost: Nno
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7,&
                                     matrix_tmp8,&
                                     1.0_dp, almo_scf_env%matrix_ks_blk(ispin),&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_release(matrix_tmp7,error=error)
       CALL cp_dbcsr_release(matrix_tmp8,error=error)

    ENDDO ! spins

    CALL timestop(handle)

  END SUBROUTINE almo_scf_ks_to_ks_blk

! *****************************************************************************
!> \brief ALMOs by diagonalizing the KS domain submatrices
!>        computes both the occupied and virtual orbitals
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_ks_xx_to_tv_xx(almo_scf_env,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_xx_to_tv_xx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_size, idomain, &
                                                info, ispin, lwork, ndomains
    LOGICAL                                  :: failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, work
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: data_copy
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_ks_xx_orthog, subm_t, &
                                                subm_tmp

    CALL timeset(routineN,handle)

    IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular .AND. &
        almo_scf_env%mat_distr_aos==almo_mat_distr_atomic) THEN
       CPErrorMessage(cp_failure_level,routineP,"a domain must be located entirely on a CPU",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ndomains=almo_scf_env%ndomains
    ALLOCATE(subm_tmp(ndomains))
    ALLOCATE(subm_ks_xx_orthog(ndomains))
    ALLOCATE(subm_t(ndomains))

    DO ispin=1,almo_scf_env%nspins

       CALL init_submatrices(subm_tmp,error)
       CALL init_submatrices(subm_ks_xx_orthog,error)
    
       ! TRY: project out T0-occupied space for each domain
       ! F=(1-R_du).F.(1-tr(R_du))
       !CALL copy_submatrices(almo_scf_env%domain_ks_xx(:,ispin),&
       !        subm_ks_xx_orthog,copy_data=.TRUE.,error=error)
       !CALL multiply_submatrices('N','N',1.0_dp,&
       !        almo_scf_env%domain_r_down_up(:,ispin),&
       !        almo_scf_env%domain_ks_xx(:,ispin),0.0_dp,subm_tmp,error)
       !CALL add_submatrices(1.0_dp,subm_ks_xx_orthog,-1.0_dp,subm_tmp,'N',error)
       !CALL add_submatrices(1.0_dp,subm_ks_xx_orthog,-1.0_dp,subm_tmp,'T',error)
       !!CALL multiply_submatrices('N','T',1.0_dp,subm_tmp,&
       !!        almo_scf_env%domain_r_down_up(:,ispin),&
       !!        1.0_dp,subm_ks_xx_orthog,error)

       ! convert blocks to the orthogonal basis set
       ! TRY: replace one multiply
       !CALL multiply_submatrices('N','N',1.0_dp,subm_ks_xx_orthog,&
       !        almo_scf_env%domain_s_sqrt_inv(:,ispin),0.0_dp,subm_tmp,error)
       CALL multiply_submatrices('N','N',1.0_dp,almo_scf_env%domain_ks_xx(:,ispin),&
               almo_scf_env%domain_s_sqrt_inv(:,ispin),0.0_dp,subm_tmp,error)
       CALL multiply_submatrices('N','N',1.0_dp,almo_scf_env%domain_s_sqrt_inv(:,ispin),&
               subm_tmp,0.0_dp,subm_ks_xx_orthog,error)
       CALL release_submatrices(subm_tmp,error)
       
       ! create temporary matrices for occupied and virtual orbitals 
       ! represented in the orthogonalized basis set
       CALL init_submatrices(subm_t,error)
   
       ! loop over domains - perform diagonalization
       DO idomain = 1, ndomains
       
          ! check if the submatrix exists
          IF (subm_ks_xx_orthog(idomain)%domain.gt.0) THEN
   
             iblock_size=subm_ks_xx_orthog(idomain)%nrows
      
             ! Prepare data
             ALLOCATE(eigenvalues(iblock_size))
             ALLOCATE(data_copy(iblock_size,iblock_size))
             data_copy(:,:)=subm_ks_xx_orthog(idomain)%mdata(:,:)
      
             ! Query the optimal workspace for dsyev
             LWORK = -1
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             LWORK = INT(WORK( 1 ))
             DEALLOCATE(WORK)
      
             ! Allocate the workspace and solve the eigenproblem
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             IF( INFO.NE.0 ) THEN
                CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF

!WRITE (*,*) "Domain", idomain ,"OCC energies", eigenvalues( 1:almo_scf_env%nocc_of_domain(idomain,ispin) )
!WRITE (*,*) "Domain", idomain ,"VIR energies", eigenvalues( almo_scf_env%nocc_of_domain(idomain,ispin)+1 : iblock_size )
   
             ! Copy occupied eigenvectors
             IF ( almo_scf_env%domain_t(idomain,ispin)%ncols.NE.&
                almo_scf_env%nocc_of_domain(idomain,ispin) ) THEN
                CPErrorMessage(cp_failure_level,routineP,"wrong domain structure",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
             CALL copy_submatrices(almo_scf_env%domain_t(idomain,ispin),&
                     subm_t(idomain),.FALSE.,error)
             CALL copy_submatrix_data(data_copy(:,1:almo_scf_env%nocc_of_domain(idomain,ispin)),&
                     subm_t(idomain),error)

             DEALLOCATE(WORK)
             DEALLOCATE(data_copy)
             DEALLOCATE(eigenvalues)
      
          ENDIF ! submatrix for the domain exists
   
       ENDDO ! loop over domains

       CALL release_submatrices(subm_ks_xx_orthog,error)

       ! convert orbitals to the AO basis set (from orthogonalized AOs)
       CALL multiply_submatrices('N','N',1.0_dp,almo_scf_env%domain_s_sqrt_inv(:,ispin),&
               subm_t,0.0_dp,almo_scf_env%domain_t(:,ispin),error)
       CALL release_submatrices(subm_t,error)
       
       ! convert domain orbitals to a dbcsr matrix
       CALL construct_dbcsr_from_submatrices(&
               almo_scf_env%matrix_t(ispin),&
               almo_scf_env%domain_t(:,ispin),&
               almo_scf_env%quench_t(ispin),&
               error)
       CALL cp_dbcsr_filter(almo_scf_env%matrix_t(ispin),&
               almo_scf_env%eps_filter,error=error)
       
       ! TRY: add T0 component
       !!CALL cp_dbcsr_add(almo_scf_env%matrix_t(ispin),&
       !!        almo_scf_env%matrix_t_blk(ispin),1.0_dp,1.0_dp,error=error)

    ENDDO ! spins

    DEALLOCATE(subm_tmp)
    DEALLOCATE(subm_ks_xx_orthog)
    DEALLOCATE(subm_t)

    CALL timestop(handle)

  END SUBROUTINE almo_scf_ks_xx_to_tv_xx

! *****************************************************************************
!> \brief computes ALMOs by diagonalizing the projected blocked KS matrix
!>        uses the diagonalization code for blocks
!>        computes both the occupied and virtual orbitals
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_ks_blk_to_tv_blk(almo_scf_env,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_blk_to_tv_blk', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iblock_col, iblock_row, iblock_size, info, ispin, &
      lwork, nocc_of_block, nvirt_of_block, orbital
    LOGICAL                                  :: block_needed, failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, work
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: data_copy
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_ks_blk_orthog, &
                                                matrix_t_blk_orthog, &
                                                matrix_tmp, &
                                                matrix_v_blk_orthog

    CALL timeset(routineN,handle)

    IF (almo_scf_env%domain_layout_aos==almo_domain_layout_molecular .AND. &
        almo_scf_env%mat_distr_aos==almo_mat_distr_atomic) THEN
       CPErrorMessage(cp_failure_level,routineP,"a domain must be located entirely on a CPU",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    DO ispin=1,almo_scf_env%nspins

       CALL cp_dbcsr_init(matrix_tmp,error=error)
       CALL cp_dbcsr_init(matrix_ks_blk_orthog,error=error)
       CALL cp_dbcsr_create(matrix_tmp,template=almo_scf_env%matrix_ks_blk(ispin),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_create(matrix_ks_blk_orthog,template=almo_scf_env%matrix_ks_blk(ispin),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
   
       ! convert blocks to the orthogonal basis set
       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_ks_blk(ispin),&
               almo_scf_env%matrix_s_blk_sqrt_inv(1),0.0_dp,matrix_tmp,&
               filter_eps=almo_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_s_blk_sqrt_inv(1),&
               matrix_tmp,0.0_dp,matrix_ks_blk_orthog,&
               filter_eps=almo_scf_env%eps_filter,error=error)
       
       CALL cp_dbcsr_release(matrix_tmp,error=error)

       ! create temporary matrices for occupied and virtual orbitals 
       ! represented in the orthogonalized AOs basis set
       CALL cp_dbcsr_init(matrix_t_blk_orthog,error=error)
       CALL cp_dbcsr_init(matrix_v_blk_orthog,error=error)
       CALL cp_dbcsr_create(matrix_t_blk_orthog,template=almo_scf_env%matrix_t_blk(ispin),&
               error=error)
       CALL cp_dbcsr_create(matrix_v_blk_orthog,template=almo_scf_env%matrix_v_full_blk(ispin),&
               error=error)
       CALL cp_dbcsr_work_create(matrix_t_blk_orthog,work_mutable=.TRUE.,&
               error=error)
       CALL cp_dbcsr_work_create(matrix_v_blk_orthog,work_mutable=.TRUE.,&
               error=error)
       
       CALL cp_dbcsr_work_create(almo_scf_env%matrix_eoo(ispin),work_mutable=.TRUE.,&
               error=error)
       CALL cp_dbcsr_work_create(almo_scf_env%matrix_evv_full(ispin),work_mutable=.TRUE.,&
               error=error)
   
       CALL cp_dbcsr_iterator_start(iter,matrix_ks_blk_orthog)
   
       DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
          CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,row_size=iblock_size)
   
          block_needed=.FALSE.
       
          IF (iblock_row==iblock_col) THEN
              block_needed=.TRUE.
          ENDIF
   
          IF (.NOT.block_needed) THEN
             CPErrorMessage(cp_failure_level,routineP,"off-diagonal block found",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ENDIF
   
          IF (block_needed) THEN

             ! Prepare data
             ALLOCATE(eigenvalues(iblock_size))
             ALLOCATE(data_copy(iblock_size,iblock_size))
             data_copy(:,:)=data_p(:,:)
   
             ! Query the optimal workspace for dsyev
             LWORK = -1
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             LWORK = INT(WORK( 1 ))
             DEALLOCATE(WORK)
   
             ! Allocate the workspace and solve the eigenproblem
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             IF( INFO.NE.0 ) THEN
                CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF

             !!! RZK-warning                                               !!!
             !!! IT IS EXTREMELY IMPORTANT THAT THE DIAGONAL BLOCKS OF THE !!!
             !!! FOLLOWING MATRICES ARE LOCATED ON THE SAME NODES WITH     !!!
             !!! THE CORRESPONDING DIAGONAL BLOCKS OF THE FOCK MATRIX:     !!!
             !!! T, V, E_o, E_v

             ! copy eigenvectors into two cp_dbcsr matrices - occupied and virtuals
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(matrix_t_blk_orthog,iblock_row,iblock_col,p_new_block)
             nocc_of_block=SIZE(p_new_block,2)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             CPPrecondition(nocc_of_block.gt.0,cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = data_copy(:,1:nocc_of_block)
             ! now virtuals
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(matrix_v_blk_orthog,iblock_row,iblock_col,p_new_block)
             nvirt_of_block=SIZE(p_new_block,2)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             CPPrecondition(nvirt_of_block.gt.0,cp_failure_level,routineP,error,failure)
             !CPPrecondition((nvirt_of_block+nocc_of_block.eq.iblock_size),cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = data_copy(:,(nocc_of_block+1):(nocc_of_block+nvirt_of_block))
   
             ! copy eigenvalues into two diagonal cp_dbcsr matrices - Eoo and Evv
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_eoo(ispin),iblock_row,iblock_col,p_new_block)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = 0.0_dp 
             DO orbital=1,nocc_of_block
                p_new_block(orbital,orbital)=eigenvalues(orbital)
             ENDDO
             ! virtual energies
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(almo_scf_env%matrix_evv_full(ispin),iblock_row,iblock_col,p_new_block)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = 0.0_dp 
             DO orbital=1,nvirt_of_block
                p_new_block(orbital,orbital)=eigenvalues(nocc_of_block+orbital)
             ENDDO

             DEALLOCATE(WORK)
             DEALLOCATE(data_copy)
             DEALLOCATE(eigenvalues)
   
          ENDIF
       
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
   
       CALL cp_dbcsr_finalize(matrix_t_blk_orthog,error=error)
       CALL cp_dbcsr_finalize(matrix_v_blk_orthog,error=error)
       CALL cp_dbcsr_finalize(almo_scf_env%matrix_eoo(ispin),error=error)
       CALL cp_dbcsr_finalize(almo_scf_env%matrix_evv_full(ispin),error=error)
      
       CALL cp_dbcsr_filter(matrix_t_blk_orthog,almo_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_filter(matrix_v_blk_orthog,almo_scf_env%eps_filter,error=error)
       
       CALL cp_dbcsr_release(matrix_ks_blk_orthog, error=error)
   
       ! convert orbitals to the AO basis set (from orthogonalized AOs)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_s_blk_sqrt_inv(1),&
               matrix_t_blk_orthog,0.0_dp,almo_scf_env%matrix_t_blk(ispin),&
               filter_eps=almo_scf_env%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,almo_scf_env%matrix_s_blk_sqrt_inv(1),&
               matrix_v_blk_orthog,0.0_dp,almo_scf_env%matrix_v_full_blk(ispin),&
               filter_eps=almo_scf_env%eps_filter,error=error)
         
       CALL cp_dbcsr_release(matrix_t_blk_orthog, error=error)
       CALL cp_dbcsr_release(matrix_v_blk_orthog, error=error)

    ENDDO ! spins

    CALL timestop(handle)

  END SUBROUTINE almo_scf_ks_blk_to_tv_blk

! *****************************************************************************
!> \brief inverts block-diagonal blocks of a cp_dbcsr_matrix
!> \param matrix_in ...
!> \param matrix_out ...
!> \param nocc ...
!> \param error ...
!> \par History
!>       2012.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE pseudo_invert_diagonal_blk(matrix_in,matrix_out,nocc,error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_in
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_out
    INTEGER, DIMENSION(:)                    :: nocc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pseudo_invert_diagonal_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, iblock_size, &
                                                methodID
    LOGICAL                                  :: failure
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: data_copy
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_create(matrix_out,template=matrix_in,&
            error=error)
    CALL cp_dbcsr_work_create(matrix_out,work_mutable=.TRUE.,&
            error=error)
    
    CALL cp_dbcsr_iterator_start(iter,matrix_in)
   
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       
       CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,row_size=iblock_size)
   
       IF (iblock_row==iblock_col) THEN

          ! Prepare data
          ALLOCATE(data_copy(iblock_size,iblock_size))
          !data_copy(:,:)=data_p(:,:)

          ! 0. Cholesky factorization
          ! 1. Diagonalization
          methodID=1
          CALL pseudo_invert_matrix(data_p,data_copy,iblock_size,&
                  methodID,&
                  range1=nocc(iblock_row),range2=nocc(iblock_row),&
                  !range1_thr,range2_thr,&
                  shift=1.0E-5_dp,&
                  error=error)
   
          !!! IT IS EXTREMELY IMPORTANT THAT THE BLOCKS OF THE "OUT"  !!!
          !!! MATRIX ARE DISTRIBUTED AS THE BLOCKS OF THE "IN" MATRIX !!!

          NULLIFY (p_new_block)
          CALL cp_dbcsr_reserve_block2d(matrix_out,iblock_row,iblock_col,p_new_block)
          CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
          p_new_block(:,:) = data_copy(:,:)
   
          DEALLOCATE(data_copy)
   
       ENDIF
    
    ENDDO
    CALL cp_dbcsr_iterator_stop(iter)
   
    CALL cp_dbcsr_finalize(matrix_out,error=error)
       
    CALL timestop(handle)

  END SUBROUTINE pseudo_invert_diagonal_blk

! *****************************************************************************
!> \brief computes occupied ALMOs from the superimposed atomic density blocks
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_p_blk_to_t_blk(almo_scf_env,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_p_blk_to_t_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(cp_dbcsr_type)                      :: matrix_t_blk_tmp
    TYPE(dbcsr_error_type)                   :: dbcsr_error

    CALL timeset(routineN,handle)

    DO ispin=1,almo_scf_env%nspins

       !! create a temporary matrix to keep the eigenvectors
       !CALL cp_dbcsr_init(matrix_t_blk_tmp, error=error)
       !CALL cp_dbcsr_create(matrix_t_blk_tmp,&
       !                     template=almo_scf_env%matrix_t_blk(ispin),&
       !                     error=error)
       !CALL cp_dbcsr_work_create(matrix_t_blk_tmp,&
       !        work_mutable=.TRUE.,&
       !        error=error)
   
       !CALL cp_dbcsr_iterator_start(iter,almo_scf_env%matrix_p_blk(ispin))
       !DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       !   CALL cp_dbcsr_iterator_next_block(iter,iblock_row,iblock_col,data_p,row_size=iblock_size)
   
       !   block_needed=.FALSE.
       !
       !   IF (iblock_row==iblock_col) THEN
       !       block_needed=.TRUE.
       !   ENDIF
   
       !   IF (.NOT.block_needed) THEN
       !      CPErrorMessage(cp_failure_level,routineP,"off-diagonal block found",error)
       !      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !   ENDIF
   
       !   IF (block_needed) THEN

       !      ! Prepare data
       !      ALLOCATE(eigenvalues(iblock_size))
       !      ALLOCATE(data_copy(iblock_size,iblock_size))
       !      data_copy(:,:)=data_p(:,:)
   
       !      ! Query the optimal workspace for dsyev
       !      LWORK = -1
       !      ALLOCATE(WORK(MAX(1,LWORK)))
       !      CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
       !      LWORK = INT(WORK( 1 ))
       !      DEALLOCATE(WORK)
   
       !      ! Allocate the workspace and solve the eigenproblem
       !      ALLOCATE(WORK(MAX(1,LWORK)))
       !      CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
       !      IF( INFO.NE.0 ) THEN
       !         CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
       !         CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       !      END IF

       !      !!! IT IS EXTREMELY IMPORTANT THAT THE DIAGONAL BLOCKS OF THE !!!
       !      !!! P AND T MATRICES ARE LOCATED ON THE SAME NODES            !!!

       !      ! copy eigenvectors into two cp_dbcsr matrices - occupied and virtuals
       !      NULLIFY (p_new_block)
       !      CALL cp_dbcsr_reserve_block2d(matrix_t_blk_tmp,&
       !              iblock_row,iblock_col,p_new_block)
       !      nocc_of_block=SIZE(p_new_block,2)
       !      CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
       !      CPPrecondition(nocc_of_block.gt.0,cp_failure_level,routineP,error,failure)
       !      p_new_block(:,:) = data_copy(:,iblock_size-nocc_of_block+1:)

       !      DEALLOCATE(WORK)
       !      DEALLOCATE(data_copy)
       !      DEALLOCATE(eigenvalues)
   
       !   ENDIF
       !
       !ENDDO
       !CALL cp_dbcsr_iterator_stop(iter)
   
       !CALL cp_dbcsr_finalize(matrix_t_blk_tmp,error=error)
       !CALL cp_dbcsr_filter(matrix_t_blk_tmp,&
       !        almo_scf_env%eps_filter,error=error)
       !CALL cp_dbcsr_copy(almo_scf_env%matrix_t_blk(ispin),&
       !        matrix_t_blk_tmp,error=error)
       !CALL cp_dbcsr_release(matrix_t_blk_tmp,error=error)
       
       !! generate a random set of ALMOs 
       !! matrix_t_blk should already be initiated to the proper domain structure
       CALL dbcsr_init_random(almo_scf_env%matrix_t_blk(ispin)%matrix,&
                              keep_sparsity=.TRUE.,error=dbcsr_error)

       CALL cp_dbcsr_init(matrix_t_blk_tmp, error=error)
       CALL cp_dbcsr_create(matrix_t_blk_tmp,&
                            template=almo_scf_env%matrix_t_blk(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)

       ! use current ALMOs in matrix_t_blk and project them onto the blocked dm
       ! compute T_new = R_blk S_blk T_random
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_blk(1),&
                              almo_scf_env%matrix_t_blk(ispin),&
                              0.0_dp, matrix_t_blk_tmp,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, &
                  almo_scf_env%matrix_p_blk(ispin), matrix_t_blk_tmp,&
                  0.0_dp, almo_scf_env%matrix_t_blk(ispin),&
                  filter_eps=almo_scf_env%eps_filter,&
                  error=error)

       CALL cp_dbcsr_release(matrix_t_blk_tmp, error=error)

    ENDDO

    CALL timestop(handle)

  END SUBROUTINE almo_scf_p_blk_to_t_blk

! *****************************************************************************
!> \brief Computes the overlap matrix of MO orbitals
!> \param bra ...
!> \param ket ...
!> \param overlap ...
!> \param metric ...
!> \param retain_overlap_sparsity ...
!> \param eps_filter ...
!> \param error ...
!> \par History
!>       2011.08 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE get_overlap(bra,ket,overlap,metric,retain_overlap_sparsity,&
    eps_filter,error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: bra, ket
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: overlap
    TYPE(cp_dbcsr_type), INTENT(IN)          :: metric
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_overlap_sparsity
    REAL(KIND=dp)                            :: eps_filter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_overlap', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: local_retain_sparsity
    TYPE(cp_dbcsr_type)                      :: tmp

    CALL timeset(routineN,handle)

    IF (.NOT.PRESENT(retain_overlap_sparsity)) THEN
       local_retain_sparsity=.FALSE.
    ELSE
       local_retain_sparsity=retain_overlap_sparsity
    ENDIF

    CALL cp_dbcsr_init(tmp,error=error)
    CALL cp_dbcsr_create(tmp,template=ket,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)

    ! TMP=metric*ket
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            metric,ket,0.0_dp,tmp,&
            filter_eps=eps_filter,error=error)

    ! OVERLAP=tr(bra)*TMP
    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
            bra,tmp,0.0_dp,overlap,&
            retain_sparsity=local_retain_sparsity,&
            filter_eps=eps_filter,error=error)
 
    CALL cp_dbcsr_release(tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE get_overlap 

!! *****************************************************************************
!!> \brief Create the overlap matrix of virtual orbitals
!!> \par History
!!>       2011.07 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! *****************************************************************************
!  SUBROUTINE almo_scf_v_to_sigma_vv(almo_scf_env,error)
!
!    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!
!    CHARACTER(LEN=*), PARAMETER :: &
!      routineN = 'almo_scf_v_to_sigma_vv', &
!      routineP = moduleN//':'//routineN
!
!    TYPE(cp_dbcsr_type)                      :: tmp 
!    INTEGER                                  :: ispin, handle
!
!    CALL timeset(routineN,handle)
!
!    DO ispin=1,almo_scf_env%nspins
!   
!       CALL cp_dbcsr_init(tmp, error=error)
!       CALL cp_dbcsr_create(tmp,&
!               template=almo_scf_env%matrix_v(ispin),&
!               matrix_type=dbcsr_type_no_symmetry,&
!               error=error)
!
!       ! TMP=S.V
!       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!               almo_scf_env%matrix_s(1),&
!               almo_scf_env%matrix_v(ispin),&
!               0.0_dp,tmp,&
!               filter_eps=almo_scf_env%eps_filter,error=error)
!
!       ! Sig_vv=tr(V).S.V - get MO overlap
!       CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!               almo_scf_env%matrix_v(ispin),&
!               tmp,&
!               0.0_dp,almo_scf_env%matrix_sigma_vv(ispin),&
!               filter_eps=almo_scf_env%eps_filter,error=error)
!
!       CALL cp_dbcsr_release(tmp,error=error)
!
!    END DO
!  
!    CALL timestop(handle)
!
!  END SUBROUTINE almo_scf_v_to_sigma_vv

!! *****************************************************************************
!!> \brief orthogonalize virtual oribitals within a domain
!!> \par History
!!>       2011.07 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! *****************************************************************************
!  SUBROUTINE almo_scf_v_to_v_orthonormal_blk(almo_scf_env,error)
!
!    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
!    TYPE(cp_error_type), INTENT(INOUT)       :: error
!
!    CHARACTER(LEN=*), PARAMETER :: &
!      routineN = 'almo_scf_v_to_v_orthonormal_blk', &
!      routineP = moduleN//':'//routineN
!
!    TYPE(cp_dbcsr_type)                      :: matrix_v_tmp,&
!                                                sigma_vv_blk_sqrt,&
!                                                sigma_vv_blk_sqrt_inv
!    INTEGER                                  :: ispin, handle
!
!
!    CALL timeset(routineN,handle)
!
!    DO ispin=1,almo_scf_env%nspins
!   
!       CALL cp_dbcsr_init(matrix_v_tmp, error=error)
!       CALL cp_dbcsr_create(matrix_v_tmp,&
!                            template=almo_scf_env%matrix_v(ispin),&
!                            error=error)
!
!       ! TMP=S.V
!       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1),&
!                              almo_scf_env%matrix_v(ispin),&
!                              0.0_dp, matrix_v_tmp,&
!                              filter_eps=almo_scf_env%eps_filter,&
!                              error=error)
!
!       ! Sig_blk=tr(V).TMP - get blocked MO overlap
!       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
!                              almo_scf_env%matrix_v(ispin),&
!                              matrix_v_tmp,&
!                              0.0_dp, almo_scf_env%matrix_sigma_vv_blk(ispin),&
!                              filter_eps=almo_scf_env%eps_filter,&
!                              retain_sparsity=.TRUE.,&
!                              error=error)
!
!       CALL cp_dbcsr_init(sigma_vv_blk_sqrt,error=error)
!       CALL cp_dbcsr_init(sigma_vv_blk_sqrt_inv,error=error)
!       CALL cp_dbcsr_create(sigma_vv_blk_sqrt,template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
!                            matrix_type=dbcsr_type_no_symmetry,error=error) 
!       CALL cp_dbcsr_create(sigma_vv_blk_sqrt_inv,template=almo_scf_env%matrix_sigma_vv_blk(ispin),&
!                            matrix_type=dbcsr_type_no_symmetry,error=error) 
!
!       ! compute sqrt and sqrt_inv of the blocked MO overlap
!       CALL matrix_sqrt_Newton_Schulz(sigma_vv_blk_sqrt,sigma_vv_blk_sqrt_inv,&
!                                      almo_scf_env%matrix_sigma_vv_blk(ispin),&
!                                      threshold=almo_scf_env%eps_filter,&
!                                      order=almo_scf_env%order_lanczos,&
!                                      eps_lanczos=almo_scf_env%eps_lancsoz,&
!                                      max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
!                                      error=error)
!       
!       ! TMP_blk=V.SigSQRTInv_blk
!       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
!                              almo_scf_env%matrix_v(ispin),&
!                              sigma_vv_blk_sqrt_inv,&
!                              0.0_dp, matrix_v_tmp,&
!                              filter_eps=almo_scf_env%eps_filter,&
!                              error=error)
!
!       ! update the orbitals with the orthonormalized MOs
!       CALL cp_dbcsr_copy(almo_scf_env%matrix_v(ispin),matrix_v_tmp,&
!                          error=error)
!
!       CALL cp_dbcsr_release (matrix_v_tmp, error=error)
!       CALL cp_dbcsr_release (sigma_vv_blk_sqrt, error=error)
!       CALL cp_dbcsr_release (sigma_vv_blk_sqrt_inv, error=error)
!
!    END DO
!  
!    CALL timestop(handle)
!
!  END SUBROUTINE almo_scf_v_to_v_orthonormal_blk

! *****************************************************************************
!> \brief orthogonalize MOs
!> \param ket ...
!> \param overlap ...
!> \param metric ...
!> \param retain_locality ...
!> \param only_normalize ...
!> \param eps_filter ...
!> \param order_lanczos ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param error ...
!> \par History
!>       2012.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE orthogonalize_mos(ket,overlap,metric,retain_locality,only_normalize,&
    eps_filter,order_lanczos,eps_lanczos,max_iter_lanczos,error)

    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: ket, overlap
    TYPE(cp_dbcsr_type), INTENT(IN)          :: metric
    LOGICAL, INTENT(IN)                      :: retain_locality, &
                                                only_normalize
    REAL(KIND=dp)                            :: eps_filter
    INTEGER, INTENT(IN)                      :: order_lanczos
    REAL(KIND=dp), INTENT(IN)                :: eps_lanczos
    INTEGER, INTENT(IN)                      :: max_iter_lanczos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'orthogonalize_mos', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dim0, handle
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: diagonal
    TYPE(cp_dbcsr_type)                      :: matrix_sigma_blk_sqrt, &
                                                matrix_sigma_blk_sqrt_inv, &
                                                matrix_t_blk_tmp

    CALL timeset(routineN,handle)

    ! create block-diagonal sparsity pattern for the overlap
    ! in case retain_locality is set to true
    ! RZK-warning this will fail if distribution blocks are smaller than domains!!! 
    CALL cp_dbcsr_set(overlap,0.0_dp,error=error)
    CALL cp_dbcsr_add_on_diag(overlap,1.0_dp,error=error)
    CALL cp_dbcsr_filter(overlap,eps_filter,error=error)

    CALL get_overlap(ket,ket,overlap,metric,retain_locality,&
            eps_filter,error)

    IF (only_normalize) THEN

       CALL cp_dbcsr_get_info(overlap, nfullrows_total=dim0 )
       ALLOCATE(diagonal(dim0))
       CALL cp_dbcsr_get_diag(overlap,diagonal,error=error)
       CALL cp_dbcsr_set(overlap,0.0_dp,error=error)
       CALL cp_dbcsr_set_diag(overlap,diagonal,error=error)
       DEALLOCATE(diagonal)
       CALL cp_dbcsr_filter(overlap,eps_filter,error=error)

    ENDIF

    CALL cp_dbcsr_init(matrix_sigma_blk_sqrt,error=error)
    CALL cp_dbcsr_init(matrix_sigma_blk_sqrt_inv,error=error)
    CALL cp_dbcsr_create(matrix_sigma_blk_sqrt,template=overlap,&
            matrix_type=dbcsr_type_no_symmetry,error=error) 
    CALL cp_dbcsr_create(matrix_sigma_blk_sqrt_inv,template=overlap,&
            matrix_type=dbcsr_type_no_symmetry,error=error) 

    ! compute sqrt and sqrt_inv of the blocked MO overlap
    CALL matrix_sqrt_Newton_Schulz(matrix_sigma_blk_sqrt,matrix_sigma_blk_sqrt_inv,&
            overlap,threshold=eps_filter,&
            order=order_lanczos,&
            eps_lanczos=eps_lanczos,&
            max_iter_lanczos=max_iter_lanczos,&
            error=error)
       
    CALL cp_dbcsr_init(matrix_t_blk_tmp, error=error)
    CALL cp_dbcsr_create(matrix_t_blk_tmp,&
            template=ket,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)

    CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
            ket,&
            matrix_sigma_blk_sqrt_inv,&
            0.0_dp, matrix_t_blk_tmp,&
            filter_eps=eps_filter,&
            !retain_sparsity=retain_locality,&
            error=error)

    ! update the orbitals with the orthonormalized MOs
    CALL cp_dbcsr_copy(ket,matrix_t_blk_tmp,&
            error=error)

    CALL cp_dbcsr_release (matrix_t_blk_tmp, error=error)
    CALL cp_dbcsr_release (matrix_sigma_blk_sqrt, error=error)
    CALL cp_dbcsr_release (matrix_sigma_blk_sqrt_inv, error=error)

    CALL timestop(handle)

  END SUBROUTINE orthogonalize_mos 

! *****************************************************************************
!> \brief orthogonalize ALMOs within a domain (obsolete, use orthogonalize_mos)
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_t_blk_to_t_blk_orthonormal(almo_scf_env,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'almo_scf_t_blk_to_t_blk_orthonormal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(cp_dbcsr_type)                      :: matrix_sigma_blk_sqrt, &
                                                matrix_sigma_blk_sqrt_inv, &
                                                matrix_t_blk_tmp

    CALL timeset(routineN,handle)

    DO ispin=1,almo_scf_env%nspins
   
       CALL cp_dbcsr_init(matrix_t_blk_tmp, error=error)
       CALL cp_dbcsr_create(matrix_t_blk_tmp,&
                            template=almo_scf_env%matrix_t_blk(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)

       ! TMP_blk=S_blk.T_blk
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_blk(1),&
                              almo_scf_env%matrix_t_blk(ispin),&
                              0.0_dp, matrix_t_blk_tmp,&
                              filter_eps=almo_scf_env%eps_filter,&
                              error=error)

       ! Sig_blk=tr(T_blk).TMP_blk - get blocked MO overlap
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                              almo_scf_env%matrix_t_blk(ispin),&
                              matrix_t_blk_tmp,&
                              0.0_dp, almo_scf_env%matrix_sigma_blk(ispin),&
                              filter_eps=almo_scf_env%eps_filter,&
                              retain_sparsity=.TRUE.,&
                              error=error)

       ! RZK-warning try to use symmetry of the sqrt and sqrt_inv matrices
       CALL cp_dbcsr_init(matrix_sigma_blk_sqrt,error=error)
       CALL cp_dbcsr_init(matrix_sigma_blk_sqrt_inv,error=error)
       CALL cp_dbcsr_create(matrix_sigma_blk_sqrt,template=almo_scf_env%matrix_sigma_blk(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,error=error) 
       CALL cp_dbcsr_create(matrix_sigma_blk_sqrt_inv,template=almo_scf_env%matrix_sigma_blk(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,error=error) 

       ! compute sqrt and sqrt_inv of the blocked MO overlap
       CALL matrix_sqrt_Newton_Schulz(matrix_sigma_blk_sqrt,matrix_sigma_blk_sqrt_inv,&
                                      almo_scf_env%matrix_sigma_blk(ispin),&
                                      threshold=almo_scf_env%eps_filter,&
                                      order=almo_scf_env%order_lanczos,&
                                      eps_lanczos=almo_scf_env%eps_lanczos,&
                                      max_iter_lanczos=almo_scf_env%max_iter_lanczos,&
                                      error=error)
       
       ! TMP_blk=T_blk.SigSQRTInv_blk
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              almo_scf_env%matrix_t_blk(ispin),&
                              matrix_sigma_blk_sqrt_inv,&
                              0.0_dp, matrix_t_blk_tmp,&
                              filter_eps=almo_scf_env%eps_filter,&
                              retain_sparsity=.TRUE.,&
                              error=error)

       ! update the orbitals with the orthonormalized ALMOs
       CALL cp_dbcsr_copy(almo_scf_env%matrix_t_blk(ispin),matrix_t_blk_tmp,&
                          keep_sparsity=.TRUE.,&
                          error=error)

       CALL cp_dbcsr_release (matrix_t_blk_tmp, error=error)
       CALL cp_dbcsr_release (matrix_sigma_blk_sqrt, error=error)
       CALL cp_dbcsr_release (matrix_sigma_blk_sqrt_inv, error=error)

    END DO
  
    CALL timestop(handle)

  END SUBROUTINE almo_scf_t_blk_to_t_blk_orthonormal

! *****************************************************************************
!> \brief computes the idempotent density matrix from MOs
!>        MOs can be either orthogonal or non-orthogonal
!> \param t ...
!> \param p ...
!> \param eps_filter ...
!> \param orthog_orbs ...
!> \param s ...
!> \param sigma ...
!> \param sigma_inv ...
!> \param use_guess ...
!> \param error ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_t_to_p(t,p,eps_filter,orthog_orbs,s,sigma,sigma_inv,&
    use_guess,error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: t
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: p
    REAL(KIND=dp), INTENT(IN)                :: eps_filter
    LOGICAL, INTENT(IN)                      :: orthog_orbs
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: s
    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      OPTIONAL                               :: sigma, sigma_inv
    LOGICAL, INTENT(IN), OPTIONAL            :: use_guess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_t_to_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, use_sigma_inv_guess
    TYPE(cp_dbcsr_type)                      :: t_tmp

    CALL timeset(routineN,handle)

    ! make sure that S, sigma and sigma_inv are present for non-orthogonal orbitals
    IF (.NOT.orthog_orbs) THEN
       IF ((.NOT.PRESENT(s)).OR.(.NOT.PRESENT(sigma)).OR.(.NOT.PRESENT(sigma_inv))) THEN
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       ENDIF
    ENDIF

    use_sigma_inv_guess=.FALSE.
    IF (PRESENT(use_guess)) THEN
       use_sigma_inv_guess=use_guess
    ENDIF

    IF (orthog_orbs) THEN
    
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp,t,t,&
                              0.0_dp,p,filter_eps=eps_filter,&
                              error=error)

    ELSE

       CALL cp_dbcsr_init(t_tmp, error=error)
       CALL cp_dbcsr_create(t_tmp,template=t,error=error)
   
       ! TMP=S.T
       CALL cp_dbcsr_multiply("N","N",1.0_dp,s,t,0.0_dp,t_tmp,&
                              filter_eps=eps_filter,&
                              error=error)
   
       ! Sig=tr(T).TMP - get MO overlap
       CALL cp_dbcsr_multiply("T","N",1.0_dp,t,t_tmp,0.0_dp,sigma,&
                              filter_eps=eps_filter,&
                              error=error)
   
       ! invert MO overlap
       CALL invert_Hotelling(&
               matrix_inverse=sigma_inv,&
               matrix=sigma,&
               use_inv_as_guess=use_sigma_inv_guess,&
               threshold=eps_filter,&
               error=error)
  
       ! TMP=T.SigInv
       CALL cp_dbcsr_multiply("N","N",1.0_dp,t,sigma_inv,0.0_dp,t_tmp,&
                              filter_eps=eps_filter,&
                              error=error)
   
       ! P=TMP.tr(T_blk)
       CALL cp_dbcsr_multiply("N","T",1.0_dp,t_tmp,t,0.0_dp,p,&
                              filter_eps=eps_filter,&
                              error=error)
   
       CALL cp_dbcsr_release (t_tmp, error=error)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_t_to_p

! *****************************************************************************
!> \brief applies projector to the orbitals
!>        |psi_out> = P |psi_in>   OR   |psi_out> = (1-P) |psi_in>,
!>        where P = |psi_proj> (<psi_proj|psi_roj>)^{-1} <psi_proj|
!> \param psi_in ...
!> \param psi_out ...
!> \param psi_projector ...
!> \param metric ...
!> \param project_out ...
!> \param psi_projector_orthogonal ...
!> \param proj_in_template ...
!> \param eps_filter ...
!> \param sig_inv_projector ...
!> \param sig_inv_template ...
!> \param error ...
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE apply_projector(psi_in,psi_out,psi_projector,metric,project_out,&
    psi_projector_orthogonal,proj_in_template,eps_filter,sig_inv_projector,&
    sig_inv_template,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: psi_in
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: psi_out
    TYPE(cp_dbcsr_type), INTENT(IN)          :: psi_projector, metric
    LOGICAL, INTENT(IN)                      :: project_out, &
                                                psi_projector_orthogonal
    TYPE(cp_dbcsr_type), INTENT(IN)          :: proj_in_template
    REAL(KIND=dp), INTENT(IN)                :: eps_filter
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: sig_inv_projector, &
                                                sig_inv_template
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_projector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: tmp_no, tmp_ov, tmp_ov2, &
                                                tmp_sig, tmp_sig_inv

    CALL timeset(routineN,handle)

    ! =S*PSI_proj
    CALL cp_dbcsr_init(tmp_no, error=error)
    CALL cp_dbcsr_create(tmp_no,template=psi_projector,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            metric,psi_projector,&
            0.0_dp,tmp_no,&
            filter_eps=eps_filter,&
            error=error)
   
    ! =tr(S.PSI_proj)*PSI_in
    CALL cp_dbcsr_init(tmp_ov,error=error)
    CALL cp_dbcsr_create(tmp_ov,template=proj_in_template,error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
            tmp_no,psi_in,&
            0.0_dp,tmp_ov,&
            filter_eps=eps_filter,&
            error=error)

    IF (.NOT.psi_projector_orthogonal) THEN
       ! =SigInv_proj*Sigma_OV
       CALL cp_dbcsr_init(tmp_sig_inv,error=error)
       CALL cp_dbcsr_init(tmp_ov2,error=error)
       CALL cp_dbcsr_create(tmp_ov2,&
               template=proj_in_template,error=error)
       IF (PRESENT(sig_inv_projector)) THEN
          CALL cp_dbcsr_create(tmp_sig_inv,&
                  template=sig_inv_projector,&
                  error=error)
          CALL cp_dbcsr_copy(tmp_sig_inv,sig_inv_projector,error=error)
       ELSE
          IF (.NOT.PRESENT(sig_inv_template)) THEN
             CPErrorMessage(cp_failure_level,routineP,"PROGRAMMING ERROR: provide either template or sig_inv",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ENDIF
          ! compute inverse overlap of the projector orbitals
          CALL cp_dbcsr_init(tmp_sig,error=error)
          CALL cp_dbcsr_create(tmp_sig,&
                  template=sig_inv_template,&
                  matrix_type=dbcsr_type_no_symmetry,&
                  error=error)
          CALL cp_dbcsr_multiply("T","N",1.0_dp,&
                  psi_projector,tmp_no,0.0_dp,tmp_sig,&
                  filter_eps=eps_filter,&
                  error=error)
          CALL cp_dbcsr_create(tmp_sig_inv,&
                  template=sig_inv_template,&
                  matrix_type=dbcsr_type_no_symmetry,&
                  error=error)
          CALL invert_Hotelling(tmp_sig_inv,tmp_sig,&
                  threshold=eps_filter,error=error)
          CALL cp_dbcsr_release(tmp_sig,error=error)
       ENDIF
       CALL cp_dbcsr_multiply("N","N",1.0_dp,&
               tmp_sig_inv,tmp_ov,0.0_dp,tmp_ov2,&
               filter_eps=eps_filter,error=error)
       CALL cp_dbcsr_release(tmp_sig_inv,error=error)
       CALL cp_dbcsr_copy(tmp_ov,tmp_ov2,error=error)
       CALL cp_dbcsr_release(tmp_ov2,error=error)
    ENDIF
    CALL cp_dbcsr_release(tmp_no,error=error)

    ! =PSI_proj*TMP_OV
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            psi_projector,tmp_ov,0.0_dp,psi_out,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_release(tmp_ov,error=error)

    ! V_out=V_in-V_out
    IF (project_out) THEN
       CALL cp_dbcsr_add(psi_out,psi_in,-1.0_dp,+1.0_dp,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE apply_projector

!! *****************************************************************************
!!> \brief projects the occupied space out from the provided orbitals
!!> \par History
!!>       2011.07 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! *****************************************************************************
!  SUBROUTINE almo_scf_p_out_from_v(v_in,v_out,ov_template,ispin,almo_scf_env,error)
!    
!    TYPE(cp_dbcsr_type), INTENT(IN)                :: v_in, ov_template
!    TYPE(cp_dbcsr_type), INTENT(INOUT)             :: v_out
!    INTEGER, INTENT(IN)                            :: ispin
!    TYPE(almo_scf_env_type), INTENT(INOUT)         :: almo_scf_env
!    TYPE(cp_error_type), INTENT(INOUT)             :: error
!
!    CHARACTER(LEN=*), PARAMETER :: &
!      routineN = 'almo_scf_p_out_from_v', &
!      routineP = moduleN//':'//routineN
!
!    TYPE(cp_dbcsr_type)                      :: tmp_on, tmp_ov, tmp_ov2
!    INTEGER                                  :: handle
!    LOGICAL                                  :: failure
!       
!    CALL timeset(routineN,handle)
!
!    ! =tr(T_blk)*S
!    CALL cp_dbcsr_init(tmp_on, error=error)
!    CALL cp_dbcsr_create(tmp_on,&
!            template=almo_scf_env%matrix_t_tr(ispin),error=error)
!    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
!            almo_scf_env%matrix_t_blk(ispin),&
!            almo_scf_env%matrix_s(1),&
!            0.0_dp,tmp_on,&
!            filter_eps=almo_scf_env%eps_filter,&
!            error=error)
!   
!    ! =tr(T_blk).S*V_in
!    CALL cp_dbcsr_init(tmp_ov,error=error)
!    CALL cp_dbcsr_create(tmp_ov,template=ov_template,error=error)
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!            tmp_on,v_in,0.0_dp,tmp_ov,&
!            filter_eps=almo_scf_env%eps_filter,&
!            error=error)
!    CALL cp_dbcsr_release(tmp_on,error=error)
!
!    ! =SigmaInv*Sigma_OV
!    CALL cp_dbcsr_init(tmp_ov2, error=error)
!    CALL cp_dbcsr_create(tmp_ov2,template=ov_template,error=error)
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!            almo_scf_env%matrix_sigma_inv(ispin),&
!            tmp_ov,0.0_dp,tmp_ov2,&
!            filter_eps=almo_scf_env%eps_filter,&
!            error=error)
!    CALL cp_dbcsr_release(tmp_ov,error=error)
!    
!    ! =T_blk*SigmaInv.Sigma_OV
!    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
!            almo_scf_env%matrix_t_blk(ispin),&
!            tmp_ov2,0.0_dp,v_out,&
!            filter_eps=almo_scf_env%eps_filter,&
!            error=error)
!    CALL cp_dbcsr_release(tmp_ov2,error=error)
!
!    ! V_out=V_in-V_out=
!    CALL cp_dbcsr_add(v_out,v_in,-1.0_dp,+1.0_dp,error=error)
!
!    CALL timestop(handle)
!
!  END SUBROUTINE almo_scf_p_out_from_v

! *****************************************************************************
!> \brief computes the idempotent density matrix from ALMOs
!> \param almo_scf_env ...
!> \param use_sigma_inv_guess ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!>       2011.07 converted into a wrapper which calls almo_scf_t_to_p
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_t_blk_to_p(almo_scf_env,use_sigma_inv_guess,error)

    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    LOGICAL, INTENT(IN), OPTIONAL            :: use_sigma_inv_guess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_t_blk_to_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: use_guess
    REAL(KIND=dp)                            :: spin_factor

    CALL timeset(routineN,handle)

    use_guess=.FALSE.
    IF (PRESENT(use_sigma_inv_guess)) THEN
       use_guess=use_sigma_inv_guess
    ENDIF

    DO ispin=1,almo_scf_env%nspins
   
       CALL almo_scf_t_to_p(t=almo_scf_env%matrix_t_blk(ispin),&
                            p=almo_scf_env%matrix_p(ispin),&
                            eps_filter=almo_scf_env%eps_filter,&
                            orthog_orbs=.FALSE.,&
                            s=almo_scf_env%matrix_s(1),&
                            sigma=almo_scf_env%matrix_sigma(ispin),&
                            sigma_inv=almo_scf_env%matrix_sigma_inv(ispin),&
                            use_guess=use_guess,&
                            error=error)
       
       IF (almo_scf_env%nspins == 1) THEN
        spin_factor = 2.0_dp
       ELSE
        spin_factor = 1.0_dp
       ENDIF
       CALL cp_dbcsr_scale(almo_scf_env%matrix_p(ispin),spin_factor,&
                              error=error)

    END DO
  
    CALL timestop(handle)

  END SUBROUTINE almo_scf_t_blk_to_p

! *****************************************************************************
!> \brief computes a unitary matrix from an arbitrary "generator" matrix
!>        U = ( 1 - X + tr(X) ) ( 1 + X - tr(X) )^(-1)
!> \param X ...
!> \param U ...
!> \param eps_filter ...
!> \param error ...
!> \par History
!>       2011.08 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE generator_to_unitary(X,U,eps_filter,error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: X
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: U
    REAL(KIND=dp), INTENT(IN)                :: eps_filter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'generator_to_unitary', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    LOGICAL                                  :: safe_mode
    REAL(KIND=dp)                            :: frob_matrix, frob_matrix_base
    TYPE(cp_dbcsr_type)                      :: delta, t1, t2, tmp1
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    safe_mode=.TRUE.

    ! get a useful output_unit
    logger => cp_error_get_logger(error)
    IF (logger%para_env%mepos==logger%para_env%source) THEN
       unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.)
    ELSE
       unit_nr=-1
    ENDIF
 
    CALL cp_dbcsr_init(t1,error=error)
    CALL cp_dbcsr_create(t1,template=X,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(t2,error=error)
    CALL cp_dbcsr_create(t2,template=X,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    
    ! create antisymmetric Delta = -X + tr(X)
    CALL cp_dbcsr_init(delta,error=error)
    CALL cp_dbcsr_create(delta,template=X,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_transposed(delta,X,error=error)
! check that transposed is added correctly
    CALL cp_dbcsr_add(delta,X,1.0_dp,-1.0_dp,error=error)
    
    ! compute (1 - Delta)^(-1)
    CALL cp_dbcsr_add_on_diag(t1,1.0_dp,error=error)
    CALL cp_dbcsr_add(t1,delta,1.0_dp,-1.0_dp,error=error)
    CALL invert_Hotelling(t2,t1,threshold=eps_filter,error=error)
    
    IF (safe_mode) THEN

       CALL cp_dbcsr_init(tmp1,error=error)
       CALL cp_dbcsr_create(tmp1,template=X,&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, t2, t1, 0.0_dp, tmp1,&
               filter_eps=eps_filter,error=error)
       frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
       IF (unit_nr>0) THEN
          WRITE(unit_nr,*) "Error for (inv(A)*A-I)",frob_matrix/frob_matrix_base
       ENDIF
       CALL cp_dbcsr_release(tmp1,error=error)
    ENDIF

    CALL cp_dbcsr_multiply("N","N",1.0_dp,delta,t2,0.0_dp,U,&
            filter_eps=eps_filter,error=error)
    CALL cp_dbcsr_add(U,t2,1.0_dp,1.0_dp,error=error)

    IF (safe_mode) THEN

       CALL cp_dbcsr_init(tmp1,error=error)
       CALL cp_dbcsr_create(tmp1,template=X,&
               matrix_type=dbcsr_type_no_symmetry,error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp, U, U, 0.0_dp, tmp1,&
               filter_eps=eps_filter,error=error)
       frob_matrix_base=cp_dbcsr_frobenius_norm(tmp1)
       CALL cp_dbcsr_add_on_diag(tmp1,-1.0_dp,error=error)
       frob_matrix=cp_dbcsr_frobenius_norm(tmp1)
       IF (unit_nr>0) THEN
          WRITE(unit_nr,*) "Error for (trn(U)*U-I)",frob_matrix/frob_matrix_base
       ENDIF
       CALL cp_dbcsr_release(tmp1,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE generator_to_unitary

! *****************************************************************************
!> \brief Parallel code for domain specific operations (my_action)
!>         0. out = op1 * in
!>         1. out = in - op2 * op1 * in
!> \param matrix_in ...
!> \param matrix_out ...
!> \param operator1 ...
!> \param operator2 ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param my_action ...
!> \param filter_eps ...
!> \param matrix_trimmer ...
!> \param use_trimmer ...
!> \param error ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE apply_domain_operators(matrix_in,matrix_out,operator1,operator2,&
    dpattern,map,node_of_domain,my_action,filter_eps,matrix_trimmer,use_trimmer,&
    error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_in
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_out
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN)               :: operator1
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN), OPTIONAL     :: operator2
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    INTEGER, INTENT(IN)                      :: my_action
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: matrix_trimmer
    LOGICAL, INTENT(IN), OPTIONAL            :: use_trimmer
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_domain_operators', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndomains
    LOGICAL                                  :: failure, &
                                                matrix_trimmer_required, &
                                                my_use_trimmer, &
                                                operator2_required
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_in, subm_out, subm_temp

    CALL timeset(routineN,handle)
    
    my_use_trimmer=.FALSE.
    IF (PRESENT(use_trimmer)) THEN
       my_use_trimmer=use_trimmer
    ENDIF
    
    operator2_required=.FALSE.
    matrix_trimmer_required=.FALSE.

    IF (my_action.eq.1) operator2_required=.TRUE.

    IF (my_use_trimmer) THEN
       matrix_trimmer_required=.TRUE.
       CPErrorMessage(cp_failure_level,routineP,"TRIMMED PROJECTOR DISABLED!",error)
    ENDIF

    IF (.NOT.PRESENT(operator2).AND.operator2_required) THEN
       CPErrorMessage(cp_failure_level,routineP,"SECOND OPERATOR IS REQUIRED",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (.NOT.PRESENT(matrix_trimmer).AND.matrix_trimmer_required) THEN
       CPErrorMessage(cp_failure_level,routineP,"TRIMMER MATRIX IS REQUIRED",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ndomains = cp_dbcsr_nblkcols_total(dpattern)

    ALLOCATE(subm_in(ndomains))
    ALLOCATE(subm_temp(ndomains))
    ALLOCATE(subm_out(ndomains))
    !!!TRIM ALLOCATE(subm_trimmer(ndomains))
    CALL init_submatrices(subm_in,error)
    CALL init_submatrices(subm_temp,error)
    CALL init_submatrices(subm_out,error)

    CALL construct_submatrices(matrix_in,subm_in,&
            dpattern,map,node_of_domain,select_row,error)
    
    !!!TRIM IF (matrix_trimmer_required) THEN
    !!!TRIM    CALL construct_submatrices(matrix_trimmer,subm_trimmer,&
    !!!TRIM            dpattern,map,node_of_domain,select_row,error)
    !!!TRIM ENDIF

    IF (my_action.eq.0) THEN
       ! for example, apply preconditioner
       CALL multiply_submatrices('N','N',1.0_dp,operator1,&
               subm_in,0.0_dp,subm_out,error)
    ELSE IF (my_action.eq.1) THEN
       ! use for projectors
       CALL copy_submatrices(subm_in,subm_out,.TRUE.,error)
       CALL multiply_submatrices('N','N',1.0_dp,operator1,&
               subm_in,0.0_dp,subm_temp,error)
       CALL multiply_submatrices('N','N',-1.0_dp,operator2,&
               subm_temp,1.0_dp,subm_out,error)
!GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
!   cp_dbcsr_distribution(dpattern)))
!CALL print_submatrices(subm_out,GroupID,error)                   
    ELSE
       CPErrorMessage(cp_failure_level,routineP,"ILLEGAL ACTION",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL construct_dbcsr_from_submatrices(matrix_out,subm_out,dpattern,error)
    CALL cp_dbcsr_filter(matrix_out,filter_eps,error=error)

    CALL release_submatrices(subm_out,error)
    CALL release_submatrices(subm_temp,error)
    CALL release_submatrices(subm_in,error)

    DEALLOCATE(subm_out)
    DEALLOCATE(subm_temp)
    DEALLOCATE(subm_in)

    CALL timestop(handle)
  
  END SUBROUTINE apply_domain_operators

! *****************************************************************************
!> \brief Constructs preconditioners for each domain
!>        -1. projected preconditioner
!>         0. simple preconditioner
!> \param matrix_main ...
!> \param subm_s_inv ...
!> \param subm_r_down ...
!> \param matrix_trimmer ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param preconditioner ...
!> \param use_trimmer ...
!> \param my_action ...
!> \param error ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_domain_preconditioner(matrix_main,subm_s_inv,&
    subm_r_down,matrix_trimmer,dpattern,map,node_of_domain,preconditioner,&
    use_trimmer,my_action,error)
    
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_main
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(IN), OPTIONAL     :: subm_s_inv, subm_r_down
    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      OPTIONAL                               :: matrix_trimmer
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: preconditioner
    LOGICAL, INTENT(IN), OPTIONAL            :: use_trimmer
    INTEGER, INTENT(IN)                      :: my_action
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'construct_domain_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, idomain, index1_end, index1_start, n_domain_mos, naos, &
      nblkrows_tot, ndomains, neighbor, row
    INTEGER, DIMENSION(:), POINTER           :: nmos
    LOGICAL :: failure, matrix_r_required, matrix_s_inv_required, &
      matrix_trimmer_required, my_use_trimmer
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Minv
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_main, subm_tmp, subm_tmp2

    CALL timeset(routineN,handle)
    
    my_use_trimmer=.FALSE.
    IF (PRESENT(use_trimmer)) THEN
       my_use_trimmer=use_trimmer
    ENDIF
    
    matrix_s_inv_required=.FALSE.
    matrix_trimmer_required=.FALSE.
    matrix_r_required=.FALSE.

    IF (my_action.eq.-1) matrix_s_inv_required=.TRUE.
    IF (my_action.eq.-1) matrix_r_required=.TRUE.
    IF (my_use_trimmer) THEN
       matrix_trimmer_required=.TRUE.
       CPErrorMessage(cp_failure_level,routineP,"TRIMMED PRECONDITIONER DISABLED!",error)
    ENDIF

    IF (.NOT.PRESENT(subm_s_inv).AND.matrix_s_inv_required) THEN
       CPErrorMessage(cp_failure_level,routineP,"S_inv SUBMATRICES ARE REQUIRED",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (.NOT.PRESENT(subm_r_down).AND.matrix_r_required) THEN
       CPErrorMessage(cp_failure_level,routineP,"R SUBMATRICES ARE REQUIRED",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (.NOT.PRESENT(matrix_trimmer).AND.matrix_trimmer_required) THEN
       CPErrorMessage(cp_failure_level,routineP,"TRIMMER MATRIX IS REQUIRED",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ndomains = cp_dbcsr_nblkcols_total(dpattern)
    nblkrows_tot = cp_dbcsr_nblkrows_total(dpattern)
    nmos => array_data(cp_dbcsr_col_block_sizes(dpattern))

    ALLOCATE(subm_main(ndomains))
    CALL init_submatrices(subm_main,error)
    !!!TRIM ALLOCATE(subm_trimmer(ndomains))

    CALL construct_submatrices(matrix_main,subm_main,&
            dpattern,map,node_of_domain,select_row_col,error)
                   
    !!!TRIM IF (matrix_trimmer_required) THEN
    !!!TRIM    CALL construct_submatrices(matrix_trimmer,subm_trimmer,&
    !!!TRIM            dpattern,map,node_of_domain,select_row,error)
    !!!TRIM ENDIF

    IF (my_action.eq.-1) THEN
       ! project out the local occupied space
       !tmp=MATMUL(subm_r(idomain)%mdata,Minv)
       !Minv=MATMUL(tmp,subm_main(idomain)%mdata)
       !subm_main(idomain)%mdata=subm_main(idomain)%mdata-&
       !   Minv-TRANSPOSE(Minv)+MATMUL(Minv,TRANSPOSE(tmp))
       ALLOCATE(subm_tmp(ndomains))
       ALLOCATE(subm_tmp2(ndomains))
       CALL init_submatrices(subm_tmp,error)
       CALL init_submatrices(subm_tmp2,error)
       CALL multiply_submatrices('N','N',1.0_dp,subm_r_down,&
               subm_s_inv,0.0_dp,subm_tmp,error)
       CALL multiply_submatrices('N','N',1.0_dp,subm_tmp,&
               subm_main,0.0_dp,subm_tmp2,error)
       CALL add_submatrices(1.0_dp,subm_main,-1.0_dp,subm_tmp2,'N',error)
       CALL add_submatrices(1.0_dp,subm_main,-1.0_dp,subm_tmp2,'T',error)
       CALL multiply_submatrices('N','T',1.0_dp,subm_tmp2,&
               subm_tmp,1.0_dp,subm_main,error)
       CALL release_submatrices(subm_tmp,error)
       CALL release_submatrices(subm_tmp2,error)
       DEALLOCATE(subm_tmp2)
       DEALLOCATE(subm_tmp)
    ENDIF

    ! loop over domains - perform inversion
    DO idomain = 1, ndomains
    
       ! check if the submatrix exists
       IF (subm_main(idomain)%domain.gt.0) THEN

          ! find sizes of MO submatrices
          IF (idomain.eq.1) THEN
             index1_start = 1
          ELSE
             index1_start = map%index1(idomain-1)
          ENDIF
          index1_end=map%index1(idomain)-1

          n_domain_mos=0
          DO row = index1_start, index1_end
             neighbor=map%pairs(row,1)
             n_domain_mos=n_domain_mos+nmos(neighbor)
          ENDDO
       
          naos=subm_main(idomain)%nrows
          !WRITE(*,*) "Domain, mo_self_and_neig, ao_domain: ", idomain, n_domain_mos, naos

          ALLOCATE(Minv(naos,naos))

          
          !!!TRIM IF (my_use_trimmer) THEN
          !!!TRIM    ! THIS IS SUPER EXPENSIVE (ELIMINATE)
          !!!TRIM    ! trim the main matrix before inverting
          !!!TRIM    ! assume that the trimmer columns are different (i.e. the main matrix is different for each MO)
          !!!TRIM    allocate(tmp(naos,nmos(idomain)))
          !!!TRIM    DO ii=1, nmos(idomain)
          !!!TRIM       ! transform the main matrix using the trimmer for the current MO
          !!!TRIM       DO jj=1, naos
          !!!TRIM          DO kk=1, naos
          !!!TRIM             Mstore(jj,kk)=sumb_main(idomain)%mdata(jj,kk)*&
          !!!TRIM                subm_trimmer(idomain)%mdata(jj,ii)*&
          !!!TRIM                subm_trimmer(idomain)%mdata(kk,ii)
          !!!TRIM          ENDDO
          !!!TRIM       ENDDO
          !!!TRIM       ! invert the main matrix (exclude some eigenvalues, shift some)
          !!!TRIM       CALL pseudo_invert_matrix(A=Mstore,Ainv=Minv,N=naos,method=1,&
          !!!TRIM               !range1_thr=1.0E-9_dp,range2_thr=1.0E-9_dp,&
          !!!TRIM               shift=1.0E-5_dp,&
          !!!TRIM               range1=nmos(idomain),range2=nmos(idomain),&
          !!!TRIM               error=error)
          !!!TRIM       ! apply the inverted matrix
          !!!TRIM       ! RZK-warning this is only possible when the preconditioner is applied
          !!!TRIM       tmp(:,ii)=MATMUL(Minv,subm_in(idomain)%mdata(:,ii))
          !!!TRIM    ENDDO
          !!!TRIM    subm_out=MATMUL(tmp,sigma)
          !!!TRIM    deallocate(tmp)
          !!!TRIM ELSE
              CALL pseudo_invert_matrix(A=subm_main(idomain)%mdata,Ainv=Minv,N=naos,method=1,&
                      range1=nmos(idomain),range2=n_domain_mos,error=error)
          !!!TRIM ENDIF
   
          CALL copy_submatrices(subm_main(idomain),preconditioner(idomain),.FALSE.,error)
          CALL copy_submatrix_data(Minv,preconditioner(idomain),error)
          
          DEALLOCATE(Minv)

       ENDIF ! submatrix for the domain exists

    ENDDO ! loop over domains

    CALL release_submatrices(subm_main,error)
    DEALLOCATE(subm_main)
    !DEALLOCATE(subm_s)
    !DEALLOCATE(subm_r)

    !IF (matrix_r_required) THEN
    !   CALL cp_dbcsr_release(m_tmp_no_1,error=error)
    !   CALL cp_dbcsr_release(m_tmp_no_2,error=error)
    !   CALL cp_dbcsr_release(matrix_r,error=error)
    !ENDIF

    !RZK-warning do we need a barrier here ?

    CALL timestop(handle)
  
  END SUBROUTINE construct_domain_preconditioner

! *****************************************************************************
!> \brief Constructs S^(+1/2) and S^(-1/2) submatrices for each domain
!> \param matrix_s ...
!> \param subm_s_sqrt ...
!> \param subm_s_sqrt_inv ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_domain_s_sqrt(matrix_s,subm_s_sqrt,subm_s_sqrt_inv,&
    dpattern,map,node_of_domain,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: subm_s_sqrt, subm_s_sqrt_inv
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_s_sqrt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, naos, &
                                                ndomains
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Ssqrt, Ssqrtinv
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_s

    CALL timeset(routineN,handle)
    
    ndomains = cp_dbcsr_nblkcols_total(dpattern)
    CPPrecondition(SIZE(subm_s_sqrt).eq.ndomains,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(subm_s_sqrt_inv).eq.ndomains,cp_failure_level,routineP,error,failure)
    ALLOCATE(subm_s(ndomains))
    CALL init_submatrices(subm_s,error)

    CALL construct_submatrices(matrix_s,subm_s,&
            dpattern,map,node_of_domain,select_row_col,error)

    ! loop over domains - perform inversion
    DO idomain = 1, ndomains
    
       ! check if the submatrix exists
       IF (subm_s(idomain)%domain.gt.0) THEN

          naos=subm_s(idomain)%nrows

          ALLOCATE(Ssqrt(naos,naos))
          ALLOCATE(Ssqrtinv(naos,naos))

          CALL matrix_sqrt(A=subm_s(idomain)%mdata,Asqrt=Ssqrt,Asqrtinv=Ssqrtinv,&
                  N=naos,error=error)
   
          CALL copy_submatrices(subm_s(idomain),subm_s_sqrt(idomain),.FALSE.,error)
          CALL copy_submatrix_data(Ssqrt,subm_s_sqrt(idomain),error)
          
          CALL copy_submatrices(subm_s(idomain),subm_s_sqrt_inv(idomain),.FALSE.,error)
          CALL copy_submatrix_data(Ssqrtinv,subm_s_sqrt_inv(idomain),error)
          
          DEALLOCATE(Ssqrtinv)
          DEALLOCATE(Ssqrt)

       ENDIF ! submatrix for the domain exists

    ENDDO ! loop over domains

    CALL release_submatrices(subm_s,error)
    DEALLOCATE(subm_s)

    CALL timestop(handle)
  
  END SUBROUTINE construct_domain_s_sqrt

! *****************************************************************************
!> \brief Constructs S_inv block for each domain
!> \param matrix_s ...
!> \param subm_s_inv ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param error ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_domain_s_inv(matrix_s,subm_s_inv,dpattern,map,&
    node_of_domain,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: subm_s_inv
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_s_inv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, naos, &
                                                ndomains
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Sinv
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_s

    CALL timeset(routineN,handle)
    
    ndomains = cp_dbcsr_nblkcols_total(dpattern)

    CPPrecondition(SIZE(subm_s_inv).eq.ndomains,cp_failure_level,routineP,error,failure)
    ALLOCATE(subm_s(ndomains))
    CALL init_submatrices(subm_s,error)

    CALL construct_submatrices(matrix_s,subm_s,&
            dpattern,map,node_of_domain,select_row_col,error)

    !GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
    !   cp_dbcsr_distribution(dpattern)))
    !CALL print_submatrices(domain_r_down(:,ispin),GroupID,error)

    ! loop over domains - perform inversion
    DO idomain = 1, ndomains
    
       ! check if the submatrix exists
       IF (subm_s(idomain)%domain.gt.0) THEN

          naos=subm_s(idomain)%nrows

          ALLOCATE(Sinv(naos,naos))

          CALL pseudo_invert_matrix(A=subm_s(idomain)%mdata,Ainv=Sinv,N=naos,&
                  method=0,error=error)
   
          CALL copy_submatrices(subm_s(idomain),subm_s_inv(idomain),.FALSE.,error)
          CALL copy_submatrix_data(Sinv,subm_s_inv(idomain),error)
          
          DEALLOCATE(Sinv)

       ENDIF ! submatrix for the domain exists

    ENDDO ! loop over domains

    CALL release_submatrices(subm_s,error)
    DEALLOCATE(subm_s)

    CALL timestop(handle)
  
  END SUBROUTINE construct_domain_s_inv

! *****************************************************************************
!> \brief Constructs subblocks of the covariant-covariant DM
!> \param matrix_t ...
!> \param matrix_sigma_inv ...
!> \param matrix_s ...
!> \param subm_r_down ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param filter_eps ...
!> \param error ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_domain_r_down(matrix_t,matrix_sigma_inv,matrix_s,&
    subm_r_down,dpattern,map,node_of_domain,filter_eps,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_t, matrix_sigma_inv, &
                                                matrix_s
    TYPE(domain_submatrix_type), &
      DIMENSION(:), INTENT(INOUT)            :: subm_r_down
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_r_down', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndomains
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: m_tmp_no_1, m_tmp_no_2, &
                                                matrix_r

    CALL timeset(routineN,handle)
    
    ! compute the density matrix in the COVARIANT representation
    CALL cp_dbcsr_init(matrix_r,error=error)
    CALL cp_dbcsr_create(matrix_r,&
            template=matrix_s,&
            matrix_type=dbcsr_type_symmetric,error=error)
    CALL cp_dbcsr_init(m_tmp_no_1,error=error)
    CALL cp_dbcsr_create(m_tmp_no_1,&
            template=matrix_t,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_init(m_tmp_no_2,error=error)
    CALL cp_dbcsr_create(m_tmp_no_2,&
            template=matrix_t,&
            matrix_type=dbcsr_type_no_symmetry,error=error)

    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_s, matrix_t,&
            0.0_dp, m_tmp_no_1, filter_eps=filter_eps, error=error)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, m_tmp_no_1, matrix_sigma_inv,&
            0.0_dp, m_tmp_no_2, filter_eps=filter_eps, error=error)
    CALL cp_dbcsr_multiply("N", "T", 1.0_dp, m_tmp_no_2, m_tmp_no_1,&
            0.0_dp, matrix_r, filter_eps=filter_eps, error=error)

    CALL cp_dbcsr_release(m_tmp_no_1,error=error)
    CALL cp_dbcsr_release(m_tmp_no_2,error=error)

    ndomains = cp_dbcsr_nblkcols_total(dpattern)
    CPPrecondition(SIZE(subm_r_down).eq.ndomains,cp_failure_level,routineP,error,failure)

    CALL construct_submatrices(matrix_r,subm_r_down,&
            dpattern,map,node_of_domain,select_row_col,error)
                   
    CALL cp_dbcsr_release(matrix_r,error=error)

    CALL timestop(handle)
  
  END SUBROUTINE construct_domain_r_down

! *****************************************************************************
!> \brief Finds the square root of a matrix and its inverse
!> \param A ...
!> \param Asqrt ...
!> \param Asqrtinv ...
!> \param N ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE matrix_sqrt(A,Asqrt,Asqrtinv,N,error)
    
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: A
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: Asqrt, Asqrtinv
    INTEGER, INTENT(IN)                      :: N
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_sqrt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, INFO, jj, LWORK
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, WORK
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: test, testN

    CALL timeset(routineN,handle)

    Asqrtinv=A
    INFO=0

    !CALL DPOTRF('L', N, Ainv, N, INFO )
    !IF( INFO.NE.0 ) THEN
    !   CPErrorMessage(cp_failure_level,routineP,"DPOTRF failed",error)
    !   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    !END IF
    !CALL DPOTRI('L', N, Ainv, N, INFO )
    !IF( INFO.NE.0 ) THEN
    !   CPErrorMessage(cp_failure_level,routineP,"DPOTRI failed",error)
    !   CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    !END IF
    !! complete the matrix
    !DO ii=1,N
    !   DO jj=ii+1,N
    !      Ainv(ii,jj)=Ainv(jj,ii)
    !   ENDDO
    !   !WRITE(*,'(100F13.9)') Ainv(ii,:)
    !ENDDO

    ! diagonalize first
    ALLOCATE(eigenvalues(N))
    ! Query the optimal workspace for dsyev
    LWORK = -1
    ALLOCATE(WORK(MAX(1,LWORK)))
    CALL DSYEV('V','L',N,Asqrtinv,N,eigenvalues,WORK,LWORK,INFO)
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ! Allocate the workspace and solve the eigenproblem
    ALLOCATE(WORK(MAX(1,LWORK)))
    CALL DSYEV('V','L',N,Asqrtinv,N,eigenvalues,WORK,LWORK,INFO)
    IF ( INFO.NE.0 ) THEN
       WRITE(*,*) 'DSYEV ERROR MESSAGE: ', INFO
       CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(WORK)
   
    ! take functions of eigenvalues and use eigenvectors to compute the matrix function 
    ! first sqrt
    ALLOCATE(test(N,N))
    DO jj=1, N
       test(jj,:)=Asqrtinv(:,jj)*SQRT(eigenvalues(jj))
    ENDDO
    ALLOCATE(testN(N,N))
    testN(:,:)=MATMUL(Asqrtinv,test)
    Asqrt=testN
    ! now, sqrt_inv
    DO jj=1, N
       test(jj,:)=Asqrtinv(:,jj)/SQRT(eigenvalues(jj))
    ENDDO
    testN(:,:)=MATMUL(Asqrtinv,test)
    Asqrtinv=testN
    DEALLOCATE(test,testN)
    
    DEALLOCATE(eigenvalues)

!!!    ! compute the error 
!!!    allocate(test(N,N))
!!!    test=MATMUL(Ainv,A)
!!!    DO ii=1,N
!!!       test(ii,ii)=test(ii,ii)-1.0_dp
!!!    ENDDO
!!!    test_error=0.0_dp
!!!    DO ii=1,N
!!!       DO jj=1,N
!!!          test_error=test_error+test(jj,ii)*test(jj,ii)
!!!       ENDDO
!!!    ENDDO
!!!    WRITE(*,*) "Inversion error: ", SQRT(test_error)
!!!    deallocate(test)

    CALL timestop(handle)
  
  END SUBROUTINE matrix_sqrt

! *****************************************************************************
!> \brief Inverts a matrix using a requested method
!>         0. Cholesky factorization
!>         1. Diagonalization
!> \param A ...
!> \param Ainv ...
!> \param N ...
!> \param method ...
!> \param range1 ...
!> \param range2 ...
!> \param range1_thr ...
!> \param range2_thr ...
!> \param shift ...
!> \param error ...
!> \par History
!>       2012.04 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE pseudo_invert_matrix(A,Ainv,N,method,range1,range2,range1_thr,range2_thr,&
    shift,error)
    
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: A
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: Ainv
    INTEGER, INTENT(IN)                      :: N, method
    INTEGER, INTENT(IN), OPTIONAL            :: range1, range2
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: range1_thr, range2_thr, shift
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pseudo_invert_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ii, INFO, jj, LWORK, &
                                                range1_eiv, range2_eiv, &
                                                range3_eiv
    LOGICAL                                  :: failure, use_ranges
    REAL(KIND=dp)                            :: my_shift
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, WORK
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: test, testN

    CALL timeset(routineN,handle)

    IF (method.eq.1) THEN
       IF (PRESENT(range1)) THEN
          use_ranges=.TRUE.
          IF (.NOT.PRESENT(range2)) THEN
             CPErrorMessage(cp_failure_level,routineP,"SPECIFY TWO RANGES",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ENDIF
       ELSE
          use_ranges=.FALSE.
          IF ((.NOT.PRESENT(range1_thr)).OR.(.NOT.PRESENT(range2_thr))) THEN
             CPErrorMessage(cp_failure_level,routineP,"SPECIFY TWO THRESHOLDS",error)
             CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          ENDIF
       ENDIF
    ENDIF

    my_shift=0.0_dp
    IF (PRESENT(shift)) THEN
       my_shift=shift
    ENDIF

    Ainv=A
    INFO=0

    SELECT CASE (method)
    CASE (0)

       CALL DPOTRF('L', N, Ainv, N, INFO )
       IF( INFO.NE.0 ) THEN
          CPErrorMessage(cp_failure_level,routineP,"DPOTRF failed",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       CALL DPOTRI('L', N, Ainv, N, INFO )
       IF( INFO.NE.0 ) THEN
          CPErrorMessage(cp_failure_level,routineP,"DPOTRI failed",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       ! complete the matrix
       DO ii=1,N
          DO jj=ii+1,N
             Ainv(ii,jj)=Ainv(jj,ii)
          ENDDO
          !WRITE(*,'(100F13.9)') Ainv(ii,:)
       ENDDO

    CASE (1)
    
       ! diagonalize first
       ALLOCATE(eigenvalues(N))
       ! Query the optimal workspace for dsyev
       LWORK = -1
       ALLOCATE(WORK(MAX(1,LWORK)))
       CALL DSYEV('V','L',N,Ainv,N,eigenvalues,WORK,LWORK,INFO)
       LWORK = INT(WORK(1))
       DEALLOCATE(WORK)
       ! Allocate the workspace and solve the eigenproblem
       ALLOCATE(WORK(MAX(1,LWORK)))
       CALL DSYEV('V','L',N,Ainv,N,eigenvalues,WORK,LWORK,INFO)
       IF ( INFO.NE.0 ) THEN
          WRITE(*,*) 'DSYEV ERROR MESSAGE: ', INFO
          CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(WORK)

       !WRITE(*,*) "EIGENVALS: "
       !WRITE(*,'(4F13.9)') eigenvalues(:)
   
       ! invert eigenvalues and use eigenvectors to compute the Hessian inverse
       ! project out zero-eigenvalue directions
       ALLOCATE(test(N,N))
       range1_eiv=0
       range2_eiv=0
       range3_eiv=0
       IF (use_ranges) THEN
          DO jj=1,N 
             IF (jj.le.range1) THEN
                test(jj,:)=Ainv(:,jj)*0.0_dp
                range1_eiv=range1_eiv+1
             ELSE IF (jj.le.range2) THEN
                test(jj,:)=Ainv(:,jj)*1.0_dp
                range2_eiv=range2_eiv+1
             ELSE
                test(jj,:)=Ainv(:,jj)/(eigenvalues(jj)+my_shift)
                range3_eiv=range3_eiv+1
             ENDIF
          ENDDO
       ELSE
          DO jj=1, N
             IF (eigenvalues(jj).lt.range1_thr) THEN
                test(jj,:)=Ainv(:,jj)*0.0_dp
                range1_eiv=range1_eiv+1
             ELSE IF (eigenvalues(jj).lt.range2_thr) THEN
                test(jj,:)=Ainv(:,jj)*1.0_dp
                range2_eiv=range2_eiv+1
             ELSE
                test(jj,:)=Ainv(:,jj)/(eigenvalues(jj)+my_shift)
                range3_eiv=range3_eiv+1
             ENDIF
          ENDDO
       ENDIF
       !WRITE(*,*) ' EIV RANGES: ', range1_eiv, range2_eiv, range3_eiv
       ALLOCATE(testN(N,N))
       testN(:,:)=MATMUL(Ainv,test)
       Ainv=testN
       DEALLOCATE(test,testN)
       DEALLOCATE(eigenvalues)

    CASE DEFAULT

       CPErrorMessage(cp_failure_level,routineP,"Illegal method selected for matrix inversion",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)

    END SELECT

!!!    ! compute the inversion error 
!!!    allocate(test(N,N))
!!!    test=MATMUL(Ainv,A)
!!!    DO ii=1,N
!!!       test(ii,ii)=test(ii,ii)-1.0_dp
!!!    ENDDO
!!!    test_error=0.0_dp
!!!    DO ii=1,N
!!!       DO jj=1,N
!!!          test_error=test_error+test(jj,ii)*test(jj,ii)
!!!       ENDDO
!!!    ENDDO
!!!    WRITE(*,*) "Inversion error: ", SQRT(test_error)
!!!    deallocate(test)

    CALL timestop(handle)
  
  END SUBROUTINE pseudo_invert_matrix

! *****************************************************************************
!> \brief computes the step matrix from the gradient and Hessian using
!>         the Newton-Raphson method
!> \param matrix_grad ...
!> \param matrix_step ...
!> \param matrix_s ...
!> \param matrix_ks ...
!> \param matrix_t ...
!> \param matrix_sigma_inv ...
!> \param quench_t ...
!> \param spin_factor ...
!> \param eps_filter ...
!> \param error ...
!> \par History
!>       2012.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE newton_grad_to_step(matrix_grad,matrix_step,matrix_s,matrix_ks,&
    matrix_t,matrix_sigma_inv,quench_t,spin_factor,eps_filter,error)
    
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_grad, matrix_step, &
                                                matrix_s
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_ks, matrix_t
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_sigma_inv, quench_t
    REAL(KIND=dp), INTENT(IN)                :: spin_factor, eps_filter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'newton_grad_to_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: res_norm
    TYPE(cp_dbcsr_type) :: m_tmp_no_1, m_tmp_no_3, m_tmp_oo_2, matrix_f_ao, &
      matrix_f_mo, matrix_f_vo, matrix_s_ao, matrix_s_mo, matrix_s_vo

    CALL timeset(routineN,handle)
    
    CALL cp_dbcsr_init(matrix_s_ao,error=error)
    CALL cp_dbcsr_create(matrix_s_ao,&
            template=matrix_s,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)
    CALL cp_dbcsr_init(matrix_f_ao,error=error)
    CALL cp_dbcsr_create(matrix_f_ao,&
            template=matrix_s,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)
    CALL cp_dbcsr_init(matrix_f_mo,error=error)
    CALL cp_dbcsr_create(matrix_f_mo,&
            template=matrix_sigma_inv,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)
    CALL cp_dbcsr_init(matrix_s_mo,error=error)
    CALL cp_dbcsr_create(matrix_s_mo,&
            template=matrix_sigma_inv,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)
    CALL cp_dbcsr_init(matrix_f_vo,error=error)
    CALL cp_dbcsr_create(matrix_f_vo,&
            template=matrix_t,error=error)
    CALL cp_dbcsr_init(matrix_s_vo,error=error)
    CALL cp_dbcsr_create(matrix_s_vo,&
            template=matrix_t,error=error)

    CALL cp_dbcsr_init(m_tmp_no_1,error=error)
    CALL cp_dbcsr_create(m_tmp_no_1,&
            template=matrix_t,error=error)
    CALL cp_dbcsr_init(m_tmp_no_3,error=error)
    CALL cp_dbcsr_create(m_tmp_no_3,&
            template=matrix_t,error=error)
    CALL cp_dbcsr_init(m_tmp_oo_2,error=error)
    CALL cp_dbcsr_create(m_tmp_oo_2,&
            template=matrix_sigma_inv,&
            matrix_type=dbcsr_type_no_symmetry,&
            error=error)
    
    ! calculate S-SRS and (1-R)F(1-R)
    ! RZK-warning some optimization is ABSOLUTELY NECESSARY
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_s,&
            matrix_t,&
            0.0_dp,m_tmp_no_1,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            m_tmp_no_1,&
            matrix_sigma_inv,&
            0.0_dp,matrix_s_vo,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_desymmetrize(matrix_s,&
            matrix_s_ao,error=error)
    CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
            m_tmp_no_1,&
            matrix_s_vo,&
            1.0_dp,matrix_s_ao,&
            filter_eps=eps_filter,&
            error=error)
    
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_ks,&
            matrix_t,&
            0.0_dp,m_tmp_no_1,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_desymmetrize(matrix_ks,matrix_f_ao,error=error)
    CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
            m_tmp_no_1,&
            matrix_s_vo,&
            1.0_dp,matrix_f_ao,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
            matrix_s_vo,&
            m_tmp_no_1,&
            1.0_dp,matrix_f_ao,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
            matrix_t,&
            m_tmp_no_1,&
            0.0_dp,matrix_f_mo,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_s_vo,&
            matrix_f_mo,&
            0.0_dp,m_tmp_no_1,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","T",1.0_dp,&
            m_tmp_no_1,&
            matrix_s_vo,&
            1.0_dp,matrix_f_ao,&
            filter_eps=eps_filter,&
            error=error)

    ! calculate F_mo
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_sigma_inv,&
            matrix_f_mo,&
            0.0_dp,m_tmp_oo_2,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            m_tmp_oo_2,&
            matrix_sigma_inv,&
            1.0_dp,matrix_f_mo,&
            filter_eps=eps_filter,&
            error=error)

    ! calculate F_vo
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_ks,&
            matrix_t,&
            0.0_dp,m_tmp_no_1,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            m_tmp_no_1,&
            matrix_sigma_inv,&
            0.0_dp,matrix_f_vo,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","T",-1.0_dp,&
            matrix_s_vo,&
            m_tmp_oo_2,&
            1.0_dp,matrix_f_vo,&
            filter_eps=eps_filter,&
            error=error)

    CALL cp_dbcsr_desymmetrize(matrix_sigma_inv,matrix_s_mo,error=error)

    !!! RZK-warning: this is HIGHLIGHTED BLOCK
    !!! check it first if the procedure does not function as it supposed to
    !CALL cp_dbcsr_desymmetrize(matrix_s,matrix_s_ao,error=error)
    CALL cp_dbcsr_add(matrix_f_ao,matrix_s_ao,1.0_dp,1.0_dp,error=error)
    
    CALL cp_dbcsr_set(matrix_f_mo,0.0_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_f_mo,1.0_dp,error=error)
    CALL cp_dbcsr_filter(matrix_f_mo,eps_filter,error=error)
    
    CALL cp_dbcsr_set(matrix_s_mo,0.0_dp,error=error)
    CALL cp_dbcsr_add_on_diag(matrix_s_mo,1.0_dp,error=error)
    CALL cp_dbcsr_filter(matrix_s_mo,eps_filter,error=error)
    
    CALL cp_dbcsr_set(matrix_s_ao,0.0_dp,error=error)
    !CALL cp_dbcsr_add_on_diag(matrix_s_ao,1.0_dp,error=error)
    CALL cp_dbcsr_filter(matrix_s_ao,eps_filter,error=error)
    !!! RZK-warning: end of HIGHLIGHTED BLOCK

    CALL cp_dbcsr_scale(matrix_f_ao,&
            2.0_dp*spin_factor,&
            error=error)
    CALL cp_dbcsr_scale(matrix_s_ao,&
            -2.0_dp*spin_factor,&
            error=error)
    CALL cp_dbcsr_scale(matrix_f_vo,&
            2.0_dp*spin_factor,&
            error=error)

    !WRITE(*,*) "INSIDE newton_grad_to_step: "
    !CALL cp_dbcsr_print(matrix_s_mo,error=error)
    !CALL cp_dbcsr_print(matrix_ks,error=error)
    !CALL cp_dbcsr_print(matrix_s,error=error)
    !CALL cp_dbcsr_print(matrix_sigma_inv,error=error)

    CALL hessian_diag_apply(matrix_grad,matrix_step,matrix_s_ao,&
            matrix_f_ao,matrix_s_mo,matrix_f_mo,&
            matrix_s_vo,matrix_f_vo,quench_t,error)

    ! check that the step satisfies H.step=-grad
    CALL cp_dbcsr_copy(m_tmp_no_3,quench_t,error=error)
    CALL cp_dbcsr_copy(m_tmp_no_1,quench_t,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_f_ao,&
            matrix_step,&
            0.0_dp,m_tmp_no_1,&
            !retain_sparsity=.TRUE.,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            m_tmp_no_1,&
            matrix_s_mo,&
            0.0_dp,m_tmp_no_3,&
            retain_sparsity=.TRUE.,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_copy(m_tmp_no_1,quench_t,error=error)
    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            matrix_s_ao,&
            matrix_step,&
            0.0_dp,m_tmp_no_1,&
            !retain_sparsity=.TRUE.,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_multiply("N","N",-1.0_dp,&
            m_tmp_no_1,&
            matrix_f_mo,&
            1.0_dp,m_tmp_no_3,&
            retain_sparsity=.TRUE.,&
            filter_eps=eps_filter,&
            error=error)
    CALL cp_dbcsr_add(m_tmp_no_3,matrix_grad,&
            1.0_dp,1.0_dp,error=error)
    CALL cp_dbcsr_norm(m_tmp_no_3,&
            dbcsr_norm_maxabsnorm, norm_scalar=res_norm, error=error)
    WRITE(*,*) "NEWTON step error: ", res_norm

    CALL cp_dbcsr_release(m_tmp_no_3,error=error)
    CALL cp_dbcsr_release(m_tmp_no_1,error=error)
    CALL cp_dbcsr_release(m_tmp_oo_2,error=error)
    CALL cp_dbcsr_release(matrix_s_ao,error=error)
    CALL cp_dbcsr_release(matrix_s_mo,error=error)
    CALL cp_dbcsr_release(matrix_f_ao,error=error)
    CALL cp_dbcsr_release(matrix_f_mo,error=error)
    CALL cp_dbcsr_release(matrix_s_vo,error=error)
    CALL cp_dbcsr_release(matrix_f_vo,error=error)
    
    CALL timestop(handle)
  
  END SUBROUTINE newton_grad_to_step

! *****************************************************************************
!> \brief Serial code that constructs an approximate Hessian
!> \param matrix_grad ...
!> \param matrix_step ...
!> \param matrix_S_ao ...
!> \param matrix_F_ao ...
!> \param matrix_S_mo ...
!> \param matrix_F_mo ...
!> \param matrix_S_vo ...
!> \param matrix_F_vo ...
!> \param quench_t ...
!> \param error ...
!> \par History
!>       2012.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE hessian_diag_apply(matrix_grad,matrix_step,matrix_S_ao,&
    matrix_F_ao,matrix_S_mo,matrix_F_mo,matrix_S_vo,matrix_F_vo,quench_t,error)
    
    TYPE(cp_dbcsr_type), INTENT(INOUT) :: matrix_grad, matrix_step, &
      matrix_S_ao, matrix_F_ao, matrix_S_mo, matrix_F_mo, matrix_S_vo, &
      matrix_F_vo, quench_t
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hessian_diag_apply', &
      routineP = moduleN//':'//routineN

    INTEGER :: ao_hori_offset, ao_vert_offset, block_col, block_row, col, &
      copy, H_size, handle, ii, INFO, jj, lev1_hori_offset, lev1_vert_offset, &
      lev2_hori_offset, lev2_vert_offset, LWORK, nblkcols_tot, nblkrows_tot, &
      orb_i, orb_j, row, zero_neg_eiv
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ao_block_sizes, &
                                                ao_domain_sizes, &
                                                mo_block_sizes
    INTEGER, DIMENSION(:), POINTER           :: ao_blk_sizes, mo_blk_sizes
    LOGICAL                                  :: failure, found, found2, &
                                                found_col, found_row
    REAL(KIND=dp)                            :: test_error
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, Grad_vec, &
                                                Step_vec, tmp, tmpr, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: F_ao_block, F_mo_block, H, &
                                                H1, H2, Hinv, S_ao_block, &
                                                S_mo_block, test, test2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: block_p, block_p2, p_new_block
    TYPE(cp_dbcsr_type)                      :: matrix_F_ao_sym, &
                                                matrix_F_mo_sym, &
                                                matrix_S_ao_sym, &
                                                matrix_S_mo_sym

    CALL timeset(routineN,handle)

    nblkrows_tot = cp_dbcsr_nblkrows_total(quench_t)
    nblkcols_tot = cp_dbcsr_nblkcols_total(quench_t)
    CPPostcondition(nblkrows_tot==nblkcols_tot,cp_failure_level,routineP,error,failure)
    mo_blk_sizes => array_data(cp_dbcsr_col_block_sizes(quench_t))
    ao_blk_sizes => array_data(cp_dbcsr_row_block_sizes(quench_t))
    ALLOCATE(mo_block_sizes(nblkcols_tot),ao_block_sizes(nblkcols_tot))
    ALLOCATE(ao_domain_sizes(nblkcols_tot))
    mo_block_sizes(:)=mo_blk_sizes(:)
    ao_block_sizes(:)=ao_blk_sizes(:)
    ao_domain_sizes(:)=0

    CALL cp_dbcsr_init(matrix_S_ao_sym,error=error)
    CALL cp_dbcsr_create(matrix_S_ao_sym,&
            template=matrix_S_ao,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_desymmetrize(matrix_S_ao,matrix_S_ao_sym,error=error)

    CALL cp_dbcsr_init(matrix_F_ao_sym,error=error)
    CALL cp_dbcsr_create(matrix_F_ao_sym,&
            template=matrix_F_ao,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_desymmetrize(matrix_F_ao,matrix_F_ao_sym,error=error)

    CALL cp_dbcsr_init(matrix_S_mo_sym,error=error)
    CALL cp_dbcsr_create(matrix_S_mo_sym,&
            template=matrix_S_mo,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_desymmetrize(matrix_S_mo,matrix_S_mo_sym,error=error)

    CALL cp_dbcsr_init(matrix_F_mo_sym,error=error)
    CALL cp_dbcsr_create(matrix_F_mo_sym,&
            template=matrix_F_mo,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_desymmetrize(matrix_F_mo,matrix_F_mo_sym,error=error)

    !CALL cp_dbcsr_print(matrix_grad,error=error)
    !CALL cp_dbcsr_print(matrix_F_ao_sym,error=error)
    !CALL cp_dbcsr_print(matrix_S_ao_sym,error=error)
    !CALL cp_dbcsr_print(matrix_F_mo_sym,error=error)
    !CALL cp_dbcsr_print(matrix_S_mo_sym,error=error)

    ! loop over domains to find the size of the Hessian
    H_size=0
    DO col = 1, nblkcols_tot

       ! find sizes of AO submatrices
       DO row = 1, nblkrows_tot
          
          CALL cp_dbcsr_get_block_p(quench_t,&
                  row, col, block_p, found)
          IF (found) THEN
            ao_domain_sizes(col)=ao_domain_sizes(col)+ao_blk_sizes(row)
          ENDIF

       ENDDO

       H_size=H_size+ao_domain_sizes(col)*mo_block_sizes(col)

    ENDDO

    ALLOCATE(H(H_size,H_size))

    ! fill the Hessian matrix
    lev1_vert_offset=0
    ! loop over all pairs of fragments
    DO row = 1, nblkcols_tot

       lev1_hori_offset=0
       DO col = 1, nblkcols_tot

          ! prepare blocks for the current row-column fragment pair
          ALLOCATE(F_ao_block(ao_domain_sizes(row),ao_domain_sizes(col)))
          ALLOCATE(S_ao_block(ao_domain_sizes(row),ao_domain_sizes(col)))
          ALLOCATE(F_mo_block(mo_block_sizes(row),mo_block_sizes(col)))
          ALLOCATE(S_mo_block(mo_block_sizes(row),mo_block_sizes(col)))
   
          F_ao_block(:,:)=0.0_dp
          S_ao_block(:,:)=0.0_dp
          F_mo_block(:,:)=0.0_dp
          S_mo_block(:,:)=0.0_dp
   
          ! fill AO submatrices
          ! loop over all blocks of the AO dbcsr matrix
          ao_vert_offset=0
          DO block_row = 1, nblkcols_tot
             
             CALL cp_dbcsr_get_block_p(quench_t,&
                     block_row, row, block_p, found_row)
             IF (found_row) THEN
      
                ao_hori_offset=0
                DO block_col = 1, nblkcols_tot
                   
                   CALL cp_dbcsr_get_block_p(quench_t,&
                           block_col, col, block_p, found_col)
                   IF (found_col) THEN
         
                      CALL cp_dbcsr_get_block_p(matrix_F_ao_sym,&
                              block_row, block_col, block_p, found)
                      IF (found) THEN
                         ! copy the block into the submatrix
                         F_ao_block(ao_vert_offset+1:ao_vert_offset+ao_block_sizes(block_row),&
                           ao_hori_offset+1:ao_hori_offset+ao_block_sizes(block_col))&
                           =block_p(:,:)
                      ENDIF
         
                      CALL cp_dbcsr_get_block_p(matrix_S_ao_sym,&
                              block_row, block_col, block_p, found)
                      IF (found) THEN
                         ! copy the block into the submatrix
                         S_ao_block(ao_vert_offset+1:ao_vert_offset+ao_block_sizes(block_row),&
                           ao_hori_offset+1:ao_hori_offset+ao_block_sizes(block_col))&
                           =block_p(:,:)
                      ENDIF
         
                      ao_hori_offset=ao_hori_offset+ao_block_sizes(block_col)
                   
                   ENDIF
      
                ENDDO

                ao_vert_offset=ao_vert_offset+ao_block_sizes(block_row)

             ENDIF
   
          ENDDO

          ! fill MO submatrices
          CALL cp_dbcsr_get_block_p(matrix_F_mo_sym, row, col, block_p, found)
          IF (found) THEN
             ! copy the block into the submatrix
             F_mo_block(1:mo_block_sizes(row),1:mo_block_sizes(col))=block_p(:,:)
          ENDIF
          CALL cp_dbcsr_get_block_p(matrix_S_mo_sym, row, col, block_p, found)
          IF (found) THEN
             ! copy the block into the submatrix
             S_mo_block(1:mo_block_sizes(row),1:mo_block_sizes(col))=block_p(:,:)
          ENDIF

          !WRITE(*,*) "F_AO_BLOCK", row, col, ao_domain_sizes(row), ao_domain_sizes(col)
          !DO ii=1,ao_domain_sizes(row)
          !  WRITE(*,'(100F13.9)') F_ao_block(ii,:)
          !ENDDO
          !WRITE(*,*) "S_AO_BLOCK", row, col
          !DO ii=1,ao_domain_sizes(row)
          !  WRITE(*,'(100F13.9)') S_ao_block(ii,:)
          !ENDDO
          !WRITE(*,*) "F_MO_BLOCK", row, col
          !DO ii=1,mo_block_sizes(row)
          !  WRITE(*,'(100F13.9)') F_mo_block(ii,:)
          !ENDDO
          !WRITE(*,*) "S_MO_BLOCK", row, col, mo_block_sizes(row), mo_block_sizes(col)
          !DO ii=1,mo_block_sizes(row)
          !  WRITE(*,'(100F13.9)') S_mo_block(ii,:)
          !ENDDO

          ! construct tensor products for the current row-column fragment pair
          lev2_vert_offset=0
          DO orb_j=1,mo_block_sizes(row)

             lev2_hori_offset=0
             DO orb_i=1,mo_block_sizes(col)

                H(lev1_vert_offset+lev2_vert_offset+1:lev1_vert_offset+lev2_vert_offset+ao_domain_sizes(row),&
                  lev1_hori_offset+lev2_hori_offset+1:lev1_hori_offset+lev2_hori_offset+ao_domain_sizes(col))&
               =S_mo_block(orb_j,orb_i)*F_ao_block(:,:)&
               -F_mo_block(orb_j,orb_i)*S_ao_block(:,:)

                !WRITE(*,*) row, col, orb_j, orb_i, lev1_vert_offset+lev2_vert_offset+1, ao_domain_sizes(row),&
                !   lev1_hori_offset+lev2_hori_offset+1, ao_domain_sizes(col), S_mo_block(orb_j,orb_i)

                lev2_hori_offset=lev2_hori_offset+ao_domain_sizes(col)

             ENDDO

             lev2_vert_offset=lev2_vert_offset+ao_domain_sizes(row)

          ENDDO

          
          lev1_hori_offset=lev1_hori_offset+ao_domain_sizes(col)*mo_block_sizes(col)

          DEALLOCATE(F_ao_block)
          DEALLOCATE(S_ao_block)
          DEALLOCATE(F_mo_block)
          DEALLOCATE(S_mo_block)

       ENDDO ! col fragment

       lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(row)*mo_block_sizes(row)

    ENDDO ! row fragment

    CALL cp_dbcsr_release(matrix_S_ao_sym,error=error)
    CALL cp_dbcsr_release(matrix_F_ao_sym,error=error)
    CALL cp_dbcsr_release(matrix_S_mo_sym,error=error)
    CALL cp_dbcsr_release(matrix_F_mo_sym,error=error)
          
    ! two more terms of the Hessian
    ALLOCATE(H1(H_size,H_size))
    ALLOCATE(H2(H_size,H_size))
    H1=0.0_dp
    H2=0.0_dp
    DO row = 1, nblkcols_tot

       lev1_hori_offset=0
       DO col = 1, nblkcols_tot

          CALL cp_dbcsr_get_block_p(matrix_F_vo,&
                  row, col, block_p, found)
          CALL cp_dbcsr_get_block_p(matrix_S_vo,&
                  row, col, block_p2, found2)

          lev1_vert_offset=0
          DO block_col = 1, nblkcols_tot
             
             CALL cp_dbcsr_get_block_p(quench_t,&
                     row, block_col, p_new_block, found_row)

             IF (found_row) THEN

                ! determine offset in this short loop
                lev2_vert_offset=0
                DO block_row=1,row-1
                   CALL cp_dbcsr_get_block_p(quench_t,&
                           block_row, block_col, p_new_block, found_col)
                   IF (found_col) lev2_vert_offset=lev2_vert_offset+ao_block_sizes(block_row)
                ENDDO
                !!!!!!!! short loop

                ! over all electrons of the block
                DO orb_i=1, mo_block_sizes(col)
                
                   ! into all possible locations
                   DO orb_j=1, mo_block_sizes(block_col)
   
                      ! column is copied several times
                      DO copy=1, ao_domain_sizes(col)
   
                         IF (found) THEN

                            !WRITE(*,*) row, col, block_col, orb_i, orb_j, copy,&
                            ! lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1,&
                            ! lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy

                            H1( lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1:&
                                lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+ao_block_sizes(row),&
                                lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy )&
                              =block_p(:,orb_i)
   
                         ENDIF ! found block in the data matrix
       
                         IF (found2) THEN
          
                            H2( lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+1:&
                                lev1_vert_offset+(orb_j-1)*ao_domain_sizes(block_col)+lev2_vert_offset+ao_block_sizes(row),&
                                lev1_hori_offset+(orb_i-1)*ao_domain_sizes(col)+copy )&
                              =block_p2(:,orb_i)
   
                         ENDIF ! found block in the data matrix

                      ENDDO
   
                   ENDDO

                ENDDO
                
                !lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)

             ENDIF ! found block in the quench matrix

             lev1_vert_offset=lev1_vert_offset+&
                ao_domain_sizes(block_col)*mo_block_sizes(block_col)

          ENDDO

          lev1_hori_offset=lev1_hori_offset+&
             ao_domain_sizes(col)*mo_block_sizes(col)

       ENDDO

       !lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)

    ENDDO
    ! add terms to the hessian
!WRITE(*,*) "F_vo"
!DO ii=1,H_size
! WRITE(*,'(100F13.9)') H1(ii,:)
!ENDDO
!WRITE(*,*) "S_vo"
!DO ii=1,H_size
! WRITE(*,'(100F13.9)') H2(ii,:)
!ENDDO
    !DO ii=1,H_size
    !   DO jj=1,H_size
    !      H(ii,jj)=H(ii,jj)-H1(ii,jj)*H2(jj,ii)-H1(jj,ii)*H2(ii,jj)
    !   ENDDO
    !ENDDO
    DEALLOCATE(H1)
    DEALLOCATE(H2)
             
    ! convert gradient from the dbcsr matrix to the vector form
    ALLOCATE(Grad_vec(H_size))
    Grad_vec(:)=0.0_dp
    lev1_vert_offset=0
    ! loop over all electron blocks 
    DO col = 1, nblkcols_tot

       ! loop over AO-rows of the dbcsr matrix
       lev2_vert_offset=0
       DO row = 1, nblkrows_tot
          
          CALL cp_dbcsr_get_block_p(quench_t,&
                  row, col, block_p, found_row)
          IF (found_row) THEN
      
             CALL cp_dbcsr_get_block_p(matrix_grad,&
                     row, col, block_p, found)
             IF (found) THEN
                ! copy the data into the vector, column by column
                DO orb_i=1, mo_block_sizes(col)
                   Grad_vec(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
                            lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))&
                            =block_p(:,orb_i)
!WRITE(*,*) "GRAD: ", row, col, orb_i, lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1, ao_block_sizes(row)
                ENDDO

             ENDIF
       
             lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)

          ENDIF
   
       ENDDO

       lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(col)*mo_block_sizes(col)

    ENDDO ! loop over electron blocks 

!WRITE(*,*) "HESSIAN"
!DO ii=1,H_size
! WRITE(*,'(100F13.9)') H(ii,:)
!ENDDO

    ! invert the Hessian
    INFO=0
    ALLOCATE(Hinv(H_size,H_size))
    Hinv(:,:)=H(:,:)

    ! before inverting diagonalize
    ALLOCATE(eigenvalues(H_size))
    ! Query the optimal workspace for dsyev
    LWORK = -1
    ALLOCATE(WORK(MAX(1,LWORK)))
    CALL DSYEV('V','L',H_size,Hinv,H_size,eigenvalues,WORK,LWORK,INFO)
    LWORK = INT(WORK( 1 ))
    DEALLOCATE(WORK)
    ! Allocate the workspace and solve the eigenproblem
    ALLOCATE(WORK(MAX(1,LWORK)))
    CALL DSYEV('V','L',H_size,Hinv,H_size,eigenvalues,WORK,LWORK,INFO)
    IF( INFO.NE.0 ) THEN
       WRITE(*,*) 'DSYEV ERROR MESSAGE: ', INFO
       CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(WORK)

    ! invert eigenvalues and use eigenvectors to compute the Hessian inverse
    ! project out zero-eigenvalue directions
    ALLOCATE(test(H_size,H_size))
    zero_neg_eiv=0
    DO jj=1, H_size
       IF (eigenvalues(jj).gt.1.0E-8) THEN
          test(jj,:)=Hinv(:,jj)/eigenvalues(jj)
       ELSE
          test(jj,:)=Hinv(:,jj)*0.0_dp
          zero_neg_eiv=zero_neg_eiv+1
       ENDIF
    ENDDO
    WRITE(*,*) 'ZERO OR NEGATIVE EIGENVALUES: ', zero_neg_eiv
    ALLOCATE(test2(H_size,H_size))
    test2(:,:)=MATMUL(Hinv,test)
    Hinv(:,:)=test2(:,:)
    DEALLOCATE(test,test2)

    !! shift to kill singularity
    !shift=0.0_dp
    !IF (eigenvalues(1).lt.0.0_dp) THEN
    !   CPErrorMessage(cp_failure_level,routineP,"Negative eigenvalue(s)",error)
    !   shift=abs(eigenvalues(1))
    !   WRITE(*,*) "Lowest eigenvalue: ", eigenvalues(1)
    !ENDIF
    !DO ii=1, H_size
    !   IF (eigenvalues(ii).gt.1.0E-6_dp) THEN
    !      shift=shift+min(1.0_dp,eigenvalues(ii))*1.0E-4_dp
    !      EXIT
    !   ENDIF
    !ENDDO
    !WRITE(*,*) "Hessian shift: ", shift
    !DO ii=1, H_size
    !   H(ii,ii)=H(ii,ii)+shift
    !ENDDO
    !! end shift
    
    DEALLOCATE(eigenvalues)
    
!!!!    Hinv=H
!!!!    INFO=0
!!!!    CALL DPOTRF('L', H_size, Hinv, H_size, INFO )
!!!!    IF( INFO.NE.0 ) THEN
!!!!       WRITE(*,*) 'DPOTRF ERROR MESSAGE: ', INFO
!!!!       CPErrorMessage(cp_failure_level,routineP,"DPOTRF failed",error)
!!!!       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!!!!    END IF
!!!!    CALL DPOTRI('L', H_size, Hinv, H_size, INFO )
!!!!    IF( INFO.NE.0 ) THEN
!!!!       WRITE(*,*) 'DPOTRI ERROR MESSAGE: ', INFO
!!!!       CPErrorMessage(cp_failure_level,routineP,"DPOTRI failed",error)
!!!!       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
!!!!    END IF
!!!!    ! complete the matrix
!!!!    DO ii=1,H_size
!!!!       DO jj=ii+1,H_size
!!!!          Hinv(ii,jj)=Hinv(jj,ii)
!!!!       ENDDO
!!!!    ENDDO

    ! compute the inversion error 
    ALLOCATE(test(H_size,H_size))
!    WRITE(*,*) "SIZE: ", H_size
    test(:,:)=MATMUL(Hinv,H)
!WRITE(*,*) "TEST"
!DO ii=1,H_size
! WRITE(*,'(100F8.4)') test(ii,:)
!ENDDO
    DO ii=1,H_size
       test(ii,ii)=test(ii,ii)-1.0_dp
    ENDDO
    test_error=0.0_dp
    DO ii=1,H_size
       DO jj=1,H_size
          test_error=test_error+test(jj,ii)*test(jj,ii)
       ENDDO
    ENDDO
    WRITE(*,*) "Hessian inversion error: ", SQRT(test_error)
    DEALLOCATE(test)

    ! prepare the output vector
    ALLOCATE(Step_vec(H_size))
    ALLOCATE(tmp(H_size))
    tmp(:)=MATMUL(Hinv,Grad_vec)
    Step_vec(:)=-1.0_dp*tmp(:)
    
ALLOCATE(tmpr(H_size))
tmpr(:)=MATMUL(H,Step_vec)
tmp(:)=tmpr(:)+Grad_vec(:)
DEALLOCATE(tmpr)
WRITE(*,*) "NEWTOV step error: ", MAXVAL(ABS(tmp))

    DEALLOCATE(tmp)

    DEALLOCATE(H)
    DEALLOCATE(Hinv)
    DEALLOCATE(Grad_vec)

    ! copy the step from the vector into the cp_dbcsr matrix
    
    ! re-create the step matrix to remove all blocks
    CALL cp_dbcsr_create(matrix_step,&
            template=matrix_grad,&
            matrix_type=dbcsr_type_no_symmetry,error=error)
    CALL cp_dbcsr_work_create(matrix_step,work_mutable=.TRUE.,&
            error=error)
    
    lev1_vert_offset=0
    ! loop over all electron blocks 
    DO col = 1, nblkcols_tot

       ! loop over AO-rows of the dbcsr matrix
       lev2_vert_offset=0
       DO row = 1, nblkrows_tot
          
          CALL cp_dbcsr_get_block_p(quench_t,&
                  row, col, block_p, found_row)
          IF (found_row) THEN
      
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(matrix_step,row,col,p_new_block)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             ! copy the data column by column
             DO orb_i=1, mo_block_sizes(col)
                p_new_block(:,orb_i) = &
                   Step_vec(lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1:&
                            lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+ao_block_sizes(row))
!WRITE(*,*) "STEP: ", row, col, orb_i, lev1_vert_offset+ao_domain_sizes(col)*(orb_i-1)+lev2_vert_offset+1, ao_block_sizes(row)
             ENDDO

             lev2_vert_offset=lev2_vert_offset+ao_block_sizes(row)

          ENDIF
   
       ENDDO

       lev1_vert_offset=lev1_vert_offset+ao_domain_sizes(col)*mo_block_sizes(col)

    ENDDO ! loop over electron blocks 

    DEALLOCATE(Step_vec)
    
    CALL cp_dbcsr_finalize(matrix_step,error=error)

    DEALLOCATE(mo_block_sizes,ao_block_sizes)
    DEALLOCATE(ao_domain_sizes)
    
    CALL timestop(handle)
  
  END SUBROUTINE hessian_diag_apply 

! *****************************************************************************
!> \brief Load balancing of the submatrix computations 
!> \param almo_scf_env ...
!> \param error ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE distribute_domains(almo_scf_env,error)
    
    TYPE(almo_scf_env_type), INTENT(INOUT)   :: almo_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribute_domains', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, &
                                                least_loaded, nao, ncpus, &
                                                ndomains
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpu_load, domain_load

    CALL timeset(routineN,handle)

    ndomains = almo_scf_env%ndomains
    ncpus = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
           cp_dbcsr_distribution(almo_scf_env%matrix_s(1))))

    ALLOCATE(domain_load(ndomains))
    DO idomain=1,ndomains
       nao=almo_scf_env%nbasis_of_domain(idomain)
       domain_load(idomain)=(nao*nao*nao)*1.0_dp
    ENDDO

    ALLOCATE(index0(ndomains))
    
    CALL sort(domain_load,ndomains,index0)

    ALLOCATE(cpu_load(ncpus))
    cpu_load(:)=0.0_dp

    DO idomain=1,ndomains
      least_loaded=MINLOC(cpu_load,1)
      cpu_load(least_loaded)=cpu_load(least_loaded)+domain_load(idomain) 
      almo_scf_env%cpu_of_domain(index0(idomain))=least_loaded-1
    ENDDO

    DEALLOCATE(cpu_load)
    DEALLOCATE(index0)
    DEALLOCATE(domain_load)

    CALL timestop(handle)
  
  END SUBROUTINE distribute_domains 

! *****************************************************************************
!> \brief Tests construction and release of domain submatrices 
!> \param matrix_no ...
!> \param matrix_nn ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param error ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE construct_test(matrix_no,matrix_nn,dpattern,map,node_of_domain,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_no, matrix_nn, dpattern
    TYPE(domain_map_type), INTENT(IN)        :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: GroupID, handle, ndomains
    TYPE(cp_dbcsr_type)                      :: copy1
    TYPE(domain_submatrix_type), &
      ALLOCATABLE, DIMENSION(:)              :: subm_nn, subm_no

    CALL timeset(routineN,handle)

    ndomains = cp_dbcsr_nblkcols_total(dpattern)
    GroupID = dbcsr_mp_group(dbcsr_distribution_mp(&
       cp_dbcsr_distribution(dpattern)))
    
    ALLOCATE(subm_no(ndomains),subm_nn(ndomains))
    CALL init_submatrices(subm_no,error)
    CALL init_submatrices(subm_nn,error)

    !CALL cp_dbcsr_print(matrix_nn,error=error)
    !CALL construct_submatrices(matrix_nn,subm_nn,dpattern,map,select_row_col,error)
    !CALL print_submatrices(subm_nn,GroupID,error)
    
    !CALL cp_dbcsr_print(matrix_no,error=error)
    CALL construct_submatrices(matrix_no,subm_no,dpattern,map,node_of_domain,select_row,error)
    CALL print_submatrices(subm_no,GroupID,error)
    
    CALL cp_dbcsr_init(copy1,error=error)
    CALL cp_dbcsr_create(copy1,template=matrix_no,error=error)
    CALL cp_dbcsr_copy(copy1,matrix_no,error=error)
    CALL cp_dbcsr_print(copy1,error=error)
    CALL construct_dbcsr_from_submatrices(copy1,subm_no,dpattern,error)
    CALL cp_dbcsr_print(copy1,error=error)
    CALL cp_dbcsr_release(copy1,error=error)
    
    CALL release_submatrices(subm_no,error)
    CALL release_submatrices(subm_nn,error)
    DEALLOCATE(subm_no,subm_nn)

    CALL timestop(handle)
  
  END SUBROUTINE construct_test 

END MODULE almo_scf_methods

