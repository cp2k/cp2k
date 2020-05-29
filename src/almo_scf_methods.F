!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Subroutines for ALMO SCF
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE almo_scf_methods
   USE almo_scf_types,                  ONLY: almo_scf_env_type,&
                                              almo_scf_history_type
   USE bibliography,                    ONLY: Kolafa2004,&
                                              Kuhne2007,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, dbcsr_frobenius_norm, &
        dbcsr_get_block_p, dbcsr_get_diag, dbcsr_get_info, dbcsr_get_stored_coordinates, &
        dbcsr_init_random, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_print, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_scale_by_vector, dbcsr_set, dbcsr_set_diag, dbcsr_transposed, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_work_create
   USE domain_submatrix_methods,        ONLY: &
        add_submatrices, construct_dbcsr_from_submatrices, construct_submatrices, &
        copy_submatrices, copy_submatrix_data, init_submatrices, multiply_submatrices, &
        print_submatrices, release_submatrices
   USE domain_submatrix_types,          ONLY: domain_map_type,&
                                              domain_submatrix_type,&
                                              select_row,&
                                              select_row_col
   USE fermi_utils,                     ONLY: FermiFixed
!! for smearing
   USE input_constants,                 ONLY: almo_domain_layout_molecular,&
                                              almo_mat_distr_atomic,&
                                              almo_scf_diag,&
                                              spd_inversion_dense_cholesky,&
                                              spd_inversion_ls_hotelling,&
                                              spd_inversion_ls_taylor
   USE iterate_matrix,                  ONLY: invert_Hotelling,&
                                              invert_Taylor,&
                                              matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: binomial
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'almo_scf_methods'

   PUBLIC almo_scf_ks_to_ks_blk, almo_scf_p_blk_to_t_blk, &
      almo_scf_t_to_proj, almo_scf_ks_blk_to_tv_blk, &
      almo_scf_ks_xx_to_tv_xx, almo_scf_t_rescaling, &
      apply_projector, get_overlap, &
      generator_to_unitary, &
      orthogonalize_mos, &
      pseudo_invert_diagonal_blk, construct_test, &
      construct_domain_preconditioner, &
      apply_domain_operators, &
      construct_domain_s_inv, &
      construct_domain_s_sqrt, &
      distribute_domains, &
      almo_scf_ks_to_ks_xx, &
      construct_domain_r_down, &
      xalmo_initial_guess, &
      fill_matrix_with_ones

CONTAINS

! **************************************************************************************************
!> \brief Fill all matrix blocks with 1.0_dp
!> \param matrix ...
!> \par History
!>       2019.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE fill_matrix_with_ones(matrix)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      INTEGER                                            :: col, hold, iblock_col, iblock_row, &
                                                            mynode, nblkcols_tot, nblkrows_tot, row
      LOGICAL                                            :: tr
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_new_block
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL dbcsr_get_info(matrix, distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)
      CALL dbcsr_work_create(matrix, work_mutable=.TRUE.)
      nblkrows_tot = dbcsr_nblkrows_total(matrix)
      nblkcols_tot = dbcsr_nblkcols_total(matrix)
      DO row = 1, nblkrows_tot
         DO col = 1, nblkcols_tot
            tr = .FALSE.
            iblock_row = row
            iblock_col = col
            CALL dbcsr_get_stored_coordinates(matrix, &
                                              iblock_row, iblock_col, hold)
            IF (hold .EQ. mynode) THEN
               NULLIFY (p_new_block)
               CALL dbcsr_reserve_block2d(matrix, &
                                          iblock_row, iblock_col, p_new_block)
               CPASSERT(ASSOCIATED(p_new_block))
               p_new_block(:, :) = 1.0_dp
            ENDIF
         ENDDO
      ENDDO
      CALL dbcsr_finalize(matrix)

   END SUBROUTINE fill_matrix_with_ones

! **************************************************************************************************
!> \brief builds projected KS matrices for the overlapping domains
!>        also computes the DIIS error vector as a by-product
!> \param almo_scf_env ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_ks_to_ks_xx(almo_scf_env)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_to_ks_xx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ndomains
      REAL(KIND=dp)                                      :: eps_multiply
      TYPE(dbcsr_type) :: matrix_tmp1, matrix_tmp2, matrix_tmp3, matrix_tmp4, matrix_tmp5, &
         matrix_tmp6, matrix_tmp7, matrix_tmp8, matrix_tmp9
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_tmp1, subm_tmp2, subm_tmp3

      CALL timeset(routineN, handle)

      eps_multiply = almo_scf_env%eps_filter

      DO ispin = 1, almo_scf_env%nspins

         ndomains = dbcsr_nblkcols_total(almo_scf_env%quench_t(ispin))

         ! 0. Create KS_xx
         CALL construct_submatrices( &
            almo_scf_env%matrix_ks(ispin), &
            almo_scf_env%domain_ks_xx(:, ispin), &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row_col)

         !!!!! RZK-warning MAKE SURE THAT YOU NEED BLOCKS OUTSIDE QUENCH_T
         !!!!! FOR ALL NO-MATRICES NOT COMPUTING THEM CAN SAVE LOTS OF TIME

         ! 1. TMP1=KS.T
         !    Cost: NOn
         !matrix_tmp1 = create NxO, full
         CALL dbcsr_create(matrix_tmp1, &
                           template=almo_scf_env%matrix_t(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_ks(ispin), &
                             almo_scf_env%matrix_t(ispin), &
                             0.0_dp, matrix_tmp1, &
                             filter_eps=eps_multiply)

         ! 2. TMP2=TMP1.SigInv=KS.T.SigInv
         !    Cost: NOO
         !matrix_tmp2 = create NxO, full
         CALL dbcsr_create(matrix_tmp2, &
                           template=almo_scf_env%matrix_t(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                             almo_scf_env%matrix_sigma_inv(ispin), &
                             0.0_dp, matrix_tmp2, &
                             filter_eps=eps_multiply)

         ! 3. TMP1=S.T
         !    Cost: NOn
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1), &
                             almo_scf_env%matrix_t(ispin), &
                             0.0_dp, matrix_tmp1, &
                             filter_eps=eps_multiply)

         ! 4. TMP4=TMP2.tr(TMP1)=KS.T.SigInv.tr(T).S
         !    Cost: NNO
         !matrix_tmp4 = create NxN
         CALL dbcsr_create(matrix_tmp4, &
                           template=almo_scf_env%matrix_s(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp2, &
                             matrix_tmp1, &
                             0.0_dp, matrix_tmp4, &
                             filter_eps=eps_multiply)

         ! 5. KS_xx=KS_xx-TMP4_xx-tr(TMP4_xx)
         ALLOCATE (subm_tmp1(ndomains))
         CALL init_submatrices(subm_tmp1)
         CALL construct_submatrices( &
            matrix_tmp4, &
            subm_tmp1, &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row_col)
         CALL add_submatrices(1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                              -1.0_dp, subm_tmp1, 'N')
         CALL add_submatrices(1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                              -1.0_dp, subm_tmp1, 'T')

         ! 6. TMP3=tr(TMP4).T=S.T.SigInv.tr(T).KS.T
         !    Cost: NOn
         !matrix_tmp3 = create NxO, full
         CALL dbcsr_create(matrix_tmp3, &
                           template=almo_scf_env%matrix_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             matrix_tmp4, &
                             almo_scf_env%matrix_t(ispin), &
                             0.0_dp, matrix_tmp3, &
                             filter_eps=eps_multiply)
         CALL dbcsr_release(matrix_tmp4)

         ! 8. TMP6=TMP3.SigInv=S.T.SigInv.tr(T).KS.T.SigInv
         !    Cost: NOO
         !matrix_tmp6 = create NxO, full
         CALL dbcsr_create(matrix_tmp6, &
                           template=almo_scf_env%matrix_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             matrix_tmp3, &
                             almo_scf_env%matrix_sigma_inv(ispin), &
                             0.0_dp, matrix_tmp6, &
                             filter_eps=eps_multiply)

         ! 8A. Use intermediate matrices to evaluate the gradient/error
         !     Err=(TMP2-TMP6)_q=(KS.T.SigInv-S.T.SigInv.tr(T).KS.T.SigInv)_q
         ! error vector in AO-MO basis
         CALL dbcsr_copy(almo_scf_env%matrix_err_xx(ispin), &
                         almo_scf_env%quench_t(ispin))
         CALL dbcsr_copy(almo_scf_env%matrix_err_xx(ispin), &
                         matrix_tmp2, keep_sparsity=.TRUE.)
         CALL dbcsr_create(matrix_tmp4, &
                           template=almo_scf_env%matrix_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(matrix_tmp4, &
                         almo_scf_env%quench_t(ispin))
         CALL dbcsr_copy(matrix_tmp4, &
                         matrix_tmp6, keep_sparsity=.TRUE.)
         CALL dbcsr_add(almo_scf_env%matrix_err_xx(ispin), &
                        matrix_tmp4, 1.0_dp, -1.0_dp)
         CALL dbcsr_release(matrix_tmp4)
         !
         ! error vector in AO-AO basis
         ! RZK-warning tmp4 can be created using the sparsity pattern,
         ! then retain_sparsity can be used to perform the multiply
         ! this will save some time
         CALL dbcsr_copy(matrix_tmp3, &
                         matrix_tmp2)
         CALL dbcsr_add(matrix_tmp3, &
                        matrix_tmp6, 1.0_dp, -1.0_dp)
         CALL dbcsr_create(matrix_tmp4, &
                           template=almo_scf_env%matrix_s(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "T", 1.0_dp, &
                             matrix_tmp3, &
                             almo_scf_env%matrix_t(ispin), &
                             0.0_dp, matrix_tmp4, &
                             filter_eps=eps_multiply)
         CALL construct_submatrices( &
            matrix_tmp4, &
            almo_scf_env%domain_err(:, ispin), &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row_col)
         CALL dbcsr_release(matrix_tmp4)
         ! domain_err submatrices are in down-up representation
         ! bring them into the orthogonalized basis
         ALLOCATE (subm_tmp2(ndomains))
         CALL init_submatrices(subm_tmp2)
         CALL multiply_submatrices('N', 'N', 1.0_dp, &
                                   almo_scf_env%domain_err(:, ispin), &
                                   almo_scf_env%domain_s_sqrt(:, ispin), 0.0_dp, subm_tmp2)
         CALL multiply_submatrices('N', 'N', 1.0_dp, &
                                   almo_scf_env%domain_s_sqrt_inv(:, ispin), &
                                   subm_tmp2, 0.0_dp, almo_scf_env%domain_err(:, ispin))

         ! 9. TMP5=TMP6.tr(TMP1)=S.T.SigInv.tr(T).KS.T.SigInv.tr(T).S
         !    Cost: NNO
         !    matrix_tmp5 = create NxN, full
         ! RZK-warning tmp5 can be created using the sparsity pattern,
         ! then retain_sparsity can be used to perform the multiply
         ! this will save some time
         CALL dbcsr_create(matrix_tmp5, &
                           template=almo_scf_env%matrix_s(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "T", 1.0_dp, &
                             matrix_tmp6, &
                             matrix_tmp1, &
                             0.0_dp, matrix_tmp5, &
                             filter_eps=eps_multiply)

         ! 10. KS_xx=KS_xx+TMP5_xx
         CALL construct_submatrices( &
            matrix_tmp5, &
            subm_tmp1, &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row_col)
         CALL dbcsr_release(matrix_tmp5)
         CALL add_submatrices(1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                              1.0_dp, subm_tmp1, 'N')

         ! 11. KS_xx=KS_xx + [S.T]_xx.[SigInv.tr(T).KS.(1-T.SigInv.tr(T).S)]_xx + transposed
         ALLOCATE (subm_tmp3(ndomains))
         CALL init_submatrices(subm_tmp3)
         CALL construct_submatrices( &
            matrix_tmp2, &
            subm_tmp2, &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row)
         CALL construct_submatrices( &
            matrix_tmp6, &
            subm_tmp3, &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row)
         CALL dbcsr_release(matrix_tmp6)
         CALL add_submatrices(1.0_dp, subm_tmp2, &
                              -1.0_dp, subm_tmp3, 'N')
         CALL construct_submatrices( &
            matrix_tmp1, &
            subm_tmp3, &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row)
         CALL multiply_submatrices('N', 'T', 1.0_dp, subm_tmp2, &
                                   subm_tmp3, 0.0_dp, subm_tmp1)
         CALL add_submatrices(1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                              1.0_dp, subm_tmp1, 'N')
         CALL add_submatrices(1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                              1.0_dp, subm_tmp1, 'T')

         ! 12. TMP7=tr(T).KS.T.SigInv
         CALL dbcsr_create(matrix_tmp7, &
                           template=almo_scf_env%matrix_sigma_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             almo_scf_env%matrix_t(ispin), &
                             matrix_tmp2, &
                             0.0_dp, matrix_tmp7, &
                             filter_eps=eps_multiply)

         ! 13. TMP8=[SigInv.tr(T).KS.T.SigInv]_xx
         CALL dbcsr_create(matrix_tmp8, &
                           template=almo_scf_env%matrix_sigma_blk(ispin), &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_copy(matrix_tmp8, almo_scf_env%matrix_sigma_blk(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_sigma_inv(ispin), &
                             matrix_tmp7, &
                             0.0_dp, matrix_tmp8, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)
         CALL dbcsr_release(matrix_tmp7)

         ! 13. TMP9=[S.T]_xx
         CALL dbcsr_create(matrix_tmp9, &
                           template=almo_scf_env%matrix_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(matrix_tmp9, almo_scf_env%quench_t(ispin))
         CALL dbcsr_copy(matrix_tmp9, matrix_tmp1, keep_sparsity=.TRUE.)

         ! 14. TMP3=TMP9.TMP8=[S.T]_xx.[SigInv.tr(T).KS.T.SigInv]_xx
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             matrix_tmp9, &
                             matrix_tmp8, &
                             0.0_dp, matrix_tmp3, &
                             filter_eps=eps_multiply)
         CALL dbcsr_release(matrix_tmp8)
         CALL dbcsr_release(matrix_tmp9)

         ! 15. KS_xx=KS_xx+[S.T]_xx.[SigInv.tr(T).KS.T.SigInv]_xx.[tr(T).S]_xx
         CALL construct_submatrices( &
            matrix_tmp3, &
            subm_tmp2, &
            almo_scf_env%quench_t(ispin), &
            almo_scf_env%domain_map(ispin), &
            almo_scf_env%cpu_of_domain, &
            select_row)
         CALL multiply_submatrices('N', 'T', 1.0_dp, subm_tmp2, &
                                   subm_tmp3, 0.0_dp, subm_tmp1)
         CALL add_submatrices(1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                              1.0_dp, subm_tmp1, 'N')

         !!!!!!! use intermediate matrices to get the error vector !!!!!!!
         !!!!!!! make sure s_blk_sqrt and its inverse exist (i.e. we use diag algorithm)
         !CPPrecondition(almo_scf_env%almo_update_algorithm.eq.almo_scf_diag,cp_failure_level,routineP,failure)
         !! tmp_err = (1-S.T_blk.SigInv.tr(T_blk)).F.T_blk.SigInv
         !CALL dbcsr_init(matrix_tmp_err)
         !CALL dbcsr_create(matrix_tmp_err,&
         !        template=almo_scf_env%matrix_t(ispin))
         !CALL dbcsr_copy(matrix_tmp_err,&
         !        matrix_tmp2)
         !CALL dbcsr_add(matrix_tmp_err,matrix_tmp3,&
         !        1.0_dp,-1.0_dp)
         !! err_blk = tmp_err.tr(T_blk)
         !CALL dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
         !        almo_scf_env%matrix_s_blk_sqrt(1))
         !CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp_err,&
         !        almo_scf_env%matrix_t(ispin),&
         !        0.0_dp, almo_scf_env%matrix_err_blk(ispin),&
         !        retain_sparsity=.TRUE.,&
         !        filter_eps=eps_multiply)
         !CALL dbcsr_release(matrix_tmp_err)
         !! bring to the orthogonal basis
         !! err_blk = (S_blk^-1/2).err_blk.(S_blk^1/2)
         !CALL dbcsr_init(matrix_tmp_err)
         !CALL dbcsr_create(matrix_tmp_err,&
         !        template=almo_scf_env%matrix_err_blk(ispin))
         !CALL dbcsr_multiply("N", "N", 1.0_dp,&
         !        almo_scf_env%matrix_err_blk(ispin),&
         !        almo_scf_env%matrix_s_blk_sqrt(1),&
         !        0.0_dp, matrix_tmp_err,&
         !        filter_eps=eps_multiply)
         !CALL dbcsr_multiply("N", "N", 1.0_dp,&
         !        almo_scf_env%matrix_s_blk_sqrt_inv(1),&
         !        matrix_tmp_err,&
         !        0.0_dp, almo_scf_env%matrix_err_blk(ispin),&
         !        filter_eps=eps_multiply)
         !! subtract transpose
         !CALL dbcsr_transposed(matrix_tmp_err,&
         !        almo_scf_env%matrix_err_blk(ispin))
         !CALL dbcsr_add(almo_scf_env%matrix_err_blk(ispin),&
         !        matrix_tmp_err,&
         !        1.0_dp,-1.0_dp)
         !CALL dbcsr_release(matrix_tmp_err)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CALL release_submatrices(subm_tmp3)
         CALL release_submatrices(subm_tmp2)
         CALL release_submatrices(subm_tmp1)
         DEALLOCATE (subm_tmp3)
         DEALLOCATE (subm_tmp2)
         DEALLOCATE (subm_tmp1)
         CALL dbcsr_release(matrix_tmp3)
         CALL dbcsr_release(matrix_tmp2)
         CALL dbcsr_release(matrix_tmp1)

      ENDDO ! spins

      CALL timestop(handle)

   END SUBROUTINE almo_scf_ks_to_ks_xx

! **************************************************************************************************
!> \brief computes the projected KS from the total KS matrix
!>        also computes the DIIS error vector as a by-product
!> \param almo_scf_env ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_ks_to_ks_blk(almo_scf_env)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_to_ks_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: eps_multiply
      TYPE(dbcsr_type) :: matrix_tmp1, matrix_tmp2, matrix_tmp3, matrix_tmp4, matrix_tmp5, &
         matrix_tmp6, matrix_tmp7, matrix_tmp8, matrix_tmp9, matrix_tmp_err

      CALL timeset(routineN, handle)

      eps_multiply = almo_scf_env%eps_filter

      DO ispin = 1, almo_scf_env%nspins

         ! 1. TMP1=KS.T_blk
         !    Cost: NOn
         !matrix_tmp1 = create NxO, full
         CALL dbcsr_create(matrix_tmp1, &
                           template=almo_scf_env%matrix_t(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_ks(ispin), &
                             almo_scf_env%matrix_t_blk(ispin), &
                             0.0_dp, matrix_tmp1, &
                             filter_eps=eps_multiply)
         ! 2. TMP2=TMP1.SigInv=KS.T_blk.SigInv
         !    Cost: NOO
         !matrix_tmp2 = create NxO, full
         CALL dbcsr_create(matrix_tmp2, &
                           template=almo_scf_env%matrix_t(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                             almo_scf_env%matrix_sigma_inv(ispin), &
                             0.0_dp, matrix_tmp2, &
                             filter_eps=eps_multiply)

         !!!!!! use intermediate matrices to get the error vector !!!!!!!
         !CALL dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
         !        almo_scf_env%matrix_t_blk(ispin))
         !CALL dbcsr_copy(almo_scf_env%matrix_err_blk(ispin),&
         !        matrix_tmp2,&
         !        keep_sparsity=.TRUE.)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! 3. TMP1=S.T_blk
         !    Cost: NOn
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s(1), &
                             almo_scf_env%matrix_t_blk(ispin), &
                             0.0_dp, matrix_tmp1, &
                             filter_eps=eps_multiply)

         ! 4. TMP4_blk=TMP2.tr(TMP1)=KS.T_blk.SigInv.tr(T_blk).S
         !    Cost: NnO
         !matrix_tmp4 = create NxN, blk
         CALL dbcsr_create(matrix_tmp4, &
                           template=almo_scf_env%matrix_s_blk(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(matrix_tmp4, almo_scf_env%matrix_s_blk(1))
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp2, &
                             matrix_tmp1, &
                             0.0_dp, matrix_tmp4, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)

         ! 5. KS_blk=KS_blk-TMP4_blk
         CALL dbcsr_copy(almo_scf_env%matrix_ks_blk(ispin), &
                         almo_scf_env%matrix_ks(ispin), keep_sparsity=.TRUE.)
         CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), &
                        matrix_tmp4, &
                        1.0_dp, -1.0_dp)

         ! 6. TMP5_blk=tr(TMP4_blk)
         !    KS_blk=KS_blk-tr(TMP4_blk)
         !matrix_tmp5 = create NxN, blk
         CALL dbcsr_create(matrix_tmp5, &
                           template=almo_scf_env%matrix_s_blk(1), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_transposed(matrix_tmp5, matrix_tmp4)
         CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), matrix_tmp5, &
                        1.0_dp, -1.0_dp)

         ! 7. TMP3=tr(T_blk).TMP2=tr(T_blk).KS.T_blk.SigInv
         !    Cost: OOn
         !matrix_tmp3 = create OxO, full
         CALL dbcsr_create(matrix_tmp3, &
                           template=almo_scf_env%matrix_sigma_inv(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             almo_scf_env%matrix_t_blk(ispin), &
                             matrix_tmp2, &
                             0.0_dp, matrix_tmp3, &
                             filter_eps=eps_multiply)

         ! 8. TMP6=SigInv.TMP3=SigInv.tr(T_blk).KS.T_blk.SigInv
         !    Cost: OOO
         !matrix_tmp6 = create OxO, full
         CALL dbcsr_create(matrix_tmp6, &
                           template=almo_scf_env%matrix_sigma_inv(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_sigma_inv(ispin), &
                             matrix_tmp3, &
                             0.0_dp, matrix_tmp6, &
                             filter_eps=eps_multiply)

         ! 9. TMP3=TMP1.TMP6=S.T_blk.SigInv.tr(T_blk).KS.T_blk.SigInv
         !    Cost: NOO
         !matrix_tmp3 = re-create NxO, full
         CALL dbcsr_release(matrix_tmp3)
         CALL dbcsr_create(matrix_tmp3, &
                           template=almo_scf_env%matrix_t(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, &
                             matrix_tmp6, &
                             0.0_dp, matrix_tmp3, &
                             filter_eps=eps_multiply)

         !!!!!! use intermediate matrices to get the error vector !!!!!!!
         !CALL dbcsr_init(matrix_tmp_err)
         !CALL dbcsr_create(matrix_tmp_err,&
         !        template=almo_scf_env%matrix_t_blk(ispin))
         !CALL dbcsr_copy(matrix_tmp_err,&
         !        almo_scf_env%matrix_t_blk(ispin))
         !CALL dbcsr_copy(matrix_tmp_err,matrix_tmp3,&
         !        keep_sparsity=.TRUE.)
         !CALL dbcsr_add(almo_scf_env%matrix_err_blk(ispin),matrix_tmp_err,&
         !        1.0_dp,-1.0_dp)
         !CALL dbcsr_release(matrix_tmp_err)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!! use intermediate matrices to get the error vector !!!!!!!
         !!!!!! make sure s_blk_sqrt and its inverse exist (i.e. we use diag algorithm)
         CPASSERT(almo_scf_env%almo_update_algorithm .EQ. almo_scf_diag)
         ! tmp_err = (1-S.T_blk.SigInv.tr(T_blk)).F.T_blk.SigInv
         CALL dbcsr_create(matrix_tmp_err, &
                           template=almo_scf_env%matrix_t_blk(ispin))
         CALL dbcsr_copy(matrix_tmp_err, &
                         matrix_tmp2)
         CALL dbcsr_add(matrix_tmp_err, matrix_tmp3, &
                        1.0_dp, -1.0_dp)
         ! err_blk = tmp_err.tr(T_blk)
         CALL dbcsr_copy(almo_scf_env%matrix_err_blk(ispin), &
                         almo_scf_env%matrix_s_blk_sqrt(1))
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp_err, &
                             almo_scf_env%matrix_t_blk(ispin), &
                             0.0_dp, almo_scf_env%matrix_err_blk(ispin), &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)
         CALL dbcsr_release(matrix_tmp_err)
         ! bring to the orthogonal basis
         ! err_blk = (S_blk^-1/2).err_blk.(S_blk^1/2)
         CALL dbcsr_create(matrix_tmp_err, &
                           template=almo_scf_env%matrix_err_blk(ispin))
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_err_blk(ispin), &
                             almo_scf_env%matrix_s_blk_sqrt(1), &
                             0.0_dp, matrix_tmp_err, &
                             filter_eps=eps_multiply)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             almo_scf_env%matrix_s_blk_sqrt_inv(1), &
                             matrix_tmp_err, &
                             0.0_dp, almo_scf_env%matrix_err_blk(ispin), &
                             filter_eps=eps_multiply)

         ! subtract transpose
         CALL dbcsr_transposed(matrix_tmp_err, &
                               almo_scf_env%matrix_err_blk(ispin))
         CALL dbcsr_add(almo_scf_env%matrix_err_blk(ispin), &
                        matrix_tmp_err, &
                        1.0_dp, -1.0_dp)
         CALL dbcsr_release(matrix_tmp_err)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! later we will need only the blk version of TMP6
         ! create it here and release TMP6
         !matrix_tmp9 = create OxO, blk
         !matrix_tmp9 = copy data from matrix_tmp6, retain sparsity
         !matrix_tmp6 = release
         CALL dbcsr_create(matrix_tmp9, &
                           template=almo_scf_env%matrix_sigma_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_copy(matrix_tmp9, almo_scf_env%matrix_sigma_blk(ispin))
         CALL dbcsr_copy(matrix_tmp9, matrix_tmp6, keep_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp6)

         !10. KS_blk=KS_blk+TMP3.tr(TMP1)=
         !          =KS_blk+S.T_blk.SigInv.tr.(T_blk).KS.T_blk.SigInv.tr(T_blk).S
         !    Cost: NnO
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp3, &
                             matrix_tmp1, &
                             1.0_dp, almo_scf_env%matrix_ks_blk(ispin), &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)

         ! 11. TMP4_blk=TMP7_blk.tr(TMP8_blk)
         !    Cost: Nnn
         !matrix_tmp7 = create NxO, blk
         !matrix_tmp7 = copy data from matrix_tmp3, retain sparsity
         !matrix_tmp3 = release
         !matrix_tmp8 = create NxO, blk
         !matrix_tmp8 = copy data from matrix_tmp1, retain sparsity
         !matrix_tmp1 = release
         CALL dbcsr_create(matrix_tmp7, &
                           template=almo_scf_env%matrix_t_blk(ispin))
         ! transfer only the ALMO blocks from tmp3 into tmp7:
         ! first, copy t_blk into tmp7 to transfer the blk structure,
         ! then copy tmp3 into tmp7 with retain_sparsity
         CALL dbcsr_copy(matrix_tmp7, almo_scf_env%matrix_t_blk(ispin))
         CALL dbcsr_copy(matrix_tmp7, matrix_tmp3, keep_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp3)
         ! do the same for tmp1->tmp8
         CALL dbcsr_create(matrix_tmp8, &
                           template=almo_scf_env%matrix_t_blk(ispin))
         CALL dbcsr_copy(matrix_tmp8, almo_scf_env%matrix_t_blk(ispin))
         CALL dbcsr_copy(matrix_tmp8, matrix_tmp1, keep_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp1)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7, &
                             matrix_tmp8, &
                             0.0_dp, matrix_tmp4, &
                             filter_eps=eps_multiply, &
                             retain_sparsity=.TRUE.)

         ! 12. KS_blk=KS_blk-TMP4_blk
         CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), matrix_tmp4, &
                        1.0_dp, -1.0_dp)

         ! 13. TMP5_blk=tr(TMP5_blk)
         !     KS_blk=KS_blk-tr(TMP4_blk)
         CALL dbcsr_transposed(matrix_tmp5, matrix_tmp4)
         CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), matrix_tmp5, &
                        1.0_dp, -1.0_dp)

         ! 14. TMP4_blk=TMP7_blk.tr(TMP8_blk)
         !     Cost: Nnn
         CALL dbcsr_copy(matrix_tmp7, matrix_tmp2, keep_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp2)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7, &
                             matrix_tmp8, &
                             0.0_dp, matrix_tmp4, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)
         ! 15. KS_blk=KS_blk+TMP4_blk
         CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), matrix_tmp4, &
                        1.0_dp, 1.0_dp)

         ! 16. KS_blk=KS_blk+tr(TMP4_blk)
         CALL dbcsr_transposed(matrix_tmp5, matrix_tmp4)
         CALL dbcsr_release(matrix_tmp4)
         CALL dbcsr_add(almo_scf_env%matrix_ks_blk(ispin), matrix_tmp5, &
                        1.0_dp, 1.0_dp)
         CALL dbcsr_release(matrix_tmp5)

         ! 17. TMP10_blk=TMP8_blk.TMP9_blk
         !    Cost: Noo
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp8, &
                             matrix_tmp9, &
                             0.0_dp, matrix_tmp7, &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)
         CALL dbcsr_release(matrix_tmp9)

         ! 18. KS_blk=TMP7_blk.tr(TMP8_blk)
         !    Cost: Nno
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7, &
                             matrix_tmp8, &
                             1.0_dp, almo_scf_env%matrix_ks_blk(ispin), &
                             retain_sparsity=.TRUE., &
                             filter_eps=eps_multiply)
         CALL dbcsr_release(matrix_tmp7)
         CALL dbcsr_release(matrix_tmp8)

      ENDDO ! spins

      CALL timestop(handle)

   END SUBROUTINE almo_scf_ks_to_ks_blk

! **************************************************************************************************
!> \brief ALMOs by diagonalizing the KS domain submatrices
!>        computes both the occupied and virtual orbitals
!> \param almo_scf_env ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_ks_xx_to_tv_xx(almo_scf_env)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_xx_to_tv_xx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iblock_size, idomain, info, &
                                                            ispin, lwork, ndomains
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, work
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: data_copy
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_ks_xx_orthog, subm_t, subm_tmp

      CALL timeset(routineN, handle)

      IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular .AND. &
          almo_scf_env%mat_distr_aos == almo_mat_distr_atomic) THEN
         CPABORT("a domain must be located entirely on a CPU")
      ENDIF

      ndomains = almo_scf_env%ndomains
      ALLOCATE (subm_tmp(ndomains))
      ALLOCATE (subm_ks_xx_orthog(ndomains))
      ALLOCATE (subm_t(ndomains))

      DO ispin = 1, almo_scf_env%nspins

         CALL init_submatrices(subm_tmp)
         CALL init_submatrices(subm_ks_xx_orthog)

         ! TRY: project out T0-occupied space for each domain
         ! F=(1-R_du).F.(1-tr(R_du))
         !CALL copy_submatrices(almo_scf_env%domain_ks_xx(:,ispin),&
         !        subm_ks_xx_orthog,copy_data=.TRUE.)
         !CALL multiply_submatrices('N','N',1.0_dp,&
         !        almo_scf_env%domain_r_down_up(:,ispin),&
         !        almo_scf_env%domain_ks_xx(:,ispin),0.0_dp,subm_tmp)
         !CALL add_submatrices(1.0_dp,subm_ks_xx_orthog,-1.0_dp,subm_tmp,'N')
         !CALL add_submatrices(1.0_dp,subm_ks_xx_orthog,-1.0_dp,subm_tmp,'T')
         !!CALL multiply_submatrices('N','T',1.0_dp,subm_tmp,&
         !!        almo_scf_env%domain_r_down_up(:,ispin),&
         !!        1.0_dp,subm_ks_xx_orthog)

         ! convert blocks to the orthogonal basis set
         ! TRY: replace one multiply
         !CALL multiply_submatrices('N','N',1.0_dp,subm_ks_xx_orthog,&
         !        almo_scf_env%domain_s_sqrt_inv(:,ispin),0.0_dp,subm_tmp)
         CALL multiply_submatrices('N', 'N', 1.0_dp, almo_scf_env%domain_ks_xx(:, ispin), &
                                   almo_scf_env%domain_s_sqrt_inv(:, ispin), 0.0_dp, subm_tmp)
         CALL multiply_submatrices('N', 'N', 1.0_dp, almo_scf_env%domain_s_sqrt_inv(:, ispin), &
                                   subm_tmp, 0.0_dp, subm_ks_xx_orthog)
         CALL release_submatrices(subm_tmp)

         ! create temporary matrices for occupied and virtual orbitals
         ! represented in the orthogonalized basis set
         CALL init_submatrices(subm_t)

         ! loop over domains - perform diagonalization
         DO idomain = 1, ndomains

            ! check if the submatrix exists
            IF (subm_ks_xx_orthog(idomain)%domain .GT. 0) THEN

               iblock_size = subm_ks_xx_orthog(idomain)%nrows

               ! Prepare data
               ALLOCATE (eigenvalues(iblock_size))
               ALLOCATE (data_copy(iblock_size, iblock_size))
               data_copy(:, :) = subm_ks_xx_orthog(idomain)%mdata(:, :)

               ! Query the optimal workspace for dsyev
               LWORK = -1
               ALLOCATE (WORK(MAX(1, LWORK)))
               CALL DSYEV('V', 'L', iblock_size, data_copy, iblock_size, eigenvalues, WORK, LWORK, INFO)
               LWORK = INT(WORK(1))
               DEALLOCATE (WORK)

               ! Allocate the workspace and solve the eigenproblem
               ALLOCATE (WORK(MAX(1, LWORK)))
               CALL DSYEV('V', 'L', iblock_size, data_copy, iblock_size, eigenvalues, WORK, LWORK, INFO)
               IF (INFO .NE. 0) THEN
                  CPABORT("DSYEV failed")
               END IF

               ! Copy occupied eigenvectors
               IF (almo_scf_env%domain_t(idomain, ispin)%ncols .NE. &
                   almo_scf_env%nocc_of_domain(idomain, ispin)) THEN
                  CPABORT("wrong domain structure")
               END IF
               CALL copy_submatrices(almo_scf_env%domain_t(idomain, ispin), &
                                     subm_t(idomain), .FALSE.)
               CALL copy_submatrix_data(data_copy(:, 1:almo_scf_env%nocc_of_domain(idomain, ispin)), &
                                        subm_t(idomain))
               !! Copy occupied eigenvalues if smearing requested
               IF (almo_scf_env%smear) THEN
                  almo_scf_env%mo_energies(1 + SUM(almo_scf_env%nocc_of_domain(:idomain - 1, ispin)) &
                                           :SUM(almo_scf_env%nocc_of_domain(:idomain, ispin)), ispin) &
                     = eigenvalues(1:almo_scf_env%nocc_of_domain(idomain, ispin))
               END IF

               DEALLOCATE (WORK)
               DEALLOCATE (data_copy)
               DEALLOCATE (eigenvalues)

            ENDIF ! submatrix for the domain exists

         ENDDO ! loop over domains

         CALL release_submatrices(subm_ks_xx_orthog)

         ! convert orbitals to the AO basis set (from orthogonalized AOs)
         CALL multiply_submatrices('N', 'N', 1.0_dp, almo_scf_env%domain_s_sqrt_inv(:, ispin), &
                                   subm_t, 0.0_dp, almo_scf_env%domain_t(:, ispin))
         CALL release_submatrices(subm_t)

         ! convert domain orbitals to a dbcsr matrix
         CALL construct_dbcsr_from_submatrices( &
            almo_scf_env%matrix_t(ispin), &
            almo_scf_env%domain_t(:, ispin), &
            almo_scf_env%quench_t(ispin))
         CALL dbcsr_filter(almo_scf_env%matrix_t(ispin), &
                           almo_scf_env%eps_filter)

         ! TRY: add T0 component
         !!CALL dbcsr_add(almo_scf_env%matrix_t(ispin),&
         !!        almo_scf_env%matrix_t_blk(ispin),1.0_dp,1.0_dp)

      ENDDO ! spins

      DEALLOCATE (subm_tmp)
      DEALLOCATE (subm_ks_xx_orthog)
      DEALLOCATE (subm_t)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_ks_xx_to_tv_xx

! **************************************************************************************************
!> \brief computes ALMOs by diagonalizing the projected blocked KS matrix
!>        uses the diagonalization code for blocks
!>        computes both the occupied and virtual orbitals
!> \param almo_scf_env ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_ks_blk_to_tv_blk(almo_scf_env)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ks_blk_to_tv_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iblock_col, iblock_row, &
                                                            iblock_size, info, ispin, lwork, &
                                                            nocc_of_block, nvirt_of_block, orbital
      LOGICAL                                            :: block_needed
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, work
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: data_copy
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_ks_blk_orthog, &
                                                            matrix_t_blk_orthog, matrix_tmp, &
                                                            matrix_v_blk_orthog

      CALL timeset(routineN, handle)

      IF (almo_scf_env%domain_layout_aos == almo_domain_layout_molecular .AND. &
          almo_scf_env%mat_distr_aos == almo_mat_distr_atomic) THEN
         CPABORT("a domain must be located entirely on a CPU")
      ENDIF

      DO ispin = 1, almo_scf_env%nspins

         CALL dbcsr_create(matrix_tmp, template=almo_scf_env%matrix_ks_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_ks_blk_orthog, template=almo_scf_env%matrix_ks_blk(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         ! convert blocks to the orthogonal basis set
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_ks_blk(ispin), &
                             almo_scf_env%matrix_s_blk_sqrt_inv(1), 0.0_dp, matrix_tmp, &
                             filter_eps=almo_scf_env%eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_blk_sqrt_inv(1), &
                             matrix_tmp, 0.0_dp, matrix_ks_blk_orthog, &
                             filter_eps=almo_scf_env%eps_filter)

         CALL dbcsr_release(matrix_tmp)

         ! create temporary matrices for occupied and virtual orbitals
         ! represented in the orthogonalized AOs basis set
         CALL dbcsr_create(matrix_t_blk_orthog, template=almo_scf_env%matrix_t_blk(ispin))
         CALL dbcsr_create(matrix_v_blk_orthog, template=almo_scf_env%matrix_v_full_blk(ispin))
         CALL dbcsr_work_create(matrix_t_blk_orthog, work_mutable=.TRUE.)
         CALL dbcsr_work_create(matrix_v_blk_orthog, work_mutable=.TRUE.)

         CALL dbcsr_work_create(almo_scf_env%matrix_eoo(ispin), work_mutable=.TRUE.)
         CALL dbcsr_work_create(almo_scf_env%matrix_evv_full(ispin), work_mutable=.TRUE.)

         CALL dbcsr_iterator_start(iter, matrix_ks_blk_orthog)

         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, row_size=iblock_size)

            IF (iblock_row .NE. iblock_col) THEN
               CPABORT("off-diagonal block found")
            ENDIF

            block_needed = .TRUE.
            IF (almo_scf_env%nocc_of_domain(iblock_col, ispin) .EQ. 0 .AND. &
                almo_scf_env%nvirt_of_domain(iblock_col, ispin) .EQ. 0) THEN
               block_needed = .FALSE.
            ENDIF

            IF (block_needed) THEN

               ! Prepare data
               ALLOCATE (eigenvalues(iblock_size))
               ALLOCATE (data_copy(iblock_size, iblock_size))
               data_copy(:, :) = data_p(:, :)

               ! Query the optimal workspace for dsyev
               LWORK = -1
               ALLOCATE (WORK(MAX(1, LWORK)))
               CALL DSYEV('V', 'L', iblock_size, data_copy, iblock_size, eigenvalues, WORK, LWORK, INFO)
               LWORK = INT(WORK(1))
               DEALLOCATE (WORK)

               ! Allocate the workspace and solve the eigenproblem
               ALLOCATE (WORK(MAX(1, LWORK)))
               CALL DSYEV('V', 'L', iblock_size, data_copy, iblock_size, eigenvalues, WORK, LWORK, INFO)
               IF (INFO .NE. 0) THEN
                  CPABORT("DSYEV failed")
               END IF

               !!! RZK-warning                                               !!!
               !!! IT IS EXTREMELY IMPORTANT THAT THE DIAGONAL BLOCKS OF THE !!!
               !!! FOLLOWING MATRICES ARE LOCATED ON THE SAME NODES WITH     !!!
               !!! THE CORRESPONDING DIAGONAL BLOCKS OF THE FOCK MATRIX:     !!!
               !!! T, V, E_o, E_v

               ! copy eigenvectors into two dbcsr matrices - occupied and virtuals
               nocc_of_block = almo_scf_env%nocc_of_domain(iblock_col, ispin)
               IF (nocc_of_block .GT. 0) THEN
                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(matrix_t_blk_orthog, iblock_row, iblock_col, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:, :) = data_copy(:, 1:nocc_of_block)
                  ! copy eigenvalues into diagonal dbcsr matrix - Eoo
                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(almo_scf_env%matrix_eoo(ispin), iblock_row, iblock_col, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:, :) = 0.0_dp
                  DO orbital = 1, nocc_of_block
                     p_new_block(orbital, orbital) = eigenvalues(orbital)
                  ENDDO
                  !! Retrieve occupied MOs energies for smearing purpose, if requested
                  !! RS-WARNING: Hack to retrieve the occupied energies, since matrix_eoo seems to be ill-defined
                  !!             for multiprocessing (any idea for fix?)
                  !! RS-WARNING: This section is not suitable for parallel run !!!
                  !!             (but usually fails less than retrieving the diagonal of matrix_eoo)
                  !! RS-WARNING: This method will likely keep the energies of the initial guess if run in parallel
                  !!             (which is still a reasonable smearing in most cases...)
                  IF (almo_scf_env%smear) THEN
                     DO orbital = 1, nocc_of_block
                        almo_scf_env%mo_energies(SUM(almo_scf_env%nocc_of_domain(:iblock_row - 1, ispin)) + orbital, &
                                                 ispin) = eigenvalues(orbital)
                     END DO
                  END IF
               ENDIF

               ! now virtuals
               nvirt_of_block = almo_scf_env%nvirt_of_domain(iblock_col, ispin)
               IF (nvirt_of_block .GT. 0) THEN
                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(matrix_v_blk_orthog, iblock_row, iblock_col, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:, :) = data_copy(:, (nocc_of_block + 1):(nocc_of_block + nvirt_of_block))
                  ! virtual energies
                  NULLIFY (p_new_block)
                  CALL dbcsr_reserve_block2d(almo_scf_env%matrix_evv_full(ispin), iblock_row, iblock_col, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(:, :) = 0.0_dp
                  DO orbital = 1, nvirt_of_block
                     p_new_block(orbital, orbital) = eigenvalues(nocc_of_block + orbital)
                  ENDDO
               ENDIF

               DEALLOCATE (WORK)
               DEALLOCATE (data_copy)
               DEALLOCATE (eigenvalues)

            ENDIF

         ENDDO
         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_finalize(matrix_t_blk_orthog)
         CALL dbcsr_finalize(matrix_v_blk_orthog)
         CALL dbcsr_finalize(almo_scf_env%matrix_eoo(ispin))
         CALL dbcsr_finalize(almo_scf_env%matrix_evv_full(ispin))

         !! RS-WARNING: When matrix_eoo will be well-defined with multiprocessing,
         !!             the following should be the preferred way to retrieve the occupied energies:
         !! Retrieve occupied MOs energies for smearing purpose, if requested
         !! IF (almo_scf_env%smear) THEN
         !!    CALL dbcsr_get_diag(almo_scf_env%matrix_eoo(ispin), almo_scf_env%mo_energies(:, ispin))
         !! END IF

         CALL dbcsr_filter(matrix_t_blk_orthog, almo_scf_env%eps_filter)
         CALL dbcsr_filter(matrix_v_blk_orthog, almo_scf_env%eps_filter)

         CALL dbcsr_release(matrix_ks_blk_orthog)

         ! convert orbitals to the AO basis set (from orthogonalized AOs)
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_blk_sqrt_inv(1), &
                             matrix_t_blk_orthog, 0.0_dp, almo_scf_env%matrix_t_blk(ispin), &
                             filter_eps=almo_scf_env%eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_blk_sqrt_inv(1), &
                             matrix_v_blk_orthog, 0.0_dp, almo_scf_env%matrix_v_full_blk(ispin), &
                             filter_eps=almo_scf_env%eps_filter)

         CALL dbcsr_release(matrix_t_blk_orthog)
         CALL dbcsr_release(matrix_v_blk_orthog)

      ENDDO ! spins

      CALL timestop(handle)

   END SUBROUTINE almo_scf_ks_blk_to_tv_blk

! **************************************************************************************************
!> \brief inverts block-diagonal blocks of a dbcsr_matrix
!> \param matrix_in ...
!> \param matrix_out ...
!> \param nocc ...
!> \par History
!>       2012.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE pseudo_invert_diagonal_blk(matrix_in, matrix_out, nocc)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_out
      INTEGER, DIMENSION(:)                              :: nocc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'pseudo_invert_diagonal_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iblock_col, iblock_row, &
                                                            iblock_size, methodID
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: data_copy
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_out, template=matrix_in)
      CALL dbcsr_work_create(matrix_out, work_mutable=.TRUE.)

      CALL dbcsr_iterator_start(iter, matrix_in)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, row_size=iblock_size)

         IF (iblock_row == iblock_col) THEN

            ! Prepare data
            ALLOCATE (data_copy(iblock_size, iblock_size))
            !data_copy(:,:)=data_p(:,:)

            ! 0. Cholesky factorization
            ! 1. Diagonalization
            methodID = 1
            CALL pseudo_invert_matrix(data_p, data_copy, iblock_size, &
                                      methodID, &
                                      range1=nocc(iblock_row), range2=nocc(iblock_row), &
                                      !range1_thr,range2_thr,&
                                      shift=1.0E-5_dp)
            !!! IT IS EXTREMELY IMPORTANT THAT THE BLOCKS OF THE "OUT"  !!!
            !!! MATRIX ARE DISTRIBUTED AS THE BLOCKS OF THE "IN" MATRIX !!!

            NULLIFY (p_new_block)
            CALL dbcsr_reserve_block2d(matrix_out, iblock_row, iblock_col, p_new_block)
            CPASSERT(ASSOCIATED(p_new_block))
            p_new_block(:, :) = data_copy(:, :)

            DEALLOCATE (data_copy)

         ENDIF

      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_finalize(matrix_out)

      CALL timestop(handle)

   END SUBROUTINE pseudo_invert_diagonal_blk

! **************************************************************************************************
!> \brief computes occupied ALMOs from the superimposed atomic density blocks
!> \param almo_scf_env ...
!> \param ionic ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_p_blk_to_t_blk(almo_scf_env, ionic)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env
      LOGICAL, INTENT(IN)                                :: ionic

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_p_blk_to_t_blk', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iblock_col, iblock_row, &
                                                            iblock_size, info, ispin, lwork, &
                                                            nocc_of_block, unit_nr
      LOGICAL                                            :: block_needed
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, work
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: data_copy
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_t_blk_tmp

      CALL timeset(routineN, handle)

      ! get a useful unit_nr
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      DO ispin = 1, almo_scf_env%nspins

         IF (ionic) THEN

            ! create a temporary matrix to keep the eigenvectors
            CALL dbcsr_create(matrix_t_blk_tmp, &
                              template=almo_scf_env%matrix_t_blk(ispin))
            CALL dbcsr_work_create(matrix_t_blk_tmp, &
                                   work_mutable=.TRUE.)

            CALL dbcsr_iterator_start(iter, almo_scf_env%matrix_p_blk(ispin))
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, row_size=iblock_size)

               block_needed = .FALSE.

               IF (iblock_row == iblock_col) THEN
                  block_needed = .TRUE.
               ENDIF

               IF (.NOT. block_needed) THEN
                  CPABORT("off-diag block found")
               ENDIF

               ! Prepare data
               ALLOCATE (eigenvalues(iblock_size))
               ALLOCATE (data_copy(iblock_size, iblock_size))
               data_copy(:, :) = data_p(:, :)

               ! Query the optimal workspace for dsyev
               LWORK = -1
               ALLOCATE (WORK(MAX(1, LWORK)))
               CALL DSYEV('V', 'L', iblock_size, data_copy, iblock_size, eigenvalues, &
                          WORK, LWORK, INFO)
               LWORK = INT(WORK(1))
               DEALLOCATE (WORK)

               ! Allocate the workspace and solve the eigenproblem
               ALLOCATE (WORK(MAX(1, LWORK)))
               CALL DSYEV('V', 'L', iblock_size, data_copy, iblock_size, eigenvalues, WORK, LWORK, INFO)
               IF (INFO .NE. 0) THEN
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, *) 'BLOCK  = ', iblock_row
                     WRITE (unit_nr, *) 'INFO   =', INFO
                     WRITE (unit_nr, *) data_p(:, :)
                  END IF
                  CPABORT("DSYEV failed")
               END IF

               !!! IT IS EXTREMELY IMPORTANT THAT THE DIAGONAL BLOCKS OF THE !!!
               !!! P AND T MATRICES ARE LOCATED ON THE SAME NODES            !!!

               ! copy eigenvectors into two dbcsr matrices - occupied and virtuals
               NULLIFY (p_new_block)
               CALL dbcsr_reserve_block2d(matrix_t_blk_tmp, &
                                          iblock_row, iblock_col, p_new_block)
               nocc_of_block = SIZE(p_new_block, 2)
               CPASSERT(ASSOCIATED(p_new_block))
               CPASSERT(nocc_of_block .GT. 0)
               p_new_block(:, :) = data_copy(:, iblock_size - nocc_of_block + 1:)

               DEALLOCATE (WORK)
               DEALLOCATE (data_copy)
               DEALLOCATE (eigenvalues)

            ENDDO
            CALL dbcsr_iterator_stop(iter)

            CALL dbcsr_finalize(matrix_t_blk_tmp)
            CALL dbcsr_filter(matrix_t_blk_tmp, &
                              almo_scf_env%eps_filter)
            CALL dbcsr_copy(almo_scf_env%matrix_t_blk(ispin), &
                            matrix_t_blk_tmp)
            CALL dbcsr_release(matrix_t_blk_tmp)

         ELSE

            !! generate a random set of ALMOs
            !! matrix_t_blk should already be initiated to the proper domain structure
            CALL dbcsr_init_random(almo_scf_env%matrix_t_blk(ispin), &
                                   keep_sparsity=.TRUE.)

            CALL dbcsr_create(matrix_t_blk_tmp, &
                              template=almo_scf_env%matrix_t_blk(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            ! use current ALMOs in matrix_t_blk and project them onto the blocked dm
            ! compute T_new = R_blk S_blk T_random
            CALL dbcsr_multiply("N", "N", 1.0_dp, almo_scf_env%matrix_s_blk(1), &
                                almo_scf_env%matrix_t_blk(ispin), &
                                0.0_dp, matrix_t_blk_tmp, &
                                filter_eps=almo_scf_env%eps_filter)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                almo_scf_env%matrix_p_blk(ispin), matrix_t_blk_tmp, &
                                0.0_dp, almo_scf_env%matrix_t_blk(ispin), &
                                filter_eps=almo_scf_env%eps_filter)

            CALL dbcsr_release(matrix_t_blk_tmp)

         ENDIF

      ENDDO

      CALL timestop(handle)

   END SUBROUTINE almo_scf_p_blk_to_t_blk

! **************************************************************************************************
!> \brief Apply an occupation-rescaling trick to ALMOs for smearing.
!>        Partially occupied orbitals are considered full and rescaled by SQRT(occupation_number)
!>        (this was designed to be used with smearing only)
!> \param matrix_t ...
!> \param mo_energies ...
!> \param mu_of_domain ...
!> \param real_ne_of_domain ...
!> \param spin_kTS ...
!> \param smear_e_temp ...
!> \param ndomains ...
!> \param nocc_of_domain ...
!> \par History
!>       2018.09 created [Ruben Staub]
!> \author Ruben Staub
! **************************************************************************************************
   SUBROUTINE almo_scf_t_rescaling(matrix_t, mo_energies, mu_of_domain, real_ne_of_domain, &
                                   spin_kTS, smear_e_temp, ndomains, nocc_of_domain)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_t
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: mo_energies
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: mu_of_domain, real_ne_of_domain
      REAL(KIND=dp), INTENT(INOUT)                       :: spin_kTS
      REAL(KIND=dp), INTENT(IN)                          :: smear_e_temp
      INTEGER, INTENT(IN)                                :: ndomains
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc_of_domain

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_t_rescaling', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idomain, neigenval_used, nmo
      REAL(KIND=dp)                                      :: kTS
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupation_numbers, rescaling_factors

      CALL timeset(routineN, handle)

      !!
      !! Initialization
      !!
      nmo = SIZE(mo_energies)
      ALLOCATE (occupation_numbers(nmo))
      ALLOCATE (rescaling_factors(nmo))

      !!
      !! Set occupation numbers for smearing
      !!
      !! occupation numbers are obtained using Fermi-Dirac smearing with orbital energies stored in mo_energies
      !! nocc_of_domain is the number of partially occupied orbitals, while real_ne_of_domain is the real number of electrons
      neigenval_used = 0 !! this is used as an offset to copy sub-arrays

      !! Reset electronic entropy term
      spin_kTS = 0.0_dp

      !! Apply Fermi-Dirac smearing for each domain and store associated occupations for the whole system
      DO idomain = 1, ndomains
         CALL FermiFixed(occupation_numbers(1 + neigenval_used:nocc_of_domain(idomain) + neigenval_used), &
                         mu_of_domain(idomain), &
                         kTS, &
                         mo_energies(1 + neigenval_used:nocc_of_domain(idomain) + neigenval_used), &
                         real_ne_of_domain(idomain), &
                         smear_e_temp, &
                         1.0_dp) !! Warning, maxocc is set to 1 since we don't want to interfere with the spin_factor rescaling
         spin_kTS = spin_kTS + kTS !! Add up electronic entropy contributions
         neigenval_used = neigenval_used + nocc_of_domain(idomain) !! Update eigenvalues index offset
      END DO
      rescaling_factors(:) = SQRT(occupation_numbers) !! scale = sqrt(occupation_number)

      !!
      !! Rescaling electronic entropy contribution by spin_factor (deprecated)
      !! (currently, entropy is rescaled by spin_factor with the density matrix)
      !!
      !!IF (almo_scf_env%nspins == 1) THEN
      !!   spin_kTS = spin_kTS*2.0_dp
      !!END IF

      !!
      !! Rescaling of T (i.e. ALMOs)
      !!
      CALL dbcsr_scale_by_vector(matrix_t, rescaling_factors, side='right') !! Apply occupation-rescaling trick

      !!
      !! Debug tools (for debug purpose only)
      !!
      !! WRITE (*,*) "occupations", occupation_numbers(:) !! debug
      !! WRITE (*,*) "eigenvalues", mo_energies(:) !! debug
      !! WRITE (*,*) "kTS (spin_factor excluded) = ", spin_kTS !! debug

      !!
      !! Cleaning up before exit
      !!
      DEALLOCATE (occupation_numbers)
      DEALLOCATE (rescaling_factors)

      CALL timestop(handle)

   END SUBROUTINE almo_scf_t_rescaling

! **************************************************************************************************
!> \brief Computes the overlap matrix of MO orbitals
!> \param bra ...
!> \param ket ...
!> \param overlap ...
!> \param metric ...
!> \param retain_overlap_sparsity ...
!> \param eps_filter ...
!> \param smear ...
!> \par History
!>       2011.08 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE get_overlap(bra, ket, overlap, metric, retain_overlap_sparsity, &
                          eps_filter, smear)

      TYPE(dbcsr_type), INTENT(IN)                       :: bra, ket
      TYPE(dbcsr_type), INTENT(INOUT)                    :: overlap
      TYPE(dbcsr_type), INTENT(IN)                       :: metric
      LOGICAL, INTENT(IN), OPTIONAL                      :: retain_overlap_sparsity
      REAL(KIND=dp)                                      :: eps_filter
      LOGICAL, INTENT(IN), OPTIONAL                      :: smear

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_overlap', routineP = moduleN//':'//routineN

      INTEGER                                            :: dim0, handle
      LOGICAL                                            :: local_retain_sparsity, smearing
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_correction
      TYPE(dbcsr_type)                                   :: tmp

      CALL timeset(routineN, handle)

      IF (.NOT. PRESENT(retain_overlap_sparsity)) THEN
         local_retain_sparsity = .FALSE.
      ELSE
         local_retain_sparsity = retain_overlap_sparsity
      ENDIF

      IF (.NOT. PRESENT(smear)) THEN
         smearing = .FALSE.
      ELSE
         smearing = smear
      ENDIF

      CALL dbcsr_create(tmp, template=ket, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! TMP=metric*ket
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          metric, ket, 0.0_dp, tmp, &
                          filter_eps=eps_filter)

      ! OVERLAP=tr(bra)*TMP
      CALL dbcsr_multiply("T", "N", 1.0_dp, &
                          bra, tmp, 0.0_dp, overlap, &
                          retain_sparsity=local_retain_sparsity, &
                          filter_eps=eps_filter)

      CALL dbcsr_release(tmp)

      !! If smearing ALMO is requested, apply correction of the occupation-rescaling trick
      !! (i.e. converting rescaled orbitals into selfish orbitals)
      !! (i.e. the diagonal blocks of the rescaled overlap must remain unscaled)
      !! Since we have here orthonormal MOs within a fragment, diagonal blocks are identity matrices
      !! Therefore, one only need to restore the diagonal to 1
      !! RS-WARNING: Assume orthonormal MOs within a fragment
      IF (smearing) THEN
         CALL dbcsr_get_info(overlap, nfullrows_total=dim0)
         ALLOCATE (diag_correction(dim0))
         diag_correction = 1.0_dp
         CALL dbcsr_set_diag(overlap, diag_correction)
         DEALLOCATE (diag_correction)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_overlap

! **************************************************************************************************
!> \brief orthogonalize MOs
!> \param ket ...
!> \param overlap ...
!> \param metric ...
!> \param retain_locality ...
!> \param only_normalize ...
!> \param nocc_of_domain ...
!> \param eps_filter ...
!> \param order_lanczos ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param overlap_sqrti ...
!> \param smear ...
!> \par History
!>       2012.03 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE orthogonalize_mos(ket, overlap, metric, retain_locality, only_normalize, &
                                nocc_of_domain, eps_filter, order_lanczos, eps_lanczos, &
                                max_iter_lanczos, overlap_sqrti, smear)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: ket, overlap
      TYPE(dbcsr_type), INTENT(IN)                       :: metric
      LOGICAL, INTENT(IN)                                :: retain_locality, only_normalize
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc_of_domain
      REAL(KIND=dp)                                      :: eps_filter
      INTEGER, INTENT(IN)                                :: order_lanczos
      REAL(KIND=dp), INTENT(IN)                          :: eps_lanczos
      INTEGER, INTENT(IN)                                :: max_iter_lanczos
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: overlap_sqrti
      LOGICAL, INTENT(IN), OPTIONAL                      :: smear

      CHARACTER(LEN=*), PARAMETER :: routineN = 'orthogonalize_mos', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dim0, handle
      LOGICAL                                            :: smearing
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diagonal
      TYPE(dbcsr_type)                                   :: matrix_sigma_blk_sqrt, &
                                                            matrix_sigma_blk_sqrt_inv, &
                                                            matrix_t_blk_tmp

      CALL timeset(routineN, handle)

      IF (.NOT. PRESENT(smear)) THEN
         smearing = .FALSE.
      ELSE
         smearing = smear
      ENDIF

      ! create block-diagonal sparsity pattern for the overlap
      ! in case retain_locality is set to true
      ! RZK-warning this will fail if distribution blocks are smaller than domains!!!
      CALL dbcsr_set(overlap, 0.0_dp)
      CALL dbcsr_add_on_diag(overlap, 1.0_dp)
      CALL dbcsr_filter(overlap, eps_filter)

      CALL get_overlap(ket, ket, overlap, metric, retain_locality, &
                       eps_filter, smear=smearing)

      IF (only_normalize) THEN

         CALL dbcsr_get_info(overlap, nfullrows_total=dim0)
         ALLOCATE (diagonal(dim0))
         CALL dbcsr_get_diag(overlap, diagonal)
         CALL dbcsr_set(overlap, 0.0_dp)
         CALL dbcsr_set_diag(overlap, diagonal)
         DEALLOCATE (diagonal)
         CALL dbcsr_filter(overlap, eps_filter)

      ENDIF

      CALL dbcsr_create(matrix_sigma_blk_sqrt, template=overlap, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_sigma_blk_sqrt_inv, template=overlap, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! compute sqrt and sqrt_inv of the blocked MO overlap
      CALL set_zero_electron_blocks_in_mo_mo_matrix(overlap, nocc_of_domain, 1.0_dp)
      CALL matrix_sqrt_Newton_Schulz(matrix_sigma_blk_sqrt, matrix_sigma_blk_sqrt_inv, &
                                     overlap, threshold=eps_filter, &
                                     order=order_lanczos, &
                                     eps_lanczos=eps_lanczos, &
                                     max_iter_lanczos=max_iter_lanczos)
      CALL set_zero_electron_blocks_in_mo_mo_matrix(overlap, nocc_of_domain, 0.0_dp)
      !CALL set_zero_electron_blocks_in_mo_mo_matrix(matrix_sigma_blk_sqrt,nocc_of_domain,0.0_dp)
      CALL set_zero_electron_blocks_in_mo_mo_matrix(matrix_sigma_blk_sqrt_inv, nocc_of_domain, 0.0_dp)

      CALL dbcsr_create(matrix_t_blk_tmp, &
                        template=ket, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          ket, &
                          matrix_sigma_blk_sqrt_inv, &
                          0.0_dp, matrix_t_blk_tmp, &
                          filter_eps=eps_filter)

      ! update the orbitals with the orthonormalized MOs
      CALL dbcsr_copy(ket, matrix_t_blk_tmp)

      ! return overlap SQRT_INV if necessary
      IF (PRESENT(overlap_sqrti)) THEN
         CALL dbcsr_copy(overlap_sqrti, &
                         matrix_sigma_blk_sqrt_inv)
      ENDIF

      CALL dbcsr_release(matrix_t_blk_tmp)
      CALL dbcsr_release(matrix_sigma_blk_sqrt)
      CALL dbcsr_release(matrix_sigma_blk_sqrt_inv)

      CALL timestop(handle)

   END SUBROUTINE orthogonalize_mos

! **************************************************************************************************
!> \brief computes the idempotent density matrix from MOs
!>        MOs can be either orthogonal or non-orthogonal
!> \param t ...
!> \param p ...
!> \param eps_filter ...
!> \param orthog_orbs ...
!> \param nocc_of_domain ...
!> \param s ...
!> \param sigma ...
!> \param sigma_inv ...
!> \param use_guess ...
!> \param smear ...
!> \param algorithm to inver sigma: 0 - Hotelling (linear), 1 - Cholesky (cubic, low prefactor)
!> \param para_env ...
!> \param blacs_env ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param inverse_accelerator ...
!> \param inv_eps_factor ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!>       2018.09 smearing support [Ruben Staub]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE almo_scf_t_to_proj(t, p, eps_filter, orthog_orbs, nocc_of_domain, s, sigma, sigma_inv, &
                                 use_guess, smear, algorithm, para_env, blacs_env, eps_lanczos, &
                                 max_iter_lanczos, inverse_accelerator, inv_eps_factor)

      TYPE(dbcsr_type), INTENT(IN)                       :: t
      TYPE(dbcsr_type), INTENT(INOUT)                    :: p
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      LOGICAL, INTENT(IN)                                :: orthog_orbs
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: nocc_of_domain
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: s
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: sigma, sigma_inv
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_guess, smear
      INTEGER, INTENT(IN), OPTIONAL                      :: algorithm
      TYPE(cp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(cp_blacs_env_type), OPTIONAL, POINTER         :: blacs_env
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_lanczos
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iter_lanczos, inverse_accelerator
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: inv_eps_factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_t_to_proj', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: dim0, handle, my_accelerator, &
                                                            my_algorithm
      LOGICAL                                            :: smearing, use_sigma_inv_guess
      REAL(KIND=dp)                                      :: my_inv_eps_factor
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_correction
      TYPE(dbcsr_type)                                   :: t_tmp

      CALL timeset(routineN, handle)

      ! make sure that S, sigma and sigma_inv are present for non-orthogonal orbitals
      IF (.NOT. orthog_orbs) THEN
         IF ((.NOT. PRESENT(s)) .OR. (.NOT. PRESENT(sigma)) .OR. &
             (.NOT. PRESENT(sigma_inv)) .OR. (.NOT. PRESENT(nocc_of_domain))) THEN
            CPABORT("Nonorthogonal orbitals need more input")
         ENDIF
      ENDIF

      my_algorithm = 0
      IF (PRESENT(algorithm)) my_algorithm = algorithm

      IF (my_algorithm == 1 .AND. (.NOT. PRESENT(para_env) .OR. .NOT. PRESENT(blacs_env))) &
         CPABORT("PARA and BLACS env are necessary for cholesky algorithm")

      use_sigma_inv_guess = .FALSE.
      IF (PRESENT(use_guess)) THEN
         use_sigma_inv_guess = use_guess
      ENDIF

      IF (.NOT. PRESENT(smear)) THEN
         smearing = .FALSE.
      ELSE
         smearing = smear
      ENDIF

      my_accelerator = 1
      IF (PRESENT(inverse_accelerator)) my_accelerator = inverse_accelerator

      my_inv_eps_factor = 10.0_dp
      IF (PRESENT(inv_eps_factor)) my_inv_eps_factor = inv_eps_factor

      IF (orthog_orbs) THEN

         CALL dbcsr_multiply("N", "T", 1.0_dp, t, t, &
                             0.0_dp, p, filter_eps=eps_filter)

      ELSE

         CALL dbcsr_create(t_tmp, template=t)

         ! TMP=S.T
         CALL dbcsr_multiply("N", "N", 1.0_dp, s, t, 0.0_dp, t_tmp, &
                             filter_eps=eps_filter)

         ! Sig=tr(T).TMP - get MO overlap
         CALL dbcsr_multiply("T", "N", 1.0_dp, t, t_tmp, 0.0_dp, sigma, &
                             filter_eps=eps_filter)

         !! If smearing ALMO is requested, apply correction of the occupation-rescaling trick
         !! (i.e. converting rescaled orbitals into selfish orbitals)
         !! (i.e. the diagonal blocks of the rescaled overlap must remain unscaled)
         !! Since we have here orthonormal MOs within a fragment, diagonal blocks are identity matrices
         !! Therefore, one only need to restore the diagonal to 1
         !! RS-WARNING: Assume orthonormal MOs within a fragment
         IF (smearing) THEN
            CALL dbcsr_get_info(sigma, nfullrows_total=dim0)
            ALLOCATE (diag_correction(dim0))
            diag_correction = 1.0_dp
            CALL dbcsr_set_diag(sigma, diag_correction)
            DEALLOCATE (diag_correction)
         END IF

         ! invert MO overlap
         CALL set_zero_electron_blocks_in_mo_mo_matrix(sigma, nocc_of_domain, 1.0_dp)
         SELECT CASE (my_algorithm)
         CASE (spd_inversion_ls_taylor)

            CALL invert_Taylor( &
               matrix_inverse=sigma_inv, &
               matrix=sigma, &
               use_inv_as_guess=use_sigma_inv_guess, &
               threshold=eps_filter*my_inv_eps_factor, &
               filter_eps=eps_filter, &
               !accelerator_order=my_accelerator, &
               !eps_lanczos=eps_lanczos, &
               !max_iter_lanczos=max_iter_lanczos, &
               silent=.FALSE.)

         CASE (spd_inversion_ls_hotelling)

            CALL invert_Hotelling( &
               matrix_inverse=sigma_inv, &
               matrix=sigma, &
               use_inv_as_guess=use_sigma_inv_guess, &
               threshold=eps_filter*my_inv_eps_factor, &
               filter_eps=eps_filter, &
               accelerator_order=my_accelerator, &
               eps_lanczos=eps_lanczos, &
               max_iter_lanczos=max_iter_lanczos, &
               silent=.FALSE.)

         CASE (spd_inversion_dense_cholesky)

            ! invert using cholesky
            CALL dbcsr_copy(sigma_inv, sigma)
            CALL cp_dbcsr_cholesky_decompose(sigma_inv, &
                                             para_env=para_env, &
                                             blacs_env=blacs_env)
            CALL cp_dbcsr_cholesky_invert(sigma_inv, &
                                          para_env=para_env, &
                                          blacs_env=blacs_env, &
                                          upper_to_full=.TRUE.)
            CALL dbcsr_filter(sigma_inv, &
                              eps=eps_filter)
         CASE DEFAULT
            CPABORT("Illegal MO overalp inversion algorithm")
         END SELECT
         CALL set_zero_electron_blocks_in_mo_mo_matrix(sigma, nocc_of_domain, 0.0_dp)
         CALL set_zero_electron_blocks_in_mo_mo_matrix(sigma_inv, nocc_of_domain, 0.0_dp)

         ! TMP=T.SigInv
         CALL dbcsr_multiply("N", "N", 1.0_dp, t, sigma_inv, 0.0_dp, t_tmp, &
                             filter_eps=eps_filter)

         ! P=TMP.tr(T_blk)
         CALL dbcsr_multiply("N", "T", 1.0_dp, t_tmp, t, 0.0_dp, p, &
                             filter_eps=eps_filter)

         CALL dbcsr_release(t_tmp)

      ENDIF

      CALL timestop(handle)

   END SUBROUTINE almo_scf_t_to_proj

! **************************************************************************************************
!> \brief self-explanatory
!> \param matrix ...
!> \param nocc_of_domain ...
!> \param value ...
!> \param
!> \par History
!>       2016.12 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE set_zero_electron_blocks_in_mo_mo_matrix(matrix, nocc_of_domain, value)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nocc_of_domain
      REAL(KIND=dp), INTENT(IN)                          :: value

      INTEGER                                            :: hold, iblock_col, iblock_row, mynode, &
                                                            nblkrows_tot, row
      LOGICAL                                            :: found, tr
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_new_block
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL dbcsr_get_info(matrix, distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)
      !mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix)))
      CALL dbcsr_work_create(matrix, work_mutable=.TRUE.)

      nblkrows_tot = dbcsr_nblkrows_total(matrix)

      DO row = 1, nblkrows_tot
         IF (nocc_of_domain(row) == 0) THEN
            tr = .FALSE.
            iblock_row = row
            iblock_col = row
            CALL dbcsr_get_stored_coordinates(matrix, iblock_row, iblock_col, hold)
            IF (hold .EQ. mynode) THEN
               NULLIFY (p_new_block)
               CALL dbcsr_get_block_p(matrix, iblock_row, iblock_col, p_new_block, found)
               IF (found) THEN
                  p_new_block(1, 1) = value
               ELSE
                  CALL dbcsr_reserve_block2d(matrix, iblock_row, iblock_col, p_new_block)
                  CPASSERT(ASSOCIATED(p_new_block))
                  p_new_block(1, 1) = value
               ENDIF
            ENDIF ! mynode
         ENDIF !zero-electron block
      ENDDO

      CALL dbcsr_finalize(matrix)

   END SUBROUTINE set_zero_electron_blocks_in_mo_mo_matrix

! **************************************************************************************************
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
!> \par History
!>       2011.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE apply_projector(psi_in, psi_out, psi_projector, metric, project_out, &
                              psi_projector_orthogonal, proj_in_template, eps_filter, sig_inv_projector, &
                              sig_inv_template)

      TYPE(dbcsr_type), INTENT(IN)                       :: psi_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: psi_out
      TYPE(dbcsr_type), INTENT(IN)                       :: psi_projector, metric
      LOGICAL, INTENT(IN)                                :: project_out, psi_projector_orthogonal
      TYPE(dbcsr_type), INTENT(IN)                       :: proj_in_template
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: sig_inv_projector, sig_inv_template

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_projector', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: tmp_no, tmp_ov, tmp_ov2, tmp_sig, &
                                                            tmp_sig_inv

      CALL timeset(routineN, handle)

      ! =S*PSI_proj
      CALL dbcsr_create(tmp_no, template=psi_projector)
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          metric, psi_projector, &
                          0.0_dp, tmp_no, &
                          filter_eps=eps_filter)

      ! =tr(S.PSI_proj)*PSI_in
      CALL dbcsr_create(tmp_ov, template=proj_in_template)
      CALL dbcsr_multiply("T", "N", 1.0_dp, &
                          tmp_no, psi_in, &
                          0.0_dp, tmp_ov, &
                          filter_eps=eps_filter)

      IF (.NOT. psi_projector_orthogonal) THEN
         ! =SigInv_proj*Sigma_OV
         CALL dbcsr_create(tmp_ov2, &
                           template=proj_in_template)
         IF (PRESENT(sig_inv_projector)) THEN
            CALL dbcsr_create(tmp_sig_inv, &
                              template=sig_inv_projector)
            CALL dbcsr_copy(tmp_sig_inv, sig_inv_projector)
         ELSE
            IF (.NOT. PRESENT(sig_inv_template)) THEN
               CPABORT("PROGRAMMING ERROR: provide either template or sig_inv")
            ENDIF
            ! compute inverse overlap of the projector orbitals
            CALL dbcsr_create(tmp_sig, &
                              template=sig_inv_template, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                psi_projector, tmp_no, 0.0_dp, tmp_sig, &
                                filter_eps=eps_filter)
            CALL dbcsr_create(tmp_sig_inv, &
                              template=sig_inv_template, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL invert_Hotelling(tmp_sig_inv, tmp_sig, &
                                  threshold=eps_filter)
            CALL dbcsr_release(tmp_sig)
         ENDIF
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             tmp_sig_inv, tmp_ov, 0.0_dp, tmp_ov2, &
                             filter_eps=eps_filter)
         CALL dbcsr_release(tmp_sig_inv)
         CALL dbcsr_copy(tmp_ov, tmp_ov2)
         CALL dbcsr_release(tmp_ov2)
      ENDIF
      CALL dbcsr_release(tmp_no)

      ! =PSI_proj*TMP_OV
      CALL dbcsr_multiply("N", "N", 1.0_dp, &
                          psi_projector, tmp_ov, 0.0_dp, psi_out, &
                          filter_eps=eps_filter)
      CALL dbcsr_release(tmp_ov)

      ! V_out=V_in-V_out
      IF (project_out) THEN
         CALL dbcsr_add(psi_out, psi_in, -1.0_dp, +1.0_dp)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE apply_projector

!! **************************************************************************************************
!!> \brief projects the occupied space out from the provided orbitals
!!> \par History
!!>       2011.07 created [Rustam Z Khaliullin]
!!> \author Rustam Z Khaliullin
!! **************************************************************************************************
!  SUBROUTINE almo_scf_p_out_from_v(v_in,v_out,ov_template,ispin,almo_scf_env)
!
!    TYPE(dbcsr_type), INTENT(IN)                :: v_in, ov_template
!    TYPE(dbcsr_type), INTENT(INOUT)             :: v_out
!    INTEGER, INTENT(IN)                            :: ispin
!    TYPE(almo_scf_env_type), INTENT(INOUT)         :: almo_scf_env
!
!    CHARACTER(LEN=*), PARAMETER :: &
!      routineN = 'almo_scf_p_out_from_v', &
!      routineP = moduleN//':'//routineN
!
!    TYPE(dbcsr_type)                      :: tmp_on, tmp_ov, tmp_ov2
!    INTEGER                                  :: handle
!    LOGICAL                                  :: failure
!
!    CALL timeset(routineN,handle)
!
!    ! =tr(T_blk)*S
!    CALL dbcsr_init(tmp_on)
!    CALL dbcsr_create(tmp_on,&
!            template=almo_scf_env%matrix_t_tr(ispin))
!    CALL dbcsr_multiply("T","N",1.0_dp,&
!            almo_scf_env%matrix_t_blk(ispin),&
!            almo_scf_env%matrix_s(1),&
!            0.0_dp,tmp_on,&
!            filter_eps=almo_scf_env%eps_filter)
!
!    ! =tr(T_blk).S*V_in
!    CALL dbcsr_init(tmp_ov)
!    CALL dbcsr_create(tmp_ov,template=ov_template)
!    CALL dbcsr_multiply("N","N",1.0_dp,&
!            tmp_on,v_in,0.0_dp,tmp_ov,&
!            filter_eps=almo_scf_env%eps_filter)
!    CALL dbcsr_release(tmp_on)
!
!    ! =SigmaInv*Sigma_OV
!    CALL dbcsr_init(tmp_ov2)
!    CALL dbcsr_create(tmp_ov2,template=ov_template)
!    CALL dbcsr_multiply("N","N",1.0_dp,&
!            almo_scf_env%matrix_sigma_inv(ispin),&
!            tmp_ov,0.0_dp,tmp_ov2,&
!            filter_eps=almo_scf_env%eps_filter)
!    CALL dbcsr_release(tmp_ov)
!
!    ! =T_blk*SigmaInv.Sigma_OV
!    CALL dbcsr_multiply("N","N",1.0_dp,&
!            almo_scf_env%matrix_t_blk(ispin),&
!            tmp_ov2,0.0_dp,v_out,&
!            filter_eps=almo_scf_env%eps_filter)
!    CALL dbcsr_release(tmp_ov2)
!
!    ! V_out=V_in-V_out=
!    CALL dbcsr_add(v_out,v_in,-1.0_dp,+1.0_dp)
!
!    CALL timestop(handle)
!
!  END SUBROUTINE almo_scf_p_out_from_v

! **************************************************************************************************
!> \brief computes a unitary matrix from an arbitrary "generator" matrix
!>        U = ( 1 - X + tr(X) ) ( 1 + X - tr(X) )^(-1)
!> \param X ...
!> \param U ...
!> \param eps_filter ...
!> \par History
!>       2011.08 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE generator_to_unitary(X, U, eps_filter)

      TYPE(dbcsr_type), INTENT(IN)                       :: X
      TYPE(dbcsr_type), INTENT(INOUT)                    :: U
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'generator_to_unitary', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      LOGICAL                                            :: safe_mode
      REAL(KIND=dp)                                      :: frob_matrix, frob_matrix_base
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: delta, t1, t2, tmp1

      CALL timeset(routineN, handle)

      safe_mode = .TRUE.

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(t1, template=X, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(t2, template=X, &
                        matrix_type=dbcsr_type_no_symmetry)

      ! create antisymmetric Delta = -X + tr(X)
      CALL dbcsr_create(delta, template=X, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_transposed(delta, X)
! check that transposed is added correctly
      CALL dbcsr_add(delta, X, 1.0_dp, -1.0_dp)

      ! compute (1 - Delta)^(-1)
      CALL dbcsr_add_on_diag(t1, 1.0_dp)
      CALL dbcsr_add(t1, delta, 1.0_dp, -1.0_dp)
      CALL invert_Hotelling(t2, t1, threshold=eps_filter)

      IF (safe_mode) THEN

         CALL dbcsr_create(tmp1, template=X, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "N", 1.0_dp, t2, t1, 0.0_dp, tmp1, &
                             filter_eps=eps_filter)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)
         IF (unit_nr > 0) WRITE (unit_nr, *) "Error for (inv(A)*A-I)", frob_matrix/frob_matrix_base
         CALL dbcsr_release(tmp1)
      ENDIF

      CALL dbcsr_multiply("N", "N", 1.0_dp, delta, t2, 0.0_dp, U, &
                          filter_eps=eps_filter)
      CALL dbcsr_add(U, t2, 1.0_dp, 1.0_dp)

      IF (safe_mode) THEN

         CALL dbcsr_create(tmp1, template=X, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("T", "N", 1.0_dp, U, U, 0.0_dp, tmp1, &
                             filter_eps=eps_filter)
         frob_matrix_base = dbcsr_frobenius_norm(tmp1)
         CALL dbcsr_add_on_diag(tmp1, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(tmp1)
         IF (unit_nr > 0) WRITE (unit_nr, *) "Error for (trn(U)*U-I)", frob_matrix/frob_matrix_base
         CALL dbcsr_release(tmp1)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE generator_to_unitary

! **************************************************************************************************
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
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE apply_domain_operators(matrix_in, matrix_out, operator1, operator2, &
                                     dpattern, map, node_of_domain, my_action, filter_eps, matrix_trimmer, use_trimmer)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_out
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN)                                      :: operator1
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN), OPTIONAL                            :: operator2
      TYPE(dbcsr_type), INTENT(IN)                       :: dpattern
      TYPE(domain_map_type), INTENT(IN)                  :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_of_domain
      INTEGER, INTENT(IN)                                :: my_action
      REAL(KIND=dp)                                      :: filter_eps
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: matrix_trimmer
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_trimmer

      CHARACTER(len=*), PARAMETER :: routineN = 'apply_domain_operators', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ndomains
      LOGICAL                                            :: matrix_trimmer_required, my_use_trimmer, &
                                                            operator2_required
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_in, subm_out, subm_temp

      CALL timeset(routineN, handle)

      my_use_trimmer = .FALSE.
      IF (PRESENT(use_trimmer)) THEN
         my_use_trimmer = use_trimmer
      ENDIF

      operator2_required = .FALSE.
      matrix_trimmer_required = .FALSE.

      IF (my_action .EQ. 1) operator2_required = .TRUE.

      IF (my_use_trimmer) THEN
         matrix_trimmer_required = .TRUE.
         CPABORT("TRIMMED PROJECTOR DISABLED!")
      ENDIF

      IF (.NOT. PRESENT(operator2) .AND. operator2_required) THEN
         CPABORT("SECOND OPERATOR IS REQUIRED")
      ENDIF
      IF (.NOT. PRESENT(matrix_trimmer) .AND. matrix_trimmer_required) THEN
         CPABORT("TRIMMER MATRIX IS REQUIRED")
      ENDIF

      ndomains = dbcsr_nblkcols_total(dpattern)

      ALLOCATE (subm_in(ndomains))
      ALLOCATE (subm_temp(ndomains))
      ALLOCATE (subm_out(ndomains))
      !!!TRIM ALLOCATE(subm_trimmer(ndomains))
      CALL init_submatrices(subm_in)
      CALL init_submatrices(subm_temp)
      CALL init_submatrices(subm_out)

      CALL construct_submatrices(matrix_in, subm_in, &
                                 dpattern, map, node_of_domain, select_row)

      !!!TRIM IF (matrix_trimmer_required) THEN
      !!!TRIM    CALL construct_submatrices(matrix_trimmer,subm_trimmer,&
      !!!TRIM            dpattern,map,node_of_domain,select_row)
      !!!TRIM ENDIF

      IF (my_action .EQ. 0) THEN
         ! for example, apply preconditioner
         CALL multiply_submatrices('N', 'N', 1.0_dp, operator1, &
                                   subm_in, 0.0_dp, subm_out)
      ELSE IF (my_action .EQ. 1) THEN
         ! use for projectors
         CALL copy_submatrices(subm_in, subm_out, .TRUE.)
         CALL multiply_submatrices('N', 'N', 1.0_dp, operator1, &
                                   subm_in, 0.0_dp, subm_temp)
         CALL multiply_submatrices('N', 'N', -1.0_dp, operator2, &
                                   subm_temp, 1.0_dp, subm_out)
      ELSE
         CPABORT("ILLEGAL ACTION")
      ENDIF

      CALL construct_dbcsr_from_submatrices(matrix_out, subm_out, dpattern)
      CALL dbcsr_filter(matrix_out, filter_eps)

      CALL release_submatrices(subm_out)
      CALL release_submatrices(subm_temp)
      CALL release_submatrices(subm_in)

      DEALLOCATE (subm_out)
      DEALLOCATE (subm_temp)
      DEALLOCATE (subm_in)

      CALL timestop(handle)

   END SUBROUTINE apply_domain_operators

! **************************************************************************************************
!> \brief Constructs preconditioners for each domain
!>        -1. projected preconditioner
!>         0. simple preconditioner
!> \param matrix_main ...
!> \param subm_s_inv ...
!> \param subm_s_inv_half ...
!> \param subm_s_half ...
!> \param subm_r_down ...
!> \param matrix_trimmer ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param preconditioner ...
!> \param bad_modes_projector_down ...
!> \param use_trimmer ...
!> \param eps_zero_eigenvalues ...
!> \param my_action ...
!> \param skip_inversion ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_domain_preconditioner(matrix_main, subm_s_inv, subm_s_inv_half, &
                                          subm_s_half, subm_r_down, matrix_trimmer, dpattern, map, node_of_domain, preconditioner, &
                                             bad_modes_projector_down, use_trimmer, eps_zero_eigenvalues, my_action, skip_inversion)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_main
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(IN), OPTIONAL                            :: subm_s_inv, subm_s_inv_half, &
                                                            subm_s_half, subm_r_down
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: matrix_trimmer
      TYPE(dbcsr_type), INTENT(IN)                       :: dpattern
      TYPE(domain_map_type), INTENT(IN)                  :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_of_domain
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: preconditioner
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT), OPTIONAL                         :: bad_modes_projector_down
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_trimmer
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_zero_eigenvalues
      INTEGER, INTENT(IN)                                :: my_action
      LOGICAL, INTENT(IN), OPTIONAL                      :: skip_inversion

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_preconditioner', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idomain, index1_end, &
                                                            index1_start, n_domain_mos, naos, &
                                                            nblkrows_tot, ndomains, neighbor, row
      INTEGER, DIMENSION(:), POINTER                     :: nmos
      LOGICAL :: eps_zero_eigenvalues_required, matrix_r_required, matrix_s_half_required, &
         matrix_s_inv_half_required, matrix_s_inv_required, matrix_trimmer_required, &
         my_skip_inversion, my_use_trimmer
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Minv, proj_array
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_main, subm_tmp, subm_tmp2

      CALL timeset(routineN, handle)

      my_use_trimmer = .FALSE.
      IF (PRESENT(use_trimmer)) THEN
         my_use_trimmer = use_trimmer
      ENDIF

      my_skip_inversion = .FALSE.
      IF (PRESENT(skip_inversion)) THEN
         my_skip_inversion = skip_inversion
      ENDIF

      matrix_s_inv_half_required = .FALSE.
      matrix_s_half_required = .FALSE.
      eps_zero_eigenvalues_required = .FALSE.
      matrix_s_inv_required = .FALSE.
      matrix_trimmer_required = .FALSE.
      matrix_r_required = .FALSE.

      IF (my_action .EQ. -1) matrix_s_inv_required = .TRUE.
      IF (my_action .EQ. -1) matrix_r_required = .TRUE.
      IF (my_use_trimmer) THEN
         matrix_trimmer_required = .TRUE.
         CPABORT("TRIMMED PRECONDITIONER DISABLED!")
      ENDIF
      ! tie the following optional arguments together to prevent bad calls
      IF (PRESENT(bad_modes_projector_down)) THEN
         matrix_s_inv_half_required = .TRUE.
         matrix_s_half_required = .TRUE.
         eps_zero_eigenvalues_required = .TRUE.
      ENDIF

      ! check if all required optional arguments are provided
      IF (.NOT. PRESENT(subm_s_inv_half) .AND. matrix_s_inv_half_required) THEN
         CPABORT("S_inv_half SUBMATRICES ARE REQUIRED")
      ENDIF
      IF (.NOT. PRESENT(subm_s_half) .AND. matrix_s_half_required) THEN
         CPABORT("S_half SUBMATRICES ARE REQUIRED")
      ENDIF
      IF (.NOT. PRESENT(eps_zero_eigenvalues) .AND. eps_zero_eigenvalues_required) THEN
         CPABORT("EPS_ZERO_EIGENVALUES IS REQUIRED")
      ENDIF
      IF (.NOT. PRESENT(subm_s_inv) .AND. matrix_s_inv_required) THEN
         CPABORT("S_inv SUBMATRICES ARE REQUIRED")
      ENDIF
      IF (.NOT. PRESENT(subm_r_down) .AND. matrix_r_required) THEN
         CPABORT("R SUBMATRICES ARE REQUIRED")
      ENDIF
      IF (.NOT. PRESENT(matrix_trimmer) .AND. matrix_trimmer_required) THEN
         CPABORT("TRIMMER MATRIX IS REQUIRED")
      ENDIF

      CALL dbcsr_get_info(dpattern, &
                          nblkcols_total=ndomains, &
                          nblkrows_total=nblkrows_tot, &
                          col_blk_size=nmos)

      ALLOCATE (subm_main(ndomains))
      CALL init_submatrices(subm_main)
      !!!TRIM ALLOCATE(subm_trimmer(ndomains))

      CALL construct_submatrices(matrix_main, subm_main, &
                                 dpattern, map, node_of_domain, select_row_col)

      !!!TRIM IF (matrix_trimmer_required) THEN
      !!!TRIM    CALL construct_submatrices(matrix_trimmer,subm_trimmer,&
      !!!TRIM            dpattern,map,node_of_domain,select_row)
      !!!TRIM ENDIF

      IF (my_action .EQ. -1) THEN
         ! project out the local occupied space
         !tmp=MATMUL(subm_r(idomain)%mdata,Minv)
         !Minv=MATMUL(tmp,subm_main(idomain)%mdata)
         !subm_main(idomain)%mdata=subm_main(idomain)%mdata-&
         !   Minv-TRANSPOSE(Minv)+MATMUL(Minv,TRANSPOSE(tmp))
         ALLOCATE (subm_tmp(ndomains))
         ALLOCATE (subm_tmp2(ndomains))
         CALL init_submatrices(subm_tmp)
         CALL init_submatrices(subm_tmp2)
         CALL multiply_submatrices('N', 'N', 1.0_dp, subm_r_down, &
                                   subm_s_inv, 0.0_dp, subm_tmp)
         CALL multiply_submatrices('N', 'N', 1.0_dp, subm_tmp, &
                                   subm_main, 0.0_dp, subm_tmp2)
         CALL add_submatrices(1.0_dp, subm_main, -1.0_dp, subm_tmp2, 'N')
         CALL add_submatrices(1.0_dp, subm_main, -1.0_dp, subm_tmp2, 'T')
         CALL multiply_submatrices('N', 'T', 1.0_dp, subm_tmp2, &
                                   subm_tmp, 1.0_dp, subm_main)
         CALL release_submatrices(subm_tmp)
         CALL release_submatrices(subm_tmp2)
         DEALLOCATE (subm_tmp2)
         DEALLOCATE (subm_tmp)
      ENDIF

      IF (my_skip_inversion) THEN
         CALL copy_submatrices(subm_main, preconditioner, .TRUE.)
      ELSE
         ! loop over domains - perform inversion
         DO idomain = 1, ndomains

            ! check if the submatrix exists
            IF (subm_main(idomain)%domain .GT. 0) THEN

               ! find sizes of MO submatrices
               IF (idomain .EQ. 1) THEN
                  index1_start = 1
               ELSE
                  index1_start = map%index1(idomain - 1)
               ENDIF
               index1_end = map%index1(idomain) - 1

               n_domain_mos = 0
               DO row = index1_start, index1_end
                  neighbor = map%pairs(row, 1)
                  n_domain_mos = n_domain_mos + nmos(neighbor)
               ENDDO

               naos = subm_main(idomain)%nrows
               !WRITE(*,*) "Domain, mo_self_and_neig, ao_domain: ", idomain, n_domain_mos, naos

               ALLOCATE (Minv(naos, naos))

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
               !!!TRIM
               !!!TRIM       ! apply the inverted matrix
               !!!TRIM       ! RZK-warning this is only possible when the preconditioner is applied
               !!!TRIM       tmp(:,ii)=MATMUL(Minv,subm_in(idomain)%mdata(:,ii))
               !!!TRIM    ENDDO
               !!!TRIM    subm_out=MATMUL(tmp,sigma)
               !!!TRIM    deallocate(tmp)
               !!!TRIM ELSE

               IF (PRESENT(bad_modes_projector_down)) THEN
                  ALLOCATE (proj_array(naos, naos))
                  CALL pseudo_invert_matrix(A=subm_main(idomain)%mdata, Ainv=Minv, N=naos, method=1, &
                                            range1=nmos(idomain), range2=n_domain_mos, &
                                            range1_thr=eps_zero_eigenvalues, &
                                            bad_modes_projector_down=proj_array, &
                                            s_inv_half=subm_s_inv_half(idomain)%mdata, &
                                            s_half=subm_s_half(idomain)%mdata &
                                            )
               ELSE
                  CALL pseudo_invert_matrix(A=subm_main(idomain)%mdata, Ainv=Minv, N=naos, method=1, &
                                            range1=nmos(idomain), range2=n_domain_mos)
               ENDIF
               !!!TRIM ENDIF

               CALL copy_submatrices(subm_main(idomain), preconditioner(idomain), .FALSE.)
               CALL copy_submatrix_data(Minv, preconditioner(idomain))
               DEALLOCATE (Minv)

               IF (PRESENT(bad_modes_projector_down)) THEN
                  CALL copy_submatrices(subm_main(idomain), bad_modes_projector_down(idomain), .FALSE.)
                  CALL copy_submatrix_data(proj_array, bad_modes_projector_down(idomain))
                  DEALLOCATE (proj_array)
               ENDIF

            ENDIF ! submatrix for the domain exists

         ENDDO ! loop over domains

      ENDIF

      CALL release_submatrices(subm_main)
      DEALLOCATE (subm_main)
      !DEALLOCATE(subm_s)
      !DEALLOCATE(subm_r)

      !IF (matrix_r_required) THEN
      !   CALL dbcsr_release(m_tmp_no_1)
      !   CALL dbcsr_release(m_tmp_no_2)
      !   CALL dbcsr_release(matrix_r)
      !ENDIF

      CALL timestop(handle)

   END SUBROUTINE construct_domain_preconditioner

! **************************************************************************************************
!> \brief Constructs S^(+1/2) and S^(-1/2) submatrices for each domain
!> \param matrix_s ...
!> \param subm_s_sqrt ...
!> \param subm_s_sqrt_inv ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_domain_s_sqrt(matrix_s, subm_s_sqrt, subm_s_sqrt_inv, &
                                      dpattern, map, node_of_domain)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: subm_s_sqrt, subm_s_sqrt_inv
      TYPE(dbcsr_type), INTENT(IN)                       :: dpattern
      TYPE(domain_map_type), INTENT(IN)                  :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_of_domain

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_s_sqrt', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idomain, naos, ndomains
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Ssqrt, Ssqrtinv
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_s

      CALL timeset(routineN, handle)

      ndomains = dbcsr_nblkcols_total(dpattern)
      CPASSERT(SIZE(subm_s_sqrt) .EQ. ndomains)
      CPASSERT(SIZE(subm_s_sqrt_inv) .EQ. ndomains)
      ALLOCATE (subm_s(ndomains))
      CALL init_submatrices(subm_s)

      CALL construct_submatrices(matrix_s, subm_s, &
                                 dpattern, map, node_of_domain, select_row_col)

      ! loop over domains - perform inversion
      DO idomain = 1, ndomains

         ! check if the submatrix exists
         IF (subm_s(idomain)%domain .GT. 0) THEN

            naos = subm_s(idomain)%nrows

            ALLOCATE (Ssqrt(naos, naos))
            ALLOCATE (Ssqrtinv(naos, naos))

            CALL matrix_sqrt(A=subm_s(idomain)%mdata, Asqrt=Ssqrt, Asqrtinv=Ssqrtinv, &
                             N=naos)

            CALL copy_submatrices(subm_s(idomain), subm_s_sqrt(idomain), .FALSE.)
            CALL copy_submatrix_data(Ssqrt, subm_s_sqrt(idomain))

            CALL copy_submatrices(subm_s(idomain), subm_s_sqrt_inv(idomain), .FALSE.)
            CALL copy_submatrix_data(Ssqrtinv, subm_s_sqrt_inv(idomain))

            DEALLOCATE (Ssqrtinv)
            DEALLOCATE (Ssqrt)

         ENDIF ! submatrix for the domain exists

      ENDDO ! loop over domains

      CALL release_submatrices(subm_s)
      DEALLOCATE (subm_s)

      CALL timestop(handle)

   END SUBROUTINE construct_domain_s_sqrt

! **************************************************************************************************
!> \brief Constructs S_inv block for each domain
!> \param matrix_s ...
!> \param subm_s_inv ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_domain_s_inv(matrix_s, subm_s_inv, dpattern, map, &
                                     node_of_domain)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: subm_s_inv
      TYPE(dbcsr_type), INTENT(IN)                       :: dpattern
      TYPE(domain_map_type), INTENT(IN)                  :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_of_domain

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_s_inv', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idomain, naos, ndomains
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Sinv
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_s

      CALL timeset(routineN, handle)

      ndomains = dbcsr_nblkcols_total(dpattern)

      CPASSERT(SIZE(subm_s_inv) .EQ. ndomains)
      ALLOCATE (subm_s(ndomains))
      CALL init_submatrices(subm_s)

      CALL construct_submatrices(matrix_s, subm_s, &
                                 dpattern, map, node_of_domain, select_row_col)

      ! loop over domains - perform inversion
      DO idomain = 1, ndomains

         ! check if the submatrix exists
         IF (subm_s(idomain)%domain .GT. 0) THEN

            naos = subm_s(idomain)%nrows

            ALLOCATE (Sinv(naos, naos))

            CALL pseudo_invert_matrix(A=subm_s(idomain)%mdata, Ainv=Sinv, N=naos, &
                                      method=0)

            CALL copy_submatrices(subm_s(idomain), subm_s_inv(idomain), .FALSE.)
            CALL copy_submatrix_data(Sinv, subm_s_inv(idomain))

            DEALLOCATE (Sinv)

         ENDIF ! submatrix for the domain exists

      ENDDO ! loop over domains

      CALL release_submatrices(subm_s)
      DEALLOCATE (subm_s)

      CALL timestop(handle)

   END SUBROUTINE construct_domain_s_inv

! **************************************************************************************************
!> \brief Constructs subblocks of the covariant-covariant projectors (i.e. DM without spin factor)
!> \param matrix_t ...
!> \param matrix_sigma_inv ...
!> \param matrix_s ...
!> \param subm_r_down ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \param filter_eps ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_domain_r_down(matrix_t, matrix_sigma_inv, matrix_s, &
                                      subm_r_down, dpattern, map, node_of_domain, filter_eps)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_t, matrix_sigma_inv, matrix_s
      TYPE(domain_submatrix_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: subm_r_down
      TYPE(dbcsr_type), INTENT(IN)                       :: dpattern
      TYPE(domain_map_type), INTENT(IN)                  :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_of_domain
      REAL(KIND=dp)                                      :: filter_eps

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_domain_r_down', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ndomains
      TYPE(dbcsr_type)                                   :: m_tmp_no_1, m_tmp_no_2, matrix_r

      CALL timeset(routineN, handle)

      ! compute the density matrix in the COVARIANT representation
      CALL dbcsr_create(matrix_r, &
                        template=matrix_s, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(m_tmp_no_1, &
                        template=matrix_t, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(m_tmp_no_2, &
                        template=matrix_t, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s, matrix_t, &
                          0.0_dp, m_tmp_no_1, filter_eps=filter_eps)
      CALL dbcsr_multiply("N", "N", 1.0_dp, m_tmp_no_1, matrix_sigma_inv, &
                          0.0_dp, m_tmp_no_2, filter_eps=filter_eps)
      CALL dbcsr_multiply("N", "T", 1.0_dp, m_tmp_no_2, m_tmp_no_1, &
                          0.0_dp, matrix_r, filter_eps=filter_eps)

      CALL dbcsr_release(m_tmp_no_1)
      CALL dbcsr_release(m_tmp_no_2)

      ndomains = dbcsr_nblkcols_total(dpattern)
      CPASSERT(SIZE(subm_r_down) .EQ. ndomains)

      CALL construct_submatrices(matrix_r, subm_r_down, &
                                 dpattern, map, node_of_domain, select_row_col)

      CALL dbcsr_release(matrix_r)

      CALL timestop(handle)

   END SUBROUTINE construct_domain_r_down

! **************************************************************************************************
!> \brief Finds the square root of a matrix and its inverse
!> \param A ...
!> \param Asqrt ...
!> \param Asqrtinv ...
!> \param N ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE matrix_sqrt(A, Asqrt, Asqrtinv, N)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: A
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Asqrt, Asqrtinv
      INTEGER, INTENT(IN)                                :: N

      CHARACTER(len=*), PARAMETER :: routineN = 'matrix_sqrt', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, INFO, jj, LWORK, unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, WORK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: test, testN
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      Asqrtinv = A
      INFO = 0

      ! get a useful unit_nr
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      !CALL DPOTRF('L', N, Ainv, N, INFO )
      !IF( INFO.NE.0 ) THEN
      !   CPErrorMessage(cp_failure_level,routineP,"DPOTRF failed")
      !   CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
      !END IF
      !CALL DPOTRI('L', N, Ainv, N, INFO )
      !IF( INFO.NE.0 ) THEN
      !   CPErrorMessage(cp_failure_level,routineP,"DPOTRI failed")
      !   CPPrecondition(.FALSE.,cp_failure_level,routineP,failure)
      !END IF
      !! complete the matrix
      !DO ii=1,N
      !   DO jj=ii+1,N
      !      Ainv(ii,jj)=Ainv(jj,ii)
      !   ENDDO
      !   !WRITE(*,'(100F13.9)') Ainv(ii,:)
      !ENDDO

      ! diagonalize first
      ALLOCATE (eigenvalues(N))
      ! Query the optimal workspace for dsyev
      LWORK = -1
      ALLOCATE (WORK(MAX(1, LWORK)))
      CALL DSYEV('V', 'L', N, Asqrtinv, N, eigenvalues, WORK, LWORK, INFO)
      LWORK = INT(WORK(1))
      DEALLOCATE (WORK)
      ! Allocate the workspace and solve the eigenproblem
      ALLOCATE (WORK(MAX(1, LWORK)))
      CALL DSYEV('V', 'L', N, Asqrtinv, N, eigenvalues, WORK, LWORK, INFO)
      IF (INFO .NE. 0) THEN
         IF (unit_nr > 0) WRITE (unit_nr, *) 'DSYEV ERROR MESSAGE: ', INFO
         CPABORT("DSYEV failed")
      END IF
      DEALLOCATE (WORK)

      ! take functions of eigenvalues and use eigenvectors to compute the matrix function
      ! first sqrt
      ALLOCATE (test(N, N))
      DO jj = 1, N
         test(jj, :) = Asqrtinv(:, jj)*SQRT(eigenvalues(jj))
      ENDDO
      ALLOCATE (testN(N, N))
      testN(:, :) = MATMUL(Asqrtinv, test)
      Asqrt = testN
      ! now, sqrt_inv
      DO jj = 1, N
         test(jj, :) = Asqrtinv(:, jj)/SQRT(eigenvalues(jj))
      ENDDO
      testN(:, :) = MATMUL(Asqrtinv, test)
      Asqrtinv = testN
      DEALLOCATE (test, testN)

      DEALLOCATE (eigenvalues)

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

! **************************************************************************************************
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
!> \param shift ...
!> \param bad_modes_projector_down ...
!> \param s_inv_half ...
!> \param s_half ...
!> \par History
!>       2012.04 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE pseudo_invert_matrix(A, Ainv, N, method, range1, range2, range1_thr, &
                                   shift, bad_modes_projector_down, s_inv_half, s_half)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: A
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Ainv
      INTEGER, INTENT(IN)                                :: N, method
      INTEGER, INTENT(IN), OPTIONAL                      :: range1, range2
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: range1_thr, shift
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: bad_modes_projector_down
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: s_inv_half, s_half

      CHARACTER(len=*), PARAMETER :: routineN = 'pseudo_invert_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ii, INFO, jj, LWORK, range1_eiv, &
                                                            range2_eiv, range3_eiv, unit_nr
      LOGICAL                                            :: use_both, use_ranges_only, use_thr_only
      REAL(KIND=dp)                                      :: my_shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, WORK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: temp1, temp2, temp3, temp4
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful unit_nr
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      IF (method .EQ. 1) THEN

         IF ((PRESENT(range1) .AND. (.NOT. PRESENT(range2))) .OR. (PRESENT(range2) .AND. (.NOT. PRESENT(range1)))) THEN
            CPABORT("range1 and range2 must be provided together")
         ENDIF

         IF (PRESENT(range1) .AND. PRESENT(range1_thr)) THEN
            use_both = .TRUE.
            use_thr_only = .FALSE.
            use_ranges_only = .FALSE.
         ELSE
            use_both = .FALSE.

            IF (PRESENT(range1)) THEN
               use_ranges_only = .TRUE.
            ELSE
               use_ranges_only = .FALSE.
            ENDIF

            IF (PRESENT(range1_thr)) THEN
               use_thr_only = .TRUE.
            ELSE
               use_thr_only = .FALSE.
            ENDIF

         ENDIF

         IF ((PRESENT(s_half) .AND. (.NOT. PRESENT(s_inv_half))) .OR. (PRESENT(s_inv_half) .AND. (.NOT. PRESENT(s_half)))) THEN
            CPABORT("Domain overlap matrix missing")
         ENDIF
      ENDIF

      my_shift = 0.0_dp
      IF (PRESENT(shift)) THEN
         my_shift = shift
      ENDIF

      Ainv = A
      INFO = 0

      SELECT CASE (method)
      CASE (0) ! Inversion via cholesky factorization

         CALL DPOTRF('L', N, Ainv, N, INFO)
         IF (INFO .NE. 0) THEN
            CPABORT("DPOTRF failed")
         END IF
         CALL DPOTRI('L', N, Ainv, N, INFO)
         IF (INFO .NE. 0) THEN
            CPABORT("DPOTRI failed")
         END IF
         ! complete the matrix
         DO ii = 1, N
            DO jj = ii + 1, N
               Ainv(ii, jj) = Ainv(jj, ii)
            ENDDO
            !WRITE(*,'(100F13.9)') Ainv(ii,:)
         ENDDO

      CASE (1)

         ! diagonalize first
         ALLOCATE (eigenvalues(N))
         ALLOCATE (temp1(N, N))
         ALLOCATE (temp4(N, N))
         IF (PRESENT(s_inv_half)) THEN
            CALL DSYMM('L', 'U', N, N, 1.0_dp, s_inv_half, N, A, N, 0.0_dp, temp1, N)
            CALL DSYMM('R', 'U', N, N, 1.0_dp, s_inv_half, N, temp1, N, 0.0_dp, Ainv, N)
         ENDIF
         ! Query the optimal workspace for dsyev
         LWORK = -1
         ALLOCATE (WORK(MAX(1, LWORK)))
         CALL DSYEV('V', 'L', N, Ainv, N, eigenvalues, WORK, LWORK, INFO)

         LWORK = INT(WORK(1))
         DEALLOCATE (WORK)
         ! Allocate the workspace and solve the eigenproblem
         ALLOCATE (WORK(MAX(1, LWORK)))
         CALL DSYEV('V', 'L', N, Ainv, N, eigenvalues, WORK, LWORK, INFO)

         IF (INFO .NE. 0) THEN
            IF (unit_nr > 0) WRITE (unit_nr, *) 'EIGENSYSTEM ERROR MESSAGE: ', INFO
            CPABORT("Eigenproblem routine failed")
         END IF
         DEALLOCATE (WORK)

         !WRITE(*,*) "EIGENVALS: "
         !WRITE(*,'(4F13.9)') eigenvalues(:)

         ! invert eigenvalues and use eigenvectors to compute pseudo Ainv
         ! project out near-zero eigenvalue modes
         ALLOCATE (temp2(N, N))
         IF (PRESENT(bad_modes_projector_down)) ALLOCATE (temp3(N, N))
         temp2(1:N, 1:N) = Ainv(1:N, 1:N)

         range1_eiv = 0
         range2_eiv = 0
         range3_eiv = 0

         IF (use_both) THEN
            DO jj = 1, N
               IF ((jj .LE. range2) .AND. (eigenvalues(jj) .LT. range1_thr)) THEN
                  temp1(jj, :) = temp2(:, jj)*0.0_dp
                  IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*1.0_dp
                  range1_eiv = range1_eiv + 1
               ELSE
                  temp1(jj, :) = temp2(:, jj)/(eigenvalues(jj) + my_shift)
                  IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*0.0_dp
                  range2_eiv = range2_eiv + 1
               ENDIF
            ENDDO
         ELSE
            IF (use_ranges_only) THEN
               DO jj = 1, N
                  IF (jj .LE. range1) THEN
                     temp1(jj, :) = temp2(:, jj)*0.0_dp
                     IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*1.0_dp
                     range1_eiv = range1_eiv + 1
                  ELSE IF (jj .LE. range2) THEN
                     temp1(jj, :) = temp2(:, jj)*1.0_dp
                     IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*1.0_dp
                     range2_eiv = range2_eiv + 1
                  ELSE
                     temp1(jj, :) = temp2(:, jj)/(eigenvalues(jj) + my_shift)
                     IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*0.0_dp
                     range3_eiv = range3_eiv + 1
                  ENDIF
               ENDDO
            ELSE IF (use_thr_only) THEN
               DO jj = 1, N
                  IF (eigenvalues(jj) .LT. range1_thr) THEN
                     temp1(jj, :) = temp2(:, jj)*0.0_dp
                     IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*1.0_dp
                     range1_eiv = range1_eiv + 1
                  ELSE
                     temp1(jj, :) = temp2(:, jj)/(eigenvalues(jj) + my_shift)
                     IF (PRESENT(bad_modes_projector_down)) temp3(jj, :) = Ainv(:, jj)*0.0_dp
                     range2_eiv = range2_eiv + 1
                  ENDIF
               ENDDO
            ELSE ! no ranges, no thresholds
               CPABORT("Invert using Cholesky. It would be faster.")
            ENDIF
         ENDIF
         !WRITE(*,*) ' EIV RANGES: ', range1_eiv, range2_eiv, range3_eiv
         IF (PRESENT(bad_modes_projector_down)) THEN
            IF (PRESENT(s_half)) THEN
               CALL DSYMM('L', 'U', N, N, 1.0_dp, s_half, N, temp2, N, 0.0_dp, Ainv, N)
               CALL DSYMM('R', 'U', N, N, 1.0_dp, s_half, N, temp3, N, 0.0_dp, temp4, N)
               CALL DGEMM('N', 'N', N, N, N, 1.0_dp, Ainv, N, temp4, N, 0.0_dp, bad_modes_projector_down, N)
            ELSE
               CALL DGEMM('N', 'N', N, N, N, 1.0_dp, temp2, N, temp3, N, 0.0_dp, bad_modes_projector_down, N)
            ENDIF
         ENDIF

         IF (PRESENT(s_inv_half)) THEN
            CALL DSYMM('L', 'U', N, N, 1.0_dp, s_inv_half, N, temp2, N, 0.0_dp, temp4, N)
            CALL DSYMM('R', 'U', N, N, 1.0_dp, s_inv_half, N, temp1, N, 0.0_dp, temp2, N)
            CALL DGEMM('N', 'N', N, N, N, 1.0_dp, temp4, N, temp2, N, 0.0_dp, Ainv, N)
         ELSE
            CALL DGEMM('N', 'N', N, N, N, 1.0_dp, temp2, N, temp1, N, 0.0_dp, Ainv, N)
         ENDIF
         DEALLOCATE (temp1, temp2, temp4)
         IF (PRESENT(bad_modes_projector_down)) DEALLOCATE (temp3)
         DEALLOCATE (eigenvalues)

      CASE DEFAULT

         CPABORT("Illegal method selected for matrix inversion")

      END SELECT

      !! compute the inversion error
      !allocate(temp1(N,N))
      !temp1=MATMUL(Ainv,A)
      !DO ii=1,N
      !   temp1(ii,ii)=temp1(ii,ii)-1.0_dp
      !ENDDO
      !temp1_error=0.0_dp
      !DO ii=1,N
      !   DO jj=1,N
      !      temp1_error=temp1_error+temp1(jj,ii)*temp1(jj,ii)
      !   ENDDO
      !ENDDO
      !WRITE(*,*) "Inversion error: ", SQRT(temp1_error)
      !deallocate(temp1)

      CALL timestop(handle)

   END SUBROUTINE pseudo_invert_matrix

! **************************************************************************************************
!> \brief Find matrix power using diagonalization
!> \param A ...
!> \param Apow ...
!> \param power ...
!> \param N ...
!> \param range1 ...
!> \param range1_thr ...
!> \param shift ...
!> \par History
!>       2012.04 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************

   SUBROUTINE pseudo_matrix_power(A, Apow, power, N, range1, range1_thr, shift)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: A
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: Apow
      REAL(KIND=dp), INTENT(IN)                          :: power
      INTEGER, INTENT(IN)                                :: N
      INTEGER, INTENT(IN), OPTIONAL                      :: range1
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: range1_thr, shift

      CHARACTER(len=*), PARAMETER :: routineN = 'pseudo_matrix_power', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, INFO, jj, LWORK, range1_eiv, &
                                                            range2_eiv, unit_nr
      LOGICAL                                            :: use_both, use_ranges, use_thr
      REAL(KIND=dp)                                      :: my_shift
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, WORK
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: temp1, temp2
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      ! get a useful unit_nr
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      IF (PRESENT(range1) .AND. PRESENT(range1_thr)) THEN
         use_both = .TRUE.
      ELSE
         use_both = .FALSE.
         IF (PRESENT(range1)) THEN
            use_ranges = .TRUE.
         ELSE
            use_ranges = .FALSE.
            IF (PRESENT(range1_thr)) THEN
               use_thr = .TRUE.
            ELSE
               use_thr = .FALSE.
            ENDIF
         ENDIF
      ENDIF

      my_shift = 0.0_dp
      IF (PRESENT(shift)) THEN
         my_shift = shift
      ENDIF

      Apow = A
      INFO = 0

      ! diagonalize first
      ALLOCATE (eigenvalues(N))
      ALLOCATE (temp1(N, N))

      ! Query the optimal workspace for dsyev
      LWORK = -1
      ALLOCATE (WORK(MAX(1, LWORK)))
      CALL DSYEV('V', 'L', N, Apow, N, eigenvalues, WORK, LWORK, INFO)

      LWORK = INT(WORK(1))
      DEALLOCATE (WORK)
      ! Allocate the workspace and solve the eigenproblem
      ALLOCATE (WORK(MAX(1, LWORK)))

      CALL DSYEV('V', 'L', N, Apow, N, eigenvalues, WORK, LWORK, INFO)

      IF (INFO .NE. 0) THEN
         IF (unit_nr > 0) WRITE (unit_nr, *) 'EIGENSYSTEM ERROR MESSAGE: ', INFO
         CPABORT("Eigenproblem routine failed")
      END IF
      DEALLOCATE (WORK)

      !WRITE(*,*) "EIGENVALS: "
      !WRITE(*,'(4F13.9)') eigenvalues(:)

      ! invert eigenvalues and use eigenvectors to compute pseudo Ainv
      ! project out near-zero eigenvalue modes
      ALLOCATE (temp2(N, N))

      temp2(1:N, 1:N) = Apow(1:N, 1:N)

      range1_eiv = 0
      range2_eiv = 0

      IF (use_both) THEN
         DO jj = 1, N
            IF ((jj .LE. range1) .AND. (eigenvalues(jj) .LT. range1_thr)) THEN
               temp1(jj, :) = temp2(:, jj)*0.0_dp
               range1_eiv = range1_eiv + 1
            ELSE
               temp1(jj, :) = temp2(:, jj)*((eigenvalues(jj) + my_shift)**power)
            ENDIF
         ENDDO
      ELSE
         IF (use_ranges) THEN
            DO jj = 1, N
               IF (jj .LE. range1) THEN
                  temp1(jj, :) = temp2(:, jj)*0.0_dp
                  range1_eiv = range1_eiv + 1
               ELSE
                  temp1(jj, :) = temp2(:, jj)*((eigenvalues(jj) + my_shift)**power)
               ENDIF
            ENDDO
         ELSE
            IF (use_thr) THEN
               DO jj = 1, N
                  IF (eigenvalues(jj) .LT. range1_thr) THEN
                     temp1(jj, :) = temp2(:, jj)*0.0_dp

                     range1_eiv = range1_eiv + 1
                  ELSE
                     temp1(jj, :) = temp2(:, jj)*((eigenvalues(jj) + my_shift)**power)
                  ENDIF
               ENDDO
            ELSE
               DO jj = 1, N
                  temp1(jj, :) = temp2(:, jj)*((eigenvalues(jj) + my_shift)**power)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      !WRITE(*,*) ' EIV RANGES: ', range1_eiv, range2_eiv, range3_eiv
      Apow = MATMUL(temp2, temp1)
      DEALLOCATE (temp1, temp2)
      DEALLOCATE (eigenvalues)

      CALL timestop(handle)

   END SUBROUTINE pseudo_matrix_power

! **************************************************************************************************
!> \brief Load balancing of the submatrix computations
!> \param almo_scf_env ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE distribute_domains(almo_scf_env)

      TYPE(almo_scf_env_type), INTENT(INOUT)             :: almo_scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'distribute_domains', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, idomain, least_loaded, nao, &
                                                            ncpus, ndomains
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: index0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpu_load, domain_load
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      ndomains = almo_scf_env%ndomains
      CALL dbcsr_get_info(almo_scf_env%matrix_s(1), distribution=dist)
      CALL dbcsr_distribution_get(dist, numnodes=ncpus)

      ALLOCATE (domain_load(ndomains))
      DO idomain = 1, ndomains
         nao = almo_scf_env%nbasis_of_domain(idomain)
         domain_load(idomain) = (nao*nao*nao)*1.0_dp
      ENDDO

      ALLOCATE (index0(ndomains))

      CALL sort(domain_load, ndomains, index0)

      ALLOCATE (cpu_load(ncpus))
      cpu_load(:) = 0.0_dp

      DO idomain = 1, ndomains
         least_loaded = MINLOC(cpu_load, 1)
         cpu_load(least_loaded) = cpu_load(least_loaded) + domain_load(idomain)
         almo_scf_env%cpu_of_domain(index0(idomain)) = least_loaded - 1
      ENDDO

      DEALLOCATE (cpu_load)
      DEALLOCATE (index0)
      DEALLOCATE (domain_load)

      CALL timestop(handle)

   END SUBROUTINE distribute_domains

! **************************************************************************************************
!> \brief Tests construction and release of domain submatrices
!> \param matrix_no ...
!> \param dpattern ...
!> \param map ...
!> \param node_of_domain ...
!> \par History
!>       2013.01 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_test(matrix_no, dpattern, map, node_of_domain)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_no, dpattern
      TYPE(domain_map_type), INTENT(IN)                  :: map
      INTEGER, DIMENSION(:), INTENT(IN)                  :: node_of_domain

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_test', routineP = moduleN//':'//routineN

      INTEGER                                            :: GroupID, handle, ndomains
      TYPE(dbcsr_type)                                   :: copy1
      TYPE(domain_submatrix_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: subm_nn, subm_no

      CALL timeset(routineN, handle)

      ndomains = dbcsr_nblkcols_total(dpattern)
      CALL dbcsr_get_info(dpattern, group=GroupID)

      ALLOCATE (subm_no(ndomains), subm_nn(ndomains))
      CALL init_submatrices(subm_no)
      CALL init_submatrices(subm_nn)

      !CALL dbcsr_print(matrix_nn)
      !CALL construct_submatrices(matrix_nn,subm_nn,dpattern,map,select_row_col)
      !CALL print_submatrices(subm_nn,GroupID)

      !CALL dbcsr_print(matrix_no)
      CALL construct_submatrices(matrix_no, subm_no, dpattern, map, node_of_domain, select_row)
      CALL print_submatrices(subm_no, GroupID)

      CALL dbcsr_create(copy1, template=matrix_no)
      CALL dbcsr_copy(copy1, matrix_no)
      CALL dbcsr_print(copy1)
      CALL construct_dbcsr_from_submatrices(copy1, subm_no, dpattern)
      CALL dbcsr_print(copy1)
      CALL dbcsr_release(copy1)

      CALL release_submatrices(subm_no)
      CALL release_submatrices(subm_nn)
      DEALLOCATE (subm_no, subm_nn)

      CALL timestop(handle)

   END SUBROUTINE construct_test

! **************************************************************************************************
!> \brief create the initial guess for XALMOs
!> \param m_guess ...
!> \param m_t_in ...
!> \param m_t0 ...
!> \param m_quench_t ...
!> \param m_overlap ...
!> \param m_sigma_tmpl ...
!> \param nspins ...
!> \param xalmo_history ...
!> \param assume_t0_q0x ...
!> \param optimize_theta ...
!> \param envelope_amplitude ...
!> \param eps_filter ...
!> \param order_lanczos ...
!> \param eps_lanczos ...
!> \param max_iter_lanczos ...
!> \param nocc_of_domain ...
!> \par History
!>       2016.11 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE xalmo_initial_guess(m_guess, m_t_in, m_t0, m_quench_t, &
                                  m_overlap, m_sigma_tmpl, nspins, xalmo_history, assume_t0_q0x, &
                                  optimize_theta, envelope_amplitude, eps_filter, order_lanczos, eps_lanczos, &
                                  max_iter_lanczos, nocc_of_domain)

      TYPE(dbcsr_type), DIMENSION(:), INTENT(INOUT)      :: m_guess
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_t_in, m_t0, m_quench_t
      TYPE(dbcsr_type), INTENT(IN)                       :: m_overlap
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: m_sigma_tmpl
      INTEGER, INTENT(IN)                                :: nspins
      TYPE(almo_scf_history_type), INTENT(IN)            :: xalmo_history
      LOGICAL, INTENT(IN)                                :: assume_t0_q0x, optimize_theta
      REAL(KIND=dp), INTENT(IN)                          :: envelope_amplitude, eps_filter
      INTEGER, INTENT(IN)                                :: order_lanczos
      REAL(KIND=dp), INTENT(IN)                          :: eps_lanczos
      INTEGER, INTENT(IN)                                :: max_iter_lanczos
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: nocc_of_domain

      CHARACTER(len=*), PARAMETER :: routineN = 'xalmo_initial_guess', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iaspc, ispin, istore, naspc, &
                                                            unit_nr
      LOGICAL                                            :: aspc_guess
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: m_extrapolated, m_sigma_tmp

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      IF (optimize_theta) THEN
         CPWARN("unused option")
         ! just not to trigger unused variable
         alpha = envelope_amplitude
      END IF

      ! if extrapolation order is zero then the standard guess is used
      ! ... the number of stored history points will remain zero if extrapolation order is zero
      IF (xalmo_history%istore == 0) THEN
         aspc_guess = .FALSE.
      ELSE
         aspc_guess = .TRUE.
      ENDIF

      ! create initial guess
      IF (.NOT. aspc_guess) THEN

         DO ispin = 1, nspins

            ! zero initial guess for the delocalization amplitudes
            ! or the supplied guess for orbitals
            IF (assume_t0_q0x) THEN
               CALL dbcsr_set(m_guess(ispin), 0.0_dp)
            ELSE
               ! copy coefficients to m_guess
               CALL dbcsr_copy(m_guess(ispin), m_t_in(ispin))
            ENDIF

         ENDDO !ispins

      ELSE !aspc_guess

         CALL cite_reference(Kolafa2004)
         CALL cite_reference(Kuhne2007)

         naspc = MIN(xalmo_history%istore, xalmo_history%nstore)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(/,T2,A,/,/,T3,A,I0)") &
               "Parameters for the always stable predictor-corrector (ASPC) method:", &
               "ASPC order: ", naspc
         END IF

         DO ispin = 1, nspins

            CALL dbcsr_create(m_extrapolated, &
                              template=m_quench_t(ispin), matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(m_sigma_tmp, &
                              template=m_sigma_tmpl(ispin), matrix_type=dbcsr_type_no_symmetry)

            ! set to zero before accumulation
            CALL dbcsr_set(m_guess(ispin), 0.0_dp)

            ! extrapolation
            DO iaspc = 1, naspc

               istore = MOD(xalmo_history%istore - iaspc, xalmo_history%nstore) + 1
               alpha = (-1.0_dp)**(iaspc + 1)*REAL(iaspc, KIND=dp)* &
                       binomial(2*naspc, naspc - iaspc)/binomial(2*naspc - 2, naspc - 1)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, FMT="(T3,A2,I0,A4,F10.6)") &
                     "B(", iaspc, ") = ", alpha
               END IF

               ! m_extrapolated - initialize the correct sparsity pattern
               ! it must be kept throughout extrapolation
               CALL dbcsr_copy(m_extrapolated, m_quench_t(ispin))

               ! project t0 onto the previous DMs
               ! note that t0 is projected instead of any other matrix (e.g.
               ! t_SCF from the prev step or random t)
               ! this is done to keep orbitals phase (i.e. sign) the same as in
               ! t0. if this is not done then subtracting t0 on the next step
               ! will produce a terrible guess and extrapolation will fail
               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   xalmo_history%matrix_p_up_down(ispin, istore), &
                                   m_t0(ispin), &
                                   0.0_dp, m_extrapolated, &
                                   retain_sparsity=.TRUE.)
               ! normalize MOs
               CALL orthogonalize_mos(ket=m_extrapolated, &
                                      overlap=m_sigma_tmp, &
                                      metric=m_overlap, &
                                      retain_locality=.TRUE., &
                                      only_normalize=.FALSE., &
                                      nocc_of_domain=nocc_of_domain(:, ispin), &
                                      eps_filter=eps_filter, &
                                      order_lanczos=order_lanczos, &
                                      eps_lanczos=eps_lanczos, &
                                      max_iter_lanczos=max_iter_lanczos)

               ! now accumulate. correct sparsity is ensured
               CALL dbcsr_add(m_guess(ispin), m_extrapolated, &
                              1.0_dp, (1.0_dp*alpha)/naspc)

            ENDDO !iaspc

            CALL dbcsr_release(m_extrapolated)

            ! normalize MOs
            CALL orthogonalize_mos(ket=m_guess(ispin), &
                                   overlap=m_sigma_tmp, &
                                   metric=m_overlap, &
                                   retain_locality=.TRUE., &
                                   only_normalize=.FALSE., &
                                   nocc_of_domain=nocc_of_domain(:, ispin), &
                                   eps_filter=eps_filter, &
                                   order_lanczos=order_lanczos, &
                                   eps_lanczos=eps_lanczos, &
                                   max_iter_lanczos=max_iter_lanczos)

            CALL dbcsr_release(m_sigma_tmp)

            ! project the t0 space out from the extrapolated state
            ! this can be done outside this subroutine
            IF (assume_t0_q0x) THEN
               CALL dbcsr_add(m_guess(ispin), m_t0(ispin), &
                              1.0_dp, -1.0_dp)
            ENDIF !assume_t0_q0x

         ENDDO !ispin

      ENDIF !aspc_guess?

      CALL timestop(handle)

   END SUBROUTINE xalmo_initial_guess

END MODULE almo_scf_methods

