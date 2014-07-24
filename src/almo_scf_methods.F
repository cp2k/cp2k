!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines for ALMO SCF
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
MODULE almo_scf_methods
  USE almo_scf_aux2_methods,           ONLY: &
       construct_object01, construct_original_form_object01, &
       copy_object01_data, copy_object01_gen, ij_exists, init_object01_gen, &
       op1_object01_gen, op2_object01_gen, release_object01_gen
  USE almo_scf_aux2_types,             ONLY: object01_type,&
                                             object02_type,&
                                             select_row,&
                                             select_row_col
  USE almo_scf_types,                  ONLY: almo_objectM1_type
  USE cp_dbcsr_interface,              ONLY: &
       array_data, cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_filter, &
       cp_dbcsr_finalize, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_nblkcols_total, cp_dbcsr_nblkrows_total, cp_dbcsr_release, &
       cp_dbcsr_reserve_block2d, cp_dbcsr_scale, cp_dbcsr_transposed, &
       cp_dbcsr_type, cp_dbcsr_work_create, dbcsr_distribution_mp, &
       dbcsr_mp_numnodes, dbcsr_type_no_symmetry, dbcsr_type_symmetric
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

  PUBLIC :: almo_scf_loc_hp_blk, almo_scf_get_t_blk,&
         almo_scf_p_get_t_blk, almo_scf_ortho_blk,&
         almo_scf_p_get_t, almo_scf_get_hp_blk_and_tv_blk,&
         almo_scf_get_hp_xx_and_tv_xx,&
         get_sigma_or_s,&
         invert_blk_once,& 
         get_group_complex,&
         operations_on_sets,&
         get_group_inv,&
         get_group_sqrt,&
         get_group_distr,&
         almo_scf_loc_hp_xx,&
         get_group_rdown

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_loc_hp_xx(main_objectM1,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_loc_hp_xx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, ndomains
    REAL(KIND=dp)                            :: eps_multiply
    TYPE(cp_dbcsr_type) :: matrix_tmp1, matrix_tmp2, matrix_tmp3, &
      matrix_tmp4, matrix_tmp5, matrix_tmp6, matrix_tmp7, matrix_tmp8, &
      matrix_tmp9
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: subm_tmp1, subm_tmp2, &
                                                subm_tmp3

    CALL timeset(routineN,handle)

    eps_multiply=main_objectM1%eps_filter

    DO ispin=1,main_objectM1%nspins

       ndomains = cp_dbcsr_nblkcols_total(main_objectM1%quench_t(ispin))

       CALL construct_object01(&
               main_objectM1%hfh(ispin),&
               main_objectM1%concept_hfh_xx(:,ispin),&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row_col,&
               error)

       CALL cp_dbcsr_init(matrix_tmp1, error=error)
       CALL cp_dbcsr_create(matrix_tmp1,&
                            template=main_objectM1%enter(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, main_objectM1%hfh(ispin),&
                                     main_objectM1%enter(ispin),&
                                     0.0_dp, matrix_tmp1,&
                                     filter_eps=eps_multiply,&
                                     error=error)

       CALL cp_dbcsr_init(matrix_tmp2, error=error)
       CALL cp_dbcsr_create(matrix_tmp2,&
                            template=main_objectM1%enter(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1,&
                                     main_objectM1%rem_i(ispin),&
                                     0.0_dp, matrix_tmp2,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, main_objectM1%matrix_s(1),&
               main_objectM1%enter(ispin),&
               0.0_dp, matrix_tmp1,&
               filter_eps=eps_multiply,&
               error=error)

       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
               template=main_objectM1%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp2,&
               matrix_tmp1,&
               0.0_dp, matrix_tmp4,&
               filter_eps=eps_multiply,&
               error=error)
       
       ALLOCATE(subm_tmp1(ndomains))
       CALL init_object01_gen(subm_tmp1,error)
       CALL construct_object01(&
               matrix_tmp4,&
               subm_tmp1,&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row_col,&
               error)
       CALL op2_object01_gen(1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               -1.0_dp,subm_tmp1,'N',error)
       CALL op2_object01_gen(1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               -1.0_dp,subm_tmp1,'T',error)

       CALL cp_dbcsr_init(matrix_tmp3, error=error)
       CALL cp_dbcsr_create(matrix_tmp3,&
                            template=main_objectM1%enter(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                                     matrix_tmp4,&
                                     main_objectM1%enter(ispin),&
                                     0.0_dp, matrix_tmp3,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_release(matrix_tmp4,error=error)
       
       CALL cp_dbcsr_init(matrix_tmp6, error=error)
       CALL cp_dbcsr_create(matrix_tmp6,&
                            template=main_objectM1%enter(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              matrix_tmp3,&
                              main_objectM1%rem_i(ispin),&
                              0.0_dp, matrix_tmp6,&
                              filter_eps=eps_multiply,&
                              error=error)
       
       CALL cp_dbcsr_copy(main_objectM1%gst_xx(ispin),&
               main_objectM1%quench_t(ispin),error=error)
       CALL cp_dbcsr_copy(main_objectM1%gst_xx(ispin),&
               matrix_tmp2,keep_sparsity=.TRUE.,error=error)
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
               template=main_objectM1%enter(ispin),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_copy(matrix_tmp4,&
               main_objectM1%quench_t(ispin),error=error)
       CALL cp_dbcsr_copy(matrix_tmp4,&
               matrix_tmp6,keep_sparsity=.TRUE.,error=error)
       CALL cp_dbcsr_add(main_objectM1%gst_xx(ispin),&
               matrix_tmp4,1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_release(matrix_tmp4, error=error)
       CALL cp_dbcsr_copy(matrix_tmp3,&
               matrix_tmp2,error=error)
       CALL cp_dbcsr_add(matrix_tmp3,&
               matrix_tmp6,1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
               template=main_objectM1%matrix_s(1),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp,&
               matrix_tmp3,&
               main_objectM1%enter(ispin),&
               0.0_dp, matrix_tmp4,&
               filter_eps=eps_multiply,&
               error=error)
       CALL construct_object01(&
               matrix_tmp4,&
               main_objectM1%concept_gst(:,ispin),&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row_col,&
               error)
       CALL cp_dbcsr_release(matrix_tmp4, error=error)
       ALLOCATE(subm_tmp2(ndomains))
       CALL init_object01_gen(subm_tmp2,error)
       CALL op1_object01_gen('N','N',1.0_dp,&
               main_objectM1%concept_gst(:,ispin),&
               main_objectM1%concept_ss(:,ispin),0.0_dp,subm_tmp2,error)
       CALL op1_object01_gen('N','N',1.0_dp,&
               main_objectM1%concept_ss_inv(:,ispin),&
               subm_tmp2,0.0_dp,main_objectM1%concept_gst(:,ispin),error)
       
       CALL cp_dbcsr_init(matrix_tmp5, error=error)
       CALL cp_dbcsr_create(matrix_tmp5,&
                            template=main_objectM1%matrix_s(1),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp,&
                              matrix_tmp6,&
                              matrix_tmp1,&
                              0.0_dp, matrix_tmp5,&
                              filter_eps=eps_multiply,&
                              error=error)

       CALL construct_object01(&
               matrix_tmp5,&
               subm_tmp1,&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row_col,&
               error)
       CALL cp_dbcsr_release(matrix_tmp5,error=error)
       CALL op2_object01_gen(1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               1.0_dp,subm_tmp1,'N',error)

       ALLOCATE(subm_tmp3(ndomains))
       CALL init_object01_gen(subm_tmp3,error)
       CALL construct_object01(&
               matrix_tmp2,&
               subm_tmp2,&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row,&
               error)
       CALL construct_object01(&
               matrix_tmp6,&
               subm_tmp3,&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row,&
               error)
       CALL cp_dbcsr_release(matrix_tmp6,error=error)
       CALL op2_object01_gen(1.0_dp,subm_tmp2,&
               -1.0_dp,subm_tmp3,'N',error)
       CALL construct_object01(&
               matrix_tmp1,&
               subm_tmp3,&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row,&
               error)
       CALL op1_object01_gen('N','T',1.0_dp,subm_tmp2,&
               subm_tmp3,0.0_dp,subm_tmp1,error)
       CALL op2_object01_gen(1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               1.0_dp,subm_tmp1,'N',error)
       CALL op2_object01_gen(1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               1.0_dp,subm_tmp1,'T',error)
      
       CALL cp_dbcsr_init(matrix_tmp7, error=error)
       CALL cp_dbcsr_create(matrix_tmp7,&
                            template=main_objectM1%rem_b(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                              main_objectM1%enter(ispin),&
                              matrix_tmp2,&
                              0.0_dp, matrix_tmp7,&
                              filter_eps=eps_multiply,&
                              error=error)

       CALL cp_dbcsr_init(matrix_tmp8, error=error)
       CALL cp_dbcsr_create(matrix_tmp8,&
                            template=main_objectM1%rem_b(ispin),&
                            matrix_type=dbcsr_type_symmetric,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp8,main_objectM1%rem_b(ispin),&
                          error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              main_objectM1%rem_i(ispin),&
                              matrix_tmp7,&
                              0.0_dp, matrix_tmp8,&
                              retain_sparsity=.TRUE.,&
                              filter_eps=eps_multiply,&
                              error=error)
       CALL cp_dbcsr_release(matrix_tmp7,error=error)

       CALL cp_dbcsr_init(matrix_tmp9, error=error)
       CALL cp_dbcsr_create(matrix_tmp9,&
                            template=main_objectM1%enter(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,main_objectM1%quench_t(ispin),error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,matrix_tmp1,keep_sparsity=.TRUE.,&
                          error=error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              matrix_tmp9,&
                              matrix_tmp8,&
                              0.0_dp, matrix_tmp3,&
                              filter_eps=eps_multiply,&
                              error=error)
       CALL cp_dbcsr_release(matrix_tmp8,error=error)
       CALL cp_dbcsr_release(matrix_tmp9,error=error)

       CALL construct_object01(&
               matrix_tmp3,&
               subm_tmp2,&
               main_objectM1%quench_t(ispin),&
               main_objectM1%domain_map(ispin),&
               main_objectM1%cpu_of_domain,&
               select_row,&
               error)
       CALL op1_object01_gen('N','T',1.0_dp,subm_tmp2,&
               subm_tmp3,0.0_dp,subm_tmp1,error)
       CALL op2_object01_gen(1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               1.0_dp,subm_tmp1,'N',error)
       
       CALL release_object01_gen(subm_tmp3,error)
       CALL release_object01_gen(subm_tmp2,error)
       CALL release_object01_gen(subm_tmp1,error)
       DEALLOCATE(subm_tmp3)
       DEALLOCATE(subm_tmp2)
       DEALLOCATE(subm_tmp1)
       CALL cp_dbcsr_release(matrix_tmp3,error=error)
       CALL cp_dbcsr_release(matrix_tmp2,error=error)
       CALL cp_dbcsr_release(matrix_tmp1,error=error)

    ENDDO ! spins

    CALL timestop(handle)

  END SUBROUTINE almo_scf_loc_hp_xx

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_loc_hp_blk(main_objectM1,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_loc_hp_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: eps_multiply
    TYPE(cp_dbcsr_type) :: matrix_tmp1, matrix_tmp2, matrix_tmp3, &
      matrix_tmp4, matrix_tmp5, matrix_tmp6, matrix_tmp7, matrix_tmp8, &
      matrix_tmp9, matrix_tmp_err

    CALL timeset(routineN,handle)

    eps_multiply=main_objectM1%eps_filter

    DO ispin=1,main_objectM1%nspins

       CALL cp_dbcsr_init(matrix_tmp1, error=error)
       CALL cp_dbcsr_create(matrix_tmp1,&
                            template=main_objectM1%enter(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, main_objectM1%hfh(ispin),&
                                     main_objectM1%enter_b(ispin),&
                                     0.0_dp, matrix_tmp1,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_init(matrix_tmp2, error=error)
       CALL cp_dbcsr_create(matrix_tmp2,&
                            template=main_objectM1%enter(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1,&
                                     main_objectM1%rem_i(ispin),&
                                     0.0_dp, matrix_tmp2,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, main_objectM1%matrix_s(1),&
                                     main_objectM1%enter_b(ispin),&
                                     0.0_dp, matrix_tmp1,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       CALL cp_dbcsr_init(matrix_tmp4, error=error)
       CALL cp_dbcsr_create(matrix_tmp4,&
                            template=main_objectM1%so_b(1),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp4,main_objectM1%so_b(1),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp2,&
                                     matrix_tmp1,&
                                     0.0_dp, matrix_tmp4,&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),&
                                     matrix_tmp4,&
                                     1.0_dp,-1.0_dp,error=error)
       
       CALL cp_dbcsr_init(matrix_tmp5, error=error)
       CALL cp_dbcsr_create(matrix_tmp5,&
                            template=main_objectM1%so_b(1),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_transposed(matrix_tmp5,matrix_tmp4,error=error)
       CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),matrix_tmp5,&
                                     1.0_dp,-1.0_dp,error=error)
       
       CALL cp_dbcsr_init(matrix_tmp3, error=error)
       CALL cp_dbcsr_create(matrix_tmp3,&
                            template=main_objectM1%rem_i(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                                     main_objectM1%enter_b(ispin),&
                                     matrix_tmp2,&
                                     0.0_dp, matrix_tmp3,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       CALL cp_dbcsr_init(matrix_tmp6, error=error)
       CALL cp_dbcsr_create(matrix_tmp6,&
                            template=main_objectM1%rem_i(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              main_objectM1%rem_i(ispin),&
                              matrix_tmp3,&
                              0.0_dp, matrix_tmp6,&
                              filter_eps=eps_multiply,&
                              error=error)
       
       CALL cp_dbcsr_release(matrix_tmp3,error=error)
       CALL cp_dbcsr_init(matrix_tmp3, error=error)
       CALL cp_dbcsr_create(matrix_tmp3,&
                            template=main_objectM1%enter(ispin),&
                            error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1,&
                                     matrix_tmp6,&
                                     0.0_dp, matrix_tmp3,&
                                     filter_eps=eps_multiply,&
                                     error=error)

       CPPrecondition(main_objectM1%almo_update_algorithm.eq.almo_scf_diag,cp_failure_level,routineP,error,failure)
       CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       CALL cp_dbcsr_create(matrix_tmp_err,&
               template=main_objectM1%enter_b(ispin),&
               error=error)
       CALL cp_dbcsr_copy(matrix_tmp_err,&
               matrix_tmp2,&
               error=error)
       CALL cp_dbcsr_add(matrix_tmp_err,matrix_tmp3,&
               1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_copy(main_objectM1%gst_b(ispin),&
               main_objectM1%so_bs(1),&
               error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp_err,&
               main_objectM1%enter_b(ispin),&
               0.0_dp, main_objectM1%gst_b(ispin),&
               retain_sparsity=.TRUE.,&
               filter_eps=eps_multiply,&
               error=error)
       CALL cp_dbcsr_release(matrix_tmp_err,error=error)
       CALL cp_dbcsr_init(matrix_tmp_err,error=error)
       CALL cp_dbcsr_create(matrix_tmp_err,&
               template=main_objectM1%gst_b(ispin),&
               error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
               main_objectM1%gst_b(ispin),&
               main_objectM1%so_bs(1),&
               0.0_dp, matrix_tmp_err,&
               filter_eps=eps_multiply,&
               error=error)
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
               main_objectM1%so_bsi(1),&
               matrix_tmp_err,&
               0.0_dp, main_objectM1%gst_b(ispin),&
               filter_eps=eps_multiply,&
               error=error)
       CALL cp_dbcsr_transposed(matrix_tmp_err,&
               main_objectM1%gst_b(ispin),error=error)
       CALL cp_dbcsr_add(main_objectM1%gst_b(ispin),&
               matrix_tmp_err,&
               1.0_dp,-1.0_dp,error=error)
       CALL cp_dbcsr_release(matrix_tmp_err,error=error)

       CALL cp_dbcsr_init(matrix_tmp9, error=error)
       CALL cp_dbcsr_create(matrix_tmp9,&
                            template=main_objectM1%rem_b(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,main_objectM1%rem_b(ispin),&
                          error=error)
       CALL cp_dbcsr_copy(matrix_tmp9,matrix_tmp6,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp6,error=error)
       
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp3,&
                                     matrix_tmp1,&
                                     1.0_dp, main_objectM1%hfh_b(ispin),&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       
       CALL cp_dbcsr_init(matrix_tmp7, error=error)
       CALL cp_dbcsr_create(matrix_tmp7,&
                            template=main_objectM1%enter_b(ispin),&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp7,main_objectM1%enter_b(ispin),&
                          error=error)
       CALL cp_dbcsr_copy(matrix_tmp7,matrix_tmp3,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp3,error=error)
       CALL cp_dbcsr_init(matrix_tmp8, error=error)
       CALL cp_dbcsr_create(matrix_tmp8,&
                            template=main_objectM1%enter_b(ispin),&
                            error=error)
       CALL cp_dbcsr_copy(matrix_tmp8,main_objectM1%enter_b(ispin),&
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
       
       CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),matrix_tmp4,&
                                     1.0_dp,-1.0_dp,error=error)
       
       CALL cp_dbcsr_transposed(matrix_tmp5,matrix_tmp4,error=error)
       CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),matrix_tmp5,&
                                     1.0_dp,-1.0_dp,error=error)
       
       CALL cp_dbcsr_copy(matrix_tmp7,matrix_tmp2,keep_sparsity=.TRUE.,&
                          error=error)
       CALL cp_dbcsr_release(matrix_tmp2,error=error)
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7,&
                                     matrix_tmp8,&
                                     0.0_dp, matrix_tmp4,&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),matrix_tmp4,&
                                     1.0_dp,1.0_dp,error=error)
       
       CALL cp_dbcsr_transposed(matrix_tmp5,matrix_tmp4,error=error)
       CALL cp_dbcsr_release(matrix_tmp4,error=error)
       CALL cp_dbcsr_add(main_objectM1%hfh_b(ispin),matrix_tmp5,&
                                     1.0_dp,1.0_dp,error=error)
       CALL cp_dbcsr_release(matrix_tmp5,error=error)
                                     
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp8,&
                                     matrix_tmp9,&
                                     0.0_dp, matrix_tmp7,&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_release(matrix_tmp9,error=error)
       
       CALL cp_dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp7,&
                                     matrix_tmp8,&
                                     1.0_dp, main_objectM1%hfh_b(ispin),&
                                     retain_sparsity=.TRUE.,&
                                     filter_eps=eps_multiply,&
                                     error=error)
       CALL cp_dbcsr_release(matrix_tmp7,error=error)
       CALL cp_dbcsr_release(matrix_tmp8,error=error)

    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE almo_scf_loc_hp_blk

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_get_hp_xx_and_tv_xx(main_objectM1,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_get_hp_xx_and_tv_xx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_size, idomain, &
                                                info, ispin, lwork, ndomains
    LOGICAL                                  :: failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, work
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: data_copy
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: subm_ks_xx_orthog, subm_t, &
                                                subm_tmp

    CALL timeset(routineN,handle)

    IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular .AND. &
        main_objectM1%mat_distr_aos==almo_mat_distr_atomic) THEN
       CPErrorMessage(cp_failure_level,routineP,"a domain must be located entirely on a CPU",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    ndomains=main_objectM1%ndomains
    ALLOCATE(subm_tmp(ndomains))
    ALLOCATE(subm_ks_xx_orthog(ndomains))
    ALLOCATE(subm_t(ndomains))

    DO ispin=1,main_objectM1%nspins

       CALL init_object01_gen(subm_tmp,error)
       CALL init_object01_gen(subm_ks_xx_orthog,error)
    
       CALL op1_object01_gen('N','N',1.0_dp,main_objectM1%concept_hfh_xx(:,ispin),&
               main_objectM1%concept_ss_inv(:,ispin),0.0_dp,subm_tmp,error)
       CALL op1_object01_gen('N','N',1.0_dp,main_objectM1%concept_ss_inv(:,ispin),&
               subm_tmp,0.0_dp,subm_ks_xx_orthog,error)
       CALL release_object01_gen(subm_tmp,error)
       
       CALL init_object01_gen(subm_t,error)
   
       DO idomain = 1, ndomains
       
          IF (subm_ks_xx_orthog(idomain)%domain.gt.0) THEN
   
             iblock_size=subm_ks_xx_orthog(idomain)%nrows
      
             ALLOCATE(eigenvalues(iblock_size))
             ALLOCATE(data_copy(iblock_size,iblock_size))
             data_copy(:,:)=subm_ks_xx_orthog(idomain)%mdata(:,:)
      
             LWORK = -1
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             LWORK = INT(WORK( 1 ))
             DEALLOCATE(WORK)
      
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             IF( INFO.NE.0 ) THEN
                CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
   
             IF ( main_objectM1%concept_enter(idomain,ispin)%ncols.NE.&
                main_objectM1%nocc_of_domain(idomain,ispin) ) THEN
                CPErrorMessage(cp_failure_level,routineP,"wrong domain structure",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF
             CALL copy_object01_gen(main_objectM1%concept_enter(idomain,ispin),&
                     subm_t(idomain),.FALSE.,error)
             CALL copy_object01_data(data_copy(:,1:main_objectM1%nocc_of_domain(idomain,ispin)),&
                     subm_t(idomain),error)

             DEALLOCATE(WORK)
             DEALLOCATE(data_copy)
             DEALLOCATE(eigenvalues)
      
          ENDIF 
   
       ENDDO 

       CALL release_object01_gen(subm_ks_xx_orthog,error)

       CALL op1_object01_gen('N','N',1.0_dp,main_objectM1%concept_ss_inv(:,ispin),&
               subm_t,0.0_dp,main_objectM1%concept_enter(:,ispin),error)
       CALL release_object01_gen(subm_t,error)
       
       CALL construct_original_form_object01(&
               main_objectM1%enter(ispin),&
               main_objectM1%concept_enter(:,ispin),&
               main_objectM1%quench_t(ispin),&
               error)
       CALL cp_dbcsr_filter(main_objectM1%enter(ispin),&
               main_objectM1%eps_filter,error=error)
       
    ENDDO 

    DEALLOCATE(subm_tmp)
    DEALLOCATE(subm_ks_xx_orthog)
    DEALLOCATE(subm_t)

    CALL timestop(handle)

  END SUBROUTINE almo_scf_get_hp_xx_and_tv_xx

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_get_hp_blk_and_tv_blk(main_objectM1,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'almo_scf_get_hp_blk_and_tv_blk', &
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

    IF (main_objectM1%domain_layout_aos==almo_domain_layout_molecular .AND. &
        main_objectM1%mat_distr_aos==almo_mat_distr_atomic) THEN
       CPErrorMessage(cp_failure_level,routineP,"a domain must be located entirely on a CPU",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    DO ispin=1,main_objectM1%nspins

       CALL cp_dbcsr_init(matrix_tmp,error=error)
       CALL cp_dbcsr_init(matrix_ks_blk_orthog,error=error)
       CALL cp_dbcsr_create(matrix_tmp,template=main_objectM1%hfh_b(ispin),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
       CALL cp_dbcsr_create(matrix_ks_blk_orthog,template=main_objectM1%hfh_b(ispin),&
               matrix_type=dbcsr_type_no_symmetry,&
               error=error)
   
       CALL cp_dbcsr_multiply("N","N",1.0_dp,main_objectM1%hfh_b(ispin),&
               main_objectM1%so_bsi(1),0.0_dp,matrix_tmp,&
               filter_eps=main_objectM1%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,main_objectM1%so_bsi(1),&
               matrix_tmp,0.0_dp,matrix_ks_blk_orthog,&
               filter_eps=main_objectM1%eps_filter,error=error)
       
       CALL cp_dbcsr_release(matrix_tmp,error=error)

       CALL cp_dbcsr_init(matrix_t_blk_orthog,error=error)
       CALL cp_dbcsr_init(matrix_v_blk_orthog,error=error)
       CALL cp_dbcsr_create(matrix_t_blk_orthog,template=main_objectM1%enter_b(ispin),&
               error=error)
       CALL cp_dbcsr_create(matrix_v_blk_orthog,template=main_objectM1%v_fb(ispin),&
               error=error)
       CALL cp_dbcsr_work_create(matrix_t_blk_orthog,work_mutable=.TRUE.,&
               error=error)
       CALL cp_dbcsr_work_create(matrix_v_blk_orthog,work_mutable=.TRUE.,&
               error=error)
       
       CALL cp_dbcsr_work_create(main_objectM1%eoo(ispin),work_mutable=.TRUE.,&
               error=error)
       CALL cp_dbcsr_work_create(main_objectM1%evv(ispin),work_mutable=.TRUE.,&
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

             ALLOCATE(eigenvalues(iblock_size))
             ALLOCATE(data_copy(iblock_size,iblock_size))
             data_copy(:,:)=data_p(:,:)
   
             LWORK = -1
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             LWORK = INT(WORK( 1 ))
             DEALLOCATE(WORK)
   
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             IF( INFO.NE.0 ) THEN
                CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF

             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(matrix_t_blk_orthog,iblock_row,iblock_col,p_new_block)
             nocc_of_block=SIZE(p_new_block,2)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             CPPrecondition(nocc_of_block.gt.0,cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = data_copy(:,1:nocc_of_block)
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(matrix_v_blk_orthog,iblock_row,iblock_col,p_new_block)
             nvirt_of_block=SIZE(p_new_block,2)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             CPPrecondition(nvirt_of_block.gt.0,cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = data_copy(:,(nocc_of_block+1):(nocc_of_block+nvirt_of_block))
   
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(main_objectM1%eoo(ispin),iblock_row,iblock_col,p_new_block)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = 0.0_dp 
             DO orbital=1,nocc_of_block
                p_new_block(orbital,orbital)=eigenvalues(orbital)
             ENDDO
             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(main_objectM1%evv(ispin),iblock_row,iblock_col,p_new_block)
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
       CALL cp_dbcsr_finalize(main_objectM1%eoo(ispin),error=error)
       CALL cp_dbcsr_finalize(main_objectM1%evv(ispin),error=error)
      
       CALL cp_dbcsr_filter(matrix_t_blk_orthog,main_objectM1%eps_filter,error=error)
       CALL cp_dbcsr_filter(matrix_v_blk_orthog,main_objectM1%eps_filter,error=error)
       
       CALL cp_dbcsr_release(matrix_ks_blk_orthog, error=error)
   
       CALL cp_dbcsr_multiply("N","N",1.0_dp,main_objectM1%so_bsi(1),&
               matrix_t_blk_orthog,0.0_dp,main_objectM1%enter_b(ispin),&
               filter_eps=main_objectM1%eps_filter,error=error)
       CALL cp_dbcsr_multiply("N","N",1.0_dp,main_objectM1%so_bsi(1),&
               matrix_v_blk_orthog,0.0_dp,main_objectM1%v_fb(ispin),&
               filter_eps=main_objectM1%eps_filter,error=error)
         
       CALL cp_dbcsr_release(matrix_t_blk_orthog, error=error)
       CALL cp_dbcsr_release(matrix_v_blk_orthog, error=error)

    ENDDO 

    CALL timestop(handle)

  END SUBROUTINE almo_scf_get_hp_blk_and_tv_blk

! *****************************************************************************
!> \brief ...
!> \param matrix_in ...
!> \param matrix_out ...
!> \param error ...
!> \par History
!>       2012.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE invert_blk_once(matrix_in,matrix_out,error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_in
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_out
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_blk_once', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, iblock_size
    LOGICAL                                  :: block_needed, failure
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
   
       block_needed=.FALSE.
    
       IF (iblock_row==iblock_col) THEN
           block_needed=.TRUE.
       ENDIF
   
       IF (block_needed) THEN

          ALLOCATE(data_copy(iblock_size,iblock_size))

          CALL get_gen_inv(data_p,data_copy,iblock_size,1,&
                  range1=0,range2=0,&
                  shift=1.0E-5_dp,&
                  error=error)
   
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

  END SUBROUTINE invert_blk_once

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_get_t_blk(main_objectM1,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_get_t_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iblock_col, &
                                                iblock_row, iblock_size, &
                                                info, ispin, lwork, &
                                                nocc_of_block
    LOGICAL                                  :: block_needed, failure
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, work
    REAL(kind=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: data_copy
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: data_p, p_new_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type)                      :: matrix_t_blk_tmp

    CALL timeset(routineN,handle)

    DO ispin=1,main_objectM1%nspins
       CALL cp_dbcsr_init(matrix_t_blk_tmp, error=error)
       CALL cp_dbcsr_create(matrix_t_blk_tmp,&
                            template=main_objectM1%enter_b(ispin),&
                            error=error)
       CALL cp_dbcsr_work_create(matrix_t_blk_tmp,&
               work_mutable=.TRUE.,&
               error=error)
   
       CALL cp_dbcsr_iterator_start(iter,main_objectM1%dpp_b(ispin))
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

             ALLOCATE(eigenvalues(iblock_size))
             ALLOCATE(data_copy(iblock_size,iblock_size))
             data_copy(:,:)=data_p(:,:)
   
             LWORK = -1
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             LWORK = INT(WORK( 1 ))
             DEALLOCATE(WORK)
   
             ALLOCATE(WORK(MAX(1,LWORK)))
             CALL DSYEV('V','L',iblock_size,data_copy,iblock_size,eigenvalues,WORK,LWORK,INFO)
             IF( INFO.NE.0 ) THEN
                CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END IF

             NULLIFY (p_new_block)
             CALL cp_dbcsr_reserve_block2d(matrix_t_blk_tmp,&
                     iblock_row,iblock_col,p_new_block)
             nocc_of_block=SIZE(p_new_block,2)
             CPPostcondition(ASSOCIATED(p_new_block),cp_failure_level,routineP,error,failure)
             CPPrecondition(nocc_of_block.gt.0,cp_failure_level,routineP,error,failure)
             p_new_block(:,:) = data_copy(:,iblock_size-nocc_of_block+1:)

             DEALLOCATE(WORK)
             DEALLOCATE(data_copy)
             DEALLOCATE(eigenvalues)
   
          ENDIF
       
       ENDDO
       CALL cp_dbcsr_iterator_stop(iter)
   
       CALL cp_dbcsr_finalize(matrix_t_blk_tmp,error=error)
       CALL cp_dbcsr_filter(matrix_t_blk_tmp,&
               main_objectM1%eps_filter,error=error)
       CALL cp_dbcsr_copy(main_objectM1%enter_b(ispin),&
               matrix_t_blk_tmp,error=error)
       CALL cp_dbcsr_release(matrix_t_blk_tmp,error=error)
       
    ENDDO

    CALL timestop(handle)

  END SUBROUTINE almo_scf_get_t_blk

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE get_sigma_or_s(bra,ket,overlap,metric,retain_overlap_sparsity,&
    eps_filter,error)

    TYPE(cp_dbcsr_type), INTENT(IN)          :: bra, ket
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: overlap
    TYPE(cp_dbcsr_type), INTENT(IN)          :: metric
    LOGICAL, INTENT(IN), OPTIONAL            :: retain_overlap_sparsity
    REAL(KIND=dp)                            :: eps_filter
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_sigma_or_s', &
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

    CALL cp_dbcsr_multiply("N","N",1.0_dp,&
            metric,ket,0.0_dp,tmp,&
            filter_eps=eps_filter,error=error)

    CALL cp_dbcsr_multiply("T","N",1.0_dp,&
            bra,tmp,0.0_dp,overlap,&
            retain_sparsity=local_retain_sparsity,&
            filter_eps=eps_filter,error=error)
 
    CALL cp_dbcsr_release(tmp,error=error)

    CALL timestop(handle)

  END SUBROUTINE get_sigma_or_s 

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_ortho_blk(main_objectM1,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_ortho_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    TYPE(cp_dbcsr_type)                      :: matrix_sigma_blk_sqrt, &
                                                matrix_sigma_blk_sqrt_inv, &
                                                matrix_t_blk_tmp

    CALL timeset(routineN,handle)

    DO ispin=1,main_objectM1%nspins
   
       CALL cp_dbcsr_init(matrix_t_blk_tmp, error=error)
       CALL cp_dbcsr_create(matrix_t_blk_tmp,&
                            template=main_objectM1%enter_b(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,&
                            error=error)

       CALL cp_dbcsr_multiply("N", "N", 1.0_dp, main_objectM1%so_b(1),&
                              main_objectM1%enter_b(ispin),&
                              0.0_dp, matrix_t_blk_tmp,&
                              filter_eps=main_objectM1%eps_filter,&
                              error=error)

       CALL cp_dbcsr_multiply("T", "N", 1.0_dp,&
                              main_objectM1%enter_b(ispin),&
                              matrix_t_blk_tmp,&
                              0.0_dp, main_objectM1%rem_b(ispin),&
                              filter_eps=main_objectM1%eps_filter,&
                              retain_sparsity=.TRUE.,&
                              error=error)

       CALL cp_dbcsr_init(matrix_sigma_blk_sqrt,error=error)
       CALL cp_dbcsr_init(matrix_sigma_blk_sqrt_inv,error=error)
       CALL cp_dbcsr_create(matrix_sigma_blk_sqrt,template=main_objectM1%rem_b(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,error=error) 
       CALL cp_dbcsr_create(matrix_sigma_blk_sqrt_inv,template=main_objectM1%rem_b(ispin),&
                            matrix_type=dbcsr_type_no_symmetry,error=error) 

       CALL matrix_sqrt_Newton_Schulz(matrix_sigma_blk_sqrt,matrix_sigma_blk_sqrt_inv,&
                                      main_objectM1%rem_b(ispin),&
                                      main_objectM1%eps_filter, &
                                      3, 1.0E-4_dp, 40, error=error)
       
       CALL cp_dbcsr_multiply("N", "N", 1.0_dp,&
                              main_objectM1%enter_b(ispin),&
                              matrix_sigma_blk_sqrt_inv,&
                              0.0_dp, matrix_t_blk_tmp,&
                              filter_eps=main_objectM1%eps_filter,&
                              retain_sparsity=.TRUE.,&
                              error=error)

       CALL cp_dbcsr_copy(main_objectM1%enter_b(ispin),matrix_t_blk_tmp,&
                          keep_sparsity=.TRUE.,&
                          error=error)

       CALL cp_dbcsr_release (matrix_t_blk_tmp, error=error)
       CALL cp_dbcsr_release (matrix_sigma_blk_sqrt, error=error)
       CALL cp_dbcsr_release (matrix_sigma_blk_sqrt_inv, error=error)

    END DO
  
    CALL timestop(handle)

  END SUBROUTINE almo_scf_ortho_blk

! *****************************************************************************
!> \brief ...
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
!> \note MOs can be either orthogonal or non-orthogonal
! *****************************************************************************
  SUBROUTINE almo_scf_p_get_t(t,p,eps_filter,orthog_orbs,s,sigma,sigma_inv,&
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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_p_get_t', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, use_sigma_inv_guess
    TYPE(cp_dbcsr_type)                      :: t_tmp

    CALL timeset(routineN,handle)

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
   
       CALL cp_dbcsr_multiply("N","N",1.0_dp,s,t,0.0_dp,t_tmp,&
                              filter_eps=eps_filter,&
                              error=error)
   
       CALL cp_dbcsr_multiply("T","N",1.0_dp,t,t_tmp,0.0_dp,sigma,&
                              filter_eps=eps_filter,&
                              error=error)
   
       CALL invert_Hotelling(&
               matrix_inverse=sigma_inv,&
               matrix=sigma,&
               use_inv_as_guess=use_sigma_inv_guess,&
               threshold=eps_filter,&
               error=error)
  
       CALL cp_dbcsr_multiply("N","N",1.0_dp,t,sigma_inv,0.0_dp,t_tmp,&
                              filter_eps=eps_filter,&
                              error=error)
   
       CALL cp_dbcsr_multiply("N","T",1.0_dp,t_tmp,t,0.0_dp,p,&
                              filter_eps=eps_filter,&
                              error=error)
   
       CALL cp_dbcsr_release (t_tmp, error=error)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE almo_scf_p_get_t

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param use_sigma_inv_guess ...
!> \param error ...
!> \par History
!>       2011.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! *****************************************************************************
  SUBROUTINE almo_scf_p_get_t_blk(main_objectM1,use_sigma_inv_guess,error)

    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    LOGICAL, INTENT(IN), OPTIONAL            :: use_sigma_inv_guess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'almo_scf_p_get_t_blk', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: use_guess
    REAL(KIND=dp)                            :: spin_factor

    CALL timeset(routineN,handle)

    use_guess=.FALSE.
    IF (PRESENT(use_sigma_inv_guess)) THEN
       use_guess=use_sigma_inv_guess
    ENDIF

    DO ispin=1,main_objectM1%nspins
   
       CALL almo_scf_p_get_t(t=main_objectM1%enter_b(ispin),&
                            p=main_objectM1%dpp(ispin),&
                            eps_filter=main_objectM1%eps_filter,&
                            orthog_orbs=.FALSE.,&
                            s=main_objectM1%matrix_s(1),&
                            sigma=main_objectM1%rem(ispin),&
                            sigma_inv=main_objectM1%rem_i(ispin),&
                            use_guess=use_guess,&
                            error=error)
       
       IF (main_objectM1%nspins == 1) THEN
        spin_factor = 2.0_dp
       ELSE
        spin_factor = 1.0_dp
       ENDIF
       CALL cp_dbcsr_scale(main_objectM1%dpp(ispin),spin_factor,&
                              error=error)

    END DO
  
    CALL timestop(handle)

  END SUBROUTINE almo_scf_p_get_t_blk

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE operations_on_sets(matrix_in,matrix_out,operator1,operator2,&
    dpattern,map,node_of_domain,my_action,filter_eps,matrix_trimmer,use_trimmer,&
    error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_in
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_out
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN)                             :: operator1
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: operator2
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(object02_type), INTENT(IN)          :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    INTEGER, INTENT(IN)                      :: my_action
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_dbcsr_type), INTENT(IN), &
      OPTIONAL                               :: matrix_trimmer
    LOGICAL, INTENT(IN), OPTIONAL            :: use_trimmer
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'operations_on_sets', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndomains
    LOGICAL                                  :: failure, &
                                                matrix_trimmer_required, &
                                                my_use_trimmer, &
                                                operator2_required
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: subm_in, subm_out, subm_temp

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
    CALL init_object01_gen(subm_in,error)
    CALL init_object01_gen(subm_temp,error)
    CALL init_object01_gen(subm_out,error)

    CALL construct_object01(matrix_in,subm_in,&
            dpattern,map,node_of_domain,select_row,error)
    
    IF (my_action.eq.0) THEN
       CALL op1_object01_gen('N','N',1.0_dp,operator1,&
               subm_in,0.0_dp,subm_out,error)
    ELSE IF (my_action.eq.1) THEN
       CALL copy_object01_gen(subm_in,subm_out,.TRUE.,error)
       CALL op1_object01_gen('N','N',1.0_dp,operator1,&
               subm_in,0.0_dp,subm_temp,error)
       CALL op1_object01_gen('N','N',-1.0_dp,operator2,&
               subm_temp,1.0_dp,subm_out,error)
    ELSE
       CPErrorMessage(cp_failure_level,routineP,"ILLEGAL ACTION",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL construct_original_form_object01(matrix_out,subm_out,dpattern,error)
    CALL cp_dbcsr_filter(matrix_out,filter_eps,error=error)

    CALL release_object01_gen(subm_out,error)
    CALL release_object01_gen(subm_temp,error)
    CALL release_object01_gen(subm_in,error)

    DEALLOCATE(subm_out)
    DEALLOCATE(subm_temp)
    DEALLOCATE(subm_in)

    CALL timestop(handle)
  
  END SUBROUTINE operations_on_sets

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE get_group_complex(matrix_main,subm_s_inv,&
    subm_r_down,matrix_trimmer,dpattern,map,node_of_domain,preconditioner,&
    use_trimmer,my_action,error)
    
    TYPE(cp_dbcsr_type), INTENT(INOUT)       :: matrix_main
    TYPE(object01_type), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: subm_s_inv, subm_r_down
    TYPE(cp_dbcsr_type), INTENT(INOUT), &
      OPTIONAL                               :: matrix_trimmer
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(object02_type), INTENT(IN)          :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: preconditioner
    LOGICAL, INTENT(IN), OPTIONAL            :: use_trimmer
    INTEGER, INTENT(IN)                      :: my_action
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_group_complex', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, &
                                                n_domain_mos, naos, &
                                                nblkrows_tot, ndomains, row
    INTEGER, DIMENSION(:), POINTER           :: nmos
    LOGICAL :: failure, matrix_r_required, matrix_s_inv_required, &
      matrix_trimmer_required, my_use_trimmer
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Minv
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: subm_main, subm_tmp, subm_tmp2

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
       CPErrorMessage(cp_failure_level,routineP,"S_inv IS REQUIRED",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (.NOT.PRESENT(subm_r_down).AND.matrix_r_required) THEN
       CPErrorMessage(cp_failure_level,routineP,"R IS REQUIRED",error)
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
    CALL init_object01_gen(subm_main,error)
    CALL construct_object01(matrix_main,subm_main,&
            dpattern,map,node_of_domain,select_row_col,error)
                   

    IF (my_action.eq.-1) THEN
       ALLOCATE(subm_tmp(ndomains))
       ALLOCATE(subm_tmp2(ndomains))
       CALL init_object01_gen(subm_tmp,error)
       CALL init_object01_gen(subm_tmp2,error)
       CALL op1_object01_gen('N','N',1.0_dp,subm_r_down,&
               subm_s_inv,0.0_dp,subm_tmp,error)
       CALL op1_object01_gen('N','N',1.0_dp,subm_tmp,&
               subm_main,0.0_dp,subm_tmp2,error)
       CALL op2_object01_gen(1.0_dp,subm_main,-1.0_dp,subm_tmp2,'N',error)
       CALL op2_object01_gen(1.0_dp,subm_main,-1.0_dp,subm_tmp2,'T',error)
       CALL op1_object01_gen('N','T',1.0_dp,subm_tmp2,&
               subm_tmp,1.0_dp,subm_main,error)
       CALL release_object01_gen(subm_tmp,error)
       CALL release_object01_gen(subm_tmp2,error)
       DEALLOCATE(subm_tmp2)
       DEALLOCATE(subm_tmp)
    ENDIF

    DO idomain = 1, ndomains
    
       IF (subm_main(idomain)%domain.gt.0) THEN

          n_domain_mos=0
          DO row = 1, nblkrows_tot
             IF (ij_exists(map,row,idomain,error)) THEN
               n_domain_mos=n_domain_mos+nmos(idomain)
             ENDIF
          ENDDO
       
          naos=subm_main(idomain)%nrows

          ALLOCATE(Minv(naos,naos))

              CALL get_gen_inv(A=subm_main(idomain)%mdata,Ainv=Minv,N=naos,method=1,&
                      range1=nmos(idomain),range2=n_domain_mos,error=error)
   
          CALL copy_object01_gen(subm_main(idomain),preconditioner(idomain),.FALSE.,error)
          CALL copy_object01_data(Minv,preconditioner(idomain),error)
          
          DEALLOCATE(Minv)

       ENDIF

    ENDDO 

    CALL release_object01_gen(subm_main,error)
    DEALLOCATE(subm_main)

    CALL timestop(handle)
  
  END SUBROUTINE get_group_complex

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE get_group_sqrt(matrix_s,subm_s_sqrt,subm_s_sqrt_inv,&
    dpattern,map,node_of_domain,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: subm_s_sqrt, subm_s_sqrt_inv
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(object02_type), INTENT(IN)          :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_group_sqrt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, naos, &
                                                ndomains
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Ssqrt, Ssqrtinv
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: subm_s

    CALL timeset(routineN,handle)
    
    ndomains = cp_dbcsr_nblkcols_total(dpattern)
    CPPrecondition(SIZE(subm_s_sqrt).eq.ndomains,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(subm_s_sqrt_inv).eq.ndomains,cp_failure_level,routineP,error,failure)
    ALLOCATE(subm_s(ndomains))
    CALL init_object01_gen(subm_s,error)

    CALL construct_object01(matrix_s,subm_s,&
            dpattern,map,node_of_domain,select_row_col,error)

    DO idomain = 1, ndomains
    
       IF (subm_s(idomain)%domain.gt.0) THEN

          naos=subm_s(idomain)%nrows

          ALLOCATE(Ssqrt(naos,naos))
          ALLOCATE(Ssqrtinv(naos,naos))

          CALL get_gen_sqrt(A=subm_s(idomain)%mdata,Asqrt=Ssqrt,Asqrtinv=Ssqrtinv,&
                  N=naos,error=error)
   
          CALL copy_object01_gen(subm_s(idomain),subm_s_sqrt(idomain),.FALSE.,error)
          CALL copy_object01_data(Ssqrt,subm_s_sqrt(idomain),error)
          
          CALL copy_object01_gen(subm_s(idomain),subm_s_sqrt_inv(idomain),.FALSE.,error)
          CALL copy_object01_data(Ssqrtinv,subm_s_sqrt_inv(idomain),error)
          
          DEALLOCATE(Ssqrtinv)
          DEALLOCATE(Ssqrt)

       ENDIF 

    ENDDO 

    CALL release_object01_gen(subm_s,error)
    DEALLOCATE(subm_s)

    CALL timestop(handle)
  
  END SUBROUTINE get_group_sqrt

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE get_group_inv(matrix_s,subm_s_inv,dpattern,map,&
    node_of_domain,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_s
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: subm_s_inv
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(object02_type), INTENT(IN)          :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_group_inv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, naos, &
                                                ndomains
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Sinv
    TYPE(object01_type), ALLOCATABLE, &
      DIMENSION(:)                           :: subm_s

    CALL timeset(routineN,handle)
    
    ndomains = cp_dbcsr_nblkcols_total(dpattern)

    CPPrecondition(SIZE(subm_s_inv).eq.ndomains,cp_failure_level,routineP,error,failure)
    ALLOCATE(subm_s(ndomains))
    CALL init_object01_gen(subm_s,error)

    CALL construct_object01(matrix_s,subm_s,&
            dpattern,map,node_of_domain,select_row_col,error)

    DO idomain = 1, ndomains
    
       IF (subm_s(idomain)%domain.gt.0) THEN

          naos=subm_s(idomain)%nrows

          ALLOCATE(Sinv(naos,naos))

          CALL get_gen_inv(A=subm_s(idomain)%mdata,Ainv=Sinv,N=naos,&
                  method=0,error=error)
   
          CALL copy_object01_gen(subm_s(idomain),subm_s_inv(idomain),.FALSE.,error)
          CALL copy_object01_data(Sinv,subm_s_inv(idomain),error)
          
          DEALLOCATE(Sinv)

       ENDIF 

    ENDDO 

    CALL release_object01_gen(subm_s,error)
    DEALLOCATE(subm_s)

    CALL timestop(handle)
  
  END SUBROUTINE get_group_inv

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE get_group_rdown(matrix_t,matrix_sigma_inv,matrix_s,&
    subm_r_down,dpattern,map,node_of_domain,filter_eps,error)
    
    TYPE(cp_dbcsr_type), INTENT(IN)          :: matrix_t, matrix_sigma_inv, &
                                                matrix_s
    TYPE(object01_type), DIMENSION(:), &
      INTENT(INOUT)                          :: subm_r_down
    TYPE(cp_dbcsr_type), INTENT(IN)          :: dpattern
    TYPE(object02_type), INTENT(IN)          :: map
    INTEGER, DIMENSION(:), INTENT(IN)        :: node_of_domain
    REAL(KIND=dp)                            :: filter_eps
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_group_rdown', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndomains
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: m_tmp_no_1, m_tmp_no_2, &
                                                matrix_r

    CALL timeset(routineN,handle)
    
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

    CALL construct_object01(matrix_r,subm_r_down,&
            dpattern,map,node_of_domain,select_row_col,error)
                   
    CALL cp_dbcsr_release(matrix_r,error=error)

    CALL timestop(handle)
  
  END SUBROUTINE get_group_rdown

! *****************************************************************************
!> \brief ...
!> \param A ...
!> \param Asqrt ...
!> \param Asqrtinv ...
!> \param N ...
!> \param error ...
!> \par History
!>       2013.03 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE get_gen_sqrt(A,Asqrt,Asqrtinv,N,error)
    
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: A
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: Asqrt, Asqrtinv
    INTEGER, INTENT(IN)                      :: N
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_gen_sqrt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, INFO, jj, LWORK
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues, WORK
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: test, testN

    CALL timeset(routineN,handle)

    Asqrtinv=A
    INFO=0

    ALLOCATE(eigenvalues(N))
    LWORK = -1
    ALLOCATE(WORK(MAX(1,LWORK)))
    CALL DSYEV('V','L',N,Asqrtinv,N,eigenvalues,WORK,LWORK,INFO)
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(MAX(1,LWORK)))
    CALL DSYEV('V','L',N,Asqrtinv,N,eigenvalues,WORK,LWORK,INFO)
    IF ( INFO.NE.0 ) THEN
       WRITE(*,*) 'DSYEV ERROR MESSAGE: ', INFO
       CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(WORK)
   
    ALLOCATE(test(N,N))
    DO jj=1, N
       test(jj,:)=Asqrtinv(:,jj)*SQRT(eigenvalues(jj))
    ENDDO
    ALLOCATE(testN(N,N))
    testN(:,:)=MATMUL(Asqrtinv,test)
    Asqrt=testN
    DO jj=1, N
       test(jj,:)=Asqrtinv(:,jj)/SQRT(eigenvalues(jj))
    ENDDO
    testN(:,:)=MATMUL(Asqrtinv,test)
    Asqrtinv=testN
    DEALLOCATE(test,testN)
    
    DEALLOCATE(eigenvalues)

    CALL timestop(handle)
  
  END SUBROUTINE get_gen_sqrt

! *****************************************************************************
!> \brief ...
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
  SUBROUTINE get_gen_inv(A,Ainv,N,method,range1,range2,range1_thr,range2_thr,&
    shift,error)
    
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: A
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: Ainv
    INTEGER, INTENT(IN)                      :: N, method
    INTEGER, INTENT(IN), OPTIONAL            :: range1, range2
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: range1_thr, range2_thr, shift
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_gen_inv', &
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
       DO ii=1,N
          DO jj=ii+1,N
             Ainv(ii,jj)=Ainv(jj,ii)
          ENDDO
       ENDDO

    CASE (1)
    
       ALLOCATE(eigenvalues(N))
       LWORK = -1
       ALLOCATE(WORK(MAX(1,LWORK)))
       CALL DSYEV('V','L',N,Ainv,N,eigenvalues,WORK,LWORK,INFO)
       LWORK = INT(WORK(1))
       DEALLOCATE(WORK)
       ALLOCATE(WORK(MAX(1,LWORK)))
       CALL DSYEV('V','L',N,Ainv,N,eigenvalues,WORK,LWORK,INFO)
       IF ( INFO.NE.0 ) THEN
          WRITE(*,*) 'DSYEV ERROR MESSAGE: ', INFO
          CPErrorMessage(cp_failure_level,routineP,"DSYEV failed",error)
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE(WORK)
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
       ALLOCATE(testN(N,N))
       testN(:,:)=MATMUL(Ainv,test)
       Ainv=testN
       DEALLOCATE(test,testN)
       DEALLOCATE(eigenvalues)

    CASE DEFAULT

       CPErrorMessage(cp_failure_level,routineP,"Illegal method selected for matrix inversion",error)
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)

    END SELECT

    CALL timestop(handle)
  
  END SUBROUTINE get_gen_inv

! *****************************************************************************
!> \brief ...
!> \param main_objectM1 ...
!> \param error ...
!> \par History
!>       2013.02 created [Rustam Z. Khaliullin]
!> \author Rustam Z. Khaliullin
! *****************************************************************************
  SUBROUTINE get_group_distr(main_objectM1,error)
    
    TYPE(almo_objectM1_type), INTENT(INOUT)  :: main_objectM1
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_group_distr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idomain, &
                                                least_loaded, nao, ncpus, &
                                                ndomains
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: index0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpu_load, domain_load

    CALL timeset(routineN,handle)

    ndomains = main_objectM1%ndomains
    ncpus = dbcsr_mp_numnodes(dbcsr_distribution_mp(&
           cp_dbcsr_distribution(main_objectM1%matrix_s(1))))

    ALLOCATE(domain_load(ndomains))
    DO idomain=1,ndomains
       nao=main_objectM1%nbasis_of_domain(idomain)
       domain_load(idomain)=(nao*nao*nao)*1.0_dp
    ENDDO

    ALLOCATE(index0(ndomains))
    
    CALL sort(domain_load,ndomains,index0)

    ALLOCATE(cpu_load(ncpus))
    cpu_load(:)=0.0_dp

    DO idomain=1,ndomains
      least_loaded=MINLOC(cpu_load,1)
      cpu_load(least_loaded)=cpu_load(least_loaded)+domain_load(idomain) 
      main_objectM1%cpu_of_domain(index0(idomain))=least_loaded-1
    ENDDO

    DEALLOCATE(cpu_load)
    DEALLOCATE(index0)
    DEALLOCATE(domain_load)

    CALL timestop(handle)
  
  END SUBROUTINE get_group_distr 

END MODULE almo_scf_methods

