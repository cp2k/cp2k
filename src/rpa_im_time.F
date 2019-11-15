!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for RPA with imaginary time
!> \par History
!>      10.2015 created [Jan Wilhelm]
! **************************************************************************************************
MODULE rpa_im_time

   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_clear, dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_get_occupation, dbcsr_get_stored_coordinates, &
        dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release_p, &
        dbcsr_reserve_all_blocks, dbcsr_reserve_blocks, dbcsr_scalar, dbcsr_scale, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_clear, dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, &
        dbcsr_t_copy_tensor_to_matrix, dbcsr_t_create, dbcsr_t_destroy, &
        dbcsr_t_distribution_destroy, dbcsr_t_distribution_new, dbcsr_t_distribution_type, &
        dbcsr_t_filter, dbcsr_t_get_info, dbcsr_t_get_mapping_info, dbcsr_t_mp_environ_pgrid, &
        dbcsr_t_need_contract, dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, dbcsr_t_type, &
        dbcsr_t_write_split_info
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_alltoall,&
                                              mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_sync,&
                                              mp_waitall
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type,&
                                              two_dim_int_array
   USE qs_3c_tensors,                   ONLY: cyclic_tensor_dist
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time'

   PUBLIC :: compute_mat_P_omega, &
             setup_mat_for_mem_cut_3c, &
             zero_mat_P_omega, &
             compute_transl_dm, &
             init_cell_index_rpa

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param mat_P_global ...
!> \param matrix_s ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param ispin ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param starts_array_mc_t ...
!> \param ends_array_mc_t ...
!> \param weights_cos_tf_t_to_w ...
!> \param tj ...
!> \param tau_tj ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param alpha ...
!> \param eps_filter_im_time ...
!> \param Eigenval ...
!> \param nmo ...
!> \param num_integ_points ...
!> \param cut_memory ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param para_env ...
!> \param stabilize_exp ...
!> \param qs_env ...
!> \param index_to_cell_3c ...
!> \param cell_to_index_3c ...
!> \param has_mat_P_blocks ...
!> \param do_ri_sos_laplace_mp2 ...
! **************************************************************************************************
   SUBROUTINE compute_mat_P_omega_t(mat_P_omega, fm_scaled_dm_occ_tau, &
                                    fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                    fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                    mat_P_global, &
                                    matrix_s, mao_coeff_occ, &
                                    mao_coeff_virt, ispin, &
                                    t_3c_M, t_3c_O, &
                                    starts_array_mc_t, ends_array_mc_t, &
                                    weights_cos_tf_t_to_w, &
                                    tj, tau_tj, e_fermi, eps_filter, &
                                    alpha, eps_filter_im_time, Eigenval, nmo, &
                                    num_integ_points, cut_memory, unit_nr, &
                                    mp2_env, para_env, &
                                    stabilize_exp, qs_env, index_to_cell_3c, cell_to_index_3c, &
                                    has_mat_P_blocks, do_ri_sos_laplace_mp2)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega
      TYPE(cp_fm_type), POINTER :: fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, fm_mo_coeff_occ, &
         fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, mao_coeff_occ, mao_coeff_virt
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_M
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t_3c_O
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc_t, ends_array_mc_t
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: weights_cos_tf_t_to_w
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tj
      INTEGER, INTENT(IN)                                :: num_integ_points, nmo
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time, alpha, eps_filter, &
                                                            e_fermi
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      INTEGER, INTENT(IN)                                :: cut_memory, unit_nr
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(IN)                          :: stabilize_exp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: cell_to_index_3c
      LOGICAL, DIMENSION(:, :, :, :, :), INTENT(INOUT)   :: has_mat_P_blocks
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_P_omega_t', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, handle5, handle6, i, i_cell, i_cell_R_1, i_cell_R_1_minus_S, &
         i_cell_R_1_minus_T, i_cell_R_2, i_cell_R_2_minus_S_minus_T, i_cell_S, i_cell_T, i_mem, &
         iquad, j, j_mem, jquad, num_3c_repl, num_cells_dm, unit_nr_prv
      INTEGER(KIND=int_8)                                :: flops_1_max_occ, flops_1_max_virt, &
                                                            flops_1_occ, flops_1_virt, flops_2, &
                                                            flops_2_max
      INTEGER, DIMENSION(2, 1)                           :: ibounds_2, jbounds_2
      INTEGER, DIMENSION(2, 2)                           :: ibounds_1, jbounds_1
      INTEGER, DIMENSION(3)                              :: bounds_3c, pdims
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL :: do_Gamma_RPA, do_kpoints_cubic_RPA, do_opt_pgrid, first_cycle_im_time, &
         first_cycle_omega_loop, memory_info, pgrid_1_init_occ, pgrid_1_init_virt, pgrid_2_init, &
         R_1_minus_S_needed, R_1_minus_T_needed, R_2_minus_S_minus_T_needed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: does_mat_P_T_tau_have_blocks
      REAL(KIND=dp)                                      :: omega, omega_old, t1, t2, tau, weight, &
                                                            weight_old
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ_global, mat_dm_virt_global
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid_1_use_occ, pgrid_1_use_virt, &
                                                            pgrid_2_use
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid_1_opt_occ, pgrid_1_opt_virt, &
                                                            pgrid_2_opt
      TYPE(dbcsr_t_type)                                 :: t_3c_M_occ, t_3c_M_occ_tmp, t_3c_M_virt, &
                                                            t_3c_M_virt_tmp, t_P
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:)      :: t_dm_occ, t_dm_virt
      TYPE(dbcsr_t_type), &
         DIMENSION(SIZE(t_3c_O, 1), SIZE(t_3c_O, 2))     :: t_3c_O_occ, t_3c_O_virt

      NULLIFY (pgrid_1_opt_occ, pgrid_1_opt_virt, pgrid_2_opt)

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints
      do_Gamma_RPA = .NOT. do_kpoints_cubic_RPA
      num_3c_repl = MAXVAL(cell_to_index_3c)
      do_opt_pgrid = qs_env%mp2_env%ri_rpa_im_time%group_size_internal

      CALL timeset(routineN, handle)

      first_cycle_im_time = .TRUE.
      DO i = 1, SIZE(t_3c_O, 1)
         DO j = 1, SIZE(t_3c_O, 2)
            CALL dbcsr_t_create(t_3c_O(i, j), t_3c_O_occ(i, j))
            CALL dbcsr_t_copy(t_3c_O(i, j), t_3c_O_occ(i, j))
            CALL dbcsr_t_create(t_3c_O(i, j), t_3c_O_virt(i, j))
            CALL dbcsr_t_copy(t_3c_O(i, j), t_3c_O_virt(i, j))
            !CALL dbcsr_t_clear(t_3c_O(i, j)) ! clearing t_3c_O is not safe because it may be used later
         ENDDO
      ENDDO

      pgrid_1_init_occ = .FALSE.; pgrid_1_init_virt = .FALSE.; pgrid_2_init = .FALSE.
      DO jquad = 1, num_integ_points

         flops_1_max_virt = 0; flops_1_max_occ = 0; flops_2_max = 0

         IF (pgrid_1_init_occ) THEN
            DO i = 1, SIZE(t_3c_O, 1)
               DO j = 1, SIZE(t_3c_O, 2)
                  CALL tensor_change_pgrid(t_3c_O_occ(i, j), pgrid_1_use_occ)
               ENDDO
            ENDDO
         ENDIF

         IF (pgrid_1_init_virt) THEN
            DO i = 1, SIZE(t_3c_O, 1)
               DO j = 1, SIZE(t_3c_O, 2)
                  CALL tensor_change_pgrid(t_3c_O_virt(i, j), pgrid_1_use_virt)
               ENDDO
            ENDDO
         ENDIF

         IF (pgrid_2_init) THEN
            CALL tensor_change_pgrid(t_3c_M, pgrid_2_use, nodata=.TRUE.)
         ENDIF

         ! loop over T for chi^T(it)
         IF (memory_info) THEN
            unit_nr_prv = unit_nr
         ELSE
            unit_nr_prv = 0
         ENDIF

         CALL mp_sync(para_env%group)
         t1 = m_walltime()

         CALL compute_mat_dm_global(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, nmo, &
                                    fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                    fm_mo_coeff_virt_scaled, mat_dm_occ_global, mat_dm_virt_global, &
                                    matrix_s, mao_coeff_occ, mao_coeff_virt, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, &
                                    unit_nr, para_env, &
                                    jquad, .FALSE., stabilize_exp, do_kpoints_cubic_RPA, qs_env, &
                                    num_cells_dm, index_to_cell_dm, &
                                    does_mat_P_T_tau_have_blocks)

         IF (memory_info) THEN
            IF (pgrid_1_init_occ) THEN
               CALL dbcsr_t_get_info(t_3c_O_occ(1, 1), pdims=pdims)
               IF (unit_nr_prv > 0) THEN
                  WRITE (unit_nr_prv, "(T2,A)") "OPTIMIZED PGRID INFO M = O x D_occ"
                  WRITE (unit_nr_prv, "(T4,A,1X,3I6)") "process grid dimensions:", pdims
               ENDIF
               CALL dbcsr_t_write_split_info(pgrid_1_use_occ, unit_nr_prv)
            ENDIF

            IF (pgrid_1_init_virt) THEN
               CALL dbcsr_t_get_info(t_3c_O_virt(1, 1), pdims=pdims)
               IF (unit_nr_prv > 0) THEN
                  WRITE (unit_nr_prv, "(T2,A)") "OPTIMIZED PGRID INFO M = O x D_virt"
                  WRITE (unit_nr_prv, "(T4,A,1X,3I6)") "process grid dimensions:", pdims
               ENDIF
               CALL dbcsr_t_write_split_info(pgrid_1_use_virt, unit_nr_prv)
            ENDIF

            IF (pgrid_2_init) THEN
               CALL dbcsr_t_get_info(t_3c_M, pdims=pdims)
               IF (unit_nr_prv > 0) THEN
                  WRITE (unit_nr_prv, "(T2,A)") "OPTIMIZED PGRID INFO P = M x M"
                  WRITE (unit_nr_prv, "(T4,A,1X,3I6)") "process grid dimensions:", pdims
               ENDIF
               IF (pgrid_2_init) CALL dbcsr_t_write_split_info(pgrid_2_use, unit_nr_prv)
            ENDIF
         ENDIF

         CALL timeset(routineN//"_dbcsr_t", handle6)
         ALLOCATE (t_dm_virt(num_cells_dm))
         ALLOCATE (t_dm_occ(num_cells_dm))
         CALL dbcsr_t_create(mat_P_global%matrix, t_P, name="P (RI | RI)")
         DO i_cell = 1, num_cells_dm
            CALL dbcsr_t_create(mat_dm_virt_global(jquad, i_cell)%matrix, t_dm_virt(i_cell), name="D virt (AO | AO)")
            CALL dbcsr_t_copy_matrix_to_tensor(mat_dm_virt_global(jquad, i_cell)%matrix, t_dm_virt(i_cell))
            CALL dbcsr_clear(mat_dm_virt_global(jquad, i_cell)%matrix)
            CALL dbcsr_t_create(mat_dm_occ_global(jquad, i_cell)%matrix, t_dm_occ(i_cell), name="D occ (AO | AO)")
            CALL dbcsr_t_copy_matrix_to_tensor(mat_dm_occ_global(jquad, i_cell)%matrix, t_dm_occ(i_cell))
            CALL dbcsr_clear(mat_dm_occ_global(jquad, i_cell)%matrix)
         ENDDO
         CALL dbcsr_t_create(t_3c_O_occ(1, 1), t_3c_M_occ_tmp, name="M (RI AO | AO)")
         CALL dbcsr_t_create(t_3c_O_virt(1, 1), t_3c_M_virt_tmp, name="M (RI AO | AO)")
         CALL dbcsr_t_create(t_3c_M, t_3c_M_occ, name="M occ (RI | AO AO)")
         CALL dbcsr_t_create(t_3c_M, t_3c_M_virt, name="M virt (RI | AO AO)")

         DO i_cell_T = 1, num_cells_dm/2 + 1

            IF (.NOT. does_mat_P_T_tau_have_blocks(i_cell_T)) CYCLE

            DO j_mem = 1, cut_memory

               CALL dbcsr_t_get_info(t_3c_O_occ(1, 1), nfull_total=bounds_3c)

               jbounds_1(:, 1) = [1, bounds_3c(1)]
               jbounds_1(:, 2) = [starts_array_mc_t(j_mem), ends_array_mc_t(j_mem)]

               jbounds_2(:, 1) = [starts_array_mc_t(j_mem), ends_array_mc_t(j_mem)]

               DO i_mem = 1, cut_memory

                  ibounds_1(:, 1) = [1, bounds_3c(1)]
                  ibounds_1(:, 2) = [starts_array_mc_t(i_mem), ends_array_mc_t(i_mem)]

                  ibounds_2(:, 1) = [starts_array_mc_t(i_mem), ends_array_mc_t(i_mem)]

                  IF (unit_nr_prv > 0) WRITE (UNIT=unit_nr_prv, FMT="(T3,A,I2,1X,I2)") &
                     "RPA_IM_TIME_INFO| Memory Cut iteration", i_mem, j_mem

                  IF (do_Gamma_RPA) THEN
                     IF (.NOT. dbcsr_t_need_contract(t_3c_O_occ(1, 1), &
                                                     t_dm_occ(1), &
                                                     contract_1=[3], notcontract_1=[1, 2], &
                                                     contract_2=[2], notcontract_2=[1], &
                                                     bounds_2=jbounds_1, bounds_3=ibounds_2)) CYCLE
                     IF (.NOT. dbcsr_t_need_contract(t_3c_O_occ(1, 1), &
                                                     t_dm_virt(1), &
                                                     contract_1=[3], notcontract_1=[1, 2], &
                                                     contract_2=[2], notcontract_2=[1], &
                                                     bounds_2=ibounds_1, bounds_3=jbounds_2)) CYCLE
                  ENDIF

                  DO i_cell_R_1 = 1, num_3c_repl

                     DO i_cell_R_2 = 1, num_3c_repl

                        IF (.NOT. has_mat_P_blocks(i_cell_T, i_mem, j_mem, i_cell_R_1, i_cell_R_2)) CYCLE

                        CALL get_diff_index_3c(i_cell_R_1, i_cell_T, i_cell_R_1_minus_T, &
                                               index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                               R_1_minus_T_needed, do_kpoints_cubic_RPA)
                        DO i_cell_S = 1, num_cells_dm
                           CALL get_diff_index_3c(i_cell_R_1, i_cell_S, i_cell_R_1_minus_S, index_to_cell_3c, &
                                                  cell_to_index_3c, index_to_cell_dm, R_1_minus_S_needed, &
                                                  do_kpoints_cubic_RPA)
                           IF (R_1_minus_S_needed) THEN
                           IF (dbcsr_t_need_contract(t_3c_O_occ(i_cell_R_1_minus_S, i_cell_R_2), &
                                                     t_dm_occ(i_cell_S), &
                                                     contract_1=[3], notcontract_1=[1, 2], &
                                                     contract_2=[2], notcontract_2=[1], &
                                                     bounds_2=jbounds_1, bounds_3=ibounds_2)) THEN

                              CALL timeset(routineN//"_calc_M_occ_t", handle5)

                              CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), &
                                                    tensor_1=t_3c_O_occ(i_cell_R_1_minus_S, i_cell_R_2), &
                                                    tensor_2=t_dm_occ(i_cell_S), &
                                                    beta=dbcsr_scalar(1.0_dp), &
                                                    tensor_3=t_3c_M_occ_tmp, &
                                                    contract_1=[3], notcontract_1=[1, 2], &
                                                    contract_2=[2], notcontract_2=[1], &
                                                    map_1=[1, 2], map_2=[3], &
                                                    bounds_2=jbounds_1, bounds_3=ibounds_2, &
                                                    pgrid_opt_1=pgrid_1_opt_occ, &
                                                    filter_eps=eps_filter, unit_nr=unit_nr_prv, &
                                                    flop=flops_1_occ)

                              CALL timestop(handle5)
                              IF (do_opt_pgrid) THEN
                                 CPASSERT(ASSOCIATED(pgrid_1_opt_occ))
                                 IF (flops_1_occ .GT. flops_1_max_occ) THEN
                                    IF (pgrid_1_init_occ) CALL dbcsr_t_pgrid_destroy(pgrid_1_use_occ)
                                    pgrid_1_use_occ = pgrid_1_opt_occ
                                    DEALLOCATE (pgrid_1_opt_occ)
                                    pgrid_1_init_occ = .TRUE.
                                    flops_1_max_occ = flops_1_occ
                                 ELSE
                                    CALL dbcsr_t_pgrid_destroy(pgrid_1_opt_occ)
                                    DEALLOCATE (pgrid_1_opt_occ)
                                 ENDIF
                              ENDIF
                           ENDIF
                           ENDIF
                        ENDDO

                        CALL timeset(routineN//"_copy_M_occ_t", handle5)

                        ! copy matrix to optimal contraction layout - copy is done manually in order
                        ! to better control memory allocations (we can release data of previous
                        ! representation)
                        CALL dbcsr_t_copy(t_3c_M_occ_tmp, t_3c_M_occ, order=[1, 3, 2], move_data=.TRUE.)
                        CALL dbcsr_t_clear(t_3c_M_occ_tmp)
                        CALL dbcsr_t_filter(t_3c_M_occ, eps_filter)
                        CALL timestop(handle5)

                        DO i_cell_S = 1, num_cells_dm
                           CALL get_diff_diff_index_3c(i_cell_R_2, i_cell_S, i_cell_T, i_cell_R_2_minus_S_minus_T, &
                                                       index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                                       R_2_minus_S_minus_T_needed, do_kpoints_cubic_RPA)

                           IF (R_1_minus_T_needed .AND. R_2_minus_S_minus_T_needed) THEN

                              IF (dbcsr_t_need_contract(t_3c_O_virt(i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T), &
                                                        t_dm_virt(i_cell_S), &
                                                        contract_1=[3], notcontract_1=[1, 2], &
                                                        contract_2=[2], notcontract_2=[1], &
                                                        bounds_2=ibounds_1, bounds_3=jbounds_2)) THEN

                                 CALL timeset(routineN//"_calc_M_virt_t", handle5)
                                 CALL dbcsr_t_contract(alpha=dbcsr_scalar(alpha/2.0_dp), &
                                                       tensor_1=t_3c_O_virt( &
                                                       i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T), &
                                                       tensor_2=t_dm_virt(i_cell_S), &
                                                       beta=dbcsr_scalar(1.0_dp), &
                                                       tensor_3=t_3c_M_virt_tmp, &
                                                       contract_1=[3], notcontract_1=[1, 2], &
                                                       contract_2=[2], notcontract_2=[1], &
                                                       map_1=[1, 2], map_2=[3], &
                                                       bounds_2=ibounds_1, bounds_3=jbounds_2, &
                                                       pgrid_opt_1=pgrid_1_opt_virt, &
                                                       filter_eps=eps_filter, unit_nr=unit_nr_prv, &
                                                       flop=flops_1_virt)
                                 CALL timestop(handle5)
                                 IF (do_opt_pgrid) THEN
                                    CPASSERT(ASSOCIATED(pgrid_1_opt_virt))
                                    IF (flops_1_virt .GT. flops_1_max_virt) THEN
                                       IF (pgrid_1_init_virt) CALL dbcsr_t_pgrid_destroy(pgrid_1_use_virt)
                                       pgrid_1_use_virt = pgrid_1_opt_virt
                                       DEALLOCATE (pgrid_1_opt_virt)
                                       pgrid_1_init_virt = .TRUE.
                                       flops_1_max_virt = flops_1_virt
                                    ELSE
                                       CALL dbcsr_t_pgrid_destroy(pgrid_1_opt_virt)
                                       DEALLOCATE (pgrid_1_opt_virt)
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDDO

                        CALL timeset(routineN//"_copy_M_virt_t", handle5)
                        CALL dbcsr_t_copy(t_3c_M_virt_tmp, t_3c_M_virt, order=[1, 2, 3], move_data=.TRUE.)
                        CALL dbcsr_t_clear(t_3c_M_virt_tmp)
                        CALL dbcsr_t_filter(t_3c_M_virt, eps_filter)
                        CALL timestop(handle5)

                        flops_2 = 0
                        IF (dbcsr_t_need_contract(t_3c_M_occ, t_3c_M_virt, &
                                                  contract_1=[2, 3], notcontract_1=[1], &
                                                  contract_2=[2, 3], notcontract_2=[1])) THEN

                           CALL timeset(routineN//"_calc_P_t", handle5)

                           CALL dbcsr_t_contract(alpha=dbcsr_scalar(1.0_dp), tensor_1=t_3c_M_occ, &
                                                 tensor_2=t_3c_M_virt, &
                                                 beta=dbcsr_scalar(0.0_dp), &
                                                 tensor_3=t_P, &
                                                 contract_1=[2, 3], notcontract_1=[1], &
                                                 contract_2=[2, 3], notcontract_2=[1], &
                                                 map_1=[1], map_2=[2], &
                                                 pgrid_opt_2=pgrid_2_opt, &
                                                 filter_eps=eps_filter_im_time/REAL(cut_memory**2, KIND=dp), &
                                                 flop=flops_2, &
                                                 move_data=.TRUE., &
                                                 unit_nr=unit_nr_prv)

                           CALL dbcsr_t_copy_tensor_to_matrix(t_P, mat_P_global%matrix)

                           CALL timestop(handle5)

                           IF (do_opt_pgrid) THEN
                              CPASSERT(ASSOCIATED(pgrid_2_opt))
                              IF (flops_2 .GT. flops_2_max) THEN
                                 IF (pgrid_2_init) CALL dbcsr_t_pgrid_destroy(pgrid_2_use)
                                 pgrid_2_use = pgrid_2_opt
                                 DEALLOCATE (pgrid_2_opt)
                                 pgrid_2_init = .TRUE.
                                 flops_2_max = flops_2
                              ELSE
                                 CALL dbcsr_t_pgrid_destroy(pgrid_2_opt)
                                 DEALLOCATE (pgrid_2_opt)
                              ENDIF
                           ENDIF
                        ENDIF

                        CALL dbcsr_t_clear(t_3c_M_occ)
                        CALL dbcsr_t_clear(t_3c_M_virt)
                        CALL dbcsr_t_clear(t_P)

                        IF (do_ri_sos_laplace_mp2) THEN
                           ! For RI-SOS-Laplace-MP2 we do not perform a cosine transform,
                           ! but we have to copy P_local to the output matrix

                           CALL dbcsr_add(mat_P_omega(jquad, i_cell_T)%matrix, mat_P_global%matrix, 1.0_dp, 1.0_dp)
                        ELSE
                           CALL timeset(routineN//"_Fourier_transform", handle5)

                           ! Fourier transform of P(it) to P(iw)
                           first_cycle_omega_loop = .TRUE.

                           tau = tau_tj(jquad)

                           DO iquad = 1, num_integ_points

                              omega = tj(iquad)
                              weight = weights_cos_tf_t_to_w(iquad, jquad)

                              IF (first_cycle_omega_loop) THEN
                                 ! no multiplication with 2.0 as in Kresses paper (Kaltak, JCTC 10, 2498 (2014), Eq. 12)
                                 ! because this factor is already absorbed in the weight w_j
                                 CALL dbcsr_scale(mat_P_global%matrix, COS(omega*tau)*weight)
                              ELSE
                                 CALL dbcsr_scale(mat_P_global%matrix, COS(omega*tau)/COS(omega_old*tau)*weight/weight_old)
                              END IF

                              CALL dbcsr_add(mat_P_omega(iquad, i_cell_T)%matrix, mat_P_global%matrix, 1.0_dp, 1.0_dp)

                              first_cycle_omega_loop = .FALSE.

                              omega_old = omega
                              weight_old = weight

                           END DO

                           CALL timestop(handle5)

                           first_cycle_im_time = .FALSE.

                           CALL check_if_mat_P_T_tau_has_blocks(does_mat_P_T_tau_have_blocks, mat_P_global, i_cell_T, &
                                                                jquad, i_mem, j_mem, i_cell_R_1, i_cell_R_2, &
                                                                para_env, has_mat_P_blocks, flops_2)

                        END IF ! do_ri_sos_laplace_mp2

                        CALL dbcsr_clear(mat_P_global%matrix)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

            CALL sync_does_mat_P_T_tau_have_blocks(does_mat_P_T_tau_have_blocks, para_env, unit_nr, &
                                                   i_cell_T, jquad, index_to_cell_dm)
         ENDDO

         CALL dbcsr_t_destroy(t_P)
         DO i_cell = 1, num_cells_dm
            CALL dbcsr_t_destroy(t_dm_virt(i_cell))
            CALL dbcsr_t_destroy(t_dm_occ(i_cell))
         ENDDO

         CALL dbcsr_t_destroy(t_3c_M_occ_tmp)
         CALL dbcsr_t_destroy(t_3c_M_virt_tmp)
         CALL dbcsr_t_destroy(t_3c_M_occ)
         CALL dbcsr_t_destroy(t_3c_M_virt)
         DEALLOCATE (t_dm_virt)
         DEALLOCATE (t_dm_occ)

         CALL timestop(handle6)

         CALL mp_sync(para_env%group)
         t2 = m_walltime()

         IF (unit_nr_prv > 0) WRITE (unit_nr_prv, '(T3,A,1X,I3,A,11X,F25.6)') &
            'RPA_IM_TIME_INFO| Time for time point', jquad, ':', t2 - t1

      END DO ! time points

      DO i = 1, SIZE(t_3c_O, 1)
         DO j = 1, SIZE(t_3c_O, 2)
            CALL dbcsr_t_destroy(t_3c_O_occ(i, j))
            CALL dbcsr_t_destroy(t_3c_O_virt(i, j))
         ENDDO
      ENDDO

      IF (pgrid_1_init_virt) CALL dbcsr_t_pgrid_destroy(pgrid_1_use_virt)
      IF (pgrid_1_init_occ) CALL dbcsr_t_pgrid_destroy(pgrid_1_use_occ)
      IF (pgrid_2_init) CALL dbcsr_t_pgrid_destroy(pgrid_2_use)

      CALL clean_up(mat_dm_occ_global, mat_dm_virt_global, does_mat_P_T_tau_have_blocks)

      CALL timestop(handle)

   END SUBROUTINE compute_mat_P_omega_t

! **************************************************************************************************
!> \brief compute the matrix Q(it) (intermediate) and Fourier transform it
!>        directly to fm_mat_P_omega(iw) (output)
!> \param mat_P_omega ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param mat_P_local ...
!> \param mat_P_global ...
!> \param mat_P_global_copy ...
!> \param mat_M_mu_Pnu_occ ...
!> \param mat_M_mu_Pnu_virt ...
!> \param matrix_s ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param ispin ...
!> \param mat_M_P_munu_occ ...
!> \param mat_M_P_munu_virt ...
!> \param mat_3c_overl_int_cut ...
!> \param mat_3c_overl_int_mao_for_occ_cut ...
!> \param mat_3c_overl_int_mao_for_virt_cut ...
!> \param do_dbcsr_t ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param starts_array_mc_t ...
!> \param ends_array_mc_t ...
!> \param mat_dm_loc_occ_cut ...
!> \param mat_dm_loc_virt_cut ...
!> \param weights_cos_tf_t_to_w ...
!> \param tj ...
!> \param tau_tj ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param alpha ...
!> \param eps_filter_im_time ...
!> \param Eigenval ...
!> \param nmo ...
!> \param n_group_col ...
!> \param group_size_P ...
!> \param num_integ_points ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param unit_nr ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param starts_array_prim_col ...
!> \param ends_array_prim_col ...
!> \param starts_array_prim_row ...
!> \param ends_array_prim_row ...
!> \param starts_array_prim_fullcol ...
!> \param ends_array_prim_fullcol ...
!> \param starts_array_prim_fullrow ...
!> \param ends_array_prim_fullrow ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param offset_combi_block ...
!> \param starts_array_cm_mao_occ ...
!> \param starts_array_cm_mao_virt ...
!> \param ends_array_cm_mao_occ ...
!> \param ends_array_cm_mao_virt ...
!> \param mepos_P_from_RI_row ...
!> \param row_from_LLL ...
!> \param cycle_due_to_sparse_dm ...
!> \param multiply_needed_occ ...
!> \param multiply_needed_virt ...
!> \param non_zero_blocks_3c ...
!> \param non_zero_blocks_3c_cut_col ...
!> \param buffer_mat_M ...
!> \param do_mao ...
!> \param stabilize_exp ...
!> \param qs_env ...
!> \param index_to_cell_3c ...
!> \param cell_to_index_3c ...
!> \param needed_cutRI_mem_R1vec_R2vec_for_kp ...
!> \param has_mat_P_blocks ...
!> \param num_3c_repl ...
!> \param do_ri_sos_laplace_mp2 whether we perform ri-laplace-sos mp2
! **************************************************************************************************
   SUBROUTINE compute_mat_P_omega(mat_P_omega, fm_scaled_dm_occ_tau, &
                                  fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                  fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                  mat_P_local, mat_P_global, mat_P_global_copy, mat_M_mu_Pnu_occ, &
                                  mat_M_mu_Pnu_virt, matrix_s, mao_coeff_occ, &
                                  mao_coeff_virt, ispin, mat_M_P_munu_occ, &
                                  mat_M_P_munu_virt, mat_3c_overl_int_cut, mat_3c_overl_int_mao_for_occ_cut, &
                                  mat_3c_overl_int_mao_for_virt_cut, do_dbcsr_t, t_3c_M, t_3c_O, &
                                  starts_array_mc_t, ends_array_mc_t, &
                                  mat_dm_loc_occ_cut, mat_dm_loc_virt_cut, weights_cos_tf_t_to_w, &
                                  tj, tau_tj, e_fermi, eps_filter, &
                                  alpha, eps_filter_im_time, Eigenval, nmo, n_group_col, &
                                  group_size_P, num_integ_points, cut_memory, cut_RI, unit_nr, &
                                  mp2_env, para_env, para_env_sub, &
                                  starts_array_prim_col, ends_array_prim_col, &
                                  starts_array_prim_row, ends_array_prim_row, &
                                  starts_array_prim_fullcol, ends_array_prim_fullcol, &
                                  starts_array_prim_fullrow, ends_array_prim_fullrow, &
                                  my_group_L_starts_im_time, my_group_L_sizes_im_time, &
                                  offset_combi_block, starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                  ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                  mepos_P_from_RI_row, row_from_LLL, &
                                  cycle_due_to_sparse_dm, multiply_needed_occ, multiply_needed_virt, &
                                  non_zero_blocks_3c, non_zero_blocks_3c_cut_col, buffer_mat_M, do_mao, &
                                  stabilize_exp, qs_env, index_to_cell_3c, cell_to_index_3c, &
                                  needed_cutRI_mem_R1vec_R2vec_for_kp, &
                                  has_mat_P_blocks, num_3c_repl, do_ri_sos_laplace_mp2)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega
      TYPE(cp_fm_type), POINTER :: fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, fm_mo_coeff_occ, &
         fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type)                                 :: mat_P_local, mat_P_global, &
                                                            mat_P_global_copy
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_M_mu_Pnu_occ, mat_M_mu_Pnu_virt, &
                                                            matrix_s, mao_coeff_occ, mao_coeff_virt
      INTEGER                                            :: ispin
      TYPE(dbcsr_p_type)                                 :: mat_M_P_munu_occ, mat_M_P_munu_virt
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut, &
         mat_3c_overl_int_mao_for_occ_cut, mat_3c_overl_int_mao_for_virt_cut
      LOGICAL, INTENT(IN)                                :: do_dbcsr_t
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_M
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t_3c_O
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc_t, ends_array_mc_t
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_dm_loc_occ_cut, mat_dm_loc_virt_cut
      REAL(KIND=dp), DIMENSION(:, :)                     :: weights_cos_tf_t_to_w
      REAL(KIND=dp), DIMENSION(:)                        :: tj
      INTEGER, INTENT(IN)                                :: num_integ_points, group_size_P, &
                                                            n_group_col, nmo
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      REAL(KIND=dp)                                      :: eps_filter_im_time, alpha, eps_filter, &
                                                            e_fermi
      REAL(KIND=dp), DIMENSION(0:num_integ_points)       :: tau_tj
      INTEGER, INTENT(IN)                                :: cut_memory, cut_RI, unit_nr
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: starts_array_prim_col, ends_array_prim_col, &
         starts_array_prim_row, ends_array_prim_row, starts_array_prim_fullcol, &
         ends_array_prim_fullcol, starts_array_prim_fullrow, ends_array_prim_fullrow
      INTEGER, DIMENSION(:)                              :: my_group_L_starts_im_time, &
                                                            my_group_L_sizes_im_time
      TYPE(two_dim_int_array), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: offset_combi_block
      INTEGER, DIMENSION(:), POINTER                     :: starts_array_cm_mao_occ, &
                                                            starts_array_cm_mao_virt, &
                                                            ends_array_cm_mao_occ, &
                                                            ends_array_cm_mao_virt
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mepos_P_from_RI_row, row_from_LLL
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: cycle_due_to_sparse_dm, &
                                                            multiply_needed_occ, &
                                                            multiply_needed_virt
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: non_zero_blocks_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)        :: non_zero_blocks_3c_cut_col
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: buffer_mat_M
      LOGICAL                                            :: do_mao
      REAL(KIND=dp)                                      :: stabilize_exp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index_3c
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: needed_cutRI_mem_R1vec_R2vec_for_kp
      LOGICAL, DIMENSION(:, :, :, :, :)                  :: has_mat_P_blocks
      INTEGER                                            :: num_3c_repl
      LOGICAL                                            :: do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_P_omega', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, handle4, handle5, i_cell, i_cell_R_1, i_cell_R_1_minus_S, &
         i_cell_R_1_minus_T, i_cell_R_2, i_cell_R_2_minus_S_minus_T, i_cell_S, i_cell_T, i_cut_RI, &
         i_mem, iquad, j_mem, jquad, my_group_L_size, num_cells_dm
      INTEGER(KIND=int_8)                                :: flop_occ, flop_virt, num_flops_mat_P
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL :: do_kpoints_cubic_RPA, first_cycle_im_time, first_cycle_omega_loop, memory_info, &
         R_1_minus_S_needed, R_1_minus_T_needed, R_2_minus_S_minus_T_needed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: does_mat_P_T_tau_have_blocks
      REAL(KIND=dp)                                      :: omega, omega_old, tau, weight, weight_old
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ_global, mat_dm_virt_global

      IF (do_dbcsr_t) THEN
         CALL compute_mat_P_omega_t(mat_P_omega, fm_scaled_dm_occ_tau, &
                                    fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                    fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                    mat_P_global, &
                                    matrix_s, mao_coeff_occ, &
                                    mao_coeff_virt, ispin, &
                                    t_3c_M, t_3c_O, &
                                    starts_array_mc_t, ends_array_mc_t, &
                                    weights_cos_tf_t_to_w, &
                                    tj, tau_tj, e_fermi, eps_filter, &
                                    alpha, eps_filter_im_time, Eigenval, nmo, &
                                    num_integ_points, cut_memory, unit_nr, &
                                    mp2_env, para_env, &
                                    stabilize_exp, qs_env, index_to_cell_3c, cell_to_index_3c, &
                                    has_mat_P_blocks, do_ri_sos_laplace_mp2)
         RETURN
      ENDIF

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      CALL timeset(routineN, handle)

      first_cycle_im_time = .TRUE.

      ! iteration over time points
      DO jquad = 1, num_integ_points

         CALL compute_mat_dm_global(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, nmo, &
                                    fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                    fm_mo_coeff_virt_scaled, mat_dm_occ_global, mat_dm_virt_global, &
                                    matrix_s, mao_coeff_occ, mao_coeff_virt, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, &
                                    unit_nr, para_env, &
                                    jquad, do_mao, &
                                    stabilize_exp, do_kpoints_cubic_RPA, qs_env, num_cells_dm, index_to_cell_dm, &
                                    does_mat_P_T_tau_have_blocks)

         CALL get_cycle_due_to_sparse_dm(cycle_due_to_sparse_dm, mat_dm_occ_global, mat_dm_virt_global, &
                                         mat_3c_overl_int_cut, num_integ_points, cut_memory, cut_RI, &
                                         starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                         ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                         my_group_L_sizes_im_time, para_env, jquad)

         CALL timeset(routineN//"_im_time_repl_subgr", handle4)

         CALL mp_sync(para_env%group)

         CALL timestop(handle4)

         ! loop over T for chi^T(it)
         DO i_cell_T = 1, num_cells_dm/2 + 1

            IF (does_mat_P_T_tau_have_blocks(i_cell_T) .EQV. .FALSE.) CYCLE

            DO j_mem = 1, cut_memory

               CALL timeset(routineN//"_calc_M", handle5)
               CALL replicate_dm_to_subgroup(para_env, para_env_sub, mat_dm_virt_global, nmo, jquad, &
                                             mat_dm_loc_virt_cut, starts_array_cm_mao_virt, ends_array_cm_mao_virt, &
                                             j_mem, cut_memory, cut_RI, non_zero_blocks_3c, &
                                             non_zero_blocks_3c_cut_col, multiply_needed_virt, do_kpoints_cubic_RPA, &
                                             cell_to_index_3c, index_to_cell_dm)
               CALL timestop(handle5)

               DO i_mem = 1, cut_memory

                  IF (cycle_due_to_sparse_dm(i_mem, j_mem, jquad) .AND. (.NOT. do_kpoints_cubic_RPA)) CYCLE

                  CALL timeset(routineN//"_calc_M", handle5)
                  CALL replicate_dm_to_subgroup(para_env, para_env_sub, mat_dm_occ_global, nmo, jquad, &
                                                mat_dm_loc_occ_cut, starts_array_cm_mao_occ, ends_array_cm_mao_occ, &
                                                i_mem, cut_memory, cut_RI, non_zero_blocks_3c, &
                                                non_zero_blocks_3c_cut_col, multiply_needed_occ, do_kpoints_cubic_RPA, &
                                                cell_to_index_3c, index_to_cell_dm)
                  CALL timestop(handle5)

                  DO i_cell_R_1 = 1, num_3c_repl

                     DO i_cell_R_2 = 1, num_3c_repl

                        IF (has_mat_P_blocks(i_cell_T, i_mem, j_mem, i_cell_R_1, i_cell_R_2) .EQV. .FALSE.) CYCLE

                        CALL timeset(routineN//"_calc_M", handle5)

                        DO i_cut_RI = 1, cut_RI

                           my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

                           CALL get_diff_index_3c(i_cell_R_1, i_cell_T, i_cell_R_1_minus_T, &
                                                  index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                                  R_1_minus_T_needed, do_kpoints_cubic_RPA)

                           DO i_cell_S = 1, num_cells_dm

                              CALL get_diff_index_3c(i_cell_R_1, i_cell_S, i_cell_R_1_minus_S, index_to_cell_3c, &
                                                     cell_to_index_3c, index_to_cell_dm, R_1_minus_S_needed, &
                                                     do_kpoints_cubic_RPA)

                              ! check whether there are common mat_dm_loc_occ_cut rows and mat_3c_overl_int_cut cols
                              IF (multiply_needed_occ(j_mem, i_cut_RI, i_cell_S) .AND. R_1_minus_S_needed .AND. &
                                  needed_cutRI_mem_R1vec_R2vec_for_kp(i_cut_RI, j_mem, i_cell_R_1_minus_S, &
                                                                      i_cell_R_2)) THEN

                                 ! D^occ*(munuP)
                                 CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                                     mat_dm_loc_occ_cut(i_cut_RI, i_mem, i_cell_S)%matrix, &
                                                     mat_3c_overl_int_mao_for_occ_cut(i_cut_RI, j_mem, i_cell_R_1_minus_S, &
                                                                                      i_cell_R_2)%matrix, &
                                                     1.0_dp, mat_M_mu_Pnu_occ(i_cut_RI)%matrix, &
                                                     ! first_row=starts_array_cm(i_mem), &
                                                     ! last_row=ends_array_cm(i_mem), &
                                                     ! first_column=(starts_array_cm(j_mem)-1)*my_group_L_size+1, &
                                                     ! last_column=ends_array_cm(j_mem)*my_group_L_size, &
                                                     filter_eps=eps_filter, &
                                                     flop=flop_occ)

                              END IF

                           END DO ! i_cell_S

                           DO i_cell_S = 1, num_cells_dm

                              CALL get_diff_diff_index_3c(i_cell_R_2, i_cell_S, i_cell_T, i_cell_R_2_minus_S_minus_T, &
                                                          index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                                          R_2_minus_S_minus_T_needed, do_kpoints_cubic_RPA)

                              ! check whether there are common mat_dm_loc_virt_cut rows and mat_3c_overl_int_cut cols
                              IF (multiply_needed_virt(i_mem, i_cut_RI, i_cell_S) .AND. &
                                  R_1_minus_T_needed .AND. R_2_minus_S_minus_T_needed .AND. &
                                  needed_cutRI_mem_R1vec_R2vec_for_kp(i_cut_RI, i_mem, i_cell_R_2_minus_S_minus_T, &
                                                                      i_cell_R_1_minus_T)) THEN

                                 ! D^virt*(munuP)
                                 CALL dbcsr_multiply("N", "N", alpha/2.0_dp, &
                                                     mat_dm_loc_virt_cut(i_cut_RI, j_mem, i_cell_S)%matrix, &
                                                     mat_3c_overl_int_mao_for_virt_cut(i_cut_RI, i_mem, &
                                                                                       i_cell_R_2_minus_S_minus_T, &
                                                                                       i_cell_R_1_minus_T)%matrix, &
                                                     1.0_dp, mat_M_mu_Pnu_virt(i_cut_RI)%matrix, &
                                                     ! first_row=starts_array_cm(j_mem), &
                                                     ! last_row=ends_array_cm(j_mem), &
                                                     ! first_column=(starts_array_cm(i_mem)-1)*my_group_L_size+1, &
                                                     ! last_column=ends_array_cm(i_mem)*my_group_L_size, &
                                                     filter_eps=eps_filter, &
                                                     flop=flop_virt)

                              END IF

                           END DO ! i_cell_S

                        END DO ! i_cut_RI

                        ! free the occupied density matrix
                        ! that if is needed since we need the density matrices until the whole loops of R1 and R2
                        ! are done for kpoints
                        IF (i_cell_R_1 == num_3c_repl .AND. i_cell_R_2 == num_3c_repl) THEN
                           DO i_cut_RI = 1, cut_RI
                              DO i_cell_S = 1, num_cells_dm
                                 CALL dbcsr_set(mat_dm_loc_occ_cut(i_cut_RI, i_mem, i_cell_S)%matrix, 0.0_dp)
                                 CALL dbcsr_filter(mat_dm_loc_occ_cut(i_cut_RI, i_mem, i_cell_S)%matrix, 1.0_dp)
                              END DO
                           END DO
                        END IF

                        IF (first_cycle_im_time .AND. memory_info) THEN
                           CALL print_occupation_3c(mat_M_mu_Pnu_occ, unit_nr, "M occ before comm", para_env, cut_RI)
                           CALL print_occupation_3c(mat_M_mu_Pnu_virt, unit_nr, "M virt before comm", para_env, cut_RI)
                        END IF

                        CALL mp_sync(para_env%group)

                        CALL timestop(handle5)

                        CALL timeset(routineN//"_fill_M_tot", handle5)

                        CALL fill_mat_M_P_munu_from_M_mu_Pnu(mat_M_P_munu_occ, mat_M_mu_Pnu_occ, para_env, &
                                                             my_group_L_starts_im_time, &
                                                             my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                                             n_group_col, cut_RI, row_from_LLL, &
                                                             starts_array_prim_row, ends_array_prim_row, &
                                                             starts_array_prim_col, ends_array_prim_col, &
                                                             offset_combi_block(i_mem, j_mem)%array, &
                                                             starts_array_prim_fullcol, &
                                                             ends_array_prim_fullcol, starts_array_prim_fullrow, &
                                                             ends_array_prim_fullrow, group_size_P, &
                                                             i_mem, j_mem, buffer_mat_M, eps_filter_im_time, &
                                                             do_occ=.TRUE., do_virt=.FALSE.)

                        CALL fill_mat_M_P_munu_from_M_mu_Pnu(mat_M_P_munu_virt, mat_M_mu_Pnu_virt, para_env, &
                                                             my_group_L_starts_im_time, &
                                                             my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                                             n_group_col, cut_RI, row_from_LLL, &
                                                             starts_array_prim_row, ends_array_prim_row, &
                                                             starts_array_prim_col, ends_array_prim_col, &
                                                             offset_combi_block(i_mem, j_mem)%array, &
                                                             starts_array_prim_fullcol, &
                                                             ends_array_prim_fullcol, starts_array_prim_fullrow, &
                                                             ends_array_prim_fullrow, group_size_P, &
                                                             i_mem, j_mem, buffer_mat_M, eps_filter_im_time, &
                                                             do_occ=.FALSE., do_virt=.TRUE.)

                        IF (first_cycle_im_time .AND. memory_info) THEN
                           CALL print_occupation_2c(mat_M_P_munu_occ%matrix, unit_nr, "M occ after comm", para_env)
                           CALL print_occupation_2c(mat_M_P_munu_virt%matrix, unit_nr, "M virt after comm", para_env)
                        END IF

                        CALL mp_sync(para_env%group)

                        CALL timestop(handle5)

                        CALL timeset(routineN//"_calc_P", handle5)

                        ! P_RT = sum_mu sigma M^occ_P_mu_sigma M^virt_R_mu_sigma
                        CALL dbcsr_multiply("N", "T", 1.0_dp, mat_M_P_munu_occ%matrix, mat_M_P_munu_virt%matrix, &
                                            0.0_dp, mat_P_local%matrix, &
                                            filter_eps= &
                                            eps_filter_im_time/REAL(cut_memory**2*(para_env%num_pe/group_size_P), KIND=dp), &
                                            flop=num_flops_mat_P)

                        IF (first_cycle_im_time .AND. memory_info) THEN
                           CALL print_occupation_2c(mat_P_local%matrix, unit_nr, "local P matrix", para_env)
                        END IF

                        ! release memory
                        CALL dbcsr_set(mat_M_P_munu_occ%matrix, 0.0_dp)
                        CALL dbcsr_filter(mat_M_P_munu_occ%matrix, 1.0_dp)

                        CALL dbcsr_set(mat_M_P_munu_virt%matrix, 0.0_dp)
                        CALL dbcsr_filter(mat_M_P_munu_virt%matrix, 1.0_dp)

                        IF (first_cycle_im_time .AND. memory_info) THEN
                           CALL print_occupation_2c(mat_P_global%matrix, unit_nr, "global P matrix", para_env)
                        END IF

                        CALL mp_sync(para_env%group)

                        CALL fill_mat_P_global_from_mat_P_local(mat_P_global, mat_P_global_copy, mat_P_local, para_env, &
                                                                eps_filter_im_time)

                        CALL timestop(handle5)

                        IF (do_ri_sos_laplace_mp2) THEN
                           ! For RI-SOS-Laplace-MP2 we do not perform a cosine transform,
                           ! but we have to copy P_local to the output matrix

                           CALL dbcsr_add(mat_P_omega(jquad, i_cell_T)%matrix, mat_P_global%matrix, 1.0_dp, 1.0_dp)
                        ELSE
                           CALL timeset(routineN//"_Fourier_transform", handle5)

                           ! Fourier transform of P(it) to P(iw)
                           first_cycle_omega_loop = .TRUE.

                           tau = tau_tj(jquad)

                           DO iquad = 1, num_integ_points

                              omega = tj(iquad)
                              weight = weights_cos_tf_t_to_w(iquad, jquad)

                              IF (first_cycle_omega_loop) THEN
                                 ! no multiplication with 2.0 as in Kresses paper (Kaltak, JCTC 10, 2498 (2014), Eq. 12)
                                 ! because this factor is already absorbed in the weight w_j
                                 CALL dbcsr_scale(mat_P_global%matrix, COS(omega*tau)*weight)
                              ELSE
                                 CALL dbcsr_scale(mat_P_global%matrix, COS(omega*tau)/COS(omega_old*tau)*weight/weight_old)
                              END IF

                              CALL dbcsr_add(mat_P_omega(iquad, i_cell_T)%matrix, mat_P_global%matrix, 1.0_dp, 1.0_dp)

                              first_cycle_omega_loop = .FALSE.

                              omega_old = omega
                              weight_old = weight

                           END DO

                           CALL timestop(handle5)

                           first_cycle_im_time = .FALSE.

                           CALL check_if_mat_P_T_tau_has_blocks(does_mat_P_T_tau_have_blocks, mat_P_global, i_cell_T, &
                                                                jquad, i_mem, j_mem, i_cell_R_1, i_cell_R_2, &
                                                                para_env, has_mat_P_blocks, num_flops_mat_P)
                        END IF ! do_ri_sos_laplace_mp2

                     END DO ! i_cell_R2

                  END DO ! i_cell_R1

               END DO ! cut memory i

               ! free the virtual density matrix
               DO i_cut_RI = 1, cut_RI
                  DO i_cell = 1, num_cells_dm
                     CALL dbcsr_set(mat_dm_loc_virt_cut(i_cut_RI, j_mem, i_cell)%matrix, 0.0_dp)
                     CALL dbcsr_filter(mat_dm_loc_virt_cut(i_cut_RI, j_mem, i_cell)%matrix, 1.0_dp)
                  END DO
               END DO

            END DO ! cut memory j

            CALL sync_does_mat_P_T_tau_have_blocks(does_mat_P_T_tau_have_blocks, para_env, unit_nr, &
                                                   i_cell_T, jquad, index_to_cell_dm)

         END DO ! neighbor cell T

      END DO ! time points

      CALL clean_up(mat_dm_occ_global, mat_dm_virt_global, does_mat_P_T_tau_have_blocks)

      CALL timestop(handle)

   END SUBROUTINE compute_mat_P_omega

! **************************************************************************************************
!> \brief ...
!> \param does_mat_P_T_tau_have_blocks ...
!> \param mat_P_global ...
!> \param i_cell_T ...
!> \param jquad ...
!> \param i_mem ...
!> \param j_mem ...
!> \param i_cell_R_1 ...
!> \param i_cell_R_2 ...
!> \param para_env ...
!> \param has_mat_P_blocks ...
!> \param num_flops_mat_P ...
! **************************************************************************************************
   SUBROUTINE check_if_mat_P_T_tau_has_blocks(does_mat_P_T_tau_have_blocks, mat_P_global, i_cell_T, &
                                              jquad, i_mem, j_mem, i_cell_R_1, i_cell_R_2, &
                                              para_env, has_mat_P_blocks, num_flops_mat_P)

      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: does_mat_P_T_tau_have_blocks
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_P_global
      INTEGER, INTENT(IN)                                :: i_cell_T, jquad, i_mem, j_mem, &
                                                            i_cell_R_1, i_cell_R_2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, DIMENSION(:, :, :, :, :), INTENT(INOUT)   :: has_mat_P_blocks
      INTEGER(KIND=int_8), INTENT(INOUT)                 :: num_flops_mat_P

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_if_mat_P_T_tau_has_blocks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nblks

      CALL timeset(routineN, handle)

      IF (jquad == 1 .AND. i_mem == 1 .AND. j_mem == 1 .AND. i_cell_R_1 == 1 .AND. i_cell_R_2 == 1) THEN
         does_mat_P_T_tau_have_blocks(i_cell_T) = .FALSE.
      END IF

      nblks = dbcsr_get_num_blocks(mat_P_global%matrix)

      IF (nblks > 0) THEN
         does_mat_P_T_tau_have_blocks(i_cell_T) = .TRUE.
      END IF

      IF (jquad == 1) THEN

         CALL mp_sum(num_flops_mat_P, para_env%group)
         IF (num_flops_mat_P == 0) has_mat_P_blocks(i_cell_T, i_mem, j_mem, i_cell_R_1, i_cell_R_2) = .FALSE.

      END IF

      CALL timestop(handle)

   END SUBROUTINE check_if_mat_P_T_tau_has_blocks

! **************************************************************************************************
!> \brief ...
!> \param does_mat_P_T_tau_have_blocks ...
!> \param para_env ...
!> \param unit_nr ...
!> \param i_cell_T ...
!> \param jquad ...
!> \param index_to_cell_dm ...
! **************************************************************************************************
   SUBROUTINE sync_does_mat_P_T_tau_have_blocks(does_mat_P_T_tau_have_blocks, para_env, unit_nr, i_cell_T, &
                                                jquad, index_to_cell_dm)

      LOGICAL, DIMENSION(:), INTENT(INOUT)               :: does_mat_P_T_tau_have_blocks
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr, i_cell_T, jquad
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sync_does_mat_P_T_tau_have_blocks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, j_cell_T, j_cell_T_new
      INTEGER, ALLOCATABLE, DIMENSION(:) :: integ_does_mat_P_T_tau_have_blocks

      CALL timeset(routineN, handle)

      ALLOCATE (integ_does_mat_P_T_tau_have_blocks(SIZE(does_mat_P_T_tau_have_blocks)))
      integ_does_mat_P_T_tau_have_blocks(:) = 0

      DO j_cell_T = 1, SIZE(does_mat_P_T_tau_have_blocks)

         IF (does_mat_P_T_tau_have_blocks(j_cell_T)) THEN
            integ_does_mat_P_T_tau_have_blocks(j_cell_T) = 1
         END IF

      END DO

      CALL mp_sum(integ_does_mat_P_T_tau_have_blocks, para_env%group)

      DO j_cell_T = 1, SIZE(does_mat_P_T_tau_have_blocks)

         IF (integ_does_mat_P_T_tau_have_blocks(j_cell_T) .GE. 1) THEN
            does_mat_P_T_tau_have_blocks(j_cell_T) = .TRUE.
         ELSE IF (integ_does_mat_P_T_tau_have_blocks(j_cell_T) == 0) THEN
            does_mat_P_T_tau_have_blocks(j_cell_T) = .FALSE.
         ELSE
            CPABORT("Something is wrong when checking whether chi^T is zero.")
         END IF

      END DO

      IF ((jquad == 1) .AND. (does_mat_P_T_tau_have_blocks(i_cell_T) .EQV. .FALSE.)) THEN
         DO j_cell_T_new = i_cell_T + 1, SIZE(does_mat_P_T_tau_have_blocks)

            ! check if there is a cell which is closer to the 0-cell where the P matrix is already zero
            IF (ABS(index_to_cell_dm(1, i_cell_T)) .LE. ABS(index_to_cell_dm(1, j_cell_T_new)) .AND. &
                ABS(index_to_cell_dm(2, i_cell_T)) .LE. ABS(index_to_cell_dm(2, j_cell_T_new)) .AND. &
                ABS(index_to_cell_dm(3, i_cell_T)) .LE. ABS(index_to_cell_dm(3, j_cell_T_new))) THEN

               does_mat_P_T_tau_have_blocks(j_cell_T_new) = .FALSE.

            END IF

         END DO

      END IF

      integ_does_mat_P_T_tau_have_blocks(:) = 0
      DO j_cell_T = 1, SIZE(does_mat_P_T_tau_have_blocks)

         IF (does_mat_P_T_tau_have_blocks(j_cell_T)) THEN
            integ_does_mat_P_T_tau_have_blocks(j_cell_T) = 1
         END IF

      END DO

      IF (SUM(integ_does_mat_P_T_tau_have_blocks(i_cell_T + 1:SIZE(integ_does_mat_P_T_tau_have_blocks))) == 0) THEN
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
            "RPA_IM_TIME_INFO| Number of periodic images T for Chi_PQ^T(it):", &
            SUM(integ_does_mat_P_T_tau_have_blocks)
      END IF

      DEALLOCATE (integ_does_mat_P_T_tau_have_blocks)

      CALL timestop(handle)

   END SUBROUTINE sync_does_mat_P_T_tau_have_blocks

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param num_integ_points ...
!> \param nkp ...
! **************************************************************************************************
   SUBROUTINE zero_mat_P_omega(mat_P_omega, num_integ_points, nkp)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_P_omega
      INTEGER, INTENT(IN)                                :: num_integ_points, nkp

      INTEGER                                            :: i_kp, jquad

      DO jquad = 1, num_integ_points
         DO i_kp = 1, nkp

            CALL dbcsr_set(mat_P_omega(jquad, i_kp)%matrix, 0.0_dp)

         END DO
      END DO

   END SUBROUTINE zero_mat_P_omega

! **************************************************************************************************
!> \brief ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param mat_dm_global ...
!> \param nmo ...
!> \param jquad ...
!> \param mat_dm_loc_cut ...
!> \param starts_array_cm ...
!> \param ends_array_cm ...
!> \param i_mem ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param non_zero_blocks_3c ...
!> \param non_zero_blocks_3c_cut ...
!> \param multiply_needed ...
!> \param do_kpoints_cubic_RPA ...
!> \param cell_to_index_3c ...
!> \param index_to_cell_dm ...
! **************************************************************************************************
   SUBROUTINE replicate_dm_to_subgroup(para_env, para_env_sub, mat_dm_global, nmo, jquad, &
                                       mat_dm_loc_cut, starts_array_cm, ends_array_cm, &
                                       i_mem, cut_memory, cut_RI, non_zero_blocks_3c, &
                                       non_zero_blocks_3c_cut, multiply_needed, do_kpoints_cubic_RPA, &
                                       cell_to_index_3c, index_to_cell_dm)

      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_global
      INTEGER, INTENT(IN)                                :: nmo, jquad
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_dm_loc_cut
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_cm, ends_array_cm
      INTEGER, INTENT(IN)                                :: i_mem, cut_memory, cut_RI
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: non_zero_blocks_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :), &
         INTENT(IN)                                      :: non_zero_blocks_3c_cut
      LOGICAL, DIMENSION(:, :, :), INTENT(OUT)           :: multiply_needed
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: cell_to_index_3c
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_dm_to_subgroup', &
         routineP = moduleN//':'//routineN

      INTEGER :: block_counter, block_offset, block_size, col, col_counter, col_from_buffer, &
         col_offset, col_size, handle, handle1, i_block, i_cell, i_cell_S, i_cut_RI, i_entry, &
         i_mepos, igroup, imepos, imepos_sub, j_mem, msg_offset, nblkrows_total, ngroup, &
         num_blocks, num_cells_dm, num_pe_sub, offset, row, row_from_buffer, row_offset, row_size, &
         total_num_entries
      INTEGER, ALLOCATABLE, DIMENSION(:) :: blk_counter, cols_to_allocate_all, &
         counter_blk_to_alloc, entry_counter, num_entries_blocks_rec, num_entries_blocks_send, &
         row_block_from_index, rows_to_allocate_all, sizes_rec, sizes_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: cols_to_allocate, &
                                                            multiply_needed_tmp_int, &
                                                            rows_to_allocate
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_offset, row_blk_size
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      num_cells_dm = SIZE(mat_dm_global, 2)

      multiply_needed(:, :, :) = .FALSE.

      DO i_cell = 1, num_cells_dm

         CALL timeset("get_sizes_D", handle1)

         NULLIFY (data_block)

         ALLOCATE (num_entries_blocks_send(0:2*para_env%num_pe - 1))
         num_entries_blocks_send(:) = 0

         ALLOCATE (num_entries_blocks_rec(0:2*para_env%num_pe - 1))
         num_entries_blocks_rec(:) = 0

         ngroup = para_env%num_pe/para_env_sub%num_pe

         num_pe_sub = para_env_sub%num_pe

         CALL dbcsr_iterator_start(iter, mat_dm_global(jquad, i_cell)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            ! check whether block is in the range of the memory cutoff
            IF (row_offset + row_size - 1 >= starts_array_cm(i_mem) .AND. row_offset < ends_array_cm(i_mem)) THEN

               CALL dbcsr_get_stored_coordinates(mat_dm_loc_cut(1, 1, 1)%matrix, row, col, imepos_sub)

               DO igroup = 0, ngroup - 1

                  imepos = imepos_sub + igroup*num_pe_sub

                  num_entries_blocks_send(2*imepos) = num_entries_blocks_send(2*imepos) + row_size*col_size
                  num_entries_blocks_send(2*imepos + 1) = num_entries_blocks_send(2*imepos + 1) + 1

               END DO

            END IF

         END DO

         CALL dbcsr_iterator_stop(iter)

         CALL timestop(handle1)

         CALL timeset("send_sizes_D_1", handle1)

         total_num_entries = SUM(num_entries_blocks_send)
         CALL mp_sum(total_num_entries, para_env%group)

         CALL timestop(handle1)

         CALL timeset("send_sizes_D_2", handle1)

         IF (para_env%num_pe > 1) THEN
            CALL mp_alltoall(num_entries_blocks_send, num_entries_blocks_rec, 2, para_env%group)
         ELSE
            num_entries_blocks_rec(0:1) = num_entries_blocks_send(0:1)
         END IF

         CALL timestop(handle1)

         CALL timeset("get_data_D", handle1)

         ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
         ALLOCATE (buffer_send(0:para_env%num_pe - 1))

         ! allocate data message and corresponding indices
         DO imepos = 0, para_env%num_pe - 1

            ALLOCATE (buffer_rec(imepos)%msg(num_entries_blocks_rec(2*imepos)))
            buffer_rec(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_send(imepos)%msg(num_entries_blocks_send(2*imepos)))
            buffer_send(imepos)%msg = 0.0_dp

            ALLOCATE (buffer_rec(imepos)%indx(num_entries_blocks_rec(2*imepos + 1), 3))
            buffer_rec(imepos)%indx = 0

            ALLOCATE (buffer_send(imepos)%indx(num_entries_blocks_send(2*imepos + 1), 3))
            buffer_send(imepos)%indx = 0

         END DO

         ALLOCATE (entry_counter(0:para_env%num_pe - 1))
         entry_counter(:) = 0

         ALLOCATE (blk_counter(0:para_env%num_pe - 1))
         blk_counter = 0

         CALL dbcsr_iterator_start(iter, mat_dm_global(jquad, i_cell)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            ! check whether block is in the range of the memory cutoff
            IF (row_offset + row_size - 1 >= starts_array_cm(i_mem) .AND. row_offset < ends_array_cm(i_mem)) THEN

               CALL dbcsr_get_stored_coordinates(mat_dm_loc_cut(1, 1, 1)%matrix, row, col, imepos_sub)

               DO igroup = 0, ngroup - 1

                  imepos = imepos_sub + igroup*num_pe_sub

                  msg_offset = entry_counter(imepos)

                  block_size = row_size*col_size

                  buffer_send(imepos)%msg(msg_offset + 1:msg_offset + block_size) = &
                     RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

                  entry_counter(imepos) = entry_counter(imepos) + block_size

                  blk_counter(imepos) = blk_counter(imepos) + 1

                  block_offset = blk_counter(imepos)

                  buffer_send(imepos)%indx(block_offset, 1) = row
                  buffer_send(imepos)%indx(block_offset, 2) = col
                  buffer_send(imepos)%indx(block_offset, 3) = msg_offset

               END DO

            END IF

         END DO

         CALL dbcsr_iterator_stop(iter)

         CALL timestop(handle1)

         CALL timeset("send_data_D", handle1)

         ALLOCATE (req_array(1:para_env%num_pe, 4))

         ALLOCATE (sizes_rec(0:para_env%num_pe - 1))
         ALLOCATE (sizes_send(0:para_env%num_pe - 1))

         DO imepos = 0, para_env%num_pe - 1

            sizes_send(imepos) = num_entries_blocks_send(2*imepos)
            sizes_rec(imepos) = num_entries_blocks_rec(2*imepos)

         END DO

         CALL communicate_buffer(para_env, sizes_rec, sizes_send, buffer_rec, buffer_send, req_array)

         DEALLOCATE (req_array, sizes_rec, sizes_send)

         CALL timestop(handle1)

         CALL timeset("row_block_from_index", handle1)

         CALL dbcsr_get_info(mat_dm_loc_cut(1, 1, 1)%matrix, &
                             nblkrows_total=nblkrows_total, &
                             row_blk_offset=row_blk_offset, &
                             row_blk_size=row_blk_size)

         ALLOCATE (row_block_from_index(nmo))
         row_block_from_index = 0

         DO i_entry = 1, nmo
            DO i_block = 1, nblkrows_total

               IF (i_entry >= row_blk_offset(i_block) .AND. &
                   i_entry <= row_blk_offset(i_block) + row_blk_size(i_block) - 1) THEN

                  row_block_from_index(i_entry) = i_block

               END IF

            END DO
         END DO

         CALL timestop(handle1)

         CALL timeset("reserve_blocks_D", handle1)

         num_blocks = 0

         ! get the number of blocks, which have to be allocated
         DO imepos = 0, para_env%num_pe - 1
            num_blocks = num_blocks + num_entries_blocks_rec(2*imepos + 1)
         END DO

         ALLOCATE (rows_to_allocate_all(num_blocks))
         rows_to_allocate_all = 0

         ALLOCATE (cols_to_allocate_all(num_blocks))
         cols_to_allocate_all = 0

         block_counter = 0

         DO i_mepos = 0, para_env%num_pe - 1

            DO i_block = 1, num_entries_blocks_rec(2*i_mepos + 1)

               block_counter = block_counter + 1

               rows_to_allocate_all(block_counter) = buffer_rec(i_mepos)%indx(i_block, 1)
               cols_to_allocate_all(block_counter) = buffer_rec(i_mepos)%indx(i_block, 2)

            END DO

         END DO

         ALLOCATE (rows_to_allocate(num_blocks, cut_RI))
         rows_to_allocate = 0
         ALLOCATE (cols_to_allocate(num_blocks, cut_RI))
         cols_to_allocate = 0

         ALLOCATE (counter_blk_to_alloc(cut_RI))
         counter_blk_to_alloc = 0

         DO i_cut_RI = 1, cut_RI

            DO i_block = 1, num_blocks

               row = rows_to_allocate_all(i_block)
               col = cols_to_allocate_all(i_block)

               IF (ANY(non_zero_blocks_3c(i_cut_RI, :, :) == col)) THEN

                  counter_blk_to_alloc(i_cut_RI) = counter_blk_to_alloc(i_cut_RI) + 1
                  block_counter = counter_blk_to_alloc(i_cut_RI)

                  rows_to_allocate(block_counter, i_cut_RI) = row
                  cols_to_allocate(block_counter, i_cut_RI) = col

               END IF

            END DO

            CALL dbcsr_set(mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix, 0.0_dp)
            CALL dbcsr_filter(mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix, 1.0_dp)
            CALL dbcsr_reserve_blocks(mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix, &
                                      rows=rows_to_allocate(1:counter_blk_to_alloc(i_cut_RI), i_cut_RI), &
                                      cols=cols_to_allocate(1:counter_blk_to_alloc(i_cut_RI), i_cut_RI))
            CALL dbcsr_finalize(mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix)
            CALL dbcsr_set(mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix, 0.0_dp)

         END DO

         CALL timestop(handle1)

         CALL timeset("fill_multiply_needed", handle1)

         ! fill multiply_needed
         ALLOCATE (multiply_needed_tmp_int(cut_memory, cut_RI))
         multiply_needed_tmp_int(:, :) = 0

         ! rows from 3c matrix and cols for virt density matrix, see mat-mat-mul D^virt*(munuP)
         DO i_cut_RI = 1, cut_RI
            DO j_mem = 1, cut_memory
               DO col_counter = 1, counter_blk_to_alloc(i_cut_RI)
                  col = cols_to_allocate(col_counter, i_cut_RI)
                  IF (ANY(non_zero_blocks_3c_cut(i_cut_RI, :, :, j_mem) == col)) THEN
                     multiply_needed_tmp_int(j_mem, i_cut_RI) = 1
                  END IF
               END DO
            END DO
         END DO

         CALL mp_sum(multiply_needed_tmp_int, para_env_sub%group)

         DO i_cut_RI = 1, cut_RI
            DO j_mem = 1, cut_memory
               IF (multiply_needed_tmp_int(j_mem, i_cut_RI) > 0) THEN
                  multiply_needed(j_mem, i_cut_RI, i_cell) = .TRUE.
               END IF
            END DO
         END DO

         DEALLOCATE (rows_to_allocate_all, cols_to_allocate_all)

         CALL timestop(handle1)

         CALL timeset("fill_mat_dm_loc_cut", handle1)

         ! fill the dbcsr matrix
         DO i_cut_RI = 1, cut_RI

            CALL dbcsr_iterator_start(iter, mat_dm_loc_cut(i_cut_RI, i_mem, i_cell)%matrix)

            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                              row_size=row_size, col_size=col_size)

               DO imepos = 0, para_env%num_pe - 1

                  DO i_block = 1, num_entries_blocks_rec(2*imepos + 1)

                     row_from_buffer = buffer_rec(imepos)%indx(i_block, 1)
                     col_from_buffer = buffer_rec(imepos)%indx(i_block, 2)
                     offset = buffer_rec(imepos)%indx(i_block, 3)

                     IF (row == row_from_buffer .AND. col == col_from_buffer) THEN

                        data_block(1:row_size, 1:col_size) = &
                           RESHAPE(buffer_rec(imepos)%msg(offset + 1:offset + row_size*col_size), &
                                   (/row_size, col_size/))

                     END IF

                  END DO

               END DO

            END DO ! blocks

            CALL dbcsr_iterator_stop(iter)

         END DO ! cut_RI

         CALL timestop(handle1)

         DO imepos = 0, para_env%num_pe - 1
            DEALLOCATE (buffer_rec(imepos)%msg)
            DEALLOCATE (buffer_rec(imepos)%indx)
            DEALLOCATE (buffer_send(imepos)%msg)
            DEALLOCATE (buffer_send(imepos)%indx)
         END DO

         DEALLOCATE (buffer_rec, buffer_send)
         DEALLOCATE (entry_counter)
         DEALLOCATE (blk_counter)
         DEALLOCATE (row_block_from_index)
         DEALLOCATE (rows_to_allocate, cols_to_allocate)
         DEALLOCATE (counter_blk_to_alloc)
         DEALLOCATE (multiply_needed_tmp_int)
         DEALLOCATE (num_entries_blocks_send, num_entries_blocks_rec)

      END DO ! loop over cells

      IF (do_kpoints_cubic_RPA) THEN

         DO i_cut_RI = 1, cut_RI

            DO j_mem = 1, cut_memory

               IF (ANY(multiply_needed(j_mem, i_cut_RI, :) .EQV. .TRUE.)) THEN

                  multiply_needed(j_mem, i_cut_RI, :) = .TRUE.

               END IF

            END DO

         END DO

         DO i_cell_S = 1, num_cells_dm

            IF (index_to_cell_dm(1, i_cell_S) < LBOUND(cell_to_index_3c, 1) .OR. &
                index_to_cell_dm(1, i_cell_S) > UBOUND(cell_to_index_3c, 1) .OR. &
                index_to_cell_dm(2, i_cell_S) < LBOUND(cell_to_index_3c, 2) .OR. &
                index_to_cell_dm(2, i_cell_S) > UBOUND(cell_to_index_3c, 2) .OR. &
                index_to_cell_dm(3, i_cell_S) < LBOUND(cell_to_index_3c, 3) .OR. &
                index_to_cell_dm(3, i_cell_S) > UBOUND(cell_to_index_3c, 3)) THEN

               multiply_needed(:, :, i_cell_S) = .FALSE.

            END IF

         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE replicate_dm_to_subgroup

! **************************************************************************************************
!> \brief ...
!> \param mat_P_global ...
!> \param mat_P_global_copy ...
!> \param mat_P_local ...
!> \param para_env ...
!> \param eps_filter_im_time ...
! **************************************************************************************************
   SUBROUTINE fill_mat_P_global_from_mat_P_local(mat_P_global, mat_P_global_copy, mat_P_local, para_env, &
                                                 eps_filter_im_time)

      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global, mat_P_global_copy, &
                                                            mat_P_local
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_mat_P_global_from_mat_P_local', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: block, block_size, col, col_size, &
                                                            handle, handle1, imepos, offset, &
                                                            rec_counter, row, row_size, &
                                                            send_counter
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, entry_counter, num_blocks_rec, &
         num_blocks_send, num_entries_rec, num_entries_send, sizes_rec, sizes_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      CALL timeset("fill_Q_2_get_coord", handle1)

      CALL dbcsr_set(mat_P_global%matrix, 0.0_dp)
      CALL dbcsr_set(mat_P_global_copy%matrix, 0.0_dp)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      ALLOCATE (num_blocks_send(0:para_env%num_pe - 1))
      num_entries_rec = 0
      num_blocks_rec = 0
      num_entries_send = 0
      num_blocks_send = 0

      CALL dbcsr_iterator_start(iter, mat_P_local%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         CALL dbcsr_get_stored_coordinates(mat_P_global%matrix, row, col, imepos)

         num_entries_send(imepos) = num_entries_send(imepos) + row_size*col_size
         num_blocks_send(imepos) = num_blocks_send(imepos) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_comm_size", handle1)

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe - 1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe - 1))

         DO imepos = 0, para_env%num_pe - 1

            sizes_send(2*imepos) = num_entries_send(imepos)
            sizes_send(2*imepos + 1) = num_blocks_send(imepos)

         END DO

         CALL mp_alltoall(sizes_send, sizes_rec, 2, para_env%group)

         DO imepos = 0, para_env%num_pe - 1
            num_entries_rec(imepos) = sizes_rec(2*imepos)
            num_blocks_rec(imepos) = sizes_rec(2*imepos + 1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send(0)
         num_blocks_rec(0) = num_blocks_send(0)

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_fill_buffer", handle1)

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 5))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send(imepos), 5))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = 0

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ! fill buffer_send
      CALL dbcsr_iterator_start(iter, mat_P_local%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         CALL dbcsr_get_stored_coordinates(mat_P_global%matrix, row, col, imepos)

         block_size = row_size*col_size

         offset = entry_counter(imepos)

         buffer_send(imepos)%msg(offset + 1:offset + block_size) = &
            RESHAPE(data_block(1:row_size, 1:col_size), (/block_size/))

         block = block_counter(imepos) + 1

         buffer_send(imepos)%indx(block, 1) = row
         buffer_send(imepos)%indx(block, 2) = col
         buffer_send(imepos)%indx(block, 3) = offset

         entry_counter(imepos) = entry_counter(imepos) + block_size

         block_counter(imepos) = block_counter(imepos) + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_comm_data", handle1)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter + 1
               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
            END IF
            IF (num_entries_rec(imepos) > 0) THEN
               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter + 1
               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
            END IF
            IF (num_entries_send(imepos) > 0) THEN
               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         CALL mp_waitall(req_array(1:send_counter, 1:2))
         CALL mp_waitall(req_array(1:rec_counter, 3:4))

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_Q_2_set_blocks", handle1)

      ! fill mat_P_global_copy
      CALL dbcsr_iterator_start(iter, mat_P_global_copy%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size)

         block_size = row_size*col_size

         DO imepos = 0, para_env%num_pe - 1

            DO block = 1, num_blocks_rec(imepos)

               IF (row == buffer_rec(imepos)%indx(block, 1)) THEN

                  IF (col == buffer_rec(imepos)%indx(block, 2)) THEN

                     offset = buffer_rec(imepos)%indx(block, 3)

                     data_block(1:row_size, 1:col_size) = data_block(1:row_size, 1:col_size) + &
                                                          RESHAPE(buffer_rec(imepos)%msg(offset + 1:offset + &
                                                                                         row_size*col_size), &
                                                                  (/row_size, col_size/))

                  END IF

               END IF

            END DO

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_copy(mat_P_global%matrix, mat_P_global_copy%matrix)

      ! just remove the blocks which are exactly zero from mat_P_global
      CALL dbcsr_filter(mat_P_global%matrix, eps_filter_im_time)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      DEALLOCATE (buffer_rec, buffer_send)

      DEALLOCATE (block_counter, entry_counter)

      DEALLOCATE (req_array)

      CALL dbcsr_set(mat_P_local%matrix, 0.0_dp)
      CALL dbcsr_filter(mat_P_local%matrix, 1.0_dp)

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE fill_mat_P_global_from_mat_P_local

! **************************************************************************************************
!> \brief ...
!> \param mat_munu_array ...
!> \param unit_nr ...
!> \param matrix_name ...
!> \param para_env ...
!> \param cut_RI ...
! **************************************************************************************************
   SUBROUTINE print_occupation_3c(mat_munu_array, unit_nr, matrix_name, para_env, cut_RI)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_munu_array
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN)                       :: matrix_name
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: cut_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_occupation_3c', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cut_RI, imepos
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      REAL(KIND=dp)                                      :: local_occupation, max_occupation, &
                                                            min_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupation
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      local_occupation = 0.0_dp

      DO i_cut_RI = 1, cut_RI

         local_occupation = local_occupation + dbcsr_get_occupation(mat_munu_array(i_cut_RI)%matrix)

      END DO

      local_occupation = local_occupation/REAL(cut_RI, KIND=dp)

      ALLOCATE (buffer_send(0:para_env%num_pe - 1))
      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))

      DO imepos = 0, para_env%num_pe - 1
         ALLOCATE (buffer_send(imepos)%msg(1))
         buffer_send(imepos)%msg(1) = local_occupation

         ALLOCATE (buffer_rec(imepos)%msg(1))
      END DO

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (req_array(0:para_env%num_pe - 1, 2))

         DO imepos = 0, para_env%num_pe - 1
            CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(imepos, 1), tag=2)
            CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(imepos, 2), tag=2)
         END DO

         CALL mp_waitall(req_array(:, 1:2))

         DEALLOCATE (req_array)

      ELSE

         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      ALLOCATE (occupation(0:para_env%num_pe - 1))

      DO imepos = 0, para_env%num_pe - 1

         occupation(imepos) = buffer_rec(imepos)%msg(1)

      END DO

      max_occupation = MAXVAL(occupation)

      min_occupation = MINVAL(occupation)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%msg)
      END DO

      DEALLOCATE (buffer_send, buffer_rec, occupation)

      ! print fraction of non-zero blocks
      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T64,ES7.1,A,T74,ES7.1)") &
         "MATRIX_INFO| Frac. of non-zero bl. in ", matrix_name, ":", min_occupation, ' -', max_occupation

      CALL timestop(handle)

   END SUBROUTINE print_occupation_3c

! **************************************************************************************************
!> \brief ...
!> \param mat_munu ...
!> \param unit_nr ...
!> \param matrix_name ...
!> \param para_env ...
!> \param one_number ...
! **************************************************************************************************
   SUBROUTINE print_occupation_2c(mat_munu, unit_nr, matrix_name, para_env, one_number)
      TYPE(dbcsr_type), POINTER                          :: mat_munu
      INTEGER, INTENT(IN)                                :: unit_nr
      CHARACTER(len=*), INTENT(IN)                       :: matrix_name
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: one_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_occupation_2c', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imepos
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL                                            :: my_one_number
      REAL(KIND=dp)                                      :: local_occupation, max_occupation, &
                                                            min_occupation
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: occupation
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CALL timeset(routineN, handle)

      IF (PRESENT(one_number)) THEN
         my_one_number = one_number
      ELSE
         my_one_number = .FALSE.
      END IF

      local_occupation = dbcsr_get_occupation(mat_munu)

      ALLOCATE (buffer_send(0:para_env%num_pe - 1))
      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))

      DO imepos = 0, para_env%num_pe - 1
         ALLOCATE (buffer_send(imepos)%msg(1))
         buffer_send(imepos)%msg(1) = local_occupation

         ALLOCATE (buffer_rec(imepos)%msg(1))
      END DO

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (req_array(0:para_env%num_pe - 1, 2))

         DO imepos = 0, para_env%num_pe - 1
            CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(imepos, 1), tag=2)
            CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(imepos, 2), tag=2)
         END DO

         CALL mp_waitall(req_array(:, 1:2))

         DEALLOCATE (req_array)

      ELSE

         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      ALLOCATE (occupation(0:para_env%num_pe - 1))

      DO imepos = 0, para_env%num_pe - 1

         occupation(imepos) = buffer_rec(imepos)%msg(1)

      END DO

      max_occupation = MAXVAL(occupation)

      min_occupation = MINVAL(occupation)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%msg)
      END DO

      DEALLOCATE (buffer_send, buffer_rec, occupation)

      IF (my_one_number) THEN

         ! print fraction of non-zero blocks
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T74,ES7.1)") &
            "MATRIX_INFO| Frac. of non-zero bl. in ", matrix_name, ":", max_occupation

      ELSE

         ! print fraction of non-zero blocks
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T64,ES7.1,A,T74,ES7.1)") &
            "MATRIX_INFO| Frac. of non-zero bl. in ", matrix_name, ":", min_occupation, ' -', max_occupation

      END IF

      CALL timestop(handle)

   END SUBROUTINE print_occupation_2c

! **************************************************************************************************
!> \brief ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param tau_tj ...
!> \param num_integ_points ...
!> \param nmo ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param mat_dm_occ_global ...
!> \param mat_dm_virt_global ...
!> \param matrix_s ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param ispin ...
!> \param Eigenval ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param memory_info ...
!> \param unit_nr ...
!> \param para_env ...
!> \param jquad ...
!> \param do_mao ...
!> \param stabilize_exp ...
!> \param do_kpoints_cubic_RPA ...
!> \param qs_env ...
!> \param num_cells_dm ...
!> \param index_to_cell_dm ...
!> \param does_mat_P_T_tau_have_blocks ...
! **************************************************************************************************
   SUBROUTINE compute_mat_dm_global(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, nmo, &
                                    fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                    fm_mo_coeff_virt_scaled, mat_dm_occ_global, mat_dm_virt_global, &
                                    matrix_s, mao_coeff_occ, mao_coeff_virt, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, &
                                    unit_nr, para_env, &
                                    jquad, do_mao, stabilize_exp, do_kpoints_cubic_RPA, qs_env, &
                                    num_cells_dm, index_to_cell_dm, does_mat_P_T_tau_have_blocks)

      TYPE(cp_fm_type), POINTER                          :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      INTEGER, INTENT(IN)                                :: num_integ_points
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      INTEGER, INTENT(IN)                                :: nmo
      TYPE(cp_fm_type), POINTER                          :: fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                                            fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_virt_scaled
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ_global, mat_dm_virt_global
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, mao_coeff_occ, mao_coeff_virt
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, eps_filter
      LOGICAL, INTENT(IN)                                :: memory_info
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: jquad
      LOGICAL, INTENT(IN)                                :: do_mao
      REAL(KIND=dp), INTENT(IN)                          :: stabilize_exp
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: num_cells_dm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: does_mat_P_T_tau_have_blocks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_dm_global', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, iquad, jjB, &
                                                            ncol_local, nmao_occ, nmao_virt, &
                                                            nrow_local, size_dm_occ, size_dm_virt
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_mao_occ, col_indices, &
                                                            row_indices
      REAL(KIND=dp)                                      :: tau

      CALL timeset(routineN, handle)

      IF (memory_info .AND. unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RPA_IM_TIME_INFO| Started with time point: ", jquad

      tau = tau_tj(jquad)

      IF (do_kpoints_cubic_RPA) THEN

         ! no MAOs for kpoints
         CPASSERT(.NOT. do_mao)

         CALL compute_transl_dm(mat_dm_occ_global, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                stabilize_exp, eps_filter, num_cells_dm, index_to_cell_dm, &
                                remove_occ=.FALSE., remove_virt=.TRUE., first_jquad=1)

         CALL compute_transl_dm(mat_dm_virt_global, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                stabilize_exp, eps_filter, num_cells_dm, index_to_cell_dm, &
                                remove_occ=.TRUE., remove_virt=.FALSE., first_jquad=1)

      ELSE

         num_cells_dm = 1

         ! get info of fm_mo_coeff_occ
         CALL cp_fm_get_info(matrix=fm_mo_coeff_occ, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         ! Multiply the occupied and the virtual MO coefficients with the factor exp((-e_i-e_F)*tau/2).
         ! Then, we simply get the sum over all occ states and virt. states by a simple matrix-matrix
         ! multiplication.

         ! first, the occ
         DO jjB = 1, nrow_local
            DO iiB = 1, ncol_local
               i_global = col_indices(iiB)

               ! hard coded: exponential function gets NaN if argument is negative with large absolute value
               ! use 69, since e^(-69) = 10^(-30) which should be sufficiently small that it does not matter
               IF (ABS(tau*0.5_dp*(Eigenval(i_global) - e_fermi)) < stabilize_exp) THEN
                  fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = &
                     fm_mo_coeff_occ%local_data(jjB, iiB)*EXP(tau*0.5_dp*(Eigenval(i_global) - e_fermi))
               ELSE
                  fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = 0.0_dp
               END IF

            END DO
         END DO

         ! get info of fm_mo_coeff_virt
         CALL cp_fm_get_info(matrix=fm_mo_coeff_virt, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         ! the same for virt
         DO jjB = 1, nrow_local
            DO iiB = 1, ncol_local
               i_global = col_indices(iiB)

               IF (ABS(tau*0.5_dp*(Eigenval(i_global) - e_fermi)) < stabilize_exp) THEN
                  fm_mo_coeff_virt_scaled%local_data(jjB, iiB) = &
                     fm_mo_coeff_virt%local_data(jjB, iiB)*EXP(-tau*0.5_dp*(Eigenval(i_global) - e_fermi))
               ELSE
                  fm_mo_coeff_virt_scaled%local_data(jjB, iiB) = 0.0_dp
               END IF

            END DO
         END DO

         size_dm_occ = nmo
         size_dm_virt = nmo

         ! per default, dont do MAOs
         IF (do_mao) THEN

            CALL dbcsr_get_info(mao_coeff_occ(ispin)%matrix, &
                                nfullcols_total=nmao_occ)
            size_dm_occ = nmao_occ

            CALL dbcsr_get_info(mao_coeff_virt(ispin)%matrix, &
                                nfullcols_total=nmao_virt)
            size_dm_virt = nmao_virt

         END IF

         CALL cp_gemm(transa="N", transb="T", m=size_dm_occ, n=size_dm_occ, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mo_coeff_occ_scaled, matrix_b=fm_mo_coeff_occ_scaled, beta=0.0_dp, &
                      matrix_c=fm_scaled_dm_occ_tau)

         CALL cp_gemm(transa="N", transb="T", m=size_dm_virt, n=size_dm_virt, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mo_coeff_virt_scaled, matrix_b=fm_mo_coeff_virt_scaled, beta=0.0_dp, &
                      matrix_c=fm_scaled_dm_virt_tau)

         IF (jquad == 1) THEN

            ! transfer fm density matrices to dbcsr matrix
            NULLIFY (mat_dm_occ_global)
            CALL dbcsr_allocate_matrix_set(mat_dm_occ_global, num_integ_points, 1)

            DO iquad = 1, num_integ_points

               ALLOCATE (mat_dm_occ_global(iquad, 1)%matrix)
               ! per default, dont do MAOs
               IF (do_mao) THEN
                  CALL dbcsr_get_info(mao_coeff_occ(ispin)%matrix, &
                                      col_blk_size=blk_sizes_mao_occ)
                  CALL dbcsr_create(matrix=mat_dm_occ_global(iquad, 1)%matrix, &
                                    template=matrix_s(1)%matrix, &
                                    row_blk_size=blk_sizes_mao_occ, &
                                    col_blk_size=blk_sizes_mao_occ, nze=0, &
                                    matrix_type=dbcsr_type_no_symmetry)
               ELSE
                  CALL dbcsr_create(matrix=mat_dm_occ_global(iquad, 1)%matrix, &
                                    template=matrix_s(1)%matrix, &
                                    matrix_type=dbcsr_type_no_symmetry)
               END IF

            END DO

         END IF

         CALL copy_fm_to_dbcsr(fm_scaled_dm_occ_tau, &
                               mat_dm_occ_global(jquad, 1)%matrix, &
                               keep_sparsity=.FALSE.)

         CALL dbcsr_filter(mat_dm_occ_global(jquad, 1)%matrix, eps_filter)

         IF (jquad == 1) THEN

            NULLIFY (mat_dm_virt_global)
            CALL dbcsr_allocate_matrix_set(mat_dm_virt_global, num_integ_points, 1)

         END IF

         ALLOCATE (mat_dm_virt_global(jquad, 1)%matrix)
         CALL dbcsr_create(matrix=mat_dm_virt_global(jquad, 1)%matrix, &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_scaled_dm_virt_tau, &
                               mat_dm_virt_global(jquad, 1)%matrix, &
                               keep_sparsity=.FALSE.)

         CALL dbcsr_filter(mat_dm_virt_global(jquad, 1)%matrix, eps_filter)

         ! release memory
         IF (jquad > 1) THEN
            CALL dbcsr_set(mat_dm_occ_global(jquad - 1, 1)%matrix, 0.0_dp)
            CALL dbcsr_set(mat_dm_virt_global(jquad - 1, 1)%matrix, 0.0_dp)
            CALL dbcsr_filter(mat_dm_occ_global(jquad - 1, 1)%matrix, 0.0_dp)
            CALL dbcsr_filter(mat_dm_virt_global(jquad - 1, 1)%matrix, 0.0_dp)
         END IF

         IF (memory_info) THEN
            CALL print_occupation_2c(mat_dm_occ_global(jquad, 1)%matrix, unit_nr, &
                                     "D_occ(it)", para_env, one_number=.TRUE.)
            CALL print_occupation_2c(mat_dm_virt_global(jquad, 1)%matrix, unit_nr, &
                                     "D_virt(it)", para_env, one_number=.TRUE.)
         END IF

      END IF ! do kpoints

      IF (jquad == 1) THEN
         ALLOCATE (does_mat_P_T_tau_have_blocks(num_cells_dm/2 + 1))
         does_mat_P_T_tau_have_blocks(:) = .TRUE.
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_mat_dm_global

! **************************************************************************************************
!> \brief ...
!> \param cycle_due_to_sparse_dm ...
!> \param mat_dm_occ_global ...
!> \param mat_dm_virt_global ...
!> \param mat_3c_overl_int_cut ...
!> \param num_integ_points ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param starts_array_cm_mao_occ ...
!> \param starts_array_cm_mao_virt ...
!> \param ends_array_cm_mao_occ ...
!> \param ends_array_cm_mao_virt ...
!> \param my_group_L_sizes_im_time ...
!> \param para_env ...
!> \param jquad ...
! **************************************************************************************************
   SUBROUTINE get_cycle_due_to_sparse_dm(cycle_due_to_sparse_dm, mat_dm_occ_global, mat_dm_virt_global, &
                                         mat_3c_overl_int_cut, num_integ_points, cut_memory, cut_RI, &
                                         starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
                                         ends_array_cm_mao_occ, ends_array_cm_mao_virt, &
                                         my_group_L_sizes_im_time, para_env, jquad)

      LOGICAL, DIMENSION(:, :, :), INTENT(OUT)           :: cycle_due_to_sparse_dm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ_global, mat_dm_virt_global
      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut
      INTEGER, INTENT(IN)                                :: num_integ_points, cut_memory, cut_RI
      INTEGER, DIMENSION(:), INTENT(IN) :: starts_array_cm_mao_occ, starts_array_cm_mao_virt, &
         ends_array_cm_mao_occ, ends_array_cm_mao_virt, my_group_L_sizes_im_time
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: jquad

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_cycle_due_to_sparse_dm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, col_offset, col_size, handle, &
                                                            i_cut_RI, i_mem, j_mem, &
                                                            my_group_L_size, natom, row, &
                                                            row_offset, row_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cycle_combined, cycle_due_to_sparse_dm_occ_tmp, &
         cycle_due_to_sparse_dm_virt_tmp, non_zero_blocks_in_dm_occ, non_zero_blocks_in_dm_virt
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(mat_dm_occ_global(1, 1)%matrix, &
                          nblkrows_total=natom)

      ALLOCATE (non_zero_blocks_in_dm_occ(num_integ_points, cut_memory, natom))
      non_zero_blocks_in_dm_occ = 0

      ALLOCATE (non_zero_blocks_in_dm_virt(num_integ_points, cut_memory, natom))
      non_zero_blocks_in_dm_virt = 0

      ! for occ
      CALL get_non_zero_blocks_dm(mat_dm_occ_global, non_zero_blocks_in_dm_occ, &
                                  cut_memory, natom, jquad, &
                                  starts_array_cm_mao_occ, ends_array_cm_mao_occ, para_env)

      ! for virt
      CALL get_non_zero_blocks_dm(mat_dm_virt_global, non_zero_blocks_in_dm_virt, &
                                  cut_memory, natom, jquad, &
                                  starts_array_cm_mao_virt, ends_array_cm_mao_virt, para_env)

      ALLOCATE (cycle_due_to_sparse_dm_occ_tmp(cut_memory, cut_memory, num_integ_points))
      ALLOCATE (cycle_due_to_sparse_dm_virt_tmp(cut_memory, cut_memory, num_integ_points))
      cycle_due_to_sparse_dm_occ_tmp = 0
      cycle_due_to_sparse_dm_virt_tmp = 0

      DO i_cut_RI = 1, cut_RI

         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         DO i_mem = 1, cut_memory

            CALL dbcsr_iterator_start(iter, mat_3c_overl_int_cut(i_cut_RI, i_mem, 1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                              row_offset=row_offset, col_offset=col_offset, &
                                              row_size=row_size, col_size=col_size)

               DO j_mem = 1, cut_memory

                  ! check whether row index is restricted due to mem_cut and col index due to sparse virt dm
                  IF (row_offset + row_size - 1 >= starts_array_cm_mao_virt(i_mem) .AND. &
                      row_offset < ends_array_cm_mao_virt(i_mem) .AND. &
                      ANY(non_zero_blocks_in_dm_virt(jquad, j_mem, :) == col)) THEN

                     cycle_due_to_sparse_dm_virt_tmp(j_mem, i_mem, jquad) = 1

                  END IF

                  ! check whether row index is restricted due to sparse occ dm and col index due to j_mem
                  ! use the fact that mat_3c_overl_int_cut is symmetric (neglecting the combined col index)
                  IF (row_offset + row_size - 1 >= starts_array_cm_mao_occ(i_mem) .AND. &
                      row_offset < ends_array_cm_mao_occ(i_mem) .AND. &
                      ANY(non_zero_blocks_in_dm_occ(jquad, j_mem, :) == col)) THEN

                     cycle_due_to_sparse_dm_occ_tmp(i_mem, j_mem, jquad) = 1

                  END IF

               END DO ! j_mem

            END DO ! block

            CALL dbcsr_iterator_stop(iter)

         END DO ! i_mem

      END DO ! cut_RI

      CALL mp_sum(cycle_due_to_sparse_dm_occ_tmp, para_env%group)
      CALL mp_sum(cycle_due_to_sparse_dm_virt_tmp, para_env%group)

      DO i_mem = 1, cut_memory
         DO j_mem = 1, cut_memory

            IF (cycle_due_to_sparse_dm_occ_tmp(i_mem, j_mem, jquad) > 0) THEN
               cycle_due_to_sparse_dm_occ_tmp(i_mem, j_mem, jquad) = 1
            END IF

            IF (cycle_due_to_sparse_dm_virt_tmp(i_mem, j_mem, jquad) > 0) THEN
               cycle_due_to_sparse_dm_virt_tmp(i_mem, j_mem, jquad) = 1
            END IF

         END DO
      END DO

      ALLOCATE (cycle_combined(cut_memory, cut_memory, num_integ_points))
      cycle_combined = 0

      cycle_combined(:, :, :) = cycle_due_to_sparse_dm_occ_tmp(:, :, :) + &
                                cycle_due_to_sparse_dm_virt_tmp(:, :, :)

      DO i_mem = 1, cut_memory
         DO j_mem = 1, cut_memory

            IF (cycle_combined(i_mem, j_mem, jquad) == 2) THEN
               cycle_due_to_sparse_dm(i_mem, j_mem, jquad) = .FALSE.
            ELSE
               cycle_due_to_sparse_dm(i_mem, j_mem, jquad) = .TRUE.
            END IF

         END DO
      END DO

      DEALLOCATE (non_zero_blocks_in_dm_occ, non_zero_blocks_in_dm_virt, &
                  cycle_due_to_sparse_dm_occ_tmp, &
                  cycle_due_to_sparse_dm_virt_tmp, cycle_combined)

      CALL timestop(handle)

   END SUBROUTINE get_cycle_due_to_sparse_dm

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_global ...
!> \param non_zero_blocks_in_dm ...
!> \param cut_memory ...
!> \param natom ...
!> \param jquad ...
!> \param starts_array_cm ...
!> \param ends_array_cm ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE get_non_zero_blocks_dm(mat_dm_global, non_zero_blocks_in_dm, &
                                     cut_memory, natom, jquad, &
                                     starts_array_cm, ends_array_cm, para_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_global
      INTEGER, DIMENSION(:, :, :), INTENT(INOUT)         :: non_zero_blocks_in_dm
      INTEGER, INTENT(IN)                                :: cut_memory, natom, jquad
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_cm, ends_array_cm
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_non_zero_blocks_dm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, handle, i_mem, row, row_offset, &
                                                            row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      CALL dbcsr_iterator_start(iter, mat_dm_global(jquad, 1)%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, row_offset=row_offset)

         DO i_mem = 1, cut_memory

            ! check whether block is in the range of the memory cutoff
            IF (row_offset + row_size - 1 >= starts_array_cm(i_mem) .AND. row_offset < ends_array_cm(i_mem)) THEN

               ! 1 means that the block is there
               non_zero_blocks_in_dm(jquad, i_mem, col) = 1

            END IF

         END DO

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL mp_sum(non_zero_blocks_in_dm, para_env%group)

      ! write in non_zero_blocks_in_dm the actual row number
      DO i_mem = 1, cut_memory

         DO row = 1, natom

            IF (non_zero_blocks_in_dm(jquad, i_mem, row) > 0) THEN

               non_zero_blocks_in_dm(jquad, i_mem, row) = row

            END IF

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE get_non_zero_blocks_dm

! **************************************************************************************************
!> \brief ...
!> \param mat_M_P_munu ...
!> \param mat_M_mu_Pnu ...
!> \param para_env ...
!> \param my_group_L_starts_im_time ...
!> \param my_group_L_sizes_im_time ...
!> \param mepos_P_from_RI_row ...
!> \param n_group_col ...
!> \param cut_RI ...
!> \param row_from_LLL ...
!> \param starts_array_prim_row ...
!> \param ends_array_prim_row ...
!> \param starts_array_prim_col ...
!> \param ends_array_prim_col ...
!> \param offset_combi_block ...
!> \param starts_array_prim_fullcol ...
!> \param ends_array_prim_fullcol ...
!> \param starts_array_prim_fullrow ...
!> \param ends_array_prim_fullrow ...
!> \param group_size_P ...
!> \param i_mem ...
!> \param j_mem ...
!> \param buffer_mat_M ...
!> \param eps_filter_im_time ...
!> \param do_occ ...
!> \param do_virt ...
! **************************************************************************************************
   SUBROUTINE fill_mat_M_P_munu_from_M_mu_Pnu(mat_M_P_munu, mat_M_mu_Pnu, para_env, my_group_L_starts_im_time, &
                                              my_group_L_sizes_im_time, mepos_P_from_RI_row, &
                                              n_group_col, cut_RI, row_from_LLL, &
                                              starts_array_prim_row, ends_array_prim_row, &
                                              starts_array_prim_col, ends_array_prim_col, &
                                              offset_combi_block, starts_array_prim_fullcol, &
                                              ends_array_prim_fullcol, starts_array_prim_fullrow, &
                                              ends_array_prim_fullrow, group_size_P, &
                                              i_mem, j_mem, buffer_mat_M, eps_filter_im_time, &
                                              do_occ, do_virt)

      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_M_P_munu
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_M_mu_Pnu
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: my_group_L_starts_im_time, &
                                                            my_group_L_sizes_im_time, &
                                                            mepos_P_from_RI_row
      INTEGER, INTENT(IN)                                :: n_group_col, cut_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_from_LLL
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN) :: starts_array_prim_row, &
         ends_array_prim_row, starts_array_prim_col, ends_array_prim_col, offset_combi_block, &
         starts_array_prim_fullcol, ends_array_prim_fullcol, starts_array_prim_fullrow, &
         ends_array_prim_fullrow
      INTEGER, INTENT(IN)                                :: group_size_P, i_mem, j_mem
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: buffer_mat_M
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter_im_time
      LOGICAL, INTENT(IN)                                :: do_occ, do_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_mat_M_P_munu_from_M_mu_Pnu', &
         routineP = moduleN//':'//routineN

      INTEGER :: block, block_counter_int, block_size, col, col_end_in_data_block, col_offset, &
         col_offset_data_block, col_offset_orig, col_offset_prim, col_offset_rec_prim, &
         col_rec_prim, col_rec_prim_rel, col_size, col_size_in_data_block, col_size_orig, &
         col_size_to_send, col_start_in_data_block, color_sub_col, color_sub_P, color_sub_row, &
         end_col, end_row, handle, handle1, i_cut_RI, imepos, LLL, mepos_P, my_group_L_size, &
         my_group_L_start, n_entries_rec, num_blocks, offset, offset_rec, old_block, rec_counter, &
         row, row_end_in_data_block, row_offset, row_offset_data_block, row_offset_prim
      INTEGER :: row_offset_rec_prim, row_rec_prim, row_rec_prim_rel, row_size, &
         row_size_in_data_block, row_size_to_send, row_start_in_data_block, send_counter, &
         start_col, start_row
      INTEGER, ALLOCATABLE, DIMENSION(:) :: block_counter, cols_to_allocate, entry_counter, &
         num_blocks_rec, num_blocks_send, num_entries_rec, num_entries_send, rows_to_allocate, &
         sizes_rec, sizes_send
      INTEGER, DIMENSION(:, :), POINTER                  :: req_array
      LOGICAL                                            :: is_new_block
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send

      CPASSERT(do_occ .NEQV. do_virt)

      CALL timeset(routineN, handle)

      CALL timeset("allocate_stuff_M", handle1)

      CALL dbcsr_set(mat_M_P_munu%matrix, 0.0_dp)

      NULLIFY (data_block)

      ALLOCATE (num_entries_send(0:para_env%num_pe - 1))
      num_entries_send(:) = 0

      ALLOCATE (num_blocks_send(0:para_env%num_pe - 1))
      num_blocks_send(:) = 0

      CALL timestop(handle1)

      CALL timeset("get_sizes_M", handle1)

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         CALL dbcsr_iterator_start(iter, mat_M_mu_Pnu(i_cut_RI)%matrix)

         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            ! process to send to
            DO imepos = 0, para_env%num_pe - 1

               color_sub_P = imepos/group_size_P
               mepos_P = MODULO(imepos, group_size_P)
               color_sub_row = color_sub_P/n_group_col
               color_sub_col = MODULO(color_sub_P, n_group_col)

               IF (do_occ) THEN

                  IF (row >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      row <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      col >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      col <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     IF (row == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_row = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_row = ends_array_prim_fullrow(color_sub_row, i_mem)
                        row_size_to_send = end_row - start_row + 1

                     ELSE IF (row == starts_array_prim_row(color_sub_row, i_mem)) THEN

                        start_row = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_row = row_offset + row_size - 1
                        row_size_to_send = end_row - start_row + 1

                     ELSE IF (row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_row = row_offset
                        end_row = ends_array_prim_fullrow(color_sub_row, i_mem)
                        row_size_to_send = end_row - start_row + 1

                     ELSE

                        row_size_to_send = row_size

                     END IF

                     IF (col == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_col = (starts_array_prim_fullcol(color_sub_col, j_mem) - 1)*my_group_L_size + 1
                        end_col = ends_array_prim_fullcol(color_sub_col, j_mem)*my_group_L_size
                        col_size_to_send = end_col - start_col + 1

                     ELSE IF (col == starts_array_prim_col(color_sub_col, j_mem)) THEN

                        start_col = (starts_array_prim_fullcol(color_sub_col, j_mem) - 1)*my_group_L_size + 1
                        end_col = col_offset + col_size - 1
                        col_size_to_send = end_col - start_col + 1

                     ELSE IF (col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_col = col_offset
                        end_col = ends_array_prim_fullcol(color_sub_col, j_mem)*my_group_L_size
                        col_size_to_send = end_col - start_col + 1

                     ELSE

                        col_size_to_send = col_size

                     END IF

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P == mepos_P_from_RI_row(row_from_LLL(LLL + my_group_L_start - 1))) THEN

                           num_entries_send(imepos) = num_entries_send(imepos) + &
                                                      row_size_to_send*col_size_to_send/my_group_L_size

                           num_blocks_send(imepos) = num_blocks_send(imepos) + 1

                        END IF

                     END DO

                  END IF

                  ! everything transposed for virtuals
               ELSE IF (do_virt) THEN

                  col_size_orig = col_size/my_group_L_size

                  !             col_offset_orig = col_offset/my_group_L_size

                  col_offset_orig = (col_offset - 1)/my_group_L_size + 1

                  IF (col >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      col <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      row >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      row <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     IF (col == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_col = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_col = ends_array_prim_fullrow(color_sub_row, i_mem)
                        col_size_to_send = (end_col - start_col + 1)*my_group_L_size

                     ELSE IF (col == starts_array_prim_row(color_sub_row, i_mem)) THEN

                        start_col = starts_array_prim_fullrow(color_sub_row, i_mem)
                        end_col = col_offset_orig + col_size_orig - 1
                        col_size_to_send = (end_col - start_col + 1)*my_group_L_size

                     ELSE IF (col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        start_col = col_offset_orig
                        end_col = ends_array_prim_fullrow(color_sub_row, i_mem)
                        col_size_to_send = (end_col - start_col + 1)*my_group_L_size

                     ELSE

                        col_size_to_send = col_size

                     END IF

                     IF (row == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_row = starts_array_prim_fullcol(color_sub_col, j_mem)
                        end_row = ends_array_prim_fullcol(color_sub_col, j_mem)
                        row_size_to_send = end_row - start_row + 1

                     ELSE IF (row == starts_array_prim_col(color_sub_col, j_mem)) THEN

                        start_row = starts_array_prim_fullcol(color_sub_col, j_mem)
                        end_row = row_offset + row_size - 1
                        row_size_to_send = end_row - start_row + 1

                     ELSE IF (row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        start_row = row_offset
                        end_row = ends_array_prim_fullcol(color_sub_col, j_mem)
                        row_size_to_send = end_row - start_row + 1

                     ELSE

                        row_size_to_send = row_size

                     END IF

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P == mepos_P_from_RI_row(row_from_LLL(LLL + my_group_L_start - 1))) THEN

                           num_entries_send(imepos) = num_entries_send(imepos) + &
                                                      row_size_to_send*col_size_to_send/my_group_L_size

                           num_blocks_send(imepos) = num_blocks_send(imepos) + 1

                        END IF

                     END DO

!                 num_entries_send(imepos) = num_entries_send(imepos) + row_size_to_send*col_size_to_send
!                 num_blocks_send(imepos)  = num_blocks_send(imepos) + my_group_L_size

                  END IF

               END IF

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

      END DO ! i_cut_RI

      CALL timestop(handle1)

      CALL timeset("send_sizes_M", handle1)

      ALLOCATE (buffer_rec(0:para_env%num_pe - 1))
      ALLOCATE (buffer_send(0:para_env%num_pe - 1))

      ALLOCATE (num_entries_rec(0:para_env%num_pe - 1))
      ALLOCATE (num_blocks_rec(0:para_env%num_pe - 1))

      IF (para_env%num_pe > 1) THEN

         ALLOCATE (sizes_rec(0:2*para_env%num_pe - 1))
         ALLOCATE (sizes_send(0:2*para_env%num_pe - 1))

         DO imepos = 0, para_env%num_pe - 1
            sizes_send(2*imepos) = num_entries_send(imepos)
            sizes_send(2*imepos + 1) = num_blocks_send(imepos)
         END DO

         CALL mp_alltoall(sizes_send, sizes_rec, 2, para_env%group)

         DO imepos = 0, para_env%num_pe - 1
            num_entries_rec(imepos) = sizes_rec(2*imepos)
            num_blocks_rec(imepos) = sizes_rec(2*imepos + 1)
         END DO

         DEALLOCATE (sizes_rec, sizes_send)

      ELSE

         num_entries_rec(0) = num_entries_send(0)
         num_blocks_rec(0) = num_blocks_send(0)

      END IF

      CALL timestop(handle1)

      CALL timeset("fill_buffer_send_M", handle1)

      ! allocate data message and corresponding indices
      DO imepos = 0, para_env%num_pe - 1

         ALLOCATE (buffer_rec(imepos)%msg(num_entries_rec(imepos)))
         buffer_rec(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_send(imepos)%msg(num_entries_send(imepos)))
         buffer_send(imepos)%msg = 0.0_dp

         ALLOCATE (buffer_rec(imepos)%indx(num_blocks_rec(imepos), 10))
         buffer_rec(imepos)%indx = 0

         ALLOCATE (buffer_send(imepos)%indx(num_blocks_send(imepos), 10))
         buffer_send(imepos)%indx = 0

      END DO

      ALLOCATE (entry_counter(0:para_env%num_pe - 1))
      entry_counter(:) = 0

      ALLOCATE (block_counter(0:para_env%num_pe - 1))
      block_counter(:) = 0

      DO i_cut_RI = 1, cut_RI

         my_group_L_start = my_group_L_starts_im_time(i_cut_RI)
         my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

         ! fill buffer_send
         CALL dbcsr_iterator_start(iter, mat_M_mu_Pnu(i_cut_RI)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))

            CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)

            DO imepos = 0, para_env%num_pe - 1

!             color_sub_row = imepos/n_group_col
!             color_sub_col = MODULO(imepos,n_group_col)

               color_sub_P = imepos/group_size_P
               mepos_P = MODULO(imepos, group_size_P)
               color_sub_row = color_sub_P/n_group_col
               color_sub_col = MODULO(color_sub_P, n_group_col)

               IF (do_occ) THEN

                  IF (row >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      row <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      col >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      col <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     col_size_orig = col_size/my_group_L_size

                     ! For terminal blocks, we have to compute the sizes
                     IF (row == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem) - row_offset + 1
                        row_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem) - row_offset + 1

                     ELSE IF (row == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              row .NE. ends_array_prim_row(color_sub_row, i_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem) - row_offset + 1
                        row_end_in_data_block = row_size

                     ELSE IF (row .NE. starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              row == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        row_start_in_data_block = 1
                        row_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem) - row_offset + 1

                     ELSE

                        row_start_in_data_block = 1
                        row_end_in_data_block = row_size
                        row_size_in_data_block = row_size

                     END IF

                     row_size_in_data_block = row_end_in_data_block - row_start_in_data_block + 1

                     col_offset_orig = (col_offset - 1)/my_group_L_size + 1

                     ! For terminal blocks, we have to compute the sizes
                     IF (col == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem) - col_offset_orig + 1
                        col_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem) - col_offset_orig + 1

                     ELSE IF (col == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              col .NE. ends_array_prim_col(color_sub_col, j_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem) - col_offset_orig + 1
                        col_end_in_data_block = col_size_orig

                     ELSE IF (col .NE. starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              col == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        col_start_in_data_block = 1
                        col_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem) - col_offset_orig + 1

                     ELSE

                        col_start_in_data_block = 1
                        col_end_in_data_block = col_size_orig

                     END IF

                     col_size_in_data_block = col_end_in_data_block - col_start_in_data_block + 1

                     block_size = row_size_in_data_block*col_size_in_data_block

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P .NE. mepos_P_from_RI_row(row_from_LLL(LLL + my_group_L_start - 1))) CYCLE

                        offset = entry_counter(imepos)

                        col_offset_data_block = (LLL - 1)*col_size_orig + col_start_in_data_block

                        buffer_send(imepos)%msg(offset + 1:offset + block_size) = &
                           RESHAPE(data_block(row_start_in_data_block:row_end_in_data_block, &
                                              col_offset_data_block:col_offset_data_block + col_size_in_data_block - 1), &
                                   (/block_size/))

                        block = block_counter(imepos) + 1

                        buffer_send(imepos)%indx(block, 1) = LLL + my_group_L_start - 1
                        buffer_send(imepos)%indx(block, 2) = row
                        buffer_send(imepos)%indx(block, 3) = row_offset
                        buffer_send(imepos)%indx(block, 4) = row_size
                        buffer_send(imepos)%indx(block, 5) = col
                        buffer_send(imepos)%indx(block, 6) = col_offset
                        buffer_send(imepos)%indx(block, 7) = col_size_orig
                        buffer_send(imepos)%indx(block, 8) = offset
                        buffer_send(imepos)%indx(block, 9) = block_size

                        entry_counter(imepos) = entry_counter(imepos) + block_size

                        block_counter(imepos) = block_counter(imepos) + 1

                     END DO

                  END IF

               ELSE IF (do_virt) THEN

                  IF (col >= starts_array_prim_row(color_sub_row, i_mem) .AND. &
                      col <= ends_array_prim_row(color_sub_row, i_mem) .AND. &
                      row >= starts_array_prim_col(color_sub_col, j_mem) .AND. &
                      row <= ends_array_prim_col(color_sub_col, j_mem)) THEN

                     col_size_orig = col_size/my_group_L_size

                     col_offset_orig = (col_offset - 1)/my_group_L_size + 1

                     ! For terminal blocks, we have to compute the sizes
                     IF (col == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                         col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem) - col_offset_orig + 1
                        col_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem) - col_offset_orig + 1

                     ELSE IF (col == starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              col .NE. ends_array_prim_row(color_sub_row, i_mem)) THEN

                        col_start_in_data_block = starts_array_prim_fullrow(color_sub_row, i_mem) - col_offset_orig + 1
                        col_end_in_data_block = col_size_orig

                     ELSE IF (col .NE. starts_array_prim_row(color_sub_row, i_mem) .AND. &
                              col == ends_array_prim_row(color_sub_row, i_mem)) THEN

                        col_start_in_data_block = 1
                        col_end_in_data_block = ends_array_prim_fullrow(color_sub_row, i_mem) - col_offset_orig + 1

                     ELSE

                        col_start_in_data_block = 1
                        col_end_in_data_block = col_size_orig

                     END IF

                     col_size_in_data_block = col_end_in_data_block - col_start_in_data_block + 1

                     ! For terminal blocks, we have to compute the sizes
                     IF (row == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                         row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem) - row_offset + 1
                        row_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem) - row_offset + 1

                     ELSE IF (row == starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              row .NE. ends_array_prim_col(color_sub_col, j_mem)) THEN

                        row_start_in_data_block = starts_array_prim_fullcol(color_sub_col, j_mem) - row_offset + 1
                        row_end_in_data_block = row_size

                     ELSE IF (row .NE. starts_array_prim_col(color_sub_col, j_mem) .AND. &
                              row == ends_array_prim_col(color_sub_col, j_mem)) THEN

                        row_start_in_data_block = 1
                        row_end_in_data_block = ends_array_prim_fullcol(color_sub_col, j_mem) - row_offset + 1

                     ELSE

                        row_start_in_data_block = 1
                        row_end_in_data_block = row_size
                        row_size_in_data_block = row_size

                     END IF

                     row_size_in_data_block = row_end_in_data_block - row_start_in_data_block + 1

                     block_size = row_size_in_data_block*col_size_in_data_block

                     DO LLL = 1, my_group_L_size

                        IF (mepos_P .NE. mepos_P_from_RI_row(row_from_LLL(LLL + my_group_L_start - 1))) CYCLE

                        offset = entry_counter(imepos)

                        col_offset_data_block = (LLL - 1)*col_size_orig + col_start_in_data_block

                        buffer_send(imepos)%msg(offset + 1:offset + block_size) = &
                           RESHAPE(TRANSPOSE(data_block(row_start_in_data_block:row_end_in_data_block, &
                                                        col_offset_data_block:col_offset_data_block + &
                                                        col_size_in_data_block - 1)), &
                                   (/block_size/))

                        block = block_counter(imepos) + 1

                        buffer_send(imepos)%indx(block, 1) = LLL + my_group_L_start - 1
                        buffer_send(imepos)%indx(block, 2) = col
                        buffer_send(imepos)%indx(block, 3) = col_offset
                        buffer_send(imepos)%indx(block, 4) = col_size_orig
                        buffer_send(imepos)%indx(block, 5) = row
                        buffer_send(imepos)%indx(block, 6) = row_offset
                        buffer_send(imepos)%indx(block, 7) = row_size
                        buffer_send(imepos)%indx(block, 8) = offset
                        buffer_send(imepos)%indx(block, 9) = block_size

                        entry_counter(imepos) = entry_counter(imepos) + block_size

                        block_counter(imepos) = block_counter(imepos) + 1

                     END DO

                  END IF

               END IF

            END DO

         END DO

         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_set(mat_M_mu_Pnu(i_cut_RI)%matrix, 0.0_dp)
         CALL dbcsr_filter(mat_M_mu_Pnu(i_cut_RI)%matrix, 1.0_dp)

      END DO ! i_cut_RI

      CALL timestop(handle1)

      CALL timeset("comm_data_M", handle1)

      DEALLOCATE (entry_counter, block_counter)

      ALLOCATE (req_array(1:para_env%num_pe, 4))

      IF (para_env%num_pe > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_rec(imepos) > 0) THEN
               rec_counter = rec_counter + 1
               CALL mp_irecv(buffer_rec(imepos)%indx, imepos, para_env%group, req_array(rec_counter, 3), tag=4)
            END IF
            IF (num_entries_rec(imepos) > 0) THEN
               CALL mp_irecv(buffer_rec(imepos)%msg, imepos, para_env%group, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO imepos = 0, para_env%num_pe - 1
            IF (num_entries_send(imepos) > 0) THEN
               send_counter = send_counter + 1
               CALL mp_isend(buffer_send(imepos)%indx, imepos, para_env%group, req_array(send_counter, 1), tag=4)
            END IF
            IF (num_entries_send(imepos) > 0) THEN
               CALL mp_isend(buffer_send(imepos)%msg, imepos, para_env%group, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         CALL mp_waitall(req_array(1:send_counter, 1:2))
         CALL mp_waitall(req_array(1:rec_counter, 3:4))

      ELSE

         buffer_rec(0)%indx = buffer_send(0)%indx
         buffer_rec(0)%msg = buffer_send(0)%msg

      END IF

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_send(imepos)%msg)
         DEALLOCATE (buffer_send(imepos)%indx)
      END DO

      CALL timestop(handle1)

      CALL timeset("reserve_blocks_M", handle1)

      num_blocks = 0

      ! get the number of blocks, which have to be allocated
      DO imepos = 0, para_env%num_pe - 1

         num_blocks = num_blocks + num_blocks_rec(imepos)

      END DO

      ALLOCATE (rows_to_allocate(num_blocks))
      rows_to_allocate = 0
      ALLOCATE (cols_to_allocate(num_blocks))
      cols_to_allocate = 1

      block_counter_int = 0

      DO imepos = 0, para_env%num_pe - 1

         DO block = 1, num_blocks_rec(imepos)

            is_new_block = .TRUE.

            ! check whether block is already there
            DO old_block = 1, block_counter_int

               IF (row_from_LLL(buffer_rec(imepos)%indx(block, 1)) == rows_to_allocate(old_block)) THEN

                  is_new_block = .FALSE.

               END IF

            END DO

            IF (is_new_block .AND. buffer_rec(imepos)%indx(block, 1) .NE. 0) THEN

               block_counter_int = block_counter_int + 1

               rows_to_allocate(block_counter_int) = row_from_LLL(buffer_rec(imepos)%indx(block, 1))

            END IF

         END DO

      END DO

      CALL dbcsr_set(mat_M_P_munu%matrix, 0.0_dp)

      CALL dbcsr_filter(mat_M_P_munu%matrix, 1.0_dp)

      CALL dbcsr_reserve_blocks(mat_M_P_munu%matrix, rows=rows_to_allocate(1:block_counter_int), &
                                cols=cols_to_allocate(1:block_counter_int))

      CALL dbcsr_finalize(mat_M_P_munu%matrix)

      CALL timestop(handle1)

      CALL dbcsr_set(mat_M_P_munu%matrix, 0.0_dp)

      n_entries_rec = SUM(num_entries_rec)

      ! Fill the dbcsr matrix
      CALL timeset("fill_dbcsr_mat_M", handle1)

      color_sub_P = para_env%mepos/group_size_P

      color_sub_row = color_sub_P/n_group_col
      color_sub_col = MODULO(para_env%mepos, n_group_col)

      row_offset_prim = starts_array_prim_row(color_sub_row, i_mem)
      col_offset_prim = starts_array_prim_col(color_sub_col, j_mem)

      CALL dbcsr_iterator_start(iter, mat_M_P_munu%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)

         buffer_mat_M(:, :) = 0.0_dp

         DO imepos = 0, para_env%num_pe - 1

            DO block = 1, num_blocks_rec(imepos)

               LLL = buffer_rec(imepos)%indx(block, 1)

               IF (LLL >= row_offset .AND. LLL < row_offset + row_size) THEN

                  row_rec_prim = buffer_rec(imepos)%indx(block, 2)
                  row_offset_rec_prim = buffer_rec(imepos)%indx(block, 3)
                  col_rec_prim = buffer_rec(imepos)%indx(block, 5)
                  col_offset_rec_prim = buffer_rec(imepos)%indx(block, 6)
                  offset_rec = buffer_rec(imepos)%indx(block, 8)
                  block_size = buffer_rec(imepos)%indx(block, 9)

                  row_rec_prim_rel = row_rec_prim - row_offset_prim + 1
                  col_rec_prim_rel = col_rec_prim - col_offset_prim + 1

                  row_offset_data_block = LLL - row_offset

                  col_offset_data_block = offset_combi_block(row_rec_prim, col_rec_prim)

                  buffer_mat_M(row_offset_data_block + 1, col_offset_data_block + 1:col_offset_data_block + block_size) = &
                     buffer_rec(imepos)%msg(offset_rec + 1:offset_rec + block_size)

               END IF

            END DO

         END DO

         data_block(1:row_size, 1:col_size) = buffer_mat_M(1:row_size, 1:col_size)

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_filter(mat_M_P_munu%matrix, eps_filter_im_time)

      DO imepos = 0, para_env%num_pe - 1
         DEALLOCATE (buffer_rec(imepos)%msg)
         DEALLOCATE (buffer_rec(imepos)%indx)
      END DO

      DEALLOCATE (req_array)

      DEALLOCATE (buffer_rec, buffer_send)

      DEALLOCATE (num_entries_send, num_entries_rec, num_blocks_send, num_blocks_rec)

      CALL timestop(handle1)

      CALL timestop(handle)

   END SUBROUTINE fill_mat_M_P_munu_from_M_mu_Pnu

! **************************************************************************************************
!> \brief ...
!> \param mat_3c_overl_int_cut ...
!> \param mat_3c_overl_int ...
!> \param cut_memory ...
!> \param cut_RI ...
!> \param starts_array_cm ...
!> \param ends_array_cm ...
!> \param my_group_L_sizes_im_time ...
!> \param eps_filter ...
!> \param do_kpoints_cubic_RPA ...
!> \param do_gw_im_time ...
!> \param num_3c_repl ...
! **************************************************************************************************
   SUBROUTINE setup_mat_for_mem_cut_3c(mat_3c_overl_int_cut, mat_3c_overl_int, cut_memory, cut_RI, &
                                       starts_array_cm, ends_array_cm, my_group_L_sizes_im_time, &
                                       eps_filter, do_kpoints_cubic_RPA, do_gw_im_time, num_3c_repl)

      TYPE(dbcsr_p_type), DIMENSION(:, :, :, :), POINTER :: mat_3c_overl_int_cut
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int
      INTEGER, INTENT(IN)                                :: cut_memory, cut_RI
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_cm, ends_array_cm, &
                                                            my_group_L_sizes_im_time
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA, do_gw_im_time
      INTEGER, INTENT(IN)                                :: num_3c_repl

      CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_mat_for_mem_cut_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, col, col_end_in_data_block, col_offset, col_size, col_start_in_data_block, &
         handle, i_cell, i_cut_RI, i_mem, j_cell, my_group_L_size, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      NULLIFY (mat_3c_overl_int_cut)
      CALL dbcsr_allocate_matrix_set(mat_3c_overl_int_cut, cut_RI, cut_memory, num_3c_repl, num_3c_repl)

      DO i_mem = 1, cut_memory
         DO i_cut_RI = 1, cut_RI
            DO i_cell = 1, num_3c_repl
               DO j_cell = 1, num_3c_repl

                  ALLOCATE (mat_3c_overl_int_cut(i_cut_RI, i_mem, i_cell, j_cell)%matrix)
                  CALL dbcsr_create(matrix=mat_3c_overl_int_cut(i_cut_RI, i_mem, i_cell, j_cell)%matrix, &
                                    template=mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix)

                  CALL dbcsr_copy(mat_3c_overl_int_cut(i_cut_RI, i_mem, i_cell, j_cell)%matrix, &
                                  mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix)

                  my_group_L_size = my_group_L_sizes_im_time(i_cut_RI)

                  CALL dbcsr_iterator_start(iter, mat_3c_overl_int_cut(i_cut_RI, i_mem, i_cell, j_cell)%matrix)
                  DO WHILE (dbcsr_iterator_blocks_left(iter))
                     CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                                    col_size=col_size, col_offset=col_offset)

                     ! set the block to zero if it is outside the range for cutting the mu-sigma combined index
                     IF (col_offset > ends_array_cm(i_mem)*my_group_L_size .OR. &
                         col_offset + col_size - 1 < (starts_array_cm(i_mem) - 1)*my_group_L_size + 1) THEN

                        data_block = 0.0_dp

                     END IF

                     IF (col_offset >= ends_array_cm(i_mem)*my_group_L_size .AND. &
                         col_offset + col_size - 1 < ends_array_cm(i_mem)*my_group_L_size) THEN

                        col_end_in_data_block = ends_array_cm(i_mem)*my_group_L_size - col_offset + 1

                        data_block(:, col_end_in_data_block + 1:col_size) = 0.0_dp

                     END IF

                     IF (col_offset > (starts_array_cm(i_mem) - 1)*my_group_L_size + 1 .AND. &
                         col_offset + col_size - 1 <= (starts_array_cm(i_mem) - 1)*my_group_L_size + 1) THEN

                        col_start_in_data_block = (starts_array_cm(i_mem) - 1)*my_group_L_size + 1

                        data_block(:, 1:col_start_in_data_block - 1) = 0.0_dp

                     END IF

                  END DO

                  CALL dbcsr_iterator_stop(iter)

                  ! remove the zeroed blocks
                  CALL dbcsr_filter(mat_3c_overl_int_cut(i_cut_RI, i_mem, i_cell, j_cell)%matrix, eps_filter)

               END DO
            END DO

         END DO
      END DO

      ! we need to keep the 3c integrals in case we are doing k-points for GW
      IF (.NOT. (do_gw_im_time .AND. do_kpoints_cubic_RPA)) THEN

         ! to be 100 % safe, set original three center overlap integrals to zero and filter
         DO i_cut_RI = 1, cut_RI
            DO i_cell = 1, num_3c_repl
               DO j_cell = 1, num_3c_repl
                  CALL dbcsr_set(mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix, 0.0_dp)
                  CALL dbcsr_filter(mat_3c_overl_int(i_cut_RI, i_cell, j_cell)%matrix, 1.0_dp)
               END DO
            END DO
         END DO

         CALL dbcsr_deallocate_matrix_set(mat_3c_overl_int)

      END IF

      CALL timestop(handle)

   END SUBROUTINE setup_mat_for_mem_cut_3c

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_occ_global ...
!> \param mat_dm_virt_global ...
!> \param does_mat_P_T_tau_have_blocks ...
! **************************************************************************************************
   SUBROUTINE clean_up(mat_dm_occ_global, mat_dm_virt_global, does_mat_P_T_tau_have_blocks)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ_global, mat_dm_virt_global
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: does_mat_P_T_tau_have_blocks

      CALL dbcsr_deallocate_matrix_set(mat_dm_occ_global)
      CALL dbcsr_deallocate_matrix_set(mat_dm_virt_global)
      DEALLOCATE (does_mat_P_T_tau_have_blocks)

   END SUBROUTINE clean_up

! **************************************************************************************************
!> \brief Calculate kpoint density matrices (rho(k), owned by kpoint groups)
!> \param kpoint    kpoint environment
!> \param tau ...
!> \param e_fermi ...
!> \param stabilize_exp ...
!> \param remove_occ ...
!> \param remove_virt ...
! **************************************************************************************************
   SUBROUTINE kpoint_density_matrices_rpa(kpoint, tau, e_fermi, stabilize_exp, remove_occ, remove_virt)

      TYPE(kpoint_type), POINTER                         :: kpoint
      REAL(KIND=dp), INTENT(IN)                          :: tau, e_fermi, stabilize_exp
      LOGICAL, INTENT(IN)                                :: remove_occ, remove_virt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_density_matrices_rpa', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_mo, ikpgr, ispin, kplocal, &
                                                            nao, nmo, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, exp_scaling, occupation
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: cpmat, fwork, rpmat
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      ! only imaginary wavefunctions supported in kpoint cubic scaling RPA
      CPASSERT(kpoint%use_real_wfn .EQV. .FALSE.)

      ! work matrix
      mo_set => kpoint%kp_env(1)%kpoint_env%mos(1, 1)%mo_set
      CALL get_mo_set(mo_set, nao=nao, nmo=nmo)

      ! if this CPASSERT is triggered, please add all virtual MOs to SCF section,
      ! e.g. ADDED_MOS 1000000
      CPASSERT(nao == nmo)

      ALLOCATE (exp_scaling(nmo))

      CALL cp_fm_get_info(mo_set%mo_coeff, matrix_struct=matrix_struct)
      CALL cp_fm_create(fwork, matrix_struct)

      CALL get_kpoint_info(kpoint, kp_range=kp_range)
      kplocal = kp_range(2) - kp_range(1) + 1

      DO ikpgr = 1, kplocal
         kp => kpoint%kp_env(ikpgr)%kpoint_env
         nspin = SIZE(kp%mos, 2)
         DO ispin = 1, nspin
            mo_set => kp%mos(1, ispin)%mo_set
            CALL get_mo_set(mo_set, eigenvalues=eigenvalues)
            rpmat => kp%wmat(1, ispin)%matrix
            cpmat => kp%wmat(2, ispin)%matrix
            CALL get_mo_set(mo_set, occupation_numbers=occupation)
            CALL cp_fm_to_fm(mo_set%mo_coeff, fwork)

            IF (remove_virt) THEN
               CALL cp_fm_column_scale(fwork, occupation)
            END IF
            IF (remove_occ) THEN
               CALL cp_fm_column_scale(fwork, 2.0_dp/REAL(nspin, KIND=dp) - occupation)
            END IF

            ! proper spin
            IF (nspin == 1) THEN
               CALL cp_fm_scale(0.5_dp, fwork)
            END IF

            DO i_mo = 1, nmo
               IF (ABS(tau*0.5_dp*(eigenvalues(i_mo) - e_fermi)) < stabilize_exp) THEN
                  exp_scaling(i_mo) = EXP(-ABS(tau*(eigenvalues(i_mo) - e_fermi)))
               ELSE
                  exp_scaling(i_mo) = 0.0_dp
               END IF
            END DO

            CALL cp_fm_column_scale(fwork, exp_scaling)

            ! Re(c)*Re(c)
            CALL cp_gemm("N", "T", nao, nao, nmo, 1.0_dp, mo_set%mo_coeff, fwork, 0.0_dp, rpmat)
            mo_set => kp%mos(2, ispin)%mo_set
            ! Im(c)*Re(c)
            CALL cp_gemm("N", "T", nao, nao, nmo, -1.0_dp, mo_set%mo_coeff, fwork, 0.0_dp, cpmat)
            ! Re(c)*Im(c)
            CALL cp_gemm("N", "T", nao, nao, nmo, 1.0_dp, fwork, mo_set%mo_coeff, 1.0_dp, cpmat)

            CALL cp_fm_to_fm(mo_set%mo_coeff, fwork)

            IF (remove_virt) THEN
               CALL cp_fm_column_scale(fwork, occupation)
            END IF
            IF (remove_occ) THEN
               CALL cp_fm_column_scale(fwork, 2.0_dp/REAL(nspin, KIND=dp) - occupation)
            END IF

            ! proper spin
            IF (nspin == 1) THEN
               CALL cp_fm_scale(0.5_dp, fwork)
            END IF

            DO i_mo = 1, nmo
               IF (ABS(tau*0.5_dp*(eigenvalues(i_mo) - e_fermi)) < stabilize_exp) THEN
                  exp_scaling(i_mo) = EXP(-ABS(tau*(eigenvalues(i_mo) - e_fermi)))
               ELSE
                  exp_scaling(i_mo) = 0.0_dp
               END IF
            END DO

            CALL cp_fm_column_scale(fwork, exp_scaling)
            ! Im(c)*Im(c)
            CALL cp_gemm("N", "T", nao, nao, nmo, 1.0_dp, mo_set%mo_coeff, fwork, 1.0_dp, rpmat)

         END DO

      END DO

      CALL cp_fm_release(fwork)
      DEALLOCATE (exp_scaling)

      CALL timestop(handle)

   END SUBROUTINE kpoint_density_matrices_rpa

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_global ...
!> \param qs_env ...
!> \param ispin ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param e_fermi ...
!> \param tau ...
!> \param stabilize_exp ...
!> \param eps_filter ...
!> \param num_cells_dm ...
!> \param index_to_cell_dm ...
!> \param remove_occ ...
!> \param remove_virt ...
!> \param first_jquad ...
! **************************************************************************************************
   SUBROUTINE compute_transl_dm(mat_dm_global, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                stabilize_exp, eps_filter, num_cells_dm, index_to_cell_dm, remove_occ, remove_virt, &
                                first_jquad)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_global
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ispin, num_integ_points, jquad
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, tau, stabilize_exp, eps_filter
      INTEGER, INTENT(OUT)                               :: num_cells_dm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL, INTENT(IN)                                :: remove_occ, remove_virt
      INTEGER, INTENT(IN)                                :: first_jquad

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_transl_dm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_dim, i_img, iquad, jspin, nspin
      INTEGER, DIMENSION(3)                              :: cell_grid_dm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_global_work, matrix_s_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      matrix_s_kp=matrix_s_kp, &
                      sab_orb=sab_nl, &
                      mos=mos, &
                      dft_control=dft_control, &
                      cell=cell, &
                      kpoints=kpoints)

      nspin = SIZE(mos)

      ! we always use an odd number of image cells
      ! CAUTION: also at another point, cell_grid_dm is defined, these definitions have to be identical
      DO i_dim = 1, 3
         cell_grid_dm(i_dim) = (kpoints%nkp_grid(i_dim)/2)*2 - 1
      END DO

      num_cells_dm = cell_grid_dm(1)*cell_grid_dm(2)*cell_grid_dm(3)

      NULLIFY (mat_dm_global_work)
      CALL dbcsr_allocate_matrix_set(mat_dm_global_work, nspin, num_cells_dm)

      DO jspin = 1, nspin

         DO i_img = 1, num_cells_dm

            ALLOCATE (mat_dm_global_work(jspin, i_img)%matrix)
            CALL dbcsr_create(matrix=mat_dm_global_work(jspin, i_img)%matrix, &
                              template=matrix_s_kp(1, 1)%matrix, &
                              !                              matrix_type=dbcsr_type_symmetric)
                              matrix_type=dbcsr_type_no_symmetry)

            CALL dbcsr_reserve_all_blocks(mat_dm_global_work(jspin, i_img)%matrix)

            CALL dbcsr_set(mat_dm_global_work(ispin, i_img)%matrix, 0.0_dp)

         END DO

      END DO

      ! density matrices in k-space weighted with EXP(-|e_i-e_F|*t) for occupied orbitals
      CALL kpoint_density_matrices_rpa(kpoints, tau, e_fermi, stabilize_exp, &
                                       remove_occ=remove_occ, remove_virt=remove_virt)

      ! overwrite the cell indices in kpoints
      CALL init_cell_index_rpa(cell_grid_dm, kpoints%cell_to_index, kpoints%index_to_cell, cell)

      ! density matrices in real space, the cell vectors T for transforming are taken from kpoints%index_to_cell
      ! (custom made for RPA) and not from sab_nl (which is symmetric and from SCF)
      CALL density_matrix_from_kp_to_transl(kpoints, mat_dm_global_work, kpoints%index_to_cell)

      ! we need the index to cell for the density matrices later
      index_to_cell_dm => kpoints%index_to_cell

      ! normally, jquad = 1 to allocate the matrix set, but for GW jquad = 0 is the exchange self-energy
      IF (jquad == first_jquad) THEN

         NULLIFY (mat_dm_global)
!         CALL dbcsr_allocate_matrix_set(mat_dm_global, jquad:num_integ_points, num_cells_dm)
         ALLOCATE (mat_dm_global(first_jquad:num_integ_points, num_cells_dm))

         DO iquad = first_jquad, num_integ_points
            DO i_img = 1, num_cells_dm
               NULLIFY (mat_dm_global(iquad, i_img)%matrix)
               ALLOCATE (mat_dm_global(iquad, i_img)%matrix)
               CALL dbcsr_create(matrix=mat_dm_global(iquad, i_img)%matrix, &
                                 template=matrix_s_kp(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)

            END DO
         END DO

      END IF

      DO i_img = 1, num_cells_dm

         ! filter to get rid of the blocks full with zeros on the lower half, otherwise blocks doubled
         CALL dbcsr_filter(mat_dm_global_work(ispin, i_img)%matrix, eps_filter)

         CALL dbcsr_copy(mat_dm_global(jquad, i_img)%matrix, &
                         mat_dm_global_work(ispin, i_img)%matrix)

      END DO

      CALL dbcsr_deallocate_matrix_set(mat_dm_global_work)

      CALL timestop(handle)

   END SUBROUTINE compute_transl_dm

! **************************************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param mat_dm_global_work ...
!> \param index_to_cell ...
! **************************************************************************************************
   SUBROUTINE density_matrix_from_kp_to_transl(kpoints, mat_dm_global_work, index_to_cell)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_global_work
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: index_to_cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_from_kp_to_transl', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icell, ik, ispin, nkp, nspin, &
                                                            xcell, ycell, zcell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp)                                      :: arg, coskl, sinkl
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cp_fm_type), POINTER                          :: cpmat, rpmat
      TYPE(dbcsr_type), POINTER                          :: mat_work_im, mat_work_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL timeset(routineN, handle)

      NULLIFY (cell_to_index, xkp, wkp)

      NULLIFY (mat_work_re)
      CALL dbcsr_init_p(mat_work_re)
      CALL dbcsr_create(matrix=mat_work_re, &
                        template=mat_dm_global_work(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_work_im)
      CALL dbcsr_init_p(mat_work_im)
      CALL dbcsr_create(matrix=mat_work_im, &
                        template=mat_dm_global_work(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, wkp=wkp, &
                           cell_to_index=cell_to_index)

      nspin = SIZE(mat_dm_global_work, 1)

      CPASSERT(SIZE(mat_dm_global_work, 2) == SIZE(index_to_cell, 2))

      DO ispin = 1, nspin

         DO icell = 1, SIZE(mat_dm_global_work, 2)

            CALL dbcsr_set(mat_dm_global_work(ispin, icell)%matrix, 0.0_dp)

         END DO

      END DO

      DO ispin = 1, nspin

         DO ik = 1, nkp

            kp => kpoints%kp_env(ik)%kpoint_env
            rpmat => kp%wmat(1, ispin)%matrix
            cpmat => kp%wmat(2, ispin)%matrix

            CALL copy_fm_to_dbcsr(rpmat, mat_work_re, keep_sparsity=.FALSE.)
            CALL copy_fm_to_dbcsr(cpmat, mat_work_im, keep_sparsity=.FALSE.)

            DO icell = 1, SIZE(mat_dm_global_work, 2)

               xcell = index_to_cell(1, icell)
               ycell = index_to_cell(2, icell)
               zcell = index_to_cell(3, icell)

               arg = REAL(xcell, dp)*xkp(1, ik) + REAL(ycell, dp)*xkp(2, ik) + REAL(zcell, dp)*xkp(3, ik)
               coskl = wkp(ik)*COS(twopi*arg)
               sinkl = wkp(ik)*SIN(twopi*arg)

               CALL dbcsr_add(mat_dm_global_work(ispin, icell)%matrix, mat_work_re, 1.0_dp, coskl)
               CALL dbcsr_add(mat_dm_global_work(ispin, icell)%matrix, mat_work_im, 1.0_dp, sinkl)

            END DO

         END DO
      END DO

      CALL dbcsr_release_p(mat_work_re)
      CALL dbcsr_release_p(mat_work_im)

      CALL timestop(handle)

   END SUBROUTINE density_matrix_from_kp_to_transl

! **************************************************************************************************
!> \brief ...
!> \param cell_grid ...
!> \param cell_to_index ...
!> \param index_to_cell ...
!> \param cell ...
! **************************************************************************************************
   SUBROUTINE init_cell_index_rpa(cell_grid, cell_to_index, index_to_cell, cell)
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell_grid
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      TYPE(cell_type), POINTER                           :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_cell_index_rpa', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: cell_counter, handle, i_cell, &
                                                            index_min_dist, num_cells, xcell, &
                                                            ycell, zcell
      INTEGER, DIMENSION(3)                              :: itm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_unsorted
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_unsorted
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: abs_cell_vectors
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat

      CALL timeset(routineN, handle)

      CALL get_cell(cell=cell, h=hmat)

      num_cells = cell_grid(1)*cell_grid(2)*cell_grid(3)
      itm(:) = cell_grid(:)/2

      ! check that real space super lattice is a (2n+1)x(2m+1)x(2k+1) super lattice with the unit cell
      ! in the middle
      CPASSERT(cell_grid(1) .NE. itm(1)*2)
      CPASSERT(cell_grid(2) .NE. itm(2)*2)
      CPASSERT(cell_grid(3) .NE. itm(3)*2)

      IF (ASSOCIATED(cell_to_index)) DEALLOCATE (cell_to_index)
      IF (ASSOCIATED(index_to_cell)) DEALLOCATE (index_to_cell)

      ALLOCATE (cell_to_index_unsorted(-itm(1):itm(1), -itm(2):itm(2), -itm(3):itm(3)))
      cell_to_index_unsorted(:, :, :) = 0

      ALLOCATE (index_to_cell_unsorted(3, num_cells))
      index_to_cell_unsorted(:, :) = 0

      ALLOCATE (cell_to_index(-itm(1):itm(1), -itm(2):itm(2), -itm(3):itm(3)))
      cell_to_index(:, :, :) = 0

      ALLOCATE (index_to_cell(3, num_cells))
      index_to_cell(:, :) = 0

      ALLOCATE (abs_cell_vectors(1:num_cells))

      cell_counter = 0

      DO xcell = -itm(1), itm(1)
         DO ycell = -itm(2), itm(2)
            DO zcell = -itm(3), itm(3)

               cell_counter = cell_counter + 1
               cell_to_index_unsorted(xcell, ycell, zcell) = cell_counter

               index_to_cell_unsorted(1, cell_counter) = xcell
               index_to_cell_unsorted(2, cell_counter) = ycell
               index_to_cell_unsorted(3, cell_counter) = zcell

               cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell_unsorted(1:3, cell_counter), dp))

               abs_cell_vectors(cell_counter) = SQRT(cell_vector(1)**2 + cell_vector(2)**2 + cell_vector(3)**2)

            END DO
         END DO
      END DO

      ! first only do all symmetry non-equivalent cells, we need that because chi^T is computed for
      ! cell indices T from index_to_cell(:,1:num_cells/2+1)
      DO i_cell = 1, num_cells/2 + 1

         index_min_dist = MINLOC(abs_cell_vectors(1:num_cells/2 + 1), DIM=1)

         xcell = index_to_cell_unsorted(1, index_min_dist)
         ycell = index_to_cell_unsorted(2, index_min_dist)
         zcell = index_to_cell_unsorted(3, index_min_dist)

         index_to_cell(1, i_cell) = xcell
         index_to_cell(2, i_cell) = ycell
         index_to_cell(3, i_cell) = zcell

         cell_to_index(xcell, ycell, zcell) = i_cell

         abs_cell_vectors(index_min_dist) = 10000000000.0_dp

      END DO

      ! now all the remaining cells
      DO i_cell = num_cells/2 + 2, num_cells

         index_min_dist = MINLOC(abs_cell_vectors(1:num_cells), DIM=1)

         xcell = index_to_cell_unsorted(1, index_min_dist)
         ycell = index_to_cell_unsorted(2, index_min_dist)
         zcell = index_to_cell_unsorted(3, index_min_dist)

         index_to_cell(1, i_cell) = xcell
         index_to_cell(2, i_cell) = ycell
         index_to_cell(3, i_cell) = zcell

         cell_to_index(xcell, ycell, zcell) = i_cell

         abs_cell_vectors(index_min_dist) = 10000000000.0_dp

      END DO

      DEALLOCATE (index_to_cell_unsorted, cell_to_index_unsorted, abs_cell_vectors)

      CALL timestop(handle)

   END SUBROUTINE init_cell_index_rpa

! **************************************************************************************************
!> \brief ...
!> \param i_cell_R ...
!> \param i_cell_S ...
!> \param i_cell_R_minus_S ...
!> \param index_to_cell_3c ...
!> \param cell_to_index_3c ...
!> \param index_to_cell_dm ...
!> \param R_minus_S_needed ...
!> \param do_kpoints_cubic_RPA ...
! **************************************************************************************************
   SUBROUTINE get_diff_index_3c(i_cell_R, i_cell_S, i_cell_R_minus_S, index_to_cell_3c, &
                                cell_to_index_3c, index_to_cell_dm, R_minus_S_needed, &
                                do_kpoints_cubic_RPA)

      INTEGER, INTENT(IN)                                :: i_cell_R, i_cell_S
      INTEGER, INTENT(OUT)                               :: i_cell_R_minus_S
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: cell_to_index_3c
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL, INTENT(OUT)                               :: R_minus_S_needed
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_diff_index_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, x_cell_R, x_cell_R_minus_S, x_cell_S, y_cell_R, y_cell_R_minus_S, &
         y_cell_S, z_cell_R, z_cell_R_minus_S, z_cell_S

      CALL timeset(routineN, handle)

      IF (do_kpoints_cubic_RPA) THEN

         x_cell_R = index_to_cell_3c(1, i_cell_R)
         y_cell_R = index_to_cell_3c(2, i_cell_R)
         z_cell_R = index_to_cell_3c(3, i_cell_R)

         x_cell_S = index_to_cell_dm(1, i_cell_S)
         y_cell_S = index_to_cell_dm(2, i_cell_S)
         z_cell_S = index_to_cell_dm(3, i_cell_S)

         x_cell_R_minus_S = x_cell_R - x_cell_S
         y_cell_R_minus_S = y_cell_R - y_cell_S
         z_cell_R_minus_S = z_cell_R - z_cell_S

         IF (x_cell_R_minus_S .GE. LBOUND(cell_to_index_3c, 1) .AND. &
             x_cell_R_minus_S .LE. UBOUND(cell_to_index_3c, 1) .AND. &
             y_cell_R_minus_S .GE. LBOUND(cell_to_index_3c, 2) .AND. &
             y_cell_R_minus_S .LE. UBOUND(cell_to_index_3c, 2) .AND. &
             z_cell_R_minus_S .GE. LBOUND(cell_to_index_3c, 3) .AND. &
             z_cell_R_minus_S .LE. UBOUND(cell_to_index_3c, 3)) THEN

            i_cell_R_minus_S = cell_to_index_3c(x_cell_R_minus_S, y_cell_R_minus_S, z_cell_R_minus_S)

            ! 0 means that there is no 3c index with this R-S vector because R-S is too big and the 3c integral is 0
            IF (i_cell_R_minus_S == 0) THEN

               R_minus_S_needed = .FALSE.
               i_cell_R_minus_S = 0

            ELSE

               R_minus_S_needed = .TRUE.

            END IF

         ELSE

            i_cell_R_minus_S = 0
            R_minus_S_needed = .FALSE.

         END IF

      ELSE ! no k-points

         R_minus_S_needed = .TRUE.
         i_cell_R_minus_S = 1

      END IF

      CALL timestop(handle)

   END SUBROUTINE get_diff_index_3c

! **************************************************************************************************
!> \brief ...
!> \param i_cell_R ...
!> \param i_cell_S ...
!> \param i_cell_T ...
!> \param i_cell_R_minus_S_minus_T ...
!> \param index_to_cell_3c ...
!> \param cell_to_index_3c ...
!> \param index_to_cell_dm ...
!> \param R_minus_S_minus_T_needed ...
!> \param do_kpoints_cubic_RPA ...
! **************************************************************************************************
   SUBROUTINE get_diff_diff_index_3c(i_cell_R, i_cell_S, i_cell_T, i_cell_R_minus_S_minus_T, &
                                     index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                     R_minus_S_minus_T_needed, &
                                     do_kpoints_cubic_RPA)

      INTEGER, INTENT(IN)                                :: i_cell_R, i_cell_S, i_cell_T
      INTEGER, INTENT(OUT)                               :: i_cell_R_minus_S_minus_T
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: cell_to_index_3c
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL, INTENT(OUT)                               :: R_minus_S_minus_T_needed
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_diff_diff_index_3c', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, x_cell_R, x_cell_R_minus_S_minus_T, x_cell_S, x_cell_T, y_cell_R, &
         y_cell_R_minus_S_minus_T, y_cell_S, y_cell_T, z_cell_R, z_cell_R_minus_S_minus_T, &
         z_cell_S, z_cell_T

      CALL timeset(routineN, handle)

      IF (do_kpoints_cubic_RPA) THEN

         x_cell_R = index_to_cell_3c(1, i_cell_R)
         y_cell_R = index_to_cell_3c(2, i_cell_R)
         z_cell_R = index_to_cell_3c(3, i_cell_R)

         x_cell_S = index_to_cell_dm(1, i_cell_S)
         y_cell_S = index_to_cell_dm(2, i_cell_S)
         z_cell_S = index_to_cell_dm(3, i_cell_S)

         x_cell_T = index_to_cell_dm(1, i_cell_T)
         y_cell_T = index_to_cell_dm(2, i_cell_T)
         z_cell_T = index_to_cell_dm(3, i_cell_T)

         x_cell_R_minus_S_minus_T = x_cell_R - x_cell_S - x_cell_T
         y_cell_R_minus_S_minus_T = y_cell_R - y_cell_S - y_cell_T
         z_cell_R_minus_S_minus_T = z_cell_R - z_cell_S - z_cell_T

         IF (x_cell_R_minus_S_minus_T .GE. LBOUND(cell_to_index_3c, 1) .AND. &
             x_cell_R_minus_S_minus_T .LE. UBOUND(cell_to_index_3c, 1) .AND. &
             y_cell_R_minus_S_minus_T .GE. LBOUND(cell_to_index_3c, 2) .AND. &
             y_cell_R_minus_S_minus_T .LE. UBOUND(cell_to_index_3c, 2) .AND. &
             z_cell_R_minus_S_minus_T .GE. LBOUND(cell_to_index_3c, 3) .AND. &
             z_cell_R_minus_S_minus_T .LE. UBOUND(cell_to_index_3c, 3)) THEN

            i_cell_R_minus_S_minus_T = cell_to_index_3c(x_cell_R_minus_S_minus_T, &
                                                        y_cell_R_minus_S_minus_T, &
                                                        z_cell_R_minus_S_minus_T)

            ! index 0 means that there are only no 3c matrix elements because R-S-T is too big
            IF (i_cell_R_minus_S_minus_T == 0) THEN

               R_minus_S_minus_T_needed = .FALSE.

            ELSE

               R_minus_S_minus_T_needed = .TRUE.

            END IF

         ELSE

            i_cell_R_minus_S_minus_T = 0
            R_minus_S_minus_T_needed = .FALSE.

         END IF

         !  no k-kpoints
      ELSE

         R_minus_S_minus_T_needed = .TRUE.
         i_cell_R_minus_S_minus_T = 1

      END IF

      CALL timestop(handle)

   END SUBROUTINE get_diff_diff_index_3c

! **************************************************************************************************
!> \brief ...
!> \param t_3c ...
!> \param pgrid ...
!> \param nodata ...
! **************************************************************************************************
   SUBROUTINE tensor_change_pgrid(t_3c, pgrid, nodata)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)               :: pgrid
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata

      CHARACTER(default_string_length)                   :: name
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: bs1, bs2, bs3, dist1, dist2, dist3, &
                                                            map1, map2
      INTEGER, DIMENSION(3)                              :: pcoord, pdims, tdims
      TYPE(dbcsr_t_distribution_type)                    :: dist
      TYPE(dbcsr_t_type)                                 :: t_tmp

      CALL dbcsr_t_mp_environ_pgrid(pgrid, pdims, pcoord)
      CALL dbcsr_t_get_info(t_3c, nblks_total=tdims, blk_size_1=bs1, blk_size_2=bs2, blk_size_3=bs3, &
                            name=name)

      ALLOCATE (dist1(tdims(1)))
      ALLOCATE (dist2(tdims(2)))
      ALLOCATE (dist3(tdims(3)))

      CALL cyclic_tensor_dist(tdims(1), pdims(1), bs1, dist1)
      CALL cyclic_tensor_dist(tdims(2), pdims(2), bs2, dist2)
      CALL cyclic_tensor_dist(tdims(3), pdims(3), bs3, dist3)

      CALL dbcsr_t_get_mapping_info(t_3c%nd_index_blk, map1_2d=map1, map2_2d=map2)
      CALL dbcsr_t_distribution_new(dist, pgrid, map1, map2, dist1, dist2, dist3)
      CALL dbcsr_t_create(t_tmp, name, dist, map1, map2, dbcsr_type_real_8, bs1, bs2, bs3)
      CALL dbcsr_t_distribution_destroy(dist)

      IF (PRESENT(nodata)) THEN
         IF (.NOT. nodata) CALL dbcsr_t_copy(t_3c, t_tmp, move_data=.TRUE.)
      ELSE
         CALL dbcsr_t_copy(t_3c, t_tmp, move_data=.TRUE.)
      ENDIF

      CALL dbcsr_t_destroy(t_3c)
      t_3c = t_tmp
   END SUBROUTINE

END MODULE rpa_im_time
