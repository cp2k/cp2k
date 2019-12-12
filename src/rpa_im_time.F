!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
        dbcsr_add, dbcsr_clear, dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_get_num_blocks, &
        dbcsr_get_occupation, dbcsr_init_p, dbcsr_p_type, dbcsr_release_p, &
        dbcsr_reserve_all_blocks, dbcsr_scalar, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_clear, dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, &
        dbcsr_t_copy_tensor_to_matrix, dbcsr_t_create, dbcsr_t_destroy, dbcsr_t_filter, &
        dbcsr_t_get_info, dbcsr_t_need_contract, dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, &
        dbcsr_t_type, dbcsr_t_write_split_info
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_sync,&
                                              mp_waitall
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_tensors,                      ONLY: tensor_change_pgrid
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time'

   PUBLIC :: compute_mat_P_omega, &
             compute_transl_dm, &
             init_cell_index_rpa, &
             zero_mat_P_omega

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
!> \param ispin ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
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
   SUBROUTINE compute_mat_P_omega(mat_P_omega, fm_scaled_dm_occ_tau, &
                                  fm_scaled_dm_virt_tau, fm_mo_coeff_occ, fm_mo_coeff_virt, &
                                  fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, &
                                  mat_P_global, &
                                  matrix_s, &
                                  ispin, &
                                  t_3c_M, t_3c_O, &
                                  starts_array_mc, ends_array_mc, &
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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: t_3c_M
      TYPE(dbcsr_t_type), DIMENSION(:, :), INTENT(INOUT) :: t_3c_O
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_P_omega', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, handle5, handle6, i, i_cell, i_cell_R_1, i_cell_R_1_minus_S, &
         i_cell_R_1_minus_T, i_cell_R_2, i_cell_R_2_minus_S_minus_T, i_cell_S, i_cell_T, i_mem, &
         iquad, j, j_mem, jquad, num_3c_repl, num_cells_dm, unit_nr_2, unit_nr_prv
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
      IF (memory_info) THEN
         unit_nr_prv = unit_nr
      ELSE
         unit_nr_prv = 0
      ENDIF

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

         unit_nr_2 = unit_nr_prv
         IF (pgrid_1_init_occ) THEN
            DO i = 1, SIZE(t_3c_O, 1)
               DO j = 1, SIZE(t_3c_O, 2)
                  CALL tensor_change_pgrid(t_3c_O_occ(i, j), pgrid_1_use_occ, unit_nr=unit_nr_2)
                  unit_nr_2 = 0
               ENDDO
            ENDDO
         ENDIF

         unit_nr_2 = unit_nr_prv
         IF (pgrid_1_init_virt) THEN
            DO i = 1, SIZE(t_3c_O, 1)
               DO j = 1, SIZE(t_3c_O, 2)
                  CALL tensor_change_pgrid(t_3c_O_virt(i, j), pgrid_1_use_virt, unit_nr=unit_nr_2)
                  unit_nr_2 = 0
               ENDDO
            ENDDO
         ENDIF

         IF (pgrid_2_init) THEN
            CALL tensor_change_pgrid(t_3c_M, pgrid_2_use, nodata=.TRUE., unit_nr=unit_nr_prv)
         ENDIF

         ! loop over T for chi^T(it)

         CALL mp_sync(para_env%group)
         t1 = m_walltime()

         CALL compute_mat_dm_global(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, nmo, &
                                    fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                    fm_mo_coeff_virt_scaled, mat_dm_occ_global, mat_dm_virt_global, &
                                    matrix_s, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, &
                                    unit_nr, para_env, &
                                    jquad, stabilize_exp, do_kpoints_cubic_RPA, qs_env, &
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
               jbounds_1(:, 2) = [starts_array_mc(j_mem), ends_array_mc(j_mem)]

               jbounds_2(:, 1) = [starts_array_mc(j_mem), ends_array_mc(j_mem)]

               DO i_mem = 1, cut_memory

                  ibounds_1(:, 1) = [1, bounds_3c(1)]
                  ibounds_1(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

                  ibounds_2(:, 1) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

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

   END SUBROUTINE compute_mat_P_omega

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
!> \param ispin ...
!> \param Eigenval ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param memory_info ...
!> \param unit_nr ...
!> \param para_env ...
!> \param jquad ...
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
                                    matrix_s, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, &
                                    unit_nr, para_env, &
                                    jquad, stabilize_exp, do_kpoints_cubic_RPA, qs_env, &
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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, eps_filter
      LOGICAL, INTENT(IN)                                :: memory_info
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: jquad
      REAL(KIND=dp), INTENT(IN)                          :: stabilize_exp
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: num_cells_dm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL, ALLOCATABLE, DIMENSION(:), INTENT(INOUT)  :: does_mat_P_T_tau_have_blocks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_dm_global', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, iquad, jjB, &
                                                            ncol_local, nrow_local, size_dm_occ, &
                                                            size_dm_virt
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: tau

      CALL timeset(routineN, handle)

      IF (memory_info .AND. unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RPA_IM_TIME_INFO| Started with time point: ", jquad

      tau = tau_tj(jquad)

      IF (do_kpoints_cubic_RPA) THEN

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
               CALL dbcsr_create(matrix=mat_dm_occ_global(iquad, 1)%matrix, &
                                 template=matrix_s(1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)

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

END MODULE rpa_im_time
