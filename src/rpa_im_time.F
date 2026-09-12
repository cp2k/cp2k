!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for low-scaling RPA/GW with imaginary time
!> \par History
!>      10.2015 created [Jan Wilhelm]
! **************************************************************************************************
MODULE rpa_im_time
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_scale
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_clear, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_type, dbcsr_filter, dbcsr_get_info, dbcsr_init_p, dbcsr_p_type, &
        dbcsr_release_p, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE dbt_api,                         ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_contract, dbt_copy, &
        dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_create, dbt_destroy, dbt_filter, &
        dbt_get_info, dbt_nblks_total, dbt_nd_mp_comm, dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE hfx_types,                       ONLY: block_ind_type,&
                                              hfx_compression_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_tensors,                      ONLY: decompress_tensor,&
                                              get_tensor_occupancy
   USE qs_tensors_types,                ONLY: create_2c_tensor
   USE rpa_gw_im_time_util,             ONLY: compute_weight_re_im,&
                                              get_atom_index_from_basis_function_index
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time'

   INTEGER, PARAMETER, PUBLIC :: propagator_sector_occupied = 1
   INTEGER, PARAMETER, PUBLIC :: propagator_sector_virtual = 2
   INTEGER, PARAMETER, PUBLIC :: propagator_number_of_sectors = 2

   PUBLIC :: compute_mat_P_omega, &
             compute_transl_dm, &
             init_cell_index_rpa, &
             zero_mat_P_omega, &
             compute_periodic_dm, &
             compute_mat_dm_global, &
             compute_gamma_propagator, &
             create_propagator_matrix_set

CONTAINS

! **************************************************************************************************
!> \brief Creates the sector-, time-, and cell-resolved propagator matrix set.
!> \param propagator propagator matrix set to create
!> \param ntime number of time or frequency points
!> \param matrix_template DBCSR matrix template for every matrix
!> \param index_to_cell cell coordinates for the third matrix-set index
! **************************************************************************************************
   SUBROUTINE create_propagator_matrix_set(propagator, ntime, matrix_template, index_to_cell)
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: propagator
      INTEGER, INTENT(IN)                                :: ntime
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_template
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell

      INTEGER                                            :: icell, isector, jtime, ncell

      CPASSERT(ntime > 0)
      CPASSERT(SIZE(index_to_cell, 1) == 3)
      ncell = SIZE(index_to_cell, 2)
      CPASSERT(ncell > 0)

      CALL dbcsr_allocate_matrix_set(propagator, propagator_number_of_sectors, ntime, ncell)

      DO icell = 1, ncell
         DO jtime = 1, ntime
            DO isector = 1, propagator_number_of_sectors
               ALLOCATE (propagator(isector, jtime, icell)%matrix)
               CALL dbcsr_create(matrix=propagator(isector, jtime, icell)%matrix, &
                                 template=matrix_template, &
                                 matrix_type=dbcsr_type_no_symmetry)
            END DO
         END DO
      END DO

   END SUBROUTINE create_propagator_matrix_set

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
!> \param cfm_mo_coeff ...
!> \param homo ...
!> \param mat_P_global ...
!> \param matrix_s ...
!> \param ispin ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param starts_array_mc_block ...
!> \param ends_array_mc_block ...
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
!> \param qs_env ...
!> \param do_kpoints_from_Gamma ...
!> \param index_to_cell_3c ...
!> \param cell_to_index_3c ...
!> \param has_mat_P_blocks ...
!> \param do_ri_sos_laplace_mp2 ...
!> \param dbcsr_time ...
!> \param dbcsr_nflop ...
! **************************************************************************************************
   SUBROUTINE compute_mat_P_omega(mat_P_omega, cfm_mo_coeff, homo, &
                                  mat_P_global, &
                                  matrix_s, &
                                  ispin, &
                                  t_3c_M, t_3c_O, t_3c_O_compressed, t_3c_O_ind, &
                                  starts_array_mc, ends_array_mc, &
                                  starts_array_mc_block, ends_array_mc_block, &
                                  weights_cos_tf_t_to_w, &
                                  tj, tau_tj, e_fermi, eps_filter, &
                                  alpha, eps_filter_im_time, Eigenval, nmo, &
                                  num_integ_points, cut_memory, unit_nr, &
                                  mp2_env, para_env, &
                                  qs_env, do_kpoints_from_Gamma, &
                                  index_to_cell_3c, cell_to_index_3c, &
                                  has_mat_P_blocks, do_ri_sos_laplace_mp2, &
                                  dbcsr_time, dbcsr_nflop)
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN)    :: mat_P_omega
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mo_coeff
      INTEGER, INTENT(IN)                                :: homo
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_M
      TYPE(dbt_type), DIMENSION(:, :), INTENT(INOUT)     :: t_3c_O
      TYPE(hfx_compression_type), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: t_3c_O_compressed
      TYPE(block_ind_type), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: t_3c_O_ind
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_t_to_w
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tau_tj
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, eps_filter, alpha, &
                                                            eps_filter_im_time
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: nmo, num_integ_points, cut_memory, &
                                                            unit_nr
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: do_kpoints_from_Gamma
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: cell_to_index_3c
      LOGICAL, DIMENSION(:, :, :, :, :), INTENT(INOUT)   :: has_mat_P_blocks
      LOGICAL, INTENT(IN)                                :: do_ri_sos_laplace_mp2
      REAL(dp), INTENT(INOUT)                            :: dbcsr_time
      INTEGER(int_8), INTENT(INOUT)                      :: dbcsr_nflop

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_P_omega'

      INTEGER :: comm_2d_handle, handle, handle2, handle3, i, i_cell, i_cell_R_1, &
         i_cell_R_1_minus_S, i_cell_R_1_minus_T, i_cell_R_2, i_cell_R_2_minus_S_minus_T, i_cell_S, &
         i_cell_T, i_mem, iquad, j, j_mem, jquad, num_3c_repl, num_cells_dm, unit_nr_dbcsr
      INTEGER(int_8)                                     :: nze, nze_dm_occ, nze_dm_virt, nze_M_occ, &
                                                            nze_M_virt, nze_O
      INTEGER(KIND=int_8)                                :: flops_1_occ, flops_1_virt, flops_2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2, mc_ranges, size_dm, &
                                                            size_P
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(2, 1)                           :: ibounds_2, jbounds_2
      INTEGER, DIMENSION(2, 2)                           :: ibounds_1, jbounds_1
      INTEGER, DIMENSION(3)                              :: bounds_3c
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL :: do_Gamma_RPA, do_kpoints_cubic_RPA, first_cycle_im_time, first_cycle_omega_loop, &
         memory_info, R_1_minus_S_needed, R_1_minus_T_needed, R_2_minus_S_minus_T_needed
      REAL(dp)                                           :: occ, occ_dm_occ, occ_dm_virt, occ_M_occ, &
                                                            occ_M_virt, occ_O, t1_flop
      REAL(KIND=dp)                                      :: omega, omega_old, t1, t2, tau, weight, &
                                                            weight_old
      TYPE(dbcsr_distribution_type)                      :: dist_P
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: propagator
      TYPE(dbt_pgrid_type)                               :: pgrid_2d
      TYPE(dbt_type)                                     :: t_3c_M_occ, t_3c_M_occ_tmp, t_3c_M_virt, &
                                                            t_3c_M_virt_tmp, t_dm, t_dm_tmp, t_P, &
                                                            t_P_tmp
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_dm_occ, t_dm_virt
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_O_occ, t_3c_O_virt
      TYPE(mp_comm_type)                                 :: comm_2d

      CALL timeset(routineN, handle)

      NULLIFY (propagator)

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_dbcsr = unit_nr
      ELSE
         unit_nr_dbcsr = 0
      END IF

      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints
      do_Gamma_RPA = .NOT. do_kpoints_cubic_RPA
      num_3c_repl = MAXVAL(cell_to_index_3c)

      first_cycle_im_time = .TRUE.
      ALLOCATE (t_3c_O_occ(SIZE(t_3c_O, 1), SIZE(t_3c_O, 2)), t_3c_O_virt(SIZE(t_3c_O, 1), SIZE(t_3c_O, 2)))
      DO i = 1, SIZE(t_3c_O, 1)
         DO j = 1, SIZE(t_3c_O, 2)
            CALL dbt_create(t_3c_O(i, j), t_3c_O_occ(i, j))
            CALL dbt_create(t_3c_O(i, j), t_3c_O_virt(i, j))
         END DO
      END DO

      CALL dbt_create(t_3c_M, t_3c_M_occ, name="M occ (RI | AO AO)")
      CALL dbt_create(t_3c_M, t_3c_M_virt, name="M virt (RI | AO AO)")

      ALLOCATE (mc_ranges(cut_memory + 1))
      mc_ranges(:cut_memory) = starts_array_mc_block(:)
      mc_ranges(cut_memory + 1) = ends_array_mc_block(cut_memory) + 1

      DO jquad = 1, num_integ_points

         CALL para_env%sync()
         t1 = m_walltime()

         CALL compute_mat_dm_global(tau_tj, num_integ_points, nmo, cfm_mo_coeff, homo, propagator, &
                                    matrix_s, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, unit_nr, &
                                    jquad, do_kpoints_cubic_RPA, do_kpoints_from_Gamma, qs_env, &
                                    num_cells_dm, index_to_cell_dm, para_env)

         ALLOCATE (t_dm_virt(num_cells_dm))
         ALLOCATE (t_dm_occ(num_cells_dm))
         CALL dbcsr_get_info(mat_P_global%matrix, distribution=dist_P)
         CALL dbcsr_distribution_get(dist_P, group=comm_2d_handle, nprows=pdims_2d(1), npcols=pdims_2d(2))
         CALL comm_2d%set_handle(comm_2d_handle)

         pgrid_2d = dbt_nd_mp_comm(comm_2d, [1], [2], pdims_2d=pdims_2d)
         ALLOCATE (size_P(dbt_nblks_total(t_3c_M, 1)))
         CALL dbt_get_info(t_3c_M, blk_size_1=size_P)

         ALLOCATE (size_dm(dbt_nblks_total(t_3c_O(1, 1), 3)))
         CALL dbt_get_info(t_3c_O(1, 1), blk_size_3=size_dm)
         CALL create_2c_tensor(t_dm, dist_1, dist_2, pgrid_2d, size_dm, size_dm, name="D (AO | AO)")
         DEALLOCATE (size_dm)
         DEALLOCATE (dist_1, dist_2)
         CALL create_2c_tensor(t_P, dist_1, dist_2, pgrid_2d, size_P, size_P, name="P (RI | RI)")
         DEALLOCATE (size_P)
         DEALLOCATE (dist_1, dist_2)
         CALL dbt_pgrid_destroy(pgrid_2d)

         DO i_cell = 1, num_cells_dm
            CALL dbt_create(t_dm, t_dm_virt(i_cell), name="D virt (AO | AO)")
            CALL dbt_create(propagator(propagator_sector_virtual, jquad, i_cell)%matrix, t_dm_tmp)
            CALL dbt_copy_matrix_to_tensor(propagator(propagator_sector_virtual, jquad, i_cell)%matrix, t_dm_tmp)
            CALL dbt_copy(t_dm_tmp, t_dm_virt(i_cell), move_data=.TRUE.)
            CALL dbcsr_clear(propagator(propagator_sector_virtual, jquad, i_cell)%matrix)

            CALL dbt_create(t_dm, t_dm_occ(i_cell), name="D occ (AO | AO)")
            CALL dbt_copy_matrix_to_tensor(propagator(propagator_sector_occupied, jquad, i_cell)%matrix, t_dm_tmp)
            CALL dbt_copy(t_dm_tmp, t_dm_occ(i_cell), move_data=.TRUE.)
            CALL dbt_destroy(t_dm_tmp)
            CALL dbcsr_clear(propagator(propagator_sector_occupied, jquad, i_cell)%matrix)
         END DO

         CALL get_tensor_occupancy(t_dm_occ(1), nze_dm_occ, occ_dm_occ)
         CALL get_tensor_occupancy(t_dm_virt(1), nze_dm_virt, occ_dm_virt)

         CALL dbt_destroy(t_dm)

         CALL dbt_create(t_3c_O_occ(1, 1), t_3c_M_occ_tmp, name="M (RI AO | AO)")
         CALL dbt_create(t_3c_O_virt(1, 1), t_3c_M_virt_tmp, name="M (RI AO | AO)")

         CALL timeset(routineN//"_contract", handle2)

         CALL para_env%sync()
         t1_flop = m_walltime()

         DO i = 1, SIZE(t_3c_O_occ, 1)
            DO j = 1, SIZE(t_3c_O_occ, 2)
               CALL dbt_batched_contract_init(t_3c_O_occ(i, j), batch_range_2=mc_ranges, batch_range_3=mc_ranges)
            END DO
         END DO
         DO i = 1, SIZE(t_3c_O_virt, 1)
            DO j = 1, SIZE(t_3c_O_virt, 2)
               CALL dbt_batched_contract_init(t_3c_O_virt(i, j), batch_range_2=mc_ranges, batch_range_3=mc_ranges)
            END DO
         END DO
         CALL dbt_batched_contract_init(t_3c_M_occ_tmp, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
         CALL dbt_batched_contract_init(t_3c_M_virt_tmp, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
         CALL dbt_batched_contract_init(t_3c_M_occ, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
         CALL dbt_batched_contract_init(t_3c_M_virt, batch_range_2=mc_ranges, batch_range_3=mc_ranges)

         DO i_cell_T = 1, num_cells_dm/2 + 1

            IF (.NOT. ANY(has_mat_P_blocks(i_cell_T, :, :, :, :))) CYCLE

            CALL dbt_batched_contract_init(t_P)

            IF (do_Gamma_RPA) THEN
               nze_O = 0
               nze_M_virt = 0
               nze_M_occ = 0
               occ_M_virt = 0.0_dp
               occ_M_occ = 0.0_dp
               occ_O = 0.0_dp
            END IF

            DO j_mem = 1, cut_memory

               CALL dbt_get_info(t_3c_O_occ(1, 1), nfull_total=bounds_3c)

               jbounds_1(:, 1) = [1, bounds_3c(1)]
               jbounds_1(:, 2) = [starts_array_mc(j_mem), ends_array_mc(j_mem)]

               jbounds_2(:, 1) = [starts_array_mc(j_mem), ends_array_mc(j_mem)]

               IF (do_Gamma_RPA) CALL dbt_batched_contract_init(t_dm_virt(1))

               DO i_mem = 1, cut_memory

                  IF (.NOT. ANY(has_mat_P_blocks(i_cell_T, i_mem, j_mem, :, :))) CYCLE

                  ibounds_1(:, 1) = [1, bounds_3c(1)]
                  ibounds_1(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

                  ibounds_2(:, 1) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

                  IF (unit_nr_dbcsr > 0) WRITE (UNIT=unit_nr_dbcsr, FMT="(T3,A,I3,1X,I3)") &
                     "RPA_LOW_SCALING_INFO| Memory Cut iteration", i_mem, j_mem

                  DO i_cell_R_1 = 1, num_3c_repl

                     DO i_cell_R_2 = 1, num_3c_repl

                        IF (.NOT. has_mat_P_blocks(i_cell_T, i_mem, j_mem, i_cell_R_1, i_cell_R_2)) CYCLE

                        CALL get_diff_index_3c(i_cell_R_1, i_cell_T, i_cell_R_1_minus_T, &
                                               index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                               R_1_minus_T_needed, do_kpoints_cubic_RPA)

                        IF (do_Gamma_RPA) CALL dbt_batched_contract_init(t_dm_occ(1))
                        DO i_cell_S = 1, num_cells_dm
                           CALL get_diff_index_3c(i_cell_R_1, i_cell_S, i_cell_R_1_minus_S, index_to_cell_3c, &
                                                  cell_to_index_3c, index_to_cell_dm, R_1_minus_S_needed, &
                                                  do_kpoints_cubic_RPA)
                           IF (R_1_minus_S_needed) THEN

                              CALL timeset(routineN//"_calc_M_occ_t", handle3)
                              CALL decompress_tensor(t_3c_O(i_cell_R_1_minus_S, i_cell_R_2), &
                                                     t_3c_O_ind(i_cell_R_1_minus_S, i_cell_R_2, j_mem)%ind, &
                                                     t_3c_O_compressed(i_cell_R_1_minus_S, i_cell_R_2, j_mem), &
                                                     qs_env%mp2_env%ri_rpa_im_time%eps_compress)

                              IF (do_Gamma_RPA .AND. i_mem == 1) THEN
                                 CALL get_tensor_occupancy(t_3c_O(1, 1), nze, occ)
                                 nze_O = nze_O + nze
                                 occ_O = occ_O + occ
                              END IF

                              CALL dbt_copy(t_3c_O(i_cell_R_1_minus_S, i_cell_R_2), &
                                            t_3c_O_occ(i_cell_R_1_minus_S, i_cell_R_2), move_data=.TRUE.)

                              CALL dbt_contract(alpha=1.0_dp, &
                                                tensor_1=t_3c_O_occ(i_cell_R_1_minus_S, i_cell_R_2), &
                                                tensor_2=t_dm_occ(i_cell_S), &
                                                beta=1.0_dp, &
                                                tensor_3=t_3c_M_occ_tmp, &
                                                contract_1=[3], notcontract_1=[1, 2], &
                                                contract_2=[2], notcontract_2=[1], &
                                                map_1=[1, 2], map_2=[3], &
                                                bounds_2=jbounds_1, bounds_3=ibounds_2, &
                                                filter_eps=eps_filter, unit_nr=unit_nr_dbcsr, &
                                                flop=flops_1_occ)
                              CALL timestop(handle3)

                              dbcsr_nflop = dbcsr_nflop + flops_1_occ

                           END IF
                        END DO

                        IF (do_Gamma_RPA) CALL dbt_batched_contract_finalize(t_dm_occ(1))

                        ! copy matrix to optimal contraction layout - copy is done manually in order
                        ! to better control memory allocations (we can release data of previous
                        ! representation)
                        CALL timeset(routineN//"_copy_M_occ_t", handle3)
                        CALL dbt_copy(t_3c_M_occ_tmp, t_3c_M_occ, order=[1, 3, 2], move_data=.TRUE.)
                        CALL dbt_filter(t_3c_M_occ, eps_filter)
                        CALL timestop(handle3)

                        IF (do_Gamma_RPA) THEN
                           CALL get_tensor_occupancy(t_3c_M_occ, nze, occ)
                           nze_M_occ = nze_M_occ + nze
                           occ_M_occ = occ_M_occ + occ
                        END IF

                        DO i_cell_S = 1, num_cells_dm
                           CALL get_diff_diff_index_3c(i_cell_R_2, i_cell_S, i_cell_T, i_cell_R_2_minus_S_minus_T, &
                                                       index_to_cell_3c, cell_to_index_3c, index_to_cell_dm, &
                                                       R_2_minus_S_minus_T_needed, do_kpoints_cubic_RPA)

                           IF (R_1_minus_T_needed .AND. R_2_minus_S_minus_T_needed) THEN
                              CALL decompress_tensor(t_3c_O(i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T), &
                                                     t_3c_O_ind(i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T, i_mem)%ind, &
                                                     t_3c_O_compressed(i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T, i_mem), &
                                                     qs_env%mp2_env%ri_rpa_im_time%eps_compress)

                              CALL dbt_copy(t_3c_O(i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T), &
                                            t_3c_O_virt(i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T), move_data=.TRUE.)

                              CALL timeset(routineN//"_calc_M_virt_t", handle3)
                              CALL dbt_contract(alpha=alpha/2.0_dp, &
                                                tensor_1=t_3c_O_virt( &
                                                i_cell_R_2_minus_S_minus_T, i_cell_R_1_minus_T), &
                                                tensor_2=t_dm_virt(i_cell_S), &
                                                beta=1.0_dp, &
                                                tensor_3=t_3c_M_virt_tmp, &
                                                contract_1=[3], notcontract_1=[1, 2], &
                                                contract_2=[2], notcontract_2=[1], &
                                                map_1=[1, 2], map_2=[3], &
                                                bounds_2=ibounds_1, bounds_3=jbounds_2, &
                                                filter_eps=eps_filter, unit_nr=unit_nr_dbcsr, &
                                                flop=flops_1_virt)
                              CALL timestop(handle3)

                              dbcsr_nflop = dbcsr_nflop + flops_1_virt

                           END IF
                        END DO

                        CALL timeset(routineN//"_copy_M_virt_t", handle3)
                        CALL dbt_copy(t_3c_M_virt_tmp, t_3c_M_virt, move_data=.TRUE.)
                        CALL dbt_filter(t_3c_M_virt, eps_filter)
                        CALL timestop(handle3)

                        IF (do_Gamma_RPA) THEN
                           CALL get_tensor_occupancy(t_3c_M_virt, nze, occ)
                           nze_M_virt = nze_M_virt + nze
                           occ_M_virt = occ_M_virt + occ
                        END IF

                        flops_2 = 0

                        CALL timeset(routineN//"_calc_P_t", handle3)

                        CALL dbt_contract(alpha=1.0_dp, tensor_1=t_3c_M_occ, &
                                          tensor_2=t_3c_M_virt, &
                                          beta=1.0_dp, &
                                          tensor_3=t_P, &
                                          contract_1=[2, 3], notcontract_1=[1], &
                                          contract_2=[2, 3], notcontract_2=[1], &
                                          map_1=[1], map_2=[2], &
                                          filter_eps=eps_filter_im_time/REAL(cut_memory**2, KIND=dp), &
                                          flop=flops_2, &
                                          move_data=.TRUE., &
                                          unit_nr=unit_nr_dbcsr)

                        CALL timestop(handle3)

                        first_cycle_im_time = .FALSE.

                        IF (jquad == 1 .AND. flops_2 == 0) THEN
                           has_mat_P_blocks(i_cell_T, i_mem, j_mem, i_cell_R_1, i_cell_R_2) = .FALSE.
                        END IF

                     END DO
                  END DO
               END DO
               IF (do_Gamma_RPA) CALL dbt_batched_contract_finalize(t_dm_virt(1))
            END DO

            CALL dbt_batched_contract_finalize(t_P, unit_nr=unit_nr_dbcsr)

            CALL dbt_create(mat_P_global%matrix, t_P_tmp)
            CALL dbt_copy(t_P, t_P_tmp, move_data=.TRUE.)
            CALL dbt_copy_tensor_to_matrix(t_P_tmp, mat_P_global%matrix)
            CALL dbt_destroy(t_P_tmp)

            IF (do_ri_sos_laplace_mp2) THEN
               ! For RI-SOS-Laplace-MP2 we do not perform a cosine transform,
               ! but we have to copy P_local to the output matrix

               CALL dbcsr_add(mat_P_omega(jquad, i_cell_T)%matrix, mat_P_global%matrix, 1.0_dp, 1.0_dp)
            ELSE
               CALL timeset(routineN//"_Fourier_transform", handle3)

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

               CALL timestop(handle3)
            END IF

         END DO

         CALL timestop(handle2)

         CALL dbt_batched_contract_finalize(t_3c_M_occ_tmp)
         CALL dbt_batched_contract_finalize(t_3c_M_virt_tmp)
         CALL dbt_batched_contract_finalize(t_3c_M_occ)
         CALL dbt_batched_contract_finalize(t_3c_M_virt)

         DO i = 1, SIZE(t_3c_O_occ, 1)
            DO j = 1, SIZE(t_3c_O_occ, 2)
               CALL dbt_batched_contract_finalize(t_3c_O_occ(i, j))
            END DO
         END DO

         DO i = 1, SIZE(t_3c_O_virt, 1)
            DO j = 1, SIZE(t_3c_O_virt, 2)
               CALL dbt_batched_contract_finalize(t_3c_O_virt(i, j))
            END DO
         END DO

         CALL dbt_destroy(t_P)
         DO i_cell = 1, num_cells_dm
            CALL dbt_destroy(t_dm_virt(i_cell))
            CALL dbt_destroy(t_dm_occ(i_cell))
         END DO

         CALL dbt_destroy(t_3c_M_occ_tmp)
         CALL dbt_destroy(t_3c_M_virt_tmp)
         DEALLOCATE (t_dm_virt)
         DEALLOCATE (t_dm_occ)

         CALL para_env%sync()
         t2 = m_walltime()

         dbcsr_time = dbcsr_time + t2 - t1_flop

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/T3,A,1X,I3)') &
               'RPA_LOW_SCALING_INFO| Info for time point', jquad
            WRITE (unit_nr, '(T6,A,T56,F25.1)') &
               'Execution time (s):', t2 - t1
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of D occ:', REAL(nze_dm_occ, dp), '/', occ_dm_occ*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of D virt:', REAL(nze_dm_virt, dp), '/', occ_dm_virt*100, '%'
            IF (do_Gamma_RPA) THEN
               WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
                  'Occupancy of 3c ints:', REAL(nze_O, dp), '/', occ_O*100, '%'
               WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
                  'Occupancy of M occ:', REAL(nze_M_occ, dp), '/', occ_M_occ*100, '%'
               WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
                  'Occupancy of M virt:', REAL(nze_M_virt, dp), '/', occ_M_virt*100, '%'
            END IF
            WRITE (unit_nr, *)
            CALL m_flush(unit_nr)
         END IF

      END DO ! time points

      CALL dbt_destroy(t_3c_M_occ)
      CALL dbt_destroy(t_3c_M_virt)

      DO i = 1, SIZE(t_3c_O, 1)
         DO j = 1, SIZE(t_3c_O, 2)
            CALL dbt_destroy(t_3c_O_occ(i, j))
            CALL dbt_destroy(t_3c_O_virt(i, j))
         END DO
      END DO

      CALL dbcsr_deallocate_matrix_set(propagator)

      CALL timestop(handle)

   END SUBROUTINE compute_mat_P_omega

! **************************************************************************************************
!> \brief ...
!> \param mat_P_omega ...
! **************************************************************************************************
   SUBROUTINE zero_mat_P_omega(mat_P_omega)
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN)    :: mat_P_omega

      INTEGER                                            :: i_kp, jquad

      DO jquad = 1, SIZE(mat_P_omega, 1)
         DO i_kp = 1, SIZE(mat_P_omega, 2)

            CALL dbcsr_set(mat_P_omega(jquad, i_kp)%matrix, 0.0_dp)

         END DO
      END DO

   END SUBROUTINE zero_mat_P_omega

! **************************************************************************************************
!> \brief ...
!> \param tau_tj ...
!> \param num_integ_points ...
!> \param nmo ...
!> \param cfm_mo_coeff ...
!> \param homo ...
!> \param propagator ...
!> \param matrix_s ...
!> \param ispin ...
!> \param Eigenval ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param memory_info ...
!> \param unit_nr ...
!> \param jquad ...
!> \param do_kpoints_cubic_RPA ...
!> \param do_kpoints_from_Gamma ...
!> \param qs_env ...
!> \param num_cells_dm ...
!> \param index_to_cell_dm ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_mat_dm_global(tau_tj, num_integ_points, nmo, cfm_mo_coeff, homo, propagator, &
                                    matrix_s, ispin, &
                                    Eigenval, e_fermi, eps_filter, memory_info, unit_nr, &
                                    jquad, do_kpoints_cubic_RPA, do_kpoints_from_Gamma, qs_env, &
                                    num_cells_dm, index_to_cell_dm, para_env)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: tau_tj
      INTEGER, INTENT(IN)                                :: num_integ_points, nmo
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mo_coeff
      INTEGER, INTENT(IN)                                :: homo
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), &
         INTENT(INOUT), POINTER                          :: propagator
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: matrix_s
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, eps_filter
      LOGICAL, INTENT(IN)                                :: memory_info
      INTEGER, INTENT(IN)                                :: unit_nr, jquad
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA, &
                                                            do_kpoints_from_Gamma
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(OUT)                               :: num_cells_dm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_mat_dm_global'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(3, 1)                           :: index_to_cell_zero
      REAL(KIND=dp)                                      :: tau

      index_to_cell_zero = 0

      CALL timeset(routineN, handle)

      IF (memory_info .AND. unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "RPA_LOW_SCALING_INFO| Started with time point: ", jquad

      tau = tau_tj(jquad)

      IF (do_kpoints_cubic_RPA) THEN

         CALL compute_transl_dm(propagator, qs_env, &
                                ispin, num_integ_points, jquad, e_fermi, tau, &
                                eps_filter, num_cells_dm, index_to_cell_dm, &
                                sector=propagator_sector_occupied)

         CALL compute_transl_dm(propagator, qs_env, &
                                ispin, num_integ_points, jquad, e_fermi, tau, &
                                eps_filter, num_cells_dm, index_to_cell_dm, &
                                sector=propagator_sector_virtual)

      ELSE IF (do_kpoints_from_Gamma) THEN

         CALL compute_periodic_dm(propagator, qs_env, &
                                  ispin, num_integ_points, jquad, e_fermi, tau, &
                                  sector=propagator_sector_occupied)

         CALL compute_periodic_dm(propagator, qs_env, &
                                  ispin, num_integ_points, jquad, e_fermi, tau, &
                                  sector=propagator_sector_virtual)

         num_cells_dm = 1

      ELSE

         num_cells_dm = 1

         IF (jquad == 1) THEN
            CALL create_propagator_matrix_set(propagator, num_integ_points, matrix_s(1)%matrix, &
                                              index_to_cell_zero)
         END IF

         CALL compute_gamma_propagator(propagator, jquad, cfm_mo_coeff, homo, Eigenval, nmo, &
                                       eps_filter, e_fermi, tau, para_env)

         ! release memory
         IF (jquad > 1) THEN
            CALL dbcsr_set(propagator(propagator_sector_occupied, jquad - 1, 1)%matrix, 0.0_dp)
            CALL dbcsr_set(propagator(propagator_sector_virtual, jquad - 1, 1)%matrix, 0.0_dp)
            CALL dbcsr_filter(propagator(propagator_sector_occupied, jquad - 1, 1)%matrix, 0.0_dp)
            CALL dbcsr_filter(propagator(propagator_sector_virtual, jquad - 1, 1)%matrix, 0.0_dp)
         END IF

      END IF ! do kpoints

      CALL timestop(handle)

   END SUBROUTINE compute_mat_dm_global

! **************************************************************************************************
!> \brief Builds the Gamma-point occupied and virtual propagators.
!> \param propagator sector/time/cell propagator matrix set
!> \param jquad time-point index
!> \param cfm_mo_coeff complete MO coefficient matrix
!> \param homo number of occupied orbitals
!> \param Eigenval orbital energies
!> \param nmo number of orbitals
!> \param eps_filter DBCSR filtering threshold
!> \param e_fermi Fermi energy
!> \param tau imaginary time
!> \param para_env parallel environment
! **************************************************************************************************
   SUBROUTINE compute_gamma_propagator(propagator, jquad, cfm_mo_coeff, homo, Eigenval, nmo, &
                                       eps_filter, e_fermi, tau, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), &
         INTENT(INOUT), POINTER                          :: propagator
      INTEGER, INTENT(IN)                                :: jquad
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_mo_coeff
      INTEGER, INTENT(IN)                                :: homo
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: Eigenval
      INTEGER, INTENT(IN)                                :: nmo
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter, e_fermi, tau
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_gamma_propagator'
      REAL(KIND=dp), PARAMETER                           :: stabilize_exp = 70.0_dp

      INTEGER                                            :: handle, i_global, iiB, jjB, nao, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      TYPE(cp_cfm_type)                                  :: cfm_mo_coeff_occ_scaled, &
                                                            cfm_mo_coeff_virt_scaled, &
                                                            cfm_scaled_dm_occ_tau, &
                                                            cfm_scaled_dm_virt_tau
      TYPE(cp_fm_type)                                   :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(dbcsr_type), POINTER                          :: mat_occ, mat_virt

      CALL timeset(routineN, handle)
      CALL para_env%sync()

      mat_occ => propagator(propagator_sector_occupied, jquad, 1)%matrix
      mat_virt => propagator(propagator_sector_virtual, jquad, 1)%matrix

      CALL cp_cfm_get_info(matrix=cfm_mo_coeff, &
                           nrow_global=nao, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      CALL cp_cfm_create(cfm_scaled_dm_occ_tau, cfm_mo_coeff%matrix_struct, nrow=nao, ncol=nao)
      CALL cp_cfm_create(cfm_scaled_dm_virt_tau, cfm_mo_coeff%matrix_struct, nrow=nao, ncol=nao)
      CALL cp_cfm_create(cfm_mo_coeff_occ_scaled, cfm_mo_coeff%matrix_struct)
      CALL cp_cfm_create(cfm_mo_coeff_virt_scaled, cfm_mo_coeff%matrix_struct)

      DO jjB = 1, nrow_local
         DO iiB = 1, ncol_local
            i_global = col_indices(iiB)
            IF (i_global <= homo .AND. ABS(tau*0.5_dp*(Eigenval(i_global) - e_fermi)) < stabilize_exp) THEN
               cfm_mo_coeff_occ_scaled%local_data(jjB, iiB) = &
                  cfm_mo_coeff%local_data(jjB, iiB)*EXP(tau*0.5_dp*(Eigenval(i_global) - e_fermi))
            ELSE
               cfm_mo_coeff_occ_scaled%local_data(jjB, iiB) = (0.0_dp, 0.0_dp)
            END IF
         END DO
      END DO

      DO jjB = 1, nrow_local
         DO iiB = 1, ncol_local
            i_global = col_indices(iiB)
            IF (i_global > homo .AND. ABS(tau*0.5_dp*(Eigenval(i_global) - e_fermi)) < stabilize_exp) THEN
               cfm_mo_coeff_virt_scaled%local_data(jjB, iiB) = &
                  cfm_mo_coeff%local_data(jjB, iiB)*EXP(-tau*0.5_dp*(Eigenval(i_global) - e_fermi))
            ELSE
               cfm_mo_coeff_virt_scaled%local_data(jjB, iiB) = (0.0_dp, 0.0_dp)
            END IF
         END DO
      END DO

      CALL para_env%sync()

      CALL parallel_gemm(transa="N", transb="C", m=nao, n=nao, k=nmo, alpha=(1.0_dp, 0.0_dp), &
                         matrix_a=cfm_mo_coeff_occ_scaled, matrix_b=cfm_mo_coeff_occ_scaled, beta=(0.0_dp, 0.0_dp), &
                         matrix_c=cfm_scaled_dm_occ_tau)
      CALL parallel_gemm(transa="N", transb="C", m=nao, n=nao, k=nmo, alpha=(1.0_dp, 0.0_dp), &
                         matrix_a=cfm_mo_coeff_virt_scaled, matrix_b=cfm_mo_coeff_virt_scaled, beta=(0.0_dp, 0.0_dp), &
                         matrix_c=cfm_scaled_dm_virt_tau)

      CALL cp_fm_create(fm_scaled_dm_occ_tau, cfm_scaled_dm_occ_tau%matrix_struct)
      CALL cp_fm_create(fm_scaled_dm_virt_tau, cfm_scaled_dm_virt_tau%matrix_struct)
      CALL cp_cfm_to_fm(cfm_scaled_dm_occ_tau, fm_scaled_dm_occ_tau)
      CALL cp_cfm_to_fm(cfm_scaled_dm_virt_tau, fm_scaled_dm_virt_tau)

      CALL dbcsr_set(mat_occ, 0.0_dp)
      CALL copy_fm_to_dbcsr(fm_scaled_dm_occ_tau, mat_occ, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(mat_occ, eps_filter)
      CALL dbcsr_set(mat_virt, 0.0_dp)
      CALL copy_fm_to_dbcsr(fm_scaled_dm_virt_tau, mat_virt, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(mat_virt, eps_filter)

      CALL cp_fm_release(fm_scaled_dm_occ_tau)
      CALL cp_fm_release(fm_scaled_dm_virt_tau)
      CALL cp_cfm_release(cfm_scaled_dm_occ_tau)
      CALL cp_cfm_release(cfm_scaled_dm_virt_tau)
      CALL cp_cfm_release(cfm_mo_coeff_occ_scaled)
      CALL cp_cfm_release(cfm_mo_coeff_virt_scaled)
      CALL timestop(handle)

   END SUBROUTINE compute_gamma_propagator

! **************************************************************************************************
!> \brief Calculate one kpoint density matrix in the complex AO representation.
!> \param kpoint    kpoint environment
!> \param ikpgr     local kpoint index
!> \param ispin     spin index
!> \param tau ...
!> \param e_fermi ...
!> \param sector ...
!> \param density complex AO density matrix
! **************************************************************************************************
   SUBROUTINE build_kpoint_density_matrix_rpa(kpoint, ikpgr, ispin, tau, e_fermi, sector, density)

      TYPE(kpoint_type), POINTER                         :: kpoint
      INTEGER, INTENT(IN)                                :: ikpgr, ispin
      REAL(KIND=dp), INTENT(IN)                          :: tau, e_fermi
      INTEGER, INTENT(IN)                                :: sector
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: density

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_kpoint_density_matrix_rpa'
      REAL(KIND=dp), PARAMETER                           :: stabilize_exp = 70.0_dp

      INTEGER                                            :: handle, i_mo, nao, nmo, nspin
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, exp_scaling, occupation
      TYPE(cp_cfm_type)                                  :: cfm_mo_coeff, cfm_work
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mo_set, mo_set_im

      CALL timeset(routineN, handle)

      CPASSERT(sector == propagator_sector_occupied .OR. sector == propagator_sector_virtual)
      CPASSERT(kpoint%use_real_wfn .EQV. .FALSE.)

      kp => kpoint%kp_env(ikpgr)%kpoint_env
      nspin = SIZE(kp%mos, 2)
      CPASSERT(ispin >= 1 .AND. ispin <= nspin)
      mo_set => kp%mos(1, ispin)
      mo_set_im => kp%mos(2, ispin)
      CALL get_mo_set(mo_set, nao=nao, nmo=nmo, eigenvalues=eigenvalues, occupation_numbers=occupation)

      ! if this CPASSERT is triggered, please add all virtual MOs to SCF section,
      ! e.g. ADDED_MOS 1000000
      CPASSERT(nao == nmo)

      ALLOCATE (exp_scaling(nmo))
      CALL cp_cfm_create(cfm_mo_coeff, mo_set%mo_coeff%matrix_struct)
      CALL cp_cfm_create(cfm_work, mo_set%mo_coeff%matrix_struct)
      CALL cp_fm_to_cfm(msourcer=mo_set%mo_coeff, msourcei=mo_set_im%mo_coeff, &
                        mtarget=cfm_mo_coeff)
      CALL cp_cfm_to_cfm(cfm_mo_coeff, cfm_work)

      IF (sector == propagator_sector_occupied) THEN
         CALL cp_cfm_column_scale(cfm_work, CMPLX(occupation, 0.0_dp, KIND=dp))
      ELSE
         CALL cp_cfm_column_scale(cfm_work, &
                                  CMPLX(2.0_dp/REAL(nspin, KIND=dp) - occupation, 0.0_dp, KIND=dp))
      END IF

      IF (nspin == 1) THEN
         CALL cp_cfm_scale(0.5_dp, cfm_work)
      END IF

      DO i_mo = 1, nmo
         IF (ABS(tau*0.5_dp*(eigenvalues(i_mo) - e_fermi)) < stabilize_exp) THEN
            exp_scaling(i_mo) = EXP(-ABS(tau*(eigenvalues(i_mo) - e_fermi)))
         ELSE
            exp_scaling(i_mo) = 0.0_dp
         END IF
      END DO

      CALL cp_cfm_column_scale(cfm_work, CMPLX(exp_scaling, 0.0_dp, KIND=dp))
      CALL parallel_gemm("N", "C", nao, nao, nmo, (1.0_dp, 0.0_dp), cfm_mo_coeff, cfm_work, &
                         (0.0_dp, 0.0_dp), density)

      CALL cp_cfm_release(cfm_work)
      CALL cp_cfm_release(cfm_mo_coeff)
      DEALLOCATE (exp_scaling)

      CALL timestop(handle)

   END SUBROUTINE build_kpoint_density_matrix_rpa

! **************************************************************************************************
!> \brief ...
!> \param propagator ...
!> \param qs_env ...
!> \param ispin ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param e_fermi ...
!> \param tau ...
!> \param eps_filter ...
!> \param num_cells_dm ...
!> \param index_to_cell_dm ...
!> \param sector ...
! **************************************************************************************************
   SUBROUTINE compute_transl_dm(propagator, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                eps_filter, num_cells_dm, index_to_cell_dm, &
                                sector)
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), &
         INTENT(INOUT), POINTER                          :: propagator
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ispin, num_integ_points, jquad
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, tau, eps_filter
      INTEGER, INTENT(OUT)                               :: num_cells_dm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      INTEGER, INTENT(IN)                                :: sector

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_transl_dm'

      INTEGER                                            :: handle, i_dim, i_img, nao
      INTEGER, DIMENSION(3)                              :: cell_grid_dm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type)                                  :: cfm_density
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_dm_global_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      matrix_s_kp=matrix_s_kp, &
                      kpoints=kpoints, &
                      cell=cell)

      ! we always use an odd number of image cells
      ! CAUTION: also at another point, cell_grid_dm is defined, these definitions have to be identical
      DO i_dim = 1, 3
         cell_grid_dm(i_dim) = (kpoints%nkp_grid(i_dim)/2)*2 - 1
      END DO

      num_cells_dm = cell_grid_dm(1)*cell_grid_dm(2)*cell_grid_dm(3)

      NULLIFY (mat_dm_global_work)
      ! Only the requested spin is built.
      CALL dbcsr_allocate_matrix_set(mat_dm_global_work, num_cells_dm)

      DO i_img = 1, num_cells_dm

         ALLOCATE (mat_dm_global_work(i_img)%matrix)
         CALL dbcsr_create(matrix=mat_dm_global_work(i_img)%matrix, &
                           template=matrix_s_kp(1, 1)%matrix, &
                           !                           matrix_type=dbcsr_type_symmetric)
                           matrix_type=dbcsr_type_no_symmetry)

         CALL dbcsr_reserve_all_blocks(mat_dm_global_work(i_img)%matrix)

      END DO

      ! Reusable complex AO scratch for the k-point density source.
      CALL get_mo_set(kpoints%kp_env(1)%kpoint_env%mos(1, 1), nao=nao)
      CALL cp_cfm_create(cfm_density, kpoints%kp_env(1)%kpoint_env%mos(1, 1)%mo_coeff%matrix_struct, &
                         nrow=nao, ncol=nao)

      ! overwrite the cell indices in kpoints
      CALL init_cell_index_rpa(cell_grid_dm, kpoints%cell_to_index, kpoints%index_to_cell, cell)

      ! density matrices in real space, the cell vectors T for transforming are taken from kpoints%index_to_cell
      ! (custom made for RPA) and not from sab_nl (which is symmetric and from SCF)
      CALL density_matrix_from_kp_to_transl(kpoints, cfm_density, mat_dm_global_work, kpoints%index_to_cell, &
                                            ispin, tau, e_fermi, sector)

      CALL cp_cfm_release(cfm_density)

      ! we need the index to cell for the density matrices later
      index_to_cell_dm => kpoints%index_to_cell

      ! The first request owns allocation of the complete time-indexed propagator set.
      IF (.NOT. ASSOCIATED(propagator)) THEN
         CALL create_propagator_matrix_set(propagator, num_integ_points, matrix_s_kp(1, 1)%matrix, &
                                           kpoints%index_to_cell)
      END IF

      DO i_img = 1, num_cells_dm

         ! filter to get rid of the blocks full with zeros on the lower half, otherwise blocks doubled
         CALL dbcsr_filter(mat_dm_global_work(i_img)%matrix, eps_filter)

         CALL dbcsr_copy(propagator(sector, jquad, i_img)%matrix, &
                         mat_dm_global_work(i_img)%matrix)

      END DO

      CALL dbcsr_deallocate_matrix_set(mat_dm_global_work)

      CALL timestop(handle)

   END SUBROUTINE compute_transl_dm

! **************************************************************************************************
!> \brief ...
!> \param propagator ...
!> \param qs_env ...
!> \param ispin ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param e_fermi ...
!> \param tau ...
!> \param sector ...
! **************************************************************************************************
   SUBROUTINE compute_periodic_dm(propagator, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                  sector)
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), &
         INTENT(INOUT), POINTER                          :: propagator
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ispin, num_integ_points, jquad
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi, tau
      INTEGER, INTENT(IN)                                :: sector

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_periodic_dm'

      INTEGER                                            :: handle, nao, num_cells_dm
      INTEGER, DIMENSION(3, 1)                           :: index_to_cell_zero
      TYPE(cp_cfm_type)                                  :: cfm_density
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_dm_global_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      TYPE(kpoint_type), POINTER                         :: kpoints_G

      CALL timeset(routineN, handle)

      NULLIFY (matrix_s_kp)

      CALL get_qs_env(qs_env, &
                      matrix_s_kp=matrix_s_kp)

      kpoints_G => qs_env%mp2_env%ri_rpa_im_time%kpoints_G

      num_cells_dm = 1
      index_to_cell_zero = 0

      NULLIFY (mat_dm_global_work)
      ! Only the requested spin is built.
      CALL dbcsr_allocate_matrix_set(mat_dm_global_work, num_cells_dm)

      ! The first request owns allocation of the complete time-indexed propagator set.
      IF (.NOT. ASSOCIATED(propagator)) THEN
         CALL create_propagator_matrix_set(propagator, num_integ_points, matrix_s_kp(1, 1)%matrix, &
                                           index_to_cell_zero)
      END IF

      ALLOCATE (mat_dm_global_work(1)%matrix)
      CALL dbcsr_create(matrix=mat_dm_global_work(1)%matrix, &
                        template=matrix_s_kp(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_reserve_all_blocks(mat_dm_global_work(1)%matrix)

      CALL dbcsr_set(mat_dm_global_work(1)%matrix, 0.0_dp)

      ! Reusable complex AO scratch for the k-point density source.
      CALL get_mo_set(kpoints_G%kp_env(1)%kpoint_env%mos(1, 1), nao=nao)
      CALL cp_cfm_create(cfm_density, kpoints_G%kp_env(1)%kpoint_env%mos(1, 1)%mo_coeff%matrix_struct, &
                         nrow=nao, ncol=nao)

      CALL density_matrix_from_kp_to_mic(kpoints_G, cfm_density, mat_dm_global_work, qs_env, ispin, tau, e_fermi, sector)

      CALL cp_cfm_release(cfm_density)

      CALL dbcsr_copy(propagator(sector, jquad, 1)%matrix, &
                      mat_dm_global_work(1)%matrix)

      CALL dbcsr_deallocate_matrix_set(mat_dm_global_work)

      CALL timestop(handle)

   END SUBROUTINE compute_periodic_dm

   ! **************************************************************************************************
!> \brief ...
!> \param kpoints_G ...
!> \param density complex AO density scratch matrix
!> \param mat_dm_global_work ...
!> \param qs_env ...
!> \param ispin ...
!> \param tau ...
!> \param e_fermi ...
!> \param sector ...
! **************************************************************************************************
   SUBROUTINE density_matrix_from_kp_to_mic(kpoints_G, density, mat_dm_global_work, qs_env, ispin, tau, e_fermi, sector)

      TYPE(kpoint_type), POINTER                         :: kpoints_G
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: density
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_dm_global_work
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), INTENT(IN)                          :: tau, e_fermi
      INTEGER, INTENT(IN)                                :: sector

      CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_from_kp_to_mic'

      COMPLEX(KIND=dp), CONTIGUOUS, DIMENSION(:, :), POINTER :: cfm_data
      INTEGER                                            :: handle, iatom, iatom_old, ik, ik_global, ik_old, &
                                                            irow, jatom, jatom_old, jcol, &
                                                            nao, ncol_local, nrow_local, &
                                                            num_cells
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_ao_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp)                                      :: contribution, weight_im, weight_re
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_mat_work
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (xkp, wkp)

      CALL get_kpoint_info(kpoints_G, xkp=xkp, wkp=wkp)
      index_to_cell => kpoints_G%index_to_cell
      num_cells = SIZE(index_to_cell, 2)

      CALL cp_cfm_get_info(density, nrow_global=nao, nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)
      CALL cp_fm_create(fm_mat_work, density%matrix_struct, set_zero=.TRUE.)

      ALLOCATE (atom_from_ao_index(nao))

      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_ao_index, nao, "ORB")

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      iatom_old = 0
      jatom_old = 0
      ik_old = 0

      CALL dbcsr_set(mat_dm_global_work(1)%matrix, 0.0_dp)

      DO ik = 1, SIZE(kpoints_G%kp_env)

         kp => kpoints_G%kp_env(ik)%kpoint_env
         ik_global = kp%nkpoint
         CALL build_kpoint_density_matrix_rpa(kpoints_G, ik, ispin, tau, e_fermi, sector, density)
         CALL cp_cfm_get_info(density, local_data=cfm_data)

         DO irow = 1, nrow_local
            DO jcol = 1, ncol_local

               iatom = atom_from_ao_index(row_indices(irow))
               jatom = atom_from_ao_index(col_indices(jcol))

               IF (ik_global /= ik_old .OR. iatom /= iatom_old .OR. jatom /= jatom_old) THEN

                  CALL compute_weight_re_im(weight_re, weight_im, &
                                            num_cells, iatom, jatom, xkp(1:3, ik_global), wkp(ik_global), &
                                            cell, index_to_cell, hmat, particle_set)

                  iatom_old = iatom
                  jatom_old = jatom
                  ik_old = ik_global

               END IF

               contribution = REAL(CMPLX(weight_re, -weight_im, KIND=dp)*cfm_data(irow, jcol), KIND=dp)

               fm_mat_work%local_data(irow, jcol) = fm_mat_work%local_data(irow, jcol) + contribution

            END DO
         END DO

      END DO ! ik

      CALL copy_fm_to_dbcsr(fm_mat_work, mat_dm_global_work(1)%matrix, keep_sparsity=.FALSE.)

      CALL cp_fm_set_all(fm_mat_work, 0.0_dp)

      CALL cp_fm_release(fm_mat_work)
      DEALLOCATE (atom_from_ao_index)

      CALL timestop(handle)

   END SUBROUTINE density_matrix_from_kp_to_mic

! **************************************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param density complex AO density scratch matrix
!> \param mat_dm_global_work ...
!> \param index_to_cell ...
!> \param ispin ...
!> \param tau ...
!> \param e_fermi ...
!> \param sector ...
! **************************************************************************************************
   SUBROUTINE density_matrix_from_kp_to_transl(kpoints, density, mat_dm_global_work, index_to_cell, &
                                               ispin, tau, e_fermi, sector)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: density
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN)       :: mat_dm_global_work
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp), INTENT(IN)                          :: tau, e_fermi
      INTEGER, INTENT(IN)                                :: sector

      CHARACTER(LEN=*), PARAMETER :: routineN = 'density_matrix_from_kp_to_transl'

      INTEGER                                            :: handle, icell, ik, ik_global, xcell, ycell, zcell
      REAL(KIND=dp)                                      :: arg, coskl, sinkl
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cp_fm_type)                                   :: fm_mat_im, fm_mat_re
      TYPE(dbcsr_type), POINTER                          :: mat_work_im, mat_work_re
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL timeset(routineN, handle)

      NULLIFY (xkp, wkp)

      NULLIFY (mat_work_re)
      CALL dbcsr_init_p(mat_work_re)
      CALL dbcsr_create(matrix=mat_work_re, &
                        template=mat_dm_global_work(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      NULLIFY (mat_work_im)
      CALL dbcsr_init_p(mat_work_im)
      CALL dbcsr_create(matrix=mat_work_im, &
                        template=mat_dm_global_work(1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      CALL cp_fm_create(fm_mat_re, density%matrix_struct)
      CALL cp_fm_create(fm_mat_im, density%matrix_struct)

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp)

      CPASSERT(SIZE(mat_dm_global_work) == SIZE(index_to_cell, 2))

      DO icell = 1, SIZE(mat_dm_global_work)

         CALL dbcsr_set(mat_dm_global_work(icell)%matrix, 0.0_dp)

      END DO

      DO ik = 1, SIZE(kpoints%kp_env)

         kp => kpoints%kp_env(ik)%kpoint_env
         ik_global = kp%nkpoint
         CALL build_kpoint_density_matrix_rpa(kpoints, ik, ispin, tau, e_fermi, sector, density)
         CALL cp_cfm_to_fm(density, fm_mat_re, fm_mat_im)

         CALL copy_fm_to_dbcsr(fm_mat_re, mat_work_re, keep_sparsity=.FALSE.)
         CALL copy_fm_to_dbcsr(fm_mat_im, mat_work_im, keep_sparsity=.FALSE.)

         DO icell = 1, SIZE(mat_dm_global_work)

            xcell = index_to_cell(1, icell)
            ycell = index_to_cell(2, icell)
            zcell = index_to_cell(3, icell)

            arg = REAL(xcell, dp)*xkp(1, ik_global) + REAL(ycell, dp)*xkp(2, ik_global) + &
                  REAL(zcell, dp)*xkp(3, ik_global)
            coskl = wkp(ik_global)*COS(twopi*arg)
            sinkl = wkp(ik_global)*SIN(twopi*arg)

            CALL dbcsr_add(mat_dm_global_work(icell)%matrix, mat_work_re, 1.0_dp, coskl)
            ! The real-space convention is Re(exp(i k R) G(k)).
            CALL dbcsr_add(mat_dm_global_work(icell)%matrix, mat_work_im, 1.0_dp, -sinkl)

         END DO

      END DO

      CALL dbcsr_release_p(mat_work_re)
      CALL dbcsr_release_p(mat_work_im)
      CALL cp_fm_release(fm_mat_re)
      CALL cp_fm_release(fm_mat_im)

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
      TYPE(cell_type), INTENT(IN), POINTER               :: cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_cell_index_rpa'

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
      CPASSERT(cell_grid(1) /= itm(1)*2)
      CPASSERT(cell_grid(2) /= itm(2)*2)
      CPASSERT(cell_grid(3) /= itm(3)*2)

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
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: index_to_cell_3c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: cell_to_index_3c
      INTEGER, DIMENSION(:, :), INTENT(IN), POINTER      :: index_to_cell_dm
      LOGICAL, INTENT(OUT)                               :: R_minus_S_needed
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_diff_index_3c'

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

         IF (x_cell_R_minus_S >= LBOUND(cell_to_index_3c, 1) .AND. &
             x_cell_R_minus_S <= UBOUND(cell_to_index_3c, 1) .AND. &
             y_cell_R_minus_S >= LBOUND(cell_to_index_3c, 2) .AND. &
             y_cell_R_minus_S <= UBOUND(cell_to_index_3c, 2) .AND. &
             z_cell_R_minus_S >= LBOUND(cell_to_index_3c, 3) .AND. &
             z_cell_R_minus_S <= UBOUND(cell_to_index_3c, 3)) THEN

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
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell_dm
      LOGICAL, INTENT(OUT)                               :: R_minus_S_minus_T_needed
      LOGICAL, INTENT(IN)                                :: do_kpoints_cubic_RPA

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_diff_diff_index_3c'

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

         IF (x_cell_R_minus_S_minus_T >= LBOUND(cell_to_index_3c, 1) .AND. &
             x_cell_R_minus_S_minus_T <= UBOUND(cell_to_index_3c, 1) .AND. &
             y_cell_R_minus_S_minus_T >= LBOUND(cell_to_index_3c, 2) .AND. &
             y_cell_R_minus_S_minus_T <= UBOUND(cell_to_index_3c, 2) .AND. &
             z_cell_R_minus_S_minus_T >= LBOUND(cell_to_index_3c, 3) .AND. &
             z_cell_R_minus_S_minus_T <= UBOUND(cell_to_index_3c, 3)) THEN

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
