!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines from paper [Graml2024]
!> \author Jan Wilhelm
!> \date 07.2023
! **************************************************************************************************
MODULE gw_large_cell_gamma
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_block_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_diag, cp_fm_get_info, cp_fm_read_unformatted, cp_fm_release, &
        cp_fm_set_all, cp_fm_to_fm, cp_fm_type, cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE dbt_api,                         ONLY: dbt_clear,&
                                              dbt_contract,&
                                              dbt_copy,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE gw_communication,                ONLY: fm_to_local_tensor,&
                                              local_dbt_to_global_mat
   USE gw_utils,                        ONLY: analyt_conti_and_print,&
                                              de_init_bs_env,&
                                              time_to_freq
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_file_delete
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: MIC_contribution_from_ikp,&
                                              cfm_ikp_from_fm_Gamma,&
                                              get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_tensors,                      ONLY: build_3c_integrals
   USE rpa_gw_kpoints_util,             ONLY: cp_cfm_power,&
                                              cp_cfm_upper_to_full
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_large_cell_gamma'

   PUBLIC :: gw_calc_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief Perform GW band structure calculation
!> \param qs_env ...
!> \param bs_env ...
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE gw_calc_large_cell_Gamma(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_large_cell_Gamma'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_MIC_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      ! G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      ! χ_PQ(iτ,k=0) = sum_λν [sum_µ (µν|P) G^occ_µλ(i|τ|)] [sum_σ (σλ|Q) G^vir_σν(i|τ|)]
      CALL get_mat_chi_Gamma_tau(bs_env, qs_env, bs_env%mat_chi_Gamma_tau)

      ! χ_PQ(iτ,k=0) -> χ_PQ(iω,k) -> ε_PQ(iω,k) -> W_PQ(iω,k) -> W^MIC_PQ(iτ) -> M^-1*W^MIC*M^-1
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_MIC_time)

      ! D_µν = sum_n^occ C_µn(k=0) C_νn(k=0), V^trunc_PQ = sum_cell_R <phi_P,0|V^trunc|phi_Q,R>
      ! Σ^x_λσ(k=0) = sum_νQ [sum_P (νσ|P) V^trunc_PQ] [sum_µ (λµ|Q) D_µν)]
      CALL get_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma)

      ! Σ^c_λσ(iτ,k=0) = sum_νQ [sum_P (νσ|P) W^MIC_PQ(iτ)] [sum_µ (λµ|Q) G^occ_µν(i|τ|)], τ < 0
      ! Σ^c_λσ(iτ,k=0) = sum_νQ [sum_P (νσ|P) W^MIC_PQ(iτ)] [sum_µ (λµ|Q) G^vir_µν(i|τ|)], τ > 0
      CALL get_Sigma_c(bs_env, qs_env, fm_W_MIC_time, fm_Sigma_c_Gamma_time)

      ! Σ^c_λσ(iτ,k=0) -> Σ^c_nn(ϵ,k); ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ,k) + Σ^x_nn(k) - v^xc_nn(k)
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_large_cell_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
! **************************************************************************************************
   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, qs_env, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER :: handle, i_intval_idx, i_t, inner_loop_atoms_interval_index, ispin, j_intval_idx
      INTEGER, DIMENSION(2)                              :: i_atoms, IL_atoms, j_atoms
      LOGICAL                                            :: dist_too_long_i, dist_too_long_j
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, &
                                                            t_3c_for_Gvir, t_3c_x_Gocc, &
                                                            t_3c_x_Gocc_2, t_3c_x_Gvir, &
                                                            t_3c_x_Gvir_2

      CALL timeset(routineN, handle)

      DO i_t = 1, bs_env%num_time_freq_points

         t1 = m_walltime()

         IF (bs_env%read_chi(i_t)) THEN

            CALL fm_read(bs_env%fm_RI_RI, bs_env, bs_env%chi_name, i_t)

            CALL copy_fm_to_dbcsr(bs_env%fm_RI_RI, mat_chi_Gamma_tau(i_t)%matrix, &
                                  keep_sparsity=.FALSE.)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I5,A,I3,A,F7.1,A)') &
                  'Read χ(iτ,k=0) from file for time point  ', i_t, ' /', &
                  bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

            CYCLE

         END IF

         IF (.NOT. bs_env%calc_chi(i_t)) CYCLE

         CALL create_tensors_chi(t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, t_3c_for_Gvir, &
                                 t_3c_x_Gocc, t_3c_x_Gvir, t_3c_x_Gocc_2, t_3c_x_Gvir_2, bs_env)

         ! 1. compute G^occ and G^vir
         !    Background: G^σ(iτ) = G^occ,σ(iτ) * Θ(-τ) + G^vir,σ(iτ) * Θ(τ), σ ∈ {↑,↓}
         !    G^occ,σ_µλ(i|τ|,k=0) = sum_n^occ C^σ_µn(k=0) e^(-|(ϵ^σ_nk=0-ϵ_F)τ|) C^σ_λn(k=0)
         !    G^vir,σ_µλ(i|τ|,k=0) = sum_n^vir C^σ_µn(k=0) e^(-|(ϵ^σ_nk=0-ϵ_F)τ|) C^σ_λn(k=0)
         tau = bs_env%imag_time_points(i_t)

         DO ispin = 1, bs_env%n_spin
            CALL G_occ_vir(bs_env, tau, bs_env%fm_Gocc, ispin, occ=.TRUE., vir=.FALSE.)
            CALL G_occ_vir(bs_env, tau, bs_env%fm_Gvir, ispin, occ=.FALSE., vir=.TRUE.)

            CALL fm_to_local_tensor(bs_env%fm_Gocc, bs_env%mat_ao_ao%matrix, &
                                    bs_env%mat_ao_ao_tensor%matrix, t_2c_Gocc, bs_env, &
                                    bs_env%atoms_j_t_group)
            CALL fm_to_local_tensor(bs_env%fm_Gvir, bs_env%mat_ao_ao%matrix, &
                                    bs_env%mat_ao_ao_tensor%matrix, t_2c_Gvir, bs_env, &
                                    bs_env%atoms_i_t_group)

            ! every group has its own range of i_atoms and j_atoms; only deal with a
            ! limited number of i_atom-j_atom pairs simultaneously in a group to save memory
            DO i_intval_idx = 1, bs_env%n_intervals_i
               DO j_intval_idx = 1, bs_env%n_intervals_j
                  i_atoms = bs_env%i_atom_intervals(1:2, i_intval_idx)
                  j_atoms = bs_env%j_atom_intervals(1:2, j_intval_idx)

                  DO inner_loop_atoms_interval_index = 1, bs_env%n_intervals_inner_loop_atoms

                     IL_atoms = bs_env%inner_loop_atom_intervals(1:2, inner_loop_atoms_interval_index)

                     CALL check_dist(i_atoms, IL_atoms, qs_env, bs_env, dist_too_long_i)
                     CALL check_dist(j_atoms, IL_atoms, qs_env, bs_env, dist_too_long_j)
                     IF (dist_too_long_i .OR. dist_too_long_j) CYCLE

                     ! 2. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
                     CALL compute_3c_integrals(qs_env, bs_env, t_3c_for_Gocc, i_atoms, IL_atoms)

                     ! 3. tensor operation M_λνP(iτ) = sum_µ (µν|P) G^occ_µλ(i|τ|,k=0)
                     CALL G_times_3c(t_3c_for_Gocc, t_2c_Gocc, t_3c_x_Gocc, bs_env, &
                                     j_atoms, i_atoms, IL_atoms)

                     ! 4. compute 3-center integrals (σλ|Q) ("|": truncated Coulomb operator)
                     CALL compute_3c_integrals(qs_env, bs_env, t_3c_for_Gvir, j_atoms, IL_atoms)

                     ! 5. tensor operation N_νλQ(iτ) = sum_σ (σλ|Q) G^vir_σν(i|τ|,k=0)
                     CALL G_times_3c(t_3c_for_Gvir, t_2c_Gvir, t_3c_x_Gvir, bs_env, &
                                     i_atoms, j_atoms, IL_atoms)

                  END DO ! IL_atoms

                  ! 6. reorder tensors
                  CALL dbt_copy(t_3c_x_Gocc, t_3c_x_Gocc_2, move_data=.TRUE., order=[1, 3, 2])
                  CALL dbt_copy(t_3c_x_Gvir, t_3c_x_Gvir_2, move_data=.TRUE.)

                  ! 7. tensor operation χ_PQ(iτ,k=0) = sum_λν M_λνP(iτ) N_νλQ(iτ),
                  CALL dbt_contract(alpha=bs_env%spin_degeneracy, &
                                    tensor_1=t_3c_x_Gocc_2, tensor_2=t_3c_x_Gvir_2, &
                                    beta=1.0_dp, tensor_3=bs_env%t_chi, &
                                    contract_1=[2, 3], notcontract_1=[1], map_1=[1], &
                                    contract_2=[2, 3], notcontract_2=[1], map_2=[2], &
                                    filter_eps=bs_env%eps_filter, move_data=.TRUE.)

               END DO ! j_atoms
            END DO ! i_atoms
         END DO ! ispin

         ! 8. communicate data of χ_PQ(iτ,k=0) in tensor bs_env%t_chi (which local in the
         !    subgroup) to the global dbcsr matrix mat_chi_Gamma_tau (which stores
         !    χ_PQ(iτ,k=0) for all time points)
         CALL local_dbt_to_global_mat(bs_env%t_chi, bs_env%mat_RI_RI_tensor%matrix, &
                                      mat_chi_Gamma_tau(i_t)%matrix, bs_env%para_env)

         CALL write_matrix(mat_chi_Gamma_tau(i_t)%matrix, i_t, bs_env%chi_name, &
                           bs_env%fm_RI_RI, qs_env)

         CALL destroy_tensors_chi(t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, t_3c_for_Gvir, &
                                  t_3c_x_Gocc, t_3c_x_Gvir, t_3c_x_Gocc_2, t_3c_x_Gvir_2)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ(iτ,k=0) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE get_mat_chi_Gamma_tau

! **************************************************************************************************
!> \brief ...
!> \param fm ...
!> \param bs_env ...
!> \param mat_name ...
!> \param idx ...
! **************************************************************************************************
   SUBROUTINE fm_read(fm, bs_env, mat_name, idx)
      TYPE(cp_fm_type)                                   :: fm
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*)                                   :: mat_name
      INTEGER                                            :: idx

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_read'

      CHARACTER(LEN=default_string_length)               :: f_chi
      INTEGER                                            :: handle, unit_nr

      CALL timeset(routineN, handle)

      unit_nr = -1
      IF (bs_env%para_env%is_source()) THEN

         IF (idx < 10) THEN
            WRITE (f_chi, '(3A,I1,A)') TRIM(bs_env%prefix), TRIM(mat_name), "_0", idx, ".matrix"
         ELSE IF (idx < 100) THEN
            WRITE (f_chi, '(3A,I2,A)') TRIM(bs_env%prefix), TRIM(mat_name), "_", idx, ".matrix"
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         CALL open_file(file_name=TRIM(f_chi), file_action="READ", file_form="UNFORMATTED", &
                        file_position="REWIND", file_status="OLD", unit_number=unit_nr)

      END IF

      CALL cp_fm_read_unformatted(fm, unit_nr)

      IF (bs_env%para_env%is_source()) CALL close_file(unit_number=unit_nr)

      CALL timestop(handle)

   END SUBROUTINE fm_read

! **************************************************************************************************
!> \brief ...
!> \param t_2c_Gocc ...
!> \param t_2c_Gvir ...
!> \param t_3c_for_Gocc ...
!> \param t_3c_for_Gvir ...
!> \param t_3c_x_Gocc ...
!> \param t_3c_x_Gvir ...
!> \param t_3c_x_Gocc_2 ...
!> \param t_3c_x_Gvir_2 ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE create_tensors_chi(t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, t_3c_for_Gvir, &
                                 t_3c_x_Gocc, t_3c_x_Gvir, t_3c_x_Gocc_2, t_3c_x_Gvir_2, bs_env)

      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, &
                                                            t_3c_for_Gvir, t_3c_x_Gocc, &
                                                            t_3c_x_Gvir, t_3c_x_Gocc_2, &
                                                            t_3c_x_Gvir_2
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_tensors_chi'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_create(bs_env%t_G, t_2c_Gocc)
      CALL dbt_create(bs_env%t_G, t_2c_Gvir)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_for_Gocc)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_for_Gvir)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_Gocc)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_Gvir)
      CALL dbt_create(bs_env%t_RI__AO_AO, t_3c_x_Gocc_2)
      CALL dbt_create(bs_env%t_RI__AO_AO, t_3c_x_Gvir_2)

      CALL timestop(handle)

   END SUBROUTINE create_tensors_chi

! **************************************************************************************************
!> \brief ...
!> \param t_2c_Gocc ...
!> \param t_2c_Gvir ...
!> \param t_3c_for_Gocc ...
!> \param t_3c_for_Gvir ...
!> \param t_3c_x_Gocc ...
!> \param t_3c_x_Gvir ...
!> \param t_3c_x_Gocc_2 ...
!> \param t_3c_x_Gvir_2 ...
! **************************************************************************************************
   SUBROUTINE destroy_tensors_chi(t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, t_3c_for_Gvir, &
                                  t_3c_x_Gocc, t_3c_x_Gvir, t_3c_x_Gocc_2, t_3c_x_Gvir_2)
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, &
                                                            t_3c_for_Gvir, t_3c_x_Gocc, &
                                                            t_3c_x_Gvir, t_3c_x_Gocc_2, &
                                                            t_3c_x_Gvir_2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_tensors_chi'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_destroy(t_2c_Gocc)
      CALL dbt_destroy(t_2c_Gvir)
      CALL dbt_destroy(t_3c_for_Gocc)
      CALL dbt_destroy(t_3c_for_Gvir)
      CALL dbt_destroy(t_3c_x_Gocc)
      CALL dbt_destroy(t_3c_x_Gvir)
      CALL dbt_destroy(t_3c_x_Gocc_2)
      CALL dbt_destroy(t_3c_x_Gvir_2)

      CALL timestop(handle)

   END SUBROUTINE destroy_tensors_chi

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param matrix_index ...
!> \param matrix_name ...
!> \param fm ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE write_matrix(matrix, matrix_index, matrix_name, fm, qs_env)
      TYPE(dbcsr_type)                                   :: matrix
      INTEGER                                            :: matrix_index
      CHARACTER(LEN=*)                                   :: matrix_name
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: fm
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'write_matrix'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_set_all(fm, 0.0_dp)

      CALL copy_dbcsr_to_fm(matrix, fm)

      CALL fm_write(fm, matrix_index, matrix_name, qs_env)

      CALL timestop(handle)

   END SUBROUTINE write_matrix

! **************************************************************************************************
!> \brief ...
!> \param fm ...
!> \param matrix_index ...
!> \param matrix_name ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE fm_write(fm, matrix_index, matrix_name, qs_env)
      TYPE(cp_fm_type)                                   :: fm
      INTEGER                                            :: matrix_index
      CHARACTER(LEN=*)                                   :: matrix_name
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'fm_write'

      CHARACTER(LEN=default_string_length)               :: filename
      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input)

      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN

         IF (matrix_index < 10) THEN
            WRITE (filename, '(3A,I1)') "RESTART_", matrix_name, "_0", matrix_index
         ELSE IF (matrix_index < 100) THEN
            WRITE (filename, '(3A,I2)') "RESTART_", matrix_name, "_", matrix_index
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         unit_nr = cp_print_key_unit_nr(logger, input, key, extension=".matrix", &
                                        file_form="UNFORMATTED", middle_name=TRIM(filename), &
                                        file_position="REWIND", file_action="WRITE")

         CALL cp_fm_write_unformatted(fm, unit_nr)
         IF (unit_nr > 0) THEN
            CALL close_file(unit_nr)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE fm_write

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param tau ...
!> \param fm_G_Gamma ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
! **************************************************************************************************
   SUBROUTINE G_occ_vir(bs_env, tau, fm_G_Gamma, ispin, occ, vir)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: tau
      TYPE(cp_fm_type)                                   :: fm_G_Gamma
      INTEGER                                            :: ispin
      LOGICAL                                            :: occ, vir

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_occ_vir'

      INTEGER                                            :: handle, homo, i_row_local, j_col, &
                                                            j_col_local, n_mo, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp)                                      :: tau_E

      CALL timeset(routineN, handle)

      CPASSERT(occ .NEQV. vir)

      CALL cp_fm_get_info(matrix=bs_env%fm_work_mo(1), &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          col_indices=col_indices)

      n_mo = bs_env%n_ao
      homo = bs_env%n_occ(ispin)

      CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(ispin), bs_env%fm_work_mo(1))

      DO i_row_local = 1, nrow_local
         DO j_col_local = 1, ncol_local

            j_col = col_indices(j_col_local)

            tau_E = ABS(tau*0.5_dp*(bs_env%eigenval_scf_Gamma(j_col, ispin) - bs_env%e_fermi(ispin)))

            IF (tau_E < bs_env%stabilize_exp) THEN
               bs_env%fm_work_mo(1)%local_data(i_row_local, j_col_local) = &
                  bs_env%fm_work_mo(1)%local_data(i_row_local, j_col_local)*EXP(-tau_E)
            ELSE
               bs_env%fm_work_mo(1)%local_data(i_row_local, j_col_local) = 0.0_dp
            END IF

            IF ((occ .AND. j_col > homo) .OR. (vir .AND. j_col <= homo)) THEN
               bs_env%fm_work_mo(1)%local_data(i_row_local, j_col_local) = 0.0_dp
            END IF

         END DO
      END DO

      CALL parallel_gemm(transa="N", transb="T", m=n_mo, n=n_mo, k=n_mo, alpha=1.0_dp, &
                         matrix_a=bs_env%fm_work_mo(1), matrix_b=bs_env%fm_work_mo(1), &
                         beta=0.0_dp, matrix_c=fm_G_Gamma)

      CALL timestop(handle)

   END SUBROUTINE G_occ_vir

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param t_3c ...
!> \param atoms_AO_1 ...
!> \param atoms_AO_2 ...
!> \param atoms_RI ...
! **************************************************************************************************
   SUBROUTINE compute_3c_integrals(qs_env, bs_env, t_3c, atoms_AO_1, atoms_AO_2, atoms_RI)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbt_type)                                     :: t_3c
      INTEGER, DIMENSION(2), OPTIONAL                    :: atoms_AO_1, atoms_AO_2, atoms_RI

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_integrals'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(2)                              :: my_atoms_AO_1, my_atoms_AO_2, my_atoms_RI
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_array

      CALL timeset(routineN, handle)

      ! free memory (not clear whether memory has been freed previously)
      CALL dbt_clear(t_3c)

      ALLOCATE (t_3c_array(1, 1))
      CALL dbt_create(t_3c, t_3c_array(1, 1))

      IF (PRESENT(atoms_AO_1)) THEN
         my_atoms_AO_1 = atoms_AO_1
      ELSE
         my_atoms_AO_1 = [1, bs_env%n_atom]
      END IF
      IF (PRESENT(atoms_AO_2)) THEN
         my_atoms_AO_2 = atoms_AO_2
      ELSE
         my_atoms_AO_2 = [1, bs_env%n_atom]
      END IF
      IF (PRESENT(atoms_RI)) THEN
         my_atoms_RI = atoms_RI
      ELSE
         my_atoms_RI = [1, bs_env%n_atom]
      END IF

      CALL build_3c_integrals(t_3c_array, &
                              bs_env%eps_filter, &
                              qs_env, &
                              bs_env%nl_3c, &
                              int_eps=bs_env%eps_filter, &
                              basis_i=bs_env%basis_set_RI, &
                              basis_j=bs_env%basis_set_AO, &
                              basis_k=bs_env%basis_set_AO, &
                              potential_parameter=bs_env%ri_metric, &
                              bounds_i=atoms_RI, &
                              bounds_j=atoms_AO_1, &
                              bounds_k=atoms_AO_2, &
                              desymmetrize=.FALSE.)

      CALL dbt_copy(t_3c_array(1, 1), t_3c, move_data=.TRUE.)

      CALL dbt_destroy(t_3c_array(1, 1))
      DEALLOCATE (t_3c_array)

      CALL timestop(handle)

   END SUBROUTINE compute_3c_integrals

! **************************************************************************************************
!> \brief ...
!> \param t_3c_for_G ...
!> \param t_G ...
!> \param t_M ...
!> \param bs_env ...
!> \param atoms_AO_1 ...
!> \param atoms_AO_2 ...
!> \param atoms_IL ...
! **************************************************************************************************
   SUBROUTINE G_times_3c(t_3c_for_G, t_G, t_M, bs_env, atoms_AO_1, atoms_AO_2, atoms_IL)
      TYPE(dbt_type)                                     :: t_3c_for_G, t_G, t_M
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(2)                              :: atoms_AO_1, atoms_AO_2, atoms_IL

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_times_3c'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(2)                              :: bounds_IL, bounds_l
      INTEGER, DIMENSION(2, 2)                           :: bounds_k

      CALL timeset(routineN, handle)

      ! JW bounds_IL and bounds_k do not safe any operations, but maybe communication
      !    maybe remove "bounds_1=bounds_IL, &" and "bounds_2=bounds_k, &" later and
      !    check whether performance improves

      bounds_IL(1:2) = [bs_env%i_ao_start_from_atom(atoms_IL(1)), &
                        bs_env%i_ao_end_from_atom(atoms_IL(2))]
      bounds_k(1:2, 1) = [1, bs_env%n_RI]
      bounds_k(1:2, 2) = [bs_env%i_ao_start_from_atom(atoms_AO_2(1)), &
                          bs_env%i_ao_end_from_atom(atoms_AO_2(2))]
      bounds_l(1:2) = [bs_env%i_ao_start_from_atom(atoms_AO_1(1)), &
                       bs_env%i_ao_end_from_atom(atoms_AO_1(2))]

      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=t_3c_for_G, &
                        tensor_2=t_G, &
                        beta=1.0_dp, &
                        tensor_3=t_M, &
                        contract_1=[3], notcontract_1=[1, 2], map_1=[1, 2], &
                        contract_2=[2], notcontract_2=[1], map_2=[3], &
                        bounds_1=bounds_IL, &
                        bounds_2=bounds_k, &
                        bounds_3=bounds_l, &
                        filter_eps=bs_env%eps_filter)

      CALL dbt_clear(t_3c_for_G)

      CALL timestop(handle)

   END SUBROUTINE G_times_3c

! **************************************************************************************************
!> \brief ...
!> \param atoms_1 ...
!> \param atoms_2 ...
!> \param qs_env ...
!> \param bs_env ...
!> \param dist_too_long ...
! **************************************************************************************************
   SUBROUTINE check_dist(atoms_1, atoms_2, qs_env, bs_env, dist_too_long)
      INTEGER, DIMENSION(2)                              :: atoms_1, atoms_2
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL                                            :: dist_too_long

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_dist'

      INTEGER                                            :: atom_1, atom_2, handle
      REAL(dp)                                           :: abs_rab, min_dist_AO_atoms
      REAL(KIND=dp), DIMENSION(3)                        :: rab
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)

      min_dist_AO_atoms = 1.0E5_dp
      DO atom_1 = atoms_1(1), atoms_1(2)
         DO atom_2 = atoms_2(1), atoms_2(2)

            rab = pbc(particle_set(atom_1)%r(1:3), particle_set(atom_2)%r(1:3), cell)

            abs_rab = SQRT(rab(1)**2 + rab(2)**2 + rab(3)**2)

            min_dist_AO_atoms = MIN(min_dist_AO_atoms, abs_rab)

         END DO
      END DO

      dist_too_long = (min_dist_AO_atoms > bs_env%max_dist_AO_atoms)

      CALL timestop(handle)

   END SUBROUTINE check_dist

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE get_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_W_MIC'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (bs_env%all_W_exist) THEN
         CALL read_W_MIC_time(bs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      ELSE
         CALL compute_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V_kp ...
!> \param ikp_batch ...
! **************************************************************************************************
   SUBROUTINE compute_V_k_by_lattice_sum(bs_env, qs_env, fm_V_kp, ikp_batch)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_kp
      INTEGER                                            :: ikp_batch

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_k_by_lattice_sum'

      INTEGER                                            :: handle, ikp, ikp_end, ikp_start, &
                                                            nkp_chi_eps_W_batch, re_im
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      nkp_chi_eps_W_batch = bs_env%nkp_chi_eps_W_batch

      ikp_start = (ikp_batch - 1)*bs_env%nkp_chi_eps_W_batch + 1
      ikp_end = MIN(ikp_batch*bs_env%nkp_chi_eps_W_batch, bs_env%kpoints_chi_eps_W%nkp)

      NULLIFY (mat_V_kp)
      ALLOCATE (mat_V_kp(ikp_start:ikp_end, 2))

      DO ikp = ikp_start, ikp_end
         DO re_im = 1, 2
            NULLIFY (mat_V_kp(ikp, re_im)%matrix)
            ALLOCATE (mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_create(mat_V_kp(ikp, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
            CALL dbcsr_reserve_all_blocks(mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_set(mat_V_kp(ikp, re_im)%matrix, 0.0_dp)

         END DO ! re_im
      END DO ! ikp

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (ikp_end .LE. bs_env%nkp_chi_eps_W_orig) THEN

         ! 1. 2c Coulomb integrals for the first "original" k-point grid
         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ELSE IF (ikp_start > bs_env%nkp_chi_eps_W_orig .AND. &
               ikp_end .LE. bs_env%nkp_chi_eps_W_orig_plus_extra) THEN

         ! 2. 2c Coulomb integrals for the second "extrapolation" k-point grid
         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_extra

      ELSE

         CPABORT("Error with k-point parallelization.")

      END IF

      CALL build_2c_coulomb_matrix_kp(mat_V_kp, &
                                      bs_env%kpoints_chi_eps_W, &
                                      basis_type="RI_AUX", &
                                      cell=cell, &
                                      particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=ikp_start, &
                                      ikp_end=ikp_end)

      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ALLOCATE (fm_V_kp(ikp_start:ikp_end, 2))
      DO ikp = ikp_start, ikp_end
         DO re_im = 1, 2
            CALL cp_fm_create(fm_V_kp(ikp, re_im), bs_env%fm_RI_RI%matrix_struct)
            CALL copy_dbcsr_to_fm(mat_V_kp(ikp, re_im)%matrix, fm_V_kp(ikp, re_im))
            CALL dbcsr_deallocate_matrix(mat_V_kp(ikp, re_im)%matrix)
         END DO
      END DO
      DEALLOCATE (mat_V_kp)

      CALL timestop(handle)

   END SUBROUTINE compute_V_k_by_lattice_sum

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V_kp ...
!> \param cfm_V_sqrt_ikp ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE compute_MinvVsqrt_Vsqrt(bs_env, qs_env, fm_V_kp, cfm_V_sqrt_ikp, &
                                      cfm_M_inv_V_sqrt_ikp, ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_kp
      TYPE(cp_cfm_type)                                  :: cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_MinvVsqrt_Vsqrt'

      INTEGER                                            :: handle, info, n_RI
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_ikp, cfm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M_ikp

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      ! get here M(k) and write it to fm_M_ikp
      CALL RI_2c_integral_mat(qs_env, fm_M_ikp, fm_V_kp(ikp, 1), &
                              n_RI, bs_env%ri_metric, do_kpoints=.TRUE., &
                              kpoints=bs_env%kpoints_chi_eps_W, &
                              regularization_RI=bs_env%regularization_RI, ikp_ext=ikp, &
                              do_build_cell_index=(ikp == 1))

      IF (ikp == 1) THEN
         CALL cp_cfm_create(cfm_V_sqrt_ikp, fm_V_kp(ikp, 1)%matrix_struct)
         CALL cp_cfm_create(cfm_M_inv_V_sqrt_ikp, fm_V_kp(ikp, 1)%matrix_struct)
      END IF
      CALL cp_cfm_create(cfm_M_inv_ikp, fm_V_kp(ikp, 1)%matrix_struct)

      CALL cp_fm_to_cfm(fm_M_ikp(1, 1), fm_M_ikp(1, 2), cfm_M_inv_ikp)
      CALL cp_fm_to_cfm(fm_V_kp(ikp, 1), fm_V_kp(ikp, 2), cfm_V_sqrt_ikp)

      CALL cp_fm_release(fm_M_ikp)

      CALL cp_cfm_create(cfm_work, fm_V_kp(ikp, 1)%matrix_struct)

      ! M(k) -> M^-1(k)
      CALL cp_cfm_to_cfm(cfm_M_inv_ikp, cfm_work)
      CALL cp_cfm_cholesky_decompose(matrix=cfm_M_inv_ikp, n=n_RI, info_out=info)
      IF (info == 0) THEN
         ! successful Cholesky decomposition
         CALL cp_cfm_cholesky_invert(cfm_M_inv_ikp)
         ! symmetrize the result
         CALL cp_cfm_upper_to_full(cfm_M_inv_ikp)
      ELSE
         ! Cholesky decomposition not successful: use expensive diagonalization
         CALL cp_cfm_power(cfm_work, threshold=bs_env%eps_eigval_mat_RI, exponent=-1.0_dp)
         CALL cp_cfm_to_cfm(cfm_work, cfm_M_inv_ikp)
      END IF

      ! V(k) -> L(k) with L^H(k)*L(k) = V(k) [L(k) can be just considered to be V^0.5(k)]
      CALL cp_cfm_to_cfm(cfm_V_sqrt_ikp, cfm_work)
      CALL cp_cfm_cholesky_decompose(matrix=cfm_V_sqrt_ikp, n=n_RI, info_out=info)
      IF (info == 0) THEN
         ! successful Cholesky decomposition
         CALL clean_lower_part(cfm_V_sqrt_ikp)
      ELSE
         ! Cholesky decomposition not successful: use expensive diagonalization
         CALL cp_cfm_power(cfm_work, threshold=0.0_dp, exponent=0.5_dp)
         CALL cp_cfm_to_cfm(cfm_work, cfm_V_sqrt_ikp)
      END IF
      CALL cp_cfm_release(cfm_work)

      ! get M^-1(k)*V^0.5(k)
      CALL parallel_gemm("N", "C", n_RI, n_RI, n_RI, z_one, cfm_M_inv_ikp, cfm_V_sqrt_ikp, &
                         z_zero, cfm_M_inv_V_sqrt_ikp)

      CALL cp_cfm_release(cfm_M_inv_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_MinvVsqrt_Vsqrt

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE read_W_MIC_time(bs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_W_MIC_time'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)
      CALL create_fm_W_MIC_time(bs_env, fm_W_MIC_time)

      DO i_t = 1, bs_env%num_time_freq_points

         t1 = m_walltime()

         CALL fm_read(fm_W_MIC_time(i_t), bs_env, bs_env%W_time_name, i_t)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I5,A,I3,A,F7.1,A)') &
               'Read W^MIC(iτ) from file for time point  ', i_t, ' /', bs_env%num_time_freq_points, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE read_W_MIC_time

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE compute_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W_MIC'

      INTEGER                                            :: handle, i_t, ikp, ikp_batch, &
                                                            ikp_in_batch, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_kp

      CALL timeset(routineN, handle)

      CALL create_fm_W_MIC_time(bs_env, fm_W_MIC_time)

      DO ikp_batch = 1, bs_env%num_chi_eps_W_batches

         t1 = m_walltime()

         ! Compute V_PQ(k) = sum_R e^(ikR) <phi_P, cell 0 | 1/r | phi_Q, cell R>
         CALL compute_V_k_by_lattice_sum(bs_env, qs_env, fm_V_kp, ikp_batch)

         DO ikp_in_batch = 1, bs_env%nkp_chi_eps_W_batch

            ikp = (ikp_batch - 1)*bs_env%nkp_chi_eps_W_batch + ikp_in_batch

            IF (ikp > bs_env%nkp_chi_eps_W_orig_plus_extra) CYCLE

            CALL compute_MinvVsqrt_Vsqrt(bs_env, qs_env, fm_V_kp, &
                                         cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp, ikp)

            CALL bs_env%para_env%sync()

            DO j_w = 1, bs_env%num_time_freq_points

               ! check if we need this (ikp, ω_j) combination for approximate k-point extrapolation
               IF (bs_env%approx_kp_extrapol .AND. j_w > 1 .AND. &
                   ikp > bs_env%nkp_chi_eps_W_orig) CYCLE

               CALL compute_fm_W_MIC_freq_j(bs_env, qs_env, bs_env%fm_W_MIC_freq, j_w, ikp, &
                                            mat_chi_Gamma_tau, cfm_M_inv_V_sqrt_ikp, &
                                            cfm_V_sqrt_ikp)

               ! Fourier trafo from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
               CALL Fourier_transform_w_to_t(bs_env, fm_W_MIC_time, bs_env%fm_W_MIC_freq, j_w)

            END DO ! ω_j

            CALL cp_fm_release(fm_V_kp(ikp, 1))
            CALL cp_fm_release(fm_V_kp(ikp, 2))

         END DO ! ikp_in_batch

         DEALLOCATE (fm_V_kp)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I12,A,I3,A,F7.1,A)') &
               'Computed W(iτ,k) for k-point batch', &
               ikp_batch, ' /', bs_env%num_chi_eps_W_batches, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! ikp_batch

      IF (bs_env%approx_kp_extrapol) THEN
         CALL apply_extrapol_factor(bs_env, fm_W_MIC_time)
      END IF

      ! M^-1(k=0)*W^MIC(iτ)*M^-1(k=0)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_W_MIC_time)

      DO i_t = 1, bs_env%num_time_freq_points
         CALL fm_write(fm_W_MIC_time(i_t), i_t, bs_env%W_time_name, qs_env)
      END DO

      CALL cp_cfm_release(cfm_M_inv_V_sqrt_ikp)
      CALL cp_cfm_release(cfm_V_sqrt_ikp)
      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
!> \param ikp ...
!> \param mat_chi_Gamma_tau ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param cfm_V_sqrt_ikp ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_MIC_freq_j(bs_env, qs_env, fm_W_MIC_freq_j, j_w, ikp, mat_chi_Gamma_tau, &
                                      cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w, ikp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_W_MIC_freq_j'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_chi_ikp_freq_j, cfm_W_ikp_freq_j

      CALL timeset(routineN, handle)

      ! 1. Fourier transformation of χ_PQ(iτ,k=0) to χ_PQ(iω_j,k=0)
      CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

      CALL cp_fm_set_all(fm_W_MIC_freq_j, 0.0_dp)

      ! 2. Get χ_PQ(iω_j,k_i) from χ_PQ(iω_j,k=0) using the minimum image convention
      CALL cfm_ikp_from_fm_Gamma(cfm_chi_ikp_freq_j, bs_env%fm_chi_Gamma_freq, &
                                 ikp, qs_env, bs_env%kpoints_chi_eps_W, "RI_AUX")

      ! 3. Remove all negative eigenvalues from χ_PQ(iω_j,k_i)
      CALL cp_cfm_power(cfm_chi_ikp_freq_j, threshold=0.0_dp, exponent=1.0_dp)

      ! 4. ε(iω_j,k_i) = Id - V^0.5(k_i)*M^-1(k_i)*χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)
      !    W(iω_j,k_i) = V^0.5(k_i)*(ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)
      CALL compute_cfm_W_ikp_freq_j(bs_env, cfm_chi_ikp_freq_j, cfm_V_sqrt_ikp, &
                                    cfm_M_inv_V_sqrt_ikp, cfm_W_ikp_freq_j)

      ! 5. k-point integration W_PQ(iω_j, k_i) to W_PQ^MIC(iω_j)
      SELECT CASE (bs_env%approx_kp_extrapol)
      CASE (.FALSE.)
         ! default: standard k-point extrapolation
         CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, cfm_W_ikp_freq_j, ikp, &
                                        bs_env%kpoints_chi_eps_W, "RI_AUX")
      CASE (.TRUE.)
         ! for approximate kpoint extrapolation: get W_PQ^MIC(iω_1) with and without k-point
         ! extrapolation to compute the extrapolation factor f_PQ for every PQ-matrix element,
         ! f_PQ = (W_PQ^MIC(iω_1) with extrapolation) / (W_PQ^MIC(iω_1) without extrapolation)

         ! for ω_1, we compute the k-point extrapolated result using all k-points
         IF (j_w == 1) THEN

            ! k-point extrapolated
            CALL MIC_contribution_from_ikp(bs_env, qs_env, bs_env%fm_W_MIC_freq_1_extra, &
                                           cfm_W_ikp_freq_j, ikp, bs_env%kpoints_chi_eps_W, &
                                           "RI_AUX")
            ! non-kpoint extrapolated
            IF (ikp .LE. bs_env%nkp_chi_eps_W_orig) THEN
               CALL MIC_contribution_from_ikp(bs_env, qs_env, bs_env%fm_W_MIC_freq_1_no_extra, &
                                              cfm_W_ikp_freq_j, ikp, bs_env%kpoints_chi_eps_W, &
                                              "RI_AUX", wkp_ext=bs_env%wkp_orig)
            END IF

         END IF

         ! for all ω_j, we need to compute W^MIC without k-point extrpolation
         IF (ikp .LE. bs_env%nkp_chi_eps_W_orig) THEN
            CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, cfm_W_ikp_freq_j, &
                                           ikp, bs_env%kpoints_chi_eps_W, "RI_AUX", &
                                           wkp_ext=bs_env%wkp_orig)
         END IF
      END SELECT

      CALL cp_cfm_release(cfm_W_ikp_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_MIC_freq_j

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat ...
! **************************************************************************************************
   SUBROUTINE clean_lower_part(cfm_mat)
      TYPE(cp_cfm_type)                                  :: cfm_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'clean_lower_part'

      INTEGER                                            :: handle, i_global, i_row, j_col, &
                                                            j_global, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=cfm_mat, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      DO i_row = 1, nrow_local
         DO j_col = 1, ncol_local
            i_global = row_indices(i_row)
            j_global = col_indices(j_col)
            IF (j_global < i_global) cfm_mat%local_data(i_row, j_col) = z_zero
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE clean_lower_part

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE apply_extrapol_factor(bs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_extrapol_factor'

      INTEGER                                            :: handle, i, i_t, j, ncol_local, nrow_local
      REAL(KIND=dp)                                      :: extrapol_factor, W_extra_1, W_no_extra_1

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_W_MIC_time(1), nrow_local=nrow_local, ncol_local=ncol_local)

      DO i_t = 1, bs_env%num_time_freq_points
         DO i = 1, nrow_local
            DO j = 1, ncol_local

               W_extra_1 = bs_env%fm_W_MIC_freq_1_extra%local_data(i, j)
               W_no_extra_1 = bs_env%fm_W_MIC_freq_1_no_extra%local_data(i, j)

               IF (ABS(W_no_extra_1) > 1.0E-13) THEN
                  extrapol_factor = W_extra_1/W_no_extra_1
               ELSE
                  extrapol_factor = 1.0_dp
               END IF

               ! reset extrapolation factor if it is very large
               IF (ABS(extrapol_factor) > 10.0_dp) extrapol_factor = 1.0_dp

               fm_W_MIC_time(i_t)%local_data(i, j) = fm_W_MIC_time(i_t)%local_data(i, j) &
                                                     *extrapol_factor
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE apply_extrapol_factor

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_chi_Gamma_freq ...
!> \param j_w ...
!> \param mat_chi_Gamma_tau ...
! **************************************************************************************************
   SUBROUTINE compute_fm_chi_Gamma_freq(bs_env, fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type)                                   :: fm_chi_Gamma_freq
      INTEGER                                            :: j_w
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_chi_Gamma_freq'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: freq_j, time_i, weight_ij

      CALL timeset(routineN, handle)

      CALL dbcsr_set(bs_env%mat_RI_RI%matrix, 0.0_dp)

      freq_j = bs_env%imag_freq_points(j_w)

      DO i_t = 1, bs_env%num_time_freq_points

         time_i = bs_env%imag_time_points(i_t)
         weight_ij = bs_env%weights_cos_t_to_w(j_w, i_t)

         ! actual Fourier transform
         CALL dbcsr_add(bs_env%mat_RI_RI%matrix, mat_chi_Gamma_tau(i_t)%matrix, &
                        1.0_dp, COS(time_i*freq_j)*weight_ij)

      END DO

      CALL copy_dbcsr_to_fm(bs_env%mat_RI_RI%matrix, fm_chi_Gamma_freq)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_chi_Gamma_freq

! **************************************************************************************************
!> \brief ...
!> \param mat_ikp_re ...
!> \param mat_ikp_im ...
!> \param mat_Gamma ...
!> \param kpoints ...
!> \param ikp ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE mat_ikp_from_mat_Gamma(mat_ikp_re, mat_ikp_im, mat_Gamma, kpoints, ikp, qs_env)
      TYPE(dbcsr_type)                                   :: mat_ikp_re, mat_ikp_im, mat_Gamma
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mat_ikp_from_mat_Gamma'

      INTEGER                                            :: col, handle, i_cell, j_cell, num_cells, &
                                                            row
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL :: f, i_cell_is_the_minimum_image_cell
      REAL(KIND=dp)                                      :: abs_rab_cell_i, abs_rab_cell_j, arg
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, cell_vector_j, rab_cell_i, &
                                                            rab_cell_j
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_im, block_re, data_block
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      ! get the same blocks in mat_ikp_re and mat_ikp_im as in mat_Gamma
      CALL dbcsr_copy(mat_ikp_re, mat_Gamma)
      CALL dbcsr_copy(mat_ikp_im, mat_Gamma)
      CALL dbcsr_set(mat_ikp_re, 0.0_dp)
      CALL dbcsr_set(mat_ikp_im, 0.0_dp)

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      index_to_cell => kpoints%index_to_cell

      num_cells = SIZE(index_to_cell, 2)

      DO i_cell = 1, num_cells

         CALL dbcsr_iterator_start(iter, mat_Gamma)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_block)

            cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, i_cell), dp))

            rab_cell_i(1:3) = pbc(particle_set(row)%r(1:3), cell) - &
                              (pbc(particle_set(col)%r(1:3), cell) + cell_vector(1:3))
            abs_rab_cell_i = SQRT(rab_cell_i(1)**2 + rab_cell_i(2)**2 + rab_cell_i(3)**2)

            ! minimum image convention
            i_cell_is_the_minimum_image_cell = .TRUE.
            DO j_cell = 1, num_cells
               cell_vector_j(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, j_cell), dp))
               rab_cell_j(1:3) = pbc(particle_set(row)%r(1:3), cell) - &
                                 (pbc(particle_set(col)%r(1:3), cell) + cell_vector_j(1:3))
               abs_rab_cell_j = SQRT(rab_cell_j(1)**2 + rab_cell_j(2)**2 + rab_cell_j(3)**2)

               IF (abs_rab_cell_i > abs_rab_cell_j + 1.0E-6_dp) THEN
                  i_cell_is_the_minimum_image_cell = .FALSE.
               END IF
            END DO

            IF (i_cell_is_the_minimum_image_cell) THEN
               NULLIFY (block_re, block_im)
               CALL dbcsr_get_block_p(matrix=mat_ikp_re, row=row, col=col, block=block_re, found=f)
               CALL dbcsr_get_block_p(matrix=mat_ikp_im, row=row, col=col, block=block_im, found=f)
               CPASSERT(ALL(ABS(block_re) < 1.0E-10_dp))
               CPASSERT(ALL(ABS(block_im) < 1.0E-10_dp))

               arg = REAL(index_to_cell(1, i_cell), dp)*kpoints%xkp(1, ikp) + &
                     REAL(index_to_cell(2, i_cell), dp)*kpoints%xkp(2, ikp) + &
                     REAL(index_to_cell(3, i_cell), dp)*kpoints%xkp(3, ikp)

               block_re(:, :) = COS(twopi*arg)*data_block(:, :)
               block_im(:, :) = SIN(twopi*arg)*data_block(:, :)
            END IF

         END DO
         CALL dbcsr_iterator_stop(iter)

      END DO

      CALL timestop(handle)

   END SUBROUTINE mat_ikp_from_mat_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param cfm_chi_ikp_freq_j ...
!> \param cfm_V_sqrt_ikp ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param cfm_W_ikp_freq_j ...
! **************************************************************************************************
   SUBROUTINE compute_cfm_W_ikp_freq_j(bs_env, cfm_chi_ikp_freq_j, cfm_V_sqrt_ikp, &
                                       cfm_M_inv_V_sqrt_ikp, cfm_W_ikp_freq_j)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_chi_ikp_freq_j, cfm_V_sqrt_ikp, &
                                                            cfm_M_inv_V_sqrt_ikp, cfm_W_ikp_freq_j

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cfm_W_ikp_freq_j'

      INTEGER                                            :: handle, info, n_RI
      TYPE(cp_cfm_type)                                  :: cfm_eps_ikp_freq_j, cfm_work

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_work, cfm_chi_ikp_freq_j%matrix_struct)
      n_RI = bs_env%n_RI

      ! 1. ε(iω_j,k) = Id - V^0.5(k)*M^-1(k)*χ(iω_j,k)*M^-1(k)*V^0.5(k)

      ! 1. a) work = χ(iω_j,k)*M^-1(k)*V^0.5(k)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, z_one, &
                         cfm_chi_ikp_freq_j, cfm_M_inv_V_sqrt_ikp, z_zero, cfm_work)
      CALL cp_cfm_release(cfm_chi_ikp_freq_j)

      ! 1. b) eps_work = V^0.5(k)*M^-1(k)*work
      CALL cp_cfm_create(cfm_eps_ikp_freq_j, cfm_work%matrix_struct)
      CALL parallel_gemm('C', 'N', n_RI, n_RI, n_RI, z_one, &
                         cfm_M_inv_V_sqrt_ikp, cfm_work, z_zero, cfm_eps_ikp_freq_j)

      ! 1. c) ε(iω_j,k) = eps_work - Id
      CALL cfm_add_on_diag(cfm_eps_ikp_freq_j, z_one)

      ! 2. W(iω_j,k) = V^0.5(k)*(ε^-1(iω_j,k)-Id)*V^0.5(k)

      ! 2. a) Cholesky decomposition of ε(iω_j,k) as preparation for inversion
      CALL cp_cfm_cholesky_decompose(matrix=cfm_eps_ikp_freq_j, n=n_RI, info_out=info)
      CPASSERT(info == 0)

      ! 2. b) Inversion of ε(iω_j,k) using its Cholesky decomposition
      CALL cp_cfm_cholesky_invert(cfm_eps_ikp_freq_j)
      CALL cp_cfm_upper_to_full(cfm_eps_ikp_freq_j)

      ! 2. c) ε^-1(iω_j,k)-Id
      CALL cfm_add_on_diag(cfm_eps_ikp_freq_j, -z_one)

      ! 2. d) work = (ε^-1(iω_j,k)-Id)*V^0.5(k)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, z_one, cfm_eps_ikp_freq_j, cfm_V_sqrt_ikp, &
                         z_zero, cfm_work)

      ! 2. e) W(iw,k) = V^0.5(k)*work
      CALL cp_cfm_create(cfm_W_ikp_freq_j, cfm_work%matrix_struct)
      CALL parallel_gemm('C', 'N', n_RI, n_RI, n_RI, z_one, cfm_V_sqrt_ikp, cfm_work, &
                         z_zero, cfm_W_ikp_freq_j)

      CALL cp_cfm_release(cfm_work)
      CALL cp_cfm_release(cfm_eps_ikp_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_cfm_W_ikp_freq_j

! **************************************************************************************************
!> \brief ...
!> \param cfm ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE cfm_add_on_diag(cfm, alpha)

      TYPE(cp_cfm_type)                                  :: cfm
      COMPLEX(KIND=dp)                                   :: alpha

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cfm_add_on_diag'

      INTEGER                                            :: handle, i_global, i_row, j_col, &
                                                            j_global, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=cfm, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ! add 1 on the diagonal
      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            i_global = row_indices(i_row)
            IF (j_global == i_global) THEN
               cfm%local_data(i_row, j_col) = cfm%local_data(i_row, j_col) + alpha
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE cfm_add_on_diag

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE create_fm_W_MIC_time(bs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fm_W_MIC_time'

      INTEGER                                            :: handle, i_t

      CALL timeset(routineN, handle)

      ALLOCATE (fm_W_MIC_time(bs_env%num_time_freq_points))
      DO i_t = 1, bs_env%num_time_freq_points
         CALL cp_fm_create(fm_W_MIC_time(i_t), bs_env%fm_RI_RI%matrix_struct)
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_fm_W_MIC_time

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
! **************************************************************************************************
   SUBROUTINE Fourier_transform_w_to_t(bs_env, fm_W_MIC_time, fm_W_MIC_freq_j, j_w)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Fourier_transform_w_to_t'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: freq_j, time_i, weight_ij

      CALL timeset(routineN, handle)

      freq_j = bs_env%imag_freq_points(j_w)

      DO i_t = 1, bs_env%num_time_freq_points

         time_i = bs_env%imag_time_points(i_t)
         weight_ij = bs_env%weights_cos_w_to_t(i_t, j_w)

         ! actual Fourier transform
         CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_W_MIC_time(i_t), &
                                  beta=weight_ij*COS(time_i*freq_j), matrix_b=fm_W_MIC_freq_j)

      END DO

      CALL timestop(handle)

   END SUBROUTINE Fourier_transform_w_to_t

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'multiply_fm_W_MIC_time_with_Minv_Gamma'

      INTEGER                                            :: handle, i_t, n_RI, ndep
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Minv_Gamma

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_work, fm_W_MIC_time(1)%matrix_struct)

      ! compute Gamma-only RI-metric matrix M(k=0); no regularization
      CALL RI_2c_integral_mat(qs_env, fm_Minv_Gamma, fm_W_MIC_time(1), n_RI, &
                              bs_env%ri_metric, do_kpoints=.FALSE.)

      CALL cp_fm_power(fm_Minv_Gamma(1, 1), fm_work, -1.0_dp, 0.0_dp, ndep)

      ! M^-1(k=0)*W^MIC(iτ)*M^-1(k=0)
      DO i_t = 1, SIZE(fm_W_MIC_time)

         CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_Minv_Gamma(1, 1), &
                            fm_W_MIC_time(i_t), 0.0_dp, fm_work)

         CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_work, &
                            fm_Minv_Gamma(1, 1), 0.0_dp, fm_W_MIC_time(i_t))

      END DO

      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_Minv_Gamma)

      CALL timestop(handle)

   END SUBROUTINE multiply_fm_W_MIC_time_with_Minv_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE get_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_Sigma_x'

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      IF (bs_env%Sigma_x_exists) THEN
         DO ispin = 1, bs_env%n_spin
            CALL fm_read(fm_Sigma_x_Gamma(ispin), bs_env, bs_env%Sigma_x_name, ispin)
         END DO
      ELSE
         CALL compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_Sigma_x

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, i_intval_idx, ispin, j_intval_idx
      INTEGER, DIMENSION(2)                              :: i_atoms, j_atoms
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma
      TYPE(dbt_type)                                     :: t_2c_D, t_2c_Sigma_x, t_2c_V, t_3c_x_V

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL dbt_create(bs_env%t_G, t_2c_D)
      CALL dbt_create(bs_env%t_W, t_2c_V)
      CALL dbt_create(bs_env%t_G, t_2c_Sigma_x)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_V)
      CALL dbcsr_create(mat_Sigma_x_Gamma, template=bs_env%mat_ao_ao%matrix)

      ! 1. Compute truncated Coulomb operator matrix V^tr(k=0) (cutoff rad: cellsize/2)
      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI, bs_env%n_RI, &
                              bs_env%trunc_coulomb, do_kpoints=.FALSE.)

      ! 2. Compute M^-1(k=0) and get M^-1(k=0)*V^tr(k=0)*M^-1(k=0)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(bs_env, qs_env, fm_Vtr_Gamma(:, 1))

      DO ispin = 1, bs_env%n_spin

         ! 3. Compute density matrix D_µν
         CALL G_occ_vir(bs_env, 0.0_dp, bs_env%fm_work_mo(2), ispin, occ=.TRUE., vir=.FALSE.)

         CALL fm_to_local_tensor(bs_env%fm_work_mo(2), bs_env%mat_ao_ao%matrix, &
                                 bs_env%mat_ao_ao_tensor%matrix, t_2c_D, bs_env, &
                                 bs_env%atoms_i_t_group)

         CALL fm_to_local_tensor(fm_Vtr_Gamma(1, 1), bs_env%mat_RI_RI%matrix, &
                                 bs_env%mat_RI_RI_tensor%matrix, t_2c_V, bs_env, &
                                 bs_env%atoms_j_t_group)

         ! every group has its own range of i_atoms and j_atoms; only deal with a
         ! limited number of i_atom-j_atom pairs simultaneously in a group to save memory
         DO i_intval_idx = 1, bs_env%n_intervals_i
            DO j_intval_idx = 1, bs_env%n_intervals_j
               i_atoms = bs_env%i_atom_intervals(1:2, i_intval_idx)
               j_atoms = bs_env%j_atom_intervals(1:2, j_intval_idx)

               ! 4. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
               ! 5. M_νσQ(iτ) = sum_P (νσ|P) (M^-1(k=0)*V^tr(k=0)*M^-1(k=0))_PQ(iτ)
               CALL compute_3c_and_contract_W(qs_env, bs_env, i_atoms, j_atoms, t_3c_x_V, t_2c_V)

               ! 6. tensor operations with D and computation of Σ^x
               !    Σ^x_λσ(k=0) = sum_νQ M_νσQ(iτ) sum_µ (λµ|Q) D_µν
               CALL contract_to_Sigma(t_2c_D, t_3c_x_V, t_2c_Sigma_x, i_atoms, j_atoms, &
                                      qs_env, bs_env, occ=.TRUE., vir=.FALSE., clear_W=.TRUE.)

            END DO ! j_atoms
         END DO ! i_atoms

         CALL local_dbt_to_global_mat(t_2c_Sigma_x, bs_env%mat_ao_ao_tensor%matrix, &
                                      mat_Sigma_x_Gamma, bs_env%para_env)

         CALL write_matrix(mat_Sigma_x_Gamma, ispin, bs_env%Sigma_x_name, &
                           bs_env%fm_work_mo(1), qs_env)

         CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma(ispin))

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL dbcsr_release(mat_Sigma_x_Gamma)
      CALL dbt_destroy(t_2c_D)
      CALL dbt_destroy(t_2c_V)
      CALL dbt_destroy(t_2c_Sigma_x)
      CALL dbt_destroy(t_3c_x_V)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_time ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE get_Sigma_c(bs_env, qs_env, fm_W_MIC_time, fm_Sigma_c_Gamma_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_Sigma_c'

      INTEGER                                            :: handle, i_intval_idx, i_t, ispin, &
                                                            j_intval_idx, read_write_index
      INTEGER, DIMENSION(2)                              :: i_atoms, j_atoms
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, &
                                                            t_2c_Sigma_neg_tau, &
                                                            t_2c_Sigma_pos_tau, t_2c_W, t_3c_x_W

      CALL timeset(routineN, handle)

      CALL create_mat_for_Sigma_c(bs_env, t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                                  t_2c_Sigma_pos_tau, t_3c_x_W, &
                                  mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      DO i_t = 1, bs_env%num_time_freq_points

         DO ispin = 1, bs_env%n_spin

            t1 = m_walltime()

            read_write_index = i_t + (ispin - 1)*bs_env%num_time_freq_points

            ! read self-energy from restart
            IF (bs_env%Sigma_c_exists(i_t, ispin)) THEN
               CALL fm_read(bs_env%fm_work_mo(1), bs_env, bs_env%Sigma_p_name, read_write_index)
               CALL copy_fm_to_dbcsr(bs_env%fm_work_mo(1), mat_Sigma_pos_tau(i_t, ispin)%matrix, &
                                     keep_sparsity=.FALSE.)
               CALL fm_read(bs_env%fm_work_mo(1), bs_env, bs_env%Sigma_n_name, read_write_index)
               CALL copy_fm_to_dbcsr(bs_env%fm_work_mo(1), mat_Sigma_neg_tau(i_t, ispin)%matrix, &
                                     keep_sparsity=.FALSE.)
               IF (bs_env%unit_nr > 0) THEN
                  WRITE (bs_env%unit_nr, '(T2,2A,I3,A,I3,A,F7.1,A)') 'Read Σ^c(iτ,k=0) ', &
                     'from file for time point  ', i_t, ' /', bs_env%num_time_freq_points, &
                     ',    Execution time', m_walltime() - t1, ' s'
               END IF

               CYCLE

            END IF

            tau = bs_env%imag_time_points(i_t)

            CALL G_occ_vir(bs_env, tau, bs_env%fm_Gocc, ispin, occ=.TRUE., vir=.FALSE.)
            CALL G_occ_vir(bs_env, tau, bs_env%fm_Gvir, ispin, occ=.FALSE., vir=.TRUE.)

            ! fm G^occ, G^vir and W to local tensor
            CALL fm_to_local_tensor(bs_env%fm_Gocc, bs_env%mat_ao_ao%matrix, &
                                    bs_env%mat_ao_ao_tensor%matrix, t_2c_Gocc, bs_env, &
                                    bs_env%atoms_i_t_group)
            CALL fm_to_local_tensor(bs_env%fm_Gvir, bs_env%mat_ao_ao%matrix, &
                                    bs_env%mat_ao_ao_tensor%matrix, t_2c_Gvir, bs_env, &
                                    bs_env%atoms_i_t_group)
            CALL fm_to_local_tensor(fm_W_MIC_time(i_t), bs_env%mat_RI_RI%matrix, &
                                    bs_env%mat_RI_RI_tensor%matrix, t_2c_W, bs_env, &
                                    bs_env%atoms_j_t_group)

            ! every group has its own range of i_atoms and j_atoms; only deal with a
            ! limited number of i_atom-j_atom pairs simultaneously in a group to save memory
            DO i_intval_idx = 1, bs_env%n_intervals_i
               DO j_intval_idx = 1, bs_env%n_intervals_j
                  i_atoms = bs_env%i_atom_intervals(1:2, i_intval_idx)
                  j_atoms = bs_env%j_atom_intervals(1:2, j_intval_idx)

                  IF (bs_env%skip_Sigma_occ(i_intval_idx, j_intval_idx) .AND. &
                      bs_env%skip_Sigma_vir(i_intval_idx, j_intval_idx)) CYCLE

                  ! 1. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
                  ! 2. tensor operation M_νσQ(iτ) = sum_P (νσ|P) W^MIC_PQ(iτ)
                  CALL compute_3c_and_contract_W(qs_env, bs_env, i_atoms, j_atoms, t_3c_x_W, t_2c_W)

                  ! 3. Σ_λσ(iτ,k=0) = sum_νQ M_νσQ(iτ) sum_µ (λµ|Q) G^occ_µν(i|τ|) for τ < 0
                  !    (recall M_νσQ(iτ) = M_νσQ(-iτ) because W^MIC_PQ(iτ) = W^MIC_PQ(-iτ) )
                  CALL contract_to_Sigma(t_2c_Gocc, t_3c_x_W, t_2c_Sigma_neg_tau, i_atoms, j_atoms, &
                                         qs_env, bs_env, occ=.TRUE., vir=.FALSE., clear_W=.FALSE., &
                                         can_skip=bs_env%skip_Sigma_occ(i_intval_idx, j_intval_idx))

                  !    Σ_λσ(iτ,k=0) = sum_νQ M_νσQ(iτ) sum_µ (λµ|Q) G^vir_µν(i|τ|) for τ > 0
                  CALL contract_to_Sigma(t_2c_Gvir, t_3c_x_W, t_2c_Sigma_pos_tau, i_atoms, j_atoms, &
                                         qs_env, bs_env, occ=.FALSE., vir=.TRUE., clear_W=.TRUE., &
                                         can_skip=bs_env%skip_Sigma_vir(i_intval_idx, j_intval_idx))

               END DO ! j_atoms
            END DO ! i_atoms

            ! 4. communicate data tensor t_2c_Sigma (which is local in the subgroup)
            !    to the global dbcsr matrix mat_Sigma_pos/neg_tau (which stores Σ for all iτ)
            CALL local_dbt_to_global_mat(t_2c_Sigma_neg_tau, bs_env%mat_ao_ao_tensor%matrix, &
                                         mat_Sigma_neg_tau(i_t, ispin)%matrix, bs_env%para_env)
            CALL local_dbt_to_global_mat(t_2c_Sigma_pos_tau, bs_env%mat_ao_ao_tensor%matrix, &
                                         mat_Sigma_pos_tau(i_t, ispin)%matrix, bs_env%para_env)

            CALL write_matrix(mat_Sigma_pos_tau(i_t, ispin)%matrix, read_write_index, &
                              bs_env%Sigma_p_name, bs_env%fm_work_mo(1), qs_env)
            CALL write_matrix(mat_Sigma_neg_tau(i_t, ispin)%matrix, read_write_index, &
                              bs_env%Sigma_n_name, bs_env%fm_work_mo(1), qs_env)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I10,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ,k=0) for time point ', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL print_skipping(bs_env)

      CALL destroy_mat_Sigma_c(t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                               t_2c_Sigma_pos_tau, t_3c_x_W, fm_W_MIC_time, &
                               mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      CALL delete_unnecessary_files(bs_env)

      CALL timestop(handle)

   END SUBROUTINE get_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param t_2c_Gocc ...
!> \param t_2c_Gvir ...
!> \param t_2c_W ...
!> \param t_2c_Sigma_neg_tau ...
!> \param t_2c_Sigma_pos_tau ...
!> \param t_3c_x_W ...
!> \param mat_Sigma_neg_tau ...
!> \param mat_Sigma_pos_tau ...
! **************************************************************************************************
   SUBROUTINE create_mat_for_Sigma_c(bs_env, t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                                     t_2c_Sigma_pos_tau, t_3c_x_W, &
                                     mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_2c_W, &
                                                            t_2c_Sigma_neg_tau, &
                                                            t_2c_Sigma_pos_tau, t_3c_x_W
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_for_Sigma_c'

      INTEGER                                            :: handle, i_t, ispin

      CALL timeset(routineN, handle)

      CALL dbt_create(bs_env%t_G, t_2c_Gocc)
      CALL dbt_create(bs_env%t_G, t_2c_Gvir)
      CALL dbt_create(bs_env%t_W, t_2c_W)
      CALL dbt_create(bs_env%t_G, t_2c_Sigma_neg_tau)
      CALL dbt_create(bs_env%t_G, t_2c_Sigma_pos_tau)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_W)

      NULLIFY (mat_Sigma_neg_tau, mat_Sigma_pos_tau)
      ALLOCATE (mat_Sigma_neg_tau(bs_env%num_time_freq_points, bs_env%n_spin))
      ALLOCATE (mat_Sigma_pos_tau(bs_env%num_time_freq_points, bs_env%n_spin))

      DO i_t = 1, bs_env%num_time_freq_points
         DO ispin = 1, bs_env%n_spin
            ALLOCATE (mat_Sigma_neg_tau(i_t, ispin)%matrix)
            ALLOCATE (mat_Sigma_pos_tau(i_t, ispin)%matrix)
            CALL dbcsr_create(mat_Sigma_neg_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
            CALL dbcsr_create(mat_Sigma_pos_tau(i_t, ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_mat_for_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param i_atoms ...
!> \param j_atoms ...
!> \param t_3c_x_W ...
!> \param t_2c_W ...
! **************************************************************************************************
   SUBROUTINE compute_3c_and_contract_W(qs_env, bs_env, i_atoms, j_atoms, t_3c_x_W, t_2c_W)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, DIMENSION(2)                              :: i_atoms, j_atoms
      TYPE(dbt_type)                                     :: t_3c_x_W, t_2c_W

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_and_contract_W'

      INTEGER                                            :: handle, RI_intval_idx
      INTEGER, DIMENSION(2)                              :: bounds_j, RI_atoms
      TYPE(dbt_type)                                     :: t_3c_for_W, t_3c_x_W_tmp

      CALL timeset(routineN, handle)

      CALL dbt_create(bs_env%t_RI__AO_AO, t_3c_x_W_tmp)
      CALL dbt_create(bs_env%t_RI__AO_AO, t_3c_for_W)

      bounds_j(1:2) = [bs_env%i_RI_start_from_atom(j_atoms(1)), &
                       bs_env%i_RI_end_from_atom(j_atoms(2))]

      DO RI_intval_idx = 1, bs_env%n_intervals_inner_loop_atoms
         RI_atoms = bs_env%inner_loop_atom_intervals(1:2, RI_intval_idx)

         ! 1. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
         CALL compute_3c_integrals(qs_env, bs_env, t_3c_for_W, &
                                   atoms_AO_1=i_atoms, atoms_RI=RI_atoms)

         ! 2. tensor operation M_νσQ(iτ) = sum_P (νσ|P) W^MIC_PQ(iτ)
         CALL dbt_contract(alpha=1.0_dp, &
                           tensor_1=t_2c_W, &
                           tensor_2=t_3c_for_W, &
                           beta=1.0_dp, &
                           tensor_3=t_3c_x_W_tmp, &
                           contract_1=[2], notcontract_1=[1], map_1=[1], &
                           contract_2=[1], notcontract_2=[2, 3], map_2=[2, 3], &
                           bounds_2=bounds_j, &
                           filter_eps=bs_env%eps_filter)

      END DO ! RI_atoms

      ! 3. reorder tensor
      CALL dbt_copy(t_3c_x_W_tmp, t_3c_x_W, order=[1, 2, 3], move_data=.TRUE.)

      CALL dbt_destroy(t_3c_x_W_tmp)
      CALL dbt_destroy(t_3c_for_W)

      CALL timestop(handle)

   END SUBROUTINE compute_3c_and_contract_W

! **************************************************************************************************
!> \brief ...
!> \param t_2c_G ...
!> \param t_3c_x_W ...
!> \param t_2c_Sigma ...
!> \param i_atoms ...
!> \param j_atoms ...
!> \param qs_env ...
!> \param bs_env ...
!> \param occ ...
!> \param vir ...
!> \param clear_W ...
!> \param can_skip ...
! **************************************************************************************************
   SUBROUTINE contract_to_Sigma(t_2c_G, t_3c_x_W, t_2c_Sigma, i_atoms, j_atoms, qs_env, bs_env, &
                                occ, vir, clear_W, can_skip)
      TYPE(dbt_type)                                     :: t_2c_G, t_3c_x_W, t_2c_Sigma
      INTEGER, DIMENSION(2)                              :: i_atoms, j_atoms
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL                                            :: occ, vir, clear_W
      LOGICAL, OPTIONAL                                  :: can_skip

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_to_Sigma'

      INTEGER :: handle, inner_loop_atoms_interval_index
      INTEGER(KIND=int_8)                                :: flop
      INTEGER, DIMENSION(2)                              :: bounds_i, IL_atoms
      REAL(KIND=dp)                                      :: sign_Sigma
      TYPE(dbt_type)                                     :: t_3c_for_G, t_3c_x_G, t_3c_x_G_2

      CALL timeset(routineN, handle)

      CPASSERT(occ .EQV. (.NOT. vir))
      IF (occ) sign_Sigma = -1.0_dp
      IF (vir) sign_Sigma = 1.0_dp

      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_for_G)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_G)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_G_2)

      bounds_i(1:2) = [bs_env%i_ao_start_from_atom(i_atoms(1)), &
                       bs_env%i_ao_end_from_atom(i_atoms(2))]

      DO inner_loop_atoms_interval_index = 1, bs_env%n_intervals_inner_loop_atoms
         IL_atoms = bs_env%inner_loop_atom_intervals(1:2, inner_loop_atoms_interval_index)

         CALL compute_3c_integrals(qs_env, bs_env, t_3c_for_G, &
                                   atoms_RI=j_atoms, atoms_AO_2=IL_atoms)

         CALL dbt_contract(alpha=1.0_dp, &
                           tensor_1=t_2c_G, &
                           tensor_2=t_3c_for_G, &
                           beta=1.0_dp, &
                           tensor_3=t_3c_x_G, &
                           contract_1=[2], notcontract_1=[1], map_1=[3], &
                           contract_2=[3], notcontract_2=[1, 2], map_2=[1, 2], &
                           bounds_2=bounds_i, &
                           filter_eps=bs_env%eps_filter)

      END DO ! IL_atoms

      CALL dbt_copy(t_3c_x_G, t_3c_x_G_2, order=[1, 3, 2], move_data=.TRUE.)

      CALL dbt_contract(alpha=sign_Sigma, &
                        tensor_1=t_3c_x_W, &
                        tensor_2=t_3c_x_G_2, &
                        beta=1.0_dp, &
                        tensor_3=t_2c_Sigma, &
                        contract_1=[1, 2], notcontract_1=[3], map_1=[1], &
                        contract_2=[1, 2], notcontract_2=[3], map_2=[2], &
                        filter_eps=bs_env%eps_filter, move_data=clear_W, flop=flop)

      IF (PRESENT(can_skip)) THEN
         IF (flop == 0_int_8) can_skip = .TRUE.
      END IF

      CALL dbt_destroy(t_3c_for_G)
      CALL dbt_destroy(t_3c_x_G)
      CALL dbt_destroy(t_3c_x_G_2)

      CALL timestop(handle)

   END SUBROUTINE contract_to_Sigma

! **************************************************************************************************
!> \brief ...
!> \param fm_Sigma_c_Gamma_time ...
!> \param bs_env ...
!> \param mat_Sigma_pos_tau ...
!> \param mat_Sigma_neg_tau ...
! **************************************************************************************************
   SUBROUTINE fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                         mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_pos_tau, mat_Sigma_neg_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_Sigma_c_Gamma_time'

      INTEGER                                            :: handle, i_t, ispin, pos_neg

      CALL timeset(routineN, handle)

      ALLOCATE (fm_Sigma_c_Gamma_time(bs_env%num_time_freq_points, 2, bs_env%n_spin))
      DO i_t = 1, bs_env%num_time_freq_points
         DO ispin = 1, bs_env%n_spin
            DO pos_neg = 1, 2
               CALL cp_fm_create(fm_Sigma_c_Gamma_time(i_t, pos_neg, ispin), &
                                 bs_env%fm_s_Gamma%matrix_struct)
            END DO
            CALL copy_dbcsr_to_fm(mat_Sigma_pos_tau(i_t, ispin)%matrix, &
                                  fm_Sigma_c_Gamma_time(i_t, 1, ispin))
            CALL copy_dbcsr_to_fm(mat_Sigma_neg_tau(i_t, ispin)%matrix, &
                                  fm_Sigma_c_Gamma_time(i_t, 2, ispin))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE fill_fm_Sigma_c_Gamma_time

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE print_skipping(bs_env)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'print_skipping'

      INTEGER                                            :: handle, i_intval_idx, j_intval_idx, &
                                                            n_skip

      CALL timeset(routineN, handle)

      n_skip = 0

      DO i_intval_idx = 1, bs_env%n_intervals_i
         DO j_intval_idx = 1, bs_env%n_intervals_j
            IF (bs_env%skip_Sigma_occ(i_intval_idx, j_intval_idx) .AND. &
                bs_env%skip_Sigma_vir(i_intval_idx, j_intval_idx)) THEN
               n_skip = n_skip + 1
            END IF
         END DO
      END DO

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T74,F7.1,A)') &
            'Sparsity of Σ^c(iτ,k=0): Percentage of skipped atom pairs:', &
            REAL(100*n_skip, KIND=dp)/REAL(i_intval_idx*j_intval_idx, KIND=dp), ' %'
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_skipping

! **************************************************************************************************
!> \brief ...
!> \param t_2c_Gocc ...
!> \param t_2c_Gvir ...
!> \param t_2c_W ...
!> \param t_2c_Sigma_neg_tau ...
!> \param t_2c_Sigma_pos_tau ...
!> \param t_3c_x_W ...
!> \param fm_W_MIC_time ...
!> \param mat_Sigma_neg_tau ...
!> \param mat_Sigma_pos_tau ...
! **************************************************************************************************
   SUBROUTINE destroy_mat_Sigma_c(t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                                  t_2c_Sigma_pos_tau, t_3c_x_W, fm_W_MIC_time, &
                                  mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_2c_W, &
                                                            t_2c_Sigma_neg_tau, &
                                                            t_2c_Sigma_pos_tau, t_3c_x_W
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_mat_Sigma_c'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_destroy(t_2c_Gocc)
      CALL dbt_destroy(t_2c_Gvir)
      CALL dbt_destroy(t_2c_W)
      CALL dbt_destroy(t_2c_Sigma_neg_tau)
      CALL dbt_destroy(t_2c_Sigma_pos_tau)
      CALL dbt_destroy(t_3c_x_W)
      CALL cp_fm_release(fm_W_MIC_time)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL timestop(handle)

   END SUBROUTINE destroy_mat_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE delete_unnecessary_files(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'delete_unnecessary_files'

      CHARACTER(LEN=default_string_length)               :: f_chi, f_W_t, prefix
      INTEGER                                            :: handle, i_t

      CALL timeset(routineN, handle)

      prefix = bs_env%prefix

      DO i_t = 1, bs_env%num_time_freq_points

         IF (i_t < 10) THEN
            WRITE (f_chi, '(3A,I1,A)') TRIM(prefix), bs_env%chi_name, "_00", i_t, ".matrix"
            WRITE (f_W_t, '(3A,I1,A)') TRIM(prefix), bs_env%W_time_name, "_00", i_t, ".matrix"
         ELSE IF (i_t < 100) THEN
            WRITE (f_chi, '(3A,I2,A)') TRIM(prefix), bs_env%chi_name, "_0", i_t, ".matrix"
            WRITE (f_W_t, '(3A,I2,A)') TRIM(prefix), bs_env%W_time_name, "_0", i_t, ".matrix"
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         CALL safe_delete(f_chi, bs_env)
         CALL safe_delete(f_W_t, bs_env)

      END DO

      CALL timestop(handle)

   END SUBROUTINE delete_unnecessary_files

! **************************************************************************************************
!> \brief ...
!> \param filename ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE safe_delete(filename, bs_env)
      CHARACTER(LEN=*)                                   :: filename
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'safe_delete'

      INTEGER                                            :: handle
      LOGICAL                                            :: file_exists

      CALL timeset(routineN, handle)

      IF (bs_env%para_env%mepos == 0) THEN

         INQUIRE (file=TRIM(filename), exist=file_exists)
         IF (file_exists) CALL mp_file_delete(TRIM(filename))

      END IF

      CALL timestop(handle)

   END SUBROUTINE safe_delete

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER                                            :: handle, ikp, ispin, j_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n, V_xc_ikp_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_freq, Sigma_c_ikp_n_time
      TYPE(cp_cfm_type)                                  :: cfm_ks_ikp, cfm_mos_ikp, cfm_s_ikp, &
                                                            cfm_Sigma_x_ikp, cfm_work_ikp

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mos_ikp, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_work_ikp, bs_env%fm_s_Gamma%matrix_struct)
      ! JW TODO: fully distribute these arrays at given time; also eigenvalues in bs_env
      ALLOCATE (V_xc_ikp_n(bs_env%n_ao), Sigma_x_ikp_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_ikp_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_ikp_n_freq(bs_env%n_ao, bs_env%num_time_freq_points, 2))

      DO ispin = 1, bs_env%n_spin

         DO ikp = 1, bs_env%nkp_bs_and_DOS

            ! 1. get H^KS_µν(k_i) from H^KS_µν(k=0)
            CALL cfm_ikp_from_fm_Gamma(cfm_ks_ikp, bs_env%fm_ks_Gamma(ispin), &
                                       ikp, qs_env, bs_env%kpoints_DOS, "ORB")

            ! 2. get S_µν(k_i) from S_µν(k=0)
            CALL cfm_ikp_from_fm_Gamma(cfm_s_ikp, bs_env%fm_s_Gamma, &
                                       ikp, qs_env, bs_env%kpoints_DOS, "ORB")

            ! 3. Diagonalize (Roothaan-Hall): H_KS(k_i)*C(k_i) = S(k_i)*C(k_i)*ϵ(k_i)
            CALL cp_cfm_geeig(cfm_ks_ikp, cfm_s_ikp, cfm_mos_ikp, &
                              bs_env%eigenval_scf(:, ikp, ispin), cfm_work_ikp)

            ! 4. V^xc_µν(k=0) -> V^xc_µν(k_i) -> V^xc_nn(k_i)
            CALL to_ikp_and_mo(V_xc_ikp_n, bs_env%fm_V_xc_Gamma(ispin), &
                               ikp, qs_env, bs_env, cfm_mos_ikp)

            ! 5. Σ^x_µν(k=0) -> Σ^x_µν(k_i) -> Σ^x_nn(k_i)
            CALL to_ikp_and_mo(Sigma_x_ikp_n, fm_Sigma_x_Gamma(ispin), &
                               ikp, qs_env, bs_env, cfm_mos_ikp)

            ! 6. Σ^c_µν(k=0,+/-i|τ_j|) -> Σ^c_µν(k_i,+/-i|τ_j|) -> Σ^c_nn(k_i,+/-i|τ_j|)
            DO j_t = 1, bs_env%num_time_freq_points
               CALL to_ikp_and_mo(Sigma_c_ikp_n_time(:, j_t, 1), &
                                  fm_Sigma_c_Gamma_time(j_t, 1, ispin), &
                                  ikp, qs_env, bs_env, cfm_mos_ikp)
               CALL to_ikp_and_mo(Sigma_c_ikp_n_time(:, j_t, 2), &
                                  fm_Sigma_c_Gamma_time(j_t, 2, ispin), &
                                  ikp, qs_env, bs_env, cfm_mos_ikp)
            END DO

            ! 7. Σ^c_nn(k_i,iτ) -> Σ^c_nn(k_i,iω)
            CALL time_to_freq(bs_env, Sigma_c_ikp_n_time, Sigma_c_ikp_n_freq, ispin)

            ! 8. Analytic continuation Σ^c_nn(k_i,iω) -> Σ^c_nn(k_i,ϵ) and
            !    ϵ_nk_i^GW = ϵ_nk_i^DFT + Σ^c_nn(k_i,ϵ) + Σ^x_nn(k_i) - v^xc_nn(k_i)
            CALL analyt_conti_and_print(bs_env, Sigma_c_ikp_n_freq, Sigma_x_ikp_n, V_xc_ikp_n, &
                                        bs_env%eigenval_scf(:, ikp, ispin), ikp, ispin)

         END DO ! ikp_DOS

      END DO ! ispin

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      CALL cp_fm_release(fm_Sigma_x_Gamma)
      CALL cp_fm_release(fm_Sigma_c_Gamma_time)
      CALL cp_cfm_release(cfm_ks_ikp)
      CALL cp_cfm_release(cfm_s_ikp)
      CALL cp_cfm_release(cfm_mos_ikp)
      CALL cp_cfm_release(cfm_work_ikp)
      CALL cp_cfm_release(cfm_Sigma_x_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief ...
!> \param array_ikp_n ...
!> \param fm_Gamma ...
!> \param ikp ...
!> \param qs_env ...
!> \param bs_env ...
!> \param cfm_mos_ikp ...
! **************************************************************************************************
   SUBROUTINE to_ikp_and_mo(array_ikp_n, fm_Gamma, ikp, qs_env, bs_env, cfm_mos_ikp)

      REAL(KIND=dp), DIMENSION(:)                        :: array_ikp_n
      TYPE(cp_fm_type)                                   :: fm_Gamma
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_mos_ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'to_ikp_and_mo'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fm_ikp_mo_re

      CALL timeset(routineN, handle)

      CALL cp_fm_create(fm_ikp_mo_re, fm_Gamma%matrix_struct)

      CALL fm_Gamma_ao_to_cfm_ikp_mo(fm_Gamma, fm_ikp_mo_re, ikp, qs_env, bs_env, cfm_mos_ikp)

      CALL cp_fm_get_diag(fm_ikp_mo_re, array_ikp_n)

      CALL cp_fm_release(fm_ikp_mo_re)

      CALL timestop(handle)

   END SUBROUTINE to_ikp_and_mo

! **************************************************************************************************
!> \brief ...
!> \param fm_Gamma ...
!> \param fm_ikp_mo_re ...
!> \param ikp ...
!> \param qs_env ...
!> \param bs_env ...
!> \param cfm_mos_ikp ...
! **************************************************************************************************
   SUBROUTINE fm_Gamma_ao_to_cfm_ikp_mo(fm_Gamma, fm_ikp_mo_re, ikp, qs_env, bs_env, cfm_mos_ikp)
      TYPE(cp_fm_type)                                   :: fm_Gamma, fm_ikp_mo_re
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_mos_ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_Gamma_ao_to_cfm_ikp_mo'

      INTEGER                                            :: handle, nmo
      TYPE(cp_cfm_type)                                  :: cfm_ikp_ao, cfm_ikp_mo, cfm_tmp

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_ikp_ao, fm_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_ikp_mo, fm_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_tmp, fm_Gamma%matrix_struct)

      ! get cfm_µν(k_i) from fm_µν(k=0)
      CALL cfm_ikp_from_fm_Gamma(cfm_ikp_ao, fm_Gamma, ikp, qs_env, bs_env%kpoints_DOS, "ORB")

      nmo = bs_env%n_ao
      CALL parallel_gemm('N', 'N', nmo, nmo, nmo, z_one, cfm_ikp_ao, cfm_mos_ikp, z_zero, cfm_tmp)
      CALL parallel_gemm('C', 'N', nmo, nmo, nmo, z_one, cfm_mos_ikp, cfm_tmp, z_zero, cfm_ikp_mo)

      CALL cp_cfm_to_fm(cfm_ikp_mo, fm_ikp_mo_re)

      CALL cp_cfm_release(cfm_ikp_mo)
      CALL cp_cfm_release(cfm_ikp_ao)
      CALL cp_cfm_release(cfm_tmp)

      CALL timestop(handle)

   END SUBROUTINE fm_Gamma_ao_to_cfm_ikp_mo

END MODULE gw_large_cell_gamma
