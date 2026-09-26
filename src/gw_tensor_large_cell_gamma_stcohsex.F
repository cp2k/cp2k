!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief static COHSEX band structure calculation
!> \author Ritaj Tyagi
!> \date Jan 2025
! **************************************************************************************************
MODULE gw_tensor_large_cell_Gamma_stcohsex
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_transpose
   USE cp_fm_diag,                      ONLY: cp_fm_syevx
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbt_api,                         ONLY: dbt_contract,&
                                              dbt_copy,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE gw_tensor_large_cell_gamma,      ONLY: &
        apply_extrapol_factor_single, compute_3c_and_contract_W, compute_3c_integrals, &
        compute_MinvVsqrt_Vsqrt, compute_V_k_by_lattice_sum, compute_fm_W_MIC_freq_j, &
        delete_unnecessary_files, fm_read, fm_write, get_Sigma_x, get_mat_chi_Gamma_tau, &
        print_skipping
   USE gw_utils,                        ONLY: de_init_bs_env
   USE gw_utils_communication,          ONLY: fm_to_local_tensor,&
                                              local_dbt_to_global_mat
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: cfm_ikp_from_fm_Gamma,&
                                              get_all_VBM_CBM_bandgaps
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_io,                        ONLY: write_mo_set_to_output_unit
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_tensor_large_cell_Gamma_stcohsex'

   PUBLIC :: gw_calc_tensor_large_cell_Gamma_stcohsex

CONTAINS

! **************************************************************************************************
!> \brief Perform st-COHSEX band structure calculation
!> \param qs_env ...
!> \param bs_env ...
!> \par History
!>    * 12.2024 created [Ritaj Tyagi]
! **************************************************************************************************
   SUBROUTINE gw_calc_tensor_large_cell_Gamma_stcohsex(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_tensor_large_cell_Gamma_stcohsex'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma

      CALL timeset(routineN, handle)

      ! G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      ! χ_PQ(iτ,k=0) = sum_λν [sum_µ (µν|P) G^occ_µλ(i|τ|)] [sum_σ (σλ|Q) G^vir_σν(i|τ|)]
      CALL get_mat_chi_Gamma_tau(bs_env, qs_env, bs_env%mat_chi_Gamma_tau)

      ! χ_PQ(iτ,k=0) -> χ_PQ(0,k) -> ε_PQ(0,k) -> W_PQ(0,k) -> W^MIC_PQ(0) -> M^-1*W^MIC(0)*M^-1
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_MIC_freq_zero)

      ! D_µν = sum_n^occ C_µn(k=0) C_νn(k=0), V^trunc_PQ = sum_cell_R <phi_P,0|V^trunc|phi_Q,R>
      ! Σ^x_λσ(k=0) = sum_νQ [sum_P (νσ|P) V^trunc_PQ] [sum_µ (λµ|Q) D_µν)]
      CALL get_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma)

      ! D'_µν(k=0) = sum_n C_µn(k=0) C_νn(k=0)
      ! Σ^c_λσ(k=0) = sum_νQ [sum_P (νσ|P) W^MIC_PQ(0)] [sum_µ (λµ|Q) D'_µν]
      CALL get_Sigma_c(bs_env, qs_env, fm_W_MIC_freq_zero, fm_Sigma_c_Gamma)

      ! Σ^c_λσ(k=0) -> Σ^c_nn(ϵ,k); ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ,k) + Σ^x_nn(k) - v^xc_nn(k)
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_tensor_large_cell_Gamma_stcohsex

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_freq_zero ...
! **************************************************************************************************
   SUBROUTINE get_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_freq_zero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_W_MIC'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      IF (bs_env%W_freq_zero_exists) THEN
         CALL read_W_MIC(bs_env, mat_chi_Gamma_tau, fm_W_MIC_freq_zero)
      ELSE
         CALL compute_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_freq_zero)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_freq_zero ...
! **************************************************************************************************
   SUBROUTINE read_W_MIC(bs_env, mat_chi_Gamma_tau, fm_W_MIC_freq_zero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_W_MIC'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ALLOCATE (fm_W_MIC_freq_zero)
      CALL cp_fm_create(fm_W_MIC_freq_zero, bs_env%fm_RI_RI%matrix_struct)

      t1 = m_walltime()

      CALL fm_read(fm_W_MIC_freq_zero, bs_env, bs_env%W_freq_zero_name, 0)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,A,F7.1,A)') &
            'Read W(0,k) from file  ', 'Execution time', m_walltime() - t1, ' s'
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE read_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_freq_zero ...
! **************************************************************************************************
   SUBROUTINE compute_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_freq_zero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W_MIC'

      INTEGER                                            :: handle, ikp, ikp_batch, ikp_in_batch
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_kp

      CALL timeset(routineN, handle)

      ALLOCATE (fm_W_MIC_freq_zero)
      CALL cp_fm_create(fm_W_MIC_freq_zero, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(fm_W_MIC_freq_zero, 0.0_dp)

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

            CALL compute_fm_W_MIC_freq_j(bs_env, qs_env, bs_env%fm_W_MIC_freq, 1, ikp, &
                                         mat_chi_Gamma_tau, cfm_M_inv_V_sqrt_ikp, &
                                         cfm_V_sqrt_ikp, at_zero_freq=.TRUE.)

            CALL cp_fm_scale_and_add(1.0_dp, fm_W_MIC_freq_zero, &
                                     1.0_dp, bs_env%fm_W_MIC_freq)

            CALL cp_fm_release(fm_V_kp(ikp, 1))
            CALL cp_fm_release(fm_V_kp(ikp, 2))

         END DO ! ikp_in_batch

         DEALLOCATE (fm_V_kp)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I12,A,I3,A,F7.1,A)') &
               'Computed W(0,k) for k-point batch', &
               ikp_batch, ' /', bs_env%num_chi_eps_W_batches, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! ikp_batch

      IF (bs_env%approx_kp_extrapol) THEN
         CALL apply_extrapol_factor_single(bs_env, fm_W_MIC_freq_zero)
      END IF

      ! M^-1(k=0)*W^MIC(0)*M^-1(k=0)
      CALL multiply_fm_W_MIC_zero_with_Minv_Gamma(bs_env, fm_W_MIC_freq_zero)

      CALL fm_write(fm_W_MIC_freq_zero, 0, "W_freq_zero", qs_env)

      CALL cp_cfm_release(cfm_M_inv_V_sqrt_ikp)
      CALL cp_cfm_release(cfm_V_sqrt_ikp)
      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_freq_zero ...
! **************************************************************************************************
   SUBROUTINE multiply_fm_W_MIC_zero_with_Minv_Gamma(bs_env, fm_W_MIC_freq_zero)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero

      CHARACTER(LEN=*), PARAMETER :: routineN = 'multiply_fm_W_MIC_zero_with_Minv_Gamma'

      INTEGER                                            :: handle, n_RI
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_work, bs_env%fm_RI_RI%matrix_struct)

      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, bs_env%fm_Minv_Gamma, &
                         fm_W_MIC_freq_zero, 0.0_dp, fm_work)

      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, fm_work, &
                         bs_env%fm_Minv_Gamma, 0.0_dp, fm_W_MIC_freq_zero)

      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE multiply_fm_W_MIC_zero_with_Minv_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_zero ...
!> \param fm_Sigma_c_Gamma ...
! **************************************************************************************************
   SUBROUTINE get_Sigma_c(bs_env, qs_env, fm_W_MIC_freq_zero, fm_Sigma_c_Gamma)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_Sigma_c'

      INTEGER                                            :: handle, i_intval_idx, ispin, j_intval_idx
      INTEGER, DIMENSION(2)                              :: i_atoms, j_atoms
      REAL(KIND=dp)                                      :: t1
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_Sigma
      TYPE(dbt_type)                                     :: t_2c_G, t_2c_Sigma, t_2c_W, t_3c_x_W

      CALL timeset(routineN, handle)

      IF (bs_env%SigmaC_freq_zero_exists) THEN

         CALL read_SigmaC_zero(bs_env, fm_W_MIC_freq_zero, fm_Sigma_c_Gamma)

      ELSE

         CALL dbt_create(bs_env%t_W, t_2c_W)

         NULLIFY (mat_Sigma)
         ALLOCATE (mat_Sigma(bs_env%n_spin))

         ! fm W to local tensor
         CALL fm_to_local_tensor(fm_W_MIC_freq_zero, bs_env%mat_RI_RI%matrix, &
                                 bs_env%mat_RI_RI_tensor%matrix, t_2c_W, bs_env, &
                                 bs_env%atoms_j_t_group)

         CALL cp_fm_release(fm_W_MIC_freq_zero)

         DO ispin = 1, bs_env%n_spin

            t1 = m_walltime()

            CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_x_W)

            CALL dbt_create(bs_env%t_G, t_2c_G)
            CALL dbt_create(bs_env%t_G, t_2c_Sigma)

            CALL D_all(bs_env, bs_env%fm_work_mo(2), ispin)

            ! fm G^sp to local tensor
            CALL fm_to_local_tensor(bs_env%fm_work_mo(2), bs_env%mat_ao_ao%matrix, &
                                    bs_env%mat_ao_ao_tensor%matrix, t_2c_G, bs_env, &
                                    bs_env%atoms_i_t_group)

            ! every group has its own range of i_atoms and j_atoms; only deal with a
            ! limited number of i_atom-j_atom pairs simultaneously in a group to save memory
            DO i_intval_idx = 1, bs_env%n_intervals_i
               DO j_intval_idx = 1, bs_env%n_intervals_j

                  i_atoms = bs_env%i_atom_intervals(1:2, i_intval_idx)
                  j_atoms = bs_env%j_atom_intervals(1:2, j_intval_idx)

                  IF (bs_env%skip_Sigma_occ(i_intval_idx, j_intval_idx) .AND. &
                      bs_env%skip_Sigma_vir(i_intval_idx, j_intval_idx)) CYCLE

                  ! 1. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
                  ! 2. tensor operation M_νσQ(0) = sum_P (νσ|P) W^MIC_PQ(0)
                  CALL compute_3c_and_contract_W(qs_env, bs_env, i_atoms, j_atoms, t_3c_x_W, t_2c_W)

                  ! 3. Σ_λσ(k=0) = sum_νQ M_νσQ(0) sum_µ (λµ|Q) G^sp_µν
                  CALL contract_to_Sigma(t_2c_G, t_3c_x_W, t_2c_Sigma, i_atoms, j_atoms, &
                                         qs_env, bs_env)

               END DO ! j_atoms
            END DO ! i_atoms

            CALL dbt_destroy(t_3c_x_W)

            CALL dbt_destroy(t_2c_G)

            ! 4. communicate data tensor t_2c_Sigma (which is local in the subgroup)
            !    to the global dbcsr matrix mat_Sigma(which stores Σ)

            ALLOCATE (mat_Sigma(ispin)%matrix)
            CALL dbcsr_create(mat_Sigma(ispin)%matrix, template=bs_env%mat_ao_ao%matrix)
            CALL local_dbt_to_global_mat(t_2c_Sigma, bs_env%mat_ao_ao_tensor%matrix, &
                                         mat_Sigma(ispin)%matrix, bs_env%para_env)
            CALL dbt_destroy(t_2c_Sigma)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,A,F7.1,A)') &
                  'Computed Σ^c(k=0) ', &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO

         CALL dbt_destroy(t_2c_W)

         CALL fill_fm_Sigma_c_Gamma(fm_Sigma_c_Gamma, bs_env, mat_Sigma)

         CALL fm_write(fm_Sigma_c_Gamma(1, 1), 0, "SigmaC_freq_zero", qs_env)

         CALL print_skipping(bs_env)

         CALL delete_unnecessary_files(bs_env)

      END IF

      CALL timestop(handle)

   END SUBROUTINE get_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_freq_zero ...
!> \param fm_Sigma_c_Gamma ...
! **************************************************************************************************

   SUBROUTINE read_SigmaC_zero(bs_env, fm_W_MIC_freq_zero, fm_Sigma_c_Gamma)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE                      :: fm_W_MIC_freq_zero
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_SigmaC_zero'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      CALL cp_fm_release(fm_W_MIC_freq_zero)

      ALLOCATE (fm_Sigma_c_Gamma(1, 1))
      CALL cp_fm_create(fm_Sigma_c_Gamma(1, 1), bs_env%fm_s_Gamma%matrix_struct)

      t1 = m_walltime()

      CALL fm_read(fm_Sigma_c_Gamma(1, 1), bs_env, bs_env%SigmaC_freq_zero_name, 0)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,A,F7.1,A)') &
            'Read Σ^c(k=0) from file  ', 'Execution time', m_walltime() - t1, ' s'
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE read_SigmaC_zero

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_G_Gamma ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE D_all(bs_env, fm_G_Gamma, ispin)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type)                                   :: fm_G_Gamma
      INTEGER                                            :: ispin

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'D_all'

      INTEGER                                            :: handle, homo, j_col, n_mo, n_retained
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: col_scale

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao
      homo = bs_env%n_occ(ispin)

      n_retained = MIN(bs_env%n_mo_retained, n_mo)
      IF (n_retained < n_mo .AND. bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,I0,A,I0,A)') &
            'st-COHSEX: Coulomb-hole completeness sum restricted to ', n_retained, ' of ', &
            n_mo, ' basis modes (rest removed by EPS_EIGVAL).'
      END IF

      CALL cp_fm_to_fm(bs_env%fm_mo_coeff_Gamma(ispin), bs_env%fm_work_mo(1))

      ! scale each MO column n by -0.5 (occupied) or +0.5 (virtual); dropped modes by 0
      ALLOCATE (col_scale(n_mo))
      DO j_col = 1, n_mo
         IF (j_col > n_retained) THEN
            col_scale(j_col) = 0.0_dp
         ELSE IF (j_col <= homo) THEN
            col_scale(j_col) = -0.5_dp
         ELSE
            col_scale(j_col) = 0.5_dp
         END IF
      END DO
      CALL cp_fm_column_scale(bs_env%fm_work_mo(1), col_scale)
      DEALLOCATE (col_scale)

      ! G_µν = sum_n C_µn (factor_n C_νn) ; beta=0 overwrites fm_G_Gamma
      CALL parallel_gemm(transa="N", transb="T", m=n_mo, n=n_mo, k=n_mo, alpha=1.0_dp, &
                         matrix_a=bs_env%fm_mo_coeff_Gamma(ispin), matrix_b=bs_env%fm_work_mo(1), &
                         beta=0.0_dp, matrix_c=fm_G_Gamma)

      CALL timestop(handle)

   END SUBROUTINE D_all

! *************************************************************************************************
!> \brief ...
!> \param t_2c_G ...
!> \param t_3c_x_W ...
!> \param t_2c_Sigma ...
!> \param i_atoms ...
!> \param j_atoms ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************

   SUBROUTINE contract_to_Sigma(t_2c_G, t_3c_x_W, t_2c_Sigma, i_atoms, j_atoms, qs_env, bs_env)
      TYPE(dbt_type)                                     :: t_2c_G, t_3c_x_W, t_2c_Sigma
      INTEGER, DIMENSION(2)                              :: i_atoms, j_atoms
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_to_Sigma'

      INTEGER :: handle, inner_loop_atoms_interval_index
      INTEGER, DIMENSION(2)                              :: bounds_i, IL_atoms
      TYPE(dbt_type)                                     :: t_3c_for_G, t_3c_x_G, t_3c_x_G_2

      CALL timeset(routineN, handle)

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
      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=t_3c_x_W, &
                        tensor_2=t_3c_x_G_2, &
                        beta=1.0_dp, &
                        tensor_3=t_2c_Sigma, &
                        contract_1=[1, 2], notcontract_1=[3], map_1=[1], &
                        contract_2=[1, 2], notcontract_2=[3], map_2=[2], &
                        filter_eps=bs_env%eps_filter)

      CALL dbt_destroy(t_3c_for_G)

      CALL dbt_destroy(t_3c_x_G)
      CALL dbt_destroy(t_3c_x_G_2)

      CALL timestop(handle)

   END SUBROUTINE contract_to_Sigma

! **************************************************************************************************
!> \brief ...
!> \param fm_Sigma_c_Gamma ...
!> \param bs_env ...
!> \param mat_Sigma ...
! **************************************************************************************************
   SUBROUTINE fill_fm_Sigma_c_Gamma(fm_Sigma_c_Gamma, bs_env, mat_Sigma)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_Sigma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_Sigma_c_Gamma'

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      ALLOCATE (fm_Sigma_c_Gamma(1, bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_c_Gamma(1, ispin), &
                           bs_env%fm_s_Gamma%matrix_struct)
         CALL copy_dbcsr_to_fm(mat_Sigma(ispin)%matrix, &
                               fm_Sigma_c_Gamma(1, ispin))
      END DO

      CALL dbcsr_deallocate_matrix_set(mat_Sigma)

      CALL timestop(handle)

   END SUBROUTINE fill_fm_Sigma_c_Gamma

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param fm_Sigma_c_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER                                            :: handle, ikp, ispin
      TYPE(cp_cfm_type)                                  :: cfm_ks_ikp, cfm_mos_ikp, cfm_s_ikp, &
                                                            cfm_work_ikp
      TYPE(cp_fm_type)                                   :: fm_Sigma_c_ikp_mo, fm_Sigma_x_ikp_mo, &
                                                            fm_V_xc_ikp_mo

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mos_ikp, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_work_ikp, bs_env%fm_s_Gamma%matrix_struct)

      ! MO-basis operator matrices at k_i (n_ao x n_ao), reused for every k-point
      CALL cp_fm_create(fm_V_xc_ikp_mo, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_Sigma_x_ikp_mo, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(fm_Sigma_c_ikp_mo, bs_env%fm_s_Gamma%matrix_struct)

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

            ! 4. V^xc_µν(k=0) -> V^xc_µν(k_i) -> V^xc_mn(k_i)
            CALL fm_Gamma_ao_to_cfm_ikp_mo(bs_env%fm_V_xc_Gamma(ispin), fm_V_xc_ikp_mo, &
                                           ikp, qs_env, bs_env, cfm_mos_ikp)

            ! 5. Σ^x_µν(k=0) -> Σ^x_µν(k_i) -> Σ^x_mn(k_i)
            CALL fm_Gamma_ao_to_cfm_ikp_mo(fm_Sigma_x_Gamma(ispin), fm_Sigma_x_ikp_mo, &
                                           ikp, qs_env, bs_env, cfm_mos_ikp)

            ! 6. Σ^c_µν(ϵ(k=0),k=0) -> Σ^c_mn(k_i), symmetrized
            CALL get_Sigma_c_ikp_mo(bs_env, qs_env, fm_Sigma_c_Gamma(:, ispin), ikp, &
                                    cfm_mos_ikp, fm_Sigma_c_ikp_mo)

            ! 7. build st-COHSEX QP Hamiltonian, diagonalize, and print
            CALL diagonalize_and_print(bs_env, qs_env, fm_Sigma_x_ikp_mo, fm_V_xc_ikp_mo, &
                                       fm_Sigma_c_ikp_mo, ikp, ispin)

         END DO ! ikp_DOS

      END DO ! ispin

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      CALL cp_fm_release(fm_Sigma_x_Gamma)
      CALL cp_fm_release(fm_Sigma_c_Gamma)
      CALL cp_fm_release(fm_V_xc_ikp_mo)
      CALL cp_fm_release(fm_Sigma_x_ikp_mo)
      CALL cp_fm_release(fm_Sigma_c_ikp_mo)
      CALL cp_cfm_release(cfm_ks_ikp)
      CALL cp_cfm_release(cfm_s_ikp)
      CALL cp_cfm_release(cfm_mos_ikp)
      CALL cp_cfm_release(cfm_work_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Gamma ...
!> \param ikp ...
!> \param cfm_mos_ikp ...
!> \param fm_Sigma_c_ikp_mo ...
! **************************************************************************************************
   SUBROUTINE get_Sigma_c_ikp_mo(bs_env, qs_env, fm_Gamma, ikp, cfm_mos_ikp, fm_Sigma_c_ikp_mo)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_Gamma
      INTEGER                                            :: ikp
      TYPE(cp_cfm_type)                                  :: cfm_mos_ikp
      TYPE(cp_fm_type)                                   :: fm_Sigma_c_ikp_mo

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_Sigma_c_ikp_mo'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fm_Sigma_c_t

      CALL timeset(routineN, handle)

      ! Σ^c_µν(k=0) -> Σ^c_mn(k_i) (full MO-basis matrix)
      CALL fm_Gamma_ao_to_cfm_ikp_mo(fm_Gamma(1), fm_Sigma_c_ikp_mo, &
                                     ikp, qs_env, bs_env, cfm_mos_ikp)

      ! symmetrize: Σ^c <- 1/2 (Σ^c + Σ^c^T)
      CALL cp_fm_create(fm_Sigma_c_t, fm_Sigma_c_ikp_mo%matrix_struct)
      CALL cp_fm_transpose(fm_Sigma_c_ikp_mo, fm_Sigma_c_t)
      CALL cp_fm_scale_and_add(0.5_dp, fm_Sigma_c_ikp_mo, 0.5_dp, fm_Sigma_c_t)
      CALL cp_fm_release(fm_Sigma_c_t)

      CALL timestop(handle)

   END SUBROUTINE get_Sigma_c_ikp_mo

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
      ! get cfm_µν(k_i) from fm_µν(k=0)
      CALL cfm_ikp_from_fm_Gamma(cfm_ikp_ao, fm_Gamma, ikp, qs_env, bs_env%kpoints_DOS, "ORB")
      nmo = bs_env%n_ao

      CALL cp_cfm_create(cfm_ikp_mo, fm_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_tmp, fm_Gamma%matrix_struct)
      CALL parallel_gemm('N', 'N', nmo, nmo, nmo, z_one, cfm_ikp_ao, cfm_mos_ikp, z_zero, cfm_tmp)
      CALL parallel_gemm('C', 'N', nmo, nmo, nmo, z_one, cfm_mos_ikp, cfm_tmp, z_zero, cfm_ikp_mo)

      CALL cp_cfm_to_fm(cfm_ikp_mo, fm_ikp_mo_re)

      CALL cp_cfm_release(cfm_ikp_mo)
      CALL cp_cfm_release(cfm_ikp_ao)
      CALL cp_cfm_release(cfm_tmp)

      CALL timestop(handle)

   END SUBROUTINE fm_Gamma_ao_to_cfm_ikp_mo

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_ikp_mo ...
!> \param fm_V_xc_ikp_mo ...
!> \param fm_Sigma_c_ikp_mo ...
!> \param ikp ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE diagonalize_and_print(bs_env, qs_env, fm_Sigma_x_ikp_mo, fm_V_xc_ikp_mo, &
                                    fm_Sigma_c_ikp_mo, ikp, ispin)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_Sigma_x_ikp_mo, fm_V_xc_ikp_mo, &
                                                            fm_Sigma_c_ikp_mo
      INTEGER                                            :: ikp, ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_and_print'

      CHARACTER(len=3)                                   :: occ_vir
      CHARACTER(len=default_string_length)               :: fname
      INTEGER                                            :: handle, i_col, i_row, ikp_for_print, &
                                                            iunit, j_col, j_row, m, n_mo, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: is_bandstruc_kpoint, print_DOS_kpoints, &
                                                            print_ikp
      REAL(KIND=dp)                                      :: val
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag_Sigma_c, diag_Sigma_x, diag_SP, &
                                                            diag_V_xc, SP_eigval
      TYPE(cp_fm_type)                                   :: fm_SP_mat, SP_eigvec, SP_eigvec_ao

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao

      CALL cp_fm_create(fm_SP_mat, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(SP_eigvec, bs_env%fm_s_Gamma%matrix_struct)
      CALL cp_fm_create(SP_eigvec_ao, bs_env%fm_s_Gamma%matrix_struct)

      ALLOCATE (SP_eigval(n_mo))
      SP_eigval(:) = 0.0_dp

      ! Build the st-COHSEX QP Hamiltonian directly in the distributed fm by looping
      ! over local elements only (global indices from row_indices/col_indices):
      !   occ-occ block : Σ^x - v^xc + stfacip*Σ^c
      !   vir-vir block : Σ^x - v^xc + stfacea*Σ^c
      !   occ-vir / vir-occ : 0
      !   diagonal additionally gets + ϵ^DFT_nk
      CALL cp_fm_get_info(fm_SP_mat, nrow_local=nrow_local, ncol_local=ncol_local, &
                          row_indices=row_indices, col_indices=col_indices)

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
!$OMP             SHARED(ncol_local, nrow_local, col_indices, row_indices, bs_env, ispin, ikp, &
!$OMP                    fm_Sigma_x_ikp_mo, fm_V_xc_ikp_mo, fm_Sigma_c_ikp_mo, fm_SP_mat) &
!$OMP             PRIVATE(i_col, i_row, j_col, j_row, val)
      DO i_col = 1, ncol_local
         DO i_row = 1, nrow_local
            j_col = col_indices(i_col)
            j_row = row_indices(i_row)
            IF (j_row <= bs_env%n_occ(ispin) .AND. j_col <= bs_env%n_occ(ispin)) THEN
               val = fm_Sigma_x_ikp_mo%local_data(i_row, i_col) &
                     - fm_V_xc_ikp_mo%local_data(i_row, i_col) &
                     + bs_env%stfacip*fm_Sigma_c_ikp_mo%local_data(i_row, i_col)
            ELSE IF (j_row > bs_env%n_occ(ispin) .AND. j_col > bs_env%n_occ(ispin)) THEN
               val = fm_Sigma_x_ikp_mo%local_data(i_row, i_col) &
                     - fm_V_xc_ikp_mo%local_data(i_row, i_col) &
                     + bs_env%stfacea*fm_Sigma_c_ikp_mo%local_data(i_row, i_col)
            ELSE
               val = 0.0_dp
            END IF
            IF (j_row == j_col) val = val + bs_env%eigenval_scf(j_row, ikp, ispin)
            fm_SP_mat%local_data(i_row, i_col) = val
         END DO
      END DO
!$OMP END PARALLEL DO

      ! diagonal elements (for the per-MO QP printout) before diagonalization
      ALLOCATE (diag_SP(n_mo), diag_Sigma_c(n_mo), diag_Sigma_x(n_mo), diag_V_xc(n_mo))
      CALL cp_fm_get_diag(fm_SP_mat, diag_SP)
      CALL cp_fm_get_diag(fm_Sigma_c_ikp_mo, diag_Sigma_c)
      CALL cp_fm_get_diag(fm_Sigma_x_ikp_mo, diag_Sigma_x)
      CALL cp_fm_get_diag(fm_V_xc_ikp_mo, diag_V_xc)

      CALL cp_fm_syevx(matrix=fm_SP_mat, eigenvectors=SP_eigvec, eigenvalues=SP_eigval)
      IF (ikp == bs_env%nkp_bs_and_DOS) THEN
         CALL fm_write_formatted(SP_eigvec, 0, "Singlepole_Eigenvector_Gamma", qs_env)
      END IF

      ! transform QP eigenvectors from the KS-MO basis to the AO basis
      CALL parallel_gemm(transa="N", transb="N", m=n_mo, n=n_mo, k=n_mo, alpha=1.0_dp, &
                         matrix_a=bs_env%fm_mo_coeff_Gamma(ispin), matrix_b=SP_eigvec, &
                         beta=0.0_dp, matrix_c=SP_eigvec_ao)

      bs_env%eigenval_GW(:, ikp, ispin) = SP_eigval(:)
      IF (ALLOCATED(bs_env%eigenval_G0W0)) THEN
         bs_env%eigenval_G0W0(:, ikp, ispin) = SP_eigval(:)
      END IF
      DO m = 1, n_mo
         bs_env%eigenval_HF(m, ikp, ispin) = bs_env%eigenval_scf(m, ikp, ispin) + &
                                             diag_Sigma_x(m) - diag_V_xc(m)
      END DO

      ! only print eigenvalues of DOS k-points in case no bandstructure path has been given;
      ! in kpoints_DOS the last nkp_only_bs entries are the bandstructure-path k-points
      print_DOS_kpoints = (bs_env%nkp_only_bs <= 0)
      is_bandstruc_kpoint = (ikp > bs_env%nkp_only_DOS)
      print_ikp = print_DOS_kpoints .OR. is_bandstruc_kpoint

      !! Printing Quasiparticle energies to a file
      IF (bs_env%para_env%is_source() .AND. print_ikp) THEN

         IF (print_DOS_kpoints) THEN
            ikp_for_print = ikp
         ELSE
            ikp_for_print = ikp - bs_env%nkp_only_DOS
         END IF

         fname = "Quasiparticle_energies_singlepole"
         IF (ikp_for_print == 1 .AND. ispin == 1) THEN
            CALL open_file(TRIM(fname), unit_number=iunit, file_status="REPLACE", file_action="WRITE")
         ELSE
            CALL open_file(TRIM(fname), unit_number=iunit, file_status="OLD", file_action="WRITE", &
                           file_position="APPEND")
         END IF

         WRITE (iunit, "(A)") " "
         WRITE (iunit, "(A10,I7,A25,3F10.4)") "kpoint: ", ikp_for_print, "coordinate: ", &
            bs_env%kpoints_DOS%xkp(:, ikp)
         WRITE (iunit, "(A)") " "
         WRITE (iunit, "(A5,A24,2A17,A16,A18,A18)") "n", "ϵ_nk^DFT (eV)", "Σ^c_nk (eV)", &
            "Σ^x_nk (eV)", "v_nk^xc (eV)", "ϵ_nk^DG (eV)", "ϵ_nk^CN (eV)"
         WRITE (iunit, "(A)") " "

         DO m = 1, n_mo
            IF (m <= bs_env%n_occ(ispin)) occ_vir = 'occ'
            IF (m > bs_env%n_occ(ispin)) occ_vir = 'vir'
            WRITE (iunit, "(I5,3A,4F16.3,F17.3,F17.3)") m, ' (', occ_vir, ') ', &
               bs_env%eigenval_scf(m, ikp, ispin)*evolt, &
               diag_Sigma_c(m)*evolt, &
               diag_Sigma_x(m)*evolt, &
               diag_V_xc(m)*evolt, &
               diag_SP(m)*evolt, &
               SP_eigval(m)*evolt
         END DO

         WRITE (iunit, "(A)") " "

         CALL close_file(iunit)
      END IF

      ! MO energies/coefficients and MO cubes (DFT%PRINT%MO, DFT%PRINT%MO_CUBES), Gamma point only
      IF (ikp == 1) THEN
         CALL print_stcohsex_mos(bs_env, qs_env, SP_eigvec_ao, SP_eigval, ispin)
      END IF

      CALL cp_fm_release(fm_SP_mat)
      CALL cp_fm_release(SP_eigvec)
      CALL cp_fm_release(SP_eigvec_ao)
      DEALLOCATE (SP_eigval, diag_SP, diag_Sigma_c, diag_Sigma_x, diag_V_xc)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_and_print

! **************************************************************************************************
!> \brief Print st-COHSEX QP MO energies/occupations/coefficients and MO cubes at the Gamma point,
!>        reusing the standard DFT%PRINT%MO / DFT%PRINT%MO_CUBES print keys.
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_qp_mo_ao ...
!> \param qp_eigval ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE print_stcohsex_mos(bs_env, qs_env, fm_qp_mo_ao, qp_eigval, ispin)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_qp_mo_ao
      REAL(KIND=dp), DIMENSION(:)                        :: qp_eigval
      INTEGER                                            :: ispin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'print_stcohsex_mos'

      CHARACTER(LEN=5)                                   :: spin_label
      CHARACTER(LEN=default_path_length)                 :: fname, title
      INTEGER                                            :: first_mo, handle, ivector, last_mo, &
                                                            n_mo, nelectron, nhomo, nlumo, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL                                            :: print_coeff, print_cubes, &
                                                            print_energies, print_mo, print_occ
      REAL(KIND=dp)                                      :: maxocc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type)                                  :: qp_mo_set
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_c1d_gs_type)                               :: wf_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      CALL timeset(routineN, handle)

      n_mo = bs_env%n_ao

      CALL get_qs_env(qs_env, input=input)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      logger => cp_get_default_logger()

      ! ---------------------------------------------------------------------------------------
      ! MO energies / occupation numbers / coefficients
      ! ---------------------------------------------------------------------------------------
      CALL section_vals_val_get(dft_section, "PRINT%MO%ENERGIES", l_val=print_energies)
      CALL section_vals_val_get(dft_section, "PRINT%MO%COEFFICIENTS", l_val=print_coeff)
      CALL section_vals_val_get(dft_section, "PRINT%MO%OCCUPATION_NUMBERS", l_val=print_occ)
      print_mo = (cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO") /= 0)

      IF (print_mo .AND. (print_energies .OR. print_occ .OR. print_coeff)) THEN

         ! occupation numbers of the QP states: the lowest n_occ are occupied
         IF (bs_env%n_spin == 1) THEN
            maxocc = 2.0_dp
         ELSE
            maxocc = 1.0_dp
         END IF
         nelectron = NINT(maxocc*REAL(bs_env%n_occ(ispin), KIND=dp))

         CALL allocate_mo_set(qp_mo_set, nao=n_mo, nmo=n_mo, nelectron=nelectron, &
                              n_el_f=REAL(nelectron, KIND=dp), maxocc=maxocc, &
                              flexible_electron_count=0.0_dp)
         CALL init_mo_set(qp_mo_set, fm_struct=fm_qp_mo_ao%matrix_struct, name="stcohsex")
         qp_mo_set%eigenvalues(:) = qp_eigval(:)
         CALL cp_fm_to_fm(fm_qp_mo_ao, qp_mo_set%mo_coeff)

         CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)

         IF (bs_env%n_spin == 1) THEN
            CALL write_mo_set_to_output_unit(qp_mo_set, qs_kind_set, particle_set, dft_section, &
                                             4, 0, final_mos=.TRUE.)
         ELSE
            IF (ispin == 1) THEN
               spin_label = "ALPHA"
            ELSE
               spin_label = "BETA"
            END IF
            CALL write_mo_set_to_output_unit(qp_mo_set, qs_kind_set, particle_set, dft_section, &
                                             4, 0, final_mos=.TRUE., spin=TRIM(spin_label))
         END IF

         CALL deallocate_mo_set(qp_mo_set)
      END IF

      ! ---------------------------------------------------------------------------------------
      ! MO cubes
      ! ---------------------------------------------------------------------------------------
      print_cubes = (cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO_CUBES") /= 0)

      IF (print_cubes) THEN
         CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%NHOMO", i_val=nhomo)
         CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%NLUMO", i_val=nlumo)
         CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%STRIDE", i_vals=stride)

         IF (nhomo == -1) nhomo = bs_env%n_occ(ispin)
         IF (nlumo == -1) nlumo = n_mo - bs_env%n_occ(ispin)
         first_mo = MAX(1, bs_env%n_occ(ispin) - nhomo + 1)
         last_mo = MIN(n_mo, bs_env%n_occ(ispin) + nlumo)

         CALL get_qs_env(qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         cell=cell, &
                         dft_control=dft_control, &
                         particle_set=particle_set, &
                         pw_env=pw_env, &
                         subsys=subsys)
         CALL qs_subsys_get(subsys, particles=particles)

         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL auxbas_pw_pool%create_pw(wf_r)
         CALL auxbas_pw_pool%create_pw(wf_g)

         DO ivector = first_mo, last_mo
            CALL calculate_wavefunction(fm_qp_mo_ao, ivector, wf_r, wf_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env)

            WRITE (fname, '(A,I5.5,A,I1)') "stcohsex_WFN_", ivector, "_", ispin
            unit_nr = cp_print_key_unit_nr(logger, dft_section, "PRINT%MO_CUBES", &
                                           extension=".cube", middle_name=TRIM(fname), &
                                           file_position="REWIND", log_filename=.FALSE.)
            WRITE (title, *) "st-COHSEX WAVEFUNCTION ", ivector, " spin ", ispin, &
               " i.e. HOMO ", ivector - bs_env%n_occ(ispin)
            CALL cp_pw_to_cube(wf_r, unit_nr, title, particles, stride=stride)
            CALL cp_print_key_finished_output(unit_nr, logger, dft_section, "PRINT%MO_CUBES")
         END DO

         CALL auxbas_pw_pool%give_back_pw(wf_r)
         CALL auxbas_pw_pool%give_back_pw(wf_g)
      END IF

      CALL timestop(handle)

   END SUBROUTINE print_stcohsex_mos

! **************************************************************************************************
!> \brief ...
!> \param fm ...
!> \param matrix_index ...
!> \param matrix_name ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE fm_write_formatted(fm, matrix_index, matrix_name, qs_env)
      TYPE(cp_fm_type)                                   :: fm
      INTEGER                                            :: matrix_index
      CHARACTER(LEN=*)                                   :: matrix_name
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'fm_write_formatted'

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
                                        file_form="FORMATTED", middle_name=TRIM(filename), &
                                        file_position="REWIND", file_action="WRITE")

         CALL cp_fm_write_formatted(fm, unit_nr)
         IF (unit_nr > 0) THEN
            CALL close_file(unit_nr)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE fm_write_formatted

END MODULE gw_tensor_large_cell_Gamma_stcohsex

