!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 07.2023
! **************************************************************************************************
MODULE gw_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Wilhelm2021,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert
   USE cp_cfm_diag,                     ONLY: cp_cfm_geeig
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_get_block_p, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_type
   USE dbt_api,                         ONLY: dbt_clear,&
                                              dbt_contract,&
                                              dbt_copy,&
                                              dbt_copy_matrix_to_tensor,&
                                              dbt_copy_tensor_to_matrix,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE gw_communication,                ONLY: global_matrix_to_local_matrix,&
                                              local_matrix_to_global_matrix
   USE gw_types,                        ONLY: gw_type
   USE gw_utils,                        ONLY: create_and_init_gw_env,&
                                              de_init_gw_env
   USE input_constants,                 ONLY: ri_rpa_g0w0_crossing_newton,&
                                              xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat,&
                                              inversion_of_M_and_mult_with_chol_dec_of_V
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_tensors,                      ONLY: build_3c_integrals
   USE rpa_gw,                          ONLY: continuation_pade
   USE rpa_gw_im_time_util,             ONLY: compute_weight_re_im,&
                                              get_atom_index_from_basis_function_index
   USE rpa_gw_kpoints_util,             ONLY: cp_cfm_power,&
                                              cp_cfm_upper_to_full
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_methods'

   PUBLIC :: gw

CONTAINS

! **************************************************************************************************
!> \brief Perform GW band structure calculation
!> \param qs_env  Quickstep environment
!>        gw_input_section
!> \param gw_env ...
!> \param gw_input_section ...
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE gw(qs_env, gw_env, gw_input_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(section_vals_type), POINTER                   :: gw_input_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gw'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma_time
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CALL timeset(routineN, handle)

      CALL cite_reference(Wilhelm2021)

      CALL create_and_init_gw_env(qs_env, gw_env, gw_input_section)

      ! G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      ! χ_PQ(iτ,k=0) = sum_λν [sum_µ (µν|P) G^occ_µλ(i|τ|)] [sum_σ (σλ|Q) G^vir_σν(i|τ|)]
      CALL compute_chi_Gamma(gw_env, qs_env, mat_chi_Gamma_tau)

      ! χ_PQ(iτ,k=0) -> χ_PQ(iω,k) -> ε_PQ(iω,k) -> W_PQ(iω,k) -> W^MIC_PQ(iτ)
      CALL compute_W_MIC(gw_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)

      ! D_µν = sum_n^occ C_µn(k=0) C_νn(k=0), V^trunc_PQ = sum_cell_R <phi_P,0|V^trunc|phi_Q,R>
      ! Σ^x_λσ(k=0) = sum_νQ [sum_P (νσ|P) V^trunc_PQ] [sum_µ (λµ|Q) D_µν)]
      CALL compute_Sigma_x(gw_env, qs_env, fm_W_MIC_time, fm_Sigma_x_Gamma)

      ! Σ^c_λσ(iτ,k=0) = sum_νQ [sum_P (νσ|P) W^MIC_PQ(iτ)] [sum_µ (λµ|Q) G^occ_µν(i|τ|)], τ < 0
      ! Σ^c_λσ(iτ,k=0) = sum_νQ [sum_P (νσ|P) W^MIC_PQ(iτ)] [sum_µ (λµ|Q) G^vir_µν(i|τ|)], τ > 0
      CALL compute_Sigma_c(gw_env, qs_env, fm_W_MIC_time, fm_Sigma_c_Gamma_time)

      ! Σ^c_λσ(iτ,k=0) -> Σ^c_nn(ϵ,k); ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ,k) + Σ^x_nn(k) - v^xc_nn(k)
      CALL compute_QP_energies(gw_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_gw_env(gw_env)

      CALL timestop(handle)

   END SUBROUTINE gw

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param tau ...
!> \param fm_GGamma ...
!> \param occ ...
!> \param vir ...
! **************************************************************************************************
   SUBROUTINE G_occ_vir(gw_env, tau, fm_GGamma, occ, vir)
      TYPE(gw_type), POINTER                             :: gw_env
      REAL(KIND=dp)                                      :: tau
      TYPE(cp_fm_type)                                   :: fm_GGamma
      LOGICAL                                            :: occ, vir

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_occ_vir'

      INTEGER                                            :: handle, homo, i_row_local, j_col, &
                                                            j_col_local, n_mo, ncol_local, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp)                                      :: tau_E

      CALL timeset(routineN, handle)

      CPASSERT(occ .NEQV. vir)

      CALL cp_fm_get_info(matrix=gw_env%fm_work_mo(1), &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          col_indices=col_indices)

      n_mo = gw_env%n_ao
      homo = gw_env%n_occ(1)

      CALL cp_fm_to_fm(gw_env%fm_mo_coeff_Gamma(1), gw_env%fm_work_mo(1))

      DO i_row_local = 1, nrow_local
         DO j_col_local = 1, ncol_local

            j_col = col_indices(j_col_local)

            tau_E = ABS(tau*0.5_dp*(gw_env%eigenvalues_scf_Gamma(j_col) - gw_env%e_fermi(1)))

            IF (tau_E < gw_env%stabilize_exp) THEN
               gw_env%fm_work_mo(1)%local_data(i_row_local, j_col_local) = &
                  gw_env%fm_work_mo(1)%local_data(i_row_local, j_col_local)*EXP(-tau_E)
            ELSE
               gw_env%fm_work_mo(1)%local_data(i_row_local, j_col_local) = 0.0_dp
            END IF

            IF ((occ .AND. j_col > homo) .OR. (vir .AND. j_col <= homo)) THEN
               gw_env%fm_work_mo(1)%local_data(i_row_local, j_col_local) = 0.0_dp
            END IF

         END DO
      END DO

      CALL parallel_gemm(transa="N", transb="T", m=n_mo, n=n_mo, k=n_mo, alpha=1.0_dp, &
                         matrix_a=gw_env%fm_work_mo(1), matrix_b=gw_env%fm_work_mo(1), &
                         beta=0.0_dp, matrix_c=fm_GGamma)

      CALL timestop(handle)

   END SUBROUTINE G_occ_vir

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
! **************************************************************************************************
   SUBROUTINE compute_chi_Gamma(gw_env, qs_env, mat_chi_Gamma_tau)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_chi_Gamma'

      INTEGER                                            :: handle, i_atom, i_t, j_atom
      REAL(KIND=dp)                                      :: tau
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, &
                                                            t_3c_for_Gvir, t_3c_x_Gocc, &
                                                            t_3c_x_Gocc_2, t_3c_x_Gvir, &
                                                            t_3c_x_Gvir_2

      CALL timeset(routineN, handle)

      CALL create_mat_for_chi(t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, t_3c_for_Gvir, t_3c_x_Gocc, &
                              t_3c_x_Gvir, t_3c_x_Gocc_2, t_3c_x_Gvir_2, mat_chi_Gamma_tau, gw_env)

      DO i_t = 1, gw_env%num_time_freq_points

         gw_env%t1 = m_walltime()

         ! 1. compute G^occ and G^vir
         !    Background: G(iτ) = G^occ(iτ) * Θ(-τ) + G^vir(iτ) * Θ(τ)
         !    G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
         !    G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
         tau = gw_env%imag_time_points(i_t)
         CALL G_occ_vir(gw_env, tau, gw_env%fm_Gocc(i_t), occ=.TRUE., vir=.FALSE.)
         CALL G_occ_vir(gw_env, tau, gw_env%fm_Gvir(i_t), occ=.FALSE., vir=.TRUE.)
         CALL fm_to_local_tensor(gw_env%fm_Gocc(i_t), gw_env%mat_ao_ao%matrix, &
                                 gw_env%mat_ao_ao_tensor%matrix, t_2c_Gocc, gw_env)
         CALL fm_to_local_tensor(gw_env%fm_Gvir(i_t), gw_env%mat_ao_ao%matrix, &
                                 gw_env%mat_ao_ao_tensor%matrix, t_2c_Gvir, gw_env)

         ! every group has its own range of i_atoms and j_atoms; only deal with a
         ! single i_atom-j_atom pair simultaneously in a group to save memory
         DO i_atom = gw_env%atom_i_start, gw_env%atom_i_end, gw_env%n_atom_ij
            DO j_atom = gw_env%atom_j_start, gw_env%atom_j_end, gw_env%n_atom_ij

               ! 2. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
               CALL compute_3c_integrals(qs_env, gw_env, t_3c_for_Gocc, i_atom)
               CALL compute_3c_integrals(qs_env, gw_env, t_3c_for_Gvir, j_atom)

               ! 3. tensor operation M_λνP(iτ) = sum_µ (µν|P) G^occ_µλ(i|τ|,k=0)
               !                     N_νλQ(iτ) = sum_σ (σλ|Q) G^vir_σν(i|τ|,k=0)
               CALL G_times_3c(t_3c_for_Gocc, t_2c_Gocc, t_3c_x_Gocc, gw_env%eps_filter)
               CALL G_times_3c(t_3c_for_Gvir, t_2c_Gvir, t_3c_x_Gvir, gw_env%eps_filter)

               ! 4. reorder tensors
               CALL dbt_copy(t_3c_x_Gocc, t_3c_x_Gocc_2, move_data=.TRUE., order=[1, 3, 2])
               CALL dbt_copy(t_3c_x_Gvir, t_3c_x_Gvir_2, move_data=.TRUE.)

               ! 5. tensor operation χ_PQ(iτ,k=0) = sum_λν M_λνP(iτ) N_νλQ(iτ),
               CALL dbt_contract(alpha=gw_env%spin_degeneracy, &
                                 tensor_1=t_3c_x_Gocc_2, tensor_2=t_3c_x_Gvir_2, &
                                 beta=1.0_dp, tensor_3=gw_env%t_chi, &
                                 contract_1=[2, 3], notcontract_1=[1], map_1=[1], &
                                 contract_2=[2, 3], notcontract_2=[1], map_2=[2], &
                                 filter_eps=gw_env%eps_filter, move_data=.TRUE.)

            END DO ! j_atom
         END DO ! i_atom

         ! 6. communicate data of χ_PQ(iτ,k=0) in tensor gw_env%t_chi (which local in the
         !    subgroup) to the global dbcsr matrix mat_chi_Gamma_tau (which stores
         !    χ_PQ(iτ,k=0) for all time points)
         CALL local_dbt_to_global_mat(gw_env%t_chi, gw_env%mat_RI_RI_tensor%matrix, &
                                      mat_chi_Gamma_tau(i_t)%matrix, gw_env%para_env)

         IF (gw_env%unit_nr > 0) THEN
            WRITE (gw_env%unit_nr, '(T2,A,I7,A,I3,A,F12.1,A)') &
               'Computed χ(iτ,k=0) for time point', i_t, ' /', gw_env%num_time_freq_points, &
               ',     Execution time', m_walltime() - gw_env%t1, ' s'
         END IF

      END DO ! i_t

      IF (gw_env%unit_nr > 0) WRITE (gw_env%unit_nr, '(A)') ' '

      CALL dbt_destroy(t_2c_Gocc)
      CALL dbt_destroy(t_2c_Gvir)
      CALL dbt_destroy(t_3c_for_Gocc)
      CALL dbt_destroy(t_3c_for_Gvir)
      CALL dbt_destroy(t_3c_x_Gocc)
      CALL dbt_destroy(t_3c_x_Gvir)
      CALL dbt_destroy(t_3c_x_Gocc_2)
      CALL dbt_destroy(t_3c_x_Gvir_2)

      CALL timestop(handle)

   END SUBROUTINE compute_chi_Gamma

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
!> \param mat_chi_Gamma_tau ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE create_mat_for_chi(t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, t_3c_for_Gvir, t_3c_x_Gocc, &
                                 t_3c_x_Gvir, t_3c_x_Gocc_2, t_3c_x_Gvir_2, mat_chi_Gamma_tau, gw_env)
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_3c_for_Gocc, &
                                                            t_3c_for_Gvir, t_3c_x_Gocc, &
                                                            t_3c_x_Gvir, t_3c_x_Gocc_2, &
                                                            t_3c_x_Gvir_2
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_for_chi'

      INTEGER                                            :: handle, i_t

      CALL timeset(routineN, handle)

      CALL dbt_create(gw_env%t_G, t_2c_Gocc)
      CALL dbt_create(gw_env%t_G, t_2c_Gvir)
      CALL dbt_create(gw_env%t_RI_AO__AO, t_3c_for_Gocc)
      CALL dbt_create(gw_env%t_RI_AO__AO, t_3c_for_Gvir)
      CALL dbt_create(gw_env%t_RI_AO__AO, t_3c_x_Gocc)
      CALL dbt_create(gw_env%t_RI_AO__AO, t_3c_x_Gvir)
      CALL dbt_create(gw_env%t_RI__AO_AO, t_3c_x_Gocc_2)
      CALL dbt_create(gw_env%t_RI__AO_AO, t_3c_x_Gvir_2)

      NULLIFY (mat_chi_Gamma_tau)
      CALL dbcsr_allocate_matrix_set(mat_chi_Gamma_tau, gw_env%num_time_freq_points)

      DO i_t = 1, gw_env%num_time_freq_points
         ALLOCATE (mat_chi_Gamma_tau(i_t)%matrix)
         CALL dbcsr_create(mat_chi_Gamma_tau(i_t)%matrix, template=gw_env%mat_RI_RI%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_mat_for_chi

! **************************************************************************************************
!> \brief ...
!> \param fm_global ...
!> \param mat_global ...
!> \param mat_local ...
!> \param tensor ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE fm_to_local_tensor(fm_global, mat_global, mat_local, tensor, gw_env)

      TYPE(cp_fm_type)                                   :: fm_global
      TYPE(dbcsr_type)                                   :: mat_global, mat_local
      TYPE(dbt_type)                                     :: tensor
      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_to_local_tensor'

      INTEGER                                            :: handle
      TYPE(dbt_type)                                     :: tensor_tmp

      CALL timeset(routineN, handle)

      ! JW 2do: only distribute the matrix elements of tensor to subgroups which are actually
      !         needed in the subgroup because of iatom, jatom range

!      CALL dbcsr_set(gw_env%mat_ao_ao%matrix, 0.0_dp)
!      CALL dbt_clear(tensor)
!      CALL copy_fm_to_dbcsr(fm_global, gw_env%mat_ao_ao%matrix, keep_sparsity=.FALSE.)
!      CALL dbcsr_filter(gw_env%mat_ao_ao%matrix, gw_env%eps_filter)
!      CALL global_matrix_to_local_matrix(gw_env%mat_ao_ao%matrix, &
!                                         gw_env%mat_ao_ao_tensor%matrix, &
!                                         gw_env%para_env, gw_env%para_env_tensor%num_pe)
!      CALL dbt_create(gw_env%mat_ao_ao_tensor%matrix, tensor_tmp)
!      CALL dbt_copy_matrix_to_tensor(gw_env%mat_ao_ao_tensor%matrix, tensor_tmp)
!      CALL dbt_copy(tensor_tmp, tensor, move_data=.TRUE.)
!      CALL dbt_destroy(tensor_tmp)

      CALL dbt_clear(tensor)
      CALL copy_fm_to_dbcsr(fm_global, mat_global, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(mat_global, gw_env%eps_filter)
      CALL global_matrix_to_local_matrix(mat_global, mat_local, &
                                         gw_env%para_env, gw_env%para_env_tensor%num_pe)
      CALL dbt_create(mat_local, tensor_tmp)
      CALL dbt_copy_matrix_to_tensor(mat_local, tensor_tmp)
      CALL dbt_copy(tensor_tmp, tensor, move_data=.TRUE.)
      CALL dbt_destroy(tensor_tmp)

      CALL timestop(handle)

   END SUBROUTINE fm_to_local_tensor

! **************************************************************************************************
!> \brief ...
!> \param time ...
! **************************************************************************************************
   SUBROUTINE my_sleep(time)

      REAL(KIND=dp)                                      :: time

      INTEGER                                            :: a
      REAL(KIND=dp)                                      :: t_1

      t_1 = m_walltime()
      DO WHILE (m_walltime() - t_1 < time)
         a = 1 + 1
      END DO

   END SUBROUTINE my_sleep

! **************************************************************************************************
!> \brief ...
!> \param tensor ...
!> \param mat_tensor ...
!> \param mat_global ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE local_dbt_to_global_mat(tensor, mat_tensor, mat_global, para_env)

      TYPE(dbt_type)                                     :: tensor
      TYPE(dbcsr_type)                                   :: mat_tensor, mat_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'local_dbt_to_global_mat'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_copy_tensor_to_matrix(tensor, mat_tensor)
      CALL dbt_clear(tensor)
      CALL local_matrix_to_global_matrix(mat_tensor, mat_global, para_env)

      CALL timestop(handle)

   END SUBROUTINE local_dbt_to_global_mat

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
!> \param t_3c ...
!> \param i_atom ...
! **************************************************************************************************
   SUBROUTINE compute_3c_integrals(qs_env, gw_env, t_3c, i_atom)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(dbt_type)                                     :: t_3c
      INTEGER                                            :: i_atom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_integrals'

      INTEGER                                            :: handle
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_array

      CALL timeset(routineN, handle)

      ALLOCATE (t_3c_array(1, 1))
      CALL dbt_create(t_3c, t_3c_array(1, 1))

      CALL build_3c_integrals(t_3c_array, &
                              gw_env%eps_filter, &
                              qs_env, &
                              gw_env%nl_3c, &
                              int_eps=gw_env%eps_filter, &
                              basis_i=gw_env%basis_set_RI, &
                              basis_j=gw_env%basis_set_AO, &
                              basis_k=gw_env%basis_set_AO, &
                              potential_parameter=gw_env%ri_metric, &
                              bounds_j=[i_atom, i_atom], &
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
!> \param eps_filter ...
! **************************************************************************************************
   SUBROUTINE G_times_3c(t_3c_for_G, t_G, t_M, eps_filter)
      TYPE(dbt_type)                                     :: t_3c_for_G, t_G, t_M
      REAL(KIND=dp)                                      :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_times_3c'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=t_3c_for_G, &
                        tensor_2=t_G, &
                        beta=0.0_dp, &
                        tensor_3=t_M, &
                        contract_1=[3], notcontract_1=[1, 2], map_1=[1, 2], &
                        contract_2=[2], notcontract_2=[1], map_2=[3], &
                        filter_eps=eps_filter)

      CALL dbt_clear(t_3c_for_G)

      CALL timestop(handle)

   END SUBROUTINE G_times_3c

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_V_kp ...
! **************************************************************************************************
   SUBROUTINE compute_V_k_by_lattice_sum(gw_env, qs_env, fm_V_kp)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_k_by_lattice_sum'

      INTEGER                                            :: handle, ikp, re_im
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_RI_global
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      gw_env%t1 = m_walltime()

      NULLIFY (mat_V_kp)
      CALL dbcsr_allocate_matrix_set(mat_V_kp, gw_env%kpoints_chi_eps_W%nkp, 2)

      DO ikp = 1, gw_env%kpoints_chi_eps_W%nkp
         DO re_im = 1, 2

            ALLOCATE (mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_create(mat_V_kp(ikp, re_im)%matrix, template=gw_env%mat_RI_RI%matrix)
            CALL dbcsr_reserve_all_blocks(mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_set(mat_V_kp(ikp, re_im)%matrix, 0.0_dp)

         END DO ! re_im
      END DO ! ikp

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      ! 1. 2c Coulomb integrals for the first "original" k-point grid
      gw_env%kpoints_chi_eps_W%nkp_grid = gw_env%nkp_grid_chi_eps_W_orig
      CALL build_2c_coulomb_matrix_kp(mat_V_kp, &
                                      gw_env%kpoints_chi_eps_W, &
                                      basis_type="RI_AUX", &
                                      cell=cell, &
                                      particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=gw_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=1, &
                                      ikp_end=gw_env%nkp_chi_eps_W_orig)

      ! 2. 2c Coulomb integrals for the second "extrapolation" k-point grid
      gw_env%kpoints_chi_eps_W%nkp_grid = gw_env%nkp_grid_chi_eps_W_extra
      CALL build_2c_coulomb_matrix_kp(mat_V_kp, &
                                      gw_env%kpoints_chi_eps_W, &
                                      basis_type="RI_AUX", &
                                      cell=cell, &
                                      particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=gw_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=gw_env%nkp_chi_eps_W_orig + 1, &
                                      ikp_end=gw_env%kpoints_chi_eps_W%nkp)

      gw_env%kpoints_chi_eps_W%nkp_grid = gw_env%nkp_grid_chi_eps_W_orig

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      ! JW: THIS IS NOT GOOD AT ALL, IT WOULD BE MUCH BETTER TO ALLOCATE fm_V_kp ALREADY IN
      !     allocate_and_fill_matrices_and_arrays. THE PROBLEM WITH ALLOCATING fm_V_kp
      !     ALREADY IN allocate_and_fill_matrices_and_arrays IS THAT THERE IS SOME ERROR
      !     APPEARING AFTER THE CP2K EXECUTION.
      NULLIFY (fm_struct_RI_global)
      CALL cp_fm_struct_create(fm_struct_RI_global, context=blacs_env, &
                               nrow_global=gw_env%n_RI, ncol_global=gw_env%n_RI, &
                               para_env=para_env)

      ALLOCATE (fm_V_kp(gw_env%kpoints_chi_eps_W%nkp, 2))
      DO ikp = 1, gw_env%kpoints_chi_eps_W%nkp
         DO re_im = 1, 2
            CALL cp_fm_create(fm_V_kp(ikp, re_im), fm_struct_RI_global)
            CALL copy_dbcsr_to_fm(mat_V_kp(ikp, re_im)%matrix, fm_V_kp(ikp, re_im))
         END DO
      END DO

      CALL cp_fm_struct_release(fm_struct_RI_global)

      CALL dbcsr_deallocate_matrix_set(mat_V_kp)

      IF (gw_env%unit_nr > 0) THEN
         WRITE (gw_env%unit_nr, '(T2,A,A34,F12.1,A)') 'Computed V(k) from lattice sum,', &
            'Execution time', m_walltime() - gw_env%t1, ' s'
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_V_k_by_lattice_sum

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_V_kp ...
!> \param fm_V_sqrt_kp ...
!> \param fm_M_inv_V_sqrt_kp ...
! **************************************************************************************************
   SUBROUTINE compute_RI_metric_matrices(gw_env, qs_env, fm_V_kp, fm_V_sqrt_kp, fm_M_inv_V_sqrt_kp)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_kp, fm_V_sqrt_kp, fm_M_inv_V_sqrt_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_RI_metric_matrices'

      INTEGER                                            :: handle, ikp, re_im

      CALL timeset(routineN, handle)

      gw_env%t1 = m_walltime()

      ! get here M(k) and write it to fm_M_inv_V_sqrt_kp
      CALL RI_2c_integral_mat(qs_env, fm_M_inv_V_sqrt_kp, fm_V_kp(1, 1), &
                              gw_env%n_RI, gw_env%ri_metric, do_kpoints=.TRUE., &
                              kpoints=gw_env%kpoints_chi_eps_W, &
                              regularization_RI=gw_env%regularization_RI)

      ! calculate here inv(M(k)), sqrt(V(k)) and write result to fm_M_inv_V_sqrt_kp
      ! JW for later: do this in subgroups; this is also good because we need fm_M_inv_V_sqrt_kp
      !               in the subgroup for calculation eps(iw,k), W(iw,k)
      CALL inversion_of_M_and_mult_with_chol_dec_of_V(fm_M_inv_V_sqrt_kp, fm_V_kp, gw_env%n_RI, &
                                                      gw_env%kpoints_chi_eps_W, &
                                                      gw_env%eps_eigval_mat_RI)

      ! The routine inversion_of_M_and_mult_with_chol_dec_of_V overwrites fm_V_kp
      ! with its square root. Make this clear by the following renaming:
      ALLOCATE (fm_V_sqrt_kp(gw_env%kpoints_chi_eps_W%nkp, 2))
      DO ikp = 1, gw_env%kpoints_chi_eps_W%nkp
         DO re_im = 1, 2
            CALL cp_fm_create(fm_V_sqrt_kp(ikp, re_im), fm_V_kp(ikp, re_im)%matrix_struct)
            CALL cp_fm_to_fm(fm_V_kp(ikp, re_im), fm_V_sqrt_kp(ikp, re_im))
         END DO
      END DO

      ! JW 2do for later: Write restart for fm_M_inv_V_sqrt_kp

      CALL cp_fm_release(fm_V_kp)

      IF (gw_env%unit_nr > 0) THEN
         WRITE (gw_env%unit_nr, '(T2,A,A35,F12.1,A)') 'Computed V^0.5(k) and M^-1(k),', &
            'Execution time', m_walltime() - gw_env%t1, ' s'
         WRITE (gw_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_RI_metric_matrices

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE compute_W_MIC(gw_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W_MIC'

      INTEGER                                            :: handle, ikp, ikp_end, ikp_start, j_w
      TYPE(cp_cfm_type)                                  :: cfm_chi_ikp_freq_j, cfm_W_ikp_freq_j
      TYPE(cp_fm_type)                                   :: fm_chi_Gamma_freq_j, fm_W_MIC_freq_j
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M_inv_V_sqrt_kp, fm_V_kp, fm_V_sqrt_kp

      CALL timeset(routineN, handle)

      ! Compute V_PQ(k) = sum_R e^(ikR) <phi_P, cell 0 | 1/r | phi_Q, cell R>
      CALL compute_V_k_by_lattice_sum(gw_env, qs_env, fm_V_kp)

      CALL compute_RI_metric_matrices(gw_env, qs_env, fm_V_kp, fm_V_sqrt_kp, fm_M_inv_V_sqrt_kp)

      CALL create_fm_W_MIC_time(gw_env, fm_W_MIC_time, fm_V_sqrt_kp(1, 1))

      DO j_w = 1, gw_env%num_time_freq_points

         gw_env%t1 = m_walltime()

         ! 1. Fourier transformation of χ_PQ(iτ,k=0) to χ_PQ(iω_j,k=0)
         CALL compute_fm_chi_Gamma_freq_j(gw_env, fm_chi_Gamma_freq_j, j_w, &
                                          fm_M_inv_V_sqrt_kp(1, 1), mat_chi_Gamma_tau)

         ! 2. Replicate χ_PQ(iω_j,k=0) into kp_diag group (later, best with fm_copy_general)
         !    then fm_chi_Gamma_freq_j -> fm_chi_Gamma_freq_j_local

         ! JW needs to be adjusted for subgroups (every subgroup will have a single k-point
         !                                        or all kpoints)
         ! currently gw_env%ikp_local_start = 1, ikp_local_end = gw_env%kpoints_chi_eps_W%nkp
         ikp_start = gw_env%ikp_local_start
         ikp_end = gw_env%ikp_local_end

         DO ikp = ikp_start, ikp_end

            ! 3. Get χ_PQ(iω_j,k) from χ_PQ(iω_j,k=0) using the minimum image convention
            CALL cfm_ikp_from_fm_Gamma(cfm_chi_ikp_freq_j, fm_chi_Gamma_freq_j, &
                                       ikp, qs_env, gw_env%kpoints_chi_eps_W, "RI_AUX")

            ! 4. Remove all negative eigenvalues from χ_PQ(iω_j,k)
            CALL cp_cfm_power(cfm_chi_ikp_freq_j, threshold=0.0_dp, exponent=1.0_dp)

            ! 5. ε(iω_j,k) = Id - V^0.5(k)*M^-1(k)*χ(iω_j,k)*M^-1(k)*V^0.5(k)
            !    W(iω_j,k) = V^0.5(k)*(ε^-1(iω_j,k)-Id)*V^0.5(k)
            CALL compute_cfm_W_ikp_freq_j(gw_env, ikp, cfm_chi_ikp_freq_j, fm_V_sqrt_kp, &
                                          fm_M_inv_V_sqrt_kp, cfm_W_ikp_freq_j)

            ! 6. K-point integration (with extrapolation): W_PQ(iω_j,k) to W_PQ^MIC(iω_j)
            CALL compute_fm_W_MIC_freq_j(gw_env, qs_env, fm_W_MIC_freq_j, &
                                         cfm_W_ikp_freq_j, ikp, ikp_start)

         END DO ! ikp

         ! JW 2do: Communicate fm_W_MIC_freq_j from subgroup to global

         ! 7. Fourier trafo from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
         CALL Fourier_transform_w_to_t(gw_env, fm_W_MIC_time, fm_W_MIC_freq_j, j_w)

         ! 8. JW 2do: Sync and write restart file for W_PQ^MIC(iτ)

         ! 9. Print info
         IF (gw_env%unit_nr > 0) THEN
            WRITE (gw_env%unit_nr, '(T2,A,I7,A,I3,A,F12.1,A)') &
               'Computed W^MIC(iω) for freq point', j_w, ' /', gw_env%num_time_freq_points, &
               ',     Execution time', m_walltime() - gw_env%t1, ' s'
         END IF

         CALL cp_fm_release(fm_chi_Gamma_freq_j)

      END DO ! j_w

      ! 9. M^-1(k=0)*W^MIC(iτ)*M^-1(k=0)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(gw_env, qs_env, fm_W_MIC_time)

      CALL cp_fm_release(fm_M_inv_V_sqrt_kp)
      CALL cp_fm_release(fm_V_sqrt_kp)
      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      IF (gw_env%unit_nr > 0) WRITE (gw_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param fm_chi_Gamma_freq_j ...
!> \param j_w ...
!> \param fm_RI ...
!> \param mat_chi_Gamma_tau ...
! **************************************************************************************************
   SUBROUTINE compute_fm_chi_Gamma_freq_j(gw_env, fm_chi_Gamma_freq_j, j_w, fm_RI, mat_chi_Gamma_tau)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(cp_fm_type)                                   :: fm_chi_Gamma_freq_j
      INTEGER                                            :: j_w
      TYPE(cp_fm_type)                                   :: fm_RI
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_chi_Gamma_freq_j'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: freq_j, time_i, weight_ij

      CALL timeset(routineN, handle)

      CALL dbcsr_set(gw_env%mat_RI_RI%matrix, 0.0_dp)

      freq_j = gw_env%imag_freq_points(j_w)

      DO i_t = 1, gw_env%num_time_freq_points

         time_i = gw_env%imag_time_points(i_t)
         weight_ij = gw_env%weights_cos_t_to_w(j_w, i_t)

         ! actual Fourier transform
         CALL dbcsr_add(gw_env%mat_RI_RI%matrix, mat_chi_Gamma_tau(i_t)%matrix, &
                        1.0_dp, COS(time_i*freq_j)*weight_ij)

      END DO

      CALL cp_fm_create(fm_chi_Gamma_freq_j, fm_RI%matrix_struct)
      CALL copy_dbcsr_to_fm(gw_env%mat_RI_RI%matrix, fm_chi_Gamma_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_chi_Gamma_freq_j

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
!> \param cfm_ikp ...
!> \param fm_Gamma ...
!> \param ikp ...
!> \param qs_env ...
!> \param kpoints ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE cfm_ikp_from_fm_Gamma(cfm_ikp, fm_Gamma, ikp, qs_env, kpoints, basis_type)
      TYPE(cp_cfm_type)                                  :: cfm_ikp
      TYPE(cp_fm_type)                                   :: fm_Gamma
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*)                                   :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cfm_ikp_from_fm_Gamma'

      INTEGER :: col_global, handle, i, i_atom, i_atom_old, i_cell, i_mic_cell, i_row, j, j_atom, &
         j_atom_old, j_cell, j_col, n_bf, ncol_local, nrow_local, num_cells, row_global
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_bf
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL :: i_cell_is_the_minimum_image_cell
      REAL(KIND=dp)                                      :: abs_rab_cell_i, abs_rab_cell_j, arg
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, cell_vector_j, rab_cell_i, &
                                                            rab_cell_j
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      IF (.NOT. ASSOCIATED(cfm_ikp%local_data)) THEN
         CALL cp_cfm_create(cfm_ikp, fm_Gamma%matrix_struct)
      END IF
      CALL cp_cfm_set_all(cfm_ikp, z_zero)

      CALL cp_fm_get_info(matrix=fm_Gamma, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      ! get number of basis functions (bf) for different basis sets
      IF (basis_type == "ORB") THEN
         n_bf = qs_env%gw_env%n_ao
      ELSE IF (basis_type == "RI_AUX") THEN
         n_bf = qs_env%gw_env%n_RI
      ELSE
         CPABORT("Only ORB and RI_AUX basis implemented.")
      END IF

      ALLOCATE (atom_from_bf(n_bf))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_bf, n_bf, basis_type)

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      index_to_cell => kpoints%index_to_cell

      num_cells = SIZE(index_to_cell, 2)

      DO i_row = 1, nrow_local
         DO j_col = 1, ncol_local

            row_global = row_indices(i_row)
            col_global = col_indices(j_col)

            i_atom = atom_from_bf(row_global)
            j_atom = atom_from_bf(col_global)

            ! we only need to check for new MIC cell for new i_atom-j_atom pair
            IF (i_atom .NE. i_atom_old .OR. j_atom .NE. j_atom_old) THEN
               DO i_cell = 1, num_cells

                  ! only check nearest neigbors
                  IF (ANY(ABS(index_to_cell(1:3, i_cell)) > 1)) CYCLE

                  cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, i_cell), dp))

                  rab_cell_i(1:3) = pbc(particle_set(i_atom)%r(1:3), cell) - &
                                    (pbc(particle_set(j_atom)%r(1:3), cell) + cell_vector(1:3))
                  abs_rab_cell_i = SQRT(rab_cell_i(1)**2 + rab_cell_i(2)**2 + rab_cell_i(3)**2)

                  ! minimum image convention
                  i_cell_is_the_minimum_image_cell = .TRUE.
                  DO j_cell = 1, num_cells
                     cell_vector_j(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, j_cell), dp))
                     rab_cell_j(1:3) = pbc(particle_set(i_atom)%r(1:3), cell) - &
                                       (pbc(particle_set(j_atom)%r(1:3), cell) + cell_vector_j(1:3))
                     abs_rab_cell_j = SQRT(rab_cell_j(1)**2 + rab_cell_j(2)**2 + rab_cell_j(3)**2)

                     IF (abs_rab_cell_i > abs_rab_cell_j + 1.0E-6_dp) THEN
                        i_cell_is_the_minimum_image_cell = .FALSE.
                     END IF
                  END DO

                  IF (i_cell_is_the_minimum_image_cell) THEN
                     i_mic_cell = i_cell
                  END IF

               END DO ! i_cell
            END IF

            arg = REAL(index_to_cell(1, i_mic_cell), dp)*kpoints%xkp(1, ikp) + &
                  REAL(index_to_cell(2, i_mic_cell), dp)*kpoints%xkp(2, ikp) + &
                  REAL(index_to_cell(3, i_mic_cell), dp)*kpoints%xkp(3, ikp)

            i = i_row
            j = j_col

            cfm_ikp%local_data(i, j) = COS(twopi*arg)*fm_Gamma%local_data(i, j)*z_one + &
                                       SIN(twopi*arg)*fm_Gamma%local_data(i, j)*gaussi

            j_atom_old = j_atom
            i_atom_old = i_atom

         END DO ! j_col
      END DO ! i_row

      CALL timestop(handle)

   END SUBROUTINE cfm_ikp_from_fm_Gamma

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param ikp ...
!> \param cfm_chi_ikp_freq_j ...
!> \param fm_V_sqrt_kp ...
!> \param fm_M_inv_V_sqrt_kp ...
!> \param cfm_W_ikp_freq_j ...
! **************************************************************************************************
   SUBROUTINE compute_cfm_W_ikp_freq_j(gw_env, ikp, cfm_chi_ikp_freq_j, fm_V_sqrt_kp, &
                                       fm_M_inv_V_sqrt_kp, cfm_W_ikp_freq_j)

      TYPE(gw_type), POINTER                             :: gw_env
      INTEGER                                            :: ikp
      TYPE(cp_cfm_type)                                  :: cfm_chi_ikp_freq_j
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_V_sqrt_kp, fm_M_inv_V_sqrt_kp
      TYPE(cp_cfm_type)                                  :: cfm_W_ikp_freq_j

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cfm_W_ikp_freq_j'

      INTEGER                                            :: handle, info, n_RI
      TYPE(cp_cfm_type)                                  :: cfm_eps_ikp_freq_j, &
                                                            cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp, &
                                                            cfm_work

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_work, cfm_chi_ikp_freq_j%matrix_struct)
      CALL cp_cfm_create(cfm_M_inv_V_sqrt_ikp, cfm_chi_ikp_freq_j%matrix_struct)
      CALL cp_fm_to_cfm(fm_M_inv_V_sqrt_kp(ikp, 1), fm_M_inv_V_sqrt_kp(ikp, 2), cfm_M_inv_V_sqrt_ikp)
      n_RI = gw_env%n_RI

      ! 1. ε(iω_j,k) = Id - V^0.5(k)*M^-1(k)*χ(iω_j,k)*M^-1(k)*V^0.5(k)

      ! 1. a) work = χ(iω_j,k)*M^-1(k)*V^0.5(k)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, z_one, &
                         cfm_chi_ikp_freq_j, cfm_M_inv_V_sqrt_ikp, z_zero, cfm_work)
      CALL cp_cfm_release(cfm_chi_ikp_freq_j)

      ! 1. b) eps_work = V^0.5(k)*M^-1(k)*work
      CALL cp_cfm_create(cfm_eps_ikp_freq_j, cfm_work%matrix_struct)
      CALL parallel_gemm('C', 'N', n_RI, n_RI, n_RI, z_one, &
                         cfm_M_inv_V_sqrt_ikp, cfm_work, z_zero, cfm_eps_ikp_freq_j)
      CALL cp_cfm_release(cfm_M_inv_V_sqrt_ikp)

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
      CALL cp_cfm_create(cfm_V_sqrt_ikp, cfm_eps_ikp_freq_j%matrix_struct)
      CALL cp_fm_to_cfm(fm_V_sqrt_kp(ikp, 1), fm_V_sqrt_kp(ikp, 2), cfm_V_sqrt_ikp)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, z_one, cfm_eps_ikp_freq_j, cfm_V_sqrt_ikp, &
                         z_zero, cfm_work)

      ! 2. e) W(iw,k) = V^1/2(k)*work
      CALL cp_cfm_create(cfm_W_ikp_freq_j, cfm_work%matrix_struct)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, z_one, cfm_V_sqrt_ikp, cfm_work, &
                         z_zero, cfm_W_ikp_freq_j)

      CALL cp_cfm_release(cfm_work)
      CALL cp_cfm_release(cfm_eps_ikp_freq_j)
      CALL cp_cfm_release(cfm_V_sqrt_ikp)

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

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=cfm, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ! add 1 on the diagonal
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global) THEN
               cfm%local_data(iiB, jjB) = cfm%local_data(iiB, jjB) + alpha
            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE cfm_add_on_diag

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_j ...
!> \param cfm_W_ikp_freq_j ...
!> \param ikp ...
!> \param ikp_start ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_MIC_freq_j(gw_env, qs_env, fm_W_MIC_freq_j, cfm_W_ikp_freq_j, &
                                      ikp, ikp_start)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      TYPE(cp_cfm_type)                                  :: cfm_W_ikp_freq_j
      INTEGER                                            :: ikp, ikp_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_W_MIC_freq_j'

      INTEGER                                            :: handle, iatom, iatom_old, irow, jatom, &
                                                            jatom_old, jcol, n_RI, ncol_local, &
                                                            nrow_local, num_cells
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(KIND=dp)                                      :: contribution, weight_im, weight_re
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      IF (ikp == ikp_start) THEN
         CALL cp_fm_create(fm_W_MIC_freq_j, cfm_W_ikp_freq_j%matrix_struct)
         CALL cp_fm_set_all(fm_W_MIC_freq_j, 0.0_dp)
      END IF

      n_RI = gw_env%n_RI
      ALLOCATE (atom_from_RI_index(n_RI))
      CALL get_atom_index_from_basis_function_index(qs_env, atom_from_RI_index, n_RI, "RI_AUX")

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env, cell=cell, particle_set=particle_set)
      CALL get_cell(cell=cell, h=hmat)

      CALL cp_cfm_get_info(matrix=cfm_W_ikp_freq_j, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      CALL get_kpoint_info(gw_env%kpoints_chi_eps_W, xkp=xkp, wkp=wkp)
      index_to_cell => gw_env%kpoints_chi_eps_W%index_to_cell
      num_cells = SIZE(index_to_cell, 2)

      iatom_old = 0
      jatom_old = 0

      DO irow = 1, nrow_local
         DO jcol = 1, ncol_local

            iatom = atom_from_RI_index(row_indices(irow))
            jatom = atom_from_RI_index(col_indices(jcol))

            IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN

               ! symmetrize=.FALSE. necessary since we already have a symmetrized index_to_cell
               CALL compute_weight_re_im(weight_re, weight_im, &
                                         num_cells, iatom, jatom, xkp(1:3, ikp), wkp(ikp), &
                                         cell, index_to_cell, hmat, particle_set)

               iatom_old = iatom
               jatom_old = jatom

            END IF

            contribution = weight_re*REAL(cfm_W_ikp_freq_j%local_data(irow, jcol)) + &
                           weight_im*AIMAG(cfm_W_ikp_freq_j%local_data(irow, jcol))

            fm_W_MIC_freq_j%local_data(irow, jcol) = fm_W_MIC_freq_j%local_data(irow, jcol) &
                                                     + contribution

         END DO
      END DO

      CALL cp_cfm_release(cfm_W_ikp_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_MIC_freq_j

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param fm_W_MIC_time ...
!> \param fm_template ...
! **************************************************************************************************
   SUBROUTINE create_fm_W_MIC_time(gw_env, fm_W_MIC_time, fm_template)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type)                                   :: fm_template

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fm_W_MIC_time'

      INTEGER                                            :: handle, i_t

      CALL timeset(routineN, handle)

      ALLOCATE (fm_W_MIC_time(gw_env%num_time_freq_points))
      DO i_t = 1, gw_env%num_time_freq_points
         CALL cp_fm_create(fm_W_MIC_time(i_t), fm_template%matrix_struct)
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_fm_W_MIC_time

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param fm_W_MIC_time ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
! **************************************************************************************************
   SUBROUTINE Fourier_transform_w_to_t(gw_env, fm_W_MIC_time, fm_W_MIC_freq_j, j_w)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Fourier_transform_w_to_t'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: freq_j, time_i, weight_ij

      CALL timeset(routineN, handle)

      freq_j = gw_env%imag_freq_points(j_w)

      DO i_t = 1, gw_env%num_time_freq_points

         time_i = gw_env%imag_time_points(i_t)
         weight_ij = gw_env%weights_cos_w_to_t(i_t, j_w)

         ! actual Fourier transform
         CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_W_MIC_time(i_t), &
                                  beta=weight_ij*COS(time_i*freq_j), matrix_b=fm_W_MIC_freq_j)

      END DO

      CALL cp_fm_release(fm_W_MIC_freq_j)

      CALL timestop(handle)

   END SUBROUTINE Fourier_transform_w_to_t

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE multiply_fm_W_MIC_time_with_Minv_Gamma(gw_env, qs_env, fm_W_MIC_time)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'multiply_fm_W_MIC_time_with_Minv_Gamma'

      INTEGER                                            :: handle, i_t, n_RI, ndep
      TYPE(cp_fm_type)                                   :: fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Minv_Gamma

      CALL timeset(routineN, handle)

      n_RI = gw_env%n_RI

      CALL cp_fm_create(fm_work, fm_W_MIC_time(1)%matrix_struct)

      ! compute Gamma-only RI-metric matrix M(k=0)
      ! JW COMMENT: in the old GW implementation, there was no regularization_RI
      !             for M^-1(k=0)... Maybe this is an improvement here?
      CALL RI_2c_integral_mat(qs_env, fm_Minv_Gamma, fm_W_MIC_time(1), n_RI, &
                              gw_env%ri_metric, do_kpoints=.FALSE., &
                              regularization_RI=gw_env%regularization_RI)

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

!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_W_MIC_time ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(gw_env, qs_env, fm_W_MIC_time, fm_Sigma_x_Gamma)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type)                                   :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, i_atom, j_atom
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma
      TYPE(dbt_type)                                     :: t_2c_D, t_2c_Sigma_x, t_2c_V, &
                                                            t_3c_for_D, t_3c_x_V

      CALL timeset(routineN, handle)

      gw_env%t1 = m_walltime()

      CALL dbt_create(gw_env%t_G, t_2c_D)
      CALL dbt_create(gw_env%t_W, t_2c_V)
      CALL dbt_create(gw_env%t_G, t_2c_Sigma_x)
      CALL dbt_create(gw_env%t_RI_AO__AO, t_3c_for_D)
      CALL dbt_create(gw_env%t_RI__AO_AO, t_3c_x_V)
      CALL dbcsr_create(mat_Sigma_x_Gamma, template=gw_env%mat_ao_ao%matrix)

      ! 1. Compute density matrix D_µν
      CALL G_occ_vir(gw_env, 0.0_dp, gw_env%fm_work_mo(2), occ=.TRUE., vir=.FALSE.)
      CALL fm_to_local_tensor(gw_env%fm_work_mo(2), gw_env%mat_ao_ao%matrix, &
                              gw_env%mat_ao_ao_tensor%matrix, t_2c_D, gw_env)

      ! 2. Compute truncated Coulomb operator matrix V^tr(k=0) (cutoff rad: cellsize/2)
      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, fm_W_MIC_time(1), gw_env%n_RI, &
                              gw_env%trunc_coulomb, do_kpoints=.FALSE., &
                              regularization_RI=gw_env%regularization_RI)

      ! 3. Compute M^-1(k=0) and get M^-1(k=0)*V^tr(k=0)*M^-1(k=0)
      CALL multiply_fm_W_MIC_time_with_Minv_Gamma(gw_env, qs_env, fm_Vtr_Gamma(:, 1))

      CALL fm_to_local_tensor(fm_Vtr_Gamma(1, 1), gw_env%mat_RI_RI%matrix, &
                              gw_env%mat_RI_RI_tensor%matrix, t_2c_V, gw_env)

      ! i_atom, j_atom loop as in compute_Sigma_c
      DO i_atom = gw_env%atom_i_start, gw_env%atom_i_end, gw_env%n_atom_ij
         DO j_atom = gw_env%atom_j_start, gw_env%atom_j_end, gw_env%n_atom_ij

            ! 4. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
            ! 5. M_νσQ(iτ) = sum_P (νσ|P) (M^-1(k=0)*V^tr(k=0)*M^-1(k=0))_PQ(iτ)
            CALL compute_3c_and_contract_W(qs_env, gw_env, i_atom, j_atom, &
                                           t_3c_for_D, t_3c_x_V, t_2c_V)

            ! 6. tensor operations with D and computation of Σ^x
            !    Σ^x_λσ(k=0) = sum_νQ M_νσQ(iτ) sum_µ (λµ|Q) D_µν
            CALL contract_to_Sigma(t_3c_for_D, t_2c_D, t_3c_x_V, t_2c_Sigma_x, &
                                   gw_env%eps_filter, do_occ=.TRUE., do_virt=.FALSE.)

         END DO ! j_atom
      END DO ! i_atom

      CALL local_dbt_to_global_mat(t_2c_Sigma_x, gw_env%mat_ao_ao_tensor%matrix, &
                                   mat_Sigma_x_Gamma, gw_env%para_env)

      CALL cp_fm_create(fm_Sigma_x_Gamma, gw_env%fm_s_Gamma%matrix_struct)
      CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma)

      IF (gw_env%unit_nr > 0) THEN
         WRITE (gw_env%unit_nr, '(T2,A,T53,A,F12.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - gw_env%t1, ' s'
         WRITE (gw_env%unit_nr, '(A)') ' '
      END IF

      CALL dbcsr_release(mat_Sigma_x_Gamma)
      CALL dbt_destroy(t_2c_D)
      CALL dbt_destroy(t_2c_V)
      CALL dbt_destroy(t_2c_Sigma_x)
      CALL dbt_destroy(t_3c_x_V)
      CALL dbt_destroy(t_3c_for_D)
      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_W_MIC_time ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c(gw_env, qs_env, fm_W_MIC_time, fm_Sigma_c_Gamma_time)
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_atom, i_t, j_atom
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, &
                                                            t_2c_Sigma_neg_tau, &
                                                            t_2c_Sigma_pos_tau, t_2c_W, &
                                                            t_3c_for_G, t_3c_x_W

      CALL timeset(routineN, handle)

      CALL create_mat_for_Sigma_c(gw_env, t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                                  t_2c_Sigma_pos_tau, t_3c_for_G, t_3c_x_W, &
                                  mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      DO i_t = 1, gw_env%num_time_freq_points

         gw_env%t1 = m_walltime()

         ! fm G^occ, G^vir and W to local tensor
         CALL fm_to_local_tensor(gw_env%fm_Gocc(i_t), gw_env%mat_ao_ao%matrix, &
                                 gw_env%mat_ao_ao_tensor%matrix, t_2c_Gocc, gw_env)
         CALL fm_to_local_tensor(gw_env%fm_Gvir(i_t), gw_env%mat_ao_ao%matrix, &
                                 gw_env%mat_ao_ao_tensor%matrix, t_2c_Gvir, gw_env)
         CALL fm_to_local_tensor(fm_W_MIC_time(i_t), gw_env%mat_RI_RI%matrix, &
                                 gw_env%mat_RI_RI_tensor%matrix, t_2c_W, gw_env)

         ! every group has its own range of i_atoms and j_atoms; only deal with a
         ! single i_atom-j_atom pair simultaneously in a group to save memory
         DO i_atom = gw_env%atom_i_start, gw_env%atom_i_end, gw_env%n_atom_ij
            DO j_atom = gw_env%atom_j_start, gw_env%atom_j_end, gw_env%n_atom_ij

               ! 1. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
               ! 2. tensor operation M_νσQ(iτ) = sum_P (νσ|P) W^MIC_PQ(iτ)
               CALL compute_3c_and_contract_W(qs_env, gw_env, i_atom, j_atom, &
                                              t_3c_for_G, t_3c_x_W, t_2c_W)

               ! 3. Σ_λσ(iτ,k=0) = sum_νQ M_νσQ(iτ) sum_µ (λµ|Q) G^occ_µν(i|τ|) for τ < 0
               !    (recall M_νσQ(iτ) = M_νσQ(-iτ) because W^MIC_PQ(iτ) = W^MIC_PQ(-iτ) )
               CALL contract_to_Sigma(t_3c_for_G, t_2c_Gocc, t_3c_x_W, t_2c_Sigma_neg_tau, &
                                      gw_env%eps_filter, do_occ=.TRUE., do_virt=.FALSE.)
               !    Σ_λσ(iτ,k=0) = sum_νQ M_νσQ(iτ) sum_µ (λµ|Q) G^vir_µν(i|τ|) for τ > 0
               CALL contract_to_Sigma(t_3c_for_G, t_2c_Gvir, t_3c_x_W, t_2c_Sigma_pos_tau, &
                                      gw_env%eps_filter, do_occ=.FALSE., do_virt=.TRUE.)

            END DO ! j_atom
         END DO ! i_atom

         ! 4. communicate data tensor t_2c_Sigma (which is local in the subgroup)
         !    to the global dbcsr matrix mat_Sigma_pos/neg_tau (which stores Σ for all iτ)
         CALL local_dbt_to_global_mat(t_2c_Sigma_neg_tau, gw_env%mat_ao_ao_tensor%matrix, &
                                      mat_Sigma_neg_tau(i_t)%matrix, gw_env%para_env)
         CALL local_dbt_to_global_mat(t_2c_Sigma_pos_tau, gw_env%mat_ao_ao_tensor%matrix, &
                                      mat_Sigma_pos_tau(i_t)%matrix, gw_env%para_env)

         IF (gw_env%unit_nr > 0) THEN
            WRITE (gw_env%unit_nr, '(T2,A,I5,A,I3,A,F12.1,A)') &
               'Computed Σ^c(iτ,k=0) for time point', i_t, ' /', gw_env%num_time_freq_points, &
               ',     Execution time', m_walltime() - gw_env%t1, ' s'
         END IF

      END DO ! i_t

      IF (gw_env%unit_nr > 0) WRITE (gw_env%unit_nr, '(A)') ' '

      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, gw_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL destroy_mat_Sigma_c(t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                               t_2c_Sigma_pos_tau, t_3c_for_G, t_3c_x_W, fm_W_MIC_time, &
                               mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param t_2c_Gocc ...
!> \param t_2c_Gvir ...
!> \param t_2c_W ...
!> \param t_2c_Sigma_neg_tau ...
!> \param t_2c_Sigma_pos_tau ...
!> \param t_3c_for_G ...
!> \param t_3c_x_W ...
!> \param mat_Sigma_neg_tau ...
!> \param mat_Sigma_pos_tau ...
! **************************************************************************************************
   SUBROUTINE create_mat_for_Sigma_c(gw_env, t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                                     t_2c_Sigma_pos_tau, t_3c_for_G, t_3c_x_W, &
                                     mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_2c_W, &
                                                            t_2c_Sigma_neg_tau, &
                                                            t_2c_Sigma_pos_tau, t_3c_for_G, &
                                                            t_3c_x_W
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_Sigma_neg_tau, mat_Sigma_pos_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_for_Sigma_c'

      INTEGER                                            :: handle, i_t

      CALL timeset(routineN, handle)

      CALL dbt_create(gw_env%t_G, t_2c_Gocc)
      CALL dbt_create(gw_env%t_G, t_2c_Gvir)
      CALL dbt_create(gw_env%t_W, t_2c_W)
      CALL dbt_create(gw_env%t_G, t_2c_Sigma_neg_tau)
      CALL dbt_create(gw_env%t_G, t_2c_Sigma_pos_tau)
      CALL dbt_create(gw_env%t_RI_AO__AO, t_3c_for_G)
      CALL dbt_create(gw_env%t_RI__AO_AO, t_3c_x_W)

      NULLIFY (mat_Sigma_neg_tau, mat_Sigma_pos_tau)
      CALL dbcsr_allocate_matrix_set(mat_Sigma_neg_tau, gw_env%num_time_freq_points)
      CALL dbcsr_allocate_matrix_set(mat_Sigma_pos_tau, gw_env%num_time_freq_points)

      DO i_t = 1, gw_env%num_time_freq_points
         ALLOCATE (mat_Sigma_neg_tau(i_t)%matrix)
         ALLOCATE (mat_Sigma_pos_tau(i_t)%matrix)
         CALL dbcsr_create(mat_Sigma_neg_tau(i_t)%matrix, template=gw_env%mat_ao_ao%matrix)
         CALL dbcsr_create(mat_Sigma_pos_tau(i_t)%matrix, template=gw_env%mat_ao_ao%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_mat_for_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
!> \param i_atom ...
!> \param j_atom ...
!> \param t_3c_for_G ...
!> \param t_3c_x_W ...
!> \param t_2c_W ...
! **************************************************************************************************
   SUBROUTINE compute_3c_and_contract_W(qs_env, gw_env, i_atom, j_atom, &
                                        t_3c_for_G, t_3c_x_W, t_2c_W)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      INTEGER                                            :: i_atom, j_atom
      TYPE(dbt_type)                                     :: t_3c_for_G, t_3c_x_W, t_2c_W

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_3c_and_contract_W'

      INTEGER                                            :: handle
      TYPE(dbt_type)                                     :: t_3c_for_W, t_3c_x_W_tmp

      CALL timeset(routineN, handle)

      CALL dbt_create(gw_env%t_RI__AO_AO, t_3c_x_W_tmp)
      CALL dbt_create(gw_env%t_RI__AO_AO, t_3c_for_W)

      ! 1. compute 3-center integrals (µν|P) ("|": truncated Coulomb operator)
      CALL compute_3c_integrals(qs_env, gw_env, t_3c_for_W, i_atom)
      CALL compute_3c_integrals(qs_env, gw_env, t_3c_for_G, j_atom)

      ! 2. tensor operation M_νσQ(iτ) = sum_P (νσ|P) W^MIC_PQ(iτ)
      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=t_2c_W, &
                        tensor_2=t_3c_for_W, &
                        beta=0.0_dp, &
                        tensor_3=t_3c_x_W_tmp, &
                        contract_1=[2], notcontract_1=[1], map_1=[1], &
                        contract_2=[1], notcontract_2=[2, 3], map_2=[2, 3], &
                        filter_eps=gw_env%eps_filter)

      ! 3. reorder tensor
      CALL dbt_copy(t_3c_x_W_tmp, t_3c_x_W, order=[1, 2, 3], move_data=.TRUE.)

      CALL dbt_destroy(t_3c_x_W_tmp)
      CALL dbt_destroy(t_3c_for_W)

      CALL timestop(handle)

   END SUBROUTINE compute_3c_and_contract_W

! **************************************************************************************************
!> \brief ...
!> \param t_3c_for_G ...
!> \param t_2c_G ...
!> \param t_3c_x_W_2 ...
!> \param t_2c_Sigma ...
!> \param eps_filter ...
!> \param do_occ ...
!> \param do_virt ...
! **************************************************************************************************
   SUBROUTINE contract_to_Sigma(t_3c_for_G, t_2c_G, t_3c_x_W_2, t_2c_Sigma, eps_filter, &
                                do_occ, do_virt)
      TYPE(dbt_type)                                     :: t_3c_for_G, t_2c_G, t_3c_x_W_2, &
                                                            t_2c_Sigma
      REAL(KIND=dp)                                      :: eps_filter
      LOGICAL                                            :: do_occ, do_virt

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_to_Sigma'

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: sign_Sigma
      TYPE(dbt_type)                                     :: t_3c_x_G, t_3c_x_G_2

      CALL timeset(routineN, handle)

      CPASSERT(do_occ .EQV. (.NOT. do_virt))
      IF (do_occ) sign_Sigma = -1.0_dp
      IF (do_virt) sign_Sigma = 1.0_dp

      CALL dbt_create(t_3c_for_G, t_3c_x_G)
      CALL dbt_create(t_3c_for_G, t_3c_x_G_2)

      CALL dbt_contract(alpha=1.0_dp, &
                        tensor_1=t_2c_G, &
                        tensor_2=t_3c_for_G, &
                        beta=0.0_dp, &
                        tensor_3=t_3c_x_G, &
                        contract_1=[2], notcontract_1=[1], map_1=[3], &
                        contract_2=[3], notcontract_2=[1, 2], map_2=[1, 2], &
                        filter_eps=eps_filter)

      CALL dbt_copy(t_3c_x_G, t_3c_x_G_2, order=[1, 3, 2], move_data=.TRUE.)

      CALL dbt_contract(alpha=sign_Sigma, &
                        tensor_1=t_3c_x_W_2, &
                        tensor_2=t_3c_x_G_2, &
                        beta=1.0_dp, &
                        tensor_3=t_2c_Sigma, &
                        contract_1=[1, 2], notcontract_1=[3], map_1=[1], &
                        contract_2=[1, 2], notcontract_2=[3], map_2=[2], &
                        filter_eps=eps_filter)

      CALL dbt_destroy(t_3c_x_G)
      CALL dbt_destroy(t_3c_x_G_2)

      CALL timestop(handle)

   END SUBROUTINE contract_to_Sigma

! **************************************************************************************************
!> \brief ...
!> \param fm_Sigma_c_Gamma_time ...
!> \param gw_env ...
!> \param mat_Sigma_pos_tau ...
!> \param mat_Sigma_neg_tau ...
! **************************************************************************************************
   SUBROUTINE fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, gw_env, &
                                         mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma_time
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_Sigma_pos_tau, mat_Sigma_neg_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_fm_Sigma_c_Gamma_time'

      INTEGER                                            :: handle, i_t, pos_neg

      CALL timeset(routineN, handle)

      ALLOCATE (fm_Sigma_c_Gamma_time(gw_env%num_time_freq_points, 2))
      DO i_t = 1, gw_env%num_time_freq_points
         DO pos_neg = 1, 2
            CALL cp_fm_create(fm_Sigma_c_Gamma_time(i_t, pos_neg), gw_env%fm_s_Gamma%matrix_struct)
         END DO
         CALL copy_dbcsr_to_fm(mat_Sigma_pos_tau(i_t)%matrix, fm_Sigma_c_Gamma_time(i_t, 1))
         CALL copy_dbcsr_to_fm(mat_Sigma_neg_tau(i_t)%matrix, fm_Sigma_c_Gamma_time(i_t, 2))
      END DO

      CALL timestop(handle)

   END SUBROUTINE fill_fm_Sigma_c_Gamma_time

! **************************************************************************************************
!> \brief ...
!> \param t_2c_Gocc ...
!> \param t_2c_Gvir ...
!> \param t_2c_W ...
!> \param t_2c_Sigma_neg_tau ...
!> \param t_2c_Sigma_pos_tau ...
!> \param t_3c_for_G ...
!> \param t_3c_x_W ...
!> \param fm_W_MIC_time ...
!> \param mat_Sigma_neg_tau ...
!> \param mat_Sigma_pos_tau ...
! **************************************************************************************************
   SUBROUTINE destroy_mat_Sigma_c(t_2c_Gocc, t_2c_Gvir, t_2c_W, t_2c_Sigma_neg_tau, &
                                  t_2c_Sigma_pos_tau, t_3c_for_G, t_3c_x_W, fm_W_MIC_time, &
                                  mat_Sigma_neg_tau, mat_Sigma_pos_tau)

      TYPE(dbt_type)                                     :: t_2c_Gocc, t_2c_Gvir, t_2c_W, &
                                                            t_2c_Sigma_neg_tau, &
                                                            t_2c_Sigma_pos_tau, t_3c_for_G, &
                                                            t_3c_x_W
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_Sigma_neg_tau, mat_Sigma_pos_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_mat_Sigma_c'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_destroy(t_2c_Gocc)
      CALL dbt_destroy(t_2c_Gvir)
      CALL dbt_destroy(t_2c_W)
      CALL dbt_destroy(t_2c_Sigma_neg_tau)
      CALL dbt_destroy(t_2c_Sigma_pos_tau)
      CALL dbt_destroy(t_3c_for_G)
      CALL dbt_destroy(t_3c_x_W)
      CALL cp_fm_release(fm_W_MIC_time)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL timestop(handle)

   END SUBROUTINE destroy_mat_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param Sigma_c_n_time ...
!> \param Sigma_c_n_freq ...
! **************************************************************************************************
   SUBROUTINE time_to_freq(gw_env, Sigma_c_n_time, Sigma_c_n_freq)
      TYPE(gw_type), POINTER                             :: gw_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_n_time, Sigma_c_n_freq

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'time_to_freq'

      INTEGER                                            :: handle, i_t, j_w, n_occ
      REAL(KIND=dp)                                      :: freq_j, time_i, w_cos_ij, w_sin_ij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Sigma_c_n_cos_time, Sigma_c_n_sin_time

      CALL timeset(routineN, handle)

      ALLOCATE (Sigma_c_n_cos_time(gw_env%n_ao, gw_env%num_time_freq_points))
      ALLOCATE (Sigma_c_n_sin_time(gw_env%n_ao, gw_env%num_time_freq_points))

      Sigma_c_n_cos_time(:, :) = 0.5_dp*(Sigma_c_n_time(:, :, 1) + Sigma_c_n_time(:, :, 2))
      Sigma_c_n_sin_time(:, :) = 0.5_dp*(Sigma_c_n_time(:, :, 1) - Sigma_c_n_time(:, :, 2))

      Sigma_c_n_freq(:, :, :) = 0.0_dp

      DO i_t = 1, gw_env%num_time_freq_points

         DO j_w = 1, gw_env%num_time_freq_points

            freq_j = gw_env%imag_freq_points(j_w)
            time_i = gw_env%imag_time_points(i_t)
            ! integration weights for cosine and sine transform
            w_cos_ij = gw_env%weights_cos_t_to_w(j_w, i_t)*COS(freq_j*time_i)
            w_sin_ij = gw_env%weights_sin_t_to_w(j_w, i_t)*SIN(freq_j*time_i)

            ! 1. Re(Σ^c_nn(k_i,iω)) from cosine transform
            Sigma_c_n_freq(:, j_w, 1) = Sigma_c_n_freq(:, j_w, 1) + &
                                        w_cos_ij*Sigma_c_n_cos_time(:, i_t)

            ! 2. Im(Σ^c_nn(k_i,iω)) from sine transform
            Sigma_c_n_freq(:, j_w, 2) = Sigma_c_n_freq(:, j_w, 2) + &
                                        w_sin_ij*Sigma_c_n_sin_time(:, i_t)

         END DO

      END DO

      ! for occupied levels, we need the correlation self-energy for negative omega.
      ! Therefore, weight_sin should be computed with -omega, which results in an
      ! additional minus for the imaginary part:
      n_occ = gw_env%n_occ(1)
      Sigma_c_n_freq(1:n_occ, :, 2) = -Sigma_c_n_freq(1:n_occ, :, 2)

      CALL timestop(handle)

   END SUBROUTINE time_to_freq

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(gw_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_Sigma_x_Gamma
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER                                            :: handle, ikp, j_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n, V_xc_ikp_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_freq, Sigma_c_ikp_n_time
      TYPE(cp_cfm_type)                                  :: cfm_ks_ikp, cfm_mos_ikp, cfm_s_ikp, &
                                                            cfm_Sigma_x_ikp, cfm_work_ikp
      TYPE(cp_fm_type)                                   :: fm_V_xc_Gamma

      CALL timeset(routineN, handle)

      ! JW TODO: kpoints for kpoints_Sigma also from bandstructure !!!
      CALL cp_cfm_create(cfm_mos_ikp, gw_env%fm_s_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_work_ikp, gw_env%fm_s_Gamma%matrix_struct)
      CALL compute_fm_V_xc_Gamma(qs_env, gw_env, fm_V_xc_Gamma)
      ! JW TODO: fully distribute these arrays at given time; also eigenvalues in gw_env
      ALLOCATE (V_xc_ikp_n(gw_env%n_ao), Sigma_x_ikp_n(gw_env%n_ao))
      ALLOCATE (Sigma_c_ikp_n_time(gw_env%n_ao, gw_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_ikp_n_freq(gw_env%n_ao, gw_env%num_time_freq_points, 2))

      DO ikp = 1, gw_env%kpoints_Sigma%nkp

         ! 1. get S_µν(k_i) from S_µν(k=0)
         CALL cfm_ikp_from_fm_Gamma(cfm_ks_ikp, gw_env%fm_ks_Gamma(1), &
                                    ikp, qs_env, gw_env%kpoints_Sigma, "ORB")

         ! 2. get H^KS_µν(k_i) from H^KS_µν(k=0)
         CALL cfm_ikp_from_fm_Gamma(cfm_s_ikp, gw_env%fm_s_Gamma, &
                                    ikp, qs_env, gw_env%kpoints_Sigma, "ORB")

         ! 3. Diagonalize (Roothaan-Hall): H_KS(k_i)*C(k_i) = S(k_i)*C(k_i)*ϵ(k_i)
         CALL cp_cfm_geeig(cfm_ks_ikp, cfm_s_ikp, cfm_mos_ikp, &
                           gw_env%eigenvalues_scf(:, ikp), cfm_work_ikp)

         ! 4. V^xc_µν(k=0) -> V^xc_µν(k_i) -> V^xc_nn(k_i)
         CALL to_ikp_and_mo(V_xc_ikp_n, fm_V_xc_Gamma, ikp, qs_env, gw_env, cfm_mos_ikp)

         ! 5. Σ^x_µν(k=0) -> Σ^x_µν(k_i) -> Σ^x_nn(k_i)
         CALL to_ikp_and_mo(Sigma_x_ikp_n, fm_Sigma_x_Gamma, ikp, qs_env, gw_env, cfm_mos_ikp)

         ! 6. Σ^c_µν(k=0,+/-i|τ_j|) -> Σ^c_µν(k_i,+/-i|τ_j|) -> Σ^c_nn(k_i,+/-i|τ_j|)
         DO j_t = 1, gw_env%num_time_freq_points
            CALL to_ikp_and_mo(Sigma_c_ikp_n_time(:, j_t, 1), fm_Sigma_c_Gamma_time(j_t, 1), &
                               ikp, qs_env, gw_env, cfm_mos_ikp)
            CALL to_ikp_and_mo(Sigma_c_ikp_n_time(:, j_t, 2), fm_Sigma_c_Gamma_time(j_t, 2), &
                               ikp, qs_env, gw_env, cfm_mos_ikp)
         END DO

         ! 7. Σ^c_nn(k_i,iτ) -> Σ^c_nn(k_i,iω)
         CALL time_to_freq(gw_env, Sigma_c_ikp_n_time, Sigma_c_ikp_n_freq)

         ! 8. Analytic continuation Σ^c_nn(k_i,iω) -> Σ^c_nn(k_i,ϵ) and
         !    ϵ_nk_i^GW = ϵ_nk_i^DFT + Σ^c_nn(k_i,ϵ) + Σ^x_nn(k_i) - v^xc_nn(k_i)
         CALL analyt_conti_and_print(gw_env, Sigma_c_ikp_n_freq, Sigma_x_ikp_n, V_xc_ikp_n, ikp)

      END DO

      CALL get_VBM_CBM_bandgaps_and_print(gw_env)

      CALL cp_fm_release(fm_Sigma_x_Gamma)
      CALL cp_fm_release(fm_Sigma_c_Gamma_time)
      CALL cp_fm_release(fm_V_xc_Gamma)
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
!> \param gw_env ...
!> \param cfm_mos_ikp ...
! **************************************************************************************************
   SUBROUTINE to_ikp_and_mo(array_ikp_n, fm_Gamma, ikp, qs_env, gw_env, cfm_mos_ikp)

      REAL(KIND=dp), DIMENSION(:)                        :: array_ikp_n
      TYPE(cp_fm_type)                                   :: fm_Gamma
      INTEGER                                            :: ikp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(cp_cfm_type)                                  :: cfm_mos_ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'to_ikp_and_mo'

      INTEGER                                            :: handle, nmo
      TYPE(cp_cfm_type)                                  :: cfm_ikp_ao, cfm_ikp_mo, cfm_tmp
      TYPE(cp_fm_type)                                   :: fm_ikp_mo_re

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_ikp_ao, fm_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_tmp, fm_Gamma%matrix_struct)
      CALL cp_cfm_create(cfm_ikp_mo, fm_Gamma%matrix_struct)
      CALL cp_fm_create(fm_ikp_mo_re, fm_Gamma%matrix_struct)

      ! get cfm_µν(k_i) from fm_µν(k=0)
      CALL cfm_ikp_from_fm_Gamma(cfm_ikp_ao, fm_Gamma, ikp, qs_env, gw_env%kpoints_Sigma, "ORB")

      nmo = gw_env%n_ao
      CALL parallel_gemm('N', 'N', nmo, nmo, nmo, z_one, cfm_ikp_ao, cfm_mos_ikp, z_zero, cfm_tmp)
      CALL parallel_gemm('C', 'N', nmo, nmo, nmo, z_one, cfm_mos_ikp, cfm_tmp, z_zero, cfm_ikp_mo)

      CALL cp_cfm_to_fm(cfm_ikp_mo, fm_ikp_mo_re)
      CALL cp_fm_get_diag(fm_ikp_mo_re, array_ikp_n)

      CALL cp_cfm_release(cfm_ikp_ao)
      CALL cp_cfm_release(cfm_tmp)
      CALL cp_cfm_release(cfm_ikp_mo)
      CALL cp_fm_release(fm_ikp_mo_re)

      CALL timestop(handle)

   END SUBROUTINE to_ikp_and_mo

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gw_env ...
!> \param fm_V_xc_Gamma ...
! **************************************************************************************************
   SUBROUTINE compute_fm_V_xc_Gamma(qs_env, gw_env, fm_V_xc_Gamma)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(gw_type), POINTER                             :: gw_env
      TYPE(cp_fm_type)                                   :: fm_V_xc_Gamma

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_V_xc_Gamma'

      INTEGER                                            :: handle, ispin, myfun
      REAL(KIND=dp)                                      :: energy_ex, energy_exc, energy_total
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_ks_without_v_xc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input, energy=energy, dft_control=dft_control)

      ! previously, dft_control%nimages set to # neighbor cells, revert for Γ-only KS matrix
      dft_control%nimages = 1

      ! we need to reset XC functional, therefore, get XC input
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=myfun)
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=xc_none)

      ! save the energy before the energy gets updated
      energy_total = energy%total
      energy_exc = energy%exc
      energy_ex = energy%ex

      NULLIFY (mat_ks_without_v_xc)
      CALL dbcsr_allocate_matrix_set(mat_ks_without_v_xc, gw_env%n_spin)

      DO ispin = 1, gw_env%n_spin
         ALLOCATE (mat_ks_without_v_xc(ispin)%matrix)
         CALL dbcsr_create(mat_ks_without_v_xc(ispin)%matrix, template=gw_env%mat_ao_ao%matrix)
      END DO

      ! calculate KS-matrix without XC
      CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE., &
                                        ext_ks_matrix=mat_ks_without_v_xc)

      ! transfer dbcsr matrix to fm
      CALL cp_fm_create(fm_V_xc_Gamma, gw_env%fm_s_Gamma%matrix_struct)
      CALL copy_dbcsr_to_fm(mat_ks_without_v_xc(1)%matrix, fm_V_xc_Gamma)

      ! finally compute the xc potential in the AO basis
      CALL cp_fm_scale_and_add(alpha=-1.0_dp, matrix_a=fm_V_xc_Gamma, &
                               beta=1.0_dp, matrix_b=gw_env%fm_ks_Gamma(1))

      ! set back the energy
      energy%total = energy_total
      energy%exc = energy_exc
      energy%ex = energy_ex

      ! set the DFT functional and HF fraction back
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)

      CALL dbcsr_deallocate_matrix_set(mat_ks_without_v_xc)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_V_xc_Gamma

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
!> \param Sigma_c_ikp_n_freq ...
!> \param Sigma_x_ikp_n ...
!> \param V_xc_ikp_n ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE analyt_conti_and_print(gw_env, Sigma_c_ikp_n_freq, Sigma_x_ikp_n, V_xc_ikp_n, ikp)

      TYPE(gw_type), POINTER                             :: gw_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_freq
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n, V_xc_ikp_n
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'analyt_conti_and_print'

      CHARACTER(len=3)                                   :: occ_vir
      CHARACTER(len=default_path_length)                 :: fname
      INTEGER                                            :: handle, i_mo, iunit, n_mo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dummy, Sigma_c_ikp_n_qp

      CALL timeset(routineN, handle)

      n_mo = gw_env%n_ao
      ALLOCATE (dummy(n_mo), Sigma_c_ikp_n_qp(n_mo))
      Sigma_c_ikp_n_qp(:) = 0.0_dp

      DO i_mo = 1, n_mo

         IF (MODULO(i_mo, gw_env%para_env%num_pe) /= gw_env%para_env%mepos) CYCLE

         CALL continuation_pade(Sigma_c_ikp_n_qp, &
                                gw_env%imag_freq_points_fit, dummy, dummy, &
                                Sigma_c_ikp_n_freq(:, 1:gw_env%num_freq_points_fit, 1)*z_one + &
                                Sigma_c_ikp_n_freq(:, 1:gw_env%num_freq_points_fit, 2)*gaussi, &
                                Sigma_x_ikp_n(:) - V_xc_ikp_n(:), &
                                gw_env%eigenvalues_scf(:, ikp), &
                                gw_env%eigenvalues_scf(:, ikp), &
                                i_mo, gw_env%n_occ(1), gw_env%nparam_pade, &
                                gw_env%num_freq_points_fit, &
                                ri_rpa_g0w0_crossing_newton, gw_env%n_occ(1), &
                                0.0_dp, .TRUE., .FALSE., 1)

      END DO

      CALL gw_env%para_env%sum(Sigma_c_ikp_n_qp)

      gw_env%eigenvalues_G0W0(:, ikp) = gw_env%eigenvalues_scf(:, ikp) + &
                                        Sigma_c_ikp_n_qp(:) + &
                                        Sigma_x_ikp_n(:) - &
                                        V_xc_ikp_n(:)

      IF (ikp < 10) THEN
         WRITE (fname, "(A42,I1)") "DFT_and_G0W0_band_structure_for_kpoint_000", ikp
      ELSE IF (ikp < 100) THEN
         WRITE (fname, "(A41,I2)") "DFT_and_G0W0_band_structure_for_kpoint_00", ikp
      ELSE IF (ikp < 1000) THEN
         WRITE (fname, "(A40,I3)") "DFT_and_G0W0_band_structure_for_kpoint_0", ikp
      ELSE
         WRITE (fname, "(A39,I4)") "DFT_and_G0W0_band_structure_for_kpoint_", ikp
      END IF

      CALL open_file(TRIM(fname), unit_number=iunit, file_status="UNKNOWN", file_action="WRITE")

      WRITE (iunit, "(A)") " "
      WRITE (iunit, "(A10,3F10.4)") "kpoint: ", gw_env%kpoints_Sigma%xkp(:, ikp)
      WRITE (iunit, "(A)") " "
      WRITE (iunit, "(A5,A24,2A17,A16,A18)") "n", "ϵ_nk^DFT (eV)", "Σ^c_nk (eV)", &
         "Σ^x_nk (eV)", "v_n^xc (eV)", "ϵ_nk^G0W0 (eV)"
      WRITE (iunit, "(A)") " "
      DO i_mo = 1, n_mo
         IF (i_mo .LE. gw_env%n_occ(1)) occ_vir = 'occ'
         IF (i_mo > gw_env%n_occ(1)) occ_vir = 'vir'
         WRITE (iunit, "(I5,3A,4F16.3,F17.3)") i_mo, ' (', occ_vir, ') ', &
            gw_env%eigenvalues_scf(i_mo, ikp)*evolt, &
            Sigma_c_ikp_n_qp(i_mo)*evolt, &
            Sigma_x_ikp_n(i_mo)*evolt, &
            V_xc_ikp_n(i_mo)*evolt, &
            gw_env%eigenvalues_G0W0(i_mo, ikp)*evolt
      END DO

      CALL close_file(iunit)

      CALL timestop(handle)

   END SUBROUTINE analyt_conti_and_print

! **************************************************************************************************
!> \brief ...
!> \param gw_env ...
! **************************************************************************************************
   SUBROUTINE get_VBM_CBM_bandgaps_and_print(gw_env)

      TYPE(gw_type), POINTER                             :: gw_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_VBM_CBM_bandgaps_and_print'

      INTEGER                                            :: handle, homo, ikp, lumo, n_mo, u
      REAL(KIND=dp)                                      :: E_DBG_G0W0_at_ikp, E_DBG_scf_at_ikp

      CALL timeset(routineN, handle)

      homo = gw_env%n_occ(1)
      lumo = homo + 1
      n_mo = gw_env%n_ao
      gw_env%E_VBM_scf = MAXVAL(gw_env%eigenvalues_scf(1:homo, :))
      gw_env%E_CBM_scf = MINVAL(gw_env%eigenvalues_scf(homo + 1:n_mo, :))
      gw_env%E_IDBG_scf = gw_env%E_CBM_scf - gw_env%E_VBM_scf
      gw_env%E_VBM_G0W0 = MAXVAL(gw_env%eigenvalues_G0W0(1:homo, :))
      gw_env%E_CBM_G0W0 = MINVAL(gw_env%eigenvalues_G0W0(homo + 1:n_mo, :))
      gw_env%E_IDBG_G0W0 = gw_env%E_CBM_G0W0 - gw_env%E_VBM_G0W0

      gw_env%E_DBG_scf = 1000.0_dp
      gw_env%E_DBG_G0W0 = 1000.0_dp
      DO ikp = 1, gw_env%kpoints_Sigma%nkp
         E_DBG_scf_at_ikp = -MAXVAL(gw_env%eigenvalues_scf(1:homo, ikp)) + &
                            MINVAL(gw_env%eigenvalues_scf(homo + 1:n_mo, ikp))
         IF (E_DBG_scf_at_ikp < gw_env%E_DBG_scf) gw_env%E_DBG_scf = E_DBG_scf_at_ikp

         E_DBG_G0W0_at_ikp = -MAXVAL(gw_env%eigenvalues_G0W0(1:homo, ikp)) + &
                             MINVAL(gw_env%eigenvalues_G0W0(homo + 1:n_mo, ikp))
         IF (E_DBG_G0W0_at_ikp < gw_env%E_DBG_G0W0) gw_env%E_DBG_G0W0 = E_DBG_G0W0_at_ikp
      END DO

      u = gw_env%unit_nr
      IF (u > 0) THEN
         WRITE (u, '(T2, 2A)') 'Quasiparticle energies written to ', &
            'DFT_and_G0W0_band_structure_for_kpoint_<ikp>'
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,F49.3)') 'DFT valence band maximum (eV):', gw_env%E_VBM_scf*evolt
         WRITE (u, '(T2,A,F46.3)') 'DFT conduction band minimum (eV):', gw_env%E_CBM_scf*evolt
         WRITE (u, '(T2,A,F52.3)') 'DFT indirect band gap (eV):', gw_env%E_IDBG_scf*evolt
         WRITE (u, '(T2,A,F54.3)') 'DFT direct band gap (eV):', gw_env%E_DBG_scf*evolt
         WRITE (u, '(T2,A)') ''
         WRITE (u, '(T2,A,F48.3)') 'G0W0 valence band maximum (eV):', gw_env%E_VBM_G0W0*evolt
         WRITE (u, '(T2,A,F45.3)') 'G0W0 conduction band minimum (eV):', gw_env%E_CBM_G0W0*evolt
         WRITE (u, '(T2,A,F51.3)') 'G0W0 indirect band gap (eV):', gw_env%E_IDBG_G0W0*evolt
         WRITE (u, '(T2,A,F53.3)') 'G0W0 direct band gap (eV):', gw_env%E_DBG_G0W0*evolt
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_VBM_CBM_bandgaps_and_print

END MODULE gw_methods
