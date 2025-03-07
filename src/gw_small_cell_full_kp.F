!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 05.2024
! **************************************************************************************************
MODULE gw_small_cell_full_kp
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_diag,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE dbt_api,                         ONLY: dbt_clear,&
                                              dbt_contract,&
                                              dbt_copy,&
                                              dbt_create,&
                                              dbt_destroy,&
                                              dbt_type
   USE gw_communication,                ONLY: fm_to_local_array,&
                                              fm_to_local_tensor,&
                                              local_array_to_fm,&
                                              local_dbt_to_global_fm
   USE gw_kp_to_real_space_and_back,    ONLY: add_ikp_to_all_rs,&
                                              fm_add_ikp_to_rs,&
                                              fm_trafo_rs_to_ikp,&
                                              trafo_rs_to_ikp
   USE gw_utils,                        ONLY: add_R,&
                                              analyt_conti_and_print,&
                                              de_init_bs_env,&
                                              get_V_tr_R,&
                                              is_cell_in_index_to_cell,&
                                              power,&
                                              time_to_freq
   USE kinds,                           ONLY: dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp_small_cell
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: get_all_VBM_CBM_bandgaps
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_small_cell_full_kp'

   PUBLIC :: gw_calc_small_cell_full_kp

CONTAINS

! **************************************************************************************************
!> \brief Perform GW band structure calculation
!> \param qs_env ...
!> \param bs_env ...
!> \par History
!>    * 05.2024 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE gw_calc_small_cell_full_kp(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_small_cell_full_kp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! G^occ_µλ(i|τ|,k) = sum_n^occ C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
      ! G^vir_µλ(i|τ|,k) = sum_n^vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
      ! k-point k -> cell S: G^occ/vir_µλ^S(i|τ|) = sum_k w_k G^occ/vir_µλ(i|τ|,k) e^(ikS)
      ! χ_PQ^R(iτ) = sum_λR1νR2 [ sum_µS (µR1-S νR2 | P0) G^vir_λµ^S(i|τ|) ]
      !                         [ sum_σS (σR2-S λR1 | QR) G^occ_νσ^S(i|τ|) ]
      CALL compute_chi(bs_env)

      ! χ_PQ^R(iτ) -> χ_PQ(iω,k) -> ε_PQ(iω,k) -> W_PQ(iω,k) -> Ŵ(iω,k) = M^-1(k)*W(iω,k)*M^-1(k)
      !            -> Ŵ_PQ^R(iτ)
      CALL compute_W_real_space(bs_env, qs_env)

      ! D_µν(k) = sum_n^occ C^*_µn(k) C_νn(k), V^tr_PQ^R = <phi_P,0|V^tr|phi_Q,R>
      ! V^tr(k) = sum_R e^ikR V^tr^R, M(k) = sum_R e^ikR M^R, M(k) -> M^-1(k)
      ! -> Ṽ^tr(k) = M^-1(k) * V^tr(k) * M^-1(k) -> Ṽ^tr_PQ^R = sum_k w_k e^-ikR Ṽ^tr_PQ(k)
      ! Σ^x_λσ^R = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) D_µν^S2    ]
      !                       [ sum_QR2 (σR νS1    | QR1-R2) Ṽ^tr_PQ^R2 ]
      CALL compute_Sigma_x(bs_env, qs_env)

      ! Σ^c_λσ^R(iτ) = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) G^occ/vir_µν^S2(i|τ|) ]
      !                           [ sum_QR2 (σR νS1    | QR1-R2) Ŵ_PQ^R2(iτ)           ]
      CALL compute_Sigma_c(bs_env)

      ! Σ^c_λσ^R(iτ,k=0) -> Σ^c_nn(ϵ,k); ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ,k) + Σ^x_nn(k) - v^xc_nn(k)
      CALL compute_QP_energies(bs_env)

      CALL de_init_bs_env(bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_small_cell_full_kp

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_chi(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_chi'

      INTEGER                                            :: cell_DR(3), cell_R1(3), cell_R2(3), &
                                                            handle, i_cell_Delta_R, i_cell_R1, &
                                                            i_cell_R2, i_t, i_task_Delta_R_local, &
                                                            ispin
      LOGICAL                                            :: cell_found
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: Gocc_S, Gvir_S, t_chi_R
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_Gocc, t_Gvir

      CALL timeset(routineN, handle)

      DO i_t = 1, bs_env%num_time_freq_points

         CALL dbt_create_2c_R(Gocc_S, bs_env%t_G, bs_env%nimages_scf_desymm)
         CALL dbt_create_2c_R(Gvir_S, bs_env%t_G, bs_env%nimages_scf_desymm)
         CALL dbt_create_2c_R(t_chi_R, bs_env%t_chi, bs_env%nimages_scf_desymm)
         CALL dbt_create_3c_R1_R2(t_Gocc, bs_env%t_RI_AO__AO, bs_env%nimages_3c, bs_env%nimages_3c)
         CALL dbt_create_3c_R1_R2(t_Gvir, bs_env%t_RI_AO__AO, bs_env%nimages_3c, bs_env%nimages_3c)

         t1 = m_walltime()
         tau = bs_env%imag_time_points(i_t)

         DO ispin = 1, bs_env%n_spin

            ! 1. compute G^occ,S(iτ) and G^vir^S(iτ) in imaginary time for cell S
            !    Background: G^σ,S(iτ) = G^occ,S,σ(iτ) * Θ(-τ) + G^vir,S,σ(iτ) * Θ(τ), σ ∈ {↑,↓}
            !    G^occ_µλ(i|τ|,k) = sum_n^occ C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
            !    G^vir_µλ(i|τ|,k) = sum_n^vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
            !    k-point k -> cell S: G^occ/vir_µλ^S(i|τ|) = sum_k w_k G^occ/vir_µλ(i|τ|,k) e^(ikS)
            CALL G_occ_vir(bs_env, tau, Gocc_S, ispin, occ=.TRUE., vir=.FALSE.)
            CALL G_occ_vir(bs_env, tau, Gvir_S, ispin, occ=.FALSE., vir=.TRUE.)

            ! loop over ΔR = R_1 - R_2 which are local in the tensor subgroup
            DO i_task_Delta_R_local = 1, bs_env%n_tasks_Delta_R_local

               i_cell_Delta_R = bs_env%task_Delta_R(i_task_Delta_R_local)

               DO i_cell_R2 = 1, bs_env%nimages_3c

                  cell_R2(1:3) = bs_env%index_to_cell_3c(i_cell_R2, 1:3)
                  cell_DR(1:3) = bs_env%index_to_cell_Delta_R(i_cell_Delta_R, 1:3)

                  ! R_1 = R_2 + ΔR (from ΔR = R_2 - R_1)
                  CALL add_R(cell_R2, cell_DR, bs_env%index_to_cell_3c, cell_R1, &
                             cell_found, bs_env%cell_to_index_3c, i_cell_R1)

                  ! 3-cells check because in M^vir_νR2,λR1,QR (step 3.): R2 is index on ν
                  IF (.NOT. cell_found) CYCLE
                  ! 2. M^occ/vir_λR1,νR2,P0 = sum_µS (λR1 µR2-S | P0) G^occ/vir_νµ^S(iτ)
                  CALL G_times_3c(Gocc_S, t_Gocc, bs_env, i_cell_R1, i_cell_R2)
                  CALL G_times_3c(Gvir_S, t_Gvir, bs_env, i_cell_R2, i_cell_R1)

               END DO ! i_cell_R2

               ! 3. χ_PQ^R(iτ) = sum_λR1,νR2 M^occ_λR1,νR2,P0 M^vir_νR2,λR1,QR
               CALL contract_M_occ_vir_to_chi(t_Gocc, t_Gvir, t_chi_R, &
                                              bs_env, i_cell_Delta_R)

            END DO ! i_cell_Delta_R_local

         END DO ! ispin

         CALL bs_env%para_env%sync()

         CALL local_dbt_to_global_fm(t_chi_R, bs_env%fm_chi_R_t(:, i_t), bs_env%mat_RI_RI, &
                                     bs_env%mat_RI_RI_tensor, bs_env)

         CALL destroy_t_1d(Gocc_S)
         CALL destroy_t_1d(Gvir_S)
         CALL destroy_t_1d(t_chi_R)
         CALL destroy_t_2d(t_Gocc)
         CALL destroy_t_2d(t_Gvir)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ^R(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',      Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! i_t

      CALL timestop(handle)

   END SUBROUTINE compute_chi

! *************************************************************************************************
!> \brief ...
!> \param R ...
!> \param template ...
!> \param nimages ...
! **************************************************************************************************
   SUBROUTINE dbt_create_2c_R(R, template, nimages)

      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: R
      TYPE(dbt_type)                                     :: template
      INTEGER                                            :: nimages

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dbt_create_2c_R'

      INTEGER                                            :: handle, i_cell_S

      CALL timeset(routineN, handle)

      ALLOCATE (R(nimages))
      DO i_cell_S = 1, nimages
         CALL dbt_create(template, R(i_cell_S))
      END DO

      CALL timestop(handle)

   END SUBROUTINE dbt_create_2c_R

! **************************************************************************************************
!> \brief ...
!> \param t_3c_R1_R2 ...
!> \param t_3c_template ...
!> \param nimages_1 ...
!> \param nimages_2 ...
! **************************************************************************************************
   SUBROUTINE dbt_create_3c_R1_R2(t_3c_R1_R2, t_3c_template, nimages_1, nimages_2)

      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_3c_R1_R2
      TYPE(dbt_type)                                     :: t_3c_template
      INTEGER                                            :: nimages_1, nimages_2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbt_create_3c_R1_R2'

      INTEGER                                            :: handle, i_cell, j_cell

      CALL timeset(routineN, handle)

      ALLOCATE (t_3c_R1_R2(nimages_1, nimages_2))
      DO i_cell = 1, nimages_1
         DO j_cell = 1, nimages_2
            CALL dbt_create(t_3c_template, t_3c_R1_R2(i_cell, j_cell))
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE dbt_create_3c_R1_R2

! **************************************************************************************************
!> \brief ...
!> \param t_G_S ...
!> \param t_M ...
!> \param bs_env ...
!> \param i_cell_R1 ...
!> \param i_cell_R2 ...
! **************************************************************************************************
   SUBROUTINE G_times_3c(t_G_S, t_M, bs_env, i_cell_R1, i_cell_R2)
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_G_S
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_M
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: i_cell_R1, i_cell_R2

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_times_3c'

      INTEGER                                            :: handle, i_cell_R1_p_S, i_cell_S
      INTEGER, DIMENSION(3)                              :: cell_R1, cell_R1_plus_cell_S, cell_R2, &
                                                            cell_S
      LOGICAL                                            :: cell_found
      TYPE(dbt_type)                                     :: t_3c_int

      CALL timeset(routineN, handle)

      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_int)

      cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
      cell_R2(1:3) = bs_env%index_to_cell_3c(i_cell_R2, 1:3)

      DO i_cell_S = 1, bs_env%nimages_scf_desymm

         cell_S(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_S, 1:3)
         cell_R1_plus_cell_S(1:3) = cell_R1(1:3) + cell_S(1:3)

         CALL is_cell_in_index_to_cell(cell_R1_plus_cell_S, bs_env%index_to_cell_3c, cell_found)

         IF (.NOT. cell_found) CYCLE

         i_cell_R1_p_S = bs_env%cell_to_index_3c(cell_R1_plus_cell_S(1), cell_R1_plus_cell_S(2), &
                                                 cell_R1_plus_cell_S(3))

         IF (bs_env%nblocks_3c(i_cell_R2, i_cell_R1_p_S) == 0) CYCLE

         CALL get_t_3c_int(t_3c_int, bs_env, i_cell_R2, i_cell_R1_p_S)

         CALL dbt_contract(alpha=1.0_dp, &
                           tensor_1=t_3c_int, &
                           tensor_2=t_G_S(i_cell_S), &
                           beta=1.0_dp, &
                           tensor_3=t_M(i_cell_R1, i_cell_R2), &
                           contract_1=[3], notcontract_1=[1, 2], map_1=[1, 2], &
                           contract_2=[2], notcontract_2=[1], map_2=[3], &
                           filter_eps=bs_env%eps_filter)
      END DO

      CALL dbt_destroy(t_3c_int)

      CALL timestop(handle)

   END SUBROUTINE G_times_3c

! **************************************************************************************************
!> \brief ...
!> \param t_3c_int ...
!> \param bs_env ...
!> \param j_cell ...
!> \param k_cell ...
! **************************************************************************************************
   SUBROUTINE get_t_3c_int(t_3c_int, bs_env, j_cell, k_cell)

      TYPE(dbt_type)                                     :: t_3c_int
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: j_cell, k_cell

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_t_3c_int'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbt_clear(t_3c_int)
      IF (j_cell < k_cell) THEN
         CALL dbt_copy(bs_env%t_3c_int(k_cell, j_cell), t_3c_int, order=[1, 3, 2])
      ELSE
         CALL dbt_copy(bs_env%t_3c_int(j_cell, k_cell), t_3c_int)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_t_3c_int

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param tau ...
!> \param G_S ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
! **************************************************************************************************
   SUBROUTINE G_occ_vir(bs_env, tau, G_S, ispin, occ, vir)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp)                                      :: tau
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: G_S
      INTEGER                                            :: ispin
      LOGICAL                                            :: occ, vir

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'G_occ_vir'

      INTEGER                                            :: handle, homo, i_cell_S, ikp, j, &
                                                            j_col_local, n_mo, ncol_local, &
                                                            nimages, nkp
      INTEGER, DIMENSION(:), POINTER                     :: col_indices
      REAL(KIND=dp)                                      :: tau_E

      CALL timeset(routineN, handle)

      CPASSERT(occ .NEQV. vir)

      CALL cp_cfm_get_info(matrix=bs_env%cfm_work_mo, &
                           ncol_local=ncol_local, &
                           col_indices=col_indices)

      nkp = bs_env%nkp_scf_desymm
      nimages = bs_env%nimages_scf_desymm
      n_mo = bs_env%n_ao
      homo = bs_env%n_occ(ispin)

      DO i_cell_S = 1, bs_env%nimages_scf_desymm
         CALL cp_fm_set_all(bs_env%fm_G_S(i_cell_S), 0.0_dp)
      END DO

      DO ikp = 1, nkp

         ! get C_µn(k)
         CALL cp_cfm_to_cfm(bs_env%cfm_mo_coeff_kp(ikp, ispin), bs_env%cfm_work_mo)

         ! G^occ/vir_µλ(i|τ|,k) = sum_n^occ/vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k)
         DO j_col_local = 1, ncol_local

            j = col_indices(j_col_local)

            ! 0.5 * |(ϵ_nk-ϵ_F)τ|
            tau_E = ABS(tau*0.5_dp*(bs_env%eigenval_scf(j, ikp, ispin) - bs_env%e_fermi(ispin)))

            IF (tau_E < bs_env%stabilize_exp) THEN
               bs_env%cfm_work_mo%local_data(:, j_col_local) = &
                  bs_env%cfm_work_mo%local_data(:, j_col_local)*EXP(-tau_E)
            ELSE
               bs_env%cfm_work_mo%local_data(:, j_col_local) = z_zero
            END IF

            IF ((occ .AND. j > homo) .OR. (vir .AND. j <= homo)) THEN
               bs_env%cfm_work_mo%local_data(:, j_col_local) = z_zero
            END IF

         END DO

         CALL parallel_gemm(transa="N", transb="C", m=n_mo, n=n_mo, k=n_mo, alpha=z_one, &
                            matrix_a=bs_env%cfm_work_mo, matrix_b=bs_env%cfm_work_mo, &
                            beta=z_zero, matrix_c=bs_env%cfm_work_mo_2)

         ! trafo k-point k -> cell S:  G^occ/vir_µλ(i|τ|,k) -> G^occ/vir,S_µλ(i|τ|)
         CALL fm_add_ikp_to_rs(bs_env%cfm_work_mo_2, bs_env%fm_G_S, &
                               bs_env%kpoints_scf_desymm, ikp)

      END DO ! ikp

      ! replicate to tensor from local tensor group
      DO i_cell_S = 1, bs_env%nimages_scf_desymm
         CALL fm_to_local_tensor(bs_env%fm_G_S(i_cell_S), bs_env%mat_ao_ao%matrix, &
                                 bs_env%mat_ao_ao_tensor%matrix, G_S(i_cell_S), bs_env)
      END DO

      CALL timestop(handle)

   END SUBROUTINE G_occ_vir

! **************************************************************************************************
!> \brief ...
!> \param t_Gocc ...
!> \param t_Gvir ...
!> \param t_chi_R ...
!> \param bs_env ...
!> \param i_cell_Delta_R ...
! **************************************************************************************************
   SUBROUTINE contract_M_occ_vir_to_chi(t_Gocc, t_Gvir, t_chi_R, bs_env, i_cell_Delta_R)
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_Gocc, t_Gvir
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_chi_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: i_cell_Delta_R

      CHARACTER(LEN=*), PARAMETER :: routineN = 'contract_M_occ_vir_to_chi'

      INTEGER                                            :: handle, i_cell_R, i_cell_R1, &
                                                            i_cell_R1_minus_R, i_cell_R2, &
                                                            i_cell_R2_minus_R
      INTEGER, DIMENSION(3)                              :: cell_DR, cell_R, cell_R1, &
                                                            cell_R1_minus_R, cell_R2, &
                                                            cell_R2_minus_R
      LOGICAL                                            :: cell_found
      TYPE(dbt_type)                                     :: t_Gocc_2, t_Gvir_2

      CALL timeset(routineN, handle)

      CALL dbt_create(bs_env%t_RI__AO_AO, t_Gocc_2)
      CALL dbt_create(bs_env%t_RI__AO_AO, t_Gvir_2)

      ! χ_PQ^R(iτ) = sum_λR1,νR2 M^occ_λR1,νR2,P0 M^vir_νR2,λR1,QR
      DO i_cell_R = 1, bs_env%nimages_scf_desymm

         DO i_cell_R2 = 1, bs_env%nimages_3c

            cell_R(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)
            cell_R2(1:3) = bs_env%index_to_cell_3c(i_cell_R2, 1:3)
            cell_DR(1:3) = bs_env%index_to_cell_Delta_R(i_cell_Delta_R, 1:3)

            ! R_1 = R_2 + ΔR (from ΔR = R_2 - R_1)
            CALL add_R(cell_R2, cell_DR, bs_env%index_to_cell_3c, cell_R1, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_R1)
            IF (.NOT. cell_found) CYCLE

            ! R_1 - R
            CALL add_R(cell_R1, -cell_R, bs_env%index_to_cell_3c, cell_R1_minus_R, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_R1_minus_R)
            IF (.NOT. cell_found) CYCLE

            ! R_2 - R
            CALL add_R(cell_R2, -cell_R, bs_env%index_to_cell_3c, cell_R2_minus_R, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_R2_minus_R)
            IF (.NOT. cell_found) CYCLE

            ! reorder tensors for efficient contraction to χ_PQ^R
            CALL dbt_copy(t_Gocc(i_cell_R1, i_cell_R2), t_Gocc_2, order=[1, 3, 2])
            CALL dbt_copy(t_Gvir(i_cell_R2_minus_R, i_cell_R1_minus_R), t_Gvir_2)

            ! χ_PQ^R(iτ) = sum_λR1,νR2 M^occ_λR1,νR2,P0 M^vir_νR2,λR1,QR
            CALL dbt_contract(alpha=bs_env%spin_degeneracy, &
                              tensor_1=t_Gocc_2, tensor_2=t_Gvir_2, &
                              beta=1.0_dp, tensor_3=t_chi_R(i_cell_R), &
                              contract_1=[2, 3], notcontract_1=[1], map_1=[1], &
                              contract_2=[2, 3], notcontract_2=[1], map_2=[2], &
                              filter_eps=bs_env%eps_filter, move_data=.TRUE.)
         END DO ! i_cell_R2

      END DO ! i_cell_R

      ! remove all data from t_Gocc and t_Gvir to safe memory
      DO i_cell_R1 = 1, bs_env%nimages_3c
         DO i_cell_R2 = 1, bs_env%nimages_3c
            CALL dbt_clear(t_Gocc(i_cell_R1, i_cell_R2))
            CALL dbt_clear(t_Gvir(i_cell_R1, i_cell_R2))
         END DO
      END DO

      CALL dbt_destroy(t_Gocc_2)
      CALL dbt_destroy(t_Gvir_2)

      CALL timestop(handle)

   END SUBROUTINE contract_M_occ_vir_to_chi

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_W_real_space(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_W_real_space'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: chi_k_w, eps_k_w, W_k_w, work
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: M_inv, M_inv_V_sqrt, V_sqrt
      INTEGER                                            :: handle, i_t, ikp, ikp_local, j_w, n_RI, &
                                                            nimages_scf_desymm
      REAL(KIND=dp)                                      :: freq_j, t1, time_i, weight_ij
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: chi_R, MWM_R, W_R

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      nimages_scf_desymm = bs_env%nimages_scf_desymm

      ALLOCATE (chi_k_w(n_RI, n_RI), work(n_RI, n_RI), eps_k_w(n_RI, n_RI), W_k_w(n_RI, n_RI))
      ALLOCATE (chi_R(n_RI, n_RI, nimages_scf_desymm), W_R(n_RI, n_RI, nimages_scf_desymm), &
                MWM_R(n_RI, n_RI, nimages_scf_desymm))

      t1 = m_walltime()

      CALL compute_Minv_and_Vsqrt(bs_env, qs_env, M_inv_V_sqrt, M_inv, V_sqrt)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed V_PQ(k),', 'Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      t1 = m_walltime()

      DO j_w = 1, bs_env%num_time_freq_points

         ! χ_PQ^R(iτ) -> χ_PQ^R(iω_j) (which is stored in chi_R, single ω_j from j_w loop)
         chi_R(:, :, :) = 0.0_dp
         DO i_t = 1, bs_env%num_time_freq_points
            freq_j = bs_env%imag_freq_points(j_w)
            time_i = bs_env%imag_time_points(i_t)
            weight_ij = bs_env%weights_cos_t_to_w(j_w, i_t)*COS(time_i*freq_j)

            CALL fm_to_local_array(bs_env%fm_chi_R_t(:, i_t), chi_R, weight_ij, add=.TRUE.)
         END DO

         ikp_local = 0
         W_R(:, :, :) = 0.0_dp
         DO ikp = 1, bs_env%nkp_chi_eps_W_orig_plus_extra

            ! trivial parallelization over k-points
            IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

            ikp_local = ikp_local + 1

            ! 1. χ_PQ^R(iω_j) -> χ_PQ(iω_j,k)
            CALL trafo_rs_to_ikp(chi_R, chi_k_w, bs_env%kpoints_scf_desymm%index_to_cell, &
                                 bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))

            ! 2. remove negative eigenvalues from χ_PQ(iω,k)
            CALL power(chi_k_w, 1.0_dp, bs_env%eps_eigval_mat_RI)

            ! 3. ε(iω_j,k_i) = Id - V^0.5(k_i)*M^-1(k_i)*χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)

            ! 3. a) work = χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)
            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, chi_k_w, n_RI, &
                       M_inv_V_sqrt(:, :, ikp_local), n_RI, z_zero, work, n_RI)

            ! 3. b) eps_work = V^0.5(k_i)*M^-1(k_i)*work
            CALL ZGEMM('C', 'N', n_RI, n_RI, n_RI, z_one, M_inv_V_sqrt(:, :, ikp_local), n_RI, &
                       work, n_RI, z_zero, eps_k_w, n_RI)

            ! 3. c) ε(iω_j,k_i) = eps_work - Id
            CALL add_on_diag(eps_k_w, z_one)

            ! 4. W(iω_j,k_i) = M^-1(k_i)*V^0.5(k_i)*(ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)*M^-1(k_i)

            ! 4. a) Inversion of ε(iω_j,k_i) using its Cholesky decomposition
            CALL power(eps_k_w, -1.0_dp, 0.0_dp)

            ! 4. b) ε^-1(iω_j,k_i)-Id
            CALL add_on_diag(eps_k_w, -z_one)

            ! 4. c) work = (ε^-1(iω_j,k_i)-Id)*V^0.5(k_i)
            CALL ZGEMM('N', 'C', n_RI, n_RI, n_RI, z_one, eps_k_w, n_RI, &
                       V_sqrt(:, :, ikp_local), n_RI, z_zero, work, n_RI)

            ! 4. d) W(iω,k_i) = V^0.5(k_i)*work
            CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, V_sqrt(:, :, ikp_local), n_RI, &
                       work, n_RI, z_zero, W_k_w, n_RI)

            ! 5. W(iω,k_i) -> W^R(iω) = sum_k w_k e^(-ikR) W(iω,k) (k-point extrapolation here)
            CALL add_ikp_to_all_rs(W_k_w, W_R, bs_env%kpoints_chi_eps_W, ikp, &
                                   index_to_cell_ext=bs_env%kpoints_scf_desymm%index_to_cell)

         END DO ! ikp

         CALL bs_env%para_env%sync()
         CALL bs_env%para_env%sum(W_R)

         ! 6. W^R(iω) -> W(iω,k) [k-mesh is not extrapolated for stable mult. with M^-1(k) ]
         !            -> M^-1(k)*W(iω,k)*M^-1(k) =: Ŵ(iω,k) -> Ŵ^R(iω) (stored in MWM_R)
         CALL mult_W_with_Minv(W_R, MWM_R, bs_env, qs_env)

         ! 7. Ŵ^R(iω) -> Ŵ^R(iτ) and to fully distributed fm matrix bs_env%fm_MWM_R_t
         DO i_t = 1, bs_env%num_time_freq_points
            freq_j = bs_env%imag_freq_points(j_w)
            time_i = bs_env%imag_time_points(i_t)
            weight_ij = bs_env%weights_cos_w_to_t(i_t, j_w)*COS(time_i*freq_j)
            CALL local_array_to_fm(MWM_R, bs_env%fm_MWM_R_t(:, i_t), weight_ij, add=.TRUE.)
         END DO ! i_t

      END DO ! j_w

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T60,A,F7.1,A)') &
            'Computed W_PQ(k,iω) for all k and τ,', 'Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_W_real_space

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param M_inv_V_sqrt ...
!> \param M_inv ...
!> \param V_sqrt ...
! **************************************************************************************************
   SUBROUTINE compute_Minv_and_Vsqrt(bs_env, qs_env, M_inv_V_sqrt, M_inv, V_sqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)  :: M_inv_V_sqrt, M_inv, V_sqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Minv_and_Vsqrt'

      INTEGER                                            :: handle, ikp, ikp_local, n_RI, nkp, &
                                                            nkp_local, nkp_orig
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: M_R

      CALL timeset(routineN, handle)

      nkp = bs_env%nkp_chi_eps_W_orig_plus_extra
      nkp_orig = bs_env%nkp_chi_eps_W_orig
      n_RI = bs_env%n_RI

      nkp_local = 0
      DO ikp = 1, nkp
         ! trivial parallelization over k-points
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE
         nkp_local = nkp_local + 1
      END DO

      ALLOCATE (M_inv_V_sqrt(n_RI, n_RI, nkp_local), M_inv(n_RI, n_RI, nkp_local), &
                V_sqrt(n_RI, n_RI, nkp_local))

      M_inv_V_sqrt(:, :, :) = z_zero
      M_inv(:, :, :) = z_zero
      V_sqrt(:, :, :) = z_zero

      ! 1. 2c Coulomb integrals for the first "original" k-point grid
      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
      CALL build_2c_coulomb_matrix_kp_small_cell(V_sqrt, qs_env, bs_env%kpoints_chi_eps_W, &
                                                 bs_env%size_lattice_sum_V, basis_type="RI_AUX", &
                                                 ikp_start=1, ikp_end=nkp_orig)

      ! 2. 2c Coulomb integrals for the second "extrapolation" k-point grid
      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_extra
      CALL build_2c_coulomb_matrix_kp_small_cell(V_sqrt, qs_env, bs_env%kpoints_chi_eps_W, &
                                                 bs_env%size_lattice_sum_V, basis_type="RI_AUX", &
                                                 ikp_start=nkp_orig + 1, ikp_end=nkp)

      ! now get M^-1(k) and M^-1(k)*V^0.5(k)

      ! compute M^R_PQ = <phi_P,0|V^tr(rc=3Å)|phi_Q,R> for RI metric
      CALL get_V_tr_R(M_R, bs_env%ri_metric, bs_env%regularization_RI, bs_env, qs_env)

      ikp_local = 0
      DO ikp = 1, nkp

         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

         ikp_local = ikp_local + 1

         ! M(k) = sum_R e^ikR M^R
         CALL trafo_rs_to_ikp(M_R, M_inv(:, :, ikp_local), &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_chi_eps_W%xkp(1:3, ikp))

         ! invert M_PQ(k)
         CALL power(M_inv(:, :, ikp_local), -1.0_dp, 0.0_dp)

         ! V^0.5(k)
         CALL power(V_sqrt(:, :, ikp_local), 0.5_dp, 0.0_dp)

         ! M^-1(k)*V^0.5(k)
         CALL ZGEMM("N", "C", n_RI, n_RI, n_RI, z_one, M_inv(:, :, ikp_local), n_RI, &
                    V_sqrt(:, :, ikp_local), n_RI, z_zero, M_inv_V_sqrt(:, :, ikp_local), n_RI)

      END DO ! ikp

      CALL timestop(handle)

   END SUBROUTINE compute_Minv_and_Vsqrt

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param alpha ...
! **************************************************************************************************
   SUBROUTINE add_on_diag(matrix, alpha)
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: matrix
      COMPLEX(KIND=dp)                                   :: alpha

      CHARACTER(len=*), PARAMETER                        :: routineN = 'add_on_diag'

      INTEGER                                            :: handle, i, n

      CALL timeset(routineN, handle)

      n = SIZE(matrix, 1)
      CPASSERT(n == SIZE(matrix, 2))

      DO i = 1, n
         matrix(i, i) = matrix(i, i) + alpha
      END DO

      CALL timestop(handle)

   END SUBROUTINE add_on_diag

! **************************************************************************************************
!> \brief ...
!> \param W_R ...
!> \param MWM_R ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE mult_W_with_Minv(W_R, MWM_R, bs_env, qs_env)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: W_R, MWM_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mult_W_with_Minv'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: M_inv, W_k, work
      INTEGER                                            :: handle, ikp, n_RI
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: M_R

      CALL timeset(routineN, handle)

      ! compute M^R again
      CALL get_V_tr_R(M_R, bs_env%ri_metric, bs_env%regularization_RI, bs_env, qs_env)

      n_RI = bs_env%n_RI
      ALLOCATE (M_inv(n_RI, n_RI), W_k(n_RI, n_RI), work(n_RI, n_RI))
      MWM_R(:, :, :) = 0.0_dp

      DO ikp = 1, bs_env%nkp_scf_desymm

         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE

         ! M(k) = sum_R e^ikR M^R
         CALL trafo_rs_to_ikp(M_R, M_inv, &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))

         ! invert M_PQ(k)
         CALL power(M_inv, -1.0_dp, 0.0_dp)

         ! W(k) = sum_R e^ikR W^R [only R in the supercell that is determined by the SCF k-mesh]
         CALL trafo_rs_to_ikp(W_R, W_k, &
                              bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))

         ! 2e. M^-1(k) W^trunc(k)
         CALL ZGEMM("N", "N", n_RI, n_RI, n_RI, z_one, M_inv, n_RI, W_k, n_RI, z_zero, work, n_RI)

         ! 2f. Ŵ(k) = M^-1(k) W^trunc(k) M^-1(k)
         CALL ZGEMM("N", "N", n_RI, n_RI, n_RI, z_one, work, n_RI, M_inv, n_RI, z_zero, W_k, n_RI)

         ! 2g. Ŵ^R = sum_k w_k e^(-ikR) Ŵ^(k)
         CALL add_ikp_to_all_rs(W_k, MWM_R, bs_env%kpoints_scf_desymm, ikp)

      END DO ! ikp

      CALL bs_env%para_env%sync()
      CALL bs_env%para_env%sum(MWM_R)

      CALL timestop(handle)

   END SUBROUTINE mult_W_with_Minv

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, qs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, i_cell_Delta_R, &
                                                            i_task_Delta_R_local, ispin
      REAL(KIND=dp)                                      :: t1
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: D_S, Mi_Vtr_Mi_R, Sigma_x_R
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_V

      CALL timeset(routineN, handle)

      CALL dbt_create_2c_R(Mi_Vtr_Mi_R, bs_env%t_W, bs_env%nimages_scf_desymm)
      CALL dbt_create_2c_R(D_S, bs_env%t_G, bs_env%nimages_scf_desymm)
      CALL dbt_create_2c_R(Sigma_x_R, bs_env%t_G, bs_env%nimages_scf_desymm)
      CALL dbt_create_3c_R1_R2(t_V, bs_env%t_RI_AO__AO, bs_env%nimages_3c, bs_env%nimages_3c)

      t1 = m_walltime()

      ! V^tr_PQ^R = <phi_P,0|V^tr|phi_Q,R>, V^tr(k) = sum_R e^ikR V^tr^R
      ! M(k) = sum_R e^ikR M^R, M(k) -> M^-1(k) -> Ṽ^tr(k) = M^-1(k) * V^tr(k) * M^-1(k)
      !                                         -> Ṽ^tr_PQ^R = sum_k w_k e^-ikR Ṽ^tr_PQ(k)
      CALL get_Minv_Vtr_Minv_R(Mi_Vtr_Mi_R, bs_env, qs_env)

      ! Σ^x_λσ^R = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) D_µν^S2       ]
      !                       [ sum_QR2 (σR νS1    | QR1-R2) Ṽ^tr_PQ^R2 ]
      DO ispin = 1, bs_env%n_spin

         ! compute D^S(iτ) for cell S from D_µν(k) = sum_n^occ C^*_µn(k) C_νn(k):
         ! trafo k-point k -> cell S: D_µν^S = sum_k w_k D_µν(k) e^(ikS)
         CALL G_occ_vir(bs_env, 0.0_dp, D_S, ispin, occ=.TRUE., vir=.FALSE.)

         ! loop over ΔR = S_1 - R_1 which are local in the tensor subgroup
         DO i_task_Delta_R_local = 1, bs_env%n_tasks_Delta_R_local

            i_cell_Delta_R = bs_env%task_Delta_R(i_task_Delta_R_local)

            ! M^V_σ0,νS1,PR1 = sum_QR2 ( σ0 νS1 | QR1-R2 ) Ṽ^tr_QP^R2 for i_task_local
            CALL contract_W(t_V, Mi_Vtr_Mi_R, bs_env, i_cell_Delta_R)

            ! M^D_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) D_µν^S2
            ! Σ^x_λσ^R = sum_PR1νS1 M^D_λ0,νS1,PR1 * M^V_σR,νS1,PR1 for i_task_local, where
            !                                        M^V_σR,νS1,PR1 = M^V_σ0,νS1-R,PR1-R
            CALL contract_to_Sigma(Sigma_x_R, t_V, D_S, i_cell_Delta_R, bs_env, &
                                   occ=.TRUE., vir=.FALSE., clear_t_W=.TRUE.)

         END DO ! i_cell_Delta_R_local

         CALL bs_env%para_env%sync()

         CALL local_dbt_to_global_fm(Sigma_x_R, bs_env%fm_Sigma_x_R, bs_env%mat_ao_ao, &
                                     bs_env%mat_ao_ao_tensor, bs_env)

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x,', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      CALL destroy_t_1d(Mi_Vtr_Mi_R)
      CALL destroy_t_1d(D_S)
      CALL destroy_t_1d(Sigma_x_R)
      CALL destroy_t_2d(t_V)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_cell_Delta_R, i_t, &
                                                            i_task_Delta_R_local, ispin
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: Gocc_S, Gvir_S, Sigma_c_R_neg_tau, &
                                                            Sigma_c_R_pos_tau, W_R
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_W

      CALL timeset(routineN, handle)

      CALL dbt_create_2c_R(Gocc_S, bs_env%t_G, bs_env%nimages_scf_desymm)
      CALL dbt_create_2c_R(Gvir_S, bs_env%t_G, bs_env%nimages_scf_desymm)
      CALL dbt_create_2c_R(W_R, bs_env%t_W, bs_env%nimages_scf_desymm)
      CALL dbt_create_3c_R1_R2(t_W, bs_env%t_RI_AO__AO, bs_env%nimages_3c, bs_env%nimages_3c)
      CALL dbt_create_2c_R(Sigma_c_R_neg_tau, bs_env%t_G, bs_env%nimages_scf_desymm)
      CALL dbt_create_2c_R(Sigma_c_R_pos_tau, bs_env%t_G, bs_env%nimages_scf_desymm)

      ! Σ^c_λσ^R(iτ) = sum_PR1νS1 [ sum_µS2 (λ0 µS1-S2 | PR1   ) G^occ/vir_µν^S2(i|τ|) ]
      !                           [ sum_QR2 (σR νS1    | QR1-R2) Ŵ_PQ^R2(iτ)           ]
      DO i_t = 1, bs_env%num_time_freq_points

         DO ispin = 1, bs_env%n_spin

            t1 = m_walltime()

            tau = bs_env%imag_time_points(i_t)

            ! G^occ_µλ(i|τ|,k) = sum_n^occ C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k), τ < 0
            ! G^vir_µλ(i|τ|,k) = sum_n^vir C_µn(k)^* e^(-|(ϵ_nk-ϵ_F)τ|) C_λn(k), τ > 0
            ! k-point k -> cell S: G^occ/vir_µλ^S(i|τ|) = sum_k w_k G^occ/vir_µλ(i|τ|,k) e^(ikS)
            CALL G_occ_vir(bs_env, tau, Gocc_S, ispin, occ=.TRUE., vir=.FALSE.)
            CALL G_occ_vir(bs_env, tau, Gvir_S, ispin, occ=.FALSE., vir=.TRUE.)

            ! write data of W^R_PQ(iτ) to W_R 2-index tensor
            CALL fm_MWM_R_t_to_local_tensor_W_R(bs_env%fm_MWM_R_t(:, i_t), W_R, bs_env)

            ! loop over ΔR = S_1 - R_1 which are local in the tensor subgroup
            DO i_task_Delta_R_local = 1, bs_env%n_tasks_Delta_R_local

               i_cell_Delta_R = bs_env%task_Delta_R(i_task_Delta_R_local)

               ! for i_task_local (i.e. fixed ΔR = S_1 - R_1) and for all τ (W(iτ) = W(-iτ)):
               ! M^W_σ0,νS1,PR1 = sum_QR2 ( σ0 νS1 | QR1-R2 ) W(iτ)_QP^R2
               CALL contract_W(t_W, W_R, bs_env, i_cell_Delta_R)

               ! for τ < 0 and for i_task_local (i.e. fixed ΔR = S_1 - R_1):
               ! M^G_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) G^occ(i|τ|)_µν^S2
               ! Σ^c_λσ^R(iτ) = sum_PR1νS1 M^G_λ0,νS1,PR1 * M^W_σR,νS1,PR1
               !                                      where M^W_σR,νS1,PR1 = M^W_σ0,νS1-R,PR1-R
               CALL contract_to_Sigma(Sigma_c_R_neg_tau, t_W, Gocc_S, i_cell_Delta_R, bs_env, &
                                      occ=.TRUE., vir=.FALSE., clear_t_W=.FALSE.)

               ! for τ > 0: same as for τ < 0, but G^occ -> G^vir
               CALL contract_to_Sigma(Sigma_c_R_pos_tau, t_W, Gvir_S, i_cell_Delta_R, bs_env, &
                                      occ=.FALSE., vir=.TRUE., clear_t_W=.TRUE.)

            END DO ! i_cell_Delta_R_local

            CALL bs_env%para_env%sync()

            CALL local_dbt_to_global_fm(Sigma_c_R_pos_tau, &
                                        bs_env%fm_Sigma_c_R_pos_tau(:, i_t, ispin), &
                                        bs_env%mat_ao_ao, bs_env%mat_ao_ao_tensor, bs_env)

            CALL local_dbt_to_global_fm(Sigma_c_R_neg_tau, &
                                        bs_env%fm_Sigma_c_R_neg_tau(:, i_t, ispin), &
                                        bs_env%mat_ao_ao, bs_env%mat_ao_ao_tensor, bs_env)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I10,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point   ', i_t, ' /', bs_env%num_time_freq_points, &
                  ',      Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

      END DO ! i_t

      CALL destroy_t_1d(Gocc_S)
      CALL destroy_t_1d(Gvir_S)
      CALL destroy_t_1d(W_R)
      CALL destroy_t_1d(Sigma_c_R_neg_tau)
      CALL destroy_t_1d(Sigma_c_R_pos_tau)
      CALL destroy_t_2d(t_W)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief ...
!> \param Mi_Vtr_Mi_R ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_Minv_Vtr_Minv_R(Mi_Vtr_Mi_R, bs_env, qs_env)
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: Mi_Vtr_Mi_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_Minv_Vtr_Minv_R'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: M_inv_V_tr_kp, M_kp, Mi_Vtr_Mi_kp, &
                                                            V_tr_kp
      INTEGER                                            :: handle, i_cell_R, ikp, n_RI, &
                                                            nimages_scf, nkp_scf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: M_R, Mi_Vtr_Mi_R_arr, V_tr_R

      CALL timeset(routineN, handle)

      nimages_scf = bs_env%nimages_scf_desymm
      nkp_scf = bs_env%kpoints_scf_desymm%nkp
      n_RI = bs_env%n_RI

      CALL get_V_tr_R(V_tr_R, bs_env%trunc_coulomb, 0.0_dp, bs_env, qs_env)
      CALL get_V_tr_R(M_R, bs_env%ri_metric, bs_env%regularization_RI, bs_env, qs_env)

      ALLOCATE (V_tr_kp(n_RI, n_RI), M_kp(n_RI, n_RI), M_inv_V_tr_kp(n_RI, n_RI), &
                Mi_Vtr_Mi_kp(n_RI, n_RI), Mi_Vtr_Mi_R_arr(n_RI, n_RI, nimages_scf))
      Mi_Vtr_Mi_R_arr(:, :, :) = 0.0_dp

      DO ikp = 1, nkp_scf
         ! trivial parallelization
         IF (MODULO(ikp, bs_env%para_env%num_pe) .NE. bs_env%para_env%mepos) CYCLE
         ! V_tr(k) = sum_R e^ikR V_tr^R
         CALL trafo_rs_to_ikp(V_tr_R, V_tr_kp, bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))
         ! M(k)    = sum_R e^ikR M^R
         CALL trafo_rs_to_ikp(M_R, M_kp, bs_env%kpoints_scf_desymm%index_to_cell, &
                              bs_env%kpoints_scf_desymm%xkp(1:3, ikp))
         ! M(k) -> M^-1(k)
         CALL power(M_kp, -1.0_dp, 0.0_dp)
         ! M^-1(k) * V_tr(k)
         CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, M_kp, n_RI, &
                    V_tr_kp, n_RI, z_zero, M_inv_V_tr_kp, n_RI)
         ! Ṽ(k) = M^-1(k) * V_tr(k) * M^-1(k)
         CALL ZGEMM('N', 'N', n_RI, n_RI, n_RI, z_one, M_inv_V_tr_kp, n_RI, &
                    M_kp, n_RI, z_zero, Mi_Vtr_Mi_kp, n_RI)
         ! Ṽ^R = sum_k w_k e^-ikR Ṽ(k)
         CALL add_ikp_to_all_rs(Mi_Vtr_Mi_kp, Mi_Vtr_Mi_R_arr, bs_env%kpoints_scf_desymm, ikp)
      END DO
      CALL bs_env%para_env%sync()
      CALL bs_env%para_env%sum(Mi_Vtr_Mi_R_arr)

      ! use bs_env%fm_chi_R_t for temporary storage
      CALL local_array_to_fm(Mi_Vtr_Mi_R_arr, bs_env%fm_chi_R_t(:, 1))

      ! communicate Mi_Vtr_Mi_R to tensor format; full replication in tensor group
      DO i_cell_R = 1, nimages_scf
         CALL fm_to_local_tensor(bs_env%fm_chi_R_t(i_cell_R, 1), bs_env%mat_RI_RI%matrix, &
                                 bs_env%mat_RI_RI_tensor%matrix, Mi_Vtr_Mi_R(i_cell_R), bs_env)
      END DO

      CALL timestop(handle)

   END SUBROUTINE get_Minv_Vtr_Minv_R

! **************************************************************************************************
!> \brief ...
!> \param t_1d ...
! **************************************************************************************************
   SUBROUTINE destroy_t_1d(t_1d)
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_1d

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'destroy_t_1d'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      DO i = 1, SIZE(t_1d)
         CALL dbt_destroy(t_1d(i))
      END DO
      DEALLOCATE (t_1d)

      CALL timestop(handle)

   END SUBROUTINE destroy_t_1d

! **************************************************************************************************
!> \brief ...
!> \param t_2d ...
! **************************************************************************************************
   SUBROUTINE destroy_t_2d(t_2d)
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_2d

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'destroy_t_2d'

      INTEGER                                            :: handle, i, j

      CALL timeset(routineN, handle)

      DO i = 1, SIZE(t_2d, 1)
      DO j = 1, SIZE(t_2d, 2)
         CALL dbt_destroy(t_2d(i, j))
      END DO
      END DO
      DEALLOCATE (t_2d)

      CALL timestop(handle)

   END SUBROUTINE destroy_t_2d

! **************************************************************************************************
!> \brief ...
!> \param t_W ...
!> \param W_R ...
!> \param bs_env ...
!> \param i_cell_Delta_R ...
! **************************************************************************************************
   SUBROUTINE contract_W(t_W, W_R, bs_env, i_cell_Delta_R)
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_W
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: W_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: i_cell_Delta_R

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_W'

      INTEGER                                            :: handle, i_cell_R1, i_cell_R2, &
                                                            i_cell_R2_m_R1, i_cell_S1, &
                                                            i_cell_S1_m_R1_p_R2
      INTEGER, DIMENSION(3)                              :: cell_DR, cell_R1, cell_R2, cell_R2_m_R1, &
                                                            cell_S1, cell_S1_m_R2_p_R1
      LOGICAL                                            :: cell_found
      TYPE(dbt_type)                                     :: t_3c_int, t_W_tmp

      CALL timeset(routineN, handle)

      CALL dbt_create(bs_env%t_RI__AO_AO, t_W_tmp)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_int)

      DO i_cell_R1 = 1, bs_env%nimages_3c

         cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
         cell_DR(1:3) = bs_env%index_to_cell_Delta_R(i_cell_Delta_R, 1:3)

         ! S_1 = R_1 + ΔR (from ΔR = S_1 - R_1)
         CALL add_R(cell_R1, cell_DR, bs_env%index_to_cell_3c, cell_S1, &
                    cell_found, bs_env%cell_to_index_3c, i_cell_S1)
         IF (.NOT. cell_found) CYCLE

         DO i_cell_R2 = 1, bs_env%nimages_scf_desymm

            cell_R2(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R2, 1:3)

            ! R_2 - R_1
            CALL add_R(cell_R2, -cell_R1, bs_env%index_to_cell_3c, cell_R2_m_R1, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_R2_m_R1)
            IF (.NOT. cell_found) CYCLE

            ! S_1 - R_1 + R_2
            CALL add_R(cell_S1, cell_R2_m_R1, bs_env%index_to_cell_3c, cell_S1_m_R2_p_R1, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_S1_m_R1_p_R2)
            IF (.NOT. cell_found) CYCLE

            CALL get_t_3c_int(t_3c_int, bs_env, i_cell_S1_m_R1_p_R2, i_cell_R2_m_R1)

            ! M^W_σ0,νS1,PR1 = sum_QR2 ( σ0     νS1       | QR1-R2 ) W_QP^R2
            !                = sum_QR2 ( σR2-R1 νS1-R1+R2 | Q0     ) W_QP^R2
            ! for ΔR = S_1 - R_1
            CALL dbt_contract(alpha=1.0_dp, &
                              tensor_1=W_R(i_cell_R2), &
                              tensor_2=t_3c_int, &
                              beta=0.0_dp, &
                              tensor_3=t_W_tmp, &
                              contract_1=[1], notcontract_1=[2], map_1=[1], &
                              contract_2=[1], notcontract_2=[2, 3], map_2=[2, 3], &
                              filter_eps=bs_env%eps_filter)

            ! reorder tensor
            CALL dbt_copy(t_W_tmp, t_W(i_cell_S1, i_cell_R1), order=[1, 2, 3], &
                          move_data=.TRUE., summation=.TRUE.)

         END DO ! i_cell_R2

      END DO ! i_cell_R1

      CALL dbt_destroy(t_W_tmp)
      CALL dbt_destroy(t_3c_int)

      CALL timestop(handle)

   END SUBROUTINE contract_W

! **************************************************************************************************
!> \brief ...
!> \param Sigma_R ...
!> \param t_W ...
!> \param G_S ...
!> \param i_cell_Delta_R ...
!> \param bs_env ...
!> \param occ ...
!> \param vir ...
!> \param clear_t_W ...
! **************************************************************************************************
   SUBROUTINE contract_to_Sigma(Sigma_R, t_W, G_S, i_cell_Delta_R, bs_env, occ, vir, clear_t_W)
      TYPE(dbt_type), DIMENSION(:)                       :: Sigma_R
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :)       :: t_W
      TYPE(dbt_type), DIMENSION(:)                       :: G_S
      INTEGER                                            :: i_cell_Delta_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      LOGICAL                                            :: occ, vir, clear_t_W

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'contract_to_Sigma'

      INTEGER :: handle, handle2, i_cell_m_R1, i_cell_R, i_cell_R1, i_cell_R1_minus_R, i_cell_S1, &
         i_cell_S1_minus_R, i_cell_S1_p_S2_m_R1, i_cell_S2
      INTEGER, DIMENSION(3)                              :: cell_DR, cell_m_R1, cell_R, cell_R1, &
                                                            cell_R1_minus_R, cell_S1, &
                                                            cell_S1_minus_R, cell_S1_p_S2_m_R1, &
                                                            cell_S2
      LOGICAL                                            :: cell_found
      REAL(KIND=dp)                                      :: sign_Sigma
      TYPE(dbt_type)                                     :: t_3c_int, t_G, t_G_2

      CALL timeset(routineN, handle)

      CPASSERT(occ .EQV. (.NOT. vir))
      IF (occ) sign_Sigma = -1.0_dp
      IF (vir) sign_Sigma = 1.0_dp

      CALL dbt_create(bs_env%t_RI_AO__AO, t_G)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_G_2)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c_int)

      DO i_cell_R1 = 1, bs_env%nimages_3c

         cell_R1(1:3) = bs_env%index_to_cell_3c(i_cell_R1, 1:3)
         cell_DR(1:3) = bs_env%index_to_cell_Delta_R(i_cell_Delta_R, 1:3)

         ! S_1 = R_1 + ΔR (from ΔR = S_1 - R_1)
         CALL add_R(cell_R1, cell_DR, bs_env%index_to_cell_3c, cell_S1, cell_found, &
                    bs_env%cell_to_index_3c, i_cell_S1)
         IF (.NOT. cell_found) CYCLE

         DO i_cell_S2 = 1, bs_env%nimages_scf_desymm

            cell_S2(1:3) = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_S2, 1:3)
            cell_m_R1(1:3) = -cell_R1(1:3)
            cell_S1_p_S2_m_R1(1:3) = cell_S1(1:3) + cell_S2(1:3) - cell_R1(1:3)

            CALL is_cell_in_index_to_cell(cell_m_R1, bs_env%index_to_cell_3c, cell_found)
            IF (.NOT. cell_found) CYCLE

            CALL is_cell_in_index_to_cell(cell_S1_p_S2_m_R1, bs_env%index_to_cell_3c, cell_found)
            IF (.NOT. cell_found) CYCLE

            i_cell_m_R1 = bs_env%cell_to_index_3c(cell_m_R1(1), cell_m_R1(2), cell_m_R1(3))
            i_cell_S1_p_S2_m_R1 = bs_env%cell_to_index_3c(cell_S1_p_S2_m_R1(1), &
                                                          cell_S1_p_S2_m_R1(2), &
                                                          cell_S1_p_S2_m_R1(3))

            CALL timeset(routineN//"_3c_x_G", handle2)

            CALL get_t_3c_int(t_3c_int, bs_env, i_cell_m_R1, i_cell_S1_p_S2_m_R1)

            ! M_λ0,νS1,PR1 = sum_µS2 ( λ0   µS1-S2    | PR1 ) G^occ/vir_µν^S2(i|τ|)
            !              = sum_µS2 ( λ-R1 µS1-S2-R1 | P0  ) G^occ/vir_µν^S2(i|τ|)
            ! for ΔR = S_1 - R_1
            CALL dbt_contract(alpha=1.0_dp, &
                              tensor_1=G_S(i_cell_S2), &
                              tensor_2=t_3c_int, &
                              beta=1.0_dp, &
                              tensor_3=t_G, &
                              contract_1=[2], notcontract_1=[1], map_1=[3], &
                              contract_2=[3], notcontract_2=[1, 2], map_2=[1, 2], &
                              filter_eps=bs_env%eps_filter)

            CALL timestop(handle2)

         END DO ! i_cell_S2

         CALL dbt_copy(t_G, t_G_2, order=[1, 3, 2], move_data=.TRUE.)

         CALL timeset(routineN//"_contract", handle2)

         DO i_cell_R = 1, bs_env%nimages_scf_desymm

            cell_R = bs_env%kpoints_scf_desymm%index_to_cell(i_cell_R, 1:3)

            ! R_1 - R
            CALL add_R(cell_R1, -cell_R, bs_env%index_to_cell_3c, cell_R1_minus_R, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_R1_minus_R)
            IF (.NOT. cell_found) CYCLE

            ! S_1 - R
            CALL add_R(cell_S1, -cell_R, bs_env%index_to_cell_3c, cell_S1_minus_R, &
                       cell_found, bs_env%cell_to_index_3c, i_cell_S1_minus_R)
            IF (.NOT. cell_found) CYCLE

            ! Σ_λσ^R = sum_PR1νS1 M^G_λ0,νS1,PR1 M^W_σR,νS1,PR1, where
            ! M^G_λ0,νS1,PR1 = sum_µS2 (λ0 µS1-S2 | PR1) G_µν^S2
            ! M^W_σR,νS1,PR1 = sum_QR2 (σR νS1 | QR1-R2) W_PQ^R2 = M^W_σ0,νS1-R,PR1-R
            CALL dbt_contract(alpha=sign_Sigma, &
                              tensor_1=t_G_2, &
                              tensor_2=t_W(i_cell_S1_minus_R, i_cell_R1_minus_R), &
                              beta=1.0_dp, &
                              tensor_3=Sigma_R(i_cell_R), &
                              contract_1=[1, 2], notcontract_1=[3], map_1=[1], &
                              contract_2=[1, 2], notcontract_2=[3], map_2=[2], &
                              filter_eps=bs_env%eps_filter)

         END DO ! i_cell_R

         CALL dbt_clear(t_G_2)

         CALL timestop(handle2)

      END DO ! i_cell_R1

      ! release memory
      IF (clear_t_W) THEN
         DO i_cell_S1 = 1, bs_env%nimages_3c
            DO i_cell_R1 = 1, bs_env%nimages_3c
               CALL dbt_clear(t_W(i_cell_S1, i_cell_R1))
            END DO
         END DO
      END IF

      CALL dbt_destroy(t_G)
      CALL dbt_destroy(t_G_2)
      CALL dbt_destroy(t_3c_int)

      CALL timestop(handle)

   END SUBROUTINE contract_to_Sigma

! **************************************************************************************************
!> \brief ...
!> \param fm_W_R ...
!> \param W_R ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE fm_MWM_R_t_to_local_tensor_W_R(fm_W_R, W_R, bs_env)
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_W_R
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: W_R
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_MWM_R_t_to_local_tensor_W_R'

      INTEGER                                            :: handle, i_cell_R

      CALL timeset(routineN, handle)

      ! communicate fm_W_R to tensor W_R; full replication in tensor group
      DO i_cell_R = 1, bs_env%nimages_scf_desymm
         CALL fm_to_local_tensor(fm_W_R(i_cell_R), bs_env%mat_RI_RI%matrix, &
                                 bs_env%mat_RI_RI_tensor%matrix, W_R(i_cell_R), bs_env)
      END DO

      CALL timestop(handle)

   END SUBROUTINE fm_MWM_R_t_to_local_tensor_W_R

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE compute_QP_energies(bs_env)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_QP_energies'

      INTEGER                                            :: handle, ikp, ispin, j_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Sigma_x_ikp_n
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Sigma_c_ikp_n_freq, Sigma_c_ikp_n_time
      TYPE(cp_cfm_type)                                  :: cfm_mo_coeff

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_mo_coeff, bs_env%fm_s_Gamma%matrix_struct)
      ALLOCATE (Sigma_x_ikp_n(bs_env%n_ao))
      ALLOCATE (Sigma_c_ikp_n_time(bs_env%n_ao, bs_env%num_time_freq_points, 2))
      ALLOCATE (Sigma_c_ikp_n_freq(bs_env%n_ao, bs_env%num_time_freq_points, 2))

      DO ispin = 1, bs_env%n_spin

         DO ikp = 1, bs_env%nkp_bs_and_DOS

            ! 1. get C_µn(k)
            CALL cp_cfm_to_cfm(bs_env%cfm_mo_coeff_kp(ikp, ispin), cfm_mo_coeff)

            ! 2. Σ^x_µν(k) = sum_R Σ^x_µν^R e^ikR
            !    Σ^x_nn(k) = sum_µν C^*_µn(k) Σ^x_µν(k) C_νn(k)
            CALL trafo_to_k_and_nn(bs_env%fm_Sigma_x_R, Sigma_x_ikp_n, cfm_mo_coeff, bs_env, ikp)

            ! 3. Σ^c_µν(k,+/-i|τ_j|) = sum_R Σ^c_µν^R(+/-i|τ_j|) e^ikR
            !    Σ^c_nn(k,+/-i|τ_j|) = sum_µν C^*_µn(k) Σ^c_µν(k,+/-i|τ_j|) C_νn(k)
            DO j_t = 1, bs_env%num_time_freq_points
               CALL trafo_to_k_and_nn(bs_env%fm_Sigma_c_R_pos_tau(:, j_t, ispin), &
                                      Sigma_c_ikp_n_time(:, j_t, 1), cfm_mo_coeff, bs_env, ikp)
               CALL trafo_to_k_and_nn(bs_env%fm_Sigma_c_R_neg_tau(:, j_t, ispin), &
                                      Sigma_c_ikp_n_time(:, j_t, 2), cfm_mo_coeff, bs_env, ikp)
            END DO

            ! 4. Σ^c_nn(k_i,iω) = ∫ from -∞ to ∞ dτ e^-iωτ Σ^c_nn(k_i,iτ)
            CALL time_to_freq(bs_env, Sigma_c_ikp_n_time, Sigma_c_ikp_n_freq, ispin)

            ! 5. Analytic continuation Σ^c_nn(k_i,iω) -> Σ^c_nn(k_i,ϵ) and
            !    ϵ_nk_i^GW = ϵ_nk_i^DFT + Σ^c_nn(k_i,ϵ) + Σ^x_nn(k_i) - v^xc_nn(k_i)
            CALL analyt_conti_and_print(bs_env, Sigma_c_ikp_n_freq, Sigma_x_ikp_n, &
                                        bs_env%v_xc_n(:, ikp, ispin), &
                                        bs_env%eigenval_scf(:, ikp, ispin), ikp, ispin)

         END DO ! ikp

      END DO ! ispin

      CALL get_all_VBM_CBM_bandgaps(bs_env)

      CALL cp_cfm_release(cfm_mo_coeff)

      CALL timestop(handle)

   END SUBROUTINE compute_QP_energies

! **************************************************************************************************
!> \brief ...
!> \param fm_rs ...
!> \param array_ikp_n ...
!> \param cfm_mo_coeff ...
!> \param bs_env ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE trafo_to_k_and_nn(fm_rs, array_ikp_n, cfm_mo_coeff, bs_env, ikp)
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      REAL(KIND=dp), DIMENSION(:)                        :: array_ikp_n
      TYPE(cp_cfm_type)                                  :: cfm_mo_coeff
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'trafo_to_k_and_nn'

      INTEGER                                            :: handle, n_ao
      TYPE(cp_cfm_type)                                  :: cfm_ikp, cfm_tmp
      TYPE(cp_fm_type)                                   :: fm_ikp_re

      CALL timeset(routineN, handle)

      CALL cp_cfm_create(cfm_ikp, cfm_mo_coeff%matrix_struct)
      CALL cp_cfm_create(cfm_tmp, cfm_mo_coeff%matrix_struct)
      CALL cp_fm_create(fm_ikp_re, cfm_mo_coeff%matrix_struct)

      ! Σ_µν(k_i) = sum_R e^ik_iR Σ_µν^R
      CALL fm_trafo_rs_to_ikp(cfm_ikp, fm_rs, bs_env%kpoints_DOS, ikp)

      n_ao = bs_env%n_ao

      ! Σ_nm(k_i) = sum_µν C^*_µn(k_i) Σ_µν(k_i) C_νn(k_i)
      CALL parallel_gemm('N', 'N', n_ao, n_ao, n_ao, z_one, cfm_ikp, cfm_mo_coeff, z_zero, cfm_tmp)
      CALL parallel_gemm('C', 'N', n_ao, n_ao, n_ao, z_one, cfm_mo_coeff, cfm_tmp, z_zero, cfm_ikp)

      ! get Σ_nn(k_i) which is a real quantity as Σ^x and Σ^c(iτ) is Hermitian
      CALL cp_cfm_to_fm(cfm_ikp, fm_ikp_re)
      CALL cp_fm_get_diag(fm_ikp_re, array_ikp_n)

      CALL cp_cfm_release(cfm_ikp)
      CALL cp_cfm_release(cfm_tmp)
      CALL cp_fm_release(fm_ikp_re)

      CALL timestop(handle)

   END SUBROUTINE trafo_to_k_and_nn

END MODULE gw_small_cell_full_kp
