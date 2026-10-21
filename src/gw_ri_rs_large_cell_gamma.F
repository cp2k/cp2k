!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GW using the RI-RS approximation for large periodic cells at Gamma
!> \par History
!>      04.2026 created [Ritaj Tyagi]
! **************************************************************************************************
MODULE gw_ri_rs_large_cell_Gamma
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_get_info, dbcsr_p_type, &
        dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE gw_compute_V_and_W,              ONLY: get_W_MIC
   USE gw_tensor_large_cell_Gamma,      ONLY: G_occ_vir,&
                                              compute_QP_energies,&
                                              delete_unnecessary_files,&
                                              fill_fm_Sigma_c_Gamma_time
   USE gw_utils,                        ONLY: de_init_bs_env
   USE gw_utils_dbcsr,                  ONLY: dbcsr_contract_ABA,&
                                              hadamard_product
   USE gw_utils_fm,                     ONLY: fm_contract_ABA
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_matrices,                  ONLY: release_dbcsr_topology_and_matrices,&
                                              setup_square_topology
   USE ri_rs_setup,                     ONLY: ensure_ri_rs_representation
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_large_cell_Gamma'

   PUBLIC :: gw_calc_ri_rs_large_cell_Gamma

CONTAINS

! **************************************************************************************************
!> \brief GW calculation using RI-RS formalism for large periodic cells at Gamma
!> \param qs_env ...
!> \param bs_env Band-structure environment containing GW parameters.
! **************************************************************************************************

   SUBROUTINE gw_calc_ri_rs_large_cell_Gamma(qs_env, bs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'gw_calc_ri_rs_large_cell_Gamma'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma, fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CALL timeset(routineN, handle)

      ! Build the common grid, phi and Z representation with the system's boundary conditions.
      CALL ensure_ri_rs_representation(qs_env, bs_env)

       !!========================================================================
       !! 4. Compute Independent-Particle Polarizability (χ)
       !!    G^occ_µλ(i|τ|)  = sum_n^occ C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^vir_µλ(i|τ|)  = sum_n^vir C_µn e^(-|(ϵ_n-ϵ_F)τ|) C_λn
       !!    G^occ_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')
       !!    G^vir_ll'(i|τ|) = sum_µν ϕ_µ(r_l) G^vir_µν ϕ_ν(r_l')
       !!    χ_ll'(iτ)       = G^occ_ll'(i|τ|) * G^vir_ll'(i|τ|)
       !!    χ_PQ(iτ)        = sum_ll' Z_lP χ_ll'(iτ) Z_l'Q
       !!========================================================================
      CALL get_mat_chi_Gamma_tau(bs_env, bs_env%mat_chi_Gamma_tau, &
                                 bs_env%ri_rs%mat_phi_mu_l, bs_env%ri_rs%mat_Z_lP)

       !!========================================================================
       !! 5. Compute Screened Interaction (W^MIC)
       !!    χ_PQ(iτ) -> χ_PQ(iω) -> ε_PQ(iω) -> W_PQ(iω) -> W^MIC_PQ(iτ)
       !!========================================================================
      CALL get_W_MIC(bs_env, qs_env, bs_env%mat_chi_Gamma_tau, fm_W_time)

       !!========================================================================
       !! 6. Compute Exact Exchange Self-Energy (Σ^x)
       !!    D_µν          = sum_n^occ C_µn C_νn
       !!    D_ll'         = sum_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
       !!    V^trunc_ll'   = sum_PQ Z_lP V^trunc_PQ Z_l'Q
       !!    Σ^x_ll'       = D_ll' * V^trunc_ll'
       !!    Σ^x_λσ(k=0)   = -sum_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_x(bs_env, qs_env, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_x_Gamma)

       !!========================================================================
       !! 7. Compute Correlation Self-Energy (Σ^c)
       !!    W^MIC_ll'(iτ) =  sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
       !!    Σ^c_ll'(iτ)   = -G^occ_ll'(i|τ|) * W^MIC_ll'(iτ), for τ < 0
       !!    Σ^c_ll'(iτ)   =  G^vir_ll'(i|τ|) * W^MIC_ll'(iτ), for τ > 0
       !!    Σ^c_λσ(iτ)    =  sum_ll' ϕ_λ(r_l) Σ^c_ll'(iτ) ϕ_σ(r_l')
       !!========================================================================
      CALL compute_Sigma_c(bs_env, fm_W_time, bs_env%ri_rs%mat_phi_mu_l, &
                           bs_env%ri_rs%mat_Z_lP, fm_Sigma_c_Gamma_time)

       !!========================================================================
       !! 8. Compute Quasiparticle Energies
       !!    Σ^c_λσ(iτ) -> Σ^c_nn(ϵ)
       !!    ϵ_nk^GW = ϵ_nk^DFT + Σ^c_nn(ϵ) + Σ^x_nn - v^xc_nn
       !!========================================================================
      CALL compute_QP_energies(bs_env, qs_env, fm_Sigma_x_Gamma, fm_Sigma_c_Gamma_time)

      CALL de_init_bs_env(qs_env, bs_env)

      CALL timestop(handle)

   END SUBROUTINE gw_calc_ri_rs_large_cell_Gamma

! **************************************************************************************************
!> \brief Computes the χ(iτ, k=0) matrix
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
! **************************************************************************************************

   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, mat_phi_mu_l, mat_Z_lP)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i, i_t, ispin, npcol
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, blk_grid, dist_col_grid, &
                                                            dist_row_grid
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_grid_grid, dist_phi
      TYPE(dbcsr_type)                                   :: matrix_chi_grid, matrix_chi_grid_spin, &
                                                            matrix_G_occ_grid, matrix_G_vir_grid

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES
      ! =========================================================================
      CALL dbcsr_get_info(mat_phi_mu_l, distribution=dist_phi, row_blk_size=blk_grid, col_blk_size=blk_ao)
      CALL dbcsr_distribution_get(dist_phi, row_dist=dist_row_grid, npcols=npcol)

      ! Build a perfectly safe column distribution for the Grid dimension
      ALLOCATE (dist_col_grid(SIZE(blk_grid)))
      DO i = 1, SIZE(blk_grid)
         dist_col_grid(i) = MOD(i - 1, npcol)
      END DO

      CALL dbcsr_distribution_new(dist_grid_grid, template=dist_phi, &
                                  row_dist=dist_row_grid, col_dist=dist_col_grid)

      CALL dbcsr_create(matrix_G_occ_grid, "G_occ_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_G_vir_grid, "G_vir_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_chi_grid, "chi_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
      CALL dbcsr_create(matrix_chi_grid_spin, "chi_grid_spin", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

      ! =========================================================================
      ! 2. MAIN IMAGINARY TIME LOOP
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         t1 = m_walltime()

         tau = bs_env%time_frequency_grid%imaginary_time(i_t)
         CALL dbcsr_set(matrix_chi_grid, 0.0_dp)

         ! ----------------------------------------------------------------------
         ! A. SPIN LOOP (Allocations safely encapsulated in wrappers)
         ! ----------------------------------------------------------------------
         DO ispin = 1, bs_env%n_spin

            ! G^occ_µλ(i|τ|,k=0) = sum_n^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^occ_ll'(i|τ|,k=0) = sum_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, &
                              matrix_G_occ_grid, bs_env%eps_filter)

            ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^vir_ll'(i|τ|,k=0) = sum_µν ϕ_µ(r_l) G^vir_µν ϕ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, &
                              matrix_G_vir_grid, bs_env%eps_filter)

            ! -------------------------------------------------------------------
            ! B. ELEMENT-WISE HADAMARD PRODUCT
            ! -------------------------------------------------------------------
            ! χ_ll'(iτ,k=0) = G^occ_ll'(i|τ|,k=0) * G^vir_ll'(i|τ|,k=0)
            CALL hadamard_product(matrix_G_occ_grid, matrix_G_vir_grid, matrix_chi_grid_spin, bs_env%spin_degeneracy)

            ! Accumulate spin contributions
            CALL dbcsr_add(matrix_chi_grid, matrix_chi_grid_spin, 1.0_dp, 1.0_dp)

         END DO ! ispin

         ! ----------------------------------------------------------------------
         ! C. TRANSFORM TO AUXILIARY BASIS & EXPORT DIRECTLY
         ! χ_aux        = Z^T * χ_grid * Z
         ! χ_PQ(iτ,k=0) = sum_ll' Z_lP χ_ll'(iτ,k=0) Z_l'Q
         ! Result is dumped directly into the final array mat_chi_Gamma_tau!
         ! ----------------------------------------------------------------------
         CALL dbcsr_contract_ABA("T", "N", mat_Z_lP, matrix_chi_grid, &
                                 mat_chi_Gamma_tau(i_t)%matrix, bs_env%eps_filter)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I13,A,I3,A,F7.1,A)') &
               'Computed χ(iτ,k=0) for time point', i_t, ' /', bs_env%num_time_freq_points, &
               ',    Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! i_t

      ! =========================================================================
      ! 3. FINAL CLEANUP
      ! =========================================================================
      CALL dbcsr_release(matrix_G_occ_grid)
      CALL dbcsr_release(matrix_G_vir_grid)
      CALL dbcsr_release(matrix_chi_grid)
      CALL dbcsr_release(matrix_chi_grid_spin)
      CALL dbcsr_distribution_release(dist_grid_grid)
      DEALLOCATE (dist_col_grid)

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE get_mat_chi_Gamma_tau

! **************************************************************************************************
!> \brief Computes Green's Function in grid basis
!> \param bs_env ...
!> \param tau ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
!> \param mat_phi_mu_l ...
!> \param matrix_G_grid ...
!> \param eps_filter ...
! **************************************************************************************************

   SUBROUTINE build_G_grid(bs_env, tau, ispin, occ, vir, mat_phi_mu_l, matrix_G_grid, eps_filter)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occ, vir
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, matrix_G_grid
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_G_grid'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      TYPE(cp_fm_type), POINTER                          :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao
      TYPE(dbcsr_type)                                   :: matrix_G_ao

      CALL timeset(routineN, handle)

      ! 1. Select the correct FM matrix based on occ/vir flags
      IF (occ) THEN
         fm_G => bs_env%fm_Gocc
      ELSE
         fm_G => bs_env%fm_Gvir
      END IF

      ! 2. Compute Dense FM Green's Function
      ! G^occ/vir_µλ(i|τ|,k=0) = sum_G^occ/vir_µλn^occ/vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
      CALL G_occ_vir(bs_env, tau, fm_G, ispin, occ=occ, vir=vir)

      ! 3. Setup AO DBCSR Topology and Create Matrix dynamically
      CALL setup_square_topology(mat_phi_mu_l, 'COL', dist_ao_ao, blk_ao, dist_row_ao)

      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      ! 4. Convert FM to Sparse DBCSR
      CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)

      ! 5. Transform to Grid Basis: G_grid = phi * G_ao * phi^T
      ! G^occ/vir_ll'(i|τ|,k=0) = sum_µν ϕ_µ(r_l) G^occ/vir_µν ϕ_ν(r_l')
      CALL dbcsr_contract_ABA("N", "T", mat_phi_mu_l, matrix_G_ao, matrix_G_grid, &
                              eps_filter)

      ! 6. Release AO matrix and topology
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao, m1=matrix_G_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_grid

! **************************************************************************************************
!> \brief Computes the exact exchange part of the GW self-energy
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_x_Gamma ...
! **************************************************************************************************

   SUBROUTINE compute_Sigma_x(bs_env, qs_env, mat_phi_mu_l, mat_Z_lP, fm_Sigma_x_Gamma)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, ispin
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_Vtr_Gamma
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_type)                                   :: mat_Sigma_x_Gamma, matrix_D_grid, &
                                                            matrix_Sigma_x_grid, matrix_V_aux, &
                                                            matrix_V_grid

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      CALL dbcsr_create(mat_Sigma_x_Gamma, template=bs_env%mat_ao_ao%matrix)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES
      ! =========================================================================
      CALL setup_square_topology(mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, dist_col_grid)
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      ! =========================================================================
      ! 2. COMPUTE V^tr_ll'
      ! =========================================================================
      CALL RI_2c_integral_mat(qs_env, fm_Vtr_Gamma, bs_env%fm_RI_RI%matrix_struct, bs_env%n_RI, &
                              bs_env%trunc_coulomb)

      ! M^-1(k=0) V^tr(τ) M^-1(k=0) -> fm_Vtr_Gamma
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_Vtr_Gamma(:, 1))

      CALL dbcsr_create(matrix_V_aux, "V_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
      CALL dbcsr_create(matrix_V_grid, "V_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

      CALL copy_fm_to_dbcsr(fm_Vtr_Gamma(1, 1), matrix_V_aux, keep_sparsity=.FALSE.)

      ! V^tr_ll' = sum_PQ Z_lP V^trunc_PQ Z_l'Q
      CALL dbcsr_contract_ABA("N", "T", mat_Z_lP, matrix_V_aux, matrix_V_grid, &
                              bs_env%eps_filter)
      CALL dbcsr_release(matrix_V_aux)

      ! =========================================================================
      ! 3. SPIN LOOP FOR EXACT EXCHANGE
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! Density matrix on grid is essentially G_occ at tau = 0.0
         ! D_µν  = sum_n^occ C_µn(k=0) C_νn(k=0)
         ! D_ll' = sum_µν ϕ_µ(r_l) D_µν ϕ_ν(r_l')
         CALL dbcsr_create(matrix_D_grid, "D_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
         CALL build_G_grid(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_D_grid, bs_env%eps_filter)

         ! Element-wise Hadamard product: Σ^x_grid = D_grid ◦ V_grid
         ! Σ^x_ll' = D_ll' * V^tr_ll'
         CALL dbcsr_create(matrix_Sigma_x_grid, template=matrix_V_grid)
         CALL hadamard_product(matrix_D_grid, matrix_V_grid, matrix_Sigma_x_grid, 1.0_dp)

         CALL dbcsr_release(matrix_D_grid)

         ! Transform back to AO basis: Σ^x_ao = -1.0 * phi^T * Σ^x_grid * phi
         ! Σ^x_λσ(k=0)   = -sum_ll' ϕ_λ(r_l) Σ^x_ll' ϕ_σ(r_l')
         CALL dbcsr_contract_ABA("T", "N", mat_phi_mu_l, matrix_Sigma_x_grid, &
                                 mat_Sigma_x_Gamma, bs_env%eps_filter)
         CALL dbcsr_scale(mat_Sigma_x_Gamma, -1.0_dp)

         CALL dbcsr_release(matrix_Sigma_x_grid)

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(mat_Sigma_x_Gamma, fm_Sigma_x_Gamma(ispin))

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      ! =========================================================================
      ! 4. CLEANUP
      ! =========================================================================
      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid, &
                                               m1=mat_Sigma_x_Gamma, m2=matrix_V_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL cp_fm_release(fm_Vtr_Gamma)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy
!> \param bs_env ...
!> \param fm_W_time ...
!> \param mat_phi_mu_l ...
!> \param mat_Z_lP ...
!> \param fm_Sigma_c_Gamma_time ...
! **************************************************************************************************

   SUBROUTINE compute_Sigma_c(bs_env, fm_W_time, mat_phi_mu_l, mat_Z_lP, fm_Sigma_c_Gamma_time)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mat_phi_mu_l, mat_Z_lP
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      INTEGER                                            :: handle, i_t, ispin
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, blk_grid, dist_col_grid, &
                                                            dist_row_aux
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux, dist_grid_grid
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_Sigma_neg_tau, mat_Sigma_pos_tau
      TYPE(dbcsr_type) :: matrix_G_occ_grid, matrix_G_vir_grid, matrix_Sigma_neg_grid, &
         matrix_Sigma_pos_grid, matrix_W_aux, matrix_W_grid

      CALL timeset(routineN, handle)

      ! =========================================================================
      ! 1. SETUP CORE TOPOLOGIES AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(mat_phi_mu_l, 'ROW', dist_grid_grid, blk_grid, dist_col_grid)
      CALL setup_square_topology(mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      ! Pre-allocate local DBCSR matrices to act as targets for final output
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

      ! =========================================================================
      ! 2. MAIN IMAGINARY TIME LOOP
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)

         ! -------------------------------------------------------------------
         ! Compute W_grid = Z * W_aux * Z^T
         ! -------------------------------------------------------------------
         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, blk_aux, blk_aux)
         CALL dbcsr_create(matrix_W_grid, "W_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

         CALL copy_fm_to_dbcsr(fm_W_time(i_t), matrix_W_aux, keep_sparsity=.FALSE.)

         ! W^MIC_ll'(iτ,k=0) = sum_PQ Z_lP W^MIC_PQ(iτ) Z_l'Q
         CALL dbcsr_contract_ABA("N", "T", mat_Z_lP, matrix_W_aux, matrix_W_grid, &
                                 bs_env%eps_filter)

         CALL dbcsr_release(matrix_W_aux) ! Clean up aux basis immediately

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! -------------------------------------------------------------------
            ! A. Transform Green's Functions to the Grid
            ! -------------------------------------------------------------------
            CALL dbcsr_create(matrix_G_occ_grid, "G_occ_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)
            CALL dbcsr_create(matrix_G_vir_grid, "G_vir_grid", dist_grid_grid, dbcsr_type_no_symmetry, blk_grid, blk_grid)

            ! G^occ_µλ(i|τ|,k=0) = sum_G^occ_µλn^occ C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^occ_ll'(i|τ|,k=0) = sum_µν ϕ_µ(r_l) G^occ_µν ϕ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .TRUE., .FALSE., mat_phi_mu_l, matrix_G_occ_grid, bs_env%eps_filter)

            ! G^vir_µλ(i|τ|,k=0) = sum_n^vir C_µn(k=0) e^(-|(ϵ_nk=0-ϵ_F)τ|) C_λn(k=0)
            ! G^vir_ll'(i|τ|,k=0) = sum_µν ϕ_µ(r_l) G^vir_µν ϕ_ν(r_l')
            CALL build_G_grid(bs_env, tau, ispin, .FALSE., .TRUE., mat_phi_mu_l, matrix_G_vir_grid, bs_env%eps_filter)

            ! -------------------------------------------------------------------
            ! B. Element-wise Hadamard Products for Sigma_c on Grid
            ! Σ_neg_grid = G_occ_grid ◦ W_grid
            ! Σ_pos_grid = G_vir_grid ◦ W_grid
            ! -------------------------------------------------------------------
            CALL dbcsr_create(matrix_Sigma_neg_grid, template=matrix_W_grid)
            CALL dbcsr_create(matrix_Sigma_pos_grid, template=matrix_W_grid)

            ! Σ^c_ll'(iτ,k=0)   = -G^occ_ll'(i|τ|,k=0) * W^MIC_ll'(iτ,k=0), for τ < 0
            CALL hadamard_product(matrix_G_occ_grid, matrix_W_grid, matrix_Sigma_neg_grid, 1.0_dp)

            ! Σ^c_ll'(iτ,k=0)   =  G^vir_ll'(i|τ|,k=0) * W^MIC_ll'(iτ,k=0), for τ > 0
            CALL hadamard_product(matrix_G_vir_grid, matrix_W_grid, matrix_Sigma_pos_grid, 1.0_dp)

            ! Instantly purge massive G_grid arrays to save memory
            CALL dbcsr_release(matrix_G_occ_grid)
            CALL dbcsr_release(matrix_G_vir_grid)

            ! -------------------------------------------------------------------
            ! C. Transform Sigma back to AO Basis
            ! Σ_AO = phi^T * Σ_grid * phi
            ! -------------------------------------------------------------------

            ! Σ^c_λσ(iτ,k=0)    = sum_ll' ϕ_λ(r_l) Σ^c_ll'(iτ,k=0) ϕ_σ(r_l'), for τ < 0
            CALL dbcsr_contract_ABA("T", "N", mat_phi_mu_l, matrix_Sigma_neg_grid, &
                                    mat_Sigma_neg_tau(i_t, ispin)%matrix, bs_env%eps_filter)
            CALL dbcsr_scale(mat_Sigma_neg_tau(i_t, ispin)%matrix, -1.0_dp)

            ! Σ^c_λσ(iτ,k=0)    = sum_ll' ϕ_λ(r_l) Σ^c_ll'(iτ,k=0) ϕ_σ(r_l'), for τ > 0
            CALL dbcsr_contract_ABA("T", "N", mat_phi_mu_l, matrix_Sigma_pos_grid, &
                                    mat_Sigma_pos_tau(i_t, ispin)%matrix, bs_env%eps_filter)

            ! Purge Grid Sigma arrays
            CALL dbcsr_release(matrix_Sigma_neg_grid)
            CALL dbcsr_release(matrix_Sigma_pos_grid)

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I10,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ,k=0) for time point ', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         ! Release the W_grid for this time point
         CALL dbcsr_release(matrix_W_grid)

      END DO ! i_t

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      ! -------------------------------------------------------------------------
      ! 3. FINALIZE AND CLEANUP
      ! -------------------------------------------------------------------------
      CALL fill_fm_Sigma_c_Gamma_time(fm_Sigma_c_Gamma_time, bs_env, &
                                      mat_Sigma_pos_tau, mat_Sigma_neg_tau)

      CALL cp_fm_release(fm_W_time)

      CALL dbcsr_deallocate_matrix_set(mat_Sigma_neg_tau)
      CALL dbcsr_deallocate_matrix_set(mat_Sigma_pos_tau)

      CALL release_dbcsr_topology_and_matrices(dist=dist_grid_grid, mapped_dist=dist_col_grid)
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      CALL delete_unnecessary_files(bs_env)
      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

END MODULE gw_ri_rs_large_cell_Gamma
