!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared RI-RS polarizability and self-energy assembly.
! **************************************************************************************************
MODULE gw_ri_rs_chi_sigma_kernels
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_distribution_type,&
                                              dbcsr_filter,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_scale,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE gw_tensor_large_cell_Gamma,      ONLY: delete_unnecessary_files
   USE gw_utils,                        ONLY: G_occ_vir
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE ri_rs_matrices,                  ONLY: build_ri_rs_bare_aux_interaction,&
                                              contract_grid_panels,&
                                              print_matrix_occupation,&
                                              release_dbcsr_topology_and_matrices,&
                                              resolve_grid_panels,&
                                              setup_square_topology
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_ri_rs_chi_sigma_kernels'
   PUBLIC :: get_mat_chi_Gamma_tau, compute_Sigma_x, compute_Sigma_c
CONTAINS

! **************************************************************************************************
!> \brief Computes χ_PQ(iτ) from the occupied and virtual Green's functions.
!> \param bs_env ...
!> \param mat_chi_Gamma_tau Response matrices χ_PQ(iτ) to be computed
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
! **************************************************************************************************
   SUBROUTINE get_mat_chi_Gamma_tau(bs_env, mat_chi_Gamma_tau, cutoff)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_mat_chi_Gamma_tau'

      INTEGER                                            :: handle, i_t, ispin, n_panels
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      REAL(KIND=dp)                                      :: first_grid_occ, t1, tau
      TYPE(dbcsr_type)                                   :: matrix_G_occ_ao(1), matrix_G_vir_ao

      CALL timeset(routineN, handle)

      ! Panel boundaries for the grid-streaming contraction.
      ! The panels are identical for χ, Σ^x and Σ^c, so the count is reported once here for all three stages.

      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, cutoff)

      n_panels = SIZE(pan_first)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T74,I9)') &
            'Number of batches for χ, Σ matrices', n_panels
         WRITE (bs_env%unit_nr, '(A)') ' '
         CALL m_flush(bs_env%unit_nr)
      END IF

      ! =========================================================================
      ! IMAGINARY TIME LOOP
      ! χ_PQ(iτ) = Σ_s g_s · Z^T ( (ϕ G^occ_s ϕ^T) ∘ (ϕ G^vir_s ϕ^T) ) Z
      ! (g_s = spin degeneracy)
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         t1 = m_walltime()
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)

         DO ispin = 1, bs_env%n_spin

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_occ_ao(1))
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_vir_ao)

            ! χ_PQ += g_s · Z^T ( (ϕ G^occ ϕ^T) ∘ (ϕ G^vir ϕ^T) ) Z
            IF (i_t == 1) THEN
               CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_G_occ_ao, &
                                         L_B=bs_env%ri_rs%mat_phi_mu_l, M_B=matrix_G_vir_ao, &
                                         L_out=bs_env%ri_rs%mat_Z_lP, mat_out=mat_chi_Gamma_tau(i_t:i_t), &
                                         scale=bs_env%spin_degeneracy, eps=bs_env%eps_filter, &
                                         para_env=bs_env%para_env, &
                                         pan_first=pan_first, pan_last=pan_last, &
                                         lb_eq_la=.TRUE., lout_eq_la=.FALSE., &
                                         zero_out=(ispin == 1), &
                                         centroids=bs_env%ri_rs%chunk_centroids, &
                                         cutoff=cutoff, grid_occupation=first_grid_occ)
            ELSE
               CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_G_occ_ao, &
                                         L_B=bs_env%ri_rs%mat_phi_mu_l, M_B=matrix_G_vir_ao, &
                                         L_out=bs_env%ri_rs%mat_Z_lP, mat_out=mat_chi_Gamma_tau(i_t:i_t), &
                                         scale=bs_env%spin_degeneracy, eps=bs_env%eps_filter, &
                                         para_env=bs_env%para_env, &
                                         pan_first=pan_first, pan_last=pan_last, &
                                         lb_eq_la=.TRUE., lout_eq_la=.FALSE., &
                                         zero_out=(ispin == 1), &
                                         centroids=bs_env%ri_rs%chunk_centroids, &
                                         cutoff=cutoff)
            END IF

            CALL dbcsr_release(matrix_G_occ_ao(1))
            CALL dbcsr_release(matrix_G_vir_ao)

         END DO ! ispin

         ! Sparsity reports
         IF (i_t == 1) THEN
            CALL print_matrix_occupation(bs_env%ri_rs%mat_Z_lP, 'Z_lP', bs_env%para_env, bs_env%unit_nr, &
                                         bs_env%dbcsr_msg_elem_limit)
            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,T73,F7.2,A)') &
                  'Percentage of non-zero matrix elements in G_ll'', χ_ll'', W_ll''', &
                  first_grid_occ*100.0_dp, ' %'
               CALL m_flush(bs_env%unit_nr)
            END IF
            CALL print_matrix_occupation(mat_chi_Gamma_tau(i_t)%matrix, 'χ_PQ', &
                                         bs_env%para_env, bs_env%unit_nr, bs_env%dbcsr_msg_elem_limit, &
                                         suffix=' for time point 1')
            IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '
         END IF

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
!> \brief Computes the exact-exchange part of the GW self-energy:
!>          D_μν        = Σ_n^occ C_μn C_νn                     (density matrix = G^occ at τ=0)
!>          V^tr_PQ     = M^-1 (P|Q)_trunc M^-1                 (truncated Coulomb, RI basis)
!>          Σ^x_λσ(k=0) = -Σ_ll' ϕ_λ(r_l) [ (ϕ D ϕ^T)_ll' ∘ (Z V^tr Z^T)_ll' ] ϕ_σ(r_l')
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_Sigma_x_Gamma ...
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
! **************************************************************************************************
   SUBROUTINE compute_Sigma_x(bs_env, qs_env, fm_Sigma_x_Gamma, cutoff)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_Sigma_x_Gamma
      REAL(KIND=dp), INTENT(IN)                          :: cutoff

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_x'

      INTEGER                                            :: handle, ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      REAL(KIND=dp)                                      :: t1
      TYPE(dbcsr_p_type)                                 :: mat_Sigma_x_Gamma(1)
      TYPE(dbcsr_type)                                   :: matrix_D_ao(1), matrix_V_aux
      TYPE(dbcsr_type), TARGET                           :: matrix_Sigma_x

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      ALLOCATE (fm_Sigma_x_Gamma(bs_env%n_spin))
      DO ispin = 1, bs_env%n_spin
         CALL cp_fm_create(fm_Sigma_x_Gamma(ispin), bs_env%fm_s_Gamma%matrix_struct)
      END DO

      mat_Sigma_x_Gamma(1)%matrix => matrix_Sigma_x
      CALL dbcsr_create(matrix_Sigma_x, template=bs_env%mat_ao_ao%matrix)

      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, cutoff)

      ! =========================================================================
      ! 1. COMPUTE V^tr_PQ (RI x RI)
      ! =========================================================================
      CALL build_ri_rs_bare_aux_interaction(qs_env, bs_env%fm_Minv_Gamma, bs_env%trunc_coulomb, &
                                            bs_env%auto_ri, bs_env%ri_rs%mat_Z_lP, matrix_V_aux)
      CALL dbcsr_filter(matrix_V_aux, bs_env%eps_filter)

      ! =========================================================================
      ! 2. SPIN LOOP FOR EXACT EXCHANGE
      ! Σ^x_λσ = -Σ_ll' ϕ_λ(r_l) ( D_ll' V^tr_ll' ) ϕ_σ(r_l')
      !        = -ϕ^T ( (ϕ D ϕ^T) ∘ (Z V^tr Z^T) ) ϕ
      ! =========================================================================
      DO ispin = 1, bs_env%n_spin

         ! AO-space density matrix D_µν = G^occ at τ = 0
         CALL build_G_ao(bs_env, 0.0_dp, ispin, .TRUE., .FALSE., &
                         bs_env%ri_rs%mat_phi_mu_l, matrix_D_ao(1))

         CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_D_ao, &
                                   L_B=bs_env%ri_rs%mat_Z_lP, M_B=matrix_V_aux, &
                                   L_out=bs_env%ri_rs%mat_phi_mu_l, mat_out=mat_Sigma_x_Gamma, &
                                   scale=1.0_dp, eps=bs_env%eps_filter, &
                                   para_env=bs_env%para_env, &
                                   pan_first=pan_first, pan_last=pan_last, &
                                   lb_eq_la=.FALSE., lout_eq_la=.TRUE., zero_out=.TRUE., &
                                   centroids=bs_env%ri_rs%chunk_centroids, &
                                   cutoff=cutoff)
         CALL dbcsr_scale(matrix_Sigma_x, -1.0_dp)

         CALL dbcsr_release(matrix_D_ao(1))

         ! Data I/O and Export to CP2K Full Matrices
         CALL copy_dbcsr_to_fm(matrix_Sigma_x, fm_Sigma_x_Gamma(ispin))

      END DO ! ispin

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed Σ^x(k=0),', ' Execution time', m_walltime() - t1, ' s'
         WRITE (bs_env%unit_nr, '(A)') ' '
      END IF

      ! =========================================================================
      ! 3. CLEANUP
      ! =========================================================================
      CALL dbcsr_release(matrix_V_aux)
      CALL dbcsr_release(matrix_Sigma_x)

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_x

! **************************************************************************************************
!> \brief Computes the correlation part of the GW self-energy on the imaginary-time grid:
!>          Σ^c_λσ(iτ<0) = -Σ_ll' ϕ_λ(r_l) [ (ϕ G^occ ϕ^T)_ll' ∘ (Z W^MIC Z^T)_ll' ] ϕ_σ(r_l')
!>          Σ^c_λσ(iτ>0) = +Σ_ll' ϕ_λ(r_l) [ (ϕ G^vir ϕ^T)_ll' ∘ (Z W^MIC Z^T)_ll' ] ϕ_σ(r_l')
!> \param bs_env ...
!> \param fm_W_time ...
!> \param fm_Sigma_c_Gamma_time ...
!> \param cutoff grid-centroid distance cutoff; zero disables geometric truncation
!> \param keep_w_and_files keep W and scratch files for another Sigma_c evaluation
! **************************************************************************************************
   SUBROUTINE compute_Sigma_c(bs_env, fm_W_time, fm_Sigma_c_Gamma_time, cutoff, keep_w_and_files)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :, :)  :: fm_Sigma_c_Gamma_time
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      LOGICAL, INTENT(IN), OPTIONAL                      :: keep_w_and_files

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_Sigma_c'

      CHARACTER(LEN=32)                                  :: sigma_label
      INTEGER                                            :: handle, i_sector, i_t, ispin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pan_first, pan_last
      INTEGER, DIMENSION(:), POINTER                     :: blk_aux, dist_row_aux
      LOGICAL                                            :: retain_work
      REAL(KIND=dp)                                      :: t1, tau
      TYPE(dbcsr_distribution_type)                      :: dist_aux_aux
      TYPE(dbcsr_p_type)                                 :: mat_Sigma_tau(2)
      TYPE(dbcsr_type)                                   :: matrix_G_ao(2), matrix_W_aux
      TYPE(dbcsr_type), TARGET                           :: matrix_Sigma_tau(2)

      CALL timeset(routineN, handle)

      retain_work = .FALSE.
      IF (PRESENT(keep_w_and_files)) retain_work = keep_w_and_files

      ! =========================================================================
      ! 1. SETUP AUXILIARY TOPOLOGY AND PRE-ALLOCATE OUTPUT ARRAYS
      ! =========================================================================
      CALL setup_square_topology(bs_env%ri_rs%mat_Z_lP, 'COL', dist_aux_aux, blk_aux, dist_row_aux)

      CALL resolve_grid_panels(bs_env%ri_rs, bs_env%para_env, bs_env%dbcsr_msg_elem_limit, &
                               pan_first, pan_last, cutoff)

      ! Only the final FM time array persists; DBCSR outputs are streamed one time/spin at a time.
      ALLOCATE (fm_Sigma_c_Gamma_time(bs_env%num_time_freq_points, 2, bs_env%n_spin))
      DO i_sector = 1, 2
         mat_Sigma_tau(i_sector)%matrix => matrix_Sigma_tau(i_sector)
      END DO

      ! =========================================================================
      ! 2. IMAGINARY TIME LOOP
      ! Σ^c_neg_λσ(iτ) = -ϕ^T ( (ϕ G^occ ϕ^T) ∘ (Z W^MIC Z^T) ) ϕ
      ! Σ^c_pos_λσ(iτ) =  ϕ^T ( (ϕ G^vir ϕ^T) ∘ (Z W^MIC Z^T) ) ϕ
      ! =========================================================================
      DO i_t = 1, bs_env%num_time_freq_points
         tau = bs_env%time_frequency_grid%imaginary_time(i_t)

         CALL dbcsr_create(matrix_W_aux, "W_aux", dist_aux_aux, dbcsr_type_no_symmetry, &
                           blk_aux, blk_aux)
         CALL copy_fm_to_dbcsr(fm_W_time(i_t), matrix_W_aux, keep_sparsity=.FALSE.)
         CALL dbcsr_filter(matrix_W_aux, bs_env%eps_filter)

         DO ispin = 1, bs_env%n_spin
            t1 = m_walltime()

            ! AO-space Green's functions G^occ_µν, G^vir_µν (dense AO x AO, small)
            CALL build_G_ao(bs_env, tau, ispin, .TRUE., .FALSE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_ao(1))
            CALL build_G_ao(bs_env, tau, ispin, .FALSE., .TRUE., &
                            bs_env%ri_rs%mat_phi_mu_l, matrix_G_ao(2))

            DO i_sector = 1, 2
               CALL dbcsr_create(matrix_Sigma_tau(i_sector), template=bs_env%mat_ao_ao%matrix)
            END DO

            ! Σ^c_neg and Σ^c_pos in a single panel loop: W_pan = Z_panel × W × Z^T built once
            CALL contract_grid_panels(L_A=bs_env%ri_rs%mat_phi_mu_l, M_A=matrix_G_ao, &
                                      L_B=bs_env%ri_rs%mat_Z_lP, M_B=matrix_W_aux, &
                                      L_out=bs_env%ri_rs%mat_phi_mu_l, &
                                      mat_out=mat_Sigma_tau, &
                                      scale=1.0_dp, eps=bs_env%eps_filter, para_env=bs_env%para_env, &
                                      pan_first=pan_first, pan_last=pan_last, &
                                      lb_eq_la=.FALSE., lout_eq_la=.TRUE., zero_out=.TRUE., &
                                      centroids=bs_env%ri_rs%chunk_centroids, cutoff=cutoff)
            CALL dbcsr_scale(matrix_Sigma_tau(1), -1.0_dp)

            ! FM time sector 1 is positive, sector 2 is negative.
            DO i_sector = 1, 2
               CALL dbcsr_release(matrix_G_ao(i_sector))
               CALL cp_fm_create(fm_Sigma_c_Gamma_time(i_t, 3 - i_sector, ispin), &
                                 bs_env%fm_s_Gamma%matrix_struct)
               CALL copy_dbcsr_to_fm(matrix_Sigma_tau(i_sector), &
                                     fm_Sigma_c_Gamma_time(i_t, 3 - i_sector, ispin))
            END DO

            ! Report the first-time outputs before releasing the current DBCSR pair.
            IF (i_t == 1) THEN
               sigma_label = 'Σ^c'
               IF (bs_env%n_spin == 2 .AND. ispin == 1) sigma_label = 'Σ^c_α'
               IF (bs_env%n_spin == 2 .AND. ispin == 2) sigma_label = 'Σ^c_β'
               CALL print_matrix_occupation(matrix_Sigma_tau(2), &
                                            TRIM(sigma_label)//'(iτ)', bs_env%para_env, bs_env%unit_nr, &
                                            bs_env%dbcsr_msg_elem_limit, suffix=' for time point 1')
               CALL print_matrix_occupation(matrix_Sigma_tau(1), &
                                            TRIM(sigma_label)//'(-iτ)', bs_env%para_env, bs_env%unit_nr, &
                                            bs_env%dbcsr_msg_elem_limit, suffix=' for time point 1')
            END IF
            DO i_sector = 1, 2
               CALL dbcsr_release(matrix_Sigma_tau(i_sector))
            END DO

            IF (bs_env%unit_nr > 0) THEN
               WRITE (bs_env%unit_nr, '(T2,A,I15,A,I3,A,F7.1,A)') &
                  'Computed Σ^c(iτ) for time point', i_t, ' /', bs_env%num_time_freq_points, &
                  ',    Execution time', m_walltime() - t1, ' s'
            END IF

         END DO ! ispin

         CALL dbcsr_release(matrix_W_aux)

      END DO ! i_t

      ! The DBCSR topology is local to this contraction.
      CALL release_dbcsr_topology_and_matrices(dist=dist_aux_aux, mapped_dist=dist_row_aux)

      ! The non-periodic evGW0 path keeps these across Sigma_c evaluations.
      IF (.NOT. retain_work) THEN
         CALL cp_fm_release(fm_W_time)
         CALL delete_unnecessary_files(bs_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_Sigma_c

! **************************************************************************************************
!> \brief Builds the AO-space Green's function operator G^occ/vir_µν (AO x AO DBCSR)
!> \param bs_env ...
!> \param tau ...
!> \param ispin ...
!> \param occ ...
!> \param vir ...
!> \param template ...
!> \param matrix_G_ao ...
! **************************************************************************************************
   SUBROUTINE build_G_ao(bs_env, tau, ispin, occ, vir, template, matrix_G_ao)

      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), INTENT(IN)                          :: tau
      INTEGER, INTENT(IN)                                :: ispin
      LOGICAL, INTENT(IN)                                :: occ, vir
      TYPE(dbcsr_type), INTENT(INOUT)                    :: template
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix_G_ao

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_G_ao'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_ao, dist_row_ao
      TYPE(cp_fm_type), POINTER                          :: fm_G
      TYPE(dbcsr_distribution_type)                      :: dist_ao_ao

      CALL timeset(routineN, handle)

      IF (occ) THEN
         fm_G => bs_env%fm_Gocc
      ELSE
         fm_G => bs_env%fm_Gvir
      END IF

      CALL G_occ_vir(bs_env, tau, fm_G, ispin, occ=occ, vir=vir)

      CALL setup_square_topology(template, 'COL', dist_ao_ao, blk_ao, dist_row_ao)
      CALL dbcsr_create(matrix_G_ao, name="G_ao", dist=dist_ao_ao, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=blk_ao, col_blk_size=blk_ao)

      CALL copy_fm_to_dbcsr(fm_G, matrix_G_ao, keep_sparsity=.FALSE.)
      CALL dbcsr_filter(matrix_G_ao, bs_env%eps_filter)

      ! release only the topology; keep matrix_G_ao for the caller
      CALL release_dbcsr_topology_and_matrices(dist=dist_ao_ao, mapped_dist=dist_row_ao)

      CALL timestop(handle)

   END SUBROUTINE build_G_ao

END MODULE gw_ri_rs_chi_sigma_kernels
